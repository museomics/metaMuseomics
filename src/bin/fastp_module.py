from concurrent.futures import ThreadPoolExecutor
import os
import sys
import glob
import subprocess
import argparse
import shutil
from pathlib import Path
from seqpy_tools import clean_and_tar, run_command, setup_logging, get_read_ids2

## Functions for fastp processing.
## This module performs initial trimming (with trimmed reads output),
## merging of paired-end reads, generation of overlap plots and final summary statistics.
### NOTE: This uses an R script with a relative filepath to be updated for package release 

def run_fastp_trim(sample_ids, r1_path, r2_path, output_dir, extra_args=None):
    """Run initial fastp trimming and filtering on paired-end reads."""

    tmp_dir = os.path.join(output_dir, "tmp")
    os.makedirs(tmp_dir, exist_ok=True)

    fastp_command = [
        "fastp", "-i", r1_path, "-I", r2_path,
        "--unpaired1", f"{tmp_dir}/{sample_ids}_unpaired_1.fq",
        "--unpaired2", f"{tmp_dir}/{sample_ids}_unpaired_2.fq",
        "--out1", f"{output_dir}/{sample_ids}_trimmed_1.fq",
        "--out2", f"{output_dir}/{sample_ids}_trimmed_2.fq",
        "--detect_adapter_for_pe",
        "--qualified_quality_phred=30", "--trim_poly_g",
        "--correction", "--dedup",
        "--html", f"{output_dir}/{sample_ids}_trim.html",
        "--json", f"{output_dir}/{sample_ids}_trim.json"
    ]

    if extra_args:
        fastp_command += extra_args

    run_command(fastp_command, f"{sample_ids}_trim")

def run_fastp_merge(sample_ids, output_dir):
    """Run fastp merge on trimmed paired-end reads
    This is designed to rely on the output from run_fastp_trim.
    Assumes run_fastp_trim naming convention for input pairs."""

    fastp_command = [
        "fastp", "--in1", f"{output_dir}/{sample_ids}_trimmed_1.fq",
        "--in2", f"{output_dir}/{sample_ids}_trimmed_2.fq",
        "--merge", "--merged_out", f"{output_dir}/{sample_ids}_merged.fq",
        "--out1", f"{output_dir}/{sample_ids}_unmerged_1.fq",
        "--out2", f"{output_dir}/{sample_ids}_unmerged_2.fq",
        "--length_required", "30",
        "--html", f"{output_dir}/{sample_ids}_merge.html",
        "--json", f"{output_dir}/{sample_ids}_merge.json"
    ]

    run_command(fastp_command, f"{sample_ids}_merge")

def generate_seqkit_stats(output_dir, stats_output):
    """Generate seqkit stats summary of all output files in output_dir"""

    fastq_files = glob.glob(os.path.join(output_dir, "*.f*q"))
    seqkit_command = ["seqkit", "stats", *fastq_files]
    with open(stats_output, "w", encoding="utf-8") as stats_file:
        subprocess.run(seqkit_command, stdout=stats_file, check=True)

def run_fastp_overlap_plot(sample_ids, r1_path, r2_path, output_dir):
    """Uses fastp merge without merge output to generate overlap plots as an html file."""

    command = [
        "fastp", "--in1", r1_path, "--in2", r2_path,
        "--stdout", "--merge", "-A", "-G", "-Q", "-L",
        "--json", "/dev/null",
        "--html", os.path.join(output_dir, f"{sample_ids}_overlaps.html")
    ]
    with open(os.devnull, 'w', encoding="utf-8") as devnull:
        subprocess.run(command, stdout=devnull)


def run_fastp_json_merge(
    json_dir,
    r_script,
    logger,
    output_file="combined_fastp_json_out.csv"
):
    """Run the R fastp JSON parser on all JSON files in a directory."""

    json_dir = Path(json_dir)
    r_script = Path(r_script)

    # Confirm script exists
    if not r_script.exists():
        logger.error(f"R script not found: {r_script}")
        raise FileNotFoundError(f"R script not found: {r_script}")

    # Collect input JSONs
    json_files = sorted(json_dir.glob("*_trim.json"))
    if not json_files:
        logger.error(f"No JSON files found in: {json_dir}")
        raise RuntimeError(f"No JSON files found in: {json_dir}")

    cmd = ["Rscript", str(r_script)] + [str(f) for f in json_files]

    logger.info(f"Running fastp JSON parser: {' '.join(cmd)}")

    subprocess.run(cmd, check=True)

    logger.info(f"Combined JSON output written to: {output_file}")

def process_sample(sample_id, r1, r2, output_dir, logger, extra_args=None):
    """The function for order of fastp operations per sample."""

    logger.info("Processing %s", sample_id)
    run_fastp_trim(sample_id, r1, r2, output_dir, extra_args)
    run_fastp_merge(sample_id, output_dir)
    run_fastp_overlap_plot(sample_id, r1, r2, output_dir)
    logger.info("Finished processing %s", sample_id)

def main(args):
    """The main function to handle argument parsing and workflow of module."""

    logger = setup_logging(log_dir="./logs", log_file="fastp_processing.log")
    os.makedirs(args.output_dir, exist_ok=True)

    if args.tracking_sheet:
        input_source = args.tracking_sheet
    else:
        input_source = args.input_dir

    if args.ids is None:
        ids = get_read_ids2(
            input_source,
            logger=logger,
            paired=True,
            prefix=args.prefix,
            suffix=args.suffix,
            column_name=args.column_name,
            sheet=args.sheet
        )
    else:
        ids = args.ids.split(",")

    jobs = []
    paired_files = {}

    ### uses ids to find pairs in the input dir but redundant because get_read_ids2 uses os.walk
    for sample_id in ids:
        id_pair = glob.glob(os.path.join(args.input_dir, "**", f"{sample_id}*"), recursive=True)
        if len(id_pair) == 2:
            r1_path, r2_path = sorted(id_pair)
            logger.info("Running files: %s and %s", r1_path, r2_path)
            paired_files[sample_id] = [r1_path, r2_path]
            jobs.append((sample_id, r1_path, r2_path))
        else:
            logger.error(
                "Error: ID %s has %d matching files.",
                sample_id,
                len(id_pair)
            )

    # Run in parallel after jobs are built
    with ThreadPoolExecutor(max_workers=args.threads) as executor:
        futures = [
            executor.submit(
                process_sample,
                sample_id,
                r1,
                r2,
                args.output_dir,
                logger,
                args.fastp_extra_args
            )
            for sample_id, r1, r2 in jobs
        ]
        for future in futures:
            future.result()

    # Generate summary stats
    logger.info("Generating summary statistics...")
    generate_seqkit_stats(args.output_dir, f"{args.output_dir}/summary_stats.txt")
    generate_seqkit_stats(f"{args.output_dir}/tmp", f"{args.output_dir}/summary_stats_intermediate.txt")

    pattern = os.path.join(args.output_dir, "**", "*_trim.json")
    json_paths = glob.glob(pattern, recursive=True)
    if not json_paths:
        raise FileNotFoundError(f"No *_trim.json files found under {args.output_dir}")
    run_fastp_json_merge(
        os.path.dirname(json_paths[0]),
        r_script="/mnt/shared/scratch/mkamouyi/private/defra-fungi/PRJEB81712/parse_fastp_json.R",
        logger=logger,
        output_file=os.path.join(args.output_dir, "combined_fastp_json_out.csv")
    )
    
    # Tar and clean
    # Remove tmp directory if it exists
    tmp_dir = os.path.join(args.output_dir, "tmp")
    if os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir)
        print(f"Removed directory: {tmp_dir}")
    logger.info("All samples processed!")

if __name__ == "__main__":
    if len(sys.argv) == 1:
        parser = argparse.ArgumentParser()
        parser.print_help()
        sys.exit(1)

    parser = argparse.ArgumentParser(
        description="Pre-process and process raw read data using fastp."
    )
    parser.add_argument(
        "--input_dir",
        type=str,
        required=True,
        help="Path to input directory containing FASTQ files."
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        required=True,
        help="Path to output directory."
    )
    parser.add_argument(
        "--prefix",
        required=False,
        help="Prefix to use to find specific files in the input directory."
    )
    parser.add_argument(
        "--ids",
        required=False,
        help="List of specific file IDs in the input directory to use. " \
        "Comma-separated. If not provided, all IDs will be used."
    )
    parser.add_argument(
        "--tracking_sheet",
        required=False,
        help="CSV file with tracking metadata."
    )
    parser.add_argument(
        "--column_name",
        required='--tracking_sheet' in sys.argv,
        help="Column name in tracking sheet with library names."
    )
    parser.add_argument(
        "--sheet",
        type=int,
        default=0,
        help="(For XLSX input only) Optional sheet index to be used as a tracking sheet. " \
        "Default is the first sheet (--sheet 0)."
    )
    parser.add_argument(
        "--suffix",
        required=False,
        help="Suffix to use to find specific files in the input directory."
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=4,
        help="Number of threads for parallel processing."
    )
    parser.add_argument(
        "--fastp_extra_args",
        nargs=argparse.REMAINDER,
        help="Extra arguments for fastp to happen before merging."
    )

    args = parser.parse_args()
    main(args)
