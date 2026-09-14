import os
import sys
import subprocess
import pathlib
import shutil
import argparse
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed
from metahist_tools import run_subprocess, repair_reads, get_read_ids2, clean_and_tar, find_paired_files2, find_single_reads, setup_logging

def run_bbduk(file_path, output_dir, temp_dir, phix_ref, logger):
    ''' Run BBDuk to remove PhiX contamination from a FASTQ file. '''

    if file_path.endswith(".gz"):
        file_path_no_gz = file_path[:-3]
        filename = pathlib.Path(file_path_no_gz).stem  # Remove .gz
        output_file = os.path.join(temp_dir, f"{filename}_nophiX.fq.gz")
    else:
        filename = pathlib.Path(file_path).stem
        output_file = os.path.join(temp_dir, f"{filename}_nophiX.fq")

    command = [
        "bbduk.sh",
        f"in={file_path}",
        f"out={output_file}",
        f"ref={phix_ref}",
        "k=31",
        "hdist=1",
        "-Xmx2g",
        f"stats={output_dir}/{filename}_nophiX_stats.txt"
    ]
    run_subprocess(command, filename)
    logger.info("Processed %s for PhiX contamination", file_path)
    return output_file

def run_bwa_mem_and_samtools(sample_id, input_file, output_dir, temp_dir, human_ref, logger):
    ''' Run BWA MEM to map reads to human reference and SAMtools to extract unmapped reads. '''

    bam_file = os.path.join(temp_dir, f"{sample_id}.bam")
    sorted_bam_file = os.path.join(temp_dir, f"{sample_id}_sorted.bam")
    unmapped_fastq = os.path.join(output_dir, f"{sample_id}.fastq")
    stats_file = os.path.join(output_dir, f"{sample_id}_human_mapping_flagstats.txt")

    try:
        logger.info("Running BWA MEM and SAMtools for %s", input_file)
        bwa_cmd = ["bwa", "mem", "-M", "-t", "8", human_ref, input_file]
        with open(bam_file, "wb") as bam_out:
            bwa_proc = subprocess.Popen(bwa_cmd, stdout=subprocess.PIPE)
            view_proc = subprocess.Popen(["samtools", "view", "-b", "-"],
                stdin=bwa_proc.stdout, stdout=bam_out)
            bwa_proc.wait()
            view_proc.communicate()

        subprocess.run(["samtools", "sort", "-o", sorted_bam_file, bam_file], check=True)
        logger.info(f"Saved {sorted_bam_file}")

        with open(unmapped_fastq, "w", encoding="utf-8") as fastq_out:
            p1 = subprocess.Popen(["samtools", "view", "-f4", sorted_bam_file],
                stdout=subprocess.PIPE)
            subprocess.run(["samtools", "fastq"], stdin=p1.stdout, stdout=fastq_out, check=True)
        with open(stats_file, "w",  encoding="utf-8") as out:
            subprocess.run(["samtools", "flagstat", "-O", "tsv", sorted_bam_file],
                stdout=out, check=True)
            logger.info(f"Finished BWA/SAMtools for {sample_id}")
    except Exception as e:
        logger.error(f"Error processing {input_file}: {e}")
        raise

def cleanup_temp_dir(temp_dir, logger):
    ''' Remove temporary directory and its contents. '''

    if os.path.exists(temp_dir):
        try:
            shutil.rmtree(temp_dir)
            logger.info(f"Temporary directory {temp_dir} cleaned up.")
        except Exception as e:
            logger.error(f"Failed to clean up temporary directory {temp_dir}: {e}")


def main(args):
    logger = setup_logging(log_dir="./logs", log_file="fastp_processing.log")

    # Setup
    os.makedirs(args.output_dir, exist_ok=True)
    temp_dir = os.path.join(args.output_dir, "tmp")
    os.makedirs(temp_dir, exist_ok=True)
    os.makedirs("./logs", exist_ok=True)

    # Step 1: Get sample prefixes from tracking sheet or prefix argument
    if args.prefix:
        prefixes = args.prefix.split(",") if args.prefix else []
    else:
        if args.tracking_sheet:
            prefixes = get_read_ids2(args.tracking_sheet,
                paired=True,
                prefix=args.prefix,
                suffix=args.suffix,
                column_name=args.column_name,
                sheet=args.tracking_sheet,
                logger=logger
            )
        else:
            prefixes = get_read_ids2(args.input_dir,
                paired=True,
                prefix=args.prefix,
                suffix=args.suffix,
                column_name=args.column_name,
                sheet=args.tracking_sheet,
                logger=logger
            )

    # Step 2: Gather input files
    if args.paired:
        logger.info("Using paired-end read mode")
        all_pairs = []
        for prefix in prefixes:
            if args.suffix:
                r1, r2 = find_paired_files2(args.input_dir, prefix=prefix, suffix=args.suffix)
            else:
                r1, r2 = find_paired_files2(args.input_dir, prefix=prefix)
            logger.info("Found paired files for prefix %s: %s , %s", prefix, r1, r2)
            if r1 and r2:
                all_pairs.append((r1, r2))
            else:
                logger.error("No valid pair found for prefix: %s", prefix)
        if not all_pairs:
            logger.error("No valid read pairs found in input directory. Exiting.")
            sys.exit(1)
        files = [r for r1, r2 in all_pairs for r in (r1, r2)]
    elif args.merged:
        logger.info("Using merged read mode")
        logger.info("Using merged read mode")
        files = find_single_reads(args.input_dir, prefixes, args.suffix)
    else:
        logger.error("Specify either --paired or --merged")
        sys.exit(1)

    # Step 3: Run BBDuk
    max_workers = args.threads or os.cpu_count()
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        bbduk_futures = [executor.submit(run_bbduk, f, args.output_dir,
            temp_dir, args.phiX_ref, logger) for f in files]
        bbduk_results = []
        for future in as_completed(bbduk_futures):
            try:
                result = future.result()
                bbduk_results.append(result)
                logger.info("Finished BBDuk for %s", result)
            except Exception as e:
                logger.error("BBDuk error: %s", e)

        # Step 4: Run BWA/SAMtools
        bwa_futures = []
        for result in bbduk_results:
            sid = pathlib.Path(result).stem.split("_nophiX")[0]
            bwa_futures.append(executor.submit(run_bwa_mem_and_samtools, sid,
                result, args.output_dir, temp_dir, args.human_ref, logger))

        for future in as_completed(bwa_futures):
            try:
                future.result()
            except Exception as e:
                logger.error("BWA/SAMtools future failed: %s", e)

    # Step 5: Repair pairs if --paired
    if args.paired:
        logger.info("Finding post-BWA paired files in output directory")

        decontaminated_pairs = find_paired_files2(
            args.output_dir,
            prefixes=prefixes,
            suffix=args.suffix
        )
        if not decontaminated_pairs:
            logger.error("No valid read pairs found in output directory. Exiting.")
            sys.exit(1)

        max_workers = min(max_workers, len(decontaminated_pairs))
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = {
                executor.submit(repair_reads, r1, r2, args.output_dir,
                    logger): os.path.basename(r1).split("_")[0]
                for r1, r2 in decontaminated_pairs
            }
            for future in as_completed(futures):
                sid = futures[future]
                try:
                    future.result()
                    logger.info("Successfully repaired ID: %s", sid)
                except Exception as e:
                    logger.error("Error repairing ID %s : %s", sid, e)

    # Step 6: Clean up
    logger.info("Compressing and cleaning up data directories")
    cleanup_temp_dir(temp_dir, logger)
    logger.info("All samples processed!")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Pre-process and decontaminate raw reads",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--input_dir", required=True, help="Directory with FASTQ files")
    parser.add_argument("--output_dir", required=True, help="Directory for outputs")
    parser.add_argument("--tracking_sheet", help="CSV or XLSX tracking sheet")
    parser.add_argument("--prefix",
        help="Sample prefix(es) to match files if not providing a sample sheet.")
    parser.add_argument("--suffix", help="Suffix pattern between sample ID and 1/2.f*q.* \
        (e.g., 'unmerged' in 'SampleA_unmerged_1.fq.gz'). \
        Can be used in conjunction with --prefix or on its own.")
    parser.add_argument("--column_name", help="Column with sample IDs in tracking sheet")
    parser.add_argument("--sheet", type=int, default=0, help="Sheet index if XLSX file is used")
    parser.add_argument("--paired", action="store_true", help="Use paired-end reads")
    parser.add_argument("--merged", action="store_true", help="Use merged reads")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads")
    parser.add_argument("--phiX_ref", default="../ref/GCA_000819615.1_ViralProj14015_genomic.fna",
        help="Relative filepath for reference PhiX genome if not using default")
    parser.add_argument("--human_ref", default="../ref/GCF_000001405.40_GRCh38.p14_genomic.fna",
        help="Relative filepath for reference Human genome if not using default")
    args = parser.parse_args()
    main(args)
