import subprocess
import os
import glob
import pgzip
import shutil
import logging
import pandas as pd
import re
import argparse
from pathlib import Path
from seqpy_tools import clean_and_tar
from seqpy_tools import pair_input_files
from seqpy_tools import check_and_handle_gunzipped
from concurrent.futures import ThreadPoolExecutor, as_completed

# Set up logger
logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

def run_cutadapt(input_files, output_dir, cutadapt_error_rate, i7_barcodes, i5_barcodes, extra_args=None):
    """Run the cutadapt command."""
    logging.info(f"Running cutadapt for pair: {input_files}...")

    basename = os.path.basename(input_files[0])
    prefix_match = re.findall(r"^(.*?)(?:_[Rr]?[12])(?:_\d{3})?\..+$", basename)
    prefix = prefix_match[0]
    
    cutadapt_command = [
        "cutadapt",
        "-e", str(cutadapt_error_rate),
        "--pair-adapters",
        "-g", f"^file:{i5_barcodes}",
        "-G", f"^file:{i7_barcodes}",
        "-o", f"{output_dir}/{{name}}_{prefix}_1.fq.gz",
        "-p", f"{output_dir}/{{name}}_{prefix}_2.fq.gz",
        "--action=trim",
        "--pair-filter=both",
        *input_files
    ]

    if extra_args:
        cutadapt_command += extra_args, 
    subprocess.run(cutadapt_command, check=True)

def seqkit_sanitize(input_file, output_dir):
    """Sanitize FASTQ files using SeqKit."""
    
    # Extract the base filename and append .sanitised.fq.gz

    basename = os.path.basename(input_file)
    extensions = ".fq.gz"
    #prefix_match = re.findall(r"^(.*?)(?:_[Rr]?[12])(?:_\d{3})?\..+$", basename)
    prefix = str(basename).removesuffix(extensions)
    output_file = os.path.join(output_dir, prefix + ".sanitised.fq.gz")

    # Run SeqKit sanitize
    command = ["seqkit", "sana", input_file, "-o", output_file]
    subprocess.run(command, check=True)

    return output_file  # Return sanitized filename if needed

def find_files(output_dir):
    """Find and pair sanitized FASTQ files based on common prefixes."""
    results = {}
    file_dict = {}

    # Regex to capture the common prefix before _1.fastq or _2.fastq
    file_pattern = re.compile(r"(.+?)[._](R?[12])\.sanitised\.fq.gz$")

    # Group files by their prefix
    for file in glob.glob(f"{output_dir}/*.sanitised.fq.gz"):
        match = file_pattern.search(os.path.basename(file))
        if match:
            prefix, read_pair = match.groups()
            file_dict.setdefault(prefix, {})[read_pair] = file

    # Store only complete (_1, _2) pairs
    for prefix, files in file_dict.items():
        if "1" in files and "2" in files:
            results[prefix] = (files["1"], files["2"])

    return results

def seqkit_pair(r1_path, r2_path, output_dir):
    """Pair reads using SeqKit."""
    command = ["seqkit", "pair", "-1", r1_path, "-2", r2_path, "-O", output_dir]
    subprocess.run(command, check=True)

def generate_seqkit_stats(output_dir, stats_output):
    """Generate statistics with SeqKit."""
    fastq_files = glob.glob(os.path.join(output_dir, "*.f*q.gz"))
    seqkit_command = ["seqkit", "stats", *fastq_files]
    with open(stats_output, "w") as stats_file:
        subprocess.run(seqkit_command, stdout=stats_file, check=True)


def main(cutadapt_error_rate, i7_barcodes, i5_barcodes, input_dir, output_dir, extra_cutadapt_args=None):
    os.makedirs(output_dir, exist_ok=True)
    stats_output = os.path.join(output_dir, "cutadapt_summary.stats")
    
    os.makedirs(os.path.join(output_dir, "raw_demultiplexed"), exist_ok=True)
    demux_dir = os.path.join(output_dir, "raw_demultiplexed")

    os.makedirs(os.path.join(output_dir, "tmp"), exist_ok=True)
    tmp_dir = os.path.join(output_dir, "tmp")

    #Step0: Check to make sure all input files are gzipped
    logging.info("Checking to make sure all input files in the input directory are gzipped")     
    check_and_handle_gunzipped(input_dir)

    # Step 1: Find paired input files
    logging.info("Pairing input files")     
    paired_files = pair_input_files(input_dir, "Undetermined")

    # Step 2: Run Cutadapt in parallel
    logging.info("Starting Cutadapt processing...")
    with ThreadPoolExecutor() as executor:
        futures = [
            executor.submit(run_cutadapt, pair, demux_dir, cutadapt_error_rate, i7_barcodes, i5_barcodes, extra_cutadapt_args)
            for pair in paired_files
        ]
        for future in as_completed(futures):
            future.result()
    logging.info("Cutadapt processing complete.")

    # Step 3: Sanitize reads in parallel
    for f in glob.glob(os.path.join(demux_dir,"unknown*.f*q.gz")):
        os.remove(f)

    logging.info("Sanitizing files.")
    fastq_files = glob.glob(os.path.join(demux_dir, "*.f*q.gz"))

    with ThreadPoolExecutor() as executor:
        futures = {executor.submit(seqkit_sanitize, file, tmp_dir): file for file in fastq_files}
        for future in as_completed(futures):
            future.result()

    # Step 4: Find paired sanitized files
    logging.info("Pairing sanitised files.")
    paired_files = find_files(tmp_dir)

    # Step 5: Pair sanitised files in parallel
    with ThreadPoolExecutor() as executor:
        for r1_path, r2_path in paired_files.values():
            executor.submit(seqkit_pair, r1_path, r2_path, output_dir)

    paired_files = glob.glob(os.path.join(output_dir, "*.f*q.gz"))
    for f in paired_files:
        new_name = re.sub(r'\.sanitised', '', f)
        os.rename(f, new_name)

    # Step 6: Generate SeqKit statistics
    logging.info("Generating summary stats.")
    generate_seqkit_stats(output_dir, stats_output)

    # Step 7: Check and handle gunzipped files
    check_and_handle_gunzipped(output_dir)

    # Step 8: Remove tmp dir and clean up file names
    clean_and_tar(output_dir, input_dir)

    shutil.rmtree(demux_dir)


if __name__ == "__main__":
    if len(os.sys.argv) == 1:
        print("No arguments provided. Showing help:")
        parser = argparse.ArgumentParser(description="Concatenate paired reads and compress input directories.")
        parser.print_help()
        os.sys.exit(1)

    parser = argparse.ArgumentParser(description="Demultiplex, sanitize, and summarize FASTQ files.")
    parser.add_argument("--cutadapt_error_rate", type=float, required=True, help="Cutadapt error rate (e.g., 0.2)")
    parser.add_argument("--i7_barcodes", required=True, help="Path to i7 barcodes fasta file")
    parser.add_argument("--i5_barcodes", required=True, help="Path to i5 barcodes fasta file")
    parser.add_argument("--input_dir", required=True, help="Input directory with raw FASTQ files")
    parser.add_argument("--output_dir", required=True, help="Output directory for demultiplexed results")
    parser.add_argument("--cutadapt_extra_args", nargs=argparse.REMAINDER, help="Extra arguments (e.g., --minimum-length 30) to pass to cutadapt.")

    args = parser.parse_args()

    main(
        cutadapt_error_rate=args.cutadapt_error_rate,
        i7_barcodes=args.i7_barcodes,
        i5_barcodes=args.i5_barcodes,
        input_dir=args.input_dir,
        output_dir=args.output_dir,
        extra_cutadapt_args=args.cutadapt_extra_args
    )

# Notes on cutadapt_error_rate: cutadapt variable -e = rounded down % of adapter. e.g, 0.1 = 10%; 0.2 = 20%; For example, an adapter match of length 8 containing 1 error has an error rate of 1/8=0.125. At the default maximum error rate 0.1, it would not be found, but a match of length
