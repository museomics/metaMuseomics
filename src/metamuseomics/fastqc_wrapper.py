import os
import sys
import csv
import subprocess
import logging
import pathlib
import pandas as pd
from concurrent.futures import ThreadPoolExecutor
from seqpy_tools import pair_input_files, xlsx2csv

# Set up logging
log_dir = "./logs"
os.makedirs(log_dir, exist_ok=True)
logger = logging.getLogger()
logger.addHandler(logging.StreamHandler(sys.stdout))
logger.addHandler(logging.FileHandler(os.path.join(log_dir, "fastqc_raw.log")))
logger.setLevel(logging.INFO)


# Function to run FastQC
def run_fastqc(r1_path, r2_path, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    command = ["fastqc", "-o", output_dir, r1_path, r2_path]
    result = subprocess.run(command, capture_output=True, text=True)
    print(f"Processed {r1_path} and {r2_path}")


# Main logic function
def main(input_dir, tracking_sheet=None, column_name=None, prefix=None, sheet=None, output_dir="fastqc_output", max_workers=4):
    # Get sample prefixes from tracking sheet, prefix, or all files
    if args.tracking_sheet:
        if args.column_name:
            logger.info(f"Getting prefixes from column '{args.column_name}' in '{args.tracking_sheet}'")
            filetype = pathlib.Path(args.tracking_sheet).suffix.lower()
            if filetype == ".xlsx":
                df = xlsx2csv(args.tracking_sheet, sheet=args.sheet if hasattr(args, 'sheet') else None)
            elif filetype == ".csv":
                df = pd.read_csv(args.tracking_sheet)
            else:
                logger.error(f"Unsupported tracking sheet format: {args.tracking_sheet}")
                sys.exit(1)

            prefix_list = df[args.column_name].dropna().astype(str).tolist()
            pairs = pair_input_files(args.input_dir, prefix_list)
        else:
            logger.error("No column name specified. Use '--column_name' with '--tracking_sheet'")
            sys.exit(1)

    elif args.prefix:
        logger.info(f"Finding all files beginning with prefix: {args.prefix}")
        pairs = pair_input_files(args.input_dir, args.prefix)

    elif args.suffix:
        logger.info(f"Finding all files ending with suffix: {args.suffix}")
        pairs = pair_input_files(args.input_dir, args.suffix)

    else:
        logger.info(f"No prefix or tracking sheet specified. Using all files in {args.input_dir}")
        pairs = pair_input_files(args.input_dir)

        # Filter out pairs where either file starts with 'Undetermined'
        filtered_pairs = [
            pair for pair in pairs
            if not any(os.path.basename(f).startswith("Undetermined") for f in pair)
        ]
        logger.info(f"Filtered out 'Undetermined' files. {len(filtered_pairs)} valid pairs remain.")
        pairs = filtered_pairs

    # Run FastQC in parallel
    futures = []
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        for file_pair in pairs:
            r1_path, r2_path = file_pair
            ids = pathlib.Path(r1_path).stem.split("_")[0]  # extract sample id
            futures.append(executor.submit(run_fastqc, r1_path, r2_path, output_dir))

        for future in futures:
            future.result()


# Example usage in __main__
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Run FastQC on paired-end reads.")
    parser.add_argument("--input_dir", required=True, help="Directory containing input FASTQ files")
    parser.add_argument("--tracking_sheet", help="CSV or XLSX file with sample prefixes")
    parser.add_argument("--column_name", help="Column in tracking sheet to extract prefixes from")
    parser.add_argument("--sheet", help="(for XLSX input only) The name of the sheet in XLSX file. Default will take the first sheet.")
    parser.add_argument("--prefix", nargs="+", help="Sample prefix(es) to match files if not providing a sample sheet.")
    parser.add_argument("--suffix", help="Suffix pattern between sample ID and 1/2.f*q.* (e.g., 'unmerged' in 'SampleA_unmerged_1.fq.gz'). Can be used in conjunction with --prefix or on its own.")
    parser.add_argument("--output_dir", default="fastqc_output", help="Directory to store FastQC results")
    parser.add_argument("--max_workers", type=int, default=4, help="Maximum number of threads")

    args = parser.parse_args()

    main(
        input_dir=args.input_dir,
        tracking_sheet=args.tracking_sheet,
        column_name=args.column_name,
        prefix=args.prefix,
        sheet=args.sheet,
        output_dir=args.output_dir,
        max_workers=args.max_workers
    )

##Example usage: python3.12 fastqc_wrapper.py --input_dir ./raw_data
