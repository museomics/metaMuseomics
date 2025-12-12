#!/usr/bin/env python3
import os
import subprocess
import argparse
import glob
import csv

def run_nonpareil(input_path, output_dir, threads=4, reads_to_sample=100000):
    ''' Run Nonpareil on a single FASTQ file.'''
    sample_name = os.path.splitext(os.path.basename(input_path))[0]
    output_prefix = os.path.join(output_dir, sample_name)

    cmd = [
        "nonpareil",
        "-s", input_path,
        "-T", "kmer",
        "-f", "fastq",
        "-b", output_prefix,
        "-R", str(reads_to_sample),
        "-t", str(threads)
    ]

    print(f"\n Running Nonpareil on {sample_name}...")
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running Nonpareil on {sample_name}: {e}")
        return None

    return output_prefix + ".npo"

def parse_npo(npo_file):
    ''' Parse nonpareil .npo file to extract metrics.'''

    metrics = {
        'sample': os.path.basename(npo_file).replace('.npo', ''),
        'coverage': None,
        'redundancy': None,
        'diversity': None,
        'effort95Gbp': None
    }
    if not os.path.exists(npo_file):
        print(f"Missing .npo file: {npo_file}")
        return metrics

    with open(npo_file, encoding="utf-8") as f:
        for line in f:
            if line.startswith("Coverage"):
                metrics['coverage'] = float(line.strip().split(":")[1])
            elif line.startswith("Redundancy"):
                metrics['redundancy'] = float(line.strip().split(":")[1])
            elif line.startswith("Diversity"):
                metrics['diversity'] = float(line.strip().split(":")[1])
            elif "Effort for 95% coverage" in line:
                effort = line.strip().split(":")[1].strip().split()[0]
                metrics['effort95Gbp'] = float(effort)

    return metrics

def batch_run_nonpareil(fastq_dir, output_dir, threads=4, reads_to_sample=100000):
    ''' Run Nonpareil on all FASTQ files in a directory and summarize results.'''

    os.makedirs(output_dir, exist_ok=True)
    fastq_files = glob.glob(os.path.join(fastq_dir, "*min*.fq"))

    results = []
    for fq in sorted(fastq_files):
        npo_file = run_nonpareil(fq, output_dir, threads=threads, reads_to_sample=reads_to_sample)
        if npo_file:
            metrics = parse_npo(npo_file)
            results.append(metrics)

    # Write summary
    output_csv = os.path.join(output_dir, "nonpareil_summary.csv")
    with open(output_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=results[0].keys())
        writer.writeheader()
        writer.writerows(results)
    print(f"\n Summary saved to {output_csv}")

def main():
    ''' Main function to parse arguments and run batch Nonpareil.'''

    parser = argparse.ArgumentParser(
        description="Run Nonpareil on multiple FASTQ files for metagenomic complexity estimation."
        )
    parser.add_argument("-i", "--input_dir", required=True, help="Directory with FASTQ(.gz) files")
    parser.add_argument(
        "-o", "--output_dir", required=True, help="Directory to store Nonpareil output"
        )
    parser.add_argument(
        "-t", "--threads", type=int, default=4, help="Number of threads per run"
        )
    parser.add_argument(
        "-R", "--reads_to_sample", type=int, default=100000,
        help="Number of reads to sample per file"
        )
    args = parser.parse_args()

    batch_run_nonpareil(
        args.input_dir,
        args.output_dir,
        threads=args.threads,
        reads_to_sample=args.reads_to_sample
    )

if __name__ == "__main__":
    main()


## Input directory structure:
#data/
#+-- sample1.fq.gz
#+-- sample2.fq
#+-- sample3.fastq.gz
##python3 run_nonpareil_batch.py -i data/ -o nonpareil_out/ -t 8 -R 100000