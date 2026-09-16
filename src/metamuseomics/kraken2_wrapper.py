import argparse
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from collections import defaultdict
from typing import Optional


# -----------------------------
# CONFIG
# -----------------------------

TAXDUMP = Path("taxdump")
THREADS = 16


# UTILS

def run(cmd):
    print(" ".join(cmd))
    subprocess.run(cmd, check=True)


# STEP 1: RUN KRAKEN2

def run_kraken(reads_fasta, out_tsv, out_report, kraken_db):

    out_tsv.parent.mkdir(parents=True, exist_ok=True)

    run([
        "kraken2",
        "--db", kraken_db,
        "--threads", str(THREADS),
        "--use-names",
        "--output", str(out_tsv),
        "--report", str(out_report),
        "--memory-mapping",
        str(reads_fasta)
    ])


# TAXONOMY HELPER FUNCTIONS

def load_nodes(nodes_file):
    parent = {}
    rank = {}

    with nodes_file.open() as fh:
        for line in fh:
            fields = [f.strip() for f in line.split("|")]
            taxid = fields[0]
            parent[taxid] = fields[1]
            rank[taxid] = fields[2]

    return parent, rank


def load_names(names_file):
    names = {}

    with names_file.open() as fh:
        for line in fh:
            fields = [f.strip() for f in line.split("|")]
            if fields[3] == "scientific name":
                names[fields[0]] = fields[1]

    return names


# STEP 3: TAXID TO FAMILY

def taxid_to_family(taxid, parent, rank, names):
    """
    Walk up taxonomy until family or root
    """
    while taxid != "1" and taxid in parent:
        if rank.get(taxid) == "family":
            return names.get(taxid, "unidentified_to_family")
        taxid = parent[taxid]

    return "unidentified_to_family"


# STEP 4: PARSE KRAKEN OUTPUT

def parse_kraken(
    kraken_tsv: Path,
    parent: dict,
    rank: dict,
    names: dict
):
    families = defaultdict(list)

    with kraken_tsv.open() as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            status, read_id, taxid = fields[:3]

            if status != "C":
                family = "unidentified_to_family"
            else:
                family = taxid_to_family(taxid, parent, rank, names)

            families[family].append(read_id)

    return families


# STEP 5: WRITE FAMILY FASTAS

def load_fasta(fasta: Path):
    seqs = {}
    cur = None
    buff = []

    with fasta.open() as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if cur:
                    seqs[cur] = "".join(buff)
                cur = line[1:].split()[0]
                buff = []
            else:
                buff.append(line)
        if cur:
            seqs[cur] = "".join(buff)

    return seqs


def write_family_fastas(
    reads_fasta: Path,
    family_map: dict,
    out_dir: Path
):
    out_dir.mkdir(parents=True, exist_ok=True)

    seqs = load_fasta(reads_fasta)

    for family, reads in family_map.items():
        out_fa = out_dir / f"{family}.fasta"

        with out_fa.open("w") as fh:
            for rid in reads:
                if rid in seqs:
                    fh.write(f">{rid}\n{seqs[rid]}\n")


# PER-SAMPLE PIPELINE

def process_sample(reads_fasta: Path, sample_name: str, base_out_dir: Path, kraken_db: str):
    """
    Run the full pipeline for a single sample. Outputs are written to
    base_out_dir/<sample_name>/.
    """
    print(f"\n{'='*60}")
    print(f"Processing sample: {sample_name}")
    print(f"Input: {reads_fasta}")
    print(f"{'='*60}\n")

    sample_out  = base_out_dir / sample_name
    kraken_out  = sample_out / "kraken" / "reads.kraken"
    kraken_report = sample_out / "kraken" / "reads.report"
    family_dir  = sample_out / "families"

    run_kraken(reads_fasta, kraken_out, kraken_report, kraken_db)

    #parent, rank = load_nodes(TAXDUMP / "nodes.dmp")
    #names = load_names(TAXDUMP / "names.dmp")

    #family_map = parse_kraken(kraken_out, parent, rank, names)
    #write_family_fastas(reads_fasta, family_map, family_dir)

    print(f"[DONE] {sample_name}")
    return sample_name


# SAMPLE DISCOVERY

FASTA_EXTENSIONS = {".fa", ".fasta", ".fna"}

def discover_samples(input_dir: Path, filename: Optional[str]) -> list:
    """
    Returns a list of (sample_name, fasta_path) tuples.

    If --filename is given, search for input_dir/*/<filename>.
    Otherwise, collect FASTA files directly in input_dir or one level of
    subdirectories (only when a subdir contains exactly one FASTA file).
    """
    samples = []

    if filename:
        for subdir in sorted(input_dir.iterdir()):
            if not subdir.is_dir():
                continue
            candidate = subdir / filename
            if candidate.is_file():
                samples.append((subdir.name, candidate))
            else:
                print(f"[WARN] {filename} not found in {subdir}, skipping.", file=sys.stderr)
    else:
        # Direct FASTA files in input_dir
        direct = [
            (f.stem, f)
            for f in sorted(input_dir.iterdir())
            if f.is_file() and f.suffix.lower() in FASTA_EXTENSIONS
        ]
        # One FASTA per subdirectory
        subdir_fastas = []
        for subdir in sorted(input_dir.iterdir()):
            if not subdir.is_dir():
                continue
            found = [
                f for f in sorted(subdir.iterdir())
                if f.is_file() and f.suffix.lower() in FASTA_EXTENSIONS
            ]
            if len(found) == 1:
                subdir_fastas.append((subdir.name, found[0]))
            elif len(found) > 1:
                print(
                    f"[WARN] Multiple FASTA files found in {subdir}; "
                    "use --filename to specify one. Skipping.",
                    file=sys.stderr
                )

        samples = direct + subdir_fastas

    if not samples:
        sys.exit(f"[ERROR] No samples found in {input_dir}. Check --input_dir / --filename.")

    return samples


# MAIN

def parse_args():
    parser = argparse.ArgumentParser(
        description="Run Kraken2 LCA pipeline across multiple samples in parallel."
    )
    parser.add_argument(
        "--input_dir", required=True, type=Path,
        help="Directory containing sample FASTA files or subdirectories."
    )
    parser.add_argument(
        "--filename", default=None,
        help=(
            "If samples share a common filename inside subdirectories "
            "(e.g. final.contigs.fa), specify it here. "
            "Each subdirectory of --input_dir that contains this file becomes one sample."
        )
    )
    parser.add_argument(
        "--kraken_db", required=True,
        help="Path to the Kraken2 database directory."
    )
    parser.add_argument(
        "--output_dir", default="output", type=Path,
        help="Base output directory (default: output/)."
    )
    parser.add_argument(
        "--workers", type=int, default=4,
        help="Number of samples to process in parallel (default: 4)."
    )
    return parser.parse_args()


def main():
    args = parse_args()

    if not args.input_dir.is_dir():
        sys.exit(f"[ERROR] --input_dir {args.input_dir} does not exist or is not a directory.")

    if not Path(args.kraken_db).is_dir():
        sys.exit(f"[ERROR] --kraken_db {args.kraken_db} does not exist or is not a directory.")

    samples = discover_samples(args.input_dir, args.filename)

    print(f"Found {len(samples)} sample(s):")
    for name, path in samples:
        print(f"  {name}: {path}")
    print()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    failed = []

    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        futures = {
            executor.submit(process_sample, fasta, name, args.output_dir, args.kraken_db): name
            for name, fasta in samples
        }

        for future in as_completed(futures):
            sample_name = futures[future]
            try:
                future.result()
            except Exception as exc:
                print(f"[ERROR] Sample {sample_name} failed: {exc}", file=sys.stderr)
                failed.append(sample_name)

    if failed:
        print(f"\n[WARN] {len(failed)} sample(s) failed: {', '.join(failed)}", file=sys.stderr)
        sys.exit(1)
    else:
        print(f"\nAll {len(samples)} sample(s) completed successfully.")


if __name__ == "__main__":
    main()
