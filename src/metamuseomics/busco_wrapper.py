import subprocess
import logging
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import argparse
import json 

import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


def assemblies_exist_for_all_samples(assembly_dir, assembler):
    assembly_dir = Path(assembly_dir)

    if assembler == "megahit":
        target_file = "final.contigs.fa"
    elif assembler == "metaspades":
        target_file = "scaffolds.fasta"
    elif assembler == "idba_ud":
        target_file = "final.contig.fa"
    else:
        logger.error("Unsupported assembler: %s", assembler)
        return []

    contig_files = []

    for subdir in sorted(assembly_dir.iterdir()):
        if not subdir.is_dir():
            continue

        contigs_path = subdir / target_file

        if not contigs_path.is_file():
            logger.warning(
                "Skipping %s - contig file missing (%s)",
                subdir.name,
                contigs_path
            )
            continue

        contig_files.append(contigs_path)

    if not contig_files:
        logger.warning("No valid assemblies found in %s", assembly_dir)

    return contig_files


def busco_completed(sample_name, busco_outdir):
    """
    Check whether BUSCO has already completed for a sample.
    """
    sample_dir = busco_outdir / sample_name

    if not sample_dir.exists():
        return False

    summary_files = list(sample_dir.glob("short_summary*.txt"))

    return len(summary_files) > 0


def run_busco(contigs_path, busco_outdir, lineage, mode="genome", threads=8):

    sample_name = contigs_path.parent.name
    outdir = busco_outdir / sample_name
    outdir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "busco",
        "-i", str(contigs_path),
        "-l", lineage,
        "-m", mode,
        "-c", str(threads),
        "-o", str(outdir),
        "-f",
        "--metaeuk"
    ]

    logger.info("Running BUSCO for %s", sample_name)

    subprocess.run(cmd, check=True)

    logger.info("Finished BUSCO for %s", sample_name)

    return sample_name


def ensure_busco_lineage(lineage):

    try:
        subprocess.run(
            ["busco", "--download", lineage],
            check=True
        )
        logger.info("BUSCO lineage %s is available", lineage)

    except subprocess.CalledProcessError as e:
        logger.error("Failed to download BUSCO lineage %s", lineage)
        raise RuntimeError(
            f"Unable to download BUSCO lineage dataset: {lineage}"
        ) from e


def collect_busco_summary(busco_outdir):
    """
    Collect the 'results' section from all BUSCO JSON files
    and combine them into a single CSV file.

    Every sample directory in busco_outdir is included (NA if BUSCO has failed)
    """

    busco_outdir = Path(busco_outdir)

    sample_dirs = sorted(
        d for d in busco_outdir.iterdir()
        if d.is_dir()
    )

    data = []
    samples_with_results = set()

    # Search recursively for all JSON files
    for json_file in busco_outdir.rglob("*.json"):
        try:
            with open(json_file, "r") as f:
                content = json.load(f)

            results = content.get("results", {})

            if results:
                # Make a copy so the original JSON data is not modified
                results = results.copy()

                # Add filename and sample as a column
                results["filename"] = json_file.name
                results["sample"] = json_file.parent.name

                sample_name = json_file.parent.name

                data.append(results)
                samples_with_results.add(sample_name)

        except (json.JSONDecodeError, OSError) as e:
            logger.warning(
                "Could not read BUSCO JSON file %s: %s",
                json_file,
                e
            )

    if data:
        df = pd.DataFrame(data)
    else:
        df = pd.DataFrame()


    # Add samples with no JSON results
    for sample_dir in sample_dirs:

        sample_name = sample_dir.name
        if sample_name not in samples_with_results:
            logger.info(
                "No BUSCO JSON results found for %s - adding NA row",
                sample_name
            )

            # Create a row containing NA for every existing column
            missing_row = {
                column: pd.NA
                for column in df.columns
            }

            # Always retain the sample name
            missing_row["sample"] = sample_name

            # Filename is NA because no JSON exists
            missing_row["filename"] = pd.NA

            # Add the row
            df = pd.concat(
                [df, pd.DataFrame([missing_row])],
                ignore_index=True
            )


    # Save summary
    if df.empty:
        logger.warning(
            "No BUSCO sample directories or JSON results found in %s",
            busco_outdir
        )
        return None

    output_csv = busco_outdir / "busco_summary.csv"
    df.to_csv(output_csv, index=False)

    logger.info(
        "BUSCO summary written to %s (%d samples)",
        output_csv,
        len(df)
    )

    print("\nBUSCO summary:")
    print(df.head())

    return df


def main(
    assembly_dir,
    assembler,
    busco_outdir,
    lineage,
    threads,
    max_workers,
    mode="genome"
):

    assembly_dir = Path(assembly_dir)
    busco_outdir = Path(busco_outdir)
    busco_outdir.mkdir(parents=True, exist_ok=True)

    ensure_busco_lineage(lineage)

    contig_files = assemblies_exist_for_all_samples(
        assembly_dir,
        assembler
    )

    if not contig_files:
        logger.error("No valid contig files found. BUSCO will not be run.")
        return

    ## Resume if previously broken

    remaining_contigs = []

    for contigs in contig_files:
        sample = contigs.parent.name

        if busco_completed(sample, busco_outdir):
            logger.info("Skipping %s - BUSCO already completed", sample)
        else:
            remaining_contigs.append(contigs)

    if remaining_contigs:

        logger.info(
            "Running BUSCO on %d remaining samples using %d workers",
            len(remaining_contigs),
            max_workers
        )

        with ProcessPoolExecutor(max_workers=max_workers) as executor:

            futures = {
                executor.submit(
                    run_busco,
                    contigs_path,
                    busco_outdir,
                    lineage,
                    mode,
                    threads
                ): contigs_path
                for contigs_path in remaining_contigs
            }

            for future in as_completed(futures):
                contigs_path = futures[future]
                try:
                    sample = future.result()
                    logger.info("BUSCO completed successfully for %s", sample)

                except Exception as e:
                    logger.error(
                        "BUSCO failed for %s: %s",
                        contigs_path.parent.name, e
                    )

    else:
        logger.info("All BUSCO runs already completed.")

    # Collect BUSCO JSON results after jobs have finished.
    logger.info("Collecting BUSCO results.")
    collect_busco_summary(busco_outdir=busco_outdir)

    logger.info("Job done.")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Check assemblies and run BUSCO in parallel"
    )

    parser.add_argument(
        "--assembly-dir",
        required=True
    )

    parser.add_argument(
        "--assembler",
        required=True,
        choices=["megahit", "metaspades", "idba_ud"]
    )

    parser.add_argument(
        "--busco_outdir",
        default="busco_out"
    )

    parser.add_argument(
        "--lineage",
        required=True
    )

    parser.add_argument(
        "--threads",
        type=int,
        default=8
    )

    parser.add_argument(
        "--max-workers",
        type=int,
        default=4
    )

    parser.add_argument(
        "--mode",
        default="genome",
        choices=["genome", "proteins", "transcriptome"]
    )

    args = parser.parse_args()

    main(
        assembly_dir=args.assembly_dir,
        assembler=args.assembler,
        busco_outdir=args.busco_outdir,
        lineage=args.lineage,
        threads=args.threads,
        max_workers=args.max_workers,
        mode=args.mode
    )
