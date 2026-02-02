import os
import sys
import argparse
import re
from pathlib import Path
import pandas as pd
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from typing import Optional
import json
import csv
import gzip
import shutil


from metahist_tools import (
    find_single_reads,
    find_paired_files2,
    get_read_ids,
    get_read_ids2,
    find_program,
    setup_logging
    )

# Module for metagenomic assembly using MetaSPAdes, IDBA-UD or MEGAHIT.
# Module will run assemblies based on user parameters, check assemblies,
# and generate summary statistics of contigs.

@dataclass
class Params:
    ''' Dataclass to hold assembly parameters.'''
    assembler: str
    merged: bool
    paired: bool
    both: bool
    input_dir: str
    output_dir: str
    prefix: Optional[str] = None
    suffix: Optional[str] = None
    ids: Optional[str] = None
    tracking_sheet: Optional[str] = None
    column_name: Optional[str] = None
    output_suffix: Optional[str] = None
    unmerged_dir: Optional[str] = None
    threads: Optional[int] = 4
    extra_args: Optional[list] = None
    busco_lineage: Optional[str] = "fungi_odb10"
    correction: Optional[bool] = False

def build_params(args) -> Params:
    ''' Function to build parameters from argparse args.'''

    return Params(
        assembler=args.assembler,
        merged=args.merged,
        paired=args.paired,
        both=args.both,
        prefix=args.prefix,
        suffix=args.suffix,
        ids=args.ids,
        tracking_sheet=args.tracking_sheet,
        input_dir=args.input_dir,
        column_name=args.column_name,
        output_dir=args.output_dir,
        output_suffix=args.output_suffix,
        unmerged_dir=args.unmerged_dir,
        threads=args.threads,
        extra_args=args.extra_args,
        correction=args.correction
    )

# Set up logging
logger = setup_logging("./logs", "assembly.log")

##################### Helper functions #####################

def concatenate_fastas(fasta_a, fasta_b, output_fasta):
    """Concatenate two FASTA files in order."""
    with open(output_fasta, "w") as out:
        for fa in (fasta_a, fasta_b):
            if fa and Path(fa).exists():
                with open(fa) as fh:
                    out.write(fh.read())

def check_dependencies(params):
    """
    Check required external programs based on selected pipeline options.
    Fails fast if anything is missing.
    """

    required = set()

    # Assembler-specific requirements
    if params.assembler == "megahit":
        required.add("megahit")
    elif params.assembler == "metaspades":
        required.add("metaspades.py")
    elif params.assembler == "idba_ud":
        required.update({"idba_ud", "fq2fa"})

    logger.info("Assembler specified: %s", params.assembler)


    # metaMIC correction requirements
    if params.correction:
        required.update({
            "metaMIC",
            "bwa",
            "samtools",
            "seqkit",
        })
        logger.info("Correction selected, checking dependencies...")


    # BUSCO (always run downstream)
    required.add("busco")

    for prog in sorted(required):
        try:
            find_program(prog)
        except Exception:
            logger.error("Missing required program: %s", prog)
            raise RuntimeError(
                f"Required program not found in PATH: {prog}"
            )

##################### Assembly-specific functions #####################

def run_idba_ud(
    sample_id,
    fasta_path,
    output_dir,
    output_dir_suffix=None,
    extra_args=None,
):
    """Run IDBA-UD assembly."""

    logger.info("Starting IDBA-UD for %s", sample_id)

    assembly_dir = Path(output_dir) / f"{sample_id}_{output_dir_suffix}"

    cmd = [
        "idba_ud",
        "-r", str(fasta_path),
        "--num_threads", "1",
        "-o", str(assembly_dir),
    ]

    if extra_args:
        cmd += extra_args

    logger.debug("Running command: %s", " ".join(cmd))

    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        logger.info("IDBA-UD completed for %s", sample_id)
    except subprocess.CalledProcessError as e:
        logger.error("IDBA-UD failed for %s: %s", sample_id, e.stderr)

def process_idba_samples(params):
    """Run IDBA-UD assemblies."""

    params.output_dir.mkdir(exist_ok=True)
    fq2fa = find_program("fq2fa")

    # Resolve sample prefixes
    if params.prefix:
        prefixes = params.prefix
    elif params.tracking_sheet:
        prefixes = get_read_ids(
            params.tracking_sheet,
            mode="single",
            prefix=None,
            suffix=params.suffix,
            column_name=params.column_name,
            sheet=params.sheet,
        )
    else:
        prefixes = get_read_ids(
            params.input_dir,
            mode="single",
            prefix=None,
            suffix=params.suffix,
        )

    jobs = []

    for sample in prefixes:
        fasta_path = None

        if params.merged:
            merged = find_single_reads(params.input_dir, sample)
            if not merged:
                logger.error("[%s] No merged reads found", sample)
                continue
            merged = merged[0]

            if merged.suffix in {".fq", ".fastq"}:
                fasta_path = merged.with_suffix(".fa")
                subprocess.run([fq2fa, merged, fasta_path], check=True)
            else:
                fasta_path = merged

        elif params.paired:
            r1, r2 = find_paired_files2(params.input_dir, sample)
            if not (r1 and r2):
                logger.error("[%s] Missing paired reads", sample)
                continue

            fasta_path = params.input_dir / f"{sample}_interleaved.fa"
            subprocess.run(
                [fq2fa, "--merge", "--filter", r1, r2, fasta_path],
                check=True,
            )

        jobs.append(
            (
                sample,
                fasta_path,
                params.output_dir,
                params.output_suffix,
                params.extra_args,
            )
        )

    # Parallel execution
    with ProcessPoolExecutor(max_workers=params.threads) as exe:
        futures = [exe.submit(run_idba_ud, *job) for job in jobs]
        for fut in as_completed(futures):
            fut.result()

    check_assemblies(
        params.output_dir,
        assembler="idba_ud",
        params=params,
    )

def run_metaspades(
    sample_id,
    merged_path=None,
    r1_path=None,
    r2_path=None,
    output_dir=None,
    output_dir_suffix=None,
    extra_args=None,
):
    """Run MetaSPAdes."""

    logger.info("Starting MetaSPAdes for %s", sample_id)

    assembly_dir = Path(output_dir) / f"{sample_id}_{output_dir_suffix}"

    cmd = [
        "metaspades.py",
        "--phred-offset", "33",
        "-o", str(assembly_dir),
    ]

    if extra_args:
        cmd += extra_args

    if merged_path and r1_path and r2_path:
        cmd += ["--merged", merged_path, "-1", r1_path, "-2", r2_path]
    elif r1_path and r2_path:
        cmd += ["-1", r1_path, "-2", r2_path]
    else:
        logger.error("[%s] Missing reads for MetaSPAdes", sample_id)
        return

    logger.debug("Running command: %s", " ".join(cmd))

    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        logger.info("MetaSPAdes completed for %s", sample_id)
    except subprocess.CalledProcessError as e:
        logger.error("MetaSPAdes failed for %s: %s", sample_id, e.stderr)

def process_metaspades_samples(params):
    """Run MetaSPAdes assemblies."""

    params.output_dir.mkdir(exist_ok=True)

    if params.prefix:
        prefixes = params.prefix
    elif params.tracking_sheet:
        prefixes = get_read_ids(
            params.tracking_sheet,
            mode="single",
            prefix=None,
            suffix=params.suffix,
            column_name=params.column_name,
            sheet=params.sheet,
        )
    else:
        prefixes = get_read_ids(
            params.input_dir,
            mode="single",
            prefix=None,
            suffix=params.suffix,
        )

    jobs = []

    for sample in prefixes:
        merged = r1 = r2 = None

        if params.both:
            merged_files = find_single_reads(params.input_dir, sample)
            if not merged_files:
                logger.error("[%s] No merged reads", sample)
                continue
            merged = merged_files[0]

        if params.paired or params.both:
            r1, r2 = find_paired_files2(
                params.unmerged_dir or params.input_dir,
                sample,
            )
            if not (r1 and r2):
                logger.error("[%s] Missing paired reads", sample)
                continue

        jobs.append(
            (
                sample,
                str(merged) if merged else None,
                str(r1) if r1 else None,
                str(r2) if r2 else None,
                params.output_dir,
                params.output_suffix,
                params.extra_args,
            )
        )

    with ProcessPoolExecutor(max_workers=params.threads) as exe:
        futures = [exe.submit(run_metaspades, *job) for job in jobs]
        for fut in as_completed(futures):
            fut.result()

    check_assemblies(
        params.output_dir,
        assembler="metaspades",
        params=params,
    )

def run_megahit(sample_id, merged_path=None, r1_path=None, r2_path=None, output_dir=None, output_dir_suffix=None, extra_args=None):
    '''Function to run megahit.'''
         
    logger.info("Starting Megahit for %s.", sample_id)

    assembly_dir = os.path.join(output_dir, f"{sample_id}_{output_dir_suffix}")

    # Form the command based on the available files
    command = ["megahit", "-o", assembly_dir]

    if extra_args:
        command += extra_args
    if merged_path and not (r1_path or r2_path):
        command.extend(["-r", merged_path])
    elif not merged_path and r1_path and r2_path:
        command.extend(["-1", r1_path, "-2", r2_path])
    elif merged_path and r1_path and r2_path:
        command.extend(["-r", merged_path, "-1", r1_path, "-2", r2_path])
    else:
        logger.error(
            "ID %s: Missing necessary files to run Megahit. Skipping...",
            sample_id
        )
        return

    logger.debug("Running command: %s", " ".join(command))

    try:
        subprocess.run(command, capture_output=True, text=True, check=True)
        logger.info("Megahit completed for %s.", sample_id)
    except subprocess.CalledProcessError as e:
        logger.error("Megahit failed for %s : %s.", sample_id, e.stderr)
        return

def process_megahit_samples(params):
    ''' 
    Function that runs megahit assembly on samples based on user arguments,
    checks assemblies and generated summary statistics of contigs.
    '''

    assembly_dir = params.output_dir
    unmerged_dir = params.unmerged_dir if params.unmerged_dir else params.input_dir

    if os.path.exists(assembly_dir):
        if assemblies_exist_for_all_samples(assembly_dir, assembler= params.assembler):
            logger.info(
                "All %s assemblies already exist. Skipping assembly step.", params.assembler
            )
            check_assemblies(
                assembly_dir,
                assembler=params.assembler,
                params=params
            )
            return

    if os.path.exists(params.output_dir):
        pass
    else:
        os.mkdir(params.output_dir)

    if params.merged:
        if params.prefix:
            prefixes = [params.prefix]
        else:
            if params.tracking_sheet:
                prefixes = get_read_ids(params.tracking_sheet,
                    mode="single",
                    prefix=params.prefix,
                    suffix=params.suffix,
                    column_name=params.column_name
                )
            else:
                prefixes = get_read_ids(params.input_dir,
                    mode="single",
                    prefix=params.prefix,
                    suffix=params.suffix
                )
                print(prefixes, "\n")
    elif params.paired:
        ### uses ids to find pairs in the input dir
        if params.prefix:
            prefixes = [params.prefix]
        else:
            if params.tracking_sheet:
                prefixes = get_read_ids2(params.tracking_sheet,
                    paired=True,
                    prefix=params.prefix,
                    suffix=params.suffix,
                    column_name=params.column_name,
                    sheet=params.tracking_sheet, logger=None
                )
            elif params.ids:
                prefixes = params.ids.split(",")
            else:
                prefixes = get_read_ids2(params.input_dir,
                    paired=True,
                    prefix=params.prefix,
                    suffix=params.suffix,
                    column_name=params.column_name,
                    sheet=params.tracking_sheet, logger=None
                )
        paired_files = {}
    else:
        logger.error("Either --paired or --merged must be specified.")
        return

    # Ensure paired_files exists for all modes (paired, both, merged)
    paired_files = {}
    logger.info("Found files for %s", prefixes)

    unmerged_dir = params.unmerged_dir if params.unmerged_dir else params.input_dir
    logger.info(unmerged_dir)
    extra_args = params.extra_args if hasattr(params, "extra_args") else None

    jobs = []
    for prefix in prefixes:
        logger.info("Starting Megahit for %s", prefix)

        # Merged reads
        merged_file_path = None
        if params.merged or params.both:
            merged_files = find_single_reads(params.input_dir, prefix=prefix)
            merged_file_path = merged_files[0]
            logger.info("Using %s for prefix: %s", merged_file_path, prefix)
            if not merged_file_path and (params.merged or params.both):
                raise ValueError(f"No merged file found for prefix: {prefix}")

        # Add job parameters based on mode
        if params.paired:
            r1_path, r2_path = find_paired_files2(params.input_dir, prefix)
            if r1_path and r2_path:
                paired_files[prefix] = [r1_path, r2_path]
            paired_set = paired_files.get(prefix, (None, None))
            jobs.append((prefix, None, str(paired_set[0]), str(paired_set[1]),
                        params.output_dir, params.output_suffix, extra_args))
        elif params.both:
            r1_path, r2_path = find_paired_files2(params.input_dir, prefix)
            if r1_path and r2_path:
                paired_files[prefix] = [r1_path, r2_path]
            paired_set = paired_files.get(prefix, (None, None))
            jobs.append((prefix, str(merged_file_path), str(paired_set[0]),
                        str(paired_set[1]), params.output_dir, params.output_suffix, extra_args))
        elif params.merged:
            jobs.append((prefix, str(merged_file_path), None, None,
                        params.output_dir, params.output_suffix, extra_args))


    # Run jobs in parallel
    with ProcessPoolExecutor(max_workers=params.threads) as executor:
        futures = [executor.submit(run_megahit, *job) for job in jobs]

        for future in as_completed(futures):
            try:
                future.result()
            except subprocess.CalledProcessError as e:
                logger.error("Megahit job failed: %s.", str(e.stderr))

    ## Check assembly
    assembly_dir = params.output_dir
    check_assemblies(assembly_dir, assembler="megahit", params=params)


##################### Assembly checks #####################

def assemblies_exist_for_all_samples(assembly_dir, assembler):
    """Check if assemblies exist for all samples in the assembly directory."""

    assembly_dir = Path(assembly_dir)

    if assembler == "megahit":
        target_file = "final.contigs.fa"
    elif assembler == "metaspades":
        target_file = "scaffolds.fasta"
    elif assembler == "idba_ud":
        target_file = "final.contig.fa"
    else:
        logger.error("Unsupported assembler: %s", assembler)
        return False

    found_any = False

    for subdir in assembly_dir.iterdir():
        if not subdir.is_dir():
            continue

        found_any = True
        contigs_path = subdir / target_file

        if not contigs_path.is_file():
            logger.info(
                "Missing contig file for %s (%s)",
                subdir.name,
                contigs_path
            )
            return False
        else:
            logger.info(
                "A valid contig file exists for %s: %s",
                subdir.name,
                contigs_path
            )

    return found_any

def check_assemblies(
    assembly_dir,
    assembler,
    params
    ):
    """
    Check assemblies
    Returns a list of paths to contig FASTA files to be evaluated downstream
    """

    assembly_dir = Path(assembly_dir)

    if assembler == "megahit":
        target_file = "final.contigs.fa"
        run_func = run_megahit_restart
    elif assembler == "metaspades":
        target_file = "scaffolds.fasta"
        run_func = run_metaspades_restart
    elif assembler == "idba_ud":
        target_file = "final.contig.fa"
        run_func = None
    else:
        logger.error("Unsupported assembler: %s", assembler)
        return []

    # Validate / rerun assemblies

    for subdir in assembly_dir.iterdir():
        if not subdir.is_dir():
            continue

        sample_id = subdir.name
        contigs_path = subdir / target_file

        if contigs_path.is_file():
            logger.info(
                "Valid contig file exists for %s: %s",
                sample_id,
                contigs_path
            )
            continue

        if assembler == "idba_ud":
            scaffold_path = subdir / "scaffold.fa"
            if scaffold_path.is_file():
                logger.info(
                    "Using IDBA-UD scaffold fallback for %s: %s",
                    sample_id,
                    scaffold_path
                )
            else:
                logger.error(
                    "IDBA-UD failed for %s. Consider re-running assembly.",
                    sample_id
                )
            continue

        logger.info(
            "No contigs found for %s, rerunning %s",
            sample_id,
            assembler
        )

        if run_func is not None:
            run_func(subdir, sample_id)

        if not contigs_path.is_file():
            logger.error(
                "Assembler %s failed for sample %s",
                assembler,
                sample_id
            )

    # Collect successful assemblies
    successful_assemblies = []

    for subdir in assembly_dir.iterdir():
        if not subdir.is_dir():
            continue

        contigs_path = subdir / target_file
        if contigs_path.is_file():
            successful_assemblies.append(contigs_path)

    if not successful_assemblies:
        logger.warning("No successful assemblies found")
        return []

    logger.info(
        "Found %d successful assemblies",
        len(successful_assemblies)
    )

    return successful_assemblies

def run_metaspades_restart(output_dir, sample_id):
    '''Rerun metaspades with --continue flag in case failure was due to interruption.'''

    cmd = [
        "metaspades",
        "--continue", 
        "-o", output_dir
    ]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode == 0:
            logger.info("MetaSPAdes succeeded for %s.", sample_id)
        else:
            logger.error("MetaSPAdes failed for %s. Error: %s.", sample_id, result.stderr)
    except OSError as e:
        logger.error("OS error running MetaSPAdes (is it installed?): %s.", e)
    except Exception as e:
        logger.error("MetaSPAdes execution error for %s : %s.", sample_id, e)

def run_megahit_restart(output_dir, sample_id):
    '''Rerun megahit with --continue flag in case failure was due to interruption.'''

    cmd = [
        "megahit",
        "-o", output_dir,
        "--continue"
    ]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode == 0:
            logger.info("MEGAHIT succeeded for %s.", sample_id)
        else:
            logger.error(
                "MEGAHIT failed for %s. Error: %s.",
                sample_id,
                result.stderr
            )
    except OSError as e:
        logger.error("OS error running MEGAHIT (is it installed?): %s.", e)
    except Exception as e:
        logger.error("MEGAHIT execution error for %s: %s.", sample_id, e)


##################### Assembly summaries #####################

def get_coverage_and_correct(contig_file, params):
    """
    Run metaMIC on a contig file.
    Returns: success : bool; corrected_fasta : path to metaMIC corrected_contigs.fa if successful
    """

    contig_file = Path(contig_file)
    sample = contig_file.parent.name.replace("_merged_assembly", "")

    unmerged_dir = (
        Path(params.unmerged_dir)
        if params.unmerged_dir
        else Path(params.input_dir)
    )

    output_dir = contig_file.parent / "metaMIC_correction"
    output_dir.mkdir(parents=True, exist_ok=True)

    bam_file = output_dir / f"{sample}_reads.bam"
    pileup_file = output_dir / f"{sample}_reads.pileup"
    corrected_fasta = output_dir / "corrected_contigs.fa"

    logger.info("[%s] Locating paired reads", sample)

    try:
        if params.suffix is not None:
            r1, r2 = find_paired_files2(unmerged_dir, sample, suffix=params.suffix)
        else:
            r1, r2 = find_paired_files2(unmerged_dir, sample)
    except Exception as e:
        logger.error("[%s] Error finding paired reads: %s", sample, e)
        return False, None

    if not (r1 and r2):
        logger.warning("[%s] No paired reads found", sample)
        return False, None

    try:
        subprocess.run(["bwa", "index", contig_file], check=True)

        subprocess.run(
            f"bwa mem -a -t {params.threads} {contig_file} {r1} {r2} | "
            f"samtools view -h -q 10 -m 50 -F 4 -b | "
            f"samtools sort -o {bam_file}",
            shell=True,
            check=True,
        )

        subprocess.run(
            f"samtools mpileup -C 50 -A -f {contig_file} {bam_file} | "
            f"awk '$3 != \"N\"' > {pileup_file}",
            shell=True,
            check=True,
        )

        subprocess.run(
            [
                "metaMIC", "extract_feature",
                "--bam", bam_file,
                "-c", contig_file,
                "-o", output_dir,
                "--pileup", pileup_file,
                "-m", "meta",
            ],
            check=True,
        )

        subprocess.run(
            [
                "metaMIC", "predict",
                "-c", contig_file,
                "-o", output_dir,
                "-a", params.assembler.upper(),
                "-m", "meta",
            ],
            check=True,
        )

    except subprocess.CalledProcessError as e:
        logger.error("[%s] metaMIC failed: %s", sample, e)
        return False, None

    if not corrected_fasta.exists():
        logger.error("[%s] metaMIC finished but corrected_contigs.fa missing", sample)
        return False, None

    logger.info("[%s] metaMIC correction successful", sample)
    return True, corrected_fasta

def run_metamic_correction(assembler, summary, params):
    """Run metaMIC correction across multiple hDNA assemblies 
        This run metaMIC on each assembly, and if it fails, splits the assembly
        into contigs >=1000 bp and <1000 bp, corrects the >=1000 bp contigs,
        and re-merges them.\n\n
        Also outputs a summary TSV that includes the correction level achieved per sample,
        and output filepath.\n\n"""

    results = []

    if assembler == "megahit":
        target_file = "final.contigs.fa"
    elif assembler == "metaspades":
        target_file = "scaffolds.fasta"
    elif assembler == "idba_ud":
        target_file = "final.contig.fa"
    else:
        logger.error("Unsupported assembler: %s", assembler)
        return

    assembly_paths = sorted(params.input_dir.glob("**/" + target_file))

    for assembly in assembly_paths:
        sample = assembly.parent.name
        logger.info("Processing: %s", sample)

        final_fa = assembly.parent / "assembly.metaMIC_corrected.fa"

        # Attempt full assembly correction
        ok, corrected = get_coverage_and_correct(assembly, params)

        if ok:
            concatenate_fastas(corrected, None, final_fa)
            results.append([sample, "full", str(final_fa)])
            continue

        # Fallback: Split by length and retry

        logger.info("[%s]: Splitting assembly by contig length", sample)

        long_contigs = assembly.parent / "assembly.gt1000.fa"
        short_contigs = assembly.parent / "assembly.lt1000.fa"

        try:
            subprocess.run(
                ["seqkit", "seq", "-m", "1000", assembly],
                stdout=open(long_contigs, "w"),
                check=True,
            )
            subprocess.run(
                ["seqkit", "seq", "-M", "999", assembly],
                stdout=open(short_contigs, "w"),
                check=True,
            )
        except subprocess.CalledProcessError:
            logger.error("[%s] Contig splitting failed", sample)
            results.append([sample, "failed", str(assembly)])
            continue

        ok, corrected = get_coverage_and_correct(long_contigs, params)

        if ok:
            concatenate_fastas(corrected, short_contigs, final_fa)
            results.append([sample, "split_gt1000_merged", str(final_fa)])
        else:
            logger.warning("[%s] metaMIC failed after splitting", sample)
            results.append([sample, "failed", str(assembly)])

    # Write summary
    with open(summary, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["sample", "correction_level", "filepath"])
        writer.writerows(results)

    logger.info("Summary written to %s", summary)

def metamic_worker(contig_path: Path, params):
    sample = contig_path.parent.name.replace("_merged_assembly", "")
    
    final_fa = contig_path.parent / "assembly.metaMIC_corrected.fa"

    ok, corrected = get_coverage_and_correct(contig_path, params)

    if ok:
        concatenate_fastas(corrected, None, final_fa)
        return sample, "full", final_fa

    # Fallback: split by length
    long_contigs = contig_path.parent / "assembly.gt1000.fa"
    short_contigs = contig_path.parent / "assembly.lt1000.fa"

    try:
        subprocess.run(
            ["seqkit", "seq", "-m", "1000", contig_path],
            stdout=open(long_contigs, "w"),
            check=True,
        )
        subprocess.run(
            ["seqkit", "seq", "-M", "999", contig_path],
            stdout=open(short_contigs, "w"),
            check=True,
        )
    except subprocess.CalledProcessError:
        return sample, "failed", contig_path

    ok, corrected = get_coverage_and_correct(long_contigs, params)

    if ok:
        concatenate_fastas(corrected, short_contigs, final_fa)
        return sample, "split_gt1000_merged", final_fa

    return sample, "failed", contig_path

def run_metamic_parallel(contig_paths, params):
    summary_path = params.output_dir / "metaMIC_summary.tsv"
    results = []

    with ProcessPoolExecutor(max_workers=params.threads) as exe:
        futures = {
            exe.submit(metamic_worker, p, params): p
            for p in contig_paths
        }

        for fut in as_completed(futures):
            try:
                results.append(fut.result())
            except Exception as e:
                logger.error("metaMIC worker failed: %s", e)

    with open(summary_path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["sample", "correction_level", "filepath"])
        writer.writerows(
            [r for r in results if r is not None]
        )

    logger.info("metaMIC summary written to %s", summary_path)

    # Return corrected assemblies where available
    return [
        r[2] for r in results
        if r and Path(r[2]).exists()
    ]

# BUSCO summary
def busco_worker(contig_path, busco_outdir, lineage, threads):
    sample = contig_path.parent.name
    outdir = busco_outdir / sample
    outdir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "busco",
        "-i", str(contig_path),
        "-l", lineage,
        "-m", "genome",
        "-c", str(threads),
        "-o", sample,
        "--out_path", str(busco_outdir),
        "-f",
        "--metaeuk"
    ]

    subprocess.run(cmd, check=True)
    return sample

def run_busco_parallel(contig_paths, params):
    ensure_busco_lineage(params.busco_lineage)

    busco_dir = params.output_dir / "busco"
    busco_dir.mkdir(exist_ok=True)

    with ProcessPoolExecutor(max_workers=params.threads) as exe:
        futures = [
            exe.submit(
                busco_worker,
                p,
                busco_dir,
                params.busco_lineage,
                max(1, params.threads // 2),
            )
            for p in contig_paths
        ]

        for fut in as_completed(futures):
            try:
                fut.result()
            except Exception as e:
                logger.error("BUSCO failed: %s", e)

    return busco_dir

def ensure_busco_lineage(lineage):
    """
    Ensure that a BUSCO lineage dataset is available locally.
    Downloads it once if missing.
    """
    try:
        logger.info("Checking availability of BUSCO lineage: %s", lineage)
        subprocess.run(
            ["busco", "--list-datasets"],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
    except subprocess.CalledProcessError:
        logger.warning("Could not list BUSCO datasets Ã¢â‚¬â€ attempting download")

    # Try a no-op download check
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

def summarize_busco_json(base_dir, output_csv=None):
    """
    Summarize BUSCO results from JSON files in a directory.

    Parameters
    ----------
    base_dir : str or Path
        Directory containing BUSCO output JSON files.
    output_csv : str or Path, optional
        Path to write a summary CSV. If None, CSV is not written.

    Returns
    -------
    pd.DataFrame
        Summary table with BUSCO results and filenames.
    """
    base_dir = Path(base_dir)
    if not base_dir.exists():
        logger.error("Base directory does not exist: %s", base_dir)
        return pd.DataFrame()

    data = []

    # Iterate recursively over JSON files
    for json_file in base_dir.rglob("*.json"):
        try:
            with open(json_file, "r") as f:
                content = json.load(f)
                results = content.get("results", None)
                if results:
                    results["filename"] = json_file.name
                    results["sample"] = json_file.parent.name
                    data.append(results)
                else:
                    logger.warning("No 'results' found in %s", json_file)
        except json.JSONDecodeError:
            logger.error("Failed to parse JSON file: %s", json_file)
        except Exception as e:
            logger.error("Error processing %s: %s", json_file, e)

    if not data:
        logger.warning("No BUSCO results found in %s", base_dir)
        return pd.DataFrame()

    df = pd.DataFrame(data)

    if output_csv:
        df.to_csv(output_csv, index=False)
        logger.info("BUSCO summary written to: %s", output_csv)

    return df

# Assembly summary
def generate_seqfu_summary(contig_paths, output_dir, assembler):
    '''Run seqfu stats on all successfully assembled contigs.'''

    summary_file = Path(output_dir) / f"{assembler}_assembly_summary.tsv"

    try:
        cmd = ["seqfu", "stats", "--gc", "--csv", "-a"] + [str(p) for p in contig_paths]
        logger.info("Running seqfu on %d files", len(contig_paths))

        with summary_file.open("w") as f_out:
            subprocess.run(cmd, stdout=f_out, stderr=subprocess.PIPE, check=True, text=True)

        logger.info("Summary stats written to: %s", summary_file)

    except subprocess.CalledProcessError as e:
        logger.error("seqfu failed: %s", e.stderr)
    except OSError as e:
        logger.error("OS error running seqfu (is it installed?):  %s", e)
    except Exception as e:
        logger.error("IO error writing to summary file: %s", e)


##################### Main function #####################

def main(args):
    '''Main function to run assembly based on input arguments.'''

    params = build_params(args)

    # Normalise paths
    params.input_dir = Path(params.input_dir)
    params.output_dir = Path(params.output_dir)
    if params.unmerged_dir:
        params.unmerged_dir = Path(params.unmerged_dir)

    # Dependency check
    logger.info("Checking dependencies")
    check_dependencies(params)

    logger.info("Starting assembly pipeline")

    # Assembly
    if params.assembler == "megahit":
        process_megahit_samples(params)
    elif params.assembler == "metaspades":
        process_metaspades_samples(params)
    elif params.assembler == "idba_ud":
        process_idba_samples(params)

    # Check assemblies
    assemblies = check_assemblies(
        assembly_dir=params.output_dir,
        assembler=params.assembler,
        params=params
    )

    if not assemblies:
        logger.error("No successful assemblies found")
        sys.exit(1)

    # metaMIC correction (optional)
    if params.correction:
        logger.info("Running metaMIC correction in parallel")
        assemblies = run_metamic_parallel(assemblies, params)

    # BUSCO
    logger.info("Running BUSCO in parallel")
    busco_dir = run_busco_parallel(assemblies, params)

    summarize_busco_json(
        base_dir=busco_dir,
        output_csv=params.output_dir / "busco_summary.csv"
    )

    # SeqFu summary
    generate_seqfu_summary(
        contig_paths=assemblies,
        output_dir=params.output_dir,
        assembler=params.assembler
    )

    logger.info("Pipeline completed successfully")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "A wrapper module to run a metagenome assembler "
            "(MEGAHIT v1.0.3, MetaSPAdes genome assembler v4.0.0, or IDBA-UD v1.1.3) "
            "and output corrected contigs (using MetaMIC) and summary stats."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument("--input_dir", required=True, help="Directory with FASTQ files")
    parser.add_argument(
        "--output_dir",
        required=True,
        help="Directory for outputs." \
        "Individual assembly outputs will appear as subdirectories within this.")
    parser.add_argument(
        "--assembler",
        default="megahit",
        choices=["megahit", "metaspades", "idba_ud"],
        help="Metagenome assembler to use."
    )
    parser.add_argument(
        "--output_suffix",
        default="assembly",
        help="Suffix for output directories"
    )
    parser.add_argument(
        "--unmerged_dir",
        default=None,
        help="Directory containing unmerged (single/merged) " \
        "reads (default: input_dir).")
    parser.add_argument(
        "--correction",
        default=False, action="store_true",
        help="Run assembly correction using metaMIC after assembly. \
        This requires specifying --unmerged_dir for read files if paired \
        files are not found in input_dir."
    )
    parser.add_argument(
        "--busco_lineage",
        default="fungi_odb10",
        help="BUSCO lineage dataset to use for assembly evaluation."
    )

    p = parser.add_mutually_exclusive_group(required=False)
    p.add_argument("--tracking_sheet", help="CSV or XLSX tracking sheet")
    p.add_argument(
        "--prefix",
        nargs="+",
        help="Sample prefix(es) to match files if not providing a sample sheet."
    )
    parser.add_argument(
        "--suffix",
        help="[FOR PAIRED READS ONLY]. " \
        "Suffix pattern between sample ID and 1/2.f*q.* " \
        "(e.g., 'unmerged' in 'SampleA_unmerged_1.fq.gz')." \
        "Can be used in conjunction with --prefix or on its own."
    )
    parser.add_argument(
        "--ids",
        help="Comma-separated list of sample IDs " \
        "(used for paired reads if not using a sheet or prefix)"
    )
    parser.add_argument("--column_name", help="Column with sample IDs in tracking sheet")
    parser.add_argument("--sheet", type=int, default=0, help="Sheet index if XLSX file is used")

    g = parser.add_mutually_exclusive_group(required=True)
    g.add_argument("-p", "--paired", action="store_true", help="Use paired-end reads")
    g.add_argument("-m", "--merged", action="store_true", help="Use merged reads")
    g.add_argument("-b", "--both", action="store_true", help="Use both merged and " \
        "paired-end unmerged reads")

    parser.add_argument("--threads", type=int, default=4, help="Number of threads")
    parser.add_argument("--extra_args", nargs=argparse.REMAINDER,
        help="Extra arguments for assembly.")

    # Parse arguments first
    args = parser.parse_args()

    # Now do validation that depends on args
    if (args.column_name or args.sheet) and (args.tracking_sheet is None):
        parser.error("--column_name and --sheet require specifying an input file \
            using --tracking_sheet.")

    main(args)
