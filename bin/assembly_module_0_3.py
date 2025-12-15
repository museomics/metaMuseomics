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
        extra_args=args.extra_args
    )

# Set up logging
logger = setup_logging("./logs", "assembly.log")

#  Find required executables find_program(program_name):
## "metaspades.py"
## "idba_ud"
## "fq2fa"
## "metaMIC"

# Add in optional --correction flag for metaMIC and stats after correction

# Assembly-specific functions
def run_idba_ud(sample_id, fasta_file, output_dir, output_dir_suffix=None, extra_args=None):
    ''' Function to run IDBA-UD assembly.'''

    logger.info("Starting idba-ud for %s.", fasta_file)

    assembly_dir = os.path.join(output_dir, f"{sample_id}_{output_dir_suffix}_idba_ud")
    command = [
        "idba_ud",
        "-r", fasta_file,
        "--num_threads", "1",
        "-o", assembly_dir
    ]

    if extra_args:
        command += extra_args
    try:
        subprocess.run(command, capture_output=True, text=True, check=True)
        logger.info("IDBA-UD completed for %s.", sample_id)
    except subprocess.CalledProcessError as e:
        logger.error("IDBA-UD failed for %s : %s.", sample_id, e.stderr)
        return

    logger.debug("Running command: %s", " ".join(command))

def process_idba_samples(params):
    """ Function that runs IDBA-UD assembly on samples based on user arguments,
    checks assemblies and generated summary statistics of contigs.
    """

    if os.path.exists(params.output_dir):
        pass
    else:
        os.mkdir(params.output_dir)

    fq2fa_path = find_program("fq2fa")

    if params.merged:
        if params.prefix:
            prefixes = [params.prefix]
        else:
            if params.tracking_sheet:
                prefixes = get_read_ids(params.tracking_sheet,
                    mode="single",
                    prefix=params.prefix,
                    suffix=params.suffix,
                    column_name=params.column_name,
                    sheet=params.sheet
                )
            else:
                prefixes = get_read_ids(params.input_dir,
                    mode="single",
                    prefix=params.prefix,
                    suffix=params.suffix,
                    column_name=None,
                    sheet=None
                    )
                print(prefixes, "\n")
    elif params.paired:
        ### uses ids to find pairs in the input dir but redundant because get_read_ids2 uses os.walk
        if params.prefix:
            prefixes = [params.prefix]
        else:
            if params.tracking_sheet:
                prefixes = get_read_ids2(
                    params.tracking_sheet,
                    paired=True,
                    prefix=params.prefix,
                    suffix=params.suffix,
                    column_name=params.column_name,
                    sheet=params.sheet
                )
            elif params.ids:
                prefixes = params.ids.split(",")
            else:
                prefixes = get_read_ids2(
                    params.input_dir,
                    paired=True,
                    prefix=params.prefix,
                    suffix=params.suffix,
                    column_name=None,
                    sheet=None
                )
        paired_files = {}
    else:
        logger.error("Either --paired or --merged must be specified.")
        return

    # Ensure paired_files exists for all modes (paired, both, merged)
    paired_files = {}

    unmerged_dir = params.unmerged_dir if params.unmerged_dir else params.input_dir
    print(unmerged_dir)
    extra_args = params.extra_args if hasattr(params, "extra_args") else None

    jobs = []
    for prefix in prefixes:
        logger.info("Starting IDBA-UD for %s.", prefix)
        fasta_file = None

        if params.paired:
            r1, r2 = find_paired_files2(params.input_dir, prefix)
            paired_files[prefix] = [r1, r2]
            if r1 and r2:
                if r1.endswith((".fq", ".fastq")) and r2.endswith((".fq", ".fastq")):
                    interleaved_out = os.path.join(params.input_dir, f"{prefix}_interleaved.fa")
                    cmd = [fq2fa_path, "--merge", "--filter", r1, r2, interleaved_out]
                    subprocess.run(cmd, check=True)
                    fasta_file = interleaved_out
                elif r1.endswith((".fa", ".fas", ".fasta", ".fna")) and \
                    r2.endswith((".fa", ".fas", ".fasta", ".fna")):
                    logger.error(
                        "Found pair for %s. Convert to interleaved not supported, skipping.",
                        prefix
                        )
                    continue
                elif r1.endswith((".fq", ".fastq")) and \
                    r2.endswith((".fa", ".fas", ".fasta", ".fna")):
                    r1_fasta = os.path.join(params.input_dir, f"{os.path.basename(r1)}.fa")
                    cmd = [fq2fa_path, r1, r1_fasta]
                    subprocess.run(cmd, check=True)
                    logger.warning(
                        "Mixed file formats for %s, converted R1 to FASTA. Please verify.",
                        prefix
                    )
                    fasta_file = r1_fasta
                elif r1.endswith((".fa", ".fas", ".fasta", ".fna")) and \
                    r2.endswith((".fq", ".fastq")):
                    r2_fasta = os.path.join(params.input_dir, f"{os.path.basename(r2)}.fa")
                    cmd = [fq2fa_path, r2, r2_fasta]
                    subprocess.run(cmd, check=True)
                    logger.warning(
                        "Mixed file formats for %s, converted R2 to FASTA. Please verify.",
                        prefix
                    )
                    fasta_file = r2_fasta
                else:
                    logger.error("Unrecognized file extension for: %s", r1 or r2)
                    continue
            else:
                logger.error("Could not find both paired files for %s, skipping...", prefix)
                continue

        elif params.merged:
            single_file = find_single_reads(params.input_dir, prefix)
            single_path = None
            # find_single_reads returns a list; ensure we use a single path (or skip if none)
            if not single_file:
                logger.error("No single-end reads found for %s, skipping...", prefix)
                continue
            if isinstance(single_file, (list, tuple)):
                if len(single_file) > 1:
                    logger.error(
                        "Multiple single-end files found for %s, using first: %s",
                        prefix,
                        single_file
                    )
                    continue
            else:
                single_path = single_file[0]

            if single_path.endswith((".fq", ".fastq")):
                outfile = os.path.join(params.input_dir, f"{os.path.basename(single_file)}.fa")
                cmd = [fq2fa_path, single_file, outfile]
                subprocess.run(cmd, check=True)
                fasta_file = outfile
            elif single_path.endswith((".fa", ".fasta", ".fna", ".fas")):
                fasta_file = single_path
                logger.info("Found existing FASTA file: %s", fasta_file)
            else:
                logger.error("Unrecognized file format: %s", single_path)
                continue
        else:
            logger.error("Either --paired or --merged must be specified.")
            continue

        jobs.append((prefix, fasta_file, params.output_dir, params.output_suffix, extra_args))

    # Run jobs in parallel
    with ProcessPoolExecutor(max_workers=params.threads) as executor:
        futures = [executor.submit(run_idba_ud, *job) for job in jobs]

        for future in as_completed(futures):
            try:
                future.result()
            except subprocess.CalledProcessError as e:
                logger.exception("IDBA-UD job failed: %s", str(e))

    ## Check assembly
    assembly_dir = params.output_dir
    check_assemblies(assembly_dir, assembler="idba_ud", unmerged_dir=unmerged_dir, params=params)

def process_metaspades_samples(params):
    ''' 
    Function that runs MetaSPAdes assembly on samples based on user arguments,
    checks assemblies and generated summary statistics of contigs.
    '''

    ## Set up directory
    if os.path.exists(params.output_dir):
        pass
    else:
        os.mkdir(params.output_dir)

    ## Find metaspades.py script on system
    metaspades_path = find_program("metaspades.py")

    if params.both:
        if params.prefix:
            prefixes = [params.prefix]
        else:
            if params.tracking_sheet:
                prefixes = get_read_ids(
                    params.tracking_sheet,
                    mode="single",
                    prefix=params.prefix,
                    suffix=params.suffix,
                    column_name=params.column_name,
                    sheet=params.sheet
                )
            else:
                prefixes = get_read_ids(
                    params.input_dir,
                    mode="single",
                    prefix=params.prefix,
                    suffix=params.suffix,
                    column_name=None,
                    sheet=None
                )
                print(prefixes, "\n")
    elif params.paired:
        ### uses ids to find pairs in the input dir
        if params.prefix:
            prefixes = [params.prefix]
        else:
            if params.tracking_sheet:
                prefixes = get_read_ids2(
                    params.tracking_sheet,
                    paired=True,
                    prefix=params.prefix,
                    suffix=params.suffix,
                    column_name=params.column_name,
                    sheet=params.sheet
                )
            elif params.ids:
                prefixes = params.ids.split(",")
            else:
                prefixes = get_read_ids2(
                    params.input_dir,
                    paired=True,
                    prefix=params.prefix,
                    suffix=params.suffix,
                    column_name=None,
                    sheet=None)
        paired_files = {}
    else:
        logger.error("Either --paired or --merged must be specified.")
        return
    # Ensure paired_files exists for all modes (paired, both, merged)
    paired_files = {}
    logger.info("Found files for %s", prefixes)

    unmerged_dir = params.unmerged_dir if params.unmerged_dir else params.input_dir
    print(unmerged_dir)
    extra_args = params.extra_args if hasattr(params, "extra_args") else None

    jobs = []
    for prefix in prefixes:
        logger.info("Starting Megahit for %s", prefix)

        # Merged reads
        merged_file_path = None
        if params.both:
            # Get merged input file
            merged_files = find_single_reads(params.input_dir, prefix=prefix)
            merged_file_path = merged_files[0]
            logger.info("Using %s for prefix: %s.", merged_file_path, prefix)
            if not merged_file_path and (params.merged or params.both):
                raise ValueError(f"No merged file found for prefix: {prefix}")
            # Get paired input files
            r1_path, r2_path = find_paired_files2(unmerged_dir, prefix)
            if r1_path and r2_path:
                paired_files[prefix] = [r1_path, r2_path]
            paired_set = paired_files.get(prefix, (None, None))
        elif params.paired:
            # Get paired input files
            r1_path, r2_path = find_paired_files2(unmerged_dir, prefix)
            if r1_path and r2_path:
                paired_files[prefix] = [r1_path, r2_path]
            paired_set = paired_files.get(prefix, (None, None))
        else:
            logger.error("Either --paired or --merged must be specified.")
            return

        # Add job parameters based on mode
        if params.paired:
            jobs.append(
                (
                    prefix,
                    None,
                    str(paired_set[0]),
                    str(paired_set[1]),
                    metaspades_path,
                    params.output_dir,
                    params.output_suffix,
                    extra_args,
                    False,  # both=False represented as positional flag here
                )
            )
        elif params.both:
            jobs.append(
                (
                    prefix,
                    str(merged_file_path),
                    str(paired_set[0]),
                    str(paired_set[1]),
                    metaspades_path,
                    params.output_dir,
                    params.output_suffix,
                    extra_args,
                    True,  # both=True
                )
            )


    # Run jobs in parallel
    with ProcessPoolExecutor(max_workers=params.threads) as executor:
        futures = [executor.submit(run_metaspades, *job) for job in jobs]

        for future in as_completed(futures):
            try:
                future.result()
            except subprocess.CalledProcessError as e:
                logger.error("MetaSPAdes job failed: %s.", str(e))

    ## Check assembly
    assembly_dir = params.output_dir
    check_assemblies(assembly_dir, assembler="metaspades", unmerged_dir=unmerged_dir, params=params)

def run_metaspades(sample_id, merged_path, r1_path, r2_path, metaspades_path, output_dir, output_dir_suffix, extra_args=None, both=True):
    '''Function to run MetaSPAdes with paired input.'''

    logger.info("Starting MetaSPAdes for %s", sample_id)
    assembly_dir = os.path.join(output_dir, f"{sample_id}_{output_dir_suffix}")

    # Form the command based on the available files
    command = [
        "python", 
        metaspades_path,
        "--phred-offset", "33",
        "-o", assembly_dir
        ]

    if extra_args:
        command += extra_args
    if both and (r1_path or r2_path):
        command.extend(["--merged", merged_path, "-1", r1_path, "-2", r2_path])
    elif not both and (r1_path and r2_path):
        command.extend(["-1", r1_path, "-2", r2_path])
    else:
        logger.error("ID %s: Missing necessary files to run Megahit. Skipping...", sample_id)
        return

    logger.debug("Running command: %s", " ".join(command))

    try:
        subprocess.run(command, capture_output=True, text=True, check=True)
        logger.info("MetaSPAdes completed for %s.", sample_id)
    except subprocess.CalledProcessError as e:
        logger.error("MetaSPAdes failed for %s : %s.", sample_id, e.stderr)
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
                unmerged_dir=unmerged_dir,
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
                    sheet=params.tracking_sheet
                )
            elif params.ids:
                prefixes = params.ids.split(",")
            else:
                prefixes = get_read_ids2(params.input_dir,
                    paired=True,
                    prefix=params.prefix,
                    suffix=params.suffix,
                    column_name=params.column_name,
                    sheet=params.tracking_sheet
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
    check_assemblies(assembly_dir, assembler="megahit", unmerged_dir=unmerged_dir, params=params)

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


# Assembly checks
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

def check_assemblies(assembly_dir, assembler, unmerged_dir, params, correction=False):
    '''Function to check assemblies and rerun if necessary.'''
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
        return

    for subdir in assembly_dir.iterdir():
        if not subdir.is_dir():
            continue  # skip non-directories

        sample_id = subdir.name
        contigs_path = subdir / target_file

        if contigs_path.is_file():
            logger.info(
                "A valid contig file exists for %s: %s.",
                sample_id,
                contigs_path
            )
        elif assembler == "idba_ud":
            # Try fallback scaffold.fa
            alt_file = "scaffold.fa"
            scaffold_path = subdir / alt_file
            if scaffold_path.is_file():
                logger.info(
                    "Found scaffold fallback for %s; using IDBA-UD: %s.",
                    sample_id,
                    scaffold_path
                )
            else:
                logger.error(
                    "IDBA-UD failed for %s. Consider using a different assembler or starting again",
                    sample_id
                )
        else:
            logger.info(
                "No valid contigs found for %s. Rerunning %s in case of interruption.",
                sample_id,
                assembler)
            if run_func is not None:
                run_func(subdir, sample_id)

            # Check again
            if (subdir / target_file).is_file():
                logger.info(
                    "Assembler %s successfully created contigs for %s.",
                    assembler,
                    sample_id
                )
            else:
                logger.error(
                    "Assembler %s failed for sample %s.",
                    assembler,
                    sample_id
                )

    # After all reruns, list which assemblies worked
    successful_assemblies = []

    for paths in assembly_dir.iterdir():
        if not paths.is_dir():
            continue
        contigs_path = paths / target_file
        if contigs_path.is_file():
            successful_assemblies.append(contigs_path)

    if successful_assemblies:
        logger.info(
            "Found %d successful assemblies.",
            len(successful_assemblies)
        )

        if correction:
            logger.info("Running assembly correction using metaMIC...")
            for contig_file in successful_assemblies:
                if contig_file.is_file():
                    get_coverage_and_correct(contig_file, unmerged_dir, params)
            successful_assemblies = []
            target_file = "metaMIC_corrected_contigs.fa"
            for paths in assembly_dir.iterdir():
                if not paths.is_dir():
                    continue
                contigs_path = paths / target_file
                if contigs_path.is_file():
                    successful_assemblies.append(contigs_path)

        logger.info(
            "Running BUSCO on %d assemblies in parallel...",
            len(successful_assemblies)
        )

        busco_jobs = []
        with ProcessPoolExecutor(max_workers=params.threads) as executor:
            for contig_file in successful_assemblies:
                if contig_file.is_file():
                    busco_jobs.append(
                        executor.submit(
                            run_busco_on_assemblies,
                            contig_file,
                            lineage="fungi_odb10",
                            mode="genome"
                        )
                    )

            for future in as_completed(busco_jobs):
                try:
                    future.result()
                except Exception as e:
                    logger.error("BUSCO failed for one assembly: %s", e)

        logger.info("Summarising BUSCO results...")
        summarize_busco_results(assembly_dir, assembler)
        logger.info("Generating a contig summary with seqfu...")
        generate_seqfu_summary(successful_assemblies, assembly_dir, assembler)
    else:
        logger.warning("No successful assemblies found to summarize.")


#######################################################################################
##################### Assembly summaries

# Coverage stats & Assembly polishing

def get_coverage_and_correct(contig_file, unmerged_dir, params):
    ''' Function to estimate coverage and generate an mpileup summary file for each contig file'''

    #define variables
    unmerged_dir = params.unmerged_dir if params.unmerged_dir else params.input_dir
    prefix = contig_file.parent.name.replace("_assembly", "")
    current_dir = contig_file.parent.name
    output_dir= f"{current_dir}/metaMIC_correction"
    bam_file = f"{output_dir}/{prefix}_reads.bam"
    pileup_file = f"{output_dir}/{prefix}_reads.pileup"

    # Find paired files for each contig
    logger.info("Checking files for %s", prefix)
    r1_path, r2_path = find_paired_files2(unmerged_dir, prefix)


    if not (r1_path and r2_path):
        logger.warning("No paired reads found for %s â€” skipping.", prefix)
        return

    # Define commands
    bwa_index = ["bwa", "index", str(contig_file)]
    full_cmd = (
        f"bwa mem -a -t {params.threads} {contig_file} {r1_path} {r2_path} | "
        f"samtools view -h -q 10 -m 50 -F 4 -b | "
        f"samtools sort -o {bam_file}"
    )
    mpileup_cmd = f"samtools mpileup -C 50 -A -f {contig_file} {bam_file} \
        awk '$3 != \"N\"' > {pileup_file}"
    mic_extract = [
        "metaMIC",
        "extract_feature",
        "--bam", str(bam_file),
        "-c", str(contig_file),
        "-o", str(output_dir),
        "--pileup", str(pileup_file),
        "-m", "meta"
    ]
    mic_cmd = [
        "metaMIC",
        "predict",
        "-c", str(contig_file),
        "-o", str(output_dir),
        "-a", params.assembler.upper(),
        "-m", "meta"
    ]

    logger.info("Files found for %s: %s, %s", prefix, r1_path, r2_path)

    # Index contig file
    try:
        subprocess.run(bwa_index, check=True)
    except subprocess.CalledProcessError:
        logger.error("BWA index failed for %s", {' '.join(bwa_index)})

    # run bwa mem, samtools view/filter and sort
    try:
        subprocess.run(full_cmd, shell=True, check=True)
    except subprocess.CalledProcessError:
        logger.error("Read alignment failed for %s", {' '.join(full_cmd)})

    # Samtools mpileup
    try:
        subprocess.run(mic_extract, shell=True, check=True)
    except subprocess.CalledProcessError:
        logger.error("Samtools mpileup failed for %s", {' '.join(mic_extract)})

    # metaMIC
    try:
        subprocess.run(mpileup_cmd, check=True)
    except subprocess.CalledProcessError:
        logger.error("metaMIC extract_feature failed for %s", {' '.join(mpileup_cmd)})

    try:
        subprocess.run(mic_cmd, check=True)
    except subprocess.CalledProcessError:
        logger.error("metaMIC correction failed for %s", {' '.join(mic_cmd)})


# BUSCO summary
def run_busco_on_assemblies(fasta_path, lineage="fungi_odb10", mode="genome"):
    """
    Run BUSCO on a list of assembly FASTA files.

    Parameters:
    - fasta_path: Path to the assembly FASTA file.
    - lineage: BUSCO lineage dataset (default: fungi_odb10)
    - mode: BUSCO mode (default: genome)
    """

    output_dir = fasta_path.parent / "busco_results"
    output_dir.mkdir(exist_ok=True)

    cmd = [
            "busco",
            "-i", str(fasta_path),
            "-m", mode,
            "-l", lineage,
            "--metaeuk",
            "-f",
            "-o", str(output_dir)
        ]

    find_program_out = find_program("busco")
    if not find_program_out:
        logger.error("BUSCO not found on system. Please install BUSCO to run this function.")
        return
    else:
        logger.info("Found BUSCO at: %s", find_program_out)
        logger.info("Running BUSCO for %s : %s", fasta_path, {' '.join(cmd)})
        try:
            subprocess.run(cmd, check=True)
        except subprocess.CalledProcessError as e:
            logger.error("BUSCO failed for %s: %s", fasta_path, e.stderr)

def summarize_busco_results(assembly_dir, assembler):
    '''Summarize BUSCO results from multiple assemblies into a single TSV file.'''

    busco_files = list(assembly_dir.rglob("short_summary_*.txt"))
    if output_file is None:
        output_file = assembly_dir / f"{assembler}_busco_summary.tsv"
    parse_busco_to_tsv(busco_files, output_file)


def parse_busco_to_tsv(filepaths, output_tsv):
    '''Parse multiple BUSCO short-summary files and output a TSV.'''

    parsed_rows = []

    # Regex definitions
    version_re = re.compile(r"BUSCO version is:\s*([\d\.]+)")
    lineage_re = re.compile(r"The lineage dataset is:\s*(\S+)\s*\
        (Creation date:\s*([\d\-]+),\s*number of genomes:\s*(\d+),\s*number of BUSCOs:\s*(\d+)\)"
    )
    input_re = re.compile(r"Summarized benchmarking .* for file (.+)")
    mode_re = re.compile(r"BUSCO was run in mode:\s*(\S+)")
    predictor_re = re.compile(r"Gene predictor used:\s*(\S+)")
    metrics_re = re.compile(
        r"C:(\d+\.?\d*)%\[S:(\d+\.?\d*)%,D:(\d+\.?\d*)%\],F:(\d+\.?\d*)%,M:(\d+\.?\d*)%,n:(\d+)"
    )
    count_re = {
        "complete": re.compile(r"(\d+)\s+Complete BUSCOs"),
        "single_copy": re.compile(r"(\d+)\s+Complete and single-copy BUSCOs"),
        "duplicated": re.compile(r"(\d+)\s+Complete and duplicated BUSCOs"),
        "fragmented": re.compile(r"(\d+)\s+Fragmented BUSCOs"),
        "missing": re.compile(r"(\d+)\s+Missing BUSCOs"),
        "total_groups": re.compile(r"(\d+)\s+Total BUSCO groups searched"),
    }
    assembly_re = {
        "num_scaffolds": re.compile(r"(\d+)\s+Number of scaffolds"),
        "num_contigs": re.compile(r"(\d+)\s+Number of contigs"),
        "total_length": re.compile(r"(\d+)\s+Total length"),
        "percent_gaps": re.compile(r"([\d\.]+%)\s+Percent gaps"),
        "scaffold_n50": re.compile(r"(.+?)\s+Scaffold N50"),
        "contig_n50": re.compile(r"(.+?)\s+Contigs N50"),
    }
    dep_re = re.compile(r"(\S+):\s+(\S+)")  # dependencies

    # Parse each file
    for fp in filepaths:
        fp = Path(fp)
        row = {
            "file": str(fp),
            "busco_version": None,
            "lineage_dataset": None,
            "lineage_creation_date": None,
            "lineage_num_genomes": None,
            "lineage_num_buscos": None,
            "input_file": None,
            "mode": None,
            "gene_predictor": None,
            "C": None,
            "S": None,
            "D": None,
            "F": None,
            "M": None,
            "n": None,
            "complete": None,
            "single_copy": None,
            "duplicated": None,
            "fragmented": None,
            "missing": None,
            "total_groups": None,
            "num_scaffolds": None,
            "num_contigs": None,
            "total_length": None,
            "percent_gaps": None,
            "scaffold_n50": None,
            "contig_n50": None
        }

        dependencies = {}

        with open(fp, "r") as f:
            for line in f:
                line = line.strip()

                # main metadata
                if m := version_re.search(line):
                    row["busco_version"] = m.group(1)
                if m := lineage_re.search(line):
                    row["lineage_dataset"] = m.group(1)
                    row["lineage_creation_date"] = m.group(2)
                    row["lineage_num_genomes"] = int(m.group(3))
                    row["lineage_num_buscos"] = int(m.group(4))
                if m := input_re.search(line):
                    row["input_file"] = m.group(1)
                if m := mode_re.search(line):
                    row["mode"] = m.group(1)
                if m := predictor_re.search(line):
                    row["gene_predictor"] = m.group(1)

                # metrics line
                if m := metrics_re.search(line):
                    row["C"] = float(m.group(1))
                    row["S"] = float(m.group(2))
                    row["D"] = float(m.group(3))
                    row["F"] = float(m.group(4))
                    row["M"] = float(m.group(5))
                    row["n"] = int(m.group(6))

                # BUSCO counts
                for key, regex in count_re.items():
                    if m := regex.search(line):
                        row[key] = int(m.group(1))

                # assembly stats
                for key, regex in assembly_re.items():
                    if m := regex.search(line):
                        row[key] = m.group(1)

                # dependencies
                if m := dep_re.search(line):
                    tool, ver = m.groups()
                    dependencies[tool] = ver

        # flatten dependencies into columns
        for tool, ver in dependencies.items():
            row[f"dep_{tool}"] = ver

        parsed_rows.append(row)

    # Convert to DataFrame and save as TSV
    df = pd.DataFrame(parsed_rows)
    df.to_csv(output_tsv, sep="\t", index=False)

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


def main(args):
    '''Main function to run assembly based on input arguments.'''

    params = build_params(args)
    logger.debug("Parameters: %s", params)

    logger.info("Starting assembly pipeline...")

    # Check that assembler is installed before assembly
    find_program(params.assembler)

    try:
        # Start Assembly (depending on assembler)
        if params.assembler == "idba_ud":
            process_idba_samples(params)
        elif params.assembler == "metaspades":
            process_metaspades_samples(params)
        elif params.assembler == "megahit":
            process_megahit_samples(params)

    except ValueError as e:
        logger.error("ValueError: %s", e)
        sys.exit(1)

    logger.info("All samples assembled!")


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
