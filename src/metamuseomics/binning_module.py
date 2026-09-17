import argparse
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from metahist_tools import setup_logging, run_subprocess, find_program, find_paired_files2

"""
Module for multi-sample metagenomic binning that uses fairy, samtools and bedtools for coverage estimation and 
runs multiple binners (MaxBin2, MetaBAT2, SemiBin2, MetaDecoder, CONCOCT) in parallel. 
Finally, DAS Tool is used to integrate/refine bins from the different binners. 

Requirements:
  - MaxBin2 (run_MaxBin.pl)
  - metabat2 (CLI)
  - SemiBin2 (SemiBin2 CLI)
  - MetaDecorder (metadecoder CLI)
  - CONCOCT (and its scripts: cut_up_fasta.py, concoct_coverage_table.py) ### use find_program?
  - fairy (for coverage estimation)
  - DAS_Tool (DAS_Tool.py)
  - samtools
  - bedtools

MaxBin2: Binning assembled metagenomic sequences based on an Expectation-Maximization algorithm.
MetaBAT 2: Constructs a similarity graph and utilizes LPA for partitioning
SemiBin 2: Constructs must-link and cannot-link constraints and combines them with contrastive learning to derive feature embeddings of the contigs.
MetaDecoder: Uses DPGMM for initial clustering and a semi-supervised probability model along with a modified GMM for subsequent clustering
CONCOCT: Uses a Gaussian mixture model that integrates sequence composition and coverage across multiple samples to cluster contigs into bins.
"""

##### Coverage functions #####
def get_coverage(contig_file, paired_dir, logger, threads=8):
    """
    Estimate coverages across contig file and outputs a BAM and mpileup file.
    """
    contig_file = Path(contig_file)
    unmerged_dir = Path(paired_dir)

    prefix = contig_file.parent.name.replace("_assembly", "")
    prefix = prefix.replace("_merged", "")

    current_dir = contig_file.parent
    output_dir = current_dir / "coverage_files"
    output_dir.mkdir(parents=True, exist_ok=True)

    bam_file = output_dir / f"{prefix}_reads.bam"
    pileup_file = output_dir / f"{prefix}_reads.pileup"

    logger.info("Checking files for %s", prefix)
    r1, r2 = find_paired_files2(unmerged_dir, prefix)
    if not r1 or not r2:
        logger.warning("No paired reads found for %s; skipping.", prefix)
        return None

    r1 = Path(r1)
    r2 = Path(r2)
    logger.info("Files found for %s: %s, %s", prefix, r1, r2)

    # Commands
    bwa_index = ["bwa", "index", str(contig_file)]
    align_cmd = (
        f"bwa mem -a -t {threads} {contig_file} {r1} {r2} | "
        f"samtools view -h -q 10 -m 50 -F 4 -b | "
        f"samtools sort -o {bam_file}"
    )
    mpileup_cmd = (
        f"samtools mpileup -C 50 -A -f {contig_file} {bam_file} | "
        f"awk '$3 != \"N\"' > {pileup_file}"
    )

    # Check for existing BAM and pileup files
    if bam_file.exists() and pileup_file.exists():
        logger.info("BAM and pileup files already exist for %s, skipping coverage estimation.", prefix)
        return bam_file, pileup_file
    else: 
        logger.info("Running coverage estimation for %s", prefix)
        # Run steps
        try:
            subprocess.run(bwa_index, check=True)
        except subprocess.CalledProcessError:
            logger.error("BWA index failed for %s", contig_file)
            return
        try:
            subprocess.run(align_cmd, shell=True, check=True)
        except subprocess.CalledProcessError:
            logger.error("Read alignment failed for %s", prefix)
            return
        try:
            subprocess.run(mpileup_cmd, shell=True, check=True)
        except subprocess.CalledProcessError:
            logger.error("Samtools mpileup failed for %s", prefix)
            return

        logger.info("BWA mem/Samtools coverage estimation completed for %s", prefix)
        return bam_file, pileup_file

def concoct_coverage_file(contig_file, bed_file, output, logger):
    ''' Function to run coverage file for concoct '''

    contig_file = Path(contig_file)
    prefix = contig_file.parent.name.replace("_assembly", "")
    prefix = prefix.replace("_merged", "")

    current_dir = contig_file.parent
    output_dir = current_dir / "coverage_files"
    output_dir.mkdir(parents=True, exist_ok=True)

    bam_file = output_dir / f"{prefix}_reads.bam"
    bed_file = Path(bed_file)
    concoct_file = Path(output)
    concoct_file.parent.mkdir(parents=True, exist_ok=True)
    sample_names = output_dir / "sample_names.txt"

    # Create a sample_names text file
    logger.info("Checking files for %s", prefix)
    if bam_file.is_file():
        logger.info("File found for %s: %s", prefix, bam_file)
    else:
        return

    sample_names.write_text(prefix + "\n")

    #Run CONCOCT coverage table python script
    coverage_table_program = find_program("concoct_coverage_table.py")
    if not coverage_table_program:
        logger.error("concoct_coverage_table.py not found")
        return None
    logger.info("Generating CONCOCT coverage table...")

    concoct_cmd = [
        str(coverage_table_program),
        "--samplenames",
        str(sample_names),
        str(bed_file),
        str(bam_file),
    ]
    try:
        result = subprocess.run(
            concoct_cmd,
            check=True,
            capture_output=True,
            text=True,
        )
        concoct_file.write_text(result.stdout)
    except subprocess.CalledProcessError:
        logger.error("CONCOCT coverage table generation failed for %s", prefix)
        return

    logger.info("CONCOCT coverage table generated: %s", concoct_file)

    # Clean up concoct file
    logger.info("Standardizing line endings in %s...", concoct_file)
    file_path = Path(concoct_file)
    if not file_path.exists():
        raise FileNotFoundError(f"File not found: {file_path}")

    cleaned_lines = []
    with file_path.open("r", newline="") as f:
        for line in f:
            line = line.replace("\r", "")   # remove CR
            line = line.rstrip()            # strip trailing whitespace
            if line:                        # drop empty lines
                cleaned_lines.append(line)

    # Write back in-place
    with file_path.open("w", newline="\n") as f:
        for line in cleaned_lines:
            f.write(line + "\n")

    # Validate file
    expected_fields = None
    for i, line in enumerate(cleaned_lines, start=1):
        nf = len(line.split("\t"))
        if expected_fields is None:
            expected_fields = nf
        elif nf != expected_fields:
            logger.error("ERROR: Inconsistent field count on line %d. Expected %s, but found %d",
                i, expected_fields, nf)
            return

    return concoct_file

def fairy_coverage(contig_file, paired_dir, binner, output, logger):
    ''' Function for generating coverage files for metabat2, maxbin2, semibin2'''

    contig_file = Path(contig_file)
    unmerged_dir = Path(paired_dir)

    prefix = contig_file.parent.name.replace("_assembly", "")
    prefix = prefix.replace("_merged", "")

    current_dir = contig_file.parent
    output_dir = current_dir / "coverage_files"
    output_dir.mkdir(parents=True, exist_ok=True)
    abundance_file = Path(output)
    abundance_file.parent.mkdir(parents=True, exist_ok=True)

    logger.info("Checking files for %s", prefix)
    r1, r2 = find_paired_files2(unmerged_dir, prefix)
    if not r1 or not r2:
        logger.warning("No paired reads found for %s; skipping.", prefix)
        return None

    r1 = Path(r1)
    r2 = Path(r2)
    logger.info("Found files for prefix %s: %s, %s", prefix, r1, r2)

    # Sketch files
    sketch_dir = output_dir / "sketch_dir"
    sketch_cmd = (f"fairy sketch -1 {r1} -2 {r2} -d {sketch_dir}")
    try:
        subprocess.run(sketch_cmd, shell=True, check=True)
    except subprocess.CalledProcessError:
        logger.error("Fairy coverage table generation failed for %s", prefix)
        return

    # Prepare fairy coverage matrix for this assembly using ALL sketches in sketch_dir
    sketch_files = sorted(Path(sketch_dir).glob("*.bcsp"))
    if not sketch_files:
        raise RuntimeError(f"No fairy sketch files (*.bcsp) found in {sketch_dir} -- run sketching first")

    # produce MaxBin2 abundance file
    if binner == "MaxBin2" or binner == "maxbin2":
        cmd_maxbin_cov = ["fairy", "coverage"] + [str(p) for p in sketch_files] + [str(contig_file), "--maxbin-format", "-o", str(abundance_file)]
        logger.info("Generating MaxBin2 abundance file for %s", prefix)
        run_subprocess(cmd_maxbin_cov, log_prefix = f"maxbin2_cov_{prefix}")
        return abundance_file

    # produce semibin2 abundance files (aemb format)
    if binner == "Semibin2" or binner == "semibin2":
        cmd_semibin_cov = ["fairy", "coverage"] + [str(p) for p in sketch_files] + [str(contig_file), "--aemb-format", "-o", str(abundance_file)]
        logger.info("Generating Semibin2 abundance file for %s", prefix)
        run_subprocess(cmd_semibin_cov, log_prefix = f"semibin2_cov_{prefix}")
        return abundance_file


##### Helper functions ######
def get_assembly_filename(assembler, logger):
    """
    Determine the expected assembly fasta filename based on the assembler used.
    """
    if assembler == "megahit":
        target_file = "final.contigs.fa"
    elif assembler == "metaspades":
        target_file = "scaffolds.fasta"
    elif assembler == "idba_ud":
        target_file = "final.contig.fa"
    else:
        logger.error("Unsupported assembler: %s", assembler)
        return False

    return target_file

##### Binning functions #####
def metadecoder_binning(contig_file, bam_file, logger):
    ''' Function for coverage and running metadecoder'''

    contig_file = Path(contig_file)

    prefix = contig_file.parent.name.replace("_assembly", "")
    prefix = prefix.replace("_merged", "")

    current_dir = contig_file.parent
    output_dir = current_dir / "bins_dir"
    output_dir.mkdir(parents=True, exist_ok=True)
    output_dir_metadecoder = output_dir/ "metadecoder_bins"

    outfile = output_dir / f"{prefix}_METADECODER_gsa.COVERAGE"
    seedout = output_dir / f"{prefix}_METADECODER_gsa.SEED"

    logger.info("Checking files for %s", prefix)
    if not bam_file.exists():
        raise FileNotFoundError(f"File not found: {bam_file}")

    # Coverage
    cov_cmd = (f"metadecoder coverage --threads 16 -s {bam_file} -o {outfile}")
    try:
        subprocess.run(cov_cmd, shell=True, check=True)
        logger.info("Coverage file generated for %s", bam_file)
    except subprocess.CalledProcessError:
        logger.error("Coverage file generation failed for %s", bam_file)
        return

    # Seed
    seed_cmd = (f"metadecoder seed --threads 4 -f {contig_file} -o {seedout}")
    try:
        subprocess.run(seed_cmd, shell=True, check=True)
        logger.info("Seed set for %s", contig_file)
    except subprocess.CalledProcessError:
        logger.error("Setting a seed failed for %s", contig_file)
        return

    # Cluster
    metadecoder_cmd = (f"metadecoder cluster -f {contig_file} -c {outfile} -s {seedout} -o {output_dir_metadecoder}")
    try:
        subprocess.run(metadecoder_cmd, shell=True, check=True)
        logger.info("MetaDecoder ran for %s", contig_file)
    except subprocess.CalledProcessError:
        logger.error("MetaDecoder failed for %s", contig_file)
        return

def metabat2_binning(contig_file, bam_file, output_dir, logger):
    ''' Function for coverage and running MetaBat2; requires jgi_summarize_bam_contig_depths and metabat2 in PATH'''

    contig_file = Path(contig_file)

    prefix = contig_file.parent.name.replace("_assembly", "")
    prefix = prefix.replace("_merged", "")

    current_dir = contig_file.parent
    output_dir = current_dir / "bins_dir"
    output_dir.mkdir(parents=True, exist_ok=True)
    output_dir_metabat2 = output_dir/ "metabat2_bins"

    depth_file = output_dir / f"{prefix}_depth_file.txt"

    logger.info("Checking files for %s", prefix)
    if not bam_file.exists():
        raise FileNotFoundError(f"File not found: {bam_file}")

    # Generate coverage file
    contig_depths_cmd=f"jgi_summarize_bam_contig_depths –outputDepth {depth_file} {bam_file}"

    try:
        subprocess.run(contig_depths_cmd, shell=True, check=True)
    except subprocess.CalledProcessError:
        logger.error("Depth file generation failed for %s", prefix)
        return

    logger.info("Depth file generated: %s", depth_file)

    #Run metabat
    logger.info("Running MetaBat2 on %s, using %s.", contig_file, depth_file)
    metabat_cmd=(
        f"metabat2 -i {contig_file} -a {depth_file} -o {output_dir_metabat2}/{prefix} -minContig 300 --seed 42"
        )
    try:
        subprocess.run(metabat_cmd, shell=True, check=True)
    except subprocess.CalledProcessError:
        logger.error("MetaBat2 failed for %s", prefix)
        return

    logger.info("MetaBat2 finished for: %s", depth_file)

def process_sample(sample_dir, out_dir, threads, binner, logger, assembler):
    """
    Process a single sample:
    1. Determine the assembly file based on the assembler.
    2. Estimate coverage if required by the binner.
    3. Run the specified binner(s) and collect the produced bins.
    4. Run DAS Tool to integrate bins from different binners if applicable.
    """

    binner = binner.lower()
    binner_options = ["metabat2", "maxbin2", "semibin2", "concoct", "metadecoder", "all"]
    if binner not in binner_options:
        logger.error("Unsupported binner option: %s", binner)
        return {}

    target_file = get_assembly_filename(assembler, logger)
    if not target_file:
        return {}

    sample_dir = Path(sample_dir)
    sample_id = sample_dir.name
    sample_out = out_dir / sample_id
    sample_out.mkdir(parents=True, exist_ok=True)
    logger.info("Processing sample %s -> %s", sample_id, sample_out)

    logger.info("Searching for assembly file: %s", target_file)
    assembly_candidates = next(sample_dir.glob(target_file), None)
    if assembly_candidates is None:
        raise RuntimeError(f"No assembly fasta found in {sample_dir}")

    fasta_to_contigs2bin_path = find_program("Fasta_to_Contigs2Bin.sh")

    produced_bins = {}

    # Conditions for metadecoder, metabat 2, or concocct (all require coverage estimation via BAM) 
    if binner in {"metadecoder", "concoct", "metabat2", "all"}:
        # Estimate coverage (BAM) for this sample
        bam_file = get_coverage(assembly_candidates, sample_dir, logger, threads)
        if not bam_file:
            logger.error("Coverage estimation failed for %s; skipping binners that require coverage.", sample_id)
            return produced_bins

        # --- MetaDecoder ---
        if binner in {"metadecoder", "all"} and bam_file:
            metadecoder_outdir = sample_out / "metadecoder_bins"
            metadecoder_outdir.mkdir(exist_ok=True)
            metadecoder_path = find_program("metadecoder")
            if metadecoder_path:
                logger.info("Running MetaDecoder for %s", sample_id)
                metadecoder_binning(assembly_candidates, bam_file[0], logger)
                produced_bins['metadecoder'] = str(metadecoder_outdir)
            else:
                logger.error("Coverage estimation failed; skipping MetaDecoder for %s", sample_id)

        # --- MetaBAT2 ---
        if binner in {"metabat2", "all"} and bam_file:
            metabat2_outdir = sample_out / "metabat2_bins"
            metabat2_outdir.mkdir(exist_ok=True)
            metabat2_path = find_program("metabat2")
            if metabat2_path:
                logger.info("Running MetaBAT2 for %s", sample_id)
                metabat2_binning(assembly_candidates, bam_file[0], metabat2_outdir, logger)
                produced_bins['metabat2'] = str(metabat2_outdir)
            else:
                logger.error("Coverage estimation failed; skipping MetaBAT2 for %s", sample_id)

        # --- CONCOCT ---
        # CONCOCT generally relies on per-sample mapping (BAM) + cut_up_fasta.py + concoct_coverage_table.py.
        if binner in {"concoct", "all"} and bam_file:
            bam_for_concoct = bam_file[0]
            logger.info("Running CONCOCT for %s using BAM %s", sample_id, bam_for_concoct)

            # cut up contigs into chunk size 10000 used in CONCOCT docs
            cut_fa = sample_out / f"{sample_id}.cut10k.fa"
            bed_file = sample_out / f"{sample_id}.cut10k.bed"
            logger.info("Cutting up contigs for CONCOCT")
            cmd_cut = ["cut_up_fasta.py", str(assembly_candidates), "-c", "10000", "-o", "0", "--merge_last", "-b", str(bed_file)]

            # write stdout to cut_fa
            with cut_fa.open("w") as fh:
                proc = subprocess.run(cmd_cut, stdout=fh, stderr=subprocess.PIPE, text=True)
            if proc.returncode != 0:
                logger.error("cut_up_fasta.py failed: %s", proc.stderr)
                raise subprocess.CalledProcessError(proc.returncode, cmd_cut, stderr=proc.stderr)

            # generate coverage table via concoct_coverage_table.py (requires bam index .bai)
            cov_table = sample_out / f"{sample_id}.concoct.coverage.tsv"
            coverage_file = concoct_coverage_file(assembly_candidates, bed_file, cov_table, logger)
            if coverage_file is None or not coverage_file.exists():
                logger.error("CONCOCT coverage table was not created for %s", sample_id)
                return produced_bins

            # run concoct
            concoct_outdir = sample_out / "concoct_out"
            concoct_outdir.mkdir(exist_ok=True)
            cmd_concoct = ["concoct", "-c", "400", "--composition", "--coverage_file", str(cov_table), str(cut_fa), "-b", str(concoct_outdir)]
            logger.info("Running CONCOCT for %s", sample_id)
            run_subprocess(cmd_concoct, log_prefix = f"concoct_{sample_id}")
            produced_bins['concoct'] = str(concoct_outdir)
        else:
            logger.warning("Skipping CONCOCT for %s: no BAM provided (needed for CONCOCT coverage).", sample_id)

    # Conditions for maxbin 2
    # --- MaxBin2 ---
    if binner in {"maxbin2", "all"}:
        maxbin_outdir = sample_out / "maxbin2_bins"
        maxbin_outdir.mkdir(exist_ok=True)
        maxbin_path = find_program("run_MaxBin.pl")
        if maxbin_path:
            maxbin_abund = sample_out / f"{sample_id}.maxbin.abund.tsv"
            fairy_coverage(assembly_candidates, sample_dir, "maxbin2", maxbin_abund, logger)
            if not maxbin_abund.exists():
                logger.error("MaxBin2 abundance file not found for %s; skipping MaxBin2", sample_id)
                return produced_bins
            else:
                cmd = ["run_MaxBin.pl", "-contig", str(assembly_candidates), "-abund", str(maxbin_abund), "-out", str(maxbin_outdir / "maxbin"), "-min_contig_length", "300"]
                logger.info("Running MaxBin2 for %s", sample_id)
                run_subprocess(cmd, log_prefix = f"maxbin2_{sample_id}")
                produced_bins['maxbin2'] = str(maxbin_outdir)
                logger.info("MaxBin2 finished for %s; output in %s", sample_id, maxbin_outdir)
        else:
            logger.error("MaxBin2 executable not found; skipping MaxBin2 for %s", sample_id)


    # Conditions for SemiBin2 (also uses fairy coverage but with aemb format)
    # --- SemiBin2 ---
    if binner in {"semibin2", "all"}:
        semibin_out = sample_out / "semibin2_out"
        semibin_out.mkdir(exist_ok=True)

        semibin_path = find_program("SemiBin2")
        if not semibin_path:
            logger.error("SemiBin2 was not found")
            return produced_bins

        if semibin_path:
            aemb_file = fairy_coverage(assembly_candidates, sample_dir, "semibin2", sample_out / f"{sample_id}.semibin.aemb", logger)
            if aemb_file is None or not aemb_file.exists():
                logger.error("SemiBin2 aemb file not found for %s; skipping SemiBin2", sample_id)
                return produced_bins

    # Example SemiBin2 single_easy_bin usage from Fairy README:
    cmd = [str(semibin_path), "single_easy_bin", "-i", str(assembly_candidates), str(aemb_file), "-o", str(semibin_out), "-t", str(threads)]
    logger.info("Running SemiBin2 for %s", sample_id)

    run_subprocess(cmd, log_prefix = f"SemiBin2_{sample_id}")
    produced_bins['semibin2'] = str(semibin_out)


    # --- DAS Tool ---
    # Requires: -i list of tab delimited txts of contig name and their corresponding bin
    # Requires: -l list of labels for each binner used
    # FASTA of contigs - one file for each binner used
    # -i <tsv1,tsv2,...> -l <label1,label2,...> -o prefix
    if binner == "all" or len(produced_bins) > 1:
        logger.info("Running DAS Tool to integrate bins for %s", sample_id)

        das_outdir = sample_out / "DAS_Tool_out"
        das_outdir.mkdir(exist_ok=True)

        # 1. Labels: Collect all directories that contain any of the binners in produced_bins.keys()
        labels = []
        tsv_files = []
        logger.info("DAS Tool will attempt to use bins from: %s", produced_bins.keys())

        for tool_label, bin_directory in produced_bins.items():
            tsv_file = das_outdir / f"{tool_label}_contig2bin.tsv"
            bin_files = list(Path(bin_directory).glob("*.fa"))

            if not bin_files:
                logger.warning("No bin FASTA files found in %s for DAS Tool; skipping this binner.", bin_directory)
                continue

            labels.append(tool_label)


            with tsv_file.open("w") as output_handle:
                if not fasta_to_contigs2bin_path:
                    logger.error("Fasta_to_Contigs2Bin.sh was not found")
                    return produced_bins
                subprocess.run(
                    [
                        str(fasta_to_contigs2bin_path),
                        *[str(path) for path in bin_files],
                        "-e",
                        "fasta",
                    ],
                    stdout=output_handle,
                    check=True,
                    text=True,
                )

            tsv_files.append(tsv_file)

        if not tsv_files:
            logger.warning("No DAS Tool input files were created")
            return produced_bins

        tsv_str = ",".join(str(path) for path in tsv_files)
        labels_str = ",".join(labels)

        logger.info("DAS Tool will integrate bins from: %s", labels_str)
        logger.info("TSV files for DAS Tool found: %s", tsv_str)

        # DAS_Tool usage:
        cmd = ["DAS_Tool", "-i", tsv_str, "-l", labels_str, "-o", str(das_outdir), "-c", str(assembly_candidates)]
        logger.info("Running DAS_Tool for %s (integrating: %s)", sample_id, labels_str)
        try:
            run_subprocess(cmd, log_prefix = f"DAS_Tool_{sample_id}")
            produced_bins['DAS_Tool'] = str(das_outdir)
        except Exception as e:
            logger.error("DAS_Tool failed for %s: %s", sample_id, str(e))

    return produced_bins

def main():
    parser = argparse.ArgumentParser(description=("Run multi-sample metagenomic binning using MaxBin2, MetaBAT2, SemiBin2, MetaDecoder, CONCOCT, and optionally DAS Tool."))
    parser.add_argument("-i","--input", dest="input_dir", required=True, type=Path,
        help="Directory containing the sample directories.")
    parser.add_argument("-o", "--output", dest="output_dir", required=True, type=Path,
        help="Directory where binning results will be written.")
    parser.add_argument("-a", "--assembler", choices=["megahit", "metaspades", "idba_ud"],
        required=True, help="Assembler used to generate the assemblies.")
    parser.add_argument("-b","--binner",
        choices=["metabat2", "maxbin2", "semibin2", "concoct", "metadecoder", "all"],
        default="all", help="Binner to run. Default: all."
    )
    parser.add_argument(
        "-t", "--threads",
        type=int, default=8,
        help="Number of threads passed to coverage/binners. Default: 8."
    )
    parser.add_argument(
        "-p", "--processes",
        type=int, default=1,
        help="Number of samples to process in parallel. Default: 1.",
    )
    parser.add_argument(
        "--log-dir", type=Path, default="logs",
        help="Directory where log files will be written. Default: logs.",
    )

    args = parser.parse_args()

# Validate input directory
    if not args.input_dir.is_dir():
        parser.error(f"Input directory does not exist: {args.input_dir}")

    if args.threads < 1:
        parser.error("--threads must be >= 1")

    if args.processes < 1:
        parser.error("--processes must be >= 1")

    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Set up logging.
    logger = setup_logging(log_dir = args.log_dir, log_file = "binning_module.log")

    logger.info("Input directory: %s", args.input_dir)
    logger.info("Output directory: %s", args.output_dir)
    logger.info("Assembler: %s", args.assembler)
    logger.info("Binner: %s", args.binner)
    logger.info("Threads per sample: %d", args.threads)
    logger.info("Parallel processes: %d", args.processes)

    # Find sample directories.
    sample_dirs = sorted(
        path for path in args.input_dir.iterdir()
        if path.is_dir()
    )

    if not sample_dirs:
        logger.error("No sample directories found in %s", args.input_dir)
        return 1

    logger.info("Found %d sample directories", len(sample_dirs))

    # Run samples in parallel.
    futures = {}

    with ProcessPoolExecutor(max_workers=args.processes) as executor:
        for sample_dir in sample_dirs:
            future = executor.submit(
                process_sample,
                sample_dir,
                args.output_dir,
                args.threads,
                args.binner,
                logger,
                args.assembler,
            )
            futures[future] = sample_dir

        failed = 0

        for future in as_completed(futures):
            sample_dir = futures[future]

            try:
                produced_bins = future.result()

                if produced_bins:
                    logger.info(
                        "Finished sample %s. Produced: %s",
                        sample_dir.name,
                        ", ".join(produced_bins.keys()),
                    ) 
                else:
                    logger.warning(
                        "Sample %s completed without producing bins.",
                        sample_dir.name
                    )

            except Exception:
                failed += 1
                logger.exception(
                    "Sample %s failed.",
                    sample_dir.name,
                )

    if failed:
        logger.error(
            "%d of %d samples failed.",
            failed,
            len(sample_dirs),
        )
        return 1

    logger.info("All samples completed successfully.")
    return 0

if __name__ == "__main__":
    sys.exit(main())
