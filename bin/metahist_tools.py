import shutil
import os
import pgzip
import glob
import re
import pandas as pd
import pathlib
import subprocess
import logging
import sys
from collections import Counter
from typing import Optional, List, Dict, Tuple, Union


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)

def check_and_handle_gunzipped(input_dir): 
    """Check for gunzipped files and either gzip them or delete them if they already exist as gzipped."""
    gunzipped_files = glob.glob(os.path.join(input_dir, "*.f*q"))

    for file in gunzipped_files:
        gzipped_file = file + ".gz"
        if os.path.exists(gzipped_file):
            os.remove(file)
        else:
            with open(file, "rb") as f_in, pgzip.open(gzipped_file, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
            os.remove(file)  # Remove the original gunzipped file

def clean_and_tar(output_directory, input_directory):
    """
    Removes {output_directory}/tmp directory if it exists,
    and tars and gzips the {input_directory} into its parent directory as {input_directory}.tar.gz.
    Example: /foo/bar/raw_reads/ to /foo/bar/raw_reads.tar.gz
    """
    # Remove tmp directory if it exists
    tmp_dir = os.path.join(output_directory, "tmp")
    if os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir)
        print(f"Removed directory: {tmp_dir}")

    input_dir_abs = os.path.abspath(input_directory)
    parent_dir = os.path.dirname(input_dir_abs)
    input_dir_name = os.path.basename(input_dir_abs)
    base_name = os.path.join(parent_dir, input_dir_name)

    # Ensure the directory exists
    if not os.path.isdir(input_dir_abs):
        raise ValueError(f"Input directory does not exist: {input_dir_abs}")

    # Archive including the folder itself
    archive_path = shutil.make_archive(
        base_name,
        "gztar",
        root_dir=parent_dir,      # parent of the folder
        base_dir=input_dir_name   # folder name itself
    )

    print(f"Archived {input_directory} to {archive_path}")

def compress_file(filepath):
    if not filepath.endswith('.gz') and os.path.isfile(filepath):
        gz_path = filepath + '.gz'
        with open(filepath, 'rb') as f_in, pgzip.open(gz_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
        os.remove(filepath)
        print(f"Compressed: {os.path.basename(filepath)}")


def find_program(program_name):
    """Check if program_name is installed and return its path, else exit gracefully."""
    program_path = shutil.which(program_name)
    if program_path:
        logger.info(f"Found {program_name} at: {program_path}")
        return program_path
    else:
        logger.error(f"{program_name} executable not found in PATH. Please ensure it is installed and accessible.")
        sys.exit(1)

def find_single_reads(input_dir, prefix=None, suffix=None):
    """
    Find single-end (unpaired) read files in input_dir.

    Parameters:
    - input_dir: directory with input FASTQ files
    - prefix: optional list of prefixes to filter filenames
    - suffix: optional keyword expected before the file extension (e.g., 'merged')

    Returns:
    - List of matching single-end read FASTQ file paths. These files must not end in *1.f*q or *2.f*q
    """
    if isinstance(prefix, str):
        prefix = prefix.split(",")

    single_files = []

    fastq_pattern = re.compile(r'\.f(ast)?q(\.gz)?$', re.IGNORECASE)
    paired_pattern = re.compile(r'_[Rr]?[12](?:\.f(ast)?q(\.gz)?)?$', re.IGNORECASE)

    for fpath in glob.glob(os.path.join(input_dir, "**", "*"), recursive=True):
        if not os.path.isfile(fpath):
            continue
        fname = os.path.basename(fpath)

        if not fastq_pattern.search(fname):
            continue

        if prefix and not any(pfx in fname for pfx in prefix):
            continue

        if paired_pattern.search(fname):
            continue

        if suffix:
            suffix_pattern = re.compile(rf'{re.escape(suffix)}\.f(ast)?q(\.gz)?$', re.IGNORECASE)
            if not suffix_pattern.search(fname):
                continue

        single_files.append(fpath)  # use fpath, not os.path.join(input_dir, fname)

    # Check for duplicates
    duplicates = [f for f in single_files if single_files.count(f) > 1]
    if duplicates:
        raise ValueError(f"Duplicate files detected: {set(duplicates)}")

    return single_files

def get_ids(input_source, column_name=None, sheet=None):
    """
    Extract IDs from a spreadsheet (Excel/CSV) or from filenames in a directory.

    Parameters:
    - input_source: path to an Excel/CSV file *or* a directory of files
    - column_name: name of the column to extract IDs from (required if input_source is a file)
    - sheet: optional sheet name or index for Excel files

    Returns:
    - List of IDs (strings)
    """
    input_path = pathlib.Path(input_source)

    if input_path.is_file():
        # Handle Excel conversion if needed
        if input_source.endswith(".xlsx"):
            input_source = xlsx2csv(input_source, sheet)
            if input_source is None:
                raise FileNotFoundError("Failed to convert Excel to CSV.")

        df = pd.read_csv(input_source)

        if column_name is None:
            raise ValueError("You must provide a column_name when using a spreadsheet input.")

        if column_name not in df.columns:
            raise ValueError(f"Column '{column_name}' not found in file '{input_source}'")

        return df[column_name].dropna().astype(str).tolist()

    elif input_path.is_dir():
        # Derive IDs from filenames (like in pair_input_files)
        all_files = os.listdir(input_source)
        id_set = set()

        for file in all_files:
            match = re.match(r"(.+?)(?:_[Rr]?[12])(?:\..+)?$", file)
            if match:
                id_set.add(match.group(1))

        return sorted(id_set)

    else:
        raise FileNotFoundError(f"'{input_source}' is neither a valid file nor a directory")

def get_read_ids(
    input_source: str,
    mode: Optional[str] = None,
    prefix: Optional[List[str]] = None,
    suffix: Optional[str] = None,
    column_name: Optional[str] = None,
    sheet: Optional[str] = None
) -> List[str]:
    """
    Returns sample IDs from a directory or spreadsheet.

    Modes:
    - mode="paired": returns list of sample IDs that have both R1 and R2
    - mode="single": returns list of sample IDs for unpaired reads
    - mode=None: returns all sample IDs found (unique)

    Parameters:
    - input_source: directory or spreadsheet (CSV/XLSX)
    - mode: 'single', 'paired', or None
    - prefix: optional prefixes to filter filenames
    - suffix: optional string that must appear in the filename
    - column_name: required if input_source is a spreadsheet
    - sheet: optional sheet name for XLSX files
    """

    def keep_file(fname: str) -> bool:
        if prefix and not any(fname.startswith(pfx) for pfx in prefix):
            return False
        if suffix and suffix not in fname:
            return False
        return bool(re.match(r"(.+?)(?:_R?[12])?\.f(?:ast)?q(?:\.gz)?$", fname))

    input_path = pathlib.Path(input_source)
    ids_to_include = None

    # Spreadsheet input
    if input_path.is_file():
        if input_source.endswith(".xlsx"):
            df = pd.read_excel(input_source, sheet_name=sheet)
        else:
            df = pd.read_csv(input_source)

        if column_name is None or column_name not in df.columns:
            raise ValueError("You must provide a valid column_name for spreadsheet input.")
        ids_to_include = set(df[column_name].dropna().astype(str).tolist())
        return sorted(ids_to_include)

    # Directory input
    elif input_path.is_dir():
        all_files = sorted(os.listdir(input_source))

        if mode is None:
            id_set = set()
            for fname in all_files:
                match = re.match(r"(.+?)(?:_[Rr]?[12])(?:\..+)?$", fname)
                if match and keep_file(fname):
                    id_set.add(match.group(1))
            return sorted(id_set)

        elif mode == "paired":
            pairs: Dict[str, List[Optional[str]]] = {}
            for fname in all_files:
                if not keep_file(fname):
                    continue
                sample_id = re.match(r"(.+?)(?:_R?[12])?\.f(?:ast)?q(?:\.gz)?$", fname).group(1)
                pairs.setdefault(sample_id, [None, None])
                if "_1" in fname or "_R1" in fname:
                    pairs[sample_id][0] = fname
                elif "_2" in fname or "_R2" in fname:
                    pairs[sample_id][1] = fname
            paired_ids = [k for k, v in pairs.items() if all(v)]
            if not paired_ids:
                print("No complete read pairs found.")
            else:
                print(f"Found {len(paired_ids)} paired samples.")
            return sorted(paired_ids)

        elif mode == "single":
            single_ids = []
            for fname in all_files:
                if keep_file(fname) and not re.search(r'_[Rr]?[12]\.', fname):
                    match = re.match(r"(.+?)\.f(?:ast)?q(?:\.gz)?$", fname)
                    if match:
                        single_ids.append(match.group(1))
            if not single_ids:
                print("No single reads found.")
            else:
                print(f"Found {len(single_ids)} single reads.")
            return sorted(single_ids)

        else:
            raise ValueError("mode must be 'single', 'paired', or None")

    else:
        raise FileNotFoundError(f"'{input_source}' is neither a valid file nor directory")

def get_read_ids2(input_source, logger, paired=False, prefix=None, suffix=None, column_name=None, sheet=0):
    """
    Function to extract read IDs from a directory or spreadsheet.

    input_source : str
        Path to a directory or CSV/XLSX file.
    paired : bool, default=True
        If True, return only paired IDs (occurring twice); otherwise, return singletons.
    prefix : str, optional
        Optional filename prefix to match.
    suffix : str, optional
        Optional filename suffix to match before file extension.
    column_name : str, optional
        Required for spreadsheet input to specify the column of IDs.
    sheet : int or str, default=0
        Sheet index for Excel files.
    """

    input_path = pathlib.Path(input_source)

    # Spreadsheet input
    if input_path.is_file():
        if input_source.endswith(".xlsx"):
            xlsx_read = xlsx2csv(input_source, sheet=sheet)
            df = pd.read_csv(xlsx_read)
            logger.info(f"Using sheet: {xlsx_read}")
        else:
            df = pd.read_csv(input_source)

        if prefix:
            logger.error("Ignoring prefix due to spreadsheet input.")

        if column_name is None or column_name not in df.columns:
            logger.error("You must provide a valid column_name for spreadsheet input.")
        ids_to_include = sorted(df[column_name].dropna().astype(str).tolist())
        return ids_to_include

    # Directory input
    if input_path.is_dir():
        all_files = []
        for root, _, files in os.walk(input_source):
            for fname in files:
                all_files.append(os.path.join(root, fname))

        logger.info(f"Total files found: {len(all_files)}")

        # Filter all_files to keep FASTQ and based on prefix/suffix
        matched_files = []
        matched_ids = []
        names = [os.path.basename(f) for f in all_files]
        all_files_ids = [re.sub(r"(?:_[Rr]?[12])?\.f(?:ast)?q(?:\.gz)?$", "", n) for n in names]
        logger.info(f"File names before applying prefix/suffix filter: {len(names)}")
        if prefix or suffix:
            # Escape prefix/suffix for regex if present
            pre_pattern = prefix if prefix else ""
            suf_pattern = suffix if suffix else ""
            pattern = rf"({pre_pattern}.+?{suf_pattern})(?:_[Rr]?[12])?\.f(?:ast)?q(?:\.gz)?$"
            my_list = [re.match(pattern, n) is not None for n in names]
            logger.info(f"Files matching prefix/suffix pattern: {sum(my_list)}")
            if sum(my_list) == 0:
                logger.error("No files matched the given prefix/suffix pattern.")

            for i, match in enumerate(my_list):
                if match:
                    matched_files.append(names[i])
                    matched_ids.append(all_files_ids[i])
        else:
            matched_files = names
            matched_ids = all_files_ids
        logger.info(f"Files after applying prefix/suffix filter: {len(matched_files)}")
        if len(matched_files) == 0:
            raise ValueError("No files matched the given prefix/suffix criteria.")
        logger.info(f"Sample matched files: {matched_files} ...")
        logger.info(f"Sample matched IDs: {matched_ids} ...")
        logger.info(f"Total matched files: {len(matched_files)}")
        logger.info(f"Total matched IDs: {len(matched_ids)}")
        if len(matched_files) != len(matched_ids):
            raise ValueError("Mismatch between matched files and IDs.")
        logger.info("Counting occurrences of matched filenames...")

        # Count occurrences of full matched filenames (including extension)
        id_counts = Counter(matched_ids)
        # Build sets of paired and single files
        paired_ids = set(id_ for id_, c in id_counts.items() if c == 2)
        single_ids = set(id_ for id_, c in id_counts.items() if c == 1)
        for id_, c in id_counts.items():
            if c > 2:
                logger.info(f"Warning: File '{id_}' appears {c} times. Consider adding a suffix to filter duplicates.")

        # Return match_ids for paired or single files
        if paired:
            ids_to_include = paired_ids
        else:
            ids_to_include = single_ids
        logger.info(f"Final IDs: {ids_to_include} ...")
        return ids_to_include

    raise FileNotFoundError(f"'{input_source}' is neither a valid file nor directory")

def find_paired_files2(input_dir, prefix, recursive=True):
    """Find paired files for a given prefix in input directory.
    
    Args:
        input_dir (str): Directory to search in
        prefix (str): Sample prefix to match
        recursive (bool): Search recursively in subdirectories
        
    Returns:
        tuple: (r1_path, r2_path) if exactly two files found, (None, None) otherwise
        
    Logs:
        Info: When files are found
        Error: When incorrect number of files found
    """

    id_pair = glob.glob(os.path.join(input_dir, "**" if recursive else "", f"{prefix}*"), 
                       recursive=recursive)
    if len(id_pair) == 2:
        r1_path, r2_path = sorted(id_pair)
        logger.info(f"Running files: {r1_path} and {r2_path}")
        return r1_path, r2_path
    else:
        logger.error(f"Error: ID '{prefix}' has {len(id_pair)} matching files.")
        return None, None

def gzip_files_in_dir(input_dir, max_workers=4):
    filepaths = [
        os.path.join(input_dir, f)
        for f in os.listdir(input_dir)
        if os.path.isfile(os.path.join(input_dir, f)) and not f.endswith('.gz')
    ]
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        executor.map(compress_file, filepaths)
    print(f"All files in {input_dir} compressed")

def pair_input_files(input_directory, prefixes=None, suffix=None):
    """
    Scan input_directory for files and pair _1/_2, _r1/r2, or _R1/_R2 files.
    """
    if isinstance(prefixes, str):
        prefixes = [prefixes]

    all_files = glob.glob(os.path.join(input_directory, "**", "*"), recursive=True)
    all_files = [f for f in all_files if os.path.isfile(f)]
    file_dict = {}

    # Build regex pattern dynamically based on suffix
    if suffix:
        pattern = re.compile(
            rf"^(.*?){re.escape(suffix)}_[Rr]?[12](?:\.f(?:ast)?q(?:\.gz)?)?$"
        )
    else:
        pattern = re.compile(
            r"^(.*?)(?:_[Rr]?[12])(?:\.f(?:ast)?q(?:\.gz)?)?$"
        )

    for file in all_files:
        filename = os.path.basename(file)
        if prefixes and not any(filename.startswith(prefix) for prefix in prefixes):
            continue

        match = pattern.match(filename)
        if match:
            base = match.group(1)
            file_dict.setdefault(base, []).append(file)

    # Build pairs after collecting everything
    paired_files = []
    for base, files in file_dict.items():
        if len(files) == 2:
            paired_files.append(sorted(files))

    return paired_files



def pair_input_files_deprecated(input_directory, prefixes=None):
    """
    Scan input_directory for files and pair _1/_2, _r1/r2, or _R1/_R2 files for processing.

    Parameters:
    - input_directory: path to the directory containing FASTA/FASTQ files
    - prefixes: (optional) a single prefix string or a list of prefix strings

    Returns:
    - List of paired file lists (each with two items)
    """
    if isinstance(prefixes, str):
        prefixes = [prefixes]

    all_files = glob.glob(os.path.join(input_directory, "**", "*"), recursive=True)
    all_files = [f for f in all_files if os.path.isfile(f)]    
    paired_files = []

    file_dict = {}

    for file in all_files:
        if prefixes:
            if not any(file.startswith(prefix) for prefix in prefixes):
                continue

        # Match ending in _1, _2, _r1, _r2, _R1, or _R2 before file extension
        match = re.match(r"(.+?)(?:_[Rr]?[12])(?:\..+)?$", file)
        if match:
            base = match.group(1)
            full_path = os.path.join(input_directory, file)
            file_dict.setdefault(base, []).append(full_path)

    for base, files in file_dict.items():
        if len(files) == 2:
            files.sort()
            paired_files.append(files)

    return paired_files

def find_paired_ids(input_dir, prefixes=None, suffix=None):
    """
    Finds and pairs read files ending in fastq/fq (optionally gzipped).
    Returns a dict of {sample_id: (read1, read2)}.

    Parameters:
    - input_dir: directory with FASTQ files
    - prefixes: optional list of sample prefixes to include
    - suffix: optional string that must appear before _1/_2 (e.g. "unmerged")
    """
    dir_path = pathlib.Path(input_dir)
    read_files = sorted(dir_path.glob("*.f*q*"))

    pairs = {}
    for file in read_files:
        fname = file.name

        if prefixes and not any(fname.startswith(pfx) for pfx in prefixes):
            continue

        # Construct regex: match base ID + optional suffix + _1/_2
        if suffix:
            pattern = rf"^(.*?){re.escape(suffix)}_[Rr]?[12](?:\.f(?:ast)?q(?:\.gz)?)?$"
        else:
            pattern = r"^(.*?)[_-][Rr]?[12](?:\.f(?:ast)?q(?:\.gz)?)?$"

        match = re.match(pattern, fname)
        if not match:
            continue

        sample_id = match.group(1)  # Strip suffix from sample ID

        if suffix and not suffix in fname:
            continue  # if suffix specified but not found (shouldn't happen if regex matches, but defensive)

        pairs.setdefault(sample_id, [None, None])
        if "_1" in fname or "_R1" in fname:
            pairs[sample_id][0] = file
        elif "_2" in fname or "_R2" in fname:
            pairs[sample_id][1] = file

    paired_reads = {k: tuple(v) for k, v in pairs.items() if all(v)}

    if not paired_reads:
        logger.warning("No complete read pairs found.")
    else:
        logger.info(f"Found {len(paired_reads)} paired samples.")
    return paired_reads


def repair_reads(read1, read2, input_directory):
    """
    Repairs a pair of read files using bbmap/repair.sh.
    """
    sample_id = re.match(r"(.*?)(_R?[12])", read1.name).group(1)
    out1 = os.path.join(input_directory, f"{sample_id}_paired_2.fastq")
    out2 = os.path.join(input_directory, f"{sample_id}_paired_2.fastq")
    outsingle = os.path.join(input_directory, f"{sample_id}_single.fastq")

    cmd = [
        "repair.sh",
        f"in1={read1}",
        f"in2={read2}",
        f"out1={out1}",
        f"out2={out2}",
        f"outsingle={outsingle}",
    ]

    logger.info(f"Repairing: {sample_id}")
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        logger.info(f"Finished repairing: {sample_id}")
        if os.path.exists(outsingle):
            os.remove(outsingle)
            logger.info(f"Removed: {outsingle}")
    except subprocess.CalledProcessError as e:
        logger.error(f"repair.sh failed for {sample_id}: {e.stderr}")
    except Exception as e:
        logger.error(f"Unexpected error for {sample_id}: {e}")


def run_command(command, log_prefix, log_dir="./logs"):
    os.makedirs(log_dir, exist_ok=True)
    stdout_log = os.path.join(log_dir, f"{log_prefix}_output.log")
    stderr_log = os.path.join(log_dir, f"{log_prefix}_error.log")

    with open(stdout_log, "w") as f_out, open(stderr_log, "w") as f_err:
        f_out.write("COMMAND: " + " ".join(command) + "\n\n")  # log command at top
        f_err.write("COMMAND: " + " ".join(command) + "\n\n")
        subprocess.run(command, stdout=f_out, stderr=f_err, text=True, check=True)

def run_subprocess(command, log_prefix, log_dir="./logs", logger=None):
    os.makedirs(log_dir, exist_ok=True)
    stdout_path = os.path.join(log_dir, f"{log_prefix}_output.log")
    stderr_path = os.path.join(log_dir, f"{log_prefix}_error.log")

    with open(stdout_path, "w") as f_out, open(stderr_path, "w") as f_err:
        result = subprocess.run(command, stdout=f_out, stderr=f_err, text=True)

    if result.returncode != 0:
        if logger:
            logger.error(f"Error running {log_prefix}. See {stderr_path} for details.")
        raise subprocess.CalledProcessError(result.returncode, command)


def setup_logging(log_dir="./logs", log_file="log.out"):
    os.makedirs(log_dir, exist_ok=True)
    log_path = os.path.join(log_dir, log_file)

    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)  # Capture all levels DEBUG and above

    # Remove any existing handlers (to avoid duplicates if re-running)
    if logger.hasHandlers():
        logger.handlers.clear()

    # File handler - logs everything DEBUG and above to the file
    file_handler = logging.FileHandler(log_path)
    file_handler.setLevel(logging.DEBUG)
    file_formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(file_formatter)
    logger.addHandler(file_handler)

    # Console (stream) handler - INFO and above to stdout
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)
    console_formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    console_handler.setFormatter(console_formatter)
    logger.addHandler(console_handler)

    return logger


def xlsx2csv(file_path, sheet=0):
    if not pathlib.Path(file_path).is_file():
        logger.error(f"Error: The file '{file_path}' does not exist.")
        return None
    try:
        # If `sheet` is specified, read only that sheet; otherwise default to the first
        xlsx_read = pd.read_excel(file_path, sheet_name=sheet)
    except Exception as e:
        logger.error(f"Failed to read Excel file: {e}")
        return None

    csv_name = pathlib.Path(file_path).stem 
    dirname = os.path.dirname(file_path)
    csv_file_path = os.path.join(dirname, f"{csv_name}.csv")

    xlsx_read.to_csv(csv_file_path, index=False, header=True)
    logger.info(f"Converted '{file_path}' (sheet={sheet}) to '{csv_file_path}'")
    return csv_file_path


