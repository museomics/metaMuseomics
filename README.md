# metaMuseomics
Python package of scripts for use in metagenomic assembly and analysis of degraded, metagenomic samples common in museum specimen

Each module wraps commonly used tools into individual chunks that can be used in isolation where needed or as part of an automated pipeline. The main purpose of each wrapper is to have optimised parameters specifically for dealing with hDNA. In some cases (e.g., the assembly module) this has involved changing some source code in assembler installs (IDBA-UD), so make sure to keep this in an isolated environment if you use IDBA-UD for DNA that is not highly degraded.  

## Contents

1. [Getting started](#getting-started)
2. [Individual Modules: Use & Details](#metamuseomics-modules)
3. [Extra tools](#extra-utility-tools)
4. [Tutorials and Use Cases](https://github.com/Kamouyiaraki/metaMuseome/blob/main/Tutorial.md)

   
## Getting started 

### Installation

pip install

### The dependencies: 




## metaMuseomics Modules

![Modules flowchart](https://github.com/museomics/metaMuseomics/blob/main/img/flowchart.svg)

Module | Main function / role | Key functions | Dependencies |
|---|---|---|---|
| `setup_library_dir.py` | Library/project data setup. Downloads sequencing files from URLs, validates MD5 checksums, extracts TAR archives, moves FASTQs into `raw_data`, and cleans up the resulting directory. | `download_file()` – downloads individual files with `wget`; top-level script handles URL parsing, directory creation, checksum validation, TAR extraction and file organisation.| **External:** `wget`, `md5sum`, `tar`. **Python:** standard library|
| `nonpareil_module.py`  | Metagenomic sequencing-complexity estimation. Runs Nonpareil on FASTQ files and summarises coverage, redundancy, diversity and sequencing effort required for 95% coverage.            | `run_nonpareil()` – runs Nonpareil per FASTQ; `parse_npo()` – extracts metrics from `.npo`; `batch_run_nonpareil()` – processes all matching FASTQs and creates summary CSV.| **External:** `nonpareil`. **Python:** standard library|
| `fastp_module.py` | Read preprocessing/QC. Trims and filters paired-end FASTQs, merges overlapping pairs, produces overlap plots, and creates summary statistics. Samples can be processed in parallel. | `run_fastp_trim()`: adapter/quality trimming, poly-G trimming and deduplication; `run_fastp_merge()`: merges paired reads; `run_fastp_overlap_plot()`: creates overlap HTML; `generate_seqkit_stats()`: FASTQ statistics; `run_fastp_json_merge()`: combines fastp JSON results via R; `process_sample()`: per-sample workflow | **External**: `fastp`, `seqkit`, `R`, `seqpy-tools` (run_command, setup_logging, get_read_ids2). **Python:** standard library.|
| `decontam_module.py` | Host/contamination removal. Removes PhiX contamination, maps reads against a human reference, retains unmapped reads, and repairs paired-end files. Supports paired or merged reads and parallel processing. | `run_bbduk()` – removes PhiX with BBDuk; `run_bwa_mem_and_samtools()` – maps to human reference and extracts unmapped reads | **References**: PhiX genome and human GRCh38 reference FASTAs; **External**: `bbduk.sh`/`BBMap`, `bwa`, `samtools`, `seqpy-tools` (run_subprocess, repair_reads, get_read_ids2, find_paired_files2, find_single_reads, setup_logging). **Python**: `pandas` + standard library.|
| `assembly_module.py` | Metagenomic assembly + assembly evaluation/correction. Runs one of MEGAHIT, MetaSPAdes or IDBA-UD, validates/restarts failed assemblies, optionally applies metaMIC correction, then evaluates assemblies with BUSCO and produces contig statistics with SeqFu. | Assembly: `run_megahit()`, `run_metaspades()`, `run_idba_ud()`. **Validation**: `assemblies_exist_for_all_samples()`, `check_assemblies()`, restart functions. **Correction:** `get_coverage_and_correct()`, `run_metamic_correction()`. **Evaluation**: `run_busco_parallel()`, `ensure_busco_lineage()`, `summarize_busco_json()`, `generate_seqfu_summary()`. | **External**: `megahit`/`metaspades`/`idba_ud`, `busco`, `seqfu`, `metaMIC`, `bwa`, `samtools`, `seqkit`, `awk`, `seqpy-tools` (find_single_reads, find_paired_files2, get_read_ids, get_read_ids2, find_program, setup_logging). **Python**: `pandas` + standard library.|                                                                                                                                                                                                                 || `cutadapt_module.py`   | Barcode demultiplexing and FASTQ sanitisation. Uses i5/i7 barcodes to demultiplex raw paired-end reads, sanitises reads, repairs pairing, and generates read statistics.               | `run_cutadapt()` – barcode-based demultiplexing; `seqkit_sanitize()` – sanitises FASTQs; `find_files()` – identifies R1/R2 pairs; `seqkit_pair()` – repairs/pairs reads; `generate_seqkit_stats()` – QC statistics. | **External:** `cutadapt`, `seqkit`, `seqpy-tools` (clean_and_tar, pair_input_files, check_and_handle_gunzipped). **Python:** `pandas`, `pgzip` and standard library. |


### Modules

### Quick run modules

### Outputs 

## Extra utility tools 

Module | Main role | Key functions | Dependencies |
|---|---|---|---|
|`parse_fastp_json.R` | Parse JSON files into a summary CSV file after fastp trimming. This is coded into the fastp module but can be used as a standalone tool. | `json_parse()`: takes all JSON files, reads, extracts primary information and outputs into a dataframe that is then exported as a CSV| `R` (jsonlite)|
| `ids2csv.py`| Sample metadata to FASTQ path mapping. Takes sample IDs from a CSV, searches a project directory for corresponding trimmed paired FASTQs, and adds `forward`/`reverse` path columns.    | `get_ids()`: extracts IDs from CSV; `find_files()`: locates matching R1/R2 FASTQs; `write_to_csv()`: adds paths and writes the output CSV.| **Python:**  standard library. |
| `fastqc_wrapper.py`    | Raw-read quality control. Finds paired FASTQs using a tracking sheet, sample prefix/suffix, or all files, then runs FastQC on each pair.| `run_fastqc()`: executes FastQC.| **External:** `fastqc`, `seqpy-tools` (`pair_input_files`, `xlsx2csv`). **Python:**  standard library + `pandas`.|            

## References
- Simon Andrews, 2010. FastQC:  A Quality Control Tool for High Throughput Sequence Data [Online]. Available online at: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
- Shifu Chen, 2023. Ultrafast one-pass FASTQ data preprocessing, quality control, and deduplication using fastp. iMeta 2: e107. https://doi.org/10.1002/imt2.107
- Shifu Chen, Yanqing Zhou, Yaru Chen, Jia Gu, 2018. fastp: an ultra-fast all-in-one FASTQ preprocessor, Bioinformatics, 34, 17, i884–i890. https://doi.org/10.1093/bioinformatics/bty560
- Bushnell B. – sourceforge.net/projects/bbmap/
- Heng Li, Richard Durbin, 2009. Fast and accurate short read alignment with Burrows–Wheeler transform, Bioinformatics, 25, 14, 1754–1760, https://doi.org/10.1093/bioinformatics/btp324
- Heng Li, 2013 Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. arXiv:1303.3997v2 [q-bio.GN].
- Li, D., Liu, C-M., Luo, R., Sadakane, K., and Lam, T-W., 2015. MEGAHIT: An ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph. Bioinformatics, doi: 10.1093/bioinformatics/btv033 [PMID: 25609793].
- Li, D., Luo, R., Liu, C.M., Leung, C.M., Ting, H.F., Sadakane, K., Yamashita, H. and Lam, T.W., 2016. MEGAHIT v1.0: A Fast and Scalable Metagenome Assembler driven by Advanced Methodologies and Community Practices. Methods.
- Nurk S, Meleshko D, Korobeynikov A, Pevzner PA, 2017. metaSPAdes: a new versatile metagenomic assembler. Genome Res.;27(5):824-834. doi: 10.1101/gr.213959.116. Epub 2017 Mar 15. PMID: 28298430; PMCID: PMC5411777.
- Yu Peng, Henry C. M. Leung, S. M. Yiu, Francis Y. L. Chin, 2012. IDBA-UD: a de novo assembler for single-cell and metagenomic sequencing data with highly uneven depth, Bioinformatics, 28, 11, 1420–1428, https://doi.org/10.1093/bioinformatics/bts174
- Hofmeyr, S., Egan, R., Georganas, E. et al. 2020. Terabase-scale metagenome coassembly with MetaHipMer. Sci Rep 10, 10689.


*MetaMuseomics was written by Maria Kamouyiaros @ NHMUK (2025).*
