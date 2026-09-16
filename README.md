# metaMuseomics

*Python tools and modules for use in metagenomic assembly and analysis of degraded, metagenomic samples common in museum specimen.*

## Contents

1. [Background](#background)
2. [Getting started](#getting-started)
4. [metaMuseomics Modules: Use & Details](#metamuseomics-modules)
5. [Some extra utility tools](#extra-utility-tools)
6. [Example usage](#example-usage)
7. [References](#references)

## Background

There are 3 main modules, plus a suite of supplementary modules that can be used in a comprehensive pipeline. Each module wraps commonly used metagenomics tools into individual chunks that can be used in isolation where needed or as part of an automated pipeline. The main purpose of each module is to have optimised parameters specifically for dealing with hDNA. 

### 1. `fastp_module.py`

Performs initial adapter/quality trimming, poly-G trimming, correction and deduplication in a first step. This trimming step outputs separate stats summaries and an overlap summary. The output files of this first step are then used to merge reads, before a last round of summaries are generated. 

The poly-G trimmming, adapter, quality, correction and deduplication are all defaults, but further fastp parameters can be included in the first step by the user. 

### 2. `decontam_module.py`
The decontam module is designed to use references to remove any potential reads sourced from known contaminants. In this particular instance it is set up for PhiX and human read decontamination. 

The important point is that PhiX removal happens before human decontamination, and the two stages are separately parallelised across files, so that there are separate outputs for each stage. 

### 3. `assembly_module.py`
This is substantially larger than the other two modules: it combines assembly (with 3 assembler options: MEGAHIT, MetaSPADEs, IDBA-UD), validation/restarts, optional metaMIC correction, BUSCO assessment and SeqFu statistics. In essence `assembly_module.py` is really three pipelines in one: **assembly → optional correction → evaluation/reporting.**

In the first step, assembly happens with the chosen assembler. IDBA-UD first converts paired FASTQs into an interleaved FASTA using fq2fa; merged FASTQ can also be converted to FASTA. MEGAHIT supports three modes in the wrapper: merged reads only, paired reads only, or merged + paired reads.

Assemblies are then checked with `check_assemblies()`. This is particularly important because it acts as the gatekeeper before downstream correction/QC: it collects only assemblies with the expected contig file. IDBA-UD has a special fallback to scaffold.fa, whereas MEGAHIT and MetaSPAdes are restarted if their expected output is missing.

A correction strategy can then be implemented (with the flag `--correction`), using metaMIC, that utilises a fallback for contigs generated <1000 bp:
```
Assembly
   │
   ▼
BWA map reads → contigs
   │
   ▼
SAMtools filtering/sorting
   │
   ▼
SAMtools mpileup + AWK
   │
   ▼
metaMIC extract_feature
   │
   ▼
metaMIC predict
   │
   ├── success → corrected_contigs.fa
   │
   └── failure
         │
         ▼
   split contigs at 1000 bp
       ├── >=1000 bp → metaMIC
       └── <1000 bp → retain unchanged
              │
              ▼
       concatenate corrected + short
```

Finally, assemblies are assessed using BUSCO. BUSCO uses the selected lineage, defaults to fungi_odb10, runs in genome mode with --metaeuk, and is parallelised across assemblies. The resulting JSON files are then converted into "busco_summary.csv".

SeqFu is run once across all successful contig FASTAs with GC and CSV output enabled.


*Note: The assembly module uses process-level parallelism, whereas the preprocessing/decontamination modules use thread pools. This is mainly due to the relatively heavyweight assembler/BUSCO processes versus the lighter orchestration around external commands.*

## Getting started 

### Installation

There are two installation options: 

Option A) To install the python package only (assuming you have bioinformatics software already installed) run: 

```
python3.12 -m pip install metaMuseomics
```

Option B) To install a complete conda environment, run:  

```
conda env create -f environment.yml
conda activate metamuseomics
```

If you are planning on using the decontaminate module, then run the following in the directory above your working directory to download the reference sequences (HG38 *Homo sapiens* reference genome GRCh38.p14 (GCF_000001405.40) & PhiX NC_001422.1) used: 

```
mkdir ref
cd ./ref

# PhiX
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/viral/Sinsheimervirus_phiX174/latest_assembly_versions/GCA_000819615.1_ViralProj14015/GCA_000819615.1_ViralProj14015_genomic.fna.gz

# Homo sapiens 
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz
gunzip GCF_000001405.40_GRCh38.p14_genomic.fna.gz
bwa index GCF_000001405.40_GRCh38.p14_genomic.fna
```

Alternatively, you can specify your reference sequences using the arguments `--phiX_ref` and `--human-ref`. 

In theory any reference can be used in place of either, however you should keep in mind that as the PhiX genome is so small, `bbduk.sh` from BBTools is used, while BWA MEM is used for the much larger (and more resource intense) `--human-ref`. 

**Important note on IDBA-UD:**
[IDBA-UD](https://github.com/loneknightpy/idba) may throw an error in the presence of short insert sizes. A solution to this has been made available [here](https://www.seqanswers.com/forum/bioinformatics/bioinformatics-aa/24625-250bp-reads-in-idba_ud). 


## metaMuseomics Modules

The modules are designed to be used either as independent tools or as part of a metagenomic pipeline for historical (maybe museum-derived) data, as shown below.

![Modules flowchart](https://github.com/museomics/metaMuseomics/blob/main/img/flowchart.svg)

Module | Main function / role | Key functions | Dependencies |
|---|---|---|---|
| `nonpareil_module.py`  | Metagenomic sequencing-complexity estimation. Runs Nonpareil on FASTQ files and summarises coverage, redundancy, diversity and sequencing effort required for 95% coverage. | `run_nonpareil()` – runs Nonpareil per FASTQ; `parse_npo()` – extracts metrics from `.npo`; `batch_run_nonpareil()` – processes all matching FASTQs and creates summary CSV.| **External:** `nonpareil`. **Python:** standard library|
| `fastp_module.py` | Read preprocessing/QC. Trims and filters paired-end FASTQs, merges overlapping pairs, produces overlap plots, and creates summary statistics. Samples can be processed in parallel. | `run_fastp_trim()`: adapter/quality trimming, poly-G trimming and deduplication; `run_fastp_merge()`: merges paired reads; `run_fastp_overlap_plot()`: creates overlap HTML; `generate_seqkit_stats()`: FASTQ statistics; `run_fastp_json_merge()`: combines fastp JSON results via R; `process_sample()`: per-sample workflow | **External**: `fastp`, `seqkit`, `R`, `seqpy-tools` (run_command, setup_logging, get_read_ids2). **Python:** standard library.|
| `decontam_module.py` | Host/contamination removal. Removes PhiX contamination, maps reads against a human reference, retains unmapped reads, and repairs paired-end files. Supports paired or merged reads and parallel processing. | `run_bbduk()` – removes PhiX with BBDuk; `run_bwa_mem_and_samtools()` – maps to human reference and extracts unmapped reads | **References**: PhiX genome and human GRCh38 reference FASTAs; **External**: `bbduk.sh`/`BBMap`, `bwa`, `samtools`, `seqpy-tools` (run_subprocess, repair_reads, get_read_ids2, find_paired_files2, find_single_reads, setup_logging). **Python**: `pandas` + standard library.|
| `assembly_module.py` | Metagenomic assembly + assembly evaluation/correction. Runs one of MEGAHIT, MetaSPAdes or IDBA-UD, validates/restarts failed assemblies, optionally applies metaMIC correction, then evaluates assemblies with BUSCO and produces contig statistics with SeqFu. | Assembly: `run_megahit()`, `run_metaspades()`, `run_idba_ud()`. **Validation**: `assemblies_exist_for_all_samples()`, `check_assemblies()`, restart functions. **Correction:** `get_coverage_and_correct()`, `run_metamic_correction()`. **Evaluation**: `run_busco_parallel()`, `ensure_busco_lineage()`, `summarize_busco_json()`, `generate_seqfu_summary()`. | **External**: `megahit`/`metaspades`/`idba_ud`, `busco`, `seqfu`, `metaMIC`, `bwa`, `samtools`, `seqkit`, `awk`, `seqpy-tools` (find_single_reads, find_paired_files2, get_read_ids, get_read_ids2, find_program, setup_logging). **Python**: `pandas`, `gzip`, `json` + standard library.|
| `cutadapt_module.py`   | Barcode demultiplexing and FASTQ sanitisation. Uses i5/i7 barcodes to demultiplex raw paired-end reads, sanitises reads, repairs pairing, and generates read statistics.| `run_cutadapt()` – barcode-based demultiplexing; `seqkit_sanitize()` – sanitises FASTQs; `find_files()` – identifies R1/R2 pairs; `seqkit_pair()` – repairs/pairs reads; `generate_seqkit_stats()` – QC statistics. | **External:** `cutadapt`, `seqkit`, `seqpy-tools` (clean_and_tar, pair_input_files, check_and_handle_gunzipped). **Python:** `pandas`, `pgzip` and standard library. |


## Extra utility tools and wrappers

Module | Main role | Key functions | Dependencies |
|---|---|---|---|
| `busco_wrapper.py` | BUSCO summary. Finds existing assemblies, runs BUSCO, summarises BUSCO results in a spreadhseet. This is the same process that occurs in the assembly module, but made available for independent use outside of the assembly step. | `run_busco`: runs BUSCO on found assemblies. `collect_busco_summary()`: collects BUSCO outputs and summarises into a spreadsheet. | **External:** `BUSCO`. **Python:** standard library, `json` + `pandas` |
| `fastqc_wrapper.py` | Raw-read quality control. Finds paired FASTQs using a tracking sheet, sample prefix/suffix, or all files, then runs FastQC on each pair.| `run_fastqc()`: executes FastQC.| **External:** `fastqc`, `seqpy-tools` (`pair_input_files`, `xlsx2csv`). **Python:**  standard library + `pandas`.|  
| `kraken2_wrapper.py` | Kraken2 taxonomic classification on multiple FASTA samples in parallel. It includes functions for processing samples, running Kraken2, and outputs FASTAs grouped by Family (other ranks for future development). |  `run_kraken()`: wrapper to run kraken2 with any specified DB; `taxid_to_family()`: uses tax IDs to find Family rank, `write_family_fastas()`: outputs FASTAs of all contigs belonging to the same Family| **External**: `Kraken2`, **Python**: standard library. |
| `ids2csv.py`| Sample metadata to FASTQ path mapping. Takes sample IDs from a CSV, searches a project directory for corresponding trimmed paired FASTQs, and adds `forward`/`reverse` path columns.    | `get_ids()`: extracts IDs from CSV; `find_files()`: locates matching R1/R2 FASTQs; `write_to_csv()`: adds paths and writes the output CSV.| **Python:**  standard library. |
|`parse_fastp_json.R` | Parse JSON files into a summary CSV file after fastp trimming. This is coded into the fastp module but can be used as a standalone tool. | `json_parse()`: takes all JSON files, reads, extracts primary information and outputs into a dataframe that is then exported as a CSV| `R` (jsonlite)|
| `setup_library_dir.py` | Library/project data setup. Downloads sequencing files from URLs, validates MD5 checksums, extracts TAR archives, moves FASTQs into `raw_data`, and cleans up the resulting directory. | `download_file()` – downloads individual files with `wget`; top-level script handles URL parsing, directory creation, checksum validation, TAR extraction and file organisation.| **External:** `wget`, `md5sum`, `tar`. **Python:** standard library|

## Example usage

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
