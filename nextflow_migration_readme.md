# Nextflow Migration Guide for metaMuseomics

This document outlines the strategy for migrating the metaMuseomics Python pipeline to Nextflow, providing a robust, scalable workflow for processing degraded metagenomic DNA samples from museum specimens.

---

## Table of Contents

1. [Current Pipeline Architecture](#current-pipeline-architecture)
2. [Nextflow Implementation Strategy](#nextflow-implementation-strategy)
3. [Process Definitions](#process-definitions)
4. [Configuration](#configuration)
5. [Key Implementation Considerations](#key-implementation-considerations)
6. [Advanced Features](#advanced-features)
7. [Migration Benefits](#migration-benefits)
8. [Testing Strategy](#testing-strategy)

---

## Current Pipeline Architecture

The metaMuseomics pipeline consists of 9 Python wrapper scripts that process degraded metagenomic DNA samples through the following stages:

### Workflow Stages

1. **Data Acquisition** (`setup_library_dir.py`) - Download and organize sequencing data
2. **Initial QC** (`fastqc_wrapper.py`) - Quality control of raw reads
3. **Demultiplexing** (`cutadapt_2.3.py`) - Separate samples by barcodes
4. **Read Processing** (`fastp_module_3.py`) - Trim, filter, and merge reads
5. **Decontamination** (`decontam_module_v2.py`) - Remove PhiX and human contamination
6. **Complexity Estimation** (`nonpareil_parallel.py`) - Estimate metagenomic complexity
7. **Assembly** (`assembly_module_0_3.py`) - Assemble contigs (MEGAHIT, MetaSPAdes, or IDBA-UD)
8. **Assembly QC** - Quality assessment with BUSCO and seqfu statistics

### Key Features
- Optimized parameters for highly degraded DNA (hDNA)
- Parallel processing with ThreadPoolExecutor
- Support for paired-end and merged reads
- Multiple assembler options
- Comprehensive logging and error handling

---

## Nextflow Implementation Strategy

### Main Workflow Structure

**File: `main.nf`**

```groovy
#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Import modules
include { FASTQC_RAW } from './modules/fastqc'
include { DEMULTIPLEX } from './modules/demultiplex'
include { FASTP_PROCESS } from './modules/fastp'
include { DECONTAMINATE } from './modules/decontaminate'
include { NONPAREIL } from './modules/nonpareil'
include { ASSEMBLE } from './modules/assemble'
include { BUSCO_QC } from './modules/busco'
include { MULTIQC } from './modules/multiqc'

workflow {
    // Input channels
    ch_raw_reads = Channel
        .fromFilePairs(params.input_dir + '/*_{1,2}.fq.gz', checkIfExists: true)
        .ifEmpty { error "No read pairs found in ${params.input_dir}" }

    ch_samplesheet = params.tracking_sheet ?
        Channel.fromPath(params.tracking_sheet) :
        Channel.empty()

    // Process chain
    FASTQC_RAW(ch_raw_reads)

    DEMULTIPLEX(
        ch_raw_reads,
        params.i5_barcodes,
        params.i7_barcodes
    )

    FASTP_PROCESS(DEMULTIPLEX.out.reads)

    DECONTAMINATE(
        FASTP_PROCESS.out.trimmed_reads,
        params.phix_ref,
        params.human_ref
    )

    // Parallel complexity estimation
    NONPAREIL(DECONTAMINATE.out.clean_reads)

    // Combine merged and unmerged reads for assembly
    ch_assembly_input = FASTP_PROCESS.out.merged_reads
        .join(FASTP_PROCESS.out.unmerged_reads)

    ASSEMBLE(ch_assembly_input, params.assembler)

    BUSCO_QC(ASSEMBLE.out.contigs)

    // Aggregate all QC reports
    ch_multiqc = FASTQC_RAW.out.zip
        .mix(FASTP_PROCESS.out.json)
        .mix(BUSCO_QC.out.summary)
        .collect()

    MULTIQC(ch_multiqc)
}

workflow.onComplete {
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    log.info "Duration: $workflow.duration"
}
```

---

## Process Definitions

### Directory Structure

```
nextflow_pipeline/
├── main.nf
├── nextflow.config
├── modules/
│   ├── fastqc.nf
│   ├── demultiplex.nf
│   ├── fastp.nf
│   ├── decontaminate.nf
│   ├── assemble.nf
│   ├── busco.nf
│   ├── nonpareil.nf
│   └── multiqc.nf
├── subworkflows/
│   ├── preprocessing.nf
│   └── assembly_qc.nf
├── bin/
│   ├── parse_busco.py
│   └── parse_fastp_json.R
└── conf/
    ├── base.config
    ├── docker.config
    └── test.config
```

### A. FASTQC Process

**File: `modules/fastqc.nf`**

```groovy
process FASTQC_RAW {
    tag "$sample_id"
    publishDir "${params.outdir}/fastqc_raw", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    path "*.html", emit: html
    path "*.zip", emit: zip

    script:
    """
    fastqc -t ${task.cpus} -o . ${reads[0]} ${reads[1]}
    """
}
```

### B. DEMULTIPLEX Process

**File: `modules/demultiplex.nf`**

```groovy
process DEMULTIPLEX {
    tag "$sample_id"
    publishDir "${params.outdir}/demultiplexed", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    path i5_barcodes
    path i7_barcodes

    output:
    tuple val(sample_id), path("*_1.fq.gz"), path("*_2.fq.gz"), emit: reads
    path "*.stats", emit: stats

    script:
    def prefix = reads[0].baseName.replaceAll(/_R?[12].*/, '')
    """
    # Demultiplex with cutadapt
    cutadapt \\
        -e ${params.cutadapt_error_rate} \\
        --pair-adapters \\
        -g ^file:${i5_barcodes} \\
        -G ^file:${i7_barcodes} \\
        -o {name}_${prefix}_1.fq.gz \\
        -p {name}_${prefix}_2.fq.gz \\
        --action=trim \\
        --pair-filter=both \\
        ${reads[0]} ${reads[1]} \\
        > ${sample_id}_cutadapt.stats

    # Remove unknown reads
    rm -f unknown*.fq.gz

    # Sanitize read names
    for file in *_1.fq.gz; do
        base=\${file%_1.fq.gz}
        seqkit sana \${file} -o \${base}_1.sanitised.fq.gz
        seqkit sana \${base}_2.fq.gz -o \${base}_2.sanitised.fq.gz
    done

    # Pair reads
    for file in *_1.sanitised.fq.gz; do
        base=\${file%_1.sanitised.fq.gz}
        seqkit pair -1 \${file} -2 \${base}_2.sanitised.fq.gz -O .
        mv \${file} \${base}_1.fq.gz
        mv \${base}_2.sanitised.fq.gz \${base}_2.fq.gz
    done
    """
}
```

### C. FASTP_PROCESS

**File: `modules/fastp.nf`**

```groovy
process FASTP_PROCESS {
    tag "$sample_id"
    publishDir "${params.outdir}/fastp", mode: 'copy'

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    tuple val(sample_id), path("*_merged.fq"), emit: merged_reads
    tuple val(sample_id), path("*_unmerged_1.fq"), path("*_unmerged_2.fq"), emit: unmerged_reads
    tuple val(sample_id), path("*_trimmed_1.fq"), path("*_trimmed_2.fq"), emit: trimmed_reads
    path "*.html", emit: html
    path "*.json", emit: json

    script:
    """
    # Trimming and filtering step
    fastp \\
        -i ${r1} \\
        -I ${r2} \\
        --unpaired1 ${sample_id}_unpaired_1.fq \\
        --unpaired2 ${sample_id}_unpaired_2.fq \\
        --out1 ${sample_id}_trimmed_1.fq \\
        --out2 ${sample_id}_trimmed_2.fq \\
        --detect_adapter_for_pe \\
        --qualified_quality_phred=30 \\
        --trim_poly_g \\
        --correction \\
        --dedup \\
        --thread ${task.cpus} \\
        --html ${sample_id}_trim.html \\
        --json ${sample_id}_trim.json

    # Merging step
    fastp \\
        --in1 ${sample_id}_trimmed_1.fq \\
        --in2 ${sample_id}_trimmed_2.fq \\
        --merge \\
        --merged_out ${sample_id}_merged.fq \\
        --out1 ${sample_id}_unmerged_1.fq \\
        --out2 ${sample_id}_unmerged_2.fq \\
        --length_required 30 \\
        --thread ${task.cpus} \\
        --html ${sample_id}_merge.html \\
        --json ${sample_id}_merge.json

    # Generate overlap statistics
    fastp \\
        --in1 ${r1} \\
        --in2 ${r2} \\
        --stdout \\
        --merge \\
        -A -G -Q -L \\
        --thread ${task.cpus} \\
        --json /dev/null \\
        --html ${sample_id}_overlaps.html \\
        > /dev/null

    # Cleanup unpaired reads
    rm -f ${sample_id}_unpaired_*.fq
    """
}
```

### D. DECONTAMINATE Process

**File: `modules/decontaminate.nf`**

```groovy
process DECONTAMINATE {
    tag "$sample_id"
    publishDir "${params.outdir}/decontaminated", mode: 'copy',
        saveAs: { filename ->
            if (filename.endsWith('.fastq')) filename
            else "stats/$filename"
        }

    input:
    tuple val(sample_id), path(r1), path(r2)
    path phix_ref
    path human_ref

    output:
    tuple val(sample_id), path("${sample_id}_R1.fastq"), path("${sample_id}_R2.fastq"), emit: clean_reads
    path "*_stats.txt", emit: stats
    path "*_flagstats.txt", emit: flagstats

    script:
    """
    # PhiX removal with BBDuk
    bbduk.sh \\
        in1=${r1} \\
        in2=${r2} \\
        out1=${sample_id}_nophiX_R1.fq.gz \\
        out2=${sample_id}_nophiX_R2.fq.gz \\
        ref=${phix_ref} \\
        k=31 \\
        hdist=1 \\
        threads=${task.cpus} \\
        -Xmx${task.memory.toGiga()}g \\
        stats=${sample_id}_nophiX_stats.txt

    # Human decontamination with BWA + SAMtools
    bwa mem \\
        -M \\
        -t ${task.cpus} \\
        ${human_ref} \\
        ${sample_id}_nophiX_R1.fq.gz \\
        ${sample_id}_nophiX_R2.fq.gz | \\
        samtools view -b - | \\
        samtools sort -@ ${task.cpus} -o ${sample_id}_sorted.bam

    # Extract unmapped reads (non-human)
    samtools view -f 4 ${sample_id}_sorted.bam | \\
        samtools fastq \\
            -1 ${sample_id}_R1.fastq \\
            -2 ${sample_id}_R2.fastq \\
            -s ${sample_id}_singleton.fastq

    # Generate mapping statistics
    samtools flagstat -O tsv ${sample_id}_sorted.bam > ${sample_id}_human_mapping_flagstats.txt

    # Cleanup intermediate files
    rm -f ${sample_id}_nophiX_*.fq.gz ${sample_id}_sorted.bam ${sample_id}_singleton.fastq
    """
}
```

### E. ASSEMBLE Process

**File: `modules/assemble.nf`**

```groovy
process ASSEMBLE {
    tag "$sample_id"
    publishDir "${params.outdir}/assemblies", mode: 'copy'
    cpus params.assembly_threads
    memory params.assembly_memory
    time params.assembly_time
    errorStrategy 'retry'
    maxRetries 1

    input:
    tuple val(sample_id), path(merged), path(r1), path(r2)
    val assembler

    output:
    tuple val(sample_id), path("${sample_id}_assembly/"), emit: assembly_dir
    tuple val(sample_id), path("${sample_id}_assembly/${contig_file}"), emit: contigs
    path "${sample_id}_assembly/*.log", emit: logs optional true

    script:
    contig_file = assembler == 'megahit' ? 'final.contigs.fa' :
                  assembler == 'metaspades' ? 'scaffolds.fasta' :
                  'contig.fa'

    if (assembler == 'megahit')
        """
        megahit \\
            -r ${merged} \\
            -1 ${r1} \\
            -2 ${r2} \\
            -o ${sample_id}_assembly \\
            -t ${task.cpus} \\
            --min-contig-len 500
        """
    else if (assembler == 'metaspades')
        """
        # Check if continuing from previous run
        if [ -d ${sample_id}_assembly ]; then
            metaspades.py \\
                --continue \\
                -o ${sample_id}_assembly
        else
            metaspades.py \\
                --merged ${merged} \\
                -1 ${r1} \\
                -2 ${r2} \\
                --phred-offset 33 \\
                -o ${sample_id}_assembly \\
                -t ${task.cpus} \\
                -m ${task.memory.toGiga()}
        fi
        """
    else if (assembler == 'idba_ud')
        """
        # Convert FASTQ to FASTA and interleave
        fq2fa \\
            --merge \\
            --filter \\
            ${r1} \\
            ${r2} \\
            ${sample_id}_interleaved.fa

        # Run IDBA-UD
        idba_ud \\
            -r ${sample_id}_interleaved.fa \\
            --num_threads ${task.cpus} \\
            -o ${sample_id}_assembly \\
            --mink 20 \\
            --maxk 100 \\
            --step 10

        # Rename output for consistency
        if [ ! -f ${sample_id}_assembly/contig.fa ] && [ -f ${sample_id}_assembly/scaffold.fa ]; then
            ln -s scaffold.fa ${sample_id}_assembly/contig.fa
        fi
        """
    else
        error "Unknown assembler: ${assembler}. Choose from: megahit, metaspades, idba_ud"
}
```

### F. BUSCO_QC Process

**File: `modules/busco.nf`**

```groovy
process BUSCO_QC {
    tag "$sample_id"
    publishDir "${params.outdir}/busco", mode: 'copy'

    input:
    tuple val(sample_id), path(contigs)

    output:
    path "busco_${sample_id}/", emit: busco_dir
    path "busco_${sample_id}/short_summary*.txt", emit: summary
    path "${sample_id}_assembly_stats.tsv", emit: stats

    script:
    """
    # Run BUSCO
    busco \\
        -i ${contigs} \\
        -m genome \\
        -l ${params.busco_lineage} \\
        --metaeuk \\
        -f \\
        -c ${task.cpus} \\
        -o busco_${sample_id}

    # Generate assembly statistics with seqfu
    seqfu stats \\
        --gc \\
        --csv \\
        -a ${contigs} \\
        > ${sample_id}_assembly_stats.tsv
    """
}
```

### G. NONPAREIL Process

**File: `modules/nonpareil.nf`**

```groovy
process NONPAREIL {
    tag "$sample_id"
    publishDir "${params.outdir}/nonpareil", mode: 'copy'

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    path "${sample_id}.npo", emit: npo
    path "${sample_id}.npa", emit: npa
    path "${sample_id}.npc", emit: npc
    path "${sample_id}.npl", emit: npl

    script:
    """
    # Run Nonpareil on merged reads
    nonpareil \\
        -s ${r1} \\
        -T kmer \\
        -f fastq \\
        -b ${sample_id} \\
        -R ${params.nonpareil_reads} \\
        -t ${task.cpus}
    """
}
```

### H. MULTIQC Process

**File: `modules/multiqc.nf`**

```groovy
process MULTIQC {
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path('*')

    output:
    path "multiqc_report.html"
    path "multiqc_data/"

    script:
    """
    multiqc \\
        . \\
        --title "metaMuseomics Pipeline Report" \\
        --comment "Metagenomic analysis of degraded DNA samples" \\
        --config ${projectDir}/conf/multiqc_config.yaml
    """
}
```

---

## Configuration

### Main Configuration File

**File: `nextflow.config`**

```groovy
// metaMuseomics Pipeline Configuration

params {
    // Input/Output
    input_dir = "./raw_data"
    outdir = "./results"
    tracking_sheet = null

    // Reference genomes
    phix_ref = "./ref/GCA_000819615.1_ViralProj14015_genomic.fna"
    human_ref = "./ref/GCF_000001405.40_GRCh38.p14_genomic.fna"

    // Demultiplexing
    i5_barcodes = "./barcodes/i5_barcodes.fasta"
    i7_barcodes = "./barcodes/i7_barcodes.fasta"
    cutadapt_error_rate = 0.2

    // Assembly
    assembler = "megahit"  // Options: megahit, metaspades, idba_ud
    assembly_threads = 16
    assembly_memory = '64 GB'
    assembly_time = '48h'

    // BUSCO
    busco_lineage = "fungi_odb10"

    // Nonpareil
    nonpareil_reads = 100000

    // Workflow control
    skip_preprocessing = false
    skip_assembly = false
    skip_qc = false

    // Max resources
    max_cpus = 32
    max_memory = '128 GB'
    max_time = '72h'
}

// Process-specific resource allocation
process {
    // Default resources
    cpus = 2
    memory = '4 GB'
    time = '2h'

    withName: FASTQC_RAW {
        cpus = 2
        memory = '4 GB'
        time = '1h'
    }

    withName: DEMULTIPLEX {
        cpus = 4
        memory = '8 GB'
        time = '4h'
    }

    withName: FASTP_PROCESS {
        cpus = 4
        memory = '8 GB'
        time = '4h'
    }

    withName: DECONTAMINATE {
        cpus = 8
        memory = '16 GB'
        time = '8h'
    }

    withName: NONPAREIL {
        cpus = 4
        memory = '8 GB'
        time = '2h'
    }

    withName: ASSEMBLE {
        cpus = { check_max( 16 * task.attempt, 'cpus' ) }
        memory = { check_max( 64.GB * task.attempt, 'memory' ) }
        time = { check_max( 48.h * task.attempt, 'time' ) }
        errorStrategy = 'retry'
        maxRetries = 2
    }

    withName: BUSCO_QC {
        cpus = 8
        memory = '32 GB'
        time = '12h'
    }

    withName: MULTIQC {
        cpus = 1
        memory = '4 GB'
        time = '30m'
    }
}

// Execution profiles
profiles {
    standard {
        process.executor = 'local'
    }

    docker {
        docker.enabled = true
        docker.runOptions = '-u $(id -u):$(id -g)'
        process.container = 'metamuseomics:latest'
    }

    singularity {
        singularity.enabled = true
        singularity.autoMounts = true
        process.container = 'metamuseomics.sif'
    }

    conda {
        conda.enabled = true
        process.conda = "${projectDir}/environment.yml"
    }

    slurm {
        process.executor = 'slurm'
        process.queue = 'normal'
        process.clusterOptions = '--account=myaccount'
    }

    test {
        params.input_dir = "${projectDir}/test_data"
        params.outdir = "${projectDir}/test_results"
    }
}

// Report configuration
timeline {
    enabled = true
    file = "${params.outdir}/timeline.html"
}

report {
    enabled = true
    file = "${params.outdir}/report.html"
}

trace {
    enabled = true
    file = "${params.outdir}/trace.txt"
}

dag {
    enabled = true
    file = "${params.outdir}/dag.svg"
}

// Manifest
manifest {
    name = 'metaMuseomics'
    author = 'Maria Kamouyiaros'
    homePage = 'https://github.com/Kamouyiaraki/metaMuseomics'
    description = 'Nextflow pipeline for metagenomic assembly and analysis of degraded DNA samples'
    mainScript = 'main.nf'
    version = '1.0.0'
    nextflowVersion = '>=21.10.0'
}

// Function to check resource limits
def check_max(obj, type) {
    if (type == 'memory') {
        try {
            if (obj.compareTo(params.max_memory as nextflow.util.MemoryUnit) == 1)
                return params.max_memory as nextflow.util.MemoryUnit
            else
                return obj
        } catch (all) {
            println "   ### ERROR ###   Max memory '${params.max_memory}' is not valid! Using default value: $obj"
            return obj
        }
    } else if (type == 'time') {
        try {
            if (obj.compareTo(params.max_time as nextflow.util.Duration) == 1)
                return params.max_time as nextflow.util.Duration
            else
                return obj
        } catch (all) {
            println "   ### ERROR ###   Max time '${params.max_time}' is not valid! Using default value: $obj"
            return obj
        }
    } else if (type == 'cpus') {
        try {
            return Math.min( obj, params.max_cpus as int )
        } catch (all) {
            println "   ### ERROR ###   Max cpus '${params.max_cpus}' is not valid! Using default value: $obj"
            return obj
        }
    }
}
```

---

## Key Implementation Considerations

### A. Channel Management

**1. Paired-End Read Handling**
```groovy
// Create channel from paired files
ch_reads = Channel
    .fromFilePairs("${params.input_dir}/*_{1,2}.fq.gz")
    .ifEmpty { error "No paired files found" }

// Join channels by sample ID
ch_combined = ch_merged
    .join(ch_unmerged, by: 0)
```

**2. Broadcasting Channels**
```groovy
// Use same reference for all samples
ch_reference = Channel.value(file(params.human_ref))

// Broadcast to all samples
DECONTAMINATE(ch_reads, ch_reference)
```

**3. Conditional Channel Flow**
```groovy
// Split based on assembly success
ASSEMBLE.out.contigs
    .branch {
        success: it[1].size() > 0
        failed: it[1].size() == 0
    }
    .set { ch_assemblies }

// Process only successful assemblies
BUSCO_QC(ch_assemblies.success)
```

### B. Parallelization Strategy

**1. Sample-Level Parallelization**
- Each sample processed independently
- Automatic parallel execution up to available resources
- No manual ThreadPoolExecutor management needed

**2. Process-Level Threading**
```groovy
process ASSEMBLE {
    cpus 16

    script:
    """
    megahit -t ${task.cpus} ...
    """
}
```

**3. Dynamic Resource Scaling**
```groovy
process ASSEMBLE {
    cpus = { 8 * task.attempt }
    memory = { 32.GB * task.attempt }

    errorStrategy = 'retry'
    maxRetries = 3
}
```

### C. Error Handling

**1. Process Retry Logic**
```groovy
process ASSEMBLE {
    errorStrategy = { task.exitStatus in [143,137,104,134,139] ? 'retry' : 'finish' }
    maxRetries = 3

    script:
    """
    # Check for previous run
    if [ -d output ] && [ -f output/continue.txt ]; then
        megahit --continue -o output
    else
        megahit -o output ...
    fi
    """
}
```

**2. Graceful Failure Handling**
```groovy
process BUSCO_QC {
    errorStrategy = 'ignore'

    script:
    """
    busco ... || echo "BUSCO failed for ${sample_id}" > busco_failed.txt
    """
}
```

### D. Data Flow Diagram

```
┌─────────────┐
│  Raw Reads  │
└──────┬──────┘
       │
       ├──────────────────────────────────┐
       │                                  │
       v                                  v
┌─────────┐                        ┌──────────┐
│ FastQC  │                        │ MultiQC  │
└─────────┘                        │ (final)  │
       │                           └──────────┘
       v                                  ^
┌──────────────┐                         │
│ Demultiplex  │                         │
└──────┬───────┘                         │
       │                                 │
       v                                 │
┌─────────┐                              │
│  Fastp  │                              │
│ (trim + │                              │
│  merge) │                              │
└────┬────┘                              │
     │                                   │
     ├──────────────┬────────────────────┤
     │              │                    │
     v              v                    │
┌────────┐   ┌──────────┐               │
│Decontam│   │Nonpareil │───────────────┤
└───┬────┘   └──────────┘               │
    │                                    │
    v                                    │
┌─────────┐                              │
│ Assemble│                              │
│ (MEGAHIT│                              │
│MetaSPAdes│                             │
│ IDBA-UD)│                              │
└────┬────┘                              │
     │                                   │
     v                                   │
┌──────────┐                             │
│  BUSCO   │─────────────────────────────┘
│  QC      │
└──────────┘
```

### E. Module Organization

**Best Practices:**

1. **Separation of Concerns**
   - One process per module file
   - Reusable subworkflows for common patterns
   - Keep configuration separate from logic

2. **Naming Conventions**
   ```
   modules/local/    - Custom processes
   modules/nf-core/ - Community processes
   subworkflows/    - Multi-process workflows
   lib/             - Groovy helper functions
   bin/             - Executable scripts
   ```

3. **Version Control**
   - Tag module versions
   - Document parameter changes
   - Maintain backwards compatibility

### F. Utility Functions Integration

**Convert Python utilities to Nextflow:**

**File: `lib/Utils.groovy`**

```groovy
class Utils {

    // Parse sample ID from filename
    static String getSampleId(path) {
        path.baseName.replaceAll(/_R?[12].*$/, '')
    }

    // Validate input files
    static void validateInputs(params) {
        if (!params.input_dir) {
            error "Input directory not specified"
        }

        def inputDir = new File(params.input_dir)
        if (!inputDir.exists()) {
            error "Input directory does not exist: ${params.input_dir}"
        }
    }

    // Check reference files
    static void checkReferences(params) {
        ['phix_ref', 'human_ref'].each { ref ->
            def refFile = new File(params[ref])
            if (!refFile.exists()) {
                error "${ref} not found: ${params[ref]}"
            }
        }
    }
}
```

**Usage in main.nf:**
```groovy
// Import utilities
import Utils

workflow {
    // Validate inputs
    Utils.validateInputs(params)
    Utils.checkReferences(params)

    // Continue with workflow...
}
```

---

## Advanced Features

### A. Resume Capability

**1. Built-in Resume**
```bash
# Run workflow
nextflow run main.nf

# Resume after failure
nextflow run main.nf -resume
```

**2. Process-Level Checkpointing**
```groovy
process ASSEMBLE {
    storeDir "${params.outdir}/assemblies/${sample_id}"

    // Output stored permanently, won't rerun unless inputs change
}
```

### B. Conditional Workflows

**1. Skip Steps**
```groovy
workflow {
    ch_reads = Channel.fromFilePairs(params.input_dir + '/*_{1,2}.fq.gz')

    if (!params.skip_preprocessing) {
        FASTP_PROCESS(ch_reads)
        ch_processed = FASTP_PROCESS.out.trimmed_reads
    } else {
        ch_processed = ch_reads
    }

    if (!params.skip_assembly) {
        ASSEMBLE(ch_processed, params.assembler)
    }
}
```

**2. Assembler Selection**
```groovy
// Choose assembler dynamically
process ASSEMBLE {
    script:
    def assembler_cmd = params.assembler == 'megahit' ? 'megahit ...' :
                       params.assembler == 'metaspades' ? 'metaspades.py ...' :
                       params.assembler == 'idba_ud' ? 'idba_ud ...' :
                       error("Unknown assembler: ${params.assembler}")

    """
    ${assembler_cmd}
    """
}
```

### C. Assembly Correction Integration

**File: `modules/correct_assembly.nf`**

```groovy
process CORRECT_ASSEMBLY {
    tag "$sample_id"
    publishDir "${params.outdir}/corrected_assemblies", mode: 'copy'

    when:
    params.run_correction

    input:
    tuple val(sample_id), path(contigs)
    tuple val(sample_id), path(r1), path(r2)

    output:
    tuple val(sample_id), path("${sample_id}_corrected.fa"), emit: corrected
    path "${sample_id}_correction.log", emit: log

    script:
    """
    # Index contigs
    bwa index ${contigs}

    # Map reads to contigs
    bwa mem -a -t ${task.cpus} ${contigs} ${r1} ${r2} | \\
        samtools view -h -q 10 -m 50 -F 4 -b | \\
        samtools sort -o ${sample_id}.bam

    # Generate pileup
    samtools mpileup -C 50 -A -f ${contigs} ${sample_id}.bam | \\
        awk '\$3 != "N"' > ${sample_id}.pileup

    # Extract features for metaMIC
    metaMIC extract_feature \\
        --bam ${sample_id}.bam \\
        -c ${contigs} \\
        -o . \\
        --pileup ${sample_id}.pileup \\
        -m meta

    # Predict and correct errors
    metaMIC predict \\
        -c ${contigs} \\
        -o . \\
        -a ${params.assembler.toUpperCase()} \\
        -m meta

    # Rename output
    mv metaMIC_corrected_contigs.fa ${sample_id}_corrected.fa

    echo "Correction complete for ${sample_id}" > ${sample_id}_correction.log
    """
}
```

### D. MultiQC Custom Configuration

**File: `conf/multiqc_config.yaml`**

```yaml
title: "metaMuseomics Pipeline Report"
subtitle: "Analysis of degraded metagenomic DNA"
intro_text: "Quality control and assembly metrics for museum specimen samples"

report_header_info:
    - Pipeline: 'metaMuseomics'
    - Version: '1.0.0'
    - Optimized for: 'Highly degraded DNA (hDNA)'

module_order:
    - fastqc
    - fastp
    - cutadapt
    - bbduk
    - samtools
    - busco
    - custom_content

custom_data:
    nonpareil:
        file_format: 'csv'
        section_name: 'Nonpareil Complexity'
        description: 'Metagenomic complexity estimation'
        plot_type: 'linegraph'
        pconfig:
            title: 'Coverage vs Sequencing Effort'
            xlab: 'Sequencing Effort (Gbp)'
            ylab: 'Average Coverage (%)'

table_columns_visible:
    FastQC:
        percent_duplicates: True
        percent_gc: True
        total_sequences: True
    fastp:
        pct_adapter: True
        pct_duplication: True
        after_filtering_q30_rate: True
    BUSCO:
        Complete: True
        Single: True
        Duplicated: True
        Fragmented: True
        Missing: True

extra_fn_clean_exts:
    - '_trimmed'
    - '_merged'
    - '_unmerged'
    - '.sanitised'
```

### E. Subworkflow Example

**File: `subworkflows/preprocessing.nf`**

```groovy
include { FASTQC_RAW } from '../modules/fastqc'
include { DEMULTIPLEX } from '../modules/demultiplex'
include { FASTP_PROCESS } from '../modules/fastp'
include { DECONTAMINATE } from '../modules/decontaminate'

workflow PREPROCESSING {
    take:
    reads         // channel: [ val(sample_id), [ r1, r2 ] ]
    i5_barcodes   // path: i5 barcodes file
    i7_barcodes   // path: i7 barcodes file
    phix_ref      // path: PhiX reference
    human_ref     // path: Human reference

    main:
    FASTQC_RAW(reads)

    DEMULTIPLEX(reads, i5_barcodes, i7_barcodes)

    FASTP_PROCESS(DEMULTIPLEX.out.reads)

    DECONTAMINATE(
        FASTP_PROCESS.out.trimmed_reads,
        phix_ref,
        human_ref
    )

    emit:
    clean_reads = DECONTAMINATE.out.clean_reads
    qc_html = FASTQC_RAW.out.html.mix(FASTP_PROCESS.out.html)
    qc_stats = DEMULTIPLEX.out.stats
}
```

**Usage in main.nf:**
```groovy
include { PREPROCESSING } from './subworkflows/preprocessing'

workflow {
    ch_reads = Channel.fromFilePairs(params.input_dir + '/*_{1,2}.fq.gz')

    PREPROCESSING(
        ch_reads,
        params.i5_barcodes,
        params.i7_barcodes,
        params.phix_ref,
        params.human_ref
    )

    // Continue with assembly using clean reads
    ASSEMBLE(PREPROCESSING.out.clean_reads, params.assembler)
}
```

---

## Migration Benefits

### Advantages of Nextflow vs. Python Scripts

| Feature | Python Scripts | Nextflow |
|---------|---------------|----------|
| **Parallelization** | Manual ThreadPoolExecutor | Automatic per-sample parallelization |
| **Resume** | Custom checkpointing | Built-in `-resume` flag |
| **Resource Management** | Manual specification | Dynamic allocation & retry |
| **Portability** | Environment-dependent | Docker/Singularity/Conda support |
| **Scalability** | Limited to single machine | Cloud-ready (AWS, GCP, Azure) |
| **Error Handling** | Try-except blocks | Process-level strategies |
| **Monitoring** | Custom logging | Built-in reports, timeline, DAG |
| **Workflow Visualization** | None | Automatic DAG generation |
| **Job Scheduling** | Manual | SLURM, PBS, SGE integration |
| **Dependency Management** | requirements.txt | Environment per process |

### What to Keep from Python Scripts

**1. Complex Parsing Logic**
- BUSCO summary parsing (`parse_busco_to_tsv`)
- JSON aggregation functions
- Custom statistics calculations

**Solution:** Include as helper scripts in `bin/` directory

**2. R Script Integration**
- fastp JSON parser (`parse_fastp_json.R`)

**Solution:** Call from Nextflow process:
```groovy
script:
"""
Rscript ${projectDir}/bin/parse_fastp_json.R ${json_files}
"""
```

**3. Utility Functions**
- File pairing logic → Use Nextflow channel operators
- ID extraction → Groovy functions in `lib/`
- Logging → Nextflow native logging

### Performance Improvements

**Expected Gains:**

1. **Scalability:** 10-100x throughput on HPC clusters
2. **Reliability:** Automatic retry reduces manual intervention by ~80%
3. **Resource Efficiency:** Dynamic allocation prevents over/under-provisioning
4. **Development Time:** Modular processes reduce code duplication by ~60%
5. **Reproducibility:** Container support ensures consistent results across platforms

---

## Testing Strategy

### 1. Unit Testing (Individual Processes)

**Test FastQC:**
```bash
nextflow run modules/fastqc.nf \
    --input_dir test_data/raw \
    --outdir test_results/fastqc \
    -profile test
```

**Test Assembly:**
```bash
nextflow run modules/assemble.nf \
    --input_dir test_data/clean \
    --assembler megahit \
    --outdir test_results/assembly \
    -profile test
```

### 2. Integration Testing (Full Pipeline)

**Small Test Dataset:**
```bash
# Create test profile in nextflow.config
profiles {
    test {
        params.input_dir = "${projectDir}/test_data/mini_dataset"
        params.outdir = "${projectDir}/test_results"
        params.assembly_threads = 2
        params.max_cpus = 4
        params.max_memory = '8 GB'
    }
}

# Run test
nextflow run main.nf -profile test
```

**Medium Dataset (10 samples):**
```bash
nextflow run main.nf \
    --input_dir test_data/medium \
    --outdir results/medium_test \
    -resume
```

### 3. Validation Testing

**Compare to Python Pipeline:**
```bash
# Run Python pipeline
python bin/fastp_module_3.py --input_dir raw_data --output_dir py_results

# Run Nextflow pipeline
nextflow run main.nf --input_dir raw_data --outdir nf_results

# Compare outputs
diff -r py_results/sample1_merged.fq nf_results/fastp/sample1_merged.fq
md5sum py_results/*.fq nf_results/fastp/*.fq
```

### 4. Performance Benchmarking

**Benchmark Script:**
```bash
#!/bin/bash
# benchmark.sh

echo "Sample Count,Python Time,Nextflow Time,Speedup" > benchmark.csv

for n in 1 5 10 20 50; do
    # Python pipeline
    start=$(date +%s)
    python bin/assembly_module_0_3.py --input_dir data_${n} --output_dir py_out_${n}
    python_time=$(($(date +%s) - start))

    # Nextflow pipeline
    start=$(date +%s)
    nextflow run main.nf --input_dir data_${n} --outdir nf_out_${n}
    nf_time=$(($(date +%s) - start))

    speedup=$(echo "scale=2; $python_time / $nf_time" | bc)
    echo "$n,$python_time,$nf_time,$speedup" >> benchmark.csv
done
```

### 5. Continuous Integration

**GitHub Actions Workflow (`.github/workflows/ci.yml`):**
```yaml
name: Nextflow CI

on: [push, pull_request]

jobs:
  test:
    runs-on: ubuntu-latest

    steps:
    - uses: actions/checkout@v2

    - name: Install Nextflow
      run: |
        wget -qO- https://get.nextflow.io | bash
        sudo mv nextflow /usr/local/bin/

    - name: Pull Docker image
      run: docker pull metamuseomics:latest

    - name: Run pipeline tests
      run: |
        nextflow run main.nf -profile test,docker

    - name: Validate outputs
      run: |
        python tests/validate_outputs.py
```

### 6. Test Data Generation

**Create Minimal Test Dataset:**
```bash
# Extract subset of reads for testing
seqtk sample -s100 sample1_R1.fq.gz 10000 | gzip > test_data/sample1_R1.fq.gz
seqtk sample -s100 sample1_R2.fq.gz 10000 | gzip > test_data/sample1_R2.fq.gz

# Create minimal references
samtools faidx human_ref.fa chr1:1-1000000 > test_data/human_mini.fa
```

---

## Migration Checklist

### Phase 1: Setup (Week 1)
- [ ] Install Nextflow (>=21.10.0)
- [ ] Create project directory structure
- [ ] Set up test datasets (small, medium, full)
- [ ] Configure reference genomes and barcodes
- [ ] Build/pull Docker containers

### Phase 2: Core Processes (Weeks 2-3)
- [ ] Implement FASTQC module
- [ ] Implement DEMULTIPLEX module
- [ ] Implement FASTP_PROCESS module
- [ ] Implement DECONTAMINATE module
- [ ] Test each module independently

### Phase 3: Assembly (Week 4)
- [ ] Implement ASSEMBLE module (MEGAHIT)
- [ ] Add MetaSPAdes support
- [ ] Add IDBA-UD support
- [ ] Implement retry logic
- [ ] Test all three assemblers

### Phase 4: QC and Reporting (Week 5)
- [ ] Implement BUSCO_QC module
- [ ] Implement NONPAREIL module
- [ ] Implement MULTIQC module
- [ ] Configure custom MultiQC settings
- [ ] Generate test reports

### Phase 5: Integration (Week 6)
- [ ] Connect all modules in main workflow
- [ ] Implement channel logic
- [ ] Add conditional workflows
- [ ] Test full pipeline end-to-end

### Phase 6: Optimization (Week 7)
- [ ] Optimize resource allocation
- [ ] Implement dynamic retry strategies
- [ ] Add subworkflows for common patterns
- [ ] Profile performance bottlenecks

### Phase 7: Validation (Week 8)
- [ ] Compare outputs with Python pipeline
- [ ] Validate assembly quality
- [ ] Benchmark performance
- [ ] Document edge cases

### Phase 8: Documentation (Week 9)
- [ ] Write usage documentation
- [ ] Create parameter guide
- [ ] Document troubleshooting steps
- [ ] Prepare example datasets

### Phase 9: Deployment (Week 10)
- [ ] Set up HPC profile (SLURM/PBS)
- [ ] Configure cloud execution (AWS/GCP)
- [ ] Create CI/CD pipeline
- [ ] Release v1.0.0

---

## Usage Examples

### Basic Usage

```bash
# Run with default parameters
nextflow run main.nf --input_dir raw_data/

# Specify output directory
nextflow run main.nf \
    --input_dir raw_data/ \
    --outdir results/

# Resume failed run
nextflow run main.nf \
    --input_dir raw_data/ \
    --outdir results/ \
    -resume
```

### Advanced Usage

```bash
# Use MetaSPAdes assembler
nextflow run main.nf \
    --input_dir raw_data/ \
    --assembler metaspades \
    --assembly_threads 32

# Skip preprocessing (start from clean reads)
nextflow run main.nf \
    --input_dir clean_reads/ \
    --skip_preprocessing \
    --outdir results/

# Use custom references
nextflow run main.nf \
    --input_dir raw_data/ \
    --phix_ref refs/custom_phix.fa \
    --human_ref refs/custom_human.fa

# Run with Singularity on HPC
nextflow run main.nf \
    --input_dir raw_data/ \
    -profile singularity,slurm \
    --outdir results/
```

### Docker Execution

```bash
# Build Docker image
docker build -t metamuseomics:latest .

# Run pipeline in Docker
nextflow run main.nf \
    --input_dir raw_data/ \
    -profile docker
```

### Cloud Execution (AWS)

```bash
# Configure AWS Batch
nextflow run main.nf \
    --input_dir s3://mybucket/raw_data/ \
    --outdir s3://mybucket/results/ \
    -profile awsbatch \
    -bucket-dir s3://mybucket/scratch/
```

---

## Troubleshooting

### Common Issues

**1. Out of Memory Errors**
```bash
# Increase memory for assembly
nextflow run main.nf --assembly_memory '128 GB'

# Or edit nextflow.config:
process.withName:ASSEMBLE.memory = '128 GB'
```

**2. Assembly Failures**
```groovy
// Check process logs
cat .nextflow/log

// Inspect work directory
ls -lh work/a1/b2c3d4e5f6.../

// Rerun with increased resources
nextflow run main.nf -resume
```

**3. Channel Cardinality Errors**
```groovy
// Ensure channels have matching sample IDs
ch_merged.view()  // Debug channel contents
ch_unmerged.view()

// Use .join() to combine by sample ID
ch_combined = ch_merged.join(ch_unmerged, by: 0)
```

**4. Container Issues**
```bash
# Test container manually
docker run -it metamuseomics:latest bash

# Check container paths
docker run metamuseomics:latest which megahit

# Update container
docker pull metamuseomics:latest
```

---

## Resources

### Documentation
- [Nextflow Documentation](https://www.nextflow.io/docs/latest/)
- [nf-core Best Practices](https://nf-co.re/developers/guidelines)
- [Nextflow Patterns](https://nextflow-io.github.io/patterns/)

### Training
- [Nextflow Training](https://training.nextflow.io/)
- [nf-core Tutorials](https://nf-co.re/usage/tutorials)

### Community
- [Nextflow Slack](https://nextflow.io/slack-invite.html)
- [nf-core Slack](https://nf-co.re/join)

---

## Authors

**Original Python Pipeline:**
- Maria Kamouyiaros @ NHMUK (2025)

**Nextflow Migration:**
- [To be updated]

---

## License

[To be specified - should match main repository license]

---

## Citation

If you use metaMuseomics in your research, please cite:

```
Kamouyiaros, M. (2025). metaMuseomics: A pipeline for metagenomic
assembly and analysis of degraded DNA samples from museum specimens.
GitHub: https://github.com/Kamouyiaraki/metaMuseomics
```

---

## Changelog

### Version 1.0.0 (Planned)
- Initial Nextflow implementation
- Support for MEGAHIT, MetaSPAdes, and IDBA-UD assemblers
- Automated QC with FastQC, BUSCO, and MultiQC
- Docker and Singularity support
- HPC execution profiles (SLURM, PBS)
- Cloud execution support (AWS Batch)
