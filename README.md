# BactogenMiner

[![Nextflow](https://img.shields.io/badge/Nextflow-%E2%89%A522.10.0-23AA62.svg)](https://www.nextflow.io/)
[![DSL2](https://img.shields.io/badge/DSL-2-23AA62.svg)](https://www.nextflow.io/)
[![Docker](https://img.shields.io/badge/Docker-enabled-2496ED.svg)](https://www.docker.com/)
[![Conda](https://img.shields.io/badge/Conda-enabled-44A833.svg)](https://docs.conda.io/)

**BactogenMiner** is a scalable, reproducible Nextflow (DSL2) workflow designed for end-to-end bacterial whole-genome sequencing (WGS) data analysis. It automates pre-assembly read quality screening, read trimming, *de novo* assembly, assembly evaluation, genome annotation, taxonomic identification, sequence typing, virulence/plasmid screening, and phylogeny-ready manifest generation.

---

## Table of Contents

- [Overview & Key Features](#overview--key-features)
- [Pipeline Architecture](#pipeline-architecture)
- [Workflow Steps](#workflow-steps)
- [Prerequisites & Dependencies](#prerequisites--dependencies)
- [Reference Databases](#reference-databases)
- [Installation](#installation)
- [Usage](#usage)
  - [1. Single-sample Mode](#1-single-sample-mode)
  - [2. Multi-sample Mode (Samplesheet)](#2-multi-sample-mode-samplesheet)
  - [3. Execution Profiles](#3-execution-profiles)
- [Parameters Reference](#parameters-reference)
  - [Input / Output](#input--output)
  - [Pre-Assembly Read Screening & Quality Gating](#pre-assembly-read-screening--quality-gating)
  - [Taxonomic Classification](#taxonomic-classification)
  - [Genome Annotation](#genome-annotation)
  - [Assembly, QC & Typing](#assembly-qc--typing)
  - [Compute Resources](#compute-resources)
- [Output Directory Structure](#output-directory-structure)
- [Tools and Citations](#tools-and-citations)
- [Author and Contact](#author-and-contact)

---

## Overview & Key Features

- **Pre-Assembly Read Screening & Quality Gating**:
  - Automatically assesses read count, total basepairs, R1/R2 forward-reverse balance, estimated genome length, and estimated sequencing coverage prior to assembly using **fastq-scan** and **Mash** (adapted from Theiagen / TheiaProk screening standards).
  - Flags and skips degraded, under-covered, or contaminated samples early to avoid wasted compute.
- **Automated QC & Preprocessing**: Performs **FastQC** quality checks and **fastp** adapter trimming, poly-G clipping, and quality filtering on paired-end or single-end Illumina reads.
- **Robust *De Novo* Assembly**: Employs **SPAdes** for contig/scaffold reconstruction and assembly graph generation.
- **Multifaceted Assembly QC**:
  - Contig continuity and N50 statistics using **QUAST**.
  - Assembly graph visual exploration using **Bandage**.
  - Core conserved gene completeness benchmarking using **BUSCO** (`bacteria_odb12` lineage).
- **Flexible Genome Annotation**:
  - Rapid, standard bacterial annotation via **Prokka** (default).
  - High-precision annotation via **Bakta** (with automated database download or local database support).
- **Dual Taxonomic Classification Support**:
  - Ultra-fast genomic distance-based classification via **GAMBIT** with automated top-species resolution.
  - K-mer based metagenomic classification with **Kraken2** and abundance re-estimation with **Bracken**.
- **In Silico Sequence Typing**: Multi-Locus Sequence Typing (**MLST**) against PubMLST typing schemes with automatic schema detection.
- **Virulence & Plasmid Screening**: Screen assemblies for virulence factors (**VFDB**) and plasmid replicons (**PlasmidFinder**) via **ABRICATE**.
- **Phylogeny Manifest Generation**: Centralizes all generated scaffold assemblies into an assembly manifest (`assembly_sample_sheet.txt`) ready for downstream phylogenetic workflows (e.g. SKA, RAxML).
- **Reproducible Execution**: Built-in containerized profiles (**Docker** via BioContainers & custom images) and environment recipes (**Conda**).

---

## Pipeline Architecture

```mermaid
flowchart TD
    subgraph Inputs ["Input Data"]
        R1["Raw FASTQ Reads (Single / Paired)"]
        SS["CSV Samplesheet"]
    end

    subgraph Screen ["1. Read Screening & Quality Gating"]
        R1 --> CHECK_READS["CHECK_READS\n(fastq-scan & Mash)\nLength, Coverage, R1/R2 Balance"]
        SS --> CHECK_READS
        CHECK_READS -->|FAIL| DROP["Filtered Out\n(Warning Logged)"]
        CHECK_READS -->|PASS| PASS_READS["Validated Reads"]
    end

    subgraph Preprocessing ["2. Quality Control & Trimming"]
        PASS_READS --> FASTQC["FastQC\n(Raw Read QC)"]
        PASS_READS --> FASTP["fastp\n(Adapter Trimming & Quality Filtering)"]
    end

    subgraph Assembly ["3. De Novo Assembly"]
        FASTP --> SPADES["SPAdes\n(De Novo Genome Assembly)"]
    end

    subgraph AssemblyQC ["4. Assembly Quality Control"]
        SPADES --> BANDAGE["Bandage\n(Graph Visualization)"]
        SPADES --> QUAST["QUAST\n(Assembly Metrics & N50)"]
        SPADES --> BUSCO["BUSCO\n(Genome Completeness)"]
    end

    subgraph Annotation ["5. Genome Annotation"]
        SPADES --> PROKKA["Prokka\n(Default Annotation)"]
        SPADES --> BAKTA["Bakta\n(Optional Annotation)"]
        BAKTA_DB["Bakta DB Download (if needed)"] -.-> BAKTA
    end

    subgraph Taxonomy ["6. Species Classification"]
        SPADES --> GAMBIT["GAMBIT\n(Genomic Distance)"]
        GAMBIT --> GS_GAMBIT["Top Species Extraction"]
        SPADES --> KRAKEN2["Kraken2\n(K-mer Classification)"]
        KRAKEN2 --> BRACKEN["Bracken\n(Abundance Estimation)"]
        BRACKEN --> GS_BRACKEN["Top Species Extraction"]
    end

    subgraph Typing ["7. Sequence Typing & Profiling"]
        SPADES --> MLST["MLST\n(PubMLST Scheme Typing)"]
        SPADES --> VFDB["ABRICATE (VFDB)\n(Virulence Factors)"]
        SPADES --> PLASMID["ABRICATE (PlasmidFinder)\n(Plasmid Replicons)"]
    end

    subgraph Phylogeny ["8. Phylogeny Preparation"]
        SPADES --> MAKE_SHEET["MAKE_ASSEMBLY_SHEET\n(Aggregates all scaffolds)"]
    end

    subgraph Outputs ["Results Directory"]
        CHECK_READS --> OUT_QC["results/sample_id/ReadScreen/"]
        FASTQC --> OUT["results/sample_id/"]
        FASTP --> OUT
        BANDAGE --> OUT
        QUAST --> OUT
        BUSCO --> OUT
        PROKKA --> OUT
        BAKTA --> OUT
        GAMBIT --> OUT
        BRACKEN --> OUT
        MLST --> OUT
        VFDB --> OUT
        PLASMID --> OUT
        MAKE_SHEET --> PHYLO_OUT["results/phylogeny/assembly_sample_sheet.txt"]
    end
```

---

## Workflow Steps

1. **Pre-Assembly Read Screening & Quality Gating (`CHECK_READS`)**:
   - Computes read count and basepair metrics via `fastq-scan`.
   - Sketches reads via `Mash` to estimate genome length and sequencing coverage before performing *de novo* assembly.
   - Verifies paired-end read balance (`--min_proportion`), minimum read count (`--min_reads`), minimum total basepairs (`--min_basepairs`), estimated genome length bounds (`--min_genome_length`, `--max_genome_length`), and minimum coverage (`--min_coverage`).
   - Samples failing quality criteria are flagged in `<sample_id>_read_screen.tsv` and filtered out.
2. **Read QC and Trimming (`FASTQC` & `FASTP`)**:
   - **FastQC**: Evaluates per-base sequence quality, GC content distribution, and duplication levels.
   - **fastp**: Performs automated adapter detection, quality filtering, poly-G clipping, and front/tail trimming (`--cut_front --cut_tail`).
3. **De Novo Assembly (`SPADES`)**:
   - **SPAdes**: Assembles high-quality trimmed reads into scaffolds (`<sample_id>.fasta`) and assembly graphs (`<sample_id>.fastg`). Supports single-end and paired-end datasets.
4. **Assembly Quality Control (`BANDAGE`, `QUAST`, `BUSCO`)**:
   - **Bandage**: Generates visual representations of the assembly graph topology (`<sample_id>_Bandage_Img.jpg`).
   - **QUAST**: Computes assembly continuity metrics (N50, L50, contig count, GC content) and generates an interactive HTML summary.
   - **BUSCO**: Quantifies genome completeness using universal bacterial single-copy orthologs (`bacteria_odb12`).
5. **Genome Annotation (`PROKKA` / `BAKTA`)**:
   - **Prokka** *(Default)*: Rapidly annotates CDS, tRNA, rRNA, and signal peptides, outputting GFF3, GenBank, and FASTA files.
   - **Bakta** *(Optional)*: High-precision bacterial genome annotation. Downloads required databases automatically (`BAKTA_BD`) if a local database path is not specified.
6. **Taxonomic Classification (`GAMBIT` / `KRAKEN2` + `BRACKEN`)**:
   - **GAMBIT** *(Default)*: Rapid species identification using genomic signatures against curated reference databases.
   - **Kraken2** & **Bracken** *(Optional)*: Exact k-mer matching with Bayesian abundance re-estimation and top-species resolution.
7. **Sequence Typing & Virulence/Plasmid Profiling (`MLST` & `ABRICATE`)**:
   - **MLST**: Scans contigs against PubMLST databases to determine sequence type (ST) and allele profiles.
   - **ABRICATE**: Mass screening of contigs for:
     - **VFDB**: Bacterial virulence factors.
     - **PlasmidFinder**: Plasmid replicon typing.
8. **Phylogeny Sample Sheet Compilation (`MAKE_ASSEMBLY_SHEET`)**:
   - Compiles real file paths of all assembled scaffolds into `results/phylogeny/assembly_sample_sheet.txt` for downstream phylogenetic analyses (e.g. core-genome SNP alignment or phylogenetics).

---

## Prerequisites & Dependencies

- **Nextflow** (>= 22.10.0)
- **Container Engine or Package Manager**:
  - **Docker** (recommended for seamless reproducibility; uses official BioContainers and `mmcj/check_reads:1.1`)
  - **Conda / Mamba** (using environment definitions in `envs/`)

Hardware requirements depend on the number and size of genomes analyzed. By default, processes use available CPUs and up to 28 GB RAM (configurable via parameters).

---

## Reference Databases

Ensure required databases are accessible before running species classification or advanced annotation:

| Tool | Purpose | Source / Download Link |
|---|---|---|
| **GAMBIT** | Genomic distance classification | [GAMBIT Databases](https://gambit-genomics.readthedocs.io/en/latest/databases.html) |
| **Kraken2 / Bracken** | K-mer classification & abundance | [Ben Langmead AWS Indexes](https://benlangmead.github.io/aws-indexes/k2) |
| **Bakta DB** *(Optional)* | Full or light annotation database | Automatically downloaded if `--use_bakta` or via [Bakta Documentation](https://github.com/oschwengers/bakta#database) |
| **BUSCO Lineage** | Gene completeness benchmarking | Auto-downloaded during run (`bacteria_odb12`) |
| **VFDB / PlasmidFinder** | Virulence and plasmid typing | Bundled within ABRicate |

---

## Installation

Clone the repository:

```bash
git clone https://github.com/MOISECHRIST/bacteria_typing-phylogeny.git
cd bacteria_typing-phylogeny
```

Ensure Nextflow is installed and functioning:

```bash
nextflow -version
```

---

## Usage

### 1. Single-sample Mode

To analyze a single sample with paired-end reads:

```bash
nextflow run main.nf \
  -profile docker \
  --reads "data/sample01_{1,2}.fastq.gz" \
  --sample_name "sample01" \
  --use_gambit \
  --gambit_db "/path/to/gambit/db" \
  --outdir "results"
```

For single-end data:

```bash
nextflow run main.nf \
  -profile conda \
  --reads "data/sample01.fastq.gz" \
  --sample_name "sample01" \
  --gambit_db "/path/to/gambit/db" \
  --outdir "results"
```

---

### 2. Multi-sample Mode (Samplesheet)

To process a batch of samples in parallel, provide a CSV samplesheet formatted as follows:

```csv
sample_id,fastq_1,fastq_2
ERR4422341,test/data/ERR4422341_1.fastq.gz,test/data/ERR4422341_2.fastq.gz
ERR4422727,test/data/ERR4422727_1.fastq.gz,test/data/ERR4422727_2.fastq.gz
```

Run the pipeline:

```bash
nextflow run main.nf \
  -profile docker \
  --samplesheet_csv "samplesheet.csv" \
  --gambit_db "/path/to/gambit/db" \
  --outdir "results"
```

---

### 3. Execution Profiles

Pre-configured profiles in `nextflow.config` can be combined using comma-separated arguments:

- `docker`: Executes all processes inside official BioContainers and custom Docker images (`mmcj/check_reads:1.1`).
- `conda`: Automatically configures environments from YAML specifications in `envs/`.
- `test`: Executes pipeline with bundled test data (`test/data/sample_sheet.csv`).
- `using_gambit`: Activates GAMBIT-based taxonomy with bundled test reference database.
- `using_kraken2`: Activates Kraken2 + Bracken taxonomy with bundled test reference database.

#### Example Commands

Run test dataset with Docker and GAMBIT:
```bash
nextflow run main.nf -profile test,docker,using_gambit
```

Run test dataset with Conda and Kraken2:
```bash
nextflow run main.nf -profile test,conda,using_kraken2
```

Run with Bakta annotation enabled instead of Prokka:
```bash
nextflow run main.nf \
  -profile docker \
  --samplesheet_csv "samplesheet.csv" \
  --gambit_db "/path/to/gambit/db" \
  --use_bakta \
  --bakta_db_type "light" \
  --outdir "results"
```

---

## Parameters Reference

### Input / Output

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--reads` | String | `null` | File path pattern for single or paired-end FASTQ reads (e.g. `"data/*_{1,2}.fastq.gz"`) |
| `--samplesheet_csv` | String | `null` | Path to CSV samplesheet with header `sample_id,fastq_1,fastq_2` |
| `--sample_name` | String | `'sample01'` | Identifier used when running in single-sample mode with `--reads` |
| `--outdir` | String | `'results'` | Output directory where final results will be organized |

### Pre-Assembly Read Screening & Quality Gating

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--min_reads` | Integer | `50` | Minimum total number of reads required per sample |
| `--min_basepairs` | Integer | `2000000` | Minimum total base pairs required (2 Mb default) |
| `--min_proportion` | Integer | `40` | Minimum percentage of total bases in R1 or R2 (flags PE imbalance) |
| `--min_genome_length` | Integer | `100000` | Minimum estimated genome length (bp) via Mash (100 kb default) |
| `--max_genome_length` | Integer | `17000000` | Maximum estimated genome length (bp) via Mash (17 Mb default) |
| `--min_coverage` | Integer | `10` | Minimum estimated sequencing coverage depth (10x default) |

### Taxonomic Classification

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--use_gambit` | Boolean | `true` | Enable GAMBIT species identification |
| `--gambit_db` | String | `null` | Directory path to GAMBIT reference database |
| `--use_kraken2` | Boolean | `false` | Enable Kraken2 taxonomic classification and Bracken re-estimation |
| `--kraken_db` | String | `null` | Path to Kraken2 reference database directory |
| `--bracken_class_level` | String | `'S'` | Taxonomic rank for Bracken abundance estimation (`D`, `P`, `C`, `O`, `F`, `G`, `S`) |
| `--kraken_read_len` | String | `'100'` | Read length parameter for Bracken estimation |
| `--bracken_threshold` | String | `'0'` | Minimum read threshold required prior to abundance re-estimation |

### Genome Annotation

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--use_bakta` | Boolean | `false` | Enable Bakta annotation (defaults to Prokka if `false`) |
| `--bakta_db` | String | `null` | Path to existing Bakta database directory |
| `--bakta_db_type` | String | `'full'` | Database type to auto-download if `--bakta_db` is omitted (`'light'` or `'full'`) |

### Assembly, QC & Typing

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--busco_lineage` | String | `'bacteria_odb12'` | BUSCO lineage dataset for genome completeness assessment |
| `--mlst_scheme` | String | `null` | Specific MLST scheme (e.g. `ecoli`, `saureus`). Automatically inferred if omitted |

### Compute Resources

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--max_cpus` | Integer | System cores - 2 | Maximum CPUs allocated to high-resource processes |
| `--max_mem` | String | `'28.GB'` | Maximum memory allocated to high-resource processes |

---

## Output Directory Structure

Pipeline results are organized hierarchically by sample identifier, along with a shared phylogeny directory:

```text
results/
├── <sample_id>/
│   ├── ReadScreen/
│   │   ├── <sample_id>_read_screen.tsv         # Pre-assembly QC screening metrics & PASS/FAIL flag
│   │   ├── <sample_id>_1.fastq.gz              # Validated forward reads (if passed screening)
│   │   └── <sample_id>_2.fastq.gz              # Validated reverse reads (if passed screening)
│   ├── fastqc/
│   │   ├── <sample_id>_fastqc.html             # Raw read FastQC report
│   │   └── <sample_id>_fastqc.zip              # FastQC data metrics archive
│   ├── fastp/
│   │   ├── <sample_id>_trimmed_R1.fastq.gz     # Filtered & trimmed forward reads
│   │   ├── <sample_id>_trimmed_R2.fastq.gz     # Filtered & trimmed reverse reads
│   │   ├── <sample_id>_report.fastp.html       # fastp visual HTML summary
│   │   └── <sample_id>_report.fastp.json       # fastp machine-readable stats
│   ├── spades/
│   │   ├── <sample_id>.fasta                   # Assembled genome scaffolds
│   │   └── <sample_id>.fastg                   # Assembly graph file
│   ├── assembly_QC/
│   │   ├── <sample_id>_Bandage_Img.jpg         # Visual rendering of assembly graph
│   │   ├── quast/                              # QUAST contig metrics and HTML reports
│   │   └── Busco/                              # BUSCO completeness summaries (JSON/TXT)
│   ├── Annotation/
│   │   ├── Prokka/                             # Prokka annotation files (.gff, .gbk, .faa, .fna)
│   │   └── Bakta/                              # Bakta annotation outputs (if --use_bakta)
│   ├── Gambit/
│   │   └── <sample_id>_gambit_report.csv       # GAMBIT taxonomic predictions
│   ├── kraken2_Bracken/
│   │   ├── <sample_id>_kraken2_report.txt      # Kraken2 hierarchical classification report
│   │   ├── <sample_id>_kraken2_output.txt      # Per-contig classification output
│   │   ├── <sample_id>_bracken_report.txt      # Bracken adjusted abundance report
│   │   └── <sample_id>_bracken_output.txt      # Bracken abundance estimation table
│   ├── mlst/
│   │   └── <sample_id>_mlst_report.csv         # PubMLST sequence typing and allele profiles
│   └── Abricate/
│       ├── <sample_id>_vfdb_report.txt         # Virulence factor screening results
│       └── <sample_id>_plasmidfinder_report.txt# Plasmid replicon typing results
├── phylogeny/
│   └── assembly_sample_sheet.txt               # Aggregated sample manifest for downstream phylogenetic analysis
└── Bakta_DB/                                   # Downloaded Bakta database (if generated)
```

---

## Tools and Citations

If you use BactogenMiner in your research, please cite the tools utilized:

- **Nextflow**: Di Tommaso, P., et al. (2017). Nextflow enables reproducible computational workflows. *Nature Biotechnology*, 35(4), 316–319.
- **fastq-scan**: Brown, D. C. `fastq-scan` (GitHub repository: [https://github.com/rpetit3/fastq-scan](https://github.com/rpetit3/fastq-scan)).
- **Mash**: Ondov, B. D., et al. (2016). Mash: fast genome and metagenome distance estimation using MinHash. *Genome Biology*, 17(1), 132.
- **TheiaProk / PHB Screening Concept**: Theiagen Genomics (GitHub repository: [https://github.com/theiagen/public_health_bioinformatics](https://github.com/theiagen/public_health_bioinformatics)).
- **FastQC**: Andrews, S. (2010). FastQC: a quality control tool for high throughput sequence data.
- **fastp**: Chen, S., et al. (2018). fastp: an ultra-fast all-in-one FASTQ preprocessor. *Bioinformatics*, 34(17), i884–i890.
- **SPAdes**: Prjibelski, A., et al. (2020). Using SPAdes de novo assembler. *Current Protocols in Bioinformatics*, 70(1), e102.
- **Bandage**: Wick, R. R., et al. (2015). Bandage: interactive visualization of de novo genome assemblies. *Bioinformatics*, 31(20), 3350–3352.
- **QUAST**: Gurevich, A., et al. (2013). QUAST: quality assessment tool for genome assemblies. *Bioinformatics*, 29(8), 1072–1075.
- **BUSCO**: Manni, M., et al. (2021). BUSCO: Assessing genomic data quality and completeness. *Current Protocols in Bioinformatics*, 71(1), e104.
- **Prokka**: Seemann, T. (2014). Prokka: rapid prokaryotic genome annotation. *Bioinformatics*, 30(14), 2068–2069.
- **Bakta**: Schwengers, O., et al. (2021). Bakta: rapid and standardized annotation of bacterial genomes via alignment-free sequence identification. *Microbial Genomics*, 7(11), 000685.
- **GAMBIT**: Lumpe, J., et al. (2023). GAMBIT: rapid and accurate bacterial identification from whole-genome sequencing data. *GigaScience*, 12, giad034.
- **Kraken2**: Wood, D. E., et al. (2019). Improved metagenomic analysis with Kraken 2. *Genome Biology*, 20(1), 257.
- **Bracken**: Lu, J., et al. (2017). Bracken: estimating species abundance in metagenomics data. *PeerJ Computer Science*, 3, e104.
- **MLST**: Seemann, T. `mlst` (GitHub repository: [https://github.com/tseemann/mlst](https://github.com/tseemann/mlst)).
- **ABRICATE**: Seemann, T. `abricate` (GitHub repository: [https://github.com/tseemann/abricate](https://github.com/tseemann/abricate)).

---

## Author and Contact

- **MEKA Moise**
- Email: `moise.meka@students.unibe.ch`
- MSc Bioinformatics and Computational Biology, University of Bern