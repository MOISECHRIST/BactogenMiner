# BactogenMiner

A reproducible Nextflow (DSL2) workflow for bacterial whole-genome sequencing (WGS) data analysis, covering read preprocessing, de novo assembly, assembly quality assessment, taxonomic classification, and multi-locus sequence typing (MLST).

---

## Table of Contents

- [Overview](#overview)
- [Pipeline Architecture](#pipeline-architecture)
- [Prerequisites and Dependencies](#prerequisites-and-dependencies)
- [Databases](#databases)
- [Installation](#installation)
- [Usage](#usage)
  - [Single-sample Mode](#single-sample-mode)
  - [Multi-sample Mode (Samplesheet)](#multi-sample-mode-samplesheet)
  - [Execution Profiles](#execution-profiles)
- [Parameters](#parameters)
- [Output Structure](#output-structure)
- [Tools and Citations](#tools-and-citations)
- [Author and Contact](#author-and-contact)

---

## Overview

This pipeline automates the downstream processing of raw bacterial paired-end or single-end Illumina sequencing reads:

1. **Read QC and Trimming**: Quality evaluation with FastQC and automated adapter trimming/filtering with fastp.
2. **De Novo Assembly**: Genome reconstruction using SPAdes.
3. **Assembly Quality Control**:
   - Contig metrics and assembly statistics via QUAST.
   - Assembly graph visualization via Bandage.
   - Genome completeness assessment via BUSCO (`bacteria_odb12` lineage).
4. **Species Identification and Taxonomy**:
   - Rapid genomic distance classification via GAMBIT.
   - Optional k-mer based taxonomic classification with Kraken2 and abundance re-estimation with Bracken.
   - Optional MinHash sketch and distance estimation via Mash.
5. **Sequence Typing**: Multi-locus sequence typing via MLST (PubMLST schemes).

---

## Pipeline Architecture

```
Raw FastQ Reads
      │
      ├──> FastQC (Raw Quality Assessment)
      └──> fastp  (Adapter Trimming & Quality Filtering)
             │
             └──> Cleaned Reads
                    │
                    └──> SPAdes (De Novo Assembly)
                           │
                           ├──> Bandage (Assembly Graph Image)
                           ├──> QUAST   (Assembly Metrics)
                           ├──> BUSCO   (Lineage Completeness)
                           ├──> GAMBIT  (Species Classification)
                           ├──> Kraken2 / Bracken (Taxonomy & Abundance - Optional)
                           └──> MLST    (Sequence Typing)
```

---

## Prerequisites and Dependencies

- **Nextflow** (>= 22.10.0)
- **Container or Environment Manager** :
  - **Docker**
  - **Conda** 

---

## Databases

Before running taxonomic classification, download the required reference databases:

- **GAMBIT Database**: Available on the [GAMBIT Documentation](https://gambit-genomics.readthedocs.io/en/latest/databases.html).
- **Kraken2 / Bracken Database**: Available on the [Ben Langmead AWS Indexes](https://benlangmead.github.io/aws-indexes/k2).
- **Mash Database**: Available on the [Mash Documentation](https://mash.readthedocs.io/en/latest/tutorials.html).
- **BUSCO Lineage**: Automatically downloaded upon execution or configured locally (default lineage: `bacteria_odb12`).

---

## Installation

Clone the repository to your local system:

```bash
git clone https://github.com/MOISECHRIST/bacteria_typing-phylogeny.git
cd bacteria_typing-phylogeny
```

---

## Usage

### Single-sample Mode

To analyze a single sample with paired-end reads:

```bash
nextflow run main.nf \
  -profile conda \
  --reads "data/sample01_{1,2}.fastq.gz" \
  --sample_name "sample01" \
  --gambit_db "/path/to/gambit/db" \
  --outdir "results"
```

### Multi-sample Mode (Samplesheet)

To analyze multiple samples simultaneously, create a CSV samplesheet with the header `sample_id,fastq_1,fastq_2`:

Example `samplesheet.csv`:
```csv
sample_id,fastq_1,fastq_2
ERR4422341,test/data/ERR4422341_1.fastq.gz,test/data/ERR4422341_2.fastq.gz
ERR4422727,test/data/ERR4422727_1.fastq.gz,test/data/ERR4422727_2.fastq.gz
```

Run the pipeline:

```bash
nextflow run main.nf \
  -profile conda \
  --samplesheet_csv "samplesheet.csv" \
  --gambit_db "/path/to/gambit/db" \
  --outdir "results"
```

### Execution Profiles

The pipeline provides pre-configured execution profiles via the `-profile` option:

- `docker`: Executes pipeline processes inside BioContainers images.
- `conda`: Manages process environments using Conda yaml recipes in `envs/`.
- `test`: Executes the pipeline with bundled test data (`test/data/sample_sheet.csv`).
- `using_gambit`: Activates GAMBIT-based taxonomic classification.
- `using_kraken2`: Activates Kraken2 and Bracken classification workflows.

Example running the built-in test with Docker and GAMBIT:

```bash
nextflow run main.nf -profile test,docker,using_gambit
```

---

## Parameters

### Input and Output

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--reads` | String | `null` | Path pattern to paired-end/single-end fastq files |
| `--samplesheet_csv` | String | `null` | Path to CSV samplesheet (`sample_id,fastq_1,fastq_2`) |
| `--sample_name` | String | `sample01` | Sample identifier for single sample run |
| `--outdir` | String | `results` | Output directory for pipeline results |

### Taxonomy and Classification

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--use_gambit` | Boolean | `true` | Enable GAMBIT species identification |
| `--gambit_db` | String | `null` | Path to GAMBIT reference database |
| `--use_kraken2` | Boolean | `false` | Enable Kraken2 and Bracken classification |
| `--kraken_db` | String | `null` | Path to Kraken2 database |
| `--bracken_class_level` | String | `S` | Taxonomic level for Bracken estimation (`D`, `P`, `C`, `O`, `F`, `G`, `S`) |
| `--kraken_read_len` | String | `100` | Read length for Bracken estimation |
| `--bracken_threshold` | String | `0` | Minimum read threshold for abundance re-estimation |

### Assembly, QC and Typing

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--busco_lineage` | String | `bacteria_odb12` | BUSCO dataset lineage for completeness estimation |
| `--mlst_scheme` | String | `null` | Specific MLST scheme (e.g. `ecoli`, `saureus`). Automatically inferred if omitted |

### Resource Allocation

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--max_cpus` | Integer | System available - 2 | Maximum CPUs allocated to high-resource processes |
| `--max_mem` | String | `28.GB` | Maximum memory limit for high-resource processes |

---

## Output Structure

Results are organized by sample identifier under the specified output directory (`results/`):

```
results/
└── <sample_id>/
    ├── fastqc/
    │   ├── <sample_id>_fastqc.html
    │   └── <sample_id>_fastqc.zip
    ├── fastp/
    │   ├── <sample_id>_trimmed_R1.fastq.gz
    │   ├── <sample_id>_trimmed_R2.fastq.gz
    │   ├── <sample_id>_report.fastp.html
    │   └── <sample_id>_report.fastp.json
    ├── spades/
    │   ├── <sample_id>.fasta              <- Assembled scaffolds
    │   └── <sample_id>.fastg              <- Assembly graph
    ├── assembly_QC/
    │   ├── <sample_id>_Bandage_Img.jpg    <- Assembly graph visual
    │   ├── quast/                         <- QUAST metrics and HTML reports
    │   └── Busco/                         <- BUSCO summaries (JSON / TXT)
    ├── Gambit/
    │   └── <sample_id>_gambit_report.csv  <- GAMBIT species classification
    ├── kraken2_Bracken/
    │   ├── <sample_id>_kraken2_report.txt
    │   ├── <sample_id>_kraken2_output.txt
    │   ├── <sample_id>_bracken_report.txt
    │   └── <sample_id>_bracken_output.txt
    └── mlst/
        └── <sample_id>_mlst_report.csv    <- MLST sequence typing output
```

---

## Tools and Citations

- **FastQC**: Andrews, S. (2010). FastQC: a quality control tool for high throughput sequence data.
- **fastp**: Chen, S., et al. (2018). fastp: an ultra-fast all-in-one FASTQ preprocessor. *Bioinformatics*, 34(17), i884–i890.
- **SPAdes**: Prjibelski, A., et al. (2020). Using SPAdes de novo assembler. *Current Protocols in Bioinformatics*, 70(1), e102.
- **Bandage**: Wick, R. R., et al. (2015). Bandage: interactive visualization of de novo genome assemblies. *Bioinformatics*, 31(20), 3350-3352.
- **QUAST**: Gurevich, A., et al. (2013). QUAST: quality assessment tool for genome assemblies. *Bioinformatics*, 29(8), 1072-1075.
- **BUSCO**: Manni, M., et al. (2021). BUSCO: Assessing genomic data quality and completeness. *Current Protocols in Bioinformatics*, 71(1), e104.
- **GAMBIT**: Lumpe, J., et al. (2023). GAMBIT: rapid and accurate bacterial identification from whole-genome sequencing data. *GigaScience*, 12, giad034.
- **Kraken2**: Wood, D. E., et al. (2019). Improved metagenomic analysis with Kraken 2. *Genome Biology*, 20(1), 257.
- **Bracken**: Lu, J., et al. (2017). Bracken: estimating species abundance in metagenomics data. *PeerJ Computer Science*, 3, e104.
- **Mash**: Ondov, B. D., et al. (2016). Mash: fast genome and metagenome distance estimation using MinHash. *Genome Biology*, 17(1), 132.
- **MLST**: Seemann, T. mlst (GitHub repository: https://github.com/tseemann/mlst).

---

## Author and Contact

- **MEKA Moise**
- Email: `moise.meka@students.unibe.ch`
- MSc Bioinformatics and Computational Biology, University of Bern