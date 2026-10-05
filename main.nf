#!/usr/bin/env nextflow

include { SHORT_READS_SINGLE_PROCESSING }      from "./workflows/single_sample_pipeline.nf"
include { MAKE_ASSEMBLY_SHEET      }      from "./modules/make_assembly_sheet.nf"

def startMessage(){
        def c_reset  = "\033[0m"
        def c_dim    = "\033[2m"
        def c_bold   = "\033[1m"
        def c_cyan   = "\033[36m"
        def c_yellow = "\033[33m"
        def c_green  = "\033[32m"
    log.info """
    ${c_bold}================================================================================${c_reset}
    ${c_cyan}${c_bold}  B A C T O G E N M I N E R   v${workflow.manifest.version}${c_reset}
      Bacterial WGS Typing, Annotation & Phylogeny Preparation Pipeline
    ${c_bold}================================================================================${c_reset}
    """
}

def helpMessage() {
        def c_reset  = "\033[0m"
        def c_dim    = "\033[2m"
        def c_bold   = "\033[1m"
        def c_cyan   = "\033[36m"
        def c_yellow = "\033[33m"
        def c_green  = "\033[32m"

        log.info """
    ${c_yellow}Usage:${c_reset}
      nextflow run main.nf -profile <docker|conda>[,test,using_gambit|using_kraken2] [options]

      Input: provide EITHER --reads (one sample) OR --samplesheet_csv (one or many samples).
      Taxonomy: enable at least one of --use_gambit / --use_kraken2, otherwise samples
      cannot be routed to species-specific tools (Kleborate, ECTyper, SeqSero2, LisSero, Pasty).

    ${c_yellow}Examples:${c_reset}
      ${c_dim}# 1. Single sample (paired-end) with GAMBIT classification${c_reset}
      nextflow run main.nf -profile docker \\
          --reads "data/sample01_{1,2}.fastq.gz" \\
          --sample_name "sample01" \\
          --use_gambit --gambit_db "/path/to/gambit/db"

      ${c_dim}# 2. Multi-sample batch via CSV samplesheet${c_reset}
      nextflow run main.nf -profile docker \\
          --samplesheet_csv "samplesheet.csv" \\
          --use_gambit --gambit_db "/path/to/gambit/db"

      ${c_dim}# 3. Kraken2/Bracken classification and Bakta annotation${c_reset}
      nextflow run main.nf -profile docker \\
          --samplesheet_csv "samplesheet.csv" \\
          --use_kraken2 --kraken_db "/path/to/k2_db" \\
          --use_bakta --bakta_db_type "light"

      ${c_dim}# 4. Bundled test dataset${c_reset}
      nextflow run main.nf -profile test,docker,using_gambit

    --------------------------------------------------------------------------------

    ${c_green}General:${c_reset}
      --help                 Print this help page and exit
      -resume                (Nextflow option) Resume from cached results

    ${c_green}Input / Output Options:${c_reset}
      --reads                Glob pattern for the FASTQ files of ONE sample (single or paired-end)
      --samplesheet_csv      CSV samplesheet for one or many samples (header: sample_id,fastq_1,fastq_2;
                             leave fastq_2 empty for single-end reads)
      --sample_name          Sample identifier used with --reads [default: '${params.sample_name}']
      --outdir               Output results directory [default: '${params.outdir}']

    ${c_green}Pre-Assembly Read Screening & Gating (CHECK_READS):${c_reset}
      --min_reads            Minimum total read count [default: ${params.min_reads}]
      --min_basepairs        Minimum total sequenced base pairs [default: ${params.min_basepairs}]
      --min_proportion       Minimum percentage of bases in R1/R2 [default: ${params.min_proportion}%]
      --min_genome_length    Minimum estimated genome length in bp (Mash) [default: ${params.min_genome_length}]
      --max_genome_length    Maximum estimated genome length in bp (Mash) [default: ${params.max_genome_length}]
      --min_coverage         Minimum estimated sequencing depth [default: ${params.min_coverage}x]

    ${c_green}Taxonomic Classification:${c_reset}
      --use_gambit           Enable GAMBIT species identification [default: ${params.use_gambit}]
      --gambit_db            Path to GAMBIT reference database directory [default: ${params.gambit_db}]
      --use_kraken2          Enable Kraken2 + Bracken classification [default: ${params.use_kraken2}]
      --kraken_db            Path to Kraken2 reference database directory [default: ${params.kraken_db}]
      --bracken_class_level  Taxonomic rank for Bracken (D|P|C|O|F|G|S) [default: '${params.bracken_class_level}']
      --kraken_read_len      Read length for Bracken estimation [default: '${params.kraken_read_len}']
      --bracken_threshold    Minimum read threshold before re-estimation [default: '${params.bracken_threshold}']

    ${c_green}Genome Annotation:${c_reset}
      --use_bakta            Use Bakta instead of the default Prokka [default: ${params.use_bakta}]
      --bakta_db             Path to a pre-downloaded Bakta database [default: ${params.bakta_db}]
      --bakta_db_type        Bakta database to download if --bakta_db is omitted ('light'|'full') [default: '${params.bakta_db_type}']

    ${c_green}Assembly QC & Typing:${c_reset}
      --busco_lineage        BUSCO lineage dataset for completeness [default: '${params.busco_lineage}']
      --mlst_scheme          PubMLST scheme name [default: ${params.mlst_scheme ?: 'auto-inferred'}]
      --ectyper_db           Path to a pre-downloaded ECTyper database (Escherichia/Shigella) [default: ${params.ectyper_db ?: 'auto'}]

    ${c_green}Compute Resources:${c_reset}
      --max_cpus             Maximum CPU cores per task [default: ${params.max_cpus}]
      --max_mem              Maximum memory per task [default: '${params.max_mem}']

    ${c_green}Execution Profiles (-profile, comma-separated):${c_reset}
      docker                 Run processes inside Docker containers
      conda                  Run processes using local Conda environments
      test                   Use the bundled test samplesheet (test/data/sample_sheet.csv)
      using_gambit           GAMBIT with the bundled test reference database
      using_kraken2          Kraken2 with the bundled test reference database

    ${c_bold}================================================================================${c_reset}
    """
    }

workflow {
    main:
        startMessage()
        if (params.help) {
                helpMessage()
                exit 0
            }

        if (params.reads) {
            reads_ch = channel.fromFilePairs(params.reads, checkIfExists: true)
                .map { prefix, files -> [params.sample_name, files] }
        } else if (params.samplesheet_csv) {
            reads_ch = channel.fromPath(params.samplesheet_csv, checkIfExists: true)
                .splitCsv(header: true)
                .map { row ->
                        def files = row.fastq_2 ? [file(row.fastq_1), file(row.fastq_2)] : [file(row.fastq_1)]
                        [row.sample_id, files]
                    }
        } else {
            error "Please provide either --reads or --samplesheet_csv"
        }
        SHORT_READS_SINGLE_PROCESSING(reads_ch)
        MAKE_ASSEMBLY_SHEET(SHORT_READS_SINGLE_PROCESSING.out.scafolds.collect(flat: false))
}
