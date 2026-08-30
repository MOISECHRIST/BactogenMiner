include { SINGLE_SAMPLE_PROCESSING }      from "./workflows/single_sample_pipeline.nf"
include { MAKE_ASSEMBLY_SHEET      }      from "./modules/make_assembly_sheet.nf"
workflow {
    main:
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
        SINGLE_SAMPLE_PROCESSING(reads_ch)
        MAKE_ASSEMBLY_SHEET(SINGLE_SAMPLE_PROCESSING.out.scafolds.collect())
}
