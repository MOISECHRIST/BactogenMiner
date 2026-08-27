include { SINGLE_SAMPLE_PROCESSING } from "./workflows/single_sample_pipeline.nf"

workflow {
    main:
        if (params.reads) {
            reads_ch = channel.fromFilePairs(params.reads, checkIfExists: true)
                .map { prefix, files -> [params.sample_name, files] }
            SINGLE_SAMPLE_PROCESSING(reads_ch)
        } else if (params.samplesheet_csv) {
            reads_ch = channel.fromPath(params.samplesheet_csv, checkIfExists: true)
                .splitCsv(header: true)
                .map { row -> [row.sample_id, [file(row.fastq_1), file(row.fastq_2)]] }
            
            SINGLE_SAMPLE_PROCESSING(reads_ch)
        } else {
            error "Please provide either --reads or --samplesheet_csv"
        }
}
