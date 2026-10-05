include { SAMPLE_QC  as SHORT_READS_SAMPLE_QC                     } from "./short_reads_workflows.nf"
include { ASSEMBLY_ANNOTATION as SHORT_READS_ASSEMBLY_ANNOTATION  } from "./short_reads_workflows.nf"
include { SPECIES_CLASSIFICATION                                  } from "./post_assembly_workflow.nf"
include { SEROTYPING                                              } from "./post_assembly_workflow.nf"


workflow SHORT_READS_SINGLE_PROCESSING {
    take:
        reads_ch  // channel: [ val(sample_name), [ path(reads) ] ]

    main:
        // 01 - Quality Control 
        SHORT_READS_SAMPLE_QC(reads_ch)
        // 02 - Assembly and Annotation 
        SHORT_READS_ASSEMBLY_ANNOTATION(SHORT_READS_SAMPLE_QC.out.cleaned_reads)
        SHORT_READS_ASSEMBLY_ANNOTATION.out.assembly_status
                .map { sample_name, status -> "${sample_name}\t${status}" }
                .collectFile(
                    name: 'assembly_status.tsv',
                    storeDir: params.outdir,
                    seed: "sample_name\tassembly_status",
                    newLine: true,
                    sort: true
                )
        // 03 - Species Classification
        SPECIES_CLASSIFICATION(SHORT_READS_ASSEMBLY_ANNOTATION.out.scafolds)
        // 04 - Serotyping
        SEROTYPING(SHORT_READS_ASSEMBLY_ANNOTATION.out.scafolds, SPECIES_CLASSIFICATION.out.species_name)
    emit:
        scafolds = SHORT_READS_ASSEMBLY_ANNOTATION.out.scafolds
}
