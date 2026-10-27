include { ABRICATE as ABRICATE_AMR  } from "../../modules/abricate.nf"
include { ABRICATE as ABRICATE_VFDB } from "../../modules/abricate.nf"

workflow DEFAULT_SEROTYPING {
    take:
        scafolds_ch // channel: [ val(sample_name), path(scafolds) ]

    main:
        ABRICATE_AMR(scafolds_ch, "resfinder")
        ABRICATE_VFDB(scafolds_ch, "vfdb")

    emit:
        abricate_amr  = ABRICATE_AMR.out.abricate_results
        abricate_vfdb = ABRICATE_VFDB.out.abricate_results
}