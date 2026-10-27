include { PASTY } from "../../modules/pasty.nf"
include { ABRICATE as ABRICATE_AMR  } from "../../modules/abricate.nf"
include { ABRICATE as ABRICATE_VFDB } from "../../modules/abricate.nf"

workflow PSEUDOMONAS_SEROTYPING {
    take:
        scafolds_ch // channel: [ val(sample_name), path(scafolds) ]

    main:
        PASTY(scafolds_ch)
        ABRICATE_VFDB(scafolds_ch, "vfdb")
        ABRICATE_AMR(scafolds_ch, "resfinder")

    emit:
        pasty = PASTY.out.pasty_results
        abricate_vfdb_pasty = ABRICATE_VFDB.out.abricate_results
        abricate_amr_pasty  = ABRICATE_AMR.out.abricate_results
}