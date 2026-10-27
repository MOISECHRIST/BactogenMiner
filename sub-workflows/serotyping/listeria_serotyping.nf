include { LISSERO                   } from "../../modules/lissero.nf"
include { ABRICATE as ABRICATE_VFDB } from "../../modules/abricate.nf"
include { ABRICATE as ABRICATE_AMR  } from "../../modules/abricate.nf"

workflow LISTERIA_SEROTYPING {
    take:
        scafolds_ch // channel: [ val(sample_name), path(scafolds) ]

    main:
        LISSERO(scafolds_ch)
        ABRICATE_VFDB(scafolds_ch, "vfdb")
        ABRICATE_AMR(scafolds_ch, "resfinder")

    emit:
        lissero               = LISSERO.out.lissero_results
        abricate_vfdb_lissero = ABRICATE_VFDB.out.abricate_results
        abricate_amr_lissero  = ABRICATE_AMR.out.abricate_results
}