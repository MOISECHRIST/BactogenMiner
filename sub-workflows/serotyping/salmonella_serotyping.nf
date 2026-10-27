include { SEQSERO2                  } from "../../modules/seqsero2.nf"
include { ABRICATE as ABRICATE_AMR  } from "../../modules/abricate.nf"
include { ABRICATE as ABRICATE_VFDB } from "../../modules/abricate.nf"

workflow SALMONELLA_SEROTYPING {
    take:
        scafolds_ch  // channel: [ val(sample_name), path(scafolds) ]

    main:
        SEQSERO2(scafolds_ch, "4")   // "5" for long-read assemblies
        ABRICATE_VFDB(scafolds_ch, "vfdb")
        ABRICATE_AMR(scafolds_ch, "resfinder")


    emit:
        seqsero2 = SEQSERO2.out.seqsero2_results
        abricate_vfdb_seqsero2 = ABRICATE_VFDB.out.abricate_results
        abricate_amr_seqsero2  = ABRICATE_AMR.out.abricate_results
}