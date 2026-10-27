include { KLEBORATE                  } from "../../modules/kleborate.nf"
include { ECTYPER                    } from "../../modules/ectyper.nf"
include { ECTYPER_GET_DB             } from "../../modules/ectyper_get_db.nf"
include { ABRICATE as ABRICATE_ECOLI } from "../../modules/abricate.nf"

workflow KLEB_ECOLI_SEROTYPING {
    take:
        kleb_ch  // channel: [ val(sample_name), path(scafolds), val(sp_group) ]

    main:
        KLEBORATE(kleb_ch)

        ecoli_only = kleb_ch.filter { sn, sc, sp -> sp == "escherichia" }
        ecoli_fa   = ecoli_only.map { sn, sc, sp -> tuple(sn, sc) }

        ABRICATE_ECOLI(ecoli_fa, "ecoli_vf")

        if (params.ectyper_db) {
            ectyper_db_ch = Channel.value(file(params.ectyper_db))
        } else {
            ECTYPER_GET_DB(ecoli_only.first().map { true })
            ectyper_db_ch = ECTYPER_GET_DB.out.ectyper_db
        }
        ECTYPER(ecoli_fa, ectyper_db_ch)

    emit:
        kleborate      = KLEBORATE.out.kleborate_results
        ectyper        = ECTYPER.out.ectyper_results
        abricate_ecoli = ABRICATE_ECOLI.out.abricate_results
}