#!/usr/bin/env nextflow

include { KRAKEN2                           } from "../modules/kraken2.nf"
include { BRACKEN                           } from "../modules/Bracken.nf"
include { MLST                              } from "../modules/mlst.nf"
include { GAMBIT                            } from "../modules/gambit.nf"
include { GET_SPECIES_BRAKEN                } from "../modules/get_species_braken.nf"
include { GET_SPECIES_GAMBIT                } from "../modules/get_species_gambit.nf"
include { UNTAR_DB as UNTAR_GAMBIT          } from "../modules/processing_db"
include { UNTAR_DB as UNTAR_KRAKEN          } from "../modules/processing_db"
include { PREPARE_DB as PREPARE_ECTYPER_DB  } from "../modules/processing_db.nf"
include { PREPARE_DB as PREPARE_KRAKEN_DB   } from "../modules/processing_db.nf"
include { PREPARE_DB as PREPARE_GAMBIT_DB   } from "../modules/processing_db.nf"
include { GET_DB as GET_GAMBIT_DB           } from '../modules/processing_db.nf'
include { GET_DB as GET_KRAKEN_DB           } from '../modules/processing_db.nf'
include { JOIN_SPECIES_TOOLS                } from "../modules/join_species_with_tools.nf"
include { ABRICATE as ABRICATE_PLASMID      } from "../modules/abricate.nf"
include { KLEB_ECOLI_SEROTYPING             } from "../sub-workflows/serotyping/kleb_ecoli_serotyping.nf"
include { SALMONELLA_SEROTYPING             } from "../sub-workflows/serotyping/salmonella_serotyping.nf"
include { LISTERIA_SEROTYPING               } from "../sub-workflows/serotyping/listeria_serotyping.nf"
include { PSEUDOMONAS_SEROTYPING            } from "../sub-workflows/serotyping/pseudomonas_serotyping.nf"
include { DEFAULT_SEROTYPING                } from "../sub-workflows/serotyping/default_serotyping.nf"

workflow SPECIES_CLASSIFICATION {
    take:
        scafolds_ch  // channel: [ val(sample_name), path(scafolds) ]

    main:
        // Sequence typing
        MLST(scafolds_ch)

        // Species classification 
        if (params.use_gambit) {
            if (params.gambit_db) {
                PREPARE_GAMBIT_DB(Channel.value(params.gambit_db))
                gambit_db_ch = PREPARE_GAMBIT_DB.out.db
            } else {
                GET_GAMBIT_DB(Channel.of(tuple('gambit_db', params.gambit_db_url, params.gambit_db_md5)))
                UNTAR_GAMBIT(GET_GAMBIT_DB.out.db)
                gambit_db_ch = UNTAR_GAMBIT.out.db
            }
            
            GAMBIT(scafolds_ch, gambit_db_ch)
            GET_SPECIES_GAMBIT(GAMBIT.out.gambit_report)
            species_name = GET_SPECIES_GAMBIT.out.species_name
            gambit_report_ch  = GAMBIT.out.gambit_report
            bracken_report_ch = channel.empty()
            bracken_output_ch = channel.empty()
            kraken2_output_ch = channel.empty()
            kraken2_report_ch = channel.empty()
        } else {
            if (params.kraken_db) {
                PREPARE_KRAKEN_DB(Channel.value(params.kraken_db))
                kraken_db_ch = PREPARE_KRAKEN_DB.out.db
            } else {
                GET_KRAKEN_DB(Channel.of(tuple('k2_standard_08_GB_20260626', params.kraken_db_url, params.kraken_db_md5)))
                UNTAR_KRAKEN(GET_KRAKEN_DB.out.db)
                kraken_db_ch = UNTAR_KRAKEN.out.db
            }

            KRAKEN2(scafolds_ch, kraken_db_ch)
            BRACKEN(KRAKEN2.out.kraken_out, kraken_db_ch)
            GET_SPECIES_BRAKEN(BRACKEN.out.bracken_output)
            species_name = GET_SPECIES_BRAKEN.out.species_name
            gambit_report_ch  = channel.empty()
            bracken_report_ch = BRACKEN.out.bracken_report
            bracken_output_ch = BRACKEN.out.bracken_output
            kraken2_output_ch = KRAKEN2.out.kraken2_output
            kraken2_report_ch = KRAKEN2.out.kraken2_report
        }
    
    emit:
        species_name        = species_name
        sequence_type       = MLST.out.mlst_output
        gambit_report       = gambit_report_ch
        bracken_report      = bracken_report_ch
        bracken_output      = bracken_output_ch
        kraken2_output      = kraken2_output_ch
        kraken2_report      = kraken2_report_ch
}


workflow SEROTYPING {
    take:
        scafolds_ch   // channel: [ val(sample_name), path(scafolds) ]
        species_name  // channel: [ val(sample_name), val(species_name) ]

    main:
        JOIN_SPECIES_TOOLS(species_name)

        routed = scafolds_ch
            .join(JOIN_SPECIES_TOOLS.out.sp_group)
            .join(JOIN_SPECIES_TOOLS.out.tools)
            .branch { sn, sc, sp, tl ->
                kleborate: tl == "kleborate"
                seqsero:   tl == "SeqSero2"
                lissero:   tl == "LisSero"
                pasty:     tl == "Pasty"
                other:     true
            }

        KLEB_ECOLI_SEROTYPING(routed.kleborate.map { sn, sc, sp, tl -> tuple(sn, sc, sp) })
        SALMONELLA_SEROTYPING(routed.seqsero.map   { sn, sc, sp, tl -> tuple(sn, sc) })
        LISTERIA_SEROTYPING(routed.lissero.map     { sn, sc, sp, tl -> tuple(sn, sc) })
        PSEUDOMONAS_SEROTYPING(routed.pasty.map    { sn, sc, sp, tl -> tuple(sn, sc) })
        DEFAULT_SEROTYPING(routed.other.map        { sn, sc, sp, tl -> tuple(sn, sc) })

        // Universal: all samples
        ABRICATE_PLASMID(scafolds_ch, "plasmidfinder")

    emit:
        //KpSC, KoSC, Ecoli & Shigella
        kleborate              = KLEB_ECOLI_SEROTYPING.out.kleborate
        ectyper                = KLEB_ECOLI_SEROTYPING.out.ectyper
        abricate_ecoli         = KLEB_ECOLI_SEROTYPING.out.abricate_ecoli
        //Salmonella 
        seqsero2               = SALMONELLA_SEROTYPING.out.seqsero2
        abricate_vfdb_seqsero2 = SALMONELLA_SEROTYPING.out.abricate_vfdb_seqsero2
        abricate_amr_seqsero2  = SALMONELLA_SEROTYPING.out.abricate_amr_seqsero2
        //Listeria 
        lissero                = LISTERIA_SEROTYPING.out.lissero
        abricate_vfdb_lissero  = LISTERIA_SEROTYPING.out.abricate_vfdb_lissero
        abricate_amr_lissero   = LISTERIA_SEROTYPING.out.abricate_amr_lissero
        //Pseudomonas
        pasty                  = PSEUDOMONAS_SEROTYPING.out.pasty
        abricate_vfdb_pasty    = PSEUDOMONAS_SEROTYPING.out.abricate_vfdb_pasty
        abricate_amr_pasty     = PSEUDOMONAS_SEROTYPING.out.abricate_amr_pasty
        //Other species
        abricate_amr_default   = DEFAULT_SEROTYPING.out.abricate_amr
        abricate_vfdb_default  = DEFAULT_SEROTYPING.out.abricate_vfdb
        //Plasmids
        abricate_plasmid_all   = ABRICATE_PLASMID.out.abricate_results
}