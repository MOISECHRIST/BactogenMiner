#!/usr/bin/env nextflow

include { KRAKEN2                           }      from "../modules/kraken2.nf"
include { BRACKEN                           }      from "../modules/Bracken.nf"
include { MLST                              }      from "../modules/mlst.nf"
include { GAMBIT                            }      from "../modules/gambit.nf"
include { ABRICATE as ABRICATE_VFDB         }      from "../modules/abricate.nf"
include { ABRICATE as ABRICATE_AMR          }      from "../modules/abricate.nf"
include { ABRICATE as ABRICATE_PLASMID      }      from "../modules/abricate.nf"
include { ABRICATE as ABRICATE_ECOLI        }      from "../modules/abricate.nf"
include { GET_SPECIES_BRAKEN                }      from "../modules/get_species_braken.nf"
include { GET_SPECIES_GAMBIT                }      from "../modules/get_species_gambit.nf"
include { JOIN_SPECIES_TOOLS                }      from "../modules/join_species_with_tools.nf"
include { KLEBORATE                         }      from "../modules/kleborate.nf"
include { ECTYPER                           }      from "../modules/ectyper.nf"
include { ECTYPER_GET_DB                    }      from "../modules/ectyper_get_db.nf"
include { SEQSERO2                          }      from "../modules/seqsero2.nf"
include { LISSERO                           }      from "../modules/lissero.nf"
include { PASTY                             }      from "../modules/pasty.nf"
include { UNTAR_DB as UNTAR_GAMBIT          }      from "../modules/processing_db"
include { UNTAR_DB as UNTAR_KRAKEN          }      from "../modules/processing_db"
include { PREPARE_DB as PREPARE_ECTYPER_DB  }      from "../modules/processing_db.nf"
include { PREPARE_DB as PREPARE_KRAKEN_DB   }      from "../modules/processing_db.nf"
include { PREPARE_DB as PREPARE_GAMBIT_DB   }      from "../modules/processing_db.nf"
include { GET_DB as GET_GAMBIT_DB           }      from '../modules/processing_db.nf'
include { GET_DB as GET_KRAKEN_DB           }      from '../modules/processing_db.nf'

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
        scafolds_ch  // channel: [ val(sample_name), path(scafolds) ]
        species_name // channel: [ val(sample_name), val(species_name)  ]

    main:
        //Species specific screening genome assemblies
        JOIN_SPECIES_TOOLS(species_name)
        sp_group_ch = JOIN_SPECIES_TOOLS.out.sp_group   
        tools_ch    = JOIN_SPECIES_TOOLS.out.tools      

        routed = scafolds_ch
            .join(sp_group_ch)      
            .join(tools_ch)         
            .branch { sample_name, scafolds, sp_group, tools ->
                kleborate: tools == "kleborate"
                seqsero:   tools == "SeqSero2"
                lissero:   tools == "LisSero"
                pasty:     tools == "Pasty"
                other:     true
            }

        KLEBORATE(
            routed.kleborate.map { sn, sc, sp, tl -> tuple(sn, sc, sp) }
        )
        ecoli_only = routed.kleborate.filter { sn, sc, sp, tl -> sp == "escherichia" }
        ecoli_fa   = ecoli_only.map { sn, sc, sp, tl -> tuple(sn, sc) }

        ABRICATE_ECOLI(ecoli_fa, "ecoli_vf")

        if (params.ectyper_db) {
            ectyper_db_ch = Channel.value(file(params.ectyper_db))
        } else {
            ECTYPER_GET_DB(ecoli_only.first().map { true })
            ectyper_db_ch = ECTYPER_GET_DB.out.ectyper_db
        }

        ECTYPER(ecoli_fa, ectyper_db_ch)

        SEQSERO2(routed.seqsero.map {sn, sc, sp, tl -> tuple(sn, sc)}, "4")//for illumina assembly fasta
        //SEQSERO2(routed.seqsero.map {sn, sc, sp, tl -> tuple(sn, sc)}, "5")//for Long Reads assembly fasta

        LISSERO(routed.lissero.map {sn, sc, sp, tl -> tuple(sn, sc)})
        PASTY(routed.pasty.map {sn, sc, sp, tl -> tuple(sn, sc)})

        ABRICATE_AMR(routed.other.map { sn, sc, sp, tl -> tuple(sn, sc) }, "resfinder")
        ABRICATE_VFDB(routed.other.map { sn, sc, sp, tl -> tuple(sn, sc) }, "vfdb")
        ABRICATE_PLASMID(scafolds_ch, "plasmidfinder")
    
    emit:
        kleborate             = KLEBORATE.out.kleborate_results
        ectyper               = ECTYPER.out.ectyper_results
        seqsero2              = SEQSERO2.out.seqsero2_results
        lissero               = LISSERO.out.lissero_results
        pasty                 = PASTY.out.pasty_results
        abricate_ecoli        = ABRICATE_ECOLI.out.abricate_results
        abricate_amr          = ABRICATE_AMR.out.abricate_results
        abricate_vfdb         = ABRICATE_VFDB.out.abricate_results
        abricate_plasmid      = ABRICATE_PLASMID.out.abricate_results
}