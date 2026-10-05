#!/usr/bin/env nextflow

include { KRAKEN2                      }      from "../modules/kraken2.nf"
include { BRACKEN                      }      from "../modules/Bracken.nf"
include { MLST                         }      from "../modules/mlst.nf"
include { GAMBIT                       }      from "../modules/gambit.nf"
include { ABRICATE as ABRICATE_VFDB    }      from "../modules/abricate.nf"
include { ABRICATE as ABRICATE_AMR     }      from "../modules/abricate.nf"
include { ABRICATE as ABRICATE_PLASMID }      from "../modules/abricate.nf"
include { ABRICATE as ABRICATE_ECOLI   }      from "../modules/abricate.nf"
include { GET_SPECIES_BRAKEN           }      from "../modules/get_species_braken.nf"
include { GET_SPECIES_GAMBIT           }      from "../modules/get_species_gambit.nf"
include { JOIN_SPECIES_TOOLS           }      from "../modules/join_species_with_tools.nf"
include { KLEBORATE                    }      from "../modules/kleborate.nf"
include { ECTYPER                      }      from "../modules/ectyper.nf"
include { ECTYPER_GET_DB               }      from "../modules/ectyper_get_db.nf"
include { SEQSERO2                     }      from "../modules/seqsero2.nf"
include { LISSERO                      }      from "../modules/lissero.nf"
include { PASTY                        }      from "../modules/pasty.nf"


workflow SPECIES_CLASSIFICATION {
    take:
        scafolds_ch  // channel: [ val(sample_name), path(scafolds) ]

    main:
        // Sequence typing
        MLST(scafolds_ch)

        // Species classification 
        if (params.use_gambit) {
            GAMBIT(scafolds_ch, file(params.gambit_db))
            GET_SPECIES_GAMBIT(GAMBIT.out.gambit_report)
            species_name = GET_SPECIES_GAMBIT.out.species_name
            gambit_report_ch  = GAMBIT.out.gambit_report
            bracken_report_ch = channel.empty()
            bracken_output_ch = channel.empty()
            kraken2_output_ch = channel.empty()
            kraken2_report_ch = channel.empty()
        } else {
            KRAKEN2(scafolds_ch, file(params.kraken_db))
            BRACKEN(KRAKEN2.out.kraken_out, file(params.kraken_db))
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
        ABRICATE_ECOLI(ecoli_only.map { sn, sc, sp, tl -> tuple(sn, sc) }, "ecoli_vf")
        if(params.ectyper_db){
            ECTYPER(ecoli_only.map { sn, sc, sp, tl -> tuple(sn, sc) }, file(params.ectyper_db))
        } else {
            ECTYPER_GET_DB()
            ECTYPER(ecoli_only.map { sn, sc, sp, tl -> tuple(sn, sc) }, ECTYPER_GET_DB.out.ectyper_db)
        }

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