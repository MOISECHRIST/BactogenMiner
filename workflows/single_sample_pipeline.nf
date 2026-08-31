include { CHECK_READS                  }      from "../modules/check_reads.nf"
include { FASTQC                       }      from "../modules/fastqc.nf"
include { FASTP                        }      from "../modules/fastp.nf"
include { SPADES                       }      from "../modules/spades.nf"
include { QUAST                        }      from "../modules/quast.nf"
include { BANDAGE                      }      from "../modules/bandage.nf"
include { KRAKEN2                      }      from "../modules/kraken2.nf"
include { BRACKEN                      }      from "../modules/Bracken.nf"
include { MLST                         }      from "../modules/mlst.nf"
include { BUSCO                        }      from "../modules/busco.nf"
include { GAMBIT                       }      from "../modules/gambit.nf"
include { PROKKA                       }      from "../modules/prokka.nf"
include { BAKTA                        }      from "../modules/bakta.nf"
include { BAKTA_BD                     }      from "../modules/bakta_db.nf"
include { ABRICATE as ABRICATE_VFDB    }      from "../modules/abricate.nf"
include { ABRICATE as ABRICATE_AMR     }      from "../modules/abricate.nf"
include { ABRICATE as ABRICATE_PLASMID }      from "../modules/abricate.nf"
include { ABRICATE as ABRICATE_ECOLI   }      from "../modules/abricate.nf"
include { GET_SPECIES_BRAKEN           }      from "../modules/get_species_braken.nf"
include { GET_SPECIES_GAMBIT           }      from "../modules/get_species_gambit.nf"
include { JOIN_SPECIES_TOOLS           }      from "../modules/join_species_with_tools.nf"
include { KLEBORATE                    }      from "../modules/kleborate.nf"


workflow SINGLE_SAMPLE_PROCESSING {
    take:
        reads_ch  // channel: [ val(sample_name), [ path(reads) ] ]

    main:
        // Raw reads quality control
        CHECK_READS(reads_ch)
        passed_reads_ch = CHECK_READS.out.validated_reads
            .join(CHECK_READS.out.read_screen)
            .filter { sample_name, reads, flag ->
                if (flag != "PASS") {
                    log.warn "Sample ${sample_name} failed read screening: ${flag}"
                }
                flag == "PASS"
            }
            .map { sample_name, reads, flag -> tuple(sample_name, reads) }
        FASTQC(passed_reads_ch)
        FASTP(passed_reads_ch)

        // De novo assembly 
        SPADES(FASTP.out.cleaned_reads)

        // Assembly quality control
        BANDAGE(SPADES.out.assembly_graph)
        QUAST(SPADES.out.scafolds)
        BUSCO(SPADES.out.scafolds)

        // Genome Annotation
        if (params.use_bakta) {
            if (!params.bakta_db){
                BAKTA_BD()
                bakta_db_ch = BAKTA_BD.out.bakta_bd
            } else {
                bakta_db_ch =file(params.bakta_db)
            }
            BAKTA(SPADES.out.scafolds, bakta_db_ch)
        } else {
            PROKKA(SPADES.out.scafolds)
        }

        // Species classification 
        if (params.use_gambit) {
            GAMBIT(SPADES.out.scafolds, file(params.gambit_db))
            GET_SPECIES_GAMBIT(GAMBIT.out.gambit_report)
            species_name = GET_SPECIES_GAMBIT.out.species_name
            gambit_report_ch  = GAMBIT.out.gambit_report
            bracken_report_ch = channel.empty()
            bracken_output_ch = channel.empty()
            kraken2_output_ch = channel.empty()
            kraken2_report_ch = channel.empty()
        } else {
            KRAKEN2(SPADES.out.scafolds, file(params.kraken_db))
            BRACKEN(KRAKEN2.out.kraken_out, file(params.kraken_db))
            GET_SPECIES_BRAKEN(BRACKEN.out.bracken_output)
            species_name = GET_SPECIES_BRAKEN.out.species_name
            gambit_report_ch  = channel.empty()
            bracken_report_ch = BRACKEN.out.bracken_report
            bracken_output_ch = BRACKEN.out.bracken_output
            kraken2_output_ch = KRAKEN2.out.kraken2_output
            kraken2_report_ch = KRAKEN2.out.kraken2_report
        }

        // Sequence typing
        MLST(SPADES.out.scafolds)

        //Species specific screening genome assemblies
        JOIN_SPECIES_TOOLS(species_name)
        sp_group_ch = JOIN_SPECIES_TOOLS.out.sp_group.map {sp -> sp[1]}
        tools_ch    = JOIN_SPECIES_TOOLS.out.tools.map {tools -> tools[1]}

        routed = SPADES.out.scafolds
            .combine(sp_group_ch)
            .combine(tools_ch)
            .branch { sample_name, scafolds, sp_group, tools ->
                kleborate: tools == "kleborate"
                seqsero:   tools == "SeqSero2"
                lissero:   tools == "LisSero"
                other:     true
            }

        KLEBORATE(
            routed.kleborate.map { sn, sc, sp, tl -> tuple(sn, sc) },
            routed.kleborate.map { sn, sc, sp, tl -> sp }
        )
        ecoli_only = routed.kleborate.filter { sn, sc, sp, tl -> sp == "escherichia" }
        ABRICATE_ECOLI(ecoli_only.map { sn, sc, sp, tl -> tuple(sn, sc) }, "ecoli_vf")

        // SEQSERO2
        // LISSERO

        ABRICATE_AMR(routed.other.map { sn, sc, sp, tl -> tuple(sn, sc) }, "resfinder")
        ABRICATE_VFDB(routed.other.map { sn, sc, sp, tl -> tuple(sn, sc) }, "vfdb")
        ABRICATE_PLASMID(SPADES.out.scafolds, "plasmidfinder")
    
    emit:
        fastqc_zip = FASTQC.out.fastqc_zip
        fastp_json = FASTP.out.json
        scafolds = SPADES.out.scafolds
        quast_report = QUAST.out.quast_report
        busco_json_summary = BUSCO.out.busco_json_summary
        busco_txt_summary = BUSCO.out.busco_txt_summary
        gambit_report       = gambit_report_ch
        bracken_report      = bracken_report_ch
        bracken_output      = bracken_output_ch
        kraken2_output      = kraken2_output_ch
        kraken2_report      = kraken2_report_ch
}
