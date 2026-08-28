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
include { ABRICATE as ABRICATE_PLASMID }      from "../modules/abricate.nf"
include { GET_SPECIES_BRAKEN           }      from "../modules/get_species_braken.nf"
include { GET_SPECIES_GAMBIT           }      from "../modules/get_species_gambit.nf"


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
        }

        if (params.use_kraken2) {
            KRAKEN2(SPADES.out.scafolds, file(params.kraken_db))
            BRACKEN(KRAKEN2.out.kraken_out, file(params.kraken_db))
            GET_SPECIES_BRAKEN(BRACKEN.out.bracken_output)

        }

        // Sequence typing
        MLST(SPADES.out.scafolds)

        // Plasmid & Virulence Genes detection 
        ABRICATE_VFDB(SPADES.out.scafolds, "vfdb")
        ABRICATE_PLASMID(SPADES.out.scafolds, "plasmidfinder")
    
    emit:
        scafolds = SPADES.out.scafolds
}
