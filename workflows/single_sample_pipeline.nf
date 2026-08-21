include { FASTQC  } from "../modules/fastqc.nf"
include { FASTP   } from "../modules/fastp.nf"
include { SPADES  } from "../modules/spades.nf"
include { QUAST   } from "../modules/quast.nf"
include { BANDAGE } from "../modules/bandage.nf"
include { KRAKEN2 } from "../modules/kraken2.nf"
include { BRACKEN } from "../modules/Bracken.nf"
include { MLST    } from "../modules/mlst.nf"
include { BUSCO   } from "../modules/busco.nf"
include { GAMBIT    } from "../modules/gambit.nf"


workflow SINGLE_SAMPLE_PROCESSING {
    take:
        reads_ch  // channel: [ val(sample_name), [ path(reads) ] ]

    main:
        // Raw reads quality control
        FASTQC(reads_ch)
        FASTP(reads_ch)

        // De novo assembly 
        SPADES(FASTP.out.cleaned_reads)

        // Assembly quality control
        BANDAGE(SPADES.out.assembly_graph)
        QUAST(SPADES.out.scafolds)
        BUSCO(SPADES.out.scafolds)

        // Species classification 
        if (params.use_gambit) {
            GAMBIT(SPADES.out.scafolds, file(params.gambit_db))
        }

        if (params.use_kraken2) {
            KRAKEN2(SPADES.out.scafolds, file(params.kraken_db))
            BRACKEN(KRAKEN2.out.kraken_out, file(params.kraken_db))
        }

        // Sequence typing
        MLST(SPADES.out.scafolds)
}
