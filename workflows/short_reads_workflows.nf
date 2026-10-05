#!/usr/bin/env nextflow

include { CHECK_READS               }      from "../modules/check_reads.nf"
include { CHECK_READS as CHECK_TRIM }      from "../modules/check_reads.nf"
include { FASTQC                    }      from "../modules/fastqc.nf"
include { FASTP                     }      from "../modules/fastp.nf"
include { SPADES                    }      from "../modules/spades.nf"
include { QUAST                     }      from "../modules/quast.nf"
include { BANDAGE                   }      from "../modules/bandage.nf"
include { BUSCO                     }      from "../modules/busco.nf"
include { PROKKA                    }      from "../modules/prokka.nf"
include { BAKTA                     }      from "../modules/bakta.nf"
include { BAKTA_BD                  }      from "../modules/bakta_db.nf"

workflow SAMPLE_QC {
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
        // Quality control
        FASTQC(passed_reads_ch)

        // Cleaning and Trimming 
        FASTP(passed_reads_ch)
        CHECK_TRIM(FASTP.out.cleaned_reads)
        passed_trim_ch = CHECK_TRIM.out.validated_reads
            .join(CHECK_TRIM.out.read_screen)
            .filter { sample_name, reads, flag ->
                if (flag != "PASS") {
                    log.warn "Sample ${sample_name} failed read screening: ${flag}"
                }
                flag == "PASS"
            }
            .map { sample_name, reads, flag -> tuple(sample_name, reads) }
    
    emit:
        fastqc_zip = FASTQC.out.fastqc_zip
        fastp_json = FASTP.out.json
        cleaned_reads = passed_trim_ch
}

workflow ASSEMBLY_ANNOTATION {
    take:
        cleaned_reads // channel: [ val(sample_name), [ path(reads) ] ]
    
    main:
        // De novo assembly
        SPADES(cleaned_reads)

        SPADES.out.assembly_status
            .filter { sample_name, status -> status != 'PASS' }
            .subscribe { sample_name, status ->
                log.warn "Assembly failed for ${sample_name}: skipped downstream (check spades/spades.log in its work dir)"
            }
        
        passed_status_ch = SPADES.out.assembly_status
            .filter { sample_name, status -> status == 'PASS' }

        passed_scaffolds_ch = SPADES.out.scafolds
            .join(passed_status_ch)
            .map { sample_name, fasta, status -> tuple(sample_name, fasta) }

        passed_graph_ch = SPADES.out.assembly_graph
            .join(passed_status_ch)
            .map { sample_name, gfa, status -> tuple(sample_name, gfa) }

        // Assembly QC
        BANDAGE(passed_graph_ch)
        QUAST(passed_scaffolds_ch)
        BUSCO(passed_scaffolds_ch)

        // Genome Annotation
        if (params.use_bakta) {
            if (!params.bakta_db){
                BAKTA_BD()
                bakta_db_ch = BAKTA_BD.out.bakta_bd
            } else {
                bakta_db_ch =file(params.bakta_db)
            }
            BAKTA(passed_scaffolds_ch, bakta_db_ch)
            genome_annotation = BAKTA.out.bakta_results
        } else {
            PROKKA(passed_scaffolds_ch)
            genome_annotation = PROKKA.out.prokka_results
        }
    emit:
        scafolds = passed_scaffolds_ch
        assembly_status    = SPADES.out.assembly_status
        quast_report = QUAST.out.quast_report
        busco_json_summary = BUSCO.out.busco_json_summary
        busco_txt_summary = BUSCO.out.busco_txt_summary
        genome_annotation = genome_annotation
}