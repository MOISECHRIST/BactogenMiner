#!/home/mmeka/.local/bin/nextflow

/*
AUTHOR : MEKA Moise
EMAIL : moise.meka@students.unibe.ch
CHECK_READS : Screen reads (read count, bp count, R1/R2 balance if paired,
              estimated genome length & coverage via mash) before assembly.
              Handles both paired-end (2 files) and single-end (1 file).
              Adapted from Theiagen TheiaProk/PHB's check_reads.wdl / check_reads_se.wdl
    INPUT :
        sample_name : string with the name of the sample
        reads : path to the input fastq.gz file(s) - either [R1, R2] or [R1] only
    OUTPUT :
        reads : the same reads, only if PASS
        read_screen : "PASS" or "FAIL; <reasons>"
        read_screen_tsv : per-sample metrics table
        est_genome_length : estimated genome length (bp)
*/

process CHECK_READS {
    tag "${sample_name}"
    label 'high'
    publishDir "${params.outdir}/${sample_name}/ReadScreen", mode: 'copy'

    input:
    tuple val(sample_name), path(reads)

    output:
    tuple val(sample_name), path(reads),          emit: validated_reads
    tuple val(sample_name), env('FLAG'),          emit: read_screen
    path("${sample_name}_read_screen.tsv"),       emit: read_screen_tsv
    env('EST_GENOME_LENGTH'),                     emit: est_genome_length

    script:
    def reads_list = reads instanceof List ? reads : [reads]
    def is_paired  = reads_list.size() == 2
    def read1 = reads_list[0]
    def read2 = is_paired ? reads_list[1] : null

    """
    set -eu

    fail_log=""
    estimated_genome_length=0
    estimated_coverage=0

    read1_num=\$(zcat '${read1}' | fastq-scan | grep 'read_total' | sed 's/[^0-9]*\\([0-9]\\+\\).*/\\1/')
    read1_bp=\$(zcat '${read1}' | fastq-scan | grep 'total_bp' | sed 's/[^0-9]*\\([0-9]\\+\\).*/\\1/')

    if [ "${is_paired}" == "true" ]; then
        ## ---------- PAIRED-END ----------
        metrics="read1_count\tread2_count\tread_bp\test_genome_length\tcomment"

        read2_num=\$(zcat '${read2}' | fastq-scan | grep 'read_total' | sed 's/[^0-9]*\\([0-9]\\+\\).*/\\1/')
        read2_bp=\$(zcat '${read2}' | fastq-scan | grep 'total_bp' | sed 's/[^0-9]*\\([0-9]\\+\\).*/\\1/')

        reads_total=\$(( read1_num + read2_num ))
        bp_total=\$(( read1_bp + read2_bp ))

        if [ "\${reads_total}" -le "${params.min_reads}" ]; then
            fail_log+="; total reads (\${reads_total}) below minimum of ${params.min_reads}"
        fi

        percent_read1=\$(python3 -c "print(round((\$read1_bp / (\$read1_bp + \$read2_bp))*100))")
        percent_read2=\$(python3 -c "print(round((\$read2_bp / (\$read1_bp + \$read2_bp))*100))")

        if [ "\${percent_read1}" -lt "${params.min_proportion}" ]; then
            fail_log+="; R1 underrepresented (\${percent_read1}% vs \${percent_read2}%)"
        fi
        if [ "\${percent_read2}" -lt "${params.min_proportion}" ]; then
            fail_log+="; R2 underrepresented (\${percent_read2}% vs \${percent_read1}%)"
        fi

        if [ "\${bp_total}" -le "${params.min_basepairs}" ]; then
            fail_log+="; total basepairs (\${bp_total}) below minimum of ${params.min_basepairs}"
        fi

        mash sketch -o test -k 31 -m 3 -r '${read1}' '${read2}' > mash-output.txt 2>&1 || true
        if [ ! -f test.msh ]; then
            fail_log+="; mash failed - cannot estimate genome length/coverage"
        else
            estimated_genome_length=\$(grep "Estimated genome size:" mash-output.txt | awk '{printf("%.0f", \$4)}')
            estimated_coverage=\$(grep "Estimated coverage:" mash-output.txt | awk '{printf("%d", \$3)}')
            rm -f test.msh mash-output.txt

            if [ "\${estimated_genome_length}" -gt "${params.max_genome_length}" ] || [ "\${estimated_genome_length}" -lt "${params.min_genome_length}" ]; then
                M="-m 10"
                [ "\${estimated_genome_length}" -lt "${params.min_genome_length}" ] && M="-m 1"
                mash sketch -o test -k 31 \${M} -r '${read1}' '${read2}' > mash-output.txt 2>&1
                estimated_genome_length=\$(grep "Estimated genome size:" mash-output.txt | awk '{printf("%.0f", \$4)}')
                estimated_coverage=\$(grep "Estimated coverage:" mash-output.txt | awk '{printf("%d", \$3)}')
                rm -f test.msh mash-output.txt
            fi
        fi

        metrics+="\\n\${read1_num}\\t\${read2_num}\\t\${bp_total}\\t\${estimated_genome_length}"

    else
        ## ---------- SINGLE-END ----------
        metrics="read1_count\tread_bp\test_genome_length\tcomment"

        if [ "\${read1_num}" -le "${params.min_reads}" ]; then
            fail_log+="; number of reads (\${read1_num}) below minimum of ${params.min_reads}"
        fi
        if [ "\${read1_bp}" -le "${params.min_basepairs}" ]; then
            fail_log+="; number of basepairs (\${read1_bp}) below minimum of ${params.min_basepairs}"
        fi

        mash sketch -o test -k 31 -m 3 -r '${read1}' > mash-output.txt 2>&1 || true
        if [ ! -f test.msh ]; then
            fail_log+="; mash failed - cannot estimate genome length/coverage"
        else
            estimated_genome_length=\$(grep "Estimated genome size:" mash-output.txt | awk '{printf("%.0f", \$4)}')
            estimated_coverage=\$(grep "Estimated coverage:" mash-output.txt | awk '{printf("%d", \$3)}')
            rm -f test.msh mash-output.txt

            if [ "\${estimated_genome_length}" -gt "${params.max_genome_length}" ] || [ "\${estimated_genome_length}" -lt "${params.min_genome_length}" ]; then
                M="-m 10"
                [ "\${estimated_genome_length}" -lt "${params.min_genome_length}" ] && M="-m 1"
                mash sketch -o test -k 31 \${M} -r '${read1}' > mash-output.txt 2>&1
                estimated_genome_length=\$(grep "Estimated genome size:" mash-output.txt | awk '{printf("%.0f", \$4)}')
                estimated_coverage=\$(grep "Estimated coverage:" mash-output.txt | awk '{printf("%d", \$3)}')
                rm -f test.msh mash-output.txt
            fi
        fi

        metrics+="\\n\${read1_num}\\t\${read1_bp}\\t\${estimated_genome_length}"
    fi

    ## ---------- Common genome length / coverage checks ----------
    if [ "\${estimated_genome_length}" -ge "${params.max_genome_length}" ]; then
        fail_log+="; estimated genome length (\${estimated_genome_length}) above max of ${params.max_genome_length} bp"
    elif [ "\${estimated_genome_length}" -le "${params.min_genome_length}" ]; then
        fail_log+="; estimated genome length (\${estimated_genome_length}) below min of ${params.min_genome_length} bp"
    fi
    if [ "\${estimated_coverage}" -lt "${params.min_coverage}" ]; then
        fail_log+="; estimated coverage (\${estimated_coverage}x) below min of ${params.min_coverage}x"
    fi

    if [ -z "\${fail_log}" ]; then
        FLAG="PASS"
    else
        FLAG="FAIL\${fail_log}"
    fi

    echo -e "\${metrics}\\t\${FLAG}" > '${sample_name}_read_screen.tsv'
    echo "\${FLAG}"
    EST_GENOME_LENGTH=\${estimated_genome_length}
    """
}