#!/home/mmeka/.local/bin/nextflow

/*
AUTHOR : MEKA Moise
EMAIL : moise.meka@students.unibe.ch
FASTP : Process to apply trimming, filtering and cleanning on fastq reads using fastp
    INPUT : 
        reads : path to the input fastq file
        sample_name : string with the name of the sample 
    OUTPUT :
        fastp_results : path to a folder containing all fastp results including the html report file, the json data file and cleaned reads
*/

process FASTP{
    tag "${sample_name}"
    label 'high'
    publishDir "${params.outdir}/${sample_name}/fastp", mode: 'copy'

    input: 
    tuple val(sample_name), path(reads)

    output:
    tuple val(sample_name), path("${sample_name}_trimmed_R*.fastq.gz"), emit: cleaned_reads
    tuple val(sample_name), path("${sample_name}_report.fastp.json"), emit: json
    tuple val(sample_name), path("${sample_name}_report.fastp.html"), emit: html
    path("${sample_name}_unpaired*.fastq.gz"), emit: unpaired_reads, optional: true

    script:
    if (reads instanceof List && reads.size() == 2) {
        """
        fastp -i '${reads[0]}' -I '${reads[1]}' \\
        -o '${sample_name}_trimmed_R1.fastq.gz' \\
        -O '${sample_name}_trimmed_R2.fastq.gz' \\
        --unpaired1 '${sample_name}_unpaired_R1.fastq.gz' \\
        --unpaired2 '${sample_name}_unpaired_R2.fastq.gz' \\
        --json '${sample_name}_report.fastp.json' \\
        --html '${sample_name}_report.fastp.html' \\
        --detect_adapter_for_pe \\
        --thread ${task.cpus} --cut_front --cut_tail
        """
    } else {
        """ 
        fastp -i '${reads}' \\
        -o '${sample_name}_trimmed_R1.fastq.gz' \\
        --json '${sample_name}_report.fastp.json' \\
        --html '${sample_name}_report.fastp.html' \\
        --thread ${task.cpus} --cut_front --cut_tail
        """
    }
}