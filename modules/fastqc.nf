#!/home/mmeka/.local/bin/nextflow

/*
AUTHOR : MEKA Moise
EMAIL : moise.meka@students.unibe.ch
FASTQC : Process to assess the quality of reads in fastq format using fastqc tool
    INPUT : 
        reads : path to the input fastq file
        sample_name : string with the name of the sample 
    OUTPUT :
        fastqc_results : path to a folder containing a html report file and a .zip archive file containing data
*/

process FASTQC{
    tag "${sample_name}"
    label 'high'
    publishDir "${params.outdir}/${sample_name}/fastqc", mode: 'copy'


    input:
    tuple val(sample_name), path(reads)

    output:
    tuple val(sample_name), path("*.html"), emit: fastqc_html
    tuple val(sample_name), path("*.zip"), emit: fastqc_zip

    script:
    """
    fastqc -t '${task.cpus}' ${reads} 
    """
}