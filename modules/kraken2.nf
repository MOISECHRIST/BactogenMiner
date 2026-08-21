#!/home/mmeka/.local/bin/nextflow

/*
AUTHOR : MEKA Moise
EMAIL : moise.meka@students.unibe.ch
KRAKEN2 : Process for specie classification using Kraken2
    INPUT : 
        scafolds : path to the input scafolds fasta file
        sample_name : string with the name of the sample 
    OUTPUT :
        kraken2_report : path to the kraken2 report
        kraken2_output : path to the kraken2 report
*/


process KRAKEN2{
    tag "${sample_name}"
    label 'high'
    publishDir "${params.outdir}/${sample_name}/kraken2_Bracken", mode: 'copy'

    input:
    tuple val(sample_name), path(scafolds)
    path(kraken_db)

    output:
    tuple val(sample_name), path("${sample_name}_kraken2_report.txt"), path("${sample_name}_kraken2_output.txt"), emit: kraken_out
    tuple val(sample_name), path("${sample_name}_kraken2_report.txt"), emit: kraken2_report
    tuple val(sample_name), path("${sample_name}_kraken2_output.txt"), emit: kraken2_output

    script:
    """
    kraken2 --db '${kraken_db}' --threads '${task.cpus}' \\
    --output '${sample_name}_kraken2_output.txt' \\
    --report '${sample_name}_kraken2_report.txt' \\
    '${scafolds}' --memory-mapping
    """
}