#!/home/mmeka/.local/bin/nextflow

/*
AUTHOR : MEKA Moise
EMAIL : moise.meka@students.unibe.ch
GAMBIT : Process to assign taxonomy using Gambit
    INPUT : 
        assembly_graph : path to the input scafolds fasta file
        sample_name : string with the name of the sample 
    OUTPUT :
        gambit_report : path to the mash report
*/


process GAMBIT{
    tag "${sample_name}"
    label 'high'
    publishDir "${params.outdir}/${sample_name}/Gambit", mode: 'copy'

    input:
    tuple val(sample_name), path(scafolds)
    path gambit_db

    output:
    tuple val(sample_name), path("${sample_name}_gambit_report.csv"), emit: gambit_report

    script:
    """
    gambit --db '${gambit_db}' query -o '${sample_name}_gambit_report.csv' '${scafolds}'
    """
}