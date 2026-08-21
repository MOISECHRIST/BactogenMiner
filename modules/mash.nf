#!/home/mmeka/.local/bin/nextflow

/*
AUTHOR : MEKA Moise
EMAIL : moise.meka@students.unibe.ch
MASH : Process to assign taxonomy using mash
    INPUT : 
        assembly_graph : path to the input scafolds fasta file
        sample_name : string with the name of the sample 
    OUTPUT :
        mash_report : path to the mash report
*/


process MASH{
    tag "${sample_name}"
    label 'high'
    publishDir "${params.outdir}/${sample_name}/mash", mode: 'copy'

    input:
    tuple val(sample_name), path(scafolds)
    path mash_db

    output:
    tuple val(sample_name), path("${sample_name}_mash_report.txt"), emit: mash_report

    script:
    """
    mash sketch -m 2 -o '${sample_name}_sketch' '${scafolds}'
    mash dist '${mash_db}' '${sample_name}_sketch.msh' > '${sample_name}_mash_report.txt'
    """
}