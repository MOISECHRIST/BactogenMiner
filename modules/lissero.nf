#!/home/mmeka/.local/bin/nextflow

/*
AUTHOR : MEKA Moise
EMAIL : moise.meka@students.unibe.ch
LISSERO : Serotyping of Listeria monocytogenes using LisSero
    INPUT : 
        assembly_scafold : path to the input scafolds fasta file
        sample_name : string with the name of the sample 
    OUTPUT :
        lissero_results : path to the pasty results
*/


process LISSERO{
    tag "${sample_name}"
    label 'high'
    publishDir "${params.outdir}/${sample_name}/LisSero", mode: 'copy'

    input:
    tuple val(sample_name), path(scafolds)

    output:
    tuple val(sample_name), path("${sample_name}_lissero_report.txt"), emit: lissero_results

    script:
    """
    lissero ${scafolds} > ${sample_name}_lissero_report.txt
    """
}