#!/home/mmeka/.local/bin/nextflow

/*
AUTHOR : MEKA Moise
EMAIL : moise.meka@students.unibe.ch
QUAST : Process to do the assembly results quality control using quast
    INPUT : 
        scafolds : path to the input scafolds fasta file
        sample_name : string with the name of the sample 
    OUTPUT :
        quast_results : path to a folder containing all quast results
*/


process QUAST{
    tag "${sample_name}"
    label 'high'
    publishDir "${params.outdir}/${sample_name}/assembly_QC", mode: 'copy'

    input:
    tuple val(sample_name), path(scafolds)

    output:
    tuple val(sample_name), path("quast"), emit: quast_report

    script:
    """
    quast.py '${scafolds}' -o 'quast' --threads '${task.cpus}' 
    """
}