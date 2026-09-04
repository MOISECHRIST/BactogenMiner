#!/home/mmeka/.local/bin/nextflow

/*
AUTHOR : MEKA Moise
EMAIL : moise.meka@students.unibe.ch
PASTY : Serotyping of Pseudomonas aeruginosa using Pasty
    INPUT : 
        assembly_scafold : path to the input scafolds fasta file
        sample_name : string with the name of the sample 
    OUTPUT :
        pasty_results : path to the pasty results
*/


process PASTY{
    tag "${sample_name}"
    label 'high'
    publishDir "${params.outdir}/${sample_name}", mode: 'copy'

    input:
    tuple val(sample_name), path(scafolds)

    output:
    tuple val(sample_name), path("Pasty"), emit: pasty_results

    script:
    """
    pasty --input ${scafolds} --outdir Pasty \\
     --pathotype --prefix ${sample_name}
    """
}