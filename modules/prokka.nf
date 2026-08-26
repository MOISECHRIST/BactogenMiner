#!/home/mmeka/.local/bin/nextflow

/*
AUTHOR : MEKA Moise
EMAIL : moise.meka@students.unibe.ch
PROKKA : Genome annotation using Prokka
    INPUT : 
        assembly_scafold : path to the input scafolds fasta file
        sample_name : string with the name of the sample 
    OUTPUT :
        prokka_results : path to the prokka results
*/


process PROKKA{
    tag "${sample_name}"
    label 'high'
    publishDir "${params.outdir}/${sample_name}/Annotation", mode: 'copy'

    input:
    tuple val(sample_name), path(scafolds)

    output:
    tuple val(sample_name), path("Prokka"), emit: prokka_results

    script:
    """
    prokka --outdir Prokka --prefix '${sample_name}' \\
     --cpus '${task.cpus}' '${scafolds}'
    """
}