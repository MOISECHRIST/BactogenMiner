#!/home/mmeka/.local/bin/nextflow

/*
AUTHOR : MEKA Moise
EMAIL : moise.meka@students.unibe.ch
BAKTA : Genome annotation using Bakta
    INPUT : 
        assembly_scafold : path to the input scafolds fasta file
        sample_name : string with the name of the sample 
    OUTPUT :
        bakta_results : path to the bakta results
*/


process BAKTA{
    tag "${sample_name}"
    label 'high'
    publishDir "${params.outdir}/${sample_name}/Annotation", mode: 'copy'

    input:
    tuple val(sample_name), path(scafolds)
    path(bakta_db)

    output:
    tuple val(sample_name), path("Bakta"), emit: bakta_results

    script:
    """
    bakta --db '${bakta_db}' --output 'Bakta/db-light' --prefix '${sample_name}' \\
     --threads '${task.cpus}' '${scafolds}'
    """
}