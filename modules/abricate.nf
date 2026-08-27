#!/home/mmeka/.local/bin/nextflow

/*
AUTHOR : MEKA Moise
EMAIL : moise.meka@students.unibe.ch
ABRICATE : Finding of virulence factors, resistance genes and plasmid using abricate 
    INPUT : 
        assembly_scafold : path to the input scafolds fasta file
        sample_name : string with the name of the sample 
        database_name : string with the name of the database name
    OUTPUT :
        abricate_results : path to the abricate results
*/


process ABRICATE{
    tag "${sample_name}"
    label 'high'
    publishDir "${params.outdir}/${sample_name}/Abricate", mode: 'copy'

    input:
    tuple val(sample_name), path(scafolds)
    val(database_name)

    output:
    tuple val(sample_name), path("${sample_name}_${database_name}_report.txt"), emit: abricate_results

    script:
    """
    abricate --db '${database_name}' '${scafolds}' > '${sample_name}_${database_name}_report.txt'
    """
}