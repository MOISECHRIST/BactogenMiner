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


process KRAKEN_BRACKEN{
    tag "${sample_name}"
    label 'high'

    input:
    tuple val(sample_name), path(scafolds)

    output:
    tuple val(sample_name), path("${sample_name}/kraken2_Bracken/kraken2_output.txt"), emit: kraken2_output
    tuple val(sample_name), path("${sample_name}/kraken2_Bracken/kraken2_report.txt"), emit: kraken2_report

    script:
    """
    #Make the result directory
    mkdir -p '${sample_name}/kraken2_Bracken'

    #Now run Kraken2
    kraken2 --db '${params.kraken_db}' --threads '${task.cpus}' \\
    --output '${sample_name}/kraken2_Bracken/kraken2_output.txt' \\
    --report '${sample_name}/kraken2_Bracken/kraken2_report.txt' \\
    '${scafolds}' --memory-mapping

    ## NOTE : 
    ## -> The flag --memory-mapping is for someone who do not need to load all the database in the RAM
    """
}