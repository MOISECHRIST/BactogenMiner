#!/home/mmeka/.local/bin/nextflow

/*
AUTHOR : MEKA Moise
EMAIL : moise.meka@students.unibe.ch
SPADES : Process to do de novo assembly of our genome using spades
    INPUT : 
        reads : path to the input fastq file
        sample_name : string with the name of the sample 
    OUTPUT :
        spades_results : path to a folder containing all spades results
*/


process SPADES{
    tag "${sample_name}"
    label 'high'

    input:
    tuple val(sample_name), path(reads)

    output:
    tuple val(sample_name), path("${sample_name}/spades/scaffolds.fasta"), emit: spades_scafolds
    tuple val(sample_name), path("${sample_name}/spades/assembly_graph.fastg"), emit: spades_assembly_graph

    script:
    """
    #Make the result directory
    mkdir -p '${sample_name}/spades'

    #Now run the assembly with spades
    spades.py -1 '${reads[0]}' -2 '${reads[1]}' -o '${sample_name}/spades' \\
    --threads '${task.cpus}'
    """
}