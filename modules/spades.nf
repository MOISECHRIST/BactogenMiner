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
    publishDir "${params.outdir}/${sample_name}", pattern: "spades/{scaffolds.fasta,assembly_graph.fastg}", mode: 'copy'

    input:
    tuple val(sample_name), path(reads)

    output:
    tuple val(sample_name), path("spades/${sample_name}.fasta"), emit: scafolds
    tuple val(sample_name), path("spades/${sample_name}.fastg"), emit: assembly_graph

    script:
    if (reads instanceof List && reads.size() == 2) {
        """
        mkdir spades

        spades.py -1 '${reads[0]}' -2 '${reads[1]}' -o 'spades' \\
        --threads '${task.cpus}'

        mv spades/scaffolds.fasta 'spades/${sample_name}.fasta'
        mv spades/assembly_graph.fastg 'spades/${sample_name}.fastg'
        """
    } else {
        """
        mkdir spades
        
        spades.py -s '${reads}' -o 'spades' \\
        --threads '${task.cpus}'

        mv spades/scaffolds.fasta 'spades/${sample_name}.fasta'
        mv spades/assembly_graph.fastg 'spades/${sample_name}.fastg'
        """
    }
}
