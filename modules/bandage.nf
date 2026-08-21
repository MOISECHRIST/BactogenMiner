#!/home/mmeka/.local/bin/nextflow

/*
AUTHOR : MEKA Moise
EMAIL : moise.meka@students.unibe.ch
BANDAGE : Process to do the assembly results quality control using Bandage
    INPUT : 
        assembly_graph : path to the input scafolds fasta file
        sample_name : string with the name of the sample 
    OUTPUT :
        bandage_img : path to a image of the graph
*/


process BANDAGE{
    tag "${sample_name}"
    label 'low'
    publishDir "${params.outdir}/${sample_name}/assembly_QC", mode: 'copy'

    input:
    tuple val(sample_name), path(assembly_graph)

    output:
    tuple val(sample_name), path("${sample_name}_Bandage_Img.jpg"), emit: bandage_img

    script:
    """
    Bandage image '${assembly_graph}' '${sample_name}_Bandage_Img.jpg'
    """
}