#!/home/mmeka/.local/bin/nextflow

/*
AUTHOR : MEKA Moise
EMAIL : moise.meka@students.unibe.ch
BUSCO : Process to do the assembly results quality control using Busco
    INPUT : 
        assembly_graph : path to the input scafolds fasta file
        sample_name : string with the name of the sample 
    OUTPUT :
        busco_json_summary : path to the busco summary in json format
        busco_txt_summary : path to the busco summary in txt format
*/


process BUSCO{
    tag "${sample_name}"
    label 'high'
    publishDir "${params.outdir}/${sample_name}/assembly_QC", mode: 'copy'

    input:
    tuple val(sample_name), path(scafolds)

    output:
    tuple val(sample_name), path("Busco/short_summary.*.txt"), emit: busco_txt_summary
    tuple val(sample_name), path("Busco/short_summary.*.json"), emit: busco_json_summary

    script:
    """
    busco -i '${scafolds}' -m genome -l '${params.busco_lineage}' --cpu '${task.cpus}' \\
    --out  Busco
    """
}