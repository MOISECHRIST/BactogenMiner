#!/home/mmeka/.local/bin/nextflow

/*
AUTHOR : MEKA Moise
EMAIL : moise.meka@students.unibe.ch
KLEBORATE : Process to screen genome assemblies of Klebsiella pneumoniae species complex,
 the Klebsiella oxytoca species complex and Escherichia 
    INPUT : 
        scafolds : path to the input scafolds fasta file
        sample_name : string with the name of the sample 
        kleborate_preset : string with the presept for kleborate [kpsc|kosc|escherichia]
    OUTPUT :
        kleborate_results : path to the kleborate result directory
*/


process KLEBORATE{
    tag "${sample_name}"
    label 'high'
    publishDir "${params.outdir}/${sample_name}", mode: 'copy'

    input:
    tuple val(sample_name), path(scafolds)
    val(kleborate_preset)

    output:
    tuple val(sample_name), path("Kleborate"), emit: kleborate_results

    script:
    """
    kleborate --assemblies ${scafolds} --outdir Kleborate --preset ${kleborate_preset}  --trim_headers
    """
}