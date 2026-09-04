#!/home/mmeka/.local/bin/nextflow

/*
AUTHOR : MEKA Moise
EMAIL : moise.meka@students.unibe.ch
ECTYPER : Serotyping and prediction of pathotype of Escherichia coli and Shigella using ECTyper
    INPUT : 
        assembly_scafold : path to the input scafolds fasta file
        sample_name : string with the name of the sample 
    OUTPUT :
        ectyper_results : path to the ectyper results
*/


process ECTYPER{
    tag "${sample_name}"
    label 'high'
    publishDir "${params.outdir}/${sample_name}", mode: 'copy'

    input:
    tuple val(sample_name), path(scafolds)

    output:
    tuple val(sample_name), path("ECTyper"), emit: ectyper_results

    script:
    """
    ectyper --input ${scafolds} --output ECTyper \\
     --pathotype --cores ${task.cpus} 
    """
}