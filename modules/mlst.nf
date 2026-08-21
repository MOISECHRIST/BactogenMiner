#!/home/mmeka/.local/bin/nextflow

/*
AUTHOR : MEKA Moise
EMAIL : moise.meka@students.unibe.ch
MLST : Process to compute sequence typing using MLST
    INPUT : 
        scafolds : path to the input scafolds fasta file
        sample_name : string with the name of the sample  
    OUTPUT :
        mlst_output : path to the mslt output text file
*/


process MLST{
    tag "${sample_name}"
    label 'high'
    publishDir "${params.outdir}/${sample_name}/mlst", mode: 'copy'

    input:
    tuple val(sample_name), path(scafolds)

    output:
    tuple val(sample_name), path("${sample_name}_mlst_report.csv"), emit: mlst_output

    script:
    if(params.mlst_scheme){
        """
        mlst --legacy --scheme '${params.mlst_scheme}' '${scafolds}' --csv > '${sample_name}_mlst_report.csv'
        """
    } else {
        """
        mlst '${scafolds}' --csv > '${sample_name}_mlst_report.csv'
        """
    } 
}