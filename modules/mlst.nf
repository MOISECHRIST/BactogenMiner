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

    input:
    tuple val(sample_name), path(scafolds)

    output:
    tuple val(sample_name), path("${sample_name}/mlst/report.csv"), emit: mlst_output

    script:
    if(params.mlst_scheme){
        """
        #Make the result directory
        mkdir -p '${sample_name}/mlst'

        #Then run MLST for sequence typing 
        mlst --legacy --scheme '${params.mlst_scheme}' '${scafolds}' --csv > '${sample_name}/mlst/report.csv'
        """
    } else {
        """
        #Make the result directory
        mkdir -p '${sample_name}/mlst'

        #Then run MLST for sequence typing 
        mlst --full '${scafolds}' --csv > '${sample_name}/mlst/report.csv'
        """
    }
    
}