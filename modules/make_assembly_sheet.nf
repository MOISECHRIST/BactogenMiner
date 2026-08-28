#!/home/mmeka/.local/bin/nextflow

/*
AUTHOR : MEKA Moise
EMAIL : moise.meka@students.unibe.ch
MAKE_ASSEMBLY_SHEET : Compile all scafolds assembly fasta into a sample sheet
    INPUT : 
        assembly_scafold : path to the input scafolds fasta file
        sample_name : string with the name of the sample 
    OUTPUT :
        abricate_results : path to the abricate results
*/


process MAKE_ASSEMBLY_SHEET {
    label 'low'
    publishDir "${params.outdir}/phylogeny", mode: 'copy'

    input:
    tuple val(sample_name), path(scafolds)

    output:
    path("assembly_sample_sheet.txt"), emit: assembly_sheet

    script:
    """
    echo "${sample_name}\t${scafolds.toRealPath()}" >> assembly_sample_sheet.txt
    """
}