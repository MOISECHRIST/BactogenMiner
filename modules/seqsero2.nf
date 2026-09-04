#!/home/mmeka/.local/bin/nextflow

/*
AUTHOR : MEKA Moise
EMAIL : moise.meka@students.unibe.ch
SEQSERO2 : Serotyping of Salmonella Spp using SeqSero2
    INPUT : 
        assembly_scafold : path to the input scafolds fasta file
        sample_name : string with the name of the sample 
        sample_type : string with '4' for genome assembly and '5' for nanopore reads (fasta/fastq)
    OUTPUT :
        seqsero2_results : path to the pasty results
*/


process SEQSERO2{
    tag "${sample_name}"
    label 'high'
    publishDir "${params.outdir}/${sample_name}/SeqSero2", mode: 'copy'

    input:
    tuple val(sample_name), path(scafolds)
    val(sample_type) //4 or 5

    output:
    tuple val(sample_name), path("SeqSero_result_*/SeqSero_result.txt"), emit: seqsero2_results

    script:
    """
    SeqSero2_package.py -m k -t ${sample_type} -i ${scafolds}
    """
}