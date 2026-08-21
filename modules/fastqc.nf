#!/home/mmeka/.local/bin/nextflow

/*
AUTHOR : MEKA Moise
EMAIL : moise.meka@students.unibe.ch
FASTQC : Process to assess the quality of reads in fastq format using fastqc tool
    INPUT : 
        reads : path to the input fastq file
        sample_name : string with the name of the sample 
    OUTPUT :
        fastqc_results : path to a folder containing a html report file and a .zip archive file containing data
*/

process FASTQC{
    tag "${sample_name}"
    label 'low'


    input:
    tuple val(sample_name), path(reads)

    output:
    path "${sample_name}/fastqc", emit: fastqc_results

    script:
    """
    #Make the result directory 
    mkdir -p '${sample_name}/fastqc'

    #Now run fastqc on our fastq file
    fastqc -t '${task.cpus}' '${reads}' -o '${sample_name}/fastqc'
    """
}