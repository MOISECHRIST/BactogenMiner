#!/home/mmeka/.local/bin/nextflow

/*
AUTHOR : MEKA Moise
EMAIL : moise.meka@students.unibe.ch
BRACKEN : Process to compute specie abundance using Bracken
    INPUT : 
        kraken2_report : path to the kraken2 report
        sample_name : string with the name of the sample  
    OUTPUT :
        bracken_output : path to the bracken output text file
        bracken_report : path to the bracken report text file
*/


process BRACKEN{
    tag "${sample_name}"
    label 'high'

    input:
    tuple val(sample_name), path(kraken2_report)

    output:
    tuple val(sample_name), path("${sample_name}/kraken2_Bracken/bracken_output.txt"), emit: bracken_output
    tuple val(sample_name), path("${sample_name}/kraken2_Bracken/bracken_report.txt"), emit: bracken_report

    script:
    """
    #Make the result directory
    mkdir -p '${sample_name}/kraken2_Bracken'

    #To identify all species abundance in your sample you can use bracken 
    bracken -d '${params.kraken_db}' \\
    -i '${kraken2_report}' \\
    -o '${sample_name}/kraken2_Bracken/bracken_output.txt' \\
    -w '${sample_name}/kraken2_Bracken/bracken_report.txt' \\
    -r '${params.kraken_read_len}'' -l '${params.bracken_class_level}' -t '${params.bracken_threshold}'

    ## NOTE :
    ## -> The flag -r <READ_LEN> is for the read length to get all classifications for (default: 100)
    ## -> The flag -l S is for the level to estimate abundance at [options: D,P,C,O,F,G,S,S1,etc] (default: S)
    ## -> The flag -t <THRESHOLD> is for the number of reads required PRIOR to abundance estimation to perform reestimation (default: 0)
    """
}