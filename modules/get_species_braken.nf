#!/home/mmeka/.local/bin/nextflow

/*
AUTHOR : MEKA Moise
EMAIL : moise.meka@students.unibe.ch
GET_SPECIES_BRAKEN : Get from braken output the species name
    INPUT : 
        braken_output : path to the braken output
        sample_name : string with the name of the sample 
    OUTPUT :
        species_name : the species name
*/


process GET_SPECIES_BRAKEN{
    tag "${sample_name}"
    label 'low'

    input:
    tuple val(sample_name), path(braken_output)

    output:
    tuple val(sample_name), env('TOP_SPECIES'), emit: species_name

    script:
    """
    TOP_SPECIES=\$(tail -n+2 '${braken_output}' | cut -f1,7   | sort -t \$'\\t' -k2,2 -nr | head -1 | cut -f1)
    """
}