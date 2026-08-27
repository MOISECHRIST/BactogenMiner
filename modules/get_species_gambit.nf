#!/home/mmeka/.local/bin/nextflow

/*
AUTHOR : MEKA Moise
EMAIL : moise.meka@students.unibe.ch
GET_SPECIES_GAMBIT : Get from gambit results the species name
    INPUT : 
        gambit_report : path to the gambit report
        sample_name : string with the name of the sample 
    OUTPUT :
        species_name : the species name
*/


process GET_SPECIES_GAMBIT{
    tag "${sample_name}"
    label 'low'

    input:
    tuple val(sample_name), path(gambit_report)

    output:
    tuple val(sample_name), env('TOP_SPECIES'), emit: species_name

    script:
    """
    TOP_SPECIES=\$(tail -n1 '${gambit_report}' | cut -d, -f2)
    """
}