#!/home/mmeka/.local/bin/nextflow

/*
AUTHOR : MEKA Moise
EMAIL : moise.meka@students.unibe.ch
JOIN_SPECIES_TOOLS : Process to assign for each species the right tools 
    INPUT : 
        species_name : path to the input scafolds fasta file
        sample_name : string with the name of the sample  
    OUTPUT :
        tools : string with the name of the tools
        sp_group : string with the species group 
*/


process JOIN_SPECIES_TOOLS{
    tag "${sample_name}"
    label 'low'
    publishDir "${params.outdir}/${sample_name}", mode: 'copy'

    input:
    tuple val(sample_name), val(species_name)

    output:
    tuple val(sample_name), env("sp_group"), emit: sp_group
    tuple val(sample_name), env("tools"), emit: tools

    script:
     """
    set -eu

    result=\$(awk -F',' -v species="${species_name}" '
        match(species, \$1) > 0 { print \$2","\$3; found=1; exit }
        END { if (!found) exit 1 }
    ' '${projectDir}/meta_data/species_tools.csv') || {
        echo "ERROR: no matching tools/group found for species '${species_name}'" >&2
        exit 1
    }

    sp_group=\$(echo "\${result}" | cut -d',' -f1)
    tools=\$(echo "\${result}" | cut -d',' -f2)
    """
}