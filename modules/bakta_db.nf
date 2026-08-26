#!/home/mmeka/.local/bin/nextflow

/*
AUTHOR : MEKA Moise
EMAIL : moise.meka@students.unibe.ch
BAKTA_BD : Genome annotation using Bakta
    OUTPUT :
        bakta_db : path to the bakta database 
*/


process BAKTA_BD{
    label 'low'
    publishDir "${params.outdir}", mode: 'copy'

    output:
    path("Bakta_DB/${params.bakta_db_type == 'light' ? 'db-light' : 'db'}"), emit: bakta_bd

    script:
    """
    bakta_db download --output Bakta_DB --type '${params.bakta_db_type}'
    """
}