#!/usr/bin/env nextflow

/*
AUTHOR : MEKA Moise
EMAIL : moise.meka@students.unibe.ch
ECTYPER_GET_DB : Doawnload ECTyper database 
    INPUT:
        trigger : Trigger the database dowloading 
    OUTPUT :
        ectyper_db : path to the ectyper database 
*/


process ECTYPER_GET_DB{
    label 'high'
    publishDir "${params.outdir}/ECTyper_DB", mode: 'copy'

    input:
    val(trigger)

    output:
    path("EnteroRef_GTDBSketch_20231003_V2.msh"), emit: ectyper_db

    script:
    """
    wget https://zenodo.org/records/13969103/files/EnteroRef_GTDBSketch_20231003_V2.msh
    mash info -t EnteroRef_GTDBSketch_20231003_V2.msh > EnteroRef_GTDBSketch_20231003_V2.msh.txt
    """
}