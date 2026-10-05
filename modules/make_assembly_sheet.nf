#!/usr/bin/env nextflow

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
    val(samples)   // list of [sample_name, scafolds]

    output:
    path("assembly_sample_sheet.txt"), emit: assembly_sheet

    script:
    def rows = samples.collect { sn, sc -> "${sn}\t${sc.toRealPath()}" }.join("\n")
    """
    cat > assembly_sample_sheet.txt <<'EOF'
${rows}
EOF
    """
}