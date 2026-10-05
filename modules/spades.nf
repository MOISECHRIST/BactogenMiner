#!/home/mmeka/.local/bin/nextflow

/*
AUTHOR : MEKA Moise
EMAIL : moise.meka@students.unibe.ch
SPADES : Process to do de novo assembly of our genome using spades
    INPUT : 
        reads : path to the input fastq file
        sample_name : string with the name of the sample 
    OUTPUT :
        spades_results : path to a folder containing all spades results
*/

process SPADES {
    tag "${sample_name}"
    label 'high'

    publishDir "${params.outdir}/${sample_name}", pattern: "spades/{${sample_name}.fasta,${sample_name}.fastg}", mode: 'copy'

    input:
    tuple val(sample_name), path(reads)

    output:
    tuple val(sample_name), path("spades/${sample_name}.fasta"), emit: scafolds, optional: true
    tuple val(sample_name), path("spades/${sample_name}.fastg"), emit: assembly_graph, optional: true
    tuple val(sample_name), env('ASSEMBLY_END'), emit: assembly_status

    script:
    def input_reads = (reads instanceof List && reads.size() == 2) ?
        "-1 '${reads[0]}' -2 '${reads[1]}'" : "-s '${reads}'"
    """
    mkdir -p spades

    set +e
    spades.py ${input_reads} -o spades --threads ${task.cpus} --phred-offset 33
    status=\$?
    set -e

    if [ \$status -eq 0 ] && [ -s spades/scaffolds.fasta ]; then
        mv spades/scaffolds.fasta 'spades/${sample_name}.fasta'
        mv spades/assembly_graph.fastg 'spades/${sample_name}.fastg'
        ASSEMBLY_END=PASS
    else
        echo "SPAdes failed for ${sample_name} (exit \$status), see spades/spades.log" >&2
        ASSEMBLY_END=FAIL
    fi
    """
}