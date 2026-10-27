#!/usr/bin/env nextflow

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
    spades.py ${input_reads} -o spades --threads ${task.cpus}
    status=\$?
    set -e

    if [ \$status -eq 0 ] && [ -s spades/scaffolds.fasta ]; then
        mv spades/scaffolds.fasta 'spades/${sample_name}.fasta'
        mv spades/assembly_graph.fastg 'spades/${sample_name}.fastg'
        ASSEMBLY_END=PASS
    else
        log=spades/spades.log
        reason=""

        if [ -s "\$log" ]; then
            # 1) Wrapper-level error: "== Error ==  <message>"
            reason=\$(grep -m1 -E '== Error ==' "\$log" | sed -E 's/.*== Error ==[[:space:]]*//')

            # 2) Fallback: first C++ stage ERROR line
            [ -z "\$reason" ] && reason=\$(grep -m1 -E '\\bERROR\\b' "\$log" | sed -E 's/^.*ERROR[[:space:]]*//')

            # 3) Fallback: last non-empty line of the log
            [ -z "\$reason" ] && reason=\$(grep -v '^[[:space:]]*\$' "\$log" | tail -n1)
        fi

        # Fallbacks if the log is missing/empty
        if [ -z "\$reason" ]; then
            [ \$status -ne 0 ] && reason="spades_exit_\$status" || reason="empty_or_missing_scaffolds"
        fi

        echo "SPAdes failed for ${sample_name}: \$reason" >&2
        ASSEMBLY_END="FAIL:\$reason"
    fi
    """
}