// COMPUTE_COVERAGE: per-position depth for one BAM via samtools depth, one
// category per run (vector, helper, rep_cap, spike_in).

process COMPUTE_COVERAGE {
    tag "${meta.id}:${category}"

    publishDir { "${params.outdir}/${meta.id}/stats/coverage" }, mode: 'copy', pattern: "coverage_*.json"

    cpus   2
    memory 4.GB
    time   1.h

    input:
    tuple val(meta), val(category), path(bam), path(bai)

    output:
    tuple val(meta), val(category), path("coverage_${category}.json"), emit: coverage

    script:
    """
    python /app/aav_analyze.py compute-coverage \\
        --bam ${bam} \\
        --category ${category} \\
        --out-json coverage_${category}.json
    """
}
