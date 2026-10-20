// COMPUTE_TRUNCATION: truncation-hotspot enrichment along the vector genome
// from the vector BAM. The .bai is staged because pysam .fetch needs it.

process COMPUTE_TRUNCATION {
    tag "${meta.id}"

    publishDir { "${params.outdir}/${meta.id}/stats" }, mode: 'copy', pattern: "truncation_hotspots.*"

    cpus   2
    memory 4.GB
    time   1.h

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("truncation_hotspots.csv"), path("truncation_hotspots.json"), emit: truncation

    script:
    """
    python /app/aav_analyze.py compute-truncation \\
        --bam ${bam} \\
        --out-csv truncation_hotspots.csv \\
        --out-json truncation_hotspots.json
    """
}
