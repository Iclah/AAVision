// ALIGN_VECTOR: align prepared reads to the vector reference (minimap2 -> SAM).

process ALIGN_VECTOR {
    tag "${meta.id}"

    publishDir { "${params.outdir}/${meta.id}/bam" }, mode: 'copy'

    cpus   8
    memory 12.GB
    time   2.h

    input:
    tuple val(meta), path(reads), path(vector)

    output:
    tuple val(meta), path("vector.sam"), emit: sam

    script:
    """
    python /app/aav_analyze.py align-vector \\
        --reads ${reads} \\
        --vector ${vector} \\
        --out-sam vector.sam \\
        --threads ${task.cpus}
    """
}
