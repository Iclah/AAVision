// CLASSIFY_HQ_LQ: split mapped vector reads into high- and low-quality by
// per-read alignment coverage. Needs the vector SAM and vector.fastq.

process CLASSIFY_HQ_LQ {
    tag "${meta.id}"

    publishDir { "${params.outdir}/${meta.id}/fastq" }, mode: 'copy', pattern: "vector_*.fastq"

    cpus   2
    memory 4.GB
    time   1.h

    input:
    tuple val(meta), path(sam), path(vector_fastq)
    val min_coverage

    output:
    tuple val(meta), path("vector_HQ.fastq"), emit: hq
    tuple val(meta), path("vector_LQ.fastq"), emit: lq

    script:
    """
    python /app/aav_analyze.py classify-hq-lq \\
        --sam ${sam} \\
        --vector-fastq ${vector_fastq} \\
        --out-hq vector_HQ.fastq \\
        --out-lq vector_LQ.fastq \\
        --min-coverage ${min_coverage}
    """
}
