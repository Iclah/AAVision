// EXTRACT_VECTOR: from the vector SAM, emit mapped reads (vector.fastq), a
// sorted indexed BAM, and the unmapped reads. Needs the SAM and the prepared
// reads (unmapped = input minus mapped).

process EXTRACT_VECTOR {
    tag "${meta.id}"

    publishDir { "${params.outdir}/${meta.id}/bam" },   mode: 'copy', pattern: "vector.bam*"
    publishDir { "${params.outdir}/${meta.id}/fastq" }, mode: 'copy', pattern: "*.fastq"

    cpus   2
    memory 4.GB
    time   1.h

    input:
    tuple val(meta), path(sam), path(reads)

    output:
    tuple val(meta), path("vector.bam"), path("vector.bam.bai"), emit: bam
    tuple val(meta), path("vector.fastq"),                       emit: vector_fastq
    tuple val(meta), path("vector_unmap.fastq"),                 emit: unmapped

    script:
    """
    python /app/aav_analyze.py extract-vector \\
        --sam ${sam} \\
        --reads ${reads} \\
        --out-fastq vector.fastq \\
        --out-bam vector.bam \\
        --out-unmap vector_unmap.fastq
    """
}
