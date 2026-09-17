// ALIGN_REF: align reads to one reference, emit a sorted indexed BAM.
// Generic — used for both the hybrid pass (LQ vector reads) and the impurity
// pass (unmapped reads); the mode tag drives the output name.

process ALIGN_REF {
    tag "${meta.id}:${mode}:${ref_name}"

    publishDir { "${params.outdir}/${meta.id}/bam" }, mode: 'copy', pattern: "alignments_*.bam*"

    cpus   4
    memory 8.GB
    time   1.h

    input:
    tuple val(meta), path(reads), val(ref_name), path(ref), val(mode)

    output:
    tuple val(meta), val(mode), val(ref_name),
          path("alignments_${mode}_${ref_name}.bam"),
          path("alignments_${mode}_${ref_name}.bam.bai"), emit: bam

    script:
    """
    python /app/aav_analyze.py align-ref \\
        --reads ${reads} \\
        --ref ${ref} \\
        --ref-name ${ref_name} \\
        --out-bam alignments_${mode}_${ref_name}.bam \\
        --threads ${task.cpus}
    """
}
