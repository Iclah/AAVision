// CLASSIFY_BEST: assign each read to its best-scoring reference across all the
// per-reference BAMs, writing one FASTQ per reference plus <mode>_unmapped.fastq.
// The reference name is taken from each BAM's filename, so the result does not
// depend on channel ordering. .bai files are staged because pysam needs them.

process CLASSIFY_BEST {
    tag "${meta.id}:${mode}"

    publishDir { "${params.outdir}/${meta.id}/fastq" }, mode: 'copy', pattern: "*.fastq"

    cpus   2
    memory 8.GB
    time   1.h

    input:
    tuple val(meta), val(mode), path(bams), path(bais), path(source_fastq)

    output:
    tuple val(meta), val(mode), path("${mode}_*.fastq"), emit: classified

    script:
    """
    # Build the 'name=path' pairs for classify-best by deriving each
    # reference name from its BAM filename: alignments_<mode>_<ref>.bam -> <ref>
    bam_args=""
    for bam in ${bams}; do
        base=\$(basename "\$bam" .bam)
        # strip leading 'alignments_<mode>_'
        refname=\${base#alignments_${mode}_}
        bam_args="\$bam_args \$refname=\$bam"
    done

    python /app/aav_analyze.py classify-best \\
        --bams \$bam_args \\
        --source-fastq ${source_fastq} \\
        --out-prefix ${mode} \\
        --mode ${mode}
    """
}
