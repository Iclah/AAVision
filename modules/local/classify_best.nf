// CLASSIFY_BEST: assign each read to its best-scoring reference across all the
// per-reference BAMs, writing one FASTQ per reference plus <mode>_unmapped.fastq.
// The reference name is taken from each BAM's filename, so the result does not
// depend on channel ordering. .bai files are staged because pysam needs them.
//
// In hybrid mode, shared_beds are the sample's REF_HOMOLOGY BEDs: alignments
// inside vector-shared regions do not count, and the rejected reads are listed
// in hybrid_shared_only.tsv with counts in hybrid_filter_stats.json. In
// impurity mode shared_beds is the NO_FILE placeholder and is not used.

process CLASSIFY_BEST {
    tag "${meta.id}:${mode}"

    publishDir { "${params.outdir}/${meta.id}/fastq" }, mode: 'copy',
        pattern: "{*.fastq,*_shared_only.tsv,*_filter_stats.json}"

    cpus   2
    memory 8.GB
    time   1.h

    input:
    tuple val(meta), val(mode), path(bams), path(bais), path(source_fastq), path(shared_beds, stageAs: 'homology/*')
    val min_unique_bases

    output:
    tuple val(meta), val(mode), path("${mode}_*.fastq"), emit: classified
    tuple val(meta), path("${mode}_filter_stats.json"), optional: true, emit: filter_stats
    tuple val(meta), path("${mode}_shared_only.tsv"), optional: true, emit: shared_only

    script:
    def shared = mode == 'hybrid'
        ? "--shared-bed shared_regions.bed --min-unique-bases ${min_unique_bases}"
        : ""
    def concat = mode == 'hybrid' ? "cat homology/*.bed > shared_regions.bed" : ""
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

    ${concat}

    python /app/aav_analyze.py classify-best \\
        --bams \$bam_args \\
        --source-fastq ${source_fastq} \\
        --out-prefix ${mode} \\
        --mode ${mode} ${shared}
    """
}
