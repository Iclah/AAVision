// PREPARE_FASTQ: merge a sample's FASTQ files and optionally length-filter them
// into a single prepared FASTQ.

process PREPARE_FASTQ {
    tag "${meta.id}"

    // Publish only the canonical prepared FASTQ and the filter stats. Without
    // an explicit pattern, the merge/filter intermediate (merged.filtered.fastq.gz)
    // left in the work dir would also be published, duplicating prepared.fastq.gz.
    publishDir { "${params.outdir}/${meta.id}/fastq" }, mode: 'copy',
                pattern: "{prepared.fastq.gz,length_filter_stats.json}"

    cpus   2
    memory 4.GB
    time   1.h

    input:
    tuple val(meta), path(fastq_dir)
    val min_length

    output:
    tuple val(meta), path("prepared.fastq.gz"), emit: reads
    path "length_filter_stats.json",            optional: true, emit: filter_stats

    script:
    """
    # prepare-fastq discovers FASTQ files under the sample directory
    # (one level deep), merges them if there are several, and applies the
    # optional length filter.
    #
    # --work-dir .   : write merge/filter intermediates into the task work
    #                  dir (writable), never into the staged input directory.
    # --out          : copy the final prepared FASTQ to a predictable name
    #                  for the output channel.
    # --out-filter-stats : surface the length-filter stats as an output file.

    python /app/aav_analyze.py prepare-fastq \\
        --sample-dir ${fastq_dir} \\
        --sample-name ${meta.id} \\
        --min-length ${min_length} \\
        --work-dir . \\
        --out prepared.fastq.gz \\
        --out-filter-stats length_filter_stats.json
    """
}
