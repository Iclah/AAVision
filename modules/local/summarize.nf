// SUMMARIZE: aggregate classification counts from the categorized FASTQs into
// classification_stats.json and summary.txt. All FASTQs are staged into the
// task dir and counted from there.

process SUMMARIZE {
    tag "${meta.id}"

    publishDir { "${params.outdir}/${meta.id}/stats" }, mode: 'copy', pattern: "{classification_stats.json,summary.txt}"

    cpus   2
    memory 4.GB
    time   1.h

    input:
    tuple val(meta), path(hq), path(lq), path(hybrid_fastqs), path(impurity_fastqs), val(optional_refs)

    output:
    tuple val(meta), path("classification_stats.json"), emit: stats
    path "summary.txt", emit: summary_txt

    script:
    def opt = optional_refs ? "--optional-refs ${optional_refs.join(' ')}" : ""
    """
    python /app/aav_analyze.py summarize \\
        --fastq-dir . \\
        ${opt} \\
        --out-json classification_stats.json \\
        --out-txt summary.txt
    """
}
