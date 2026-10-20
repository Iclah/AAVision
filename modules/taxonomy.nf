// TAXONOMY (optional; gated on params.kraken_db): classify unmapped reads with
// Kraken2 + Bracken and write taxonomy_summary.json for the report.
// Kraken2 loads the whole database into RAM, so the memory request is sized to
// the database in the workflow and passed in as `mem`.

process TAXONOMY {
    tag "${meta.id}"

    publishDir { "${params.outdir}/${meta.id}/stats/taxonomy" }, mode: 'copy'

    cpus   params.kraken_cpus
    memory { mem }
    time   params.kraken_time

    input:
    tuple val(meta), path(unmapped)
    path db
    val mem

    output:
    tuple val(meta), path("taxonomy_summary.json"), emit: summary
    path "kraken2_output.txt", optional: true
    path "kraken2_report.txt", optional: true
    path "bracken_*.txt",      optional: true

    script:
    // Only pass --read-length when a Bracken read length is configured; when
    // null the Python side falls back to the database maximum.
    def rl_arg = params.bracken_read_length ? "--read-length ${params.bracken_read_length}" : ""
    """
    python /app/aav_analyze.py taxonomy \\
        --reads ${unmapped} \\
        --db ${db} \\
        --out-dir . \\
        --out-summary taxonomy_summary.json \\
        --confidence ${params.kraken_confidence} \\
        --level ${params.bracken_level} \\
        --threshold ${params.bracken_threshold} \\
        ${rl_arg} \\
        --threads ${task.cpus}
    """
}
