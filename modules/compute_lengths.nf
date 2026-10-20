// COMPUTE_LENGTHS: read-length distributions per category, as one combined JSON.
// Categories are derived from the FASTQ filenames; backbone is excluded and
// empty files are skipped.

process COMPUTE_LENGTHS {
    tag "${meta.id}"

    publishDir { "${params.outdir}/${meta.id}/stats" }, mode: 'copy', pattern: "read_length_distributions.json"

    cpus   2
    memory 4.GB
    time   1.h

    input:
    tuple val(meta), path(vector_fastq), path(impurity_fastqs)

    output:
    tuple val(meta), path("read_length_distributions.json"), emit: lengths

    script:
    """
    # Build name=path pairs. Always include vector. For each impurity FASTQ,
    # derive its category from the filename (impurity_<cat>.fastq -> <cat>),
    # skip backbone, and skip empty files.
    fastq_args="vector=${vector_fastq}"
    for fq in ${impurity_fastqs}; do
        base=\$(basename "\$fq" .fastq)
        cat=\${base#impurity_}
        # exclude backbone from length analysis
        if [ "\$cat" = "backbone" ]; then
            continue
        fi
        # skip empty files
        if [ ! -s "\$fq" ]; then
            continue
        fi
        fastq_args="\$fastq_args \$cat=\$fq"
    done

    python /app/aav_analyze.py compute-lengths \\
        --fastqs \$fastq_args \\
        --out-json read_length_distributions.json
    """
}
