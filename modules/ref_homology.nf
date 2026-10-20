// REF_HOMOLOGY: find the regions of one hybrid reference (backbone, helper,
// host, rep_cap) that share sequence with the vector, as BED (0-based
// half-open, with a '# searched:' header). CLASSIFY_BEST concatenates a
// sample's BEDs and ignores alignments inside these regions in the hybrid pass.
// The vector is staged in its own folder so it cannot clash with a reference
// that has the same file name.

process REF_HOMOLOGY {
    tag "${meta.id}:${ref_name}"

    publishDir { "${params.outdir}/${meta.id}/homology" }, mode: 'copy'

    // Same footprint as ALIGN_REF: indexing the host genome dominates
    cpus   4
    memory 8.GB
    time   1.h

    input:
    tuple val(meta), path(vector, stageAs: 'query/*'), val(ref_name), path(ref)
    val preset

    output:
    tuple val(meta), val(ref_name), path("shared_regions_${ref_name}.bed"), emit: bed

    script:
    """
    python /app/aav_analyze.py ref-homology \\
        --vector ${vector} \\
        --ref ${ref_name}=${ref} \\
        --preset ${preset} \\
        --threads ${task.cpus} \\
        --out-bed shared_regions_${ref_name}.bed
    """
}
