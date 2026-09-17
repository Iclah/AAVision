// GENERATE_REPORT: build the HTML report from the computed inputs. The per-
// category coverage JSONs are merged into one object in the script before the
// report is generated. Taxonomy is included only when a summary is present.

process GENERATE_REPORT {
    tag "${meta.id}"

    publishDir { "${params.outdir}/${meta.id}" }, mode: 'copy', pattern: "analysis_report.html"

    cpus   2
    memory 4.GB
    time   1.h

    input:
    tuple val(meta), path(stats), path(truncation_json), path(coverage_jsons), path(lengths), path(taxonomy_summary)

    output:
    tuple val(meta), path("analysis_report.html"), emit: report

    script:
    // taxonomy_summary may be a real file or a placeholder (NO_FILE). Only
    // pass --taxonomy when a genuine summary is present.
    def tax_arg = (taxonomy_summary.name != 'NO_FILE') ? "--taxonomy ${taxonomy_summary}" : ""
    """
    # Merge per-category coverage JSONs into one {category: {pos: depth}} object.
    # Category name is derived from each filename: coverage_<cat>.json -> <cat>.
    python3 - <<'PY'
import json, glob, os
combined = {}
for path in sorted(glob.glob('coverage_*.json')):
    cat = os.path.basename(path)[len('coverage_'):-len('.json')]
    with open(path) as fh:
        combined[cat] = json.load(fh)
with open('combined_coverage.json', 'w') as fh:
    json.dump(combined, fh)
PY

    python /app/aav_analyze.py generate-report \\
        --stats ${stats} \\
        --truncation ${truncation_json} \\
        --coverage combined_coverage.json \\
        --lengths ${lengths} \\
        ${tax_arg} \\
        --out analysis_report.html
    """
}
