#!/usr/bin/env nextflow

/*
============================================================
 AAVision — Nextflow pipeline
============================================================
 Samplesheet-driven analysis of AAV sequencing data.

 Usage:
   nextflow run main.nf -profile docker --input samplesheet.csv --outdir results

 Samplesheet columns (header required):
   sample    - unique sample name
   fastq     - path to a DIRECTORY containing the sample's FASTQ file(s)
   vector    - path to vector reference FASTA
   backbone  - path to backbone reference FASTA
   helper    - path to helper reference FASTA
   host      - path to host genome reference FASTA
   rep_cap   - (optional) path to rep/cap reference FASTA; leave blank if absent
   spike_in  - (optional) path to spike-in reference FASTA; leave blank if absent
============================================================
*/

include { PREPARE_FASTQ }  from './modules/local/prepare_fastq.nf'
include { ALIGN_VECTOR }   from './modules/local/align_vector.nf'
include { EXTRACT_VECTOR } from './modules/local/extract_vector.nf'
include { CLASSIFY_HQ_LQ } from './modules/local/classify_hq_lq.nf'
include { ALIGN_REF }      from './modules/local/align_ref.nf'
include { CLASSIFY_BEST }  from './modules/local/classify_best.nf'
include { COMPUTE_TRUNCATION } from './modules/local/compute_truncation.nf'
include { COMPUTE_COVERAGE }   from './modules/local/compute_coverage.nf'
include { COMPUTE_LENGTHS }    from './modules/local/compute_lengths.nf'
include { TAXONOMY }           from './modules/local/taxonomy.nf'
include { SUMMARIZE }          from './modules/local/summarize.nf'
include { GENERATE_REPORT }    from './modules/local/generate_report.nf'

// ------------------------------------------------------------
// Parse one samplesheet row into [ meta, fastq_dir, refs ]
// (a plain function — no statements at script top level)
// ------------------------------------------------------------
def parse_row(row) {
    if (!row.sample?.trim()) {
        error("Samplesheet row missing 'sample': ${row}")
    }
    if (!row.fastq?.trim()) {
        error("Sample ${row.sample}: missing 'fastq' directory")
    }

    // Validate required references using each (no for-loops in strict syntax)
    ['vector', 'backbone', 'helper', 'host'].each { req ->
        if (!row[req]?.trim()) {
            error("Sample ${row.sample}: missing required reference '${req}'")
        }
    }

    def meta = [ id: row.sample.trim() ]
    def fastq_dir = file(row.fastq.trim(), checkIfExists: true)

    def refs = [:]
    refs.vector   = file(row.vector.trim(),   checkIfExists: true)
    refs.backbone = file(row.backbone.trim(), checkIfExists: true)
    refs.helper   = file(row.helper.trim(),   checkIfExists: true)
    refs.host     = file(row.host.trim(),     checkIfExists: true)
    if (row.rep_cap?.trim()) {
        refs.rep_cap = file(row.rep_cap.trim(), checkIfExists: true)
    }
    if (row.spike_in?.trim()) {
        refs.spike_in = file(row.spike_in.trim(), checkIfExists: true)
    }

    return [ meta, fastq_dir, refs ]
}

// ------------------------------------------------------------
// Expand a [meta, reads, refs] tuple into one [meta, reads, ref_name,
// ref_path, mode] item per impurity reference, for the ALIGN_REF fan-out.
//
// hybrid  mode uses: backbone, helper, host, rep_cap   (NOT spike_in)
// impurity mode uses: backbone, helper, host, rep_cap, spike_in
// Optional refs (rep_cap, spike_in) are included only when present in refs.
// hybrid uses backbone/helper/host/rep_cap; impurity adds spike_in.
// ------------------------------------------------------------
def expand_refs(meta, reads, refs, mode) {
    def names = (mode == 'impurity')
        ? ['backbone', 'helper', 'host', 'rep_cap', 'spike_in']
        : ['backbone', 'helper', 'host', 'rep_cap']
    return names
        .findAll { name -> refs.containsKey(name) }
        .collect { name -> tuple(meta, reads, name, refs[name], mode) }
}

// ------------------------------------------------------------
// Size the TAXONOMY memory request from the Kraken2 DB on the HOST filesystem
// (computed in the workflow, before staging — the staged path inside the task
// isn't reliably enumerable at directive-evaluation time). Kraken2 loads the
// .k2d files (hash.k2d dominates) into RAM; sum them, x1.1 + 2.GB margin, floor
// 8.GB. Returns a memory string like '114 GB'. db is a Nextflow path object
// (file(params.kraken_db)) — use its native path API, not java.io.File, so
// non-native mounts (e.g. /mnt/c under WSL) resolve correctly.
// ------------------------------------------------------------
def compute_kraken_memory(db) {
    long bytes = 0L
    def found = []
    if (db.isDirectory()) {
        db.listFiles().each { f ->
            if (f.getName().endsWith('.k2d')) {
                bytes += f.size()
                found << "${f.getName()}=${f.size()}"
            }
        }
    }
    def floorBytes = 8L * 1024L * 1024L * 1024L
    def wanted = (bytes * 1.1) + (2L * 1024L * 1024L * 1024L)
    long finalBytes = (wanted < floorBytes ? floorBytes : wanted) as long
    def gb = Math.ceil(finalBytes / (1024.0 * 1024.0 * 1024.0)) as int
    log.info "TAXONOMY memory sizing: db=${db} .k2d files=[${found.join(', ')}] " +
             "total=${bytes} bytes -> request=${gb} GB" +
             (bytes == 0L ? "  (WARNING: no .k2d files found — check --kraken_db path)" : "")
    return "${gb} GB"
}

// ------------------------------------------------------------
// Main workflow
// ------------------------------------------------------------
workflow {
    main:

    if (!params.input) {
        error("ERROR: --input samplesheet.csv is required")
    }

    // ---- Reconcile parameter naming ----
    // Nextflow converts a kebab-case CLI flag (--min-length) to camelCase
    // (params.minLength), which does NOT match the snake_case config default
    // (params.min_length). Prefer an explicitly-supplied CLI value (camelCase
    // variant) and otherwise fall back to the snake_case config value, so that
    // --min-length, --minLength, and --min_length all work. snake_case is the
    // canonical name used internally.
    def min_length = params.containsKey('minLength') ? params.minLength : params.min_length

    // Coverage threshold for HQ/LQ split. CLI --quality is a percentage (e.g. 90);
    // the classifier wants a 0-1 fraction, so divide by 100. 'quality' is a
    // single-word param so it needs no kebab/camel reconciliation.
    def min_coverage = (params.quality as String).toBigDecimal() / 100.0

    samples_ch = channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)
        .map { row -> parse_row(row) }

    // Keep a [meta, refs] channel so reference paths can be rejoined to each
    // sample after PREPARE_FASTQ. meta is the grouping key throughout.
    refs_ch = samples_ch.map { meta, fastq_dir, refs -> tuple(meta, refs) }

    // --- Step 1: prepare FASTQ (merge + optional length filter) ---
    // PREPARE_FASTQ needs only meta + fastq_dir.
    PREPARE_FASTQ(
        samples_ch.map { meta, fastq_dir, refs -> tuple(meta, fastq_dir) },
        min_length
    )

    // --- Step 2: align prepared reads to the vector reference ---
    // Join prepared reads with the per-sample refs, then pass
    // [meta, reads, vector] to ALIGN_VECTOR.
    align_vector_in = PREPARE_FASTQ.out.reads
        .join(refs_ch)
        .map { meta, reads, refs -> tuple(meta, reads, refs.vector) }

    ALIGN_VECTOR(align_vector_in)

    // --- Step 3: extract mapped/unmapped reads and build the vector BAM ---
    // EXTRACT_VECTOR needs the SAM (from ALIGN_VECTOR) AND the prepared reads
    // (from PREPARE_FASTQ); join both on meta -> [meta, sam, reads].
    extract_vector_in = ALIGN_VECTOR.out.sam
        .join(PREPARE_FASTQ.out.reads)
        .map { meta, sam, reads -> tuple(meta, sam, reads) }

    EXTRACT_VECTOR(extract_vector_in)

    // --- Step 4: split vector reads into HQ / LQ by coverage ---
    // Needs the vector SAM (from ALIGN_VECTOR) and vector.fastq (from
    // EXTRACT_VECTOR); join both on meta -> [meta, sam, vector_fastq].
    classify_hq_lq_in = ALIGN_VECTOR.out.sam
        .join(EXTRACT_VECTOR.out.vector_fastq)
        .map { meta, sam, vector_fastq -> tuple(meta, sam, vector_fastq) }

    CLASSIFY_HQ_LQ(classify_hq_lq_in, min_coverage)

    // --- Steps 5+7: combined fan-out — align LQ reads (hybrid) AND unmapped
    // reads (impurity) to each impurity reference. A process can only be
    // invoked ONCE per workflow in DSL2, so we merge both passes into one
    // channel (each item already tagged with its mode) and call ALIGN_REF
    // once. This also maximizes parallelism: all hybrid + impurity alignments
    // run as independent tasks.
    hybrid_align_in = CLASSIFY_HQ_LQ.out.lq
        .join(refs_ch)
        .flatMap { meta, lq, refs -> expand_refs(meta, lq, refs, 'hybrid') }

    impurity_align_in = EXTRACT_VECTOR.out.unmapped
        .join(refs_ch)
        .flatMap { meta, unmapped, refs -> expand_refs(meta, unmapped, refs, 'impurity') }

    align_ref_in = hybrid_align_in.mix(impurity_align_in)

    ALIGN_REF(align_ref_in)

    // --- Steps 6+8: combined gather + best-alignment classification ---
    // Group BAMs per [id, mode] (mode keeps hybrid and impurity separate),
    // then join each group with its source reads. The source differs by mode:
    //   hybrid   -> LQ reads (CLASSIFY_HQ_LQ.out.lq)
    //   impurity -> unmapped reads (EXTRACT_VECTOR.out.unmapped)
    // We build a [meta_with_mode, source] channel for both and join on that.
    grouped_bams = ALIGN_REF.out.bam
        .map { meta, mode, ref_name, bam, bai ->
            tuple([id: meta.id, mode: mode], bam, bai)
        }
        .groupTuple()

    // Source-reads channel keyed by [id, mode] so it can be joined to the
    // grouped BAMs regardless of pass.
    hybrid_src = CLASSIFY_HQ_LQ.out.lq
        .map { meta, lq -> tuple([id: meta.id, mode: 'hybrid'], lq) }
    impurity_src = EXTRACT_VECTOR.out.unmapped
        .map { meta, unmapped -> tuple([id: meta.id, mode: 'impurity'], unmapped) }
    source_reads = hybrid_src.mix(impurity_src)

    classify_best_in = grouped_bams
        .join(source_reads)
        .map { key, bams, bais, src ->
            tuple([id: key.id], key.mode, bams, bais, src)
        }

    CLASSIFY_BEST(classify_best_in)

    CLASSIFY_BEST.out.classified.view { meta, mode, fastqs ->
        "CLASSIFY_BEST [${meta.id}:${mode}] -> ${fastqs}"
    }

    // --- Step 9: truncation hotspots on the vector BAM ---
    // EXTRACT_VECTOR.out.bam is [meta, bam, bai] — pass straight through
    // (pysam .fetch needs the .bai staged alongside).
    COMPUTE_TRUNCATION(EXTRACT_VECTOR.out.bam)

    COMPUTE_TRUNCATION.out.truncation.view { meta, csv, json ->
        "TRUNCATION [${meta.id}] -> ${csv}"
    }

    // --- Step 10: per-category coverage ---
    // Coverage is computed for rAAV (vector BAM) plus the helper/rep_cap/
    // spike_in IMPURITY BAMs — NOT backbone, NOT host.
    // Build a [meta, category, bam, bai] channel from two sources and run the
    // generic COMPUTE_COVERAGE once per category.

    // (a) vector / rAAV coverage from EXTRACT_VECTOR
    coverage_vector = EXTRACT_VECTOR.out.bam
        .map { meta, bam, bai -> tuple(meta, 'vector', bam, bai) }

    // (b) impurity coverage for helper, rep_cap, spike_in only
    coverage_impurity = ALIGN_REF.out.bam
        .filter { meta, mode, ref_name, bam, bai ->
            mode == 'impurity' && ref_name in ['helper', 'rep_cap', 'spike_in']
        }
        .map { meta, mode, ref_name, bam, bai -> tuple(meta, ref_name, bam, bai) }

    coverage_in = coverage_vector.mix(coverage_impurity)

    COMPUTE_COVERAGE(coverage_in)

    COMPUTE_COVERAGE.out.coverage.view { meta, category, json ->
        "COVERAGE [${meta.id}:${category}] -> ${json}"
    }

    // --- Step 11: read-length distributions ---
    // One combined JSON across categories: vector + impurity helper/rep_cap/
    // host/spike_in/unmapped (backbone excluded; handled in the process).
    // Pass vector.fastq and the impurity FASTQ list as separate inputs;
    // COMPUTE_LENGTHS derives impurity category labels from filenames.
    impurity_fastqs = CLASSIFY_BEST.out.classified
        .filter { meta, mode, fastqs -> mode == 'impurity' }
        .map { meta, mode, fastqs -> tuple(meta, fastqs) }

    lengths_in = EXTRACT_VECTOR.out.vector_fastq
        .join(impurity_fastqs)
        .map { meta, vector_fq, impurity_fqs ->
            tuple(meta, vector_fq, impurity_fqs)
        }

    COMPUTE_LENGTHS(lengths_in)

    COMPUTE_LENGTHS.out.lengths.view { meta, json ->
        "LENGTHS [${meta.id}] -> ${json}"
    }

    // --- Step 12: taxonomy (OPTIONAL, gated on --kraken_db) ---
    // Only runs when a Kraken2 DB is provided. Classifies the unmapped
    // impurity reads. The unmapped FASTQ is selected from the impurity
    // CLASSIFY_BEST output by filename (impurity_unmapped.fastq).
    //
    // For the report we need a taxonomy entry for EVERY sample (real summary
    // or a NO_FILE placeholder). We do NOT use join(remainder:true) against an
    // empty channel — that has a documented bug (nextflow #1015) where an empty
    // right channel makes join emit only the key, collapsing the tuple. Instead
    // we build a placeholder channel keyed by meta and (when taxonomy ran)
    // override it with the real summary via mix + groupTuple-free join.
    no_file = file("${projectDir}/assets/NO_FILE")

    // Placeholder baseline: every sample gets NO_FILE. If taxonomy runs, real
    // summaries override the placeholder. This works whether taxonomy is off,
    // on-for-all-samples, or on-but-skipped-for-some — and avoids the empty-
    // channel join bug (nextflow #1015) entirely.
    tax_placeholder = refs_ch.map { meta, refs -> tuple(meta, no_file) }

    if (params.kraken_db) {
        unmapped_ch = CLASSIFY_BEST.out.classified
            .filter { meta, mode, fastqs -> mode == 'impurity' }
            .map { meta, mode, fastqs ->
                def flist = fastqs instanceof List ? fastqs : [fastqs]
                def unmapped = flist.find { f -> f.name == 'impurity_unmapped.fastq' }
                tuple(meta, unmapped)
            }
            .filter { meta, unmapped -> unmapped != null }

        // Resolve the memory request ONCE, on the host DB path, before staging.
        // 'auto' -> size from the DB .k2d files; otherwise use the explicit value.
        def kraken_mem = (params.kraken_memory == 'auto')
            ? compute_kraken_memory(file(params.kraken_db))
            : params.kraken_memory

        TAXONOMY(unmapped_ch, file(params.kraken_db), kraken_mem)

        // Real summaries take precedence; placeholder fills any sample TAXONOMY
        // didn't emit for. mix + groupTuple + pick-real keeps it deterministic.
        taxonomy_ch = TAXONOMY.out.summary
            .mix(tax_placeholder)
            .groupTuple()
            .map { meta, files ->
                def real = files.find { f -> f.name != 'NO_FILE' }
                tuple(meta, real != null ? real : no_file)
            }
    } else {
        taxonomy_ch = tax_placeholder
    }

    // --- Step 13: SUMMARIZE — aggregate classification stats ---
    // Stage HQ, LQ, all hybrid FASTQs and all impurity FASTQs into one task.
    // optional_refs (which of rep_cap/spike_in this sample has) is derived
    // per-sample from its refs map and threaded through the channel.
    hybrid_fastqs_ch = CLASSIFY_BEST.out.classified
        .filter { meta, mode, fastqs -> mode == 'hybrid' }
        .map { meta, mode, fastqs -> tuple(meta, fastqs) }

    impurity_fastqs_for_sum = CLASSIFY_BEST.out.classified
        .filter { meta, mode, fastqs -> mode == 'impurity' }
        .map { meta, mode, fastqs -> tuple(meta, fastqs) }

    optional_refs_ch = refs_ch
        .map { meta, refs ->
            def opt = ['rep_cap', 'spike_in'].findAll { name -> refs.containsKey(name) }
            tuple(meta, opt)
        }

    summarize_in = CLASSIFY_HQ_LQ.out.hq
        .join(CLASSIFY_HQ_LQ.out.lq)
        .join(hybrid_fastqs_ch)
        .join(impurity_fastqs_for_sum)
        .join(optional_refs_ch)
        .map { meta, hq, lq, hybrid_fqs, impurity_fqs, opt ->
            tuple(meta, hq, lq, hybrid_fqs, impurity_fqs, opt)
        }

    SUMMARIZE(summarize_in)

    // --- Step 14: GENERATE_REPORT — the fan-in ---
    // Gather all per-category coverage JSONs into a list per sample.
    coverage_grouped = COMPUTE_COVERAGE.out.coverage
        .map { meta, category, json -> tuple(meta, json) }
        .groupTuple()

    // truncation: keep only the JSON (output is [meta, csv, json]).
    truncation_json_ch = COMPUTE_TRUNCATION.out.truncation
        .map { meta, csv, json -> tuple(meta, json) }

    // Assemble [meta, stats, truncation, coverage_list, lengths, taxonomy].
    // taxonomy_ch is complete (one entry per sample, real or NO_FILE), so a
    // plain join is safe.
    report_in = SUMMARIZE.out.stats
        .join(truncation_json_ch)
        .join(coverage_grouped)
        .join(COMPUTE_LENGTHS.out.lengths)
        .join(taxonomy_ch)
        .map { meta, stats, trunc, cov_list, lengths, tax ->
            tuple(meta, stats, trunc, cov_list, lengths, tax)
        }

    GENERATE_REPORT(report_in)

    GENERATE_REPORT.out.report.view { meta, html ->
        "REPORT [${meta.id}] -> ${html}"
    }
}
