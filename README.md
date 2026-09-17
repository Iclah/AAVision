# AAVision

A Nextflow workflow for quality control and characterization of recombinant
adeno-associated virus (rAAV) preparations from Oxford Nanopore sequencing.

It aligns reads to the vector and reference sequences, classifies them
(rAAV / hybrid / impurity), quantifies truncation hotspots and read-length
distributions, optionally performs taxonomic classification of unmapped reads
with Kraken2/Bracken, and produces an interactive HTML report.

## Requirements

- [Nextflow](https://www.nextflow.io/) (DSL2; developed against 26.04)
- [Docker](https://www.docker.com/) (the workflow runs each process in a
  container; the image bundles minimap2, samtools, seqkit, kraken2, bracken and
  the Python analysis code)

No other local installation is required — the tools live in the container.

## Quick start

You need [Nextflow](https://www.nextflow.io/) and [Docker](https://www.docker.com/)
installed and running (see [Requirements](#requirements)). Everything else — the
pipeline code and the analysis container — is fetched automatically.

Create a samplesheet (see [Samplesheet](#samplesheet)), then run the pipeline
straight from GitHub:

```bash
nextflow run Iclah/AAVision -r v1.0.0 -profile docker \
  --input samplesheet.csv \
  --outdir results \
  --min-length 300
```

`-r v1.0.0` pins a specific release. On the first run, Nextflow downloads the
pipeline and pulls the container image (`ghcr.io/iclah/aavision:1.0.0`)
automatically; there is no build or manual install step. The interactive report
is written to `results/<sample>/analysis_report.html`.

### Running from a local clone

If you prefer to clone the repository (or want to modify the workflow), run
`main.nf` directly instead:

```bash
git clone https://github.com/Iclah/AAVision.git
cd AAVision
nextflow run main.nf -profile docker \
  --input samplesheet.csv --outdir results --min-length 300
```

### Developing against a local image

By default the workflow uses the published image. To build and use a local
image instead, override the container:

```bash
docker build -t aavision:dev -f AAV_analyzer.dockerfile .
nextflow run main.nf -profile docker --input samplesheet.csv \
  --outdir results --min-length 300 --container aavision:dev
```

## Samplesheet

A comma-separated file with one row per sample. The `fastq` column points at a
**directory** of FASTQ files (e.g. an ONT `fastq_pass` folder); they are merged
automatically. Reference columns are explicit per row.

```csv
sample,fastq,vector,backbone,helper,host,rep_cap,spike_in
vKUL156,/data/run1/vKUL156/fastq_pass,/data/refs/vector.fasta,/data/refs/backbone.fasta,/data/refs/helper.fasta,/data/refs/host.fasta,/data/refs/rep_cap.fasta,
```

| Column | Required | Description |
|--------|----------|-------------|
| `sample` | yes | Unique sample identifier |
| `fastq` | yes | Directory containing the sample's FASTQ files |
| `vector` | yes | rAAV vector reference FASTA |
| `backbone` | yes | Plasmid backbone reference FASTA |
| `helper` | yes | Helper plasmid reference FASTA |
| `host` | yes | Host genome reference FASTA |
| `rep_cap` | no | Rep/Cap reference FASTA (leave blank to skip) |
| `spike_in` | no | Spike-in reference FASTA (leave blank to skip) |

The file must be comma-separated. Some spreadsheet applications export
semicolon-separated CSVs on non-US locales — check the delimiter if the
samplesheet fails to parse.

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | (required) | Path to the samplesheet CSV |
| `--outdir` | `results` | Directory for published outputs |
| `--quality` | `90` | Minimum alignment coverage (%) for high-quality vector reads |
| `--min-length` | `0` | Minimum read length filter in bp (`0` = off) |
| `--kraken_db` | `null` | Path to a Kraken2 database directory. If unset, taxonomy is skipped |
| `--kraken_confidence` | `0.05` | Kraken2 confidence threshold |
| `--bracken_level` | `S` | Bracken taxonomic level (D/P/C/O/F/G/S) |
| `--bracken_threshold` | `10` | Bracken minimum read count per taxon |
| `--bracken_read_length` | `null` | Bracken read length (`-r`). Must match a length the Bracken DB was built for. If unset, the DB maximum is used |
| `--kraken_memory` | `auto` | TAXONOMY memory request. `auto` sizes it from the DB `.k2d` files; override with an explicit value (e.g. `128.GB`) |
| `--kraken_cpus` | `8` | Threads for Kraken2 |
| `--kraken_time` | `12.h` | Time limit for the taxonomy step |

### Resource limits

By default the workflow uses whatever the machine provides. To cap resources
(e.g. on a laptop or a memory-limited VM), uncomment and set `resourceLimits` in
`nextflow.config`. This is documented inline in that file.

### Taxonomy (optional)

Taxonomic classification of unmapped reads runs only when `--kraken_db` is
supplied. Kraken2 memory-maps the entire database into RAM, so this step
requires a machine sized to the database (a standard database is 100+ GB). The
memory request is sized automatically from the database; on a machine that
cannot provide it the step fails clearly rather than silently under-provisioning.

To reproduce a specific read length (e.g. matching the read-length peak of your
data), set `--bracken_read_length` to a value the Bracken database was built for.

## Output

```
results/<sample>/
├── analysis_report.html          # interactive report — start here
├── fastq/                        # classified reads
│   ├── prepared.fastq.gz         # merged + length-filtered input
│   ├── vector_HQ.fastq           # high-quality vector reads
│   ├── vector_LQ.fastq           # low-quality vector reads
│   ├── hybrid_*.fastq            # hybrid reads by reference
│   └── impurity_*.fastq          # impurity reads by reference
├── bam/                          # alignments (IGV-compatible)
└── stats/                        # JSON/CSV statistics, coverage, truncation,
                                  # read lengths, and taxonomy (if enabled)
```

## The report

The HTML report includes overall classification statistics (rAAV / hybrid /
impurity), a read-distribution table and pie chart, read-length histograms,
per-reference coverage, a truncation-hotspot plot, and — when taxonomy is
enabled — a taxonomic breakdown of unmapped reads. Charts render client-side
with Plotly, so viewing the report requires an internet connection.

## Method

1. **Vector alignment** — reads are aligned to the vector reference with
   minimap2; mapped and unmapped reads are separated.
2. **Quality classification** — mapped reads are split into high- and
   low-quality by alignment coverage (default ≥90%).
3. **Hybrid detection** — low-quality vector reads are aligned to the other
   references and assigned to the best-matching one.
4. **Impurity detection** — unmapped reads are aligned to all references and
   classified by best alignment.
5. **Truncation analysis** — read termini are counted, normalized by coverage,
   smoothed, and reported as enrichment along the vector genome.
6. **Read-length analysis** — length distributions are binned per category.
7. **Taxonomy (optional)** — unmapped reads are classified with Kraken2 and
   abundances re-estimated with Bracken.

## Repository layout

```
main.nf                  # workflow entry point
nextflow.config          # parameters and configuration
modules/local/           # one Nextflow process per file
assets/                  # samplesheet example, placeholder files
env.yaml                 # pinned container dependencies
AAV_analyzer.dockerfile  # container definition
aav_analyze.py           # CLI the processes invoke inside the container
aav_analyzer/            # analysis library (alignment, classification,
                         # truncation, taxonomy, report)
```

## License

MIT — see the LICENSE file.
