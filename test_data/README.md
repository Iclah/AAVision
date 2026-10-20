# Test dataset for `-profile test`

Used by `nextflow run ... -profile test,docker` to check an installation.
It is a smoke test, not a validation: the host reference is a subset, so a
few reads can align differently than against the full genome, and the
numbers differ slightly from a full-host run.

| File | Content |
|------|---------|
| `samplesheet.csv` | one sample; paths use `${projectDir}`, which `main.nf` replaces with the pipeline directory |
| `reads/` | test_reads.fastq.gz, copied unchanged |
| `references/vector,backbone,helper,rep_cap.fasta` | copied unchanged |
| `references/host_subset.fasta` | regions of GRCh38 that the reads or the vector align to |
| `references/host_subset_regions.bed` | those regions (0-based, half-open); column 4 is the contig name in `host_subset.fasta` |

## How the host subset was made

All reads and the vector were aligned to the full host; every alignment
(primary, secondary, supplementary) was padded by 5000 bp per side,
clipped to contig ends and merged.

10439 host alignment records on 25 contigs -> 164 regions, 1,813,321 bp (pad 5000 bp)

Commands (inside `ghcr.io/iclah/aavision:1.0.2`):

```
minimap2 -x map-ont -I 4G host.fasta test_reads.fastq.gz (PAF from an earlier run of this command)
minimap2 -x asm20 --secondary=yes -I 4G host.fasta vector.fasta (PAF from an earlier run of this command)
samtools faidx -r regions.txt host.fasta
```

Tools: minimap2 2.31-r1302; samtools 1.23.1
Generated 2026-10-20 by `scripts/make_test_profile_data.py` at commit 9d2a6f1.
