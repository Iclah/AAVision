# Design decisions

Non-obvious design decisions, newest first. Check here before changing
behaviour that one of these entries covers.

---

## 2026-10-16 — HQ uses primary + supplementary read coverage; snapbacks stay LQ

**Decision.** A vector read is HQ if both hold:
1. Its primary and supplementary vector records together cover at least
   `min_coverage` (`-q`, default 90%) of the full read. Each record's covered
   read interval (leading H+S clip to clip + M/I/=/X bases, flipped for
   minus-strand records) is merged with the others, so overlapping pieces
   count once. The denominator is the full read length including hard clips,
   so the result is the same with or without minimap2 `-Y`. Secondary records
   are ignored.
2. It is not a snapback: no two of its records lie on opposite strands with
   vector reference spans overlapping by at least `min_snapback_overlap`
   (default 200 bp).

HQ reads are copied from vector.fastq (the primary records, full read, original
orientation). LQ is still vector.fastq minus HQ.

**Supersedes** the ecc6b65 rule (HQ decided from the primary record alone).

**Why.** Primary-only made reads whose pieces are all vector, split over two
records, LQ, while the vector fully explains them. Merging the records fixes
that. Without an extra rule, merging would make snapbacks (a vector region
followed by its own reverse complement) HQ, which is what ecc6b65 set out to
prevent; snapbacks are not intact genomes and belong in LQ.

**Alternatives considered.**
- *Primary only (ecc6b65).* Simple, but calls all-vector split reads LQ.
- *Merged coverage with snapbacks HQ* (the original task spec). Treats a read
  as HQ whenever vector explains it, whatever its structure; rejected so that
  snapbacks keep counting as LQ.
- *Overlap on any strand.* Would also make same-strand tandem repeats
  (head-to-tail concatemers) LQ; kept to opposite strands to match the
  snapback definition.
- *Any overlap (≥ 1 bp).* A short flip/flop ITR supplementary (~145 bp) on a
  full genome would make it LQ; 200 bp clears the ITR length.

**Consequences (flagged).**
- Rearranged genomes whose pieces don't overlap on the vector (inversions,
  permutations, junctions between distant vector regions) are HQ when they
  cover ≥ 90% of the read.
- Same-strand tandem repeats are HQ when coverage ≥ 90%.
- For scAAV, vector.fasta must be the full self-complementary construct; if it
  holds only one half, every full sc read is a snapback and becomes LQ.
- Foldbacks overlapping < 200 bp do not make a read LQ.

---

## 2026-10-16 — Hybrid calls ignore reference regions shared with the vector

**Decision.** Before the hybrid pass, vector.fasta is aligned to each other
reference (backbone, helper, host, rep_cap) with `minimap2 -x asm20
--secondary=yes`. The reference intervals it hits are the *shared regions*
(BED, 0-based half-open, strand ignored, merged per contig). A low-quality
(LQ) vector read is assigned to a reference only if its alignments to that
reference cover at least `min_unique_bases` (default 100) reference bases
outside the shared regions, summed over the read's primary and supplementary
records. Among eligible references the highest AS wins, as before. Reads with
no eligible reference stay pure LQ (`hybrid_unmapped`). The impurity pass is
unchanged.

**Why.** Before this change any alignment to another reference made an LQ
read a hybrid. When the vector shares sequence with another reference (a
backbone plasmid that includes the cassette, a promoter or polyA also present
in the helper, a human-derived transgene in the host genome), pure truncated
or snapback vector reads align there and were counted as hybrids, deflating
pure LQ and inflating hybrids.

**Alternatives considered.**
- *N-mask the shared regions in copies of the references.* Leaves the
  classifier untouched, but writes a masked copy of the host genome (~3 GB)
  per sample, and the hybrid BAMs no longer show the true reference in IGV.
- *Compare read-side intervals instead* (count only the part of the read not
  covered by its vector alignment). Needs no reference-vs-reference step, but
  ties the hybrid call to the vector alignment's clipping, and a read part
  that is shared sequence would still be called a hybrid.
- *Threshold 50 bp / 200 bp.* 50 bp is close to how far an alignment can
  extend past a homology boundary by chance; 200 bp misses short non-vector
  segments. 100 bp chosen as the middle ground.
- *`asm5` / `map-ont` presets.* `asm5` misses near-identical variant copies of
  an element; `map-ont` is tuned for noisy reads, not reference sequence.

**Known limitations.**
- A real hybrid whose non-vector part lies entirely inside a shared region
  cannot be told apart from pure vector. It is counted as vector: HQ if the
  vector alignments cover ≥ 90% of the read, otherwise pure LQ.
- `asm20` does not splice, so a cDNA transgene is matched to the host genome
  exon by exon; exons shorter than ~50–100 bp may be missed.
- Homology between two non-vector references (e.g. ori/AmpR in backbone and
  helper) is not handled; it only affects which non-vector reference wins.

Every read rejected by the filter is listed in `hybrid_shared_only.tsv`, and
per-reference counts are written to `hybrid_filter_stats.json` and the
classification stats.

---

## (back-filled) HQ coverage ignores reference deletions

**Decision.** HQ/LQ coverage is read bases in M/I/=/X operations divided by
the full read length (including hard clips). Deletions and skips (D/N) do not
lower coverage, so a read with a large internal deletion bridged by one
alignment is HQ.

**Why.** The HQ metric measures how much of the *read* the primary alignment
explains (read-side clipping), not whether the genome is intact; deletions
are characterised separately.

*Back-filled from the code comments in `alignment.calculate_alignment_coverage`;
the original date and alternatives were not recorded.*
