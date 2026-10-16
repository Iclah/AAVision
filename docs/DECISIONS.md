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
