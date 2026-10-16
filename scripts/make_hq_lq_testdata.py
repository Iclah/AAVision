"""
Generate synthetic reads for the HQ/LQ vector-read classification tests.

Writes three files into the output directory:
  vector.fa    - a random 4,700 bp "vector" reference
  reads.fastq  - reads with known structure, in their ORIGINAL orientation
  truth.tsv    - read_id, structure, strand, expected_class, read_length

expected_class is one of:
  HQ        primary + supplementary vector alignments cover >= 90% of the read
            (merged in read coordinates) and the read is not a snapback
  LQ        mapped to the vector, but < 90% covered, or a snapback (two records
            on opposite strands overlapping >= 200 bp on the vector)
  unmapped  no alignment to the vector

All coordinates below are 0-based, half-open Python slices of the vector.
Reads labelled strand '-' are the reverse complement of the '+' construct.

The structures are chosen so that each expected class is well clear of the
0.9 threshold even with the simulated substitution errors, EXCEPT where the
structure is the point of the test (primary92_small_supp, below_threshold).
Internal deletions bridged by a single alignment are expected HQ: the HQ
metric measures read-side clipping only (see docs/DECISIONS.md).

Usage:
    python scripts/make_hq_lq_testdata.py --outdir tests/data/hq_lq [--seed 20261016]
"""

import argparse
import os
import random

VECTOR_LEN = 4700
SUB_RATE = 0.01  # substitution-only errors keep read lengths exact
READS_PER_CLASS = 15

_COMP = str.maketrans("ACGT", "TGCA")


def revcomp(seq: str) -> str:
    """Reverse complement of an uppercase ACGT sequence."""
    return seq[::-1].translate(_COMP)


def _random_seq(rng: random.Random, n: int) -> str:
    return "".join(rng.choice("ACGT") for _ in range(n))


def _mutate(rng: random.Random, seq: str, rate: float) -> str:
    out = []
    for base in seq:
        if rng.random() < rate:
            out.append(rng.choice([b for b in "ACGT" if b != base]))
        else:
            out.append(base)
    return "".join(out)


def _constructs(rng: random.Random, v: str):
    """
    Yield (structure, expected_class, '+'-strand sequence) per read.

    Each class draws a small random jitter so reads are not identical.
    """
    for _ in range(READS_PER_CLASS):
        j = rng.randint(0, 100)

        # Full-length genome
        yield "full_length", "HQ", v

        # Snapback, 50/50: first half of a 4 kb read is vector[j:j+2000],
        # second half folds back on itself. Primary ~0.5, supplementary
        # (hard-clipped without -Y) covers the other half.
        a = v[j : j + 2000]
        yield "snapback_50_50", "LQ", a + revcomp(a)

        # Snapback, 70/30: primary ~0.72, supplementary ~0.28
        yield "snapback_70_30", "LQ", v[: 3300 - j] + revcomp(v[2000 : 3300 - j])

        # Truncated genome + 2 kb non-vector sequence: primary ~0.5, no supplementary
        yield "trunc_plus_junk", "LQ", v[j : j + 2000] + _random_seq(rng, 2000)

        # 1.7 kb internal deletion, bridged by one alignment (D op): expected HQ
        yield "deletion_bridged", "HQ", v[: 1500 + j] + v[3200 + j :]

        # Primary ~0.97 of the read, plus a small 150 bp foldback (overlap
        # < 200 bp, so not a snapback even if it gets its own record): HQ
        yield "primary92_small_supp", "HQ", v[:4400] + revcomp(v[4250:4400])

        # 400 bp foldback: if it aligns it overlaps the primary by 400 bp on the
        # opposite strand (snapback); if not, coverage 4400/4800 < 0.9. LQ either way
        yield "foldback_400", "LQ", v[:4400] + revcomp(v[4000:4400])

        # 2 kb vector + 1 kb junk + 1 kb vector: merged coverage ~0.75 -> LQ.
        # The second vector piece lies upstream of the first on the vector, so
        # minimap2 must split the read instead of bridging the junk as one
        # alignment with a 1 kb insertion (insertions count as covered bases).
        yield (
            "split_vec_junk_vec",
            "LQ",
            v[2500 + j : 4500 + j] + _random_seq(rng, 1000) + v[0:1000],
        )

        # All vector, permuted (second half first): two non-overlapping
        # same-strand records cover the whole read -> HQ
        yield "split_all_vector", "HQ", v[2400:] + v[: 2300 - j]

        # 3.5 kb vector + 500 bp junk: primary ~0.875 < 0.9 -> LQ
        yield "below_threshold", "LQ", v[j : j + 3500] + _random_seq(rng, 500)

        # Pure junk: does not map to the vector
        yield "unmapped_junk", "unmapped", _random_seq(rng, 3000)


def generate(outdir: str, seed: int = 20261016) -> dict:
    """
    Write vector.fa, reads.fastq and truth.tsv into outdir.

    Args:
        outdir: Output directory (created if missing)
        seed: Random seed; the same seed always gives identical files

    Returns:
        Dict with keys 'vector', 'reads', 'truth' mapping to the file paths
    """
    rng = random.Random(seed)
    os.makedirs(outdir, exist_ok=True)

    vector = _random_seq(rng, VECTOR_LEN)
    paths = {
        "vector": os.path.join(outdir, "vector.fa"),
        "reads": os.path.join(outdir, "reads.fastq"),
        "truth": os.path.join(outdir, "truth.tsv"),
    }

    with open(paths["vector"], "w") as fa:
        fa.write(">vector\n" + vector + "\n")

    with open(paths["reads"], "w") as fq, open(paths["truth"], "w") as tt:
        tt.write("read_id\tstructure\tstrand\texpected_class\tread_length\n")
        n = 0
        for structure, expected, seq in _constructs(rng, vector):
            for strand in ("+", "-"):
                read = seq if strand == "+" else revcomp(seq)
                read = _mutate(rng, read, SUB_RATE)
                read_id = f"r{n:04d}_{structure}_{'fwd' if strand == '+' else 'rev'}"
                n += 1
                fq.write(f"@{read_id}\n{read}\n+\n{'I' * len(read)}\n")
                tt.write(f"{read_id}\t{structure}\t{strand}\t{expected}\t{len(read)}\n")

    return paths


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--seed", type=int, default=20261016, help="Random seed")
    args = parser.parse_args()
    paths = generate(args.outdir, args.seed)
    for name, path in paths.items():
        print(f"{name}: {path}")


if __name__ == "__main__":
    main()
