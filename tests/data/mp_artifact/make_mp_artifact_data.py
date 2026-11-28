#!/usr/bin/env python3
"""Generate the reference and read pairs for tests/mp_sequencing_artifact_mates.

Deterministic (seeded); rerun to regenerate the committed files:
    python3 tests/data/mp_artifact/make_mp_artifact_data.py tests/data/mp_artifact

Reference (mp_artifact.fna, 20 kb of random sequence). The reads come from a SAMPLE genome that is
the reference plus a novel 1,500-bp insertion of random sequence after reference position 14000, so
the mates that fall inside the insertion align nowhere and make a real MP.

Reads: 2x150 FR pairs, fragment length 500, q30-q40 with 0.5% q15 errors, at about 40x of the sample
genome. Then ARTIFACT_PAIRS extra fragments starting just left of reference position 6000 whose mate 2
has failed into G-rich read-through -- the base it starts on and then 70-100% G with scattered other
bases, at mixed quality, as two-colour chemistry produces after a long homopolymer -- so it aligns
nowhere. Without the artifact gate these make an MP at 6000 with nothing inserted there.
"""
import gzip
import os
import random
import sys

SEED = 20261128
LENGTH = 20000
INSERT_AFTER = 14000
INSERT_LENGTH = 1500
ARTIFACT_SITE = 6000
READ, FRAGMENT = 150, 500
COVERAGE = 40
ARTIFACT_PAIRS = 60
ERROR_RATE = 0.005

COMP = str.maketrans("ACGT", "TGCA")


def rc(s):
    return s.translate(COMP)[::-1]


def noisy(rng, read):
    bases, quals = list(read), []
    for i in range(len(bases)):
        if rng.random() < ERROR_RATE:
            bases[i] = rng.choice([b for b in "ACGT" if b != bases[i]])
            quals.append(chr(33 + 15))
        else:
            quals.append(chr(33 + rng.randint(30, 40)))
    return "".join(bases), "".join(quals)


def g_rich_garbage(rng):
    """A failed read: G-rich throughout (70-100% G), some calls at high quality."""
    g_fraction = rng.uniform(0.70, 1.0)
    bases = "".join("G" if rng.random() < g_fraction else rng.choice("ACT") for _ in range(READ))
    quals = "".join(chr(33 + rng.randint(15, 37)) for _ in range(READ))
    return bases, quals


def main(out_dir):
    rng = random.Random(SEED)
    ref = "".join(rng.choice("ACGT") for _ in range(LENGTH))
    insertion = "".join(rng.choice("ACGT") for _ in range(INSERT_LENGTH))
    sample = ref[:INSERT_AFTER] + insertion + ref[INSERT_AFTER:]

    with open(os.path.join(out_dir, "mp_artifact.fna"), "w") as f:
        f.write(">mp-artifact\n")
        for i in range(0, LENGTH, 60):
            f.write(ref[i:i + 60] + "\n")

    pairs = []
    n_tiled = len(sample) * COVERAGE // (2 * READ)
    for i in range(n_tiled):
        start = rng.randint(1, len(sample) - FRAGMENT + 1)
        frag = sample[start - 1:start - 1 + FRAGMENT]
        s1, q1 = noisy(rng, frag[:READ])
        s2f, q2f = noisy(rng, frag[-READ:])
        pairs.append((s1, q1, rc(s2f), q2f[::-1]))

    for i in range(ARTIFACT_PAIRS):
        # Mate 1 forward, ending just before the site; mate 2 is the failed read of the other end.
        start = rng.randint(ARTIFACT_SITE - READ - 60, ARTIFACT_SITE - READ)
        s1, q1 = noisy(rng, ref[start - 1:start - 1 + READ])
        s2, q2 = g_rich_garbage(rng)
        pairs.append((s1, q1, s2, q2))

    rng.shuffle(pairs)
    with gzip.open(os.path.join(out_dir, "mp_artifact_R1.fastq.gz"), "wt", compresslevel=9) as r1, \
         gzip.open(os.path.join(out_dir, "mp_artifact_R2.fastq.gz"), "wt", compresslevel=9) as r2:
        for i, (s1, q1, s2, q2) in enumerate(pairs):
            r1.write("@p%d/1\n%s\n+\n%s\n" % (i, s1, q1))
            r2.write("@p%d/2\n%s\n+\n%s\n" % (i, s2, q2))


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__)))
