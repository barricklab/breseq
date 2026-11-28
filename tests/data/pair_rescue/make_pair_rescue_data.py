#!/usr/bin/env python3
"""Generate the reference and read pairs for tests/pair_rescue_repeat_snp.

Deterministic (seeded); rerun to regenerate the committed files:
    python3 tests/data/pair_rescue/make_pair_rescue_data.py tests/data/pair_rescue

Reference (pair_rescue.fna, 20 kb): random sequence with two copies of a 1-kb repeat, A at 5001-6000
and B at 15001-16000, identical except at the repeat's 500th base (5500 / 15500).

Reads: 2x150 FR pairs, fragment length 250 (the mates overlap by 50 bases), base qualities q30-q40
with ERROR_RATE of bases replaced by a wrong base at q15 (an error-free library leaves the error-rate
plot with nothing to draw), tiled across the whole reference at about 40x. The error never falls on
the differing base. Then RESCUE_PAIRS extra fragments from copy A whose two
mates BOTH cover 5500: mate 1 reads copy A's base there at q40, mate 2 reads copy B's base at q10.
Each mate then has a unique best alignment, on different copies, so without rescue the pair is
discordant; placing mate 2 on copy A costs one q10 mismatch.
"""
import gzip
import os
import random
import sys

SEED = 20261128
LENGTH = 20000
REPEAT = 1000
A_START, B_START = 5001, 15001
DIFF_OFFSET = 500           # 1-based offset in the repeat of the one differing base
READ, FRAGMENT = 150, 250
COVERAGE = 40
RESCUE_PAIRS = 30
ERROR_RATE = 0.005

COMP = str.maketrans("ACGT", "TGCA")


def rc(s):
    return s.translate(COMP)[::-1]


def noisy(rng, read, protect=-1):
    """q30-q40 qualities, with ERROR_RATE of bases (never index `protect`) changed at q15."""
    bases, quals = list(read), []
    for i in range(len(bases)):
        if (i != protect) and (rng.random() < ERROR_RATE):
            bases[i] = rng.choice([b for b in "ACGT" if b != bases[i]])
            quals.append(chr(33 + 15))
        else:
            quals.append(chr(33 + rng.randint(30, 40)))
    return "".join(bases), "".join(quals)


def main(out_dir):
    rng = random.Random(SEED)
    seq = [rng.choice("ACGT") for _ in range(LENGTH)]
    repeat = [rng.choice("ACGT") for _ in range(REPEAT)]
    seq[A_START - 1:A_START - 1 + REPEAT] = repeat
    seq[B_START - 1:B_START - 1 + REPEAT] = repeat
    a_pos, b_pos = A_START + DIFF_OFFSET - 1, B_START + DIFF_OFFSET - 1
    seq[b_pos - 1] = rng.choice([b for b in "ACGT" if b != seq[a_pos - 1]])
    seq = "".join(seq)

    with open(os.path.join(out_dir, "pair_rescue.fna"), "w") as f:
        f.write(">pair-rescue\n")
        for i in range(0, LENGTH, 60):
            f.write(seq[i:i + 60] + "\n")

    def protect_index(start, read_start):
        """Index of a_pos or b_pos within a read starting at read_start, or -1."""
        for p in (a_pos, b_pos):
            if read_start <= p < read_start + READ:
                return p - read_start
        return -1

    pairs = []
    n_tiled = LENGTH * COVERAGE // (2 * READ)
    for i in range(n_tiled):
        start = rng.randint(1, LENGTH - FRAGMENT + 1)
        frag = seq[start - 1:start - 1 + FRAGMENT]
        s1, q1 = noisy(rng, frag[:READ], protect_index(start, start))
        s2_fwd, q2_fwd = noisy(rng, frag[-READ:], protect_index(start, start + FRAGMENT - READ))
        pairs.append((s1, q1, rc(s2_fwd), q2_fwd[::-1]))

    b_base = seq[b_pos - 1]
    for i in range(RESCUE_PAIRS):
        # Both mates cover a_pos: mate 1 covers [start, start+149], mate 2 [start+100, start+249].
        start = rng.randint(a_pos - READ + 1, a_pos - (FRAGMENT - READ))
        frag = seq[start - 1:start - 1 + FRAGMENT]
        k = a_pos - (start + FRAGMENT - READ)          # index of a_pos within mate 2, forward strand
        s1, q1 = noisy(rng, frag[:READ], a_pos - start)
        s2_fwd, q2_fwd = noisy(rng, frag[-READ:], k)
        s2_fwd = s2_fwd[:k] + b_base + s2_fwd[k + 1:]
        q2_fwd = q2_fwd[:k] + chr(33 + 10) + q2_fwd[k + 1:]   # copy B's base, at q10
        pairs.append((s1, q1, rc(s2_fwd), q2_fwd[::-1]))

    rng.shuffle(pairs)
    with gzip.open(os.path.join(out_dir, "pair_rescue_R1.fastq.gz"), "wt", compresslevel=9) as r1, \
         gzip.open(os.path.join(out_dir, "pair_rescue_R2.fastq.gz"), "wt", compresslevel=9) as r2:
        for i, (s1, q1, s2, q2) in enumerate(pairs):
            r1.write("@p%d/1\n%s\n+\n%s\n" % (i, s1, q1))
            r2.write("@p%d/2\n%s\n+\n%s\n" % (i, s2, q2))


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__)))
