#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
TEST_CORES=4
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/data/annotated.gd"
EXPECTED_OUTPUTS[0]="${SELF}/expected.gd"

# A junction made by the sequencer, not the sample: reads that skipped a hairpin.
#
# Real data, because the signature is in the base qualities and simulated qualities would not be
# evidence of anything. The reference is 3 kb of plasmid CP049750 of Rhodococcoides fascians A21d2
# (GCA_000760905.2, positions 15801-18850) and the reads are every read of Illumina run SRR1259109
# (2x101) that breseq mapped there in a whole-genome run, as single-end reads. Position 1505-1549
# of the fragment (17305-17349 of the plasmid) is the loop side of a near-perfect inverted repeat,
#
#     ...GGTGGCCGTGTCGGGGATGGACTCAGCATCCCCGACACGGCCACC...
#
# and about a third of the reads sequenced in one direction jump it. They pick up again on the far
# side, so they still match the reference there and split into a clean 45-bp deletion junction,
# JC 1504 (-1) / 1550 (+1) -- but the base caller has lost phase and calls the bases after the jump
# at low quality (median about Q6, against Q37 before it). 31 of the 32 reads breseq counts for the
# junction are on one strand.
#
# Expected: the JC is reported with reject=LOW_QUALITY_SIDE, low_quality_side_read_count=21 (of 32)
# and low_quality_side_expected_read_count=3.2.
#
# A read is counted when the mean quality of its bases on one side of the breakpoint is under
# --junction-side-quality-cutoff (15) and at least 10 below the other side, over at least a quarter
# of the read (up to 20 bases). Read by read that is indistinguishable from an ordinary read whose
# quality collapsed part-way through, which this library produces too -- so the junction is judged
# against its library. breseq applies the same test to ordinary reads of the read file at every
# cut position, and the expected count is what this junction's 32 reads, cut where they are cut,
# would give by that rate: 3.2. Seeing 21 is rejected by a one-sided binomial test at
# --junction-low-quality-side-p-value-cutoff (0.001). Reads like this also do not seed candidate
# junctions, nor does a read whose far piece has no base above --base-quality-cutoff at all.
#
# --junction-score-cutoff 0 turns the coverage-evenness test off. On a 3 kb reference that test
# rejects this junction on its own, which would hide whether the quality test did anything; in
# the whole-genome run it does not (score 2.8 against a cutoff of 3), which is how this junction
# came to be reported at 29% frequency.
TESTCMD=" \
    ${BRESEQ} \
    ${BRESEQ_TEST_THREAD_ARG} \
    --junction-score-cutoff 0 \
    -o ${SELF} \
    -r ${DATADIR}/jc_low_quality_side/CP049750_15801-18850.fna \
    ${DATADIR}/jc_low_quality_side/SRR1259109_CP049750_15801-18850.fastq.gz \
    && test \$(grep '^JC' ${SELF}/data/annotated.gd | grep -c 'reject=LOW_QUALITY_SIDE') -eq 1 \
	"

do_test $1 ${SELF}
