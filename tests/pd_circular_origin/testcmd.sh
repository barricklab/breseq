#!/bin/bash
#
# Pair distance (PD) evidence within one fragment length of a circular origin is ignored.
#
# Pair handling is linear: a pair whose fragment runs across the origin of a circular sequence is
# never measured (a mate crossing it aligns to the origin junction and the pair is not classified;
# mates on either side of it are measured about a sequence length apart). Those are the LONG
# fragments covering positions near the origin, so the surviving pairs there look short and PD saw a
# size shift with no event -- PD 5603/5604, -170 bp, on a 6,114-bp Salmonella plasmid at 620x.
# CLAUDE.md ("Pair handling is linear") records what a real fix would take.
#
# Simulated 2x150 pairs, 450-bp fragments, 600x on the unmutated 10,411-bp circular TMV plasmid.
# Before the fix this predicted an unrejected PD 562/563 (size_shift -37) and nothing else. check.sh
# requires every PD within the library's distance cutoff of the origin to be
# ignore=CIRCULAR_CHROMOSOME and no mutation at all. PDs away from an origin are covered by
# lambda_pair_distances, which this change leaves as it was.

SELF=`dirname ${BASH_SOURCE}`
TEST_CORES=4
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/data/annotated.gd"
EXPECTED_OUTPUTS[0]="${SELF}/expected.gd"
REFERENCE_ARG="-r ${DATADIR}/tmv_plasmid/tmv-plasmid.gbk"

TESTCMD=" \
    ${BRESEQ} SIMULATE-READS \
        -m paired-end \
        ${REFERENCE_ARG} \
        -o ${SELF}/output.simulated.fastq \
        -l 150 \
        -c 600 \
        --mean 450 \
        --stdev 60 \
        --seed 7 \
    && ${BRESEQ} \
        ${BRESEQ_TEST_THREAD_ARG} \
        -o ${SELF} \
        ${REFERENCE_ARG} \
        ${SELF}/output.simulated_1.fastq \
        ${SELF}/output.simulated_2.fastq \
    && bash ${SELF}/check.sh ${SELF}/data/annotated.gd 10411 \
    "

do_test $1 ${SELF}
