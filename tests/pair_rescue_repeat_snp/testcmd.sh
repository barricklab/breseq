#!/bin/bash
#
# Concordant rescue of read pairs between near-identical repeat copies (rescue_concordant_pair(),
# resolve_alignments.cpp; --pair-rescue-max-log10-likelihood-loss).
#
# tests/data/pair_rescue (make_pair_rescue_data.py) is a 20-kb reference with two 1-kb copies, A at
# 5001-6000 and B at 15001-16000, that differ at one base (5500 / 15500), and 2x150 pairs at about
# 40x with overlapping mates. In 30 of the pairs from copy A both mates cover 5500: mate 1 reads A's
# base there at q40, mate 2 reads B's base at q10. Each mate then has a unique best alignment, on
# different copies, so without rescue every one of those pairs is discordant and together they make a
# DP joining the two copies -- the shape of four spurious DPs on a Thermincola CM166461 Illumina run,
# where nanopore reads showed the reference to be right.
#
# The same reads are run twice:
#   output_baseline/  rescue OFF: the DP between the copies is predicted (the case is live)
#   data/             rescue ON (3 log10 units): no DP; mate 2 costs one q10 mismatch on copy A
# check.sh states both outcomes. The golden is the rescue run's.

SELF=`dirname ${BASH_SOURCE}`
TEST_CORES=2
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/data/annotated.gd"
EXPECTED_OUTPUTS[0]="${SELF}/expected.gd"
REFERENCE_ARG="-r ${DATADIR}/pair_rescue/pair_rescue.fna"
READS="${DATADIR}/pair_rescue/pair_rescue_R1.fastq.gz ${DATADIR}/pair_rescue/pair_rescue_R2.fastq.gz"

TESTCMD=" \
    ${BRESEQ} \
        ${BRESEQ_TEST_THREAD_ARG} \
        -o ${SELF}/output_baseline \
        ${REFERENCE_ARG} \
        --pair-rescue-max-log10-likelihood-loss 0 \
        --pair-ambiguity-max-log10-likelihood-loss 0 \
        ${READS} \
    && ${BRESEQ} \
        ${BRESEQ_TEST_THREAD_ARG} \
        -o ${SELF} \
        ${REFERENCE_ARG} \
        --pair-rescue-max-log10-likelihood-loss 3 \
        ${READS} \
    && bash ${SELF}/check.sh ${SELF}/output_baseline/output/output.gd ${SELF}/data/annotated.gd \
    "

do_test $1 ${SELF}
