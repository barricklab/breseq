#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
TEST_CORES=4
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/data/annotated.gd"
EXPECTED_OUTPUTS[0]="${SELF}/expected.gd"

# Regression test for the starvation guard.
#
# --no-read-alignment-prediction skips the stage-08 pileup, which is also what writes the discordant-
# pair/missing-pair/pair-distance candidate-region CSVs and the per-position coverage table CN reads.
# The stages that CONSUME those files sit outside the guard that skips the pileup, so once CN/DP/MP/PD
# became default-on this combination ran them with their inputs absent and died on a missing file.
# Settings now turns all four off when RA prediction is off.
#
# Nothing here asserts the evidence is absent -- the golden covers that -- but the run completing at
# all is the assertion, since without the guard it aborts.
TESTCMD="\
    ${BRESEQ} \
    ${BRESEQ_TEST_THREAD_ARG} \
    -o ${SELF} \
    --no-read-alignment-prediction \
    -r ${DATADIR}/lambda/lambda.gbk \
    ${DATADIR}/lambda/lambda_mixed_population.fastq.gz \
    "

do_test $1 ${SELF}
