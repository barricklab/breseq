#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
TEST_CORES=4
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/data/annotated.gd"
EXPECTED_OUTPUTS[0]="${SELF}/expected.gd"

# The four deprecated opt-in flags must still be accepted, must warn, and must change nothing --
# they are now what the run does anyway. A deprecated option that errors out breaks every
# existing script and pipeline that passes it, which is the whole point of keeping it registered
# (DEPRECATED_OPTION in anyoption.h: parsed normally, never shown in help).
#
# The grep -c asserts all four warnings fire. Without it the test would still pass if a flag were
# silently dropped from the option table, since AnyOption exits on an UNKNOWN option but this
# golden would be identical either way.
TESTCMD="\
    ${BRESEQ} \
    ${BRESEQ_TEST_THREAD_ARG} \
    -o ${SELF} \
    --predict-copy-number \
    --predict-discordant-pairs \
    --predict-missing-pairs \
    --predict-pair-distance \
    -r ${DATADIR}/lambda/lambda.gbk \
    ${DATADIR}/lambda/lambda_mixed_population.fastq.gz \
    2> ${SELF}/output.stderr.txt \
    && test \$(grep -c 'option is DEPRECATED' ${SELF}/output.stderr.txt) -eq 4 \
    "

do_test $1 ${SELF}
