#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
TEST_CORES=4
. ${SELF}/../common.sh

# Re-analyzes the lambda_mixed_pop reads against the reference with some of that test's own
# mutations applied first (--apply-check). Those mutations should no longer be predicted, the rest
# shift by the applied indels, and every entry carries original_* fields whose values are the
# positions lambda_mixed_pop reports. The input also applies a 10-bp INS the reads do not have, so
# the DEL/MC/UN/JC that call it back sit inside inserted sequence and exercise the
# original_*_offset fields (anchor 30000, offsets 1-10) and the "30,000+1" HTML rendering.
CURRENT_OUTPUTS[0]="${SELF}/data/annotated.gd"
EXPECTED_OUTPUTS[0]="${SELF}/expected.gd"
CURRENT_OUTPUTS[1]="${SELF}/data/original_coordinates.tsv"
EXPECTED_OUTPUTS[1]="${SELF}/expected_original_coordinates.tsv"
REFERENCE_ARG="-r ${DATADIR}/lambda/lambda.gbk"
COMPARE_ARG="--genbank-field-for-seq-id VERSION"

TESTCMD="\
    ${BRESEQ} \
    ${BRESEQ_TEST_THREAD_ARG} \
    -o ${SELF} \
    --genbank-field-for-seq-id VERSION \
    --apply-check ${SELF}/input.gd \
    ${REFERENCE_ARG} \
    ${DATADIR}/lambda/lambda_mixed_population.fastq.gz \
    "

do_test $1 ${SELF}
