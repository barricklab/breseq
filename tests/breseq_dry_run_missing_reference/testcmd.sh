#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

# Failure test: --dry-run must exit non-zero when a -r reference file does not exist.
EXPECTED_EXIT_CODE=1

TESTCMD="\
    ${BRESEQ} \
        --dry-run \
        -o ${SELF}/output \
        -r ${SELF}/does_not_exist.gbk \
        ${DATADIR}/lambda/lambda_mixed_population.fastq.gz \
    "

do_test $1 ${SELF}
