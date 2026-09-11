#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

# A clean --dry-run must exit 0 AND must leave the filesystem untouched.
#
# The trailing [ ! -e ] is the real regression guard here: Settings::log() calls
# create_path(output_path), so before --dry-run existed, merely constructing a
# Settings created the output directory. If that logging is ever un-gated, this
# test fails -- whereas the exit code on its own would still look fine.
TESTCMD="\
    ${BRESEQ} \
        --dry-run \
        -o ${SELF}/output \
        -r ${DATADIR}/lambda/lambda.gbk \
        ${DATADIR}/lambda/lambda_mixed_population.fastq.gz \
    && [ ! -e ${SELF}/output ] \
    "

do_test $1 ${SELF}
