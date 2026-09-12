#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

# Failure test: --dry-run must exit non-zero when a READ file does not exist.
#
# The read files are positional arguments rather than named options, so this is the
# test for AnyOption::setPositionalArgumentsRole() specifically -- every named option
# could be checked correctly and this case still be missed.
EXPECTED_EXIT_CODE=1

TESTCMD="\
    ${BRESEQ} \
        --dry-run \
        -o ${SELF}/output \
        -r ${DATADIR}/lambda/lambda.gbk \
        ${SELF}/does_not_exist.fastq \
    "

do_test $1 ${SELF}
