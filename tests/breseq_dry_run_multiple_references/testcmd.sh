#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

# A clean --dry-run with -r repeated and with -c and -s also supplied.
#
# Two things are pinned here that a single-reference test cannot see:
#
#  1. Repeated options are stored newline-joined by AnyOption::setValue, so the path
#     checker has to split that value rather than check the literal path "a\nb".
#  2. -c and -s used to be declared with NULL as their default. NULL resolves to the
#     TEMPLATE operator() overload (decltype(NULL) is an integer type, an exact match
#     for const T&), not the void* "takes no argument" one, so they registered a
#     default VALUE of the string "0" -- which `breseq -h` printed as "(DEFAULT=0)".
#     Checking their effective value therefore reported a missing file named "0" on
#     every run that omitted them.
TESTCMD="\
    ${BRESEQ} \
        --dry-run \
        -o ${SELF}/output \
        -r ${DATADIR}/lambda/lambda.1-2.gbk \
        -r ${DATADIR}/lambda/lambda.3.gbk \
        -c ${DATADIR}/lambda/lambda-contig.gbk \
        -s ${DATADIR}/lambda/lambda.4.gbk \
        ${DATADIR}/lambda/lambda_mixed_population.fastq.gz \
    && [ ! -e ${SELF}/output ] \
    "

do_test $1 ${SELF}
