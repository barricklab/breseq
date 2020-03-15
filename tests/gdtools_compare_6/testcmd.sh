#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/output.fa"
EXPECTED_OUTPUTS[0]="${SELF}/expected.fa"

TESTCMD="\
    ${GDTOOLS} \
        COMPARE \
        -f FASTA \
        -o ${SELF}/output.fa \
        -r ${DATADIR}/lambda/lambda.1-2.gbk \
        -r ${DATADIR}/lambda/lambda.3.gbk \
        -r ${DATADIR}/lambda/lambda.4.gbk \
        -r ${DATADIR}/lambda/lambda.5.gbk \
        ${SELF}/input_1.gd \
        ${SELF}/input_2.gd \
        ${SELF}/input_3.gd \
    "

do_test $1 ${SELF}
