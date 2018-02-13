#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/output.csv"
EXPECTED_OUTPUTS[0]="${SELF}/expected.csv"

TESTCMD="\
    ${GDTOOLS} \
        COUNT \
        -b \
        -o ${SELF}/output.csv \
        -r ${DATADIR}/lambda/lambda.1-2.gbk \
        -r ${DATADIR}/lambda/lambda.3.gbk \
        -r ${DATADIR}/lambda/lambda.4.gbk \
        -r ${DATADIR}/lambda/lambda.5.gbk \
        ${SELF}/input.gd \
    "

do_test $1 ${SELF}
