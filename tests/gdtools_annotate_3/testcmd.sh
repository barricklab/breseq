#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/output.json"
EXPECTED_OUTPUTS[0]="${SELF}/expected.json"

TESTCMD="\
    ${GDTOOLS} \
        COMPARE \
        -f JSON \
        --preserve-evidence \
        -o ${SELF}/output.json \
        -r ${DATADIR}/lambda/lambda.1-2.gbk \
        -r ${DATADIR}/lambda/lambda.3.gbk \
        -r ${DATADIR}/lambda/lambda.4.gbk \
        -r ${DATADIR}/lambda/lambda.5.gbk \
        ${SELF}/input.gd \
    "

do_test $1 ${SELF}
