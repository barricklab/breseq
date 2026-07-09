#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/output.tre"
EXPECTED_OUTPUTS[0]="${SELF}/expected.tre"
CURRENT_OUTPUTS[1]="${SELF}/output.mutation.key.txt"
EXPECTED_OUTPUTS[1]="${SELF}/expected.mutation.key.txt"

TESTCMD="\
    ${GDTOOLS} \
        PHYLOGENY \
        --missing-as-ancestral \
        -o ${SELF}/output \
        -r ${DATADIR}/lambda/lambda.1-2.gbk \
        -r ${DATADIR}/lambda/lambda.3.gbk \
        -r ${DATADIR}/lambda/lambda.4.gbk \
        -r ${DATADIR}/lambda/lambda.5.gbk \
        ${SELF}/input_1.gd \
        ${SELF}/input_2.gd \
        ${SELF}/input_3.gd \
    "

do_test $1 ${SELF}
