#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="output.gd"
EXPECTED_OUTPUTS[0]="expected.gd"

TESTCMD="\
    ${GDTOOLS} \
        CHECK
        -o ${SELF}/output.gd \
        -r ${DATADIR}/lambda/lambda.1-2.gbk \
        -r ${DATADIR}/lambda/lambda.3.gbk \
        -r ${DATADIR}/lambda/lambda.4.gbk \
        -r ${DATADIR}/lambda/lambda.5.gbk \
        ${SELF}/control.gd \
        ${SELF}/test.gd \
    "

do_test $1 ${SELF}
