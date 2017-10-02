#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="output.json"
EXPECTED_OUTPUTS[0]="expected.json"

TESTCMD="\
    ${GDTOOLS} \
        CONVERT
        -a
        -f JSON
        -o ${SELF}/output.json \
        -r ${DATADIR}/bull/bull_2.gbk \
        ${SELF}/input.gd \
    "

do_test $1 ${SELF}
