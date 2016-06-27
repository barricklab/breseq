#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="output.html"
EXPECTED_OUTPUTS[0]="expected.html"

TESTCMD="\
    ${GDTOOLS} \
        COMPARE
        -f HTML
        -o ${SELF}/output.html \
        -r ${DATADIR}/pDCAF3/pDCAF3.gbk \
        ${SELF}/input.gd \
    "

do_test $1 ${SELF}
