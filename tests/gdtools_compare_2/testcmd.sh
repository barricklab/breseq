#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/output.html"
EXPECTED_OUTPUTS[0]="${SELF}/expected.html"

TESTCMD="\
    ${GDTOOLS} \
        COMPARE \
        -f HTML \
        --repeat-header 10 \
        -o ${SELF}/output.html \
        -r ${DATADIR}/pDCAF3/pDCAF3.gbk \
        ${SELF}/input.gd \
    "

do_test $1 ${SELF}
