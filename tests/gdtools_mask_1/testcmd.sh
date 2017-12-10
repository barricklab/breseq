#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="output.gd"
EXPECTED_OUTPUTS[0]="expected.gd"

TESTCMD="\
    ${GDTOOLS} \
        MASK \
        -o ${SELF}/output.gd \
        ${SELF}/input.gd \
        ${SELF}/mask.gd \
    "

do_test $1 ${SELF}
