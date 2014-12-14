#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="output.gd"
EXPECTED_OUTPUTS[0]="expected.gd"

TESTCMD="\
    ${GDTOOLS} \
        INTERSECT
        -o ${SELF}/output.gd \
        ${SELF}/input_1.gd \
        ${SELF}/input_2.gd \
    "

do_test $1 ${SELF}
