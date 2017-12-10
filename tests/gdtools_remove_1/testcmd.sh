#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="output.gd"
EXPECTED_OUTPUTS[0]="expected.gd"

TESTCMD="\
    ${GDTOOLS} \
        REMOVE \
        -c 'frequency<1' \
        -c label==remove_me \
        -o ${SELF}/output.gd \
        ${SELF}/input.gd \
    "

do_test $1 ${SELF}
