#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/output.gd"
EXPECTED_OUTPUTS[0]="${SELF}/expected.gd"

TESTCMD="\
    ${GDTOOLS} \
        REMOVE \
        --clear-ids \
        -o ${SELF}/output.gd \
        ${SELF}/input.gd \
    "

do_test $1 ${SELF}
