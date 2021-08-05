#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/output.gd"
EXPECTED_OUTPUTS[0]="${SELF}/expected.gd"

TESTCMD="\
    ${GDTOOLS} \
        REMOVE \
        -e \
        --comment \
        -c 'frequency<1' \
        -c label==remove_me \
        -o ${SELF}/output.gd \
        ${SELF}/input.gd \
    "

do_test $1 ${SELF}
