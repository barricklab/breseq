#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

# Failure test: SUBTRACT must exit non-zero when an input file does not exist.
# The first operand path is intentionally missing; other.gd is a valid second operand.
EXPECTED_EXIT_CODE=1

TESTCMD="\
    ${GDTOOLS} \
        SUBTRACT \
        -o ${SELF}/output.gd \
        ${SELF}/does_not_exist.gd \
        ${SELF}/other.gd \
    "

do_test $1 ${SELF}
