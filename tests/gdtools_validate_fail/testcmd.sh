#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

# Failure test: VALIDATE must exit non-zero when a file has a formatting error.
# input.gd contains a garbled (non key=value) trailing field, a fatal parse error.
EXPECTED_EXIT_CODE=1

TESTCMD="\
    ${GDTOOLS} \
        VALIDATE \
        ${SELF}/input.gd \
    "

do_test $1 ${SELF}
