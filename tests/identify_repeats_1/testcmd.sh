#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/output.csv"
EXPECTED_OUTPUTS[0]="${SELF}/expected.csv"

TESTCMD="\
    ${BRESEQ} \
        IDENTIFY-REPEATS \
        -r ${DATADIR}/REL606/REL606.fragment.fna \
        -o ${SELF}/output.csv \
        -l 15 \
        -i 90 \
    "

do_test $1 ${SELF}
