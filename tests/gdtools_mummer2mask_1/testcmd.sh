#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="output.gd"
EXPECTED_OUTPUTS[0]="expected.gd"

TESTCMD="\
    ${GDTOOLS} \
        MUMMER2MASK \
        -p 10 \
        -o ${SELF}/output.gd \
        -r ${DATADIR}/REL606/REL606.fragment.fna \
        ${SELF}/input.coords \
    "

do_test $1 ${SELF}
