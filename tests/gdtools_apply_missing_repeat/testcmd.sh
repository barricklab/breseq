#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

# Failure test: APPLY must exit non-zero when a MOB references a repeat family that does not
# exist in the reference sequences (a reference feature that isn't there).
EXPECTED_EXIT_CODE=1

TESTCMD="\
    ${GDTOOLS} \
        APPLY \
        -f FASTA \
        -o ${SELF}/output.fasta \
        -r ${DATADIR}/REL606/REL606.fragment.gbk \
        -r ${DATADIR}/REL606/REL606.fragment.2.gbk \
        ${SELF}/input.gd \
    "

do_test $1 ${SELF}
