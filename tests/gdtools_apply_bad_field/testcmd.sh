#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

# Failure test: a malformed field in an otherwise-valid line must abort a pipeline subcommand.
# input.gd has a garbled (non key=value) trailing field, now a fatal parse error on read().
EXPECTED_EXIT_CODE=1

TESTCMD="\
    ${GDTOOLS} \
        APPLY \
        -f FASTA \
        -o ${SELF}/output.fasta \
        -r ${DATADIR}/REL606/REL606.fragment.2.gbk \
        ${SELF}/input.gd \
    "

do_test $1 ${SELF}
