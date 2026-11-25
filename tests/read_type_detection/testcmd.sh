#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

# Read type detection from read names and base qualities, and parsing of the <read type>:<path>
# prefix on read file arguments -- see run_checks.sh for what it pins down. It builds its read
# files from the lambda reads and never runs the pipeline, so it takes a second or two.

CURRENT_OUTPUTS[0]="${SELF}/output.txt"
EXPECTED_OUTPUTS[0]="${SELF}/expected.txt"

TESTCMD="\
    ${SELF}/run_checks.sh ${SELF} ${BRESEQ} ${DATADIR} \
    "

do_test $1 ${SELF}
