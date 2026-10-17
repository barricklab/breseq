#!/bin/bash

# Two sequences whose names collapse onto one seq id (`__x__` and `[x]` are both `x`) must be
# refused, not merged: the FASTA reader used to let the second record silently replace the
# first, because it passed a string where set_sequence_loaded_from_file's bool went.

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

EXPECTED_EXIT_CODE=1

TESTCMD="\
    ${BRESEQ} \
        CONVERT-REFERENCE \
        --list-seq-ids \
        -o ${SELF}/output.tsv \
        ${DATADIR}/unsafe_seq_ids/colliding_seq_ids.fasta \
    "

do_test $1 ${SELF}
