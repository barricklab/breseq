#!/bin/bash

# CONVERT-REFERENCE --list-seq-ids prints each sequence's name in the file beside the seq id
# breseq uses for it, which is the rule safe_seq_id_name applies stated as data: a tool that
# reads breseq's output can be checked against it rather than against a reading of the C++.

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/output.tsv"
EXPECTED_OUTPUTS[0]="${SELF}/expected.tsv"

TESTCMD="\
    ${BRESEQ} \
        CONVERT-REFERENCE \
        --list-seq-ids \
        -o ${SELF}/output.tsv \
        ${DATADIR}/unsafe_seq_ids/unsafe_seq_ids.fasta \
    "

do_test $1 ${SELF}
