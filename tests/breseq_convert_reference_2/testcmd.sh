#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/output.gff3"
EXPECTED_OUTPUTS[0]="${SELF}/expected.gff3"

TESTCMD="\
    ${BRESEQ} \
        CONVERT-REFERENCE \
        --genbank-field-for-seq-id ACCESSION \
        -f GFF \
        -o ${SELF}/output.gff3 \
        ${DATADIR}/pDCAF3/KC619530.gbk \
    "

do_test $1 ${SELF}
