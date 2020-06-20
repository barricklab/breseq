#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/output.gff3"
EXPECTED_OUTPUTS[0]="${SELF}/expected.gff3"

TESTCMD="\
    ${BRESEQ} \
        CONVERT-REFERENCE \
        -f GFF \
        -o ${SELF}/output.gff3 \
        ${DATADIR}/pDCAF3/KC619530_variant_LOCUS_no_ACCESSION_no_VERSION.gbk \
    "

do_test $1 ${SELF}
