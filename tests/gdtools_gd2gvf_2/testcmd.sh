#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

# The date and the breseq version are in the GVF header
DIFF_IGNORE='^##(file-date|source-method)'

CURRENT_OUTPUTS[0]="${SELF}/output.gvf"
EXPECTED_OUTPUTS[0]="${SELF}/expected.gvf"

# Annotated (adds Variant_effect and the codon attributes), and with every sequence in full
TESTCMD="\
    ${GDTOOLS} \
        CONVERT \
        -f GVF \
        --annotate \
        --gvf-max-sequence-length 0 \
        -o ${SELF}/output.gvf \
        -r ${DATADIR}/REL606/REL606.fragment.gbk \
        ${SELF}/input.gd \
    "

do_test $1 ${SELF}
