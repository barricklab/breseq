#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

# The date and the breseq version are in the VCF header
DIFF_IGNORE='^##(fileDate|source)='

CURRENT_OUTPUTS[0]="${SELF}/output.vcf"
EXPECTED_OUTPUTS[0]="${SELF}/expected.vcf"

TESTCMD="\
    ${GDTOOLS} \
        CONVERT \
        -f VCF \
        -o ${SELF}/output.vcf \
        -r ${DATADIR}/REL606/REL606.fragment.gbk \
        ${SELF}/input.gd \
    "

do_test $1 ${SELF}
