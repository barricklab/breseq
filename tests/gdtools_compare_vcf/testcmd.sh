#!/bin/bash
#
# gdtools COMPARE -f VCF writes one GT:AF column per sample (#238). These three samples (titles all
# "gd", renamed gd/gd_1/gd_2) cover mutations that are fixed, polymorphic (GT 1, AF below 1), absent
# (0:0), and uncallable in a sample ('?' from a UN region -> .:.).
#
SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

# The date and the breseq version are in the VCF header
DIFF_IGNORE='^##(fileDate|source)='

CURRENT_OUTPUTS[0]="${SELF}/output.vcf"
EXPECTED_OUTPUTS[0]="${SELF}/expected.vcf"

TESTCMD="\
    ${GDTOOLS} \
        COMPARE \
        -f VCF \
        -o ${SELF}/output.vcf \
        -r ${DATADIR}/lambda/lambda.1-2.gbk \
        -r ${DATADIR}/lambda/lambda.3.gbk \
        -r ${DATADIR}/lambda/lambda.4.gbk \
        -r ${DATADIR}/lambda/lambda.5.gbk \
        ${SELF}/../gdtools_compare_3/input_1.gd \
        ${SELF}/../gdtools_compare_3/input_2.gd \
        ${SELF}/../gdtools_compare_3/input_3.gd \
    "

do_test $1 ${SELF}
