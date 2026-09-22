#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

# Every type of mutation crossing the origin of a circular sequence: what gdtools APPLY makes of
# it (GFF3 has the features and the sequence), and the VCF and GVF that describe the same change.
# The IS element for the MOBs comes from the second reference file.

# The date and the breseq version are in the VCF and GVF headers
DIFF_IGNORE='^##(fileDate=|source=|file-date|source-method)'

CASES="sub del inv amp con mob mob_del"
REFERENCE_ARG="-r ${DATADIR}/tmv_plasmid/tmv-plasmid.gbk -r ${DATADIR}/REL606/REL606.fragment.gbk"

TESTCMD="true"
i=0
for CASE in ${CASES}; do
  for FORMAT in gff3 vcf gvf; do
    CURRENT_OUTPUTS[$i]="${SELF}/output.${CASE}.${FORMAT}"
    EXPECTED_OUTPUTS[$i]="${SELF}/expected.${CASE}.${FORMAT}"
    i=$((i+1))
  done
  TESTCMD="${TESTCMD} \
    && ${GDTOOLS} APPLY -f GFF3 -s TMV-plasmid -o ${SELF}/output.${CASE}.gff3 ${REFERENCE_ARG} ${SELF}/${CASE}.gd \
    && ${GDTOOLS} CONVERT -f VCF -o ${SELF}/output.${CASE}.vcf ${REFERENCE_ARG} ${SELF}/${CASE}.gd \
    && ${GDTOOLS} CONVERT -f GVF -o ${SELF}/output.${CASE}.gvf ${REFERENCE_ARG} ${SELF}/${CASE}.gd \
    "
done

do_test $1 ${SELF}
