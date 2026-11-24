#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
TEST_CORES=4
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/data/annotated.gd"
EXPECTED_OUTPUTS[0]="${SELF}/expected.gd"
REFERENCE_ARG="-r ${DATADIR}/lambda/lambda.gbk"

# lambda_polymorphism_aligned_sam with every base quality in the SAM's QUAL column rewritten to '!'
# (Q0), as for PacBio CLR reads aligned with minimap2 after downloading them from the SRA. With
# --aligned-sam the scores arrive in the SAM records rather than through FASTQ conversion: breseq
# must find that the file has one score during sequence conversion, warn once, and replace the
# scores as the records are loaded. The mutations predicted, and their frequencies to within a few
# percent, are the same as lambda_polymorphism_aligned_sam's from the original scores.
#
# The rewritten SAM is written under output.reads/, which "clean" removes with the rest.
TESTCMD="\
    mkdir -p ${SELF}/output.reads \
    && awk 'BEGIN { FS = OFS = \"\\t\" } /^@/ { print; next } { if (\$11 != \"*\") gsub(/./, \"!\", \$11); print }' \
       ${DATADIR}/lambda/lambda_mixed_population.sam > ${SELF}/output.reads/lambda_q0.sam \
    && ${BRESEQ} \
    ${BRESEQ_TEST_THREAD_ARG} \
    -p \
    --aligned-sam \
    -o ${SELF} \
    ${REFERENCE_ARG} \
    ${SELF}/output.reads/lambda_q0.sam \
    2> ${SELF}/output.stderr.txt \
    && test \$(grep -c 'has the same base quality score (Q0)' ${SELF}/output.stderr.txt) -eq 1 \
    && rm -f ${SELF}/output.reads/lambda_q0.sam \
    "

do_test $1 ${SELF}
