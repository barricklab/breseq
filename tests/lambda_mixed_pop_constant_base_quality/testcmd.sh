#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
TEST_CORES=4
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/data/annotated.gd"
EXPECTED_OUTPUTS[0]="${SELF}/expected.gd"
REFERENCE_ARG="-r ${DATADIR}/lambda/lambda.gbk"

# Reads whose bases all have one quality score, as PacBio CLR reads downloaded from the SRA do
# (every score '!', Q0): the lambda_mixed_pop reads with every quality character rewritten to '!'.
# breseq must warn about the file, once, and give every base the replacement score (Q20 by default),
# whose error rate is then calibrated from the reads. The mutations predicted are the same as
# lambda_mixed_pop's from the original scores. Before this was handled, every base was below
# --base-quality-cutoff, and drawing the error rate plot stopped breseq with a gnuplot error.
#
# The rewritten reads are written under output.reads/, which "clean" removes with the rest.
TESTCMD="\
    mkdir -p ${SELF}/output.reads \
    && gunzip -c ${DATADIR}/lambda/lambda_mixed_population.fastq.gz \
       | awk 'NR % 4 == 0 { gsub(/./, \"!\") } { print }' > ${SELF}/output.reads/lambda_q0.fastq \
    && ${BRESEQ} \
    ${BRESEQ_TEST_THREAD_ARG} \
    -o ${SELF} \
    ${REFERENCE_ARG} \
    ${SELF}/output.reads/lambda_q0.fastq \
    2> ${SELF}/output.stderr.txt \
    && test \$(grep -c 'has the same base quality score (Q0)' ${SELF}/output.stderr.txt) -eq 1 \
    && rm -f ${SELF}/output.reads/lambda_q0.fastq \
    "

do_test $1 ${SELF}
