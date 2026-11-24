#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
TEST_CORES=4
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/data/annotated.gd"
EXPECTED_OUTPUTS[0]="${SELF}/expected.gd"
REFERENCE_ARG="-r ${DATADIR}/lambda/lambda.gbk"

# As lambda_mixed_pop_constant_base_quality (every base Q0), but with
# --no-constant-base-quality-replacement, so the scores are kept as they are. The error rate table
# then has a single quality score, which gave gnuplot an x range of zero width ("Can't plot with an
# empty x range!") and stopped breseq. breseq must finish and draw the plot (checked in the
# unzipped evidence directory), and must not warn about the scores. Every base is at or below
# --base-quality-cutoff, so no mutations are predicted.
TESTCMD="\
    mkdir -p ${SELF}/output.reads \
    && gunzip -c ${DATADIR}/lambda/lambda_mixed_population.fastq.gz \
       | awk 'NR % 4 == 0 { gsub(/./, \"!\") } { print }' > ${SELF}/output.reads/lambda_q0.fastq \
    && ${BRESEQ} \
    ${BRESEQ_TEST_THREAD_ARG} \
    --no-constant-base-quality-replacement \
    --unzipped-html \
    -o ${SELF} \
    ${REFERENCE_ARG} \
    ${SELF}/output.reads/lambda_q0.fastq \
    2> ${SELF}/output.stderr.txt \
    && test -s ${SELF}/output/evidence/lambda_q0.error_rates.svg \
    && test \$(grep -c 'has the same base quality score' ${SELF}/output.stderr.txt) -eq 0 \
    && rm -f ${SELF}/output.reads/lambda_q0.fastq \
    "

do_test $1 ${SELF}
