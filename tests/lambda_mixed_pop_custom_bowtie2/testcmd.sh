#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
TEST_CORES=4
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/data/annotated.gd"
EXPECTED_OUTPUTS[0]="${SELF}/expected.gd"
REFERENCE_ARG="-r ${DATADIR}/lambda/lambda.gbk"

TESTCMD="\
    ${BRESEQ} \
    ${BRESEQ_TEST_THREAD_ARG} \
    -o ${SELF} \
    --bowtie2-stage1 '--ma 1 --mp 2 --np 1 --rdg 4,1 --rfg 4,1 --local --score-min L,10,0.4 -k 20' \
    --bowtie2-stage2 '' \
    --bowtie2-junction '--ma 1 --mp 2 --np 1 --rdg 4,1 --rfg 4,1 --local --score-min L,10,0.30  -k 20' \
    ${REFERENCE_ARG} \
    ${DATADIR}/lambda/lambda_mixed_population.fastq.gz \
    "

do_test $1 ${SELF}
