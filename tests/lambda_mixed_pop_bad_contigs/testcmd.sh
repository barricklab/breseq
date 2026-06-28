#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
TEST_CORES=4
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/data/annotated.gd"
EXPECTED_OUTPUTS[0]="${SELF}/expected.gd"
REFERENCE_ARG="-r ${DATADIR}/lambda/lambda.gbk -r ${DATADIR}/lambda/other.gbk -r ${DATADIR}/REL606/REL606.fragment.gbk"

TESTCMD="\
    ${BRESEQ} \
    ${BRESEQ_TEST_THREAD_ARG} \
    --output-unmapped-reads \
    -o ${SELF} \
    ${REFERENCE_ARG} \
    ${DATADIR}/lambda/lambda_mixed_population.fastq.gz ; \
    if [ ! -f "${SELF}/data/unmapped_reads.fastq.gz" ]; then \
        echo 'Unmapped read fastq file not created!'; \
        exit 1; \
    fi \
    "

do_test $1 ${SELF}
