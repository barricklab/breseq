#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
TEST_CORES=4
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/data/annotated.gd"
EXPECTED_OUTPUTS[0]="${SELF}/expected.gd"
REFERENCE_ARG="-r ${DATADIR}/lambda/lambda-contig.gbk"

# check_reference_groups.sh pins that the five contigs reached CNery as one reference group. It is
# chained here rather than expressed in expected.gd because it cannot be: these contigs are equal
# slices of one genome at one depth, so grouped and ungrouped runs produce an identical .gd.
TESTCMD="\
    ${BRESEQ} \
    ${BRESEQ_TEST_THREAD_ARG} \
    -o ${SELF} \
    -c ${DATADIR}/lambda/lambda-contig.gbk \
    ${DATADIR}/lambda/lambda_mixed_population.fastq.gz \
    && ${SELF}/check_reference_groups.sh ${SELF} \
    "

do_test $1 ${SELF}
