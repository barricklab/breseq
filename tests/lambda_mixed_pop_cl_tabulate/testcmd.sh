#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
# This test must be run AFTER lambda_mixed_pop, which generates its input
# (data/reference.bam, data/reference.fasta, data/reference.gff3). The
# Snakefile parses TEST_DEPENDS to enforce this ordering when running tests
# in parallel; common.sh ignores it (serial runners rely on sorted discovery
# order, where lambda_mixed_pop already comes first).
TEST_DEPENDS=lambda_mixed_pop
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/contingency_loci.csv"
EXPECTED_OUTPUTS[0]="${SELF}/expected.csv"

TESTCMD="\
    ${BRESEQ} CL-TABULATE \
    -b ${SELF}/../lambda_mixed_pop/data/reference.bam \
    -f ${SELF}/../lambda_mixed_pop/data/reference.fasta \
    -r ${SELF}/../lambda_mixed_pop/data/reference.gff3 \
    -m 5 \
    -s \
    -o ${SELF}/contingency_loci.csv \
    "

do_test $1 ${SELF}
