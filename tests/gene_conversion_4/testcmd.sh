#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
TEST_CORES=4
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/data/annotated.gd"
EXPECTED_OUTPUTS[0]="${SELF}/expected.gd"

# NEGATIVE control for gene conversion prediction: a whole paralog copy is deleted outright.
#
# Same two-copy reference as gene_conversion_1. input.gd deletes the entire first copy, 4001-5200,
# so the recipient loses its coverage just as a conversion would, and every read from the surviving
# copy still says "donor allele" at every column where the copies differ. The one thing that
# distinguishes this from a conversion is depth: the donor column carries its normal single-copy
# depth, not double. The gene conversion step must leave this alone, and the deletion must be
# predicted as the DEL it is.
TESTCMD=" \
    ${GDTOOLS} APPLY \
        -f FASTA \
        -s NC_001416 \
        -r ${DATADIR}/lambda/lambda.gbk \
        -o ${SELF}/output.reference.fna \
        ${SELF}/reference.gd \
    && ${GDTOOLS} APPLY \
        -f FASTA \
        -s NC_001416 \
        -r ${SELF}/output.reference.fna \
        -o ${SELF}/output.mutated.fna \
        ${SELF}/input.gd \
    && ${BRESEQ} SIMULATE-READS \
        -r ${SELF}/output.mutated.fna \
        -o ${SELF}/output.simulated.fastq \
        -l 36 \
        -c 40 \
        --seed 1 \
    && ${BRESEQ} \
        ${BRESEQ_TEST_THREAD_ARG} \
        -o ${SELF} \
        -r ${SELF}/output.reference.fna \
        ${SELF}/output.simulated.fastq \
    && test \$(cut -f1 ${SELF}/data/annotated.gd | grep -c '^CON\$') -eq 0 \
	"

do_test $1 ${SELF}
