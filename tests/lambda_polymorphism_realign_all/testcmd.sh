#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
TEST_CORES=4
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/data/annotated.gd"
EXPECTED_OUTPUTS[0]="${SELF}/expected.gd"

# The lambda_polymorphism_linkage mixture (see its testcmd.sh for what each simulated mutation is
# for) run with --realign all, so that EVERY RA column is re-scored against its candidate
# haplotypes, lone SNPs included, not only the linked runs and polymorphic indels. The same
# check.sh must still hold: the re-scoring may move a frequency a little but must not lose a call,
# split the SUB, or merge the two insertion haplotypes at 10000. The .gd inputs and check.sh are
# shared with that test rather than copied.
SIBLING=${SELF}/../lambda_polymorphism_linkage
TESTCMD=" \
    ${GDTOOLS} APPLY \
        -f FASTA \
        -s NC_001416 \
        -r ${DATADIR}/lambda/lambda.gbk \
        -o ${SELF}/output.mutA.fna \
        ${SIBLING}/mutA.gd \
    && ${GDTOOLS} APPLY \
        -f FASTA \
        -s NC_001416 \
        -r ${DATADIR}/lambda/lambda.gbk \
        -o ${SELF}/output.mutB.fna \
        ${SIBLING}/mutB.gd \
    && ${BRESEQ} SIMULATE-READS \
        -r ${DATADIR}/lambda/lambda.gbk \
        -o ${SELF}/output.ref.fastq \
        -l 100 \
        -c 25 \
        --seed 1 \
    && ${BRESEQ} SIMULATE-READS \
        -r ${SELF}/output.mutA.fna \
        -o ${SELF}/output.mutA.fastq \
        -l 100 \
        -c 15 \
        --seed 2 \
    && ${BRESEQ} SIMULATE-READS \
        -r ${SELF}/output.mutB.fna \
        -o ${SELF}/output.mutB.fastq \
        -l 100 \
        -c 10 \
        --seed 3 \
    && ${BRESEQ} \
        ${BRESEQ_TEST_THREAD_ARG} \
        --polymorphism-prediction \
        --realign all \
        -o ${SELF} \
        -r ${DATADIR}/lambda/lambda.gbk \
        ${SELF}/output.ref.fastq \
        ${SELF}/output.mutA.fastq \
        ${SELF}/output.mutB.fastq \
    && bash ${SIBLING}/check.sh ${SELF}/data/annotated.gd \
	"

do_test $1 ${SELF}
