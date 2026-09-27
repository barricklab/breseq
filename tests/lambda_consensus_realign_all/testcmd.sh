#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
TEST_CORES=4
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/data/annotated.gd"
EXPECTED_OUTPUTS[0]="${SELF}/expected.gd"

# Read linkage (LN evidence) and local realignment in CONSENSUS mode, which only --realign all
# turns on there. A clone carrying the lambda_polymorphism_linkage mutA.gd genotype (see that
# testcmd.sh for what each mutation exercises: a 2-base insertion, one more unit of a GCA repeat,
# a T>CA substitution, two SNPs 11 bp apart, a lone SNP, one unit of an ACG repeat deleted, one A
# of a homopolymer deleted) is simulated at 50x with 100-bp reads and called as a clone.
#
# check.sh states what the golden must mean: every simulated mutation is called exactly once, as a
# fixed mutation with no frequency, the multi-column ones licensed by a contiguous LN whose
# refined fit calls the all-variant haplotype consensus, and every RA column carries realigned=1.
TESTCMD=" \
    ${GDTOOLS} APPLY \
        -f FASTA \
        -s NC_001416 \
        -r ${DATADIR}/lambda/lambda.gbk \
        -o ${SELF}/output.mutA.fna \
        ${SELF}/../lambda_polymorphism_linkage/mutA.gd \
    && ${BRESEQ} SIMULATE-READS \
        -r ${SELF}/output.mutA.fna \
        -o ${SELF}/output.mutA.fastq \
        -l 100 \
        -c 50 \
        --seed 4 \
    && ${BRESEQ} \
        ${BRESEQ_TEST_THREAD_ARG} \
        --realign all \
        -o ${SELF} \
        -r ${DATADIR}/lambda/lambda.gbk \
        ${SELF}/output.mutA.fastq \
    && bash ${SELF}/check.sh ${SELF}/data/annotated.gd \
	"

do_test $1 ${SELF}
