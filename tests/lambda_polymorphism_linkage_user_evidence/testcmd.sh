#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
TEST_CORES=4
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/data/annotated.gd"
EXPECTED_OUTPUTS[0]="${SELF}/expected.gd"

# The lambda_polymorphism_linkage mixture (see that test: 50% lambda, 30% mutA, 20% mutB, with
# INS AC at 10000 in mutA and INS A at 10000 in mutB), run with --user-evidence-gd naming the A
# insertion at 10000. The user-named allele is reported on its own at the column's frequency AND the
# linked haplotype INS AC is reported at its own; the haplotype INS A is not emitted twice.
# check.sh states what the golden must mean.
TESTCMD=" \
    ${GDTOOLS} APPLY \
        -f FASTA \
        -s NC_001416 \
        -r ${DATADIR}/lambda/lambda.gbk \
        -o ${SELF}/output.mutA.fna \
        ${SELF}/mutA.gd \
    && ${GDTOOLS} APPLY \
        -f FASTA \
        -s NC_001416 \
        -r ${DATADIR}/lambda/lambda.gbk \
        -o ${SELF}/output.mutB.fna \
        ${SELF}/mutB.gd \
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
        --user-evidence-gd ${SELF}/user_evidence.gd \
        -o ${SELF} \
        -r ${DATADIR}/lambda/lambda.gbk \
        ${SELF}/output.ref.fastq \
        ${SELF}/output.mutA.fastq \
        ${SELF}/output.mutB.fastq \
    && bash ${SELF}/check.sh ${SELF}/data/annotated.gd \
	"

do_test $1 ${SELF}
