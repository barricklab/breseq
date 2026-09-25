#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
TEST_CORES=4
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/data/annotated.gd"
EXPECTED_OUTPUTS[0]="${SELF}/expected.gd"

# The same gene conversion as gene_conversion_1 seen through PAIRED 100 bp reads.
#
# A read inside the converted tract whose mate lies in the recipient's unique flank is pulled to
# the recipient by pair rescue, so it is placed uniquely there and the columns nearest the tract
# edges get called: the SNP-driven step sees a cluster of donor-base SNPs at each edge and makes a
# CON out of them, but only as far as the calls reach. The reads then extend that tract across the
# interior, where every read went to the donor, and the CON must come out the same as in
# gene_conversion_1, now resting on the edge RA evidence, the MC and any discordant pair that joins
# the flank to the donor.
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
        -m paired-end \
        -r ${SELF}/output.mutated.fna \
        -o ${SELF}/output.simulated.fastq \
        -l 100 \
        -c 40 \
        --mean 300 \
        --stdev 30 \
        --seed 1 \
    && ${BRESEQ} \
        ${BRESEQ_TEST_THREAD_ARG} \
        -o ${SELF} \
        -r ${SELF}/output.reference.fna \
        ${SELF}/output.simulated_1.fastq \
        ${SELF}/output.simulated_2.fastq \
    && test \$(cut -f1 ${SELF}/data/annotated.gd | grep -c '^CON\$') -eq 1 \
	"

do_test $1 ${SELF}
