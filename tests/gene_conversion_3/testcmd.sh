#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
TEST_CORES=4
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/data/annotated.gd"
EXPECTED_OUTPUTS[0]="${SELF}/expected.gd"

# The gene conversion of gene_conversion_1 with the donor on the OTHER strand.
#
# reference.gd here writes the diverged copy at 12001-13200 as the reverse complement of the one
# homologous_deletion_1 uses, so the two paralogs are 96% identical but antiparallel. input.gd
# converts the same 400 bp of the first copy, taking its sequence from the second copy read
# backwards (region start > end). Everything about the reads is as in gene_conversion_1; what is
# being checked is that the donor is found on the reverse strand, that the tract is extended along
# it in the right direction, and that the CON's region is written backwards.
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
    && test \$(cut -f1 ${SELF}/data/annotated.gd | grep -c '^CON\$') -eq 1 \
	"

do_test $1 ${SELF}
