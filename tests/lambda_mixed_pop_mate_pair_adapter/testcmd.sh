#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
TEST_CORES=4
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/data/annotated.gd"
EXPECTED_OUTPUTS[0]="${SELF}/expected.gd"
REFERENCE_ARG="-r ${DATADIR}/lambda/lambda.gbk"

# Reads as an Illumina Nextera mate-pair library has them (see MATE-PAIR JUNCTION ADAPTERS in
# src/breseq/fastq.h): every fourth lambda_mixed_pop read is rewritten as its own first 22 bases,
# the 38-base junction adapter, and the first 22 bases of the read before it -- genome, adapter,
# then genome from somewhere else, which is what a read through the junction is. Every eighth has
# one base of the adapter changed, as a sequencing error would (mate_pair_reads.awk).
#
# breseq must recognize the library, warn about it once, cut those reads at the adapter (keeping
# the 22 bases before it), and predict what lambda_mixed_pop predicts: in particular no junction
# made of adapter, which is what these reads produced before, at about 100 times the usual cost of
# the junction re-alignment step.
#
# The rewritten reads are written under output.reads/, which "clean" removes with the rest.
TESTCMD="\
    mkdir -p ${SELF}/output.reads \
    && gunzip -c ${DATADIR}/lambda/lambda_mixed_population.fastq.gz \
       | awk -f ${SELF}/mate_pair_reads.awk > ${SELF}/output.reads/lambda_mate_pair.fastq \
    && ${BRESEQ} \
    ${BRESEQ_TEST_THREAD_ARG} \
    -o ${SELF} \
    ${REFERENCE_ARG} \
    ${SELF}/output.reads/lambda_mate_pair.fastq \
    2> ${SELF}/output.stderr.txt \
    && test \$(grep -c 'looks like a MATE-PAIR library' ${SELF}/output.stderr.txt) -eq 1 \
    && grep -q 'Mate-pair junction adapter trimmed from reads' ${SELF}/output.stderr.txt \
    && rm -f ${SELF}/output.reads/lambda_mate_pair.fastq \
    "

do_test $1 ${SELF}
