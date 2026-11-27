#!/bin/bash
#
# Polymorphism-mode deletions and amplifications, each carried by only part of the population, in
# one simulated mixture. This was two tests (lambda_polymorphism_del and lambda_polymorphism_amp)
# with two breseq runs; the four mutations are far apart on lambda and do not interact, so one run
# asserts all of it. check.sh states what the golden must mean, mutation by mutation.
#
# A simulated population of three genotypes, so that each mutation keeps the frequency its own test
# gave it:
#   50%  carry del_amp.gd  (both deletions and both amplifications)   reads at 40x
#   25%  carry amp.gd      (the amplifications only)                  reads at 20x
#   25%  are the reference                                            reads at 20x
#
# Deletions at 50% that leave no missing coverage, because only part of the population carries
# them, predicted from their junction (and, when CN can resolve them, the relative coverage):
#   DEL 10001 3000   CN sees it at half copy; frequency from jc_frequency and cn_frequency.
#   DEL 30001 80     too small for CN's 100 bp windows; frequency from jc_frequency alone.
#
# Amplifications at 75% with a calculated frequency (barricklab/conda#2: the AMPs built from JC +
# CN evidence carried no frequency, and the predictor asserted):
#   AMP 20001 5000 x2   new_copy_number 2.
#   AMP 35001 3000 x3   new_copy_number 3 -- the smallest whole copy number that keeps the
#                       frequency <= 1 for a mean copy number of 2.5.
# Each with the jc_frequency and cn_frequency they were combined from.
#
# The deletions are NOT simply put at 75% alongside the amplifications in a two-genotype mixture:
# that was tried, and at 75% the 80 bp deletion empties enough of its one 100 bp window for CN to
# call it, which loses the junction-only case this test exists to cover.
#
# The reads are generated at test time, into output_sim/ (removed by do_clean), so no read data is
# committed. The seeds are not arbitrary: an AMP's jc_frequency is x / (1 - x) of its junction's
# frequency x, which near x = 0.43 turns a few reads of sampling noise into a large swing (0.62,
# 0.89 and 1.02 across three seed choices for the 5 kb AMP). With seeds 2 and 3 for the two mutated
# genomes the junction drew high, jc_frequency passed 1 and the combined frequency clamped to 1;
# seeds 4 and 5 give a draw near the expected value. If a rebuild moves an AMP frequency out of
# check.sh's tolerance, look at that junction's read counts before suspecting the predictor.
#
# The reference-genome reads are deliberately named reference_*: stage 02 used to remove its bowtie2
# index with "rm -f 02_reference_alignment/reference*", which also deleted such a read file's
# reference_reads.stage1.bam (see track_bowtie2_index_files in breseq_cmdline.cpp).
#

SELF=`dirname ${BASH_SOURCE}`
TEST_CORES=4
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/data/annotated.gd"
EXPECTED_OUTPUTS[0]="${SELF}/expected.gd"
REFERENCE_ARG="-r ${DATADIR}/lambda/lambda.gbk"
SIMDIR="${SELF}/output_sim"

TESTCMD="\
    mkdir -p ${SIMDIR} \
    && ${GDTOOLS} APPLY \
        -f FASTA \
        ${REFERENCE_ARG} \
        -o ${SIMDIR}/amp.fasta \
        ${SELF}/amp.gd \
    && ${GDTOOLS} APPLY \
        -f FASTA \
        ${REFERENCE_ARG} \
        -o ${SIMDIR}/del_amp.fasta \
        ${SELF}/del_amp.gd \
    && ${BRESEQ} SIMULATE-READS \
        ${REFERENCE_ARG} \
        -c 20 -l 100 --seed 1 \
        -o ${SIMDIR}/reference_reads.fastq \
    && ${BRESEQ} SIMULATE-READS \
        -r ${SIMDIR}/amp.fasta \
        -c 20 -l 100 --seed 4 \
        -o ${SIMDIR}/amp_reads.fastq \
    && ${BRESEQ} SIMULATE-READS \
        -r ${SIMDIR}/del_amp.fasta \
        -c 40 -l 100 --seed 5 \
        -o ${SIMDIR}/del_amp_reads.fastq \
    && ${BRESEQ} \
        ${BRESEQ_TEST_THREAD_ARG} \
        -p \
        -o ${SELF} \
        ${REFERENCE_ARG} \
        ${SIMDIR}/reference_reads.fastq \
        ${SIMDIR}/amp_reads.fastq \
        ${SIMDIR}/del_amp_reads.fastq \
    && bash ${SELF}/check.sh ${SELF}/data/annotated.gd \
    "

do_test $1 ${SELF}
