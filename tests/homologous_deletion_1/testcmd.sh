#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
TEST_CORES=4
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/data/annotated.gd"
EXPECTED_OUTPUTS[0]="${SELF}/expected.gd"

# POSITIVE control for predicting a large deletion from MC evidence ALONE, with no JC and no DP.
#
# reference.gd builds the reference: lambda, with a 1200 bp block at 4001-5200 copied to 12001-13200
# and 4% of its bases changed, so the genome carries two 96%-identical paralogs exactly 8000 bp apart.
# The copy is written as a plain SUB, so there is NO repeat_region annotation on it -- which is the
# whole point. predictMCplusJCtoDEL's case (2) handles deletions between ANNOTATED repeats; this is
# the case where nothing in the reference file says the two copies are related.
#
# input.gd then deletes exactly the interval between them, crossing over at 4593 -- deliberately one
# of the 48 columns where the two copies differ, so the answer is a single position rather than an
# ambiguity interval, and the test asserts base resolution rather than "close enough".
#
# The deletion leaves no junction to find: a 36 bp read crossing the crossover aligns to either copy,
# so candidate junction identification has nothing to split, and reads are single-end so there is no
# discordant pair either. Coverage is all that is left, and MC alone is what has to carry it. What
# locates the crossover is the reads themselves -- at each column where the copies differ, coverage
# has moved from one copy to the other, and the place where that flips is the breakpoint.
#
# Expected: DEL NC_001416 4593 8000, resting on the single MC.
#
# Reference and reads are both generated here rather than committed; SIMULATE-READS is exactly
# reproducible given --seed, so the same bytes are produced on Linux and macOS. That also means any
# change to the simulator or its defaults requires rebuilding this expected.gd.
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
	"

do_test $1 ${SELF}
