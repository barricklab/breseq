#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
TEST_CORES=1
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/output.reads_1.fastq"
CURRENT_OUTPUTS[1]="${SELF}/output.reads_2.fastq"
EXPECTED_OUTPUTS[0]="${SELF}/expected.reads_1.fastq"
EXPECTED_OUTPUTS[1]="${SELF}/expected.reads_2.fastq"

# Asserts that SIMULATE-READS is exactly reproducible for a given --seed, by comparing the generated
# FASTQ itself rather than anything downstream of it.
#
# This matters because other tests now generate their read files at run time instead of committing
# them, and CI runs on both Linux and macOS. If the generator ever stops being platform-independent
# -- someone reintroduces rand(), a <random> distribution, or a libm call into the read path -- this
# test fails directly and says so, whereas the tests that go on to run breseq would only show an
# unexplained annotated.gd difference that looks like an aligner change.
#
# Kept deliberately tiny (25 pairs, ~6 KB per file) so it runs in well under a second.
TESTCMD=" \
    ${BRESEQ} SIMULATE-READS \
        -m paired-end \
        -r ${DATADIR}/lambda/lambda.gbk \
        -o ${SELF}/output.reads.fastq \
        -l 100 \
        -n 25 \
        --mean 350 \
        --stdev 40 \
        --seed 1 \
	"

do_test $1 ${SELF}
