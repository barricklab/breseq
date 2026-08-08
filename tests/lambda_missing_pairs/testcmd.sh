#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
TEST_CORES=4
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/data/annotated.gd"
EXPECTED_OUTPUTS[0]="${SELF}/expected.gd"
REFERENCE_ARG="-r ${DATADIR}/lambda/lambda.gbk"

# POSITIVE control for --predict-missing-pairs.
#
# A 3274 bp sequence with no homology to lambda is inserted into the reference, paired reads are
# simulated from the MUTATED genome, and breseq is then run against the UNMUTATED lambda. Mates that
# land inside the insertion have nowhere to align -- the sequence is in neither the reference nor any
# candidate junction -- so they become singletons, and MP must fire on both flanks of the insertion
# point. This is the case no other evidence type can see: JC needs a split read whose halves both
# map, DP needs both mates mapped.
#
# Compare with tests/tmv_plasmid_missing_pairs, the negative control, which checks that MP does NOT
# fire on real data whose unmappable mates are just a flat background.
#
# The read files are generated here rather than committed. SIMULATE-READS is exactly reproducible
# given --seed -- integer-only arithmetic, no libc rand(), no libm -- so the same bytes are produced
# on Linux and macOS. expected.gd asserts that directly: DIFF_IGNORE does not strip the
# #=CONVERTED-BASES / #=INPUT-READS header lines, so they are compared. That also means any change to
# the simulator or its defaults requires rebuilding this expected.gd.
TESTCMD=" \
    ${GDTOOLS} APPLY \
        -f FASTA \
        -s NC_001416 \
        -r ${DATADIR}/lambda/lambda.gbk \
        -r ${DATADIR}/tmv_plasmid/tmv-plasmid-truncate-start.gbk \
        -o ${SELF}/output.mutated.fna \
        ${SELF}/input.gd \
    && ${BRESEQ} SIMULATE-READS \
        -m paired-end \
        -r ${SELF}/output.mutated.fna \
        -o ${SELF}/output.simulated.fastq \
        -l 100 \
        -c 40 \
        --mean 350 \
        --stdev 40 \
        --seed 1 \
    && ${BRESEQ} \
        ${BRESEQ_TEST_THREAD_ARG} \
        --predict-missing-pairs \
        -o ${SELF} \
        ${REFERENCE_ARG} \
        ${SELF}/output.simulated_1.fastq \
        ${SELF}/output.simulated_2.fastq \
	"

do_test $1 ${SELF}
