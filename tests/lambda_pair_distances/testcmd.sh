#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
TEST_CORES=4
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/data/annotated.gd"
EXPECTED_OUTPUTS[0]="${SELF}/expected.gd"
REFERENCE_ARG="-r ${DATADIR}/lambda/lambda.gbk"

# POSITIVE control for --predict-pair-distance.
#
# A 200 bp deletion and a 200 bp insertion are applied to lambda, paired reads are simulated from the
# MUTATED genome, and breseq is then run against the UNMUTATED lambda. Read pairs whose unsequenced
# middle gap spans the deletion map 200 bp farther apart than the library says they should; pairs
# spanning the insertion map 200 bp closer. No individual pair is anomalous either way -- both shifts
# sit well inside the discordant distance cutoff, and there is no short-side cutoff at all -- so this
# is exactly the case DP cannot see. PD tests the collective shift instead.
#
# The library parameters are not arbitrary. A pair can only place a 200 bp event in its gap if the
# fragment is longer than 2*read_length + 200 = 400, so --mean 600 leaves the great majority of pairs
# informative. At --mean 350 (what lambda_missing_pairs uses) almost no fragment could span either
# event with both mates still aligned, and PD would correctly find nothing. That is the sensitivity
# limit worth remembering: PD reaches insertions up to roughly (median distance - 2*read length).
#
# --predict-discordant-pairs is on as well, deliberately. The deletion pushes a small tail of its
# spanning pairs past the discordant cutoff, so DP does seed there, and expected.gd is where the rule
# that PD supersedes DP on a shared breakpoint gets asserted.
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
        -c 60 \
        --mean 600 \
        --stdev 60 \
        --seed 1 \
    && ${BRESEQ} \
        ${BRESEQ_TEST_THREAD_ARG} \
        --predict-pair-distance \
        --predict-discordant-pairs \
        -o ${SELF} \
        ${REFERENCE_ARG} \
        ${SELF}/output.simulated_1.fastq \
        ${SELF}/output.simulated_2.fastq \
	"

do_test $1 ${SELF}
