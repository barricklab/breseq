#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
TEST_CORES=4
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/data/annotated.gd"
EXPECTED_OUTPUTS[0]="${SELF}/expected.gd"
REFERENCE_ARG="-r ${DATADIR}/lambda/lambda.gbk"

# POSITIVE control for pair distance (PD) AND missing pair (MP) evidence, in one run. They were two
# tests (this one and lambda_missing_pairs) with two simulated libraries; the events are 12 kb apart
# and each evidence type looks at a different set of pairs, so one library serves both. check.sh
# states what the golden must mean for each.
#
# PD: a 200 bp deletion and a 200 bp insertion are applied to lambda, paired reads are simulated
# from the MUTATED genome, and breseq is then run against the UNMUTATED lambda. Read pairs whose
# unsequenced middle gap spans the deletion map 200 bp farther apart than the library says they
# should; pairs spanning the insertion map 200 bp closer. No individual pair is anomalous either way
# -- both shifts sit well inside the discordant distance cutoff, and there is no short-side cutoff
# at all -- so this is exactly the case DP cannot see. PD tests the collective shift instead.
#
# MP: a 3074 bp sequence with no homology to lambda is inserted as well. Mates that land inside it
# have nowhere to align -- the sequence is in neither the reference nor any candidate junction -- so
# they become singletons, and MP must fire on both flanks of the insertion point. This is the case
# no other evidence type can see: JC needs a split read whose halves both map, DP needs both mates
# mapped. The 200 bp insertion leaves singletons too, but only from the few pairs with a mate wholly
# inside 200 bp, so its two MP items score right at --missing-pair-score-cutoff (3) and can land on
# either side of it; check.sh asserts nothing about them. That event is PD's to call.
# tests/tmv_plasmid_circular_deletion is the negative control on real data, where the unmappable
# mates are just a flat background.
#
# The library parameters are not arbitrary. A pair can only place a 200 bp event in its gap if the
# fragment is longer than 2*read_length + 200 = 400, so --mean 600 leaves the great majority of pairs
# informative. At --mean 350 (what lambda_missing_pairs used) almost no fragment could span either
# event with both mates still aligned, and PD would correctly find nothing. That is the sensitivity
# limit worth remembering: PD reaches insertions up to roughly (median distance - 2*read length).
# MP has no such limit, which is why it could move to this library and PD could not move to that one.
#
# --predict-discordant-pairs is on as well, deliberately. The deletion pushes a small tail of its
# spanning pairs past the discordant cutoff, so DP does seed there, and expected.gd is where the rule
# that PD supersedes DP on a shared breakpoint gets asserted. (All three --predict-* flags here are
# deprecated and inert: these evidence types are predicted by default.)
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
        --predict-missing-pairs \
        -o ${SELF} \
        ${REFERENCE_ARG} \
        ${SELF}/output.simulated_1.fastq \
        ${SELF}/output.simulated_2.fastq \
    && bash ${SELF}/check.sh ${SELF}/data/annotated.gd \
	"

do_test $1 ${SELF}
