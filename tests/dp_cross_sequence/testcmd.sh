#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
TEST_CORES=4
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/data/annotated.gd"
EXPECTED_OUTPUTS[0]="${SELF}/expected.gd"

# POSITIVE control for CROSS-SEQUENCE discordant pair (DP) evidence.
#
# Paired reads are simulated from the intact, single-sequence lambda genome, and breseq is then run
# against lambda-contig.gbk -- that same genome split in order into 5 contigs (NC_001416-0..4, ~9700
# bp each; they are exact contiguous slices of lambda.gbk, differing only where two N's were resolved
# to real bases). Every read pair whose unsequenced middle gap straddles one of the 4 internal contig
# boundaries therefore has its two mates on DIFFERENT reference sequences.
#
# That is the one thing this test is about, and it is what breseq could not see at all: a pair whose
# mates land on different reference sequences has no within-sequence orientation, so mark_pair_info
# stamps it XP:Z:NA with TLEN 0, and the DP seeding step in identify_mutations.cpp dropped every read
# carrying that tag on its `else continue`. No chromosome-to-plasmid integration and no translocation
# between two contigs could produce DP evidence, however many pairs supported it. Such pairs now get
# their own orientation slot (kDPnOrientations) instead of being discarded.
#
# Each internal boundary should therefore yield a DP item with side_1_seq_id != side_2_seq_id: side_1
# at the end of contig i (strand -1) and side_2 at the start of contig i+1 (strand +1). Before the
# change this test produced ZERO DP items.
#
# The library parameters follow lambda_pair_distances: at -l 100 --mean 600 and 60x, roughly
# (60/200)*(600-2*100) = 120 pairs bridge each boundary -- far above --discordant-pair-seed (3) and
# --discordant-pair-minimum-pairs (2). A contig boundary is a sequence END, so no concordant pair can
# bracket either placed side and the local frequency is ~1.0; both the frequency and skew gates are
# exercised on their merits rather than bypassed.
#
# -c rather than -r, matching lambda_contig_ref: these five sequences are contigs of one genome, as a
# draft assembly would present them. That also makes this the test that pins DP's deliberate absence
# from ignore_evidence_near_contig_ends -- RA/MC/JC/SC are ignored near a contig end because that is
# where their statistic stops being measurable, and PD/MP set their own equivalent, but a legitimate
# cross-sequence DP sits at two contig ends BY CONSTRUCTION. If that rule is ever extended to DP,
# every DP item here disappears and this test says so.
#
# The exact diff against expected.gd doubles as the false-positive control: a spurious DP elsewhere in
# the 5 contigs, or a cross-sequence pair mis-binned into an intra-contig region, shows up as an extra
# or altered line. The single-sequence DP tests (tmv_plasmid_circular_deletion, lambda_pair_distances
# and the long LTEE ones) are the "no same-sequence behavior changed" control and must not move.
#
# The read files are generated here rather than committed. SIMULATE-READS is exactly reproducible
# given --seed -- integer-only arithmetic, no libc rand(), no libm -- so the same bytes are produced
# on Linux and macOS. expected.gd asserts that directly: DIFF_IGNORE does not strip the
# #=CONVERTED-BASES / #=INPUT-READS header lines, so they are compared. Any change to the simulator or
# its defaults therefore requires rebuilding this expected.gd.
TESTCMD=" \
    ${BRESEQ} SIMULATE-READS \
        -m paired-end \
        -r ${DATADIR}/lambda/lambda.gbk \
        -o ${SELF}/output.simulated.fastq \
        -l 100 \
        -c 60 \
        --mean 600 \
        --stdev 60 \
        --seed 1 \
    && ${BRESEQ} \
        ${BRESEQ_TEST_THREAD_ARG} \
        --predict-discordant-pairs \
        -o ${SELF} \
        -c ${DATADIR}/lambda/lambda-contig.gbk \
        ${SELF}/output.simulated_1.fastq \
        ${SELF}/output.simulated_2.fastq \
	"

do_test $1 ${SELF}
