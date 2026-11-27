#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
TEST_CORES=4
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/data/annotated.gd"
EXPECTED_OUTPUTS[0]="${SELF}/expected.gd"

# Gene conversion (CON) and deletion predicted from missing coverage ALONE, between paralogs that
# nothing in the reference file says are related. Four cases in one run; they were four tests
# (gene_conversion_1, gene_conversion_3, gene_conversion_4 and homologous_deletion_1), each a whole
# breseq run on the same kind of reads to make one call. check.sh states what the golden must mean,
# case by case.
#
# reference.gd builds the reference: lambda with four 1200 bp blocks each copied elsewhere with 4%
# of its bases changed (48 columns), so the genome carries four pairs of 96%-identical paralogs.
# The copies are written as plain SUBs, so there is NO repeat_region annotation on any of them --
# which is the whole point. No pair shares sequence with another pair.
#
#   pair A   4001-5200  and 12001-13200           (the copy every earlier version of this test used)
#   pair B  16001-17200 and 24001-25200, ANTIPARALLEL: the copy is written reverse-complemented
#   pair C  28001-29200 and 32001-33200
#   pair D  37001-38200 and 45001-46200
#
# input.gd then makes one event in each pair, and 36 bp single-end reads are simulated at 40x.
#
# (A) POSITIVE control for a gene conversion: CON 4401 400, taking pair A's second copy.
#     Every read from the converted tract that spans a column where the copies differ now matches
#     the second copy exactly and the first copy with a mismatch, so it is placed uniquely at the
#     DONOR. The recipient loses its coverage over the tract (an MC), the donor's doubles, and
#     nothing inside the tract is ever called as a mutation. What says this is a conversion and not
#     a deletion is the reads themselves: at each column where the copies differ, the donor column
#     carries twice the depth it should, all of it the donor allele, while a deleted recipient copy
#     would leave the donor's depth unchanged.
#     Expected: one CON whose tract is the whole identical stretch between the discriminating
#     columns that flank the converted 400 bp, with the donor region 8000 bp downstream, resting on
#     the MC.
#
# (B) The same with the donor on the OTHER strand: CON 16401 400, taking pair B's second copy read
#     backwards (region start > end). What is being checked is that the donor is found on the
#     reverse strand, that the tract is extended along it in the right direction, and that the
#     CON's region is written backwards.
#
# (C) NEGATIVE control for gene conversion: DEL 28001 1200, the whole first copy of pair C.
#     The recipient loses its coverage just as a conversion would, and every read from the
#     surviving copy still says "donor allele" at every column where the copies differ. The one
#     thing that distinguishes this from a conversion is depth: the donor column carries its normal
#     single-copy depth, not double. The gene conversion step must leave this alone, and the
#     deletion must be predicted as the DEL it is.
#
# (D) POSITIVE control for a large deletion from MC evidence ALONE, with no JC and no DP:
#     DEL 37582 8000, exactly the interval between pair D's copies. predictMCplusJCtoDEL's case (2)
#     handles deletions between ANNOTATED repeats; this is the case where nothing in the reference
#     file says the two copies are related. The crossover, 37582, is deliberately one of the 48
#     columns where the two copies differ, so the answer is a single position rather than an
#     ambiguity interval, and the test asserts base resolution rather than "close enough".
#     The deletion leaves no junction to find: a 36 bp read crossing the crossover aligns to either
#     copy, so candidate junction identification has nothing to split, and reads are single-end so
#     there is no discordant pair either. Coverage is all that is left, and MC alone is what has to
#     carry it. What locates the crossover is the reads themselves -- at each column where the
#     copies differ, coverage has moved from one copy to the other, and the place where that flips
#     is the breakpoint.
#
# tests/gene_conversion_2 is case (A) seen through paired 100 bp reads.
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
    && bash ${SELF}/check.sh ${SELF}/data/annotated.gd \
	"

do_test $1 ${SELF}
