#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
TEST_CORES=4
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/data/annotated.gd"
EXPECTED_OUTPUTS[0]="${SELF}/expected.gd"

# Read linkage (LN evidence) and cluster-local realignment on a mixture of KNOWN composition.
#
# Three genomes are simulated and mixed 50/30/20 at 50x total with 100-bp reads: lambda itself,
# lambda carrying mutA.gd, and lambda carrying mutB.gd. Nothing is committed but the two small .gd
# files; the reads are regenerated on every run (SIMULATE-READS is byte-reproducible for a seed).
#
# mutA (30%), each chosen to exercise one thing the column-by-column RA caller could not do:
#   INS AC at 10000        a 2-base insertion at a non-repeat site: two adjacent RA columns
#                          (10000.1, 10000.2) that must merge into ONE INS at ~30%, not two.
#   INS GCA at 11701       one more unit of the GCAGCAGCA repeat at 11693-11701: where the aligner
#                          puts the gap differs between reads, so the CIGAR shifter and the
#                          realignment against the haplotype sequence must agree on one call.
#   SNP 12000 A>G, 12011 A>C   11 bp apart in cis: reported as a phase=cis LN.
#   SUB 13020 T>CA         a substitution plus an insertion in one read: the INS + SNP pair that
#                          polymorphism mode used to report as two events.
#   SNP 15000 G>A          a lone control.
#   DEL 17240 3            one unit of the ACGACGACG repeat at 17240-17248.
#   DEL 22368 1            one A of the AAAAAAA homopolymer at 22368-22374.
# mutB (20%):
#   SNP 12005 T>C          5 bp from mutA's 12000, on the OTHER genotype: a phase=trans LN.
#   SNP 16000 T>G          a lone control.
#
# expected.gd pins the exact fitted values; check.sh states what the numbers must mean, so a rebuild
# after a change to the simulator or to a default cannot silently lose them. Coordinates in check.sh
# are the right-normalized ones breseq reports (11702, 17246, 22375), not the ones mutA.gd applies.
TESTCMD=" \
    ${GDTOOLS} APPLY \
        -f FASTA \
        -s NC_001416 \
        -r ${DATADIR}/lambda/lambda.gbk \
        -o ${SELF}/output.mutA.fna \
        ${SELF}/mutA.gd \
    && ${GDTOOLS} APPLY \
        -f FASTA \
        -s NC_001416 \
        -r ${DATADIR}/lambda/lambda.gbk \
        -o ${SELF}/output.mutB.fna \
        ${SELF}/mutB.gd \
    && ${BRESEQ} SIMULATE-READS \
        -r ${DATADIR}/lambda/lambda.gbk \
        -o ${SELF}/output.ref.fastq \
        -l 100 \
        -c 25 \
        --seed 1 \
    && ${BRESEQ} SIMULATE-READS \
        -r ${SELF}/output.mutA.fna \
        -o ${SELF}/output.mutA.fastq \
        -l 100 \
        -c 15 \
        --seed 2 \
    && ${BRESEQ} SIMULATE-READS \
        -r ${SELF}/output.mutB.fna \
        -o ${SELF}/output.mutB.fastq \
        -l 100 \
        -c 10 \
        --seed 3 \
    && ${BRESEQ} \
        ${BRESEQ_TEST_THREAD_ARG} \
        --polymorphism-prediction \
        -o ${SELF} \
        -r ${DATADIR}/lambda/lambda.gbk \
        ${SELF}/output.ref.fastq \
        ${SELF}/output.mutA.fastq \
        ${SELF}/output.mutB.fastq \
    && bash ${SELF}/check.sh ${SELF}/data/annotated.gd \
	"

do_test $1 ${SELF}
