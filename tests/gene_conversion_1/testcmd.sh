#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
TEST_CORES=4
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/data/annotated.gd"
EXPECTED_OUTPUTS[0]="${SELF}/expected.gd"

# POSITIVE control for predicting a gene conversion (CON) from missing coverage ALONE.
#
# reference.gd is the one homologous_deletion_1 uses: lambda with a 1200 bp block at 4001-5200
# copied to 12001-13200 and 4% of its bases changed, so the genome carries two 96%-identical
# paralogs 8000 bp apart with no repeat_region annotation. input.gd converts 400 bp in the middle of
# the first copy to the second copy's sequence.
#
# With 36 bp single-end reads, every read from the converted tract that spans a column where the
# copies differ now matches the second copy exactly and the first copy with a mismatch, so it is
# placed uniquely at the DONOR. The recipient loses its coverage over the tract (an MC), the donor's
# doubles, and nothing inside the tract is ever called as a mutation. What says this is a
# conversion and not a deletion is the reads themselves: at each column where the copies differ, the
# donor column carries twice the depth it should, all of it the donor allele, while a deleted
# recipient copy would leave the donor's depth unchanged.
#
# Expected: one CON whose tract is the whole identical stretch between the discriminating columns
# that flank the converted 400 bp, with the donor region 8000 bp downstream, resting on the MC.
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
    && test \$(cut -f1 ${SELF}/data/annotated.gd | grep -c '^CON\$') -eq 1 \
	"

do_test $1 ${SELF}
