#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
TEST_CORES=4
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/data/annotated.gd"
EXPECTED_OUTPUTS[0]="${SELF}/expected.gd"
REFERENCE_ARG="-r ${DATADIR}/REL606/REL606.fragment.gbk"

# The read file is GENERATED here rather than committed: input.gd is applied to the reference to make
# the mutated genome, and breseq's own tiled simulator turns that into reads. Tiled mode uses no
# randomness at all -- fixed spacing, both strands at every start, constant quality, no simulated
# errors -- so the FASTQ is identical on every platform and run. This replaces a committed 1.1 MB
# REL606.tiled_reads_continuation.fastq.gz, which these exact commands reproduce byte for byte.
#
# -c 100 gives spacing 1 (ceil(2 * read_length / coverage)), which is what the committed reads used.
TESTCMD=" \
    ${GDTOOLS} APPLY \
        -f FASTA \
        -r ${DATADIR}/REL606/REL606.fragment.gbk \
        -o ${SELF}/output.mutated.fna \
        ${SELF}/input.gd \
    && ${BRESEQ} SIMULATE-READS \
        -m tiled \
        -l 50 \
        -c 100 \
        -r ${SELF}/output.mutated.fna \
        -o ${SELF}/output.tiled_reads.fastq \
    && ${BRESEQ} \
        ${BRESEQ_TEST_THREAD_ARG} \
        -p \
        -o ${SELF} \
        ${REFERENCE_ARG} \
        -s ${DATADIR}/REL606/REL606.is.gbk \
        ${SELF}/output.tiled_reads.fastq \
	"

do_test $1 ${SELF}
