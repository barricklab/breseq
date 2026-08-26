#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
TEST_CORES=4
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/data/annotated.gd"
EXPECTED_OUTPUTS[0]="${SELF}/expected.gd"

# REGRESSION test for an IS-mediated deletion whose MC boundary stops SHORT of the IS element.
#
# predictMCplusJCtoDEL's case (3) pairs an MC boundary lying inside a repeat with a junction from
# that repeat family to the other MC boundary. It used to locate the repeat with within_repeat(),
# strict containment -- so an MC boundary that fell a few bases short of the element made both
# boundaries look like unique sequence, case (2b) claimed the MC, found only one of the two junctions
# it needs, and the deletion was dropped entirely. This was seen on real data (an IS150-mediated
# deletion in E. coli REL606 missed by 4 bp).
#
# input.gd is the same four mutations as REL606_fragment_is_mediated_dels. The second of them,
# DEL 28183 3119, ends exactly at the start of the IS1 copy at 31302..32069, so its MC normally runs
# back INTO that element (see that test's MC 6 ... 28183 31540 0 238).
#
# spike.gd then supplies the shortfall the way real data does: reads from an excised-element
# transposition intermediate, which carries the element's 9 bp target-site duplication and therefore
# maps uniquely for 9 bases past the annotated boundary. That holds coverage over 31293..31301 and
# pushes the MC end back to 31292 -- 9 bp outside IS1, where within_repeat() sees nothing.
#
# Expected: DEL 28183 3119 mediated=IS1 is still predicted, its right boundary taken from the
# junction (base-pair resolution) rather than from the coverage. The other three mutations are
# carried along unchanged, so this also asserts the relaxation did not disturb them.
#
# Reference stays the stock fragment; both read sets are generated here rather than committed.
# SIMULATE-READS is exactly reproducible given --seed, so the same bytes are produced on Linux and
# macOS -- which also means any change to the simulator requires rebuilding this expected.gd.
TESTCMD=" \
    ${GDTOOLS} APPLY \
        -f FASTA \
        -r ${DATADIR}/REL606/REL606.fragment.gbk \
        -o ${SELF}/output.mutated.fna \
        ${SELF}/input.gd \
    && ${GDTOOLS} APPLY \
        -f FASTA \
        -r ${DATADIR}/REL606/REL606.fragment.gbk \
        -o ${SELF}/output.spike.fna \
        ${SELF}/spike.gd \
    && ${BRESEQ} SIMULATE-READS \
        -r ${SELF}/output.mutated.fna \
        -o ${SELF}/output.simulated.fastq \
        -l 36 \
        -c 40 \
        --seed 1 \
    && ${BRESEQ} SIMULATE-READS \
        -r ${SELF}/output.spike.fna \
        -o ${SELF}/output.spike.fastq \
        -l 36 \
        -c 400 \
        --seed 2 \
    && ${BRESEQ} \
        ${BRESEQ_TEST_THREAD_ARG} \
        -b 0 \
        -o ${SELF} \
        -r ${DATADIR}/REL606/REL606.fragment.gbk \
        ${SELF}/output.simulated.fastq \
        ${SELF}/output.spike.fastq \
	"

do_test $1 ${SELF}
