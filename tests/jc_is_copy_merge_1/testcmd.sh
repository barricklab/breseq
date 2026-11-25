#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
TEST_CORES=3
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/data/annotated.gd"
EXPECTED_OUTPUTS[0]="${SELF}/expected.gd"

# Two candidate junctions that resolve to ONE JC entry, both PASSING.
#
# One unique side (AP025630:21294, -1) is joined to the 5' end of IS630_443 through two different
# copies (ending at 33786 and 36808; see tests/data/annamia_is630/README). The copies differ near
# that end, so the candidate junction sequences differ and each is scored on its own reads. Then
# junction_to_diff_entry() normalizes each redundant IS side onto the consensus copy (ending at
# 3353), and both come out as JC AP025630 21294 -1 AP025630 3353 -1 0. Before the fix this stopped
# breseq in resolve_alignments with "Fatal formatting error in GD file being written:
# .../05_alignment_correction/jc_evidence.gd".
#
# Now one entry is kept per identity. Here both candidates pass, so the reads resolution credited to
# the weaker one count toward the survivor: its key appears in the survivor's merged_keys, and
# new_junction_read_count and the frequency include those reads. jc_is_copy_merge_2 covers the case
# where the weaker candidate is rejected.
TESTCMD="\
    ${BRESEQ} \
    ${BRESEQ_TEST_THREAD_ARG} \
    -o ${SELF} \
    -r ${DATADIR}/annamia_is630/AP025630_is630_windows.gbk.gz \
    ${DATADIR}/annamia_is630/DRR285493_subset_1.fastq.gz \
    ${DATADIR}/annamia_is630/DRR285493_subset_2.fastq.gz \
    "

do_test $1 ${SELF}
