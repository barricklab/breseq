#!/bin/bash
#
# Post-run checks for lambda_contig_ref: that the five contigs of the -c reference actually reached
# CNery as ONE reference group.
#
# This cannot be an expected.gd comparison, and that is the whole reason the file exists. The five
# contigs here are equal slices of one lambda genome at one depth, so a shared baseline and five
# per-contig baselines call exactly the same copy number and produce an identical .gd. The failure
# being pinned is silent by construction: if the `group` column went blank -- a seq_id that no longer
# matches a coverage table's name, a grouping derived from the wrong set -- CNery would fall back to
# refitting a baseline per contig, every contig of a real draft assembly would come out at copy
# number 1 however amplified it was, and every golden in the suite would still pass.
#
# Checked through breseq's own data/summary.json rather than CNery's output, because that file
# survives the run (09_copy_number_variation/ does not, absent -k) and because it makes this a test
# of the handoff end to end: breseq writes the group table, CNery groups on it, breseq reads the
# result back.
#
# Exits non-zero with a message naming what went wrong; testcmd.sh chains this after breseq with &&.

set -u

SELF="$1"
SUMMARY="${SELF}/data/summary.json"

fail() { echo "CHECK FAILED: $*"; exit 1; }

[ -f "${SUMMARY}" ] || fail "breseq wrote no summary: ${SUMMARY}"

# CNery declines the ori-ter correction for a group of two or more -- contig order and orientation
# are unknown, so there is no coordinate for a replication ramp to run along -- and says so, naming
# the group and its size. That string is therefore the direct evidence that the grouping arrived:
# an ungrouped contig would report a fit or "No ori-ter bias detected" instead.
#
# The group is named for the -c file, by basename, so this does not depend on where the run started.
expected='"correction_type": "No OTR correction (reference group '\''lambda-contig.gbk'\'': 5 sequences, order and orientation unknown)"'

#
# The seq_id keys of the copy_number block sit at four spaces; the same ids appear again, at six,
# under "references" -- so the indentation in the pattern is what selects the right section, not
# decoration.
for seq_id in NC_001416-0 NC_001416-1 NC_001416-2 NC_001416-3 NC_001416-4; do
    block=`sed -n "/^    \"${seq_id}\": {/,/^    }/p" "${SUMMARY}"`
    [ -n "${block}" ] || fail "no copy number summary for ${seq_id} in ${SUMMARY}"

    echo "${block}" | grep -qF "${expected}" \
        || fail "${seq_id} was not grouped with the other contigs of lambda-contig.gbk -- CNery reported: `echo "${block}" | sed -n 's/.*"correction_type": //p'`"
done

echo "Reference group checks passed."
