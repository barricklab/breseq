#!/bin/bash
#
# Post-run checks for lambda_mixed_pop_cn_no_coverage. See testcmd.sh for why these live here rather
# than in the expected .gd.
#
# Exits non-zero with a message naming what went wrong; testcmd.sh chains this after breseq with &&.

set -u

SELF="$1"
COVDIR="${SELF}/08_mutation_identification"

fail() { echo "CHECK FAILED: $*"; exit 1; }

# The reference that got the reads, and the one that got none. Both must have a table, and both must
# have a row per reference position -- a header-only table is the bug this test exists for.
#   NC_001416   lambda.gbk, 48502 bp
#   AF322221    other.gbk,    687 bp
for entry in "NC_001416 48502" "AF322221 687"; do
    set -- $entry
    seq_id="$1"
    length="$2"
    table="${COVDIR}/${seq_id}.coverage.tsv"

    [ -f "${table}" ] || fail "no coverage table for ${seq_id}: ${table}"

    lines=`grep -c "" "${table}"`
    expected=$((length + 1))          # one header row plus one row per position
    [ "${lines}" = "${expected}" ] \
        || fail "${seq_id}.coverage.tsv has ${lines} line(s), expected ${expected} (header + ${length} positions)"
done

# AF322221 got no reads, so every one of its rows must read zero coverage. This is what distinguishes
# "the table was written" from "the table was written with the wrong sequence's numbers in it".
nonzero=`awk 'NR > 1 && ($3 != 0 || $4 != 0)' "${COVDIR}/AF322221.coverage.tsv" | grep -c ""`
[ "${nonzero}" = "0" ] || fail "AF322221 has no reads but ${nonzero} of its coverage rows are non-zero"

# The junction-only reference is deliberately not analyzed, so it must have NO table -- and CNery
# must not have been asked for one, which is what the run reaching this point at all proves.
[ ! -f "${COVDIR}/REL606-5.coverage.tsv" ] \
    || fail "junction-only reference REL606-5 should not have a coverage table"

# CNery must have produced its per-reference output for both analyzed sequences, and for neither
# more nor fewer. A missing OTR record costs breseq all of its ori-ter reporting for that sequence.
OTRDIR="${SELF}/09_copy_number_variation/cnery_out/OTR_corr"
for seq_id in NC_001416 AF322221; do
    [ -f "${OTRDIR}/cnery_out${seq_id}_otr_results.json" ] \
        || fail "no CNery OTR results for ${seq_id}"
done
[ ! -f "${OTRDIR}/cnery_outREL606-5_otr_results.json" ] \
    || fail "CNery analyzed the junction-only reference REL606-5"

# ...and CNery must have recognized AF322221 as having nothing to measure, rather than fitting noise.
grep -q '"No usable coverage reason": "' "${OTRDIR}/cnery_outAF322221_otr_results.json" \
    || fail "CNery did not report AF322221 as having no usable coverage"

# breseq must actually READ that ori-ter fit, which nothing else in the suite checks: no .gd file
# carries an origin or a terminus, so when CNery renamed "Origin-to-Termius/Bias Ratio" to
# "...Terminus..." every golden kept passing while breseq silently reported "no ori-ter bias
# detected" for every reference sequence. Assert breseq's summary against CNery's own numbers rather
# than against a hard-coded verdict, so this stays a test of the reader and not of the fit.
CNJSON="${OTRDIR}/cnery_outNC_001416_otr_results.json"
SUMMARY="${SELF}/09_copy_number_variation/copy_number_summary.json"

[ -f "${SUMMARY}" ] || fail "breseq wrote no copy number summary: ${SUMMARY}"

# The summary holds one object per reference sequence, keyed by seq_id and SORTED -- so AF322221
# comes first and a whole-file grep would answer for the wrong sequence. Slice out NC_001416's block
# and ask only that.
BLOCK=`sed -n '/^  "NC_001416": {/,/^  }/p' "${SUMMARY}"`
[ -n "${BLOCK}" ] || fail "no NC_001416 block in ${SUMMARY}"

# "key": value  ->  value, with the trailing comma removed. Reads stdin. The sed delimiter is '|'
# rather than '/' because one of the keys asked for is "Origin-to-Terminus/Bias Ratio".
json_value() { sed -n "s|.*\"$1\"[[:space:]]*:[[:space:]]*\(.*\)|\1|p" | sed 's/,$//' | head -1; }

cnery_ratio=`json_value "Origin-to-Terminus/Bias Ratio" < "${CNJSON}"`
breseq_detected=`echo "${BLOCK}" | json_value "otr_detected"`

case "${cnery_ratio}" in
    '"Not detected"'|null|'')
        [ "${breseq_detected}" = "false" ] \
            || fail "CNery reported no ori-ter fit for NC_001416 but breseq claims one"
        ;;
    *)
        [ "${breseq_detected}" = "true" ] \
            || fail "CNery fit an ori-ter ramp for NC_001416 (ratio ${cnery_ratio}) but breseq did not read it -- has the ratio key been renamed again?"

        # The coordinates are integers, so they compare exactly; the ratio is a float breseq
        # reformats, so it is deliberately not string-compared here.
        cnery_ori=`json_value "Origin window" < "${CNJSON}"`
        cnery_ter=`json_value "Terminus window" < "${CNJSON}"`
        breseq_ori=`echo "${BLOCK}" | json_value "origin"`
        breseq_ter=`echo "${BLOCK}" | json_value "terminus"`
        [ "${cnery_ori}" = "${breseq_ori}" ] \
            || fail "breseq read origin ${breseq_ori}, CNery wrote 'Origin window' ${cnery_ori}"
        [ "${cnery_ter}" = "${breseq_ter}" ] \
            || fail "breseq read terminus ${breseq_ter}, CNery wrote 'Terminus window' ${cnery_ter}"
        ;;
esac

# The fields added alongside the rename. relative_copy_number is 1.0 for the longest (here, the only
# analyzed) sequence by construction, so a zero means breseq did not read the key at all.
[ "`echo "${BLOCK}" | json_value relative_copy_number`" = "1.0" ] \
    || fail "breseq did not read CNery's relative copy number for NC_001416"
echo "${BLOCK}" | grep -q '"correction_type": "Ori-ter coordinates' \
    || fail "breseq did not read CNery's correction type for NC_001416"

# ...and the same three for the sequence with no coverage, where they are the ONLY thing reported:
# there is no fit, so if these do not come through, breseq says nothing at all about that sequence.
NOCOV=`sed -n '/^  "AF322221": {/,/^  }/p' "${SUMMARY}"`
echo "${NOCOV}" | grep -q '"no_coverage_reason": "..*"' \
    || fail "breseq did not read CNery's no-usable-coverage reason for AF322221"
echo "${NOCOV}" | grep -q '"correction_type": "No usable coverage"' \
    || fail "breseq did not read CNery's correction type for AF322221"

echo "Coverage table checks passed."
