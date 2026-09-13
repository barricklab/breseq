#!/bin/bash
#
# Post-run checks for lambda_polymorphism_cn_evidence. See testcmd.sh for why these live here rather
# than in the expected .gd: on this read set polymorphism mode reaches the same segmentation the
# integer grid does, so the .gd cannot distinguish "breseq passed --polymorphism-mode" from "breseq
# did not". The dtype CNery wrote can.
#
# Exits non-zero with a message naming what went wrong; testcmd.sh chains this after breseq with &&.

set -u

SELF="$1"
CNVDIR="${SELF}/09_copy_number_variation/cnery_out/CNV_csv"
BREAKS="${CNVDIR}/cnery_outNC_001416_break_pts.csv"

fail() { echo "CHECK FAILED: $*"; exit 1; }

[ -f "${BREAKS}" ] || fail "no CNery break points file: ${BREAKS} (was CN prediction skipped?)"

# Startpos,State,Segment_Size -- State is the copy number. CNery writes it as a bare integer in
# consensus mode and as a float under --polymorphism-mode, because the two decode over different
# grids and the polymorphism one keeps its column a float dtype even where a level happens to be
# whole. So a State of "1.0" IS the assertion that breseq passed the flag, and a State of "1" is the
# assertion that it did not.
states=`awk -F, 'NR > 1 { print $2 }' "${BREAKS}"`
[ -n "${states}" ] || fail "${BREAKS} has a header and no segments"

for state in ${states}; do
    case "${state}" in
        *.*) ;;
        *)   fail "CNery wrote copy number '${state}' as an integer in ${BREAKS};
       breseq did not pass --polymorphism-mode, or CNery ignored it" ;;
    esac
done

# ...and the other half of the round trip: an integral level must come back OUT of breseq without
# its decimal point, because copy_number has read "0" and "2" in every consensus .gd ever written
# and a polymorphism run has no business spelling the same value differently. This is the one place
# cn_copy_number_string()'s trailing-zero strip is checked end to end.
ANNOTATED="${SELF}/data/annotated.gd"
[ -f "${ANNOTATED}" ] || fail "no annotated .gd: ${ANNOTATED}"

bad=`awk -F'\t' '$1 == "CN" && $7 ~ /\./ && $7 ~ /\.0*$/ { print $7 }' "${ANNOTATED}"`
[ -z "${bad}" ] \
    || fail "CN entry in ${ANNOTATED} reports an integral copy number with a decimal point: ${bad}"

echo "Polymorphism-mode copy number checks passed."
