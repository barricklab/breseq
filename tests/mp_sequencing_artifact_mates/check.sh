#!/bin/bash
#
# What mp_sequencing_artifact_mates's two runs must mean. See testcmd.sh.
#
# Usage: check.sh <baseline output.gd> <gated annotated.gd> <gated reference.bam>.

BASELINE="$1"
GATED="$2"
BAM="$3"

fail() { echo "CHECK FAILED: $*"; exit 1; }

for f in "${BASELINE}" "${GATED}" "${BAM}"; do [ -f "${f}" ] || fail "missing ${f}"; done

# Unrejected MP items within 400 bases of a position.
mp_near() {
    awk -F'\t' -v p="$2" '$1 == "MP" && !/reject=/ && !/ignore=/ && ($5 > p - 400) && ($5 < p + 400)' "$1"
}

[ -n "`mp_near "${BASELINE}" 6000`" ] || fail "baseline (gate off) has no MP at the artifact site: the case is not exercised"
[ -n "`mp_near "${BASELINE}" 14000`" ] || fail "baseline has no MP at the real insertion"

[ -z "`mp_near "${GATED}" 6000`" ] || fail "gated run still reports an MP at the artifact site: `mp_near "${GATED}" 6000`"
[ -n "`mp_near "${GATED}" 14000`" ] || fail "gated run lost the MP at the real insertion"

# The reads whose mates were set aside carry X9:i:2 (MP counts only X9:i:1).
if command -v samtools > /dev/null 2>&1; then
    n=`samtools view "${BAM}" | grep -c 'X9:i:2'`
    [ "${n}" -ge 30 ] || fail "only ${n} reads with an artifact mate (X9:i:2), expected most of the 60 planted"
    echo "check.sh: MP artifact-mate checks passed (${n} reads with an artifact mate)"
else
    echo "check.sh: MP artifact-mate checks passed (samtools not found: X9 count not checked)"
fi
exit 0
