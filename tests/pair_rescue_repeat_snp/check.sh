#!/bin/bash
#
# What pair_rescue_repeat_snp's two runs must mean. See testcmd.sh.
#
# Usage: check.sh <baseline output.gd> <rescue annotated.gd>. Exits non-zero naming what went wrong.

BASELINE="$1"
RESCUE="$2"

fail() { echo "CHECK FAILED: $*"; exit 1; }

[ -f "${BASELINE}" ] || fail "no baseline .gd: ${BASELINE}"
[ -f "${RESCUE}" ] || fail "no rescue .gd: ${RESCUE}"

# Unrejected DPs with one side in each repeat copy (allowing a read length or so of placement slop).
between_copies() {
    awk -F'\t' '$1 == "DP" && !/reject=/ {
        p1 = $5; p2 = $8
        if (p1 > 4800 && p1 < 6200 && p2 > 14800 && p2 < 16200) print
    }' "$1"
}
mutations() { awk -F'\t' '$1 ~ /^(SNP|SUB|DEL|INS|MOB|AMP|INV|CON|INT)$/' "$1"; }

b=`between_copies "${BASELINE}"`
[ -n "${b}" ] || fail "baseline (rescue off) predicts no DP between the copies: the case is not exercised"
echo "${b}" | tr '\t' '\n' | grep -qx "discordant_count=30" || fail "baseline DP does not have the 30 planted pairs: ${b}"

r=`between_copies "${RESCUE}"`
[ -z "${r}" ] || fail "rescue run still predicts a DP between the copies: ${r}"

[ -z "`mutations "${BASELINE}"`" ] || fail "baseline predicts a mutation"
[ -z "`mutations "${RESCUE}"`" ] || fail "rescue run predicts a mutation (a rescued mate's q10 base must not call one)"

echo "check.sh: pair rescue checks passed"
exit 0
