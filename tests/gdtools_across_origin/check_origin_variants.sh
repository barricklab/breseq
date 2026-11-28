#!/bin/bash
#
# What gdtools_across_origin's origin_variants output must mean. See the origin_variants comment in
# testcmd.sh for the layout of each case.
#
# Usage: check_origin_variants.sh <output .gd>. Exits non-zero naming what went wrong.

GD="$1"

fail() { echo "CHECK FAILED (${GD}): $*"; exit 1; }

[ -f "${GD}" ] || fail "no output .gd"

# item <type> <side/seq fields...>: lines of that type whose fields 4.. start with the given values.
item() {
    local t="$1"; shift
    awk -F'\t' -v t="${t}" -v want="$*" '
        $1 != t { next }
        { n = split(want, w, " "); ok = 1; for (i = 1; i <= n; i++) if ($(i + 3) != w[i]) ok = 0; if (ok) print }' "${GD}"
}
count() { if [ -z "$1" ]; then echo 0; else echo "$1" | wc -l | tr -d ' '; fi }
has() { echo "$1" | tr '\t' '\n' | grep -qx "$2"; }

[ `grep -c 'ignore=CIRCULAR_CHROMOSOME' "${GD}"` -eq 3 ] || fail "expected 3 CIRCULAR_CHROMOSOME items"
[ `grep -c 'ignore=SEQUENCE_END' "${GD}"` -eq 1 ] || fail "expected 1 SEQUENCE_END item"

# DP 1: re-placed from the inverted copy to the copy that ends 30 bases short of the end.
dp1=`item DP origin-inverted 1 1 origin-inverted 2970 -1`
[ `count "${dp1}"` -eq 1 ] || fail "DP origin-inverted 1/1001 was not re-placed at 1 (+1) / 2970 (-1)"
has "${dp1}" "ignore=CIRCULAR_CHROMOSOME" || fail "DP origin-inverted 1/2970 is not ignored"
has "${dp1}" "side_2_placed_at_origin_from=origin-inverted:1001:1" || fail "DP origin-inverted 1/2970 lost where it was reported"

# DP 2: unique sequence, left alone.
dp2=`item DP origin-inverted 1 1 origin-inverted 1600 1`
[ `count "${dp2}"` -eq 1 ] || fail "DP origin-inverted 1/1600 moved"
echo "${dp2}" | grep -q 'ignore=' && fail "DP origin-inverted 1/1600 was ignored"

# DP 3: re-placed from another sequence, sides back in coordinate order.
dp3=`item DP plasmid-end 1 1 plasmid-end 1600 -1`
[ `count "${dp3}"` -eq 1 ] || fail "DP plasmid-carrier 1360/plasmid-end 1 was not re-placed at plasmid-end 1 (+1) / 1600 (-1)"
has "${dp3}" "ignore=CIRCULAR_CHROMOSOME" || fail "DP plasmid-end 1/1600 is not ignored"
has "${dp3}" "side_2_placed_at_origin_from=plasmid-carrier:1360:-1" || fail "DP plasmid-end 1/1600 lost where it was reported"
has "${dp3}" "side_2_redundant=1" || fail "DP plasmid-end 1/1600: the redundant flag did not move with its side"

# JC 4: inverted copy of the start, re-placed at the origin.
jc4=`item JC origin-jc-inverted 1 1 origin-jc-inverted 2000 -1`
[ `count "${jc4}"` -eq 1 ] || fail "JC origin-jc-inverted 1300/2000 was not re-placed at 1 (+1) / 2000 (-1)"
has "${jc4}" "ignore=CIRCULAR_CHROMOSOME" || fail "JC origin-jc-inverted 1/2000 is not ignored"

# MC 5 / MC 6 on the linear sequence.
mc5=`item MC linear-dtr 2501 3000`
has "${mc5}" "ignore=SEQUENCE_END" || fail "MC linear-dtr 2501-3000 is not ignored"
mc6=`item MC linear-dtr 1 60`
[ `count "${mc6}"` -eq 1 ] || fail "MC linear-dtr 1-60 is missing"
echo "${mc6}" | grep -q 'ignore=' && fail "MC linear-dtr 1-60 was ignored"

# Nothing here is a mutation except possibly a deletion from MC 6.
n=`grep -v '^#' "${GD}" | cut -f1 | grep -cE '^(SNP|SUB|INS|MOB|AMP|INV|CON|INT)$'`
[ "${n}" -eq 0 ] || fail "unexpected mutations"

echo "check_origin_variants.sh: all checks passed for ${GD}"
exit 0
