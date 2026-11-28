#!/bin/bash
#
# What gdtools_cn_amp_between_repeats's output must mean. See testcmd.sh.
#
# Usage: check.sh <output.gd>. Exits non-zero naming what went wrong.

GD="$1"

fail() { echo "CHECK FAILED: $*"; exit 1; }

[ -f "${GD}" ] || fail "no output .gd: ${GD}"

# amp <position> <size>: the AMP line at that position and size, or nothing.
amp() { awk -F'\t' -v p="$1" -v s="$2" '$1 == "AMP" && $5 == p && $6 == s' "${GD}"; }
has() { echo "$1" | tr '\t' '\n' | grep -qx "$2"; }

n=`awk -F'\t' '$1 ~ /^(SNP|SUB|DEL|INS|MOB|AMP|INV|CON|INT)$/' "${GD}" | wc -l | tr -d ' '`
[ "${n}" -eq 2 ] || fail "expected 2 mutations, found ${n}"

a1=`amp 6501 20000`
[ -n "${a1}" ] || fail "no AMP 6501 20000 between the ISX copies"
echo "${a1}" | awk -F'\t' '$7 == 2 { ok = 1 } END { exit !ok }' || fail "AMP 6501 20000 is not 2 copies"
has "${a1}" "between=ISX" || fail "AMP 6501 20000 is not between=ISX"

a2=`amp 32802 13000`
[ -n "${a2}" ] || fail "no AMP 32802 13000 between the unannotated copies"
echo "${a2}" | awk -F'\t' '$7 == 2 { ok = 1 } END { exit !ok }' || fail "AMP 32802 13000 is not 2 copies"
has "${a2}" "between=repeat" || fail "AMP 32802 13000 is not between=repeat"
has "${a2}" "repeat_length=801" || fail "AMP 32802 13000 does not have repeat_length=801"

echo "check.sh: all CN-between-repeats checks passed"
exit 0
