#!/bin/bash
#
# What lambda_polymorphism_del_amp's golden must mean, so that a rebuild after a change to the
# simulator or to a default cannot silently lose one of the four calls. See testcmd.sh.
#
# Usage: check.sh <annotated.gd>. Exits non-zero naming what went wrong.

GD="$1"

fail() { echo "CHECK FAILED: $*"; exit 1; }

[ -f "${GD}" ] || fail "no annotated .gd: ${GD}"

# field <type> <position> <size> <key>: the value of key= on the one mutation of that type,
# position and size, or nothing.
field() {
    awk -F'\t' -v t="$1" -v p="$2" -v s="$3" -v k="$4" '
        $1 == t && $5 == p && $6 == s {
            for (i = 7; i <= NF; i++) { split($i, kv, "="); if (kv[1] == k) print kv[2] }
        }' "${GD}"
}

# near <value> <target> <tolerance>
near() {
    awk -v v="$1" -v t="$2" -v tol="$3" 'BEGIN { d = v - t; if (d < 0) d = -d; exit !(v != "" && d <= tol) }'
}

# Exactly the four simulated mutations, and nothing else.
n=`awk -F'\t' '$1 ~ /^(SNP|SUB|DEL|INS|MOB|AMP|INV|CON|INT)$/' "${GD}" | wc -l | tr -d ' '`
[ "${n}" -eq 4 ] || fail "expected 4 mutations, found ${n}"

# The 3 kb deletion: junction and coverage both, combined.
f=`field DEL 10001 3000 frequency`
near "${f}" 0.5 0.08 || fail "DEL 10001 3000 frequency '${f}' is not ~0.5"
[ -n "`field DEL 10001 3000 jc_frequency`" ] || fail "DEL 10001 3000 has no jc_frequency"
[ -n "`field DEL 10001 3000 cn_frequency`" ] || fail "DEL 10001 3000 has no cn_frequency"

# The 80 bp deletion: below CN's window size, so the junction alone.
f=`field DEL 30001 80 frequency`
near "${f}" 0.5 0.08 || fail "DEL 30001 80 frequency '${f}' is not ~0.5"
[ -n "`field DEL 30001 80 jc_frequency`" ] || fail "DEL 30001 80 has no jc_frequency"
[ -z "`field DEL 30001 80 cn_frequency`" ] || fail "DEL 30001 80 has a cn_frequency; CN should not resolve 80 bp"

# The amplifications: copy number 2 and 3, each with a frequency and both of its sources.
for amp in "20001 5000 2" "35001 3000 3"; do
    set -- ${amp}
    copies=`awk -F'\t' -v p="$1" -v s="$2" '$1 == "AMP" && $5 == p && $6 == s { print $7 }' "${GD}"`
    [ "${copies}" = "$3" ] || fail "AMP $1 $2 has new_copy_number '${copies}', expected $3"
    f=`field AMP $1 $2 frequency`
    near "${f}" 0.75 0.08 || fail "AMP $1 $2 frequency '${f}' is not ~0.75"
    [ -n "`field AMP $1 $2 jc_frequency`" ] || fail "AMP $1 $2 has no jc_frequency"
    [ -n "`field AMP $1 $2 cn_frequency`" ] || fail "AMP $1 $2 has no cn_frequency"
done

echo "Polymorphism deletion and amplification checks passed."
