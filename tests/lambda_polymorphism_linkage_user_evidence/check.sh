#!/bin/bash
# Checks on the annotated.gd of the simulated mixture run with --user-evidence-gd naming the A
# insertion at 10000 (see testcmd.sh). Portable awk only.
#
# Usage: check.sh <annotated.gd>

GD="$1"

fail() { echo "CHECK FAILED: $1"; exit 1; }

# The user asked about INS A at 10000. It is reported on its own, as user evidence always is, at the
# COLUMN's frequency (A is the first inserted base of both the A-only and the AC lineage, so about
# 50%), and the linked haplotypes are reported too: INS AC at about 30%. The haplotype INS A would
# duplicate the user's mutation and is not emitted a second time.
n=$(awk -F'\t' '$1=="INS" && $4=="NC_001416" && $5=="10000" && $6=="A"' "$GD" | wc -l)
[ "$n" -eq 1 ] || fail "expected exactly one INS A at 10000, found $n"
n=$(awk -F'\t' '$1=="INS" && $4=="NC_001416" && $5=="10000" && $6=="AC"' "$GD" | wc -l)
[ "$n" -eq 1 ] || fail "expected exactly one INS AC at 10000, found $n"
n=$(awk -F'\t' '$1=="INS" && $4=="NC_001416" && $5=="10000" && $6=="C"' "$GD" | wc -l)
[ "$n" -eq 0 ] || fail "found $n nested INS C at 10000"

f=$(awk -F'\t' '$1=="INS" && $4=="NC_001416" && $5=="10000" && $6=="A" { for (i = 6; i <= NF; i++) if ($i ~ /^frequency=/) { split($i, a, "="); print a[2] + 0 } }' "$GD")
awk -v f="$f" 'BEGIN { exit !(f >= 0.35 && f <= 0.65) }' || fail "user INS A frequency $f outside [0.35, 0.65]"
f=$(awk -F'\t' '$1=="INS" && $4=="NC_001416" && $5=="10000" && $6=="AC" { for (i = 6; i <= NF; i++) if ($i ~ /^frequency=/) { split($i, a, "="); print a[2] + 0 } }' "$GD")
awk -v f="$f" 'BEGIN { exit !(f >= 0.20 && f <= 0.40) }' || fail "INS AC frequency $f outside [0.20, 0.40]"

# The user's RA is marked as such.
n=$(awk -F'\t' '$1=="RA" && $4=="NC_001416" && $5=="10000" && $6=="1" { for (i = 9; i <= NF; i++) if ($i == "user_defined=1") c++ } END { print c+0 }' "$GD")
[ "$n" -ge 1 ] || fail "expected the RA at 10000.1 to be user_defined"

exit 0
