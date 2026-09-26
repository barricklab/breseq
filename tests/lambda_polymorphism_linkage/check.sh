#!/bin/bash
# Checks on the annotated.gd of the simulated mixture (see testcmd.sh). Portable awk only:
# the harness runs under whatever grep the platform ships.
#
# Usage: check.sh <annotated.gd>
# Exits non-zero, naming the failed check, if what the test is for is not in the file.

GD="$1"

fail() { echo "CHECK FAILED: $1"; exit 1; }

# The 30% genotype carries INS AC at 10000 and the 20% genotype INS A at the same site. Linkage
# resolves them as two haplotypes: one INS AC and one INS A, each at its own frequency -- not one
# INS A at the summed 50% plus a nested INS C at insert position 2.
n=$(awk -F'\t' '$1=="INS" && $4=="NC_001416" && $5=="10000" && $6=="AC"' "$GD" | wc -l)
[ "$n" -eq 1 ] || fail "expected exactly one INS AC at 10000, found $n"
n=$(awk -F'\t' '$1=="INS" && $4=="NC_001416" && $5=="10000" && $6=="A"' "$GD" | wc -l)
[ "$n" -eq 1 ] || fail "expected exactly one INS A at 10000, found $n"
n=$(awk -F'\t' '$1=="INS" && $4=="NC_001416" && $5=="10000" && $6=="C"' "$GD" | wc -l)
[ "$n" -eq 0 ] || fail "found $n nested INS C at 10000"
f=$(awk -F'\t' '$1=="INS" && $4=="NC_001416" && $5=="10000" && $6=="AC" { for (i = 6; i <= NF; i++) if ($i ~ /^frequency=/) { split($i, a, "="); print a[2] + 0 } }' "$GD")
awk -v f="$f" 'BEGIN { exit !(f >= 0.20 && f <= 0.40) }' || fail "INS AC frequency $f outside [0.20, 0.40]"
f=$(awk -F'\t' '$1=="INS" && $4=="NC_001416" && $5=="10000" && $6=="A" { for (i = 6; i <= NF; i++) if ($i ~ /^frequency=/) { split($i, a, "="); print a[2] + 0 } }' "$GD")
awk -v f="$f" 'BEGIN { exit !(f >= 0.10 && f <= 0.30) }' || fail "INS A frequency $f outside [0.10, 0.30]"

# One SUB T>CA at 13020, and no leftover SNP or INS from calling its columns separately.
n=$(awk -F'\t' '$1=="SUB" && $4=="NC_001416" && $5=="13020" && $6=="1" && $7=="CA"' "$GD" | wc -l)
[ "$n" -eq 1 ] || fail "expected exactly one SUB CA at 13020, found $n"
n=$(awk -F'\t' '($1=="SNP" || $1=="INS") && $4=="NC_001416" && ($5=="13019" || $5=="13020")' "$GD" | wc -l)
[ "$n" -eq 0 ] || fail "found $n separate SNP/INS calls at 13019-13020"

# Every simulated mutation is called once, at a frequency consistent with its 30% (or 20%)
# genotype given read sampling at ~35 spanning reads.
awk -F'\t' '
  BEGIN { want["11702"]; want["12000"]; want["12005"]; want["12011"]; want["13020"]; want["15000"]; want["16000"]; want["17246"]; want["22375"] }
  ($1=="INS" || $1=="SUB" || $1=="DEL" || $1=="SNP") && $4=="NC_001416" && ($5 in want) {
    seen[$5]++
    # The first key=value column differs by type (an INS has fewer fixed columns than a SUB).
    f = -1
    for (i = 6; i <= NF; i++) if ($i ~ /^frequency=/) { split($i, a, "="); f = a[2] + 0 }
    if (f < 0.15 || f > 0.45) { print "CHECK FAILED: frequency " f " at " $5 " outside [0.15, 0.45]"; bad = 1 }
  }
  END {
    for (p in want) if (seen[p] != 1) { print "CHECK FAILED: expected exactly one mutation at " p ", found " seen[p]+0; bad = 1 }
    exit bad
  }' "$GD" || exit 1

# The two cis SNPs and the trans pair are reported as such.
n=$(awk -F'\t' '$1=="LN" { for (i = 8; i <= NF; i++) if ($i == "phase=cis") c++ } END { print c+0 }' "$GD")
[ "$n" -ge 1 ] || fail "expected a phase=cis LN"
n=$(awk -F'\t' '$1=="LN" { for (i = 8; i <= NF; i++) if ($i == "phase=trans") c++ } END { print c+0 }' "$GD")
[ "$n" -ge 1 ] || fail "expected a phase=trans LN"

exit 0
