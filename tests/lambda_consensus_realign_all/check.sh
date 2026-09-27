#!/bin/bash
# Checks on the annotated.gd of the simulated clone (see testcmd.sh). Portable awk only: the
# harness runs under whatever grep the platform ships.
#
# Usage: check.sh <annotated.gd>
# Exits non-zero, naming the failed check, if what the test is for is not in the file.

GD="$1"

fail() { echo "CHECK FAILED: $1"; exit 1; }

# Every simulated mutation is called exactly once, at the right-normalized coordinate breseq
# reports, with the right sequence, and as a FIXED mutation: consensus-mode mutations carry no
# frequency key, whether they came from the per-column join or from an LN.
awk -F'\t' '
  BEGIN {
    want["INS 10000 AC"]; want["INS 11702 CAG"]; want["SNP 12000 G"]; want["SNP 12011 C"];
    want["SUB 13020 1 CA"]; want["SNP 15000 A"]; want["DEL 17246 3"]; want["DEL 22375 1"]
  }
  ($1=="INS" || $1=="SNP" || $1=="DEL") && $4=="NC_001416" { key = $1 " " $5 " " $6; n_first = 7 }
  $1=="SUB" && $4=="NC_001416" { key = $1 " " $5 " " $6 " " $7; n_first = 8 }
  ($1=="INS" || $1=="SNP" || $1=="DEL" || $1=="SUB") && $4=="NC_001416" {
    total++
    if (key in want) seen[key]++
    else { print "CHECK FAILED: unexpected mutation " key; bad = 1 }
    for (i = n_first; i <= NF; i++) if ($i ~ /^frequency=/) { print "CHECK FAILED: consensus-mode mutation carries " $i " at " $5; bad = 1 }
  }
  END {
    for (k in want) if (seen[k] != 1) { print "CHECK FAILED: expected exactly one " k ", found " seen[k]+0; bad = 1 }
    exit bad
  }' "$GD" || exit 1

# The multi-column mutations were taken over by a contiguous LN whose refined (realigned=1) fit
# calls the all-variant haplotype consensus: INS AC at 10000, the T>CA SUB at 13019-13020, and
# the 3-base repeat-unit deletion at 17246. (INS CAG at 11702 is a run too, but its columns all
# sit at one position, so it is identified by its start alone.)
for start in 10000 11702 13019 17246; do
  n=$(awk -F'\t' -v s="$start" '$1=="LN" && $4=="NC_001416" && $5==s {
      c = 0; r = 0; p = 0
      for (i = 8; i <= NF; i++) { if ($i == "contiguous=1") c = 1; if ($i == "realigned=1") r = 1; if ($i ~ /^haplotype_predictions=/ && $i ~ /consensus/) p = 1 }
      if (c && r && p) n++
    } END { print n+0 }' "$GD")
  [ "$n" -ge 1 ] || fail "expected a contiguous, realigned LN starting at $start with a consensus haplotype"
done

# Under --realign all every RA column is re-scored, lone SNPs included. A column inside a linked
# run is re-scored as part of the run, and the run's LN carries realigned=1 for it; every other
# column must carry realigned=1 itself.
awk -F'\t' '
  $1=="LN" { c = 0; for (i = 8; i <= NF; i++) if ($i == "contiguous=1") c = 1; if (c) { n_ln++; ln_start[n_ln] = $5 + 0; ln_end[n_ln] = $7 + 0 } }
  $1=="RA" { r = 0; for (i = 9; i <= NF; i++) if ($i == "realigned=1") r = 1; if (!r) { n_ra++; ra_pos[n_ra] = $5 + 0 } }
  END {
    for (k = 1; k <= n_ra; k++) {
      inside = 0
      for (l = 1; l <= n_ln; l++) if ((ra_pos[k] >= ln_start[l]) && (ra_pos[k] <= ln_end[l])) inside = 1
      if (!inside) { print "CHECK FAILED: RA column at " ra_pos[k] " was not realigned and is in no linked run"; bad = 1 }
    }
    exit bad
  }' "$GD" || exit 1

exit 0
