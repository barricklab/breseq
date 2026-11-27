#!/bin/bash
#
# What lambda_pair_distances' golden must mean, so that a rebuild after a change to the simulator
# or to a default cannot silently lose the PD or the MP control. See testcmd.sh.
#
# Usage: check.sh <annotated.gd>. Exits non-zero naming what went wrong.

GD="$1"

fail() { echo "CHECK FAILED: $*"; exit 1; }

[ -f "${GD}" ] || fail "no annotated .gd: ${GD}"

# count <awk condition>: how many lines satisfy it.
count() { awk -F'\t' "$1 { n++ } END { print n + 0 }" "${GD}"; }

# The deletion is the only event with a junction, so it is the only mutation.
n=`count '$1 ~ /^(SNP|SUB|DEL|INS|MOB|AMP|INV|CON|INT)$/'`
[ "${n}" -eq 1 ] || fail "expected 1 mutation, found ${n}"
n=`count '$1 == "DEL" && $5 == 12000 && $6 == 200'`
[ "${n}" -eq 1 ] || fail "DEL 12000 200 was not predicted"

# PD: one accepted item at each 200 bp event, and nowhere else.
n=`count '$1 == "PD" && $0 !~ /reject=/'`
[ "${n}" -eq 2 ] || fail "expected 2 accepted PD items, found ${n}"
# (PD places a breakpoint from the pairs, to within a few bases; the JC is what makes it exact.)
n=`count '$1 == "PD" && $0 !~ /reject=/ && $5 >= 11989 && $5 <= 12009 && $8 >= 12190 && $8 <= 12210'`
[ "${n}" -eq 1 ] || fail "no accepted PD at the 200 bp deletion (11999 / 12200)"
n=`count '$1 == "PD" && $0 !~ /reject=/ && $5 >= 23980 && $5 <= 24020'`
[ "${n}" -eq 1 ] || fail "no accepted PD at the 200 bp insertion (24000)"

# MP: accepted on both flanks of the 3 kb insertion, one facing each way. The items at the 200 bp
# insertion score at the cutoff and are deliberately not asserted either way (see testcmd.sh), but
# nothing may be accepted anywhere else.
n=`count '$1 == "MP" && $0 !~ /reject=/ && $5 >= 35980 && $5 <= 36020 && $6 == -1'`
[ "${n}" -eq 1 ] || fail "expected 1 accepted MP with strand -1 at the 3 kb insertion (36000), found ${n}"
n=`count '$1 == "MP" && $0 !~ /reject=/ && $5 >= 35980 && $5 <= 36020 && $6 == 1'`
[ "${n}" -eq 1 ] || fail "expected 1 accepted MP with strand 1 at the 3 kb insertion (36000), found ${n}"
n=`count '$1 == "MP" && $0 !~ /reject=/ && !($5 >= 35980 && $5 <= 36020) && !($5 >= 23800 && $5 <= 24200)'`
[ "${n}" -eq 0 ] || fail "${n} accepted MP item(s) away from both insertions"

# The two ends of this LINEAR reference are joined in the reads. Both the junction and the
# discordant pairs must say so, visibly: neither rejected nor ignored (CIRCULAR_CHROMOSOME is for
# sequences declared circular), and no mutation made out of them.
n=`count '$1 == "JC" && $5 == 1 && $6 == 1 && $8 == 48502 && $9 == -1 && $0 !~ /reject=/ && $0 !~ /ignore=/'`
[ "${n}" -eq 1 ] || fail "no visible JC joining the two ends of the linear reference (1 / 48502)"
n=`count '$1 == "DP" && $5 == 1 && $6 == 1 && $8 == 48502 && $9 == -1 && $0 !~ /reject=/ && $0 !~ /ignore=/'`
[ "${n}" -eq 1 ] || fail "no visible DP joining the two ends of the linear reference (1 / 48502)"

echo "Pair distance, missing pair and joined-ends checks passed."
