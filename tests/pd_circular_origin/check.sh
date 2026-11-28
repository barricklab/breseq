#!/bin/bash
#
# What pd_circular_origin's output must mean. See testcmd.sh.
#
# Usage: check.sh <annotated.gd> <sequence length>. Exits non-zero naming what went wrong.

GD="$1"
LENGTH="$2"
# Generous bound on the library's distance cutoff (450-bp fragments, sd 60): a PD this close to
# either end of the plasmid is inside the origin window.
WINDOW=1000

fail() { echo "CHECK FAILED: $*"; exit 1; }

[ -f "${GD}" ] || fail "no annotated .gd: ${GD}"

n=`awk -F'\t' '$1 ~ /^(SNP|SUB|DEL|INS|MOB|AMP|INV|CON|INT)$/' "${GD}" | wc -l | tr -d ' '`
[ "${n}" -eq 0 ] || fail "expected no mutations on the unmutated plasmid, found ${n}"

near=`awk -F'\t' -v L="${LENGTH}" -v w="${WINDOW}" '$1 == "PD" && ($5 <= w || $5 >= L - w)' "${GD}"`
[ -n "${near}" ] || fail "no PD near the origin: the case is not exercised"
bad=`echo "${near}" | grep -v 'ignore=CIRCULAR_CHROMOSOME'`
[ -z "${bad}" ] || fail "PD near the origin not ignored: ${bad}"

open=`awk -F'\t' '$1 == "PD" && !/reject=/ && !/ignore=/' "${GD}"`
[ -z "${open}" ] || fail "unrejected, unignored PD on the unmutated plasmid: ${open}"

# The -M1/-M2 pieces of a read split across the origin junction must carry the read's own strand on
# the reference: in this FR library, the opposite strand to its mate. write_moved_alignment used to
# reverse-complement such a piece without toggling 0x10, so every pair here read as FF/RR.
BAM="`dirname "${GD}"`/reference.bam"
if [ -f "${BAM}" ] && command -v samtools > /dev/null 2>&1; then
    counts=`samtools view "${BAM}" | awk '{
        rev = int($2 / 16) % 2
        if ($1 ~ /-M[12]$/) { b = $1; sub(/-M[12]$/, "", b); split(b, p, ":"); m = ((p[1] == 1) ? 2 : 1) ":" p[2]; piece[m] = rev }
        else mate[$1] = rev
      } END { for (m in piece) if (m in mate) { tot++; if (piece[m] != mate[m]) opp++ } print tot + 0, opp + 0 }'`
    set -- ${counts}
    [ "$1" -gt 0 ] || fail "no origin-crossing -M pieces with a mate in ${BAM}"
    [ "$1" -eq "$2" ] || fail "$(( $1 - $2 )) of $1 origin -M pieces are on the same strand as their mate"
fi

echo "check.sh: PD circular-origin checks passed"
exit 0
