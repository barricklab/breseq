#!/bin/bash
#
# What gene_conversion_1's golden must mean, so that a rebuild after a change to the simulator or
# to a default cannot silently lose one of its four cases. See testcmd.sh for cases (A)-(D).
#
# Usage: check.sh <annotated.gd>. Exits non-zero naming what went wrong.

GD="$1"

fail() { echo "CHECK FAILED: $*"; exit 1; }

[ -f "${GD}" ] || fail "no annotated .gd: ${GD}"

# count <awk condition>: how many lines satisfy it.
count() { awk -F'\t' "$1 { n++ } END { print n + 0 }" "${GD}"; }

# Exactly four mutations: two conversions, two deletions.
n=`count '$1 ~ /^(SNP|SUB|DEL|INS|MOB|AMP|INV|CON|INT)$/'`
[ "${n}" -eq 4 ] || fail "expected 4 mutations, found ${n}"
n=`count '$1 == "CON"'`
[ "${n}" -eq 2 ] || fail "expected 2 CON mutations, found ${n}"

# con <lo> <hi>: "position size region" of the one CON starting in lo..hi.
con() { awk -F'\t' -v lo="$1" -v hi="$2" '$1 == "CON" && $5 >= lo && $5 <= hi { print $5, $6, $7 }' "${GD}"; }

# (A) Forward donor. The tract must contain the converted 4401-4800 and stop inside the paralog,
#     and the donor is the same stretch 8000 bp downstream, written forwards.
set -- `con 4001 4401`
[ -n "$1" ] || fail "(A) no CON at pair A (converted 4401-4800)"
end=$(( $1 + $2 - 1 ))
[ "${end}" -ge 4800 -a "${end}" -le 5200 ] || fail "(A) CON $1 $2 does not cover 4401-4800 within the paralog"
[ "$3" = "NC_001416:$(( $1 + 8000 ))-$(( end + 8000 ))" ] \
    || fail "(A) CON $1 $2 has donor region '$3', expected the same stretch 8000 bp downstream"

# (B) Reverse-strand donor. Offset o of the block lies at 25200 - o in the antiparallel copy, so
#     the region must run backwards from 41201 - start to 41201 - end.
set -- `con 16001 16401`
[ -n "$1" ] || fail "(B) no CON at pair B (converted 16401-16800)"
end=$(( $1 + $2 - 1 ))
[ "${end}" -ge 16800 -a "${end}" -le 17200 ] || fail "(B) CON $1 $2 does not cover 16401-16800 within the paralog"
[ "$3" = "NC_001416:$(( 41201 - $1 ))-$(( 41201 - end ))" ] \
    || fail "(B) CON $1 $2 has donor region '$3', expected the antiparallel copy read backwards"

# (C) A deleted copy is a DEL, not a CON. The deletion's ends are identical for a base or two past
#     the block in some registers, so the position may shift by a few bases; the size may not.
n=`count '$1 == "DEL" && $5 >= 27996 && $5 <= 28006 && $6 == 1200'`
[ "${n}" -eq 1 ] || fail "(C) the deleted paralog copy (28001, 1200 bp) was not predicted as a DEL"

# (D) The deletion between two unannotated paralogs, to the base.
n=`count '$1 == "DEL" && $5 == 37582 && $6 == 8000'`
[ "${n}" -eq 1 ] || fail "(D) DEL 37582 8000 between the pair D copies was not predicted exactly"

echo "Gene conversion and paralog deletion checks passed."
