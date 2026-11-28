#!/bin/bash
#
# What gdtools_apply_feature_copy_strand's output must mean. See testcmd.sh.
#
# Usage: check.sh <output.gff3>. Exits non-zero naming what went wrong.

GFF="$1"

fail() { echo "CHECK FAILED: $*"; exit 1; }

[ -f "${GFF}" ] || fail "no output GFF3: ${GFF}"

# cds <gene>: "start end strand" for every CDS of that gene, sorted by start.
cds() {
    awk -F'\t' -v g="$1" '
        $3 == "CDS" && $9 ~ ("(^|;)(gene|Name)=" g "(;|$)") { print $4, $5, $7 }' "${GFF}" \
    | sort -n -k1,1 | uniq
}

# The AMP unit is 4001-6300 and is inserted before itself, so the first copy keeps the original
# coordinates and the second is shifted by 2300.
mltd=`cds mltD`
glob=`cds gloB`

echo "${glob}" | grep -qx "5461 6216 -" || fail "gloB first AMP copy is not complement(5461..6216): ${glob}"
echo "${glob}" | grep -qx "7761 8516 -" || fail "gloB second AMP copy is not complement(7761..8516): ${glob}"
[ `echo "${glob}" | wc -l | tr -d ' '` -eq 2 ] || fail "expected 2 gloB CDS: ${glob}"

echo "${mltd}" | grep -qx "4031 5389 -" || fail "mltD first AMP copy is not complement(4031..5389): ${mltd}"
echo "${mltd}" | grep -qx "6331 7689 -" || fail "mltD second AMP copy is not complement(6331..7689): ${mltd}"

# Two more mltD copies from the INTs: the forward donor keeps the minus strand, the reverse donor
# (an inverted copy) is plus. Each spans the gene's length, 1359 bases.
[ `echo "${mltd}" | wc -l | tr -d ' '` -eq 4 ] || fail "expected 4 mltD CDS: ${mltd}"
[ `echo "${mltd}" | awk '$3 == "-"' | wc -l | tr -d ' '` -eq 3 ] || fail "expected 3 minus-strand mltD: ${mltd}"
[ `echo "${mltd}" | awk '$3 == "+"' | wc -l | tr -d ' '` -eq 1 ] || fail "expected 1 plus-strand mltD: ${mltd}"
[ `echo "${mltd}" | awk '$2 - $1 + 1 != 1359' | wc -l | tr -d ' '` -eq 0 ] || fail "an mltD copy changed length: ${mltd}"

# Within its 1500-base donor the gene starts 30 bases in when copied forward (4031 - 4001) and 111
# bases in when the donor is reversed (5500 - 5389), and the two INTs are inserted 5000 reference
# bases apart plus the first INT's 1500 bases. So the plus copy starts 5000 + 1500 + 111 - 30 = 6581
# bases after the forward one.
fwd_start=`echo "${mltd}" | awk '$1 > 9000 && $3 == "-" { print $1 }'`
rev_start=`echo "${mltd}" | awk '$3 == "+" { print $1 }'`
[ -n "${fwd_start}" ] && [ -n "${rev_start}" ] || fail "missing an INT copy of mltD: ${mltd}"
[ $(( rev_start - fwd_start )) -eq 6581 ] || fail "inverted INT copy of mltD is not mirrored in its donor: ${mltd}"

echo "check.sh: all feature-copy strand checks passed"
exit 0
