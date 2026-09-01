#!/bin/bash
#
# Verify that every GenomeDiff example quoted in docs/ actually exists in the golden it names.
#
# The evidence documentation (docs/evidence-*.md) quotes real .gd lines from the long tests'
# committed expected.gd files, inside fenced blocks tagged with the file they came from:
#
#     ```text title="tests/long_ltee_ara_m1_40k_pe36/expected.gd"
#     MP	1951	.	REL606	2039143	1	score=13.0	...
#     ```
#
# Long lines are wrapped for display and fields are trimmed to those under discussion, so we cannot
# compare the block verbatim. Instead we take the identifying prefix of each quoted entry -- the
# type, the entry id, and the first two positional fields -- and require a line in the named golden
# that starts with exactly those. That is enough to catch an example that was invented, mistyped, or
# left behind by a golden rebuild, which is the failure this guards against.
#
# Exits non-zero if any example cannot be found. Reads only; never re-runs a test.

set -u

SELF=$(dirname "${BASH_SOURCE[0]}")
DOCSDIR="${SELF}/../docs"
ROOTDIR="${SELF}/.."

fail=0
checked=0

# Evidence types whose .gd lines we know how to identify.
TYPES='RA|MC|JC|CN|UN|SC|DP|MP|PD|SNP|SUB|DEL|INS|MOB|AMP|INV|CON|INT'

for doc in "${DOCSDIR}"/*.md; do
    [ -e "$doc" ] || continue

    # Walk the file, tracking which golden the current fenced block is tagged with.
    golden=""
    while IFS= read -r line; do
        case "$line" in
            '```'*title=*)
                golden=$(printf '%s\n' "$line" | sed -n 's/.*title="\([^"]*\)".*/\1/p')
                continue
                ;;
            '```'*)
                golden=""
                continue
                ;;
        esac

        [ -n "$golden" ] || continue

        # Only look at lines that open a GenomeDiff entry: TYPE <tab> id <tab> ...
        printf '%s\n' "$line" | grep -qE "^(${TYPES})"$'\t' || continue

        if [ ! -e "${ROOTDIR}/${golden}" ]; then
            echo "MISSING GOLDEN: ${doc}: names ${golden}, which does not exist"
            fail=1
            continue
        fi

        # Identify the entry by type, id, and the first two positional fields. Field 3 is
        # parent-ids, which is '.' for evidence; fields 4 and 5 are seq_id and position/start.
        key=$(printf '%s\n' "$line" | cut -f1,2,4,5)
        if ! cut -f1,2,4,5 "${ROOTDIR}/${golden}" | grep -qxF "$key"; then
            echo "NOT FOUND: ${doc}"
            echo "    golden: ${golden}"
            echo "    entry:  $(printf '%s' "$key" | tr '\t' ' ')"
            fail=1
        fi
        checked=$((checked + 1))
    done < "$doc"
done

if [ "$fail" -eq 0 ]; then
    echo "check_doc_examples: ${checked} quoted GenomeDiff example(s) verified against their goldens"
else
    echo "check_doc_examples: FAILED -- see above"
fi

exit "$fail"
