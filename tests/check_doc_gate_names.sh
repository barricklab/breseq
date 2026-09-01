#!/bin/bash
#
# Verify that the metric and gate row names emitted into summary.html match the "###" headings
# that document them, one for one and in the same order.
#
# The four experimental evidence types (DP, PD, SC, MP) each print a metrics table and a gates
# table in summary.html, and each row of those tables has a matching section on the type's
# documentation page. Renaming a row in output.cpp without renaming the heading -- or adding a
# gate and forgetting to document it -- silently leaves the manual describing a report that no
# longer exists. That has to be caught mechanically, because both sides look fine on their own.
#
# Reads only; never builds or runs anything. Exits non-zero on a mismatch.

set -u

SELF=$(dirname "${BASH_SOURCE[0]}")
ROOTDIR="${SELF}/.."

exec python3 - "$ROOTDIR" <<'PYEOF'
import re, sys, os

root = sys.argv[1]
src = open(os.path.join(root, 'src/breseq/output.cpp')).read()

FUNCS = [
    ('evidence-dp.md', 'html_discordant_pair_gates_string'),
    ('evidence-pd.md', 'html_pair_distance_gates_string'),
    ('evidence-sc.md', 'html_soft_clipping_gates_string'),
    ('evidence-mp.md', 'html_missing_pair_gates_string'),
]

fail = 0
for doc, fn in FUNCS:
    i = src.index('string ' + fn)
    j = src.index('\n}\n', i)
    body = src[i:j]

    # Row names are the first cell of each emitted row: tr(td("<name>") ...
    rows = re.findall(r'tr\(td\("([a-z][a-z0-9 \-]*)"\)', body)

    path = os.path.join(root, 'docs', doc)
    heads = re.findall(r'^### ([a-z].*)$', open(path).read(), re.M)

    missing = [r for r in rows if r not in heads]
    extra = [h for h in heads if h not in rows]

    if missing or extra:
        fail = 1
        print("MISMATCH: docs/%s vs %s()" % (doc, fn))
        for r in missing:
            print("    row emitted in summary.html but not documented: '%s'" % r)
        for h in extra:
            print("    documented but no longer emitted:               '%s'" % h)
    else:
        print("  docs/%-20s %2d rows match" % (doc, len(rows)))

if fail:
    print("check_doc_gate_names: FAILED -- see above")
else:
    print("check_doc_gate_names: metric and gate names match their documentation")
sys.exit(fail)
PYEOF
