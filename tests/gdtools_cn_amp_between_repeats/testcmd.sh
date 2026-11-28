#!/bin/bash
#
# A copy-number gain between two copies of a repeat that are both in the reference is a tandem
# duplication (predictCNbetweenRepeatsToAMP), even with no junction: recombination between the
# copies makes none, since unit end -> repeat and repeat -> unit start both read as reference.
# Seen as a bare CN 1528001-1700900 between two IS6_292 copies of Thermincola CM166461, which
# dcasm confirmed with nanopore reads as AMP 1527539 175246 2.
#
# tests/data/cn_amp_repeats/cn_amp_repeats.gff3 is 60 kb of seeded random sequence with:
#   ISX  annotated repeat_region, 1500 bp, direct copies at 5001-6500 and 25001-26500
#   R    800 bp, NOT annotated, direct copies at 32001-32800 and 45001-45800
#   ISY  annotated repeat_region, 600 bp, at 50001-50600 and INVERTED at 56001-56600
#
# evidence.gd has CN gains, none with a junction. A gain has to cover EVERYTHING between the two
# copies; an edge may stop short of a copy only over windows CNery censored for redundant coverage,
# which the CN item records as redundant_before= / redundant_after= (bases), plus one tile.
#   CN 6801-24900   between the ISX copies (unique stretch 6501-25000), reaching 6501 and 25000
#                   over censored windows: AMP 6501 20000 2 between=ISX
#   CN 33001-44900  between the copies of R (unique stretch 32802-45000), likewise:
#                   AMP 32802 13000 2 between=repeat repeat_length=801
#                   (the base after each copy happens to agree, so the unit starts one base later)
#   CN 50701-55900  between INVERTED copies: recombination would invert, not duplicate. Left as CN.
#   CN 27001-31000  a repeat at one edge only: left as CN.
#   CN 6601-8000    near both ISX copies, but covers only 1.4 kb of the 18.5 kb between them, so it
#                   is not that unit's gain: left as CN.
#   CN 33501-44900  between the copies of R, but its first 700 bases are uncovered and NOT censored:
#                   left as CN.
#
# check.sh states those outcomes, so a golden rebuild cannot quietly drop one.

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/output.gd"
EXPECTED_OUTPUTS[0]="${SELF}/expected.gd"

TESTCMD="\
    ${GDTOOLS} MUTATIONS \
        -o ${SELF}/output.gd \
        -r ${DATADIR}/cn_amp_repeats/cn_amp_repeats.gff3 \
        ${SELF}/evidence.gd \
    && bash ${SELF}/check.sh ${SELF}/output.gd \
    "

do_test $1 ${SELF}
