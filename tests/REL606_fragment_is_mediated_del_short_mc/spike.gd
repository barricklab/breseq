#=GENOME_DIFF	1.0
#
# Stand-in for an IS-element transposition intermediate. The excised circle carries the element's
# target-site duplication, so its reads match the chromosome for TSD-many bases PAST the annotated
# element boundary -- which is how a real MC boundary ends up stopping short of the element that
# mediated the deletion.
#
# This keeps reference 31293..32078: the IS1 copy at 31302..32069 together with both copies of its
# real 9 bp TSD (CCGCGACAA at 31293..31301 and 32070..32078). Every base of it is contiguous
# reference sequence, so it introduces no junction of its own. Its only effect is to hold coverage
# over 31293..31301, the 9 bases between DEL 28183 3119's true endpoint and the IS1 boundary.
#
# The shortfall it creates is therefore bounded by the TSD length, 9 bp -- the same bound real data
# obeys, and comfortably inside kMaxMCBoundaryDistanceToRepeat. Coverage ramps in over the contig's
# first read length, so the exact MC end lands a base or two inside the TSD; the test asserts that
# the deletion survives, not the precise shortfall.
DEL	.	.	REL606-5	1	31292
DEL	.	.	REL606-5	32079	14220
