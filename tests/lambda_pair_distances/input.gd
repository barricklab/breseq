#=GENOME_DIFF	1.0
#=TITLE	A 200 bp deletion, a 200 bp insertion and a 3 kb insertion in lambda, for the pair-distance (PD) and missing-pair (MP) controls
#=COMMENT	These are the two events PD exists to find, one in each direction. Read pairs whose
#=COMMENT	unsequenced middle gap spans the DEL map 200 bp FARTHER apart on the reference than the
#=COMMENT	library says they should, because the reference still has the 200 bases the sample lost.
#=COMMENT	Pairs spanning the insertion map 200 bp CLOSER, for the mirror-image reason.
#=COMMENT	Neither shift is large enough to make any single pair discordant, which is precisely the
#=COMMENT	gap PD fills: DP tests each pair on its own and has no short-side cutoff at all.
#=COMMENT	A tandem duplication would NOT belong here. Pairs crossing its junction map in reversed
#=COMMENT	(everted) order rather than at a shifted distance, which is DP's signature, not PD's.
#=COMMENT	INT (not INS) takes the inserted sequence from a SECOND loaded reference, so the 200 bp of
#=COMMENT	novel sequence do not have to be committed literally here. Size 0 makes it a pure
#=COMMENT	insertion -- no reference bases are removed. The donor is a TMV plasmid fragment purely
#=COMMENT	because it is already in tests/data and shares no 24-mer with lambda on either strand.
#=COMMENT	The two events are 12 kb apart, and both are far from the contig-end guard bands PD
#=COMMENT	ignores, so neither can interfere with the other's covering pairs.
#=COMMENT	The third event is the missing-pair (MP) control: 3074 bp of the same TMV fragment (the
#=COMMENT	part the 200 bp insertion does not use), another 12 kb along. It is far too long for any
#=COMMENT	pair to span, so PD has nothing to say about it; mates that land inside it have nowhere
#=COMMENT	to align and become singletons on both flanks of the insertion point.
DEL	1	.	NC_001416	12000	200
INT	2	.	NC_001416	24000	0	TMV-plasmid-truncate-start:1-200
INT	3	.	NC_001416	36000	0	TMV-plasmid-truncate-start:201-3274
