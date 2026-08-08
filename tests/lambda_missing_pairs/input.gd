#=GENOME_DIFF	1.0
#=TITLE	Novel sequence inserted into lambda, for the missing-pair (MP) positive control
#=COMMENT	INT (not INS) takes the inserted sequence from a SECOND loaded reference, so the 3274 bp
#=COMMENT	of novel sequence do not have to be committed literally here. Size 0 makes it a pure
#=COMMENT	insertion -- no reference bases are removed.
#=COMMENT	The donor is a TMV plasmid fragment purely because it is already in tests/data and shares
#=COMMENT	no 24-mer with lambda on either strand; the test is about read geometry, not biology. Its
#=COMMENT	length matters: at ~5x the largest simulated fragment, no read pair can bridge it, so every
#=COMMENT	mate landing inside is unmappable and both flanks get clean MP support.
#=COMMENT	Position 24000 is mid-genome, far from the contig-end guard bands that MP ignores.
INT	1	.	NC_001416	24000	0	TMV-plasmid-truncate-start:1-3274
