#=GENOME_DIFF	1.0
#=TITLE	Some of the mutations that tests/lambda_mixed_pop predicts, to apply before re-analyzing its reads
DEL	1	.	NC_001416.1	139	1
INS	2	.	NC_001416.1	14266	G
SNP	3	.	NC_001416.1	20661	G
DEL	4	.	NC_001416.1	21738	5996
SUB	5	.	NC_001416.1	47977	2	AC
# The reads do not carry this insertion, so breseq deletes it again: that DEL and its RA
# evidence sit inside inserted sequence and report original_position=30000 with an offset.
INS	6	.	NC_001416.1	30000	GGATCCTAGG
