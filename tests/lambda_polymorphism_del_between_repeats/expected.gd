#=GENOME_DIFF	1.0
#=CREATED	16:34:54 16 Nov 2026
#=PROGRAM	breseq 0.50.0 revision 2db563856b58
#=COMMAND	./src/breseq/breseq -j 4 -p -o ./tests/lambda_polymorphism_del_between_repeats -r ./tests/lambda_polymorphism_del_between_repeats/output_sim/lambda_two_is150.gbk ./tests/lambda_polymorphism_del_between_repeats/output_sim/ancestor_reads.fastq ./tests/lambda_polymorphism_del_between_repeats/output_sim/deleted_reads.fastq
#=REFSEQ	./tests/lambda_polymorphism_del_between_repeats/output_sim/lambda_two_is150.gbk
#=READSEQ	./tests/lambda_polymorphism_del_between_repeats/output_sim/ancestor_reads.fastq
#=READSEQ	./tests/lambda_polymorphism_del_between_repeats/output_sim/deleted_reads.fastq
#=CONVERTED-BASES	3893600
#=CONVERTED-READS	38936
#=INPUT-BASES	3893600
#=INPUT-READS	38936
#=MAPPED-BASES	3893600
#=MAPPED-READS	38936
DEL	1	4	NC_001416	16446	5446	between=IS150	cn_frequency=4.78000000e-01	frequency=4.78000000e-01	gene_name=[I]–[lom]	gene_product=[I],J,insJ-4,insK-4,[lom]	genes_inactivated=I,J,insJ-4,insK-4,lom	locus_tag=[lambdap20]–[lambdap26]	locus_tags_inactivated=lambdap20,lambdap21,ECB_03408,ECB_03409,lambdap26	mutation_category=large_deletion	position_end=21891	position_start=16446	ref_seq=5446-bp
MC	2	.	NC_001416	23553	23553	0	0	gene_name=orf-401|orf206b	gene_position=coding (1012/1206 nt)|coding (107/621 nt)	gene_product=Tail fiber protein|hypothetical protein	gene_strand=>|<	left_inside_cov=0	left_outside_cov=75	locus_tag=lambdap27|lambdap90	right_inside_cov=0	right_outside_cov=74	snp_type=|
MC	3	.	NC_001416	50209	50209	0	0	gene_name=lambdap78	gene_position=coding (259/534 nt)	gene_product=putative envelope protein	gene_strand=<	left_inside_cov=0	left_outside_cov=78	locus_tag=lambdap78	right_inside_cov=0	right_outside_cov=79
CN	4	.	NC_001416	16501	20400	0.5	gene_name=[I]–J	gene_product=[I],J	locus_tag=[lambdap20]–[lambdap21]	relative_coverage=0.522	relative_coverage_se=0.0126	tile_size=100
UN	5	.	NC_001416	1	8
UN	6	.	NC_001416	15095	16362
UN	7	.	NC_001416	16366	16374
UN	8	.	NC_001416	20541	21797
UN	9	.	NC_001416	23553	23553
UN	10	.	NC_001416	50209	50209
UN	11	.	NC_001416	51386	51394
