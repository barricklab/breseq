#=GENOME_DIFF	1.0
#=CREATED	21:06:11 13 Aug 2026
#=PROGRAM	breseq 0.50.0 revision d3ec44338e89
#=COMMAND	./src/breseq/breseq -j 4 --predict-pair-distance --predict-discordant-pairs -o ./tests/lambda_pair_distances -r ./tests/lambda_pair_distances/../data/lambda/lambda.gbk ./tests/lambda_pair_distances/output.simulated_1.fastq ./tests/lambda_pair_distances/output.simulated_2.fastq
#=REFSEQ	./tests/lambda_pair_distances/../data/lambda/lambda.gbk
#=READSEQ	./tests/lambda_pair_distances/output.simulated_1.fastq
#=READSEQ	./tests/lambda_pair_distances/output.simulated_2.fastq
#=CONVERTED-BASES	2910000
#=CONVERTED-READS	29100
#=INPUT-BASES	2910000
#=INPUT-READS	29100
#=MAPPED-BASES	2885290
#=MAPPED-READS	28853
DEL	1	3,7,15	NC_001416	12000	200	gene_name=H	gene_position=coding (1459-1658/2562 nt)	gene_product=tail component	gene_strand=>	genes_inactivated=H	locus_tag=lambdap16	locus_tags_inactivated=lambdap16	mutation_category=large_deletion	position_end=12199	position_start=12000	ref_seq=200-bp
MC	2	.	NC_001416	1	373	0	0	gene_name=[nu1]	gene_product=[nu1]	left_inside_cov=0	left_outside_cov=NA	locus_tag=[lambdap01]	right_inside_cov=35	right_outside_cov=36
MC	3	.	NC_001416	12000	12199	0	0	gene_name=H	gene_position=coding (1459-1658/2562 nt)	gene_product=tail component	gene_strand=>	left_inside_cov=0	left_outside_cov=65	locus_tag=lambdap16	right_inside_cov=0	right_outside_cov=66
MC	4	.	NC_001416	20661	20661	0	0	gene_name=orf-401|orf206b	gene_position=coding (1012/1206 nt)|coding (107/621 nt)	gene_product=Tail fiber protein|hypothetical protein	gene_strand=>|<	left_inside_cov=0	left_outside_cov=61	locus_tag=lambdap27|lambdap90	right_inside_cov=0	right_outside_cov=62	snp_type=|
MC	5	.	NC_001416	47317	47317	0	0	gene_name=lambdap78	gene_position=coding (259/534 nt)	gene_product=putative envelope protein	gene_strand=<	left_inside_cov=0	left_outside_cov=67	locus_tag=lambdap78	right_inside_cov=0	right_outside_cov=70
MC	6	.	NC_001416	48002	48502	0	0	gene_name=lambdap79/–	gene_position=intergenic (+58/–)	gene_product=hypothetical protein/–	gene_strand=>/–	left_inside_cov=35	left_outside_cov=36	locus_tag=lambdap79/–	right_inside_cov=0	right_outside_cov=NA
JC	7	.	NC_001416	11999	-1	NC_001416	12200	1	0	alignment_overlap=0	coverage_minus=33	coverage_plus=31	flanking_left=100	flanking_right=100	frequency=1.000e+00	frequency_lower=9.521e-01	frequency_upper=1.000e+00	junction_effective_depth=61.00	junction_mixture_iterations=1	junction_possible_overlap_registers=94	junction_possible_overlap_registers_before_trimming=99	key=NC_001416__11999__-1__NC_001416__12200__1__0____100__100__0__0	max_left=99	max_left_minus=99	max_left_plus=97	max_min_left=49	max_min_left_minus=49	max_min_left_plus=48	max_min_right=49	max_min_right_minus=48	max_min_right_plus=49	max_pos_hash_score=198	max_right=95	max_right_minus=95	max_right_plus=90	neg_log10_pos_hash_p_value=0.1	new_junction_coverage=1.08	new_junction_read_count=61	new_junction_reference_weighted_read_count=0.00	new_junction_weighted_read_count=61.00	pos_hash_score=57	prediction=consensus	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=0.00	side_1_gene_name=H	side_1_gene_position=coding (1458/2562 nt)	side_1_gene_product=tail component	side_1_gene_strand=>	side_1_locus_tag=lambdap16	side_1_overlap=0	side_1_possible_overlap_registers=94	side_1_possible_overlap_registers_before_trimming=99	side_1_read_count=0	side_1_redundant=0	side_1_weighted_read_count=0.00	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=0.00	side_2_gene_name=H	side_2_gene_position=coding (1659/2562 nt)	side_2_gene_product=tail component	side_2_gene_strand=>	side_2_locus_tag=lambdap16	side_2_overlap=0	side_2_possible_overlap_registers=94	side_2_possible_overlap_registers_before_trimming=99	side_2_read_count=0	side_2_redundant=0	side_2_weighted_read_count=0.00	total_non_overlap_reads=64
UN	8	.	NC_001416	1	11
UN	9	.	NC_001416	17	17
UN	10	.	NC_001416	12000	12199
UN	11	.	NC_001416	20661	20661
UN	12	.	NC_001416	24000	24002
UN	13	.	NC_001416	47317	47317
UN	14	.	NC_001416	48496	48502
PD	15	.	NC_001416	11999	-1	NC_001416	12200	1	ambiguous_pair_count=8	candidate_covering_count=129	distinct_pair_count=120	frequency=0.9917	frequency_lower=0.9614	frequency_upper=0.9996	normal_pair_count=1	position_range=1	score=53.0	seed_z_score=19.09	shifted_pair_count=120	side_1_annotate_key=gene	side_1_gene_name=H	side_1_gene_position=coding (1458/2562 nt)	side_1_gene_product=tail component	side_1_gene_strand=>	side_1_locus_tag=lambdap16	side_2_annotate_key=gene	side_2_gene_name=H	side_2_gene_position=coding (1659/2562 nt)	side_2_gene_product=tail component	side_2_gene_strand=>	side_2_locus_tag=lambdap16	size_shift=200	size_shift_lower=195	size_shift_upper=208	snapped_to_junction=1	total_pair_count=129
PD	16	.	NC_001416	23996	-1	NC_001416	23997	1	ambiguous_pair_count=3	candidate_covering_count=59	distinct_pair_count=50	frequency=0.9444	frequency_lower=0.8626	frequency_upper=0.9847	normal_pair_count=3	position_range=7	score=22.6	seed_z_score=-12.86	shifted_pair_count=51	side_1_annotate_key=gene	side_1_gene_name=ea47/ea31	side_1_gene_position=intergenic (-78/+513)	side_1_gene_product=ea47/ea31	side_1_gene_strand=</<	side_1_locus_tag=lambdap80/lambdap81	side_2_annotate_key=gene	side_2_gene_name=ea47/ea31	side_2_gene_position=intergenic (-79/+512)	side_2_gene_product=ea47/ea31	side_2_gene_strand=</<	side_2_locus_tag=lambdap80/lambdap81	size_shift=-208	size_shift_lower=-233	size_shift_upper=-197	total_pair_count=57
