#=GENOME_DIFF	1.0
#=CREATED	09:41:06 28 Nov 2026
#=PROGRAM	breseq 0.50.0 revision 43c19b63ebc3
#=COMMAND	./src/breseq/breseq -j 4 -o ./tests/gene_conversion_1 -r ./tests/gene_conversion_1/output.reference.fna ./tests/gene_conversion_1/output.simulated.fastq
#=REFSEQ	./tests/gene_conversion_1/output.reference.fna
#=READSEQ	./tests/gene_conversion_1/output.simulated.fastq
#=CONVERTED-BASES	1572048
#=CONVERTED-READS	43668
#=INPUT-BASES	1572048
#=INPUT-READS	43668
#=MAPPED-BASES	1570641
#=MAPPED-READS	43629
CON	1	7	NC_001416	4391	425	NC_001416:12391-12815	gene_name=–/–	gene_position=intergenic (–/–)	gene_product=–/–	gene_strand=–/–	locus_tag=–/–	mutation_category=gene_conversion	position_end=4815	position_start=4391	ref_seq=425-bp
CON	2	8	NC_001416	16348	459	NC_001416:24853-24395	gene_name=–/–	gene_position=intergenic (–/–)	gene_product=–/–	gene_strand=–/–	locus_tag=–/–	mutation_category=gene_conversion	position_end=16806	position_start=16348	ref_seq=459-bp
DEL	3	10,13,14	NC_001416	28001	1200	gene_name=–/–	gene_position=intergenic (–/–)	gene_product=–/–	gene_strand=–/–	locus_tag=–/–	mutation_category=large_deletion	position_end=29200	position_start=28001	ref_seq=1200-bp
DEL	4	11,15	NC_001416	37582	8000	gene_name=–/–	gene_position=intergenic (–/–)	gene_product=–/–	gene_strand=–/–	locus_tag=–/–	mutation_category=large_deletion	position_end=45581	position_start=37582	ref_seq=8000-bp
RA	5	.	NC_001416	37582	0	C	T	allele_frequencies=T:1.000e+00	deleted=1	fisher_strand_p_value=1.00000e+00	frequency=1.000e+00	frequency_lower=9.209e-01	frequency_upper=1.000e+00	ks_quality_p_value=1.00000e+00	major_base=T	major_cov=9/9	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=9/9	prediction=consensus	ref_cov=0/0	score=25.9	total_cov=9/9
MC	6	.	NC_001416	1	20	0	0	gene_name=–/–	gene_position=intergenic (–/–)	gene_product=–/–	gene_strand=–/–	left_inside_cov=0	left_outside_cov=NA	locus_tag=–/–	right_inside_cov=20	right_outside_cov=22
MC	7	.	NC_001416	4385	4801	40	20	gene_name=–/–	gene_position=intergenic (–/–)	gene_product=–/–	gene_strand=–/–	left_inside_cov=20	left_outside_cov=23	locus_tag=–/–	right_inside_cov=19	right_outside_cov=23
MC	8	.	NC_001416	16339	16812	0	191	gene_name=–/–	gene_position=intergenic (–/–)	gene_product=–/–	gene_strand=–/–	left_inside_cov=20	left_outside_cov=22	locus_tag=–/–	right_inside_cov=20	right_outside_cov=21
MC	9	.	NC_001416	20661	20661	0	0	gene_name=–/–	gene_position=intergenic (–/–)	gene_product=–/–	gene_strand=–/–	left_inside_cov=0	left_outside_cov=42	locus_tag=–/–	right_inside_cov=0	right_outside_cov=44
MC	10	.	NC_001416	28001	29200	0	0	gene_name=–/–	gene_position=intergenic (–/–)	gene_product=–/–	gene_strand=–/–	left_inside_cov=0	left_outside_cov=40	locus_tag=–/–	right_inside_cov=0	right_outside_cov=39
MC	11	.	NC_001416	37582	45609	50	42	gene_name=–/–	gene_position=intergenic (–/–)	gene_product=–/–	gene_strand=–/–	left_inside_cov=18	left_outside_cov=23	locus_tag=–/–	right_inside_cov=20	right_outside_cov=21
MC	12	.	NC_001416	47317	47317	0	0	gene_name=–/–	gene_position=intergenic (–/–)	gene_product=–/–	gene_strand=–/–	left_inside_cov=0	left_outside_cov=46	locus_tag=–/–	right_inside_cov=0	right_outside_cov=47
JC	13	.	NC_001416	28000	-1	NC_001416	29201	1	0	alignment_overlap=0	coverage_minus=21	coverage_plus=18	flanking_left=36	flanking_right=36	frequency=1.000e+00	frequency_lower=9.222e-01	frequency_upper=1.000e+00	junction_effective_depth=37.00	junction_mixture_iterations=1	junction_possible_overlap_registers=32	junction_possible_overlap_registers_before_trimming=35	key=NC_001416__28000__-1__NC_001416__29201__1__0____36__36__0__0	max_left=35	max_left_minus=35	max_left_plus=31	max_min_left=17	max_min_left_minus=15	max_min_left_plus=17	max_min_right=18	max_min_right_minus=18	max_min_right_plus=16	max_pos_hash_score=70	max_right=35	max_right_minus=35	max_right_plus=34	neg_log10_pos_hash_p_value=0.0	new_junction_coverage=1.03	new_junction_read_count=37	new_junction_reference_weighted_read_count=0.00	new_junction_weighted_read_count=37.00	pos_hash_score=33	prediction=consensus	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=0.00	side_1_gene_name=–/–	side_1_gene_position=intergenic (–/–)	side_1_gene_product=–/–	side_1_gene_strand=–/–	side_1_locus_tag=–/–	side_1_overlap=0	side_1_possible_overlap_registers=30	side_1_possible_overlap_registers_before_trimming=35	side_1_read_count=0	side_1_redundant=0	side_1_weighted_read_count=0.00	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=0.00	side_2_gene_name=–/–	side_2_gene_position=intergenic (–/–)	side_2_gene_product=–/–	side_2_gene_strand=–/–	side_2_locus_tag=–/–	side_2_overlap=0	side_2_possible_overlap_registers=32	side_2_possible_overlap_registers_before_trimming=35	side_2_read_count=0	side_2_redundant=0	side_2_weighted_read_count=0.00	total_non_overlap_reads=39
CN	14	.	NC_001416	29001	29200	0	gene_name=–/–	gene_position=intergenic (–/–)	gene_product=–/–	gene_strand=–/–	locus_tag=–/–	redundant_before=1000	relative_coverage=0	tile_size=100
CN	15	.	NC_001416	38201	45200	0	gene_name=–/–	gene_position=intergenic (–/–)	gene_product=–/–	gene_strand=–/–	locus_tag=–/–	redundant_after=1000	redundant_before=1000	relative_coverage=0	tile_size=100
UN	16	.	NC_001416	1	8
UN	17	.	NC_001416	4399	4791
UN	18	.	NC_001416	16033	16034
UN	19	.	NC_001416	16140	16148
UN	20	.	NC_001416	16252	16252
UN	21	.	NC_001416	16370	16797
UN	22	.	NC_001416	16868	16882
UN	23	.	NC_001416	17147	17151
UN	24	.	NC_001416	20661	20661
UN	25	.	NC_001416	24046	24056
UN	26	.	NC_001416	24316	24334
UN	27	.	NC_001416	24445	24548
UN	28	.	NC_001416	24814	24821
UN	29	.	NC_001416	24949	24952
UN	30	.	NC_001416	25053	25058
UN	31	.	NC_001416	28001	29200
UN	32	.	NC_001416	32236	32242
UN	33	.	NC_001416	32410	32427
UN	34	.	NC_001416	32680	32722
UN	35	.	NC_001416	37374	37385
UN	36	.	NC_001416	37525	37538
UN	37	.	NC_001416	37589	45583
UN	38	.	NC_001416	46055	46069
UN	39	.	NC_001416	47317	47317
UN	40	.	NC_001416	48500	48502
