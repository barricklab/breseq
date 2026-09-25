#=GENOME_DIFF	1.0
#=CREATED	09:43:26 25 Sep 2026
#=PROGRAM	breseq 0.50.0 revision f30111c975ee
#=COMMAND	./src/breseq/breseq -j 4 -o ./tests/gene_conversion_4 -r ./tests/gene_conversion_4/output.reference.fna ./tests/gene_conversion_4/output.simulated.fastq
#=REFSEQ	./tests/gene_conversion_4/output.reference.fna
#=READSEQ	./tests/gene_conversion_4/output.simulated.fastq
#=CONVERTED-BASES	1892052
#=CONVERTED-READS	52557
#=INPUT-BASES	1892052
#=INPUT-READS	52557
#=MAPPED-BASES	1890648
#=MAPPED-READS	52518
DEL	1	2,5	NC_001416	4002	1200	gene_name=–/–	gene_position=intergenic (–/–)	gene_product=–/–	gene_strand=–/–	locus_tag=–/–	mutation_category=large_deletion	position_end=5201	position_start=4002	ref_seq=1200-bp
MC	2	.	NC_001416	4002	5201	0	0	gene_name=–/–	gene_position=intergenic (–/–)	gene_product=–/–	gene_strand=–/–	left_inside_cov=0	left_outside_cov=40	locus_tag=–/–	right_inside_cov=2	right_outside_cov=42
MC	3	.	NC_001416	20661	20661	0	0	gene_name=–/–	gene_position=intergenic (–/–)	gene_product=–/–	gene_strand=–/–	left_inside_cov=0	left_outside_cov=55	locus_tag=–/–	right_inside_cov=0	right_outside_cov=54
MC	4	.	NC_001416	47317	47317	0	0	gene_name=–/–	gene_position=intergenic (–/–)	gene_product=–/–	gene_strand=–/–	left_inside_cov=0	left_outside_cov=40	locus_tag=–/–	right_inside_cov=0	right_outside_cov=45
JC	5	.	NC_001416	4001	-1	NC_001416	5202	1	0	alignment_overlap=4	coverage_minus=19	coverage_plus=20	flanking_left=36	flanking_right=36	frequency=1.000e+00	frequency_lower=9.180e-01	frequency_upper=1.000e+00	junction_effective_depth=35.00	junction_mixture_iterations=1	junction_possible_overlap_registers=29	junction_possible_overlap_registers_before_trimming=31	key=NC_001416__4001__-1__NC_001416__5198__1__4____36__36__0__0	max_left=31	max_left_minus=29	max_left_plus=31	max_min_left=15	max_min_left_minus=15	max_min_left_plus=15	max_min_right=16	max_min_right_minus=16	max_min_right_plus=16	max_pos_hash_score=62	max_right=31	max_right_minus=31	max_right_plus=30	neg_log10_pos_hash_p_value=0.1	new_junction_coverage=1.09	new_junction_read_count=35	new_junction_reference_weighted_read_count=0.00	new_junction_weighted_read_count=35.00	pos_hash_score=29	prediction=consensus	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=0.00	side_1_gene_name=–/–	side_1_gene_position=intergenic (–/–)	side_1_gene_product=–/–	side_1_gene_strand=–/–	side_1_locus_tag=–/–	side_1_overlap=4	side_1_possible_overlap_registers=33	side_1_possible_overlap_registers_before_trimming=35	side_1_read_count=0	side_1_redundant=0	side_1_weighted_read_count=0.00	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=0.00	side_2_gene_name=–/–	side_2_gene_position=intergenic (–/–)	side_2_gene_product=–/–	side_2_gene_strand=–/–	side_2_locus_tag=–/–	side_2_overlap=0	side_2_possible_overlap_registers=28	side_2_possible_overlap_registers_before_trimming=31	side_2_read_count=0	side_2_redundant=0	side_2_weighted_read_count=0.00	total_non_overlap_reads=39
UN	6	.	NC_001416	1	8
UN	7	.	NC_001416	4002	5201
UN	8	.	NC_001416	20661	20661
UN	9	.	NC_001416	47317	47317
UN	10	.	NC_001416	48501	48502
