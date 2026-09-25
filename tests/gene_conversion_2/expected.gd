#=GENOME_DIFF	1.0
#=CREATED	09:43:06 25 Sep 2026
#=PROGRAM	breseq 0.50.0 revision f30111c975ee
#=COMMAND	./src/breseq/breseq -j 4 -o ./tests/gene_conversion_2 -r ./tests/gene_conversion_2/output.reference.fna ./tests/gene_conversion_2/output.simulated_1.fastq ./tests/gene_conversion_2/output.simulated_2.fastq
#=REFSEQ	./tests/gene_conversion_2/output.reference.fna
#=READSEQ	./tests/gene_conversion_2/output.simulated_1.fastq
#=READSEQ	./tests/gene_conversion_2/output.simulated_2.fastq
#=CONVERTED-BASES	1940000
#=CONVERTED-READS	19400
#=INPUT-BASES	1940000
#=INPUT-READS	19400
#=MAPPED-BASES	1940000
#=MAPPED-READS	19400
CON	1	2,3,6,10,11,22,23	NC_001416	4391	425	NC_001416:12391-12815	gene_name=–/–	gene_position=intergenic (–/–)	gene_product=–/–	gene_strand=–/–	locus_tag=–/–	mutation_category=gene_conversion	position_end=4815	position_start=4391	ref_seq=425-bp
RA	2	.	NC_001416	4401	0	A	C	allele_frequencies=C:1.000e+00	deleted=1	fisher_strand_p_value=1.00000e+00	frequency=1.000e+00	frequency_lower=8.657e-01	frequency_upper=1.000e+00	ks_quality_p_value=1.00000e+00	major_base=C	major_cov=3/7	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=3/7	prediction=consensus	ref_cov=0/0	score=11.6	total_cov=3/7
RA	3	.	NC_001416	4785	0	T	A	allele_frequencies=A:1.000e+00	deleted=1	fisher_strand_p_value=1.00000e+00	frequency=1.000e+00	frequency_lower=8.940e-01	frequency_upper=1.000e+00	ks_quality_p_value=1.00000e+00	major_base=A	major_cov=10/4	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=10/4	prediction=consensus	ref_cov=0/0	score=22.4	total_cov=10/4
RA	4	.	NC_001416	12390	0	G	C	allele_frequencies=C:2.444e-01,G:7.556e-01	consensus_reject=FREQUENCY_CUTOFF	fisher_strand_p_value=2.23778e-01	frequency=2.444e-01	frequency_lower=1.525e-01	frequency_upper=3.537e-01	gene_name=–/–	gene_position=intergenic (–/–)	gene_product=–/–	gene_strand=–/–	ks_quality_p_value=9.97490e-02	locus_tag=–/–	major_base=G	major_cov=17/21	major_frequency=7.556e-01	minor_base=C	minor_cov=10/5	new_cov=10/5	prediction=polymorphism	ref_cov=17/21	score=12.0	snp_type=intergenic	total_cov=27/26
MC	5	.	NC_001416	1	132	0	0	gene_name=–/–	gene_position=intergenic (–/–)	gene_product=–/–	gene_strand=–/–	left_inside_cov=0	left_outside_cov=NA	locus_tag=–/–	right_inside_cov=20	right_outside_cov=21
MC	6	.	NC_001416	4376	4791	0	0	gene_name=–/–	gene_position=intergenic (–/–)	gene_product=–/–	gene_strand=–/–	left_inside_cov=20	left_outside_cov=22	locus_tag=–/–	right_inside_cov=19	right_outside_cov=21
MC	7	.	NC_001416	20661	20661	0	0	gene_name=–/–	gene_position=intergenic (–/–)	gene_product=–/–	gene_strand=–/–	left_inside_cov=0	left_outside_cov=39	locus_tag=–/–	right_inside_cov=0	right_outside_cov=41
MC	8	.	NC_001416	47317	47317	0	0	gene_name=–/–	gene_position=intergenic (–/–)	gene_product=–/–	gene_strand=–/–	left_inside_cov=0	left_outside_cov=36	locus_tag=–/–	right_inside_cov=0	right_outside_cov=36
MC	9	.	NC_001416	48367	48502	0	0	gene_name=–/–	gene_position=intergenic (–/–)	gene_product=–/–	gene_strand=–/–	left_inside_cov=20	left_outside_cov=21	locus_tag=–/–	right_inside_cov=0	right_outside_cov=NA
CN	10	.	NC_001416	4401	4800	0	gene_name=–/–	gene_position=intergenic (–/–)	gene_product=–/–	gene_strand=–/–	locus_tag=–/–	relative_coverage=0	tile_size=100
CN	11	.	NC_001416	12401	12800	2	gene_name=–/–	gene_position=intergenic (–/–)	gene_product=–/–	gene_strand=–/–	locus_tag=–/–	relative_coverage=1.72	tile_size=100
CN	12	.	NC_001416	48401	48500	0	gene_name=–/–	gene_position=intergenic (–/–)	gene_product=–/–	gene_strand=–/–	locus_tag=–/–	relative_coverage=0.152	tile_size=100
UN	13	.	NC_001416	1	37
UN	14	.	NC_001416	4405	4775
UN	15	.	NC_001416	20661	20661
UN	16	.	NC_001416	47317	47317
UN	17	.	NC_001416	48457	48457
UN	18	.	NC_001416	48465	48465
UN	19	.	NC_001416	48468	48468
UN	20	.	NC_001416	48470	48470
UN	21	.	NC_001416	48473	48502
DP	22	.	NC_001416	4431	-1	NC_001416	12357	1	candidate_discordant_count=31	concordant_count=9.0	discordant_count=31	distinct_discordant_count=31	expected_concordant_count=20.0	frequency=0.7750	frequency_lower=0.6402	frequency_upper=0.8773	neg_log10_discordance_p_value=0.0	new_junction_coverage=1.55	side_1_annotate_key=gene	side_1_concordant_count=0	side_1_coverage=0.00	side_1_discordant_count=31	side_1_gene_name=–/–	side_1_gene_position=intergenic (–/–)	side_1_gene_product=–/–	side_1_gene_strand=–/–	side_1_locus_tag=–/–	side_1_unpaired_count=0	side_2_annotate_key=gene	side_2_concordant_count=18	side_2_coverage=0.90	side_2_discordant_count=31	side_2_gene_name=–/–	side_2_gene_position=intergenic (–/–)	side_2_gene_product=–/–	side_2_gene_strand=–/–	side_2_locus_tag=–/–	side_2_unpaired_count=0
DP	23	.	NC_001416	4775	1	NC_001416	12827	-1	candidate_discordant_count=35	concordant_count=8.0	discordant_count=35	distinct_discordant_count=35	expected_concordant_count=20.0	frequency=0.8140	frequency_lower=0.6893	frequency_upper=0.9039	neg_log10_discordance_p_value=0.0	new_junction_coverage=1.75	side_1_annotate_key=gene	side_1_concordant_count=0	side_1_coverage=0.00	side_1_discordant_count=35	side_1_gene_name=–/–	side_1_gene_position=intergenic (–/–)	side_1_gene_product=–/–	side_1_gene_strand=–/–	side_1_locus_tag=–/–	side_1_unpaired_count=0	side_2_annotate_key=gene	side_2_concordant_count=16	side_2_coverage=0.80	side_2_discordant_count=35	side_2_gene_name=–/–	side_2_gene_position=intergenic (–/–)	side_2_gene_product=–/–	side_2_gene_strand=–/–	side_2_locus_tag=–/–	side_2_unpaired_count=0
