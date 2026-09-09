#=GENOME_DIFF	1.0
#=CREATED	07:32:59 09 Sep 2026
#=PROGRAM	breseq 0.50.0 revision 1a30aed3f53b
#=COMMAND	./src/breseq/breseq -j 4 --predict-discordant-pairs -o tests/dp_cross_sequence -c tests/dp_cross_sequence/../data/lambda/lambda-contig.gbk tests/dp_cross_sequence/output.simulated_1.fastq tests/dp_cross_sequence/output.simulated_2.fastq
#=REFSEQ	tests/dp_cross_sequence/../data/lambda/lambda-contig.gbk
#=READSEQ	tests/dp_cross_sequence/output.simulated_1.fastq
#=READSEQ	tests/dp_cross_sequence/output.simulated_2.fastq
#=CONVERTED-BASES	2910000
#=CONVERTED-READS	29100
#=INPUT-BASES	2910000
#=INPUT-READS	29100
#=MAPPED-BASES	2903618
#=MAPPED-READS	29091
MC	1	.	NC_001416-0	1	373	0	0	gene_name=[nu1]	gene_product=[nu1]	ignore=CONTIG_END	left_inside_cov=0	left_outside_cov=NA	locus_tag=[lambdap01]	right_inside_cov=35	right_outside_cov=36
MC	2	.	NC_001416-2	1261	1261	0	0	gene_name=orf-401|orf206b	gene_position=coding (1012/1206 nt)|coding (107/621 nt)	gene_product=Tail fiber protein|hypothetical protein	gene_strand=>|<	left_inside_cov=0	left_outside_cov=53	locus_tag=lambdap27|lambdap90	right_inside_cov=0	right_outside_cov=56	snp_type=|
MC	3	.	NC_001416-4	8516	8516	0	0	gene_name=lambdap78	gene_position=coding (259/534 nt)	gene_product=putative envelope protein	gene_strand=<	left_inside_cov=0	left_outside_cov=65	locus_tag=lambdap78	right_inside_cov=0	right_outside_cov=68
MC	4	.	NC_001416-4	9201	9701	0	0	gene_name=lambdap79/–	gene_position=intergenic (+58/–)	gene_product=hypothetical protein/–	gene_strand=>/–	ignore=CONTIG_END	left_inside_cov=35	left_outside_cov=36	locus_tag=lambdap79/–	right_inside_cov=0	right_outside_cov=NA
CN	5	.	NC_001416-4	9201	9700	0	gene_name=lambdap79/–	gene_position=intergenic (+58/–)	gene_product=hypothetical protein/–	gene_strand=>/–	locus_tag=lambdap79/–	relative_coverage=0.365	tile_size=100
UN	6	.	NC_001416-0	1	11
UN	7	.	NC_001416-0	17	17
UN	8	.	NC_001416-0	9699	9700
UN	9	.	NC_001416-1	1	4
UN	10	.	NC_001416-1	9700	9700
UN	11	.	NC_001416-2	1	2
UN	12	.	NC_001416-2	1261	1261
UN	13	.	NC_001416-2	9701	9701
UN	14	.	NC_001416-3	1	1
UN	15	.	NC_001416-3	9700	9700
UN	16	.	NC_001416-4	1	2
UN	17	.	NC_001416-4	8516	8516
UN	18	.	NC_001416-4	9695	9701
DP	19	.	NC_001416-0	9700	-1	NC_001416-1	1	1	candidate_discordant_count=139	concordant_count=0.0	discordant_count=139	distinct_discordant_count=137	expected_concordant_count=122.4	frequency=1.0000	frequency_lower=0.9784	frequency_upper=1.0000	neg_log10_discordance_p_value=0.1	new_junction_coverage=1.12	side_1_annotate_key=gene	side_1_concordant_count=0	side_1_coverage=0.00	side_1_discordant_count=139	side_1_gene_name=V/–	side_1_gene_position=intergenic (+5/–)	side_1_gene_product=tail component/–	side_1_gene_strand=>/–	side_1_locus_tag=lambdap13/–	side_1_unpaired_count=0	side_2_annotate_key=gene	side_2_concordant_count=0	side_2_coverage=0.00	side_2_discordant_count=139	side_2_gene_name=–/G	side_2_gene_position=intergenic (–/-10)	side_2_gene_product=–/tail component	side_2_gene_strand=–/>	side_2_locus_tag=–/lambdap14	side_2_unpaired_count=1
PD	20	.	NC_001416-1	8982	-1	NC_001416-1	8983	1	ambiguous_pair_count=84	candidate_covering_count=25	distinct_pair_count=11	frequency=0.2683	frequency_lower=0.1583	frequency_upper=0.4048	normal_pair_count=30	position_range=28	reject=PAIR_DISTANCE_SCORE	score=0.1	seed_z_score=-4.51	shifted_pair_count=11	side_1_annotate_key=gene	side_1_gene_name=J	side_1_gene_position=coding (3178/3399 nt)	side_1_gene_product=tail:host specificity protein	side_1_gene_strand=>	side_1_locus_tag=lambdap21	side_2_annotate_key=gene	side_2_gene_name=J	side_2_gene_position=coding (3179/3399 nt)	side_2_gene_product=tail:host specificity protein	side_2_gene_strand=>	side_2_locus_tag=lambdap21	size_shift=-58	size_shift_lower=-90	size_shift_upper=-30	total_pair_count=125
DP	21	.	NC_001416-1	9700	-1	NC_001416-2	1	1	candidate_discordant_count=153	concordant_count=0.0	discordant_count=153	distinct_discordant_count=153	expected_concordant_count=122.4	frequency=1.0000	frequency_lower=0.9806	frequency_upper=1.0000	neg_log10_discordance_p_value=0.0	new_junction_coverage=1.25	side_1_annotate_key=gene	side_1_concordant_count=0	side_1_coverage=0.00	side_1_discordant_count=153	side_1_gene_name=lom	side_1_gene_position=coding (436/436 nt)	side_1_gene_product=outer host membrane	side_1_gene_strand=>	side_1_locus_tag=lambdap26	side_1_unpaired_count=2	side_2_annotate_key=gene	side_2_concordant_count=0	side_2_coverage=0.00	side_2_discordant_count=153	side_2_gene_name=–/orf-401	side_2_gene_position=intergenic (–/-249)	side_2_gene_product=–/Tail fiber protein	side_2_gene_strand=–/>	side_2_locus_tag=–/lambdap27	side_2_unpaired_count=1
DP	22	.	NC_001416-2	9701	-1	NC_001416-3	1	1	candidate_discordant_count=146	concordant_count=0.0	discordant_count=146	distinct_discordant_count=146	expected_concordant_count=122.4	frequency=1.0000	frequency_lower=0.9797	frequency_upper=1.0000	neg_log10_discordance_p_value=0.0	new_junction_coverage=1.19	side_1_annotate_key=gene	side_1_concordant_count=0	side_1_coverage=0.00	side_1_discordant_count=146	side_1_gene_name=xis/–	side_1_gene_position=intergenic (-23/–)	side_1_gene_product=Excisionase/–	side_1_gene_strand=</–	side_1_locus_tag=lambdap34/–	side_1_unpaired_count=2	side_2_annotate_key=gene	side_2_concordant_count=0	side_2_coverage=0.00	side_2_discordant_count=146	side_2_gene_name=–/lambdap35	side_2_gene_position=intergenic (–/+16)	side_2_gene_product=–/hypothetical protein	side_2_gene_strand=–/<	side_2_locus_tag=–/lambdap35	side_2_unpaired_count=0
DP	23	.	NC_001416-3	9700	-1	NC_001416-4	1	1	candidate_discordant_count=148	concordant_count=0.0	discordant_count=148	distinct_discordant_count=148	expected_concordant_count=122.4	frequency=1.0000	frequency_lower=0.9800	frequency_upper=1.0000	neg_log10_discordance_p_value=0.0	new_junction_coverage=1.21	side_1_annotate_key=gene	side_1_concordant_count=0	side_1_coverage=0.00	side_1_discordant_count=148	side_1_gene_name=O	side_1_gene_position=coding (116/116 nt)	side_1_gene_product=DNA replication protein	side_1_gene_strand=>	side_1_locus_tag=lambdap89	side_1_unpaired_count=1	side_2_annotate_key=gene	side_2_concordant_count=0	side_2_coverage=0.00	side_2_discordant_count=148	side_2_gene_name=–/P	side_2_gene_position=intergenic (–/-780)	side_2_gene_product=–/DNA replication protein	side_2_gene_strand=–/>	side_2_locus_tag=–/lambdap61	side_2_unpaired_count=2
