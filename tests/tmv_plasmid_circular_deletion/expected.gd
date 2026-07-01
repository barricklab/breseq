#=GENOME_DIFF	1.0
#=CREATED	20:36:54 30 Jun 2026
#=PROGRAM	breseq 0.50.0
#=COMMAND	./src/c/breseq/breseq -j 4 -o ./tests/tmv_plasmid_circular_deletion -r ./tests/tmv_plasmid_circular_deletion/../data/tmv_plasmid/tmv-plasmid.gbk ./tests/tmv_plasmid_circular_deletion/../data/tmv_plasmid/D3-9_1P.fastq.gz ./tests/tmv_plasmid_circular_deletion/../data/tmv_plasmid/D3-9_2P.fastq.gz
#=REFSEQ	./tests/tmv_plasmid_circular_deletion/../data/tmv_plasmid/tmv-plasmid.gbk
#=READSEQ	./tests/tmv_plasmid_circular_deletion/../data/tmv_plasmid/D3-9_1P.fastq.gz
#=READSEQ	./tests/tmv_plasmid_circular_deletion/../data/tmv_plasmid/D3-9_2P.fastq.gz
#=CONVERTED-BASES	2733978
#=CONVERTED-READS	18356
#=INPUT-BASES	2738037
#=INPUT-READS	18386
#=MAPPED-BASES	2286649
#=MAPPED-READS	15438
DEL	1	4,5,7	TMV-plasmid	10205	7344	gene_name=–/–	gene_position=intergenic (–/–)	gene_product=–/–	gene_strand=–/–	locus_tag=–/–	mutation_category=large_deletion	position_end=17548	position_start=10205	ref_seq=7344-bp
RA	2	.	TMV-plasmid	1565	0	C	T	consensus_score=10.3	deleted=1	frequency=1	major_base=T	major_cov=1/3	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=1/3	polymorphism_frequency=1.000e+00	polymorphism_score=NA	prediction=consensus	ref_cov=0/0	total_cov=1/3
RA	3	.	TMV-plasmid	5198	0	A	G	consensus_score=18.2	deleted=1	frequency=1	major_base=G	major_cov=3/3	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=3/3	polymorphism_frequency=1.000e+00	polymorphism_score=NA	prediction=consensus	ref_cov=0/0	total_cov=3/3
MC	4	.	TMV-plasmid	1	7137	0	0	gene_name=–/–	gene_position=intergenic (–/–)	gene_product=–/–	gene_strand=–/–	left_inside_cov=1	left_outside_cov=NA	locus_tag=–/–	right_inside_cov=41	right_outside_cov=812
MC	5	.	TMV-plasmid	10205	10411	0	0	gene_name=–/–	gene_position=intergenic (–/–)	gene_product=–/–	gene_strand=–/–	left_inside_cov=85	left_outside_cov=858	locus_tag=–/–	right_inside_cov=1	right_outside_cov=NA
JC	6	.	TMV-plasmid	7033	1	TMV-plasmid	10253	-1	0	alignment_overlap=0	coverage_minus=22	coverage_plus=17	flanking_left=151	flanking_right=151	frequency=1	junction_possible_overlap_registers=147	key=TMV-plasmid__7033__1__TMV-plasmid__10253__-1__0____151__151__0__0	max_left=149	max_left_minus=149	max_left_plus=147	max_min_left=72	max_min_left_minus=71	max_min_left_plus=72	max_min_right=75	max_min_right_minus=75	max_min_right_plus=60	max_pos_hash_score=296	max_right=150	max_right_minus=150	max_right_plus=138	neg_log10_pos_hash_p_value=0.9	new_junction_coverage=0.19	new_junction_read_count=39	polymorphism_frequency=8.571e-01	pos_hash_score=28	prediction=consensus	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=0.05	side_1_gene_name=–/–	side_1_gene_position=intergenic (–/–)	side_1_gene_product=–/–	side_1_gene_strand=–/–	side_1_locus_tag=–/–	side_1_overlap=0	side_1_possible_overlap_registers=147	side_1_read_count=10	side_1_redundant=0	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=0.01	side_2_gene_name=–/–	side_2_gene_position=intergenic (–/–)	side_2_gene_product=–/–	side_2_gene_strand=–/–	side_2_locus_tag=–/–	side_2_overlap=0	side_2_possible_overlap_registers=147	side_2_read_count=3	side_2_redundant=0	total_non_overlap_reads=39
JC	7	.	TMV-plasmid	7138	1	TMV-plasmid	10204	-1	0	alignment_overlap=9	coverage_minus=388	coverage_plus=373	flanking_left=151	flanking_right=151	frequency=1	junction_possible_overlap_registers=138	key=TMV-plasmid__7138__1__TMV-plasmid__10213__-1__9____151__151__0__0	max_left=141	max_left_minus=141	max_left_plus=141	max_min_left=69	max_min_left_minus=69	max_min_left_plus=69	max_min_right=71	max_min_right_minus=71	max_min_right_plus=71	max_pos_hash_score=278	max_right=141	max_right_minus=141	max_right_plus=141	neg_log10_pos_hash_p_value=0.0	new_junction_coverage=3.97	new_junction_read_count=769	polymorphism_frequency=9.728e-01	pos_hash_score=203	prediction=consensus	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=0.16	side_1_gene_name=–/–	side_1_gene_position=intergenic (–/–)	side_1_gene_product=–/–	side_1_gene_strand=–/–	side_1_locus_tag=–/–	side_1_overlap=9	side_1_possible_overlap_registers=147	side_1_read_count=33	side_1_redundant=0	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=0.06	side_2_gene_name=–/–	side_2_gene_position=intergenic (–/–)	side_2_gene_product=–/–	side_2_gene_strand=–/–	side_2_locus_tag=–/–	side_2_overlap=0	side_2_possible_overlap_registers=138	side_2_read_count=12	side_2_redundant=0	total_non_overlap_reads=761
JC	8	.	TMV-plasmid	7503	1	TMV-plasmid	8682	-1	0	alignment_overlap=1	coverage_minus=2	coverage_plus=1	flanking_left=151	flanking_right=151	frequency=4.130e-03	junction_possible_overlap_registers=146	key=TMV-plasmid__7503__1__TMV-plasmid__8683__-1__1____151__151__0__0	max_left=148	max_left_minus=148	max_left_plus=78	max_min_left=0	max_min_left_minus=0	max_min_left_plus=0	max_min_right=72	max_min_right_minus=72	max_min_right_plus=70	max_pos_hash_score=294	max_right=72	max_right_minus=72	max_right_plus=70	neg_log10_pos_hash_p_value=5.9	new_junction_coverage=0.01	new_junction_read_count=3	no_show=1	polymorphism_frequency=4.130e-03	pos_hash_score=3	prediction=polymorphism	reject=COVERAGE_EVENNESS_SKEW,FREQUENCY_CUTOFF	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=3.72	side_1_overlap=1	side_1_possible_overlap_registers=147	side_1_read_count=768	side_1_redundant=0	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=3.33	side_2_overlap=0	side_2_possible_overlap_registers=146	side_2_read_count=684	side_2_redundant=0	total_non_overlap_reads=3
JC	9	.	TMV-plasmid	7508	-1	TMV-plasmid	7607	1	0	alignment_overlap=0	coverage_minus=3	coverage_plus=1	flanking_left=151	flanking_right=151	frequency=5.686e-03	junction_possible_overlap_registers=147	key=TMV-plasmid__7508__-1__TMV-plasmid__7607__1__0____151__151__0__0	max_left=150	max_left_minus=150	max_left_plus=82	max_min_left=0	max_min_left_minus=0	max_min_left_plus=0	max_min_right=70	max_min_right_minus=55	max_min_right_plus=70	max_pos_hash_score=296	max_right=70	max_right_minus=55	max_right_plus=70	neg_log10_pos_hash_p_value=5.9	new_junction_coverage=0.02	new_junction_read_count=4	no_show=1	polymorphism_frequency=5.686e-03	pos_hash_score=3	prediction=polymorphism	reject=COVERAGE_EVENNESS_SKEW,FREQUENCY_CUTOFF	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=3.62	side_1_overlap=0	side_1_possible_overlap_registers=147	side_1_read_count=748	side_1_redundant=0	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=3.15	side_2_overlap=0	side_2_possible_overlap_registers=147	side_2_read_count=651	side_2_redundant=0	total_non_overlap_reads=4
JC	10	.	TMV-plasmid	7512	1	TMV-plasmid	8986	-1	0	alignment_overlap=0	coverage_minus=2	coverage_plus=1	flanking_left=151	flanking_right=151	frequency=4.149e-03	junction_possible_overlap_registers=147	key=TMV-plasmid__7512__1__TMV-plasmid__8986__-1__0____151__151__0__0	max_left=95	max_left_minus=95	max_left_plus=46	max_min_left=46	max_min_left_minus=0	max_min_left_plus=46	max_min_right=75	max_min_right_minus=75	max_min_right_plus=0	max_pos_hash_score=296	max_right=105	max_right_minus=75	max_right_plus=105	neg_log10_pos_hash_p_value=5.9	new_junction_coverage=0.01	new_junction_read_count=3	no_show=1	polymorphism_frequency=4.149e-03	pos_hash_score=3	prediction=polymorphism	reject=COVERAGE_EVENNESS_SKEW,FREQUENCY_CUTOFF	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=3.61	side_1_overlap=0	side_1_possible_overlap_registers=147	side_1_read_count=745	side_1_redundant=0	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=3.36	side_2_overlap=0	side_2_possible_overlap_registers=147	side_2_read_count=695	side_2_redundant=0	total_non_overlap_reads=3
JC	11	.	TMV-plasmid	7526	1	TMV-plasmid	7778	1	0	alignment_overlap=0	coverage_minus=2	coverage_plus=2	flanking_left=151	flanking_right=151	frequency=5.682e-03	junction_possible_overlap_registers=147	key=TMV-plasmid__7526__1__TMV-plasmid__7778__1__0____151__151__0__0	max_left=150	max_left_minus=150	max_left_plus=69	max_min_left=69	max_min_left_minus=0	max_min_left_plus=69	max_min_right=53	max_min_right_minus=53	max_min_right_plus=0	max_pos_hash_score=296	max_right=146	max_right_minus=53	max_right_plus=146	neg_log10_pos_hash_p_value=5.4	new_junction_coverage=0.02	new_junction_read_count=4	polymorphism_frequency=5.682e-03	pos_hash_score=4	prediction=polymorphism	reject=COVERAGE_EVENNESS_SKEW,FREQUENCY_CUTOFF	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=3.39	side_1_gene_name=–/–	side_1_gene_position=intergenic (–/–)	side_1_gene_product=–/–	side_1_gene_strand=–/–	side_1_locus_tag=–/–	side_1_overlap=0	side_1_possible_overlap_registers=147	side_1_read_count=701	side_1_redundant=0	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=3.38	side_2_gene_name=–/–	side_2_gene_position=intergenic (–/–)	side_2_gene_product=–/–	side_2_gene_strand=–/–	side_2_locus_tag=–/–	side_2_overlap=0	side_2_possible_overlap_registers=147	side_2_read_count=699	side_2_redundant=0	total_non_overlap_reads=4
JC	12	.	TMV-plasmid	7579	1	TMV-plasmid	9496	1	0	alignment_overlap=0	coverage_minus=1	coverage_plus=6	flanking_left=151	flanking_right=151	frequency=8.866e-03	junction_possible_overlap_registers=147	key=TMV-plasmid__7579__1__TMV-plasmid__9496__1__0____151__151__0__0	max_left=68	max_left_minus=68	max_left_plus=48	max_min_left=68	max_min_left_minus=68	max_min_left_plus=48	max_min_right=0	max_min_right_minus=0	max_min_right_plus=0	max_pos_hash_score=296	max_right=149	max_right_minus=82	max_right_plus=149	neg_log10_pos_hash_p_value=4.9	new_junction_coverage=0.03	new_junction_read_count=7	polymorphism_frequency=8.866e-03	pos_hash_score=5	prediction=polymorphism	reject=COVERAGE_EVENNESS_SKEW,FREQUENCY_CUTOFF	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=3.11	side_1_gene_name=–/–	side_1_gene_position=intergenic (–/–)	side_1_gene_product=–/–	side_1_gene_strand=–/–	side_1_locus_tag=–/–	side_1_overlap=0	side_1_possible_overlap_registers=147	side_1_read_count=643	side_1_redundant=0	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=4.46	side_2_gene_name=–/–	side_2_gene_position=intergenic (–/–)	side_2_gene_product=–/–	side_2_gene_strand=–/–	side_2_locus_tag=–/–	side_2_overlap=0	side_2_possible_overlap_registers=147	side_2_read_count=922	side_2_redundant=0	total_non_overlap_reads=7
JC	13	.	TMV-plasmid	7616	1	TMV-plasmid	8679	-1	0	alignment_overlap=2	coverage_minus=2	coverage_plus=2	flanking_left=151	flanking_right=151	frequency=7.381e-03	junction_possible_overlap_registers=145	key=TMV-plasmid__7616__1__TMV-plasmid__8681__-1__2____151__151__0__0	max_left=148	max_left_minus=148	max_left_plus=75	max_min_left=0	max_min_left_minus=0	max_min_left_plus=0	max_min_right=74	max_min_right_minus=54	max_min_right_plus=74	max_pos_hash_score=292	max_right=74	max_right_minus=54	max_right_plus=74	neg_log10_pos_hash_p_value=5.8	new_junction_coverage=0.02	new_junction_read_count=5	no_show=1	polymorphism_frequency=7.381e-03	pos_hash_score=3	prediction=polymorphism	reject=COVERAGE_EVENNESS_SKEW,FREQUENCY_CUTOFF	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=3.29	side_1_overlap=2	side_1_possible_overlap_registers=147	side_1_read_count=679	side_1_redundant=0	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=3.31	side_2_overlap=0	side_2_possible_overlap_registers=145	side_2_read_count=675	side_2_redundant=0	total_non_overlap_reads=4
JC	14	.	TMV-plasmid	7657	1	TMV-plasmid	9272	1	0	alignment_overlap=2	coverage_minus=2	coverage_plus=1	flanking_left=151	flanking_right=151	frequency=4.166e-03	junction_possible_overlap_registers=145	key=TMV-plasmid__7657__1__TMV-plasmid__9270__1__2____151__151__0__0	max_left=148	max_left_minus=148	max_left_plus=77	max_min_left=0	max_min_left_minus=0	max_min_left_plus=0	max_min_right=72	max_min_right_minus=54	max_min_right_plus=72	max_pos_hash_score=292	max_right=72	max_right_minus=54	max_right_plus=72	neg_log10_pos_hash_p_value=5.8	new_junction_coverage=0.01	new_junction_read_count=3	polymorphism_frequency=4.166e-03	pos_hash_score=3	prediction=polymorphism	reject=COVERAGE_EVENNESS_SKEW,FREQUENCY_CUTOFF	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=3.17	side_1_gene_name=–/–	side_1_gene_position=intergenic (–/–)	side_1_gene_product=–/–	side_1_gene_strand=–/–	side_1_locus_tag=–/–	side_1_overlap=2	side_1_possible_overlap_registers=147	side_1_read_count=654	side_1_redundant=0	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=3.87	side_2_gene_name=–/–	side_2_gene_position=intergenic (–/–)	side_2_gene_product=–/–	side_2_gene_strand=–/–	side_2_locus_tag=–/–	side_2_overlap=0	side_2_possible_overlap_registers=145	side_2_read_count=789	side_2_redundant=0	total_non_overlap_reads=3
JC	15	.	TMV-plasmid	7707	1	TMV-plasmid	8742	1	0	alignment_overlap=0	coverage_minus=4	coverage_plus=2	flanking_left=151	flanking_right=151	frequency=1.070e-02	junction_possible_overlap_registers=147	key=TMV-plasmid__7707__1__TMV-plasmid__8742__1__0____151__151__0__0	max_left=150	max_left_minus=150	max_left_plus=59	max_min_left=59	max_min_left_minus=0	max_min_left_plus=59	max_min_right=57	max_min_right_minus=57	max_min_right_plus=0	max_pos_hash_score=296	max_right=92	max_right_minus=57	max_right_plus=92	neg_log10_pos_hash_p_value=5.4	new_junction_coverage=0.03	new_junction_read_count=7	polymorphism_frequency=1.070e-02	pos_hash_score=4	prediction=polymorphism	reject=COVERAGE_EVENNESS_SKEW,FREQUENCY_CUTOFF	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=3.14	side_1_gene_name=–/–	side_1_gene_position=intergenic (–/–)	side_1_gene_product=–/–	side_1_gene_strand=–/–	side_1_locus_tag=–/–	side_1_overlap=0	side_1_possible_overlap_registers=147	side_1_read_count=649	side_1_redundant=0	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=3.12	side_2_gene_name=–/–	side_2_gene_position=intergenic (–/–)	side_2_gene_product=–/–	side_2_gene_strand=–/–	side_2_locus_tag=–/–	side_2_overlap=0	side_2_possible_overlap_registers=147	side_2_read_count=645	side_2_redundant=0	total_non_overlap_reads=6
JC	16	.	TMV-plasmid	7709	1	TMV-plasmid	9354	1	0	alignment_overlap=0	coverage_minus=2	coverage_plus=3	flanking_left=151	flanking_right=151	frequency=6.562e-03	junction_possible_overlap_registers=147	key=TMV-plasmid__7709__1__TMV-plasmid__9354__1__0____151__151__0__0	max_left=87	max_left_minus=87	max_left_plus=46	max_min_left=46	max_min_left_minus=0	max_min_left_plus=46	max_min_right=64	max_min_right_minus=64	max_min_right_plus=0	max_pos_hash_score=296	max_right=150	max_right_minus=64	max_right_plus=150	neg_log10_pos_hash_p_value=5.9	new_junction_coverage=0.02	new_junction_read_count=5	no_show=1	polymorphism_frequency=6.562e-03	pos_hash_score=3	prediction=polymorphism	reject=COVERAGE_EVENNESS_SKEW,FREQUENCY_CUTOFF	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=3.23	side_1_overlap=0	side_1_possible_overlap_registers=147	side_1_read_count=668	side_1_redundant=0	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=4.09	side_2_overlap=0	side_2_possible_overlap_registers=147	side_2_read_count=846	side_2_redundant=0	total_non_overlap_reads=5
JC	17	.	TMV-plasmid	7729	1	TMV-plasmid	9723	1	0	alignment_overlap=1	coverage_minus=3	coverage_plus=1	flanking_left=151	flanking_right=151	frequency=5.339e-03	junction_possible_overlap_registers=146	key=TMV-plasmid__7729__1__TMV-plasmid__9722__1__1____151__151__0__0	max_left=149	max_left_minus=149	max_left_plus=60	max_min_left=61	max_min_left_minus=61	max_min_left_plus=60	max_min_right=1	max_min_right_minus=1	max_min_right_plus=0	max_pos_hash_score=294	max_right=90	max_right_minus=89	max_right_plus=90	neg_log10_pos_hash_p_value=5.9	new_junction_coverage=0.02	new_junction_read_count=4	no_show=1	polymorphism_frequency=5.339e-03	pos_hash_score=3	prediction=polymorphism	reject=COVERAGE_EVENNESS_SKEW,FREQUENCY_CUTOFF	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=3.32	side_1_overlap=1	side_1_possible_overlap_registers=147	side_1_read_count=686	side_1_redundant=0	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=3.94	side_2_overlap=0	side_2_possible_overlap_registers=146	side_2_read_count=809	side_2_redundant=0	total_non_overlap_reads=4
JC	18	.	TMV-plasmid	7846	-1	TMV-plasmid	9812	1	0	alignment_overlap=0	coverage_minus=1	coverage_plus=2	flanking_left=151	flanking_right=151	frequency=3.630e-03	junction_possible_overlap_registers=147	key=TMV-plasmid__7846__-1__TMV-plasmid__9812__1__0____151__151__0__0	max_left=89	max_left_minus=89	max_left_plus=36	max_min_left=36	max_min_left_minus=0	max_min_left_plus=36	max_min_right=62	max_min_right_minus=62	max_min_right_plus=0	max_pos_hash_score=296	max_right=150	max_right_minus=62	max_right_plus=150	neg_log10_pos_hash_p_value=5.9	new_junction_coverage=0.01	new_junction_read_count=3	no_show=1	polymorphism_frequency=3.630e-03	pos_hash_score=3	prediction=polymorphism	reject=COVERAGE_EVENNESS_SKEW,FREQUENCY_CUTOFF	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=3.89	side_1_overlap=0	side_1_possible_overlap_registers=147	side_1_read_count=804	side_1_redundant=0	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=4.08	side_2_overlap=0	side_2_possible_overlap_registers=147	side_2_read_count=843	side_2_redundant=0	total_non_overlap_reads=3
JC	19	.	TMV-plasmid	7847	1	TMV-plasmid	8069	-1	0	alignment_overlap=0	coverage_minus=1	coverage_plus=3	flanking_left=151	flanking_right=151	frequency=5.145e-03	junction_possible_overlap_registers=147	key=TMV-plasmid__7847__1__TMV-plasmid__8069__-1__0____151__151__0__0	max_left=109	max_left_minus=109	max_left_plus=51	max_min_left=51	max_min_left_minus=0	max_min_left_plus=51	max_min_right=42	max_min_right_minus=42	max_min_right_plus=0	max_pos_hash_score=296	max_right=150	max_right_minus=42	max_right_plus=150	neg_log10_pos_hash_p_value=5.9	new_junction_coverage=0.02	new_junction_read_count=4	no_show=1	polymorphism_frequency=5.145e-03	pos_hash_score=3	prediction=polymorphism	reject=COVERAGE_EVENNESS_SKEW,FREQUENCY_CUTOFF	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=3.89	side_1_overlap=0	side_1_possible_overlap_registers=147	side_1_read_count=804	side_1_redundant=0	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=3.60	side_2_overlap=0	side_2_possible_overlap_registers=147	side_2_read_count=743	side_2_redundant=0	total_non_overlap_reads=4
JC	20	.	TMV-plasmid	7958	1	TMV-plasmid	8183	-1	0	alignment_overlap=0	coverage_minus=1	coverage_plus=2	flanking_left=151	flanking_right=151	frequency=4.405e-03	junction_possible_overlap_registers=147	key=TMV-plasmid__7958__1__TMV-plasmid__8183__-1__0____151__151__0__0	max_left=65	max_left_minus=65	max_left_plus=39	max_min_left=65	max_min_left_minus=65	max_min_left_plus=39	max_min_right=0	max_min_right_minus=0	max_min_right_plus=0	max_pos_hash_score=296	max_right=149	max_right_minus=85	max_right_plus=149	neg_log10_pos_hash_p_value=5.9	new_junction_coverage=0.01	new_junction_read_count=3	no_show=1	polymorphism_frequency=4.405e-03	pos_hash_score=3	prediction=polymorphism	reject=COVERAGE_EVENNESS_SKEW,FREQUENCY_CUTOFF	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=3.26	side_1_overlap=0	side_1_possible_overlap_registers=147	side_1_read_count=673	side_1_redundant=0	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=3.31	side_2_overlap=0	side_2_possible_overlap_registers=147	side_2_read_count=683	side_2_redundant=0	total_non_overlap_reads=3
JC	21	.	TMV-plasmid	8013	1	TMV-plasmid	8947	-1	0	alignment_overlap=7	coverage_minus=2	coverage_plus=2	flanking_left=151	flanking_right=151	frequency=6.513e-03	junction_possible_overlap_registers=140	key=TMV-plasmid__8013__1__TMV-plasmid__8954__-1__7____151__151__0__0	max_left=143	max_left_minus=143	max_left_plus=43	max_min_left=43	max_min_left_minus=0	max_min_left_plus=43	max_min_right=50	max_min_right_minus=50	max_min_right_plus=0	max_pos_hash_score=282	max_right=141	max_right_minus=50	max_right_plus=141	neg_log10_pos_hash_p_value=5.2	new_junction_coverage=0.02	new_junction_read_count=4	polymorphism_frequency=6.513e-03	pos_hash_score=4	prediction=polymorphism	reject=COVERAGE_EVENNESS_SKEW,FREQUENCY_CUTOFF	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=3.32	side_1_gene_name=–/–	side_1_gene_position=intergenic (–/–)	side_1_gene_product=–/–	side_1_gene_strand=–/–	side_1_locus_tag=–/–	side_1_overlap=7	side_1_possible_overlap_registers=147	side_1_read_count=686	side_1_redundant=0	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=2.88	side_2_gene_name=–/–	side_2_gene_position=intergenic (–/–)	side_2_gene_product=–/–	side_2_gene_strand=–/–	side_2_locus_tag=–/–	side_2_overlap=0	side_2_possible_overlap_registers=140	side_2_read_count=567	side_2_redundant=0	total_non_overlap_reads=4
JC	22	.	TMV-plasmid	8112	-1	TMV-plasmid	9262	1	0	alignment_overlap=3	coverage_minus=2	coverage_plus=4	flanking_left=151	flanking_right=151	frequency=7.599e-03	junction_possible_overlap_registers=144	key=TMV-plasmid__8112__-1__TMV-plasmid__9259__1__3____151__151__0__0	max_left=101	max_left_minus=101	max_left_plus=86	max_min_left=1	max_min_left_minus=0	max_min_left_plus=1	max_min_right=61	max_min_right_minus=47	max_min_right_plus=61	max_pos_hash_score=290	max_right=147	max_right_minus=47	max_right_plus=147	neg_log10_pos_hash_p_value=5.8	new_junction_coverage=0.03	new_junction_read_count=6	polymorphism_frequency=7.599e-03	pos_hash_score=3	prediction=polymorphism	reject=COVERAGE_EVENNESS_SKEW,FREQUENCY_CUTOFF	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=3.78	side_1_gene_name=–/–	side_1_gene_position=intergenic (–/–)	side_1_gene_product=–/–	side_1_gene_strand=–/–	side_1_locus_tag=–/–	side_1_overlap=3	side_1_possible_overlap_registers=147	side_1_read_count=781	side_1_redundant=0	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=3.96	side_2_gene_name=–/–	side_2_gene_position=intergenic (–/–)	side_2_gene_product=–/–	side_2_gene_strand=–/–	side_2_locus_tag=–/–	side_2_overlap=0	side_2_possible_overlap_registers=144	side_2_read_count=802	side_2_redundant=0	total_non_overlap_reads=6
JC	23	.	TMV-plasmid	8127	1	TMV-plasmid	9785	-1	0	alignment_overlap=1	coverage_minus=2	coverage_plus=1	flanking_left=151	flanking_right=151	frequency=3.950e-03	junction_possible_overlap_registers=146	key=TMV-plasmid__8127__1__TMV-plasmid__9786__-1__1____151__151__0__0	max_left=149	max_left_minus=149	max_left_plus=101	max_min_left=0	max_min_left_minus=0	max_min_left_plus=0	max_min_right=49	max_min_right_minus=46	max_min_right_plus=49	max_pos_hash_score=294	max_right=49	max_right_minus=46	max_right_plus=49	neg_log10_pos_hash_p_value=5.9	new_junction_coverage=0.01	new_junction_read_count=3	no_show=1	polymorphism_frequency=3.950e-03	pos_hash_score=3	prediction=polymorphism	reject=COVERAGE_EVENNESS_SKEW,FREQUENCY_CUTOFF	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=3.60	side_1_overlap=1	side_1_possible_overlap_registers=147	side_1_read_count=743	side_1_redundant=0	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=3.78	side_2_overlap=0	side_2_possible_overlap_registers=146	side_2_read_count=775	side_2_redundant=0	total_non_overlap_reads=3
JC	24	.	TMV-plasmid	8626	1	TMV-plasmid	8784	1	0	alignment_overlap=0	coverage_minus=1	coverage_plus=2	flanking_left=151	flanking_right=151	frequency=4.785e-03	junction_possible_overlap_registers=147	key=TMV-plasmid__8626__1__TMV-plasmid__8784__1__0____151__151__0__0	max_left=102	max_left_minus=102	max_left_plus=102	max_min_left=1	max_min_left_minus=0	max_min_left_plus=1	max_min_right=46	max_min_right_minus=41	max_min_right_plus=46	max_pos_hash_score=296	max_right=150	max_right_minus=41	max_right_plus=150	neg_log10_pos_hash_p_value=5.9	new_junction_coverage=0.01	new_junction_read_count=3	no_show=1	polymorphism_frequency=4.785e-03	pos_hash_score=3	prediction=polymorphism	reject=COVERAGE_EVENNESS_SKEW,FREQUENCY_CUTOFF	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=3.09	side_1_overlap=0	side_1_possible_overlap_registers=147	side_1_read_count=638	side_1_redundant=0	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=2.95	side_2_overlap=0	side_2_possible_overlap_registers=147	side_2_read_count=610	side_2_redundant=0	total_non_overlap_reads=3
JC	25	.	TMV-plasmid	8658	1	TMV-plasmid	9800	1	0	alignment_overlap=2	coverage_minus=1	coverage_plus=2	flanking_left=151	flanking_right=151	frequency=4.104e-03	junction_possible_overlap_registers=145	key=TMV-plasmid__8658__1__TMV-plasmid__9798__1__2____151__151__0__0	max_left=89	max_left_minus=89	max_left_plus=35	max_min_left=35	max_min_left_minus=0	max_min_left_plus=35	max_min_right=60	max_min_right_minus=60	max_min_right_plus=0	max_pos_hash_score=292	max_right=143	max_right_minus=60	max_right_plus=143	neg_log10_pos_hash_p_value=5.8	new_junction_coverage=0.01	new_junction_read_count=3	polymorphism_frequency=4.104e-03	pos_hash_score=3	prediction=polymorphism	reject=COVERAGE_EVENNESS_SKEW,FREQUENCY_CUTOFF	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=3.21	side_1_gene_name=–/–	side_1_gene_position=intergenic (–/–)	side_1_gene_product=–/–	side_1_gene_strand=–/–	side_1_locus_tag=–/–	side_1_overlap=2	side_1_possible_overlap_registers=147	side_1_read_count=664	side_1_redundant=0	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=3.93	side_2_gene_name=–/–	side_2_gene_position=intergenic (–/–)	side_2_gene_product=–/–	side_2_gene_strand=–/–	side_2_locus_tag=–/–	side_2_overlap=0	side_2_possible_overlap_registers=145	side_2_read_count=801	side_2_redundant=0	total_non_overlap_reads=3
JC	26	.	TMV-plasmid	8694	1	TMV-plasmid	9377	1	0	alignment_overlap=1	coverage_minus=1	coverage_plus=3	flanking_left=151	flanking_right=151	frequency=5.197e-03	junction_possible_overlap_registers=146	key=TMV-plasmid__8694__1__TMV-plasmid__9376__1__1____151__151__0__0	max_left=60	max_left_minus=60	max_left_plus=53	max_min_left=60	max_min_left_minus=60	max_min_left_plus=53	max_min_right=0	max_min_right_minus=0	max_min_right_plus=0	max_pos_hash_score=294	max_right=149	max_right_minus=90	max_right_plus=149	neg_log10_pos_hash_p_value=5.9	new_junction_coverage=0.02	new_junction_read_count=4	no_show=1	polymorphism_frequency=5.197e-03	pos_hash_score=3	prediction=polymorphism	reject=COVERAGE_EVENNESS_SKEW,FREQUENCY_CUTOFF	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=3.27	side_1_overlap=1	side_1_possible_overlap_registers=147	side_1_read_count=675	side_1_redundant=0	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=4.20	side_2_overlap=0	side_2_possible_overlap_registers=146	side_2_read_count=861	side_2_redundant=0	total_non_overlap_reads=4
JC	27	.	TMV-plasmid	8991	1	TMV-plasmid	9092	-1	0	alignment_overlap=0	coverage_minus=2	coverage_plus=1	flanking_left=151	flanking_right=151	frequency=3.856e-03	junction_possible_overlap_registers=147	key=TMV-plasmid__8991__1__TMV-plasmid__9092__-1__0____151__151__0__0	max_left=148	max_left_minus=148	max_left_plus=77	max_min_left=0	max_min_left_minus=0	max_min_left_plus=0	max_min_right=73	max_min_right_minus=54	max_min_right_plus=73	max_pos_hash_score=296	max_right=73	max_right_minus=54	max_right_plus=73	neg_log10_pos_hash_p_value=5.9	new_junction_coverage=0.01	new_junction_read_count=3	no_show=1	polymorphism_frequency=3.856e-03	pos_hash_score=3	prediction=polymorphism	reject=COVERAGE_EVENNESS_SKEW,FREQUENCY_CUTOFF	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=3.47	side_1_overlap=0	side_1_possible_overlap_registers=147	side_1_read_count=717	side_1_redundant=0	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=4.03	side_2_overlap=0	side_2_possible_overlap_registers=147	side_2_read_count=833	side_2_redundant=0	total_non_overlap_reads=3
JC	28	.	TMV-plasmid	9280	-1	TMV-plasmid	9388	1	0	alignment_overlap=0	coverage_minus=3	coverage_plus=2	flanking_left=151	flanking_right=151	frequency=5.893e-03	junction_possible_overlap_registers=147	key=TMV-plasmid__9280__-1__TMV-plasmid__9388__1__0____151__151__0__0	max_left=146	max_left_minus=146	max_left_plus=37	max_min_left=37	max_min_left_minus=0	max_min_left_plus=37	max_min_right=57	max_min_right_minus=57	max_min_right_plus=0	max_pos_hash_score=296	max_right=113	max_right_minus=57	max_right_plus=113	neg_log10_pos_hash_p_value=5.9	new_junction_coverage=0.02	new_junction_read_count=5	no_show=1	polymorphism_frequency=5.893e-03	pos_hash_score=3	prediction=polymorphism	reject=COVERAGE_EVENNESS_SKEW,FREQUENCY_CUTOFF	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=3.83	side_1_overlap=0	side_1_possible_overlap_registers=147	side_1_read_count=791	side_1_redundant=0	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=4.34	side_2_overlap=0	side_2_possible_overlap_registers=147	side_2_read_count=896	side_2_redundant=0	total_non_overlap_reads=5
JC	29	.	TMV-plasmid	9443	1	TMV-plasmid	9519	1	0	alignment_overlap=8	coverage_minus=1	coverage_plus=4	flanking_left=151	flanking_right=151	frequency=6.766e-03	junction_possible_overlap_registers=139	key=TMV-plasmid__9443__1__TMV-plasmid__9511__1__8____151__151__0__0	max_left=83	max_left_minus=71	max_left_plus=83	max_min_left=71	max_min_left_minus=71	max_min_left_plus=61	max_min_right=56	max_min_right_minus=0	max_min_right_plus=56	max_pos_hash_score=280	max_right=140	max_right_minus=72	max_right_plus=140	neg_log10_pos_hash_p_value=4.7	new_junction_coverage=0.03	new_junction_read_count=6	polymorphism_frequency=6.766e-03	pos_hash_score=5	prediction=polymorphism	reject=COVERAGE_EVENNESS_SKEW,FREQUENCY_CUTOFF	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=4.40	side_1_gene_name=–/–	side_1_gene_position=intergenic (–/–)	side_1_gene_product=–/–	side_1_gene_strand=–/–	side_1_locus_tag=–/–	side_1_overlap=8	side_1_possible_overlap_registers=147	side_1_read_count=909	side_1_redundant=0	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=4.62	side_2_gene_name=–/–	side_2_gene_position=intergenic (–/–)	side_2_gene_product=–/–	side_2_gene_strand=–/–	side_2_locus_tag=–/–	side_2_overlap=0	side_2_possible_overlap_registers=139	side_2_read_count=902	side_2_redundant=0	total_non_overlap_reads=5
JC	30	.	TMV-plasmid	9544	-1	TMV-plasmid	9840	1	0	alignment_overlap=3	coverage_minus=2	coverage_plus=1	flanking_left=151	flanking_right=151	frequency=3.272e-03	junction_possible_overlap_registers=144	key=TMV-plasmid__9544__-1__TMV-plasmid__9837__1__3____151__151__0__0	max_left=84	max_left_minus=84	max_left_plus=43	max_min_left=51	max_min_left_minus=51	max_min_left_plus=43	max_min_right=71	max_min_right_minus=71	max_min_right_plus=0	max_pos_hash_score=290	max_right=105	max_right_minus=97	max_right_plus=105	neg_log10_pos_hash_p_value=5.8	new_junction_coverage=0.01	new_junction_read_count=3	polymorphism_frequency=3.272e-03	pos_hash_score=3	prediction=polymorphism	reject=COVERAGE_EVENNESS_SKEW,FREQUENCY_CUTOFF	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=4.79	side_1_gene_name=–/–	side_1_gene_position=intergenic (–/–)	side_1_gene_product=–/–	side_1_gene_strand=–/–	side_1_locus_tag=–/–	side_1_overlap=3	side_1_possible_overlap_registers=147	side_1_read_count=990	side_1_redundant=0	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=4.24	side_2_gene_name=–/–	side_2_gene_position=intergenic (–/–)	side_2_gene_product=–/–	side_2_gene_strand=–/–	side_2_locus_tag=–/–	side_2_overlap=0	side_2_possible_overlap_registers=144	side_2_read_count=858	side_2_redundant=0	total_non_overlap_reads=3
JC	31	.	TMV-plasmid	9676	1	TMV-plasmid	9753	1	0	alignment_overlap=2	coverage_minus=2	coverage_plus=1	flanking_left=151	flanking_right=151	frequency=3.734e-03	junction_possible_overlap_registers=145	key=TMV-plasmid__9676__1__TMV-plasmid__9751__1__2____151__151__0__0	max_left=148	max_left_minus=148	max_left_plus=65	max_min_left=65	max_min_left_minus=0	max_min_left_plus=65	max_min_right=61	max_min_right_minus=61	max_min_right_plus=0	max_pos_hash_score=292	max_right=84	max_right_minus=61	max_right_plus=84	neg_log10_pos_hash_p_value=5.8	new_junction_coverage=0.01	new_junction_read_count=3	polymorphism_frequency=3.734e-03	pos_hash_score=3	prediction=polymorphism	reject=COVERAGE_EVENNESS_SKEW,FREQUENCY_CUTOFF	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=3.95	side_1_gene_name=–/–	side_1_gene_position=intergenic (–/–)	side_1_gene_product=–/–	side_1_gene_strand=–/–	side_1_locus_tag=–/–	side_1_overlap=2	side_1_possible_overlap_registers=147	side_1_read_count=817	side_1_redundant=0	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=3.90	side_2_gene_name=–/–	side_2_gene_position=intergenic (–/–)	side_2_gene_product=–/–	side_2_gene_strand=–/–	side_2_locus_tag=–/–	side_2_overlap=0	side_2_possible_overlap_registers=145	side_2_read_count=795	side_2_redundant=0	total_non_overlap_reads=3
UN	32	.	TMV-plasmid	1	185
UN	33	.	TMV-plasmid	191	191
UN	34	.	TMV-plasmid	210	211
UN	35	.	TMV-plasmid	213	213
UN	36	.	TMV-plasmid	219	219
UN	37	.	TMV-plasmid	231	231
UN	38	.	TMV-plasmid	238	238
UN	39	.	TMV-plasmid	241	241
UN	40	.	TMV-plasmid	268	268
UN	41	.	TMV-plasmid	274	274
UN	42	.	TMV-plasmid	283	283
UN	43	.	TMV-plasmid	294	294
UN	44	.	TMV-plasmid	298	298
UN	45	.	TMV-plasmid	301	301
UN	46	.	TMV-plasmid	304	304
UN	47	.	TMV-plasmid	306	772
UN	48	.	TMV-plasmid	779	1544
UN	49	.	TMV-plasmid	1564	1564
UN	50	.	TMV-plasmid	1566	1567
UN	51	.	TMV-plasmid	1569	1571
UN	52	.	TMV-plasmid	1576	1576
UN	53	.	TMV-plasmid	1579	1579
UN	54	.	TMV-plasmid	1582	1583
UN	55	.	TMV-plasmid	1586	1586
UN	56	.	TMV-plasmid	1590	1592
UN	57	.	TMV-plasmid	1598	1598
UN	58	.	TMV-plasmid	1600	1754
UN	59	.	TMV-plasmid	1764	1764
UN	60	.	TMV-plasmid	1774	1774
UN	61	.	TMV-plasmid	1777	1777
UN	62	.	TMV-plasmid	1779	1779
UN	63	.	TMV-plasmid	1781	1781
UN	64	.	TMV-plasmid	1785	1785
UN	65	.	TMV-plasmid	1958	1958
UN	66	.	TMV-plasmid	1964	1964
UN	67	.	TMV-plasmid	1969	2427
UN	68	.	TMV-plasmid	2430	2430
UN	69	.	TMV-plasmid	2441	2441
UN	70	.	TMV-plasmid	2540	2540
UN	71	.	TMV-plasmid	2547	2548
UN	72	.	TMV-plasmid	2552	2552
UN	73	.	TMV-plasmid	2563	2563
UN	74	.	TMV-plasmid	2570	2879
UN	75	.	TMV-plasmid	2882	2882
UN	76	.	TMV-plasmid	2888	2888
UN	77	.	TMV-plasmid	2935	2936
UN	78	.	TMV-plasmid	2938	2938
UN	79	.	TMV-plasmid	2944	2944
UN	80	.	TMV-plasmid	2952	2956
UN	81	.	TMV-plasmid	2989	2989
UN	82	.	TMV-plasmid	3024	3031
UN	83	.	TMV-plasmid	3037	3037
UN	84	.	TMV-plasmid	3044	3044
UN	85	.	TMV-plasmid	3056	3056
UN	86	.	TMV-plasmid	3062	3062
UN	87	.	TMV-plasmid	3070	3070
UN	88	.	TMV-plasmid	3072	3072
UN	89	.	TMV-plasmid	3076	3076
UN	90	.	TMV-plasmid	3109	3109
UN	91	.	TMV-plasmid	3181	3181
UN	92	.	TMV-plasmid	3183	3189
UN	93	.	TMV-plasmid	3194	3195
UN	94	.	TMV-plasmid	3199	3199
UN	95	.	TMV-plasmid	3202	3202
UN	96	.	TMV-plasmid	3204	3204
UN	97	.	TMV-plasmid	3206	3206
UN	98	.	TMV-plasmid	3213	3214
UN	99	.	TMV-plasmid	3216	3216
UN	100	.	TMV-plasmid	3218	3219
UN	101	.	TMV-plasmid	3222	3222
UN	102	.	TMV-plasmid	3224	3225
UN	103	.	TMV-plasmid	3228	3229
UN	104	.	TMV-plasmid	3235	3235
UN	105	.	TMV-plasmid	3239	3292
UN	106	.	TMV-plasmid	3298	3322
UN	107	.	TMV-plasmid	3408	3408
UN	108	.	TMV-plasmid	3416	3416
UN	109	.	TMV-plasmid	3438	3576
UN	110	.	TMV-plasmid	3582	3582
UN	111	.	TMV-plasmid	3739	3739
UN	112	.	TMV-plasmid	3758	3758
UN	113	.	TMV-plasmid	3763	3802
UN	114	.	TMV-plasmid	3837	3837
UN	115	.	TMV-plasmid	3842	3842
UN	116	.	TMV-plasmid	3846	3846
UN	117	.	TMV-plasmid	3864	3864
UN	118	.	TMV-plasmid	3866	3877
UN	119	.	TMV-plasmid	3879	3883
UN	120	.	TMV-plasmid	3887	3887
UN	121	.	TMV-plasmid	3891	3892
UN	122	.	TMV-plasmid	3894	3894
UN	123	.	TMV-plasmid	3896	3896
UN	124	.	TMV-plasmid	3900	3900
UN	125	.	TMV-plasmid	3903	3904
UN	126	.	TMV-plasmid	3906	3907
UN	127	.	TMV-plasmid	3909	3913
UN	128	.	TMV-plasmid	3915	3915
UN	129	.	TMV-plasmid	3917	3921
UN	130	.	TMV-plasmid	3926	3927
UN	131	.	TMV-plasmid	3929	3930
UN	132	.	TMV-plasmid	3932	3932
UN	133	.	TMV-plasmid	3934	3934
UN	134	.	TMV-plasmid	3937	3938
UN	135	.	TMV-plasmid	3940	4236
UN	136	.	TMV-plasmid	4239	4241
UN	137	.	TMV-plasmid	4391	4391
UN	138	.	TMV-plasmid	4393	4393
UN	139	.	TMV-plasmid	4397	4398
UN	140	.	TMV-plasmid	4401	4401
UN	141	.	TMV-plasmid	4403	4407
UN	142	.	TMV-plasmid	4409	4409
UN	143	.	TMV-plasmid	4416	4416
UN	144	.	TMV-plasmid	4418	4432
UN	145	.	TMV-plasmid	4434	4434
UN	146	.	TMV-plasmid	4436	4682
UN	147	.	TMV-plasmid	4686	4686
UN	148	.	TMV-plasmid	4831	4831
UN	149	.	TMV-plasmid	4846	4846
UN	150	.	TMV-plasmid	4850	5173
UN	151	.	TMV-plasmid	5299	5299
UN	152	.	TMV-plasmid	5303	5400
UN	153	.	TMV-plasmid	5406	5406
UN	154	.	TMV-plasmid	5483	5483
UN	155	.	TMV-plasmid	5493	5493
UN	156	.	TMV-plasmid	5509	5542
UN	157	.	TMV-plasmid	5545	5551
UN	158	.	TMV-plasmid	5557	5558
UN	159	.	TMV-plasmid	5562	5562
UN	160	.	TMV-plasmid	5564	5564
UN	161	.	TMV-plasmid	5586	5586
UN	162	.	TMV-plasmid	5595	5595
UN	163	.	TMV-plasmid	5598	5598
UN	164	.	TMV-plasmid	5625	5625
UN	165	.	TMV-plasmid	5634	5634
UN	166	.	TMV-plasmid	5643	5643
UN	167	.	TMV-plasmid	5673	5673
UN	168	.	TMV-plasmid	5677	5677
UN	169	.	TMV-plasmid	5688	6662
UN	170	.	TMV-plasmid	6800	6800
UN	171	.	TMV-plasmid	6802	6802
UN	172	.	TMV-plasmid	6804	6888
UN	173	.	TMV-plasmid	10254	10411
