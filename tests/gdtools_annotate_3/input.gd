#=GENOME_DIFF	1.0
#=CREATED	23:54:33 19 Jan 2020
#=PROGRAM	breseq 0.35.0 revision fe39cd8f5880
#=COMMAND	./src/c/breseq/breseq -j 4 -o tests/lambda_mult_ref_read -r tests/lambda_mult_ref_read/../data/lambda/lambda.1-2.gbk -r tests/lambda_mult_ref_read/../data/lambda/lambda.3.gbk -r tests/lambda_mult_ref_read/../data/lambda/lambda.4.gbk -r tests/lambda_mult_ref_read/../data/lambda/lambda.5.gbk -l 50 tests/lambda_mult_ref_read/../data/lambda/empty.fastq tests/lambda_mult_ref_read/../data/lambda/only_bad.fastq tests/lambda_mult_ref_read/../data/lambda/lambda_mixed_population.1.fastq tests/lambda_mult_ref_read/../data/lambda/lambda_mixed_population.2.fastq tests/lambda_mult_ref_read/../data/lambda/lambda_mixed_population.3.fastq tests/lambda_mult_ref_read/../data/lambda/lambda_mixed_population.4.fastq tests/lambda_mult_ref_read/../data/lambda/lambda_mixed_population.5.fastq
#=REFSEQ	tests/lambda_mult_ref_read/../data/lambda/lambda.1-2.gbk
#=REFSEQ	tests/lambda_mult_ref_read/../data/lambda/lambda.3.gbk
#=REFSEQ	tests/lambda_mult_ref_read/../data/lambda/lambda.4.gbk
#=REFSEQ	tests/lambda_mult_ref_read/../data/lambda/lambda.5.gbk
#=READSEQ	tests/lambda_mult_ref_read/../data/lambda/empty.fastq
#=READSEQ	tests/lambda_mult_ref_read/../data/lambda/only_bad.fastq
#=READSEQ	tests/lambda_mult_ref_read/../data/lambda/lambda_mixed_population.1.fastq
#=READSEQ	tests/lambda_mult_ref_read/../data/lambda/lambda_mixed_population.2.fastq
#=CONVERTED-BASES	2425115
#=CONVERTED-READS	69289
#=INPUT-BASES	2800140
#=INPUT-READS	80004
#=MAPPED-BASES	1783840
#=MAPPED-READS	51266
DEL	1	29	NC_001416-0	139	1
INS	2	30	NC_001416-1	4566	G
SNP	3	31	NC_001416-2	1261	G
INS	4	32	NC_001416-2	1435	C
SNP	5	33	NC_001416-2	2314	A
DEL	6	60,65	NC_001416-2	2338	5996
SNP	7	34	NC_001416-3	1915	C
SNP	8	35	NC_001416-3	5833	G
DEL	9	36	NC_001416-3	8717	1
SNP	10	37	NC_001416-4	6817	C
INS	11	39	NC_001416-4	8156	A
SNP	12	40	NC_001416-4	8184	T
SNP	13	41	NC_001416-4	8191	T
SNP	14	42	NC_001416-4	8203	A
SNP	15	43	NC_001416-4	8328	G
SNP	16	44	NC_001416-4	8342	T
SNP	17	45	NC_001416-4	8442	A
SNP	18	46	NC_001416-4	8514	A
SNP	19	47	NC_001416-4	8559	A
SNP	20	48	NC_001416-4	8597	T
SNP	21	49	NC_001416-4	8708	C
SNP	22	50	NC_001416-4	8728	T
SNP	23	51	NC_001416-4	8774	A
SNP	24	52	NC_001416-4	8868	C
SNP	25	53	NC_001416-4	9077	G
SNP	26	54	NC_001416-4	9172	C
SUB	27	55,56	NC_001416-4	9176	2	AC
SNP	28	57	NC_001416-4	9359	C
RA	29	.	NC_001416-0	138	0	G	.	consensus_score=66.5	frequency=1	major_base=.	major_cov=8/11	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=8/11	polymorphism_frequency=1.000e+00	polymorphism_score=NA	prediction=consensus	ref_cov=0/0	total_cov=8/11
RA	30	.	NC_001416-1	4566	1	.	G	consensus_score=66.6	frequency=1	major_base=G	major_cov=13/13	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=13/13	polymorphism_frequency=1.000e+00	polymorphism_score=NA	prediction=consensus	ref_cov=0/0	total_cov=13/13
RA	31	.	NC_001416-2	1261	0	A	G	consensus_score=66.5	frequency=1	major_base=G	major_cov=8/18	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=8/18	polymorphism_frequency=1.000e+00	polymorphism_score=NA	prediction=consensus	ref_cov=0/0	total_cov=8/18
RA	32	.	NC_001416-2	1432	1	.	C	consensus_score=84.6	frequency=1	major_base=C	major_cov=9/21	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=9/21	polymorphism_frequency=1.000e+00	polymorphism_score=NA	prediction=consensus	ref_cov=0/0	total_cov=9/21
RA	33	.	NC_001416-2	2314	0	G	A	bias_e_value=48501.9	bias_p_value=0.999998	consensus_score=90.2	fisher_strand_p_value=1	frequency=1	ks_quality_p_value=0.998219	major_base=A	major_cov=12/21	major_frequency=9.706e-01	minor_base=G	minor_cov=0/1	new_cov=12/21	polymorphism_frequency=9.706e-01	polymorphism_score=-4.4	prediction=consensus	ref_cov=0/1	total_cov=13/22
RA	34	.	NC_001416-3	1915	0	T	C	bias_e_value=48502	bias_p_value=1	consensus_score=108.4	fisher_strand_p_value=1	frequency=1	ks_quality_p_value=1	major_base=C	major_cov=25/17	major_frequency=9.767e-01	minor_base=A	minor_cov=1/0	new_cov=25/17	polymorphism_frequency=9.767e-01	polymorphism_score=NA	prediction=consensus	ref_cov=0/0	total_cov=26/17
RA	35	.	NC_001416-3	5833	0	A	G	consensus_score=59.4	frequency=1	major_base=G	major_cov=6/19	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=6/19	polymorphism_frequency=1.000e+00	polymorphism_score=NA	prediction=consensus	ref_cov=0/0	total_cov=6/19
RA	36	.	NC_001416-3	8714	0	C	.	consensus_score=71.5	frequency=1	major_base=.	major_cov=13/7	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=13/7	polymorphism_frequency=1.000e+00	polymorphism_score=NA	prediction=consensus	ref_cov=0/0	total_cov=13/7
RA	37	.	NC_001416-4	6817	0	T	C	bias_e_value=45093	bias_p_value=0.929714	consensus_score=129.4	fisher_strand_p_value=1	frequency=1	ks_quality_p_value=0.649315	major_base=C	major_cov=24/25	major_frequency=9.800e-01	minor_base=T	minor_cov=0/1	new_cov=24/25	polymorphism_frequency=9.800e-01	polymorphism_score=-3.0	prediction=consensus	ref_cov=0/1	total_cov=25/26
RA	38	.	NC_001416-4	7820	0	A	C	bias_e_value=22792.1	bias_p_value=0.469921	consensus_reject=FREQUENCY_CUTOFF	consensus_score=71.9	fisher_strand_p_value=1	frequency=2.105e-01	ks_quality_p_value=0.169261	major_base=A	major_cov=17/13	major_frequency=7.895e-01	minor_base=C	minor_cov=5/3	new_cov=5/3	polymorphism_frequency=2.105e-01	polymorphism_score=12.5	prediction=polymorphism	ref_cov=17/13	total_cov=22/16
RA	39	.	NC_001416-4	8152	1	.	A	bias_e_value=47141.3	bias_p_value=0.971946	consensus_score=50.7	fisher_strand_p_value=1	frequency=1	ks_quality_p_value=0.772677	major_base=A	major_cov=8/11	major_frequency=9.500e-01	minor_base=.	minor_cov=0/1	new_cov=8/11	polymorphism_frequency=9.500e-01	polymorphism_score=-1.1	prediction=consensus	ref_cov=0/1	total_cov=8/12
RA	40	.	NC_001416-4	8184	0	C	T	consensus_score=105.0	frequency=1	major_base=T	major_cov=21/16	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=21/16	polymorphism_frequency=1.000e+00	polymorphism_score=NA	prediction=consensus	ref_cov=0/0	total_cov=21/16
RA	41	.	NC_001416-4	8191	0	C	T	consensus_score=69.9	frequency=1	major_base=T	major_cov=12/14	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=12/14	polymorphism_frequency=1.000e+00	polymorphism_score=NA	prediction=consensus	ref_cov=0/0	total_cov=12/14
RA	42	.	NC_001416-4	8203	0	G	A	consensus_score=77.1	frequency=1	major_base=A	major_cov=12/18	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=12/18	polymorphism_frequency=1.000e+00	polymorphism_score=NA	prediction=consensus	ref_cov=0/0	total_cov=12/18
RA	43	.	NC_001416-4	8328	0	A	G	consensus_score=90.6	frequency=1	major_base=G	major_cov=19/14	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=19/14	polymorphism_frequency=1.000e+00	polymorphism_score=NA	prediction=consensus	ref_cov=0/0	total_cov=19/14
RA	44	.	NC_001416-4	8342	0	C	T	bias_e_value=48502	bias_p_value=1	consensus_score=85.1	fisher_strand_p_value=1	frequency=1	ks_quality_p_value=1	major_base=T	major_cov=18/13	major_frequency=9.688e-01	minor_base=C	minor_cov=1/0	new_cov=18/13	polymorphism_frequency=9.688e-01	polymorphism_score=-3.7	prediction=consensus	ref_cov=1/0	total_cov=19/14
RA	45	.	NC_001416-4	8442	0	G	A	consensus_score=59.0	frequency=1	major_base=A	major_cov=7/15	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=7/15	polymorphism_frequency=1.000e+00	polymorphism_score=NA	prediction=consensus	ref_cov=0/0	total_cov=7/15
RA	46	.	NC_001416-4	8514	0	G	A	consensus_score=119.0	frequency=1	major_base=A	major_cov=23/19	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=23/19	polymorphism_frequency=1.000e+00	polymorphism_score=NA	prediction=consensus	ref_cov=0/0	total_cov=23/19
RA	47	.	NC_001416-4	8559	0	G	A	consensus_score=98.0	frequency=1	major_base=A	major_cov=13/22	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=13/22	polymorphism_frequency=1.000e+00	polymorphism_score=NA	prediction=consensus	ref_cov=0/0	total_cov=13/22
RA	48	.	NC_001416-4	8597	0	C	T	consensus_score=77.3	frequency=1	major_base=T	major_cov=15/13	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=15/13	polymorphism_frequency=1.000e+00	polymorphism_score=NA	prediction=consensus	ref_cov=0/0	total_cov=15/13
RA	49	.	NC_001416-4	8708	0	T	C	consensus_score=51.4	frequency=1	major_base=C	major_cov=12/8	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=12/8	polymorphism_frequency=1.000e+00	polymorphism_score=NA	prediction=consensus	ref_cov=0/0	total_cov=12/8
RA	50	.	NC_001416-4	8728	0	C	T	consensus_score=81.5	frequency=1	major_base=T	major_cov=10/19	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=10/19	polymorphism_frequency=1.000e+00	polymorphism_score=NA	prediction=consensus	ref_cov=0/0	total_cov=10/19
RA	51	.	NC_001416-4	8774	0	C	A	bias_e_value=38186.6	bias_p_value=0.78732	consensus_score=87.2	fisher_strand_p_value=0.424242	frequency=1	ks_quality_p_value=0.998108	major_base=A	major_cov=13/19	major_frequency=9.697e-01	minor_base=C	minor_cov=1/0	new_cov=13/19	polymorphism_frequency=9.697e-01	polymorphism_score=NA	prediction=consensus	ref_cov=1/0	total_cov=14/19
RA	52	.	NC_001416-4	8868	0	T	C	consensus_score=105.8	frequency=1	major_base=C	major_cov=23/17	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=23/17	polymorphism_frequency=1.000e+00	polymorphism_score=NA	prediction=consensus	ref_cov=0/0	total_cov=23/17
RA	53	.	NC_001416-4	9077	0	A	G	consensus_score=80.2	frequency=1	major_base=G	major_cov=13/18	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=13/18	polymorphism_frequency=1.000e+00	polymorphism_score=NA	prediction=consensus	ref_cov=0/0	total_cov=13/18
RA	54	.	NC_001416-4	9172	0	T	C	consensus_score=97.5	frequency=1	major_base=C	major_cov=18/19	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=18/19	polymorphism_frequency=1.000e+00	polymorphism_score=NA	prediction=consensus	ref_cov=0/0	total_cov=18/19
RA	55	.	NC_001416-4	9176	0	G	A	consensus_score=95.5	frequency=1	major_base=A	major_cov=18/17	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=18/17	polymorphism_frequency=1.000e+00	polymorphism_score=NA	prediction=consensus	ref_cov=0/0	total_cov=18/17
RA	56	.	NC_001416-4	9177	0	T	C	bias_e_value=40570.2	bias_p_value=0.836465	consensus_score=91.4	fisher_strand_p_value=0.485714	frequency=1	ks_quality_p_value=1	major_base=C	major_cov=18/16	major_frequency=9.714e-01	minor_base=G	minor_cov=0/1	new_cov=18/16	polymorphism_frequency=9.714e-01	polymorphism_score=-4.3	prediction=consensus	ref_cov=0/0	total_cov=18/17
RA	57	.	NC_001416-4	9359	0	T	C	bias_e_value=42433.8	bias_p_value=0.874889	consensus_score=58.6	fisher_strand_p_value=0.625341	frequency=1	ks_quality_p_value=0.869132	major_base=C	major_cov=12/16	major_frequency=8.485e-01	minor_base=T	minor_cov=1/4	new_cov=12/16	polymorphism_frequency=8.485e-01	polymorphism_score=5.9	prediction=consensus	ref_cov=1/4	total_cov=13/20
RA	58	.	NC_001416-4	9494	0	C	A	bias_e_value=21103.3	bias_p_value=0.435101	consensus_reject=FREQUENCY_CUTOFF	consensus_score=61.3	fisher_strand_p_value=0.25343	frequency=2.250e-01	ks_quality_p_value=0.592988	major_base=C	major_cov=19/12	major_frequency=7.750e-01	minor_base=A	minor_cov=3/6	new_cov=3/6	polymorphism_frequency=2.250e-01	polymorphism_score=18.3	prediction=polymorphism	ref_cov=19/12	total_cov=22/18
MC	59	.	NC_001416-0	1	2	0	0	left_inside_cov=0	left_outside_cov=NA	right_inside_cov=0	right_outside_cov=49
MC	60	.	NC_001416-2	2338	8333	0	0	left_inside_cov=0	left_outside_cov=34	right_inside_cov=1	right_outside_cov=36
MC	61	.	NC_001416-3	1	13	0	0	left_inside_cov=0	left_outside_cov=NA	right_inside_cov=16	right_outside_cov=20
MC	62	.	NC_001416-4	9654	9701	0	0	left_inside_cov=12	left_outside_cov=13	right_inside_cov=0	right_outside_cov=NA
JC	63	.	NC_001416-0	9700	-1	NC_001416-1	1	1	0	alignment_overlap=0	coverage_minus=28	coverage_plus=14	flanking_left=35	flanking_right=35	frequency=1	junction_possible_overlap_registers=34	key=NC_001416-0__9700__-1__NC_001416-1__1__1__0____35__35__0__0	max_left=34	max_left_minus=34	max_left_plus=32	max_min_left=17	max_min_left_minus=17	max_min_left_plus=15	max_min_right=16	max_min_right_minus=16	max_min_right_plus=12	max_pos_hash_score=68	max_right=34	max_right_minus=34	max_right_plus=31	neg_log10_pos_hash_p_value=0.1	new_junction_coverage=1.31	new_junction_read_count=46	polymorphism_frequency=1.000e+00	pos_hash_score=29	prediction=consensus	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=0.00	side_1_overlap=0	side_1_possible_overlap_registers=34	side_1_read_count=0	side_1_redundant=0	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=0.00	side_2_overlap=0	side_2_possible_overlap_registers=34	side_2_read_count=0	side_2_redundant=0	total_non_overlap_reads=42
JC	64	.	NC_001416-1	9700	-1	NC_001416-2	1	1	0	alignment_overlap=0	coverage_minus=12	coverage_plus=16	flanking_left=35	flanking_right=35	frequency=1	junction_possible_overlap_registers=34	key=NC_001416-1__9700__-1__NC_001416-2__1__1__0____35__35__0__0	max_left=34	max_left_minus=33	max_left_plus=34	max_min_left=16	max_min_left_minus=16	max_min_left_plus=16	max_min_right=17	max_min_right_minus=17	max_min_right_plus=14	max_pos_hash_score=68	max_right=34	max_right_minus=34	max_right_plus=33	neg_log10_pos_hash_p_value=0.0	new_junction_coverage=0.79	new_junction_read_count=29	polymorphism_frequency=1.000e+00	pos_hash_score=21	prediction=consensus	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=0.00	side_1_overlap=0	side_1_possible_overlap_registers=34	side_1_read_count=0	side_1_redundant=0	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=0.00	side_2_overlap=0	side_2_possible_overlap_registers=34	side_2_read_count=0	side_2_redundant=0	total_non_overlap_reads=28
JC	65	.	NC_001416-2	2337	-1	NC_001416-2	8334	1	0	alignment_overlap=5	coverage_minus=10	coverage_plus=11	flanking_left=35	flanking_right=35	frequency=1	junction_possible_overlap_registers=29	key=NC_001416-2__2337__-1__NC_001416-2__8329__1__5____35__35__0__0	max_left=27	max_left_minus=26	max_left_plus=27	max_min_left=9	max_min_left_minus=9	max_min_left_plus=6	max_min_right=15	max_min_right_minus=15	max_min_right_plus=14	max_pos_hash_score=58	max_right=29	max_right_minus=29	max_right_plus=27	neg_log10_pos_hash_p_value=0.0	new_junction_coverage=0.76	new_junction_read_count=25	polymorphism_frequency=9.804e-01	pos_hash_score=16	prediction=consensus	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=0.00	side_1_overlap=5	side_1_possible_overlap_registers=34	side_1_read_count=0	side_1_redundant=0	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=0.03	side_2_overlap=0	side_2_possible_overlap_registers=29	side_2_read_count=1	side_2_redundant=0	total_non_overlap_reads=21
JC	66	.	NC_001416-3	9700	-1	NC_001416-4	1	1	0	alignment_overlap=0	coverage_minus=25	coverage_plus=21	flanking_left=35	flanking_right=35	frequency=1	junction_possible_overlap_registers=34	key=NC_001416-3__9700__-1__NC_001416-4__1__1__0____35__35__0__0	max_left=34	max_left_minus=34	max_left_plus=31	max_min_left=16	max_min_left_minus=16	max_min_left_plus=16	max_min_right=17	max_min_right_minus=17	max_min_right_plus=17	max_pos_hash_score=68	max_right=33	max_right_minus=33	max_right_plus=33	neg_log10_pos_hash_p_value=0.1	new_junction_coverage=1.15	new_junction_read_count=54	polymorphism_frequency=1.000e+00	pos_hash_score=33	prediction=consensus	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=0.00	side_1_overlap=0	side_1_possible_overlap_registers=34	side_1_read_count=0	side_1_redundant=0	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=0.00	side_2_overlap=0	side_2_possible_overlap_registers=34	side_2_read_count=0	side_2_redundant=0	total_non_overlap_reads=46
UN	67	.	NC_001416-0	1	8
UN	68	.	NC_001416-0	34	37
UN	69	.	NC_001416-0	8045	8053
UN	70	.	NC_001416-1	4676	4676
UN	71	.	NC_001416-1	5406	5406
UN	72	.	NC_001416-1	7291	7292
UN	73	.	NC_001416-1	7295	7296
UN	74	.	NC_001416-2	2338	8333
UN	75	.	NC_001416-2	9697	9701
UN	76	.	NC_001416-3	1	10
UN	77	.	NC_001416-4	9666	9666
UN	78	.	NC_001416-4	9668	9701
