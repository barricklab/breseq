#=GENOME_DIFF	1.0
#=COMMAND	./src/c/breseq/breseq -j 4 -o tests/tmv_plasmid_circular_deletion_start_only -r tests/tmv_plasmid_circular_deletion_start_only/../data/tmv_plasmid/tmv-plasmid-truncate-end.gbk tests/tmv_plasmid_circular_deletion_start_only/../data/tmv_plasmid/D3-9_1P.fastq.gz tests/tmv_plasmid_circular_deletion_start_only/../data/tmv_plasmid/D3-9_2P.fastq.gz
#=REFSEQ	tests/tmv_plasmid_circular_deletion_start_only/../data/tmv_plasmid/tmv-plasmid-truncate-end.gbk
#=READSEQ	tests/tmv_plasmid_circular_deletion_start_only/../data/tmv_plasmid/D3-9_1P.fastq.gz
#=READSEQ	tests/tmv_plasmid_circular_deletion_start_only/../data/tmv_plasmid/D3-9_2P.fastq.gz
#=CONVERTED-BASES	2735772
#=CONVERTED-READS	18371
#=INPUT-BASES	2738037
#=INPUT-READS	18386
#=MAPPED-BASES	2284581
#=MAPPED-READS	15427
DEL	1	4,6	TMV-plasmid-truncate-end	1	7137	gene_name=–/–	gene_position=intergenic (–/–)	gene_product=–/–	gene_strand=–/–	locus_tag=–/–	mutation_category=large_deletion	position_end=7137	position_start=1	ref_seq=7137-bp
RA	2	.	TMV-plasmid-truncate-end	1565	0	C	T	consensus_score=10.1	deleted=1	frequency=1	major_base=T	major_cov=1/3	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=1/3	polymorphism_frequency=1.000e+00	polymorphism_score=NA	prediction=consensus	ref_cov=0/0	total_cov=1/3
RA	3	.	TMV-plasmid-truncate-end	5198	0	A	G	consensus_score=18.2	deleted=1	frequency=1	major_base=G	major_cov=3/3	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=3/3	polymorphism_frequency=1.000e+00	polymorphism_score=NA	prediction=consensus	ref_cov=0/0	total_cov=3/3
MC	4	.	TMV-plasmid-truncate-end	1	7137	0	0	gene_name=–/–	gene_position=intergenic (–/–)	gene_product=–/–	gene_strand=–/–	left_inside_cov=0	left_outside_cov=NA	locus_tag=–/–	right_inside_cov=41	right_outside_cov=852
JC	5	.	TMV-plasmid-truncate-end	7033	1	TMV-plasmid-truncate-end	10204	-1	-49	alignment_overlap=-49	coverage_minus=19	coverage_plus=12	flanking_left=151	flanking_right=151	frequency=1	junction_possible_overlap_registers=98	key=TMV-plasmid-truncate-end__7033__1__TMV-plasmid-truncate-end__10204__-1__-49__CAACACCCGCTGACGCGCCCTGACGGGCTTGTCTGCTCCCGGCATCCGC__151__151__0__0	max_left=87	max_left_minus=86	max_left_plus=87	max_min_left=50	max_min_left_minus=33	max_min_left_plus=50	max_min_right=39	max_min_right_minus=34	max_min_right_plus=39	max_pos_hash_score=198	max_right=101	max_right_minus=101	max_right_plus=89	neg_log10_pos_hash_p_value=1.4	new_junction_coverage=0.17	new_junction_read_count=31	polymorphism_frequency=9.029e-01	pos_hash_score=20	prediction=consensus	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=0.04	side_1_gene_name=–/–	side_1_gene_position=intergenic (–/–)	side_1_gene_product=–/–	side_1_gene_strand=–/–	side_1_locus_tag=–/–	side_1_overlap=0	side_1_possible_overlap_registers=147	side_1_read_count=10	side_1_redundant=0	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=0.00	side_2_gene_name=–/–	side_2_gene_position=intergenic (–/–)	side_2_gene_product=–/–	side_2_gene_strand=–/–	side_2_locus_tag=–/–	side_2_overlap=0	side_2_possible_overlap_registers=147	side_2_read_count=0	side_2_redundant=0	total_non_overlap_reads=31	unique_read_sequence=CAACACCCGCTGACGCGCCCTGACGGGCTTGTCTGCTCCCGGCATCCGC
JC	6	.	TMV-plasmid-truncate-end	7138	1	TMV-plasmid-truncate-end	10204	-1	0	alignment_overlap=0	coverage_minus=405	coverage_plus=394	flanking_left=151	flanking_right=151	frequency=1	junction_possible_overlap_registers=147	key=TMV-plasmid-truncate-end__7138__1__TMV-plasmid-truncate-end__10204__-1__0____151__151__0__0	max_left=150	max_left_minus=150	max_left_plus=150	max_min_left=75	max_min_left_minus=75	max_min_left_plus=74	max_min_right=75	max_min_right_minus=75	max_min_right_plus=75	max_pos_hash_score=296	max_right=149	max_right_minus=147	max_right_plus=149	neg_log10_pos_hash_p_value=-0.0	new_junction_coverage=2.92	new_junction_read_count=809	polymorphism_frequency=9.800e-01	pos_hash_score=214	prediction=consensus	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=0.12	side_1_gene_name=–/–	side_1_gene_position=intergenic (–/–)	side_1_gene_product=–/–	side_1_gene_strand=–/–	side_1_locus_tag=–/–	side_1_overlap=0	side_1_possible_overlap_registers=147	side_1_read_count=33	side_1_redundant=0	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=0.00	side_2_gene_name=–/–	side_2_gene_position=intergenic (–/–)	side_2_gene_product=–/–	side_2_gene_strand=–/–	side_2_locus_tag=–/–	side_2_overlap=0	side_2_possible_overlap_registers=147	side_2_read_count=0	side_2_redundant=0	total_non_overlap_reads=799
JC	7	.	TMV-plasmid-truncate-end	7496	1	TMV-plasmid-truncate-end	7789	-1	0	alignment_overlap=0	coverage_minus=7	coverage_plus=7	flanking_left=151	flanking_right=151	frequency=1.574e-02	junction_possible_overlap_registers=147	key=TMV-plasmid-truncate-end__7496__1__TMV-plasmid-truncate-end__7789__-1__0____151__151__0__0	max_left=144	max_left_minus=144	max_left_plus=118	max_min_left=6	max_min_left_minus=3	max_min_left_plus=6	max_min_right=70	max_min_right_minus=70	max_min_right_plus=67	max_pos_hash_score=296	max_right=150	max_right_minus=149	max_right_plus=150	neg_log10_pos_hash_p_value=5.6	new_junction_coverage=0.05	new_junction_read_count=14	polymorphism_frequency=1.574e-02	pos_hash_score=9	prediction=polymorphism	reject=COVERAGE_EVENNESS_SKEW,FREQUENCY_CUTOFF	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=3.27	side_1_gene_name=–/–	side_1_gene_position=intergenic (–/–)	side_1_gene_product=–/–	side_1_gene_strand=–/–	side_1_locus_tag=–/–	side_1_overlap=0	side_1_possible_overlap_registers=147	side_1_read_count=906	side_1_redundant=0	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=3.05	side_2_gene_name=–/–	side_2_gene_position=intergenic (–/–)	side_2_gene_product=–/–	side_2_gene_strand=–/–	side_2_locus_tag=–/–	side_2_overlap=0	side_2_possible_overlap_registers=147	side_2_read_count=845	side_2_redundant=0	total_non_overlap_reads=14
JC	8	.	TMV-plasmid-truncate-end	9443	1	TMV-plasmid-truncate-end	9519	1	0	alignment_overlap=8	coverage_minus=1	coverage_plus=4	flanking_left=151	flanking_right=151	frequency=6.755e-03	junction_possible_overlap_registers=139	key=TMV-plasmid-truncate-end__9443__1__TMV-plasmid-truncate-end__9511__1__8____151__151__0__0	max_left=83	max_left_minus=71	max_left_plus=83	max_min_left=71	max_min_left_minus=71	max_min_left_plus=61	max_min_right=56	max_min_right_minus=0	max_min_right_plus=56	max_pos_hash_score=280	max_right=140	max_right_minus=72	max_right_plus=140	neg_log10_pos_hash_p_value=7.0	new_junction_coverage=0.02	new_junction_read_count=6	polymorphism_frequency=6.755e-03	pos_hash_score=5	prediction=polymorphism	reject=COVERAGE_EVENNESS_SKEW,FREQUENCY_CUTOFF	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=3.29	side_1_gene_name=–/–	side_1_gene_position=intergenic (–/–)	side_1_gene_product=–/–	side_1_gene_strand=–/–	side_1_locus_tag=–/–	side_1_overlap=8	side_1_possible_overlap_registers=147	side_1_read_count=912	side_1_redundant=0	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=3.44	side_2_gene_name=–/–	side_2_gene_position=intergenic (–/–)	side_2_gene_product=–/–	side_2_gene_strand=–/–	side_2_locus_tag=–/–	side_2_overlap=0	side_2_possible_overlap_registers=139	side_2_read_count=902	side_2_redundant=0	total_non_overlap_reads=5
UN	9	.	TMV-plasmid-truncate-end	1	185
UN	10	.	TMV-plasmid-truncate-end	191	191
UN	11	.	TMV-plasmid-truncate-end	204	204
UN	12	.	TMV-plasmid-truncate-end	210	211
UN	13	.	TMV-plasmid-truncate-end	213	213
UN	14	.	TMV-plasmid-truncate-end	217	217
UN	15	.	TMV-plasmid-truncate-end	219	219
UN	16	.	TMV-plasmid-truncate-end	231	231
UN	17	.	TMV-plasmid-truncate-end	238	238
UN	18	.	TMV-plasmid-truncate-end	241	241
UN	19	.	TMV-plasmid-truncate-end	268	268
UN	20	.	TMV-plasmid-truncate-end	274	274
UN	21	.	TMV-plasmid-truncate-end	283	283
UN	22	.	TMV-plasmid-truncate-end	294	294
UN	23	.	TMV-plasmid-truncate-end	298	298
UN	24	.	TMV-plasmid-truncate-end	301	301
UN	25	.	TMV-plasmid-truncate-end	304	304
UN	26	.	TMV-plasmid-truncate-end	306	772
UN	27	.	TMV-plasmid-truncate-end	779	1544
UN	28	.	TMV-plasmid-truncate-end	1564	1564
UN	29	.	TMV-plasmid-truncate-end	1566	1567
UN	30	.	TMV-plasmid-truncate-end	1569	1571
UN	31	.	TMV-plasmid-truncate-end	1576	1576
UN	32	.	TMV-plasmid-truncate-end	1579	1579
UN	33	.	TMV-plasmid-truncate-end	1582	1583
UN	34	.	TMV-plasmid-truncate-end	1586	1586
UN	35	.	TMV-plasmid-truncate-end	1590	1592
UN	36	.	TMV-plasmid-truncate-end	1598	1598
UN	37	.	TMV-plasmid-truncate-end	1600	1754
UN	38	.	TMV-plasmid-truncate-end	1764	1764
UN	39	.	TMV-plasmid-truncate-end	1774	1774
UN	40	.	TMV-plasmid-truncate-end	1777	1777
UN	41	.	TMV-plasmid-truncate-end	1779	1779
UN	42	.	TMV-plasmid-truncate-end	1781	1781
UN	43	.	TMV-plasmid-truncate-end	1785	1785
UN	44	.	TMV-plasmid-truncate-end	1958	1958
UN	45	.	TMV-plasmid-truncate-end	1964	1964
UN	46	.	TMV-plasmid-truncate-end	1969	2427
UN	47	.	TMV-plasmid-truncate-end	2430	2430
UN	48	.	TMV-plasmid-truncate-end	2435	2435
UN	49	.	TMV-plasmid-truncate-end	2441	2441
UN	50	.	TMV-plasmid-truncate-end	2540	2540
UN	51	.	TMV-plasmid-truncate-end	2547	2548
UN	52	.	TMV-plasmid-truncate-end	2552	2552
UN	53	.	TMV-plasmid-truncate-end	2563	2563
UN	54	.	TMV-plasmid-truncate-end	2570	2879
UN	55	.	TMV-plasmid-truncate-end	2882	2882
UN	56	.	TMV-plasmid-truncate-end	2888	2888
UN	57	.	TMV-plasmid-truncate-end	2935	2936
UN	58	.	TMV-plasmid-truncate-end	2938	2938
UN	59	.	TMV-plasmid-truncate-end	2944	2944
UN	60	.	TMV-plasmid-truncate-end	2952	2956
UN	61	.	TMV-plasmid-truncate-end	2989	2989
UN	62	.	TMV-plasmid-truncate-end	3024	3031
UN	63	.	TMV-plasmid-truncate-end	3037	3037
UN	64	.	TMV-plasmid-truncate-end	3044	3044
UN	65	.	TMV-plasmid-truncate-end	3056	3056
UN	66	.	TMV-plasmid-truncate-end	3062	3062
UN	67	.	TMV-plasmid-truncate-end	3070	3070
UN	68	.	TMV-plasmid-truncate-end	3072	3072
UN	69	.	TMV-plasmid-truncate-end	3076	3076
UN	70	.	TMV-plasmid-truncate-end	3109	3109
UN	71	.	TMV-plasmid-truncate-end	3181	3181
UN	72	.	TMV-plasmid-truncate-end	3183	3189
UN	73	.	TMV-plasmid-truncate-end	3194	3195
UN	74	.	TMV-plasmid-truncate-end	3199	3199
UN	75	.	TMV-plasmid-truncate-end	3202	3202
UN	76	.	TMV-plasmid-truncate-end	3204	3204
UN	77	.	TMV-plasmid-truncate-end	3206	3206
UN	78	.	TMV-plasmid-truncate-end	3213	3214
UN	79	.	TMV-plasmid-truncate-end	3216	3216
UN	80	.	TMV-plasmid-truncate-end	3218	3219
UN	81	.	TMV-plasmid-truncate-end	3222	3222
UN	82	.	TMV-plasmid-truncate-end	3224	3225
UN	83	.	TMV-plasmid-truncate-end	3228	3229
UN	84	.	TMV-plasmid-truncate-end	3235	3235
UN	85	.	TMV-plasmid-truncate-end	3239	3292
UN	86	.	TMV-plasmid-truncate-end	3298	3322
UN	87	.	TMV-plasmid-truncate-end	3408	3408
UN	88	.	TMV-plasmid-truncate-end	3416	3416
UN	89	.	TMV-plasmid-truncate-end	3438	3576
UN	90	.	TMV-plasmid-truncate-end	3582	3582
UN	91	.	TMV-plasmid-truncate-end	3739	3739
UN	92	.	TMV-plasmid-truncate-end	3758	3758
UN	93	.	TMV-plasmid-truncate-end	3763	3802
UN	94	.	TMV-plasmid-truncate-end	3837	3837
UN	95	.	TMV-plasmid-truncate-end	3842	3842
UN	96	.	TMV-plasmid-truncate-end	3846	3846
UN	97	.	TMV-plasmid-truncate-end	3864	3864
UN	98	.	TMV-plasmid-truncate-end	3866	3877
UN	99	.	TMV-plasmid-truncate-end	3879	3883
UN	100	.	TMV-plasmid-truncate-end	3887	3887
UN	101	.	TMV-plasmid-truncate-end	3891	3892
UN	102	.	TMV-plasmid-truncate-end	3894	3894
UN	103	.	TMV-plasmid-truncate-end	3896	3896
UN	104	.	TMV-plasmid-truncate-end	3900	3900
UN	105	.	TMV-plasmid-truncate-end	3903	3904
UN	106	.	TMV-plasmid-truncate-end	3906	3907
UN	107	.	TMV-plasmid-truncate-end	3909	3913
UN	108	.	TMV-plasmid-truncate-end	3915	3915
UN	109	.	TMV-plasmid-truncate-end	3917	3921
UN	110	.	TMV-plasmid-truncate-end	3926	3927
UN	111	.	TMV-plasmid-truncate-end	3929	3930
UN	112	.	TMV-plasmid-truncate-end	3932	3932
UN	113	.	TMV-plasmid-truncate-end	3934	3934
UN	114	.	TMV-plasmid-truncate-end	3937	3938
UN	115	.	TMV-plasmid-truncate-end	3940	4236
UN	116	.	TMV-plasmid-truncate-end	4239	4241
UN	117	.	TMV-plasmid-truncate-end	4393	4393
UN	118	.	TMV-plasmid-truncate-end	4397	4398
UN	119	.	TMV-plasmid-truncate-end	4401	4402
UN	120	.	TMV-plasmid-truncate-end	4404	4407
UN	121	.	TMV-plasmid-truncate-end	4416	4416
UN	122	.	TMV-plasmid-truncate-end	4418	4432
UN	123	.	TMV-plasmid-truncate-end	4434	4434
UN	124	.	TMV-plasmid-truncate-end	4436	4682
UN	125	.	TMV-plasmid-truncate-end	4686	4686
UN	126	.	TMV-plasmid-truncate-end	4831	4831
UN	127	.	TMV-plasmid-truncate-end	4846	4846
UN	128	.	TMV-plasmid-truncate-end	4850	5173
UN	129	.	TMV-plasmid-truncate-end	5299	5299
UN	130	.	TMV-plasmid-truncate-end	5303	5400
UN	131	.	TMV-plasmid-truncate-end	5406	5406
UN	132	.	TMV-plasmid-truncate-end	5483	5483
UN	133	.	TMV-plasmid-truncate-end	5493	5493
UN	134	.	TMV-plasmid-truncate-end	5509	5542
UN	135	.	TMV-plasmid-truncate-end	5545	5551
UN	136	.	TMV-plasmid-truncate-end	5557	5558
UN	137	.	TMV-plasmid-truncate-end	5562	5562
UN	138	.	TMV-plasmid-truncate-end	5564	5564
UN	139	.	TMV-plasmid-truncate-end	5586	5586
UN	140	.	TMV-plasmid-truncate-end	5595	5595
UN	141	.	TMV-plasmid-truncate-end	5598	5598
UN	142	.	TMV-plasmid-truncate-end	5625	5625
UN	143	.	TMV-plasmid-truncate-end	5634	5634
UN	144	.	TMV-plasmid-truncate-end	5643	5643
UN	145	.	TMV-plasmid-truncate-end	5673	5673
UN	146	.	TMV-plasmid-truncate-end	5677	5677
UN	147	.	TMV-plasmid-truncate-end	5688	6662
UN	148	.	TMV-plasmid-truncate-end	6800	6800
UN	149	.	TMV-plasmid-truncate-end	6802	6802
UN	150	.	TMV-plasmid-truncate-end	6804	6888
