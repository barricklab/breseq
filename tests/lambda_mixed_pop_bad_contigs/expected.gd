#=GENOME_DIFF	1.0
#=AUTHOR	
#=REFSEQ	tests/lambda_mixed_pop_bad_contigs/../data/lambda/lambda.gbk
#=REFSEQ	tests/lambda_mixed_pop_bad_contigs/../data/lambda/other.gbk
#=REFSEQ	tests/lambda_mixed_pop_bad_contigs/../data/REL606/REL606.fragment.gbk
#=READSEQ	tests/lambda_mixed_pop_bad_contigs/../data/lambda/lambda_mixed_population.fastq
DEL	60	57	AF322221	1	687
DEL	63	3	NC_001416	139	1
INS	64	6	NC_001416	14266	G
SNP	65	7	NC_001416	20661	G
INS	66	8	NC_001416	20835	C
SNP	67	9	NC_001416	21714	A
DEL	61	17,1	NC_001416	21738	5996
SNP	68	19	NC_001416	31016	C
SNP	69	20	NC_001416	34934	G
DEL	70	21	NC_001416	37818	1
SNP	71	22	NC_001416	45618	C
INS	72	33	NC_001416	46957	A
SNP	73	34	NC_001416	46985	T
SNP	74	35	NC_001416	46992	T
SNP	75	36	NC_001416	47004	A
SNP	76	37	NC_001416	47129	G
SNP	77	38	NC_001416	47143	T
SNP	78	39	NC_001416	47243	A
SNP	79	40	NC_001416	47315	A
SNP	80	41	NC_001416	47317	T
SNP	81	42	NC_001416	47360	A
SNP	82	43	NC_001416	47398	T
SNP	83	44	NC_001416	47509	C
SNP	84	45	NC_001416	47529	T
SNP	85	46	NC_001416	47575	A
SNP	86	47	NC_001416	47669	C
SNP	87	48	NC_001416	47878	G
SNP	88	49	NC_001416	47973	C
SUB	89	50,51	NC_001416	47977	2	AC
DEL	62	58	REL606-5	1	46298
MC	57		AF322221	1	687	0	0	left_inside_cov=0	left_outside_cov=NA	right_inside_cov=0	right_outside_cov=NA
MC	2		NC_001416	1	2	0	0	left_inside_cov=0	left_outside_cov=NA	right_inside_cov=0	right_outside_cov=157
RA	3		NC_001416	139	0	G	.	frequency=1	new_cov=34/40	quality=299.3	ref_cov=0/0	tot_cov=34/40
RA	4		NC_001416	176	0	T	G	frequency=0.0800	new_cov=6/0	polymorphism_quality=11.3002	quality=216.3	ref_cov=48/21	reject=POLYMORPHISM_FREQUENCY_CUTOFF	tot_cov=54/21
RA	5		NC_001416	10551	0	C	A	frequency=0.1485	new_cov=14/1	polymorphism_quality=10.5445	quality=253.7	ref_cov=21/65	tot_cov=35/67
RA	6		NC_001416	14266	1	.	G	frequency=1	new_cov=41/31	quality=203.6	ref_cov=0/0	tot_cov=41/31
RA	7		NC_001416	20661	0	N	G	frequency=1	new_cov=22/58	quality=222.5	ref_cov=0/0	tot_cov=22/60
RA	8		NC_001416	20835	1	.	C	frequency=1	new_cov=32/63	quality=279.1	ref_cov=0/0	tot_cov=32/63
RA	9		NC_001416	21714	0	G	A	frequency=1	new_cov=32/61	quality=286.3	ref_cov=0/1	tot_cov=33/62
JC	1		NC_001416	21737	-1	NC_001416	27734	1	0	alignment_overlap=5	coverage_minus=35	coverage_plus=41	flanking_left=35	flanking_right=35	key=NC_001416__21737__-1__NC_001416__27729__1__5____35__35__0__0	max_left=29	max_left_minus=27	max_left_plus=29	max_min_left=13	max_min_left_minus=13	max_min_left_plus=11	max_min_right=15	max_min_right_minus=15	max_min_right_plus=15	max_pos_hash_score=58	max_right=29	max_right_minus=29	max_right_plus=28	neg_log10_pos_hash_p_value=0.1	pos_hash_score=40	side_1_annotate_key=gene	side_1_overlap=5	side_1_redundant=0	side_2_annotate_key=gene	side_2_overlap=0	side_2_redundant=0	total_non_overlap_reads=76
MC	17		NC_001416	21738	27733	0	0	left_inside_cov=0	left_outside_cov=102	right_inside_cov=0	right_outside_cov=102
RA	10		NC_001416	22444	0	A	G	deleted=1	frequency=1	new_cov=1/0	quality=3.1	ref_cov=0/0	reject=EVALUE	tot_cov=1/0
RA	11		NC_001416	22470	0	C	A	deleted=1	frequency=1	new_cov=1/0	quality=2.4	ref_cov=0/0	reject=EVALUE	tot_cov=1/0
RA	12		NC_001416	27555	0	A	G	deleted=1	frequency=1	new_cov=1/0	quality=3.1	ref_cov=0/0	reject=EVALUE	tot_cov=1/0
RA	13		NC_001416	27563	1	.	A	deleted=1	frequency=1	new_cov=1/0	quality=3.0	ref_cov=0/0	reject=EVALUE	tot_cov=1/0
RA	14		NC_001416	27601	0	C	T	deleted=1	frequency=1	new_cov=1/0	quality=3.0	ref_cov=0/0	reject=EVALUE	tot_cov=1/0
RA	15		NC_001416	27653	0	A	T	deleted=1	frequency=1	new_cov=0/1	quality=3.0	ref_cov=0/0	reject=EVALUE	tot_cov=0/1
RA	16		NC_001416	27668	0	T	G	deleted=1	frequency=1	new_cov=0/1	quality=2.5	ref_cov=0/0	reject=EVALUE	tot_cov=0/1
RA	19		NC_001416	31016	0	T	C	frequency=1	new_cov=57/47	quality=294.9	ref_cov=0/0	tot_cov=58/47
RA	20		NC_001416	34934	0	A	G	frequency=1	new_cov=21/44	quality=182.8	ref_cov=0/0	tot_cov=21/45
RA	21		NC_001416	37818	0	C	.	frequency=1	new_cov=29/17	quality=180.5	ref_cov=1/0	tot_cov=30/17
RA	22		NC_001416	45618	0	T	C	frequency=1	new_cov=58/78	quality=393.1	ref_cov=0/1	tot_cov=60/79
RA	23		NC_001416	46154	0	G	A	frequency=0.0759	new_cov=8/3	polymorphism_quality=13.6230	quality=358.7	ref_cov=87/47	reject=POLYMORPHISM_FREQUENCY_CUTOFF	tot_cov=96/50
RA	24		NC_001416	46157	0	C	G	frequency=0.0633	new_cov=7/3	polymorphism_quality=17.0458	quality=408.3	ref_cov=102/46	reject=POLYMORPHISM_FREQUENCY_CUTOFF	tot_cov=109/49
RA	25		NC_001416	46162	0	C	A	frequency=0.0500	new_cov=5/3	polymorphism_quality=11.2923	quality=429.0	ref_cov=97/55	reject=POLYMORPHISM_FREQUENCY_CUTOFF	tot_cov=102/58
RA	26		NC_001416	46185	0	G	A	frequency=0.0761	new_cov=6/8	polymorphism_quality=18.9822	quality=458.1	ref_cov=97/73	reject=POLYMORPHISM_FREQUENCY_CUTOFF	tot_cov=103/82
RA	27		NC_001416	46190	0	T	G	frequency=0.1084	new_cov=8/10	polymorphism_quality=27.8527	quality=462.7	ref_cov=81/67	tot_cov=90/77
RA	28		NC_001416	46430	0	T	C	frequency=0.2353	new_cov=16/4	polymorphism_quality=30.8289	quality=148.1	ref_cov=28/37	tot_cov=44/41
RA	29		NC_001416	46597	0	G	A	frequency=0.2000	new_cov=6/9	polymorphism_quality=26.3887	quality=132.9	ref_cov=39/21	tot_cov=45/31
RA	30		NC_001416	46633	0	A	G	frequency=0.1458	new_cov=9/5	polymorphism_quality=20.7638	quality=213.6	ref_cov=43/39	tot_cov=52/44
RA	31		NC_001416	46659	0	C	T	frequency=0.1102	new_cov=11/2	polymorphism_quality=19.0931	quality=257.4	ref_cov=66/39	tot_cov=78/42
RA	32		NC_001416	46680	0	G	C	frequency=0.0902	new_cov=8/3	polymorphism_quality=13.3108	quality=317.9	ref_cov=53/58	reject=POLYMORPHISM_FREQUENCY_CUTOFF	tot_cov=61/61
RA	33		NC_001416	46957	1	.	A	frequency=1	new_cov=25/29	quality=170.4	ref_cov=2/0	tot_cov=27/29
RA	34		NC_001416	46985	0	C	T	frequency=1	new_cov=35/20	quality=171.1	ref_cov=0/0	tot_cov=35/20
RA	35		NC_001416	46992	0	C	T	frequency=1	new_cov=32/17	quality=142.1	ref_cov=2/0	tot_cov=34/17
RA	36		NC_001416	47004	0	G	A	frequency=1	new_cov=12/21	quality=92.9	ref_cov=0/0	tot_cov=12/21
RA	37		NC_001416	47129	0	A	G	frequency=1	new_cov=32/26	quality=171.5	ref_cov=0/0	tot_cov=32/26
RA	38		NC_001416	47143	0	C	T	frequency=1	new_cov=29/23	quality=161.4	ref_cov=1/0	tot_cov=30/24
RA	39		NC_001416	47243	0	G	A	frequency=1	new_cov=25/41	quality=204.0	ref_cov=0/1	tot_cov=25/42
RA	40		NC_001416	47315	0	G	A	frequency=1	new_cov=42/51	quality=290.9	ref_cov=1/0	tot_cov=45/51
RA	41		NC_001416	47317	0	N	T	frequency=1	new_cov=50/51	quality=316.0	ref_cov=0/0	tot_cov=50/52
RA	42		NC_001416	47360	0	G	A	frequency=1	new_cov=24/45	quality=216.4	ref_cov=0/0	tot_cov=24/45
RA	43		NC_001416	47398	0	C	T	frequency=1	new_cov=33/44	quality=238.1	ref_cov=1/0	tot_cov=34/44
RA	44		NC_001416	47509	0	T	C	frequency=1	new_cov=34/25	quality=168.7	ref_cov=0/0	tot_cov=34/25
RA	45		NC_001416	47529	0	C	T	frequency=1	new_cov=41/41	quality=258.1	ref_cov=0/0	tot_cov=41/42
RA	46		NC_001416	47575	0	C	A	frequency=1	new_cov=32/45	quality=241.4	ref_cov=1/0	tot_cov=33/45
RA	47		NC_001416	47669	0	T	C	frequency=1	new_cov=45/44	quality=256.3	ref_cov=0/0	tot_cov=45/44
RA	48		NC_001416	47878	0	A	G	frequency=1	new_cov=35/51	quality=241.0	ref_cov=0/0	tot_cov=35/51
RA	49		NC_001416	47973	0	T	C	frequency=1	new_cov=44/48	quality=269.4	ref_cov=0/0	tot_cov=44/48
RA	50		NC_001416	47977	0	G	A	frequency=1	new_cov=36/39	quality=224.5	ref_cov=0/0	tot_cov=36/39
RA	51		NC_001416	47978	0	T	C	frequency=1	new_cov=36/38	quality=214.1	ref_cov=0/0	tot_cov=36/39
RA	52		NC_001416	48160	0	T	C	frequency=0.8161	new_cov=35/36	polymorphism_quality=28.9175	quality=151.8	ref_cov=6/10	tot_cov=41/46
RA	53		NC_001416	48202	0	G	A	frequency=0.1310	new_cov=2/9	polymorphism_quality=16.5709	quality=178.2	ref_cov=43/30	tot_cov=45/39
RA	54		NC_001416	48295	0	C	A	frequency=0.1667	new_cov=10/10	polymorphism_quality=36.9422	quality=271.4	ref_cov=53/47	tot_cov=63/57
MC	55		NC_001416	48467	48502	0	0	left_inside_cov=28	left_outside_cov=29	right_inside_cov=0	right_outside_cov=NA
MC	58		REL606-5	1	46298	0	0	left_inside_cov=0	left_outside_cov=NA	right_inside_cov=0	right_outside_cov=NA
UN	18		NC_001416	1	8
UN	56		NC_001416	21738	27733
UN	59		NC_001416	48494	48502
