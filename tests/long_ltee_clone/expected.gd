#=GENOME_DIFF	1.0
#=CREATED	19:35:46 12 Aug 2026
#=PROGRAM	breseq 0.50.0 revision c333ad77bdea
#=COMMAND	./src/breseq/breseq -j 7 -o tests/long_ltee_clone -r tests/long_ltee_clone/../data/downloads/ltee_REL606/REL606.gbk --predict-copy-number tests/long_ltee_clone/../data/downloads/ena_SRR2589061/SRR2589061_1.fastq.gz tests/long_ltee_clone/../data/downloads/ena_SRR2589061/SRR2589061_2.fastq.gz
#=REFSEQ	tests/long_ltee_clone/../data/downloads/ltee_REL606/REL606.gbk
#=READSEQ	tests/long_ltee_clone/../data/downloads/ena_SRR2589061/SRR2589061_1.fastq.gz
#=READSEQ	tests/long_ltee_clone/../data/downloads/ena_SRR2589061/SRR2589061_2.fastq.gz
#=CONVERTED-BASES	481720712
#=CONVERTED-READS	4769512
#=INPUT-BASES	481722328
#=INPUT-READS	4769528
#=MAPPED-BASES	473713701
#=MAPPED-READS	4690619
MOB	1	47,48	REL606	467507	IS150	-1	3	gene_name=htpG	gene_position=coding (321-323/1875 nt)	gene_product=heat shock protein 90	gene_strand=>	genes_inactivated=htpG	locus_tag=ECB_00424	locus_tags_inactivated=ECB_00424	mutation_category=mobile_element_insertion	position_end=467509	position_start=467507	ref_seq=ATC	repeat_size=1443
DEL	2	43,49	REL606	474383	9	gene_name=ybaL	gene_position=coding (915-923/1677 nt)	gene_product=predicted transporter with NAD(P)-binding Rossmann-fold domain	gene_strand=<	genes_overlapping=ybaL	locus_tag=ECB_00429	locus_tags_overlapping=ECB_00429	mutation_category=small_indel	position_end=474391	position_start=474383	ref_seq=GTCGCCAGC
SNP	3	28	REL606	651601	A	aa_new_seq=L	aa_position=7	aa_ref_seq=Q	codon_new_seq=CTG	codon_number=7	codon_position=2	codon_ref_seq=CAG	gene_name=ybeB	gene_position=20	gene_product=hypothetical protein	gene_strand=<	genes_overlapping=ybeB	locus_tag=ECB_00606	locus_tags_overlapping=ECB_00606	mutation_category=snp_nonsynonymous	position_end=651601	position_start=651601	ref_seq=T	snp_type=nonsynonymous	transl_table=11
MOB	4	56,62	REL606	968326	IS150	-1	3	gene_name=pflA/pflB	gene_position=intergenic (-153/+37)	gene_product=pyruvate formate lyase activating enzyme 1/pyruvate formate lyase I	gene_strand=</<	locus_tag=ECB_00906/ECB_00907	mutation_category=mobile_element_insertion	position_end=968328	position_start=968326	ref_seq=ATT	repeat_size=1443
SNP	5	29	REL606	1089367	A	aa_new_seq=V	aa_position=310	aa_ref_seq=A	codon_new_seq=GTC	codon_number=310	codon_position=2	codon_ref_seq=GCC	gene_name=ycdM	gene_position=929	gene_product=predicted monooxygenase	gene_strand=<	genes_overlapping=ycdM	locus_tag=ECB_01015	locus_tags_overlapping=ECB_01015	mutation_category=snp_nonsynonymous	position_end=1089367	position_start=1089367	ref_seq=G	snp_type=nonsynonymous	transl_table=11
SNP	6	30	REL606	1166950	T	aa_new_seq=I	aa_position=148	aa_ref_seq=T	codon_new_seq=ATT	codon_number=148	codon_position=2	codon_ref_seq=ACT	gene_name=fabF	gene_position=443	gene_product=3-oxoacyl-(acyl carrier protein) synthase	gene_strand=>	genes_overlapping=fabF	locus_tag=ECB_01091	locus_tags_overlapping=ECB_01091	mutation_category=snp_nonsynonymous	position_end=1166950	position_start=1166950	ref_seq=C	snp_type=nonsynonymous	transl_table=11
MOB	7	63,64	REL606	1305934	IS1	-1	9	gene_name=cls	gene_position=coding (298-306/1461 nt)	gene_product=cardiolipin synthetase	gene_strand=<	genes_inactivated=cls	locus_tag=ECB_01223	locus_tags_inactivated=ECB_01223	mutation_category=mobile_element_insertion	position_end=1305942	position_start=1305934	ref_seq=TGGCGCAGC	repeat_size=768
SNP	8	31	REL606	1330018	A	aa_new_seq=E	aa_position=200	aa_ref_seq=A	codon_new_seq=GAG	codon_number=200	codon_position=2	codon_ref_seq=GCG	gene_name=topA	gene_position=599	gene_product=DNA topoisomerase I	gene_strand=>	genes_overlapping=topA	locus_tag=ECB_01250	locus_tags_overlapping=ECB_01250	mutation_category=snp_nonsynonymous	position_end=1330018	position_start=1330018	ref_seq=C	snp_type=nonsynonymous	transl_table=11
MOB	9	50,57	REL606	1462251	IS150	-1	3	gene_name=mokB/trg	gene_position=intergenic (-13/-326)	gene_product=regulatory peptide/methyl-accepting chemotaxis protein III, ribose and galactose sensor receptor	gene_strand=</>	genes_promoter=mokB	locus_tag=ECB_01377/ECB_01378	locus_tags_promoter=ECB_01377	mutation_category=mobile_element_insertion	position_end=1462253	position_start=1462251	ref_seq=GTT	repeat_size=1443
SNP	10	32	REL606	1587635	G	aa_new_seq=C	aa_position=99	aa_ref_seq=Y	codon_new_seq=TGC	codon_number=99	codon_position=2	codon_ref_seq=TAC	gene_name=marA	gene_position=296	gene_product=DNA-binding transcriptional dual activator of multiple antibiotic resistance	gene_strand=>	genes_overlapping=marA	locus_tag=ECB_01490	locus_tags_overlapping=ECB_01490	mutation_category=snp_nonsynonymous	position_end=1587635	position_start=1587635	ref_seq=A	snp_type=nonsynonymous	transl_table=11
SNP	11	33	REL606	1733865	T	aa_new_seq=S	aa_position=301	aa_ref_seq=A	codon_new_seq=TCC	codon_number=301	codon_position=1	codon_ref_seq=GCC	gene_name=pykF	gene_position=901	gene_product=pyruvate kinase	gene_strand=>	genes_overlapping=pykF	locus_tag=ECB_01645	locus_tags_overlapping=ECB_01645	mutation_category=snp_nonsynonymous	position_end=1733865	position_start=1733865	ref_seq=G	snp_type=nonsynonymous	transl_table=11
DEL	12	44,65	REL606	2100308	22146	gene_name=ogrK–[ECB_02013]	gene_product=ogrK,yegZ,ECB_01989,ECB_01990,ECB_01991,ECB_01992,ECB_01993,ECB_01994,ECB_01995,ECB_01996,ECB_01997,ECB_01998,ECB_01999,ECB_02000,ECB_02001,ECB_02002,ECB_02003,ECB_02004,ECB_02005,ECB_02006,ECB_02007,ECB_02008,ECB_02009,ECB_02010,ECB_02011,ECB_02012,[ECB_02013]	genes_inactivated=ogrK,yegZ,ECB_01989,ECB_01990,ECB_01991,ECB_01992,ECB_01993,ECB_01994,ECB_01995,ECB_01996,ECB_01997,ECB_01998,ECB_01999,ECB_02000,ECB_02001,ECB_02002,ECB_02003,ECB_02004,ECB_02005,ECB_02006,ECB_02007,ECB_02008,ECB_02009,ECB_02010,ECB_02011,ECB_02012	genes_overlapping=ECB_02013	locus_tag=[ECB_01987]–[ECB_02013]	locus_tags_inactivated=ECB_01987,ECB_01988,ECB_01989,ECB_01990,ECB_01991,ECB_01992,ECB_01993,ECB_01994,ECB_01995,ECB_01996,ECB_01997,ECB_01998,ECB_01999,ECB_02000,ECB_02001,ECB_02002,ECB_02003,ECB_02004,ECB_02005,ECB_02006,ECB_02007,ECB_02008,ECB_02009,ECB_02010,ECB_02011,ECB_02012	locus_tags_overlapping=ECB_02013	mutation_category=large_deletion	position_end=2122453	position_start=2100308	ref_seq=22146-bp
MOB	13	51,58	REL606	2157881	IS150	1	3	gene_name=yehM	gene_position=coding (991-993/2280 nt)	gene_product=hypothetical protein	gene_strand=>	genes_inactivated=yehM	locus_tag=ECB_02049	locus_tags_inactivated=ECB_02049	mutation_category=mobile_element_insertion	position_end=2157883	position_start=2157881	ref_seq=AAC	repeat_size=1443
SNP	14	34	REL606	2275076	G	aa_new_seq=H	aa_position=474	aa_ref_seq=Q	codon_new_seq=CAC	codon_number=474	codon_position=3	codon_ref_seq=CAA	gene_name=yfaQ	gene_position=1422	gene_product=hypothetical protein	gene_strand=<	genes_overlapping=yfaQ	locus_tag=ECB_02153	locus_tags_overlapping=ECB_02153	mutation_category=snp_nonsynonymous	position_end=2275076	position_start=2275076	ref_seq=T	snp_type=nonsynonymous	transl_table=11
MOB	15	66,67	REL606	2283472	IS150	-1	3	gene_name=yfaA/gyrA	gene_position=intergenic (-128/+19)	gene_product=hypothetical protein/DNA gyrase subunit A	gene_strand=</<	genes_promoter=yfaA	locus_tag=ECB_02156/ECB_02157	locus_tags_promoter=ECB_02156	mutation_category=mobile_element_insertion	position_end=2283474	position_start=2283472	ref_seq=TTT	repeat_size=1443
MOB	16	68,69	REL606	2448493	IS186	1	6	gene_name=nupC/yfeA	gene_position=intergenic (+24/+21)	gene_product=nucleoside (except guanosine) transporter/predicted diguanylate cyclase	gene_strand=>/<	locus_tag=ECB_02302/ECB_02303	mutation_category=mobile_element_insertion	position_end=2448498	position_start=2448493	ref_seq=ATAATT	repeat_size=1343
AMP	17	70,72	REL606	2749275	26603	2	gene_name=[fhlA]–insJ-3	gene_product=[fhlA],ygbA,mutS,pphB,ygbI,ygbJ,ygbK,ygbL,ygbM,ygbN,rpoS,nlpD,pcm,surE,ygbO,ispF,ispD,ftsB,ygbE,cysC,cysN,cysD,iap,insL-5,insK-3,insJ-3	genes_overlapping=fhlA,ygbA,mutS,pphB,ygbI,ygbJ,ygbK,ygbL,ygbM,ygbN,rpoS,nlpD,pcm,surE,ygbO,ispF,ispD,ftsB,ygbE,cysC,cysN,cysD,iap,insL-5,insK-3,insJ-3	locus_tag=[ECB_02581]–[ECB_02606]	locus_tags_overlapping=ECB_02581,ECB_02582,ECB_02583,ECB_02584,ECB_02585,ECB_02586,ECB_02587,ECB_02588,ECB_02589,ECB_02590,ECB_02591,ECB_02592,ECB_02593,ECB_02594,ECB_02595,ECB_02596,ECB_02597,ECB_02598,ECB_02599,ECB_02600,ECB_02601,ECB_02602,ECB_02603,ECB_02604,ECB_02605,ECB_02606	mutation_category=large_amplification	position_end=2775877	position_start=2749275	ref_seq=26603-bp
SNP	18	36	REL606	3107610	A	aa_new_seq=I	aa_position=99	aa_ref_seq=I	codon_new_seq=ATA	codon_number=99	codon_position=3	codon_ref_seq=ATT	gene_name=ygiN	gene_position=297	gene_product=quinol monooxygenase	gene_strand=>	genes_overlapping=ygiN	locus_tag=ECB_02901	locus_tags_overlapping=ECB_02901	mutation_category=snp_synonymous	position_end=3107610	position_start=3107610	ref_seq=T	snp_type=synonymous	transl_table=11
SNP	19	37	REL606	3249544	T	aa_new_seq=I	aa_position=569	aa_ref_seq=V	codon_new_seq=ATT	codon_number=569	codon_position=1	codon_ref_seq=GTT	gene_name=infB	gene_position=1705	gene_product=translation initiation factor IF-2	gene_strand=<	genes_overlapping=infB	locus_tag=ECB_03035	locus_tags_overlapping=ECB_03035	mutation_category=snp_nonsynonymous	position_end=3249544	position_start=3249544	ref_seq=C	snp_type=nonsynonymous	transl_table=11
SNP	20	38	REL606	3289548	A	gene_name=yhcC/gltB	gene_position=intergenic (-263/-319)	gene_product=predicted Fe-S oxidoreductase/glutamate synthase, large subunit	gene_strand=</>	locus_tag=ECB_03076/ECB_03077	mutation_category=snp_intergenic	position_end=3289548	position_start=3289548	ref_seq=G	snp_type=intergenic
SNP	21	39	REL606	3340296	A	aa_new_seq=I	aa_position=251	aa_ref_seq=M	codon_new_seq=ATA	codon_number=251	codon_position=3	codon_ref_seq=ATG	gene_name=yhdJ	gene_position=753	gene_product=predicted methyltransferase	gene_strand=>	genes_overlapping=yhdJ	locus_tag=ECB_03120	locus_tags_overlapping=ECB_03120	mutation_category=snp_nonsynonymous	position_end=3340296	position_start=3340296	ref_seq=G	snp_type=nonsynonymous	transl_table=11
DEL	22	45,59	REL606	3894997	5627	gene_name=rbsD–[rbsR]	gene_product=rbsD,rbsA,rbsC,rbsB,rbsK,[rbsR]	genes_inactivated=rbsD,rbsA,rbsC,rbsB,rbsK,rbsR	locus_tag=[ECB_03634]–[ECB_03639]	locus_tags_inactivated=ECB_03634,ECB_03635,ECB_03636,ECB_03637,ECB_03638,ECB_03639	mediated=IS150	mutation_category=large_deletion	position_end=3900623	position_start=3894997	ref_seq=5627-bp
SNP	23	40	REL606	4100211	T	aa_new_seq=L	aa_position=340	aa_ref_seq=F	codon_new_seq=TTA	codon_number=340	codon_position=3	codon_ref_seq=TTC	gene_name=hslU	gene_position=1020	gene_product=ATP-dependent protease ATP-binding subunit	gene_strand=<	genes_overlapping=hslU	locus_tag=ECB_03816	locus_tags_overlapping=ECB_03816	mutation_category=snp_nonsynonymous	position_end=4100211	position_start=4100211	ref_seq=G	snp_type=nonsynonymous	transl_table=11
SNP	24	41	REL606	4201937	T	aa_new_seq=D	aa_position=208	aa_ref_seq=G	codon_new_seq=GAT	codon_number=208	codon_position=2	codon_ref_seq=GGT	gene_name=iclR	gene_position=623	gene_product=DNA-binding transcriptional repressor	gene_strand=<	genes_overlapping=iclR	locus_tag=ECB_03890	locus_tags_overlapping=ECB_03890	mutation_category=snp_nonsynonymous	position_end=4201937	position_start=4201937	ref_seq=C	snp_type=nonsynonymous	transl_table=11
MOB	25	52,60	REL606	4252526	IS150	1	3	gene_name=uvrA	gene_position=coding (276-278/2823 nt)	gene_product=excinuclease ABC subunit A	gene_strand=<	genes_inactivated=uvrA	locus_tag=ECB_03930	locus_tags_inactivated=ECB_03930	mutation_category=mobile_element_insertion	position_end=4252528	position_start=4252526	ref_seq=TTA	repeat_size=1443
SNP	26	42	REL606	4270684	T	aa_new_seq=L	aa_position=80	aa_ref_seq=L	codon_new_seq=CTT	codon_number=80	codon_position=3	codon_ref_seq=CTG	gene_name=nrfE	gene_position=240	gene_product=heme lyase (NrfEFG) for insertion of heme into c552, subunit NrfE	gene_strand=>	genes_overlapping=nrfE	locus_tag=ECB_03946	locus_tags_overlapping=ECB_03946	mutation_category=snp_synonymous	position_end=4270684	position_start=4270684	ref_seq=G	snp_type=synonymous	transl_table=11
MOB	27	53,61	REL606	4588156	IS150	1	3	gene_name=mdoB	gene_position=coding (518-520/2292 nt)	gene_product=phosphoglycerol transferase I	gene_strand=<	genes_inactivated=mdoB	locus_tag=ECB_04236	locus_tags_inactivated=ECB_04236	mutation_category=mobile_element_insertion	position_end=4588158	position_start=4588156	ref_seq=GCT	repeat_size=1443
RA	28	.	REL606	651601	0	T	A	aa_new_seq=L	aa_position=7	aa_ref_seq=Q	allele_frequencies=A:1.000e+00	codon_new_seq=CTG	codon_number=7	codon_position=2	codon_ref_seq=CAG	fisher_strand_p_value=1.00000e+00	frequency=1.000e+00	frequency_lower=9.852e-01	frequency_upper=1.000e+00	gene_name=ybeB	gene_position=20	gene_product=hypothetical protein	gene_strand=<	ks_quality_p_value=1.00000e+00	locus_tag=ECB_00606	major_base=A	major_cov=52/39	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=52/39	new_seq=A	prediction=consensus	ref_cov=0/0	ref_seq=T	score=329.9	snp_type=nonsynonymous	total_cov=52/39	transl_table=11
RA	29	.	REL606	1089367	0	G	A	aa_new_seq=V	aa_position=310	aa_ref_seq=A	allele_frequencies=A:1.000e+00	codon_new_seq=GTC	codon_number=310	codon_position=2	codon_ref_seq=GCC	fisher_strand_p_value=1.00000e+00	frequency=1.000e+00	frequency_lower=9.742e-01	frequency_upper=1.000e+00	gene_name=ycdM	gene_position=929	gene_product=predicted monooxygenase	gene_strand=<	ks_quality_p_value=1.00000e+00	locus_tag=ECB_01015	major_base=A	major_cov=35/36	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=35/36	new_seq=A	prediction=consensus	ref_cov=0/0	ref_seq=G	score=250.8	snp_type=nonsynonymous	total_cov=36/36	transl_table=11
RA	30	.	REL606	1166950	0	C	T	aa_new_seq=I	aa_position=148	aa_ref_seq=T	allele_frequencies=T:1.000e+00	codon_new_seq=ATT	codon_number=148	codon_position=2	codon_ref_seq=ACT	fisher_strand_p_value=1.00000e+00	frequency=1.000e+00	frequency_lower=9.806e-01	frequency_upper=1.000e+00	gene_name=fabF	gene_position=443	gene_product=3-oxoacyl-(acyl carrier protein) synthase	gene_strand=>	ks_quality_p_value=1.00000e+00	locus_tag=ECB_01091	major_base=T	major_cov=35/34	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=35/34	new_seq=T	prediction=consensus	ref_cov=0/0	ref_seq=C	score=247.1	snp_type=nonsynonymous	total_cov=35/34	transl_table=11
RA	31	.	REL606	1330018	0	C	A	aa_new_seq=E	aa_position=200	aa_ref_seq=A	allele_frequencies=A:1.000e+00	codon_new_seq=GAG	codon_number=200	codon_position=2	codon_ref_seq=GCG	fisher_strand_p_value=1.00000e+00	frequency=1.000e+00	frequency_lower=9.800e-01	frequency_upper=1.000e+00	gene_name=topA	gene_position=599	gene_product=DNA topoisomerase I	gene_strand=>	ks_quality_p_value=1.00000e+00	locus_tag=ECB_01250	major_base=A	major_cov=29/38	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=29/38	new_seq=A	prediction=consensus	ref_cov=0/0	ref_seq=C	score=235.1	snp_type=nonsynonymous	total_cov=29/38	transl_table=11
RA	32	.	REL606	1587635	0	A	G	aa_new_seq=C	aa_position=99	aa_ref_seq=Y	allele_frequencies=G:1.000e+00	codon_new_seq=TGC	codon_number=99	codon_position=2	codon_ref_seq=TAC	fisher_strand_p_value=1.00000e+00	frequency=1.000e+00	frequency_lower=9.828e-01	frequency_upper=1.000e+00	gene_name=marA	gene_position=296	gene_product=DNA-binding transcriptional dual activator of multiple antibiotic resistance	gene_strand=>	ks_quality_p_value=1.00000e+00	locus_tag=ECB_01490	major_base=G	major_cov=38/40	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=38/40	new_seq=G	prediction=consensus	ref_cov=0/0	ref_seq=A	score=321.2	snp_type=nonsynonymous	total_cov=38/40	transl_table=11
RA	33	.	REL606	1733865	0	G	T	aa_new_seq=S	aa_position=301	aa_ref_seq=A	allele_frequencies=G:2.353e-02,T:9.765e-01	codon_new_seq=TCC	codon_number=301	codon_position=1	codon_ref_seq=GCC	fisher_strand_p_value=1.00000e+00	frequency=9.765e-01	frequency_lower=9.388e-01	frequency_upper=9.944e-01	gene_name=pykF	gene_position=901	gene_product=pyruvate kinase	gene_strand=>	ks_quality_p_value=7.25210e-01	locus_tag=ECB_01645	major_base=T	major_cov=36/47	major_frequency=9.765e-01	minor_base=G	minor_cov=1/1	new_cov=36/47	new_seq=T	prediction=consensus	ref_cov=1/1	ref_seq=G	score=291.7	snp_type=nonsynonymous	total_cov=37/48	transl_table=11
RA	34	.	REL606	2275076	0	T	G	aa_new_seq=H	aa_position=474	aa_ref_seq=Q	allele_frequencies=G:1.000e+00	codon_new_seq=CAC	codon_number=474	codon_position=3	codon_ref_seq=CAA	fisher_strand_p_value=1.00000e+00	frequency=1.000e+00	frequency_lower=9.777e-01	frequency_upper=1.000e+00	gene_name=yfaQ	gene_position=1422	gene_product=hypothetical protein	gene_strand=<	ks_quality_p_value=1.00000e+00	locus_tag=ECB_02153	major_base=G	major_cov=29/31	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=29/31	new_seq=G	prediction=consensus	ref_cov=0/0	ref_seq=T	score=246.3	snp_type=nonsynonymous	total_cov=29/31	transl_table=11
RA	35	.	REL606	2766603	0	G	A	aa_new_seq=G	aa_position=13	aa_ref_seq=G	allele_frequencies=A:4.399e-01,G:5.601e-01	codon_new_seq=GGT	codon_number=13	codon_position=3	codon_ref_seq=GGC	consensus_reject=FREQUENCY_CUTOFF	fisher_strand_p_value=8.86816e-01	frequency=4.399e-01	frequency_lower=3.829e-01	frequency_upper=4.980e-01	gene_name=ispF	gene_position=39	gene_product=2-C-methyl-D-erythritol 2,4-cyclodiphosphate synthase	gene_strand=<	ks_quality_p_value=9.60450e-01	locus_tag=ECB_02596	major_base=G	major_cov=54/58	major_frequency=5.601e-01	minor_base=A	minor_cov=44/44	new_cov=44/44	new_seq=A	prediction=polymorphism	ref_cov=54/58	ref_seq=G	score=282.6	snp_type=synonymous	total_cov=98/102	transl_table=11
RA	36	.	REL606	3107610	0	T	A	aa_new_seq=I	aa_position=99	aa_ref_seq=I	allele_frequencies=A:1.000e+00	codon_new_seq=ATA	codon_number=99	codon_position=3	codon_ref_seq=ATT	fisher_strand_p_value=1.00000e+00	frequency=1.000e+00	frequency_lower=9.847e-01	frequency_upper=1.000e+00	gene_name=ygiN	gene_position=297	gene_product=quinol monooxygenase	gene_strand=>	ks_quality_p_value=1.00000e+00	locus_tag=ECB_02901	major_base=A	major_cov=42/46	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=42/46	new_seq=A	prediction=consensus	ref_cov=0/0	ref_seq=T	score=313.3	snp_type=synonymous	total_cov=42/46	transl_table=11
RA	37	.	REL606	3249544	0	C	T	aa_new_seq=I	aa_position=569	aa_ref_seq=V	allele_frequencies=T:1.000e+00	codon_new_seq=ATT	codon_number=569	codon_position=1	codon_ref_seq=GTT	fisher_strand_p_value=1.00000e+00	frequency=1.000e+00	frequency_lower=9.849e-01	frequency_upper=1.000e+00	gene_name=infB	gene_position=1705	gene_product=translation initiation factor IF-2	gene_strand=<	ks_quality_p_value=1.00000e+00	locus_tag=ECB_03035	major_base=T	major_cov=40/49	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=40/49	new_seq=T	prediction=consensus	ref_cov=0/0	ref_seq=C	score=305.7	snp_type=nonsynonymous	total_cov=40/49	transl_table=11
RA	38	.	REL606	3289548	0	G	A	allele_frequencies=A:1.000e+00	fisher_strand_p_value=1.00000e+00	frequency=1.000e+00	frequency_lower=9.896e-01	frequency_upper=1.000e+00	gene_name=yhcC/gltB	gene_position=intergenic (-263/-319)	gene_product=predicted Fe-S oxidoreductase/glutamate synthase, large subunit	gene_strand=</>	ks_quality_p_value=1.00000e+00	locus_tag=ECB_03076/ECB_03077	major_base=A	major_cov=56/73	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=56/73	prediction=consensus	ref_cov=0/0	score=463.9	snp_type=intergenic	total_cov=56/73
RA	39	.	REL606	3340296	0	G	A	aa_new_seq=I	aa_position=251	aa_ref_seq=M	allele_frequencies=A:1.000e+00	codon_new_seq=ATA	codon_number=251	codon_position=3	codon_ref_seq=ATG	fisher_strand_p_value=1.00000e+00	frequency=1.000e+00	frequency_lower=9.846e-01	frequency_upper=1.000e+00	gene_name=yhdJ	gene_position=753	gene_product=predicted methyltransferase	gene_strand=>	ks_quality_p_value=1.00000e+00	locus_tag=ECB_03120	major_base=A	major_cov=38/49	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=38/49	new_seq=A	prediction=consensus	ref_cov=0/0	ref_seq=G	score=313.9	snp_type=nonsynonymous	total_cov=38/49	transl_table=11
RA	40	.	REL606	4100211	0	G	T	aa_new_seq=L	aa_position=340	aa_ref_seq=F	allele_frequencies=T:1.000e+00	codon_new_seq=TTA	codon_number=340	codon_position=3	codon_ref_seq=TTC	fisher_strand_p_value=1.00000e+00	frequency=1.000e+00	frequency_lower=9.880e-01	frequency_upper=1.000e+00	gene_name=hslU	gene_position=1020	gene_product=ATP-dependent protease ATP-binding subunit	gene_strand=<	ks_quality_p_value=1.00000e+00	locus_tag=ECB_03816	major_base=T	major_cov=50/62	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=50/62	new_seq=T	prediction=consensus	ref_cov=0/0	ref_seq=G	score=409.7	snp_type=nonsynonymous	total_cov=50/62	transl_table=11
RA	41	.	REL606	4201937	0	C	T	aa_new_seq=D	aa_position=208	aa_ref_seq=G	allele_frequencies=T:1.000e+00	codon_new_seq=GAT	codon_number=208	codon_position=2	codon_ref_seq=GGT	fisher_strand_p_value=1.00000e+00	frequency=1.000e+00	frequency_lower=9.877e-01	frequency_upper=1.000e+00	gene_name=iclR	gene_position=623	gene_product=DNA-binding transcriptional repressor	gene_strand=<	ks_quality_p_value=1.00000e+00	locus_tag=ECB_03890	major_base=T	major_cov=55/54	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=55/54	new_seq=T	prediction=consensus	ref_cov=0/0	ref_seq=C	score=377.1	snp_type=nonsynonymous	total_cov=55/54	transl_table=11
RA	42	.	REL606	4270684	0	G	T	aa_new_seq=L	aa_position=80	aa_ref_seq=L	allele_frequencies=T:1.000e+00	codon_new_seq=CTT	codon_number=80	codon_position=3	codon_ref_seq=CTG	fisher_strand_p_value=1.00000e+00	frequency=1.000e+00	frequency_lower=9.811e-01	frequency_upper=1.000e+00	gene_name=nrfE	gene_position=240	gene_product=heme lyase (NrfEFG) for insertion of heme into c552, subunit NrfE	gene_strand=>	ks_quality_p_value=1.00000e+00	locus_tag=ECB_03946	major_base=T	major_cov=34/37	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=34/37	new_seq=T	prediction=consensus	ref_cov=0/0	ref_seq=G	score=251.3	snp_type=synonymous	total_cov=34/37	transl_table=11
MC	43	.	REL606	474383	474391	0	0	gene_name=ybaL	gene_position=coding (915-923/1677 nt)	gene_product=predicted transporter with NAD(P)-binding Rossmann-fold domain	gene_strand=<	left_inside_cov=0	left_outside_cov=99	locus_tag=ECB_00429	right_inside_cov=0	right_outside_cov=101
MC	44	.	REL606	2100308	2122453	0	0	gene_name=ogrK–[ECB_02013]	gene_product=ogrK,yegZ,ECB_01989,ECB_01990,ECB_01991,ECB_01992,ECB_01993,ECB_01994,ECB_01995,ECB_01996,ECB_01997,ECB_01998,ECB_01999,ECB_02000,ECB_02001,ECB_02002,ECB_02003,ECB_02004,ECB_02005,ECB_02006,ECB_02007,ECB_02008,ECB_02009,ECB_02010,ECB_02011,ECB_02012,[ECB_02013]	left_inside_cov=0	left_outside_cov=90	locus_tag=[ECB_01987]–[ECB_02013]	right_inside_cov=23	right_outside_cov=111
MC	45	.	REL606	3894931	3900623	65	0	gene_name=[insK-5]–[rbsR]	gene_product=[insK-5],rbsD,rbsA,rbsC,rbsB,rbsK,[rbsR]	left_inside_cov=42	left_outside_cov=43	locus_tag=[ECB_03633]–[ECB_03639]	right_inside_cov=1	right_outside_cov=135
JC	46	.	REL606	1	1	REL606	4629812	-1	0	alignment_overlap=0	coverage_minus=51	coverage_plus=74	flanking_left=101	flanking_right=101	frequency=NA	frequency_lower=NA	frequency_upper=NA	ignore=CIRCULAR_CHROMOSOME	junction_effective_depth=116.00	junction_possible_overlap_registers=98	junction_possible_overlap_registers_before_trimming=100	key=REL606__1__1__REL606__4629812__-1__0____101__101__0__0	max_left=100	max_left_minus=100	max_left_plus=100	max_min_left=50	max_min_left_minus=50	max_min_left_plus=50	max_min_right=49	max_min_right_minus=48	max_min_right_plus=49	max_pos_hash_score=200	max_right=100	max_right_minus=100	max_right_plus=100	neg_log10_pos_hash_p_value=0.1	new_junction_coverage=1.17	new_junction_read_count=116	new_junction_reference_weighted_read_count=0.01	new_junction_weighted_read_count=115.99	pos_hash_score=88	prediction=unknown	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=NA	side_1_gene_name=–/thrL	side_1_gene_position=intergenic (–/-189)	side_1_gene_product=–/thr operon leader peptide	side_1_gene_strand=–/>	side_1_locus_tag=–/ECB_00001	side_1_overlap=0	side_1_possible_overlap_registers=0	side_1_possible_overlap_registers_before_trimming=100	side_1_read_count=0	side_1_redundant=0	side_1_weighted_read_count=0.00	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=NA	side_2_gene_name=lasT/–	side_2_gene_position=intergenic (+24/–)	side_2_gene_product=predicted rRNA methyltransferase/–	side_2_gene_strand=>/–	side_2_locus_tag=ECB_04279/–	side_2_overlap=0	side_2_possible_overlap_registers=0	side_2_possible_overlap_registers_before_trimming=100	side_2_read_count=0	side_2_redundant=0	side_2_weighted_read_count=0.00	total_non_overlap_reads=125
JC	47	.	REL606	467507	1	REL606	664688	1	0	alignment_overlap=1	coverage_minus=54	coverage_plus=65	flanking_left=101	flanking_right=101	frequency=1.000e+00	frequency_lower=9.751e-01	frequency_upper=1.000e+00	junction_effective_depth=119.00	junction_mixture_iterations=1	junction_possible_overlap_registers=97	junction_possible_overlap_registers_before_trimming=99	key=REL606__467506__1__REL606__2775877__-1__1____101__101__0__1	max_left=98	max_left_minus=95	max_left_plus=98	max_min_left=49	max_min_left_minus=45	max_min_left_plus=49	max_min_right=50	max_min_right_minus=50	max_min_right_plus=47	max_pos_hash_score=198	max_right=97	max_right_minus=97	max_right_plus=97	neg_log10_pos_hash_p_value=0.1	new_junction_coverage=1.21	new_junction_read_count=119	new_junction_reference_weighted_read_count=0.00	new_junction_weighted_read_count=119.00	pos_hash_score=82	prediction=consensus	read_count_offset=3	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=0.00	side_1_gene_name=htpG	side_1_gene_position=coding (321/1875 nt)	side_1_gene_product=heat shock protein 90	side_1_gene_strand=>	side_1_locus_tag=ECB_00424	side_1_overlap=0	side_1_possible_overlap_registers=93	side_1_possible_overlap_registers_before_trimming=96	side_1_read_count=0	side_1_redundant=0	side_1_weighted_read_count=0.00	side_2_annotate_key=repeat	side_2_continuation=0	side_2_coverage=NA	side_2_gene_name=IS150	side_2_gene_position=noncoding (1/1443 nt)	side_2_gene_product=repeat region	side_2_gene_strand=>	side_2_overlap=1	side_2_possible_overlap_registers=NA	side_2_possible_overlap_registers_before_trimming=NA	side_2_read_count=NA	side_2_redundant=1	side_2_weighted_read_count=NA	total_non_overlap_reads=119
JC	48	.	REL606	467509	-1	REL606	666130	-1	0	alignment_overlap=0	coverage_minus=67	coverage_plus=70	flanking_left=101	flanking_right=101	frequency=1.000e+00	frequency_lower=9.779e-01	frequency_upper=1.000e+00	junction_effective_depth=134.00	junction_mixture_iterations=1	junction_possible_overlap_registers=98	junction_possible_overlap_registers_before_trimming=100	key=REL606__467509__-1__REL606__590471__-1__0____101__101__0__1	max_left=100	max_left_minus=99	max_left_plus=100	max_min_left=49	max_min_left_minus=48	max_min_left_plus=49	max_min_right=50	max_min_right_minus=50	max_min_right_plus=50	max_pos_hash_score=200	max_right=98	max_right_minus=98	max_right_plus=97	neg_log10_pos_hash_p_value=0.0	new_junction_coverage=1.35	new_junction_read_count=134	new_junction_reference_weighted_read_count=0.00	new_junction_weighted_read_count=134.00	pos_hash_score=92	prediction=consensus	read_count_offset=3	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=0.00	side_1_gene_name=htpG	side_1_gene_position=coding (323/1875 nt)	side_1_gene_product=heat shock protein 90	side_1_gene_strand=>	side_1_locus_tag=ECB_00424	side_1_overlap=0	side_1_possible_overlap_registers=94	side_1_possible_overlap_registers_before_trimming=97	side_1_read_count=0	side_1_redundant=0	side_1_weighted_read_count=0.00	side_2_annotate_key=repeat	side_2_continuation=0	side_2_coverage=NA	side_2_gene_name=IS150	side_2_gene_position=noncoding (1443/1443 nt)	side_2_gene_product=repeat region	side_2_gene_strand=>	side_2_overlap=0	side_2_possible_overlap_registers=NA	side_2_possible_overlap_registers_before_trimming=NA	side_2_read_count=NA	side_2_redundant=1	side_2_weighted_read_count=NA	total_non_overlap_reads=137
JC	49	.	REL606	474382	-1	REL606	474392	1	0	alignment_overlap=8	coverage_minus=49	coverage_plus=48	flanking_left=101	flanking_right=101	frequency=1.000e+00	frequency_lower=9.676e-01	frequency_upper=1.000e+00	junction_effective_depth=91.00	junction_mixture_iterations=1	junction_possible_overlap_registers=86	junction_possible_overlap_registers_before_trimming=92	key=REL606__474382__-1__REL606__474384__1__8____101__101__0__0	max_left=91	max_left_minus=86	max_left_plus=91	max_min_left=46	max_min_left_minus=46	max_min_left_plus=43	max_min_right=46	max_min_right_minus=42	max_min_right_plus=46	max_pos_hash_score=184	max_right=91	max_right_minus=90	max_right_plus=91	neg_log10_pos_hash_p_value=0.2	new_junction_coverage=1.04	new_junction_read_count=91	new_junction_reference_weighted_read_count=0.00	new_junction_weighted_read_count=91.00	pos_hash_score=72	prediction=consensus	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=0.00	side_1_gene_name=ybaL	side_1_gene_position=coding (924/1677 nt)	side_1_gene_product=predicted transporter with NAD(P)-binding Rossmann-fold domain	side_1_gene_strand=<	side_1_locus_tag=ECB_00429	side_1_overlap=8	side_1_possible_overlap_registers=98	side_1_possible_overlap_registers_before_trimming=100	side_1_read_count=0	side_1_redundant=0	side_1_weighted_read_count=0.00	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=0.00	side_2_gene_name=ybaL	side_2_gene_position=coding (914/1677 nt)	side_2_gene_product=predicted transporter with NAD(P)-binding Rossmann-fold domain	side_2_gene_strand=<	side_2_locus_tag=ECB_00429	side_2_overlap=0	side_2_possible_overlap_registers=90	side_2_possible_overlap_registers_before_trimming=92	side_2_read_count=0	side_2_redundant=0	side_2_weighted_read_count=0.00	total_non_overlap_reads=97
JC	50	.	REL606	664688	1	REL606	1462251	1	0	alignment_overlap=0	coverage_minus=46	coverage_plus=49	flanking_left=101	flanking_right=101	frequency=1.000e+00	frequency_lower=9.686e-01	frequency_upper=1.000e+00	junction_effective_depth=94.00	junction_mixture_iterations=1	junction_possible_overlap_registers=98	junction_possible_overlap_registers_before_trimming=100	key=REL606__664688__1__REL606__1462251__1__0____101__101__1__0	max_left=98	max_left_minus=98	max_left_plus=98	max_min_left=49	max_min_left_minus=49	max_min_left_plus=49	max_min_right=50	max_min_right_minus=50	max_min_right_plus=49	max_pos_hash_score=200	max_right=100	max_right_minus=97	max_right_plus=100	neg_log10_pos_hash_p_value=0.5	new_junction_coverage=0.95	new_junction_read_count=94	new_junction_reference_weighted_read_count=0.00	new_junction_weighted_read_count=94.00	pos_hash_score=67	prediction=consensus	read_count_offset=3	side_1_annotate_key=repeat	side_1_continuation=0	side_1_coverage=NA	side_1_gene_name=IS150	side_1_gene_position=noncoding (1/1443 nt)	side_1_gene_product=repeat region	side_1_gene_strand=>	side_1_overlap=0	side_1_possible_overlap_registers=NA	side_1_possible_overlap_registers_before_trimming=NA	side_1_read_count=NA	side_1_redundant=1	side_1_weighted_read_count=NA	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=0.00	side_2_gene_name=mokB/trg	side_2_gene_position=intergenic (-13/-328)	side_2_gene_product=regulatory peptide/methyl-accepting chemotaxis protein III, ribose and galactose sensor receptor	side_2_gene_strand=</>	side_2_locus_tag=ECB_01377/ECB_01378	side_2_overlap=0	side_2_possible_overlap_registers=95	side_2_possible_overlap_registers_before_trimming=97	side_2_read_count=0	side_2_redundant=0	side_2_weighted_read_count=0.00	total_non_overlap_reads=95
JC	51	.	REL606	664688	1	REL606	2157883	-1	0	alignment_overlap=0	coverage_minus=47	coverage_plus=35	flanking_left=101	flanking_right=101	frequency=1.000e+00	frequency_lower=9.632e-01	frequency_upper=1.000e+00	junction_effective_depth=80.00	junction_mixture_iterations=1	junction_possible_overlap_registers=98	junction_possible_overlap_registers_before_trimming=100	key=REL606__664688__1__REL606__2157883__-1__0____101__101__1__0	max_left=95	max_left_minus=95	max_left_plus=91	max_min_left=48	max_min_left_minus=48	max_min_left_plus=47	max_min_right=49	max_min_right_minus=49	max_min_right_plus=48	max_pos_hash_score=200	max_right=100	max_right_minus=100	max_right_plus=100	neg_log10_pos_hash_p_value=0.7	new_junction_coverage=0.80	new_junction_read_count=80	new_junction_reference_weighted_read_count=0.00	new_junction_weighted_read_count=80.00	pos_hash_score=63	prediction=consensus	read_count_offset=3	side_1_annotate_key=repeat	side_1_continuation=0	side_1_coverage=NA	side_1_gene_name=IS150	side_1_gene_position=noncoding (1/1443 nt)	side_1_gene_product=repeat region	side_1_gene_strand=>	side_1_overlap=0	side_1_possible_overlap_registers=NA	side_1_possible_overlap_registers_before_trimming=NA	side_1_read_count=NA	side_1_redundant=1	side_1_weighted_read_count=NA	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=0.00	side_2_gene_name=yehM	side_2_gene_position=coding (993/2280 nt)	side_2_gene_product=hypothetical protein	side_2_gene_strand=>	side_2_locus_tag=ECB_02049	side_2_overlap=0	side_2_possible_overlap_registers=92	side_2_possible_overlap_registers_before_trimming=97	side_2_read_count=0	side_2_redundant=0	side_2_weighted_read_count=0.00	total_non_overlap_reads=82
JC	52	.	REL606	664688	1	REL606	4252528	-1	0	alignment_overlap=2	coverage_minus=56	coverage_plus=49	flanking_left=101	flanking_right=101	frequency=1.000e+00	frequency_lower=9.719e-01	frequency_upper=1.000e+00	junction_effective_depth=105.00	junction_mixture_iterations=1	junction_possible_overlap_registers=96	junction_possible_overlap_registers_before_trimming=98	key=REL606__3893554__1__REL606__4252530__-1__2____101__101__1__0	max_left=96	max_left_minus=94	max_left_plus=96	max_min_left=49	max_min_left_minus=47	max_min_left_plus=49	max_min_right=49	max_min_right_minus=49	max_min_right_plus=48	max_pos_hash_score=196	max_right=97	max_right_minus=97	max_right_plus=97	neg_log10_pos_hash_p_value=0.2	new_junction_coverage=1.08	new_junction_read_count=105	new_junction_reference_weighted_read_count=0.00	new_junction_weighted_read_count=105.00	pos_hash_score=77	prediction=consensus	read_count_offset=3	side_1_annotate_key=repeat	side_1_continuation=0	side_1_coverage=NA	side_1_gene_name=IS150	side_1_gene_position=noncoding (1/1443 nt)	side_1_gene_product=repeat region	side_1_gene_strand=>	side_1_overlap=2	side_1_possible_overlap_registers=NA	side_1_possible_overlap_registers_before_trimming=NA	side_1_read_count=NA	side_1_redundant=1	side_1_weighted_read_count=NA	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=0.00	side_2_gene_name=uvrA	side_2_gene_position=coding (276/2823 nt)	side_2_gene_product=excinuclease ABC subunit A	side_2_gene_strand=<	side_2_locus_tag=ECB_03930	side_2_overlap=0	side_2_possible_overlap_registers=89	side_2_possible_overlap_registers_before_trimming=95	side_2_read_count=0	side_2_redundant=0	side_2_weighted_read_count=0.00	total_non_overlap_reads=105
JC	53	.	REL606	664688	1	REL606	4588158	-1	0	alignment_overlap=0	coverage_minus=57	coverage_plus=48	flanking_left=101	flanking_right=101	frequency=1.000e+00	frequency_lower=9.716e-01	frequency_upper=1.000e+00	junction_effective_depth=104.00	junction_mixture_iterations=1	junction_possible_overlap_registers=98	junction_possible_overlap_registers_before_trimming=100	key=REL606__3893554__1__REL606__4588158__-1__0____101__101__1__0	max_left=99	max_left_minus=96	max_left_plus=99	max_min_left=50	max_min_left_minus=49	max_min_left_plus=50	max_min_right=49	max_min_right_minus=49	max_min_right_plus=49	max_pos_hash_score=200	max_right=100	max_right_minus=98	max_right_plus=100	neg_log10_pos_hash_p_value=0.2	new_junction_coverage=1.05	new_junction_read_count=104	new_junction_reference_weighted_read_count=0.00	new_junction_weighted_read_count=104.00	pos_hash_score=79	prediction=consensus	read_count_offset=3	side_1_annotate_key=repeat	side_1_continuation=0	side_1_coverage=NA	side_1_gene_name=IS150	side_1_gene_position=noncoding (1/1443 nt)	side_1_gene_product=repeat region	side_1_gene_strand=>	side_1_overlap=0	side_1_possible_overlap_registers=NA	side_1_possible_overlap_registers_before_trimming=NA	side_1_read_count=NA	side_1_redundant=1	side_1_weighted_read_count=NA	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=0.00	side_2_gene_name=mdoB	side_2_gene_position=coding (518/2292 nt)	side_2_gene_product=phosphoglycerol transferase I	side_2_gene_strand=<	side_2_locus_tag=ECB_04236	side_2_overlap=0	side_2_possible_overlap_registers=94	side_2_possible_overlap_registers_before_trimming=97	side_2_read_count=0	side_2_redundant=0	side_2_weighted_read_count=0.00	total_non_overlap_reads=105
JC	54	.	REL606	665214	1	REL606	15784	-1	0	alignment_overlap=1	coverage_minus=10	coverage_plus=16	flanking_left=101	flanking_right=101	frequency=NA	frequency_lower=NA	frequency_upper=NA	junction_effective_depth=26.00	junction_possible_overlap_registers=97	junction_possible_overlap_registers_before_trimming=99	key=REL606__589555__1__REL606__591114__-1__1____101__101__1__1	max_left=95	max_left_minus=93	max_left_plus=95	max_min_left=43	max_min_left_minus=31	max_min_left_plus=43	max_min_right=45	max_min_right_minus=35	max_min_right_plus=45	max_pos_hash_score=198	max_right=97	max_right_minus=92	max_right_plus=97	neg_log10_pos_hash_p_value=3.8	new_junction_coverage=0.26	new_junction_read_count=26	new_junction_reference_weighted_read_count=0.00	new_junction_weighted_read_count=26.00	pos_hash_score=24	prediction=unknown	reject=COVERAGE_EVENNESS_SKEW	side_1_annotate_key=repeat	side_1_continuation=0	side_1_coverage=NA	side_1_gene_name=IS150	side_1_gene_position=noncoding (527/1443 nt)	side_1_gene_product=repeat region	side_1_gene_strand=>	side_1_overlap=1	side_1_possible_overlap_registers=NA	side_1_possible_overlap_registers_before_trimming=NA	side_1_read_count=NA	side_1_redundant=1	side_1_weighted_read_count=NA	side_2_annotate_key=repeat	side_2_continuation=0	side_2_coverage=NA	side_2_gene_name=IS186	side_2_gene_position=noncoding (399/1343 nt)	side_2_gene_product=repeat region	side_2_gene_strand=>	side_2_overlap=0	side_2_possible_overlap_registers=NA	side_2_possible_overlap_registers_before_trimming=NA	side_2_read_count=NA	side_2_redundant=1	side_2_weighted_read_count=NA	total_non_overlap_reads=26
JC	55	.	REL606	665981	-1	REL606	15777	1	0	alignment_overlap=1	coverage_minus=10	coverage_plus=4	flanking_left=101	flanking_right=101	frequency=NA	frequency_lower=NA	frequency_upper=NA	junction_effective_depth=14.00	junction_possible_overlap_registers=97	junction_possible_overlap_registers_before_trimming=99	key=REL606__590322__-1__REL606__591105__1__1____101__101__1__1	max_left=86	max_left_minus=86	max_left_plus=61	max_min_left=46	max_min_left_minus=46	max_min_left_plus=45	max_min_right=44	max_min_right_minus=44	max_min_right_plus=36	max_pos_hash_score=198	max_right=80	max_right_minus=69	max_right_plus=80	neg_log10_pos_hash_p_value=4.8	new_junction_coverage=0.14	new_junction_read_count=14	new_junction_reference_weighted_read_count=0.00	new_junction_weighted_read_count=14.00	pos_hash_score=14	prediction=unknown	reject=COVERAGE_EVENNESS_SKEW	side_1_annotate_key=repeat	side_1_continuation=0	side_1_coverage=NA	side_1_gene_name=IS150	side_1_gene_position=noncoding (1294/1443 nt)	side_1_gene_product=repeat region	side_1_gene_strand=>	side_1_overlap=1	side_1_possible_overlap_registers=NA	side_1_possible_overlap_registers_before_trimming=NA	side_1_read_count=NA	side_1_redundant=1	side_1_weighted_read_count=NA	side_2_annotate_key=repeat	side_2_continuation=0	side_2_coverage=NA	side_2_gene_name=IS186	side_2_gene_position=noncoding (392/1343 nt)	side_2_gene_product=repeat region	side_2_gene_strand=>	side_2_overlap=0	side_2_possible_overlap_registers=NA	side_2_possible_overlap_registers_before_trimming=NA	side_2_read_count=NA	side_2_redundant=1	side_2_weighted_read_count=NA	total_non_overlap_reads=14
JC	56	.	REL606	666130	-1	REL606	968328	-1	0	alignment_overlap=0	coverage_minus=55	coverage_plus=55	flanking_left=101	flanking_right=101	frequency=1.000e+00	frequency_lower=9.729e-01	frequency_upper=1.000e+00	junction_effective_depth=109.00	junction_mixture_iterations=1	junction_possible_overlap_registers=97	junction_possible_overlap_registers_before_trimming=100	key=REL606__666130__-1__REL606__968328__-1__0____101__101__1__0	max_left=98	max_left_minus=98	max_left_plus=98	max_min_left=49	max_min_left_minus=49	max_min_left_plus=49	max_min_right=50	max_min_right_minus=49	max_min_right_plus=50	max_pos_hash_score=200	max_right=99	max_right_minus=99	max_right_plus=99	neg_log10_pos_hash_p_value=0.1	new_junction_coverage=1.11	new_junction_read_count=109	new_junction_reference_weighted_read_count=0.00	new_junction_weighted_read_count=109.00	pos_hash_score=81	prediction=consensus	read_count_offset=3	side_1_annotate_key=repeat	side_1_continuation=0	side_1_coverage=NA	side_1_gene_name=IS150	side_1_gene_position=noncoding (1443/1443 nt)	side_1_gene_product=repeat region	side_1_gene_strand=>	side_1_overlap=0	side_1_possible_overlap_registers=NA	side_1_possible_overlap_registers_before_trimming=NA	side_1_read_count=NA	side_1_redundant=1	side_1_weighted_read_count=NA	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=0.00	side_2_gene_name=pflA/pflB	side_2_gene_position=intergenic (-155/+37)	side_2_gene_product=pyruvate formate lyase activating enzyme 1/pyruvate formate lyase I	side_2_gene_strand=</<	side_2_locus_tag=ECB_00906/ECB_00907	side_2_overlap=0	side_2_possible_overlap_registers=93	side_2_possible_overlap_registers_before_trimming=97	side_2_read_count=0	side_2_redundant=0	side_2_weighted_read_count=0.00	total_non_overlap_reads=110
JC	57	.	REL606	666130	-1	REL606	1462253	-1	0	alignment_overlap=0	coverage_minus=51	coverage_plus=60	flanking_left=101	flanking_right=101	frequency=1.000e+00	frequency_lower=9.729e-01	frequency_upper=1.000e+00	junction_effective_depth=109.00	junction_mixture_iterations=1	junction_possible_overlap_registers=97	junction_possible_overlap_registers_before_trimming=100	key=REL606__666130__-1__REL606__1462253__-1__0____101__101__1__0	max_left=97	max_left_minus=97	max_left_plus=97	max_min_left=48	max_min_left_minus=48	max_min_left_plus=43	max_min_right=50	max_min_right_minus=50	max_min_right_plus=50	max_pos_hash_score=200	max_right=100	max_right_minus=100	max_right_plus=100	neg_log10_pos_hash_p_value=0.1	new_junction_coverage=1.11	new_junction_read_count=109	new_junction_reference_weighted_read_count=0.00	new_junction_weighted_read_count=109.00	pos_hash_score=83	prediction=consensus	read_count_offset=3	side_1_annotate_key=repeat	side_1_continuation=0	side_1_coverage=NA	side_1_gene_name=IS150	side_1_gene_position=noncoding (1443/1443 nt)	side_1_gene_product=repeat region	side_1_gene_strand=>	side_1_overlap=0	side_1_possible_overlap_registers=NA	side_1_possible_overlap_registers_before_trimming=NA	side_1_read_count=NA	side_1_redundant=1	side_1_weighted_read_count=NA	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=0.00	side_2_gene_name=mokB/trg	side_2_gene_position=intergenic (-15/-326)	side_2_gene_product=regulatory peptide/methyl-accepting chemotaxis protein III, ribose and galactose sensor receptor	side_2_gene_strand=</>	side_2_locus_tag=ECB_01377/ECB_01378	side_2_overlap=0	side_2_possible_overlap_registers=95	side_2_possible_overlap_registers_before_trimming=97	side_2_read_count=0	side_2_redundant=0	side_2_weighted_read_count=0.00	total_non_overlap_reads=111
JC	58	.	REL606	666130	-1	REL606	2157881	1	0	alignment_overlap=0	coverage_minus=61	coverage_plus=61	flanking_left=101	flanking_right=101	frequency=1.000e+00	frequency_lower=9.753e-01	frequency_upper=1.000e+00	junction_effective_depth=120.00	junction_mixture_iterations=1	junction_possible_overlap_registers=97	junction_possible_overlap_registers_before_trimming=100	key=REL606__666130__-1__REL606__2157881__1__0____101__101__1__0	max_left=97	max_left_minus=97	max_left_plus=97	max_min_left=50	max_min_left_minus=50	max_min_left_plus=50	max_min_right=50	max_min_right_minus=50	max_min_right_plus=48	max_pos_hash_score=200	max_right=100	max_right_minus=100	max_right_plus=98	neg_log10_pos_hash_p_value=0.1	new_junction_coverage=1.22	new_junction_read_count=120	new_junction_reference_weighted_read_count=0.00	new_junction_weighted_read_count=120.00	pos_hash_score=82	prediction=consensus	read_count_offset=3	side_1_annotate_key=repeat	side_1_continuation=0	side_1_coverage=NA	side_1_gene_name=IS150	side_1_gene_position=noncoding (1443/1443 nt)	side_1_gene_product=repeat region	side_1_gene_strand=>	side_1_overlap=0	side_1_possible_overlap_registers=NA	side_1_possible_overlap_registers_before_trimming=NA	side_1_read_count=NA	side_1_redundant=1	side_1_weighted_read_count=NA	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=0.00	side_2_gene_name=yehM	side_2_gene_position=coding (991/2280 nt)	side_2_gene_product=hypothetical protein	side_2_gene_strand=>	side_2_locus_tag=ECB_02049	side_2_overlap=0	side_2_possible_overlap_registers=92	side_2_possible_overlap_registers_before_trimming=97	side_2_read_count=0	side_2_redundant=0	side_2_weighted_read_count=0.00	total_non_overlap_reads=122
JC	59	.	REL606	666130	-1	REL606	3900624	1	0	alignment_overlap=1	coverage_minus=56	coverage_plus=74	flanking_left=101	flanking_right=101	frequency=9.922e-01	frequency_lower=9.639e-01	frequency_upper=9.996e-01	junction_effective_depth=130.00	junction_mixture_iterations=1	junction_possible_overlap_registers=97	junction_possible_overlap_registers_before_trimming=99	key=REL606__3894996__-1__REL606__3900623__1__1____101__101__1__0	max_left=98	max_left_minus=98	max_left_plus=97	max_min_left=49	max_min_left_minus=48	max_min_left_plus=49	max_min_right=50	max_min_right_minus=50	max_min_right_plus=49	max_pos_hash_score=198	max_right=98	max_right_minus=85	max_right_plus=98	neg_log10_pos_hash_p_value=0.0	new_junction_coverage=1.31	new_junction_read_count=129	new_junction_reference_weighted_read_count=0.00	new_junction_weighted_read_count=129.00	pos_hash_score=97	prediction=consensus	side_1_annotate_key=repeat	side_1_continuation=0	side_1_coverage=NA	side_1_gene_name=IS150	side_1_gene_position=noncoding (1443/1443 nt)	side_1_gene_product=repeat region	side_1_gene_strand=>	side_1_overlap=1	side_1_possible_overlap_registers=NA	side_1_possible_overlap_registers_before_trimming=NA	side_1_read_count=NA	side_1_redundant=1	side_1_weighted_read_count=NA	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=0.01	side_2_gene_name=rbsR	side_2_gene_position=coding (591/993 nt)	side_2_gene_product=DNA-binding transcriptional repressor of ribose metabolism	side_2_gene_strand=>	side_2_locus_tag=ECB_03639	side_2_overlap=0	side_2_possible_overlap_registers=96	side_2_possible_overlap_registers_before_trimming=99	side_2_read_count=1	side_2_redundant=0	side_2_weighted_read_count=1.00	total_non_overlap_reads=130
JC	60	.	REL606	666130	-1	REL606	4252526	1	0	alignment_overlap=0	coverage_minus=53	coverage_plus=74	flanking_left=101	flanking_right=101	frequency=1.000e+00	frequency_lower=9.751e-01	frequency_upper=1.000e+00	junction_effective_depth=119.00	junction_mixture_iterations=1	junction_possible_overlap_registers=97	junction_possible_overlap_registers_before_trimming=100	key=REL606__3894996__-1__REL606__4252526__1__0____101__101__1__0	max_left=98	max_left_minus=95	max_left_plus=98	max_min_left=49	max_min_left_minus=49	max_min_left_plus=46	max_min_right=50	max_min_right_minus=48	max_min_right_plus=50	max_pos_hash_score=200	max_right=100	max_right_minus=100	max_right_plus=94	neg_log10_pos_hash_p_value=0.1	new_junction_coverage=1.21	new_junction_read_count=119	new_junction_reference_weighted_read_count=0.00	new_junction_weighted_read_count=119.00	pos_hash_score=81	prediction=consensus	read_count_offset=3	side_1_annotate_key=repeat	side_1_continuation=0	side_1_coverage=NA	side_1_gene_name=IS150	side_1_gene_position=noncoding (1443/1443 nt)	side_1_gene_product=repeat region	side_1_gene_strand=>	side_1_overlap=0	side_1_possible_overlap_registers=NA	side_1_possible_overlap_registers_before_trimming=NA	side_1_read_count=NA	side_1_redundant=1	side_1_weighted_read_count=NA	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=0.00	side_2_gene_name=uvrA	side_2_gene_position=coding (278/2823 nt)	side_2_gene_product=excinuclease ABC subunit A	side_2_gene_strand=<	side_2_locus_tag=ECB_03930	side_2_overlap=0	side_2_possible_overlap_registers=93	side_2_possible_overlap_registers_before_trimming=97	side_2_read_count=0	side_2_redundant=0	side_2_weighted_read_count=0.00	total_non_overlap_reads=127
JC	61	.	REL606	666130	-1	REL606	4588156	1	0	alignment_overlap=0	coverage_minus=53	coverage_plus=85	flanking_left=101	flanking_right=101	frequency=1.000e+00	frequency_lower=9.785e-01	frequency_upper=1.000e+00	junction_effective_depth=138.00	junction_mixture_iterations=1	junction_possible_overlap_registers=98	junction_possible_overlap_registers_before_trimming=100	key=REL606__3652533__-1__REL606__4588156__1__0____101__101__1__0	max_left=99	max_left_minus=99	max_left_plus=98	max_min_left=50	max_min_left_minus=49	max_min_left_plus=50	max_min_right=50	max_min_right_minus=44	max_min_right_plus=50	max_pos_hash_score=200	max_right=98	max_right_minus=98	max_right_plus=95	neg_log10_pos_hash_p_value=0.1	new_junction_coverage=1.39	new_junction_read_count=138	new_junction_reference_weighted_read_count=0.00	new_junction_weighted_read_count=138.00	pos_hash_score=87	prediction=consensus	read_count_offset=3	side_1_annotate_key=repeat	side_1_continuation=0	side_1_coverage=NA	side_1_gene_name=IS150	side_1_gene_position=noncoding (1443/1443 nt)	side_1_gene_product=repeat region	side_1_gene_strand=>	side_1_overlap=0	side_1_possible_overlap_registers=NA	side_1_possible_overlap_registers_before_trimming=NA	side_1_read_count=NA	side_1_redundant=1	side_1_weighted_read_count=NA	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=0.00	side_2_gene_name=mdoB	side_2_gene_position=coding (520/2292 nt)	side_2_gene_product=phosphoglycerol transferase I	side_2_gene_strand=<	side_2_locus_tag=ECB_04236	side_2_overlap=0	side_2_possible_overlap_registers=94	side_2_possible_overlap_registers_before_trimming=97	side_2_read_count=0	side_2_redundant=0	side_2_weighted_read_count=0.00	total_non_overlap_reads=138
JC	62	.	REL606	968326	1	REL606	664688	1	0	alignment_overlap=0	coverage_minus=60	coverage_plus=64	flanking_left=101	flanking_right=101	frequency=1.000e+00	frequency_lower=9.759e-01	frequency_upper=1.000e+00	junction_effective_depth=123.00	junction_mixture_iterations=1	junction_possible_overlap_registers=98	junction_possible_overlap_registers_before_trimming=100	key=REL606__968326__1__REL606__2775877__-1__0____101__101__0__1	max_left=100	max_left_minus=100	max_left_plus=99	max_min_left=50	max_min_left_minus=50	max_min_left_plus=42	max_min_right=50	max_min_right_minus=50	max_min_right_plus=49	max_pos_hash_score=200	max_right=99	max_right_minus=99	max_right_plus=99	neg_log10_pos_hash_p_value=0.1	new_junction_coverage=1.24	new_junction_read_count=123	new_junction_reference_weighted_read_count=0.00	new_junction_weighted_read_count=123.00	pos_hash_score=88	prediction=consensus	read_count_offset=3	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=0.00	side_1_gene_name=pflA/pflB	side_1_gene_position=intergenic (-153/+39)	side_1_gene_product=pyruvate formate lyase activating enzyme 1/pyruvate formate lyase I	side_1_gene_strand=</<	side_1_locus_tag=ECB_00906/ECB_00907	side_1_overlap=0	side_1_possible_overlap_registers=93	side_1_possible_overlap_registers_before_trimming=97	side_1_read_count=0	side_1_redundant=0	side_1_weighted_read_count=0.00	side_2_annotate_key=repeat	side_2_continuation=0	side_2_coverage=NA	side_2_gene_name=IS150	side_2_gene_position=noncoding (1/1443 nt)	side_2_gene_product=repeat region	side_2_gene_strand=>	side_2_overlap=0	side_2_possible_overlap_registers=NA	side_2_possible_overlap_registers_before_trimming=NA	side_2_read_count=NA	side_2_redundant=1	side_2_weighted_read_count=NA	total_non_overlap_reads=124
JC	63	.	REL606	1305934	1	REL606	1193615	1	0	alignment_overlap=0	coverage_minus=44	coverage_plus=41	flanking_left=101	flanking_right=101	frequency=1.000e+00	frequency_lower=9.650e-01	frequency_upper=1.000e+00	junction_effective_depth=84.00	junction_mixture_iterations=1	junction_possible_overlap_registers=97	junction_possible_overlap_registers_before_trimming=100	key=REL606__1305934__1__REL606__1322721__1__0____101__101__0__1	max_left=99	max_left_minus=93	max_left_plus=99	max_min_left=50	max_min_left_minus=50	max_min_left_plus=45	max_min_right=50	max_min_right_minus=49	max_min_right_plus=50	max_pos_hash_score=200	max_right=97	max_right_minus=96	max_right_plus=97	neg_log10_pos_hash_p_value=0.5	new_junction_coverage=0.85	new_junction_read_count=84	new_junction_reference_weighted_read_count=0.00	new_junction_weighted_read_count=84.00	pos_hash_score=69	prediction=consensus	read_count_offset=9	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=0.00	side_1_gene_name=cls	side_1_gene_position=coding (306/1461 nt)	side_1_gene_product=cardiolipin synthetase	side_1_gene_strand=<	side_1_locus_tag=ECB_01223	side_1_overlap=0	side_1_possible_overlap_registers=85	side_1_possible_overlap_registers_before_trimming=91	side_1_read_count=0	side_1_redundant=0	side_1_weighted_read_count=0.00	side_2_annotate_key=repeat	side_2_continuation=0	side_2_coverage=NA	side_2_gene_name=IS1	side_2_gene_position=noncoding (1/768 nt)	side_2_gene_product=repeat region	side_2_gene_strand=>	side_2_overlap=0	side_2_possible_overlap_registers=NA	side_2_possible_overlap_registers_before_trimming=NA	side_2_read_count=NA	side_2_redundant=1	side_2_weighted_read_count=NA	total_non_overlap_reads=85
JC	64	.	REL606	1305942	-1	REL606	1194382	-1	0	alignment_overlap=1	coverage_minus=44	coverage_plus=41	flanking_left=101	flanking_right=101	frequency=1.000e+00	frequency_lower=9.654e-01	frequency_upper=1.000e+00	junction_effective_depth=85.00	junction_mixture_iterations=1	junction_possible_overlap_registers=97	junction_possible_overlap_registers_before_trimming=99	key=REL606__1305943__-1__REL606__1323488__-1__1____101__101__0__1	max_left=98	max_left_minus=98	max_left_plus=96	max_min_left=49	max_min_left_minus=46	max_min_left_plus=49	max_min_right=49	max_min_right_minus=49	max_min_right_plus=49	max_pos_hash_score=198	max_right=97	max_right_minus=97	max_right_plus=88	neg_log10_pos_hash_p_value=0.8	new_junction_coverage=0.86	new_junction_read_count=85	new_junction_reference_weighted_read_count=0.00	new_junction_weighted_read_count=85.00	pos_hash_score=61	prediction=consensus	read_count_offset=9	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=0.00	side_1_gene_name=cls	side_1_gene_position=coding (298/1461 nt)	side_1_gene_product=cardiolipin synthetase	side_1_gene_strand=<	side_1_locus_tag=ECB_01223	side_1_overlap=0	side_1_possible_overlap_registers=81	side_1_possible_overlap_registers_before_trimming=90	side_1_read_count=0	side_1_redundant=0	side_1_weighted_read_count=0.00	side_2_annotate_key=repeat	side_2_continuation=0	side_2_coverage=NA	side_2_gene_name=IS1	side_2_gene_position=noncoding (768/768 nt)	side_2_gene_product=repeat region	side_2_gene_strand=>	side_2_overlap=1	side_2_possible_overlap_registers=NA	side_2_possible_overlap_registers_before_trimming=NA	side_2_read_count=NA	side_2_redundant=1	side_2_weighted_read_count=NA	total_non_overlap_reads=85
JC	65	.	REL606	2100307	-1	REL606	2122454	1	0	alignment_overlap=23	coverage_minus=56	coverage_plus=31	flanking_left=101	flanking_right=101	frequency=1.000e+00	frequency_lower=9.650e-01	frequency_upper=1.000e+00	junction_effective_depth=84.00	junction_mixture_iterations=1	junction_possible_overlap_registers=74	junction_possible_overlap_registers_before_trimming=77	key=REL606__2100307__-1__REL606__2122431__1__23____101__101__0__0	max_left=77	max_left_minus=74	max_left_plus=77	max_min_left=36	max_min_left_minus=36	max_min_left_plus=36	max_min_right=38	max_min_right_minus=38	max_min_right_plus=24	max_pos_hash_score=154	max_right=77	max_right_minus=77	max_right_plus=76	neg_log10_pos_hash_p_value=0.1	new_junction_coverage=1.12	new_junction_read_count=84	new_junction_reference_weighted_read_count=0.00	new_junction_weighted_read_count=84.00	pos_hash_score=66	prediction=consensus	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=0.00	side_1_gene_name=yegQ/ogrK	side_1_gene_position=intergenic (+171/+102)	side_1_gene_product=predicted peptidase/phage DNA-binding transcriptional regulator	side_1_gene_strand=>/<	side_1_locus_tag=ECB_01986/ECB_01987	side_1_overlap=23	side_1_possible_overlap_registers=93	side_1_possible_overlap_registers_before_trimming=100	side_1_read_count=0	side_1_redundant=0	side_1_weighted_read_count=0.00	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=0.00	side_2_gene_name=ECB_02013	side_2_gene_position=coding (172/216 nt)	side_2_gene_product=conserved hypothetical protein; putative exported protein	side_2_gene_strand=<	side_2_locus_tag=ECB_02013	side_2_overlap=0	side_2_possible_overlap_registers=73	side_2_possible_overlap_registers_before_trimming=77	side_2_read_count=0	side_2_redundant=0	side_2_weighted_read_count=0.00	total_non_overlap_reads=87
JC	66	.	REL606	2283472	1	REL606	664688	1	0	alignment_overlap=0	coverage_minus=55	coverage_plus=57	flanking_left=101	flanking_right=101	frequency=1.000e+00	frequency_lower=9.716e-01	frequency_upper=1.000e+00	junction_effective_depth=104.00	junction_mixture_iterations=1	junction_possible_overlap_registers=95	junction_possible_overlap_registers_before_trimming=100	key=REL606__2283472__1__REL606__2775877__-1__0____101__101__0__1	max_left=100	max_left_minus=99	max_left_plus=100	max_min_left=50	max_min_left_minus=48	max_min_left_plus=50	max_min_right=47	max_min_right_minus=47	max_min_right_plus=47	max_pos_hash_score=200	max_right=98	max_right_minus=98	max_right_plus=98	neg_log10_pos_hash_p_value=0.3	new_junction_coverage=1.08	new_junction_read_count=104	new_junction_reference_weighted_read_count=0.00	new_junction_weighted_read_count=104.00	pos_hash_score=75	prediction=consensus	read_count_offset=3	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=0.00	side_1_gene_name=yfaA/gyrA	side_1_gene_position=intergenic (-128/+21)	side_1_gene_product=hypothetical protein/DNA gyrase subunit A	side_1_gene_strand=</<	side_1_locus_tag=ECB_02156/ECB_02157	side_1_overlap=0	side_1_possible_overlap_registers=93	side_1_possible_overlap_registers_before_trimming=97	side_1_read_count=0	side_1_redundant=0	side_1_weighted_read_count=0.00	side_2_annotate_key=repeat	side_2_continuation=0	side_2_coverage=NA	side_2_gene_name=IS150	side_2_gene_position=noncoding (1/1443 nt)	side_2_gene_product=repeat region	side_2_gene_strand=>	side_2_overlap=0	side_2_possible_overlap_registers=NA	side_2_possible_overlap_registers_before_trimming=NA	side_2_read_count=NA	side_2_redundant=1	side_2_weighted_read_count=NA	total_non_overlap_reads=112
JC	67	.	REL606	2283474	-1	REL606	666130	-1	0	alignment_overlap=2	coverage_minus=58	coverage_plus=46	flanking_left=101	flanking_right=101	frequency=1.000e+00	frequency_lower=9.716e-01	frequency_upper=1.000e+00	junction_effective_depth=104.00	junction_mixture_iterations=1	junction_possible_overlap_registers=94	junction_possible_overlap_registers_before_trimming=98	key=REL606__2283476__-1__REL606__3894996__-1__2____101__101__0__1	max_left=95	max_left_minus=92	max_left_plus=95	max_min_left=49	max_min_left_minus=48	max_min_left_plus=49	max_min_right=49	max_min_right_minus=49	max_min_right_plus=49	max_pos_hash_score=196	max_right=95	max_right_minus=95	max_right_plus=95	neg_log10_pos_hash_p_value=0.3	new_junction_coverage=1.09	new_junction_read_count=104	new_junction_reference_weighted_read_count=0.00	new_junction_weighted_read_count=104.00	pos_hash_score=74	prediction=consensus	read_count_offset=3	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=0.00	side_1_gene_name=yfaA/gyrA	side_1_gene_position=intergenic (-130/+19)	side_1_gene_product=hypothetical protein/DNA gyrase subunit A	side_1_gene_strand=</<	side_1_locus_tag=ECB_02156/ECB_02157	side_1_overlap=0	side_1_possible_overlap_registers=91	side_1_possible_overlap_registers_before_trimming=95	side_1_read_count=0	side_1_redundant=0	side_1_weighted_read_count=0.00	side_2_annotate_key=repeat	side_2_continuation=0	side_2_coverage=NA	side_2_gene_name=IS150	side_2_gene_position=noncoding (1443/1443 nt)	side_2_gene_product=repeat region	side_2_gene_strand=>	side_2_overlap=2	side_2_possible_overlap_registers=NA	side_2_possible_overlap_registers_before_trimming=NA	side_2_read_count=NA	side_2_redundant=1	side_2_weighted_read_count=NA	total_non_overlap_reads=104
JC	68	.	REL606	2448493	1	REL606	16728	-1	0	alignment_overlap=3	coverage_minus=51	coverage_plus=43	flanking_left=101	flanking_right=101	frequency=1.000e+00	frequency_lower=9.680e-01	frequency_upper=1.000e+00	junction_effective_depth=92.00	junction_mixture_iterations=1	junction_possible_overlap_registers=95	junction_possible_overlap_registers_before_trimming=97	key=REL606__2448490__1__REL606__2774196__-1__3____101__101__0__1	max_left=97	max_left_minus=96	max_left_plus=97	max_min_left=48	max_min_left_minus=48	max_min_left_plus=47	max_min_right=49	max_min_right_minus=49	max_min_right_plus=48	max_pos_hash_score=194	max_right=93	max_right_minus=93	max_right_plus=92	neg_log10_pos_hash_p_value=0.5	new_junction_coverage=0.95	new_junction_read_count=92	new_junction_reference_weighted_read_count=0.00	new_junction_weighted_read_count=92.00	pos_hash_score=67	prediction=consensus	read_count_offset=6	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=0.00	side_1_gene_name=nupC/yfeA	side_1_gene_position=intergenic (+24/+26)	side_1_gene_product=nucleoside (except guanosine) transporter/predicted diguanylate cyclase	side_1_gene_strand=>/<	side_1_locus_tag=ECB_02302/ECB_02303	side_1_overlap=0	side_1_possible_overlap_registers=86	side_1_possible_overlap_registers_before_trimming=91	side_1_read_count=0	side_1_redundant=0	side_1_weighted_read_count=0.00	side_2_annotate_key=repeat	side_2_continuation=0	side_2_coverage=NA	side_2_gene_name=IS186	side_2_gene_position=noncoding (1343/1343 nt)	side_2_gene_product=repeat region	side_2_gene_strand=>	side_2_overlap=3	side_2_possible_overlap_registers=NA	side_2_possible_overlap_registers_before_trimming=NA	side_2_read_count=NA	side_2_redundant=1	side_2_weighted_read_count=NA	total_non_overlap_reads=94
JC	69	.	REL606	2448498	-1	REL606	15386	1	0	alignment_overlap=3	coverage_minus=36	coverage_plus=60	flanking_left=101	flanking_right=101	frequency=1.000e+00	frequency_lower=9.683e-01	frequency_upper=1.000e+00	junction_effective_depth=93.00	junction_mixture_iterations=1	junction_possible_overlap_registers=94	junction_possible_overlap_registers_before_trimming=97	key=REL606__2448501__-1__REL606__2772854__1__3____101__101__0__1	max_left=97	max_left_minus=97	max_left_plus=97	max_min_left=48	max_min_left_minus=48	max_min_left_plus=44	max_min_right=48	max_min_right_minus=44	max_min_right_plus=48	max_pos_hash_score=194	max_right=90	max_right_minus=90	max_right_plus=90	neg_log10_pos_hash_p_value=0.3	new_junction_coverage=0.97	new_junction_read_count=93	new_junction_reference_weighted_read_count=0.00	new_junction_weighted_read_count=93.00	pos_hash_score=70	prediction=consensus	read_count_offset=6	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=0.00	side_1_gene_name=nupC/yfeA	side_1_gene_position=intergenic (+29/+21)	side_1_gene_product=nucleoside (except guanosine) transporter/predicted diguanylate cyclase	side_1_gene_strand=>/<	side_1_locus_tag=ECB_02302/ECB_02303	side_1_overlap=0	side_1_possible_overlap_registers=86	side_1_possible_overlap_registers_before_trimming=91	side_1_read_count=0	side_1_redundant=0	side_1_weighted_read_count=0.00	side_2_annotate_key=repeat	side_2_continuation=0	side_2_coverage=NA	side_2_gene_name=IS186	side_2_gene_position=noncoding (1/1343 nt)	side_2_gene_product=repeat region	side_2_gene_strand=>	side_2_overlap=3	side_2_possible_overlap_registers=NA	side_2_possible_overlap_registers_before_trimming=NA	side_2_read_count=NA	side_2_redundant=1	side_2_weighted_read_count=NA	total_non_overlap_reads=96
JC	70	.	REL606	2749275	1	REL606	664688	1	0	alignment_overlap=1	coverage_minus=52	coverage_plus=63	flanking_left=101	flanking_right=101	frequency=6.108e-01	frequency_lower=5.481e-01	frequency_upper=6.708e-01	junction_effective_depth=185.00	junction_mixture_iterations=2	junction_possible_overlap_registers=97	junction_possible_overlap_registers_before_trimming=99	key=REL606__2749274__1__REL606__2775877__-1__1____101__101__0__1	max_left=99	max_left_minus=99	max_left_plus=98	max_min_left=49	max_min_left_minus=47	max_min_left_plus=49	max_min_right=50	max_min_right_minus=49	max_min_right_plus=50	max_pos_hash_score=198	max_right=98	max_right_minus=97	max_right_plus=98	neg_log10_pos_hash_p_value=0.1	new_junction_coverage=1.15	new_junction_read_count=113	new_junction_reference_weighted_read_count=0.00	new_junction_weighted_read_count=113.00	pos_hash_score=84	prediction=consensus	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=0.73	side_1_gene_name=fhlA	side_1_gene_position=coding (78/2079 nt)	side_1_gene_product=DNA-binding transcriptional activator	side_1_gene_strand=>	side_1_locus_tag=ECB_02581	side_1_overlap=0	side_1_possible_overlap_registers=97	side_1_possible_overlap_registers_before_trimming=99	side_1_read_count=72	side_1_redundant=0	side_1_weighted_read_count=72.00	side_2_annotate_key=repeat	side_2_continuation=0	side_2_coverage=NA	side_2_gene_name=IS150	side_2_gene_position=noncoding (1/1443 nt)	side_2_gene_product=repeat region	side_2_gene_strand=>	side_2_overlap=1	side_2_possible_overlap_registers=NA	side_2_possible_overlap_registers_before_trimming=NA	side_2_read_count=NA	side_2_redundant=1	side_2_weighted_read_count=NA	total_non_overlap_reads=115
CN	71	.	REL606	2100301	2122500	0	gene_name=ogrK–[ECB_02013]	gene_product=ogrK,yegZ,ECB_01989,ECB_01990,ECB_01991,ECB_01992,ECB_01993,ECB_01994,ECB_01995,ECB_01996,ECB_01997,ECB_01998,ECB_01999,ECB_02000,ECB_02001,ECB_02002,ECB_02003,ECB_02004,ECB_02005,ECB_02006,ECB_02007,ECB_02008,ECB_02009,ECB_02010,ECB_02011,ECB_02012,[ECB_02013]	locus_tag=[ECB_01987]–[ECB_02013]	relative_coverage=0.0047	tile_size=200
CN	72	.	REL606	2749301	2776100	2	gene_name=[fhlA]–[cysH]	gene_product=[fhlA],ygbA,mutS,pphB,ygbI,ygbJ,ygbK,ygbL,ygbM,ygbN,rpoS,nlpD,pcm,surE,ygbO,ispF,ispD,ftsB,ygbE,cysC,cysN,cysD,iap,insL-5,insK-3,insJ-3,[cysH]	locus_tag=[ECB_02581]–[ECB_02607]	relative_coverage=2.21	tile_size=200
CN	73	.	REL606	3895001	3900700	0	gene_name=rbsD–[rbsR]	gene_product=rbsD,rbsA,rbsC,rbsB,rbsK,[rbsR]	locus_tag=[ECB_03634]–[ECB_03639]	relative_coverage=0.0149	tile_size=200
UN	74	.	REL606	15623	15623
UN	75	.	REL606	15625	15625
UN	76	.	REL606	15628	15628
UN	77	.	REL606	15634	16462
UN	78	.	REL606	16465	16465
UN	79	.	REL606	23597	23597
UN	80	.	REL606	23599	23602
UN	81	.	REL606	23604	23604
UN	82	.	REL606	23606	23608
UN	83	.	REL606	23611	23613
UN	84	.	REL606	23615	23615
UN	85	.	REL606	23619	23622
UN	86	.	REL606	23624	23797
UN	87	.	REL606	227283	227286
UN	88	.	REL606	227289	228274
UN	89	.	REL606	228276	228277
UN	90	.	REL606	228279	228279
UN	91	.	REL606	228282	228282
UN	92	.	REL606	228284	228284
UN	93	.	REL606	228286	228286
UN	94	.	REL606	228291	228292
UN	95	.	REL606	228295	228300
UN	96	.	REL606	228303	228303
UN	97	.	REL606	228308	228313
UN	98	.	REL606	228315	228316
UN	99	.	REL606	228318	228318
UN	100	.	REL606	228323	228323
UN	101	.	REL606	229019	229019
UN	102	.	REL606	229023	231514
UN	103	.	REL606	231519	231520
UN	104	.	REL606	231528	231530
UN	105	.	REL606	231535	231537
UN	106	.	REL606	231547	231547
UN	107	.	REL606	241544	241544
UN	108	.	REL606	241546	241546
UN	109	.	REL606	241549	241688
UN	110	.	REL606	241690	241690
UN	111	.	REL606	241692	241694
UN	112	.	REL606	241696	241696
UN	113	.	REL606	241698	241698
UN	114	.	REL606	241700	241700
UN	115	.	REL606	241702	241702
UN	116	.	REL606	263096	263103
UN	117	.	REL606	263105	263105
UN	118	.	REL606	263120	263221
UN	119	.	REL606	353535	353535
UN	120	.	REL606	353537	353575
UN	121	.	REL606	353578	353578
UN	122	.	REL606	353580	353582
UN	123	.	REL606	353584	353586
UN	124	.	REL606	353588	353658
UN	125	.	REL606	353664	353667
UN	126	.	REL606	353670	353698
UN	127	.	REL606	353701	353702
UN	128	.	REL606	353716	353716
UN	129	.	REL606	377033	377033
UN	130	.	REL606	377035	377035
UN	131	.	REL606	377038	377038
UN	132	.	REL606	377050	377051
UN	133	.	REL606	377053	377056
UN	134	.	REL606	377058	377059
UN	135	.	REL606	377062	377177
UN	136	.	REL606	429787	429787
UN	137	.	REL606	429789	429894
UN	138	.	REL606	429948	429953
UN	139	.	REL606	429955	429962
UN	140	.	REL606	429964	429970
UN	141	.	REL606	429972	429972
UN	142	.	REL606	429974	430553
UN	143	.	REL606	430555	430555
UN	144	.	REL606	430557	430558
UN	145	.	REL606	430560	430560
UN	146	.	REL606	430564	430564
UN	147	.	REL606	430567	430567
UN	148	.	REL606	430569	430574
UN	149	.	REL606	430578	430578
UN	150	.	REL606	430581	430582
UN	151	.	REL606	474383	474391
UN	152	.	REL606	495699	495699
UN	153	.	REL606	495701	498075
UN	154	.	REL606	498087	498090
UN	155	.	REL606	498092	498096
UN	156	.	REL606	498098	498098
UN	157	.	REL606	498100	498100
UN	158	.	REL606	498102	498102
UN	159	.	REL606	498104	498106
UN	160	.	REL606	498108	498136
UN	161	.	REL606	498573	498573
UN	162	.	REL606	498575	498664
UN	163	.	REL606	498666	498666
UN	164	.	REL606	498668	498669
UN	165	.	REL606	498673	498681
UN	166	.	REL606	498687	498687
UN	167	.	REL606	547191	547191
UN	168	.	REL606	547219	547383
UN	169	.	REL606	553535	553536
UN	170	.	REL606	553613	553613
UN	171	.	REL606	553615	553619
UN	172	.	REL606	553621	553621
UN	173	.	REL606	553623	554444
UN	174	.	REL606	588716	588817
UN	175	.	REL606	588820	588821
UN	176	.	REL606	588823	588823
UN	177	.	REL606	588826	588827
UN	178	.	REL606	588829	588829
UN	179	.	REL606	588831	588832
UN	180	.	REL606	588834	588834
UN	181	.	REL606	588836	588836
UN	182	.	REL606	588839	588839
UN	183	.	REL606	588843	588843
UN	184	.	REL606	588845	588846
UN	185	.	REL606	588849	588849
UN	186	.	REL606	588851	588853
UN	187	.	REL606	588855	589253
UN	188	.	REL606	589255	589257
UN	189	.	REL606	589259	589260
UN	190	.	REL606	589262	589263
UN	191	.	REL606	589304	589304
UN	192	.	REL606	589306	589310
UN	193	.	REL606	589897	590016
UN	194	.	REL606	590924	590936
UN	195	.	REL606	590965	590965
UN	196	.	REL606	591020	591020
UN	197	.	REL606	591022	591023
UN	198	.	REL606	591025	591031
UN	199	.	REL606	591033	591767
UN	200	.	REL606	619392	619392
UN	201	.	REL606	619399	619403
UN	202	.	REL606	619405	619408
UN	203	.	REL606	619410	619410
UN	204	.	REL606	619412	619564
UN	205	.	REL606	634237	634387
UN	206	.	REL606	634389	634392
UN	207	.	REL606	634394	634397
UN	208	.	REL606	634399	634399
UN	209	.	REL606	634401	634406
UN	210	.	REL606	634409	634465
UN	211	.	REL606	664909	664911
UN	212	.	REL606	664914	665786
UN	213	.	REL606	665789	665791
UN	214	.	REL606	665794	665819
UN	215	.	REL606	665821	665823
UN	216	.	REL606	665827	665827
UN	217	.	REL606	1110769	1110769
UN	218	.	REL606	1110772	1110772
UN	219	.	REL606	1110774	1110776
UN	220	.	REL606	1110778	1110780
UN	221	.	REL606	1110782	1110782
UN	222	.	REL606	1110789	1111468
UN	223	.	REL606	1111472	1111472
UN	224	.	REL606	1111475	1111475
UN	225	.	REL606	1111477	1111478
UN	226	.	REL606	1111480	1111496
UN	227	.	REL606	1193811	1193811
UN	228	.	REL606	1193813	1193813
UN	229	.	REL606	1193815	1193816
UN	230	.	REL606	1193822	1193823
UN	231	.	REL606	1193826	1194115
UN	232	.	REL606	1194118	1194118
UN	233	.	REL606	1322975	1322983
UN	234	.	REL606	1322989	1322993
UN	235	.	REL606	1322997	1322997
UN	236	.	REL606	1322999	1323000
UN	237	.	REL606	1323005	1323010
UN	238	.	REL606	1323013	1323221
UN	239	.	REL606	1442872	1442872
UN	240	.	REL606	1442875	1442875
UN	241	.	REL606	1442877	1442878
UN	242	.	REL606	1442880	1442880
UN	243	.	REL606	1442888	1442889
UN	244	.	REL606	1442891	1442891
UN	245	.	REL606	1442894	1442894
UN	246	.	REL606	1442896	1442897
UN	247	.	REL606	1442899	1442903
UN	248	.	REL606	1442906	1442910
UN	249	.	REL606	1442912	1442916
UN	250	.	REL606	1442918	1442918
UN	251	.	REL606	1442920	1442922
UN	252	.	REL606	1442925	1442925
UN	253	.	REL606	1442927	1442927
UN	254	.	REL606	1442930	1442930
UN	255	.	REL606	1442932	1443651
UN	256	.	REL606	1443656	1443659
UN	257	.	REL606	1460440	1460441
UN	258	.	REL606	1460443	1460444
UN	259	.	REL606	1460446	1460645
UN	260	.	REL606	1500483	1500483
UN	261	.	REL606	1500485	1500488
UN	262	.	REL606	1500490	1500490
UN	263	.	REL606	1500493	1500499
UN	264	.	REL606	1500505	1500505
UN	265	.	REL606	1500508	1502876
UN	266	.	REL606	1502878	1502879
UN	267	.	REL606	1502884	1502887
UN	268	.	REL606	1502889	1502889
UN	269	.	REL606	1502891	1502891
UN	270	.	REL606	1502893	1502894
UN	271	.	REL606	1502900	1502900
UN	272	.	REL606	1502903	1502903
UN	273	.	REL606	1502905	1502907
UN	274	.	REL606	1502909	1502909
UN	275	.	REL606	1502919	1502919
UN	276	.	REL606	1502922	1502922
UN	277	.	REL606	1502925	1502925
UN	278	.	REL606	1502927	1502927
UN	279	.	REL606	1502929	1502929
UN	280	.	REL606	1502936	1502939
UN	281	.	REL606	1502941	1502941
UN	282	.	REL606	1502944	1502957
UN	283	.	REL606	1503407	1503410
UN	284	.	REL606	1503412	1503412
UN	285	.	REL606	1503415	1503415
UN	286	.	REL606	1503426	1503430
UN	287	.	REL606	1503433	1503435
UN	288	.	REL606	1503438	1503439
UN	289	.	REL606	1503462	1503463
UN	290	.	REL606	1503466	1503467
UN	291	.	REL606	1503469	1503470
UN	292	.	REL606	1503474	1503479
UN	293	.	REL606	1503483	1503484
UN	294	.	REL606	1503541	1503549
UN	295	.	REL606	1606842	1607647
UN	296	.	REL606	1607649	1607649
UN	297	.	REL606	1607651	1607653
UN	298	.	REL606	1607684	1607684
UN	299	.	REL606	1607688	1607688
UN	300	.	REL606	1607691	1607691
UN	301	.	REL606	1608190	1608191
UN	302	.	REL606	1608194	1608194
UN	303	.	REL606	1608197	1608199
UN	304	.	REL606	1608201	1608203
UN	305	.	REL606	1608212	1608212
UN	306	.	REL606	1608227	1608229
UN	307	.	REL606	1608232	1608892
UN	308	.	REL606	1608894	1608894
UN	309	.	REL606	1608903	1608903
UN	310	.	REL606	1615737	1615737
UN	311	.	REL606	1615740	1615760
UN	312	.	REL606	1615762	1615763
UN	313	.	REL606	1615765	1616457
UN	314	.	REL606	1616459	1616470
UN	315	.	REL606	1616472	1616472
UN	316	.	REL606	1616476	1616476
UN	317	.	REL606	1857382	1858140
UN	318	.	REL606	1973522	1973522
UN	319	.	REL606	1973530	1973531
UN	320	.	REL606	1973533	1973533
UN	321	.	REL606	1973565	1973744
UN	322	.	REL606	1973748	1973748
UN	323	.	REL606	1973750	1973751
UN	324	.	REL606	1973755	1973755
UN	325	.	REL606	1973757	1973759
UN	326	.	REL606	1973764	1973764
UN	327	.	REL606	2034623	2034777
UN	328	.	REL606	2100308	2122436
UN	329	.	REL606	2128878	2129110
UN	330	.	REL606	2143313	2143476
UN	331	.	REL606	2143478	2143478
UN	332	.	REL606	2143480	2143480
UN	333	.	REL606	2143482	2143484
UN	334	.	REL606	2143488	2143488
UN	335	.	REL606	2143490	2143539
UN	336	.	REL606	2143542	2143542
UN	337	.	REL606	2143544	2143545
UN	338	.	REL606	2143547	2143547
UN	339	.	REL606	2254540	2254549
UN	340	.	REL606	2254820	2254820
UN	341	.	REL606	2254822	2254990
UN	342	.	REL606	2254993	2254993
UN	343	.	REL606	2263011	2263011
UN	344	.	REL606	2263015	2263015
UN	345	.	REL606	2263020	2263215
UN	346	.	REL606	2263218	2263218
UN	347	.	REL606	2263220	2263221
UN	348	.	REL606	2407737	2407738
UN	349	.	REL606	2407743	2407743
UN	350	.	REL606	2407748	2407752
UN	351	.	REL606	2407755	2407909
UN	352	.	REL606	2618357	2618360
UN	353	.	REL606	2618362	2618363
UN	354	.	REL606	2618368	2618372
UN	355	.	REL606	2618376	2618509
UN	356	.	REL606	2647761	2647762
UN	357	.	REL606	2647764	2647764
UN	358	.	REL606	2647766	2647766
UN	359	.	REL606	2647768	2647769
UN	360	.	REL606	2647771	2647773
UN	361	.	REL606	2647776	2647779
UN	362	.	REL606	2647783	2647783
UN	363	.	REL606	2647785	2652451
UN	364	.	REL606	2683017	2683022
UN	365	.	REL606	2773152	2773152
UN	366	.	REL606	2773157	2773878
UN	367	.	REL606	2773881	2773881
UN	368	.	REL606	2773884	2773885
UN	369	.	REL606	2773888	2773890
UN	370	.	REL606	2773894	2773895
UN	371	.	REL606	2773897	2773899
UN	372	.	REL606	2773902	2773903
UN	373	.	REL606	2773906	2773906
UN	374	.	REL606	2773908	2773909
UN	375	.	REL606	2773911	2773911
UN	376	.	REL606	2774779	2775578
UN	377	.	REL606	2775580	2775580
UN	378	.	REL606	2775582	2775586
UN	379	.	REL606	2822509	2822511
UN	380	.	REL606	2822513	2822513
UN	381	.	REL606	2822515	2822634
UN	382	.	REL606	2822637	2822639
UN	383	.	REL606	2822658	2822661
UN	384	.	REL606	2882974	2882975
UN	385	.	REL606	2882979	2882979
UN	386	.	REL606	2882981	2882983
UN	387	.	REL606	2882986	2882986
UN	388	.	REL606	2882989	2883641
UN	389	.	REL606	2883643	2883645
UN	390	.	REL606	2883653	2883653
UN	391	.	REL606	2883657	2883659
UN	392	.	REL606	2883665	2883665
UN	393	.	REL606	3024251	3024252
UN	394	.	REL606	3024254	3024255
UN	395	.	REL606	3024259	3024420
UN	396	.	REL606	3024423	3024423
UN	397	.	REL606	3024425	3024425
UN	398	.	REL606	3351910	3351911
UN	399	.	REL606	3351923	3351944
UN	400	.	REL606	3351946	3351946
UN	401	.	REL606	3351948	3354021
UN	402	.	REL606	3354023	3354026
UN	403	.	REL606	3354030	3354030
UN	404	.	REL606	3354032	3354052
UN	405	.	REL606	3354680	3354694
UN	406	.	REL606	3354696	3354698
UN	407	.	REL606	3354701	3354701
UN	408	.	REL606	3354705	3354705
UN	409	.	REL606	3354707	3354714
UN	410	.	REL606	3354716	3354716
UN	411	.	REL606	3354718	3354718
UN	412	.	REL606	3354748	3354750
UN	413	.	REL606	3354756	3354757
UN	414	.	REL606	3354759	3354759
UN	415	.	REL606	3354761	3354761
UN	416	.	REL606	3354763	3354765
UN	417	.	REL606	3354767	3354772
UN	418	.	REL606	3354774	3354774
UN	419	.	REL606	3354776	3354778
UN	420	.	REL606	3354780	3354781
UN	421	.	REL606	3354784	3354784
UN	422	.	REL606	3354786	3356401
UN	423	.	REL606	3356405	3356406
UN	424	.	REL606	3433750	3433944
UN	425	.	REL606	3433950	3433951
UN	426	.	REL606	3433957	3433958
UN	427	.	REL606	3550261	3551066
UN	428	.	REL606	3551072	3551072
UN	429	.	REL606	3551075	3551075
UN	430	.	REL606	3551636	3551636
UN	431	.	REL606	3551638	3551638
UN	432	.	REL606	3551656	3551656
UN	433	.	REL606	3551660	3551661
UN	434	.	REL606	3551663	3551666
UN	435	.	REL606	3551668	3551699
UN	436	.	REL606	3551705	3551705
UN	437	.	REL606	3551708	3552007
UN	438	.	REL606	3552015	3552015
UN	439	.	REL606	3552017	3552017
UN	440	.	REL606	3552031	3552031
UN	441	.	REL606	3651352	3651352
UN	442	.	REL606	3651354	3651356
UN	443	.	REL606	3651359	3652226
UN	444	.	REL606	3652230	3652230
UN	445	.	REL606	3652232	3652234
UN	446	.	REL606	3697313	3697313
UN	447	.	REL606	3697316	3698046
UN	448	.	REL606	3698048	3698078
UN	449	.	REL606	3698080	3698084
UN	450	.	REL606	3698601	3698933
UN	451	.	REL606	3698938	3698938
UN	452	.	REL606	3698942	3698942
UN	453	.	REL606	3698945	3698948
UN	454	.	REL606	3698952	3698952
UN	455	.	REL606	3698954	3698955
UN	456	.	REL606	3698957	3699036
UN	457	.	REL606	3699038	3699038
UN	458	.	REL606	3699041	3699042
UN	459	.	REL606	3699045	3699045
UN	460	.	REL606	3699049	3699049
UN	461	.	REL606	3699053	3699053
UN	462	.	REL606	3741469	3741469
UN	463	.	REL606	3741504	3741650
UN	464	.	REL606	3741652	3741652
UN	465	.	REL606	3741654	3741659
UN	466	.	REL606	3893806	3893809
UN	467	.	REL606	3893812	3893812
UN	468	.	REL606	3893821	3893821
UN	469	.	REL606	3893853	3894618
UN	470	.	REL606	3894620	3894620
UN	471	.	REL606	3894999	3900623
UN	472	.	REL606	3903743	3903745
UN	473	.	REL606	3903747	3903747
UN	474	.	REL606	3903751	3908359
UN	475	.	REL606	4013911	4015254
UN	476	.	REL606	4015257	4015263
UN	477	.	REL606	4015271	4015271
UN	478	.	REL606	4015275	4015275
UN	479	.	REL606	4015880	4015881
UN	480	.	REL606	4015883	4015899
UN	481	.	REL606	4015901	4015901
UN	482	.	REL606	4015903	4015905
UN	483	.	REL606	4016328	4016328
UN	484	.	REL606	4016330	4016330
UN	485	.	REL606	4016332	4016332
UN	486	.	REL606	4016334	4017105
UN	487	.	REL606	4017107	4017107
UN	488	.	REL606	4017117	4017117
UN	489	.	REL606	4017121	4017122
UN	490	.	REL606	4018395	4018395
UN	491	.	REL606	4018397	4018793
UN	492	.	REL606	4018803	4018803
UN	493	.	REL606	4146527	4146527
UN	494	.	REL606	4146529	4146529
UN	495	.	REL606	4146531	4146533
UN	496	.	REL606	4146535	4146536
UN	497	.	REL606	4146538	4146539
UN	498	.	REL606	4146543	4146543
UN	499	.	REL606	4146545	4146545
UN	500	.	REL606	4146547	4146549
UN	501	.	REL606	4146551	4146552
UN	502	.	REL606	4146556	4147794
UN	503	.	REL606	4147799	4147799
UN	504	.	REL606	4147802	4147802
UN	505	.	REL606	4147804	4147808
UN	506	.	REL606	4147810	4147810
UN	507	.	REL606	4147814	4147814
UN	508	.	REL606	4147817	4147818
UN	509	.	REL606	4147820	4147820
UN	510	.	REL606	4147822	4147822
UN	511	.	REL606	4147834	4147834
UN	512	.	REL606	4147837	4147838
UN	513	.	REL606	4147840	4147840
UN	514	.	REL606	4149048	4149050
UN	515	.	REL606	4149074	4149075
UN	516	.	REL606	4149078	4149078
UN	517	.	REL606	4149082	4149889
UN	518	.	REL606	4149891	4149891
UN	519	.	REL606	4149893	4149893
UN	520	.	REL606	4149901	4149901
UN	521	.	REL606	4150385	4150385
UN	522	.	REL606	4150388	4150390
UN	523	.	REL606	4150395	4151223
UN	524	.	REL606	4151226	4151226
UN	525	.	REL606	4151228	4151229
UN	526	.	REL606	4151232	4151233
UN	527	.	REL606	4151236	4151242
UN	528	.	REL606	4151245	4151246
UN	529	.	REL606	4151248	4151250
UN	530	.	REL606	4151264	4151264
UN	531	.	REL606	4187767	4187767
UN	532	.	REL606	4187769	4192570
UN	533	.	REL606	4192572	4192574
UN	534	.	REL606	4192576	4192576
UN	535	.	REL606	4192578	4192578
UN	536	.	REL606	4274870	4274870
UN	537	.	REL606	4505633	4505764
UN	538	.	REL606	4505776	4505776
UN	539	.	REL606	4505857	4505857
UN	540	.	REL606	4521795	4521795
UN	541	.	REL606	4521819	4521824
UN	542	.	REL606	4521830	4521842
UN	543	.	REL606	4521846	4521846
UN	544	.	REL606	4521848	4521849
UN	545	.	REL606	4521851	4521853
UN	546	.	REL606	4521855	4522018
UN	547	.	REL606	4522021	4522022
UN	548	.	REL606	4550992	4551148
