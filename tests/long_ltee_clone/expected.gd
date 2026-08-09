#=GENOME_DIFF	1.0
#=CREATED	16:39:16 08 Aug 2026
#=PROGRAM	breseq 0.50.0 revision 791a79d27ad5
#=COMMAND	./src/breseq/breseq -j 8 -o ./tests/long_ltee_clone -r ./tests/long_ltee_clone/../data/downloads/ncbi_REL606/GCF_000017985.1_ASM1798v1_genomic.gbff.gz ./tests/long_ltee_clone/../data/downloads/ena_SRR2589061/SRR2589061_1.fastq.gz ./tests/long_ltee_clone/../data/downloads/ena_SRR2589061/SRR2589061_2.fastq.gz
#=REFSEQ	./tests/long_ltee_clone/../data/downloads/ncbi_REL606/GCF_000017985.1_ASM1798v1_genomic.gbff.gz
#=READSEQ	./tests/long_ltee_clone/../data/downloads/ena_SRR2589061/SRR2589061_1.fastq.gz
#=READSEQ	./tests/long_ltee_clone/../data/downloads/ena_SRR2589061/SRR2589061_2.fastq.gz
#=CONVERTED-BASES	481720712
#=CONVERTED-READS	4769512
#=INPUT-BASES	481722328
#=INPUT-READS	4769528
#=MAPPED-BASES	473713701
#=MAPPED-READS	4690619
DEL	1	32,38	NC_012967	474383	9	gene_name=ybaL	gene_position=coding (915-923/1677 nt)	gene_product=YbaL family putative K(+) efflux transporter	gene_strand=<	genes_overlapping=ybaL	locus_tag=ECB_RS02245	locus_tags_overlapping=ECB_RS02245	mutation_category=small_indel	position_end=474391	position_start=474383	ref_seq=GTCGCCAGC
SNP	2	17	NC_012967	651601	A	aa_new_seq=L	aa_position=7	aa_ref_seq=Q	codon_new_seq=CTG	codon_number=7	codon_position=2	codon_ref_seq=CAG	gene_name=rsfS	gene_position=20	gene_product=ribosome silencing factor	gene_strand=<	genes_overlapping=rsfS	locus_tag=ECB_RS03155	locus_tags_overlapping=ECB_RS03155	mutation_category=snp_nonsynonymous	position_end=651601	position_start=651601	ref_seq=T	snp_type=nonsynonymous	transl_table=11
SNP	3	18	NC_012967	1089367	A	aa_new_seq=V	aa_position=291	aa_ref_seq=A	codon_new_seq=GTC	codon_number=291	codon_position=2	codon_ref_seq=GCC	gene_name=rutA	gene_position=872	gene_product=pyrimidine utilization protein A	gene_strand=<	genes_overlapping=rutA	locus_tag=ECB_RS05365	locus_tags_overlapping=ECB_RS05365	mutation_category=snp_nonsynonymous	position_end=1089367	position_start=1089367	ref_seq=G	snp_type=nonsynonymous	transl_table=11
SNP	4	19	NC_012967	1166950	T	aa_new_seq=I	aa_position=148	aa_ref_seq=T	codon_new_seq=ATT	codon_number=148	codon_position=2	codon_ref_seq=ACT	gene_name=fabF	gene_position=443	gene_product=beta-ketoacyl-ACP synthase II	gene_strand=>	genes_overlapping=fabF	locus_tag=ECB_RS05790	locus_tags_overlapping=ECB_RS05790	mutation_category=snp_nonsynonymous	position_end=1166950	position_start=1166950	ref_seq=C	snp_type=nonsynonymous	transl_table=11
SNP	5	20	NC_012967	1330018	A	aa_new_seq=E	aa_position=200	aa_ref_seq=A	codon_new_seq=GAG	codon_number=200	codon_position=2	codon_ref_seq=GCG	gene_name=topA	gene_position=599	gene_product=type I DNA topoisomerase	gene_strand=>	genes_overlapping=topA	locus_tag=ECB_RS06635	locus_tags_overlapping=ECB_RS06635	mutation_category=snp_nonsynonymous	position_end=1330018	position_start=1330018	ref_seq=C	snp_type=nonsynonymous	transl_table=11
SNP	6	21	NC_012967	1587635	G	aa_new_seq=C	aa_position=99	aa_ref_seq=Y	codon_new_seq=TGC	codon_number=99	codon_position=2	codon_ref_seq=TAC	gene_name=marA	gene_position=296	gene_product=MDR efflux pump AcrAB transcriptional activator MarA	gene_strand=>	genes_overlapping=marA	locus_tag=ECB_RS07885	locus_tags_overlapping=ECB_RS07885	mutation_category=snp_nonsynonymous	position_end=1587635	position_start=1587635	ref_seq=A	snp_type=nonsynonymous	transl_table=11
SNP	7	22	NC_012967	1733865	T	aa_new_seq=S	aa_position=301	aa_ref_seq=A	codon_new_seq=TCC	codon_number=301	codon_position=1	codon_ref_seq=GCC	gene_name=pykF	gene_position=901	gene_product=pyruvate kinase PykF	gene_strand=>	genes_overlapping=pykF	locus_tag=ECB_RS08685	locus_tags_overlapping=ECB_RS08685	mutation_category=snp_nonsynonymous	position_end=1733865	position_start=1733865	ref_seq=G	snp_type=nonsynonymous	transl_table=11
DEL	8	33,49	NC_012967	2100308	22146	gene_name=ogrK–ECB_RS10640	gene_product=ogrK,ECB_RS25980,ECB_RS25985,ECB_RS10535,ECB_RS23815,ECB_RS23820,ECB_RS10540,ECB_RS10545,ECB_RS10550,ECB_RS10555,ECB_RS10560,ECB_RS10565,ECB_RS10570,ECB_RS25495,ECB_RS10580,ECB_RS10585,ECB_RS10590,ECB_RS25990,ECB_RS10595,ECB_RS10600,ECB_RS10605,ECB_RS10610,ECB_RS10615,ECB_RS10620,ECB_RS10625,ECB_RS25745,ECB_RS10630,ECB_RS10635,ECB_RS10640	genes_inactivated=ogrK,ECB_RS25980,ECB_RS25985,ECB_RS10535,ECB_RS23815,ECB_RS23820,ECB_RS10540,ECB_RS10545,ECB_RS10550,ECB_RS10555,ECB_RS10560,ECB_RS10565,ECB_RS10570,ECB_RS25495,ECB_RS10580,ECB_RS10585,ECB_RS10590,ECB_RS25990,ECB_RS10595,ECB_RS10600,ECB_RS10605,ECB_RS10610,ECB_RS10615,ECB_RS10620,ECB_RS10625,ECB_RS25745,ECB_RS10630,ECB_RS10635,ECB_RS10640	locus_tag=[ECB_RS23805]–[ECB_RS10640]	locus_tags_inactivated=ECB_RS23805,ECB_RS25980,ECB_RS25985,ECB_RS10535,ECB_RS23815,ECB_RS23820,ECB_RS10540,ECB_RS10545,ECB_RS10550,ECB_RS10555,ECB_RS10560,ECB_RS10565,ECB_RS10570,ECB_RS25495,ECB_RS10580,ECB_RS10585,ECB_RS10590,ECB_RS25990,ECB_RS10595,ECB_RS10600,ECB_RS10605,ECB_RS10610,ECB_RS10615,ECB_RS10620,ECB_RS10625,ECB_RS25745,ECB_RS10630,ECB_RS10635,ECB_RS10640	mutation_category=large_deletion	position_end=2122453	position_start=2100308	ref_seq=22146-bp
SNP	9	23	NC_012967	2275076	G	aa_new_seq=H	aa_position=474	aa_ref_seq=Q	codon_new_seq=CAC	codon_number=474	codon_position=3	codon_ref_seq=CAA	gene_name=yfaQ	gene_position=1422	gene_product=YfaQ family protein	gene_strand=<	genes_overlapping=yfaQ	locus_tag=ECB_RS11345	locus_tags_overlapping=ECB_RS11345	mutation_category=snp_nonsynonymous	position_end=2275076	position_start=2275076	ref_seq=T	snp_type=nonsynonymous	transl_table=11
SNP	10	25	NC_012967	3107610	A	aa_new_seq=I	aa_position=99	aa_ref_seq=I	codon_new_seq=ATA	codon_number=99	codon_position=3	codon_ref_seq=ATT	gene_name=ygiN	gene_position=297	gene_product=putative quinol monooxygenase	gene_strand=>	genes_overlapping=ygiN	locus_tag=ECB_RS15330	locus_tags_overlapping=ECB_RS15330	mutation_category=snp_synonymous	position_end=3107610	position_start=3107610	ref_seq=T	snp_type=synonymous	transl_table=11
SNP	11	26	NC_012967	3249544	T	aa_new_seq=I	aa_position=569	aa_ref_seq=V	codon_new_seq=ATT	codon_number=569	codon_position=1	codon_ref_seq=GTT	gene_name=infB	gene_position=1705	gene_product=translation initiation factor IF-2	gene_strand=<	genes_overlapping=infB	locus_tag=ECB_RS16035	locus_tags_overlapping=ECB_RS16035	mutation_category=snp_nonsynonymous	position_end=3249544	position_start=3249544	ref_seq=C	snp_type=nonsynonymous	transl_table=11
SNP	12	27	NC_012967	3289548	A	gene_name=yhcC/gltB	gene_position=intergenic (-263/-412)	gene_product=TIGR01212 family radical SAM protein/glutamate synthase large subunit	gene_strand=</>	locus_tag=ECB_RS16250/ECB_RS16255	mutation_category=snp_intergenic	position_end=3289548	position_start=3289548	ref_seq=G	snp_type=intergenic
SNP	13	28	NC_012967	3340296	A	aa_new_seq=I	aa_position=251	aa_ref_seq=M	codon_new_seq=ATA	codon_number=251	codon_position=3	codon_ref_seq=ATG	gene_name=yhdJ	gene_position=753	gene_product=adenine-specific DNA-methyltransferase	gene_strand=>	genes_overlapping=yhdJ	locus_tag=ECB_RS16490	locus_tags_overlapping=ECB_RS16490	mutation_category=snp_nonsynonymous	position_end=3340296	position_start=3340296	ref_seq=G	snp_type=nonsynonymous	transl_table=11
SNP	14	29	NC_012967	4100211	T	aa_new_seq=L	aa_position=340	aa_ref_seq=F	codon_new_seq=TTA	codon_number=340	codon_position=3	codon_ref_seq=TTC	gene_name=hslU	gene_position=1020	gene_product=HslU--HslV peptidase ATPase subunit	gene_strand=<	genes_overlapping=hslU	locus_tag=ECB_RS20235	locus_tags_overlapping=ECB_RS20235	mutation_category=snp_nonsynonymous	position_end=4100211	position_start=4100211	ref_seq=G	snp_type=nonsynonymous	transl_table=11
SNP	15	30	NC_012967	4201937	T	aa_new_seq=D	aa_position=208	aa_ref_seq=G	codon_new_seq=GAT	codon_number=208	codon_position=2	codon_ref_seq=GGT	gene_name=iclR	gene_position=623	gene_product=glyoxylate bypass operon transcriptional repressor IclR	gene_strand=<	genes_overlapping=iclR	locus_tag=ECB_RS20685	locus_tags_overlapping=ECB_RS20685	mutation_category=snp_nonsynonymous	position_end=4201937	position_start=4201937	ref_seq=C	snp_type=nonsynonymous	transl_table=11
SNP	16	31	NC_012967	4270684	T	aa_new_seq=L	aa_position=100	aa_ref_seq=L	codon_new_seq=CTT	codon_number=100	codon_position=3	codon_ref_seq=CTG	gene_name=nrfE	gene_position=300	gene_product=heme lyase CcmF/NrfE family subunit	gene_strand=>	genes_overlapping=nrfE	locus_tag=ECB_RS20990	locus_tags_overlapping=ECB_RS20990	mutation_category=snp_synonymous	position_end=4270684	position_start=4270684	ref_seq=G	snp_type=synonymous	transl_table=11
RA	17	.	NC_012967	651601	0	T	A	aa_new_seq=L	aa_position=7	aa_ref_seq=Q	allele_frequencies=A:1.000e+00	codon_new_seq=CTG	codon_number=7	codon_position=2	codon_ref_seq=CAG	fisher_strand_p_value=1.00000e+00	frequency=1.000e+00	frequency_lower=9.852e-01	frequency_upper=1.000e+00	gene_name=rsfS	gene_position=20	gene_product=ribosome silencing factor	gene_strand=<	ks_quality_p_value=1.00000e+00	locus_tag=ECB_RS03155	major_base=A	major_cov=52/39	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=52/39	new_seq=A	prediction=consensus	ref_cov=0/0	ref_seq=T	score=329.9	snp_type=nonsynonymous	total_cov=52/39	transl_table=11
RA	18	.	NC_012967	1089367	0	G	A	aa_new_seq=V	aa_position=291	aa_ref_seq=A	allele_frequencies=A:1.000e+00	codon_new_seq=GTC	codon_number=291	codon_position=2	codon_ref_seq=GCC	fisher_strand_p_value=1.00000e+00	frequency=1.000e+00	frequency_lower=9.742e-01	frequency_upper=1.000e+00	gene_name=rutA	gene_position=872	gene_product=pyrimidine utilization protein A	gene_strand=<	ks_quality_p_value=1.00000e+00	locus_tag=ECB_RS05365	major_base=A	major_cov=35/36	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=35/36	new_seq=A	prediction=consensus	ref_cov=0/0	ref_seq=G	score=250.8	snp_type=nonsynonymous	total_cov=36/36	transl_table=11
RA	19	.	NC_012967	1166950	0	C	T	aa_new_seq=I	aa_position=148	aa_ref_seq=T	allele_frequencies=T:1.000e+00	codon_new_seq=ATT	codon_number=148	codon_position=2	codon_ref_seq=ACT	fisher_strand_p_value=1.00000e+00	frequency=1.000e+00	frequency_lower=9.806e-01	frequency_upper=1.000e+00	gene_name=fabF	gene_position=443	gene_product=beta-ketoacyl-ACP synthase II	gene_strand=>	ks_quality_p_value=1.00000e+00	locus_tag=ECB_RS05790	major_base=T	major_cov=35/34	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=35/34	new_seq=T	prediction=consensus	ref_cov=0/0	ref_seq=C	score=247.1	snp_type=nonsynonymous	total_cov=35/34	transl_table=11
RA	20	.	NC_012967	1330018	0	C	A	aa_new_seq=E	aa_position=200	aa_ref_seq=A	allele_frequencies=A:1.000e+00	codon_new_seq=GAG	codon_number=200	codon_position=2	codon_ref_seq=GCG	fisher_strand_p_value=1.00000e+00	frequency=1.000e+00	frequency_lower=9.800e-01	frequency_upper=1.000e+00	gene_name=topA	gene_position=599	gene_product=type I DNA topoisomerase	gene_strand=>	ks_quality_p_value=1.00000e+00	locus_tag=ECB_RS06635	major_base=A	major_cov=29/38	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=29/38	new_seq=A	prediction=consensus	ref_cov=0/0	ref_seq=C	score=235.1	snp_type=nonsynonymous	total_cov=29/38	transl_table=11
RA	21	.	NC_012967	1587635	0	A	G	aa_new_seq=C	aa_position=99	aa_ref_seq=Y	allele_frequencies=G:1.000e+00	codon_new_seq=TGC	codon_number=99	codon_position=2	codon_ref_seq=TAC	fisher_strand_p_value=1.00000e+00	frequency=1.000e+00	frequency_lower=9.828e-01	frequency_upper=1.000e+00	gene_name=marA	gene_position=296	gene_product=MDR efflux pump AcrAB transcriptional activator MarA	gene_strand=>	ks_quality_p_value=1.00000e+00	locus_tag=ECB_RS07885	major_base=G	major_cov=38/40	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=38/40	new_seq=G	prediction=consensus	ref_cov=0/0	ref_seq=A	score=321.2	snp_type=nonsynonymous	total_cov=38/40	transl_table=11
RA	22	.	NC_012967	1733865	0	G	T	aa_new_seq=S	aa_position=301	aa_ref_seq=A	allele_frequencies=G:2.353e-02,T:9.765e-01	codon_new_seq=TCC	codon_number=301	codon_position=1	codon_ref_seq=GCC	fisher_strand_p_value=1.00000e+00	frequency=9.765e-01	frequency_lower=9.388e-01	frequency_upper=9.944e-01	gene_name=pykF	gene_position=901	gene_product=pyruvate kinase PykF	gene_strand=>	ks_quality_p_value=7.25210e-01	locus_tag=ECB_RS08685	major_base=T	major_cov=36/47	major_frequency=9.765e-01	minor_base=G	minor_cov=1/1	new_cov=36/47	new_seq=T	prediction=consensus	ref_cov=1/1	ref_seq=G	score=291.7	snp_type=nonsynonymous	total_cov=37/48	transl_table=11
RA	23	.	NC_012967	2275076	0	T	G	aa_new_seq=H	aa_position=474	aa_ref_seq=Q	allele_frequencies=G:1.000e+00	codon_new_seq=CAC	codon_number=474	codon_position=3	codon_ref_seq=CAA	fisher_strand_p_value=1.00000e+00	frequency=1.000e+00	frequency_lower=9.777e-01	frequency_upper=1.000e+00	gene_name=yfaQ	gene_position=1422	gene_product=YfaQ family protein	gene_strand=<	ks_quality_p_value=1.00000e+00	locus_tag=ECB_RS11345	major_base=G	major_cov=29/31	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=29/31	new_seq=G	prediction=consensus	ref_cov=0/0	ref_seq=T	score=246.3	snp_type=nonsynonymous	total_cov=29/31	transl_table=11
RA	24	.	NC_012967	2766603	0	G	A	aa_new_seq=G	aa_position=13	aa_ref_seq=G	allele_frequencies=A:4.399e-01,G:5.601e-01	codon_new_seq=GGT	codon_number=13	codon_position=3	codon_ref_seq=GGC	consensus_reject=FREQUENCY_CUTOFF	fisher_strand_p_value=8.86816e-01	frequency=4.399e-01	frequency_lower=3.829e-01	frequency_upper=4.980e-01	gene_name=ispF	gene_position=39	gene_product=2-C-methyl-D-erythritol 2,4-cyclodiphosphate synthase	gene_strand=<	ks_quality_p_value=9.60450e-01	locus_tag=ECB_RS13720	major_base=G	major_cov=54/58	major_frequency=5.601e-01	minor_base=A	minor_cov=44/44	new_cov=44/44	new_seq=A	prediction=polymorphism	ref_cov=54/58	ref_seq=G	score=282.6	snp_type=synonymous	total_cov=98/102	transl_table=11
RA	25	.	NC_012967	3107610	0	T	A	aa_new_seq=I	aa_position=99	aa_ref_seq=I	allele_frequencies=A:1.000e+00	codon_new_seq=ATA	codon_number=99	codon_position=3	codon_ref_seq=ATT	fisher_strand_p_value=1.00000e+00	frequency=1.000e+00	frequency_lower=9.847e-01	frequency_upper=1.000e+00	gene_name=ygiN	gene_position=297	gene_product=putative quinol monooxygenase	gene_strand=>	ks_quality_p_value=1.00000e+00	locus_tag=ECB_RS15330	major_base=A	major_cov=42/46	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=42/46	new_seq=A	prediction=consensus	ref_cov=0/0	ref_seq=T	score=313.3	snp_type=synonymous	total_cov=42/46	transl_table=11
RA	26	.	NC_012967	3249544	0	C	T	aa_new_seq=I	aa_position=569	aa_ref_seq=V	allele_frequencies=T:1.000e+00	codon_new_seq=ATT	codon_number=569	codon_position=1	codon_ref_seq=GTT	fisher_strand_p_value=1.00000e+00	frequency=1.000e+00	frequency_lower=9.849e-01	frequency_upper=1.000e+00	gene_name=infB	gene_position=1705	gene_product=translation initiation factor IF-2	gene_strand=<	ks_quality_p_value=1.00000e+00	locus_tag=ECB_RS16035	major_base=T	major_cov=40/49	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=40/49	new_seq=T	prediction=consensus	ref_cov=0/0	ref_seq=C	score=305.7	snp_type=nonsynonymous	total_cov=40/49	transl_table=11
RA	27	.	NC_012967	3289548	0	G	A	allele_frequencies=A:1.000e+00	fisher_strand_p_value=1.00000e+00	frequency=1.000e+00	frequency_lower=9.896e-01	frequency_upper=1.000e+00	gene_name=yhcC/gltB	gene_position=intergenic (-263/-412)	gene_product=TIGR01212 family radical SAM protein/glutamate synthase large subunit	gene_strand=</>	ks_quality_p_value=1.00000e+00	locus_tag=ECB_RS16250/ECB_RS16255	major_base=A	major_cov=56/73	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=56/73	prediction=consensus	ref_cov=0/0	score=463.9	snp_type=intergenic	total_cov=56/73
RA	28	.	NC_012967	3340296	0	G	A	aa_new_seq=I	aa_position=251	aa_ref_seq=M	allele_frequencies=A:1.000e+00	codon_new_seq=ATA	codon_number=251	codon_position=3	codon_ref_seq=ATG	fisher_strand_p_value=1.00000e+00	frequency=1.000e+00	frequency_lower=9.846e-01	frequency_upper=1.000e+00	gene_name=yhdJ	gene_position=753	gene_product=adenine-specific DNA-methyltransferase	gene_strand=>	ks_quality_p_value=1.00000e+00	locus_tag=ECB_RS16490	major_base=A	major_cov=38/49	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=38/49	new_seq=A	prediction=consensus	ref_cov=0/0	ref_seq=G	score=313.9	snp_type=nonsynonymous	total_cov=38/49	transl_table=11
RA	29	.	NC_012967	4100211	0	G	T	aa_new_seq=L	aa_position=340	aa_ref_seq=F	allele_frequencies=T:1.000e+00	codon_new_seq=TTA	codon_number=340	codon_position=3	codon_ref_seq=TTC	fisher_strand_p_value=1.00000e+00	frequency=1.000e+00	frequency_lower=9.880e-01	frequency_upper=1.000e+00	gene_name=hslU	gene_position=1020	gene_product=HslU--HslV peptidase ATPase subunit	gene_strand=<	ks_quality_p_value=1.00000e+00	locus_tag=ECB_RS20235	major_base=T	major_cov=50/62	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=50/62	new_seq=T	prediction=consensus	ref_cov=0/0	ref_seq=G	score=409.7	snp_type=nonsynonymous	total_cov=50/62	transl_table=11
RA	30	.	NC_012967	4201937	0	C	T	aa_new_seq=D	aa_position=208	aa_ref_seq=G	allele_frequencies=T:1.000e+00	codon_new_seq=GAT	codon_number=208	codon_position=2	codon_ref_seq=GGT	fisher_strand_p_value=1.00000e+00	frequency=1.000e+00	frequency_lower=9.877e-01	frequency_upper=1.000e+00	gene_name=iclR	gene_position=623	gene_product=glyoxylate bypass operon transcriptional repressor IclR	gene_strand=<	ks_quality_p_value=1.00000e+00	locus_tag=ECB_RS20685	major_base=T	major_cov=55/54	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=55/54	new_seq=T	prediction=consensus	ref_cov=0/0	ref_seq=C	score=377.1	snp_type=nonsynonymous	total_cov=55/54	transl_table=11
RA	31	.	NC_012967	4270684	0	G	T	aa_new_seq=L	aa_position=100	aa_ref_seq=L	allele_frequencies=T:1.000e+00	codon_new_seq=CTT	codon_number=100	codon_position=3	codon_ref_seq=CTG	fisher_strand_p_value=1.00000e+00	frequency=1.000e+00	frequency_lower=9.811e-01	frequency_upper=1.000e+00	gene_name=nrfE	gene_position=300	gene_product=heme lyase CcmF/NrfE family subunit	gene_strand=>	ks_quality_p_value=1.00000e+00	locus_tag=ECB_RS20990	major_base=T	major_cov=34/37	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=34/37	new_seq=T	prediction=consensus	ref_cov=0/0	ref_seq=G	score=251.3	snp_type=synonymous	total_cov=34/37	transl_table=11
MC	32	.	NC_012967	474383	474391	0	0	gene_name=ybaL	gene_position=coding (915-923/1677 nt)	gene_product=YbaL family putative K(+) efflux transporter	gene_strand=<	left_inside_cov=0	left_outside_cov=99	locus_tag=ECB_RS02245	right_inside_cov=0	right_outside_cov=101
MC	33	.	NC_012967	2100308	2122453	0	0	gene_name=ogrK–ECB_RS10640	gene_product=ogrK,ECB_RS25980,ECB_RS25985,ECB_RS10535,ECB_RS23815,ECB_RS23820,ECB_RS10540,ECB_RS10545,ECB_RS10550,ECB_RS10555,ECB_RS10560,ECB_RS10565,ECB_RS10570,ECB_RS25495,ECB_RS10580,ECB_RS10585,ECB_RS10590,ECB_RS25990,ECB_RS10595,ECB_RS10600,ECB_RS10605,ECB_RS10610,ECB_RS10615,ECB_RS10620,ECB_RS10625,ECB_RS25745,ECB_RS10630,ECB_RS10635,ECB_RS10640	left_inside_cov=0	left_outside_cov=90	locus_tag=[ECB_RS23805]–[ECB_RS10640]	right_inside_cov=23	right_outside_cov=111
MC	34	.	NC_012967	3894931	3900622	65	0	gene_name=[ECB_RS25215]–[rbsR]	gene_product=[ECB_RS25215],rbsD,rbsA,rbsC,rbsB,rbsK,[rbsR]	left_inside_cov=42	left_outside_cov=43	locus_tag=[ECB_RS25215]–[ECB_RS19230]	right_inside_cov=1	right_outside_cov=134
JC	35	.	NC_012967	1	1	NC_012967	4629812	-1	0	alignment_overlap=0	coverage_minus=51	coverage_plus=74	flanking_left=101	flanking_right=101	frequency=NA	frequency_lower=NA	frequency_upper=NA	ignore=CIRCULAR_CHROMOSOME	junction_effective_depth=116.00	junction_possible_overlap_registers=98	junction_possible_overlap_registers_before_trimming=100	key=NC_012967__1__1__NC_012967__4629812__-1__0____101__101__0__0	max_left=100	max_left_minus=100	max_left_plus=100	max_min_left=50	max_min_left_minus=50	max_min_left_plus=50	max_min_right=49	max_min_right_minus=48	max_min_right_plus=49	max_pos_hash_score=200	max_right=100	max_right_minus=100	max_right_plus=100	neg_log10_pos_hash_p_value=0.1	new_junction_coverage=1.17	new_junction_read_count=116	new_junction_reference_weighted_read_count=0.01	new_junction_weighted_read_count=115.99	pos_hash_score=88	prediction=unknown	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=NA	side_1_gene_name=–/thrL	side_1_gene_position=intergenic (–/-189)	side_1_gene_product=–/thr operon leader peptide	side_1_gene_strand=–/>	side_1_locus_tag=–/ECB_RS00005	side_1_overlap=0	side_1_possible_overlap_registers=0	side_1_possible_overlap_registers_before_trimming=100	side_1_read_count=0	side_1_redundant=0	side_1_weighted_read_count=0.00	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=NA	side_2_gene_name=yjtD/–	side_2_gene_position=intergenic (+24/–)	side_2_gene_product=tRNA/rRNA methyltransferase/–	side_2_gene_strand=>/–	side_2_locus_tag=ECB_RS22810/–	side_2_overlap=0	side_2_possible_overlap_registers=0	side_2_possible_overlap_registers_before_trimming=100	side_2_read_count=0	side_2_redundant=0	side_2_weighted_read_count=0.00	total_non_overlap_reads=125
JC	36	.	NC_012967	467506	1	NC_012967	2775876	-1	0	alignment_overlap=1	coverage_minus=54	coverage_plus=65	flanking_left=101	flanking_right=101	frequency=9.754e-01	frequency_lower=9.377e-01	frequency_upper=9.933e-01	junction_effective_depth=122.00	junction_mixture_iterations=1	junction_possible_overlap_registers=97	junction_possible_overlap_registers_before_trimming=99	key=NC_012967__467506__1__NC_012967__2775877__-1__1____101__101__0__1	max_left=98	max_left_minus=95	max_left_plus=98	max_min_left=49	max_min_left_minus=45	max_min_left_plus=49	max_min_right=50	max_min_right_minus=50	max_min_right_plus=47	max_pos_hash_score=198	max_right=97	max_right_minus=97	max_right_plus=97	neg_log10_pos_hash_p_value=0.1	new_junction_coverage=1.21	new_junction_read_count=119	new_junction_reference_weighted_read_count=0.00	new_junction_weighted_read_count=119.00	pos_hash_score=82	prediction=consensus	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=0.03	side_1_gene_name=htpG	side_1_gene_position=coding (320/1875 nt)	side_1_gene_product=molecular chaperone HtpG	side_1_gene_strand=>	side_1_locus_tag=ECB_RS02220	side_1_overlap=1	side_1_possible_overlap_registers=97	side_1_possible_overlap_registers_before_trimming=100	side_1_read_count=3	side_1_redundant=0	side_1_weighted_read_count=3.00	side_2_annotate_key=repeat	side_2_continuation=0	side_2_coverage=NA	side_2_gene_name=ECB_RS25125/cysH	side_2_gene_position=intergenic (-46/+203)	side_2_gene_product=IS3-like element IS150 family transposase/phosphoadenosine phosphosulfate reductase	side_2_gene_strand=</<	side_2_locus_tag=ECB_RS25125/ECB_RS13780	side_2_overlap=0	side_2_possible_overlap_registers=NA	side_2_possible_overlap_registers_before_trimming=NA	side_2_read_count=NA	side_2_redundant=1	side_2_weighted_read_count=NA	total_non_overlap_reads=119
JC	37	.	NC_012967	467509	-1	NC_012967	590471	-1	0	alignment_overlap=0	coverage_minus=67	coverage_plus=70	flanking_left=101	flanking_right=101	frequency=9.925e-01	frequency_lower=9.652e-01	frequency_upper=9.996e-01	junction_effective_depth=135.00	junction_mixture_iterations=1	junction_possible_overlap_registers=98	junction_possible_overlap_registers_before_trimming=100	key=NC_012967__467509__-1__NC_012967__590471__-1__0____101__101__0__1	max_left=100	max_left_minus=99	max_left_plus=100	max_min_left=49	max_min_left_minus=48	max_min_left_plus=49	max_min_right=50	max_min_right_minus=50	max_min_right_plus=50	max_pos_hash_score=200	max_right=98	max_right_minus=98	max_right_plus=97	neg_log10_pos_hash_p_value=0.0	new_junction_coverage=1.35	new_junction_read_count=134	new_junction_reference_weighted_read_count=0.00	new_junction_weighted_read_count=134.00	pos_hash_score=92	prediction=consensus	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=0.01	side_1_gene_name=htpG	side_1_gene_position=coding (323/1875 nt)	side_1_gene_product=molecular chaperone HtpG	side_1_gene_strand=>	side_1_locus_tag=ECB_RS02220	side_1_overlap=0	side_1_possible_overlap_registers=97	side_1_possible_overlap_registers_before_trimming=100	side_1_read_count=1	side_1_redundant=0	side_1_weighted_read_count=1.00	side_2_annotate_key=repeat	side_2_continuation=0	side_2_coverage=NA	side_2_gene_name=ECB_RS24845/hokE	side_2_gene_position=intergenic (+26/-72)	side_2_gene_product=IS3 family transposase/type I toxin-antitoxin system toxin HokE	side_2_gene_strand=>/>	side_2_locus_tag=ECB_RS24845/ECB_RS02860	side_2_overlap=0	side_2_possible_overlap_registers=NA	side_2_possible_overlap_registers_before_trimming=NA	side_2_read_count=NA	side_2_redundant=1	side_2_weighted_read_count=NA	total_non_overlap_reads=137
JC	38	.	NC_012967	474382	-1	NC_012967	474392	1	0	alignment_overlap=8	coverage_minus=49	coverage_plus=48	flanking_left=101	flanking_right=101	frequency=1.000e+00	frequency_lower=9.676e-01	frequency_upper=1.000e+00	junction_effective_depth=91.00	junction_mixture_iterations=1	junction_possible_overlap_registers=86	junction_possible_overlap_registers_before_trimming=92	key=NC_012967__474382__-1__NC_012967__474384__1__8____101__101__0__0	max_left=91	max_left_minus=86	max_left_plus=91	max_min_left=46	max_min_left_minus=46	max_min_left_plus=43	max_min_right=46	max_min_right_minus=42	max_min_right_plus=46	max_pos_hash_score=184	max_right=91	max_right_minus=90	max_right_plus=91	neg_log10_pos_hash_p_value=0.2	new_junction_coverage=1.04	new_junction_read_count=91	new_junction_reference_weighted_read_count=0.00	new_junction_weighted_read_count=91.00	pos_hash_score=72	prediction=consensus	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=0.00	side_1_gene_name=ybaL	side_1_gene_position=coding (924/1677 nt)	side_1_gene_product=YbaL family putative K(+) efflux transporter	side_1_gene_strand=<	side_1_locus_tag=ECB_RS02245	side_1_overlap=8	side_1_possible_overlap_registers=98	side_1_possible_overlap_registers_before_trimming=100	side_1_read_count=0	side_1_redundant=0	side_1_weighted_read_count=0.00	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=0.00	side_2_gene_name=ybaL	side_2_gene_position=coding (914/1677 nt)	side_2_gene_product=YbaL family putative K(+) efflux transporter	side_2_gene_strand=<	side_2_locus_tag=ECB_RS02245	side_2_overlap=0	side_2_possible_overlap_registers=90	side_2_possible_overlap_registers_before_trimming=92	side_2_read_count=0	side_2_redundant=0	side_2_weighted_read_count=0.00	total_non_overlap_reads=97
JC	39	.	NC_012967	589555	1	NC_012967	591113	-1	0	alignment_overlap=1	coverage_minus=10	coverage_plus=16	flanking_left=101	flanking_right=101	frequency=NA	frequency_lower=NA	frequency_upper=NA	junction_effective_depth=26.00	junction_possible_overlap_registers=97	junction_possible_overlap_registers_before_trimming=99	key=NC_012967__589555__1__NC_012967__591114__-1__1____101__101__1__1	max_left=95	max_left_minus=93	max_left_plus=95	max_min_left=43	max_min_left_minus=31	max_min_left_plus=43	max_min_right=45	max_min_right_minus=35	max_min_right_plus=45	max_pos_hash_score=198	max_right=97	max_right_minus=92	max_right_plus=97	neg_log10_pos_hash_p_value=3.8	new_junction_coverage=0.26	new_junction_read_count=26	new_junction_reference_weighted_read_count=0.00	new_junction_weighted_read_count=26.00	pos_hash_score=24	prediction=unknown	reject=COVERAGE_EVENNESS_SKEW	side_1_annotate_key=repeat	side_1_continuation=0	side_1_coverage=NA	side_1_gene_name=ECB_RS24840/ECB_RS23205	side_1_gene_position=intergenic (+1/+15)	side_1_gene_product=IS3-like element IS150 family transposase/IS1-like element IS1A family transposase	side_1_gene_strand=>/<	side_1_locus_tag=ECB_RS24840/ECB_RS23205	side_1_overlap=1	side_1_possible_overlap_registers=NA	side_1_possible_overlap_registers_before_trimming=NA	side_1_read_count=NA	side_1_redundant=1	side_1_weighted_read_count=NA	side_2_annotate_key=repeat	side_2_continuation=0	side_2_coverage=NA	side_2_gene_name=ECB_RS02865	side_2_gene_position=coding (342/1113 nt)	side_2_gene_product=IS4-like element IS421 family transposase	side_2_gene_strand=>	side_2_locus_tag=ECB_RS02865	side_2_overlap=0	side_2_possible_overlap_registers=NA	side_2_possible_overlap_registers_before_trimming=NA	side_2_read_count=NA	side_2_redundant=1	side_2_weighted_read_count=NA	total_non_overlap_reads=26
JC	40	.	NC_012967	590322	-1	NC_012967	591106	1	0	alignment_overlap=1	coverage_minus=10	coverage_plus=4	flanking_left=101	flanking_right=101	frequency=NA	frequency_lower=NA	frequency_upper=NA	junction_effective_depth=14.00	junction_possible_overlap_registers=97	junction_possible_overlap_registers_before_trimming=99	key=NC_012967__590322__-1__NC_012967__591105__1__1____101__101__1__1	max_left=86	max_left_minus=86	max_left_plus=61	max_min_left=46	max_min_left_minus=46	max_min_left_plus=45	max_min_right=44	max_min_right_minus=44	max_min_right_plus=36	max_pos_hash_score=198	max_right=80	max_right_minus=69	max_right_plus=80	neg_log10_pos_hash_p_value=4.8	new_junction_coverage=0.14	new_junction_read_count=14	new_junction_reference_weighted_read_count=0.00	new_junction_weighted_read_count=14.00	pos_hash_score=14	prediction=unknown	reject=COVERAGE_EVENNESS_SKEW	side_1_annotate_key=repeat	side_1_continuation=0	side_1_coverage=NA	side_1_gene_name=ECB_RS23205/ECB_RS24845	side_1_gene_position=intergenic (-55/-1)	side_1_gene_product=IS1-like element IS1A family transposase/IS3 family transposase	side_1_gene_strand=</>	side_1_locus_tag=ECB_RS23205/ECB_RS24845	side_1_overlap=1	side_1_possible_overlap_registers=NA	side_1_possible_overlap_registers_before_trimming=NA	side_1_read_count=NA	side_1_redundant=1	side_1_weighted_read_count=NA	side_2_annotate_key=repeat	side_2_continuation=0	side_2_coverage=NA	side_2_gene_name=ECB_RS02865	side_2_gene_position=coding (335/1113 nt)	side_2_gene_product=IS4-like element IS421 family transposase	side_2_gene_strand=>	side_2_locus_tag=ECB_RS02865	side_2_overlap=0	side_2_possible_overlap_registers=NA	side_2_possible_overlap_registers_before_trimming=NA	side_2_read_count=NA	side_2_redundant=1	side_2_weighted_read_count=NA	total_non_overlap_reads=14
JC	41	.	NC_012967	664688	1	NC_012967	1462251	1	0	alignment_overlap=0	coverage_minus=46	coverage_plus=49	flanking_left=101	flanking_right=101	frequency=1.000e+00	frequency_lower=9.686e-01	frequency_upper=1.000e+00	junction_effective_depth=94.00	junction_mixture_iterations=1	junction_possible_overlap_registers=98	junction_possible_overlap_registers_before_trimming=100	key=NC_012967__664688__1__NC_012967__1462251__1__0____101__101__1__0	max_left=98	max_left_minus=98	max_left_plus=98	max_min_left=49	max_min_left_minus=49	max_min_left_plus=49	max_min_right=50	max_min_right_minus=50	max_min_right_plus=49	max_pos_hash_score=200	max_right=100	max_right_minus=97	max_right_plus=100	neg_log10_pos_hash_p_value=0.5	new_junction_coverage=0.95	new_junction_read_count=94	new_junction_reference_weighted_read_count=0.00	new_junction_weighted_read_count=94.00	pos_hash_score=67	prediction=consensus	side_1_annotate_key=repeat	side_1_continuation=0	side_1_coverage=NA	side_1_gene_name=gltL/ECB_RS24860	side_1_gene_position=intergenic (-1/-47)	side_1_gene_product=glutamate/aspartate ABC transporter ATP-binding protein GltL/IS3-like element IS150 family transposase	side_1_gene_strand=</>	side_1_locus_tag=ECB_RS03215/ECB_RS24860	side_1_overlap=0	side_1_possible_overlap_registers=NA	side_1_possible_overlap_registers_before_trimming=NA	side_1_read_count=NA	side_1_redundant=1	side_1_weighted_read_count=NA	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=0.00	side_2_gene_name=hokB/trg	side_2_gene_position=intergenic (-71/-328)	side_2_gene_product=type I toxin-antitoxin system toxin HokB/methyl-accepting chemotaxis protein Trg	side_2_gene_strand=</>	side_2_locus_tag=ECB_RS07305/ECB_RS07310	side_2_overlap=0	side_2_possible_overlap_registers=98	side_2_possible_overlap_registers_before_trimming=100	side_2_read_count=0	side_2_redundant=0	side_2_weighted_read_count=0.00	total_non_overlap_reads=95
JC	42	.	NC_012967	664688	1	NC_012967	2157883	-1	0	alignment_overlap=0	coverage_minus=47	coverage_plus=35	flanking_left=101	flanking_right=101	frequency=9.873e-01	frequency_lower=9.422e-01	frequency_upper=9.993e-01	junction_effective_depth=81.00	junction_mixture_iterations=1	junction_possible_overlap_registers=98	junction_possible_overlap_registers_before_trimming=100	key=NC_012967__664688__1__NC_012967__2157883__-1__0____101__101__1__0	max_left=95	max_left_minus=95	max_left_plus=91	max_min_left=48	max_min_left_minus=48	max_min_left_plus=47	max_min_right=49	max_min_right_minus=49	max_min_right_plus=48	max_pos_hash_score=200	max_right=100	max_right_minus=100	max_right_plus=100	neg_log10_pos_hash_p_value=0.7	new_junction_coverage=0.80	new_junction_read_count=80	new_junction_reference_weighted_read_count=0.00	new_junction_weighted_read_count=80.00	pos_hash_score=63	prediction=consensus	side_1_annotate_key=repeat	side_1_continuation=0	side_1_coverage=NA	side_1_gene_name=gltL/ECB_RS24860	side_1_gene_position=intergenic (-1/-47)	side_1_gene_product=glutamate/aspartate ABC transporter ATP-binding protein GltL/IS3-like element IS150 family transposase	side_1_gene_strand=</>	side_1_locus_tag=ECB_RS03215/ECB_RS24860	side_1_overlap=0	side_1_possible_overlap_registers=NA	side_1_possible_overlap_registers_before_trimming=NA	side_1_read_count=NA	side_1_redundant=1	side_1_weighted_read_count=NA	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=0.01	side_2_gene_name=yehM	side_2_gene_position=coding (993/2280 nt)	side_2_gene_product=DUF5682 family protein	side_2_gene_strand=>	side_2_locus_tag=ECB_RS10800	side_2_overlap=0	side_2_possible_overlap_registers=95	side_2_possible_overlap_registers_before_trimming=100	side_2_read_count=1	side_2_redundant=0	side_2_weighted_read_count=1.00	total_non_overlap_reads=82
JC	43	.	NC_012967	666130	-1	NC_012967	968328	-1	0	alignment_overlap=0	coverage_minus=55	coverage_plus=55	flanking_left=101	flanking_right=101	frequency=9.820e-01	frequency_lower=9.444e-01	frequency_upper=9.968e-01	junction_effective_depth=111.00	junction_mixture_iterations=1	junction_possible_overlap_registers=97	junction_possible_overlap_registers_before_trimming=100	key=NC_012967__666130__-1__NC_012967__968328__-1__0____101__101__1__0	max_left=98	max_left_minus=98	max_left_plus=98	max_min_left=49	max_min_left_minus=49	max_min_left_plus=49	max_min_right=50	max_min_right_minus=49	max_min_right_plus=50	max_pos_hash_score=200	max_right=99	max_right_minus=99	max_right_plus=99	neg_log10_pos_hash_p_value=0.1	new_junction_coverage=1.11	new_junction_read_count=109	new_junction_reference_weighted_read_count=0.00	new_junction_weighted_read_count=109.00	pos_hash_score=81	prediction=consensus	side_1_annotate_key=repeat	side_1_continuation=0	side_1_coverage=NA	side_1_gene_name=ECB_RS24860/gltL	side_1_gene_position=intergenic (+26/+1)	side_1_gene_product=IS3-like element IS150 family transposase/glutamate/aspartate ABC transporter ATP-binding protein GltL	side_1_gene_strand=>/<	side_1_locus_tag=ECB_RS24860/ECB_RS03230	side_1_overlap=0	side_1_possible_overlap_registers=NA	side_1_possible_overlap_registers_before_trimming=NA	side_1_read_count=NA	side_1_redundant=1	side_1_weighted_read_count=NA	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=0.02	side_2_gene_name=pflA/pflB	side_2_gene_position=intergenic (-155/+37)	side_2_gene_product=pyruvate formate lyase 1-activating protein/formate C-acetyltransferase	side_2_gene_strand=</<	side_2_locus_tag=ECB_RS04790/ECB_RS04795	side_2_overlap=0	side_2_possible_overlap_registers=97	side_2_possible_overlap_registers_before_trimming=100	side_2_read_count=2	side_2_redundant=0	side_2_weighted_read_count=2.00	total_non_overlap_reads=110
JC	44	.	NC_012967	666130	-1	NC_012967	1462253	-1	0	alignment_overlap=0	coverage_minus=51	coverage_plus=60	flanking_left=101	flanking_right=101	frequency=1.000e+00	frequency_lower=9.729e-01	frequency_upper=1.000e+00	junction_effective_depth=109.00	junction_mixture_iterations=1	junction_possible_overlap_registers=97	junction_possible_overlap_registers_before_trimming=100	key=NC_012967__666130__-1__NC_012967__1462253__-1__0____101__101__1__0	max_left=97	max_left_minus=97	max_left_plus=97	max_min_left=48	max_min_left_minus=48	max_min_left_plus=43	max_min_right=50	max_min_right_minus=50	max_min_right_plus=50	max_pos_hash_score=200	max_right=100	max_right_minus=100	max_right_plus=100	neg_log10_pos_hash_p_value=0.1	new_junction_coverage=1.11	new_junction_read_count=109	new_junction_reference_weighted_read_count=0.00	new_junction_weighted_read_count=109.00	pos_hash_score=83	prediction=consensus	side_1_annotate_key=repeat	side_1_continuation=0	side_1_coverage=NA	side_1_gene_name=ECB_RS24860/gltL	side_1_gene_position=intergenic (+26/+1)	side_1_gene_product=IS3-like element IS150 family transposase/glutamate/aspartate ABC transporter ATP-binding protein GltL	side_1_gene_strand=>/<	side_1_locus_tag=ECB_RS24860/ECB_RS03230	side_1_overlap=0	side_1_possible_overlap_registers=NA	side_1_possible_overlap_registers_before_trimming=NA	side_1_read_count=NA	side_1_redundant=1	side_1_weighted_read_count=NA	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=0.00	side_2_gene_name=hokB/trg	side_2_gene_position=intergenic (-73/-326)	side_2_gene_product=type I toxin-antitoxin system toxin HokB/methyl-accepting chemotaxis protein Trg	side_2_gene_strand=</>	side_2_locus_tag=ECB_RS07305/ECB_RS07310	side_2_overlap=0	side_2_possible_overlap_registers=97	side_2_possible_overlap_registers_before_trimming=100	side_2_read_count=0	side_2_redundant=0	side_2_weighted_read_count=0.00	total_non_overlap_reads=111
JC	45	.	NC_012967	666130	-1	NC_012967	2157881	1	0	alignment_overlap=0	coverage_minus=61	coverage_plus=61	flanking_left=101	flanking_right=101	frequency=1.000e+00	frequency_lower=9.753e-01	frequency_upper=1.000e+00	junction_effective_depth=120.00	junction_mixture_iterations=1	junction_possible_overlap_registers=97	junction_possible_overlap_registers_before_trimming=100	key=NC_012967__666130__-1__NC_012967__2157881__1__0____101__101__1__0	max_left=97	max_left_minus=97	max_left_plus=97	max_min_left=50	max_min_left_minus=50	max_min_left_plus=50	max_min_right=50	max_min_right_minus=50	max_min_right_plus=48	max_pos_hash_score=200	max_right=100	max_right_minus=100	max_right_plus=98	neg_log10_pos_hash_p_value=0.1	new_junction_coverage=1.22	new_junction_read_count=120	new_junction_reference_weighted_read_count=0.00	new_junction_weighted_read_count=120.00	pos_hash_score=82	prediction=consensus	side_1_annotate_key=repeat	side_1_continuation=0	side_1_coverage=NA	side_1_gene_name=ECB_RS24860/gltL	side_1_gene_position=intergenic (+26/+1)	side_1_gene_product=IS3-like element IS150 family transposase/glutamate/aspartate ABC transporter ATP-binding protein GltL	side_1_gene_strand=>/<	side_1_locus_tag=ECB_RS24860/ECB_RS03230	side_1_overlap=0	side_1_possible_overlap_registers=NA	side_1_possible_overlap_registers_before_trimming=NA	side_1_read_count=NA	side_1_redundant=1	side_1_weighted_read_count=NA	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=0.00	side_2_gene_name=yehM	side_2_gene_position=coding (991/2280 nt)	side_2_gene_product=DUF5682 family protein	side_2_gene_strand=>	side_2_locus_tag=ECB_RS10800	side_2_overlap=0	side_2_possible_overlap_registers=97	side_2_possible_overlap_registers_before_trimming=100	side_2_read_count=0	side_2_redundant=0	side_2_weighted_read_count=0.00	total_non_overlap_reads=122
JC	46	.	NC_012967	968326	1	NC_012967	2775877	-1	0	alignment_overlap=0	coverage_minus=60	coverage_plus=64	flanking_left=101	flanking_right=101	frequency=1.000e+00	frequency_lower=9.759e-01	frequency_upper=1.000e+00	junction_effective_depth=123.00	junction_mixture_iterations=1	junction_possible_overlap_registers=98	junction_possible_overlap_registers_before_trimming=100	key=NC_012967__968326__1__NC_012967__2775877__-1__0____101__101__0__1	max_left=100	max_left_minus=100	max_left_plus=99	max_min_left=50	max_min_left_minus=50	max_min_left_plus=42	max_min_right=50	max_min_right_minus=50	max_min_right_plus=49	max_pos_hash_score=200	max_right=99	max_right_minus=99	max_right_plus=99	neg_log10_pos_hash_p_value=0.1	new_junction_coverage=1.24	new_junction_read_count=123	new_junction_reference_weighted_read_count=0.00	new_junction_weighted_read_count=123.00	pos_hash_score=88	prediction=consensus	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=0.00	side_1_gene_name=pflA/pflB	side_1_gene_position=intergenic (-153/+39)	side_1_gene_product=pyruvate formate lyase 1-activating protein/formate C-acetyltransferase	side_1_gene_strand=</<	side_1_locus_tag=ECB_RS04790/ECB_RS04795	side_1_overlap=0	side_1_possible_overlap_registers=96	side_1_possible_overlap_registers_before_trimming=100	side_1_read_count=0	side_1_redundant=0	side_1_weighted_read_count=0.00	side_2_annotate_key=repeat	side_2_continuation=0	side_2_coverage=NA	side_2_gene_name=ECB_RS25125/cysH	side_2_gene_position=intergenic (-47/+202)	side_2_gene_product=IS3-like element IS150 family transposase/phosphoadenosine phosphosulfate reductase	side_2_gene_strand=</<	side_2_locus_tag=ECB_RS25125/ECB_RS13780	side_2_overlap=0	side_2_possible_overlap_registers=NA	side_2_possible_overlap_registers_before_trimming=NA	side_2_read_count=NA	side_2_redundant=1	side_2_weighted_read_count=NA	total_non_overlap_reads=124
JC	47	.	NC_012967	1305934	1	NC_012967	1322721	1	0	alignment_overlap=0	coverage_minus=44	coverage_plus=41	flanking_left=101	flanking_right=101	frequency=9.105e-01	frequency_lower=8.455e-01	frequency_upper=9.542e-01	junction_effective_depth=92.00	junction_mixture_iterations=1	junction_possible_overlap_registers=97	junction_possible_overlap_registers_before_trimming=100	key=NC_012967__1305934__1__NC_012967__1322721__1__0____101__101__0__1	max_left=99	max_left_minus=93	max_left_plus=99	max_min_left=50	max_min_left_minus=50	max_min_left_plus=45	max_min_right=50	max_min_right_minus=49	max_min_right_plus=50	max_pos_hash_score=200	max_right=97	max_right_minus=96	max_right_plus=97	neg_log10_pos_hash_p_value=0.5	new_junction_coverage=0.85	new_junction_read_count=84	new_junction_reference_weighted_read_count=0.00	new_junction_weighted_read_count=84.00	pos_hash_score=69	prediction=consensus	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=0.08	side_1_gene_name=cls	side_1_gene_position=coding (306/1461 nt)	side_1_gene_product=cardiolipin synthase	side_1_gene_strand=<	side_1_locus_tag=ECB_RS06490	side_1_overlap=0	side_1_possible_overlap_registers=94	side_1_possible_overlap_registers_before_trimming=100	side_1_read_count=8	side_1_redundant=0	side_1_weighted_read_count=8.00	side_2_annotate_key=repeat	side_2_continuation=0	side_2_coverage=NA	side_2_gene_name=ECB_RS06585/ECB_RS06590	side_2_gene_position=intergenic (+1/-55)	side_2_gene_product=DUF2207 domain-containing protein/IS1-like element IS1A family transposase	side_2_gene_strand=>/>	side_2_locus_tag=ECB_RS06585/ECB_RS06590	side_2_overlap=0	side_2_possible_overlap_registers=NA	side_2_possible_overlap_registers_before_trimming=NA	side_2_read_count=NA	side_2_redundant=1	side_2_weighted_read_count=NA	total_non_overlap_reads=85
JC	48	.	NC_012967	1305943	-1	NC_012967	1323487	-1	0	alignment_overlap=1	coverage_minus=44	coverage_plus=41	flanking_left=101	flanking_right=101	frequency=9.328e-01	frequency_lower=8.724e-01	frequency_upper=9.700e-01	junction_effective_depth=91.00	junction_mixture_iterations=1	junction_possible_overlap_registers=97	junction_possible_overlap_registers_before_trimming=99	key=NC_012967__1305943__-1__NC_012967__1323488__-1__1____101__101__0__1	max_left=98	max_left_minus=98	max_left_plus=96	max_min_left=49	max_min_left_minus=46	max_min_left_plus=49	max_min_right=49	max_min_right_minus=49	max_min_right_plus=49	max_pos_hash_score=198	max_right=97	max_right_minus=97	max_right_plus=88	neg_log10_pos_hash_p_value=0.8	new_junction_coverage=0.86	new_junction_read_count=85	new_junction_reference_weighted_read_count=0.00	new_junction_weighted_read_count=85.00	pos_hash_score=61	prediction=consensus	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=0.06	side_1_gene_name=cls	side_1_gene_position=coding (297/1461 nt)	side_1_gene_product=cardiolipin synthase	side_1_gene_strand=<	side_1_locus_tag=ECB_RS06490	side_1_overlap=1	side_1_possible_overlap_registers=95	side_1_possible_overlap_registers_before_trimming=100	side_1_read_count=6	side_1_redundant=0	side_1_weighted_read_count=6.00	side_2_annotate_key=repeat	side_2_continuation=0	side_2_coverage=NA	side_2_gene_name=ECB_RS06590/ECB_RS06595	side_2_gene_position=intergenic (+14/-1)	side_2_gene_product=IS1-like element IS1A family transposase/DUF2207 domain-containing protein	side_2_gene_strand=>/>	side_2_locus_tag=ECB_RS06590/ECB_RS06595	side_2_overlap=0	side_2_possible_overlap_registers=NA	side_2_possible_overlap_registers_before_trimming=NA	side_2_read_count=NA	side_2_redundant=1	side_2_weighted_read_count=NA	total_non_overlap_reads=85
JC	49	.	NC_012967	2100307	-1	NC_012967	2122454	1	0	alignment_overlap=23	coverage_minus=56	coverage_plus=31	flanking_left=101	flanking_right=101	frequency=1.000e+00	frequency_lower=9.650e-01	frequency_upper=1.000e+00	junction_effective_depth=84.00	junction_mixture_iterations=1	junction_possible_overlap_registers=74	junction_possible_overlap_registers_before_trimming=77	key=NC_012967__2100307__-1__NC_012967__2122431__1__23____101__101__0__0	max_left=77	max_left_minus=74	max_left_plus=77	max_min_left=36	max_min_left_minus=36	max_min_left_plus=36	max_min_right=38	max_min_right_minus=38	max_min_right_plus=24	max_pos_hash_score=154	max_right=77	max_right_minus=77	max_right_plus=76	neg_log10_pos_hash_p_value=0.1	new_junction_coverage=1.12	new_junction_read_count=84	new_junction_reference_weighted_read_count=0.00	new_junction_weighted_read_count=84.00	pos_hash_score=66	prediction=consensus	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=0.00	side_1_gene_name=yegQ/ogrK	side_1_gene_position=intergenic (+171/+102)	side_1_gene_product=tRNA 5-hydroxyuridine modification protein YegQ/prophage transcriptional regulator OgrK	side_1_gene_strand=>/<	side_1_locus_tag=ECB_RS10515/ECB_RS23805	side_1_overlap=23	side_1_possible_overlap_registers=93	side_1_possible_overlap_registers_before_trimming=100	side_1_read_count=0	side_1_redundant=0	side_1_weighted_read_count=0.00	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=0.00	side_2_gene_name=ECB_RS10640/yegR	side_2_gene_position=intergenic (+107/+159)	side_2_gene_product=site-specific integrase/protein YegR	side_2_gene_strand=>/<	side_2_locus_tag=ECB_RS10640/ECB_RS10650	side_2_overlap=0	side_2_possible_overlap_registers=73	side_2_possible_overlap_registers_before_trimming=77	side_2_read_count=0	side_2_redundant=0	side_2_weighted_read_count=0.00	total_non_overlap_reads=87
JC	50	.	NC_012967	2283472	1	NC_012967	2775877	-1	0	alignment_overlap=0	coverage_minus=55	coverage_plus=57	flanking_left=101	flanking_right=101	frequency=1.000e+00	frequency_lower=9.716e-01	frequency_upper=1.000e+00	junction_effective_depth=104.00	junction_mixture_iterations=1	junction_possible_overlap_registers=95	junction_possible_overlap_registers_before_trimming=100	key=NC_012967__2283472__1__NC_012967__2775877__-1__0____101__101__0__1	max_left=100	max_left_minus=99	max_left_plus=100	max_min_left=50	max_min_left_minus=48	max_min_left_plus=50	max_min_right=47	max_min_right_minus=47	max_min_right_plus=47	max_pos_hash_score=200	max_right=98	max_right_minus=98	max_right_plus=98	neg_log10_pos_hash_p_value=0.3	new_junction_coverage=1.08	new_junction_read_count=104	new_junction_reference_weighted_read_count=0.00	new_junction_weighted_read_count=104.00	pos_hash_score=75	prediction=consensus	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=0.00	side_1_gene_name=yfaA/gyrA	side_1_gene_position=intergenic (-128/+21)	side_1_gene_product=DUF2138 domain-containing protein/DNA topoisomerase (ATP-hydrolyzing) subunit A	side_1_gene_strand=</<	side_1_locus_tag=ECB_RS11360/ECB_RS11365	side_1_overlap=0	side_1_possible_overlap_registers=93	side_1_possible_overlap_registers_before_trimming=100	side_1_read_count=0	side_1_redundant=0	side_1_weighted_read_count=0.00	side_2_annotate_key=repeat	side_2_continuation=0	side_2_coverage=NA	side_2_gene_name=ECB_RS25125/cysH	side_2_gene_position=intergenic (-47/+202)	side_2_gene_product=IS3-like element IS150 family transposase/phosphoadenosine phosphosulfate reductase	side_2_gene_strand=</<	side_2_locus_tag=ECB_RS25125/ECB_RS13780	side_2_overlap=0	side_2_possible_overlap_registers=NA	side_2_possible_overlap_registers_before_trimming=NA	side_2_read_count=NA	side_2_redundant=1	side_2_weighted_read_count=NA	total_non_overlap_reads=112
JC	51	.	NC_012967	2283476	-1	NC_012967	3894994	-1	0	alignment_overlap=2	coverage_minus=58	coverage_plus=46	flanking_left=101	flanking_right=101	frequency=9.644e-01	frequency_lower=9.192e-01	frequency_upper=9.881e-01	junction_effective_depth=108.00	junction_mixture_iterations=1	junction_possible_overlap_registers=94	junction_possible_overlap_registers_before_trimming=98	key=NC_012967__2283476__-1__NC_012967__3894996__-1__2____101__101__0__1	max_left=95	max_left_minus=92	max_left_plus=95	max_min_left=49	max_min_left_minus=48	max_min_left_plus=49	max_min_right=49	max_min_right_minus=49	max_min_right_plus=49	max_pos_hash_score=196	max_right=95	max_right_minus=95	max_right_plus=95	neg_log10_pos_hash_p_value=0.3	new_junction_coverage=1.09	new_junction_read_count=104	new_junction_reference_weighted_read_count=0.00	new_junction_weighted_read_count=104.00	pos_hash_score=74	prediction=consensus	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=0.04	side_1_gene_name=yfaA/gyrA	side_1_gene_position=intergenic (-132/+17)	side_1_gene_product=DUF2138 domain-containing protein/DNA topoisomerase (ATP-hydrolyzing) subunit A	side_1_gene_strand=</<	side_1_locus_tag=ECB_RS11360/ECB_RS11365	side_1_overlap=2	side_1_possible_overlap_registers=98	side_1_possible_overlap_registers_before_trimming=100	side_1_read_count=4	side_1_redundant=0	side_1_weighted_read_count=4.00	side_2_annotate_key=repeat	side_2_continuation=0	side_2_coverage=NA	side_2_gene_name=ECB_RS25215/rbsD	side_2_gene_position=intergenic (+24/-164)	side_2_gene_product=IS3-like element IS150 family transposase/D-ribose pyranase	side_2_gene_strand=>/>	side_2_locus_tag=ECB_RS25215/ECB_RS19205	side_2_overlap=0	side_2_possible_overlap_registers=NA	side_2_possible_overlap_registers_before_trimming=NA	side_2_read_count=NA	side_2_redundant=1	side_2_weighted_read_count=NA	total_non_overlap_reads=104
JC	52	.	NC_012967	2448490	1	NC_012967	2774193	-1	0	alignment_overlap=3	coverage_minus=51	coverage_plus=43	flanking_left=101	flanking_right=101	frequency=9.587e-01	frequency_lower=9.079e-01	frequency_upper=9.859e-01	junction_effective_depth=96.32	junction_mixture_iterations=5	junction_possible_overlap_registers=95	junction_possible_overlap_registers_before_trimming=97	key=NC_012967__2448490__1__NC_012967__2774196__-1__3____101__101__0__1	max_left=97	max_left_minus=96	max_left_plus=97	max_min_left=48	max_min_left_minus=48	max_min_left_plus=47	max_min_right=49	max_min_right_minus=49	max_min_right_plus=48	max_pos_hash_score=194	max_right=93	max_right_minus=93	max_right_plus=92	neg_log10_pos_hash_p_value=0.5	new_junction_coverage=0.95	new_junction_read_count=92	new_junction_reference_weighted_read_count=0.00	new_junction_weighted_read_count=92.00	pos_hash_score=67	prediction=consensus	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=0.08	side_1_gene_name=nupC/pdeA	side_1_gene_position=intergenic (+21/+29)	side_1_gene_product=nucleoside permease NupC/bifunctional diguanylate cyclase/phosphodiesterase	side_1_gene_strand=>/<	side_1_locus_tag=ECB_RS12115/ECB_RS12120	side_1_overlap=3	side_1_possible_overlap_registers=96	side_1_possible_overlap_registers_before_trimming=100	side_1_read_count=8	side_1_redundant=0	side_1_weighted_read_count=4.17	side_2_annotate_key=repeat	side_2_continuation=0	side_2_coverage=NA	side_2_gene_name=ECB_RS13760/ECB_RS13765	side_2_gene_position=intergenic (+170/+23)	side_2_gene_product=IS4-like element IS421 family transposase/Hok/Gef family protein	side_2_gene_strand=>/<	side_2_locus_tag=ECB_RS13760/ECB_RS13765	side_2_overlap=0	side_2_possible_overlap_registers=NA	side_2_possible_overlap_registers_before_trimming=NA	side_2_read_count=NA	side_2_redundant=1	side_2_weighted_read_count=NA	total_non_overlap_reads=94
JC	53	.	NC_012967	2448501	-1	NC_012967	2772857	1	0	alignment_overlap=3	coverage_minus=36	coverage_plus=60	flanking_left=101	flanking_right=101	frequency=8.878e-01	frequency_lower=8.243e-01	frequency_upper=9.342e-01	junction_effective_depth=106.06	junction_mixture_iterations=5	junction_possible_overlap_registers=94	junction_possible_overlap_registers_before_trimming=97	key=NC_012967__2448501__-1__NC_012967__2772854__1__3____101__101__0__1	max_left=97	max_left_minus=97	max_left_plus=97	max_min_left=48	max_min_left_minus=48	max_min_left_plus=44	max_min_right=48	max_min_right_minus=44	max_min_right_plus=48	max_pos_hash_score=194	max_right=90	max_right_minus=90	max_right_plus=90	neg_log10_pos_hash_p_value=0.3	new_junction_coverage=0.97	new_junction_read_count=93	new_junction_reference_weighted_read_count=0.00	new_junction_weighted_read_count=93.00	pos_hash_score=70	prediction=consensus	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=0.17	side_1_gene_name=nupC/pdeA	side_1_gene_position=intergenic (+32/+18)	side_1_gene_product=nucleoside permease NupC/bifunctional diguanylate cyclase/phosphodiesterase	side_1_gene_strand=>/<	side_1_locus_tag=ECB_RS12115/ECB_RS12120	side_1_overlap=3	side_1_possible_overlap_registers=96	side_1_possible_overlap_registers_before_trimming=100	side_1_read_count=17	side_1_redundant=0	side_1_weighted_read_count=12.56	side_2_annotate_key=repeat	side_2_continuation=0	side_2_coverage=NA	side_2_gene_name=iap/ECB_RS13760	side_2_gene_position=intergenic (+378/-54)	side_2_gene_product=alkaline phosphatase isozyme conversion aminopeptidase/IS4-like element IS421 family transposase	side_2_gene_strand=>/>	side_2_locus_tag=ECB_RS13755/ECB_RS13760	side_2_overlap=0	side_2_possible_overlap_registers=NA	side_2_possible_overlap_registers_before_trimming=NA	side_2_read_count=NA	side_2_redundant=1	side_2_weighted_read_count=NA	total_non_overlap_reads=96
JC	54	.	NC_012967	2749274	1	NC_012967	2775876	-1	0	alignment_overlap=1	coverage_minus=52	coverage_plus=63	flanking_left=101	flanking_right=101	frequency=6.035e-01	frequency_lower=5.412e-01	frequency_upper=6.633e-01	junction_effective_depth=188.00	junction_mixture_iterations=2	junction_possible_overlap_registers=97	junction_possible_overlap_registers_before_trimming=99	key=NC_012967__2749274__1__NC_012967__2775877__-1__1____101__101__0__1	max_left=99	max_left_minus=99	max_left_plus=98	max_min_left=49	max_min_left_minus=47	max_min_left_plus=49	max_min_right=50	max_min_right_minus=49	max_min_right_plus=50	max_pos_hash_score=198	max_right=98	max_right_minus=97	max_right_plus=98	neg_log10_pos_hash_p_value=0.1	new_junction_coverage=1.15	new_junction_read_count=113	new_junction_reference_weighted_read_count=0.00	new_junction_weighted_read_count=113.00	pos_hash_score=84	prediction=consensus	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=0.75	side_1_gene_name=flhA	side_1_gene_position=coding (77/2079 nt)	side_1_gene_product=formate hydrogenlyase transcriptional activator FlhA	side_1_gene_strand=>	side_1_locus_tag=ECB_RS13645	side_1_overlap=1	side_1_possible_overlap_registers=98	side_1_possible_overlap_registers_before_trimming=100	side_1_read_count=75	side_1_redundant=0	side_1_weighted_read_count=75.00	side_2_annotate_key=repeat	side_2_continuation=0	side_2_coverage=NA	side_2_gene_name=ECB_RS25125/cysH	side_2_gene_position=intergenic (-46/+203)	side_2_gene_product=IS3-like element IS150 family transposase/phosphoadenosine phosphosulfate reductase	side_2_gene_strand=</<	side_2_locus_tag=ECB_RS25125/ECB_RS13780	side_2_overlap=0	side_2_possible_overlap_registers=NA	side_2_possible_overlap_registers_before_trimming=NA	side_2_read_count=NA	side_2_redundant=1	side_2_weighted_read_count=NA	total_non_overlap_reads=115
JC	55	.	NC_012967	3652533	-1	NC_012967	4588156	1	0	alignment_overlap=0	coverage_minus=53	coverage_plus=85	flanking_left=101	flanking_right=101	frequency=9.927e-01	frequency_lower=9.662e-01	frequency_upper=9.996e-01	junction_effective_depth=139.00	junction_mixture_iterations=1	junction_possible_overlap_registers=98	junction_possible_overlap_registers_before_trimming=100	key=NC_012967__3652533__-1__NC_012967__4588156__1__0____101__101__1__0	max_left=99	max_left_minus=99	max_left_plus=98	max_min_left=50	max_min_left_minus=49	max_min_left_plus=50	max_min_right=50	max_min_right_minus=44	max_min_right_plus=50	max_pos_hash_score=200	max_right=98	max_right_minus=98	max_right_plus=95	neg_log10_pos_hash_p_value=0.1	new_junction_coverage=1.39	new_junction_read_count=138	new_junction_reference_weighted_read_count=0.00	new_junction_weighted_read_count=138.00	pos_hash_score=87	prediction=consensus	side_1_annotate_key=repeat	side_1_continuation=0	side_1_coverage=NA	side_1_gene_name=ECB_RS25190/glyS	side_1_gene_position=intergenic (+26/+253)	side_1_gene_product=IS3-like element IS150 family transposase/glycine--tRNA ligase subunit beta	side_1_gene_strand=>/<	side_1_locus_tag=ECB_RS25190/ECB_RS18060	side_1_overlap=0	side_1_possible_overlap_registers=NA	side_1_possible_overlap_registers_before_trimming=NA	side_1_read_count=NA	side_1_redundant=1	side_1_weighted_read_count=NA	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=0.01	side_2_gene_name=opgB	side_2_gene_position=coding (520/2292 nt)	side_2_gene_product=phosphatidylglycerol--membrane-oligosaccharide glycerophosphotransferase	side_2_gene_strand=<	side_2_locus_tag=ECB_RS22570	side_2_overlap=0	side_2_possible_overlap_registers=97	side_2_possible_overlap_registers_before_trimming=100	side_2_read_count=1	side_2_redundant=0	side_2_weighted_read_count=1.00	total_non_overlap_reads=138
JC	56	.	NC_012967	3893554	1	NC_012967	4588158	-1	0	alignment_overlap=0	coverage_minus=57	coverage_plus=48	flanking_left=101	flanking_right=101	frequency=9.905e-01	frequency_lower=9.556e-01	frequency_upper=9.995e-01	junction_effective_depth=105.00	junction_mixture_iterations=1	junction_possible_overlap_registers=98	junction_possible_overlap_registers_before_trimming=100	key=NC_012967__3893554__1__NC_012967__4588158__-1__0____101__101__1__0	max_left=99	max_left_minus=96	max_left_plus=99	max_min_left=50	max_min_left_minus=49	max_min_left_plus=50	max_min_right=49	max_min_right_minus=49	max_min_right_plus=49	max_pos_hash_score=200	max_right=100	max_right_minus=98	max_right_plus=100	neg_log10_pos_hash_p_value=0.2	new_junction_coverage=1.05	new_junction_read_count=104	new_junction_reference_weighted_read_count=0.00	new_junction_weighted_read_count=104.00	pos_hash_score=79	prediction=consensus	side_1_annotate_key=repeat	side_1_continuation=0	side_1_coverage=NA	side_1_gene_name=kup/ECB_RS25215	side_1_gene_position=intergenic (+9/-47)	side_1_gene_product=low affinity potassium transporter Kup/IS3-like element IS150 family transposase	side_1_gene_strand=>/>	side_1_locus_tag=ECB_RS19190/ECB_RS25215	side_1_overlap=0	side_1_possible_overlap_registers=NA	side_1_possible_overlap_registers_before_trimming=NA	side_1_read_count=NA	side_1_redundant=1	side_1_weighted_read_count=NA	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=0.01	side_2_gene_name=opgB	side_2_gene_position=coding (518/2292 nt)	side_2_gene_product=phosphatidylglycerol--membrane-oligosaccharide glycerophosphotransferase	side_2_gene_strand=<	side_2_locus_tag=ECB_RS22570	side_2_overlap=0	side_2_possible_overlap_registers=98	side_2_possible_overlap_registers_before_trimming=100	side_2_read_count=1	side_2_redundant=0	side_2_weighted_read_count=1.00	total_non_overlap_reads=105
JC	57	.	NC_012967	3893556	1	NC_012967	4252530	-1	0	alignment_overlap=2	coverage_minus=56	coverage_plus=49	flanking_left=101	flanking_right=101	frequency=9.545e-01	frequency_lower=9.068e-01	frequency_upper=9.819e-01	junction_effective_depth=110.00	junction_mixture_iterations=2	junction_possible_overlap_registers=96	junction_possible_overlap_registers_before_trimming=98	key=NC_012967__3893554__1__NC_012967__4252530__-1__2____101__101__1__0	max_left=96	max_left_minus=94	max_left_plus=96	max_min_left=49	max_min_left_minus=47	max_min_left_plus=49	max_min_right=49	max_min_right_minus=49	max_min_right_plus=48	max_pos_hash_score=196	max_right=97	max_right_minus=97	max_right_plus=97	neg_log10_pos_hash_p_value=0.2	new_junction_coverage=1.08	new_junction_read_count=105	new_junction_reference_weighted_read_count=0.00	new_junction_weighted_read_count=105.00	pos_hash_score=77	prediction=consensus	side_1_annotate_key=repeat	side_1_continuation=0	side_1_coverage=NA	side_1_gene_name=kup/ECB_RS25215	side_1_gene_position=intergenic (+11/-45)	side_1_gene_product=low affinity potassium transporter Kup/IS3-like element IS150 family transposase	side_1_gene_strand=>/>	side_1_locus_tag=ECB_RS19190/ECB_RS25215	side_1_overlap=0	side_1_possible_overlap_registers=NA	side_1_possible_overlap_registers_before_trimming=NA	side_1_read_count=NA	side_1_redundant=1	side_1_weighted_read_count=NA	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=0.05	side_2_gene_name=uvrA	side_2_gene_position=coding (274/2823 nt)	side_2_gene_product=excinuclease ABC subunit UvrA	side_2_gene_strand=<	side_2_locus_tag=ECB_RS20910	side_2_overlap=2	side_2_possible_overlap_registers=96	side_2_possible_overlap_registers_before_trimming=100	side_2_read_count=5	side_2_redundant=0	side_2_weighted_read_count=5.00	total_non_overlap_reads=105
JC	58	.	NC_012967	3894995	-1	NC_012967	3900623	1	0	alignment_overlap=1	coverage_minus=56	coverage_plus=74	flanking_left=101	flanking_right=101	frequency=9.922e-01	frequency_lower=9.639e-01	frequency_upper=9.996e-01	junction_effective_depth=130.00	junction_mixture_iterations=1	junction_possible_overlap_registers=97	junction_possible_overlap_registers_before_trimming=99	key=NC_012967__3894996__-1__NC_012967__3900623__1__1____101__101__1__0	max_left=98	max_left_minus=98	max_left_plus=97	max_min_left=49	max_min_left_minus=48	max_min_left_plus=49	max_min_right=50	max_min_right_minus=50	max_min_right_plus=49	max_pos_hash_score=198	max_right=98	max_right_minus=85	max_right_plus=98	neg_log10_pos_hash_p_value=0.0	new_junction_coverage=1.31	new_junction_read_count=129	new_junction_reference_weighted_read_count=0.00	new_junction_weighted_read_count=129.00	pos_hash_score=97	prediction=consensus	side_1_annotate_key=repeat	side_1_continuation=0	side_1_coverage=NA	side_1_gene_name=ECB_RS25215/rbsD	side_1_gene_position=intergenic (+25/-163)	side_1_gene_product=IS3-like element IS150 family transposase/D-ribose pyranase	side_1_gene_strand=>/>	side_1_locus_tag=ECB_RS25215/ECB_RS19205	side_1_overlap=0	side_1_possible_overlap_registers=NA	side_1_possible_overlap_registers_before_trimming=NA	side_1_read_count=NA	side_1_redundant=1	side_1_weighted_read_count=NA	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=0.01	side_2_gene_name=rbsR	side_2_gene_position=coding (590/993 nt)	side_2_gene_product=ribose operon transcriptional repressor RbsR	side_2_gene_strand=>	side_2_locus_tag=ECB_RS19230	side_2_overlap=1	side_2_possible_overlap_registers=96	side_2_possible_overlap_registers_before_trimming=100	side_2_read_count=1	side_2_redundant=0	side_2_weighted_read_count=1.00	total_non_overlap_reads=130
JC	59	.	NC_012967	3894996	-1	NC_012967	4252526	1	0	alignment_overlap=0	coverage_minus=53	coverage_plus=74	flanking_left=101	flanking_right=101	frequency=9.831e-01	frequency_lower=9.484e-01	frequency_upper=9.969e-01	junction_effective_depth=121.00	junction_mixture_iterations=1	junction_possible_overlap_registers=97	junction_possible_overlap_registers_before_trimming=100	key=NC_012967__3894996__-1__NC_012967__4252526__1__0____101__101__1__0	max_left=98	max_left_minus=95	max_left_plus=98	max_min_left=49	max_min_left_minus=49	max_min_left_plus=46	max_min_right=50	max_min_right_minus=48	max_min_right_plus=50	max_pos_hash_score=200	max_right=100	max_right_minus=100	max_right_plus=94	neg_log10_pos_hash_p_value=0.1	new_junction_coverage=1.21	new_junction_read_count=119	new_junction_reference_weighted_read_count=0.00	new_junction_weighted_read_count=119.00	pos_hash_score=81	prediction=consensus	side_1_annotate_key=repeat	side_1_continuation=0	side_1_coverage=NA	side_1_gene_name=ECB_RS25215/rbsD	side_1_gene_position=intergenic (+26/-162)	side_1_gene_product=IS3-like element IS150 family transposase/D-ribose pyranase	side_1_gene_strand=>/>	side_1_locus_tag=ECB_RS25215/ECB_RS19205	side_1_overlap=0	side_1_possible_overlap_registers=NA	side_1_possible_overlap_registers_before_trimming=NA	side_1_read_count=NA	side_1_redundant=1	side_1_weighted_read_count=NA	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=0.02	side_2_gene_name=uvrA	side_2_gene_position=coding (278/2823 nt)	side_2_gene_product=excinuclease ABC subunit UvrA	side_2_gene_strand=<	side_2_locus_tag=ECB_RS20910	side_2_overlap=0	side_2_possible_overlap_registers=95	side_2_possible_overlap_registers_before_trimming=100	side_2_read_count=2	side_2_redundant=0	side_2_weighted_read_count=2.00	total_non_overlap_reads=127
UN	60	.	NC_012967	15623	15623
UN	61	.	NC_012967	15625	15625
UN	62	.	NC_012967	15628	15628
UN	63	.	NC_012967	15634	16462
UN	64	.	NC_012967	16465	16465
UN	65	.	NC_012967	23597	23597
UN	66	.	NC_012967	23599	23602
UN	67	.	NC_012967	23604	23604
UN	68	.	NC_012967	23606	23608
UN	69	.	NC_012967	23611	23613
UN	70	.	NC_012967	23615	23615
UN	71	.	NC_012967	23619	23622
UN	72	.	NC_012967	23624	23797
UN	73	.	NC_012967	227283	227286
UN	74	.	NC_012967	227289	228274
UN	75	.	NC_012967	228276	228277
UN	76	.	NC_012967	228279	228279
UN	77	.	NC_012967	228282	228282
UN	78	.	NC_012967	228284	228284
UN	79	.	NC_012967	228286	228286
UN	80	.	NC_012967	228291	228292
UN	81	.	NC_012967	228295	228300
UN	82	.	NC_012967	228303	228303
UN	83	.	NC_012967	228308	228313
UN	84	.	NC_012967	228315	228316
UN	85	.	NC_012967	228318	228318
UN	86	.	NC_012967	228323	228323
UN	87	.	NC_012967	229019	229019
UN	88	.	NC_012967	229023	231514
UN	89	.	NC_012967	231519	231520
UN	90	.	NC_012967	231528	231530
UN	91	.	NC_012967	231535	231537
UN	92	.	NC_012967	231547	231547
UN	93	.	NC_012967	241544	241544
UN	94	.	NC_012967	241546	241546
UN	95	.	NC_012967	241549	241688
UN	96	.	NC_012967	241690	241690
UN	97	.	NC_012967	241692	241694
UN	98	.	NC_012967	241696	241696
UN	99	.	NC_012967	241698	241698
UN	100	.	NC_012967	241700	241700
UN	101	.	NC_012967	241702	241702
UN	102	.	NC_012967	263096	263103
UN	103	.	NC_012967	263105	263105
UN	104	.	NC_012967	263120	263221
UN	105	.	NC_012967	353535	353535
UN	106	.	NC_012967	353537	353575
UN	107	.	NC_012967	353578	353578
UN	108	.	NC_012967	353580	353582
UN	109	.	NC_012967	353584	353586
UN	110	.	NC_012967	353588	353658
UN	111	.	NC_012967	353664	353667
UN	112	.	NC_012967	353670	353698
UN	113	.	NC_012967	353701	353702
UN	114	.	NC_012967	353716	353716
UN	115	.	NC_012967	377033	377033
UN	116	.	NC_012967	377035	377035
UN	117	.	NC_012967	377038	377038
UN	118	.	NC_012967	377050	377051
UN	119	.	NC_012967	377053	377056
UN	120	.	NC_012967	377058	377059
UN	121	.	NC_012967	377062	377177
UN	122	.	NC_012967	429787	429787
UN	123	.	NC_012967	429789	429894
UN	124	.	NC_012967	429948	429953
UN	125	.	NC_012967	429955	429962
UN	126	.	NC_012967	429964	429970
UN	127	.	NC_012967	429972	429972
UN	128	.	NC_012967	429974	430553
UN	129	.	NC_012967	430555	430555
UN	130	.	NC_012967	430557	430558
UN	131	.	NC_012967	430560	430560
UN	132	.	NC_012967	430564	430564
UN	133	.	NC_012967	430567	430567
UN	134	.	NC_012967	430569	430574
UN	135	.	NC_012967	430578	430578
UN	136	.	NC_012967	430581	430582
UN	137	.	NC_012967	474383	474391
UN	138	.	NC_012967	495699	495699
UN	139	.	NC_012967	495701	498075
UN	140	.	NC_012967	498087	498090
UN	141	.	NC_012967	498092	498096
UN	142	.	NC_012967	498098	498098
UN	143	.	NC_012967	498100	498100
UN	144	.	NC_012967	498102	498102
UN	145	.	NC_012967	498104	498106
UN	146	.	NC_012967	498108	498136
UN	147	.	NC_012967	498573	498573
UN	148	.	NC_012967	498575	498664
UN	149	.	NC_012967	498666	498666
UN	150	.	NC_012967	498668	498669
UN	151	.	NC_012967	498673	498681
UN	152	.	NC_012967	498687	498687
UN	153	.	NC_012967	547191	547191
UN	154	.	NC_012967	547219	547383
UN	155	.	NC_012967	553535	553536
UN	156	.	NC_012967	553613	553613
UN	157	.	NC_012967	553615	553619
UN	158	.	NC_012967	553621	553621
UN	159	.	NC_012967	553623	554444
UN	160	.	NC_012967	588716	588817
UN	161	.	NC_012967	588820	588821
UN	162	.	NC_012967	588823	588823
UN	163	.	NC_012967	588826	588827
UN	164	.	NC_012967	588829	588829
UN	165	.	NC_012967	588831	588832
UN	166	.	NC_012967	588834	588834
UN	167	.	NC_012967	588836	588836
UN	168	.	NC_012967	588839	588839
UN	169	.	NC_012967	588843	588843
UN	170	.	NC_012967	588845	588846
UN	171	.	NC_012967	588849	588849
UN	172	.	NC_012967	588851	588853
UN	173	.	NC_012967	588855	589253
UN	174	.	NC_012967	589255	589257
UN	175	.	NC_012967	589259	589260
UN	176	.	NC_012967	589262	589263
UN	177	.	NC_012967	589304	589304
UN	178	.	NC_012967	589306	589310
UN	179	.	NC_012967	589897	590016
UN	180	.	NC_012967	590924	590936
UN	181	.	NC_012967	590965	590965
UN	182	.	NC_012967	591020	591020
UN	183	.	NC_012967	591022	591023
UN	184	.	NC_012967	591025	591031
UN	185	.	NC_012967	591033	591767
UN	186	.	NC_012967	619392	619392
UN	187	.	NC_012967	619399	619403
UN	188	.	NC_012967	619405	619408
UN	189	.	NC_012967	619410	619410
UN	190	.	NC_012967	619412	619564
UN	191	.	NC_012967	634237	634387
UN	192	.	NC_012967	634389	634392
UN	193	.	NC_012967	634394	634397
UN	194	.	NC_012967	634399	634399
UN	195	.	NC_012967	634401	634406
UN	196	.	NC_012967	634409	634465
UN	197	.	NC_012967	664909	664911
UN	198	.	NC_012967	664914	665786
UN	199	.	NC_012967	665789	665791
UN	200	.	NC_012967	665794	665819
UN	201	.	NC_012967	665821	665823
UN	202	.	NC_012967	665827	665827
UN	203	.	NC_012967	1110769	1110769
UN	204	.	NC_012967	1110772	1110772
UN	205	.	NC_012967	1110774	1110776
UN	206	.	NC_012967	1110778	1110780
UN	207	.	NC_012967	1110782	1110782
UN	208	.	NC_012967	1110789	1111468
UN	209	.	NC_012967	1111472	1111472
UN	210	.	NC_012967	1111475	1111475
UN	211	.	NC_012967	1111477	1111478
UN	212	.	NC_012967	1111480	1111496
UN	213	.	NC_012967	1193811	1193811
UN	214	.	NC_012967	1193813	1193813
UN	215	.	NC_012967	1193815	1193816
UN	216	.	NC_012967	1193822	1193823
UN	217	.	NC_012967	1193826	1194115
UN	218	.	NC_012967	1194118	1194118
UN	219	.	NC_012967	1322975	1322983
UN	220	.	NC_012967	1322989	1322993
UN	221	.	NC_012967	1322997	1322997
UN	222	.	NC_012967	1322999	1323000
UN	223	.	NC_012967	1323005	1323010
UN	224	.	NC_012967	1323013	1323221
UN	225	.	NC_012967	1442872	1442872
UN	226	.	NC_012967	1442875	1442875
UN	227	.	NC_012967	1442877	1442878
UN	228	.	NC_012967	1442880	1442880
UN	229	.	NC_012967	1442888	1442889
UN	230	.	NC_012967	1442891	1442891
UN	231	.	NC_012967	1442894	1442894
UN	232	.	NC_012967	1442896	1442897
UN	233	.	NC_012967	1442899	1442903
UN	234	.	NC_012967	1442906	1442910
UN	235	.	NC_012967	1442912	1442916
UN	236	.	NC_012967	1442918	1442918
UN	237	.	NC_012967	1442920	1442922
UN	238	.	NC_012967	1442925	1442925
UN	239	.	NC_012967	1442927	1442927
UN	240	.	NC_012967	1442930	1442930
UN	241	.	NC_012967	1442932	1443651
UN	242	.	NC_012967	1443656	1443659
UN	243	.	NC_012967	1460440	1460441
UN	244	.	NC_012967	1460443	1460444
UN	245	.	NC_012967	1460446	1460645
UN	246	.	NC_012967	1500483	1500483
UN	247	.	NC_012967	1500485	1500488
UN	248	.	NC_012967	1500490	1500490
UN	249	.	NC_012967	1500493	1500499
UN	250	.	NC_012967	1500505	1500505
UN	251	.	NC_012967	1500508	1502876
UN	252	.	NC_012967	1502878	1502879
UN	253	.	NC_012967	1502884	1502887
UN	254	.	NC_012967	1502889	1502889
UN	255	.	NC_012967	1502891	1502891
UN	256	.	NC_012967	1502893	1502894
UN	257	.	NC_012967	1502900	1502900
UN	258	.	NC_012967	1502903	1502903
UN	259	.	NC_012967	1502905	1502907
UN	260	.	NC_012967	1502909	1502909
UN	261	.	NC_012967	1502919	1502919
UN	262	.	NC_012967	1502922	1502922
UN	263	.	NC_012967	1502925	1502925
UN	264	.	NC_012967	1502927	1502927
UN	265	.	NC_012967	1502929	1502929
UN	266	.	NC_012967	1502936	1502939
UN	267	.	NC_012967	1502941	1502941
UN	268	.	NC_012967	1502944	1502957
UN	269	.	NC_012967	1503407	1503410
UN	270	.	NC_012967	1503412	1503412
UN	271	.	NC_012967	1503415	1503415
UN	272	.	NC_012967	1503426	1503430
UN	273	.	NC_012967	1503433	1503435
UN	274	.	NC_012967	1503438	1503439
UN	275	.	NC_012967	1503462	1503463
UN	276	.	NC_012967	1503466	1503467
UN	277	.	NC_012967	1503469	1503470
UN	278	.	NC_012967	1503474	1503479
UN	279	.	NC_012967	1503483	1503484
UN	280	.	NC_012967	1503541	1503549
UN	281	.	NC_012967	1606842	1607647
UN	282	.	NC_012967	1607649	1607649
UN	283	.	NC_012967	1607651	1607653
UN	284	.	NC_012967	1607684	1607684
UN	285	.	NC_012967	1607688	1607688
UN	286	.	NC_012967	1607691	1607691
UN	287	.	NC_012967	1608190	1608191
UN	288	.	NC_012967	1608194	1608194
UN	289	.	NC_012967	1608197	1608199
UN	290	.	NC_012967	1608201	1608203
UN	291	.	NC_012967	1608212	1608212
UN	292	.	NC_012967	1608227	1608229
UN	293	.	NC_012967	1608232	1608892
UN	294	.	NC_012967	1608894	1608894
UN	295	.	NC_012967	1608903	1608903
UN	296	.	NC_012967	1615737	1615737
UN	297	.	NC_012967	1615740	1615760
UN	298	.	NC_012967	1615762	1615763
UN	299	.	NC_012967	1615765	1616457
UN	300	.	NC_012967	1616459	1616470
UN	301	.	NC_012967	1616472	1616472
UN	302	.	NC_012967	1616476	1616476
UN	303	.	NC_012967	1857382	1858140
UN	304	.	NC_012967	1973522	1973522
UN	305	.	NC_012967	1973530	1973531
UN	306	.	NC_012967	1973533	1973533
UN	307	.	NC_012967	1973565	1973744
UN	308	.	NC_012967	1973748	1973748
UN	309	.	NC_012967	1973750	1973751
UN	310	.	NC_012967	1973755	1973755
UN	311	.	NC_012967	1973757	1973759
UN	312	.	NC_012967	1973764	1973764
UN	313	.	NC_012967	2034623	2034777
UN	314	.	NC_012967	2100308	2122436
UN	315	.	NC_012967	2128878	2129110
UN	316	.	NC_012967	2143313	2143476
UN	317	.	NC_012967	2143478	2143478
UN	318	.	NC_012967	2143480	2143480
UN	319	.	NC_012967	2143482	2143484
UN	320	.	NC_012967	2143488	2143488
UN	321	.	NC_012967	2143490	2143539
UN	322	.	NC_012967	2143542	2143542
UN	323	.	NC_012967	2143544	2143545
UN	324	.	NC_012967	2143547	2143547
UN	325	.	NC_012967	2254540	2254549
UN	326	.	NC_012967	2254820	2254820
UN	327	.	NC_012967	2254822	2254990
UN	328	.	NC_012967	2254993	2254993
UN	329	.	NC_012967	2263011	2263011
UN	330	.	NC_012967	2263015	2263015
UN	331	.	NC_012967	2263020	2263215
UN	332	.	NC_012967	2263218	2263218
UN	333	.	NC_012967	2263220	2263221
UN	334	.	NC_012967	2407737	2407738
UN	335	.	NC_012967	2407743	2407743
UN	336	.	NC_012967	2407748	2407752
UN	337	.	NC_012967	2407755	2407909
UN	338	.	NC_012967	2618357	2618360
UN	339	.	NC_012967	2618362	2618363
UN	340	.	NC_012967	2618368	2618372
UN	341	.	NC_012967	2618376	2618509
UN	342	.	NC_012967	2647761	2647762
UN	343	.	NC_012967	2647764	2647764
UN	344	.	NC_012967	2647766	2647766
UN	345	.	NC_012967	2647768	2647769
UN	346	.	NC_012967	2647771	2647773
UN	347	.	NC_012967	2647776	2647779
UN	348	.	NC_012967	2647783	2647783
UN	349	.	NC_012967	2647785	2652451
UN	350	.	NC_012967	2683017	2683022
UN	351	.	NC_012967	2773152	2773152
UN	352	.	NC_012967	2773157	2773878
UN	353	.	NC_012967	2773881	2773881
UN	354	.	NC_012967	2773884	2773885
UN	355	.	NC_012967	2773888	2773890
UN	356	.	NC_012967	2773894	2773895
UN	357	.	NC_012967	2773897	2773899
UN	358	.	NC_012967	2773902	2773903
UN	359	.	NC_012967	2773906	2773906
UN	360	.	NC_012967	2773908	2773909
UN	361	.	NC_012967	2773911	2773911
UN	362	.	NC_012967	2774779	2775578
UN	363	.	NC_012967	2775580	2775580
UN	364	.	NC_012967	2775582	2775586
UN	365	.	NC_012967	2822509	2822511
UN	366	.	NC_012967	2822513	2822513
UN	367	.	NC_012967	2822515	2822634
UN	368	.	NC_012967	2822637	2822639
UN	369	.	NC_012967	2822658	2822661
UN	370	.	NC_012967	2882974	2882975
UN	371	.	NC_012967	2882979	2882979
UN	372	.	NC_012967	2882981	2882983
UN	373	.	NC_012967	2882986	2882986
UN	374	.	NC_012967	2882989	2883641
UN	375	.	NC_012967	2883643	2883645
UN	376	.	NC_012967	2883653	2883653
UN	377	.	NC_012967	2883657	2883659
UN	378	.	NC_012967	2883665	2883665
UN	379	.	NC_012967	3024251	3024252
UN	380	.	NC_012967	3024254	3024255
UN	381	.	NC_012967	3024259	3024420
UN	382	.	NC_012967	3024423	3024423
UN	383	.	NC_012967	3024425	3024425
UN	384	.	NC_012967	3351910	3351911
UN	385	.	NC_012967	3351923	3351944
UN	386	.	NC_012967	3351946	3351946
UN	387	.	NC_012967	3351948	3354021
UN	388	.	NC_012967	3354023	3354026
UN	389	.	NC_012967	3354030	3354030
UN	390	.	NC_012967	3354032	3354052
UN	391	.	NC_012967	3354680	3354694
UN	392	.	NC_012967	3354696	3354698
UN	393	.	NC_012967	3354701	3354701
UN	394	.	NC_012967	3354705	3354705
UN	395	.	NC_012967	3354707	3354714
UN	396	.	NC_012967	3354716	3354716
UN	397	.	NC_012967	3354718	3354718
UN	398	.	NC_012967	3354748	3354750
UN	399	.	NC_012967	3354756	3354757
UN	400	.	NC_012967	3354759	3354759
UN	401	.	NC_012967	3354761	3354761
UN	402	.	NC_012967	3354763	3354765
UN	403	.	NC_012967	3354767	3354772
UN	404	.	NC_012967	3354774	3354774
UN	405	.	NC_012967	3354776	3354778
UN	406	.	NC_012967	3354780	3354781
UN	407	.	NC_012967	3354784	3354784
UN	408	.	NC_012967	3354786	3356401
UN	409	.	NC_012967	3356405	3356406
UN	410	.	NC_012967	3433750	3433944
UN	411	.	NC_012967	3433950	3433951
UN	412	.	NC_012967	3433957	3433958
UN	413	.	NC_012967	3550261	3551066
UN	414	.	NC_012967	3551072	3551072
UN	415	.	NC_012967	3551075	3551075
UN	416	.	NC_012967	3551636	3551636
UN	417	.	NC_012967	3551638	3551638
UN	418	.	NC_012967	3551656	3551656
UN	419	.	NC_012967	3551660	3551661
UN	420	.	NC_012967	3551663	3551666
UN	421	.	NC_012967	3551668	3551699
UN	422	.	NC_012967	3551705	3551705
UN	423	.	NC_012967	3551708	3552007
UN	424	.	NC_012967	3552015	3552015
UN	425	.	NC_012967	3552017	3552017
UN	426	.	NC_012967	3552031	3552031
UN	427	.	NC_012967	3651352	3651352
UN	428	.	NC_012967	3651354	3651356
UN	429	.	NC_012967	3651359	3652226
UN	430	.	NC_012967	3652230	3652230
UN	431	.	NC_012967	3652232	3652234
UN	432	.	NC_012967	3697313	3697313
UN	433	.	NC_012967	3697316	3698046
UN	434	.	NC_012967	3698048	3698078
UN	435	.	NC_012967	3698080	3698084
UN	436	.	NC_012967	3698601	3698933
UN	437	.	NC_012967	3698938	3698938
UN	438	.	NC_012967	3698942	3698942
UN	439	.	NC_012967	3698945	3698948
UN	440	.	NC_012967	3698952	3698952
UN	441	.	NC_012967	3698954	3698955
UN	442	.	NC_012967	3698957	3699036
UN	443	.	NC_012967	3699038	3699038
UN	444	.	NC_012967	3699041	3699042
UN	445	.	NC_012967	3699045	3699045
UN	446	.	NC_012967	3699049	3699049
UN	447	.	NC_012967	3699053	3699053
UN	448	.	NC_012967	3741469	3741469
UN	449	.	NC_012967	3741504	3741650
UN	450	.	NC_012967	3741652	3741652
UN	451	.	NC_012967	3741654	3741659
UN	452	.	NC_012967	3893806	3893809
UN	453	.	NC_012967	3893812	3893812
UN	454	.	NC_012967	3893821	3893821
UN	455	.	NC_012967	3893853	3894618
UN	456	.	NC_012967	3894620	3894620
UN	457	.	NC_012967	3894999	3900622
UN	458	.	NC_012967	3903743	3903745
UN	459	.	NC_012967	3903747	3903747
UN	460	.	NC_012967	3903751	3908359
UN	461	.	NC_012967	4013911	4015254
UN	462	.	NC_012967	4015257	4015263
UN	463	.	NC_012967	4015271	4015271
UN	464	.	NC_012967	4015275	4015275
UN	465	.	NC_012967	4015880	4015881
UN	466	.	NC_012967	4015883	4015899
UN	467	.	NC_012967	4015901	4015901
UN	468	.	NC_012967	4015903	4015905
UN	469	.	NC_012967	4016328	4016328
UN	470	.	NC_012967	4016330	4016330
UN	471	.	NC_012967	4016332	4016332
UN	472	.	NC_012967	4016334	4017105
UN	473	.	NC_012967	4017107	4017107
UN	474	.	NC_012967	4017117	4017117
UN	475	.	NC_012967	4017121	4017122
UN	476	.	NC_012967	4018395	4018395
UN	477	.	NC_012967	4018397	4018793
UN	478	.	NC_012967	4018803	4018803
UN	479	.	NC_012967	4146527	4146527
UN	480	.	NC_012967	4146529	4146529
UN	481	.	NC_012967	4146531	4146533
UN	482	.	NC_012967	4146535	4146536
UN	483	.	NC_012967	4146538	4146539
UN	484	.	NC_012967	4146543	4146543
UN	485	.	NC_012967	4146545	4146545
UN	486	.	NC_012967	4146547	4146549
UN	487	.	NC_012967	4146551	4146552
UN	488	.	NC_012967	4146556	4147794
UN	489	.	NC_012967	4147799	4147799
UN	490	.	NC_012967	4147802	4147802
UN	491	.	NC_012967	4147804	4147808
UN	492	.	NC_012967	4147810	4147810
UN	493	.	NC_012967	4147814	4147814
UN	494	.	NC_012967	4147817	4147818
UN	495	.	NC_012967	4147820	4147820
UN	496	.	NC_012967	4147822	4147822
UN	497	.	NC_012967	4147834	4147834
UN	498	.	NC_012967	4147837	4147838
UN	499	.	NC_012967	4147840	4147840
UN	500	.	NC_012967	4149048	4149050
UN	501	.	NC_012967	4149074	4149075
UN	502	.	NC_012967	4149078	4149078
UN	503	.	NC_012967	4149082	4149889
UN	504	.	NC_012967	4149891	4149891
UN	505	.	NC_012967	4149893	4149893
UN	506	.	NC_012967	4149901	4149901
UN	507	.	NC_012967	4150385	4150385
UN	508	.	NC_012967	4150388	4150390
UN	509	.	NC_012967	4150395	4151223
UN	510	.	NC_012967	4151226	4151226
UN	511	.	NC_012967	4151228	4151229
UN	512	.	NC_012967	4151232	4151233
UN	513	.	NC_012967	4151236	4151242
UN	514	.	NC_012967	4151245	4151246
UN	515	.	NC_012967	4151248	4151250
UN	516	.	NC_012967	4151264	4151264
UN	517	.	NC_012967	4187767	4187767
UN	518	.	NC_012967	4187769	4192570
UN	519	.	NC_012967	4192572	4192574
UN	520	.	NC_012967	4192576	4192576
UN	521	.	NC_012967	4192578	4192578
UN	522	.	NC_012967	4274870	4274870
UN	523	.	NC_012967	4505633	4505764
UN	524	.	NC_012967	4505776	4505776
UN	525	.	NC_012967	4505857	4505857
UN	526	.	NC_012967	4521795	4521795
UN	527	.	NC_012967	4521819	4521824
UN	528	.	NC_012967	4521830	4521842
UN	529	.	NC_012967	4521846	4521846
UN	530	.	NC_012967	4521848	4521849
UN	531	.	NC_012967	4521851	4521853
UN	532	.	NC_012967	4521855	4522018
UN	533	.	NC_012967	4522021	4522022
UN	534	.	NC_012967	4550992	4551148
