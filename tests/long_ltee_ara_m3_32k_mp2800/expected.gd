#=GENOME_DIFF	1.0
#=CREATED	07:47:11 31 Aug 2026
#=PROGRAM	breseq 0.50.0 revision 4ff75f7e1fd2
#=COMMAND	./src/breseq/breseq -j 7 -o tests/long_ltee_ara_m3_32k_mp2800 -r tests/long_ltee_ara_m3_32k_mp2800/../data/downloads/ltee_REL606/REL606.gbk -l 80 --junction-alignment-pair-limit 2000000 --predict-copy-number --predict-discordant-pairs --predict-missing-pairs --predict-pair-distance tests/long_ltee_ara_m3_32k_mp2800/../data/downloads/ena_SRR098033/SRR098033_1.fastq.gz tests/long_ltee_ara_m3_32k_mp2800/../data/downloads/ena_SRR098033/SRR098033_2.fastq.gz
#=REFSEQ	tests/long_ltee_ara_m3_32k_mp2800/../data/downloads/ltee_REL606/REL606.gbk
#=READSEQ	tests/long_ltee_ara_m3_32k_mp2800/../data/downloads/ena_SRR098033/SRR098033_1.fastq.gz
#=READSEQ	tests/long_ltee_ara_m3_32k_mp2800/../data/downloads/ena_SRR098033/SRR098033_2.fastq.gz
#=CONVERTED-BASES	370384980
#=CONVERTED-READS	10582428
#=INPUT-BASES	1778009800
#=INPUT-READS	50800280
#=MAPPED-BASES	335105570
#=MAPPED-READS	9583197
SNP	1	46	REL606	9972	G	aa_new_seq=T	aa_position=174	aa_ref_seq=N	codon_new_seq=ACC	codon_number=174	codon_position=2	codon_ref_seq=AAC	gene_name=yaaH	gene_position=521	gene_product=conserved inner membrane protein associated with acetate transport	gene_strand=<	genes_overlapping=yaaH	locus_tag=ECB_00010	locus_tags_overlapping=ECB_00010	mutation_category=snp_nonsynonymous	position_end=9972	position_start=9972	ref_seq=T	snp_type=nonsynonymous	transl_table=11
SNP	2	47	REL606	81158	C	aa_new_seq=A	aa_position=245	aa_ref_seq=E	codon_new_seq=GCG	codon_number=245	codon_position=2	codon_ref_seq=GAG	gene_name=setA	gene_position=734	gene_product=broad specificity sugar efflux system	gene_strand=>	genes_overlapping=setA	locus_tag=ECB_00072	locus_tags_overlapping=ECB_00072	mutation_category=snp_nonsynonymous	position_end=81158	position_start=81158	ref_seq=A	snp_type=nonsynonymous	transl_table=11
SNP	3	48	REL606	216480	T	gene_name=tilS/rof	gene_position=intergenic (+9/+40)	gene_product=tRNA(Ile)-lysidine synthetase/modulator of Rho-dependent transcription termination	gene_strand=>/<	locus_tag=ECB_00186/ECB_00187	mutation_category=snp_intergenic	position_end=216480	position_start=216480	ref_seq=C	snp_type=intergenic
SNP	4	49	REL606	247796	C	aa_new_seq=G	aa_position=27	aa_ref_seq=S	codon_new_seq=GGT	codon_number=27	codon_position=1	codon_ref_seq=AGT	gene_name=fadE	gene_position=79	gene_product=acyl-CoA dehydrogenase	gene_strand=<	genes_overlapping=fadE	locus_tag=ECB_00216	locus_tags_overlapping=ECB_00216	mutation_category=snp_nonsynonymous	position_end=247796	position_start=247796	ref_seq=T	snp_type=nonsynonymous	transl_table=11
SNP	5	50	REL606	281923	T	gene_name=ykgK/ykgL	gene_position=intergenic (-406/-369)	gene_product=predicted regulator/hypothetical protein	gene_strand=</>	locus_tag=ECB_00253/ECB_00254	mutation_category=snp_intergenic	position_end=281923	position_start=281923	ref_seq=G	snp_type=intergenic
SNP	6	51	REL606	398061	C	aa_new_seq=A	aa_position=45	aa_ref_seq=G	codon_new_seq=GCC	codon_number=45	codon_position=2	codon_ref_seq=GGC	gene_name=secF	gene_position=134	gene_product=protein export protein SecF	gene_strand=>	genes_overlapping=secF	locus_tag=ECB_00357	locus_tags_overlapping=ECB_00357	mutation_category=snp_nonsynonymous	position_end=398061	position_start=398061	ref_seq=G	snp_type=nonsynonymous	transl_table=11
INS	7	52	REL606	433366	T	gene_name=lon/hupB	gene_position=intergenic (+69/-140)	gene_product=DNA-binding ATP-dependent protease La/HU, DNA-binding transcriptional regulator, beta subunit	gene_strand=>/>	genes_promoter=hupB	locus_tag=ECB_00391/ECB_00392	locus_tags_promoter=ECB_00392	mutation_category=small_indel	position_end=433366	position_start=433366	ref_seq=T	repeat_length=1	repeat_new_copies=8	repeat_ref_copies=7	repeat_seq=T
INS	8	53,54	REL606	473904	GC	gene_name=ybaL	gene_position=coding (1402/1677 nt)	gene_product=predicted transporter with NAD(P)-binding Rossmann-fold domain	gene_strand=<	genes_overlapping=ybaL	locus_tag=ECB_00429	locus_tags_overlapping=ECB_00429	mutation_category=small_indel	position_end=473904	position_start=473904	ref_seq=C
SNP	9	55	REL606	648692	T	aa_new_seq=H	aa_position=69	aa_ref_seq=R	codon_new_seq=CAC	codon_number=69	codon_position=2	codon_ref_seq=CGC	gene_name=mrdB	gene_position=206	gene_product=cell wall shape-determining protein	gene_strand=<	genes_overlapping=mrdB	locus_tag=ECB_00603	locus_tags_overlapping=ECB_00603	mutation_category=snp_nonsynonymous	position_end=648692	position_start=648692	ref_seq=C	snp_type=nonsynonymous	transl_table=11
SNP	10	56	REL606	734488	T	aa_new_seq=T	aa_position=258	aa_ref_seq=A	codon_new_seq=ACT	codon_number=258	codon_position=1	codon_ref_seq=GCT	gene_name=gltA	gene_position=772	gene_product=citrate synthase	gene_strand=<	genes_overlapping=gltA	locus_tag=ECB_00680	locus_tags_overlapping=ECB_00680	mutation_category=snp_nonsynonymous	position_end=734488	position_start=734488	ref_seq=C	snp_type=nonsynonymous	transl_table=11
SNP	11	57	REL606	943475	T	aa_new_seq=I	aa_position=9	aa_ref_seq=M	codon_new_seq=ATA	codon_number=9	codon_position=3	codon_ref_seq=ATG	gene_name=infA	gene_position=27	gene_product=translation initiation factor IF-1	gene_strand=<	genes_overlapping=infA	locus_tag=ECB_00888	locus_tags_overlapping=ECB_00888	mutation_category=snp_nonsynonymous	position_end=943475	position_start=943475	ref_seq=C	snp_type=nonsynonymous	transl_table=11
SNP	12	58	REL606	1331794	A	aa_new_seq=K	aa_position=792	aa_ref_seq=T	codon_new_seq=AAG	codon_number=792	codon_position=2	codon_ref_seq=ACG	gene_name=topA	gene_position=2375	gene_product=DNA topoisomerase I	gene_strand=>	genes_overlapping=topA	locus_tag=ECB_01250	locus_tags_overlapping=ECB_01250	mutation_category=snp_nonsynonymous	position_end=1331794	position_start=1331794	ref_seq=C	snp_type=nonsynonymous	transl_table=11
SNP	13	59	REL606	1383157	G	aa_new_seq=G	aa_position=224	aa_ref_seq=G	codon_new_seq=GGG	codon_number=224	codon_position=3	codon_ref_seq=GGT	gene_name=ycjX	gene_position=672	gene_product=conserved protein with nucleoside triphosphate hydrolase domain	gene_strand=>	genes_overlapping=ycjX	locus_tag=ECB_01298	locus_tags_overlapping=ECB_01298	mutation_category=snp_synonymous	position_end=1383157	position_start=1383157	ref_seq=T	snp_type=synonymous	transl_table=11
MOB	14	149,152,1917	REL606	1451972	IS150	-1	3	gene_name=acpD	gene_position=coding (218-220/606 nt)	gene_product=acyl carrier protein phosphodiesterase	gene_strand=<	genes_inactivated=acpD	locus_tag=ECB_01367	locus_tags_inactivated=ECB_01367	mutation_category=mobile_element_insertion	position_end=1451974	position_start=1451972	ref_seq=GAG	repeat_size=1443
DEL	15	60	REL606	1607917	1	gene_name=ECB_01510	gene_position=coding (16/1320 nt)	gene_product=putative tail component of prophage	gene_strand=<	genes_inactivated=ECB_01510	locus_tag=ECB_01510	locus_tags_inactivated=ECB_01510	mutation_category=small_indel	position_end=1607917	position_start=1607917	ref_seq=T
SNP	16	61	REL606	1733343	A	aa_new_seq=N	aa_position=127	aa_ref_seq=D	codon_new_seq=AAT	codon_number=127	codon_position=1	codon_ref_seq=GAT	gene_name=pykF	gene_position=379	gene_product=pyruvate kinase	gene_strand=>	genes_overlapping=pykF	locus_tag=ECB_01645	locus_tags_overlapping=ECB_01645	mutation_category=snp_nonsynonymous	position_end=1733343	position_start=1733343	ref_seq=G	snp_type=nonsynonymous	transl_table=11
DEL	17	102,181	REL606	2032711	23293	between=manB,cpsG	gene_name=[manB]–[cpsG]	gene_product=[manB],manC,insB-14,insA-14,wbbD,wbbC,wzy,wbbB,wbbA,vioB,vioA,wzx,rmlC,rfbA,rfbD,rfbB,galF,wcaM,wcaL,wcaK,wzxC,wcaJ,[cpsG]	genes_inactivated=manB,manC,insB-14,insA-14,wbbD,wbbC,wzy,wbbB,wbbA,vioB,vioA,wzx,rmlC,rfbA,rfbD,rfbB,galF,wcaM,wcaL,wcaK,wzxC,wcaJ,cpsG	locus_tag=[ECB_01932]–[ECB_01954]	locus_tags_inactivated=ECB_01932,ECB_01933,ECB_01934,ECB_01935,ECB_01936,ECB_01937,ECB_01938,ECB_01939,ECB_01940,ECB_01941,ECB_01942,ECB_01943,ECB_01944,ECB_01945,ECB_01946,ECB_01947,ECB_01948,ECB_01949,ECB_01950,ECB_01951,ECB_01952,ECB_01953,ECB_01954	mutation_category=large_deletion	position_end=2056003	position_start=2032711	ref_seq=23293-bp
MOB	18	167,168,1925	REL606	2170662	IS150	-1	3	gene_name=yehZ	gene_position=coding (15-17/918 nt)	gene_product=predicted transporter subunit: periplasmic-binding component of ABC superfamily	gene_strand=<	genes_inactivated=yehZ	locus_tag=ECB_02061	locus_tags_inactivated=ECB_02061	mutation_category=mobile_element_insertion	position_end=2170664	position_start=2170662	ref_seq=AGC	repeat_size=1443
INS	19	63	REL606	2333539	T	gene_name=ECB_02200/nuoN	gene_position=intergenic (+35/+33)	gene_product=YadA protein/NADH dehydrogenase subunit N	gene_strand=>/<	locus_tag=ECB_02200/ECB_02201	mutation_category=small_indel	position_end=2333539	position_start=2333539	ref_seq=T
SNP	20	64	REL606	2389878	T	aa_new_seq=R	aa_position=33	aa_ref_seq=S	codon_new_seq=AGA	codon_number=33	codon_position=3	codon_ref_seq=AGC	gene_name=mepA	gene_position=99	gene_product=penicillin-insensitive murein endopeptidase	gene_strand=<	genes_overlapping=mepA	locus_tag=ECB_02253	locus_tags_overlapping=ECB_02253	mutation_category=snp_nonsynonymous	position_end=2389878	position_start=2389878	ref_seq=G	snp_type=nonsynonymous	transl_table=11
SNP	21	65	REL606	2446984	C	gene_name=mntH/nupC	gene_position=intergenic (-53/-283)	gene_product=manganese transport protein MntH/nucleoside (except guanosine) transporter	gene_strand=</>	genes_promoter=mntH	locus_tag=ECB_02301/ECB_02302	locus_tags_promoter=ECB_02301	mutation_category=snp_intergenic	position_end=2446984	position_start=2446984	ref_seq=A	snp_type=intergenic
SNP	22	66	REL606	2665639	T	aa_new_seq=Q	aa_position=100	aa_ref_seq=L	codon_new_seq=CAG	codon_number=100	codon_position=2	codon_ref_seq=CTG	gene_name=rplS	gene_position=299	gene_product=50S ribosomal protein L19	gene_strand=<	genes_overlapping=rplS	locus_tag=ECB_02495	locus_tags_overlapping=ECB_02495	mutation_category=snp_nonsynonymous	position_end=2665639	position_start=2665639	ref_seq=A	snp_type=nonsynonymous	transl_table=11
SNP	23	67	REL606	2843036	A	aa_new_seq=G	aa_position=198	aa_ref_seq=G	codon_new_seq=GGT	codon_number=198	codon_position=3	codon_ref_seq=GGC	gene_name=amiC	gene_position=594	gene_product=N-acetylmuramoyl-L-alanine amidase	gene_strand=<	genes_overlapping=amiC	locus_tag=ECB_02665	locus_tags_overlapping=ECB_02665	mutation_category=snp_synonymous	position_end=2843036	position_start=2843036	ref_seq=G	snp_type=synonymous	transl_table=11
SNP	24	68	REL606	2999330	A	gene_name=yeeP/flu	gene_position=intergenic (+138/-234)	gene_product=b1999; predicted GTP-binding protein/antigen 43 (Ag43) phase-variable biofilm formation autotransporter	gene_strand=>/>	locus_tag=ECB_02799/ECB_02800	mutation_category=snp_intergenic	position_end=2999330	position_start=2999330	ref_seq=G	snp_type=intergenic
DEL	25	117,174	REL606	3260572	6	gene_name=hflB	gene_position=coding (1593-1598/1935 nt)	gene_product=protease, ATP-dependent zinc-metallo	gene_strand=<	genes_overlapping=hflB	locus_tag=ECB_03043	locus_tags_overlapping=ECB_03043	mutation_category=small_indel	position_end=3260577	position_start=3260572	ref_seq=TTCGCT	repeat_length=6	repeat_new_copies=1	repeat_ref_copies=2	repeat_seq=TTCGCT
SNP	26	69	REL606	3288025	A	aa_new_seq=L	aa_position=79	aa_ref_seq=Q	codon_new_seq=CTA	codon_number=79	codon_position=2	codon_ref_seq=CAA	gene_name=arcB	gene_position=236	gene_product=hybrid sensory histidine kinase in two-component regulatory system with ArcA	gene_strand=<	genes_overlapping=arcB	locus_tag=ECB_03075	locus_tags_overlapping=ECB_03075	mutation_category=snp_nonsynonymous	position_end=3288025	position_start=3288025	ref_seq=T	snp_type=nonsynonymous	transl_table=11
SNP	27	70	REL606	3339313	C	aa_new_seq=S	aa_position=51	aa_ref_seq=Y	codon_new_seq=TCT	codon_number=51	codon_position=2	codon_ref_seq=TAT	gene_name=fis	gene_position=152	gene_product=DNA-binding protein Fis	gene_strand=>	genes_overlapping=fis	locus_tag=ECB_03119	locus_tags_overlapping=ECB_03119	mutation_category=snp_nonsynonymous	position_end=3339313	position_start=3339313	ref_seq=A	snp_type=nonsynonymous	transl_table=11
SNP	28	71	REL606	3401754	A	aa_new_seq=L	aa_position=78	aa_ref_seq=R	codon_new_seq=CTC	codon_number=78	codon_position=2	codon_ref_seq=CGC	gene_name=rpsG	gene_position=233	gene_product=30S ribosomal protein S7	gene_strand=<	genes_overlapping=rpsG	locus_tag=ECB_03192	locus_tags_overlapping=ECB_03192	mutation_category=snp_nonsynonymous	position_end=3401754	position_start=3401754	ref_seq=C	snp_type=nonsynonymous	transl_table=11
SNP	29	72	REL606	3463900	T	aa_new_seq=H	aa_position=218	aa_ref_seq=R	codon_new_seq=CAC	codon_number=218	codon_position=2	codon_ref_seq=CGC	gene_name=envZ	gene_position=653	gene_product=osmolarity sensor protein	gene_strand=<	genes_overlapping=envZ	locus_tag=ECB_03256	locus_tags_overlapping=ECB_03256	mutation_category=snp_nonsynonymous	position_end=3463900	position_start=3463900	ref_seq=C	snp_type=nonsynonymous	transl_table=11
SNP	30	73	REL606	3481820	G	aa_new_seq=A	aa_position=46	aa_ref_seq=T	codon_new_seq=GCC	codon_number=46	codon_position=1	codon_ref_seq=ACC	gene_name=malT	gene_position=136	gene_product=transcriptional regulator MalT	gene_strand=>	genes_overlapping=malT	locus_tag=ECB_03270	locus_tags_overlapping=ECB_03270	mutation_category=snp_nonsynonymous	position_end=3481820	position_start=3481820	ref_seq=A	snp_type=nonsynonymous	transl_table=11
SNP	31	74	REL606	3488669	C	aa_new_seq=A	aa_position=180	aa_ref_seq=S	codon_new_seq=GCG	codon_number=180	codon_position=1	codon_ref_seq=TCG	gene_name=glpR	gene_position=538	gene_product=DNA-binding transcriptional repressor	gene_strand=<	genes_overlapping=glpR	locus_tag=ECB_03274	locus_tags_overlapping=ECB_03274	mutation_category=snp_nonsynonymous	position_end=3488669	position_start=3488669	ref_seq=A	snp_type=nonsynonymous	transl_table=11
SNP	32	75	REL606	3613023	C	gene_name=dctA/yhjK	gene_position=intergenic (-84/+99)	gene_product=C4-dicarboxylate transport protein/predicted diguanylate cyclase	gene_strand=</<	genes_promoter=dctA	locus_tag=ECB_03376/ECB_03377	locus_tags_promoter=ECB_03376	mutation_category=snp_intergenic	position_end=3613023	position_start=3613023	ref_seq=A	snp_type=intergenic
DEL	33	76,77	REL606	3893549	2	gene_name=kup/insJ-5	gene_position=intergenic (+4/-51)	gene_product=potassium transporter/IS150 hypothetical protein	gene_strand=>/>	genes_promoter=insJ-5	locus_tag=ECB_03631/ECB_03632	locus_tags_promoter=ECB_03632	mutation_category=small_indel	position_end=3893550	position_start=3893549	ref_seq=CA
SNP	34	78	REL606	3909807	T	gene_name=hdfR/yifE	gene_position=intergenic (-40/-79)	gene_product=transcriptional regulator HdfR/hypothetical protein	gene_strand=</>	genes_promoter=hdfR,yifE	locus_tag=ECB_03642/ECB_03643	locus_tags_promoter=ECB_03642,ECB_03643	mutation_category=snp_intergenic	position_end=3909807	position_start=3909807	ref_seq=G	snp_type=intergenic
DEL	35	79	REL606	3911972	1	gene_name=yifB/ilvL	gene_position=intergenic (-173/-150)	gene_product=predicted bifunctional enzyme and transcriptional regulator/ilvG operon leader peptide	gene_strand=</>	genes_promoter=ilvL	locus_tag=ECB_03644/ECB_03645	locus_tags_promoter=ECB_03645	mutation_category=small_indel	position_end=3911972	position_start=3911972	ref_seq=T	repeat_length=1	repeat_new_copies=6	repeat_ref_copies=7	repeat_seq=T
SNP	36	80	REL606	4100183	G	aa_new_seq=P	aa_position=350	aa_ref_seq=S	codon_new_seq=CCT	codon_number=350	codon_position=1	codon_ref_seq=TCT	gene_name=hslU	gene_position=1048	gene_product=ATP-dependent protease ATP-binding subunit	gene_strand=<	genes_overlapping=hslU	locus_tag=ECB_03816	locus_tags_overlapping=ECB_03816	mutation_category=snp_nonsynonymous	position_end=4100183	position_start=4100183	ref_seq=A	snp_type=nonsynonymous	transl_table=11
SNP	37	81	REL606	4107018	A	aa_new_seq=P	aa_position=123	aa_ref_seq=P	codon_new_seq=CCT	codon_number=123	codon_position=3	codon_ref_seq=CCA	gene_name=ECB_03822	gene_position=369	gene_product=putative outer membrane lipoprotein	gene_strand=<	genes_overlapping=ECB_03822	locus_tag=ECB_03822	locus_tags_overlapping=ECB_03822	mutation_category=snp_synonymous	position_end=4107018	position_start=4107018	ref_seq=T	snp_type=synonymous	transl_table=11
SNP	38	82	REL606	4180200	G	aa_new_seq=G	aa_position=8	aa_ref_seq=V	codon_new_seq=GGG	codon_number=8	codon_position=2	codon_ref_seq=GTG	gene_name=yjaH	gene_position=23	gene_product=hypothetical protein	gene_strand=>	genes_overlapping=yjaH	locus_tag=ECB_03878	locus_tags_overlapping=ECB_03878	mutation_category=snp_nonsynonymous	position_end=4180200	position_start=4180200	ref_seq=T	snp_type=nonsynonymous	transl_table=11
SNP	39	83	REL606	4201958	C	aa_new_seq=R	aa_position=201	aa_ref_seq=L	codon_new_seq=CGC	codon_number=201	codon_position=2	codon_ref_seq=CTC	gene_name=iclR	gene_position=602	gene_product=DNA-binding transcriptional repressor	gene_strand=<	genes_overlapping=iclR	locus_tag=ECB_03890	locus_tags_overlapping=ECB_03890	mutation_category=snp_nonsynonymous	position_end=4201958	position_start=4201958	ref_seq=A	snp_type=nonsynonymous	transl_table=11
SNP	40	84	REL606	4363338	A	aa_new_seq=G	aa_position=263	aa_ref_seq=G	codon_new_seq=GGA	codon_number=263	codon_position=3	codon_ref_seq=GGC	gene_name=yjeM	gene_position=789	gene_product=predicted transporter	gene_strand=>	genes_overlapping=yjeM	locus_tag=ECB_04028	locus_tags_overlapping=ECB_04028	mutation_category=snp_synonymous	position_end=4363338	position_start=4363338	ref_seq=C	snp_type=synonymous	transl_table=11
DEL	41	85,86	REL606	4431394	2	gene_name=ytfN	gene_position=coding (3526-3527/3780 nt)	gene_product=hypothetical protein	gene_strand=>	genes_overlapping=ytfN	locus_tag=ECB_04092	locus_tags_overlapping=ECB_04092	mutation_category=small_indel	position_end=4431395	position_start=4431394	ref_seq=GG
SNP	42	87	REL606	4433347	G	aa_new_seq=E	aa_position=46	aa_ref_seq=K	codon_new_seq=GAA	codon_number=46	codon_position=1	codon_ref_seq=AAA	gene_name=ytfQ	gene_position=136	gene_product=predicted sugar transporter subunit: periplasmic-binding component of ABC superfamily	gene_strand=>	genes_overlapping=ytfQ	locus_tag=ECB_04096	locus_tags_overlapping=ECB_04096	mutation_category=snp_nonsynonymous	position_end=4433347	position_start=4433347	ref_seq=A	snp_type=nonsynonymous	transl_table=11
DEL	43	142,160,186	REL606	4522339	38944	gene_name=[fimB]–[hsdR]	gene_product=[fimB],fimE,fimA,fimI,fimC,fimD,fimF,fimG,fimH,gntP,uxuA,uxuB,uxuR,yjiC,yjiD,yjiE,iadA,yjiG,yjiH,kptA,yjiJ,yjiK,yjiL,yjiM,yjiN,yjiO,yjiPQ,insA-28,insB-28,yjiV,mcrC,mcrB,yjiW,hsdS,hsdM,[hsdR]	genes_inactivated=fimB,fimE,fimA,fimI,fimC,fimD,fimF,fimG,fimH,gntP,uxuA,uxuB,uxuR,yjiC,yjiD,yjiE,iadA,yjiG,yjiH,kptA,yjiJ,yjiK,yjiL,yjiM,yjiN,yjiO,yjiPQ,insA-28,insB-28,yjiV,mcrC,mcrB,yjiW,hsdS,hsdM	genes_overlapping=hsdR	locus_tag=[ECB_04179]–[ECB_04216]	locus_tags_inactivated=ECB_04179,ECB_04182,ECB_04183,ECB_04184,ECB_04185,ECB_04186,ECB_04187,ECB_04188,ECB_04189,ECB_04190,ECB_04191,ECB_04192,ECB_04193,ECB_04194,ECB_04195,ECB_04196,ECB_04197,ECB_04198,ECB_04199,ECB_04200,ECB_04201,ECB_04202,ECB_04203,ECB_04204,ECB_04205,ECB_04206,ECB_04207,ECB_04208,ECB_04209,ECB_04210,ECB_04211,ECB_04212,ECB_04213,ECB_04214,ECB_04215	locus_tags_overlapping=ECB_04216	mediated=IS1	mutation_category=large_deletion	position_end=4561282	position_start=4522339	ref_seq=38944-bp
MOB	44	148,157,1937	REL606	4613019	IS150	1	0	gene_name=smp/serB	gene_position=intergenic (-16/-90)	gene_product=hypothetical protein/3-phosphoserine phosphatase	gene_strand=</>	genes_promoter=smp,serB	ins_start=TT	locus_tag=ECB_04263/ECB_04264	locus_tags_promoter=ECB_04263,ECB_04264	mutation_category=mobile_element_insertion	position_end=4613019	position_start=4613019	ref_seq=A	repeat_size=1443
SNP	45	88	REL606	4616538	C	aa_new_seq=S	aa_position=337	aa_ref_seq=Y	codon_new_seq=TCC	codon_number=337	codon_position=2	codon_ref_seq=TAC	gene_name=nadR	gene_position=1010	gene_product=nicotinamide-nucleotide adenylyltransferase	gene_strand=>	genes_overlapping=nadR	locus_tag=ECB_04266	locus_tags_overlapping=ECB_04266	mutation_category=snp_nonsynonymous	position_end=4616538	position_start=4616538	ref_seq=A	snp_type=nonsynonymous	transl_table=11
RA	46	.	REL606	9972	0	T	G	aa_new_seq=T	aa_position=174	aa_ref_seq=N	allele_frequencies=G:1.000e+00	codon_new_seq=ACC	codon_number=174	codon_position=2	codon_ref_seq=AAC	fisher_strand_p_value=1.00000e+00	frequency=1.000e+00	frequency_lower=9.573e-01	frequency_upper=1.000e+00	gene_name=yaaH	gene_position=521	gene_product=conserved inner membrane protein associated with acetate transport	gene_strand=<	ks_quality_p_value=1.00000e+00	locus_tag=ECB_00010	major_base=G	major_cov=26/5	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=26/5	new_seq=G	prediction=consensus	ref_cov=0/0	ref_seq=T	score=126.0	snp_type=nonsynonymous	total_cov=26/5	transl_table=11
RA	47	.	REL606	81158	0	A	C	aa_new_seq=A	aa_position=245	aa_ref_seq=E	allele_frequencies=C:1.000e+00	codon_new_seq=GCG	codon_number=245	codon_position=2	codon_ref_seq=GAG	fisher_strand_p_value=1.00000e+00	frequency=1.000e+00	frequency_lower=9.757e-01	frequency_upper=1.000e+00	gene_name=setA	gene_position=734	gene_product=broad specificity sugar efflux system	gene_strand=>	ks_quality_p_value=1.00000e+00	locus_tag=ECB_00072	major_base=C	major_cov=24/31	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=24/31	new_seq=C	prediction=consensus	ref_cov=0/0	ref_seq=A	score=227.1	snp_type=nonsynonymous	total_cov=24/31	transl_table=11
RA	48	.	REL606	216480	0	C	T	allele_frequencies=T:1.000e+00	fisher_strand_p_value=1.00000e+00	frequency=1.000e+00	frequency_lower=9.621e-01	frequency_upper=1.000e+00	gene_name=tilS/rof	gene_position=intergenic (+9/+40)	gene_product=tRNA(Ile)-lysidine synthetase/modulator of Rho-dependent transcription termination	gene_strand=>/<	ks_quality_p_value=1.00000e+00	locus_tag=ECB_00186/ECB_00187	major_base=T	major_cov=4/31	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=4/31	prediction=consensus	ref_cov=0/0	score=135.8	snp_type=intergenic	total_cov=4/31
RA	49	.	REL606	247796	0	T	C	aa_new_seq=G	aa_position=27	aa_ref_seq=S	allele_frequencies=C:1.000e+00	codon_new_seq=GGT	codon_number=27	codon_position=1	codon_ref_seq=AGT	fisher_strand_p_value=1.00000e+00	frequency=1.000e+00	frequency_lower=9.716e-01	frequency_upper=1.000e+00	gene_name=fadE	gene_position=79	gene_product=acyl-CoA dehydrogenase	gene_strand=<	ks_quality_p_value=1.00000e+00	locus_tag=ECB_00216	major_base=C	major_cov=12/35	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=12/35	new_seq=C	prediction=consensus	ref_cov=0/0	ref_seq=T	score=202.2	snp_type=nonsynonymous	total_cov=12/35	transl_table=11
RA	50	.	REL606	281923	0	G	T	allele_frequencies=T:1.000e+00	fisher_strand_p_value=1.00000e+00	frequency=1.000e+00	frequency_lower=9.641e-01	frequency_upper=1.000e+00	gene_name=ykgK/ykgL	gene_position=intergenic (-406/-369)	gene_product=predicted regulator/hypothetical protein	gene_strand=</>	ks_quality_p_value=1.00000e+00	locus_tag=ECB_00253/ECB_00254	major_base=T	major_cov=0/37	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=0/37	prediction=consensus	ref_cov=0/0	score=150.2	snp_type=intergenic	total_cov=0/37
RA	51	.	REL606	398061	0	G	C	aa_new_seq=A	aa_position=45	aa_ref_seq=G	allele_frequencies=C:1.000e+00	codon_new_seq=GCC	codon_number=45	codon_position=2	codon_ref_seq=GGC	fisher_strand_p_value=1.00000e+00	frequency=1.000e+00	frequency_lower=9.544e-01	frequency_upper=1.000e+00	gene_name=secF	gene_position=134	gene_product=protein export protein SecF	gene_strand=>	ks_quality_p_value=1.00000e+00	locus_tag=ECB_00357	major_base=C	major_cov=17/12	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=17/12	new_seq=C	prediction=consensus	ref_cov=0/0	ref_seq=G	score=114.1	snp_type=nonsynonymous	total_cov=17/12	transl_table=11
RA	52	.	REL606	433366	1	.	T	allele_frequencies=T:1.000e+00	fisher_strand_p_value=1.00000e+00	frequency=1.000e+00	frequency_lower=9.429e-01	frequency_upper=1.000e+00	gene_name=lon/hupB	gene_position=intergenic (+69/-140)	gene_product=DNA-binding ATP-dependent protease La/HU, DNA-binding transcriptional regulator, beta subunit	gene_strand=>/>	ks_quality_p_value=1.00000e+00	locus_tag=ECB_00391/ECB_00392	major_base=T	major_cov=15/8	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=15/8	prediction=consensus	ref_cov=0/0	score=77.2	snp_type=intergenic	total_cov=15/8
RA	53	.	REL606	473901	1	.	C	allele_frequencies=C:9.846e-01,.:1.538e-02	fisher_strand_p_value=2.03390e-01	frequency=9.846e-01	frequency_lower=9.410e-01	frequency_upper=9.998e-01	gene_name=ybaL	gene_position=coding (1405/1677 nt)	gene_product=predicted transporter with NAD(P)-binding Rossmann-fold domain	gene_strand=<	ks_quality_p_value=9.83051e-01	locus_tag=ECB_00429	major_base=C	major_cov=47/11	major_frequency=9.846e-01	minor_base=.	minor_cov=0/1	new_cov=47/11	prediction=consensus	ref_cov=0/1	score=152.2	total_cov=47/12
RA	54	.	REL606	473901	2	.	G	allele_frequencies=G:9.847e-01,.:1.526e-02	fisher_strand_p_value=2.03390e-01	frequency=9.847e-01	frequency_lower=9.411e-01	frequency_upper=9.999e-01	gene_name=ybaL	gene_position=coding (1405/1677 nt)	gene_product=predicted transporter with NAD(P)-binding Rossmann-fold domain	gene_strand=<	ks_quality_p_value=1.00000e+00	locus_tag=ECB_00429	major_base=G	major_cov=47/11	major_frequency=9.847e-01	minor_base=.	minor_cov=0/1	new_cov=47/11	prediction=consensus	ref_cov=0/1	score=153.9	total_cov=47/12
RA	55	.	REL606	648692	0	C	T	aa_new_seq=H	aa_position=69	aa_ref_seq=R	allele_frequencies=T:1.000e+00	codon_new_seq=CAC	codon_number=69	codon_position=2	codon_ref_seq=CGC	fisher_strand_p_value=1.00000e+00	frequency=1.000e+00	frequency_lower=9.404e-01	frequency_upper=1.000e+00	gene_name=mrdB	gene_position=206	gene_product=cell wall shape-determining protein	gene_strand=<	ks_quality_p_value=1.00000e+00	locus_tag=ECB_00603	major_base=T	major_cov=13/9	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=13/9	new_seq=T	prediction=consensus	ref_cov=0/0	ref_seq=C	score=87.4	snp_type=nonsynonymous	total_cov=13/9	transl_table=11
RA	56	.	REL606	734488	0	C	T	aa_new_seq=T	aa_position=258	aa_ref_seq=A	allele_frequencies=T:1.000e+00	codon_new_seq=ACT	codon_number=258	codon_position=1	codon_ref_seq=GCT	fisher_strand_p_value=1.00000e+00	frequency=1.000e+00	frequency_lower=9.452e-01	frequency_upper=1.000e+00	gene_name=gltA	gene_position=772	gene_product=citrate synthase	gene_strand=<	ks_quality_p_value=1.00000e+00	locus_tag=ECB_00680	major_base=T	major_cov=9/15	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=9/15	new_seq=T	prediction=consensus	ref_cov=0/0	ref_seq=C	score=95.5	snp_type=nonsynonymous	total_cov=9/15	transl_table=11
RA	57	.	REL606	943475	0	C	T	aa_new_seq=I	aa_position=9	aa_ref_seq=M	allele_frequencies=T:1.000e+00	codon_new_seq=ATA	codon_number=9	codon_position=3	codon_ref_seq=ATG	fisher_strand_p_value=1.00000e+00	frequency=1.000e+00	frequency_lower=9.854e-01	frequency_upper=1.000e+00	gene_name=infA	gene_position=27	gene_product=translation initiation factor IF-1	gene_strand=<	ks_quality_p_value=1.00000e+00	locus_tag=ECB_00888	major_base=T	major_cov=64/28	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=64/28	new_seq=T	prediction=consensus	ref_cov=0/0	ref_seq=C	score=380.3	snp_type=nonsynonymous	total_cov=64/28	transl_table=11
RA	58	.	REL606	1331794	0	C	A	aa_new_seq=K	aa_position=792	aa_ref_seq=T	allele_frequencies=A:1.000e+00	codon_new_seq=AAG	codon_number=792	codon_position=2	codon_ref_seq=ACG	fisher_strand_p_value=1.00000e+00	frequency=1.000e+00	frequency_lower=9.733e-01	frequency_upper=1.000e+00	gene_name=topA	gene_position=2375	gene_product=DNA topoisomerase I	gene_strand=>	ks_quality_p_value=1.00000e+00	locus_tag=ECB_01250	major_base=A	major_cov=44/6	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=44/6	new_seq=A	prediction=consensus	ref_cov=0/0	ref_seq=C	score=203.5	snp_type=nonsynonymous	total_cov=44/6	transl_table=11
RA	59	.	REL606	1383157	0	T	G	aa_new_seq=G	aa_position=224	aa_ref_seq=G	allele_frequencies=G:7.148e-01,T:2.852e-01	codon_new_seq=GGG	codon_number=224	codon_position=3	codon_ref_seq=GGT	fisher_strand_p_value=9.38899e-04	frequency=7.148e-01	frequency_lower=5.975e-01	frequency_upper=8.151e-01	gene_name=ycjX	gene_position=672	gene_product=conserved protein with nucleoside triphosphate hydrolase domain	gene_strand=>	ks_quality_p_value=4.93189e-08	locus_tag=ECB_01298	major_base=G	major_cov=33/0	major_frequency=7.148e-01	minor_base=T	minor_cov=8/5	new_cov=33/0	new_seq=G	prediction=consensus	ref_cov=8/5	ref_seq=T	score=70.8	snp_type=synonymous	total_cov=41/5	transl_table=11
RA	60	.	REL606	1607917	0	T	.	allele_frequencies=.:1.000e+00	fisher_strand_p_value=1.00000e+00	frequency=1.000e+00	frequency_lower=9.667e-01	frequency_upper=1.000e+00	gene_name=ECB_01510	gene_position=coding (16/1320 nt)	gene_product=putative tail component of prophage	gene_strand=<	ks_quality_p_value=1.00000e+00	locus_tag=ECB_01510	major_base=.	major_cov=32/8	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=32/8	prediction=consensus	ref_cov=0/0	score=148.1	total_cov=32/8
RA	61	.	REL606	1733343	0	G	A	aa_new_seq=N	aa_position=127	aa_ref_seq=D	allele_frequencies=A:1.000e+00	codon_new_seq=AAT	codon_number=127	codon_position=1	codon_ref_seq=GAT	fisher_strand_p_value=1.00000e+00	frequency=1.000e+00	frequency_lower=9.621e-01	frequency_upper=1.000e+00	gene_name=pykF	gene_position=379	gene_product=pyruvate kinase	gene_strand=>	ks_quality_p_value=1.00000e+00	locus_tag=ECB_01645	major_base=A	major_cov=17/18	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=17/18	new_seq=A	prediction=consensus	ref_cov=0/0	ref_seq=G	score=135.6	snp_type=nonsynonymous	total_cov=17/18	transl_table=11
RA	62	.	REL606	2054830	0	T	A	allele_frequencies=A:1.000e+00	deleted=1	fisher_strand_p_value=1.00000e+00	frequency=1.000e+00	frequency_lower=7.981e-01	frequency_upper=1.000e+00	ks_quality_p_value=1.00000e+00	major_base=A	major_cov=0/12	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=0/12	prediction=consensus	ref_cov=0/0	score=18.3	total_cov=0/12
RA	63	.	REL606	2333539	1	.	T	allele_frequencies=T:1.000e+00	fisher_strand_p_value=1.00000e+00	frequency=1.000e+00	frequency_lower=9.586e-01	frequency_upper=1.000e+00	gene_name=ECB_02200/nuoN	gene_position=intergenic (+35/+33)	gene_product=YadA protein/NADH dehydrogenase subunit N	gene_strand=>/<	ks_quality_p_value=1.00000e+00	locus_tag=ECB_02200/ECB_02201	major_base=T	major_cov=17/15	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=17/15	prediction=consensus	ref_cov=0/0	score=101.2	snp_type=intergenic	total_cov=17/15
RA	64	.	REL606	2389878	0	G	T	aa_new_seq=R	aa_position=33	aa_ref_seq=S	allele_frequencies=T:1.000e+00	codon_new_seq=AGA	codon_number=33	codon_position=3	codon_ref_seq=AGC	fisher_strand_p_value=1.00000e+00	frequency=1.000e+00	frequency_lower=9.559e-01	frequency_upper=1.000e+00	gene_name=mepA	gene_position=99	gene_product=penicillin-insensitive murein endopeptidase	gene_strand=<	ks_quality_p_value=1.00000e+00	locus_tag=ECB_02253	major_base=T	major_cov=21/9	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=21/9	new_seq=T	prediction=consensus	ref_cov=0/0	ref_seq=G	score=119.8	snp_type=nonsynonymous	total_cov=21/9	transl_table=11
RA	65	.	REL606	2446984	0	A	C	allele_frequencies=C:1.000e+00	fisher_strand_p_value=1.00000e+00	frequency=1.000e+00	frequency_lower=9.838e-01	frequency_upper=1.000e+00	gene_name=mntH/nupC	gene_position=intergenic (-53/-283)	gene_product=manganese transport protein MntH/nucleoside (except guanosine) transporter	gene_strand=</>	ks_quality_p_value=1.00000e+00	locus_tag=ECB_02301/ECB_02302	major_base=C	major_cov=51/32	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=51/32	prediction=consensus	ref_cov=0/0	score=350.2	snp_type=intergenic	total_cov=51/32
RA	66	.	REL606	2665639	0	A	T	aa_new_seq=Q	aa_position=100	aa_ref_seq=L	allele_frequencies=A:2.000e-02,T:9.800e-01	codon_new_seq=CAG	codon_number=100	codon_position=2	codon_ref_seq=CTG	fisher_strand_p_value=1.00000e+00	frequency=9.800e-01	frequency_lower=9.290e-01	frequency_upper=9.979e-01	gene_name=rplS	gene_position=299	gene_product=50S ribosomal protein L19	gene_strand=<	ks_quality_p_value=8.40000e-01	locus_tag=ECB_02495	major_base=T	major_cov=41/8	major_frequency=9.800e-01	minor_base=A	minor_cov=1/0	new_cov=41/8	new_seq=T	prediction=consensus	ref_cov=1/0	ref_seq=A	score=191.5	snp_type=nonsynonymous	total_cov=42/8	transl_table=11
RA	67	.	REL606	2843036	0	G	A	aa_new_seq=G	aa_position=198	aa_ref_seq=G	allele_frequencies=A:1.000e+00	codon_new_seq=GGT	codon_number=198	codon_position=3	codon_ref_seq=GGC	fisher_strand_p_value=1.00000e+00	frequency=1.000e+00	frequency_lower=9.823e-01	frequency_upper=1.000e+00	gene_name=amiC	gene_position=594	gene_product=N-acetylmuramoyl-L-alanine amidase	gene_strand=<	ks_quality_p_value=1.00000e+00	locus_tag=ECB_02665	major_base=A	major_cov=18/67	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=18/67	new_seq=A	prediction=consensus	ref_cov=0/1	ref_seq=G	score=323.6	snp_type=synonymous	total_cov=18/68	transl_table=11
RA	68	.	REL606	2999330	0	G	A	allele_frequencies=A:1.000e+00	fisher_strand_p_value=1.00000e+00	frequency=1.000e+00	frequency_lower=7.981e-01	frequency_upper=1.000e+00	gene_name=yeeP/flu	gene_position=intergenic (+138/-234)	gene_product=b1999; predicted GTP-binding protein/antigen 43 (Ag43) phase-variable biofilm formation autotransporter	gene_strand=>/>	ks_quality_p_value=1.00000e+00	locus_tag=ECB_02799/ECB_02800	major_base=A	major_cov=6/0	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=6/0	prediction=consensus	ref_cov=0/0	score=19.1	snp_type=intergenic	total_cov=6/0
RA	69	.	REL606	3288025	0	T	A	aa_new_seq=L	aa_position=79	aa_ref_seq=Q	allele_frequencies=A:1.000e+00	codon_new_seq=CTA	codon_number=79	codon_position=2	codon_ref_seq=CAA	fisher_strand_p_value=1.00000e+00	frequency=1.000e+00	frequency_lower=9.728e-01	frequency_upper=1.000e+00	gene_name=arcB	gene_position=236	gene_product=hybrid sensory histidine kinase in two-component regulatory system with ArcA	gene_strand=<	ks_quality_p_value=1.00000e+00	locus_tag=ECB_03075	major_base=A	major_cov=37/12	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=37/12	new_seq=A	prediction=consensus	ref_cov=0/0	ref_seq=T	score=197.5	snp_type=nonsynonymous	total_cov=37/12	transl_table=11
RA	70	.	REL606	3339313	0	A	C	aa_new_seq=S	aa_position=51	aa_ref_seq=Y	allele_frequencies=C:1.000e+00	codon_new_seq=TCT	codon_number=51	codon_position=2	codon_ref_seq=TAT	fisher_strand_p_value=1.00000e+00	frequency=1.000e+00	frequency_lower=9.743e-01	frequency_upper=1.000e+00	gene_name=fis	gene_position=152	gene_product=DNA-binding protein Fis	gene_strand=>	ks_quality_p_value=1.00000e+00	locus_tag=ECB_03119	major_base=C	major_cov=22/30	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=22/30	new_seq=C	prediction=consensus	ref_cov=0/0	ref_seq=A	score=220.3	snp_type=nonsynonymous	total_cov=22/30	transl_table=11
RA	71	.	REL606	3401754	0	C	A	aa_new_seq=L	aa_position=78	aa_ref_seq=R	allele_frequencies=A:1.000e+00	codon_new_seq=CTC	codon_number=78	codon_position=2	codon_ref_seq=CGC	fisher_strand_p_value=1.00000e+00	frequency=1.000e+00	frequency_lower=9.753e-01	frequency_upper=1.000e+00	gene_name=rpsG	gene_position=233	gene_product=30S ribosomal protein S7	gene_strand=<	ks_quality_p_value=1.00000e+00	locus_tag=ECB_03192	major_base=A	major_cov=15/39	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=15/39	new_seq=A	prediction=consensus	ref_cov=0/0	ref_seq=C	score=220.9	snp_type=nonsynonymous	total_cov=15/39	transl_table=11
RA	72	.	REL606	3463900	0	C	T	aa_new_seq=H	aa_position=218	aa_ref_seq=R	allele_frequencies=T:1.000e+00	codon_new_seq=CAC	codon_number=218	codon_position=2	codon_ref_seq=CGC	fisher_strand_p_value=1.00000e+00	frequency=1.000e+00	frequency_lower=9.710e-01	frequency_upper=1.000e+00	gene_name=envZ	gene_position=653	gene_product=osmolarity sensor protein	gene_strand=<	ks_quality_p_value=1.00000e+00	locus_tag=ECB_03256	major_base=T	major_cov=13/33	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=13/33	new_seq=T	prediction=consensus	ref_cov=0/0	ref_seq=C	score=188.4	snp_type=nonsynonymous	total_cov=13/33	transl_table=11
RA	73	.	REL606	3481820	0	A	G	aa_new_seq=A	aa_position=46	aa_ref_seq=T	allele_frequencies=G:1.000e+00	codon_new_seq=GCC	codon_number=46	codon_position=1	codon_ref_seq=ACC	fisher_strand_p_value=1.00000e+00	frequency=1.000e+00	frequency_lower=9.346e-01	frequency_upper=1.000e+00	gene_name=malT	gene_position=136	gene_product=transcriptional regulator MalT	gene_strand=>	ks_quality_p_value=1.00000e+00	locus_tag=ECB_03270	major_base=G	major_cov=9/11	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=9/11	new_seq=G	prediction=consensus	ref_cov=0/0	ref_seq=A	score=81.5	snp_type=nonsynonymous	total_cov=9/11	transl_table=11
RA	74	.	REL606	3488669	0	A	C	aa_new_seq=A	aa_position=180	aa_ref_seq=S	allele_frequencies=C:1.000e+00	codon_new_seq=GCG	codon_number=180	codon_position=1	codon_ref_seq=TCG	fisher_strand_p_value=1.00000e+00	frequency=1.000e+00	frequency_lower=8.604e-01	frequency_upper=1.000e+00	gene_name=glpR	gene_position=538	gene_product=DNA-binding transcriptional repressor	gene_strand=<	ks_quality_p_value=1.00000e+00	locus_tag=ECB_03274	major_base=C	major_cov=1/8	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=1/8	new_seq=C	prediction=consensus	ref_cov=0/0	ref_seq=A	score=33.4	snp_type=nonsynonymous	total_cov=1/8	transl_table=11
RA	75	.	REL606	3613023	0	A	C	allele_frequencies=C:1.000e+00	fisher_strand_p_value=1.00000e+00	frequency=1.000e+00	frequency_lower=9.757e-01	frequency_upper=1.000e+00	gene_name=dctA/yhjK	gene_position=intergenic (-84/+99)	gene_product=C4-dicarboxylate transport protein/predicted diguanylate cyclase	gene_strand=</<	ks_quality_p_value=1.00000e+00	locus_tag=ECB_03376/ECB_03377	major_base=C	major_cov=2/53	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=2/53	prediction=consensus	ref_cov=0/0	score=232.4	snp_type=intergenic	total_cov=2/53
RA	76	.	REL606	3893549	0	C	.	allele_frequencies=.:1.000e+00	fisher_strand_p_value=1.00000e+00	frequency=1.000e+00	frequency_lower=9.450e-01	frequency_upper=1.000e+00	gene_name=kup/insJ-5	gene_position=intergenic (+4/-52)	gene_product=potassium transporter/IS150 hypothetical protein	gene_strand=>/>	ks_quality_p_value=1.00000e+00	locus_tag=ECB_03631/ECB_03632	major_base=.	major_cov=7/17	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=7/17	prediction=consensus	ref_cov=0/0	score=62.2	snp_type=intergenic	total_cov=7/17
RA	77	.	REL606	3893550	0	A	.	allele_frequencies=.:1.000e+00	fisher_strand_p_value=1.00000e+00	frequency=1.000e+00	frequency_lower=9.450e-01	frequency_upper=1.000e+00	gene_name=kup/insJ-5	gene_position=intergenic (+5/-51)	gene_product=potassium transporter/IS150 hypothetical protein	gene_strand=>/>	ks_quality_p_value=1.00000e+00	locus_tag=ECB_03631/ECB_03632	major_base=.	major_cov=7/17	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=7/17	prediction=consensus	ref_cov=0/0	score=62.2	snp_type=intergenic	total_cov=7/17
RA	78	.	REL606	3909807	0	G	T	allele_frequencies=T:1.000e+00	fisher_strand_p_value=1.00000e+00	frequency=1.000e+00	frequency_lower=9.886e-01	frequency_upper=1.000e+00	gene_name=hdfR/yifE	gene_position=intergenic (-40/-79)	gene_product=transcriptional regulator HdfR/hypothetical protein	gene_strand=</>	ks_quality_p_value=1.00000e+00	locus_tag=ECB_03642/ECB_03643	major_base=T	major_cov=68/50	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=68/50	prediction=consensus	ref_cov=0/0	score=464.6	snp_type=intergenic	total_cov=68/50
RA	79	.	REL606	3911972	0	T	.	allele_frequencies=.:1.000e+00	fisher_strand_p_value=1.00000e+00	frequency=1.000e+00	frequency_lower=9.733e-01	frequency_upper=1.000e+00	gene_name=yifB/ilvL	gene_position=intergenic (-173/-150)	gene_product=predicted bifunctional enzyme and transcriptional regulator/ilvG operon leader peptide	gene_strand=</>	ks_quality_p_value=1.00000e+00	locus_tag=ECB_03644/ECB_03645	major_base=.	major_cov=30/20	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=30/20	prediction=consensus	ref_cov=0/0	score=238.9	snp_type=intergenic	total_cov=30/20
RA	80	.	REL606	4100183	0	A	G	aa_new_seq=P	aa_position=350	aa_ref_seq=S	allele_frequencies=G:1.000e+00	codon_new_seq=CCT	codon_number=350	codon_position=1	codon_ref_seq=TCT	fisher_strand_p_value=1.00000e+00	frequency=1.000e+00	frequency_lower=9.803e-01	frequency_upper=1.000e+00	gene_name=hslU	gene_position=1048	gene_product=ATP-dependent protease ATP-binding subunit	gene_strand=<	ks_quality_p_value=1.00000e+00	locus_tag=ECB_03816	major_base=G	major_cov=29/39	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=29/39	new_seq=G	prediction=consensus	ref_cov=0/0	ref_seq=A	score=286.4	snp_type=nonsynonymous	total_cov=29/39	transl_table=11
RA	81	.	REL606	4107018	0	T	A	aa_new_seq=P	aa_position=123	aa_ref_seq=P	allele_frequencies=A:1.000e+00	codon_new_seq=CCT	codon_number=123	codon_position=3	codon_ref_seq=CCA	fisher_strand_p_value=1.00000e+00	frequency=1.000e+00	frequency_lower=9.832e-01	frequency_upper=1.000e+00	gene_name=ECB_03822	gene_position=369	gene_product=putative outer membrane lipoprotein	gene_strand=<	ks_quality_p_value=1.00000e+00	locus_tag=ECB_03822	major_base=A	major_cov=34/46	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=34/46	new_seq=A	prediction=consensus	ref_cov=0/0	ref_seq=T	score=326.5	snp_type=synonymous	total_cov=34/46	transl_table=11
RA	82	.	REL606	4180200	0	T	G	allele_frequencies=G:8.974e-01,.:1.026e-01	fisher_strand_p_value=6.01816e-03	frequency=8.974e-01	frequency_lower=7.999e-01	frequency_upper=9.591e-01	gene_name=yjaH	gene_position=coding (23/696 nt)	gene_product=hypothetical protein	gene_strand=>	ks_quality_p_value=8.02969e-01	locus_tag=ECB_03878	major_base=G	major_cov=27/8	major_frequency=8.974e-01	minor_base=.	minor_cov=0/4	new_cov=27/8	prediction=consensus	ref_cov=0/0	score=144.5	total_cov=27/12
RA	83	.	REL606	4201958	0	A	C	aa_new_seq=R	aa_position=201	aa_ref_seq=L	allele_frequencies=C:1.000e+00	codon_new_seq=CGC	codon_number=201	codon_position=2	codon_ref_seq=CTC	fisher_strand_p_value=1.00000e+00	frequency=1.000e+00	frequency_lower=9.800e-01	frequency_upper=1.000e+00	gene_name=iclR	gene_position=602	gene_product=DNA-binding transcriptional repressor	gene_strand=<	ks_quality_p_value=1.00000e+00	locus_tag=ECB_03890	major_base=C	major_cov=37/30	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=37/30	new_seq=C	prediction=consensus	ref_cov=0/0	ref_seq=A	score=282.8	snp_type=nonsynonymous	total_cov=37/30	transl_table=11
RA	84	.	REL606	4363338	0	C	A	aa_new_seq=G	aa_position=263	aa_ref_seq=G	allele_frequencies=A:1.000e+00	codon_new_seq=GGA	codon_number=263	codon_position=3	codon_ref_seq=GGC	fisher_strand_p_value=1.00000e+00	frequency=1.000e+00	frequency_lower=9.777e-01	frequency_upper=1.000e+00	gene_name=yjeM	gene_position=789	gene_product=predicted transporter	gene_strand=>	ks_quality_p_value=1.00000e+00	locus_tag=ECB_04028	major_base=A	major_cov=37/23	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=37/23	new_seq=A	prediction=consensus	ref_cov=0/0	ref_seq=C	score=241.8	snp_type=synonymous	total_cov=37/23	transl_table=11
RA	85	.	REL606	4431394	0	G	.	allele_frequencies=.:1.000e+00	fisher_strand_p_value=1.00000e+00	frequency=1.000e+00	frequency_lower=9.621e-01	frequency_upper=1.000e+00	gene_name=ytfN	gene_position=coding (3526/3780 nt)	gene_product=hypothetical protein	gene_strand=>	ks_quality_p_value=1.00000e+00	locus_tag=ECB_04092	major_base=.	major_cov=21/14	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=21/14	prediction=consensus	ref_cov=0/0	score=125.8	total_cov=21/14
RA	86	.	REL606	4431395	0	G	.	allele_frequencies=.:1.000e+00	fisher_strand_p_value=1.00000e+00	frequency=1.000e+00	frequency_lower=9.621e-01	frequency_upper=1.000e+00	gene_name=ytfN	gene_position=coding (3527/3780 nt)	gene_product=hypothetical protein	gene_strand=>	ks_quality_p_value=1.00000e+00	locus_tag=ECB_04092	major_base=.	major_cov=21/14	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=21/14	prediction=consensus	ref_cov=0/0	score=125.8	total_cov=21/14
RA	87	.	REL606	4433347	0	A	G	aa_new_seq=E	aa_position=46	aa_ref_seq=K	allele_frequencies=G:1.000e+00	codon_new_seq=GAA	codon_number=46	codon_position=1	codon_ref_seq=AAA	fisher_strand_p_value=1.00000e+00	frequency=1.000e+00	frequency_lower=9.757e-01	frequency_upper=1.000e+00	gene_name=ytfQ	gene_position=136	gene_product=predicted sugar transporter subunit: periplasmic-binding component of ABC superfamily	gene_strand=>	ks_quality_p_value=1.00000e+00	locus_tag=ECB_04096	major_base=G	major_cov=14/41	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=14/41	new_seq=G	prediction=consensus	ref_cov=0/0	ref_seq=A	score=232.5	snp_type=nonsynonymous	total_cov=14/41	transl_table=11
RA	88	.	REL606	4616538	0	A	C	aa_new_seq=S	aa_position=337	aa_ref_seq=Y	allele_frequencies=C:1.000e+00	codon_new_seq=TCC	codon_number=337	codon_position=2	codon_ref_seq=TAC	fisher_strand_p_value=1.00000e+00	frequency=1.000e+00	frequency_lower=9.710e-01	frequency_upper=1.000e+00	gene_name=nadR	gene_position=1010	gene_product=nicotinamide-nucleotide adenylyltransferase	gene_strand=>	ks_quality_p_value=1.00000e+00	locus_tag=ECB_04266	major_base=C	major_cov=23/23	major_frequency=1.000e+00	minor_base=N	minor_cov=0/0	new_cov=23/23	new_seq=C	prediction=consensus	ref_cov=0/0	ref_seq=A	score=194.5	snp_type=nonsynonymous	total_cov=23/23	transl_table=11
MC	89	.	REL606	167429	167438	5	0	gene_name=hrpB/mrcB	gene_position=intergenic (+44/-143)	gene_product=predicted ATP-dependent helicase/penicillin-binding protein 1b	gene_strand=>/>	left_inside_cov=2	left_outside_cov=6	locus_tag=ECB_00147/ECB_00148	right_inside_cov=0	right_outside_cov=6
MC	90	.	REL606	214460	214468	0	0	gene_name=ldcC	gene_position=coding (1940-1948/2142 nt)	gene_product=lysine decarboxylase 2, constitutive	gene_strand=>	left_inside_cov=3	left_outside_cov=16	locus_tag=ECB_00184	right_inside_cov=2	right_outside_cov=8
MC	91	.	REL606	328242	328255	0	0	gene_name=codB	gene_position=coding (765-778/1260 nt)	gene_product=cytosine transporter	gene_strand=>	left_inside_cov=2	left_outside_cov=6	locus_tag=ECB_00290	right_inside_cov=0	right_outside_cov=6
MC	92	.	REL606	547684	550350	16	0	gene_name=[nmpC]–[ybcU]	gene_product=[nmpC],ybcR,ybcS,ybcT,[ybcU]	left_inside_cov=4	left_outside_cov=6	locus_tag=[ECB_00503]–[ECB_00509]	right_inside_cov=0	right_outside_cov=10
MC	93	.	REL606	551250	551255	0	0	gene_name=ECB_00510/nohB	gene_position=intergenic (-333/-51)	gene_product=conserved hypothetical protein/DNA packaging protein	gene_strand=</>	left_inside_cov=1	left_outside_cov=17	locus_tag=ECB_00510/ECB_00511	right_inside_cov=1	right_outside_cov=9
MC	94	.	REL606	558072	558075	0	0	gene_name=ECB_00516	gene_position=coding (9-12/294 nt)	gene_product=conserved hypothetical protein	gene_strand=>	left_inside_cov=4	left_outside_cov=19	locus_tag=ECB_00516	right_inside_cov=0	right_outside_cov=5
MC	95	.	REL606	1029734	1029744	0	0	gene_name=pqiA	gene_position=coding (643-653/1254 nt)	gene_product=paraquat-inducible membrane protein A	gene_strand=>	left_inside_cov=3	left_outside_cov=15	locus_tag=ECB_00954	right_inside_cov=0	right_outside_cov=7
MC	96	.	REL606	1177332	1177340	0	0	gene_name=ycfM	gene_position=coding (126-134/642 nt)	gene_product=predicted outer membrane lipoprotein	gene_strand=>	left_inside_cov=4	left_outside_cov=12	locus_tag=ECB_01101	right_inside_cov=4	right_outside_cov=7
MC	97	.	REL606	1270055	1270662	100	84	gene_name=ldrC	gene_product=ldrC	left_inside_cov=4	left_outside_cov=28	locus_tag=[ECB_01193]	right_inside_cov=3	right_outside_cov=18
MC	98	.	REL606	1315434	1315434	0	0	gene_name=trpB	gene_position=coding (576/1194 nt)	gene_product=tryptophan synthase subunit beta	gene_strand=<	left_inside_cov=0	left_outside_cov=5	locus_tag=ECB_01235	right_inside_cov=0	right_outside_cov=7
MC	99	.	REL606	1462254	1462324	0	0	gene_name=mokB/trg	gene_position=intergenic (-16/-255)	gene_product=regulatory peptide/methyl-accepting chemotaxis protein III, ribose and galactose sensor receptor	gene_strand=</>	left_inside_cov=0	left_outside_cov=126	locus_tag=ECB_01377/ECB_01378	right_inside_cov=3	right_outside_cov=17
MC	100	.	REL606	1609173	1615472	3	1	gene_name=ybcW–ECB_01523	gene_product=ybcW,gnsB,ynfN,ECB_01516,cspI,ydfP,ydfQ,ydfR,essQ,ECB_01522,ECB_01523	left_inside_cov=2	left_outside_cov=6	locus_tag=[ECB_01513]–[ECB_01523]	right_inside_cov=4	right_outside_cov=7
MC	101	.	REL606	1881605	1881607	0	0	gene_name=manX/manY	gene_position=intergenic (+7/-54)	gene_product=fused mannose-specific PTS enzymes: IIA component/IIB component/mannose-specific enzyme IIC component of PTS	gene_strand=>/>	left_inside_cov=0	left_outside_cov=8	locus_tag=ECB_01787/ECB_01788	right_inside_cov=0	right_outside_cov=8
MC	102	.	REL606	2032710	2055364	0	0	gene_name=[manB]–[cpsG]	gene_product=[manB],manC,insB-14,insA-14,wbbD,wbbC,wzy,wbbB,wbbA,vioB,vioA,wzx,rmlC,rfbA,rfbD,rfbB,galF,wcaM,wcaL,wcaK,wzxC,wcaJ,[cpsG]	left_inside_cov=0	left_outside_cov=9	locus_tag=[ECB_01932]–[ECB_01954]	right_inside_cov=0	right_outside_cov=8
MC	103	.	REL606	2055399	2055606	0	32	gene_name=cpsG	gene_position=coding (491-698/1371 nt)	gene_product=phosphomannomutase	gene_strand=<	left_inside_cov=0	left_outside_cov=8	locus_tag=ECB_01954	right_inside_cov=4	right_outside_cov=13
MC	104	.	REL606	2055662	2055778	0	4	gene_name=cpsG	gene_position=coding (319-435/1371 nt)	gene_product=phosphomannomutase	gene_strand=<	left_inside_cov=0	left_outside_cov=11	locus_tag=ECB_01954	right_inside_cov=0	right_outside_cov=7
MC	105	.	REL606	2055826	2055902	0	12	gene_name=cpsG	gene_position=coding (195-271/1371 nt)	gene_product=phosphomannomutase	gene_strand=<	left_inside_cov=4	left_outside_cov=5	locus_tag=ECB_01954	right_inside_cov=4	right_outside_cov=11
MC	106	.	REL606	2086518	2122425	93	0	gene_name=yegM–[ECB_02013]	gene_product=yegM,yegN,yegO,yegB,baeS,baeR,yegP,yegQ,ogrK,yegZ,ECB_01989,ECB_01990,ECB_01991,ECB_01992,ECB_01993,ECB_01994,ECB_01995,ECB_01996,ECB_01997,ECB_01998,ECB_01999,ECB_02000,ECB_02001,ECB_02002,ECB_02003,ECB_02004,ECB_02005,ECB_02006,ECB_02007,ECB_02008,ECB_02009,ECB_02010,ECB_02011,ECB_02012,[ECB_02013]	left_inside_cov=2	left_outside_cov=12	locus_tag=[ECB_01979]–[ECB_02013]	right_inside_cov=0	right_outside_cov=8
MC	107	.	REL606	2157884	2157890	0	0	gene_name=yehM	gene_position=coding (994-1000/2280 nt)	gene_product=hypothetical protein	gene_strand=>	left_inside_cov=0	left_outside_cov=96	locus_tag=ECB_02049	right_inside_cov=4	right_outside_cov=24
MC	108	.	REL606	2240635	2240635	0	0	gene_name=narP	gene_position=coding (147/648 nt)	gene_product=DNA-binding response regulator in two-component regulatory system with NarQ or NarX	gene_strand=>	left_inside_cov=0	left_outside_cov=5	locus_tag=ECB_02120	right_inside_cov=0	right_outside_cov=5
MC	109	.	REL606	2496803	2496808	0	0	gene_name=ypfE/maeB	gene_position=intergenic (-170/+118)	gene_product=predicted carboxysome structural protein with predicted role in ethanol utilization/malic enzyme	gene_strand=</<	left_inside_cov=0	left_outside_cov=8	locus_tag=ECB_02353/ECB_02354	right_inside_cov=3	right_outside_cov=29
MC	110	.	REL606	2628950	2628952	0	0	gene_name=rseC	gene_position=coding (343-345/480 nt)	gene_product=RseC protein involved in reduction of the SoxR iron-sulfur cluster	gene_strand=<	left_inside_cov=4	left_outside_cov=15	locus_tag=ECB_02464	right_inside_cov=0	right_outside_cov=8
MC	111	.	REL606	2815226	2815230	0	0	gene_name=yqcD/ygdH	gene_position=intergenic (+101/-7)	gene_product=hypothetical protein/hypothetical protein	gene_strand=>/>	left_inside_cov=0	left_outside_cov=19	locus_tag=ECB_02639/ECB_02640	right_inside_cov=0	right_outside_cov=19
MC	112	.	REL606	2990103	2990112	0	0	gene_name=mltC	gene_position=coding (10-19/1080 nt)	gene_product=membrane-bound lytic murein transglycosylase C	gene_strand=>	left_inside_cov=0	left_outside_cov=8	locus_tag=ECB_02793	right_inside_cov=0	right_outside_cov=6
MC	113	.	REL606	3001861	3001867	0	0	gene_name=flu	gene_position=coding (2298-2304/2847 nt)	gene_product=antigen 43 (Ag43) phase-variable biofilm formation autotransporter	gene_strand=>	left_inside_cov=4	left_outside_cov=6	locus_tag=ECB_02800	right_inside_cov=3	right_outside_cov=16
MC	114	.	REL606	3015775	3035122	0	0	gene_name=[ECB_02816]–[ECB_02836]	gene_product=[ECB_02816],ECB_02817,ECB_02818,ECB_02819,ECB_02820,ECB_02821,ECB_02822,insB-22,insA-22,ECB_02825,ECB_02826,ECB_02827,ECB_02828,yghD,yghE,ECB_02831,ECB_02832,ECB_02833,ECB_02834,ECB_02835,[ECB_02836]	left_inside_cov=0	left_outside_cov=7	locus_tag=[ECB_02816]–[ECB_02836]	right_inside_cov=0	right_outside_cov=72
MC	115	.	REL606	3042912	3042917	0	0	gene_name=yghJ	gene_position=coding (3087-3092/4554 nt)	gene_product=predicted inner membrane lipoprotein	gene_strand=<	left_inside_cov=4	left_outside_cov=5	locus_tag=ECB_02842	right_inside_cov=4	right_outside_cov=20
MC	116	.	REL606	3093628	3093628	0	0	gene_name=ygiQ	gene_position=coding (1691/2220 nt)	gene_product=hypothetical protein	gene_strand=<	left_inside_cov=0	left_outside_cov=5	locus_tag=ECB_02888	right_inside_cov=0	right_outside_cov=5
MC	117	.	REL606	3260565	3260570	0	0	gene_name=hflB	gene_position=coding (1600-1605/1935 nt)	gene_product=protease, ATP-dependent zinc-metallo	gene_strand=<	left_inside_cov=2	left_outside_cov=59	locus_tag=ECB_03043	right_inside_cov=0	right_outside_cov=60
MC	118	.	REL606	3317226	3317226	0	0	gene_name=yhcR	gene_position=coding (9/204 nt)	gene_product=membrane protein of efflux system	gene_strand=<	left_inside_cov=0	left_outside_cov=6	locus_tag=ECB_03102	right_inside_cov=0	right_outside_cov=15
MC	119	.	REL606	3401857	3401866	0	0	gene_name=rpsG	gene_position=coding (121-130/471 nt)	gene_product=30S ribosomal protein S7	gene_strand=<	left_inside_cov=2	left_outside_cov=10	locus_tag=ECB_03192	right_inside_cov=4	right_outside_cov=9
MC	120	.	REL606	3511244	3511244	0	0	gene_name=yhhX	gene_position=coding (1015/1038 nt)	gene_product=predicted oxidoreductase with NAD(P)-binding Rossmann-fold domain	gene_strand=<	left_inside_cov=0	left_outside_cov=48	locus_tag=ECB_03291	right_inside_cov=0	right_outside_cov=58
MC	121	.	REL606	3549937	3550087	11	124	gene_name=rhsB	gene_position=coding (58-208/4236 nt)	gene_product=rhsB element core protein RshB	gene_strand=>	left_inside_cov=3	left_outside_cov=5	locus_tag=ECB_03331	right_inside_cov=0	right_outside_cov=9
MC	122	.	REL606	3551294	3551296	0	0	gene_name=rhsB	gene_position=coding (1415-1417/4236 nt)	gene_product=rhsB element core protein RshB	gene_strand=>	left_inside_cov=0	left_outside_cov=7	locus_tag=ECB_03331	right_inside_cov=0	right_outside_cov=7
MC	123	.	REL606	3551332	3551351	0	12	gene_name=rhsB	gene_position=coding (1453-1472/4236 nt)	gene_product=rhsB element core protein RshB	gene_strand=>	left_inside_cov=0	left_outside_cov=7	locus_tag=ECB_03331	right_inside_cov=4	right_outside_cov=17
MC	124	.	REL606	3551390	3551405	2	0	gene_name=rhsB	gene_position=coding (1511-1526/4236 nt)	gene_product=rhsB element core protein RshB	gene_strand=>	left_inside_cov=2	left_outside_cov=5	locus_tag=ECB_03331	right_inside_cov=0	right_outside_cov=17
MC	125	.	REL606	3552258	3552313	0	25	gene_name=rhsB	gene_position=coding (2379-2434/4236 nt)	gene_product=rhsB element core protein RshB	gene_strand=>	left_inside_cov=1	left_outside_cov=15	locus_tag=ECB_03331	right_inside_cov=4	right_outside_cov=10
MC	126	.	REL606	3552373	3552392	0	0	gene_name=rhsB	gene_position=coding (2494-2513/4236 nt)	gene_product=rhsB element core protein RshB	gene_strand=>	left_inside_cov=0	left_outside_cov=14	locus_tag=ECB_03331	right_inside_cov=0	right_outside_cov=7
MC	127	.	REL606	3552431	3552447	0	9	gene_name=rhsB	gene_position=coding (2552-2568/4236 nt)	gene_product=rhsB element core protein RshB	gene_strand=>	left_inside_cov=2	left_outside_cov=7	locus_tag=ECB_03331	right_inside_cov=4	right_outside_cov=5
MC	128	.	REL606	3552535	3552542	0	3	gene_name=rhsB	gene_position=coding (2656-2663/4236 nt)	gene_product=rhsB element core protein RshB	gene_strand=>	left_inside_cov=1	left_outside_cov=12	locus_tag=ECB_03331	right_inside_cov=4	right_outside_cov=17
MC	129	.	REL606	3552685	3552713	9	6	gene_name=rhsB	gene_position=coding (2806-2834/4236 nt)	gene_product=rhsB element core protein RshB	gene_strand=>	left_inside_cov=3	left_outside_cov=11	locus_tag=ECB_03331	right_inside_cov=0	right_outside_cov=15
MC	130	.	REL606	3553222	3553254	28	0	gene_name=rhsB	gene_position=coding (3343-3375/4236 nt)	gene_product=rhsB element core protein RshB	gene_strand=>	left_inside_cov=4	left_outside_cov=10	locus_tag=ECB_03331	right_inside_cov=0	right_outside_cov=7
MC	131	.	REL606	3577451	3577460	0	0	gene_name=gor	gene_position=coding (463-472/1353 nt)	gene_product=glutathione reductase	gene_strand=>	left_inside_cov=0	left_outside_cov=15	locus_tag=ECB_03349	right_inside_cov=2	right_outside_cov=5
MC	132	.	REL606	3741958	3742148	5	0	gene_name=[waaT]	gene_product=[waaT]	left_inside_cov=2	left_outside_cov=7	locus_tag=[ECB_03483]	right_inside_cov=0	right_outside_cov=5
MC	133	.	REL606	3868895	3868904	0	0	gene_name=pstB/pstA	gene_position=intergenic (-162/+12)	gene_product=phosphate transporter subunit/phosphate transporter subunit	gene_strand=</<	left_inside_cov=0	left_outside_cov=5	locus_tag=ECB_03609/ECB_03610	right_inside_cov=0	right_outside_cov=18
MC	134	.	REL606	3894999	3901455	0	0	gene_name=rbsD–[yieO]	gene_product=rbsD,rbsA,rbsC,rbsB,rbsK,rbsR,[yieO]	left_inside_cov=3	left_outside_cov=8	locus_tag=[ECB_03634]–[ECB_03640]	right_inside_cov=3	right_outside_cov=20
MC	135	.	REL606	3951678	3951693	0	0	gene_name=hemC	gene_position=coding (904-919/942 nt)	gene_product=porphobilinogen deaminase	gene_strand=<	left_inside_cov=0	left_outside_cov=8	locus_tag=ECB_03680	right_inside_cov=3	right_outside_cov=22
MC	136	.	REL606	4017747	4017769	8	5	gene_name=rrlA	gene_position=noncoding (1856-1878/2904 nt)	gene_product=23S ribosomal RNA	gene_strand=>	left_inside_cov=0	left_outside_cov=5	locus_tag=ECB_r00015	right_inside_cov=0	right_outside_cov=20
MC	137	.	REL606	4046102	4046104	0	0	gene_name=yihQ	gene_position=coding (1394-1396/2004 nt)	gene_product=alpha-glucosidase	gene_strand=<	left_inside_cov=0	left_outside_cov=6	locus_tag=ECB_03763	right_inside_cov=0	right_outside_cov=6
MC	138	.	REL606	4108227	4108231	0	0	gene_name=metJ	gene_position=coding (272-276/318 nt)	gene_product=transcriptional repressor protein MetJ	gene_strand=<	left_inside_cov=0	left_outside_cov=7	locus_tag=ECB_03824	right_inside_cov=0	right_outside_cov=5
MC	139	.	REL606	4110233	4110237	0	0	gene_name=metL	gene_position=coding (292-296/2433 nt)	gene_product=bifunctional aspartate kinase II/homoserine dehydrogenase II	gene_strand=>	left_inside_cov=0	left_outside_cov=13	locus_tag=ECB_03826	right_inside_cov=3	right_outside_cov=14
MC	140	.	REL606	4357682	4357688	0	0	gene_name=ampC/frdD	gene_position=intergenic (-27/+30)	gene_product=beta-lactamase/D-alanine carboxypeptidase/fumarate reductase subunit D	gene_strand=</<	left_inside_cov=3	left_outside_cov=6	locus_tag=ECB_04022/ECB_04023	right_inside_cov=0	right_outside_cov=6
MC	141	.	REL606	4507891	4507893	0	0	gene_name=sgcA	gene_position=coding (226-228/432 nt)	gene_product=predicted phosphotransferase enzyme IIA component	gene_strand=<	left_inside_cov=0	left_outside_cov=25	locus_tag=ECB_04167	right_inside_cov=2	right_outside_cov=17
MC	142	.	REL606	4522279	4561282	59	0	gene_name=[insA-27]–[hsdR]	gene_product=[insA-27],fimB,fimE,fimA,fimI,fimC,fimD,fimF,fimG,fimH,gntP,uxuA,uxuB,uxuR,yjiC,yjiD,yjiE,iadA,yjiG,yjiH,kptA,yjiJ,yjiK,yjiL,yjiM,yjiN,yjiO,yjiPQ,insA-28,insB-28,yjiV,mcrC,mcrB,yjiW,hsdS,hsdM,[hsdR]	left_inside_cov=3	left_outside_cov=6	locus_tag=[ECB_04181]–[ECB_04216]	right_inside_cov=0	right_outside_cov=45
JC	143	.	REL606	1	1	REL606	4629812	-1	0	alignment_overlap=0	coverage_minus=59	coverage_plus=73	flanking_left=36	flanking_right=36	frequency=NA	frequency_lower=NA	frequency_upper=NA	ignore=CIRCULAR_CHROMOSOME	junction_effective_depth=131.94	junction_possible_overlap_registers=32	junction_possible_overlap_registers_before_trimming=34	key=REL606__1__1__REL606__4629812__-1__0____36__36__0__0	max_left=33	max_left_minus=29	max_left_plus=33	max_min_left=17	max_min_left_minus=0	max_min_left_plus=17	max_min_right=15	max_min_right_minus=15	max_min_right_plus=15	max_pos_hash_score=68	max_right=19	max_right_minus=15	max_right_plus=19	neg_log10_pos_hash_p_value=0.1	new_junction_coverage=1.98	new_junction_read_count=132	new_junction_reference_weighted_read_count=0.89	new_junction_weighted_read_count=131.11	pos_hash_score=16	prediction=unknown	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=NA	side_1_gene_name=–/thrL	side_1_gene_position=intergenic (–/-189)	side_1_gene_product=–/thr operon leader peptide	side_1_gene_strand=–/>	side_1_locus_tag=–/ECB_00001	side_1_overlap=0	side_1_possible_overlap_registers=0	side_1_possible_overlap_registers_before_trimming=34	side_1_read_count=0	side_1_redundant=0	side_1_weighted_read_count=0.00	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=NA	side_2_gene_name=lasT/–	side_2_gene_position=intergenic (+24/–)	side_2_gene_product=predicted rRNA methyltransferase/–	side_2_gene_strand=>/–	side_2_locus_tag=ECB_04279/–	side_2_overlap=0	side_2_possible_overlap_registers=0	side_2_possible_overlap_registers_before_trimming=34	side_2_read_count=0	side_2_redundant=0	side_2_weighted_read_count=0.00	total_non_overlap_reads=132
JC	144	.	REL606	15386	1	REL606	3517313	-1	0	alignment_overlap=3	coverage_minus=27	coverage_plus=37	flanking_left=36	flanking_right=36	frequency=8.657e-01	frequency_lower=7.773e-01	frequency_upper=9.281e-01	junction_effective_depth=67.00	junction_mixture_iterations=1	junction_possible_overlap_registers=29	junction_possible_overlap_registers_before_trimming=31	key=REL606__590715__1__REL606__3517316__-1__3____36__36__1__0	max_left=28	max_left_minus=28	max_left_plus=25	max_min_left=15	max_min_left_minus=15	max_min_left_plus=8	max_min_right=14	max_min_right_minus=5	max_min_right_plus=14	max_pos_hash_score=62	max_right=31	max_right_minus=22	max_right_plus=31	neg_log10_pos_hash_p_value=0.2	new_junction_coverage=0.96	new_junction_read_count=58	new_junction_reference_weighted_read_count=0.00	new_junction_weighted_read_count=58.00	pos_hash_score=11	prediction=consensus	side_1_annotate_key=repeat	side_1_continuation=0	side_1_coverage=NA	side_1_gene_name=IS186	side_1_gene_position=noncoding (1/1343 nt)	side_1_gene_product=repeat region	side_1_gene_strand=>	side_1_overlap=3	side_1_possible_overlap_registers=NA	side_1_possible_overlap_registers_before_trimming=NA	side_1_read_count=NA	side_1_redundant=1	side_1_weighted_read_count=NA	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=0.15	side_2_gene_name=ggt	side_2_gene_position=coding (187/1743 nt)	side_2_gene_product=gamma-glutamyltranspeptidase periplasmic precursor	side_2_gene_strand=<	side_2_locus_tag=ECB_03296	side_2_overlap=0	side_2_possible_overlap_registers=29	side_2_possible_overlap_registers_before_trimming=31	side_2_read_count=9	side_2_redundant=0	side_2_weighted_read_count=9.00	total_non_overlap_reads=64
JC	145	.	REL606	16992	-1	REL606	664688	1	0	alignment_overlap=0	coverage_minus=15	coverage_plus=14	flanking_left=36	flanking_right=36	frequency=5.385e-01	frequency_lower=3.958e-01	frequency_upper=6.765e-01	junction_effective_depth=39.00	junction_mixture_iterations=1	junction_possible_overlap_registers=32	junction_possible_overlap_registers_before_trimming=34	key=REL606__16992__-1__REL606__2775877__-1__0____36__36__0__1	max_left=34	max_left_minus=34	max_left_plus=26	max_min_left=14	max_min_left_minus=14	max_min_left_plus=12	max_min_right=17	max_min_right_minus=17	max_min_right_plus=10	max_pos_hash_score=68	max_right=30	max_right_minus=23	max_right_plus=30	neg_log10_pos_hash_p_value=0.4	new_junction_coverage=0.31	new_junction_read_count=21	new_junction_reference_weighted_read_count=0.00	new_junction_weighted_read_count=21.00	pos_hash_score=9	prediction=consensus	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=0.27	side_1_gene_name=mokC/nhaA	side_1_gene_position=intergenic (-34/-496)	side_1_gene_product=regulatory protein for HokC, overlaps CDS of hokC/pH-dependent sodium/proton antiporter	side_1_gene_strand=</>	side_1_locus_tag=ECB_00017/ECB_00018	side_1_overlap=0	side_1_possible_overlap_registers=32	side_1_possible_overlap_registers_before_trimming=34	side_1_read_count=18	side_1_redundant=0	side_1_weighted_read_count=18.00	side_2_annotate_key=repeat	side_2_continuation=0	side_2_coverage=NA	side_2_gene_name=IS150	side_2_gene_position=noncoding (1/1443 nt)	side_2_gene_product=repeat region	side_2_gene_strand=>	side_2_overlap=0	side_2_possible_overlap_registers=NA	side_2_possible_overlap_registers_before_trimming=NA	side_2_read_count=NA	side_2_redundant=1	side_2_weighted_read_count=NA	total_non_overlap_reads=29
JC	146	.	REL606	546191	1	REL606	1194382	-1	0	alignment_overlap=1	coverage_minus=32	coverage_plus=18	flanking_left=36	flanking_right=36	frequency=8.495e-01	frequency_lower=7.502e-01	frequency_upper=9.201e-01	junction_effective_depth=58.00	junction_mixture_iterations=1	junction_possible_overlap_registers=31	junction_possible_overlap_registers_before_trimming=33	key=REL606__546190__1__REL606__546932__1__1____36__36__0__1	max_left=30	max_left_minus=29	max_left_plus=30	max_min_left=16	max_min_left_minus=16	max_min_left_plus=13	max_min_right=13	max_min_right_minus=9	max_min_right_plus=13	max_pos_hash_score=66	max_right=22	max_right_minus=21	max_right_plus=22	neg_log10_pos_hash_p_value=0.2	new_junction_coverage=0.77	new_junction_read_count=50	new_junction_reference_weighted_read_count=0.00	new_junction_weighted_read_count=50.00	pos_hash_score=12	prediction=consensus	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=0.14	side_1_gene_name=ybcQ	side_1_gene_position=coding (346/384 nt)	side_1_gene_product=predicted antitermination protein	side_1_gene_strand=>	side_1_locus_tag=ECB_00502	side_1_overlap=0	side_1_possible_overlap_registers=28	side_1_possible_overlap_registers_before_trimming=33	side_1_read_count=8	side_1_redundant=0	side_1_weighted_read_count=8.00	side_2_annotate_key=repeat	side_2_continuation=0	side_2_coverage=NA	side_2_gene_name=IS1	side_2_gene_position=noncoding (768/768 nt)	side_2_gene_product=repeat region	side_2_gene_strand=>	side_2_overlap=1	side_2_possible_overlap_registers=NA	side_2_possible_overlap_registers_before_trimming=NA	side_2_read_count=NA	side_2_redundant=1	side_2_weighted_read_count=NA	total_non_overlap_reads=50
JC	147	.	REL606	643364	-1	REL606	664688	1	0	alignment_overlap=1	coverage_minus=20	coverage_plus=51	flanking_left=36	flanking_right=36	frequency=1.000e+00	frequency_lower=9.050e-01	frequency_upper=1.000e+00	junction_effective_depth=30.00	junction_mixture_iterations=1	junction_possible_overlap_registers=25	junction_possible_overlap_registers_before_trimming=33	key=REL606__643365__-1__REL606__664688__1__1____36__36__0__1	max_left=31	max_left_minus=31	max_left_plus=31	max_min_left=16	max_min_left_minus=16	max_min_left_plus=12	max_min_right=16	max_min_right_minus=16	max_min_right_plus=8	max_pos_hash_score=66	max_right=31	max_right_minus=23	max_right_plus=31	neg_log10_pos_hash_p_value=0.2	new_junction_coverage=0.57	new_junction_read_count=30	new_junction_reference_weighted_read_count=0.00	new_junction_weighted_read_count=30.00	pos_hash_score=13	prediction=consensus	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=0.00	side_1_gene_name=ybeF	side_1_gene_position=coding (598/954 nt)	side_1_gene_product=predicted DNA-binding transcriptional regulator	side_1_gene_strand=<	side_1_locus_tag=ECB_00598	side_1_overlap=0	side_1_possible_overlap_registers=24	side_1_possible_overlap_registers_before_trimming=33	side_1_read_count=0	side_1_redundant=0	side_1_weighted_read_count=0.00	side_2_annotate_key=repeat	side_2_continuation=0	side_2_coverage=NA	side_2_gene_name=IS150	side_2_gene_position=noncoding (1/1443 nt)	side_2_gene_product=repeat region	side_2_gene_strand=>	side_2_overlap=1	side_2_possible_overlap_registers=NA	side_2_possible_overlap_registers_before_trimming=NA	side_2_read_count=NA	side_2_redundant=1	side_2_weighted_read_count=NA	total_non_overlap_reads=71
JC	148	.	REL606	664686	1	REL606	4613019	-1	0	alignment_overlap=0	coverage_minus=16	coverage_plus=95	flanking_left=36	flanking_right=36	frequency=1.000e+00	frequency_lower=9.680e-01	frequency_upper=1.000e+00	junction_effective_depth=92.00	junction_mixture_iterations=1	junction_possible_overlap_registers=28	junction_possible_overlap_registers_before_trimming=34	key=REL606__664686__1__REL606__4613019__-1__0____36__36__1__0	max_left=35	max_left_minus=31	max_left_plus=35	max_min_left=14	max_min_left_minus=14	max_min_left_plus=10	max_min_right=15	max_min_right_minus=8	max_min_right_plus=15	max_pos_hash_score=68	max_right=33	max_right_minus=24	max_right_plus=33	neg_log10_pos_hash_p_value=0.1	new_junction_coverage=1.57	new_junction_read_count=92	new_junction_reference_weighted_read_count=0.00	new_junction_weighted_read_count=92.00	pos_hash_score=17	prediction=consensus	read_count_offset=0	side_1_annotate_key=repeat	side_1_continuation=0	side_1_coverage=NA	side_1_gene_name=gltL/insJ-2	side_1_gene_position=intergenic (-77/-49)	side_1_gene_product=ATP-binding protein of glutamate/aspartate transport system; b0652_1/IS150 hypothetical protein	side_1_gene_strand=</>	side_1_locus_tag=ECB_00618/ECB_00619	side_1_overlap=0	side_1_possible_overlap_registers=NA	side_1_possible_overlap_registers_before_trimming=NA	side_1_read_count=NA	side_1_redundant=1	side_1_weighted_read_count=NA	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=0.00	side_2_gene_name=smp/serB	side_2_gene_position=intergenic (-16/-90)	side_2_gene_product=hypothetical protein/3-phosphoserine phosphatase	side_2_gene_strand=</>	side_2_locus_tag=ECB_04263/ECB_04264	side_2_overlap=0	side_2_possible_overlap_registers=30	side_2_possible_overlap_registers_before_trimming=34	side_2_read_count=0	side_2_redundant=0	side_2_weighted_read_count=0.00	total_non_overlap_reads=111
JC	149	.	REL606	664688	1	REL606	1451972	1	0	alignment_overlap=2	coverage_minus=35	coverage_plus=21	flanking_left=36	flanking_right=36	frequency=1.000e+00	frequency_lower=9.202e-01	frequency_upper=1.000e+00	junction_effective_depth=36.00	junction_mixture_iterations=1	junction_possible_overlap_registers=28	junction_possible_overlap_registers_before_trimming=32	key=REL606__664688__1__REL606__1451970__1__2____36__36__1__0	max_left=30	max_left_minus=29	max_left_plus=30	max_min_left=13	max_min_left_minus=10	max_min_left_plus=13	max_min_right=12	max_min_right_minus=12	max_min_right_plus=4	max_pos_hash_score=64	max_right=32	max_right_minus=31	max_right_plus=32	neg_log10_pos_hash_p_value=0.3	new_junction_coverage=0.62	new_junction_read_count=36	new_junction_reference_weighted_read_count=0.00	new_junction_weighted_read_count=36.00	pos_hash_score=11	prediction=consensus	read_count_offset=3	side_1_annotate_key=repeat	side_1_continuation=0	side_1_coverage=NA	side_1_gene_name=IS150	side_1_gene_position=noncoding (1/1443 nt)	side_1_gene_product=repeat region	side_1_gene_strand=>	side_1_overlap=2	side_1_possible_overlap_registers=NA	side_1_possible_overlap_registers_before_trimming=NA	side_1_read_count=NA	side_1_redundant=1	side_1_weighted_read_count=NA	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=0.00	side_2_gene_name=acpD	side_2_gene_position=coding (220/606 nt)	side_2_gene_product=acyl carrier protein phosphodiesterase	side_2_gene_strand=<	side_2_locus_tag=ECB_01367	side_2_overlap=0	side_2_possible_overlap_registers=19	side_2_possible_overlap_registers_before_trimming=29	side_2_read_count=0	side_2_redundant=0	side_2_weighted_read_count=0.00	total_non_overlap_reads=56
JC	150	.	REL606	664688	1	REL606	3035123	1	0	alignment_overlap=2	coverage_minus=41	coverage_plus=24	flanking_left=36	flanking_right=36	frequency=1.000e+00	frequency_lower=9.418e-01	frequency_upper=1.000e+00	junction_effective_depth=50.00	junction_mixture_iterations=1	junction_possible_overlap_registers=29	junction_possible_overlap_registers_before_trimming=32	key=REL606__2775877__-1__REL606__3035121__1__2____36__36__1__0	max_left=31	max_left_minus=31	max_left_plus=30	max_min_left=13	max_min_left_minus=10	max_min_left_plus=13	max_min_right=14	max_min_right_minus=13	max_min_right_plus=14	max_pos_hash_score=64	max_right=31	max_right_minus=31	max_right_plus=22	neg_log10_pos_hash_p_value=0.1	new_junction_coverage=0.83	new_junction_read_count=50	new_junction_reference_weighted_read_count=0.00	new_junction_weighted_read_count=50.00	pos_hash_score=14	prediction=consensus	side_1_annotate_key=repeat	side_1_continuation=0	side_1_coverage=NA	side_1_gene_name=IS150	side_1_gene_position=noncoding (1/1443 nt)	side_1_gene_product=repeat region	side_1_gene_strand=>	side_1_overlap=2	side_1_possible_overlap_registers=NA	side_1_possible_overlap_registers_before_trimming=NA	side_1_read_count=NA	side_1_redundant=1	side_1_weighted_read_count=NA	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=0.00	side_2_gene_name=ECB_02836	side_2_gene_position=coding (319/1224 nt)	side_2_gene_product=hypothetical type II secretion protein GspF	side_2_gene_strand=<	side_2_locus_tag=ECB_02836	side_2_overlap=0	side_2_possible_overlap_registers=26	side_2_possible_overlap_registers_before_trimming=32	side_2_read_count=0	side_2_redundant=0	side_2_weighted_read_count=0.00	total_non_overlap_reads=65
JC	151	.	REL606	666130	-1	REL606	1047818	1	0	alignment_overlap=1	coverage_minus=65	coverage_plus=23	flanking_left=36	flanking_right=36	frequency=9.670e-01	frequency_lower=9.170e-01	frequency_upper=9.910e-01	junction_effective_depth=91.00	junction_mixture_iterations=1	junction_possible_overlap_registers=30	junction_possible_overlap_registers_before_trimming=33	key=REL606__666130__-1__REL606__1047817__1__1____36__36__1__0	max_left=30	max_left_minus=30	max_left_plus=28	max_min_left=17	max_min_left_minus=17	max_min_left_plus=4	max_min_right=16	max_min_right_minus=12	max_min_right_plus=16	max_pos_hash_score=66	max_right=29	max_right_minus=20	max_right_plus=29	neg_log10_pos_hash_p_value=0.1	new_junction_coverage=1.40	new_junction_read_count=88	new_junction_reference_weighted_read_count=0.00	new_junction_weighted_read_count=88.00	pos_hash_score=15	prediction=consensus	side_1_annotate_key=repeat	side_1_continuation=0	side_1_coverage=NA	side_1_gene_name=IS150	side_1_gene_position=noncoding (1443/1443 nt)	side_1_gene_product=repeat region	side_1_gene_strand=>	side_1_overlap=1	side_1_possible_overlap_registers=NA	side_1_possible_overlap_registers_before_trimming=NA	side_1_read_count=NA	side_1_redundant=1	side_1_weighted_read_count=NA	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=0.05	side_2_gene_name=yccK/yccA	side_2_gene_position=intergenic (-59/+32)	side_2_gene_product=predicted sulfite reductase subunit/inner membrane protein	side_2_gene_strand=</<	side_2_locus_tag=ECB_00973/ECB_00974	side_2_overlap=0	side_2_possible_overlap_registers=30	side_2_possible_overlap_registers_before_trimming=33	side_2_read_count=3	side_2_redundant=0	side_2_weighted_read_count=3.00	total_non_overlap_reads=88
JC	152	.	REL606	666130	-1	REL606	1451974	-1	0	alignment_overlap=0	coverage_minus=43	coverage_plus=16	flanking_left=36	flanking_right=36	frequency=1.000e+00	frequency_lower=9.470e-01	frequency_upper=1.000e+00	junction_effective_depth=55.00	junction_mixture_iterations=1	junction_possible_overlap_registers=27	junction_possible_overlap_registers_before_trimming=34	key=REL606__666130__-1__REL606__1451974__-1__0____36__36__1__0	max_left=31	max_left_minus=31	max_left_plus=19	max_min_left=16	max_min_left_minus=16	max_min_left_plus=13	max_min_right=18	max_min_right_minus=18	max_min_right_plus=15	max_pos_hash_score=68	max_right=31	max_right_minus=21	max_right_plus=31	neg_log10_pos_hash_p_value=0.1	new_junction_coverage=0.98	new_junction_read_count=55	new_junction_reference_weighted_read_count=0.00	new_junction_weighted_read_count=55.00	pos_hash_score=16	prediction=consensus	read_count_offset=3	side_1_annotate_key=repeat	side_1_continuation=0	side_1_coverage=NA	side_1_gene_name=IS150	side_1_gene_position=noncoding (1443/1443 nt)	side_1_gene_product=repeat region	side_1_gene_strand=>	side_1_overlap=0	side_1_possible_overlap_registers=NA	side_1_possible_overlap_registers_before_trimming=NA	side_1_read_count=NA	side_1_redundant=1	side_1_weighted_read_count=NA	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=0.00	side_2_gene_name=acpD	side_2_gene_position=coding (218/606 nt)	side_2_gene_product=acyl carrier protein phosphodiesterase	side_2_gene_strand=<	side_2_locus_tag=ECB_01367	side_2_overlap=0	side_2_possible_overlap_registers=24	side_2_possible_overlap_registers_before_trimming=31	side_2_read_count=0	side_2_redundant=0	side_2_weighted_read_count=0.00	total_non_overlap_reads=59
JC	153	.	REL606	666130	-1	REL606	1462253	-1	0	alignment_overlap=0	coverage_minus=44	coverage_plus=81	flanking_left=36	flanking_right=36	frequency=1.000e+00	frequency_lower=9.763e-01	frequency_upper=1.000e+00	junction_effective_depth=125.00	junction_mixture_iterations=1	junction_possible_overlap_registers=31	junction_possible_overlap_registers_before_trimming=34	key=REL606__666130__-1__REL606__1462253__-1__0____36__36__1__0	max_left=30	max_left_minus=30	max_left_plus=25	max_min_left=17	max_min_left_minus=17	max_min_left_plus=16	max_min_right=16	max_min_right_minus=11	max_min_right_plus=16	max_pos_hash_score=68	max_right=28	max_right_minus=23	max_right_plus=28	neg_log10_pos_hash_p_value=0.1	new_junction_coverage=1.93	new_junction_read_count=125	new_junction_reference_weighted_read_count=0.00	new_junction_weighted_read_count=125.00	pos_hash_score=17	prediction=consensus	side_1_annotate_key=repeat	side_1_continuation=0	side_1_coverage=NA	side_1_gene_name=IS150	side_1_gene_position=noncoding (1443/1443 nt)	side_1_gene_product=repeat region	side_1_gene_strand=>	side_1_overlap=0	side_1_possible_overlap_registers=NA	side_1_possible_overlap_registers_before_trimming=NA	side_1_read_count=NA	side_1_redundant=1	side_1_weighted_read_count=NA	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=0.00	side_2_gene_name=mokB/trg	side_2_gene_position=intergenic (-15/-326)	side_2_gene_product=regulatory peptide/methyl-accepting chemotaxis protein III, ribose and galactose sensor receptor	side_2_gene_strand=</>	side_2_locus_tag=ECB_01377/ECB_01378	side_2_overlap=0	side_2_possible_overlap_registers=31	side_2_possible_overlap_registers_before_trimming=34	side_2_read_count=0	side_2_redundant=0	side_2_weighted_read_count=0.00	total_non_overlap_reads=125
JC	154	.	REL606	666130	-1	REL606	1651966	1	0	alignment_overlap=0	coverage_minus=78	coverage_plus=43	flanking_left=36	flanking_right=36	frequency=9.274e-01	frequency_lower=8.768e-01	frequency_upper=9.616e-01	junction_effective_depth=124.00	junction_mixture_iterations=1	junction_possible_overlap_registers=30	junction_possible_overlap_registers_before_trimming=34	key=REL606__666130__-1__REL606__1651966__1__0____36__36__1__0	max_left=31	max_left_minus=31	max_left_plus=31	max_min_left=17	max_min_left_minus=17	max_min_left_plus=8	max_min_right=18	max_min_right_minus=15	max_min_right_plus=18	max_pos_hash_score=68	max_right=31	max_right_minus=31	max_right_plus=26	neg_log10_pos_hash_p_value=0.1	new_junction_coverage=1.84	new_junction_read_count=115	new_junction_reference_weighted_read_count=0.00	new_junction_weighted_read_count=115.00	pos_hash_score=18	prediction=consensus	side_1_annotate_key=repeat	side_1_continuation=0	side_1_coverage=NA	side_1_gene_name=IS150	side_1_gene_position=noncoding (1443/1443 nt)	side_1_gene_product=repeat region	side_1_gene_strand=>	side_1_overlap=0	side_1_possible_overlap_registers=NA	side_1_possible_overlap_registers_before_trimming=NA	side_1_read_count=NA	side_1_redundant=1	side_1_weighted_read_count=NA	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=0.14	side_2_gene_name=ydgG	side_2_gene_position=coding (859/1035 nt)	side_2_gene_product=predicted inner membrane protein	side_2_gene_strand=>	side_2_locus_tag=ECB_01570	side_2_overlap=0	side_2_possible_overlap_registers=30	side_2_possible_overlap_registers_before_trimming=34	side_2_read_count=9	side_2_redundant=0	side_2_weighted_read_count=9.00	total_non_overlap_reads=121
JC	155	.	REL606	666130	-1	REL606	2086596	-1	0	alignment_overlap=0	coverage_minus=54	coverage_plus=76	flanking_left=36	flanking_right=36	frequency=NA	frequency_lower=NA	frequency_upper=NA	junction_effective_depth=130.00	junction_possible_overlap_registers=32	junction_possible_overlap_registers_before_trimming=34	key=REL606__666130__-1__REL606__2086596__-1__0____36__36__1__1	max_left=31	max_left_minus=31	max_left_plus=31	max_min_left=17	max_min_left_minus=17	max_min_left_plus=4	max_min_right=12	max_min_right_minus=12	max_min_right_plus=10	max_pos_hash_score=68	max_right=33	max_right_minus=21	max_right_plus=33	neg_log10_pos_hash_p_value=0.1	new_junction_coverage=1.95	new_junction_read_count=130	new_junction_reference_weighted_read_count=0.00	new_junction_weighted_read_count=130.00	pos_hash_score=15	prediction=unknown	side_1_annotate_key=repeat	side_1_continuation=0	side_1_coverage=NA	side_1_gene_name=IS150	side_1_gene_position=noncoding (1443/1443 nt)	side_1_gene_product=repeat region	side_1_gene_strand=>	side_1_overlap=0	side_1_possible_overlap_registers=NA	side_1_possible_overlap_registers_before_trimming=NA	side_1_read_count=NA	side_1_redundant=1	side_1_weighted_read_count=NA	side_2_annotate_key=repeat	side_2_continuation=0	side_2_coverage=NA	side_2_gene_name=yegL/yegM	side_2_gene_position=intergenic (-1017/-528)	side_2_gene_product=hypothetical protein/multidrug efflux system, subunit A	side_2_gene_strand=</>	side_2_locus_tag=ECB_01978/ECB_01979	side_2_overlap=0	side_2_possible_overlap_registers=NA	side_2_possible_overlap_registers_before_trimming=NA	side_2_read_count=NA	side_2_redundant=1	side_2_weighted_read_count=NA	total_non_overlap_reads=130
JC	156	.	REL606	666130	-1	REL606	4415712	-1	0	alignment_overlap=1	coverage_minus=30	coverage_plus=17	flanking_left=36	flanking_right=36	frequency=1.000e+00	frequency_lower=9.327e-01	frequency_upper=1.000e+00	junction_effective_depth=43.00	junction_mixture_iterations=1	junction_possible_overlap_registers=29	junction_possible_overlap_registers_before_trimming=33	key=REL606__3894996__-1__REL606__4415713__-1__1____36__36__1__0	max_left=30	max_left_minus=29	max_left_plus=30	max_min_left=16	max_min_left_minus=16	max_min_left_plus=14	max_min_right=8	max_min_right_minus=8	max_min_right_plus=6	max_pos_hash_score=66	max_right=31	max_right_minus=31	max_right_plus=29	neg_log10_pos_hash_p_value=0.1	new_junction_coverage=0.71	new_junction_read_count=43	new_junction_reference_weighted_read_count=0.00	new_junction_weighted_read_count=43.00	pos_hash_score=14	prediction=consensus	side_1_annotate_key=repeat	side_1_continuation=0	side_1_coverage=NA	side_1_gene_name=IS150	side_1_gene_position=noncoding (1443/1443 nt)	side_1_gene_product=repeat region	side_1_gene_strand=>	side_1_overlap=1	side_1_possible_overlap_registers=NA	side_1_possible_overlap_registers_before_trimming=NA	side_1_read_count=NA	side_1_redundant=1	side_1_weighted_read_count=NA	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=0.00	side_2_gene_name=cycA	side_2_gene_position=coding (850/1413 nt)	side_2_gene_product=D-alanine/D-serine/glycine transporter	side_2_gene_strand=>	side_2_locus_tag=ECB_04080	side_2_overlap=0	side_2_possible_overlap_registers=28	side_2_possible_overlap_registers_before_trimming=33	side_2_read_count=0	side_2_redundant=0	side_2_weighted_read_count=0.00	total_non_overlap_reads=47
JC	157	.	REL606	666130	-1	REL606	4613020	1	0	alignment_overlap=1	coverage_minus=71	coverage_plus=25	flanking_left=36	flanking_right=36	frequency=1.000e+00	frequency_lower=9.693e-01	frequency_upper=1.000e+00	junction_effective_depth=96.00	junction_mixture_iterations=1	junction_possible_overlap_registers=31	junction_possible_overlap_registers_before_trimming=33	key=REL606__3894996__-1__REL606__4613019__1__1____36__36__1__0	max_left=31	max_left_minus=31	max_left_plus=21	max_min_left=16	max_min_left_minus=16	max_min_left_plus=13	max_min_right=13	max_min_right_minus=13	max_min_right_plus=12	max_pos_hash_score=66	max_right=29	max_right_minus=22	max_right_plus=29	neg_log10_pos_hash_p_value=0.1	new_junction_coverage=1.48	new_junction_read_count=96	new_junction_reference_weighted_read_count=0.00	new_junction_weighted_read_count=96.00	pos_hash_score=18	prediction=consensus	read_count_offset=0	side_1_annotate_key=repeat	side_1_continuation=0	side_1_coverage=NA	side_1_gene_name=IS150	side_1_gene_position=noncoding (1443/1443 nt)	side_1_gene_product=repeat region	side_1_gene_strand=>	side_1_overlap=1	side_1_possible_overlap_registers=NA	side_1_possible_overlap_registers_before_trimming=NA	side_1_read_count=NA	side_1_redundant=1	side_1_weighted_read_count=NA	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=0.00	side_2_gene_name=smp/serB	side_2_gene_position=intergenic (-17/-89)	side_2_gene_product=hypothetical protein/3-phosphoserine phosphatase	side_2_gene_strand=</>	side_2_locus_tag=ECB_04263/ECB_04264	side_2_overlap=0	side_2_possible_overlap_registers=30	side_2_possible_overlap_registers_before_trimming=33	side_2_read_count=0	side_2_redundant=0	side_2_weighted_read_count=0.00	total_non_overlap_reads=96
JC	158	.	REL606	816050	1	REL606	818967	-1	0	alignment_overlap=0	coverage_minus=7	coverage_plus=0	flanking_left=36	flanking_right=36	frequency=1.153e-01	frequency_lower=5.594e-02	frequency_upper=2.046e-01	junction_effective_depth=62.00	junction_mixture_iterations=1	junction_possible_overlap_registers=29	junction_possible_overlap_registers_before_trimming=34	key=REL606__816050__1__REL606__818967__-1__0____36__36__0__0	max_left=14	max_left_minus=14	max_left_plus=0	max_min_left=14	max_min_left_minus=14	max_min_left_plus=0	max_min_right=0	max_min_right_minus=0	max_min_right_plus=0	max_pos_hash_score=68	max_right=22	max_right_minus=22	max_right_plus=0	neg_log10_pos_hash_p_value=1.4	new_junction_coverage=0.12	new_junction_read_count=7	new_junction_reference_weighted_read_count=0.00	new_junction_weighted_read_count=7.00	pos_hash_score=3	prediction=polymorphism	reject=FREQUENCY_CUTOFF	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=1.35	side_1_gene_name=ybhO	side_1_gene_position=coding (315/1242 nt)	side_1_gene_product=cardiolipin synthase 2	side_1_gene_strand=<	side_1_locus_tag=ECB_00756	side_1_overlap=0	side_1_possible_overlap_registers=29	side_1_possible_overlap_registers_before_trimming=34	side_1_read_count=82	side_1_redundant=0	side_1_weighted_read_count=82.00	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=0.42	side_2_gene_name=ybhS	side_2_gene_position=coding (911/1134 nt)	side_2_gene_product=predicted transporter subunit: membrane component of ABC superfamily	side_2_gene_strand=<	side_2_locus_tag=ECB_00760	side_2_overlap=0	side_2_possible_overlap_registers=32	side_2_possible_overlap_registers_before_trimming=34	side_2_read_count=28	side_2_redundant=0	side_2_weighted_read_count=28.00	total_non_overlap_reads=7
JC	159	.	REL606	1009100	-1	REL606	3511243	-1	0	alignment_overlap=4	coverage_minus=31	coverage_plus=17	flanking_left=36	flanking_right=36	frequency=6.545e-01	frequency_lower=5.430e-01	frequency_upper=7.545e-01	junction_effective_depth=62.00	junction_mixture_iterations=1	junction_possible_overlap_registers=25	junction_possible_overlap_registers_before_trimming=30	key=REL606__1009100__-1__REL606__3511247__-1__4____36__36__0__0	max_left=13	max_left_minus=9	max_left_plus=13	max_min_left=13	max_min_left_minus=9	max_min_left_plus=13	max_min_right=0	max_min_right_minus=0	max_min_right_plus=0	max_pos_hash_score=60	max_right=28	max_right_minus=28	max_right_plus=23	neg_log10_pos_hash_p_value=0.4	new_junction_coverage=0.71	new_junction_read_count=37	new_junction_reference_weighted_read_count=0.00	new_junction_weighted_read_count=37.00	pos_hash_score=8	prediction=consensus	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=0.75	side_1_gene_name=pepN	side_1_gene_position=coding (1388/2613 nt)	side_1_gene_product=aminopeptidase N	side_1_gene_strand=>	side_1_locus_tag=ECB_00936	side_1_overlap=4	side_1_possible_overlap_registers=32	side_1_possible_overlap_registers_before_trimming=34	side_1_read_count=50	side_1_redundant=0	side_1_weighted_read_count=50.00	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=0.00	side_2_gene_name=yhhX	side_2_gene_position=coding (1016/1038 nt)	side_2_gene_product=predicted oxidoreductase with NAD(P)-binding Rossmann-fold domain	side_2_gene_strand=<	side_2_locus_tag=ECB_03291	side_2_overlap=0	side_2_possible_overlap_registers=24	side_2_possible_overlap_registers_before_trimming=30	side_2_read_count=0	side_2_redundant=0	side_2_weighted_read_count=0.00	total_non_overlap_reads=48
JC	160	.	REL606	1193615	1	REL606	4561283	1	0	alignment_overlap=0	coverage_minus=20	coverage_plus=25	flanking_left=36	flanking_right=36	frequency=1.000e+00	frequency_lower=9.356e-01	frequency_upper=1.000e+00	junction_effective_depth=45.00	junction_mixture_iterations=1	junction_possible_overlap_registers=31	junction_possible_overlap_registers_before_trimming=34	key=REL606__4550679__1__REL606__4561283__1__0____36__36__1__0	max_left=30	max_left_minus=14	max_left_plus=30	max_min_left=14	max_min_left_minus=14	max_min_left_plus=0	max_min_right=18	max_min_right_minus=0	max_min_right_plus=18	max_pos_hash_score=68	max_right=30	max_right_minus=30	max_right_plus=18	neg_log10_pos_hash_p_value=0.5	new_junction_coverage=0.70	new_junction_read_count=45	new_junction_reference_weighted_read_count=0.00	new_junction_weighted_read_count=45.00	pos_hash_score=8	prediction=consensus	side_1_annotate_key=repeat	side_1_continuation=0	side_1_coverage=NA	side_1_gene_name=IS1	side_1_gene_position=noncoding (1/768 nt)	side_1_gene_product=repeat region	side_1_gene_strand=>	side_1_overlap=0	side_1_possible_overlap_registers=NA	side_1_possible_overlap_registers_before_trimming=NA	side_1_read_count=NA	side_1_redundant=1	side_1_weighted_read_count=NA	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=0.00	side_2_gene_name=hsdR	side_2_gene_position=coding (3454/3513 nt)	side_2_gene_product=endonuclease R	side_2_gene_strand=<	side_2_locus_tag=ECB_04216	side_2_overlap=0	side_2_possible_overlap_registers=32	side_2_possible_overlap_registers_before_trimming=34	side_2_read_count=0	side_2_redundant=0	side_2_weighted_read_count=0.00	total_non_overlap_reads=45
JC	161	.	REL606	1194382	-1	REL606	1428800	-1	0	alignment_overlap=0	coverage_minus=9	coverage_plus=76	flanking_left=36	flanking_right=36	frequency=9.426e-01	frequency_lower=8.844e-01	frequency_upper=9.767e-01	junction_effective_depth=90.00	junction_mixture_iterations=1	junction_possible_overlap_registers=30	junction_possible_overlap_registers_before_trimming=34	key=REL606__1323488__-1__REL606__1428800__-1__0____36__36__1__0	max_left=25	max_left_minus=25	max_left_plus=22	max_min_left=13	max_min_left_minus=13	max_min_left_plus=12	max_min_right=15	max_min_right_minus=11	max_min_right_plus=15	max_pos_hash_score=68	max_right=31	max_right_minus=31	max_right_plus=29	neg_log10_pos_hash_p_value=0.1	new_junction_coverage=1.36	new_junction_read_count=85	new_junction_reference_weighted_read_count=0.00	new_junction_weighted_read_count=85.00	pos_hash_score=16	prediction=consensus	side_1_annotate_key=repeat	side_1_continuation=0	side_1_coverage=NA	side_1_gene_name=IS1	side_1_gene_position=noncoding (768/768 nt)	side_1_gene_product=repeat region	side_1_gene_strand=>	side_1_overlap=0	side_1_possible_overlap_registers=NA	side_1_possible_overlap_registers_before_trimming=NA	side_1_read_count=NA	side_1_redundant=1	side_1_weighted_read_count=NA	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=0.08	side_2_gene_name=ompN/ydbK	side_2_gene_position=intergenic (-301/+66)	side_2_gene_product=outer membrane pore protein N, non-specific/fused predicted pyruvate-flavodoxin oxidoreductase: conserved protein/conserved protein/FeS binding protein	side_2_gene_strand=</<	side_2_locus_tag=ECB_01348/ECB_01349	side_2_overlap=0	side_2_possible_overlap_registers=29	side_2_possible_overlap_registers_before_trimming=34	side_2_read_count=5	side_2_redundant=0	side_2_weighted_read_count=5.00	total_non_overlap_reads=85
JC	162	.	REL606	1420712	-1	REL606	1111786	-1	0	alignment_overlap=1	coverage_minus=52	coverage_plus=26	flanking_left=36	flanking_right=36	frequency=9.515e-01	frequency_lower=8.912e-01	frequency_upper=9.836e-01	junction_effective_depth=80.00	junction_mixture_iterations=1	junction_possible_overlap_registers=31	junction_possible_overlap_registers_before_trimming=33	key=REL606__1420712__-1__REL606__1607917__1__1____36__36__0__1	max_left=32	max_left_minus=28	max_left_plus=32	max_min_left=14	max_min_left_minus=14	max_min_left_plus=11	max_min_right=16	max_min_right_minus=8	max_min_right_plus=16	max_pos_hash_score=66	max_right=28	max_right_minus=28	max_right_plus=24	neg_log10_pos_hash_p_value=0.2	new_junction_coverage=1.17	new_junction_read_count=76	new_junction_reference_weighted_read_count=0.00	new_junction_weighted_read_count=76.00	pos_hash_score=13	prediction=consensus	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=0.06	side_1_gene_name=ynaK/ECB_01341	side_1_gene_position=intergenic (+24/-9)	side_1_gene_product=conserved protein/hypothetical protein	side_1_gene_strand=>/>	side_1_locus_tag=ECB_01340/ECB_01341	side_1_overlap=1	side_1_possible_overlap_registers=32	side_1_possible_overlap_registers_before_trimming=34	side_1_read_count=4	side_1_redundant=0	side_1_weighted_read_count=4.00	side_2_annotate_key=repeat	side_2_continuation=0	side_2_coverage=NA	side_2_gene_name=insE-1/serX	side_2_gene_position=intergenic (-55/+347)	side_2_gene_product=putative transposase-related protein/tRNA-Ser	side_2_gene_strand=</<	side_2_locus_tag=ECB_01029/ECB_t00024	side_2_overlap=0	side_2_possible_overlap_registers=NA	side_2_possible_overlap_registers_before_trimming=NA	side_2_read_count=NA	side_2_redundant=1	side_2_weighted_read_count=NA	total_non_overlap_reads=78
JC	163	.	REL606	1480958	1	REL606	3923283	1	0	alignment_overlap=5	coverage_minus=5	coverage_plus=0	flanking_left=36	flanking_right=36	frequency=7.773e-02	frequency_lower=3.376e-02	frequency_upper=1.493e-01	junction_effective_depth=74.50	junction_mixture_iterations=1	junction_possible_overlap_registers=25	junction_possible_overlap_registers_before_trimming=29	key=REL606__1480958__1__REL606__3923278__1__5____36__36__0__0	max_left=23	max_left_minus=23	max_left_plus=0	max_min_left=0	max_min_left_minus=0	max_min_left_plus=0	max_min_right=8	max_min_right_minus=8	max_min_right_plus=0	max_pos_hash_score=58	max_right=8	max_right_minus=8	max_right_plus=0	neg_log10_pos_hash_p_value=1.2	new_junction_coverage=0.10	new_junction_read_count=5	new_junction_reference_weighted_read_count=0.00	new_junction_weighted_read_count=5.00	pos_hash_score=3	prediction=polymorphism	reject=FREQUENCY_CUTOFF	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=1.39	side_1_gene_name=ydcR	side_1_gene_position=coding (846/1407 nt)	side_1_gene_product=fused predicted DNA-binding transcriptional regulator/predicted amino transferase	side_1_gene_strand=>	side_1_locus_tag=ECB_01396	side_1_overlap=5	side_1_possible_overlap_registers=32	side_1_possible_overlap_registers_before_trimming=34	side_1_read_count=93	side_1_redundant=0	side_1_weighted_read_count=93.00	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=0.88	side_2_gene_name=rep	side_2_gene_position=coding (805/2022 nt)	side_2_gene_product=DNA helicase and single-stranded DNA-dependent ATPase	side_2_gene_strand=>	side_2_locus_tag=ECB_03656	side_2_overlap=0	side_2_possible_overlap_registers=25	side_2_possible_overlap_registers_before_trimming=29	side_2_read_count=46	side_2_redundant=0	side_2_weighted_read_count=46.00	total_non_overlap_reads=5
JC	164	.	REL606	1729054	-1	REL606	666130	-1	0	alignment_overlap=0	coverage_minus=24	coverage_plus=40	flanking_left=36	flanking_right=36	frequency=6.154e-01	frequency_lower=5.261e-01	frequency_upper=6.991e-01	junction_effective_depth=95.00	junction_mixture_iterations=1	junction_possible_overlap_registers=30	junction_possible_overlap_registers_before_trimming=34	key=REL606__1729054__-1__REL606__2774435__1__0____36__36__0__1	max_left=31	max_left_minus=30	max_left_plus=31	max_min_left=14	max_min_left_minus=14	max_min_left_plus=10	max_min_right=17	max_min_right_minus=11	max_min_right_plus=17	max_pos_hash_score=68	max_right=31	max_right_minus=22	max_right_plus=31	neg_log10_pos_hash_p_value=0.2	new_junction_coverage=0.91	new_junction_read_count=57	new_junction_reference_weighted_read_count=0.00	new_junction_weighted_read_count=57.00	pos_hash_score=13	prediction=consensus	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=0.57	side_1_gene_name=ydhV	side_1_gene_position=coding (2044/2103 nt)	side_1_gene_product=predicted oxidoreductase	side_1_gene_strand=<	side_1_locus_tag=ECB_01642	side_1_overlap=0	side_1_possible_overlap_registers=32	side_1_possible_overlap_registers_before_trimming=34	side_1_read_count=38	side_1_redundant=0	side_1_weighted_read_count=38.00	side_2_annotate_key=repeat	side_2_continuation=0	side_2_coverage=NA	side_2_gene_name=IS150	side_2_gene_position=noncoding (1443/1443 nt)	side_2_gene_product=repeat region	side_2_gene_strand=>	side_2_overlap=0	side_2_possible_overlap_registers=NA	side_2_possible_overlap_registers_before_trimming=NA	side_2_read_count=NA	side_2_redundant=1	side_2_weighted_read_count=NA	total_non_overlap_reads=64
JC	165	.	REL606	1776437	-1	REL606	666130	-1	0	alignment_overlap=0	coverage_minus=49	coverage_plus=84	flanking_left=36	flanking_right=36	frequency=8.625e-01	frequency_lower=7.939e-01	frequency_upper=9.146e-01	junction_effective_depth=103.00	junction_mixture_iterations=1	junction_possible_overlap_registers=29	junction_possible_overlap_registers_before_trimming=34	key=REL606__1776437__-1__REL606__2774435__1__0____36__36__0__1	max_left=35	max_left_minus=35	max_left_plus=31	max_min_left=11	max_min_left_minus=7	max_min_left_plus=11	max_min_right=17	max_min_right_minus=15	max_min_right_plus=17	max_pos_hash_score=68	max_right=32	max_right_minus=32	max_right_plus=32	neg_log10_pos_hash_p_value=0.1	new_junction_coverage=1.45	new_junction_read_count=88	new_junction_reference_weighted_read_count=0.00	new_junction_weighted_read_count=88.00	pos_hash_score=18	prediction=consensus	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=0.23	side_1_gene_name=pheS/pheM	side_1_gene_position=intergenic (-228/+56)	side_1_gene_product=phenylalanyl-tRNA synthetase alpha subunit/phenylalanyl-tRNA synthetase operon leader peptide	side_1_gene_strand=</<	side_1_locus_tag=ECB_01683/ECB_01684	side_1_overlap=0	side_1_possible_overlap_registers=31	side_1_possible_overlap_registers_before_trimming=34	side_1_read_count=15	side_1_redundant=0	side_1_weighted_read_count=15.00	side_2_annotate_key=repeat	side_2_continuation=0	side_2_coverage=NA	side_2_gene_name=IS150	side_2_gene_position=noncoding (1443/1443 nt)	side_2_gene_product=repeat region	side_2_gene_strand=>	side_2_overlap=0	side_2_possible_overlap_registers=NA	side_2_possible_overlap_registers_before_trimming=NA	side_2_read_count=NA	side_2_redundant=1	side_2_weighted_read_count=NA	total_non_overlap_reads=133
JC	166	.	REL606	2157883	-1	REL606	666130	-1	0	alignment_overlap=0	coverage_minus=20	coverage_plus=76	flanking_left=36	flanking_right=36	frequency=1.000e+00	frequency_lower=9.693e-01	frequency_upper=1.000e+00	junction_effective_depth=96.00	junction_mixture_iterations=1	junction_possible_overlap_registers=32	junction_possible_overlap_registers_before_trimming=34	key=REL606__2157883__-1__REL606__2774435__1__0____36__36__0__1	max_left=31	max_left_minus=31	max_left_plus=21	max_min_left=17	max_min_left_minus=17	max_min_left_plus=12	max_min_right=18	max_min_right_minus=11	max_min_right_plus=18	max_pos_hash_score=68	max_right=33	max_right_minus=24	max_right_plus=33	neg_log10_pos_hash_p_value=0.1	new_junction_coverage=1.44	new_junction_read_count=96	new_junction_reference_weighted_read_count=0.00	new_junction_weighted_read_count=96.00	pos_hash_score=18	prediction=consensus	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=0.00	side_1_gene_name=yehM	side_1_gene_position=coding (993/2280 nt)	side_1_gene_product=hypothetical protein	side_1_gene_strand=>	side_1_locus_tag=ECB_02049	side_1_overlap=0	side_1_possible_overlap_registers=29	side_1_possible_overlap_registers_before_trimming=34	side_1_read_count=0	side_1_redundant=0	side_1_weighted_read_count=0.00	side_2_annotate_key=repeat	side_2_continuation=0	side_2_coverage=NA	side_2_gene_name=IS150	side_2_gene_position=noncoding (1443/1443 nt)	side_2_gene_product=repeat region	side_2_gene_strand=>	side_2_overlap=0	side_2_possible_overlap_registers=NA	side_2_possible_overlap_registers_before_trimming=NA	side_2_read_count=NA	side_2_redundant=1	side_2_weighted_read_count=NA	total_non_overlap_reads=96
JC	167	.	REL606	2170662	1	REL606	664688	1	0	alignment_overlap=0	coverage_minus=13	coverage_plus=38	flanking_left=36	flanking_right=36	frequency=1.000e+00	frequency_lower=9.430e-01	frequency_upper=1.000e+00	junction_effective_depth=51.00	junction_mixture_iterations=1	junction_possible_overlap_registers=32	junction_possible_overlap_registers_before_trimming=34	key=REL606__2170662__1__REL606__2775877__-1__0____36__36__0__1	max_left=34	max_left_minus=16	max_left_plus=34	max_min_left=16	max_min_left_minus=16	max_min_left_plus=12	max_min_right=12	max_min_right_minus=0	max_min_right_plus=12	max_pos_hash_score=68	max_right=32	max_right_minus=20	max_right_plus=32	neg_log10_pos_hash_p_value=0.4	new_junction_coverage=0.76	new_junction_read_count=51	new_junction_reference_weighted_read_count=0.00	new_junction_weighted_read_count=51.00	pos_hash_score=10	prediction=consensus	read_count_offset=3	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=0.00	side_1_gene_name=yehZ	side_1_gene_position=coding (17/918 nt)	side_1_gene_product=predicted transporter subunit: periplasmic-binding component of ABC superfamily	side_1_gene_strand=<	side_1_locus_tag=ECB_02061	side_1_overlap=0	side_1_possible_overlap_registers=26	side_1_possible_overlap_registers_before_trimming=31	side_1_read_count=0	side_1_redundant=0	side_1_weighted_read_count=0.00	side_2_annotate_key=repeat	side_2_continuation=0	side_2_coverage=NA	side_2_gene_name=IS150	side_2_gene_position=noncoding (1/1443 nt)	side_2_gene_product=repeat region	side_2_gene_strand=>	side_2_overlap=0	side_2_possible_overlap_registers=NA	side_2_possible_overlap_registers_before_trimming=NA	side_2_read_count=NA	side_2_redundant=1	side_2_weighted_read_count=NA	total_non_overlap_reads=51
JC	168	.	REL606	2170664	-1	REL606	666130	-1	0	alignment_overlap=1	coverage_minus=36	coverage_plus=72	flanking_left=36	flanking_right=36	frequency=1.000e+00	frequency_lower=9.724e-01	frequency_upper=1.000e+00	junction_effective_depth=107.00	junction_mixture_iterations=1	junction_possible_overlap_registers=31	junction_possible_overlap_registers_before_trimming=33	key=REL606__2170665__-1__REL606__2774435__1__1____36__36__0__1	max_left=32	max_left_minus=30	max_left_plus=32	max_min_left=17	max_min_left_minus=17	max_min_left_plus=12	max_min_right=16	max_min_right_minus=11	max_min_right_plus=16	max_pos_hash_score=66	max_right=31	max_right_minus=24	max_right_plus=31	neg_log10_pos_hash_p_value=0.1	new_junction_coverage=1.65	new_junction_read_count=107	new_junction_reference_weighted_read_count=0.00	new_junction_weighted_read_count=107.00	pos_hash_score=17	prediction=consensus	read_count_offset=3	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=0.00	side_1_gene_name=yehZ	side_1_gene_position=coding (15/918 nt)	side_1_gene_product=predicted transporter subunit: periplasmic-binding component of ABC superfamily	side_1_gene_strand=<	side_1_locus_tag=ECB_02061	side_1_overlap=0	side_1_possible_overlap_registers=26	side_1_possible_overlap_registers_before_trimming=30	side_1_read_count=0	side_1_redundant=0	side_1_weighted_read_count=0.00	side_2_annotate_key=repeat	side_2_continuation=0	side_2_coverage=NA	side_2_gene_name=IS150	side_2_gene_position=noncoding (1443/1443 nt)	side_2_gene_product=repeat region	side_2_gene_strand=>	side_2_overlap=1	side_2_possible_overlap_registers=NA	side_2_possible_overlap_registers_before_trimming=NA	side_2_read_count=NA	side_2_redundant=1	side_2_weighted_read_count=NA	total_non_overlap_reads=108
JC	169	.	REL606	2175935	1	REL606	2178784	-1	0	alignment_overlap=3	coverage_minus=0	coverage_plus=9	flanking_left=36	flanking_right=36	frequency=1.228e-01	frequency_lower=6.498e-02	frequency_upper=2.058e-01	junction_effective_depth=71.50	junction_mixture_iterations=1	junction_possible_overlap_registers=29	junction_possible_overlap_registers_before_trimming=31	key=REL606__2175935__1__REL606__2178787__-1__3____36__36__0__0	max_left=20	max_left_minus=0	max_left_plus=20	max_min_left=0	max_min_left_minus=0	max_min_left_plus=0	max_min_right=13	max_min_right_minus=0	max_min_right_plus=13	max_pos_hash_score=62	max_right=13	max_right_minus=0	max_right_plus=13	neg_log10_pos_hash_p_value=1.3	new_junction_coverage=0.15	new_junction_read_count=9	new_junction_reference_weighted_read_count=0.00	new_junction_weighted_read_count=9.00	pos_hash_score=3	prediction=polymorphism	reject=FREQUENCY_CUTOFF	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=0.85	side_1_gene_name=pbpG	side_1_gene_position=coding (174/933 nt)	side_1_gene_product=D-alanyl-D-alanine endopeptidase	side_1_gene_strand=<	side_1_locus_tag=ECB_02064	side_1_overlap=3	side_1_possible_overlap_registers=27	side_1_possible_overlap_registers_before_trimming=34	side_1_read_count=48	side_1_redundant=0	side_1_weighted_read_count=48.00	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=1.27	side_2_gene_name=yohGH	side_2_gene_position=coding (1142/1365 nt)	side_2_gene_product=predicted outer membrane protein	side_2_gene_strand=<	side_2_locus_tag=ECB_02068	side_2_overlap=0	side_2_possible_overlap_registers=29	side_2_possible_overlap_registers_before_trimming=31	side_2_read_count=77	side_2_redundant=0	side_2_weighted_read_count=77.00	total_non_overlap_reads=9
JC	170	.	REL606	2482732	1	REL606	2485952	-1	0	alignment_overlap=0	coverage_minus=14	coverage_plus=0	flanking_left=36	flanking_right=36	frequency=2.165e-01	frequency_lower=1.352e-01	frequency_upper=3.186e-01	junction_effective_depth=63.50	junction_mixture_iterations=1	junction_possible_overlap_registers=31	junction_possible_overlap_registers_before_trimming=34	key=REL606__2482732__1__REL606__2485952__-1__0____36__36__0__0	max_left=19	max_left_minus=19	max_left_plus=0	max_min_left=0	max_min_left_minus=0	max_min_left_plus=0	max_min_right=18	max_min_right_minus=18	max_min_right_plus=0	max_pos_hash_score=68	max_right=18	max_right_minus=18	max_right_plus=0	neg_log10_pos_hash_p_value=1.4	new_junction_coverage=0.22	new_junction_read_count=14	new_junction_reference_weighted_read_count=0.00	new_junction_weighted_read_count=14.00	pos_hash_score=3	prediction=polymorphism	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=0.56	side_1_gene_name=yfeG/yffI	side_1_gene_position=intergenic (-31/+15)	side_1_gene_product=predicted DNA-binding transcriptional regulator/predicted carboxysome structural protein with predicted role in ethanolamine utilization	side_1_gene_strand=</<	side_1_locus_tag=ECB_02337/ECB_02338	side_1_overlap=0	side_1_possible_overlap_registers=29	side_1_possible_overlap_registers_before_trimming=34	side_1_read_count=34	side_1_redundant=0	side_1_weighted_read_count=34.00	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=1.00	side_2_gene_name=eutB	side_2_gene_position=coding (247/1362 nt)	side_2_gene_product=ethanolamine ammonia-lyase, large subunit, heavy chain	side_2_gene_strand=<	side_2_locus_tag=ECB_02341	side_2_overlap=0	side_2_possible_overlap_registers=31	side_2_possible_overlap_registers_before_trimming=34	side_2_read_count=65	side_2_redundant=0	side_2_weighted_read_count=65.00	total_non_overlap_reads=14
JC	171	.	REL606	2655783	-1	REL606	666130	-1	0	alignment_overlap=1	coverage_minus=39	coverage_plus=54	flanking_left=36	flanking_right=36	frequency=8.232e-01	frequency_lower=7.480e-01	frequency_upper=8.833e-01	junction_effective_depth=99.00	junction_mixture_iterations=2	junction_possible_overlap_registers=29	junction_possible_overlap_registers_before_trimming=33	key=REL606__2655784__-1__REL606__3894996__-1__1____36__36__0__1	max_left=33	max_left_minus=33	max_left_plus=24	max_min_left=11	max_min_left_minus=10	max_min_left_plus=11	max_min_right=16	max_min_right_minus=2	max_min_right_plus=16	max_pos_hash_score=66	max_right=32	max_right_minus=32	max_right_plus=31	neg_log10_pos_hash_p_value=0.1	new_junction_coverage=1.35	new_junction_read_count=82	new_junction_reference_weighted_read_count=0.00	new_junction_weighted_read_count=82.00	pos_hash_score=14	prediction=consensus	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=0.29	side_1_gene_name=yfiH	side_1_gene_position=coding (660/732 nt)	side_1_gene_product=hypothetical protein	side_1_gene_strand=<	side_1_locus_tag=ECB_02483	side_1_overlap=0	side_1_possible_overlap_registers=28	side_1_possible_overlap_registers_before_trimming=33	side_1_read_count=17	side_1_redundant=0	side_1_weighted_read_count=17.00	side_2_annotate_key=repeat	side_2_continuation=0	side_2_coverage=NA	side_2_gene_name=IS150	side_2_gene_position=noncoding (1443/1443 nt)	side_2_gene_product=repeat region	side_2_gene_strand=>	side_2_overlap=1	side_2_possible_overlap_registers=NA	side_2_possible_overlap_registers_before_trimming=NA	side_2_read_count=NA	side_2_redundant=1	side_2_weighted_read_count=NA	total_non_overlap_reads=93
JC	172	.	REL606	2714060	1	REL606	2717533	-1	0	alignment_overlap=0	coverage_minus=8	coverage_plus=0	flanking_left=36	flanking_right=36	frequency=1.113e-01	frequency_lower=5.729e-02	frequency_upper=1.904e-01	junction_effective_depth=74.00	junction_mixture_iterations=1	junction_possible_overlap_registers=30	junction_possible_overlap_registers_before_trimming=34	key=REL606__2714060__1__REL606__2717533__-1__0____36__36__0__0	max_left=20	max_left_minus=20	max_left_plus=0	max_min_left=0	max_min_left_minus=0	max_min_left_plus=0	max_min_right=17	max_min_right_minus=17	max_min_right_plus=0	max_pos_hash_score=68	max_right=17	max_right_minus=17	max_right_plus=0	neg_log10_pos_hash_p_value=1.4	new_junction_coverage=0.13	new_junction_read_count=8	new_junction_reference_weighted_read_count=0.00	new_junction_weighted_read_count=8.00	pos_hash_score=3	prediction=polymorphism	reject=FREQUENCY_CUTOFF	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=1.27	side_1_gene_name=csrA/alaS	side_1_gene_position=intergenic (-55/+180)	side_1_gene_product=carbon storage regulator/alanyl-tRNA synthetase	side_1_gene_strand=</<	side_1_locus_tag=ECB_02546/ECB_02547	side_1_overlap=0	side_1_possible_overlap_registers=31	side_1_possible_overlap_registers_before_trimming=34	side_1_read_count=82	side_1_redundant=0	side_1_weighted_read_count=82.00	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=0.77	side_2_gene_name=recX/recA	side_2_gene_position=intergenic (-35/+34)	side_2_gene_product=RecA regulator RecX/recombinase A	side_2_gene_strand=</<	side_2_locus_tag=ECB_02548/ECB_02549	side_2_overlap=0	side_2_possible_overlap_registers=31	side_2_possible_overlap_registers_before_trimming=34	side_2_read_count=50	side_2_redundant=0	side_2_weighted_read_count=50.00	total_non_overlap_reads=8
JC	173	.	REL606	3015773	-1	REL606	666130	-1	0	alignment_overlap=1	coverage_minus=15	coverage_plus=38	flanking_left=36	flanking_right=36	frequency=1.000e+00	frequency_lower=9.450e-01	frequency_upper=1.000e+00	junction_effective_depth=53.00	junction_mixture_iterations=1	junction_possible_overlap_registers=30	junction_possible_overlap_registers_before_trimming=33	key=REL606__3015774__-1__REL606__3894996__-1__1____36__36__0__1	max_left=29	max_left_minus=29	max_left_plus=22	max_min_left=16	max_min_left_minus=16	max_min_left_plus=10	max_min_right=15	max_min_right_minus=9	max_min_right_plus=15	max_pos_hash_score=66	max_right=29	max_right_minus=17	max_right_plus=29	neg_log10_pos_hash_p_value=0.3	new_junction_coverage=0.85	new_junction_read_count=53	new_junction_reference_weighted_read_count=0.00	new_junction_weighted_read_count=53.00	pos_hash_score=10	prediction=consensus	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=0.00	side_1_gene_name=ECB_02816	side_1_gene_position=coding (1194/1677 nt)	side_1_gene_product=KpsD protein	side_1_gene_strand=>	side_1_locus_tag=ECB_02816	side_1_overlap=0	side_1_possible_overlap_registers=30	side_1_possible_overlap_registers_before_trimming=33	side_1_read_count=0	side_1_redundant=0	side_1_weighted_read_count=0.00	side_2_annotate_key=repeat	side_2_continuation=0	side_2_coverage=NA	side_2_gene_name=IS150	side_2_gene_position=noncoding (1443/1443 nt)	side_2_gene_product=repeat region	side_2_gene_strand=>	side_2_overlap=1	side_2_possible_overlap_registers=NA	side_2_possible_overlap_registers_before_trimming=NA	side_2_read_count=NA	side_2_redundant=1	side_2_weighted_read_count=NA	total_non_overlap_reads=53
JC	174	.	REL606	3260571	1	REL606	3260564	-1	0	alignment_overlap=7	coverage_minus=21	coverage_plus=36	flanking_left=36	flanking_right=36	frequency=1.000e+00	frequency_lower=9.019e-01	frequency_upper=1.000e+00	junction_effective_depth=29.00	junction_mixture_iterations=1	junction_possible_overlap_registers=18	junction_possible_overlap_registers_before_trimming=27	key=REL606__3260571__1__REL606__3260571__-1__7____36__36__0__0	max_left=18	max_left_minus=18	max_left_plus=13	max_min_left=13	max_min_left_minus=0	max_min_left_plus=13	max_min_right=12	max_min_right_minus=12	max_min_right_plus=0	max_pos_hash_score=54	max_right=21	max_right_minus=12	max_right_plus=21	neg_log10_pos_hash_p_value=0.3	new_junction_coverage=0.77	new_junction_read_count=29	new_junction_reference_weighted_read_count=0.00	new_junction_weighted_read_count=29.00	pos_hash_score=8	prediction=consensus	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=0.00	side_1_gene_name=hflB	side_1_gene_position=coding (1599/1935 nt)	side_1_gene_product=protease, ATP-dependent zinc-metallo	side_1_gene_strand=<	side_1_locus_tag=ECB_03043	side_1_overlap=7	side_1_possible_overlap_registers=21	side_1_possible_overlap_registers_before_trimming=34	side_1_read_count=0	side_1_redundant=0	side_1_weighted_read_count=0.00	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=0.00	side_2_gene_name=hflB	side_2_gene_position=coding (1606/1935 nt)	side_2_gene_product=protease, ATP-dependent zinc-metallo	side_2_gene_strand=<	side_2_locus_tag=ECB_03043	side_2_overlap=0	side_2_possible_overlap_registers=20	side_2_possible_overlap_registers_before_trimming=27	side_2_read_count=0	side_2_redundant=0	side_2_weighted_read_count=0.00	total_non_overlap_reads=57
JC	175	.	REL606	3511245	1	REL606	666130	-1	0	alignment_overlap=1	coverage_minus=14	coverage_plus=44	flanking_left=36	flanking_right=36	frequency=1.000e+00	frequency_lower=9.497e-01	frequency_upper=1.000e+00	junction_effective_depth=58.00	junction_mixture_iterations=1	junction_possible_overlap_registers=29	junction_possible_overlap_registers_before_trimming=33	key=REL606__3511244__1__REL606__3652533__-1__1____36__36__0__1	max_left=31	max_left_minus=30	max_left_plus=31	max_min_left=11	max_min_left_minus=4	max_min_left_plus=11	max_min_right=17	max_min_right_minus=17	max_min_right_plus=16	max_pos_hash_score=66	max_right=29	max_right_minus=29	max_right_plus=29	neg_log10_pos_hash_p_value=0.1	new_junction_coverage=0.96	new_junction_read_count=58	new_junction_reference_weighted_read_count=0.00	new_junction_weighted_read_count=58.00	pos_hash_score=16	prediction=consensus	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=0.00	side_1_gene_name=yhhX	side_1_gene_position=coding (1014/1038 nt)	side_1_gene_product=predicted oxidoreductase with NAD(P)-binding Rossmann-fold domain	side_1_gene_strand=<	side_1_locus_tag=ECB_03291	side_1_overlap=0	side_1_possible_overlap_registers=26	side_1_possible_overlap_registers_before_trimming=33	side_1_read_count=0	side_1_redundant=0	side_1_weighted_read_count=0.00	side_2_annotate_key=repeat	side_2_continuation=0	side_2_coverage=NA	side_2_gene_name=IS150	side_2_gene_position=noncoding (1443/1443 nt)	side_2_gene_product=repeat region	side_2_gene_strand=>	side_2_overlap=1	side_2_possible_overlap_registers=NA	side_2_possible_overlap_registers_before_trimming=NA	side_2_read_count=NA	side_2_redundant=1	side_2_weighted_read_count=NA	total_non_overlap_reads=58
JC	176	.	REL606	3630716	-1	REL606	666130	-1	0	alignment_overlap=0	coverage_minus=27	coverage_plus=54	flanking_left=36	flanking_right=36	frequency=NA	frequency_lower=NA	frequency_upper=NA	junction_effective_depth=69.00	junction_possible_overlap_registers=29	junction_possible_overlap_registers_before_trimming=34	key=REL606__3630716__-1__REL606__3652533__-1__0____36__36__1__1	max_left=32	max_left_minus=15	max_left_plus=32	max_min_left=15	max_min_left_minus=15	max_min_left_plus=14	max_min_right=17	max_min_right_minus=0	max_min_right_plus=17	max_pos_hash_score=68	max_right=31	max_right_minus=28	max_right_plus=31	neg_log10_pos_hash_p_value=0.1	new_junction_coverage=1.14	new_junction_read_count=69	new_junction_reference_weighted_read_count=0.03	new_junction_weighted_read_count=68.97	pos_hash_score=18	prediction=unknown	side_1_annotate_key=repeat	side_1_continuation=0	side_1_coverage=NA	side_1_gene_name=ldrD/yhjV	side_1_gene_position=intergenic (-171/-305)	side_1_gene_product=toxic polypeptide, small/predicted transporter	side_1_gene_strand=</>	side_1_locus_tag=ECB_03389/ECB_03390	side_1_overlap=0	side_1_possible_overlap_registers=NA	side_1_possible_overlap_registers_before_trimming=NA	side_1_read_count=NA	side_1_redundant=1	side_1_weighted_read_count=NA	side_2_annotate_key=repeat	side_2_continuation=0	side_2_coverage=NA	side_2_gene_name=IS150	side_2_gene_position=noncoding (1443/1443 nt)	side_2_gene_product=repeat region	side_2_gene_strand=>	side_2_overlap=0	side_2_possible_overlap_registers=NA	side_2_possible_overlap_registers_before_trimming=NA	side_2_read_count=NA	side_2_redundant=1	side_2_weighted_read_count=NA	total_non_overlap_reads=81
JC	177	.	REL606	3728853	1	REL606	3731968	-1	0	alignment_overlap=1	coverage_minus=0	coverage_plus=7	flanking_left=36	flanking_right=36	frequency=1.627e-01	frequency_lower=7.976e-02	frequency_upper=2.824e-01	junction_effective_depth=44.00	junction_mixture_iterations=1	junction_possible_overlap_registers=29	junction_possible_overlap_registers_before_trimming=33	key=REL606__3728853__1__REL606__3731969__-1__1____36__36__0__0	max_left=19	max_left_minus=0	max_left_plus=19	max_min_left=0	max_min_left_minus=0	max_min_left_plus=0	max_min_right=16	max_min_right_minus=0	max_min_right_plus=16	max_pos_hash_score=66	max_right=16	max_right_minus=0	max_right_plus=16	neg_log10_pos_hash_p_value=1.3	new_junction_coverage=0.12	new_junction_read_count=7	new_junction_reference_weighted_read_count=0.00	new_junction_weighted_read_count=7.00	pos_hash_score=3	prediction=polymorphism	reject=FREQUENCY_CUTOFF	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=0.31	side_1_gene_name=yibP	side_1_gene_position=coding (1250/1260 nt)	side_1_gene_product=protease with a role in cell division	side_1_gene_strand=>	side_1_locus_tag=ECB_03471	side_1_overlap=1	side_1_possible_overlap_registers=32	side_1_possible_overlap_registers_before_trimming=34	side_1_read_count=21	side_1_redundant=0	side_1_weighted_read_count=21.00	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=0.88	side_2_gene_name=tdh	side_2_gene_position=coding (144/1026 nt)	side_2_gene_product=L-threonine 3-dehydrogenase	side_2_gene_strand=<	side_2_locus_tag=ECB_03474	side_2_overlap=0	side_2_possible_overlap_registers=29	side_2_possible_overlap_registers_before_trimming=33	side_2_read_count=53	side_2_redundant=0	side_2_weighted_read_count=53.00	total_non_overlap_reads=7
JC	178	.	REL606	4388552	1	REL606	4391717	-1	0	alignment_overlap=0	coverage_minus=0	coverage_plus=7	flanking_left=36	flanking_right=36	frequency=8.244e-02	frequency_lower=3.836e-02	frequency_upper=1.516e-01	junction_effective_depth=80.50	junction_mixture_iterations=1	junction_possible_overlap_registers=31	junction_possible_overlap_registers_before_trimming=34	key=REL606__4388552__1__REL606__4391717__-1__0____36__36__0__0	max_left=20	max_left_minus=0	max_left_plus=20	max_min_left=0	max_min_left_minus=0	max_min_left_plus=0	max_min_right=17	max_min_right_minus=0	max_min_right_plus=17	max_pos_hash_score=68	max_right=17	max_right_minus=0	max_right_plus=17	neg_log10_pos_hash_p_value=1.4	new_junction_coverage=0.11	new_junction_read_count=7	new_junction_reference_weighted_read_count=0.00	new_junction_weighted_read_count=7.00	pos_hash_score=3	prediction=polymorphism	reject=FREQUENCY_CUTOFF	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=0.75	side_1_gene_name=yjfI	side_1_gene_position=coding (265/402 nt)	side_1_gene_product=hypothetical protein	side_1_gene_strand=>	side_1_locus_tag=ECB_04048	side_1_overlap=0	side_1_possible_overlap_registers=32	side_1_possible_overlap_registers_before_trimming=34	side_1_read_count=50	side_1_redundant=0	side_1_weighted_read_count=50.00	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=1.66	side_2_gene_name=yjfC	side_2_gene_position=coding (535/1164 nt)	side_2_gene_product=predicted synthetase/amidase	side_2_gene_strand=>	side_2_locus_tag=ECB_04053	side_2_overlap=0	side_2_possible_overlap_registers=28	side_2_possible_overlap_registers_before_trimming=34	side_2_read_count=97	side_2_redundant=0	side_2_weighted_read_count=97.00	total_non_overlap_reads=7
CN	179	.	REL606	547701	550400	0	gene_name=[nmpC]–ybcU	gene_product=[nmpC],ybcR,ybcS,ybcT,ybcU	locus_tag=[ECB_00503]–[ECB_00509]	relative_coverage=0	tile_size=100
CN	180	.	REL606	1609201	1615400	0	gene_name=ybcW–ECB_01523	gene_product=ybcW,gnsB,ynfN,ECB_01516,cspI,ydfP,ydfQ,ydfR,essQ,ECB_01522,ECB_01523	locus_tag=[ECB_01513]–[ECB_01523]	relative_coverage=0	tile_size=100
CN	181	.	REL606	2032801	2055500	0	gene_name=manC–[cpsG]	gene_product=manC,insB-14,insA-14,wbbD,wbbC,wzy,wbbB,wbbA,vioB,vioA,wzx,rmlC,rfbA,rfbD,rfbB,galF,wcaM,wcaL,wcaK,wzxC,wcaJ,[cpsG]	locus_tag=[ECB_01933]–[ECB_01954]	relative_coverage=0	tile_size=100
CN	182	.	REL606	2087001	2122400	0	gene_name=yegM–ECB_02012	gene_product=yegM,yegN,yegO,yegB,baeS,baeR,yegP,yegQ,ogrK,yegZ,ECB_01989,ECB_01990,ECB_01991,ECB_01992,ECB_01993,ECB_01994,ECB_01995,ECB_01996,ECB_01997,ECB_01998,ECB_01999,ECB_02000,ECB_02001,ECB_02002,ECB_02003,ECB_02004,ECB_02005,ECB_02006,ECB_02007,ECB_02008,ECB_02009,ECB_02010,ECB_02011,ECB_02012	locus_tag=[ECB_01979]–[ECB_02012]	relative_coverage=0	tile_size=100
CN	183	.	REL606	3015801	3035100	0	gene_name=[ECB_02816]–[ECB_02836]	gene_product=[ECB_02816],ECB_02817,ECB_02818,ECB_02819,ECB_02820,ECB_02821,ECB_02822,insB-22,insA-22,ECB_02825,ECB_02826,ECB_02827,ECB_02828,yghD,yghE,ECB_02831,ECB_02832,ECB_02833,ECB_02834,ECB_02835,[ECB_02836]	locus_tag=[ECB_02816]–[ECB_02836]	relative_coverage=0	tile_size=100
CN	184	.	REL606	3742001	3742200	0	gene_name=waaT	gene_position=pseudogene (166-365/405 nt)	gene_product=UDP-galactose:(Glucosyl) LPS alpha1,2-galactosyltransferase WaaT	gene_strand=<	locus_tag=ECB_03483	relative_coverage=0.0631	tile_size=100
CN	185	.	REL606	3895001	3901500	0	gene_name=rbsD–[yieO]	gene_product=rbsD,rbsA,rbsC,rbsB,rbsK,rbsR,[yieO]	locus_tag=[ECB_03634]–[ECB_03640]	relative_coverage=0	tile_size=100
CN	186	.	REL606	4522401	4561300	0	gene_name=[fimB]–[hsdR]	gene_product=[fimB],fimE,fimA,fimI,fimC,fimD,fimF,fimG,fimH,gntP,uxuA,uxuB,uxuR,yjiC,yjiD,yjiE,iadA,yjiG,yjiH,kptA,yjiJ,yjiK,yjiL,yjiM,yjiN,yjiO,yjiPQ,insA-28,insB-28,yjiV,mcrC,mcrB,yjiW,hsdS,hsdM,[hsdR]	locus_tag=[ECB_04179]–[ECB_04216]	relative_coverage=0	tile_size=100
UN	187	.	REL606	274	283
UN	188	.	REL606	5245	5252
UN	189	.	REL606	5480	5480
UN	190	.	REL606	7954	7954
UN	191	.	REL606	15245	15246
UN	192	.	REL606	15248	15250
UN	193	.	REL606	15538	15551
UN	194	.	REL606	15788	15802
UN	195	.	REL606	15910	15911
UN	196	.	REL606	16144	16165
UN	197	.	REL606	16263	16264
UN	198	.	REL606	16300	16302
UN	199	.	REL606	16550	16556
UN	200	.	REL606	16677	16686
UN	201	.	REL606	16879	16884
UN	202	.	REL606	34913	34915
UN	203	.	REL606	37330	37332
UN	204	.	REL606	39693	39698
UN	205	.	REL606	53680	53694
UN	206	.	REL606	56880	56882
UN	207	.	REL606	56884	56886
UN	208	.	REL606	56890	56891
UN	209	.	REL606	61581	61582
UN	210	.	REL606	66488	66488
UN	211	.	REL606	69398	69401
UN	212	.	REL606	69430	69458
UN	213	.	REL606	69524	69576
UN	214	.	REL606	73898	73903
UN	215	.	REL606	77248	77251
UN	216	.	REL606	81779	81779
UN	217	.	REL606	81786	81786
UN	218	.	REL606	93684	93685
UN	219	.	REL606	95839	95839
UN	220	.	REL606	96552	96556
UN	221	.	REL606	104606	104607
UN	222	.	REL606	111735	111742
UN	223	.	REL606	112014	112018
UN	224	.	REL606	112020	112023
UN	225	.	REL606	122968	122968
UN	226	.	REL606	122972	122972
UN	227	.	REL606	122974	122980
UN	228	.	REL606	123835	123835
UN	229	.	REL606	123839	123839
UN	230	.	REL606	123842	123842
UN	231	.	REL606	123971	123971
UN	232	.	REL606	123975	123981
UN	233	.	REL606	124303	124310
UN	234	.	REL606	124376	124377
UN	235	.	REL606	125225	125230
UN	236	.	REL606	128619	128621
UN	237	.	REL606	130408	130414
UN	238	.	REL606	133933	133933
UN	239	.	REL606	135638	135638
UN	240	.	REL606	144346	144346
UN	241	.	REL606	144348	144349
UN	242	.	REL606	146558	146564
UN	243	.	REL606	146843	146843
UN	244	.	REL606	146845	146848
UN	245	.	REL606	152164	152167
UN	246	.	REL606	152492	152492
UN	247	.	REL606	152498	152498
UN	248	.	REL606	152500	152500
UN	249	.	REL606	152655	152655
UN	250	.	REL606	154055	154055
UN	251	.	REL606	157977	157983
UN	252	.	REL606	159097	159097
UN	253	.	REL606	161869	161874
UN	254	.	REL606	167426	167441
UN	255	.	REL606	168092	168092
UN	256	.	REL606	168154	168155
UN	257	.	REL606	169320	169327
UN	258	.	REL606	173100	173106
UN	259	.	REL606	179178	179178
UN	260	.	REL606	179184	179184
UN	261	.	REL606	179186	179186
UN	262	.	REL606	188112	188112
UN	263	.	REL606	196310	196314
UN	264	.	REL606	197575	197580
UN	265	.	REL606	211769	211769
UN	266	.	REL606	213485	213487
UN	267	.	REL606	214458	214469
UN	268	.	REL606	215239	215239
UN	269	.	REL606	216186	216186
UN	270	.	REL606	216188	216188
UN	271	.	REL606	226685	226693
UN	272	.	REL606	226814	226817
UN	273	.	REL606	226954	226975
UN	274	.	REL606	227192	227222
UN	275	.	REL606	227291	227333
UN	276	.	REL606	227359	227361
UN	277	.	REL606	227363	227367
UN	278	.	REL606	227499	227499
UN	279	.	REL606	227504	227513
UN	280	.	REL606	227768	227787
UN	281	.	REL606	227789	227795
UN	282	.	REL606	227889	227900
UN	283	.	REL606	228268	228274
UN	284	.	REL606	228369	228371
UN	285	.	REL606	228670	228673
UN	286	.	REL606	229874	229874
UN	287	.	REL606	230395	230397
UN	288	.	REL606	230467	230471
UN	289	.	REL606	230829	230829
UN	290	.	REL606	231223	231224
UN	291	.	REL606	231226	231230
UN	292	.	REL606	231233	231233
UN	293	.	REL606	231353	231367
UN	294	.	REL606	231517	231517
UN	295	.	REL606	231889	231898
UN	296	.	REL606	231903	231907
UN	297	.	REL606	234718	234718
UN	298	.	REL606	241753	241760
UN	299	.	REL606	241815	241815
UN	300	.	REL606	242193	242193
UN	301	.	REL606	250261	250262
UN	302	.	REL606	251528	251531
UN	303	.	REL606	258217	258218
UN	304	.	REL606	263233	263236
UN	305	.	REL606	265320	265332
UN	306	.	REL606	265607	265617
UN	307	.	REL606	266571	266572
UN	308	.	REL606	272703	272704
UN	309	.	REL606	280721	280721
UN	310	.	REL606	291667	291675
UN	311	.	REL606	297459	297472
UN	312	.	REL606	300323	300323
UN	313	.	REL606	302156	302162
UN	314	.	REL606	304909	304913
UN	315	.	REL606	309938	309938
UN	316	.	REL606	309942	309946
UN	317	.	REL606	312065	312065
UN	318	.	REL606	312126	312126
UN	319	.	REL606	312437	312552
UN	320	.	REL606	312661	312697
UN	321	.	REL606	316546	316546
UN	322	.	REL606	316548	316549
UN	323	.	REL606	321736	321738
UN	324	.	REL606	322400	322451
UN	325	.	REL606	322483	322547
UN	326	.	REL606	326685	326691
UN	327	.	REL606	326888	326888
UN	328	.	REL606	327279	327284
UN	329	.	REL606	328238	328257
UN	330	.	REL606	330057	330111
UN	331	.	REL606	330257	330311
UN	332	.	REL606	332927	332931
UN	333	.	REL606	346904	346908
UN	334	.	REL606	349755	349756
UN	335	.	REL606	351359	351362
UN	336	.	REL606	357804	357805
UN	337	.	REL606	370826	370826
UN	338	.	REL606	376992	377021
UN	339	.	REL606	380996	380997
UN	340	.	REL606	391296	391300
UN	341	.	REL606	391306	391310
UN	342	.	REL606	391444	391453
UN	343	.	REL606	393137	393139
UN	344	.	REL606	401793	401797
UN	345	.	REL606	402541	402544
UN	346	.	REL606	404216	404224
UN	347	.	REL606	417994	418001
UN	348	.	REL606	422767	422767
UN	349	.	REL606	423091	423091
UN	350	.	REL606	423093	423100
UN	351	.	REL606	434996	434996
UN	352	.	REL606	447063	447071
UN	353	.	REL606	453286	453290
UN	354	.	REL606	456246	456246
UN	355	.	REL606	460092	460102
UN	356	.	REL606	461980	461986
UN	357	.	REL606	463011	463029
UN	358	.	REL606	464735	464740
UN	359	.	REL606	465329	465350
UN	360	.	REL606	467550	467555
UN	361	.	REL606	469594	469604
UN	362	.	REL606	469942	469951
UN	363	.	REL606	471812	471812
UN	364	.	REL606	480763	480964
UN	365	.	REL606	483251	483251
UN	366	.	REL606	484074	484082
UN	367	.	REL606	489546	489547
UN	368	.	REL606	489549	489551
UN	369	.	REL606	493455	493460
UN	370	.	REL606	493464	493464
UN	371	.	REL606	495337	495337
UN	372	.	REL606	495456	495472
UN	373	.	REL606	495593	495602
UN	374	.	REL606	495625	495636
UN	375	.	REL606	495678	495696
UN	376	.	REL606	495735	495738
UN	377	.	REL606	495771	495796
UN	378	.	REL606	495841	495851
UN	379	.	REL606	495883	495893
UN	380	.	REL606	495973	495993
UN	381	.	REL606	496018	496021
UN	382	.	REL606	496106	496157
UN	383	.	REL606	496206	496214
UN	384	.	REL606	496401	496432
UN	385	.	REL606	496521	496522
UN	386	.	REL606	496769	496771
UN	387	.	REL606	496782	496785
UN	388	.	REL606	496853	496856
UN	389	.	REL606	496939	496949
UN	390	.	REL606	496953	496953
UN	391	.	REL606	496956	496969
UN	392	.	REL606	497021	497036
UN	393	.	REL606	497433	497433
UN	394	.	REL606	498182	498203
UN	395	.	REL606	498483	498486
UN	396	.	REL606	498581	498581
UN	397	.	REL606	498741	498751
UN	398	.	REL606	498801	498811
UN	399	.	REL606	504089	504089
UN	400	.	REL606	505946	505946
UN	401	.	REL606	514477	514489
UN	402	.	REL606	517557	517563
UN	403	.	REL606	517853	517854
UN	404	.	REL606	519982	519990
UN	405	.	REL606	523156	523156
UN	406	.	REL606	523158	523159
UN	407	.	REL606	527209	527215
UN	408	.	REL606	540575	540586
UN	409	.	REL606	541292	541296
UN	410	.	REL606	547438	547441
UN	411	.	REL606	547683	547684
UN	412	.	REL606	547686	547687
UN	413	.	REL606	547689	550352
UN	414	.	REL606	551247	551259
UN	415	.	REL606	558068	558077
UN	416	.	REL606	569234	569241
UN	417	.	REL606	572335	572335
UN	418	.	REL606	572617	572622
UN	419	.	REL606	572728	572736
UN	420	.	REL606	572797	572807
UN	421	.	REL606	573199	573201
UN	422	.	REL606	575010	575010
UN	423	.	REL606	582060	582061
UN	424	.	REL606	582063	582063
UN	425	.	REL606	582065	582066
UN	426	.	REL606	589903	589904
UN	427	.	REL606	600633	600636
UN	428	.	REL606	601521	601524
UN	429	.	REL606	603318	603324
UN	430	.	REL606	607326	607326
UN	431	.	REL606	612713	612719
UN	432	.	REL606	616027	616027
UN	433	.	REL606	618556	618563
UN	434	.	REL606	619555	619555
UN	435	.	REL606	619654	619657
UN	436	.	REL606	619659	619663
UN	437	.	REL606	619666	619666
UN	438	.	REL606	625044	625044
UN	439	.	REL606	627114	627116
UN	440	.	REL606	627348	627353
UN	441	.	REL606	628722	628722
UN	442	.	REL606	629973	629982
UN	443	.	REL606	632816	632816
UN	444	.	REL606	638295	638295
UN	445	.	REL606	643365	643368
UN	446	.	REL606	643371	643371
UN	447	.	REL606	648058	648069
UN	448	.	REL606	653558	653558
UN	449	.	REL606	654140	654141
UN	450	.	REL606	658128	658128
UN	451	.	REL606	659719	659721
UN	452	.	REL606	659723	659735
UN	453	.	REL606	664967	664967
UN	454	.	REL606	664972	664972
UN	455	.	REL606	665277	665283
UN	456	.	REL606	676516	676517
UN	457	.	REL606	678554	678575
UN	458	.	REL606	678662	678682
UN	459	.	REL606	678789	678810
UN	460	.	REL606	678880	678901
UN	461	.	REL606	678989	679011
UN	462	.	REL606	679181	679202
UN	463	.	REL606	689209	689209
UN	464	.	REL606	694038	694045
UN	465	.	REL606	703260	703261
UN	466	.	REL606	711472	711478
UN	467	.	REL606	711859	711874
UN	468	.	REL606	711876	711885
UN	469	.	REL606	712048	712072
UN	470	.	REL606	712105	712135
UN	471	.	REL606	712193	712194
UN	472	.	REL606	712227	712271
UN	473	.	REL606	712323	712362
UN	474	.	REL606	712403	712437
UN	475	.	REL606	712543	712549
UN	476	.	REL606	712587	712646
UN	477	.	REL606	712676	712677
UN	478	.	REL606	712736	712768
UN	479	.	REL606	712867	712948
UN	480	.	REL606	713921	713943
UN	481	.	REL606	713976	714003
UN	482	.	REL606	714064	714065
UN	483	.	REL606	714127	714141
UN	484	.	REL606	714195	714210
UN	485	.	REL606	714217	714231
UN	486	.	REL606	714276	714298
UN	487	.	REL606	714332	714339
UN	488	.	REL606	714448	714556
UN	489	.	REL606	714603	714634
UN	490	.	REL606	714727	714759
UN	491	.	REL606	714788	714813
UN	492	.	REL606	715400	715407
UN	493	.	REL606	717914	717934
UN	494	.	REL606	718115	718122
UN	495	.	REL606	719178	719194
UN	496	.	REL606	719382	719399
UN	497	.	REL606	719685	719686
UN	498	.	REL606	720786	720786
UN	499	.	REL606	735940	735943
UN	500	.	REL606	736307	736309
UN	501	.	REL606	737311	737315
UN	502	.	REL606	739240	739243
UN	503	.	REL606	740388	740391
UN	504	.	REL606	755343	755345
UN	505	.	REL606	755874	755880
UN	506	.	REL606	761461	761466
UN	507	.	REL606	761468	761485
UN	508	.	REL606	761671	761692
UN	509	.	REL606	761745	761768
UN	510	.	REL606	761871	761906
UN	511	.	REL606	761956	761974
UN	512	.	REL606	762179	762201
UN	513	.	REL606	762386	762406
UN	514	.	REL606	764665	764666
UN	515	.	REL606	765460	765460
UN	516	.	REL606	790432	790433
UN	517	.	REL606	802418	802418
UN	518	.	REL606	808214	808215
UN	519	.	REL606	812773	812775
UN	520	.	REL606	814922	814927
UN	521	.	REL606	816089	816091
UN	522	.	REL606	826103	826108
UN	523	.	REL606	826221	826221
UN	524	.	REL606	826223	826224
UN	525	.	REL606	833662	833671
UN	526	.	REL606	836635	836636
UN	527	.	REL606	839268	839271
UN	528	.	REL606	839812	839818
UN	529	.	REL606	846942	846942
UN	530	.	REL606	848479	848482
UN	531	.	REL606	850055	850057
UN	532	.	REL606	850215	850227
UN	533	.	REL606	851747	851747
UN	534	.	REL606	865232	865235
UN	535	.	REL606	865842	865842
UN	536	.	REL606	866264	866264
UN	537	.	REL606	872424	872429
UN	538	.	REL606	880730	880733
UN	539	.	REL606	883394	883403
UN	540	.	REL606	888295	888301
UN	541	.	REL606	888303	888303
UN	542	.	REL606	888305	888305
UN	543	.	REL606	888309	888309
UN	544	.	REL606	907450	907450
UN	545	.	REL606	909162	909162
UN	546	.	REL606	913794	913801
UN	547	.	REL606	913959	913961
UN	548	.	REL606	914341	914347
UN	549	.	REL606	918099	918099
UN	550	.	REL606	918109	918110
UN	551	.	REL606	918329	918329
UN	552	.	REL606	918390	918390
UN	553	.	REL606	924282	924293
UN	554	.	REL606	925092	925094
UN	555	.	REL606	928061	928062
UN	556	.	REL606	929335	929341
UN	557	.	REL606	931793	931810
UN	558	.	REL606	931864	931872
UN	559	.	REL606	931911	931922
UN	560	.	REL606	935370	935376
UN	561	.	REL606	942562	942680
UN	562	.	REL606	942735	942861
UN	563	.	REL606	942912	942943
UN	564	.	REL606	950979	950980
UN	565	.	REL606	952703	952709
UN	566	.	REL606	953833	953833
UN	567	.	REL606	953835	953835
UN	568	.	REL606	967413	967413
UN	569	.	REL606	967415	967415
UN	570	.	REL606	985783	985783
UN	571	.	REL606	988923	988923
UN	572	.	REL606	991165	991167
UN	573	.	REL606	999060	999060
UN	574	.	REL606	1019998	1019998
UN	575	.	REL606	1020000	1020000
UN	576	.	REL606	1025263	1025274
UN	577	.	REL606	1025524	1025530
UN	578	.	REL606	1025893	1025893
UN	579	.	REL606	1029733	1029745
UN	580	.	REL606	1030388	1030393
UN	581	.	REL606	1036729	1036730
UN	582	.	REL606	1047817	1047817
UN	583	.	REL606	1048282	1048283
UN	584	.	REL606	1049935	1049937
UN	585	.	REL606	1049941	1049943
UN	586	.	REL606	1062335	1062335
UN	587	.	REL606	1072043	1072045
UN	588	.	REL606	1075702	1075704
UN	589	.	REL606	1082312	1082314
UN	590	.	REL606	1084776	1084781
UN	591	.	REL606	1084783	1084783
UN	592	.	REL606	1085145	1085146
UN	593	.	REL606	1089538	1089538
UN	594	.	REL606	1089542	1089542
UN	595	.	REL606	1090044	1090045
UN	596	.	REL606	1093283	1093286
UN	597	.	REL606	1099873	1099874
UN	598	.	REL606	1100694	1100699
UN	599	.	REL606	1102959	1102962
UN	600	.	REL606	1111282	1111287
UN	601	.	REL606	1111813	1111838
UN	602	.	REL606	1111910	1111917
UN	603	.	REL606	1111960	1111966
UN	604	.	REL606	1111994	1112029
UN	605	.	REL606	1113371	1113371
UN	606	.	REL606	1114527	1114536
UN	607	.	REL606	1128101	1128105
UN	608	.	REL606	1129799	1129799
UN	609	.	REL606	1130119	1130122
UN	610	.	REL606	1139813	1139821
UN	611	.	REL606	1142646	1142654
UN	612	.	REL606	1144604	1144605
UN	613	.	REL606	1146779	1146779
UN	614	.	REL606	1165636	1165636
UN	615	.	REL606	1165639	1165639
UN	616	.	REL606	1167551	1167551
UN	617	.	REL606	1167553	1167554
UN	618	.	REL606	1167556	1167565
UN	619	.	REL606	1176821	1176833
UN	620	.	REL606	1177331	1177340
UN	621	.	REL606	1177400	1177410
UN	622	.	REL606	1182403	1182405
UN	623	.	REL606	1184648	1184649
UN	624	.	REL606	1185099	1185103
UN	625	.	REL606	1185424	1185439
UN	626	.	REL606	1188059	1188060
UN	627	.	REL606	1188063	1188069
UN	628	.	REL606	1192807	1192811
UN	629	.	REL606	1194109	1194115
UN	630	.	REL606	1196404	1196404
UN	631	.	REL606	1200831	1200832
UN	632	.	REL606	1207587	1207598
UN	633	.	REL606	1210183	1210188
UN	634	.	REL606	1210193	1210194
UN	635	.	REL606	1218265	1218265
UN	636	.	REL606	1219357	1219363
UN	637	.	REL606	1222784	1222785
UN	638	.	REL606	1222787	1222792
UN	639	.	REL606	1228734	1228738
UN	640	.	REL606	1228741	1228741
UN	641	.	REL606	1228744	1228744
UN	642	.	REL606	1247402	1247405
UN	643	.	REL606	1254350	1254360
UN	644	.	REL606	1263505	1263505
UN	645	.	REL606	1263509	1263509
UN	646	.	REL606	1269250	1269361
UN	647	.	REL606	1269406	1269443
UN	648	.	REL606	1269554	1269581
UN	649	.	REL606	1269788	1269885
UN	650	.	REL606	1269945	1270000
UN	651	.	REL606	1270054	1270663
UN	652	.	REL606	1270724	1270739
UN	653	.	REL606	1270853	1270862
UN	654	.	REL606	1279430	1279434
UN	655	.	REL606	1280228	1280228
UN	656	.	REL606	1280233	1280238
UN	657	.	REL606	1287067	1287093
UN	658	.	REL606	1287361	1287388
UN	659	.	REL606	1303010	1303016
UN	660	.	REL606	1309413	1309416
UN	661	.	REL606	1310740	1310744
UN	662	.	REL606	1315348	1315349
UN	663	.	REL606	1315432	1315435
UN	664	.	REL606	1316063	1316063
UN	665	.	REL606	1318184	1318193
UN	666	.	REL606	1320103	1320113
UN	667	.	REL606	1326066	1326070
UN	668	.	REL606	1328140	1328140
UN	669	.	REL606	1333962	1333967
UN	670	.	REL606	1334473	1334475
UN	671	.	REL606	1336302	1336305
UN	672	.	REL606	1340295	1340299
UN	673	.	REL606	1341166	1341170
UN	674	.	REL606	1349559	1349561
UN	675	.	REL606	1349564	1349564
UN	676	.	REL606	1350064	1350064
UN	677	.	REL606	1350520	1350520
UN	678	.	REL606	1353176	1353182
UN	679	.	REL606	1353192	1353195
UN	680	.	REL606	1374774	1374781
UN	681	.	REL606	1375510	1375512
UN	682	.	REL606	1375514	1375516
UN	683	.	REL606	1397533	1397540
UN	684	.	REL606	1398085	1398088
UN	685	.	REL606	1402579	1402584
UN	686	.	REL606	1415909	1415912
UN	687	.	REL606	1418078	1418086
UN	688	.	REL606	1419877	1419877
UN	689	.	REL606	1421993	1422004
UN	690	.	REL606	1422006	1422013
UN	691	.	REL606	1422043	1422051
UN	692	.	REL606	1422205	1422224
UN	693	.	REL606	1422226	1422226
UN	694	.	REL606	1422423	1422427
UN	695	.	REL606	1425991	1425993
UN	696	.	REL606	1436154	1436154
UN	697	.	REL606	1436156	1436159
UN	698	.	REL606	1454430	1454430
UN	699	.	REL606	1456486	1456486
UN	700	.	REL606	1456566	1456568
UN	701	.	REL606	1460278	1460278
UN	702	.	REL606	1460284	1460290
UN	703	.	REL606	1462254	1462326
UN	704	.	REL606	1478578	1478585
UN	705	.	REL606	1478587	1478588
UN	706	.	REL606	1478590	1478590
UN	707	.	REL606	1478593	1478593
UN	708	.	REL606	1490211	1490222
UN	709	.	REL606	1492600	1492601
UN	710	.	REL606	1498365	1498368
UN	711	.	REL606	1498398	1498440
UN	712	.	REL606	1499157	1499157
UN	713	.	REL606	1500331	1500332
UN	714	.	REL606	1500388	1500393
UN	715	.	REL606	1500442	1500451
UN	716	.	REL606	1500467	1500468
UN	717	.	REL606	1500532	1500534
UN	718	.	REL606	1500620	1500639
UN	719	.	REL606	1500701	1500701
UN	720	.	REL606	1500703	1500724
UN	721	.	REL606	1500733	1500734
UN	722	.	REL606	1500736	1500736
UN	723	.	REL606	1500831	1500837
UN	724	.	REL606	1500872	1500873
UN	725	.	REL606	1500876	1500903
UN	726	.	REL606	1500974	1500974
UN	727	.	REL606	1500976	1500996
UN	728	.	REL606	1501244	1501266
UN	729	.	REL606	1501288	1501288
UN	730	.	REL606	1501628	1501648
UN	731	.	REL606	1501724	1501727
UN	732	.	REL606	1501803	1501805
UN	733	.	REL606	1501992	1501999
UN	734	.	REL606	1502966	1502966
UN	735	.	REL606	1503088	1503089
UN	736	.	REL606	1503203	1503218
UN	737	.	REL606	1503221	1503221
UN	738	.	REL606	1503615	1503620
UN	739	.	REL606	1507073	1507081
UN	740	.	REL606	1514337	1514351
UN	741	.	REL606	1526377	1526379
UN	742	.	REL606	1526381	1526382
UN	743	.	REL606	1540171	1540175
UN	744	.	REL606	1542867	1542868
UN	745	.	REL606	1544737	1544742
UN	746	.	REL606	1545715	1545728
UN	747	.	REL606	1549084	1549087
UN	748	.	REL606	1549089	1549089
UN	749	.	REL606	1556360	1556360
UN	750	.	REL606	1572152	1572154
UN	751	.	REL606	1600398	1600400
UN	752	.	REL606	1602489	1602493
UN	753	.	REL606	1605198	1605199
UN	754	.	REL606	1606855	1606866
UN	755	.	REL606	1606936	1606944
UN	756	.	REL606	1606948	1606961
UN	757	.	REL606	1607213	1607218
UN	758	.	REL606	1607397	1607398
UN	759	.	REL606	1607402	1607408
UN	760	.	REL606	1607734	1607747
UN	761	.	REL606	1607775	1607792
UN	762	.	REL606	1607874	1607898
UN	763	.	REL606	1607963	1607974
UN	764	.	REL606	1607995	1607995
UN	765	.	REL606	1608092	1608126
UN	766	.	REL606	1608156	1608161
UN	767	.	REL606	1608224	1608262
UN	768	.	REL606	1608348	1608348
UN	769	.	REL606	1608365	1608391
UN	770	.	REL606	1608505	1608509
UN	771	.	REL606	1608576	1608582
UN	772	.	REL606	1608596	1608597
UN	773	.	REL606	1608630	1608631
UN	774	.	REL606	1608643	1608647
UN	775	.	REL606	1608678	1608715
UN	776	.	REL606	1608719	1608720
UN	777	.	REL606	1608746	1608746
UN	778	.	REL606	1608748	1608748
UN	779	.	REL606	1608943	1608944
UN	780	.	REL606	1609049	1609067
UN	781	.	REL606	1609107	1609107
UN	782	.	REL606	1609109	1609138
UN	783	.	REL606	1609172	1615473
UN	784	.	REL606	1615618	1615619
UN	785	.	REL606	1615889	1615889
UN	786	.	REL606	1615971	1615977
UN	787	.	REL606	1616278	1616280
UN	788	.	REL606	1616283	1616283
UN	789	.	REL606	1616285	1616285
UN	790	.	REL606	1616287	1616287
UN	791	.	REL606	1616327	1616329
UN	792	.	REL606	1616343	1616350
UN	793	.	REL606	1616385	1616389
UN	794	.	REL606	1616642	1616642
UN	795	.	REL606	1616687	1616688
UN	796	.	REL606	1619587	1619587
UN	797	.	REL606	1619646	1619651
UN	798	.	REL606	1625003	1625003
UN	799	.	REL606	1626419	1626422
UN	800	.	REL606	1632921	1632927
UN	801	.	REL606	1634633	1634643
UN	802	.	REL606	1634775	1634775
UN	803	.	REL606	1634778	1634778
UN	804	.	REL606	1634781	1634781
UN	805	.	REL606	1636034	1636040
UN	806	.	REL606	1652860	1652860
UN	807	.	REL606	1655397	1655397
UN	808	.	REL606	1665560	1665561
UN	809	.	REL606	1668763	1668763
UN	810	.	REL606	1675336	1675343
UN	811	.	REL606	1682451	1682466
UN	812	.	REL606	1685887	1686075
UN	813	.	REL606	1688990	1688990
UN	814	.	REL606	1700396	1700396
UN	815	.	REL606	1710311	1710311
UN	816	.	REL606	1710827	1710834
UN	817	.	REL606	1710903	1710906
UN	818	.	REL606	1710911	1710914
UN	819	.	REL606	1716941	1716943
UN	820	.	REL606	1733139	1733148
UN	821	.	REL606	1755500	1755505
UN	822	.	REL606	1764888	1764890
UN	823	.	REL606	1769996	1769996
UN	824	.	REL606	1777753	1777754
UN	825	.	REL606	1788318	1788318
UN	826	.	REL606	1790686	1790687
UN	827	.	REL606	1810681	1810682
UN	828	.	REL606	1811325	1811332
UN	829	.	REL606	1815476	1815478
UN	830	.	REL606	1819579	1819585
UN	831	.	REL606	1836453	1836453
UN	832	.	REL606	1848198	1848199
UN	833	.	REL606	1851049	1851049
UN	834	.	REL606	1857235	1857235
UN	835	.	REL606	1860496	1860496
UN	836	.	REL606	1864076	1864083
UN	837	.	REL606	1865050	1865056
UN	838	.	REL606	1871662	1871663
UN	839	.	REL606	1872448	1872459
UN	840	.	REL606	1875006	1875006
UN	841	.	REL606	1879112	1879114
UN	842	.	REL606	1881603	1881609
UN	843	.	REL606	1883343	1883346
UN	844	.	REL606	1888146	1888153
UN	845	.	REL606	1893942	1893945
UN	846	.	REL606	1895194	1895202
UN	847	.	REL606	1896677	1896698
UN	848	.	REL606	1908944	1908948
UN	849	.	REL606	1918152	1918152
UN	850	.	REL606	1961624	1961624
UN	851	.	REL606	1984817	1984821
UN	852	.	REL606	1985248	1985261
UN	853	.	REL606	1991320	1991325
UN	854	.	REL606	2005385	2005390
UN	855	.	REL606	2006503	2006506
UN	856	.	REL606	2018411	2018414
UN	857	.	REL606	2019262	2019270
UN	858	.	REL606	2020124	2020129
UN	859	.	REL606	2020742	2020754
UN	860	.	REL606	2024366	2024374
UN	861	.	REL606	2031327	2031328
UN	862	.	REL606	2031961	2031977
UN	863	.	REL606	2032099	2032099
UN	864	.	REL606	2032347	2032348
UN	865	.	REL606	2032351	2032360
UN	866	.	REL606	2032709	2045242
UN	867	.	REL606	2045275	2045533
UN	868	.	REL606	2045543	2045543
UN	869	.	REL606	2045546	2045546
UN	870	.	REL606	2045565	2054822
UN	871	.	REL606	2054873	2054947
UN	872	.	REL606	2054978	2055045
UN	873	.	REL606	2055047	2055047
UN	874	.	REL606	2055084	2055124
UN	875	.	REL606	2055179	2055247
UN	876	.	REL606	2055276	2055307
UN	877	.	REL606	2055342	2055607
UN	878	.	REL606	2055646	2055781
UN	879	.	REL606	2055815	2055815
UN	880	.	REL606	2055818	2055903
UN	881	.	REL606	2055950	2055959
UN	882	.	REL606	2056469	2056469
UN	883	.	REL606	2058837	2058837
UN	884	.	REL606	2071997	2071998
UN	885	.	REL606	2077218	2077226
UN	886	.	REL606	2081511	2081514
UN	887	.	REL606	2081516	2081516
UN	888	.	REL606	2085980	2085985
UN	889	.	REL606	2086104	2086107
UN	890	.	REL606	2086229	2086231
UN	891	.	REL606	2086308	2086317
UN	892	.	REL606	2086512	2122428
UN	893	.	REL606	2122900	2122920
UN	894	.	REL606	2125178	2125183
UN	895	.	REL606	2131939	2131947
UN	896	.	REL606	2138815	2138820
UN	897	.	REL606	2143700	2143713
UN	898	.	REL606	2156987	2156994
UN	899	.	REL606	2157884	2157891
UN	900	.	REL606	2161799	2161805
UN	901	.	REL606	2165281	2165281
UN	902	.	REL606	2170250	2170252
UN	903	.	REL606	2173460	2173460
UN	904	.	REL606	2174160	2174164
UN	905	.	REL606	2180234	2180234
UN	906	.	REL606	2180236	2180236
UN	907	.	REL606	2188546	2188547
UN	908	.	REL606	2204484	2204499
UN	909	.	REL606	2207463	2207465
UN	910	.	REL606	2215477	2215477
UN	911	.	REL606	2215597	2215599
UN	912	.	REL606	2216083	2216086
UN	913	.	REL606	2229512	2229518
UN	914	.	REL606	2233420	2233420
UN	915	.	REL606	2240633	2240636
UN	916	.	REL606	2242418	2242422
UN	917	.	REL606	2245716	2245719
UN	918	.	REL606	2247097	2247103
UN	919	.	REL606	2247369	2247369
UN	920	.	REL606	2247499	2247505
UN	921	.	REL606	2252683	2252684
UN	922	.	REL606	2254468	2254638
UN	923	.	REL606	2254692	2254692
UN	924	.	REL606	2254694	2255049
UN	925	.	REL606	2256138	2256138
UN	926	.	REL606	2260782	2260782
UN	927	.	REL606	2262584	2262586
UN	928	.	REL606	2268783	2268788
UN	929	.	REL606	2269464	2269465
UN	930	.	REL606	2270463	2270468
UN	931	.	REL606	2270471	2270478
UN	932	.	REL606	2274829	2274835
UN	933	.	REL606	2279689	2279689
UN	934	.	REL606	2280059	2280067
UN	935	.	REL606	2281331	2281334
UN	936	.	REL606	2283303	2283303
UN	937	.	REL606	2283306	2283307
UN	938	.	REL606	2286830	2286839
UN	939	.	REL606	2289040	2289041
UN	940	.	REL606	2293925	2293930
UN	941	.	REL606	2293933	2293936
UN	942	.	REL606	2298310	2298310
UN	943	.	REL606	2306959	2306965
UN	944	.	REL606	2309627	2309631
UN	945	.	REL606	2318033	2318034
UN	946	.	REL606	2318038	2318038
UN	947	.	REL606	2324508	2324513
UN	948	.	REL606	2333113	2333113
UN	949	.	REL606	2335517	2335521
UN	950	.	REL606	2346932	2346938
UN	951	.	REL606	2367463	2367470
UN	952	.	REL606	2371369	2371369
UN	953	.	REL606	2375046	2375046
UN	954	.	REL606	2378330	2378337
UN	955	.	REL606	2382766	2382770
UN	956	.	REL606	2389208	2389208
UN	957	.	REL606	2390865	2390866
UN	958	.	REL606	2393519	2393519
UN	959	.	REL606	2393675	2393679
UN	960	.	REL606	2393681	2393682
UN	961	.	REL606	2395224	2395228
UN	962	.	REL606	2398666	2398669
UN	963	.	REL606	2403653	2403662
UN	964	.	REL606	2404173	2404175
UN	965	.	REL606	2404260	2404265
UN	966	.	REL606	2407936	2407936
UN	967	.	REL606	2407939	2407939
UN	968	.	REL606	2407941	2407951
UN	969	.	REL606	2410144	2410151
UN	970	.	REL606	2450942	2450973
UN	971	.	REL606	2451057	2451078
UN	972	.	REL606	2453831	2453859
UN	973	.	REL606	2453863	2453863
UN	974	.	REL606	2453952	2453979
UN	975	.	REL606	2454074	2454103
UN	976	.	REL606	2461637	2461641
UN	977	.	REL606	2472454	2472454
UN	978	.	REL606	2472456	2472456
UN	979	.	REL606	2472458	2472458
UN	980	.	REL606	2484102	2484102
UN	981	.	REL606	2484104	2484104
UN	982	.	REL606	2486515	2486515
UN	983	.	REL606	2486518	2486518
UN	984	.	REL606	2488882	2488918
UN	985	.	REL606	2488985	2489013
UN	986	.	REL606	2491101	2491102
UN	987	.	REL606	2491620	2491622
UN	988	.	REL606	2496797	2496812
UN	989	.	REL606	2498235	2498241
UN	990	.	REL606	2498757	2498761
UN	991	.	REL606	2502471	2502498
UN	992	.	REL606	2503338	2503338
UN	993	.	REL606	2503340	2503341
UN	994	.	REL606	2509561	2509561
UN	995	.	REL606	2525388	2525388
UN	996	.	REL606	2525978	2525983
UN	997	.	REL606	2529863	2529863
UN	998	.	REL606	2531911	2531916
UN	999	.	REL606	2538596	2538599
UN	1000	.	REL606	2540216	2540221
UN	1001	.	REL606	2552147	2552149
UN	1002	.	REL606	2553410	2553415
UN	1003	.	REL606	2560611	2560612
UN	1004	.	REL606	2560962	2560964
UN	1005	.	REL606	2563737	2563746
UN	1006	.	REL606	2570970	2570970
UN	1007	.	REL606	2571124	2571132
UN	1008	.	REL606	2573100	2573112
UN	1009	.	REL606	2583073	2583263
UN	1010	.	REL606	2583390	2583392
UN	1011	.	REL606	2584955	2584984
UN	1012	.	REL606	2597735	2597740
UN	1013	.	REL606	2612982	2612982
UN	1014	.	REL606	2612985	2612985
UN	1015	.	REL606	2612989	2612990
UN	1016	.	REL606	2612992	2612993
UN	1017	.	REL606	2612995	2612995
UN	1018	.	REL606	2612997	2613001
UN	1019	.	REL606	2614605	2614608
UN	1020	.	REL606	2615844	2615847
UN	1021	.	REL606	2618447	2618474
UN	1022	.	REL606	2628950	2628953
UN	1023	.	REL606	2632949	2632950
UN	1024	.	REL606	2642080	2642093
UN	1025	.	REL606	2647670	2647670
UN	1026	.	REL606	2647678	2647683
UN	1027	.	REL606	2647896	2647912
UN	1028	.	REL606	2648009	2648017
UN	1029	.	REL606	2648026	2648042
UN	1030	.	REL606	2648179	2648179
UN	1031	.	REL606	2648436	2648443
UN	1032	.	REL606	2648594	2648606
UN	1033	.	REL606	2649578	2649583
UN	1034	.	REL606	2650242	2650250
UN	1035	.	REL606	2651200	2651224
UN	1036	.	REL606	2651359	2651372
UN	1037	.	REL606	2651664	2651664
UN	1038	.	REL606	2651667	2651680
UN	1039	.	REL606	2651803	2651808
UN	1040	.	REL606	2651850	2651855
UN	1041	.	REL606	2652013	2652047
UN	1042	.	REL606	2652092	2652094
UN	1043	.	REL606	2652208	2652215
UN	1044	.	REL606	2652451	2652451
UN	1045	.	REL606	2652475	2652487
UN	1046	.	REL606	2652772	2652783
UN	1047	.	REL606	2658318	2658322
UN	1048	.	REL606	2669890	2669907
UN	1049	.	REL606	2670782	2670785
UN	1050	.	REL606	2685197	2685197
UN	1051	.	REL606	2685306	2685312
UN	1052	.	REL606	2686127	2686127
UN	1053	.	REL606	2689145	2689146
UN	1054	.	REL606	2693114	2693115
UN	1055	.	REL606	2693117	2693118
UN	1056	.	REL606	2693121	2693121
UN	1057	.	REL606	2693124	2693124
UN	1058	.	REL606	2697950	2697950
UN	1059	.	REL606	2698131	2698132
UN	1060	.	REL606	2698898	2698914
UN	1061	.	REL606	2703750	2703751
UN	1062	.	REL606	2706346	2706362
UN	1063	.	REL606	2710013	2710014
UN	1064	.	REL606	2710016	2710016
UN	1065	.	REL606	2710018	2710019
UN	1066	.	REL606	2712664	2712698
UN	1067	.	REL606	2712793	2712817
UN	1068	.	REL606	2712939	2712993
UN	1069	.	REL606	2713078	2713132
UN	1070	.	REL606	2713208	2713238
UN	1071	.	REL606	2713353	2713375
UN	1072	.	REL606	2721329	2721330
UN	1073	.	REL606	2721332	2721339
UN	1074	.	REL606	2721341	2721341
UN	1075	.	REL606	2724945	2724948
UN	1076	.	REL606	2729187	2729195
UN	1077	.	REL606	2738460	2738462
UN	1078	.	REL606	2738464	2738474
UN	1079	.	REL606	2738907	2738918
UN	1080	.	REL606	2743151	2743152
UN	1081	.	REL606	2743936	2743944
UN	1082	.	REL606	2751570	2751581
UN	1083	.	REL606	2756303	2756303
UN	1084	.	REL606	2760115	2760118
UN	1085	.	REL606	2760120	2760120
UN	1086	.	REL606	2760122	2760122
UN	1087	.	REL606	2760125	2760125
UN	1088	.	REL606	2764827	2764834
UN	1089	.	REL606	2769333	2769340
UN	1090	.	REL606	2770153	2770156
UN	1091	.	REL606	2770807	2770815
UN	1092	.	REL606	2773497	2773498
UN	1093	.	REL606	2773731	2773748
UN	1094	.	REL606	2773764	2773764
UN	1095	.	REL606	2774014	2774025
UN	1096	.	REL606	2778869	2778880
UN	1097	.	REL606	2779193	2779196
UN	1098	.	REL606	2786189	2786197
UN	1099	.	REL606	2790045	2790045
UN	1100	.	REL606	2792991	2793002
UN	1101	.	REL606	2800013	2800024
UN	1102	.	REL606	2800549	2800560
UN	1103	.	REL606	2810103	2810103
UN	1104	.	REL606	2810106	2810106
UN	1105	.	REL606	2810109	2810110
UN	1106	.	REL606	2810378	2810385
UN	1107	.	REL606	2811108	2811108
UN	1108	.	REL606	2815225	2815235
UN	1109	.	REL606	2822710	2822719
UN	1110	.	REL606	2826645	2826645
UN	1111	.	REL606	2828174	2828181
UN	1112	.	REL606	2828384	2828386
UN	1113	.	REL606	2842026	2842058
UN	1114	.	REL606	2842124	2842168
UN	1115	.	REL606	2842236	2842278
UN	1116	.	REL606	2843732	2843734
UN	1117	.	REL606	2844584	2844584
UN	1118	.	REL606	2845381	2845384
UN	1119	.	REL606	2845387	2845387
UN	1120	.	REL606	2845390	2845390
UN	1121	.	REL606	2845392	2845392
UN	1122	.	REL606	2846254	2846262
UN	1123	.	REL606	2848978	2848978
UN	1124	.	REL606	2858932	2858939
UN	1125	.	REL606	2875068	2875072
UN	1126	.	REL606	2879417	2879427
UN	1127	.	REL606	2899617	2899622
UN	1128	.	REL606	2900708	2900709
UN	1129	.	REL606	2900712	2900714
UN	1130	.	REL606	2900922	2900924
UN	1131	.	REL606	2906302	2906310
UN	1132	.	REL606	2918565	2918572
UN	1133	.	REL606	2920271	2920275
UN	1134	.	REL606	2933226	2933226
UN	1135	.	REL606	2933230	2933231
UN	1136	.	REL606	2935487	2935488
UN	1137	.	REL606	2942553	2942553
UN	1138	.	REL606	2945515	2945520
UN	1139	.	REL606	2952366	2952377
UN	1140	.	REL606	2956013	2956013
UN	1141	.	REL606	2957190	2957190
UN	1142	.	REL606	2957192	2957192
UN	1143	.	REL606	2969996	2970000
UN	1144	.	REL606	2972145	2972147
UN	1145	.	REL606	2973706	2973707
UN	1146	.	REL606	2974073	2974084
UN	1147	.	REL606	2974485	2974497
UN	1148	.	REL606	2990097	2990113
UN	1149	.	REL606	2990963	2990968
UN	1150	.	REL606	2992238	2992240
UN	1151	.	REL606	2992244	2992244
UN	1152	.	REL606	2992247	2992247
UN	1153	.	REL606	2996107	2996122
UN	1154	.	REL606	2998752	2998757
UN	1155	.	REL606	3001862	3001869
UN	1156	.	REL606	3008651	3008654
UN	1157	.	REL606	3010377	3010377
UN	1158	.	REL606	3015774	3035122
UN	1159	.	REL606	3038035	3038044
UN	1160	.	REL606	3042913	3042921
UN	1161	.	REL606	3045671	3045671
UN	1162	.	REL606	3046069	3046069
UN	1163	.	REL606	3053113	3053116
UN	1164	.	REL606	3061161	3061162
UN	1165	.	REL606	3063209	3063209
UN	1166	.	REL606	3064549	3064550
UN	1167	.	REL606	3070444	3070444
UN	1168	.	REL606	3070451	3070451
UN	1169	.	REL606	3074441	3074445
UN	1170	.	REL606	3075089	3075089
UN	1171	.	REL606	3077736	3077737
UN	1172	.	REL606	3077740	3077740
UN	1173	.	REL606	3077742	3077743
UN	1174	.	REL606	3080187	3080223
UN	1175	.	REL606	3080269	3080298
UN	1176	.	REL606	3090015	3090017
UN	1177	.	REL606	3091827	3091836
UN	1178	.	REL606	3093624	3093629
UN	1179	.	REL606	3109530	3109530
UN	1180	.	REL606	3109533	3109533
UN	1181	.	REL606	3110527	3110527
UN	1182	.	REL606	3116658	3116668
UN	1183	.	REL606	3117243	3117250
UN	1184	.	REL606	3127683	3127688
UN	1185	.	REL606	3128341	3128341
UN	1186	.	REL606	3131645	3131652
UN	1187	.	REL606	3136026	3136040
UN	1188	.	REL606	3142674	3142679
UN	1189	.	REL606	3148534	3148534
UN	1190	.	REL606	3154911	3154914
UN	1191	.	REL606	3158099	3158106
UN	1192	.	REL606	3161589	3161597
UN	1193	.	REL606	3164188	3164191
UN	1194	.	REL606	3172554	3172555
UN	1195	.	REL606	3172645	3172653
UN	1196	.	REL606	3173889	3173889
UN	1197	.	REL606	3185360	3185365
UN	1198	.	REL606	3186531	3186532
UN	1199	.	REL606	3192260	3192263
UN	1200	.	REL606	3193567	3193584
UN	1201	.	REL606	3193733	3193733
UN	1202	.	REL606	3193736	3193736
UN	1203	.	REL606	3194418	3194425
UN	1204	.	REL606	3194510	3194515
UN	1205	.	REL606	3194518	3194519
UN	1206	.	REL606	3196683	3196683
UN	1207	.	REL606	3196685	3196692
UN	1208	.	REL606	3203381	3203383
UN	1209	.	REL606	3204176	3204176
UN	1210	.	REL606	3204180	3204180
UN	1211	.	REL606	3204182	3204182
UN	1212	.	REL606	3210526	3210528
UN	1213	.	REL606	3212029	3212030
UN	1214	.	REL606	3212549	3212549
UN	1215	.	REL606	3226605	3226606
UN	1216	.	REL606	3239081	3239082
UN	1217	.	REL606	3241242	3241242
UN	1218	.	REL606	3241244	3241245
UN	1219	.	REL606	3241247	3241247
UN	1220	.	REL606	3245437	3245438
UN	1221	.	REL606	3253206	3253206
UN	1222	.	REL606	3257717	3257732
UN	1223	.	REL606	3260565	3260570
UN	1224	.	REL606	3261728	3261738
UN	1225	.	REL606	3266101	3266101
UN	1226	.	REL606	3267951	3267954
UN	1227	.	REL606	3283428	3283429
UN	1228	.	REL606	3285201	3285201
UN	1229	.	REL606	3288953	3288961
UN	1230	.	REL606	3292666	3292667
UN	1231	.	REL606	3297607	3297614
UN	1232	.	REL606	3297875	3297875
UN	1233	.	REL606	3297878	3297881
UN	1234	.	REL606	3301393	3301393
UN	1235	.	REL606	3301395	3301395
UN	1236	.	REL606	3304150	3304152
UN	1237	.	REL606	3307590	3307593
UN	1238	.	REL606	3310387	3310398
UN	1239	.	REL606	3310401	3310410
UN	1240	.	REL606	3311036	3311038
UN	1241	.	REL606	3314094	3314097
UN	1242	.	REL606	3315720	3315739
UN	1243	.	REL606	3316499	3316500
UN	1244	.	REL606	3317223	3317227
UN	1245	.	REL606	3317974	3317981
UN	1246	.	REL606	3318115	3318115
UN	1247	.	REL606	3320034	3320037
UN	1248	.	REL606	3320073	3320132
UN	1249	.	REL606	3320181	3320223
UN	1250	.	REL606	3324551	3324552
UN	1251	.	REL606	3324835	3324840
UN	1252	.	REL606	3333567	3333570
UN	1253	.	REL606	3337125	3337125
UN	1254	.	REL606	3337127	3337127
UN	1255	.	REL606	3351163	3351175
UN	1256	.	REL606	3351345	3351391
UN	1257	.	REL606	3351590	3351590
UN	1258	.	REL606	3351593	3351635
UN	1259	.	REL606	3351685	3351690
UN	1260	.	REL606	3351898	3351912
UN	1261	.	REL606	3352026	3352045
UN	1262	.	REL606	3352075	3352075
UN	1263	.	REL606	3352179	3352179
UN	1264	.	REL606	3352436	3352436
UN	1265	.	REL606	3352467	3352467
UN	1266	.	REL606	3352469	3352471
UN	1267	.	REL606	3352541	3352542
UN	1268	.	REL606	3354413	3354428
UN	1269	.	REL606	3354834	3354834
UN	1270	.	REL606	3355003	3355008
UN	1271	.	REL606	3355014	3355015
UN	1272	.	REL606	3355052	3355056
UN	1273	.	REL606	3355367	3355367
UN	1274	.	REL606	3355369	3355370
UN	1275	.	REL606	3355376	3355394
UN	1276	.	REL606	3355840	3355846
UN	1277	.	REL606	3355950	3355985
UN	1278	.	REL606	3356074	3356079
UN	1279	.	REL606	3356311	3356318
UN	1280	.	REL606	3356333	3356344
UN	1281	.	REL606	3356450	3356450
UN	1282	.	REL606	3356459	3356464
UN	1283	.	REL606	3356466	3356467
UN	1284	.	REL606	3356469	3356469
UN	1285	.	REL606	3356494	3356501
UN	1286	.	REL606	3356731	3356747
UN	1287	.	REL606	3356804	3356818
UN	1288	.	REL606	3356824	3356843
UN	1289	.	REL606	3364756	3364756
UN	1290	.	REL606	3364758	3364758
UN	1291	.	REL606	3368191	3368199
UN	1292	.	REL606	3369527	3369532
UN	1293	.	REL606	3379947	3379952
UN	1294	.	REL606	3379954	3379955
UN	1295	.	REL606	3380977	3380979
UN	1296	.	REL606	3391425	3391434
UN	1297	.	REL606	3394108	3394108
UN	1298	.	REL606	3398152	3398157
UN	1299	.	REL606	3398488	3398490
UN	1300	.	REL606	3398664	3398665
UN	1301	.	REL606	3398917	3398917
UN	1302	.	REL606	3400481	3400484
UN	1303	.	REL606	3401856	3401866
UN	1304	.	REL606	3404494	3404496
UN	1305	.	REL606	3405923	3405931
UN	1306	.	REL606	3406333	3406340
UN	1307	.	REL606	3409798	3409803
UN	1308	.	REL606	3419607	3419612
UN	1309	.	REL606	3426630	3426632
UN	1310	.	REL606	3429362	3429362
UN	1311	.	REL606	3430728	3430734
UN	1312	.	REL606	3432750	3432754
UN	1313	.	REL606	3440987	3440992
UN	1314	.	REL606	3447472	3447473
UN	1315	.	REL606	3448180	3448187
UN	1316	.	REL606	3456070	3456081
UN	1317	.	REL606	3456915	3456915
UN	1318	.	REL606	3456926	3456929
UN	1319	.	REL606	3459929	3459931
UN	1320	.	REL606	3463181	3463185
UN	1321	.	REL606	3465001	3465010
UN	1322	.	REL606	3474629	3474634
UN	1323	.	REL606	3474855	3474855
UN	1324	.	REL606	3488466	3488470
UN	1325	.	REL606	3489749	3489758
UN	1326	.	REL606	3497656	3497659
UN	1327	.	REL606	3510159	3510160
UN	1328	.	REL606	3511244	3511244
UN	1329	.	REL606	3517089	3517094
UN	1330	.	REL606	3518574	3518576
UN	1331	.	REL606	3518578	3518581
UN	1332	.	REL606	3518583	3518591
UN	1333	.	REL606	3518597	3518601
UN	1334	.	REL606	3522991	3522992
UN	1335	.	REL606	3525462	3525473
UN	1336	.	REL606	3525616	3525622
UN	1337	.	REL606	3525625	3525625
UN	1338	.	REL606	3526226	3526230
UN	1339	.	REL606	3538295	3538304
UN	1340	.	REL606	3539408	3539420
UN	1341	.	REL606	3542497	3542500
UN	1342	.	REL606	3545006	3545008
UN	1343	.	REL606	3545087	3545092
UN	1344	.	REL606	3546876	3546880
UN	1345	.	REL606	3549178	3549180
UN	1346	.	REL606	3549930	3550088
UN	1347	.	REL606	3550225	3550248
UN	1348	.	REL606	3550335	3550431
UN	1349	.	REL606	3550445	3550445
UN	1350	.	REL606	3550448	3550448
UN	1351	.	REL606	3550452	3550452
UN	1352	.	REL606	3550454	3550454
UN	1353	.	REL606	3550460	3550487
UN	1354	.	REL606	3550616	3550617
UN	1355	.	REL606	3550667	3550672
UN	1356	.	REL606	3550748	3550756
UN	1357	.	REL606	3550844	3550850
UN	1358	.	REL606	3551037	3551049
UN	1359	.	REL606	3551281	3551297
UN	1360	.	REL606	3551330	3551352
UN	1361	.	REL606	3551389	3551407
UN	1362	.	REL606	3551590	3551594
UN	1363	.	REL606	3551666	3551671
UN	1364	.	REL606	3552257	3552319
UN	1365	.	REL606	3552372	3552395
UN	1366	.	REL606	3552429	3552449
UN	1367	.	REL606	3552453	3552453
UN	1368	.	REL606	3552483	3552502
UN	1369	.	REL606	3552534	3552543
UN	1370	.	REL606	3552615	3552623
UN	1371	.	REL606	3552684	3552714
UN	1372	.	REL606	3552749	3552791
UN	1373	.	REL606	3552821	3552843
UN	1374	.	REL606	3552845	3552845
UN	1375	.	REL606	3552885	3552907
UN	1376	.	REL606	3552940	3552951
UN	1377	.	REL606	3552983	3553025
UN	1378	.	REL606	3553087	3553158
UN	1379	.	REL606	3553222	3553257
UN	1380	.	REL606	3553458	3553458
UN	1381	.	REL606	3553463	3553463
UN	1382	.	REL606	3553504	3553524
UN	1383	.	REL606	3554643	3554660
UN	1384	.	REL606	3555138	3555140
UN	1385	.	REL606	3555317	3555323
UN	1386	.	REL606	3555513	3555521
UN	1387	.	REL606	3555765	3555768
UN	1388	.	REL606	3555911	3555917
UN	1389	.	REL606	3559762	3559766
UN	1390	.	REL606	3559977	3559983
UN	1391	.	REL606	3560268	3560270
UN	1392	.	REL606	3560690	3560692
UN	1393	.	REL606	3563577	3563592
UN	1394	.	REL606	3571973	3571973
UN	1395	.	REL606	3575166	3575166
UN	1396	.	REL606	3575974	3575974
UN	1397	.	REL606	3577450	3577461
UN	1398	.	REL606	3577471	3577471
UN	1399	.	REL606	3578807	3578817
UN	1400	.	REL606	3579038	3579043
UN	1401	.	REL606	3581363	3581363
UN	1402	.	REL606	3584278	3584279
UN	1403	.	REL606	3586725	3586725
UN	1404	.	REL606	3589516	3589518
UN	1405	.	REL606	3592534	3592536
UN	1406	.	REL606	3594276	3594281
UN	1407	.	REL606	3595987	3595992
UN	1408	.	REL606	3596039	3596042
UN	1409	.	REL606	3596044	3596044
UN	1410	.	REL606	3596046	3596046
UN	1411	.	REL606	3596049	3596050
UN	1412	.	REL606	3596237	3596238
UN	1413	.	REL606	3596454	3596470
UN	1414	.	REL606	3605616	3605616
UN	1415	.	REL606	3607815	3607822
UN	1416	.	REL606	3607924	3607925
UN	1417	.	REL606	3608011	3608012
UN	1418	.	REL606	3611421	3611421
UN	1419	.	REL606	3619952	3619953
UN	1420	.	REL606	3619955	3619955
UN	1421	.	REL606	3620197	3620201
UN	1422	.	REL606	3629406	3629408
UN	1423	.	REL606	3629440	3629458
UN	1424	.	REL606	3629499	3629538
UN	1425	.	REL606	3629676	3629757
UN	1426	.	REL606	3629824	3629887
UN	1427	.	REL606	3630034	3630041
UN	1428	.	REL606	3630159	3630174
UN	1429	.	REL606	3630320	3630371
UN	1430	.	REL606	3630406	3630421
UN	1431	.	REL606	3630482	3630529
UN	1432	.	REL606	3630642	3630723
UN	1433	.	REL606	3631054	3631058
UN	1434	.	REL606	3637433	3637433
UN	1435	.	REL606	3638825	3638825
UN	1436	.	REL606	3638827	3638899
UN	1437	.	REL606	3639014	3639067
UN	1438	.	REL606	3645600	3645602
UN	1439	.	REL606	3655200	3655200
UN	1440	.	REL606	3657613	3657614
UN	1441	.	REL606	3658172	3658178
UN	1442	.	REL606	3663902	3663910
UN	1443	.	REL606	3669123	3669128
UN	1444	.	REL606	3671079	3671079
UN	1445	.	REL606	3678254	3678254
UN	1446	.	REL606	3681255	3681255
UN	1447	.	REL606	3683677	3683682
UN	1448	.	REL606	3688148	3688154
UN	1449	.	REL606	3695323	3695330
UN	1450	.	REL606	3695408	3695408
UN	1451	.	REL606	3695410	3695410
UN	1452	.	REL606	3697231	3697245
UN	1453	.	REL606	3697422	3697428
UN	1454	.	REL606	3697485	3697485
UN	1455	.	REL606	3697490	3697490
UN	1456	.	REL606	3697501	3697501
UN	1457	.	REL606	3697601	3697605
UN	1458	.	REL606	3697656	3697669
UN	1459	.	REL606	3697792	3697792
UN	1460	.	REL606	3697798	3697807
UN	1461	.	REL606	3697883	3697884
UN	1462	.	REL606	3698034	3698046
UN	1463	.	REL606	3698197	3698207
UN	1464	.	REL606	3698576	3698581
UN	1465	.	REL606	3699311	3699316
UN	1466	.	REL606	3699480	3699493
UN	1467	.	REL606	3699743	3699744
UN	1468	.	REL606	3700140	3700147
UN	1469	.	REL606	3700190	3700196
UN	1470	.	REL606	3700226	3700235
UN	1471	.	REL606	3700281	3700379
UN	1472	.	REL606	3700463	3700515
UN	1473	.	REL606	3700564	3700570
UN	1474	.	REL606	3701943	3701988
UN	1475	.	REL606	3702052	3702132
UN	1476	.	REL606	3702228	3702239
UN	1477	.	REL606	3702314	3702314
UN	1478	.	REL606	3703026	3703028
UN	1479	.	REL606	3706098	3706105
UN	1480	.	REL606	3714884	3714894
UN	1481	.	REL606	3715118	3715133
UN	1482	.	REL606	3719707	3719728
UN	1483	.	REL606	3720128	3720135
UN	1484	.	REL606	3720928	3720931
UN	1485	.	REL606	3720936	3720939
UN	1486	.	REL606	3725917	3725920
UN	1487	.	REL606	3741546	3741550
UN	1488	.	REL606	3741955	3742149
UN	1489	.	REL606	3749564	3749564
UN	1490	.	REL606	3752118	3752122
UN	1491	.	REL606	3764952	3764955
UN	1492	.	REL606	3788834	3788834
UN	1493	.	REL606	3788838	3788839
UN	1494	.	REL606	3792917	3792922
UN	1495	.	REL606	3794639	3794642
UN	1496	.	REL606	3797261	3797261
UN	1497	.	REL606	3799600	3799602
UN	1498	.	REL606	3804278	3804279
UN	1499	.	REL606	3806050	3806051
UN	1500	.	REL606	3809018	3809018
UN	1501	.	REL606	3809020	3809021
UN	1502	.	REL606	3809025	3809025
UN	1503	.	REL606	3809504	3809514
UN	1504	.	REL606	3812481	3812486
UN	1505	.	REL606	3813996	3813998
UN	1506	.	REL606	3819941	3819944
UN	1507	.	REL606	3820910	3820911
UN	1508	.	REL606	3827984	3827986
UN	1509	.	REL606	3827988	3827988
UN	1510	.	REL606	3834450	3834453
UN	1511	.	REL606	3837927	3837929
UN	1512	.	REL606	3840687	3840687
UN	1513	.	REL606	3842250	3842258
UN	1514	.	REL606	3842580	3842586
UN	1515	.	REL606	3863903	3863907
UN	1516	.	REL606	3868893	3868907
UN	1517	.	REL606	3869188	3869197
UN	1518	.	REL606	3876963	3876971
UN	1519	.	REL606	3889897	3889902
UN	1520	.	REL606	3890952	3890954
UN	1521	.	REL606	3893793	3893800
UN	1522	.	REL606	3893830	3893836
UN	1523	.	REL606	3893931	3893935
UN	1524	.	REL606	3894095	3894095
UN	1525	.	REL606	3894097	3894100
UN	1526	.	REL606	3894138	3894149
UN	1527	.	REL606	3894197	3894200
UN	1528	.	REL606	3894233	3894238
UN	1529	.	REL606	3894302	3894312
UN	1530	.	REL606	3894451	3894451
UN	1531	.	REL606	3894858	3894858
UN	1532	.	REL606	3894861	3894865
UN	1533	.	REL606	3894991	3901456
UN	1534	.	REL606	3903318	3903433
UN	1535	.	REL606	3903502	3903502
UN	1536	.	REL606	3903509	3903509
UN	1537	.	REL606	3903513	3903513
UN	1538	.	REL606	3903515	3903773
UN	1539	.	REL606	3903778	3903778
UN	1540	.	REL606	3903805	3903918
UN	1541	.	REL606	3904005	3904053
UN	1542	.	REL606	3904073	3904073
UN	1543	.	REL606	3904075	3904075
UN	1544	.	REL606	3904080	3904080
UN	1545	.	REL606	3904082	3904082
UN	1546	.	REL606	3904085	3904087
UN	1547	.	REL606	3904122	3904124
UN	1548	.	REL606	3904154	3904166
UN	1549	.	REL606	3904210	3904210
UN	1550	.	REL606	3904213	3904224
UN	1551	.	REL606	3904279	3904279
UN	1552	.	REL606	3904285	3904286
UN	1553	.	REL606	3904290	3904292
UN	1554	.	REL606	3904297	3904320
UN	1555	.	REL606	3904324	3904337
UN	1556	.	REL606	3904567	3904571
UN	1557	.	REL606	3904891	3904891
UN	1558	.	REL606	3904895	3904898
UN	1559	.	REL606	3904963	3904970
UN	1560	.	REL606	3905514	3905520
UN	1561	.	REL606	3905522	3905522
UN	1562	.	REL606	3905584	3905593
UN	1563	.	REL606	3905759	3905760
UN	1564	.	REL606	3905763	3905763
UN	1565	.	REL606	3905765	3905766
UN	1566	.	REL606	3906250	3906254
UN	1567	.	REL606	3906785	3906789
UN	1568	.	REL606	3907309	3907311
UN	1569	.	REL606	3907314	3907316
UN	1570	.	REL606	3907573	3907595
UN	1571	.	REL606	3908000	3908000
UN	1572	.	REL606	3908116	3908123
UN	1573	.	REL606	3908163	3908168
UN	1574	.	REL606	3908267	3908281
UN	1575	.	REL606	3908437	3908442
UN	1576	.	REL606	3908489	3908494
UN	1577	.	REL606	3910976	3910983
UN	1578	.	REL606	3910985	3910985
UN	1579	.	REL606	3920272	3920275
UN	1580	.	REL606	3926906	3926907
UN	1581	.	REL606	3928784	3928785
UN	1582	.	REL606	3928787	3928789
UN	1583	.	REL606	3931640	3931642
UN	1584	.	REL606	3934355	3934357
UN	1585	.	REL606	3937983	3938000
UN	1586	.	REL606	3940239	3940240
UN	1587	.	REL606	3941167	3941168
UN	1588	.	REL606	3941777	3941786
UN	1589	.	REL606	3942764	3942764
UN	1590	.	REL606	3946126	3946139
UN	1591	.	REL606	3947537	3947542
UN	1592	.	REL606	3949676	3949676
UN	1593	.	REL606	3949731	3949769
UN	1594	.	REL606	3950052	3950052
UN	1595	.	REL606	3951427	3951427
UN	1596	.	REL606	3951675	3951696
UN	1597	.	REL606	3954410	3954424
UN	1598	.	REL606	3956239	3956242
UN	1599	.	REL606	3957253	3957259
UN	1600	.	REL606	3957667	3957668
UN	1601	.	REL606	3958485	3958493
UN	1602	.	REL606	3958648	3958658
UN	1603	.	REL606	3965057	3965085
UN	1604	.	REL606	3965605	3965607
UN	1605	.	REL606	3966324	3966324
UN	1606	.	REL606	3966326	3966326
UN	1607	.	REL606	3966751	3966751
UN	1608	.	REL606	3966755	3966755
UN	1609	.	REL606	3966757	3966782
UN	1610	.	REL606	3977496	3977498
UN	1611	.	REL606	3979130	3979134
UN	1612	.	REL606	3985691	3985699
UN	1613	.	REL606	3992403	3992403
UN	1614	.	REL606	3992411	3992412
UN	1615	.	REL606	4003069	4003070
UN	1616	.	REL606	4004872	4004877
UN	1617	.	REL606	4007536	4007536
UN	1618	.	REL606	4009136	4009139
UN	1619	.	REL606	4010834	4010837
UN	1620	.	REL606	4010933	4010934
UN	1621	.	REL606	4013648	4013659
UN	1622	.	REL606	4013755	4013760
UN	1623	.	REL606	4014017	4014017
UN	1624	.	REL606	4014019	4014036
UN	1625	.	REL606	4014098	4014112
UN	1626	.	REL606	4014255	4014256
UN	1627	.	REL606	4014258	4014259
UN	1628	.	REL606	4014261	4014262
UN	1629	.	REL606	4014265	4014271
UN	1630	.	REL606	4014422	4014434
UN	1631	.	REL606	4014439	4014444
UN	1632	.	REL606	4014446	4014446
UN	1633	.	REL606	4014492	4014494
UN	1634	.	REL606	4014696	4014703
UN	1635	.	REL606	4014803	4014808
UN	1636	.	REL606	4014840	4014840
UN	1637	.	REL606	4014916	4014927
UN	1638	.	REL606	4015009	4015019
UN	1639	.	REL606	4015102	4015110
UN	1640	.	REL606	4015184	4015195
UN	1641	.	REL606	4015435	4015450
UN	1642	.	REL606	4015563	4015568
UN	1643	.	REL606	4015638	4015645
UN	1644	.	REL606	4015650	4015652
UN	1645	.	REL606	4015964	4015984
UN	1646	.	REL606	4016016	4016022
UN	1647	.	REL606	4016146	4016168
UN	1648	.	REL606	4017178	4017179
UN	1649	.	REL606	4017746	4017773
UN	1650	.	REL606	4018131	4018131
UN	1651	.	REL606	4018525	4018558
UN	1652	.	REL606	4018671	4018680
UN	1653	.	REL606	4033085	4033089
UN	1654	.	REL606	4033091	4033091
UN	1655	.	REL606	4046100	4046105
UN	1656	.	REL606	4054645	4054645
UN	1657	.	REL606	4056342	4056342
UN	1658	.	REL606	4058549	4058552
UN	1659	.	REL606	4060834	4060851
UN	1660	.	REL606	4063061	4063062
UN	1661	.	REL606	4063064	4063065
UN	1662	.	REL606	4066391	4066391
UN	1663	.	REL606	4066393	4066393
UN	1664	.	REL606	4073508	4073512
UN	1665	.	REL606	4080298	4080301
UN	1666	.	REL606	4091676	4091680
UN	1667	.	REL606	4092235	4092235
UN	1668	.	REL606	4094070	4094073
UN	1669	.	REL606	4094655	4094655
UN	1670	.	REL606	4095117	4095126
UN	1671	.	REL606	4095128	4095137
UN	1672	.	REL606	4102215	4102219
UN	1673	.	REL606	4102222	4102223
UN	1674	.	REL606	4103723	4103724
UN	1675	.	REL606	4108226	4108233
UN	1676	.	REL606	4110232	4110238
UN	1677	.	REL606	4110679	4110683
UN	1678	.	REL606	4110722	4110722
UN	1679	.	REL606	4128916	4128919
UN	1680	.	REL606	4132426	4132426
UN	1681	.	REL606	4132429	4132430
UN	1682	.	REL606	4135209	4135209
UN	1683	.	REL606	4137350	4137359
UN	1684	.	REL606	4141658	4141670
UN	1685	.	REL606	4146249	4146249
UN	1686	.	REL606	4146371	4146377
UN	1687	.	REL606	4146424	4146456
UN	1688	.	REL606	4146768	4146769
UN	1689	.	REL606	4146771	4146772
UN	1690	.	REL606	4146775	4146775
UN	1691	.	REL606	4146779	4146782
UN	1692	.	REL606	4146836	4146881
UN	1693	.	REL606	4147059	4147071
UN	1694	.	REL606	4147191	4147191
UN	1695	.	REL606	4147229	4147237
UN	1696	.	REL606	4147276	4147312
UN	1697	.	REL606	4147474	4147476
UN	1698	.	REL606	4147679	4147699
UN	1699	.	REL606	4147706	4147706
UN	1700	.	REL606	4147806	4147811
UN	1701	.	REL606	4148037	4148047
UN	1702	.	REL606	4148050	4148051
UN	1703	.	REL606	4149850	4149853
UN	1704	.	REL606	4149886	4149913
UN	1705	.	REL606	4150225	4150227
UN	1706	.	REL606	4150595	4150595
UN	1707	.	REL606	4150602	4150603
UN	1708	.	REL606	4150759	4150763
UN	1709	.	REL606	4151405	4151410
UN	1710	.	REL606	4160259	4160261
UN	1711	.	REL606	4160263	4160270
UN	1712	.	REL606	4161500	4161502
UN	1713	.	REL606	4161702	4161706
UN	1714	.	REL606	4161708	4161715
UN	1715	.	REL606	4161956	4161958
UN	1716	.	REL606	4163015	4163015
UN	1717	.	REL606	4163018	4163018
UN	1718	.	REL606	4173811	4173812
UN	1719	.	REL606	4187503	4187508
UN	1720	.	REL606	4187510	4187510
UN	1721	.	REL606	4187512	4187513
UN	1722	.	REL606	4187515	4187516
UN	1723	.	REL606	4187556	4187591
UN	1724	.	REL606	4187638	4187645
UN	1725	.	REL606	4187819	4187823
UN	1726	.	REL606	4187952	4187952
UN	1727	.	REL606	4187998	4188002
UN	1728	.	REL606	4188006	4188007
UN	1729	.	REL606	4188010	4188020
UN	1730	.	REL606	4188071	4188101
UN	1731	.	REL606	4188192	4188197
UN	1732	.	REL606	4188294	4188295
UN	1733	.	REL606	4188365	4188365
UN	1734	.	REL606	4188374	4188374
UN	1735	.	REL606	4188377	4188386
UN	1736	.	REL606	4188420	4188420
UN	1737	.	REL606	4188422	4188422
UN	1738	.	REL606	4188424	4188425
UN	1739	.	REL606	4188427	4188427
UN	1740	.	REL606	4188430	4188458
UN	1741	.	REL606	4188521	4188527
UN	1742	.	REL606	4188531	4188531
UN	1743	.	REL606	4188565	4188574
UN	1744	.	REL606	4188637	4188654
UN	1745	.	REL606	4188688	4188690
UN	1746	.	REL606	4188892	4188893
UN	1747	.	REL606	4188895	4188896
UN	1748	.	REL606	4188900	4188900
UN	1749	.	REL606	4188904	4188904
UN	1750	.	REL606	4188908	4188912
UN	1751	.	REL606	4188921	4188923
UN	1752	.	REL606	4189082	4189089
UN	1753	.	REL606	4189239	4189239
UN	1754	.	REL606	4189327	4189362
UN	1755	.	REL606	4189683	4189717
UN	1756	.	REL606	4189878	4189885
UN	1757	.	REL606	4190904	4190908
UN	1758	.	REL606	4191536	4191541
UN	1759	.	REL606	4191600	4191601
UN	1760	.	REL606	4191603	4191604
UN	1761	.	REL606	4192275	4192292
UN	1762	.	REL606	4192386	4192402
UN	1763	.	REL606	4192437	4192438
UN	1764	.	REL606	4192621	4192629
UN	1765	.	REL606	4192691	4192708
UN	1766	.	REL606	4192972	4192979
UN	1767	.	REL606	4193025	4193026
UN	1768	.	REL606	4193028	4193030
UN	1769	.	REL606	4197902	4197935
UN	1770	.	REL606	4198000	4198042
UN	1771	.	REL606	4205151	4205161
UN	1772	.	REL606	4216108	4216108
UN	1773	.	REL606	4216110	4216111
UN	1774	.	REL606	4216513	4216515
UN	1775	.	REL606	4216583	4216592
UN	1776	.	REL606	4226332	4226332
UN	1777	.	REL606	4231485	4231500
UN	1778	.	REL606	4231972	4231980
UN	1779	.	REL606	4236676	4236681
UN	1780	.	REL606	4253543	4253546
UN	1781	.	REL606	4256760	4256761
UN	1782	.	REL606	4256883	4256885
UN	1783	.	REL606	4258405	4258415
UN	1784	.	REL606	4262616	4262616
UN	1785	.	REL606	4269615	4269622
UN	1786	.	REL606	4272602	4272603
UN	1787	.	REL606	4274781	4274980
UN	1788	.	REL606	4275016	4275199
UN	1789	.	REL606	4275831	4275835
UN	1790	.	REL606	4276238	4276240
UN	1791	.	REL606	4281586	4281592
UN	1792	.	REL606	4287752	4287753
UN	1793	.	REL606	4288252	4288253
UN	1794	.	REL606	4293144	4293145
UN	1795	.	REL606	4294573	4294576
UN	1796	.	REL606	4296278	4296291
UN	1797	.	REL606	4298937	4298942
UN	1798	.	REL606	4301799	4301799
UN	1799	.	REL606	4304678	4304738
UN	1800	.	REL606	4304779	4304824
UN	1801	.	REL606	4304886	4304937
UN	1802	.	REL606	4304979	4305041
UN	1803	.	REL606	4308525	4308525
UN	1804	.	REL606	4311777	4311785
UN	1805	.	REL606	4313550	4313550
UN	1806	.	REL606	4314765	4314771
UN	1807	.	REL606	4319772	4319772
UN	1808	.	REL606	4321564	4321564
UN	1809	.	REL606	4333499	4333504
UN	1810	.	REL606	4334727	4334727
UN	1811	.	REL606	4336356	4336356
UN	1812	.	REL606	4340305	4340306
UN	1813	.	REL606	4340405	4340405
UN	1814	.	REL606	4341246	4341246
UN	1815	.	REL606	4341309	4341313
UN	1816	.	REL606	4341730	4341732
UN	1817	.	REL606	4344178	4344179
UN	1818	.	REL606	4344181	4344181
UN	1819	.	REL606	4345473	4345478
UN	1820	.	REL606	4346152	4346156
UN	1821	.	REL606	4346158	4346158
UN	1822	.	REL606	4355404	4355434
UN	1823	.	REL606	4357681	4357689
UN	1824	.	REL606	4370544	4370569
UN	1825	.	REL606	4370634	4370681
UN	1826	.	REL606	4370743	4370788
UN	1827	.	REL606	4373056	4373068
UN	1828	.	REL606	4375992	4375994
UN	1829	.	REL606	4375996	4375996
UN	1830	.	REL606	4375998	4376002
UN	1831	.	REL606	4378957	4378961
UN	1832	.	REL606	4385479	4385479
UN	1833	.	REL606	4387959	4387959
UN	1834	.	REL606	4397337	4397340
UN	1835	.	REL606	4400797	4400797
UN	1836	.	REL606	4414287	4414287
UN	1837	.	REL606	4415713	4415715
UN	1838	.	REL606	4419299	4419301
UN	1839	.	REL606	4429786	4429787
UN	1840	.	REL606	4432008	4432012
UN	1841	.	REL606	4432227	4432232
UN	1842	.	REL606	4439326	4439326
UN	1843	.	REL606	4439328	4439328
UN	1844	.	REL606	4439723	4439723
UN	1845	.	REL606	4440220	4440226
UN	1846	.	REL606	4442647	4442648
UN	1847	.	REL606	4444557	4444560
UN	1848	.	REL606	4445710	4445717
UN	1849	.	REL606	4449594	4449602
UN	1850	.	REL606	4453735	4453738
UN	1851	.	REL606	4458957	4458957
UN	1852	.	REL606	4458959	4458959
UN	1853	.	REL606	4463980	4463982
UN	1854	.	REL606	4467621	4467622
UN	1855	.	REL606	4470088	4470098
UN	1856	.	REL606	4471774	4471776
UN	1857	.	REL606	4472886	4472887
UN	1858	.	REL606	4475598	4475601
UN	1859	.	REL606	4475678	4475680
UN	1860	.	REL606	4475726	4475726
UN	1861	.	REL606	4479960	4479967
UN	1862	.	REL606	4489789	4489798
UN	1863	.	REL606	4494154	4494163
UN	1864	.	REL606	4494169	4494169
UN	1865	.	REL606	4494227	4494232
UN	1866	.	REL606	4496789	4496794
UN	1867	.	REL606	4499588	4499591
UN	1868	.	REL606	4502287	4502289
UN	1869	.	REL606	4502937	4502939
UN	1870	.	REL606	4505173	4505173
UN	1871	.	REL606	4507500	4507501
UN	1872	.	REL606	4507890	4507894
UN	1873	.	REL606	4510382	4510384
UN	1874	.	REL606	4516169	4516170
UN	1875	.	REL606	4521780	4521780
UN	1876	.	REL606	4521783	4521783
UN	1877	.	REL606	4521786	4521788
UN	1878	.	REL606	4521790	4521790
UN	1879	.	REL606	4521792	4521794
UN	1880	.	REL606	4521800	4521801
UN	1881	.	REL606	4521814	4521814
UN	1882	.	REL606	4521919	4521925
UN	1883	.	REL606	4522063	4522080
UN	1884	.	REL606	4522082	4522084
UN	1885	.	REL606	4522277	4561282
UN	1886	.	REL606	4561801	4561801
UN	1887	.	REL606	4568562	4568569
UN	1888	.	REL606	4573169	4573172
UN	1889	.	REL606	4574051	4574055
UN	1890	.	REL606	4574854	4574857
UN	1891	.	REL606	4575574	4575577
UN	1892	.	REL606	4578947	4578948
UN	1893	.	REL606	4579661	4579662
UN	1894	.	REL606	4585660	4585660
UN	1895	.	REL606	4588671	4588671
UN	1896	.	REL606	4589273	4589274
UN	1897	.	REL606	4595342	4595371
UN	1898	.	REL606	4595580	4595606
UN	1899	.	REL606	4596651	4596651
UN	1900	.	REL606	4599134	4599137
UN	1901	.	REL606	4599589	4599589
UN	1902	.	REL606	4603520	4603526
UN	1903	.	REL606	4603588	4603629
UN	1904	.	REL606	4604402	4604404
UN	1905	.	REL606	4605476	4605476
UN	1906	.	REL606	4607334	4607346
UN	1907	.	REL606	4610901	4610906
UN	1908	.	REL606	4612516	4612518
UN	1909	.	REL606	4612520	4612520
PD	1910	.	REL606	546190	-1	REL606	546191	1	ambiguous_pair_count=634	candidate_covering_count=2023	distinct_pair_count=253	frequency=0.7792	frequency_lower=0.7564	frequency_upper=0.8009	normal_pair_count=217	position_range=21	repeat_name=IS1	repeat_name_evidence=junction	repeat_size_difference=-8	score=30.8	seed_z_score=-35.15	shifted_pair_count=766	side_1_annotate_key=gene	side_1_gene_name=ybcQ	side_1_gene_position=coding (345/384 nt)	side_1_gene_product=predicted antitermination protein	side_1_gene_strand=>	side_1_locus_tag=ECB_00502	side_2_annotate_key=gene	side_2_gene_name=ybcQ	side_2_gene_position=coding (346/384 nt)	side_2_gene_product=predicted antitermination protein	side_2_gene_strand=>	side_2_locus_tag=ECB_00502	size_shift=-760	size_shift_lower=-762	size_shift_upper=-760	total_pair_count=1617
PD	1911	.	REL606	547517	-1	REL606	549789	1	ambiguous_pair_count=0	candidate_covering_count=1096	distinct_pair_count=393	frequency=1.0000	frequency_lower=0.9973	frequency_upper=1.0000	normal_pair_count=0	position_range=1127	score=87.3	seed_z_score=56.91	shifted_pair_count=1096	side_1_annotate_key=gene	side_1_gene_name=IS1	side_1_gene_position=noncoding (183/768 nt)	side_1_gene_product=repeat region	side_1_gene_strand=<	side_2_annotate_key=gene	side_2_gene_name=ybcT	side_2_gene_position=coding (214/462 nt)	side_2_gene_product=predicted murein endopeptidase	side_2_gene_strand=>	side_2_locus_tag=ECB_00508	size_shift=2271	size_shift_lower=2267	size_shift_upper=2276	total_pair_count=1096
PD	1912	.	REL606	643364	-1	REL606	643365	1	ambiguous_pair_count=42	candidate_covering_count=1390	distinct_pair_count=460	frequency=0.9941	frequency_lower=0.9893	frequency_upper=0.9970	normal_pair_count=8	position_range=13	repeat_name=IS150	repeat_name_evidence=junction	repeat_size_difference=64	score=111.7	seed_z_score=-64.07	shifted_pair_count=1340	side_1_annotate_key=gene	side_1_gene_name=ybeF	side_1_gene_position=coding (598/954 nt)	side_1_gene_product=predicted DNA-binding transcriptional regulator	side_1_gene_strand=<	side_1_locus_tag=ECB_00598	side_2_annotate_key=gene	side_2_gene_name=ybeF	side_2_gene_position=coding (597/954 nt)	side_2_gene_product=predicted DNA-binding transcriptional regulator	side_2_gene_strand=<	side_2_locus_tag=ECB_00598	size_shift=-1507	size_shift_lower=-1518	size_shift_upper=-1494	total_pair_count=1390
PD	1913	.	REL606	1047828	-1	REL606	1047829	1	ambiguous_pair_count=50	candidate_covering_count=1104	distinct_pair_count=385	frequency=0.9896	frequency_lower=0.9828	frequency_upper=0.9941	normal_pair_count=11	position_range=17	repeat_name=IS150	repeat_name_evidence=junction	repeat_size_difference=-14	score=87.0	seed_z_score=-56.83	shifted_pair_count=1043	side_1_annotate_key=gene	side_1_gene_name=yccK/yccA	side_1_gene_position=intergenic (-69/+22)	side_1_gene_product=predicted sulfite reductase subunit/inner membrane protein	side_1_gene_strand=</<	side_1_locus_tag=ECB_00973/ECB_00974	side_2_annotate_key=gene	side_2_gene_name=yccK/yccA	side_2_gene_position=intergenic (-70/+21)	side_2_gene_product=predicted sulfite reductase subunit/inner membrane protein	side_2_gene_strand=</<	side_2_locus_tag=ECB_00973/ECB_00974	size_shift=-1429	size_shift_lower=-1432	size_shift_upper=-1420	total_pair_count=1104
PD	1914	.	REL606	1270361	-1	REL606	1270362	1	ambiguous_pair_count=209	candidate_covering_count=624	distinct_pair_count=144	frequency=0.9159	frequency_lower=0.8900	frequency_upper=0.9372	normal_pair_count=35	position_range=612	repeat_name=IS1	repeat_name_evidence=size	repeat_size_candidates=IS1	repeat_size_difference=31	score=33.5	seed_z_score=-36.52	shifted_pair_count=381	side_1_annotate_key=gene	side_1_gene_name=ldrB/ldrC	side_1_gene_position=intergenic (-403/+25)	side_1_gene_product=toxic polypeptide, small/toxic polypeptide, small	side_1_gene_strand=</<	side_1_locus_tag=ECB_01192/ECB_01193	side_2_annotate_key=gene	side_2_gene_name=ldrB/ldrC	side_2_gene_position=intergenic (-404/+24)	side_2_gene_product=toxic polypeptide, small/toxic polypeptide, small	side_2_gene_strand=</<	side_2_locus_tag=ECB_01192/ECB_01193	size_shift=-799	size_shift_lower=-845	size_shift_upper=-777	total_pair_count=625
PD	1915	.	REL606	1420708	-1	REL606	1420709	1	ambiguous_pair_count=74	candidate_covering_count=1151	distinct_pair_count=345	frequency=0.9759	frequency_lower=0.9667	frequency_upper=0.9830	normal_pair_count=26	position_range=4	repeat_name=IS3	repeat_name_evidence=junction	repeat_size_difference=69	score=89.5	seed_z_score=-57.61	shifted_pair_count=1051	side_1_annotate_key=gene	side_1_gene_name=ynaK/ECB_01341	side_1_gene_position=intergenic (+20/-13)	side_1_gene_product=conserved protein/hypothetical protein	side_1_gene_strand=>/>	side_1_locus_tag=ECB_01340/ECB_01341	side_2_annotate_key=gene	side_2_gene_name=ynaK/ECB_01341	side_2_gene_position=intergenic (+21/-12)	side_2_gene_product=conserved protein/hypothetical protein	side_2_gene_strand=>/>	side_2_locus_tag=ECB_01340/ECB_01341	size_shift=-1324	size_shift_lower=-1334	size_shift_upper=-1320	total_pair_count=1151
PD	1916	.	REL606	1428810	-1	REL606	1428811	1	ambiguous_pair_count=612	candidate_covering_count=2175	distinct_pair_count=485	frequency=0.8932	frequency_lower=0.8794	frequency_upper=0.9058	normal_pair_count=167	position_range=23	repeat_name=IS1	repeat_name_evidence=junction	repeat_size_difference=-3	score=118.3	seed_z_score=-65.87	shifted_pair_count=1396	side_1_annotate_key=gene	side_1_gene_name=ompN/ydbK	side_1_gene_position=intergenic (-311/+56)	side_1_gene_product=outer membrane pore protein N, non-specific/fused predicted pyruvate-flavodoxin oxidoreductase: conserved protein/conserved protein/FeS binding protein	side_1_gene_strand=</<	side_1_locus_tag=ECB_01348/ECB_01349	side_2_annotate_key=gene	side_2_gene_name=ompN/ydbK	side_2_gene_position=intergenic (-312/+55)	side_2_gene_product=outer membrane pore protein N, non-specific/fused predicted pyruvate-flavodoxin oxidoreductase: conserved protein/conserved protein/FeS binding protein	side_2_gene_strand=</<	side_2_locus_tag=ECB_01348/ECB_01349	size_shift=-765	size_shift_lower=-766	size_shift_upper=-765	total_pair_count=2175
PD	1917	.	REL606	1451970	-1	REL606	1451971	1	ambiguous_pair_count=54	candidate_covering_count=1229	distinct_pair_count=394	frequency=0.9540	frequency_lower=0.9427	frequency_upper=0.9637	normal_pair_count=54	position_range=13	repeat_name=IS150	repeat_name_evidence=junction	repeat_size_difference=40	score=97.0	seed_z_score=-59.86	shifted_pair_count=1121	side_1_annotate_key=gene	side_1_gene_name=acpD	side_1_gene_position=coding (222/606 nt)	side_1_gene_product=acyl carrier protein phosphodiesterase	side_1_gene_strand=<	side_1_locus_tag=ECB_01367	side_2_annotate_key=gene	side_2_gene_name=acpD	side_2_gene_position=coding (221/606 nt)	side_2_gene_product=acyl carrier protein phosphodiesterase	side_2_gene_strand=<	side_2_locus_tag=ECB_01367	size_shift=-1483	size_shift_lower=-1487	size_shift_upper=-1481	total_pair_count=1229
PD	1918	.	REL606	1462295	-1	REL606	1462296	1	ambiguous_pair_count=34	candidate_covering_count=1235	distinct_pair_count=448	frequency=0.9850	frequency_lower=0.9779	frequency_upper=0.9903	normal_pair_count=18	position_range=96	repeat_name=IS150	repeat_name_evidence=junction	repeat_size_difference=-37	score=97.8	seed_z_score=-60.11	shifted_pair_count=1183	side_1_annotate_key=gene	side_1_gene_name=mokB/trg	side_1_gene_position=intergenic (-57/-284)	side_1_gene_product=regulatory peptide/methyl-accepting chemotaxis protein III, ribose and galactose sensor receptor	side_1_gene_strand=</>	side_1_locus_tag=ECB_01377/ECB_01378	side_2_annotate_key=gene	side_2_gene_name=mokB/trg	side_2_gene_position=intergenic (-58/-283)	side_2_gene_product=regulatory peptide/methyl-accepting chemotaxis protein III, ribose and galactose sensor receptor	side_2_gene_strand=</>	side_2_locus_tag=ECB_01377/ECB_01378	size_shift=-1406	size_shift_lower=-1406	size_shift_upper=-1403	total_pair_count=1235
PD	1919	.	REL606	1506815	-1	REL606	1507093	1	ambiguous_pair_count=2970	candidate_covering_count=645	distinct_pair_count=17	frequency=0.0949	frequency_lower=0.0732	frequency_upper=0.1207	normal_pair_count=410	position_range=1210	reject=PAIR_DISTANCE_FREQUENCY,PAIR_DISTANCE_SCORE	score=1.5	seed_z_score=14.11	shifted_pair_count=43	side_1_annotate_key=gene	side_1_gene_name=ydcC	side_1_gene_position=coding (438/792 nt)	side_1_gene_product=hypothetical protein	side_1_gene_strand=>	side_1_locus_tag=ECB_01418	side_2_annotate_key=gene	side_2_gene_name=ydcC	side_2_gene_position=coding (716/792 nt)	side_2_gene_product=hypothetical protein	side_2_gene_strand=>	side_2_locus_tag=ECB_01418	size_shift=277	size_shift_lower=265	size_shift_upper=303	total_pair_count=3423
PD	1920	.	REL606	1608483	-1	REL606	1616153	1	ambiguous_pair_count=0	candidate_covering_count=261	distinct_pair_count=95	frequency=1.0000	frequency_lower=0.9886	frequency_upper=1.0000	normal_pair_count=0	position_range=1086	score=18.0	seed_z_score=27.96	shifted_pair_count=261	side_1_annotate_key=gene	side_1_gene_name=IS3	side_1_gene_position=noncoding (564/1255 nt)	side_1_gene_product=repeat region	side_1_gene_strand=>	side_2_annotate_key=gene	side_2_gene_name=IS3	side_2_gene_position=noncoding (680/1255 nt)	side_2_gene_product=repeat region	side_2_gene_strand=>	size_shift=7669	size_shift_lower=7614	size_shift_upper=7689	total_pair_count=261
PD	1921	.	REL606	1651969	-1	REL606	1651970	1	ambiguous_pair_count=39	candidate_covering_count=1261	distinct_pair_count=430	frequency=0.9845	frequency_lower=0.9773	frequency_upper=0.9898	normal_pair_count=19	position_range=2	repeat_name=IS150	repeat_name_evidence=junction	repeat_size_difference=-53	score=99.9	seed_z_score=-60.72	shifted_pair_count=1203	side_1_annotate_key=gene	side_1_gene_name=ydgG	side_1_gene_position=coding (862/1035 nt)	side_1_gene_product=predicted inner membrane protein	side_1_gene_strand=>	side_1_locus_tag=ECB_01570	side_2_annotate_key=gene	side_2_gene_name=ydgG	side_2_gene_position=coding (863/1035 nt)	side_2_gene_product=predicted inner membrane protein	side_2_gene_strand=>	side_2_locus_tag=ECB_01570	size_shift=-1390	size_shift_lower=-1397	size_shift_upper=-1383	total_pair_count=1261
PD	1922	.	REL606	1729045	-1	REL606	1729046	1	ambiguous_pair_count=114	candidate_covering_count=1449	distinct_pair_count=421	frequency=0.8577	frequency_lower=0.8410	frequency_upper=0.8732	normal_pair_count=190	position_range=5	repeat_name=IS150	repeat_name_evidence=junction	repeat_size_difference=4	score=54.8	seed_z_score=-45.70	shifted_pair_count=1145	side_1_annotate_key=gene	side_1_gene_name=ydhV	side_1_gene_position=coding (2053/2103 nt)	side_1_gene_product=predicted oxidoreductase	side_1_gene_strand=<	side_1_locus_tag=ECB_01642	side_2_annotate_key=gene	side_2_gene_name=ydhV	side_2_gene_position=coding (2052/2103 nt)	side_2_gene_product=predicted oxidoreductase	side_2_gene_strand=<	side_2_locus_tag=ECB_01642	size_shift=-1447	size_shift_lower=-1447	size_shift_upper=-1447	total_pair_count=1449
PD	1923	.	REL606	1776433	-1	REL606	1776434	1	ambiguous_pair_count=65	candidate_covering_count=1400	distinct_pair_count=498	frequency=0.9888	frequency_lower=0.9828	frequency_upper=0.9931	normal_pair_count=15	position_range=0	repeat_name=IS150	repeat_name_evidence=junction	repeat_size_difference=-23	score=111.5	seed_z_score=-64.01	shifted_pair_count=1320	side_1_annotate_key=gene	side_1_gene_name=pheS/pheM	side_1_gene_position=intergenic (-224/+60)	side_1_gene_product=phenylalanyl-tRNA synthetase alpha subunit/phenylalanyl-tRNA synthetase operon leader peptide	side_1_gene_strand=</<	side_1_locus_tag=ECB_01683/ECB_01684	side_2_annotate_key=gene	side_2_gene_name=pheS/pheM	side_2_gene_position=intergenic (-225/+59)	side_2_gene_product=phenylalanyl-tRNA synthetase alpha subunit/phenylalanyl-tRNA synthetase operon leader peptide	side_2_gene_strand=</<	side_2_locus_tag=ECB_01683/ECB_01684	size_shift=-1420	size_shift_lower=-1429	size_shift_upper=-1417	total_pair_count=1400
PD	1924	.	REL606	2157883	-1	REL606	2157884	1	ambiguous_pair_count=62	candidate_covering_count=913	distinct_pair_count=355	frequency=0.9812	frequency_lower=0.9716	frequency_upper=0.9882	normal_pair_count=16	position_range=4	repeat_name=IS150	repeat_name_evidence=junction	repeat_size_difference=-72	score=70.4	seed_z_score=-51.38	shifted_pair_count=835	side_1_annotate_key=gene	side_1_gene_name=yehM	side_1_gene_position=coding (993/2280 nt)	side_1_gene_product=hypothetical protein	side_1_gene_strand=>	side_1_locus_tag=ECB_02049	side_2_annotate_key=gene	side_2_gene_name=yehM	side_2_gene_position=coding (994/2280 nt)	side_2_gene_product=hypothetical protein	side_2_gene_strand=>	side_2_locus_tag=ECB_02049	size_shift=-1371	size_shift_lower=-1381	size_shift_upper=-1366	total_pair_count=913
PD	1925	.	REL606	2170658	-1	REL606	2170659	1	ambiguous_pair_count=43	candidate_covering_count=1040	distinct_pair_count=404	frequency=0.9769	frequency_lower=0.9675	frequency_upper=0.9842	normal_pair_count=23	position_range=13	repeat_name=IS150	repeat_name_evidence=junction	repeat_size_difference=-22	score=81.8	seed_z_score=-55.18	shifted_pair_count=974	side_1_annotate_key=gene	side_1_gene_name=yehZ	side_1_gene_position=coding (21/918 nt)	side_1_gene_product=predicted transporter subunit: periplasmic-binding component of ABC superfamily	side_1_gene_strand=<	side_1_locus_tag=ECB_02061	side_2_annotate_key=gene	side_2_gene_name=yehZ	side_2_gene_position=coding (20/918 nt)	side_2_gene_product=predicted transporter subunit: periplasmic-binding component of ABC superfamily	side_2_gene_strand=<	side_2_locus_tag=ECB_02061	size_shift=-1421	size_shift_lower=-1447	size_shift_upper=-1413	total_pair_count=1040
PD	1926	.	REL606	2448497	-1	REL606	2448498	1	ambiguous_pair_count=43	candidate_covering_count=1408	distinct_pair_count=445	frequency=0.9780	frequency_lower=0.9703	frequency_upper=0.9841	normal_pair_count=30	position_range=0	repeat_size_candidates=IS150,IS186	score=111.8	seed_z_score=-64.10	shifted_pair_count=1335	side_1_annotate_key=gene	side_1_gene_name=nupC/yfeA	side_1_gene_position=intergenic (+28/+22)	side_1_gene_product=nucleoside (except guanosine) transporter/predicted diguanylate cyclase	side_1_gene_strand=>/<	side_1_locus_tag=ECB_02302/ECB_02303	side_2_annotate_key=gene	side_2_gene_name=nupC/yfeA	side_2_gene_position=intergenic (+29/+21)	side_2_gene_product=nucleoside (except guanosine) transporter/predicted diguanylate cyclase	side_2_gene_strand=>/<	side_2_locus_tag=ECB_02302/ECB_02303	size_shift=-1391	size_shift_lower=-1397	size_shift_upper=-1380	total_pair_count=1408
PD	1927	.	REL606	2655777	-1	REL606	2655778	1	ambiguous_pair_count=54	candidate_covering_count=993	distinct_pair_count=324	frequency=0.9585	frequency_lower=0.9461	frequency_upper=0.9686	normal_pair_count=39	position_range=6	repeat_name=IS150	repeat_name_evidence=junction	repeat_size_difference=-46	score=73.8	seed_z_score=-52.55	shifted_pair_count=900	side_1_annotate_key=gene	side_1_gene_name=yfiH	side_1_gene_position=coding (666/732 nt)	side_1_gene_product=hypothetical protein	side_1_gene_strand=<	side_1_locus_tag=ECB_02483	side_2_annotate_key=gene	side_2_gene_name=yfiH	side_2_gene_position=coding (665/732 nt)	side_2_gene_product=hypothetical protein	side_2_gene_strand=<	side_2_locus_tag=ECB_02483	size_shift=-1397	size_shift_lower=-1414	size_shift_upper=-1391	total_pair_count=993
PD	1928	.	REL606	2899618	-1	REL606	2899619	1	ambiguous_pair_count=52	candidate_covering_count=1227	distinct_pair_count=399	frequency=0.9702	frequency_lower=0.9607	frequency_upper=0.9779	normal_pair_count=35	position_range=5	repeat_size_candidates=IS150,IS186	score=95.4	seed_z_score=-59.37	shifted_pair_count=1140	side_1_annotate_key=gene	side_1_gene_name=yqeB	side_1_gene_position=coding (377/1626 nt)	side_1_gene_product=conserved protein with NAD(P)-binding Rossman fold	side_1_gene_strand=<	side_1_locus_tag=ECB_02708	side_2_annotate_key=gene	side_2_gene_name=yqeB	side_2_gene_position=coding (376/1626 nt)	side_2_gene_product=conserved protein with NAD(P)-binding Rossman fold	side_2_gene_strand=<	side_2_locus_tag=ECB_02708	size_shift=-1391	size_shift_lower=-1396	size_shift_upper=-1386	total_pair_count=1227
PD	1929	.	REL606	3268046	-1	REL606	3268047	1	ambiguous_pair_count=58	candidate_covering_count=1430	distinct_pair_count=503	frequency=0.9679	frequency_lower=0.9590	frequency_upper=0.9754	normal_pair_count=44	position_range=3	repeat_name=IS150	repeat_name_evidence=size	repeat_size_candidates=IS150	repeat_size_difference=36	score=112.5	seed_z_score=-64.30	shifted_pair_count=1328	side_1_annotate_key=gene	side_1_gene_name=yhbE/rpmA	side_1_gene_position=intergenic (-77/+50)	side_1_gene_product=conserved inner membrane protein/50S ribosomal protein L27	side_1_gene_strand=</<	side_1_locus_tag=ECB_03049/ECB_03050	side_2_annotate_key=gene	side_2_gene_name=yhbE/rpmA	side_2_gene_position=intergenic (-78/+49)	side_2_gene_product=conserved inner membrane protein/50S ribosomal protein L27	side_2_gene_strand=</<	side_2_locus_tag=ECB_03049/ECB_03050	size_shift=-1479	size_shift_lower=-1484	size_shift_upper=-1473	total_pair_count=1430
PD	1930	.	REL606	3511244	-1	REL606	3511245	1	ambiguous_pair_count=31	candidate_covering_count=1124	distinct_pair_count=378	frequency=0.9799	frequency_lower=0.9714	frequency_upper=0.9863	normal_pair_count=22	position_range=17	repeat_name=IS150	repeat_name_evidence=junction	repeat_size_difference=23	score=88.7	seed_z_score=-57.37	shifted_pair_count=1071	side_1_annotate_key=gene	side_1_gene_name=yhhX	side_1_gene_position=coding (1015/1038 nt)	side_1_gene_product=predicted oxidoreductase with NAD(P)-binding Rossmann-fold domain	side_1_gene_strand=<	side_1_locus_tag=ECB_03291	side_2_annotate_key=gene	side_2_gene_name=yhhX	side_2_gene_position=coding (1014/1038 nt)	side_2_gene_product=predicted oxidoreductase with NAD(P)-binding Rossmann-fold domain	side_2_gene_strand=<	side_2_locus_tag=ECB_03291	size_shift=-1466	size_shift_lower=-1482	size_shift_upper=-1455	total_pair_count=1124
PD	1931	.	REL606	3598266	-1	REL606	3598500	1	ambiguous_pair_count=1693	candidate_covering_count=1812	distinct_pair_count=13	frequency=0.2339	frequency_lower=0.1815	frequency_upper=0.2934	normal_pair_count=131	position_range=76	reject=PAIR_DISTANCE_SCORE	score=0.7	seed_z_score=13.01	shifted_pair_count=40	side_1_annotate_key=gene	side_1_gene_name=yhjA	side_1_gene_position=coding (414/1398 nt)	side_1_gene_product=predicted cytochrome C peroxidase	side_1_gene_strand=<	side_1_locus_tag=ECB_03366	side_2_annotate_key=gene	side_2_gene_name=yhjA	side_2_gene_position=coding (180/1398 nt)	side_2_gene_product=predicted cytochrome C peroxidase	side_2_gene_strand=<	side_2_locus_tag=ECB_03366	size_shift=233	size_shift_lower=223	size_shift_upper=238	total_pair_count=1864
PD	1932	.	REL606	3741577	-1	REL606	3741856	1	ambiguous_pair_count=1825	candidate_covering_count=1995	distinct_pair_count=42	frequency=0.5354	frequency_lower=0.4744	frequency_upper=0.5955	normal_pair_count=92	position_range=844	score=13.8	seed_z_score=25.16	shifted_pair_count=106	side_1_annotate_key=gene	side_1_gene_name=IS1	side_1_gene_position=noncoding (389/768 nt)	side_1_gene_product=repeat region	side_1_gene_strand=<	side_2_annotate_key=gene	side_2_gene_name=IS1	side_2_gene_position=noncoding (110/768 nt)	side_2_gene_product=repeat region	side_2_gene_strand=<	size_shift=278	size_shift_lower=268	size_shift_upper=280	total_pair_count=2023
PD	1933	.	REL606	3894257	-1	REL606	3900791	1	ambiguous_pair_count=0	candidate_covering_count=988	distinct_pair_count=358	frequency=0.9970	frequency_lower=0.9922	frequency_upper=0.9992	normal_pair_count=3	position_range=1346	score=78.7	seed_z_score=54.17	shifted_pair_count=985	side_1_annotate_key=gene	side_1_gene_name=IS150	side_1_gene_position=noncoding (704/1443 nt)	side_1_gene_product=repeat region	side_1_gene_strand=>	side_2_annotate_key=gene	side_2_gene_name=rbsR	side_2_gene_position=coding (758/993 nt)	side_2_gene_product=DNA-binding transcriptional repressor of ribose metabolism	side_2_gene_strand=>	side_2_locus_tag=ECB_03639	size_shift=6533	size_shift_lower=6525	size_shift_upper=6538	total_pair_count=988
PD	1934	.	REL606	4110229	-1	REL606	4110230	1	ambiguous_pair_count=50	candidate_covering_count=1000	distinct_pair_count=378	frequency=0.9895	frequency_lower=0.9822	frequency_upper=0.9943	normal_pair_count=10	position_range=13	repeat_size_candidates=IS150,IS186	score=78.5	seed_z_score=-54.12	shifted_pair_count=940	side_1_annotate_key=gene	side_1_gene_name=metL	side_1_gene_position=coding (288/2433 nt)	side_1_gene_product=bifunctional aspartate kinase II/homoserine dehydrogenase II	side_1_gene_strand=>	side_1_locus_tag=ECB_03826	side_2_annotate_key=gene	side_2_gene_name=metL	side_2_gene_position=coding (289/2433 nt)	side_2_gene_product=bifunctional aspartate kinase II/homoserine dehydrogenase II	side_2_gene_strand=>	side_2_locus_tag=ECB_03826	size_shift=-1449	size_shift_lower=-1450	size_shift_upper=-1428	total_pair_count=1000
PD	1935	.	REL606	4415705	-1	REL606	4415706	1	ambiguous_pair_count=66	candidate_covering_count=1004	distinct_pair_count=357	frequency=0.9925	frequency_lower=0.9860	frequency_upper=0.9965	normal_pair_count=7	position_range=7	repeat_name=IS150	repeat_name_evidence=junction	repeat_size_difference=35	score=79.1	seed_z_score=-54.31	shifted_pair_count=931	side_1_annotate_key=gene	side_1_gene_name=cycA	side_1_gene_position=coding (843/1413 nt)	side_1_gene_product=D-alanine/D-serine/glycine transporter	side_1_gene_strand=>	side_1_locus_tag=ECB_04080	side_2_annotate_key=gene	side_2_gene_name=cycA	side_2_gene_position=coding (844/1413 nt)	side_2_gene_product=D-alanine/D-serine/glycine transporter	side_2_gene_strand=>	side_2_locus_tag=ECB_04080	size_shift=-1478	size_shift_lower=-1486	size_shift_upper=-1473	total_pair_count=1004
PD	1936	.	REL606	4462444	-1	REL606	4462627	1	ambiguous_pair_count=3240	candidate_covering_count=2643	distinct_pair_count=7	frequency=0.1364	frequency_lower=0.0932	frequency_upper=0.1904	normal_pair_count=133	position_range=470	reject=PAIR_DISTANCE_FREQUENCY,PAIR_DISTANCE_SCORE	score=0.7	seed_z_score=13.05	shifted_pair_count=21	side_1_annotate_key=gene	side_1_gene_name=yjgM	side_1_gene_position=coding (127/504 nt)	side_1_gene_product=predicted acetyltransferase	side_1_gene_strand=<	side_1_locus_tag=ECB_04122	side_2_annotate_key=gene	side_2_gene_name=yjgM/yjgN	side_2_gene_position=intergenic (-57/-136)	side_2_gene_product=predicted acetyltransferase/conserved inner membrane protein	side_2_gene_strand=</>	side_2_locus_tag=ECB_04122/ECB_04123	size_shift=182	size_shift_lower=160	size_shift_upper=190	total_pair_count=3394
PD	1937	.	REL606	4613019	-1	REL606	4613020	1	ambiguous_pair_count=47	candidate_covering_count=739	distinct_pair_count=276	frequency=0.9971	frequency_lower=0.9909	frequency_upper=0.9995	normal_pair_count=2	position_range=1	repeat_name=IS150	repeat_name_evidence=junction	repeat_size_difference=-65	score=56.9	seed_z_score=-46.52	shifted_pair_count=690	side_1_annotate_key=gene	side_1_gene_name=smp/serB	side_1_gene_position=intergenic (-16/-90)	side_1_gene_product=hypothetical protein/3-phosphoserine phosphatase	side_1_gene_strand=</>	side_1_locus_tag=ECB_04263/ECB_04264	side_2_annotate_key=gene	side_2_gene_name=smp/serB	side_2_gene_position=intergenic (-17/-89)	side_2_gene_product=hypothetical protein/3-phosphoserine phosphatase	side_2_gene_strand=</>	side_2_locus_tag=ECB_04263/ECB_04264	size_shift=-1378	size_shift_lower=-1408	size_shift_upper=-1367	total_pair_count=739
