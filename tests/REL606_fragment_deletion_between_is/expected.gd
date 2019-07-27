#=GENOME_DIFF	1.0
#=COMMAND	./src/c/breseq/breseq -j 4 -b 0 -o tests/REL606_fragment_deletion_between_is -r tests/REL606_fragment_deletion_between_is/../data/REL606/REL606.fragment.gbk tests/REL606_fragment_deletion_between_is/../data/REL606/REL606.fragment.2.fastq
#=REFSEQ	tests/REL606_fragment_deletion_between_is/../data/REL606/REL606.fragment.gbk
#=READSEQ	tests/REL606_fragment_deletion_between_is/../data/REL606/REL606.fragment.2.fastq
#=CONVERTED-BASES	866700
#=CONVERTED-READS	24075
#=INPUT-BASES	866700
#=INPUT-READS	24075
#=MAPPED-BASES	864799
#=MAPPED-READS	24028
DEL	1	2	REL606-5	10535	21535	between=IS1	gene_name=[ECB_00212]–[phoE]	gene_product=[ECB_00212],yhhI,yafV,ykfE,fadE,lpcA,yafJ,yafK,yafQ,dinJ,yafL,yafM,fhiA,mbhA,dinB,yafN,yafO,yafP,ykfJ,prfH,pepD,gpt,yafA,crl,insA-3,insB-3,[phoE]	genes_inactivated=ECB_00212,yhhI,yafV,ykfE,fadE,lpcA,yafJ,yafK,yafQ,dinJ,yafL,yafM,fhiA,mbhA,dinB,yafN,yafO,yafP,ykfJ,prfH,pepD,gpt,yafA,crl,insA-3,insB-3	genes_overlapping=phoE	locus_tag=[ECB_00212]–[ECB_00238]	locus_tags_inactivated=ECB_00212,ECB_00213,ECB_00214,ECB_00215,ECB_00216,ECB_00217,ECB_00218,ECB_00219,ECB_00220,ECB_00221,ECB_00222,ECB_00223,ECB_00224,ECB_00225,ECB_00226,ECB_00227,ECB_00228,ECB_00229,ECB_00230,ECB_00231,ECB_00232,ECB_00233,ECB_00234,ECB_00235,ECB_00236,ECB_00237	locus_tags_overlapping=ECB_00238	mutation_category=large_deletion	position_end=32069	position_start=10535	ref_seq=21535-bp
MC	2	.	REL606-5	10229	32044	305	298	gene_name=[insB-2]–[insB-3]	gene_product=[insB-2],ECB_00212,yhhI,yafV,ykfE,fadE,lpcA,yafJ,yafK,yafQ,dinJ,yafL,yafM,fhiA,mbhA,dinB,yafN,yafO,yafP,ykfJ,prfH,pepD,gpt,yafA,crl,insA-3,[insB-3]	left_inside_cov=15	left_outside_cov=18	locus_tag=[ECB_00211]–[ECB_00237]	right_inside_cov=15	right_outside_cov=17
JC	3	.	REL606-5	1	1	REL606-5	46298	-1	0	alignment_overlap=0	circular_chromosome=1	coverage_minus=14	coverage_plus=13	flanking_left=36	flanking_right=36	frequency=1	junction_possible_overlap_registers=35	key=REL606-5__1__1__REL606-5__46298__-1__0____36__36__0__0	max_left=34	max_left_minus=33	max_left_plus=34	max_min_left=16	max_min_left_minus=15	max_min_left_plus=16	max_min_right=18	max_min_right_minus=18	max_min_right_plus=17	max_pos_hash_score=70	max_right=30	max_right_minus=29	max_right_plus=30	neg_log10_pos_hash_p_value=0.0	new_junction_coverage=0.79	new_junction_read_count=27	polymorphism_frequency=1.000e+00	pos_hash_score=21	prediction=consensus	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=0.00	side_1_overlap=0	side_1_possible_overlap_registers=35	side_1_read_count=0	side_1_redundant=0	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=0.00	side_2_overlap=0	side_2_possible_overlap_registers=35	side_2_read_count=0	side_2_redundant=0	total_non_overlap_reads=27
UN	4	.	REL606-5	389	408
UN	5	.	REL606-5	8392	8416
UN	6	.	REL606-5	9797	10001
UN	7	.	REL606-5	10059	10181
UN	8	.	REL606-5	10244	32036
