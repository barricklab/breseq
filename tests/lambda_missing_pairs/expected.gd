#=GENOME_DIFF	1.0
#=CREATED	11:08:22 08 Aug 2026
#=PROGRAM	breseq 0.50.0 revision 9365b7ce37ec
#=COMMAND	./src/breseq/breseq -j 4 --predict-missing-pairs -o ./tests/lambda_missing_pairs -r ./tests/lambda_missing_pairs/../data/lambda/lambda.gbk ./tests/lambda_missing_pairs/output.simulated_1.fastq ./tests/lambda_missing_pairs/output.simulated_2.fastq
#=REFSEQ	./tests/lambda_missing_pairs/../data/lambda/lambda.gbk
#=READSEQ	./tests/lambda_missing_pairs/output.simulated_1.fastq
#=READSEQ	./tests/lambda_missing_pairs/output.simulated_2.fastq
#=CONVERTED-BASES	2071000
#=CONVERTED-READS	20710
#=INPUT-BASES	2071000
#=INPUT-READS	20710
#=MAPPED-BASES	1940700
#=MAPPED-READS	19407
MC	1	.	NC_001416	1	97	0	0	gene_name=–/nu1	gene_position=intergenic (–/-94)	gene_product=–/DNA packaging protein	gene_strand=–/>	left_inside_cov=0	left_outside_cov=NA	locus_tag=–/lambdap01	right_inside_cov=20	right_outside_cov=21
MC	2	.	NC_001416	20661	20661	0	0	gene_name=orf-401|orf206b	gene_position=coding (1012/1206 nt)|coding (107/621 nt)	gene_product=Tail fiber protein|hypothetical protein	gene_strand=>|<	left_inside_cov=0	left_outside_cov=43	locus_tag=lambdap27|lambdap90	right_inside_cov=0	right_outside_cov=44	snp_type=|
MC	3	.	NC_001416	47317	47317	0	0	gene_name=lambdap78	gene_position=coding (259/534 nt)	gene_product=putative envelope protein	gene_strand=<	left_inside_cov=0	left_outside_cov=37	locus_tag=lambdap78	right_inside_cov=0	right_outside_cov=37
MC	4	.	NC_001416	48422	48502	0	0	gene_name=lambdap79/–	gene_position=intergenic (+478/–)	gene_product=hypothetical protein/–	gene_strand=>/–	left_inside_cov=19	left_outside_cov=21	locus_tag=lambdap79/–	right_inside_cov=0	right_outside_cov=NA
UN	5	.	NC_001416	1	16
UN	6	.	NC_001416	20	20
UN	7	.	NC_001416	20661	20661
UN	8	.	NC_001416	23987	24002
UN	9	.	NC_001416	47317	47317
UN	10	.	NC_001416	48484	48502
MP	11	.	NC_001416	24000	-1	annotate_key=gene	candidate_unpaired_count=49	distinct_read_count=44	frequency=1.0000	frequency_lower=0.9407	frequency_upper=1.0000	gene_name=ea47/ea31	gene_position=intergenic (-82/+509)	gene_product=ea47/ea31	gene_strand=</<	locus_tag=lambdap80/lambdap81	spanning_pair_count=0	total_read_count=49	unpaired_read_count=49
MP	12	.	NC_001416	24002	1	annotate_key=gene	candidate_unpaired_count=42	distinct_read_count=40	frequency=1.0000	frequency_lower=0.9312	frequency_upper=1.0000	gene_name=ea47/ea31	gene_position=intergenic (-84/+507)	gene_product=ea47/ea31	gene_strand=</<	locus_tag=lambdap80/lambdap81	spanning_pair_count=0	total_read_count=42	unpaired_read_count=42
