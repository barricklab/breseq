BRESEQ=./src/perl/breseq.pl
TEST_INPUT=./tests/data
TEST_OUTPUT=./tests/output

all:

check: test

test:  lambda_test \
       lambda_test_multiple_ref_seqs \
       lambda_test_multiple_read_files \
       lambda_test_multiple_ref_seq_and_read_files

clean : 
	rm -r ${TEST_OUTPUT}

lambda_test:
	${BRESEQ} -a --no-junction \
		-o ${TEST_OUTPUT}/lambda \
		-r ${TEST_INPUT}/lambda/lambda.gbk \
		${TEST_INPUT}/lambda/lambda_mixed_population.fastq

lambda_test_multiple_ref_seqs:
	${BRESEQ} -a \
		-o ${TEST_OUTPUT}/lambda_multiple_ref_seqs \
		-r ${TEST_INPUT}/lambda/lambda.1.gbk \
		-r ${TEST_INPUT}/lambda/lambda.2.gbk \
		-r ${TEST_INPUT}/lambda/lambda.3.gbk \
		-r ${TEST_INPUT}/lambda/lambda.4.gbk \
		-r ${TEST_INPUT}/lambda/lambda.5.gbk \
		${TEST_INPUT}/lambda/lambda_mixed_population.fastq

lambda_test_multiple_read_files:
	${BRESEQ} -a --no-junction \
		-o ${TEST_OUTPUT}/lambda_multiple_read_files \
		-r ${TEST_INPUT}/lambda/lambda.gbk \
		${TEST_INPUT}/lambda/lambda_mixed_population.1.fastq \
		${TEST_INPUT}/lambda/lambda_mixed_population.2.fastq \
		${TEST_INPUT}/lambda/lambda_mixed_population.3.fastq \
		${TEST_INPUT}/lambda/lambda_mixed_population.4.fastq \
		${TEST_INPUT}/lambda/lambda_mixed_population.5.fastq

lambda_test_multiple_ref_seq_and_read_files:
	${BRESEQ} -a \
		-o ${TEST_OUTPUT}/lambda_multiple_ref_seq_and_read_files \
		-r ${TEST_INPUT}/lambda/lambda.1.gbk \
		-r ${TEST_INPUT}/lambda/lambda.2.gbk \
		-r ${TEST_INPUT}/lambda/lambda.3.gbk \
		-r ${TEST_INPUT}/lambda/lambda.4.gbk \
		-r ${TEST_INPUT}/lambda/lambda.5.gbk \
		${TEST_INPUT}/lambda/lambda_mixed_population.1.fastq \
		${TEST_INPUT}/lambda/lambda_mixed_population.2.fastq \
		${TEST_INPUT}/lambda/lambda_mixed_population.3.fastq \
		${TEST_INPUT}/lambda/lambda_mixed_population.4.fastq \
		${TEST_INPUT}/lambda/lambda_mixed_population.5.fastq