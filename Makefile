BRESEQ=./src/perl/breseq.pl
TEST_INPUT=./tests/data
TEST_OUTPUT=./tests/output

all:

check: lambda_test

lambda_test:
	${BRESEQ} --no-junction \
		-o ${TEST_OUTPUT}/lambda \
		-r ${TEST_INPUT}/lambda/lambda.gbk \
		${TEST_INPUT}/lambda/lambda_mixed_population.fastq
