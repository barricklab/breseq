BRESEQ=./src/perl/breseq.pl

all:

check: lambda_test
#	python tests/_testrunner/testrunner.py $@ ---testrunner-name=./run_tests

lambda_test:
	${BRESEQ} --no-junction \
		-o ${TEST_OUTPUT}/lambda \
		-r ${TEST_INPUT}/lambda/lambda.gbk \
		${TEST_INPUT}/lambda/lambda_mixed_population.fastq
