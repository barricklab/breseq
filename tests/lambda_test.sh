#rm -r output/lambda
#breseq.pl --no-junction -o output/lambda -r data/lambda/lambda.gbk data/lambda/lambda_mixed_population.fastq 

breseq.pl -o output/lambda_multiple_reference \
-r data/lambda/lambda-0.gbk \
-r data/lambda/lambda-1.gbk \
-r data/lambda/lambda-2.gbk \
-r data/lambda/lambda-3.gbk \
-r data/lambda/lambda-4.gbk \
data/lambda/lambda_mixed_population.fastq 