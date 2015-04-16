BRESEQVERSIONSTRING=`perl -ne 's/AC_INIT\(\[(.+?)\].+?\[(.+?)\].+/\1-\2/ && print' configure.ac`

./bootstrap.sh
BINARYPLATFORM=`uname`
BINARYARCH=`arch`
MYCFLAGS="$CFLAGS"
BINARYNAME=${BINARYPLATFORM}-${BINARYARCH}
if [ "$BINARYPLATFORM" == "Darwin" ]; then
	BINARYARCH="universal"
	MYCFLAGS="-arch i386 -arch x86_64"
	BINARYNAME="MacOSX"
fi

BINARYLOCALDIR=${BRESEQVERSIONSTRING}-${BINARYNAME}
BINARYDIR=${PWD}/${BINARYLOCALDIR}
rm -r ${BINARYDIR} ${BINARYDIR}.tgz
./configure CFLAGS="${MYCFLAGS}" CXXFLAGS="${MYCFLAGS}" LDFLAGS="${MYCFLAGS}" --prefix="${BINARYDIR}"
make clean
make -j 6 install

#Documentation and information
make docs
cp -r src/doc/_build/html ${BINARYDIR}/documentation
cp -r LICENSE ${BINARYDIR}
cp -r README-BINARY ${BINARYDIR}/README

#Test
mkdir -p ${BINARYDIR}/tests/lambda_mult_ref_read
cp tests/common.sh ${BINARYDIR}/tests
cp tests/test.sh ${BINARYDIR}/tests
cp tests/lambda_mult_ref_read/expected.gd ${BINARYDIR}/tests/lambda_mult_ref_read
cp tests/lambda_mult_ref_read/testcmd.sh ${BINARYDIR}/tests/lambda_mult_ref_read
mkdir -p ${BINARYDIR}/tests/data/lambda
cp tests/data/lambda/lambda_mixed_population.1.fastq ${BINARYDIR}/tests/data/lambda
cp tests/data/lambda/lambda_mixed_population.2.fastq ${BINARYDIR}/tests/data/lambda
cp tests/data/lambda/lambda_mixed_population.3.fastq ${BINARYDIR}/tests/data/lambda
cp tests/data/lambda/lambda_mixed_population.4.fastq ${BINARYDIR}/tests/data/lambda
cp tests/data/lambda/lambda_mixed_population.5.fastq ${BINARYDIR}/tests/data/lambda
cp tests/data/lambda/lambda.1-2.gbk ${BINARYDIR}/tests/data/lambda
cp tests/data/lambda/lambda.3.gbk ${BINARYDIR}/tests/data/lambda
cp tests/data/lambda/lambda.4.gbk ${BINARYDIR}/tests/data/lambda
cp tests/data/lambda/lambda.5.gbk ${BINARYDIR}/tests/data/lambda

echo "export TESTBINPREFIX=bin" > ${BINARYDIR}/tests/test.config;
echo "export BRESEQ_DATA_PATH=share/breseq" >> ${BINARYDIR}/tests/test.config;
echo "export BRESEQ_SAMTOOLS_PATH=bin" >> ${BINARYDIR}/tests/test.config;

echo "tests/test.sh clean tests" > ${BINARYDIR}/run_tests.sh
echo "tests/test.sh test tests" >> ${BINARYDIR}/run_tests.sh
chmod a+x ${BINARYDIR}/run_tests.sh

tar -czf ${BINARYLOCALDIR}.tar.gz ${BINARYLOCALDIR}
rm -r ${BINARYLOCALDIR}
