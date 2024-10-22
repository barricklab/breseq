BRESEQVERSIONSTRING=`perl -ne 's/AC_INIT\(\[(.+?)\].+?\[(.+?)\].+/\1-\2/ && print' configure.ac`

./bootstrap.sh
BINARYPLATFORM=`uname`
BINARYARCH=`arch`
ARCHFLAGS=""
CONFIGFLAGS=""
if [ "$BINARYPLATFORM" == "Darwin" ]; then
	BINARYARCH="universal"
	ARCHFLAGS="-mmacosx-version-min=10.13 -arch arm64 -arch x86_64"
	#ARCHFLAGS="-mmacosx-version-min=10.13 -arch x86_64"

	BINARYPLATFORM="MacOSX-10.13+"
	
	#Libz scenarios
	
	# 1) Using macports
	# You need to make sure this is a +universal build of libz for this to work
  	# $ port install libz +univsersal
	
	#For static zlib installed with macports
  	#CONFIGFLAGS="--with-static-libz=/opt/local"
  	
	# 2) Using own install
  	
  	# Even better is installing it from source with the same options used here
  	#
  	# From within the libz source code directory:
  	# 
  	# $ CFLAGS="-mmacosx-version-min=10.13" ./configure --static --prefix=$HOME/local --archs="-arch arm64 -arch x86_64" 
  	# $ make install
  	
  	CONFIGFLAGS="--with-static-libz=$HOME/local"
  	
fi
BINARYNAME=${BINARYPLATFORM}-${BINARYARCH}


BINARYLOCALDIR=${BRESEQVERSIONSTRING}-${BINARYNAME}
BINARYDIR=${PWD}/${BINARYLOCALDIR}
rm -rf ${BINARYDIR} ${BINARYDIR}.tar.gz ${BINARYDIR}.tgz

echo "${BINARYLOCALDIR}"
echo "${BINARYDIR}"
echo "./configure --without-libunwind --prefix=\"${BINARYDIR}\" --disable-shared --enable-static CFLAGS=\"${ARCHFLAGS} ${CFLAGS}\" CXXFLAGS=\"${ARCHFLAGS} ${CXXFLAGS}\" LDFLAGS=\"${ARCHFLAGS} ${LDFLAGS}\" $CONFIGFLAGS"
./configure --without-libunwind --prefix="${BINARYDIR}" CFLAGS="${ARCHFLAGS} ${CFLAGS}" CXXFLAGS="${ARCHFLAGS} ${CXXFLAGS}" LDFLAGS="${ARCHFLAGS} ${LDFLAGS}" $CONFIGFLAGS


make clean
make -j 12
#make test
make install

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

#need to update the #COMMAND line of the expected GenomeDiff to match path this is run from
sed -i -e 's/.\/src\/c\/breseq/bin/g' ${BINARYDIR}/tests/*/expected.gd

mkdir -p ${BINARYDIR}/tests/data/lambda
cp tests/data/lambda/empty.fastq ${BINARYDIR}/tests/data/lambda
cp tests/data/lambda/only_bad.fastq ${BINARYDIR}/tests/data/lambda
cp tests/data/lambda/lambda_mixed_population.1.fastq ${BINARYDIR}/tests/data/lambda
cp tests/data/lambda/lambda_mixed_population.2.fastq ${BINARYDIR}/tests/data/lambda
cp tests/data/lambda/lambda_mixed_population.3.fastq ${BINARYDIR}/tests/data/lambda
cp tests/data/lambda/lambda_mixed_population.4.fastq ${BINARYDIR}/tests/data/lambda
cp tests/data/lambda/lambda_mixed_population.5.fastq ${BINARYDIR}/tests/data/lambda
cp tests/data/lambda/lambda.1-2.gbk ${BINARYDIR}/tests/data/lambda
cp tests/data/lambda/lambda.3.gbk ${BINARYDIR}/tests/data/lambda
cp tests/data/lambda/lambda.4.gbk ${BINARYDIR}/tests/data/lambda
cp tests/data/lambda/lambda.5.gbk ${BINARYDIR}/tests/data/lambda

#options here need to match those in Makefile.am for test to pass
echo "export TESTBINPREFIX=bin" > ${BINARYDIR}/tests/test.config;
echo "export BRESEQ_DATA_PATH=share/breseq" >> ${BINARYDIR}/tests/test.config;
echo "export BRESEQ_TEST_THREAD_ARG=\"-j 4\"" >> ${BINARYDIR}/tests/test.config

echo "tests/test.sh clean tests" > ${BINARYDIR}/run_tests.sh
echo "tests/test.sh test tests" >> ${BINARYDIR}/run_tests.sh

#Fix permssions to 644 and add back executables
chmod 644 $(find ${BINARYDIR} -type f) 
chmod a+x ${BINARYDIR}/run_tests.sh
chmod a+x ${BINARYDIR}/bin/*
chmod a+x ${BINARYDIR}/tests/*.sh

tar -czf ${BINARYLOCALDIR}.tar.gz ${BINARYLOCALDIR}
#rm -r ${BINARYLOCALDIR}
