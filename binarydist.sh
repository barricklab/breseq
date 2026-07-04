BRESEQVERSIONSTRING=`perl -ne 's/AC_INIT\(\[(.+?)\].+?\[(.+?)\].+/\1-\2/ && print' configure.ac`

# Pass --rebuild-prereqs to force a clean rebuild of the static-link
# prerequisite libraries (zlib, miniz, htslib, and htslib's CRAM-codec
# dependencies bzip2/xz) in binarydist_build/, even if already-built copies
# matching the pinned versions in dev-environment.yml are present.
REBUILD_PREREQS=0
for ARG in "$@"; do
	if [ "$ARG" == "--rebuild-prereqs" ]; then
		REBUILD_PREREQS=1
	fi
done

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

	# Discard any CFLAGS/CXXFLAGS/CPPFLAGS/LDFLAGS inherited from an active
	# conda environment (e.g. breseq-dev). They point into that env's
	# include/lib directories (-isystem .../env/include, -L.../env/lib,
	# -Wl,-rpath,.../env/lib), whose dylibs (e.g. libc++.dylib) are
	# single-architecture. The linker picks those up ahead of the system's
	# universal libc++, "ignores" them for the non-matching slice, and then
	# fails with missing C++ runtime symbols (___cxa_throw, etc.) for that
	# architecture — and even when it links, the resulting binary would carry
	# an rpath dependency on the conda env, defeating the point of a portable
	# static binary distribution.
	CFLAGS=""
	CXXFLAGS=""
	CPPFLAGS=""
	LDFLAGS=""

	# Force the system Apple clang for the prerequisite builds below: it
	# can produce true universal (arm64 + x86_64) binaries from a single
	# invocation via multiple -arch flags. Conda dev environments commonly
	# export CC/CXX pointing at single-architecture cross-compilers (e.g.
	# arm64-apple-darwin20.0.0-clang), which silently build only one slice
	# and cause "ld: symbol(s) not found for architecture ..." at link time.
	PREREQ_CC=/usr/bin/clang
	PREREQ_CXX=/usr/bin/clang++

	# Static-link prerequisites (zlib, miniz, htslib, plus htslib's CRAM
	# codec dependencies bzip2/xz, built so htslib retains CRAM support for
	# future development even though breseq itself only uses BAM) are built as
	# static-universal (arm64 + x86_64) libraries into binarydist_build/ under
	# the source root, so that the resulting breseq/gdtools binaries have
	# no runtime dependency on these libraries. Each is (re)built only when
	# missing or when its built version doesn't match the pin in
	# dev-environment.yml (or when --rebuild-prereqs is passed).
	PREREQDIR="${PWD}/binarydist_build"
	mkdir -p "${PREREQDIR}"

	get_pinned_version () {
		grep -E "^[[:space:]]*-[[:space:]]*$1=" dev-environment.yml | sed -E 's/.*=([0-9.]+).*/\1/'
	}

	needs_build () {
		LIBDIR="$1"
		PINNED_VERSION="$2"
		VERSION_STAMP="${PREREQDIR}/${LIBDIR}/.version"
		if [ "$REBUILD_PREREQS" -eq 1 ]; then
			return 0
		fi
		if [ ! -f "$VERSION_STAMP" ]; then
			return 0
		fi
		if [ "$(cat "$VERSION_STAMP")" != "$PINNED_VERSION" ]; then
			return 0
		fi
		return 1
	}

	# zlib
	ZLIB_VERSION=`get_pinned_version zlib`
	if needs_build zlib "$ZLIB_VERSION"; then
		echo "Building static universal zlib ${ZLIB_VERSION} into ${PREREQDIR}/zlib"
		rm -rf "${PREREQDIR}/zlib" "${PREREQDIR}/zlib-src"
		git clone --branch "v${ZLIB_VERSION}" --depth 1 https://github.com/madler/zlib "${PREREQDIR}/zlib-src"
		( cd "${PREREQDIR}/zlib-src" \
			&& CC="${PREREQ_CC}" CFLAGS="-mmacosx-version-min=10.13" ./configure --static --prefix="${PREREQDIR}/zlib" --archs="-arch arm64 -arch x86_64" \
			&& make && make install ) || exit 1
		rm -rf "${PREREQDIR}/zlib-src"
		echo "$ZLIB_VERSION" > "${PREREQDIR}/zlib/.version"
	fi

	# miniz
	MINIZ_VERSION=`get_pinned_version miniz`
	if needs_build miniz "$MINIZ_VERSION"; then
		echo "Building static universal miniz ${MINIZ_VERSION} into ${PREREQDIR}/miniz"
		rm -rf "${PREREQDIR}/miniz" "${PREREQDIR}/miniz-src"
		git clone --branch "${MINIZ_VERSION}" --depth 1 https://github.com/richgel999/miniz "${PREREQDIR}/miniz-src"
		( mkdir -p "${PREREQDIR}/miniz-src/build" && cd "${PREREQDIR}/miniz-src/build" \
			&& cmake .. -DBUILD_SHARED_LIBS=OFF -DINSTALL_PROJECT=ON \
				-DBUILD_TESTS=OFF -DBUILD_EXAMPLES=OFF \
				-DCMAKE_C_COMPILER="${PREREQ_CC}" \
				-DCMAKE_CXX_COMPILER="${PREREQ_CXX}" \
				-DCMAKE_INSTALL_PREFIX="${PREREQDIR}/miniz" \
				-DCMAKE_OSX_ARCHITECTURES="arm64;x86_64" \
				-DCMAKE_OSX_DEPLOYMENT_TARGET=10.13 \
			&& make && make install ) || exit 1
		rm -rf "${PREREQDIR}/miniz-src"
		echo "$MINIZ_VERSION" > "${PREREQDIR}/miniz/.version"
	fi

	# bzip2 (CRAM codec dependency of htslib)
	BZ2_VERSION=`get_pinned_version bzip2`
	if needs_build bz2 "$BZ2_VERSION"; then
		echo "Building static universal bzip2 ${BZ2_VERSION} into ${PREREQDIR}/bz2"
		rm -rf "${PREREQDIR}/bz2" "${PREREQDIR}/bz2-src"
		git clone --branch "bzip2-${BZ2_VERSION}" --depth 1 https://github.com/libarchive/bzip2 "${PREREQDIR}/bz2-src"
		( cd "${PREREQDIR}/bz2-src" \
			&& make CC="${PREREQ_CC}" CFLAGS="-mmacosx-version-min=10.13 -arch arm64 -arch x86_64 -O2" libbz2.a \
			&& mkdir -p "${PREREQDIR}/bz2/include" "${PREREQDIR}/bz2/lib" \
			&& cp bzlib.h "${PREREQDIR}/bz2/include" \
			&& cp libbz2.a "${PREREQDIR}/bz2/lib" ) || exit 1
		rm -rf "${PREREQDIR}/bz2-src"
		echo "$BZ2_VERSION" > "${PREREQDIR}/bz2/.version"
	fi

	# xz / liblzma (CRAM codec dependency of htslib)
	XZ_VERSION=`get_pinned_version xz`
	if needs_build xz "$XZ_VERSION"; then
		echo "Building static universal xz/liblzma ${XZ_VERSION} into ${PREREQDIR}/xz"
		rm -rf "${PREREQDIR}/xz" "${PREREQDIR}/xz-src"
		mkdir -p "${PREREQDIR}/xz-src"
		# Use the release tarball (not a git clone) because it ships a
		# pre-generated ./configure — the git checkout requires running
		# autogen.sh, which depends on gettext/autopoint that aren't part
		# of dev-environment.yml.
		curl -sL "https://github.com/tukaani-project/xz/releases/download/v${XZ_VERSION}/xz-${XZ_VERSION}.tar.gz" \
			| tar -xz -C "${PREREQDIR}/xz-src" --strip-components=1
		( cd "${PREREQDIR}/xz-src" \
			&& CC="${PREREQ_CC}" CFLAGS="-mmacosx-version-min=10.13 -arch arm64 -arch x86_64" \
			   ./configure --disable-shared --enable-static --prefix="${PREREQDIR}/xz" \
			       --disable-doc --disable-scripts --disable-nls \
			       --disable-xz --disable-xzdec --disable-lzmadec --disable-lzmainfo \
			&& make && make install ) || exit 1
		rm -rf "${PREREQDIR}/xz-src"
		echo "$XZ_VERSION" > "${PREREQDIR}/xz/.version"
	fi

	# htslib
	HTSLIB_VERSION=`get_pinned_version htslib`
	if needs_build htslib "$HTSLIB_VERSION"; then
		echo "Building static universal htslib ${HTSLIB_VERSION} into ${PREREQDIR}/htslib"
		rm -rf "${PREREQDIR}/htslib" "${PREREQDIR}/htslib-src"
		git clone --branch "${HTSLIB_VERSION}" --depth 1 --shallow-submodules --recurse-submodules https://github.com/samtools/htslib "${PREREQDIR}/htslib-src"
		( cd "${PREREQDIR}/htslib-src" \
			&& autoreconf -i \
			&& CC="${PREREQ_CC}" \
			   CFLAGS="-mmacosx-version-min=10.13 -arch arm64 -arch x86_64" \
			   CPPFLAGS="-I${PREREQDIR}/zlib/include -I${PREREQDIR}/bz2/include -I${PREREQDIR}/xz/include" \
			   LDFLAGS="-L${PREREQDIR}/zlib/lib -L${PREREQDIR}/bz2/lib -L${PREREQDIR}/xz/lib" \
			   ./configure --disable-libcurl --prefix="${PREREQDIR}/htslib" \
			&& make && make install ) || exit 1
		rm -rf "${PREREQDIR}/htslib-src"
		# htslib's ./configure has no --disable-shared (it warns "unrecognized
		# options") and always builds+installs both libhts.a and libhts.dylib;
		# the latter comes out single-architecture (htslib's dylib link rule
		# doesn't pass through our -arch flags), so the linker would pick it
		# over our universal .a and fail with "symbol(s) not found for
		# architecture x86_64". Removing it leaves only the static archive,
		# matching every other prereq's lib/ directory.
		rm -f "${PREREQDIR}/htslib/lib/"libhts.*dylib
		echo "$HTSLIB_VERSION" > "${PREREQDIR}/htslib/.version"
	fi

	# Force the same system Apple clang for the main breseq/gdtools build:
	# conda dev environments commonly export CC/CXX as single-architecture
	# cross-compilers (e.g. arm64-apple-darwin20.0.0-clang++), which silently
	# emit incomplete x86_64 object code (missing libc++ runtime symbols like
	# ___cxa_throw) and fail to link a true universal binary even when given
	# -arch arm64 -arch x86_64.
	CONFIGFLAGS="CC=${PREREQ_CC} CXX=${PREREQ_CXX}"

	# Point the build at our just-(re)built static-universal prerequisites
	# instead of whatever libz/miniz/htslib the active conda environment (or
	# system) provides — configure.ac always emits -I${prefix}/include and
	# -lz/-lminiz/-lhts (it no longer has --with-static-* options); putting
	# our -I/-L paths first makes those resolve to our universal .a files
	# (each of these lib/ dirs holds only a .a, no .dylib, so -lfoo can't
	# pick up a wrong-architecture dynamic library by mistake). bzip2/xz
	# are purely a transitive CRAM-codec dependency of htslib with no
	# -lfoo reference anywhere in Makefile.am, so their static archives are
	# passed via the dedicated CRAM_CODEC_LIBS configure variable, which
	# Makefile.am appends to breseq_LDADD/gdtools_LDADD (after $(HTSLIB_LINK),
	# i.e. exactly where libhts.a's undefined bz2/lzma symbols need them) —
	# kept out of the global LIBS so libtool doesn't also see these static
	# archives while linking the libbreseq.la convenience library and emit
	# "Linking the shared library ... against the static library ... is not
	# portable" warnings.
	STATIC_PREREQ_CPPFLAGS="-I${PREREQDIR}/zlib/include -I${PREREQDIR}/miniz/include -I${PREREQDIR}/htslib/include -I${PREREQDIR}/bz2/include -I${PREREQDIR}/xz/include"
	STATIC_PREREQ_LDFLAGS="-L${PREREQDIR}/zlib/lib -L${PREREQDIR}/miniz/lib -L${PREREQDIR}/htslib/lib"
	CRAM_CODEC_LIBS="${PREREQDIR}/bz2/lib/libbz2.a ${PREREQDIR}/xz/lib/liblzma.a"

fi
BINARYNAME=${BINARYPLATFORM}-${BINARYARCH}


BINARYLOCALDIR=${BRESEQVERSIONSTRING}-${BINARYNAME}
BINARYDIR=${PWD}/${BINARYLOCALDIR}
rm -rf ${BINARYDIR} ${BINARYDIR}.tar.gz ${BINARYDIR}.tgz

echo "${BINARYLOCALDIR}"
echo "${BINARYDIR}"
echo "./configure --prefix=\"${BINARYDIR}\" --disable-shared --enable-static CFLAGS=\"${ARCHFLAGS} ${CFLAGS}\" CXXFLAGS=\"${ARCHFLAGS} ${CXXFLAGS}\" CPPFLAGS=\"${STATIC_PREREQ_CPPFLAGS} ${CPPFLAGS}\" LDFLAGS=\"${ARCHFLAGS} ${STATIC_PREREQ_LDFLAGS} ${LDFLAGS}\" LIBS=\"${LIBS}\" CRAM_CODEC_LIBS=\"${CRAM_CODEC_LIBS}\" $CONFIGFLAGS"
./configure --prefix="${BINARYDIR}" CFLAGS="${ARCHFLAGS} ${CFLAGS}" CXXFLAGS="${ARCHFLAGS} ${CXXFLAGS}" CPPFLAGS="${STATIC_PREREQ_CPPFLAGS} ${CPPFLAGS}" LDFLAGS="${ARCHFLAGS} ${STATIC_PREREQ_LDFLAGS} ${LDFLAGS}" LIBS="${LIBS}" CRAM_CODEC_LIBS="${CRAM_CODEC_LIBS}" $CONFIGFLAGS


make clean
make -j 12
#make test
make install

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

echo "tests/test.sh clean tests" > ${BINARYDIR}/run_tests.sh
echo "tests/test.sh test tests" >> ${BINARYDIR}/run_tests.sh

#Fix permissions to 644 and add back executables
chmod 644 $(find ${BINARYDIR} -type f) 
chmod a+x ${BINARYDIR}/run_tests.sh
chmod a+x ${BINARYDIR}/bin/*
chmod a+x ${BINARYDIR}/tests/*.sh

tar -czf ${BINARYLOCALDIR}.tar.gz ${BINARYLOCALDIR}
rm -r ${BINARYLOCALDIR}
