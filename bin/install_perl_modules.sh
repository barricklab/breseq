## helper script to install Perl modules that require compilation

THISDIR=`dirname ${BASH_SOURCE}`
CURRENTDIR=`pwd`

## Bio::SamTools

## we need this header file to be in the same directory as libbam.a to install
cp lib/libbam.a src/c/samtools-0.1.7a

BIOSAMTOOLS="$CURRENTDIR/$THISDIR/../extern/Bio-SamTools-1.19"
cd $BIOSAMTOOLS
export SAMTOOLS=../../src/c/samtools-0.1.7a
echo $SAMTOOLS
perl Build.PL --install_base=../../src/perl/extern
./Build
./Build install

## Math::CDF
MATHCDF="$CURRENTDIR/$THISDIR/../extern/Math-CDF-0.1"
echo $MATHCDF
cd $MATHCDF
perl Makefile.PL PREFIX=../../src/perl/extern
make
make install

## Math::Pari
MATHPARI="$CURRENTDIR/$THISDIR/../extern/Math-Pari-2.01080604"
echo $MATHPARI
cd $MATHPARI
perl Makefile.PL PREFIX=../../src/perl/extern
make
make install
