## helper script to install Perl modules that require compilation

THISDIR=`dirname ${BASH_SOURCE}`
CURRENTDIR=`pwd`

## Bio::SamTools
BIOSAMTOOLS="$CURRENTDIR/$THISDIR/../extern/Bio-SamTools-1.19"
cd $BIOSAMTOOLS
export SAMTOOLS=../../src/c/samtools-0.1.7a
echo $SAMTOOLS
perl Build.PL --install_base=../..
./Build
./Build install

## Math::CDF
MATHCDF="$CURRENTDIR/$THISDIR/../extern/Math-CDF-0.1"
echo $MATHCDF
cd $MATHCDF
perl Makefile.PL PREFIX=../..
make
make install
