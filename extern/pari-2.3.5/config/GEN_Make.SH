dir=examples
name=Makefile.$osname-$arch
file=$dir/$name
lnfile=Makefile

echo Extracting $file
rm -f $file $dir/$lnfile

if test -z "$DLLD"; then static=y; fi
case "$static" in
  n) dft=dyn ;;
  y) dft=sta ;;
esac

RUNPTH=
if test -n "$runpathprefix"; then
  RUNPTH=`config/ldflags "$LDneedsWl" "$runpathprefix $runpath"`
fi

cat > $file << EOT
# Generic Makefile for PARI programs -- $pretty
#
#  This file was created by Configure. Any change made to it will be
#  lost when Configure is run.
#
# make all will create
#  extgcd-dyn (linked dynamically with libpari)
#  extgcd-sta (linked statically)
#  libextgcd.so (to be used by "install" under GP)
#
# Under GP: install("extgcd", "GG&&", "gcdex", "./libextgcd.so") enables
# you to subsequently use gcdex to call extgcd (see the reference manual).
#

# change this TARGET to compile your own programs
TARGET = extgcd
SHELL  = $make_sh

DBGFLAGS   = $DBGFLAGS $CFLAGS
CFLAGS     = $OPTFLAGS $CFLAGS
#CFLAGS    = \$(DBGFLAGS)

# Various linkers use different flags to force static compilation. Choose
# the one which is relevant for your installation.
#
# Solaris ld (global)
#STATIC    = -dn

# Solaris ld (toggle: no shared object accepted until -B dynamic is seen
#STATIC    = -B static

# gcc
#STATIC    = -static

CC         = $CC
CPPFLAGS   = -I. -I$includedir
LD         = $LD
LDFLAGS    = $LDFLAGS
MODLD       = $MODLD
MODLDFLAGS  = $MODLDFLAGS
EXTRAMODLDFLAGS = $EXTRAMODLDFLAGS
EXTRALIBS  =

RUNPTH     = $RUNPTH
DLCFLAGS   = $DLCFLAGS
LIBS       = $LIBS -L$libdir -lpari

RM = rm -f


OBJS = \$(TARGET).o
DYN = lib\$(TARGET).$DLSUFFIX
ALL = \$(TARGET)-sta \$(TARGET)-dyn \$(DYN)

dft: \$(TARGET)-$dft

all: \$(ALL)

sta: \$(TARGET)-sta

dyn: \$(TARGET)-dyn

dynlib: \$(DYN)

\$(DYN): \$(OBJS)
	\$(MODLD) -o \$@ \$(MODLDFLAGS) \$(OBJS) \$(EXTRAMODLDFLAGS)

\$(TARGET)-sta: \$(OBJS)
	\$(LD) -o \$@ \$(LDFLAGS) \$< \$(EXTRALIBS) \$(STATIC) \$(LIBS)

\$(TARGET)-dyn: \$(OBJS)
	\$(LD) -o \$@ \$(LDFLAGS) \$< \$(RUNPTH) \$(LIBS)

%.o: %.c
	\$(CC) -c \$(CFLAGS) \$(CPPFLAGS) \$(DLCFLAGS) \$<
clean:
	-\$(RM) *.o \$(ALL)
EOT

( cd $dir ; $ln_s $name $lnfile )
