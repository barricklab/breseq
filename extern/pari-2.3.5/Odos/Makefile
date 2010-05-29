# Makefile for Pari/GP -- ix86 running DOS (ix86 kernel)
# This file was created by hand. Update it when config/ evolves
# $Id: Makefile 1320 2001-03-26 22:20:24Z karim $
SHELL    = /bin/sh

CFLAGS   = -O3 -DGCC_INLINE -Wall
#  Uncomment for debugging version
#CFLAGS   = -g -Wall
#XTRAOINL = $(OINLINE)

LDOBJS  =$(OBJS)
LDOBJSGP=$(OBJSGP)
#  Uncomment if the shell complains about long lines (bash will not, DOS will)
#LDOBJS  =@objs
#LDOBJSGP=@objsgp

#  Correct according to your local configuration. We expect to find
#    standard libraries    in $PREFIX/lib
#    standard headers      in $PREFIX/include
#    gcc                   in $PREFIXBIN/bin
PREFIX=/usr/local/cross-tools/i386-mingw32msvc
PREFIX=/usr/local/i386-pc-cygwin32/i386-pc-cygwin32
#PREFIX=/usr/local
PREFIXBIN=/usr/local/cross-tools
PREFIXBIN=$(PREFIX)

CC        = $(PREFIXBIN)/bin/gcc
LD        = $(CC)
CPP       = $(CC) -E -I.
CPPFLAGS  = -I. -I../src/headers -I$(PREFIX)/include
LDFLAGS   = -L$(PREFIX)/lib

#  Uncomment (and correct) if you use EMX
#EMXBIND   = emxbind -q C:/emx/bin/emxl.exe  gp gp.exe -C0
#LDFLAGS   = -Zsysv-signals

RLLIBS    = -lreadline
#  Uncomment if you use GPRL.ZIP or correct if you installed readline in a
#  non-standard place
#RLDIR     = ../readline
#RLINCLUDE = -I$(RLDIR)
#RLLIBS    = -L$(RLDIR) -lreadline -ltermcap

# use del if RM is not found.
RM = rm -f
#RM = del

# Graphic library.
PLOTFILE   = plotnull.c
PLOTCFLAGS =
PLOTLIBS   =

OBJS   =  kernel.o mp.o es.o alglin1.o alglin2.o arith1.o arith2.o base1.o base2.o base3.o base4.o base5.o bibli1.o bibli2.o buch1.o buch2.o buch3.o buch4.o galconj.o gen1.o gen2.o gen3.o ifactor1.o mpqs.o polarit1.o polarit2.o polarit3.o rootpol.o trans1.o trans2.o trans3.o elliptic.o galois.o kummer.o nffactor.o stark.o subfield.o subgroup.o thue.o anal.o compat.o errmsg.o helpmsg.o highlvl.o init.o sumiter.o
OBJSGP =  gp.o gp_init.o gp_rl.o plot.o plotport.o whatnow.o
OINLINE = mpinl.o

dft: gp.exe

lib: libpari.a

all: dft lib

libpari.a: $(OBJS) $(OINLINE)
	-$(RM) $@
	ar r $@ $(OINLINE) $(LDOBJS)
	ar s $@

gp: $(OBJS) $(OBJSGP) $(XTRAOINL)
	$(LD) -o gp $(LDFLAGS) $(XTRAOINL) $(LDOBJS) $(LDOBJSGP) $(RLLIBS) $(PLOTLIBS) $(LIBS)

gp.exe: gp
	$(EMXBIND)

cleantest:
	$(RM) *.out *.dif pari.ps

clean: cleantest
	-$(RM) *.o *.s libpari.* gp gp.exe

HEADERS = pariinl.h ../src/headers/pari.h ../src/headers/paritype.h ../src/headers/paricast.h ../src/headers/paricom.h ../src/headers/parierr.h ../src/headers/parigen.h ../src/headers/paristio.h ../src/headers/parisys.h ../src/headers/pariport.h
level0.s: ../src/kernel/ix86/l0asm.c ../src/kernel/ix86/l0asm.h
	$(CPP) -o $@ ../src/kernel/ix86/l0asm.c
kernel.o: level0.s
	$(CC) -c $(CFLAGS) $(CPPFLAGS) -o kernel.o level0.s
mp.o: $(HEADERS) ../src/kernel/none/mp.c
	$(CC) -c $(CFLAGS) $(CPPFLAGS) -o mp.o ../src/kernel/none/mp.c
mpinl.o: $(HEADERS) ../src/kernel/none/level1.h
	$(CC) -c $(CFLAGS) $(CPPFLAGS) -o mpinl.o ../src/kernel/none/level1.c

alglin1.o: $(HEADERS)  ../src/basemath/alglin1.c
	$(CC) -c $(CFLAGS) $(CPPFLAGS)  -o alglin1.o ../src/basemath/alglin1.c
alglin2.o: $(HEADERS)  ../src/basemath/alglin2.c
	$(CC) -c $(CFLAGS) $(CPPFLAGS)  -o alglin2.o ../src/basemath/alglin2.c
arith1.o: $(HEADERS)  ../src/basemath/arith1.c
	$(CC) -c $(CFLAGS) $(CPPFLAGS)  -o arith1.o ../src/basemath/arith1.c
arith2.o: $(HEADERS)  ../src/basemath/arith2.c
	$(CC) -c $(CFLAGS) $(CPPFLAGS)  -o arith2.o ../src/basemath/arith2.c
base1.o: $(HEADERS) ../src/headers/parinf.h ../src/basemath/base1.c
	$(CC) -c $(CFLAGS) $(CPPFLAGS)  -o base1.o ../src/basemath/base1.c
base2.o: $(HEADERS)  ../src/basemath/base2.c
	$(CC) -c $(CFLAGS) $(CPPFLAGS)  -o base2.o ../src/basemath/base2.c
base3.o: $(HEADERS)  ../src/basemath/base3.c
	$(CC) -c $(CFLAGS) $(CPPFLAGS)  -o base3.o ../src/basemath/base3.c
base4.o: $(HEADERS)  ../src/basemath/base4.c
	$(CC) -c $(CFLAGS) $(CPPFLAGS)  -o base4.o ../src/basemath/base4.c
base5.o: $(HEADERS)  ../src/basemath/base5.c
	$(CC) -c $(CFLAGS) $(CPPFLAGS)  -o base5.o ../src/basemath/base5.c
bibli1.o: $(HEADERS) ../src/headers/parinf.h ../src/basemath/bibli1.c
	$(CC) -c $(CFLAGS) $(CPPFLAGS)  -o bibli1.o ../src/basemath/bibli1.c
bibli2.o: $(HEADERS)  ../src/basemath/bibli2.c
	$(CC) -c $(CFLAGS) $(CPPFLAGS)  -o bibli2.o ../src/basemath/bibli2.c
buch1.o: $(HEADERS)  ../src/basemath/buch1.c
	$(CC) -c $(CFLAGS) $(CPPFLAGS)  -o buch1.o ../src/basemath/buch1.c
buch2.o: $(HEADERS) ../src/headers/parinf.h ../src/basemath/buch2.c
	$(CC) -c $(CFLAGS) $(CPPFLAGS)  -o buch2.o ../src/basemath/buch2.c
buch3.o: $(HEADERS)  ../src/basemath/buch3.c
	$(CC) -c $(CFLAGS) $(CPPFLAGS)  -o buch3.o ../src/basemath/buch3.c
buch4.o: $(HEADERS)  ../src/basemath/buch4.c
	$(CC) -c $(CFLAGS) $(CPPFLAGS)  -o buch4.o ../src/basemath/buch4.c
galconj.o: $(HEADERS)  ../src/basemath/galconj.c
	$(CC) -c $(CFLAGS) $(CPPFLAGS)  -o galconj.o ../src/basemath/galconj.c
gen1.o: $(HEADERS)  ../src/basemath/gen1.c
	$(CC) -c $(CFLAGS) $(CPPFLAGS)  -o gen1.o ../src/basemath/gen1.c
gen2.o: $(HEADERS)  ../src/basemath/gen2.c
	$(CC) -c $(CFLAGS) $(CPPFLAGS)  -o gen2.o ../src/basemath/gen2.c
gen3.o: $(HEADERS)  ../src/basemath/gen3.c
	$(CC) -c $(CFLAGS) $(CPPFLAGS)  -o gen3.o ../src/basemath/gen3.c
ifactor1.o: $(HEADERS)  ../src/basemath/ifactor1.c
	$(CC) -c $(CFLAGS) $(CPPFLAGS)  -o ifactor1.o ../src/basemath/ifactor1.c
mpqs.o: $(HEADERS)  ../src/modules/mpqs.c
	$(CC) -c $(CFLAGS) $(CPPFLAGS)  -o mpqs.o ../src/modules/mpqs.c
polarit1.o: $(HEADERS)  ../src/basemath/polarit1.c
	$(CC) -c $(CFLAGS) $(CPPFLAGS)  -o polarit1.o ../src/basemath/polarit1.c
polarit2.o: $(HEADERS)  ../src/basemath/polarit2.c
	$(CC) -c $(CFLAGS) $(CPPFLAGS)  -o polarit2.o ../src/basemath/polarit2.c
polarit3.o: $(HEADERS)  ../src/basemath/polarit3.c
	$(CC) -c $(CFLAGS) $(CPPFLAGS)  -o polarit3.o ../src/basemath/polarit3.c
rootpol.o: $(HEADERS)  ../src/basemath/rootpol.c
	$(CC) -c $(CFLAGS) $(CPPFLAGS)  -o rootpol.o ../src/basemath/rootpol.c
subgroup.o: $(HEADERS)  ../src/basemath/subgroup.c
	$(CC) -c $(CFLAGS) $(CPPFLAGS)  -o subgroup.o ../src/basemath/subgroup.c
trans1.o: $(HEADERS)  ../src/basemath/trans1.c
	$(CC) -c $(CFLAGS) $(CPPFLAGS)  -o trans1.o ../src/basemath/trans1.c
trans2.o: $(HEADERS)  ../src/basemath/trans2.c
	$(CC) -c $(CFLAGS) $(CPPFLAGS)  -o trans2.o ../src/basemath/trans2.c
trans3.o: $(HEADERS)  ../src/basemath/trans3.c
	$(CC) -c $(CFLAGS) $(CPPFLAGS)  -o trans3.o ../src/basemath/trans3.c
elliptic.o: $(HEADERS)  ../src/modules/elliptic.c
	$(CC) -c $(CFLAGS) $(CPPFLAGS)  -o elliptic.o ../src/modules/elliptic.c
galois.o: $(HEADERS)  ../src/modules/galois.c
	$(CC) -c $(CFLAGS) $(CPPFLAGS)  -o galois.o ../src/modules/galois.c
kummer.o: $(HEADERS)  ../src/modules/kummer.c
	$(CC) -c $(CFLAGS) $(CPPFLAGS)  -o kummer.o ../src/modules/kummer.c
nffactor.o: $(HEADERS)  ../src/modules/nffactor.c
	$(CC) -c $(CFLAGS) $(CPPFLAGS)  -o nffactor.o ../src/modules/nffactor.c
stark.o: $(HEADERS)  ../src/modules/stark.c
	$(CC) -c $(CFLAGS) $(CPPFLAGS)  -o stark.o ../src/modules/stark.c
subfield.o: $(HEADERS)  ../src/modules/subfield.c
	$(CC) -c $(CFLAGS) $(CPPFLAGS)  -o subfield.o ../src/modules/subfield.c
thue.o: $(HEADERS)  ../src/modules/thue.c
	$(CC) -c $(CFLAGS) $(CPPFLAGS)  -o thue.o ../src/modules/thue.c
anal.o: $(HEADERS) ../src/language/anal.h ../src/headers/parinf.h ../src/language/anal.c
	$(CC) -c $(CFLAGS) $(CPPFLAGS)  -o anal.o ../src/language/anal.c
compat.o: $(HEADERS)  ../src/language/compat.c
	$(CC) -c $(CFLAGS) $(CPPFLAGS)  -o compat.o ../src/language/compat.c
errmsg.o: $(HEADERS)  ../src/language/errmsg.c
	$(CC) -c $(CFLAGS) $(CPPFLAGS)  -o errmsg.o ../src/language/errmsg.c
es.o: $(HEADERS) ../src/language/anal.h ../src/language/es.c
	$(CC) -c $(CFLAGS) $(CPPFLAGS)  -o es.o ../src/language/es.c
helpmsg.o: $(HEADERS)  ../src/language/helpmsg.c
	$(CC) -c $(CFLAGS) $(CPPFLAGS)  -o helpmsg.o ../src/language/helpmsg.c
highlvl.o: $(HEADERS)  ../src/gp/highlvl.c
	$(CC) -c $(CFLAGS) $(CPPFLAGS)  -o highlvl.o ../src/gp/highlvl.c
init.o: $(HEADERS) ../src/language/anal.h ../src/language/init.c
	$(CC) -c $(CFLAGS) $(CPPFLAGS)  -o init.o ../src/language/init.c
sumiter.o: $(HEADERS) ../src/language/anal.h ../src/language/sumiter.c
	$(CC) -c $(CFLAGS) $(CPPFLAGS)  -o sumiter.o ../src/language/sumiter.c
gp.o: $(HEADERS) ../src/language/anal.h ../src/gp/gp.h ./paricfg.h ../src/gp/gp.c
	$(CC) -c $(CFLAGS) $(CPPFLAGS) -I../src/language $(RLINCLUDE) -o gp.o ../src/gp/gp.c
gp_init.o: $(HEADERS) ../src/graph/rect.h ../src/gp/gp_init.c
	$(CC) -c $(CFLAGS) $(CPPFLAGS) -I../src/graph -o gp_init.o ../src/gp/gp_init.c
gp_rl.o: $(HEADERS) ../src/language/anal.h ../src/gp/gp.h ./paricfg.h ../src/gp/gp_rl.c
	$(CC) -c $(CFLAGS) $(CPPFLAGS) -I../src/language $(RLINCLUDE) -o gp_rl.o ../src/gp/gp_rl.c
plot.o: $(HEADERS) ../src/graph/rect.h ../src/graph/$(PLOTFILE)
	$(CC) -c $(CFLAGS) $(CPPFLAGS) $(PLOTCFLAGS) -o plot.o ../src/graph/$(PLOTFILE)
plotport.o: $(HEADERS) ../src/graph/rect.h ../src/graph/plotport.c
	$(CC) -c $(CFLAGS) $(CPPFLAGS) -I../src/graph -o plotport.o ../src/graph/plotport.c
whatnow.o: $(HEADERS) ../src/language/anal.h ./paricfg.h ../src/gp/whatnow.c
	$(CC) -c $(CFLAGS) $(CPPFLAGS) -I../src/language -o whatnow.o ../src/gp/whatnow.c
