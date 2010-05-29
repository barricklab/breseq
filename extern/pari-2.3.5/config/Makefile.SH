file=$objdir/Makefile

echo Extracting $file
rm -f $file

# Path to directories
cfg=../$config_dir 
data=../$data_dir 
doc=../$doc_dir 
desc=../$desc_dir 
emx=../$emacs_dir 
ex=../$examples_dir 
misc=../$misc_dir 
src=../$src_dir 

knone=$src/kernel/none
kern0=$src/kernel/$kernlvl0
kern1=$src/kernel/$kernlvl1

#
# File lists
#
kernel="mp mpinl"
gp=`ls $src_dir/gp/*.c | sed 's,.*/\(.*\)\.c,\1,'`
language=`ls $src_dir/language/*.c | sed 's,.*/\(.*\)\.c,\1,'`
basemath=`ls $src_dir/basemath/*.c | sed 's,.*/\(.*\)\.c,\1,'`
modules=`ls $src_dir/modules/*.c   | sed 's,.*/\(.*\)\.c,\1,'`
systems=
if test -d "$src_dir/systems/$osname"; then
  systems=`ls $src_dir/systems/$osname/*.c | sed 's,.*/\(.*\)\.c,\1,'`
fi
# remove paridecl: no general recompilation when only changing a prototype
hlist=`ls $src_dir/headers/*.h | grep -v paridecl | grep -v paripriv |\
       sed 's,.*/\(.*\)\.h,\1,'`
# for installation put paridecl back. Remove \n otherwise make will choke.
headerlist="paridecl paripriv `echo $hlist| tr '\n' ' '`"

# special systems (OS/2 for now)
shortlib_prefix=lib
lib_prefix=lib
dllib_prefix=lib
case "$osname" in
  os2)  shortlib_prefix= # 8.3 filenames
        dllib_prefix= # != lib_prefix to allow gp-dyn link with DLL
	export_lib_create=emximp; export_lib='$(DLLIBPARI)$(_A)'
	export_file=pari.def; export_create="emxexp -u"
	# Actually, the build will fail until the switch to -Zomf
	dlld_ignore=- ;;
esac

PLOTCFLAGS=
PLOTLIBS=
postconfig="@: "
case "$which_graphic_lib" in
none)
  graph=plotnull;;
qt)
  PLOTCFLAGS='-D__FANCY_WIN__ -I$(QTDIR)/include'
  PLOTLIBS="-L\$(QTDIR)/lib $QTLIB"
  graph=plotQt;;
fltk)
  PLOTCFLAGS="-I\$(FLTKDIR)/include $X11_INC"
  PLOTLIBS="-L\$(FLTKDIR)/lib -lfltk $FLTK_LIBS"
  postconfig='-fltk-config --post '
  graph=plotfltk;;
esac
graph="plotport $graph"

plotrunpath=
case "$which_graphic_lib" in
  *X11*)
    PLOTCFLAGS="$PLOTCFLAGS $X11_INC"
    PLOTLIBS="$PLOTLIBS $X11_LIBS"
    plotrunpath=$X11
    graph="plotX $graph"
  ;;
esac

KERNOBJS=
for f in $kernel; do
  KERNOBJS="$KERNOBJS $f\$(_O)"
done
OBJS=$KERNOBJS
for f in $basemath $language $modules $systems; do
  OBJS="$OBJS $f\$(_O)"
done
OBJSGP=
for f in $gp $graph; do
  OBJSGP="$OBJSGP $f\$(_O)"
done
HEADERS=pariinl.h
for f in $hlist; do
  HEADERS="$HEADERS $src/headers/$f.h";
done

# runpath
tmp=$runpath
for d in "$plotrunpath" "$gmp" "$readline" "/usr/lib"; do
  if test -n "$d"; then
    case "$tmp" in
      $d|*$dir_sep$d|*$dir_sep$d$dir_sep*);;
      *) tmp="$tmp$dir_sep$d";;
    esac
  fi
done
RUNPTH_FINAL=
RUNPTH=
if test -n "$runpathprefix"; then
  RUNPTH_FINAL=`$config_dir/ldflags "$LDneedsWl" "$runpathprefix$tmp"`
  RUNPTH=`$config_dir/ldflags "$LDneedsWl" "$runpathprefix$TOP/$objdir$dir_sep$tmp"`
fi

if test -z "$DLLD"; then
  exec="gp-sta"; static=y
else
  exec="gp-sta gp-dyn"
fi
case "$sizeof_long" in
  4) numbits=32;;
  8) numbits=64;;
esac
dotest="env \"LD_LIBRARY_PATH=.$dir_sep\$\$LD_LIBRARY_PATH\" \$(SHELL) $src/test/dotest $numbits"
if test -z "$emacs"; then
  install_emacs_files=
else
  install_emacs_files=install_emacs_files
fi
case "$static" in
  n) dft=dyn; libdft=lib-dyn;;
  y) dft=sta; libdft= ;;
esac

CPPFLAGS="-I. -I$src/headers"
if test "$has_dlopen" = builtin; then
  CPPFLAGS="$CPPFLAGS -I$src/systems/$osname"
fi

case "$ln_s" in
  *cp*) ln_objdir=".";;
  *)    ln_objdir="$objdir"
esac

cat > $file << EOT
# Makefile for Pari/GP -- $pretty
#
#  This file was created by Configure. Any change made to it will be
#  lost when Configure is run.
#
SHELL      = $make_sh
PERL       = $perl
AR	   = ar

SIZEOF_LONG= $sizeof_long

CC_FLAVOR  =
CC         = $CC \$(CC_FLAVOR)
CPPFLAGS   = $CPPFLAGS
CFLAGS     = $CFLAGS
DLCFLAGS   = $DLCFLAGS
LD_FLAVOR  = \$(CC_FLAVOR)
LD         = $LD \$(LD_FLAVOR)
LDFLAGS    = $LDFLAGS
DLLD_FLAVOR  = \$(LD_FLAVOR)
DLLD       = $DLLD \$(DLLD_FLAVOR)
DLLDFLAGS  = $DLLDFLAGS
EXTRADLLDFLAGS = $EXTRADLLDFLAGS

# HIGHLY EXPERIMENTAL (only tested with gmp-4.0 on ix86 and Ultra).
# If you've configured and compiled GMP and would like to tune PARI using
# the nice cycle counting functions in GMP, uncomment the 4 lines below
# (correct the first one to the path to your gmp source tree).
#
#GMP     = /some/directory/gmp-4.0/tune
#GMPFLAGS= -DGMP_TIMER
#GMPO1 = \$(GMP)/time.o \$(GMP)/freq.o
#GMPOBJS=\$(GMPO1) \$(GMPO2)
#
#You may need to add a few object files to GMPOBJS. On UltraSparc, uncomment
#the following line
#GMPO2 = \$(GMP)/sparcv9.o

_O	   = .o
_A	   = .a
LIB_PREFIX = $lib_prefix
DLLIB_PREFIX = $dllib_prefix
LIBPARI_BASE = $libpari_base
LIBPARI      = \$(LIB_PREFIX)pari
DLLIBPARI    = \$(DLLIB_PREFIX)pari
LIBPARI_STA  = \$(LIBPARI)\$(_A)
LIBPARI_SO   = \$(DLLIBPARI).$DLSUFFIX
SOLIBPARI    = \$(DLLIB_PREFIX)\$(LIBPARI_BASE).$DLSUFFIX
LIBPARI_DYN   = \$(SOLIBPARI)$sodest
LIBPARI_SONAME= \$(SOLIBPARI)$soname
DL_DFLT_NAME = $DL_DFLT_NAME

EXPORT_FILE   = $export_file
EXPORT_CREATE = $export_create
EXPORT_LIB    = $export_lib
EXPORT_LIB_CREATE = $export_lib_create
DLLD_IGNORE   = $dlld_ignore

RUNPTH       = $RUNPTH
RUNPTH_FINAL = $RUNPTH_FINAL
LDDYN        = $LDDYN
LIBS         = $LIBS

RM = rm -f
MV = mv -f
LN = $ln_s
CP_F = cp -f
STRIP      = strip
STRIPFLAGS =

# Change these installation directories to suit your needs.
# DESTDIR is used to install to a false hierachy (to build a Debian package)
INCLUDEDIR= "\$(DESTDIR)$includedir"
LIBDIR    = "\$(DESTDIR)$libdir"
BINDIR    = "\$(DESTDIR)$bindir"
MANDIR    = "\$(DESTDIR)$mandir"
EMACSDIR  = "\$(DESTDIR)$emacsdir"
DATADIR   = "\$(DESTDIR)$datadir"
SYSDATADIR= "\$(DESTDIR)$sysdatadir"

EXDIR     = \$(DATADIR)/examples
MISCDIR   = \$(DATADIR)/misc
DOCDIR    = \$(DATADIR)/doc

INSTALL = $cfg/install
INSTALL_PROGRAM = \$(INSTALL)
INSTALL_DATA = \$(INSTALL) -m 644

# Readline
RLINCLUDE = $RLINCLUDE
RLLIBS    = $RLLIBS
# GMP
GMPINCLUDE = $GMPINCLUDE
# Graphic library.
QTDIR      = "$QTDIR"
MOC        = \$(QTDIR)/bin/moc
PLOTCFLAGS = $PLOTCFLAGS
PLOTLIBS   = $PLOTLIBS
CPLUSPLUS  = g++

FLTKDIR    = "$FLTKDIR"

TOPLDDYN   = -L"$TOP/$objdir" \$(RLLIBS)
# Description system
DESC          = pari.desc
DESC_HELP_GEN = $desc/gen_proto $desc/PARI/822.pm
TMPSUF        = -$osname-$arch-$$.tmp

DOTEST=$dotest

OBJS   = $OBJS
OBJSGP = $OBJSGP

gp: gp-$dft ../gp$suffix$exe_suff $libdft

../gp$suffix$exe_suff: gp-$dft
	-\$(RM) \$@
	-\$(LN) $ln_objdir/gp-$dft$exe_suff \$@

all: $exec lib-sta

lib-sta: \$(LIBPARI_STA)
lib-dyn: \$(LIBPARI_DYN)

\$(LIBPARI_STA): $add_funclist \$(OBJS)
	-\$(RM) \$@
	\$(AR) r \$@ \$(OBJS)
	-$ranlib \$@

kerntest\$(_O): $src/test/kerntest.c
	\$(CC) -c \$(CFLAGS) -I$src/language \$(CPPFLAGS) -o \$@ $src/test/kerntest.c

kerntest: $KERNOBJS kerntest\$(_O)
	\$(CC) \$(CFLAGS) \$(CPPFLAGS) -o \$@ $KERNOBJS kerntest\$(_O) \$(LIBS)

mpinl\$(_O): .headers parilvl0.h parilvl1.h $knone/mpinl.c
	\$(CC) -c \$(CFLAGS) \$(CPPFLAGS) \$(DLCFLAGS) -o mpinl\$(_O) $knone/mpinl.c

test-kernel:: kerntest
	@./kerntest > gp.out;\
	diff -c gp.out $src/test/$numbits/kernel > kern.dif;\
	if test -s kern.dif; then echo "KERNEL BUG"; else echo OK; fi

tune.o: $src/test/tune.c
	\$(CC) \$(CFLAGS) \$(GMPFLAGS) \$(CPPFLAGS) -o \$@ \$< -c

tune: tune-sta
	-\$(RM) tune ../tune
	-\$(LN) tune-sta tune
	-\$(LN) $ln_objdir/tune ../tune

tune-dyn: tune.o \$(LIBPARI_DYN)
	\$(LD) \$(LDFLAGS) \$(RUNPTH) -o \$@ \$< \$(GMPOBJS) -L. \$(LDDYN) \$(LIBS)
tune-sta: tune.o \$(LIBPARI_STA)
	\$(LD) \$(LDFLAGS) -o \$@ \$< \$(GMPOBJS) ./\$(LIBPARI_STA) \$(LIBS)

gp-sta: $add_funclist \$(OBJS) \$(OBJSGP)
	\$(RM) \$@
	\$(LD) -o \$@ \$(LDFLAGS) \$(OBJS) \$(OBJSGP) \$(RUNPTH) \$(RLLIBS) \$(PLOTLIBS) \$(LIBS)
	$postconfig gp-sta

cleanemacs:
	\$(RM) $emx/*.elc
cleantest:
	\$(RM) *.dif gp.out pari.ps
cleanobj:
	-\$(RM) *\$(_O) mp.c *.s pariinl.h parilvl0.h parilvl1.h libpari* tune* $exec
cleandesc:
	-\$(RM) $desc/\$(DESC) *\$(TMPSUF)

clean: cleanobj cleantest cleanemacs cleandesc

.headers: $HEADERS
	@touch \$@

install: install-lib-$dft install-include install-bin install-man install-misc install-doc install-examples install-cfg install-data

install-include:
	-mkdir -p \$(INCLUDEDIR)/pari
	-for i in paricfg.h pariinl.h; do \\
	  \$(INSTALL_DATA) \$\$i \$(INCLUDEDIR)/pari; done
	-for i in $headerlist; do \\
	   \$(INSTALL_DATA) $src/headers/\$\$i.h  \$(INCLUDEDIR)/pari; done
	-\$(RM) \$(INCLUDEDIR)/pari/genpari.h
	-\$(LN) pari.h \$(INCLUDEDIR)/pari/genpari.h

install-bin: install-bin-$dft

install-bin-dyn: gp-dyn
	-mkdir -p \$(BINDIR)
	-\$(RM) \$(BINDIR)/gp-$version$exe_suff \$(BINDIR)/gp$exe_suff
	\$(LD) -o \$(BINDIR)/gp-$version$exe_suff \$(LDFLAGS) \$(OBJSGP) \$(RUNPTH_FINAL) \$(TOPLDDYN) \$(LDDYN) \$(PLOTLIBS) \$(LIBS)
	-\$(STRIP) \$(STRIPFLAGS) \$(BINDIR)/gp-$version$exe_suff
	$postconfig \$(BINDIR)/gp-$version$exe_suff
	-\$(LN) gp-$version$exe_suff \$(BINDIR)/gp$exe_suff

install-bin-sta: gp-sta
	-mkdir -p \$(BINDIR)
	-\$(RM) \$(BINDIR)/gp-$version$exe_suff \$(BINDIR)/gp$exe_suff
	\$(INSTALL_PROGRAM) gp-sta$exe_suff \$(BINDIR)/gp-$version$exe_suff
	-\$(LN) gp-$version$exe_suff \$(BINDIR)/gp$exe_suff
# Can't strip it if we want install() to work on OSF.
#	-\$(STRIP) \$(STRIPFLAGS)  \$(BINDIR)/gp-$version
	$postconfig \$(BINDIR)/gp-$version$exe_suff

install-man::
	-mkdir -p \$(MANDIR)
	-\$(RM) \$(MANDIR)/pari.1 \$(MANDIR)/gp.1 \$(MANDIR)/gp-$version.1
	\$(INSTALL_DATA) $doc/gphelp.1 \$(MANDIR)
	\$(INSTALL_DATA) $doc/gp.1 \$(MANDIR)/gp-$version.1
	-\$(LN) gp.1 \$(MANDIR)/pari.1
	-\$(LN) gp-$version.1 \$(MANDIR)/gp.1

install-misc: $install_emacs_files
	-mkdir -p \$(MISCDIR) \$(BINDIR)
	\$(INSTALL_PROGRAM) $misc/tex2mail \$(BINDIR)
	\$(INSTALL_DATA) $doc/tex2mail.1 \$(MANDIR)
	\$(INSTALL_DATA) $misc/README    \$(MISCDIR)
	\$(INSTALL_DATA) $misc/color.dft \$(MISCDIR)
	\$(INSTALL_DATA) $misc/gpalias   \$(MISCDIR)
	\$(INSTALL_PROGRAM) $misc/gpflog \$(MISCDIR)
	\$(INSTALL_DATA) $misc/gprc.dft  \$(MISCDIR)
	\$(INSTALL_DATA) $misc/new.dic   \$(MISCDIR)
	\$(INSTALL_DATA) $misc/pari.xpm  \$(MISCDIR)
	\$(INSTALL_PROGRAM) $misc/xgp    \$(MISCDIR)

install-cfg::
	-mkdir -p \$(SYSDATADIR)
	-\$(INSTALL_DATA) pari.cfg \$(SYSDATADIR)
	-if test -n "$add_funclist"; then\
	   mkdir -p \$(DATADIR)/PARI;\
	   \$(INSTALL_DATA) $desc/PARI/822.pm \$(DATADIR)/PARI;\
	   \$(INSTALL_DATA) $desc/\$(DESC) \$(DATADIR); fi

install-doc::
	-mkdir -p \$(BINDIR) \$(DOCDIR)
	-cd $doc; \$(MAKE) all
	-\$(INSTALL_PROGRAM) $doc/gphelp    \$(BINDIR)
	-\$(INSTALL_DATA) $doc/Makefile     \$(DOCDIR)
	-\$(INSTALL_DATA) $doc/translations \$(DOCDIR)
	-\$(INSTALL_DATA) $doc/appa.tex     \$(DOCDIR)
	-\$(INSTALL_DATA) $doc/appb.tex     \$(DOCDIR)
	-\$(INSTALL_DATA) $doc/appc.tex     \$(DOCDIR)
	-\$(INSTALL_DATA) $doc/parimacro.tex \$(DOCDIR)
	-\$(INSTALL_DATA) $doc/pdfmacs.tex  \$(DOCDIR)
	-\$(INSTALL_DATA) $doc/refcard.tex  \$(DOCDIR)
	-\$(INSTALL_DATA) $doc/tutorial.tex \$(DOCDIR)
	-\$(INSTALL_DATA) $doc/users.tex    \$(DOCDIR)
	-\$(INSTALL_DATA) $doc/usersch1.tex \$(DOCDIR)
	-\$(INSTALL_DATA) $doc/usersch2.tex \$(DOCDIR)
	-\$(INSTALL_DATA) $doc/usersch3.tex \$(DOCDIR)
	-\$(INSTALL_DATA) $doc/usersch4.tex \$(DOCDIR)
	-\$(INSTALL_DATA) $doc/usersch5.tex \$(DOCDIR)
	-\$(INSTALL_DATA) $doc/paricfg.tex  \$(DOCDIR)
	-\$(INSTALL_DATA) $doc/libpari.dvi    \$(DOCDIR)
	-\$(INSTALL_DATA) $doc/users.dvi    \$(DOCDIR)
	-\$(INSTALL_DATA) $doc/tutorial.dvi \$(DOCDIR)
	-\$(INSTALL_DATA) $doc/refcard.dvi  \$(DOCDIR)
	-\$(INSTALL_DATA) $doc/refcard.ps   \$(DOCDIR)

install-examples:
	-mkdir -p \$(EXDIR)
	-\$(INSTALL_DATA) $ex/EXPLAIN     \$(EXDIR)
	-\$(INSTALL_DATA) $ex/Inputrc     \$(EXDIR)
	-\$(INSTALL_DATA) $ex/Makefile    \$(EXDIR)
	-\$(INSTALL_DATA) $ex/bench.gp    \$(EXDIR)
	-\$(INSTALL_DATA) $ex/cl.gp       \$(EXDIR)
	-\$(INSTALL_DATA) $ex/classno.gp  \$(EXDIR)
	-\$(INSTALL_DATA) $ex/contfrac.gp \$(EXDIR)
	-\$(INSTALL_DATA) $ex/lucas.gp    \$(EXDIR)
	-\$(INSTALL_DATA) $ex/extgcd.c    \$(EXDIR)
	-\$(INSTALL_DATA) $ex/rho.gp      \$(EXDIR)
	-\$(INSTALL_DATA) $ex/squfof.gp   \$(EXDIR)
	-\$(INSTALL_DATA) $ex/taylor.gp   \$(EXDIR)

install-data:
	-if test -d $data; then cd $data; \
	   for d in \`ls\`; do \
	     mkdir -p \$(DATADIR)/\$\$d && \
	     for f in \`ls \$\$d\`; do \
	       \$(INSTALL_DATA) \$\$d/\$\$f \$(DATADIR)/\$\$d; \
	     done >/dev/null;\
	   done; \
	 fi

install-lib-sta: \$(LIBPARI_STA)
	-mkdir -p \$(LIBDIR)
	-\$(RM) \$(LIBDIR)/\$(LIBPARI_STA)
	\$(INSTALL_DATA) \$(LIBPARI_STA) \$(LIBDIR)/\$(LIBPARI_STA)

install-lib-dyn-base:
	-mkdir -p \$(LIBDIR)
	-\$(RM) \$(LIBDIR)/\$(LIBPARI_DYN) \$(LIBDIR)/\$(LIBPARI_SONAME) \$(LIBDIR)/\$(LIBPARI_SO)
	\$(INSTALL_PROGRAM) \$(LIBPARI_DYN) \$(LIBDIR)/\$(LIBPARI_DYN)
install-lib-dyn-link: install-lib-dyn-base
	-if test "\$(LIBPARI_SONAME)" != "\$(LIBPARI_SO)"; then \
	  \$(LN) \$(LIBPARI_DYN) \$(LIBDIR)/\$(LIBPARI_SONAME);\
	  \$(LN) \$(LIBPARI_DYN) \$(LIBDIR)/\$(LIBPARI_SO); fi
install-lib-dyn: \$(LIBPARI_DYN) install-lib-dyn-link

pariinl.h: parilvl0.h parilvl1.h
	cat parilvl0.h parilvl1.h > pariinl.h

bench: $exec
	@\$(DOTEST) $test_basic
dobench::
	@\$(DOTEST) $test_basic
test-all: $exec
	@\$(DOTEST) $test_basic $test_extra
dotest-all::
	@\$(DOTEST) $test_basic $test_extra
EOT
for i in $test_extra $test_extra_out; do
  cat >> $file << EOT
test-$i: $exec
	@\$(DOTEST) $i
dotest-$i::
	@\$(DOTEST) $i
EOT
done

if test -z "$DLLD"; then
  cat >> $file << EOT

\$(LIBPARI_DYN)::
	@echo "Configure could not find a way to build a shared library on this machine"

EOT
else
  if test -n "$export_file"; then
    cat >> $file << EOT

EXPORT_FILE_BASE = $src/systems/os2/pari.def.base
VERSION_VERBOSE = $pari_release_verbose

\$(EXPORT_FILE): \$(OBJS) \$(EXPORT_FILE_BASE)
	cat \$(EXPORT_FILE_BASE) | sed 's/<DLL_BASE>/\$(DLLIBPARI)/' | sed 's/<VENDOR>/pari.math.u-bordeaux.fr\//' | sed 's/<VERSION>/\$(VERSION_VERBOSE)/' | sed 's/<DESCR>/GP\/PARI compiled with \$(CFLAGS)/' > \$@
	\$(EXPORT_CREATE) \$(OBJS) >> \$@

\$(DLLIBPARI)\$(_A): \$(EXPORT_FILE)
	\$(EXPORT_LIB_CREATE) -o \$@ \$(EXPORT_FILE)

EOT
  fi
  cat >> $file << EOT

gp-dyn: $add_funclist \$(OBJSGP) \$(LIBPARI_DYN) \$(EXPORT_LIB)
	\$(RM) \$@
	\$(LD) -o \$@ \$(LDFLAGS) \$(OBJSGP) \$(RUNPTH) \$(TOPLDDYN) \$(LDDYN) \$(PLOTLIBS) \$(LIBS)
	$postconfig gp-dyn

\$(LIBPARI_DYN): $add_funclist \$(OBJS) \$(EXPORT_FILE)
	-\$(RM) \$(LIBPARI_DYN)
	\$(DLLD_IGNORE)\$(DLLD) -o \$(LIBPARI_DYN) \$(DLLDFLAGS) \$(OBJS) \$(EXTRADLLDFLAGS) \$(EXPORT_FILE)
	-if test "\$(LIBPARI_DYN)" != "\$(LIBPARI_SO)"; then \
	  \$(RM) \$(LIBPARI_SO);\
	  \$(LN) \$(LIBPARI_DYN) \$(LIBPARI_SO); fi
	-if test "\$(LIBPARI_DYN)" != "\$(LIBPARI_SONAME)"; then \
	  \$(RM) \$(LIBPARI_SONAME);\
	  \$(LN) \$(LIBPARI_DYN) \$(LIBPARI_SONAME); fi
EOT
fi

if test -n "$emacs"; then
  cat >> $file << EOT
elc:
	@echo "Byte-Compiling elisp files..."
	-cd $emx; $emacs -batch \
	   --eval '(setq load-path (append (list ".") load-path))' \
	   -f batch-byte-compile *.el

$install_emacs_files: elc
	-mkdir -p \$(EMACSDIR)
	-cd $emx; for f in *.elc *.el pariemacs.txt; do \
	  \$(INSTALL_DATA) \$\$f \$(EMACSDIR); \
	done
EOT
fi

if test -n "$add_funclist"; then
  cat >> $file << EOT
desc:: $src/language/init.h $src/gp/gp_init.h $src/gp/highlvl.h $src/language/members.h

$src/funclist::
	@-cksum \`env LANG= ls -d $src/functions/*/* | grep -v CVS\` > funclist\$(TMPSUF)
	@-sh $cfg/mv-if-diff funclist\$(TMPSUF) $src/funclist

$desc/\$(DESC): $src/funclist $desc/merge_822 $desc/PARI/822.pm
	cd $desc && \$(PERL) merge_822 ../functions/*/* > def\$(TMPSUF)
	mv $desc/def\$(TMPSUF) \$@

$src/language/members.h: $desc/\$(DESC) $desc/gen_member
	cd $desc && \$(PERL) gen_member \$(DESC) > members\$(TMPSUF)
	mv $desc/members\$(TMPSUF) \$@

$src/language/init.h: $desc/\$(DESC) \$(DESC_HELP_GEN)
	cd $desc && \$(PERL) gen_proto basic \$(DESC) > init\$(TMPSUF)
	mv $desc/init\$(TMPSUF) \$@

$src/gp/gp_init.h: $desc/\$(DESC) \$(DESC_HELP_GEN)
	cd $desc && \$(PERL) gen_proto gp \$(DESC) > gp_init\$(TMPSUF)
	mv $desc/gp_init\$(TMPSUF) \$@

$src/gp/highlvl.h: $desc/\$(DESC) \$(DESC_HELP_GEN)
	cd $desc && \$(PERL) gen_proto highlevel \$(DESC) > highlvl\$(TMPSUF)
	mv $desc/highlvl\$(TMPSUF) \$@

L0MODS=$kern0/asm0.h $knone/addll.h $knone/mulll.h $knone/bfffo.h $knone/divll.h
EOT
fi

# Level 0
f=$src_dir/kernel/$kernlvl0/MakeLVL0.SH
if test -s $f; then
  . $f
else
  cat >> $file <<EOT
parilvl0.h: \$(L0MODS)
	$cfg/genkernel $kern0/asm0.h > parilvl0.h
EOT
fi

# Level 1
f=$src_dir/kernel/$kernlvl1/MakeLVL1.SH
if test -s $f; then . $f; fi

HUGELINE=
for dir in basemath modules language gp graph graphdyn systems; do
  eval list='$'$dir
  case "$dir" in
    systems)  SRC=$src/$dir/$osname;;
    *)        SRC=$src/$dir;;
  esac
  for f in $list; do

  source="$SRC/$f.c"
  f=`basename $f`
  depend=
  cflags="\$(CFLAGS) \$(CPPFLAGS)"
  compile='$(CC)'
  case "$f" in
  gp|gp_rl)
    cflags="$cflags -I$src/language \$(RLINCLUDE)"
    depend="$src/language/anal.h $src/gp/gp.h ./paricfg.h"
    ;;
  gp_init)
    cflags="$cflags -I$src/graph"
    depend="$src/graph/rect.h $src/gp/gp_init.h"
    ;;
  default)
    depend="$src/language/anal.h ./paricfg.h"
    cflags="$cflags \$(DLCFLAGS)"
    ;;
  plotport)
    cflags="$cflags -I$src/graph"
    depend="$src/graph/rect.h"
    ;;
  highlvl)
    cflags="$cflags -DDL_DFLT_NAME=\$(DL_DFLT_NAME)"
    depend="$src/gp/highlvl.h"
    ;;
  plotQt)
    cflags="$cflags \$(PLOTCFLAGS)"
    depend="$src/graph/rect.h"
    compile="\$(MOC) -o plotQt.moc.cpp $src/graph/plotQt.c && \$(CPLUSPLUS)"
    ;;
  plotfltk)
    cflags="$cflags \$(PLOTCFLAGS)"
    depend="$src/graph/rect.h"
    compile="\$(CPLUSPLUS)"
    ;;
  plot*)
    cflags="$cflags \$(PLOTCFLAGS)"
    depend="$src/graph/rect.h"
    ;;
  es|sumiter)
    depend="$src/language/anal.h"
    cflags="$cflags \$(DLCFLAGS)"
    ;;
  whatnow);;
  init)
    depend="$src/language/anal.h $src/language/init.h"
    cflags="$cflags \$(DLCFLAGS)"
    ;;
  anal)
    depend="$src/language/anal.h $src/language/members.h"
    cflags="$cflags \$(DLCFLAGS)"
    ;;
  *)
    cflags="$cflags \$(DLCFLAGS)"
    ;;
  esac
HUGELINE="$HUGELINE
$f\$(_O): .headers $depend $source
	$compile -c $cflags -o $f\$(_O) $source"
  done
done
# avoid execing too many "cat". MS-DOS can't handle it
echo "$HUGELINE" >> $file

echo ".PHONY: gp dft clean"  >> $file
