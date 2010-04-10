CC=			gcc
CFLAGS=		-g -Wall -O2 #-m64 #-arch ppc
DFLAGS=		-D_FILE_OFFSET_BITS=64 -D_USE_KNETFILE -D_CURSES_LIB=1
LOBJS=		bgzf.o kstring.o bam_aux.o bam.o bam_import.o sam.o bam_index.o	\
			bam_pileup.o bam_lpileup.o bam_md.o glf.o razf.o faidx.o knetfile.o	\
			bam_sort.o sam_header.o
AOBJS=		bam_tview.o bam_maqcns.o bam_plcmd.o sam_view.o	\
			bam_rmdup.o bam_rmdupse.o bam_mate.o bam_stat.o bam_color.o	\
			bamtk.o kaln.o
PROG=		samtools
INCLUDES=
SUBDIRS=	. misc
LIBPATH=
LIBCURSES=	-lcurses # -lXCurses

.SUFFIXES:.c .o

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

all-recur lib-recur clean-recur cleanlocal-recur install-recur:
		@target=`echo $@ | sed s/-recur//`; \
		wdir=`pwd`; \
		list='$(SUBDIRS)'; for subdir in $$list; do \
			cd $$subdir; \
			$(MAKE) CC="$(CC)" DFLAGS="$(DFLAGS)" CFLAGS="$(CFLAGS)" \
				INCLUDES="$(INCLUDES)" LIBPATH="$(LIBPATH)" $$target || exit 1; \
			cd $$wdir; \
		done;

all:$(PROG)

lib:libbam.a

libbam.a:$(LOBJS)
		$(AR) -cru $@ $(LOBJS)

samtools:lib $(AOBJS)
		$(CC) $(CFLAGS) -o $@ $(AOBJS) -lm $(LIBPATH) $(LIBCURSES) -lz -L. -lbam

razip:razip.o razf.o knetfile.o
		$(CC) $(CFLAGS) -o $@ razf.o razip.o knetfile.o -lz

bgzip:bgzip.o bgzf.o
		$(CC) $(CFLAGS) -o $@ bgzf.o bgzip.o -lz

razip.o:razf.h
bam.o:bam.h razf.h bam_endian.h kstring.h sam_header.h
sam.o:sam.h bam.h
bam_import.o:bam.h kseq.h khash.h razf.h
bam_pileup.o:bam.h razf.h ksort.h
bam_plcmd.o:bam.h faidx.h bam_maqcns.h glf.h
bam_index.o:bam.h khash.h ksort.h razf.h bam_endian.h
bam_lpileup.o:bam.h ksort.h
bam_tview.o:bam.h faidx.h bam_maqcns.h
bam_maqcns.o:bam.h ksort.h bam_maqcns.h
bam_sort.o:bam.h ksort.h razf.h
bam_md.o:bam.h faidx.h
glf.o:glf.h
sam_header.o:sam_header.h khash.h

faidx.o:faidx.h razf.h khash.h
faidx_main.o:faidx.h razf.h

cleanlocal:
		rm -fr gmon.out *.o a.out *.dSYM razip $(PROG) *~ *.a

clean:cleanlocal-recur
