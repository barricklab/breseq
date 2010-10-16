## Remember to drag along the LICENSE, NEWS, unused headers, etc.
EXTRA_DIST = \
	AUTHORS \
	ChangeLog \
	COPYING \
	INSTALL \
	Makefile.orig \
	Makefile.mingw \
	NEWS \
	win32/xcurses.h \
	win32/zconf.h \
	win32/zlib.h

noinst_LTLIBRARIES = libbam.la
libbam_la_SOURCES = \
	bgzf.c \
	kstring.c \
	bam_aux.c \
	bam.c \
	bam_import.c \
	sam.c \
	bam_index.c \
	bam_pileup.c \
	bam_lpileup.c \
	bam_md.c \
	glf.c \
	razf.c \
	faidx.c \
	knetfile.c \
	bam_sort.c \
	sam_header.c \
	bam_endian.h \
	bam_maqcns.h \
	bgzf.h \
	faidx.h \
	glf.h \
	kaln.h \
	khash.h \
	klist.h \
	knetfile.h \
	kseq.h \
	ksort.h \
	kstring.h \
	razf.h \
	sam_header.h \
	sam.h 

bin_PROGRAMS = breseq_samtools
breseq_samtools_SOURCES = \
	bam_tview.c \
	bam_maqcns.c \
	bam_plcmd.c \
	sam_view.c	\
	bam_rmdup.c \
	bam_rmdupse.c \
	bam_mate.c \
	bam_stat.c \
	bam_color.c	\
	bamtk.c \
	kaln.c 
breseq_samtools_CFLAGS = -D_CURSES_LIB=1
breseq_samtools_LDADD = libbam.la
