###
##
## AUTHORS
##
## Jeffrey E. Barrick <jeffrey.e.barrick@gmail.com>
## David B. Knoester
##
## LICENSE AND COPYRIGHT
##
## Copyright (c) 2010 Michigan State University
##
## breseq is free software; you can redistribute it and/or modify it under the terms the
## GNU General Public License as published by the Free Software Foundation; either
## version 1, or (at your option) any later version.
##
###

SUBDIRS = htslib-1.3.1

BUILT_SOURCES = version.h
DISTCLEANFILES =  version.h
CLEANFILES = libbam.a

PACKAGE_VERSION  = 1.3.1

version.h:
	echo '#define SAMTOOLS_VERSION "$(PACKAGE_VERSION)"' > $@

AM_CFLAGS = -Wall -O2 -I$(top_srcdir)/$(HTSLIB) $(SAMTOOLS_ADDITIONAL_CFLAGS)

LIBADD = $(top_builddir)/$(HTSLIB)/libhts.la

noinst_LTLIBRARIES = libbam.la
libbam_la_SOURCES = \
	bam_aux.c \
	bam.c \
	bam_import.c \
	sam.c \
	sam_header.c \
	bam_plbuf.c \
	bam_index.c \
	bam_plcmd.c \
	sam_view.c \
	bam_cat.c \
	bam_md.c \
	bam_reheader.c \
	bam_sort.c \
	bedidx.c \
	kprobaln.c \
	bam_rmdup.c \
	bam_rmdupse.c \
	bam_mate.c \
	bam_stat.c \
	bam_color.c \
	bamtk.c \
	bam2bcf.c \
	bam2bcf_indel.c \
	errmod.c \
	sample.c \
	sam_opts.c \
	cut_target.c \
	phase.c \
	bam2depth.c \
	padding.c \
	bedcov.c \
	bamshuf.c \
	faidx.c \
	dict.c \
	stats.c \
	stats_isize.c \
	bam_flags.c \
	bam_split.c \
	bam_tview.c \
	bam_tview_html.c \
	bam_lpileup.c \
	bam_quickcheck.c \
	bam_addrprg.c

EXTRA_DIST = \
	AUTHORS \
	ChangeLog.old \
	INSTALL \
	LICENSE \
	Makefile \
	Makefile.mingw \
	Makefile.original \
	NEWS \
	README \
	aclocal.m4 \
	bam.h \
	bam2bcf.h \
	bam_endian.h \
	bam_lpileup.h \
	bam_plbuf.h \
	bam_tview.h \
	bam_tview_curses.c \
	config.h \
	config.h.in \
	config.mk \
	config.mk.in \
	configure \
	configure.ac \
	errmod.h \
	examples/00README.txt \
	examples/ex1.fa \
	examples/ex1.sam.gz \
	examples/toy.fa \
	examples/toy.sam \
	install-sh \
	kprobaln.h \
	misc/HmmGlocal.java \
	misc/ace2sam.c \
	misc/blast2sam.pl \
	misc/bowtie2sam.pl \
	misc/export2sam.pl \
	misc/interpolate_sam.pl \
	misc/maq2sam.c \
	misc/md5fa.c \
	misc/md5sum-lite.c \
	misc/novo2sam.pl \
	misc/plot-bamstats \
	misc/psl2sam.pl \
	misc/r2plot.lua \
	misc/sam2vcf.pl \
	misc/samtools.pl \
	misc/seq_cache_populate.pl \
	misc/soap2sam.pl \
	misc/varfilter.py \
	misc/vcfutils.lua \
	misc/wgsim.1 \
	misc/wgsim.c \
	misc/wgsim_eval.pl \
	misc/zoom2sam.pl \
	sam.h \
	sam_header.h \
	sam_opts.h \
	sample.h \
	samtools.1 \
	samtools.h \
	stats_isize.h \
	win32/libcurses.a \
	win32/libz.a \
	win32/xcurses.h \
	win32/zconf.h \
	win32/zlib.h
