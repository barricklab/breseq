###
##
## MODIFICATIONS
##
## This automake Makefile.am is used to configure, compile, and install SAMtools for
##   use by the breseq computational pipeline.
## 
## The original Makefile distributed with SAMtools has been renamed Makefile.orig
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
## version 1, or (at your option) any later version. Or under any terms
##
###

SUBDIRS = bcftools

AM_LDFLAGS = -static 

## Remember to drag along the LICENSE, NEWS, unused headers, etc.
EXTRA_DIST = \
	../samtools_modifications.txt \
	AUTHORS \
	ChangeLog \
	COPYING \
	INSTALL \
	Makefile.original \
	Makefile.mingw \
	NEWS \
	samtools.1 \
	win32/xcurses.h \
	win32/zconf.h \
	win32/zlib.h \
	examples/00README.txt \
	examples/bam2bed.c \
	examples/calDepth.c \
	examples/ex1.fa \
	examples/ex1.sam.gz \
	examples/Makefile \
	examples/toy.fa \
	examples/toy.sam \
	misc/blast2sam.pl \
	misc/bowtie2sam.pl \
	misc/export2sam.pl \
	misc/HmmGlocal.java \
	misc/interpolate_sam.pl \
	misc/Makefile \
	misc/maq2sam.c \
	misc/md5.c \
	misc/md5.h \
	misc/md5fa.c \
	misc/novo2sam.pl \
	misc/psl2sam.pl \
	misc/sam2vcf.pl \
	misc/samtools.pl \
	misc/seqtk.c \
	misc/soap2sam.pl \
	misc/varfilter.py \
	misc/wgsim_eval.pl \
	misc/wgsim.c \
	misc/zoom2sam.pl

noinst_LTLIBRARIES = libbam.la

libbam_la_SOURCES = \
	bam_aux.c \
	bam_cat.c \
	bam_endian.h \
	bam_import.c \
	bam_index.c \
	bam_lpileup.c \
	bam_md.c \
	bam_pileup.c \
	bam_reheader.c \
	bam_sort.c \
	bam.c \
	bam.h \
	bam2bcf.h \
	bam2bcf_indel.c \
	bam2depth.c \
	bedidx.c \
	bgzf.c \
	bgzf.h \
	bgzip.c \
	faidx.c \
	faidx.h \
	kaln.h \
	khash.h \
	klist.h \
	knetfile.c \
	knetfile.h \
	kprobaln.c \
	kprobaln.h \
	kseq.h \
	ksort.h \
	kstring.c \
	kstring.h \
	razf.c \
	razf.h \
	razip.c \
	sam_header.c \
	sam_header.h \
	sam.c \
	sam.h
libbam_la_LDFLAGS = -fPIC
libbam_la_CFLAGS = -std=gnu89 -D_CURSES_LIB=1 -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -D_USE_KNETFILE $(SAMTOOLS_ADDITIONAL_CFLAGS)

bin_PROGRAMS = samtools
samtools_SOURCES = \
	bam_tview.c \
	bam_plcmd.c \
	sam_view.c	\
	bam_rmdup.c \
	bam_rmdupse.c \
	bam_mate.c \
	bam_stat.c \
	bam_color.c	\
	bamtk.c \
	cut_target.c \
	kaln.c \
	bam2bcf.c \
	bam2bcf_indel.c \
	errmod.c \
	errmod.h \
	sample.c \
	sample.h \
	phase.c 
samtools_CFLAGS = -std=gnu89 -D_CURSES_LIB=1 -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -D_USE_KNETFILE $(SAMTOOLS_ADDITIONAL_CFLAGS)
samtools_LDADD = libbam.la bcftools/libbcf.la

noinst_HEADERS = bam.h


