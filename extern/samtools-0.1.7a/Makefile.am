###
##
## MODIFICATIONS
##
## This automake Makefile.ac is used to configure, compile, and install SAMtools for
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

lib_LTLIBRARIES = libbam.la
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
libbam_la_LDFLAGS = -static

bin_PROGRAMS = samtools
samtools_SOURCES = \
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
samtools_CFLAGS = -D_CURSES_LIB=1
samtools_LDADD = libbam.la

include_HEADERS = bam.h


