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

AM_LDFLAGS = -static 

EXTRA_DIST = \
	Makefile.original \
	README \
	vcfutils.pl

noinst_LTLIBRARIES = libbcf.la

libbcf_la_SOURCES = \
	bcf.c \
	bcf.h \
	bcf.tex \
	bcf2qcall.c \
	vcf.c \
	bcfutils.c \
	prob1.c \
	prob1.h \
	em.c \
	kfunc.c \
	kmin.c \
	kmin.h \
	index.c \
	fet.c \
	mut.c \
	bcf2qcall.c
libbcf_la_LDFLAGS = -fPIC
libbcf_la_CFLAGS = -I$(top_srcdir)/$(LIBBAM) $(SAMTOOLS_ADDITIONAL_CFLAGS)

bin_PROGRAMS = bcftools
bcftools_SOURCES = \
	call1.c \
	main.c \
	../bgzf.c \
	../kstring.c \
	../knetfile.c \
	../bedidx.c
bcftools_CFLAGS = -g -Wall -O2 -D_FILE_OFFSET_BITS=64 -D_USE_KNETFILE -I$(top_srcdir)/$(LIBBAM) $(SAMTOOLS_ADDITIONAL_CFLAGS)
bcftools_LDADD = libbcf.la
