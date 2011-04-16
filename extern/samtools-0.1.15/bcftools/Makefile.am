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

## Remember to drag along the LICENSE, NEWS, unused headers, etc.
noinst_LTLIBRARIES = libbcf.la



##	call1.h \
##	main.h \ 

libbcf_la_SOURCES = \
	bcf.c \
	vcf.c \
	bcfutils.c \
	prob1.c \
	ld.c \
	kfunc.c \
	index.c \
	fet.c \
	bcf2qcall.c \
	bcf.h \
	vcf.h \
	bcfutils.h \
	prob1.h \
	ld.h \
	kfunc.h \
	index.h \
	fet.h \
	bcf2qcall.h \
	kstring.h \
	bgzf.h \
	knetfile.h \
	bedidx.h
libbcf_la_LDFLAGS = -fPIC
libbcf_la_CPPFLAGS = -I$(top_srcdir)/$(LIBBAM)


