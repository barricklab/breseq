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

BUILT_SOURCES = version.h
DISTCLEANFILES =  version.h

PACKAGE_VERSION  = 1.3.1

version.h:
	echo '#define HTS_VERSION "$(PACKAGE_VERSION)"' > $@

AM_CFLAGS = -Wall -O2 $(SAMTOOLS_ADDITIONAL_CFLAGS)

noinst_LTLIBRARIES = libhts.la
libhts_la_SOURCES = \
	kfunc.c \
	knetfile.c \
	kstring.c \
	bgzf.c \
	faidx.c \
	hfile.c \
	hfile_net.c \
	hts.c \
	md5.c \
	regidx.c \
	sam.c \
	synced_bcf_reader.c \
	vcf_sweep.c \
	tbx.c \
	vcf.c \
	vcfutils.c \
	cram/cram_codecs.c \
	cram/cram_decode.c \
	cram/cram_encode.c \
	cram/cram_external.c \
	cram/cram_index.c \
	cram/cram_io.c \
	cram/cram_samtools.c \
	cram/cram_stats.c \
	cram/files.c \
	cram/mFILE.c \
	cram/open_trace_file.c \
	cram/pooled_alloc.c \
	cram/rANS_static.c \
	cram/sam_header.c \
	cram/string_alloc.c \
	cram/thread_pool.c \
	cram/vlen.c \
	cram/zfio.c


EXTRA_DIST = \
	INSTALL \
	LICENSE \
	Makefile.original \
	Makefile.xcode \
	Makefile.am \
	NEWS \
	README \
	bgzip.c \
	config.h \
	config.h.in \
	config.mk.in \
	configure \
	configure.ac \
	cram/cram.h \
	cram/cram_codecs.h \
	cram/cram_decode.h \
	cram/cram_encode.h \
	cram/cram_index.h \
	cram/cram_io.h \
	cram/cram_samtools.h \
	cram/cram_stats.h \
	cram/cram_structs.h \
	cram/mFILE.h \
	cram/misc.h \
	cram/open_trace_file.h \
	cram/os.h \
	cram/pooled_alloc.h \
	cram/rANS_byte.h \
	cram/rANS_static.h \
	cram/sam_header.h \
	cram/string_alloc.h \
	cram/thread_pool.h \
	cram/vlen.h \
	cram/zfio.h \
	faidx.5 \
	hfile_internal.h \
	hfile_irods.c \
	hfile_libcurl.c \
	hts_internal.h \
	htsfile.1 \
	htslib.mk \
	htslib.pc.in \
	htslib/bgzf.h \
	htslib/cram.h \
	htslib/faidx.h \
	htslib/hfile.h \
	htslib/hts.h \
	htslib/hts_defs.h \
	htslib/kbitset.h \
	htslib/kfunc.h \
	htslib/khash.h \
	htslib/khash_str2int.h \
	htslib/klist.h \
	htslib/knetfile.h \
	htslib/kseq.h \
	htslib/ksort.h \
	htslib/kstring.h \
	htslib/regidx.h \
	htslib/sam.h \
	htslib/synced_bcf_reader.h \
	htslib/tbx.h \
	htslib/vcf.h \
	htslib/vcf_sweep.h \
	htslib/vcfutils.h \
	htslib_vars.mk \
	plugin.c \
	sam.5 \
	synced_bcf_reader.c \
	tabix.1 \
	tabix.c \
	vcf.5 \
	vcf_sweep.c
