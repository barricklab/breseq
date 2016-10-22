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

version.h:
	echo '#define HTS_VERSION "$(PACKAGE_VERSION)"' > $@

PACKAGE_VERSION  = 1.3.1
AM_CFLAGS = -g -Wall -O2

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
	test/aux\#aux_java.cram \
	test/auxf\#values.sam \
	test/auxf.fa \
	test/auxf.fa.fai \
	test/c1\#bounds.sam \
	test/c1\#clip.sam \
	test/c1\#noseq.sam \
	test/c1\#pad1.sam \
	test/c1\#pad2.sam \
	test/c1\#pad3.sam \
	test/c1\#unknown.sam \
	test/c1.fa \
	test/c1.fa.fai \
	test/c2\#pad.sam \
	test/c2.fa \
	test/c2.fa.fai \
	test/ce\#1.sam \
	test/ce\#2.sam \
	test/ce\#5.sam \
	test/ce\#5b.sam \
	test/ce\#5b_java.cram \
	test/ce\#large_seq.sam \
	test/ce\#supp.sam \
	test/ce\#tag_depadded.sam \
	test/ce\#tag_padded.sam \
	test/ce\#unmap.sam \
	test/ce\#unmap1.sam \
	test/ce\#unmap2.sam \
	test/ce.fa \
	test/ce.fa.fai \
	test/compare_sam.pl \
	test/cross_validate.sh \
	test/faidx.fa \
	test/fieldarith.c \
	test/fieldarith.sam \
	test/formatcols.vcf \
	test/hfile.c \
	test/noroundtrip-out.vcf \
	test/noroundtrip.vcf \
	test/sam.c \
	test/test-regidx.c \
	test/test-vcf-api.c \
	test/test-vcf-api.out \
	test/test-vcf-hdr-in.vcf \
	test/test-vcf-hdr.out \
	test/test-vcf-sweep.c \
	test/test-vcf-sweep.out \
	test/test.pl \
	test/test_view.c \
	test/xx\#blank.sam \
	test/xx\#large_aux.sam \
	test/xx\#large_aux2.sam \
	test/xx\#large_aux_java.cram \
	test/xx\#minimal.sam \
	test/xx\#pair.sam \
	test/xx\#repeated.sam \
	test/xx\#rg.sam \
	test/xx\#tlen.sam \
	test/xx\#tlen2.sam \
	test/xx\#triplet.sam \
	test/xx\#unsorted.sam \
	test/xx.fa \
	test/xx.fa.fai \
	vcf.5 \
	vcf_sweep.c
