###
##
## AUTHORS
##
##   Jeffrey E. Barrick <jeffrey.e.barrick@gmail.com> and other contributors
##
## LICENSE AND COPYRIGHT
##
##   Copyright (c) 2008-2010 Michigan State University
##   Copyright (c) 2011-2025 The University of Texas at Austin
##   Copyright (c) 2025-     Michigan State University
##
##   breseq is free software; you can redistribute it and/or modify it under the
##   terms of the GNU General Public License as published by the Free Software
##   Foundation; either version 2, or (at your option) any later version.
##
##   SPDX-License-Identifier: GPL-2.0-or-later
##
###

##This is the most compatible way of passing through the data directory to the C++ code
BRESEQDATADIR= $(DESTDIR)$(datadir)/$(PACKAGE_NAME)
BRESEQRELATIVEDATADIR= $(BIN_TO_DATA_RELATIVE_PATH)/$(PACKAGE_NAME)

# $(STATIC_LIBZ_INC) and $(MINIZ_INC) add include paths for static builds.
# -DMINIZ_NO_ZLIB_COMPATIBLE_NAMES prevents miniz aliases clashing with libz.
AM_CPPFLAGS = -O3 $(STATIC_LIBZ_INC) $(MINIZ_INC) $(HTSLIB_INC) -DMINIZ_NO_ZLIB_COMPATIBLE_NAMES \
              -I$(top_srcdir)/aux_build -DDATADIR=\"$(BRESEQRELATIVEDATADIR)/\"

AM_CXXFLAGS = -std=c++11

noinst_LTLIBRARIES = libbreseq.la
libbreseq_la_SOURCES = \
	alignment.cpp \
	alignment_output.cpp \
	anyoption.cpp \
	calculate_trims.cpp \
	candidate_junctions.cpp \
	contingency_loci.cpp \
	coverage_output.cpp \
	coverage_distribution.cpp \
	error_count.cpp \
	fasta.cpp \
	fastq.cpp \
	flagged_regions.cpp \
	genome_diff.cpp \
	genome_diff_entry.cpp \
	identify_mutations.cpp \
	mutation_predictor.cpp \
	nw.cpp \
	output.cpp \
	pgzstream.cpp \
	pileup.cpp \
	pileup_base.cpp \
	reference_sequence.cpp \
	resolve_alignments.cpp \
	samtools_commands.cpp \
	settings.cpp \
	soft_clipping.cpp \
	stats.cpp \
	summary.cpp \
	libbreseq/alignment.h \
	libbreseq/alignment_output.h \
	libbreseq/anyoption.h \
	libbreseq/calculate_trims.h \
	libbreseq/candidate_junctions.h \
	libbreseq/common.h \
	libbreseq/counted_ptr.h \
	libbreseq/ctpl_stl.h \
	libbreseq/contingency_loci.h \
	libbreseq/coverage_distribution.h \
	libbreseq/coverage_output.h \
	libbreseq/error_count.h \
	libbreseq/fasta.h \
	libbreseq/fastq.h \
	libbreseq/file_parse_errors.h \
	libbreseq/flagged_regions.h \
	libbreseq/genome_diff.h \
	libbreseq/genome_diff_entry.h \
	libbreseq/identify_mutations.h \
	libbreseq/json.hpp \
	libbreseq/mutation_predictor.h \
	libbreseq/nw.h \
	libbreseq/output.h \
	libbreseq/pgzstream.h \
	libbreseq/pileup.h \
	libbreseq/pileup_base.h \
	libbreseq/reference_sequence.h \
	libbreseq/resolve_alignments.h \
	libbreseq/samtools_commands.h \
	libbreseq/settings.h \
	libbreseq/soft_clipping.h \
	libbreseq/stats.h \
	libbreseq/storable.h \
	libbreseq/summary.h
libbreseq_la_LIBADD =

bin_PROGRAMS = breseq gdtools

# $(CRAM_CODEC_LIBS) carries full .a paths to htslib's CRAM-codec
# dependencies (bzip2/xz) when binarydist.sh links statically — kept out of
# libbreseq_la_LIBADD/LIBS deliberately, since libtool warns "Linking the
# shared library libbreseq.la against the static library ... is not portable"
# whenever a noinst_LTLIBRARIES link sees a bare .a path; these are only
# needed (and only valid) on the final program links, after libhts.a.
breseq_SOURCES = breseq_cmdline.cpp
breseq_LDADD   = ./libbreseq.la $(HTSLIB_LINK) $(MINIZ_LINK) $(LIBZ_LINK) $(CRAM_CODEC_LIBS)

gdtools_SOURCES = gdtools_cmdline.cpp
gdtools_LDADD   = ./libbreseq.la $(HTSLIB_LINK) $(MINIZ_LINK) $(LIBZ_LINK) $(CRAM_CODEC_LIBS)

EXTRA_DIST = breseq.png \
	breseq_icon.png \
	list.min.js \
	jszip.min.js

install-data-local :

	##copy our data files
	mkdir -p $(BRESEQDATADIR)
	cp $(srcdir)/breseq.png $(BRESEQDATADIR)
	cp $(srcdir)/breseq_icon.png $(BRESEQDATADIR)
	cp $(srcdir)/list.min.js $(BRESEQDATADIR)
	cp $(srcdir)/jszip.min.js $(BRESEQDATADIR)

uninstall-local :

	##remove our data files
	rm -rf $(BRESEQDATADIR)

