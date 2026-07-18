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
	cn_evidence.cpp \
	contingency_loci.cpp \
	coverage_output.cpp \
	coverage_distribution.cpp \
	dp_evidence.cpp \
	error_count.cpp \
	fasta.cpp \
	fastq.cpp \
	flagged_regions.cpp \
	genome_diff.cpp \
	genome_diff_entry.cpp \
	identify_mutations.cpp \
	mutation_predictor.cpp \
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
	alignment.h \
	alignment_output.h \
	anyoption.h \
	calculate_trims.h \
	candidate_junctions.h \
	cn_evidence.h \
	common.h \
	counted_ptr.h \
	ctpl_stl.h \
	contingency_loci.h \
	coverage_distribution.h \
	coverage_output.h \
	dp_evidence.h \
	error_count.h \
	fasta.h \
	fastq.h \
	file_parse_errors.h \
	flagged_regions.h \
	genome_diff.h \
	genome_diff_entry.h \
	identify_mutations.h \
	json.hpp \
	mutation_predictor.h \
	output.h \
	pgzstream.h \
	pileup.h \
	pileup_base.h \
	reference_sequence.h \
	resolve_alignments.h \
	samtools_commands.h \
	settings.h \
	soft_clipping.h \
	stats.h \
	storable.h \
	summary.h
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

## These live in ../share/breseq/ (a sibling of this directory) so that
## Settings::set_global_paths()'s DATADIR-relative fallback -- which assumes
## a bin/ and share/breseq/ laid out side by side, matching the offset
## between --bindir and --datadir -- finds them automatically for an
## uninstalled in-tree build, not just an installed one.
EXTRA_DIST = ../share/breseq/breseq.png \
	../share/breseq/breseq_icon.png \
	../share/breseq/list.min.js \
	../share/breseq/jszip.min.js

install-data-local :

	##copy our data files
	mkdir -p $(BRESEQDATADIR)
	cp $(srcdir)/../share/breseq/breseq.png $(BRESEQDATADIR)
	cp $(srcdir)/../share/breseq/breseq_icon.png $(BRESEQDATADIR)
	cp $(srcdir)/../share/breseq/list.min.js $(BRESEQDATADIR)
	cp $(srcdir)/../share/breseq/jszip.min.js $(BRESEQDATADIR)

uninstall-local :

	##remove our data files
	rm -rf $(BRESEQDATADIR)

