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

##This is the most compatible way of passing through the data directory to the C++ code
BRESEQDATADIR= $(DESTDIR)$(datadir)/$(PACKAGE_NAME)
BRESEQRELATIVEDATADIR= $(BIN_TO_DATA_RELATIVE_PATH)/$(PACKAGE_NAME)

# If configured for static libz linking, add include path here
AM_CPPFLAGS = -O3 -std=c++11 $(STATIC_LIBZ_INC) -I$(top_srcdir)/$(LIBBAM) -I$(top_srcdir)/$(HTSLIB) -I$(top_srcdir)/aux_build -DDATADIR=\"$(BRESEQRELATIVEDATADIR)/\"

noinst_LTLIBRARIES = libbreseq.la
libbreseq_la_SOURCES = \
	alignment.cpp \
	alignment_output.cpp \
	anyoption.cpp \
	calculate_trims.cpp \
	candidate_junctions.cpp \
	chisquare.cpp \
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
	pileup.cpp \
	pileup_base.cpp \
	reference_sequence.cpp \
	resolve_alignments.cpp \
	samtools_commands.cpp \
	settings.cpp \
	soft_clipping.cpp \
	summary.cpp \
	libbreseq/alignment.h \
	libbreseq/alignment_output.h \
	libbreseq/anyoption.h \
	libbreseq/calculate_trims.h \
	libbreseq/candidate_junctions.h \
	libbreseq/chisquare.h \
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
	libbreseq/gzstream.h \
	libbreseq/identify_mutations.h \
	libbreseq/json.hpp \
	libbreseq/mutation_predictor.h \
	libbreseq/nw.h \
	libbreseq/output.h \
	libbreseq/pileup.h \
	libbreseq/pileup_base.h \
	libbreseq/reference_sequence.h \
	libbreseq/resolve_alignments.h \
	libbreseq/samtools_commands.h \
	libbreseq/settings.h \
	libbreseq/soft_clipping.h \
	libbreseq/storable.h \
	libbreseq/summary.h
libbreseq_la_LIBADD = $(UNWIND_LIBS) $(top_builddir)/$(HTSLIB)/libhts.la $(top_builddir)/$(LIBBAM)/libbam.la

# If configured for static libz linking, add full path to libz.a to force static linking
LDADD = libbreseq.la $(STATIC_LIBZ_LIB)

bin_PROGRAMS = breseq gdtools
breseq_SOURCES = breseq_cmdline.cpp
gdtools_SOURCES = gdtools_cmdline.cpp

EXTRA_DIST = breseq.png \
	breseq_small.png \
	coverage_distribution.r \
	plot_coverage.r \
	plot_error_rate.r \
	plot_jc_scores.r \
	polymorphism_statistics.r \
	list.min.js \
	run_circos.sh \
	indels.conf \
	mutations.conf \
	mobs.conf \
	combined_circos.conf

install-data-local :

	##copy our data files
	mkdir -p $(BRESEQDATADIR)
	cp $(srcdir)/breseq.png $(BRESEQDATADIR)
	cp $(srcdir)/breseq_small.png $(BRESEQDATADIR)
	cp $(srcdir)/coverage_distribution.r $(BRESEQDATADIR)
	cp $(srcdir)/plot_coverage.r $(BRESEQDATADIR)
	cp $(srcdir)/plot_error_rate.r $(BRESEQDATADIR)
	cp $(srcdir)/plot_jc_scores.r $(BRESEQDATADIR)
	cp $(srcdir)/polymorphism_statistics.r $(BRESEQDATADIR)
	cp $(srcdir)/list.min.js $(BRESEQDATADIR)
	cp $(srcdir)/run_circos.sh $(BRESEQDATADIR)
	cp $(srcdir)/indels.conf $(BRESEQDATADIR)
	cp $(srcdir)/mutations.conf $(BRESEQDATADIR)
	cp $(srcdir)/mobs.conf $(BRESEQDATADIR)
	cp $(srcdir)/combined_circos.conf $(BRESEQDATADIR)

uninstall-local :

	##remove our data files
	rm -rf $(BRESEQDATADIR)

