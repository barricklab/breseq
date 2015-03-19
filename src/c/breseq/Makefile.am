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
BRESEQDATADIR= $(DESTDIR)/$(datadir)/$(PACKAGE_NAME)
BRESEQRELATIVEDATADIR= $(BIN_TO_DATA_RELATIVE_PATH)/$(PACKAGE_NAME)
AM_CPPFLAGS = -I$(top_srcdir)/$(LIBBAM) -I$(top_srcdir)/aux_build  -DDATADIR=\"$(BRESEQRELATIVEDATADIR)/\"
AM_LDFLAGS = -static
LDADD =   $(top_builddir)/$(LIBBAM)/libbam.la \
					libbreseq.la
LIBS = 

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
	identify_mutations.cpp \
	mutation_predictor.cpp \
	output.cpp \
	pileup.cpp \
	pileup_base.cpp \
	reference_sequence.cpp \
	resolve_alignments.cpp \
	rna_seq.cpp \
	settings.cpp \
	summary.cpp \
	libbreseq/alignment.h \
	libbreseq/alignment_output.h \
	libbreseq/anyoption.h \
	libbreseq/calculate_trims.h \
	libbreseq/candidate_junctions.h \
	libbreseq/chisquare.h \
	libbreseq/common.h \
	libbreseq/contingency_loci.h \
	libbreseq/coverage_distribution.h \
	libbreseq/coverage_output.h \
	libbreseq/error_count.h \
	libbreseq/fasta.h \
	libbreseq/fastq.h \
	libbreseq/flagged_regions.h \
	libbreseq/genome_diff.h \
	libbreseq/identify_mutations.h \
	libbreseq/mutation_predictor.h \
	libbreseq/output.h \
	libbreseq/pileup.h \
	libbreseq/pileup_base.h \
	libbreseq/reference_sequence.h \
	libbreseq/resolve_alignments.h \
	libbreseq/rna_seq.h \
	libbreseq/settings.h \
	libbreseq/storable.h \
	libbreseq/summary.h 
	
libbreseq_la_LIBADD = $(top_builddir)/$(LIBBAM)/libbam.la

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
	cp $(srcdir)/run_circos.sh $(BRESEQDATADIR)
	cp $(srcdir)/indels.conf $(BRESEQDATADIR)
	cp $(srcdir)/mutations.conf $(BRESEQDATADIR)
	cp $(srcdir)/mobs.conf $(BRESEQDATADIR)
	cp $(srcdir)/combined_circos.conf $(BRESEQDATADIR)

uninstall-local :

	##remove our data files
	rm -rf $(BRESEQDATADIR)

