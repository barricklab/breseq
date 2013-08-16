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
BRESEQDATADIR= $(datadir)/$(PACKAGE_NAME)
AM_CPPFLAGS = -I$(top_srcdir)/$(LIBBAM) -I$(top_srcdir)/aux_build  -DDATADIR=\"$(BRESEQDATADIR)\"
AM_LDFLAGS = -static
LDADD =   $(top_srcdir)/$(LIBBAM)/libbam.la \
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
	
libbreseq_la_LIBADD = $(top_srcdir)/$(LIBBAM)/libbam.la

bin_PROGRAMS = breseq gdtools 
breseq_SOURCES = breseq_cmdline.cpp
gdtools_SOURCES = gdtools_cmdline.cpp

EXTRA_DIST = \
	breseq.png \
	breseq_small.png \
	coverage_distribution.r \
	plot_coverage.r \
	plot_error_rate.r \
	plot_jc_scores.r \
	polymorphism_statistics.r \
	run_circos.sh \
	ideogram.conf \
	indels.conf \
	karyotype.and.layout.conf \
	mobs.conf \
	mutations.conf \
	combined_circos.conf
	

install-data-local :

	##copy our data files
	mkdir -p $(DESTDIR)$(BRESEQDATADIR)
	cp breseq.png $(DESTDIR)$(BRESEQDATADIR)
	cp breseq_small.png $(DESTDIR)$(BRESEQDATADIR)
	cp coverage_distribution.r $(DESTDIR)$(BRESEQDATADIR)
	cp plot_coverage.r $(DESTDIR)$(BRESEQDATADIR)
	cp plot_error_rate.r $(DESTDIR)$(BRESEQDATADIR)
	cp plot_jc_scores.r $(DESTDIR)$(BRESEQDATADIR)
	cp polymorphism_statistics.r $(DESTDIR)$(BRESEQDATADIR)
	cp run_circos.sh $(DESTDIR)$(BRESEQDATADIR)
	cp indels.conf $(DESTDIR)$(BRESEQDATADIR)
	cp mutations.conf $(DESTDIR)$(BRESEQDATADIR)
	cp karyotype.and.layout.conf $(DESTDIR)$(BRESEQDATADIR)
	cp mobs.conf $(DESTDIR)$(BRESEQDATADIR)
	cp ideogram.conf $(DESTDIR)$(BRESEQDATADIR)
	cp combined_circos.conf $(DESTDIR)$(BRESEQDATADIR)
