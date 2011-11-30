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

##is there a variable for this?
BRESEQDATADIR=$(datadir)/$(PACKAGE_NAME)
AM_CPPFLAGS = -I$(top_srcdir)/$(LIBBAM) -I$(top_srcdir)/aux_build  -DDATADIR=\"$(BRESEQDATADIR)\"
AM_LDFLAGS = -static
LDADD =   $(top_srcdir)/$(LIBBAM)/libbam.la \
					libbreseq.la
LIBS = 

noinst_LTLIBRARIES = libbreseq.la
libbreseq_la_SOURCES = \
	alignment.cpp \
	alignment_output.cpp \
	annotated_sequence.cpp \
	anyoption.cpp \
	calculate_trims.cpp \
	candidate_junctions.cpp \
	contingency_loci.cpp \
	coverage_output.cpp \
	coverage_distribution.cpp \
	chisquare.cpp \
	error_count.cpp \
	fasta.cpp \
	fastq.cpp \
	genome_diff.cpp \
	identify_mutations.cpp \
	mutation_predictor.cpp \
	output.cpp \
	pileup_base.cpp \
	pileup.cpp \
	resolve_alignments.cpp \
	settings.cpp \
	simulated_read.cpp \
	libbreseq/alignment.h \
	libbreseq/alignment_output.h \
	libbreseq/annotated_sequence.h \
	libbreseq/anyoption.h \
	libbreseq/calculate_trims.h \
	libbreseq/candidate_junctions.h \
	libbreseq/chisquare.h \
	libbreseq/contingency_loci.h \
	libbreseq/common.h \
	libbreseq/coverage_distribution.h \
	libbreseq/coverage_output.h \
	libbreseq/error_count.h \
	libbreseq/fasta.h \
	libbreseq/fastq.h \
	libbreseq/genome_diff.h \
	libbreseq/identify_mutations.h \
	libbreseq/mutation_predictor.h \
	libbreseq/output.h \
	libbreseq/pileup_base.h \
	libbreseq/pileup.h \
	libbreseq/resolve_alignments.h \
	libbreseq/settings.h \
	libbreseq/storable.h \
  libbreseq/simulated_read.h	
	
libbreseq_la_LIBADD = $(top_srcdir)/$(LIBBAM)/libbam.la

bin_PROGRAMS = breseq 
breseq_SOURCES = breseq_cmdline.cpp

EXTRA_DIST = \
	breseq.png \
	breseq_small.png \
	coverage_distribution.r \
	plot_coverage.r \
	plot_error_rate.r \
	polymorphism_statistics.r

install-data-local :

	##copy our data files
	mkdir -p $(BRESEQDATADIR)
	cp breseq.png $(BRESEQDATADIR)
	cp breseq_small.png $(BRESEQDATADIR)
	cp coverage_distribution.r $(BRESEQDATADIR)
	cp plot_coverage.r $(BRESEQDATADIR)
	cp plot_error_rate.r $(BRESEQDATADIR)
	cp polymorphism_statistics.r $(BRESEQDATADIR)
