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

AM_CPPFLAGS = -I$(top_srcdir)/$(LIBBAM) $(BOOST_CPPFLAGS)
AM_LDFLAGS = $(BOOST_PROGRAM_OPTIONS_LDFLAGS) -static
LDADD = $(BOOST_PROGRAM_OPTIONS_LIBS) \
		$(top_srcdir)/$(LIBBAM)/libbam.la \
		libbreseq.la
LIBS = $(BOOST_PROGRAM_OPTIONS_LIBS)

noinst_LTLIBRARIES = libbreseq.la
libbreseq_la_SOURCES = \
	alignment.cpp \
	alignment_output.cpp \
	annotated_sequence.cpp \
	calculate_trims.cpp \
	candidate_junctions.cpp \
	error_count.cpp \
	fasta.cpp \
	fastq.cpp \
	genome_diff.cpp \
	identify_mutations.cpp \
	pileup_base.cpp \
	pileup.cpp \
	resolve_alignments.cpp \
	tabulate_coverage.cpp \
	breseq/alignment.h \
	breseq/alignment_output.h \
	breseq/annotated_sequence.h \
	breseq/calculate_trims.h \
	breseq/candidate_junctions.h \
	breseq/common.h \
	breseq/error_count.h \
	breseq/fasta.h \
	breseq/fastq.h \
	breseq/genome_diff.h \
	breseq/identify_mutations.h \
	breseq/pileup_base.h \
	breseq/pileup.h \
	breseq/resolve_alignments.h \
	breseq/tabulate_coverage.h 
libbreseq_la_LIBADD = $(top_srcdir)/$(LIBBAM)/libbam.la

bin_PROGRAMS = cbreseq cbam2aln

cbreseq_SOURCES = breseq_cmdline.cpp
cbam2aln_SOURCES = bam2aln_cmdline.cpp