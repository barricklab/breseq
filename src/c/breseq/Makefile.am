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

noinst_LTLIBRARIES = libbreseq.la
libbreseq_la_SOURCES = \
	alignment.cpp \
	annotated_sequence.cpp \
	error_count.cpp \
	fastq.cpp \
	genome_diff.cpp \
	identify_mutations.cpp \
	pileup_base.cpp \
	pileup.cpp \
	tabulate_coverage.cpp \
	breseq/alignment.h \
	breseq/annotated_sequence.h \
	breseq/common.h \
	breseq/error_count.h \
	breseq/fastq.h \
	breseq/genome_diff.h \
	breseq/identify_mutations.h \
	breseq/pileup_base.h \
	breseq/pileup.h \
	breseq/tabulate_coverage.h 
libbreseq_la_LDFLAGS = -static
libbreseq_la_LIBADD = $(top_srcdir)/$(LIBBAM)/libbam.la

bin_PROGRAMS = error_count identify_mutations calc_trims tabulate_coverage convert_genbank cfastq_utils

error_count_SOURCES = error_count_cmdline.cpp

identify_mutations_SOURCES = identify_mutations_cmdline.cpp

calc_trims_SOURCES = calc_trims_cmdline.cpp

tabulate_coverage_SOURCES = tabulate_coverage_cmdline.cpp

convert_genbank_SOURCES = convert_genbank_cmdline.cpp

cfastq_utils_SOURCES = fastq_utils_cmdline.cpp

