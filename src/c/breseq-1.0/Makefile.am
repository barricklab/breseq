AM_CPPFLAGS = -I$(top_srcdir)/$(LIBBAM) $(BOOST_CPPFLAGS)
AM_LDFLAGS = $(BOOST_PROGRAM_OPTIONS_LDFLAGS)
LDADD = $(top_srcdir)/$(LIBBAM)/libbam.la \
	libbreseq.la \
	$(BOOST_PROGRAM_OPTIONS_LIBS)

noinst_LTLIBRARIES = libbreseq.la
libbreseq_la_SOURCES = \
	alignment.cpp \
	error_count.cpp \
	genome_diff.cpp \
	identify_mutations.cpp \
	pileup.cpp \
	pileup_base.cpp \
	breseq/alignment.h \
	breseq/common.h \
	breseq/error_count.h \
	breseq/genome_diff.h \
	breseq/identify_mutations.h \
	breseq/pileup_base.h \
	breseq/pileup.h 
libbreseq_la_LIBADD = $(top_srcdir)/$(LIBBAM)/libbam.la

bin_PROGRAMS = error_count identify_mutations calc_trims
error_count_SOURCES = error_count_cmdline.cpp
identify_mutations_SOURCES = identify_mutations_cmdline.cpp
calc_trims_SOURCES = calc_trims_cmdline.cpp
