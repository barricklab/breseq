lib_LTLIBRARIES = libbreseq.la
libbreseq_la_SOURCES = \
	alignment.cpp \
	error_count.cpp \
	genome_diff.cpp \
	identify_mutations.cpp \
	pileup.cpp \
	pileup_base.cpp
include_HEADERS = \
	breseq/alignment.h \
	breseq/common.h \
	breseq/error_count.h \
	breseq/genome_diff.h \
	breseq/identify_mutations.h \
	breseq/pileup.h \
	breseq/pileup_base.h

bin_PROGRAMS = error_count identify_mutations calc_trims
error_count_SOURCES = error_count_cmdline.cpp
error_count_LDADD = libbreseq.la

identify_mutations_SOURCES = identify_mutations_cmdline.cpp
identify_mutations_LDADD = libbreseq.la

calc_trims_SOURCES = calc_trims_cmdline.cpp
calc_trims_LDADD = libbreseq.la
