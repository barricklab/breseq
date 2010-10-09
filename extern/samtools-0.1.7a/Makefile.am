#lib_LTLIBRARIES = libbam.la
#libbam_la_SOURCES = \
#	bgzf.c \
#	kstring.c \
#	bam_aux.c \
#	bam.c \
#	bam_import.c \
#	sam.c \
#	bam_index.c \
#	bam_pileup.c \
#	bam_lpileup.c \
#	bam_md.c \
#	glf.c \
#	razf.c \
#	faidx.c \
#	knetfile.c \
#	bam_sort.c \
#	sam_header.c
#libbam_la_CFLAGS = -fPIC -static

lib_LIBRARIES = libbam.a
libbam_a_SOURCES = \
	bgzf.c \
	kstring.c \
	bam_aux.c \
	bam.c \
	bam_import.c \
	sam.c \
	bam_index.c \
	bam_pileup.c \
	bam_lpileup.c \
	bam_md.c \
	glf.c \
	razf.c \
	faidx.c \
	knetfile.c \
	bam_sort.c \
	sam_header.c
libbam_a_CFLAGS = -fPIC -static

bin_PROGRAMS = samtools
samtools_SOURCES = \
	bam_tview.c \
	bam_maqcns.c \
	bam_plcmd.c \
	sam_view.c	\
	bam_rmdup.c \
	bam_rmdupse.c \
	bam_mate.c \
	bam_stat.c \
	bam_color.c	\
	bamtk.c \
	kaln.c
samtools_CFLAGS = -D_CURSES_LIB=1
samtools_LDADD = libbam.la

include_HEADERS = bam.h


