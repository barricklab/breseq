### breseq Perl Module
#EXTRA_DIST = \
#	src/perl/Breseq/bin \
#	src/perl/Breseq/Build.PL \
#	src/perl/Breseq/Changes \
#	src/perl/Breseq/lib \
#	src/perl/Breseq/MANIFEST \
#	src/perl/Breseq/README \
#	src/perl/Breseq/t

### Bio::DB::Sam Perl Module 
#EXTRA_DIST += \
#	extern/Bio-SamTools-1.19/bin \
#	extern/Bio-SamTools-1.19/Build.PL \
#	extern/Bio-SamTools-1.19/Changes \
#	extern/Bio-SamTools-1.19/DISCLAIMER \
#	extern/Bio-SamTools-1.19/lib \
#	extern/Bio-SamTools-1.19/LICENSE \
#	extern/Bio-SamTools-1.19/MANIFEST \
#	extern/Bio-SamTools-1.19/README \
#	extern/Bio-SamTools-1.19/t \
#	extern/Bio-SamTools-1.19/typemap

### Tests!
EXTRA_DIST = tests

ACLOCAL_AMFLAGS = -I aux_build/m4

SUBDIRS = extern/samtools-0.1.7a src/c

ROOTDIR=$(PWD)
DOCDIR=$(ROOTDIR)/src/doc
BRESEQDIR=$(ROOTDIR)/src/perl/Breseq
BIOSAMTOOLS=$(ROOTDIR)/extern/Bio-SamTools-1.19
SAMTOOLSDIR=$(ROOTDIR)/extern/samtools-0.1.7a
BIOPERL=$(ROOTDIR)/extern/bioperl-live

STAGEDIR=$(ROOTDIR)/build

## Main hooks called after normal make steps

all-local :

## We should build Perl modules here (and not during install) because
## otherwise a sudo make install creates superuser owned files in the
## distribution, that are not removed by make clean.

	cd $(BRESEQDIR) ; \
	perl Build.PL --install_base=$(prefix) ; \
	./Build ;

	## Bio::SamTools
	cd $(BIOSAMTOOLS) ; \
	export SAMTOOLS=$(ROOTDIR)/extern/samtools-0.1.7a ; \
	perl Build.PL --install_base=$(prefix) ; \
	./Build ;

	## Bioperl
#	cd $(BIOPERL) ; \
#	perl Build.PL --accept --install_base=$(prefix) ; \
#	./Build ; 

install-exec-local :

	## Breseq
	cd $(BRESEQDIR) ; \
	./Build install;

#	cd $(BIOPERL) ; \
#	./Build install ;

	## Bio::SamTools
	cd $(BIOSAMTOOLS) ; \
	./Build install ;
	
clean-local :

	cd $(BRESEQDIR) ; \
	./Build clean

#	cd $(BIOPERL) ; \
#	./Build clean
	
	cd $(BIOSAMTOOLS) ; \
	./Build clean
	
	rm -rf $(STAGEDIR)
	
	tests/test.sh clean tests/	

dist-hook:

	cd $(DOCDIR) ; \
	make html ;
	
	cp -r $(DOCDIR)/_build/html $(distdir) ;

	cd $(BRESEQDIR) ; \
	./Build distdir ; \
	echo $(distdir)/src/perl/Breseq ;
	
	mkdir -p $(distdir)/src/perl
	cp -r $(BRESEQDIR)/Breseq-* $(distdir)/src/perl/Breseq ; 

	cd $(BIOSAMTOOLS) ; \
	./Build distdir ;
	
	cp -r $(BIOSAMTOOLS)/Bio-SamTools-* $(distdir)/extern ; 
	
	
## breseq only	
make-breseq:
	cd $(BRESEQDIR) ; \
	perl Build.PL --install_base=$(prefix) ; \
	./Build ; 

install-breseq:
	cd $(BRESEQDIR) ; \
	./Build install
	
clean-breseq :
	cd $(BRESEQDIR) ; \
	./Build clean
	
#utility to avoid everything else
breseq :

	cd $(BRESEQDIR) ; \
	perl Build.PL --install_base=$(prefix) ; \
	./Build

	cd $(BRESEQDIR) ; \
	./Build install
	
remake-breseq :: clean-breseq make-breseq install-breseq	

## tests

test:
	tests/lambda_mult_ref_read/testcmd.sh  test
	tests/lambda_mixed_pop/testcmd.sh  test
	tests/REL606_fragment_manipulated/testcmd.sh  test
	
test-all:
	tests/test.sh test tests

clean-tests:
	tests/test.sh clean tests
	
