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

ACLOCAL_AMFLAGS = -I aux_build/m4

SUBDIRS = $(LIBBAM) $(BRESEQ)

### Tests! -- only one test included for now in distribution.
EXTRA_DIST = \
	tests/common.sh \
	tests/test.sh \
	tests/lambda_mult_ref_read/expected.gd \
	tests/lambda_mult_ref_read/testcmd.sh \
	tests/data/lambda/lambda_mixed_population.1.fastq \
	tests/data/lambda/lambda_mixed_population.2.fastq \
	tests/data/lambda/lambda_mixed_population.3.fastq \
	tests/data/lambda/lambda_mixed_population.4.fastq \
	tests/data/lambda/lambda_mixed_population.5.fastq \
	tests/data/lambda/lambda.1.gbk \
	tests/data/lambda/lambda.2.gbk \
	tests/data/lambda/lambda.3.gbk \
	tests/data/lambda/lambda.4.gbk \
	tests/data/lambda/lambda.5.gbk

ROOTDIR=$(PWD)
DOCDIR=$(ROOTDIR)/src/doc
PERLBRESEQ=$(ROOTDIR)/src/perl/Breseq
BIOSAMTOOLS=$(ROOTDIR)/extern/Bio-SamTools-1.19
SAMTOOLSDIR=$(ROOTDIR)/$(LIBBAM)

## Main hooks called after normal make steps

all-local :

## We should build Perl modules here (and not during install) because
## otherwise a sudo make install creates superuser owned files in the
## distribution, that are not removed by make clean.

	## Bio::SamTools
	cd $(BIOSAMTOOLS) ; \
	export SAMTOOLSLIBRARY=$(SAMTOOLSDIR)/.libs; \
	export SAMTOOLSHEADERS=$(SAMTOOLSDIR); \
	perl Build.PL --install_base=$(prefix) --install_path bin=$(bindir) --install_path script=$(bindir) ; \
	./Build ;

	## Breseq
	cd $(PERLBRESEQ) ; \
	perl Build.PL --install_base=$(prefix) --install_path bin=$(bindir)  --install_path script=$(bindir) ; \
	./Build ;

install-exec-local :

	echo $(prefix)/lib

### Alternate way of doing things
##	## Bio::SamTools - compiled here because libbam must be in final location for static link to work!
##	cd $(BIOSAMTOOLS) ; \
##	export SAMTOOLSLIBRARY=$(libdir); \
##	export SAMTOOLSHEADERS=$(SAMTOOLSDIR); \
##	perl Build.PL --install_base=$(prefix) --install_path bin=$(bindir) --install_path script=$(bindir) ; \
##	./Build ; \
##	./Build install ;

	## Bio::SamTools
	cd $(BIOSAMTOOLS) ; \
	./Build install ;
	
	## Breseq
	cd $(PERLBRESEQ) ; \
	./Build install;
	
clean-local :
	
	cd $(BIOSAMTOOLS) ; \
	./Build clean
		
	cd $(PERLBRESEQ) ; \
	./Build clean
	
	tests/test.sh clean tests/	
	
dist-hook:

	cd $(DOCDIR) ; \
	make html ;
	
	cp -r $(DOCDIR)/_build/html $(distdir)/documentation ;

	cd $(BIOSAMTOOLS) ; \
	./Build distdir ;
	cp -r $(BIOSAMTOOLS)/Bio-SamTools-* $(distdir)/extern ; 

	cd $(PERLBRESEQ) ; \
	./Build distdir ; \
	echo $(distdir)/src/perl/Breseq ;
	
	mkdir -p $(distdir)/src/perl
	cp -r $(PERLBRESEQ)/Breseq-* $(distdir)/src/perl/Breseq ; 
	
	## Remove ._ files
	rm -rf `find $(distdir) -name ".*"`
		
## breseq only	
make-breseq:
	cd $(PERLBRESEQ) ; \
	perl Build.PL --install_base=$(prefix) ; \
	./Build ; 

install-breseq:
	cd $(PERLBRESEQ) ; \
	./Build install
	
clean-breseq :
	cd $(PERLBRESEQ) ; \
	./Build clean
	
#utility to avoid everything else
breseq :

	cd $(PERLBRESEQ) ; \
	perl Build.PL --install_base=$(prefix) ; \
	./Build

	cd $(PERLBRESEQ) ; \
	./Build install
	
remake-breseq :: clean-breseq make-breseq install-breseq	

## tests

test:
	export TESTBINPREFIX=$(bindir); tests/test.sh test tests

clean-tests:
	tests/test.sh clean tests
