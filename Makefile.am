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

PACKAGEDIR = $(ROOTDIR)/package
MACPACKAGEDIR = $(PACKAGEDIR)/macosx
MACPACKAGEBUILDDIR = $(MACPACKAGEDIR)/build

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
BIOSAMTOOLS=$(ROOTDIR)/extern/Bio-SamTools-1.28
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
	./Build install --destdir=${DESTDIR};
	
	## Breseq
	cd $(PERLBRESEQ) ; \
	./Build install --destdir=${DESTDIR};
	
clean-local :: clean-tests clean-breseq
	
	cd $(BIOSAMTOOLS) ; \
	./Build clean
	
		
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
make-breseq :
	cd $(PERLBRESEQ) ; \
	perl Build.PL --install_base=$(prefix) ; \
	./Build ; 

install-breseq :
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

##

macosx-pkg-download : 

	mkdir -p $(PACKAGEDIR)
	mkdir -p $(MACPACKAGEDIR)

	cd $(MACPACKAGEDIR) ; \
	curl -O http://cran.r-project.org/src/base/R-2/R-2.13.0.tar.gz ; \
	tar -xvzf R-2.13.0.tar.gz ;

	cd $(MACPACKAGEDIR) ; \
	curl -O http://mirror.datapipe.net/CPAN/authors/id/C/CJ/CJFIELDS/BioPerl-1.6.900.tar.gz ; \
	tar -xvzf BioPerl-1.6.900.tar.gz;
	
	cd $(MACPACKAGEDIR) ; \
	curl -O ftp://ftp.sanger.ac.uk/pub4/resources/software/ssaha2/ssaha2_mac.tgz ; \
	tar -xvzf ssaha2_mac.tgz;

macosx-pkg-compile :

	mkdir -p $(MACPACKAGEBUILDDIR)

	## We have a problem that when building an XS module
	## the architecture differs between the system Perl /usr/bin/perl
	## and MacPorts Perl, so we rebuild using the system Perl here
	cd $(BIOSAMTOOLS) ; \
	export SAMTOOLSLIBRARY=$(SAMTOOLSDIR)/.libs; \
	export SAMTOOLSHEADERS=$(SAMTOOLSDIR); \
	/usr/bin/perl Build.PL --install_base=$(MACPACKAGEBUILDDIR) --install_path bin=$(MACPACKAGEBUILDDIR)/bin --install_path script=$(MACPACKAGEBUILDDIR)/bin ; \
	./Build clean; \
	./Build ; \
	./Build install ;

	cd $(MACPACKAGEDIR)/R-2.13.0; \
	./configure --prefix=$(MACPACKAGEBUILDDIR) G77=gfortran --without-aqua --disable-R-framework ; \
	make ; \
	make install ;

	cd $(MACPACKAGEDIR)/BioPerl-1.6.900; \
	perl Build.PL --accept --install_base=$(MACPACKAGEBUILDDIR) ; \
	./Build install ;	
	
	cp $(MACPACKAGEDIR)/ssaha2*/ssaha2* $(MACPACKAGEBUILDDIR)/bin ;
	
	## Install outselves into this dir
	./configure --prefix=$(MACPACKAGEBUILDDIR) ; \
	make install ;

	
clean-pkg :

	rm -rf package

## tests

test:
	export TESTBINPREFIX=$(bindir); tests/test.sh test tests

clean-tests:
	export TESTBINPREFIX=$(bindir); tests/test.sh clean tests
