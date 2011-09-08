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

## tests - does not work yet.
##TESTS_ENVIRONMENT = PROGRAMDATAPATH=$(BRESEQ) SAMTOOLSPATH=$(SAMTOOLSDIR)
##TESTS = tests/lambda_mult_ref_read/testcmd.sh

ROOTDIR=$(PWD)
DOCDIR=$(ROOTDIR)/src/doc
PERLGENOMEDIFF=$(ROOTDIR)/src/perl/GenomeDiff
SAMTOOLSDIR=$(ROOTDIR)/$(LIBBAM)

## Main hooks called after normal make steps

all-local :

## We should build Perl modules here (and not during install) because
## otherwise a sudo make install creates superuser owned files in the
## distribution, that are not removed by make clean.
	
	## This saves the install path for the test script
	echo "export TESTBINPREFIX=$(bindir)" > tests/test.config;

install-exec-local :

	
clean-local :: clean-tests
	
	## This saves the install path for the test script	
	rm -f tests/test.config

## These two targets are optional for Perl genomediff tools.

genomediff :

	## GenomeDiff
	cd $(PERLGENOMEDIFF) ; \
	perl Build.PL --install_base=$(prefix) --install_path bin=$(bindir)  --install_path script=$(bindir) ; \
	./Build install;
	
clean-genomediff :

	## GenomeDiff
	cd $(PERLGENOMEDIFF) ; \
	./Build clean;	
	
dist-hook:

## This environmental variable is used by python Sphinx
	echo $(PACKAGE_VERSION); \
	export BRESEQ_VERSION=$(PACKAGE_VERSION); \
	cd $(DOCDIR) ; \
	make html ;
	
	cp -r $(DOCDIR)/_build/html $(distdir)/documentation ;

	## GenomeDiff
	cd $(PERLGENOMEDIFF) ; \
	./Build distdir ; \
	echo $(distdir)/src/perl/GenomeDiff ;
	
	mkdir -p $(distdir)/src/perl
	cp -r $(PERLGENOMEDIFF)/GenomeDiff-* $(distdir)/src/perl/GenomeDiff ; 
	
	## Remove ._ files
	rm -rf `find $(distdir) -name ".*"`
			

macosx-pkg-download : 

	mkdir -p $(PACKAGEDIR)
	mkdir -p $(MACPACKAGEDIR)
	
	cd $(MACPACKAGEDIR) ; \
	curl -O ftp://ftp.sanger.ac.uk/pub4/resources/software/ssaha2/ssaha2_mac.tgz ; \
	tar -xvzf ssaha2_mac.tgz;


macosx-pkg-compile :

	mkdir -p $(MACPACKAGEBUILDDIR)
	
	cp $(MACPACKAGEDIR)/ssaha2*/ssaha2* $(MACPACKAGEBUILDDIR)/bin ;
	
	## Install outselves into this dir
	./configure --prefix=$(MACPACKAGEBUILDDIR) ; \
	make install ;
	
	
clean-pkg :
	rm -rf package


test:
	tests/test.sh test tests


clean-tests:
	tests/test.sh clean tests
