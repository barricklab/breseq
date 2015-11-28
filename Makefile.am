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

PACKAGEDIR = $(top_srcdir)/package
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
	tests/data/lambda/lambda.1-2.gbk \
	tests/data/lambda/lambda.3.gbk \
	tests/data/lambda/lambda.4.gbk \
	tests/data/lambda/lambda.5.gbk \
	src/doc/Makefile \
	src/doc/conf.py \
	src/doc/acknowledgments.rst \
	src/doc/faq.rst \
	src/doc/gd_format.rst \
	src/doc/gd_usage.rst \
	src/doc/hidden_code_block.py \
	src/doc/index.rst \
	src/doc/installation.rst \
	src/doc/introduction.rst \
	src/doc/methods.rst \
	src/doc/output.rst \
	src/doc/references.rst \
	src/doc/test_drive.rst \
	src/doc/tutorial_barcoded_targeted.rst \
	src/doc/tutorial_clones.rst \
	src/doc/tutorial_populations.rst \
	src/doc/usage.rst \
	src/doc/images \
	src/doc/images/Ara-3.tree.png \
	src/doc/images/ZDB83.ECB02816-ECB02836.JC.png \
	src/doc/images/ZDB83.ECB02816-ECB02836.MC.png \
	src/doc/images/ZDB83.rnk-citG.JC.png \
	src/doc/images/ZDB83.rnk-citG.JC.schematic.png \
	src/doc/images/ZDB83.rnk-citG.zoom-in.png \
	src/doc/images/ZDB83.rrlA.MC.png \
	src/doc/images/ZDB83.waaT.IS-mediated-deletion.png \
	src/doc/images/amp_1.png \
	src/doc/images/breseq.png \
	src/doc/images/coverage_distribution.png \
	src/doc/images/del_1.png \
	src/doc/images/del_2.png \
	src/doc/images/end_trimming_example.docx \
	src/doc/images/end_trimming_example.png \
	src/doc/images/error_rates.png \
	src/doc/images/ins_1.png \
	src/doc/images/ins_2.png \
	src/doc/images/jc_1.png \
	src/doc/images/jc_2.png \
	src/doc/images/mc_1.png \
	src/doc/images/missed_mutation_no_trimming.png \
	src/doc/images/mob_1.png \
	src/doc/images/mob_2.png \
	src/doc/images/mob_3.png \
	src/doc/images/ra_1.png \
	src/doc/images/region_coverage_example.png \
	src/doc/images/snp_1.png \
	src/doc/images/snp_2.png \
	src/doc/images/sub_1.png


## tests - does not work yet.
##TESTS_ENVIRONMENT = PROGRAMDATAPATH=$(BRESEQ) SAMTOOLSPATH=$(SAMTOOLSDIR)
##TESTS = tests/lambda_mult_ref_read/testcmd.sh

DOCSRCDIR=$(top_srcdir)/src/doc
DOCBUILDDIR=$(abs_top_builddir)/src/doc/_build

PERLGENOMEDIFF=$(top_srcdir)/src/perl/GenomeDiff
SAMTOOLSDIR=$(top_srcdir)/$(LIBBAM)

## Main hooks called after normal make steps

all-local :

## We should build Perl modules here (and not during install) because
## otherwise a sudo make install creates superuser owned files in the
## distribution, that are not removed by make clean.
	
	## OLD method: 'make test' after 'make install'
	## This saves the install path for the test script
	##echo "export TESTBINPREFIX=$(bindir)" > tests/test.config;
	
	
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

docs :
## This environmental variable is used by python Sphinx
	echo $(PACKAGE_VERSION); \
	export BRESEQ_VERSION=$(PACKAGE_VERSION); \
	cd $(DOCSRCDIR) ; \
	make -f Makefile BUILDDIR=$(DOCBUILDDIR) html ;
	
dist-hook: docs

## This environmental variable is used by python Sphinx
	
	cp -r $(DOCBUILDDIR)/html $(distdir)/documentation ;
	
	## Uncomment below to add Perl back to dist
	
	## GenomeDiff
	##cd $(abs_top_srcdir)/src/perl/GenomeDiff ; \
	##chmod u+w . ; \
	##perl Build.PL ; \
	##./Build distdir ; 
	
	##mkdir -p $(distdir)/src/perl
	##cp -r $(abs_top_srcdir)/src/perl/GenomeDiff/GenomeDiff-* $(distdir)/src/perl/GenomeDiff ; 
	
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
	## NEW method: 'make test' after 'make'
	## This must be located here so that it is run after EXTRA_DIST files are included
	echo "export TESTBINPREFIX=$(top_srcdir)/src/c/breseq" > $(top_srcdir)/tests/test.config;
	echo "export BRESEQ_DATA_PATH=$(top_srcdir)/src/c/breseq" >> $(top_srcdir)/tests/test.config;
	echo "export BRESEQ_SAMTOOLS_PATH=$(top_srcdir)/extern/samtools-0.1.18" >> $(top_srcdir)/tests/test.config;
	$(top_srcdir)/tests/test.sh test tests


clean-tests:
	##This is commented out to not break make distcheck, which requires we not modify the source in this way?
	##echo "export TESTBINPREFIX=$(top_srcdir)/src/c/breseq" > $(top_srcdir)/tests/test.config;
	$(top_srcdir)/tests/test.sh clean tests

distclean-local:
	rm -rf $(DOCBUILDDIR);

	
