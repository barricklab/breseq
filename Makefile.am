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
	tests/data/lambda/empty.fastq \
	tests/data/lambda/only_bad.fastq \
	tests/data/lambda/lambda_mixed_population.1.fastq \
	tests/data/lambda/lambda_mixed_population.2.fastq \
	tests/data/lambda/lambda_mixed_population.3.fastq \
	tests/data/lambda/lambda_mixed_population.4.fastq \
	tests/data/lambda/lambda_mixed_population.5.fastq \
	tests/data/lambda/lambda.1-2.gbk \
	tests/data/lambda/lambda.3.gbk \
	tests/data/lambda/lambda.4.gbk \
	tests/data/lambda/lambda.5.gbk

## tests - does not work yet.
##TESTS_ENVIRONMENT = PROGRAMDATAPATH=$(BRESEQ) SAMTOOLSPATH=$(SAMTOOLSDIR)
##TESTS = tests/lambda_mult_ref_read/testcmd.sh

DOCSRCDIR=$(top_srcdir)/src/doc
DOCBUILDDIR=$(abs_top_builddir)/src/doc/_build

SAMTOOLSDIR=$(top_srcdir)/$(LIBBAM)

## Main hooks called after normal make steps

all-local :


install-exec-local :


clean-local :: clean-tests

	## This saves the install path for the test script
	rm -f $(abs_top_builddir)/tests/test.config
	rm -rf $(DOCBUILDDIR);
	find $(abs_top_builddir) -name '.DS_Store' -type f -delete
	rm -rf $(abs_top_builddir)/src/c/breseq/*.dSYM

docs :
## This environmental variable is used by python Sphinx
	echo $(PACKAGE_VERSION); \
	export BRESEQ_VERSION=$(PACKAGE_VERSION); \
	cd $(DOCSRCDIR) ; \
	make -f Makefile BUILDDIR=$(DOCBUILDDIR) html ;

dist-hook:

	if [ -d "$(DOCSRCDIR)" ]; then \
		make docs ; \
		cp -r $(DOCBUILDDIR)/html $(distdir)/documentation ; \
	fi ;

	find $(distdir) -name '.DS_Store' -type f -delete

	## Remove ._ files
	rm -rf `find $(distdir) -name ".*"`
	rm -rf ./src/c/breseq/gdtools.dSYM

	## Set up run_tests.sh and environmental variables
	echo "export TESTBINPREFIX=$(top_srcdir)/src/c/breseq" > $(distdir)/tests/test.config;
	echo "export BRESEQ_DATA_PATH=$(top_srcdir)/src/c/breseq" >> $(distdir)/tests/test.config;
	echo "export BRESEQ_TEST_THREAD_ARG=\"-j 4\"" >> $(distdir)/tests/test.config;

	echo "tests/test.sh clean tests" > ${distdir}/run_tests.sh
	echo "tests/test.sh test tests" >> ${distdir}/run_tests.sh
	chmod a+x ${distdir}/run_tests.sh

test:
	## NEW method: 'make test' after 'make'
	## This must be located here so that it is run after EXTRA_DIST files are included
	echo "export TESTBINPREFIX=$(top_srcdir)/src/c/breseq" > $(top_srcdir)/tests/test.config;
	echo "export BRESEQ_DATA_PATH=$(top_srcdir)/src/c/breseq" >> $(top_srcdir)/tests/test.config;
	echo "export BRESEQ_TEST_THREAD_ARG=\"-j 4\"" >> $(top_srcdir)/tests/test.config;
	$(top_srcdir)/tests/test.sh test tests

test-long:
	## NEW method: 'make test' after 'make'
	## This must be located here so that it is run after EXTRA_DIST files are included
	echo "export TESTBINPREFIX=$(top_srcdir)/src/c/breseq" > $(top_srcdir)/tests/test.config;
	echo "export BRESEQ_DATA_PATH=$(top_srcdir)/src/c/breseq" >> $(top_srcdir)/tests/test.config;
	echo "export BRESEQ_TEST_THREAD_ARG=\"-j 4\"" >> $(top_srcdir)/tests/test.config;
	$(top_srcdir)/tests/test.sh test tests long


clean-tests:
	##This is commented out to not break make distcheck, which requires we not modify the source in this way?
	##echo "export TESTBINPREFIX=$(top_srcdir)/src/c/breseq" > $(top_srcdir)/tests/test.config;
	$(top_srcdir)/tests/test.sh clean tests long

