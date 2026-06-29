###
##
## AUTHORS
##
##   Jeffrey E. Barrick <jeffrey.e.barrick@gmail.com> and other contributors
##
## LICENSE AND COPYRIGHT
##
##   Copyright (c) 2008-2010 Michigan State University
##   Copyright (c) 2011-2025 The University of Texas at Austin
##   Copyright (c) 2025-     Michigan State University
##
##   breseq is free software; you can redistribute it and/or modify it under the
##   terms of the GNU General Public License as published by the Free Software
##   Foundation; either version 2, or (at your option) any later version.
##
##   SPDX-License-Identifier: GPL-2.0-or-later
##
###

ACLOCAL_AMFLAGS = -I aux_build/m4

SUBDIRS = $(BRESEQ)

### Tests! -- only one test included for now in distribution.

## Top level files included (and required) by default
## AUTHORS, ChangeLog, COPYING, INSTALL, NEWS, README
EXTRA_DIST = \
	DEVELOPER \
	dev-environment.yml \
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

DOCSRCDIR=$(top_srcdir)/src/doc
DOCBUILDDIR=$(abs_top_builddir)/src/doc/_build

## Total core budget for parallel test execution.
## Auto-detected from make's -j flag by inspecting make's own command line via
## ps (works with all make versions, including macOS's bundled GNU Make 3.81
## which does not expose the job count via MAKEFLAGS). Falls back to 'all'
## (all available cores) when no -j was given or ps is unavailable.
## Override explicitly on the command line: 'make test TEST_CORES=16'.
## Individual tests declare how many of these cores they use via 'TEST_CORES=N'
## in their testcmd.sh (see tests/common.sh and tests/Snakefile).
TEST_CORES ?= $(or $(shell ps -p $$PPID -o args= 2>/dev/null | grep -oE -- '(-j *[0-9]+|--jobs=[0-9]+)' | grep -oE '[0-9]+' | head -1),all)

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

	echo "tests/test.sh clean tests" > ${distdir}/run_tests.sh
	echo "tests/test.sh test tests" >> ${distdir}/run_tests.sh
	chmod a+x ${distdir}/run_tests.sh

## 'make test'/'make test-long' run the consistency tests in parallel via
## Snakemake (see tests/Snakefile), bounded by a total core budget
## (TEST_CORES, auto-detected from make -j or 'all' if unset;
## override with e.g. 'make test TEST_CORES=16').
## Each test writes its output to tests/<name>/test.log (removed by
## 'make clean-tests') and a one-line tests/<name>/test.result that
## tests/print_test_summary.sh aggregates into the final pass/fail + timing
## summary. The overall exit status reflects whether every test passed, so
## this is suitable for driving CI (e.g. GitHub Actions).
test:
	## NEW method: 'make test' after 'make'
	## This must be located here so that it is run after EXTRA_DIST files are included
	## NOTE: keep these as paths relative to top_srcdir (not absolute) -- the
	## literal command line (including these paths) gets embedded as the
	## '#=COMMAND' line in each test's output .gd file, and expected.gd files
	## were generated this way. The Snakefile runs tests from this directory
	## (it does not change the working directory), so the paths still resolve.
	echo "export TESTBINPREFIX=$(top_srcdir)/src/c/breseq" > $(top_srcdir)/tests/test.config;
	echo "export BRESEQ_DATA_PATH=$(top_srcdir)/src/c/breseq" >> $(top_srcdir)/tests/test.config;
	START=$$(date '+%s'); \
	snakemake -s $(top_srcdir)/tests/Snakefile --cores $(TEST_CORES) --keep-going; \
	SNAKEMAKE_STATUS=$$?; \
	END=$$(date '+%s'); \
	$(top_srcdir)/tests/print_test_summary.sh $(top_srcdir)/tests $$((END - START)); \
	SUMMARY_STATUS=$$?; \
	test $$SNAKEMAKE_STATUS -eq 0 -a $$SUMMARY_STATUS -eq 0

test-long:
	## NEW method: 'make test' after 'make'
	## This must be located here so that it is run after EXTRA_DIST files are included
	## NOTE: see comment in the 'test' target above about why these paths are relative.
	echo "export TESTBINPREFIX=$(top_srcdir)/src/c/breseq" > $(top_srcdir)/tests/test.config;
	echo "export BRESEQ_DATA_PATH=$(top_srcdir)/src/c/breseq" >> $(top_srcdir)/tests/test.config;
	START=$$(date '+%s'); \
	snakemake -s $(top_srcdir)/tests/Snakefile --cores $(TEST_CORES) --keep-going --config long=1; \
	SNAKEMAKE_STATUS=$$?; \
	END=$$(date '+%s'); \
	$(top_srcdir)/tests/print_test_summary.sh $(top_srcdir)/tests $$((END - START)); \
	SUMMARY_STATUS=$$?; \
	test $$SNAKEMAKE_STATUS -eq 0 -a $$SUMMARY_STATUS -eq 0


clean-tests:
	##This is commented out to not break make distcheck, which requires we not modify the source in this way?
	##echo "export TESTBINPREFIX=$(top_srcdir)/src/c/breseq" > $(top_srcdir)/tests/test.config;
	$(top_srcdir)/tests/test.sh clean tests long
	rm -rf $(top_srcdir)/tests/.snakemake

