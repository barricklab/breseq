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
	tests/data/lambda/lambda_mixed_population.A.fastq.gz \
	tests/data/lambda/lambda_mixed_population.B.fastq.gz \
	tests/data/lambda/lambda_mixed_population.3.fastq.gz \
	tests/data/lambda/lambda_mixed_population.4.fastq.gz \
	tests/data/lambda/lambda_mixed_population.5.fastq.gz \
	tests/data/lambda/lambda.1-2.gbk \
	tests/data/lambda/lambda.3.gbk \
	tests/data/lambda/lambda.4.gbk \
	tests/data/lambda/lambda.5.gbk

## Total core budget for parallel test execution.
## Auto-detected from make's -j flag by inspecting make's own command line via
## ps (works with all make versions, including macOS's bundled GNU Make 3.81
## which does not expose the job count via MAKEFLAGS). Falls back to 'all'
## (all available cores) when no -j was given or ps is unavailable.
## Override explicitly on the command line: 'make test TEST_CORES=16'.
## Individual tests declare how many of these cores they use via 'TEST_CORES=N'
## in their testcmd.sh (see tests/common.sh and tests/Snakefile).
TEST_CORES ?= $(or $(shell ps -p $$PPID -o args= 2>/dev/null | grep -oE -- '-j *[0-9]+|--jobs=[0-9]+' | grep -oE '[0-9]+' | head -1),all)

## Writes the TESTBINPREFIX/BRESEQ_DATA_PATH exports that tests/common.sh
## sources, at the destination path passed as $(1). Shared by all-local,
## test, test-long, and dist-hook so there's a single place to edit.
## NOTE: kept as a single printf (not multiple echo lines) because GNU
## Make treats embedded newlines from macro expansion as separate recipe
## lines / shell invocations, which breaks callers that embed this inside
## a larger shell construct (e.g. all-local's 'if ... fi' guard below).
define write_test_config
	printf 'export TESTBINPREFIX=$(top_srcdir)/src/breseq\nexport BRESEQ_DATA_PATH=$(top_srcdir)/src/share/breseq\n' > $(1);
endef

## Documentation (docs/, built with mkdocs into DOCBUILDDIR -- kept outside
## docs/ itself so generated HTML never mixes with the committed Markdown
## source). Soft/optional: skipped with a warning (not a failure) if mkdocs
## isn't installed, so contributors without it aren't blocked from building
## or from 'make dist'/'./binarydist.sh'.
DOCBUILDDIR = $(abs_top_builddir)/docs_site

## Without this, make treats the 'docs' target as satisfied by the existing
## docs/ source directory of the same name and never runs the recipe.
.PHONY: docs

docs:
if HAVE_MKDOCS
	cd $(top_srcdir) && $(MKDOCS) build --site-dir $(DOCBUILDDIR)
else
	@echo "mkdocs not found -- skipping documentation build (see dev-environment.yml)"
endif

## Main hooks called after normal make steps

## Regenerates tests/test.config after a plain 'make'/'make all' so that
## running a single test directly (./tests/run.sh, build.sh, rebuild.sh)
## doesn't require having run 'make test' first. Guarded by a writability
## check because this also runs during 'make distcheck', which builds from
## a read-only copy of the source tree -- distcheck never invokes
## 'test'/'test-long', so skipping generation there is harmless.
all-local :
	@if [ -w "$(top_srcdir)/tests" ]; then \
		$(call write_test_config,$(top_srcdir)/tests/test.config) \
	fi


install-exec-local :


clean-local :: clean-tests

	## This saves the install path for the test script
	rm -f $(abs_top_builddir)/tests/test.config
	find $(abs_top_builddir) -name '.DS_Store' -type f -delete
	rm -rf $(abs_top_builddir)/src/breseq/*.dSYM
	rm -rf $(DOCBUILDDIR)

dist-hook:
	find $(distdir) -name '.DS_Store' -type f -delete

	## Remove ._ files
	rm -rf `find $(distdir) -name ".*"`
	rm -rf ./src/breseq/gdtools.dSYM

	## Set up run_tests.sh and environmental variables
	$(call write_test_config,$(distdir)/tests/test.config)

	echo "tests/test.sh clean tests" > ${distdir}/run_tests.sh
	echo "tests/test.sh test tests" >> ${distdir}/run_tests.sh
	chmod a+x ${distdir}/run_tests.sh

	## Ship the built documentation site (not the raw Markdown source -- anyone
	## wanting to edit docs/ should clone the git repo instead).
if HAVE_MKDOCS
	$(MAKE) docs
	cp -r $(DOCBUILDDIR) $(distdir)/docs
else
	@echo "mkdocs not found -- distribution tarball will not include built documentation HTML"
endif

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
	$(call write_test_config,$(top_srcdir)/tests/test.config)
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
	$(call write_test_config,$(top_srcdir)/tests/test.config)
	START=$$(date '+%s'); \
	snakemake -s $(top_srcdir)/tests/Snakefile --cores $(TEST_CORES) --keep-going --config long=1; \
	SNAKEMAKE_STATUS=$$?; \
	END=$$(date '+%s'); \
	$(top_srcdir)/tests/print_test_summary.sh $(top_srcdir)/tests $$((END - START)); \
	SUMMARY_STATUS=$$?; \
	test $$SNAKEMAKE_STATUS -eq 0 -a $$SUMMARY_STATUS -eq 0


clean-tests:
	## 'all' cleans every test dir's output, including 'long'-named and
	## '*_disabled' ones (a disabled test can still leave output behind when
	## run by name). It is a clean-only escalation -- there is no 'make'
	## target that RUNS disabled tests.
	$(top_srcdir)/tests/test.sh clean tests all
	rm -rf $(top_srcdir)/tests/.snakemake

