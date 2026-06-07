# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What is breseq

breseq is a computational pipeline for finding mutations in haploid microbial genomes by comparing short-read sequencing data to a reference sequence. It uses bowtie2 for alignment, samtools for BAM manipulation, and R for statistical steps. The companion tool `gdtools` manipulates Genome Diff (`.gd`) files produced by breseq.

## Build System

breseq uses GNU autotools. The build flow is:

```bash
# First time or after adding source files:
./bootstrap.sh
./configure

# Normal rebuild:
make

# Install:
make install
```

**Dependencies for building** (use conda with `dev-environment.yml`):
```bash
conda env create -f dev-environment.yml
conda activate breseq-dev
# Set env vars as shown in DEVELOPER file, then re-activate
./bootstrap.sh && ./configure && make
```

**Dependencies for running tests** (separate env — the run env installs a compiler that conflicts with building):
```bash
conda env create -f run-environment.yml
conda activate breseq-run
```

Runtime dependencies: `bowtie2`, `R`/`Rscript`, optionally `phylip`.

When adding a new `.cpp` source file, re-run `./bootstrap.sh` so autotools picks it up, then `./configure && make`.

## Testing

Tests are consistency tests that run breseq/gdtools and compare output `.gd` files against expected results.

```bash
# Run all tests
make test

# Run all tests including long tests
make test-long

# Run a single test (or 'all' for every test) directly, without Snakemake
./tests/run.sh lambda_mixed_pop

# Clean test outputs
make clean-tests

# Rebuild expected output after intentional changes (re-runs the test; 'all' rebuilds everything)
./tests/rebuild.sh lambda_mixed_pop

# Promote a failed run's existing output to be the new expected output (no re-run)
./tests/build.sh lambda_mixed_pop
```

`make test`/`make test-long` write `tests/test.config` (sets `TESTBINPREFIX`, `BRESEQ_DATA_PATH`)
and then run all discovered tests **in parallel** via Snakemake (`tests/Snakefile`), bounded by a
total core budget. The shared test infra (the `do_build`/`do_check`/`do_clean`/... conventions
that every `testcmd.sh` dispatches into) lives in `tests/common.sh`; `Snakefile` only discovers
and launches `testcmd.sh test` for each test, it doesn't change those conventions.

- **Core budget**: control the total number of cores used across all parallel tests with
  `make test TEST_CORES=16` (default: all available cores). Each test declares how many of those cores *it* uses via
  a `TEST_CORES=N` line near the top of its `testcmd.sh`, set *before* sourcing `common.sh` (see
  `tests/lambda_mixed_pop/testcmd.sh`); tests that don't set it (typically `gdtools`-only tests)
  default to `TEST_CORES=1`. `common.sh` derives `BRESEQ_TEST_THREAD_ARG="-j ${TEST_CORES}"` from
  this, so the same value drives both scheduling and the actual thread count passed to `breseq`.
- **Logs**: each test's output is captured to `tests/<test_name>/test.log`; `make clean-tests`
  removes these (along with `test.result`, `.test_done`, and Snakemake's `.snakemake/` directory).
- **Summary & CI**: after the run, `tests/print_test_summary.sh` prints a PASS/FAIL + timing table
  per test plus an overall total, and exits non-zero if any test failed — `make test`/`make
  test-long` propagate that status, so they're suitable for driving CI (e.g. GitHub Actions).

Test directories prefixed with `_` are skipped. Directories with `long` in their name are only run with `make test-long`.

**Creating a new test:**
1. Add data files under `tests/data/` (reuse existing data when possible; don't add large files to git).
2. Create `tests/<test_name>/testcmd.sh` by copying an existing one and editing the `TESTCMD=` and `CURRENT_OUTPUTS`/`EXPECTED_OUTPUTS` lines. If the test runs `breseq` with multiple threads, set `TEST_CORES=N` (matching the `-j` value) before sourcing `common.sh`; otherwise leave it unset (defaults to 1).
3. Run `./tests/<test_name>/testcmd.sh rebuild` to generate `expected.gd`.

## Code Architecture

### Two binaries, one library

- `src/c/breseq/breseq` — the main pipeline binary
- `src/c/breseq/gdtools` — genome diff manipulation tool
- `src/c/breseq/libbreseq.la` — static library shared by both

All headers live in `src/c/breseq/libbreseq/`. Implementation files are in `src/c/breseq/`.

Entry points: `breseq_cmdline.cpp` (dispatches breseq subcommands) and `gdtools_cmdline.cpp` (dispatches gdtools subcommands).

### breseq pipeline stages (in order)

The main `breseq` default action runs these stages sequentially, using done-files to skip completed steps on restart:

1. **Read and reference sequence file input** — parse FASTQ/FASTA inputs, convert reference GenBank → FASTA
2. **Read alignment to reference genome** — runs `bowtie2-build` + `bowtie2` (optionally 2-stage: stringent then relaxed)
3. **Preprocessing alignments for candidate junction identification**
4. **Preliminary analysis of coverage distribution** — uses R scripts
5. **Identifying junction candidates** — `candidate_junctions.cpp`
6. **Re-alignment to junction candidates** — second bowtie2 alignment pass
7. **Resolving best read alignments** — `resolve_alignments.cpp`
8. **Creating BAM files** — samtools sort/index
9. **Tabulating error counts** — `error_count.cpp`
10. **Re-calibrating base error rates** — R scripts
11. **Examining read alignment evidence** — `identify_mutations.cpp`
12. **Polymorphism statistics** — R scripts
13. **Predicting copy number variation**
14. **Output** — mutation prediction, annotation, HTML report (`output.cpp`)

### Genome Diff (`.gd`) format

The central data structure is `cGenomeDiff` (`genome_diff.h`/`genome_diff.cpp`) containing a list of `cDiffEntry` objects (`genome_diff_entry.h`/`genome_diff_entry.cpp`).

`cDiffEntry` is a string→string map. Entries are tab-delimited in files; columns after the fixed positional fields are `key=value` pairs.

Mutation types (enum `gd_entry_type`): `SNP`, `SUB`, `DEL`, `INS`, `MOB`, `AMP`, `INV`, `CON`, `INT`

Evidence types: `RA` (read alignment), `MC` (missing coverage), `JC` (junction candidate), `CN` (copy number), `UN` (unmatched)

Validation/annotation types: `CURA`, `FPOS`, `PHYL`, `TSEQ`, `PFLP`, `RFLP`, `PFGE`, `NOTE`, `MASK`

### Key modules

- `settings.h/cpp` — `Settings` class holds all configuration; `cReadFiles` manages input read files; static `Settings::pool` is the thread pool
- `reference_sequence.h/cpp` — `cReferenceSequences`/`cAnnotatedSequence` hold parsed reference genome data
- `output.h/cpp` — generates HTML output and handles zip (uses bundled miniz for `.zip` output)
- `alignment_output.h/cpp` — renders HTML alignment views (used by `BAM2ALN` subcommand)
- `coverage_distribution.h/cpp` — coverage analysis
- `mutation_predictor.h/cpp` — final mutation prediction logic
- `anyoption.h/cpp` — command-line option parsing (custom library)

### Bundled dependencies

- `extern/samtools-1.3.1` — modified samtools used as a library (see `extern/samtools_modifications.txt`)
- `miniz*.{c,h}` — miniz zip library (split into multiple files; compiled with `-DMINIZ_NO_ZLIB_COMPATIBLE_NAMES` to avoid conflicts with system libz)

### R scripts

R scripts embedded in the binary's data directory: `coverage_distribution.r`, `plot_coverage.r`, `plot_error_rate.r`, `plot_jc_scores.r`, `polymorphism_statistics.r`. These are invoked as subprocesses via `Rscript`.

## Output Structure

breseq writes to `-o <output_dir>` (default: current directory):
- `output/index.html` — main HTML report
- `output/annotated.gd` — final Genome Diff with mutations
- `data/reference.bam` — read alignments
- `data/reference.fasta` — reference sequences

## gdtools Subcommands

`VALIDATE`, `APPLY`, `ANNOTATE`/`COMPARE`, `MUTATIONS`, `CHECK`, `NORMALIZE`, `SUBTRACT`, `INTERSECT`, `UNION`/`MERGE`, `FILTER`/`REMOVE`, `MASK`, `NOT-EVIDENCE`, `GD2VCF`, `VCF2GD`, `GD2GVF`, `GD2CIRCOS`, `MUMMER2MASK`, `COUNT`, `PHYLOGENY`

## Distribution

```bash
# Source archive
make distcheck

# macOS universal binary
./binarydist.sh
```
