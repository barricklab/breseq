# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What is breseq

breseq is a computational pipeline for finding mutations in haploid microbial genomes by comparing short-read sequencing data to a reference sequence. It uses bowtie2 for alignment, samtools for BAM manipulation, and gnuplot for diagnostic plots. The companion tool `gdtools` manipulates Genome Diff (`.gd`) files produced by breseq.

## Build System

breseq uses GNU autotools. **A fresh checkout or worktree has no `configure` script and no
`Makefile` — you MUST run `./bootstrap.sh` then `./configure` before the first `make`, or `make`
fails with "No targets specified and no makefile found".** This bites us repeatedly; do not skip it.

```bash
# REQUIRED first build in any new checkout/worktree (generates configure + Makefile):
./bootstrap.sh                        # regenerates ./configure (also re-run after adding a .cpp source file)
./configure --prefix=$CONDA_PREFIX    # generates Makefile / config.h from configure.ac + *.am

# Normal rebuild (only works AFTER bootstrap + configure have been run once):
make

# Install (rarely wanted — see the warning below):
make install
```

Quick check before building: if `./configure` or `Makefile` is missing (`ls configure Makefile`),
you still need the bootstrap + configure step above.

**`--prefix` is not optional.** It is both where breseq installs *and* where `./configure` looks for
the zlib, miniz and htslib it links against, so it must point at the conda env holding them; it also
becomes the binaries' runtime library search path (`-rpath`), which is what lets a freshly built
`breseq` run without the env activated. A bare `./configure` can appear to work inside an activated
env — conda's compiler-activation scripts export the same paths on their own — but that is precisely
what silently broke CI when conda-forge moved to gcc 15 (`miniz/miniz.h not found at NONE/include`).
Do not rely on it. Note `make install` with this prefix installs into the *shared* conda env, so it
overwrites the `breseq`/`gdtools` every other worktree sees; building and testing in place (`make`,
`make test`) never installs anything.

**Dependencies for building and running tests** (use conda with `dev-environment.yml`):

When working in a worktree, reuse the pre-built conda env in the main repo (`breseq/env`, which is at `../../../env` relative to a worktree). Use `conda run` to invoke commands inside it without activating:
```bash
conda run -p ../../../env ./bootstrap.sh
conda run -p ../../../env sh -c './configure --prefix="$CONDA_PREFIX"'
conda run -p ../../../env make
```
(`conda run -p` sets `CONDA_PREFIX` inside the command it runs, but not in your shell — hence the
`sh -c` wrapper, which keeps the env's absolute path out of the `--prefix` you have to type.)

If that env doesn't exist yet (e.g., fresh clone or working directly in the main repo), create it first as described in `DEVELOPER`:
```bash
conda env create -f dev-environment.yml --prefix=$PWD/env
# then use: conda run -p $PWD/env <command>
```

Runtime dependencies: `bowtie2`, `gnuplot`, optionally `phylip`.

When adding a new `.cpp` source file, re-run `./bootstrap.sh` so autotools picks it up, then
`./configure --prefix=$CONDA_PREFIX && make`.

## Testing

Tests are consistency tests that run breseq/gdtools and compare output `.gd` files against expected results.

```bash
# Run all tests
make test

# Run all tests including long tests
make test-long

# Clean test outputs (always run this between testing cycles to clear stale results)
make clean-tests
```

**CRITICAL — `make test` does NOT re-run tests whose output already exists; it reports STALE results.**
Test runs are gated by per-test done-files (`tests/<name>/test.done`, `test.result`, and Snakemake's
`.snakemake/`). If those are present from a previous run, `make test` skips the test and re-prints its
**cached** PASS/FAIL — it does **not** rebuild `data/annotated.gd` or re-compare it. So after **any**
code change (or `git` operation that changes the binary), running `make test` alone can report "all
passed" even though the affected tests never actually ran against your new build. This has caused real
regressions to be merged. Rules:
- **Always run `make clean-tests` immediately before `make test`** after changing any `src/` file or
  rebuilding — every time, not just "between cycles". For a single test, run
  `./tests/test.sh clean <name>` (or `make clean-tests`) before `./tests/test.sh test <name>`.
  (`./tests/test.sh rebuild <name>` cleans + re-runs on its own; a bare `test` action does not.)
- **Sanity-check the timing.** `tests/print_test_summary.sh` prints a per-test time and an overall
  total. A real full-suite run takes on the order of minutes; if the summary shows a near-zero total
  (e.g. `total time 00:00:01` for 70+ tests), **nothing actually re-ran** — the results are stale.
  Treat a suspiciously fast run as a failed run: `make clean-tests` and run again.
- Do not commit/merge/push based on a `make test` result you did not first force to re-run this way.

**To run a single test, use `tests/test.sh <action> <name>`** — this is the canonical single-test
entry point (streams output live, no Snakemake). Prefer it over the thin `run.sh`/`rebuild.sh`/
`build.sh` wrappers, which all just dispatch to the same `tests/<name>/testcmd.sh <action>`:

```bash
# Run one test:
./tests/test.sh test lambda_mixed_pop

# Rebuild its expected output after intentional changes (re-runs the test):
./tests/test.sh rebuild lambda_mixed_pop

# Promote a failed run's existing output to be the new expected output (no re-run):
./tests/test.sh build lambda_mixed_pop

# Clean just one test's output:
./tests/test.sh clean lambda_mixed_pop
```

(Use `tests` as the name — e.g. `./tests/test.sh test tests` — to batch every discovered test
serially; `make test` is preferred for the full suite since it runs them in parallel.)

A plain `make`/`make all` writes `tests/test.config` (sets `TESTBINPREFIX`, `BRESEQ_DATA_PATH`) as
part of `all-local`, so running a single test directly (`./tests/test.sh <action> <name>`)
never requires a prior `make test` (but it does require a successful `make` first — see Build System). `make test`/`make test-long` also (re)write it before running
all discovered tests **in parallel** via Snakemake (`tests/Snakefile`), bounded by a total core
budget. The shared test infra (the `do_build`/`do_check`/`do_clean`/... conventions that every
`testcmd.sh` dispatches into) lives in `tests/common.sh`; `Snakefile` only discovers and launches
`testcmd.sh test` for each test, it doesn't change those conventions.

- **Core budget**: control the total number of cores used across all parallel tests with
  `make -j8 test` (auto-detected) or `make test TEST_CORES=16` (explicit override; default: all
  available cores when no `-j` is given). Each test declares how many of those cores *it* uses via
  a `TEST_CORES=N` line near the top of its `testcmd.sh`, set *before* sourcing `common.sh` (see
  `tests/lambda_mixed_pop/testcmd.sh`); tests that don't set it (typically `gdtools`-only tests)
  default to `TEST_CORES=1`. `common.sh` derives `BRESEQ_TEST_THREAD_ARG="-j ${TEST_CORES}"` from
  this, so the same value drives both scheduling and the actual thread count passed to `breseq`.
- **Ordering dependent tests**: a test that needs another test's output to exist first (e.g.
  `lambda_mixed_pop_cl_tabulate` reads `data/reference.{bam,fasta,gff3}` generated by
  `lambda_mixed_pop`) declares this with `TEST_DEPENDS=<name>[,<name>...]` next to `TEST_CORES`
  in its `testcmd.sh`. `Snakefile` reads this to add a real dependency edge so the prerequisite
  test finishes successfully first — `common.sh` ignores the variable, and the serial
  `./tests/test.sh <action> tests` runner (and its `run.sh`/`build.sh`/`rebuild.sh` wrappers)
  doesn't need it since their sorted discovery order already happens to run dependencies first.
- **External test data (long tests only)**: a `long`-named test that needs large real-world data
  declares `TEST_DATA=<dataset>[,<dataset>...]` next to `TEST_CORES` in its `testcmd.sh`. Datasets
  are defined in `tests/test_data_manifest.tsv` (`dataset` / `file` / `md5` / `url`, tab-separated)
  and fetched by `tests/fetch_test_data.sh`, which downloads with `curl` to a temp file in the
  destination directory, md5-verifies it, and only then renames it into place — so a failed or
  corrupt download can never leave a plausible-looking file behind. Reference files are used as the
  gzipped `.gbff.gz`/`.fna.gz` that NCBI serves (breseq reads gzipped references transparently), so
  the md5 covers exactly the bytes breseq parses.
  - **A checksum mismatch is a hard failure, never a warning.** It means the remote data drifted, so
    every `expected.gd` built from it is suspect. Recovery: re-query the source for the new md5
    (the manifest header has the exact `curl` commands), update it, `./tests/test.sh rebuild
    <name>`, and *review the resulting diff* — don't just paste in the new checksum.
  - **Cache**: `tests/data/downloads/` (gitignored), or set `BRESEQ_TEST_DATA_DIR` to share one
    cache between worktrees. It lives under `tests/data/` but deliberately outside any test
    directory, because `do_clean` does `rm -Rf <testdir>/data`. Budget roughly 0.5–2 GB per dataset
    *plus* a similar amount of `tests/<name>/data/` output while a long test runs.
  - **`make test` never touches the network.** Long tests aren't discovered without
    `--config long=1`, and `Snakefile` raises a `WorkflowError` if a non-long test declares
    `TEST_DATA`. `Snakefile` parses `TEST_DATA` with the same anchored regex it uses for
    `TEST_CORES`/`TEST_DEPENDS` and makes each dataset a job-graph node, so a dataset shared by
    several long tests downloads once; `common.sh`'s `do_test` runs the same fetch (a stat-only
    no-op when already cached) so `./tests/test.sh test <name>`, which bypasses Snakemake, works too.
  - **Targets**: `make test-data` pre-fetches everything, `make verify-test-data` forces a full
    re-hash of the cache (to detect local corruption — it re-hashes, it does not re-download), and
    `make clean-test-data` deletes it. `make clean-tests` and `make clean` deliberately do **not**
    delete downloads; they only remove `tests/.test_data_stamps/`.
  - When rebuilding a long test's `expected.gd`, do it with `BRESEQ_TEST_DATA_DIR` **unset** so the
    committed file keeps repo-relative paths. (`DIFF_IGNORE` strips `#=COMMAND`/`#=REFSEQ`/
    `#=READSEQ` before comparing, so an absolute path doesn't *break* the test — it just bakes
    someone's home directory into a committed file.) Note `#=MAPPED-BASES`/`#=MAPPED-READS` *are*
    compared, so a bowtie2 version bump fails a long test even with breseq unchanged.
  - `tests/fetch_test_data_1` is an offline self-test of the fetcher (it uses `curl`'s `file://`
    support against small committed files), so the download, checksum-mismatch and failure paths are
    covered by every `make test` run without any network.
  - Beware: test-name matching for "long" is a **substring** check, so a directory named e.g.
    `nonlong_foo` is treated as a long test and skipped by `make test`.
- **Logs**: each test's output is captured to `tests/<test_name>/test.log`; `make clean-tests`
  removes these (along with `test.result`, `test.done`, Snakemake's `.snakemake/` directory, and
  `tests/.test_data_stamps/`).
- **Summary & CI**: after the run, `tests/print_test_summary.sh` prints a PASS/FAIL + timing table
  per test plus an overall total, and exits non-zero if any test failed — `make test`/`make
  test-long` propagate that status, so they're suitable for driving CI (e.g. GitHub Actions).

Test directories prefixed with `_` are skipped. Directories with `long` in their name are only run with `make test-long`.

A test directory whose name ends in `_disabled` is skipped by both `make test` and `make test-long` — it runs *only* when invoked explicitly by name (e.g. `./tests/test.sh test <name>_disabled`, `./tests/test.sh rebuild <name>_disabled`). There is deliberately no `make` target that runs disabled tests in bulk (the serial `all` runners skip them too). Use this to park a temporarily-broken test in the tree without failing CI. `make clean-tests` still cleans disabled dirs' output.

A test directory whose name matches `tests/*_local*/` (e.g. `tests/<name>_local`) is gitignored and understood to be an intentionally untracked, personal/in-progress test — it's still discovered and run normally by `make test`/`make test-long`, since test discovery is filesystem-based (`tests/Snakefile`'s glob over `tests/*/testcmd.sh`, and equivalently `tests/test.sh`) and doesn't care about git tracking status. The suffixes compose: a `tests/<name>_local_disabled` directory is both untracked (matches `*_local*`) and disabled (ends in `_disabled`), i.e. an uncommitted test that only runs when named explicitly. The matching gitignore patterns for untracked data are `tests/*_local*/` and `tests/data/*_local*/`.

**Creating a new test:**
1. Add data files under `tests/data/` (reuse existing data when possible; don't add large files to git).
2. Create `tests/<test_name>/testcmd.sh` by copying an existing one and editing the `TESTCMD=` and `CURRENT_OUTPUTS`/`EXPECTED_OUTPUTS` lines. If the test runs `breseq` with multiple threads, set `TEST_CORES=N` (matching the `-j` value) before sourcing `common.sh`; otherwise leave it unset (defaults to 1).
3. Run `./tests/<test_name>/testcmd.sh rebuild` to generate `expected.gd`.
4. For a test on large real-world data, name the directory `long_<something>`, add its files to
   `tests/test_data_manifest.tsv` (with checksums taken from the source — see the manifest header),
   and declare `TEST_DATA=<dataset>[,...]` in the `testcmd.sh`. Reference it as
   `${DOWNLOADDIR}/<dataset>/<file>`. It then runs only under `make test-long`.

## Code Architecture

### Two binaries, one library

- `src/breseq/breseq` — the main pipeline binary
- `src/breseq/gdtools` — genome diff manipulation tool
- `src/breseq/libbreseq.la` — static library shared by both

Headers and implementation files both live directly in `src/breseq/`.

Entry points: `breseq_cmdline.cpp` (dispatches breseq subcommands) and `gdtools_cmdline.cpp` (dispatches gdtools subcommands).

### breseq pipeline stages (in order)

The main `breseq` default action runs these stages sequentially, using done-files to skip completed steps on restart:

1. **Read and reference sequence file input** — parse FASTQ/FASTA inputs, convert reference GenBank → FASTA
2. **Read alignment to reference genome** — runs `bowtie2-build` + `bowtie2` (optionally 2-stage: stringent then relaxed)
3. **Preprocessing alignments for candidate junction identification**
4. **Preliminary analysis of coverage distribution** — `coverage_distribution.cpp`
5. **Identifying junction candidates** — `candidate_junctions.cpp`
6. **Re-alignment to junction candidates** — second bowtie2 alignment pass
7. **Resolving best read alignments** — `resolve_alignments.cpp`
8. **Creating BAM files** — samtools sort/index
9. **Tabulating error counts** — `error_count.cpp`
10. **Re-calibrating base error rates** — `error_count.cpp`
11. **Examining read alignment evidence** — `identify_mutations.cpp`
12. **Polymorphism statistics** — `reference_sequence.cpp`
13. **Predicting copy number variation**
14. **Output** — mutation prediction, annotation, HTML report (`output.cpp`)

### Intermediate files and cleanup

Each stage writes a done-file when it completes; on restart, `Settings::do_step()` skips a stage whose
done-file already exists. Large intermediate files (BAMs, per-read sidecars, etc.) are deleted once they
are no longer needed — unless the user passes `-k`/`--keep-intermediates` (`settings.keep_all_intermediates`).

**When you add a new file that only some stages need, register it for cleanup** with
`settings.track_intermediate_file(done_key, file_path)` (`settings.cpp`). The `done_key` is the
**done-file of the last stage that consumes the file** — cleanup happens inside `Settings::done_step(done_key)`,
which deletes every file tracked under that key (respecting `keep_all_intermediates`). A file may be
*written* during an earlier stage than the one it is keyed to; e.g. the candidate-junction split BAMs
(`#.split.bam`) and their paired-mapping sidecars (`#.split_pair_positions.csv`) are written during
preprocessing (stage 3) but keyed to `candidate_junction_done_file_name` because stage 5
(`identify_candidate_junctions`) is the last consumer. The tracking map is serialized into the done-file
and restored on restart, so cleanup still happens even when the producing stage was skipped. Match this
pattern (same `done_key`, same `do_preprocess`/gating guards) whenever you emit a companion file next to
an existing tracked one, so it shares that file's lifecycle instead of leaking.

**Re-running a single pipeline stage for verification** — the trick of "delete a stage's done-file and
re-run breseq so only that stage re-executes" **only works if the original run used `-k`/`--keep-intermediates`**.
Without `-k`, the intermediate files that stage *consumes* have already been deleted by cleanup, so the
re-run crashes (e.g. "Could not open file for reading ..."). Normal test runs (`make test`,
`tests/test.sh`) do **not** pass `-k`, so their output dirs can't be poked this way. To iterate on, say,
the Output/HTML stage: run breseq once with `-k` into a scratch `-o` dir, then (with intermediates kept)
edit an intermediate `.gd` and remove just that stage's done-file to re-render. The last-stage
`data/annotated.gd` may linger even without `-k`, but the upstream files the stage needs will not.

**Also remove the top-level `output/output.done`.** Before any per-stage `do_step` check, breseq
short-circuits the *entire* pipeline if `settings.output_done_file_name` (`output/output.done`) exists
(`breseq_cmdline.cpp` ~line 1437: prints only "ALREADY COMPLETE Output" and returns). So deleting a
mid-pipeline stage's done-file alone does nothing — the run exits immediately. To re-run e.g. stage 08,
remove **both** `output/output.done` and `08_mutation_identification/mutation_identification.done`;
breseq then reports "ALREADY COMPLETE" for stages 01–07 (cached) and re-runs 08 onward.

**Delete only the `*.done` file, never the whole stage directory.** Do NOT `rm -rf` a numbered stage
directory (`0*_.../`) to "reset" a step: some output artifacts are produced once by an early stage and
are not regenerated when a *later* stage re-runs — e.g. the paired-mapping-distance plot
(`evidence/#.pair_distance.svg`) is drawn only in stage 3, so deleting its directory (or the whole
`evidence/` dir) leaves the Output step unable to recreate it, and the `summary.html` link silently
drops. Deleting only the done-file preserves the on-disk intermediates that skipped-but-upstream stages
already wrote.

### Genome Diff (`.gd`) format

The central data structure is `cGenomeDiff` (`genome_diff.h`/`genome_diff.cpp`) containing a list of `cDiffEntry` objects (`genome_diff_entry.h`/`genome_diff_entry.cpp`).

`cDiffEntry` is a string→string map. Entries are tab-delimited in files; columns after the fixed positional fields are `key=value` pairs.

Mutation types (enum `gd_entry_type`): `SNP`, `SUB`, `DEL`, `INS`, `MOB`, `AMP`, `INV`, `CON`, `INT`

Evidence types: `RA` (read alignment), `MC` (missing coverage), `JC` (new junction), `CN` (copy
number), `UN` (unknown), `SC` (soft clipping), `DP` (discordant pair), `MP` (missing pair),
`PD` (pair distance)

The last four are experimental and opt-in (`--predict-soft-clipping`, `--predict-discordant-pairs`,
`--predict-missing-pairs`, `--predict-pair-distance`). The three pair-based types divide the space by
what is anomalous about a read pair: `DP` fires on pairs that are individually discordant (wrong
orientation, or distance past the cutoff), `MP` on reads whose mate did not map anywhere, and `PD` on
pairs that are individually unremarkable but collectively shifted in mapping distance — the deletions
and insertions of a few hundred bases the other two cannot see. Where a PD and a DP describe the same
breakpoint the DP is removed, since PD uses the whole pair population where DP uses only its tail.

`MP` counts only reads whose mate produced **no alignment at all** — not merely reads flagged
`BAM_FMUNMAP`, which also covers mates the aligner placed and breseq then rejected on
`--require-match-fraction` (see `kBreseqMateNeverAlignedBAMTag`). `SC`, `MP` and `PD` are each
decided by a genome-wide `score` (`--soft-clipping-score-cutoff`, `--missing-pair-score-cutoff`,
`--pair-distance-score-cutoff`, all defaulting to 3) computed against a null **fitted to the run**,
not against a fixed cutoff; each reports that null in a gates table in `summary.html`. A fixed
frequency cutoff cannot substitute: it implicitly assumes a zero background, which holds in a
simulation and in nothing else.

For `SC` the score is not what does most of the work on real data — the **strand** test is
(`--soft-clipping-fisher-strand-p-value-cutoff`, reject `FISHER_STRAND`). The dominant SC false
positive is a dark-cycle poly-G read tail, which is always the read's 3' end and therefore appears
on exactly one strand for a given clip direction, whereas reads clipped at a real breakpoint come
from both. Measured over 29 LTEE clones, 95% of accepted SC calls had every clipped read on one
strand. The tail-consensus test cannot see this at all, because poly-G tails agree with each other
perfectly; the companion `LOW_COMPLEXITY_TAIL` gate
(`--soft-clipping-maximum-tail-homopolymer-fraction`) catches the low-count positions where the
strand test has no power. The SC gates table in `summary.html` reports what fraction of the run's
clip events were one-strand — a value near 100% means the run's clipping is artifact, and it is the
first thing to check when a library predicts implausibly many SC items.

Validation/annotation types: `CURA`, `FPOS`, `PHYL`, `TSEQ`, `PFLP`, `RFLP`, `PFGE`, `NOTE`, `MASK`

### Key modules

- `settings.h/cpp` — `Settings` class holds all configuration; `cReadFiles` manages input read files; static `Settings::pool` is the thread pool
- `reference_sequence.h/cpp` — `cReferenceSequences`/`cAnnotatedSequence` hold parsed reference genome data
- `output.h/cpp` — generates HTML output and handles zip (uses bundled miniz for `.zip` output)
- `alignment_output.h/cpp` — renders HTML alignment views (used by `BAM2ALN` subcommand)
- `coverage_distribution.h/cpp` — coverage analysis
- `mutation_predictor.h/cpp` — final mutation prediction logic
- `anyoption.h/cpp` — command-line option parsing (custom library)

### Plotting

All statistics and plotting are native C++; no R dependency remains. Diagnostic plots (coverage distribution, coverage overview, error rates, JC precision/sensitivity) are drawn by writing a gnuplot script to a temp file and invoking `gnuplot` as a subprocess via the `run_gnuplot_script()` helper in `common.h`.

## Output Structure

breseq writes to `-o <output_dir>` (default: current directory):
- `output/index.html` — main HTML report
- `output/annotated.gd` — final Genome Diff with mutations
- `data/reference.bam` — read alignments
- `data/reference.fasta` — reference sequences

## gdtools Subcommands

`VALIDATE`, `APPLY`, `ANNOTATE`/`COMPARE`, `MUTATIONS`, `CHECK`, `NORMALIZE`, `SUBTRACT`, `INTERSECT`, `UNION`/`MERGE`, `FILTER`/`REMOVE`, `MASK`, `NOT-EVIDENCE`, `GD2VCF`, `VCF2GD`, `GD2GVF`, `GD2CIRCOS`, `MUMMER2MASK`, `COUNT`, `PHYLOGENY`

## Claude Code Workflow Rules

- **Always ask before committing, merging, or pushing** — never do these automatically.
- **Work in the current worktree** — all edits and builds happen here. Do not modify files in the main repository except:
  1. When running inside a git worktree, use the shared conda env at `../../../env` (i.e., `breseq/env` in the main repo) via `conda run -p ../../../env <command>`. If that env is absent, fall back to creating one locally as described in `DEVELOPER`.
  2. When explicitly told to, you may merge changes back to the `master` branch in the main repo.
- **To land a worktree branch on `master`, push it into the repo from the worktree** — a session running
  in a worktree is refused any git command that redirects at the shared checkout (`git -C
  /Users/jbarrick/src/breseq ...`, `--git-dir`, `cd` into it), and that refusal applies to `!`-prefixed
  commands too, so `git -C <main repo> merge <branch>` is not available. This is:

  ```bash
  git push . HEAD:master          # from the worktree root, no -C
  ```

  It moves the `master` ref *and* updates the main checkout's working files in one step, because the
  repo sets `receive.denyCurrentBranch = updateInstead` (in the shared `.git/config`; note
  `extensions.worktreeConfig` is on, so never set this with `--worktree` or the main checkout won't see
  it). Without that setting the push is refused outright, which is the safe default — the dangerous
  values are `warn`/`ignore`, which move the ref while leaving the main checkout's files stale so every
  merged change shows up there as an uncommitted deletion.

  Two limits: it only fast-forwards, so a branch that has diverged still needs a real `git merge` run
  from a shell in the main checkout; and it refuses (`Working directory has unstaged changes`) if the
  main checkout is dirty, leaving the ref untouched. Both failures are safe — nothing moves.
- **NEVER edit files under `/Users/jbarrick/src/breseq/src/...` when working in a worktree — that is the MAIN repo, not your worktree.** This mistake has happened repeatedly and is silent: the edit "succeeds," but you built and tested the *unmodified* worktree copy, so nothing changes and the diff hides in the main repo. Concrete rules to prevent it:
  - Your worktree root is the `Primary working directory` in the environment block (e.g. `.../breseq/.claude/worktrees/<name>/`). **Every file you Read, Edit, or Write must have that worktree root as a prefix.** A path that starts with `/Users/jbarrick/src/breseq/src/` or `/Users/jbarrick/src/breseq/tests/` (i.e. the repo root *without* the `.claude/worktrees/<name>/` segment) is the wrong copy — do not edit it.
  - **Explore/Plan/general-purpose subagents report absolute paths in the MAIN repo** (they resolve symlinks / run from the repo root). Treat every `file:line` a subagent returns as main-repo-relative and **translate it to the worktree path before editing** — strip the `/Users/jbarrick/src/breseq/` prefix and re-root it under your worktree directory. The line numbers still apply; only the path prefix changes.
  - Before your first Edit in a session, sanity-check with `git rev-parse --show-toplevel` (or confirm the path against the `Primary working directory`) so you know which tree you're touching.
  - If you discover edits landed in the main repo, re-apply them under the worktree, then get them
    reverted in the main repo. A worktree session cannot do that itself — `git -C
    /Users/jbarrick/src/breseq checkout -- <files>` is refused by the same isolation guard described
    above, `!`-prefixed included — so name the files and ask for `git checkout -- <files>` to be run
    from a shell in the main checkout.
- **Run `make test` as part of the development cycle** — run the full test suite before considering a feature or fix complete.
- **Run tests from the worktree root** — `./tests/test.sh <action> <name>` (single test) and `make test` (full suite) must be invoked from the top of the worktree (not from a subdirectory). `make` auto-generates this worktree's `tests/test.config` as part of `all-local` and points `TESTBINPREFIX` at this worktree's `src/breseq` — you don't create or edit it by hand. (It also exports `BRESEQ_DATA_PATH=./src/share/breseq`, but that export is redundant for in-tree runs: a binary built at `src/breseq/` already finds its data at the sibling `src/share/breseq/` on its own — see the next bullet — so `BRESEQ_DATA_PATH` only actually matters if the executable is moved away from that sibling.)
- **A freshly built `breseq`/`gdtools` runs standalone with no env vars needed** — the runtime data files (`breseq_icon.png`, `jszip.min.js`, etc.) live in `src/share/breseq/`, a sibling of `src/breseq/` where the binaries are built, matching the relative `bin/`↔`share/<pkg>/` layout autotools already computes between `--bindir`/`--datadir`. `BRESEQ_DATA_PATH` is only needed to override this (e.g. a relocated/custom install).

## Distribution

```bash
# Source archive
make distcheck

# macOS universal binary
./binarydist.sh
```
