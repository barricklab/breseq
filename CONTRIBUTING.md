# Contributing to breseq

Thanks for your interest in contributing to _breseq_! This document covers how to
get involved in the project. For build/test/distribution instructions, see the
[`DEVELOPER`](DEVELOPER) file. Please also read our
[Code of Conduct](CODE_OF_CONDUCT.md) before participating.

## Contributing code

We welcome your contributions! You can check out suggested
[coding project issues](https://github.com/barricklab/breseq/labels/coding-project)
or [look for bugs to fix](https://github.com/barricklab/breseq/issues). Please
comment on an issue if you are working on it so that we can avoid duplicating
efforts and provide guidance.

To contribute your code, you should:

1. [Fork the _breseq_ repository](https://help.github.com/en/articles/fork-a-repo).
2. Work on your code, create tests for the code, have others try it out...
3. [Submit a pull request to the _breseq_ repository](https://help.github.com/en/articles/creating-a-pull-request-from-a-fork), so that it can be integrated into the project.

## Developing breseq extensions

If you are creating a Python or R package that operates on the output of
_breseq_ runs, these tools may require the installation of other programs
(e.g., Circos, de novo assemblers).

1. Install _breseq_ as described in [`DEVELOPER`](DEVELOPER), or download the
   most recent [release version](https://github.com/barricklab/breseq/releases)
   and follow the [installation instructions](https://breseq.barricklab.org/installation/).

### Developing Python modules

If you are developing code that reads in `genomediff` files output by
_breseq_, we recommend that you use the official
[`genomediff` Python module](https://github.com/barricklab/genomediff-python)
for this purpose rather than re-inventing the wheel. If you need to update the
functionality of the `genomediff` Python parser, please contribute to that
GitHub project.

## Developing breseq core code

If you are coding in C++ to change the core `breseq` and/or `gdtools`
commands:

1. Clone the [GitHub repository](https://github.com/barricklab/breseq). (Do
   not download a release version.)
2. Follow the "Developing with Command-Line Compilation" instructions in
   [`DEVELOPER`](DEVELOPER).

## All developers

You should familiarize yourself with how _breseq_ works and its output files
by completing these activities from the
[_breseq_ documentation](https://breseq.barricklab.org/):

1. [Test Drive](https://breseq.barricklab.org/test-drive/)
2. Tutorials: [Clonal Samples](https://breseq.barricklab.org/tutorial-clones/) | [Mixed Populations](https://breseq.barricklab.org/tutorial-populations/)

## Coding style

Our philosophy is that the _breseq_ core executables should have minimal
external requirements and maximum compatibility for installation and
compilation on a wide range of systems.

In order to ensure this is the case and for maintainability of the code,
please follow these coding conventions:

* C++ code
  * In general, follow the [Google C++ Style Guide](https://google.github.io/styleguide/cppguide.html).
  * Only utilize features available in the **C++11** standard.
  * Do not link to any external convenience libraries (e.g., Boost).
* R code
  * In general, follow the [Google R Style Guide](https://google.github.io/styleguide/Rguide.html).
  * Use only core libraries/functionality included with a base R install.

If you are writing a _breseq_ extension in Python or R, you can utilize
whatever other modules are needed for functionality. In these cases, we
recommend releasing your code via [PyPI](https://pypi.org/) or
[CRAN](https://cran.r-project.org/) to manage requirements.
