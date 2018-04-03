Installation
==============

|breseq| is a command line tool implemented in C++ and R. It is compatible with a variety of UNIX-like platforms, including Linux, MacOSX, and Cygwin. 

The most recent |breseq| binary distributions and source code packages are available for download from `GitHub <https://github.com/barricklab/breseq/releases>`_.

Install external dependencies
++++++++++++++++++++++++++++++

|breseq| requires these software programs to be installed on your system:

* `Bowtie2 <http://bowtie-bio.sourceforge.net/bowtie2>`_ (version 2.1.0 or higher) read mapping program
* `R <http://www.r-project.org>`_ (version 2.1.4 or higher) statistical programming language 

To install each dependency, visit the respective web pages linked above and follow the instructions for your platform. You must make sure that the executables for |Bowtie2| and :program:`R` are in your environment's ``$PATH`` for |breseq| to function.

Method 1. Binary download
+++++++++++++++++++++++++++++

Linux and MacOSX packages with precompiled executables are available for download. Using these is the quickest and easiest install option that should be used by most users.

You should be able to immediately run |breseq| from within the unarchived directory structure.

For example, run this command to get breseq help:

.. code-block:: bash

   $ bin/breseq

You may also run this command to test the |breseq| install:

.. code-block:: bash

   $ ./run_tests.sh

.. NOTE::
   If you relocate the executable files in the ``bin`` directory, then you must also relocate the files in the ``share`` directory to the same location relative to the binaries (e.g., ``bin/../share/breseq``).

.. _installing-source-code-package:

Method 2. Source code download
++++++++++++++++++++++++++++++

|breseq| installation from the source code requires some basic familiarity with UNIX commands and environments.

In addition to the normal dependencies, you must also have a C++ compiler installed on your system. For example: `GCC <http://gcc.gnu.org>`_.

MacOSX does not have a C++ compiler installed by default. You must first download and install the `Apple Developer tools <http://developer.apple.com/tools/>`_ (available from the App store). Be sure that you also complete any additional steps that are necessary to install the "command-line tools".

If you have admin privileges and want to install |breseq| in a standard location accessible to all users of a computer, then see :ref:`installing-in-a-system-wide-location`. If you do not have admin privileges on your computer, then see :ref:`installing-in-the-source-directory` or :ref:`installing-in-a-custom-location`. 

The |breseq| source distribution relies on open-source code developed by others. This code is included in the |breseq| distribution under /extern:

* `SAMtools <http://samtools.sourceforge.net>`_ 

.. _installing-in-a-system-wide-location:

Installing in a system-wide location
************************************

This method requires that you have admin privileges on your machine. After installation, all users of the machine will be able to run |breseq|.

Open a terminal window and change directory to the root of the |breseq| source distribution. Then, run these commands:

.. code-block:: bash

   $ ./configure
   $ make
   $ make test
   $ sudo make install

``make test`` is optional, but recommended. It should take less than 5 minutes to run and report success at the end if everything is operating correctly.

.. _installing-in-the-source-directory:

Installing in the source directory
**********************************

This is the most robust way to compile and install |breseq| if you do not have admin privileges on a system. All of the compiled programs and libraries will be self-contained in the original source tree.

Open a terminal window and change directory to the root of the |breseq| source distribution. Then, run these commands:

.. code-block:: bash

   $ ./configure --prefix=${PWD}
   $ make
   $ make test
   $ make install

After installation, if you want to be able to call |breseq| commands without specifying the entire path to them, you will need to add the newly created "bin" directory within the |breseq| source to your $PATH.

For a :program:`bash` shell you can usually use a command like this:

.. code-block:: bash

   $ echo "export PATH=\$PATH:${PWD}/bin" >> ~/.bashrc

But the exact way to do this may depend on your system. Once you open a new terminal window so that it registers this change to your $PATH, you should be able to invoke |breseq| commands.

.. _installing-in-a-custom-location:

Installing in a custom location
*******************************

We'll assume that you've chosen to install |breseq| in ``/mnt/home/me/local``. Open a terminal window and change directory to the root of the |breseq| source distribution. Then, run these commands:

.. code-block:: bash

   $ ./configure --prefix=/mnt/home/me/local
   $ make
   $ make test
   $ make install

This will create a usual UNIX grouping of program directories (with sub-directories like ``bin``, ``lib``, ``man``, etc). 

After installation, if you want to be able to call |breseq| commands without specifying the entire path to them, you will need to add the newly created "bin" directory within the |breseq| source to your $PATH.

For a :program:`bash` shell you can usually use a command like this:

.. code-block:: bash

   $ echo "export PATH=\$PATH:/mnt/home/me/local/bin" >> ~/.bashrc

But the exact way to do this may depend on your system. You may also want to similarly update your $MANPATH, $CPPFLAGS, $LD_FLAGS, etc. Now you should be able to invoke |breseq| commands once you open a new terminal window.

Method 3. GitHub source code
+++++++++++++++++++++++++++++++++

If you are working with a development version of |breseq| cloned from the `GitHub code repository <https://github.com/barricklab/breseq>`_, then you will need to run some additional commands and have other tools installed on your system in order to get it to compile or work with the XCode project.

These requirements and commands are detailed in the DEVELOPER text file located in the main directory of the source code.

Installing on Cygwin (Windows)
+++++++++++++++++++++++++++++++++

It is possible to compile and install |breseq| and all of its dependencies in the Cygwin environment on a Windows computer. We do not currently provide a binary installer for Cygwin and are unable to help troubleshoot these installs, but here is what has worked for other users.

Before you start, use the Cygwin package manager to install these packages (which provide libraries needed to compile |breseq| and |Bowtie2|). When prompted whether to install further dependencies of a package, answer yes.

.. code-block:: bash

   R                    libncurses-devel
   gcc-core             zlib-devel
   gcc-g++              byacc
   gcc-objc++           bool
   python               pkg-config
   m4                   perl-File-Copy-Recursive
   make                 perl-Config-AutoConf
   automake             perl-ExtUtils-PkgConfig
   autoconf             mingw-pthreads
   diffutils            mingw64-x86_64-pthreads
   libiconv             mingw64-x86_64-winpthreads

Now, compile and install |Bowtie2| from source code and use the :ref:`installing-source-code-package` instructions to install |breseq|.

If the configure or make steps in either install fail, try to diagnose which dependencies are missing from the warnings and install further packages as necessary.

Troubleshooting installation
+++++++++++++++++++++++++++++++++
If you have a problem installing |breseq|, please send a detailed report to jeffrey.e.barrick@gmail.com.
