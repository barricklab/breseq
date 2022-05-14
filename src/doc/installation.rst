Installation
==============

|breseq| is a command line tool implemented in C++ and R. It is compatible with a variety of UNIX-like platforms, including Linux, MacOSX, and Windows Subsystem for Linux (WSL).

The most recent |breseq| binary distributions and source code packages are available for download from `GitHub <https://github.com/barricklab/breseq/releases>`_.
The instructions in the following sections explain how to install |breseq| using these files.

.. image:: https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat
   :target: http://bioconda.github.io/recipes/breseq/README.html

**New:** Another installation option is to use the `Conda package manager <https://docs.conda.io/en/latest/>`_ to install |breseq| and all of the programs it requires. Make sure you have
`Bioconda <https://bioconda.github.io/user/install.html>`_ set up, then follow the directions for the `breseq package <http://bioconda.github.io/recipes/breseq/README.html>`_.


.. image:: images/galaxy.png
   :target: https://usegalaxy.org

**New:** If you are not comfortable with running commands in a terminal, it is also possible to install and use |breseq| on the web-based `Galaxy platform <https://usegalaxy.org>`_ (See `Installing on Galaxy`_).

Install external dependencies
++++++++++++++++++++++++++++++

|breseq| requires these software programs to be installed on your system:

* `Bowtie2 <http://bowtie-bio.sourceforge.net/bowtie2>`_ (version 2.1.0 or higher) read mapping program
* `R <http://www.r-project.org>`_ (version 2.1.4 or higher) statistical programming language

To install each dependency, visit the respective web pages linked above and follow the instructions for your platform. You must make sure that the executables for |Bowtie2| and :program:`R` are in your environment's ``$PATH`` for |breseq| to function.

.. warning::
   Certain versions of |Bowtie2|, including 2.3.1 and 2.4.5, produce SAM output that can cause |breseq| to crash. If you are using one of these versions you will get an error message when you run |breseq|.

   Output of |breseq| may vary slightly depending on your version of |Bowtie2|. When reporting |breseq| results in a publication you should always report the version of |breseq| you used and also the version of |Bowtie2|. Currently, we recommend using |Bowtie2| version 2.4.4. This is the version used by the consistency tests.

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

Installing on Windows (using WSL)
+++++++++++++++++++++++++++++++++

Download and install `Windows Subsystem for Linux (WSL) <https://docs.microsoft.com/en-us/windows/wsl/about>`_ on your machine. In the WSL terminal, you should be able to use any of the methods described above for installation. For example, you can install Conda and then use it to install the |breseq|.

Installing on Galaxy
+++++++++++++++++++++++++++++++++

If you administer a Galaxy server, |breseq| is available to install from the `Main Tool Shed <https://toolshed.g2.bx.psu.edu/>`_. See also, the directions for `Installing Tools into Galaxy <https://galaxyproject.org/admin/tools/add-tool-from-toolshed-tutorial/>`_.

If you would like to run |breseq| through the Galaxy web interface on your own computer, you can follow these steps:

1. Install a local copy of Galaxy using `planemo <https://planemo.readthedocs.io/en/latest/installation.html>`_.

2. Clone a copy of the Galaxy Toolshed (requires `git <https://git-scm.com/>`_).

.. code-block:: bash

   git clone https://github.com/galaxyproject/tools-iuc.git

3. Start the local Galaxy server

.. code-block:: bash

   cd tools-iuc/tools/breseq
   planemo serve

.. |br| raw:: html

   <br />

.. warning::
   In either case, you need to go to the settings of your Galaxy install and choose to "Whitelist" |breseq| so that it can return HTML output to the web browser. |br|

   .. image:: images/galaxy_select_whitelist.png
      :width: 300
      :target: _images/galaxy_select_whitelist.png

   .. image:: images/galaxy_select_breseq.png
      :width: 300
      :target: _images/galaxy_select_breseq.png

Troubleshooting installation
+++++++++++++++++++++++++++++++++
If you have a problem installing |breseq|, please post a detailed report as an `issue on GitHub <https://github.com/barricklab/breseq/issues>`_.
