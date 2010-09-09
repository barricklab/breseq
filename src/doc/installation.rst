Installation
==============

:program:`breseq` is a command line tool implemented in Perl/C++/R. It should compile and function on a variety of Unix platforms, including MacOSX.

1. Install external dependencies
---------------------------------

:program:`breseq` requires several external packages and software programs to be installed.

* GCC (C++ compiler) already installed on many systems
* Perl (version 5.8 or higher) already installed on many systems
* R Statistical Programming Language (version 2.10 or higher)
* Boost

To install each missing dependency, visit the respective web pages linked above and follow the instructions for your platform.

MacOSX Instructions
********************

You must have administrator privileges to install :program:`breseq`. We recommend that you first install the package manager `MacPorts <http://www.macports.org/>`_.

* :program:`GCC`: Download the :program:`Developer tools` package from http://developer.apple.com/tools/ to install GCC. 
* :program:`Perl`: is already installed on MacOSX systems. 
* :program:`R`: install with :program:`MacPorts` command: ``%% sudo port install R``
* :program:`Boost`: install with :program:`MacPorts` command: ``%% sudo port install boost``

2. Build :program:`breseq`
----------------------------

Before compiling :program:`breseq`, move the entire archive that you downloaded to a stable location on your system. If you later move this directory, then you may need to re-compile :program:`breseq` to have it operate correctly.

Open a terminal window, navigate to the root of the source distribution and run these commands::

  %% ./configure
  %% make
  %% make install
  %% make test

These commands compile and install not only :program:`breseq`, they also install several tools developed by others that are included in the source package for ease of installation.

They include:

* `BioPerl <http://www.bioperl.org>`_ (*partial* install)
* `SAMtools <http://samtools.sourceforge.net>`_ 
* `Bio::DB::Sam <http://search.cpan.org/~lds/Bio-SamTools/lib/Bio/DB/Sam.pm>`_ 

In order to not interfere with other versions of these tools that you may have installed, these files are not copied into their common system-wide paths. All of the files required for breseq to function are created and collected in the ``stage`` directory of the installation.

3. Add :program:`breseq` to your $PATH
----------------------------------------

Breseq can now be run by invoking the executables located under ROOT/stage/bin. For convenience, you probably want to add this directory to your $PATH, so that you can invoke the commands without typing out the full path.

For a bash shell, you can run the command::

  echo "export PATH=\$PATH:BRESEQ_ROOT/stage/bin" >> ~/.profile
  
to do this, replacing [BRESEQ_ROOT] with the absolute path to the root of the :program:`breseq` source archive, e.g. "/Users/jbarrick/src/breseq".  
  
If you have other versions of SAMtools installed on your system, be careful about the order of paths in your $PATH variable. If you include the :program:`breseq` path *last*, then it will not override your commands going to the version you are normally using.

Common installation problems
---------------------------------

None known yet. If you have a problem, please contact breseq@barricklab.org.