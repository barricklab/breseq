Installation
==============

:program:`breseq` is a command line tool implemented in Perl/C++/R. It should compile and function on a variety of Unix platforms, including MacOSX.

1. Install external dependencies
---------------------------------

:program:`breseq` requires several external packages and software programs to be installed.

* `GCC <http://gcc.gnu.org>`_ (or other C++ compiler) already installed on many systems
* `Perl <http://www.perl.org>`_ (version 5.8 or higher) already installed on many systems
* `SSAHA2 <http://www.sanger.ac.uk/resources/software/ssaha2/>`_ read mapping program
* `R <http://www.r-project.org>`_ (version 2.10 or higher) statistical programming language 

To install each missing dependency, use your system's package manager or visit the respective web pages linked above and follow the instructions for your platform.

MacOSX Instructions
********************

You must have administrator privileges to install :program:`breseq` using these instructions. We recommend that you install and use the package manager `MacPorts <http://www.macports.org/>`_ to simplfy some installation steps.

* :program:`GCC`: download and install `Apple Developer tools <http://developer.apple.com/tools/>`_. 
* :program:`Perl`: is already installed on MacOSX systems. 
* :program:`SSAHA2`: Download and install the appropriate package from the `Sanger Center <http://www.sanger.ac.uk/resources/software/ssaha2/>`_.
* :program:`R`: install with :program:`MacPorts` command: ``%% sudo port install R``

2. Build :program:`breseq`
----------------------------

Before compiling :program:`breseq`, move the entire archive that you downloaded to a stable location on your system. If you later move this directory, then you may need to re-compile :program:`breseq` to have it operate correctly.

Open a terminal window, navigate to the root of the source distribution and run these commands::

  %% ./configure
  %% make
  %% make install
  %% make test (`this step is optional`)

These commands compile and install not only :program:`breseq`, they also install several tools developed by others that are included in the source package in the /extern directory for ease of installation.

They include:

* `BioPerl <http://www.bioperl.org>`_
* `SAMtools <http://samtools.sourceforge.net>`_ 
* `Bio::DB::Sam <http://search.cpan.org/~lds/Bio-SamTools/lib/Bio/DB/Sam.pm>`_ 
* `Boost <http://www.boost.org>`_

In order to not interfere with other versions of these tools that you may have installed, these files are not copied into their common system-wide paths. All of the files required for breseq to function are created and collected in the ``stage`` directory of the installation.

3. Add :program:`breseq` to your $PATH
----------------------------------------

Breseq can now be run by invoking the executables located under ROOT/stage/bin. For convenience, you probably want to add this directory to your $PATH, so that you can invoke the commands without typing the full path.

For a bash shell, you can run the command::

  echo "export PATH=\$PATH:BRESEQ_ROOT/stage/bin" >> ~/.profile
  
to do this, replacing [BRESEQ_ROOT] with the absolute path to the root of the :program:`breseq` source archive, e.g. "/Users/jbarrick/src/breseq".  
  
If you have other versions of SAMtools installed on your system, be careful about the order of paths in your $PATH variable. If you include the :program:`breseq` path *last*, then it will not override your commands going to the version you are normally using.

Common installation problems
---------------------------------

None known yet. If you have a problem, please contact breseq@barricklab.org.