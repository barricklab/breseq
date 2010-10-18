Installation
==============

|breseq| is a command line tool implemented in Perl, C++, and R. It will compile and function on a variety of Unix platforms, including MacOSX.

1. Install external dependencies
---------------------------------

Several external packages and software programs need to be installed to compile and use |breseq|:

* `GCC <http://gcc.gnu.org>`_ (or other C++ compiler) already installed on some systems
* `Perl <http://www.perl.org>`_ (version 5.8 or higher) already installed on most systems
* `BioPerl <http://www.bioperl.org>`_ (version 1.4 or higher)
* `SSAHA2 <http://www.sanger.ac.uk/resources/software/ssaha2/>`_ read mapping program
* `R <http://www.r-project.org>`_ (version 2.10 or higher) statistical programming language 
* `Boost <http://www.boost.org>`_ (version 1.42 or higher) C++ libraries

To install each missing dependency, use your system's package manager or visit the respective web pages linked above and follow the instructions for your platform. More specific directions are available below for some platforms. You must make sure that the executables for |SSAHA2| and :program:`R` are in your environment's $PATH for |breseq| to function.

MacOS X Instructions
********************

You will need administrator privileges to install |breseq| using these instructions. We recommend that you install and use the package manager `MacPorts <http://www.macports.org/>`_ to simplfy some installation steps.

* :program:`GCC`: download and install `Apple Developer tools <http://developer.apple.com/tools/>`_.
* :program:`Perl`: is already installed on MacOSX systems. 
* :program:`BioPerl`: download and install according to the directions at `BioPerl.org <http://www.bioperl.org>`_ 
* :program:`SSAHA2`: download and install the MacOSX package from the `Sanger Center <http://www.sanger.ac.uk/resources/software/ssaha2/>`_.
* :program:`R`: install with :program:`MacPorts` terminal command: ``sudo port install R``
* :program:`Boost`: install with :program:`MacPorts` terminal command: ``sudo port install boost``

2. Build :program:`breseq`
----------------------------

Now, open a terminal window, change directory to the root of the |breseq| source distribution, and run these commands::

  ./configure
  make
  make install
  
Optionally, you can test your |breseq| installation with this command::

  make test

These commands compile and install not only |breseq|, but also some open-source code developed by others. These packages are included in the |breseq| source distribution under ``extern``:

* `SAMtools <http://samtools.sourceforge.net>`_ 
* `Bio::DB::Sam <http://search.cpan.org/~lds/Bio-SamTools/lib/Bio/DB/Sam.pm>`_ 

In order to not interfere with other versions of these tools that you may have installed, these files are not copied into common system-wide paths. All of the files required for |breseq| to function will be in the newly created ``breseq`` directory within the source path after installation. After you have successfully built |breseq|, you can move the entirety of this directory to any location on your system and |breseq| will continue to function.

3. Add |breseq| to your $PATH
----------------------------------------

|breseq| can now be run by invoking the executables located under ``breseq/bin``. For convenience, you may want to add this directory to your environment's $PATH variable so that you can invoke |breseq| commands from a shell without typing the full path.

For a bash shell, you can do this by running the command::

  echo "export PATH=\$PATH:[LOCATION]/breseq/bin" >> ~/.profile
  
replacing ``[LOCATION]`` with the absolute path to the root of the |breseq| source archive, for example, to make ``/Users/me/my_programs/breseq/bin``.  You will need to open a new shell after you do this for the change to take effect.
  
If you already have |SAMtools| installed on your system, be careful about the order of paths in your $PATH environmental variable. Include the |breseq| path *last*, so that you will not override your previously installed version.

Common installation problems
---------------------------------

None known yet. If you have a problem, please contact jeffrey.e.barrick@gmail.com.


