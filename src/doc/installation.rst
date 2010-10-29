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

MacOSX Instructions
********************

You will need administrator privileges to install |breseq| dependencies using these instructions. We recommend that you install and use the package manager `MacPorts <http://www.macports.org/>`_ to simplfy some installation steps.

* :program:`GCC`: download and install `Apple Developer tools <http://developer.apple.com/tools/>`_.
* :program:`Perl`: is already installed on MacOSX systems. 
* :program:`BioPerl`: download and install according to the directions at `BioPerl.org <http://www.bioperl.org>`_ 
* :program:`SSAHA2`: download and install the MacOSX package from the `Sanger Center <http://www.sanger.ac.uk/resources/software/ssaha2/>`_.
* :program:`R`: install with :program:`MacPorts` terminal command: ``sudo port install R``
* :program:`Boost`: install with :program:`MacPorts` terminal command: ``sudo port install boost``

2. Build :program:`breseq`
----------------------------

Now, open a terminal window, change directory to the root of the |breseq| source distribution. If you have admin privileges on your computer, you can run these commands to install |breseq|::

  ./configure
  make
  sudo make install

If you do not have admin privileges on your computer, then see :ref:`installing-in-a-user-location`.

These commands compile and install not only |breseq|, but also some open-source code developed by others. These packages are included in the |breseq| source distribution under ``extern``:

* `SAMtools <http://samtools.sourceforge.net>`_ 
* `Bio::DB::Sam <http://search.cpan.org/~lds/Bio-SamTools/lib/Bio/DB/Sam.pm>`_ 

.. WARNING::
   These commands will overwrite any other versions of :program:`SAMtools` or the Perl module :program:`Bio::DB::Sam` that you have in the default ./configure install locations. To avoid this, you can follow the instructions in :ref:`installing-in-a-user-location` to safely install |breseq| elsewhere.

Finally, we recommend that you test that your |breseq| installation functions with this command::

  make test
  
This should take 5-10 minutes to run and report success at the end if everything is operating correctly.

.. _installing-in-a-user-location:

Installing in a user location
*****************************

If you do not have admin privileges on your computer, then you need to specify a location in your home directory to install in. We'll assume that you've chosen to install |breseq| in ``/mnt/home/me/local``, in which case you would use these commands::

  ./configure --prefix=/mnt/home/me/local
  make
  make install

Before you test or use this kind of installation, you will need to tell your shell that ``/mnt/home/me/local`` contains a usual UNIX grouping of program directories (with sub-directories like ``bin``, ``lib``, ``man``, etc). To do this you can use these commands, if you are using a bash shell::

  echo "export PATH=\$PATH:/mnt/home/me/local/bin" >> ~/.profile

Now, you should be able to invoke |breseq| commands if you open a new terminal window.

Common installation problems
---------------------------------

None known yet. If you have a problem, please contact jeffrey.e.barrick@gmail.com.


