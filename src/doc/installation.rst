Installation
==============

|breseq| is a command line tool implemented in C++, and R. It will compile and function on a variety of UNIX platforms, including MacOSX. |breseq| installation from the source code requires some basic familiarity with UNIX commands and environments.

1. Install external dependencies
---------------------------------

Several external packages and software programs need to be installed to compile and use |breseq|:

* `GCC <http://gcc.gnu.org>`_ (or other C++ compiler) already installed on many systems
* `SSAHA2 <http://www.sanger.ac.uk/resources/software/ssaha2/>`_ read mapping program
* `R <http://www.r-project.org>`_ (version 2.1.0 or higher) statistical programming language 

To install each missing dependency, use your system's package manager or visit the respective web pages linked above and follow the instructions for your platform. More specific directions are available below for some platforms. You must make sure that the executables for |SSAHA2| and :program:`R` are in your environment's $PATH for |breseq| to function.

MacOSX Instructions
********************

You will need administrator privileges to install |breseq| dependencies using these instructions. We recommend that you install and use the package manager `MacPorts <http://www.macports.org/>`_ to simplify these installation steps. Using  :program:`MacPorts` will generally take longer than downloading and installing the packages in other ways, but it greatly simplifies the searching you might otherwise have to do to track down all the prerequisites.

* :program:`GCC`: Download and install the `Apple Developer tools <http://developer.apple.com/tools/>`_. You can use either version 3 (which is free) or download version 4 from the App Store for a small price.
* :program:`SSAHA2`: download and install the MacOSX package from the `Sanger Center <http://www.sanger.ac.uk/resources/software/ssaha2/>`_. You will need to move the executables to where your system can use them. If you change into the ssaha2_v2.5.1_MacOS directory, you can use this command:

>>> sudo cp ssaha2* /usr/local/bin

You should now get a message like this, telling you that the system can find your :program:`SSAHA2`: executable, if you type this command in the terminal:

>>> which ssaha2
/users/local/bin/ssaha2

* :program:`R`: Download the installer package from http://www.r-project.org/.

2. Compile and install |breseq|
-------------------------------

If you have admin privileges and want to install |breseq| in a standard location accessible to all users of a computer, then see :ref:`installing-in-a-system-wide-location`. If you do not have admin privileges on your computer, then see :ref:`installing-in-the-source-directory` or :ref:`installing-in-a-custom-location`. 

.. NOTE::
   If you encounter problems with one of the other install methods, we recommend that you try :ref:`installing-in-the-source-directory`.   
   
.. _installing-in-a-system-wide-location:

Installing in a system-wide location
************************************

This method requires that you have admin privileges on your machine. After installation, all users of the machine will be able to run |breseq|.

Open a terminal window and change directory to the root of the |breseq| source distribution. Then, run these commands::

  ./configure
  make
  sudo make install

These commands compile and install not only |breseq|, but also some open-source code developed by others. These packages are included in the |breseq| source distribution under /extern:

* `SAMtools <http://samtools.sourceforge.net>`_ 

.. WARNING::
   Installing |breseq| will overwrite any other versions of :program:`SAMtools` that you have in the default ./configure install locations. To avoid this, you can follow the instructions in :ref:`installing-in-the-source-directory` or :ref:`installing-in-a-custom-location` to safely install |breseq| elsewhere.

Finally, we recommend that you test that your |breseq| installation functions with this command after installatin::

  make test
  
This should take less than 5 minutes to run and report success at the end if everything is operating correctly.

.. _installing-in-the-source-directory:

Installing in the source directory
**********************************

This is the most robust way to install |breseq| if you do not have admin privileges on a system. All of the compiled programs and libraries will be self-contained in the original source tree.

Open a terminal window and change directory to the root of the |breseq| source distribution. Then, run these commands::

  ./configure --prefix=${PWD}
  make
  make install
  make test

After installation, if you want to be able to call |breseq| commands without specifying the entire path to them, you will need to add the newly created "bin" directory within the |breseq| source to your $PATH.

For a :program:`bash` shell you can usually use a command like this::

  echo "export PATH=\$PATH:${PWD}/bin" >> ~/.bashrc

But the exact way to do this may depend on your system. Once you open a new terminal window so that it registers this change to your $PATH, you should be able to invoke |breseq| commands.

.. _installing-in-a-custom-location:

Installing in a custom location
*******************************

We'll assume that you've chosen to install |breseq| in ``/mnt/home/me/local``. Open a terminal window and change directory to the root of the |breseq| source distribution. Then, run these commands::

  ./configure --prefix=/mnt/home/me/local
  make
  make install

This will create a usual UNIX grouping of program directories (with sub-directories like ``bin``, ``lib``, ``man``, etc). 

After installation, if you want to be able to call |breseq| commands without specifying the entire path to them, you will need to add the newly created "bin" directory within the |breseq| source to your $PATH.

For a :program:`bash` shell you can usually use a command like this::

  echo "export PATH=\$PATH:/mnt/home/me/local/bin" >> ~/.bashrc

But the exact way to do this may depend on your system. You may also want to similarly update your $MANPATH, $CPPFLAGS, $LD_FLAGS, etc. Now you should be able to invoke |breseq| commands once you open a new terminal window.

Common installation problems
---------------------------------

Dependencies installed in custom locations
******************************************

In general, you will need to be sure that your environment is set up correctly to find and use each dependency. This will likely be taken care of for you if you use an installer package. If you install some dependencies from source or in custom locations, and run into problems with |breseq| installation, be sure to check that:

#. :program:`R` is in your $PATH.
#. :program:`SSAHA2` is in your $PATH.

.. note::
   You may need to use absolute paths (i.e. ``/absolute/path``) rather than paths relative to your home directory (i.e ``~/path/relative/to/home``) for these settings.

Other problems
***************

If you have a problem installing |breseq|, please send a detailed report to jeffrey.e.barrick@gmail.com.


