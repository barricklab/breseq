We use the Subversion open-source revision control system. For us developers,
it provides network-transparent source control. For ordinary users it provides
a convenient way to obtain patched versions in between releases. Subversion
clients are available for all major platforms: Unix, MacOS, Windows, see

  http://subversion.tigris.org/project_packages.html

In particular, the svn command-line client is readily available in all Linux
distributions.

This file documents access to the PARI Subversion server, which is intended for
PARI lovers who want the very latest bleeding edge release. These sources may
contain severe bugs, they may not even compile, benches may fail and so on.
Stable releases are made available on a regular basis using the customary
method: a message to pari-announce.

Note that in order to use the sources fetched via Subversion, you will need a
working perl installation as well as the regular build system.


1) First connection to the anonymous Subversion server:
=======================================================

Type the following command from the shell
  
  svn checkout svn://pari.math.u-bordeaux.fr/pari/trunk

It creates a local copy of the distribution from the distant repository in
local directory trunk, which you may move or rename if you wish, say to 'pari'.
From now on, you can go to this pari directory and use any svn command directly
(without the svn: argument), as long as you remain there, or in a subdirectory.

2) What can I do now ?
======================

* You can build pari in the usual way (see INSTALL) as if this 'pari' directory
had been created by fetching, then extracting, an archive on an FTP server.

* You can update your local copy at any time using 'svn update', which puts you
in synch with the repository.

* You can check exact differences between successive versions of a given file
by using 'svn diff'. If you modify some files on your local copy, this also
enables you to track down your changes, and produce a patch. You will not be
able to commit your changes using anonymous access. Send the output of svn
diff, to the pari-dev mailing list with a short description of what you have
done, or to pari@math.u-bordeaux.fr if you are not subscribed to pari-dev.

3) Version tags and branches:
=============================

Official releases (starting from version 2.0.17) are 'tagged' so that all files
pertaining to a given release can be simultaneously accessed without tracking
version numbers. Tag names are release-version with dots replaced by dashes,
e.g. release-2-0-20 for 2.0.20.

To fetch a specific version of pari (2.0.17 or more recent), type for instance

  svn checkout svn://pari.math.u-bordeaux.fr/pari/tags/release-2-0-20

The branch release-2-3-patches denotes the stable branch 2.3 as a whole, and
can be used to checkout up to date sources from that branch in between
releases. For instance:

  svn checkout svn://pari.math.u-bordeaux.fr/pari/branches/release-2-3-patches

produces the latest stable distribution with all relevant patches [the ones not
affecting stability] backported. 

Tips and Caveats:
=================

* svn diff gives you the difference between your local copy and the sources
they were based on, not with the current state of the repository. Use 

  svn diff -r HEAD

for that.

* To see the log message associated to the last commit leading to the current
state of your local repository, type 

  svn log -r COMMITTED

You may add a file or directory name to get the log message for the last commit
which modified it.
