# ======================================================================================
#  http://www.gnu.org/software/autoconf-archive/ax_compute_standard_relative_paths.html
# ======================================================================================
#
# SYNOPSIS
#
#   AX_COMPUTE_STANDARD_RELATIVE_PATHS
#
# DESCRIPTION
#
#   Here is the standard hierarchy of paths, as defined by the GNU Coding
#   Standards:
#
#     prefix
#        exec_prefix
#           bindir
#           libdir
#           libexecdir
#           sbindir
#        datadir
#        sysconfdir
#        sharestatedir
#        localstatedir
#        infodir
#        lispdir
#        includedir
#        oldincludedir
#        mandir
#
#   This macro will setup a set of variables of the form
#   'xxx_forward_relative_path' and 'xxx_backward_relative_path' where xxx
#   is one of the above directories. The latter variable is set to the
#   relative path to go from xxx to its parent directory, while the former
#   hold the other way.
#
#   For instance `bindir_relative_path' will contains the value to add to
#   $exec_prefix to reach the $bindir directory (usually 'bin'), and
#   `bindir_backward_relative_path' the value to append to $bindir to reach
#   the $exec_prefix directory (usually '..').
#
#   This macro requires AX_COMPUTE_RELATIVE_PATHS which itself requires
#   AX_NORMALIZE_PATH.
#
# LICENSE
#
#   Copyright (c) 2008 Alexandre Duret-Lutz <adl@gnu.org>
#
#   This program is free software; you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation; either version 2 of the License, or (at your
#   option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#   Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program. If not, see <http://www.gnu.org/licenses/>.
#
#   As a special exception, the respective Autoconf Macro's copyright owner
#   gives unlimited permission to copy, distribute and modify the configure
#   scripts that are the output of Autoconf when processing the Macro. You
#   need not follow the terms of the GNU General Public License when using
#   or distributing such scripts, even though portions of the text of the
#   Macro appear in them. The GNU General Public License (GPL) does govern
#   all other use of the material that constitutes the Autoconf Macro.
#
#   This special exception to the GPL applies to versions of the Autoconf
#   Macro released by the Autoconf Archive. When you make and distribute a
#   modified version of the Autoconf Macro, you may extend this special
#   exception to the GPL to apply to your modified version as well.

#serial 7

AU_ALIAS([ADL_COMPUTE_STANDARD_RELATIVE_PATHS], [AX_COMPUTE_STANDARD_RELATIVE_PATHS])
AC_DEFUN([AX_COMPUTE_STANDARD_RELATIVE_PATHS],
## These calls need to be on separate lines for aclocal to work!
[AX_COMPUTE_RELATIVE_PATHS(dnl
AX_STANDARD_RELATIVE_PATH_LIST)])

dnl AX_STANDARD_RELATIVE_PATH_LIST
dnl ===============================
dnl A list of standard paths, ready to supply to AX_COMPUTE_RELATIVE_PATHS.
AC_DEFUN([AX_STANDARD_RELATIVE_PATH_LIST],
[pushdef([TRIPLET],
[$][1:$][2:$][2_forward_relative_path $]dnl
[2:$][1:$][2_backward_relative_path])dnl
TRIPLET(prefix, exec_prefix) dnl
TRIPLET(exec_prefix, bindir) dnl
TRIPLET(exec_prefix, libdir) dnl
TRIPLET(exec_prefix, libexecdir) dnl
TRIPLET(exec_prefix, sbindir) dnl
TRIPLET(prefix, datadir) dnl
TRIPLET(prefix, sysconfdir) dnl
TRIPLET(prefix, sharestatedir) dnl
TRIPLET(prefix, localstatedir) dnl
TRIPLET(prefix, infodir) dnl
TRIPLET(prefix, lispdir) dnl
TRIPLET(prefix, includedir) dnl
TRIPLET(prefix, oldincludedir) dnl
TRIPLET(prefix, mandir) dnl
popdef([TRIPLET])])
