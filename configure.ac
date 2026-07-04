###
##
## AUTHORS
##
##   Jeffrey E. Barrick <jeffrey.e.barrick@gmail.com> and other contributors
##
## LICENSE AND COPYRIGHT
##
##   Copyright (c) 2008-2010 Michigan State University
##   Copyright (c) 2011-2025 The University of Texas at Austin
##   Copyright (c) 2025-     Michigan State University
##
##   breseq is free software; you can redistribute it and/or modify it under the
##   terms of the GNU General Public License as published by the Free Software
##   Foundation; either version 2, or (at your option) any later version.
##
##   SPDX-License-Identifier: GPL-2.0-or-later
##
###


#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
AC_PREREQ([2.65])
AC_INIT([breseq], [0.50.0], [https://github.com/barricklab/breseq/issues], [breseq], [http://barricklab.org/breseq])
AC_CONFIG_AUX_DIR(aux_build)
AC_CONFIG_MACRO_DIR([aux_build/m4])
AC_CONFIG_HEADERS([aux_build/config.h])
AC_CONFIG_SRCDIR([])

AC_PREFIX_DEFAULT([$prefix])

LT_INIT
AM_INIT_AUTOMAKE([subdir-objects])

AM_PROG_CC_C_O
AC_PROG_CXX
AC_PROG_AWK
AC_PROG_CPP
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_PROG_MKDIR_P

AX_COMPUTE_RELATIVE_PATHS([bindir:datadir:bin_to_data])
AC_MSG_NOTICE([Local path from bin to data...])
AC_MSG_NOTICE([$bin_to_data])
AC_SUBST([BIN_TO_DATA_RELATIVE_PATH], [$bin_to_data])


AC_CHECK_HEADERS([arpa/inet.h fcntl.h inttypes.h netdb.h stddef.h stdint.h stdlib.h string.h sys/socket.h unistd.h wchar.h])
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_INT16_T
AC_TYPE_INT32_T
AC_TYPE_INT64_T
AC_TYPE_INT8_T
AC_TYPE_OFF_T
AC_TYPE_SIZE_T
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T
AC_TYPE_UINT8_T
AC_FUNC_FSEEKO
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([bzero gethostbyaddr gethostbyname memmove memset pow select socket sqrt strchr strdup strerror strstr strtol ftello])

#AC_SEARCH_LIBS([wclear], [curses ncurses], [],
#   [AC_MSG_ERROR(["lib(n)curses is required."])])

AC_SEARCH_LIBS([acosf], [m], [], [AC_MSG_ERROR(["libm is required."])])
AC_SEARCH_LIBS([pthread_create], [pthread], [], [AC_MSG_ERROR(["libpthread is required."])])


#### libz

# Dynamic linking to libz (bioconda install) — verify at ${prefix}/lib and set
# LIBZ_LINK/STATIC_LIBZ_INC. (binarydist.sh links statically against its own
# universal build of zlib by passing -I/-L paths and LIBS directly to
# ./configure, ahead of these -lfoo flags resolving against ${prefix}/lib.)

save_LDFLAGS="$LDFLAGS"
LDFLAGS="$LDFLAGS -L${prefix}/lib"
AC_CHECK_LIB([z], [deflateSetHeader], [],
  [AC_MSG_ERROR([libz >=1.2.2.1 not found at ${prefix}/lib. Install zlib or set --prefix=\$CONDA_PREFIX])])
LDFLAGS="$save_LDFLAGS"
AC_SUBST([LIBZ_LINK], ["-lz"])
AC_SUBST([STATIC_LIBZ_INC], [""])


#### miniz (required for --zip-html evidence ZIP archive support)

# Dynamic linking to libminiz (bioconda install) — verify at ${prefix} and set
# MINIZ_LINK/MINIZ_INC. (binarydist.sh links statically against its own
# universal build of miniz by passing -I/-L paths and LIBS directly to
# ./configure, ahead of these -lfoo flags resolving against ${prefix}/lib.)

save_CPPFLAGS="$CPPFLAGS"
save_LDFLAGS="$LDFLAGS"
CPPFLAGS="$CPPFLAGS -I${prefix}/include"
LDFLAGS="$LDFLAGS -L${prefix}/lib"
AC_CHECK_HEADER([miniz/miniz.h], [],
  [AC_MSG_ERROR([miniz/miniz.h not found at ${prefix}/include. Install: conda install -c conda-forge miniz and configure with --prefix=\$CONDA_PREFIX])])
AC_CHECK_LIB([miniz], [mz_zip_writer_init_file], [],
  [AC_MSG_ERROR([libminiz not found at ${prefix}/lib. Install: conda install -c conda-forge miniz and configure with --prefix=\$CONDA_PREFIX])])
CPPFLAGS="$save_CPPFLAGS"
LDFLAGS="$save_LDFLAGS"
AC_SUBST([MINIZ_INC],  [""])
AC_SUBST([MINIZ_LINK], ["-lminiz"])


# On Linux, export symbols so dladdr() can resolve names from the main binary
case "$host_os" in
  darwin*) ;;
  *) LDFLAGS="$LDFLAGS -rdynamic" ;;
esac
# On older Linux glibc (< 2.34), dladdr() is in libdl; newer glibc and macOS have it in libc
AC_SEARCH_LIBS([dladdr], [dl], [], [AC_MSG_ERROR(["libdl is required."])])

#### htslib (required for BAM/SAM reading and pileup)

# Dynamic linking to libhts (bioconda install) — verify at ${prefix} and set
# HTSLIB_LINK/HTSLIB_INC. (binarydist.sh links statically against its own
# universal build of htslib by passing -I/-L paths and LIBS directly to
# ./configure, ahead of these -lfoo flags resolving against ${prefix}/lib.)

save_CPPFLAGS="$CPPFLAGS"
save_LDFLAGS="$LDFLAGS"
save_LIBS="$LIBS"
CPPFLAGS="$CPPFLAGS -I${prefix}/include"
LDFLAGS="$LDFLAGS -L${prefix}/lib"
# $CRAM_CODEC_LIBS (set by binarydist.sh to full .a paths for bzip2/xz) must be
# present for this link check when statically linking: libhts.a has undefined
# CRAM-codec symbols that only those archives resolve, and AC_CHECK_LIB's probe
# program needs to link cleanly to detect -lhts. Dynamic/conda builds leave it
# unset, so this is a no-op there (libhts.so resolves its own dependencies).
LIBS="$LIBS $CRAM_CODEC_LIBS"
AC_CHECK_HEADER([htslib/sam.h], [],
  [AC_MSG_ERROR([htslib/sam.h not found at ${prefix}/include. Install: conda install -c bioconda htslib and configure with --prefix=\$CONDA_PREFIX])])
AC_CHECK_LIB([hts], [hts_open], [],
  [AC_MSG_ERROR([libhts not found at ${prefix}/lib. Install: conda install -c bioconda htslib and configure with --prefix=\$CONDA_PREFIX])])
# Check whether this htslib version has sam_hdr_destroy (added in 1.10)
LIBS="$LIBS -lhts"
AC_CHECK_FUNC([sam_hdr_destroy],
  [AC_DEFINE([HAVE_SAM_HDR_DESTROY], [1], [htslib has sam_hdr_destroy (>= 1.10)])])
CPPFLAGS="$save_CPPFLAGS"
LDFLAGS="$save_LDFLAGS"
LIBS="$save_LIBS"
AC_SUBST([HTSLIB_INC],  [""])
AC_SUBST([HTSLIB_LINK], ["-lhts"])

# Lets binarydist.sh hand the final breseq/gdtools program links full .a paths
# to htslib's CRAM-codec dependencies (bzip2/xz) — e.g.
# ./configure CRAM_CODEC_LIBS="/path/to/libbz2.a /path/to/liblzma.a". Routed
# through breseq_LDADD/gdtools_LDADD (not the global LIBS) so libtool doesn't
# also see these static archives while linking the libbreseq.la convenience
# library and emit "Linking the shared library ... against the static
# library ... is not portable" warnings; conda/dynamic builds leave it unset.
AC_ARG_VAR([CRAM_CODEC_LIBS], [Static libraries (e.g. full paths to libbz2.a and liblzma.a) to append to the breseq/gdtools link, after libhts, for CRAM codec support])


AC_SUBST([LIBZ], [extern/zlib-1.3.1])
AC_SUBST([BRESEQ], [src/breseq])


# Reading .git/refs/heads/master directly only worked for an in-place checkout
# of the master branch with unpacked refs. Ask git itself instead, so this
# also works from a worktree (where .git is a file, not a directory), on any
# branch/detached HEAD, and after `git gc` has packed the refs.
AC_PATH_PROG([GIT], [git])
GITHUB_REVISION_HASH=
GITHUB_REVISION_PREFIX=
if test -n "$GIT" && $GIT -C "$srcdir" rev-parse --git-dir >/dev/null 2>&1; then
	GITHUB_REVISION_HASH=`$GIT -C "$srcdir" rev-parse --short=12 HEAD 2>/dev/null`
	test -n "$GITHUB_REVISION_HASH" && GITHUB_REVISION_PREFIX="revision "
fi
AC_SUBST(GITHUB_REVISION_HASH)
AC_SUBST(GITHUB_REVISION_PREFIX)

AC_DEFINE_UNQUOTED([GITHUB_REVISION_STRING], ["${GITHUB_REVISION_PREFIX}${GITHUB_REVISION_HASH}"], [Defining GitHub repository version])


AC_CANONICAL_HOST
# Check for which host we are on and setup a few things
# specifically based on the host
case $host_os in
  cygwin* )
  		AC_MSG_NOTICE([Detected host 'cygwin'.])
  		AC_SUBST(SAMTOOLS_ADDITIONAL_CFLAGS, ["-Dexpl=exp -Dlogl=log"])
        ;;
  darwin* )
    	AC_MSG_NOTICE([Detected host 'darwin'.])
        ;;
esac


AC_CONFIG_FILES([
    src/breseq/Makefile
    Makefile
])

AC_OUTPUT
