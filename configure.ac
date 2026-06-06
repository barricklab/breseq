###
##
## AUTHORS
##
## Jeffrey E. Barrick <jeffrey.e.barrick@gmail.com>
## David B. Knoester
##
## LICENSE AND COPYRIGHT
##
## Copyright (c) 2010 Michigan State University
##
## breseq is free software; you can redistribute it and/or modify it under the terms the
## GNU General Public License as published by the Free Software Foundation; either
## version 1, or (at your option) any later version.
##
###


#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
AC_PREREQ([2.65])
AC_INIT([breseq], [0.40.1], [https://github.com/barricklab/breseq/issues], [breseq], [http://barricklab.org/breseq])
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

# We have options for static linking (precompiled binaries) and dynamic linking (bioconda install)
# Check for command line option and environmental variable versions of --with-static-lib and WITH_STATIC_LIBZ
AC_ARG_WITH([static-libz],
  [AS_HELP_STRING([--with-static-libz=PATH], [Prefix to static install of libz. Directory that has lib and include within it.])],
  [with_static_libz="$withval"],
  [with_static_libz="${WITH_STATIC_LIBZ:-}"]
)
AC_MSG_NOTICE([with_static_libz = $with_static_libz])

# 1) static linking to libz if command line option / environmental variable provided and library is found
# 2) dynamic linking to libz otherwise — verify at ${prefix}/lib and set LIBZ_LINK

AS_IF([test "x$with_static_libz" != "x"], [
  AS_IF([test -f "$with_static_libz/lib/libz.a"], [
    AC_MSG_NOTICE([Static linking to libz.a located in "$with_static_libz/lib"])
    AC_SUBST([STATIC_LIBZ_INC], ["-I$with_static_libz/lib/include"])
    AC_SUBST([LIBZ_LINK], ["$with_static_libz/lib/libz.a"])
  ], [
    AC_MSG_ERROR([Could not find libz.a library for static linking in specified path: "$with_static_libz/lib"])
  ])
], [
  save_LDFLAGS="$LDFLAGS"
  LDFLAGS="$LDFLAGS -L${prefix}/lib"
  AC_CHECK_LIB([z], [deflateSetHeader], [],
    [AC_MSG_ERROR([libz >=1.2.2.1 not found at ${prefix}/lib. Install zlib or set --prefix=\$CONDA_PREFIX])])
  LDFLAGS="$save_LDFLAGS"
  AC_SUBST([LIBZ_LINK], ["-lz"])
])


#### miniz (required for --zip-html evidence ZIP archive support)

# We have options for static linking (precompiled binaries) and dynamic linking (bioconda install)
# Check for command line option and environmental variable versions of --with-static-miniz and WITH_STATIC_MINIZ
AC_ARG_WITH([static-miniz],
  [AS_HELP_STRING([--with-static-miniz=PATH], [Prefix to static install of miniz. Directory that has lib/ and include/ within it.])],
  [with_static_miniz="$withval"],
  [with_static_miniz="${WITH_STATIC_MINIZ:-}"]
)
AC_MSG_NOTICE([with_static_miniz = $with_static_miniz])

# 1) static linking to libminiz if command line option / environmental variable provided and library is found
# 2) dynamic linking to libminiz — verify at ${prefix} and set MINIZ_LINK/MINIZ_INC

AS_IF([test "x$with_static_miniz" != "x"], [
  AS_IF([test -f "$with_static_miniz/lib/libminiz.a"], [
    AC_MSG_NOTICE([Static linking to libminiz.a located in "$with_static_miniz/lib"])
    AC_SUBST([MINIZ_INC],  ["-I$with_static_miniz/include"])
    AC_SUBST([MINIZ_LINK], ["$with_static_miniz/lib/libminiz.a"])
  ], [
    AC_MSG_ERROR([Could not find libminiz.a for static linking in: "$with_static_miniz/lib"])
  ])
], [
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
])


#Find libunwind
AC_ARG_WITH([libunwind],
	AS_HELP_STRING([--without-libunwind], [Disable preferring to use libunwind for debug info]))

AS_IF([test "x$with_libunwind" != "xno"], [
  AC_CHECK_HEADER(libunwind.h, FOUND_LIBUNWIND_H=1, FOUND_LIBUNWIND_H=0)
  AC_SEARCH_LIBS([unw_getcontext], [unwind], FOUND_LIBUNWIND=1,  FOUND_LIBUNWIND=0)

  if [[ $FOUND_LIBUNWIND == 0 ]] && [[ $FOUND_LIBUNWIND_H == 1 ]]; then
    AC_MSG_CHECKING([libunwind])
    AC_MSG_CHECKING([for target architecture])
    case x"$target" in
    xNONE | x)
    target_or_host="$host" ;;
    *)
    target_or_host="$target" ;;
    esac
    AC_MSG_RESULT([$target_or_host])

    case "$target_or_host" in
    i*86-*)
    TARGET_PLATFORM="x86"
    ;;

    x86_64-*)
    TARGET_PLATFORM="x86_64"
    ;;

    arm-*)
    TARGET_PLATFORM="arm"
    ;;
    esac

    AC_SEARCH_LIBS([_UL${TARGET_PLATFORM}_step], [unwind], FOUND_LIBUNWIND=1, [])

    if [[ $FOUND_LIBUNWIND == 1 ]]; then
        UNWIND_LIBS="-lunwind-$TARGET_PLATFORM"
        AC_SUBST(UNWIND_LIBS)
    fi
  fi
],
[
  FOUND_LIBUNWIND_H=0
  FOUND_LIBUNWIND=0
])


if [[ $FOUND_LIBUNWIND == 1 ]] && [[ $FOUND_LIBUNWIND_H == 1 ]]; then
    AC_MSG_RESULT([Using libunwind.h for backtrace.])
    AC_DEFINE([HAVE_LIBUNWIND],[1],[Lib unwind])
else
    AC_MSG_RESULT([Using execinfo.h for backtrace.])
fi



#### htslib (required for BAM/SAM reading and pileup)

# We have options for static linking (precompiled binaries) and dynamic linking (bioconda install)
AC_ARG_WITH([static-htslib],
  [AS_HELP_STRING([--with-static-htslib=PATH], [Prefix to static install of htslib. Directory that has lib/ and include/ within it.])],
  [with_static_htslib="$withval"],
  [with_static_htslib="${WITH_STATIC_HTSLIB:-}"]
)
AC_MSG_NOTICE([with_static_htslib = $with_static_htslib])

# 1) static linking to libhts if command line option / environmental variable provided and library is found
# 2) dynamic linking to libhts via prefix otherwise

AS_IF([test "x$with_static_htslib" != "x"], [
  AS_IF([test -f "$with_static_htslib/lib/libhts.a"], [
    AC_MSG_NOTICE([Static linking to libhts.a located in "$with_static_htslib/lib"])
    AC_SUBST([HTSLIB_INC],  ["-I$with_static_htslib/include"])
    AC_SUBST([HTSLIB_LINK], ["$with_static_htslib/lib/libhts.a"])
    # Check whether this htslib version has sam_hdr_destroy (added in 1.10)
    save_CPPFLAGS="$CPPFLAGS"; save_LIBS="$LIBS"
    CPPFLAGS="$CPPFLAGS -I$with_static_htslib/include"
    LIBS="$LIBS $with_static_htslib/lib/libhts.a"
    CPPFLAGS="$save_CPPFLAGS"; LIBS="$save_LIBS"
  ], [
    AC_MSG_ERROR([Could not find libhts.a for static linking in: "$with_static_htslib/lib"])
  ])
], [
  save_CPPFLAGS="$CPPFLAGS"
  save_LDFLAGS="$LDFLAGS"
  save_LIBS="$LIBS"
  CPPFLAGS="$CPPFLAGS -I${prefix}/include"
  LDFLAGS="$LDFLAGS -L${prefix}/lib"
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
])


AC_SUBST([LIBZ], [extern/zlib-1.3.1])
AC_SUBST([BRESEQ], [src/c/breseq])
AC_SUBST([GDTOOLS], [src/c/gdtools])


AC_CHECK_FILE([.git/refs/heads/master],
	AC_SUBST(GITHUB_REVISION_HASH, [`eval head -n 1 .git/refs/heads/master | sed 's/^\(.\{12\}\).*/\1/'`]),
	AC_SUBST(GITHUB_REVISION_HASH, [])
)

AC_CHECK_FILE([.git/refs/heads/master],
	AC_SUBST(GITHUB_REVISION_PREFIX, "revision "),
	AC_SUBST(GITHUB_REVISION_PREFIX, [])
)

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
    src/c/breseq/Makefile
    Makefile
])

AC_OUTPUT
