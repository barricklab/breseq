#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

#initialiazation
AC_PREREQ([2.65])

AC_INIT([FULL-PACKAGE-NAME], [VERSION], [BUG-REPORT-ADDRESS])

AC_CONFIG_AUX_DIR(aux_build)
AC_CONFIG_MACRO_DIR([aux_build/m4])
AC_CONFIG_HEADERS([aux_build/config.h])
AC_CONFIG_SRCDIR([])

AC_PREFIX_DEFAULT($PWD/stage)

LT_INIT
AM_INIT_AUTOMAKE

## Compiling with per-target flags
AM_PROG_CC_C_O

# FIXME: Replace `main' with a function in `-lz':
AC_CHECK_LIB([boost_serialization], [main])

# FIXME: Replace `main' with a function in `-lz':
AC_CHECK_LIB([boost_program_options], [main])

##########################################
### BEGIN content from  autoscan

# Checks for programs.
AC_PROG_CXX
AC_PROG_AWK
AC_PROG_CPP
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_PROG_MKDIR_P

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lbam':
AC_CHECK_LIB([bam], [main])
# FIXME: Replace `main' with a function in `-lcurses':
AC_CHECK_LIB([curses], [main])
# FIXME: Replace `main' with a function in `-lm':
AC_CHECK_LIB([m], [main])
# FIXME: Replace `main' with a function in `-lz':
AC_CHECK_LIB([z], [main])

# Checks for header files.
AC_CHECK_HEADERS([arpa/inet.h fcntl.h inttypes.h netdb.h stddef.h stdint.h stdlib.h string.h sys/socket.h unistd.h wchar.h])

# Checks for typedefs, structures, and compiler characteristics.
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

# Checks for library functions.
AC_FUNC_FSEEKO
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([bzero gethostbyaddr gethostbyname memmove memset pow select socket sqrt strchr strdup strerror strstr strtol])

### END content from  autoscan
##########################################


AC_CONFIG_FILES([extern/samtools-0.1.7a/Makefile
				 src/c/Makefile
                 Makefile ])

                 
AC_OUTPUT
