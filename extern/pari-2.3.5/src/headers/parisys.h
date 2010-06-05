/* $Id: parisys.h 7873 2006-04-14 15:56:31Z kb $

Copyright (C) 2000  The PARI group.

This file is part of the PARI/GP package.

PARI/GP is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation. It is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY WHATSOEVER.

Check the License for details. You should have received a copy of it, along
with the package; see the file 'COPYING'. If not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA. */

/* This files contains macros depending on system and compiler    */

#ifndef LITTLE_ENDIAN_64
#  define   LITTLE_ENDIAN_64 12345678
#endif
#ifndef BIG_ENDIAN_64
#  define   BIG_ENDIAN_64    87654321
#endif
#ifndef LITTLE_ENDIAN
#  define   LITTLE_ENDIAN 1234
#endif
#ifndef BIG_ENDIAN
#  define   BIG_ENDIAN    4321
#endif
#ifndef PDP_ENDIAN
#  define   PDP_ENDIAN    3412
#endif

#ifdef __cplusplus
#  define ANYARG ...
#  define BEGINEXTERN extern "C" {
#  define ENDEXTERN }
#else
#  define ANYARG
#  define BEGINEXTERN
#  define ENDEXTERN
#endif

#ifdef DISABLE_INLINE
#  undef ASMINLINE
#else
#  ifdef __cplusplus
#    define INLINE inline static
#  elif defined(__GNUC__)
#    define INLINE __inline__ static
#  endif
#endif

#ifndef DISABLE_VOLATILE
#  ifdef __GNUC__
#    define VOLATILE volatile
#  endif
#endif

#ifndef VOLATILE
#  define VOLATILE
#endif
#ifndef INLINE
#  define INLINE static
#endif

#if defined(_WIN32) || defined(__CYGWIN32__)
/* ANSI C does not allow to longjmp() out of a signal handler, in particular,
 * the SIGINT handler. On Win32, the handler is executed in another thread, and
 * longjmp'ing into another thread's stack will utterly confuse the system.
 * Instead, we check whether win32ctrlc is set in new_chunk(). */
BEGINEXTERN
  extern int win32ctrlc;
  void dowin32ctrlc();
ENDEXTERN
#endif
