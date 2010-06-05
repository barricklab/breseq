#line 2 "../src/kernel/none/divll.h"
/* $Id: divll.h 7867 2006-04-14 15:26:51Z kb $

Copyright (C) 2003  The PARI group.

This file is part of the PARI/GP package.

PARI/GP is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation. It is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY WHATSOEVER.

Check the License for details. You should have received a copy of it, along
with the package; see the file 'COPYING'. If not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA. */

/* This file originally adapted from gmp-3.1.1 (from T. Granlund), files
 * longlong.h and gmp-impl.h

  Copyright (C) 2000 Free Software Foundation, Inc. */

#undef  LOCAL_HIREMAINDER
#define LOCAL_HIREMAINDER
extern ulong hiremainder;

#if !defined(INLINE)
extern long divll(ulong x, ulong y);
#else

#define __GLUE(hi, lo) (((hi) << BITS_IN_HALFULONG) | (lo))
#define __SPLIT(a, b, c) b = HIGHWORD(a); c = LOWWORD(a)
#define __LDIV(a, b, q, r) q = a / b; r = a - q*b
extern ulong hiremainder;

/* divide (hiremainder * 2^BITS_IN_LONG + n0) by d; assume hiremainder < d.
 * Return quotient, set hiremainder to remainder */

#if defined(__GNUC__) && !defined(DISABLE_INLINE)
#undef LOCAL_HIREMAINDER
#define LOCAL_HIREMAINDER register ulong hiremainder

#define divll(n0, d) 							\
({ 									\
  ulong __d1, __d0, __q1, __q0, __r1, __r0, __m, __n1, __n0;	        \
  ulong __k, __d;                                                       \
                                                                        \
  __n1 = hiremainder; __n0 = n0; __d = d;                               \
  if (__n1 == 0)                                                        \
  { /* Only one division needed */                                      \
    __LDIV(__n0, __d, __q1, hiremainder);                               \
  }                                                                     \
  else if (__d < LOWMASK)                                               \
  { /* Two half-word divisions  */                                      \
    __n1 = __GLUE(__n1, HIGHWORD(__n0));                                \
    __LDIV(__n1, __d, __q1, __r1);                                      \
    __n1 = __GLUE(__r1,  LOWWORD(__n0));                                \
    __LDIV(__n1, __d, __q0, hiremainder);                               \
    __q1 = __GLUE(__q1, __q0);                                          \
  }                                                                     \
  else                                                                  \
  { /* General case */                                                  \
    if (__d & HIGHBIT)                                                  \
    {                                                                   \
      __k = 0; __SPLIT(__d, __d1, __d0);                                \
    }                                                                   \
    else                                                                \
    {                                                                   \
      __k = bfffo(__d);                                                 \
      __n1 = (__n1 << __k) | (__n0 >> (BITS_IN_LONG - __k));            \
      __n0 <<= __k;                                                     \
      __d = __d << __k; __SPLIT(__d, __d1, __d0);                       \
    }                                                                   \
    __LDIV(__n1, __d1, __q1, __r1);                                     \
    __m =  __q1 * __d0; 					        \
    __r1 = __GLUE(__r1, HIGHWORD(__n0));  				\
    if (__r1 < __m)							\
    {									\
      __q1--, __r1 += __d;						\
      if (__r1 >= __d) /* we didn't get carry when adding to __r1 */    \
        if (__r1 < __m)	__q1--, __r1 += __d;				\
    }									\
    __r1 -= __m;							\
    __LDIV(__r1, __d1, __q0, __r0);                                     \
    __m =  __q0 * __d0;  					        \
    __r0 = __GLUE(__r0, LOWWORD(__n0));   				\
    if (__r0 < __m)							\
    {									\
      __q0--, __r0 += __d;						\
      if (__r0 >= __d)			        			\
        if (__r0 < __m)	__q0--, __r0 += __d;				\
    }									\
    hiremainder = (__r0 - __m) >> __k;		                        \
    __q1 = __GLUE(__q1, __q0);                 				\
  }                                   				        \
  __q1;					                                \
})

#else /* __GNUC__ */

INLINE long
divll(ulong n0, ulong d) 							
{
  ulong __d1, __d0, __q1, __q0, __r1, __r0, __m, __n1, __n0;
  ulong __k, __d;

  __n1 = hiremainder; __n0 = n0; __d = d;

  if (__n1 == 0)
  { /* Only one division needed */
    __LDIV(__n0, __d, __q1, hiremainder);
  }
  else if (__d < LOWMASK)
  { /* Two half-word divisions  */
    __n1 = __GLUE(__n1, HIGHWORD(__n0));
    __LDIV(__n1, __d, __q1, __r1);
    __n1 = __GLUE(__r1,  LOWWORD(__n0));
    __LDIV(__n1, __d, __q0, hiremainder);
    __q1 = __GLUE(__q1, __q0);
  }
  else
  { /* General case */
    if (__d & HIGHBIT)
    {
      __k = 0; __SPLIT(__d, __d1, __d0);
    }
    else
    {
      __k = bfffo(__d);
      __n1 = (__n1 << __k) | (__n0 >> (BITS_IN_LONG - __k));
      __n0 = __n0 << __k;
      __d = __d << __k; __SPLIT(__d, __d1, __d0);
    }
    __LDIV(__n1, __d1, __q1, __r1);
    __m =  __q1 * __d0;
    __r1 = __GLUE(__r1, HIGHWORD(__n0));
    if (__r1 < __m)
      {
        __q1--, __r1 += __d;
        if (__r1 >= __d) /* we didn't get carry when adding to __r1 */
          if (__r1 < __m) __q1--, __r1 += __d;
      }
    __r1 -= __m;
    __LDIV(__r1, __d1, __q0, __r0);
    __m =  __q0 * __d0;
    __r0 = __GLUE(__r0, LOWWORD(__n0));
    if (__r0 < __m)
      {
        __q0--, __r0 += __d;
        if (__r0 >= __d)
          if (__r0 < __m) __q0--, __r0 += __d;
      }
    hiremainder = (__r0 - __m) >> __k;
    __q1 = __GLUE(__q1, __q0);
  }
  return __q1;
}

#endif /* __GNUC__ */

#endif
