#line 2 "../src/kernel/hppa64/asm0.h"
/* $Id: asm0.h 11631 2009-02-28 22:03:38Z bill $

Copyright (C) 2004  The PARI group.

This file is part of the PARI/GP package.

PARI/GP is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation. It is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY WHATSOEVER.

Check the License for details. You should have received a copy of it, along
with the package; see the file 'COPYING'. If not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA. */

/* This file was made using idea from Bruno Haible ix86 asm inline kernel
 * and code from Nigel Smart hppa asm kernel. mulll was inspired from
 * longlong.h from the GNU MP package.*/

/*
ASM addll mulll
NOASM bfffo divll
*/
#ifdef ASMINLINE
#define LOCAL_HIREMAINDER  register ulong hiremainder
#define LOCAL_OVERFLOW     register ulong overflow

#define addll(a,b) \
({ ulong __value, __arg1 = (a), __arg2 = (b); \
   __asm__ ("add %2,%3,%0\n\tadd,dc %%r0,%%r0,%1" \
        : "=r" (__value), "=r" (overflow) \
        : "r" (__arg1), "r" (__arg2) \
        : "cc"); \
  __value; \
})

#define addllx(a,b) \
({ ulong __value, __arg1 = (a), __arg2 = (b); \
   __asm__ ("sub %4,%5,%%r0\n\tadd,dc %2,%3,%0\n\tadd,dc %%r0,%%r0,%1" \
        : "=r" (__value), "=r" (overflow) \
        : "r" (__arg1), "r" (__arg2), "r" (overflow), "r" ((ulong) 1)\
        : "cc"); \
  __value; \
})

#define subll(a,b) \
({ ulong __value, __arg1 = (a), __arg2 = (b); \
   __asm__ ("sub %2,%3,%0\n\tadd,dc %%r0,%%r0,%1\n\tsubi 1,%1,%1" \
        : "=r" (__value), "=r" (overflow) \
        : "r" (__arg1), "r" (__arg2) \
        : "cc"); \
  __value; \
})

#define subllx(a,b) \
({ ulong __value, __arg1 = (a), __arg2 = (b); \
   __asm__ ("sub %%r0,%4,%%r0\n\tsub,db %2,%3,%0\n\tadd,dc %%r0,%%r0,%1\n\tsubi 1,%1,%1" \
        : "=&r" (__value), "=r" (overflow) \
        : "r" (__arg1), "r" (__arg2), "r" (overflow)\
        : "cc"); \
  __value; \
})

/* z=a+b; c+= carry; return z */
#define __addllc(a,b,c) \
({ ulong __value, __arg1 = (a), __arg2 = (b); \
   __asm__ ("add %2,%3,%0\n\tadd,dc %4,%%r0,%1" \
        : "=&r" (__value), "=r" (c) \
        : "r" (__arg1), "r" (__arg2), "r" (c) \
        : "cc"); \
  __value; \
})

/* 32x32->64 multiply*/
#define __mulhh(a,b) \
({ unsigned int __arg1 = (a), __arg2 = (b); \
   ulong __value; \
   __asm__ ("xmpyu %1,%2,%0" \
        : "=f" (__value) \
        : "f" (__arg1), "f" (__arg2) \
        : "cc"); \
   __value; \
})

#define mulll(arg1,arg2) \
({ \
  const ulong __x=(arg1), __y=(arg2); \
  const ulong __xlo = LOWWORD(__x), __xhi = HIGHWORD(__x); \
  const ulong __ylo = LOWWORD(__y), __yhi = HIGHWORD(__y); \
  ulong __xylo,__xymid,__xyhi,__xymidhi,__xymidlo; \
  ulong __xylh,__xyhl; \
  __xylo = __mulhh(__xlo,__ylo); __xyhi = __mulhh(__xhi,__yhi); \
  __xylh = __mulhh(__xlo,__yhi); __xyhl = __mulhh(__xhi,__ylo); \
  __xymid = __xylh+__xyhl; \
  if (__xymid<__xylh) __xyhi += (1UL << BITS_IN_HALFULONG); \
  __xymidhi = HIGHWORD(__xymid); \
  __xymidlo = __xymid << BITS_IN_HALFULONG; \
  __xylo = __addllc(__xylo,__xymidlo,__xyhi); \
  hiremainder = __xyhi + __xymidhi; \
  __xylo; \
})

#define addmul(arg1,arg2) \
({ \
  const ulong __x=(arg1), __y=(arg2); \
  const ulong __xlo = LOWWORD(__x), __xhi = HIGHWORD(__x); \
  const ulong __ylo = LOWWORD(__y), __yhi = HIGHWORD(__y); \
  ulong __xylo,__xymid,__xyhi,__xymidhi,__xymidlo; \
  ulong __xylh,__xyhl; \
  __xylo = __mulhh(__xlo,__ylo); __xyhi = __mulhh(__xhi,__yhi); \
  __xylh = __mulhh(__xlo,__yhi); __xyhl = __mulhh(__xhi,__ylo); \
  __xymid = __xylh+__xyhl; \
  if (__xymid<__xylh) __xyhi += (1UL << BITS_IN_HALFULONG); \
  __xylo = __addllc(__xylo,hiremainder,__xyhi); \
  __xymidhi = HIGHWORD(__xymid); \
  __xymidlo = __xymid << BITS_IN_HALFULONG; \
  __xylo = __addllc(__xylo,__xymidlo,__xyhi); \
  hiremainder = __xyhi + __xymidhi; \
  __xylo; \
})

#endif
