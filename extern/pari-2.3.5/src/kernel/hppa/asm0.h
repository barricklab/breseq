#line 2 "../src/kernel/hppa/asm0.h"
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
   __asm__ ("add %2,%3,%0\n\taddc %%r0,%%r0,%1" \
        : "=r" (__value), "=r" (overflow) \
        : "r" (__arg1), "r" (__arg2) \
        : "cc"); \
  __value; \
})

#define addllx(a,b) \
({ ulong __value, __arg1 = (a), __arg2 = (b); \
   __asm__ ("sub %4,%5,%%r0\n\taddc %2,%3,%0\n\taddc %%r0,%%r0,%1" \
        : "=r" (__value), "=r" (overflow) \
        : "r" (__arg1), "r" (__arg2), "r" (overflow), "r" ((ulong) 1)\
        : "cc"); \
  __value; \
})

#define subll(a,b) \
({ ulong __value, __arg1 = (a), __arg2 = (b); \
   __asm__ ("sub %2,%3,%0\n\taddc %%r0,%%r0,%1\n\tsubi 1,%1,%1" \
        : "=r" (__value), "=r" (overflow) \
        : "r" (__arg1), "r" (__arg2) \
        : "cc"); \
  __value; \
})

#define subllx(a,b) \
({ ulong __value, __arg1 = (a), __arg2 = (b); \
   __asm__ ("sub %%r0,%4,%%r0\n\tsubb %2,%3,%0\n\taddc %%r0,%%r0,%1\n\tsubi 1,%1,%1" \
        : "=r" (__value), "=r" (overflow) \
        : "r" (__arg1), "r" (__arg2), "r" (overflow)\
        : "cc"); \
  __value; \
})

#define mulll(a,b) \
({ ulong __arg1 = (a), __arg2 = (b); \
   union {double z; ulong x[2];} __vtab; \
   __asm__ ("xmpyu %1,%2,%0" \
        : "=f" (__vtab.z) \
        : "f" (__arg1), "f" (__arg2) \
        : "cc"); \
   hiremainder=__vtab.x[0]; \
   __vtab.x[1]; \
})

#define addmul(a,b) \
({ ulong __value, __arg1 = (a), __arg2 = (b); \
    union {double z; ulong x[2];} __vtab; \
    __asm__ ("xmpyu %1,%2,%0" \
	: "=f" (__vtab.z) \
	: "f" (__arg1), "f" (__arg2) \
	: "cc"); \
    __asm__ ("add %2,%3,%0\n\taddc %%r0, %4, %1" \
        : "=&r" (__value), "=r" (hiremainder) \
        : "r" (__vtab.x[1]),"r" (hiremainder), "r" (__vtab.x[0]) \
        : "cc"); \
    __value; \
})

#endif
