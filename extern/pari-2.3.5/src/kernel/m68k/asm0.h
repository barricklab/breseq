#line 2 "../src/kernel/m68k/asm0.h"
/* $Id: asm0.h 7911 2006-05-10 20:51:17Z kb $

Copyright (C) 2006  The PARI group.

This file is part of the PARI/GP package.

PARI/GP is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation. It is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY WHATSOEVER.

Check the License for details. You should have received a copy of it, along
with the package; see the file 'COPYING'. If not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA. */

/* Written by Bill Allombert and dedicated to thoses who wrote the original 
 * m68k kernel mp.s */

/*
ASM addll mulll bfffo divll
*/

#ifdef ASMINLINE
#define LOCAL_HIREMAINDER  register ulong hiremainder
#define LOCAL_OVERFLOW     register ulong overflow

#define addll(a,b)                                              \
({ ulong __value, __arg1 = (a), __arg2 = (b);                   \
   __asm__ ("add.l %2,%0 ; addx.l %1,%1"                        \
        : "=&d" (__value), "=d" (overflow)                      \
        : "rm" (__arg1), "0" (__arg2), "1" (0UL)                \
        : "cc");                                                \
  __value;                                                      \
})

#define addllx(a,b)                                             \
({ ulong __value, __arg1 = (a), __arg2 = (b), __temp;           \
   __asm__ ("neg.l %2 ; addx.l %4,%0 ; addx.l %1,%1"            \
        : "=d" (__value), "=d" (overflow), "=d" (__temp)        \
        : "0" (__arg1), "d" (__arg2), "2" (overflow), "1" (0UL) \
        : "cc");                                                \
  __value;                                                      \
})

#define subll(a,b)                                              \
({ ulong __value, __arg1 = (a), __arg2 = (b);                   \
   __asm__ ("sub.l %3,%0 ; addx.l %1,%1"                        \
        : "=&d" (__value), "=d" (overflow)                      \
        : "0" (__arg1), "rm" (__arg2), "1" (0UL)                \
        : "cc");                                                \
  __value;                                                      \
})

#define subllx(a,b)                                             \
({ ulong __value, __arg1 = (a), __arg2 = (b), __temp;           \
   __asm__ ("neg.l %2 ; subx.l %4,%0 ; addx.l %1,%1"            \
        : "=d" (__value), "=d" (overflow), "=d" (__temp)        \
        : "0" (__arg1), "d" (__arg2), "2" (overflow), "1" (0UL) \
        : "cc");                                                \
  __value;                                                      \
})

#define mulll(a, b)                                                     \
({                                                                      \
  ulong __arg1 = (a), __arg2 = (b), __value;                            \
  __asm__ ("mulu.l %2, %0:%1"                                           \
           : "=d" (hiremainder), "=d" (__value)                         \
           : "md" (__arg1) , "1" (__arg2)                               \
           : "cc");                                                     \
  __value;                                                              \
})

#define addmul(a, b)                                                    \
({                                                                      \
  ulong __arg1 = (a), __arg2 = (b), __value;                            \
  __asm__ ("mulu.l %2, %0:%1; add.l %4,%1; addx.l %5,%0"                \
           : "=&d" (hiremainder), "=&d" (__value)                       \
           : "md" (__arg1), "1" (__arg2), "d" (hiremainder), "d" (0UL)  \
           : "cc" );                                                    \
  __value;                                                              \
})

#define bfffo(a)                                                        \
({                                                                      \
  ulong __arg1 = (a), __value;                                          \
  __asm__ ("bfffo %1{#0:#0}, %0"                                        \
           : "=d" (__value)                                             \
           : "md" (__arg1)                                              \
           : "cc" );                                                    \
  __value;                                                              \
})

#define divll(a, b)                                                     \
({                                                                      \
  ulong __arg2 = (b), __value =(a);                                     \
  __asm__ ("divu.l %2, %0:%1"                                           \
           : "+d" (hiremainder), "+d" (__value)                         \
           : "md" (__arg2)                                              \
           : "cc");                                                     \
  __value;                                                              \
})
#endif
