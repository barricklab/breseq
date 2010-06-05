/* $Id: asm0.h 7881 2006-04-18 17:23:00Z kb $

Copyright (C) 2000  The PARI group.

This file is part of the PARI/GP package.

PARI/GP is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation. It is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY WHATSOEVER.

Check the License for details. You should have received a copy of it, along
with the package; see the file 'COPYING'. If not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA. */
/*
ASM addll mulll bfffo
NOASM divll
*/
#ifdef ASMINLINE
#define LOCAL_HIREMAINDER  register ulong hiremainder
#define LOCAL_OVERFLOW     register ulong overflow

#define addll(a, b)\
({ ulong __value, __arg1 = (a), __arg2 = (b); \
   __asm__ ("addc %0,%2,%3\n\txor %1,%2,%2\n\taddze %1,%4\n\t" \
   : "=&r" (__value), "=r" (overflow) \
   : "r" (__arg1), "r" (__arg2), "1" ((ulong) 0)); \
  __value; \
})

#define addllx(a, b)\
({ ulong __value, __arg1 = (a), __arg2 = (b); \
 __asm__ ("addc %0,%3,%4\n\tli %1,0\n\taddze %1,%4\n\taddc %0,%2,%5\n\taddze %1,%4\n\t" \
   : "=&r" (__value), "=r" (overflow) \
   : "r" (__arg1), "r" (__arg2), "1" (overflow), "0" ((ulong) 0)); \
__value; \
})

#define bfffo(a) \
({ ulong __a = (a), __value; \
    __asm__ ("cntlzw %0, %1" : "=r" (__value) : "r" (__a)); \
    __value; \
})

#define subll(a, b)\
({ ulong __value, __arg1 = (a), __arg2 = (b); \
  __asm__ ("subfc %0,%3,%2\n\tli %1,0\n\taddme %1,%4\n\tneg %1,%4" \
   : "=&r" (__value), "=r" (overflow) \
   : "r" (__arg1), "r" (__arg2), "1" ((ulong)0)); \
  __value; \
})

#define subllx(a, b)\
({ ulong __value, __arg1 = (a), __arg2 = (b); \
__asm__ ("subfc %0,%5,%2\n\tli %1,0\n\taddme %1,%5\n\tsubfc %0,%3,%4\n\taddme %1,%5\n\tneg %1,%5" \
   : "=r" (__value), "=r" (overflow) \
   : "r" (__arg1), "r" (__arg2), "0" ((ulong)0), "1" (overflow)); \
 __value; \
})

#define mulll(a, b) \
({ ulong __value, __arg1 = (a), __arg2 = (b); \
 __asm__ ("mulhwu %1,%2,%3\n\tmullw %0,%2,%3\n\t" \
   : "=r" (__value), "=&r" (hiremainder) \
   : "r" (__arg1), "r" (__arg2)); \
 __value; \
})

#define addmul(a, b) \
({ ulong __value, __arg1 = (a), __arg2 = (b), __temp; \
 __asm__ ("mullw %0,%3,%4\n\tmulhwu %2,%3,%4\n\taddc %0,%5,%6\n\taddze %1,%7\n\t" \
   : "=&r" (__value), "=r" (hiremainder), "=r" (__temp) \
   : "r" (__arg1), "r" (__arg2), "0" ((ulong) 0), "1" (hiremainder), "2" ((ulong) 0)); \
 __value; \
})
#endif
