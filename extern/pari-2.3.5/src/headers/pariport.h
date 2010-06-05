/* $Id: pariport.h 6593 2005-02-13 15:12:55Z kb $

Copyright (C) 2000  The PARI group.

This file is part of the PARI/GP package.

PARI/GP is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation. It is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY WHATSOEVER.

Check the License for details. You should have received a copy of it, along
with the package; see the file 'COPYING'. If not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA. */

/*******************************************************************/
/*                                                                 */
/*           DECLARATIONS SPECIFIC TO THE PORTABLE VERSION         */
/*                                                                 */
/*******************************************************************/

/*******************************************************************/
/*                                                                 */
/*                    OPERATIONS BY VALUE                          */
/*            f is a pointer to the function called.               */
/*            result is gaffect-ed to last parameter               */
/*                                                                 */
/*******************************************************************/
#define TRgopgz(f, x, y, fstr)  STMT_START {\
  GEN __x = (x), __y = (y);\
  long prec = precision(__y);\
  pari_sp __av = avma;\
  if (!prec) err(infprecer, fstr);\
  gaffect(f(__x), __y); avma=__av; } STMT_END
#define gopgz(f, x, y)  STMT_START {\
  GEN __x = (x), __y = (y);\
  pari_sp __av = avma;\
  gaffect(f(__x), __y); avma=__av; } STMT_END
#define gopggz(f, x, y, z)  STMT_START {\
  GEN __x = (x), __y = (y), __z = (z);\
  pari_sp __av = avma;\
  gaffect(f(__x,__y), __z); avma=__av; } STMT_END
#define gopgsz(f, x, s, z)  STMT_START {\
  GEN __x = (x), __z = (z);\
  long __s = (s);\
  pari_sp __av = avma;\
  gaffect(f(__x,__s), __z); avma=__av; } STMT_END
#define gopsgz(f, s, y, z)  STMT_START {\
  GEN __y = (y), __z = (z);\
  long __s = (s);\
  pari_sp __av = avma;\
  gaffect(f(__s,__y), __z); avma=__av; } STMT_END
#define gopssz(f, x, y, z)  STMT_START {\
  GEN __z = (z);\
  long __x = (x), __y = (y);\
  pari_sp __av = avma;\
  gaffect(f(__x,__y), __z); avma=__av; } STMT_END

/* mp.c */

#define equaliu(x,y)    (equalui((y),(x)))
#define equalis(x,y)    (equalsi((y),(x)))
#define cmpiu(x,y)     (-cmpui((y),(x)))
#define cmpis(x,y)     (-cmpsi((y),(x)))
#define cmprs(x,y)     (-cmpsr((y),(x)))
#define cmpri(x,y)     (-cmpir((y),(x)))
#define subis(x,y)     (addsi(-(y),(x)))
#define subrs(x,y)     (addsr(-(y),(x)))

#define truedivii(a,b) (truedvmdii((a),(b),NULL))
#define truedivis(a,b) (truedvmdis((a),(b),NULL))
#define divii(a,b)     (dvmdii((a),(b),NULL))
#define remii(a,b)     (dvmdii((a),(b),ONLY_REM))

#define mpshift(x,s)   ((typ(x)==t_INT)?shifti((x),(s)):shiftr((x),(s)))
#define mptruncz(x,y)  gopgz(mptrunc,(x),(y))
#define mpfloorz(x,y)  gopgz(mpfloor,(x),(y))
#define mpaddz(x,y,z)  gopggz(mpadd,(x),(y),(z))
#define addsiz(s,y,z)  gopsgz(addsi,(s),(y),(z))
#define addsrz(s,y,z)  gopsgz(addsr,(s),(y),(z))
#define addiiz(x,y,z)  gopggz(addii,(x),(y),(z))
#define addirz(x,y,z)  gopggz(addir,(x),(y),(z))
#define addriz(x,y,z)  gopggz(addir,(y),(x),(z))
#define addrrz(x,y,z)  gopggz(addrr,(x),(y),(z))
#define mpsubz(x,y,z)  gopggz(mpsub,(x),(y),(z))
#define subss(x,y)     (addss(-(y),(x)))
#define subssz(x,y,z)  (addssz((x),-(y),(z)))
#define subsiz(s,y,z)  gopsgz(subsi,(s),(y),(z))
#define subsrz(s,y,z)  gopsgz(subsr,(s),(y),(z))
#define subisz(y,s,z)  gopsgz(addsi,-(s),(y),(z))
#define subrsz(y,s,z)  gopsgz(addsr,-(s),(y),(z))
#define subiiz(x,y,z)  gopggz(subii,(x),(y),(z))
#define subirz(x,y,z)  gopggz(subir,(x),(y),(z))
#define subriz(x,y,z)  gopggz(subri,(x),(y),(z))
#define subrrz(x,y,z)  gopggz(subrr,(x),(y),(z))
#define mpmulz(x,y,z)  gopggz(mpmul,(x),(y),(z))
#define mulsiz(s,y,z)  gopsgz(mulsi,(s),(y),(z))
#define mulsrz(s,y,z)  gopsgz(mulsr,(s),(y),)(z)
#define muliiz(x,y,z)  gopggz(mulii,(x),(y),(z))
#define mulirz(x,y,z)  gopggz(mulir,(x),(y),(z))
#define mulriz(x,y,z)  gopggz(mulir,(y),(x),(z))
#define mulrrz(x,y,z)  gopggz(mulrr,(x),(y),(z))
#define mpdvmdz(x,y,z,t) (dvmdiiz((x),(y),(z),(t))
#define addssz(s,y,z)  gopssz(addss,(s),(y),(z))
#define modssz(s,y,z)  gopssz(modss,(s),(y),(z))
#define mulssz(s,y,z)  gopssz(mulss,(s),(y),(z))
#define modsiz(s,y,z)  gopsgz(modsi,(s),(y),(z))
#define modisz(y,s,z)  gopgsz(modis,(y),(s),(z))
#define remsiz(s,y,z)  gopsgz(remsi,(s),(y),(z))
#define remisz(y,s,z)  gopgsz(remis,(y),(s),(z))
#define remssz(s,y,z)  gopssz(remss,(s),(y),(z))
#define diviiz(x,y,z)  gopggz(divii,(x),(y),(z))
#define divirz(x,y,z)  gopggz(divir,(x),(y),(z))
#define divisz(x,y,z)  gopgsz(divis,(x),(y),(z))
#define divriz(x,y,z)  gopggz(divri,(x),(y),(z))
#define divsrz(s,y,z)  gopsgz(divsr,(s),(y),(z))
#define divrsz(y,s,z)  gopgsz(divrs,(y),(s),(z))
