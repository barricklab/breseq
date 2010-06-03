/* $Id: parigen.h 9115 2007-10-29 10:16:34Z kb $

Copyright (C) 2000  The PARI group.

This file is part of the PARI/GP package.

PARI/GP is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation. It is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY WHATSOEVER.

Check the License for details. You should have received a copy of it, along
with the package; see the file 'COPYING'. If not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA. */

/* This file defines the parameters of the GEN type               */

typedef long *GEN;
typedef unsigned long pari_ulong;
#define ulong pari_ulong

#ifdef LONG_IS_64BIT
#  define TWOPOTBYTES_IN_LONG  3
#else
#  define TWOPOTBYTES_IN_LONG  2
#endif

#define DEFAULTPREC    2 + (8>>TWOPOTBYTES_IN_LONG)
#define MEDDEFAULTPREC 2 + (16>>TWOPOTBYTES_IN_LONG)
#define BIGDEFAULTPREC 2 + (24>>TWOPOTBYTES_IN_LONG)
#define TWOPOTBITS_IN_LONG (TWOPOTBYTES_IN_LONG+3)
#define BYTES_IN_LONG (1L<<TWOPOTBYTES_IN_LONG)
#define BITS_IN_LONG  (1L<<TWOPOTBITS_IN_LONG)
#define HIGHBIT (1UL << (BITS_IN_LONG-1))
#define BITS_IN_HALFULONG (BITS_IN_LONG>>1)
#define MAXULONG (~0x0UL)
#define MAXHALFULONG ((1UL<<BITS_IN_HALFULONG) - 1)
#define LOWMASK  (MAXHALFULONG)
#define HIGHMASK (~LOWMASK)
#define SMALL_MASK (HIGHBIT>>1)

#define HIGHWORD(a) ((a) >> BITS_IN_HALFULONG)
#define LOWWORD(a) ((a) & LOWMASK)

/* Order of bits in codewords:
 *  x[0]       TYPBITS, CLONEBIT, LGBITS
 *  x[1].real  SIGNBITS, EXPOBITS
 *       int   SIGNBITS, LGBITS
 *       pol   SIGNBITS, VARNBITS
 *       ser   SIGNBITS, VARNBITS, VALPBITS
 *       padic VALPBITS, PRECPBITS */
#define TYPnumBITS   7
#define SIGNnumBITS  2

#ifdef LONG_IS_64BIT
#  define VARNnumBITS 16 /* otherwise MAXVARN too large */
#else
#  define VARNnumBITS 14
#endif

/* no user serviceable parts below :-) */
#define   LGnumBITS (BITS_IN_LONG - 1 - TYPnumBITS)
#define VALPnumBITS (BITS_IN_LONG - SIGNnumBITS - VARNnumBITS)
#define EXPOnumBITS (BITS_IN_LONG - SIGNnumBITS)
#define PRECPSHIFT VALPnumBITS
#define  VARNSHIFT VALPnumBITS
#define   TYPSHIFT (BITS_IN_LONG - TYPnumBITS)
#define  SIGNSHIFT (BITS_IN_LONG - SIGNnumBITS)

#define EXPOBITS    ((1UL<<EXPOnumBITS)-1)
#define SIGNBITS    (~((1UL<<SIGNSHIFT) - 1))
#define  TYPBITS    (~((1UL<< TYPSHIFT) - 1))
#define PRECPBITS   (~VALPBITS)
#define LGBITS      ((1UL<<LGnumBITS)-1)
#define VALPBITS    ((1UL<<VALPnumBITS)-1)
#define VARNBITS    (MAXVARN<<VARNSHIFT)
#define MAXVARN     ((1UL<<VARNnumBITS)-1)
#define HIGHEXPOBIT (1UL<<(EXPOnumBITS-1))
#define HIGHVALPBIT (1UL<<(VALPnumBITS-1))
#define CLONEBIT    (1UL<<LGnumBITS)

#define evaltyp(x)    (((ulong)(x)) << TYPSHIFT)
#define evalvarn(x)   (((ulong)(x)) << VARNSHIFT)
#define evalsigne(x)  ((ulong)(((long)(x)) << SIGNSHIFT))
#define evalprecp(x)   (((long)(x)) << PRECPSHIFT)
#define _evalexpo(x)  (HIGHEXPOBIT + (x))
#define _evalvalp(x)  (HIGHVALPBIT + (x))
#define evallgefint(x)  (x)
#define evallgeflist(x) (x)
#define _evallg(x)    (x)

#define typ(x)        ((long)((((ulong*)(x))[0]) >> TYPSHIFT))
#define settyp(x,s)   (((ulong*)(x))[0]=\
                        (((ulong*)(x))[0]&(~TYPBITS)) | evaltyp(s))

#define isclone(x)    (((ulong*) (x))[0] & CLONEBIT)
#define setisclone(x) (((ulong*) (x))[0] |= CLONEBIT)
#define unsetisclone(x) (((ulong*) (x))[0] &= (~CLONEBIT))

#define lg(x)         ((long)(((ulong*)(x))[0] & LGBITS))
#define setlg(x,s)    (((ulong*)(x))[0]=\
                        (((ulong*)(x))[0]&(~LGBITS)) | evallg(s))

#define signe(x)      ((((long*)(x))[1]) >> SIGNSHIFT)
#define setsigne(x,s) (((ulong*)(x))[1]=\
                        (((ulong*)(x))[1]&(~SIGNBITS)) | (ulong)evalsigne(s))
#define togglesign(x) (void)((((GEN)(x))[1] & SIGNBITS) && (((GEN)(x))[1] ^= HIGHBIT))

#define lgeflist(x)      (((long*)(x))[1])
#define setlgeflist(x,l) (((ulong*)(x))[1]=(ulong)(l))

#define lgefint(x)      ((long)(((ulong*)(x))[1] & LGBITS))
#define setlgefint(x,s) (((ulong*)(x))[1]=\
                          (((ulong*)(x))[1]&(~LGBITS)) | (ulong)evallgefint(s))

#define expo(x)       ((long) ((((ulong*)(x))[1] & EXPOBITS) - HIGHEXPOBIT))
#define setexpo(x,s)  (((ulong*)(x))[1]=\
		       (((ulong*)(x))[1]&(~EXPOBITS)) | (ulong)evalexpo(s))

#define valp(x)       ((long) ((((ulong*)(x))[1] & VALPBITS) - HIGHVALPBIT))
#define setvalp(x,s)  (((ulong*)(x))[1]=\
		       (((ulong*)(x))[1]&(~VALPBITS)) | (ulong)evalvalp(s))

#define precp(x)      ((long) (((ulong*)(x))[1] >> PRECPSHIFT))
#define setprecp(x,s) (((ulong*)(x))[1]=\
		       (((ulong*)(x))[1]&(~PRECPBITS)) | (ulong)evalprecp(s))

#define varn(x)       ((long)((((ulong*)(x))[1]&VARNBITS) >> VARNSHIFT))
#define setvarn(x,s)  (((ulong*)(x))[1]=\
		       (((ulong*)(x))[1]&(~VARNBITS)) | (ulong)evalvarn(s))

#define varncmp(x,y)  ((x)-(y))

