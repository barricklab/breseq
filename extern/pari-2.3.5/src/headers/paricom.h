/* $Id: paricom.h 11924 2009-09-16 09:25:40Z bill $

Copyright (C) 2004  The PARI group.

This file is part of the PARI/GP package.

PARI/GP is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation. It is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY WHATSOEVER.

Check the License for details. You should have received a copy of it, along
with the package; see the file 'COPYING'. If not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA. */

/******************************************************************/
/*                                                                */
/*              PARI header file (common to all versions)         */
/*                                                                */
/******************************************************************/
#ifdef STMT_START /* perl headers */
#  undef STMT_START
#endif
#ifdef STMT_END
#  undef STMT_END
#endif
/* STMT_START { statements; } STMT_END;
 * can be used as a single statement, as in
 * if (x) STMT_START { ... } STMT_END; else ...
 * [ avoid "dangling else" problem in macros ] */
#define STMT_START	do
#define STMT_END	while (0)
/*=====================================================================*/
/* CATCH(numer) {
 *   recovery 
 * } TRY {
 *   code
 * } ENDCATCH
 * will execute 'code', then 'recovery' if exception 'numer' is thrown
 * [ any exception if numer == CATCH_ALL ].
 * RETRY = as TRY, but execute 'recovery', then 'code' again [still catching] */
#define CATCH(err) {         \
  VOLATILE long __err = err; \
  int pari_errno;            \
  jmp_buf __env;             \
  void *__catcherr = NULL;   \
  if ((pari_errno = setjmp(__env))) 

#define RETRY { __catcherr = err_catch(__err, &__env); {
#define TRY else RETRY

/* Take address of __catcher to prevent compiler from putting it into a register
 * (could be clobbered by longjmp otherwise) */
#define CATCH_RELEASE() err_leave(&__catcherr)
#define ENDCATCH }} CATCH_RELEASE(); }

#define CATCH_ALL -1
/*=====================================================================*/
/* VOLATILE int errorN;
 * CATCH_ERR(errorN) {
 *   code
 * } ENDCATCH_ERR
 * executes 'code', setting errorN to the number of exception thrown;
 * errorN is 0 if no error was thrown. */

#define CATCH_ERR(__err) {  \
  jmp_buf __env;            \
  __err = setjmp(__env);    \
  if (!__err) {		    \
    void *__catcherr = err_catch(CATCH_ALL, &__env);

#define ENDCATCH_ERR	    \
    CATCH_RELEASE();	    \
  }}
/*=====================================================================*/

#define LOG2   (0.6931471805599453) /* log(2) */
#define L2SL10 (0.3010299956639812) /* log_10(2) */

#ifndef  PI
#  define PI (3.141592653589)
#endif

/*3.32~log_2(10)*/
#define ndec2nlong(x) (1 + (long)((x)*(3.321928094887362/BITS_IN_LONG)))
#define ndec2prec(x) (3 + (long)((x)*(3.321928094887362/BITS_IN_LONG)))
#define nbits2prec(x) (((x)+3*BITS_IN_LONG-1) >> TWOPOTBITS_IN_LONG)
#define nbits2nlong(x) (((x)+BITS_IN_LONG-1) >> TWOPOTBITS_IN_LONG)
#define nchar2nlong(x) (((x)+BYTES_IN_LONG-1) >> TWOPOTBYTES_IN_LONG)
#define bit_accuracy(x) (((x)-2) << TWOPOTBITS_IN_LONG)
#define bit_accuracy_mul(x,y) (((x)-2) * (BITS_IN_LONG*(y)))
#define prec2ndec(x) ((long)bit_accuracy_mul((x), L2SL10))
#define GSTR(x) ((char*) (((GEN) (x)) + 1 ))

#include "pariold.h"

/* Common global variables: */
extern ulong DEBUGFILES, DEBUGLEVEL, DEBUGMEM, precdl;
extern long  *ordvar;
extern GEN   bernzone,gpi,geuler;
extern GEN   polvar,*pol_1,*pol_x,primetab;
extern GEN   gen_m1,gen_1,gen_2,ghalf,gi,gen_0,gnil;

extern const long lontyp[];
extern void* global_err_data;

extern int new_galois_format;
extern int factor_add_primes;

enum manage_var_t {
  manage_var_create,
  manage_var_delete,
  manage_var_init,
  manage_var_next,
  manage_var_max_avail,
  manage_var_pop
};

#ifdef LONG_IS_64BIT
#  define VERYBIGINT (9223372036854775807L) /* 2^63-1 */
#  define BIGINT (2147483647)               /* 2^31-1 */
#  define u_OK_ULONG(p) ((ulong)p <= 3037000493UL)
#else
#  define VERYBIGINT (2147483647L) /* 2^31-1 */
#  define BIGINT (32767)           /* 2^15-1 */
#  define u_OK_ULONG(p) ((ulong)p <= 46337UL)
#endif

/* 2p^2 < 2^BITS_IN_LONG */
#define OK_ULONG(p) (lgefint(p) == 3 && u_OK_ULONG(p[2]))

#ifndef HAS_EXP2
#  undef exp2
#  define exp2(x) (exp((double)(x)*LOG2))
#endif
#ifndef HAS_LOG2
#  undef log2
#  define log2(x) (log((double)(x))/LOG2)
#endif

#ifdef min
#  undef min
#endif
#ifdef max
#  undef max
#endif
#define min(a,b) ((a)>(b)?(b):(a))
#define max(a,b) ((a)>(b)?(a):(b))

#define gval(x,v) (ggval((x),pol_x[v]))

#define absr  mpabs
#define absi  mpabs
#define negi  mpneg
#define negr  mpneg
#define rcopy mpcopy

#define addis(x,s)  (addsi((s),(x)))
#define addrs(x,s)  (addsr((s),(x)))
#define addri(x,s)  (addir((s),(x)))
#define mulis(x,s)  (mulsi((s),(x)))
#define muliu(x,s)  (mului((s),(x)))
#define mulru(x,s)  (mulur((s),(x)))
#define mulri(x,s)  (mulir((s),(x)))
#define mulrs(x,s)  (mulsr((s),(x)))

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
  if (!prec) pari_err(infprecer, fstr);\
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

#define mpexpz(x,y)    gopgz(mpexp,(x),(y))
#define mplogz(x,y)    gopgz(mplog,(x),(y))
#define mpcosz(x,y)    gopgz(mpcos,(x),(y))
#define mpsinz(x,y)    gopgz(mpsin,(x),(y))
#define gnegz(x,y)     gopgz(gneg,(x),(y))

#define gachz(x,y)     TRgopgz(gach,(x),(y), "gachz")
#define gacosz(x,y)    TRgopgz(gacos,(x),(y), "gacosz")
#define gashz(x,y)     TRgopgz(gash,(x),(y), "gashz")
#define gasinz(x,y)    TRgopgz(gasin,(x),(y), "gasinz")
#define gatanz(x,y)    TRgopgz(gatan,(x),(y), "gatanz")
#define gathz(x,y)     TRgopgz(gath,(x),(y), "gathz")
#define gchz(x,y)      TRgopgz(gch,(x),(y), "gchz")
#define gcosz(x,y)     TRgopgz(gcos,(x),(y), "gcosz")
#define gcotanz(x,y)   TRgopgz(gcotan,(x),(y), "gcotanz")
#define gexpz(x,y)     TRgopgz(gexp,(x),(y), "gexpz")
#define ggamdz(x,y)    TRgopgz(ggamd,(x),(y), "ggamdz")
#define ggammaz(x,y)   TRgopgz(ggamma,(x),(y), "ggammaz")
#define glngammaz(x,y) TRgopgz(glngamma,(x),(y), "glngammaz")
#define glogz(x,y)     TRgopgz(glog,(x),(y), "glogz")
#define gpsiz(x,y)     TRgopgz(gpsi,(x),(y), "gpsiz")
#define gshz(x,y)      TRgopgz(gsh,(x),(y), "gshz")
#define gsinz(x,y)     TRgopgz(gsin,(x),(y), "gsinz")
#define gsqrtz(x,y)    TRgopgz(gsqrt,(x),(y), "gsqrtz")
#define gtanz(x,y)     TRgopgz(gtan,(x),(y), "gtanz")
#define gthz(x,y)      TRgopgz(gth,(x),(y), "gthz")
#define gzetaz(x,y)    TRgopgz(gzeta,(x),(y), "gzetaz")

#define gabsz(x,prec,y) gopgsz(gabs,(x),(prec),(y))
#define gmaxz(x,y,z)    gopggz(gmax,(x),(y),(z))
#define gminz(x,y,z)    gopggz(gmin,(x),(y),(z))
#define gaddz(x,y,z)    gopggz(gadd,(x),(y),(z))
#define gsubz(x,y,z)    gopggz(gsub,(x),(y),(z))
#define gmulz(x,y,z)    gopggz(gmul,(x),(y),(z))
#define gdivz(x,y,z)    gopggz(gdiv,(x),(y),(z))
#define gdiventz(x,y,z) gopggz(gdivent,(x),(y),(z))
#define gmodz(x,y,z)    gopggz(gmod,(x),(y),(z))

#define gaddgs(y,s)     gaddsg((s),(y))
#define gcmpgs(s,y)     (-gcmpsg(y,s))
#define gdiventsg(s,y)  (gopsg2(gdivent,(s),(y)))
#define gdivsg(s,y)     (gopsg2(gdiv,(s),(y)))
#define gequalgs(s,y)    (gequalsg((y),(s)))
#define gmaxsg(s,y)     (gmaxgs((y),(s)))
#define gminsg(s,y)     (gmings((y),(s)))
#define gmodsg(s,y)     (gopsg2(gmod,(s),(y)))
#define gmulgs(y,s)     (gmulsg((s),(y)))
#define gsubgs(y,s)     gaddgs((y), -(s))
#define gsubsg(s,y)     (gopsg2(gsub,(s),(y)))

#define gaddgsz(y,s,z)    gopgsz(gaddgs,(y),(s),(z))
#define gaddsgz(s,y,z)    gaddgsz((y),(s),(z))
#define gdiventgsz(y,s,z) gopgsz(gdiventgs,(y),(s),(z))
#define gdiventsgz(s,y,z) gopsg2z(gdivent,(s),y),(z)
#define gdivgsz(y,s,z)    gopgsz(gdivgs,(y),(s),(z))
#define gdivsgz(s,y,z)    gopsg2z(gdiv,(s),(y),(z))
#define gmaxgsz(y,s,z)    gopgsz(gmaxgs,(y),(s),(z))
#define gmaxsgz(s,y,z)    gmaxgsz((y),(s),(z))
#define gmingsz(y,s,z)    gopgsz(gmings,(y),(s),(z))
#define gminsgz(s,y,z)    gmingsz((y),(s),(z))
#define gmodgsz(y,s,z)    gopgsz(gmodgs,(y),(s),(z))
#define gmodsgz(s,y,z)    gopsg2z(gmod,(s),(y),(z))
#define gmulgsz(y,s,z)    gopsgz(gmulsg,(s),(y),(z))
#define gmulsgz(s,y,z)    gmulgsz((y),(s),(z))
#define gsubgsz(y,s,z)    gopgsz(gaddgs,(y),-(s),(z))
#define gsubsgz(s,y,z)    gopsg2z(gsub,(s),(y),(z))

#define gmul2nz(x,s,z)  gopgsz(gmul2n,(x),(s),(z))
#define gshiftz(x,s,z)  gopgsz(gshift,(x),(s),(z))

#define bern(i)       (bernzone + 3 + (i)*bernzone[2])

/* works only for POSITIVE integers */
#define modBIL(x)  (*int_LSW(x))
#define mod64(x)  (modBIL(x) & 63)
#define mod32(x)  (modBIL(x) & 31)
#define mod16(x)  (modBIL(x) & 15)
#define mod8(x)   (modBIL(x) & 7)
#define mod4(x)   (modBIL(x) & 3)
#define mod2(x)   (modBIL(x) & 1)
#define is_bigint_lg(n,l) ((l)>3 || ((l)==3 && (((GEN)(n))[2] & HIGHBIT)))
#define is_pm1_lg(n,l)    ((l)==3 && ((GEN)(n))[2]==1)
#define is_pm1(n)    is_pm1_lg   (n, lgefint(n))
#define is_bigint(n) is_bigint_lg(n, lgefint(n))

#define degpol(a) ((long)lg(a)-3)
#define lgpol(a) ((long)lg(a)-2)

#define odd(x) ((x) & 1)
#define mpodd(x) (signe(x) && mod2(x))

#define ONLY_REM ((GEN*)0x1L)
#define ONLY_DIVIDES ((GEN*)0x2L)
#define RgX_div(x,y)     (RgX_divrem((x),(y),NULL))
#define RgX_rem(x,y)     (RgX_divrem((x),(y),ONLY_REM))
#define RgXQX_div(x,y,T) (RgXQX_divrem((x),(y),(T),NULL))
#define RgXQX_rem(x,y,T) (RgXQX_divrem((x),(y),(T),ONLY_REM))
#define gdeuc(x,y)       (poldivrem((x),(y),NULL))
#define grem(x,y)        (poldivrem((x),(y),ONLY_REM))
#define FpX_div(x,y,p)   (FpX_divrem((x),(y),(p), NULL))
#define FpX_rem(x,y,p)   (FpX_divrem((x),(y),(p), ONLY_REM))
#define Flx_div(x,y,p)   (Flx_divrem((x),(y),(p), NULL))

#define FpX_renormalize   ZX_renormalize
#define FpXX_renormalize  ZX_renormalize
#define FpXQX_renormalize ZX_renormalize

#define ZY_ZXY_resultant(a,b) ZY_ZXY_rnfequation((a),(b),NULL)

#define RgX_add gadd
#define RgX_sub gsub
#define RgX_neg gneg

#define ZX_mul RgX_mul
#define ZX_sqr RgX_sqr

#define zv_to_ZV(x)    (vecsmall_to_vec((x)))
#define zc_to_ZC(x)    (vecsmall_to_col((x)))
#define ZV_to_zv(x)    (vec_to_vecsmall((x)))
#define zx_to_zv(x,y)  (Flx_to_Flv((x),(y)))
#define zv_to_zx(x,y)  (Flv_to_Flx((x),(y)))
#define zm_to_zxV(x,y) (Flm_to_FlxV((x),(y)))
#define zero_zx(x)     (zero_Flx((x)))
#define polx_zx(x)     (polx_Flx((x)))
#define zx_shift(x,y)  (Flx_shift((x),(y)))

#define matpascal(n) matqpascal((n),NULL)
#define sturm(x) (sturmpart((x),NULL,NULL))
#define Z_issquare(x) (Z_issquarerem((x),NULL))
#define subres(x,y) (subresall((x),(y),NULL))
/* #define subres(x,y) (resultantducos((x),(y))) */

#define invmat(a) (gauss((a),NULL))

/* output of get_nf and get_bnf */
#define typ_NULL 0
#define typ_POL  1
#define typ_Q    2
#define typ_NF   3
#define typ_BNF  4
#define typ_BNR  5
#define typ_CLA  6 /* bnfclassunit   */
#define typ_ELL  7 /* elliptic curve */
#define typ_QUA  8 /* quadclassunit  */
#define typ_GAL  9 /* galoisinit     */
#define typ_BID 10
/* for gen_sort */
#define cmp_IND 1
#define cmp_LEX 2
#define cmp_REV 4
#define cmp_C   8

#define DIFFPTR_SKIP	255		/* Skip these entries */
#define NEXT_PRIME_VIADIFF(p,d)	 STMT_START \
  { while (*(d) == DIFFPTR_SKIP) (p) += *(d)++; (p) += *(d)++; } STMT_END
#define NEXT_PRIME_VIADIFF_CHECK(p,d)  STMT_START \
  { if (!*(d)) pari_err(primer1); NEXT_PRIME_VIADIFF(p,d); } STMT_END
 
/* For use with pari_init_opts */
 
#define INIT_JMPm 1
#define INIT_SIGm 2
#define INIT_DFTm 4
