#line 2 "../src/kernel/none/mp_indep.c"
/* $Id: mp_indep.c 12097 2010-01-29 11:57:49Z bill $

Copyright (C) 2000  The PARI group.

This file is part of the PARI/GP package.

PARI/GP is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation. It is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY WHATSOEVER.

Check the License for details. You should have received a copy of it, along
with the package; see the file 'COPYING'. If not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA. */

/* Find c such that 1=c*b mod 2^BITS_IN_LONG, assuming b odd (unchecked) */
ulong
invrev(ulong b)
{
  static int tab[] = { 0, 0, 0, 8, 0, 8, 0, 0 };
  ulong x = b + tab[b & 7]; /* b^(-1) mod 2^4 */

  /* Newton applied to 1/x - b = 0 */
#ifdef LONG_IS_64BIT
  x = x*(2-b*x); /* one more pass necessary */
#endif
  x = x*(2-b*x);
  x = x*(2-b*x); return x*(2-b*x);
}

void
affrr(GEN x, GEN y)
{
  long lx,ly,i;

  y[1] = x[1]; if (!signe(x)) return;

  lx=lg(x); ly=lg(y);
  if (lx <= ly)
  {
    for (i=2; i<lx; i++) y[i]=x[i];
    for (   ; i<ly; i++) y[i]=0;
    return;
  }
  for (i=2; i<ly; i++) y[i]=x[i];
  /* lx > ly: round properly */
  if (x[ly] & HIGHBIT) roundr_up_ip(y, ly);
}

GEN
ishiftr(GEN x, long s)
{
  long ex,lx,n;
  if (!signe(x)) return gen_0;
  ex = expo(x) + s; if (ex < 0) return gen_0;
  lx = lg(x);
  n=ex - bit_accuracy(lx) + 1;
  return ishiftr_lg(x, lx, n);
}

GEN
icopy_spec(GEN x, long nx)
{
  GEN z;
  long i;
  if (!nx) return gen_0;
  z = cgeti(nx+2); z[1] = evalsigne(1)|evallgefint(nx+2);
  for (i=0; i<nx; i++) z[i+2] = x[i];
  return z;
}

GEN
mului(ulong x, GEN y)
{
  long s = signe(y);
  GEN z;

  if (!s || !x) return gen_0;
  z = muluispec(x, y+2, lgefint(y)-2);
  setsigne(z,s); return z;
}

GEN
mulsi(long x, GEN y)
{
  long s = signe(y);
  GEN z;

  if (!s || !x) return gen_0;
  if (x<0) { s = -s; x = -x; }
  z = muluispec((ulong)x, y+2, lgefint(y)-2);
  setsigne(z,s); return z;
}

/* assume x > 1, y != 0. Return u * y with sign s */
static GEN
mulur_2(ulong x, GEN y, long s)
{
  long m, sh, i, lx = lg(y), e = expo(y);
  GEN z = cgetr(lx);
  ulong garde;
  LOCAL_HIREMAINDER;

  y--; garde = mulll(x,y[lx]);
  for (i=lx-1; i>=3; i--) z[i]=addmul(x,y[i]);
  z[2]=hiremainder; /* != 0 since y normalized and |x| > 1 */

  sh = bfffo(hiremainder); m = BITS_IN_LONG-sh;
  if (sh) shift_left2(z,z, 2,lx-1, garde,sh,m);
  z[1] = evalsigne(s) | evalexpo(m+e); return z;
}

GEN
mulsr(long x, GEN y)
{
  long s, e;

  if (!x) return gen_0;
  s = signe(y);
  if (!s)
  {
    if (x < 0) x = -x;
    e = expo(y) + (BITS_IN_LONG-1)-bfffo(x);
    return real_0_bit(e);
  }
  if (x==1)  return rcopy(y);
  if (x==-1) return negr(y);
  if (x < 0)
    return mulur_2((ulong)-x, y, -s);
  else
    return mulur_2((ulong)x, y, s);
}

GEN
mulur(ulong x, GEN y)
{
  long s, e;

  if (!x) return gen_0;
  s = signe(y);
  if (!s)
  {
    e = expo(y) + (BITS_IN_LONG-1)-bfffo(x);
    return real_0_bit(e);
  }
  if (x==1) return rcopy(y);
  return mulur_2(x, y, s);
}

#ifdef KARAMULR_VARIANT
static GEN addshiftw(GEN x, GEN y, long d);
static GEN
karamulrr1(GEN y, GEN x, long ly, long lz)
{
  long i, l, lz2 = (lz+2)>>1, lz3 = lz-lz2;
  GEN lo1, lo2, hi;

  hi = muliispec_mirror(x,y, lz2,lz2);
  i = lz2; while (i<lz && !x[i]) i++;
  lo1 = muliispec_mirror(y,x+i, lz2,lz-i);
  i = lz2; while (i<ly && !y[i]) i++;
  lo2 = muliispec_mirror(x,y+i, lz2,ly-i);
  if (signe(lo1))
  {
    if (ly!=lz) { lo2 = addshiftw(lo1,lo2,1); lz3++; }
    else lo2 = addii(lo1,lo2);
  }
  l = lgefint(lo2)-(lz3+2);
  if (l > 0) hi = addiispec(hi+2,lo2+2, lgefint(hi)-2,l);
  return hi;
}
#endif

/* set z <-- x*y, floating point multiplication.
 * lz = lg(z) = lg(x) <= ly <= lg(y), sz = signe(z). flag = lg(x) < lg(y) */
INLINE void
mulrrz_i(GEN z, GEN x, GEN y, long lz, long flag, long sz)
{
  long ez = expo(x) + expo(y);
  long i, j, lzz, p1;
  ulong garde;
  GEN y1;
  LOCAL_HIREMAINDER;
  LOCAL_OVERFLOW;

  if (lz > KARATSUBA_MULR_LIMIT) 
  {
    pari_sp av = avma;
#ifdef KARAMULR_VARIANT
    GEN hi = karamulrr1(y+2, x+2, lz+flag-2, lz-2); 
#else
    GEN hi = muliispec_mirror(y+2, x+2, lz+flag-2, lz-2);
#endif
    garde = hi[lz];
    if (hi[2] < 0)
    {
      ez++;
      for (i=2; i<lz ; i++) z[i] = hi[i];
    }
    else
    {
      shift_left(z,hi,2,lz-1, garde, 1);
      garde <<= 1;
    }
    if (garde & HIGHBIT)
    { /* round to nearest */
      i = lz; do z[--i]++; while (z[i]==0 && i > 1);
      if (i == 1) { z[2] = HIGHBIT; ez++; }
    }
    z[1] = evalsigne(sz)|evalexpo(ez);
#if 0
{
GEN U;
KARATSUBA_MULR_LIMIT = 100000;
U = mulrr(x, y);
KARATSUBA_MULR_LIMIT = 4;
if (!gequal(U, z)) pari_err(talker,"toto");
}
#endif
    avma = av; return;
  }
  if (lz == 3)
  {
    if (flag)
    {
      (void)mulll(x[2],y[3]);
      garde = addmul(x[2],y[2]);
    }
    else
      garde = mulll(x[2],y[2]);
    if (hiremainder & HIGHBIT)
    {
      ez++;
      /* hiremainder < (2^BIL-1)^2 / 2^BIL, hence hiremainder+1 != 0 */
      if (garde & HIGHBIT) hiremainder++; /* round properly */
    }
    else
    {
      hiremainder = (hiremainder<<1) | (garde>>(BITS_IN_LONG-1));
      if (garde & (HIGHBIT-1))
      {
        hiremainder++; /* round properly */
        if (!hiremainder) { hiremainder = HIGHBIT; ez++; }
      }
    }
    z[1] = evalsigne(sz) | evalexpo(ez);
    z[2] = hiremainder; return;
  }

  if (flag) { (void)mulll(x[2],y[lz]); garde = hiremainder; } else garde = 0;
  lzz=lz-1; p1=x[lzz];
  if (p1)
  {
    (void)mulll(p1,y[3]);
    garde = addll(addmul(p1,y[2]), garde);
    z[lzz] = overflow+hiremainder;
  }
  else z[lzz]=0;
  for (j=lz-2, y1=y-j; j>=3; j--)
  {
    p1 = x[j]; y1++;
    if (p1)
    {
      (void)mulll(p1,y1[lz+1]);
      garde = addll(addmul(p1,y1[lz]), garde);
      for (i=lzz; i>j; i--)
      {
        hiremainder += overflow;
	z[i] = addll(addmul(p1,y1[i]), z[i]);
      }
      z[j] = hiremainder+overflow;
    }
    else z[j]=0;
  }
  p1 = x[2]; y1++;
  garde = addll(mulll(p1,y1[lz]), garde);
  for (i=lzz; i>2; i--)
  {
    hiremainder += overflow;
    z[i] = addll(addmul(p1,y1[i]), z[i]);
  }
  z[2] = hiremainder+overflow;

  if (z[2] < 0)
    ez++;
  else
  {
    shift_left(z,z,2,lzz, garde, 1);
    garde <<= 1;
  }
  if (garde & HIGHBIT)
  { /* round to nearest */
    i = lz; do z[--i]++; while (z[i]==0 && i > 2);
    if (z[i] == 0) { z[2] = (long)HIGHBIT; ez++; }
  }
  z[1] = evalsigne(sz) | evalexpo(ez);
}

GEN
mulrr(GEN x, GEN y)
{
  long flag, ly, lz, sx = signe(x), sy = signe(y);
  GEN z;

  if (!sx || !sy) return real_0_bit(expo(x) + expo(y));
  if (sy < 0) sx = -sx;
  lz = lg(x);
  ly = lg(y);
  if (lz > ly) { lz = ly; swap(x, y); flag = 1; } else flag = (lz != ly);
  z = cgetr(lz);
  mulrrz_i(z, x,y, lz,flag, sx);
  return z;
}

GEN
mulir(GEN x, GEN y)
{
  long sx = signe(x), sy, lz;
  GEN z;

  if (!sx) return gen_0;
  if (!is_bigint(x)) return mulsr(itos(x), y);
  sy = signe(y);
  if (!sy) return real_0_bit(expi(x) + expo(y));
  if (sy < 0) sx = -sx;
  lz = lg(y); z = cgetr(lz);
  mulrrz_i(z, itor(x,lz),y, lz,0, sx);
  avma = (pari_sp)z; return z;
}

/* written by Bruno Haible following an idea of Robert Harley */
long
vals(ulong z)
{
  static char tab[64]={-1,0,1,12,2,6,-1,13,3,-1,7,-1,-1,-1,-1,14,10,4,-1,-1,8,-1,-1,25,-1,-1,-1,-1,-1,21,27,15,31,11,5,-1,-1,-1,-1,-1,9,-1,-1,24,-1,-1,20,26,30,-1,-1,-1,-1,23,-1,19,29,-1,22,18,28,17,16,-1};
#ifdef LONG_IS_64BIT
  long s;
#endif

  if (!z) return -1;
#ifdef LONG_IS_64BIT
  if (! (z&0xffffffff)) { s = 32; z >>=32; } else s = 0;
#endif
  z |= ~z + 1;
  z += z << 4;
  z += z << 6;
  z ^= z << 16; /* or  z -= z<<16 */
#ifdef LONG_IS_64BIT
  return s + tab[(z&0xffffffff)>>26];
#else
  return tab[z>>26];
#endif
}

GEN
divsi(long x, GEN y)
{
  long p1, s = signe(y);
  LOCAL_HIREMAINDER;

  if (!s) pari_err(gdiver);
  if (!x || lgefint(y)>3 || ((long)y[2])<0) return gen_0;
  hiremainder=0; p1=divll(labs(x),y[2]);
  if (x<0) { hiremainder = -((long)hiremainder); p1 = -p1; }
  if (s<0) p1 = -p1;
  return stoi(p1);
}

GEN
divir(GEN x, GEN y)
{
  GEN z;
  long ly;
  pari_sp av;

  if (!signe(y)) pari_err(gdiver);
  if (!signe(x)) return gen_0;
  ly = lg(y); z = cgetr(ly); av = avma; 
  affrr(divrr(itor(x, ly+1), y), z);
  avma = av; return z;
}

void
mpdivz(GEN x, GEN y, GEN z)
{
  pari_sp av = avma;
  long tx = typ(x), ty = typ(y);
  GEN q = tx==t_INT && ty==t_INT? divii(x,y): mpdiv(x,y);

  if (typ(z) == t_REAL) affrr(q, z);
  else
  {
    if (typ(q) == t_REAL) pari_err(gdiver);
    affii(q, z); 
  }
  avma = av;
}

GEN
divsr(long x, GEN y)
{
  pari_sp av;
  long ly;
  GEN z;

  if (!signe(y)) pari_err(gdiver);
  if (!x) return gen_0;
  ly = lg(y); z = cgetr(ly); av = avma;
  affrr(divrr(stor(x,ly+1), y), z);
  avma = av; return z;
}

GEN
modii(GEN x, GEN y)
{
  switch(signe(x))
  {
    case 0: return gen_0;
    case 1: return remii(x,y);
    default:
    {
      pari_sp av = avma;
      (void)new_chunk(lgefint(y));
      x = remii(x,y); avma=av;
      if (x==gen_0) return x;
      return subiispec(y+2,x+2,lgefint(y)-2,lgefint(x)-2);
    }
  }
}

void
modiiz(GEN x, GEN y, GEN z)
{
  const pari_sp av = avma;
  affii(modii(x,y),z); avma=av;
}

GEN
divrs(GEN x, long y)
{
  long i,lx,garde,sh,s=signe(x);
  GEN z;
  LOCAL_HIREMAINDER;

  if (!y) pari_err(gdiver);
  if (!s) return real_0_bit(expo(x) - (BITS_IN_LONG-1)+bfffo(y));
  if (y<0) { s = -s; y = -y; }
  if (y==1) { z=rcopy(x); setsigne(z,s); return z; }

  z=cgetr(lx=lg(x)); hiremainder=0;
  for (i=2; i<lx; i++) z[i] = divll(x[i],y);

  /* we may have hiremainder != 0 ==> garde */
  garde=divll(0,y); sh=bfffo(z[2]);
  if (sh) shift_left(z,z, 2,lx-1, garde,sh);
  z[1] = evalsigne(s) | evalexpo(expo(x)-sh);
  if ((garde << sh) & HIGHBIT) roundr_up_ip(z, lx);
  return z;
}

GEN
truedvmdii(GEN x, GEN y, GEN *z)
{
  pari_sp av;
  GEN r, q, *gptr[2];
  if (!is_bigint(y)) return truedvmdis(x, itos(y), z);

  av = avma;
  q = dvmdii(x,y,&r); /* assume that r is last on stack */

  if (signe(r)>=0)
  {
    if (z == ONLY_REM) return gerepileuptoint(av,r);
    if (z) *z = r; else cgiv(r);
    return q;
  }

  if (z == ONLY_REM)
  { /* r += sign(y) * y, that is |y| */
    r = subiispec(y+2,r+2, lgefint(y)-2,lgefint(r)-2);
    return gerepileuptoint(av, r);
  }
  q = addis(q, -signe(y));
  if (!z) return gerepileuptoint(av, q);

  *z = subiispec(y+2,r+2, lgefint(y)-2,lgefint(r)-2);
  gptr[0]=&q; gptr[1]=z; gerepilemanysp(av,(pari_sp)r,gptr,2);
  return q;
}
GEN
truedvmdis(GEN x, long y, GEN *z)
{
  pari_sp av = avma;
  long r;
  GEN q = divis_rem(x,y,&r);

  if (r >= 0)
  {
    if (z == ONLY_REM) { avma = av; return utoi(r); }
    if (z) *z = utoi(r);
    return q;
  }
  if (z == ONLY_REM) { avma = av; return utoi(r + labs(y)); }
  q = gerepileuptoint(av, addis(q, (y < 0)? 1: -1));
  if (z) *z = utoi(r + labs(y));
  return q;
}

/* 2^n = shifti(gen_1, n) */
GEN
int2n(long n) {
  long i, m, d, l;
  GEN z;
  if (n < 0) return gen_0;
  if (n == 0) return gen_1;

  d = n>>TWOPOTBITS_IN_LONG;
  m = n & (BITS_IN_LONG-1);
  l = d + 3; z = cgeti(l);
  z[1] = evalsigne(1) | evallgefint(l);
  for (i = 2; i < l; i++) z[i] = 0;
  *int_MSW(z) = 1L << m; return z;
}
/* To avoid problems when 2^(BIL-1) < n. Overflow cleanly, where int2n
 * returns gen_0 */
GEN
int2u(ulong n) {
  ulong i, m, d, l;
  GEN z;
  if (n == 0) return gen_1;

  d = n>>TWOPOTBITS_IN_LONG;
  m = n & (BITS_IN_LONG-1);
  l = d + 3; z = cgeti(l);
  z[1] = evalsigne(1) | evallgefint(l);
  for (i = 2; i < l; i++) z[i] = 0;
  *int_MSW(z) = 1L << m; return z;
}

/* actual operations will take place on a+2 and b+2: we strip the codewords */
GEN
mulii(GEN a,GEN b)
{
  long sa,sb;
  GEN z;

  sa=signe(a); if (!sa) return gen_0;
  sb=signe(b); if (!sb) return gen_0;
  if (sb<0) sa = -sa;
  z = muliispec(a+2,b+2, lgefint(a)-2,lgefint(b)-2);
  setsigne(z,sa); return z;
}

GEN
remiimul(GEN x, GEN sy)
{
  GEN r, q, y = (GEN)sy[1], invy;
  long k;
  pari_sp av = avma;

  k = cmpii(x, y);
  if (k <= 0) return k? icopy(x): gen_0;
  invy = (GEN)sy[2];
  q = mulir(x,invy);
  q = truncr(q); /* differs from divii(x, y) at most by 1 */
  r = subii(x, mulii(y,q));
  if (signe(r) < 0)
    r = subiispec(y+2,r+2, lgefint(y)-2, lgefint(r)-2); /* y - |r| */
  else
  {
    /* remii(x,y) + y >= r >= remii(x,y) */
    k = absi_cmp(r, y);
    if (k >= 0)
    {
      if (k == 0) { avma = av; return gen_0; }
      r = subiispec(r+2, y+2, lgefint(r)-2, lgefint(y)-2);
    }
  }
#if 0
  q = subii(r,remii(x,y));
  if (signe(q))
    pari_err(talker,"bug in remiimul: x = %Z\ny = %Z\ndif = %Z", x,y,q);
#endif
  return gerepileuptoint(av, r); /* = remii(x, y) */
}

GEN
sqri(GEN a) { return sqrispec(a+2, lgefint(a)-2); }

/* Old cgiv without reference count (which was not used anyway)
 * Should be a macro.
 */
void
cgiv(GEN x)
{
  if (x == (GEN) avma)
    avma = (pari_sp) (x+lg(x));
}

/* sqrt()'s result may be off by 1 when a is not representable exactly as a
 * double [64bit machine] */
ulong
usqrtsafe(ulong a)
{
  ulong x = (ulong)sqrt((double)a);
#ifdef LONG_IS_64BIT
  if (x > LOWMASK || x*x > a) x--;
#endif
  return x;
}

/********************************************************************/
/**                                                                **/
/**                         RANDOM INTEGERS                        **/
/**                                                                **/
/********************************************************************/
static long pari_randseed = 1;

/* BSD rand gives this: seed = 1103515245*seed + 12345 */
/*Return 31 ``random'' bits.*/
long
pari_rand31(void)
{
  pari_randseed = (1000276549*pari_randseed + 12347) & 0x7fffffff;
  return pari_randseed;
}

long
setrand(long seed) { return (pari_randseed = seed); }

long
getrand(void) { return pari_randseed; }

static ulong
pari_rand(void)
{
#define GLUE2(hi, lo) (((hi) << BITS_IN_HALFULONG) | (lo))
#if !defined(LONG_IS_64BIT)
  return GLUE2((pari_rand31()>>12)&LOWMASK,
               (pari_rand31()>>12)&LOWMASK);
#else
#define GLUE4(hi1,hi2, lo1,lo2) GLUE2(((hi1)<<16)|(hi2), ((lo1)<<16)|(lo2))
#  define LOWMASK2 0xffffUL
  return GLUE4((pari_rand31()>>12)&LOWMASK2,
               (pari_rand31()>>12)&LOWMASK2,
               (pari_rand31()>>12)&LOWMASK2,
               (pari_rand31()>>12)&LOWMASK2);
#endif
}

#if 1
/* assume N > 0 */
GEN
randomi(GEN N)
{
  ulong n;
  long lx, i;
  GEN x, xMSW, xd, Nd;

  lx = lgefint(N); x = cgeti(lx);
  x[1] = evalsigne(1) | evallgefint(lx);
  xMSW = xd = int_MSW(x);
  for (i=2; i<lx; i++) { *xd = pari_rand(); xd = int_precW(xd); }

  Nd = int_MSW(N); n = (ulong)*Nd;
  xd = xMSW;
  if (lx == 3) n--;
  else
    for (i=3; i<lx; i++)
    {
      xd = int_precW(xd);
      Nd = int_precW(Nd);
      if (*xd != *Nd) {
        if ((ulong)*xd > (ulong)*Nd) n--;
        break;
      }
    }
  /* MSW needs to be generated between 0 and n */
  if (n) {
    LOCAL_HIREMAINDER;
    (void)mulll((ulong)*xMSW, n + 1);
    n = hiremainder;
  }
  *xMSW = (long)n;
  if (!n) x = int_normalize(x, 1);
  return x;
}
#else
/* assume N > 0 */
GEN
randomi(GEN N)
{
  long i, lx = lgefint(N), shift = bfffo(*int_MSW(N));
  GEN x = cgeti(lx), xMSW, xd;

  x[1] = evalsigne(1) | evallgefint(lx);
  xMSW = int_MSW(x);
  for (xd = xMSW;;) {
    for (i=2; i<lx; i++) { *xd = pari_rand(); xd = int_precW(xd); }
    *xMSW = ((ulong)*xMSW) >> shift;
    x = int_normalize(x, 0);
    if (absi_cmp(x, N) < 0) return x;
  }
}
#endif

/********************************************************************/
/**                                                                **/
/**              EXPONENT / CONVERSION t_REAL --> double           **/
/**                                                                **/
/********************************************************************/

#ifdef LONG_IS_64BIT
long
dblexpo(double x)
{
  union { double f; ulong i; } fi;
  const int mant_len = 52;  /* mantissa bits (excl. hidden bit) */
  const int exp_mid = 0x3ff;/* exponent bias */

  if (x==0.) return -exp_mid;
  fi.f = x;
  return ((fi.i & (HIGHBIT-1)) >> mant_len) - exp_mid;
}

ulong
dblmantissa(double x)
{
  union { double f; ulong i; } fi;
  const int expo_len = 11; /* number of bits of exponent */

  if (x==0.) return 0;
  fi.f = x;
  return (fi.i << expo_len) | HIGHBIT;
}

GEN
dbltor(double x)
{
  GEN z;
  long e;
  union { double f; ulong i; } fi;
  const int mant_len = 52;  /* mantissa bits (excl. hidden bit) */
  const int exp_mid = 0x3ff;/* exponent bias */
  const int expo_len = 11; /* number of bits of exponent */

  if (x==0.) return real_0_bit(-exp_mid);
  fi.f = x; z = cgetr(DEFAULTPREC);
  {
    const ulong a = fi.i;
    ulong A;
    e = ((a & (HIGHBIT-1)) >> mant_len) - exp_mid;
    if (e == exp_mid+1) pari_err(talker, "NaN or Infinity in dbltor");
    A = a << expo_len;
    if (e == -exp_mid)
    { /* unnormalized values */
      int sh = bfffo(A);
      e -= sh-1;
      z[2] = A << sh;
    }
    else
      z[2] = HIGHBIT | A;
    z[1] = evalexpo(e) | evalsigne(x<0? -1: 1);
  }
  return z;
}

double
rtodbl(GEN x)
{
  long ex,s=signe(x);
  ulong a;
  union { double f; ulong i; } fi;
  const int mant_len = 52;  /* mantissa bits (excl. hidden bit) */
  const int exp_mid = 0x3ff;/* exponent bias */
  const int expo_len = 11; /* number of bits of exponent */

  if (typ(x)==t_INT && !s) return 0.0;
  if (typ(x)!=t_REAL) pari_err(typeer,"rtodbl");
  if (!s || (ex=expo(x)) < - exp_mid) return 0.0;

  /* start by rounding to closest */
  a = (x[2] & (HIGHBIT-1)) + 0x400;
  if (a & HIGHBIT) { ex++; a=0; }
  if (ex >= exp_mid) pari_err(rtodber);
  fi.i = ((ex + exp_mid) << mant_len) | (a >> expo_len);
  if (s<0) fi.i |= HIGHBIT;
  return fi.f;
}

#else /* LONG_IS_64BIT */

#if   PARI_DOUBLE_FORMAT == 1
#  define INDEX0 1
#  define INDEX1 0
#elif PARI_DOUBLE_FORMAT == 0
#  define INDEX0 0
#  define INDEX1 1
#endif

long
dblexpo(double x)
{
  union { double f; ulong i[2]; } fi;
  const int mant_len = 52;  /* mantissa bits (excl. hidden bit) */
  const int exp_mid = 0x3ff;/* exponent bias */
  const int shift = mant_len-32;

  if (x==0.) return -exp_mid;
  fi.f = x;
  {
    const ulong a = fi.i[INDEX0];
    return ((a & (HIGHBIT-1)) >> shift) - exp_mid;
  }
}

ulong
dblmantissa(double x)
{
  union { double f; ulong i[2]; } fi;
  const int expo_len = 11; /* number of bits of exponent */

  if (x==0.) return 0;
  fi.f = x;
  {
    const ulong a = fi.i[INDEX0];
    const ulong b = fi.i[INDEX1];
    return HIGHBIT | b >> (BITS_IN_LONG-expo_len) | (a << expo_len);
  }
}

GEN
dbltor(double x)
{
  GEN z;
  long e;
  union { double f; ulong i[2]; } fi;
  const int mant_len = 52;  /* mantissa bits (excl. hidden bit) */
  const int exp_mid = 0x3ff;/* exponent bias */
  const int expo_len = 11; /* number of bits of exponent */
  const int shift = mant_len-32;

  if (x==0.) return real_0_bit(-exp_mid);
  fi.f = x; z = cgetr(DEFAULTPREC);
  {
    const ulong a = fi.i[INDEX0];
    const ulong b = fi.i[INDEX1];
    ulong A, B;
    e = ((a & (HIGHBIT-1)) >> shift) - exp_mid;
    if (e == exp_mid+1) pari_err(talker, "NaN or Infinity in dbltor");
    A = b >> (BITS_IN_LONG-expo_len) | (a << expo_len);
    B = b << expo_len;
    if (e == -exp_mid)
    { /* unnormalized values */
      int sh;
      if (A)
      {
        sh = bfffo(A);
        e -= sh-1;
        z[2] = (A << sh) | (B >> (32-sh));
        z[3] = B << sh;
      }
      else 
      {
        sh = bfffo(B); /* B != 0 */
        e -= sh-1 + 32;
        z[2] = B << sh;
        z[3] = 0;
      }
    }
    else
    {
      z[3] = B;
      z[2] = HIGHBIT | A;
    }
    z[1] = evalexpo(e) | evalsigne(x<0? -1: 1);
  }
  return z;
}

double
rtodbl(GEN x)
{
  long ex,s=signe(x),lx=lg(x);
  ulong a,b,k;
  union { double f; ulong i[2]; } fi;
  const int mant_len = 52;  /* mantissa bits (excl. hidden bit) */
  const int exp_mid = 0x3ff;/* exponent bias */
  const int expo_len = 11; /* number of bits of exponent */
  const int shift = mant_len-32;

  if (typ(x)==t_INT && !s) return 0.0;
  if (typ(x)!=t_REAL) pari_err(typeer,"rtodbl");
  if (!s || (ex=expo(x)) < - exp_mid) return 0.0;

  /* start by rounding to closest */
  a = x[2] & (HIGHBIT-1);
  if (lx > 3)
  {
    b = x[3] + 0x400UL; if (b < 0x400UL) a++;
    if (a & HIGHBIT) { ex++; a=0; }
  }
  else b = 0;
  if (ex >= exp_mid) pari_err(rtodber);
  ex += exp_mid;
  k = (a >> expo_len) | (ex << shift);
  if (s<0) k |= HIGHBIT;
  fi.i[INDEX0] = k;
  fi.i[INDEX1] = (a << (BITS_IN_LONG-expo_len)) | (b >> expo_len);
  return fi.f;
}
#endif /* LONG_IS_64BIT */

