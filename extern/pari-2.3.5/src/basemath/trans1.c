/* $Id: trans1.c 12114 2010-02-03 22:59:09Z bill $

Copyright (C) 2000  The PARI group.

This file is part of the PARI/GP package.

PARI/GP is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation. It is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY WHATSOEVER.

Check the License for details. You should have received a copy of it, along
with the package; see the file 'COPYING'. If not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA. */

/********************************************************************/
/**                                                                **/
/**                   TRANSCENDENTAL FUNCTIONS                     **/
/**                                                                **/
/********************************************************************/
#include "pari.h"
#include "paripriv.h"

#ifdef LONG_IS_64BIT
# define SQRTVERYBIGINT 3037000500   /* ceil(sqrt(VERYBIGINT)) */
# define CBRTVERYBIGINT 2097152      /* ceil(cbrt(VERYBIGINT)) */
#else
# define SQRTVERYBIGINT 46341
# define CBRTVERYBIGINT  1291
#endif

static GEN glog2;
void
pari_init_floats(void)
{
  geuler = gpi = bernzone = glog2 = NULL;
}

/********************************************************************/
/**                                                                **/
/**                               PI                               **/
/**                                                                **/
/********************************************************************/
#if 0
/* Ramanujan's formula:
 *                         ----
 *  53360 (640320)^(1/2)   \    (6n)! (545140134 n + 13591409)
 *  -------------------- = /    ------------------------------
 *        Pi               ----   (n!)^3 (3n)! (-640320)^(3n)
 *                         n>=0
 */
GEN
piold(long prec)
{
  const long k1 = 545140134, k2 = 13591409, k3 = 640320;
  const double alpha2 = 47.11041314/BITS_IN_LONG; /* 3log(k3/12) / log(2^BIL) */
  GEN p1,p2,p3,tmppi;
  long l, n, n1;
  pari_sp av = avma, av2;
  double alpha;

  tmppi = cgetr(prec);
  prec++;
  n = (long)(1 + (prec-2)/alpha2);
  n1 = 6*n - 1;
  p2 = addsi(k2, mulss(n,k1));
  p1 = itor(p2, prec);

  /* initialize mantissa length */
  if (prec>=4) l=4; else l=prec;
  setlg(p1,l); alpha = (double)l;

  av2 = avma;
  while (n)
  {
    if (n < CBRTVERYBIGINT) /* p3 = n1(n1-2)(n1-4) / n^3 */
      p3 = divrs(mulsr(n1-4,mulsr(n1*(n1-2),p1)),n*n*n);
    else
    {
      if (n1 < SQRTVERYBIGINT)
	p3 = divrs(divrs(mulsr(n1-4,mulsr(n1*(n1-2),p1)),n*n),n);
      else
	p3 = divrs(divrs(divrs(mulsr(n1-4,mulsr(n1,mulsr(n1-2,p1))),n),n),n);
    }
    p3 = divrs(divrs(p3,100100025), 327843840);
    subisz(p2,k1,p2);
    subirz(p2,p3,p1);
    alpha += alpha2; l = (long)(1+alpha); /* new mantissa length */
    if (l > prec) l = prec;
    setlg(p1,l); avma = av2;
    n--; n1-=6;
  }
  p1 = divsr(53360,p1);
  return gerepileuptoleaf(av, mulrr(p1,sqrtr_abs(stor(k3,prec))));
}
#endif
/* Gauss - Brent-Salamin AGM iteration */
void
constpi(long prec)
{
  GEN A, B, C, tmppi;
  long i, G;
  pari_sp av, av2;

  if (gpi && lg(gpi) >= prec) return;

  av = avma; tmppi = newbloc(prec);
  *tmppi = evaltyp(t_REAL) | evallg(prec);
  G = - bit_accuracy(prec);
  prec++;

  A = real_1(prec);
  B = sqrtr_abs(real2n(1,prec)); setexpo(B, -1); /* = 1/sqrt(2) */
  C = real2n(-2, prec); av2 = avma;
  for (i = 0;; i++)
  {
    GEN y, a, b, B_A = subrr(B, A);
    if (expo(B_A) < G) break;
    a = addrr(A,B); setexpo(a, expo(a)-1);
    b = sqrtr_abs( mulrr(A, B) );
    y = gsqr(B_A); setexpo(y, expo(y) + i - 2);
    affrr(subrr(C, y), C);
    affrr(a, A);
    affrr(b, B); avma = av2;
  }
  setexpo(C, expo(C)+2);
  affrr(divrr(gsqr(addrr(A,B)), C), tmppi);
  if (gpi) gunclone(gpi);
  avma = av;  gpi = tmppi;
}

GEN
mppi(long prec)
{
  GEN x = cgetr(prec);
  constpi(prec); affrr(gpi,x); return x;
}

/* Pi * 2^n */
GEN
Pi2n(long n, long prec)
{
  GEN x = mppi(prec); setexpo(x, 1+n);
  return x;
}

/* I * Pi * 2^n */
GEN
PiI2n(long n, long prec)
{
  GEN z = cgetg(3,t_COMPLEX);
  gel(z,1) = gen_0;
  gel(z,2) = Pi2n(n, prec); return z;
}

/* 2I * Pi */
GEN
PiI2(long prec) { return PiI2n(1, prec); }

/********************************************************************/
/**                                                                **/
/**                       EULER CONSTANT                           **/
/**                                                                **/
/********************************************************************/

void
consteuler(long prec)
{
  GEN u,v,a,b,tmpeuler;
  long l, n1, n, k, x;
  pari_sp av1, av2;

  if (geuler && lg(geuler) >= prec) return;

  av1 = avma; tmpeuler = newbloc(prec);
  *tmpeuler = evaltyp(t_REAL) | evallg(prec);

  prec++;

  l = prec+1; x = (long) (1 + bit_accuracy_mul(l, LOG2/4));
  a = stor(x,l); u=logr_abs(a); setsigne(u,-1); affrr(u,a);
  b = real_1(l);
  v = real_1(l);
  n = (long)(1+3.591*x); /* z=3.591: z*[ ln(z)-1 ]=1 */
  n1 = min(n, SQRTVERYBIGINT);
  if (x < SQRTVERYBIGINT)
  {
    long xx = x*x;
    av2 = avma;
    for (k=1; k<n1; k++)
    {
      divrsz(mulsr(xx,b),k*k, b);
      divrsz(addrr(divrs(mulsr(xx,a),k),b),k, a);
      addrrz(u,a,u);
      addrrz(v,b,v); avma = av2;
    }
    for (   ; k<=n; k++)
    {
      divrsz(divrs(mulsr(xx,b),k), k, b);
      divrsz(addrr(divrs(mulsr(xx,a),k),b),k, a);
      addrrz(u,a,u);
      addrrz(v,b,v); avma = av2;
    }
  }
  else
  {
    GEN xx = mulss(x,x);
    av2 = avma;
    for (k=1; k<n1; k++)
    {
      divrsz(mulir(xx,b),k*k, b);
      divrsz(addrr(divrs(mulir(xx,a),k),b),k, a);
      addrrz(u,a,u);
      addrrz(v,b,v); avma = av2;
    }
    for (   ; k<=n; k++)
    {
      divrsz(divrs(mulir(xx,b),k), k, b);
      divrsz(addrr(divrs(mulir(xx,a),k),b),k, a);
      addrrz(u,a,u);
      addrrz(v,b,v); avma = av2;
    }
  }
  divrrz(u,v,tmpeuler);
  if (geuler) gunclone(geuler);
  avma = av1; geuler = tmpeuler;
}

GEN
mpeuler(long prec)
{
  GEN x = cgetr(prec);
  consteuler(prec); affrr(geuler,x); return x;
}

/********************************************************************/
/**                                                                **/
/**          TYPE CONVERSION FOR TRANSCENDENTAL FUNCTIONS          **/
/**                                                                **/
/********************************************************************/

GEN
transc(GEN (*f)(GEN,long), GEN x, long prec)
{
  pari_sp tetpil, av = avma;
  GEN p1, y;
  long lx, i;

  if (prec < 2) pari_err(talker, "incorrect precision in transc");
  switch(typ(x))
  {
    case t_INT:
      p1 = itor(x, prec); tetpil=avma;
      return gerepile(av,tetpil,f(p1,prec));

    case t_FRAC:
      p1 = rdivii(gel(x,1), gel(x,2), prec); tetpil=avma;
      return gerepile(av,tetpil,f(p1,prec));

    case t_QUAD:
      p1 = quadtoc(x, prec); tetpil = avma;
      return gerepile(av,tetpil,f(p1,prec));

    case t_POL: case t_RFRAC:
      return gerepileupto(av, f(toser_i(x), prec));

    case t_VEC: case t_COL: case t_MAT:
      lx = lg(x); y = cgetg(lx,typ(x));
      for (i=1; i<lx; i++) gel(y,i) = f(gel(x,i),prec);
      return y;

    case t_POLMOD:
      p1 = cleanroots(gel(x,1),prec); lx = lg(p1);
      for (i=1; i<lx; i++) gel(p1,i) = poleval(gel(x,2),gel(p1,i));
      tetpil = avma; y = cgetg(lx,t_COL);
      for (i=1; i<lx; i++) gel(y,i) = f(gel(p1,i),prec);
      return gerepile(av,tetpil,y);

    default: pari_err(typeer,"a transcendental function");
  }
  return f(x,prec);
}

/*******************************************************************/
/*                                                                 */
/*                            POWERING                             */
/*                                                                 */
/*******************************************************************/
static GEN
puiss0(GEN x)
{
  long lx, i;
  GEN y;

  switch(typ(x))
  {
    case t_INT: case t_REAL: case t_FRAC: case t_COMPLEX:
    case t_PADIC: case t_QUAD:
      return gen_1;

    case t_INTMOD:
      y = cgetg(3,t_INTMOD); gel(y,1) = icopy(gel(x,1));
      gel(y,2) = gen_1; return y;

    case t_POLMOD:
      y = cgetg(3,t_POLMOD); gel(y,1) = gcopy(gel(x,1));
      gel(y,2) = pol_1[varn(x[1])]; return y;

    case t_POL: case t_SER: case t_RFRAC:
      return pol_1[gvar(x)];

    case t_MAT:
      lx=lg(x); if (lx==1) return cgetg(1,t_MAT);
      if (lx != lg(x[1])) pari_err(mattype1,"gpow");
      y = matid(lx-1);
      for (i=1; i<lx; i++) gcoeff(y,i,i) = puiss0(gcoeff(x,i,i));
      return y;
    case t_QFR: return qfr_unit(x);
    case t_QFI: return qfi_unit(x);
    case t_VECSMALL: return perm_identity(lg(x) - 1);
  }
  pari_err(typeer,"gpow");
  return NULL; /* not reached */
}

static GEN
_sqr(void *data /* ignored */, GEN x) { (void)data; return gsqr(x); }
static GEN
_mul(void *data /* ignored */, GEN x, GEN y) { (void)data; return gmul(x,y); }
static GEN
_sqri(void *data /* ignored */, GEN x) { (void)data; return sqri(x); }
static GEN
_muli(void *data /* ignored */, GEN x, GEN y) { (void)data; return mulii(x,y); }

/* INTEGER POWERING (a^n for integer a != 0 and integer n > 0)
 *
 * Use left shift binary algorithm (RS is wasteful: multiplies big numbers,
 * with LS one of them is the base, hence small). Sign of result is set
 * to s (= 1,-1). Makes life easier for caller, which otherwise might do a
 * setsigne(gen_1 / gen_m1) */
static GEN
powiu_sign(GEN a, ulong N, long s)
{
  pari_sp av;
  GEN y;

  if (lgefint(a) == 3)
  { /* easy if |a| < 3 */
    if (a[2] == 1) return (s>0)? gen_1: gen_m1;
    if (a[2] == 2) { a = int2u(N); setsigne(a,s); return a; }
  }
  if (N == 1) { a = icopy(a); setsigne(a,s); return a; }
  if (N == 2) return sqri(a);
  av = avma;
  y = leftright_pow_u(a, N, NULL, &_sqri, &_muli);
  setsigne(y,s); return gerepileuptoint(av, y);
}
/* a^N */
GEN
powiu(GEN a, ulong N)
{
  long s;
  if (!N) return gen_1;
  s = signe(a);
  if (!s) return gen_0;
  return powiu_sign(a, N, (s < 0 && odd(N))? -1: 1);
}
GEN
powuu(ulong p, ulong N)
{
  long P[] = {evaltyp(t_INT)|_evallg(3), evalsigne(1)|evallgefint(3),0};
  if (!N) return gen_1;
  if (!p) return gen_0;
  P[2] = p;
  return powiu_sign(P, N, 1);
}

/* assume p^k is SMALL */
ulong
upowuu(ulong p, ulong k)
{
  ulong i, pk;

  if (!k) return 1;
  if (p == 2) return 1UL<<k;
  pk = p; for (i=2; i<=k; i++) pk *= p;
  return pk;
}

typedef struct {
  long prec, a;
  GEN (*sqr)(GEN);
  GEN (*mulug)(ulong,GEN);
} sr_muldata;

static GEN
_rpowuu_mul(void *data, GEN x, GEN y/*unused*/)
{
  sr_muldata *D = (sr_muldata *)data;
  (void)y; return D->mulug(D->a, x);
}

static GEN
_rpowuu_sqr(void *data, GEN x)
{
  sr_muldata *D = (sr_muldata *)data;
  if (typ(x) == t_INT && lgefint(x) >= D->prec)
  { /* switch to t_REAL */
    D->sqr   = &gsqr;
    D->mulug = &mulur; x = itor(x, D->prec);
  }
  return D->sqr(x);
}

/* return a^n as a t_REAL of precision prec. Assume a > 0, n > 0 */
GEN
rpowuu(ulong a, ulong n, long prec)
{
  pari_sp av;
  GEN y;
  sr_muldata D;

  if (a == 1) return real_1(prec);
  if (a == 2) return real2n(n, prec);
  if (n == 1) return stor(a, prec);
  av = avma;
  D.sqr   = &sqri;
  D.mulug = &mului;
  D.prec = prec;
  D.a = (long)a;
  y = leftright_pow_u(utoipos(a), n, (void*)&D, &_rpowuu_sqr, &_rpowuu_mul);
  if (typ(y) == t_INT) y = itor(y, prec);
  return gerepileuptoleaf(av, y);
}

/* x^(s/2), assume x t_REAL */
GEN
powrshalf(GEN x, long s)
{
  if (s & 1) return sqrtr(gpowgs(x, s));
  return gpowgs(x, s>>1);
}
/* x^(n/d), assume x t_REAL, return t_REAL */
GEN
powrfrac(GEN x, long n, long d)
{
  long z;
  if (!n) return real_1(lg(x));
  z = cgcd(n, d); if (z > 1) { n /= z; d /= z; }
  if (d == 1) return gpowgs(x, n);
  x = gpowgs(x, n);
  if (d == 2) return sqrtr(x);
  return sqrtnr(x, d);
}

/* assume x != 0 */
static GEN
pow_monome(GEN x, long n)
{
  long i, d, dx = degpol(x);
  GEN A, b, y;

  if (n < 0) { n = -n; y = cgetg(3, t_RFRAC); } else y = NULL;

  if (HIGHWORD(dx) || HIGHWORD(n))
  {
    LOCAL_HIREMAINDER;
    d = (long)mulll((ulong)dx, (ulong)n);
    if (hiremainder || (d &~ LGBITS)) d = LGBITS; /* overflow */
    d += 2;
  }
  else
    d = dx*n + 2;
  if ((d + 1) & ~LGBITS) pari_err(talker,"degree overflow in pow_monome");
  A = cgetg(d+1, t_POL); A[1] = x[1];
  for (i=2; i < d; i++) gel(A,i) = gen_0;
  b = gpowgs(gel(x,dx+2), n); /* not memory clean if (n < 0) */
  if (!y) y = A;
  else {
    GEN c = denom(b);
    gel(y,1) = c; if (c != gen_1) b = gmul(b,c);
    gel(y,2) = A;
  }
  gel(A,d) = b; return y;
}

/* x t_PADIC */
static GEN
powps(GEN x, long n)
{
  long e = n*valp(x), v;
  GEN t, y, mod, p = gel(x,2);
  pari_sp av;

  if (!signe(x[4])) {
    if (n < 0) pari_err(gdiver);
    return zeropadic(p, e);
  }
  v = z_pval(n, p);

  y = cgetg(5,t_PADIC);
  mod = gel(x,3);
  if (v == 0) mod = icopy(mod);
  else
  {
    if (precp(x) == 1 && equaliu(p, 2)) v++;
    mod = mulii(mod, powiu(p,v));
    mod = gerepileuptoint((pari_sp)y, mod);
  }
  y[1] = evalprecp(precp(x) + v) | evalvalp(e);
  gel(y,2) = icopy(p);
  gel(y,3) = mod;

  av = avma; t = gel(x,4);
  if (n < 0) { t = Fp_inv(t, mod); n = -n; }
  t = Fp_powu(t, n, mod);
  gel(y,4) = gerepileuptoint(av, t);
  return y;
}
/* x t_PADIC */
static GEN
powp(GEN x, GEN n)
{
  long v;
  GEN y, mod, p = gel(x,2);

  if (valp(x)) pari_err(errlg);

  if (!signe(x[4])) {
    if (signe(n) < 0) pari_err(gdiver);
    return zeropadic(p, 0);
  }
  v = Z_pval(n, p);

  y = cgetg(5,t_PADIC);
  mod = gel(x,3);
  if (v == 0) mod = icopy(mod);
  else
  {
    mod = mulii(mod, powiu(p,v));
    mod = gerepileuptoint((pari_sp)y, mod);
  }
  y[1] = evalprecp(precp(x) + v) | evalvalp(0);
  gel(y,2) = icopy(p);
  gel(y,3) = mod;
  gel(y,4) = Fp_pow(gel(x,4), n, mod);
  return y;
}

GEN
gpowgs(GEN x, long n)
{
  long m;
  pari_sp av;
  GEN y;

  if (n == 0) return puiss0(x);
  if (n == 1) return gcopy(x);
  if (n ==-1) return ginv(x);
  switch(typ(x))
  {
    case t_INT:
    {
      long sx = signe(x), s;
      GEN t;
      if (!sx) {
        if (n < 0) pari_err(gdiver);
        return gen_0;
      }
      s = (sx < 0 && odd(n))? -1: 1;
      if (n > 0) return powiu_sign(x, n, s);
      t = (s > 0)? gen_1: gen_m1;
      if (is_pm1(x)) return t;
      /* n < 0, |x| > 1 */
      y = cgetg(3,t_FRAC);
      gel(y,1) = t;
      gel(y,2) = powiu_sign(x, -n, 1); /* force denominator > 0 */
      return y;
    }
    case t_INTMOD:
      y = cgetg(3,t_INTMOD); gel(y,1) = icopy(gel(x,1));
      gel(y,2) = Fp_pows(gel(x,2), n, gel(x,1));
      return y;
    case t_FRAC:
    {
      GEN a = gel(x,1), b = gel(x,2);
      long sx = signe(a), s;
      if (!sx) {
        if (n < 0) pari_err(gdiver);
        return gen_0;
      }
      s = (sx < 0 && odd(n))? -1: 1;
      if (n < 0) {
        n = -n;
        if (is_pm1(a)) return powiu_sign(b, n, s); /* +-1/x[2] inverts to t_INT */
        swap(a, b);
      }
      y = cgetg(3, t_FRAC);
      gel(y,1) = powiu_sign(a, n, s);
      gel(y,2) = powiu_sign(b, n, 1);
      return y;
    }
    case t_PADIC: return powps(x, n);
    case t_RFRAC:
    {
      av = avma; y = cgetg(3, t_RFRAC); m = labs(n);
      gel(y,1) = gpowgs(gel(x,1),m);
      gel(y,2) = gpowgs(gel(x,2),m);
      if (n < 0) y = ginv(y);
      return gerepileupto(av,y);
    }
    case t_POL:
      if (ismonome(x)) return pow_monome(x, n);
    default: {
      pari_sp av = avma;
      y = leftright_pow_u(x, (ulong)labs(n), NULL, &_sqr, &_mul);
      if (n < 0) y = ginv(y);
      return gerepileupto(av,y);
    }
  }
}

/* n a t_INT */
GEN
powgi(GEN x, GEN n)
{
  GEN y;

  if (!is_bigint(n)) return gpowgs(x, itos(n));
  /* probable overflow for non-modular types (typical exception: (X^0)^N) */
  switch(typ(x))
  {
    case t_INTMOD:
      y = cgetg(3,t_INTMOD); gel(y,1) = icopy(gel(x,1));
      gel(y,2) = Fp_pow(gel(x,2), n, gel(x,1));
      return y;
    case t_PADIC: return powp(x, n);

    case t_INT:
      if (is_pm1(x)) return (signe(x) < 0 && mpodd(n))? gen_m1: gen_1;
      if (signe(x)) pari_err(errlg);
      if (signe(n) < 0) pari_err(gdiver);
      return gen_0;
    case t_FRAC:
      if (signe(x[1])) pari_err(errlg);
      if (signe(n) < 0) pari_err(gdiver);
      return gen_0;

    case t_QFR: return qfr_pow(x,n);
    default: {
      pari_sp av = avma;
      y = leftright_pow(x, n, NULL, &_sqr, &_mul);
      if (signe(n) < 0) y = ginv(y);
      return gerepileupto(av,y);
    }
  }
}

/* we suppose n != 0, valp(x) = 0 and leading-term(x) != 0. Not stack clean */
static GEN
ser_pow(GEN x, GEN n, long prec)
{
  pari_sp av, tetpil;
  GEN y, p1, p2, lead;

  if (gvar(n) <= varn(x)) return gexp(gmul(n, glog(x,prec)), prec);
  lead = gel(x,2);
  if (gcmp1(lead))
  {
    GEN X, Y, alp = gadd(n,gen_1); /* will be left on the stack */
    long lx, mi, i, j, d;

    lx = lg(x);
    y = cgetg(lx,t_SER);
    y[1] = evalsigne(1) | evalvalp(0) | evalvarn(varn(x));
    X = x+2;
    Y = y+2;

    d = mi = lx-3; while (mi>=1 && isexactzero(gel(X,mi))) mi--;
    gel(Y,0) = gen_1;
    for (i=1; i<=d; i++)
    {
      av = avma; p1 = gen_0;
      for (j=1; j<=min(i,mi); j++)
      {
        p2 = gsubgs(gmulgs(alp,j),i);
        p1 = gadd(p1, gmul(gmul(p2,gel(X,j)),gel(Y,i-j)));
      }
      tetpil = avma; gel(Y,i) = gerepile(av,tetpil, gdivgs(p1,i));
    }
    return y;
  }
  p1 = gdiv(x,lead);
  if (typ(p1) != t_SER) pari_err(typeer, "ser_pow");
  gel(p1,2) = gen_1; /* in case it's inexact */
  if (typ(n) == t_FRAC && !isinexact(lead) && ispower(lead, gel(n,2), &p2))
    p2 = powgi(p2, gel(n,1));
  else
    p2 = gpow(lead,n, prec);
  return gmul(p2, gpow(p1,  n, prec));
}

static long
val_from_i(GEN E)
{
  if (is_bigint(E)) pari_err(talker,"valuation overflow in sqrtn");
  return itos(E);
}

/* return x^q, assume typ(x) = t_SER, typ(q) = t_FRAC and q != 0 */
static GEN
ser_powfrac(GEN x, GEN q, long prec)
{
  long e = valp(x);
  GEN y, E = gmulsg(e, q);

  if (gcmp0(x)) return zeroser(varn(x), val_from_i(gfloor(E)));
  if (typ(E) != t_INT)
    pari_err(talker,"%Z should divide valuation (= %ld) in sqrtn",q[2], e);
  e = val_from_i(E);
  y = shallowcopy(x); setvalp(y, 0);
  y = ser_pow(y, q, prec);
  if (typ(y) == t_SER) /* generic case */
    y[1] = evalsigne(1) | evalvalp(e) | evalvarn(varn(x));
  else /* e.g coeffs are POLMODs */
    y = gmul(y, monomial(gen_1, e, varn(x)));
  return y;
}

GEN
gpow(GEN x, GEN n, long prec)
{
  long i, lx, tx, tn = typ(n);
  pari_sp av;
  GEN y;

  if (tn == t_INT) return powgi(x,n);
  tx = typ(x);
  if (is_matvec_t(tx))
  {
    lx = lg(x); y = cgetg(lx,tx);
    for (i=1; i<lx; i++) gel(y,i) = gpow(gel(x,i),n,prec);
    return y;
  }
  av = avma;
  if (tx == t_POL || tx == t_RFRAC) { x = toser_i(x); tx = t_SER; }
  if (tx == t_SER)
  {
    if (tn == t_FRAC) return gerepileupto(av, ser_powfrac(x, n, prec));
    if (valp(x))
      pari_err(talker,"gpow: need integer exponent if series valuation != 0");
    if (lg(x) == 2) return gcopy(x); /* O(1) */
    return gerepileupto(av, ser_pow(x, n, prec));
  }
  if (gcmp0(x))
  {
    if (!is_scalar_t(tn) || tn == t_PADIC || tn == t_INTMOD)
      pari_err(talker,"gpow: 0 to a forbidden power");
    n = real_i(n);
    if (gsigne(n) <= 0)
      pari_err(talker,"gpow: 0 to a non positive exponent");
    if (!precision(x)) return gcopy(x);

    x = ground(gmulsg(gexpo(x),n));
    if (is_bigint(x) || (ulong)x[2] >= HIGHEXPOBIT)
      pari_err(talker,"gpow: underflow or overflow");
    avma = av; return real_0_bit(itos(x));
  }
  if (tn == t_FRAC)
  {
    GEN z, d = gel(n,2), a = gel(n,1);
    if (tx == t_INTMOD)
    {
      if (!BSW_psp(gel(x,1))) pari_err(talker,"gpow: modulus %Z is not prime",x[1]);
      y = cgetg(3,t_INTMOD); gel(y,1) = icopy(gel(x,1));
      av = avma;
      z = Fp_sqrtn(gel(x,2), d, gel(x,1), NULL);
      if (!z) pari_err(talker,"gpow: nth-root does not exist");
      gel(y,2) = gerepileuptoint(av, Fp_pow(z, a, gel(x,1)));
      return y;
    }
    else if (tx == t_PADIC)
    {
      z = equaliu(d, 2)? padic_sqrt(x): padic_sqrtn(x, d, NULL);
      if (!z) pari_err(talker, "gpow: nth-root does not exist");
      return gerepileupto(av, powgi(z, a));
    }
  }
  i = (long) precision(n); if (i) prec=i;
  y = gmul(n, glog(x,prec));
  return gerepileupto(av, gexp(y,prec));
}

/********************************************************************/
/**                                                                **/
/**                        RACINE CARREE                           **/
/**                                                                **/
/********************************************************************/

GEN
sqrtr(GEN x) {
  long s = signe(x);
  GEN y;
  if (typ(x) != t_REAL) pari_err(typeer,"sqrtr");
  if (s == 0) return real_0_bit(expo(x) >> 1);
  if (s >= 0) return sqrtr_abs(x);
  y = cgetg(3,t_COMPLEX);
  gel(y,2) = sqrtr_abs(x);
  gel(y,1) = gen_0; return y;
}

/* assume x unit, precp(x) = pp > 3 */
static GEN
sqrt_2adic(GEN x, long pp)
{
  GEN z = mod16(x)==1? gen_1: utoipos(3);
  long zp;
  pari_sp av, lim;

  if (pp == 4) return z;
  zp = 3; /* number of correct bits in z (compared to sqrt(x)) */

  av = avma; lim = stack_lim(av,2);
  for(;;)
  {
    GEN mod;
    zp = (zp<<1) - 1;
    if (zp > pp) zp = pp;
    mod = int2n(zp);
    z = addii(z, resmod2n(mulii(x, Fp_inv(z,mod)), zp));
    z = shifti(z, -1); /* (z + x/z) / 2 */
    if (pp == zp) return z;

    if (zp < pp) zp--;

    if (low_stack(lim,stack_lim(av,2)))
    {
      if (DEBUGMEM > 1) pari_warn(warnmem,"padic_sqrt");
      z = gerepileuptoint(av,z);
    }
  }
}

/* x unit defined modulo modx = p^pp, p != 2, pp > 0 */
static GEN
sqrt_padic(GEN x, GEN modx, long pp, GEN p)
{
  GEN mod, z = Fp_sqrt(x, p);
  long zp = 1;
  pari_sp av, lim;

  if (!z) pari_err(sqrter5);
  if (pp <= zp) return z;

  av = avma; lim = stack_lim(av,2);
  mod = p;
  for(;;)
  { /* could use the hensel_lift_accel idea. Not really worth it */
    GEN inv2;
    zp <<= 1;
    if (zp < pp) mod = sqri(mod); else { zp = pp; mod = modx; }
    inv2 = shifti(addis(mod,1), -1); /* = (mod + 1)/2 = 1/2 */
    z = addii(z, remii(mulii(x, Fp_inv(z,mod)), mod));
    z = mulii(z, inv2);
    z = modii(z, mod); /* (z + x/z) / 2 */
    if (pp <= zp) return z;

    if (low_stack(lim,stack_lim(av,2)))
    {
      GEN *gptr[2]; gptr[0]=&z; gptr[1]=&mod;
      if (DEBUGMEM>1) pari_warn(warnmem,"padic_sqrt");
      gerepilemany(av,gptr,2);
    }
  }
}

GEN
padic_sqrt(GEN x)
{
  long pp, e = valp(x);
  GEN z,y,mod, p = gel(x,2);
  pari_sp av;

  if (gcmp0(x)) return zeropadic(p, (e+1) >> 1);
  if (e & 1) pari_err(talker,"odd exponent in p-adic sqrt");

  y = cgetg(5,t_PADIC);
  pp = precp(x);
  mod = gel(x,3);
  x   = gel(x,4); /* lift to t_INT */
  e >>= 1; av = avma;
  if (equaliu(p,2))
  {
    long r = mod8(x);
    if (pp <= 3)
    {
      switch(pp) {
        case 1: break;
        case 2: if ((r&3) == 1) break;
        case 3: if (r == 1) break;
        default: pari_err(sqrter5);
      }
      z = gen_1;
      pp = 1;
    }
    else
    {
      if (r != 1) pari_err(sqrter5);
      z = sqrt_2adic(x, pp);
      z = gerepileuptoint(av, z);
      pp--;
    }
    mod = int2n(pp);
  }
  else /* p != 2 */
  {
    z = sqrt_padic(x, mod, pp, p);
    z = gerepileuptoint(av, z);
    mod = icopy(mod);
  }
  y[1] = evalprecp(pp) | evalvalp(e);
  gel(y,2) = icopy(p);
  gel(y,3) = mod;
  gel(y,4) = z; return y;
}

GEN
gsqrt(GEN x, long prec)
{
  pari_sp av;
  GEN y, p1, p2;

  switch(typ(x))
  {
    case t_REAL: return sqrtr(x);

    case t_INTMOD:
      y = cgetg(3,t_INTMOD); gel(y,1) = icopy(gel(x,1));
      p1 = Fp_sqrt(gel(x,2),gel(y,1));
      if (!p1) pari_err(sqrter5);
      gel(y,2) = p1; return y;

    case t_COMPLEX:
      if (isexactzero(gel(x,2))) return gsqrt(gel(x,1), prec);
      y = cgetg(3,t_COMPLEX); av = avma;

      p1 = gsqr(gel(x,1));
      p2 = gsqr(gel(x,2)); p1 = gsqrt(gadd(p1,p2), prec);
      if (gcmp0(p1)) { gel(y,1) = gel(y,2) = sqrtr(p1); return y; }
      if (gsigne(gel(x,1)) < 0)
      {
        p1 = sqrtr( gmul2n(gsub(p1,gel(x,1)), -1) );
        if (gsigne(gel(x,2)) < 0) setsigne(p1, -signe(p1));
        gel(y,2) = gerepileuptoleaf(av, p1); av = avma;
        gel(y,1) = gerepileuptoleaf(av, gdiv(gel(x,2), gmul2n(p1,1)));
      } else {
        p1 = sqrtr( gmul2n(gadd(p1,gel(x,1)), -1) );
        gel(y,1) = gerepileuptoleaf(av, p1); av = avma;
        gel(y,2) = gerepileuptoleaf(av, gdiv(gel(x,2), gmul2n(p1,1)));
      }
      return y;

    case t_PADIC:
      return padic_sqrt(x);

    default:
      av = avma; if (!(y = toser_i(x))) break;
      return gerepileupto(av, ser_powfrac(y, ghalf, prec));
  }
  return transc(gsqrt,x,prec);
}
/********************************************************************/
/**                                                                **/
/**                          N-th ROOT                             **/
/**                                                                **/
/********************************************************************/
/* exp(2Ipi/n), assume n positive t_INT */
GEN
rootsof1complex(GEN n, long prec)
{
  pari_sp av = avma;
  if (is_pm1(n)) return real_1(prec);
  if (lgefint(n)==3 && n[2]==2) return stor(-1, prec);
  return gerepileupto(av, exp_Ir( divri(Pi2n(1, prec), n) ));
}

/*Only the O() of y is used*/
GEN
rootsof1padic(GEN n, GEN y)
{
  pari_sp av0 = avma, av;
  GEN z, r = cgetp(y);

  av = avma; (void)Fp_sqrtn(gen_1,n,gel(y,2),&z);
  if (z==gen_0) { avma = av0; return gen_0; }/*should not happen*/
  z = padicsqrtnlift(gen_1, n, z, gel(y,2), precp(y));
  affii(z, gel(r,4)); avma = av; return r;
}

static GEN exp_p(GEN x);
/*compute the p^e th root of x p-adic, assume x != 0 */
GEN
padic_sqrtn_ram(GEN x, long e)
{
  pari_sp ltop=avma;
  GEN a, p = gel(x,2), n = powiu(p,e);
  long v = valp(x);
  if (v)
  {
    long z;
    v = sdivsi_rem(v, n, &z);
    if (z) return NULL;
    x = gcopy(x); setvalp(x,0);
  }
  /*If p=2 -1 is an root of unity in U1,we need an extra check*/
  if (lgefint(p)==3 && p[2]==2 && mod8(gel(x,4))!=signe(gel(x,4)))
    return NULL;
  a = exp_p(gdiv(palog(x), n));
  if (!a) return NULL;
  /*Here n=p^e and a^n=z*x where z is a (p-1)th-root of unity. Note that
      z^p=z; hence for b = a/z, then b^n=x. We say b=x/a^(n-1)*/
  a = gdiv(x, powgi(a,addis(n,-1))); if (v) setvalp(a,v);
  return gerepileupto(ltop,a);
}

/*compute the nth root of x p-adic p prime with n*/
GEN
padic_sqrtn_unram(GEN x, GEN n, GEN *zetan)
{
  pari_sp av;
  GEN Z, a, r, p = gel(x,2);
  long v = valp(x);
  if (v)
  {
    long z;
    v = sdivsi_rem(v,n,&z);
    if (z) return NULL;
  }
  r = cgetp(x); setvalp(r,v);
  Z = NULL; /* -Wall */
  if (zetan) Z = cgetp(x);
  av = avma; a = Fp_sqrtn(gel(x,4), n, p, zetan);
  if (!a) return NULL;
  affii(padicsqrtnlift(gel(x,4), n, a, p, precp(x)), gel(r,4));
  if (zetan)
  {
    affii(padicsqrtnlift(gen_1, n, *zetan, p, precp(x)), gel(Z,4));
    *zetan = Z;
  }
  avma = av; return r;
}

GEN
padic_sqrtn(GEN x, GEN n, GEN *zetan)
{
  pari_sp av = avma, tetpil;
  GEN q, p = gel(x,2);
  long e;
  if (!signe(x[4]))
  {
    long m = itos(n);
    if (zetan) *zetan = gen_1;
    return zeropadic(p, (valp(x)+m-1)/m);
  }
  /* treat the ramified part using logarithms */
  e = Z_pvalrem(n, p, &q);
  if (e) { x = padic_sqrtn_ram(x,e); if (!x) return NULL; }
  if (is_pm1(q))
  { /* finished */
    if (signe(q) < 0) x = ginv(x);
    x = gerepileupto(av, x);
    if (zetan)
      *zetan = (e && lgefint(p)==3 && p[2]==2)? gen_m1 /*-1 in Q_2*/
                                              : gen_1;
    return x;
  }
  tetpil = avma;
  /* use hensel lift for unramified case */
  x = padic_sqrtn_unram(x, q, zetan);
  if (!x) return NULL;
  if (zetan)
  {
    GEN *gptr[2];
    if (e && lgefint(p)==3 && p[2]==2)/*-1 in Q_2*/
    {
      tetpil = avma; x = gcopy(x); *zetan = gneg(*zetan);
    }
    gptr[0] = &x; gptr[1] = zetan;
    gerepilemanysp(av,tetpil,gptr,2);
    return x;
  }
  return gerepile(av,tetpil,x);
}

/* x^(1/n) */
GEN
sqrtnr(GEN x, long n)
{
  if (typ(x) != t_REAL) pari_err(typeer,"sqrtnr");
  return mpexp(divrs(mplog(x), n));
}

GEN
gsqrtn(GEN x, GEN n, GEN *zetan, long prec)
{
  long i, lx, tx;
  pari_sp av;
  GEN y, z;
  if (typ(n)!=t_INT) pari_err(talker,"second arg must be integer in gsqrtn");
  if (!signe(n)) pari_err(talker,"1/0 exponent in gsqrtn");
  if (is_pm1(n))
  {
    if (zetan) *zetan = gen_1;
    return (signe(n) > 0)? gcopy(x): ginv(x);
  }
  if (zetan) *zetan = gen_0;
  tx = typ(x);
  if (is_matvec_t(tx))
  {
    lx = lg(x); y = cgetg(lx,tx);
    for (i=1; i<lx; i++) gel(y,i) = gsqrtn(gel(x,i),n,NULL,prec);
    return y;
  }
  av = avma;
  switch(tx)
  {
  case t_INTMOD:
    z = gen_0;
    y = cgetg(3,t_INTMOD);  gel(y,1) = icopy(gel(x,1));
    if (zetan) { z = cgetg(3,t_INTMOD); gel(z,1) = gel(y,1); }
    gel(y,2) = Fp_sqrtn(gel(x,2),n,gel(x,1),zetan);
    if (!y[2]) {
      if (zetan) return gen_0;
      pari_err(talker,"nth-root does not exist in gsqrtn");
    }
    if (zetan) { gel(z,2) = *zetan; *zetan = z; }
    return y;

  case t_PADIC:
    y = padic_sqrtn(x,n,zetan);
    if (!y) {
      if (zetan) return gen_0;
      pari_err(talker,"nth-root does not exist in gsqrtn");
    }
    return y;

  case t_INT: case t_FRAC: case t_REAL: case t_COMPLEX:
    i = precision(x); if (i) prec = i;
    if (tx==t_INT && is_pm1(x) && signe(x) > 0)
     /*speed-up since there is no way to call rootsof1complex from gp*/
      y = real_1(prec);
    else if (gcmp0(x))
    {
      if (signe(n) < 0) pari_err(gdiver);
      if (isinexactreal(x))
      {
        long e = gexpo(x), junk;
        y = real_0_bit(e < 2? 0: sdivsi_rem(e, n, &junk));
      }
      else
        y = real_0(prec);
    }
    else
      y = gerepileupto(av, gexp(gdiv(glog(x,prec), n), prec));
    if (zetan) *zetan = rootsof1complex(n,prec);
    return y;

  case t_QUAD:
    return gsqrtn(quadtoc(x, prec), n, zetan, prec);

  default:
    av = avma; if (!(y = toser_i(x))) break;
    return gerepileupto(av, ser_powfrac(y, ginv(n), prec));
  }
  pari_err(typeer,"gsqrtn");
  return NULL;/* not reached */
}

/********************************************************************/
/**                                                                **/
/**                             EXP(X) - 1                         **/
/**                                                                **/
/********************************************************************/
/* exp(|x|) - 1, assume x != 0 */
GEN
exp1r_abs(GEN x)
{
  long l = lg(x), l2 = l+1, ex = expo(x), l1, i, n, m, s;
  GEN y = cgetr(l), p1, p2, p3, X, unr;
  pari_sp av2, av = avma;
  double a, b, beta, gama = 2.0 /* optimized for SUN3 */;
                                /* KB: 3.0 is better for UltraSparc */
  beta = 5. + bit_accuracy_mul(l, LOG2);
  a = sqrt(beta/(gama*LOG2));
  b = (BITS_IN_LONG-1-ex)
      + log2(a * (gama/2.718281828459045) / (double)(ulong)x[2]);
  if (a >= b)
  {
    n = (long)(1+a*gama);
    m = (long)(1+a-b);
    l2 += m>>TWOPOTBITS_IN_LONG;
  } else { /* rare ! */
    b = -1 - log((double)(ulong)x[2]) + (BITS_IN_LONG-1-ex)*LOG2; /*-1-log(x)*/
    n = (long)(1.1 + beta/b);
    m = 0;
  }
  unr=real_1(l2);
  p2 =real_1(l2); setlg(p2,3);
  X = cgetr(l2); affrr(x, X); setsigne(X, 1);
  if (m) setexpo(X, ex-m);

  s = 0; l1 = 3; av2 = avma;
  for (i=n; i>=2; i--)
  { /* compute X^(n-1)/n! + ... + X/2 + 1 */
    setlg(X,l1); p3 = divrs(X,i);
    s -= expo(p3); p1 = mulrr(p3,p2); setlg(p1,l1);
    l1 += s>>TWOPOTBITS_IN_LONG; if (l1>l2) l1=l2;
    s &= (BITS_IN_LONG-1);
    setlg(unr,l1); p1 = addrr_sign(unr,1, p1,1);
    setlg(p2,l1); affrr(p1,p2); avma = av2; /* p2 <- 1 + (X/i)*p2 */
  }
  setlg(X,l2); p2 = mulrr(X,p2);

  for (i=1; i<=m; i++)
  {
    setlg(p2,l2);
    p2 = mulrr(p2, addsr(2,p2));
  }
  affr_fixlg(p2,y); avma = av; return y;
}

GEN
mpexp1(GEN x)
{
  long sx = signe(x);
  GEN y, z;
  pari_sp av;
  if (!sx) return real_0_bit(expo(x));
  if (sx > 0) return exp1r_abs(x);
  /* compute exp(x) * (1 - exp(-x)) */
  av = avma; y = exp1r_abs(x);
  z = addsr(1, y); setsigne(z, -1);
  return gerepileupto(av, divrr(y, z));
}

/********************************************************************/
/**                                                                **/
/**                             EXP(X)                             **/
/**                                                                **/
/********************************************************************/

static GEN
mpexp_basecase(GEN x)
{
  pari_sp av = avma;
  GEN y = addsr(1, exp1r_abs(x));
  if (signe(x) < 0) y = ginv(y);
  return gerepileupto(av,y);
}

GEN
mpexp(GEN x)
{
  const long s = 6; /*Initial steps using basecase*/
  long i, n, mask, p, l, sx = signe(x), sh=0;
  GEN a, z;

  if (!sx) {
    long e = expo(x);
    return e >= 0? real_0_bit(e): real_1(nbits2prec(-e));
  }

  l = lg(x);
  if (l <= max(EXPNEWTON_LIMIT, 1<<s)) return mpexp_basecase(x);
  z = cgetr(l); /* room for result */
  if (expo(x) >= 0)
  { /* x>=1 : we do x %= log(2) to keep x small*/
    sh = (long) (rtodbl(x)/LOG2);
    x = subrr(rtor(x,l+1), mulsr(sh, mplog2(l+1)));
    if (!signe(x)) { avma = (pari_sp)(z+l); return real2n(sh, l); }
  }
  n = hensel_lift_accel(l-1,&mask);
  for(i=0, p=1; i<s; i++) { p <<= 1; if (mask&(1<<i)) p--; }
  a = mpexp_basecase(rtor(x, p+2));
  x = addrs(x,1);
  if (lg(x) < l+1) x = rtor(x, l+1);
  for(i=s; i<n; i++)
  {
    p <<= 1; if (mask&(1<<i)) p--;
    setlg(x, p+2); a = rtor(a, p+2);
    a = mulrr(a, subrr(x, logr_abs(a))); /* a := a (x - log(a)) */
  }
  affrr(a,z); 
  if (sh) setexpo(z, expo(z) + sh);
  avma = (pari_sp)z; return z;
}

static long
exp_p_prec(GEN x)
{
  long k, e = valp(x), n = e + precp(x);
  GEN p = gel(x,2);
  int is2 = equaliu(p,2);
  if (e < 1 || (e == 1 && is2)) return -1;
  if (is2)
  {
    n--; e--; k = n/e;
    if (n%e == 0) k--;
  }
  else
  {
    GEN r, t = subis(p, 1);
    k = itos(dvmdii(subis(mulis(t,n), 1), subis(mulis(t,e), 1), &r));
    if (!signe(r)) k--;
  }
  return k;
}

static GEN
exp_p(GEN x)
{
  long k;
  pari_sp av;
  GEN y;

  if (gcmp0(x)) return gaddgs(x,1);
  k = exp_p_prec(x);
  if (k < 0) return NULL;
  av = avma;
  for (y=gen_1; k; k--) y = gaddsg(1, gdivgs(gmul(y,x), k));
  return gerepileupto(av, y);
}
static GEN
cos_p(GEN x)
{
  long k;
  pari_sp av;
  GEN x2, y;

  if (gcmp0(x)) return gaddgs(x,1);
  k = exp_p_prec(x);
  if (k < 0) return NULL;
  av = avma; x2 = gsqr(x);
  if (k & 1) k--;
  for (y=gen_1; k; k-=2) 
  {
    GEN t = gdiv(gmul(y,x2), mulss(k, k-1));
    y = gsubsg(1, t);
  }
  return gerepileupto(av, y);
}
static GEN
sin_p(GEN x)
{
  long k;
  pari_sp av;
  GEN x2, y;

  if (gcmp0(x)) return gaddgs(x,1);
  k = exp_p_prec(x);
  if (k < 0) return NULL;
  av = avma; x2 = gsqr(x);
  if (k & 1) k--;
  for (y=gen_1; k; k-=2) 
  {
    GEN t = gdiv(gmul(y,x2), mulss(k, k+1));
    y = gsubsg(1, t);
  }
  return gerepileupto(av, gmul(y, x));
}

static GEN
cxexp(GEN x, long prec)
{
  GEN r,p1,p2, y = cgetg(3,t_COMPLEX);
  pari_sp av = avma, tetpil;
  r = gexp(gel(x,1),prec);
  if (gcmp0(r)) { gel(y,1) = r; gel(y,2) = r; return y; }
  gsincos(gel(x,2),&p2,&p1,prec);
  tetpil = avma;
  gel(y,1) = gmul(r,p1);
  gel(y,2) = gmul(r,p2);
  gerepilecoeffssp(av,tetpil,y+1,2);
  return y;
}

static GEN
serexp(GEN x, long prec)
{
  pari_sp av;
  long i,j,lx,ly,ex,mi;
  GEN p1,y,xd,yd;

  ex = valp(x);
  if (ex < 0) pari_err(negexper,"gexp");
  if (gcmp0(x)) return gaddsg(1,x);
  lx = lg(x);
  if (ex)
  {
    ly = lx+ex; y = cgetg(ly,t_SER);
    mi = lx-1; while (mi>=3 && isexactzero(gel(x,mi))) mi--;
    mi += ex-2;
    y[1] = evalsigne(1) | evalvalp(0) | evalvarn(varn(x));
    /* zd[i] = coefficient of X^i in z */
    xd = x+2-ex; yd = y+2; ly -= 2;
    gel(yd,0) = gen_1;
    for (i=1; i<ex; i++) gel(yd,i) = gen_0;
    for (   ; i<ly; i++)
    {
      av = avma; p1 = gen_0;
      for (j=ex; j<=min(i,mi); j++)
        p1 = gadd(p1, gmulgs(gmul(gel(xd,j),gel(yd,i-j)),j));
      gel(yd,i) = gerepileupto(av, gdivgs(p1,i));
    }
    return y;
  }
  av = avma; y = cgetg(lx, t_SER);
  y[1] = x[1]; gel(y,2) = gen_0;
  for (i=3; i <lx; i++) y[i] = x[i];
  p1 = gexp(gel(x,2),prec);
  y = gmul(p1, serexp(normalize(y),prec));
  return gerepileupto(av, y);
}

GEN
gexp(GEN x, long prec)
{
  switch(typ(x))
  {
    case t_REAL: return mpexp(x);
    case t_COMPLEX: return cxexp(x,prec);
    case t_PADIC: x = exp_p(x);
      if (!x) pari_err(talker,"p-adic argument out of range in gexp");
      return x;
    case t_INTMOD: pari_err(typeer,"gexp");
    default:
    {
      pari_sp av = avma;
      GEN y;
      if (!(y = toser_i(x))) break;
      return gerepileupto(av, serexp(y,prec));
    }
  }
  return transc(gexp,x,prec);
}

/********************************************************************/
/**                                                                **/
/**                           AGM(X, Y)                            **/
/**                                                                **/
/********************************************************************/
static int
agmr_gap(GEN a, GEN b, long L)
{
  GEN d = subrr(b, a);
  return (signe(d) && expo(d) - expo(b) >= L);
}
/* assume x > 0 */
static GEN
agm1r_abs(GEN x)
{
  long l = lg(x), L = 5-bit_accuracy(l);
  GEN a1, b1, y = cgetr(l);
  pari_sp av = avma;

  a1 = addrr(real_1(l), x); setexpo(a1, expo(a1)-1);
  b1 = sqrtr_abs(x);
  while (agmr_gap(a1,b1,L))
  {
    GEN a = a1;
    a1 = addrr(a,b1); setexpo(a1, expo(a1)-1);
    b1 = sqrtr_abs(mulrr(a,b1));
  }
  affr_fixlg(a1,y); avma = av; return y;
}

static int
agmcx_gap(GEN a, GEN b, long L)
{
  GEN d = gsub(b, a);
  return (!gcmp0(d) && gexpo(d) - gexpo(b) >= L);
}
static GEN
agm1cx(GEN x, long prec)
{
  GEN a1, b1;
  pari_sp av = avma, av2;
  long L, l = precision(x); if (!l) l = prec;

  L = 5-bit_accuracy(l);
  a1 = gtofp(gmul2n(gadd(real_1(l), x), -1), l); /* avoid loss of accuracy */
  av2 = avma;
  b1 = gsqrt(x, prec);
  while (agmcx_gap(a1,b1,L))
  {
    GEN a = a1;
    a1 = gmul2n(gadd(a,b1),-1);
    av2 = avma;
    b1 = gsqrt(gmul(a,b1), prec);
  }
  avma = av2; return gerepileupto(av,a1);
}

/* agm(1,x) */
static GEN
agm1(GEN x, long prec)
{
  GEN p1, a, a1, b1, y;
  long l, l2, ep;
  pari_sp av;

  if (gcmp0(x)) return gcopy(x);
  switch(typ(x))
  {
    case t_INT:
      if (!is_pm1(x)) break;
      return (signe(x) > 0)? real_1(prec): real_0(prec);

    case t_REAL: return signe(x) > 0? agm1r_abs(x): agm1cx(x, prec);

    case t_COMPLEX:
      if (gcmp0(gel(x,2)) && gsigne(gel(x,1)) > 0)
        return agm1(gel(x,1), prec);
      return agm1cx(x, prec);

    case t_PADIC:
      av = avma;
      a1 = x; b1 = gen_1; l = precp(x);
      do
      {
	a = a1;
	a1 = gmul2n(gadd(a,b1),-1);
        b1 = padic_sqrt(gmul(a,b1));
	p1 = gsub(b1,a1); ep = valp(p1)-valp(b1);
	if (ep<=0) { b1 = gneg_i(b1); p1 = gsub(b1,a1); ep=valp(p1)-valp(b1); }
      }
      while (ep<l && !gcmp0(p1));
      return gerepilecopy(av,a1);

    default:
      av = avma; if (!(y = toser_i(x))) break;
      a1 = y; b1 = gen_1; l = lg(y)-2;
      l2 = 5-bit_accuracy(prec);
      do
      {
	a = a1;
	a1 = gmul2n(gadd(a,b1),-1);
        b1 = ser_powfrac(gmul(a,b1), ghalf, prec);
	p1 = gsub(b1,a1); ep = valp(p1)-valp(b1);
      }
      while (ep<l && !gcmp0(p1)
                  && (!isinexactreal(p1) || gexpo(p1) - gexpo(b1) >= l2));
      return gerepilecopy(av,a1);
  }
  return transc(agm1,x,prec);
}

GEN
agm(GEN x, GEN y, long prec)
{
  long ty = typ(y);
  pari_sp av;

  if (is_matvec_t(ty))
  {
    ty = typ(x);
    if (is_matvec_t(ty)) pari_err(talker,"agm of two vector/matrices");
    swap(x, y);
  }
  if (gcmp0(y)) return gcopy(y);
  av = avma;
  return gerepileupto(av, gmul(y, agm1(gdiv(x,y), prec)));
}

/********************************************************************/
/**                                                                **/
/**                             LOG(X)                             **/
/**                                                                **/
/********************************************************************/
/* cf logagmr_abs(). Compute Pi/2agm(1, 4/2^n) ~ log(2^n) = n log(2) */
GEN
constlog2(long prec)
{
  pari_sp av;
  long l, n;
  GEN y, tmplog2;

  if (glog2 && lg(glog2) >= prec) return glog2;

  tmplog2 = newbloc(prec);
  *tmplog2 = evaltyp(t_REAL) | evallg(prec);
  av = avma;
  l = prec+1;
  n = bit_accuracy(l) >> 1;
  y = divrr(Pi2n(-1, l), agm1r_abs( real2n(2 - n, l) ));
  affrr(divrs(y,n), tmplog2);
  if (glog2) gunclone(glog2);
  glog2 = tmplog2; avma = av; return glog2;
}

GEN
mplog2(long prec)
{
  GEN x = cgetr(prec);
  affrr(constlog2(prec), x); return x;
}

/*return log(|x|), assuming x != 0 */
GEN
logr_abs(GEN X)
{
  pari_sp ltop, av;
  long EX, l1, l2, m, n, k, e, s, l = lg(X);
  double a, b;
  GEN z, x, y, y2, S, unr;
  ulong u, v;

  if (l > LOGAGM_LIMIT) return logagmr_abs(X);
  EX = expo(X);
  if (absrnz_egal2n(X)) return EX? mulsr(EX, mplog2(l)): real_0(l);

  av = avma; z = cgetr(l); ltop = avma;
  l2 = l+1; x = cgetr(l2); affrr(X,x);
  x[1] = evalsigne(1) | evalexpo(0);
  /* X = x 2^EX, 1 < x < 2 */
  av = avma; l -= 2;
  k = 2;
  u = ((ulong)x[k]) & (~HIGHBIT); /* x[2] - HIGHBIT, assuming HIGHBIT set */
  v = BITS_IN_LONG-1;
  while (!u) { v += BITS_IN_LONG; u = (ulong)x[++k]; } /* terminates: x>1 */
  a = (double)v - log2((double)u); /* ~ -log2(x - 1) */
  b = sqrt((BITS_IN_HALFULONG/3.0) * l);
  if (a <= b)
  {
    n = 1 + (long)(3*b);
    m = 1 + (long)(b-a);
    if ((ulong)m >= BITS_IN_LONG) { GEN t;
      l2 += m>>TWOPOTBITS_IN_LONG;
      t = cgetr(l2); affrr(x,t); x = t;
    }
    for (k=1; k<=m; k++) x = sqrtr_abs(x);
  }
  else
  {
    n = 1 + (long)(BITS_IN_HALFULONG*l / a);
    m = 0;
  }
  y = divrr(subrex01(x), addrex01(x)); /* = (x-1) / (x+1) ~ 0 */
  y2 = gsqr(y);
  /* log(x) = log(1+y) - log(1-y) = 2 \sum_{k odd} y^k / k */
  k = 2*n + 1;
  unr = real_1(l2); S = x; av = avma;
  setlg(S,  3);
  setlg(unr,3); affrr(divrs(unr,k), S); /* destroy x, not needed anymore */

  s = 0; e = expo(y2); l1 = 3;
  for (k -= 2; k > 0; k -= 2) /* k = 2n+1, ..., 1 */
  {
    GEN T; /* S = y^(2n+1-k)/(2n+1) + ... + 1 / k */
    setlg(y2, l1); T = mulrr(S,y2);
    setlg(unr,l1);
    s -= e; /* >= 0 */
    l1 += s>>TWOPOTBITS_IN_LONG; if (l1>l2) l1=l2;
    s &= (BITS_IN_LONG-1);
    setlg(S, l1);
    affrr(addrr(divrs(unr, k), T), S); avma = av;
  }
  setlg(S, l2); y = mulrr(y,S); /* = log(X)/2 */
  setexpo(y, expo(y)+m+1);
  if (EX) y = addrr(y, mulsr(EX, mplog2(l2)));
  affr_fixlg(y, z); avma = ltop; return z;
}

GEN
logagmr_abs(GEN q)
{
  long prec = lg(q), lim, e = expo(q);
  GEN z, y, Q;
  pari_sp av;

  if (absrnz_egal2n(q)) return e? mulsr(e, mplog2(prec)): real_0(prec);
  z = cgetr(prec); av = avma; prec++;
  lim = bit_accuracy(prec) >> 1;
  Q = cgetr(prec); affrr(q, Q);
  Q[1] = evalsigne(1) | evalexpo(lim);

  /* Pi / 2agm(1, 4/Q) ~ log(Q), q = Q * 2^(e-lim) */
  y = divrr(Pi2n(-1, prec), agm1r_abs( divsr(4, Q) ));
  y = addrr(y, mulsr(e - lim, mplog2(prec)));
  affr_fixlg(y, z); avma = av; return z;
}

/* assume Im(q) != 0 */
GEN
logagmcx(GEN q, long prec)
{
  GEN z, y, Q, a, b;
  long lim, e, ea, eb;
  pari_sp av;
  int neg = 0;
  
  z = cgetc(prec); av = avma; prec++;
  if (gsigne(gel(q,1)) < 0) { q = gneg(q); neg = 1; }
  lim = bit_accuracy(prec) >> 1;
  Q = gtofp(q, prec);
  a = gel(Q,1);
  b = gel(Q,2);
  if (gcmp0(a)) {
    affr_fixlg(logr_abs(b), gel(z,1));
    y = Pi2n(-1, prec);
    if (signe(b) < 0) setsigne(y, -1);
    affr_fixlg(y, gel(z,2)); avma = av; return z;
  }
  ea = expo(a);
  eb = expo(b);
  if (ea <= eb)
  {
    setexpo(Q[1], lim - eb + ea);
    setexpo(Q[2], lim);
    e = lim - eb;
  } else {
    setexpo(Q[1], lim);
    setexpo(Q[2], lim - ea + eb);
    e = lim - ea;
  }

  /* Pi / 2agm(1, 4/Q) ~ log(Q), q = Q * 2^e */
  y = gdiv(Pi2n(-1, prec), agm1cx( gdivsg(4, Q), prec ));
  a = gel(y,1);
  b = gel(y,2);
  a = addrr(a, mulsr(-e, mplog2(prec)));
  if (neg) b = gsigne(b) <= 0? gadd(b, mppi(prec))
                             : gsub(b, mppi(prec));
  affr_fixlg(a, gel(z,1));
  affr_fixlg(b, gel(z,2)); avma = av; return z;
}

GEN
mplog(GEN x)
{
  if (signe(x)<=0) pari_err(talker,"non positive argument in mplog");
  return logr_abs(x);
}

GEN
teich(GEN x)
{
  GEN p,q,y,z,aux,p1;
  long n, k;
  pari_sp av;

  if (typ(x)!=t_PADIC) pari_err(talker,"not a p-adic argument in teichmuller");
  if (!signe(x[4])) return gcopy(x);
  p = gel(x,2);
  q = gel(x,3);
  z = gel(x,4); y = cgetp(x); av = avma;
  if (equaliu(p,2))
    z = (mod4(z) & 2)? addsi(-1,q): gen_1;
  else
  {
    p1 = addsi(-1, p);
    z = remii(z, p);
    aux = diviiexact(addsi(-1,q),p1); n = precp(x);
    for (k=1; k<n; k<<=1)
      z = modii(mulii(z,addsi(1,mulii(aux,addsi(-1,Fp_pow(z,p1,q))))), q);
  }
  affii(z, gel(y,4)); avma = av; return y;
}

/* Let x = 1 mod p and y := (x-1)/(x+1) = 0 (p). Then
 * log(x) = log(1+y) - log(1-y) = 2 \sum_{k odd} y^k / k.
 * palogaux(x) returns the last sum (not multiplied by 2) */
static GEN
palogaux(GEN x)
{
  long k,e,pp;
  GEN y,s,y2, p = gel(x,2);

  if (equalii(gen_1, gel(x,4)))
  {
    long v = valp(x)+precp(x);
    if (equaliu(p,2)) v--;
    return zeropadic(p, v);
  }
  y = gdiv(gaddgs(x,-1), gaddgs(x,1));
  e = valp(y); /* > 0 */
  pp = e+precp(y);
  if (equaliu(p,2)) pp--;
  else
  {
    GEN p1;
    for (p1=utoipos(e); cmpui(pp,p1) > 0; pp++) p1 = mulii(p1, p);
    pp -= 2;
  }
  k = pp/e; if (!odd(k)) k--;
  y2 = gsqr(y); s = gdivgs(gen_1,k);
  while (k > 2)
  {
    k -= 2; s = gadd(gmul(y2,s), gdivgs(gen_1,k));
  }
  return gmul(s,y);
}

GEN
palog(GEN x)
{
  pari_sp av = avma;
  GEN y, p = gel(x,2);

  if (!signe(x[4])) pari_err(talker,"zero argument in palog");
  if (equaliu(p,2))
  {
    y = gsqr(x); setvalp(y,0);
    y = palogaux(y);
  }
  else
  { /* compute log(x^(p-1)) / (p-1) */
    GEN mod = gel(x,3), p1 = subis(p,1);
    y = cgetp(x);
    gel(y,4) = Fp_pow(gel(x,4), p1, mod);
    p1 = diviiexact(subis(mod,1), p1); /* 1/(1-p) */
    y = gmul(palogaux(y), mulis(p1, -2));
  }
  return gerepileupto(av,y);
}

GEN
log0(GEN x, long flag,long prec)
{
  switch(flag)
  {
    case 0:
    case 1: return glog(x,prec);
    default: pari_err(flagerr,"log");
  }
  return NULL; /* not reached */
}

GEN
glog(GEN x, long prec)
{
  pari_sp av, tetpil;
  GEN y, p1;

  switch(typ(x))
  {
    case t_REAL:
      if (signe(x) >= 0)
      {
        if (!signe(x)) pari_err(talker,"zero argument in mplog");
        return logr_abs(x);
      }
      y = cgetg(3,t_COMPLEX);
      gel(y,1) = logr_abs(x);
      gel(y,2) = mppi(lg(x)); return y;

    case t_COMPLEX:
      if (gcmp0(gel(x,2))) return glog(gel(x,1), prec);
      if (prec > LOGAGMCX_LIMIT) return logagmcx(x, prec);
      y = cgetg(3,t_COMPLEX);
      gel(y,2) = garg(x,prec);
      av = avma; p1 = glog(cxnorm(x),prec); tetpil = avma;
      gel(y,1) = gerepile(av,tetpil,gmul2n(p1,-1)); return y;

    case t_PADIC:
      return palog(x);

    case t_INTMOD: pari_err(typeer,"glog");
    default:
      av = avma; if (!(y = toser_i(x))) break;
      if (valp(y) || gcmp0(y)) pari_err(talker,"log is not meromorphic at 0");
      p1 = integ(gdiv(derivser(y), y), varn(y)); /* log(y)' = y'/y */
      if (!gcmp1(gel(y,2))) p1 = gadd(p1, glog(gel(y,2),prec));
      return gerepileupto(av, p1);
  }
  return transc(glog,x,prec);
}
/********************************************************************/
/**                                                                **/
/**                        SINE, COSINE                            **/
/**                                                                **/
/********************************************************************/

/* Reduce x0 mod Pi/2 to x in [-Pi/4, Pi/4]. Return cos(x)-1 */
static GEN
mpsc1(GEN x, long *ptmod8)
{
  long e = expo(x), l = lg(x), l1, l2, i, n, m, s;
  pari_sp av;
  double beta, a, b, d;
  GEN y, unr, p2, p1, x2;

  n = 0;
  if (e >= 0)
  {
    GEN q, z, pitemp = mppi(DEFAULTPREC + (e >> TWOPOTBITS_IN_LONG));
    setexpo(pitemp,-1);
    z = addrr(x,pitemp); /* = x + Pi/4 */
    if (expo(z) >= bit_accuracy(min(l, lg(z))) + 3) pari_err(precer,"mpsc1");
    setexpo(pitemp, 0);
    q = floorr( divrr(z,pitemp) ); /* round ( x / (Pi/2) ) */
    if (signe(q))
    {
      x = subrr(x, mulir(q, Pi2n(-1, l+1))); /* x mod Pi/2  */
      e = expo(x);
      n = mod4(q); if (n && signe(q) < 0) n = 4 - n;
    }
  }
  s = signe(x); *ptmod8 = (s < 0)? 4 + n: n;
  if (!s) return real_0_bit(-bit_accuracy(l)<<1);

  l = lg(x); l2 = l+1; y = cgetr(l);
  beta = 5. + bit_accuracy_mul(l2,LOG2);
  a = sqrt(beta / LOG2);
  d = a + 1/LOG2 - log2(a / (double)(ulong)x[2]) - (BITS_IN_LONG-1-e);
  if (d >= 0)
  {
    n = (long)((1+a) / 2.0);
    m = (long)(1+d);
    l2 += m>>TWOPOTBITS_IN_LONG;
  } else { /* rare ! */
    b = -1 - log((double)(ulong)x[2]) + (BITS_IN_LONG-1-e)*LOG2; /*-1-log(x)*/
    n = (long)(1 + beta/(2.0*b));
    m = 0;
  }
  unr= real_1(l2);
  p2 = real_1(l2);
  x2 = cgetr(l2); av = avma;
  affrr(gsqr(x), x2);
  if (m) setexpo(x2, expo(x2) - (m<<1));

  setlg(x2, 3); p1 = divrs(x2, 2*n+1);
  s = -expo(p1);
  l1 = 3 + (s>>TWOPOTBITS_IN_LONG);
  setlg(p2,l1);
  s = 0;
  for (i=n; i>=2; i--)
  {
    setlg(x2,l1); p1 = divrsns(x2, 2*i-1);
    s -= expo(p1); p1 = mulrr(p1,p2);
    l1 += s>>TWOPOTBITS_IN_LONG; if (l1>l2) l1=l2;
    s &= (BITS_IN_LONG-1);
    setlg(unr,l1); p1 = addrr_sign(unr,1, p1,-signe(p1));
    setlg(p2,l1); affrr(p1,p2); avma = av;
  }
  p2[1] = evalsigne(-signe(p2)) | evalexpo(expo(p2)-1); /* p2 := -p2/2 */
  setlg(p2,l2);
  setlg(x2,l2); p2 = mulrr(x2,p2);
  /* Now p2 = sum {1<= i <=n} (-1)^i x^(2i) / (2i)! ~ cos(x) - 1 */
  for (i=1; i<=m; i++)
  { /* p2 = cos(x)-1 --> cos(2x)-1 */
    setlg(p2,l2);
    p2 = mulrr(p2, addsr(2,p2));
    setexpo(p2, expo(p2)+1);
  }
  affr_fixlg(p2,y); return y;
}

/* sqrt (|1 - (1+x)^2|) = sqrt(|x*(x+2)|). Sends cos(x)-1 to |sin(x)| */
static GEN
mpaut(GEN x)
{
  pari_sp av = avma;
  GEN t = mulrr(x, addsr(2,x)); /* != 0 */
  if (!signe(t)) return real_0_bit(expo(t) >> 1);
  return gerepileuptoleaf(av, sqrtr_abs(t));
}

/********************************************************************/
/**                            COSINE                              **/
/********************************************************************/

GEN
mpcos(GEN x)
{
  long mod8;
  pari_sp av;
  GEN y,p1;

  if (!signe(x)) return real_1(3 + ((-expo(x)) >> TWOPOTBITS_IN_LONG));

  av = avma; p1 = mpsc1(x,&mod8);
  switch(mod8)
  {
    case 0: case 4: y = addsr(1,p1); break;
    case 1: case 7: y = mpaut(p1); setsigne(y,-signe(y)); break;
    case 2: case 6: y = subsr(-1,p1); break;
    default:        y = mpaut(p1); break; /* case 3: case 5: */
  }
  return gerepileuptoleaf(av, y);
}

/* convert INT or FRAC to REAL, which is later reduced mod 2Pi : avoid
 * cancellation */
static GEN
tofp_safe(GEN x, long prec)
{
  return (typ(x) == t_INT || gexpo(x) > 0)? gadd(x, real_0(prec)) 
                                          : fractor(x, prec);
}

GEN
gcos(GEN x, long prec)
{
  pari_sp av;
  GEN r, u, v, y, u1, v1;
  long i;

  switch(typ(x))
  {
    case t_REAL:
      return mpcos(x);

    case t_COMPLEX:
      i = precision(x); if (!i) i = prec;
      y = cgetc(i); av = avma;
      r = gexp(gel(x,2),prec);
      v1 = gmul2n(addrr(ginv(r),r), -1); /* = cos(I*Im(x)) */
      u1 = subrr(v1, r); /* = - I*sin(I*Im(x)) */
      gsincos(gel(x,1),&u,&v,prec);
      affr_fixlg(gmul(v1,v), gel(y,1));
      affr_fixlg(gmul(u1,u), gel(y,2)); return y;

    case t_INT: case t_FRAC:
      y = cgetr(prec); av = avma;
      affr_fixlg(mpcos(tofp_safe(x,prec)), y); avma = av; return y;

    case t_INTMOD: pari_err(typeer,"gcos");
    
    case t_PADIC: x = cos_p(x);
      if (!x) pari_err(talker,"p-adic argument out of range in gcos");
      return x;

    default:
      av = avma; if (!(y = toser_i(x))) break;
      if (gcmp0(y)) return gaddsg(1,y);
      if (valp(y) < 0) pari_err(negexper,"gcos");
      gsincos(y,&u,&v,prec);
      return gerepilecopy(av,v);
  }
  return transc(gcos,x,prec);
}
/********************************************************************/
/**                             SINE                               **/
/********************************************************************/

GEN
mpsin(GEN x)
{
  long mod8;
  pari_sp av;
  GEN y,p1;

  if (!signe(x)) return real_0_bit(expo(x));

  av = avma; p1 = mpsc1(x,&mod8);
  switch(mod8)
  {
    case 0: case 6: y=mpaut(p1); break;
    case 1: case 5: y=addsr(1,p1); break;
    case 2: case 4: y=mpaut(p1); setsigne(y,-signe(y)); break;
    default:        y=subsr(-1,p1); break; /* case 3: case 7: */
  }
  return gerepileuptoleaf(av, y);
}

GEN
gsin(GEN x, long prec)
{
  pari_sp av;
  GEN r, u, v, y, v1, u1;
  long i;

  switch(typ(x))
  {
    case t_REAL:
      return mpsin(x);

    case t_COMPLEX:
      i = precision(x); if (!i) i = prec;
      y = cgetc(i); av = avma;
      r = gexp(gel(x,2),prec);
      v1 = gmul2n(addrr(ginv(r),r), -1); /* = cos(I*Im(x)) */
      u1 = subrr(r, v1); /* = I*sin(I*Im(x)) */
      gsincos(gel(x,1),&u,&v,prec);
      affr_fixlg(gmul(v1,u), gel(y,1));
      affr_fixlg(gmul(u1,v), gel(y,2)); return y;

    case t_INT: case t_FRAC:
      y = cgetr(prec); av = avma;
      affr_fixlg(mpsin(tofp_safe(x,prec)), y); avma = av; return y;

    case t_INTMOD: pari_err(typeer,"gsin");

    case t_PADIC: x = sin_p(x);
      if (!x) pari_err(talker,"p-adic argument out of range in gsin");
      return x;

    default:
      av = avma; if (!(y = toser_i(x))) break;
      if (gcmp0(y)) return gcopy(y);
      if (valp(y) < 0) pari_err(negexper,"gsin");
      gsincos(y,&u,&v,prec);
      return gerepilecopy(av,u);
  }
  return transc(gsin,x,prec);
}
/********************************************************************/
/**                       SINE, COSINE together                    **/
/********************************************************************/

void
mpsincos(GEN x, GEN *s, GEN *c)
{
  long mod8;
  pari_sp av, tetpil;
  GEN p1, *gptr[2];

  if (!signe(x))
  {
    long e = expo(x);
    *s = real_0_bit(e);
    *c = e >= 0? real_0_bit(e): real_1(nbits2prec(-e));
    return;
  }

  av=avma; p1=mpsc1(x,&mod8); tetpil=avma;
  switch(mod8)
  {
    case 0: *c=addsr( 1,p1); *s=mpaut(p1); break;
    case 1: *s=addsr( 1,p1); *c=mpaut(p1); setsigne(*c,-signe(*c)); break;
    case 2: *c=subsr(-1,p1); *s=mpaut(p1); setsigne(*s,-signe(*s)); break;
    case 3: *s=subsr(-1,p1); *c=mpaut(p1); break;
    case 4: *c=addsr( 1,p1); *s=mpaut(p1); setsigne(*s,-signe(*s)); break;
    case 5: *s=addsr( 1,p1); *c=mpaut(p1); break;
    case 6: *c=subsr(-1,p1); *s=mpaut(p1); break;
    case 7: *s=subsr(-1,p1); *c=mpaut(p1); setsigne(*c,-signe(*c)); break;
  }
  gptr[0]=s; gptr[1]=c;
  gerepilemanysp(av,tetpil,gptr,2);
}

/* return exp(ix), x a t_REAL */
GEN
exp_Ir(GEN x)
{
  pari_sp av = avma;
  GEN v = cgetg(3,t_COMPLEX);
  mpsincos(x, (GEN*)(v+2), (GEN*)(v+1));
  if (!signe(x)) return gerepilecopy(av, gel(v,1));
  return v;
}

void
gsincos(GEN x, GEN *s, GEN *c, long prec)
{
  long ii, i, j, ex, ex2, lx, ly, mi;
  pari_sp av, tetpil;
  GEN y, r, u, v, u1, v1, p1, p2, p3, p4, ps, pc;
  GEN *gptr[4];

  switch(typ(x))
  {
    case t_INT: case t_FRAC:
      *s = cgetr(prec);
      *c = cgetr(prec); av = avma;
      mpsincos(tofp_safe(x, prec), &ps, &pc);
      affr_fixlg(ps,*s);
      affr_fixlg(pc,*c); avma = av; return;

    case t_REAL:
      mpsincos(x,s,c); return;

    case t_COMPLEX:
      i = precision(x); if (!i) i = prec;
      ps = cgetc(i); *s = ps;
      pc = cgetc(i); *c = pc; av = avma;
      r = gexp(gel(x,2),prec);
      v1 = gmul2n(addrr(ginv(r),r), -1); /* = cos(I*Im(x)) */
      u1 = subrr(r, v1); /* = I*sin(I*Im(x)) */
      gsincos(gel(x,1), &u,&v, prec);
      affr_fixlg(mulrr(v1,u),       gel(ps,1));
      affr_fixlg(mulrr(u1,v),       gel(ps,2));
      affr_fixlg(mulrr(v1,v),       gel(pc,1));
      affr_fixlg(mulrr(mpneg(u1),u),gel(pc,2)); return;

    case t_QUAD:
      av = avma; gsincos(quadtoc(x, prec), s, c, prec);
      gerepileall(av, 2, s, c); return;

    default:
      av = avma; if (!(y = toser_i(x))) break;
      if (gcmp0(y)) { *c = gaddsg(1,y); *s = gcopy(y); return; }

      ex = valp(y); lx = lg(y); ex2 = 2*ex+2;
      if (ex < 0) pari_err(talker,"non zero exponent in gsincos");
      if (ex2 > lx)
      {
        *s = x == y? gcopy(y): gerepilecopy(av, y); av = avma;
        *c = gerepileupto(av, gsubsg(1, gdivgs(gsqr(y),2)));
	return;
      }
      if (!ex)
      {
        p1 = shallowcopy(y); gel(p1,2) = gen_0;
        gsincos(normalize(p1),&u,&v,prec);
        gsincos(gel(y,2),&u1,&v1,prec);
        p1 = gmul(v1,v);
        p2 = gmul(u1,u);
        p3 = gmul(v1,u);
        p4 = gmul(u1,v); tetpil = avma;
        *c = gsub(p1,p2);
        *s = gadd(p3,p4);
	gptr[0]=s; gptr[1]=c;
	gerepilemanysp(av,tetpil,gptr,2);
	return;
      }

      ly = lx+2*ex;
      mi = lx-1; while (mi>=3 && isexactzero(gel(y,mi))) mi--;
      mi += ex-2;
      pc = cgetg(ly,t_SER); *c = pc;
      ps = cgetg(lx,t_SER); *s = ps;
      pc[1] = evalsigne(1) | evalvalp(0) | evalvarn(varn(y));
      gel(pc,2) = gen_1; ps[1] = y[1];
      for (i=2; i<ex+2; i++) gel(ps,i) = gcopy(gel(y,i));
      for (i=3; i< ex2; i++) gel(pc,i) = gen_0;
      for (i=ex2; i<ly; i++)
      {
	ii = i-ex; av = avma; p1 = gen_0;
	for (j=ex; j<=min(ii-2,mi); j++)
	  p1 = gadd(p1, gmulgs(gmul(gel(y,j-ex+2),gel(ps,ii-j)),j));
	gel(pc,i) = gerepileupto(av, gdivgs(p1,2-i));
	if (ii < lx)
	{
	  av = avma; p1 = gen_0;
	  for (j=ex; j<=min(i-ex2,mi); j++)
	    p1 = gadd(p1,gmulgs(gmul(gel(y,j-ex+2),gel(pc,i-j)),j));
	  p1 = gdivgs(p1,i-2);
	  gel(ps,i-ex) = gerepileupto(av, gadd(p1,gel(y,i-ex)));
	}
      }
      return;
  }
  pari_err(typeer,"gsincos");
}

/********************************************************************/
/**                                                                **/
/**                     TANGENT and COTANGENT                      **/
/**                                                                **/
/********************************************************************/
static GEN
mptan(GEN x)
{
  pari_sp av = avma;
  GEN s, c;

  mpsincos(x,&s,&c);
  if (!signe(c)) pari_err(talker, "can't compute tan(Pi/2 + kPi)");
  return gerepileuptoleaf(av, divrr(s,c));
}

GEN
gtan(GEN x, long prec)
{
  pari_sp av;
  GEN y, s, c;

  switch(typ(x))
  {
    case t_REAL: return mptan(x);

    case t_COMPLEX:
      av = avma; gsincos(x,&s,&c,prec);
      return gerepileupto(av, gdiv(s,c));

    case t_INT: case t_FRAC:
      y = cgetr(prec); av = avma;
      affr_fixlg(mptan(tofp_safe(x,prec)), y); avma = av; return y;

    case t_PADIC:
      av = avma;
      return gerepileupto(av, gdiv(gsin(x,prec), gcos(x,prec)));

    case t_INTMOD: pari_err(typeer,"gtan");

    default:
      av = avma; if (!(y = toser_i(x))) break;
      if (gcmp0(y)) return gcopy(y);
      if (valp(y) < 0) pari_err(negexper,"gtan");
      gsincos(y,&s,&c,prec);
      return gerepileupto(av, gdiv(s,c));
  }
  return transc(gtan,x,prec);
}

static GEN
mpcotan(GEN x)
{
  pari_sp av=avma, tetpil;
  GEN s,c;

  mpsincos(x,&s,&c); tetpil=avma;
  return gerepile(av,tetpil,divrr(c,s));
}

GEN
gcotan(GEN x, long prec)
{
  pari_sp av;
  GEN y, s, c;

  switch(typ(x))
  {
    case t_REAL:
      return mpcotan(x);

    case t_COMPLEX:
      av = avma; gsincos(x,&s,&c,prec);
      return gerepileupto(av, gdiv(c,s));

    case t_INT: case t_FRAC:
      y = cgetr(prec); av = avma;
      affr_fixlg(mpcotan(tofp_safe(x,prec)), y); avma = av; return y;

    case t_PADIC: 
      av = avma;
      return gerepileupto(av, gdiv(gcos(x,prec), gsin(x,prec)));

    case t_INTMOD: pari_err(typeer,"gcotan");

    default:
      av = avma; if (!(y = toser_i(x))) break;
      if (gcmp0(y)) pari_err(talker,"0 argument in cotan");
      if (valp(y) < 0) pari_err(negexper,"cotan"); /* fall through */
      gsincos(y,&s,&c,prec);
      return gerepileupto(av, gdiv(c,s));
  }
  return transc(gcotan,x,prec);
}
