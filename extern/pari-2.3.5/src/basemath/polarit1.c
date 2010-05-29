/* $Id: polarit1.c 11851 2009-07-21 21:52:51Z bill $

Copyright (C) 2000-2004  The PARI group.

This file is part of the PARI/GP package.

PARI/GP is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation. It is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY WHATSOEVER.

Check the License for details. You should have received a copy of it, along
with the package; see the file 'COPYING'. If not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA. */

/***********************************************************************/
/**                                                                   **/
/**               ARITHMETIC OPERATIONS ON POLYNOMIALS                **/
/**                         (first part)                              **/
/**                                                                   **/
/***********************************************************************/
#include "pari.h"
#include "paripriv.h"

/*******************************************************************/
/*                                                                 */
/*                           DIVISIBILITY                          */
/*                 Return 1 if y  |  x,  0 otherwise               */
/*                                                                 */
/*******************************************************************/

int
gdvd(GEN x, GEN y)
{
  pari_sp av=avma;
  x=gmod(x,y); avma=av; return gcmp0(x);
}

int
poldvd(GEN x, GEN y, GEN *z)
{
  pari_sp av = avma;
  GEN p1 = poldivrem(x,y,ONLY_DIVIDES);
  if (p1) { *z = p1; return 1; }
  avma=av; return 0;
}

/*******************************************************************/
/*                                                                 */
/*                  POLYNOMIAL EUCLIDEAN DIVISION                  */
/*                                                                 */
/*******************************************************************/
GEN
poldivrem(GEN x, GEN y, GEN *pr)
{
  long ty = typ(y), tx, vx = gvar(x), vy = gvar(y);
  GEN p1;

  if (is_scalar_t(ty) || varncmp(vx, vy) < 0)
  {
    if (pr == ONLY_REM) {
      if (gcmp0(y)) pari_err(gdiver);
      return gen_0;
    }
    if (pr && pr != ONLY_DIVIDES) *pr=gen_0;
    return gdiv(x,y);
  }
  if (ty != t_POL) pari_err(typeer,"euclidean division (poldivrem)");
  tx = typ(x);
  if (is_scalar_t(tx) || varncmp(vx, vy) > 0)
  {
    if (!signe(y)) pari_err(gdiver);
    if (!degpol(y)) /* constant */
    {
      if (pr == ONLY_REM) return zeropol(vy);
      p1 = gdiv(x, gel(y,2));
      if (pr == ONLY_DIVIDES) return p1;
      if (pr) *pr = zeropol(vy);
      return p1;
    }
    if (pr == ONLY_REM) return gcopy(x);
    if (pr == ONLY_DIVIDES) return gcmp0(x)? gen_0: NULL;
    if (pr) *pr = gcopy(x);
    return gen_0;
  }
  if (tx != t_POL) pari_err(typeer,"euclidean division (poldivrem)");

  if (varncmp(vx, vy) < 0)
  {
    if (pr && pr != ONLY_DIVIDES)
    {
      p1 = zeropol(vx); if (pr == ONLY_REM) return p1;
      *pr = p1;
    }
    return gdiv(x,y);
  }
  return RgX_divrem(x, y, pr);
}

/*******************************************************************/
/*                                                                 */
/*           ROOTS MODULO a prime p (no multiplicities)            */
/*                                                                 */
/*******************************************************************/
static ulong
init_p(GEN pp)
{
  ulong p;
  if ((ulong)expi(pp) > BITS_IN_LONG - 3) p = 0;
  else
  {
    p = itou(pp);
    if (p < 2 || signe(pp) < 0) pari_err(talker,"not a prime in factmod");
  }
  return p;
}

static long
factmod_init(GEN *F, GEN p)
{
  long d;
  if (typ(*F)!=t_POL || typ(p)!=t_INT) pari_err(typeer,"factmod");
  *F = FpX_normalize(RgX_to_FpX(*F, p), p);
  d = degpol(*F); if (d < 0) pari_err(zeropoler,"factmod");
  return d;
}

static GEN
root_mod_2(GEN f)
{
  int z1, z0 = !signe(constant_term(f));
  long i,n;
  GEN y;

  for (i=2, n=1; i < lg(f); i++)
    if (signe(f[i])) n++;
  z1 = n & 1;
  y = cgetg(z0+z1+1, t_COL); i = 1;
  if (z0) gel(y,i++) = gen_0;
  if (z1) gel(y,i) = gen_1;
  return y;
}

#define i_mod4(x) (signe(x)? mod4((GEN)(x)): 0)
static GEN
root_mod_4(GEN f)
{
  long i, no, ne;
  GEN y, t;
  int z0, z1, z2, z3;

  t = constant_term(f);
  z0 = !signe(t);
  z2 = ((i_mod4(t) + 2*i_mod4(f[3])) & 3) == 0;

  for (ne=0,i=2; i<lg(f); i+=2)
  {
    t = gel(f,i);
    if (signe(t)) ne += *int_LSW(t);
  }
  for (no=0,i=3; i<lg(f); i+=2)
  {
    t = gel(f,i);
    if (signe(t)) ne += *int_LSW(t);
  }
  no &= 3; ne &= 3;
  z3 = (no == ne);
  z1 = (no == ((4-ne)&3));
  y=cgetg(1+z0+z1+z2+z3,t_COL); i = 1;
  if (z0) gel(y,i++) = gen_0;
  if (z1) gel(y,i++) = gen_1;
  if (z2) gel(y,i++) = gen_2;
  if (z3) gel(y,i) = utoipos(3);
  return y;
}
#undef i_mod4

/* p even, accept p = 4 for p-adic stuff */
INLINE GEN
root_mod_even(GEN f, ulong p)
{
  switch(p)
  {
    case 2: return root_mod_2(f);
    case 4: return root_mod_4(f);
  }
  pari_err(talker,"not a prime in rootmod");
  return NULL; /* not reached */
}

/* by checking f(0..p-1) */
static GEN
Flx_roots_naive(GEN f, ulong p)
{
  long d = degpol(f), n = 0;
  ulong s = 1UL, r;
  GEN q, y = cgetg(d + 1, t_VECSMALL);
  pari_sp av = avma;

  if (!f[2]) y[++n] = 0;
  do
  {
    q = Flx_div_by_X_x(f, s, p, &r); /* TODO: FFT-type multi-evaluation */
    if (r) avma = av;
    else
    {
      y[++n] = s;
      f = q; av = avma;
    }
    s++;
  }
  while (n < d-1 && p > s);
  if (n == d-1 && p != s) /* -f[2]/f[3] */
    y[++n] = Fl_mul(p - Fl_inv((ulong)f[3], p), (ulong)f[2], p);
  setlg(y, n+1); return y;
}

GEN
rootmod2(GEN f, GEN pp)
{
  pari_sp av = avma;
  ulong p;
  GEN y;

  if (!factmod_init(&f, pp)) { avma = av; return cgetg(1,t_COL); }
  p = init_p(pp); if (!p) pari_err(talker,"prime too big in rootmod2");
  if (p & 1)
    y = Flc_to_ZC(Flx_roots_naive(ZX_to_Flx(f,p), p));
  else
    y = root_mod_even(f,p);
  return gerepileupto(av, FpC_to_mod(y, pp));
}

/* assume x reduced mod p, monic. */
static int
FpX_quad_factortype(GEN x, GEN p)
{
  GEN b = gel(x,3), c = gel(x,2);
  GEN D;

  if (equaliu(p, 2)) {
    if (!signe(b)) return 0;
    return signe(c)? -1: 1;
  }
  D = subii(sqri(b), shifti(c,2));
  return kronecker(D,p);
}
/* assume x reduced mod p, monic. Return one root, or NULL if irreducible */
GEN
FpX_quad_root(GEN x, GEN p, int unknown)
{
  GEN s, u, D, b = gel(x,3), c = gel(x,2);

  if (equaliu(p, 2)) {
    if (!signe(b)) return c;
    return signe(c)? NULL: gen_1;
  }
  D = subii(sqri(b), shifti(c,2));
  D = remii(D,p);
  if (unknown && kronecker(D,p) == -1) return NULL;

  s = Fp_sqrt(D,p);
  /* p is not prime, go on and give e.g. maxord a chance to recover */
  if (!s) return NULL;
  u = addis(shifti(p,-1), 1); /* = 1/2 */
  return modii(mulii(u, subii(s,b)), p);
}
static GEN
otherroot(GEN x, GEN r, GEN p)
{
  GEN s = addii(gel(x,3), r);
  if (!signe(s)) return s;
  s = subii(p, s); if (signe(s) < 0) s = addii(s,p);
  return s;
}

/* by splitting, assume p > 2 prime, deg(f) > 0, and f monic */
static GEN
FpX_roots_i(GEN f, GEN p)
{
  long n, j, da, db;
  GEN y, pol, pol0, a, b, q = shifti(p,-1);

  n = ZX_valuation(f, &f)? 1: 0;
  y = cgetg(lg(f), t_COL);
  j = 1;
  if (n) {
    gel(y, j++) = gen_0;
    if (lg(f) <= 3) { setlg(y, 2); return y; }
    n = 1;
  }
  da = degpol(f);
  if (da == 1) { gel(y,j++) = subii(p, gel(f,2)); setlg(y,j); return y; }
  if (da == 2) {
    GEN s, r = FpX_quad_root(f, p, 1);
    if (r) {
      gel(y, j++) = r;
      s = otherroot(f,r, p);
      if (!equalii(r, s)) gel(y, j++) = s;
    }
    setlg(y, j); return sort(y);
  }

  /* take gcd(x^(p-1) - 1, f) by splitting (x^q-1) * (x^q+1) */
  b = FpXQ_pow(pol_x[varn(f)],q, f,p);
  if (lg(b) < 3) pari_err(talker,"not a prime in rootmod");
  b = ZX_Z_add(b, gen_m1); /* b = x^((p-1)/2) - 1 mod f */
  a = FpX_gcd(f,b, p);
  b = ZX_Z_add(b, gen_2); /* b = x^((p-1)/2) + 1 mod f */
  b = FpX_gcd(f,b, p);
  da = degpol(a);
  db = degpol(b); n += da + db; setlg(y, n+1);
  if (db) gel(y,j)    = FpX_normalize(b,p);
  if (da) gel(y,j+db) = FpX_normalize(a,p);
  pol = gadd(pol_x[varn(f)], gen_1); pol0 = constant_term(pol);
  while (j <= n)
  { /* cf FpX_split_Berlekamp */
    a = gel(y,j); da = degpol(a);
    if (da==1)
      gel(y,j++) = subii(p, gel(a,2));
    else if (da==2)
    {
      GEN r = FpX_quad_root(a, p, 0);
      gel(y, j++) = r;
      gel(y, j++) = otherroot(a,r, p);
    }
    else for (pol0[2]=1; ; pol0[2]++)
    {
      b = ZX_Z_add(FpXQ_pow(pol,q, a,p), gen_m1); /* pol^(p-1)/2 - 1 */
      b = FpX_gcd(a,b, p); db = degpol(b);
      if (db && db < da)
      {
        b = FpX_normalize(b, p);
        gel(y,j+db) = FpX_div(a,b, p);
        gel(y,j)    = b; break;
      }
      if (pol0[2] == 100 && !BSW_psp(p))
        pari_err(talker, "not a prime in polrootsmod");
    }
  }
  return sort(y);
}

static GEN
FpX_factmod_init(GEN f, GEN p) { return FpX_normalize(FpX_red(f,p), p); }
GEN
FpX_roots(GEN f, GEN p) {
  pari_sp av = avma;
  long q = modBIL(p);
  GEN F = FpX_factmod_init(f,p);
  switch(degpol(F))
  {
    case -1: pari_err(zeropoler,"factmod");
    case 0: avma = av; return cgetg(1, t_COL);
  }
  return gerepileupto(av, odd(q)? FpX_roots_i(F, p): root_mod_even(F,q));
}

GEN
rootmod(GEN f, GEN p)
{
  ulong q;
  pari_sp av = avma;
  GEN y;

  if (!factmod_init(&f, p)) { avma=av; return cgetg(1,t_COL); }
  q = init_p(p); if (!q) q = modBIL(p);
  if (q & 1)
    y = FpX_roots_i(f, p);
  else
    y = root_mod_even(f,q);
  return gerepileupto(av, FpC_to_mod(y, p));
}

GEN
rootmod0(GEN f, GEN p, long flag)
{
  switch(flag)
  {
    case 0: return rootmod(f,p);
    case 1: return rootmod2(f,p);
    default: pari_err(flagerr,"polrootsmod");
  }
  return NULL; /* not reached */
}

/*******************************************************************/
/*                                                                 */
/*                     FACTORISATION MODULO p                      */
/*                                                                 */
/*******************************************************************/
static GEN spec_FpXQ_pow(GEN x, GEN p, GEN S);

/* Functions giving information on the factorisation. */

/*TODO: merge with FpXQ_powers */
/* u in Z[X], return kernel of (Frob - Id) over Fp[X] / u */
GEN
FpX_Berlekamp_ker(GEN u, GEN p)
{
  long j,N = degpol(u);
  GEN v,w,Q,p1;
  Q = cgetg(N+1,t_MAT); gel(Q,1) = zerocol(N);
  w = v = FpXQ_pow(pol_x[varn(u)],p,u,p);
  for (j=2; j<=N; j++)
  {
    p1 = RgX_to_RgV(w, N);
    gel(p1,j) = addis(gel(p1,j), -1);
    gel(Q,j) = p1;
    if (j < N)
    {
      pari_sp av = avma; /*FpXQ_mul is not stack clean*/
      w = gerepileupto(av, FpXQ_mul(w, v, u, p));
    }
  }
  return FpM_ker(Q,p);
}

GEN
FqX_Berlekamp_ker(GEN u, GEN T, GEN q, GEN p)
{
  pari_sp ltop=avma;
  long j,N = degpol(u);
  GEN v,w,Q,p1;
  pari_timer Ti;
  if (DEBUGLEVEL>=4) TIMERstart(&Ti);
  Q = cgetg(N+1,t_MAT); gel(Q,1) = zerocol(N);
  w = v = FpXQYQ_pow(pol_x[varn(u)], q, u, T, p);
  if (DEBUGLEVEL>=4) msgTIMER(&Ti, "FpXQYQ_pow");
  for (j=2; j<=N; j++)
  {
    p1 = RgX_to_RgV(w, N);
    gel(p1,j) = gaddgs(gel(p1,j), -1);
    gel(Q,j) = p1;
    if (j < N)
    {
      pari_sp av = avma;
      w = gerepileupto(av, FpXQX_divrem(FpXQX_mul(w,v, T,p), u,T,p,ONLY_REM));
    }
  }
  if (DEBUGLEVEL>=4) msgTIMER(&Ti, "Berlekamp_matrix");
  p1 = FqM_ker(Q,T,p);
  if (DEBUGLEVEL>=4) msgTIMER(&Ti, "Berlekamp_ker");
  return gerepileupto(ltop,p1);
}

GEN
Flx_Berlekamp_ker(GEN u, ulong p)
{
  pari_sp ltop=avma;
  long j,N = degpol(u);
  GEN v,w,Q,p1;
  pari_timer T;
  TIMERstart(&T);
  Q = cgetg(N+1,t_VEC); gel(Q,1) = const_vecsmall(N,0);
  w = v = Flxq_pow(polx_Flx(u[1]), utoipos(p), u, p);
  for (j=2; j<=N; j++)
  {
    p1 = Flx_to_Flv(w, N);
    p1[j] = Fl_sub((ulong)p1[j], 1, p);
    gel(Q,j) = p1;
    if (j < N)
    {
      pari_sp av = avma; /*Flxq_mul is not stack clean*/
      w = gerepileupto(av, Flxq_mul(w, v, u, p));
    }
  }
  if(DEBUGLEVEL>=9) msgTIMER(&T,"Berlekamp matrix");
  Q=Flm_ker_sp(Q,p,0);
  if(DEBUGLEVEL>=9) msgTIMER(&T,"kernel");
  return gerepileupto(ltop,Q);
}

GEN
ZX_deriv(GEN x)
{
  long i,lx = lg(x)-1;
  GEN y;

  if (lx<3) return zeropol(varn(x));
  y = cgetg(lx,t_POL);
  for (i=2; i<lx ; i++) gel(y,i) = mulsi(i-1,gel(x,i+1));
  y[1] = x[1]; return y;
}

GEN
FpX_deriv(GEN f, GEN p) { return FpX_red(ZX_deriv(f), p); }
GEN
FqX_deriv(GEN f, /*unused*/GEN T, GEN p) {
  (void)T; return FpXX_red(derivpol(f), p);
}

/* f in ZZ[X] and p a prime number. */
long
FpX_is_squarefree(GEN f, GEN p)
{
  pari_sp av = avma;
  GEN z;
  z = FpX_gcd(f,ZX_deriv(f),p);
  avma = av;
  return lg(z)==3;
}

/* product of terms of degree 1 in factorization of f */
static GEN
FpX_split_part(GEN f, GEN p)
{
  long n = degpol(f);
  GEN z, X = pol_x[varn(f)];
  if (n <= 1) return f;
  f = FpX_red(f, p);
  z = FpXQ_pow(X, p, f, p);
  z = FpX_sub(z, X, p);
  return FpX_gcd(z,f,p);
}

static GEN
FqX_split_part(GEN f, GEN T, GEN p)
{
  long n = degpol(f);
  GEN z, X = pol_x[varn(f)];
  if (n <= 1) return f;
  f = FpXQX_red(f, T, p);
  z = FpXQYQ_pow(X, powiu(p, degpol(T)), f, T, p);
  z = gsub(z, X);
  return FqX_gcd(z, f, T, p);
}

/* Compute the number of roots in Fq without counting multiplicity
 * return -1 for 0 polynomial. lc(f) must be prime to p. */
long
FpX_nbroots(GEN f, GEN p)
{
  pari_sp av = avma;
  GEN z = FpX_split_part(f, p);
  avma = av; return degpol(z);
}

long
FqX_nbroots(GEN f, GEN T, GEN p)
{
  pari_sp av = avma;
  GEN z = FqX_split_part(f, T, p);
  avma = av; return degpol(z);
}

long
FpX_is_totally_split(GEN f, GEN p)
{
  long n=degpol(f);
  pari_sp av = avma;
  GEN z;
  if (n <= 1) return 1;
  if (cmpui(n, p) > 0) return 0;
  f = FpX_red(f, p);
  z = FpXQ_pow(pol_x[varn(f)], p, f, p);
  avma = av;
  return degpol(z) == 1 && gcmp1(gel(z,3)) && !signe(z[2]); /* x^p = x ? */
}

/* Flv_Flx( Flm_Flc_mul(x, Flx_Flv(y), p) ) */
static GEN
Flm_Flx_mul(GEN x, GEN y, ulong p)
{
  long i,k,l, ly = lg(y)-1;
  GEN z;
  long vs=y[1];
  if (ly==1) return zero_Flx(vs);
  l = lg(x[1]);
  y++;
  z = const_vecsmall(l,0) + 1;
  if (u_OK_ULONG(p))
  {
    for (k=1; k<ly; k++)
    {
      GEN c;
      if (!y[k]) continue;
      c = gel(x,k);
      if (y[k] == 1)
        for (i=1; i<l; i++)
        {
          z[i] += c[i];
          if (z[i] & HIGHBIT) z[i] %= p;
        }
      else
        for (i=1; i<l; i++)
        {
          z[i] += c[i] * y[k];
          if (z[i] & HIGHBIT) z[i] %= p;
        }
    }
    for (i=1; i<l; i++) z[i] %= p;
  }
  else
  {
    for (k=1; k<ly; k++)
    {
      GEN c;
      if (!y[k]) continue;
      c = gel(x,k);
      if (y[k] == 1)
        for (i=1; i<l; i++)
          z[i] = Fl_add(z[i], c[i], p);
      else
        for (i=1; i<l; i++)
          z[i] = Fl_add(z[i], Fl_mul(c[i],y[k],p), p);
    }
  }
  while (--l && !z[l]);
  if (!l) return zero_Flx(vs);
  *z-- = vs; return z;
}

/* assume deg(u) > 0 */
static GEN
Flx_Frobenius(GEN u, ulong p)
{
  long j,N = degpol(u);
  GEN v,w,Q;
  pari_timer T;

  if (DEBUGLEVEL > 7) TIMERstart(&T);
  Q = cgetg(N+1,t_MAT);
  gel(Q,1) = const_vecsmall(N, 0);
  coeff(Q,1,1) = 1;
  w = v = Flxq_pow(polx_Flx(u[1]), utoipos(p), u, p);
  for (j=2; j<=N; j++)
  {
    gel(Q,j) = Flx_to_Flv(w, N);
    if (j < N)
    {
      pari_sp av = avma;
      w = gerepileupto(av, Flxq_mul(w, v, u, p));
    }
  }
  if (DEBUGLEVEL > 7) msgTIMER(&T, "frobenius");
  return Q;
}

/* z must be squarefree mod p*/
long
Flx_nbfact(GEN z, ulong p)
{
  long lgg, nfacp = 0, d = 0, e = degpol(z);
  GEN g, w, MP = Flx_Frobenius(z, p), PolX = polx_Flx(z[1]);

  w = PolX;
  while (d < (e>>1))
  { /* here e = degpol(z) */
    d++;
    w = Flm_Flx_mul(MP, w, p); /* w^p mod (z,p) */
    g = Flx_gcd(z, Flx_sub(w, PolX, p), p);
    lgg = degpol(g);
    if (!lgg) continue;

    e -= lgg;
    nfacp += lgg/d;
    if (DEBUGLEVEL>5)
      fprintferr("   %3ld fact. of degree %3ld\n", lgg/d, d);
    if (!e) break;
    z = Flx_div(z, g, p);
    w = Flx_rem(w, z, p);
  }
  if (e)
  {
    if (DEBUGLEVEL>5) fprintferr("   %3ld factor of degree %3ld\n",1,e);
    nfacp++;
  }
  return nfacp;
}

long
Flx_nbroots(GEN f, ulong p)
{
  long n = degpol(f);
  pari_sp av = avma;
  GEN z, X;
  if (n <= 1) return n;
  X = polx_Flx(f[1]);
  z = Flxq_pow(X, utoipos(p), f, p);
  z = Flx_sub(z, X, p);
  z = Flx_gcd(z, f, p);
  avma = av; return degpol(z);
}

long
FpX_nbfact(GEN u, GEN p)
{
  pari_sp av = avma;
  GEN vker = FpX_Berlekamp_ker(u, p);
  avma = av; return lg(vker)-1;
}

long
FqX_nbfact(GEN u, GEN T, GEN p)
{
  pari_sp av = avma;
  GEN vker;
  if (!T) return FpX_nbfact(u, p);
  vker = FqX_Berlekamp_ker(u, T, powiu(p, degpol(T)), p);
  avma = av; return lg(vker)-1;
}

/************************************************************/
GEN
trivfact(void)
{
  GEN y = cgetg(3,t_MAT);
  gel(y,1) = cgetg(1,t_COL);
  gel(y,2) = cgetg(1,t_COL); return y;
}

static GEN
try_pow(GEN w0, GEN pol, GEN p, GEN q, long r)
{
  GEN w2, w = FpXQ_pow(w0,q, pol,p);
  long s;
  if (gcmp1(w)) return w0;
  for (s=1; s<r; s++,w=w2)
  {
    w2 = gsqr(w);
    w2 = FpX_rem(w2, pol, p);
    if (gcmp1(w2)) break;
  }
  return gcmp_1(w)? NULL: w;
}

/* INPUT:
 *  m integer (converted to polynomial w in Z[X] by stopoly)
 *  p prime; q = (p^d-1) / 2^r
 *  t[0] polynomial of degree k*d product of k irreducible factors of degree d
 *  t[0] is expected to be normalized (leading coeff = 1)
 * OUTPUT:
 *  t[0],t[1]...t[k-1] the k factors, normalized */
static void
split(ulong m, GEN *t, long d, GEN p, GEN q, long r, GEN S)
{
  long l, v, dv;
  ulong ps;
  pari_sp av0, av;
  GEN w,w0;

  dv=degpol(*t); if (dv==d) return;
  v=varn(*t); av0=avma; ps = (ulong)p[2];
  for(av=avma;;avma=av)
  {
    if (ps==2)
    {
      w0 = w = FpXQ_pow(pol_x[v], utoi(m-1), *t, gen_2); m += 2;
      for (l=1; l<d; l++)
        w = gadd(w0, spec_FpXQ_pow(w, p, S));
    }
    else
    {
      w = FpX_rem(stopoly(m,ps,v),*t, p);
      m++; w = try_pow(w,*t,p,q,r);
      if (!w) continue;
      w = ZX_Z_add(w, gen_m1);
    }
    w = FpX_gcd(*t,w, p);
    l = degpol(w); if (l && l!=dv) break;
  }
  w = FpX_normalize(w, p);
  w = gerepileupto(av0, w);
  l /= d; t[l]=FpX_div(*t,w,p); *t=w;
  split(m,t+l,d,p,q,r,S);
  split(m,t,  d,p,q,r,S);
}

static void
splitgen(GEN m, GEN *t, long d, GEN  p, GEN q, long r)
{
  long l, v, dv = degpol(*t);
  pari_sp av;
  GEN w;

  if (dv==d) return;
  v = varn(*t);
  m = setloop(m);
  av = avma; 
  for(;; avma = av)
  {
    m = incloop(m);
    w = FpX_rem(stopoly_gen(m,p,v),*t, p);
    w = try_pow(w,*t,p,q,r);
    if (!w) continue;
    w = ZX_Z_add(w, gen_m1);
    w = FpX_gcd(*t,w, p); l=degpol(w);
    if (l && l!=dv) break;
  }
  w = FpX_normalize(w, p);
  w = gerepileupto(av, w);
  l /= d; t[l]=FpX_div(*t,w,p); *t=w;
  splitgen(m,t+l,d,p,q,r);
  splitgen(m,t,  d,p,q,r);
}

/* return S = [ x^p, x^2p, ... x^(n-1)p ] mod (p, T), n = degree(T) > 0 */
static GEN
init_spec_FpXQ_pow(GEN p, GEN T)
{
  long i, n = degpol(T), v = varn(T);
  GEN S = cgetg(n, t_VEC), x;
  if (n == 1) return S;
  x = FpXQ_pow(pol_x[v], p, T, p);
  gel(S,1) = x;
  if ((degpol(x)<<1) < degpol(T)) {
    for (i=2; i < n; i++)
      gel(S,i) = FpXQ_mul(gel(S,i-1), x, T,p);
  } else {
    for (i=2; i < n; i++)
      gel(S,i) = (i&1)? FpXQ_mul(gel(S,i-1), x, T,p)
                      : FpXQ_sqr(gel(S,i>>1), T, p);
  }
  return S;
}

/* compute x^p, S is as above */
static GEN
spec_FpXQ_pow(GEN x, GEN p, GEN S)
{
  long i, dx = degpol(x);
  pari_sp av = avma, lim = stack_lim(av, 1);
  GEN x0 = x+2, z = gel(x0,0);
  if (dx < 0) pari_err(talker, "zero polynomial in FpXQ_pow. %Z not prime", p);
  for (i = 1; i <= dx; i++)
  {
    GEN d, c = gel(x0,i); /* assume coeffs in [0, p-1] */
    if (!signe(c)) continue;
    d = gel(S,i); if (!gcmp1(c)) d = gmul(c,d);
    z = gadd(z, d);
    if (low_stack(lim, stack_lim(av,1)))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"spec_FpXQ_pow");
      z = gerepileupto(av, z);
    }
  }
  return gerepileupto(av, FpX_red(z, p));
}

static int
cmpGsGs(GEN a, GEN b) { return (long)a - (long)b; }

static GEN
FpX_is_irred_2(GEN f, GEN p, long d)
{
  if (!d) return NULL;
  if (d == 1) return gen_1;
  return FpX_quad_factortype(f, p) == -1? gen_1: NULL;
}
static GEN
FpX_degfact_2(GEN f, GEN p, long d)
{
  if (!d) return trivfact();
  if (d == 1) return mkvec2(mkvecsmall(1), mkvecsmall(1));
  switch(FpX_quad_factortype(f, p)) {
    case 1: return mkvec2(mkvecsmall2(1,1), mkvecsmall2(1,1));
    case -1:return mkvec2(mkvecsmall(2), mkvecsmall(1));
    default: return mkvec2(mkvecsmall(1), mkvecsmall(2));
  }
}
static GEN
FpX_factor_2(GEN f, GEN p, long d)
{
  GEN r, s, R, S;
  long v;
  int sgn;
  if (d < 0) pari_err(zeropoler,"FpX_factor_2");
  if (!d) return trivfact();
  if (d == 1) return mkmat2(mkcol(f), mkvecsmall(1));
  r = FpX_quad_root(f, p, 1);
  if (!r) return mkmat2(mkcol(f), mkvecsmall(1));
  v = varn(f);
  s = otherroot(f, r, p);
  if (signe(r)) r = subii(p, r);
  if (signe(s)) s = subii(p, s);
  sgn = cmpii(s, r); if (sgn < 0) swap(s,r);
  R = deg1pol_i(gen_1, r, v);
  if (!sgn) return mkmat2(mkcol(R), mkvecsmall(2));
  S = deg1pol_i(gen_1, s, v);
  return mkmat2(mkcol2(R,S), mkvecsmall2(1,1));
}

/* factor f mod pp.
 * flag = 1: return the degrees, not the factors
 * flag = 2: return NULL if f is not irreducible */
static GEN
FpX_factcantor_i(GEN f, GEN pp, long flag)
{
  long j, e, vf, nbfact, d = degpol(f);
  ulong p, k;
  GEN E,y,f2,g,g1,u,v,pd,q;
  GEN *t;

  if (d <= 2) switch(flag) {
    case 2: return FpX_is_irred_2(f, pp, d);
    case 1: return FpX_degfact_2(f, pp, d);
    default: return FpX_factor_2(f, pp, d);
  }
  p = init_p(pp);

  /* to hold factors and exponents */
  t = (GEN*)cgetg(d+1,t_VEC);
  E = cgetg(d+1, t_VECSMALL);
  vf=varn(f); e = nbfact = 1;
  for(;;)
  {
    f2 = FpX_gcd(f,ZX_deriv(f), pp);
    if (flag > 1 && lg(f2) > 3) return NULL;
    g1 = FpX_div(f,f2,pp);
    k = 0;
    while (lg(g1)>3)
    {
      long du,dg;
      GEN S;
      k++; if (p && !(k%p)) { k++; f2 = FpX_div(f2,g1,pp); }
      u = g1; g1 = FpX_gcd(f2,g1, pp);
      if (lg(g1)>3)
      {
        u = FpX_div( u,g1,pp);
        f2= FpX_div(f2,g1,pp);
      }
      du = degpol(u);
      if (du <= 0) continue;

      /* here u is square-free (product of irred. of multiplicity e * k) */
      S = init_spec_FpXQ_pow(pp, u);
      pd=gen_1; v=pol_x[vf];
      for (d=1; d <= du>>1; d++)
      {
        if (!flag) pd = mulii(pd,pp);
        v = spec_FpXQ_pow(v, pp, S);
        g = FpX_gcd(gadd(v, gneg(pol_x[vf])), u, pp);
        dg = degpol(g);
        if (dg <= 0) continue;

        /* g is a product of irred. pols, all of which have degree d */
        j = nbfact+dg/d;
        if (flag)
        {
          if (flag > 1) return NULL;
          for ( ; nbfact<j; nbfact++) { t[nbfact]=(GEN)d; E[nbfact]=e*k; }
        }
        else
        {
          long r;
          g = FpX_normalize(g, pp);
          t[nbfact]=g; q = subis(pd,1); /* also ok for p=2: unused */
          r = vali(q); q = shifti(q,-r);
         /* First parameter is an integer m, converted to polynomial w_m, whose
          * coeffs are its digits in base p (initially m = p --> w_m = X). Take
          * gcd(g, w_m^(p^d-1)/2), m++, until a factor is found. p = 2 is
          * treated separately */
          if (p)
            split(p,t+nbfact,d,pp,q,r,S);
          else
            splitgen(pp,t+nbfact,d,pp,q,r);
          for (; nbfact<j; nbfact++) E[nbfact]=e*k;
        }
        du -= dg;
        u = FpX_div(u,g,pp);
        v = FpX_rem(v,u,pp);
      }
      if (du)
      {
        t[nbfact] = flag? (GEN)du: FpX_normalize(u, pp);
        E[nbfact++]=e*k;
      }
    }
    j = lg(f2); if (j==3) break;
    e *= p; f = poldeflate_i(f2, p);
  }
  if (flag > 1) return gen_1; /* irreducible */
  setlg(t, nbfact);
  setlg(E, nbfact); y = mkvec2((GEN)t, E);
  return flag ? sort_factor_gen(y, cmpGsGs)
              : sort_factor(y, cmpii);
}
GEN
FpX_factcantor(GEN f, GEN pp, long flag)
{
  pari_sp av = avma;
  GEN z = FpX_factcantor_i(FpX_factmod_init(f,pp),pp,flag);
  if (flag == 2) { avma = av; return z; }
  return gerepileupto(av, z);
}
GEN
factcantor0(GEN f, GEN pp, long flag)
{
  pari_sp av = avma;
  long j, nbfact;
  GEN z, y, t, E, u, v;
  if (! factmod_init(&f, pp)) { avma=av; return trivfact(); }
  z = FpX_factcantor_i(f,pp,flag); t = gel(z,1); E = gel(z,2);
  y = cgetg(3, t_MAT); nbfact = lg(t);
  u = cgetg(nbfact,t_COL); gel(y,1) = u;
  v = cgetg(nbfact,t_COL); gel(y,2) = v;
  if (flag)
    for (j=1; j<nbfact; j++)
    {
      gel(u,j) = utoi((ulong)t[j]);
      gel(v,j) = utoi((ulong)E[j]);
    }
  else
    for (j=1; j<nbfact; j++)
    {
      gel(u,j) = FpX_to_mod(gel(t,j), pp);
      gel(v,j) = utoi((ulong)E[j]);
    }
  return gerepileupto(av, y);
}

/* Use this function when you think f is reducible, and that there are lots of
 * factors. If you believe f has few factors, use FpX_nbfact(f,p)==1 instead */
long
FpX_is_irred(GEN f, GEN p) {
  return !!FpX_factcantor_i(FpX_factmod_init(f,p),p,2);
}
GEN
FpX_degfact(GEN f, GEN p) {
  pari_sp av = avma;
  GEN z = FpX_factcantor_i(FpX_factmod_init(f,p),p,1);
  settyp(z[1], t_VECSMALL);
  settyp(z, t_MAT); return gerepilecopy(av, z);
}
GEN
factcantor(GEN f, GEN p) { return factcantor0(f,p,0); }
GEN
simplefactmod(GEN f, GEN p) { return factcantor0(f,p,1); }

GEN
col_to_ff(GEN x, long v)
{
  long i, k = lg(x);
  GEN p;

  while (--k && gcmp0(gel(x,k)));
  if (k <= 1) return k? gel(x,1): gen_0;
  i = k+2; p = cgetg(i,t_POL);
  p[1] = evalsigne(1) | evalvarn(v);
  x--; for (k=2; k<i; k++) p[k] = x[k];
  return p;
}

/* set x <-- x + c*y mod p */
/* x is not required to be normalized.*/
static void
Flx_addmul_inplace(GEN gx, GEN gy, ulong c, ulong p)
{
  long i, lx, ly;
  ulong *x=(ulong *)gx;
  ulong *y=(ulong *)gy;
  if (!c) return;
  lx = lg(gx);
  ly = lg(gy);
  if (lx<ly) pari_err(bugparier,"lx<ly in Flx_addmul_inplace");
  if (u_OK_ULONG(p))
    for (i=2; i<ly;  i++) x[i] = (x[i] + c*y[i]) % p;
  else
    for (i=2; i<ly;  i++) x[i] = Fl_add(x[i], Fl_mul(c,y[i],p),p);
}

static long
small_rand(ulong p)
{
  return (p==2)? ((pari_rand31() & 0x1000) == 0): pari_rand31() % p;
}

GEN
FpX_rand(long d1, long v, GEN p)
{
  long i, d = d1+2;
  GEN y = cgetg(d,t_POL); y[1] = evalsigne(1) | evalvarn(v);
  for (i=2; i<d; i++) gel(y,i) = genrand(p);
  (void)normalizepol_i(y,d); return y;
}

/* return a random polynomial in F_q[v], degree < d1 */
GEN
FqX_rand(long d1, long v, GEN T, GEN p)
{
  long i, d = d1+2, k = degpol(T), w = varn(T);
  GEN y = cgetg(d,t_POL); y[1] = evalsigne(1) | evalvarn(v);
  for (i=2; i<d; i++) gel(y,i) = FpX_rand(k, w, p);
  (void)normalizepol_i(y,d); return y;
}

#define set_irred(i) { if ((i)>ir) swap(t[i],t[ir]); ir++;}

long
FpX_split_Berlekamp(GEN *t, GEN p)
{
  GEN u = *t, a,b,po2,vker;
  long d, i, ir, L, la, lb, vu = varn(u);
  long l = lg(u);
  ulong ps = itou_or_0(p);
  if (ps)
  {
    vker = Flx_Berlekamp_ker(ZX_to_Flx(u,ps),ps);
    vker = Flm_to_FlxV(vker, u[1]);
  }
  else
  {
    vker = FpX_Berlekamp_ker(u,p);
    vker = RgM_to_RgXV(vker,vu);
  }
  d = lg(vker)-1;
  po2 = shifti(p, -1); /* (p-1) / 2 */
  ir = 0;
  /* t[i] irreducible for i <= ir, still to be treated for i < L */
  for (L=1; L<d; )
  {
    GEN polt;
    if (ps)
    {
      GEN pol = const_vecsmall(l-2,0);
      pol[1] = u[1];
      pol[2] = small_rand(ps); /*Assume vker[1]=1*/
      for (i=2; i<=d; i++)
        Flx_addmul_inplace(pol, gel(vker,i), (ulong)small_rand(ps), ps);
      (void)Flx_renormalize((GEN)pol,l-1);

      polt = Flx_to_ZX(pol);
    }
    else
    {
      GEN pol = scalarpol(genrand(p), vu);
      for (i=2; i<=d; i++)
        pol = ZX_add(pol, ZX_Z_mul(gel(vker,i), randomi(p)));
      polt = FpX_red(pol,p);
    }
    for (i=ir; i<L && L<d; i++)
    {
      a = t[i]; la = degpol(a);
      if (la == 1) { set_irred(i); }
      else if (la == 2)
      {
        GEN r = FpX_quad_root(a,p,1);
        if (r)
        {
          t[i] = deg1pol_i(gen_1, subii(p,r), vu); r = otherroot(a,r,p);
          t[L] = deg1pol_i(gen_1, subii(p,r), vu); L++;
        }
        set_irred(i);
      }
      else
      {
        pari_sp av = avma;
        b = FpX_rem(polt, a, p);
        if (degpol(b) <= 0) { avma=av; continue; }
        b = ZX_Z_add(FpXQ_pow(b,po2, a,p), gen_m1);
        b = FpX_gcd(a,b, p); lb = degpol(b);
        if (lb && lb < la)
        {
          b = FpX_normalize(b, p);
          t[L] = FpX_div(a,b,p);
          t[i]= b; L++;
        }
        else avma = av;
      }
    }
  }
  return d;
}

GEN
FqX_gcd(GEN P,GEN Q,GEN T,GEN p) {return T? FpXQX_gcd(P,Q,T,p): FpX_gcd(P,Q,p);}
long
FqX_is_squarefree(GEN P, GEN T, GEN p)
{
  pari_sp av = avma;
  GEN z = FqX_gcd(P, FqX_deriv(P, T, p), T, p);
  avma = av;
  return degpol(z)==0;
}

long
FqX_split_Berlekamp(GEN *t, GEN q, GEN T, GEN p)
{
  GEN u = *t, a,b,qo2,vker,pol;
  long N = degpol(u), vu = varn(u), vT = varn(T), dT = degpol(T);
  long d, i, ir, L, la, lb;

  vker = FqX_Berlekamp_ker(u,T,q,p);
  vker = RgM_to_RgXV(vker,vu);
  d = lg(vker)-1;
  qo2 = shifti(q, -1); /* (q-1) / 2 */
  pol = cgetg(N+3,t_POL);
  ir = 0;
  /* t[i] irreducible for i < ir, still to be treated for i < L */
  for (L=1; L<d; )
  {
    GEN polt;
    gel(pol,2) = FpX_rand(dT,vT,p);
    setlg(pol, signe(pol[2])? 3: 2);
    pol[1] = u[1];
    for (i=2; i<=d; i++)
      pol = gadd(pol, gmul(gel(vker,i), FpX_rand(dT,vT,p)));
    polt = FpXQX_red(pol,T,p);
    for (i=ir; i<L && L<d; i++)
    {
      a = t[i]; la = degpol(a);
      if (la == 1) { set_irred(i); }
      else
      {
        pari_sp av = avma;
        b = FqX_rem(polt, a, T,p);
        if (!degpol(b)) { avma=av; continue; }
        b = FpXQYQ_pow(b,qo2, a,T,p);
        if (!degpol(b)) { avma=av; continue; }
        gel(b,2) = gadd(gel(b,2), gen_1);
        b = FqX_gcd(a,b, T,p); lb = degpol(b);
        if (lb && lb < la)
        {
          b = FqX_normalize(b, T,p);
          t[L] = FqX_div(a,b,T,p);
          t[i]= b; L++;
        }
        else avma = av;
      }
    }
  }
  return d;
}

static GEN
FpX_factor_i(GEN f, GEN pp)
{
  long e, N, nbfact, val, d = degpol(f);
  ulong p, k, j;
  GEN E, f2, g1, u, *t;

  if (d <= 2) return FpX_factor_2(f, pp, d);
  p = init_p(pp);

  /* to hold factors and exponents */
  t = (GEN*)cgetg(d+1,t_COL); E = cgetg(d+1,t_VECSMALL);
  val = ZX_valuation(f, &f);
  e = nbfact = 1;
  if (val) { t[1] = pol_x[varn(f)]; E[1] = val; nbfact++; }

  for(;;)
  {
    f2 = FpX_gcd(f,ZX_deriv(f), pp);
    g1 = lg(f2)==3? f: FpX_div(f,f2,pp); /* is squarefree */
    k = 0;
    while (lg(g1)>3)
    {
      k++; if (p && !(k%p)) { k++; f2 = FpX_div(f2,g1,pp); }
      u = g1; 
      if (!degpol(f2)) g1 = pol_1[0]; /* only its degree (= 0) matters */
      else
      {
        g1 = FpX_gcd(f2,g1, pp);
        if (degpol(g1))
        {
          u = FpX_div( u,g1,pp);
          f2= FpX_div(f2,g1,pp);
        }
      }
      /* u is square-free (product of irred. of multiplicity e * k) */
      N = degpol(u);
      if (N > 0)
      {
        t[nbfact] = FpX_normalize(u,pp);
        d = (N==1)? 1: FpX_split_Berlekamp(t+nbfact, pp);
        for (j=0; j<(ulong)d; j++) E[nbfact+j] = e*k;
        nbfact += d;
      }
    }
    if (!p) break;
    j = degpol(f2); if (!j) break;
    if (j % p) pari_err(talker, "factmod: %lu is not prime", p);
    e *= p; f = poldeflate_i(f2, p);
  }
  setlg(t, nbfact);
  setlg(E, nbfact); return sort_factor(mkvec2((GEN)t,E), cmpii);
}
GEN
FpX_factor(GEN f, GEN p)
{
  pari_sp av = avma;
  return gerepilecopy(av, FpX_factor_i(FpX_factmod_init(f, p), p));
}

GEN
factmod(GEN f, GEN p)
{
  pari_sp av = avma;
  long nbfact;
  long j;
  GEN y, u, v, z, t, E;

  if (!factmod_init(&f, p)) { avma = av; return trivfact(); }
  z = FpX_factor_i(f, p); t = gel(z,1); E = gel(z,2);
  y = cgetg(3,t_MAT); nbfact = lg(t);
  u = cgetg(nbfact,t_COL); gel(y,1) = u;
  v = cgetg(nbfact,t_COL); gel(y,2) = v;
  for (j=1; j<nbfact; j++)
  {
    gel(u,j) = FpX_to_mod(gel(t,j), p);
    gel(v,j) = utoi((ulong)E[j]);
  }
  return gerepileupto(av, y);
}
GEN
factormod0(GEN f, GEN p, long flag)
{
  switch(flag)
  {
    case 0: return factmod(f,p);
    case 1: return simplefactmod(f,p);
    default: pari_err(flagerr,"factormod");
  }
  return NULL; /* not reached */
}

/*******************************************************************/
/*                                                                 */
/*                  CONVERSIONS RELATED TO p-ADICS                 */
/*                                                                 */
/*******************************************************************/
static GEN
Zp_to_Z(GEN x) {
  switch(typ(x))
  {
    case t_INT: break;
    case t_PADIC: x = gtrunc(x); break;
    default: pari_err(typeer,"QpX_to_ZX");
  }
  return x;
}
static GEN
ZpX_to_ZX(GEN f) {
  long i, l = lg(f);
  GEN F = cgetg(l, t_POL); F[1] = f[1];
  for (i=2; i<l; i++) gel(f,i) = Zp_to_Z(gel(f,i));
  return f;
}
/* make f suitable for [root|factor]padic */
static GEN
QpX_to_ZX(GEN f)
{
  GEN c = content(f);
  if (gcmp0(c)) /*  O(p^n) can occur */
  {
    if (typ(c) != t_PADIC) pari_err(typeer,"QpX_to_ZX");
    f = gdiv(f, gpowgs(gel(c,2), valp(c)));
  }
  else
    f = gdiv(f,c);
  return ZpX_to_ZX(f);
}

/* x in Z return x + O(pr), pr = p^r. Return gen_0 instead of zeropadic */
static GEN
Z_to_Zp(GEN x, GEN p, GEN pr, long r)
{
  GEN y;
  long v, sx = signe(x);

  if (!sx) return gen_0;
  v = Z_pvalrem(x,p,&x);
  r -= v; if (r <= 0) return gen_0;
  y = cgetg(5,t_PADIC);
  y[1] = evalprecp(r)|evalvalp(v);
  gel(y,2) = p;
  gel(y,3) = pr;
  gel(y,4) = modii(x,pr); return y;
}
static GEN
ZV_to_ZpV(GEN z, GEN p, long prec)
{
  long i, l = lg(z);
  GEN Z = cgetg(l, typ(z)), q = powiu(p, prec);
  for (i=1; i<lg(z); i++) gel(Z,i) = Z_to_Zp(gel(z,i),p,q,prec);
  return Z;
}
static GEN
ZX_to_ZpX(GEN z, GEN p, GEN q, long prec)
{
  long i, l = lg(z);
  GEN Z = cgetg(l, t_POL); Z[1] = z[1];
  for (i=2; i<lg(z); i++) gel(Z,i) = Z_to_Zp(gel(z,i),p,q,prec);
  return Z;
}
/* return (x + O(p^r)) normalized (multiply by a unit such that leading coeff
 * is a power of p), x in Z[X] (or Z_p[X]) */
static GEN
ZX_to_ZpX_normalized(GEN x, GEN p, GEN pr, long r)
{
  long i, lx = lg(x);
  GEN z, lead = leading_term(x);

  if (gcmp1(lead)) return ZX_to_ZpX(x, p, pr, r);
  (void)Z_pvalrem(lead, p, &lead); lead = Fp_inv(lead, pr);
  z = cgetg(lx,t_POL);
  for (i=2; i < lx; i++) gel(z,i) = Z_to_Zp(mulii(lead,gel(x,i)),p,pr,r);
  z[1] = x[1]; return z;
}
static GEN
ZXV_to_ZpXQV(GEN z, GEN T, GEN p, long prec)
{
  long i, l = lg(z);
  GEN Z = cgetg(l, typ(z)), q = powiu(p, prec);
  T = gcopy(T);
  for (i=1; i<lg(z); i++) gel(Z,i) = mkpolmod(ZX_to_ZpX(gel(z,i),p,q,prec),T);
  return Z;
}

static GEN
QpXQ_to_ZXY(GEN f)
{
  GEN c = content(f);
  long i, l = lg(f);
  if (gcmp0(c)) /*  O(p^n) can occur */
  {
    if (typ(c) != t_PADIC) pari_err(typeer,"QpXQ_to_ZXY");
    f = gdiv(f, gpowgs(gel(c,2), valp(c)));
  }
  else
    f = gdiv(f,c);
  for (i=2; i<l; i++)
  {
    GEN t = gel(f,i);
    switch(typ(t))
    {
      case t_POL: t = ZpX_to_ZX(t); break;
      default: t = Zp_to_Z(t); break;
    }
    gel(f,i) = t;
  }
  return f;
}

/*******************************************************************/
/*                                                                 */
/*                         p-ADIC ROOTS                            */
/*                                                                 */
/*******************************************************************/

/* f primitive ZX, squarefree, leading term prime to p. a in Z such that
 * f(a) = 0 mod p. Return p-adic roots of f equal to a mod p, in
 * precision >= prec */
static GEN
ZX_Zp_root(GEN f, GEN a, GEN p, long prec)
{
  GEN z, R, a0 = modii(a, p);
  long i, j, k;

  if (signe(FpX_eval(FpX_deriv(f, p), a0, p)))
  { /* simple zero mod p, go all the way to p^prec */
    if (prec > 1) a0 = ZpX_liftroot(f, a0, p, prec);
    return mkcol(a0);
  }

  f = poleval(f, gadd(a, gmul(p,pol_x[varn(f)])));
  f = gdivexact(f, powiu(p,ggval(f, p)));
  z = cgetg(degpol(f)+1,t_COL);

  R = FpX_roots(f, p);
  for (j=i=1; i<lg(R); i++)
  {
    GEN u = ZX_Zp_root(f, gel(R,i), p, prec-1);
    for (k=1; k<lg(u); k++) gel(z,j++) = gadd(a, gmul(p, gel(u,k)));
  }
  setlg(z,j); return z;
}

/* a t_PADIC, return vector of p-adic roots of f equal to a (mod p)
 * We assume 1) f(a) = 0 mod p (mod 4 if p=2).
 *           2) leading coeff prime to p. */
GEN
Zp_appr(GEN f, GEN a)
{
  pari_sp av = avma;
  long prec;
  GEN z, p;
  if (typ(f) != t_POL) pari_err(notpoler,"Zp_appr");
  if (gcmp0(f)) pari_err(zeropoler,"Zp_appr");
  if (typ(a) != t_PADIC) pari_err(typeer,"Zp_appr");
  p = gel(a,2); prec = gcmp0(a)? valp(a): precp(a);
  f = QpX_to_ZX(f);
  z = modulargcd(f, ZX_deriv(f));
  if (degpol(z) > 0) f = RgX_div(f,z);
  z = ZX_Zp_root(f, gtrunc(a), p, prec);
  return gerepilecopy(av, ZV_to_ZpV(z, p, prec));
}
/* vector of p-adic roots of the ZX f, leading term prime to p */
static GEN
ZX_Zp_roots(GEN f, GEN p, long prec)
{
  GEN y, z, rac;
  long lx, i, j, k;

  z = modulargcd(f, ZX_deriv(f));
  if (degpol(z) > 0) f = RgX_div(f,z);
  rac = FpX_roots(f, p);
  lx = lg(rac); if (lx == 1) return rac;
  y = cgetg(degpol(f)+1,t_COL);
  for (j=i=1; i<lx; i++)
  {
    z = ZX_Zp_root(f, gel(rac,i), p, prec);
    for (k=1; k<lg(z); k++,j++) gel(y,j) = gel(z,k);
  }
  setlg(y,j); return ZV_to_ZpV(y, p, prec);
}

/* f a ZX */
static GEN
pnormalize(GEN f, GEN p, long prec, long n, GEN *plead, long *pprec, int *prev)
{
  *plead = leading_term(f);
  *pprec = prec;
  *prev = 0;
  if (!is_pm1(*plead))
  {
    long v = ggval(*plead,p), v1 = ggval(constant_term(f),p);
    if (v1 < v)
    {
      *prev = 1; f = polrecip_i(f);
     /* beware loss of precision from lc(factor), whose valuation is <= v */
      *pprec += v; v = v1;
    }
    *pprec += v * n;
  }
  return pol_to_monic(f, plead);
}

/* return p-adic roots of f, precision prec */
GEN
rootpadic(GEN f, GEN p, long prec)
{
  pari_sp av = avma;
  GEN lead,y;
  long PREC,i,k;
  int reverse;

  if (typ(p)!=t_INT) pari_err(typeer,"rootpadic");
  if (typ(f)!=t_POL) pari_err(notpoler,"rootpadic");
  if (gcmp0(f)) pari_err(zeropoler,"rootpadic");
  if (prec <= 0) pari_err(talker,"non-positive precision in rootpadic");
  f = QpX_to_ZX(f);
  f = pnormalize(f, p, prec, 1, &lead, &PREC, &reverse);
  y = ZX_Zp_roots(f,p,PREC);
  k = lg(y);
  if (lead)
    for (i=1; i<k; i++) gel(y,i) = gdiv(gel(y,i), lead);
  if (reverse)
    for (i=1; i<k; i++) gel(y,i) = ginv(gel(y,i));
  return gerepilecopy(av, y);
}
/*************************************************************************/
/*                             rootpadicfast                             */
/*************************************************************************/

/* lift accelerator */
long
hensel_lift_accel(long n, long *pmask)
{
  long a, j, mask = 0;
  for(j=BITS_IN_LONG-1, a=n ;; j--)
  {
    mask |= (a&1)<<j;
    a = (a>>1) + (a&1);
    if (a==1) break;
  }
  *pmask = mask>>j;
  return BITS_IN_LONG-j;
}
/* SPEC:
 * p is a t_INT > 1, e >= 0
 * f is a ZX with leading term prime to p.
 * a is a simple root mod l for all l|p.
 * Return roots of f mod p^e, as integers (implicitly mod p^e)
 * STANDARD USE: p is a prime power */
GEN
ZpX_liftroot(GEN f, GEN a, GEN p, long e)
{
  pari_sp ltop=avma;
  GEN     qold, q, qm1;
  GEN     W, fr, ar, Wr = gen_0;
  long    i, nb, mask;
  qold = q = p; qm1 = gen_1;
  nb = hensel_lift_accel(e, &mask);
  fr = FpX_red(f,q);
  a = modii(a,q);
  W = FpX_eval(ZX_deriv(fr), a, q);
  W = Fp_inv(W,q);
  for(i=0;i<nb;i++)
  {
    qm1 = (mask&(1<<i))?sqri(qm1):mulii(qm1, q);
    q   =  mulii(qm1, p);
    fr = FpX_red(f,q);
    ar = a;
    if (i)
    {
      W = modii(mulii(Wr,FpX_eval(ZX_deriv(fr),ar,qold)), qold);
      W = modii(mulii(Wr, subsi(2,W)), qold);
    }
    Wr = W;
    a = subii(ar, mulii(Wr, FpX_eval(fr, ar,q)));
    a = modii(a,q);
    qold = q;
  }
  return gerepileupto(ltop,a);
}
GEN
ZpXQX_liftroot(GEN f, GEN a, GEN T, GEN p, long e)
{
  pari_sp ltop=avma;
  GEN     qold, q, qm1;
  GEN     W, fr, ar, Wr = gen_0;
  long    i, nb, mask;
  qold = p ;  q = p; qm1 = gen_1;
  nb=hensel_lift_accel(e, &mask);
  fr = FpXQX_red(f, T, q);
  a = Fq_red(a, T, q);
  W = FqX_eval(derivpol(fr), a, T, q);
  W = Fq_inv(W,T,q);
  for(i=0;i<nb;i++)
  {
    qm1 = (mask&(1<<i))?sqri(qm1):mulii(qm1, q);
    q   =  mulii(qm1, p);
    fr = FpXQX_red(f,T,q);
    ar = a;
    if (i)
    {
      W = Fq_red(gmul(Wr,FqX_eval(derivpol(fr),ar,T,qold)), T, qold);
      W = Fq_red(gmul(Wr, gadd(gen_2, gneg(W))), T, qold);
    }
    Wr = W;
    a = gadd(ar, gmul(gneg(Wr), FqX_eval(fr, ar, T, q)));
    a = Fq_red(a, T, q);
    qold = q;
  }
  return gerepileupto(ltop,a);
}
/* Apply ZpX_liftroot to all roots in S and trace trick.
 * Elements of S must be distinct simple roots mod p for all p|q. */
GEN
ZpX_liftroots(GEN f, GEN S, GEN q, long e)
{
  long i, d, l = lg(S), n = l-1;
  GEN y = cgetg(l, typ(S));
  if (!n) return y;
  for (i=1; i<n; i++)
    gel(y,i) = ZpX_liftroot(f, gel(S,i), q, e);
  d = degpol(f);
  if (n != d) /* not totally split*/
    gel(y,n) = ZpX_liftroot(f, gel(S,n), q, e);
  else
  { /* totally split: use trace trick */
    pari_sp av = avma;
    GEN z = gel(f, d+1);/* -trace(roots) */
    for(i=1; i<n;i++) z = addii(z, gel(y,i));
    z = modii(negi(z), powiu(q,e));
    gel(y,n) = gerepileuptoint(av,z);
  }
  return y;
}
/* p is prime
 * f in a ZX, with leading term prime to p.
 * f must have no multiple roots mod p.
 *
 * return p-adics roots of f with prec p^e, as integers (implicitly mod p^e) */
GEN
rootpadicfast(GEN f, GEN p, long e)
{
  pari_sp av = avma;
  GEN y, S = FpX_roots(f, p); /*no multiplicity*/
  if (lg(S)==1) { avma = av; return cgetg(1,t_COL); }
  S = gclone(S); avma = av;
  y = ZpX_liftroots(f,S,p,e);
  gunclone(S); return y;
}
/* Same as ZpX_liftroot for the polynomial X^n-T
 * TODO: generalize to sparse polynomials. */
GEN
padicsqrtnlift(GEN T, GEN n, GEN a, GEN p, long e)
{
  pari_sp ltop=avma;
  GEN     qold, q, qm1;
  GEN     W, ar, Wr = gen_0;
  long    i, nb, mask;
  qold = q = p; qm1 = gen_1;
  nb = hensel_lift_accel(e, &mask);
  W = modii(mulii(n,Fp_pow(a,subis(n,1),q)),q);
  W = Fp_inv(W,q);
  for(i=0;i<nb;i++)
  {
    qm1 = (mask&(1<<i))?sqri(qm1):mulii(qm1, q);
    q   =  mulii(qm1, p);
    ar = a;
    if (i)
    {
      W = modii(mulii(Wr,mulii(n,Fp_pow(ar,subis(n,1),qold))),qold);
      W = modii(mulii(Wr, subsi(2,W)),qold);
    }
    Wr = W;
    a = subii(ar, mulii(Wr, subii(Fp_pow(ar,n,q),T)));
    a = modii(a,q);
    qold = q;
  }
  return gerepileupto(ltop,a);
}
/**************************************************************************/

static void
scalar_getprec(GEN x, long *pprec, GEN *pp)
{
  if (typ(x)==t_PADIC)
  {
    long e = valp(x); if (signe(x[4])) e += precp(x);
    if (e < *pprec) *pprec = e;
    if (*pp && !equalii(*pp, gel(x,2))) pari_err(consister,"apprpadic");
    *pp = gel(x,2);
  }
}

static void
getprec(GEN x, long *pprec, GEN *pp)
{
  long i;
  if (typ(x) != t_POL) scalar_getprec(x, pprec, pp);
  else
    for (i = lg(x)-1; i>1; i--) scalar_getprec(gel(x,i), pprec, pp);
}
static GEN
ZXY_ZpQ_root(GEN f, GEN a, GEN T, GEN p, long prec)
{
  GEN z, R;
  long i, j, k, lR;
  if (signe(FqX_eval(FqX_deriv(f,T,p), a, T,p)))
  { /* simple zero mod (T,p), go all the way to p^prec */
    if (prec > 1) a = ZpXQX_liftroot(f, a, T, p, prec);
    return mkcol(a);
  }
  /* TODO: need RgX_RgYQX_compo ? */
  f = lift_intern(poleval(f, gadd(mkpolmod(a,T), gmul(p, pol_x[varn(f)]))));
  f = gdiv(f, powiu(p, ggval(f,p)));
  z = cgetg(degpol(f)+1,t_COL);

#if 1 /* TODO: need a public FqX_roots */
  lR = FqX_split_deg1(&R, FqX_red(f, T, p), powiu(p, degpol(T)), T, p) + 1;
  R = roots_from_deg1(FqX_split_roots(R, T, p, NULL));
#else
  R = FqX_factor(FqX_red(f, T, p), T, p);
  R = gel(R,1); lR = lg(R);
  for (i=1; i<lR; i++)
  {
    GEN u = gel(R,i); if (degpol(u) > 1) break;
    gel(R,i) = gneg(gel(u,2));
  }
  lR = i; /* "truncate" R to the deg 1 factors, demoted to roots above */
#endif
  for(j=i=1; i<lR; i++)
  {
    GEN u = ZXY_ZpQ_root(f, gel(R,i), T, p, prec-1);
    for (k=1; k<lg(u); k++) gel(z,j++) = gadd(a, gmul(p, gel(u,k)));
  }
  setlg(z,j); return z;
}

/* a belongs to finite extension of Q_p, return all roots of f equal to a
 * mod p. We assume f(a) = 0 (mod p) [mod 4 if p=2] */
GEN
padicappr(GEN f, GEN a)
{
  GEN p, z, T;
  long prec;
  pari_sp av = avma;

  switch(typ(a)) {
    case t_PADIC: return Zp_appr(f,a);
    case t_POLMOD: break;
    default: pari_err(typeer,"padicappr");
  }
  if (typ(f)!=t_POL) pari_err(notpoler,"padicappr");
  if (gcmp0(f)) pari_err(zeropoler,"padicappr");
  z = ggcd(f, derivpol(f));
  if (degpol(z) > 0) f = RgX_div(f,z);
  T = gel(a,1); a = gel(a,2);
  p = NULL; prec = BIGINT;
  getprec(a, &prec, &p);
  getprec(T, &prec, &p); if (!p) pari_err(typeer,"padicappr");
  f = QpXQ_to_ZXY(lift_intern(f));
  a = QpX_to_ZX(a);
  T = QpX_to_ZX(T);
  z = ZXY_ZpQ_root(f, a, T, p, prec);
  return gerepilecopy(av, ZXV_to_ZpXQV(z, T, p, prec));
}

/*******************************************************************/
/*                                                                 */
/*                      FACTORIZATION in Zp[X]                     */
/*                                                                 */
/*******************************************************************/
int
cmp_padic(GEN x, GEN y)
{
  long vx, vy;
  if (x == gen_0) return -1;
  if (y == gen_0) return  1;
  vx = valp(x);
  vy = valp(y);
  if (varncmp(vx, vy) < 0) return  1;
  if (varncmp(vx, vy) > 0) return -1;
  return cmpii(gel(x,4), gel(y,4));
}

/*******************************/
/*   Using Buchmann--Lenstra   */
/*******************************/

/* factor T = nf[1] in Zp to precision k */
static GEN
padicff2(GEN nf,GEN p,long k)
{
  GEN z, mat, D, U,fa, pk, dec_p, Ui, mulx;
  long i, l;

  mulx = eltmul_get_table(nf, gmael(nf,8,2)); /* mul by (x mod T) */
  dec_p = primedec(nf,p);
  l = lg(dec_p); fa = cgetg(l,t_COL);
  D = NULL; /* -Wall */
  for (i=1; i<l; i++)
  {
    GEN P = gel(dec_p,i);
    long e = itos(gel(P,3)), ef = e * itos(gel(P,4));
    D = smithall(idealpows(nf,P, k*e), &U, NULL);
    Ui= ginv(U); setlg(Ui, ef+1); /* cf smithrel */
    U = rowslice(U, 1, ef);
    mat = gmul(U, gmul(mulx, Ui));
    gel(fa,i) = caradj(mat,0,NULL);
  }
  pk = gcoeff(D,1,1); /* = p^k */
  z = cgetg(l,t_COL); pk = icopy(pk);
  for (i=1; i<l; i++)
    gel(z,i) = ZX_to_ZpX_normalized(gel(fa,i),p,pk,k);
  return z;
}

static GEN
padicff(GEN x,GEN p,long pr)
{
  pari_sp av = avma;
  GEN q, bas, invbas, mul, dK, nf, g, e, dx = absi(ZX_disc(x));
  long n = degpol(x), v = Z_pvalrem(dx,p,&q);

  nf = cgetg(10,t_VEC); gel(nf,1) = x;
  if (is_pm1(q)) {
    e = mkcol(utoi(v));
    g = mkcol(p);
  } else {
    e = mkcol2(utoi(v), gen_1);
    g = mkcol2(p, q);
  }
  bas = nfbasis(x, &dK, 0, mkmat2(g,e));
  gel(nf,3) = dK;
  gel(nf,4) = dvdii( diviiexact(dx, dK), p )? p: gen_1;
  invbas = QM_inv(RgXV_to_RgM(bas,n), gen_1);
  mul = get_mul_table(x,bas,invbas);
  gel(nf,7) = bas;
  gel(nf,8) = invbas;
  gel(nf,9) = mul;
  gel(nf,2) = gel(nf,5) = gel(nf,6) = gen_0;
  return gerepileupto(av,padicff2(nf,p,pr));
}

static GEN
padic_trivfact(GEN x, GEN p, long r)
{
  return mkmat2(mkcol(ZX_to_ZpX_normalized(x, p, powiu(p,r), r)),
                mkcol(gen_1));
}

GEN
factorpadic2(GEN f, GEN p, long prec)
{
  pari_sp av = avma;
  GEN fa,ex,y;
  long n,i,l;

  if (typ(f)!=t_POL) pari_err(notpoler,"factorpadic");
  if (typ(p)!=t_INT) pari_err(arither1);
  if (gcmp0(f)) pari_err(zeropoler,"factorpadic");
  if (prec <= 0) pari_err(talker,"non-positive precision in factorpadic");

  n = degpol(f);
  if (n==0) return trivfact();
  f = QpX_to_ZX(f);
  if (n==1) return gerepilecopy(av, padic_trivfact(f,p,prec));
  if (!gcmp1(leading_term(f)))
    pari_err(impl,"factorpadic2 for non-monic polynomial");

  fa = ZX_squff(f, &ex);
  l = lg(fa); n = 0;
  for (i=1; i<l; i++)
  {
    gel(fa,i) = padicff(gel(fa,i),p,prec);
    n += lg(fa[i])-1;
  }

  y = fact_from_DDF(fa,ex,n);
  return gerepileupto(av, sort_factor(y, cmp_padic));
}

/***********************/
/*   Using ROUND 4     */
/***********************/
static int
expo_is_squarefree(GEN e)
{
  long i, l = lg(e);
  for (i=1; i<l; i++)
    if (e[i] != 1) return 0;
  return 1;
}

/* assume f a ZX with leading_term 1, degree > 0 */
GEN
ZX_monic_factorpadic(GEN f, GEN p, long prec)
{
  GEN w, poly, p1, p2, ex, P, E;
  long n=degpol(f), i, k, j, pr;

  if (n==1) return mkmat2(mkcol(f), mkcol(gen_1));
  pr = prec;

  poly = ZX_squff(f,&ex);
  P = cgetg(n+1,t_COL);
  E = cgetg(n+1,t_COL); n = lg(poly);
  for (j=i=1; i<n; i++)
  {
    pari_sp av1 = avma;
    GEN fx = gel(poly,i), fa = FpX_factor(fx,p);
    w = gel(fa,1);
    if (expo_is_squarefree(gel(fa,2)))
    { /* no repeated factors: Hensel lift */
      p1 = hensel_lift_fact(fx, w, NULL, p, powiu(p,pr), pr);
      p2 = utoipos(ex[i]);
      for (k=1; k<lg(p1); k++,j++)
      {
        P[j] = p1[k];
        gel(E,j) = p2;
      }
      continue;
    }
    /* use Round 4 */
    p2 = maxord_i(p, fx, Z_pval(ZX_disc(fx),p), w, pr);
    if (p2)
    {
      p2 = gerepilecopy(av1,p2);
      p1 = gel(p2,1);
      p2 = gel(p2,2);
      for (k=1; k<lg(p1); k++,j++)
      {
        P[j] = p1[k];
        gel(E,j) = mulis(gel(p2,k),ex[i]);
      }
    }
    else
    {
      avma = av1;
      gel(P,j) = fx;
      gel(E,j) = utoipos(ex[i]); j++;
    }
  }
  setlg(P,j);
  setlg(E,j); return mkmat2(P, E);
}

GEN
factorpadic4(GEN f,GEN p,long prec)
{
  pari_sp av = avma;
  GEN y, P, ppow, lead, lt;
  long i, l, pr, n = degpol(f);
  int reverse = 0;

  if (typ(f)!=t_POL) pari_err(notpoler,"factorpadic");
  if (typ(p)!=t_INT) pari_err(arither1);
  if (gcmp0(f)) pari_err(zeropoler,"factorpadic");
  if (prec <= 0) pari_err(talker,"non-positive precision in factorpadic");
  if (n == 0) return trivfact();

  f = QpX_to_ZX(f); (void)Z_pvalrem(leading_term(f), p, &lt);
  f = pnormalize(f, p, prec, n-1, &lead, &pr, &reverse);
  y = ZX_monic_factorpadic(f, p, pr);
  P = gel(y,1); l = lg(P);
  if (lead)
    for (i=1; i<l; i++) gel(P,i) = primpart( RgX_unscale(gel(P,i), lead) );
  ppow = powiu(p,prec);
  for (i=1; i<l; i++)
  {
    GEN t = gel(P,i);
    if (reverse) t = normalizepol(polrecip_i(t));
    gel(P,i) = ZX_to_ZpX_normalized(t,p,ppow,prec);
  }
  if (!gcmp1(lt)) gel(P,1) = gmul(gel(P,1), lt);
  return gerepilecopy(av, sort_factor(y, cmp_padic));
}

GEN
factorpadic0(GEN f,GEN p,long r,long flag)
{
  switch(flag)
  {
     case 0: return factorpadic4(f,p,r);
     case 1: return factorpadic2(f,p,r);
     default: pari_err(flagerr,"factorpadic");
  }
  return NULL; /* not reached */
}

/*******************************************************************/
/*                                                                 */
/*                     FACTORIZATION IN F_q                        */
/*                                                                 */
/*******************************************************************/
static GEN spec_FqXQ_pow(GEN x, GEN S, GEN T, GEN p);

static GEN
to_Fq(GEN x, GEN T, GEN p)
{
  long i, lx, tx = typ(x);
  GEN y;

  if (tx == t_INT)
    y = mkintmod(x,p);
  else
  {
    if (tx != t_POL) pari_err(typeer,"to_Fq");
    lx = lg(x);
    y = cgetg(lx,t_POL); y[1] = x[1];
    for (i=2; i<lx; i++) gel(y,i) = mkintmod(gel(x,i), p);
  }
  return mkpolmod(y, T);
}

static GEN
to_Fq_pol(GEN x, GEN T, GEN p)
{
  long i, lx, tx = typ(x);
  if (tx != t_POL) pari_err(typeer,"to_Fq_pol");
  lx = lg(x);
  for (i=2; i<lx; i++) gel(x,i) = to_Fq(gel(x,i),T,p);
  return x;
}

static GEN
to_Fq_fact(GEN P, GEN E, GEN T, GEN p, pari_sp av)
{
  GEN y = cgetg(3,t_MAT), u, v;
  long j, l = lg(P), nbf = lg(P);

  u = cgetg(nbf,t_COL); gel(y,1) = u;
  v = cgetg(nbf,t_COL); gel(y,2) = v;
  for (j=1; j<l; j++)
  {
    gel(u,j) = simplify_i(gel(P,j)); /* may contain pols of degree 0 */
    gel(v,j) = utoi((ulong)E[j]);
  }
  y = gerepilecopy(av, y); u = gel(y,1);
  p = icopy(p);
  T = FpX_to_mod(T, p);
  for (j=1; j<nbf; j++) gel(u,j) = to_Fq_pol(gel(u,j), T,p);
  return y;
}

/* split into r factors of degree d */
static void
FqX_split(GEN *t, long d, GEN q, GEN S, GEN T, GEN p)
{
  long l, v, is2, cnt, dt = degpol(*t), dT = degpol(T);
  pari_sp av;
  GEN w,w0;

  if (dt == d) return;
  v = varn(*t);
  if (DEBUGLEVEL > 6) (void)timer2();
  av = avma; is2 = equaliu(p, 2);
  for(cnt = 1;;cnt++, avma = av)
  { /* splits *t with probability ~ 1 - 2^(1-r) */
    w = w0 = FqX_rand(dt,v, T,p);
    if (degpol(w) <= 0) continue;
    for (l=1; l<d; l++) /* sum_{0<i<d} w^(q^i), result in (F_q)^r */
      w = gadd(w0, spec_FqXQ_pow(w, S, T, p));
    w = FpXQX_red(w, T,p);
    if (is2)
    {
      w0 = w;
      for (l=1; l<dT; l++) /* sum_{0<i<k} w^(2^i), result in (F_2)^r */
      {
        w = FqX_rem(FqX_sqr(w,T,p), *t, T,p);
        w = FpXX_red(gadd(w0,w), p);
      }
    }
    else
    {
      w = FpXQYQ_pow(w, shifti(q,-1), *t, T, p);
      /* w in {-1,0,1}^r */
      if (degpol(w) <= 0) continue;
      gel(w,2) = gadd(gel(w,2), gen_1);
    }
    w = FqX_gcd(*t,w, T,p); l = degpol(w);
    if (l && l != dt) break;
  }
  w = gerepileupto(av,w);
  if (DEBUGLEVEL > 6)
    fprintferr("[FqX_split] splitting time: %ld (%ld trials)\n",timer2(),cnt);
  l /= d; t[l] = FqX_div(*t,w, T,p); *t = w;
  FqX_split(t+l,d,q,S,T,p);
  FqX_split(t  ,d,q,S,T,p);
}

/* to "compare" (real) scalars and t_INTMODs */
static int
cmp_coeff(GEN x, GEN y)
{
  if (typ(x) == t_INTMOD) x = gel(x,2);
  if (typ(y) == t_INTMOD) y = gel(y,2);
  return gcmp(x,y);
}

int
cmp_pol(GEN x, GEN y)
{
  long fx[3], fy[3];
  long i,lx,ly;
  int fl;
  if (typ(x) == t_POLMOD) x = gel(x,2);
  if (typ(y) == t_POLMOD) y = gel(y,2);
  if (typ(x) == t_POL) lx = lg(x); else { lx = 3; gel(fx,2) = x; x = fx; }
  if (typ(y) == t_POL) ly = lg(y); else { ly = 3; gel(fy,2) = y; y = fy; }
  if (lx > ly) return  1;
  if (lx < ly) return -1;
  for (i=lx-1; i>1; i--)
    if ((fl = cmp_coeff(gel(x,i), gel(y,i)))) return fl;
  return 0;
}

/*******************************************************************/
/*                                                                 */
/*                  FACTOR USING TRAGER'S TRICK                    */
/*                                                                 */
/*******************************************************************/
/* Factor polynomial a on the number field defined by polynomial T */
GEN
polfnf(GEN a, GEN T)
{
  pari_sp av = avma;
  GEN x0, P, E, u, G, fa, n, unt, dent, A;
  long lx, i, k, e;
  int sqfree, tmonic;

  if (typ(a)!=t_POL || typ(T)!=t_POL) pari_err(typeer,"polfnf");
  if (gcmp0(a)) return gcopy(a);
  A = lift(fix_relative_pol(T,a,0));
  P = content(A); if (!gcmp1(P)) A = gdiv(A, P);
  T = primpart(T);
  tmonic = is_pm1(leading_term(T));

  dent = tmonic? indexpartial(T, NULL): ZX_disc(T);
  unt = mkpolmod(gen_1,T);
  G = nfgcd(A,derivpol(A), T, dent);
  sqfree = gcmp1(G);
  u = sqfree? A: Q_primpart(RgXQX_div(A, G, T));
  k = 0; n = ZY_ZXY_rnfequation(T, u, &k);
  if (DEBUGLEVEL>4) fprintferr("polfnf: choosing k = %ld\n",k);
  if (!sqfree)
  {
    G = poleval(G, gadd(pol_x[varn(A)], gmulsg(k, pol_x[varn(T)])));
    G = ZY_ZXY_resultant(T, G);
  }
  /* n guaranteed to be squarefree */
  fa = ZX_DDF(n,0); lx = lg(fa);
  P = cgetg(lx,t_COL);
  E = cgetg(lx,t_COL);
  if (lx == 2)
  { /* P^k, k irreducible */
    gel(P,1) = gmul(unt,u);
    gel(E,1) = utoipos(degpol(A) / degpol(u));
    return gerepilecopy(av, mkmat2(P,E));
  }
  x0 = gadd(pol_x[varn(A)], gmulsg(-k, mkpolmod(pol_x[varn(T)], T)));
  for (i=lx-1; i>0; i--)
  {
    GEN f = gel(fa,i), F = lift_intern(poleval(f, x0));
    if (!tmonic) F = primpart(F);
    F = nfgcd(u, F, T, dent);
    if (typ(F) != t_POL || degpol(F) == 0)
      pari_err(talker,"reducible modulus in factornf");
    e = 1;
    if (!sqfree)
    {
      while (poldvd(G,f, &G)) e++;
      if (degpol(G) == 0) sqfree = 1;
    }
    gel(P,i) = gdiv(gmul(unt,F), leading_term(F));
    gel(E,i) = utoipos(e);
  }
  return gerepilecopy(av, sort_factor(mkmat2(P,E), cmp_pol));
}

static GEN
FqX_frob_deflate(GEN f, GEN T, GEN p)
{
  GEN F = poldeflate_i(f, itos(p)), frobinv = powiu(p, degpol(T)-1);
  long i, l = lg(F);
  for (i=2; i<l; i++) gel(F,i) = Fq_pow(gel(F,i), frobinv, T,p);
  return F;
}
/* Factor _sqfree_ polynomial a on the finite field Fp[X]/(T) */
static GEN
FqX_split_Trager(GEN A, GEN T, GEN p)
{
  GEN x0, P, u, fa, n;
  long lx, i, k;

  u = A;
  n = NULL;
  for (k = 0; cmpui(k, p) < 0; k++)
  {
    GEN U = poleval(u, gadd(pol_x[varn(A)], gmulsg(k, pol_x[varn(T)])));
    n = FpY_FpXY_resultant(T, U, p);
    if (FpX_is_squarefree(n, p)) break;
    n = NULL;
  }
  if (!n) return NULL;
  if (DEBUGLEVEL>4) fprintferr("FqX_split_Trager: choosing k = %ld\n",k);
  /* n guaranteed to be squarefree */
  fa = FpX_factor(n, p); fa = gel(fa,1); lx = lg(fa);
  P = cgetg(lx,t_COL);
  if (lx == 2)
  { /* P^k, k irreducible */
    gel(P,1) = u;
    return P;
  }
  x0 = gadd(pol_x[varn(A)], gmulsg(-k, mkpolmod(pol_x[varn(T)], T)));
  for (i=lx-1; i>1; i--)
  {
    GEN f = gel(fa,i), F = lift_intern(poleval(f, x0));
    F = FqX_gcd(u, F, T, p);
    if (typ(F) != t_POL || degpol(F) == 0)
      pari_err(talker,"reducible modulus in factornf");
    u = FqX_div(u, F, T, p);
    gel(P,i) = F;
  }
  gel(P,1) = u; return P;
}


/* assume n = deg(u) > 1, X over FqX */
/* return S = [ X^q, X^2q, ... X^(n-1)q ] mod u (in Fq[X]) in Kronecker form */
static GEN
init_spec_FqXQ_pow(GEN X, GEN q, GEN u, GEN T, GEN p)
{
  long i, n = degpol(u);
  GEN x, S = cgetg(n, t_VEC);

  if (n == 1) return S;
  x = FpXQYQ_pow(X, q, u, T, p);
  gel(S,1) = x;
  if ((degpol(x)<<1) < degpol(T)) {
    for (i=2; i < n; i++)
      gel(S,i) = FqX_rem(gmul(gel(S,i-1), x), u, T,p);
  } else {
    for (i=2; i < n; i++)
      gel(S,i) = (i&1)? FqX_rem(gmul(gel(S,i-1), x), u, T,p)
                      : FqX_rem(gsqr(gel(S,i>>1)), u, T,p);
  }
  for (i=1; i < n; i++) gel(S,i) = to_Kronecker(gel(S,i), T);
  return S;
}

/* compute x^q, S is as above */
static GEN
spec_FqXQ_pow(GEN x, GEN S, GEN T, GEN p)
{
  pari_sp av = avma, lim = stack_lim(av, 1);
  GEN x0 = x+2, z = gel(x0,0);
  long i, dx = degpol(x);

  for (i = 1; i <= dx; i++)
  {
    GEN d, c = gel(x0,i);
    if (gcmp0(c)) continue;
    d = gel(S,i); if (!gcmp1(c)) d = gmul(c,d);
    z = gadd(z, d);
    if (low_stack(lim, stack_lim(av,1)))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"spec_FqXQ_pow");
      z = gerepileupto(av, z);
    }
  }
  z = FpXQX_from_Kronecker(z, T, p);
  setvarn(z, varn(x)); return gerepileupto(av, z);
}

static long
isabsolutepol(GEN f)
{
  long i, l = lg(f);
  for(i=2; i<l; i++)
  {
    GEN c = gel(f,i);
    if (typ(c) == t_POL && degpol(c) > 0) return 0;
  }
  return 1;
}

typedef struct {
  GEN S, L, Xq;
  GEN q;    /* p^deg(T) */
  GEN p, T; /* split mod (p, T(X)) */
} FqX_split_t;

static void
add(GEN z, GEN g, long d) { appendL(z, mkvec2(utoipos(d), g)); }
/* return number of roots */
long
FqX_split_deg1(GEN *pz, GEN u, GEN q, GEN T, GEN p)
{
  long dg, N = degpol(u);
  GEN v, S, g, X, z = cget1(N+1, t_VEC);

  *pz = z;
  if (N == 1) return 1;
  v = X = pol_x[varn(u)];
  S = init_spec_FqXQ_pow(X, q, u, T, p);
  appendL(z, S);
  v = spec_FqXQ_pow(v, S, T, p);
  g = FqX_gcd(gsub(v,X),u, T,p);
  dg = degpol(g);
  if (dg > 0) add(z, g, dg);
  return dg;
}

/* return number of factors */
long
FqX_split_by_degree(GEN *pz, GEN u, GEN q, GEN T, GEN p)
{
  long nb = 0, d, dg, N = degpol(u);
  GEN v, S, g, X, z = cget1(N+1, t_VEC);

  *pz = z;
  if (N == 1) return 1;
  v = X = pol_x[varn(u)];
  S = init_spec_FqXQ_pow(X, q, u, T, p);
  appendL(z, S);
  for (d=1; d <= N>>1; d++)
  {
    v = spec_FqXQ_pow(v, S, T, p);
    g = FqX_gcd(gsub(v,X),u, T,p);
    dg = degpol(g); if (dg <= 0) continue;
    /* all factors of g have degree d */
    add(z, g, dg / d); nb += dg / d;
    N -= dg;
    if (N)
    {
      u = FqX_div(u,g, T,p);
      v = FqX_rem(v,u, T,p);
    }
  }
  if (N) { add(z, u, 1); nb++; }
  return nb;
}

static GEN
FqX_split_equal(GEN L, GEN S, GEN T, GEN p)
{
  long n = itos(gel(L,1));
  GEN u = gel(L,2), z = cgetg(n + 1, t_VEC);
  gel(z,1) = u;
  FqX_split((GEN*)(z+1), degpol(u) / n, powiu(p, degpol(T)), S, T, p);
  return z;
}
GEN
FqX_split_roots(GEN z, GEN T, GEN p, GEN pol)
{
  GEN S = gel(z,1), L = gel(z,2), rep = FqX_split_equal(L, S, T, p);
  if (pol) rep = shallowconcat(rep, FqX_div(pol, gel(L,2), T,p));
  return rep;
}
GEN
FqX_split_all(GEN z, GEN T, GEN p)
{
  GEN S = gel(z,1), rep = cgetg(1, t_VEC);
  long i, l = lg(z);
  for (i = 2; i < l; i++)
    rep = shallowconcat(rep, FqX_split_equal(gel(z,i), S, T, p));
  return rep;
}

static long
FqX_sqf_split(GEN *t0, GEN q, GEN T, GEN p)
{
  GEN *t = t0, u = *t, v, S, g, X;
  long d, dg, N = degpol(u);

  if (N == 1) return 1;
  v = X = pol_x[varn(u)];
  S = init_spec_FqXQ_pow(X, q, u, T, p);
  for (d=1; d <= N>>1; d++)
  {
    v = spec_FqXQ_pow(v, S, T, p);
    g = FqX_gcd(gsub(v,X),u, T,p);
    dg = degpol(g); if (dg <= 0) continue;

    /* all factors of g have degree d */
    *t = g;
    FqX_split(t, d, q, S, T, p);
    t += dg / d;
    N -= dg;
    if (N)
    {
      u = FqX_div(u,g, T,p);
      v = FqX_rem(v,u, T,p);
    }
  }
  if (N) *t++ = u;
  return t - t0;
}

/* not memory-clean */
/* TODO: provide a public and clean FpX_factorff */
static GEN
FpX_factorff(GEN P,GEN l, GEN Q)
{
  GEN V,E, F = FpX_factor(P,l);
  long lfact = 1, nmax = lgpol(P), n = lg(gel(F,1));
  long i;
  V = cgetg(nmax,t_VEC);
  E = cgetg(nmax,t_VECSMALL);
  for(i=1;i<n;i++)
  {
    GEN R = FpX_factorff_irred(gmael(F,1,i),Q,l);
    long j, r = lg(R);
    for (j=1;j<r;j++)
    {
      V[lfact] = R[j];
      E[lfact] = mael(F,2,i); lfact++;
    }
  }
  setlg(V,lfact);
  setlg(E,lfact); return sort_factor(mkvec2(V,E), cmp_pol);
}

static GEN
FqX_factor_i(GEN f, GEN T, GEN p)
{
  long pg, j, k, d, e, N, nbfact, pk;
  GEN E, f2, f3, df1, df2, g1, u, q, *t;

  if (!signe(f)) pari_err(zeropoler,"FqX_factor");
  d = degpol(f); if (!d) return trivfact();
  T = FpX_normalize(T, p);
  f = FqX_normalize(f, T, p);
  if (isabsolutepol(f)) return FpX_factorff(simplify_i(f), p, T);

  pg = itos_or_0(p);
  df2  = NULL; /* gcc -Wall */
  t = (GEN*)cgetg(d+1,t_VEC);
  E = cgetg(d+1, t_VECSMALL);

  q = powiu(p, degpol(T));
  e = nbfact = 1;
  pk = 1;
  f3 = NULL;
  df1 = FqX_deriv(f, T, p);
  for(;;)
  {
    long nb0;
    while (gcmp0(df1))
    { /* needs d >= p: pg = 0 can't happen  */
      pk *= pg; e = pk;
      f = FqX_frob_deflate(f, T, p);
      df1 = FqX_deriv(f, T, p); f3 = NULL;
    }
    f2 = f3? f3: FqX_gcd(f,df1, T,p);
    if (!degpol(f2)) u = f;
    else
    {
      g1 = FqX_div(f,f2, T,p);
      df2 = FqX_deriv(f2, T,p);
      if (gcmp0(df2)) { u = g1; f3 = f2; }
      else
      {
        f3 = FqX_gcd(f2,df2, T,p);
        u = degpol(f3)? FqX_div(f2, f3, T,p): f2;
        u = FqX_div(g1, u, T,p);
      }
    }
    /* u is square-free (product of irreducibles of multiplicity e) */
    nb0 = nbfact; N = degpol(u);
    t[nbfact] = FqX_normalize(u, T,p);
    if (N) {
      nb0 = nbfact;
      t[nbfact] = FqX_normalize(u, T,p);
      if (N == 1) nbfact++;
      else
      {
#if 0
        nbfact += FqX_split_Berlekamp(t+nbfact, q, T, p);
#else
        GEN P = FqX_split_Trager(t[nbfact], T, p);
        if (P) {
          for (j = 1; j < lg(P); j++) t[nbfact++] = gel(P,j);
        } else {
          if (DEBUGLEVEL) pari_warn(warner, "FqX_split_Trager failed!");
          nbfact += FqX_sqf_split(t+nbfact, q, T, p);
        }
#endif
      }
      for (j = nb0; j < nbfact; j++) E[j] = e;
    }
    if (!degpol(f2)) break;
    f = f2; df1 = df2; e += pk;
  }

  for (j=1; j<nbfact; j++)
  {
    t[j] = FqX_normalize(gel(t,j), T,p);
    for (k=1; k<j; k++)
      if (gequal(t[j],t[k]))
      {
        E[k] += E[j]; nbfact--;
        E[j] = E[nbfact];
        t[j] = t[nbfact]; break;
      }
  }
  setlg(t, nbfact);
  setlg(E, nbfact); return sort_factor(mkvec2((GEN)t, E), cmp_pol);
}
GEN
factorff(GEN f, GEN p, GEN T)
{
  pari_sp av;
  long v;
  GEN z;

  if (typ(T)!=t_POL || typ(f)!=t_POL || typ(p)!=t_INT) pari_err(typeer,"factorff");
  v = varn(T);
  if (varncmp(v, varn(f)) <= 0)
    pari_err(talker,"polynomial variable must have higher priority in factorff");
  T = RgX_to_FpX(T, p); av = avma;
  z = FqX_factor_i(RgX_to_FqX(f, T,p), T, p);
  return to_Fq_fact(gel(z,1),gel(z,2), T,p,av);
}
/* factorization of x modulo (T,p). Assume x already reduced */
GEN
FqX_factor(GEN x, GEN T, GEN p)
{
  pari_sp av = avma;
  if (!T) return FpX_factor(x, p);
  return gerepilecopy(av, FqX_factor_i(x, T, p));
}
/* See also: Isomorphisms between finite field and relative
 * factorization in polarit3.c */

/*******************************************************************/
/*                                                                 */
/*                       COMPLEX ROOTS                             */
/*                                                                 */
/*******************************************************************/
static GEN laguer(GEN pol,long N,GEN y0,long EPS,long PREC);
GEN zrhqr(GEN a,long PREC);

GEN
rootsold(GEN x, long prec)
{
  long i, j, f, real, exact, fr, deg, ln;
  pari_sp av=avma, av0, av1, av2, av3;
  long exc,expmin,m,deg0,k,ti,h,ii,e;
  GEN y,xc,xd0,xd,xdabs,p1,p2,p3,p4,p5,p6,p7;
  GEN p11,p12,p1r,p1i,pa,pax,pb,pp,pq,ps, pi;

  if (typ(x)!=t_POL) pari_err(typeer,"rootsold");
  deg0 = degpol(x); expmin = 12 - bit_accuracy(prec);
  if (!signe(x)) pari_err(zeropoler,"rootsold");
  y = cgetg(deg0+1,t_COL); if (!deg0) return y;
  for (i=1; i<=deg0; i++)
  {
    p1 = cgetc(prec); gel(y,i) = p1;
    for (j=3; j<prec; j++) (gel(p1,2))[j] = (gel(p1,1))[j] = 0;
  }
  real=1; exact=1;
  for (i=2; i<=deg0+2; i++)
  {
    ti = typ(x[i]);
    if (ti==t_REAL) exact = 0;
    else if (ti==t_QUAD)
    {
      p2 = gmael3(x,i,1,2);
      if (gsigne(p2) > 0) real = 0;
    } else if (ti != t_INT && ti != t_FRAC) real = 0;
  }
  av1 = avma;
  k = polvaluation_inexact(x, &pax);
  for (i = 1; i <= k; i++) gaffsg(0,gel(y,i));
  if (k == deg0) return y;

  pi = mppi(DEFAULTPREC);
  p2 = mkcomplex(pi, divrs(pi,10)); /* Pi * (1+I/10) */
  p11 = cgetg(4,t_POL); p11[1] = x[1];
  gel(p11,3) = gen_1;

  p12 = cgetg(5,t_POL); p12[1] = x[1];
  gel(p12,4) = gen_1;

  xd0 = derivpol(pax); pa = pax;
  pq = NULL; /* for lint */
  if (exact) { pp = ggcd(pax,xd0); h = degpol(pp); if (h) pq = RgX_div(pax,pp); }
  else{ pp = gen_1; h = 0; }
  m = 0;
  while (k != deg0)
  {
    m++;
    if (h)
    {
      pa = pp; pb = pq; pp = ggcd(pa,derivpol(pa)); h = degpol(pp);
      pq = h? RgX_div(pa,pp): pa;
      ps = RgX_div(pb,pq);
    }
    else ps = pa;
    deg = degpol(ps); if (!deg) continue;

    /* roots of exact order m */
    e = gexpo(ps) - gexpo(leading_term(ps));
    if (e < 0) e = 0; if (ps!=pax) xd0 = derivpol(ps);
    xdabs = cgetg(deg+2,t_POL); xdabs[1] = xd0[1];
    for (i=2; i<deg+2; i++)
    {
      av3 = avma; p3 = gel(xd0,i);
      gel(xdabs,i) = gerepileupto(av3, gadd(gabs(real_i(p3),prec),
                                     gabs(imag_i(p3),prec)));
    }
    av0 = avma; xc = gcopy(ps); xd = gcopy(xd0); av2 = avma;
    for (i=1; i<=deg; i++)
    {
      if (i == deg)
      {
        p1 = (GEN)y[k+m*i];
        gdivz(gneg_i(gel(xc,2)),gel(xc,3), p1);
        p1r = gel(p1,1);
        p1i = gel(p1,2);
      }
      else
      {
        p3 = gshift(p2,e);
        p4 = poleval(xc,p3);
        p5 = gnorm(p4);
        exc = 0;
        while (exc >= -20)
        {
          p7 = gneg_i(gdiv(p4, poleval(xd,p3)));
          av3 = avma;
          exc = gcmp0(p5)? -32: expo(gnorm(p7))-expo(gnorm(p3));
          avma = av3;
          for (j=1; j<=10; j++)
          {
            GEN p8, p9, p10;
            p8 = gadd(p3,p7);
            p9 = poleval(xc,p8);
            p10= gnorm(p9);
            if (exc < -20 || cmprr(p10,p5) < 0)
            {
              GEN *gptr[3];
              p3=p8; p4=p9; p5=p10;
              gptr[0]=&p3; gptr[1]=&p4; gptr[2]=&p5;
              gerepilemanysp(av2,av3,gptr,3);
              break;
            }
            gshiftz(p7,-2,p7); avma = av3;
          }
          if (j > 10)
          {
            if (DEBUGLEVEL)
              pari_warn(warner,"too many iterations in rootsold(): using roots2()");
            avma = av; return roots2(x,prec);
          }
        }
        p1 = (GEN)y[k+m*i];
        p1r = gel(p1,1); setlg(p1r, 3);
        p1i = gel(p1,2); setlg(p1i, 3); gaffect(p3, p1); avma = av2;
        for (ln = 4; ln <= prec; ln = (ln<<1)-2)
        {
          setlg(p1r,ln); if (gcmp0(p1r)) gel(p1,1) = gen_0;
          setlg(p1i,ln); if (gcmp0(p1i)) gel(p1,2) = gen_0;
          p6 = gadd(p1, gneg_i(gdiv(poleval(xc,p1), poleval(xd,p1))));
          gel(p1,1) = p1r;
          gel(p1,2) = p1i; gaffect(p6, p1); avma = av2;
        }
      }
      setlg(p1r,prec);
      setlg(p1i,prec); p7 = gcopy(p1);
      p1r = gel(p7,1); setlg(p1r,prec+1);
      p1i = gel(p7,2); setlg(p1i,prec+1);
      for (ii=1; ii<=5; ii++)
      {
        if (typ(p7) == t_COMPLEX)
        {
          if (gcmp0(gel(p7,1))) gel(p7,1) = gen_0;
          if (gcmp0(gel(p7,2))) gel(p7,2) = gen_0;
        }
        p7 = gadd(p7, gneg_i(gdiv(poleval(ps,p7), poleval(xd0,p7))));
      }
      gaffect(p7, p1);
      p6 = gdiv(poleval(ps,p7), poleval(xdabs,gabs(p7,prec)));
      if (gexpo(p6) >= expmin)
      {
        if (DEBUGLEVEL) pari_warn(warner,"error in rootsold(): using roots2()");
        avma = av; return roots2(x,prec);
      }
      avma = av2;
      if (expo(p1[2]) < expmin && real)
      {
        gaffect(gen_0, gel(p1,2));
        for (j=1; j<m; j++) gaffect(p1, (GEN)y[k+(i-1)*m+j]);
        gel(p11,2) = gneg(gel(p1,1));
        xc = gerepileupto(av0, RgX_div(xc,p11));
      }
      else
      {
        for (j=1; j<m; j++) gaffect(p1, (GEN)y[k+(i-1)*m+j]);
        if (real)
        {
          p1 = gconj(p1);
          for (j=1; j<=m; j++) gaffect(p1, (GEN)y[k+i*m+j]);
          i++;
          gel(p12,2) = gnorm(p1);
          gel(p12,3) = gmulsg(-2,gel(p1,1));
          xc = gerepileupto(av0, RgX_div(xc,p12));
        }
        else
        {
          gel(p11,2) = gneg(p1);
          xc = gerepileupto(av0, RgX_div(xc,p11));
        }
      }
      xd = derivpol(xc); av2 = avma;
    }
    k += deg*m;
  }
  avma = av1;
  for (j=2; j<=deg0; j++)
  {
    p1 = gel(y,j);
    fr = !gcmp0(gel(p1,2));
    for (k=j-1; k>=1; k--)
    {
      p2 = gel(y,k);
      f = !gcmp0(gel(p2,2));
      if (!f && fr) break;
      if (f == fr && gcmp(gel(p2,1),gel(p1,1)) <= 0) break;
      y[k+1] = y[k];
    }
    gel(y,k+1) = p1;
  }
  return y;
}

GEN
roots2(GEN pol,long PREC)
{
  pari_sp av = avma;
  long N,flagexactpol,flagrealpol,flagrealrac,ti,i,j;
  long nbpol, k, multiqol, deg, nbroot, fr, f, EPS;
  pari_sp av1;
  GEN p1,p2,rr,qol,qolbis,x,b,c,ad,v, ex, factors;

  if (typ(pol)!=t_POL) pari_err(typeer,"roots2");
  if (!signe(pol)) pari_err(zeropoler,"roots2");
  N=degpol(pol);
  if (!N) return cgetg(1,t_COL);
  if (N==1)
  {
    p1 = gmul(real_1(PREC),gel(pol,3));
    p2 = gneg_i(gdiv(gel(pol,2),p1));
    return gerepilecopy(av,p2);
  }
  EPS = 12 - bit_accuracy(PREC); /* 2^EPS is "zero" */
  flagrealpol = flagexactpol = 1;
  for (i=2; i<=N+2; i++)
  {
    c = gel(pol,i);
    switch (typ(c)) {
      case t_INT: case t_FRAC: break;

      case t_REAL: flagexactpol = 0; break;

      case t_COMPLEX: flagexactpol = flagrealpol = 0; break;

      case t_QUAD: flagexactpol = 0;
        if (gsigne(gmael(c,1,2)) > 0) flagrealpol = 0;
        break;
      default: pari_err(typeer, "roots2");
    }
  }
  rr=cgetg(N+1,t_COL);
  for (i=1; i<=N; i++)
  {
    p1 = cgetc(PREC); gel(rr,i) = p1;
    for (j=3; j<PREC; j++) mael(p1,2,j) = mael(p1,1,j) = 0;
  }
  if (flagexactpol) { pol = Q_primpart(pol); factors = ZX_squff(pol, &ex); }
  else
  {
    factors = mkcol(pol);
    ex = mkvecsmall(1);
  }
  nbpol = lg(ex)-1;
  nbroot= 0;
  for (k=1; k<=nbpol; k++)
  {
    av1=avma; qol = gel(factors,k); qolbis = gcopy(qol);
    multiqol = ex[k]; deg = degpol(qol);
    for (j=deg; j>=1; j--)
    {
      x = gen_0; flagrealrac = 0;
      if (j==1) x = gneg_i(gdiv(gel(qolbis,2),gel(qolbis,3)));
      else
      {
        x = laguer(qolbis,j,x,EPS,PREC);
        if (x == NULL) goto RLAB;
      }
      if (flagexactpol)
      {
        x = gprec_w(x, PREC+2);
        x = laguer(qol,deg,x, EPS-BITS_IN_LONG, PREC+1);
      }
      else
        x = laguer(qol,deg,x,EPS,PREC);
      if (x == NULL) goto RLAB;

      if (typ(x)==t_COMPLEX && gexpo(imag_i(x)) <= gexpo(real_i(x)) + EPS+1)
        { gel(x,2) = gen_0; flagrealrac = 1; }
      else if (j==1 && flagrealpol)
        { gel(x,2) = gen_0; flagrealrac = 1; }
      else if (typ(x)!=t_COMPLEX) flagrealrac = 1;

      for (i=1; i<=multiqol; i++) gaffect(x, gel(rr,nbroot+i));
      nbroot += multiqol;
      if (!flagrealpol || flagrealrac)
      {
        ad = new_chunk(j+1);
        for (i=0; i<=j; i++) ad[i] = qolbis[i+2];
        b = gel(ad,j);
        for (i=j-1; i>=0; i--)
        {
          c = gel(ad,i); gel(ad,i) = b;
          b = gadd(gmul(gel(rr,nbroot),b),c);
        }
        v = cgetg(j+1,t_VEC); for (i=1; i<=j; i++) v[i] = ad[j-i];
        qolbis = gtopoly(v,varn(qolbis));
        if (flagrealpol)
          for (i=2; i<=j+1; i++)
            if (typ(qolbis[i])==t_COMPLEX) gmael(qolbis,i,2)= gen_0;
      }
      else
      {
        ad = new_chunk(j-1);
        ad[j-2] = qolbis[j+2];
        p1 = gmulsg(2,real_i(gel(rr,nbroot)));
        p2 = gnorm(gel(rr,nbroot));
        gel(ad,j-3) = gadd(gel(qolbis,j+1), gmul(p1,gel(ad,j-2)));
        for (i=j-2; i>=2; i--)
          gel(ad,i-2) = gadd(gel(qolbis,i+2),
                             gsub(gmul(p1,gel(ad,i-1)),
                                  gmul(p2,gel(ad,i))));
        v = cgetg(j,t_VEC); for (i=1; i<=j-1; i++) v[i] = ad[j-1-i];
        qolbis = gtopoly(v,varn(qolbis));
        for (i=2; i<=j; i++)
          if (typ(qolbis[i])==t_COMPLEX) gmael(qolbis,i,2)= gen_0;
        for (i=1; i<=multiqol; i++)
          gaffect(gconj(gel(rr,nbroot)), gel(rr,nbroot+i));
        nbroot+=multiqol; j--;
      }
    }
    avma=av1;
  }
  for (j=2; j<=N; j++)
  {
    x=gel(rr,j); if (gcmp0(gel(x,2))) fr=0; else fr=1;
    for (i=j-1; i>=1; i--)
    {
      if (gcmp0(gmael(rr,i,2))) f=0; else f=1;
      if (f<fr) break;
      if (f==fr && gcmp(real_i(gel(rr,i)),real_i(x)) <= 0) break;
      rr[i+1]=rr[i];
    }
    gel(rr,i+1) = x;
  }
  return gerepilecopy(av,rr);

 RLAB:
  avma = av;
  for(i=2;i<=N+2;i++)
  {
    ti = typ(pol[i]);
    if (!is_intreal_t(ti)) pari_err(talker,"too many iterations in roots");
  }
  if (DEBUGLEVEL)
  {
    fprintferr("too many iterations in roots2() ( laguer() ):\n");
    fprintferr("     real coefficients polynomial, using zrhqr()\n");
  }
  return zrhqr(pol,PREC);
}

#define MR 8
#define MT 10

static GEN
laguer(GEN pol,long N,GEN y0,long EPS,long PREC)
{
  long MAXIT, iter, j;
  pari_sp av = avma, av1;
  GEN rac,erre,I,x,abx,abp,abm,dx,x1,b,d,f,g,h,sq,gp,gm,g2,*ffrac;

  MAXIT = MR*MT; rac = cgetc(PREC);
  av1 = avma;
  I = mkcomplex(gen_1,gen_1);
  ffrac = (GEN*)new_chunk(MR+1);
  ffrac[0] = dbltor(0.0);
  ffrac[1] = dbltor(0.5);
  ffrac[2] = dbltor(0.25);
  ffrac[3] = dbltor(0.75);
  ffrac[4] = dbltor(0.13);
  ffrac[5] = dbltor(0.38);
  ffrac[6] = dbltor(0.62);
  ffrac[7] = dbltor(0.88);
  ffrac[8] = dbltor(1.0);
  x=y0;
  for (iter=1; iter<=MAXIT; iter++)
  {
    b = gel(pol,N+2); d = f = gen_0;
    erre = QuickNormL1(b,PREC);
    abx  = QuickNormL1(x,PREC);
    for (j=N-1; j>=0; j--)
    {
      f = gadd(gmul(x,f), d);
      d = gadd(gmul(x,d), b);
      b = gadd(gmul(x,b), gel(pol,j+2));
      erre = gadd(QuickNormL1(b,PREC), gmul(abx,erre));
    }
    erre = gmul2n(erre, EPS);
    if (gcmp(QuickNormL1(b,PREC),erre)<=0)
    {
      gaffect(x,rac); avma = av1; return rac;
    }
    g = gdiv(d,b);
    g2 = gsqr(g); h = gsub(g2, gmul2n(gdiv(f,b),1));
    sq = gsqrt(gmulsg(N-1,gsub(gmulsg(N,h),g2)),PREC);
    gp = gadd(g,sq); abp = gnorm(gp);
    gm = gsub(g,sq); abm = gnorm(gm);
    if (gcmp(abp,abm) < 0) gp = gm;
    if (gsigne(gmax(abp,abm)) > 0)
      dx = gdivsg(N,gp);
    else
      dx = gmul(gadd(gen_1,abx), gexp(gmulgs(I,iter),PREC));
    x1 = gsub(x,dx);
    if (gexpo(QuickNormL1(gsub(x,x1),PREC)) < EPS)
    {
      gaffect(x,rac); avma = av1; return rac;
    }
    if (iter%MT) x = gcopy(x1); else x = gsub(x, gmul(ffrac[iter/MT],dx));
  }
  avma = av; return NULL;
}
#undef MR
#undef MT
/***********************************************************************/
/**                                                                   **/
/**             ROOTS of a polynomial with REAL coeffs                **/
/**                                                                   **/
/***********************************************************************/
#define RADIX 1L
#define COF 0.95

/* x t_MAT in M_n(R) : compute a symmetric matrix with the same eigenvalues */
static GEN
balanc(GEN x)
{
  pari_sp av = avma;
  long last, i, j, sqrdx = (RADIX<<1), n = lg(x);
  GEN r, c, cofgen, a = shallowcopy(x);

  last = 0; cofgen = dbltor(COF);
  while (!last)
  {
    last = 1;
    for (i=1; i<n; i++)
    {
      r = c = gen_0;
      for (j=1; j<n; j++)
        if (j!=i)
        {
          c = gadd(c, gabs(gcoeff(a,j,i),0));
          r = gadd(r, gabs(gcoeff(a,i,j),0));
        }
      if (!gcmp0(r) && !gcmp0(c))
      {
        GEN g, s = gmul(cofgen, gadd(c,r));
        long ex = 0;
        g = gmul2n(r,-RADIX); while (gcmp(c,g) < 0) {ex++; c=gmul2n(c, sqrdx);}
        g = gmul2n(r, RADIX); while (gcmp(c,g) > 0) {ex--; c=gmul2n(c,-sqrdx);}
        if (gcmp(gadd(c,r), gmul2n(s,ex)) < 0)
        {
          last = 0;
          for (j=1; j<n; j++) gcoeff(a,i,j) = gmul2n(gcoeff(a,i,j),-ex);
          for (j=1; j<n; j++) gcoeff(a,j,i) = gmul2n(gcoeff(a,j,i), ex);
        }
      }
    }
  }
  return gerepilecopy(av, a);
}

#define SIGN(a,b) ((b)>=0.0 ? fabs(a) : -fabs(a))
/* find the eigenvalues of the symmetric matrix mat */
static GEN
hqr(GEN mat)
{
  long n = lg(mat)-1, N, m, l, k, j, i, mmin, flj, flk;
  double **a, p, q, r, s, t, u, v, w, x, y, z, anorm, *wr, *wi;
  const double eps = 0.000001;
  GEN eig;

  init_dalloc();
  a = (double**)stackmalloc(sizeof(double*)*(n+1));
  for (i=1; i<=n; i++) a[i] = (double*)stackmalloc(sizeof(double)*(n+1));
  for (j=1; j<=n; j++)
    for (i=1; i<=n; i++) a[i][j] = gtodouble(gcoeff(mat,i,j));
  wr = (double*)stackmalloc(sizeof(double)*(n+1));
  wi = (double*)stackmalloc(sizeof(double)*(n+1));

  anorm = fabs(a[1][1]);
  for (i=2; i<=n; i++) for (j=(i-1); j<=n; j++) anorm += fabs(a[i][j]);
  N = n; t = 0.;
  p = q = r = 0.; /* -Wall */
  if (DEBUGLEVEL>3) { fprintferr("* Finding eigenvalues\n"); flusherr(); }
  while (N>=1)
  {
    long its = 0;
    for(;;)
    {
      for (l=N; l>=2; l--)
      {
        s = fabs(a[l-1][l-1])+fabs(a[l][l]); if (s==0.) s = anorm;
        if (fabs(a[l][l-1])+s == s) break;
      }
      x = a[N][N];
      if (l == N){ wr[N] = x+t; wi[N] = 0.; N--; break; } /* OK */

      y = a[N-1][N-1];
      w = a[N][N-1]*a[N-1][N];
      if (l == N-1)
      {
        p = 0.5*(y-x); q = p*p+w; z = sqrt(fabs(q)); x += t;
        if (q >= 0. || fabs(q) <= eps)
        {
          z = p + SIGN(z,p);
          wr[N-1] = wr[N] = x+z;
          if (fabs(z)>eps) wr[N] = x-w/z;
          wi[N-1] = wi[N] = 0.;
        }
        else { wr[N-1] = wr[N]= x+p; wi[N-1] = -z; wi[N] = z; }
        N -= 2; break; /* OK */
      }

      if (its==30) pari_err(talker,"too many iterations in hqr");
      if (its==10 || its==20)
      {
        t += x; for (i=1; i<=N; i++) a[i][i] -= x;
        s = fabs(a[N][N-1]) + fabs(a[N-1][N-2]);
        y = x = 0.75*s;
        w = -0.4375*s*s;
      }
      its++;
      for (m=N-2; m>=l; m--)
      {
        z = a[m][m]; r = x-z; s = y-z;
        p = (r*s-w)/a[m+1][m]+a[m][m+1];
        q = a[m+1][m+1]-z-r-s;
        r = a[m+2][m+1];
        s = fabs(p)+fabs(q)+fabs(r); p/=s; q/=s; r/=s;
        if (m==l) break;
        u = fabs(a[m][m-1])*(fabs(q)+fabs(r));
        v = fabs(p) * (fabs(a[m-1][m-1])+fabs(z)+fabs(a[m+1][m+1]));
        if (u+v==v) break;
      }
      for (i=m+2; i<=N; i++){ a[i][i-2]=0.; if (i!=m+2) a[i][i-3]=0.; }
      for (k=m; k<=N-1; k++)
      {
        if (k!=m)
        {
          p = a[k][k-1]; q = a[k+1][k-1];
          r = (k != N-1)? a[k+2][k-1]: 0.;
          x = fabs(p)+fabs(q)+fabs(r);
          if (x != 0.) { p/=x; q/=x; r/=x; }
        }
        s = SIGN(sqrt(p*p+q*q+r*r),p);
        if (s == 0.) continue;

        if (k==m)
          { if (l!=m) a[k][k-1] = -a[k][k-1]; }
        else
          a[k][k-1] = -s*x;
        p+=s; x=p/s; y=q/s; z=r/s; q/=p; r/=p;
        for (j=k; j<=N; j++)
        {
          p = a[k][j]+q*a[k+1][j];
          if (k != N-1) { p+=r*a[k+2][j]; a[k+2][j]-=p*z; }
          a[k+1][j] -= p*y; a[k][j] -= p*x;
        }
        mmin = (N < k+3)? N: k+3;
        for (i=l; i<=mmin; i++)
        {
          p = x*a[i][k]+y*a[i][k+1];
          if (k != N-1) { p+=z*a[i][k+2]; a[i][k+2]-=p*r; }
          a[i][k+1] -= p*q; a[i][k] -= p;
        }
      }
    }
  }
  for (j=2; j<=n; j++) /* ordering the roots */
  {
    x = wr[j];
    y = wi[j]; flj = (y != 0.);
    for (k=j-1; k>=1; k--)
    {
      flk = (wi[k] != 0.);
      if (!flk && flj) break;
      if (flk == flj && wr[k] <= x) break;
      wr[k+1] = wr[k];
      wi[k+1] = wi[k];
    }
    wr[k+1] = x;
    wi[k+1] = y;
  }
  if (DEBUGLEVEL>3) { fprintferr("* Eigenvalues computed\n"); flusherr(); }
  eig = cgetg(n+1,t_COL);
  for (i=1; i<=n; i++)
    gel(eig,i) = (wi[i] == 0.)? dbltor(wr[i])
                              : mkcomplex(dbltor(wr[i]), dbltor(wi[i]));
  return eig;
}

/* a t_POL in R[X], squarefree: give the roots of the polynomial a (real roots
 * first) in increasing order of their real parts. */
GEN
zrhqr(GEN a, long prec)
{
  pari_sp av = avma;
  long i, prec2, n = degpol(a), ex = -bit_accuracy(prec);
  GEN aa, b, rt, rr, x, dx, y, newval, oldval;

  rt = hqr(balanc(assmat(a)));
  prec2 = 2*prec; /* polishing the roots */
  aa = gprec_w(a, prec2);
  b = derivpol(aa); rr = cgetg(n+1,t_COL);
  for (i=1; i<=n; i++)
  {
    x = gprec_w(gel(rt,i), prec2);
    for (oldval=NULL;; oldval=newval, x=y)
    { /* Newton iteration */
      dx = poleval(b,x);
      if (gexpo(dx) < ex)
        pari_err(talker,"polynomial has probably multiple roots in zrhqr");
      y = gsub(x, gdiv(poleval(aa,x),dx));
      newval = gabs(poleval(aa,y),prec2);
      if (gexpo(newval) < ex || (oldval && gcmp(newval,oldval) > 0)) break;
    }
    if (DEBUGLEVEL>3) fprintferr("%ld ",i);
    gel(rr,i) = gtofp(y, prec);
  }
  if (DEBUGLEVEL>3) { fprintferr("\npolished roots = %Z",rr); flusherr(); }
  return gerepilecopy(av, rr);
}
