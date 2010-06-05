/* $Id: nffactor.c 11784 2009-06-26 09:38:07Z bill $

Copyright (C) 2000-2004  The PARI group.

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
/*            POLYNOMIAL FACTORIZATION IN A NUMBER FIELD           */
/*                                                                 */
/*******************************************************************/
#include "pari.h"
#include "paripriv.h"

static GEN nfsqff(GEN nf,GEN pol,long fl);

/* FIXME: neither nfcmbf nor LLL_cmbf can handle the non-nf case */
/* for nf_bestlift: reconstruction of algebraic integers known mod \wp^k */
typedef struct {
  long k;    /* input known mod \wp^k */
  GEN p, pk;    /* p^k */
  GEN den;   /* denom(prk^-1) = p^k [ assume pr unramified ] */
  GEN prk;   /* |.|^2 LLL-reduced basis (b_i) of \wp^k  (NOT T2-reduced) */
  GEN prkHNF;/* HNF basis of \wp^k */
  GEN iprk;  /* den * prk^-1 */
  GEN GSmin; /* min |b_i^*|^2 */

  GEN Tp; /* Tpk mod p */
  GEN Tpk;
  GEN ZqProj;/* projector to Zp / \wp^k = Z/p^k[X] / Tpk */

  GEN tozk;
  GEN topow;
  GEN topowden; /* topow x / topowden = basistoalg(x) */
} nflift_t;

typedef struct /* for use in nfsqff */
{
  nflift_t *L;
  GEN nf;
  GEN pol, polbase; 
  GEN fact; 
  GEN pr; 
  GEN Br, bound, ZC, BS_2;
  GEN dn;
  long hint;
} nfcmbf_t;

static GEN
unifpol0(GEN nf,GEN x,long flag)
{
  switch(typ(x))
  {
    case t_INT: case t_FRAC:
      return gcopy(x);

    case t_POLMOD:
      x = gel(x,2); /* fall through */
      if (typ(x) != t_POL) return gcopy(x); /* scalar */
    case t_POL:
      if (!degpol(x)) return gcopy(constant_term(x));
      return (flag == t_COL)? algtobasis(nf, x): gmodulo(x, gel(nf,1));

    default: /* t_COL */
      return (flag == t_COL)? gcopy(x): basistoalg(nf, x);
  }
}

/* Let x be a polynomial with coefficients in Z or nf (vectors or polymods)
 * return the same polynomial with coefficients expressed:
 *  if flag=t_COL: as vectors (on the integral basis).
 *  if flag=t_POLMOD: as polmods.
 */
GEN
unifpol(GEN nf, GEN x, long flag)
{
  if (typ(x)==t_POL && varncmp(varn(x), varn(nf[1])) < 0)
  {
    long i, d = lg(x);
    GEN y = cgetg(d,t_POL); y[1] = x[1];
    for (i=2; i<d; i++) gel(y,i) = unifpol0(nf, gel(x,i), flag);
    return y;
  }
  return unifpol0(nf, x, flag);
}

GEN
nffactormod(GEN nf, GEN x, GEN pr)
{
  long j, l, vx = varn(x), vn;
  pari_sp av = avma;
  GEN F, E, rep, xrd, modpr, T, p;

  nf = checknf(nf);
  vn = varn(nf[1]);
  if (typ(x)!=t_POL) pari_err(typeer,"nffactormod");
  if (varncmp(vx,vn) >= 0)
    pari_err(talker,"polynomial variable must have highest priority in nffactormod");

  modpr = nf_to_ff_init(nf, &pr, &T, &p);
  xrd = modprX(x, nf, modpr);
  rep = FqX_factor(xrd,T,p);
  settyp(rep, t_MAT);
  F = gel(rep,1); l = lg(F);
  E = gel(rep,2); settyp(E, t_COL);
  for (j = 1; j < l; j++) {
    gel(F,j) = modprX_lift(gel(F,j), modpr);
    gel(E,j) = stoi(E[j]);
  }
  return gerepilecopy(av, rep);
}

static GEN
QXQX_normalize(GEN P, GEN T)
{
  GEN P0 = leading_term(P);
  if (!gcmp1(P0))
  {
    long t = typ(P0);
    if (t == t_POL && !degpol(P0)) P0 = gel(P0,2);
    if (is_rational_t(t))
      P = gdiv(P, P0);
    else
      P = RgXQX_RgXQ_mul(P, QXQ_inv(P0,T), T);
  }
  return P;
}

/* to INT / FRAC / (POLMOD mod T), not memory clean because T not copied */
static GEN
RgXQ_to_mod(GEN x, GEN T)
{
  long d;
  switch(typ(x))
  {
    case t_INT: case t_FRAC:
      return gcopy(x);

    default:
      d = degpol(x);
      if (d < 0) return gen_0;
      if (d == 0) return gcopy(gel(x,2));
      return mkpolmod(gcopy(x), T);
  }
}
/* T a ZX, z lifted from (Q[Y]/(T(Y)))[X], apply RgXQ_to_mod to all coeffs.
 * Not memory clean because T not copied */
static GEN
RgXQX_to_mod(GEN z, GEN T)
{
  long i,l = lg(z);
  GEN x = cgetg(l,t_POL);
  for (i=2; i<l; i++) gel(x,i) = RgXQ_to_mod(gel(z,i), T);
  x[1] = z[1]; return normalizepol_i(x,l);
}
/* Apply RgXQX_to_mod to all entries. Memory-clean ! */
static GEN
RgXQXV_to_mod(GEN V, GEN T)
{
  long i, l = lg(V);
  GEN z = cgetg(l, t_VEC); T = gcopy(T);
  for (i=1;i<l; i++) gel(z,i) = RgXQX_to_mod(gel(V,i), T);
  return z;
}
/* Apply RgXQ_to_mod to all entries. Memory-clean ! */
static GEN
RgXQV_to_mod(GEN V, GEN T)
{
  long i, l = lg(V);
  GEN z = cgetg(l, t_VEC); T = gcopy(T);
  for (i=1;i<l; i++) gel(z,i) = RgXQ_to_mod(gel(V,i), T);
  return z;
}

/* return the roots of pol in nf */
GEN
nfroots(GEN nf,GEN pol)
{
  pari_sp av = avma;
  GEN A,g, T;
  long d;

  if (!nf) return nfrootsQ(pol);

  nf = checknf(nf); T = gel(nf,1);
  if (typ(pol) != t_POL) pari_err(notpoler,"nfroots");
  if (varncmp(varn(pol), varn(T)) >= 0)
    pari_err(talker,"polynomial variable must have highest priority in nfroots");
  d = degpol(pol);
  if (d == 0) return cgetg(1,t_VEC);
  if (d == 1)
  {
    A = gneg_i(gdiv(gel(pol,2),gel(pol,3)));
    return gerepilecopy(av, mkvec( basistoalg(nf,A) ));
  }
  A = fix_relative_pol(nf,pol,0);
  A = Q_primpart( lift_intern(A) );
  if (DEBUGLEVEL>3) fprintferr("test if polynomial is square-free\n");
  g = nfgcd(A, derivpol(A), T, gel(nf,4));

  if (degpol(g))
  { /* not squarefree */
    g = QXQX_normalize(g, T);
    A = RgXQX_div(A,g,T);
  }
  A = QXQX_normalize(A, T);
  A = Q_primpart(A);
  A = nfsqff(nf,A,1);
  A = RgXQV_to_mod(A, T);
  return gerepileupto(av, gen_sort(A, 0, cmp_pol));
}

/* assume x is squarefree */
int
nfissplit(GEN nf, GEN x)
{
  pari_sp av = avma;
  long l;
  if (typ(x) != t_POL) pari_err(typeer, "nfissplit");
  l = lg(nfsqff(checknf(nf), x, 2));
  avma = av; return l != 1;
}

/* nf = K[y] / (P) [P irreducible / K]. Is nf Galois over K ? */
int
nfisgalois(GEN nf, GEN P)
{
  if (typ(P) != t_POL) pari_err(typeer, "nfissplit");
  return degpol(P) <= 2 || nfissplit(nf, P);
}

/* return a minimal lift of elt modulo id */
static GEN
nf_bestlift(GEN elt, GEN bound, nflift_t *L)
{
  GEN u;
  long i,l = lg(L->prk), t = typ(elt);
  if (t != t_INT)
  {
    if (t == t_POL) elt = mulmat_pol(L->tozk, elt);
    u = gmul(L->iprk,elt);
    for (i=1; i<l; i++) gel(u,i) = diviiround(gel(u,i), L->den);
  }
  else
  {
    u = gmul(elt, gel(L->iprk,1));
    for (i=1; i<l; i++) gel(u,i) = diviiround(gel(u,i), L->den);
    elt = gscalcol(elt, l-1);
  }
  u = gsub(elt, gmul(L->prk, u));
  if (bound && gcmp(QuickNormL2(u,DEFAULTPREC), bound) > 0) u = NULL;
  return u;
}

/* Warning: return L->topowden * (best lift) */
static GEN
nf_bestlift_to_pol(GEN elt, GEN bound, nflift_t *L)
{
  pari_sp av = avma;
  GEN v, u = nf_bestlift(elt,bound,L);
  if (!u) return NULL;
  v = gclone(u); avma = av;
  u = gmul(L->topow, v); 
  gunclone(v); return u;
}

/* return the T->powden * (lift of pol with coefficients of T2-norm <= C)
 * if it exists */
static GEN
nf_pol_lift(GEN pol, GEN bound, nfcmbf_t *T)
{
  long i, l = lg(pol);
  GEN x = cgetg(l,t_POL);

  x[1] = pol[1];
  for (i=2; i<l; i++)
  {
    gel(x,i) = nf_bestlift_to_pol(gel(pol,i), bound, T->L);
    if (!x[i]) return NULL;
  }
  return x;
}

static GEN
zerofact(long v)
{
  GEN z = cgetg(3, t_MAT);
  gel(z,1) = mkcol(zeropol(v));
  gel(z,2) = mkcol(gen_1); return z;
}

/* return the factorization of x in nf */
GEN
nffactor(GEN nf,GEN pol)
{
  GEN A,g,y,p1,T, rep = cgetg(3, t_MAT);
  long l, j, dA;
  pari_sp av = avma;
  pari_timer ti;

  if (DEBUGLEVEL>2) { TIMERstart(&ti); fprintferr("\nEntering nffactor:\n"); }
  nf = checknf(nf); T = gel(nf,1);
  if (typ(pol) != t_POL) pari_err(notpoler,"nffactor");
  if (varncmp(varn(pol), varn(T)) >= 0)
    pari_err(talker,"polynomial variable must have highest priority in nffactor");

  A = fix_relative_pol(nf,pol,0);
  dA = degpol(A);
  if (dA <= 0) {
    avma = (pari_sp)(rep + 3);
    return dA == 0? trivfact(): zerofact(varn(pol));
  }
  A = Q_primpart( QXQX_normalize(A, T) );
  if (dA == 1) {
    GEN c;
    A = gerepilecopy(av, A); c = gel(A,2);
    if (typ(c) == t_POL && degpol(c) > 0) gel(A,2) = mkpolmod(c, gcopy(T));
    gel(rep,1) = mkcol(A);
    gel(rep,2) = mkcol(gen_1); return rep;
  }
  if (degpol(T) == 1)
    return gerepileupto(av, factpol(Q_primpart(simplify(pol)), 0));

  A = Q_primpart( lift_intern(A) );
  g = nfgcd(A, derivpol(A), T, gel(nf,4));

  A = QXQX_normalize(A, T);
  A = Q_primpart(A);
  if (DEBUGLEVEL>2) msgTIMER(&ti, "squarefree test");

  if (degpol(g))
  { /* not squarefree */
    pari_sp av1;
    GEN ex;
    g = QXQX_normalize(g, T);
    A = RgXQX_div(A,g, T);

    y = nfsqff(nf,A,0); av1 = avma;
    l = lg(y);
    ex=(GEN)gpmalloc(l * sizeof(long));
    for (j=l-1; j>=1; j--)
    {
      GEN fact = lift(gel(y,j)), quo = g, q;
      long e = 0;
      for(e = 1;; e++)
      {
        q = RgXQX_divrem(quo,fact,T, ONLY_DIVIDES);
        if (!q) break;
        quo = q;
      }
      ex[j] = e;
    }
    avma = av1; y = gerepileupto(av, RgXQXV_to_mod(y,T));
    p1 = cgetg(l, t_COL); for (j=l-1; j>=1; j--) gel(p1,j) = utoipos(ex[j]);
    free(ex);
  }
  else
  {
    y = gerepileupto(av, RgXQXV_to_mod(nfsqff(nf,A,0), T));
    l = lg(y);
    p1 = cgetg(l, t_COL); for (j=l-1; j>=1; j--) gel(p1,j) = gen_1;
  }
  if (DEBUGLEVEL>3)
    fprintferr("number of factor(s) found: %ld\n", lg(y)-1);
  gel(rep,1) = y;
  gel(rep,2) = p1; return sort_factor(rep, cmp_pol);
}

/* assume x scalar or t_COL */
static GEN
arch_for_T2(GEN G, GEN x)
{
  return (typ(x) == t_COL)? gmul(G,x): gmul(gel(G,1),x);
}

/* return a bound for T_2(P), P | polbase in C[X]
 * NB: Mignotte bound: A | S ==>
 *  |a_i| <= binom(d-1, i-1) || S ||_2 + binom(d-1, i) lc(S)
 *
 * Apply to sigma(S) for all embeddings sigma, then take the L_2 norm over
 * sigma, then take the sup over i.
 **/
static GEN
nf_Mignotte_bound(GEN nf, GEN polbase)
{
  GEN G = gmael(nf,5,2), lS = leading_term(polbase); /* t_INT */
  GEN p1, C, N2, matGS, binlS, bin;
  long prec, i, j, d = degpol(polbase), n = degpol(nf[1]), r1 = nf_get_r1(nf);

  binlS = bin = vecbinome(d-1);
  if (!gcmp1(lS)) binlS = gmul(lS, bin);

  N2 = cgetg(n+1, t_VEC);
  prec = gprecision(G);
  for (;;)
  {
    nffp_t F;

    matGS = cgetg(d+2, t_MAT);
    for (j=0; j<=d; j++) gel(matGS,j+1) = arch_for_T2(G, gel(polbase,j+2));
    matGS = shallowtrans(matGS);
    for (j=1; j <= r1; j++) /* N2[j] = || sigma_j(S) ||_2 */
    {
      gel(N2,j) = gsqrt( QuickNormL2(gel(matGS,j), DEFAULTPREC), DEFAULTPREC );
      if (lg(N2[j]) < DEFAULTPREC) goto PRECPB;
    }
    for (   ; j <= n; j+=2)
    {
      GEN q1 = QuickNormL2(gel(matGS,j  ), DEFAULTPREC);
      GEN q2 = QuickNormL2(gel(matGS,j+1), DEFAULTPREC);
      p1 = gmul2n(mpadd(q1, q2), -1);
      gel(N2,j) = gel(N2,j+1) = gsqrt( p1, DEFAULTPREC );
      if (lg(N2[j]) < DEFAULTPREC) goto PRECPB;
    }
    if (j > n) break; /* done */
PRECPB:
    prec = (prec<<1)-2;
    remake_GM(nf, &F, prec); G = F.G;
    if (DEBUGLEVEL>1) pari_warn(warnprec, "nf_factor_bound", prec);
  }

  /* Take sup over 0 <= i <= d of
   * sum_sigma | binom(d-1, i-1) ||sigma(S)||_2 + binom(d-1,i) lc(S) |^2 */

  /* i = 0: n lc(S)^2 */
  C = mulsi(n, sqri(lS));
  /* i = d: sum_sigma ||sigma(S)||_2^2 */
  p1 = gnorml2(N2); if (gcmp(C, p1) < 0) C = p1;
  for (i = 1; i < d; i++)
  {
    GEN s = gen_0;
    for (j = 1; j <= n; j++)
    {
      p1 = mpadd( mpmul(gel(bin,i), gel(N2,j)), gel(binlS,i+1) );
      s = mpadd(s, gsqr(p1));
    }
    if (gcmp(C, s) < 0) C = s;
  }
  return C;
}

/* return a bound for T_2(P), P | polbase
 * max |b_i|^2 <= 3^{3/2 + d} / (4 \pi d) [P]_2,
 * where [P]_2 is Bombieri's 2-norm
 * Sum over conjugates
*/
static GEN
nf_Beauzamy_bound(GEN nf, GEN polbase)
{
  GEN lt,C,run,s, G = gmael(nf,5,2), POL, bin;
  long i,prec,precnf, d = degpol(polbase), n = degpol(nf[1]);

  precnf = gprecision(G);
  prec   = MEDDEFAULTPREC;
  bin = vecbinome(d);
  POL = polbase + 2;
  /* compute [POL]_2 */
  for (;;)
  {
    run= real_1(prec);
    s = real_0(prec);
    for (i=0; i<=d; i++)
    {
      GEN p1 = gnorml2(arch_for_T2(G, gmul(run, gel(POL,i)))); /* T2(POL[i]) */
      if (!signe(p1)) continue;
      if (lg(p1) == 3) break;
      /* s += T2(POL[i]) / binomial(d,i) */
      s = addrr(s, gdiv(p1, gel(bin,i+1)));
    }
    if (i > d) break;

    prec = (prec<<1)-2;
    if (prec > precnf)
    {
      nffp_t F; remake_GM(nf, &F, prec); G = F.G;
      if (DEBUGLEVEL>1) pari_warn(warnprec, "nf_factor_bound", prec);
    }
  }
  lt = leading_term(polbase);
  s = gmul(s, mulis(sqri(lt), n));
  C = powrshalf(stor(3,DEFAULTPREC), 3 + 2*d); /* 3^{3/2 + d} */
  return gdiv(gmul(C, s), gmulsg(d, mppi(DEFAULTPREC)));
}

static GEN
nf_factor_bound(GEN nf, GEN polbase)
{
  pari_sp av = avma;
  GEN a = nf_Mignotte_bound(nf, polbase);
  GEN b = nf_Beauzamy_bound(nf, polbase);
  if (DEBUGLEVEL>2)
  {
    fprintferr("Mignotte bound: %Z\n",a);
    fprintferr("Beauzamy bound: %Z\n",b);
  }
  return gerepileupto(av, gmin(a, b));
}

static long
ZXY_get_prec(GEN P)
{
  long i, j, z, prec = 0;
  for (i=2; i<lg(P); i++)
  {
    GEN p = gel(P,i);
    if (typ(p) == t_INT)
    {
      z = lgefint(p);
      if (z > prec) prec = z;
    }
    else
    {
      for (j=2; j<lg(p); j++)
      {
        z = lgefint(p[j]);
        if (z > prec) prec = z;
      }
    }
  }
  return prec + 1;
}

long
ZM_get_prec(GEN x)
{
  long i, j, l, k = 2, lx = lg(x);

  for (j=1; j<lx; j++)
  {
    GEN c = gel(x,j);
    for (i=1; i<lx; i++) { l = lgefint(c[i]); if (l > k) k = l; }
  }
  return k;
}
long
ZX_get_prec(GEN x)
{
  long j, l, k = 2, lx = lg(x);

  for (j=2; j<lx; j++)
  {
    l = lgefint(x[j]); if (l > k) k = l;
  }
  return k;
}

/* return Bs: if r a root of sigma_i(P), |r| < Bs[i] */
static GEN
nf_root_bounds(GEN P, GEN T)
{
  long lR, i, j, l, prec;
  GEN Ps, R, V, nf;

  if (RgX_is_rational(P)) return logmax_modulus_bound(P);
  T = get_nfpol(T, &nf);

  P = Q_primpart(P);
  prec = ZXY_get_prec(P);
  l = lg(P);
  if (nf && nfgetprec(nf) >= prec)
    R = gel(nf,6);
  else
    R = roots(T, prec);
  lR = lg(R);
  V = cgetg(lR, t_VEC);
  Ps = cgetg(l, t_POL); /* sigma (P) */
  Ps[1] = P[1];
  for (j=1; j<lg(R); j++)
  {
    GEN r = gel(R,j);
    for (i=2; i<l; i++) gel(Ps,i) = poleval(gel(P,i), r);
    gel(V,j) = logmax_modulus_bound(Ps);
  }
  return V;
}

/* return B such that if x in O_K, K = Z[X]/(T), then the L2-norm of the
 * coordinates of the numerator of x [on the power, resp. integral, basis if T
 * is a polynomial, resp. an nf] is  <= B T_2(x)
 * *ptden = multiplicative bound for denom(x) */
static GEN
L2_bound(GEN T, GEN tozk, GEN *ptden)
{
  GEN nf, M, L, prep, den;
  long prec;

  T = get_nfpol(T, &nf);
  prec = ZX_get_prec(T) + ZM_get_prec(tozk);
  den = nf? gen_1: NULL;
  den = initgaloisborne(T, den, prec, &L, &prep, NULL);
  M = vandermondeinverse(L, gmul(T, real_1(prec)), den, prep);
  if (nf) M = gmul(tozk, M);
  if (gcmp1(den)) den = NULL;
  *ptden = den;
  return QuickNormL2(M, DEFAULTPREC);
}

/* || L ||_p^p in dimension n (L may be a scalar) */
static GEN
normlp(GEN L, long p, long n)
{
  long i,l, t = typ(L);
  GEN z;

  if (!is_vec_t(t)) return gmulsg(n, gpowgs(L, p));

  l = lg(L); z = gen_0;
  /* assert(n == l-1); */
  for (i=1; i<l; i++)
    z = gadd(z, gpowgs(gel(L,i), p));
  return z;
}

/* S = S0 + tS1, P = P0 + tP1 (Euclidean div. by t integer). For a true
 * factor (vS, vP), we have:
 *    | S vS + P vP |^2 < Btra
 * This implies | S1 vS + P1 vP |^2 < Bhigh, assuming t > sqrt(Btra).
 * d = dimension of low part (= [nf:Q])
 * n0 = bound for |vS|^2
 * */
static double
get_Bhigh(long n0, long d)
{
  double sqrtd = sqrt((double)d);
  double z = n0*sqrtd + sqrtd/2 * (d * (n0+1));
  z = 1. + 0.5 * z; return z * z;
}

typedef struct {
  GEN d;
  GEN dPinvS;   /* d P^(-1) S   [ integral ] */
  double **PinvSdbl; /* P^(-1) S as double */
  GEN S1, P1;   /* S = S0 + S1 q, idem P */
} trace_data;

/* S1 * u - P1 * round(P^-1 S u). K non-zero coords in u given by ind */
static GEN
get_trace(GEN ind, trace_data *T)
{
  long i, j, l, K = lg(ind)-1;
  GEN z, s, v;

  s = gel(T->S1, ind[1]);
  if (K == 1) return s;

  /* compute s = S1 u */
  for (j=2; j<=K; j++) s = gadd(s, gel(T->S1, ind[j]));

  /* compute v := - round(P^1 S u) */
  l = lg(s);
  v = cgetg(l, t_VECSMALL);
  for (i=1; i<l; i++)
  {
    double r, t = 0.;
    /* quick approximate computation */
    for (j=1; j<=K; j++) t += T->PinvSdbl[ ind[j] ][i];
    r = floor(t + 0.5);
    if (fabs(t + 0.5 - r) < 0.0001)
    { /* dubious, compute exactly */
      z = gen_0;
      for (j=1; j<=K; j++) z = addii(z, ((GEN**)T->dPinvS)[ ind[j] ][i]);
      v[i] = - itos( diviiround(z, T->d) );
    }
    else
      v[i] = - (long)r;
  }
  return gadd(s, ZM_zc_mul(T->P1, v));
} 

static trace_data *
init_trace(trace_data *T, GEN S, nflift_t *L, GEN q)
{
  long e = gexpo(S), i,j, l,h;
  GEN qgood, S1, invd;

  if (e < 0) return NULL; /* S = 0 */

  qgood = int2n(e - 32); /* single precision check */
  if (cmpii(qgood, q) > 0) q = qgood;

  S1 = gdivround(S, q);
  if (gcmp0(S1)) return NULL;

  invd = ginv(itor(L->den, DEFAULTPREC));

  T->dPinvS = gmul(L->iprk, S);
  l = lg(S);
  h = lg(T->dPinvS[1]);
  T->PinvSdbl = (double**)cgetg(l, t_MAT);
  init_dalloc();
  for (j = 1; j < l; j++)
  {
    double *t = dalloc(h * sizeof(double));
    GEN c = gel(T->dPinvS,j);
    pari_sp av = avma;
    T->PinvSdbl[j] = t;
    for (i=1; i < h; i++) t[i] = rtodbl(mpmul(invd, gel(c,i)));
    avma = av;
  }

  T->d  = L->den;
  T->P1 = gdivround(L->prk, q);
  T->S1 = S1; return T;
}

static void
update_trace(trace_data *T, long k, long i)
{
  if (!T) return;
  T->S1[k]       = T->S1[i];
  T->dPinvS[k]   = T->dPinvS[i];
  T->PinvSdbl[k] = T->PinvSdbl[i];
}

/* reduce coeffs mod (T,pk), then center mod pk */
static GEN
FqX_centermod(GEN z, GEN T, GEN pk, GEN pks2)
{
  long i, l;
  GEN y;
  if (!T) return centermod_i(z, pk, pks2);
  y = FpXQX_red(z, T, pk); l = lg(y);
  for (i = 2; i < l; i++)
  {
    GEN c = gel(y,i);
    if (typ(c) == t_INT)
      c = centermodii(c, pk, pks2);
    else
      c = FpX_center(c, pk);
    gel(y,i) = c;
  }
  return y;
}

/* Naive recombination of modular factors: combine up to maxK modular
 * factors, degree <= klim and divisible by hint
 *
 * target = polynomial we want to factor
 * famod = array of modular factors.  Product should be congruent to
 * target/lc(target) modulo p^a
 * For true factors: S1,S2 <= p^b, with b <= a and p^(b-a) < 2^31 */
static GEN
nfcmbf(nfcmbf_t *T, GEN p, long a, long maxK, long klim)
{
  GEN pol = T->pol, nf = T->nf, famod = T->fact, dn = T->dn;
  GEN bound = T->bound;
  GEN nfpol = gel(nf,1);
  long K = 1, cnt = 1, i,j,k, curdeg, lfamod = lg(famod)-1, dnf = degpol(nfpol);
  GEN res = cgetg(3, t_VEC);
  pari_sp av0 = avma;
  GEN pk = gpowgs(p,a), pks2 = shifti(pk,-1);

  GEN ind      = cgetg(lfamod+1, t_VECSMALL);
  GEN degpol   = cgetg(lfamod+1, t_VECSMALL);
  GEN degsofar = cgetg(lfamod+1, t_VECSMALL);
  GEN listmod  = cgetg(lfamod+1, t_COL);
  GEN fa       = cgetg(lfamod+1, t_COL);
  GEN lc = absi(leading_term(pol)), lt = is_pm1(lc)? NULL: lc;
  GEN C2ltpol, C = T->L->topowden, Tpk = T->L->Tpk;
  GEN Clt  = mul_content(C, lt);
  GEN C2lt = mul_content(C,Clt);
  const double Bhigh = get_Bhigh(lfamod, dnf);
  trace_data _T1, _T2, *T1, *T2;
  pari_timer ti;

  TIMERstart(&ti);

  if (maxK < 0) maxK = lfamod-1;

  C2ltpol = C2lt? gmul(C2lt,pol): pol;
  {
    GEN q = ceil_safe(sqrtr(T->BS_2));
    GEN t1,t2, ltdn, lt2dn;
    GEN trace1   = cgetg(lfamod+1, t_MAT);
    GEN trace2   = cgetg(lfamod+1, t_MAT);

    ltdn = mul_content(lt, dn);
    lt2dn= mul_content(ltdn, lt);

    for (i=1; i <= lfamod; i++)
    {
      pari_sp av = avma;
      GEN P = gel(famod,i);
      long d = degpol(P);

      degpol[i] = d; P += 2;
      t1 = gel(P,d-1);/* = - S_1 */
      t2 = gsqr(t1);
      if (d > 1) t2 = gsub(t2, gmul2n(gel(P,d-2), 1));
      /* t2 = S_2 Newton sum */
      t2 = typ(t2)!=t_INT? FpX_rem(t2, Tpk, pk): modii(t2, pk);
      if (lt)
      {
        if (typ(t2)!=t_INT) {
          t1 = FpX_red(gmul(ltdn, t1), pk);
          t2 = FpX_red(gmul(lt2dn,t2), pk);
        } else {
          t1 = remii(mulii(ltdn, t1), pk);
          t2 = remii(mulii(lt2dn,t2), pk);
        }
      }
      gel(trace1,i) = gclone( nf_bestlift(t1, NULL, T->L) );
      gel(trace2,i) = gclone( nf_bestlift(t2, NULL, T->L) ); avma = av;
    }
    T1 = init_trace(&_T1, trace1, T->L, q);
    T2 = init_trace(&_T2, trace2, T->L, q);
    for (i=1; i <= lfamod; i++) { 
      gunclone(gel(trace1,i));
      gunclone(gel(trace2,i));
    }
  }
  degsofar[0] = 0; /* sentinel */

  /* ind runs through strictly increasing sequences of length K,
   * 1 <= ind[i] <= lfamod */
nextK:
  if (K > maxK || 2*K > lfamod) goto END;
  if (DEBUGLEVEL > 3)
    fprintferr("\n### K = %d, %Z combinations\n", K,binomial(utoipos(lfamod), K));
  setlg(ind, K+1); ind[1] = 1;
  i = 1; curdeg = degpol[ind[1]];
  for(;;)
  { /* try all combinations of K factors */
    for (j = i; j < K; j++)
    {
      degsofar[j] = curdeg;
      ind[j+1] = ind[j]+1; curdeg += degpol[ind[j+1]];
    }
    if (curdeg <= klim && curdeg % T->hint == 0) /* trial divide */
    {
      GEN t, y, q, list;
      pari_sp av;

      av = avma;
      /* d - 1 test */
      if (T1)
      {
        t = get_trace(ind, T1);
        if (rtodbl(QuickNormL2(t,DEFAULTPREC)) > Bhigh)
        {
          if (DEBUGLEVEL>6) fprintferr(".");
          avma = av; goto NEXT;
        }
      }
      /* d - 2 test */
      if (T2)
      {
        t = get_trace(ind, T2);
        if (rtodbl(QuickNormL2(t,DEFAULTPREC)) > Bhigh)
        {
          if (DEBUGLEVEL>3) fprintferr("|");
          avma = av; goto NEXT;
        }
      }
      avma = av;
      y = lt; /* full computation */
      for (i=1; i<=K; i++)
      {
        GEN q = gel(famod, ind[i]);
        if (y) q = gmul(y, q);
        y = FqX_centermod(q, Tpk, pk, pks2);
      }
      y = nf_pol_lift(y, bound, T);
      if (!y)
      {
        if (DEBUGLEVEL>3) fprintferr("@");
        avma = av; goto NEXT;
      }
      /* try out the new combination: y is the candidate factor */
      q = RgXQX_divrem(C2ltpol, y, nfpol, ONLY_DIVIDES);
      if (!q)
      {
        if (DEBUGLEVEL>3) fprintferr("*");
        avma = av; goto NEXT;
      }

      /* found a factor */
      list = cgetg(K+1, t_VEC);
      gel(listmod,cnt) = list;
      for (i=1; i<=K; i++) list[i] = famod[ind[i]];

      y = Q_primpart(y);
      gel(fa,cnt++) = QXQX_normalize(y, nfpol);
      /* fix up pol */
      pol = q;
      for (i=j=k=1; i <= lfamod; i++)
      { /* remove used factors */
        if (j <= K && i == ind[j]) j++;
        else
        {
          famod[k] = famod[i];
          update_trace(T1, k, i);
          update_trace(T2, k, i);
          degpol[k] = degpol[i]; k++;
        }
      }
      lfamod -= K;
      if (lfamod < 2*K) goto END;
      i = 1; curdeg = degpol[ind[1]];

      if (C2lt) pol = Q_primpart(pol);
      if (lt) lt = absi(leading_term(pol));
      Clt  = mul_content(C, lt);
      C2lt = mul_content(C,Clt);
      C2ltpol = C2lt? gmul(C2lt,pol): pol;
      if (DEBUGLEVEL > 2)
      {
        fprintferr("\n"); msgTIMER(&ti, "to find factor %Z",y);
        fprintferr("remaining modular factor(s): %ld\n", lfamod);
      }
      continue;
    }

NEXT:
    for (i = K+1;;)
    {
      if (--i == 0) { K++; goto nextK; }
      if (++ind[i] <= lfamod - K + i)
      {
        curdeg = degsofar[i-1] + degpol[ind[i]];
        if (curdeg <= klim) break;
      }
    }
  }
END:
  if (degpol(pol) > 0)
  { /* leftover factor */
    GEN LC = leading_term(pol);
    if (signe(LC) < 0) { pol = gneg_i(pol); LC = leading_term(pol); }

    if (lfamod >= 2*K) {
      if (!equalii(LC, lc)) pol = RgX_Rg_mul(pol, gdiv(lc,LC));
    } else {
      if (C2lt) pol = QXQX_normalize(Q_primpart(pol), nfpol);
    }
    setlg(famod, lfamod+1);
    gel(listmod,cnt) = shallowcopy(famod);
    gel(fa,cnt++) = pol;
  }
  if (DEBUGLEVEL>6) fprintferr("\n");
  if (cnt == 2) { 
    avma = av0; 
    gel(res,1) = mkvec(T->pol);
    gel(res,2) = mkvec(T->fact);
  }
  else
  {
    setlg(listmod, cnt); setlg(fa, cnt);
    gel(res,1) = fa;
    gel(res,2) = listmod;
    res = gerepilecopy(av0, res);
  }
  return res;
}

static GEN
nf_chk_factors(nfcmbf_t *T, GEN P, GEN M_L, GEN famod, GEN pk)
{
  GEN nf = T->nf, bound = T->bound;
  GEN nfT = gel(nf,1);
  long i, r;
  GEN pol = P, list, piv, y;
  GEN C2ltpol, C = T->L->topowden, Tpk = T->L->Tpk;
  GEN lc = absi(leading_term(pol)), lt = is_pm1(lc)? NULL: lc;
  GEN Clt  = mul_content(C, lt);
  GEN C2lt = mul_content(C,Clt);

  piv = special_pivot(M_L);
  if (!piv) return NULL;
  if (DEBUGLEVEL>3) fprintferr("special_pivot output:\n%Z\n",piv);

  r  = lg(piv)-1;
  list = cgetg(r+1, t_COL);
  C2ltpol = C2lt? gmul(C2lt,pol): pol;
  for (i = 1;;)
  {
    pari_sp av = avma;
    if (DEBUGLEVEL)
      fprintferr("nf_LLL_cmbf: checking factor %ld (avma - bot = %lu)\n",
                 i, avma - bot);
    y = chk_factors_get(lt, famod, gel(piv,i), Tpk, pk);
    if (DEBUGLEVEL>2) fprintferr("... mod p^k (avma - bot = %lu)\n", avma-bot);
    
    if (! (y = nf_pol_lift(y, bound, T)) ) return NULL;
    if (DEBUGLEVEL>2) fprintferr("... lifted (avma - bot = %lu)\n", avma-bot);

    y = gerepilecopy(av, y);
    /* y is the candidate factor */
    pol = RgXQX_divrem(C2ltpol, y, nfT, ONLY_DIVIDES);
    if (!pol) return NULL;

    y = Q_primpart(y);
    gel(list,i) = QXQX_normalize(y, nfT);
    if (++i >= r) break;

    if (C2lt) pol = Q_primpart(pol);
    if (lt) lt = absi(leading_term(pol));
    Clt  = mul_content(C, lt);
    C2lt = mul_content(C,Clt);
    C2ltpol = C2lt? gmul(C2lt,pol): pol;
  }
  y = Q_primpart(pol);
  gel(list,i) = QXQX_normalize(y, nfT); return list;
}

static GEN
nf_to_Zq(GEN x, GEN T, GEN pk, GEN pks2, GEN proj)
{
  GEN y;
  if (typ(x) != t_COL) return centermodii(x, pk, pks2);
  y = gmul(proj, x);
  if (!T) return centermodii(y, pk, pks2);
  y = RgV_to_RgX(y, varn(T));
  return centermod_i(FpX_rem(y, T, pk), pk, pks2);
}

/* Assume P in nfX form [ unifpol(,t_COL) ], lc(P) != 0 mod p.
 * Reduce P to Zp[X]/(T) mod p^a */
static GEN
ZqX(GEN P, GEN pk, GEN T, GEN proj)
{
  long i, l = lg(P);
  GEN z, pks2 = shifti(pk,-1);

  z = cgetg(l,t_POL); z[1] = P[1];
  for (i=2; i<l; i++) gel(z,i) = nf_to_Zq(gel(P,i),T,pk,pks2,proj);
  return normalizepol(z);
}

static GEN
ZqX_normalize(GEN P, GEN lt, nflift_t *L)
{
  GEN R = lt? gmul(Fp_inv(lt, L->pk), P): P;
  return ZqX(R, L->pk, L->Tpk, L->ZqProj);
}

/* We want to be able to reconstruct x, |x|^2 < C, from x mod pr^k */
static double
bestlift_bound(GEN C, long d, double alpha, GEN Npr)
{
  const double y = 1 / (alpha - 0.25); /* = 2 if alpha = 3/4 */
  double t;
  if (typ(C) != t_REAL) C = gmul(C, real_1(DEFAULTPREC));
  setlg(C, DEFAULTPREC);
  t = rtodbl(mplog(gmul2n(divrs(C,d), 4))) * 0.5 + (d-1) * log(1.5 * sqrt(y));
  return ceil((t * d) / log(gtodouble(Npr)));
}

static GEN
get_R(GEN M)
{ 
  GEN R;
  long i, l, prec = DEFAULTPREC + (gexpo(M) >> TWOPOTBITS_IN_LONG);

  for(;;)
  {
    R = sqred1_from_QR(M, prec);
    if (R) break;
    prec = (prec-1)<<1;
  }
  l = lg(R);
  for (i=1; i<l; i++) gcoeff(R,i,i) = gen_1;
  return R;
}

static void
init_proj(nflift_t *L, GEN nfT, GEN p)
{
  if (L->Tp)
  {
    GEN z = cgetg(3, t_VEC), proj;
    gel(z,1) = L->Tp;
    gel(z,2) = FpX_div(FpX_red(nfT,p), L->Tp, p);
    z = hensel_lift_fact(nfT, z, NULL, p, L->pk, L->k);
    L->Tpk = gel(z,1);
    proj = get_proj_modT(L->topow, L->Tpk, L->pk);
    if (L->topowden)
      proj = FpM_red(gmul(Fp_inv(L->topowden, L->pk), proj), L->pk);
    L->ZqProj = proj;
  }
  else
  {
    L->Tpk = NULL;
    L->ZqProj = dim1proj(L->prkHNF);
  }
}

static void
bestlift_init(long a, GEN nf, GEN pr, GEN C, nflift_t *L)
{
  const long D = 100;
  const double alpha = ((double)D-1) / D; /* LLL parameter */
  const long d = degpol(nf[1]);
  pari_sp av = avma;
  GEN prk, PRK, B, GSmin, pk;
  pari_timer ti;

  TIMERstart(&ti);
  if (!a) a = (long)bestlift_bound(C, d, alpha, pr_norm(pr));

  for (;; avma = av, a<<=1)
  {
    if (DEBUGLEVEL>2) fprintferr("exponent: %ld\n",a);
    PRK = prk = idealpows(nf, pr, a);
    pk = gcoeff(prk,1,1);
    /* reduce size first, "scramble" matrix */
    PRK = lllintpartial_ip(PRK);
    /* now floating point reduction is fast */
    PRK = lllint_fp_ip(PRK, 4);
    PRK = lllint_i(PRK, D, 0, NULL, NULL, &B);
    if (!PRK) { PRK = prk; GSmin = pk; } /* nf = Q */
    else
    {
      pari_sp av2 = avma;
      GEN S = invmat( get_R(PRK) ), BB = GS_norms(B, DEFAULTPREC);
      GEN smax = gen_0;
      long i, j;
      for (i=1; i<=d; i++)
      {
        GEN s = gen_0;
        for (j=1; j<=d; j++)
          s = gadd(s, gdiv( gsqr(gcoeff(S,i,j)), gel(BB,j)));
        if (gcmp(s, smax) > 0) smax = s;
      }
      GSmin = gerepileupto(av2, ginv(gmul2n(smax, 2)));
    }
    if (gcmp(GSmin, C) >= 0) break;
  }
  if (DEBUGLEVEL>2)
    fprintferr("for this exponent, GSmin = %Z\nTime reduction: %ld\n",
      GSmin, TIMER(&ti));
  L->k = a;
  L->den = L->pk = pk;
  L->prk = PRK;
  L->iprk = ZM_inv(PRK, pk);
  L->GSmin= GSmin;
  L->prkHNF = prk;
  init_proj(L, gel(nf,1), gel(pr,1));
}

/* Let X = Tra * M_L, Y = bestlift(X) return V s.t Y = X - PRK V
 * and set *eT2 = gexpo(Y)  [cf nf_bestlift, but memory efficient] */
static GEN
get_V(GEN Tra, GEN M_L, GEN PRK, GEN PRKinv, GEN pk, long *eT2)
{
  long i, e = 0, l = lg(M_L);
  GEN V = cgetg(l, t_MAT);
  *eT2 = 0;
  for (i = 1; i < l; i++)
  { /* cf nf_bestlift(Tra * c) */
    pari_sp av = avma, av2;
    GEN v, T2 = gmul(Tra, gel(M_L,i));
      
    v = gdivround(gmul(PRKinv, T2), pk); /* small */
    av2 = avma;
    T2 = gsub(T2, gmul(PRK, v));
    e = gexpo(T2); if (e > *eT2) *eT2 = e;
    avma = av2;
    gel(V,i) = gerepileupto(av, v); /* small */
  }
  return V;
}

static GEN
nf_LLL_cmbf(nfcmbf_t *T, GEN p, long k, long rec)
{
  nflift_t *L = T->L;
  GEN pk = L->pk, PRK = L->prk, PRKinv = L->iprk, GSmin = L->GSmin;
  GEN Tpk = L->Tpk;

  GEN famod = T->fact, nf = T->nf, ZC = T->ZC, Br = T->Br;
  GEN Pbase = T->polbase, P = T->pol, dn = T->dn;
  GEN nfT = gel(nf,1);
  GEN Btra;
  long dnf = degpol(nfT), dP = degpol(P);

  double BitPerFactor = 0.5; /* nb bits / modular factor */
  long i, C, tmax, n0;
  GEN lP, Bnorm, Tra, T2, TT, CM_L, m, list, ZERO;
  double Bhigh;
  pari_sp av, av2, lim;
  long ti_LLL = 0, ti_CF = 0;
  pari_timer ti2, TI;

  lP = absi(leading_term(P));
  if (is_pm1(lP)) lP = NULL;

  n0 = lg(famod) - 1;
 /* Lattice: (S PRK), small vector (vS vP). To find k bound for the image,
  * write S = S1 q + S0, P = P1 q + P0
  * |S1 vS + P1 vP|^2 <= Bhigh for all (vS,vP) assoc. to true factors */
  Btra = mulrr(ZC, mulsr(dP*dP, normlp(Br, 2, dnf)));
  Bhigh = get_Bhigh(n0, dnf);
  C = (long)ceil(sqrt(Bhigh/n0)) + 1; /* C^2 n0 ~ Bhigh */
  Bnorm = dbltor( n0 * C * C + Bhigh );
  ZERO = zeromat(n0, dnf);

  av = avma; lim = stack_lim(av, 1);
  TT = cgetg(n0+1, t_VEC);
  Tra  = cgetg(n0+1, t_MAT);
  for (i=1; i<=n0; i++) TT[i] = 0;
  CM_L = gscalsmat(C, n0);
  /* tmax = current number of traces used (and computed so far) */
  for(tmax = 0;; tmax++)
  {
    long a, b, bmin, bgood, delta, tnew = tmax + 1, r = lg(CM_L)-1;
    GEN oldCM_L, M_L, q, S1, P1, VV;
    int first = 1;

    /* bound for f . S_k(genuine factor) = ZC * bound for T_2(S_tnew) */
    Btra = mulrr(ZC, mulsr(dP*dP, normlp(Br, 2*tnew, dnf)));
    bmin = logint(ceil_safe(sqrtr(Btra)), gen_2, NULL);
    if (DEBUGLEVEL>2)
      fprintferr("\nLLL_cmbf: %ld potential factors (tmax = %ld, bmin = %ld)\n",
                 r, tmax, bmin);

    /* compute Newton sums (possibly relifting first) */
    if (gcmp(GSmin, Btra) < 0)
    {
      nflift_t L1;
      GEN polred;

      bestlift_init(k<<1, nf, T->pr, Btra, &L1);
      polred = ZqX_normalize(Pbase, lP, &L1);
      k      = L1.k;
      pk     = L1.pk;
      PRK    = L1.prk;
      PRKinv = L1.iprk;
      GSmin  = L1.GSmin;
      Tpk    = L1.Tpk;
      famod = hensel_lift_fact(polred, famod, Tpk, p, pk, k);
      for (i=1; i<=n0; i++) TT[i] = 0;
    }
    for (i=1; i<=n0; i++)
    {
      GEN h, lPpow = lP? gpowgs(lP, tnew): NULL;
      GEN z = polsym_gen(gel(famod,i), gel(TT,i), tnew, Tpk, pk);
      gel(TT,i) = z;
      h = gel(z,tnew+1);
      /* make Newton sums integral */
      lPpow = mul_content(lPpow, dn);
      if (lPpow)
        h = (typ(h) == t_INT)? modii(mulii(h, lPpow), L->pk): 
                               FpX_Fp_mul(h, lPpow, L->pk);
      gel(Tra,i) = nf_bestlift(h, NULL, L); /* S_tnew(famod) */
    }

    /* compute truncation parameter */
    if (DEBUGLEVEL>2) { TIMERstart(&ti2); TIMERstart(&TI); }
    oldCM_L = CM_L;
    av2 = avma;
    b = delta = 0; /* -Wall */
AGAIN:
    M_L = Q_div_to_int(CM_L, utoipos(C));
    VV = get_V(Tra, M_L, PRK, PRKinv, pk, &a);
    if (first)
    { /* initialize lattice, using few p-adic digits for traces */
      bgood = (long)(a - max(32, BitPerFactor * r));
      b = max(bmin, bgood);
      delta = a - b;
    }
    else
    { /* add more p-adic digits and continue reduction */
      if (a < b) b = a;
      b = max(b-delta, bmin);
      if (b - delta/2 < bmin) b = bmin; /* near there. Go all the way */
    }

    /* restart with truncated entries */
    q = int2n(b);
    P1 = gdivround(PRK, q);
    S1 = gdivround(Tra, q);
    T2 = gsub(gmul(S1, M_L), gmul(P1, VV));
    m = vconcat( CM_L, T2 );
    if (first)
    {
      first = 0;
      m = shallowconcat( m, vconcat(ZERO, P1) );
      /*     [ C M_L   0  ]
       * m = [            ]   square matrix
       *     [  T2'   PRK ]   T2' = Tra * M_L  truncated
       */
    }
    CM_L = LLL_check_progress(Bnorm, n0, m, b == bmin, /*dbg:*/ &ti_LLL);
    if (DEBUGLEVEL>2)
      fprintferr("LLL_cmbf: (a,b) =%4ld,%4ld; r =%3ld -->%3ld, time = %ld\n",
                 a,b, lg(m)-1, CM_L? lg(CM_L)-1: 1, TIMER(&TI));
    if (!CM_L) { list = mkcol(QXQX_normalize(P,nfT)); break; }
    if (b > bmin)
    {
      CM_L = gerepilecopy(av2, CM_L);
      goto AGAIN;
    }
    if (DEBUGLEVEL>2) msgTIMER(&ti2, "for this trace");

    i = lg(CM_L) - 1;
    if (i == r && gequal(CM_L, oldCM_L))
    {
      CM_L = oldCM_L;
      avma = av2; continue;
    }

    if (i <= r && i*rec < n0)
    {
      pari_timer ti;
      if (DEBUGLEVEL>2) TIMERstart(&ti);
      list = nf_chk_factors(T, P, Q_div_to_int(CM_L,utoipos(C)), famod, pk);
      if (DEBUGLEVEL>2) ti_CF += TIMER(&ti);
      if (list) break;
      CM_L = gerepilecopy(av2, CM_L);
    }
    if (low_stack(lim, stack_lim(av,1)))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"nf_LLL_cmbf");
      gerepileall(av, Tpk? 9: 8,
                      &CM_L,&TT,&Tra,&famod,&pk,&GSmin,&PRK,&PRKinv,&Tpk);
    }
  }
  if (DEBUGLEVEL>2)
    fprintferr("* Time LLL: %ld\n* Time Check Factor: %ld\n",ti_LLL,ti_CF);
  return list;
}

static GEN
nf_combine_factors(nfcmbf_t *T, GEN polred, GEN p, long a, long klim)
{
  GEN res, L, listmod, famod = T->fact, nf = T->nf;
  long l, maxK = 3, nft = lg(famod)-1;
  pari_timer ti;

  if (DEBUGLEVEL>2) TIMERstart(&ti);
  T->fact = hensel_lift_fact(polred, famod, T->L->Tpk, p, T->L->pk, a);
  if (nft < 11) maxK = -1; /* few modular factors: try all posibilities */
  if (DEBUGLEVEL>2) msgTIMER(&ti, "Hensel lift");

  L = nfcmbf(T, p, a, maxK, klim);
  if (DEBUGLEVEL>2) msgTIMER(&ti, "Naive recombination");

  res     = gel(L,1);
  listmod = gel(L,2); l = lg(listmod)-1;
  famod = gel(listmod,l);
  if (maxK >= 0 && lg(famod)-1 > 2*maxK)
  {
    if (l > 1)
    {
      T->pol = gel(res,l);
      T->polbase = unifpol(nf, gel(res,l), t_COL);
    }
    L = nf_LLL_cmbf(T, p, a, maxK);
    /* remove last elt, possibly unfactored. Add all new ones. */
    setlg(res, l); res = shallowconcat(res, L);
  }
  return res;
}

static GEN
nf_DDF_roots(GEN pol, GEN polred, GEN nfpol, GEN lt, GEN init_fa, long nbf,
             long fl, nflift_t *L)
{
  long Cltx_r[] = { evaltyp(t_POL)|_evallg(4), 0,0,0 };
  long i, m;
  GEN C2ltpol, C = L->topowden;
  GEN Clt  = mul_content(C, lt);
  GEN C2lt = mul_content(C,Clt);
  GEN z;

  if (L->Tpk)
  {
    int cof = (degpol(pol) > nbf); /* non trivial cofactor ? */
    z = FqX_split_roots(init_fa, L->Tp, L->p, cof? polred: NULL);
    z = hensel_lift_fact(polred, z, L->Tpk, L->p, L->pk, L->k);
    if (cof) setlg(z, lg(z)-1); /* remove cofactor */
    z = roots_from_deg1(z);
  }
  else
    z = rootpadicfast(polred, L->p, L->k);
  Cltx_r[1] = evalsigne(1) | evalvarn(varn(pol));
  gel(Cltx_r,3) = Clt? Clt: gen_1;
  C2ltpol  = C2lt? gmul(C2lt, pol): pol;
  for (m=1,i=1; i<lg(z); i++)
  {
    GEN q, r = gel(z,i);

    r = nf_bestlift_to_pol(lt? gmul(lt,r): r, NULL, L);
    gel(Cltx_r,2) = gneg(r); /* check P(r) == 0 */
    q = RgXQX_divrem(C2ltpol, Cltx_r, nfpol, ONLY_DIVIDES); /* integral */
    if (q) { 
      C2ltpol = C2lt? gmul(Clt,q): q;
      if (Clt) r = gdiv(r, Clt);
      gel(z,m++) = r;
    }
    else if (fl == 2) return cgetg(1, t_VEC);
  }
  z[0] = evaltyp(t_VEC) | evallg(m);
  return z;
}

static long
nf_pick_prime(long ct, GEN nf, GEN polbase, long fl,
              GEN *lt, GEN *Fa, GEN *pr, GEN *Tp)
{
  GEN nfpol = gel(nf,1), dk, bad;
  long maxf, n = degpol(nfpol), dpol = degpol(polbase), nbf = 0;
  byteptr pt = diffptr;
  ulong pp = 0;

  *lt  = leading_term(polbase); /* t_INT */
  if (gcmp1(*lt)) *lt = NULL;
  dk = absi(gel(nf,3));
  bad = mulii(dk,gel(nf,4));

  /* FIXME: slow factorization of large polynomials over large Fq */
  maxf = 1;
  if (ct > 1) {
    if (dpol > 100) /* tough */
    {
      if (n >= 20) maxf = 4;
    }
    else
    {
      if (n >= 15) maxf = 4;
    }
  }
  
  for (ct = 5;;)
  {
    GEN aT, apr, ap, modpr, red;
    long anbf;
    ulong ltp = 0;
    pari_timer ti_pr;

    GEN list, r = NULL, fa = NULL;
    pari_sp av2 = avma;
    if (DEBUGLEVEL>3) TIMERstart(&ti_pr);
    for (;;)
    {
      NEXT_PRIME_VIADIFF_CHECK(pp, pt);
      if (! umodiu(bad,pp)) continue;
      if (*lt) { ltp = umodiu(*lt, pp); if (!ltp) continue; }
      ap = utoipos(pp);
      list = (GEN)FpX_factor(nfpol, ap)[1];
      if (maxf == 1)
      { /* deg.1 factors are best */
        r = gel(list,1);
        if (degpol(r) == 1) break;
      }
      else
      { /* otherwise, pick factor of largish degree */
        long i, dr;
        for (i = lg(list)-1; i > 0; i--)
        {
          r = gel(list,i); dr = degpol(r);
          if (dr <= maxf) break;
        }
        if (i > 0) break;
      }
      avma = av2;
    }
    apr = primedec_apply_kummer(nf,r,1,ap);

    modpr = zk_to_ff_init(nf,&apr,&aT,&ap);
    red = modprX(polbase, nf, modpr);
    if (!aT)
    { /* degree 1 */
      red = ZX_to_Flx(red, pp);
      if (ltp) red = Flx_normalize(red, pp);
      if (!Flx_is_squarefree(red, pp)) { avma = av2; continue; }
      anbf = fl? Flx_nbroots(red, pp): Flx_nbfact(red, pp);
    }
    else
    {
      GEN q;
      if (ltp) red = FqX_normalize(red, aT,ap);
      if (!FqX_is_squarefree(red,aT,ap)) { avma = av2; continue; }
      q = gpowgs(ap, degpol(aT));
      anbf = fl? FqX_split_deg1(&fa, red, q, aT, ap)
               : FqX_split_by_degree(&fa, red, q, aT, ap);
    }
    if (fl == 2 && anbf < dpol) return anbf;
    if (anbf <= 1)
    {
      if (!fl) return anbf; /* irreducible */
      if (!anbf) return 0; /* no root */
    }

    if (!nbf || anbf < nbf
             || (anbf == nbf && cmpii(gel(apr,4), gel(*pr,4)) > 0))
    {
      nbf = anbf;
      *pr = apr;
      *Tp = aT;
      *Fa = fa;
    }
    else avma = av2;
    if (DEBUGLEVEL>3)
      fprintferr("%3ld %s at prime\n  %Z\nTime: %ld\n",
                 anbf, fl?"roots": "factors", apr, TIMER(&ti_pr));
    if (--ct <= 0) return nbf;
  }
}

/* return the factorization of the square-free polynomial x.
   The coeffs of x are in Z_nf and its leading term is a rational integer.
   deg(x) > 1, deg(nfpol) > 1
   If fl = 1, return only the roots of x in nf
   If fl = 2, as fl=1 if pol splits, [] otherwise */
static GEN
nfsqff(GEN nf, GEN pol, long fl)
{
  long n, nbf, dpol = degpol(pol);
  GEN pr, C0, polbase, init_fa = NULL;
  GEN N2, rep, polmod, polred, lt, nfpol = gel(nf,1);
  nfcmbf_t T;
  nflift_t L;
  pari_timer ti, ti_tot;

  if (DEBUGLEVEL>2) { TIMERstart(&ti); TIMERstart(&ti_tot); }
  n = degpol(nfpol);
  polbase = unifpol(nf, pol, t_COL);
  if (typ(polbase) != t_POL) pari_err(typeer, "nfsqff");
  polmod  = lift_intern( unifpol(nf, pol, t_POLMOD) );
  if (dpol == 1) return mkvec(QXQX_normalize(polmod, nfpol));
  /* heuristic */
  if (dpol*3 < n) 
  {
    GEN z, t;
    long i;
    if (DEBUGLEVEL>2) fprintferr("Using Trager's method\n");
    z = (GEN)polfnf(polmod, nfpol)[1];
    if (fl) {
      long l = lg(z);
      for (i = 1; i < l; i++)
      {
        t = gel(z,i); if (degpol(t) > 1) break;
        gel(z,i) = gneg(gdiv(gel(t,2), gel(t,3)));
      }
      setlg(z, i);
      if (fl == 2 && i != l) return cgetg(1,t_VEC);
    }
    return z;
  }

  nbf = nf_pick_prime(5, nf, polbase, fl, &lt, &init_fa, &pr, &L.Tp);
  if (fl == 2 && nbf < dpol) return cgetg(1,t_VEC);
  if (nbf <= 1)
  {
    if (!fl) return mkvec(QXQX_normalize(polmod, nfpol)); /* irreducible */
    if (!nbf) return cgetg(1,t_VEC); /* no root */
  }

  if (DEBUGLEVEL>2) {
    msgTIMER(&ti, "choice of a prime ideal");
    fprintferr("Prime ideal chosen: %Z\n", pr);
  }

  pol = simplify_i(lift(polmod));
  L.tozk = gel(nf,8);
  L.topow= Q_remove_denom(gel(nf,7), &L.topowden);
  T.ZC = L2_bound(nf, L.tozk, &(T.dn));
  T.Br = nf_root_bounds(pol, nf); if (lt) T.Br = gmul(T.Br, lt);

  if (fl) C0 = normlp(T.Br, 2, n);
  else    C0 = nf_factor_bound(nf, polbase); /* bound for T_2(Q_i), Q | P */
  T.bound = mulrr(T.ZC, C0); /* bound for |Q_i|^2 in Z^n on chosen Z-basis */

  N2 = mulsr(dpol*dpol, normlp(T.Br, 4, n)); /* bound for T_2(lt * S_2) */
  T.BS_2 = mulrr(T.ZC, N2); /* bound for |S_2|^2 on chosen Z-basis */

  if (DEBUGLEVEL>2) {
    msgTIMER(&ti, "bound computation");
    fprintferr("  1) T_2 bound for %s: %Z\n", fl?"root":"factor", C0);
    fprintferr("  2) Conversion from T_2 --> | |^2 bound : %Z\n", T.ZC);
    fprintferr("  3) Final bound: %Z\n", T.bound);
  }

  L.p = gel(pr,1);
  if (L.Tp && degpol(L.Tp) == 1) L.Tp = NULL;
  bestlift_init(0, nf, pr, T.bound, &L);
  if (DEBUGLEVEL>2) TIMERstart(&ti);
  polred = ZqX_normalize(polbase, lt, &L); /* monic */

  if (fl) {
    GEN z = nf_DDF_roots(pol, polred, nfpol, lt, init_fa, nbf, fl, &L);
    if (lg(z) == 1) return cgetg(1, t_VEC);
    return z;
  }

  {
    pari_sp av = avma;
    if (L.Tp)
      rep = FqX_split_all(init_fa, L.Tp, L.p);
    else
    {
      long d;
      rep = cgetg(dpol + 1, t_VEC); gel(rep,1) = FpX_red(polred,L.p);
      d = FpX_split_Berlekamp((GEN*)(rep + 1), L.p);
      setlg(rep, d + 1);
    }
    T.fact  = gerepilecopy(av, sort_vecpol(rep, &cmp_pol));
  }
  if (DEBUGLEVEL>2) msgTIMER(&ti, "splitting mod %Z", pr);
  T.pr = pr;
  T.L  = &L;
  T.polbase = polbase;
  T.pol   = pol;
  T.nf    = nf;
  T.hint  = 1; /* useless */

  rep = nf_combine_factors(&T, polred, L.p, L.k, dpol-1);
  if (DEBUGLEVEL>2)
    fprintferr("Total Time: %ld\n===========\n", TIMER(&ti_tot));
  return rep;
}

GEN
nfrootsall_and_pr(GEN nf, GEN pol)
{
  GEN J1,J2, pr, T;
  pari_sp av = avma;
  GEN z = gerepileupto(av, nfsqff(checknf(nf), pol, 2));
  if (lg(z) == 1) return NULL;
  (void)nf_pick_prime(1, nf, unifpol(nf, pol, t_COL), 2,
                      &J1, &J2, &pr, &T);
  return mkvec2(z, pr);
}

/* return the characteristic polynomial of alpha over nf, where alpha
   is an element of the algebra nf[X]/(T) given as a polynomial in X */
GEN
rnfcharpoly(GEN nf, GEN T, GEN alpha, long v)
{
  long vnf, vT, lT;
  pari_sp av = avma;
  GEN p1;

  nf=checknf(nf); vnf = varn(nf[1]);
  if (v<0) v = 0;
  T = fix_relative_pol(nf,T,1);
  if (typ(alpha) == t_POLMOD) alpha = lift_to_pol(alpha);
  lT = lg(T);
  if (typ(alpha) != t_POL || varn(alpha) == vnf)
    return gerepileupto(av, gpowgs(gsub(pol_x[v], alpha), lT - 3));
  vT = varn(T);
  if (varn(alpha) != vT || varncmp(v, vnf)>=0)
    pari_err(talker,"incorrect variables in rnfcharpoly");
  if (lg(alpha) >= lT) alpha = RgX_rem(alpha, T);
  if (lT <= 4)
    return gerepileupto(av, gsub(pol_x[v], alpha));
  p1 = caract2(T, unifpol(nf,alpha, t_POLMOD), v);
  return gerepileupto(av, unifpol(nf, p1, t_POLMOD));
}
