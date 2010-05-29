/* $Id: polarit2.c 8460 2007-03-27 15:53:51Z kb $

Copyright (C) 2000  The PARI group.

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
/**                         (second part)                             **/
/**                                                                   **/
/***********************************************************************/
#include "pari.h"
#include "paripriv.h"

#define addshift(x,y) addshiftpol((x),(y),1)

/* compute Newton sums S_1(P), ... , S_n(P). S_k(P) = sum a_j^k, a_j root of P
 * If N != NULL, assume p-adic roots and compute mod N [assume integer coeffs]
 * If T != NULL, compute mod (T,N) [assume integer coeffs if N != NULL]
 * If y0!= NULL, precomputed i-th powers, i=1..m, m = length(y0).
 * Not memory clean in the latter case */
GEN
polsym_gen(GEN P, GEN y0, long n, GEN T, GEN N)
{
  long dP=degpol(P), i, k, m;
  pari_sp av1, av2;
  GEN s,y,P_lead;

  if (n<0) pari_err(impl,"polsym of a negative n");
  if (typ(P) != t_POL) pari_err(typeer,"polsym");
  if (!signe(P)) pari_err(zeropoler,"polsym");
  y = cgetg(n+2,t_COL);
  if (y0)
  {
    if (typ(y0) != t_COL) pari_err(typeer,"polsym_gen");
    m = lg(y0)-1;
    for (i=1; i<=m; i++) y[i] = y0[i]; /* not memory clean */
  }
  else
  {
    m = 1;
    gel(y,1) = stoi(dP);
  }
  P += 2; /* strip codewords */

  P_lead = gel(P,dP); if (gcmp1(P_lead)) P_lead = NULL;
  if (P_lead)
  {
    if (N) P_lead = Fq_inv(P_lead,T,N);
    else if (T) P_lead = QXQ_inv(P_lead,T);
  }
  for (k=m; k<=n; k++)
  {
    av1 = avma; s = (dP>=k)? gmulsg(k,gel(P,dP-k)): gen_0;
    for (i=1; i<k && i<=dP; i++)
      s = gadd(s, gmul(gel(y,k-i+1),gel(P,dP-i)));
    if (N)
    {
      s = Fq_red(s, T, N);
      if (P_lead) s = Fq_mul(s, P_lead, T, N);
    }
    else if (T)
    {
      s = grem(s, T);
      if (P_lead) s = grem(gmul(s, P_lead), T);
    }
    else
      if (P_lead) s = gdiv(s, P_lead);
    av2 = avma; gel(y,k+1) = gerepile(av1,av2, gneg(s));
  }
  return y;
}

GEN
polsym(GEN x, long n)
{
  return polsym_gen(x, NULL, n, NULL,NULL);
}

/* assume x and y are polynomials in the same variable whose coeffs can be
 * compared (used to sort polynomial factorizations)
 */

static int
polcmp(void *data, GEN x, GEN y)
{
  int (*coeff_cmp)(GEN,GEN)=(int(*)(GEN,GEN))data;
  long i, lx = lg(x), ly = lg(y);
  int fl;
  if (lx > ly) return  1;
  if (lx < ly) return -1;
  for (i=lx-1; i>1; i--)
    if ((fl = coeff_cmp(gel(x,i), gel(y,i)))) return fl;
  return 0;
}

/* sort generic factorisation */
GEN
sort_factor_gen_aux(GEN y, void *data, int (*cmp)(void *,GEN,GEN))
{
  long n, i;
  pari_sp av = avma;
  GEN a,b,A,B,w;
  a = gel(y,1); n = lg(a); A = new_chunk(n);
  b = gel(y,2);            B = new_chunk(n);
  w = gen_sort_aux(a, cmp_IND | cmp_C, data, cmp);
  for (i=1; i<n; i++) { A[i] = a[w[i]]; B[i] = b[w[i]]; }
  for (i=1; i<n; i++) { a[i] = A[i]; b[i] = B[i]; }
  avma = av; return y;
}

/* sort generic factorisation */
GEN
sort_factor_gen(GEN y, int (*cmp)(GEN,GEN))
{
  long n, i;
  pari_sp av = avma;
  GEN a,b,A,B,w;
  a = gel(y,1); n = lg(a); A = new_chunk(n);
  b = gel(y,2);            B = new_chunk(n);
  w = gen_sort(a, cmp_IND | cmp_C, cmp);
  for (i=1; i<n; i++) { A[i] = a[w[i]]; B[i] = b[w[i]]; }
  for (i=1; i<n; i++) { a[i] = A[i]; b[i] = B[i]; }
  avma = av; return y;
}

GEN
sort_factor(GEN y,int (*cmp)(GEN,GEN))
{
  (void)sort_factor_gen_aux(y,(void*)cmp,polcmp);
  return y;
}

GEN
sort_vecpol_gen(GEN a, int (*cmp)(GEN,GEN))
{
  long n, i;
  pari_sp av = avma;
  GEN A,w;
  n = lg(a); A = new_chunk(n);
  w = gen_sort_aux(a, cmp_IND | cmp_C,(void*)cmp, polcmp);
  for (i=1; i<n; i++) A[i] = a[w[i]];
  for (i=1; i<n; i++) a[i] = A[i];
  avma = av; return a;
}

/*In place sort*/
GEN
sort_vecpol(GEN y, int (*cmp)(GEN,GEN))
{
  (void)sort_vecpol_gen(y,cmp);
  return y;
}

/* centered residue x mod p. po2 = shifti(p, -1) or NULL (euclidean residue) */
GEN
centermodii(GEN x, GEN p, GEN po2)
{
  GEN y = remii(x, p);
  switch(signe(y))
  {
    case 0: break;
    case 1: if (po2 && absi_cmp(y,po2) > 0) y = subii(y, p);
      break;
    case -1: if (!po2 || absi_cmp(y,po2) > 0) y = addii(y, p);
      break;
  }
  return y;
}

long
s_centermod(long x, ulong pp, ulong pps2)
{
  long y = x % (long)pp;
  if (y < 0) y += pp;
  return Fl_center(y, pp,pps2);
}

/* for internal use */
GEN
centermod_i(GEN x, GEN p, GEN ps2)
{
  long i, lx;
  pari_sp av;
  GEN y;

  if (!ps2) ps2 = shifti(p,-1);
  switch(typ(x))
  {
    case t_INT: return centermodii(x,p,ps2);

    case t_POL: lx = lg(x);
      y = cgetg(lx,t_POL); y[1] = x[1];
      for (i=2; i<lx; i++)
      {
	av = avma;
	gel(y,i) = gerepileuptoint(av, centermodii(gel(x,i),p,ps2));
      }
      return normalizepol_i(y, lx);

    case t_COL: lx = lg(x);
      y = cgetg(lx,t_COL);
      for (i=1; i<lx; i++) gel(y,i) = centermodii(gel(x,i),p,ps2);
      return y;

    case t_MAT: lx = lg(x);
      y = cgetg(lx,t_MAT);
      for (i=1; i<lx; i++) gel(y,i) = centermod_i(gel(x,i),p,ps2);
      return y;
    
    case t_VECSMALL: lx = lg(x);
    {
      ulong pp = itou(p), pps2 = itou(ps2);
      y = cgetg(lx,t_VECSMALL);
      for (i=1; i<lx; i++) y[i] = s_centermod(x[i], pp, pps2);
      return y;
    }
  }
  return x;
}

GEN
centermod(GEN x, GEN p) { return centermod_i(x,p,NULL); }

/* assume same variables, y normalized, non constant */
static GEN
polidivis(GEN x, GEN y, GEN bound)
{
  long dx, dy, dz, i, j;
  pari_sp av;
  GEN z,p1,y_lead;

  dy=degpol(y);
  dx=degpol(x);
  dz=dx-dy; if (dz<0) return NULL;
  z=cgetg(dz+3,t_POL); z[1] = x[1];
  x += 2; y += 2; z += 2;
  y_lead = gel(y,dy);
  if (gcmp1(y_lead)) y_lead = NULL;

  p1 = gel(x,dx);
  gel(z,dz) = y_lead? diviiexact(p1,y_lead): icopy(p1);
  for (i=dx-1; i>=dy; i--)
  {
    av = avma; p1 = gel(x,i);
    for (j=i-dy+1; j<=i && j<=dz; j++)
      p1 = subii(p1, mulii(gel(z,j),gel(y,i-j)));
    if (y_lead) p1 = diviiexact(p1,y_lead);
    if (bound && absi_cmp(p1, bound) > 0) return NULL;
    p1 = gerepileupto(av,p1);
    gel(z,i-dy) = p1;
  }
  av = avma;
  for (;; i--)
  {
    p1 = gel(x,i);
    /* we always enter this loop at least once */
    for (j=0; j<=i && j<=dz; j++)
      p1 = subii(p1, mulii(gel(z,j),gel(y,i-j)));
    if (!gcmp0(p1)) return NULL;
    avma = av;
    if (!i) break;
  }
  return z - 2;
}

/***********************************************************************/
/**                                                                   **/
/**       QUADRATIC HENSEL LIFT (adapted from V. Shoup's NTL)         **/
/**                                                                   **/
/***********************************************************************/
/* Setup for divide/conquer quadratic Hensel lift
 *  a = set of k t_POL in Z[X] = factors over Fp (T=NULL) or Fp[Y]/(T)
 *  V = set of products of factors built as follows
 *  1) V[1..k] = initial a
 *  2) iterate:
 *      append to V the two smallest factors (minimal degree) in a, remove them
 *      from a and replace them by their product [net loss for a = 1 factor]
 *
 * W = bezout coeffs W[i]V[i] + W[i+1]V[i+1] = 1
 *
 * link[i] = -j if V[i] = a[j]
 *            j if V[i] = V[j] * V[j+1]
 * Arrays (link, V, W) pre-allocated for 2k - 2 elements */
static void
BuildTree(GEN link, GEN V, GEN W, GEN a, GEN T, GEN p)
{
  long k = lg(a)-1;
  long i, j, s, minp, mind;

  for (i=1; i<=k; i++) { V[i] = a[i]; link[i] = -i; }

  for (j=1; j <= 2*k-5; j+=2,i++)
  {
    minp = j;
    mind = degpol(V[j]);
    for (s=j+1; s<i; s++)
      if (degpol(V[s]) < mind) { minp = s; mind = degpol(V[s]); }

    lswap(V[j], V[minp]);
    lswap(link[j], link[minp]);

    minp = j+1;
    mind = degpol(V[j+1]);
    for (s=j+2; s<i; s++)
      if (degpol(V[s]) < mind) { minp = s; mind = degpol(V[s]); }

    lswap(V[j+1], V[minp]);
    lswap(link[j+1], link[minp]);

    if (T)
      gel(V,i) = FpXQX_mul(gel(V,j), gel(V,j+1), T, p);
    else
      gel(V,i) = FpX_mul(gel(V,j), gel(V,j+1), p);
    link[i] = j;
  }

  for (j=1; j <= 2*k-3; j+=2)
  {
    GEN d, u, v;
    if (T)
      d = FpXQX_extgcd(gel(V,j), gel(V,j+1), T, p, &u, &v);
    else
      d = FpX_extgcd(gel(V,j), gel(V,j+1), p, &u, &v);
    if (degpol(d) > 0) pari_err(talker, "relatively prime polynomials expected");
    d = gel(d,2);
    if (!gcmp1(d))
    {
      if (typ(d)==t_POL)
      {
        d = FpXQ_inv(d, T, p);
        u = FqX_Fq_mul(u, d, T, p);
        v = FqX_Fq_mul(v, d, T, p);
      }
      else
      {
        d = Fp_inv(d, p);
        u = FpX_Fp_mul(u, d, p);
        v = FpX_Fp_mul(v, d, p);
      }
    }
    gel(W,j) = u;
    gel(W,j+1) = v;
  }
}

/* au + bv = 1 (p0), ab = f (p0). Lift mod p1 = p0 pd (<= p0^2).
 * If noinv is set, don't lift the inverses u and v */
static void
HenselLift(GEN V, GEN W, long j, GEN f, GEN T, GEN pd, GEN p0, int noinv)
{
  pari_sp av = avma;
  long space = lg(f) * (lgefint(pd) + lgefint(p0));
  GEN a2,b2,g,z,s,t;
  GEN a = gel(V,j), b = gel(V,j+1);
  GEN u = gel(W,j), v = gel(W,j+1);

  if (T) space *= lg(T);

  (void)new_chunk(space); /* HACK */
  g = gadd(f, gneg_i(gmul(a,b)));
  if (T) g = FpXQX_red(g, T, mulii(p0,pd));
  g = gdivexact(g, p0);
  if (T)
  {
    z = FpXQX_mul(v,g, T,pd);
    t = FpXQX_divrem(z,a, T,pd, &s);
  }
  else
  {
    g = FpX_red(g, pd);
    z = FpX_mul(v,g, pd);
    t = FpX_divrem(z,a, pd, &s);
  }
  t = gadd(gmul(u,g), gmul(t,b)); t = T? FpXQX_red(t, T, pd): FpX_red(t, pd);
  t = gmul(t,p0);
  s = gmul(s,p0);
  avma = av;

  /* already reduced mod p1 = pd p0 */
  a2 = gadd(a,s); gel(V,j) = a2;
  b2 = gadd(b,t); gel(V,j+1) = b2;
  if (noinv) return;

  av = avma;
  (void)new_chunk(space); /* HACK */
  g = gadd(gmul(u,a2), gmul(v,b2));
  g = gadd(gneg_i(g), gen_1);

  if (T) g = FpXQX_red(g, T, mulii(p0,pd));
  g = gdivexact(g, p0);
  if (T)
  {
    z = FpXQX_mul(v,g, T,pd);
    t = FpXQX_divrem(z,a, T,pd, &s);
  }
  else
  {
    g = FpX_red(g, pd);
    z = FpX_mul(v,g, pd);
    t = FpX_divrem(z,a, pd, &s);
  }
  t = gadd(gmul(u,g), gmul(t,b)); t = T? FpXQX_red(t, T, pd): FpX_red(t, pd);
  t = gmul(t,p0);
  s = gmul(s,p0);
  avma = av;

  u = gadd(u,t); gel(W,j) = u;
  v = gadd(v,s); gel(W,j+1) = v;
}

/* v list of factors, w list of inverses.  f = v[j] v[j+1]
 * Lift v[j] and v[j+1] mod p0 pd (possibly mod T), then all their divisors */
static void
RecTreeLift(GEN link, GEN v, GEN w, GEN T, GEN pd, GEN p0, GEN f, long j, int noinv)
{
  if (j < 0) return;

  HenselLift(v, w, j, f, T, pd, p0, noinv);
  RecTreeLift(link, v, w, T, pd, p0, gel(v,j)  , link[j  ], noinv);
  RecTreeLift(link, v, w, T, pd, p0, gel(v,j+1), link[j+1], noinv);
}

/* lift from p^{e0} to p^{e1} */
static void
TreeLift(GEN link, GEN v, GEN w, GEN T, GEN p, long e0, long e1, GEN f, int noinv)
{
  GEN p0 = powiu(p, e0);
  GEN pd = powiu(p, e1-e0);
  RecTreeLift(link, v, w, T, pd, p0, f, lgpol(v), noinv);
}

/* Successive accuracies for a quadratic lift.
 * Eg 9 --> 9,5,3,2,1 instead of 9,8,4,2,1 */
GEN
Newton_exponents(long e)
{
  GEN E = cgetg(BITS_IN_LONG, t_VECSMALL);
  long l = 1; E[l++] = e;
  while (e > 1) { e = (e+1)>>1; E[l++] = e; }
  setlg(E, l); return E;
}

/* a = modular factors of f mod (p,T) [possibly T=NULL], lift to precision p^e0
 * flag = 0: standard.
 * flag = 1: return TreeLift structure
 * flag = 2: a = TreeLift structure, go on lifting, as flag = 1 otherwise */
static GEN
MultiLift(GEN f, GEN a, GEN T, GEN p, long e0, long flag)
{
  long l, i, e = e0, k = lg(a) - 1;
  GEN E, v, w, link;
  pari_timer Ti;

  if (k < 2 || e < 1) pari_err(talker, "MultiLift: bad args");
  if (e == 1) return a;
  if (typ(a[1]) == t_INT) flag = 2;
  else if (flag == 2) flag = 1;

  E = Newton_exponents(e);
  e = 1;
  l = lg(E)-1;

  if (DEBUGLEVEL > 3) TIMERstart(&Ti);

  if (flag != 2)
  {
    v = cgetg(2*k - 2 + 1, t_VEC);
    w = cgetg(2*k - 2 + 1, t_VEC);
    link=cgetg(2*k - 2 + 1, t_VECSMALL);
    BuildTree(link, v, w, a, T, p);
    if (DEBUGLEVEL > 3) msgTIMER(&Ti, "building tree");
  }
  else
  {
    e = itos(gel(a,1));
    link = gel(a,2);
    v    = gel(a,3);
    w    = gel(a,4);
  }

  for (i = l; i > 1; i--) {
    if (E[i-1] < e) continue;
    TreeLift(link, v, w, T, p, E[i], E[i-1], f, (flag == 0) && (i == 2));
    if (DEBUGLEVEL > 3) msgTIMER(&Ti, "lifting to prec %ld", E[i-1]);
  }

  if (flag)
    E = mkvec4(stoi(e0), link, v, w);
  else
  {
    E = cgetg(k+1, t_VEC);
    for (i = 1; i <= 2*k-2; i++)
    {
      long t = link[i];
      if (t < 0) E[-t] = v[i];
    }
  }
  return E;
}

/* Q list of (coprime, monic) factors of pol mod (T,p). Lift mod p^e = pe.
 * T may be NULL */
GEN
hensel_lift_fact(GEN pol, GEN Q, GEN T, GEN p, GEN pe, long e)
{
  pari_sp av = avma;
  if (lg(Q) == 2) return mkvec(pol);
  pol = FqX_normalize(pol, T, pe);
  return gerepilecopy(av, MultiLift(pol, Q, T, p, e, 0));
}

/* U = NULL treated as 1 */
static void
BezoutPropagate(GEN link, GEN v, GEN w, GEN pe, GEN U, GEN f, long j)
{
  GEN Q, R;
  if (j < 0) return;

  Q = FpX_mul(gel(v,j), gel(w,j), pe);
  if (U)
  {
    Q = FpXQ_mul(Q, U, f, pe);
    R = FpX_sub(U, Q, pe);
  }
  else
    R = FpX_Fp_add(FpX_neg(Q, pe), gen_1, pe);
  gel(w,j+1) = Q; /*  0 mod U v[j],  1 mod (1-U) v[j+1] */
  gel(w,j) = R; /*  1 mod U v[j],  0 mod (1-U) v[j+1] */
  BezoutPropagate(link, v, w, pe, R, f, link[j  ]);
  BezoutPropagate(link, v, w, pe, Q, f, link[j+1]);
}

/* as above, but return the Bezout coefficients for the lifted modular factors
 *   U[i] = 1 mod Qlift[i]
 *          0 mod Qlift[j], j != i */
GEN
bezout_lift_fact(GEN pol, GEN Q, GEN p, long e)
{
  pari_sp av = avma;
  GEN E, link, v, w, pe;
  long i, k = lg(Q)-1;
  if (k == 1) return mkvec(pol);
  pe = powiu(p, e);
  pol = FpX_normalize(pol, pe);
  E = MultiLift(pol, Q, NULL, p, e, 1);
  link = gel(E,2);
  v    = gel(E,3);
  w    = gel(E,4);
  BezoutPropagate(link, v, w, pe, NULL, pol, lgpol(v));
  E = cgetg(k+1, t_VEC);
  for (i = 1; i <= 2*k-2; i++)
  {
    long t = link[i];
    if (t < 0) E[-t] = w[i];
  }
  return gerepilecopy(av, E);
}

/* Front-end for hensel_lift_fact:
   lift the factorization of pol mod p given by fct to p^exp (if possible) */
GEN
polhensellift(GEN pol, GEN fct, GEN p, long exp)
{
  GEN p1, p2;
  long i, j, l;
  pari_sp av = avma;

  /* we check the arguments */
  if (typ(pol) != t_POL) pari_err(talker, "not a polynomial in polhensellift");
  if ((typ(fct) != t_COL && typ(fct) != t_VEC) || (lg(fct) < 3))
    pari_err(talker, "not a factorization in polhensellift");
  if (typ(p) != t_INT) pari_err(talker, "not a prime number in polhensellift");
  if (exp < 1) pari_err(talker, "not a positive exponent in polhensellift");

  l = lg(pol);
  for (i = 2; i < l; i++)
    if (typ(pol[i]) != t_INT)
      pari_err(talker, "not an integral polynomial in polhensellift");
  p1 = lift(fct); /* make sure the coeffs are integers and not intmods */
  l = lg(p1);
  for (i = 1; i < l; i++)
  {
    p2 = gel(p1,i);
    if (typ(p2) != t_POL)
    {
      if (typ(p2) != t_INT)
        pari_err(talker, "not an integral factorization in polhensellift");
      gel(p1,i) = scalarpol(p2, varn(pol));
    }
  }

  /* then we check that pol \equiv \prod f ; f \in fct mod p */
  p2 = gel(p1,1);
  for (j = 2; j < l; j++) p2 = FpX_mul(p2, gel(p1,j), p);
  if (!gcmp0(FpX_sub(pol, p2, p)))
    pari_err(talker, "not a correct factorization in polhensellift");

  /* finally we check that the elements of fct are coprime mod p */
  if (!FpX_is_squarefree(pol, p))
  {
    for (i = 1; i < l; i++)
      for (j = 1; j < i; j++)
        if (degpol(FpX_gcd(gel(p1,i), gel(p1,j), p)))
          pari_err(talker, "polhensellift: factors %Z and %Z are not coprime",
                     p1[i], p1[j]);
  }
  return gerepilecopy(av, hensel_lift_fact(pol,p1,NULL,p,powiu(p,exp),exp));
}

#if 0
/* cf Beauzamy et al: upper bound for
 *      lc(x) * [2^(5/8) / pi^(3/8)] e^(1/4n) 2^(n/2) sqrt([x]_2)/ n^(3/8)
 * where [x]_2 = sqrt(\sum_i=0^n x[i]^2 / binomial(n,i)). One factor has
 * all coeffs less than then bound */
static GEN
two_factor_bound(GEN x)
{
  long i, j, n = lg(x) - 3;
  pari_sp av = avma;
  GEN *invbin, c, r = cgetr(3), z;

  x += 2; invbin = (GEN*)new_chunk(n+1);
  z = real_1(3); /* invbin[i] = 1 / binomial(n, i) */
  for (i=0,j=n; j >= i; i++,j--)
  {
    invbin[i] = invbin[j] = z;
    z = divrs(mulrs(z, i+1), n-i);
  }
  z = invbin[0]; /* = 1 */
  for (i=0; i<=n; i++)
  {
    c = gel(x,i); if (!signe(c)) continue;
    affir(c, r);
    z = addrr(z, mulrr(gsqr(r), invbin[i]));
  }
  z = shiftr(sqrtr(z), n);
  z = divrr(z, dbltor(pow((double)n, 0.75)));
  z = grndtoi(sqrtr(z), &i);
  z = mulii(z, absi(gel(x,n)));
  return gerepileuptoint(av, shifti(z, 1));
}
#endif

/* A | S ==> |a_i| <= binom(d-1, i-1) || S ||_2 + binom(d-1, i) lc(S) */
static GEN
Mignotte_bound(GEN S)
{
  long i, d = degpol(S);
  GEN lS, C, binlS, bin, N2, p1;
  
  N2 = sqrtr(QuickNormL2(S,DEFAULTPREC));
  binlS = bin = vecbinome(d-1);
  lS = leading_term(S);
  if (!is_pm1(lS)) binlS = gmul(lS, bin);

  /* i = 0 */
  C = gel(binlS,1);
  /* i = d */
  p1 = N2; if (gcmp(C, p1) < 0) C = p1;
  for (i = 1; i < d; i++)
  {
    p1 = gadd(gmul(gel(bin,i), N2), gel(binlS,i+1));
    if (gcmp(C, p1) < 0) C = p1;
  }
  return C;
}
/* A | S ==> |a_i|^2 <= 3^{3/2 + d} / (4 \pi d) [P]_2^2,
 * where [P]_2 is Bombieri's 2-norm */
static GEN
Beauzamy_bound(GEN S)
{
  const long prec = DEFAULTPREC;
  long i, d = degpol(S);
  GEN bin, lS, s, C;
  bin = vecbinome(d);

  s = real_0(prec);
  for (i=0; i<=d; i++)
  {
    GEN c = gel(S,i+2);
    if (gcmp0(c)) continue;
    /* s += P_i^2 / binomial(d,i) */
    s = addrr(s, gdiv(itor(sqri(c), prec), gel(bin,i+1)));
  }
  /* s = [S]_2^2 */
  C = powrshalf(stor(3,prec), 3 + 2*d); /* 3^{3/2 + d} */
  C = gdiv(gmul(C, s), gmulsg(4*d, mppi(prec)));
  lS = absi(leading_term(S));
  return mulir(lS, sqrtr(C));
}

static GEN
factor_bound(GEN S)
{
  pari_sp av = avma;
  GEN a = Mignotte_bound(S);
  GEN b = Beauzamy_bound(S);
  if (DEBUGLEVEL>2)
  {
    fprintferr("Mignotte bound: %Z\n",a);
    fprintferr("Beauzamy bound: %Z\n",b);
  }
  return gerepileupto(av, ceil_safe(gmin(a, b)));
}

#if 0
/* all factors have coeffs less than the bound */
static GEN
all_factor_bound(GEN x)
{
  long i, n = degpol(x);
  GEN t, z = gen_0;
  for (i=2; i<=n+2; i++) z = addii(z, sqri(gel(x,i)));
  t = absi(gel(x,n+2));
  z = addii(t, addsi(1, sqrti(z)));
  z = mulii(z, binomial(stoi(n-1), n>>1));
  return shifti(mulii(t,z),1);
}
#endif

typedef ulong *uGEN;

/* Naive recombination of modular factors: combine up to maxK modular
 * factors, degree <= klim and divisible by hint
 *
 * target = polynomial we want to factor
 * famod = array of modular factors.  Product should be congruent to
 * target/lc(target) modulo p^a
 * For true factors: S1,S2 <= p^b, with b <= a and p^(b-a) < 2^31 */
static GEN
cmbf(GEN pol, GEN famod, GEN bound, GEN p, long a, long b,
     long maxK, long klim,long hint)
{
  long K = 1, cnt = 1, i,j,k, curdeg, lfamod = lg(famod)-1;
  ulong spa_b, spa_bs2, Sbound;
  GEN lc, lcpol, pa = powiu(p,a), pas2 = shifti(pa,-1);
  uGEN trace1   = (uGEN)cgetg(lfamod+1, t_VECSMALL);
  uGEN trace2   = (uGEN)cgetg(lfamod+1, t_VECSMALL);
  GEN ind      = cgetg(lfamod+1, t_VECSMALL);
  GEN degpol   = cgetg(lfamod+1, t_VECSMALL);
  GEN degsofar = cgetg(lfamod+1, t_VECSMALL);
  GEN listmod  = cgetg(lfamod+1, t_COL);
  GEN fa       = cgetg(lfamod+1, t_COL);

  if (maxK < 0) maxK = lfamod-1;

  lc = absi(leading_term(pol));
  if (is_pm1(lc)) lc = NULL;
  lcpol = lc? gmul(lc,pol): pol;

  {
    GEN pa_b,pa_bs2,pb, lc2 = lc? sqri(lc): NULL;

    pa_b = powiu(p, a-b); /* < 2^31 */
    pa_bs2 = shifti(pa_b,-1);
    pb= powiu(p, b);
    for (i=1; i <= lfamod; i++)
    {
      GEN T1,T2, P = gel(famod,i);
      long d = degpol(P);

      degpol[i] = d; P += 2;
      T1 = gel(P,d-1);/* = - S_1 */
      T2 = sqri(T1);
      if (d > 1) T2 = subii(T2, shifti(gel(P,d-2),1));
      T2 = modii(T2, pa); /* = S_2 Newton sum */
      if (lc)
      {
        T1 = modii(mulii(lc, T1), pa);
        T2 = modii(mulii(lc2,T2), pa);
      }
      trace1[i] = itou(diviiround(T1, pb));
      trace2[i] = itou(diviiround(T2, pb));
    }
    spa_b   = (ulong)  pa_b[2]; /* < 2^31 */
    spa_bs2 = (ulong)pa_bs2[2]; /* < 2^31 */
  }
  degsofar[0] = 0; /* sentinel */

  /* ind runs through strictly increasing sequences of length K,
   * 1 <= ind[i] <= lfamod */
nextK:
  if (K > maxK || 2*K > lfamod) goto END;
  if (DEBUGLEVEL > 3)
    fprintferr("\n### K = %d, %Z combinations\n", K,binomial(utoipos(lfamod), K));
  setlg(ind, K+1); ind[1] = 1;
  Sbound = (ulong) ((K+1)>>1);
  i = 1; curdeg = degpol[ind[1]];
  for(;;)
  { /* try all combinations of K factors */
    for (j = i; j < K; j++)
    {
      degsofar[j] = curdeg;
      ind[j+1] = ind[j]+1; curdeg += degpol[ind[j+1]];
    }
    if (curdeg <= klim && curdeg % hint == 0) /* trial divide */
    {
      GEN y, q, list;
      pari_sp av;
      ulong t;

      /* d - 1 test */
      for (t=trace1[ind[1]],i=2; i<=K; i++)
        t = Fl_add(t, trace1[ind[i]], spa_b);
      if (t > spa_bs2) t = spa_b - t;
      if (t > Sbound)
      {
        if (DEBUGLEVEL>6) fprintferr(".");
        goto NEXT;
      }
      /* d - 2 test */
      for (t=trace2[ind[1]],i=2; i<=K; i++)
        t = Fl_add(t, trace2[ind[i]], spa_b);
      if (t > spa_bs2) t = spa_b - t;
      if (t > Sbound)
      {
        if (DEBUGLEVEL>6) fprintferr("|");
        goto NEXT;
      }

      av = avma;
      /* check trailing coeff */
      y = lc;
      for (i=1; i<=K; i++)
      {
        GEN q = constant_term(gel(famod,ind[i]));
        if (y) q = mulii(y, q);
        y = centermod_i(q, pa, pas2);
      }
      if (!signe(y) || remii(constant_term(lcpol), y) != gen_0)
      {
        if (DEBUGLEVEL>3) fprintferr("T");
        avma = av; goto NEXT;
      }
      y = lc; /* full computation */
      for (i=1; i<=K; i++)
      {
        GEN q = gel(famod,ind[i]);
        if (y) q = gmul(y, q);
        y = centermod_i(q, pa, pas2);
      }

      /* y is the candidate factor */
      if (! (q = polidivis(lcpol,y,bound)) )
      {
        if (DEBUGLEVEL>3) fprintferr("*");
        avma = av; goto NEXT;
      }
      /* found a factor */
      list = cgetg(K+1, t_VEC);
      gel(listmod,cnt) = list;
      for (i=1; i<=K; i++) list[i] = famod[ind[i]];

      y = Q_primpart(y);
      gel(fa,cnt++) = y;
      /* fix up pol */
      pol = q;
      if (lc) pol = Q_div_to_int(pol, leading_term(y));
      for (i=j=k=1; i <= lfamod; i++)
      { /* remove used factors */
        if (j <= K && i == ind[j]) j++;
        else
        {
          famod[k] = famod[i];
          trace1[k] = trace1[i];
          trace2[k] = trace2[i];
          degpol[k] = degpol[i]; k++;
        }
      }
      lfamod -= K;
      if (lfamod < 2*K) goto END;
      i = 1; curdeg = degpol[ind[1]];
      bound = factor_bound(pol);
      if (lc) lc = absi(leading_term(pol));
      lcpol = lc? gmul(lc,pol): pol;
      if (DEBUGLEVEL>3)
        fprintferr("\nfound factor %Z\nremaining modular factor(s): %ld\n",
                   y, lfamod);
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
    if (signe(leading_term(pol)) < 0) pol = gneg_i(pol);

    setlg(famod, lfamod+1);
    gel(listmod,cnt) = shallowcopy(famod);
    gel(fa,cnt++) = pol;
  }
  if (DEBUGLEVEL>6) fprintferr("\n");
  setlg(listmod, cnt);
  setlg(fa, cnt); return mkvec2(fa, listmod);
}

void
factor_quad(GEN x, GEN res, long *ptcnt)
{
  GEN a = gel(x,4), b = gel(x,3), c = gel(x,2), d, u, z1, z2, t;
  GEN D = subii(sqri(b), shifti(mulii(a,c), 2));
  long v, cnt = *ptcnt;

  if (!Z_issquarerem(D, &d)) { gel(res,cnt++) = x; *ptcnt = cnt; return; }

  t = shifti(negi(addii(b, d)), -1);
  z1 = gdiv(t, a); u = denom(z1);
  z2 = gdiv(addii(t, d), a);
  v = varn(x);
  gel(res,cnt++) = gmul(u, gsub(pol_x[v], z1)); u = diviiexact(a, u);
  gel(res,cnt++) = gmul(u, gsub(pol_x[v], z2)); *ptcnt = cnt;
}

/* y > 1 and B integers. Let n such that y^(n-1) <= B < y^n.
 * Return e = max(n,1), set *ptq = y^e if non-NULL */
long
logint(GEN B, GEN y, GEN *ptq)
{
  pari_sp av = avma;
  long e,i,fl;
  GEN q,pow2, r = y;

  if (typ(B) != t_INT) B = ceil_safe(B);
  if (expi(B) <= (expi(y) << 6)) /* e small, be naive */
  {
    for (e=1; cmpii(r, B) < 0; e++) r = mulii(r,y);
    goto END;
  }
  /* binary splitting: compute bits of e one by one */
  /* compute pow2[i] = y^(2^i) [i < very crude upper bound for log_2(n)] */
  pow2 = new_chunk(bit_accuracy(lgefint(B)));
  gel(pow2,0) = y;
  for (i=0,q=r;; )
  {
    fl = cmpii(r,B); if (fl >= 0) break;
    q = r; r = sqri(q);
    i++; gel(pow2,i) = r;
  }
  if (i == 0) { e = 1; goto END; } /* y <= B */

  for (i--, e=1<<i;;)
  { /* y^e = q < B <= r = q * y^(2^i) */
    if (!fl) break; /* B = r */
    /* q < B < r */
    if (--i < 0) { if (fl > 0) e++; break; }
    r = mulii(q, gel(pow2,i));
    fl = cmpii(r, B);
    if (fl <= 0) { e += (1<<i); q = r; }
  }
  if (fl <= 0) { e++; r = mulii(r,y); }
END:
  if (ptq) *ptq = gerepileuptoint(av, icopy(r)); else avma = av;
  return e;
}

/* recombination of modular factors: van Hoeij's algorithm */

/* Q in Z[X], return Q(2^n) */
static GEN
shifteval(GEN Q, long n)
{
  long i, l = lg(Q);
  GEN s;

  if (!signe(Q)) return gen_0;
  s = gel(Q,l-1);
  for (i = l-2; i > 1; i--) s = addii(gel(Q,i), shifti(s, n));
  return s;
}

/* return integer y such that all |a| <= y if P(a) = 0 */
static GEN
root_bound(GEN P0)
{
  GEN Q = shallowcopy(P0), lP = absi(leading_term(Q)), x,y,z;
  long k, d = degpol(Q);

  /* P0 = lP x^d + Q, deg Q < d */
  Q = normalizepol_i(Q, d+2);
  for (k=lg(Q)-1; k>1; k--) gel(Q,k) = absi(gel(Q,k));
  k = (long)(cauchy_bound(P0) / LOG2);
  for (  ; k >= 0; k--)
  {
    pari_sp av = avma;
    /* y = 2^k; Q(y) >= lP y^d ? */
    if (cmpii(shifteval(Q,k), shifti(lP, d*k)) >= 0) break;
    avma = av;
  }
  if (k < 0) k = 0;
  x = int2n(k);
  y = int2n(k+1);
  for(k=0; ; k++)
  {
    z = shifti(addii(x,y), -1);
    if (equalii(x,z) || k > 5) break;
    if (cmpii(poleval(Q,z), mulii(lP, powiu(z, d))) < 0)
      y = z;
    else
      x = z;
  }
  return y;
}

static GEN
ZM_HNFimage(GEN x)
{
  return (lg(x) > 50)? hnflll_i(x, NULL, 1): hnfall_i(x, NULL, 1);
}

GEN
special_pivot(GEN x)
{
  GEN t, H = ZM_HNFimage(x);
  long i,j, l = lg(H), h = lg(H[1]);
  for (i=1; i<h; i++)
  {
    int fl = 0;
    for (j=1; j<l; j++)
    {
      t = gcoeff(H,i,j);
      if (signe(t))
      {
        if (!is_pm1(t) || fl) return NULL;
        fl = 1;
      }
    }
  }
  return H;
}

/* B from lllint_i: return [ |b_i^*|^2, i ] */
GEN
GS_norms(GEN B, long prec)
{
  long i, l = lg(B);
  GEN v = gmul(B, real_1(prec));
  l--; setlg(v, l);
  for (i=1; i<l; i++)
    gel(v,i) = divrr(gel(v,i+1), gel(v,i));
  return v;
}

GEN
chk_factors_get(GEN lt, GEN famod, GEN c, GEN T, GEN N)
{
  long i = 1, j, l = lg(famod);
  GEN V = cgetg(l, t_VEC);
  for (j = 1; j < l; j++)
    if (signe(c[j])) V[i++] = famod[j];
  if (lt && i > 1) gel(V,1) = gmul(lt, gel(V,1));
  setlg(V, i); 
  if (T) return FpXQXV_prod(V, T, N);
  else return FpXV_prod(V,N);
}

static GEN
chk_factors(GEN P, GEN M_L, GEN bound, GEN famod, GEN pa)
{
  long i, r;
  GEN pol = P, list, piv, y, ltpol, lt;

  piv = special_pivot(M_L);
  if (!piv) return NULL;
  if (DEBUGLEVEL>3) fprintferr("special_pivot output:\n%Z\n",piv);

  r  = lg(piv)-1;
  list = cgetg(r+1, t_COL);
  lt = absi(leading_term(pol));
  if (is_pm1(lt)) lt = NULL;
  ltpol = lt? gmul(lt, pol): pol;
  for (i = 1;;)
  {
    if (DEBUGLEVEL) fprintferr("LLL_cmbf: checking factor %ld\n",i);
    y = chk_factors_get(lt, famod, gel(piv,i), NULL, pa);
    y = FpX_center(y, pa);
    if (! (pol = polidivis(ltpol,y,bound)) ) return NULL;
    if (lt) y = Q_primpart(y);
    gel(list,i) = y;
    if (++i >= r) break;

    if (lt)
    {
      pol = gdivexact(pol, leading_term(y));
      lt = absi(leading_term(pol));
      ltpol = gmul(lt, pol);
    }
    else
      ltpol = pol;
  }
  y = Q_primpart(pol);
  gel(list,i) = y; return list;
}

GEN
LLL_check_progress(GEN Bnorm, long n0, GEN m, int final, long *ti_LLL)
{
  GEN B, norm, u;
  long i, R;
  pari_timer T;

  if (DEBUGLEVEL>2) TIMERstart(&T);
  u = lllint_i(m, final? 1000: 4, 0, NULL, NULL, &B);
  if (DEBUGLEVEL>2) *ti_LLL += TIMER(&T);
  norm = GS_norms(B, DEFAULTPREC);
  for (R=lg(m)-1; R > 0; R--)
    if (cmprr(gel(norm,R), Bnorm) < 0) break;
  for (i=1; i<=R; i++) setlg(u[i], n0+1);
  if (R <= 1)
  {
    if (!R) pari_err(bugparier,"LLL_cmbf [no factor]");
    return NULL; /* irreducible */
  }
  setlg(u, R+1); return u;
}

static ulong
next2pow(ulong a)
{
  ulong b = 1;
  while (b < a) b <<= 1;
  return b;
}

/* Recombination phase of Berlekamp-Zassenhaus algorithm using a variant of
 * van Hoeij's knapsack
 *
 * P = squarefree in Z[X].
 * famod = array of (lifted) modular factors mod p^a
 * bound = Mignotte bound for the size of divisors of P (for the sup norm)
 * previously recombined all set of factors with less than rec elts */
static GEN
LLL_cmbf(GEN P, GEN famod, GEN p, GEN pa, GEN bound, long a, long rec)
{
  const long N0 = 1; /* # of traces added at each step */
  double BitPerFactor = 0.5; /* nb bits in p^(a-b) / modular factor */
  long i,j,tmax,n0,C, dP = degpol(P);
  double logp = log((double)itos(p)), LOGp2 = LOG2/logp;
  double b0 = log((double)dP*2) / logp, logBr;
  GEN lP, Br, Bnorm, Tra, T2, TT, CM_L, m, list, ZERO;
  pari_sp av, av2, lim;
  long ti_LLL = 0, ti_CF  = 0;

  lP = absi(leading_term(P));
  if (is_pm1(lP)) lP = NULL;
  Br = root_bound(P);
  if (lP) Br = gmul(lP, Br);
  logBr = gtodouble(glog(Br, DEFAULTPREC)) / logp;

  n0 = lg(famod) - 1;
  C = (long)ceil( sqrt(N0 * n0 / 4.) ); /* > 1 */
  Bnorm = dbltor(n0 * (C*C + N0*n0/4.) * 1.00001);
  ZERO = zeromat(n0, N0);

  av = avma; lim = stack_lim(av, 1);
  TT = cgetg(n0+1, t_VEC);
  Tra  = cgetg(n0+1, t_MAT);
  for (i=1; i<=n0; i++)
  {
    TT[i]  = 0;
    gel(Tra,i) = cgetg(N0+1, t_COL);
  }
  CM_L = gscalsmat(C, n0);
  /* tmax = current number of traces used (and computed so far) */
  for (tmax = 0;; tmax += N0)
  {
    long b, bmin, bgood, delta, tnew = tmax + N0, r = lg(CM_L)-1;
    GEN M_L, q, CM_Lp, oldCM_L;
    int first = 1;
    pari_timer ti2, TI;
    
    bmin = (long)ceil(b0 + tnew*logBr);
    if (DEBUGLEVEL>2)
      fprintferr("\nLLL_cmbf: %ld potential factors (tmax = %ld, bmin = %ld)\n",
                 r, tmax, bmin);

    /* compute Newton sums (possibly relifting first) */
    if (a <= bmin)
    {
      a = (long)ceil(bmin + 3*N0*logBr) + 1; /* enough for 3 more rounds */
      a = (long)next2pow((ulong)a);

      pa = powiu(p,a);
      famod = hensel_lift_fact(P,famod,NULL,p,pa,a);
      for (i=1; i<=n0; i++) TT[i] = 0;
    }
    for (i=1; i<=n0; i++)
    {
      GEN p1 = gel(Tra,i);
      GEN p2 = polsym_gen(gel(famod,i), gel(TT,i), tnew, NULL, pa);
      gel(TT,i) = p2;
      p2 += 1+tmax; /* ignore traces number 0...tmax */
      for (j=1; j<=N0; j++) p1[j] = p2[j];
      if (lP)
      { /* make Newton sums integral */
        GEN lPpow = powiu(lP, tmax);
        for (j=1; j<=N0; j++)
        {
          lPpow = mulii(lPpow,lP);
          gel(p1,j) = mulii(gel(p1,j), lPpow);
        }
      }
    }

    /* compute truncation parameter */
    if (DEBUGLEVEL>2) { TIMERstart(&ti2); TIMERstart(&TI); }
    oldCM_L = CM_L;
    av2 = avma;
    delta = b = 0; /* -Wall */
AGAIN:
    M_L = Q_div_to_int(CM_L, utoipos(C));
    T2 = centermod( gmul(Tra, M_L), pa );
    if (first)
    { /* initialize lattice, using few p-adic digits for traces */
      double t = gexpo(T2) - max(32, BitPerFactor*r);
      bgood = (long) (t * LOGp2);
      b = max(bmin, bgood);
      delta = a - b;
    }
    else
    { /* add more p-adic digits and continue reduction */
      long b0 = (long)(gexpo(T2) * LOGp2);
      if (b0 < b) b = b0;
      b = max(b-delta, bmin);
      if (b - delta/2 < bmin) b = bmin; /* near there. Go all the way */
    }

    q = powiu(p, b);
    m = vconcat( CM_L, gdivround(T2, q) );
    if (first)
    {
      GEN P1 = gscalmat(powiu(p, a-b), N0);
      first = 0;
      m = shallowconcat( m, vconcat(ZERO, P1) );
      /*     [ C M_L        0     ]
       * m = [                    ]   square matrix
       *     [  T2'  p^(a-b) I_N0 ]   T2' = Tra * M_L  truncated
       */
    }

    CM_L = LLL_check_progress(Bnorm, n0, m, b == bmin, /*dbg:*/ &ti_LLL);
    if (DEBUGLEVEL>2)
      fprintferr("LLL_cmbf: (a,b) =%4ld,%4ld; r =%3ld -->%3ld, time = %ld\n",
                 a,b, lg(m)-1, CM_L? lg(CM_L)-1: 1, TIMER(&TI));
    if (!CM_L) { list = mkcol(P); break; }
    if (b > bmin) 
    {
      CM_L = gerepilecopy(av2, CM_L);
      goto AGAIN;
    }
    if (DEBUGLEVEL>2) msgTIMER(&ti2, "for this block of traces");

    i = lg(CM_L) - 1;
    if (i == r && gequal(CM_L, oldCM_L))
    {
      CM_L = oldCM_L;
      avma = av2; continue;
    }

    CM_Lp = FpM_image(CM_L, utoipos(27449)); /* inexpensive test */
    if (lg(CM_Lp) != lg(CM_L))
    {
      if (DEBUGLEVEL>2) fprintferr("LLL_cmbf: rank decrease\n");
      CM_L = ZM_HNFimage(CM_L);
    }

    if (i <= r && i*rec < n0)
    {
      pari_timer ti;
      if (DEBUGLEVEL>2) TIMERstart(&ti);
      list = chk_factors(P, Q_div_to_int(CM_L,utoipos(C)), bound, famod, pa);
      if (DEBUGLEVEL>2) ti_CF += TIMER(&ti);
      if (list) break;
      if (DEBUGLEVEL>2) fprintferr("LLL_cmbf: chk_factors failed");
      CM_L = gerepilecopy(av2, CM_L);
    }
    if (low_stack(lim, stack_lim(av,1)))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"LLL_cmbf");
      gerepileall(av, 5, &CM_L, &TT, &Tra, &famod, &pa);
    }
  }
  if (DEBUGLEVEL>2)
    fprintferr("* Time LLL: %ld\n* Time Check Factor: %ld\n",ti_LLL,ti_CF);
  return list;
}

/* Find a,b minimal such that A < q^a, B < q^b, 1 << q^(a-b) < 2^31 */
int
cmbf_precs(GEN q, GEN A, GEN B, long *pta, long *ptb, GEN *qa, GEN *qb)
{
  long a,b,amin,d = (long)(31 * LOG2/gtodouble(glog(q,DEFAULTPREC)) - 1e-5);
  int fl = 0;

  b = logint(B, q, qb);
  amin = b + d;
  if (gcmp(powiu(q, amin), A) <= 0)
  {
    a = logint(A, q, qa);
    b = a - d; *qb = powiu(q, b);
  }
  else
  { /* not enough room */
    a = amin;  *qa = powiu(q, a);
    fl = 1;
  }
  if (DEBUGLEVEL > 3) {
    fprintferr("S_2   bound: %Z^%ld\n", q,b);
    fprintferr("coeff bound: %Z^%ld\n", q,a);
  }
  *pta = a;
  *ptb = b; return fl;
}

/* use van Hoeij's knapsack algorithm */
GEN
combine_factors(GEN target, GEN famod, GEN p, long klim, long hint)
{
  GEN la, B, A, res, L, pa, pb, listmod;
  long a,b, l, maxK = 3, nft = lg(famod)-1, n = degpol(target);
  pari_timer T;

  A = factor_bound(target);

  la = absi(leading_term(target));
  B = mulsi(n, sqri(gmul(la, root_bound(target)))); /* = bound for S_2 */

  (void)cmbf_precs(p, A, B, &a, &b, &pa, &pb);

  if (DEBUGLEVEL>2) (void)TIMER(&T);
  famod = hensel_lift_fact(target,famod,NULL,p,pa,a);
  if (nft < 11) maxK = -1; /* few modular factors: try all posibilities */
  if (DEBUGLEVEL>2) msgTIMER(&T, "Hensel lift (mod %Z^%ld)", p,a);
  L = cmbf(target, famod, A, p, a, b, maxK, klim, hint);
  if (DEBUGLEVEL>2) msgTIMER(&T, "Naive recombination");

  res     = gel(L,1);
  listmod = gel(L,2); l = lg(listmod)-1;
  famod = gel(listmod,l);
  if (maxK >= 0 && lg(famod)-1 > 2*maxK)
  {
    if (l!=1) A = factor_bound(gel(res,l));
    if (DEBUGLEVEL > 4) fprintferr("last factor still to be checked\n");
    L = LLL_cmbf(gel(res,l), famod, p, pa, A, a, maxK);
    if (DEBUGLEVEL>2) msgTIMER(&T,"Knapsack");
    /* remove last elt, possibly unfactored. Add all new ones. */
    setlg(res, l); res = shallowconcat(res, L);
  }
  return res;
}

/* assume pol(0) != 0, polp = pol/lc(pol) mod p.
 * Return vector of rational roots of a */
GEN
DDF_roots(GEN pol, GEN polp, GEN p)
{
  GEN lc, lcpol, z, pe, pes2, bound;
  long i, m, e, lz, v = varn(pol);
  pari_sp av, lim;
  pari_timer T;

  if (DEBUGLEVEL>2) TIMERstart(&T);
  lc = absi(leading_term(pol));
  if (is_pm1(lc)) lc = NULL;
  lcpol = lc? gmul(lc,pol): pol;

  bound = root_bound(pol);
  if (lc) bound = mulii(lc, bound);
  e = logint(addis(shifti(bound, 1), 1), p, &pe);
  pes2 = shifti(pe, -1);
  if (DEBUGLEVEL>2) msgTIMER(&T, "Root bound");

  av = avma; lim = stack_lim(av,2);
  z = FpX_roots(polp, p);
  lz = lg(z)-1;
  if (lz > (degpol(pol) >> 2))
  { /* many roots */
    z = shallowconcat(deg1_from_roots(z, v),
                 FpX_div(polp, FpV_roots_to_pol(z, p, v), p));
    z = hensel_lift_fact(pol, z, NULL, p, pe, e);
  }
  else
  {
    z = ZpX_liftroots(pol, z, p, e);
    z = deg1_from_roots(z, v);
  }
  if (DEBUGLEVEL>2) msgTIMER(&T, "Hensel lift (mod %Z^%ld)", p,e);

  for (m=1, i=1; i <= lz; i++)
  {
    GEN q, r, y = gel(z,i);
    if (lc) y = gmul(y, lc);
    y = centermod_i(y, pe, pes2);
    if (! (q = polidivis(lcpol, y, NULL)) ) continue;

    lcpol = pol = q;
    r = negi( constant_term(y) );
    if (lc) {
      r = gdiv(r,lc);
      pol = Q_primpart(pol);
      lc = absi( leading_term(pol) );
      if (is_pm1(lc)) lc = NULL; else lcpol = gmul(lc, pol);
    }
    gel(z,m++) = r;
    if (low_stack(lim, stack_lim(av,2)))
    {
      if (DEBUGMEM>1) pari_warn(warnmem,"DDF_roots, m = %ld", m);
      gerepileall(av, lc? 4:2, &z, &pol, &lc, &lcpol);
    
    }
  }
  if (DEBUGLEVEL>2) msgTIMER(&T, "Recombination");
  z[0] = evaltyp(t_VEC) | evallg(m); return z;
}

/* Assume a squarefree, degree(a) > 0, a(0) != 0.
 * If fl != 0 look only for rational roots */
static GEN
DDF(GEN a, long hint, int fl)
{
  GEN lead, prime, famod, z, ap;
  const long da = degpol(a);
  long chosenp, p, nfacp, np, nmax, ti = 0;
  pari_sp av = avma, av1;
  byteptr pt = diffptr;
  const long MAXNP = 7;
  pari_timer T, T2;

  if (DEBUGLEVEL>2) { TIMERstart(&T); TIMERstart(&T2); }
  if (hint <= 0) hint = 1;
  nmax = da+1;
  chosenp = 0;
  lead = gel(a,da+2); if (gcmp1(lead)) lead = NULL;
  av1 = avma;
  for (p = np = 0; np < MAXNP; avma = av1)
  {
    NEXT_PRIME_VIADIFF_CHECK(p,pt);
    if (lead && !smodis(lead,p)) continue;
    z = ZX_to_Flx(a, p);
    if (!Flx_is_squarefree(z, p)) continue;

    nfacp = fl? Flx_nbroots(z, p): Flx_nbfact(z, p);
    if (DEBUGLEVEL>4)
      fprintferr("...tried prime %3ld (%-3ld %s). Time = %ld\n",
                  p, nfacp, fl?"roots": "factors", TIMER(&T2));
    if (nfacp < nmax)
    {
      if (nfacp <= 1)
      {
        if (!fl) { avma = av; return mkcol(a); } /* irreducible */
        if (!nfacp) return cgetg(1, t_VEC); /* no root */
      }
      nmax = nfacp; chosenp = p;
      if (da > 100 && nmax < 5) break; /* large degree, few factors. Enough */
    }
    np++;
  }
  prime = utoipos(chosenp);
  ap = lead? FpX_normalize(a, prime): FpX_red(a, prime);
  if (fl) return gerepilecopy(av, DDF_roots(a, ap, prime));

  famod = cgetg(nmax+1,t_COL);
  gel(famod,1) = ap;
  if (nmax != FpX_split_Berlekamp((GEN*)(famod+1), prime))
    pari_err(bugparier,"DDF: wrong numbers of factors");
  if (DEBUGLEVEL>2)
  {
    if (DEBUGLEVEL>4) msgTIMER(&T2, "splitting mod p = %ld", chosenp);
    ti = TIMER(&T);
    fprintferr("Time setup: %ld\n", ti);
  }
  z = combine_factors(a, famod, prime, da-1, hint);
  if (DEBUGLEVEL>2)
    fprintferr("Total Time: %ld\n===========\n", ti + TIMER(&T));
  return gerepilecopy(av, z);
}

/* A(X^d) --> A(X) */
GEN
poldeflate_i(GEN x0, long d)
{
  GEN z, y, x;
  long i,id, dy, dx = degpol(x0);
  if (d <= 1) return x0;
  if (dx < 0) return zeropol(varn(x0));
  dy = dx/d;
  y = cgetg(dy+3, t_POL); y[1] = x0[1];
  z = y + 2;
  x = x0+ 2;
  for (i=id=0; i<=dy; i++,id+=d) z[i] = x[id];
  return y;
}

long
checkdeflate(GEN x)
{
  long d = 0, i, lx = lg(x);
  for (i=3; i<lx; i++)
    if (!gcmp0(gel(x,i))) { d = cgcd(d,i-2); if (d == 1) break; }
  return d;
}

GEN
gdeflate(GEN x, long v, long d)
{
  long i, lx, tx = typ(x);
  GEN z;
  if (is_scalar_t(tx)) return gcopy(x);
  if (d <= 0) pari_err(talker,"need positive degree in gdeflate");
  if (tx == t_POL || tx == t_SER)
  {
    long vx = varn(x);
    pari_sp av;
    if (varncmp(vx, v) < 0)
    {
      lx = lg(x); z = cgetg(lx, tx); z[1] = x[1];
      for (i=2; i<lx; i++) gel(z,i) = gdeflate(gel(x,i),v,d);
      return z;
    }
    if (varncmp(vx, v) > 0) return gcopy(x);
    av = avma;
    if (tx == t_SER)
    {
      long V = valp(x);
      GEN y;

      lx = lg(x);
      if (lx == 2) return zeroser(v, V / d);
      y = ser2pol_i(x, lx);
      if (V % d != 0 || checkdeflate(y) % d != 0)
        pari_err(talker, "can't deflate this power series (d = %ld): %Z", d, x);
      y = poltoser(poldeflate_i(y, d), v, 1 + (lx-3)/d);
      setvalp(y, V/d); return gerepilecopy(av, y);
    }
    if (checkdeflate(x) % d != 0) pari_err(cant_deflate);
    return gerepilecopy(av, poldeflate_i(x,d));
  }
  if (tx == t_RFRAC)
  {
    z = cgetg(3, t_RFRAC);
    gel(z,1) = gdeflate(gel(x,1),v,d);
    gel(z,2) = gdeflate(gel(x,2),v,d);
    return z;
  }
  if (is_matvec_t(tx))
  {
    lx = lg(x); z = cgetg(lx, tx);
    for (i=1; i<lx; i++) gel(z,i) = gdeflate(gel(x,i),v,d);
    return z;
  }
  pari_err(typeer,"gdeflate");
  return NULL; /* not reached */
}

/* set *m to the largest d such that x0 = A(X^d); return A */
GEN
poldeflate(GEN x, long *m)
{
  *m = checkdeflate(x);
  return poldeflate_i(x, *m);
}

/* return x0(X^d) */
GEN
polinflate(GEN x0, long d)
{
  long i, id, dy, dx = degpol(x0);
  GEN x = x0 + 2, z, y;
  dy = dx*d;
  y = cgetg(dy+3, t_POL); y[1] = x0[1];
  z = y + 2;
  for (i=0; i<=dy; i++) gel(z,i) = gen_0;
  for (i=id=0; i<=dx; i++,id+=d) z[id] = x[i];
  return y;
}

/* Distinct Degree Factorization (deflating first)
 * Assume x squarefree, degree(x) > 0, x(0) != 0 */
GEN
ZX_DDF(GEN x, long hint)
{
  GEN L;
  long m;
  x = poldeflate(x, &m);
  L = DDF(x, hint, 0);
  if (m > 1)
  {
    GEN e, v, fa = factoru(m);
    long i,j,k, l;

    e = gel(fa,2); k = 0;
    fa= gel(fa,1); l = lg(fa);
    for (i=1; i<l; i++) k += e[i];
    v = cgetg(k+1, t_VECSMALL); k = 1;
    for (i=1; i<l; i++)
      for (j=1; j<=e[i]; j++) v[k++] = fa[i];
    for (k--; k; k--)
    {
      GEN L2 = cgetg(1,t_VEC);
      for (i=1; i < lg(L); i++)
        L2 = shallowconcat(L2, DDF(polinflate(gel(L,i), v[k]), hint, 0));
      L = L2;
    }
  }
  return L;
}

/* SquareFree Factorization. f = prod P^e, all e distinct, in Z[X] (char 0
 * would be enough, if modulargcd --> ggcd). Return (P), set *ex = (e) */
GEN
ZX_squff(GEN f, GEN *ex)
{
  GEN T, V, W, P, e;
  long i, k, dW, n, val;

  if (signe(leading_term(f)) < 0) f = gneg_i(f);
  val = ZX_valuation(f, &f);
  n = 1 + degpol(f); if (val) n++;
  e = cgetg(n,t_VECSMALL);
  P = cgetg(n,t_COL);

  T = modulargcd(derivpol(f), f);
  V = RgX_div(f,T);
  for (k=i=1;; k++)
  {
    W = modulargcd(T,V); T = RgX_div(T,W); dW = degpol(W);
    /* W = prod P^e, e > k; V = prod P^e, e >= k */
    if (dW != degpol(V)) { gel(P,i) = RgX_div(V,W); e[i] = k; i++; }
    if (dW <= 0) break;
    V = W;
  }
  if (val) { gel(P,i) = pol_x[varn(f)]; e[i] = val; i++;}
  setlg(P,i);
  setlg(e,i); *ex = e; return P;
}

GEN
fact_from_DDF(GEN fa, GEN e, long n)
{
  GEN v,w, y = cgetg(3, t_MAT);
  long i,j,k, l = lg(fa);

  v = cgetg(n+1,t_COL); gel(y,1) = v;
  w = cgetg(n+1,t_COL); gel(y,2) = w;
  for (k=i=1; i<l; i++)
  {
    GEN L = gel(fa,i), ex = utoipos(e[i]);
    long J = lg(L);
    for (j=1; j<J; j++,k++)
    {
      gel(v,k) = gcopy(gel(L,j));
      gel(w,k) = ex;
    }
  }
  return y;
}

/* Factor x in Z[t]. Assume all factors have degree divisible by hint */
GEN
factpol(GEN x, long hint)
{
  pari_sp av = avma;
  GEN fa,ex,y;
  long n,i,l;

  if (typ(x)!=t_POL) pari_err(notpoler,"factpol");
  if (!signe(x)) pari_err(zeropoler,"factpol");

  fa = ZX_squff(Q_primpart(x), &ex);
  l = lg(fa); n = 0;
  for (i=1; i<l; i++)
  {
    gel(fa,i) = ZX_DDF(gel(fa,i), hint);
    n += lg(fa[i])-1;
  }
  y = fact_from_DDF(fa,ex,n);
  return gerepileupto(av, sort_factor(y, cmpii));
}

GEN
nfrootsQ(GEN x)
{
  pari_sp av = avma;
  GEN z, d;
  long val;
  
  if (typ(x)!=t_POL) pari_err(notpoler,"nfrootsQ");
  if (!signe(x)) pari_err(zeropoler,"nfrootsQ");
  val = ZX_valuation(Q_primpart(x), &x);
  d = modulargcd(derivpol(x), x);
  if (degpol(d)) x = RgX_div(x, d);
  z = DDF(x, 1, 1);
  if (val) z = shallowconcat(z, gen_0);
  return gerepilecopy(av, z);
}

/***********************************************************************/
/**                                                                   **/
/**                          FACTORIZATION                            **/
/**                                                                   **/
/***********************************************************************/
#define LT 17
#define assign_or_fail(x,y) {\
  if (y==NULL) y=x; else if (!gequal(x,y)) return 0;\
}
#define tsh 6
#define typs(x,y) ((x << tsh) | y)
#define typ1(x) (x >> tsh)
#define typ2(x) (x & ((1<<tsh)-1))

static long
poltype(GEN x, GEN *ptp, GEN *ptpol, long *ptpa)
{
  long t[LT]; /* code for 0,1,2,3,61,62,63,67,7,81,82,83,86,87,91,93,97 */
  long tx = typ(x),lx,i,j,s,pa=BIGINT;
  GEN pcx=NULL, p=NULL,pol=NULL,p1,p2;

  if (is_scalar_t(tx))
  {
    if (tx == t_POLMOD) return 0;
    x = scalarpol(x,0);
  }
  for (i=2; i<LT; i++) t[i]=0; /* t[0..1] unused */
  lx = lg(x);
  for (i=2; i<lx; i++)
  {
    p1=gel(x,i);
    switch(typ(p1))
    {
      case t_INT: case t_FRAC:
        break;
      case t_REAL:
        s = precision(p1); if (s < pa) pa = s;
        t[2]=1; break;
      case t_INTMOD:
	assign_or_fail(gel(p1,1),p);
        t[3]=1; break;
      case t_COMPLEX:
        if (!pcx) pcx = mkpoln(3, gen_1,gen_0,gen_1); /* x^2 + 1 */
	for (j=1; j<=2; j++)
	{
	  p2 = gel(p1,j);
	  switch(typ(p2))
	  {
	    case t_INT: case t_FRAC:
	      assign_or_fail(pcx,pol);
	      t[4]=1; break;
	    case t_REAL:
              s = precision(p2); if (s < pa) pa = s;
	      t[5]=1; break;
	    case t_INTMOD:
	      assign_or_fail(gel(p2,1),p);
	      assign_or_fail(pcx,pol);
	      t[6]=1; break;
	    case t_PADIC:
	      s = precp(p2) + valp(p2); if (s < pa) pa = s;
	      assign_or_fail(gel(p2,2),p);
	      assign_or_fail(pcx,pol);
	      t[7]=1; break;
	    default: return 0;
	  }
	}
	break;
      case t_PADIC:
        s = precp(p1) + valp(p1); if (s < pa) pa = s;
	assign_or_fail(gel(p1,2),p);
        t[8]=1; break;
      case t_QUAD:
	for (j=2; j<=3; j++)
	{
	  p2 = gel(p1,j);
	  switch(typ(p2))
	  {
	    case t_INT: case t_FRAC:
	      assign_or_fail(gel(p1,1),pol);
	      t[9]=1; break;
	    case t_REAL:
	      s = precision(p2); if (s < pa) pa = s;
	      if (gsigne(discsr(gel(p1,1)))>0) t[10]=1; else t[12]=1;
	      break;
	    case t_INTMOD:
	      assign_or_fail(gel(p2,1),p);
	      assign_or_fail(gel(p1,1),pol);
	      t[11]=1; break;
	    case t_PADIC:
	      s = precp(p2) + valp(p2); if (s < pa) pa = s;
	      assign_or_fail(gel(p2,2),p);
	      assign_or_fail(gel(p1,1),pol);
	      t[13]=1; break;
	    default: return 0;
	  }
	}
	break;
      case t_POLMOD:
	assign_or_fail(gel(p1,1),pol);
	for (j=1; j<=2; j++)
	{
	  GEN pbis = NULL, polbis = NULL;
	  long pabis;
	  switch(poltype(gel(p1,j),&pbis,&polbis,&pabis))
	  {
	    case t_INT: t[14]=1; break;
	    case t_INTMOD: t[15]=1; break;
	    case t_PADIC: t[16]=1; if (pabis<pa) pa=pabis; break;
	    default: return 0;
	  }
	  if (pbis) assign_or_fail(pbis,p);
	  if (polbis) assign_or_fail(polbis,pol);
	}
	break;
      default: return 0;
    }
  }
  if (t[5]||t[12])
  {
    if (t[3]||t[6]||t[7]||t[8]||t[11]||t[13]||t[14]||t[15]||t[16]) return 0;
    *ptpa=pa; return t_COMPLEX;
  }
  if (t[2]||t[10])
  {
    if (t[3]||t[6]||t[7]||t[8]||t[11]||t[13]||t[14]||t[15]||t[16]) return 0;
    *ptpa=pa; return t[4]?t_COMPLEX:t_REAL;
  }
  if (t[6]||t[11]||t[15])
  {
    *ptpol=pol; *ptp=p;
    i = t[15]? t_POLMOD: (t[11]? t_QUAD: t_COMPLEX);
    return typs(i, t_INTMOD);
  }
  if (t[7]||t[13]||t[16])
  {
    *ptpol=pol; *ptp=p; *ptpa=pa;
    i = t[16]? t_POLMOD: (t[13]? t_QUAD: t_COMPLEX);
    return typs(i, t_PADIC);
  }
  if (t[4]||t[9]||t[14])
  {
    *ptpol=pol;
    i = t[14]? t_POLMOD: (t[9]? t_QUAD: t_COMPLEX);
    return typs(i, t_INT);
  }
  if (t[3]) { *ptp=p; return t_INTMOD; }
  if (t[8]) { *ptp=p; *ptpa=pa; return t_PADIC; }
  return t_INT;
}
#undef LT

GEN
factor0(GEN x,long flag)
{
  long tx=typ(x);

  if (flag<0) return factor(x);
  if (is_matvec_t(tx)) return gboundfact(x,flag);
  if (tx==t_INT || tx==t_FRAC) return boundfact(x,flag);
  pari_err(talker,"partial factorization is not meaningful here");
  return NULL; /* not reached */
}

GEN
concat_factor(GEN f, GEN g)
{
  if (lg(f) == 1) return g;
  if (lg(g) == 1) return f;
  return mkmat2(shallowconcat(gel(f,1), gel(g,1)),
                shallowconcat(gel(f,2), gel(g,2)));
}

/* assume f and g coprime integer factorizations */
GEN
merge_factor_i(GEN f, GEN g)
{
  if (lg(f) == 1) return g;
  if (lg(g) == 1) return f;
  return sort_factor_gen(concat_factor(f,g), cmpii);
}

GEN
deg1_from_roots(GEN L, long v)
{
  long i, l = lg(L);
  GEN z = cgetg(l,t_COL);
  for (i=1; i<l; i++)
    gel(z,i) = deg1pol_i(gen_1, gneg(gel(L,i)), v);
  return z;
}
GEN
roots_from_deg1(GEN x)
{
  long i,l = lg(x);
  GEN r = cgetg(l,t_VEC);
  for (i=1; i<l; i++) gel(r,i) = gneg(constant_term(gel(x,i)));
  return r;
}

static GEN
gauss_factor_p(GEN p)
{ 
  GEN a, b; (void)cornacchia(gen_1, p, &a,&b);
  return mkcomplex(a, b);
}

static GEN
gauss_primpart(GEN x, GEN *c)
{
  GEN y, a = gel(x,1), b = gel(x,2), n = gcdii(a, b);
  *c = n; if (n == gen_1) return x;
  y = cgetg(3, t_COMPLEX);
  gel(y,1) = diviiexact(a, n);
  gel(y,2) = diviiexact(b, n); return y;
}

static GEN
gauss_primpart_try(GEN x, GEN c)
{
  GEN r, y;
  if (typ(x) == t_INT)
  {
    y = dvmdii(x, c, &r); if (r != gen_0) return NULL;
  }
  else
  {
    GEN a = gel(x,1), b = gel(x,2); y = cgetg(3, t_COMPLEX);
    gel(y,1) = dvmdii(a, c, &r); if (r != gen_0) return NULL;
    gel(y,2) = dvmdii(b, c, &r); if (r != gen_0) return NULL;
  }
  return y;
}

static int
gauss_cmp(GEN x, GEN y)
{
  int v;
  if (typ(x) != t_COMPLEX)
    return (typ(y) == t_COMPLEX)? -1: gcmp(x, y);
  if (typ(y) != t_COMPLEX) return 1;
  v = cmpii(gel(x,2), gel(y,2));
  if (v) return v;
  return gcmp(gel(x,1), gel(y,1));
}

/* 0 or canonical representative in Z[i]^* / <i> (impose imag(x) >= 0) */
static GEN
gauss_normal(GEN x)
{
  if (typ(x) != t_COMPLEX) return (signe(x) < 0)? absi(x): x;
  if (signe(x[1]) < 0) x = gneg(x);
  if (signe(x[2]) < 0) x = mulcxI(x);
  return x;
}

static GEN
Ipow(long e) {
  switch(e & 3)
  {
    case 1: return gi;
    case 2: return gen_m1;
    case 3: return pureimag(gen_m1);
  }
  return gen_1;
}

static GEN
gauss_factor(GEN x)
{
  pari_sp av = avma;
  GEN a = gel(x,1), b = gel(x,2), d = gen_1, n, y, fa, P, E, P2, E2;
  long t1 = typ(a);
  long t2 = typ(b), i, j, l, exp = 0;
  if (t1 == t_FRAC) d = gel(a,2);
  if (t2 == t_FRAC) d = lcmii(d, gel(b,2));
  if (d == gen_1) y = x;
  else
  {
    y = gmul(x, d);
    a = gel(y,1); t1 = typ(a);
    b = gel(y,2); t2 = typ(b);
  }
  if (t1 != t_INT || t2 != t_INT) return NULL;
  y = gauss_primpart(y, &n);
  fa = factor(cxnorm(y));
  P = gel(fa,1);
  E = gel(fa,2); l = lg(P);
  P2 = cgetg(l, t_COL);
  E2 = cgetg(l, t_COL);
  for (j = 1, i = l-1; i > 0; i--) /* remove largest factors first */
  {
    GEN p = gel(P,i), w, w2, t, we, pe;
    long v, e = itos(gel(E,i));
    int is2 = equaliu(p, 2);
    w = is2? mkcomplex(gen_1,gen_1): gauss_factor_p(p);
    w2 = gauss_normal( gconj(w) );
    /* w * w2 * I^3 = p, w2 = gconj(w) * I */
    pe = powiu(p, e);
    we = gpowgs(w, e);
    t = gauss_primpart_try( gmul(y, gconj(we)), pe );
    if (t) y = t; /* y /= w^e */
    else {
      /* y /= conj(w)^e, should be y /= w2^e */
      y = gauss_primpart_try( gmul(y, we), pe );
      swap(w, w2); exp += 3 * e;
    }
    gel(P,i) = w;
    v = Z_pvalrem(n, p, &n);
    if (v) {
      exp += 3*v;
      if (is2) v <<= 1; /* 2 = w^2 I^3 */
      else {
        gel(P2,j) = w2;
        gel(E2,j) = utoipos(v); j++;
      }
      gel(E,i) = stoi(e + v);
    }
    v = Z_pvalrem(d, p, &d);
    if (v) {
      exp -= 3*v;
      if (is2) v <<= 1; /* 2 is ramified */
      else {
        gel(P2,j) = w2;
        gel(E2,j) = utoineg(v); j++;
      }
      gel(E,i) = stoi(e - v);
    }
    exp &= 3;
  }
  if (j > 1) {
    setlg(P2, j);
    setlg(E2, j);
    fa = concat_factor(fa, mkmat2(P2,E2));
  }
  if (!is_pm1(n) || !is_pm1(d))
  {
    GEN Fa = factor(gdiv(n, d));
    P = gel(Fa,1); l = lg(P);
    E = gel(Fa,2);
    for (i = j = 1; i < l; i++)
    {
      GEN w, p = gel(P,i);
      long e;
      int is2;
      switch(mod4(p))
      {
        case 3: continue;
        case 2: is2 = 1; break;
        default:is2 = 0; break;
      }
      e = itos(gel(E,i));
      w = is2? mkcomplex(gen_1,gen_1): gauss_factor_p(p);
      gel(P,i) = w;
      if (is2)
        gel(E,i) = stoi(e << 1);
      else
      {
        P = shallowconcat(P, gauss_normal( gconj(w) ));
        E = shallowconcat(E, gel(E,i));
      }
      exp += 3*e;
      exp &= 3;
    }
    gel(Fa,1) = P;
    gel(Fa,2) = E;
    fa = concat_factor(fa, Fa);
  }
  fa = sort_factor_gen(fa, &gauss_cmp);

  y = gmul(y, Ipow(exp));
  if (!gcmp1(y)) {
    gel(fa,1) = shallowconcat(mkcol(y), gel(fa,1));
    gel(fa,2) = shallowconcat(gen_1,     gel(fa,2));
  }
  return gerepilecopy(av, fa);
}

GEN
factor(GEN x)
{
  long tx=typ(x), lx, i, j, pa, v, r1;
  pari_sp av, tetpil;
  GEN  y,p,p1,p2,p3,p4,p5,pol;

  if (is_matvec_t(tx))
  {
    lx=lg(x); y=cgetg(lx,tx);
    for (i=1; i<lx; i++) gel(y,i) = factor(gel(x,i));
    return y;
  }
  if (gcmp0(x))
  {
    y = cgetg(3,t_MAT);
    gel(y,1) = mkcolcopy(x);
    gel(y,2) = mkcol(gen_1); return y;
  }
  av = avma;
  switch(tx)
  {
    case t_INT: return Z_factor(x);

    case t_FRAC:
      p1 = Z_factor(gel(x,1));
      p2 = Z_factor(gel(x,2)); gel(p2,2) = gneg_i(gel(p2,2));
      return gerepilecopy(av, merge_factor_i(p1,p2));

    case t_POL:
      tx=poltype(x,&p,&pol,&pa);
      switch(tx)
      {
        case 0: pari_err(impl,"factor for general polynomials");
	case t_INT: return factpol(x,1);
	case t_INTMOD: return factmod(x,p);

	case t_COMPLEX: y=cgetg(3,t_MAT); lx=lg(x)-2; v=varn(x);
	  av = avma; p1 = roots(x,pa); tetpil = avma;
          p2 = deg1_from_roots(p1, v);
	  gel(y,1) = gerepile(av,tetpil,p2);
	  p3=cgetg(lx,t_COL); for (i=1; i<lx; i++) gel(p3,i) = gen_1;
          gel(y,2) = p3; return y;

	case t_REAL: y=cgetg(3,t_MAT); lx=lg(x)-2; v=varn(x);
	  av=avma; p1=cleanroots(x,pa); tetpil=avma;
	  for(r1=1; r1<lx; r1++)
            if (typ(p1[r1]) == t_COMPLEX) break;
	  lx=(r1+lx)>>1; p2=cgetg(lx,t_COL);
	  for(i=1; i<r1; i++)
            gel(p2,i) = deg1pol_i(gen_1, negr(gel(p1,i)), v);
	  for(   ; i<lx; i++)
	  {
	    GEN a = gel(p1,2*i-r1);
	    p = cgetg(5, t_POL); gel(p2,i) = p;
	    p[1] = x[1];
	    gel(p,2) = gnorm(a);
	    gel(p,3) = gmul2n(gel(a,1),1); setsigne(p[3],-signe(p[3]));
	    gel(p,4) = gen_1;
	  }
	  gel(y,1) = gerepile(av,tetpil,p2);
	  p3=cgetg(lx,t_COL); for (i=1; i<lx; i++) gel(p3,i) = gen_1;
          gel(y,2) = p3; return y;

	case t_PADIC: return factorpadic4(x,p,pa);

        default:
        {
          long killv;
	  x = shallowcopy(x); lx=lg(x);
          pol = shallowcopy(pol);
          v = manage_var(manage_var_max_avail,NULL);
          for(i=2; i<lx; i++)
          {
            p1=gel(x,i);
            switch(typ(p1))
            {
              case t_QUAD: p1++;
              case t_COMPLEX:
                p2 = cgetg(3, t_POLMOD); gel(x,i) = p2;
                gel(p2,1) = pol;
                gel(p2,2) = deg1pol_i(gel(p1,2), gel(p1,1), v);
            }
          }
          killv = (avma != (pari_sp)pol);
          if (killv) setvarn(pol, fetch_var());
          switch (typ2(tx))
          {
            case t_INT: p1 = polfnf(x,pol); break;
            case t_INTMOD: p1 = factorff(x,p,pol); break;
	    default: pari_err(impl,"factor of general polynomial");
              return NULL; /* not reached */
          }
          switch (typ1(tx))
          {
            case t_POLMOD:
              if (killv) (void)delete_var();
              return gerepileupto(av,p1);
            case t_COMPLEX: p5 = gi; break;
            case t_QUAD: p5=cgetg(4,t_QUAD);
              gel(p5,1) = pol; gel(p5,2) = gen_0; gel(p5,3) = gen_1;
              break;
	    default: pari_err(impl,"factor of general polynomial");
              return NULL; /* not reached */
          }
          p2=gel(p1,1);
          for(i=1; i<lg(p2); i++)
          {
            p3=gel(p2,i);
            for(j=2; j<lg(p3); j++)
            {
              p4=gel(p3,j);
              if(typ(p4)==t_POLMOD) gel(p3,j) = gsubst(gel(p4,2),v,p5);
            }
          }
          if (killv) (void)delete_var();
          return gerepilecopy(av, p1);
        }
      }
    case t_RFRAC:
      p1 = factor(gel(x,1));
      p2 = factor(gel(x,2)); gel(p2,2) = gneg_i(gel(p2,2));
      return gerepilecopy(av, concat_factor(p1,p2));

    case t_COMPLEX:
      y = gauss_factor(x);
      if (y) return y;
  }
  pari_err(talker,"can't factor %Z",x);
  return NULL; /* not reached */
}
#undef typ1
#undef typ2
#undef typs
#undef tsh

/* assume n > 0. Compute x^n using left-right binary powering */
GEN
leftright_pow_u_fold(GEN x, ulong n, void *data, GEN  (*sqr)(void*,GEN),
                                                 GEN (*msqr)(void*,GEN))
{
  GEN y;
  long m, j;

  if (n == 1) return gcopy(x);

  m = (long)n; j = 1+bfffo(m);
  y = x;

  /* normalize, i.e set highest bit to 1 (we know m != 0) */
  m<<=j; j = BITS_IN_LONG-j;
  /* first bit is now implicit */
  for (; j; m<<=1,j--)
  {
    if (m < 0) y = msqr(data,y); /* first bit set: multiply by base */
    else y = sqr(data,y);
  }
  return y;
}


/* assume n != 0, t_INT. Compute x^|n| using left-right binary powering */
GEN
leftright_pow_fold(GEN x, GEN n, void *data, GEN (*sqr)(void*,GEN),
                                             GEN (*msqr)(void*,GEN))
{
  long ln = lgefint(n);
  if (ln == 3) return leftright_pow_u_fold(x, n[2], data, sqr, msqr);
  else
  {
    GEN nd = int_MSW(n), y = x;
    long i, m = *nd, j = 1+bfffo((ulong)m);
    pari_sp av = avma, lim = stack_lim(av, 1);

    /* normalize, i.e set highest bit to 1 (we know m != 0) */
    m<<=j; j = BITS_IN_LONG-j;
    /* first bit is now implicit */
    for (i=ln-2;;)
    {
      for (; j; m<<=1,j--)
      {
        if (m < 0) y = msqr(data,y); /* first bit set: multiply by base */
        else y = sqr(data,y);
        if (low_stack(lim, stack_lim(av,1)))
        {
          if (DEBUGMEM>1) pari_warn(warnmem,"leftright_pow");
          y = gerepilecopy(av, y);
        }
      }
      if (--i == 0) return y;
      nd=int_precW(nd);
      m = *nd; j = BITS_IN_LONG;
    }
  }
}

struct leftright_fold
{
  void *data;
  GEN x;
  GEN (*sqr)(void*,GEN);
  GEN (*mul)(void*,GEN,GEN);
};

static GEN
leftright_sqr(void* data, GEN y)
{
  struct leftright_fold *d=(struct leftright_fold*) data;
  return d->sqr(d->data,y);
}

static GEN
leftright_msqr(void* data, GEN y)
{
  struct leftright_fold *d=(struct leftright_fold*) data;
  return d->mul(d->data,d->sqr(d->data,y),d->x);
}

GEN
leftright_pow(GEN x, GEN n, void *data, GEN (*sqr)(void*,GEN),
                                        GEN (*mul)(void*,GEN,GEN))
{
  struct leftright_fold d;
  d.data=data; d.mul=mul; d.sqr=sqr; d.x=x;
  return leftright_pow_fold(x, n, (void *)&d, leftright_sqr, leftright_msqr);
}

GEN
leftright_pow_u(GEN x, ulong n, void *data, GEN (*sqr)(void*,GEN),
                                        GEN (*mul)(void*,GEN,GEN))
{
  struct leftright_fold d;
  d.data=data; d.mul=mul; d.sqr=sqr; d.x=x;
  return leftright_pow_u_fold(x, n, (void *)&d, leftright_sqr, leftright_msqr);
}

GEN
divide_conquer_assoc(GEN x, GEN (*mul)(void *,GEN,GEN),void *data)
{
  pari_sp ltop, lim;
  long i,k,lx = lg(x);

  if (lx == 1) return gen_1;
  if (lx == 2) return gcopy(gel(x,1));
  x = shallowcopy(x); k = lx;
  ltop=avma; lim = stack_lim(ltop,1);
  while (k > 2)
  {
    if (DEBUGLEVEL>7)
      fprintferr("prod: remaining objects %ld\n",k-1);
    lx = k; k = 1;
    for (i=1; i<lx-1; i+=2)
      gel(x,k++) = mul(data,gel(x,i),gel(x,i+1));
    if (i < lx) x[k++] = x[i];
    if (low_stack(lim,stack_lim(av,1)))
      gerepilecoeffs(ltop,x+1,k-1);
  }
  return gel(x,1);
}

static GEN
_domul(void *data, GEN x, GEN y)
{
  GEN (*mul)(GEN,GEN)=(GEN (*)(GEN,GEN)) data;
  return mul(x,y);
}

GEN
divide_conquer_prod(GEN x, GEN (*mul)(GEN,GEN))
{
  return divide_conquer_assoc(x, _domul, (void *)mul);
}

static GEN
idmulred(void *nf, GEN x, GEN y) { return idealmulred((GEN) nf, x, y, 0); }
static GEN
idpowred(void *nf, GEN x, GEN n) { return idealpowred((GEN) nf, x, n, 0); }
static GEN
idmul(void *nf, GEN x, GEN y) { return idealmul((GEN) nf, x, y); }
static GEN
idpow(void *nf, GEN x, GEN n) { return idealpow((GEN) nf, x, n); }
static GEN
eltmul(void *nf, GEN x, GEN y) { return element_mul((GEN) nf, x, y); }
static GEN
eltpow(void *nf, GEN x, GEN n) { return element_pow((GEN) nf, x, n); }

#if 0
static GEN
ellmul(void *ell, GEN x, GEN y) { return addell((GEN) ell, x, y); }
static GEN
ellpow(void *GEN x, GEN n) { return powell((GEN) ell, x, n); }
#endif

GEN
factorback_aux(GEN fa, GEN e, GEN (*_mul)(void*,GEN,GEN), GEN (*_pow)(void*,GEN,GEN), void *data)
{
  pari_sp av = avma;
  long k,l,lx,t = typ(fa);
  GEN p,x;

  if (e) /* supplied vector of exponents */
    p = fa;
  else /* genuine factorization */
  {
    if (t == t_MAT) {
      l = lg(fa);
      if (l == 1) return gen_1;
      if (l != 3) pari_err(talker,"not a factorisation in factorback");
    } else {
      if (!is_vec_t(t)) pari_err(talker,"not a factorisation in factorback");
      return gerepileupto(av, divide_conquer_assoc(fa, _mul,data));
    }
    p = gel(fa,1);
    e = gel(fa,2);
  }
  lx = lg(p);
  t = t_INT; /* dummy */
  /* check whether e is an integral vector of correct length */
  if (is_vec_t(typ(e)) && lx == lg(e))
  {
    for (k=1; k<lx; k++)
      if (typ(e[k]) != t_INT) break;
    if (k == lx) t = t_MAT;
  }
  if (t != t_MAT) pari_err(talker,"not a factorisation in factorback");
  if (lx == 1) return gen_1;
  x = cgetg(lx,t_VEC);
  for (l=1,k=1; k<lx; k++)
    if (signe(e[k]))
      gel(x,l++) = _pow(data,gel(p,k),gel(e,k));
  setlg(x,l);
  return gerepileupto(av, divide_conquer_assoc(x, _mul,data));
}

static GEN _agmul(void *a, GEN x, GEN y) { return gmul(x,y);}
static GEN _apowgi(void *a, GEN x, GEN y) { return powgi(x,y);}

GEN
factorback_i(GEN fa, GEN e, GEN OBJ, int red)
{
  if (!OBJ)
  {
    if (e) {
      OBJ = checknf_i(e); if (OBJ) e = NULL;
    }
    if (!OBJ) return factorback_aux(fa, e, &_agmul, &_apowgi, NULL);
  }
  if (red) return factorback_aux(fa, e, &idmulred, &idpowred, OBJ);
  else     return factorback_aux(fa, e, &idmul,    &idpow, OBJ);
}

GEN
factorbackelt(GEN fa, GEN e, GEN nf)
{
  if (!nf && e && lg(e) > 1 && typ(e[1]) != t_INT) { nf = e; e = NULL; }
  if (!nf) pari_err(talker, "missing nf in factorbackelt");

  nf = checknf(nf);
  return factorback_aux(fa, e, &eltmul, &eltpow, nf);
}

GEN
factorback0(GEN fa, GEN e, GEN nf)
{
  return factorback_i(fa,e,nf,0);
}

GEN
factorback(GEN fa, GEN nf)
{
  return factorback_i(fa,nf,NULL,0);
}

GEN
gisirreducible(GEN x)
{
  long tx = typ(x), l, i;
  pari_sp av=avma;
  GEN y;

  if (is_matvec_t(tx))
  {
    l=lg(x); y=cgetg(l,tx);
    for (i=1; i<l; i++) gel(y,i) = gisirreducible(gel(x,i));
    return y;
  }
  if (is_intreal_t(tx) || tx == t_FRAC) return gen_0;
  if (tx!=t_POL) pari_err(notpoler,"gisirreducible");
  l=lg(x); if (l<=3) return gen_0;
  y=factor(x); avma=av;
  return (lg(gcoeff(y,1,1))==l)?gen_1:gen_0;
}

/*******************************************************************/
/*                                                                 */
/*                         GENERIC GCD                             */
/*                                                                 */
/*******************************************************************/
GEN
gcd0(GEN x, GEN y, long flag)
{
  if (!y) return content(x);
  switch(flag)
  {
    case 0: return ggcd(x,y);
    case 1: return modulargcd(x,y);
    case 2: return srgcd(x,y);
    default: pari_err(flagerr,"gcd");
  }
  return NULL; /* not reached */
}

/* x is a COMPLEX or a QUAD */
static GEN
triv_cont_gcd(GEN x, GEN y)
{
  pari_sp av = avma, tetpil;
  GEN p1 = (typ(x)==t_COMPLEX)? ggcd(gel(x,1),gel(x,2))
                              : ggcd(gel(x,2),gel(x,3));
  tetpil=avma; return gerepile(av,tetpil,ggcd(p1,y));
}

/* y is a PADIC, x a rational number or an INTMOD */
static GEN
padic_gcd(GEN x, GEN y)
{
  GEN p = gel(y,2);
  long v = ggval(x,p), w = valp(y);
  if (w < v) v = w;
  return gpowgs(p, v);
}

/* x,y in Z[i], at least one of which is t_COMPLEX */
static GEN
gauss_gcd(GEN x, GEN y)
{
  pari_sp av = avma;
  GEN dx, dy;
  dx = denom(x); x = gmul(x, dx);
  dy = denom(y); y = gmul(y, dy);
  while (!gcmp0(y))
  {
    GEN z = gsub(x, gmul(ground(gdiv(x,y)), y));
    x = y; y = z;
  }
  x = gauss_normal(x);
  if (typ(x) == t_COMPLEX)
  {
    if      (gcmp0(gel(x,2))) x = gel(x,1);
    else if (gcmp0(gel(x,1))) x = gel(x,2);
  }
  return gerepileupto(av, gdiv(x, lcmii(dx, dy)));
}

static int
c_is_rational(GEN x)
{
  return (is_rational(gel(x,1)) && is_rational(gel(x,2)));
}
static GEN
c_zero_gcd(GEN c)
{
  GEN x = gel(c,1), y = gel(c,2);
  long tx = typ(x), ty = typ(y);
  if (tx == t_REAL || ty == t_REAL) return gen_1;
  if (tx == t_PADIC || tx == t_INTMOD
   || ty == t_PADIC || ty == t_INTMOD) return ggcd(x, y);
  return gauss_gcd(c, gen_0);
}

/* y == 0 */
static GEN
zero_gcd(GEN x, long tx)
{
  pari_sp av;
  switch(tx)
  {
    case t_INT: return absi(x);
    case t_FRAC: return gabs(x,0);
    case t_COMPLEX: return c_zero_gcd(x);
    case t_REAL: return gen_1;
    case t_PADIC: return gpowgs(gel(x,2), valp(x));
    case t_SER: return monomial(gen_1, valp(x), varn(x));
    case t_POLMOD: {
      GEN d = gel(x,2);
      if (typ(d) == t_POL && varn(d) == varn(gel(x,1))) return content(d);
      return isinexact(d)? zero_gcd(d, typ(d)): gcopy(d);
    }
    case t_POL:
      if (!isinexact(x)) break;
      av = avma;
      return gerepileupto(av, 
        monomialcopy(content(x), polvaluation(x,NULL), varn(x))
      );

    case t_RFRAC:
      if (!isinexact(x)) break;
      av = avma;
      return gerepileupto(av,
        gdiv(zero_gcd(gel(x,1), typ(gel(x,1))), gel(x,2))
      );
  }
  return gcopy(x);
}

/* tx = t_RFRAC, y considered as constant */
static GEN
cont_gcd_rfrac(GEN x, GEN y)
{
  pari_sp av = avma;
  GEN cx; x = primitive_part(x, &cx);
  return gerepileupto(av, gred_rfrac_simple(ggcd(cx? cx: gen_1, y), gel(x,2)));
}
/* !is_const(tx), tx != t_RFRAC, y considered as constant */
static GEN
cont_gcd_gen(GEN x, GEN y)
{
  pari_sp av = avma;
  return gerepileupto(av, ggcd(content(x),y));
}
/* !is_const(tx), y considered as constant */
static GEN
cont_gcd(GEN x, long tx, GEN y)
{
  return (tx == t_RFRAC)? cont_gcd_rfrac(x, y): cont_gcd_gen(x, y);
}

GEN
ggcd(GEN x, GEN y)
{
  long l, i, vx, vy, tx = typ(x), ty = typ(y);
  pari_sp av, tetpil;
  GEN p1,z;

  if (is_noncalc_t(tx) || is_noncalc_t(ty)) pari_err(operf,"g",x,y);
  if (is_matvec_t(ty))
  {
    l = lg(y); z = cgetg(l,ty);
    for (i=1; i<l; i++) gel(z,i) = ggcd(x,gel(y,i));
    return z;
  }
  if (is_matvec_t(tx))
  {
    l = lg(x); z = cgetg(l,tx);
    for (i=1; i<l; i++) gel(z,i) = ggcd(gel(x,i),y);
    return z;
  }
  if (tx>ty) { swap(x,y); lswap(tx,ty); }
  /* tx <= ty */
  if (isexactzero(x)) return zero_gcd(y, ty);
  if (isexactzero(y)) return zero_gcd(x, tx);
  if (is_const_t(tx))
  {
    if (ty == tx) switch(tx)
    {
      case t_INT:
        return gcdii(x,y);

      case t_INTMOD: z=cgetg(3,t_INTMOD);
        if (equalii(gel(x,1),gel(y,1)))
          gel(z,1) = gcopy(gel(x,1));
        else
          gel(z,1) = gcdii(gel(x,1),gel(y,1));
        if (gcmp1(gel(z,1))) gel(z,2) = gen_0;
        else
        {
          av = avma; p1 = gcdii(gel(z,1),gel(x,2));
          if (!is_pm1(p1)) p1 = gerepileuptoint(av, gcdii(p1,gel(y,2)));
          gel(z,2) = p1;
        }
        return z;

      case t_FRAC: z=cgetg(3,t_FRAC);
        gel(z,1) = gcdii(gel(x,1), gel(y,1));
        gel(z,2) = lcmii(gel(x,2), gel(y,2));
        return z;

      case t_COMPLEX:
        if (c_is_rational(x) && c_is_rational(y)) return gauss_gcd(x,y);
        return triv_cont_gcd(y,x);

      case t_PADIC:
        if (!equalii(gel(x,2),gel(y,2))) return gen_1;
        vx = valp(x);
        vy = valp(y); return gpowgs(gel(y,2), min(vy,vx));

      case t_QUAD:
        av=avma; p1=gdiv(x,y);
        if (gcmp0(gel(p1,3)))
        {
          p1=gel(p1,2);
          if (typ(p1)==t_INT) { avma=av; return gcopy(y); }
          tetpil=avma; return gerepile(av,tetpil, gdiv(y,gel(p1,2)));
        }
        if (typ(p1[2])==t_INT && typ(p1[3])==t_INT) {avma=av; return gcopy(y);}
        p1 = ginv(p1); avma=av;
        if (typ(p1[2])==t_INT && typ(p1[3])==t_INT) return gcopy(x);
        return triv_cont_gcd(y,x);

      default: return gen_1; /* t_REAL */
    }
    if (is_const_t(ty)) switch(tx)
    {
      case t_INT:
        switch(ty)
        {
          case t_INTMOD: z = cgetg(3,t_INTMOD);
            gel(z,1) = gcopy(gel(y,1)); av = avma;
            p1 = gcdii(gel(y,1),gel(y,2));
            if (!is_pm1(p1)) p1 = gerepileuptoint(av, gcdii(x,p1));
            gel(z,2) = p1; return z;

          case t_FRAC: z = cgetg(3,t_FRAC);
            gel(z,1) = gcdii(x,gel(y,1));
            gel(z,2) = icopy(gel(y,2)); return z;

          case t_COMPLEX:
            if (c_is_rational(y)) return gauss_gcd(x,y);
          case t_QUAD:
            return triv_cont_gcd(y,x);

          case t_PADIC:
            return padic_gcd(x,y);

          default: return gen_1; /* t_REAL */
        }

      case t_INTMOD:
        switch(ty)
        {
          case t_FRAC:
            av = avma; p1=gcdii(gel(x,1),gel(y,2)); avma = av;
            if (!is_pm1(p1)) pari_err(operi,"g",x,y);
            return ggcd(gel(y,1), x);

          case t_COMPLEX: case t_QUAD:
            return triv_cont_gcd(y,x);

          case t_PADIC:
            return padic_gcd(x,y);
        }

      case t_FRAC:
        switch(ty)
        {
          case t_COMPLEX:
            if (c_is_rational(y)) return gauss_gcd(x,y);
          case t_QUAD:
            return triv_cont_gcd(y,x);

          case t_PADIC:
            return padic_gcd(x,y);
        }

      case t_COMPLEX: /* ty = PADIC or QUAD */
        return triv_cont_gcd(x,y);

      case t_PADIC: /* ty = QUAD */
        return triv_cont_gcd(y,x);

      default: return gen_1; /* tx = t_REAL */
    }
    return cont_gcd(y,ty, x);
  }

  if (tx == t_POLMOD)
  {
    if (ty == t_POLMOD)
    {
      z = cgetg(3,t_POLMOD);
      if (gequal(gel(x,1),gel(y,1)))
        gel(z,1) = gcopy(gel(x,1));
      else
        gel(z,1) = ggcd(gel(x,1),gel(y,1));
      if (degpol(z[1])<=0) gel(z,2) = gen_0;
      else
      {
        GEN X, Y, d;
        av = avma; X = gel(x,2); Y = gel(y,2);
        d = ggcd(content(X), content(Y));
        if (!gcmp1(d)) { X = gdiv(X,d); Y = gdiv(Y,d); }
        p1 = ggcd(gel(z,1), X);
        gel(z,2) = gerepileupto(av, gmul(d, ggcd(p1, Y)));
      }
      return z;
    }
    vx = varn(x[1]);
    switch(ty)
    {
      case t_POL:
        vy = varn(y);
        if (varncmp(vy,vx) < 0) return cont_gcd_gen(y, x);
        z = cgetg(3,t_POLMOD);
        gel(z,1) = gcopy(gel(x,1));
        av = avma; p1 = ggcd(gel(x,1),gel(x,2));
        gel(z,2) = gerepileupto(av, ggcd(p1,y));
        return z;

      case t_RFRAC:
        vy = varn(y[2]);
        if (varncmp(vy,vx) < 0) return cont_gcd_rfrac(y, x);
        av = avma; 
        p1 = ggcd(gel(x,1),gel(y,2));
        if (degpol(p1)) pari_err(operi,"g",x,y);
        avma = av; return gdiv(ggcd(gel(y,1),x), content(gel(y,2)));
    }
  }

  vx = gvar(x);
  vy = gvar(y);
  if (varncmp(vy, vx) < 0) return cont_gcd(y,ty, x);
  if (varncmp(vy, vx) > 0) return cont_gcd(x,tx, y);

  /* same main variable */
  switch(tx)
  {
    case t_POL:
      switch(ty)
      {
	case t_POL: return srgcd(x,y);
	case t_SER: 
          z = ggcd(content(x), content(y));
          return monomialcopy(z, min(valp(y),gval(x,vx)), vx);
	case t_RFRAC: return cont_gcd_rfrac(y, x);
      }
      break;

    case t_SER:
      z = ggcd(content(x), content(y));
      switch(ty)
      {
	case t_SER:   return monomialcopy(z, min(valp(x),valp(y)), vx);
	case t_RFRAC: return monomialcopy(z, min(valp(x),gval(y,vx)), vx);
      }
      break;

    case t_RFRAC: z=cgetg(3,t_RFRAC);
      if (ty != t_RFRAC) pari_err(operf,"g",x,y);
      p1 = gdeuc(gel(y,2), ggcd(gel(x,2), gel(y,2)));
      tetpil = avma;
      gel(z,2) = gerepile((pari_sp)z,tetpil,gmul(p1, gel(x,2)));
      gel(z,1) = ggcd(gel(x,1), gel(y,1)); return z;
  }
  pari_err(operf,"g",x,y);
  return NULL; /* not reached */
}

/* x a t_VEC,t_COL or t_MAT */
static GEN
vec_lcm(GEN x)
{
  if (typ(x) == t_MAT)
  {
    long i, l = lg(x);
    GEN z = cgetg(l, t_VEC);
    for (i = 1; i < l; i++) gel(z,i) = glcm0(gel(x,i), NULL);
    x = z;
  }
  return glcm0(x, NULL);
}
static GEN
scal_lcm(GEN x, GEN y)
{
  pari_sp av = avma;
  long tx = typ(x), ty = typ(y);
  if (is_matvec_t(tx)) x = vec_lcm(x);
  if (is_matvec_t(ty)) y = vec_lcm(y);
  return gerepileupto(av, glcm(x, y));
}

static GEN
fix_lcm(GEN x)
{
  GEN t;
  switch(typ(x))
  {
    case t_INT: if (signe(x)<0) x = negi(x);
      break;
    case t_POL:
      if (lg(x) <= 2) break;
      t = leading_term(x);
      if (typ(t) == t_INT && signe(t) < 0) x = gneg(x);
  }
  return x;
}

GEN
glcm0(GEN x, GEN y) {
  if (!y && lg(x) == 2)
  {
    long tx = typ(x);
    if (is_vec_t(tx))
    {
      x = gel(x,1);
      tx = typ(x);
      return is_matvec_t(tx)? vec_lcm(x): fix_lcm(x);
    }
  }
  return gassoc_proto(scal_lcm,x,y);
}

GEN
glcm(GEN x, GEN y)
{
  long tx, ty, i, l;
  pari_sp av;
  GEN p1,p2,z;

  ty = typ(y);
  if (is_matvec_t(ty))
  {
    l=lg(y); z=cgetg(l,ty);
    for (i=1; i<l; i++) gel(z,i) = glcm(x,gel(y,i));
    return z;
  }
  tx = typ(x);
  if (is_matvec_t(tx))
  {
    l=lg(x); z=cgetg(l,tx);
    for (i=1; i<l; i++) gel(z,i) = glcm(gel(x,i),y);
    return z;
  }
  if (tx==t_INT && ty==t_INT) return lcmii(x,y);
  if (gcmp0(x)) return gen_0;

  av = avma;
  p1 = ggcd(x,y); if (!gcmp1(p1)) y = gdiv(y,p1);
  p2 = gmul(x,y);
  return gerepileupto(av,fix_lcm(p2));
}

/* x + r ~ x ? Assume x,r are t_POL, deg(r) <= deg(x) */
static int
pol_approx0(GEN r, GEN x, int exact)
{
  long i, lx,lr;
  if (exact) return gcmp0(r);
  lx = lg(x);
  lr = lg(r); if (lr < lx) lx = lr;
  for (i=2; i<lx; i++)
    if (!approx_0(gel(r,i), gel(x,i))) return 0;
  return 1;
}

GEN
RgX_gcd_simple(GEN x, GEN y)
{
  pari_sp av1, av = avma, lim = stack_lim(av, 1);
  GEN r, yorig = y;
  int exact = !(isinexactreal(x) || isinexactreal(y));

  for(;;)
  {
    av1 = avma; r = grem(x,y);
    if (pol_approx0(r, x, exact))
    {
      avma = av1;
      if (y == yorig) return gcopy(y);
      y = normalizepol_approx(y, lg(y));
      if (lg(y) == 3) { avma = av; return gen_1; }
      return gerepileupto(av,y);
    }
    x = y; y = r;
    if (low_stack(lim,stack_lim(av,1))) {
      if(DEBUGMEM>1) pari_warn(warnmem,"RgX_gcd_simple");
      gerepileall(av,2, &x,&y);
    }
  }
}
GEN
RgX_extgcd_simple(GEN a, GEN b, GEN *pu, GEN *pv)
{
  pari_sp av = avma;
  GEN q, r, d, d1, u, v, v1;
  int exact = !(isinexactreal(a) || isinexactreal(b));

  d = a; d1 = b; v = gen_0; v1 = gen_1;
  for(;;)
  {
    if (pol_approx0(d1, a, exact)) break;
    q = poldivrem(d,d1, &r);
    v = gadd(v, gneg_i(gmul(q,v1)));
    u=v; v=v1; v1=u;
    u=r; d=d1; d1=u;
  }
  u = gadd(d, gneg_i(gmul(b,v)));
  u = RgX_div(u,a);

  gerepileall(av, 3, &u,&v,&d);
  return d;
}

static int issimplefield(GEN x);
static int
issimplepol(GEN x)
{
  long i, lx = lg(x);
  for (i=2; i<lx; i++)
    if (issimplefield(gel(x,i))) return 1;
  return 0;
}
/* return 1 if coeff explosion is not possible */
static int
issimplefield(GEN x)
{
  switch(typ(x))
  {
    case t_REAL: case t_INTMOD: case t_PADIC: case t_SER:
      return 1;
    case t_COMPLEX:
      return issimplefield(gel(x,1)) || issimplefield(gel(x,2));
    case t_POLMOD: 
      return (typ(x[2]) == t_POL && issimplepol(gel(x,2)))
           || issimplefield(gel(x,2)) || issimplepol(gel(x,1));
  }
  return 0;
}

static int
can_use_modular_gcd(GEN x)
{
  long i;
  for (i = lg(x)-1; i > 1; i--)
  {
    long t = typ(gel(x,i));
    if (!is_rational_t(t)) return 0;
  }
  return 1;
}

static GEN
gcdmonome(GEN x, GEN y)
{
  long dx=degpol(x), v=varn(x), e=gval(y, v);
  pari_sp av = avma;
  GEN t = ggcd(leading_term(x), content(y));

  if (dx < e) e = dx;
  return gerepileupto(av, monomialcopy(t, e, v));
}

/*******************************************************************/
/*                                                                 */
/*                    CONTENT / PRIMITIVE PART                     */
/*                                                                 */
/*******************************************************************/

GEN
content(GEN x)
{
  long lx, i, t, tx = typ(x);
  pari_sp av = avma;
  GEN c;

  if (is_scalar_t(tx)) return zero_gcd(x, tx);
  switch(tx)
  {
    case t_RFRAC:
    {
      GEN n = gel(x,1), d = gel(x,2);
      /* -- varncmp(vn, vd) < 0 can't happen
       * -- if n is POLMOD, its main variable (in the sense of gvar2)
       *    has lower priority than denominator */
      if (typ(n) == t_POLMOD || varncmp(gvar(n), varn(d)) > 0)
        n = isinexact(n)? zero_gcd(n, typ(n)): gcopy(n);
      else 
        n = content(n);
      return gerepileupto(av, gdiv(n, content(d)));
    }

    case t_VEC: case t_COL:
      lx = lg(x); if (lx==1) return gen_1;
      break;

    case t_MAT:
    {
      long hx, j;
      lx = lg(x);
      if (lx == 1) return gen_1;
      hx = lg(x[1]);
      if (hx == 1) return gen_1;
      if (lx == 2) { x = gel(x,1); lx = lg(x); break; }
      if (hx == 2) { x = row_i(x, 1, 1, lx-1); break; }
      c = content(gel(x,1));
      for (j=2; j<lx; j++)
        for (i=1; i<hx; i++) c = ggcd(c,gcoeff(x,i,j));
      if (typ(c) == t_INTMOD || isinexact(c)) { avma=av; return gen_1; }
      return gerepileupto(av,c);
    }

    case t_POL: case t_SER:
      lx = lg(x); if (lx == 2) return gen_0;
      break;
    case t_QFR: case t_QFI:
      lx = 4; break;

    default: pari_err(typeer,"content");
      return NULL; /* not reached */
  }
  for (i=lontyp[tx]; i<lx; i++)
    if (typ(x[i]) != t_INT) break;
  lx--; c = gel(x,lx);
  t = typ(c); if (is_matvec_t(t)) c = content(c);
  if (i > lx)
  { /* integer coeffs */
    while (lx-- > lontyp[tx])
    {
      c = gcdii(c, gel(x,lx));
      if (is_pm1(c)) { avma=av; return gen_1; }
    }
  }
  else
  {
    if (isinexact(c)) c = zero_gcd(c, typ(c));
    while (lx-- > lontyp[tx])
    {
      GEN d = gel(x,lx);
      t = typ(d); if (is_matvec_t(t)) d = content(d);
      c = ggcd(c, d);
    }
    if (typ(c) == t_INTMOD || isinexact(c)) { avma=av; return gen_1; }
  }
  switch(typ(c))
  {
    case t_INT:
      if (signe(c) < 0) c = negi(c);
      break;
    case t_VEC: case t_COL: case t_MAT:
      pari_err(typeer, "content");
  }

  return av==avma? gcopy(c): gerepileupto(av,c);
}

GEN
primitive_part(GEN x, GEN *ptc)
{
  pari_sp av = avma;
  GEN c = content(x);
  if (gcmp1(c)) { avma = av; c = NULL; }
  else if (!gcmp0(c)) x = gdiv(x,c);
  if (ptc) *ptc = c;
  return x;
}

/* NOT MEMORY CLEAN
 * As content(), but over Q. Treats polynomial as elts of Q[x1,...xn], instead
 * of Q(x2,...,xn)[x1] */
GEN
Q_content(GEN x)
{
  long i, l = lg(x);
  GEN d;
  pari_sp av;

  switch(typ(x))
  {
    case t_INT:  return absi(x);
    case t_FRAC: return gabs(x,0);

    case t_VEC: case t_COL: case t_MAT:
      l = lg(x); if (l == 1) return gen_1;
      av = avma; d = Q_content(gel(x,1));
      for (i=2; i<l; i++)
      {
        d = ggcd(d, Q_content(gel(x,i)));
        if ((i & 255) == 0) d = gerepileupto(av, d);
      }
      return gerepileupto(av, d);

    case t_POL:
      l = lg(x); if (l == 2) return gen_0;
      av = avma; d = Q_content(gel(x,2));
      for (i=3; i<l; i++)
        d = ggcd(d, Q_content(gel(x,i)));
      return gerepileupto(av, d);
    case t_POLMOD: return Q_content(gel(x,2));
    case t_COMPLEX: return ggcd(Q_content(gel(x,1)), Q_content(gel(x,2)));
  }
  pari_err(typeer,"Q_content");
  return NULL; /* not reached */
}

GEN
Q_primitive_part(GEN x, GEN *ptc)
{
  pari_sp av = avma;
  GEN c = Q_content(x);
  if (gcmp1(c)) { avma = av; c = NULL; }
  else if (!gcmp0(c)) x = Q_div_to_int(x,c);
  if (ptc) *ptc = c;
  return x;
}

GEN
primpart(GEN x) { return primitive_part(x, NULL); }

GEN
Q_primpart(GEN x) { return Q_primitive_part(x, NULL); }

/* NOT MEMORY CLEAN (because of t_FRAC).
 * As denom(), but over Q. Treats polynomial as elts of Q[x1,...xn], instead
 * of Q(x2,...,xn)[x1] */
GEN
Q_denom(GEN x)
{
  long i, l = lg(x);
  GEN d, D;
  pari_sp av;

  switch(typ(x))
  {
    case t_INT: return gen_1;
    case t_FRAC: return gel(x,2);

    case t_VEC: case t_COL: case t_MAT:
      l = lg(x); if (l == 1) return gen_1;
      av = avma; d = Q_denom(gel(x,1));
      for (i=2; i<l; i++)
      {
        D = Q_denom(gel(x,i));
        if (D != gen_1) d = lcmii(d, D);
        if ((i & 255) == 0) d = gerepileuptoint(av, d);
      }
      return gerepileuptoint(av, d);

    case t_POL:
      l = lg(x); if (l == 2) return gen_1;
      av = avma; d = Q_denom(gel(x,2));
      for (i=3; i<l; i++)
      {
        D = Q_denom(gel(x,i));
        if (D != gen_1) d = lcmii(d, D);
      }
      return gerepileuptoint(av, d);
  }
  pari_err(typeer,"Q_denom");
  return NULL; /* not reached */
}

GEN
Q_remove_denom(GEN x, GEN *ptd)
{
  GEN d = Q_denom(x);
  if (gcmp1(d)) d = NULL; else x = Q_muli_to_int(x,d);
  if (ptd) *ptd = d;
  return x;
}

/* return y = x * d, assuming x rational, and d,y integral */
GEN
Q_muli_to_int(GEN x, GEN d)
{
  long i, l, t = typ(x);
  GEN y, xn, xd;
  pari_sp av;

  if (typ(d) != t_INT) pari_err(typeer,"Q_muli_to_int");
  switch (t)
  {
    case t_INT:
      return mulii(x,d);

    case t_FRAC:
      xn = gel(x,1);
      xd = gel(x,2); av = avma;
      y = mulii(xn, diviiexact(d, xd));
      return gerepileuptoint(av, y);

    case t_VEC: case t_COL: case t_MAT:
      l = lg(x); y = cgetg(l, t);
      for (i=1; i<l; i++) gel(y,i) = Q_muli_to_int(gel(x,i), d);
      return y;

    case t_POL:
      l = lg(x); y = cgetg(l, t_POL); y[1] = x[1];
      for (i=2; i<l; i++) gel(y,i) = Q_muli_to_int(gel(x,i), d);
      return y;

    case t_POLMOD:
      y = cgetg(3, t_POLMOD);
      gel(y,1) = gcopy(gel(x,1));
      gel(y,2) = Q_muli_to_int(gel(x,2), d);
      return y;
  }
  pari_err(typeer,"Q_muli_to_int");
  return NULL; /* not reached */
}

/* return x * n/d. x: rational; d,n,result: integral. n = NULL represents 1 */
GEN
Q_divmuli_to_int(GEN x, GEN d, GEN n)
{
  long i, l, t = typ(x);
  GEN y, xn, xd;
  pari_sp av;
  
  switch(t)
  {
    case t_INT:
      av = avma; y = diviiexact(x,d);
      if (n) y = gerepileuptoint(av, mulii(y,n));
      return y;

    case t_FRAC:
      xn = gel(x,1);
      xd = gel(x,2); av = avma;
      y = mulii(diviiexact(xn, d), diviiexact(n, xd));
      return gerepileuptoint(av, y);

    case t_VEC: case t_COL: case t_MAT:
      l = lg(x); y = cgetg(l, t);
      for (i=1; i<l; i++) gel(y,i) = Q_divmuli_to_int(gel(x,i), d,n);
      return y;

    case t_POL:
      l = lg(x); y = cgetg(l, t_POL); y[1] = x[1];
      for (i=2; i<l; i++) gel(y,i) = Q_divmuli_to_int(gel(x,i), d,n);
      return y;

    case t_POLMOD:
      y = cgetg(3, t_POLMOD);
      gel(y,1) = gcopy(gel(x,1));
      gel(y,2) = Q_divmuli_to_int(gel(x,2), d,n);
      return y;
  }
  pari_err(typeer,"Q_divmuli_to_int");
  return NULL; /* not reached */
}

/* return y = x / c, assuming x,c rational, and y integral */
GEN
Q_div_to_int(GEN x, GEN c)
{
  GEN d, n;
  switch(typ(c))
  {
    case t_INT:
      return Q_divmuli_to_int(x, c, NULL);
    case t_FRAC:
      n = gel(c,1);
      d = gel(c,2); if (gcmp1(n)) return Q_muli_to_int(x,d);
      return Q_divmuli_to_int(x, n,d);
  }
  pari_err(typeer,"Q_div_to_int");
  return NULL; /* not reached */
}

/*******************************************************************/
/*                                                                 */
/*                           SUBRESULTANT                          */
/*                                                                 */
/*******************************************************************/
/* for internal use */
GEN
gdivexact(GEN x, GEN y)
{
  long i,lx;
  GEN z;
  if (gcmp1(y)) return x;
  switch(typ(x))
  {
    case t_INT:
      if (typ(y)==t_INT) return diviiexact(x,y);
      if (!signe(x)) return gen_0;
      break;
    case t_INTMOD:
    case t_POLMOD: return gmul(x,ginv(y));
    case t_POL:
      switch(typ(y))
      {
        case t_INTMOD:
        case t_POLMOD: return gmul(x,ginv(y));
        case t_POL:
          if (varn(x)==varn(y)) return poldivrem(x,y, NULL);
      }
      lx = lg(x); z = new_chunk(lx);
      for (i=2; i<lx; i++) gel(z,i) = gdivexact(gel(x,i),y);
      z[1] = x[1]; 
      z[0] = x[0]; return z;
    case t_VEC: case t_COL: case t_MAT:
      lx = lg(x); z = new_chunk(lx);
      for (i=1; i<lx; i++) gel(z,i) = gdivexact(gel(x,i),y);
      z[0] = x[0]; return z;
  }
  if (DEBUGLEVEL) pari_warn(warner,"missing case in gdivexact");
  return gdiv(x,y);
}

static GEN
init_resultant(GEN x, GEN y)
{
  long tx,ty;
  if (gcmp0(x) || gcmp0(y)) return gen_0;
  tx = typ(x); ty = typ(y);
  if (is_scalar_t(tx) || is_scalar_t(ty))
  {
    if (tx==t_POL) return gpowgs(y,degpol(x));
    if (ty==t_POL) return gpowgs(x,degpol(y));
    return gen_1;
  }
  if (tx!=t_POL || ty!=t_POL) pari_err(typeer,"subresall");
  if (varn(x)==varn(y)) return NULL;
  return (varn(x)<varn(y))? gpowgs(y,degpol(x)): gpowgs(x,degpol(y));
}

/* return coefficients s.t x = x_0 X^n + ... + x_n */
GEN
revpol(GEN x)
{
  long i,n = degpol(x);
  GEN y = cgetg(n+3,t_POL);
  y[1] = x[1]; x += 2; y += 2;
  for (i=0; i<=n; i++) y[i] = x[n-i];
  return y;
}

/* assume dx >= dy, y non constant, mod either NULL or a t_POL. */
static GEN
pseudorem_i(GEN x, GEN y, GEN mod)
{
  long vx = varn(x), dx, dy, dz, i, lx, p;
  pari_sp av = avma, av2, lim;

  if (!signe(y)) pari_err(gdiver);
  (void)new_chunk(2);
  dx=degpol(x); x = revpol(x);
  dy=degpol(y); y = revpol(y); dz=dx-dy; p = dz+1;
  av2 = avma; lim = stack_lim(av2,1);
  for (;;)
  {
    gel(x,0) = gneg(gel(x,0)); p--;
    for (i=1; i<=dy; i++)
    {
      gel(x,i) = gadd(gmul(gel(y,0), gel(x,i)), gmul(gel(x,0),gel(y,i)));
      if (mod) gel(x,i) = RgX_rem(gel(x,i), mod);
    }
    for (   ; i<=dx; i++)
    {
      gel(x,i) = gmul(gel(y,0), gel(x,i));
      if (mod) gel(x,i) = RgX_rem(gel(x,i), mod);
    }
    do { x++; dx--; } while (dx >= 0 && gcmp0(gel(x,0)));
    if (dx < dy) break;
    if (low_stack(lim,stack_lim(av2,1)))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"pseudorem dx = %ld >= %ld",dx,dy);
      gerepilecoeffs(av2,x,dx+1);
    }
  }
  if (dx < 0) return zeropol(vx);
  lx = dx+3; x -= 2;
  x[0] = evaltyp(t_POL) | evallg(lx);
  x[1] = evalsigne(1) | evalvarn(vx);
  x = revpol(x) - 2;
  if (p)
  { /* multiply by y[0]^p   [beware dummy vars from FpY_FpXY_resultant] */
    GEN t = gel(y,0);
    if (mod)
    { /* assume p fairly small */
      for (i=1; i<p; i++)
        t = RgX_rem(gmul(t, gel(y,0)), mod);
    }
    else
      t = gpowgs(t, p);
    for (i=2; i<lx; i++)
    {
      gel(x,i) = gmul(gel(x,i), t);
      if (mod) gel(x,i) = RgX_rem(gel(x,i), mod);
    }
    if (!mod) return gerepileupto(av, x);
  }
  return gerepilecopy(av, x);
}

GEN
pseudorem(GEN x, GEN y) { return pseudorem_i(x,y, NULL); }

/* assume dx >= dy, y non constant
 * Compute z,r s.t lc(y)^(dx-dy+1) x = z y + r */
GEN
pseudodiv(GEN x, GEN y, GEN *ptr)
{
  long vx = varn(x), dx, dy, dz, i, iz, lx, lz, p;
  pari_sp av = avma, av2, lim;
  GEN z, r, ypow;

  if (!signe(y)) pari_err(gdiver);
  (void)new_chunk(2);
  dx=degpol(x); x = revpol(x);
  dy=degpol(y); y = revpol(y); dz=dx-dy; p = dz+1;
  lz = dz+3; z = cgetg(lz, t_POL) + 2;
  ypow = new_chunk(dz+1);
  gel(ypow,0) = gen_1;
  for (i=1; i<=dz; i++) gel(ypow,i) = gmul(gel(ypow,i-1), gel(y,0));
  av2 = avma; lim = stack_lim(av2,1);
  for (iz=0;;)
  {
    p--;
    gel(z,iz++) = gmul(gel(x,0), gel(ypow,p));
    gel(x,0) = gneg(gel(x,0));
    for (i=1; i<=dy; i++)
      gel(x,i) = gadd(gmul(gel(y,0), gel(x,i)), gmul(gel(x,0),gel(y,i)));
    for (   ; i<=dx; i++)
      gel(x,i) = gmul(gel(y,0), gel(x,i));
    x++; dx--;
    while (dx >= dy && gcmp0(gel(x,0))) { x++; dx--; gel(z,iz++) = gen_0; }
    if (dx < dy) break;
    if (low_stack(lim,stack_lim(av2,1)))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"pseudodiv dx = %ld >= %ld",dx,dy);
      gerepilecoeffs2(av2,x,dx+1, z,iz);
    }
  }
  while (dx >= 0 && gcmp0(gel(x,0))) { x++; dx--; }
  if (dx < 0)
    x = zeropol(vx);
  else
  {
    lx = dx+3; x -= 2;
    x[0] = evaltyp(t_POL) | evallg(lx);
    x[1] = evalsigne(1) | evalvarn(vx);
    x = revpol(x) - 2;
  }

  z -= 2;
  z[0] = evaltyp(t_POL) | evallg(lz);
  z[1] = evalsigne(1) | evalvarn(vx);
  z = revpol(z) - 2;
  r = gmul(x, gel(ypow,p));
  gerepileall(av, 2, &z, &r);
  *ptr = r; return z;
}

/* Return resultant(u,v). If sol != NULL: set *sol to the last non-zero
 * polynomial in the prs IF the sequence was computed, and gen_0 otherwise */
GEN
subresall(GEN u, GEN v, GEN *sol)
{
  pari_sp av, av2, lim;
  long degq,dx,dy,du,dv,dr,signh;
  GEN z,g,h,r,p1,p2,cu,cv;

  if (sol) *sol=gen_0;
  if ((r = init_resultant(u,v))) return r;

  if (isinexact(u) || isinexact(v)) return resultant2(u,v);
  dx=degpol(u); dy=degpol(v); signh=1;
  if (dx < dy)
  {
    swap(u,v); lswap(dx,dy);
    if (both_odd(dx, dy)) signh = -signh;
  }
  if (dy==0) return gpowgs(gel(v,2),dx);
  av = avma;
  u = primitive_part(u, &cu);
  v = primitive_part(v, &cv);
  g = h = gen_1; av2 = avma; lim = stack_lim(av2,1);
  for(;;)
  {
    r = pseudorem(u,v); dr = lg(r);
    if (dr == 2)
    {
      if (sol) { avma = (pari_sp)(r+2); *sol=gerepileupto(av,v); } else avma = av;
      return gen_0;
    }
    du = degpol(u); dv = degpol(v); degq = du-dv;
    u = v; p1 = g; g = leading_term(u);
    switch(degq)
    {
      case 0: break;
      case 1:
        p1 = gmul(h,p1); h = g; break;
      default:
        p1 = gmul(gpowgs(h,degq),p1);
        h = gdivexact(gpowgs(g,degq), gpowgs(h,degq-1));
    }
    if (both_odd(du,dv)) signh = -signh;
    v = gdivexact(r,p1);
    if (dr==3) break;
    if (low_stack(lim,stack_lim(av2,1)))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"subresall, dr = %ld",dr);
      gerepileall(av2,4, &u, &v, &g, &h);
    }
  }
  z = gel(v,2);
  if (dv > 1) z = gdivexact(gpowgs(z,dv), gpowgs(h,dv-1));
  if (signh < 0) z = gneg(z); /* z = resultant(ppart(x), ppart(y)) */
  p2 = gen_1;
  if (cu) p2 = gmul(p2, gpowgs(cu,dy));
  if (cv) p2 = gmul(p2, gpowgs(cv,dx));
  z = gmul(z, p2);

  if (sol) u = gclone(u);
  z = gerepileupto(av, z);
  if (sol) { *sol = gcopy(u); gunclone(u); }
  return z;
}

static GEN
scalar_res(GEN x, GEN y, GEN *U, GEN *V)
{
  *V=gpowgs(y,degpol(x)-1); *U=gen_0; return gmul(y,*V);
}

/* compute U, V s.t Ux + Vy = resultant(x,y) */
GEN
subresext(GEN x, GEN y, GEN *U, GEN *V)
{
  pari_sp av, av2, tetpil, lim;
  long dx, dy, signh, tx = typ(x), ty = typ(y);
  GEN z, g, h, p1, cu, cv, u, v, um1, uze, vze;
  GEN *gptr[3];

  if (!is_extscalar_t(tx) || !is_extscalar_t(ty)) pari_err(typeer,"subresext");
  if (gcmp0(x) || gcmp0(y)) { *U = *V = gen_0; return gen_0; }
  if (tx != t_POL) {
    if (ty != t_POL) { *U = ginv(x); *V = gen_0; return gen_1; }
    return scalar_res(y,x,V,U);
  }
  if (ty != t_POL) return scalar_res(x,y,U,V);
  if (varn(x) != varn(y))
    return varncmp(varn(x), varn(y)) < 0? scalar_res(x,y,U,V)
                                        : scalar_res(y,x,V,U);
  dx = degpol(x); dy = degpol(y); signh = 1;
  if (dx < dy)
  {
    pswap(U,V); lswap(dx,dy); swap(x,y);
    if (both_odd(dx, dy)) signh = -signh;
  }
  if (dy == 0)
  {
    *V = gpowgs(gel(y,2),dx-1);
    *U = gen_0; return gmul(*V,gel(y,2));
  }
  av = avma;
  u = x = primitive_part(x, &cu);
  v = y = primitive_part(y, &cv);
  g = h = gen_1; av2 = avma; lim = stack_lim(av2,1);
  um1 = gen_1; uze = gen_0;
  for(;;)
  {
    GEN r, q = pseudodiv(u,v, &r);
    long du, dv, degq, dr = lg(r);
    if (dr == 2) { *U = *V = gen_0; avma = av; return gen_0; }

    du = degpol(u); dv = degpol(v); degq = du-dv;
    /* lead(v)^(degq + 1) * um1 - q * uze */
    p1 = gsub(gmul(gpowgs(gel(v,dv+2),degq+1),um1), gmul(q,uze));
    um1 = uze; uze = p1;
    u = v; p1 = g; g = leading_term(u);
    switch(degq)
    {
      case 0: break;
      case 1: p1 = gmul(h,p1); h = g; break;
      default:
        p1 = gmul(gpowgs(h,degq),p1);
        h = gdivexact(gpowgs(g,degq), gpowgs(h,degq-1));
    }
    if (both_odd(du, dv)) signh = -signh;
    v  = gdivexact(r,p1);
    uze= gdivexact(uze,p1);
    if (dr == 3) {
      z = gel(v,2);
      if (dv > 1)
      { /* z = gdivexact(gpowgs(z,dv), gpowgs(h,dv-1)); */
        p1 = gpowgs(gdiv(z,h),dv-1);
        z = gmul(z,p1); uze = gmul(uze,p1);
      }
      break;
    }
    if (low_stack(lim,stack_lim(av2,1)))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"subresext, dr = %ld",dr);
      gerepileall(av2,6, &u,&v,&g,&h,&uze,&um1);
    }
  }
  if (signh < 0) { z = gneg_i(z); uze = gneg_i(uze); }
  p1 = gadd(z, gneg(gmul(uze,x)));
  vze = RgX_divrem(p1, y, &p1);
  if (!gcmp0(p1)) pari_warn(warner,"inexact computation in subresext");
  /* uze ppart(x) + vze ppart(y) = z = resultant(ppart(x), ppart(y)), */
  p1 = gen_1;
  if (cu) p1 = gmul(p1, gpowgs(cu,dy));
  if (cv) p1 = gmul(p1, gpowgs(cv,dx));
  cu = cu? gdiv(p1,cu): p1;
  cv = cv? gdiv(p1,cv): p1;

  tetpil = avma;
  z = gmul(z,p1);
  *U = gmul(uze,cu);
  *V = gmul(vze,cv);
  gptr[0]=&z; gptr[1]=U; gptr[2]=V;
  gerepilemanysp(av,tetpil,gptr,3);
  return z;
}

static GEN
scalar_bezout(GEN x, GEN y, GEN *U, GEN *V)
{
  *U=gen_0; *V=ginv(y); return pol_1[varn(x)];
}

static GEN
zero_bezout(GEN y, GEN *U, GEN *V)
{
  GEN x=content(y);
  *U=gen_0; *V = ginv(x); return gmul(y,*V);
}

/* compute U, V s.t Ux + Vy = GCD(x,y) using subresultant */
GEN
RgX_extgcd(GEN x, GEN y, GEN *U, GEN *V)
{
  pari_sp av, av2, tetpil, lim;
  long dx, dy, tx = typ(x), ty = typ(y);
  GEN z, g, h, p1, cu, cv, u, v, um1, uze, vze, *gptr[3];

  if (!is_extscalar_t(tx) || !is_extscalar_t(ty)) pari_err(typeer,"subresext");
  if (gcmp0(x)) {
    if (gcmp0(y)) { *U = *V = gen_0; return gen_0; }
    return zero_bezout(y,U,V);
  }
  if (gcmp0(y)) return zero_bezout(x,V,U);
  if (tx != t_POL) {
    if (ty != t_POL) { *U = ginv(x); *V = gen_0; return pol_1[0]; }
    return scalar_bezout(y,x,V,U);
  }
  if (ty != t_POL) return scalar_bezout(x,y,U,V);
  if (varn(x) != varn(y))
    return varncmp(varn(x), varn(y)) < 0? scalar_bezout(x,y,U,V)
                                        : scalar_bezout(y,x,V,U);
  dx = degpol(x); dy = degpol(y);
  if (dx < dy)
  {
    pswap(U,V); lswap(dx,dy); swap(x,y);
  }
  if (dy==0) return scalar_bezout(x,y,U,V);

  av = avma;
  u = x = primitive_part(x, &cu);
  v = y = primitive_part(y, &cv);
  g = h = gen_1; av2 = avma; lim = stack_lim(av2,1);
  um1 = gen_1; uze = gen_0;
  for(;;)
  {
    GEN r, q = pseudodiv(u,v, &r);
    long du, dv, degq, dr = lg(r);
    if (dr == 2) break;

    du = degpol(u); dv = degpol(v); degq = du-dv;
    p1 = gsub(gmul(gpowgs(gel(v,dv+2),degq+1),um1), gmul(q,uze));
    um1 = uze; uze = p1;
    u = v; p1 = g; g  = leading_term(u);
    switch(degq)
    {
      case 0: break;
      case 1:
        p1 = gmul(h,p1); h = g; break;
      default:
        p1 = gmul(gpowgs(h,degq), p1);
        h = gdiv(gpowgs(g,degq), gpowgs(h,degq-1));
    }
    v  = gdivexact(r,p1);
    uze= gdivexact(uze,p1);
    if (dr==3) break;
    if (low_stack(lim,stack_lim(av2,1)))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"RgX_extgcd, dr = %ld",dr);
      gerepileall(av2,6,&u,&v,&g,&h,&uze,&um1);
    }
  }
  p1 = gadd(v, gneg(gmul(uze,x)));
  vze = RgX_divrem(p1, y, &p1);
  if (!gcmp0(p1)) pari_warn(warner,"inexact computation in RgX_extgcd");
  if (cu) uze = gdiv(uze,cu);
  if (cv) vze = gdiv(vze,cv);
  p1 = ginv(content(v));
  
  tetpil = avma;
  *U = gmul(uze,p1);
  *V = gmul(vze,p1);
  z = gmul(v,p1);
  gptr[0]=U; gptr[1]=V; gptr[2]=&z;
  gerepilemanysp(av,tetpil,gptr,3); return z;
}

/*******************************************************************/
/*                                                                 */
/*                RESULTANT USING DUCOS VARIANT                    */
/*                                                                 */
/*******************************************************************/

static GEN
reductum(GEN P)
{
  if (signe(P)==0) return P;
  return normalizepol_i(shallowcopy(P),lg(P)-1);
}

static GEN
Lazard(GEN x, GEN y, long n)
{
  long a, b;
  GEN c;

  if (n<=1) return x;
  a=1; while (n >= (b=a+a)) a=b;
  c=x; n-=a;
  while (a>1)
  {
    a>>=1; c=gdivexact(gsqr(c),y);
    if (n>=a) { c=gdivexact(gmul(c,x),y); n -= a; }
  }
  return c;
}

static GEN
Lazard2(GEN F, GEN x, GEN y, long n)
{
  if (n<=1) return F;
  return gdivexact(gmul(Lazard(x,y,n-1), F), y);
}

/* deg(P) > deg(Q) */
static GEN
nextSousResultant(GEN P, GEN Q, GEN Z, GEN s)
{
  GEN p0, q0, z0, H, A;
  long pr, p, q, j, v = varn(P);
  pari_sp av, lim;

  z0 = leading_term(Z);
  p = degpol(P); p0 = leading_term(P); P = reductum(P);
  q = degpol(Q); q0 = leading_term(Q); Q = reductum(Q);

  av = avma; lim = stack_lim(av,1);
  H = gneg(reductum(Z));
  pr = degpol(P);
  A = (q <= pr)? gmul(gel(P,q+2),H): gen_0;
  for (j = q+1; j < p; j++)
  {
    H = (degpol(H) == q-1) ?
      addshift(reductum(H), gdivexact(gmul(gneg(gel(H,q+1)),Q), q0)) :
      addshift(H, zeropol(v));
    if (j <= pr) A = gadd(A,gmul(gel(P,j+2),H));
    if (low_stack(lim,stack_lim(av,1)))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"nextSousResultant j = %ld/%ld",j,p);
      gerepileall(av,2,&A,&H);
    }
  }
  P = normalizepol_i(P, min(pr+3,q+2));
  A = gdivexact(gadd(A,gmul(z0,P)), p0);
  A = (degpol(H) == q-1) ?
    gadd(gmul(q0,addshift(reductum(H),A)), gmul(gneg(gel(H,q+1)),Q)) :
    gmul(q0, addshift(H,A));
  return gdivexact(A, ((p-q)&1)? s: gneg(s));
}

GEN
resultantducos(GEN P, GEN Q)
{
  pari_sp av = avma, av2, lim;
  long dP,dQ,delta;
  GEN cP,cQ,Z,s;

  if ((Z = init_resultant(P,Q))) return Z;
  dP = degpol(P);
  dQ = degpol(Q);
  P = primitive_part(P, &cP);
  Q = primitive_part(Q, &cQ);
  delta = dP - dQ;
  if (delta < 0)
  {
    Z = (dP & dQ & 1)? gneg(Q): Q;
    Q = P; P = Z; delta = -delta;
  }
  s = gen_1;
  if (degpol(Q) > 0)
  {
    av2 = avma; lim = stack_lim(av2,1);
    s = gpowgs(leading_term(Q),delta);
    Z = Q;
    Q = pseudorem(P, gneg(Q));
    P = Z;
    while(degpol(Q) > 0)
    {
      if (low_stack(lim,stack_lim(av,1)))
      {
        if(DEBUGMEM>1) pari_warn(warnmem,"resultantducos, degpol Q = %ld",degpol(Q));
        gerepileall(av2,2,&P,&Q); s = leading_term(P);
      }
      delta = degpol(P) - degpol(Q);
      Z = Lazard2(Q, leading_term(Q), s, delta);
      Q = nextSousResultant(P, Q, Z, s);
      P = Z;
      s = leading_term(P);
    }
  }
  if (!signe(Q)) { avma = av; return gen_0; }
  if (!degpol(P)){ avma = av; return gen_1; }
  s = Lazard(leading_term(Q), s, degpol(P));
  if (cP) s = gmul(s, gpowgs(cP,dQ));
  if (cQ) s = gmul(s, gpowgs(cQ,dP)); else if (!cP) s = gcopy(s);
  return gerepileupto(av, s);
}

/*******************************************************************/
/*                                                                 */
/*               RESULTANT USING SYLVESTER MATRIX                  */
/*                                                                 */
/*******************************************************************/
static GEN
_zeropol(void)
{
  GEN x = cgetg(3,t_POL);
  x[1] = 0;
  gel(x,2) = gen_0; return x;
}

static GEN
sylvester_col(GEN x, long j, long d, long D)
{
  GEN c = cgetg(d+1,t_COL);
  long i;
  for (i=1; i< j; i++) gel(c,i) = gen_0;
  for (   ; i<=D; i++) c[i]=x[D-i+2];
  for (   ; i<=d; i++) gel(c,i) = gen_0;
  return c;
}

GEN
sylvestermatrix_i(GEN x, GEN y)
{
  long j,d,dx,dy;
  GEN M;

  dx = degpol(x); if (dx < 0) { dx = 0; x = _zeropol(); }
  dy = degpol(y); if (dy < 0) { dy = 0; y = _zeropol(); }
  d = dx+dy; M = cgetg(d+1,t_MAT);
  for (j=1; j<=dy; j++) gel(M,j) = sylvester_col(x,j,d,j+dx);
  for (j=1; j<=dx; j++) gel(M,j+dy) = sylvester_col(y,j,d,j+dy);
  return M;
}

GEN
sylvestermatrix(GEN x, GEN y)
{
  long i,j,lx;
  if (typ(x)!=t_POL || typ(y)!=t_POL) pari_err(typeer,"sylvestermatrix");
  if (varn(x) != varn(y))
    pari_err(talker,"not the same variables in sylvestermatrix");
  x = sylvestermatrix_i(x,y); lx = lg(x);
  for (i=1; i<lx; i++)
    for (j=1; j<lx; j++) gcoeff(x,i,j) = gcopy(gcoeff(x,i,j));
  return x;
}

GEN
resultant2(GEN x, GEN y)
{
  pari_sp av;
  GEN r;
  if ((r = init_resultant(x,y))) return r;
  av = avma; return gerepileupto(av,det(sylvestermatrix_i(x,y)));
}

static GEN
fix_pol(GEN x, long v, long *mx)
{
  long vx;
  GEN p1;

  if (typ(x) == t_POL)
  {
    vx = varn(x);
    if (vx)
    {
      if (v>=vx) return gsubst(x,v,pol_x[0]);
      p1 = cgetg(3,t_POL);
      p1[1] = evalvarn(0)|evalsigne(signe(x));
      gel(p1,2) = x; return p1;
    }
    if (v)
    {
      *mx = 1;
      return gsubst(gsubst(x,0,pol_x[MAXVARN]),v,pol_x[0]);
    }
  }
  return x;
}

/* resultant of x and y with respect to variable v, or with respect to their
 * main variable if v < 0. */
GEN
polresultant0(GEN x, GEN y, long v, long flag)
{
  long m = 0;
  pari_sp av = avma;

  if (v >= 0)
  {
    x = fix_pol(x,v, &m);
    y = fix_pol(y,v, &m);
  }
  switch(flag)
  {
    case 0: x=subresall(x,y,NULL); break;
    case 1: x=resultant2(x,y); break;
    case 2: x=resultantducos(x,y); break;
    default: pari_err(flagerr,"polresultant");
  }
  if (m) x = gsubst(x,MAXVARN,pol_x[0]);
  return gerepileupto(av,x);
}

/*******************************************************************/
/*                                                                 */
/*                  GCD USING SUBRESULTANT                         */
/*                                                                 */
/*******************************************************************/
GEN
srgcd(GEN x, GEN y)
{
  long tx = typ(x), ty = typ(y), dx, dy, vx;
  pari_sp av, av1, tetpil, lim;
  GEN d, g, h, p1, p2, u, v;

  if (!signe(y)) return gcopy(x);
  if (!signe(x)) return gcopy(y);
  if (is_scalar_t(tx) || is_scalar_t(ty)) return gen_1;
  if (tx!=t_POL || ty!=t_POL) pari_err(typeer,"srgcd");
  vx=varn(x); if (vx!=varn(y)) return gen_1;
  if (ismonome(x)) return gcdmonome(x,y);
  if (ismonome(y)) return gcdmonome(y,x);
  av = avma;
  if (can_use_modular_gcd(x) &&
      can_use_modular_gcd(y)) return modulargcd(x,y); /* Q[X] */

  if (issimplepol(x) || issimplepol(y)) x = RgX_gcd_simple(x,y);
  else
  {
    dx=lg(x); dy=lg(y);
    if (dx<dy) { swap(x,y); lswap(dx,dy); }
    p1=content(x); p2=content(y); d=ggcd(p1,p2);

    tetpil=avma; d = scalarpol(d, vx);
    if (dy==3) return gerepile(av,tetpil,d);

    av1=avma; lim=stack_lim(av1,1);
    u=gdiv(x,p1); v=gdiv(y,p2); g=h=gen_1;
    for(;;)
    {
      GEN r = pseudorem(u,v);
      long degq, du, dv, dr=lg(r);

      if (dr <= 3)
      {
        if (gcmp0(r)) break;
        avma=av1; return gerepile(av,tetpil,d);
      }
      if (DEBUGLEVEL > 9) fprintferr("srgcd: dr = %ld\n", dr);
      du=lg(u); dv=lg(v); degq=du-dv; u=v;
      switch(degq)
      {
        case 0:
          v = gdiv(r,g);
          g = leading_term(u);
          break;
        case 1:
          v = gdiv(r,gmul(h,g));
          h = g = leading_term(u);
          break;
        default:
          v = gdiv(r,gmul(gpowgs(h,degq),g));
          g = leading_term(u);
          h = gdiv(gpowgs(g,degq), gpowgs(h,degq-1));
      }
      if (low_stack(lim, stack_lim(av1,1)))
      {
        if(DEBUGMEM>1) pari_warn(warnmem,"srgcd");
        gerepileall(av1,4,&u,&v,&g,&h);
      }
    }
    p1 = content(v); if (!gcmp1(p1)) v = gdiv(v,p1);
    x = gmul(d,v);
  }

  if (typ(x)!=t_POL) x = scalarpol(x, vx);
  else
  {
    p1=leading_term(x); ty=typ(p1);
    if ((ty == t_FRAC || is_intreal_t(ty)) && gsigne(p1)<0) x = gneg(x);
  }
  return gerepileupto(av,x);
}

GEN
poldisc0(GEN x, long v)
{
  long tx=typ(x), i;
  pari_sp av;
  GEN z,p1,p2;

  switch(tx)
  {
    case t_POL:
      if (gcmp0(x)) return gen_0;
      av = avma; i = 0;
      if (v >= 0 && v != varn(x)) x = fix_pol(x,v, &i);
      p1 = subres(x, derivpol(x));
      p2 = leading_term(x); if (!gcmp1(p2)) p1 = gdiv(p1,p2);
      if (degpol(x) & 2) p1 = gneg(p1);
      if (i) p1 = gsubst(p1, MAXVARN, pol_x[0]);
      return gerepileupto(av,p1);

    case t_COMPLEX:
      return utoineg(4);

    case t_QUAD: case t_POLMOD:
      return poldisc0(gel(x,1), v);

    case t_QFR: case t_QFI:
      av = avma; return gerepileuptoint(av, qf_disc(x));

    case t_VEC: case t_COL: case t_MAT:
      i=lg(x); z=cgetg(i,tx);
      for (i--; i; i--) gel(z,i) = poldisc0(gel(x,i), v);
      return z;
  }
  pari_err(typeer,"discsr");
  return NULL; /* not reached */
}

GEN
discsr(GEN x) { return poldisc0(x, -1); }

GEN
reduceddiscsmith(GEN pol)
{
  long i, j, n;
  pari_sp av = avma;
  GEN polp, p1, m;

  if (typ(pol)!=t_POL) pari_err(typeer,"reduceddiscsmith");
  n=degpol(pol);
  if (n<=0) pari_err(constpoler,"reduceddiscsmith");
  check_ZX(pol,"poldiscreduced");
  if (!gcmp1(gel(pol,n+2)))
    pari_err(talker,"non-monic polynomial in poldiscreduced");
  polp = derivpol(pol);
  m=cgetg(n+1,t_MAT);
  for (j=1; j<=n; j++)
  {
    p1=cgetg(n+1,t_COL); gel(m,j) = p1;
    for (i=1; i<=n; i++) gel(p1,i) = truecoeff(polp,i-1);
    if (j<n) polp = grem(RgX_shift_shallow(polp, 1), pol);
  }
  return gerepileupto(av, smith(m));
}

/***********************************************************************/
/**							              **/
/**	                  STURM ALGORITHM                             **/
/**	         (number of real roots of x in ]a,b])		      **/
/**								      **/
/***********************************************************************/

/* if a (resp. b) is NULL, set it to -oo (resp. +oo) */
long
sturmpart(GEN x, GEN a, GEN b)
{
  long sl, sr, s, t, r1;
  pari_sp av = avma, lim = stack_lim(av, 1);
  GEN g,h,u,v;

  if (gcmp0(x)) pari_err(zeropoler,"sturm");
  t = typ(x);
  if (t != t_POL)
  {
    if (t == t_INT || t == t_REAL || t == t_FRAC) return 0;
    pari_err(typeer,"sturm");
  }
  s=lg(x); if (s==3) return 0;

  sl = gsigne(leading_term(x));
  if (s==4)
  {
    t = a? gsigne(poleval(x,a)): -sl;
    if (t == 0) { avma = av; return 0; }
    s = b? gsigne(poleval(x,b)):  sl;
    avma = av; return (s == t)? 0: 1;
  }
  u=gdiv(x,content(x)); v=derivpol(x);
  v=gdiv(v,content(v));
  g=gen_1; h=gen_1;
  s = b? gsigne(poleval(u,b)): sl;
  t = a? gsigne(poleval(u,a)): ((lg(u)&1)? sl: -sl);
  r1=0;
  sr = b? gsigne(poleval(v,b)): s;
  if (sr)
  {
    if (!s) s=sr;
    else if (sr!=s) { s= -s; r1--; }
  }
  sr = a? gsigne(poleval(v,a)): -t;
  if (sr)
  {
    if (!t) t=sr;
    else if (sr!=t) { t= -t; r1++; }
  }
  for(;;)
  {
    GEN p1, r = pseudorem(u,v);
    long du=lg(u), dv=lg(v), dr=lg(r), degq=du-dv;

    if (dr<=2) pari_err(talker,"not a squarefree polynomial in sturm");
    if (gsigne(leading_term(v)) > 0 || degq&1) r=gneg_i(r);
    sl = gsigne(gel(r,dr-1));
    sr = b? gsigne(poleval(r,b)): sl;
    if (sr)
    {
      if (!s) s=sr;
      else if (sr!=s) { s= -s; r1--; }
    }
    sr = a? gsigne(poleval(r,a)): ((dr&1)? sl: -sl);
    if (sr)
    {
      if (!t) t=sr;
      else if (sr!=t) { t= -t; r1++; }
    }
    if (dr==3) { avma=av; return r1; }

    u=v; p1 = g; g = gabs(leading_term(u),DEFAULTPREC);
    switch(degq)
    {
      case 0: break;
      case 1:
        p1 = gmul(h,p1); h = g; break;
      default:
        p1 = gmul(gpowgs(h,degq),p1);
        h = gdivexact(gpowgs(g,degq), gpowgs(h,degq-1));
    }
    v = gdivexact(r,p1);
    if (low_stack(lim,stack_lim(av,1)))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"polsturm, dr = %ld",dr);
      gerepileall(av,4,&u,&v,&g,&h);
    }
  }
}

/***********************************************************************/
/**                                                                   **/
/**                        GENERIC EXTENDED GCD                       **/
/**                                                                   **/
/***********************************************************************/

static GEN
RgXQ_inv_inexact(GEN x, GEN y)
{
  pari_sp av = avma;
  long i, dx = degpol(x), dy = degpol(y), dz = dx+dy;
  GEN v, z;

  if (dx < 0 || dy < 0) pari_err(talker,"non-invertible polynomial in RgXQ_inv");
  v = gauss(sylvestermatrix(y,x), col_ei(dz, dz));
  z = cgetg(dy+2,t_POL); z[1] = y[1];
  for (i=2; i<dy+2; i++) z[i] = v[dz-i+2];
  return gerepilecopy(av, normalizepol_i(z, dy+2));
}
/* assume typ(x) = t_POL */
static GEN
RgXQ_inv(GEN x, GEN y)
{
  long vx=varn(x), vy=varn(y);
  pari_sp av, av1;
  GEN u, v, d;

  while (vx != vy)
  {
    if (varncmp(vx,vy) > 0)
    {
      if (vx == BIGINT) return ginv(x);
      return gred_rfrac_simple(gen_1, x);
    }
    if (lg(x)!=3) pari_err(talker,"non-invertible polynomial in RgXQ_inv");
    x = gel(x,2); vx = gvar(x);
  }
  if (isinexact(x) || isinexact(y)) return RgXQ_inv_inexact(x,y);

  av = avma; d = subresext(x,y,&u,&v);
  if (gcmp0(d)) pari_err(talker,"non-invertible polynomial in RgXQ_inv");
  if (typ(d) == t_POL && varn(d) == vx)
  {
    if (lg(d) > 3) pari_err(talker,"non-invertible polynomial in RgXQ_inv");
    d = gel(d,2);
  }
  av1 = avma; return gerepile(av,av1, gdiv(u,d));
}

GEN
gbezout(GEN x, GEN y, GEN *u, GEN *v)
{
  if (typ(x) == t_INT && typ(y) == t_INT) return bezout(x,y,u,v);
  return RgX_extgcd(x,y,u,v);
}

GEN
vecbezout(GEN x, GEN y)
{
  GEN z=cgetg(4,t_VEC);
  gel(z,3) = gbezout(x,y,(GEN*)(z+1),(GEN*)(z+2));
  return z;
}

GEN
vecbezoutres(GEN x, GEN y)
{
  GEN z=cgetg(4,t_VEC);
  gel(z,3) = subresext(x,y,(GEN*)(z+1),(GEN*)(z+2));
  return z;
}


/*******************************************************************/
/*                                                                 */
/*                    GENERIC (modular) INVERSE                    */
/*                                                                 */
/*******************************************************************/

GEN
ginvmod(GEN x, GEN y)
{
  long tx=typ(x);

  switch(typ(y))
  {
    case t_POL:
      if (tx==t_POL) return RgXQ_inv(x,y);
      if (is_scalar_t(tx)) return ginv(x);
      break;
    case t_INT:
      if (tx==t_INT) return Fp_inv(x,y);
      if (tx==t_POL) return gen_0;
  }
  pari_err(typeer,"ginvmod");
  return NULL; /* not reached */
}

/***********************************************************************/
/**							              **/
/**		            NEWTON POLYGON                            **/
/**								      **/
/***********************************************************************/

/* assume leading coeff of x is non-zero */
GEN
newtonpoly(GEN x, GEN p)
{
  GEN y;
  long n,ind,a,b,c,u1,u2,r1,r2;
  long *vval, num[] = {evaltyp(t_INT) | _evallg(3), 0, 0};

  if (typ(x)!=t_POL) pari_err(typeer,"newtonpoly");
  n=degpol(x); if (n<=0) { y=cgetg(1,t_VEC); return y; }
  y = cgetg(n+1,t_VEC); x += 2; /* now x[i] = term of degree i */
  vval = (long *) gpmalloc(sizeof(long)*(n+1));
  for (a=0; a<=n; a++) vval[a] = ggval(gel(x,a),p);
  for (a=0, ind=1; a<n; a++)
  {
    if (vval[a] != VERYBIGINT) break;
    gel(y,ind++) = utoipos(VERYBIGINT);
  }
  for (b=a+1; b<=n; a=b, b=a+1)
  {
    while (vval[b] == VERYBIGINT) b++;
    u1=vval[a]-vval[b]; u2=b-a;
    for (c=b+1; c<=n; c++)
    {
      if (vval[c] == VERYBIGINT) continue;
      r1=vval[a]-vval[c]; r2=c-a;
      if (u1*r2<=u2*r1) { u1=r1; u2=r2; b=c; }
    }
    while (ind<=b) { affsi(u1,num); gel(y,ind++) = gdivgs(num,u2); }
  }
  free(vval); return y;
}

static GEN
lift_to_frac(GEN t, GEN mod, GEN amax, GEN bmax, GEN denom)
{
  GEN a, b;
  if (signe(t) < 0) t = addii(t, mod); /* in case t is a centerlift */
  if (!ratlift(t, mod, &a, &b, amax, bmax)
     || (denom && !dvdii(denom,b))
     || !gcmp1(gcdii(a,b))) return NULL;
  if (!is_pm1(b)) a = mkfrac(a, b);
  return a;
}

/* compute rational lifting for all the components of M modulo mod.
 * If one components fails, return NULL.
 * See ratlift.
 * If denom is not NULL, check that the denominators divide denom
 *
 * FIXME: NOT stack clean ! a & b stay on the stack.
 * If we suppose mod and denom coprime, then a and b are coprime
 * so we can do a cgetg(t_FRAC). */
GEN
matratlift(GEN M, GEN mod, GEN amax, GEN bmax, GEN denom)
{
  pari_sp ltop = avma;
  GEN N, a;
  long l2, l3;
  long i, j;
  if (typ(M)!=t_MAT) pari_err(typeer,"matratlift");
  l2 = lg(M); l3 = lg(gel(M,1));
  N = cgetg(l2, t_MAT);
  for (j = 1; j < l2; ++j)
  {
    gel(N,j) = cgetg(l3, t_COL);
    for (i = 1; i < l3; ++i)
    {
      a = lift_to_frac(gcoeff(M,i,j), mod, amax,bmax,denom);
      if (!a) { avma = ltop; return NULL; }
      gcoeff(N, i, j) = a;
    }
  }
  return N;
}

GEN
polratlift(GEN P, GEN mod, GEN amax, GEN bmax, GEN denom)
{
  pari_sp ltop = avma;
  GEN Q,a;
  long l2, j;
  if (typ(P)!=t_POL) pari_err(typeer,"polratlift");
  l2 = lg(P);
  Q = cgetg(l2, t_POL); Q[1] = P[1];
  for (j = 2; j < l2; ++j)
  {
    a = lift_to_frac(gel(P,j), mod, amax,bmax,denom);
    if (!a) { avma = ltop; return NULL; }
    gel(Q,j) = a;
  }
  return Q;
}

/* P,Q in Z[X,Y], nf in Z[Y] irreducible. compute GCD in Q[Y]/(nf)[X].
 *
 * We essentially follow M. Encarnacion "On a modular Algorithm for computing
 * GCDs of polynomials over number fields" (ISSAC'94).
 *
 * We procede as follows
 *  1:compute the gcd modulo primes discarding bad primes as they are detected.
 *  2:reconstruct the result via matratlift, stoping as soon as we get weird
 *    denominators.
 *  3:if matratlift succeeds, try the full division.
 * Suppose accuracy is insufficient to get the result right: matratlift will
 * rarely succeed, and even if it does the polynomial we get has sensible
 * coefficients, so the full division will not be too costly.
 *
 * FIXME: Handle rational coefficient for P and Q.
 * If not NULL, den must a a multiple of the denominator of the gcd,
 * for example the discriminant of nf.
 *
 * NOTE: if nf is not irreducible, nfgcd may loop forever, esp. if gcd | nf */
GEN
nfgcd(GEN P, GEN Q, GEN nf, GEN den)
{
  pari_sp ltop = avma;
  GEN sol, mod = NULL;
  long x = varn(P);
  long y = varn(nf);
  long d = degpol(nf);
  GEN lP, lQ;
  if (!signe(P) || !signe(Q)) return zeropol(x);
  /*Compute denominators*/
  if (!den) den = ZX_disc(nf);
  lP = leading_term(P);
  lQ = leading_term(Q);
  if ( !((typ(lP)==t_INT && is_pm1(lP)) || (typ(lQ)==t_INT && is_pm1(lQ))) )
    den = mulii(den, gcdii(ZX_resultant(lP, nf), ZX_resultant(lQ, nf)));
  { /*Modular GCD*/
    pari_sp btop = avma, st_lim = stack_lim(btop, 1);
    ulong p;
    long dM=0, dR;
    GEN M, dsol;
    GEN R, ax, bo;
    byteptr primepointer;
    for (p = 27449, primepointer = diffptr + 3000; ; )
    {
      NEXT_PRIME_VIADIFF_CHECK(p, primepointer);
      /*Discard primes dividing disc(T) or lc(PQ) */
      if (!smodis(den, p)) continue;
      if (DEBUGLEVEL>5) fprintferr("nfgcd: p=%d\n",p);
      /*Discard primes when modular gcd does not exist*/
      if ((R = FlxqX_safegcd(ZXX_to_FlxX(P,p,y),
                             ZXX_to_FlxX(Q,p,y),
                             ZX_to_Flx(nf,p), p)) == NULL) continue;
      dR = degpol(R);
      if (dR == 0) return scalarpol(gen_1, x);
      if (mod && dR > dM) continue; /* p divides Res(P/gcd, Q/gcd). Discard. */

      R = RgXX_to_RgM(FlxX_to_ZXX(R), d);
      /* previous primes divided Res(P/gcd, Q/gcd)? Discard them. */
      if (!mod || dR < dM) { M = R; mod = utoipos(p); dM = dR; continue; }
      if (low_stack(st_lim, stack_lim(btop, 1)))
      {
	if (DEBUGMEM>1) pari_warn(warnmem,"nfgcd");
	gerepileall(btop, 2, &M, &mod);
      }

      ax = mulis(Fp_inv(utoipos(p), mod), p);
      M = gadd(R, gmul(ax, gsub(M, R)));
      mod = mulis(mod, p);
      M = lift(FpM_to_mod(M, mod));
      /* I suspect it must be better to take amax > bmax*/
      bo = sqrti(shifti(mod, -1));
      if ((sol = matratlift(M, mod, bo, bo, den)) == NULL) continue;
      sol = RgM_to_RgXX(sol,x,y);
      dsol = primpart(sol);
      if (gcmp0(pseudorem_i(P, dsol, nf))
       && gcmp0(pseudorem_i(Q, dsol, nf))) break;
    }
  }
  return gerepilecopy(ltop, sol);
}
