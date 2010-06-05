/* $Id: bibli1.c 7977 2006-08-03 17:12:01Z kb $

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
/**                 LLL Algorithm and close friends                **/
/**                                                                **/
/********************************************************************/
#include "pari.h"
#include "paripriv.h"

/* default quality ratio for LLL: 99/100 */
#define LLLDFT 100

/* scalar product x.x */
GEN
sqscal(GEN x)
{
  long i, lx;
  pari_sp av;
  GEN z;
  lx = lg(x);
  if (lx == 1) return gen_0;
  av = avma;
  z = gsqr(gel(x,1));
  for (i=2; i<lx; i++)
    z = gadd(z, gsqr(gel(x,i)));
  return gerepileupto(av,z);
}

/* scalar product x.y */
GEN
gscal(GEN x,GEN y)
{
  long i, lx;
  pari_sp av;
  GEN z;
  if (x == y) return sqscal(x);
  lx = lg(x);
  if (lx == 1) return gen_0;
  av = avma;
  z = gmul(gel(x,1),gel(y,1));
  for (i=2; i<lx; i++)
    z = gadd(z, gmul(gel(x,i),gel(y,i)));
  return gerepileupto(av,z);
}

static GEN
sqscali(GEN x)
{
  long i, lx;
  pari_sp av;
  GEN z;
  lx = lg(x);
  if (lx == 1) return gen_0;
  av = avma;
  z = sqri(gel(x,1));
  for (i=2; i<lx; i++)
    z = addii(z, sqri(gel(x,i)));
  return gerepileuptoint(av,z);
}

static GEN
gscali(GEN x,GEN y)
{
  long i, lx;
  pari_sp av;
  GEN z;
  if (x == y) return sqscali(x);
  lx = lg(x);
  if (lx == 1) return gen_0;
  av = avma;
  z = mulii(gel(x,1),gel(y,1));
  for (i=2; i<lx; i++)
    z = addii(z, mulii(gel(x,i),gel(y,i)));
  return gerepileuptoint(av,z);
}

/********************************************************************/
/**             QR Factorization via Householder matrices          **/
/********************************************************************/
static int
no_prec_pb(GEN x)
{
  return (typ(x) != t_REAL || lg(x) >  3
                           || expo(x) < (long)BITS_IN_HALFULONG);
}
/* zero x[1..k-1], fill mu */
static int
FindApplyQ(GEN x, GEN mu, GEN B, long k, GEN Q, long prec)
{
  long i, lx = lg(x)-1, lv = lx - (k-1);
  GEN v, beta, Nx, x2, x1, xd = x + (k-1);

  x1 = gel(xd,1);
  x2 = gsqr(x1);
  if (k < lx)
  {
    for (i=2; i<=lv; i++) x2 = mpadd(x2, gsqr(gel(xd,i)));
    Nx = gsqrt(x2, prec);
    if (signe(x1) < 0) setsigne(Nx, -1);
    v = cgetg(lv+1, t_VEC);
    gel(v,1) = mpadd(x1, Nx);
    for (i=2; i<=lv; i++) v[i] = xd[i];
    if (gcmp0(x2)) return 0;

    if (gcmp0(x1))
      beta = mpmul(x2, real_1(prec)); /* make sure typ(beta) != t_INT */
    else
      beta = mpadd(x2, mpmul(Nx,x1));
    gel(Q,k) = mkvec2(ginv(beta), v);

    gcoeff(mu,k,k) = mpneg(Nx);
  }
  else
    coeff(mu,k,k) = x[k];
  if (B)
  {
    gel(B,k) = x2;
    for (i=1; i<k; i++) coeff(mu,k,i) = x[i];
  }
  else
    for (i=1; i<k; i++) coeff(mu,i,k) = x[i];
  return no_prec_pb(x2);
}

static void
ApplyQ(GEN Q, GEN r)
{
  GEN s, rd, beta = gel(Q,1), v = gel(Q,2);
  long i, l = lg(v), lr = lg(r);

  rd = r + (lr - l);
  s = mpmul(gel(v,1), gel(rd,1));
  for (i=2; i<l; i++) s = mpadd(s, mpmul(gel(v,i) ,gel(rd,i)));
  s = mpneg(mpmul(beta, s));
  for (i=1; i<l; i++) gel(rd,i) = mpadd(gel(rd,i), mpmul(s, gel(v,i)));
}

static GEN
ApplyAllQ(GEN Q, GEN r0, long k)
{
  pari_sp av = avma;
  GEN r = shallowcopy(r0);
  long j;
  for (j=1; j<k; j++) ApplyQ(gel(Q,j), r);
  return gerepilecopy(av, r);
}

/* compute B[k] = | x[k] |^2, update mu(k, 1..k-1) using Householder matrices
 * (Q = Householder(x[1..k-1]) in factored form) */
static int
incrementalQ(GEN x, GEN L, GEN B, GEN Q, long k, long prec)
{
  GEN r = ApplyAllQ(Q, gel(x,k), k);
  return FindApplyQ(r, L, B, k, Q, prec);
}

/* Q vector of Householder matrices orthogonalizing x[1..j0].
 * Q[i] = 0 means not computed yet */
static int
Householder_get_mu(GEN x, GEN L, GEN B, long k, GEN Q, long prec)
{
  GEN Nx, invNx, m;
  long i, j, j0;
  if (!Q) Q = zerovec(k);
  for (j=1; j<=k; j++)
    if (typ(Q[j]) == t_INT) break;
  j0 = j;
  for (   ; j<=k; j++)
    if (! incrementalQ(x, L, B, Q, j, prec)) return 0;
  for (j=1; j<k; j++)
  {
    m = gel(L,j); Nx = gel(m,j); /* should set gel(m,j) = gen_1; but need it later */
    invNx = ginv(Nx);
    for (i=max(j0, j+1); i<=k; i++) gel(m,i) = mpmul(invNx, gel(m,i));
  }
  return 1;
}

GEN
sqred1_from_QR(GEN x, long prec)
{
  long j, k = lg(x)-1;
  GEN L, B = zerovec(k);
  L = cgetg(k+1, t_MAT);
  for (j=1; j<=k; j++) gel(L,j) = zerocol(k);
  if (!Householder_get_mu(x, L, B, k, NULL, prec)) return NULL;
  for (j=1; j<=k; j++) coeff(L,j,j) = B[j];
  return shallowtrans(L);
}

GEN
R_from_QR(GEN x, long prec)
{
  long j, k = lg(x)-1;
  GEN L, B = zerovec(k), Q = cgetg(k+1, t_VEC);
  L = cgetg(k+1, t_MAT);
  for (j=1; j<=k; j++) gel(L,j) = zerocol(k);
  for (j=1; j<=k; j++)
    if (!incrementalQ(x, L, B, Q, j, prec)) return NULL;
  return shallowtrans(L);
}

/********************************************************************/
/**             QR Factorization via Gram-Schmidt                  **/
/********************************************************************/

/* compute B[k] = | x[k] |^2, update mu(k, 1..k-1).
 * Classical Gram-Schmidt (unstable!) */
static int
incrementalGS(GEN x, GEN mu, GEN B, long k)
{
  GEN s,A = cgetg(k+1, t_COL); /* scratch space */
  long i, j;
  pari_sp av;
  A[1] = coeff(x,k,1);
  for(j=1;j<k;)
  {
    gcoeff(mu,k,j) = mpdiv(gel(A,j), gel(B,j));
    j++; av = avma;
    /* A[j] <-- x[k,j] - sum_{i<j} mu[j,i] A[i] */
    s = mpmul(gcoeff(mu,j,1),gel(A,1));
    for (i=2; i<j; i++) s = mpadd(s, mpmul(gcoeff(mu,j,i),gel(A,i)));
    s = mpneg(s); gel(A,j) = gerepileuptoleaf(av, mpadd(gcoeff(x,k,j), s));
  }
  B[k] = A[k]; return (signe(gel(B,k)) > 0 && no_prec_pb(gel(B,k)));
}

#if 0
/* return Gram-Schmidt orthogonal basis (f) associated to (e), B is the
 * vector of the (f_i . f_i)*/
GEN
gram_schmidt(GEN e, GEN *ptB)
{
  long i,j,lx = lg(e);
  GEN f = shallowcopy(e), B, iB;

  B = cgetg(lx, t_VEC);
  iB= cgetg(lx, t_VEC);

  for (i=1;i<lx;i++)
  {
    GEN p1 = NULL;
    pari_sp av;
    gel(B,i) = sqscal(gel(f,i));
    gel(iB,i) = ginv(gel(B,i)); av = avma;
    for (j=1; j<i; j++)
    {
      GEN mu = gmul(gscal(gel(e,i),gel(f,j)), gel(iB,j));
      GEN p2 = gmul(mu, gel(f,j));
      p1 = p1? gadd(p1,p2): p2;
    }
    p1 = p1? gerepileupto(av, gsub(gel(e,i), p1)): gel(e,i);
    gel(f,i) = p1;
  }
  *ptB = B; return f;
}
#endif

/********************************************************************/
/**                                                                **/
/**                          LLL ALGORITHM                         **/
/**                                                                **/
/********************************************************************/
static GEN
lll_trivial(GEN x, long flag)
{
  GEN y;
  if (lg(x) == 1)
  { /* dim x = 0 */
    if (! (flag & lll_ALL)) return cgetg(1,t_MAT);
    y=cgetg(3,t_VEC);
    gel(y,1) = cgetg(1,t_MAT);
    gel(y,2) = cgetg(1,t_MAT); return y;
  }
  /* here dim = 1 */
  if (gcmp0(gel(x,1)))
  {
    switch(flag & (~lll_GRAM))
    {
      case lll_KER: return matid(1);
      case lll_IM : return cgetg(1,t_MAT);
      default: y=cgetg(3,t_VEC);
        gel(y,1) = matid(1);
        gel(y,2) = cgetg(1,t_MAT); return y;
    }
  }
  if (flag & lll_GRAM) flag ^= lll_GRAM; else x = NULL;
  switch (flag)
  {
    case lll_KER: return cgetg(1,t_MAT);
    case lll_IM : return matid(1);
    default: y=cgetg(3,t_VEC);
      gel(y,1) = cgetg(1,t_MAT);
      gel(y,2) = x? gcopy(x): matid(1); return y;
  }
}

static GEN
lll_finish(GEN h,GEN fl,long flag)
{
  long i, k, l = lg(fl);
  GEN g;

  k=1; while (k<l && !fl[k]) k++;
  switch(flag & (~lll_GRAM))
  {
    case lll_KER: setlg(h,k);
      return h;

    case lll_IM: h += k-1; h[0] = evaltyp(t_MAT) | evallg(l-k+1);
      return h;
  }
  g = cgetg(k, t_MAT); for (i=1; i<k; i++) g[i] = h[i];
  h += k-1; h[0] = evaltyp(t_MAT) | evallg(l-k+1);
  return mkvec2(g, h);
}

/* h[,k] += q * h[,l]. Inefficient if q = 0 */
static void
Zupdate_col(long k, long l, GEN q, long K, GEN h)
{
  GEN *hl, *hk;
  long i, qq = itos_or_0(q);

  if (!h) return;
  hl = (GEN*)h[l]; hk = (GEN*)h[k];
  if (!qq) {
    for (i=1;i<=K;i++) if (signe(hl[i])) hk[i] = addii(hk[i],mulii(q,hl[i]));
    return;
  }
  if (qq == 1) {
    for (i=1;i<=K;i++) { if (signe(hl[i])) hk[i] = addii(hk[i],hl[i]); }
  } else if (qq == -1) {
    for (i=1;i<=K;i++) { if (signe(hl[i])) hk[i] = subii(hk[i],hl[i]); }
  } else {
    for (i=1;i<=K;i++) if (signe(hl[i])) hk[i] = addii(hk[i],mulsi(qq,hl[i]));
  }
}

/* L[k,] += q * L[l,], l < k. Inefficient if q = 0 */
static void
Zupdate_row(long k, long l, GEN q, GEN L, GEN B)
{
  long i, qq = itos_or_0(q);
  if (!qq)
  {
    for(i=1;i<l;i++)  gcoeff(L,k,i) = addii(gcoeff(L,k,i),mulii(q,gcoeff(L,l,i)));
    gcoeff(L,k,l) = addii(gcoeff(L,k,l), mulii(q,B));
    return;
  }
  if (qq == 1) {
    for (i=1;i<l; i++) gcoeff(L,k,i) = addii(gcoeff(L,k,i),gcoeff(L,l,i));
    gcoeff(L,k,l) = addii(gcoeff(L,k,l), B);
  } else if (qq == -1) {
    for (i=1;i<l; i++) gcoeff(L,k,i) = subii(gcoeff(L,k,i),gcoeff(L,l,i));
    gcoeff(L,k,l) = addii(gcoeff(L,k,l), negi(B));
  } else {
    for(i=1;i<l;i++) gcoeff(L,k,i) = addii(gcoeff(L,k,i),mulsi(qq,gcoeff(L,l,i)));
    gcoeff(L,k,l) = addii(gcoeff(L,k,l), mulsi(qq,B));
  }
}

static void
update_row(long k, long l, GEN q, GEN L)
{
  long i;
  if (is_pm1(q))
  {
    if (signe(q) > 0)
    {
      for (i=1;i<l; i++) gcoeff(L,k,i) = mpadd(gcoeff(L,k,i),gcoeff(L,l,i));
    } else {
      for (i=1;i<l; i++) gcoeff(L,k,i) = mpsub(gcoeff(L,k,i),gcoeff(L,l,i));
    }
  } else {
    for(i=1;i<l;i++)  gcoeff(L,k,i) = mpadd(gcoeff(L,k,i),mpmul(q,gcoeff(L,l,i)));
  }
  gcoeff(L,k,l) = mpadd(gcoeff(L,k,l), q);
}

static void
ZRED_gram(long k, long l, GEN x, GEN h, GEN L, GEN B, long K)
{
  long i,lx;
  GEN q = truedivii(addii(B,shifti(gcoeff(L,k,l),1)), shifti(B,1));
  GEN xk,xl;
  if (!signe(q)) return;
  q = negi(q);
  xl = gel(x,l); xk = gel(x,k);
  lx = lg(xl);
  if (is_pm1(q))
  {
    if (signe(q) > 0)
    {
      gel(xk,k) = addii(gel(xk,k), gel(xl,k));
      for (i=1;i<lx;i++)
        gcoeff(x,k,i) = gel(xk,i) = addii(gel(xk,i), gel(xl,i));
    } else {
      gel(xk,k) = subii(gel(xk,k), gel(xl,k));
      for (i=1;i<lx;i++)
        gcoeff(x,k,i) = gel(xk,i) = subii(gel(xk,i), gel(xl,i));
    }
  } else { /* h[,k] += q* h[,l]. x[,k] += q * x[,l]. L[k,] += q* L[l,] */
    gel(xk,k) = addii(gel(xk,k), mulii(q,gel(xl,k)));
    for(i=1;i<lx;i++)
      gcoeff(x,k,i)=gel(xk,i) = addii(gel(xk,i),mulii(q,gel(xl,i)));
  }
  Zupdate_row(k,l,q,L,B);
  Zupdate_col(k,l,q,K,h);
}

static void
ZRED(long k, long l, GEN x, GEN h, GEN L, GEN B, long K)
{
  GEN q = truedivii(addii(B,shifti(gcoeff(L,k,l),1)), shifti(B,1));
  if (!signe(q)) return;
  q = negi(q);
  Zupdate_row(k,l,q,L,B);
  Zupdate_col(k,l,q,K,h);
  gel(x,k) = ZV_lincomb(gen_1, q, gel(x,k), gel(x,l));
}

static GEN
round_safe(GEN q)
{
  if (typ(q) == t_INT) return q;
  if (expo(q) > 30)
  {
    long e;
    q = grndtoi(q, &e);
    if (e > 0) return NULL;
  } else
    q = ground(q);
  return q;
}

/* b[k] <-- b[k] - round(L[k,l]) b[l], only b[1] ... b[K] modified so far
 * assume l < k and update x = Gram(b), L = Gram Schmidt coeffs. */
static int
RED_gram(long k, long l, GEN x, GEN h, GEN L, long K)
{
  long i,lx;
  GEN q = round_safe(gcoeff(L,k,l));
  GEN xk, xl;

  if (!q) return 0;
  if (!signe(q)) return 1;
  q = negi(q); lx = lg(x);
  xk = gel(x,k); xl = gel(x,l);
  if (is_pm1(q))
  {
    if (signe(q) > 0)
    {
      gel(xk,k) = mpadd(gel(xk,k), gel(xl,k));
      for (i=1;i<lx;i++) gcoeff(x,k,i)=gel(xk,i) = mpadd(gel(xk,i), gel(xl,i));
    } else {
      gel(xk,k) = mpsub(gel(xk,k), gel(xl,k));
      for (i=1;i<lx;i++) gcoeff(x,k,i)=gel(xk,i) = mpsub(gel(xk,i), gel(xl,i));
    }
  } else {
    gel(xk,k) = mpadd(gel(xk,k), mpmul(q,gel(xl,k)));
    for (i=1;i<lx;i++)
      gcoeff(x,k,i)=gel(xk,i) = mpadd(gel(xk,i), mpmul(q,gel(xl,i)));
  }
  update_row(k,l,q,L);
  Zupdate_col(k,l,q,K,h); return 1;
}

static void
update_col(long k, long l, GEN q, GEN x)
{
  long i,lx;
  GEN xk = gel(x,k), xl = gel(x,l);
  lx = lg(xk);
  if (is_pm1(q))
  {
    if (signe(q) > 0)
    {
      for (i=1;i<lx;i++) gel(xk,i) = mpadd(gel(xk,i), gel(xl,i));
    } else {
      for (i=1;i<lx;i++) gel(xk,i) = mpsub(gel(xk,i), gel(xl,i));
    }
  } else {
    for (i=1;i<lx;i++) gel(xk,i) = mpadd(gel(xk,i), mpmul(q,gel(xl,i)));
  }
}

static int
RED(long k, long l, GEN x, GEN h, GEN L, long K)
{
  GEN q = round_safe(gcoeff(L,k,l));
  if (!q) return 0;
  if (!signe(q)) return 1;
  q = negi(q);
  update_col(k,l,q,x);
  update_row(k,l,q,L);
  Zupdate_col(k,l,q,K,h); return 1;
}

static int
do_ZSWAP(GEN x, GEN h, GEN L, GEN B, long kmax, long k, long D, GEN fl,
         int gram)
{
  GEN la,la2,p1,Bk;
  long i, j, lx;
  pari_sp av;

  if (!fl[k-1]) return 0;
  av = avma;
  la = gcoeff(L,k,k-1); la2 = sqri(la);
  Bk = gel(B,k);
  if (fl[k])
  {
    GEN q;
    if (!D) return 0; /* only lswap non-kernel + kernel vector */
    q = addii(la2, mulii(gel(B,k-1),gel(B,k+1)));
    if (cmpii(mulsi(D-1,sqri(Bk)), mulsi(D,q)) <= 0) {
      avma = av; return 0;
    }
    gel(B,k) = diviiexact(q, Bk);
  }
  /* ZSWAP(k-1,k) */
  if (DEBUGLEVEL>3 && k==kmax)
  { /* output diagnostics associated to re-normalized rational quantities */
    pari_sp av1 = avma;
    GEN d = mulii(gel(B,k-1),gel(B,k+1));
    p1 = subii(mulsi(D-1, sqri(Bk)), mulsi(D, la2));
    fprintferr(" (%ld)", expi(p1) - expi(mulsi(D, d)));
    avma = av1;
  }
  if (h) lswap(h[k-1], h[k]);
  lswap(x[k-1], x[k]); lx = lg(x);
  if (gram)
    for (j=1; j < lx; j++) lswap(coeff(x,k-1,j), coeff(x,k,j));
  for (j=1; j<k-1; j++) lswap(coeff(L,k-1,j), coeff(L,k,j))
  if (fl[k])
  {
    av = avma;
    for (i=k+1; i<=kmax; i++)
    {
      GEN t = gcoeff(L,i,k);
      p1 = subii(mulii(gel(B,k+1),gcoeff(L,i,k-1)), mulii(la,t));
      p1 = diviiexact(p1, Bk);
      gcoeff(L,i,k) = icopy_av(p1,(GEN)av);
      av = avma = (pari_sp)coeff(L,i,k);

      p1 = addii(mulii(la,gcoeff(L,i,k-1)), mulii(gel(B,k-1),t));
      p1 = diviiexact(p1, Bk);
      gcoeff(L,i,k-1) = icopy_av(p1,(GEN)av);
      av = avma = (pari_sp)coeff(L,i,k-1);
    }
  }
  else if (signe(la))
  {
    p1 = diviiexact(la2, Bk);
    gel(B,k+1) = gel(B,k) = p1;
    for (i=k+2; i<=lx; i++) gel(B,i) = diviiexact(mulii(p1,gel(B,i)), Bk);
    for (i=k+1; i<=kmax; i++)
      gcoeff(L,i,k-1) = diviiexact(mulii(la,gcoeff(L,i,k-1)), Bk);
    for (j=k+1; j<kmax; j++)
      for (i=j+1; i<=kmax; i++)
        gcoeff(L,i,j) = diviiexact(mulii(p1,gcoeff(L,i,j)), Bk);
  }
  else
  {
    for (i=k+1; i<=kmax; i++)
    {
      gcoeff(L,i,k) = gcoeff(L,i,k-1);
      gcoeff(L,i,k-1) = gen_0;
    }
    B[k] = B[k-1]; fl[k] = 1; fl[k-1] = 0;
  }
  return 1;
}

static int
do_SWAP(GEN x, GEN h, GEN L, GEN B, long kmax, long k, GEN delta, int gram)
{
  GEN la,la2, BK,BB,q;
  long i, j, lx;
  pari_sp av;

  av = avma;
  la = gcoeff(L,k,k-1); la2 = gsqr(la);
  q = mpmul(gel(B,k-1), mpsub(delta,la2));
  if (mpcmp(q, gel(B,k)) <= 0) { avma = av; return 0; }

  /* SWAP(k-1,k) */
  if (DEBUGLEVEL>3 && k==kmax) {
    fprintferr(" (%ld)", gexpo(q) - gexpo(gel(B,k))); flusherr();
  }
  BB = mpadd(gel(B,k), mpmul(gel(B,k-1),la2));
  if (!signe(BB)) { B[k] = 0; return 1; } /* precision pb */

  gcoeff(L,k,k-1) = mpdiv(mpmul(la,gel(B,k-1)), BB);
  BK = mpdiv(gel(B,k), BB);
  gel(B,k) = mpmul(gel(B,k-1), BK);
  gel(B,k-1) = BB;

  if (h) lswap(h[k-1],h[k]);
  lswap(x[k-1],x[k]); lx = lg(x);
  if (gram)
    for (j=1; j < lx; j++) lswap(coeff(x,k-1,j), coeff(x,k,j));
  for (j=1; j<k-1; j++) lswap(coeff(L,k-1,j), coeff(L,k,j))
  for (i=k+1; i<=kmax; i++)
  {
    GEN t = gcoeff(L,i,k);
    gcoeff(L,i,k) = mpsub(gcoeff(L,i,k-1), mpmul(la,t));
    gcoeff(L,i,k-1) = mpadd(t, mpmul(gcoeff(L,k,k-1), gcoeff(L,i,k)));
  }
  return 1;
}
static void
ZincrementalGS(GEN x, GEN L, GEN B, long k, GEN fl, int gram)
{
  GEN u = NULL; /* gcc -Wall */
  long i, j, s;
  for (j=1; j<=k; j++)
    if (j==k || fl[j])
    {
      pari_sp av = avma;
      u = gram? gcoeff(x,k,j): gscali(gel(x,k), gel(x,j));
      for (i=1; i<j; i++)
        if (fl[i])
        {
          u = subii(mulii(gel(B,i+1), u), mulii(gcoeff(L,k,i), gcoeff(L,j,i)));
          u = diviiexact(u, gel(B,i));
        }
      gcoeff(L,k,j) = gerepileuptoint(av, u);
    }
  s = signe(u);
  if (s == 0) B[k+1] = B[k];
  else
  {
    if (s < 0) pari_err(lllger3);
    B[k+1] = coeff(L,k,k); gcoeff(L,k,k) = gen_1; fl[k] = 1;
  }
}

/* x integer matrix. Beware: this function can return NULL */
GEN
lllint_marked(long *pMARKED, GEN x, long D, int gram,
              GEN *pth, GEN *ptfl, GEN *ptB)
{
  long lx = lg(x), hx, i, j, k, l, n, kmax, MARKED;
  pari_sp av, lim;
  GEN B,L,h,fl;

  if (typ(x) != t_MAT) pari_err(typeer,"lllint");
  fl = cgetg(lx,t_VECSMALL);
  if (ptfl) *ptfl = fl;
  n = lx-1; if (n <= 1) return NULL;
  MARKED = pMARKED? *pMARKED: 0;
  hx = lg(x[1]);
  if (gram && hx != lx) pari_err(mattype1,"lllint");

  av = avma; lim = stack_lim(av,1); x = shallowcopy(x);
  B = gscalcol_i(gen_1, lx);
  L = cgetg(lx,t_MAT);
  for (j=1; j<lx; j++)
  {
    for (i=1; i<hx; i++)
      if (typ(gcoeff(x,i,j)) != t_INT) pari_err(typeer,"lllint_marked");
    fl[j] = 0; gel(L,j) = zerocol(n);
  }
  h = pth? matid(n): NULL;
  ZincrementalGS(x, L, B, 1, fl, gram);
  kmax = 1;
  if (DEBUGLEVEL>5) fprintferr("k = ");
  for (k=2;;)
  {
    if (k > kmax)
    {
      if (DEBUGLEVEL>3) fprintferr("K%ld ",k);
      ZincrementalGS(x, L, B, k, fl, gram);
      kmax = k;
    }
    if (k != MARKED)
    {
      if (!gram) ZRED(k,k-1, x,h,L,gel(B,k),kmax);
      else  ZRED_gram(k,k-1, x,h,L,gel(B,k),kmax);
    }
    if (do_ZSWAP(x,h,L,B,kmax,k,D,fl,gram))
    {
      if      (MARKED == k)   MARKED = k-1;
      else if (MARKED == k-1) MARKED = k;
      if (k > 2) k--;
    }
    else
    {
      if (k != MARKED)
        for (l=k-2; l; l--)
        {
          if (!gram) ZRED(k,l, x,h,L,gel(B,l+1),kmax);
          else  ZRED_gram(k,l, x,h,L,gel(B,l+1),kmax);
          if (low_stack(lim, stack_lim(av,1)))
          {
            if(DEBUGMEM>1) pari_warn(warnmem,"lllint[1], kmax = %ld", kmax);
            gerepileall(av,h?4:3,&B,&L,&x,&h);
          }
        }
      if (++k > n) break;
    }
    if (low_stack(lim, stack_lim(av,1)))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"lllint[2], kmax = %ld", kmax);
      gerepileall(av,h?4:3,&B,&L,&x,&h);
    }
  }
  if (DEBUGLEVEL>3) fprintferr("\n");
  if (ptB)  *ptB  = B;
  if (ptfl) *ptfl = fl;
  if (pth)  *pth = h;
  if (pMARKED) *pMARKED = MARKED;
  return h? h: x;
}

/* Beware: this function can return NULL (dim x <= 1) */
GEN
lllint_i(GEN x, long D, int gram, GEN *pth, GEN *ptfl, GEN *ptB)
{
  return lllint_marked(NULL, x,D,gram,pth,ptfl,ptB);
}

/* return x * lllint(x). No garbage collection */
GEN
lllint_ip(GEN x, long D)
{
  GEN fl, h = lllint_i(x, D, 0, NULL, &fl, NULL);
  if (!h) return x;
  return lll_finish(h, fl, lll_IM);
}

GEN
lllall(GEN x, long D, int gram, long flag)
{
  pari_sp av = avma;
  GEN fl, junk, h = lllint_i(x, D, gram, &junk, &fl, NULL);
  if (!h) return lll_trivial(x,flag);
  return gerepilecopy(av, lll_finish(h,fl,flag));
}

GEN
lllint(GEN x) { return lllall(x,LLLDFT,0, lll_IM); }

GEN
lllkerim(GEN x) { return lllall(x,LLLDFT,0, lll_ALL); }

GEN
lllgramint(GEN x) { return lllall(x,LLLDFT,1, lll_IM | lll_GRAM); }

GEN
lllgramkerim(GEN x) { return lllall(x,LLLDFT,1, lll_ALL | lll_GRAM); }

static int
pslg(GEN x)
{
  long tx;
  if (gcmp0(x)) return 2;
  tx = typ(x); return is_scalar_t(tx)? 3: lg(x);
}

static GEN
to_MP(GEN x, long prec)
{ return (typ(x) == t_INT && !signe(x))? gen_0: gtofp(x, prec); }
static GEN
col_to_MP(GEN x, long prec)
{
  long j, l = lg(x);
  GEN y = cgetg(l, t_COL);
  for (j=1; j<l; j++) gel(y,j) = to_MP(gel(x,j), prec);
  return y;
}
static GEN
mat_to_MP(GEN x, long prec)
{
  long j, l = lg(x);
  GEN y;
  if (typ(x) != t_MAT) return col_to_MP(x, prec);
  y = cgetg(l, t_MAT);
  for (j=1; j<l; j++) gel(y,j) = col_to_MP(gel(x,j), prec);
  return y;
}

static GEN
to_mp(GEN x, long prec)
{ return (typ(x) == t_INT)? x: gtofp(x, prec); }
static GEN
col_to_mp(GEN x, long prec)
{
  long j, l = lg(x);
  GEN y = cgetg(l, t_COL);
  for (j=1; j<l; j++) gel(y,j) = to_mp(gel(x,j), prec);
  return y;
}
static GEN
mat_to_mp(GEN x, long prec)
{
  long j, l = lg(x);
  GEN y = cgetg(l, t_MAT);
  for (j=1; j<l; j++) gel(y,j) = col_to_mp(gel(x,j), prec);
  return y;
}

static int
REDgen(long k, long l, GEN h, GEN L, GEN B)
{
  GEN q, u = gcoeff(L,k,l);
  long i;

  if (pslg(u) < pslg(B)) return 0;

  q = gneg(gdeuc(u,B));
  gel(h,k) = gadd(gel(h,k), gmul(q,gel(h,l)));
  for (i=1; i<l; i++) gcoeff(L,k,i) = gadd(gcoeff(L,k,i), gmul(q,gcoeff(L,l,i)));
  gcoeff(L,k,l) = gadd(gcoeff(L,k,l), gmul(q,B)); return 1;
}

static int
do_SWAPgen(GEN h, GEN L, GEN B, long k, GEN fl, int *flc)
{
  GEN p1, la, la2, Bk;
  long ps1, ps2, i, j, lx;

  if (!fl[k-1]) return 0;

  la = gcoeff(L,k,k-1); la2 = gsqr(la);
  Bk = gel(B,k);
  if (fl[k])
  {
    GEN q = gadd(la2, gmul(gel(B,k-1),gel(B,k+1)));
    ps1 = pslg(gsqr(Bk));
    ps2 = pslg(q);
    if (ps1 <= ps2 && (ps1 < ps2 || !*flc)) return 0;
    *flc = (ps1 != ps2);
    gel(B,k) = gdiv(q, Bk);
  }

  lswap(h[k-1], h[k]); lx = lg(L);
  for (j=1; j<k-1; j++) lswap(coeff(L,k-1,j), coeff(L,k,j));
  if (fl[k])
  {
    for (i=k+1; i<lx; i++)
    {
      GEN t = gcoeff(L,i,k);
      p1 = gsub(gmul(gel(B,k+1),gcoeff(L,i,k-1)), gmul(la,t));
      gcoeff(L,i,k) = gdiv(p1, Bk);
      p1 = gadd(gmul(la,gcoeff(L,i,k-1)), gmul(gel(B,k-1),t));
      gcoeff(L,i,k-1) = gdiv(p1, Bk);
    }
  }
  else if (!gcmp0(la))
  {
    p1 = gdiv(la2, Bk);
    gel(B,k+1) = gel(B,k) = p1;
    for (i=k+2; i<=lx; i++) gel(B,i) = gdiv(gmul(p1,gel(B,i)),Bk);
    for (i=k+1; i<lx; i++)
      gcoeff(L,i,k-1) = gdiv(gmul(la,gcoeff(L,i,k-1)), Bk);
    for (j=k+1; j<lx-1; j++)
      for (i=j+1; i<lx; i++)
        gcoeff(L,i,j) = gdiv(gmul(p1,gcoeff(L,i,j)), Bk);
  }
  else
  {
    gcoeff(L,k,k-1) = gen_0;
    for (i=k+1; i<lx; i++)
    {
      gcoeff(L,i,k) = gcoeff(L,i,k-1);
      gcoeff(L,i,k-1) = gen_0;
    }
    B[k] = B[k-1]; fl[k] = 1; fl[k-1] = 0;
  }
  return 1;
}

static void
incrementalGSgen(GEN x, GEN L, GEN B, long k, GEN fl)
{
  GEN u = NULL; /* gcc -Wall */
  long i, j, tu;
  for (j=1; j<=k; j++)
    if (j==k || fl[j])
    {
      u = gcoeff(x,k,j); tu = typ(u);
      if (! is_extscalar_t(tu)) pari_err(typeer,"incrementalGSgen");
      for (i=1; i<j; i++)
        if (fl[i])
        {
          u = gsub(gmul(gel(B,i+1),u), gmul(gcoeff(L,k,i),gcoeff(L,j,i)));
          u = gdiv(u, gel(B,i));
        }
      gcoeff(L,k,j) = u;
    }
  if (gcmp0(u)) B[k+1] = B[k];
  else
  {
    B[k+1] = coeff(L,k,k); gcoeff(L,k,k) = gen_1; fl[k] = 1;
  }
}

static GEN
lllgramallgen(GEN x, long flag)
{
  long lx = lg(x), i, j, k, l, n;
  pari_sp av0 = avma, av, lim;
  GEN B, L, h, fl;
  int flc;

  if (typ(x) != t_MAT) pari_err(typeer,"lllgramallgen");
  n = lx-1; if (n<=1) return lll_trivial(x,flag);
  if (lg(x[1]) != lx) pari_err(mattype1,"lllgramallgen");

  fl = cgetg(lx, t_VECSMALL);

  av = avma; lim = stack_lim(av,1);
  B = gscalcol_i(gen_1, lx);
  L = cgetg(lx,t_MAT);
  for (j=1; j<lx; j++) { gel(L,j) = zerocol(n); fl[j] = 0; }

  h = matid(n);
  for (i=1; i<lx; i++)
    incrementalGSgen(x, L, B, i, fl);
  flc = 0;
  for(k=2;;)
  {
    if (REDgen(k, k-1, h, L, gel(B,k))) flc = 1;
    if (do_SWAPgen(h, L, B, k, fl, &flc)) { if (k > 2) k--; }
    else
    {
      for (l=k-2; l>=1; l--)
        if (REDgen(k, l, h, L, gel(B,l+1))) flc = 1;
      if (++k > n) break;
    }
    if (low_stack(lim, stack_lim(av,1)))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"lllgramallgen");
      gerepileall(av,3,&B,&L,&h);
    }
  }
  return gerepilecopy(av0, lll_finish(h,fl,flag));
}

static GEN
_mul2n(GEN x, long e) { return e? gmul2n(x, e): x; }

/* E = scaling exponents
 * F = backward scaling exponents
 * X = lattice basis, Xs = associated scaled basis
 * Q = vector of Householder transformations
 * h = base change matrix (from X0)
 * R = from QR factorization of Xs[1..k-1] */
static int
HRS(long MARKED, long k, int prim, long kmax, GEN X, GEN Xs, GEN h, GEN R,
    GEN Q, GEN E, GEN F)
{
  long e, i, N = lg(X[k]), rounds = 0;
  const long prec = MEDDEFAULTPREC; /* 128 bits */
  GEN q, tau2, rk;
  int overf;

  E[k] = prim? E[k-1]: 0;
  F[k] = 0;
  gel(Xs,k) = E[k]? gmul2n(gel(X,k), E[k]): shallowcopy(gel(X,k));
  rounds = 0;
  if (k == MARKED) goto DONE; /* no size reduction/scaling */

UP:
  rk = ApplyAllQ(Q, col_to_MP(gel(Xs,k), prec), k);
  overf = 0;
  for (i = k-1; i > 0; i--)
  { /* size seduction of Xs[k] */
    GEN q2, mu, muint;
    long ex;
    pari_sp av = avma;

    mu = mpdiv(gel(rk,i), gcoeff(R,i,i));
    if (!signe(mu)) { avma = av; continue; }
    mu = _mul2n(mu, -F[i]); /*backward scaling*/
    ex = gexpo(mu);
    if (ex > 10) {
      overf = 1; /* large size reduction */
      muint = ex > 30? ceil_safe(mu): ground(mu);
    }
    else if (fabs(gtodouble(mu)) > 0.51)
      muint = ground(mu);
    else { avma = av; continue; }

    q = gerepileuptoint(av, negi(muint));
    q2= _mul2n(q, F[i]); /*backward scaling*/
    Zupdate_col(k, i, q2, N-1, Xs);
    rk = gadd(rk, gmul(q2, gel(R,i)));
    /* if in HRS', propagate the last SR on X (for SWAP test) */
    if (prim && (i == k-1)) {
      Zupdate_col(k, i, q, kmax, h);
       update_col(k, i, q, X);
    }
  }
  /* If large size-reduction performed, restart from exact data */
  if (overf)
  {
    if (++rounds > 100) return 0; /* detect infinite loop */
    goto UP;
  }

  if (prim || k == 1) goto DONE; /* DONE */

  rounds = 0;
  /* rescale Xs[k] ? */
  tau2 = signe(rk[k])? gsqr(gel(rk,k)): gen_0;
  for (i = k+1; i < N; i++)
    if (signe(rk[i])) tau2 = mpadd(tau2, gsqr(gel(rk,i)));
  q = mpdiv(gsqr(gcoeff(R,1,1)), tau2);
  e = gexpo(q)/2 + F[1]; /* ~ log_2 (||Xs[1]|| / tau), backward scaling*/
  if (e > 0)
  { /* rescale */
    if (e > 30) e = 30;
    gel(Xs,k) = gmul2n(gel(Xs,k), e);
    E[k] += e; goto UP; /* one more round */
  }
  else if (e < 0)
  { /* enable backward scaling */
    for (i = 1; i < k; i++) F[i] -= e;
  }

DONE:
  rk = ApplyAllQ(Q, col_to_MP(gel(Xs,k), prec), k);
  (void)FindApplyQ(rk, R, NULL, k, Q, prec);
  return 1;
}

static GEN
rescale_to_int(GEN x)
{
  long e, prec = gprecision(x);
  GEN y;
  if (!prec) return Q_primpart(x);

  y = gmul2n(x, bit_accuracy(prec) - gexpo(x));
  return gcvtoi(y, &e);
}

GEN
lll_scaled(long MARKED, GEN X0, long D)
{
  GEN delta, X, Xs, h, R, Q, E, F;
  long j, kmax = 1, k, N = lg(X0);
  pari_sp lim, av, av0 = avma;
  int retry = 0;

  delta = stor(D-1, DEFAULTPREC);
  delta = divrs(delta,D);
  E  = const_vecsmall(N-1, 0);
  F  = const_vecsmall(N-1, 0);

  av = avma; lim = stack_lim(av, 1);
  h = matid(N-1);
  X = rescale_to_int(X0);

PRECPB:
  k = 1;
  if (retry++) return mkvec(h);
  Q  = zerovec(N-1);
  Xs = zeromat(N-1, N-1);
  R  = cgetg(N, t_MAT);
  for (j=1; j<N; j++) gel(R,j) = zerocol(N-1);

  while (k < N)
  {
    pari_sp av1;
    if (k == 1) {
      HRS(MARKED, 1, 0, kmax, X, Xs, h, R, Q, E, F);
      k = 2;
    }
    if (k > kmax) {
      kmax = k;
      if (DEBUGLEVEL>3) {fprintferr(" K%ld",k);flusherr();}
    }
    if (!HRS(MARKED, k, 1, kmax, X, Xs, h, R, Q, E, F)) goto PRECPB;

    av1 = avma;
    if (mpcmp( mpmul(delta, gsqr(gcoeff(R, k-1,k-1))),
               mpadd(gsqr(gcoeff(R,k-1,k)), gsqr(gcoeff(R, k,k)))) > 0)
    {
      if (DEBUGLEVEL>3 && k == kmax)
      {
        GEN q = mpsub( mpmul(delta, gsqr(gcoeff(R, k-1,k-1))),
                       gsqr(gcoeff(R,k-1,k)));
        fprintferr(" (%ld)", gexpo(q) - gexpo(gsqr(gcoeff(R, k,k))));
      }
      lswap(X[k], X[k-1]);
      lswap(h[k], h[k-1]);
      if      (MARKED == k)   MARKED = k-1;
      else if (MARKED == k-1) MARKED = k;
      avma = av1; k--;
    }
    else
    {
      avma = av1;
      if (!HRS(MARKED, k, 0, kmax, X, Xs, h, R, Q, E, F)) goto PRECPB;
      k++;
    }
    if (low_stack(lim, stack_lim(av,1)))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"lllfp[1]");
      gerepileall(av,5,&X,&Xs, &R,&h,&Q);
    }
  }
  return gerepilecopy(av0, h);
}

static long
good_prec(GEN x, long kmax)
{
  long prec = ((kmax<<2) + gexpo(gel(x,kmax))) >> TWOPOTBITS_IN_LONG;
  if (prec < DEFAULTPREC) prec = DEFAULTPREC;
  return prec;
}

/* If gram = 1, x = Gram(b_i), x = (b_i) otherwise
 * Quality ratio = delta = (D-1)/D. Suggested values: D = 4 or D = 100
 *
 *   if (flag = 1): if precision problems, return NULL
 *   if (flag = 2): if precision problems, return mkvec ( h )
 *                  [ partial transformation matrix ], unless we could not
 *                  even start; return NULL in this case.
 *   if (flag = 3): assume x exact and !gram; return LLL-reduced basis, not
 *                  base change
 *
 * If MARKED != 0 make sure e[MARKED] is the first vector of the output basis
 * (which may then not be LLL-reduced) */
GEN
lllfp_marked(long *pMARKED, GEN x, long D, long flag, long prec, int gram)
{
  GEN xinit,L,h,B,L1,delta, Q, H = NULL;
  long retry = 2, lx = lg(x), hx, l, i, j, k, k1, n, kmax, KMAX, MARKED;
  long count, count_max = 8;
  pari_sp av0 = avma, av, lim;
  int isexact, exact_can_leave;
  const int in_place = (flag == 3);

  if (typ(x) != t_MAT) pari_err(typeer,"lllfp");
  n = lx-1; if (n <= 1) return matid(n);
#if 0 /* doesn't work yet */
  return lll_scaled(MARKED, x, D);
#endif

  hx = lg(x[1]);
  if (hx != lx)
  {
    if (gram) pari_err(mattype1,"lllfp");
    if (lx > hx) pari_err(talker,"dependent vectors in lllfp");
  }
  delta = divrs(stor(D-1, DEFAULTPREC), D);
  xinit = x;
  av = avma; lim = stack_lim(av,1);
  if (gram) {
    for (k=2,j=1; j<lx; j++)
    {
      GEN p1=gel(x,j);
      for (i=1; i<=j; i++)
        if (typ(p1[i]) == t_REAL) { l = lg(p1[i]); if (l>k) k=l; }
    }
  } else {
    for (k=2,j=1; j<lx; j++)
    {
      GEN p1=gel(x,j);
      for (i=1; i<hx; i++)
        if (typ(p1[i]) == t_REAL) { l = lg(p1[i]); if (l>k) k=l; }
    }
  }
  if (k == 2)
  {
    if (!prec) return lllint_marked(pMARKED, x, D, gram, &h, NULL, NULL);
    x = mat_to_MP(x, prec);
    isexact = 1;
  }
  else
  {
    if (prec < k) prec = k;
    x = mat_to_mp(x, prec+1);
    isexact = 0;
  }
 /* kmax = maximum column index attained during this run
  * KMAX = same over all runs (after PRECPB) */
  MARKED = pMARKED? *pMARKED: 0;
  kmax = KMAX = 1;
  h = matid(n);

#ifdef LONG_IS_64BIT
#  define PREC_THRESHOLD 32
#  define PREC_DEC_THRESHOLD 7
#else
#  define PREC_THRESHOLD 62
#  define PREC_DEC_THRESHOLD 12
#endif
PRECPB:
  switch(retry--)
  {
    case 2: break; /* entry */
    case 1:
      if (DEBUGLEVEL>3) fprintferr("\n");
      if (flag == 2) return mkvec(h);
      if (isexact || (gram && kmax > 2))
      { /* some progress but precision loss, try again */
        if (prec < PREC_THRESHOLD)
          prec = (prec<<1)-2;
        else
          prec = (long)((prec-2) * 1.25 + 2);
        if (isexact)
        {
          if (DEBUGLEVEL>2) pari_warn(warnprec,"lllfp (exact)",prec);
          if (!in_place) H = H? gmul(H, h): h;
          xinit = gram? qf_base_change(xinit, h, 1): gmul(xinit, h);
          gerepileall(av, in_place? 1: 2, &xinit, &H);
          x = mat_to_MP(xinit, prec);
          h = matid(n);
          retry = 1; /* never abort if x is exact */
          count_max = min(count_max << 1, 512);
          if (DEBUGLEVEL>3) fprintferr("count_max = %ld\n", count_max);
        }
        else
        {
          if (DEBUGLEVEL) pari_warn(warnprec,"lllfp",prec);
          x = gprec_w(xinit,prec);
          x = gram? qf_base_change(x, h, 1): gmul(x, h);
          gerepileall(av, 2, &h, &x);
        }
        kmax = 1; break;
      } /* fall through */
    case 0: /* give up */
      if (DEBUGLEVEL>3) fprintferr("\n");
      if (DEBUGLEVEL) pari_warn(warner,"lllfp giving up");
      if (flag) { avma=av; return NULL; }
      pari_err(lllger3);
  }
  exact_can_leave = 1;
  count = 0;
  Q = zerovec(n);
  L = cgetg(lx,t_MAT);
  B = cgetg(lx,t_COL);
  for (j=1; j<lx; j++) { gel(L,j) = zerocol(n); gel(B,j) = gen_0; }
  if (gram && !incrementalGS(x, L, B, 1))
  {
    if (flag) return NULL;
    pari_err(lllger3);
  }
  if (DEBUGLEVEL>5) fprintferr("k =");
  for(k=2;;)
  {
    if (k > kmax)
    {
      kmax = k;
      if (KMAX < kmax) { KMAX = kmax; count_max = 8; }
      if (DEBUGLEVEL>3) {fprintferr("K%ld ",k);flusherr();}
      if (gram) j = incrementalGS(x, L, B, k);
      else      j = Householder_get_mu(x, L, B, k, Q, prec);
      if (!j) goto PRECPB;
      count = 0;
    }
    else if (isexact && prec > PREC_DEC_THRESHOLD && k == kmax-1
                                                  && ++count > count_max)
    { /* try to reduce precision */
      count = 0;
      prec = (prec+2) >> 1;
      if (DEBUGLEVEL>3) fprintferr("\n...LLL reducing precision to %ld\n",prec);
      if (!in_place) H = H? gmul(H, h): h;
      xinit = gram? qf_base_change(xinit, h, 1): gmul(xinit, h);
      gerepileall(av, in_place? 4: 5,&B,&L,&Q,&xinit, &H);
      x = mat_to_MP(xinit, prec);
      h = matid(n);
    }
    else if (DEBUGLEVEL>5) fprintferr(" %ld",k);
    L1 = gcoeff(L,k,k-1);
    if (typ(L1) == t_REAL && expo(L1) + 20 > bit_accuracy(lg(L1)))
    {
      if (!gram) goto PRECPB;
      if (DEBUGLEVEL>3)
	fprintferr("\nRecomputing Gram-Schmidt, kmax = %ld\n", kmax);
      for (k1=1; k1<=kmax; k1++)
        if (!incrementalGS(x, L, B, k1)) goto PRECPB;
    }
    if (k != MARKED)
    {
      if (!gram) j = RED(k,k-1, x,h,L,KMAX);
      else  j = RED_gram(k,k-1, x,h,L,KMAX);
      if (!j) goto PRECPB;
    }
    if (do_SWAP(x,h,L,B,kmax,k,delta,gram))
    {
      if      (MARKED == k)   MARKED = k-1;
      else if (MARKED == k-1) MARKED = k;
      if (!B[k]) goto PRECPB;
      gel(Q,k) = gel(Q,k-1) = gen_0;
      exact_can_leave = 0;
      if (k>2) k--;
    }
    else
    {
      if (k != MARKED)
        for (l=k-2; l; l--)
        {
          if (!gram) j = RED(k,l, x,h,L,KMAX);
          else  j = RED_gram(k,l, x,h,L,KMAX);
          if (!j) goto PRECPB;
          if (low_stack(lim, stack_lim(av,1)))
          {
            if(DEBUGMEM>1) pari_warn(warnmem,"lllfp[1], kmax = %ld", kmax);
            gerepileall(av, H? 7: 6, &B,&L,&h,&x,&Q,&xinit, &H);
          }
        }
      if (++k > n)
      {
        if (isexact)
        {
          if (exact_can_leave) { if (!in_place && H) h = gmul(H, h); break; }

          if (DEBUGLEVEL>3) fprintferr("\nChecking LLL basis...");
          if (!in_place) H = H? gmul(H, h): h;
          xinit = gram? qf_base_change(xinit, h, 1): gmul(xinit, h);

          prec = good_prec(xinit, kmax);
          if (DEBUGLEVEL>3) fprintferr("in precision %ld\n", prec);
          x = mat_to_MP(xinit, prec);
          h = matid(n);
          exact_can_leave = 1;
          k = 2; kmax = 1; continue;
        }
        else if (!gram && gel(Q,n-1) == gen_0)
        {
          if (DEBUGLEVEL>3) fprintferr("\nChecking LLL basis\n");
          j = Householder_get_mu(gmul(xinit,h), L, B, n, Q, prec);
          if (!j) goto PRECPB;
          k = 2; continue;
        }
        break;
      }
    }
    if (low_stack(lim, stack_lim(av,1)))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"lllfp[2], kmax = %ld", kmax);
      gerepileall(av, H? 7: 6, &B,&L,&h,&x,&Q,&xinit, &H);
    }
  }
  if (in_place) h = gmul(xinit, h);
  if (DEBUGLEVEL>3) fprintferr("\n");
  if (pMARKED) *pMARKED = MARKED;
  return gerepilecopy(av0, h);
}

/* x integral, maximal rank, LLL-reduce in place using fp */
GEN
lllint_fp_ip(GEN x, long D)
{
  return lllfp_marked(NULL, x,D, 3,DEFAULTPREC,0);
}

GEN
lllgramintern(GEN x, long D, long flag, long prec)
{
  return lllfp_marked(NULL, x,D,flag,prec,1);
}

GEN
lllintern(GEN x, long D, long flag, long prec)
{
  return lllfp_marked(NULL, x,D,flag,prec,0);
}

GEN
lllgram(GEN x,long prec) { return lllgramintern(x,LLLDFT,0,prec); }

GEN
lll(GEN x,long prec) { return lllintern(x,LLLDFT,0,prec); }

GEN
qflll0(GEN x, long flag, long prec)
{
  switch(flag)
  {
    case 0: return lll(x,prec);
    case 1: return lllint(x);
    case 2: return lllintpartial(x);
    case 4: return lllkerim(x);
    case 5: return lllkerimgen(x);
    case 8: return lllgen(x);
    default: pari_err(flagerr,"qflll");
  }
  return NULL; /* not reached */
}

GEN
qflllgram0(GEN x, long flag, long prec)
{
  switch(flag)
  {
    case 0: return lllgram(x,prec);
    case 1: return lllgramint(x);
    case 4: return lllgramkerim(x);
    case 5: return lllgramkerimgen(x);
    case 8: return lllgramgen(x);
    default: pari_err(flagerr,"qflllgram");
  }
  return NULL; /* not reached */
}

GEN
gram_matrix(GEN b)
{
  long i,j, lx = lg(b);
  GEN g;
  if (typ(b) != t_MAT) pari_err(typeer,"gram");
  g = cgetg(lx,t_MAT);
  for (i=1; i<lx; i++)
  {
    gel(g,i) = cgetg(lx,t_COL);
    for (j=1; j<=i; j++)
      gcoeff(g,i,j) = gcoeff(g,j,i) = gscal(gel(b,i),gel(b,j));
  }
  return g;
}

GEN
lllgen(GEN x) {
  pari_sp av = avma;
  return gerepileupto(av, lllgramgen(gram_matrix(x)));
}

GEN
lllkerimgen(GEN x) {
  pari_sp av = avma;
  return gerepileupto(av, lllgramallgen(gram_matrix(x),lll_ALL));
}

GEN
lllgramgen(GEN x) { return lllgramallgen(x,lll_IM); }

GEN
lllgramkerimgen(GEN x) { return lllgramallgen(x,lll_ALL); }

/* Def: a matrix M is said to be -partially reduced- if | m1 +- m2 | >= |m1|
 * for any two columns m1 != m2, in M.
 *
 * Input: an integer matrix mat whose columns are linearly independent. Find
 * another matrix T such that mat * T is partially reduced.
 *
 * Output: mat * T if flag = 0;  T if flag != 0,
 *
 * This routine is designed to quickly reduce lattices in which one row
 * is huge compared to the other rows.  For example, when searching for a
 * polynomial of degree 3 with root a mod N, the four input vectors might
 * be the coefficients of
 *     X^3 - (a^3 mod N), X^2 - (a^2 mod N), X - (a mod N), N.
 * All four constant coefficients are O(p) and the rest are O(1). By the
 * pigeon-hole principle, the coefficients of the smallest vector in the
 * lattice are O(p^(1/4)), hence significant reduction of vector lengths
 * can be anticipated.
 *
 * An improved algorithm would look only at the leading digits of dot*.  It
 * would use single-precision calculations as much as possible.
 *
 * Original code: Peter Montgomery (1994) */
GEN
lllintpartialall(GEN m, long flag)
{
  const long ncol = lg(m)-1;
  const pari_sp av = avma;
  GEN tm1, tm2, mid;

  if (typ(m) != t_MAT) pari_err(typeer,"lllintpartial");
  if (ncol <= 1) return flag? matid(ncol): gcopy(m);

  tm1 = flag? matid(ncol): NULL;
  {
    const pari_sp av2 = avma;
    GEN dot11 = sqscali(gel(m,1));
    GEN dot22 = sqscali(gel(m,2));
    GEN dot12 = gscali(gel(m,1), gel(m,2));
    GEN tm  = matid(2); /* For first two columns only */

    int progress = 0;
    long npass2 = 0;

/* Row reduce the first two columns of m. Our best result so far is
 * (first two columns of m)*tm.
 *
 * Initially tm = 2 x 2 identity matrix.
 * Inner products of the reduced matrix are in dot11, dot12, dot22. */
    while (npass2 < 2 || progress)
    {
      GEN dot12new, q = diviiround(dot12, dot22);

      npass2++; progress = signe(q);
      if (progress)
      {/* Conceptually replace (v1, v2) by (v1 - q*v2, v2), where v1 and v2
        * represent the reduced basis for the first two columns of the matrix.
        * We do this by updating tm and the inner products. */
        q = negi(q);
        dot12new = addii(dot12, mulii(q, dot22));
        dot11 = addii(dot11, mulii(q, addii(dot12, dot12new)));
        dot12 = dot12new;
        gel(tm,1) = ZV_lincomb(gen_1,q, gel(tm,1),gel(tm,2));
      }

      /* Interchange the output vectors v1 and v2.  */
      swap(dot11,dot22); lswap(tm[1],tm[2]);

      /* Occasionally (including final pass) do garbage collection.  */
      if (npass2 % 8 == 0 || !progress)
        gerepileall(av2, 4, &dot11,&dot12,&dot22,&tm);
    } /* while npass2 < 2 || progress */

    {
      long i;
      GEN det12 = subii(mulii(dot11, dot22), mulii(dot12, dot12));

      mid = cgetg(ncol+1, t_MAT);
      for (i = 1; i <= 2; i++)
      {
        if (tm1)
        {
          coeff(tm1,1,i) = coeff(tm,1,i);
          coeff(tm1,2,i) = coeff(tm,2,i);
        }
        gel(mid,i) = ZV_lincomb(
           gcoeff(tm,1,i),gcoeff(tm,2,i), gel(m,1),gel(m,2));
      }
      for (i = 3; i <= ncol; i++)
      {
        GEN c = gel(m,i);
	GEN dot1i = gscali(gel(mid,1), c);
        GEN dot2i = gscali(gel(mid,2), c);
       /* ( dot11  dot12 ) (q1)   ( dot1i )
        * ( dot12  dot22 ) (q2) = ( dot2i )
        *
        * Round -q1 and -q2 to nearest integer. Then compute
        *   c - q1*mid[1] - q2*mid[2].
        * This will be approximately orthogonal to the first two vectors in
        * the new basis. */
	GEN q1neg = subii(mulii(dot12, dot2i), mulii(dot22, dot1i));
        GEN q2neg = subii(mulii(dot12, dot1i), mulii(dot11, dot2i));

        q1neg = diviiround(q1neg, det12);
        q2neg = diviiround(q2neg, det12);
        if (tm1)
        {
          gcoeff(tm1, 1, i) = gadd(mulii(q1neg, gcoeff(tm,1,1)),
                                     mulii(q2neg, gcoeff(tm,1,2)));
          gcoeff(tm1, 2, i) = gadd(mulii(q1neg, gcoeff(tm,2,1)),
                                     mulii(q2neg, gcoeff(tm,2,2)));
        }
        gel(mid,i) = gadd(c, ZV_lincomb(q1neg,q2neg, gel(mid,1),gel(mid,2)));
      } /* for i */
    } /* local block */
  }
  if (DEBUGLEVEL>6)
  {
    if (tm1) fprintferr("tm1 = %Z",tm1);
    fprintferr("mid = %Z",mid); /* = m * tm1 */
  }
  gerepileall(av, tm1? 2: 1, &mid, &tm1);
  {
   /* For each pair of column vectors v and w in mid * tm2,
    * try to replace (v, w) by (v, v - q*w) for some q.
    * We compute all inner products and check them repeatedly. */
    const pari_sp av3 = avma, lim = stack_lim(av3,2);
    long i, j, npass = 0, e = VERYBIGINT;
    GEN dot = cgetg(ncol+1, t_MAT); /* scalar products */

    tm2 = matid(ncol);
    for (i=1; i <= ncol; i++)
    {
      gel(dot,i) = cgetg(ncol+1,t_COL);
      for (j=1; j <= i; j++)
	gcoeff(dot,j,i) = gcoeff(dot,i,j) = gscali(gel(mid,i),gel(mid,j));
    }
    for(;;)
    {
      long reductions = 0, olde = e;
      for (i=1; i <= ncol; i++)
      {
	long ijdif;
        for (ijdif=1; ijdif < ncol; ijdif++)
	{
          long d, k1, k2;
          GEN codi, q;

          j = i + ijdif; if (j > ncol) j -= ncol;
	  /* let k1, resp. k2,  index of larger, resp. smaller, column */
          if (cmpii(gcoeff(dot,i,i),
		    gcoeff(dot,j,j)) > 0) { k1 = i; k2 = j; }
          else                            { k1 = j; k2 = i; }
	  codi = gcoeff(dot,k2,k2);
          q = signe(codi)? diviiround(gcoeff(dot,k1,k2), codi): gen_0;
          if (!signe(q)) continue;

	  /* Try to subtract a multiple of column k2 from column k1.  */
          reductions++; q = negi(q);
          gel(tm2,k1) = ZV_lincomb(gen_1,q, gel(tm2,k1), gel(tm2,k2));
          gel(dot,k1) = ZV_lincomb(gen_1,q, gel(dot,k1), gel(dot,k2));
          gcoeff(dot, k1, k1) = addii(gcoeff(dot,k1,k1),
                                      mulii(q, gcoeff(dot,k2,k1)));
          for (d = 1; d <= ncol; d++) coeff(dot,k1,d) = coeff(dot,d,k1);
        } /* for ijdif */
        if (low_stack(lim, stack_lim(av3,2)))
	{
          if(DEBUGMEM>1) pari_warn(warnmem,"lllintpartialall");
	  gerepileall(av3, 2, &dot,&tm2);
        }
      } /* for i */
      if (!reductions) break;
      e = 0;
      for (i = 1; i <= ncol; i++) e += expi( gcoeff(dot,i,i) );
      if (e == olde) break;
      if (DEBUGLEVEL>6)
      {
        npass++;
	fprintferr("npass = %ld, red. last time = %ld, log_2(det) ~ %ld\n\n",
	            npass, reductions, e);
      }
    } /* for(;;)*/

   /* Sort columns so smallest comes first in m * tm1 * tm2.
    * Use insertion sort. */
    for (i = 1; i < ncol; i++)
    {
      long j, s = i;

      for (j = i+1; j <= ncol; j++)
	if (cmpii(gcoeff(dot,s,s),gcoeff(dot,j,j)) > 0) s = j;
      if (i != s)
      { /* Exchange with proper column */
        /* Only diagonal of dot is updated */
        lswap(tm2[i], tm2[s]);
        lswap(coeff(dot,i,i), coeff(dot,s,s));
      }
    }
    i = 1;
    while (i <= ncol && !signe(gcoeff(dot,i,i))) i++;
    if (i > 1)
    {
      tm2 += (i - 1);
      tm2[0] = evaltyp(t_MAT)|evallg(ncol - i);
    }
  } /* local block */
  return gerepileupto(av, gmul(tm1? tm1: mid, tm2));
}

GEN
lllintpartial(GEN mat)
{
  return lllintpartialall(mat,1);
}

GEN
lllintpartial_ip(GEN mat)
{
  return lllintpartialall(mat,0);
}

/********************************************************************/
/**                                                                **/
/**                    COPPERSMITH ALGORITHM                       **/
/**           Finding small roots of univariate equations.         **/
/**                                                                **/
/********************************************************************/

static int
check_condition(double beta, double tau, double rho, long d, long delta, long t)
{
  long dim = d*delta + t;
  double cond = d*delta*(delta+1)/2 - beta*delta*dim
    + rho*delta*(delta - 1) / 2
    + rho * t * delta + tau*dim*(dim - 1)/2;

  if (DEBUGLEVEL >= 4) 
    fprintferr("delta = %d, t = %d, cond = %lf\n", delta, t, cond); 

  return (cond <= 0);
}

static void
choose_params(GEN P, GEN N, GEN X, GEN B, long *pdelta, long *pt)
{
  long d = degpol(P);
  GEN P0 = leading_term(P);
  double logN = gtodouble(glog(N, DEFAULTPREC));
  double tau, beta, rho;
  long delta, t;
  tau = gtodouble(glog(X, DEFAULTPREC)) / logN;
  beta = gtodouble(glog(B, DEFAULTPREC)) / logN;
  if (tau >= beta * beta / d) pari_err(talker, "bound too large");
  /* TODO : remove P0 completely ! */
  rho = gtodouble(glog(P0, DEFAULTPREC)) / logN;

  /* Enumerate (delta,t) by increasing dimension of resulting lattice.
   * Not subtle, but o(1) for computing time */
  t = d; delta = 0;
  for(;;)
  {
    t += d * delta + 1; delta = 0;
    while (t >= 0) {
      if (check_condition(beta, tau, rho, d, delta, t)) {
        *pdelta = delta; *pt = t; return;
      }
      delta++; t -= d;
    }
  }
}

static GEN 
do_exhaustive(GEN P, GEN N, long x, GEN B) 
{
  GEN tst, sol = cget1(2*x + 2, t_VECSMALL);
  long j, l;

  for (j = -x; j <= x; j++)
  {
    tst = gcdii(FpX_eval(P, stoi(j), N), N);
    
    if (cmpii(tst, B) >= 0) /* We have found a factor of N >= B */
    {
      for (l = 1; l < lg(sol) && j != sol[l]; l++) /*empty*/;
      if (l == lg(sol)) appendL(sol, (GEN)j);
    }
  }
  return zv_to_ZV(sol); 
}

#define X_SMALL 1000

/* General Coppersmith, look for a root x0 <= p, p >= B, p | N, |x0| <= X */
GEN
zncoppersmith(GEN P0, GEN N, GEN X, GEN B)
{
  GEN Q, R, N0, M, sh, short_pol, *Xpowers, z, r, sol, nsp, P, tst, Z;
  long delta, i, j, row, d, l, dim, t, bnd = 10;
  pari_sp av = avma;

  if (typ(P0) != t_POL || typ(N) != t_INT) pari_err(typeer, "zncoppersmith");
  if (typ(X) != t_INT) {
    X = gfloor(X);
    if (typ(X) != t_INT) pari_err(typeer, "zncoppersmith");
  }
  if (signe(X) < 0) pari_err(talker, "negative bound in zncoppersmith");
  if (!B) B = N;
  if (typ(B) != t_INT) B = gceil(B);

  if (cmpis(X, X_SMALL) <= 0) 
    return gerepileupto(av, do_exhaustive(P0, N, itos(X), B)); 

  /* bnd-hack is only for the case B = N */
  if (!equalii(B,N)) bnd = 1;

  P = shallowcopy(P0); d = degpol(P);
  if (d == 0) { avma = av; return cgetg(1, t_VEC); }
  if (d < 0) pari_err(talker, "zero polynomial forbidden");

  if (!gcmp1(gel(P,d+2)))
  {
    gel(P,d+2) = bezout(gel(P,d+2), N, &z, &r);
    for (j = 0; j < d; j++) gel(P,j+2) = modii(mulii(gel(P,j+2), z), N);
  }
  if (DEBUGLEVEL >= 2) fprintferr("Modified P: %Z\n", P);

  choose_params(P, N, X, B, &delta, &t);
  if (DEBUGLEVEL >= 2)
    fprintferr("Init: trying delta = %d, t = %d\n", delta, t);
  for(;;)
  {
    dim = d * delta + t;

    /* TODO: In case of failure do not recompute the full vector */
    Xpowers = (GEN*)new_chunk(dim + 1);
    Xpowers[0] = gen_1;
    for (j = 1; j <= dim; j++) Xpowers[j] = gmul(Xpowers[j-1], X);

    /* TODO: in case of failure, use the part of the matrix already computed */
    M = cgetg(dim + 1, t_MAT);
    for (j = 1; j <= dim; j++) gel(M,j) = zerocol(dim);

    /* Rows of M correspond to the polynomials
     * N^delta, N^delta Xi, ... N^delta (Xi)^d-1,
     * N^(delta-1)P(Xi), N^(delta-1)XiP(Xi), ... N^(delta-1)P(Xi)(Xi)^d-1,
     * ...
     * P(Xi)^delta, XiP(Xi)^delta, ..., P(Xi)^delta(Xi)^t-1 */
    for (j = 1; j <= d;   j++) gcoeff(M, j, j) = gel(Xpowers,j-1);

    /* P-part */
    if (delta) row = d + 1; else row = 0; 

    Q = P;
    for (i = 1; i < delta; i++)
    {
      for (j = 0; j < d; j++,row++)
        for (l = j + 1; l <= row; l++)
          gcoeff(M, l, row) = mulii(Xpowers[l-1], gel(Q,l-j+1));
      Q = RgX_mul(Q, P);
    }
    for (j = 0; j < t; row++, j++)
      for (l = j + 1; l <= row; l++)
        gcoeff(M, l, row) = mulii(Xpowers[l-1], gel(Q,l-j+1));

    /* N-part */
    row = dim - t; N0 = N;
    while (row >= 1)
    {	
      for (j = 0; j < d; j++,row--)
        for (l = 1; l <= row; l++)
          gcoeff(M, l, row) = mulii(gmael(M, row, l), N0);
      if (row >= 1) N0 = mulii(N0, N);
    }
    /* Z is the upper bound for the L^1 norm of the polynomial, 
       ie. N^delta if B = N, B^delta otherwise */ 
    if (B != N) Z = powiu(B, delta); else Z = N0; 

    if (DEBUGLEVEL >= 2)
    {
      if (DEBUGLEVEL >= 6) fprintferr("Matrix to be reduced:\n%Z\n", M);
      fprintferr("Entering LLL\nbitsize bound: %ld\n", expi(Z));
      fprintferr("expected shvector bitsize: %ld\n", expi(dethnf_i(M))/dim);
    }

    sh = lllint_fp_ip(M, 4);
    /* Take the first vector if it is non constant */
    short_pol = gel(sh,1);
    for (j = 2; j <= dim; j++)
      if (signe(gel(short_pol, j))) break;
    if (j > dim) short_pol = gel(sh, 2);

    nsp = gen_0;
    for (j = 1; j <= dim; j++) nsp = addii(nsp, absi(gel(short_pol,j)));

    if (DEBUGLEVEL >= 2)
    {
      fprintferr("Candidate: %Z\n", short_pol);
      fprintferr("bitsize Norm: %ld\n", expi(nsp));
      fprintferr("bitsize bound: %ld\n", expi(mulsi(bnd, Z)));
    }
    if (cmpii(nsp, mulsi(bnd, Z)) < 0) break; /* SUCCESS */

    /* Failed with the precomputed or supplied value */
    t++; if (t == d) { delta++; t = 1; }
    if (DEBUGLEVEL >= 2)
      fprintferr("Increasing dim, delta = %d t = %d\n", delta, t);
  }
  bnd = itos(divii(nsp, Z)) + 1;

  while (!signe(short_pol[dim])) dim--;

  R = cgetg(dim + 2, t_POL); R[1] = P[1];
  for (j = 1; j <= dim; j++)
    gel(R,j+1) = diviiexact(gel(short_pol,j), Xpowers[j-1]);
  gel(R,2) = subii(gel(R,2), mulsi(bnd - 1, N0));

  sol = cgetg(1, t_VEC);
  for (i = -bnd + 1; i < bnd; i++)
  {
    r = nfrootsQ(R);
    if (DEBUGLEVEL >= 2) fprintferr("Roots: %Z\n", r);

    for (j = 1; j < lg(r); j++)
    {
      z = gel(r,j);
      tst = gcdii(FpX_eval(P, z, N), N);
    
      if (cmpii(tst, B) >= 0) /* We have found a factor of N >= B */
      {
        for (l = 1; l < lg(sol) && !equalii(z, gel(sol,l)); l++) /*empty*/;
        if (l == lg(sol)) sol = shallowconcat(sol, z);
      }
    }
    if (i < bnd) gel(R,2) = addii(gel(R,2), Z);
  }
  return gerepilecopy(av, sol);
}

/********************************************************************/
/**                                                                **/
/**                   LINEAR & ALGEBRAIC DEPENDENCE                **/
/**                                                                **/
/********************************************************************/

static int
real_indep(GEN re, GEN im, long bitprec)
{
  GEN p1 = gsub(gmul(gel(re,1),gel(im,2)),
                gmul(gel(re,2),gel(im,1)));
  return (!gcmp0(p1) && gexpo(p1) > - bitprec);
}

GEN
lindep2(GEN x, long bit)
{
  long tx=typ(x), lx=lg(x), ly, i, j, e;
  pari_sp av = avma;
  GEN re, im, M;

  if (! is_vec_t(tx)) pari_err(typeer,"lindep2");
  if (lx<=2) return cgetg(1,t_VEC);
  if (bit < 0) pari_err(talker, "negative accuracy in lindep2");
  if (!bit)
  {
    bit = gprecision(x);
    if (!bit)
    {
      x = primpart(x);
      bit = 32 + gexpo(x);
    }
    else
      bit = (long)bit_accuracy_mul(bit, 0.8);
  }
  else
    bit = (long) (bit/L2SL10);
  re = real_i(x);
  im = imag_i(x);
  /* independent over R ? */
  if (lx == 3 && real_indep(re,im,bit))
    { avma = av; return cgetg(1, t_VEC); }

  if (gcmp0(im)) im = NULL;
  ly = im? lx+2: lx+1;
  M = cgetg(lx,t_MAT);
  for (i=1; i<lx; i++)
  {
    GEN c = cgetg(ly,t_COL); gel(M,i) = c;
    for (j=1; j<lx; j++) gel(c,j) = (i==j)? gen_1: gen_0;
    gel(c,lx)           = gcvtoi(gshift(gel(re,i),bit), &e);
    if (im) gel(c,lx+1) = gcvtoi(gshift(gel(im,i),bit), &e);
  }
  M = lllint_fp_ip(M,100);
  M = gel(M,1);
  M[0] = evaltyp(t_COL) | evallg(lx);
  return gerepilecopy(av, M);
}

#define quazero(x) (gcmp0(x) || (typ(x)==t_REAL && expo(x) < EXP))
GEN
lindep(GEN x, long prec)
{
  GEN *b,*be,*bn,**m,qzer;
  GEN c1,c2,c3,px,py,pxy,re,im,p1,p2,r,f,em;
  long i,j,fl,k, lx = lg(x), tx = typ(x), n = lx-1;
  pari_sp av = avma, lim = stack_lim(av,1), av0, av1;
  const long EXP = - bit_accuracy(prec) + 2*n;

  if (! is_vec_t(tx)) pari_err(typeer,"lindep");
  if (n <= 1) return cgetg(1,t_VEC);
  x = gmul(x, real_1(prec)); if (tx != t_COL) settyp(x,t_COL);
  re = real_i(x);
  im = imag_i(x);
  /* independent over R ? */
  if (n == 2 && real_indep(re,im,bit_accuracy(prec)))
    { avma = av; return cgetg(1, t_VEC); }
  if (EXP > -10) pari_err(precer,"lindep");

  qzer = cgetg(lx, t_VECSMALL);
  b = (GEN*)matid(n);
  be= (GEN*)new_chunk(lx);
  bn= (GEN*)new_chunk(lx);
  m = (GEN**)new_chunk(lx);
  for (i=1; i<=n; i++)
  {
    bn[i] = cgetr(prec+1);
    be[i] = cgetg(lx,t_COL);
    m[i]  = (GEN*)new_chunk(lx);
    for (j=1; j<i ; j++) m[i][j] = cgetr(prec+1);
    for (j=1; j<=n; j++) be[i][j]= (long)cgetr(prec+1);
  }
  px = sqscal(re);
  py = sqscal(im); pxy = gscal(re,im);
  p1 = mpsub(mpmul(px,py), gsqr(pxy));
  if (quazero(px)) { re = im; px = py; fl = 1; } else fl = quazero(p1);
  av0 = av1 = avma;
  for (i=1; i<=n; i++)
  {
    p2 = gscal(b[i],re);
    if (fl) p2 = gmul(gdiv(p2,px),re);
    else
    {
      GEN p5,p6,p7;
      p5 = gscal(b[i],im);
      p6 = gdiv(gsub(gmul(py,p2),gmul(pxy,p5)), p1);
      p7 = gdiv(gsub(gmul(px,p5),gmul(pxy,p2)), p1);
      p2 = gadd(gmul(p6,re), gmul(p7,im));
    }
    p2 = gsub(b[i],p2);
    for (j=1; j<i; j++)
      if (qzer[j]) affrr(bn[j], m[i][j]);
      else
      {
        gdivz(gscal(b[i],be[j]),bn[j], m[i][j]);
        p2 = gsub(p2, gmul(m[i][j],be[j]));
      }
    gaffect(p2,          be[i]);
    affrr(sqscal(be[i]), bn[i]);
    qzer[i] = quazero(bn[i]); avma = av1;
  }
  while (qzer[n])
  {
    long e;
    if (DEBUGLEVEL>8)
    {
      fprintferr("qzer[%ld]=%ld\n",n,qzer[n]);
      for (k=1; k<=n; k++)
	for (i=1; i<k; i++) output(m[k][i]);
    }
    em=bn[1]; j=1;
    for (i=2; i<n; i++)
    {
      p1=shiftr(bn[i],i);
      if (cmprr(p1,em)>0) { em=p1; j=i; }
    }
    i=j; k=i+1;
    avma = av1; r = grndtoi(m[k][i], &e);
    if (e >= 0) pari_err(precer,"lindep");
    r  = negi(r);
    p1 = ZV_lincomb(gen_1,r, b[k],b[i]);
    b[k] = b[i];
    b[i]  = p1;
    av1 = avma;
    f = addir(r,m[k][i]);
    for (j=1; j<i; j++)
      if (!qzer[j])
      {
        p1 = mpadd(m[k][j], mulir(r,m[i][j]));
        affrr(m[i][j], m[k][j]);
        affrr(p1, m[i][j]);
      }
    c1 = addrr(bn[k], mulrr(gsqr(f),bn[i]));
    if (!quazero(c1))
    {
      c2 = divrr(mulrr(bn[i],f),c1); affrr(c2, m[k][i]);
      c3 = divrr(bn[k],c1);
      mulrrz(c3,bn[i], bn[k]); qzer[k] = quazero(bn[k]);
      affrr(c1,        bn[i]); qzer[i] = 0;
      for (j=i+2; j<=n; j++)
      {
        p1 = addrr(mulrr(m[j][k],c3), mulrr(m[j][i],c2));
        subrrz(m[j][i],mulrr(f,m[j][k]), m[j][k]);
        affrr(p1, m[j][i]);
      }
    }
    else
    {
      affrr(bn[i], bn[k]); qzer[k] = qzer[i];
      affrr(c1,    bn[i]); qzer[i] = 1;
      for (j=i+2; j<=n; j++) affrr(m[j][i], m[j][k]);
    }
    if (low_stack(lim, stack_lim(av,1)))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"lindep");
      b = (GEN*)gerepilecopy(av0, (GEN)b);
      av1 = avma;
    }
  }
  p1 = cgetg(lx,t_COL); gel(p1,n) = gen_1; for (i=1; i<n; i++) gel(p1,i) = gen_0;
  return gerepileupto(av, gauss(shallowtrans((GEN)b),p1));
}

/* PSLQ Programs */

typedef struct {
  long vmind, t12, t1234, reda, fin;
  long ct;
} pslq_timer;

/* WARNING: for double ** matrices, A[i] is the i-th ROW of A */
typedef struct {
  double *y, **H, **A, **B;
  double *W; /* scratch space */
  long n;
  pslq_timer *T;
} pslqL2_M;

typedef struct {
  GEN y, H, A, B;
  long n, EXP;
  int flreal;
  pslq_timer *T;
} pslq_M;

void
init_dalloc()
{ /* correct alignment for dalloc */
  (void)new_chunk((avma % sizeof(double)) / sizeof(long));
}

double *
dalloc(size_t n)
{
  return (double*)new_chunk(n / sizeof(long));
}

char *
stackmalloc(size_t N)
{
  long n = nchar2nlong(N);
  return (char*)new_chunk(n);
}

static double
conjd(double x) { return x; }

static double
sqrd(double a) { return a*a; }

static void
redall(pslq_M *M, long i, long jsup)
{
  long j, k, n = M->n;
  GEN t,b;
  GEN A = M->A, B = M->B, H = M->H, y = M->y;
  const GEN Bi = gel(B,i);

  for (j=jsup; j>=1; j--)
  {
    pari_sp av = avma;
    t = round_safe( gdiv(gcoeff(H,i,j), gcoeff(H,j,j)) );
    if (!t || gcmp0(t)) { avma = av; continue; }

    b = gel(B,j);
    gel(y,j) = gadd(gel(y,j), gmul(t,gel(y,i)));
    for (k=1; k<=j; k++)
      gcoeff(H,i,k) = gsub(gcoeff(H,i,k), gmul(t,gcoeff(H,j,k)));
    for (k=1; k<=n; k++)
    {
      gcoeff(A,i,k) = gsub(gcoeff(A,i,k), gmul(t,gcoeff(A,j,k)));
      gel(b,k) = gadd(gel(b,k), gmul(t,gel(Bi,k)));
    }
  }
}

static void
redallbar(pslqL2_M *Mbar, long i, long jsup)
{
  long j, k, n = Mbar->n;
  double t;
  double *hi = Mbar->H[i], *ai = Mbar->A[i], *hj, *aj;

#ifdef DEBUGPSLQ
fprintferr("%ld:\n==\n",i);
#endif
  for (j=jsup; j>=1; j--)
  {
    hj = Mbar->H[j];
    t = floor(0.5 + hi[j] / hj[j]);
    if (!t) continue;
#ifdef DEBUGPSLQ
fprintferr("%15.15e ",t);
#endif
    aj = Mbar->A[j];

    Mbar->y[j] += t * Mbar->y[i];
    for (k=1; k<=j; k++) hi[k] -= t * hj[k];
    for (k=1; k<=n; k++) {
      ai[k]         -= t * aj[k];
      Mbar->B[k][j] += t * Mbar->B[k][i];
    }
#ifdef DEBUGPSLQ
fprintferr("  %ld:\n",j); dprintmat(Mbar->H,n,n-1);
#endif
  }
}

static long
vecabsminind(GEN v)
{
  long l = lg(v), m = 1, i;
  GEN t, la = mpabs(gel(v,1));
  for (i=2; i<l; i++)
  {
    t = mpabs(gel(v,i));
    if (mpcmp(t,la) < 0) { la = t; m = i; }
  }
  return m;
}

static long
vecmaxind(GEN v)
{
  long l = lg(v), m = 1, i;
  GEN t, la = gel(v,1);
  for (i=2; i<l; i++)
  {
    t = gel(v,i);
    if (mpcmp(t,la) > 0) { la = t; m = i; }
  }
  return m;
}

static long
vecmaxindbar(double *v, long n)
{
  long m = 1, i;
  double la = v[1];
  for (i=2; i<=n; i++)
    if (v[i] > la) { la = v[i]; m = i; }
  return m;
}

static GEN
maxnorml2(pslq_M *M)
{
  long n = M->n, i, j;
  GEN ma = gen_0, s;

  for (i=1; i<=n; i++)
  {
    s = gen_0;
    for (j=1; j<n; j++) s = gadd(s, gnorm(gcoeff(M->H,i,j)));
    ma = gmax(ma, s);
  }
  return sqrtr(gmul(ma, real_1(DEFAULTPREC)));
}

static void
init_timer(pslq_timer *T)
{
  T->vmind = T->t12 = T->t1234 = T->reda = T->fin = T->ct = 0;
}

static int
is_zero(GEN x, long e, long prec)
{
  if (gcmp0(x)) return 1;
  if (typ(x) == t_REAL)
  {
    long ex = expo(x);
    return (ex < e || (prec != 3 && lg(x) == 3 && ex < (e>>1)));
  }
  return gexpo(x) < e;
}

static GEN
init_pslq(pslq_M *M, GEN x, long *PREC)
{
  long tx = typ(x), lx = lg(x), n = lx-1, i, j, k, prec;
  GEN s1, s, sinv;

  if (! is_vec_t(tx)) pari_err(typeer,"pslq");
  /* check trivial cases */
  for (k = 1; k <= n; k++)
    if (gcmp0(gel(x,k))) return col_ei(n, k);
  if (n <= 1) return cgetg(1, t_COL);
  prec = gprecision(x) - 1; /* don't trust the last word */
  if (prec < 0)
  { /* exact components */
    pari_sp av = avma;
    GEN im, U = NULL;
    x = Q_primpart(x);
    im = imag_i(x);
    x = real_i(x); settyp(x, t_VEC);
    if (!gcmp0(im))
    {
      U = (GEN)extendedgcd(im)[2];
      setlg(U, lg(U)-1); /* remove last column */
      x = gmul(x, U);
      if (n == 2) /* x has a single component */
        return gcmp0(gel(x,1))? gel(U,1): cgetg(1, t_COL);
    }
    x = (GEN)extendedgcd(x)[2];
    x = gel(x,1);
    if (U) x = gmul(U, x);
    return gerepilecopy(av, x);
  }
  if (prec < DEFAULTPREC) prec = DEFAULTPREC;
  *PREC = prec;
  M->EXP = - bit_accuracy(prec) + max(n, 8);
  M->flreal = is_zero(imag_i(x), M->EXP, prec);
  if (!M->flreal)
    return lindep(x,prec); /* FIXME */
  else
    x = real_i(x);

  if (DEBUGLEVEL>=3) { (void)timer(); init_timer(M->T); }
  x = col_to_MP(x, prec); settyp(x,t_VEC);
  M->n = n;
  M->A = matid(n);
  M->B = matid(n);
  s1 = cgetg(lx,t_VEC); gel(s1,n) = gnorm(gel(x,n));
  s  = cgetg(lx,t_VEC); gel(s,n) = gabs(gel(x,n),prec);
  for (k=n-1; k>=1; k--)
  {
    gel(s1,k) = gadd(gel(s1,k+1), gnorm(gel(x,k)));
    gel(s,k) = gsqrt(gel(s1,k), prec);
  }
  sinv = ginv(gel(s,1));
  s    = gmul(sinv,s);
  M->y = gmul(sinv, x);
  M->H = cgetg(n,t_MAT);
  for (j=1; j<n; j++)
  {
    GEN d, c = cgetg(lx,t_COL);

    gel(M->H,j) = c;
    for (i=1; i<j; i++) gel(c,i) = gen_0;
    gel(c,j) = gdiv(gel(s,j+1),gel(s,j));
    d = gneg( gdiv(gel(M->y,j), gmul(gel(s,j),gel(s,j+1)) ));
    for (i=j+1; i<=n; i++) gel(c,i) = gmul(gconj(gel(M->y,i)), d);
  }
  for (i=2; i<=n; i++) redall(M, i, i-1);
  return NULL;
}

static void
SWAP(pslq_M *M, long m)
{
  long j, n = M->n;
  lswap(M->y[m], M->y[m+1]);
  lswap(M->B[m], M->B[m+1]);
  for (j=1; j<=n; j++) lswap(coeff(M->A,m,j), coeff(M->A,m+1,j));
  for (j=1; j<n;  j++) lswap(coeff(M->H,m,j), coeff(M->H,m+1,j));
}

static void
SWAPbar(pslqL2_M *M, long m)
{
  long j, n = M->n;
  dswap(M->y[m], M->y[m+1]);
  pdswap(M->A[m], M->A[m+1]);
  pdswap(M->H[m], M->H[m+1]);
  for (j=1; j<=n; j++) dswap(M->B[j][m], M->B[j][m+1]);
}

static GEN
one_step_gen(pslq_M *M, GEN tabga, long prec)
{
  GEN H = M->H, p1;
  long n = M->n, i, m;

  p1 = cgetg(n,t_VEC);
  for (i=1; i<n; i++) gel(p1,i) = gmul(gel(tabga,i), gabs(gcoeff(H,i,i),prec));
  m = vecmaxind(p1);
  if (DEBUGLEVEL>3) M->T->vmind += timer();
  SWAP(M, m);
  if (m <= n-2)
  {
    GEN tinv, t3, t4, t1c, t2c, t1 = gcoeff(H,m,m), t2 = gcoeff(H,m,m+1);
    tinv = ginv( gsqrt(gadd(gnorm(t1), gnorm(t2)), prec) );
    t1 = gmul(t1, tinv);
    t2 = gmul(t2, tinv);
    if (M->flreal) { t1c = t1; t2c = t2; }
    else           { t1c = gconj(t1); t2c = gconj(t2); }
    if (DEBUGLEVEL>3) M->T->t12 += timer();
    for (i=m; i<=n; i++)
    {
      t3 = gcoeff(H,i,m);
      t4 = gcoeff(H,i,m+1);
      gcoeff(H,i,m) = gadd(gmul(t1c,t3), gmul(t2c,t4));
      gcoeff(H,i,m+1) = gsub(gmul(t1, t4), gmul(t2, t3));
    }
    if (DEBUGLEVEL>3) M->T->t1234 += timer();
  }
  for (i=1; i<n; i++)
    if (is_zero(gcoeff(H,i,i), M->EXP, prec)) {
      m = vecabsminind(M->y); return (GEN)M->B[m];
    }
  for (i=m+1; i<=n; i++) redall(M, i, min(i-1,m+1));

  if (DEBUGLEVEL>3) M->T->reda += timer();
  if (gexpo(M->A) >= -M->EXP) return ginv(maxnorml2(M));
  m = vecabsminind(M->y);
  if (is_zero((GEN)M->y[m], M->EXP, prec)
   && gexpo(M->y) - gexpo((GEN)M->y[m]) > 20)
    return (GEN)M->B[m];

  if (DEBUGLEVEL>2)
  {
    if (DEBUGLEVEL>3) M->T->fin += timer();
    M->T->ct++;
    if ((M->T->ct&0xff) == 0)
    {
      if (DEBUGLEVEL == 3)
        fprintferr("time for ct = %ld : %ld\n",M->T->ct,timer());
      else
        fprintferr("time [max,t12,loop,reds,fin] = [%ld, %ld, %ld, %ld, %ld]\n",
                   M->T->vmind, M->T->t12, M->T->t1234, M->T->reda, M->T->fin);
    }
  }
  return NULL; /* nothing interesting */
}

static GEN
get_tabga(int flreal, long n, long prec)
{
  GEN ga = sqrtr( flreal? divrs(stor(4, prec), 3): stor(2, prec) );
  GEN tabga = cgetg(n,t_VEC);
  long i;
  gel(tabga,1) = ga;
  for (i = 2; i < n; i++) gel(tabga,i) = gmul(gel(tabga,i-1),ga);
  return tabga;
}

GEN
pslq(GEN x)
{
  GEN tabga, p1;
  long prec;
  pari_sp av0 = avma, lim = stack_lim(av0,1), av;
  pslq_M M;
  pslq_timer T; M.T = &T;

  p1 = init_pslq(&M, x, &prec);
  if (p1) return p1;

  tabga = get_tabga(M.flreal, M.n, prec);
  av = avma;
  if (DEBUGLEVEL>=3) printf("Initialization time = %ld\n",timer());
  for (;;)
  {
    if ((p1 = one_step_gen(&M, tabga, prec)))
      return gerepilecopy(av0, p1);

    if (low_stack(lim, stack_lim(av,1)))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"pslq");
      gerepileall(av,4,&M.y,&M.H,&M.A,&M.B);
    }
  }
}

/* W de longueur n-1 */

static double
dnorml2(double *W, long n, long row)
{
  long i;
  double s = 0.;

  for (i=row; i<n; i++) s += W[i]*W[i];
  return s;
}

/* Hbar *= Pbar */
static void
dmatmul(pslqL2_M *Mbar, double **Pbar, long row)
{
  const long n = Mbar->n; /* > row */
  long i, j, k;
  double s, *W = Mbar->W, **H = Mbar->H;

  for (i = row; i <= n; i++)
  {
    for (j = row; j < n; j++)
    {
      k = row; s = H[i][k] * Pbar[k][j];
      for ( ; k < n; k++) s += H[i][k] * Pbar[k][j];
      W[j] = s;
    }
    for (j = row; j < n; j++) H[i][j] = W[j];
  }
}

/* compute n-1 x n-1 matrix Pbar */
static void
dmakep(pslqL2_M *Mbar, double **Pbar, long row)
{
  long i, j, n = Mbar->n;
  double pro, nc, *C = Mbar->H[row], *W = Mbar->W;

  nc = sqrt(dnorml2(C,n,row));
  W[row] = (C[row] < 0) ? C[row] - nc : C[row] + nc;
  for (i=row; i<n; i++) W[i] = C[i];
  pro = -2.0 / dnorml2(W, n, row);
      /* must have dnorml2(W,n,row) = 2*nc*(nc+fabs(C[1])) */
  for (j=row; j<n; j++)
  {
    for (i=j+1; i<n; i++)
      Pbar[j][i] = Pbar[i][j] = pro * W[i] * W[j];
    Pbar[j][j] = 1.0 + pro * W[j] * W[j];
  }
}

static void
dLQdec(pslqL2_M *Mbar, double **Pbar)
{
  long row, j, n = Mbar->n;

  for (row=1; row<n; row++)
  {
    dmakep(Mbar, Pbar, row);
    dmatmul(Mbar, Pbar, row);
    for (j=row+1; j<n; j++) Mbar->H[row][j] = 0.;
  }
}

#ifdef DEBUGPSLQ
static void
dprintvec(double *V, long m)
{
  long i;
  fprintferr("[");
  for (i=1; i<m; i++) fprintferr("%15.15e, ",V[i]);
  fprintferr("%15.15e]\n",V[m]); pariflush();
}

static void
dprintmat(double **M, long r, long c)
{
  long i, j;
  fprintferr("[");
  for (i=1; i<r; i++)
  {
    for (j=1; j<c; j++) fprintferr("%15.15e, ",M[i][j]);
    fprintferr("%15.15e;\n ",M[i][c]);
  }
  for (j=1; j<c; j++) fprintferr("%15.15e, ",M[r][j]);
  fprintferr("%15.15e]\n",M[r][c]); pariflush();
}
#endif

static long
initializedoubles(pslqL2_M *Mbar, pslq_M *M, long prec)
{
  long i, j, n = Mbar->n;
  GEN ypro;
  pari_sp av = avma;

  ypro = gdiv(M->y, vecmax(gabs(M->y,prec)));
  for (i=1; i<=n; i++)
  {
    if (gexpo(gel(ypro,i)) < -0x3ff) return 0;
    Mbar->y[i] = rtodbl(gel(ypro,i));
  }
  avma = av;
  for (j=1; j<=n; j++)
    for (i=1; i<=n; i++)
    {
      if (i==j) Mbar->A[i][j] = Mbar->B[i][j] = 1.;
      else      Mbar->A[i][j] = Mbar->B[i][j] = 0.;
      if (j < n)
      {
        GEN h = gcoeff(M->H,i,j);
        if (!gcmp0(h) && labs(gexpo(h)) > 0x3ff) return 0;
        Mbar->H[i][j] = rtodbl(h);
      }
    }
  return 1;
}

/* T(arget) := S(ource) */
static void
storeprecdoubles(pslqL2_M *T, pslqL2_M *S)
{
  long n = T->n, i, j;

  for (i=1; i<=n; i++)
  {
    for (j=1; j<n; j++)
    {
      T->H[i][j] = S->H[i][j];
      T->A[i][j] = S->A[i][j];
      T->B[i][j] = S->B[i][j];
    }
    T->A[i][n] = S->A[i][n];
    T->B[i][n] = S->B[i][n];
    T->y[i] = S->y[i];
  }
}

static long
checkentries(pslqL2_M *Mbar)
{
  long n = Mbar->n, i, j;
  double *p1, *p2;

  for (i=1; i<=n; i++)
  {
    if (dblexpo(Mbar->y[i]) < -46) return 0;
    p1 = Mbar->A[i];
    p2 = Mbar->B[i];
    for (j=1; j<=n; j++)
      if (dblexpo(p1[j]) > 43 || dblexpo(p2[j]) > 43) return 0;
  }
  return 1;
}

static long
applybar(pslq_M *M, pslqL2_M *Mbar, GEN Abargen, GEN Bbargen)
{
  long n = Mbar->n, i, j;
  double *p1, *p2;

  for (i=1; i<=n; i++)
  {
    p1 = Mbar->A[i];
    p2 = Mbar->B[i];
    for (j=1; j<=n; j++)
    {
      if (dblexpo(p1[j]) >= 52 || dblexpo(p2[j]) >= 52) return 0;
      gcoeff(Abargen,i,j) = stoi((long)floor(p1[j]));
      gcoeff(Bbargen,i,j) = stoi((long)floor(p2[j]));
    }
  }
  M->y = gmul(M->y, Bbargen);
  M->B = gmul(M->B, Bbargen);
  M->A = gmul(Abargen, M->A);
  M->H = gmul(Abargen, M->H); return 1;
}

static GEN
checkend(pslq_M *M, long prec)
{
  long i, m, n = M->n;

  for (i=1; i<=n-1; i++)
    if (is_zero(gcoeff(M->H,i,i), M->EXP, prec))
    {
      m = vecabsminind(M->y);
      return (GEN)M->B[m];
    }
  if (gexpo(M->A) >= -M->EXP)
    return ginv( maxnorml2(M) );
  m = vecabsminind(M->y);
  if (is_zero((GEN)M->y[m], M->EXP, prec)) return (GEN)M->B[m];
  return NULL;
}

GEN
pslqL2(GEN x)
{
  GEN Abargen, Bbargen, tabga, p1;
  long lx = lg(x), n = lx-1, i, m, ctpro, flreal, flit, prec;
  pari_sp av0 = avma, lim = stack_lim(av0,1), av;
  double *tabgabar, gabar, tinvbar, t1bar, t2bar, t3bar, t4bar;
  double **Pbar, **Abar, **Bbar, **Hbar;
  pslqL2_M Mbar, Mbarst;
  pslq_M M;
  pslq_timer T; M.T = &T;

  p1 = init_pslq(&M, x, &prec);
  if (p1) return p1;

  flreal = M.flreal;
  tabga = get_tabga(flreal, n, prec);
  Abargen = matid(n);
  Bbargen = matid(n);

  Mbarst.n = Mbar.n = n;
  Mbar.A = Abar = (double**)new_chunk(n+1);
  Mbar.B = Bbar = (double**)new_chunk(n+1);
  Mbar.H = Hbar = (double**)new_chunk(n+1);
  Mbarst.A = (double**)new_chunk(n+1);
  Mbarst.B = (double**)new_chunk(n+1);
  Mbarst.H = (double**)new_chunk(n+1);
  Pbar   = (double**)new_chunk(n);

  tabgabar = dalloc((n+1)*sizeof(double));
  Mbar.y = dalloc((n+1)*sizeof(double));
  Mbarst.y = dalloc((n+1)*sizeof(double));

  Mbar.W = dalloc((n+1)*sizeof(double));
  for (i=1; i< n; i++)  Pbar[i] = dalloc((n+1)*sizeof(double));
  for (i=1; i<=n; i++)  Abar[i] = dalloc((n+1)*sizeof(double));
  for (i=1; i<=n; i++)  Bbar[i] = dalloc((n+1)*sizeof(double));
  for (i=1; i<=n; i++)  Hbar[i] = dalloc(n*sizeof(double));
  for (i=1; i<=n; i++) Mbarst.A[i] = dalloc((n+1)*sizeof(double));
  for (i=1; i<=n; i++) Mbarst.B[i] = dalloc((n+1)*sizeof(double));
  for (i=1; i<=n; i++) Mbarst.H[i] = dalloc(n*sizeof(double));

  gabar = gtodouble(gel(tabga,1)); tabgabar[1] = gabar;
  for (i=2; i<n; i++) tabgabar[i] = tabgabar[i-1]*gabar;

  av = avma;
  if (DEBUGLEVEL>=3) printf("Initialization time = %ld\n",timer());
RESTART:
  flit = initializedoubles(&Mbar, &M, prec);
  storeprecdoubles(&Mbarst, &Mbar);
  if (flit) dLQdec(&Mbar, Pbar);
  ctpro = 0;
  for (;;)
  {
    if (low_stack(lim, stack_lim(av,1)))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"pslq");
      gerepileall(av,4,&M.y,&M.H,&M.A,&M.B);
    }
    if (flit)
    {
      ctpro++;
      for (i=1; i<n; i++) Mbar.W[i] = tabgabar[i]*fabs(Hbar[i][i]);
      m = vecmaxindbar(Mbar.W, n-1);
      SWAPbar(&Mbar, m);
      if (m <= n-2)
      {
	tinvbar = 1.0 / sqrt(sqrd(Hbar[m][m]) + sqrd(Hbar[m][m+1]));
	t1bar = tinvbar*Hbar[m][m];
        t2bar = tinvbar*Hbar[m][m+1];
	if (DEBUGLEVEL>=4) T.t12 += timer();
	for (i=m; i<=n; i++)
	{
	  t3bar = Hbar[i][m];
          t4bar = Hbar[i][m+1];
	  if (flreal)
            Hbar[i][m] = t1bar*t3bar + t2bar*t4bar;
	  else
            Hbar[i][m] = conjd(t1bar)*t3bar + conjd(t2bar)*t4bar;
	  Hbar[i][m+1] = t1bar*t4bar - t2bar*t3bar;
	}
	if (DEBUGLEVEL>=4) T.t1234 += timer();
      }

      flit = checkentries(&Mbar);
      if (flit)
      {
	storeprecdoubles(&Mbarst, &Mbar);
	for (i=m+1; i<=n; i++) redallbar(&Mbar, i, min(i-1,m+1));
      }
      else
      {
	if (applybar(&M, &Mbar, Abargen,Bbargen))
	{
	  if ( (p1 = checkend(&M, prec)) ) return gerepilecopy(av0, p1);
	  goto RESTART;
	}
	else
        {
          if (ctpro == 1) goto DOGEN;
          storeprecdoubles(&Mbar, &Mbarst); /* restore */
          if (! applybar(&M, &Mbar, Abargen,Bbargen)) pari_err(bugparier,"pslqL2");
	  if ( (p1 = checkend(&M, prec)) ) return gerepilecopy(av0, p1);
          goto RESTART;
        }
      }
    }
    else
    {
DOGEN:
      if ((p1 = one_step_gen(&M, tabga, prec)))
        return gerepilecopy(av, p1);
    }
  }
}

/* x is a vector of elts of a p-adic field */
GEN
plindep(GEN x)
{
  long i, j, prec = VERYBIGINT, nx = lg(x)-1, v;
  pari_sp av = avma;
  GEN p = NULL, pn,p1,m,a;

  if (nx < 2) return cgetg(1,t_VEC);
  for (i=1; i<=nx; i++)
  {
    p1 = gel(x,i);
    if (typ(p1) != t_PADIC) continue;

    j = precp(p1); if (j < prec) prec = j;
    if (!p) p = gel(p1,2);
    else if (!equalii(p, gel(p1,2)))
      pari_err(talker,"inconsistent primes in plindep");
  }
  if (!p) pari_err(talker,"not a p-adic vector in plindep");
  v = ggval(x,p); pn = powiu(p,prec);
  if (v != 0) x = gmul(x, gpowgs(p, -v));
  x = RgV_to_FpV(x, pn);

  a = negi(gel(x,1));
  m = cgetg(nx,t_MAT);
  for (i=1; i<nx; i++)
  {
    GEN c = zerocol(nx);
    gel(c,1+i) = a;
    gel(c,1) = gel(x,i+1);
    gel(m,i) = c;
  }
  m = lllintpartial_ip( hnfmodid(m, pn) );
  m = lllint_fp_ip(m, 100);
  return gerepilecopy(av, gel(m,1));
}

GEN
lindep0(GEN x,long bit,long prec)
{
  long i, tx = typ(x);
  if (! is_vec_t(tx) && tx != t_MAT) pari_err(typeer,"lindep");
  for (i = 1; i < lg(x); i++)
    if (typ(gel(x,i)) == t_PADIC) return plindep(x);
  switch (bit)
  {
    case -1: return lindep(x,prec);
    case -2: return deplin(x);
    case -3: return pslq(x);
    case -4: return pslqL2(x);
    default: return lindep2(x, bit);
  }
}

GEN
algdep0(GEN x, long n, long bit, long prec)
{
  long tx = typ(x), i;
  pari_sp av;
  GEN y;

  if (! is_scalar_t(tx)) pari_err(typeer,"algdep0");
  if (tx==t_POLMOD) { y = gcopy(gel(x,1)); setvarn(y,0); return y; }
  if (gcmp0(x)) return pol_x[0];
  if (n <= 0)
  {
    if (!n) return gen_1;
    pari_err(talker,"negative polynomial degree in algdep");
  }

  av = avma; y = cgetg(n+2,t_COL);
  gel(y,1) = gen_1;
  gel(y,2) = x; /* n >= 1 */
  for (i=3; i<=n+1; i++) gel(y,i) = gmul(gel(y,i-1),x);
  if (typ(x) == t_PADIC)
    y = plindep(y);
  else
  {
    y = lindep0(y, bit, prec);
    if (typ(y) == t_REAL) return gerepileupto(av, y);
  }
  if (lg(y) < 2) pari_err(talker,"higher degree than expected in algdep");
  y = RgV_to_RgX(y, 0);
  if (gsigne(leading_term(y)) > 0) return gerepilecopy(av, y);
  return gerepileupto(av, gneg(y));
}

GEN
algdep2(GEN x, long n, long bit)
{
  return algdep0(x,n,bit,0);
}

GEN
algdep(GEN x, long n, long prec)
{
  return algdep0(x,n,0,prec);
}

/********************************************************************/
/**                                                                **/
/**                   INTEGRAL KERNEL (LLL REDUCED)                **/
/**                                                                **/
/********************************************************************/

GEN
matkerint0(GEN x, long flag)
{
  switch(flag)
  {
    case 0: return kerint(x);
    case 1: return kerint1(x);
    default: pari_err(flagerr,"matkerint");
  }
  return NULL; /* not reached */
}

GEN
kerint1(GEN x)
{
  pari_sp av = avma;
  return gerepilecopy(av, lllint_fp_ip(matrixqz3(ker(x)), 100));
}

GEN
kerint(GEN x)
{
  pari_sp av = avma;
  GEN fl, junk, h = lllint_i(x, 0, 0, &junk, &fl, NULL);
  if (h) h = lll_finish(h,fl, lll_KER);
  else   h = lll_trivial(x, lll_KER);
  if (lg(h)==1) { avma = av; return cgetg(1, t_MAT); }
  return gerepilecopy(av, lllint_ip(h, 100));
}

/********************************************************************/
/**                                                                **/
/**                              MINIM                             **/
/**                                                                **/
/********************************************************************/
/* x non-empty ZM, y a compatible zc (dimension > 0). */
static GEN
ZM_zc_mul_i(GEN x, GEN y, long c, long l)
{
  long i, j;
  pari_sp av;
  GEN z = cgetg(l,t_COL), s;

  for (i=1; i<l; i++)
  {
    av = avma; s = mulis(gcoeff(x,i,1),y[1]);
    for (j=2; j<c; j++)
      if (y[j]) s = addii(s, mulis(gcoeff(x,i,j),y[j]));
    gel(z,i) = gerepileuptoint(av,s);
  }
  return z;
}
GEN
ZM_zc_mul(GEN x, GEN y) {
  long l = lg(x);
  if (l == 1) return cgetg(1, t_COL);
  return ZM_zc_mul_i(x,y, l, lg(x[1]));
}

/* x ZM, y a compatible zm (dimension > 0). */
GEN
ZM_zm_mul(GEN x, GEN y)
{
  long j, c, l = lg(x), ly = lg(y);
  GEN z = cgetg(ly, t_MAT);
  if (l == 1) return z;
  c = lg(x[1]);
  for (j = 1; j < ly; j++) gel(z,j) = ZM_zc_mul_i(x, gel(y,j), l,c);
  return z;
}

void
minim_alloc(long n, double ***q, GEN *x, double **y,  double **z, double **v)
{
  long i, s;

  *x = cgetg(n, t_VECSMALL);
  *q = (double**) new_chunk(n);
  s = n * sizeof(double);
  init_dalloc();
  *y = dalloc(s);
  *z = dalloc(s);
  *v = dalloc(s);
  for (i=1; i<n; i++) (*q)[i] = dalloc(s);
}

/* If V depends linearly from the columns of the matrix, return 0.
 * Otherwise, update INVP and L and return 1. No GC. */
static int
addcolumntomatrix(GEN V, GEN invp, GEN L)
{
  GEN a = RgM_zc_mul(invp,V);
  long i,j,k, n = lg(invp);

  if (DEBUGLEVEL>6)
  {
    fprintferr("adding vector = %Z\n",V);
    fprintferr("vector in new basis = %Z\n",a);
    fprintferr("list = %Z\n",L);
    fprintferr("base change matrix =\n"); outerr(invp);
  }
  k = 1; while (k<n && (L[k] || gcmp0(gel(a,k)))) k++;
  if (k == n) return 0;
  L[k] = 1;
  for (i=k+1; i<n; i++) gel(a,i) = gdiv(gneg_i(gel(a,i)),gel(a,k));
  for (j=1; j<=k; j++)
  {
    GEN c = gel(invp,j), ck = gel(c,k);
    if (gcmp0(ck)) continue;
    gel(c,k) = gdiv(ck, gel(a,k));
    if (j==k)
      for (i=k+1; i<n; i++)
	gel(c,i) = gmul(gel(a,i), ck);
    else
      for (i=k+1; i<n; i++)
	gel(c,i) = gadd(gel(c,i), gmul(gel(a,i), ck));
  }
  return 1;
}

/* Minimal vectors for the integral definite quadratic form: a.
 * Result u:
 *   u[1]= Number of vectors of square norm <= BORNE
 *   u[2]= maximum norm found
 *   u[3]= list of vectors found (at most STOCKMAX)
 *
 *  If BORNE = gen_0: Minimal non-zero vectors.
 *  flag = min_ALL,   as above
 *  flag = min_FIRST, exits when first suitable vector is found.
 *  flag = min_PERF,  only compute rank of the family of v.v~ (v min.)
 *  flag = min_VECSMALL, return a t_VECSMALL of (half) the number of vectors for each norm
 *  flag = min_VECSMALL2, same but count only vectors with even norm, and shift the answer
 */
static GEN
minim0(GEN a, GEN BORNE, GEN STOCKMAX, long flag)
{
  GEN x,res,p1,u,r,L,gnorme,invp,V;
  long n = lg(a), i, j, k, s, maxrank;
  pari_sp av0 = avma, av1, av, lim;
  double p,maxnorm,BOUND,*v,*y,*z,**q;
  const double eps = 0.0001;

  if (!BORNE) BORNE = gen_0;
  if (!STOCKMAX) pari_err(talker,"maximal number of vectors must be provided");
  BORNE = gfloor(BORNE);
  if (typ(BORNE) != t_INT || typ(STOCKMAX) != t_INT)
    pari_err(typeer, "minim0");
  if (typ(a) != t_MAT) pari_err(typeer,"minim0");

  maxrank = 0; res = V = invp = NULL; /* gcc -Wall */
  switch(flag)
  {
    case min_FIRST:
      if (gcmp0(BORNE)) pari_err(talker,"bound = 0 in minim2");
      res = cgetg(3,t_VEC); break;
    case min_ALL: res = cgetg(4,t_VEC); break;
    case min_PERF: break;
    case min_VECSMALL:
    case min_VECSMALL2:
      maxrank = itos(BORNE);
      if (maxrank <= 0) return cgetg(1, t_VECSMALL);

      res = const_vecsmall(maxrank, 0);
      if (flag == min_VECSMALL2) BORNE = shifti(BORNE,1);
      if (gcmp0(BORNE)) return res;
      break;
    default: pari_err(flagerr, "minim0");
  }
  if (n == 1)
  {
    switch(flag)
    {
      case min_FIRST: avma=av0; return cgetg(1,t_VEC);
      case min_VECSMALL:
      case min_VECSMALL2: return res;
      case min_PERF:  avma=av0; return gen_0;
    }
    gel(res,1) = gel(res,2) = gen_0;
    gel(res,3) = cgetg(1,t_MAT); return res;
  }

  av = avma;
  minim_alloc(n, &q, &x, &y, &z, &v);
  av1 = avma;

  u = lllgramint(a);
  if (lg(u) != n) pari_err(talker,"not a definite form in minim0");
  a = qf_base_change(a,u,1);

  n--;
  a = mat_to_MP(a, DEFAULTPREC); r = sqred1(a);
  for (j=1; j<=n; j++)
  {
    v[j] = rtodbl(gcoeff(r,j,j));
    for (i=1; i<j; i++) q[i][j] = rtodbl(gcoeff(r,i,j));
  }

  if (flag==min_PERF || gcmp0(BORNE))
  {
    double c;
    p = rtodbl(gcoeff(a,1,1));
    for (i=2; i<=n; i++) { c = rtodbl(gcoeff(a,i,i)); if (c < p) p = c; }
    BORNE = roundr(dbltor(p));
    maxnorm = -1.; /* don't update maxnorm */
  }
  else
  {
    p = gtodouble(BORNE);
    maxnorm = 0.;
  }
  BOUND = p + eps;
  if (BOUND == p) pari_err(precer, "minim0");

  switch(flag)
  {
    case min_ALL:
      maxrank = itos(STOCKMAX);
      if (maxrank < 0) pari_err(talker,"negative number of vectors in minim0");
      L = new_chunk(1+maxrank);
      break;
    case min_PERF:
      BORNE = gerepileupto(av1,BORNE);
      maxrank = (n*(n+1)) >> 1;
      L = const_vecsmall(maxrank, 0);
      V = cgetg(1+maxrank, t_VECSMALL);
  }

  s = 0; av1 = avma; lim = stack_lim(av1,1);
  k = n; y[n] = z[n] = 0;
  x[n] = (long)sqrt(BOUND/v[n]);
  if (flag == min_PERF) invp = matid(maxrank);
  for(;;x[1]--)
  {
    do
    {
      if (k>1)
      {
        long l = k-1;
	z[l] = 0;
	for (j=k; j<=n; j++) z[l] += q[l][j]*x[j];
	p = (double)x[k] + z[k];
	y[l] = y[k] + p*p*v[k];
	x[l] = (long)floor(sqrt((BOUND-y[l])/v[l])-z[l]);
        k = l;
      }
      for(;;)
      {
	p = (double)x[k] + z[k];
	if (y[k] + p*p*v[k] <= BOUND) break;
	k++; x[k]--;
      }
    }
    while (k > 1);
    if (! x[1] && y[1]<=eps) break;
    p = (double)x[1] + z[1]; p = y[1] + p*p*v[1]; /* norm(x) */
    if (maxnorm >= 0)
    {
      if (flag == min_FIRST)
      {
        gel(res,2) = gerepileupto(av, ZM_zc_mul(u,x));
        av = avma;
        gel(res,1) = gerepileupto(av, ground(dbltor(p))); return res;
      }
      if (p > maxnorm) maxnorm = p;
    }
    else
    {
      pari_sp av2 = avma;
      gnorme = ground(dbltor(p));
      if (gcmp(gnorme,BORNE) >= 0) avma = av2;
      else
      {
        BOUND=gtodouble(gnorme)+eps; s=0;
        affii(gnorme,BORNE); avma = av1;
        if (flag == min_PERF) invp = matid(maxrank);
      }
    }
    s++;

    switch(flag)
    {
      case min_ALL:
        if (s<=maxrank)
        {
          p1 = new_chunk(n+1); gel(L,s) = p1;
          for (i=1; i<=n; i++) p1[i] = x[i];
        }
        break;

      case min_VECSMALL:
	{
	  ulong norm = (ulong)(p + 0.5);
	  res[norm]++;
	}
	break;

      case min_VECSMALL2:
	{
	  ulong norm = (ulong)(p + 0.5);
	  if ((norm&1) == 0) res[norm>>1]++;
	}
	break;

      case min_PERF:
      {
        long I=1;
        pari_sp av2=avma;

        for (i=1; i<=n; i++)
          for (j=i; j<=n; j++,I++) V[I] = x[i]*x[j];
        if (! addcolumntomatrix(V,invp,L))
        {
          if (DEBUGLEVEL>1) { fprintferr("."); flusherr(); }
          s--; avma=av2; continue;
        }

        if (DEBUGLEVEL>1) { fprintferr("*"); flusherr(); }
        if (s == maxrank)
        {
          if (DEBUGLEVEL>1) { fprintferr("\n"); flusherr(); }
          avma=av0; return stoi(s);
        }

        if (low_stack(lim, stack_lim(av1,1)))
        {
          if(DEBUGMEM>1) pari_warn(warnmem,"minim0, rank>=%ld",s);
          invp = gerepilecopy(av1, invp);
        }
      }
    }
  }
  switch(flag)
  {
    case min_FIRST:
      avma=av0; return cgetg(1,t_VEC);
    case min_VECSMALL:
    case min_VECSMALL2:
      avma=av; return res;
    case min_PERF:
      if (DEBUGLEVEL>1) { fprintferr("\n"); flusherr(); }
      avma=av0; return stoi(s);
  }
  k = min(s,maxrank);
  r = (maxnorm >= 0) ? ground(dbltor(maxnorm)): BORNE;

  L[0] = evaltyp(t_MAT) | evallg(k + 1);
  for (j=1; j<=k; j++) gel(L,j) = ZM_zc_mul(u, gel(L,j));

  gerepileall(av, 2, &r, &L);
  gel(res,1) = stoi(s<<1);
  gel(res,2) = r;
  gel(res,3) = L; return res;
}

GEN
qfrep0(GEN a, GEN borne, long flag)
{
  pari_sp av = avma;
  GEN g = minim0(a, borne, gen_0, (flag & 1)? min_VECSMALL2: min_VECSMALL);
  if ((flag & 2) == 0) g = gerepileupto(av, gtovec(g));
  return g;
}

GEN
qfminim0(GEN a, GEN borne, GEN stockmax, long flag, long prec)
{
  switch(flag)
  {
    case 0: return minim0(a,borne,stockmax,min_ALL);
    case 1: return minim0(a,borne,gen_0   ,min_FIRST);
    case 2:
    {
      long maxnum = stockmax? itos(stockmax): -2;
      return fincke_pohst(a,borne,maxnum,prec,NULL);
    }
    default: pari_err(flagerr,"qfminim");
  }
  return NULL; /* not reached */
}

GEN
minim(GEN a, GEN borne, GEN stockmax)
{
  return minim0(a,borne,stockmax,min_ALL);
}

GEN
minim2(GEN a, GEN borne, GEN stockmax)
{
  return minim0(a,borne,stockmax,min_FIRST);
}

GEN
perf(GEN a)
{
  return minim0(a,gen_0,gen_0,min_PERF);
}

static GEN
clonefill(GEN S, long s, long t)
{ /* initialize to dummy values */
  GEN T = S, dummy = cgetg(1, t_STR);
  long i;
  for (i = s+1; i <= t; i++) gel(S,i) = dummy;
  S = gclone(S); if (isclone(T)) gunclone(T);
  return S;
}

INLINE void
step(GEN x, GEN y, GEN inc, long k)
{
  if (!signe(y[k]))
    gel(x,k) = addis(gel(x,k), 1); /* leading coeff > 0 */
  else
  {
    long i = inc[k];
    gel(x,k) = addis(gel(x,k), i),
    inc[k] = (i > 0)? -1-i: 1-i;
  }
}
/* q is the Gauss reduction (sqred1) of the quadratic form */
/* general program for positive definit quadratic forms (real coeffs).
 * Enumerate vectors whose norm is less than BORNE, minimal vectors
 * if BORNE = NULL (implies check = NULL).
 * If (check != NULL) consider only vectors passing the check, and assumes
 *   we only want the smallest possible vectors */
static GEN
smallvectors(GEN q, GEN BORNE, long maxnum, FP_chk_fun *CHECK)
{
  long N = lg(q), n = N-1, i, j, k, s, epsbit, prec, checkcnt = 1;
  pari_sp av, av1, lim;
  GEN inc, S, x, y, z, v,  eps, p1, alpha, norms;
  GEN norme1, normax1, borne1, borne2;
  GEN (*check)(void *,GEN) = CHECK? CHECK->f: NULL;
  void *data = CHECK? CHECK->data: NULL;
  long stockmax, skipfirst = CHECK? CHECK->skipfirst: 0;
  int stockall = (maxnum < 0);

  prec = gprecision(q);
  epsbit = bit_accuracy(prec) >> 1;
  eps = real2n(-epsbit, 3);
  alpha = dbltor(0.95);
  normax1 = gen_0;
  norme1 = BORNE ? BORNE: gsqr(gcoeff(q,1,1));
  borne1 = mpadd(norme1,eps);
  if (!BORNE) borne2 = mpsub(norme1,eps);
  else        borne2 = mpmul(norme1,alpha);
  if (DEBUGLEVEL)
    fprintferr("smallvectors looking for norm < %Z\n",gprec_w(borne1,3));

  v = cgetg(N,t_VEC);
  inc = const_vecsmall(n, 1);

  av = avma; lim = stack_lim(av,2);
  stockmax = stockall? 200: maxnum;
  if (check) norms = cgetg(stockmax+1,t_VEC);
  S = cgetg(stockmax+1,t_VEC);
  x = cgetg(N,t_COL);
  y = cgetg(N,t_COL);
  z = cgetg(N,t_COL);
  for (i=1; i<N; i++) {
    gel(v,i) = gcoeff(q,i,i);
    gel(x,i) = gel(y,i) = gel(z,i) = gen_0;
  }

  gel(x,n) = gen_0; s = 0; k = n;
  for(;; step(x,y,inc,k)) /* main */
  {
    do
    {
      int fl = 0;
      if (k > 1)
      {
        k--;
        av1 = avma; p1 = mpmul(gcoeff(q,k,k+1),gel(x,k+1));
	for (j=k+2; j<N; j++)
	  p1 = mpadd(p1, mpmul(gcoeff(q,k,j),gel(x,j)));
        gel(z,k) = gerepileuptoleaf(av1,p1);

        av1 = avma; p1 = gsqr(mpadd(gel(x,k+1),gel(z,k+1)));
        p1 = mpadd(gel(y,k+1), mpmul(p1,gel(v,k+1)));
	gel(y,k) = gerepileuptoleaf(av1, p1);
        /* skip the [x_1,...,x_skipfirst,0,...,0] */
        if ((k <= skipfirst && !signe(y[skipfirst]))
         || mpcmp(borne1, gel(y,k)) < 0) fl = 1;
        else
          gel(x,k) = ground( mpneg(gel(z,k)) );
      }
      for(;; step(x,y,inc,k))
      {
	if (!fl)
	{
          av1 = avma;
	  p1 = mpmul(gel(v,k), gsqr(mpadd(gel(x,k), gel(z,k))));
	  i = mpcmp(mpsub(mpadd(p1,gel(y,k)), borne1), gmul2n(p1,-epsbit));
          avma = av1; if (i <= 0) break;

          step(x,y,inc,k);

          av1 = avma; /* same as above */
	  p1 = mpmul(gel(v,k), gsqr(mpadd(gel(x,k), gel(z,k))));
	  i = mpcmp(mpsub(mpadd(p1,gel(y,k)), borne1), gmul2n(p1,-epsbit));
          avma = av1; if (i <= 0) break;
	}
        fl = 0; inc[k] = 1;
        if (++k > n) goto END;
      }

      if (low_stack(lim, stack_lim(av,2)))
      {
	if(DEBUGMEM>1) pari_warn(warnmem,"smallvectors");
	if (stockmax) S = clonefill(S, s, stockmax);
        if (check) {
          GEN dummy = cgetg(1, t_STR);
          for (i=s+1; i<=stockmax; i++) gel(norms,i) = dummy;
        }
	gerepileall(av,check?7:6,&x,&y,&z,&normax1,&borne1,&borne2,&norms);
      }
    }
    while (k > 1);
    if (!signe(x[1]) && !signe(y[1])) continue; /* exclude 0 */

    av1 = avma; p1 = gsqr(mpadd(gel(x,1),gel(z,1)));
    norme1 = mpadd(gel(y,1), mpmul(p1, gel(v,1)));
    if (mpcmp(norme1,borne1) > 0) { avma=av1; continue; /* main */ }

    norme1 = gerepileupto(av1,norme1);
    if (check)
    {
      if (checkcnt < 5 && mpcmp(norme1, borne2) < 0)
      {
        if (!check(data,x)) { checkcnt++ ; continue; /* main */}
        if (DEBUGLEVEL>4) fprintferr("New bound: %Z", norme1);
        borne1 = mpadd(norme1, eps);
        borne2 = mpmul(borne1, alpha);
        s = 0; checkcnt = 0;
      }
    }
    else
    {
      if (!BORNE) /* find minimal vectors */
      {
        if (mpcmp(norme1, borne2) < 0)
        {
          borne1 = mpadd(norme1, eps);
          borne2 = mpsub(norme1, eps);
          s = 0; 
        }
      }
      else
        if (mpcmp(norme1,normax1) > 0) normax1 = norme1;
    }

    if (++s <= stockmax)
    {
      if (check) gel(norms,s) = norme1;
      gel(S,s) = shallowcopy(x);
      if (s == stockmax)
      { /* overflow */
        long stockmaxnew= (stockall && (stockmax < 10000L || maxnum != -1))
                          ? stockmax<<1 : stockmax;
        GEN Snew = cgetg(stockmaxnew + 1, t_VEC);
        if (check)
        {
          pari_sp av2 = avma;
          GEN per = sindexsort(norms);
          if (DEBUGLEVEL) fprintferr("sorting...\n");
          for (j = 0, i = 1; i <= s; i++)
          { /* let N be the minimal norm so far for x satisfying 'check'. Keep
             * all elements of norm N */
            long k = per[i];
            norme1 = gel(norms,k);
            if (j  && mpcmp(norme1, borne1) > 0) break;
            if (j  || check(data,gel(S,k)))
            {
              if (!j) borne1 = mpadd(norme1,eps);
              Snew[++j] = S[k];
            }
          }
          s = j; avma = av2;
          if (s)
          {
            norme1 = gel(norms, per[i-1]);
            borne1 = mpadd(norme1, eps);
            borne2 = mpmul(borne1, alpha);
            checkcnt = 0;
          }
        }
        else
        {
          if (!stockall && BORNE) goto END;
          for (i = 1; i <= s; i++) Snew[i] = S[i];
        }
        if (stockmax != stockmaxnew)
        {
          stockmax = stockmaxnew;
          norms = cgetg(stockmax+1, t_VEC);
          for (i = 1; i <= s; i++) gel(norms,i) = norme1;
          Snew = clonefill(Snew, s, stockmax);
          if (isclone(S)) gunclone(S);
          S = Snew;
        }
      }
    }
  }
END:
  if (s < stockmax) stockmax = s;
  if (check)
  {
    GEN per, alph, pols, p;
    if (DEBUGLEVEL) fprintferr("final sort & check...\n");
    setlg(norms,stockmax+1); per = sindexsort(norms);
    alph = cgetg(stockmax+1,t_VEC);
    pols = cgetg(stockmax+1,t_VEC);
    for (j=0,i=1; i<=stockmax; i++)
    {
      long t = per[i];
      norme1 = gel(norms,t);
      if (j && mpcmp(norme1, borne1) > 0) break;
      if ((p = check(data,gel(S,t))))
      {
        if (!j) borne1 = mpadd(norme1,eps);
        j++; gel(pols,j) = p; alph[j]=S[t];
      }
    }
    setlg(pols,j+1);
    setlg(alph,j+1);
    if (stockmax && isclone(S)) { alph = gcopy(alph); gunclone(S); }
    return mkvec2(pols, alph);
  }
  if (stockmax)
  {
    setlg(S,stockmax+1);
    settyp(S,t_MAT);
    if (isclone(S)) { p1 = S; S = gcopy(S); gunclone(p1); }
  }
  else
    S = cgetg(1,t_MAT);
  if (!BORNE) normax1 = mpsub(borne1, eps);
  return mkvec3(utoi(s<<1), normax1, S);
}

/* solve q(x) = x~.a.x <= bound, a > 0.
 * If check is non-NULL keep x only if check(x).
 * If a is a vector, assume a[1] is the LLL-reduced Cholesky form of q */
GEN
fincke_pohst(GEN a, GEN B0, long stockmax, long PREC, FP_chk_fun *CHECK)
{
  pari_sp av = avma;
  VOLATILE long i,j,l;
  VOLATILE GEN r,rinvtrans,u,v,res,z,vnorm,rperm,perm,uperm, bound = B0;

  if (typ(a) == t_VEC)
  {
    r = gel(a,1);
    u = NULL;
  }
  else
  {
    long prec = PREC;
    l = lg(a);
    if (l == 1)
    {
      if (CHECK) pari_err(talker, "dimension 0 in fincke_pohst");
      z = cgetg(4,t_VEC);
      gel(z,1) = gel(z,2) = gen_0;
      gel(z,3) = cgetg(1,t_MAT); return z;
    }
    i = gprecision(a); if (i) prec = i;
    if (DEBUGLEVEL>2) fprintferr("first LLL: prec = %ld\n", prec);
    u = lllgramintern(a, 4, 1, (prec<<1)-2);
    if (!u) return NULL;
    r = qf_base_change(a,u,1);
    if (!i) {
      prec = DEFAULTPREC + nbits2nlong(gexpo(r));
      if (prec < PREC) prec = PREC;
    }
    r = sqred1intern(r);
    if (!r) return NULL;
    for (i=1; i<l; i++)
    {
      GEN s = gsqrt(gcoeff(r,i,i), prec);
      gcoeff(r,i,i) = s;
      for (j=i+1; j<l; j++) gcoeff(r,i,j) = gmul(s, gcoeff(r,i,j));
    }
  }
  /* now r~ * r = a in LLL basis */
  rinvtrans = shallowtrans( invmat(r) );
  if (DEBUGLEVEL>2)
    fprintferr("Fincke-Pohst, final LLL: prec = %ld\n", gprecision(rinvtrans));
  v = lllintern(rinvtrans, 100, 1, 0);
  if (!v) return NULL;

  rinvtrans = gmul(rinvtrans, v);
  v = ZM_inv(shallowtrans(v),gen_1);
  r = gmul(r,v);
  u = u? gmul(u,v): v;

  l = lg(r);
  vnorm = cgetg(l,t_VEC);
  for (j=1; j<l; j++) gel(vnorm,j) = gnorml2(gel(rinvtrans,j));
  rperm = cgetg(l,t_MAT);
  uperm = cgetg(l,t_MAT); perm = sindexsort(vnorm);
  for (i=1; i<l; i++) { uperm[l-i] = u[perm[i]]; rperm[l-i] = r[perm[i]]; }
  u = uperm;
  r = rperm; res = NULL;
  CATCH(precer) { }
  TRY {
    if (CHECK && CHECK->f_init) bound = CHECK->f_init(CHECK, r, u);
    r = sqred1_from_QR(r, gprecision(vnorm));
    if (!r) pari_err(precer,"fincke_pohst");
    res = smallvectors(r, bound, stockmax, CHECK);
  } ENDCATCH;
  if (CHECK)
  {
    if (CHECK->f_post) res = CHECK->f_post(CHECK, res, u);
    return res;
  }
  if (!res) pari_err(precer,"fincke_pohst");

  z = cgetg(4,t_VEC);
  gel(z,1) = gcopy(gel(res,1));
  gel(z,2) = gcopy(gel(res,2));
  gel(z,3) = gmul(u, gel(res,3)); return gerepileupto(av,z);
}
