/* $Id: base5.c 7838 2006-04-08 12:11:17Z kb $

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
/*                       BASIC NF OPERATIONS                       */
/*                          (continued 2)                          */
/*                                                                 */
/*******************************************************************/
#include "pari.h"
#include "paripriv.h"

static GEN
_checkrnfeq(GEN x)
{
  if (typ(x) == t_VEC)
    switch(lg(x))
    {
      case 13: /* checkrnf(x); */ return gel(x,11);
      case  4: return x;
    }
  return NULL;
}

GEN
checkrnfeq(GEN x)
{
  x = _checkrnfeq(x);
  if (!x) pari_err(talker,"please apply rnfequation(,,1)");
  return x;
}

GEN
eltreltoabs(GEN rnfeq, GEN x)
{
  long i, k, va;
  pari_sp av = avma;
  GEN polabs, teta, alpha, s;

  rnfeq = checkrnfeq(rnfeq);
  polabs= gel(rnfeq,1);
  alpha = lift_intern(gel(rnfeq,2));
  k     = itos(gel(rnfeq,3));

  va = varn(polabs);
  if (varncmp(gvar(x), va) > 0) x = scalarpol(x,va);
  /* Mod(X - k alpha, polabs(X)), alpha root of the polynomial defining base */
  teta = gadd(pol_x[va], gmulsg(-k,alpha));
  s = gen_0;
  for (i=lg(x)-1; i>1; i--)
  {
    GEN c = gel(x,i);
    long tc = typ(c);
    switch(tc)
    {
      case t_POLMOD: c = gel(c,2); /* fall through */
      case t_POL:    c = RgX_RgXQ_compo(c, alpha, polabs); break;
      default:
        if (!is_const_t(tc)) pari_err(talker, "incorrect data in eltreltoabs");
    }
    s = RgX_rem(gadd(c, gmul(teta,s)), polabs);
  }
  return gerepileupto(av, s);
}

#if 0
static GEN
rnfmakematrices(GEN rnf)
{
  long i, j, k, n, ru, vnf;
  GEN nf, pol, sym, ro, w, ronf, z, vecM, T;

  pol = gel(rnf,1); n = degpol(pol); sym = polsym(pol, n-1);
  ro  = gel(rnf,6);
  w   = gel(rnf,7); w = gel(w,1);
  nf  = gel(rnf,10); vnf = varn(nf[1]);
  T = cgetg(n+1,t_MAT);
  for (j=1; j<=n; j++)
  {
    GEN c = cgetg(n+1,t_COL); gel(T,j) = c;
    for (i=1; i<=n; i++)
    {
      GEN d = grem(gmul(gel(w,i),gel(w,j)), pol);
      gel(c,i) = lift_if_rational( quicktrace(d, sym) );
    }
  }
  w = lift(w); ru = lg(ro)-1; ronf = gel(nf,6);
  vecM = cgetg(ru+1,t_VEC);
  for (k=1; k<=ru; k++)
  {
    GEN rok = gel(ro,k), M = cgetg(n+1,t_MAT);
    long l = lg(rok);
    gel(vecM,k) = M;
    for (j=1; j<=n; j++)
    {
      GEN a, c = cgetg(l,t_COL); gel(M,j) = c;
      a = gsubst(gel(w,j), vnf, gel(ronf,k));
      for (i=1; i<l; i++) gel(c,i) = poleval(a, gel(rok,i));
    }
  }

  z = cgetg(8,t_VEC);
  gel(z,1) = vecM;
  gel(z,4) = T;
  /* dummies */
  gel(z,2) = cgetg(1, t_VEC);
  gel(z,3) = cgetg(1, t_VEC);
  gel(z,5) = cgetg(1,t_MAT);
  gel(z,6) = cgetg(1,t_MAT);
  gel(z,7) = cgetg(1,t_MAT); return z;
}

static GEN
rnf_roots(GEN nf, GEN pol, long prec, GEN *pts)
{
  long r1, r2, j, v = varn(nf[1]), n = degpol(pol);
  GEN s, r, ro;

  nf_get_sign(nf, &r1, &r2);
  s = cgetg(r1+r2+1,t_VEC);
  r = cgetg(r1+r2+1,t_VEC);
  ro = gel(nf,6); pol = lift(pol);
  for (j=1; j<=r1; j++)
  {
    long r1j = 0;
    ro = roots(gsubst(pol,v,gel(ro,j)), prec);
    while (r1j<n && gcmp0(imag_i(gel(ro,r1j+1)))) r1j++;
    gel(s,j) = mkvec2s(r1j, (n-r1j)>>1);
    gel(r,j) = get_roots(ro, r1j, 0);
  }
  for (; j<=r1+r2; j++)
  {
    ro = roots(gsubst(pol,v,gel(ro,j)), prec);
    gel(s,j) = mkvec2s(0, n);
    gel(r,j) = ro;
  }
  *pts = s; return r;
}

#else /* dummies */
static GEN rnfmakematrices(GEN rnf) { (void)rnf; return cgetg(1, t_VEC); }
static GEN
rnf_roots(GEN nf, GEN pol, long prec, GEN *pts) {
  (void)nf; (void)pol; (void)prec;
  *pts = cgetg(1,t_VEC); return cgetg(1, t_VEC);
}
#endif

static GEN
modulereltoabs(GEN rnf, GEN x)
{
  GEN w = gel(x,1), I = gel(x,2), nf = gel(rnf,10), rnfeq = gel(rnf,11);
  GEN M, basnf, cobasnf, T = gel(nf,1), polabs = gel(rnfeq,1);
  long i, j, k, n = lg(w)-1, m = degpol(T);

  M = cgetg(n*m+1, t_VEC);
  basnf = lift_intern( gsubst(gel(nf,7), varn(T), gel(rnfeq,2)) );
  basnf = Q_primitive_part(basnf, &cobasnf); /* remove denom. --> faster */
  for (k=i=1; i<=n; i++)
  {
    GEN c0, om = gel(w,i), id = gel(I,i);

    om = Q_primitive_part(eltreltoabs(rnfeq, om), &c0);
    c0 = mul_content(c0, cobasnf);
    for (j=1; j<=m; j++)
    {
      GEN c, z = Q_primitive_part(gmul(basnf,gel(id,j)), &c);
      z = RgX_rem(gmul(om, RgX_rem(z,polabs)), polabs);
      c = mul_content(c, c0); if (c) z = gmul(c, z);
      gel(M,k++) = z;
    }
  }
  return M;
}

GEN
hnfcenter_ip(GEN M)
{
  long i, j, k, N = lg(M)-1;
  GEN a, Mj, Mk;

  for (j=N-1; j>0; j--)
  {
    Mj = gel(M,j); a = gel(Mj,j);
    if (cmpiu(a, 2) <= 0) continue;
    a = shifti(a, -1);
    for (k = j+1; k <= N; k++)
    {
      Mk = gel(M,k);
      if (cmpii(gel(Mk,j),a) > 0)
        for (i = 1; i <= j; i++) gel(Mk,i) = subii(gel(Mk,i), gel(Mj,i));
    }
  }
  return M;
}

static GEN
makenfabs(GEN rnf)
{
  GEN M, d, rnfeq, pol, nf, NF = zerovec(9);
  long n;

  rnfeq = gel(rnf,11); pol = gel(rnfeq,1);
  nf = gel(rnf,10);

  M = modulereltoabs(rnf, gel(rnf,7));
  n = lg(M)-1;
  M = RgXV_to_RgM(Q_remove_denom(M, &d), n);
  if (d) M = gdiv(hnfcenter_ip(hnfmodid(M, d)), d);
  else   M = matid(n);

  gel(NF,1) = pol;
  gel(NF,3) = mulii(powiu(gel(nf,3), degpol(rnf[1])),
                      idealnorm(nf, gel(rnf,3)));
  gel(NF,7) = RgM_to_RgXV(M,varn(pol));
  gel(NF,8) = invmat(M);
  gel(NF,9) = get_mul_table(pol, gel(NF,7), gel(NF,8));
  /* possibly wrong, but correct prime divisors [for primedec] */
  gel(NF,4) = Q_denom(gel(NF,7));
  return NF;
}

static GEN
makenorms(GEN rnf)
{
  GEN f = gel(rnf,4);
  return typ(f) == t_INT? gen_1: dethnf(f);
}

#define NFABS 1
#define NORMS 2
GEN
check_and_build_nfabs(GEN rnf) {
  return check_and_build_obj(rnf, NFABS, &makenfabs);
}
GEN
check_and_build_norms(GEN rnf) {
  return check_and_build_obj(rnf, NORMS, &makenorms);
}

GEN
rnfinitalg(GEN nf, GEN pol, long prec)
{
  pari_sp av = avma;
  long vpol;
  GEN rnf, delta, bas, D,d,f, B;

  if (typ(pol)!=t_POL) pari_err(notpoler,"rnfinitalg");
  nf = checknf(nf); vpol = varn(pol);
  pol = fix_relative_pol(nf,pol,0);
  if (vpol >= varn(nf[1]))
    pari_err(talker,"main variable must be of higher priority in rnfinitalg");

  bas = rnfallbase(nf,pol, &D,&d, &f);
  B = matbasistoalg(nf,gel(bas,1));
  gel(bas,1) = lift_if_rational( RgM_to_RgXV(B,vpol) );
  delta = mkvec2(D, d);

  rnf = cgetg(13, t_VEC);
  gel(rnf,1) = pol;
  gel(rnf,3) = delta;
  gel(rnf,4) = f;
  gel(rnf,6) = rnf_roots(nf, pol, prec, (GEN*)rnf+2);
  gel(rnf,7) = bas;
  gel(rnf,8) = lift_if_rational( invmat(B) );
  gel(rnf,9) = cgetg(1,t_VEC); /* dummy */
  gel(rnf,10) = nf;
  gel(rnf,11) = rnfequation2(nf,pol);
  gel(rnf,12) = gen_0;
  gel(rnf,5) = rnfmakematrices(rnf);
  return gerepilecopy(av, rnf);
}

GEN
rnfelementreltoabs(GEN rnf,GEN x)
{
  long i, lx, tx = typ(x);
  GEN z;

  switch(tx)
  {
    case t_VEC: case t_COL: case t_MAT:
      lx = lg(x); z = cgetg(lx,tx);
      for (i=1; i<lx; i++) gel(z,i) = rnfelementreltoabs(rnf, gel(x,i));
      return z;

    case t_POLMOD: x = lift_to_pol(x); /* fall through */
    case t_POL: return eltreltoabs(rnf, x);
    default: return gcopy(x);
  }
}

/* assume x,T,pol t_POL. T defines base field, pol defines rnf over T.
 * x an absolute element of the extension */
GEN
get_theta_abstorel(GEN T, GEN pol, GEN k)
{
  return mkpolmod(gadd(pol_x[varn(pol)],
                       gmul(k, mkpolmod(pol_x[varn(T)],T))), pol);
}
GEN
eltabstorel(GEN x, GEN T, GEN pol, GEN k)
{
  return poleval(x, get_theta_abstorel(T,pol,k));
}
static GEN
rnf_get_theta_abstorel(GEN rnf)
{
  GEN k, nf, T, pol, rnfeq;
  rnfeq  = gel(rnf,11); k = gel(rnfeq,3);
  nf = gel(rnf,10); T = gel(nf,1);
  pol = gel(rnf,1);
  return get_theta_abstorel(T, pol, k);
}

GEN
rnfelementabstorel(GEN rnf,GEN x)
{
  pari_sp av = avma;
  long tx, i, lx;
  GEN z;

  checkrnf(rnf); tx = typ(x);
  switch(tx)
  {
    case t_VEC: case t_COL: case t_MAT:
      lx = lg(x); z = cgetg(lx,tx);
      for (i=1; i<lx; i++) gel(z,i) = rnfelementabstorel(rnf,gel(x,i));
      return z;

    case t_POLMOD:
      x = lift_to_pol(x); /* fall through */
    case t_POL:
      return gerepileupto(av, poleval(x, rnf_get_theta_abstorel(rnf)));

    default: return gcopy(x);
  }
}

/* x t_POLMOD or t_POL or vector of such objects */
GEN
rnfelementup(GEN rnf,GEN x)
{
  long i, lx, tx;
  GEN z;

  checkrnf(rnf); tx = typ(x);
  switch(tx)
  {
    case t_VEC: case t_COL: case t_MAT:
      lx = lg(x); z = cgetg(lx,tx);
      for (i=1; i<lx; i++) gel(z,i) = rnfelementup(rnf,gel(x,i));
      return z;

    case t_POLMOD: x = gel(x,2); /* fall through */
    case t_POL:
      return poleval(x, gmael(rnf,11,2));

    default: return gcopy(x);
  }
}

/* x t_POLMOD or t_POL or vector of such objects */
GEN
rnfelementdown(GEN rnf,GEN x)
{
  pari_sp av;
  long i, lx, tx;
  GEN z;

  checkrnf(rnf); tx = typ(x);
  switch(tx)
  {
    case t_VEC: case t_COL: case t_MAT:
      lx = lg(x); z = cgetg(lx,tx);
      for (i=1; i<lx; i++) gel(z,i) = rnfelementdown(rnf,gel(x,i));
      return z;

    case t_POLMOD: x = gel(x,2); /* fall through */
    case t_POL:
      if (gcmp0(x)) return gen_0;
      av = avma; z = rnfelementabstorel(rnf,x);
      if (typ(z)==t_POLMOD && varn(z[1])==varn(rnf[1])) z = gel(z,2);
      if (varncmp(gvar(z), varn(rnf[1])) <= 0)
      {
        lx = lg(z);
        if (lx == 2) { avma = av; return gen_0; }
        if (lx > 3)
          pari_err(talker,"element is not in the base field in rnfelementdown");
        z = gel(z,2);
      }
      return gerepilecopy(av, z);

    default: return gcopy(x);
  }
}

static GEN
rnfid(long n, long m)
{
  return matid_intern(n, col_ei(m,1), zerocol(m));
}

/* x est exprime sur la base relative */
static GEN
rnfprincipaltohermite(GEN rnf,GEN x)
{
  pari_sp av = avma;
  GEN bas = gel(rnf,7), nf = gel(rnf,10);

  x = rnfbasistoalg(rnf,x);
  x = rnfalgtobasis(rnf, gmul(x, gmodulo(gel(bas,1), gel(rnf,1))));
  settyp(x, t_MAT);
  return gerepileupto(av, nfhermite(nf, mkvec2(x, gel(bas,2))));
}

GEN
rnfidealhermite(GEN rnf, GEN x)
{
  GEN z, nf, bas;

  checkrnf(rnf); nf = gel(rnf,10);
  switch(typ(x))
  {
    case t_INT: case t_FRAC:
      bas = gel(rnf,7); z = cgetg(3,t_VEC);
      gel(z,1) = rnfid(degpol(rnf[1]), degpol(nf[1]));
      gel(z,2) = gmul(x, gel(bas,2)); return z;

    case t_VEC:
      if (lg(x) == 3 && typ(x[1]) == t_MAT) return nfhermite(nf, x);
      return rnfidealabstorel(rnf, x);

    case t_POLMOD: case t_POL: case t_COL:
      return rnfprincipaltohermite(rnf,x);
  }
  pari_err(typeer,"rnfidealhermite");
  return NULL; /* not reached */
}

GEN
prodid(GEN nf, GEN I)
{
  long i, l = lg(I);
  GEN z;
  if (l == 1) return matid(degpol(nf[1]));
  z = gel(I,1);
  for (i=2; i<l; i++) z = idealmul(nf, z, gel(I,i));
  return z;
}

static GEN
prodidnorm(GEN I)
{
  long i, l = lg(I);
  GEN z;
  if (l == 1) return gen_1;
  z = dethnf(gel(I,1));
  for (i=2; i<l; i++) z = gmul(z, dethnf(gel(I,i)));
  return z;
}

GEN
rnfidealnormrel(GEN rnf, GEN id)
{
  pari_sp av = avma;
  GEN z, nf = gel(rnf,10);

  checkrnf(rnf);
  if (degpol(rnf[1]) == 1) return matid(degpol(nf[1]));

  z = prodid(nf, (GEN)rnfidealhermite(rnf,id)[2]);
  return gerepileupto(av, idealmul(nf,z, gel(rnf,4)));
}

GEN
rnfidealnormabs(GEN rnf, GEN id)
{
  pari_sp av = avma;
  GEN z;

  checkrnf(rnf);
  if (degpol(rnf[1]) == 1) return gen_1;

  z = prodidnorm( (GEN)rnfidealhermite(rnf,id)[2] );
  return gerepileupto(av, gmul(z, check_and_build_norms(rnf)));
}

GEN
rnfidealreltoabs(GEN rnf,GEN x)
{
  pari_sp av = avma;
  long i, l;
  GEN w;

  x = rnfidealhermite(rnf,x);
  w = gel(x,1); l = lg(w); settyp(w, t_VEC);
  for (i=1; i<l; i++) gel(w,i) = lift( rnfbasistoalg(rnf, gel(w,i)) );
  return gerepilecopy(av, modulereltoabs(rnf, x));
}

GEN
rnfidealabstorel(GEN rnf, GEN x)
{
  long N, m, j;
  pari_sp av = avma;
  GEN nf, A, I, z, id, invbas;

  checkrnf(rnf); nf = gel(rnf,10); invbas = gel(rnf,8);
  m = degpol(nf[1]);
  N = m * degpol(rnf[1]);
  if (lg(x)-1 != N) pari_err(typeer, "rnfidealabstorel");
  if (typ(x) != t_VEC) pari_err(typeer,"rnfidealabstorel");
  A = cgetg(N+1,t_MAT);
  I = cgetg(N+1,t_VEC); z = mkvec2(A,I); id = matid(m);
  for (j=1; j<=N; j++)
  {
    GEN t = lift_intern( rnfelementabstorel(rnf, gel(x,j)) );
    gel(A,j) = mulmat_pol(invbas, t);
    gel(I,j) = id;
  }
  return gerepileupto(av, nfhermite(nf,z));
}

GEN
rnfidealdown(GEN rnf,GEN x)
{
  pari_sp av = avma; x = rnfidealhermite(rnf,x);
  return gerepilecopy(av, gmael(x,2,1));
}

/* lift ideal x to the relative extension, returns a Z-basis */
GEN
rnfidealup(GEN rnf,GEN x)
{
  pari_sp av = avma;
  long i, n;
  GEN nf, bas, bas2, I, z;

  checkrnf(rnf); nf = gel(rnf,10);
  n = degpol(rnf[1]);
  bas = gel(rnf,7); bas2 = gel(bas,2);

  (void)idealtyp(&x, &z); /* z is junk */
  I = cgetg(n+1,t_VEC); z = mkvec2(gel(bas,1), I);
  for (i=1; i<=n; i++) gel(I,i) = idealmul(nf,x,gel(bas2,i));
  return gerepilecopy(av, modulereltoabs(rnf, z));
}

/* x a relative HNF ---> vector of 2 generators (relative polymods) */
GEN
rnfidealtwoelement(GEN rnf, GEN x)
{
  pari_sp av = avma;
  GEN y, z, NF;

  checkrnf(rnf);
  NF = check_and_build_nfabs(rnf);
  y = rnfidealreltoabs(rnf,x);
  y = algtobasis(NF, y); settyp(y, t_MAT);
  y = ideal_two_elt(NF, hnf(y));
  z = rnfelementabstorel(rnf, gmul(gel(NF,7), gel(y,2)));
  return gerepilecopy(av, mkvec2(gel(y,1), z));
}

GEN
rnfidealmul(GEN rnf,GEN x,GEN y) /* x et y sous HNF relative uniquement */
{
  pari_sp av = avma;
  GEN z, nf, x1, x2, p1, p2;

  z = rnfidealtwoelement(rnf,y);
  nf = gel(rnf,10);
  x = rnfidealhermite(rnf,x);
  x1 = gmodulo(gmul(gmael(rnf,7,1), matbasistoalg(nf,gel(x,1))),gel(rnf,1));
  x2 = gel(x,2);
  p1 = gmul(gel(z,1), gel(x,1));
  p2 = rnfalgtobasis(rnf, gmul(gel(z,2), x1)); settyp(p2, t_MAT);
  z = mkvec2(shallowconcat(p1, p2), shallowconcat(x2, x2));
  return gerepileupto(av, nfhermite(nf,z));
}
