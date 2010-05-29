/* $Id: Qfb.c 9113 2007-10-29 09:35:40Z kb $

Copyright (C) 2000-2005  The PARI group.

This file is part of the PARI/GP package.

PARI/GP is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation. It is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY WHATSOEVER.

Check the License for details. You should have received a copy of it, along
with the package; see the file 'COPYING'. If not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA. */
#include "pari.h"
#include "paripriv.h"
/*******************************************************************/
/*                                                                 */
/*         QUADRATIC POLYNOMIAL ASSOCIATED TO A DISCRIMINANT       */
/*                                                                 */
/*******************************************************************/

void
check_quaddisc(GEN x, long *s, long *r, char *f)
{
  if (typ(x) != t_INT) pari_err(arither1);
  *s = signe(x); if (!*s) pari_err(talker,"zero discriminant in %s", f);
  if (Z_issquare(x)) pari_err(talker,"square discriminant in %s", f);
  *r = mod4(x); if (*s < 0 && *r) *r = 4 - *r;
  if (*r > 1) pari_err(talker, "discriminant not congruent to 0,1 mod 4 in %s", f);
}
void
check_quaddisc_real(GEN x, long *r, char *f)
{
  long sx; check_quaddisc(x, &sx, r, f);
  if (sx < 0) pari_err(talker, "negative discriminant in %s", f);
}
void
check_quaddisc_imag(GEN x, long *r, char *f)
{
  long sx; check_quaddisc(x, &sx, r, f);
  if (sx > 0) pari_err(talker, "positive discriminant in %s", f);
}

static GEN
Zquadpoly(GEN x, long v)
{
  long res, sx;
  GEN y, p1;

  check_quaddisc(x, &sx, &res, "quadpoly");
  y = cgetg(5,t_POL);
  y[1] = evalsigne(1) | evalvarn(v);

  p1 = shifti(x,-2); togglesign(p1);
  /* p1 = - floor(x/4) [ = -x/4 or (1-x)/4 ] */
  if (!res) gel(y,3) = gen_0;
  else
  {
    if (sx < 0) p1 = gerepileuptoint((pari_sp)y, addsi(1,p1));
    gel(y,3) = gen_m1;
  }
  gel(y,2) = p1;
  gel(y,4) = gen_1; return y;
}

GEN
quadpoly0(GEN x, long v)
{
  long tx = typ(x);
  if (is_matvec_t(tx))
  {
    long i, l = lg(x);
    GEN y = cgetg(l, tx);
    for (i=1; i<l; i++) gel(y,i) = quadpoly0(gel(x,i),v);
    return y;
  }
  if (v < 0) v = 0;
  return Zquadpoly(x,v);
}

GEN
quadpoly(GEN x) { return quadpoly0(x, -1); }

GEN
quadgen(GEN x)
{
  GEN y = cgetg(4,t_QUAD);
  gel(y,1) = Zquadpoly(x,0); gel(y,2) = gen_0; gel(y,3) = gen_1; return y;
}

/***********************************************************************/
/**                                                                   **/
/**                      BINARY QUADRATIC FORMS                       **/
/**                                                                   **/
/***********************************************************************/

static GEN
qf_disc0(GEN x, GEN y, GEN z) { return subii(sqri(y), shifti(mulii(x,z),2)); }
GEN
qf_disc(GEN x) { return qf_disc0(gel(x,1), gel(x,2), gel(x,3)); }

GEN
qfi(GEN x, GEN y, GEN z)
{
  GEN t = cgetg(4,t_QFI);
  if (signe(x) < 0) pari_err(impl,"negative definite t_QFI");
  gel(t,1) = icopy(x);
  gel(t,2) = icopy(y);
  gel(t,3) = icopy(z); return t;
}
GEN
qfr(GEN x, GEN y, GEN z, GEN d)
{
  GEN t = cgetg(5,t_QFR);
  if (typ(d) != t_REAL) pari_err(talker,"Shanks distance must be a t_REAL in qfr");
  gel(t,1) = icopy(x);
  gel(t,2) = icopy(y);
  gel(t,3) = icopy(z);
  gel(t,4) = rcopy(d); return t;
}

GEN
Qfb0(GEN x, GEN y, GEN z, GEN d, long prec)
{
  pari_sp av = avma;
  long s;
  if (typ(x)!=t_INT || typ(y)!=t_INT || typ(z)!=t_INT) pari_err(typeer,"Qfb");
  s = signe(qf_disc0(x,y,z)); avma = av;
  if (!s) pari_err(talker,"zero discriminant in Qfb");
  if (s < 0) return qfi(x, y, z);

  d = d? gtofp(d,prec): real_0(prec);
  return qfr(x,y,z,d);
}

/* composition */
static void
qfb_sqr(GEN z, GEN x)
{
  GEN c, d1, x2, y2, v1, v2, c3, m, p1, r;

  d1 = bezout(gel(x,2),gel(x,1),&x2,&y2);
  c = gel(x,3);
  m = mulii(c,x2);
  if (is_pm1(d1))
    v1 = v2 = gel(x,1);
  else
  {
    v1 = diviiexact(gel(x,1),d1);
    v2 = mulii(v1, gcdii(d1,c)); /* = v1 iff x primitive */
    c = mulii(c, d1);
  }
  setsigne(m, -signe(m));
  r = modii(m,v2);
  p1 = mulii(r, v1);
  c3 = addii(c, mulii(r,addii(gel(x,2),p1)));
  gel(z,1) = mulii(v1,v2);
  gel(z,2) = addii(gel(x,2), shifti(p1,1));
  gel(z,3) = diviiexact(c3,v2);
}
/* z <- x * y */
void
qfb_comp(GEN z, GEN x, GEN y)
{
  GEN s, n, c, d, y1, v1, v2, c3, m, p1, r;

  if (x == y) { qfb_sqr(z,x); return; }
  s = shifti(addii(gel(x,2),gel(y,2)), -1);
  n = subii(gel(y,2),s);
  v1 = gel(x,1);
  v2 = gel(y,1);
  c  = gel(y,3);
  d = bezout(v2,v1,&y1,NULL);
  if (is_pm1(d))
    m = mulii(y1,n);
  else
  {
    GEN x2, y2, d1 = bezout(s,d,&x2,&y2); /* x2 s + y2 (x1 v1 + y1 v2) = d1 */
    if (!is_pm1(d1))
    {
      v1 = diviiexact(v1,d1);
      v2 = diviiexact(v2,d1); /* gcd = 1 iff x or y primitive */
      v1 = mulii(v1, gcdii(c,gcdii(gel(x,3),gcdii(d1,n))));
      c = mulii(c, d1);
    }
    m = addii(mulii(mulii(y1,y2),n), mulii(gel(y,3),x2));
  }
  setsigne(m, -signe(m));
  r = modii(m, v1);
  p1 = mulii(r, v2);
  c3 = addii(c, mulii(r,addii(gel(y,2),p1)));
  gel(z,1) = mulii(v1,v2);
  gel(z,2) = addii(gel(y,2), shifti(p1,1));
  gel(z,3) = dvmdii(c3,v1, &s);
  if (signe(s)) pari_err(talker,"different discriminants in qfb_comp");
}

static GEN
compimag0(GEN x, GEN y, int raw)
{
  pari_sp av = avma;
  long tx = typ(x);
  GEN z = cgetg(4,t_QFI);
  if (typ(y) != tx || tx != t_QFI) pari_err(typeer,"composition");
  if (absi_cmp(gel(x,1), gel(y,1)) > 0) swap(x, y);
  qfb_comp(z, x,y);
  if (raw) return gerepilecopy(av,z);
  return gerepileupto(av, redimag(z));
}
static GEN
compreal0(GEN x, GEN y, int raw)
{
  pari_sp av = avma;
  long tx = typ(x);
  GEN z = cgetg(5,t_QFR);
  if (typ(y) != tx || tx != t_QFR) pari_err(typeer,"composition");
  qfb_comp(z, x,y); gel(z,4) = addrr(gel(x,4),gel(y,4));
  if (raw) return gerepilecopy(av,z);
  return gerepileupto(av, redreal(z));
}
GEN
compreal(GEN x, GEN y) { return compreal0(x,y,0); }
GEN
comprealraw(GEN x, GEN y) { return compreal0(x,y,1); }
GEN
compimag(GEN x, GEN y) { return compimag0(x,y,0); }
GEN
compimagraw(GEN x, GEN y) { return compimag0(x,y,1); }
GEN
compraw(GEN x, GEN y)
{ return (typ(x)==t_QFI)? compimagraw(x,y): comprealraw(x,y); }

static GEN
sqcompimag0(GEN x, long raw)
{
  pari_sp av = avma;
  GEN z = cgetg(4,t_QFI);

  if (typ(x)!=t_QFI) pari_err(typeer,"composition");
  qfb_sqr(z,x);
  if (raw) return gerepilecopy(av,z);
  return gerepileupto(av, redimag(z));
}
static GEN
sqcompreal0(GEN x, long raw)
{
  pari_sp av = avma;
  GEN z = cgetg(5,t_QFR);

  if (typ(x)!=t_QFR) pari_err(typeer,"composition");
  qfb_sqr(z,x); gel(z,4) = shiftr(gel(x,4),1);
  if (raw) return gerepilecopy(av,z);
  return gerepileupto(av, redreal(z));
}
GEN
sqcompreal(GEN x) { return sqcompreal0(x,0); }
GEN
sqcomprealraw(GEN x) { return sqcompreal0(x,1); }
GEN
sqcompimag(GEN x) { return sqcompimag0(x,0); }
GEN
sqcompimagraw(GEN x) { return sqcompimag0(x,1); }

static GEN
qfr_unit_by_disc(GEN D, long prec)
{
  GEN y = cgetg(5,t_QFR), isqrtD;
  pari_sp av = avma;
  long r;

  check_quaddisc_real(D, /*junk*/&r, "qfr_unit_by_disc");
  gel(y,1) = gen_1; isqrtD = sqrti(D);
  if ((r & 1) != mod2(isqrtD)) /* we know isqrtD > 0 */
    isqrtD = gerepileuptoint(av, addsi(-1,isqrtD));
  gel(y,2) = isqrtD; av = avma;
  gel(y,3) = gerepileuptoint(av, shifti(subii(sqri(isqrtD), D),-2));
  gel(y,4) = real_0(prec); return y;
}
GEN
qfr_unit(GEN x)
{
  long prec;
  if (typ(x) != t_QFR) pari_err(typeer,"qfr_unit");
  prec = precision(gel(x,4));
  if (!prec) pari_err(talker,"not a t_REAL in 4th component of a t_QFR");
  return qfr_unit_by_disc(qf_disc(x), prec);
}

static GEN
qfi_unit_by_disc(GEN D)
{
  GEN y = cgetg(4,t_QFI);
  long r;

  check_quaddisc_imag(D, &r, "qfi_unit_by_disc");
  gel(y,1) = gen_1;
  gel(y,2) = r? gen_1: gen_0;
  /* upon return, y[3] = (1-D) / 4 or -D / 4, whichever is an integer */
  gel(y,3) = shifti(D,-2);
  if (r) {
    pari_sp av = avma;
    gel(y,3) = gerepileuptoint(av, addis(gel(y,3),-1));
  }
  /* at this point y[3] < 0 */
  setsigne(y[3], 1); return y;
}
GEN
qfi_unit(GEN x)
{
  if (typ(x) != t_QFI) pari_err(typeer,"qfi_unit");
  return qfi_unit_by_disc(qf_disc(x));
}

static GEN
invraw(GEN x)
{
  GEN y = gcopy(x);
  setsigne(y[2], -signe(y[2]));
  if (typ(y) == t_QFR) setsigne(y[4], -signe(y[4]));
  return y;
}
GEN
powrealraw(GEN x, long n)
{
  pari_sp av = avma;
  long m;
  GEN y;

  if (typ(x) != t_QFR) pari_err(talker,"not a t_QFR in powrealraw");
  if (!n) return qfr_unit(x);
  if (n== 1) return gcopy(x);
  if (n==-1) return invraw(x);

  y = NULL; m = labs(n);
  for (; m>1; m>>=1)
  {
    if (m&1) y = y? comprealraw(y,x): x;
    x = sqcomprealraw(x);
  }
  y = y? comprealraw(y,x): x;
  if (n < 0) y = invraw(y);
  return gerepileupto(av,y);
}
GEN
powimagraw(GEN x, long n)
{
  pari_sp av = avma;
  long m;
  GEN y;

  if (typ(x) != t_QFI) pari_err(talker,"not a t_QFI in powimag");
  if (!n) return qfi_unit(x);
  if (n== 1) return gcopy(x);
  if (n==-1) return invraw(x);

  y = NULL; m = labs(n);
  for (; m>1; m>>=1)
  {
    if (m&1) y = y? compimagraw(y,x): x;
    x = sqcompimagraw(x);
  }
  y = y? compimagraw(y,x): x;
  if (n < 0) y = invraw(y);
  return gerepileupto(av,y);
}

GEN
powraw(GEN x, long n)
{ return (typ(x)==t_QFI)? powimagraw(x,n): powrealraw(x,n); }

static long
parteucl(GEN L, GEN *d, GEN *v3, GEN *v, GEN *v2)
{
  long z;
  *v = gen_0; *v2 = gen_1;
  for (z=0; absi_cmp(*v3,L) > 0; z++)
  {
    GEN t3, t2 = subii(*v, mulii(truedvmdii(*d,*v3,&t3),*v2));
    *v = *v2; *d = *v3; *v2 = t2; *v3 = t3;
  }
  return z;
}

/* composition: Shanks' NUCOMP & NUDUPL */
/* L = floor((|d|/4)^(1/4)) */
GEN
nucomp(GEN x, GEN y, GEN L)
{
  pari_sp av = avma;
  long z;
  GEN a, a1, a2, b2, b, d, d1, g, n, p1, q1, q2, s, u, u1, v, v1, v2, v3, Q;

  if (x==y) return nudupl(x,L);
  if (typ(x) != t_QFI || typ(y) != t_QFI) pari_err(talker,"not a t_QFI in nucomp");

  if (absi_cmp(gel(x,1),gel(y,1)) < 0) swap(x, y);
  s = shifti(addii(gel(x,2),gel(y,2)), -1);
  n = subii(gel(y,2), s);
  a1 = gel(x,1);
  a2 = gel(y,1); d = bezout(a2,a1,&u,&v);
  if (is_pm1(d)) { a = negi(mulii(u,n)); d1 = d; }
  else
    if (remii(s,d) == gen_0) /* d | s */
    {
      a = negi(mulii(u,n)); d1 = d;
      a1 = diviiexact(a1, d1);
      a2 = diviiexact(a2, d1);
      s = diviiexact(s, d1);
    }
    else
    {
      GEN p2, l;
      d1 = bezout(s,d,&u1,&v1);
      if (!is_pm1(d1))
      {
        a1 = diviiexact(a1,d1);
        a2 = diviiexact(a2,d1);
        s = diviiexact(s,d1);
        d = diviiexact(d,d1);
      }
      p1 = remii(gel(x,3),d);
      p2 = remii(gel(y,3),d);
      l = modii(mulii(negi(u1), addii(mulii(u,p1),mulii(v,p2))), d);
      a = subii(mulii(l,diviiexact(a1,d)), mulii(u,diviiexact(n,d)));
    }
  a = modii(a,a1); p1 = subii(a,a1); if (absi_cmp(a,p1) > 0) a = p1;
  d = a1; v3 = a; z = parteucl(L, &d,&v3, &v,&v2);
  Q = cgetg(4,t_QFI);
  if (!z) {
    g = diviiexact(addii(mulii(v3,s),gel(y,3)), d);
    b = a2;
    b2 = gel(y,2);
    v2 = d1;
    gel(Q,1) = mulii(d,b);
  } else {
    GEN e, q3, q4;
    if (z&1) { v3 = negi(v3); v2 = negi(v2); }
    b = diviiexact(addii(mulii(a2,d), mulii(n,v)), a1);
    e = diviiexact(addii(mulii(s,d),mulii(gel(y,3),v)), a1);
    q3 = mulii(e,v2);
    q4 = subii(q3,s);
    b2 = addii(q3,q4);
    g = diviiexact(q4,v);
    if (!is_pm1(d1)) { v2 = mulii(d1,v2); v = mulii(d1,v); b2 = mulii(d1,b2); }
    gel(Q,1) = addii(mulii(d,b), mulii(e,v));
  }
  q1 = mulii(b, v3);
  q2 = addii(q1,n);
  gel(Q,2) = addii(b2, z? addii(q1,q2): shifti(q1, 1));
  gel(Q,3) = addii(mulii(v3,diviiexact(q2,d)), mulii(g,v2));
  return gerepileupto(av, redimag(Q));
}

GEN
nudupl(GEN x, GEN L)
{
  pari_sp av = avma;
  long z;
  GEN u, v, d, d1, p1, a, b, c, a2, b2, c2, Q, v2, v3, g;

  if (typ(x) != t_QFI) pari_err(talker,"not a t_QFI in nudupl");
  a = gel(x,1);
  b = gel(x,2);
  d1 = bezout(b,a, &u,&v);
  if (!is_pm1(d1))
  {
    a = diviiexact(a, d1);
    b = diviiexact(b, d1);
  }
  c = modii(negi(mulii(u,gel(x,3))), a);
  p1 = subii(c,a); if (absi_cmp(c,p1) > 0) c = p1;
  d = a; v3 = c; z = parteucl(L, &d,&v3, &v,&v2);
  a2 = sqri(d);
  c2 = sqri(v3);
  Q = cgetg(4,t_QFI);
  if (!z) {
    g = diviiexact(addii(mulii(v3,b),gel(x,3)), d);
    b2 = gel(x,2);
    v2 = d1;
    gel(Q,1) = a2;
  } else {
    GEN e;
    if (z&1) { v = negi(v); d = negi(d); }
    e = diviiexact(addii(mulii(gel(x,3),v), mulii(b,d)), a);
    g = diviiexact(subii(mulii(e,v2), b), v);
    b2 = addii(mulii(e,v2), mulii(v,g));
    if (!is_pm1(d1)) { b2 = mulii(d1,b2); v = mulii(d1,v); v2 = mulii(d1,v2); }
    gel(Q,1) = addii(a2, mulii(e,v));
  }
  gel(Q,2) = addii(b2, subii(sqri(addii(d,v3)), addii(a2,c2)));
  gel(Q,3) = addii(c2, mulii(g,v2));
  return gerepileupto(av, redimag(Q));
}

static GEN
mul_nucomp(void *l, GEN x, GEN y) { return nucomp(x, y, (GEN)l); }
static GEN
mul_nudupl(void *l, GEN x) { return nudupl(x, (GEN)l); }
GEN
nupow(GEN x, GEN n)
{
  pari_sp av;
  GEN y, l;

  if (typ(n) != t_INT) pari_err(talker,"not an integer exponent in nupow");
  if (gcmp1(n)) return gcopy(x);
  av = avma; y = qfi_unit(x);
  if (!signe(n)) return y;

  l = sqrti(shifti(sqrti(gel(y,3)),1));
  y = leftright_pow(x, n, (void*)l, &mul_nudupl, &mul_nucomp);
  if (signe(n) < 0
  && !absi_equal(gel(y,1),gel(y,2))
  && !absi_equal(gel(y,1),gel(y,3))) setsigne(y[2],-signe(y[2]));
  return gerepileupto(av, y);
}

/* Reduction */

/* reduce b mod 2*a. Update b,c */
#define REDB_i(a,b,c)\
  GEN r, a2 = shifti(a, 1), q = dvmdii(b, a2, &r);\
  if (signe(b) >= 0) {\
    if (absi_cmp(r, a) > 0) { q = addis(q,  1); r = subii(r, a2); }\
  } else { /* r <= 0 */\
    if (absi_cmp(r, a) >= 0){ q = addis(q, -1); r = addii(r, a2); }\
  }\
  c = subii(c, mulii(q, shifti(addii(b, r),-1)));\
  b = r;

#define REDBU(a,b,c, u1,u2) { REDB_i(a,b,c); u2 = subii(u2, mulii(q, u1)); }
#define REDB(a,b,c) { REDB_i(a,b,c); }

GEN
redimagsl2(GEN q, GEN *U)
{
  GEN Q = cgetg(4, t_QFI);
  pari_sp av = avma, lim = stack_lim(av, 1);
  GEN z, u1,u2,v1,v2, a = gel(q,1), b = gel(q,2), c = gel(q,3);
  long cmp;
  /* upper bound for size of final (a,b,c) */
  (void)new_chunk(lgefint(a) + lgefint(b) + lgefint(c) + 3);
  u1 = gen_1;
  u2 = gen_0; cmp = absi_cmp(a, b);
  if (cmp <= 0 && (cmp || signe(b) < 0)) REDBU(a,b,c, u1,u2);
  for(;;)
  {
    cmp = absi_cmp(a, c); if (cmp <= 0) break;
    swap(a,c); b = negi(b);
    z = u1; u1 = u2; u2 = negi(z);
    REDBU(a,b,c, u1,u2);
    if (low_stack(lim, stack_lim(av, 1))) {
      if (DEBUGMEM>1) pari_warn(warnmem, "redimagsl2");
      gerepileall(av, 5, &a,&b,&c, &u1,&u2);
    }
  }
  if (cmp == 0 && signe(b) < 0)
  {
    b = negi(b);
    z = u1; u1 = u2; u2 = negi(z);
  }
  avma = av;
  a = icopy(a); gel(Q,1) = a;
  b = icopy(b); gel(Q,2) = b;
  c = icopy(c); gel(Q,3) = c;
  u1 = icopy(u1);
  u2 = icopy(u2); av = avma;

  /* Let q = (A,B,C). q o [u1,u2; v1,v2] = Q implies
   * [v1,v2] = (1/C) [(b-B)/2 u1 - a u2, c u1 - (b+B)/2 u2] */
  z = shifti(subii(b, gel(q,2)), -1);
  v1 = subii(mulii(z, u1), mulii(a, u2)); v1 = diviiexact(v1, gel(q,3));
  z = subii(z, b);
  v2 = addii(mulii(z, u2), mulii(c, u1)); v2 = diviiexact(v2, gel(q,3));
  avma = av;
  v1 = icopy(v1);
  v2 = icopy(v2);
  *U = cgetg(3, t_MAT);
  gel(*U,1) = mkcol2(u1,v1);
  gel(*U,2) = mkcol2(u2,v2); return Q;
}

GEN
redimag(GEN q)
{
  GEN Q = cgetg(4, t_QFI);
  pari_sp av = avma, lim = stack_lim(av, 1);
  GEN a = gel(q,1), b = gel(q,2), c = gel(q,3);
  long cmp;
  /* upper bound for size of final (a,b,c) */
  (void)new_chunk(lgefint(a) + lgefint(b) + lgefint(c) + 3);
  cmp = absi_cmp(a, b);
  if (cmp <= 0 && (cmp || signe(b) < 0)) REDB(a,b,c);
  for(;;)
  {
    cmp = absi_cmp(a, c); if (cmp <= 0) break;
    swap(a,c); b = negi(b); /* apply rho */
    REDB(a,b,c);
    if (low_stack(lim, stack_lim(av, 1))) {
      if (DEBUGMEM>1) pari_warn(warnmem, "redimag");
      gerepileall(av, 3, &a,&b,&c);
    }
  }
  if (cmp == 0 && signe(b) < 0) b = negi(b);

  avma = av;
  gel(Q,1) = icopy(a);
  gel(Q,2) = icopy(b);
  gel(Q,3) = icopy(c); return Q;
}

static GEN
rhoimag(GEN x)
{
  GEN a = gel(x,1), b = gel(x,2), c = gel(x,3);
  int fl = absi_cmp(a, c);
  if (fl <= 0) {
    int fg = absi_cmp(a, b);
    if (fg >= 0) {
      x = qfi(a,b,c);
      if ((!fl || !fg) && signe(x[2]) < 0) setsigne(x[2], 1);
      return x;
    }
  }
  x = cgetg(4, t_QFI);
  (void)new_chunk(lgefint(a) + lgefint(b) + lgefint(c) + 3);
  swap(a,c); b = negi(b);
  REDB(a, b, c); avma = (pari_sp)x;
  gel(x,1) = icopy(a);
  gel(x,2) = icopy(b);
  gel(x,3) = icopy(c); return x;
}

/* qfr3 / qfr5 */

/* t_QFR are unusable: D, sqrtD, isqrtD are recomputed all the time and the
 * logarithmic Shanks's distance is costly and hard to control.
 * qfr3 / qfr5 routines take a container of t_INTs (e.g a t_VEC) as argument,
 * at least 3 (resp. 5) components [it is a feature that they do not check the
 * precise type or length of the input]. They return a vector of length 3
 * (resp. 5). A qfr3 [a,b,c] contains the form coeffs, in a qfr5 [a,b,c, e,d]
 * the t_INT e is a binary exponent, d a t_REAL, coding the distance in
 * multiplicative form: the true distance is obtained from qfr5_dist.
 * D, sqrtD, isqrtD are included in the function's arguments [sqrtD is only
 * used for distance computations].
 * All other qfr routines are obsolete (inefficient) wrappers */

/* static functions are not stack-clean. Unless mentionned otherwise, public
 * functions are. */

#define EMAX 22
static void
fix_expo(GEN x)
{
  long d = expo(x[5]) - (1 << EMAX);
  if (d >= 0) {
    gel(x,4) = addsi(1, gel(x,4));
    setexpo(x[5], d);
  }
}

/* (1/2) log (d * 2^{e * 2^EMAX}). Not stack clean if e != 0 */
GEN
qfr5_dist(GEN e, GEN d, long prec)
{
  GEN t = logr_abs(d);
  if (signe(e)) {
    GEN u = mulir(e, mplog2(prec));
    setexpo(u, expo(u+EMAX)); t = addrr(t, u);
  }
  setexpo(t, expo(t)-1); return t;
}

static void
rho_get_BC(GEN *B, GEN *C, GEN b, GEN c, GEN D, GEN isqrtD)
{
  GEN t, u;
  u = shifti(c,1); if (u == gen_0) pari_err(talker, "reducible form in qfr_rho");
  t = (absi_cmp(isqrtD,c) >= 0)? isqrtD: c;
  u = remii(addii_sign(t,1, b,signe(b)), u);
  *B = addii_sign(t, 1, u, -signe(u)); /* |t| - (|t|+b) % |2c| */
  if (*B == gen_0)
  { u = shifti(D, -2); setsigne(u, -1); }
  else
    u = shifti(addii_sign(sqri(*B),1, D,-1), -2);
  *C = diviiexact(u, c); /* = (B^2-D)/4c */
}
/* Not stack-clean */
GEN
qfr3_rho(GEN x, GEN D, GEN isqrtD)
{
  GEN B, C, y, b = gel(x,2), c = gel(x,3);

  rho_get_BC(&B, &C, b, c, D, isqrtD);
  y = cgetg(4, t_VEC);
  gel(y,1) = c;
  gel(y,2) = B;
  gel(y,3) = C; return y;
}
/* Not stack-clean */
GEN
qfr5_rho(GEN x, GEN D, GEN sqrtD, GEN isqrtD)
{
  GEN B, C, y, b = gel(x,2), c = gel(x,3);
  long sb = signe(b);

  rho_get_BC(&B, &C, b, c, D, isqrtD);
  y = cgetg(6, t_VEC);
  gel(y,1) = c;
  gel(y,2) = B;
  gel(y,3) = C;
  y[4] = x[4];
  y[5] = x[5];
  if (sb) {
    GEN t = subii(sqri(b), D);
    if (sb < 0)
      t = divir(t, gsqr(subir(b,sqrtD)));
    else
      t = divri(gsqr(addir(b,sqrtD)), t);
    /* t = (b + sqrt(D)) / (b - sqrt(D)), evaluated stably */
    gel(y,5) = mulrr(t, gel(y,5)); fix_expo(y);
  }
  return y;
}

/* Not stack-clean */
GEN
qfr_to_qfr5(GEN x, long prec)
{
  GEN y = cgetg(6,t_VEC);
  y[1] = x[1];
  y[2] = x[2];
  y[3] = x[3];
  gel(y,4) = gen_0;
  gel(y,5) = real_1(prec); return y;
}

/* d0 = initial distance, x = [a,b,c, expo(d), d], d = exp(2*distance) */
static GEN
qfr5_to_qfr(GEN x, GEN d0)
{
  GEN y;
  if (lg(x) ==  6)
  {
    GEN n = gel(x,4), d = absr(gel(x,5));
    if (signe(n))
    {
      n = addis(shifti(n, EMAX), expo(d));
      setexpo(d, 0); d = logr_abs(d);
      d = mpadd(d, mulir(n, mplog2(lg(d0))));
    }
    else
      d = gcmp1(d)? NULL: logr_abs(d); /* avoid loss of precision */
    if (d) d0 = addrr(d0, shiftr(d,-1));
  }
  y = cgetg(5, t_QFR);
  y[1] = x[1];
  y[2] = x[2];
  y[3] = x[3];
  gel(y,4) = d0; return y;
}

/* Not stack-clean */
GEN
qfr3_to_qfr(GEN x, GEN d)
{
  GEN z = cgetg(5, t_QFR);
  z[1] = x[1];
  z[2] = x[2];
  z[3] = x[3];
  gel(z,4) = d; return z;
}

static int
abi_isreduced(GEN a, GEN b, GEN isqrtD)
{
  if (signe(b) > 0 && absi_cmp(b, isqrtD) <= 0)
  {
    GEN t = addii_sign(isqrtD,1, shifti(a,1),-1);
    long l = absi_cmp(b, t); /* compare |b| and |floor(sqrt(D)) - |2a|| */
    if (l > 0 || (l == 0 && signe(t) < 0)) return 1;
  }
  return 0;
}

INLINE int
qfr_isreduced(GEN x, GEN isqrtD)
{
  return abi_isreduced(gel(x,1),gel(x,2),isqrtD);
}

/* Not stack-clean */
GEN
qfr5_red(GEN x, GEN D, GEN sqrtD, GEN isqrtD) {
  while (!qfr_isreduced(x,isqrtD)) x = qfr5_rho(x,D,sqrtD,isqrtD);
  return x;
}
/* Not stack-clean */
GEN
qfr3_red(GEN x, GEN D, GEN isqrtD) {
  while (!qfr_isreduced(x,isqrtD)) x = qfr3_rho(x,D,isqrtD);
  return x;
}

static void
get_disc(GEN x, GEN *D)
{
  if (!*D) *D = qf_disc(x);
  else if (typ(*D) != t_INT) pari_err(arither1);
  if (!signe(*D)) pari_err(talker,"reducible form in qfr_init");
}

static GEN
qfr5_init(GEN x, GEN *D, GEN *isqrtD, GEN *sqrtD)
{
  GEN d = gel(x,4);
  long prec = lg(d), l = nbits2prec(-expo(d));
  if (l > prec) prec = l;
  if (prec < 3) prec = 3;
  x = qfr_to_qfr5(x,prec);

  get_disc(x, D);
  if (!*sqrtD) *sqrtD = sqrtr(itor(*D,prec));
  else if (typ(*sqrtD) != t_REAL) pari_err(arither1);

  if (!*isqrtD) *isqrtD = truncr(*sqrtD);
  else if (typ(*isqrtD) != t_INT) pari_err(arither1);
  return x;
}
static GEN
qfr3_init(GEN x, GEN *D, GEN *isqrtD)
{
  get_disc(x, D);

  if (!*isqrtD) *isqrtD = sqrti(*D);
  else if (typ(*isqrtD) != t_INT) pari_err(arither1);
  return x;
}

#define qf_NOD  2
#define qf_STEP 1

static GEN
redreal0(GEN x, long flag, GEN D, GEN isqrtD, GEN sqrtD)
{
  pari_sp av = avma;
  GEN d = gel(x,4);
  if (typ(x) != t_QFR) pari_err(talker,"not a real quadratic form in redreal");
  x = (flag & qf_NOD)? qfr3_init(x, &D,&isqrtD)
                     : qfr5_init(x, &D,&isqrtD,&sqrtD);
  switch(flag) {
    case 0:              x = qfr5_red(x,D,sqrtD,isqrtD); break;
    case qf_NOD:         x = qfr3_red(x,D,isqrtD); break;
    case qf_STEP:        x = qfr5_rho(x,D,sqrtD,isqrtD); break;
    case qf_STEP|qf_NOD: x = qfr3_rho(x,D,isqrtD); break;
    default: pari_err(flagerr,"qfbred");
  }
  return gerepilecopy(av, qfr5_to_qfr(x,d));
}
GEN
redreal(GEN x)
{ return redreal0(x,0,NULL,NULL,NULL); }
GEN
rhoreal(GEN x)
{ return redreal0(x,qf_STEP,NULL,NULL,NULL); }
GEN
redrealnod(GEN x, GEN isqrtD)
{ return redreal0(x,qf_NOD,NULL,isqrtD,NULL); }
GEN
rhorealnod(GEN x, GEN isqrtD)
{ return redreal0(x,qf_STEP|qf_NOD,NULL,isqrtD,NULL); }
GEN
qfbred0(GEN x, long flag, GEN D, GEN isqrtD, GEN sqrtD)
{
  if (typ(x) == t_QFI)
    return (flag & qf_STEP)? rhoimag(x): redimag(x);
  return redreal0(x,flag,D,isqrtD,sqrtD);
}

GEN
qfr5_comp(GEN x, GEN y, GEN D, GEN sqrtD, GEN isqrtD)
{
  pari_sp av = avma;
  GEN z = cgetg(6,t_VEC); qfb_comp(z,x,y);
  if (x == y)
  {
    gel(z,4) = shifti(gel(x,4),1);
    gel(z,5) = gsqr(gel(x,5));
  }
  else
  {
    gel(z,4) = addii(gel(x,4),gel(y,4));
    gel(z,5) = mulrr(gel(x,5),gel(y,5));
  }
  fix_expo(z); z = qfr5_red(z,D,sqrtD,isqrtD);
  return gerepilecopy(av,z);
}
/* Not stack-clean */
GEN
qfr3_comp(GEN x, GEN y, GEN D, GEN isqrtD)
{
  GEN z = cgetg(4,t_VEC); qfb_comp(z,x,y);
  return qfr3_red(z, D, isqrtD);
}

/* assume n != 0, return x^|n|. Not stack-clean */
GEN
qfr5_pow(GEN x, GEN n, GEN D, GEN sqrtD, GEN isqrtD)
{
  GEN y = NULL;
  long i, m;
  for (i=lgefint(n)-1; i>1; i--)
  {
    m = n[i];
    for (; m; m>>=1)
    {
      if (m&1) y = y? qfr5_comp(y,x,D,sqrtD,isqrtD): x;
      if (m == 1 && i == 2) break;
      x = qfr5_comp(x,x,D,sqrtD,isqrtD);
    }
  }
  return y;
}
/* assume n != 0, return x^|n|. Not stack-clean */
GEN
qfr3_pow(GEN x, GEN n, GEN D, GEN isqrtD)
{
  GEN y = NULL;
  long i, m;
  for (i=lgefint(n)-1; i>1; i--)
  {
    m = n[i];
    for (; m; m>>=1)
    {
      if (m&1) y = y? qfr3_comp(y,x,D,isqrtD): x;
      if (m == 1 && i == 2) break;
      x = qfr3_comp(x,x,D,isqrtD);
    }
  }
  return y;
}

static GEN
qfr_inv(GEN x) {
  GEN z = cgetg(5, t_QFR);
  z[1] = x[1];
  gel(z,2) = negi(gel(x,2));
  z[3] = x[3];
  z[4] = x[4]; return z;
}
/* assume n != 0 */
GEN
qfr_pow(GEN x, GEN n)
{
  pari_sp av = avma;
  GEN D, sqrtD, isqrtD, d0;

  if (is_pm1(n)) return signe(n) > 0? gcopy(x): ginv(x);
  if (signe(n) < 0) x = qfr_inv(x);
  d0 = gel(x,4);
  D = sqrtD = isqrtD = NULL;
  if (!signe(d0)) {
    x = qfr3_init(x, &D,&isqrtD);
    x = qfr3_pow(x, n, D,isqrtD);
    x = qfr3_to_qfr(x, d0);
  } else {
    x = qfr5_init(x, &D,&isqrtD,&sqrtD);
    x = qfr5_pow(qfr_to_qfr5(x, lg(sqrtD)), n, D,sqrtD,isqrtD);
    x = qfr5_to_qfr(x, mulri(d0,n));
  }
  return gerepilecopy(av, x);
}

/* Prime forms associated to prime ideals of degree 1 */

/* assume x != 0 a t_INT, p > 0
 * Return a t_QFI, but discriminant sign is not checked: can be used for
 * real forms as well */
GEN
primeform_u(GEN x, ulong p)
{
  GEN c, y = cgetg(4, t_QFI);
  pari_sp av = avma;
  ulong b, s;

  s = mod8(x); if (signe(x) < 0 && s) s = 8-s;
  /* 2 or 3 mod 4 */
  if (s & 2) pari_err(talker,"discriminant not congruent to 0,1 mod 4 in primeform");
  if (p == 2) {
    switch(s) {
      case 0: b = 0; break;
      case 1: b = 1; break;
      case 4: b = 2; break;
      default: pari_err(sqrter5); b = 0; /* -Wall */
    }
    c = shifti(subsi(s,x), -3);
  } else {
    b = Fl_sqrt(umodiu(x,p), p); if (b == ~0UL) pari_err(sqrter5);
    /* mod(b) != mod2(x) ? */
    if ((b & 1) != (s & 1)) b = p - b;
    c = diviuexact(shifti(subii(sqru(b), x), -2), p);
  }
  gel(y,3) = gerepileuptoint(av, c);
  gel(y,2) = utoi(b);
  gel(y,1) = utoipos(p); return y;
}

/* special case: p = 1 return unit form */
GEN
primeform(GEN x, GEN p, long prec)
{
  pari_sp av;
  long s, sx = signe(x), sp = signe(p);
  GEN y, b, absp;

  if (typ(x) != t_INT || !sx) pari_err(arither1);
  if (typ(p) != t_INT || !sp) pari_err(arither1);
  if (is_pm1(p)) {
    if (sx < 0) return qfi_unit_by_disc(x);
    y = qfr_unit_by_disc(x,prec);
    if (sp < 0) { gel(y,1) = negi(gel(y,1)); gel(y,3) = negi(gel(y,3)); }
    return y;
  }
  if (sp < 0 && sx < 0) pari_err(impl,"negative definite t_QFI");
  if (lgefint(p) == 3)
  { 
    y = primeform_u(x, p[2]);
    if (sx < 0) return y;
    if (sp < 0) { gel(y,1) = negi(gel(y,1)); gel(y,3) = negi(gel(y,3)); }
    return gcopy( qfr3_to_qfr(y, real_0(prec)) );
  }
  s = mod8(x);
  if (sx < 0)
  {
    if (s) s = 8-s;
    y = cgetg(4, t_QFI);
  }
  else
  {
    y = cgetg(5, t_QFR);
    gel(y,4) = real_0(prec);
  }
  /* 2 or 3 mod 4 */
  if (s & 2) pari_err(talker,"discriminant not congruent to 0,1 mod 4 in primeform");
  absp = absi(p); av = avma;
  b = Fp_sqrt(x, absp); if (!b) pari_err(sqrter5);
  s &= 1; /* s = x mod 2 */
  /* mod(b) != mod2(x) ? [Warning: we may have b == 0] */
  if ((!signe(b) && s) || mod2(b) != s) b = gerepileuptoint(av, subii(absp,b));

  av = avma;
  gel(y,3) = gerepileuptoint(av, diviiexact(shifti(subii(sqri(b), x), -2), p));
  gel(y,2) = b;
  gel(y,1) = gcopy(p);
  return y;
}

/* Let M and N in SL_2(Z), return (N*M^-1)[,1]
 * 
*/
static GEN 
SL2_div_mul_e1(GEN N, GEN M)
{
  GEN b = gcoeff(M, 2, 1);
  GEN d = gcoeff(M, 2, 2);
  GEN p2 = cgetg(3, t_VEC);
  gel(p2,1) = subii(mulii(gcoeff(N, 1, 1), d),
                    mulii(gcoeff(N, 1, 2), b));
  gel(p2,2) = subii(mulii(gcoeff(N, 2, 1), d),
                    mulii(gcoeff(N, 2, 2), b));
  return p2;
}
/* Let M and N in SL_2(Z), return (N*[1,0;0,-1]*M^-1)[,1]
 * 
*/
static GEN 
SL2_swap_div_mul_e1(GEN N, GEN M)
{
  GEN b = gcoeff(M, 2, 1);
  GEN d = gcoeff(M, 2, 2);
  GEN p2 = cgetg(3, t_VEC);
  gel(p2,1) = addii(mulii(gcoeff(N, 1, 1), d),
                    mulii(gcoeff(N, 1, 2), b));
  gel(p2,2) = addii(mulii(gcoeff(N, 2, 1), d),
                    mulii(gcoeff(N, 2, 2), b));
  return p2;
}

/*Test equality modulo GL2*/
static int
GL2_qfb_equal(GEN a, GEN b)
{
  return equalii(gel(a,1),gel(b,1))
   && absi_equal(gel(a,2),gel(b,2))
   &&    equalii(gel(a,3),gel(b,3));
}

static GEN
qfbsolve_cornacchia(GEN c, GEN p, int swap)
{
  pari_sp av = avma;
  GEN M, N;
  if (kronecker(negi(c), p) < 0 || !cornacchia(c, p, &M,&N)) {
    avma = av; return gen_0;
  }
  return gerepilecopy(av, swap? mkvec2(N,M): mkvec2(M,N));
}

GEN
qfbimagsolvep(GEN Q, GEN p)
{
  GEN M, N, x,y, a,b,c, d;
  pari_sp av = avma;
  if (!signe(gel(Q,2)))
  {
    a = gel(Q,1);
    c = gel(Q,3); /* if principal form, use faster cornacchia */
    if (gcmp1(a)) return qfbsolve_cornacchia(c, p, 0);
    if (gcmp1(c)) return qfbsolve_cornacchia(a, p, 1);
  }
  d = qf_disc(Q); if (kronecker(d,p) < 0) return gen_0;
  a = redimagsl2(Q, &N);
  if (is_pm1(gel(a,1))) /* principal form */
  {
    long r;
    if (!signe(gel(a,2)))
    {
      a = qfbsolve_cornacchia(gel(a,3), p, 0);
      if (a == gen_0) { avma = av; return gen_0; }
      return gerepileupto(av, gmul(a, shallowtrans(N)));
    }
    /* x^2 + xy + ((1-d)/4)y^2 = p <==> (2x + y)^2 - d y^2 = 4p */
    if (!cornacchia2(negi(d), p, &x, &y)) { avma = av; return gen_0; }
    x = divis_rem(subii(x,y), 2, &r); if (r) { avma = av; return gen_0; }
    return gerepileupto(av, gmul(mkvec2(x,y), shallowtrans(N)));
  }
  b = redimagsl2(primeform(d, p, 0), &M);
  if (!GL2_qfb_equal(a,b)) { avma = av; return gen_0; }
  if (signe(gel(a,2))==signe(gel(b,2)))
    x = SL2_div_mul_e1(N,M);
  else
    x = SL2_swap_div_mul_e1(N,M);
  return gerepilecopy(av, x);
}

GEN
redrealsl2step(GEN A)
{
  pari_sp ltop = avma;
  GEN N;
  GEN V = gel(A,1); 
  GEN M = gel(A,2);
  GEN a = gel(V,1); 
  GEN b = gel(V,2);
  GEN c = gel(V,3);
  GEN d = qf_disc0(a,b,c);
  GEN rd = sqrti(d); 
  GEN ac = mpabs(c);
  GEN r = addii(b, gmax(rd, ac));
  GEN q = truedvmdii(r, shifti(ac, 1), NULL);
  r = subii(mulii(shifti(q, 1), ac), b);
  a = c; b = r;
  c = truedvmdii(subii(sqri(r), d), shifti(c,2), NULL);
  if (signe(a) < 0) q = negi(q);
  N = mkmat2(gel(M,2),
             mkcol2(subii(mulii(q, gcoeff(M, 1, 2)), gcoeff(M, 1, 1)),
                    subii(mulii(q, gcoeff(M, 2, 2)), gcoeff(M, 2, 1))));
  return gerepilecopy(ltop, mkvec2(mkvec3(a,b,c),N));
}

GEN
redrealsl2(GEN V)
{
  pari_sp ltop = avma, btop, st_lim;
  GEN u1, u2, v1, v2;
  GEN M;
  GEN a = gel(V,1);
  GEN b = gel(V,2);
  GEN c = gel(V,3);
  GEN d = qf_disc0(a,b,c);
  GEN rd = sqrti(d);
  btop = avma; st_lim = stack_lim(btop, 1);
  u1 = v2 = gen_1; v1 = u2 = gen_0; 
  while (!abi_isreduced(a,b,rd))
  {
    GEN ac = mpabs(c);
    GEN r = addii(b, gmax(rd,ac));
    GEN q = truedvmdii(r, mulsi(2, ac), NULL);
    r = subii(mulii(mulis(q, 2), ac), b);
    a = c; b = r;
    c = truedvmdii(subii(sqri(r), d), mulsi(4, c), NULL);
    q = mulis(q, signe(a));
    r = u1; u1 = v1; v1 = subii(mulii(q, v1), r);
    r = u2; u2 = v2; v2 = subii(mulii(q, v2), r);
    if (low_stack(st_lim, stack_lim(btop, 1)))
    {
      GEN *bptr[7];
      bptr[0]=&a; bptr[1]=&b; bptr[2]=&c;
      bptr[3]=&u1; bptr[4]=&u2;
      bptr[5]=&v1; bptr[6]=&v2;
      gerepilemany(ltop, bptr, 7);
    }
  }
  M = mkmat2(mkcol2(u1,u2), mkcol2(v1,v2));
  return gerepilecopy(ltop, mkvec2(mkvec3(a,b,c), M));
}

GEN
qfbrealsolvep(GEN Q, GEN p)
{
  pari_sp ltop = avma, btop, st_lim;
  GEN N, P, P1, P2, M, d = qf_disc(Q);
  if (kronecker(d, p) < 0) { avma = ltop; return gen_0; }
  M = N = redrealsl2(Q);
  P = primeform(d, p, DEFAULTPREC);
  P1 = redrealsl2(P);
  gel(P,2) = negi(gel(P,2));
  P2 = redrealsl2(P);
  btop = avma; st_lim = stack_lim(btop, 1);
  while (!gequal(gel(M,1), gel(P1,1)) && !gequal(gel(M,1), gel(P2,1)))
  {
    M = redrealsl2step(M);
    if (gequal(gel(M,1), gel(N,1))) { avma = ltop; return gen_0; }
    if (low_stack(st_lim, stack_lim(btop, 1))) M = gerepileupto(btop, M);
  }
  if (gequal(gel(M,1),gel(P1,1)))
    return gerepilecopy(ltop, SL2_div_mul_e1(gel(M,2),gel(P1,2)));
  else
    return gerepilecopy(ltop, SL2_div_mul_e1(gel(M,2),gel(P2,2)));
}

GEN
qfbsolve(GEN Q,GEN n)
{
  if (typ(n)!=t_INT) pari_err(typeer,"qfbsolve");
  switch(typ(Q))
  {
  case t_QFI: return qfbimagsolvep(Q,n);
  case t_QFR: return qfbrealsolvep(Q,n);
  default:
    pari_err(typeer,"qfbsolve");
    return NULL; /* NOT REACHED */
  }
}

/* 1 if there exists x,y such that x^2 + dy^2 = p [prime], 0 otherwise */
long
cornacchia(GEN d, GEN p, GEN *px, GEN *py)
{
  pari_sp av = avma, av2, lim;
  GEN a, b, c, L, r;

  if (typ(d) != t_INT || typ(p) != t_INT) pari_err(typeer, "cornacchia");
  if (signe(d) <= 0) pari_err(talker, "d must be positive");
  *px = *py = gen_0;
  b = subii(p, d);
  if (signe(b) < 0) return 0;
  if (signe(b) == 0) { avma = av; *py = gen_1; return 1; }
  b = Fp_sqrt(b, p); /* sqrt(-d) */
  if (!b) { avma = av; return 0; }
  if (absi_cmp(shifti(b,1), p) > 0) b = subii(b,p);
  a = p; L = sqrti(p);
  av2 = avma; lim = stack_lim(av2, 1);
  while (absi_cmp(b, L) > 0)
  {
    r = remii(a, b); a = b; b = r;
    if (low_stack(lim, stack_lim(av2, 1))) {
      if (DEBUGMEM>1) pari_warn(warnmem,"cornacchia");
      gerepileall(av2, 2, &a,&b);
    }
  }
  a = subii(p, sqri(b));
  c = dvmdii(a, d, &r);
  if (r != gen_0 || !Z_issquarerem(c, &c)) { avma = av; return 0; }
  avma = av;
  *px = icopy(b);
  *py = icopy(c); return 1;
}
/* 1 if there exists x,y such that x^2 + dy^2 = 4p [p prime], 0 otherwise */
long
cornacchia2(GEN d, GEN p, GEN *px, GEN *py)
{
  pari_sp av = avma, av2, lim;
  GEN a, b, c, L, r, px4;
  long k;

  if (typ(d) != t_INT || typ(p) != t_INT) pari_err(typeer, "cornacchia");
  if (signe(d) <= 0) pari_err(talker, "d must be positive");
  *px = *py = gen_0;
  k = mod4(d);
  if (k == 1 || k == 2) pari_err(talker,"d must be 0 or 3 mod 4");
  px4 = shifti(p,2);
  if (absi_cmp(px4, d) < 0) { avma = av; return 0; }
  if (equaliu(p, 2))
  {
    avma = av;
    switch (itou_or_0(d)) {
      case 4: *px = gen_2; break;
      case 7: *px = gen_1; break;
      default: return 0;
    }
    *py = gen_1; return 1;
  } 
  b = Fp_sqrt(negi(d), p);
  if (!b) { avma = av; return 0; }
  if (!signe(b)) { /* d = p,2p,3p,4p */
    avma = av;
    if (absi_equal(d, px4)){ *py = gen_1; return 1; }
    if (absi_equal(d, p))  { *py = gen_2; return 1; }
    return 0;
  }
  if (mod2(b) != (k & 1)) b = subii(p,b);
  a = shifti(p,1); L = sqrti(px4);
  av2 = avma; lim = stack_lim(av2, 1);
  while (cmpii(b, L) > 0)
  {
    r = remii(a, b); a = b; b = r;
    if (low_stack(lim, stack_lim(av2, 1))) {
      if (DEBUGMEM>1) pari_warn(warnmem,"cornacchia");
      gerepileall(av2, 2, &a,&b);
    }
  }
  a = subii(px4, sqri(b));
  c = dvmdii(a, d, &r);
  if (r != gen_0 || !Z_issquarerem(c, &c)) { avma = av; return 0; }
  avma = av;
  *px = icopy(b);
  *py = icopy(c); return 1;
}
