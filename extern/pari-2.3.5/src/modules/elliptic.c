/* $Id: elliptic.c 12105 2010-02-03 15:04:07Z bill $

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
/**                       ELLIPTIC CURVES                          **/
/**                                                                **/
/********************************************************************/
#include "pari.h"
#include "paripriv.h"

#define is_inf(z) (lg(z) < 3)

void
checkpt(GEN z)
{ if (typ(z)!=t_VEC) pari_err(elliper1); }
void
checksell(GEN e)
{ if (typ(e)!=t_VEC || lg(e) < 6) pari_err(elliper1); }
void
checkell(GEN e)
{ if (typ(e)!=t_VEC || lg(e) < 14) pari_err(elliper1); }
void /* check for roots, we don't want a curve over Fp */
checkbell(GEN e)
{ if (typ(e)!=t_VEC || lg(e) < 20 || typ(e[14]) != t_COL) pari_err(elliper1); }

static void
checkch(GEN z)
{ if (typ(z)!=t_VEC || lg(z) != 5) pari_err(elliper1); }

/* 4 X^3 + b2 X^2 + 2b4 X + b6 */
static GEN
RHSpol(GEN e)
{
  GEN z = cgetg(6, t_POL); z[1] = evalsigne(1);
  z[2] = e[8];
  gel(z,3) = gmul2n(gel(e,7),1);
  z[4] = e[6];
  gel(z,5) = utoipos(4); return z;
}

/* x^3 + a2 x^2 + a4 x + a6 */
static GEN
ellRHS(GEN e, GEN x)
{
  GEN z;
  z = gadd(gel(e,2),x);
  z = gadd(gel(e,4), gmul(x,z));
  z = gadd(gel(e,5), gmul(x,z));
  return z;
}

/* a1 x + a3 */
static GEN
ellLHS0(GEN e, GEN x)
{
  return gcmp0(gel(e,1))? gel(e,3): gadd(gel(e,3), gmul(x,gel(e,1)));
}

static GEN
ellLHS0_i(GEN e, GEN x)
{
  return signe(e[1])? addii(gel(e,3), mulii(x, gel(e,1))): gel(e,3);
}

/* y^2 + a1 xy + a3 y */
static GEN
ellLHS(GEN e, GEN z)
{
  GEN y = gel(z,2);
  return gmul(y, gadd(y, ellLHS0(e,gel(z,1))));
}

/* 2y + a1 x + a3 */
static GEN
d_ellLHS(GEN e, GEN z)
{
  return gadd(ellLHS0(e, gel(z,1)), gmul2n(gel(z,2),1));
}

static GEN
ell_to_small(GEN E)
{
  long i;
  GEN e;
  if (lg(E) <= 14) return E;
  e = cgetg(14,t_VEC);
  for (i = 1; i < 14; i++) e[i] = E[i];
  return e;
}

static void
smallinitell0(GEN x, GEN y)
{
  GEN b2, b4, b6, b8, D, j, a11, a13, a33, b22, c4, c6;

  checksell(x);
  y[1] = x[1];
  y[2] = x[2];
  y[3] = x[3];
  y[4] = x[4];
  y[5] = x[5];
  a11 = gsqr(gel(y,1));
  b2 = gadd(a11, gmul2n(gel(y,2),2));
  gel(y,6) = b2; /* a1^2 + 4a2 */

  a13 = gmul(gel(y,1),gel(y,3));
  b4 = gadd(a13, gmul2n(gel(y,4),1));
  gel(y,7) = b4; /* a1 a3 + 2a4 */

  a33 = gsqr(gel(y,3));
  b6 = gadd(a33, gmul2n(gel(y,5),2));
  gel(y,8) = b6; /* a3^2 + 4 a6 */
  b8 = gsub(gadd(gmul(a11,gel(y,5)), gmul(b6, gel(y,2))),
            gmul(gel(y,4), gadd(gel(y,4),a13)));
  gel(y,9) = b8; /* a1^2 a6 + 4a6 a2 + a2 a3^2 + 4 a6 - a4(a4 + a1 a3) */

  b22 = gsqr(b2);
  c4 = gadd(b22, gmulsg(-24,b4));
  gel(y,10) = c4; /* b2^2 - 24 b4 */

  c6 = gadd(gmul(b2,gsub(gmulsg(36,b4),b22)),gmulsg(-216,b6));
  gel(y,11) = c6; /* 36 b2 b4 - b2^3 - 216 b6 */

  D = gsub(gmul(b4, gadd(gmulsg(9,gmul(b2,b6)),gmulsg(-8,gsqr(b4)))),
           gadd(gmul(b22,b8),gmulsg(27,gsqr(b6))));
  gel(y,12) = D;
  if (gcmp0(D)) pari_err(talker,"singular curve in ellinit");

  j = gdiv(gmul(gsqr(c4),c4), D);
  gel(y,13) = j;
}

GEN
smallinitell(GEN x)
{
  pari_sp av = avma;
  GEN y = cgetg(14,t_VEC);
  if (typ(x)==t_STR)
    x=gel(ellsearchcurve(x),2);
  smallinitell0(x,y); return gerepilecopy(av,y);
}

GEN
ellinit0(GEN x, long flag,long prec)
{
  switch(flag)
  {
    case 0: return initell(x,prec);
    case 1: return smallinitell(x);
    default: pari_err(flagerr,"ellinit");
  }
  return NULL; /* not reached */
}

void
ellprint(GEN e)
{
  pari_sp av = avma;
  long vx, vy;
  GEN z;
  checksell(e);
  vx = fetch_var(); name_var(vx, "X");
  vy = fetch_var(); name_var(vy, "Y"); z = mkvec2(pol_x[vx], pol_x[vy]);
  fprintferr("%Z - (%Z)\n", ellLHS(e, z), ellRHS(e, pol_x[vx]));
  (void)delete_var();
  (void)delete_var(); avma = av;
}

/* compute a,b such that E1: y^2 = x(x-a)(x-b) ~ E0 */
static GEN
new_coords(GEN e, GEN x, GEN *pta, GEN *ptb, int flag, long prec)
{
  GEN a, b, p1, p2, w, e1 = gmael(e,14,1), b2 = gel(e,6);
  long ty = typ(e1);

  p2 = gmul2n(gadd(gmulsg(12,e1), b2), -2); /* = (12 e1 + b2) / 4 */
  if (ty == t_PADIC)
    w = gel(e,18);
  else
  {
    GEN b4 = gel(e,7);
    if (!is_const_t(ty)) pari_err(typeer,"zell");

    /* w^2 = 2b4 + 2b2 e1 + 12 e1^2 = 4(e1-e2)(e1-e3) */
    w = sqrtr( gmul2n(gadd(b4, gmul(e1,gadd(b2, mulur(6,e1)))),1) );
    if (gsigne(p2) > 0) setsigne(w, -1);
  }
  *pta = a = gmul2n(gsub(w,p2),-2);
  *ptb = b = gmul2n(w,-1); /* = sqrt( (e1 - e2)(e1 - e3) ) */
  if (!x) return NULL;
  if (flag)
  {
    GEN d = gsub(a,b);
    p1 = gadd(x, gmul2n(gadd(gmul2n(e1,2), b2),-3));
    p1 = gmul2n(p1,-1);
    p1 = gadd(p1, gsqrt(gsub(gsqr(p1), gmul(a,d)), prec));
    return gmul(p1, gsqr(gmul2n(gaddsg(1,gsqrt(gdiv(gadd(p1,d),p1),prec)),-1)));
  }
  x = gsub(x, e1);
  p1 = gadd(x, b);
  return gmul2n(gadd(p1, gsqrt(gsub(gsqr(p1), gmul2n(gmul(a,x),2)),prec)), -1);
}

/* a1, b1 are non-0 t_REALs */
static GEN
do_agm(GEN *ptx, GEN a1, GEN b1)
{
  const long s = signe(b1), l = min(lg(a1), lg(b1)), G = 6 - bit_accuracy(l);
  GEN p1, a, b, x;

  x = gmul2n(subrr(a1,b1),-2);
  if (!signe(x)) pari_err(precer,"initell");
  for(;;)
  {
    GEN d;
    a = a1; b = b1;
    b1 = sqrtr(mulrr(a,b)); setsigne(b1, s);
    a1 = gmul2n(addrr(addrr(a,b), gmul2n(b1,1)),-2);
    d = subrr(a1,b1);
    if (!signe(d)) break;
    p1 = sqrtr( divrr(addrr(x,d),x) );
    x = mulrr(x, gsqr(addsr(1,p1)));
    setexpo(x, expo(x)-2);
    if (expo(d) <= G + expo(b1)) break;
  }
  *ptx = x; return ginv(gmul2n(a1,2));
}
/* a1, b1 are t_PADICs */
static GEN
do_padic_agm(GEN *ptx, GEN a1, GEN b1, GEN p)
{
  GEN p1, a, b, bmod1, bmod = modii(gel(b1,4),p), x = *ptx;
  long mi;

  if (!x) x = gmul2n(gsub(a1,b1),-2);
  if (gcmp0(x)) pari_err(precer,"initell");
  mi = min(precp(a1),precp(b1));
  for(;;)
  {
    GEN d;
    a = a1; b = b1;
    b1 = gprec(padic_sqrt(gmul(a,b)),mi);
    bmod1 = modii(gel(b1,4),p);
    if (!equalii(bmod1,bmod)) b1 = gneg_i(b1);
    a1 = gprec(gmul2n(gadd(gadd(a,b),gmul2n(b1,1)),-2),mi);
    d = gsub(a1,b1);
    if (gcmp0(d)) break;
    p1 = padic_sqrt(gdiv(gadd(x,d),x));
    if (! gcmp1(modii(gel(p1,4),p))) p1 = gneg_i(p1);
    x = gmul(x, gsqr(gmul2n(gaddsg(1,p1),-1)));
  }
  *ptx = x; return ginv(gmul2n(a1,2));
}

static GEN
padic_initell(GEN y, GEN p, long prec)
{
  GEN b2, b4, c4, c6, p1, w, pv, a1, b1, x1, u2, q, e0, e1;
  long i, alpha;

  for (i=1; i<=13; i++)
    if (typ(gel(y,i)) != t_PADIC) gel(y,i) = gcvtop(gel(y,i), p, prec);
  if (gcmp0(gel(y,13)) || valp(gel(y,13)) >= 0) /* p | j */
    pari_err(talker,"valuation of j must be negative in p-adic ellinit");
  if (equaliu(p,2))
  {
    pv = utoipos(4);
    pari_err(impl,"initell for 2-adic numbers");
  }
  else
    pv = p;

  b2 = gel(y,6);
  b4 = gel(y,7);
  c4 = gel(y,10);
  c6 = gel(y,11); alpha = valp(c4) >> 1;
  setvalp(c4,0);
  setvalp(c6,0);
  e1 = gdiv(c6, gmulsg(6,c4));
  c4 = gdivgs(c4,48);
  c6 = gdivgs(c6,864);
  do
  {
    GEN e2 = gsqr(e1);
    e0 = e1;  /* (c6 + 2e1^3) / (3e1^2 - c4) */
    e1 = gdiv(gadd(gmul2n(gmul(e0,e2),1),c6), gsub(gmulsg(3,e2),c4));
  }
  while (!gequal(e0,e1));
  setvalp(e1, valp(e1)+alpha);

  e1 = gsub(e1, gdivgs(b2,12));
  w = padic_sqrt(gmul2n(gadd(b4,gmul(e1,gadd(b2,gmulsg(6,e1)))),1));

  p1 = gaddgs(gdiv(gmulsg(3,e0),w),1);
  if (valp(p1) <= 0) w = gneg_i(w);
  gel(y,18) = w;

  a1 = gmul2n(gsub(w,gadd(gmulsg(3,e1),gmul2n(b2,-2))),-2);
  b1 = gmul2n(w,-1); x1 = NULL;
  u2 = do_padic_agm(&x1,a1,b1,pv);

  p1 = ginv(gmul2n(gmul(u2,x1),1));
  w = gaddsg(1,p1);
  q = padic_sqrt(gmul(p1, gaddgs(p1,2))); /* sqrt(w^2 - 1) */
  p1 = gadd(w,q);
  q = gcmp0(p1)? gsub(w,q): p1;
  if (valp(q) < 0) q = ginv(q);

  gel(y,14) = mkvec(e1);
  gel(y,15) = u2;
  gel(y,16) = ((valp(u2)&1) || kronecker(gel(u2,4),p) <= 0)? gen_0: padic_sqrt(u2);
  gel(y,17) = q;
  gel(y,19) = gen_0; return y;
}

static int
invcmp(GEN x, GEN y) { return -gcmp(x,y); }

static void
set_dummy(GEN y) {
  gel(y,14)=gel(y,15)=gel(y,16)=gel(y,17)=gel(y,18)=gel(y,19) = gen_0;
}

/* 2iPi/x, more efficient when x pure imaginary */
static GEN
PiI2div(GEN x, long prec) { return gdiv(Pi2n(1, prec), mulcxmI(x)); }
/* exp(I x y), more efficient for x in R, y pure imaginary */
static GEN
expIxy(GEN x, GEN y, long prec) { return gexp(gmul(x, mulcxI(y)), prec); }
static GEN
check_real(GEN q)
{ return (typ(q) == t_COMPLEX && gcmp0(gel(q,2)))? gel(q,1): q; }

static GEN
initell0(GEN x, long prec)
{
  GEN D, R, T, p, w, a1, b1, x1, u2, q, pi, pi2, tau, w1, w2;
  GEN y = cgetg(20,t_VEC);
  long PREC, i, e, stop = 0;

  smallinitell0(x,y);

  e = BIGINT; p = NULL;
  for (i=1; i<=5; i++)
  {
    q = gel(y,i);
    switch(typ(q)) {
      case t_PADIC:
      {
        long e2 = signe(q[4])? precp(q)+valp(q): valp(q);
        if (e2 < e) e = e2;
        if (!p) p = gel(q,2);
        else if (!equalii(p,gel(q,2)))
          pari_err(talker,"incompatible p-adic numbers in initell");
        break;
      }
      case t_INT: case t_REAL: case t_FRAC:
        break;
      default:
        stop = 1; break;
    }
  }
  if (e < BIGINT) return padic_initell(y,p,e);
  if (!prec || stop) { set_dummy(y); return y; }

  D = gel(y,12);
  switch(typ(D))
  {
    case t_INT: e = expi(D); break;
    case t_FRAC:e = max(expi(gel(D,1)), expi(gel(D,2))); break;
    default: e = -1; break;
  }
  PREC = prec;
  if (e > 0) PREC += nbits2nlong(e >> 1);
  R = cleanroots(RHSpol(y), PREC);
  /* sort roots in decreasing order */
  if (gsigne(D) > 0) R = gen_sort(R, 0, invcmp);
  gel(y,14) = R;

  (void)new_coords(y, NULL, &a1, &b1, 0, 0);
  u2 = do_agm(&x1,a1,b1); /* 1/4M */

  w = addsr(1, ginv(gmul2n(mulrr(u2,x1),1)));
  q = sqrtr( addrs(gsqr(w),-1) );
  if (signe(real_i(w)) > 0)
    q = ginv(gadd(w, q));
  else
    q = gsub(w, q);
  if (gexpo(q) >= 0) q = ginv(q);
  pi = mppi(prec); pi2 = gmul2n(pi,1);
  tau = mulcxmI( gdiv(glog(q,prec),pi2) );

  gel(y,19) = gmul(gmul(gsqr(pi2), mpabs(u2)), imag_i(tau));
  w1 = gmul(pi2, sqrtr(mpneg(u2)));
  w2 = gmul(tau, w1);
  if (signe(b1) < 0)
    q = gsqrt(q,prec);
  else
  {
    w1= gmul2n(mpabs(gel(w2,1)), 1);
    q = mpexp(mulrr(negr(pi), divrr(gel(w2,2),w1)));
    if (signe(w2[1]) < 0) setsigne(q, -1);
    q = pureimag(q);
  }
  gel(y,15) = w1;
  gel(y,16) = w2;
  T = vecthetanullk(q, 2, prec);
  if (gcmp0(gel(T,1))) pari_err(precer,"initell");
  T = check_real(gdiv(gel(T,2), gel(T,1)));
  /* pi^2 / 6w1 * theta'''(q,0) / theta'(q,0) */
  gel(y,17) = gdiv(gmul(gsqr(pi),T), gmulsg(6,w1));
  gel(y,18) = gdiv(gadd(gmul(gel(y,17),w2), mulcxmI(pi)), w1);
  return y;
}

GEN
initell(GEN x, long prec)
{
  pari_sp av = avma;
  if (typ(x)==t_STR)
    x=gel(ellsearchcurve(x),2);
  return gerepilecopy(av, initell0(x,prec));
}

/********************************************************************/
/**                                                                **/
/**                       Coordinate Change                        **/
/**                                                                **/
/********************************************************************/
/* [1,0,0,0] */
static GEN
init_ch() { 
  GEN v = cgetg(5, t_VEC);
  gel(v,1) = gen_1;
  gel(v,2) = gel(v,3) = gel(v,4) = gen_0;
  return v;
}

static GEN
coordch4(GEN e, GEN u, GEN r, GEN s, GEN t)
{
  GEN R, y, p1, p2, v, v2, v3, v4, v6, r2, b2r, rx3 = gmulsg(3,r);
  long i, lx = lg(e);

  y = cgetg(lx,t_VEC);
  v = ginv(u); v2 = gsqr(v); v3 = gmul(v,v2); v4 = gsqr(v2); v6 = gsqr(v3);
  /* A1 = (a1 + 2s) / u */
  gel(y,1) = gmul(v,gadd(gel(e,1),gmul2n(s,1)));
  /* A2 = (a2 + 3r - (a1 s + s^2)) / u^2 */
  gel(y,2) = gmul(v2,gsub(gadd(gel(e,2),rx3),gmul(s,gadd(gel(e,1),s))));
  p2 = ellLHS0(e,r);
  p1 = gadd(gmul2n(t,1), p2);
  /* A3 = (2t + a1 r + a3) / u^3 */
  gel(y,3) = gmul(v3,p1);
  p1 = gsub(gel(e,4),gadd(gmul(t,gel(e,1)),gmul(s,p1)));
  /* A4 = (a4 - (a1 t + s (2t + a1 r + a3)) + 2a2 r + 3r^2) / u^4 */
  gel(y,4) = gmul(v4,gadd(p1,gmul(r,gadd(gmul2n(gel(e,2),1),rx3))));
  /* A6 = (r^3 + a2 r^2 + a4 r + a6 - t(t + a1 r + a3)) / u^6 */
  gel(y,5) = gmul(v6,gsub(ellRHS(e,r), gmul(t,gadd(t, p2))));
  if (lx == 6) return y;
  if (lx < 14) pari_err(elliper1);

  /* B2 = (b2 + 12r) / u^2 */
  gel(y,6) = gmul(v2,gadd(gel(e,6),gmul2n(rx3,2)));
  b2r = gmul(r, gel(e,6));
  r2 = gsqr(r);
  /* B4 = (b4 + b2 r + 6r^2) / u^4 */
  gel(y,7) = gmul(v4,gadd(gel(e,7),gadd(b2r, gmulsg(6,r2))));
  /* B6 = (b6 + 2b4 r + 2b2 r^2 + 4r^3) / u^6 */
  gel(y,8) = gmul(v6,gadd(gel(e,8),gmul(r,gadd(gmul2n(gel(e,7),1),
                                            gadd(b2r,gmul2n(r2,2))))));
  /* B8 = (b8 + 3b6r + 3b4 r^2 + b2 r^3 + 3r^4) / u^8 */
  p1 = gadd(gmulsg(3,gel(e,7)),gadd(b2r, gmulsg(3,r2)));
  gel(y,9) = gmul(gsqr(v4),
              gadd(gel(e,9), gmul(r,gadd(gmulsg(3,gel(e,8)), gmul(r,p1)))));
  gel(y,10) = gmul(v4,gel(e,10));
  gel(y,11) = gmul(v6,gel(e,11));
  gel(y,12) = gmul(gsqr(v6),gel(e,12));
  gel(y,13) = gel(e,13);
  if (lx == 14) return y;
  if (lx < 20) pari_err(elliper1);
  R = gel(e,14);
  if (typ(R) != t_COL) set_dummy(y);
  else if (typ(e[1])==t_PADIC)
  {
    gel(y,14) = mkvec( gmul(v2, gsub(gel(R,1),r)) );
    gel(y,15) = gmul(gel(e,15), gsqr(u));
    gel(y,16) = gmul(gel(e,16), u);
    gel(y,17) = gel(e,17);
    gel(y,18) = gmul(gel(e,18), v2);
    gel(y,19) = gen_0;
  }
  else
  {
    p2 = cgetg(4,t_COL);
    for (i=1; i<=3; i++) gel(p2,i) = gmul(v2, gsub(gel(R,i),r));
    gel(y,14) = p2;
    gel(y,15) = gmul(gel(e,15), u);
    gel(y,16) = gmul(gel(e,16), u);
    gel(y,17) = gdiv(gel(e,17), u);
    gel(y,18) = gdiv(gel(e,18), u);
    gel(y,19) = gmul(gel(e,19), gsqr(u));
  }
  return y;
}
static GEN
_coordch(GEN e, GEN w)
{ return coordch4(e, gel(w,1), gel(w,2), gel(w,3), gel(w,4)); }
GEN
coordch(GEN e, GEN w)
{
  pari_sp av = avma;
  checkch(w); checksell(e);
  return gerepilecopy(av, _coordch(e, w));
}

/* accumulate the effects of variable changes [u,r,s,t] and [U,R,S,T], all of
 * them integers. Frequent special cases: (U = 1) or (r = s = t = 0) */
static void
cumulev(GEN *vtotal, GEN u, GEN r, GEN s, GEN t)
{
  GEN U2,U,R,S,T, v = *vtotal;
  pari_sp av;

  U = gel(v,1);
  R = gel(v,2);
  S = gel(v,3);
  T = gel(v,4);
  if (gcmp1(U))
  {
    gel(v,1) = u;
    gel(v,2) = addii(R, r);
    gel(v,3) = addii(S, s); av = avma;
    gel(v,4) = gerepileuptoint(av, addii(T, addii(t, mulii(S, r))));
  }
  else if (!signe(r) && !signe(s) && !signe(t))
    gel(v,1) = mulii(U, u);
  else /* general case */
  {
    U2 = sqri(U);
    gel(v,1) = mulii(U, u);
    gel(v,2) = addii(R, mulii(U2, r));
    gel(v,3) = addii(S, mulii(U, s));
    gel(v,4) = addii(T, mulii(U2, addii(mulii(U, t), mulii(S, r))));
  }
}
/* as above, no assumption */
static void
gcumulev(GEN *vtotal, GEN w)
{
  GEN u,r,s,t,U2,U,R,S,T, v = *vtotal;

  u = gel(w,1);
  r = gel(w,2);
  s = gel(w,3);
  t = gel(w,4);
  U = gel(v,1); gel(v,1) = gmul(U, u); U2 = gsqr(U);
  R = gel(v,2); gel(v,2) = gadd(R, gmul(U2, r));
  S = gel(v,3); gel(v,3) = gadd(S, gmul(U, s));
  T = gel(v,4); gel(v,4) = gadd(T, gmul(U2, gadd(gmul(U, t), gmul(S, r))));
}

static void
cumule(GEN *vtotal, GEN *e, GEN u, GEN r, GEN s, GEN t)
{
  *e = coordch4(*e, u, r, s, t);
  cumulev(vtotal, u, r, s, t);
}

/* X = (x-r)/u^2
 * Y = (y - s(x-r) - t) / u^3 */
static GEN
pointch0(GEN x, GEN v2, GEN v3, GEN mor, GEN s, GEN t)
{
  GEN p1,z;

  if (is_inf(x)) return x;

  z = cgetg(3,t_VEC); p1 = gadd(gel(x,1),mor);
  gel(z,1) = gmul(v2, p1);
  gel(z,2) = gmul(v3, gsub(gel(x,2), gadd(gmul(s,p1),t)));
  return z;
}

GEN
pointch(GEN x, GEN ch)
{
  GEN y, v, v2, v3, mor, r, s, t, u;
  long tx, i, lx = lg(x);
  pari_sp av = avma;

  checkpt(x); checkch(ch);
  if (lx < 2) return gcopy(x);
  u = gel(ch,1);
  r = gel(ch,2);
  s = gel(ch,3);
  t = gel(ch,4);
  v = ginv(u); v2 = gsqr(v); v3 = gmul(v,v2);
  mor = gneg_i(r);
  tx = typ(x[1]);
  if (is_matvec_t(tx))
  {
    y = cgetg(lx,tx);
    for (i=1; i<lx; i++)
      gel(y,i) = pointch0(gel(x,i),v2,v3,mor,s,t);
  }
  else
    y = pointch0(x,v2,v3,mor,s,t);
  return gerepilecopy(av,y);
}

/* x =  u^2*X +r
 * y =  u^3*Y +s*u^2*X+t */
static GEN
pointchinv0(GEN x, GEN u2, GEN u3, GEN r, GEN s, GEN t)
{
  GEN u2X, z;
  GEN X=gel(x,1), Y=gel(x,2);
  if (is_inf(x)) return x;

  u2X = gmul(u2,X);
  z = cgetg(3, t_VEC); 
  gel(z,1) = gadd(u2X, r);
  gel(z,2) = gadd(gmul(u3, Y), gadd(gmul(s, u2X), t));
  return z;
}

GEN
pointchinv(GEN x, GEN ch)
{
  GEN y, u, r, s, t, u2, u3;
  long tx, i, lx = lg(x);
  pari_sp av = avma;

  checkpt(x); checkch(ch);
  if (lx < 2) return gcopy(x);
  u = gel(ch,1);
  r = gel(ch,2);
  s = gel(ch,3);
  t = gel(ch,4);
  tx = typ(x[1]);
  u2=gsqr(u); u3=gmul(u,u2);
  if (is_matvec_t(tx))
  {
    y = cgetg(lx,tx);
    for (i=1; i<lx; i++)
      gel(y,i) = pointchinv0(gel(x,i),u2,u3,r,s,t);
  }
  else
    y = pointchinv0(x,u2,u3,r,s,t);
  return gerepilecopy(av,y);
}


static long
ellexpo(GEN E)
{
  long i, f, e = -(long)HIGHEXPOBIT;
  for (i=1; i<=5; i++)
  {
    f = gexpo(gel(E,i));
    if (f > e) e = f;
  }
  return e;
}

/* Exactness of lhs and rhs in the following depends in non-obvious ways
 * on the coeffs of the curve as well as on the components of the point z.
 * Thus if e is exact, with a1==0, and z has exact y coordinate only, the
 * lhs will be exact but the rhs won't. */
int
oncurve(GEN e, GEN z)
{
  GEN LHS, RHS, x;
  long pl, pr, ex, expx;
  pari_sp av;

  checkpt(z); if (is_inf(z)) return 1; /* oo */
  av = avma;
  LHS = ellLHS(e,z);
  RHS = ellRHS(e,gel(z,1)); x = gsub(LHS,RHS);
  if (gcmp0(x)) { avma = av; return 1; }
  pl = precision(LHS);
  pr = precision(RHS);
  if (!pl && !pr) { avma = av; return 0; } /* both of LHS, RHS are exact */
  /* at least one of LHS,RHS is inexact */
  ex = pr? gexpo(RHS): gexpo(LHS); /* don't take exponent of exact 0 */
  if (!pr || (pl && pl < pr)) pr = pl; /* min among nonzero elts of {pl,pr} */
  expx = gexpo(x);
  pr = (expx < ex - bit_accuracy(pr) + 15
     || expx < ellexpo(e) - bit_accuracy(pr) + 5);
  avma = av; return pr;
}

GEN
ellisoncurve(GEN e, GEN a)
{
  long i, tx = typ(a), lx = lg(a);

  checksell(e); if (!is_vec_t(tx)) pari_err(elliper1);
  lx = lg(a); if (lx==1) return cgetg(1,tx);
  tx = typ(a[1]);
  if (is_vec_t(tx))
  {
    GEN z = cgetg(lx,tx);
    for (i=1; i<lx; i++) gel(z,i) = ellisoncurve(e,gel(a,i));
    return z;
  }
  return oncurve(e, a)? gen_1: gen_0;
}

GEN
addell(GEN e, GEN z1, GEN z2)
{
  GEN p1, p2, x, y, x1, x2, y1, y2;
  pari_sp av = avma, tetpil;

  checksell(e); checkpt(z1); checkpt(z2);
  if (is_inf(z1)) return gcopy(z2);
  if (is_inf(z2)) return gcopy(z1);

  x1 = gel(z1,1); y1 = gel(z1,2);
  x2 = gel(z2,1); y2 = gel(z2,2);
  if (x1 == x2 || gequal(x1,x2))
  { /* y1 = y2 or -LHS0-y2 */
    if (y1 != y2)
    {
      int eq;
      if (precision(y1) || precision(y2))
        eq = (gexpo(gadd(ellLHS0(e,x1),gadd(y1,y2))) >= gexpo(y1));
      else
        eq = gequal(y1,y2);
      if (!eq) { avma = av; return mkvec(gen_0); }
    }
    p2 = d_ellLHS(e,z1);
    if (gcmp0(p2)) { avma = av; return mkvec(gen_0); }
    p1 = gadd(gsub(gel(e,4),gmul(gel(e,1),y1)),
              gmul(x1,gadd(gmul2n(gel(e,2),1),gmulsg(3,x1))));
  }
  else {
    p1 = gsub(y2,y1);
    p2 = gsub(x2,x1);
  }
  p1 = gdiv(p1,p2);
  x = gsub(gmul(p1,gadd(p1,gel(e,1))), gadd(gadd(x1,x2),gel(e,2)));
  y = gadd(gadd(y1, ellLHS0(e,x)), gmul(p1,gsub(x,x1)));
  tetpil = avma; p1 = cgetg(3,t_VEC);
  gel(p1,1) = gcopy(x);
  gel(p1,2) = gneg(y); return gerepile(av,tetpil,p1);
}

static GEN
invell(GEN e, GEN z)
{
  GEN t;
  if (is_inf(z)) return z;
  t = cgetg(3,t_VEC);
  t[1] = z[1];
  gel(t,2) = gneg_i(gadd(gel(z,2), ellLHS0(e,gel(z,1))));
  return t;
}

GEN
subell(GEN e, GEN z1, GEN z2)
{
  pari_sp av = avma;
  checksell(e); checkpt(z2);
  return gerepileupto(av, addell(e, z1, invell(e,z2)));
}

GEN
ordell(GEN e, GEN x, long prec)
{
  long td, i, tx = typ(x);
  pari_sp av = avma;
  GEN D, a, b, d, y;

  checksell(e);
  if (is_matvec_t(tx))
  {
    long lx = lg(x); y = cgetg(lx,tx);
    for (i=1; i<lx; i++) gel(y,i) = ordell(e,gel(x,i),prec);
    return y;
  }

  a = ellRHS(e,x);
  b = ellLHS0(e,x); /* y*(y+b) = a */
  D = gadd(gsqr(b), gmul2n(a,2));
  td = typ(D);
  if (td == t_INTMOD && equaliu(gel(D,1), 2))
  { /* curve over F_2 */
    avma = av;
    if (!signe(D[2])) {
      y = cgetg(2,t_VEC);
      gel(y,1) = mkintmodu(gcmp0(a)?0:1, 2);
    } else {
      if (!gcmp0(a)) return cgetg(1,t_VEC);
      y = cgetg(3,t_VEC);
      gel(y,1) = mkintmodu(0,2);
      gel(y,2) = mkintmodu(1,2);
    }
    return y;
  }

  if (gcmp0(D)) {
    b = gneg_i(b);
    y = cgetg(2,t_VEC);
    gel(y,1) = gmul2n(b,-1);
    return gerepileupto(av,y);
  }
  switch(td)
  {
    case t_INT:
      if (!Z_issquarerem(D,&d)) { avma = av; return cgetg(1,t_VEC); }
      break;
    case t_FRAC:
      if (gissquarerem(D,&d) == gen_0) { avma = av; return cgetg(1,t_VEC); }
      break;
    case t_INTMOD:
      if (kronecker(gel(D,2),gel(D,1)) < 0) {
        avma = av; return cgetg(1,t_VEC);
      } /* fall through */
    default:
      d = gsqrt(D,prec);
  }
  a = gsub(d,b); y = cgetg(3,t_VEC);
  gel(y,1) = gmul2n(a, -1);
  gel(y,2) = gsub(gel(y,1),d);
  return gerepileupto(av,y);
}

/* n t_QUAD */
static GEN
CM_ellpow(GEN e, GEN z, GEN n)
{
  GEN x, y, p0, p1, q0, q1, z1, z2, pol, grdx, b2ov12;
  long ln, ep, vn;
  pari_sp av = avma;

  if (is_inf(z)) return gcopy(z);
  pol = gel(n,1);
  if (signe(pol[2]) < 0) pari_err(typeer,"CM_ellpow");
  if (typ(n[2]) != t_INT || typ(n[3]) != t_INT)
    pari_err(impl, "powell for nonintegral CM exponent");

  ln = itos_or_0( shifti(addsi(1, quadnorm(n)), 2) );
  if (!ln) pari_err(talker, "norm too large in CM");
  vn = (ln-4)>>2;
  z1 = weipell(e, ln);
  z2 = gsubst(z1, 0, monomial(n, 1, 0));
  b2ov12 = gdivgs(gel(e,6), 12); /* x - b2/12 */
  grdx = gadd(gel(z,1), b2ov12);
  p0 = gen_0; p1 = gen_1;
  q0 = gen_1; q1 = gen_0;
  do
  {
    GEN p2,q2, ss = gen_0;
    do
    {
      ep = (-valp(z2)) >> 1;
      ss = gadd(ss, gmul(gel(z2,2), monomial(gen_1, ep, 0)));
      z2 = gsub(z2, gmul(gel(z2,2), gpowgs(z1, ep)));
    }
    while (valp(z2) <= 0);
    p2 = gadd(p0, gmul(ss,p1)); p0 = p1; p1 = p2;
    q2 = gadd(q0, gmul(ss,q1)); q0 = q1; q1 = q2;
    if (!signe(z2)) break;
    z2 = ginv(z2);
  }
  while (degpol(p1) < vn);
  if (degpol(p1) > vn || signe(z2))
    pari_err(talker,"not a complex multiplication in powell");
  x = gdiv(p1,q1);
  y = gdiv(deriv(x,0),n);
  x = gsub(poleval(x,grdx), b2ov12);
  y = gsub( gmul(d_ellLHS(e,z), poleval(y,grdx)), ellLHS0(e,x));
  z = cgetg(3,t_VEC);
  gel(z,1) = gcopy(x);
  gel(z,2) = gmul2n(y,-1); return gerepileupto(av, z);
}

static GEN
_sqr(void *e, GEN x) { return addell((GEN)e, x, x); }
static GEN
_mul(void *e, GEN x, GEN y) { return addell((GEN)e, x, y); }

GEN
powell(GEN e, GEN z, GEN n)
{
  pari_sp av = avma;
  long s;

  checksell(e); checkpt(z);
  if (typ(n)==t_QUAD) return CM_ellpow(e,z,n);
  if (typ(n) != t_INT) pari_err(impl,"powell for non integral, non CM, exponents");
  s = signe(n);
  if (!s || is_inf(z)) return mkvec(gen_0);
  if (s < 0) z = invell(e,z);
  if (is_pm1(n)) return s < 0? gerepilecopy(av, z): gcopy(z);
  return gerepileupto(av, leftright_pow(z, n, (void*)e, &_sqr, &_mul));
}

/********************************************************************/
/**                                                                **/
/**                       ELLIPTIC FUNCTIONS                       **/
/**                                                                **/
/********************************************************************/
static GEN
quot(GEN x, GEN y)
{
  GEN z = mpdiv(x, y), q = floorr(z);
  if (gsigne(y) < 0 && !gequal(z, q)) q = addis(q, 1);
  return q;
}
GEN
zell(GEN e, GEN z, long prec)
{
  long ty, sw, fl;
  pari_sp av = avma;
  GEN t, u, p1, p2, a, b, x1, u2, D = gel(e,12);

  checkbell(e); checkpt(z);
  ty = typ(D); if (ty == t_INTMOD) pari_err(typeer,"zell");
  if (is_inf(z)) return (ty==t_PADIC)? gen_1: gen_0;

  x1 = new_coords(e,gel(z,1),&a,&b,1, prec);
  if (ty==t_PADIC)
  {
    u2 = do_padic_agm(&x1,a,b,gel(D,2));
    if (!gcmp0(gel(e,16)))
    {
      t = padic_sqrt(gaddsg(1, gdiv(x1,a)));
      t = gdiv(gaddsg(-1,t), gaddsg(1,t));
    }
    else t = gaddsg(2, ginv(gmul(u2,x1)));
    return gerepileupto(av,t);
  }

  sw = gsigne(real_i(b)); fl=0;
  for(;;) /* ~ agm */
  {
    GEN a0 = a, b0 = b, x0 = x1, d;

    b = gsqrt(gmul(a0,b0),prec);
    if (gsigne(real_i(b)) != sw) b = gneg_i(b);
    a = gmul2n(gadd(gadd(a0,b0),gmul2n(b,1)),-2);
    d = gsub(a,b);
    if (gcmp0(d) || gexpo(d) < gexpo(a) - bit_accuracy(prec) + 5) break;
    p1 = gsqrt(gdiv(gadd(x0,d),x0),prec);
    x1 = gmul(x0,gsqr(gmul2n(gaddsg(1,p1),-1)));
    d = gsub(x1,x0);
    if (gcmp0(d) || gexpo(d) < gexpo(x1) - bit_accuracy(prec) + 5)
    {
      if (fl) break;
      fl = 1;
    }
    else fl = 0;
  }
  u = gdiv(x1,a); t = gaddsg(1,u);
  if (gcmp0(t) || gexpo(t) <  5 - bit_accuracy(prec))
    t = gen_m1;
  else
    t = gdiv(u,gsqr(gaddsg(1,gsqrt(t,prec))));
  u = gsqrt(ginv(gmul2n(a,2)),prec);
  t = gmul(u, glog(t,prec));

  /* which square root? test the reciprocal function (pointell) */
  if (!gcmp0(t))
  {
    GEN z1,z2;
    int bad;

    z1 = pointell(e,gprec_w(t,3),3); /* we don't need much precision */
    /* Either z = z1 (ok: keep t), or z = z2 (bad: t <-- -t) */
    z2 = invell(e, z1);
    bad = (gexpo(gsub(z,z1)) > gexpo(gsub(z,z2)));
    if (bad) t = gneg(t);
    if (DEBUGLEVEL) {
      if (DEBUGLEVEL>4) {
        fprintferr("  z  = %Z\n",z);
        fprintferr("  z1 = %Z\n",z1);
        fprintferr("  z2 = %Z\n",z2);
      }
      fprintferr("ellpointtoz: %s square root\n", bad? "bad": "good");
      flusherr();
    }
  }
  /* send t to the fundamental domain if necessary */
  p2 = quot(imag_i(t), imag_i(gel(e,16)));
  if (signe(p2)) t = gsub(t, gmul(p2, gel(e,16)));
  p2 = quot(real_i(t), gel(e,15));
  if (signe(p2)) t = gsub(t, gmul(p2, gel(e,15)));
  return gerepileupto(av,t);
}

typedef struct {
  GEN w1,w2,tau; /* original basis for L = <w1,w2> = w2 <1,tau> */
  GEN W1,W2,Tau; /* new basis for L = <W1,W2> = W2 <1,tau> */
  GEN a,b,c,d; /* tau in F = h/Sl2, tau = g.t, g=[a,b;c,d] in SL(2,Z) */
  GEN x,y; /* z/w2 defined mod <1,tau> --> z + x tau + y reduced mod <1,tau> */
  int swap; /* 1 if we swapped w1 and w2 */
} SL2_red;

/* compute gamma in SL_2(Z) gamma(t) is in the usual
   fundamental domain. Internal function no check, no garbage. */
static void
set_gamma(SL2_red *T)
{
  GEN t = T->tau, a, b, c, d, run = dbltor(1. - 1e-8);

  a = d = gen_1;
  b = c = gen_0;
  for(;;)
  {
    GEN m, p1, n = ground(real_i(t));
    if (signe(n))
    { /* apply T^n */
      n = negi(n); t = gadd(t,n);
      a = addii(a, mulii(n,c));
      b = addii(b, mulii(n,d));
    }
    m = cxnorm(t); if (gcmp(m,run) > 0) break;
    t = gneg_i(gdiv(gconj(t), m)); /* apply S */
    p1 = negi(c); c = a; a = p1;
    p1 = negi(d); d = b; b = p1;
  }
  T->a = a; T->b = b;
  T->c = c; T->d = d;
}

/* swap w1, w2 so that Im(t := w1/w2) > 0. Set tau = representative of t in
 * the standard fundamental domain, and g in Sl_2, such that tau = g.t */
static void
red_modSL2(SL2_red *T)
{
  long s;
  T->tau = gdiv(T->w1,T->w2);
  s = gsigne(imag_i(T->tau));
  if (!s) pari_err(talker,"w1 and w2 R-linearly dependent in elliptic function");
  T->swap = (s < 0);
  if (T->swap) { swap(T->w1, T->w2); T->tau = ginv(T->tau); }
  set_gamma(T);
  /* update lattice */
  T->W1 = gadd(gmul(T->a,T->w1), gmul(T->b,T->w2));
  T->W2 = gadd(gmul(T->c,T->w1), gmul(T->d,T->w2));
  T->Tau = gdiv(T->W1, T->W2);
}

static int
get_periods(GEN e, SL2_red *T)
{
  long tx = typ(e);
  if (is_vec_t(tx))
    switch(lg(e))
    {
      case  3: T->w1 = gel(e,1);  T->w2 = gel(e,2); red_modSL2(T); return 1;
      case 20: T->w1 = gel(e,15); T->w2 = gel(e,16);red_modSL2(T); return 1;
    }
  return 0;
}

/* Return E_k(tau). Slow if tau is not in standard fundamental domain */
static GEN
trueE(GEN tau, long k, long prec)
{
  pari_sp lim, av;
  GEN p1, q, y, qn;
  long n = 1;

  q = expIxy(Pi2n(1, prec), tau, prec);
  q = check_real(q);
  y = gen_0;
  av = avma; lim = stack_lim(av,2); qn = gen_1;
  for(;; n++)
  { /* compute y := sum_{n>0} n^(k-1) q^n / (1-q^n) */
    qn = gmul(q,qn);
    p1 = gdiv(gmul(powuu(n,k-1),qn), gsub(gen_1,qn));
    if (gcmp0(p1) || gexpo(p1) <= - bit_accuracy(prec) - 5) break;
    y = gadd(y, p1);
    if (low_stack(lim, stack_lim(av,2)))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"elleisnum");
      gerepileall(av, 2, &y,&qn);
    }
  }
  return gadd(gen_1, gmul(y, gdiv(gen_2, szeta(1-k, prec))));
}

/* (2iPi/W2)^k E_k(W1/W2) */
static GEN
_elleisnum(SL2_red *T, long k, long prec)
{
  GEN y = trueE(T->Tau, k, prec);
  y = gmul(y, gpowgs(mulcxI(gdiv(Pi2n(1,prec), T->W2)),k));
  return check_real(y);
}

/* Return (2iPi)^k E_k(L) = (2iPi/w2)^k E_k(tau), with L = <w1,w2>, k > 0 even
 * E_k(tau) = 1 + 2/zeta(1-k) * sum(n>=1, n^(k-1) q^n/(1-q^n))
 * If flag is != 0 and k=4 or 6, compute g2 = E4/12 or g3 = -E6/216 resp. */
GEN
elleisnum(GEN om, long k, long flag, long prec)
{
  pari_sp av = avma;
  GEN p1, y;
  SL2_red T;

  if (k&1 || k<=0) pari_err(talker,"k not a positive even integer in elleisnum");
  if (!get_periods(om, &T)) pari_err(typeer,"elleisnum");
  y = _elleisnum(&T, k, prec);
  if (k==2 && signe(T.c))
  {
    p1 = gmul(Pi2n(1,prec), mulsi(12, T.c));
    y = gsub(y, mulcxI(gdiv(p1, gmul(T.w2, T.W2))));
  }
  else if (k==4 && flag) y = gdivgs(y,  12);
  else if (k==6 && flag) y = gdivgs(y,-216);
  return gerepileupto(av,y);
}

/* return quasi-periods associated to [w1,w2] */
static GEN
_elleta(SL2_red *T, long prec)
{
  GEN y, y1, y2, e2 = gdivgs(_elleisnum(T,2,prec), 12);
  y2 = gmul(T->W2, e2);
  y1 = gadd(PiI2div(T->W2, prec), gmul(T->W1,e2));
  y = cgetg(3,t_VEC);
  gel(y,1) = gneg(y1);
  gel(y,2) = gneg(y2); return y;
}

/* compute eta1, eta2 */
GEN
elleta(GEN om, long prec)
{
  pari_sp av = avma;
  GEN y1, y2, E2, pi = mppi(prec);
  SL2_red T;
  if (!get_periods(om, &T)) pari_err(typeer,"elleta");
  E2 = trueE(T.Tau, 2, prec); /* E_2(Tau) */
  if (signe(T.c))
  {
    GEN u = gdiv(T.w2, T.W2);
    /* E2 := u^2 E2 + 6iuc/pi = E_2(tau) */
    E2 = gadd(gmul(gsqr(u), E2), mulcxI(gdiv(gmul(mulsi(6,T.c), u), pi)));
  }
  y2 = gdiv(gmul(E2, gsqr(pi)), gmulsg(3, T.w2));
  if (T.swap)
  {
    y1 = y2;
    y2 = gadd(gmul(T.tau,y1), PiI2div(T.w2, prec));
  }
  else
    y1 = gsub(gmul(T.tau,y2), PiI2div(T.w2, prec));
  return gerepilecopy(av, mkvec2(y1,y2));
}

static GEN
reduce_z(GEN z, SL2_red *T)
{
  GEN Z = gdiv(z, T->W2);
  long t = typ(z), pr;

  if (!is_scalar_t(t) || t == t_INTMOD || t == t_PADIC || t == t_POLMOD)
    pari_err(typeer,"reduction mod SL2 (reduce_z)");
  T->x = ground(gdiv(imag_i(Z), imag_i(T->Tau)));
  Z = gsub(Z, gmul(T->x,T->Tau));
  T->y = ground(real_i(Z));
  Z = gsub(Z, T->y);
  pr = gprecision(Z);
  if (gcmp0(Z) || (pr && gexpo(Z) < 5 - bit_accuracy(pr))) Z = NULL; /*z in L*/
  return Z;
}

/* computes the numerical value of wp(z | L), L = om1 Z + om2 Z
 * return NULL if z in L.  If flall=1, compute also wp' */
static GEN
weipellnumall(SL2_red *T, GEN z, long flall, long prec)
{
  long toadd;
  pari_sp av=avma, lim, av1;
  GEN p1, pi2, q, u, y, yp, u1, u2, qn, v;

  z = reduce_z(z, T);
  if (!z) return NULL;

  /* Now L,z normalized to <1,tau>. z in fund. domain of <1, tau> */
  pi2 = Pi2n(1, prec);
  q = expIxy(pi2, T->Tau, prec);
  u = expIxy(pi2, z, prec);
  u1= gsub(gen_1,u); u2 = gsqr(u1);
  y = gadd(mkfrac(gen_1, utoipos(12)), gdiv(u,u2));
  if (flall) yp = gdiv(gadd(gen_1,u), gmul(u1,u2));
  toadd = (long)ceil(9.065*gtodouble(imag_i(z)));

  av1 = avma; lim = stack_lim(av1,1); qn = q;
  for(;;)
  {
    GEN qnu,qnu1,qnu2,qnu3,qnu4;

    qnu = gmul(qn,u);     /* q^n u */
    qnu1 = gsub(gen_1,qnu); /* 1 - q^n u */
    qnu2 = gsqr(qnu1);    /* (1 - q^n u)^2 */
    qnu3 = gsub(qn,u);    /* q^n - u */
    qnu4 = gsqr(qnu3);    /* (q^n - u)^2 */
    p1 = gsub(gmul(u, gadd(ginv(qnu2),ginv(qnu4))),
              gmul2n(ginv(gsqr(gsub(gen_1,qn))), 1));
    y = gadd(y, gmul(qn,p1));
    if (flall)
    {
      p1 = gadd(gdiv(gadd(gen_1,qnu),gmul(qnu1,qnu2)),
                gdiv(gadd(qn,u),gmul(qnu3,qnu4)));

      yp = gadd(yp, gmul(qn,p1));
    }
    qn = gmul(q,qn);
    if (gexpo(qn) <= - bit_accuracy(prec) - 5 - toadd) break;
    if (low_stack(lim, stack_lim(av1,1)))
    {
      GEN *gptr[3]; gptr[0]=&y; gptr[1]=&qn; gptr[2]=&yp;
      if(DEBUGMEM>1) pari_warn(warnmem,"weipellnum");
      gerepilemany(av1,gptr,flall?3:2);
    }
  }

  u1 = gdiv(pi2, mulcxmI(T->W2));
  u2 = gsqr(u1);
  y = gmul(u2,y); /* y *= (2i pi / w2)^2 */
  if (flall)
  {
    yp = gmul(u, gmul(gmul(u1,u2),yp));/* yp *= u (2i pi / w2)^3 */
    v = mkvec2(y, gmul2n(yp,-1));
  }
  else v = y;
  return gerepilecopy(av, v);
}

GEN
ellzeta(GEN om, GEN z, long prec)
{
  long toadd;
  pari_sp av = avma, lim, av1;
  GEN Z, pi2, q, u, y, qn, et = NULL;
  SL2_red T;

  if (!get_periods(om, &T)) pari_err(typeer,"ellzeta");
  Z = reduce_z(z, &T);
  if (!Z) pari_err(talker,"can't evaluate ellzeta at a pole");
  if (!gcmp0(T.x) || !gcmp0(T.y))
  {
    et = _elleta(&T,prec);
    et = gadd(gmul(T.x,gel(et,1)), gmul(T.y,gel(et,2)));
  }

  pi2 = Pi2n(1, prec);
  q = expIxy(pi2, T.Tau, prec);
  u = expIxy(pi2, Z, prec);

  y = mulcxmI(gdiv(gmul(gsqr(T.W2),_elleisnum(&T,2,prec)), pi2));
  y = gadd(ghalf, gdivgs(gmul(Z,y),-12));
  y = gadd(y, ginv(gsubgs(u, 1)));
  toadd = (long)ceil(9.065*gtodouble(imag_i(Z)));
  av1 = avma; lim = stack_lim(av1,1);

  /* y += sum q^n ( u/(u*q^n - 1) + 1/(u - q^n) ) */
  for (qn = q;;)
  {
    GEN p1 = gadd(gdiv(u,gsub(gmul(qn,u),gen_1)), ginv(gsub(u,qn)));
    y = gadd(y, gmul(qn,p1));
    qn = gmul(q,qn);
    if (gexpo(qn) <= - bit_accuracy(prec) - 5 - toadd) break;
    if (low_stack(lim, stack_lim(av1,1)))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"ellzeta");
      gerepileall(av1,2, &y,&qn);
    }
  }
  y = mulcxI(gmul(gdiv(pi2,T.W2), y));
  return et? gerepileupto(av, gadd(y,et)): gerepilecopy(av, y);
}

/* if flag=0, return ellsigma, otherwise return log(ellsigma) */
GEN
ellsigma(GEN w, GEN z, long flag, long prec)
{
  long toadd;
  pari_sp av = avma, lim, av1;
  GEN Z, zinit, p1, pi, pi2, q, u, y, y1, u1, qn, uinv, et, etnew, uhalf;
  int doprod = (flag >= 2), dolog = (flag & 1);
  SL2_red T;

  if (!get_periods(w, &T)) pari_err(typeer,"ellsigma");
  Z = reduce_z(z, &T);
  if (!Z)
  {
    if (!dolog) return gen_0;
    pari_err(talker,"can't evaluate log(ellsigma) at lattice point");
  }
  et = _elleta(&T, prec);
  etnew = gadd(gmul(T.x,gel(et,1)), gmul(T.y,gel(et,2)));

  pi2 = Pi2n(1,prec);
  pi  = mppi(prec);
  zinit = gmul(Z,T.W2);
  p1 = gadd(zinit, gmul2n(gadd(gmul(T.x,T.W1), gmul(T.y,T.W2)),-1));
  etnew = gmul(etnew, p1);
  if (mpodd(T.x) || mpodd(T.y)) etnew = gadd(etnew, mulcxI(pi));

  y1 = gadd(etnew, gmul2n(gmul(gmul(Z,zinit),gel(et,2)),-1));

  toadd = (long)ceil((2*PI/LOG2) * fabs(gtodouble(imag_i(Z))));
  uhalf = expIxy(pi, Z, prec); /* exp(i Pi Z) */
  u = gsqr(uhalf);
  if (doprod) { /* use product */
    q = expIxy(pi2, T.Tau, prec);
    uinv = ginv(u);
    u1 = gsub(uhalf,ginv(uhalf));
    y = mulcxmI(gdiv(gmul(T.W2,u1), pi2));
    av1 = avma; lim = stack_lim(av1,1); qn=q;
    for(;;)
    {
      p1 = gmul(gadd(gmul(qn,u),gen_m1),gadd(gmul(qn,uinv),gen_m1));
      p1 = gdiv(p1,gsqr(gadd(qn,gen_m1)));
      y = gmul(y,p1);
      qn = gmul(q,qn);
      if (gexpo(qn) <= - bit_accuracy(prec) - 5 - toadd) break;
      if (low_stack(lim, stack_lim(av1,1)))
      {
        if(DEBUGMEM>1) pari_warn(warnmem,"ellsigma");
        gerepileall(av1,2, &y,&qn);
      }
    }
  } else { /* use sum */
    GEN q8, qn2, urn, urninv;
    long n;
    q8 = expIxy(gmul2n(pi2,-3), T.Tau, prec);
    q = gpowgs(q8,8);
    u = gneg_i(u); uinv = ginv(u);
    y = gen_0;
    av1 = avma; lim = stack_lim(av1,1);
    qn = q; qn2 = gen_1;
    urn = uhalf; urninv = ginv(uhalf);
    for(n=0;;n++)
    {
      y = gadd(y,gmul(qn2,gsub(urn,urninv)));
      qn2 = gmul(qn,qn2);
      qn  = gmul(q,qn);
      urn = gmul(urn,u);
      urninv = gmul(urninv,uinv);
      if (gexpo(qn2) + n*toadd <= - bit_accuracy(prec) - 5) break;
      if (low_stack(lim, stack_lim(av1,1)))
      {
        if(DEBUGMEM>1) pari_warn(warnmem,"ellsigma");
        gerepileall(av1,5, &y,&qn,&qn2,&urn,&urninv);
      }
    }
    p1 = gmul(gmul(y,q8),
              gdiv(mulcxmI(T.W2), gmul(pi2,gpowgs(trueeta(T.Tau,prec),3))));
  }
  y1 = dolog? gadd(y1, glog(p1,prec)): gmul(p1, gexp(y1,prec));
  return gerepileupto(av, y1);
}

GEN
pointell(GEN e, GEN z, long prec)
{
  pari_sp av = avma;
  GEN v;
  SL2_red T;

  checkbell(e); (void)get_periods(e, &T);
  v = weipellnumall(&T,z,1,prec);
  if (!v) { avma = av; return mkvec(gen_0); }
  gel(v,1) = gsub(gel(v,1), gdivgs(gel(e,6),12));
  gel(v,2) = gsub(gel(v,2), gmul2n(ellLHS0(e,gel(v,1)),-1));
  return gerepilecopy(av, v);
}

static GEN
_weipell(GEN c4, GEN c6, long PREC)
{
  long i, k, l, precres = 2*PREC;
  pari_sp av;
  GEN t, res = cgetg(precres+2,t_SER), *P = (GEN*)(res + 2);

  res[1] = evalsigne(1) | evalvalp(-2) | evalvarn(0);
  if (!PREC) { setsigne(res,0); return res; }

  for (i=1; i<precres; i+=2) P[i]= gen_0;
  switch(PREC)
  {
    default:P[6] = gdivgs(c6,6048);
    case 3: P[4] = gdivgs(c4, 240);
    case 2: P[2] = gen_0;
    case 1: P[0] = gen_1;
    case 0: break;
  }
  if (PREC == 4) return res;
  av = avma;
  P[8] = gerepileupto(av, gdivgs(gsqr(P[4]), 3));
  for (k=5; k<PREC; k++)
  {
    av = avma;
    t = gmul(P[4], P[(k-2)<<1]);
    for (l=3; (l<<1) < k; l++) t = gadd(t, gmul(P[l<<1], P[(k-l)<<1]));
    t = gmul2n(t, 1);
    if ((k & 1) == 0) t = gadd(gsqr(P[k]), t);
    if (k % 3 == 2)
      t = gdivgs(gmulsg(3, t), (k-3)*(2*k+1));
    else
      t = gdivgs(t, ((k-3)*(2*k+1)) / 3);
    P[k<<1] = gerepileupto(av, t);
  }
  return res;
}

GEN
weipell(GEN e, long PREC)
{
  GEN c4 = gel(e,10);
  GEN c6 = gel(e,11);
  checkell(e); return _weipell(c4,c6,PREC);
}

GEN
weipell0(GEN e, long prec, long PREC)
{
  GEN c4,c6;

  if (lg(e) > 3) return weipell(e, PREC);
  c4 = elleisnum(e, 4, 0, prec);
  c6 = elleisnum(e, 6, 0, prec); c6 = gneg(c6);
  return _weipell(c4,c6,PREC);
}

/* assume x a t_POL */
static int
is_simple_var(GEN x)
{
  return (degpol(x) == 1 && gcmp0(gel(x,2)) && gcmp1(gel(x,3)));
}

GEN
ellwp0(GEN w, GEN z, long flag, long prec, long PREC)
{
  GEN v;
  pari_sp av = avma;
  SL2_red T;

  if (!z) return weipell0(w,prec,PREC);
  if (typ(z)==t_POL)
  {
    if (!is_simple_var(z)) pari_err(talker,"expecting a simple variable in ellwp");
    v = weipell0(w,prec,PREC); setvarn(v, varn(z));
    return v;
  }
  if (!get_periods(w, &T)) pari_err(typeer,"ellwp");
  switch(flag)
  {
    case 0: v = weipellnumall(&T,z,0,prec);
      if (!v) { avma = av; v = gpowgs(z,-2); }
      return v;
    case 1: v = weipellnumall(&T,z,1,prec);
      if (!v)
      {
        GEN p1 = gmul2n(gpowgs(z,3),1);
        pari_sp tetpil = avma;
        v = cgetg(3,t_VEC);
	gel(v,1) = gpowgs(z,-2);
	gel(v,2) = gneg(p1); return gerepile(av,tetpil,v);
      }
      return v;
    case 2: return pointell(w,z,prec);
    default: pari_err(flagerr,"ellwp"); return NULL;
  }
}

/********************************************************************/
/**                                                                **/
/**                 Tate's algorithm e (cf Anvers IV)              **/
/**               Kodaira types, global minimal model              **/
/**                                                                **/
/********************************************************************/

/* Given an integral elliptic curve in ellinit form, and a prime p, returns the
  type of the fiber at p of the Neron model, as well as the change of variables
  in the form [f, kod, v, c].

  * The integer f is the conductor's exponent.

  * The integer kod is the Kodaira type using the following notation:
    II , III , IV  -->  2, 3, 4
    I0  -->  1
    Inu --> 4+nu for nu > 0
  A '*' negates the code (e.g I* --> -2)

  * v is a quadruple [u, r, s, t] yielding a minimal model

  * c is the Tamagawa number.

  Uses Tate's algorithm (Anvers IV). Given the remarks at the bottom of
  page 46, the "long" algorithm is used for p = 2,3 only. */
static GEN
localred_result(long f, long kod, long c, GEN v)
{
  GEN z = cgetg(5, t_VEC);
  gel(z,1) = stoi(f);
  gel(z,2) = stoi(kod);
  gel(z,3) = gcopy(v);
  gel(z,4) = stoi(c); return z;
}
static GEN
localredbug(GEN p, char *s)
{
  if (BSW_psp(p)) pari_err(bugparier, s);
  pari_err(talker,"not a prime in localred");
  return NULL; /* not reached */
}

/* Here p > 3. e assumed integral */
static GEN
localred_p(GEN e, GEN p, int minim)
{
  long k, f, kod, c, nuj, nuD;
  GEN p2, v = init_ch();
  GEN c4, c6, D, tri;

  c4 = gel(e,10);
  c6 = gel(e,11);
  D  = gel(e,12);
  nuj = gcmp0(gel(e,13))? 0: - ggval(gel(e,13), p);
  nuD = Z_pval(D, p);
  k = (nuj > 0 ? nuD - nuj : nuD) / 12;
  if (k <= 0)
  {
    if (minim) return v;
  }
  else
  { /* model not minimal */
    GEN pk = powiu(p,k), p2k = sqri(pk), p4k = sqri(p2k), p6k = mulii(p4k,p2k);
    GEN r, s, t;

    s = negi(gel(e,1));
    if (mpodd(s)) s = addii(s, pk);
    s = shifti(s, -1);

    r = subii(gel(e,2), mulii(s, addii(gel(e,1), s))); /* a_2' */
    switch(umodiu(r, 3))
    {
      default: break; /* 0 */
      case 2: r = addii(r, p2k); break;
      case 1: r = subii(r, p2k); break;
    }
    r = negi( diviuexact(r, 3) );

    t = negi(ellLHS0_i(e,r)); /* - a_3' */
    if (mpodd(t)) t = addii(t, mulii(pk, p2k));
    t = shifti(t, -1);

    gel(v,1) = pk;
    gel(v,2) = r;
    gel(v,3) = s;
    gel(v,4) = t;
    if (minim) return v;

    nuD -= 12 * k;
    c4 = diviiexact(c4, p4k);
    c6 = diviiexact(c6, p6k);
    D = diviiexact(D, sqri(p6k));
  }

  if (nuj > 0) switch(nuD - nuj)
  {
    case 0: f = 1; kod = 4+nuj; /* Inu */
      switch(kronecker(negi(c6),p))
      {
	case  1: c = nuD; break;
	case -1: c = odd(nuD)? 1: 2; break;
	default: return localredbug(p,"localred (p | c6)");
      }
      break;
    case 6: f = 2; kod = -4-nuj; /* Inu* */
      if (nuj & 1)
	c = 3 + kronecker(diviiexact(mulii(c6, D),powiu(p, 9+nuj)), p);
      else
	c = 3 + kronecker(diviiexact(D, powiu(p, 6+nuj)), p);
      break;
    default: return localredbug(p,"localred (nu_D - nu_j != 0,6)");
  }
  else switch(nuD)
  {
    case  0: f = 0; kod = 1; c = 1; break; /* I0, regular */
    case  2: f = 2; kod = 2; c = 1; break; /* II   */
    case  3: f = 2; kod = 3; c = 2; break; /* III  */
    case  4: f = 2; kod = 4; /* IV   */
      c = 2 + kronecker(mulis(diviiexact(c6,sqri(p)), -6), p);
      break;
    case  6: f = 2; kod = -1; /* I0*  */
      p2 = sqri(p);
      /* x^3 - 3c4/p^2 x - 2c6/p^3 */
      tri = mkpoln(4, gen_1, gen_0,
                            negi(mulsi(3, diviiexact(c4, p2))),
                            negi(shifti(diviiexact(c6, mulii(p2,p)), 1)));
      c = 1 + FpX_nbroots(tri, p);
      break;
    case  8: f = 2; kod = -4; /* IV*  */
      c = 2 + kronecker(mulsi(-6, diviiexact(c6, sqri(sqri(p)))), p);
      break;
    case  9: f = 2; kod = -3; c = 2; break; /* III* */
    case 10: f = 2; kod = -2; c = 1; break; /* II*  */
    default: return localredbug(p,"localred");
  }
  return localred_result(f, kod, c, v);
}

/* return a_{ k,l } in Tate's notation, pl = p^l */
static ulong
aux(GEN ak, ulong q, ulong pl)
{
  return umodiu(ak, q) / pl;
}

static ulong
aux2(GEN ak, ulong p, GEN pl)
{
  pari_sp av = avma;
  ulong res = umodiu(diviiexact(ak, pl), p);
  avma = av; return res;
}

/* number of distinct roots of X^3 + aX^2 + bX + c modulo p = 2 or 3
 * assume a,b,c in {0, 1} [ p = 2] or {0, 1, 2} [ p = 3 ]
 * if there's a multiple root, put it in *mult */
static long
numroots3(long a, long b, long c, long p, long *mult)
{
  if (p == 2)
  {
    if ((c + a * b) & 1) return 3;
    *mult = b; return (a + b) & 1 ? 2 : 1;
  }
  /* p = 3 */
  if (!a) { *mult = -c; return b ? 3 : 1; }
  *mult = a * b;
  if (b == 2)
    return (a + c) == 3 ? 2 : 3;
  else
    return c ? 3 : 2;
}

/* same for aX^2 +bX + c */
static long
numroots2(long a, long b, long c, long p, long *mult)
{
  if (p == 2) { *mult = c; return b & 1 ? 2 : 1; }
  /* p = 3 */
  *mult = a * b; return (b * b - a * c) % 3 ? 2 : 1;
}

/* p = 2 or 3 */
static GEN
localred_23(GEN e, long p)
{
  long c, nu, nuD, r, s, t;
  long theroot, p2, p3, p4, p5, p6, a21, a42, a63, a32, a64;
  GEN v;

  nuD = Z_lval(gel(e,12), (ulong)p);
  v = init_ch();
  if (p == 2) { p2 = 4; p3 = 8;  p4 = 16; p5 = 32; p6 = 64;}
  else        { p2 = 9; p3 = 27; p4 = 81; p5 =243; p6 =729; }

  for (;;)
  {
    if (!nuD) return localred_result(0, 1, 1, v);
        /* I0   */
    if (umodiu(gel(e,6), p)) /* p \nmid b2 */
    {
      if (umodiu(negi(gel(e,11)), p == 2 ? 8 : 3) == 1)
        c = nuD;
      else
        c = 2 - (nuD & 1);
      return localred_result(1, 4 + nuD, c, v);
    }
        /* Inu  */
    if (p == 2)
    {
      r = umodiu(gel(e,4), 2);
      s = umodiu(gel(e,2), 2);
      t = umodiu(gel(e,5), 2);
      if (r) { t = (s + t) & 1; s = (s + 1) & 1; }
    }
    else /* p == 3 */
    {
      r = - umodiu(gel(e,8), 3);
      s = umodiu(gel(e,1), 3);
      t = umodiu(gel(e,3), 3);
      if (s) { t  = (t + r*s) % 3; if (t < 0) t += 3; }
    }
    /* p | (a1, a2, a3, a4, a6) */
    if (r || s || t) cumule(&v, &e, gen_1, stoi(r), stoi(s), stoi(t));
    if (umodiu(gel(e,5), p2))
      return localred_result(nuD, 2, 1, v);
        /* II   */
    if (umodiu(gel(e,9), p3))
      return localred_result(nuD - 1, 3, 2, v);
        /* III  */
    if (umodiu(gel(e,8), p3))
    {
      if (umodiu(gel(e,8), (p==2)? 32: 27) == (ulong)p2)
        c = 3;
      else
        c = 1;
      return localred_result(nuD - 2, 4, c, v);
    }
        /* IV   */

    if (umodiu(gel(e,5), p3))
      cumule(&v, &e, gen_1, gen_0, gen_0, p == 2? gen_2: modis(gel(e,3), 9));
        /* p | a1, a2; p^2  | a3, a4; p^3 | a6 */
    a21 = aux(gel(e,2), p2, p);
    a42 = aux(gel(e,4), p3, p2);
    a63 = aux(gel(e,5), p4, p3);
    switch (numroots3(a21, a42, a63, p, &theroot))
    {
      case 3:
        c = a63 ? 1: 2;
        if (p == 2)
          c += ((a21 + a42 + a63) & 1);
        else {
          if (((1 + a21 + a42 + a63) % 3) == 0) c++;
          if (((1 - a21 + a42 - a63) % 3) == 0) c++;
        }
        return localred_result(nuD - 4, -1, c, v);
      case 2: /* I0*  */
      { /* compute nu */
        GEN pk, pk1, p2k;
        long al, be, ga;
        if (theroot) cumule(&v, &e, gen_1, stoi(theroot * p), gen_0, gen_0);
            /* p | a1; p^2  | a2, a3; p^3 | a4; p^4 | a6 */
        nu = 1;
        pk  = utoipos(p2);
        p2k = utoipos(p4);
        for(;;)
        {
          be =  aux2(gel(e,3), p, pk);
          ga = -aux2(gel(e,5), p, p2k);
          al = 1;
          if (numroots2(al, be, ga, p, &theroot) == 2) break;
          if (theroot) cumule(&v, &e, gen_1, gen_0, gen_0, mulsi(theroot,pk));
          pk1 = pk;
          pk  = mului(p, pk);
          p2k = mului(p, p2k); nu++;

          al = a21;
          be = aux2(gel(e,4), p, pk);
          ga = aux2(gel(e,5), p, p2k);
          if (numroots2(al, be, ga, p, &theroot) == 2) break;
          if (theroot) cumule(&v, &e, gen_1, mulsi(theroot, pk1), gen_0, gen_0);
          p2k = mulsi(p, p2k); nu++;
        }
        if (p == 2)
          c = 4 - 2 * (ga & 1);
        else
          c = 3 + kross(be * be - al * ga, 3);
        return localred_result(nuD - 4 - nu, -4 - nu, c, v);
      }
      case 1: /* Inu* */
        if (theroot) cumule(&v, &e, gen_1, stoi(theroot*p), gen_0, gen_0);
            /* p | a1; p^2  | a2, a3; p^3 | a4; p^4 | a6 */
        a32 = aux(gel(e,3), p3, p2);
        a64 = aux(gel(e,5), p5, p4);
        if (numroots2(1, a32, -a64, p, &theroot) == 2)
        {
          if (p == 2)
            c = 3 - 2 * a64;
          else
            c = 2 + kross(a32 * a32 + a64, 3);
          return localred_result(nuD - 6, -4, c, v);
        }
            /* IV*  */
        if (theroot) cumule(&v, &e, gen_1, gen_0, gen_0, stoi(theroot*p2));
            /* p | a1; p^2 | a2; p^3 | a3, a4; p^5 | a6 */
        if (umodiu(gel(e,4), p4))
          return localred_result(nuD - 7, -3, 2, v);
            /* III* */

        if (umodiu(gel(e,5), p6))
          return localred_result(nuD - 8, -2, 1, v);
            /* II*  */
        cumule(&v, &e, utoipos(p), gen_0, gen_0, gen_0); /* not minimal */
        nuD -= 12;
    }
  }
}

static GEN
localred(GEN e, GEN p, int minim)
{
  if (cmpiu(p, 3) > 0) /* p != 2,3 */
    return localred_p(e,p, minim);
  else
  {
    long l = itos(p);
    GEN z;
    if (l < 2) pari_err(talker,"not a prime in localred");
    z = localred_23(e, l);
    return minim? gel(z,3): z;
  }
}

GEN
elllocalred(GEN e, GEN p)
{
  pari_sp av = avma;
  checkell(e);
  if (typ(e[12]) != t_INT)
    pari_err(talker,"not an integral curve in elllocalred");
  if (typ(p) != t_INT || signe(p) <= 0) pari_err(typeer,"elllocalred");
  return gerepileupto(av, localred(e, p, 0));
}

static GEN
ellintegralmodel(GEN e)
{
  GEN a = cgetg(6,t_VEC), v, L, u;
  long i, l, k;

  checkell(e);
  L = cgetg(1, t_VEC);
  for (i = 1; i < 6; i++)
  {
    a[i] = e[i]; u = gel(a,i);
    switch(typ(u))
    {
      case t_INT: break;
      case t_FRAC: /* partial factorization */
        L = shallowconcat(L, (GEN)auxdecomp(gel(u,2), 0)[1]);
        break;
      default: pari_err(talker, "not a rational curve in ellintegralmodel");
    }
  }
  /* a = [a1, a2, a3, a4, a6] */
  l = lg(L);
  if (l == 1) return NULL;
  L = sort(L);
  for (k = i = 2; i < l; i++)
    if (!equalii(gel(L,i), gel(L,i-1))) L[k++] = L[i];

  l = k; u = gen_1;
  for (k = 1; k < l; k++)
  {
    GEN p = gel(L,k);
    long n = 0, m;
    for (i = 1; i < 6; i++)
      if (!gcmp0(gel(a,i)))
      {
        long r = (i == 5)? 6: i; /* a5 is missing */
	m = r * n + ggval(gel(a,i), p);
	while (m < 0) { n++; m += r; }
      }
    u = mulii(u, powiu(p, n));
  }
  v = init_ch(); gel(v,1) = ginv(u); return v;
}

/* e integral model */
static void
standard_model(GEN e, GEN *pv)
{
  GEN a1 = gel(e,1), a2 = gel(e,2);
  GEN r, t, s = truedivis(a1, -2);
  r = truedivis(addis(subii(a2, mulii(s,addii(s,a1))), 1), -3);
  t = truedivis(ellLHS0_i(e,r), -2);
  cumulev(pv, gen_1, r, s, t);
}

GEN
ellminimalmodel(GEN E, GEN *ptv)
{
  pari_sp av = avma;
  GEN c4, c6, e, v, v0, P;
  long l, k;

  v0 = ellintegralmodel(E);
  e = ell_to_small(E);
  if (v0) e = _coordch(e, v0);
  v = init_ch();
  c4 = gel(e,10);
  c6 = gel(e,11);
  P = (GEN)Z_factor(gcdii(c4,c6))[1];
  l = lg(P);
  for (k = 1; k < l; k++)
  {
    GEN w = localred(e, gel(P,k), 1);
    if (!gcmp1(gel(w,1)))
      cumule(&v, &e, gel(w,1), gel(w,2), gel(w,3), gel(w,4));
  }
  standard_model(e, &v);
  if (v0) { gcumulev(&v0, v); v = v0; }
  e = _coordch(E, v);
  if (ptv) { gerepileall(av, 2, &e, &v); *ptv = v; }
  else e = gerepilecopy(av, e);
  return e;
}

/* Reduction of a rational curve E to its standard minimal model
 * (a1,a3 = 0 or 1, a2 = -1, 0 or 1).
 *
 * Return [N, [u,r,s,t], c], where
 *   N = arithmetic conductor of E
 *   c = product of the local Tamagawa numbers cp
 *   [u, r, s, t] = the change of variable reducing E to its minimal model,
 *     with u > 0 */
GEN
ellglobalred(GEN E)
{
  long k, l;
  pari_sp av = avma;
  GEN c, P, N, v, v0, e, c4, c6, D;

  v0 = ellintegralmodel(E);
  e = ell_to_small(E);
  if (v0) e = _coordch(e, v0);
  v = init_ch();
  c4 = gel(e,10);
  c6 = gel(e,11);
  D  = gel(e,12);
  P = (GEN)Z_factor(gcdii(c4,c6))[1];
  l = lg(P);
  for (k = 1; k < l; k++) (long)Z_pvalrem(D, gel(P,k), &D);
  if (!is_pm1(D)) P = shallowconcat(P, (GEN)Z_factor(absi(D))[1]);
  l = lg(P); c = N = gen_1;
  for (k = 1; k < l; k++)
  {
    GEN p = gel(P,k), q = localred(e, p, 0), w = gel(q,3);
    N = mulii(N, powgi(p, gel(q,1)));
    c = mulii(c, gel(q,4));
    if (!gcmp1(gel(w,1)))
      cumule(&v, &e, gel(w,1), gel(w,2), gel(w,3), gel(w,4));
  }
  standard_model(e, &v);
  if (v0) { gcumulev(&v0, v); v = v0; }
  return gerepilecopy(av, mkvec3(N,v,c));
}

/********************************************************************/
/**                                                                **/
/**           ROOT NUMBER (after Halberstadt at p = 2,3)           **/
/**                                                                **/
/********************************************************************/

/* p = 2 or 3 */
static long
neron(GEN e, long p, long* ptkod)
{
  long kod, v4, v6, vd;
  pari_sp av=avma;
  GEN c4, c6, d, nv;

  nv = localred_23(e,p);
  *ptkod = kod = itos(gel(nv,2));
  c4=gel(e,10); c6=gel(e,11); d=gel(e,12);
  v4 = gcmp0(c4) ? 12 : Z_lval(c4,p);
  v6 = gcmp0(c6) ? 12 : Z_lval(c6,p);
  vd = Z_lval(d,p); avma = av;
  if (p == 2) {
    if (kod > 4) return 1;
    switch(kod)
    {
      case 1: return (v6>0) ? 2 : 1;
      case 2:
        if (vd==4) return 1;
        else
        {
          if (vd==7) return 3;
          else return v4==4 ? 2 : 4;
        }
      case 3:
        switch(vd)
        {
          case 6: return 3;
          case 8: return 4;
          case 9: return 5;
          default: return v4==5 ? 2 : 1;
        }
      case 4: return v4>4 ? 2 : 1;
      case -1:
        switch(vd)
        {
          case 9: return 2;
          case 10: return 4;
          default: return v4>4 ? 3 : 1;
        }
      case -2:
        switch(vd)
        {
          case 12: return 2;
          case 14: return 3;
          default: return 1;
        }
      case -3:
        switch(vd)
        {
          case 12: return 2;
          case 14: return 3;
          case 15: return 4;
          default: return 1;
        }
      case -4: return v6==7 ? 2 : 1;
      case -5: return (v6==7 || v4==6) ? 2 : 1;
      case -6:
        switch(vd)
        {
          case 12: return 2;
          case 13: return 3;
          default: return v4==6 ? 2 : 1;
        }
      case -7: return (vd==12 || v4==6) ? 2 : 1;
      default: return v4==6 ? 2 : 1;
    }
  } else {
    if (labs(kod) > 4) return 1;
    switch(kod)
    {
      case -1: case 1: return v4&1 ? 2 : 1;
      case -3: case 3: return (2*v6>vd+3) ? 2 : 1;
      case -4: case 2:
        switch (vd%6)
        {
          case 4: return 3;
          case 5: return 4;
          default: return v6%3==1 ? 2 : 1;
        }
      default: /* kod = -2 et 4 */
        switch (vd%6)
        {
          case 0: return 2;
          case 1: return 3;
          default: return 1;
        }
    }
  }
}

static long
val_aux(GEN x, long p, long pk, long *u) {
  long v;
  GEN z;
  if (!signe(x)) { *u = 0; return 12; }
  v = Z_lvalrem(x,p,&z);
  *u = umodiu(z,pk); return v;
}
static void
val_init(GEN e, long p, long pk,
         long *v4, long *u, long *v6, long *v, long *vd, long *d1)
{
  GEN c4 = gel(e,10), c6 = gel(e,11), D = gel(e,12);
  pari_sp av = avma;
  *v4 = val_aux(c4, p,pk, u);
  *v6 = val_aux(c6, p,pk, v);
  *vd = val_aux(D , p,pk, d1); avma = av;
}

static long
ellrootno_2(GEN e)
{
  long n2, kod, u, v, x1, y1, d1, vd, v4, v6;

  val_init(e, 2, 64, &v4, &u, &v6, &v, &vd, &d1);
  if (!vd) return 1;
  n2 = neron(e,2,&kod);
  if (kod>=5)
    return odd(umodiu(gel(e,2),2) + umodiu(gel(e,3),2)) ? 1 : -1;
  if (kod<-9) return (n2==2) ? -kross(-1,v) : -1;
  x1 = u+v+v;
  switch(kod)
  {
    case 1: return 1;
    case 2:
      switch(n2)
      {
	case 1:
	  switch(v4)
	  {
	    case 4: return kross(-1,u);
	    case 5: return 1;
	    default: return -1;
	  }
	case 2: return (v6==7) ? 1 : -1;
	case 3: return (v%8==5 || (u*v)%8==5) ? 1 : -1;
	case 4: if (v4>5) return kross(-1,v);
	  return (v4==5) ? -kross(-1,u) : -1;
      }
    case 3:
      switch(n2)
      {
	case 1: return -kross(2,u*v);
	case 2: return -kross(2,v);
	case 3: y1 = (u - (v << (v6-5))) & 15;
	  return (y1==7 || y1==11) ? 1 : -1;
	case 4: return (v%8==3 || (2*u+v)%8==7) ? 1 : -1;
	case 5: return v6==8 ? kross(2,x1) : kross(-2,u);
      }
    case -1:
      switch(n2)
      {
	case 1: return -kross(2,x1);
	case 2: return (v%8==7) || (x1%32==11) ? 1 : -1;
	case 3: return v4==6 ? 1 : -1;
	case 4: if (v4>6) return kross(-1,v);
	  return v4==6 ? -kross(-1,u*v) : -1;
      }
    case -2: return n2==1 ? kross(-2,v) : kross(-1,v);
    case -3:
      switch(n2)
      {
	case 1: y1=(u-2*v)%64; if (y1<0) y1+=64;
	  return (y1==3) || (y1==19) ? 1 : -1;
	case 2: return kross(2*kross(-1,u),v);
	case 3: return -kross(-1,u)*kross(-2*kross(-1,u),u*v);
	case 4: return v6==11 ? kross(-2,x1) : -kross(-2,u);
      }
    case -5:
      if (n2==1) return x1%32==23 ? 1 : -1;
      else return -kross(2,2*u+v);
    case -6:
      switch(n2)
      {
	case 1: return 1;
	case 2: return v6==10 ? 1 : -1;
	case 3: return (u%16==11) || ((u+4*v)%16==3) ? 1 : -1;
      }
    case -7:
      if (n2==1) return 1;
      else
      {
	y1 = (u + (v << (v6-8))) & 15;
	if (v6==10) return (y1==9 || y1==13) ? 1 : -1;
	else return (y1==9 || y1==5) ? 1 : -1;
      }
    case -8: return n2==2 ? kross(-1,v*d1) : -1;
    case -9: return n2==2 ? -kross(-1,d1) : -1;
    default: return -1;
  }
}

static long
ellrootno_3(GEN e)
{
  long n2, kod, u, v, d1, r6, K4, K6, vd, v4, v6;

  val_init(e, 3, 81, &v4, &u, &v6, &v, &vd, &d1);
  if (!vd) return 1;
  n2 = neron(e,3,&kod);
  K6 = kross(v,3); if (kod>4) return K6;
  r6 = v%9; K4 = kross(u,3);
  switch(kod)
  {
    case 1: case 3: case -3: return 1;
    case 2:
      switch(n2)
      {
	case 1: return (r6==4 || r6>6) ? 1 : -1;
	case 2: return -K4*K6;
	case 3: return 1;
	case 4: return -K6;
      }
    case 4:
      switch(n2)
      {
	case 1: return K6*kross(d1,3);
	case 2: return -K4;
	case 3: return -K6;
      }
    case -2: return n2==2 ? 1 : K6;
    case -4:
      switch(n2)
      {
	case 1:
	  if (v4==4) return (r6==4 || r6==8) ? 1 : -1;
	  else return (r6==1 || r6==2) ? 1 : -1;
	case 2: return -K6;
	case 3: return (r6==2 || r6==7) ? 1 : -1;
	case 4: return K6;
      }
    default: return -1;
  }
}

/* assume p > 3, p^ex || N(E) */
static long
ellrootno_p(GEN e, GEN p, ulong ex)
{
  GEN j;
  long ep,z;

  if (!ex) return 1;
  if (ex == 1) return -kronecker(negi(gel(e,11)),p);
  j=gel(e,13);
  if (!gcmp0(j) && ggval(j,p) < 0) return krosi(-1,p);
  ep = 12 / cgcd(12, Z_pval(gel(e,12),p));
  if (ep==4) z = 2; else z = (ep&1) ? 3 : 1;
  return krosi(-z, p);
}

static long
ellrootno_global(GEN e, GEN N)
{
  long i, v, s = -1;
  GEN fa, P, E;
  
  v = Z_lvalrem(N, 2, &N); if (v) s *= ellrootno_2(e);
  v = Z_lvalrem(N, 3, &N); if (v) s *= ellrootno_3(e);
  fa = Z_factor(N);
  P = gel(fa,1);
  E = gel(fa,2);
  for (i=1; i<lg(P); i++) s *= ellrootno_p(e, gel(P,i), itou(gel(E,i)));
  return s;
}

/* local epsilon factor at p (over Q), including p=0 for the infinite place.
 * Global if p==1 or NULL. */
long
ellrootno(GEN e, GEN p)
{
  pari_sp av = avma;
  long s;
  GEN gr, N;
  checkell(e);
  e = ell_to_small(e); gr = ellglobalred(e);
  e = _coordch(e,gel(gr,2)); N = gel(gr,1);
  if (!p || gcmp1(p))
    s = ellrootno_global(e, N);
  else
  {
    if (typ(p) != t_INT || signe(p) < 0) pari_err(typeer,"ellrootno");
    if (cmpiu(p,3) > 0) s = ellrootno_p(e,p, Z_pval(N, p));
    else switch(itou(p))
    {
      case 2: s = ellrootno_2(e); break;
      case 3: s = ellrootno_3(e); break;
      default: s = -1; break; /* local factor at infinity */
    }
  }
  avma = av; return s;
}

/********************************************************************/
/**                                                                **/
/**                       TRACE OF FROBENIUS                       **/
/**                                                                **/
/********************************************************************/

/* compute a_2 */
static GEN
a2(GEN e)
{ /* solve y(1 + a1x + a3) = x (1 + a2 + a4) + a6 */
  pari_sp av = avma;
  ulong a1 = Rg_to_Fl(gel(e,1), 2);
  ulong a2 = Rg_to_Fl(gel(e,2), 2);
  ulong a3 = Rg_to_Fl(gel(e,3), 2);
  ulong a4 = Rg_to_Fl(gel(e,4), 2);
  ulong a6 = Rg_to_Fl(gel(e,5), 2);
  long N = 1; /* oo */
  if (!a3) N ++; /* x = 0, y=0 or 1 */
  else if (!a6) N += 2; /* x = 0, y arbitrary */
  if ((a3 ^ a1) == 0) N++; /* x = 1, y = 0 or 1 */
  else if (a2 ^ a4 ^ a6) N += 2; /* x = 1, y arbitrary */
  avma = av; return stoi(3 - N);
}
/* a_p using Jacobi sums */
static GEN
ap_jacobi(GEN e, ulong p)
{
  if (p == 2) return a2(e);
  else
  {
    ulong i;
    ulong e6 = Rg_to_Fl(gel(e,6), p);
    ulong e8 = Rg_to_Fl(gel(e,8), p);
    ulong e72= Rg_to_Fl(gel(e,7), p) << 1;
    long s = krouu(e8, p) + krouu((e8 + e72 + e6 + 4) % p, p); /* i = 0,1 */
    if (p < 757UL)
      for (i=2; i<p; i++)
        s += krouu((e8 + i*(e72 + i*(e6 + (i<<2)))) % p, p);
    else
      for (i=2; i<p; i++)
        s += krouu(e8 + Fl_mul(i, e72 + Fl_mul(i, e6 + (i<<2), p), p), p);
    return stoi(-s);
  }
}

GEN
apell2(GEN e, GEN pp)
{
  checkell(e); if (typ(pp)!=t_INT) pari_err(elliper1);
  if (expi(pp) > 29) pari_err(talker,"prime too large in apell2, use apell");
  return ap_jacobi(e, (ulong)pp[2]);
}

/* invert all elements of x mod p using Montgomery's trick */
GEN
multi_invmod(GEN x, GEN p)
{
  long i, lx = lg(x);
  GEN u,y = cgetg(lx, t_VEC);

  y[1] = x[1];
  for (i=2; i<lx; i++)
    gel(y,i) = remii(mulii(gel(y,i-1), gel(x,i)), p);

  u = Fp_inv(gel(y,--i), p);
  for ( ; i > 1; i--)
  {
    gel(y,i) = remii(mulii(u, gel(y,i-1)), p);
    u = remii(mulii(u, gel(x,i)), p); /* u = 1 / (x[1] ... x[i-1]) */
  }
  gel(y,1) = u; return y;
}

static GEN
addsell(GEN e, GEN z1, GEN z2, GEN p)
{
  GEN z,p1,p2,x,x1,x2,y,y1,y2;
  pari_sp av;

  if (!z1) return z2;
  if (!z2) return z1;

  x1 = gel(z1,1); y1 = gel(z1,2);
  x2 = gel(z2,1); y2 = gel(z2,2);
  z = cgetg(3, t_VEC); av = avma;
  if (x1 == x2 || equalii(x1, x2))
  {
    if (!signe(y1) || !equalii(y1,y2)) return NULL;
    p2 = shifti(y1,1);
    p1 = addii(e, mulii(x1,mulsi(3,x1)));
    p1 = remii(p1, p);
  }
  else { p1 = subii(y2,y1); p2 = subii(x2, x1); }
  p1 = mulii(p1, Fp_inv(p2, p));
  p1 = remii(p1, p);
  x = subii(sqri(p1), addii(x1,x2));
  y = negi(addii(y1, mulii(p1,subii(x,x1))));
  x = modii(x,p);
  y = modii(y,p); avma = av;
  gel(z,1) = icopy(x);
  gel(z,2) = icopy(y); return z;
}

/* z1 <-- z1 + z2 */
static void
addsell_part2(GEN e, GEN z1, GEN z2, GEN p, GEN p2inv)
{
  GEN p1,x,x1,x2,y,y1,y2;

  x1 = gel(z1,1); y1 = gel(z1,2);
  x2 = gel(z2,1); y2 = gel(z2,2);
  if (x1 == x2)
  {
    p1 = addii(e, mulii(x1,mulsi(3,x1)));
    p1 = remii(p1, p);
  }
  else p1 = subii(y2,y1);

  p1 = mulii(p1, p2inv);
  p1 = remii(p1, p);
  x = subii(sqri(p1), addii(x1,x2)); x = modii(x,p);
  y = negi(addii(y1, mulii(p1,subii(x,x1)))); y = modii(y,p);
  affii(x, x1);
  affii(y, y1);
}

static GEN
negsell(GEN f, GEN p)
{
  GEN g = cgetg(3, t_VEC), y = gel(f,2);
  gel(g,1) = gel(f,1);
  gel(g,2) = signe(y)? subii(p, y): y;
  return g;
}

typedef struct {
  GEN e, p;
} sellp;

static GEN
mul_sell(void *d, GEN x, GEN y)
{
  sellp *S = (sellp*)d;
  return addsell(S->e, x, y, S->p);
}
static GEN
sqr_sell(void *d, GEN x)
{
  sellp *S = (sellp*)d;
  return addsell(S->e, x, x, S->p);
}

static GEN
powsell(GEN e, GEN z, GEN n, GEN p)
{
  long s = signe(n);
  sellp S;

  if (!s || !z) return NULL;
  if (s < 0) z = negsell(z, p);
  if (is_pm1(n)) return z;
  S.e = e;
  S.p = p;
  return leftright_pow(z, n, &S, &sqr_sell, &mul_sell);
}

/* assume H.f = 0, return exact order of f */
static GEN
exact_order(GEN H, GEN f, GEN c4, GEN p)
{
  GEN P, e, h = H, fa = Z_factor(H);
  long i, j, l;

  P = gel(fa,1); l = lg(P);
  e = gel(fa,2);
  for (i=1; i<l; i++)
    for (j=itos(gel(e,i)); j; j--)
    {
      GEN n = diviiexact(h,gel(P,i));
      if (powsell(c4,f,n,p)) break;
      h = n;
    }
  return h;
}

/* make sure *x has lgefint >= k */
static void
_fix(GEN x, long k)
{
  GEN y = (GEN)*x;
  if (lgefint(y) < k) { GEN p1 = cgeti(k); affii(y,p1); *x = (long)p1; }
}

INLINE long safemodBIL(GEN x) { return signe(x)?modBIL(x):0; }

/* Return the lift of a (mod b), which is closest to h */
static GEN
closest_lift(GEN a, GEN b, GEN h)
{
  return addii(a, mulii(b, diviiround(gsub(h,a), b)));
}

/* compute a_p using Shanks/Mestre + Montgomery's trick. Assume p > 457 */
GEN
apell1(GEN e, GEN p)
{
  long *tx, *ty, *ti, pfinal, i, j, s, KRO, KROold, nb;
  ulong x;
  pari_sp av = avma, av2;
  GEN p1, h, mfh, F, f, fh, fg, pordmin, u, v, p1p, p2p, A, B, c4, c6, cp4, pts;
  tx = NULL;
  ty = ti = NULL; /* gcc -Wall */

  if (DEBUGLEVEL) (void)timer2();
  c4 = Rg_to_Fp(gdivgs(gel(e,10),  -48), p);
  c6 = Rg_to_Fp(gdivgs(gel(e,11), -864), p);
  /* once #E(Fp) is know mod B >= pordmin, it is completely determined */
  pordmin = addis(sqrti(gmul2n(p,4)), 1); /* ceil( 4sqrt(p) ) */
  p1p = addsi(1, p);
  p2p = shifti(p1p, 1);
  x = 0; u = c6; KRO = kronecker(u, p); KROold = - KRO;
  A = gen_0; B = gen_1; h = p1p;
  for(;;)
  {
    long CODE;
    while (!KRO || KRO == KROold)
    { /* look for points alternatively on E and its quadratic twist E' */
      x++; /* u = x^3 + c4 x + c6 */
      u = modii(addii(c6, mului(x, addii(c4, muluu(x,x)))), p);
      KRO = kronecker(u, p);
    }
    KROold = KRO;
    /* [ux, u^2] is on E_u: y^2 = x^3 + c4 u^2 x + c6 u^3
     * E_u isomorphic to E (resp. E') iff KRO = 1 (resp. -1)
     * #E(F_p) = p+1 - a_p, #E'(F_p) = p+1 + a_p
     *
     * #E_u(Fp) = A (mod B),  h is close to #E_u(Fp) */

    f = cgetg(3,t_VEC);
    gel(f,1) = modii(mului(x,u), p);
    gel(f,2) = modii(sqri(u),    p);
    cp4 = modii(mulii(c4, gel(f,2)), p); /* c4 for E_u */
    fh = powsell(cp4,f,h,p);
    if (!fh) goto FOUND;

    s = itos( gceil(gsqrt(gdiv(pordmin,B),DEFAULTPREC)) ) >> 1;
    CODE = evaltyp(t_VECSMALL) | evallg(s+1);
    /* look for h s.t f^h = 0 */
    if (!tx)
    { /* first time: initialize */
      tx = newbloc(3*(s+1));
      ty = tx + (s+1);
      ti = ty + (s+1);
    }
    F = powsell(cp4,f,B,p);
    *tx = CODE;

    /* F = B.f */
    p1 = gcopy(fh);
    if (s < 3)
    { /* we're nearly done: naive search */
      GEN q1 = p1, mF = negsell(F, p); /* -F */
      for (i=1;; i++)
      {
        p1 = addsell(cp4,p1, F,p); /* h.f + i.F */
        if (!p1) { h = addii(h, mulsi( i,B)); goto FOUND; }
        q1 = addsell(cp4,q1,mF,p); /* h.f - i.F */
        if (!q1) { h = addii(h, mulsi(-i,B)); goto FOUND; }
      }
    }
    /* Baby Step/Giant Step */
    nb = min(128, s >> 1); /* > 0. Will do nb pts at a time: faster inverse */
    pts = cgetg(nb+1, t_VEC);
    j = lgefint(p);
    for (i=1; i<=nb; i++)
    { /* baby steps */
      gel(pts,i) = p1; /* h.f + (i-1).F */
      _fix(p1+1, j); tx[i] = safemodBIL(gel(p1,1));
      _fix(p1+2, j); ty[i] = safemodBIL(gel(p1,2));
      p1 = addsell(cp4,p1,F,p); /* h.f + i.F */
      if (!p1) { h = addii(h, mulsi(i,B)); goto FOUND; }
    }
    mfh = negsell(fh, p);
    fg = addsell(cp4,p1,mfh,p); /* h.f + nb.F - h.f = nb.F */
    if (!fg) { h = mulsi(nb,B); goto FOUND; }
    u = cgetg(nb+1, t_VEC);
    av2 = avma; /* more baby steps, nb points at a time */
    while (i <= s)
    {
      long maxj;
      for (j=1; j<=nb; j++) /* adding nb.F (part 1) */
      {
        p1 = gel(pts,j); /* h.f + (i-nb-1+j-1).F */
        gel(u,j) = subii(gel(fg,1), gel(p1,1));
        if (u[j] == (long)gen_0) /* sum = 0 or doubling */
        {
          long k = i+j-2;
          if (equalii(gel(p1,2),gel(fg,2))) k -= 2*nb; /* fg == p1 */
          h = addii(h, mulsi(k,B)); goto FOUND;
        }
      }
      v = multi_invmod(u, p);
      maxj = (i-1 + nb <= s)? nb: s % nb;
      for (j=1; j<=maxj; j++,i++) /* adding nb.F (part 2) */
      {
        p1 = gel(pts,j);
        addsell_part2(cp4,p1,fg,p, gel(v,j));
        tx[i] = safemodBIL(gel(p1,1));
        ty[i] = safemodBIL(gel(p1,2));
      }
      avma = av2;
    }
    p1 = addsell(cp4,gel(pts,j-1),mfh,p); /* = (s-1).F */
    if (!p1) { h = mulsi(s-1,B); goto FOUND; }
    if (DEBUGLEVEL) msgtimer("[apell1] baby steps, s = %ld",s);

    /* giant steps: fg = s.F */
    fg = addsell(cp4,p1,F,p);
    if (!fg) { h = mulsi(s,B); goto FOUND; }
    pfinal = safemodBIL(p); av2 = avma;

    p1 = vecsmall_indexsort(tx);
    for (i=1; i<=s; i++) ti[i] = tx[p1[i]];
    for (i=1; i<=s; i++) { tx[i] = ti[i]; ti[i] = ty[p1[i]]; }
    for (i=1; i<=s; i++) { ty[i] = ti[i]; ti[i] = p1[i]; }
    if (DEBUGLEVEL) msgtimer("[apell1] sorting");
    avma = av2;

    gaffect(fg, gel(pts,1));
    for (j=2; j<=nb; j++) /* pts[j] = j.fg = (s*j).F */
    {
      p1 = addsell(cp4,gel(pts,j-1),fg,p);
      if (!p1) { h = mulii(mulss(s,j), B); goto FOUND; }
      gaffect(p1, gel(pts,j));
    }
    /* replace fg by nb.fg since we do nb points at a time */
    avma = av2;
    fg = gcopy(gel(pts,nb));
    av2 = avma;

    for (i=1,j=1; ; i++)
    {
      GEN ftest = gel(pts,j);
      ulong m, l = 1, r = s+1;
      long k, k2, j2;

      avma = av2;
      k = safemodBIL(gel(ftest,1));
      while (l<r)
      {
        m = (l+r) >> 1;
        if (tx[m] < k) l = m+1; else r = m;
      }
      if (r <= (ulong)s && tx[r] == k)
      {
        while (tx[r] == k && r) r--;
        k2 = safemodBIL(gel(ftest,2));
        for (r++; tx[r] == k && r <= (ulong)s; r++)
          if (ty[r] == k2 || ty[r] == pfinal - k2)
          { /* [h+j2] f == +/- ftest (= [i.s] f)? */
            j2 = ti[r] - 1;
            if (DEBUGLEVEL) msgtimer("[apell1] giant steps, i = %ld",i);
            p1 = addsell(cp4, powsell(cp4,F,stoi(j2),p),fh,p);
            if (equalii(gel(p1,1), gel(ftest,1)))
            {
              if (equalii(gel(p1,2), gel(ftest,2))) i = -i;
              h = addii(h, mulii(addis(mulss(s,i), j2), B));
              goto FOUND;
            }
          }
      }
      if (++j > nb)
      { /* compute next nb points */
        long save = 0; /* gcc -Wall */;
        for (j=1; j<=nb; j++)
        {
          p1 = gel(pts,j);
          gel(u,j) = subii(gel(fg,1), gel(p1,1));
          if (u[j] == (long)gen_0) /* occurs once: i = j = nb, p1 == fg */
          {
            gel(u,j) = shifti(gel(p1,2),1);
            save = fg[1]; fg[1] = p1[1];
          }
        }
        v = multi_invmod(u, p);
        for (j=1; j<=nb; j++)
          addsell_part2(cp4, gel(pts,j),fg,p, gel(v,j));
        if (i == nb) { fg[1] = save; }
        j = 1;
      }
    }
FOUND: /* found a point of exponent h on E_u */
    h = exact_order(h, f, cp4, p);
    /* h | #E_u(Fp) = A (mod B) */
    if (B == gen_1) B = h;
    else
    {
      p1 = chinese(mkintmod(A,B), mkintmod(gen_0, h));
      A = gel(p1,2);
      B = gel(p1,1);
    }

    i = (cmpii(B, pordmin) < 0);
    /* If we are not done, update A mod B for the _next_ curve, isomorphic to
     * the quadratic twist of this one */
    if (i) A = remii(subii(p2p,A), B); /* #E(Fp)+#E'(Fp) = 2p+2 */

    /* h = A mod B, closest lift to p+1 */
    h = closest_lift(A, B, p1p);
    if (!i) break;
  }
  if (tx) gunclone(tx);
  return gerepileuptoint(av, KRO==1? subii(p1p,h): subii(h,p1p));
}

typedef struct
{
  int isnull;
  long x,y;
} sellpt;

/* P <-- P + Q, safe with P = Q */
static void
s_addell(sellpt *P, sellpt *Q, long c4, long p)
{
  ulong num, den, lambda;

  if (P->isnull) { *P = *Q; return; }
  if (Q->isnull) return;
  if (P->x == Q->x)
  {
    if (! P->y || P->y != Q->y) { P->isnull = 1; return; }
    num = Fl_add(c4, Fl_mul(3, Fl_mul(P->x, P->x, p), p), p);
    den = Fl_add(P->y, P->y, p);
  }
  else
  {
    num = Fl_sub(P->y, Q->y, p);
    den = Fl_sub(P->x, Q->x, p);
  }
  lambda = Fl_div(num, den, p);
  num = Fl_sub(Fl_mul(lambda, lambda, p), Fl_add(P->x, Q->x, p), p);
  P->y = Fl_sub(Fl_mul(lambda, Fl_sub(Q->x, num, p), p), Q->y, p);
  P->x = num; /* necessary in case P = Q: we need Q->x above */
}

/* Q <-- P^n */
static void
s_powell(sellpt *Q, sellpt *P, long n, long c4, long p)
{
  sellpt R = *P;

  if (n < 0) { n = -n; if (R.y) R.y = p - R.y; }
  Q->isnull = 1;
  Q->x = Q->y = 0; /* -Wall */
  for(;;)
  {
    if (n&1) s_addell(Q, &R, c4, p);
    n >>= 1; if (!n) return;
    s_addell(&R, &R, c4, p);
  }
}

/* assume H.f = 0, return exact order of f, cf. exact_order */
static long
sexact_order(long H, sellpt *f, long c4, long p)
{
  GEN P, e, fa = factoru(H);
  long h = H, pp, i, j, l;
  sellpt fh;

  P = gel(fa,1); l = lg(P);
  e = gel(fa,2);
  for (i=1; i<l; i++)
  {
    pp = P[i];
    for (j=e[i]; j; j--)
    {
      long n = h / pp;
      s_powell(&fh, f, n, c4, p);
      if (!fh.isnull) break;
      h = n;
    }
  }
  return h;
}

typedef struct
{
  long x,y,i;
} multiple;

static int
compare_multiples(multiple *a, multiple *b) { return a->x - b->x; }

/* assume e has good reduction at p. Should use Montgomery.
 * See apell1() */
static GEN
apell0(GEN e, ulong p)
{
  sellpt f, fh, fg, ftest, F;
  ulong x, u, c4, c6, cp4, p1p, p2p, h;
  long pordmin,A,B;
  long i, s, KRO, KROold, l, r, m;
  pari_sp av;
  multiple *table;

  if (p < 99) return ap_jacobi(e,(ulong)p);
  table = NULL;

  av = avma;
  c4 = Rg_to_Fl(gdivgs(gel(e,10),  -48), p);
  c6 = Rg_to_Fl(gdivgs(gel(e,11), -864), p);
  pordmin = (long)(1 + 4*sqrt((float)p));
  p1p = p+1;
  p2p = p1p << 1;
  x = 0; u = c6; KRO = kross(u, p); KROold = -KRO;
  A = 0; B = 1; h = p1p;
  for(;;)
  {
    while (!KRO || KRO == KROold)
    {
      ulong t;
      if (++x >= p) pari_err(talker, "%lu is not prime, use ellak", p);
      t = Fl_add(c4, Fl_mul(x,x,p), p);
      u = Fl_add(c6, Fl_mul(x, t, p), p);
      KRO = kross(u,p);
    }
    KROold = KRO;
    f.isnull = 0;
    f.x = Fl_mul(x, u, p);
    f.y = Fl_mul(u, u, p);
    cp4 = Fl_mul(c4, f.y, p);
    s_powell(&fh, &f, h, cp4, p);
    s = (long) (sqrt(((float)pordmin)/B) / 2);
    if (!s) s = 1;
    if (!table)
    {
      table = (multiple *) gpmalloc((s+1) * sizeof(multiple));
      F = f;
    }
    else
      s_powell(&F, &f, B, cp4, p);
    for (i=0; i < s; i++)
    {
      if (fh.isnull) { h += B*i; goto FOUND; }
      table[i].x = fh.x;
      table[i].y = fh.y;
      table[i].i = i;
      s_addell(&fh, &F, cp4, p);
    }
    qsort(table,s,sizeof(multiple),(QSCOMP)compare_multiples);
    s_powell(&fg, &F, s, cp4, p); ftest = fg;
    for (i=1; ; i++)
    {
      if (ftest.isnull) {
        if (!uisprime(p)) pari_err(talker,"%lu is not prime, use ellak", p);
        pari_err(bugparier,"apell (f^(i*s) = 1)");
      }
      l=0; r=s;
      while (l<r)
      {
	m = (l+r) >> 1;
	if (table[m].x < ftest.x) l=m+1; else r=m;
      }
      if (r < s && table[r].x == ftest.x) break;
      s_addell(&ftest, &fg, cp4, p);
    }
    h += table[r].i * B;
    if (table[r].y == ftest.y) i = -i;
    h += s * i * B;

FOUND:
    h = sexact_order(h, &f, cp4, p);
    if (B == 1) B = h;
    else
    {
      GEN p1 = chinese(mkintmodu(smodss(A,B),B), mkintmodu(0,h));
      A = itos(gel(p1,2));
      if (is_bigint(p1[1])) { h = A; break; }
      B = itos(gel(p1,1));
    }

    i = (B < pordmin);
    if (i)
    {
      A = (p2p - A) % B;
      if ((A << 1) > B) A -= B;
    }
    /* h = A mod B, closest lift to p+1 */
    h = A + B * (((ulong)(p2p + B - (A << 1))) / (B << 1));
    avma = av; if (!i) break;
  }
  if (table) free(table);
  return stoi(KRO==1? p1p-h: h-p1p);
}

/** apell from CM (original code contributed by Mark Watkins) **/

static ulong
Mod16(GEN x) {
  long s = signe(x);
  ulong m;
  if (!s) return 0;
  m = mod16(x); if (!m) return m;
  if (s < 0) m = 16 - m;
  return m;
}
#define Mod2(x) (Mod16(x) & 1)
#define Mod4(x) (Mod16(x) & 3)
#define Mod8(x) (Mod16(x) & 7)

static GEN
ap_j0(GEN E,GEN p)
{
  GEN a, b, e, d;
  if (umodiu(p,3) != 1) return gen_0;
  (void)cornacchia2(utoipos(27),p, &a,&b);
  if (umodiu(a, 3) == 1) a = negi(a);
  d = Rg_to_Fp(gmulgs(gel(E,11), 8), p);
  e = diviuexact(shifti(p,-1), 3); /* (p-1) / 6 */
  return centermod(mulii(a, Fp_pow(d, e, p)), p);
}
static GEN
ap_j1728(GEN E,GEN p)
{
  GEN a, b, d, e;
  if (mod4(p) != 1) return gen_0;
  (void)cornacchia2(utoipos(4),p, &a,&b);
  if (Mod4(a)==0) a = b;
  if (Mod2(a)==1) a = shifti(a,1);
  if (Mod8(a)==6) a = negi(a);
  d = Rg_to_Fp(gmulgs(gel(E,10), -27), p);
  e = shifti(p,-2); /* (p-1) / 4 */
  return centermod(mulii(a, Fp_pow(d, e, p)), p);
}
static GEN
ap_j8000(GEN p)
{
  GEN a, b;
  long r = mod8(p);
  if (r != 1 && r != 3) return gen_0;
  (void)cornacchia2(utoipos(8),p, &a,&b);
  switch(Mod16(a)) {
    case 2: case 6:   if (Mod4(b)) a = negi(a);
      break;
    case 10: case 14: if (!Mod4(b)) a = negi(a);
      break;
  }
  return a;
}
static GEN
ap_j287496(GEN p)
{
  GEN a, b;
  if (mod4(p) != 1) return gen_0;
  (void)cornacchia2(utoipos(4),p, &a,&b);
  if (Mod4(a)==0) a = b;
  if (Mod2(a)==1) a = shifti(a,1);
  if (Mod8(a)==6) a = negi(a);
  if (krosi(2,p) < 0) a = negi(a);
  return a;
}
static GEN
ap_cm(int CM, GEN p)
{
  GEN a, b;
  if (krosi(CM,p) < 0) return gen_0;
  (void)cornacchia2(utoipos(-CM),p, &a, &b);
  if ((CM&3) == 0) CM >>= 2;
  if ((krois(a, -CM) > 0) ^ (CM == -7)) a = negi(a);
  return a;
}
static GEN
ec_ap_cm(GEN J,GEN C6B,GEN C6E,int CM,GEN jd,GEN jn,GEN p)
{
  GEN a;
  if (!equalii(modii(mulii(jd,J),p), jn)) return NULL;
  if      (CM == -8)  a = ap_j8000(p);
  else if (CM == -16) a = ap_j287496(p);
  else                a = ap_cm(CM,p);
  if (kronecker(mulii(C6E,C6B), p) < 0) a = negi(a);
  return a;
}

static GEN
ap_bad_red(GEN e, GEN p)
{
  pari_sp av = avma;
  GEN c6 = Rg_to_Fp(gel(e,11), p);
  long s = kronecker(c6, p);
  if (mod4(p) == 3) s = -s;
  avma = av; return stoi(s);
}
static GEN
u2tonegi(ulong a, ulong b) { GEN z = u2toi(a,b); setsigne(z, -1); return z; }

GEN
CM_ellap(GEN E, GEN p)
{
  pari_sp av = avma;
  GEN C4E, C6E, jn, jd, a, t, u;

  if (cmpiu(p, 99) < 0) return ap_jacobi(E, itou(p));
  if (!signe(Rg_to_Fp(gel(E,12), p))) { avma = av; return ap_bad_red(E,p); }
#define CHECK(CM,J,C6B) a = ec_ap_cm(J,C6B,C6E,CM,jd,jn,p); if (a) goto DONE;
  C4E = Rg_to_Fp(gel(E,10), p);
  if (!signe(C4E)) { a = ap_j0(E,p); goto DONE;}
  C6E = Rg_to_Fp(gel(E,11), p);
  if (!signe(C6E)) { a = ap_j1728(E,p); goto DONE;}
  jn = Rg_to_Fp(numer(gel(E,13)), p);
  jd = Rg_to_Fp(denom(gel(E,13)), p); /* j = jn/jd */
  CHECK(-7,  utoineg(3375),      utoipos(1323));
  CHECK(-8,  utoipos(8000),      utoineg(1792));
  CHECK(-12, utoipos(54000),     utoineg(19008));
  CHECK(-11, utoineg(32768),     utoineg(6776));
  CHECK(-16, utoipos(287496),    utoipos(12096));
  CHECK(-19, utoineg(884736),    utoineg(77976));
  CHECK(-27, utoineg(12288000),  utoineg(54648));
  CHECK(-7,  utoipos(16581375),  utoipos(75411));
  CHECK(-43, utoineg(884736000), utoineg(8387064));
  t = u2tonegi(0x00000022UL, 0x45ae8000UL); /* -27878400*5280 */
  CHECK(-67, t, utoineg(210408408));
  t = u2tonegi(0x03a4b862UL, 0xc4b40000UL); /* -640320^3 */
  u = u2tonegi(0x000000f8UL, 0x4414c858UL); /* -705220967*1512 */
  CHECK(-163, t, u);
#undef CHECK
  avma = av; return NULL;
DONE:
  return gerepileuptoint(av, icopy(a));
}

/* for ellsea() */
GEN
CM_CardEFp(GEN E, GEN p)
{
  GEN ap = CM_ellap(E, p);
  return ap? subii(addis(p,1), ap): gen_0;
}

GEN
apell(GEN e, GEN p)
{
  GEN a;
  checkell(e);
  if (typ(p)!=t_INT || signe(p) <= 0) pari_err(talker,"not a prime in apell");
  a = CM_ellap(e, p); if (a) return a;

  if (cmpiu(p, 0x3fffffff) > 0) return apell1(e, p);
  return apell0(e, itou(p));
}

GEN
ellap0(GEN e, GEN p, long flag) { return flag? apell2(e,p): apell(e,p); }

static void
checkell_int(GEN e)
{
  checkell(e);
  if (typ(e[1]) != t_INT || typ(e[2]) != t_INT || typ(e[3]) != t_INT
   || typ(e[4]) != t_INT || typ(e[5]) != t_INT)
    pari_err(talker,"not an integral model");
}

GEN
anell(GEN e, long n0)
{
  long tab[4]={0,1,1,-1}; /* p prime; (-1/p) = tab[p&3]. tab[0] not used */
  long P[3] = {evaltyp(t_INT)|_evallg(3), evalsigne(1)|evallgefint(3), 0};
  ulong p, m, SQRTn, n = (ulong)n0;
  GEN *an, D, c6;

  checkell_int(e);
  if (n0 <= 0) return cgetg(1,t_VEC);
  if (n >= LGBITS) pari_err(impl,"anell for n >= %lu", LGBITS);
  SQRTn = (ulong)sqrt(n);
  c6= gel(e,11);
  D = gel(e,12);

  an = (GEN*)cgetg(n+1,t_VEC); an[1] = gen_1;
  for (p=2; p <= n; p++) an[p] = NULL;
  for (p=2; p<=n; p++)
  {
    if (an[p]) continue; /* p not prime */

    if (!umodiu(D,p)) /* bad reduction, p | D */
      switch (tab[p&3] * krois(c6,p)) /* (-c6/p) */
      {
        case -1:  /* non deployee */
          for (m=p; m<=n; m+=p)
            if (an[m/p]) an[m] = negi(an[m/p]);
          continue;
        case 0:   /* additive */
          for (m=p; m<=n; m+=p) an[m] = gen_0;
          continue;
        case 1:   /* deployee */
          for (m=p; m<=n; m+=p)
            if (an[m/p]) an[m] = an[m/p];
          continue;
      }
    else /* good reduction */
    {
      GEN ap;
      P[2] = p; ap = apell(e, P);

      if (p <= SQRTn) {
        ulong pk, oldpk = 1;
        for (pk=p; pk <= n; oldpk=pk, pk *= p)
        {
          if (pk == p) an[pk] = ap;
          else
          {
            pari_sp av = avma;
            GEN u = mulii(ap, an[oldpk]);
            GEN v = mului(p, an[oldpk/p]);
            an[pk] = gerepileuptoint(av, subii(u,v));
          }
          for (m = n/pk; m > 1; m--)
            if (an[m] && m%p) an[m*pk] = mulii(an[m], an[pk]);
        }
      } else {
        an[p] = ap;
        for (m = n/p; m > 1; m--)
          if (an[m]) an[m*p] = mulii(an[m], ap);
      }
    }
  }
  return (GEN)an;
}

GEN
akell(GEN e, GEN n)
{
  long i, j, s, ex;
  pari_sp av = avma;
  GEN fa, P, E, D, c6, ap, u, v, w, y, p;

  checkell(e);
  if (typ(n) != t_INT) pari_err(typeer,"akell");
  if (signe(n)<= 0) return gen_0;
  if (gcmp1(n)) return gen_1;
  c6= gel(e,11);
  D = gel(e,12);
  if (typ(D) != t_INT) pari_err(talker,"not an integral model in akell");
  u = coprime_part(n, D);
  s = 1;
  if (!equalii(u, n))
  { /* bad reduction at primes dividing n/u */
    fa = Z_factor(diviiexact(n, u));
    P = gel(fa,1);
    E = gel(fa,2);
    for (i=1; i<lg(P); i++)
    {
      p = gel(P,i);
      j = kronecker(c6,p); if (!j) { avma = av; return gen_0; }
      if (mod2(gel(E,i)))
      {
        if (mod4(p) == 3) j = -j;
        if (j < 0) s = -s;
      }
    }
  }
  y = stoi(s); fa = Z_factor(u);
  P = gel(fa,1);
  E = gel(fa,2);
  for (i=1; i<lg(P); i++)
  { /* good reduction */
    p = gel(P,i);
    ex = itos(gel(E,i));
    ap = apell(e,p);
    u = ap; v = gen_1;
    for (j=2; j<=ex; j++)
    {
      w = subii(mulii(ap,u), mulii(p,v));
      v = u; u = w;
    }
    y = mulii(u,y);
  }
  return gerepileuptoint(av,y);
}

GEN
elllseries(GEN e, GEN s, GEN A, long prec)
{
  pari_sp av = avma, av1, lim;
  ulong l, n;
  long eps, flun;
  GEN z, cg, v, cga, cgb, s2, ns, gs, N, gr;

  if (!A) A = gen_1;
  else
  {
    if (gsigne(A)<=0)
      pari_err(talker,"cut-off point must be positive in lseriesell");
    if (gcmpgs(A,1) < 0) A = ginv(A);
  }
  if (isint(s, &s) && signe(s) <= 0) { avma = av; return gen_0; }
  flun = gcmp1(A) && gcmp1(s);
  checkell(e);
  e = ell_to_small(e); gr = ellglobalred(e);
  e = _coordch(e,gel(gr,2));
  N = gel(gr,1);
  eps = ellrootno_global(e, N);
  if (flun && eps < 0) { avma = av; return real_0(prec); }

  gs = ggamma(s, prec);
  cg = divrr(Pi2n(1, prec), gsqrt(N,prec));
  cga = gmul(cg, A);
  cgb = gdiv(cg, A);
  l = (ulong)((bit_accuracy_mul(prec, LOG2) +
              fabs(gtodouble(real_i(s))-1.) * log(rtodbl(cga)))
            / rtodbl(cgb) + 1);
  if ((long)l < 1) l = 1;
  v = anell(e, min(l,LGBITS-1));
  s2 = ns = NULL; /* gcc -Wall */
  if (!flun) { s2 = gsubsg(2,s); ns = gpow(cg, gsubgs(gmul2n(s,1),2),prec); }
  z = gen_0;
  av1 = avma; lim = stack_lim(av1,1);
  for (n = 1; n <= l; n++)
  {
    GEN p1, an, gn = utoipos(n);
    an = ((ulong)n<LGBITS)? gel(v,n): akell(e,gn);
    if (!signe(an)) continue;

    p1 = gdiv(incgam0(s,mulur(n,cga),gs,prec), gpow(gn,s,prec));
    if (flun)
      p1 = gmul2n(p1, 1);
    else
    {
      GEN p2 = gdiv(gmul(ns, incgam(s2,mulur(n,cgb),prec)), gpow(gn, s2,prec));
      if (eps < 0) p2 = gneg_i(p2);
      p1 = gadd(p1, p2);
    }
    z = gadd(z, gmul(p1, an));
    if (low_stack(lim, stack_lim(av1,1)))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"lseriesell");
      z = gerepilecopy(av1,z);
    }
  }
  return gerepileupto(av, gdiv(z,gs));
}

/********************************************************************/
/**                                                                **/
/**                       CANONICAL HEIGHT                         **/
/**                                                                **/
/********************************************************************/

/* h' := h_oo(a) + 1/2 log(denom(a)) */
static GEN
hell(GEN e, GEN a, long prec)
{
  long n;
  pari_sp av = avma;
  GEN p1, p2, y, z, q, pi2surw, qn, ps;

  checkbell(e);
  pi2surw = gdiv(Pi2n(1, prec), gel(e,15));
  z = gmul(real_i(zell(e,a,prec)), pi2surw);
  q = real_i( expIxy(pi2surw, gel(e,16), prec) );
  y = gsin(z,prec); qn = gen_1; ps = gneg_i(q);
  for (n = 3; ; n += 2)
  {
    qn = gmul(qn, ps);
    ps = gmul(ps, q);
    y = gadd(y, gmul(qn, gsin(gmulsg(n,z),prec)));
    if (gexpo(qn) < - bit_accuracy(prec)) break;
  }
  p1 = gmul(gsqr(gdiv(gmul2n(y,1), d_ellLHS(e,a))), pi2surw);
  p2 = gsqr(gsqr(gdiv(p1, gsqr(gsqr(denom(gel(a,1)))))));
  p1 = gdiv(gmul(p2,q), gel(e,12));
  p1 = gmul2n(glog(gabs(p1,prec),prec), -5);
  return gerepileupto(av, gneg(p1));
}

/* h' := h_oo(x) + 1/2 log(denom(x)) */
static GEN
hells(GEN e, GEN x, long prec)
{
  GEN b8 = gel(e,9), b6 = gel(e,8), b4 = gel(e,7), b2 = gel(e,6);
  GEN w, z, t, mu, b42, b62;
  long n, lim;

  t = gdiv(real_1(prec), gel(x,1));
  mu = gmul2n(glog(numer(gel(x,1)),prec),-1);
  b42 = gmul2n(b4,1);
  b62 = gmul2n(b6,1);
  lim = 15 + bit_accuracy(prec);
  for (n = 3; n < lim; n += 2)
  {
    /* 4 + b2 t + 2b4 t^2 + b6 t^3 */
    w = gmul(t, gaddsg(4, gmul(t, gadd(b2, gmul(t, gadd(b42, gmul(t, b6)))))));
    /* 1 - (b4 t^2 + 2b6 t^3 + b8 t^4) */
    z = gsub(gen_1, gmul(gsqr(t), gadd(b4, gmul(t, gadd(b62, gmul(t, b8))))));
    mu = gadd(mu, gmul2n(glog(z,prec), -n));
    t = gdiv(w, z);
  }
  return mu;
}

static GEN
hell2(GEN e, GEN x, long prec)
{
  GEN e3, ro, v, D;
  pari_sp av = avma;

  if (is_inf(x)) return gen_0;
  D = gel(e,12);
  ro= gel(e,14);
  e3 = (gsigne(D) < 0)? gel(ro,1): gel(ro,3);
  v = init_ch(); gel(v,2) = addis(gfloor(e3),-1);
  return gerepileupto(av, hells(_coordch(e,v), pointch(x,v), prec));
}

/* exp( h_oo(z) ), assume z on neutral component.
 * If flag, return exp(4 h_oo(z)) instead */
static GEN
exphellagm(GEN e, GEN z, int flag, long prec)
{
  GEN x_a, a, b, r, V = cgetg(1, t_VEC), x = gel(z,1);
  long n, ex = 5-bit_accuracy(prec);

  x = new_coords(e, x, &a,&b, 0, prec);
  x_a = gsub(x, a);
  if (gsigne(a) > 0)
  {
    GEN a0 = a;
    x = gsub(x, b);
    a = gneg(b);
    b = gsub(a0, b);
  }
  a = gsqrt(gneg(a), prec);
  b = gsqrt(gneg(b), prec);
  /* compute height on isogenous curve E1 ~ E0 */
  for(n=0; ; n++)
  {
    GEN p1, p2, ab, a0 = a;
    a = gmul2n(gadd(a0,b), -1);
    r = gsub(a, a0);
    if (gcmp0(r) || gexpo(r) < ex) break;
    ab = gmul(a0, b);
    b = gsqrt(ab, prec);

    p1 = gmul2n(gsub(x, ab), -1);
    p2 = gsqr(a);
    x = gadd(p1, gsqrt(gadd(gsqr(p1), gmul(x, p2)), prec));
    V = shallowconcat(V, gadd(x, p2));
  }
  if (n) {
    x = gel(V,n);
    while (--n > 0) x = gdiv(gsqr(x), gel(V,n));
  } else {
    x = gadd(x, gsqr(a));
  }
  /* height on E1 is log(x)/2. Go back to E0 */
  return flag? gsqr( gdiv(gsqr(x), x_a) )
             : gdiv(x, sqrtr( mpabs(x_a) ));
}
/* exp( 4h_oo(z) ) */
static GEN
exp4hellagm(GEN e, GEN z, long prec)
{
  GEN e1 = gmael(e,14,1), x = gel(z,1);
  if (gcmp(x, e1) < 0) /* z not on neutral component */
  {
    GEN eh = exphellagm(e, addell(e, z,z), 0, prec);
    /* h_oo(2P) = 4h_oo(P) - log |2y + a1x + a3| */
    return gmul(eh, gabs(d_ellLHS(e, z), prec));
  }
  return exphellagm(e, z, 1, prec);
}

GEN
ellheightoo(GEN e, GEN z, long prec)
{
  GEN e1, h, x = gel(z,1);
  pari_sp av = avma;
  checkell(e);
  e1 = gmael(e,14,1);
  if (gcmp(x, e1) < 0) /* z not on neutral component */
  {
    GEN eh = exphellagm(e, addell(e, z,z), 0, prec);
    /* h_oo(2P) = 4h_oo(P) - log |2y + a1x + a3| */
    h = gmul(eh, gabs(d_ellLHS(e, z), prec));
  }
  else
    h = exphellagm(e, z, 1, prec);
  return gerepileuptoleaf(av, gmul2n(mplog(h), -2));
}

/* Assume e integral, given by a minimal model */
GEN
ellheight0(GEN e, GEN a, long flag, long prec)
{
  long i, tx = typ(a), lx = lg(a);
  pari_sp av = avma;
  GEN Lp, x, y, z, phi2, psi2, psi3;

  if (flag > 2 || flag < 0) pari_err(flagerr,"ellheight");
  checkbell(e); if (!is_matvec_t(tx)) pari_err(elliper1);
  lx = lg(a); if (lx==1) return cgetg(1,tx);
  tx = typ(a[1]);
  if (is_matvec_t(tx))
  {
    z = cgetg(lx,tx);
    for (i=1; i<lx; i++) gel(z,i) = ellheight0(e,gel(a,i),flag,prec);
    return z;
  }
  if (is_inf(a)) return gen_0;
  if (!oncurve(e,a)) pari_err(talker, "point not on elliptic curve");

  psi2 = numer(d_ellLHS(e,a));
  if (!signe(psi2)) { avma = av; return gen_0; }
  switch(flag)
  {
    case 0:  z = hell2(e,a,prec); break; /* Tate 4^n */
    case 1:  z = hell(e,a,prec);  break; /* Silverman's log(sigma) */
    default:
    {
      GEN d = denom(gel(a,1));
      z = exp4hellagm(e,a,prec); /* = exp(4h_oo(a)), Mestre's AGM */
      if (!is_pm1(d)) z = gmul(z, sqri(d));
      z = gmul2n(mplog(z), -2); break;
    }
  }
  x = gel(a,1);
  y = gel(a,2);
  psi3 = numer( /* b8 + 3x b6 + 3x^2 b4 + x^3 b2 + 3 x^4 */
     gadd(gel(e,9), gmul(x,
     gadd(gmulsg(3,gel(e,8)), gmul(x,
     gadd(gmulsg(3,gel(e,7)), gmul(x, gadd(gel(e,6), gmulsg(3,x)))))))) );
  if (!signe(psi3)) { avma=av; return gen_0; }

  phi2 = numer( /* a4 + 2a2 x + 3x^2 - y a1*/
    gsub(gadd(gel(e,4),gmul(x,gadd(shifti(gel(e,2),1),gmulsg(3,x)))),
         gmul(gel(e,1),y)) );
  Lp = (GEN)Z_factor(gcdii(psi2,phi2))[1];
  lx = lg(Lp);
  for (i=1; i<lx; i++)
  {
    GEN p = gel(Lp,i);
    long u, v, n, n2;
    if (signe(remii(gel(e,10),p)))
    { /* p \nmid c4 */
      long N = Z_pval(gel(e,12),p);
      if (!N) continue;
      n2 = Z_pval(psi2,p); n = n2<<1;
      if (n > N) n = N;
      u = n * ((N<<1) - n);
      v = N << 3;
    }
    else
    {
      n2 = Z_pval(psi2, p);
      n  = Z_pval(psi3, p);
      if (n >= 3*n2) { u = n2; v = 3; } else { u = n; v = 8; }
    }
    /* z -= u log(p) / v */
    z = gadd(z, divrs(mulsr(-u, glog(p,prec)), v));
  }
  return gerepileupto(av, gmul2n(z, 1));
}

GEN
ghell2(GEN e, GEN a, long prec) { return ellheight0(e,a,0,prec); }

GEN
ghell(GEN e, GEN a, long prec) { return ellheight0(e,a,2,prec); }

GEN
mathell(GEN e, GEN x, long prec)
{
  GEN y, h, pdiag;
  long lx = lg(x),i,j,tx=typ(x);
  pari_sp av = avma;

  if (!is_vec_t(tx)) pari_err(elliper1);
  y = cgetg(lx,t_MAT); pdiag = new_chunk(lx);
  for (i=1; i<lx; i++)
  {
    gel(pdiag,i) = ghell(e,gel(x,i),prec);
    gel(y,i) = cgetg(lx,t_COL);
  }
  for (i=1; i<lx; i++)
  {
    coeff(y,i,i) = pdiag[i];
    for (j=i+1; j<lx; j++)
    {
      h = ghell(e, addell(e,gel(x,i),gel(x,j)), prec);
      h = gsub(h, gadd(gel(pdiag,i),gel(pdiag,j)));
      gcoeff(y,j,i) = gcoeff(y,i,j) = gmul2n(h, -1);
    }
  }
  return gerepilecopy(av,y);
}

static GEN
bilhells(GEN e, GEN z1, GEN z2, GEN h2, long prec)
{
  long lz1=lg(z1), tx, i;
  pari_sp av = avma;
  GEN y,p1,p2;

  if (lz1==1) return cgetg(1,typ(z1));

  tx = typ(z1[1]);
  if (!is_matvec_t(tx))
  {
    p1 = ghell(e, addell(e,z1,z2),prec);
    p2 = gadd(h2, ghell(e,z1,prec));
    return gerepileupto(av, gmul2n(gsub(p1,p2), -1));
  }
  y = cgetg(lz1, typ(z1));
  for (i=1; i<lz1; i++) gel(y,i) = bilhells(e,gel(z1,i),z2,h2,prec);
  return y;
}

GEN
bilhell(GEN e, GEN z1, GEN z2, long prec)
{
  GEN p1, h2;
  long tz1 = typ(z1), tz2 = typ(z2);
  pari_sp av = avma;

  if (!is_matvec_t(tz1) || !is_matvec_t(tz2)) pari_err(elliper1);
  if (lg(z1)==1) return cgetg(1,tz1);
  if (lg(z2)==1) return cgetg(1,tz2);

  tz1 = typ(z1[1]);
  tz2 = typ(z2[1]);
  if (is_matvec_t(tz2))
  {
    if (is_matvec_t(tz1)) pari_err(talker,"two vector/matrix types in bilhell");
    p1 = z1; z1 = z2; z2 = p1;
  }
  h2 = ghell(e,z2,prec);
  return gerepileupto(av, bilhells(e,z1,z2,h2,prec));
}

/********************************************************************/
/**                                                                **/
/**                    Modular Parametrization                     **/
/**                                                                **/
/********************************************************************/

GEN
elltaniyama(GEN e, long prec)
{
  GEN x, w, c, d, s1, s2, s3, X, C;
  long n, m;
  pari_sp av=avma, tetpil;

  checkell(e); x = cgetg(prec+3,t_SER);
  x[1] = evalsigne(1) | evalvalp(-2) | evalvarn(0);
  gel(x,2) = gen_1;
  d = ginv(gtoser(anell(e,prec+1), 0)); setvalp(d,-1);
  /* 2y(t) + a1x + a3 = d tx'(t). Solve for x(t),y(t):
   * 4y^2 = 4x^3 + b2 x^2 + 2b4 x + b6 */

  if (!prec) goto END;
  c = gsqr(d);
  /* 4x^3 + b2 x^2 + 2b4 x + b6 = c (t x'(t))^2; c = 1/t^2 + O(1/t) */
  C = c+4; /* C[i] = coeff(c, t^i) */
  X = x+4;
  /* n = -3 */
  gel(X,-1) = gmul2n(gmul(gel(X,-2),gel(C,-1)), -1);
  for (n=-2; n <= prec-4; n++)
  {
    if (n != 2)
    {
      s3 = gmul(gel(e,6),gel(X,n));
      if (!n) s3 = gadd(s3, gel(e,7));
      s2 = gen_0;
      for (m=-2; m<=n+1; m++)
	s2 = gadd(s2,gmulsg(m*(n+m),gmul(gel(X,m),gel(C,n-m))));
      s2 = gmul2n(s2,-1);
      s1 = gen_0;
      for (m=-1; m+m<=n; m++)
      {
	if (m+m==n)
          s1 = gadd(s1, gsqr(gel(X,m)));
	else
          s1 = gadd(s1, gmul2n(gmul(gel(X,m),gel(X,n-m)),1));
      }
      gel(X,n+2) = gdivgs(gsub(gadd(gmulsg(6,s1),s3),s2), (n+2)*(n+1)-12);
    }
    else
    {
      setlg(x, 9); gel(x,8) = pol_x[MAXVARN];
      w = derivser(x); setvalp(w,-2); /* 4v^3 + b2 x^2 + 2b4 x + b6 */
      s1 = gadd(gel(e,8), gmul(x, gadd(gmul2n(gel(e,7),1),
                                        gmul(x,gadd(gel(e,6),gmul2n(x,2))))));
      setlg(x, prec+3);
      s2 = gsub(s1, gmul(c,gsqr(w)));
      s2 = gel(s2,2);
      gel(X,n+2) = gneg(gdiv(gel(s2,2),gel(s2,3)));
    }
  }
END:
  w = gmul(d,derivser(x)); setvalp(w, valp(w)+1);
  w = gsub(w, ellLHS0(e,x));
  tetpil = avma; s1 = cgetg(3,t_VEC);
  gel(s1,1) = gcopy(x);
  gel(s1,2) = gmul2n(w,-1); return gerepile(av,tetpil,s1);
}

/********************************************************************/
/**                                                                **/
/**                       TORSION POINTS (over Q)                  **/
/**                                                                **/
/********************************************************************/
static int
smaller_x(GEN p, GEN q)
{
  int s = absi_cmp(denom(p), denom(q));
  return (s<0 || (s==0 && absi_cmp(numer(p),numer(q)) < 0));
}

/* best generator in cycle of length k */
static GEN
best_in_cycle(GEN e, GEN p, long k)
{
  GEN p0 = p,q = p;
  long i;

  for (i=2; i+i<k; i++)
  {
    q = addell(e,q,p0);
    if (cgcd(i,k)==1 && smaller_x(gel(q,1), gel(p,1))) p = q;
  }
  return (gsigne(d_ellLHS(e,p)) < 0)? invell(e,p): p;
}

/* <p,q> = E_tors, possibly NULL (= oo), p,q independent unless NULL
 * order p = k, order q = 2 unless NULL */
static GEN
tors(GEN e, long k, GEN p, GEN q, GEN v)
{
  GEN r;
  if (q)
  {
    long n = k>>1;
    GEN p1, best = q, np = powell(e,p,utoipos(n));
    if (n % 2 && smaller_x(gel(np,1), gel(best,1))) best = np;
    p1 = addell(e,q,np);
    if (smaller_x(gel(p1,1), gel(best,1))) q = p1;
    else if (best == np) { p = addell(e,p,q); q = np; }
    p = best_in_cycle(e,p,k);
    if (v)
    {
      p = pointch(p,v);
      q = pointch(q,v);
    }
    r = cgetg(4,t_VEC);
    gel(r,1) = utoipos(2*k);
    gel(r,2) = mkvec2(utoipos(k), gen_2);
    gel(r,3) = mkvec2copy(p, q);
  }
  else
  {
    if (p)
    {
      p = best_in_cycle(e,p,k);
      if (v) p = pointch(p,v);
      r = cgetg(4,t_VEC);
      gel(r,1) = utoipos(k);
      gel(r,2) = mkvec( gel(r,1) );
      gel(r,3) = mkvec( gcopy(p) );
    }
    else
    {
      r = cgetg(4,t_VEC);
      gel(r,1) = gen_1;
      gel(r,2) = cgetg(1,t_VEC);
      gel(r,3) = cgetg(1,t_VEC);
    }
  }
  return r;
}

/* assume e is defined over Q (use Mazur's theorem) */
static long
_orderell(GEN e, GEN p)
{
  pari_sp av = avma;
  GEN p1 = p;
  long k;
  for (k = 1; k < 16; k++)
  {
    if (is_inf(p1)) { avma = av; return k; }
    p1 = addell(e, p1, p);
  }
  avma = av; return 0;
}
GEN
orderell(GEN e, GEN p)
{
  long t;
  checkell(e); checkpt(p); t = typ(e[13]);
  if (!is_rational_t(t)) pari_err(impl,"orderell for nonrational elliptic curves");
  return utoi( _orderell(e, p) );
}

/* Using Lutz-Nagell */

/* p in Z[X] of degree 3. Return vector of x/4, x integral root of p */
GEN
ratroot(GEN p)
{
  GEN L, a, ld;
  long i, t, v = ZX_valuation(p, &p);

  if (v == 3) return mkvec(gen_0);
  if (v == 2) return mkvec2(gen_0, gmul2n(negi(gel(p,2)), -2));

  L = cgetg(4,t_VEC); t = 1;
  if (v == 1) gel(L,t++) = gen_0;
  ld = divisors(gel(p,2));
  for (i=1; i<lg(ld); i++)
  {
    a = gel(ld,i);
    if (!signe(poleval(p,a))) gel(L,t++) = gmul2n(a, -2);
    a = negi(a);
    if (!signe(poleval(p,a))) gel(L,t++) = gmul2n(a, -2);
  }
  setlg(L,t); return L;
}

static int
is_new_torsion(GEN e, GEN v, GEN p, long t2) {
  GEN pk = p, pkprec = NULL;
  long k,l;

  for (k=2; k<=6; k++)
  {
    pk = addell(e,pk,p); /* = [k] p */
    if (is_inf(pk)) return 1;

    for (l=2; l<=t2; l++)
      if (gequal(gel(pk,1),gmael(v,l,1))) return 1;

    if (pkprec && k<=5)
      if (gequal(gel(pk,1),gel(pkprec,1))) return 1;
    pkprec=pk;
  }
  return 0;
}

static GEN
nagelllutz(GEN e)
{
  GEN ld, pol, p1, lr, r, v, w2, w3;
  long i, j, nlr, t, t2, k, k2;
  pari_sp av=avma;

  v = ellintegralmodel(e);
  if (v) e = _coordch(e,v);
  pol = RgX_rescale(RHSpol(e), utoipos(4));
  r = cgetg(17, t_VEC);
  gel(r,1) = mkvec(gen_0);
  lr = ratroot(pol); nlr=lg(lr)-1;
  for (t=1,i=1; i<=nlr; i++)
  {
    GEN x = gel(lr,i), y = gmul2n(gneg(ellLHS0(e,x)), -1);
    gel(r,++t) = mkvec2(x, y);
  }
  ld = Z_factor(gmul2n(absi(gel(e,12)), 4));
  p1 = gel(ld,2); k = lg(p1);
  for (i=1; i<k; i++) gel(p1,i) = shifti(gel(p1,i), -1);
  ld = divisors(ld);
  for (t2=t,j=1; j<lg(ld); j++)
  {
    GEN d = gel(ld,j);
    lr = ratroot(gsub(pol, shifti(sqri(d), 6)));
    for (i=1; i<lg(lr); i++)
    {
      GEN x = gel(lr,i), y = gmul2n(gadd(d, gneg(ellLHS0(e,x))), -1);
      p1 = mkvec2(x, y);
      if (is_new_torsion(e,r,p1,t2))
      {
	gel(r,++t) = p1;
        gel(r,++t) = mkvec2(x, gsub(y, d));
      }
    }
  }
  if (t == 1) { avma = av; return tors(e,1,NULL,NULL,v); }

  if (nlr < 3)
  {
    w2 = mkvec( utoipos(t) );
    for (k=2; k<=t; k++)
      if (_orderell(e,gel(r,k)) == t) break;
    if (k>t) pari_err(bugparier,"torsell (bug1)");

    w3 = mkvec( gel(r,k) );
  }
  else
  {
    if (t&3) pari_err(bugparier,"torsell (bug2)");
    t2 = t>>1;
    w2 = mkvec2(utoipos(t2), gen_2);
    for (k=2; k<=t; k++)
      if (_orderell(e,gel(r,k)) == t2) break;
    if (k>t) pari_err(bugparier,"torsell (bug3)");

    p1 = powell(e,gel(r,k),utoipos(t>>2));
    k2 = (!is_inf(p1) && gequal(gel(r,2),p1))? 3: 2;
    w3 = mkvec2(gel(r,k), gel(r,k2));
  }
  if (v)
  {
    gel(v,1) = ginv(gel(v,1));
    w3 = pointch(w3,v);
  }
  return gerepilecopy(av, mkvec3(utoipos(t), w2,w3));
}

/* Using Doud's algorithm */

/* finds a bound for #E_tor */
static long
torsbound(GEN e)
{
  long m, b, bold, prime = 2;
  pari_sp av = avma;
  byteptr p = diffptr;
  GEN D = gel(e,12);
  long n = bit_accuracy(lgefint(D)) >> 3;
  /* n = number of primes to try ~ 1 prime every 8 bits in D */
  b = bold = 5040; /* = 2^4 * 3^2 * 5 * 7 */
  m = 0; p++;
  while (m < n)
  {
    NEXT_PRIME_VIADIFF_CHECK(prime,p);
    if (umodiu(D, prime))
    {
      b = cgcd(b, prime+1 - itos(apell0(e, (ulong)prime)));
      avma = av;
      if (b == 1) break;
      if (b == bold) m++; else { bold = b; m = 0; }
    }
  }
  return b;
}

static GEN
myround(GEN x, long *e)
{
  GEN y = grndtoi(x,e);
  if (*e > -5 && bit_accuracy(gprecision(x)) < gexpo(y) - 10)
    pari_err(talker, "ellinit data not accurate enough. Increase precision");
  return y;
}

/* E the curve, w in C/Lambda ~ E of order n, returns q = pointell(w) as a
 * rational point on the curve, or NULL if q is not rational. */
static GEN
torspnt(GEN E, GEN w, long n, long prec)
{
  GEN p = cgetg(3,t_VEC), q = pointell(E, w, prec);
  long e;
  gel(p,1) = gmul2n(myround(gmul2n(gel(q,1),2), &e),-2);
  if (e > -5 || typ(p[1]) == t_COMPLEX) return NULL;
  gel(p,2) = gmul2n(myround(gmul2n(gel(q,2),3), &e),-3);
  if (e > -5 || typ(p[2]) == t_COMPLEX) return NULL;
  return (oncurve(E,p)
      && is_inf(powell(E,p,utoipos(n)))
      && _orderell(E,p) == n)? p: NULL;
}

GEN
torsell(GEN e)
{
  long B, i, ord, pr, prec, k = 1;
  pari_sp av=avma;
  GEN v,w,w1,w22,w1j,w12,p,tor1,tor2;

  checkbell(e);
  v = ellintegralmodel(e);
  if (v) e = _coordch(e,v);

  B = torsbound(e); /* #E_tor | B */
  if (B == 1) { avma = av; return tors(e,1,NULL,NULL, v); }

  pr = DEFAULTPREC + ((lgefint(gel(e,12))-2) >> 1); /* pr >= size of sqrt(D) */
  w1 = gel(e,15);
  prec = precision(w1);
  if (prec < pr) pari_err(precer, "torsell");
  if (pr < prec) { prec = pr; e = gprec_w(e, pr); w1 = gel(e,15); }
  if (v) gel(v,1) = ginv(gel(v,1));
  w22 = gmul2n(gel(e,16),-1);
  if (B % 4)
  { /* cyclic of order 1, p, 2p, p <= 5 */
    p = NULL;
    for (i=10; i>1; i--)
    {
      if (B%i != 0) continue;
      w1j = gdivgs(w1,i);
      p = torspnt(e,w1j,i,prec);
      if (!p && i%2==0)
      {
        p = torspnt(e,gadd(w22,w1j),i,prec);
        if (!p) p = torspnt(e,gadd(w22,gmul2n(w1j,1)),i,prec);
      }
      if (p) { k = i; break; }
    }
    return gerepileupto(av, tors(e,k,p,NULL, v));
  }

  ord = 0; tor1 = tor2 = NULL;
  w12 = gmul2n(w1,-1);
  if ((p = torspnt(e,w12,2,prec)))
  {
    tor1 = p; ord++;
  }
  w = w22;
  if ((p = torspnt(e,w,2,prec)))
  {
    tor2 = p; ord += 2;
  }
  if (!ord)
  {
    w = gadd(w12,w22);
    if ((p = torspnt(e,w,2,prec)))
    {
      tor2 = p; ord += 2;
    }
  }
  p = NULL;
  switch(ord)
  {
    case 0: /* no point of order 2 */
      for (i=9; i>1; i-=2)
      {
        if (B%i != 0) continue;
        w1j = gdivgs(w1,i);
        p = torspnt(e,w1j,i,prec);
        if (p) { k = i; break; }
      }
      break;

    case 1: /* 1 point of order 2: w1 / 2 */
      for (i=12; i>2; i-=2)
      {
        if (B%i != 0) continue;
        w1j = gdivgs(w1,i);
        p = torspnt(e,w1j,i,prec);
        if (!p && i%4==0)
          p = torspnt(e,gadd(w22,w1j),i,prec);
        if (p) { k = i; break; }
      }
      if (!p) { p = tor1; k = 2; }
      break;

    case 2: /* 1 point of order 2: w = w2/2 or (w1+w2)/2 */
      for (i=5; i>1; i-=2)
      {
        if (B%i != 0) continue;
        w1j = gdivgs(w1,i);
        p = torspnt(e,gadd(w,w1j),2*i,prec);
        if (p) { k = 2*i; break; }
      }
      if (!p) { p = tor2; k = 2; }
      tor2 = NULL; break;

    case 3: /* 2 points of order 2: w1/2 and w2/2 */
      for (i=8; i>2; i-=2)
      {
        if (B%(2*i) != 0) continue;
        w1j = gdivgs(w1,i);
        p = torspnt(e,w1j,i,prec);
        if (p) { k = i; break; }
      }
      if (!p) { p = tor1; k = 2; }
      break;
  }
  return gerepileupto(av, tors(e,k,p,tor2, v));
}

GEN
elltors0(GEN e, long flag)
{
  switch(flag)
  {
    case 0: return torsell(e);
    case 1: return nagelllutz(e);
    default: pari_err(flagerr,"torsell");
  }
  return NULL; /* not reached */
}
