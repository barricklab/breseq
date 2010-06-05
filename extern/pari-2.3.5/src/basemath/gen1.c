/* $Id: gen1.c 11470 2008-12-19 19:33:27Z kb $

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
/**                      GENERIC OPERATIONS                        **/
/**                         (first part)                           **/
/**                                                                **/
/********************************************************************/
#include "pari.h"
#include "paripriv.h"

#define fix_frac(z) if (signe(z[2])<0)\
{\
  setsigne(z[1],-signe(z[1]));\
  setsigne(z[2],1);\
}

/* assume z[1] was created last */
#define fix_frac_if_int(z) if (is_pm1(z[2]))\
  z = gerepileupto((pari_sp)(z+3), gel(z,1));

/* assume z[1] was created last */
#define fix_frac_if_int_GC(z,tetpil) { if (is_pm1(z[2]))\
  z = gerepileupto((pari_sp)(z+3), gel(z,1));\
else\
  gerepilecoeffssp((pari_sp)z, tetpil, z+1, 2); }

static long
kro_quad(GEN x, GEN y)
{
  long k;
  pari_sp av=avma;

  x = subii(sqri(gel(x,3)), shifti(gel(x,2),2));
  k = kronecker(x,y); avma=av; return k;
}

/* is -1 not a square in Zp, assume p prime */
INLINE int
Zp_nosquare_m1(GEN p) { return (mod4(p) & 2); /* 2 or 3 mod 4 */ }

static GEN addpp(GEN x, GEN y);
static GEN mulpp(GEN x, GEN y);
static GEN divpp(GEN x, GEN y);
/* Argument codes for inline routines
 * c: complex, p: padic, q: quad, f: floating point (REAL, some complex)
 * R: without imaginary part (INT, REAL, INTMOD, FRAC, PADIC if -1 not square)
 * T: some type (to be converted to PADIC)
 */
static GEN
addRc(GEN x, GEN y) {
  GEN z = cgetg(3,t_COMPLEX);
  gel(z,1) = gadd(x,gel(y,1));
  gel(z,2) = gcopy(gel(y,2)); return z;
}
static GEN
mulRc(GEN x, GEN y) {
  GEN z = cgetg(3,t_COMPLEX);
  gel(z,1) = gmul(x,gel(y,1));
  gel(z,2) = gmul(x,gel(y,2)); return z;
}
static GEN
divRc(GEN x, GEN y) {
  GEN a, b, N, z = cgetg(3,t_COMPLEX);
  pari_sp tetpil, av = avma;
  a = gmul(x, gel(y,1));
  b = gmul(x, gel(y,2)); if(!gcmp0(b)) b = gneg_i(b);
  N = cxnorm(y); tetpil = avma;
  gel(z,1) = gdiv(a, N);
  gel(z,2) = gdiv(b, N); gerepilecoeffssp(av,tetpil,z+1,2); return z;
}
static GEN
divcR(GEN x, GEN y) {
  GEN z = cgetg(3,t_COMPLEX);
  gel(z,1) = gdiv(gel(x,1), y);
  gel(z,2) = gdiv(gel(x,2), y); return z;
}
static GEN
addRq(GEN x, GEN y) {
  GEN z = cgetg(4,t_QUAD);
  gel(z,1) = gcopy(gel(y,1));
  gel(z,2) = gadd(x, gel(y,2));
  gel(z,3) = gcopy(gel(y,3)); return z;
}
static GEN
mulRq(GEN x, GEN y) {
  GEN z = cgetg(4,t_QUAD);
  gel(z,1) = gcopy(gel(y,1));
  gel(z,2) = gmul(x,gel(y,2));
  gel(z,3) = gmul(x,gel(y,3)); return z;
}
static GEN
addqf(GEN x, GEN y, long prec) { pari_sp av = avma;
  long i = gexpo(x) - gexpo(y);
  if (i <= 0) i = 0; else i >>= TWOPOTBITS_IN_LONG;
  return gerepileupto(av, gadd(y, quadtoc(x, prec + i)));
}
static GEN
mulqf(GEN x, GEN y, long prec) { pari_sp av = avma;
  return gerepileupto(av, gmul(y, quadtoc(x, prec)));
}
static GEN
divqf(GEN x, GEN y, long prec) { pari_sp av = avma;
  return gerepileupto(av, gdiv(quadtoc(x,prec), y));
}
static GEN
divfq(GEN x, GEN y, long prec) { pari_sp av = avma;
  return gerepileupto(av, gdiv(x, quadtoc(y,prec)));
}
/* y PADIC, x + y by converting x to padic */
static GEN
addTp(GEN x, GEN y) { pari_sp av = avma; GEN z;

  if (!valp(y)) z = cvtop2(x,y);
  else {
    long l = signe(y[4])? valp(y) + precp(y): valp(y);
    z  = cvtop(x, gel(y,2), l);
  }
  return gerepileupto(av, addpp(z, y));
}
/* y PADIC, x * y by converting x to padic */
static GEN
mulTp(GEN x, GEN y) { pari_sp av = avma;
  return gerepileupto(av, mulpp(cvtop2(x,y), y));
}
/* y PADIC, non zero x / y by converting x to padic */
static GEN
divTp(GEN x, GEN y) { pari_sp av = avma;
  return gerepileupto(av, divpp(cvtop2(x,y), y));
}
/* x PADIC, x / y by converting y to padic */
static GEN
divpT(GEN x, GEN y) { pari_sp av = avma;
  return gerepileupto(av, divpp(x, cvtop2(y,x)));
}

/* z := Mod(x,X) + Mod(y,X) [ t_INTMOD preallocated ], x,y,X INT, 0 <= x,y < X
 * clean memory from z on */
static GEN
add_intmod_same(GEN z, GEN X, GEN x, GEN y) {
  if (lgefint(X) == 3) {
    ulong u = Fl_add(itou(x),itou(y), X[2]);
    avma = (pari_sp)z; gel(z,2) = utoi(u);
  }
  else {
    GEN u = addii(x,y); if (cmpii(u, X) >= 0) u = subii(u, X);
    gel(z,2) = gerepileuptoint((pari_sp)z, u);
  }
  gel(z,1) = icopy(X); return z;
}
/* cf add_intmod_same */
static GEN
mul_intmod_same(GEN z, GEN X, GEN x, GEN y) {
  if (lgefint(X) == 3) {
    ulong u = Fl_mul(itou(x),itou(y), X[2]);
    avma = (pari_sp)z; gel(z,2) = utoi(u);
  }
  else
    gel(z,2) = gerepileuptoint((pari_sp)z, remii(mulii(x,y), X) );
  gel(z,1) = icopy(X); return z;
}
/* cf add_intmod_same */
static GEN
div_intmod_same(GEN z, GEN X, GEN x, GEN y)
{
  if (lgefint(X) == 3) {
    ulong m = (ulong)X[2], u = Fl_div(itou(x), itou(y), m);
    avma = (pari_sp)z; gel(z,2) = utoi(u);
  }
  else
    gel(z,2) = gerepileuptoint((pari_sp)z, remii(mulii(x, Fp_inv(y,X)), X) );
  gel(z,1) = icopy(X); return z;
}

/*******************************************************************/
/*                                                                 */
/*        REDUCTION to IRREDUCIBLE TERMS (t_FRAC/t_RFRAC)          */
/*                                                                 */
/* (static routines are not memory clean, but OK for gerepileupto) */
/*******************************************************************/
/* d a t_POL, n a coprime t_POL of same var or "scalar". Not memory clean */
GEN
gred_rfrac_simple(GEN n, GEN d)
{
  GEN c, cn, cd, z;

  cd = content(d);
  cn = (typ(n) == t_POL && varn(n) == varn(d))? content(n): n;
  if (!gcmp1(cd)) {
    d = RgX_Rg_div(d,cd);
    if (!gcmp1(cn))
    {
      if (gcmp0(cn)) {
        n = (cn != n)? RgX_Rg_div(n,cd): gdiv(n, cd);
        c = gen_1;
      } else {
        n = (cn != n)? RgX_Rg_div(n,cn): gen_1;
        c = gdiv(cn,cd);
      }
    }
    else
      c = ginv(cd);
  } else {
    if (!gcmp1(cn))
    {
      if (gcmp0(cn)) {
        c = gen_1;
      } else {
        n = (cn != n)? RgX_Rg_div(n,cn): gen_1;
        c = cn;
      }
    } else {
      GEN y = cgetg(3,t_RFRAC);
      gel(y,1) = gcopy(n);
      gel(y,2) = gcopy(d); return y;
    }
  }

  if (typ(c) == t_POL)
  {
    z = c;
    do { z = content(z); } while (typ(z) == t_POL);
    cd = denom(z);
    cn = gmul(c, cd);
  }
  else
  {
    cn = numer(c);
    cd = denom(c);
  }
  z = cgetg(3,t_RFRAC);
  gel(z,1) = gmul(n, cn);
  gel(z,2) = gmul(d, cd); return z;
}

static GEN
fix_rfrac(GEN x, long d)
{
  GEN z, N, D;
  if (!d) return x;
  z = cgetg(3, t_RFRAC);
  N = gel(x,1);
  D = gel(x,2);
  if (d > 0) {
    gel(z, 1) = (typ(N)==t_POL && varn(N)==varn(D))? RgX_shift(N,d)
                                                   : monomialcopy(N,d,varn(D));
    gel(z, 2) = gcopy(D);
  } else {
    gel(z, 1) = gcopy(N);
    gel(z, 2) = RgX_shift(D, -d);
  }
  return z;
}

/* assume d != 0 */
static GEN
gred_rfrac2_i(GEN n, GEN d)
{
  GEN y, z;
  long tn, td, v;

  n = simplify_i(n);
  if (isexactzero(n)) return gcopy(n);
  d = simplify_i(d);
  td = typ(d);
  if (td!=t_POL || varncmp(varn(d), gvar(n)) > 0) return gdiv(n,d);
  tn = typ(n);
  if (tn!=t_POL)
  {
    if (varncmp(varn(d), gvar2(n)) < 0) return gred_rfrac_simple(n,d);
    pari_err(talker,"incompatible variables in gred");
  }
  if (varncmp(varn(d), varn(n)) < 0) return gred_rfrac_simple(n,d);
  if (varncmp(varn(d), varn(n)) > 0) return RgX_Rg_div(n,d);

  /* now n and d are t_POLs in the same variable */
  v = polvaluation(n, &n) - polvaluation(d, &d);
  if (!degpol(d))
  {
    n = RgX_Rg_div(n,gel(d,2));
    return v? RgX_mulXn(n,v): n;
  }

  /* X does not divide gcd(n,d), deg(d) > 0 */
  if (!isinexact(n) && !isinexact(d))
  {
    y = RgX_divrem(n, d, &z);
    if (!signe(z)) return v? RgX_mulXn(y, v): y;
    z = srgcd(d, z);
    if (degpol(z)) { n = gdeuc(n,z); d = gdeuc(d,z); }
  }
  return fix_rfrac(gred_rfrac_simple(n,d), v);
}

GEN
gred_rfrac2(GEN x1, GEN x2)
{
  pari_sp av = avma;
  return gerepileupto(av, gred_rfrac2_i(x1, x2));
}

/* x1,x2 t_INT, return x1/x2 in reduced form */
GEN
gred_frac2(GEN x1, GEN x2)
{
  GEN p1, y = dvmdii(x1,x2,&p1);
  pari_sp av;

  if (p1 == gen_0) return y; /* gen_0 intended */
  av = avma;
  p1 = gcdii(x2,p1);
  if (is_pm1(p1))
  {
    avma = av; y = cgetg(3,t_FRAC);
    gel(y,1) = icopy(x1);
    gel(y,2) = icopy(x2);
  }
  else
  {
    p1 = gclone(p1);
    avma = av; y = cgetg(3,t_FRAC);
    gel(y,1) = diviiexact(x1,p1);
    gel(y,2) = diviiexact(x2,p1);
    gunclone(p1);
  }
  fix_frac(y); return y;
}

/********************************************************************/
/**                                                                **/
/**                          SUBTRACTION                           **/
/**                                                                **/
/********************************************************************/

GEN
gsub(GEN x, GEN y)
{
  pari_sp tetpil, av = avma;
  y = gneg_i(y); tetpil = avma;
  return gerepile(av,tetpil, gadd(x,y));
}

/********************************************************************/
/**                                                                **/
/**                           ADDITION                             **/
/**                                                                **/
/********************************************************************/
/* x, y compatible PADIC */
static GEN
addpp(GEN x, GEN y)
{
  pari_sp av = avma;
  long c,d,e,r,rx,ry;
  GEN u, z, mod, p = gel(x,2);

  (void)new_chunk(5 + lgefint(x[3]) + lgefint(y[3]));
  e = valp(x);
  r = valp(y); d = r-e;
  if (d < 0) { swap(x,y); e = r; d = -d; }
  rx = precp(x);
  ry = precp(y);
  if (d) /* v(x) < v(y) */
  {
    r = d+ry; z = powiu(p,d);
    if (r < rx) mod = mulii(z,gel(y,3)); else { r = rx; mod = gel(x,3); }
    u = addii(gel(x,4), mulii(z,gel(y,4)));
    u = remii(u, mod);
  }
  else
  {
    if (ry < rx) { r=ry; mod = gel(y,3); } else { r=rx; mod = gel(x,3); }
    u = addii(gel(x,4), gel(y,4));
    if (!signe(u) || (c = Z_pvalrem(u,p,&u)) >= r)
    {
      avma = av; return zeropadic(p, e+r);
    }
    if (c)
    {
      mod = diviiexact(mod, powiu(p,c));
      r -= c;
      e += c;
    }
    u = remii(u, mod);
  }
  avma = av; z = cgetg(5,t_PADIC);
  z[1] = evalprecp(r) | evalvalp(e);
  gel(z,2) = icopy(p);
  gel(z,3) = icopy(mod);
  gel(z,4) = icopy(u); return z;
}

/* return x + y, where x is t_INT or t_FRAC(N), y t_PADIC */
static GEN
addQp(GEN x, GEN y)
{
  pari_sp av;
  long tx,vy,py,d,r,e;
  GEN z,q,p,p1,p2,mod,u;

  if (gcmp0(x)) return gcopy(y);

  av = avma; p = gel(y,2); tx = typ(x);
  e = (tx == t_INT)? Z_pvalrem(x,p,&p1)
                   : Z_pvalrem(gel(x,1),p,&p1) -
                     Z_pvalrem(gel(x,2),p,&p2);
  vy = valp(y); d = vy - e; py = precp(y); r = d + py;
  if (r <= 0) { avma = av; return gcopy(y); }
  mod = gel(y,3);
  u   = gel(y,4);
  (void)new_chunk(5 + ((lgefint(mod) + lgefint(p)*labs(d)) << 1));

  if (d > 0)
  {
    q = powiu(p,d);
    mod = mulii(mod, q);
    u   = mulii(u, q);
    if (tx != t_INT && !is_pm1(p2)) p1 = mulii(p1, Fp_inv(p2,mod));
    u = addii(u, p1);
  }
  else if (d < 0)
  {
    q = powiu(p,-d);
    if (tx != t_INT && !is_pm1(p2)) p1 = mulii(p1, Fp_inv(p2,mod));
    p1 = mulii(p1, q);
    u = addii(u, p1);
    r = py; e = vy;
  }
  else
  {
    long c;
    if (tx != t_INT && !is_pm1(p2)) p1 = mulii(p1, Fp_inv(p2,mod));
    u = addii(u, p1);
    if (!signe(u) || (c = Z_pvalrem(u,p,&u)) >= r)
    {
      avma = av; return zeropadic(p,e+r);
    }
    if (c)
    {
      mod = diviiexact(mod, powiu(p,c));
      r -= c;
      e += c;
    }
  }
  u = modii(u, mod);
  avma = av; z = cgetg(5,t_PADIC);
  z[1] = evalprecp(r) | evalvalp(e);
  gel(z,2) = icopy(p);
  gel(z,3) = icopy(mod);
  gel(z,4) = icopy(u); return z;
}

/* Mod(x,X) + Mod(y,X) */
#define add_polmod_same add_polmod_scal
/* Mod(x,X) + Mod(y,Y) */
static GEN
add_polmod(GEN X, GEN Y, GEN x, GEN y)
{
  long T[3] = { evaltyp(t_POLMOD) | _evallg(3),0,0 };
  GEN z = cgetg(3,t_POLMOD);
  long vx = varn(X), vy = varn(Y);
  if (vx==vy) {
    pari_sp av;
    gel(z,1) = srgcd(X,Y); av = avma;
    gel(z,2) = gerepileupto(av, gmod(gadd(x, y), gel(z,1))); return z;
  }
  if (varncmp(vx, vy) < 0)
  { gel(z,1) = gcopy(X); gel(T,1) = Y; gel(T,2) = y; y = T; }
  else
  { gel(z,1) = gcopy(Y); gel(T,1) = X; gel(T,2) = x; x = T; }
  gel(z,2) = gadd(x, y); return z;
}
/* Mod(y, Y) + x,  assuming x scalar or polynomial in same var and reduced degree */
static GEN
add_polmod_scal(GEN Y, GEN y, GEN x)
{
  GEN z = cgetg(3,t_POLMOD);
  gel(z,1) = gcopy(Y);
  gel(z,2) = gadd(x, y); return z;
}

/* check y[a..b-1] and set signe to 1 if one coeff is non-0, 0 otherwise
 * For t_POL and t_SER */
static GEN
NORMALIZE_i(GEN y, long a, long b)
{
  long i;
  for (i = a; i < b; i++)
    if (!gcmp0(gel(y,i))) { setsigne(y, 1); return y; }
  setsigne(y, 0); return y;
}
/* typ(y) == t_POL, varn(y) = vy, x "scalar" [e.g object in lower variable] */
static GEN
add_pol_scal(GEN y, GEN x, long vy)
{
  long i, ly = lg(y);
  GEN z;
  if (ly <= 3) {
    if (ly == 2) return isexactzero(x)? zeropol(vy): scalarpol(x,vy);
    z = cgetg(3, t_POL); z[1] = y[1];
    gel(z,2) = gadd(x,gel(y,2));
    if (gcmp0(gel(z,2))) {
      if (isexactzero(gel(z,2))) { avma = (pari_sp)(z+3); return zeropol(vy); }
      setsigne(z, 0);
    }
    return z;
  }
  z = cgetg(ly,t_POL); z[1] = y[1];
  gel(z,2) = gadd(x,gel(y,2));
  for (i = 3; i < ly; i++) gel(z,i) = gcopy(gel(y,i));
  return NORMALIZE_i(z, 2, ly);
}
/* typ(y) == t_SER, varn(y) = vy, x "scalar" [e.g object in lower variable]
 * l = valp(y) */
static GEN
add_ser_scal(GEN y, GEN x, long vy, long l)
{
  long i, j, ly;
  pari_sp av;
  GEN z;

  if (isexactzero(x)) return gcopy(y);
  ly = lg(y);
  if (l < 3-ly) return gcopy(y);
  if (l < 0)
  {
    z = cgetg(ly,t_SER); z[1] = y[1];
    for (i = 2; i <= 1-l; i++) gel(z,i) = gcopy(gel(y,i));
    gel(z,i) = gadd(x,gel(y,i)); i++;
    for (     ; i < ly; i++)   gel(z,i) = gcopy(gel(y,i));
    return NORMALIZE_i(z, 2, ly);
  }
  if (l > 0)
  {
    ly += l; y -= l; z = cgetg(ly,t_SER);
    z[1] = evalsigne(1) | evalvalp(0) | evalvarn(vy);
    gel(z,2) = gcopy(x);
    for (i=3; i<=l+1; i++) gel(z,i) = gen_0;
    for (   ; i < ly; i++) gel(z,i) = gcopy(gel(y,i));
    if (gcmp0(x)) return NORMALIZE_i(z, l+2, ly);
    return z;
  }
  /* l = 0, !isexactzero(x) */
  if (ly==2) return zeroser(vy, 0); /* 1 + O(1) --> O(1) */

  av = avma; z = cgetg(ly,t_SER);
  x = gadd(x, gel(y,2));
  if (!isexactzero(x))
  {
    z[1] = evalsigne(1) | evalvalp(0) | evalvarn(vy);
    gel(z,2) = x;
    for (i=3; i<ly; i++) gel(z,i) = gcopy(gel(y,i));
    if (gcmp0(x)) return NORMALIZE_i(z, 3, ly);
    return z;
  }
  avma = av; /* first coeff is 0 */
  i = 3; while (i < ly && isexactzero(gel(y,i))) i++;
  i -= 2; ly -= i; y += i;
  z = cgetg(ly,t_SER); z[1] = evalvalp(i) | evalvarn(vy);
  for (j = 2; j < ly; j++) gel(z,j) = gcopy(gel(y,j));
  return NORMALIZE_i(z, 2, ly);
}
/* typ(y) == RFRAC, x polynomial in same variable or "scalar" */
static GEN
add_rfrac_scal(GEN y, GEN x)
{
  pari_sp av = avma;
  GEN n = gadd(gmul(x, gel(y,2)), gel(y,1));
  return gerepileupto(av, gred_rfrac_simple(n, gel(y,2)));
}

/* x "scalar", ty != t_MAT and non-scalar */
static GEN
add_scal(GEN y, GEN x, long ty, long vy)
{
  long tx;
  switch(ty)
  {
    case t_POL: return add_pol_scal(y, x, vy);
    case t_SER: return add_ser_scal(y, x, vy, valp(y));
    case t_RFRAC: return add_rfrac_scal(y, x);
    case t_VEC: case t_COL:
      tx = typ(x);
      if (!is_matvec_t(tx) && isexactzero(x)) return gcopy(y);
      break;
  }
  pari_err(operf,"+",x,y);
  return NULL; /* not reached */
}

static GEN
addfrac(GEN x, GEN y)
{
  GEN x1 = gel(x,1), x2 = gel(x,2), z = cgetg(3,t_FRAC);
  GEN y1 = gel(y,1), y2 = gel(y,2), p1, r, n, d, delta;

  delta = gcdii(x2,y2);
  if (is_pm1(delta))
  { /* numerator is non-zero */
    gel(z,1) = gerepileuptoint((pari_sp)z, addii(mulii(x1,y2), mulii(y1,x2)));
    gel(z,2) = mulii(x2,y2); return z;
  }
  x2 = diviiexact(x2,delta);
  y2 = diviiexact(y2,delta);
  n = addii(mulii(x1,y2), mulii(y1,x2));
  if (!signe(n)) { avma = (pari_sp)(z+3); return gen_0; }
  d = mulii(x2, y2);
  p1 = dvmdii(n, delta, &r);
  if (r == gen_0)
  {
    if (is_pm1(d)) { avma = (pari_sp)(z+3); return icopy(p1); }
    avma = (pari_sp)z;
    gel(z,2) = icopy(d);
    gel(z,1) = icopy(p1); return z;
  }
  p1 = gcdii(delta, r);
  if (!is_pm1(p1))
  {
    delta = diviiexact(delta, p1);
    n     = diviiexact(n, p1);
  }
  d = mulii(d,delta); avma = (pari_sp)z;
  gel(z,1) = icopy(n);
  gel(z,2) = icopy(d); return z;
}

static GEN
add_rfrac(GEN x, GEN y)
{
  GEN x1 = gel(x,1), x2 = gel(x,2), z = cgetg(3,t_RFRAC);
  GEN y1 = gel(y,1), y2 = gel(y,2), p1, r, n, d, delta;
  pari_sp tetpil;

  delta = ggcd(x2,y2);
  if (gcmp1(delta))
  { /* numerator is non-zero */
    gel(z,1) = gerepileupto((pari_sp)z, gadd(gmul(x1,y2), gmul(y1,x2)));
    gel(z,2) = gmul(x2, y2); return z;
  }
  x2 = gdeuc(x2,delta);
  y2 = gdeuc(y2,delta);
  n = gadd(gmul(x1,y2), gmul(y1,x2));
  if (gcmp0(n)) return gerepileupto((pari_sp)(z+3), n);
  tetpil = avma; d = gmul(x2, y2);
  p1 = poldivrem(n, delta, &r); /* we want gcd(n,delta) */
  if (gcmp0(r))
  {
    if (lg(d) == 3) /* "constant" denominator */
    {
      d = gel(d,2);
           if (gcmp_1(d)) p1 = gneg(p1);
      else if (!gcmp1(d)) p1 = gdiv(p1, d);
      return gerepileupto((pari_sp)(z+3), p1);
    }
    gel(z,1) = p1; gel(z,2) = d;
    gerepilecoeffssp((pari_sp)z,tetpil,z+1,2); return z;
  }
  p1 = ggcd(delta, r);
  if (gcmp1(p1))
  {
    tetpil = avma;
    gel(z,1) = gcopy(n);
  }
  else
  {
    delta = gdeuc(delta, p1);
    tetpil = avma;
    gel(z,1) = gdeuc(n,p1);
  }
  gel(z,2) = gmul(d,delta);
  gerepilecoeffssp((pari_sp)z,tetpil,z+1,2); return z;
}

GEN
gadd(GEN x, GEN y)
{
  long tx = typ(x), ty = typ(y), vx, vy, lx, ly, i, l;
  pari_sp av, tetpil;
  GEN z, p1;

  if (tx == ty) switch(tx) /* shortcut to generic case */
  {
    case t_INT: return addii(x,y);
    case t_REAL: return addrr(x,y);
    case t_INTMOD:  { GEN X = gel(x,1), Y = gel(y,1);
      z = cgetg(3,t_INTMOD);
      if (X==Y || equalii(X,Y))
        return add_intmod_same(z, X, gel(x,2), gel(y,2));
      gel(z,1) = gcdii(X,Y);
      av = avma; p1 = addii(gel(x,2),gel(y,2));
      gel(z,2) = gerepileuptoint(av, remii(p1, gel(z,1))); return z;
    }
    case t_FRAC: return addfrac(x,y);
    case t_COMPLEX: z = cgetg(3,t_COMPLEX);
      gel(z,2) = gadd(gel(x,2),gel(y,2));
      if (isexactzero(gel(z,2)))
      {
        avma = (pari_sp)(z+3);
        return gadd(gel(x,1),gel(y,1));
      }
      gel(z,1) = gadd(gel(x,1),gel(y,1));
      return z;
    case t_PADIC:
      if (!equalii(gel(x,2),gel(y,2))) pari_err(operi,"+",x,y);
      return addpp(x,y);
    case t_QUAD: z = cgetg(4,t_QUAD);
      if (!gequal(gel(x,1),gel(y,1))) pari_err(operi,"+",x,y);
      gel(z,1) = gcopy(gel(x,1));
      gel(z,2) = gadd(gel(x,2),gel(y,2));
      gel(z,3) = gadd(gel(x,3),gel(y,3)); return z;
    case t_POLMOD:
      if (gequal(gel(x,1), gel(y,1)))
        return add_polmod_same(gel(x,1), gel(x,2), gel(y,2));
      return add_polmod(gel(x,1), gel(y,1), gel(x,2), gel(y,2));
    case t_POL:
      vx = varn(x);
      vy = varn(y);
      if (vx != vy) {
        if (varncmp(vx, vy) < 0) return add_pol_scal(x, y, vx);
        else                     return add_pol_scal(y, x, vy);
      }
      /* same variable */
      lx = lg(x);
      ly = lg(y);
      if (lx == ly) {
        z = cgetg(lx, t_POL); z[1] = x[1];
        for (i=2; i < lx; i++) gel(z,i) = gadd(gel(x,i),gel(y,i));
        return normalizepol_i(z, lx);
      }
      if (ly < lx) {
        z = cgetg(lx,t_POL); z[1] = x[1];
        for (i=2; i < ly; i++) gel(z,i) = gadd(gel(x,i),gel(y,i));
        for (   ; i < lx; i++) gel(z,i) = gcopy(gel(x,i));
        if (!signe(x)) z = normalizepol_i(z, lx);
      } else {
        z = cgetg(ly,t_POL); z[1] = y[1];
        for (i=2; i < lx; i++) gel(z,i) = gadd(gel(x,i),gel(y,i));
        for (   ; i < ly; i++) gel(z,i) = gcopy(gel(y,i));
        if (!signe(y)) z = normalizepol_i(z, ly);
      }
      return z;
    case t_SER:
      vx = varn(x);
      vy = varn(y);
      if (vx != vy) {
        if (varncmp(vx, vy) < 0) return add_ser_scal(x, y, vx, valp(x));
        else                     return add_ser_scal(y, x, vy, valp(y));
      }
      l = valp(y) - valp(x);
      if (l < 0) { l = -l; swap(x,y); }
      /* valp(x) <= valp(y) */
      lx = lg(x);
      ly = lg(y) + l; if (lx < ly) ly = lx;
      if (l)
      {
        if (l > lx-2) return gcopy(x);
        z = cgetg(ly,t_SER);
        for (i=2; i<=l+1; i++) gel(z,i) = gcopy(gel(x,i));
        for (   ; i < ly; i++) gel(z,i) = gadd(gel(x,i),gel(y,i-l));
      } else {
        z = cgetg(ly,t_SER);
        for (i=2; i < ly; i++) gel(z,i) = gadd(gel(x,i),gel(y,i));
      }
      z[1] = x[1]; return normalize(z);
    case t_RFRAC:
      vx = gvar(x);
      vy = gvar(y);
      if (vx != vy) {
        if (varncmp(vx, vy) < 0) return add_rfrac_scal(x, y);
        else                     return add_rfrac_scal(y, x);
      }
      return add_rfrac(x,y);
    case t_VEC: case t_COL: case t_MAT:
      ly = lg(y);
      if (ly != lg(x)) pari_err(operi,"+",x,y);
      z = cgetg(ly, ty);
      for (i = 1; i < ly; i++) gel(z,i) = gadd(gel(x,i),gel(y,i));
      return z;

    default: pari_err(operf,"+",x,y);
  }
  /* tx != ty */
  if (tx > ty) { swap(x,y); lswap(tx,ty); }

  if (is_const_t(ty)) switch(tx) /* tx < ty, is_const_t(tx) && is_const_t(ty) */
  {
    case t_INT:
      switch(ty)
      {
        case t_REAL: return addir(x,y);
        case t_INTMOD:
          z = cgetg(3, t_INTMOD);
          return add_intmod_same(z, gel(y,1), gel(y,2), modii(x, gel(y,1)));
        case t_FRAC: z = cgetg(3,t_FRAC);
          gel(z,1) = gerepileuptoint((pari_sp)z, addii(gel(y,1), mulii(gel(y,2),x)));
          gel(z,2) = icopy(gel(y,2)); return z;
        case t_COMPLEX: return addRc(x, y);
        case t_PADIC: return addQp(x,y);
        case t_QUAD: return addRq(x, y);
      }

    case t_REAL:
      switch(ty)
      {
        case t_FRAC:
          if (!signe(y[1])) return rcopy(x);
          if (!signe(x))
          {
            lx = expi(gel(y,1)) - expi(gel(y,2)) - expo(x);
            return lx <= 0? rcopy(x): fractor(y, nbits2prec(lx));
          }
          av=avma; z=addir(gel(y,1),mulir(gel(y,2),x)); tetpil=avma;
          return gerepile(av,tetpil,divri(z,gel(y,2)));
        case t_COMPLEX: return addRc(x, y);
        case t_QUAD: return gcmp0(y)? rcopy(x): addqf(y, x, lg(x));

        default: pari_err(operf,"+",x,y);
      }

    case t_INTMOD:
      switch(ty)
      {
        case t_FRAC: { GEN X = gel(x,1);
          z = cgetg(3, t_INTMOD);
          p1 = modii(mulii(gel(y,1), Fp_inv(gel(y,2),X)), X);
          return add_intmod_same(z, X, p1, gel(x,2));
        }
        case t_COMPLEX: return addRc(x, y);
        case t_PADIC: { GEN X = gel(x,1);
          z = cgetg(3, t_INTMOD);
          return add_intmod_same(z, X, gel(x,2), padic_to_Fp(y, X));
        }
        case t_QUAD: return addRq(x, y);
      }

    case t_FRAC:
      switch (ty)
      {
        case t_COMPLEX: return addRc(x, y);
        case t_PADIC: return addQp(x,y);
        case t_QUAD: return addRq(x, y);
      }

    case t_COMPLEX:
      switch(ty)
      {
        case t_PADIC:
          return Zp_nosquare_m1(gel(y,2))? addRc(y, x): addTp(x, y);
        case t_QUAD:
          lx = precision(x); if (!lx) pari_err(operi,"+",x,y);
          return gcmp0(y)? rcopy(x): addqf(y, x, lx);
      }

    case t_PADIC: /* ty == t_QUAD */
      return (kro_quad(gel(y,1),gel(x,2)) == -1)? addRq(x, y): addTp(y, x);
  }
  /* tx < ty, !is_const_t(y) */
  if (ty == t_MAT) {
    if (is_matvec_t(tx)) pari_err(operf,"+",x,y);
    if (isexactzero(x)) return gcopy(y);
    return gaddmat(x,y);
  }
  if (ty == t_POLMOD) /* is_const_t(tx) in this case */
    return add_polmod_scal(gel(y,1), gel(y,2), x);
  vy = gvar(y);
  if (is_scalar_t(tx))  {
    if (tx == t_POLMOD)
    {
      vx = varn(x[1]);
      if (vx == vy) y = gmod(y, gel(x,1)); /* error if ty == t_SER */
      else
        if (varncmp(vx,vy) > 0) return add_scal(y, x, ty, vy);
      return add_polmod_scal(gel(x,1), gel(x,2), y);
    }
    return add_scal(y, x, ty, vy);
  }
  /* x and y are not scalars, ty != t_MAT */
  vx = gvar(x);
  if (vx != vy) { /* x or y is treated as a scalar */
    if (is_vec_t(tx) || is_vec_t(ty)) pari_err(operf,"+",x,y);
    return (varncmp(vx, vy) < 0)? add_scal(x, y, tx, vx)
                                : add_scal(y, x, ty, vy);
  }
  /* vx = vy */
  switch(tx)
  {
    case t_POL:
      switch (ty)
      {
	case t_SER:
	  if (lg(x) == 2) return gcopy(y);
	  i = lg(y) + valp(y) - polvaluation(x, NULL);
	  if (i < 3) return gcopy(y);

	  p1 = greffe(x,i,0); y = gadd(p1,y);
          free(p1); return y;
	
        case t_RFRAC: return add_rfrac_scal(y, x);
      }
      break;

    case t_SER:
      if (ty == t_RFRAC)
      {
        GEN n, d;
        long vn, vd;
        n = gel(y,1); vn = gval(n, vy);
        d = gel(y,2); vd = polvaluation(d, &d);

	l = lg(x) + valp(x) - (vn - vd);
	if (l < 3) return gcopy(x);

	av = avma; 
        /* take advantage of y = t^n ! */
        if (degpol(d)) {
          y = gdiv(n, greffe(d,l,1));
        } else {
          y = gdiv(n, gel(d,2));
          if (gvar(y) == vy) y = greffe(y,l,1); else y = scalarser(y, vy, l);
        }
        setvalp(y, valp(y) - vd);
        return gerepileupto(av, gadd(y, x));
      }
      break;
  }
  pari_err(operf,"+",x,y);
  return NULL; /* not reached */
}

GEN
gaddsg(long x, GEN y)
{
  long ty = typ(y);
  GEN z;

  switch(ty)
  {
    case t_INT:  return addsi(x,y);
    case t_REAL: return addsr(x,y);
    case t_INTMOD:
      z = cgetg(3, t_INTMOD);
      return add_intmod_same(z, gel(y,1), gel(y,2), modsi(x, gel(y,1)));
    case t_FRAC: z = cgetg(3,t_FRAC);
      gel(z,1) = gerepileuptoint((pari_sp)z, addii(gel(y,1), mulis(gel(y,2),x)));
      gel(z,2) = icopy(gel(y,2)); return z;
    case t_COMPLEX:
      z = cgetg(3, t_COMPLEX);
      gel(z,1) = gaddsg(x, gel(y,1));
      gel(z,2) = gcopy(gel(y,2)); return z;

    default: return gadd(stoi(x), y);
  }
}

/********************************************************************/
/**                                                                **/
/**                        MULTIPLICATION                          **/
/**                                                                **/
/********************************************************************/
static GEN
mul_ser_scal(GEN y, GEN x) {
  long ly, i;
  GEN z;
  if (isexactzero(x)) { long vy = varn(y); return zeropol(vy); }
  ly = lg(y); z = cgetg(ly,t_SER); z[1] = y[1];
  for (i = 2; i < ly; i++) gel(z,i) = gmul(x,gel(y,i));
  return normalize(z);
}
/* (n/d) * x, x "scalar" or polynomial in the same variable as d
 * [n/d a valid RFRAC]  */
static GEN
mul_rfrac_scal(GEN n, GEN d, GEN x)
{
  pari_sp av = avma;
  GEN z;
  switch(typ(x))
  {
    case t_INTMOD: case t_POLMOD:
      n = gmul(n, x);
      d = gmul(d, gmodulo(gen_1, gel(x,1)));
      return gerepileupto(av, gdiv(n,d));
  }
  z = gred_rfrac2_i(x, d);
  n = simplify_i(n);
  if (typ(z) == t_RFRAC) 
    z = gred_rfrac_simple(gmul(gel(z,1), n), gel(z,2));
  else
    z = gmul(z, n);
  return gerepileupto(av, z);
}
static GEN
mul_scal(GEN y, GEN x, long ty)
{
  switch(ty)
  {
    case t_POL: return RgX_Rg_mul(y, x);
    case t_SER: return mul_ser_scal(y, x);
    case t_RFRAC: return mul_rfrac_scal(gel(y,1),gel(y,2), x);
    case t_QFI: case t_QFR:
      if (typ(x) == t_INT && gcmp1(x)) return gcopy(y); /* fall through */
  }
  pari_err(operf,"*",x,y);
  return NULL; /* not reached */
}

static GEN
mul_gen_rfrac(GEN X, GEN Y)
{
  GEN y1 = gel(Y,1), y2 = gel(Y,2);
  long vx = gvar(X), vy = varn(y2);
  return (varncmp(vx, vy) <= 0)? mul_scal(Y, X, typ(Y)):
                                 gred_rfrac_simple(gmul(y1,X), y2);
}
/* (x1/x2) * (y1/y2) */
static GEN
mul_rfrac(GEN x1, GEN x2, GEN y1, GEN y2)
{
  GEN z, X, Y;
  pari_sp av = avma;

  X = gred_rfrac2_i(x1, y2);
  Y = gred_rfrac2_i(y1, x2);
  if (typ(X) == t_RFRAC)
  {
    if (typ(Y) == t_RFRAC) {
      x1 = gel(X,1);
      x2 = gel(X,2);
      y1 = gel(Y,1);
      y2 = gel(Y,2);
      z = gred_rfrac_simple(gmul(x1,y1), gmul(x2,y2));
    } else
      z = mul_gen_rfrac(Y, X);
  }
  else if (typ(Y) == t_RFRAC)
    z = mul_gen_rfrac(X, Y);
  else
    z = gmul(X, Y);
  return gerepileupto(av, z);
}
/* (x1/x2) /y2 */
static GEN
div_rfrac_pol(GEN x1, GEN x2, GEN y2)
{
  GEN z, X;
  pari_sp av = avma;

  X = gred_rfrac2_i(x1, y2);
  if (typ(X) == t_RFRAC)
     z = gred_rfrac_simple(gel(X,1), gmul(gel(X,2),x2));
  else
    z = mul_gen_rfrac(X, mkrfrac(gen_1, x2));
  return gerepileupto(av, z);
}

/* Mod(y, Y) * x,  assuming x scalar */
static GEN
mul_polmod_scal(GEN Y, GEN y, GEN x)
{
  GEN z = cgetg(3,t_POLMOD);
  gel(z,1) = gcopy(Y);
  gel(z,2) = gmul(x,y); return z;
}
/* Mod(x,X) * Mod(y,X) */
static GEN
mul_polmod_same(GEN X, GEN x, GEN y)
{
  GEN t, z = cgetg(3,t_POLMOD);
  pari_sp av;
  long v;
  gel(z,1) = gcopy(X); av = avma;
  t = gmul(x, y);
  /* gmod(t, gel(z,1))) optimised */
  if (typ(t) == t_POL  && (v = varn(X)) == varn(t) && lg(t) >= lg(X))
    gel(z,2) = gerepileupto(av, RgX_divrem(t, X, ONLY_REM));
  else
    gel(z,2) = t;
  return z;
}
/* Mod(x,X) * Mod(y,Y) */
static GEN
mul_polmod(GEN X, GEN Y, GEN x, GEN y)
{
  long T[3] = { evaltyp(t_POLMOD) | _evallg(3),0,0 };
  long vx = varn(X), vy = varn(Y);
  GEN z = cgetg(3,t_POLMOD);
  pari_sp av;

  if (vx==vy) {
    gel(z,1) = srgcd(X,Y); av = avma;
    gel(z,2) = gerepileupto(av, gmod(gmul(x, y), gel(z,1))); return z;
  }
  if (varncmp(vx, vy) < 0)
  { gel(z,1) = gcopy(X); gel(T,1) = Y; gel(T,2) = y; y = T; }
  else
  { gel(z,1) = gcopy(Y); gel(T,1) = X; gel(T,2) = x; x = T; }
  gel(z,2) = gmul(x, y); return z;
}

/* compatible t_VEC * t_COL, l = lg(x) = lg(y) */
static GEN
VC_mul(GEN x, GEN y, long l)
{
  pari_sp av = avma;
  GEN z = gen_0;
  long i;
  for (i=1; i<l; i++)
  {
    GEN c = gel(y,i);
    if (!isexactzeroscalar(c)) z = gadd(z, gmul(gel(x,i), c));
  }
  return gerepileupto(av,z);
}
/* compatible t_MAT * t_COL, l = lg(x) = lg(y), lz = l>1? lg(x[1]): 1 */
static GEN
MC_mul(GEN x, GEN y, long l, long lz)
{
  GEN z = cgetg(lz,t_COL);
  long i, j;
  for (i=1; i<lz; i++)
  {
    pari_sp av = avma;
    GEN t = gen_0;
    for (j=1; j<l; j++)
    {
      GEN c = gel(y,j);
      if (!isexactzeroscalar(c)) t = gadd(t, gmul(gcoeff(x,i,j), c));
    }
    gel(z,i) = gerepileupto(av,t);
  }
  return z;
}
/* x,y COMPLEX */
static GEN
mulcc(GEN x, GEN y)
{
  GEN xr = gel(x,1), xi = gel(x,2);
  GEN yr = gel(y,1), yi = gel(y,2);
  GEN p1, p2, p3, p4, z = cgetg(3,t_COMPLEX);
  pari_sp tetpil, av = avma;
#if 1 /* 3M */
  p1 = gmul(xr,yr);
  p2 = gmul(xi,yi); p2 = gneg(p2);
  p3 = gmul(gadd(xr,xi), gadd(yr,yi));
  p4 = gadd(p2, gneg(p1));
#else /* standard product */
  p1 = gmul(xr,yr);
  p2 = gmul(xi,yi); p2 = gneg(p2);
  p3 = gmul(xr,yi);
  p4 = gmul(xi,yr);
#endif
  tetpil = avma;
  gel(z,1) = gadd(p1,p2);
  gel(z,2) = gadd(p3,p4);
  if (isexactzero(gel(z,2)))
  {
    cgiv(gel(z,2));
    return gerepileupto((pari_sp)(z+3), gel(z,1));
  }
  gerepilecoeffssp(av,tetpil, z+1,2); return z;
}
/* x,y PADIC */
static GEN
mulpp(GEN x, GEN y) {
  long l = valp(x) + valp(y);
  pari_sp av;
  GEN z, t;
  if (!equalii(gel(x,2),gel(y,2))) pari_err(operi,"*",x,y);
  if (!signe(x[4])) return zeropadic(gel(x,2), l);
  if (!signe(y[4])) return zeropadic(gel(x,2), l);

  t = (precp(x) > precp(y))? y: x;
  z = cgetp(t); setvalp(z,l); av = avma;
  affii(remii(mulii(gel(x,4),gel(y,4)), gel(t,3)), gel(z,4));
  avma = av; return z;
}
/* x,y QUAD */
static GEN
mulqq(GEN x, GEN y) {
  GEN p1,p2,p3,p4, z = cgetg(4,t_QUAD);
  pari_sp av, tetpil;
  p1 = gel(x,1);
  if (!gequal(p1, gel(y,1))) pari_err(operi,"*",x,y);

  gel(z,1) = gcopy(p1); av = avma;
  p2 = gmul(gel(x,2),gel(y,2));
  p3 = gmul(gel(x,3),gel(y,3));
  p4 = gmul(gneg_i(gel(p1,2)),p3);

  if (gcmp0(gel(p1,3)))
  {
    tetpil = avma;
    gel(z,2) = gerepile(av,tetpil,gadd(p4,p2));
    av = avma;
    p2 = gmul(gel(x,2),gel(y,3));
    p3 = gmul(gel(x,3),gel(y,2)); tetpil = avma;
    gel(z,3) = gerepile(av,tetpil,gadd(p2,p3)); return z;
  }

  p1 = gadd(gmul(gel(x,2),gel(y,3)), gmul(gel(x,3),gel(y,2)));
  tetpil = avma;
  gel(z,2) = gadd(p2,p4);
  gel(z,3) = gadd(p1,p3);
  gerepilecoeffssp(av,tetpil,z+2,2); return z;
}

GEN
mulcxI(GEN x)
{
  GEN z;
  switch(typ(x))
  {
    case t_INT: case t_REAL: case t_FRAC:
      return mkcomplex(gen_0, x);
    case t_COMPLEX:
      if (isexactzero(gel(x,1))) return gneg(gel(x,2));
      z = cgetg(3,t_COMPLEX);
      gel(z,1) = gneg(gel(x,2));
      z[2] = x[1]; return z;
    default:
      return gmul(gi, x);
  }
}
GEN
mulcxmI(GEN x)
{
  GEN z;
  switch(typ(x))
  {
    case t_INT: case t_REAL: case t_FRAC:
      return mkcomplex(gen_0, gneg(x));
    case t_COMPLEX:
      if (isexactzero(gel(x,1))) return gel(x,2);
      z = cgetg(3,t_COMPLEX);
      z[1] = x[2];
      gel(z,2) = gneg(gel(x,1)); return z;
    default:
      return gmul(mkcomplex(gen_0, gen_m1), x);
  }
}

GEN
gmul(GEN x, GEN y)
{
  long tx, ty, lx, ly, vx, vy, i, j, l;
  pari_sp av, tetpil;
  GEN z, p1, p2;

  if (x == y) return gsqr(x);
  tx = typ(x); ty = typ(y);
  if (tx == ty) switch(tx)
  {
    case t_INT: return mulii(x,y);
    case t_REAL: return mulrr(x,y);
    case t_INTMOD: { GEN X = gel(x,1), Y = gel(y,1);
      z = cgetg(3,t_INTMOD); 
      if (X==Y || equalii(X,Y))
        return mul_intmod_same(z, X, gel(x,2), gel(y,2));
      gel(z,1) = gcdii(X,Y); av = avma; p1 = mulii(gel(x,2),gel(y,2));
      gel(z,2) = gerepileuptoint(av, remii(p1, gel(z,1))); return z;
    }
    case t_FRAC:
    {
      GEN x1 = gel(x,1), x2 = gel(x,2);
      GEN y1 = gel(y,1), y2 = gel(y,2);
      z=cgetg(3,t_FRAC);
      p1 = gcdii(x1, y2);
      if (!is_pm1(p1)) { x1 = diviiexact(x1,p1); y2 = diviiexact(y2,p1); }
      p1 = gcdii(x2, y1);
      if (!is_pm1(p1)) { x2 = diviiexact(x2,p1); y1 = diviiexact(y1,p1); }
      tetpil = avma;
      gel(z,2) = mulii(x2,y2);
      gel(z,1) = mulii(x1,y1);
      fix_frac_if_int_GC(z,tetpil); return z;
    }
    case t_COMPLEX: return mulcc(x, y);
    case t_PADIC: return mulpp(x, y);
    case t_QUAD: return mulqq(x, y);
    case t_POLMOD:
      if (gequal(gel(x,1), gel(y,1)))
        return mul_polmod_same(gel(x,1), gel(x,2), gel(y,2));
      return mul_polmod(gel(x,1), gel(y,1), gel(x,2), gel(y,2));
    case t_POL:
      vx = varn(x);
      vy = varn(y);
      if (vx != vy) {
        if (varncmp(vx, vy) < 0) return RgX_Rg_mul(x, y);
        else                     return RgX_Rg_mul(y, x);
      }
      return RgX_mul(x, y);

    case t_SER: {
      long mix, miy;
      vx = varn(x);
      vy = varn(y);
      if (vx != vy) {
        if (varncmp(vx, vy) < 0) return mul_ser_scal(x, y);
        else                     return mul_ser_scal(y, x);
      }
      lx = lg(x); if (lx > lg(y)) { lx = lg(y); swap(x, y); }
      if (lx == 2) return zeroser(vx, valp(x)+valp(y));
      z = cgetg(lx,t_SER);
      z[1] = evalvalp(valp(x)+valp(y)) | evalvarn(vx) | evalsigne(1);
      if (lx > 200) /* threshold for 32bit coeffs: 400, 512 bits: 100 */
      {
        long ly;
        y = RgX_mul(ser2pol_i(x, lx), ser2pol_i(y, lx));
        ly = lg(y);
        if (ly >= lx) {
          for (i = 2; i < lx; i++) z[i] = y[i];
        } else {
          for (i = 2; i < ly; i++) z[i] = y[i];
          for (     ; i < lx; i++) gel(z,i) = gen_0;
        }
        z = normalize(z);
        return gerepilecopy((pari_sp)(z + lx), z);
      }
      x += 2; y += 2; z += 2; lx -= 3;
      p2 = (GEN)gpmalloc((lx+1)*sizeof(long));
      mix = miy = 0;
      for (i=0; i<=lx; i++)
      {
        p2[i] = !isexactzero(gel(y,i)); if (p2[i]) miy = i;
        if (!isexactzero(gel(x,i))) mix = i;
        p1 = gen_0; av = avma;
        for (j=i-mix; j<=min(i,miy); j++)
          if (p2[j]) p1 = gadd(p1, gmul(gel(y,j),gel(x,i-j)));
        gel(z,i) = gerepileupto(av,p1);
      }
      z -= 2; /* back to normalcy */
      free(p2); return normalize(z);
    }
    case t_QFI: return compimag(x,y);
    case t_QFR: return compreal(x,y);
    case t_RFRAC: return mul_rfrac(gel(x,1),gel(x,2), gel(y,1),gel(y,2));
    case t_MAT:
      ly = lg(y); if (ly == 1) return cgetg(1,t_MAT);
      lx = lg(x);
      if (lx != lg(y[1])) pari_err(operi,"*",x,y);
      z = cgetg(ly,t_MAT);
      l = (lx == 1)? 1: lg(x[1]);
      for (j=1; j<ly; j++) gel(z,j) = MC_mul(x, gel(y,j), lx, l);
      return z;

    case t_VECSMALL: /* multiply as permutation. cf perm_mul */
      l = lg(x); z = cgetg(l, t_VECSMALL);
      if (l != lg(y)) break;
      for (i=1; i<l; i++)
      {
        long yi = y[i];
        if (yi < 1 || yi >= l) pari_err(operf,"*",x,y);
        z[i] = x[yi];
      }
      return z;


    default:
      pari_err(operf,"*",x,y);
  }
  /* tx != ty */
  if (is_const_t(ty) && is_const_t(tx))  {
    if (tx > ty) { swap(x,y); lswap(tx,ty); }
    switch(tx) {
    case t_INT:
      switch(ty)
      {
        case t_REAL: return mulir(x,y);
        case t_INTMOD:
          z = cgetg(3, t_INTMOD);
          return mul_intmod_same(z, gel(y,1), gel(y,2), modii(x, gel(y,1)));
        case t_FRAC:
          if (!signe(x)) return gen_0;
          z=cgetg(3,t_FRAC);
          p1 = gcdii(x,gel(y,2));
          if (is_pm1(p1))
          {
            avma = (pari_sp)z;
            gel(z,2) = icopy(gel(y,2));
            gel(z,1) = mulii(gel(y,1), x);
          }
          else
          {
            x = diviiexact(x,p1); tetpil = avma;
            gel(z,2) = diviiexact(gel(y,2), p1);
            gel(z,1) = mulii(gel(y,1), x);
            fix_frac_if_int_GC(z,tetpil);
          }
          return z;
        case t_COMPLEX: return mulRc(x, y);
        case t_PADIC: return signe(x)? mulTp(x, y): gen_0;
        case t_QUAD: return mulRq(x,y);
      }

    case t_REAL:
      switch(ty)
      {
        case t_FRAC:
          if (!signe(y[1])) return gen_0;
          av = avma;
          return gerepileuptoleaf(av, divri(mulri(x,gel(y,1)), gel(y,2)));
        case t_COMPLEX: return mulRc(x, y);
        case t_QUAD: return mulqf(y, x, lg(x));
        default: pari_err(operf,"*",x,y);
      }

    case t_INTMOD:
      switch(ty)
      {
        case t_FRAC: { GEN X = gel(x,1);
          z = cgetg(3, t_INTMOD); p1 = modii(mulii(gel(y,1), gel(x,2)), X);
          return div_intmod_same(z, X, p1, remii(gel(y,2), X)); 
        }
        case t_COMPLEX: return mulRc(x, y);
        case t_PADIC: { GEN X = gel(x,1);
          z = cgetg(3, t_INTMOD);
          return mul_intmod_same(z, X, gel(x,2), padic_to_Fp(y, X));
        }
        case t_QUAD: return mulRq(x, y);
      }

    case t_FRAC:
      switch(ty)
      {
        case t_COMPLEX: return mulRc(x, y);
        case t_PADIC: return signe(x[1])? mulTp(x, y): gen_0;
        case t_QUAD: return mulRq(x, y);
      }

    case t_COMPLEX:
      switch(ty)
      {
        case t_PADIC:
          return Zp_nosquare_m1(gel(y,2))? mulRc(y, x): mulTp(x, y);
        case t_QUAD:
          lx = precision(x); if (!lx) pari_err(operi,"*",x,y);
          return mulqf(y, x, lx);
      }

    case t_PADIC: /* ty == t_QUAD */
      return (kro_quad(gel(y,1),gel(x,2))== -1)? mulRq(x, y): mulTp(y, x);
    }
  }

  if (is_matvec_t(ty))
  {
    ly = lg(y);
    if (!is_matvec_t(tx))
    {
      if (is_noncalc_t(tx)) pari_err(operf, "*",x,y); /* necessary if ly = 1 */
      z = cgetg(ly,ty);
      for (i=1; i<ly; i++) gel(z,i) = gmul(x,gel(y,i));
      return z;
    }
    lx = lg(x);

    switch(tx)
    {
      case t_VEC:
        switch(ty)
        {
          case t_COL:
            if (lx != ly) pari_err(operi,"*",x,y);
            return VC_mul(x, y, lx);

          case t_MAT:
            if (ly == 1) return cgetg(1,t_VEC);
            if (lx != lg(y[1])) pari_err(operi,"*",x,y);
            z = cgetg(ly, t_VEC);
            for (i=1; i<ly; i++) gel(z,i) = VC_mul(x, gel(y,i), lx);
            return z;
        }
        break;

      case t_COL:
        switch(ty)
        {
          case t_VEC:
            z = cgetg(ly,t_MAT);
            for (j=1; j<ly; j++)
            {
              GEN c = cgetg(lx, t_COL); gel(z,j) = c;
              for (i=1; i<lx; i++) gel(c,i) = gmul(gel(x,i), gel(y,j));
            }
            return z;

          case t_MAT:
            if (ly != 1 && lg(y[1]) != 2) pari_err(operi,"*",x,y);
            z = cgetg(ly,t_MAT);
            for (i=1; i<ly; i++) gel(z,i) = gmul(gcoeff(y,1,i),x);
            return z;
        }
        break;

      case t_MAT:
        switch(ty)
        {
          case t_VEC:
            if (lx != 2) pari_err(operi,"*",x,y);
            return gmul(gel(x,1), y);

          case t_COL:
            if (lx != ly) pari_err(operi,"*",x,y);
            return MC_mul(x, y, lx, (lx == 1)? 1: lg(x[1]));
        }
    }
  }
  if (is_matvec_t(tx))
  {
    if (is_noncalc_t(ty)) pari_err(operf, "*",x,y); /* necessary if lx = 1 */
    lx = lg(x); z = cgetg(lx,tx);
    for (i=1; i<lx; i++) gel(z,i) = gmul(y,gel(x,i));
    return z;
  }
  if (tx > ty) { swap(x,y); lswap(tx,ty); }
  /* tx < ty, !ismatvec(x and y) */

  if (ty == t_POLMOD) /* is_const_t(tx) in this case */
    return mul_polmod_scal(gel(y,1), gel(y,2), x);
  if is_scalar_t(tx) {
    if (tx == t_POLMOD) {
      vx = varn(x[1]);
      vy = gvar(y);
      if (vx != vy) {
        if (varncmp(vx,vy) > 0) return mul_scal(y, x, ty);
        return mul_polmod_scal(gel(x,1), gel(x,2), y);
      }
      /* error if ty == t_SER */
      av = avma; y = gmod(y, gel(x,1));
      return gerepileupto(av, mul_polmod_same(gel(x,1), gel(x,2), y));
    }
    return mul_scal(y, x, ty);
  }

  /* x and y are not scalars, nor matvec */
  vx = gvar(x);
  vy = gvar(y);
  if (vx != vy) { /* x or y is treated as a scalar */
    return (varncmp(vx, vy) < 0)? mul_scal(x, y, tx)
                                : mul_scal(y, x, ty);
  }
  /* vx = vy */
  switch(tx)
  {
    case t_POL:
      switch (ty)
      {
	case t_SER:
        {
          long vn;
          if (lg(x) == 2) return zeropol(vx);
          if (lg(y) == 2) return zeroser(vx, valp(y)+polvaluation(x,NULL));
          av = avma;
          vn = polvaluation(x, &x);
          avma = av;
          /* take advantage of x = t^n ! */
          if (degpol(x)) {
            p1 = greffe(x,lg(y),0);
            p2 = gmul(p1,y); free(p1);
          } else
            p2 = mul_ser_scal(y, gel(x,2));
          setvalp(p2, valp(p2) + vn);
          return p2;
        }
	
        case t_RFRAC: return mul_rfrac_scal(gel(y,1),gel(y,2), x);
      }
      break;
	
    case t_SER:
      switch (ty)
      {
	case t_RFRAC:
	  av = avma;
          return gerepileupto(av, gdiv(gmul(gel(y,1),x), gel(y,2)));
      }
      break;
  }
  pari_err(operf,"*",x,y);
  return NULL; /* not reached */
}

int
ff_poltype(GEN *x, GEN *p, GEN *pol)
{
  GEN Q, P = *x, pr,p1,p2,y;
  long i, lx;

  if (!signe(P)) return 0;
  lx = lg(P); Q = *pol;
  for (i=2; i<lx; i++)
  {
    p1 = gel(P,i); if (typ(p1) != t_POLMOD) {Q=NULL;break;}
    p2 = gel(p1,1);
    if (Q==NULL) { Q = p2; if (degpol(Q) <= 0) return 0; }
    else if (p2 != Q)
    {
      if (!gequal(p2, Q))
      {
        if (DEBUGMEM) pari_warn(warner,"different modulus in ff_poltype");
        return 0;
      }
      if (DEBUGMEM > 2) pari_warn(warner,"different pointers in ff_poltype");
    }
  }
  if (Q) {
    *x = P = to_Kronecker(P, Q);
    *pol = Q; lx = lg(P);
  }
  pr = *p; y = cgetg(lx, t_POL);
  for (i=lx-1; i>1; i--)
  {
    p1 = gel(P,i);
    switch(typ(p1))
    {
      case t_INTMOD: break;
      case t_INT:
        if (*p) p1 = modii(p1, *p);
        gel(y,i) = p1; continue;
      default:
        return (Q && !pr)? 1: 0;
    }
    p2 = gel(p1,1);
    if (pr==NULL) pr = p2;
    else if (p2 != pr)
    {
      if (!equalii(p2, pr))
      {
        if (DEBUGMEM) pari_warn(warner,"different modulus in ff_poltype");
        return 0;
      }
      if (DEBUGMEM > 2) pari_warn(warner,"different pointers in ff_poltype");
    }
    y[i] = p1[2];
  }
  y[1] = P[1];
  *x = y; *p = pr; return (Q || pr);
}

GEN
gsqr(GEN x)
{
  long tx=typ(x), lx, i, j, l;
  pari_sp av, tetpil;
  GEN z, p1, p2, p3, p4;

  if (is_scalar_t(tx))
    switch(tx)
    {
      case t_INT: return sqri(x);
      case t_REAL: return mulrr(x,x);
      case t_INTMOD: { GEN X = gel(x,1);
        z = cgetg(3,t_INTMOD);
        gel(z,2) = gerepileuptoint((pari_sp)z, remii(sqri(gel(x,2)), X));
        gel(z,1) = icopy(X); return z;
      }
      case t_FRAC: z=cgetg(3,t_FRAC);
	gel(z,1) = sqri(gel(x,1));
	gel(z,2) = sqri(gel(x,2)); return z;

      case t_COMPLEX:
        if (isexactzero(gel(x,1))) {
          av = avma;
          return gerepileupto(av, gneg(gsqr(gel(x,2))));
        }
	z = cgetg(3,t_COMPLEX); av = avma;
	p1 = gadd(gel(x,1),gel(x,2));
	p2 = gadd(gel(x,1), gneg_i(gel(x,2)));
	p3 = gmul(gel(x,1),gel(x,2));
	tetpil = avma;
	gel(z,1) = gmul(p1,p2);
        gel(z,2) = gshift(p3,1); gerepilecoeffssp(av,tetpil,z+1,2); return z;
	
      case t_PADIC:
	z = cgetg(5,t_PADIC);
	i = (equaliu(gel(x,2), 2) && signe(x[4]))? 1: 0;
        if (i && precp(x) == 1) i = 2; /* (1 + O(2))^2 = 1 + O(2^3) */
        z[1] = evalprecp(precp(x)+i) | evalvalp(valp(x) << 1);
        gel(z,2) = icopy(gel(x,2));
        gel(z,3) = shifti(gel(x,3), i); av = avma;
	gel(z,4) = gerepileuptoint(av, remii(sqri(gel(x,4)), gel(z,3)));
	return z;
	
      case t_QUAD: z = cgetg(4,t_QUAD);
	p1 = gel(x,1);
        gel(z,1) = gcopy(p1); av = avma;
	p2 = gsqr(gel(x,2));
        p3 = gsqr(gel(x,3));
	p4 = gmul(gneg_i(gel(p1,2)),p3);

	if (gcmp0(gel(p1,3)))
	{
	  tetpil = avma;
	  gel(z,2) = gerepile(av,tetpil,gadd(p4,p2));
	  av = avma;
          p2 = gmul(gel(x,2),gel(x,3)); tetpil = avma;
	  gel(z,3) = gerepile(av,tetpil,gmul2n(p2,1)); return z;
	}

	p1 = gmul2n(gmul(gel(x,2),gel(x,3)), 1);
        tetpil = avma;
	gel(z,2) = gadd(p2,p4);
        gel(z,3) = gadd(p1,p3);
	gerepilecoeffssp(av,tetpil,z+2,2); return z;

      case t_POLMOD:
        z=cgetg(3,t_POLMOD); gel(z,1) = gcopy(gel(x,1));
	av=avma; p1=gsqr(gel(x,2)); tetpil=avma;
        gel(z,2) = gerepile(av,tetpil, grem(p1,gel(z,1)));
	return z;
    }

  switch(tx)
  {
    case t_POL:
    {
      GEN a = x, p = NULL, pol = NULL;
      long vx = varn(x);
      av = avma;
      if (ff_poltype(&x,&p,&pol))
      {
        z = FpX_sqr(x, p);
        if (p) z = FpX_to_mod(z,p);
        if (pol) z = from_Kronecker(z,pol);
        z = gerepileupto(av, z);
      }
      else { avma = av; z = RgX_sqr(a); }
      setvarn(z, vx); return z;
    }

    case t_SER:
    {
      long mi;
      lx = lg(x);
      if (lx == 2) return zeroser(varn(x), 2*valp(x));
      z = cgetg(lx, t_SER);
      z[1] = evalvalp(2*valp(x)) | evalvarn(varn(x));
      x += 2; z += 2; lx -= 3;
      p2 = (GEN)gpmalloc((lx+1)*sizeof(long));
      mi = 0;
      for (i=0; i<=lx; i++)
      {
	p2[i] = !isexactzero(gel(x,i)); if (p2[i]) mi = i;
        p1=gen_0; av=avma; l=((i+1)>>1) - 1;
        for (j=i-mi; j<=min(l,mi); j++)
          if (p2[j] && p2[i-j]) p1 = gadd(p1, gmul(gel(x,j),gel(x,i-j)));
        p1 = gshift(p1,1);
        if ((i&1) == 0 && p2[i>>1])
          p1 = gadd(p1, gsqr((GEN)x[i>>1]));
        gel(z,i) = gerepileupto(av,p1);
      }
      z -= 2; free(p2); return normalize(z);
    }
    case t_RFRAC: z = cgetg(3,t_RFRAC);
      gel(z,1) = gsqr(gel(x,1));
      gel(z,2) = gsqr(gel(x,2)); return z;

    case t_MAT:
      lx = lg(x);
      if (lx!=1 && lx != lg(x[1])) pari_err(operi,"*",x,x);
      z = cgetg(lx, t_MAT);
      for (j=1; j<lx; j++) gel(z,j) = MC_mul(x, gel(x,j), lx, lx);
      return z;

    case t_QFR: return sqcompreal(x);
    case t_QFI: return sqcompimag(x);
    case t_VECSMALL:
      l = lg(x); z = cgetg(l, t_VECSMALL);
      for (i=1; i<l; i++)
      {
        long xi = x[i];
        if (xi < 1 || xi >= l) pari_err(operf,"*",x,x);
        z[i] = x[xi];
      }
      return z;
  }
  pari_err(operf,"*",x,x);
  return NULL; /* not reached */
}

/********************************************************************/
/**                                                                **/
/**                           DIVISION                             **/
/**                                                                **/
/********************************************************************/
static GEN
div_rfrac_scal(GEN x, GEN y)
{ 
  pari_sp av = avma;
  return gerepileupto(av, gred_rfrac_simple(gel(x,1),gmul(y, gel(x,2))));
}
static GEN
div_scal_rfrac(GEN x, GEN y)
{ 
  GEN y1 = gel(y,1), y2 = gel(y,2);
  pari_sp av = avma;
  if (typ(y1) == t_POL && varn(y2) == varn(y1) && degpol(y1) > 0)
    return gerepileupto(av, gred_rfrac_simple(gmul(x, y2), y1));
  return RgX_Rg_mul(y2, gdiv(x,y1));
}
static GEN
div_rfrac(GEN x, GEN y)
{ return mul_rfrac(gel(x,1),gel(x,2), gel(y,2),gel(y,1)); }

static GEN
div_ser_scal(GEN x, GEN y) {
  long i, lx = lg(x);
  GEN z = cgetg_copy(lx, x); z[1] = x[1];
  for (i=2; i<lx; i++) gel(z,i) = gdiv(gel(x,i),y);
  return normalize(z);
}
static GEN
div_T_scal(GEN x, GEN y, long tx) {
  switch(tx)
  {
    case t_POL: return RgX_Rg_div(x, y);
    case t_SER: return div_ser_scal(x, y);
    case t_RFRAC: return div_rfrac_scal(x,y);
  }
  pari_err(operf,"/",x,y);
  return NULL; /* not reached */
}

static GEN
div_scal_pol(GEN x, GEN y) {
  long ly = lg(y);
  pari_sp av;
  if (ly == 3) return gdiv(x,gel(y,2));
  if (isexactzero(x)) return zeropol(varn(y));
  av = avma;
  return gerepileupto(av, gred_rfrac_simple(x,y));
}
static GEN
div_scal_ser(GEN x, GEN y) { /* TODO: improve */
  GEN z;
  long ly, i;
  if (gcmp0(x)) { pari_sp av=avma; return gerepileupto(av, gmul(x, ginv(y))); }
  ly = lg(y); z = (GEN)gpmalloc(ly*sizeof(long));
  z[0] = evaltyp(t_SER) | evallg(ly);
  z[1] = evalsigne(1) | evalvalp(0) | evalvarn(varn(y));
  gel(z,2) = x; for (i=3; i<ly; i++) gel(z,i) = gen_0;
  y = gdiv(z,y); free(z); return y;
}
static GEN
div_scal_T(GEN x, GEN y, long ty) {
  switch(ty)
  {
    case t_POL: return div_scal_pol(x, y);
    case t_SER: return div_scal_ser(x, y);
    case t_RFRAC: return div_scal_rfrac(x, y);
  }
  pari_err(operf,"/",x,y);
  return NULL; /* not reached */
}

/* assume tx = ty = t_SER, same variable vx */
static GEN
div_ser(GEN x, GEN y, long vx)
{
  long i, j, l = valp(x) - valp(y), lx = lg(x), ly = lg(y);
  GEN y_lead, p1, p2, z;
  pari_sp av;

  if (!signe(y)) pari_err(gdiver);
  if (lx == 2) return zeroser(vx, l);
  y_lead = gel(y,2);
  if (gcmp0(y_lead)) /* normalize denominator if leading term is 0 */
  {
    pari_warn(warner,"normalizing a series with 0 leading term");
    for (i=3,y++; i<ly; i++,y++)
    {
      y_lead = gel(y,2); ly--; l--;
      if (!gcmp0(y_lead)) break;
    }
  }
  if (ly < lx) lx = ly;
  p2 = (GEN)gpmalloc(lx*sizeof(long));
  for (i=3; i<lx; i++)
  {
    p1 = gel(y,i);
    if (isexactzero(p1)) p2[i] = 0;
    else
    {
      av = avma; gel(p2,i) = gclone(gneg_i(p1));
      avma = av;
    }
  }
  z = cgetg(lx,t_SER);
  z[1] = evalvalp(l) | evalvarn(vx) | evalsigne(1);
  gel(z,2) = gdiv(gel(x,2), y_lead);
  for (i=3; i<lx; i++)
  {
    av = avma; p1 = gel(x,i);
    for (j=2; j<i; j++)
    {
      l = i-j+2;
      if (p2[l]) p1 = gadd(p1, gmul(gel(z,j), gel(p2,l)));
    }
    gel(z,i) = gerepileupto(av, gdiv(p1, y_lead));
  }
  for (i=3; i<lx; i++)
    if (p2[i]) gunclone(gel(p2,i));
  free(p2); return normalize(z);
}
/* x,y compatible PADIC */
static GEN
divpp(GEN x, GEN y) {
  pari_sp av;
  long a, b;
  GEN z, M;

  if (!signe(x[4])) return zeropadic(gel(x,2), valp(x)-valp(y));
  a = precp(x);
  b = precp(y); if (a > b) { M = gel(y,3); } else { M = gel(x,3); b = a; }
  z = cgetg(5, t_PADIC);
  z[1] = evalprecp(b) | evalvalp(valp(x) - valp(y));
  gel(z,2) = icopy(gel(x,2));
  gel(z,3) = icopy(M); av = avma;
  gel(z,4) = gerepileuptoint(av, remii(mulii(gel(x,4), Fp_inv(gel(y,4), M)), M) );
  return z;
}

GEN
gdiv(GEN x, GEN y)
{
  long tx = typ(x), ty = typ(y), lx, ly, vx, vy, i;
  pari_sp av, tetpil;
  GEN z, p1, p2;

  if (tx == ty) switch(tx)
  {
    case t_INT:
      if (is_pm1(y)) return (signe(y) < 0)? negi(x): icopy(x);
      if (is_pm1(x)) {
        long s = signe(y);
        if (!s) pari_err(gdiver);
        if (signe(x) < 0) s = -s;
        z = cgetg(3, t_FRAC);
        gel(z,1) = s<0? gen_m1: gen_1;
        gel(z,2) = absi(y); return z;
      }
      return gred_frac2(x,y);

    case t_REAL: return divrr(x,y);
    case t_INTMOD: { GEN X = gel(x,1), Y = gel(y,1);
      z = cgetg(3,t_INTMOD);
      if (X==Y || equalii(X,Y))
        return div_intmod_same(z, X, gel(x,2), gel(y,2));
      gel(z,1) = gcdii(X,Y); av = avma;
      p1 = mulii(gel(x,2), Fp_inv(gel(y,2), gel(z,1)));
      gel(z,2) = gerepileuptoint(av, remii(p1, gel(z,1))); return z;
    }
    case t_FRAC: {
      GEN x1 = gel(x,1), x2 = gel(x,2);
      GEN y1 = gel(y,1), y2 = gel(y,2);
      z = cgetg(3, t_FRAC);
      p1 = gcdii(x1, y1);
      if (!is_pm1(p1)) { x1 = diviiexact(x1,p1); y1 = diviiexact(y1,p1); }
      p1 = gcdii(x2, y2);
      if (!is_pm1(p1)) { x2 = diviiexact(x2,p1); y2 = diviiexact(y2,p1); }
      tetpil = avma;
      gel(z,2) = mulii(x2,y1);
      gel(z,1) = mulii(x1,y2);
      fix_frac(z);
      fix_frac_if_int_GC(z,tetpil);
      return z;
    }
    case t_COMPLEX:
      av=avma; p1 = cxnorm(y); p2 = mulcc(x, gconj(y)); tetpil = avma;
      return gerepile(av, tetpil, gdiv(p2,p1));

    case t_PADIC:
      if (!equalii(gel(x,2),gel(y,2))) pari_err(operi,"/",x,y);
      return divpp(x, y);

    case t_QUAD:
      if (!gequal(gel(x,1),gel(y,1))) pari_err(operi,"/",x,y);
      av = avma; p1 = quadnorm(y); p2 = mulqq(x, gconj(y)); tetpil = avma;
      return gerepile(av, tetpil, gdiv(p2,p1));

    case t_POLMOD: av = avma;
      if (gequal(gel(x,1), gel(y,1)))
      {
        GEN X = gel(x,1);
        x = gel(x,2);
        y = gel(y,2);
        if (degpol(X) == 2) { /* optimized for speed */
          z = mul_polmod_same(X, x, quad_polmod_conj(y, X));
          return gdiv(z, quad_polmod_norm(y, X));
        }
        y = ginvmod(y, X);
        z = mul_polmod_same(X, x, y);
      } else z = gmul(x, ginv(y));
      return gerepileupto(av, z);

    case t_POL:
      vx = varn(x);
      vy = varn(y);
      if (vx != vy) {
        if (varncmp(vx, vy) < 0) return RgX_Rg_div(x, y);
                            else return div_scal_pol(x, y);
      }
      if (!signe(y)) pari_err(gdiver);
      if (lg(y) == 3) return RgX_Rg_div(x,gel(y,2));
      if (isexactzero(x)) return zeropol(vy);
      return gred_rfrac2(x,y);

    case t_SER:
      vx = varn(x);
      vy = varn(y);
      if (vx != vy) {
        if (varncmp(vx, vy) < 0) return div_ser_scal(x, y);
                            else return div_scal_ser(x, y);
      }
      return div_ser(x, y, vx);
    case t_RFRAC:
      vx = gvar(x);
      vy = gvar(y);
      if (vx != vy) {
        if (varncmp(vx, vy) < 0) return div_rfrac_scal(x, y);
                            else return div_scal_rfrac(x, y);
      }
      return div_rfrac(x,y);

    case t_QFI: av = avma; return gerepileupto(av, compimag(x, ginv(y)));
    case t_QFR: av = avma; return gerepileupto(av, compreal(x, ginv(y)));

    case t_MAT:
      av = avma;
      return gerepileupto(av, gmul(x, invmat(y)));

    default: pari_err(operf,"/",x,y);
  }

  if (tx==t_INT && is_const_t(ty)) /* optimized for speed */
  {
    long s = signe(x);
    if (!s) {
      if (gcmp0(y)) pari_err(gdiver);
      if (ty != t_INTMOD) return gen_0;
      z = cgetg(3,t_INTMOD);
      gel(z,1) = icopy(gel(y,1));
      gel(z,2) = gen_0; return z;
    }
    if (is_pm1(x)) {
      if (s > 0) return ginv(y);
      av = avma; return gerepileupto(av, ginv(gneg(y)));
    }
    switch(ty)
    {
      case t_REAL: return divir(x,y);
      case t_INTMOD:
        z = cgetg(3, t_INTMOD);
        return div_intmod_same(z, gel(y,1), modii(x, gel(y,1)), gel(y,2));
      case t_FRAC:
        z = cgetg(3,t_FRAC); p1 = gcdii(x,gel(y,1));
        if (is_pm1(p1))
        {
          avma = (pari_sp)z;
          gel(z,2) = icopy(gel(y,1));
          gel(z,1) = mulii(gel(y,2), x);
          fix_frac(z);
          fix_frac_if_int(z);
        }
        else
        {
          x = diviiexact(x,p1); tetpil = avma;
          gel(z,2) = diviiexact(gel(y,1), p1);
          gel(z,1) = mulii(gel(y,2), x);
          fix_frac(z);
          fix_frac_if_int_GC(z,tetpil);
        }
        return z;

      case t_COMPLEX: return divRc(x,y);
      case t_PADIC: return divTp(x, y);
      case t_QUAD:
        av = avma; p1 = quadnorm(y); p2 = mulRq(x, gconj(y)); tetpil = avma;
        return gerepile(av, tetpil, gdiv(p2,p1));
    }
  }
  if (gcmp0(y) && ty != t_MAT) pari_err(gdiver);

  if (is_const_t(tx) && is_const_t(ty)) switch(tx)
  {
    case t_REAL:
      switch(ty)
      {
        case t_INT: return divri(x,y);
        case t_FRAC:
          av = avma; z = divri(mulri(x,gel(y,2)), gel(y,1));
          return gerepileuptoleaf(av, z);
        case t_COMPLEX: return divRc(x, y);
        case t_QUAD: return divfq(x, y, lg(x));
        default: pari_err(operf,"/",x,y);
      }

    case t_INTMOD:
      switch(ty)
      {
        case t_INT:
          z = cgetg(3, t_INTMOD);
          return div_intmod_same(z, gel(x,1), gel(x,2), modii(y, gel(x,1)));
        case t_FRAC: { GEN X = gel(x,1);
          z = cgetg(3,t_INTMOD); p1 = remii(mulii(gel(y,2), gel(x,2)), X);
          return div_intmod_same(z, X, p1, modii(gel(y,1), X));
        }
        case t_COMPLEX:
          av = avma; p1 = cxnorm(y); p2 = mulRc(x, gconj(y)); tetpil = avma;
          return gerepile(av,tetpil, gdiv(p2,p1));

        case t_QUAD:
          av = avma; p1 = quadnorm(y); p2 = gmul(x,gconj(y)); tetpil = avma;
          return gerepile(av,tetpil, gdiv(p2,p1));

        case t_PADIC: { GEN X = gel(x,1);
          z = cgetg(3, t_INTMOD);
          return div_intmod_same(z, X, gel(x,2), padic_to_Fp(y, X));
        }
        case t_REAL: pari_err(operf,"/",x,y);
      }

    case t_FRAC:
      switch(ty)
      {
        case t_INT: z = cgetg(3, t_FRAC);
        p1 = gcdii(y,gel(x,1));
        if (is_pm1(p1))
        {
          avma = (pari_sp)z; tetpil = 0;
          gel(z,1) = icopy(gel(x,1));
        }
        else
        {
          y = diviiexact(y,p1); tetpil = avma;
          gel(z,1) = diviiexact(gel(x,1), p1);
        }
        gel(z,2) = mulii(gel(x,2),y);
        fix_frac(z);
        if (tetpil) fix_frac_if_int_GC(z,tetpil);
        return z;

        case t_REAL:
          av=avma; p1=mulri(y,gel(x,2)); tetpil=avma;
          return gerepile(av, tetpil, divir(gel(x,1), p1));

        case t_INTMOD: { GEN Y = gel(y,1);
          z = cgetg(3,t_INTMOD); p1 = remii(mulii(gel(y,2),gel(x,2)), Y);
          return div_intmod_same(z, Y, gel(x,1), p1);
        }
        case t_COMPLEX: return divRc(x, y);

        case t_PADIC:
          if (!signe(x[1])) return gen_0;
          return divTp(x, y);

        case t_QUAD:
          av=avma; p1=quadnorm(y); p2=gmul(x,gconj(y)); tetpil=avma;
          return gerepile(av,tetpil,gdiv(p2,p1));
      }

    case t_COMPLEX:
      switch(ty)
      {
        case t_INT: case t_REAL: case t_INTMOD: case t_FRAC: return divcR(x,y);
        case t_PADIC:
          return Zp_nosquare_m1(gel(y,2))? divcR(x,y): divTp(x, y);
        case t_QUAD:
          lx = precision(x); if (!lx) pari_err(operi,"/",x,y);
          return divfq(x, y, lx);
      }

    case t_PADIC:
      switch(ty)
      {
        case t_INT: case t_FRAC: { GEN p = gel(x,2);
          return signe(x[4])? divpT(x, y)
                            : zeropadic(p, valp(x) - ggval(y,p));
        }
        case t_INTMOD: { GEN Y = gel(y,1);
          z = cgetg(3, t_INTMOD);
          return div_intmod_same(z, Y, padic_to_Fp(x, Y), gel(y,2));
        }
        case t_COMPLEX: case t_QUAD:
          av=avma; p1=gmul(x,gconj(y)); p2=gnorm(y); tetpil=avma;
          return gerepile(av,tetpil,gdiv(p1,p2));

        case t_REAL: pari_err(operf,"/",x,y);
      }

    case t_QUAD:
      switch (ty)
      {
        case t_INT: case t_INTMOD: case t_FRAC:
          z = cgetg(4,t_QUAD);
          gel(z,1) = gcopy(gel(x,1));
          gel(z,2) = gdiv(gel(x,2), y);
          gel(z,3) = gdiv(gel(x,3), y); return z;
        case t_REAL: return divqf(x, y, lg(y));
        case t_PADIC: return divTp(x, y);
        case t_COMPLEX:
          ly = precision(y); if (!ly) pari_err(operi,"/",x,y);
          return divqf(x, y, ly);
      }
  }
  switch(ty) {
    case t_REAL: case t_INTMOD: case t_PADIC: case t_POLMOD:
      return gmul(x, ginv(y)); /* missing gerepile, for speed */
    case t_MAT:
      av = avma; return gerepileupto(av, gmul(x, invmat(y)));
  }
  if (is_matvec_t(tx)) {
    lx = lg(x); z = cgetg_copy(lx, x);
    for (i=1; i<lx; i++) gel(z,i) = gdiv(gel(x,i),y);
    return z;
  }

  vy = gvar(y);
  if (tx == t_POLMOD) { GEN X = gel(x,1);
    vx = varn(X);
    if (vx != vy) {
      if (varncmp(vx, vy) > 0) return div_scal_T(x, y, ty);
      z = cgetg(3,t_POLMOD);
      gel(z,1) = gcopy(X);
      gel(z,2) = gdiv(gel(x,2), y); return z;
    }
    /* y is POL, SER or RFRAC */
    av = avma; y = ginvmod(gmod(y,X), X);
    return gerepileupto(av, mul_polmod_same(X, gel(x,2), y));
  }
  /* x and y are not both is_scalar_t. If one of them is scalar, it's not a
   * POLMOD (done already), hence its variable is BIGINT. If the other has
   * variable BIGINT, then the operation is incorrect */
  vx = gvar(x);
  if (vx != vy) { /* includes cases where one is scalar */
    if (varncmp(vx, vy) < 0) return div_T_scal(x, y, tx);
                        else return div_scal_T(x, y, ty);
  }
  switch(tx)
  {
    case t_POL:
      switch(ty)
      {
	case t_SER:
          if (lg(y) == 2)
            return zeroser(vx, polvaluation(x,NULL) - valp(y));
	  p1 = greffe(x,lg(y),0);
          p2 = div_ser(p1, y, vx);
          free(p1); return p2;

        case t_RFRAC:
        {
          GEN y1 = gel(y,1), y2 = gel(y,2);
          if (typ(y1) == t_POL && varn(y1) == vx)
            return mul_rfrac_scal(y2, y1, x);
          av = avma;
          return gerepileupto(av, RgX_Rg_div(RgX_mul(y2, x), y1));
        }
      }
      break;

    case t_SER:
      switch(ty)
      {
	case t_POL:
          if (lg(x) == 2)
            return zeroser(vx, valp(x) - polvaluation(y,NULL));
	  p1 = greffe(y,lg(x),0);
          p2 = div_ser(x, p1, vx);
          free(p1); return p2;
	case t_RFRAC:
	  av = avma;
	  return gerepileupto(av, gdiv(gmul(x,gel(y,2)), gel(y,1)));
      }
      break;

    case t_RFRAC:
      switch(ty)
      {
	case t_POL: return div_rfrac_pol(gel(x,1),gel(x,2), y);
	case t_SER:
	  av = avma;
	  return gerepileupto(av, gdiv(gel(x,1), gmul(gel(x,2),y)));
      }
      break;
  }
  pari_err(operf,"/",x,y);
  return NULL; /* not reached */
}

/********************************************************************/
/**                                                                **/
/**                     SIMPLE MULTIPLICATION                      **/
/**                                                                **/
/********************************************************************/
GEN
gmulsg(long s, GEN y)
{
  long ty = typ(y), ly, i;
  pari_sp av;
  GEN z;

  switch(ty)
  {
    case t_INT:  return mulsi(s,y);
    case t_REAL: return mulsr(s,y);
    case t_INTMOD: { GEN p = gel(y,1);
      z = cgetg(3,t_INTMOD);
      gel(z,2) = gerepileuptoint((pari_sp)z, modii(mulsi(s,gel(y,2)), p));
      gel(z,1) = icopy(p); return z;
    }
    case t_FRAC:
      if (!s) return gen_0;
      z = cgetg(3,t_FRAC);
      i = cgcd(s, smodis(gel(y,2), s));
      if (i == 1)
      {
        gel(z,2) = icopy(gel(y,2));
        gel(z,1) = mulis(gel(y,1), s);
      }
      else
      {
        gel(z,2) = divis(gel(y,2), i);
        gel(z,1) = mulis(gel(y,1), s/i);
        fix_frac_if_int(z);
      }
      return z;

    case t_COMPLEX: z = cgetg(3, t_COMPLEX);
      gel(z,1) = gmulsg(s,gel(y,1));
      gel(z,2) = gmulsg(s,gel(y,2)); return z;

    case t_PADIC:
      if (!s) return gen_0;
      av = avma; return gerepileupto(av, mulpp(cvtop2(stoi(s),y), y));

    case t_QUAD: z = cgetg(4, t_QUAD);
      gel(z,1) = gcopy(gel(y,1));
      gel(z,2) = gmulsg(s,gel(y,2));
      gel(z,3) = gmulsg(s,gel(y,3)); return z;

    case t_POLMOD: z = cgetg(3, t_POLMOD);
      gel(z,1) = gcopy(gel(y,1));
      gel(z,2) = gmulsg(s,gel(y,2)); return z;

    case t_POL:
      if (!s || !signe(y)) return zeropol(varn(y));
      ly = lg(y); z = cgetg(ly,t_POL); z[1]=y[1];
      for (i=2; i<ly; i++) gel(z,i) = gmulsg(s,gel(y,i));
      return normalizepol_i(z, ly);

    case t_SER:
      if (!s) return zeropol(varn(y));
      ly = lg(y); z = cgetg(ly,t_SER); z[1] = y[1];
      for (i=2; i<ly; i++) gel(z,i) = gmulsg(s,gel(y,i));
      return normalize(z);

    case t_RFRAC:
      if (!s) return zeropol(gvar(y));
      z = cgetg(3, t_RFRAC);
      i = itos( ggcd(stoi(s),gel(y,2)) );
      avma = (pari_sp)z;
      if (i == 1)
      {
        gel(z,1) = gmulgs(gel(y,1), s);
        gel(z,2) = gcopy(gel(y,2));
      }
      else
      {
        gel(z,1) = gmulgs(gel(y,1), s/i);
        gel(z,2) = gdivgs(gel(y,2), i);
      }
      return z;

    case t_VEC: case t_COL: case t_MAT:
      ly = lg(y); z = cgetg(ly,ty);
      for (i=1; i<ly; i++) gel(z,i) = gmulsg(s,gel(y,i));
      return z;
  }
  pari_err(typeer,"gmulsg");
  return NULL; /* not reached */
}

/********************************************************************/
/**                                                                **/
/**                       SIMPLE DIVISION                          **/
/**                                                                **/
/********************************************************************/

GEN
gdivgs(GEN x, long s)
{
  long tx = typ(x), lx, i;
  pari_sp av;
  GEN z, y, p1;

  if (!s) pari_err(gdiver);
  switch(tx)
  {
    case t_INT:
      av = avma; z = divis_rem(x,s,&i);
      if (!i) return z;

      i = cgcd(s, i);
      avma=av; z = cgetg(3,t_FRAC);
      if (i == 1)
        y = icopy(x);
      else
      {
        s /= i; y = diviuexact(x, i);
        if (signe(x) < 0) setsigne(y, -1);
      }
      gel(z,1) = y;
      gel(z,2) = stoi(s); fix_frac(z); return z;

    case t_REAL:
      return divrs(x,s);

    case t_INTMOD:
      z = cgetg(3, t_INTMOD);
      return div_intmod_same(z, gel(x,1), gel(x,2), modsi(s, gel(x,1)));

    case t_FRAC: z = cgetg(3, t_FRAC);
      i = cgcd(s, smodis(gel(x,1), s));
      if (i == 1)
      {
        gel(z,2) = mulsi(s, gel(x,2));
        gel(z,1) = icopy(gel(x,1));
      }
      else
      {
        gel(z,2) = mulsi(s/i, gel(x,2));
        gel(z,1) = divis(gel(x,1), i);
      }
      fix_frac(z);
      fix_frac_if_int(z); return z;

    case t_COMPLEX: z = cgetg(3, t_COMPLEX);
      gel(z,1) = gdivgs(gel(x,1),s);
      gel(z,2) = gdivgs(gel(x,2),s); return z;

    case t_PADIC: return gdiv(x, stoi(s));

    case t_QUAD: z = cgetg(4, t_QUAD);
      gel(z,1) = gcopy(gel(x,1));
      gel(z,2) = gdivgs(gel(x,2),s);
      gel(z,3) = gdivgs(gel(x,3),s); return z;

    case t_POLMOD: z = cgetg(3, t_POLMOD);
      gel(z,1) = gcopy(gel(x,1));
      gel(z,2) = gdivgs(gel(x,2),s); return z;

    case t_RFRAC:
      av = avma;
      p1 = ggcd(stoi(s),gel(x,1));
      if (typ(p1) == t_INT)
      {
        avma = av;
        z = cgetg(3, t_RFRAC);
        i = p1[2];
        if (i == 1)
        {
          gel(z,1) = gcopy(gel(x,1));
          gel(z,2) = gmulsg(s,gel(x,2));
        }
        else
        {
          gel(z,1) = gdivgs(gel(x,1), i);
          gel(z,2) = gmulgs(gel(x,2), s/i);
        }
      }
      else /* t_FRAC */
      {
        z = cgetg(3, t_RFRAC);
        gel(z,1) = gdiv(gel(x,1), p1);
        gel(z,2) = gmul(gel(x,2), gdivsg(s,p1));
        z = gerepilecopy(av, z);
      }
      return z;

    case t_POL: case t_SER: case t_VEC: case t_COL: case t_MAT:
      z = init_gen_op(x, tx, &lx, &i);
      for (; i<lx; i++) gel(z,i) = gdivgs(gel(x,i),s);
      return z;
    
  }
  pari_err(operf,"/",x, stoi(s));
  return NULL; /* not reached */
}

/* True shift (exact multiplication by 2^n) */
GEN
gmul2n(GEN x, long n)
{
  long tx = typ(x), lx, i, k, l;
  GEN z, a, b;

  switch(tx)
  {
    case t_INT:
      if (n>=0) return shifti(x,n);
      if (!signe(x)) return gen_0;
      l = vali(x); n = -n;
      if (n<=l) return shifti(x,-n);
      z = cgetg(3,t_FRAC);
      gel(z,1) = shifti(x,-l);
      gel(z,2) = int2n(n-l); return z;
	
    case t_REAL:
      return shiftr(x,n);

    case t_INTMOD: b = gel(x,1); a = gel(x,2);
      z = cgetg(3,t_INTMOD);
      if (n <= 0) return div_intmod_same(z, b, a, modii(int2n(-n), b));
      gel(z,2) = gerepileuptoint((pari_sp)z, modii(shifti(a,n), b));
      gel(z,1) = icopy(b); return z;

    case t_FRAC: a = gel(x,1); b = gel(x,2);
      l = vali(a);
      k = vali(b);
      if (n+l >= k)
      {
        if (expi(b) == k) return shifti(a,n-k); /* b power of 2 */
        l = n-k; k = -k;
      }
      else
      {
        k = -(l+n); l = -l;
      }
      z = cgetg(3,t_FRAC);
      gel(z,1) = shifti(a,l);
      gel(z,2) = shifti(b,k); return z;

    case t_QUAD: z = cgetg(4,t_QUAD);
      gel(z,1) = gcopy(gel(x,1));
      gel(z,2) = gmul2n(gel(x,2),n);
      gel(z,3) = gmul2n(gel(x,3),n); return z;

    case t_POLMOD: z = cgetg(3,t_POLMOD);
      gel(z,1) = gcopy(gel(x,1));
      gel(z,2) = gmul2n(gel(x,2),n); return z;

    case t_POL:
      z = init_gen_op(x, tx, &lx, &i);
      for (; i<lx; i++) gel(z,i) = gmul2n(gel(x,i),n);
      return normalizepol_i(z, lx); /* needed if char = 2 */
    case t_SER:
      z = init_gen_op(x, tx, &lx, &i);
      for (; i<lx; i++) gel(z,i) = gmul2n(gel(x,i),n);
      return normalize(z); /* needed if char = 2 */
    case t_COMPLEX: case t_VEC: case t_COL: case t_MAT:
      z = init_gen_op(x, tx, &lx, &i);
      for (; i<lx; i++) gel(z,i) = gmul2n(gel(x,i),n);
      return z;

    case t_RFRAC: /* int2n wrong if n < 0 */
      return mul_rfrac_scal(gel(x,1),gel(x,2), gmul2n(gen_1,n));

    case t_PADIC: /* int2n wrong if n < 0 */
      return gmul(gmul2n(gen_1,n),x);
  }
  pari_err(typeer,"gmul2n");
  return NULL; /* not reached */
}
