/* $Id: rootpol.c 7857 2006-04-11 17:28:55Z kb $

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
/*                ROOTS OF COMPLEX POLYNOMIALS                     */
/*  (original code contributed by Xavier Gourdon, INRIA RR 1852)   */
/*                                                                 */
/*******************************************************************/
#include "pari.h"
#include "paripriv.h"

static const double pariINFINITY = 100000.;

/********************************************************************/
/**                                                                **/
/**                   FAST ARITHMETIC over Z[i]                    **/
/**                                                                **/
/********************************************************************/
static long KARASQUARE_LIMIT, COOKSQUARE_LIMIT;

/* fast sum of x,y: t_INT or t_COMPLEX(t_INT) */
static GEN
addCC(GEN x, GEN y)
{
  GEN z;

  if (typ(x) == t_INT)
  {
    if (typ(y) == t_INT) return addii(x,y);
    /* ty == t_COMPLEX */
    z = cgetg(3,t_COMPLEX);
    gel(z,1) = addii(x, gel(y,1));
    gel(z,2) = icopy(gel(y,2)); return z;
  }
  /* tx == t_COMPLEX */
  z = cgetg(3,t_COMPLEX);
  if (typ(y) == t_INT)
  {
    gel(z,1) = addii(gel(x,1),y);
    gel(z,2) = icopy(gel(x,2)); return z;
  }
  /* ty == t_COMPLEX */
  gel(z,1) = addii(gel(x,1),gel(y,1));
  gel(z,2) = addii(gel(x,2),gel(y,2)); return z;
}
/* fast product of x,y: t_INT or t_COMPLEX(t_INT) */
static GEN
mulCC(GEN x, GEN y)
{
  GEN z;

  if (typ(x) == t_INT)
  {
    if (typ(y) == t_INT) return mulii(x,y);
    /* ty == t_COMPLEX */
    z = cgetg(3,t_COMPLEX);
    gel(z,1) = mulii(x, gel(y,1));
    gel(z,2) = mulii(x, gel(y,2)); return z;
  }
  /* tx == t_COMPLEX */
  z = cgetg(3,t_COMPLEX);
  if (typ(y) == t_INT)
  {
    gel(z,1) = mulii(gel(x,1),y);
    gel(z,2) = mulii(gel(x,2),y); return z;
  }
  /* ty == t_COMPLEX */
  {
    pari_sp av = avma, tetpil;
    GEN p1, p2;

    p1 = mulii(gel(x,1),gel(y,1));
    p2 = mulii(gel(x,2),gel(y,2));
    y = mulii(addii(gel(x,1),gel(x,2)),
              addii(gel(y,1),gel(y,2)));
    x = addii(p1,p2); tetpil = avma;
    gel(z,1) = subii(p1,p2);
    gel(z,2) = subii(y,x); gerepilecoeffssp(av,tetpil,z+1,2);
    return z;
  }
}
/* fast squaring x: t_INT or t_COMPLEX(t_INT) */
static GEN
sqrCC(GEN x)
{
  GEN z;

  if (typ(x) == t_INT) return sqri(x);
  /* tx == t_COMPLEX */
  z = cgetg(3,t_COMPLEX);
  {
    pari_sp av = avma, tetpil;
    GEN y, p1, p2;

    p1 = sqri(gel(x,1));
    p2 = sqri(gel(x,2));
    y = sqri(addii(gel(x,1),gel(x,2)));
    x = addii(p1,p2); tetpil = avma;
    gel(z,1) = subii(p1,p2);
    gel(z,2) = subii(y,x); gerepilecoeffssp(av,tetpil,z+1,2);
    return z;
  }
}

static void
set_karasquare_limit(long bit)
{
  if (bit<600)       { KARASQUARE_LIMIT=8; COOKSQUARE_LIMIT=400; }
  else if (bit<2000) { KARASQUARE_LIMIT=4; COOKSQUARE_LIMIT=200; }
  else if (bit<3000) { KARASQUARE_LIMIT=4; COOKSQUARE_LIMIT=125; }
  else if (bit<5000) { KARASQUARE_LIMIT=2; COOKSQUARE_LIMIT= 75; }
  else               { KARASQUARE_LIMIT=1; COOKSQUARE_LIMIT= 50; }
}

/* assume lP > 0, lP = lgpol(P) */
static GEN
CX_square_spec(GEN P, long lP)
{
  GEN s, t;
  long i, j, l, nn, n = lP - 1;
  pari_sp av;

  nn = n<<1; s = cgetg(nn+3,t_POL); s[1] = evalsigne(1)|evalvarn(0);
  gel(s,2) = sqrCC(gel(P,0)); /* i = 0 */
  for (i=1; i<=n; i++)
  {
    av = avma; l = (i+1)>>1;
    t = mulCC(gel(P,0), gel(P,i)); /* j = 0 */
    for (j=1; j<l; j++) t = addCC(t, mulCC(gel(P,j), gel(P,i-j)));
    t = gmul2n(t,1);
    if ((i & 1) == 0) t = addCC(t, sqrCC((GEN)P[i>>1]));
    gel(s,i+2) = gerepileupto(av, t);
  }
  gel(s,nn+2) = sqrCC(gel(P,n)); /* i = nn */
  for (   ; i<nn; i++)
  {
    av = avma; l = (i+1)>>1;
    t = mulCC(gel(P,i-n),gel(P,n)); /* j = i-n */
    for (j=i-n+1; j<l; j++) t = addCC(t, mulCC(gel(P,j),gel(P,i-j)));
    t = gmul2n(t,1);
    if ((i & 1) == 0) t = addCC(t, sqrCC((GEN)P[i>>1]));
    gel(s,i+2) = gerepileupto(av, t);
  }
  return normalizepol_i(s, nn+3);
}
/* not stack clean */
static GEN
RgX_addspec(GEN x, long nx, GEN y, long ny)
{
  GEN z, t;
  long i;
  if (nx == ny) {
    z = cgetg(nx+2,t_POL); z[1] = evalsigne(1)|evalvarn(0); t = z+2;
    for (i=0; i < nx; i++) gel(t,i) = gadd(gel(x,i),gel(y,i));
    return normalizepol_i(z, nx+2);
  }
  if (ny < nx) {
    z = cgetg(nx+2,t_POL); z[1] = evalsigne(1)|evalvarn(0); t = z+2;
    for (i=0; i < ny; i++) gel(t,i) = gadd(gel(x,i),gel(y,i));
    for (   ; i < nx; i++) t[i] = x[i];
    return normalizepol_i(z, nx+2);
  } else {
    z = cgetg(ny+2,t_POL); z[1] = evalsigne(1)|evalvarn(0); t = z+2;
    for (i=0; i < nx; i++) gel(t,i) = gadd(gel(x,i),gel(y,i));
    for (   ; i < ny; i++) t[i] = y[i];
    return normalizepol_i(z, ny+2);
  }
}
/* nx = lgpol(x) */
static GEN
RgX_s_mulspec(GEN x, long nx, long s)
{
  GEN z, t;
  long i;
  if (!s || !nx) return zeropol(0);
  z = cgetg(nx+2, t_POL); z[1] = evalsigne(1)|evalvarn(0); t = z + 2;
  for (i=0; i < nx; i++) gel(t,i) = gmulgs(gel(x,i), s);
  return z;
}
/* nx = lgpol(x), return x << s. Inefficient if s = 0... */
static GEN
RgX_shiftspec(GEN x, long nx, long s)
{
  GEN z, t;
  long i;
  if (!nx) return zeropol(0);
  z = cgetg(nx+2, t_POL); z[1] = evalsigne(1)|evalvarn(0); t = z + 2;
  for (i=0; i < nx; i++) gel(t,i) = gmul2n(gel(x,i), s);
  return z;
}

/* spec function. nP = lgpol(P) */
static GEN
karasquare(GEN P, long nP)
{
  GEN Q, s0, s1, s2, a, t;
  long n0, n1, i, l, N, N0, N1, n = nP - 1; /* degree(P) */
  pari_sp av;

  if (n <= KARASQUARE_LIMIT) return nP? CX_square_spec(P, nP): zeropol(0);
  av = avma;
  n0 = (n>>1) + 1; n1 = nP - n0;
  s0 = karasquare(P, n0); Q = P + n0;
  s2 = karasquare(Q, n1);
  s1 = RgX_addspec(P, n0, Q, n1);
  s1 = gadd(karasquare(s1+2, lgpol(s1)), gneg(gadd(s0,s2)));
  N = (n<<1) + 1;
  a = cgetg(N + 2, t_POL); a[1] = evalsigne(1)|evalvarn(0);
  t = a+2; l = lgpol(s0); s0 += 2; N0 = n0<<1;
  for (i=0; i < l;  i++) t[i] = s0[i];
  for (   ; i < N0; i++) gel(t,i) = gen_0;
  t = a+2 + N0; l = lgpol(s2); s2 += 2; N1 = N - N0;
  for (i=0; i < l;  i++) t[i] = s2[i];
  for (   ; i < N1; i++) gel(t,i) = gen_0;
  t = a+2 + n0; l = lgpol(s1); s1 += 2;
  for (i=0; i < l; i++)  gel(t,i) = gadd(gel(t,i), gel(s1,i));
  return gerepilecopy(av, normalizepol_i(a, N+2));
}
/* spec function. nP = lgpol(P) */
static GEN
cook_square(GEN P, long nP)
{
  GEN Q, p0, p1, p2, p3, q, r, t, vp, vm;
  long n0, n3, i, j, n = nP - 1;
  pari_sp av;

  if (n <= COOKSQUARE_LIMIT) return  nP? karasquare(P, nP): zeropol(0);
  av = avma;

  n0 = (n+1) >> 2; n3 = n+1 - 3*n0;
  p0 = P;
  p1 = p0+n0;
  p2 = p1+n0;
  p3 = p2+n0; /* lgpol(p0,p1,p2) = n0, lgpol(p3) = n3 */

  q = cgetg(8,t_VEC) + 4;
  Q = cook_square(p0, n0);
  r = RgX_addspec(p0,n0, p2,n0);
  t = RgX_addspec(p1,n0, p3,n3);
  gel(q,-1) = gadd(r,gneg(t));
  gel(q,1)  = gadd(r,t);
  r = RgX_addspec(p0,n0, RgX_shiftspec(p2,n0, 2)+2,n0);
  t = gmul2n(RgX_addspec(p1,n0, RgX_shiftspec(p3,n3, 2)+2,n3), 1);
  gel(q,-2) = gadd(r,gneg(t));
  gel(q,2)  = gadd(r,t);
  r = RgX_addspec(p0,n0, RgX_s_mulspec(p2,n0, 9)+2,n0);
  t = gmulsg(3, RgX_addspec(p1,n0, RgX_s_mulspec(p3,n3, 9)+2,n3));
  gel(q,-3) = gadd(r,gneg(t));
  gel(q,3)  = gadd(r,t);

  r = new_chunk(7);
  vp = cgetg(4,t_VEC);
  vm = cgetg(4,t_VEC);
  for (i=1; i<=3; i++)
  {
    GEN a = gel(q,i), b = gel(q,-i);
    a = cook_square(a+2, lgpol(a));
    b = cook_square(b+2, lgpol(b));
    gel(vp,i) = gadd(b, a);
    gel(vm,i) = gadd(b, gneg(a));
  }
  gel(r,0) = Q;
  gel(r,1) = gdivgs(gsub(gsub(gmulgs(gel(vm,2),9),gel(vm,3)),
                     gmulgs(gel(vm,1),45)),
                60);
  gel(r,2) = gdivgs(gadd(gadd(gmulgs(gel(vp,1),270),gmulgs(Q,-490)),
                     gadd(gmulgs(gel(vp,2),-27),gmulgs(gel(vp,3),2))),
                360);
  gel(r,3) = gdivgs(gadd(gadd(gmulgs(gel(vm,1),13),gmulgs(gel(vm,2),-8)),
                    gel(vm,3)),
                48);
  gel(r,4) = gdivgs(gadd(gadd(gmulgs(Q,56),gmulgs(gel(vp,1),-39)),
                     gsub(gmulgs(gel(vp,2),12),gel(vp,3))),
                144);
  gel(r,5) = gdivgs(gsub(gadd(gmulgs(gel(vm,1),-5),gmulgs(gel(vm,2),4)),
                     gel(vm,3)),
                240);
  gel(r,6) = gdivgs(gadd(gadd(gmulgs(Q,-20),gmulgs(gel(vp,1),15)),
                     gadd(gmulgs(gel(vp,2),-6),gel(vp,3))),
                720);
  q = cgetg(2*n+3,t_POL); q[1] = evalsigne(1)|evalvarn(0);
  t = q+2;
  for (i=0; i<=2*n; i++) gel(t,i) = gen_0;
  for (i=0; i<=6; i++,t += n0)
  {
    GEN h = gel(r,i);
    long d = lgpol(h);
    h += 2;
    for (j=0; j<d; j++) gel(t,j) = gadd(gel(t,j), gel(h,j));
  }
  return gerepilecopy(av, normalizepol_i(q, 2*n+3));
}

static GEN
graeffe(GEN p)
{
  GEN p0, p1, s0, s1, t;
  long n = degpol(p), n0, n1, i, ns1;

  if (!n) return gcopy(p);
  n0 = (n>>1)+1; n1 = n+1 - n0;
  p0 = new_chunk(n0); for (i=0; i<n0; i++) p0[i] = p[2+(i<<1)];
  p1 = new_chunk(n1); for (i=0; i<n1; i++) p1[i] = p[3+(i<<1)];

  s0 = cook_square(p0, n0);
  s1 = cook_square(p1, n1); ns1 = degpol(s1);
  t = cgetg(ns1+4, t_POL);
  t[1] = evalsigne(1)|evalvarn(0);
  gel(t,2) = gen_0;
  for (i=0; i<=ns1; i++) gel(t,3+i) = gneg(gel(s1,2+i));
  return gadd(s0,t); /* now t contains -x * s1 */
}

/********************************************************************/
/**                                                                **/
/**                       MODULUS OF ROOTS                         **/
/**                                                                **/
/********************************************************************/

static double
log2ir(GEN x)
{
  double l;

  if (!signe(x)) return -pariINFINITY;
  if (typ(x) == t_INT)
  {
    GEN m = int_MSW(x);
    l = (double)(ulong)*m;
    if (lgefint(x)==3) return log2(l);
#ifndef LONG_IS_64BIT /* overkill ? The first word should be enough... */
    l += ((double)(ulong)*int_precW(m)) / 4294967296.; /* 2^32 */
#endif
    return log2(l) + (double)(BITS_IN_LONG*(lgefint(x)-3));
  }
  /* else t_REAL */
  l = (double)(ulong)x[2];
  return log2(l) + (double)(expo(x) - (long)(BITS_IN_LONG-1));
}
/* return log(|x|) */
static double
dblogr(GEN x) {
  double l;
  if (!signe(x)) return -pariINFINITY;
  l = (double)(ulong)x[2];
  return log(l) + LOG2 * (expo(x) - (long)(BITS_IN_LONG-1));
}
static GEN /* beware overflow */
dblexp(double x) { return fabs(x) < 100.? dbltor(exp(x)): mpexp(dbltor(x)); }

double
dbllog2(GEN z)
{
  double x, y;

  if (typ(z) != t_COMPLEX) return log2ir(z);
  x = log2ir(gel(z,1));
  y = log2ir(gel(z,2));
  if (fabs(x-y) > 10) return max(x,y);
  return x + 0.5*log2(1 + exp2(2*(y-x)));
}

/* find s such that  A_h <= 2^s <= 2 A_i  for one h and all i < n = deg(p),
 * with  A_i := (binom(n,i) lc(p) / p_i) ^ 1/(n-i), and  p = sum p_i X^i */
static long
findpower(GEN p)
{
  double x, L, mins = pariINFINITY;
  long n = degpol(p),i;

  L = dbllog2(gel(p,n+2)); /* log2(lc * binom(n,i)) */
  for (i=n-1; i>=0; i--)
  {
    L += log2((double)(i+1) / (double)(n-i));
    x = dbllog2(gel(p,i+2));
    if (x != -pariINFINITY)
    {
      double s = (L - x) / (double)(n-i);
      if (s < mins) mins = s;
    }
  }
  i = (long)ceil(mins);
  if (i - mins > 1 - 1e-12) i--;
  return i;
}

/* returns the exponent for logmodulus(), from the Newton diagram */
static long
newton_polygon(GEN p, long k)
{
  pari_sp av = avma;
  double *logcoef, slope;
  long n = degpol(p), i, j, h, l, *vertex;

  init_dalloc();
  logcoef = (double*)stackmalloc((n+1)*sizeof(double));
  vertex = (long*)new_chunk(n+1);

  /* vertex[i] = 1 if i a vertex of convex hull, 0 otherwise */
  for (i=0; i<=n; i++) { logcoef[i] = dbllog2(gel(p,2+i)); vertex[i] = 0; }
  vertex[0] = 1; /* sentinel */
  for (i=0; i < n; i=h)
  {
    slope = logcoef[i+1]-logcoef[i];
    for (j = h = i+1; j<=n; j++)
    {
      double pij = (logcoef[j]-logcoef[i])/(double)(j-i);
      if (slope < pij) { slope = pij; h = j; }
    }
    vertex[h] = 1;
  }
  h = k;   while (!vertex[h]) h++;
  l = k-1; while (!vertex[l]) l--;
  avma = av;
  return (long)floor((logcoef[h]-logcoef[l])/(double)(h-l) + 0.5);
}

/* change z into z*2^e, where z is real or complex of real */
static void
myshiftrc(GEN z, long e)
{
  if (typ(z)==t_COMPLEX)
  {
    if (signe(z[1])) setexpo(z[1], expo(z[1])+e);
    if (signe(z[2])) setexpo(z[2], expo(z[2])+e);
  }
  else
    if (signe(z)) setexpo(z,expo(z)+e);
}

/* return z*2^e, where z is integer or complex of integer (destroy z) */
static GEN
myshiftic(GEN z, long e)
{
  if (typ(z)==t_COMPLEX)
  {
    gel(z,1) = signe(z[1])? mpshift(gel(z,1),e): gen_0;
    gel(z,2) = mpshift(gel(z,2),e);
    return z;
  }
  return signe(z)? mpshift(z,e): gen_0;
}

/* as real_1 with precision in bits, not in words */
static GEN
myreal_1(long bit)
{
  if (bit < 0) bit = 0;
  return real_1(nbits2prec(bit));
}

static GEN
mygprecrc(GEN x, long prec, long e)
{
  GEN y;
  switch(typ(x))
  {
    case t_REAL: return signe(x)? rtor(x, prec): real_0_bit(e);
    case t_COMPLEX:
      y = cgetg(3,t_COMPLEX);
      gel(y,1) = mygprecrc(gel(x,1),prec,e);
      gel(y,2) = mygprecrc(gel(x,2),prec,e);
      return y;
    default: return gcopy(x);
  }
}

/* gprec behaves badly with the zero for polynomials.
The second parameter in mygprec is the precision in base 2 */
static GEN
mygprec(GEN x, long bit)
{
  long lx, i, e, prec;
  GEN y;

  if (bit < 0) bit = 0; /* should rarely happen */
  e = gexpo(x) - bit;
  prec = nbits2prec(bit);
  switch(typ(x))
  {
    case t_POL:
      lx = lg(x); y = cgetg(lx, t_POL); y[1] = x[1];
      for (i=2; i<lx; i++) gel(y,i) = mygprecrc(gel(x,i),prec,e);
      break;

    default: y = mygprecrc(x,prec,e);
  }
  return y;
}

/* normalize a polynomial p, that is change it with coefficients in Z[i],
after making product by 2^shift */
static GEN
pol_to_gaussint(GEN p, long shift)
{
  long i, l = lg(p);
  GEN q = cgetg(l, t_POL); q[1] = p[1];
  for (i=2; i<l; i++) gel(q,i) = gfloor2n(gel(p,i), shift);
  return q;
}

/* returns a polynomial q in Z[i][x] keeping bit bits of p */
static GEN
eval_rel_pol(GEN p, long bit)
{
  long i;
  for (i = 2; i < lg(p); i++)
    if (gcmp0(gel(p,i))) gel(p,i) = gen_0; /* bad behaviour of gexpo */
  return pol_to_gaussint(p, bit-gexpo(p)+1);
}

/* returns p(R*x)/R^n (in R or R[i]), R = exp(lrho), bit bits of precision */
static GEN
homothetie(GEN p, double lrho, long bit)
{
  GEN q, r, t, iR;
  long n = degpol(p), i;

  iR = mygprec(dblexp(-lrho),bit);
  q = mygprec(p, bit);
  r = cgetg(n+3,t_POL); r[1] = p[1];
  t = iR; r[n+2] = q[n+2];
  for (i=n-1; i>0; i--)
  {
    gel(r,i+2) = gmul(t, gel(q,i+2));
    t = mulrr(t, iR);
  }
  gel(r,2) = gmul(t, gel(q,2)); return r;
}

/* change q in 2^(n*e) p(x*2^(-e)), n=deg(q)  [ ~as above with R = 2^-e ]*/
static void
homothetie2n(GEN p, long e)
{
  if (e)
  {
    long i,n = lg(p)-1;
    for (i=2; i<=n; i++) myshiftrc(gel(p,i), (n-i)*e);
  }
}

/* return 2^f * 2^(n*e) p(x*2^(-e)), n=deg(q) */
static void
homothetie_gauss(GEN p, long e, long f)
{
  if (e || f)
  {
    long i, n = lg(p)-1;
    for (i=2; i<=n; i++) gel(p,i) = myshiftic(gel(p,i), f+(n-i)*e);
  }
}

/* Lower bound on the modulus of the largest root z_0
 * k is set to an upper bound for #{z roots, |z-z_0| < eps} */
static double
lower_bound(GEN p, long *k, double eps)
{
  long n = degpol(p), i, j;
  pari_sp ltop = avma;
  GEN a, s, S, ilc;
  double r, R, rho;

  if (n < 4) { *k = n; return 0.; }
  S = cgetg(5,t_VEC);
  a = cgetg(5,t_VEC); ilc = gdiv(real_1(DEFAULTPREC), gel(p,n+2));
  for (i=1; i<=4; i++) gel(a,i) = gmul(ilc,gel(p,n+2-i));
  /* i = 1 split out from next loop for efficiency and initialization */
  s = gel(a,1);
  gel(S,1) = gneg(s); /* Newton sum S_i */
  rho = r = gtodouble(gabs(s,3));
  R = r / n;
  for (i=2; i<=4; i++)
  {
    s = gmulsg(i,gel(a,i));
    for (j=1; j<i; j++) s = gadd(s, gmul(gel(S,j),gel(a,i-j)));
    gel(S,i) = gneg(s); /* Newton sum S_i */
    r = gtodouble(gabs(s,3));
    if (r > 0.)
    {
      r = exp(log(r/n) / (double)i);
      if (r > R) R = r;
    }
  }
  if (R > 0. && eps < 1.2)
    *k = (long)floor((rho/R + n) / (1 + exp(-eps)*cos(eps)));
  else
    *k = n;
  avma = ltop; return R;
}

/* log of modulus of the largest root of p with relative error tau */
static double
logmax_modulus(GEN p, double tau)
{
  GEN r,q,aux,gunr;
  pari_sp av, ltop = avma;
  long i,k,n=degpol(p),nn,bit,M,e;
  double rho,eps, tau2 = (tau > 3.0)? 0.5: tau/6.;

  r = cgeti(BIGDEFAULTPREC);
  av = avma;

  eps = - 1/log(1.5*tau2); /* > 0 */
  bit = (long) ((double) n*log2(1./tau2)+3*log2((double) n))+1;
  gunr = myreal_1(bit+2*n);
  aux = gdiv(gunr, gel(p,2+n));
  q = gmul(aux,p); gel(q,2+n) = gunr;
  e = findpower(q);
  homothetie2n(q,e);
  affsi(e, r);
  q = pol_to_gaussint(q, bit);
  M = (long) (log2( log(4.*n) / (2*tau2) )) + 2;
  nn = n;
  for (i=0,e=0;;)
  { /* nn = deg(q) */
    rho = lower_bound(q, &k, eps);
    if (rho > exp2(-(double)e)) e = (long)-floor(log2(rho));
    affii(shifti(addis(r,e), 1), r);
    if (++i == M) break;

    bit = (long) ((double)k * log2(1./tau2) +
                     (double)(nn-k)*log2(1./eps) + 3*log2((double)nn)) + 1;
    homothetie_gauss(q, e, bit-(long)floor(dbllog2(gel(q,2+nn))+0.5));
    nn -= polvaluation(q, &q);
    set_karasquare_limit(gexpo(q));
    q = gerepileupto(av, graeffe(q));
    tau2 *= 1.5; if (tau2 > 0.9) tau2 = 0.5;
    eps = -1/log(tau2); /* > 0 */
    e = findpower(q);
  }
  if (!signe(r)) { avma = ltop; return 0.; }
  r = itor(r, DEFAULTPREC); setexpo(r, expo(r) - M);
  avma = ltop; return -rtodbl(r) * LOG2; /* -log(2) sum e_i 2^-i */
}
GEN
logmax_modulus_bound(GEN P)
{
  return dblexp(logmax_modulus(P, 0.01) + 0.01);
}

/* log of modulus of the smallest root of p, with relative error tau */
static double
logmin_modulus(GEN p, double tau)
{
  pari_sp av = avma;
  double r;

  if (gcmp0(gel(p,2))) return -pariINFINITY;
  r = - logmax_modulus(polrecip_i(p),tau);
  avma = av; return r;
}

/* return the log of the k-th modulus (ascending order) of p, rel. error tau*/
static double
logmodulus(GEN p, long k, double tau)
{
  GEN q, gunr;
  long i, kk = k, imax, n = degpol(p), nn, bit, e;
  pari_sp av, ltop=avma;
  double r, tau2 = tau/6;

  bit = (long)(n * (2. + log2(3.*n) + log2(1./tau2)));
  gunr = myreal_1(bit);
  av = avma;
  q = gprec_w(p, nbits2prec(bit));
  q = gmul(gunr, q);
  e = newton_polygon(q,k);
  r = (double)e;
  homothetie2n(q,e);
  imax = (long)(log2(3./tau) + log2(log(4.*n)))+1;
  for (i=1; i<imax; i++)
  {
    q = eval_rel_pol(q,bit);
    kk -= polvaluation(q, &q);
    nn = degpol(q);

    set_karasquare_limit(bit);
    q = gerepileupto(av, graeffe(q));
    e = newton_polygon(q,kk);
    r += e / exp2((double)i);
    q = gmul(gunr, q);
    homothetie2n(q,e);

    tau2 *= 1.5; if (tau2 > 1.) tau2 = 1.;
    bit = 1 + (long)(nn*(2. + log2(3.*nn) + log2(1./tau2)));
  }
  avma = ltop; return -r * LOG2;
}

/* return the log of the k-th modulus r_k of p, rel. error tau, knowing that
 * rmin < r_k < rmax. This information helps because we may reduce precision
 * quicker */
static double
logpre_modulus(GEN p, long k, double tau, double lrmin, double lrmax)
{
  GEN q;
  long n = degpol(p), i, imax, imax2, bit;
  pari_sp ltop = avma, av;
  double lrho, aux, tau2 = tau/6.;

  aux = (lrmax - lrmin) / 2. + 4*tau2;
  imax = (long) log2(log((double)n)/ aux);
  if (imax <= 0) return logmodulus(p,k,tau);

  lrho  = (lrmin + lrmax) / 2;
  av = avma;
  bit = (long)(n*(2. + aux / LOG2 - log2(tau2)));
  q = homothetie(p, lrho, bit);
  imax2 = (long)(log2(3./tau) + log2(log(4.*n))) + 1;
  if (imax > imax2) imax = imax2;

  for (i=0; i<imax; i++)
  {
    q = eval_rel_pol(q,bit);
    set_karasquare_limit(bit);
    q = gerepileupto(av, graeffe(q));
    aux = 2*aux + 2*tau2;
    tau2 *= 1.5;
    bit = (long)(n*(2. + aux / LOG2 - log2(1-exp(-tau2))));
    q = gmul(myreal_1(bit),q);
  }
  aux = exp2((double)imax);
  aux = logmodulus(q,k, aux*tau/3.) / aux;
  avma = ltop; return lrho + aux;
}

static double
ind_maxlog2(GEN q)
{
  long i, k = -1;
  double L = - pariINFINITY;
  for (i=0; i<=degpol(q); i++)
  {
    double d = dbllog2(gel(q,2+i));
    if (d > L) { L = d; k = i; }
  }
  return k;
}

/* Returns k such that r_k e^(-tau) < R < r_{k+1} e^tau.
 * Assume that l <= k <= n-l */
static long
dual_modulus(GEN p, double lrho, double tau, long l)
{
  long i, imax, delta_k = 0, n = degpol(p), nn, v2, v, bit, ll = l;
  double tau2 = tau * 7./8.;
  pari_sp av = avma;
  GEN q;

  bit = 6*n - 5*l + (long)(n*(log2(1/tau2) + tau2 * 8./7.));
  q = homothetie(p, lrho, bit);
  imax = (long)(log(log(2.*n)/tau2)/log(7./4.)+1);

  for (i=0; i<imax; i++)
  {
    q = eval_rel_pol(q,bit); v2 = n - degpol(q);
    v = polvaluation(q, &q);
    ll -= max(v, v2); if (ll < 0) ll = 0;

    nn = degpol(q); delta_k += v;
    if (!nn) return delta_k;

    set_karasquare_limit(bit);
    q = gerepileupto(av, graeffe(q));
    tau2 *= 7./4.;
    bit = 6*nn - 5*ll + (long)(nn*(log2(1/tau2) + tau2 * 8./7.));
  }
  avma = av; return delta_k + (long)ind_maxlog2(q);
}

/********************************************************************/
/**                                                                **/
/**              FACTORS THROUGH CIRCLE INTEGRATION                **/
/**                                                                **/
/********************************************************************/
/* l power of 2 */
static void
fft(GEN Omega, GEN p, GEN f, long step, long l)
{
  pari_sp ltop;
  long i, l1, l2, l3, rapi, step4;
  GEN f1, f2, f3, f02, f13, g02, g13, ff;

  if (l == 2)
  {
    gel(f,0) = gadd(gel(p,0),gel(p,step));
    gel(f,1) = gadd(gel(p,0), gneg(gel(p,step))); return;
  }
  if (l == 4)
  {
    f1 = gadd(gel(p,0),   (GEN)p[step<<1]);
    f2 = gadd(gel(p,0),   gneg((GEN)p[step<<1]));
    f3 = gadd(gel(p,step),(GEN)p[3*step]);
    f02= gadd(gel(p,step),gneg((GEN)p[3*step]));
    f02 = mulcxI(f02);
    gel(f,0) = gadd(f1, f3);
    gel(f,1) = gadd(f2, f02);
    gel(f,2) = gadd(f1, gneg(f3));
    gel(f,3) = gadd(f2, gneg(f02)); return;
  }

  ltop = avma;
  l1 = l>>2; l2 = l1<<1; l3 = l1+l2; step4 = step<<2;
  fft(Omega,p,          f,   step4,l1);
  fft(Omega,p+step,     f+l1,step4,l1);
  fft(Omega,p+(step<<1),f+l2,step4,l1);
  fft(Omega,p+3*step,   f+l3,step4,l1);

  ff = cgetg(l+1,t_VEC);
  for (i=0; i<l1; i++)
  {
    rapi = step*i;
    f1 = gmul(gel(Omega,rapi),    gel(f,i+l1));
    f2 = gmul((GEN)Omega[rapi<<1], gel(f,i+l2));
    f3 = gmul((GEN)Omega[3*rapi],  gel(f,i+l3));

    f02 = gadd(gel(f,i),f2);
    g02 = gadd(gel(f,i),gneg(f2));
    f13 = gadd(f1,f3);
    g13 = mulcxI(gadd(f1,gneg(f3)));

    gel(ff,i+1)    = gadd(f02, f13);
    gel(ff,i+l1+1) = gadd(g02, g13);
    gel(ff,i+l2+1) = gadd(f02, gneg(f13));
    gel(ff,i+l3+1) = gadd(g02, gneg(g13));
  }
  ff = gerepilecopy(ltop,ff);
  for (i=0; i<l; i++) f[i] = ff[i+1];
}

/* e(1/N) */
static GEN
RUgen(long N, long bit)
{
  if (N == 2) return mpneg(real_1(nbits2prec(bit)));
  if (N == 4) return gi;
  return exp_Ir(divrs(Pi2n(1, nbits2prec(bit)), N));
}

/* N=2^k. returns a vector RU which contains exp(2*i*k*Pi/N), k=0..N-1 */
static GEN
initRU(long N, long bit)
{
  GEN *RU, z = RUgen(N, bit);
  long i, N2 = (N>>1), N4 = (N>>2), N8 = (N>>3);

  RU = (GEN*)cgetg(N+1,t_VEC); RU++;

  RU[0] = myreal_1(bit);
  RU[1] = z;
  for (i=1; i<N8; i++)
  {
    GEN t = RU[i];
    RU[i+1] = gmul(z, t);
    RU[N4-i] = mkcomplex(gel(t,2), gel(t,1));
  }
  for (i=0; i<N4; i++) RU[i+N4] = mulcxI(RU[i]);
  for (i=0; i<N2; i++) RU[i+N2] = gneg(RU[i]);
  return (GEN)RU;
}

/* as above, N arbitrary */
static GEN
initRUgen(long N, long bit)
{
  GEN *RU = (GEN*)cgetg(N+1,t_VEC), z = RUgen(N,bit);
  long i, k = (N+3)>>1;
  RU[0] = gen_1;
  RU[1] = z;
  for (i=2; i<k; i++) RU[i] = gmul(z, RU[i-1]);
  for (   ; i<N; i++) RU[i] = gconj(RU[N-i]);
  return (GEN)RU;
}

GEN
FFTinit(long k, long prec)
{
  if (k <= 0) pari_err(typeer,"FFTinit");
  return initRU(1 << k, bit_accuracy(prec)) - 1;
}

GEN
FFT(GEN x, GEN Omega)
{
  long i, l = lg(Omega), n = lg(x);
  GEN y, z;
  if (n > l || !is_vec_t(typ(x)) || typ(Omega) != t_VEC) pari_err(typeer,"FFT");
  if (n < l) {
    z = cgetg(l, t_VECSMALL); /* cf stackdummy */
    for (i = 1; i < n; i++) z[i] = x[i];
    for (     ; i < l; i++) gel(z,i) = gen_0;
  }
  else z = x;
  y = cgetg(l, t_VEC);
  fft(Omega+1, z+1, y+1, 1, l-1);
  return y;
}

/* returns 1 if p has only real coefficients, 0 else */
static int
isreal(GEN p)
{
  long n=degpol(p),i=0;

  while (i<=n && typ(p[i+2])!=t_COMPLEX) i++;
  return (i>n);
}

/* x non complex */
static GEN
abs_update_r(GEN x, double *mu) {
  GEN y = gabs(gprec_w(x, DEFAULTPREC), DEFAULTPREC);
  double ly = dblogr(y); if (ly < *mu) *mu = ly;
  return y;
}
/* return |x|, low accuracy. Set *mu = min(log(y), *mu) */
static GEN
abs_update(GEN x, double *mu) {
  GEN y, xr, yr;
  double ly;
  if (typ(x) != t_COMPLEX) return abs_update_r(x, mu);
  xr = gel(x,1);
  yr = gel(x,2);
  if (gcmp0(xr)) return abs_update_r(yr,mu);
  if (gcmp0(yr)) return abs_update_r(xr,mu);
  /* have to treat 0 specially: 0E-10 + 1e-20 = 0E-10 */
  xr = gprec_w(xr, DEFAULTPREC);
  yr = gprec_w(yr, DEFAULTPREC);
  y = gsqrt(gadd(gsqr(xr), gsqr(yr)), DEFAULTPREC);
  ly = dblogr(y); if (ly < *mu) *mu = ly;
  return y;
}

static void
parameters(GEN p, long *LMAX, double *mu, double *gamma,
           int polreal, double param, double param2)
{
  GEN q, pc, Omega, A, RU, prim, g, ONE,TWO;
  long n = degpol(p), bit, NN, K, i, j, Lmax;
  pari_sp av2, av = avma, lim = stack_lim(av, 1);

  bit = gexpo(p) + (long)param2+8;
  Lmax = 4; while (Lmax <= n) Lmax <<= 1;
  NN = (long)(param*3.14)+1; if (NN < Lmax) NN = Lmax;
  K = NN/Lmax; if (K & 1) K++;
  NN = Lmax*K;
  if (polreal) K = K/2+1;

  Omega = initRU(Lmax,bit);
  prim = RUgen(NN, bit);

  q = mygprec(p,bit) + 2;
  A = cgetg(Lmax+1,t_VEC); A++;
  pc= cgetg(Lmax+1,t_VEC); pc++;
  for (i=0; i <= n; i++) pc[i]= q[i];
  for (   ; i<Lmax; i++) gel(pc,i) = gen_0;

  *mu = pariINFINITY;
  g = real_0_bit(-bit);
  ONE = real_1(DEFAULTPREC);
  TWO = real2n(1, DEFAULTPREC);
  av2 = avma;
  RU = myreal_1(bit);
  for (i=0; i<K; i++)
  {
    if (i) {
      GEN z = RU;
      for (j=1; j<n; j++)
      {
        gel(pc,j) = gmul(gel(q,j),z);
        z = gmul(z,RU); /* RU = prim^i, z=prim^(ij) */
      }
      gel(pc,n) = gmul(gel(q,n),z);
    }

    fft(Omega,pc,A,1,Lmax);
    if (polreal && i>0 && i<K-1)
      for (j=0; j<Lmax; j++) g = addrr(g, divrr(TWO, abs_update(gel(A,j),mu)));
    else
      for (j=0; j<Lmax; j++) g = addrr(g, divrr(ONE, abs_update(gel(A,j),mu)));
    RU = gmul(RU, prim);
    if (low_stack(lim, stack_lim(av,1)))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"parameters");
      gerepileall(av2,2, &g,&RU);
    }
  }
  *gamma = dblogr(divrs(g,NN)) / LOG2;
  *LMAX = Lmax; avma = av;
}

/* NN is a multiple of Lmax */
static void
dft(GEN p, long k, long NN, long Lmax, long bit, GEN F, GEN H, long polreal)
{
  GEN Omega, q, qd, pc, pd, A, B, C, RU, aux, U, W, prim, prim2;
  long n = degpol(p), i, j, K;
  pari_sp ltop;

  Omega = initRU(Lmax,bit);
  prim = RUgen(NN, bit);
  RU = cgetg(n+2,t_VEC); RU++;

  K = NN/Lmax; if (polreal) K = K/2+1;
  q = mygprec(p,bit);
  qd = derivpol(q);

  A = cgetg(Lmax+1,t_VEC); A++;
  B = cgetg(Lmax+1,t_VEC); B++;
  C = cgetg(Lmax+1,t_VEC); C++;
  pc = cgetg(Lmax+1,t_VEC); pc++;
  pd = cgetg(Lmax+1,t_VEC); pd++;
  pc[0] = q[2];  for (i=n+1; i<Lmax; i++) gel(pc,i) = gen_0;
  pd[0] = qd[2]; for (i=n;   i<Lmax; i++) gel(pd,i) = gen_0;

  ltop = avma;
  W = cgetg(k+1,t_VEC);
  U = cgetg(k+1,t_VEC);
  for (i=1; i<=k; i++) gel(W,i) = gel(U,i) = gen_0;

  gel(RU,0) = gen_1;
  prim2 = myreal_1(bit);
  for (i=0; i<K; i++)
  {
    gel(RU,1) = prim2;
    for (j=1; j<n; j++) gel(RU,j+1) = gmul(gel(RU,j),prim2);
    /* RU[j] = prim^(ij)= prim2^j */

    for (j=1; j<n; j++) gel(pd,j) = gmul(gel(qd,j+2),gel(RU,j));
    fft(Omega,pd,A,1,Lmax);
    for (j=1; j<=n; j++) gel(pc,j) = gmul(gel(q,j+2),gel(RU,j));
    fft(Omega,pc,B,1,Lmax);
    for (j=0; j<Lmax; j++) gel(C,j) = ginv(gel(B,j));
    for (j=0; j<Lmax; j++) gel(B,j) = gmul(gel(A,j),gel(C,j));
    fft(Omega,B,A,1,Lmax);
    fft(Omega,C,B,1,Lmax);

    if (polreal) /* p has real coefficients */
    {
      if (i>0 && i<K-1)
      {
        for (j=1; j<=k; j++)
        {
          gel(W,j) = gadd(gel(W,j), gshift(real_i(gmul(gel(A,j+1),gel(RU,j+1))),1));
          gel(U,j) = gadd(gel(U,j), gshift(real_i(gmul(gel(B,j),gel(RU,j))),1));
        }
      }
      else
      {
        for (j=1; j<=k; j++)
        {
          gel(W,j) = gadd(gel(W,j), real_i(gmul(gel(A,j+1),gel(RU,j+1))));
          gel(U,j) = gadd(gel(U,j), real_i(gmul(gel(B,j),gel(RU,j))));
        }
      }
    }
    else
    {
      for (j=1; j<=k; j++)
      {
        gel(W,j) = gadd(gel(W,j), gmul(gel(A,j+1),gel(RU,j+1)));
        gel(U,j) = gadd(gel(U,j), gmul(gel(B,j),gel(RU,j)));
      }
    }
    prim2 = gmul(prim2,prim);
    gerepileall(ltop,3, &W,&U,&prim2);
  }

  for (i=1; i<=k; i++)
  {
    aux=gel(W,i);
    for (j=1; j<i; j++) aux = gadd(aux, gmul(gel(W,i-j),gel(F,k+2-j)));
    gel(F,k+2-i) = gdivgs(aux,-i*NN);
  }
  for (i=0; i<k; i++)
  {
    aux=gel(U,k-i);
    for (j=1+i; j<k; j++) aux = gadd(aux,gmul(gel(F,2+j),gel(U,j-i)));
    gel(H,i+2) = gdivgs(aux,NN);
  }
}

#define NEWTON_MAX 10
static GEN
refine_H(GEN F, GEN G, GEN HH, long bit, long Sbit)
{
  GEN H = HH, D, aux;
  pari_sp ltop = avma, lim = stack_lim(ltop, 1);
  long error, i, bit1, bit2;

  D = gsub(gen_1, grem(gmul(H,G),F)); error = gexpo(D);
  bit2 = bit + Sbit;
  for (i=0; error>-bit && i<NEWTON_MAX && error<=0; i++)
  {
    if (low_stack(lim, stack_lim(ltop,1)))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"refine_H");
      gerepileall(ltop,2, &D,&H);
    }
    bit1 = -error + Sbit;
    aux = gmul(mygprec(H,bit1), mygprec(D,bit1));
    aux = grem(mygprec(aux,bit1), mygprec(F,bit1));

    bit1 = -error*2 + Sbit; if (bit1 > bit2) bit1 = bit2;
    H = gadd(mygprec(H,bit1), aux);
    D = gsub(gen_1, grem(gmul(H,G),F));
    error = gexpo(D); if (error < -bit1) error = -bit1;
  }
  if (error > -bit/2) return NULL; /* FAIL */
  return gerepilecopy(ltop,H);
}

/* return 0 if fails, 1 else */
static long
refine_F(GEN p, GEN *F, GEN *G, GEN H, long bit, double gamma)
{
  GEN f0, FF, GG, r, HH = H;
  long error, i, bit1 = 0, bit2, Sbit, Sbit2,  enh, normF, normG, n = degpol(p);
  pari_sp av = avma, lim = stack_lim(av, 1);

  FF = *F; GG = poldivrem(p, FF, &r);
  error = gexpo(r); if (error <= -bit) error = 1-bit;
  normF = gexpo(FF);
  normG = gexpo(GG);
  enh = gexpo(H); if (enh < 0) enh = 0;
  Sbit = normF + 2*normG + enh + (long)(4.*log2((double)n)+gamma) + 1;
  Sbit2 = enh + 2*(normF+normG) + (long)(2.*gamma+5.*log2((double)n)) + 1;
  bit2 = bit + Sbit;
  for (i=0; error>-bit && i<NEWTON_MAX && error<=0; i++)
  {
    if (bit1 == bit2 && i >= 2) { Sbit += n; Sbit2 += n; bit2 += n; }
    if (low_stack(lim, stack_lim(av,1)))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"refine_F");
      gerepileall(av,4, &FF,&GG,&r,&HH);
    }

    bit1 = -error + Sbit2;
    HH = refine_H(mygprec(FF,bit1), mygprec(GG,bit1), mygprec(HH,bit1),
                  1-error, Sbit2);
    if (!HH) return 0; /* FAIL */

    bit1 = -error + Sbit;
    r = gmul(mygprec(HH,bit1), mygprec(r,bit1));
    f0 = grem(mygprec(r,bit1), mygprec(FF,bit1));

    bit1 = -2*error + Sbit; if (bit1 > bit2) bit1 = bit2;
    FF = gadd(mygprec(FF,bit1),f0);

    bit1 = -3*error + Sbit; if (bit1 > bit2) bit1 = bit2;
    GG = poldivrem(mygprec(p,bit1), mygprec(FF,bit1), &r);
    error = gexpo(r); if (error < -bit1) error = -bit1;
  }
  if (error>-bit) return 0; /* FAIL */
  *F = FF; *G = GG; return 1;
}

/* returns F and G from the unit circle U such that |p-FG|<2^(-bit) |cd|,
where cd is the leading coefficient of p */
static void
split_fromU(GEN p, long k, double delta, long bit,
            GEN *F, GEN *G, double param, double param2)
{
  GEN pp, FF, GG, H;
  long n = degpol(p), NN, bit2, Lmax;
  int polreal = isreal(p);
  pari_sp ltop;
  double mu, gamma;

  pp = gdiv(p, gel(p,2+n));
  parameters(pp, &Lmax,&mu,&gamma, polreal,param,param2);

  H  = cgetg(k+2,t_POL); H[1] = p[1];
  FF = cgetg(k+3,t_POL); FF[1]= p[1];
  gel(FF,k+2) = gen_1;

  NN = (long)(0.5/delta); NN |= 1; if (NN < 2) NN = 2;
  NN *= Lmax; ltop = avma;
  for(;;)
  {
    bit2 = (long)(((double)NN*delta-mu)/LOG2) + gexpo(pp) + 8;
    dft(pp, k, NN, Lmax, bit2, FF, H, polreal);
    if (refine_F(pp,&FF,&GG,H,bit,gamma)) break;
    NN <<= 1; avma = ltop;
  }
  *G = gmul(GG,gel(p,2+n)); *F = FF;
}

static void
optimize_split(GEN p, long k, double delta, long bit,
            GEN *F, GEN *G, double param, double param2)
{
  long n = degpol(p);
  GEN FF, GG;

  if (k <= n/2)
    split_fromU(p,k,delta,bit,F,G,param,param2);
  else
  {
    split_fromU(polrecip_i(p),n-k,delta,bit,&FF,&GG,param,param2);
    *F = polrecip(GG);
    *G = polrecip(FF);
  }
}

/********************************************************************/
/**                                                                **/
/**               SEARCH FOR SEPARATING CIRCLE                     **/
/**                                                                **/
/********************************************************************/

/* return p(2^e*x) *2^(-n*e) */
static void
scalepol2n(GEN p, long e)
{
  long i,n=lg(p)-1;
  for (i=2; i<=n; i++) gel(p,i) = gmul2n(gel(p,i),(i-n)*e);
}

/* returns p(x/R)*R^n */
static GEN
scalepol(GEN p, GEN R, long bit)
{
  GEN q,aux,gR;
  long i;

  aux = gR = mygprec(R,bit); q = mygprec(p,bit);
  for (i=lg(p)-2; i>=2; i--)
  {
    gel(q,i) = gmul(aux,gel(q,i));
    aux = gmul(aux,gR);
  }
  return q;
}

/* return (conj(a)X-1)^n * p[ (X-a) / (conj(a)X-1) ] */
static GEN
conformal_pol(GEN p, GEN a, long bit)
{
  GEN z, r, ma = gneg(a), ca = gconj(a);
  long n = degpol(p), i;
  pari_sp av = avma, lim = stack_lim(av,2);
  
  z = mkpoln(2, ca, negr(myreal_1(bit)));
  r = scalarpol(gel(p,2+n), 0);
  for (i=n-1; ; i--)
  {
    r = addmulXn(r, gmul(ma,r), 1); /* r *= (X - a) */
    r = gadd(r, gmul(z, gel(p,2+i)));
    if (i == 0) return gerepileupto(av, r);
    z = addmulXn(gmul(z,ca), gneg(z), 1); /* z *= conj(a)X - 1 */
    if (low_stack(lim, stack_lim(av,2)))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"conformal_pol");
      gerepileall(av,2, &r,&z);
    }
  }
}

static const double UNDEF = -100000.;

static double
logradius(double *radii, GEN p, long k, double aux, double *delta)
{
  long i, n = degpol(p);
  double lrho, lrmin, lrmax;
  if (k > 1)
  {
    i = k-1; while (i>0 && radii[i] == UNDEF) i--;
    lrmin = logpre_modulus(p,k,aux, radii[i], radii[k]);
  }
  else /* k=1 */
    lrmin = logmin_modulus(p,aux);
  radii[k] = lrmin;

  if (k+1<n)
  {
    i = k+2; while (i<=n && radii[i] == UNDEF) i++;
    lrmax = logpre_modulus(p,k+1,aux, radii[k+1], radii[i]);
  }
  else /* k+1=n */
    lrmax = logmax_modulus(p,aux);
  radii[k+1] = lrmax;

  lrho = radii[k];
  for (i=k-1; i>=1; i--)
  {
    if (radii[i] == UNDEF || radii[i] > lrho)
      radii[i] = lrho;
    else
      lrho = radii[i];
  }
  lrho = radii[k+1];
  for (i=k+1; i<=n; i++)
  {
    if (radii[i] == UNDEF || radii[i] < lrho)
      radii[i] = lrho;
    else
      lrho = radii[i];
  }
  *delta = (lrmax - lrmin) / 2;
  if (*delta > 1.) *delta = 1.;
  return (lrmin + lrmax) / 2;
}

static void
update_radius(long n, double *radii, double lrho, double *par, double *par2)
{
  double t, param = 0., param2 = 0.;
  long i;
  for (i=1; i<=n; i++)
  {
    radii[i] -= lrho;
    t = fabs(rtodbl( ginv(subsr(1, dblexp(radii[i]))) ));
    param += t; if (t > 1.) param2 += log2(t);
  }
  *par = param; *par2 = param2;
}

/* apply the conformal mapping then split from U */
static void
conformal_mapping(double *radii, GEN ctr, GEN p, long k, long bit,
                  double aux, GEN *F,GEN *G)
{
  long bit2, n = degpol(p), i;
  pari_sp ltop = avma, av;
  GEN q, FF, GG, a, R;
  double lrho, delta, param, param2;
  /* n * (2.*log2(2.732)+log2(1.5)) + 1 */
  bit2 = bit + (long)(n*3.4848775) + 1;
  a = sqrtr_abs( stor(3, 2*MEDDEFAULTPREC - 2) );
  a = divrs(a, -6);
  a = gmul(mygprec(a,bit2), mygprec(ctr,bit2)); /* a = -ctr/2sqrt(3) */

  av = avma;
  q = conformal_pol(mygprec(p,bit2), a, bit2);
  for (i=1; i<=n; i++)
    if (radii[i] != UNDEF) /* update array radii */
    {
      pari_sp av2 = avma;
      GEN t, r = dblexp(radii[i]), r2 = gsqr(r);
      /* 2(r^2 - 1) / (r^2 - 3(r-1)) */
      t = divrr(shiftr((subrs(r2,1)),1), subrr(r2, mulsr(3,subrs(r,1))));
      radii[i] = dblogr(addsr(1,t)) / 2;
      avma = av2;
    }
  lrho = logradius(radii, q,k,aux/10., &delta);
  update_radius(n, radii, lrho, &param, &param2);

  bit2 += (long)(n * fabs(lrho)/LOG2 + 1.);
  R = mygprec(dblexp(-lrho), bit2);
  q = scalepol(q,R,bit2);
  gerepileall(av,2, &q,&R);

  optimize_split(q,k,delta,bit2,&FF,&GG,param,param2);
  bit2 += n; R = ginv(R);
  FF = scalepol(FF,R,bit2);
  GG = scalepol(GG,R,bit2);

  a = mygprec(a,bit2);
  FF = conformal_pol(FF,a,bit2);
  GG = conformal_pol(GG,a,bit2);

  a = ginv(gsub(gen_1, gnorm(a)));
  FF = gmul(FF, gpowgs(a,k));
  GG = gmul(GG, gpowgs(a,n-k));

  *F = mygprec(FF,bit+n);
  *G = mygprec(GG,bit+n); gerepileall(ltop,2, F,G);
}

/* split p, this time without scaling. returns in F and G two polynomials
 * such that |p-FG|< 2^(-bit)|p| */
static void
split_2(GEN p, long bit, GEN ctr, double thickness, GEN *F, GEN *G)
{
  GEN q, FF, GG, R;
  double aux, delta, param, param2;
  long n = degpol(p), i, j, k, bit2;
  double lrmin, lrmax, lrho, *radii;

  init_dalloc();
  radii = (double*) stackmalloc((n+1) * sizeof(double));

  for (i=2; i<n; i++) radii[i] = UNDEF;
  aux = thickness/(double)(4 * n);
  lrmin = logmin_modulus(p, aux);
  lrmax = logmax_modulus(p, aux);
  radii[1] = lrmin;
  radii[n] = lrmax;
  i = 1; j = n;
  lrho = (lrmin + lrmax) / 2;
  k = dual_modulus(p, lrho, aux, 1);
  if (5*k < n || (n < 2*k && 5*k < 4*n))
    { lrmax = lrho; j=k+1; radii[j] = lrho; }
  else
    { lrmin = lrho; i=k;   radii[i] = lrho; }
  while (j > i+1)
  {
    if (i+j == n+1)
      lrho = (lrmin + lrmax) / 2;
    else
    {
      double kappa = 2. - log(1. + min(i,n-j)) / log(1. + min(j,n-i));
      if (i+j < n+1) lrho = lrmax * kappa + lrmin;
      else           lrho = lrmin * kappa + lrmax;
      lrho /= 1+kappa;
    }
    aux = (lrmax - lrmin) / (4*(j-i));
    k = dual_modulus(p, lrho, aux, min(i,n+1-j));
    if (k-i < j-k-1 || (k-i == j-k-1 && 2*k > n))
      { lrmax = lrho; j=k+1; radii[j] = lrho - aux; }
    else
      { lrmin = lrho; i=k;   radii[i] = lrho + aux; }
  }
  aux = lrmax - lrmin;

  if (ctr)
  {
    lrho = (lrmax + lrmin) / 2;
    for (i=1; i<=n; i++)
      if (radii[i] != UNDEF) radii[i] -= lrho;

    bit2 = bit + (long)(n * fabs(lrho)/LOG2 + 1.);
    R = mygprec(dblexp(-lrho), bit2);
    q = scalepol(p,R,bit2);
    conformal_mapping(radii, ctr, q, k, bit2, aux, &FF, &GG);
  }
  else
  {
    lrho = logradius(radii, p, k, aux/10., &delta);
    update_radius(n, radii, lrho, &param, &param2);

    bit2 = bit + (long)(n * fabs(lrho)/LOG2 + 1.);
    R = mygprec(dblexp(-lrho), bit2);
    q = scalepol(p,R,bit2);
    optimize_split(q, k, delta, bit2, &FF, &GG, param, param2);
  }
  bit  += n;
  bit2 += n; R = ginv(mygprec(R,bit2));
  *F = mygprec(scalepol(FF,R,bit2), bit);
  *G = mygprec(scalepol(GG,R,bit2), bit);
}

/* procedure corresponding to steps 5,6,.. page 44 in RR n. 1852 */
/* put in F and G two polynomial such that |p-FG|<2^(-bit)|p|
 * where the maximum modulus of the roots of p is <=1.
 * Assume sum of roots is 0. */
static void
split_1(GEN p, long bit, GEN *F, GEN *G)
{
  long i, imax, n = degpol(p), polreal = isreal(p), ep = gexpo(p), bit2 = bit+n;
  GEN TWO, ctr, q, qq, FF, GG, v, gr, r, newq;
  double lrmin, lrmax, lthick;
  const double LOG3 = 1.098613;

  lrmax = logmax_modulus(p, 0.01);
  gr = mygprec(dblexp(-lrmax), bit2);
  q = scalepol(p,gr,bit2);

  bit2 = bit + gexpo(q) - ep + (long)((double)n*2.*log2(3.)+1);
  TWO = myreal_1(bit2); setexpo(TWO,1);
  v = cgetg(5,t_VEC);
  gel(v,1) = TWO;
  gel(v,2) = mpneg(TWO);
  gel(v,3) = pureimag(gel(v,1));
  gel(v,4) = pureimag(gel(v,2));
  q = mygprec(q,bit2); lthick = 0;
  newq = ctr = NULL; /* -Wall */
  imax = polreal? 3: 4;
  for (i=1; i<=imax; i++)
  {
    qq = translate_pol(q, gel(v,i));
    lrmin = logmin_modulus(qq,0.05);
    if (LOG3 > lrmin + lthick)
    {
      double lquo = logmax_modulus(qq,0.05) - lrmin;
      if (lquo > lthick) { lthick = lquo; newq = qq; ctr = gel(v,i); }
    }
    if (lthick > LOG2) break;
    if (polreal && i==2 && lthick > LOG3 - LOG2) break;
  }
  bit2 = bit + gexpo(newq) - ep + (long)(n*LOG3/LOG2 + 1);
  split_2(newq, bit2, ctr, lthick, &FF, &GG);
  r = gneg(mygprec(ctr,bit2));
  FF = translate_pol(FF,r);
  GG = translate_pol(GG,r);

  gr = ginv(gr); bit2 = bit - ep + gexpo(FF)+gexpo(GG);
  *F = scalepol(FF,gr,bit2);
  *G = scalepol(GG,gr,bit2);
}

/* put in F and G two polynomials such that |P-FG|<2^(-bit)|P|,
where the maximum modulus of the roots of p is < 0.5 */
static int
split_0_2(GEN p, long bit, GEN *F, GEN *G)
{
  GEN q, b, FF, GG;
  long n = degpol(p), k, bit2, eq;
  double aux = dbllog2(gel(p,n+1)) - dbllog2(gel(p,n+2));

  /* beware double overflow */
  if (aux >= 0 && (aux > 1e4 || exp2(aux) > 2.5*n)) return 0;

  aux = (aux < -300)? 0.: n*log2(1 + exp2(aux)/(double)n);
  bit2 = bit+1 + (long)(log2((double)n) + aux);

  q = mygprec(p,bit2);
  b = gdivgs(gdiv(gel(q,n+1),gel(q,n+2)),-n);
  q = translate_pol(q,b); gel(q,n+1) = gen_0; eq = gexpo(q);
  k = 0;
  while (k <= n/2 && (- gexpo(gel(q,k+2)) > bit2 + 2*(n-k) + eq
                      || gcmp0(gel(q,k+2)))) k++;
  if (k > 0)
  {
    if (k > n/2) k = n/2;
    bit2 += k<<1;
    FF = monomial(myreal_1(bit2), k, 0);
    GG = RgX_shift_shallow(q, -k);
  }
  else
  {
    split_1(q,bit2,&FF,&GG);
    bit2 = bit + gexpo(FF) + gexpo(GG) - gexpo(p) + (long)aux+1;
    FF = mygprec(FF,bit2);
  }
  GG = mygprec(GG,bit2); b = mygprec(gneg(b),bit2);
  *F = translate_pol(FF, b);
  *G = translate_pol(GG, b); return 1;
}

/* put in F and G two polynomials such that |P-FG|<2^(-bit)|P|.
 * Assume max_modulus(p) < 2 */
static void
split_0_1(GEN p, long bit, GEN *F, GEN *G)
{
  GEN FF, GG;
  long n, bit2, normp;

  if  (split_0_2(p,bit,F,G)) return;

  normp = gexpo(p);
  scalepol2n(p,2); /* p := 4^(-n) p(4*x) */
  n = degpol(p); bit2 = bit + 2*n + gexpo(p) - normp;
  split_1(mygprec(p,bit2), bit2,&FF,&GG);
  scalepol2n(FF,-2);
  scalepol2n(GG,-2); bit2 = bit + gexpo(FF) + gexpo(GG) - normp;
  *F = mygprec(FF,bit2);
  *G = mygprec(GG,bit2);
}

/* put in F and G two polynomials such that |P-FG|<2^(-bit)|P| */
static void
split_0(GEN p, long bit, GEN *F, GEN *G)
{
  const double LOG1_9 = 0.6418539;
  long n = degpol(p), k = 0;
  GEN q;

  while (gexpo(gel(p,k+2)) < -bit && k <= n/2) k++;
  if (k > 0)
  {
    if (k > n/2) k = n/2;
    *F = monomial(myreal_1(bit), k, 0);
    *G = RgX_shift_shallow(p, -k);
  }
  else
  {
    double lr = logmax_modulus(p, 0.05);
    if (lr < LOG1_9) split_0_1(p, bit, F, G);
    else
    {
      q = polrecip_i(p);
      lr = logmax_modulus(q,0.05);
      if (lr < LOG1_9)
      {
        split_0_1(q, bit, F, G);
        *F = polrecip(*F);
        *G = polrecip(*G);
      }
      else
        split_2(p,bit,NULL, 1.2837,F,G);
    }
  }
}

/********************************************************************/
/**                                                                **/
/**                ERROR ESTIMATE FOR THE ROOTS                    **/
/**                                                                **/
/********************************************************************/

static GEN
root_error(long n, long k, GEN roots_pol, long pari_err, GEN shatzle)
{
  GEN rho, d, eps, epsbis, eps2, prod, aux, rap = NULL;
  long i, j, m;

  d = cgetg(n+1,t_VEC);
  for (i=1; i<=n; i++)
  {
    if (i!=k)
    {
      aux = gsub(gel(roots_pol,i), gel(roots_pol,k));
      gel(d,i) = gabs(mygprec(aux,31), DEFAULTPREC);
    }
  }
  rho = gabs(mygprec(gel(roots_pol,k),31), DEFAULTPREC);
  if (expo(rho) < 0) rho = real_1(DEFAULTPREC);
  eps = mulrr(rho, shatzle);
  aux = shiftr(gpowgs(rho,n), pari_err);

  for (j=1; j<=2 || (j<=5 && gcmp(rap, dbltor(1.2)) > 0); j++)
  {
    m = n; prod = real_1(DEFAULTPREC);
    epsbis = mulrr(eps, dbltor(1.25));
    for (i=1; i<=n; i++)
    {
      if (i != k && cmprr(gel(d,i),epsbis) > 0)
      {
        m--;
        prod = mulrr(prod, subrr(gel(d,i),eps));
      }
    }
    eps2 = sqrtnr(mpdiv(shiftr(aux,2*m-2), prod), m);
    rap = divrr(eps,eps2); eps = eps2;
  }
  return eps;
}

/* round a complex or real number x to an absolute value of 2^(-bit) */
static GEN
mygprec_absolute(GEN x, long bit)
{
  long e;
  GEN y;

  switch(typ(x))
  {
    case t_REAL:
      e = expo(x) + bit;
      return (e <= 0 || !signe(x))? real_0_bit(-bit): rtor(x, nbits2prec(e));
    case t_COMPLEX:
      if (gexpo(gel(x,2)) < -bit) return mygprec_absolute(gel(x,1),bit);
      y = cgetg(3,t_COMPLEX);
      gel(y,1) = mygprec_absolute(gel(x,1),bit);
      gel(y,2) = mygprec_absolute(gel(x,2),bit);
      return y;
    default: return x;
  }
}

static long
a_posteriori_errors(GEN p, GEN roots_pol, long pari_err)
{
  long i, e, n = degpol(p), e_max = -(long)EXPOBITS;
  GEN sigma, shatzle, x;

  pari_err += (long)log2((double)n) + 1;
  if (pari_err > -2) return 0;
  sigma = real2n(-pari_err, 3);
  /*  2 / ((s - 1)^(1/n) - 1) */
  shatzle = divsr(2, subrs(sqrtnr(subrs(sigma,1),n), 1));
  for (i=1; i<=n; i++)
  {
    x = root_error(n,i,roots_pol,pari_err,shatzle);
    e = gexpo(x); if (e > e_max) e_max = e;
    gel(roots_pol,i) = mygprec_absolute(gel(roots_pol,i), -e);
  }
  return e_max;
}

/********************************************************************/
/**                                                                **/
/**                           MAIN                                 **/
/**                                                                **/
/********************************************************************/
static GEN
append_clone(GEN r, GEN a) { a = gclone(a); appendL(r, a); return a; }

/* put roots in placeholder roots_pol so that |P - L_1...L_n| < 2^(-bit)|P|
 * returns prod (x-roots_pol[i]) */
static GEN
split_complete(GEN p, long bit, GEN roots_pol)
{
  long n = degpol(p);
  pari_sp ltop;
  GEN p1, F, G, a, b, m1, m2;

  if (n == 1)
  {
    a = gneg_i(gdiv(gel(p,2), gel(p,3)));
    (void)append_clone(roots_pol,a); return p;
  }
  ltop = avma;
  if (n == 2)
  {
    F = gsub(gsqr(gel(p,3)), gmul2n(gmul(gel(p,2),gel(p,4)), 2));
    F = gsqrt(F, nbits2prec(bit));
    p1 = ginv(gmul2n(gel(p,4),1));
    a = gneg_i(gmul(gadd(F,gel(p,3)), p1));
    b =        gmul(gsub(F,gel(p,3)), p1);
    a = append_clone(roots_pol,a);
    b = append_clone(roots_pol,b); avma = ltop;
    a = mygprec(a, 3*bit);
    b = mygprec(b, 3*bit);
    return gmul(gel(p,4), mkpoln(3, gen_1, gneg(gadd(a,b)), gmul(a,b)));
  }
  split_0(p,bit,&F,&G);
  m1 = split_complete(F,bit,roots_pol);
  m2 = split_complete(G,bit,roots_pol);
  return gerepileupto(ltop, gmul(m1,m2));
}

/* bound the log of the largest root of p */
double
cauchy_bound(GEN p)
{
  pari_sp av = avma;
  long i, n = degpol(p), prec = DEFAULTPREC;
  GEN lc, y, q = gmul(p, real_1(prec));
  double L = 0, Lmax = -pariINFINITY;

  if (n <= 0) pari_err(constpoler,"cauchy_bound");

  lc = gabs(gel(q,n+2),prec); /* leading coefficient */
  lc = ginv(lc);
  for (i=0; i<n; i++)
  {
    y = gel(q,i+2); if (gcmp0(y)) continue;
    L = dblogr(gmul(gabs(y,prec), lc)) / (n-i);
    if (L > Lmax) Lmax = L;
  }
  avma = av; return Lmax + LOG2;
}

static GEN
mygprecrc_special(GEN x, long prec, long e)
{
  GEN y;
  switch(typ(x))
  {
    case t_REAL:
      if (!signe(x)) return real_0_bit(min(e, expo(x)));
      return (prec > lg(x))? rtor(x, prec): x;
    case t_COMPLEX:
      y = cgetg(3,t_COMPLEX);
      gel(y,1) = mygprecrc_special(gel(x,1),prec,e);
      gel(y,2) = mygprecrc_special(gel(x,2),prec,e);
      return y;
    default: return x;
  }
}

/* like mygprec but keep at least the same precision as before */
static GEN
mygprec_special(GEN x, long bit)
{
  long lx, i, e, prec;
  GEN y;

  if (bit < 0) bit = 0; /* should not happen */
  e = gexpo(x) - bit;
  prec = nbits2prec(bit);
  switch(typ(x))
  {
    case t_POL:
      lx = lg(x); y = cgetg(lx,t_POL); y[1] = x[1];
      for (i=2; i<lx; i++) gel(y,i) = mygprecrc_special(gel(x,i),prec,e);
      break;

    default: y = mygprecrc_special(x,prec,e);
  }
  return y;
}

static GEN
fix_roots1(GEN r, GEN *m, long bit)
{
  long i, l = lg(r);
  GEN allr = cgetg(l, t_VEC);
  for (i=1; i<l; i++)
  {
    GEN t = gel(r,i);
    gel(allr,i) = gcopy(t); gunclone(t);
  }
  return allr;
}
static GEN
fix_roots(GEN r, GEN *m, long h, long bit)
{
  long i, j, k, l, prec;
  GEN allr, ro1;
  if (h == 1) return fix_roots1(r, m, bit);
  ro1 = initRUgen(h, bit);
  prec = nbits2prec(bit);
  l = lg(r)-1;
  allr = cgetg(h*l+1, t_VEC);
  for (k=1,i=1; i<=l; i++)
  {
    GEN p2, p1 = gel(r,i);
    p2 = (h == 2)? gsqrt(p1, prec): gsqrtn(p1, utoipos(h), NULL, prec);
    for (j=0; j<h; j++) gel(allr,k++) = gmul(p2, gel(ro1,j));
    gunclone(p1);
  }
  *m = roots_to_pol(allr, 0);
  return allr;
}

static GEN
all_roots(GEN p, long bit)
{
  GEN lc, pd, q, roots_pol, m;
  long bit0,  bit2, n = degpol(p), i, e, h;
  pari_sp av;

  pd = poldeflate(p, &h); lc = leading_term(pd);
  e = (long)((2/LOG2) * cauchy_bound(pd)); if (e < 0) e = 0;
  bit0 = bit + gexpo(pd) - gexpo(lc) + (long)log2(n/h)+1+e;
  bit2 = bit0; e = 0;
  for (av=avma,i=1;; i++,avma=av)
  {
    roots_pol = cget1(n+1,t_VEC);
    bit2 += e + (n << i);
    q = gmul(myreal_1(bit2), mygprec(pd,bit2));
    q[1] = evalsigne(1)|evalvarn(0);
    m = split_complete(q,bit2,roots_pol);
    roots_pol = fix_roots(roots_pol, &m, h, bit2);
    q = mygprec_special(p,bit2); lc = leading_term(q);
    q[1] = evalsigne(1)|evalvarn(0);
    if (h > 1) m = gmul(m,lc);

    e = gexpo(gsub(q, m)) - gexpo(lc) + (long)log2((double)n) + 1;
    if (e < -2*bit2) e = -2*bit2; /* avoid e = -pariINFINITY */
    if (e < 0)
    {
      e = bit + a_posteriori_errors(p,roots_pol,e);
      if (e < 0) return roots_pol;
    }
    if (DEBUGLEVEL > 7)
      fprintferr("all_roots: restarting, i = %ld, e = %ld\n", i,e);
  }
}

INLINE int
isexactscalar(GEN x) { long tx = typ(x); return is_rational_t(tx); }

static int
isexactpol(GEN p)
{
  long i,n = degpol(p);
  for (i=0; i<=n; i++)
    if (!isexactscalar(gel(p,i+2))) return 0;
  return 1;
}

static long
isvalidcoeff(GEN x)
{
  switch (typ(x))
  {
    case t_INT: case t_REAL: case t_FRAC: return 1;
    case t_COMPLEX:
      if (isvalidcoeff(gel(x,1)) && isvalidcoeff(gel(x,2))) return 1;
  }
  return 0;
}

static long
isvalidpol(GEN p)
{
  long i,n = lg(p);
  for (i=2; i<n; i++)
    if (!isvalidcoeff(gel(p,i))) return 0;
  return 1;
}

static GEN
solve_exact_pol(GEN p, long bit)
{
  long i, j, k, m, n = degpol(p), iroots = 0;
  GEN ex, factors, v = zerovec(n);

  factors = ZX_squff(Q_primpart(p), &ex);
  for (i=1; i<lg(factors); i++)
  {
    GEN roots_fact = all_roots(gel(factors,i), bit);
    n = degpol(factors[i]);
    m = ex[i];
    for (j=1; j<=n; j++)
      for (k=1; k<=m; k++) v[++iroots] = roots_fact[j];
  }
  return v;
}

/* return the roots of p with absolute error bit */
static GEN
roots_com(GEN q, long bit)
{
  GEN L, p;
  long v = polvaluation_inexact(q, &p);
  if (lg(p) == 3) L = cgetg(1,t_VEC); /* constant polynomial */
  else L = isexactpol(p)? solve_exact_pol(p,bit): all_roots(p,bit);
  if (v)
  {
    GEN M, z, t = gel(q,2);
    long i, x, y, l, n;

    if (isexactzero(t)) x = -bit;
    else
    {
      n = gexpo(t);
      x = n / v; l = degpol(q);
      for (i = v; i <= l; i++)
      {
        t  = gel(q,i+2);
        if (isexactzero(t)) continue;
        y = (n - gexpo(t)) / i;
        if (y < x) x = y;
      }
    }
    z = real_0_bit(x); l = v + lg(L);
    M = cgetg(l, t_VEC); L -= v;
    for (i = 1; i <= v; i++) gel(M,i) = z;
    for (     ; i <  l; i++) M[i] = L[i];
    L = M;
  }
  return L;
}

static GEN
tocomplex(GEN x, long l)
{
  GEN y = cgetg(3,t_COMPLEX);

  gel(y,1) = cgetr(l);
  if (typ(x) == t_COMPLEX)
    { gel(y,2) = cgetr(l); gaffect(x,y); }
  else
    { gaffect(x,gel(y,1)); gel(y,2) = real_0(l); }
 return y;
}

/* Check if x is approximately real with precision e */
int
isrealappr(GEN x, long e)
{
  long tx=typ(x),lx,i;
  switch(tx)
  {
    case t_INT: case t_REAL: case t_FRAC:
      return 1;
    case t_COMPLEX:
      return (gexpo(gel(x,2)) < e);
    case t_POL: case t_SER: case t_RFRAC: case t_VEC: case t_COL: case t_MAT:
      lx = lg(x);
      for (i=lontyp[tx]; i<lx; i++)
        if (! isrealappr(gel(x,i),e)) return 0;
      return 1;
    default: pari_err(typeer,"isrealappr"); return 0;
  }
}

/* x,y are t_COMPLEX */
static int
isconj(GEN x, GEN y, long e)
{
  pari_sp av = avma;
  long i= (gexpo( gsub(gel(x,1),gel(y,1)) ) < e
        && gexpo( gadd(gel(x,2),gel(y,2)) ) < e);
  avma = av; return i;
}

/* the vector of roots of p, with absolute error 2^(- bit_accuracy(l)) */
GEN
roots(GEN p, long l)
{
  pari_sp av = avma;
  long n, i, k, s, t, e;
  GEN c, L, p1, res, rea, com;

  if (gcmp0(p)) pari_err(zeropoler,"roots");
  if (typ(p) != t_POL)
  {
    if (!isvalidcoeff(p)) pari_err(typeer,"roots");
    return cgetg(1,t_VEC); /* constant polynomial */
  }
  if (!isvalidpol(p)) pari_err(talker,"invalid coefficients in roots");
  if (lg(p) == 3) return cgetg(1,t_VEC); /* constant polynomial */

  if (l < 3) l = 3;
  L = roots_com(p, bit_accuracy(l)); n = lg(L);
  if (!isreal(p))
  {
    res = cgetg(n,t_COL);
    for (i=1; i<n; i++) gel(res,i) = tocomplex(gel(L,i),l);
    return gerepileupto(av,res);
  }
  e = 5 - bit_accuracy(l);
  rea = cgetg(n,t_COL); s = 0;
  com = cgetg(n,t_COL); t = 0;
  for (i=1; i<n; i++)
  {
    p1 = gel(L,i);
    if (isrealappr(p1,e)) {
      if (typ(p1) == t_COMPLEX) p1 = gel(p1,1);
      gel(rea,++s) = p1;
    }
    else gel(com,++t) = p1;
  }
  setlg(rea,s+1); rea = sort(rea);
  res = cgetg(n,t_COL);
  for (i=1; i<=s; i++) gel(res,i) = tocomplex(gel(rea,i),l);
  for (i=1; i<=t; i++)
  {
    c = gel(com,i); if (!c) continue;
    gel(res,++s) = tocomplex(c,l);
    for (k=i+1; k<=t; k++)
    {
      p1 = gel(com,k); if (!p1) continue;
      if (isconj(c,p1,e))
      {
        gel(res,++s) = tocomplex(p1,l);
        com[k] = 0; break;
      }
    }
    if (k==n) pari_err(bugparier,"roots (conjugates)");
  }
  return gerepileupto(av,res);
}

GEN
roots0(GEN p, long flag,long l)
{
  switch(flag)
  {
    case 0: return roots(p,l);
    case 1: return rootsold(p,l);
    default: pari_err(flagerr,"polroots");
  }
  return NULL; /* not reached */
}

/* clean up roots. If root is real replace it by its real part */
GEN
cleanroots(GEN p, long prec)
{
  GEN s, r = roots(p,prec);
  long i, l = lg(r);
  for (i=1; i<l; i++)
  {
    s = gel(r,i);
    if (signe(s[2])) break; /* remaining roots are complex */
    r[i] = s[1]; /* root is real; take real part */
  }
  return r;
}
