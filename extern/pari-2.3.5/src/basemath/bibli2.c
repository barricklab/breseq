/* $Id: bibli2.c 10293 2008-06-10 15:57:06Z kb $

Copyright (C) 2000  The PARI group.

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
/**                                                               **/
/**                      SPECIAL POLYNOMIALS                      **/
/**                                                               **/
/*******************************************************************/
#ifdef LONG_IS_64BIT
# define SQRTVERYBIGINT 3037000500   /* ceil(sqrt(VERYBIGINT)) */
#else
# define SQRTVERYBIGINT 46341
#endif

/* Tchebichev polynomial: T0=1; T1=X; T(n)=2*X*T(n-1)-T(n-2)
 * T(n) = (n/2) sum_{k=0}^{n/2} a_k x^(n-2k)
 *   where a_k = (-1)^k 2^(n-2k) (n-k-1)! / k!(n-2k)! is an integer
 *   and a_0 = 2^(n-1), a_k / a_{k-1} = - (n-2k+2)(n-2k+1) / 4k(n-k) */
GEN
tchebi(long n, long v) /* Assume 4*n < VERYBIGINT */
{
  long k, l;
  pari_sp av;
  GEN q,a,r;

  if (v<0) v = 0;
  if (n < 0) n = -n;
  if (n==0) return pol_1[v];
  if (n==1) return pol_x[v];

  q = cgetg(n+3, t_POL); r = q + n+2;
  a = int2n(n-1);
  gel(r--,0) = a;
  gel(r--,0) = gen_0;
  if (n < SQRTVERYBIGINT)
    for (k=1,l=n; l>1; k++,l-=2)
    {
      av = avma;
      a = divis(mulis(a, l*(l-1)), 4*k*(n-k));
      a = gerepileuptoint(av, negi(a));
      gel(r--,0) = a;
      gel(r--,0) = gen_0;
    }
  else
    for (k=1,l=n; l>1; k++,l-=2)
    {
      av = avma;
      a = mulis(mulis(a, l), l-1);
      a = divis(divis(a, 4*k), n-k);
      a = gerepileuptoint(av, negi(a));
      gel(r--,0) = a;
      gel(r--,0) = gen_0;
    }
  q[1] = evalsigne(1) | evalvarn(v);
  return q;
}

GEN addmulXn(GEN x, GEN y, long d);
/* Legendre polynomial */
/* L0=1; L1=X; (n+1)*L(n+1)=(2*n+1)*X*L(n)-n*L(n-1) */
GEN
legendre(long n, long v)
{
  long m;
  pari_sp av, tetpil, lim;
  GEN p0,p1,p2;

  if (v<0) v = 0;
  if (n < 0) pari_err(talker,"negative degree in legendre");
  if (n==0) return pol_1[v];
  if (n==1) return pol_x[v];

  p0=pol_1[v]; av=avma; lim=stack_lim(av,2);
  p1=gmul2n(pol_x[v],1);
  for (m=1; m<n; m++)
  {
    p2 = addmulXn(gmulsg(4*m+2,p1), gmulsg(-4*m,p0), 1);
    setvarn(p2,v);
    p0 = p1; tetpil=avma; p1 = gdivgs(p2,m+1);
    if (low_stack(lim, stack_lim(av,2)))
    {
      GEN *gptr[2];
      if(DEBUGMEM>1) pari_warn(warnmem,"legendre");
      p0=gcopy(p0); gptr[0]=&p0; gptr[1]=&p1;
      gerepilemanysp(av,tetpil,gptr,2);
    }
  }
  tetpil=avma; return gerepile(av,tetpil,gmul2n(p1,-n));
}

/* cyclotomic polynomial */
GEN
cyclo(long n, long v)
{
  long d, q, m;
  pari_sp av=avma, tetpil;
  GEN yn,yd;

  if (n <= 0) pari_err(talker, "argument must be positive in polcyclo");
  if (v<0) v = 0;
  yn = yd = pol_1[0];
  for (d=1; d*d<=n; d++)
  {
    if (n%d) continue;
    q = n/d;
    m = mu(utoipos(q));
    if (m)
    { /* y *= (x^d - 1) */
      if (m>0) yn = addmulXn(yn, gneg(yn), d);
      else     yd = addmulXn(yd, gneg(yd), d);
    }
    if (q==d) break;
    m = mu(utoipos(d));
    if (m)
    { /* y *= (x^q - 1) */
      if (m>0) yn = addmulXn(yn, gneg(yn), q);
      else     yd = addmulXn(yd, gneg(yd), q);
    }
  }
  tetpil=avma; yn = gerepile(av,tetpil,RgX_div(yn,yd));
  setvarn(yn,v); return yn;
}

/* compute prod (L*x +/- a[i]) */
GEN
roots_to_pol_intern(GEN L, GEN a, long v, int plus)
{
  long i,k,lx = lg(a), code;
  GEN p1,p2;
  if (lx == 1) return pol_1[v];
  p1 = cgetg(lx, t_VEC);
  code = evalsigne(1)|evalvarn(v);
  for (k=1,i=1; i<lx-1; i+=2)
  {
    p2 = cgetg(5,t_POL); gel(p1,k++) = p2;
    gel(p2,2) = gmul(gel(a,i),gel(a,i+1));
    gel(p2,3) = gadd(gel(a,i),gel(a,i+1));
    if (plus == 0) gel(p2,3) = gneg(gel(p2,3));
    gel(p2,4) = L; p2[1] = code;
  }
  if (i < lx)
  {
    p2 = cgetg(4,t_POL); gel(p1,k++) = p2;
    p2[1] = code = evalsigne(1)|evalvarn(v);
    gel(p2,2) = plus? gel(a,i): gneg(gel(a,i));
    gel(p2,3) = L;
  }
  setlg(p1, k); return divide_conquer_prod(p1, gmul);
}

GEN
roots_to_pol(GEN a, long v)
{
  return roots_to_pol_intern(gen_1,a,v,0);
}

/* prod_{i=1..r1} (x - a[i]) prod_{i=1..r2} (x - a[i])(x - conj(a[i]))*/
GEN
roots_to_pol_r1r2(GEN a, long r1, long v)
{
  long i,k,lx = lg(a), code;
  GEN p1;
  if (lx == 1) return pol_1[v];
  p1 = cgetg(lx, t_VEC);
  code = evalsigne(1)|evalvarn(v);
  for (k=1,i=1; i<r1; i+=2)
  {
    GEN p2 = cgetg(5,t_POL); gel(p1,k++) = p2;
    gel(p2,2) = gmul(gel(a,i),gel(a,i+1));
    gel(p2,3) = gneg(gadd(gel(a,i),gel(a,i+1)));
    gel(p2,4) = gen_1; p2[1] = code;
  }
  if (i < r1+1)
    gel(p1,k++) = gadd(pol_x[v], gneg(gel(a,i)));
  for (i=r1+1; i<lx; i++)
  {
    GEN p2 = cgetg(5,t_POL); gel(p1,k++) = p2;
    gel(p2,2) = gnorm(gel(a,i));
    gel(p2,3) = gneg(gtrace(gel(a,i)));
    gel(p2,4) = gen_1; p2[1] = code;
  }
  setlg(p1, k); return divide_conquer_prod(p1, gmul);
}

/********************************************************************/
/**                                                                **/
/**                  HILBERT & PASCAL MATRICES                     **/
/**                                                                **/
/********************************************************************/
GEN
mathilbert(long n) /* Hilbert matrix of order n */
{
  long i,j;
  GEN p;

  if (n < 0) n = 0;
  p = cgetg(n+1,t_MAT);
  for (j=1; j<=n; j++)
  {
    gel(p,j) = cgetg(n+1,t_COL);
    for (i=1+(j==1); i<=n; i++)
      gcoeff(p,i,j) = mkfrac(gen_1, utoipos(i+j-1));
  }
  if (n) gcoeff(p,1,1) = gen_1;
  return p;
}

/* q-Pascal triangle = (choose(i,j)_q) (ordinary binomial if q = NULL) */
GEN
matqpascal(long n, GEN q)
{
  long i, j, I;
  pari_sp av = avma;
  GEN m, *qpow = NULL; /* gcc -Wall */

  if (n<0) n = -1;
  n++; m = cgetg(n+1,t_MAT);
  for (j=1; j<=n; j++) gel(m,j) = cgetg(n+1,t_COL);
  if (q)
  {
    I = (n+1)/2;
    if (I > 1) { qpow = (GEN*)new_chunk(I+1); qpow[2]=q; }
    for (j=3; j<=I; j++) qpow[j] = gmul(q, qpow[j-1]);
  }
  for (i=1; i<=n; i++)
  {
    I = (i+1)/2; gcoeff(m,i,1)= gen_1;
    if (q)
    {
      for (j=2; j<=I; j++)
        gcoeff(m,i,j) = gadd(gmul(qpow[j],gcoeff(m,i-1,j)), gcoeff(m,i-1,j-1));
    }
    else
    {
      for (j=2; j<=I; j++)
        gcoeff(m,i,j) = addii(gcoeff(m,i-1,j), gcoeff(m,i-1,j-1));
    }
    for (   ; j<=i; j++) gcoeff(m,i,j) = gcoeff(m,i,i+1-j);
    for (   ; j<=n; j++) gcoeff(m,i,j) = gen_0;
  }
  return gerepilecopy(av, m);
}

/********************************************************************/
/**                                                                **/
/**                  LAPLACE TRANSFORM (OF A SERIES)               **/
/**                                                                **/
/********************************************************************/

GEN
laplace(GEN x)
{
  pari_sp av = avma;
  long i, l = lg(x), e = valp(x);
  GEN y, t;

  if (typ(x) != t_SER) pari_err(talker,"not a series in laplace");
  if (e < 0) pari_err(talker,"negative valuation in laplace");
  y = cgetg(l,t_SER);
  t = mpfact(e); y[1] = x[1];
  for (i=2; i<l; i++)
  {
    gel(y,i) = gmul(t, gel(x,i));
    e++; t = mulsi(e,t);
  }
  return gerepilecopy(av,y);
}

/********************************************************************/
/**                                                                **/
/**              CONVOLUTION PRODUCT (OF TWO SERIES)               **/
/**                                                                **/
/********************************************************************/

GEN
convol(GEN x, GEN y)
{
  long j, lx, ly, ex, ey, vx = varn(x);
  GEN z;

  if (typ(x) != t_SER || typ(y) != t_SER) pari_err(talker,"not a series in convol");
  if (varn(y) != vx) pari_err(talker,"different variables in convol");
  ex = valp(x); lx = lg(x) + ex; x -= ex;
  ey = valp(y); ly = lg(y) + ey; y -= ey;
  /* inputs shifted: x[i] and y[i] now correspond to monomials of same degree */
  if (ly < lx) lx = ly; /* min length */
  if (ex < ey) ex = ey; /* max valuation */
  if (lx - ex < 3) return zeroser(vx, lx-2);

  z = cgetg(lx - ex, t_SER);
  z[1] = evalvalp(ex) | evalvarn(vx);
  for (j = ex+2; j<lx; j++) gel(z,j-ex) = gmul(gel(x,j),gel(y,j));
  return normalize(z);
}

/******************************************************************/
/**                                                              **/
/**                       PRECISION CHANGES                      **/
/**                                                              **/
/******************************************************************/

GEN
gprec(GEN x, long l)
{
  long tx = typ(x), lx, i;
  GEN y;

  if (l <= 0) pari_err(talker,"precision<=0 in gprec");
  switch(tx)
  {
    case t_REAL:
      return rtor(x, ndec2prec(l));

    case t_PADIC:
      if (!signe(x[4])) return zeropadic(gel(x,2), l+precp(x));
      y=cgetg(5,t_PADIC);
      y[1]=x[1]; setprecp(y,l);
      gel(y,2) = gcopy(gel(x,2));
      gel(y,3) = gpowgs(gel(x,2),l);
      gel(y,4) = modii(gel(x,4), gel(y,3));
      break;

    case t_SER:
      if (lg(x) == 2) return zeroser(varn(x), l);
      y=cgetg(l+2,t_SER); y[1]=x[1]; l++; i=l;
      lx = lg(x);
      if (l>=lx)
	for ( ; i>=lx; i--) gel(y,i) = gen_0;
      for ( ; i>=2; i--) gel(y,i) = gcopy(gel(x,i));
      break;

    case t_COMPLEX: case t_POLMOD: case t_POL: case t_RFRAC:
    case t_VEC: case t_COL: case t_MAT:
      y = init_gen_op(x, tx, &lx, &i);
      for (; i<lx; i++) gel(y,i) = gprec(gel(x,i),l);
      break;
    default: y = gcopy(x);
  }
  return y;
}

/* internal: precision given in word length (including codewords) */
GEN
gprec_w(GEN x, long pr)
{
  long tx = typ(x), lx, i;
  GEN y;

  switch(tx)
  {
    case t_REAL:
      if (signe(x)) return rtor(x,pr);
      i = -bit_accuracy(pr);
      if (i < expo(x)) return real_0_bit(i);
      y = cgetr(2); y[1] = x[1]; return y;
    case t_COMPLEX: case t_POLMOD: case t_POL: case t_RFRAC:
    case t_VEC: case t_COL: case t_MAT:
      y = init_gen_op(x, tx, &lx, &i);
      for (; i<lx; i++) gel(y,i) = gprec_w(gel(x,i),pr);
      break;
    default: return x;
  }
  return y;
}

/* internal: precision given in word length (including codewords), truncate
 * mantissa to precision 'pr' but never _increase_ it */
GEN
gprec_wtrunc(GEN x, long pr)
{
  long tx = typ(x), lx, i;
  GEN y;

  switch(tx)
  {
    case t_REAL:
      return (signe(x) && lg(x) > pr)? rtor(x,pr): x;
    case t_COMPLEX: case t_POLMOD: case t_POL: case t_RFRAC:
    case t_VEC: case t_COL: case t_MAT:
      y = init_gen_op(x, tx, &lx, &i);
      for (; i<lx; i++) gel(y,i) = gprec_wtrunc(gel(x,i),pr);
      break;
    default: return x;
  }
  return y;
}

/*******************************************************************/
/**                                                               **/
/**                     RECIPROCAL POLYNOMIAL                     **/
/**                                                               **/
/*******************************************************************/
/* return coefficients s.t x = x_0 X^n + ... + x_n */
GEN
polrecip(GEN x)
{
  long lx = lg(x), i, j;
  GEN y = cgetg(lx,t_POL);

  if (typ(x) != t_POL) pari_err(typeer,"polrecip");
  y[1] = x[1]; for (i=2,j=lx-1; i<lx; i++,j--) gel(y,i) = gcopy(gel(x,j));
  return normalizepol_i(y,lx);
}

/* as above. Internal (don't copy or normalize) */
GEN
polrecip_i(GEN x)
{
  long lx = lg(x), i, j;
  GEN y = cgetg(lx,t_POL);
  y[1] = x[1]; for (i=2,j=lx-1; i<lx; i++,j--) y[i] = x[j];
  return y;
}

/*******************************************************************/
/**                                                               **/
/**                      BINOMIAL COEFFICIENTS                    **/
/**                                                               **/
/*******************************************************************/

GEN
binomial(GEN n, long k)
{
  long i;
  pari_sp av;
  GEN y;

  if (k <= 1)
  {
    if (is_noncalc_t(typ(n))) pari_err(typeer,"binomial");
    if (k < 0) return gen_0;
    if (k == 0) return gen_1;
    return gcopy(n);
  }
  av = avma;
  if (typ(n) == t_INT)
  {
    if (signe(n) > 0)
    {
      GEN z = subis(n,k);
      if (cmpis(z,k) < 0) 
      {
        k = itos(z); avma = av;
        if (k <= 1)
        {
          if (k < 0) return gen_0;
          if (k == 0) return gen_1;
          return icopy(n);
        }
      }
    }
    /* k > 1 */
    if (lgefint(n) == 3 && signe(n) > 0)
    {
      ulong N = itou(n);
      y = seq_umul(N-(ulong)k+1, N);
    }
    else
    {
      y = cgetg(k+1,t_VEC);
      for (i=1; i<=k; i++) gel(y,i) = subis(n,i-1);
      y = divide_conquer_prod(y,mulii);
    }
    y = diviiexact(y, mpfact(k));
  }
  else
  {
    y = cgetg(k+1,t_VEC);
    for (i=1; i<=k; i++) gel(y,i) = gsubgs(n,i-1);
    y = divide_conquer_prod(y,gmul);
    y = gdiv(y, mpfact(k));
  }
  return gerepileupto(av, y);
}

/* Assume n >= 1, return bin, bin[k+1] = binomial(n, k) */
GEN
vecbinome(long n)
{
  long d = (n + 1)/2, k;
  GEN bin = cgetg(n+2, t_VEC), *C;
  C = (GEN*)(bin + 1); /* C[k] = binomial(n, k) */
  C[0] = gen_1;
  for (k=1; k <= d; k++)
  {
    pari_sp av = avma;
    C[k] = gerepileuptoint(av, diviiexact(mulsi(n-k+1, C[k-1]), utoipos(k)));
  }
  for (   ; k <= n; k++) C[k] = C[n - k];
  return bin;
}

/********************************************************************/
/**                                                                **/
/**                  POLYNOMIAL INTERPOLATION                      **/
/**                                                                **/
/********************************************************************/
/* assume n > 1 */
GEN
polint_i(GEN xa, GEN ya, GEN x, long n, GEN *ptdy)
{
  long i, m, ns = 0, tx = typ(x);
  pari_sp av = avma, tetpil;
  GEN y, c, d, dy;

  if (!xa)
  {
    xa = cgetg(n+1, t_VEC);
    for (i=1; i<=n; i++) gel(xa,i) = utoipos(i);
    xa++;
  }
  if (is_scalar_t(tx) && tx != t_INTMOD && tx != t_PADIC && tx != t_POLMOD)
  {
    GEN dif = NULL, dift;
    for (i=0; i<n; i++)
    {
      dift = gabs(gsub(x,gel(xa,i)), MEDDEFAULTPREC);
      if (!dif || gcmp(dift,dif)<0) { ns = i; dif = dift; }
    }
  }
  c = new_chunk(n);
  d = new_chunk(n); for (i=0; i<n; i++) c[i] = d[i] = ya[i];
  y = gel(d,ns--);
  dy = NULL; tetpil = 0; /* gcc -Wall */
  for (m=1; m<n; m++)
  {
    for (i=0; i<n-m; i++)
    {
      GEN ho = gsub(gel(xa,i),x);
      GEN hp = gsub(gel(xa,i+m),x), den = gsub(ho,hp);
      if (gcmp0(den)) pari_err(talker,"two abcissas are equal in polint");
      den = gdiv(gsub(gel(c,i+1),gel(d,i)), den);
      gel(c,i) = gmul(ho,den);
      gel(d,i) = gmul(hp,den);
    }
    dy = (2*(ns+1) < n-m)? gel(c,ns+1): gel(d,ns--);
    tetpil = avma; y = gadd(y,dy);
  }
  if (!ptdy) y = gerepile(av,tetpil,y);
  else
  {
    GEN *gptr[2];
    *ptdy=gcopy(dy); gptr[0]=&y; gptr[1]=ptdy;
    gerepilemanysp(av,tetpil,gptr,2);
  }
  return y;
}

GEN
polint(GEN xa, GEN ya, GEN x, GEN *ptdy)
{
  long tx=typ(xa), ty, lx=lg(xa);

  if (ya) ty = typ(ya); else { ya = xa; ty = tx; xa = NULL; }

  if (! is_vec_t(tx) || ! is_vec_t(ty))
    pari_err(talker,"not vectors in polinterpolate");
  if (lx != lg(ya))
    pari_err(talker,"different lengths in polinterpolate");
  if (lx <= 2)
  {
    if (lx == 1) pari_err(talker,"no data in polinterpolate");
    ya=gcopy(gel(ya,1)); if (ptdy) *ptdy = ya;
    return ya;
  }
  if (!x) x = pol_x[0];
  return polint_i(xa? xa+1: xa,ya+1,x,lx-1,ptdy);
}

/***********************************************************************/
/*                                                                     */
/*                          SET OPERATIONS                             */
/*                                                                     */
/***********************************************************************/
GEN
gtoset(GEN x)
{
  pari_sp av;
  long i,c,tx,lx;
  GEN y;

  if (!x) return cgetg(1, t_VEC);
  tx = typ(x); lx = lg(x);
  if (!is_vec_t(tx))
  {
    if (tx != t_LIST)
      { y=cgetg(2,t_VEC); gel(y,1) = GENtocanonicalstr(x); return y; }
    lx = lgeflist(x)-1; x++;
  }
  if (lx==1) return cgetg(1,t_VEC);
  av=avma; y=cgetg(lx,t_VEC);
  for (i=1; i<lx; i++) gel(y,i) = GENtocanonicalstr(gel(x,i));
  y = sort(y);
  c=1;
  for (i=2; i<lx; i++)
    if (!gequal(gel(y,i), gel(y,c))) y[++c] = y[i];
  setlg(y,c+1); return gerepilecopy(av,y);
}

long
setisset(GEN x)
{
  long lx,i;

  if (typ(x)!=t_VEC) return 0;
  lx=lg(x)-1; if (!lx) return 1;
  for (i=1; i<lx; i++)
    if (typ(x[i]) != t_STR || gcmp(gel(x,i+1),gel(x,i))<=0) return 0;
  return typ(x[i]) == t_STR;
}

/* looks if y belongs to the set x and returns the index if yes, 0 if no */
long
gen_search_aux(GEN x, GEN y, long flag, void *data, int (*cmp)(void*,GEN,GEN))
{
  long lx,j,li,ri,fl, tx = typ(x);

  if (tx==t_VEC) lx = lg(x);
  else
  {
    if (tx!=t_LIST) pari_err(talker,"not a set in setsearch");
    lx=lgeflist(x)-1; x++;
  }
  if (lx==1) return flag? 1: 0;

  li=1; ri=lx-1;
  do
  {
    j = (ri+li)>>1; fl = cmp(data,gel(x,j),y);
    if (!fl) return flag? 0: j;
    if (fl<0) li=j+1; else ri=j-1;
  } while (ri>=li);
  if (!flag) return 0;
  return (fl<0)? j+1: j;
}

static int cmp_nodata(void *data, GEN x, GEN y)
{
  int (*cmp)(GEN,GEN)=(int (*)(GEN,GEN)) data;
  return cmp(x,y);
}

long
gen_search(GEN x, GEN y, long flag, int (*cmp)(GEN,GEN))
{
  return gen_search_aux(x, y, flag, (void *)cmp, cmp_nodata);
}

long
setsearch(GEN x, GEN y, long flag)
{
  pari_sp av = avma;
  long res;
  if (typ(y) != t_STR) y = GENtocanonicalstr(y);
  res=gen_search(x,y,flag,gcmp);
  avma=av;
  return res;
}
long 
ZV_search(GEN x, GEN y) { return gen_search(x, y, 0, cmpii); }

GEN
ZV_sort_uniq(GEN L)
{
  long i, c, l = lg(L);
  pari_sp av = avma;
  GEN perm;

  if (l < 2) return cgetg(1, typ(L));
  perm = gen_sort(L, cmp_C, &cmpii);
  L = vecpermute(L, perm);
  c = 1;
  for (i = 2; i < l; i++)
    if (!equalii(gel(L,i), gel(L,c))) L[++c] = L[i];
  setlg(L, c+1); return gerepilecopy(av, L);
}

#if 0
GEN
gen_union(GEN x, GEN y, int (*cmp)(GEN,GEN))
{
  if (typ(x) != t_VEC || typ(y) != t_VEC) pari_err(talker,"not a set in setunion");

}
#endif

GEN
setunion(GEN x, GEN y)
{
  pari_sp av=avma, tetpil;
  GEN z;

  if (typ(x) != t_VEC || typ(y) != t_VEC) pari_err(talker,"not a set in setunion");
  z=shallowconcat(x,y); tetpil=avma; return gerepile(av,tetpil,gtoset(z));
}

GEN
setintersect(GEN x, GEN y)
{
  long i, lx, c;
  pari_sp av=avma;
  GEN z;

  if (!setisset(x) || !setisset(y)) pari_err(talker,"not a set in setintersect");
  lx=lg(x); z=cgetg(lx,t_VEC); c=1;
  for (i=1; i<lx; i++)
    if (setsearch(y, gel(x,i), 0)) z[c++] = x[i];
  setlg(z,c); return gerepilecopy(av,z);
}

GEN
gen_setminus(GEN set1, GEN set2, int (*cmp)(GEN,GEN))
{
  pari_sp ltop=avma;
  long find;
  long i,j,k;
  GEN  diff=cgetg(lg(set1),t_VEC);
  for(i=1,j=1,k=1; i < lg(set1); i++)
  {
    for(find=0; j < lg(set2); j++)
    {
      int s=cmp(gel(set1,i),gel(set2,j));
      if (s<0)  break ;
      if (s>0)  continue;
      find=1;
    }
    if (!find)
      diff[k++]=set1[i];
  }
  setlg(diff,k);
  return gerepilecopy(ltop,diff);
}

GEN
setminus(GEN x, GEN y)
{
  if (!setisset(x) || !setisset(y)) pari_err(talker,"not a set in setminus");
  return gen_setminus(x,y,gcmp);
}

/***********************************************************************/
/*                                                                     */
/*               OPERATIONS ON DIRICHLET SERIES                        */
/*                                                                     */
/***********************************************************************/

/* Addition, subtraction and scalar multiplication of Dirichlet series
   are done on the corresponding vectors */

static long
dirval(GEN x)
{
  long i = 1, lx = lg(x);
  while (i < lx && gcmp0(gel(x,i))) i++;
  return i;
}

GEN
dirmul(GEN x, GEN y)
{
  pari_sp av = avma, lim = stack_lim(av, 1);
  long lx, ly, lz, dx, dy, i, j, k;
  GEN z;

  if (typ(x)!=t_VEC || typ(y)!=t_VEC) pari_err(typeer,"dirmul");
  dx = dirval(x); lx = lg(x);
  dy = dirval(y); ly = lg(y);
  if (ly-dy < lx-dx) { swap(x,y); lswap(lx,ly); lswap(dx,dy); }
  lz = min(lx*dy,ly*dx);
  z = zerovec(lz-1);
  for (j=dx; j<lx; j++)
  {
    GEN c = gel(x,j);
    if (gcmp0(c)) continue;
    if (gcmp1(c))
      for (k=dy,i=j*dy; i<lz; i+=j,k++) gel(z,i) = gadd(gel(z,i),gel(y,k));
    else
    {
      if (gcmp_1(c))
        for (k=dy,i=j*dy; i<lz; i+=j,k++) gel(z,i) = gsub(gel(z,i),gel(y,k));
      else
        for (k=dy,i=j*dy; i<lz; i+=j,k++) gel(z,i) = gadd(gel(z,i),gmul(c,gel(y,k)));
    }
    if (low_stack(lim, stack_lim(av,1)))
    {
      if (DEBUGLEVEL) fprintferr("doubling stack in dirmul\n");
      z = gerepilecopy(av,z);
    }
  }
  return gerepilecopy(av,z);
}

GEN
dirdiv(GEN x, GEN y)
{
  pari_sp av = avma;
  long lx,ly,lz,dx,dy,i,j;
  GEN z,p1;

  if (typ(x)!=t_VEC || typ(y)!=t_VEC) pari_err(typeer,"dirmul");
  dx = dirval(x); lx = lg(x);
  dy = dirval(y); ly = lg(y);
  if (dy != 1 || ly == 1) pari_err(talker,"not an invertible dirseries in dirdiv");
  lz = min(lx,ly*dx); p1 = gel(y,1);
  if (!gcmp1(p1)) { y = gdiv(y,p1); x = gdiv(x,p1); } else x = shallowcopy(x);
  z = zerovec(lz-1);
  for (j=dx; j<lz; j++)
  {
    p1=gel(x,j); gel(z,j) = p1;
    if (gcmp0(p1)) continue;
    if (gcmp1(p1))
      for (i=j+j; i<lz; i+=j) gel(x,i) = gsub(gel(x,i),gel(y,i/j));
    else
    {
      if (gcmp_1(p1))
        for (i=j+j; i<lz; i+=j) gel(x,i) = gadd(gel(x,i),gel(y,i/j));
      else
        for (i=j+j; i<lz; i+=j) gel(x,i) = gsub(gel(x,i),gmul(p1,gel(y,i/j)));
    }
  }
  return gerepilecopy(av,z);
}

/*************************************************************************/
/**									**/
/**			         RANDOM					**/
/**									**/
/*************************************************************************/

GEN
genrand(GEN N)
{
  if (!N) return stoi( pari_rand31() );
  if (typ(N)!=t_INT || signe(N)<=0) pari_err(talker,"invalid bound in random");
  return randomi(N);
}

long
getstack(void) { return top-avma; }

long
gettime(void) { return timer(); }

/***********************************************************************/
/**							              **/
/**       		     PERMUTATIONS                             **/
/**								      **/
/***********************************************************************/

GEN
numtoperm(long n, GEN x)
{
  pari_sp av;
  long i, r;
  GEN v;

  if (n < 0) pari_err(talker,"n too small (%ld) in numtoperm",n);
  if (typ(x) != t_INT) pari_err(arither1);
  v = cgetg(n+1, t_VEC);
  v[1] = 1; av = avma;
  if (signe(x) <= 0) x = modii(x, mpfact(n));
  for (r=2; r<=n; r++)
  {
    long a;
    x = divis_rem(x, r,&a);
    for (i=r; i>=a+2; i--) v[i] = v[i-1];
    v[i] = r;
    if ((r & 0x1f) == 0) x = gerepileuptoint(av, x);
  }
  avma = av;
  for (i=1; i<=n; i++) gel(v,i) = stoi(v[i]);
  return v;
}

GEN
permtonum(GEN x)
{
  long lx=lg(x)-1, n=lx, last, ind, tx = typ(x);
  pari_sp av=avma;
  GEN ary,res;

  if (!is_vec_t(tx)) pari_err(talker,"not a vector in permtonum");
  ary = cgetg(lx+1,t_VECSMALL);
  for (ind=1; ind<=lx; ind++)
  {
    res = gel(++x, 0);
    if (typ(res) != t_INT) pari_err(typeer,"permtonum");
    ary[ind] = itos(res);
  }
  ary++; res = gen_0;
  for (last=lx; last>0; last--)
  {
    lx--; ind = lx;
    while (ind>0 && ary[ind] != last) ind--;
    res = addis(mulis(res,last), ind);
    while (ind++<lx) ary[ind-1] = ary[ind];
  }
  if (!signe(res)) res = mpfact(n);
  return gerepileuptoint(av, res);
}

/********************************************************************/
/**                                                                **/
/**                       MODREVERSE                               **/
/**                                                                **/
/********************************************************************/
/* return y such that Mod(y, charpoly(Mod(a,T)) = Mod(a,T) */
GEN
modreverse_i(GEN a, GEN T)
{
  pari_sp av = avma;
  long n = degpol(T);
  GEN y;

  if (n <= 0) return gcopy(a);
  if (n == 1)
    return gerepileupto(av, gneg(gdiv(gel(T,2), gel(T,3))));
  if (gcmp0(a) || typ(a) != t_POL) pari_err(talker,"reverse polmod does not exist");

  y = RgXV_to_RgM(RgX_powers(a,T,n-1), n);
  y = gauss(y, col_ei(n, 2));
  return gerepilecopy(av, RgV_to_RgX(y, varn(T)));
}

GEN
polymodrecip(GEN x)
{
  long v, n;
  GEN T, a, y;

  if (typ(x)!=t_POLMOD) pari_err(talker,"not a polmod in modreverse");
  T = gel(x,1);
  a = gel(x,2);
  n = degpol(T); if (n <= 0) return gcopy(x);
  v = varn(T);
  y = cgetg(3,t_POLMOD);
  gel(y,1) = (n==1)? gsub(pol_x[v], a): caract2(T, a, v);
  gel(y,2) = modreverse_i(a, T); return y;
}

/********************************************************************/
/**                                                                **/
/**                           HEAPSORT                             **/
/**                                                                **/
/********************************************************************/

#define icmp(a,b) ((a)>(b)?1:(a)<(b)?-1:0)

int
pari_compare_lg(GEN x, GEN y)
{
  return icmp(lg(x),lg(y));
}

int
pari_compare_long(long *a,long *b)
{
  return icmp(*a,*b);
}

static int
pari_compare_small(GEN x, GEN y)
{
  return icmp((long)x,(long)y);
}

#undef icmp

struct veccmp_s
{
  long lk;
  GEN k;
  int (*cmp)(GEN,GEN);
};

static int
veccmp(void *data, GEN x, GEN y)
{
  struct veccmp_s *v=(struct veccmp_s *) data;
  long i,s;

  for (i=1; i<v->lk; i++)
  {
    s = v->cmp(gel(x,v->k[i]), gel(y,v->k[i]));
    if (s) return s;
  }
  return 0;
}

static GEN
gen_sortspec(GEN v, long n, void *data, int (*cmp)(void*,GEN,GEN))
{
  long nx=n>>1, ny=n-nx;
  long m, ix, iy;
  GEN x, y;
  GEN w=cgetg(n+1,t_VECSMALL);
  if (n<=2)
  {
    if (n==1) 
      w[1]=1;
    else if (n==2)
    {
      if (cmp(data,gel(v,1),gel(v,2))<=0) { w[1]=1; w[2]=2; }
      else { w[1]=2; w[2]=1; }
    }
    return w;
  }
  x=gen_sortspec(v,nx,data,cmp);
  y=gen_sortspec(v+nx,ny,data,cmp);
  for (m=1, ix=1, iy=1; ix<=nx && iy<=ny; )
    if (cmp(data, gel(v,x[ix]), gel(v,y[iy]+nx))<=0)
      w[m++]=x[ix++];
    else
      w[m++]=y[iy++]+nx;
  for(;ix<=nx;) w[m++]=x[ix++];
  for(;iy<=ny;) w[m++]=y[iy++]+nx;
  avma = (pari_sp) w;
  return w;
}

/* Sort x = vector of elts, using cmp to compare them.
 *  flag & cmp_IND: indirect sort: return permutation that would sort x
 *  flag & cmp_C  : as cmp_IND, but return permutation as vector of C-longs
 */
GEN
gen_sort_aux(GEN x, long flag, void *data, int (*cmp)(void*,GEN,GEN))
{
  long i, j;
  long tx = typ(x), lx = lg(x);
  GEN y;

  if (tx == t_LIST) { lx = lgeflist(x)-1; tx = t_VEC; x++; }
  if (!is_matvec_t(tx) && tx != t_VECSMALL) pari_err(typeer,"gen_sort");
  if      (flag & cmp_C)   tx = t_VECSMALL;
  else if (flag & cmp_IND) tx = t_VEC;
  if (lx<=2)
  {
    y=cgetg(lx,tx);
    if (lx==1) return y;
    if (lx==2)
    {
      if      (flag & cmp_C)   y[1] = 1;
      else if (flag & cmp_IND) gel(y,1) = gen_1;
      else if (tx == t_VECSMALL) y[1] = x[1];
      else gel(y,1) = gcopy(gel(x,1)); 
      return y;
    }
  }

  y = gen_sortspec(x,lx-1,data,cmp);

  if (flag & cmp_REV)
  { /* reverse order */
    for (j=1; j<=((lx-1)>>1); j++) 
    {
      long z=y[j];
      y[j]=y[lx-j];
      y[lx-j]=z;
    }
  }
  if (flag & cmp_C) return y;
  settyp(y,tx);
  if (flag & cmp_IND)
    for (i=1; i<lx; i++) gel(y,i) = utoipos(y[i]);
  else if (tx == t_VECSMALL)
    for (i=1; i<lx; i++) y[i] = x[y[i]];
  else 
    for (i=1; i<lx; i++) gel(y,i) = gcopy(gel(x,y[i]));
  return y;
}

GEN
gen_sort(GEN x, long flag, int (*cmp)(GEN,GEN))
{
  return gen_sort_aux(x, flag, (void *)cmp, cmp_nodata);
}

#define sort_fun(flag) ((flag & cmp_LEX)? &lexcmp: &gcmp)

GEN
gen_vecsort(GEN x, GEN k, long flag)
{
  long i,j,l,t, lx = lg(x), tmp[2];
  struct veccmp_s v;

  if (lx<=2) return gen_sort(x,flag,sort_fun(flag));
  t = typ(k); v.cmp = sort_fun(flag);
  if (t==t_INT)
  {
    gel(tmp,1) = k; k = tmp;
    v.lk = 2;
  }
  else
  {
    if (! is_vec_t(t)) pari_err(talker,"incorrect lextype in vecsort");
    v.lk = lg(k);
  }
  l = 0;
  v.k = (GEN)gpmalloc(v.lk * sizeof(long));
  for (i=1; i<v.lk; i++)
  {
    j = itos(gel(k,i));
    if (j<=0) pari_err(talker,"negative index in vecsort");
    v.k[i]=j; if (j>l) l=j;
  }
  t = typ(x);
  if (! is_matvec_t(t)) pari_err(typeer,"vecsort");
  for (j=1; j<lx; j++)
  {
    t = typ(x[j]);
    if (! is_vec_t(t)) pari_err(typeer,"vecsort");
    if (lg(gel(x,j)) <= l) pari_err(talker,"index too large in vecsort");
  }
  x = gen_sort_aux(x, flag, (void *) &v, veccmp);
  free(v.k); return x;
}

GEN
vecsort0(GEN x, GEN k, long flag)
{
  if (flag < 0 || flag >= cmp_C) pari_err(flagerr,"vecsort");
  if (k) return gen_vecsort(x, k, flag);
  return gen_sort(x, flag, (typ(x) == t_VECSMALL)?
                       pari_compare_small:sort_fun(flag));
}

GEN
vecsort(GEN x, GEN k)
{
  return gen_vecsort(x,k, 0);
}

GEN
sindexsort(GEN x)
{
  return gen_sort(x, cmp_IND | cmp_C, gcmp);
}

GEN
sindexlexsort(GEN x)
{
  return gen_sort(x, cmp_IND | cmp_C, lexcmp);
}

GEN
indexsort(GEN x)
{
  return gen_sort(x, cmp_IND, gcmp);
}

GEN
indexlexsort(GEN x)
{
  return gen_sort(x, cmp_IND, lexcmp);
}

GEN
sort(GEN x)
{
  return gen_sort(x, 0, gcmp);
}

GEN
lexsort(GEN x)
{
  return gen_sort(x, 0, lexcmp);
}

/* index of x in table T, 0 otherwise */
long
tablesearch(GEN T, GEN x, int (*cmp)(GEN,GEN))
{
  long l=1,u=lg(T)-1,i,s;

  while (u>=l)
  {
    i = (l+u)>>1; s = cmp(x,gel(T,i));
    if (!s) return i;
    if (s<0) u=i-1; else l=i+1;
  }
  return 0;
}

/* assume lg(x) = lg(y), x,y in Z^n */
int
cmp_vecint(GEN x, GEN y)
{
  long fl,i, lx = lg(x);
  for (i=1; i<lx; i++)
    if (( fl = cmpii(gel(x,i), gel(y,i)) )) return fl;
  return 0;
}

/* assume x and y come from the same primedec call (uniformizer unique) */
int
cmp_prime_over_p(GEN x, GEN y)
{
  long k = mael(x,4,2) - mael(y,4,2); /* diff. between residue degree */
  return k? ((k > 0)? 1: -1)
          : cmp_vecint(gel(x,2), gel(y,2));
}

int
cmp_prime_ideal(GEN x, GEN y)
{
  int k = cmpii(gel(x,1), gel(y,1));
  return k? k: cmp_prime_over_p(x,y);
}
