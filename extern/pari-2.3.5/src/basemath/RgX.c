/* $Id: RgX.c 11396 2008-12-08 11:04:43Z kb $

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

int 
is_rational(GEN x) { long t = typ(x); return is_rational_t(t); }
int
RgX_is_rational(GEN x)
{
  long i;
  for (i = lg(x)-1; i>1; i--)
    if (!is_rational(gel(x,i))) return 0;
  return 1;
}

/********************************************************************/
/**                                                                **/
/**                        LINEAR ALGEBRA                          **/
/**                                                                **/
/********************************************************************/

/* x non-empty t_MAT, y a compatible zc (dimension > 0). */
static GEN
RgM_zc_mul_i(GEN x, GEN y, long c, long l)
{
  long i, j;
  pari_sp av;
  GEN z = cgetg(l,t_COL), s;

  for (i=1; i<l; i++)
  {
    av = avma; s = gmulgs(gcoeff(x,i,1),y[1]);
    for (j=2; j<c; j++)
       if (y[j]) s = gadd(s, gmulgs(gcoeff(x,i,j),y[j]));
    gel(z,i) = gerepileupto(av,s);
  }
  return z;
}
GEN
RgM_zc_mul(GEN x, GEN y) { return RgM_zc_mul_i(x,y, lg(x), lg(x[1])); }

/* x t_MAT, y a compatible zm (dimension > 0). */
GEN 
RgM_zm_mul(GEN x, GEN y)
{
  long j, c, l = lg(x), ly = lg(y);
  GEN z = cgetg(ly, t_MAT);
  if (l == 1) return z;
  c = lg(x[1]);
  for (j = 1; j < ly; j++) gel(z,j) = RgM_zc_mul_i(x, gel(y,j), l,c);
  return z;
}
static GEN
RgV_zc_mul_i(GEN x, GEN y, long l)
{
  long i;
  GEN z = gen_0;
  pari_sp av = avma;
  for (i = 1; i < l; i++) z = gadd(z, gmulgs(gel(x,i), y[i]));
  return gerepileupto(av, z);
}
GEN
RgV_zc_mul(GEN x, GEN y) { return RgV_zc_mul_i(x, y, lg(x)); }

GEN 
RgV_zm_mul(GEN x, GEN y)
{
  long j, l = lg(x), ly = lg(y);
  GEN z = cgetg(ly, t_VEC);
  for (j = 1; j < ly; j++) gel(z,j) = RgV_zc_mul_i(x, gel(y,j), l);
  return z;
}

/********************************************************************/
/**                                                                **/
/**                         COMPOSITION                            **/
/**                                                                **/
/********************************************************************/

/* evaluate f(x) mod T */
GEN
RgX_RgXQ_compo(GEN f, GEN x, GEN T)
{
  pari_sp av = avma, limit;
  long l;
  GEN y;

  if (typ(f) != t_POL) return gcopy(f);
  l = lg(f)-1; if (l == 1) return zeropol(varn(T));
  y = gel(f,l);
  limit = stack_lim(av, 1);
  for (l--; l>=2; l--)
  {
    y = grem(gadd(gmul(y,x), gel(f,l)), T);
    if (low_stack(limit,stack_lim(av,1)))
    {
      if (DEBUGMEM > 1) pari_warn(warnmem, "RgX_RgXQ_compo");
      y = gerepileupto(av, y);
    }
  }
  return gerepileupto(av, y);
}

/* return (1,...a^l) mod T. Not memory clean */
GEN
RgX_powers(GEN a, GEN T, long l)
{
  long i;
  GEN v;

  if (typ(a) != t_POL) pari_err(typeer,"RgX_powers");
  l += 2;
  v = cgetg(l,t_VEC);
  gel(v,1) = gen_1; if (l == 2) return v;

  if (degpol(a) >= degpol(T)) a = grem(a, T);
  gel(v,2) = a;
  for (i=3; i<l; i++) gel(v,i) = grem(gmul(gel(v,i-1), a), T);
  return v;
}

/* Return P(h * x), not memory clean */
GEN
RgX_unscale(GEN P, GEN h)
{
  long i, l = lg(P);
  GEN hi = gen_1, Q = cgetg(l, t_POL);
  Q[1] = P[1];
  gel(Q,2) = gcopy(gel(P,2));
  for (i=3; i<l; i++)
  {
    hi = gmul(hi,h);
    gel(Q,i) = gmul(gel(P,i), hi);
  }
  return Q;
}

GEN
RgXV_unscale(GEN v, GEN h)
{
  long i, l;
  GEN w;
  if (!h) return v;
  l = lg(v); w = cgetg(l, typ(v));
  for (i=1; i<l; i++) gel(w,i) = RgX_unscale(gel(v,i), h);
  return w;
}

/* Return h^degpol(P) P(x / h), not memory clean */
GEN
RgX_rescale(GEN P, GEN h)
{
  long i, l = lg(P);
  GEN Q = cgetg(l,t_POL), hi = h;
  Q[l-1] = P[l-1];
  for (i=l-2; i>=2; i--)
  {
    gel(Q,i) = gmul(gel(P,i), hi);
    if (i == 2) break;
    hi = gmul(hi,h);
  }
  Q[1] = P[1]; return Q;
}

/********************************************************************/
/**                                                                **/
/**                          CONVERSIONS                           **/
/**                       (not memory clean)                       **/
/**                                                                **/
/********************************************************************/
GEN
RgX_renormalize(GEN x)
{
  long i, lx = lg(x);
  for (i = lx-1; i>1; i--)
    if (! gcmp0(gel(x,i))) break; /* _not_ isexactzero */
  stackdummy((pari_sp)(x + lg(x)), (pari_sp)(x + i+1));
  setlg(x, i+1); setsigne(x, i != 1); return x;
}

GEN
RgV_to_RgX(GEN x, long v)
{
  long i, k = lg(x);
  GEN p;

  while (--k && gcmp0(gel(x,k)));
  if (!k) return zeropol(v);
  i = k+2; p = cgetg(i,t_POL);
  p[1] = evalsigne(1) | evalvarn(v);
  x--; for (k=2; k<i; k++) p[k] = x[k];
  return p;
}

/* return the (N-dimensional) vector of coeffs of p */
GEN
RgX_to_RgV(GEN x, long N)
{
  long i, l;
  GEN z = cgetg(N+1,t_COL);
  if (typ(x) != t_POL)
  {
    gel(z,1) = x;
    for (i=2; i<=N; i++) gel(z,i) = gen_0;
    return z;
  }
  l = lg(x)-1; x++;
  for (i=1; i<l ; i++) z[i]=x[i];
  for (   ; i<=N; i++) gel(z,i) = gen_0;
  return z;
}

/* vector of polynomials (in v) whose coeffs are given by the columns of x */
GEN
RgM_to_RgXV(GEN x, long v)
{
  long j, lx = lg(x);
  GEN y = cgetg(lx, t_VEC);
  for (j=1; j<lx; j++) gel(y,j) = RgV_to_RgX(gel(x,j), v);
  return y;
}

/* matrix whose entries are given by the coeffs of the polynomials in
 * vector v (considered as degree n-1 polynomials) */
GEN
RgXV_to_RgM(GEN v, long n)
{
  long j, N = lg(v);
  GEN y = cgetg(N, t_MAT);
  for (j=1; j<N; j++) gel(y,j) = RgX_to_RgV(gel(v,j), n);
  return y;
}

/* polynomial (in v) of polynomials (in w) whose coeffs are given by the columns of x */
GEN
RgM_to_RgXX(GEN x, long v,long w)
{
  long j, lx = lg(x);
  GEN y = cgetg(lx+1, t_POL);
  y[1]=evalsigne(1) | evalvarn(v);
  y++;
  for (j=1; j<lx; j++) gel(y,j) = RgV_to_RgX(gel(x,j), w);
  return normalizepol_i(--y, lx+1);
}

/* matrix whose entries are given by the coeffs of the polynomial v in
 * two variables (considered as degree n polynomials) */
GEN
RgXX_to_RgM(GEN v, long n)
{
  long j, N = lg(v)-1;
  GEN y = cgetg(N, t_MAT);
  v++;
  for (j=1; j<N; j++) gel(y,j) = RgX_to_RgV(gel(v,j), n);
  return y;
}

/* P(X,Y) --> P(Y,X), n is an upper bound for deg_Y(P) */
GEN
RgXY_swap(GEN x, long n, long w)
{
  long j, lx = lg(x), ly = n+3;
  long v=varn(x);
  GEN y = cgetg(ly, t_POL);
  y[1]=evalsigne(1) | evalvarn(v);
  for (j=2; j<ly; j++)
  {
    long k;
    GEN p1=cgetg(lx,t_POL);
    p1[1] = evalsigne(1) | evalvarn(w);
    for (k=2; k<lx; k++)
      if (j < lg(x[k]))
        gel(p1,k) = gmael(x,k,j);
      else
        gel(p1,k) = gen_0;
    gel(y,j) = normalizepol_i(p1,lx);
  }
  return normalizepol_i(y,ly);
}

/* return (x * X^n). Shallow */
GEN
RgX_shift_shallow(GEN a, long n)
{
  long i, l = lg(a);
  GEN  b;
  if (lg(a) == 2 || !n) return a;
  l += n;
  if (n < 0)
  {
    if (l <= 2) return zeropol(varn(a));
    b = cgetg(l, t_POL); b[1] = a[1];
    a -= n;
    for (i=2; i<l; i++) gel(b,i) = gel(a,i);
  } else {
    b = cgetg(l, t_POL); b[1] = a[1];
    a -= n; n += 2;
    for (i=2; i<n; i++) gel(b,i) = gen_0;
    for (   ; i<l; i++) gel(b,i) = gel(a,i);
  }
  return b;
}
/* return (x * X^n). */
GEN
RgX_shift(GEN a, long n)
{
  long i, l = lg(a);
  GEN  b;
  if (lg(a) == 2 || !n) return gcopy(a);
  l += n;
  if (n < 0)
  {
    if (l <= 2) return zeropol(varn(a));
    b = cgetg(l, t_POL); b[1] = a[1];
    a -= n;
    for (i=2; i<l; i++) gel(b,i) = gcopy(gel(a,i));
  } else {
    b = cgetg(l, t_POL); b[1] = a[1];
    a -= n; n += 2;
    for (i=2; i<n; i++) gel(b,i) = gen_0;
    for (   ; i<l; i++) gel(b,i) = gcopy(gel(a,i));
  }
  return b;
}
GEN
RgX_mulXn(GEN x, long d)
{
  pari_sp av;
  GEN z;
  long v;
  if (d >= 0) return RgX_shift(x, d);
  d = -d;
  v = polvaluation(x, NULL);
  if (v >= d) return RgX_shift(x, -d);
  av = avma;
  z = gred_rfrac_simple( RgX_shift(x, -v),
                         monomial(gen_1, d - v, varn(x)));
  return gerepileupto(av, z);
}

GEN
RgXQC_red(GEN P, GEN T)
{
  long i, l = lg(P);
  GEN Q = cgetg(l, t_COL);
  for (i=1; i<l; i++) gel(Q,i) = grem(gel(P,i), T);
  return Q;
}

GEN
RgXQV_red(GEN P, GEN T)
{
  long i, l = lg(P);
  GEN Q = cgetg(l, t_VEC);
  for (i=1; i<l; i++) gel(Q,i) = grem(gel(P,i), T);
  return Q;
}

GEN
RgXQX_red(GEN P, GEN T)
{
  long i, l = lg(P);
  GEN Q = cgetg(l, t_POL);
  Q[1] = P[1];
  for (i=2; i<l; i++) gel(Q,i) = grem(gel(P,i), T);
  return normalizepol_i(Q, l);
}

/*******************************************************************/
/*                                                                 */
/*                  KARATSUBA (for polynomials)                    */
/*                                                                 */
/*******************************************************************/
GEN
zx_copy_spec(GEN x, long nx)
{
  GEN z = cgetg(nx+2,t_POL);
  long i;
  for (i=0; i<nx; i++) gel(z,i+2) = stoi(x[i]);
  z[1] = evalsigne(1); return z;
}

GEN
RgX_copy_spec(GEN x, long nx)
{
  GEN z = cgetg(nx+2,t_POL);
  long i;
  for (i=0; i<nx; i++) z[i+2] = x[i];
  z[1] = evalsigne(1); return z;
}

/* generic multiplication */

static GEN
addpol(GEN x, GEN y, long lx, long ly)
{
  long i,lz;
  GEN z;

  if (ly>lx) swapspec(x,y, lx,ly);
  lz = lx+2; z = cgetg(lz,t_POL) + 2;
  for (i=0; i<ly; i++) gel(z,i) = gadd(gel(x,i),gel(y,i));
  for (   ; i<lx; i++) z[i]=x[i];
  z -= 2; z[1]=0; return normalizepol_i(z, lz);
}

static GEN
addpolcopy(GEN x, GEN y, long lx, long ly)
{
  long i,lz;
  GEN z;

  if (ly>lx) swapspec(x,y, lx,ly);
  lz = lx+2; z = cgetg(lz,t_POL) + 2;
  for (i=0; i<ly; i++) gel(z,i) = gadd(gel(x,i),gel(y,i));
  for (   ; i<lx; i++) gel(z,i) = gcopy(gel(x,i));
  z -= 2; z[1]=0; return normalizepol_i(z, lz);
}

INLINE GEN
mulpol_limb(GEN x, GEN y, char *ynonzero, long a, long b)
{
  GEN p1 = NULL;
  long i;
  pari_sp av = avma;
  for (i=a; i<b; i++)
    if (ynonzero[i])
    {
      GEN p2 = gmul(gel(y,i),gel(x,-i));
      p1 = p1 ? gadd(p1, p2): p2;
    }
  return p1 ? gerepileupto(av, p1): gen_0;
}

/* assume nx >= ny > 0 */
static GEN
mulpol(GEN x, GEN y, long nx, long ny)
{
  long i,lz,nz;
  GEN z;
  char *p1;

  lz = nx+ny+1; nz = lz-2;
  z = cgetg(lz, t_POL) + 2; /* x:y:z [i] = term of degree i */
  p1 = gpmalloc(ny);
  for (i=0; i<ny; i++)
  {
    p1[i] = !isexactzero(gel(y,i));
    gel(z,i) = mulpol_limb(x+i,y,p1,0,i+1);
  }
  for (  ; i<nx; i++) gel(z,i) = mulpol_limb(x+i,y,p1,0,ny);
  for (  ; i<nz; i++) gel(z,i) = mulpol_limb(x+i,y,p1,i-nx+1,ny);
  free(p1); z -= 2; z[1]=0; return normalizepol_i(z, lz);
}

/* return (x * X^d) + y. Assume d > 0, y != 0 */
GEN
addmulXn(GEN x, GEN y, long d)
{
  GEN xd, yd, zd;
  long a, lz, nx, ny;
  
  if (!signe(x)) return y;
  ny = lgpol(y);
  nx = lgpol(x);
  zd = (GEN)avma;
  x += 2; y += 2; a = ny-d;
  if (a <= 0)
  {
    lz = (a>nx)? ny+2: nx+d+2;
    (void)new_chunk(lz); xd = x+nx; yd = y+ny;
    while (xd > x) gel(--zd,0) = gel(--xd,0);
    x = zd + a;
    while (zd > x) gel(--zd,0) = gen_0;
  }
  else
  {
    xd = new_chunk(d); yd = y+d;
    x = addpol(x,yd, nx,a);
    lz = (a>nx)? ny+2: lg(x)+d;
    x += 2; while (xd > x) *--zd = *--xd;
  }
  while (yd > y) *--zd = *--yd;
  *--zd = evalsigne(1);
  *--zd = evaltyp(t_POL) | evallg(lz); return zd;
}

GEN
addshiftpol(GEN x, GEN y, long d)
{
  long v = varn(x);
  x = addmulXn(x,y,d);
  setvarn(x,v); return x;
}

/* as above, producing a clean malloc */
static GEN
addmulXncopy(GEN x, GEN y, long d)
{
  GEN xd, yd, zd;
  long a, lz, nx, ny;
  
  if (!signe(x)) return gcopy(y);
  nx = lgpol(x);
  ny = lgpol(y);
  zd = (GEN)avma;
  x += 2; y += 2; a = ny-d;
  if (a <= 0)
  {
    lz = nx+d+2;
    (void)new_chunk(lz); xd = x+nx; yd = y+ny;
    while (xd > x) gel(--zd,0) = gcopy(gel(--xd,0));
    x = zd + a;
    while (zd > x) gel(--zd,0) = gen_0;
  }
  else
  {
    xd = new_chunk(d); yd = y+d;
    x = addpolcopy(x,yd, nx,a);
    lz = (a>nx)? ny+2: lg(x)+d;
    x += 2; while (xd > x) *--zd = *--xd;
  }
  while (yd > y) gel(--zd,0) = gcopy(gel(--yd,0));
  *--zd = evalsigne(1);
  *--zd = evaltyp(t_POL) | evallg(lz); return zd;
}

/* shift polynomial in place. assume v free cells have been left before x */
static GEN
shiftpol_ip(GEN x, long v)
{
  long i, lx;
  GEN y, z;
  if (!v) return x;
  lx = lg(x);
  if (lx == 2) return x;
  y = x + v;
  z = x + lx;
  /* stackdummy from normalizepol: move it up */
  if (lg(z) != v) x[lx + v] = z[0];
  for (i = lx-1; i >= 2; i--) y[i] = x[i];
  for (i = v+1;  i >= 2; i--) gel(x,i) = gen_0;
  /* leave x[1] alone: it is correct */
  x[0] = evaltyp(t_POL) | evallg(lx+v); return x;
}

/* fast product (Karatsuba) of polynomials a,b. These are not real GENs, a+2,
 * b+2 were sent instead. na, nb = number of terms of a, b.
 * Only c, c0, c1, c2 are genuine GEN.
 */
GEN
RgX_mulspec(GEN a, GEN b, long na, long nb)
{
  GEN a0, c, c0;
  long n0, n0a, i, v = 0;
  pari_sp av;

  while (na && isexactzero(gel(a,0))) { a++; na--; v++; }
  while (nb && isexactzero(gel(b,0))) { b++; nb--; v++; }
  if (na < nb) swapspec(a,b, na,nb);
  if (!nb) return zeropol(0);

  if (v) (void)cgetg(v,t_VECSMALL); /* v gerepile-safe cells for shiftpol_ip */
  if (nb < RgX_MUL_LIMIT)
    return shiftpol_ip(mulpol(a,b,na,nb), v);
  i = (na>>1); n0 = na-i; na = i;
  av = avma; a0 = a+n0; n0a = n0;
  while (n0a && isexactzero(gel(a,n0a-1))) n0a--;

  if (nb > n0)
  {
    GEN b0,c1,c2;
    long n0b;

    nb -= n0; b0 = b+n0; n0b = n0;
    while (n0b && isexactzero(gel(b,n0b-1))) n0b--;
    c = RgX_mulspec(a,b,n0a,n0b);
    c0 = RgX_mulspec(a0,b0, na,nb);

    c2 = addpol(a0,a, na,n0a);
    c1 = addpol(b0,b, nb,n0b);

    c1 = RgX_mulspec(c1+2,c2+2, lgpol(c1),lgpol(c2));
    c2 = gadd(c1, gneg_i(gadd(c0,c)));
    c0 = addmulXn(c0, c2, n0);
  }
  else
  {
    c = RgX_mulspec(a,b,n0a,nb);
    c0 = RgX_mulspec(a0,b,na,nb);
  }
  c0 = addmulXncopy(c0,c,n0);
  return shiftpol_ip(gerepileupto(av,c0), v);
}

static GEN
sqrpol(GEN x, long nx)
{
  long i, j, l, lz, nz;
  pari_sp av;
  GEN p1,z;
  char *p2;

  if (!nx) return zeropol(0);
  lz = (nx << 1) + 1, nz = lz-2;
  z = cgetg(lz,t_POL) + 2;
  p2 = gpmalloc(nx);
  for (i=0; i<nx; i++)
  {
    p2[i] = !isexactzero(gel(x,i));
    p1=gen_0; av=avma; l=(i+1)>>1;
    for (j=0; j<l; j++)
      if (p2[j] && p2[i-j])
        p1 = gadd(p1, gmul(gel(x,j),gel(x,i-j)));
    p1 = gshift(p1,1);
    if ((i&1) == 0 && p2[i>>1])
      p1 = gadd(p1, gsqr((GEN)x[i>>1]));
    gel(z,i) = gerepileupto(av,p1);
  }
  for (  ; i<nz; i++)
  {
    p1=gen_0; av=avma; l=(i+1)>>1;
    for (j=i-nx+1; j<l; j++)
      if (p2[j] && p2[i-j])
        p1 = gadd(p1, gmul(gel(x,j),gel(x,i-j)));
    p1 = gshift(p1,1);
    if ((i&1) == 0 && p2[i>>1])
      p1 = gadd(p1, gsqr((GEN)x[i>>1]));
    gel(z,i) = gerepileupto(av,p1);
  }
  free(p2); z -= 2; z[1]=0; return normalizepol_i(z,lz);
}

GEN
RgX_sqrspec(GEN a, long na)
{
  GEN a0, c, c0, c1;
  long n0, n0a, i, v = 0;
  pari_sp av;

  while (na && isexactzero(gel(a,0))) { a++; na--; v += 2; }
  if (v) (void)cgetg(v, t_VECSMALL);
  if (na<RgX_SQR_LIMIT) return shiftpol_ip(sqrpol(a,na), v);
  i = (na>>1); n0 = na-i; na = i;
  av = avma; a0 = a+n0; n0a = n0;
  while (n0a && isexactzero(gel(a,n0a-1))) n0a--;

  c = RgX_sqrspec(a,n0a);
  c0 = RgX_sqrspec(a0,na);
  c1 = gmul2n(RgX_mulspec(a0,a, na,n0a), 1);
  c0 = addmulXn(c0,c1, n0);
  c0 = addmulXncopy(c0,c,n0);
  return shiftpol_ip(gerepileupto(av,c0), v);
}
GEN
RgX_mul(GEN x,GEN y)
{
  GEN z = RgX_mulspec(y+2, x+2, lgpol(y), lgpol(x));
  setvarn(z,varn(x)); return z;
}
GEN
RgX_sqr(GEN x)
{
  GEN z = RgX_sqrspec(x+2, lgpol(x));
  setvarn(z,varn(x)); return z;
}

/*******************************************************************/
/*                                                                 */
/*                               DIVISION                          */
/*                                                                 */
/*******************************************************************/
GEN
RgX_Rg_div(GEN x, GEN y) {
  long i, lx = lg(x);
  GEN z = cgetg_copy(lx, x); z[1] = x[1];
  for (i=2; i<lx; i++) gel(z,i) = gdiv(gel(x,i),y);
  return normalizepol_i(z, lx);
}
GEN
RgX_div_by_X_x(GEN a, GEN x, GEN *r)
{
  long l = lg(a), i;
  GEN a0, z0, z = cgetg(l-1, t_POL);
  z[1] = a[1];
  a0 = a + l-1;
  z0 = z + l-2; *z0 = *a0--;
  for (i=l-3; i>1; i--) /* z[i] = a[i+1] + x*z[i+1] */
  {
    GEN t = gadd(gel(a0--,0), gmul(x, gel(z0--,0)));
    gel(z0,0) = t;
  }
  if (r) *r = gadd(gel(a0,0), gmul(x, gel(z0,0)));
  return z;
}
/* Polynomial division x / y:
 *   if z = ONLY_REM  return remainder, otherwise return quotient
 *   if z != NULL set *z to remainder
 *   *z is the last object on stack (and thus can be disposed of with cgiv
 *   instead of gerepile) */
/* assume, typ(x) = typ(y) = t_POL, same variable */
GEN
RgX_divrem(GEN x, GEN y, GEN *pr)
{
  pari_sp avy, av, av1;
  long dx,dy,dz,i,j,sx,lr;
  GEN z,p1,rem,y_lead,mod;
  GEN (*f)(GEN,GEN);

  if (!signe(y)) pari_err(gdiver);

  dy = degpol(y);
  y_lead = gel(y,dy+2);
  if (gcmp0(y_lead)) /* normalize denominator if leading term is 0 */
  {
    pari_warn(warner,"normalizing a polynomial with 0 leading term");
    for (dy--; dy>=0; dy--)
    {
      y_lead = gel(y,dy+2);
      if (!gcmp0(y_lead)) break;
    }
  }
  if (!dy) /* y is constant */
  {
    if (pr && pr != ONLY_DIVIDES)
    {
      if (pr == ONLY_REM) return zeropol(varn(x));
      *pr = zeropol(varn(x));
    }
    return gdiv(x, y_lead);
  }
  dx = degpol(x);
  if (dx < dy)
  {
    if (pr)
    {
      if (pr == ONLY_DIVIDES) return gcmp0(x)? gen_0: NULL;
      if (pr == ONLY_REM) return gcopy(x);
      *pr = gcopy(x);
    }
    return zeropol(varn(x));
  }

  /* x,y in R[X], y non constant */
  dz = dx-dy; av = avma;
  p1 = new_chunk(dy+3);
  for (i=2; i<dy+3; i++)
  {
    GEN p2 = gel(y,i);
    /* gneg to avoid gsub's later on */
    gel(p1,i) = isexactzero(p2)? NULL: gneg_i(p2);
  }
  y = p1;
  switch(typ(y_lead))
  {
    case t_INTMOD:
    case t_POLMOD: y_lead = ginv(y_lead);
      f = gmul; mod = gmodulo(gen_1, gel(y_lead,1));
      break;
    default: if (gcmp1(y_lead)) y_lead = NULL;
      f = gdiv; mod = NULL;
  }
  avy=avma;
  z = cgetg(dz+3,t_POL); z[1] = x[1];
  x += 2; y += 2; z += 2;

  p1 = gel(x,dx);
  gel(z,dz) = y_lead? f(p1,y_lead): gcopy(p1);
  for (i=dx-1; i>=dy; i--)
  {
    av1=avma; p1=gel(x,i);
    for (j=i-dy+1; j<=i && j<=dz; j++)
      if (y[i-j] && gel(z,j) != gen_0) p1 = gadd(p1, gmul(gel(z,j),gel(y,i-j)));
    if (y_lead) p1 = f(p1,y_lead);

    if (isexactzero(p1)) { avma=av1; p1 = gen_0; }
    else
      p1 = avma==av1? gcopy(p1): gerepileupto(av1,p1);
    gel(z,i-dy) = p1;
  }
  if (!pr) return gerepileupto(av,z-2);

  rem = (GEN)avma; av1 = (pari_sp)new_chunk(dx+3);
  for (sx=0; ; i--)
  {
    p1 = gel(x,i);
    /* we always enter this loop at least once */
    for (j=0; j<=i && j<=dz; j++)
      if (y[i-j] && gel(z,j) != gen_0) p1 = gadd(p1, gmul(gel(z,j),gel(y,i-j)));
    if (mod && avma==av1) p1 = gmul(p1,mod);
    if (!gcmp0(p1)) { sx = 1; break; } /* remainder is non-zero */
    if (!isinexactreal(p1) && !isexactzero(p1)) break;
    if (!i) break;
    avma=av1;
  }
  if (pr == ONLY_DIVIDES)
  {
    if (sx) { avma=av; return NULL; }
    avma = (pari_sp)rem;
    return gerepileupto(av,z-2);
  }
  lr=i+3; rem -= lr;
  if (avma==av1) { avma = (pari_sp)rem; p1 = gcopy(p1); }
  else p1 = gerepileupto((pari_sp)rem,p1);
  rem[0] = evaltyp(t_POL) | evallg(lr);
  rem[1] = z[-1];
  rem += 2;
  gel(rem,i) = p1;
  for (i--; i>=0; i--)
  {
    av1=avma; p1 = gel(x,i);
    for (j=0; j<=i && j<=dz; j++)
      if (y[i-j] && gel(z,j) != gen_0) p1 = gadd(p1, gmul(gel(z,j),gel(y,i-j)));
    if (mod && avma==av1) p1 = gmul(p1,mod);
    gel(rem,i) = avma==av1? gcopy(p1):gerepileupto(av1,p1);
  }
  rem -= 2;
  if (!sx) (void)normalizepol_i(rem, lr);
  if (pr == ONLY_REM) return gerepileupto(av,rem);
  z -= 2;
  {
    GEN *gptr[2]; gptr[0]=&z; gptr[1]=&rem;
    gerepilemanysp(av,avy,gptr,2); *pr = rem; return z;
  }
}

/* x and y in (R[Y]/T)[X]  (lifted), T in R[Y]. y preferably monic */
GEN
RgXQX_divrem(GEN x, GEN y, GEN T, GEN *pr)
{
  long vx, dx, dy, dz, i, j, sx, lr;
  pari_sp av0, av, tetpil;
  GEN z,p1,rem,lead;

  if (!signe(y)) pari_err(gdiver);
  vx = varn(x);
  dx = degpol(x);
  dy = degpol(y);
  if (dx < dy)
  {
    if (pr)
    {
      av0 = avma; x = RgXQX_red(x, T);
      if (pr == ONLY_DIVIDES) { avma=av0; return signe(x)? NULL: gen_0; }
      if (pr == ONLY_REM) return x;
      *pr = x;
    }
    return zeropol(vx);
  }
  lead = leading_term(y);
  if (!dy) /* y is constant */
  {
    if (pr && pr != ONLY_DIVIDES)
    {
      if (pr == ONLY_REM) return zeropol(vx);
      *pr = zeropol(vx);
    }
    if (gcmp1(lead)) return gcopy(x);
    av0 = avma; x = gmul(x, ginvmod(lead,T)); tetpil = avma;
    return gerepile(av0,tetpil,RgXQX_red(x,T));
  }
  av0 = avma; dz = dx-dy;
  lead = gcmp1(lead)? NULL: gclone(ginvmod(lead,T));
  avma = av0;
  z = cgetg(dz+3,t_POL); z[1] = x[1];
  x += 2; y += 2; z += 2;

  p1 = gel(x,dx); av = avma;
  gel(z,dz) = lead? gerepileupto(av, grem(gmul(p1,lead), T)): gcopy(p1);
  for (i=dx-1; i>=dy; i--)
  {
    av=avma; p1=gel(x,i);
    for (j=i-dy+1; j<=i && j<=dz; j++)
      p1 = gsub(p1, gmul(gel(z,j),gel(y,i-j)));
    if (lead) p1 = gmul(grem(p1, T), lead);
    tetpil=avma; gel(z,i-dy) = gerepile(av,tetpil, grem(p1, T));
  }
  if (!pr) { if (lead) gunclone(lead); return z-2; }

  rem = (GEN)avma; av = (pari_sp)new_chunk(dx+3);
  for (sx=0; ; i--)
  {
    p1 = gel(x,i);
    for (j=0; j<=i && j<=dz; j++)
      p1 = gsub(p1, gmul(gel(z,j),gel(y,i-j)));
    tetpil=avma; p1 = grem(p1, T); if (!gcmp0(p1)) { sx = 1; break; }
    if (!i) break;
    avma=av;
  }
  if (pr == ONLY_DIVIDES)
  {
    if (lead) gunclone(lead);
    if (sx) { avma=av0; return NULL; }
    avma = (pari_sp)rem; return z-2;
  }
  lr=i+3; rem -= lr;
  rem[0] = evaltyp(t_POL) | evallg(lr);
  rem[1] = z[-1];
  p1 = gerepile((pari_sp)rem,tetpil,p1);
  rem += 2; gel(rem,i) = p1;
  for (i--; i>=0; i--)
  {
    av=avma; p1 = gel(x,i);
    for (j=0; j<=i && j<=dz; j++)
      p1 = gsub(p1, gmul(gel(z,j),gel(y,i-j)));
    tetpil=avma; gel(rem,i) = gerepile(av,tetpil, grem(p1, T));
  }
  rem -= 2;
  if (lead) gunclone(lead);
  if (!sx) (void)normalizepol_i(rem, lr);
  if (pr == ONLY_REM) return gerepileupto(av0,rem);
  *pr = rem; return z-2;
}

GEN
RgXQX_mul(GEN x, GEN y, GEN T)
{
  return RgXQX_red(RgX_mul(x,y), T);
}
GEN
RgX_Rg_mul(GEN y, GEN x) {
  long i, ly;
  GEN z;
  if (isexactzero(x)) { long vy = varn(y); return zeropol(vy); }
  ly = lg(y);
  z = cgetg(ly,t_POL); z[1] = y[1];
  if (ly == 2) return z;
  for (i = 2; i < ly; i++) gel(z,i) = gmul(x,gel(y,i));
  return normalizepol_i(z,ly);
}
GEN
RgXQX_RgXQ_mul(GEN x, GEN y, GEN T)
{
  return RgXQX_red(RgX_Rg_mul(x,y), T);
}
GEN
RgXQX_sqr(GEN x, GEN T)
{
  return RgXQX_red(RgX_sqr(x), T);
}

GEN
RgXQ_mul(GEN y, GEN x, GEN T) { return RgX_rem(RgX_mul(y, x), T); }
GEN
RgXQ_sqr(GEN x, GEN T) { return RgX_rem(RgX_sqr(x), T); }

static GEN
_sqr(void *data, GEN x) { return RgXQ_sqr(x, (GEN)data); }
static GEN
_mul(void *data, GEN x, GEN y) { return RgXQ_mul(x,y, (GEN)data); }

/* x,T in Rg[X], n in N, compute lift(x^n mod T)) */
GEN
RgXQ_u_pow(GEN x, ulong n, GEN T)
{
  pari_sp av;
  GEN y;

  if (!n) return pol_1[varn(x)];
  if (n == 1) return gcopy(x);
  av = avma;
  y = leftright_pow_u(x, n, (void*)T, &_sqr, &_mul);
  return gerepileupto(av, y);
}

/* generates the list of powers of x of degree 0,1,2,...,l*/
GEN
RgXQ_powers(GEN x, long l, GEN T)
{
  GEN V=cgetg(l+2,t_VEC);
  long i;
  gel(V,1) = pol_1[varn(T)]; if (l==0) return V;
  gel(V,2) = gcopy(x);       if (l==1) return V;
  gel(V,3) = RgXQ_sqr(x,T);
  if ((degpol(x)<<1) < degpol(T)) {
    for(i = 4; i < l+2; i++)
      gel(V,i) = RgXQ_mul(gel(V,i-1),x,T);
  } else { /* use squarings if degree(x) is large */
    for(i = 4; i < l+2; i++)
      gel(V,i) = (i&1)? RgXQ_sqr(gel(V, (i+1)>>1),T)
                      : RgXQ_mul(gel(V, i-1),x,T);
  }
  return V;
}

GEN
RgXQ_matrix_pow(GEN y, long n, long m, GEN P)
{
  return RgXV_to_RgM(RgXQ_powers(y,m-1,P),n);
}

GEN
RgXQ_minpoly_naive(GEN y, GEN P)
{
  pari_sp ltop=avma;
  long n=lgpol(P);
  GEN M=ker(RgXQ_matrix_pow(y,n,n,P));
  M=content(RgM_to_RgXV(M,varn(P)));
  return gerepileupto(ltop,M);
}

