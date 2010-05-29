/* $Id: Flx.c 8433 2007-03-22 22:46:18Z bill $

Copyright (C) 2004  The PARI group.

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

/* Not so fast arithmetic with polynomials with small coefficients. */

/***********************************************************************/
/**                                                                   **/
/**               Flx                                                 **/
/**                                                                   **/
/***********************************************************************/
/* Flx objects are defined as follows:
   Let l an ulong. An Flx is a t_VECSMALL:
   x[0] = codeword
   x[1] = evalvarn(variable number)  (signe is not stored).
   x[2] = a_0 x[3] = a_1, etc.
   With 0 <= a_i < l
   
   signe(x) is not valid. Use degpol(x)>=0 instead.
*/

#define both_odd(x,y) ((x)&(y)&1)

/***********************************************************************/
/**                                                                   **/
/**          Conversion from Flx                                      **/
/**                                                                   **/
/***********************************************************************/

GEN
Flx_to_ZX(GEN z)
{
  long i, l = lg(z);
  GEN x = cgetg(l,t_POL);
  for (i=2; i<l; i++) gel(x,i) = utoi(z[i]);
  x[1] = evalsigne(l-2!=0)| z[1]; return x;
}

GEN
Flv_to_ZV(GEN z)
{
  long i, l = lg(z);
  GEN x = cgetg(l, t_VEC);
  for (i=1; i<l; i++) gel(x,i) = utoi(z[i]);
  return x;
}

GEN
Flc_to_ZC(GEN z)
{
  long i, l = lg(z);
  GEN x = cgetg(l,t_COL);
  for (i=1; i<l; i++) gel(x,i) = utoi(z[i]);
  return x;
}

GEN
Flm_to_ZM(GEN z)
{
  long i, l = lg(z);
  GEN x = cgetg(l,t_MAT);
  for (i=1; i<l; i++) gel(x,i) = Flc_to_ZC(gel(z,i));
  return x;
}

/* same as Flx_to_ZX, in place */
GEN
Flx_to_ZX_inplace(GEN z)
{
  long i, l = lg(z);
  for (i=2; i<l; i++) gel(z,i) = utoi(z[i]);
  settyp(z, t_POL); z[1]=evalsigne(l-2!=0)|z[1]; return z;
}

/*Flx_to_Flv=zx_to_zv*/
GEN
Flx_to_Flv(GEN x, long N)
{
  long i, l;
  GEN z = cgetg(N+1,t_VECSMALL);
  if (typ(x) != t_VECSMALL) pari_err(typeer,"Flx_to_Flv");
  l = lg(x)-1; x++;
  for (i=1; i<l ; i++) z[i]=x[i];
  for (   ; i<=N; i++) z[i]=0;
  return z;
}

/*Flv_to_Flx=zv_to_zx*/
GEN
Flv_to_Flx(GEN x, long vs)
{
  long i, l=lg(x)+1;
  GEN z = cgetg(l,t_VECSMALL); z[1]=vs;
  x--; 
  for (i=2; i<l ; i++) z[i]=x[i];
  return Flx_renormalize(z,l);
}

/*Flm_to_FlxV=zm_to_zxV*/
GEN
Flm_to_FlxV(GEN x, long sv)
{
  long j, lx = lg(x);
  GEN y = cgetg(lx, t_VEC);
  for (j=1; j<lx; j++) gel(y,j) = Flv_to_Flx(gel(x,j), sv);
  return y;
}

/*FlxC_to_ZXC=zxC_to_ZXC*/
GEN
FlxC_to_ZXC(GEN x)
{
  long i, l=lg(x);
  GEN z = cgetg(l,t_COL); 
  for (i=1; i<l ; i++) gel(z,i) = Flx_to_ZX(gel(x,i));
  return z;
}

/*FlxM_to_ZXM=zxM_to_ZXM*/
GEN
FlxM_to_ZXM(GEN z)
{
  long i, l = lg(z);
  GEN x = cgetg(l,t_MAT);
  for (i=1; i<l; i++) gel(x,i) = FlxC_to_ZXC(gel(z,i));
  return x;
}

/***********************************************************************/
/**                                                                   **/
/**          Conversion to Flx                                        **/
/**                                                                   **/
/***********************************************************************/

/*sv is a evalvarn()*/
/*zero_Flx=zero_zx*/
GEN
zero_Flx(long sv)
{
  GEN x = cgetg(2, t_VECSMALL);
  x[1] = sv; return x;
}

/* polx_Flx=polx_zx*/

GEN
polx_Flx(long sv)
{
  GEN z = cgetg(4, t_VECSMALL);
  z[1] = sv;
  z[2] = 0;
  z[3] = 1;
  return z;
}

/* Take an integer and return a scalar polynomial mod p,
 * with evalvarn=vs */

GEN
Fl_to_Flx(ulong x, long sv)
{
  GEN z;
  if (!x) return zero_Flx(sv);
  z = cgetg(3, t_VECSMALL);
  z[1] = sv;
  z[2] = x; return z;
}

GEN
Z_to_Flx(GEN x, ulong p, long v)
{
  GEN z;
  long sv=evalvarn(v);
  z = cgetg(3, t_VECSMALL);
  z[1] = sv;
  z[2] = umodiu(x,p); 
  if (!z[2]) {avma = (pari_sp)(z + lg(z)); z = zero_Flx(sv);}
  return z;
}

/* return x[0 .. dx] mod p as t_VECSMALL. Assume x a t_POL*/
GEN
ZX_to_Flx(GEN x, ulong p)
{
  long i;
  long lx = lg(x); 
  GEN a = cgetg(lx, t_VECSMALL);
  a[1]=((ulong)x[1])&VARNBITS;
  for (i=2; i<lx; i++) a[i] = umodiu(gel(x,i), p);
  return Flx_renormalize(a,lx);
}

GEN
ZV_to_Flv(GEN x, ulong p)
{
  long i, n = lg(x);
  GEN y = cgetg(n,t_VECSMALL);
  for (i=1; i<n; i++) y[i] = umodiu(gel(x,i), p);
  return y;
}

GEN
ZM_to_Flm(GEN x, ulong p)
{
  long j,n = lg(x);
  GEN y = cgetg(n,t_MAT);
  if (n == 1) return y;
  for (j=1; j<n; j++) gel(y,j) = ZV_to_Flv(gel(x,j), p);
  return y;
}

/***********************************************************************/
/**                                                                   **/
/**          Basic operation on Flx                                   **/
/**                                                                   **/
/***********************************************************************/
/*Similar to normalizepol, in place*/
/*FIXME: should be zx_renormalize*/
GEN
Flx_renormalize(GEN /*in place*/ x, long lx)
{
  long i;
  for (i = lx-1; i>1; i--)
    if (x[i]) break;
  stackdummy((pari_sp)(x + lg(x)), (pari_sp)(x + i+1));
  setlg(x, i+1); return x;
}

/*Do not renormalize. Must not use z[0],z[1]*/
static GEN
Flx_red_lg_i(GEN z, long l, ulong p)
{
  long i;
  ulong *y=(ulong *)z;
  GEN x = cgetg(l, t_VECSMALL);
  for (i=2; i<l; i++) x[i] = y[i]%p;
  return x; 
}

GEN
Flx_red(GEN z, ulong p)
{
  long l = lg(z);
  GEN x = Flx_red_lg_i(z,l,p); 
  x[1] = z[1]; 
  return Flx_renormalize(x,l);
}

GEN
Flx_addspec(GEN x, GEN y, ulong p, long lx, long ly)
{
  long i,lz;
  GEN z;

  if (ly>lx) swapspec(x,y, lx,ly);
  lz = lx+2; z = cgetg(lz, t_VECSMALL) + 2;
  for (i=0; i<ly; i++) z[i] = Fl_add(x[i], y[i], p);
  for (   ; i<lx; i++) z[i] = x[i];
  z -= 2; return Flx_renormalize(z, lz);
}

GEN
Flx_add(GEN x, GEN y, ulong p)
{
  long i,lz;
  GEN z; 
  long lx=lg(x);
  long ly=lg(y);
  if (ly>lx) swapspec(x,y, lx,ly);
  lz = lx; z = cgetg(lz, t_VECSMALL); z[1]=x[1];
  for (i=2; i<ly; i++) z[i] = Fl_add(x[i], y[i], p);
  for (   ; i<lx; i++) z[i] = x[i];
  return Flx_renormalize(z, lz);
}

GEN
Flx_subspec(GEN x, GEN y, ulong p, long lx, long ly)
{
  long i,lz;
  GEN z;

  if (ly <= lx)
  {
    lz = lx+2; z = cgetg(lz, t_VECSMALL)+2;
    for (i=0; i<ly; i++) z[i] = Fl_sub(x[i],y[i],p);
    for (   ; i<lx; i++) z[i] = x[i];
  }
  else
  {
    lz = ly+2; z = cgetg(lz, t_VECSMALL)+2;
    for (i=0; i<lx; i++) z[i] = Fl_sub(x[i],y[i],p);
    for (   ; i<ly; i++) z[i] = y[i]? (long)(p - y[i]): y[i];
  }
 return Flx_renormalize(z-2, lz);
}

GEN
Flx_sub(GEN x, GEN y, ulong p)
{
  long i,lz,lx = lg(x), ly = lg(y);
  GEN z;

  if (ly <= lx)
  {
    lz = lx; z = cgetg(lz, t_VECSMALL);
    for (i=2; i<ly; i++) z[i] = Fl_sub(x[i],y[i],p);
    for (   ; i<lx; i++) z[i] = x[i];
  }
  else
  {
    lz = ly; z = cgetg(lz, t_VECSMALL);
    for (i=2; i<lx; i++) z[i] = Fl_sub(x[i],y[i],p);
    for (   ; i<ly; i++) z[i] = y[i]? (long)(p - y[i]): y[i];
  }
  z[1]=x[1]; return Flx_renormalize(z, lz);
}

GEN
Flx_negspec(GEN x, ulong p, long l)
{
  long i;
  GEN z = cgetg(l+2, t_VECSMALL) + 2;
  for (i=0; i<l; i++) z[i] = Fl_neg(x[i], p);
  return z-2;
}


GEN
Flx_neg(GEN x, ulong p)
{
  GEN z = Flx_negspec(x+2, p, lgpol(x));
  z[1] = x[1];
  return z;
}

GEN
Flx_neg_inplace(GEN x, ulong p)
{
  long i, l = lg(x);
  for (i=2; i<l; i++)
    if (x[i]) x[i] = p - x[i];
  return x;
}

GEN
Flx_Fl_mul(GEN y, ulong x, ulong p)
{
  GEN z;
  long i, l;
  if (!x) return zero_Flx(y[1]);
  l = lg(y); z = cgetg(l, t_VECSMALL); z[1]=y[1]; 
  if (HIGHWORD(x | p))
    for(i=2; i<l; i++) z[i] = Fl_mul(y[i], x, p);
  else
    for(i=2; i<l; i++) z[i] = (y[i] * x) % p;
  return z;
}

/*
 * Return a*x^n
 */

GEN
Flx_shift(GEN a, long n)
{
  long i, l = lg(a);
  GEN  b;
  if (l==2) return vecsmall_copy(a);
  b = cgetg(l+n, t_VECSMALL);
  b[1] = a[1];
  for (i=0; i<n; i++) b[2+i] = 0;
  for (i=2; i<l; i++) b[i+n] = a[i];
  return b;
}

GEN
Flx_normalize(GEN z, ulong p)
{
  long l = lg(z)-1;
  ulong p1 = z[l]; /* leading term */
  if (p1 == 1) return z;
  return Flx_Fl_mul(z, Fl_inv(p1,p), p);
}

/* return (x * X^d) + y. Assume d > 0, x > 0 and y >= 0 */
GEN
Flx_addshift(GEN x, GEN y, ulong p, long d)
{
  GEN xd,yd,zd = (GEN)avma;
  long a,lz,ny = lgpol(y), nx = lgpol(x);
  long vs = x[1];

  x += 2; y += 2; a = ny-d;
  if (a <= 0)
  {
    lz = (a>nx)? ny+2: nx+d+2;
    (void)new_chunk(lz); xd = x+nx; yd = y+ny;
    while (xd > x) *--zd = *--xd;
    x = zd + a;
    while (zd > x) *--zd = 0;
  }
  else
  {
    xd = new_chunk(d); yd = y+d;
    x = Flx_addspec(x,yd,p, nx,a);
    lz = (a>nx)? ny+2: lg(x)+d;
    x += 2; while (xd > x) *--zd = *--xd;
  }
  while (yd > y) *--zd = *--yd;
  *--zd = vs;
  *--zd = evaltyp(t_VECSMALL) | evallg(lz); return zd;
}

/* shift polynomial + gerepile */
/* Do not set evalvarn*/
static GEN
Flx_shiftip(pari_sp av, GEN x, long v)
{
  long i, lx = lg(x), ly;
  GEN y;
  if (!v || lx==2) return gerepileuptoleaf(av, x);
  avma = av; ly = lx + v;
  x += lx; y = new_chunk(ly) + ly; /*cgetg could overwrite x!*/
  for (i = 2; i<lx; i++) *--y = *--x;
  for (i = 0; i< v; i++) *--y = 0;
  y -= 2; y[0] = evaltyp(t_VECSMALL) | evallg(ly);
  return y;
}

INLINE long
maxlenghtcoeffpol (ulong p, long n)
{
  pari_sp ltop=avma;
  long l;
  GEN z=muluu(p-1,p-1);
  z=mulis(z,n);
  l=lgefint(z)-2;
  avma=ltop;
  return l;
}

INLINE ulong
Flx_mullimb_ok(GEN x, GEN y, ulong p, long a, long b)
{ /* Assume OK_ULONG*/
  ulong p1 = 0;
  long i;
  for (i=a; i<b; i++)
    if (y[i])
    {
      p1 += y[i] * x[-i];
      if (p1 & HIGHBIT) p1 %= p;
    }
  return p1 % p;
}

INLINE ulong
Flx_mullimb(GEN x, GEN y, ulong p, long a, long b)
{
  ulong p1 = 0;
  long i;
  for (i=a; i<b; i++)
    if (y[i])
      p1 = Fl_add(p1, Fl_mul(y[i],x[-i],p), p);
  return p1;
}

/* assume nx >= ny > 0 */
static GEN
Flx_mulspec_basecase(GEN x, GEN y, ulong p, long nx, long ny)
{
  long i,lz,nz;
  GEN z;

  lz = nx+ny+1; nz = lz-2;
  z = cgetg(lz, t_VECSMALL) + 2; /* x:y:z [i] = term of degree i */
  if (u_OK_ULONG(p))
  {
    for (i=0; i<ny; i++)z[i] = Flx_mullimb_ok(x+i,y,p,0,i+1);
    for (  ; i<nx; i++) z[i] = Flx_mullimb_ok(x+i,y,p,0,ny);
    for (  ; i<nz; i++) z[i] = Flx_mullimb_ok(x+i,y,p,i-nx+1,ny);
  }
  else
  {
    for (i=0; i<ny; i++)z[i] = Flx_mullimb(x+i,y,p,0,i+1);
    for (  ; i<nx; i++) z[i] = Flx_mullimb(x+i,y,p,0,ny);
    for (  ; i<nz; i++) z[i] = Flx_mullimb(x+i,y,p,i-nx+1,ny);
  }
  z -= 2; return z; 
}

INLINE GEN
Flx_mulspec_mulii(GEN a, GEN b, ulong p, long na, long nb)
{
  GEN z=muliispec(a,b,na,nb);
  return Flx_red_lg_i(z,lgefint(z),p);
}

/* fast product (Karatsuba) of polynomials a,b. These are not real GENs, a+2,
 * b+2 were sent instead. na, nb = number of terms of a, b.
 * Only c, c0, c1, c2 are genuine GEN.
 */
GEN
Flx_mulspec(GEN a, GEN b, ulong p, long na, long nb)
{
  GEN a0,c,c0;
  long n0, n0a, i, v = 0;
  pari_sp av;

  while (na && !a[0]) { a++; na--; v++; }
  while (nb && !b[0]) { b++; nb--; v++; }
  if (na < nb) swapspec(a,b, na,nb);
  if (!nb) return zero_Flx(0);

  av = avma;
  if (na>30 && maxlenghtcoeffpol(p,nb)==1)
    return Flx_shiftip(av,Flx_mulspec_mulii(a,b,p,na,nb), v);
  if (nb < Flx_MUL_LIMIT)
    return Flx_shiftip(av,Flx_mulspec_basecase(a,b,p,na,nb), v);
  i=(na>>1); n0=na-i; na=i;
  a0=a+n0; n0a=n0;
  while (n0a && !a[n0a-1]) n0a--;

  if (nb > n0)
  {
    GEN b0,c1,c2;
    long n0b;

    nb -= n0; b0 = b+n0; n0b = n0;
    while (n0b && !b[n0b-1]) n0b--;
    c =  Flx_mulspec(a,b,p,n0a,n0b);
    c0 = Flx_mulspec(a0,b0,p,na,nb);

    c2 = Flx_addspec(a0,a,p,na,n0a);
    c1 = Flx_addspec(b0,b,p,nb,n0b);

    c1 = Flx_mul(c1,c2,p);
    c2 = Flx_add(c0,c,p);

    c2 = Flx_neg_inplace(c2,p);
    c2 = Flx_add(c1,c2,p);
    c0 = Flx_addshift(c0,c2 ,p, n0);
  }
  else
  {
    c  = Flx_mulspec(a,b,p,n0a,nb);
    c0 = Flx_mulspec(a0,b,p,na,nb);
  }
  c0 = Flx_addshift(c0,c,p,n0);
  return Flx_shiftip(av,c0, v);
}


GEN
Flx_mul(GEN x, GEN y, ulong p)
{
 GEN z = Flx_mulspec(x+2,y+2,p, lgpol(x),lgpol(y));
 z[1] = x[1]; return z;
}

static GEN
Flx_sqrspec_basecase(GEN x, ulong p, long nx)
{
  long i, lz, nz;
  ulong p1;
  GEN z;

  if (!nx) return zero_Flx(0);
  lz = (nx << 1) + 1, nz = lz-2;
  z = cgetg(lz, t_VECSMALL) + 2;
  if (u_OK_ULONG(p))
  {
    for (i=0; i<nx; i++)
    {
      p1 = Flx_mullimb_ok(x+i,x,p,0, (i+1)>>1);
      p1 <<= 1; 
      if ((i&1) == 0) p1 += x[i>>1] * x[i>>1];
      z[i] = (p1 % p);
    }
    for (  ; i<nz; i++)
    {
      p1 = Flx_mullimb_ok(x+i,x,p,i-nx+1, (i+1)>>1);
      p1 <<= 1; 
      if ((i&1) == 0) p1 += x[i>>1] * x[i>>1];
      z[i] = (p1 % p);
    }
  }
  else
  {
    for (i=0; i<nx; i++)
    {
      p1 = Flx_mullimb(x+i,x,p,0, (i+1)>>1);
      p1 = Fl_add(p1, p1, p);
      if ((i&1) == 0) p1 = Fl_add(p1, Fl_mul(x[i>>1], x[i>>1], p), p);
      z[i] = p1;
    }
    for (  ; i<nz; i++)
    {
      p1 = Flx_mullimb(x+i,x,p,i-nx+1, (i+1)>>1);
      p1 = Fl_add(p1, p1, p);
      if ((i&1) == 0) p1 = Fl_add(p1, Fl_mul(x[i>>1], x[i>>1], p), p);
      z[i] = p1;
    }
  }
  z -= 2; return z;
}

#if 0
/* used only by Flx_sqrspec #if 0 code.*/
static GEN
Flx_2_mul(GEN x, ulong p)
{
  long i, l = lg(x);
  GEN z = cgetg(l, t_VECSMALL);
  for (i=2; i<l; i++) z[i] = Fl_add(x[i], x[i], p);
  z[1] = x[1]; return z;
}
#endif

static GEN
Flx_sqrspec_sqri(GEN a, ulong p, long na)
{
  GEN z=sqrispec(a,na);
  return Flx_red_lg_i(z,lgefint(z),p);
}

GEN
Flx_sqrspec(GEN a, ulong p, long na)
{
  GEN a0,c,c0,c1;
  long n0, n0a, i, v = 0;
  pari_sp av;

  while (na && !a[0]) { a++; na--; v += 2; }
  av = avma;
  if (na > 30 && maxlenghtcoeffpol(p,na)==1)
    return Flx_shiftip(av, Flx_sqrspec_sqri(a,p,na), v);
  if (na < Flx_SQR_LIMIT) 
    return Flx_shiftip(av, Flx_sqrspec_basecase(a,p,na), v);
  i=(na>>1); n0=na-i; na=i;
  a0=a+n0; n0a=n0;
  while (n0a && !a[n0a-1]) n0a--;

  c = Flx_sqrspec(a,p,n0a);
  c0= Flx_sqrspec(a0,p,na);
  if (p == 2) n0 *= 2;
  else
  {
#if 0
    c1 = Flx_2_mul(Flx_mulspec(a0,a,p, na,n0a), p);
#else
    GEN  t = Flx_addspec(a0,a,p,na,n0a);
    t = Flx_sqr(t,p);
    c1= Flx_add(c0,c, p);
    c1= Flx_sub(t, c1, p);
#endif
    c0 = Flx_addshift(c0,c1,p,n0);
  }
  c0 = Flx_addshift(c0,c,p,n0);
  return Flx_shiftip(av,c0,v);
}

GEN
Flx_sqr(GEN x, ulong p)
{
  GEN z = Flx_sqrspec(x+2,p, lgpol(x));
  z[1] = x[1]; return z;
}

GEN
Flx_pow(GEN x, long n, ulong p)
{
  GEN y = Fl_to_Flx(1,x[1]), z;
  long m;
  if (n == 0) return y;
  m = n; z = x;
  for (;;)
  {
    if (m&1) y = Flx_mul(y,z, p);
    m >>= 1; if (!m) return y;
    z = Flx_sqr(z, p);
  }
}

/* separate from Flx_divrem for maximal speed. */
GEN
Flx_rem(GEN x, GEN y, ulong p)
{ 
  pari_sp av;
  GEN z, c;
  long dx,dy,dz,i,j;
  ulong p1,inv;
  long vs=x[1];

  dy = degpol(y); if (!dy) return zero_Flx(x[1]);
  dx = degpol(x);
  dz = dx-dy; if (dz < 0) return vecsmall_copy(x);
  x += 2; y += 2;
  inv = y[dy];
  if (inv != 1UL) inv = Fl_inv(inv,p);

  c = cgetg(dy+3, t_VECSMALL); c[1]=vs; c += 2; av=avma;
  z = cgetg(dz+3, t_VECSMALL); z[1]=vs; z += 2;

  if (u_OK_ULONG(p))
  {
    z[dz] = (inv*x[dx]) % p;
    for (i=dx-1; i>=dy; --i)
    {
      p1 = p - x[i]; /* compute -p1 instead of p1 (pb with ulongs otherwise) */
      for (j=i-dy+1; j<=i && j<=dz; j++)
      {
        p1 += z[j]*y[i-j];
        if (p1 & HIGHBIT) p1 %= p;
      }
      p1 %= p;
      z[i-dy] = p1? ((p - p1)*inv) % p: 0;
    }
    for (i=0; i<dy; i++)
    {
      p1 = z[0]*y[i];
      for (j=1; j<=i && j<=dz; j++)
      {
        p1 += z[j]*y[i-j];
        if (p1 & HIGHBIT) p1 %= p;
      }
      c[i] = Fl_sub(x[i], p1%p, p);
    }
  }
  else
  {
    z[dz] = Fl_mul(inv, x[dx], p);
    for (i=dx-1; i>=dy; --i)
    {
      p1 = p - x[i]; /* compute -p1 instead of p1 (pb with ulongs otherwise) */
      for (j=i-dy+1; j<=i && j<=dz; j++)
        p1 = Fl_add(p1, Fl_mul(z[j],y[i-j],p), p);
      z[i-dy] = p1? Fl_mul(p - p1, inv, p): 0;
    }
    for (i=0; i<dy; i++)
    {
      p1 = Fl_mul(z[0],y[i],p);
      for (j=1; j<=i && j<=dz; j++)
        p1 = Fl_add(p1, Fl_mul(z[j],y[i-j],p), p);
      c[i] = Fl_sub(x[i], p1, p);
    }
  }
  i = dy-1; while (i>=0 && !c[i]) i--;
  avma=av;
  return Flx_renormalize(c-2, i+3);
}

/* as FpX_divrem but working only on ulong types. ASSUME pr != ONLY_DIVIDES
 * if relevant, *pr is the last object on stack */
GEN
Flx_divrem(GEN x, GEN y, ulong p, GEN *pr)
{
  GEN z,q,c;
  long dx,dy,dz,i,j;
  ulong p1,inv;
  long sv=x[1];

  if (pr == ONLY_REM) return Flx_rem(x, y, p);
  dy = degpol(y);
  if (!dy)
  {
    if (y[2] == 1UL)
      q = vecsmall_copy(x);
    else
      q = Flx_Fl_mul(x, Fl_inv(y[2], p), p);
    if (pr) *pr = zero_Flx(sv);
    return q;
  }
  dx = degpol(x);
  dz = dx-dy;
  if (dz < 0)
  {
    q = zero_Flx(sv);
    if (pr) *pr = vecsmall_copy(x);
    return q;
  }
  x += 2;
  y += 2;
  z = cgetg(dz + 3, t_VECSMALL); z[1] = sv; z += 2;
  inv = (ulong)y[dy];
  if (inv != 1UL) inv = Fl_inv(inv,p);

  if (u_OK_ULONG(p))
  {
    z[dz] = (inv*x[dx]) % p;
    for (i=dx-1; i>=dy; --i)
    {
      p1 = p - x[i]; /* compute -p1 instead of p1 (pb with ulongs otherwise) */
      for (j=i-dy+1; j<=i && j<=dz; j++)
      {
        p1 += z[j]*y[i-j];
        if (p1 & HIGHBIT) p1 %= p;
      }
      p1 %= p;
      z[i-dy] = p1? (long) ((p - p1)*inv) % p: 0;
    }
  }
  else
  {
    z[dz] = Fl_mul(inv, x[dx], p);
    for (i=dx-1; i>=dy; --i)
    { /* compute -p1 instead of p1 (pb with ulongs otherwise) */
      p1 = p - (ulong)x[i];
      for (j=i-dy+1; j<=i && j<=dz; j++)
        p1 = Fl_add(p1, Fl_mul(z[j],y[i-j],p), p);
      z[i-dy] = p1? Fl_mul(p - p1, inv, p): 0;
    }
  }
  q = Flx_renormalize(z-2, dz+3);
  if (!pr) return q;

  c = cgetg(dy + 3, t_VECSMALL); c[1] = sv; c += 2;
  if (u_OK_ULONG(p))
  {
    for (i=0; i<dy; i++)
    {
      p1 = (ulong)z[0]*y[i];
      for (j=1; j<=i && j<=dz; j++)
      {
        p1 += (ulong)z[j]*y[i-j];
        if (p1 & HIGHBIT) p1 %= p;
      }
      c[i] = Fl_sub(x[i], p1%p, p);
    }
  }
  else
  {
    for (i=0; i<dy; i++)
    {
      p1 = Fl_mul(z[0],y[i],p);
      for (j=1; j<=i && j<=dz; j++)
        p1 = Fl_add(p1, Fl_mul(z[j],y[i-j],p), p);
      c[i] = Fl_sub(x[i], p1, p);
    }
  }
  i=dy-1; while (i>=0 && !c[i]) i--;
  c = Flx_renormalize(c-2, i+3);
  *pr = c; return q;
}

long 
Flx_valuation(GEN x)
{
  long i, l=lg(x);
  if (l==2)  return VERYBIGINT;
  for (i=2; i<l && x[i]==0; i++);
  return i-2;
}

GEN
Flx_recipspec(GEN x, long l, long n)
{
  long i;
  GEN z=cgetg(n+2,t_VECSMALL)+2;
  for(i=0; i<l; i++)
    z[n-i-1] = x[i];
  for(   ; i<n; i++)
    z[n-i-1] = 0;
  return Flx_renormalize(z-2,n+2);
}

GEN
Flx_recip(GEN x)
{
  GEN z=Flx_recipspec(x+2,lgpol(x),lgpol(x));
  z[1]=x[1];
  return z;
}

/*
 * x/polrecip(P)+O(x^n)
 */
static GEN 
Flx_invmontgomery_basecase(GEN T, ulong p)
{
  long i, l=lg(T)-1, k;
  GEN r=cgetg(l,t_VECSMALL); r[1]=T[1];
  r[2]=0; r[3]=1;
  if (u_OK_ULONG(p)) {
    for (i=4;i<l;i++)
    {
      long u = 0;
      for (k=3;k<i;k++) { u += T[l-i+k] * r[k]; if (u & HIGHBIT) u %= p; }
      u %= p;
      r[i] = Fl_neg(u, p);
    }
  }
  else {
    for (i=4;i<l;i++)
    {
      ulong u = 0;
      for (k=3;k<i;k++) u = Fl_sub(u, Fl_mul(T[l-i+k],r[k],p),p);
      r[i] = u;
    }
  }
  r = Flx_renormalize(r,l);
  return r;
}

/* Return new lgpol */
static long
Flx_lgrenormalizespec(GEN x, long lx)
{
  long i;
  for (i = lx-1; i>=0; i--)
    if (x[i]) break;
  return i+1;
}
static GEN
Flx_invmontgomery_newton(GEN T, ulong p)
{
  long i, lx, lz, lq, e, l = degpol(T);
  GEN E, q, y, z, x = const_vecsmall(l+1, 0) + 2;
  pari_sp av;

  y = T+2;
  q = Flx_recipspec(y,l,l+1) + 2;
  E = Newton_exponents(l-2); /* assume l > 3 */
  av = avma;
  /* We work on _spec_ Flx's, all the l[xzq12] below are lgpol's */

  q[0] = y[l];
  /* initialize */
  x[0] = Fl_inv(q[0], p);
  if (q[1])
  {
    ulong u = q[1];
    if (x[0] != 1) u = Fl_mul(u, Fl_sqr(x[0],p), p);
    x[1] = p - u; lx = 2;
  }
  else
    lx = 1;
  for (e = lg(E)-1; e > 1; e--, avma = av)
  { /* set x -= x(x*q - 1) + O(t^l2), knowing x*q = 1 + O(t^l1) */
    long l2 = E[e-1] + 1, l1 = E[e] + 1;

    lq = Flx_lgrenormalizespec(q, l2);
    z = Flx_mulspec(x, q, p, lx, lq); /* FIXME: high product */
    lz = lgpol(z); if (lz > l2) lz = l2;
    z += 2;
    /* subtract 1 [=>first l1-1 words are 0]: renormalize so that z(0) != 0 */
    for (i = l1-1; i < lz; i++) if (z[i]) break;
    if (i >= lz) continue; /* z-1 = 0(t^l2) */

    /* z + i represents (x*q - 1) / t^i */
    z = Flx_mulspec(x, z+i, p, lx, lz-i); /* FIXME: low product */
    lz = lgpol(z); z += 2;
    if (lz > l2-i) lz = Flx_lgrenormalizespec(z, l2-i);

    lx = lz+ i;
    y  = x + i; /* x -= z * t^i, in place */
    for (i = 0; i < lz; i++) y[i] = Fl_neg(z[i], p);
  }
  x -= 2; setlg(x, lx + 2); x[1] = T[1];
  return Flx_shift(x,1);
}

/* x/polrecip(T)+O(x^deg(T)) */
GEN
Flx_invmontgomery(GEN T, ulong p)
{
  pari_sp ltop=avma;
  long l=lg(T);
  GEN r;
  if (l<5) return zero_Flx(T[1]);
  if (l<Flx_INVMONTGOMERY_LIMIT)
  {
    ulong c=T[l-1], ci=1;
    if (c!=1)
    {
      ci=Fl_inv(c,p);
      T=Flx_Fl_mul(T, ci, p);
    }
    r=Flx_invmontgomery_basecase(T,p);
    if (c!=1) r=Flx_Fl_mul(r,ci,p);
  }
  else
    r=Flx_invmontgomery_newton(T,p);
  return gerepileuptoleaf(ltop, r);
}

/* Compute x mod T where lg(x)<=2*lg(T)-2
 * and mg is the Montgomery inverse of T. 
 */
GEN
Flx_rem_montgomery(GEN x, GEN mg, GEN T, ulong p)
{
  pari_sp ltop=avma;
  GEN z;
  long l=lgpol(x);
  long lt=degpol(T); /*We discard the leading term*/
  long lead=lt-1;
  long ld=l-lt+1;
  long lm=min(ld,lgpol(mg));
  if (l<=lt)
    return vecsmall_copy(x);
  (void)new_chunk(lt);
  z=Flx_recipspec(x+2+lead,ld,ld);             /* z = rec(x)      lz<=ld*/
  z=Flx_mulspec(z+2,mg+2,p,lgpol(z),lm);       /* z = rec(x) * mg lz<=ld+lm*/
  z=Flx_recipspec(z+2,min(ld,lgpol(z)),ld);    /* z = rec (rec(x) * mg) lz<=ld*/
  z=Flx_mulspec(z+2,T+2,p,lgpol(z),lt);        /* z *= pol        lz<=ld+lt*/
  avma=ltop;
  z=Flx_subspec(x+2,z+2,p,lt,min(lt,lgpol(z)));/* z = x - z       lz<=lt */
  z[1]=T[1];
  return z;
}

GEN
Flx_deriv(GEN z, ulong p)
{
  long i,l = lg(z)-1;
  GEN x;
  if (l < 2) l = 2;
  x = cgetg(l, t_VECSMALL); x[1] = z[1]; z++;
  if (HIGHWORD(l | p))
    for (i=2; i<l; i++) x[i] = Fl_mul((ulong)i-1, z[i], p);
  else
    for (i=2; i<l; i++) x[i] = ((i-1) * z[i]) % p;
  return Flx_renormalize(x,l);
}

/*Do not garbage collect*/
GEN
Flx_gcd_i(GEN a, GEN b, ulong p)
{
  GEN c;
  if (lg(b) > lg(a)) swap(a, b);
  while (lgpol(b))
  {
    c = Flx_rem(a,b,p);
    a = b; b = c;
  }
  return a;
}

GEN
Flx_gcd(GEN a, GEN b, ulong p)
{
  pari_sp av = avma;
  return gerepileuptoleaf(av, Flx_gcd_i(a,b,p));
}

int
Flx_is_squarefree(GEN z, ulong p)
{
  pari_sp av = avma;
  GEN d = Flx_gcd_i(z, Flx_deriv(z,p) , p);
  long res= (degpol(d) == 0);
  avma = av; return res;
}

GEN
Flx_extgcd(GEN a, GEN b, ulong p, GEN *ptu, GEN *ptv)
{
  GEN q,z,u,v, x = a, y = b;

  u = zero_Flx(a[1]);
  v = Fl_to_Flx(1,a[1]); /* v = 1 */
  while (lgpol(y))
  {
    q = Flx_divrem(x,y,p,&z);
    x = y; y = z; /* (x,y) = (y, x - q y) */
    z = Flx_sub(u, Flx_mul(q,v, p), p);
    u = v; v = z; /* (u,v) = (v, u - q v) */
  }
  z = Flx_sub(x, Flx_mul(b,u,p), p);
  z = Flx_div(z,a,p);
  *ptu = z;
  *ptv = u; return x;
}

ulong
Flx_resultant(GEN a, GEN b, ulong p)
{
  long da,db,dc,cnt;
  ulong lb, res = 1UL;
  pari_sp av;
  GEN c;

  if (lgpol(a)==0 || lgpol(b)==0) return 0;
  da = degpol(a);
  db = degpol(b);
  if (db > da)
  {
    swapspec(a,b, da,db);
    if (both_odd(da,db)) res = p-res;
  }
  if (!da) return 1; /* = res * a[2] ^ db, since 0 <= db <= da = 0 */
  cnt = 0; av = avma;
  while (db)
  {
    lb = b[db+2];
    c = Flx_rem(a,b, p);
    a = b; b = c; dc = degpol(c);
    if (dc < 0) { avma = av; return 0; }

    if (both_odd(da,db)) res = p - res;
    if (lb != 1) res = Fl_mul(res, Fl_pow(lb, da - dc, p), p);
    if (++cnt == 4) { cnt = 0; avma = av; }
    da = db; /* = degpol(a) */
    db = dc; /* = degpol(b) */
  }
  avma = av; return Fl_mul(res, Fl_pow(b[2], da, p), p);
}

/* If resultant is 0, *ptU and *ptU are not set */
ulong
Flx_extresultant(GEN a, GEN b, ulong p, GEN *ptU, GEN *ptV)
{
  GEN z,q,u,v, x = a, y = b;
  ulong lb, res = 1UL;
  pari_sp av = avma;
  long dx, dy, dz;
  long vs=a[1];

  dx = degpol(x);
  dy = degpol(y);
  if (dy > dx)
  {
    swap(x,y); lswap(dx,dy); pswap(ptU, ptV);
    a = x; b = y;
    if (both_odd(dx,dy)) res = p-res;
  }
  /* dx <= dy */
  if (dx < 0) return 0;

  u = zero_Flx(vs);
  v = Fl_to_Flx(1,vs); /* v = 1 */
  while (dy)
  { /* b u = x (a), b v = y (a) */
    lb = y[dy+2];
    q = Flx_divrem(x,y, p, &z);
    x = y; y = z; /* (x,y) = (y, x - q y) */
    dz = degpol(z); if (dz < 0) { avma = av; return 0; }
    z = Flx_sub(u, Flx_mul(q,v, p), p);
    u = v; v = z; /* (u,v) = (v, u - q v) */

    if (both_odd(dx,dy)) res = p - res;
    if (lb != 1) res = Fl_mul(res, Fl_pow(lb, dx-dz, p), p);
    dx = dy; /* = degpol(x) */
    dy = dz; /* = degpol(y) */
  }
  res = Fl_mul(res, Fl_pow(y[2], dx, p), p);
  lb = Fl_mul(res, Fl_inv(y[2],p), p);
  v = gerepileuptoleaf(av, Flx_Fl_mul(v, lb, p));
  av = avma;
  u = Flx_sub(Fl_to_Flx(res,vs), Flx_mul(b,v,p), p);
  u = gerepileuptoleaf(av, Flx_div(u,a,p)); /* = (res - b v) / a */
  *ptU = u;
  *ptV = v; return res;
}

ulong
Flx_eval(GEN x, ulong y, ulong p)
{
  ulong p1,r;
  long j, i=lg(x)-1;
  if (i<=2)
    return (i==2)? x[2]: 0;
  p1 = x[i];
  /* specific attention to sparse polynomials (see poleval)*/
  if (u_OK_ULONG(p))
  {
    for (i--; i>=2; i=j-1)
    {
      for (j=i; !x[j]; j--)
        if (j==2)
        {
          if (i != j) y = Fl_pow(y, i-j+1, p);
          return (p1 * y) % p;
        }
      r = (i==j)? y: Fl_pow(y, i-j+1, p);
      p1 = ((p1*r) + x[j]) % p;
    }
  }
  else
  {
    for (i--; i>=2; i=j-1)
    {
      for (j=i; !x[j]; j--)
        if (j==2)
        {
          if (i != j) y = Fl_pow(y, i-j+1, p);
          return Fl_mul(p1, y, p);
        }
      r = (i==j)? y: Fl_pow(y, i-j+1, p);
      p1 = Fl_add((ulong)x[j], Fl_mul(p1,r,p), p);
    }
  }
  return p1;
}

static GEN
_Flx_mul(void *p, GEN a, GEN b)
{
  return Flx_mul(a,b, *(ulong*)p);
}

/* compute prod (x - a[i]) */
GEN
Flv_roots_to_pol(GEN a, ulong p, long vs)
{
  long i,k,lx = lg(a);
  GEN p1,p2;
  if (lx == 1) return Fl_to_Flx(1,vs);
  p1 = cgetg(lx, t_VEC); 
  for (k=1,i=1; i<lx-1; i+=2)
  {
    p2 = cgetg(5,t_VECSMALL); gel(p1,k++) = p2;
    p2[1] = vs;
    p2[2] = Fl_mul(a[i], a[i+1], p);
    p2[3] = a[i] + a[i+1];
    if ((ulong)p2[3] >= p) p2[3] -= p;
    if (p2[3]) p2[3] = p - p2[3]; /* - (a[i] + a[i+1]) mod p */
    p2[4] = 1; 
  }
  if (i < lx)
  {
    p2 = cgetg(4,t_VECSMALL); gel(p1,k++) = p2;
    p2[1] = vs;
    p2[2] = a[i]?p - a[i]:0;
    p2[3] = 1;
  }
  setlg(p1, k); return divide_conquer_assoc(p1, _Flx_mul,(void *)&p);
}

GEN
Flx_div_by_X_x(GEN a, ulong x, ulong p, ulong *rem)
{
  long l = lg(a), i;
  GEN a0, z0;
  GEN z = cgetg(l-1,t_VECSMALL);
  z[1] = a[1];
  a0 = a + l-1;
  z0 = z + l-2; *z0 = *a0--;
  if (u_OK_ULONG(p))
  {
    for (i=l-3; i>1; i--) /* z[i] = (a[i+1] + x*z[i+1]) % p */
    {
      ulong t = (*a0-- + x *  *z0--) % p;
      *z0 = (long)t;
    }
    if (rem) *rem = (*a0 + x *  *z0) % p;
  }
  else
  {
    for (i=l-3; i>1; i--)
    {
      ulong t = Fl_add((ulong)*a0--, Fl_mul(x, *z0--, p), p);
      *z0 = (long)t;
    }
    if (rem) *rem = Fl_add((ulong)*a0, Fl_mul(x, *z0, p), p);
  }
  return z;
}

/* u P(X) + v P(-X) */
GEN
Flx_even_odd_comb(GEN P, ulong u, ulong v, ulong p)
{
  long i, l = lg(P);
  GEN y = cgetg(l,t_VECSMALL);
  y[1]=P[1];
  for (i=2; i<l; i++)
  {
    ulong t = P[i];
    y[i] = (t == 0)? 0:
                     (i&1)? Fl_mul(t, u + (p - v), p)
                          : Fl_mul(t, u + v, p);
  }
  return Flx_renormalize(y,l);
}

/* xa, ya = t_VECSMALL */
GEN
Flv_polint(GEN xa, GEN ya, ulong p, long vs)
{
  long i, j, n = lg(xa);
  GEN T,dP, P = cgetg(n+1, t_VECSMALL);
  GEN Q = Flv_roots_to_pol(xa, p, vs);
  ulong inv;
  P[1] = vs;
  for (j=2; j<=n; j++) P[j] = 0UL;
  for (i=1; i<n; i++)
  {
    if (!ya[i]) continue;
    T = Flx_div_by_X_x(Q, xa[i], p, NULL);
    inv = Fl_inv(Flx_eval(T,xa[i], p), p);
    if (i < n-1 && (ulong)(xa[i] + xa[i+1]) == p)
    {
      dP = Flx_even_odd_comb(T, Fl_mul(ya[i],inv,p), Fl_mul(ya[i+1],inv,p), p);
      i++; /* x_i = -x_{i+1} */
    }
    else
      dP = Flx_Fl_mul(T, Fl_mul(ya[i],inv,p), p);
    for (j=2; j<lg(dP); j++) P[j] = Fl_add(P[j], dP[j], p);
    avma = (pari_sp)Q;
  }
  avma = (pari_sp)P;
  return Flx_renormalize(P,n+1);
}

/***********************************************************************/
/**                                                                   **/
/**               Flxq                                                **/
/**                                                                   **/
/***********************************************************************/
/* Flxq objects are defined as follows:
   They are Flx modulo another Flx called q.
*/

/* Product of y and x in Z/pZ[X]/(pol), as t_VECSMALL. */
GEN
Flxq_mul(GEN y,GEN x,GEN pol,ulong p)
{
  GEN z = Flx_mul(y,x,p);
  return Flx_rem(z,pol,p);
}

/* Square of y in Z/pZ[X]/(pol), as t_VECSMALL. */
GEN
Flxq_sqr(GEN y,GEN pol,ulong p)
{
  GEN z = Flx_sqr(y,p);
  return Flx_rem(z,pol,p);
}

typedef struct {
  GEN pol;
  GEN mg;
  ulong p;
} Flxq_muldata;


static GEN
_sqr_montgomery(void *data, GEN x)
{
  Flxq_muldata *D = (Flxq_muldata*)data;
  return Flx_rem_montgomery(Flx_sqr(x,D->p),D->mg, D->pol, D->p);
}
static GEN
_mul_montgomery(void *data, GEN x, GEN y)
{
  Flxq_muldata *D = (Flxq_muldata*)data;
  return Flx_rem_montgomery(Flx_mul(x,y,D->p),D->mg, D->pol, D->p);
}

static GEN
_Flxq_sqr(void *data, GEN x)
{
  Flxq_muldata *D = (Flxq_muldata*)data;
  return Flxq_sqr(x, D->pol, D->p);
}
static GEN
_Flxq_mul(void *data, GEN x, GEN y)
{
  Flxq_muldata *D = (Flxq_muldata*)data;
  return Flxq_mul(x,y, D->pol, D->p);
}

/* n-Power of x in Z/pZ[X]/(pol), as t_VECSMALL. */
GEN
Flxq_pow(GEN x, GEN n, GEN pol, ulong p)
{
  pari_sp av = avma;
  Flxq_muldata D;
  GEN y;
  if (!signe(n)) return Fl_to_Flx(1,pol[1]);
  if (signe(n) < 0)
    x=Flxq_inv(x,pol,p);
  else
    x=Flx_rem(x, pol, p);
  if (is_pm1(n)) return x;
  D.pol = pol;
  D.p   = p;
  /* not tuned*/
  if (pol[2] && degpol(pol) >= Flx_POW_MONTGOMERY_LIMIT)
  {
    /* We do not handle polynomials multiple of x yet */
    D.mg  = Flx_invmontgomery(pol,p);
    y = leftright_pow(x, n, (void*)&D, &_sqr_montgomery, &_mul_montgomery);
  }
  else
    y = leftright_pow(x, n, (void*)&D, &_Flxq_sqr, &_Flxq_mul);
  return gerepileuptoleaf(av, y);
}

/* Inverse of x in Z/pZ[X]/(pol) or NULL if inverse doesn't exist
 * return lift(1 / (x mod (p,pol)))
 * not stack clean.
 * */
GEN
Flxq_invsafe(GEN x, GEN T, ulong p)
{
  GEN U, V;
  GEN z = Flx_extgcd(x, T, p, &U, &V);
  ulong iz;
  if (degpol(z)) return NULL;
  iz = Fl_inv ((ulong)z[2], p);
  return Flx_Fl_mul(U, iz, p);
}

GEN
Flxq_inv(GEN x,GEN T,ulong p)
{
  pari_sp av=avma;
  GEN U = Flxq_invsafe(x, T, p);
  if (!U) pari_err(talker,"non invertible polynomial in Flxq_inv");
  return gerepileuptoleaf(av, U);
}

/* generates the list of powers of x of degree 0,1,2,...,l*/
GEN
Flxq_powers(GEN x, long l, GEN T, ulong p)
{
  GEN V = cgetg(l+2,t_VEC);
  long i, v = T[1];
  gel(V,1) = Fl_to_Flx(1, v);  if (l==0) return V;
  gel(V,2) = vecsmall_copy(x); if (l==1) return V;
  gel(V,3) = Flxq_sqr(x,T,p);
  if ((degpol(x)<<1) < degpol(T)) {
    for(i = 4; i < l+2; i++)
      gel(V,i) = Flxq_mul(gel(V,i-1),x,T,p);
  } else {
    for(i = 4; i < l+2; i++) {
      gel(V,i) = (i&1)? Flxq_sqr(gel(V, (i+1)>>1),T,p)
                      : Flxq_mul(gel(V, i-1),x,T,p);
    }
  }
  return V;
}

GEN
FlxV_Flc_mul(GEN V, GEN W, ulong p)
{
  pari_sp ltop=avma;
  long i;
  GEN z = Flx_Fl_mul(gel(V,1),W[1],p);
  for(i=2;i<lg(V);i++)
    z=Flx_add(z,Flx_Fl_mul(gel(V,i),W[i],p),p);
  return gerepileuptoleaf(ltop,z);
}

GEN
ZXV_to_FlxV(GEN v, ulong p)
{
  long j, N = lg(v);
  GEN y = cgetg(N, t_VEC);
  for (j=1; j<N; j++) gel(y,j) = ZX_to_Flx(gel(v,j), p);
  return y;
}

GEN
FlxV_to_Flm(GEN v, long n)
{
  long j, N = lg(v);
  GEN y = cgetg(N, t_MAT);
  for (j=1; j<N; j++) gel(y,j) = Flx_to_Flv(gel(v,j), n);
  return y;
}


/**************************************************************
 **                 FlxX                                    **
 **                                                          **
 **************************************************************/
/*Similar to normalizepol, in place*/
/*FlxX_renormalize=zxX_renormalize */
GEN
FlxX_renormalize(GEN /*in place*/ x, long lx)
{
  long i;
  for (i = lx-1; i>1; i--)
    if (lgpol(x[i])) break;
  stackdummy((pari_sp)(x + lg(x)), (pari_sp)(x + i+1));
  setlg(x, i+1); setsigne(x, i!=1); return x;
}

/* FlxX are t_POL with Flx coefficients.
 * Normally the variable ordering should be respected.*/
GEN 
FlxX_to_ZXX(GEN B)
{
  long lb=lg(B);
  long i;
  GEN b=cgetg(lb,t_POL);
  for (i=2; i<lb; i++) 
    if (lgpol(B[i]))
      gel(b,i) = Flx_to_ZX(gel(B,i));
    else
      gel(b,i) = gen_0;
  b[1] = B[1]; return b;
}

/* Note: v is used _only_ for the t_INT. It must match
 * the variable of any t_POL coefficients. */
GEN 
ZXX_to_FlxX(GEN B, ulong p, long v)
{
  long lb=lg(B);
  long i;
  GEN b=cgetg(lb,t_POL);
  b[1]=evalsigne(1)|(((ulong)B[1])&VARNBITS);
  for (i=2; i<lb; i++) 
    switch (typ(B[i]))
    {
    case t_INT:  
      gel(b,i) = Z_to_Flx(gel(B,i), p, v);
      break;
    case t_POL:
      gel(b,i) = ZX_to_Flx(gel(B,i), p);
      break;
    }
  return FlxX_renormalize(b, lb);
}

GEN
ZXXV_to_FlxXV(GEN V, ulong p, long v)
{
  long j, N = lg(V);
  GEN y = cgetg(N, t_VEC);
  for (j=1; j<N; j++) gel(y,j) = ZXX_to_FlxX(gel(V,j), p, v);
  return y;
}

/* matrix whose entries are given by the coeffs of the polynomial v in
 * two variables (considered as degree n polynomials) */
GEN
FlxX_to_Flm(GEN v, long n)
{
  long j, N = lg(v)-1;
  GEN y = cgetg(N, t_MAT);
  v++;
  for (j=1; j<N; j++) gel(y,j) = Flx_to_Flv(gel(v,j), n);
  return y;
}

GEN
Flm_to_FlxX(GEN x, long v,long w)
{
  long j, lx = lg(x);
  GEN y = cgetg(lx+1, t_POL);
  y[1]=evalsigne(1) | v;
  y++;
  for (j=1; j<lx; j++) gel(y,j) = Flv_to_Flx(gel(x,j), w);
  return FlxX_renormalize(--y, lx+1);
}

/* P(X,Y) --> P(Y,X), n-1 is the degree in Y */
GEN
FlxX_swap(GEN x, long n, long ws)
{
  long j, lx = lg(x), ly = n+3;
  GEN y = cgetg(ly, t_POL);
  y[1] = x[1];
  for (j=2; j<ly; j++)
  {
    long k;
    GEN p1 = cgetg(lx, t_VECSMALL);
    p1[1] = ws;
    for (k=2; k<lx; k++)
      if( j<lg(x[k]))
        p1[k] = mael(x,k,j);
      else
        p1[k] = 0;
    gel(y,j) = Flx_renormalize(p1,lx);
  }
  return FlxX_renormalize(y,ly);
}

/*Fix me should be  zxX_to_Kronecker since it does not use l*/
GEN
FlxX_to_Kronecker_spec(GEN P, long lp, GEN Q)
{
  /* P(X) = sum Mod(.,Q(Y)) * X^i, lift then set X := Y^(2n-1) */
  long i,j,k,l;
  long N = (degpol(Q)<<1) + 1;
  GEN p1;
  GEN y = cgetg((N-2)*lp + 2, t_VECSMALL)+2;
  for (k=i=0; i<lp; i++)
  {
    p1 = gel(P,i);
    l = lg(p1);
    for (j=2; j < l; j++) y[k++] = p1[j];
    if (i == lp-1) break;
    for (   ; j < N; j++) y[k++] = 0;
  }
  y-=2;
  setlg(y, k+2); return y;
}

/*Fix me should be  zxX_to_Kronecker since it does not use l*/
GEN
FlxX_to_Kronecker(GEN P, GEN Q)
{
  /* P(X) = sum Mod(.,Q(Y)) * X^i, lift then set X := Y^(2n-1) */
  long i,j,k,l;
  long lx = lg(P), N = (degpol(Q)<<1) + 1;
  GEN p1;
  GEN y = cgetg((N-2)*(lx-2) + 2, t_VECSMALL);
  y[1] = P[1];
  for (k=i=2; i<lx; i++)
  {
    p1 = gel(P,i);
    l = lg(p1);
    for (j=2; j < l; j++) y[k++] = p1[j];
    if (i == lx-1) break;
    for (   ; j < N; j++) y[k++] = 0;
  }
  setlg(y, k); return y;
}

GEN
FlxX_add(GEN x, GEN y, ulong p)
{
  long i,lz;
  GEN z; 
  long lx=lg(x);
  long ly=lg(y);
  if (ly>lx) swapspec(x,y, lx,ly);
  lz = lx; z = cgetg(lz, t_POL); z[1]=x[1];
  for (i=2; i<ly; i++) gel(z,i) = Flx_add(gel(x,i), gel(y,i), p);
  for (   ; i<lx; i++) gel(z,i) = vecsmall_copy(gel(x,i));
  return FlxX_renormalize(z, lz);
}

GEN
FlxX_subspec(GEN x, GEN y, ulong p, long lx, long ly)
{
  long i,lz;
  GEN z;

  if (ly <= lx)
  {
    lz = lx+2; z = cgetg(lz, t_POL)+2;
    for (i=0; i<ly; i++) gel(z,i) = Flx_sub(gel(x,i),gel(y,i),p);
    for (   ; i<lx; i++) gel(z,i) = vecsmall_copy(gel(x,i));
  }
  else
  {
    lz = ly+2; z = cgetg(lz, t_POL)+2;
    for (i=0; i<lx; i++) gel(z,i) = Flx_sub(gel(x,i),gel(y,i),p);
    for (   ; i<ly; i++) gel(z,i) = Flx_neg(gel(x,i),p);
  }
 return FlxX_renormalize(z-2, lz);
}


/*Unused/untested*/
GEN
FlxX_sub(GEN x, GEN y, ulong p)
{
  long lx,ly,i,lz;
  GEN z;
  lx = lg(x); ly = lg(y);
  lz=max(lx,ly);
  z = cgetg(lz,t_POL);
  if (lx >= ly)
  {
    z[1] = x[1];
    for (i=2; i<ly; i++) gel(z,i) = Flx_sub(gel(x,i),gel(y,i),p);
    for (   ; i<lx; i++) gel(z,i) = vecsmall_copy(gel(x,i));
    if (lx==ly) z = FlxX_renormalize(z, lz);
  }
  else
  {
    z[1] = y[1];
    for (i=2; i<lx; i++) gel(z,i) = Flx_sub(gel(x,i),gel(y,i),p);
    for (   ; i<ly; i++) gel(z,i) = Flx_neg(gel(y,i),p);
  }
  if (!lgpol(z)) { avma = (pari_sp)(z + lz); z = zeropol(varn(x)); }
  return z;
}

GEN
FlxX_shift(GEN a, long n)
{
  long i, l = lg(a);
  GEN  b;
  long vs;
  if (!signe(a)) return a;
  vs = mael(a,2,1);
  b = cgetg(l+n, t_POL);
  b[1] = a[1];
  for (i=0; i<n; i++) gel(b,2+i) = zero_Flx(vs);
  for (i=2; i<l; i++) b[i+n] = a[i];
  return b;
}

GEN
FlxX_recipspec(GEN x, long l, long n, long vs)
{
  long i;
  GEN z=cgetg(n+2,t_POL)+2;
  for(i=0; i<l; i++)
    gel(z,n-i-1) = vecsmall_copy(gel(x,i));
  for(   ; i<n; i++)
    gel(z,n-i-1) = zero_Flx(vs);
  return FlxX_renormalize(z-2,n+2);
}

/**************************************************************
 **                 FlxqX                                    **
 **                                                          **
 **************************************************************/

/* FlxqX are t_POL with Flxq coefficients.
 * Normally the variable ordering should be respected.*/

/*Not stack clean.*/
GEN
FlxqX_from_Kronecker(GEN z, GEN T, ulong p)
{
  long i,j,lx,l, N = (degpol(T)<<1) + 1;
  GEN x, t = cgetg(N,t_VECSMALL);
  t[1] = T[1];
  l = lg(z); lx = (l-2) / (N-2);
  x = cgetg(lx+3,t_POL); x[1]=z[1];
  for (i=2; i<lx+2; i++)
  {
    for (j=2; j<N; j++) t[j] = z[j];
    z += (N-2);
    gel(x,i) = Flx_rem(Flx_renormalize(t,N), T,p);
  }
  N = (l-2) % (N-2) + 2;
  for (j=2; j<N; j++) t[j] = z[j];
  gel(x,i) = Flx_rem(Flx_renormalize(t,N), T,p);
  return FlxX_renormalize(x, i+1);
}

GEN
FlxqX_red(GEN z, GEN T, ulong p)
{
  GEN res;
  long i, l = lg(z);
  res = cgetg(l,t_POL); res[1] = z[1];
  for(i=2;i<l;i++)
    gel(res,i) = Flx_rem(gel(z,i),T,p);
  return FlxX_renormalize(res,lg(res));
}

GEN
FlxqX_mulspec(GEN x, GEN y, GEN T, ulong p, long lx, long ly)
{
  pari_sp ltop=avma;
  GEN z,kx,ky;
  kx= FlxX_to_Kronecker_spec(x,lx,T);
  ky= FlxX_to_Kronecker_spec(y,ly,T);
  z = Flx_mul(ky, kx, p);
  z = FlxqX_from_Kronecker(z,T,p);
  return gerepileupto(ltop,z);
}

GEN
FlxqX_mul(GEN x, GEN y, GEN T, ulong p)
{
  pari_sp ltop=avma;
  GEN z,kx,ky;
  kx= FlxX_to_Kronecker(x,T);
  ky= FlxX_to_Kronecker(y,T);
  z = Flx_mul(ky, kx, p);
  z = FlxqX_from_Kronecker(z,T,p);
  return gerepileupto(ltop,z);
}

GEN
FlxqX_sqr(GEN x, GEN T, ulong p)
{
  pari_sp ltop=avma;
  GEN z,kx;
  kx= FlxX_to_Kronecker(x,T);
  z = Flx_sqr(kx, p); 
  z = FlxqX_from_Kronecker(z,T,p);
  return gerepileupto(ltop,z);
}

GEN
FlxqX_Flxq_mul(GEN P, GEN U, GEN T, ulong p)
{
  long i, lP = lg(P);
  GEN res = cgetg(lP,t_POL);
  res[1] = P[1];
  for(i=2; i<lP; i++)
    gel(res,i) = Flxq_mul(U,gel(P,i), T,p);
  return FlxX_renormalize(res,lg(res));
}

GEN
FlxqX_normalize(GEN z, GEN T, ulong p)
{
  GEN p1 = leading_term(z);
  if (!lgpol(z) || (!degpol(p1) && p1[1] == 1)) return z;
  return FlxqX_Flxq_mul(z, Flxq_inv(p1,T,p), T,p);
}

/* x and y in Z[Y][X]. Assume T irreducible mod p */
GEN
FlxqX_divrem(GEN x, GEN y, GEN T, ulong p, GEN *pr)
{
  long vx, dx, dy, dz, i, j, sx, lr;
  pari_sp av0, av, tetpil;
  GEN z,p1,rem,lead;

  if (!signe(y)) pari_err(gdiver);
  vx=varn(x); dy=degpol(y); dx=degpol(x);
  if (dx < dy)
  {
    if (pr)
    {
      av0 = avma; x = FlxqX_red(x, T, p);
      if (pr == ONLY_DIVIDES) 
      { 
        avma=av0; 
        return signe(x)? NULL: zeropol(vx); 
      }
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
    av0 = avma; x = FlxqX_normalize(x,T,p); tetpil = avma;
    return gerepile(av0,tetpil,FlxqX_red(x,T,p));
  }
  av0 = avma; dz = dx-dy;
  lead = (!degpol(lead) && lead[2]==1)? NULL: gclone(Flxq_inv(lead,T,p));
  avma = av0;
  z = cgetg(dz+3,t_POL); z[1] = x[1];
  x += 2; y += 2; z += 2;

  p1 = gel(x,dx); av = avma;
  gel(z,dz) = lead? gerepileupto(av, Flxq_mul(p1,lead, T, p)): gcopy(p1);
  for (i=dx-1; i>=dy; i--)
  {
    av=avma; p1=gel(x,i);
    for (j=i-dy+1; j<=i && j<=dz; j++)
      p1 = Flx_sub(p1, Flx_mul(gel(z,j),gel(y,i-j),p),p);
    if (lead) p1 = Flx_mul(p1, lead,p);
    tetpil=avma; gel(z,i-dy) = gerepile(av,tetpil,Flx_rem(p1,T,p));
  }
  if (!pr) { if (lead) gunclone(lead); return z-2; }

  rem = (GEN)avma; av = (pari_sp)new_chunk(dx+3);
  for (sx=0; ; i--)
  {
    p1 = gel(x,i);
    for (j=0; j<=i && j<=dz; j++)
      p1 = Flx_sub(p1, Flx_mul(gel(z,j),gel(y,i-j),p),p);
    tetpil=avma; p1 = Flx_rem(p1, T, p); if (lgpol(p1)) { sx = 1; break; }
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
      p1 = Flx_sub(p1, Flx_mul(gel(z,j),gel(y,i-j),p), p);
    tetpil=avma; gel(rem,i) = gerepile(av,tetpil, Flx_rem(p1, T, p));
  }
  rem -= 2;
  if (lead) gunclone(lead);
  if (!sx) (void)FlxX_renormalize(rem, lr);
  if (pr == ONLY_REM) return gerepileupto(av0,rem);
  *pr = rem; return z-2;
}

static GEN 
FlxqX_invmontgomery_basecase(GEN T, GEN Q, ulong p)
{
  long i, l=lg(T)-1, k;
  long sv=Q[1];
  GEN r=cgetg(l,t_POL); r[1]=T[1];
  gel(r,2) = zero_Flx(sv); 
  gel(r,3) = Fl_to_Flx(1,sv);
  for (i=4;i<l;i++)
  {
    pari_sp ltop=avma;
    GEN z = zero_Flx(sv);
    for (k=3;k<i;k++)
      z = Flx_sub(z,Flxq_mul(gel(T,l-i+k),gel(r,k),Q,p),p);
    gel(r,i) = gerepileupto(ltop, z);
  }
  r = FlxX_renormalize(r,l);
  return r;
}

/* x/polrecip(P)+O(x^n) */
GEN
FlxqX_invmontgomery(GEN T, GEN Q, ulong p)
{
  pari_sp ltop=avma;
  long l=lg(T);
  GEN r;
  GEN c=gel(T,l-1), ci=NULL;
  if (l<5) return zero_Flx(T[1]);
  if (degpol(c) || c[2]!=1)
  {
    ci=Flxq_inv(c,Q,p);
    T=FlxqX_Flxq_mul(T, ci, Q, p);
  }
  r=FlxqX_invmontgomery_basecase(T,Q,p);
  if (ci) r=FlxqX_Flxq_mul(r,ci,Q,p);
  return gerepileupto(ltop, r);
}

GEN
FlxqX_rem_montgomery(GEN x, GEN mg, GEN T, GEN Q, ulong p)
{
  pari_sp ltop=avma;
  GEN z;
  long vs=Q[1];
  long l=lgpol(x);
  long lt=degpol(T); /*We discard the leading term*/
  long lead=lt-1;
  long ld=l-lt+1;
  long lm=min(ld,lgpol(mg));
  if (l<=lt)
    return gcopy(x);
  z=FlxX_recipspec(x+2+lead,ld,ld,vs);         /* z = rec(x)      lz<=ld*/
  z=FlxqX_mulspec(z+2,mg+2,Q,p,lgpol(z),lm);   /* z = rec(x) * mg lz<=ld+lm*/
  z=FlxX_recipspec(z+2,min(ld,lgpol(z)),ld,vs);/* z = rec (rec(x) * mg) lz<=ld*/
  z=FlxqX_mulspec(z+2,T+2,Q,p,lgpol(z),lt);    /* z*= pol         lz<=ld+lt*/
  z=FlxX_subspec(x+2,z+2,p,lt,min(lt,lgpol(z)));/*z = x - z       lz<=lt */
  z[1]=T[1];
  return gerepileupto(ltop,z);
}

GEN
FlxqX_safegcd(GEN P, GEN Q, GEN T, ulong p)
{
  pari_sp btop, ltop = avma, st_lim;
  long dg;
  GEN U, q;
  if (!signe(P)) return gcopy(Q);
  if (!signe(Q)) return gcopy(P);
  btop = avma; st_lim = stack_lim(btop, 1);
  dg = lg(P)-lg(Q);
  if (dg < 0) { swap(P, Q); dg = -dg; }
  for(;;)
  {
    U = Flxq_invsafe(leading_term(Q), T, p);
    if (!U) { avma = ltop; return NULL; }
    do /* set P := P % Q */
    {
      q = Flxq_mul(U, Flx_neg(leading_term(P), p), T, p);
      P = FlxX_add(P, FlxqX_Flxq_mul(FlxX_shift(Q, dg), q, T, p), p);
      dg = lg(P)-lg(Q);
    } while (dg >= 0);
    if (!signe(P)) break;

    if (low_stack(st_lim, stack_lim(btop, 1)))
    {
      if (DEBUGMEM>1) pari_warn(warnmem,"FlxqX_safegcd");
      gerepileall(btop, 2, &P,&Q);
    }
    swap(P, Q); dg = -dg;
  }
  Q = FlxqX_Flxq_mul(Q, U, T, p); /* normalize GCD */
  return gerepileupto(ltop, Q);
}

/*******************************************************************/
/*                                                                 */
/*                       (Fl[X]/T(X))[Y] / S(Y)                    */
/*                                                                 */
/*******************************************************************/

/*Preliminary implementation to speed up FpX_ffisom*/
typedef struct {
  GEN S, T, mg;
  ulong p;
} FlxYqQ_muldata;

/* reduce x in Fl[X, Y] in the algebra Fl[X, Y]/ (P(X),Q(Y)) */
static GEN
FlxYqQ_redswap(GEN x, GEN S, GEN mg, GEN T, ulong p)
{
  pari_sp ltop=avma;
  long n=degpol(S);
  long m=degpol(T);
  long w = S[1];
  GEN V = FlxX_swap(x,n,w);
  V = FlxqX_red(V,T,p);
  V = FlxX_swap(V,m,w);
  return gerepilecopy(ltop,V); 
}
static GEN
FlxYqQ_sqr(void *data, GEN x)
{
  FlxYqQ_muldata *D = (FlxYqQ_muldata*)data;
  return FlxYqQ_redswap(FlxqX_sqr(x, D->S, D->p),D->S,D->mg,D->T,D->p);
}

static GEN
FlxYqQ_mul(void *data, GEN x, GEN y)
{
  FlxYqQ_muldata *D = (FlxYqQ_muldata*)data;
  return FlxYqQ_redswap(FlxqX_mul(x,y, D->S, D->p),D->S,D->mg,D->T,D->p);
}

/* x in Z[X,Y], S in Z[X] over Fq = Z[Y]/(p,T); compute lift(x^n mod (S,T,p)) */
GEN
FlxYqQ_pow(GEN x, GEN n, GEN S, GEN T, ulong p)
{
  pari_sp av = avma;
  FlxYqQ_muldata D;
  GEN y;
  D.S = S;
  D.T = T;
  D.p = p;
  y = leftright_pow(x, n, (void*)&D, &FlxYqQ_sqr, &FlxYqQ_mul);
  return gerepileupto(av, y);
}

typedef struct {
  GEN T, S;
  GEN mg;
  ulong p;
} kronecker_muldata;

static GEN
FlxqXQ_red(void *data, GEN x)
{
  kronecker_muldata *D = (kronecker_muldata*)data;
  GEN t = FlxqX_from_Kronecker(x, D->T,D->p);
  t = FlxqX_divrem(t, D->S,D->T,D->p, ONLY_REM);
  return FlxX_to_Kronecker(t,D->T);
}
static GEN
FlxqXQ_mul(void *data, GEN x, GEN y) {
  return FlxqXQ_red(data, Flx_mul(x,y,((kronecker_muldata*) data)->p));
}
static GEN
FlxqXQ_sqr(void *data, GEN x) {
  return FlxqXQ_red(data, Flx_sqr(x,((kronecker_muldata*) data)->p));
}

#if 0
static GEN
FlxqXQ_red_montgomery(void *data, GEN x)
{
  kronecker_muldata *D = (kronecker_muldata*)data;
  GEN t = FlxqX_from_Kronecker(x, D->T,D->p);
  t = FlxqX_rem_montgomery(t, D->S,D->mg,D->T,D->p);
  return FlxX_to_Kronecker(t,D->T);
}
static GEN
FlxqXQ_mul_montgomery(void *data, GEN x, GEN y) {
  return FlxqXQ_red_montgomery(data, Flx_mul(x,y,((kronecker_muldata*) data)->p));
}
static GEN
FlxqXQ_sqr_montgomery(void *data, GEN x) {
  return FlxqXQ_red_montgomery(data, Flx_sqr(x,((kronecker_muldata*) data)->p));
}
#endif

long FlxqXQ_POW_MONTGOMERY_LIMIT = 0;

/* x over Fq, return lift(x^n) mod S */
GEN
FlxqXQ_pow(GEN x, GEN n, GEN S, GEN T, ulong p)
{
  pari_sp av0 = avma;
  GEN y;
  kronecker_muldata D;
  D.S = S;
  D.T = T;
  D.p = p;
#if 0
  if (lgpol(S[2]) && degpol(S) >= FlxqXQ_POW_MONTGOMERY_LIMIT)
  {
    /* We do not handle polynomials multiple of x yet */
    D.mg  = FlxqX_invmontgomery(S,T,p);
    y = leftright_pow(FlxX_to_Kronecker(x,T), n, 
        (void*)&D, &FlxqXQ_sqr_montgomery, &FlxqXQ_mul_montgomery);
    y = leftright_pow(FlxX_to_Kronecker(x,T), n, 
        (void*)&D, &FlxqXQ_sqr, &FlxqXQ_mul);
  }
  else
#endif
    y = leftright_pow(FlxX_to_Kronecker(x,T), n, 
        (void*)&D, &FlxqXQ_sqr, &FlxqXQ_mul);
  y = FlxqX_from_Kronecker(y, T,p);
  return gerepileupto(av0, y);
}

struct _FlxqX {ulong p; GEN T;};
static GEN _FlxqX_mul(void *data,GEN a,GEN b)
{
  struct _FlxqX *d=(struct _FlxqX*)data;
  return FlxqX_mul(a,b,d->T,d->p);
}

GEN 
FlxqXV_prod(GEN V, GEN T, ulong p)
{
  struct _FlxqX d; d.p=p; d.T=T;
  return divide_conquer_assoc(V, &_FlxqX_mul, (void*)&d);
}

GEN
FlxqV_roots_to_pol(GEN V, GEN T, ulong p, long v)
{
  pari_sp ltop = avma;
  long k;
  GEN W = cgetg(lg(V),t_VEC);
  for(k=1; k < lg(V); k++)
    gel(W,k) = deg1pol_i(Fl_to_Flx(1,T[1]),Flx_neg(gel(V,k),p),v);
  return gerepileupto(ltop, FlxqXV_prod(W, T, p));
}

