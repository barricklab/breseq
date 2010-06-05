/* $Id: polarit3.c 8320 2007-03-04 22:14:36Z bill $

Copyright (C) 2000-2005  The PARI group.

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
/**                         (third part)                              **/
/**                                                                   **/
/***********************************************************************/
#include "pari.h"
#include "paripriv.h"

/*Renormalize (in place) polynomial with t_INT or t_POL coefficients.*/

GEN
ZX_renormalize(GEN x, long lx)
{
  long i;
  for (i = lx-1; i>1; i--)
    if (signe(gel(x,i))) break;
  stackdummy((pari_sp)(x + lg(x)), (pari_sp)(x + (i+1)));
  setlg(x, i+1); setsigne(x, i!=1); return x;
}

GEN
ZX_add(GEN x, GEN y)
{
  long lx,ly,i;
  GEN z;
  lx = lg(x); ly = lg(y); if (lx < ly) swapspec(x,y, lx,ly);
  z = cgetg(lx,t_POL); z[1] = x[1];
  for (i=2; i<ly; i++) gel(z,i) = addii(gel(x,i),gel(y,i));
  for (   ; i<lx; i++) gel(z,i) = icopy(gel(x,i));
  z = ZX_renormalize(z, lx);
  if (!lgpol(z)) { avma = (pari_sp)(z + lx); return zeropol(varn(x)); }
  return z;
}

GEN
ZX_sub(GEN x,GEN y)
{
  long lx,ly,i,lz;
  GEN z;
  lx = lg(x); ly = lg(y);
  lz=max(lx,ly);
  z = cgetg(lz,t_POL);
  if (lx >= ly)
  {
    z[1] = x[1];
    for (i=2; i<ly; i++) gel(z,i) = subii(gel(x,i),gel(y,i));
    for (   ; i<lx; i++) gel(z,i) = icopy(gel(x,i));
    if (lx == ly) z = ZX_renormalize(z, lz);
  }
  else
  {
    z[1] = y[1];
    for (i=2; i<lx; i++) gel(z,i) = subii(gel(x,i),gel(y,i));
    for (   ; i<ly; i++) gel(z,i) = negi(gel(y,i));
  }
  if (!lgpol(z)) { avma = (pari_sp)(z + lz); z = zeropol(varn(x)); }
  return z;
}

GEN
ZX_neg(GEN x)
{
  long i,d=lg(x);
  GEN y;
  y=cgetg(d,t_POL); y[1]=x[1];
  for(i=2;i<d;i++)
    gel(y,i) = negi(gel(x,i));
  return y;
}

GEN
ZX_Z_add(GEN y, GEN x)
{
  GEN z;
  long lz, i;
  if (!signe(y))
    return scalarpol(x,varn(y));
  lz=lg(y);
  z=cgetg(lz,t_POL);
  z[1]=y[1];
  gel(z,2) = addii(gel(y,2),x);
  for(i=3;i<lz;i++)
    gel(z,i) = icopy(gel(y,i));
  if (lz==3) z = ZX_renormalize(z,lz);
  return z;
}

/*ZX_mul and ZX_sqr are alias for RgX_mul and Rgx_sqr currently*/

GEN
ZX_Z_mul(GEN y,GEN x)
{
  GEN z;
  long i;
  if (!signe(x)) 
    return zeropol(varn(y));
  z=cgetg(lg(y),t_POL);
  z[1]=y[1];
  for(i=2;i<lg(y);i++)
    gel(z,i) = mulii(gel(y,i),x);
  return z;
}

/************************************************************************
 **                                                                    ** 
 **                      Arithmetic in Z/pZ[X]                         **
 **                                                                    **
 ************************************************************************/

/* In practice, p is not assumed to be prime. */

/* p > 0 a t_INT, return lift(x * Mod(1,p)).
 * If x is an INTMOD, assume modulus is a multiple of p. */
GEN
Rg_to_Fp(GEN x, GEN p)
{
  if (lgefint(p) == 3) return utoi(Rg_to_Fl(x, (ulong)p[2]));
  switch(typ(x))
  {
    case t_INT: return modii(x, p);
    case t_FRAC: {
      pari_sp av = avma;
      GEN z = modii(gel(x,1), p);
      if (z == gen_0) return gen_0;
      return gerepileuptoint(av, remii(mulii(z, Fp_inv(gel(x,2), p)), p));
    }
    case t_PADIC: return padic_to_Fp(x, p);
    case t_INTMOD: {
      GEN q = gel(x,1), a = gel(x,2);
      if (equalii(q, p)) return icopy(a);
      return remii(a, p);
    }
    default: pari_err(typeer, "Rg_to_Fp");
      return NULL; /* not reached */
  }
}
ulong
Rg_to_Fl(GEN x, ulong p)
{
  switch(typ(x))
  {
    case t_INT: return umodiu(x, p);
    case t_FRAC: {
      ulong z = umodiu(gel(x,1), p);
      if (!z) return 0;
      return Fl_div(z, umodiu(gel(x,2), p), p);
    }
    case t_PADIC: return padic_to_Fl(x, p);
    case t_INTMOD: {
      GEN q = gel(x,1), a = gel(x,2);
      if (equaliu(q, p)) return itou(a);
      return umodiu(a, p);
    }
    default: pari_err(typeer, "Rg_to_Fl");
      return 0; /* not reached */
  }
}
/* If x is a POLMOD, assume modulus is a multiple of T. */
GEN
Rg_to_FpXQ(GEN x, GEN T, GEN p)
{
  long ta, tx = typ(x), v = varn(T);
  GEN a, b;
  if (is_const_t(tx)) return scalarpol(Rg_to_Fp(x, p), v);
  switch(tx)
  {
    case t_POLMOD:
      b = gel(x,1);
      a = gel(x,2); ta = typ(a);
      if (is_const_t(ta)) return Rg_to_Fp(a, p);
      b = RgX_to_FpX(b, p); if (varn(b) != v) break;
      a = RgX_to_FpX(a, p); if (gequal(b,T)) return a;
      return FpX_rem(a, T, p);
    case t_POL:
      if (varn(x) != v) break;
      return FpX_rem(RgX_to_FpX(x,p), T, p);
    case t_RFRAC:
      a = Rg_to_FpXQ(gel(x,1), T,p);
      b = Rg_to_FpXQ(gel(x,2), T,p);
      return FpXQ_div(a,b, T,p);
  }
  pari_err(typeer,"Rg_to_FpXQ");
  return NULL; /* not reached */
}
GEN
RgX_to_FpX(GEN x, GEN p)
{
  long i, l = lg(x);
  GEN z = cgetg(l, t_POL); z[1] = x[1];
  for (i = 2; i < l; i++) gel(z,i) = Rg_to_Fp(gel(x,i), p);
  return normalizepol_i(z, l);
}

GEN
RgV_to_FpV(GEN x, GEN p)
{
  long i, l = lg(x);
  GEN z = cgetg(l, t_VEC);
  for (i = 1; i < l; i++) gel(z,i) = Rg_to_Fp(gel(x,i), p);
  return z;
}

GEN
RgC_to_FpC(GEN x, GEN p)
{
  long i, l = lg(x);
  GEN z = cgetg(l, t_COL);
  for (i = 1; i < l; i++) gel(z,i) = Rg_to_Fp(gel(x,i), p);
  return z;
}

GEN
RgX_to_FpXQX(GEN x, GEN T, GEN p)
{
  long i, l = lg(x);
  GEN z = cgetg(l, t_POL); z[1] = x[1];
  for (i = 2; i < l; i++) gel(z,i) = Rg_to_FpXQ(gel(x,i), T,p);
  return normalizepol_i(z, l);
}
GEN
RgX_to_FqX(GEN x, GEN T, GEN p)
{
  long i, l = lg(x);
  GEN z = cgetg(l, t_POL); z[1] = x[1];
  for (i = 2; i < l; i++) gel(z,i) = simplify_i(Rg_to_FpXQ(gel(x,i), T,p));
  return normalizepol_i(z, l);
}

/*********************************************************************
These functions suppose polynomials to be already reduced.
They are clean and memory efficient.
**********************************************************************/

GEN
FpX_center(GEN T,GEN mod)
{/*OK centermod exists, but is not so clean*/
  pari_sp av;
  long i, l=lg(T);
  GEN P,mod2;
  P=cgetg(l,t_POL);
  P[1]=T[1];
  av=avma;
  mod2=gclone(shifti(mod,-1));/*clone*/
  avma=av;
  for(i=2;i<l;i++)
    gel(P,i) = cmpii(gel(T,i),mod2)<=0? icopy(gel(T,i)): subii(gel(T,i),mod);
  gunclone(mod2);/*unclone*/
  return P;
}

GEN
FpX_neg(GEN x,GEN p)
{
  long i,d=lg(x);
  GEN y;
  y=cgetg(d,t_POL); y[1]=x[1];
  for(i=2;i<d;i++)
    if (signe(x[i])) gel(y,i) = subii(p,gel(x,i));
    else gel(y,i) = gen_0;
  return y;
}
/**********************************************************************
Unclean functions, do not garbage collect.
This is a feature: The stack is corrupted only by the call to FpX_red
so garbage collecting so often is not desirable.
FpX_red can sometime be avoided by passing NULL for p.
In this case the function is usually clean (see below for detail)
Added to help not using POLMOD of INTMOD which are deadly slow.
gerepileupto of the result is legible.   Bill.
I don't like C++.  I am wrong.
**********************************************************************/
/*
 *If p is NULL no reduction is performed and the function is clean.
 * for FpX_add,FpX_mul,FpX_sqr,FpX_Fp_mul
 */
GEN
FpX_add(GEN x,GEN y,GEN p)
{
  GEN z = ZX_add(x,y);
  return p? FpX_red(z, p): z;
}

GEN
FpX_sub(GEN x,GEN y,GEN p)
{
  GEN z = ZX_sub(x,y);
  return p? FpX_red(z, p): z;
}

GEN
FpX_mul(GEN x,GEN y,GEN p)
{
  GEN z = ZX_mul(x, y);
  return p? FpX_red(z, p): z;
}

GEN
FpX_sqr(GEN x,GEN p)
{
  GEN z = ZX_sqr(x);
  return p? FpX_red(z, p): z;
}

GEN
FpX_Fp_mul(GEN x,GEN y,GEN p)
{
  GEN z = ZX_Z_mul(x,y);
  return p? FpX_red(z, p): z;
}

/* Inverse of x in Z/pZ[X]/(pol) or NULL if inverse doesn't exist
 * return lift(1 / (x mod (p,pol))) */
GEN
FpXQ_invsafe(GEN x, GEN T, GEN p)
{
  GEN z, U, V;

  z = FpX_extgcd(x, T, p, &U, &V);
  if (degpol(z)) return NULL;
  z = Fp_invsafe(gel(z,2), p);
  if (!z) return NULL;
  return FpX_Fp_mul(U, z, p);
}

/* Product of y and x in Z/pZ[X]/(T)
 * return lift(lift(Mod(x*y*Mod(1,p),T*Mod(1,p)))) */
/* x and y must be polynomials in the same var as T.
 * t_INT are not allowed. Use Fq_mul instead.
 */
GEN
FpXQ_mul(GEN y,GEN x,GEN T,GEN p)
{
  GEN z = RgX_mulspec(y+2, x+2, lgpol(y), lgpol(x)); setvarn(z,varn(y));
  z = FpX_red(z, p); return FpX_rem(z,T, p);
}

/* Square of y in Z/pZ[X]/(pol)
 * return lift(lift(Mod(y^2*Mod(1,p),pol*Mod(1,p)))) */
GEN
FpXQ_sqr(GEN y,GEN pol,GEN p)
{
  GEN z = RgX_sqrspec(y+2,lgpol(y)); setvarn(z,varn(y));
  z = FpX_red(z, p); return FpX_rem(z,pol, p);
}
/*Modify y[2].
 *No reduction if p is NULL
 */
GEN
FpX_Fp_add(GEN y,GEN x,GEN p)
{
  if (!signe(x)) return y;
  if (!signe(y))
    return scalarpol(x,varn(y));
  gel(y,2) = addii(gel(y,2),x);
  if (p) gel(y,2) = modii(gel(y,2),p);
  if (!signe(y[2]) && degpol(y) == 0) return zeropol(varn(y));
  return y;
}
/* as above over Fp[X] */
GEN
FpX_rescale(GEN P, GEN h, GEN p)
{
  long i, l = lg(P);
  GEN Q = cgetg(l,t_POL), hi = h;
  Q[l-1] = P[l-1];
  for (i=l-2; i>=2; i--)
  {
    gel(Q,i) = modii(mulii(gel(P,i), hi), p);
    if (i == 2) break;
    hi = modii(mulii(hi,h), p);
  }
  Q[1] = P[1]; return Q;
}
/*****************************************************************
 *                 End of unclean functions.                     *
 *****************************************************************/

/*****************************************************************
 Clean and with no reduced hypothesis.  Beware that some operations
 will be much slower with big unreduced coefficient
*****************************************************************/
/* Inverse of x in Z[X] / (p,T)
 * return lift(lift(Mod(x*Mod(1,p), T*Mod(1,p))^-1)); */
GEN
FpXQ_inv(GEN x,GEN T,GEN p)
{
  pari_sp av = avma;
  GEN U = FpXQ_invsafe(x, T, p);
  if (!U) pari_err(talker,"non invertible polynomial in FpXQ_inv");
  return gerepileupto(av, U);
}

GEN
FpXQ_div(GEN x,GEN y,GEN T,GEN p)
{
  pari_sp av = avma;
  return gerepileupto(av, FpXQ_mul(x,FpXQ_inv(y,T,p),T,p));
}

GEN
FpXV_FpC_mul(GEN V, GEN W, GEN p)
{
  pari_sp ltop=avma;
  long i;
  GEN z = ZX_Z_mul(gel(V,1),gel(W,1));
  for(i=2;i<lg(V);i++)
    z=ZX_add(z,ZX_Z_mul(gel(V,i),gel(W,i)));
  return gerepileupto(ltop,FpX_red(z,p));
}

/* generates the list of powers of x of degree 0,1,2,...,l*/
GEN
FpXQ_powers(GEN x, long l, GEN T, GEN p)
{
  GEN V=cgetg(l+2,t_VEC);
  long i;
  gel(V,1) = pol_1[varn(T)]; if (l==0) return V;
  gel(V,2) = gcopy(x);       if (l==1) return V;
  if (lgefint(p) == 3) {
    long pp = p[2];
    return FlxC_to_ZXC(Flxq_powers(ZX_to_Flx(x, pp), l, ZX_to_Flx(T,pp), pp));
  }
  gel(V,3) = FpXQ_sqr(x,T,p);
  if ((degpol(x)<<1) < degpol(T)) {
    for(i = 4; i < l+2; i++)
      gel(V,i) = FpXQ_mul(gel(V,i-1),x,T,p);
  } else { /* use squarings if degree(x) is large */
    for(i = 4; i < l+2; i++)
      gel(V,i) = (i&1)? FpXQ_sqr(gel(V, (i+1)>>1),T,p)
                      : FpXQ_mul(gel(V, i-1),x,T,p);
  }
  return V;
}
#if 0
static long brent_kung_nbmul(long d, long n, long p)
{
  return p+n*((d+p-1)/p);
}
  TODO: This the the optimal parameter for the purpose of reducing
  multiplications, but profiling need to be done to ensure it is not slower 
  than the other option in practice
/*Return optimal parameter l for the evaluation of n polynomials of degree d*/
long brent_kung_optpow(long d, long n)
{
  double r;
  long f,c,pr;
  if (n>=d ) return d;
  pr=n*d;
  if (pr<=1) return 1;
  r=d/sqrt(pr);
  c=(long)ceil(d/ceil(r));
  f=(long)floor(d/floor(r));
  return (brent_kung_nbmul(d, n, c) <= brent_kung_nbmul(d, n, f))?c:f;
}
#endif 
/*Return optimal parameter l for the evaluation of n polynomials of degree d*/
long
brent_kung_optpow(long d, long n)
{
  long l, pr;
  if (n >= d) return d;
  pr = n*d; if (pr <= 1) return 1;
  l = (long) ((double)d / sqrt(pr));
  return (d+l-1) / l;
}

/*Close to FpXV_FpC_mul*/

static GEN
spec_compo_powers(GEN P, GEN V, long a, long n)
{
  long i;
  GEN z;
  z = scalarpol(gel(P,2+a),varn(P));
  for(i=1;i<=n;i++)
    z=ZX_add(z,ZX_Z_mul(gel(V,i+1),gel(P,2+a+i)));
  return z;
}
/*Try to implement algorithm in Brent & Kung (Fast algorithms for
 *manipulating formal power series, JACM 25:581-595, 1978)
 
 V must be as output by FpXQ_powers.
 For optimal performance, l (of FpXQ_powers) must be as output by
 brent_kung_optpow
 */

GEN
FpX_FpXQV_compo(GEN P, GEN V, GEN T, GEN p)
{
  pari_sp ltop=avma;
  long l=lg(V)-1;
  GEN z,u;
  long d=degpol(P),cnt=0;
  if (d==-1) return zeropol(varn(T));
  if (d<l)
  {
    z=spec_compo_powers(P,V,0,d);
    return gerepileupto(ltop,FpX_red(z,p));
  }
  if (l<=1)
    pari_err(talker,"powers is only [] or [1] in FpX_FpXQV_compo");
  z=spec_compo_powers(P,V,d-l+1,l-1);
  d-=l;
  while(d>=l-1)
  {
    u=spec_compo_powers(P,V,d-l+2,l-2);
    z=ZX_add(u,FpXQ_mul(z,gel(V,l),T,p));
    d-=l-1;
    cnt++;
  }
  u=spec_compo_powers(P,V,0,d);
  z=ZX_add(u,FpXQ_mul(z,gel(V,d+2),T,p));
  cnt++;
  if (DEBUGLEVEL>=8) fprintferr("FpX_FpXQV_compo: %d FpXQ_mul [%d]\n",cnt,l-1);
  return gerepileupto(ltop,FpX_red(z,p));
}

/* T in Z[X] and  x in Z/pZ[X]/(pol)
 * return lift(lift(subst(T,variable(T),Mod(x*Mod(1,p),pol*Mod(1,p)))));
 */
GEN
FpX_FpXQ_compo(GEN T,GEN x,GEN pol,GEN p)
{
  pari_sp ltop=avma;
  GEN z;
  long d=degpol(T),rtd;
  if (!signe(T)) return zeropol(varn(T));
  rtd = (long) sqrt((double)d);
  z = FpX_FpXQV_compo(T,FpXQ_powers(x,rtd,pol,p),pol,p);
  return gerepileupto(ltop,z);
}

/* Evaluation in Fp
 * x a ZX and y an Fp, return x(y) mod p
 *
 * If p is very large (several longs) and x has small coefficients(<<p),
 * then Brent & Kung algorithm is faster. */
GEN
FpX_eval(GEN x,GEN y,GEN p)
{
  pari_sp av;
  GEN p1,r,res;
  long j, i=lg(x)-1;
  if (i<=2)
    return (i==2)? modii(gel(x,2),p): gen_0;
  res=cgeti(lgefint(p));
  av=avma; p1=gel(x,i);
  /* specific attention to sparse polynomials (see poleval)*/
  /*You've guessed it! It's a copy-paste(tm)*/
  for (i--; i>=2; i=j-1)
  {
    for (j=i; !signe(gel(x,j)); j--)
      if (j==2)
      {
	if (i!=j) y = Fp_powu(y,i-j+1,p);
	p1=mulii(p1,y);
	goto fppoleval;/*sorry break(2) no implemented*/
      }
    r = (i==j)? y: Fp_powu(y,i-j+1,p);
    p1 = modii(addii(mulii(p1,r), gel(x,j)),p);
  }
 fppoleval:
  modiiz(p1,p,res);
  avma=av;
  return res;
}
GEN
FqX_eval(GEN x, GEN y, GEN T, GEN p)
{
  pari_sp av;
  GEN p1, r;
  long j, i=lg(x)-1;
  if (i<=2)
    return (i==2)? Fq_red(gel(x,2), T, p): gen_0;
  av=avma; p1=gel(x,i);
  /* specific attention to sparse polynomials (see poleval)*/
  /*You've guessed it! It's a copy-paste(tm)*/
  for (i--; i>=2; i=j-1)
  {
    for (j=i; !signe(gel(x,j)); j--)
      if (j==2)
      {
	if (i!=j) y = Fq_pow(y,utoipos(i-j+1), T, p);
        return gerepileupto(av, gmul(p1,y));
      }
    r = (i==j)? y: Fq_pow(y, utoipos(i-j+1), T, p);
    p1 = Fq_red(gadd(gmul(p1,r), gel(x,j)), T, p);
  }
  return gerepileupto(av, p1);
}
/* Tz=Tx*Ty where Tx and Ty coprime
 * return lift(chinese(Mod(x*Mod(1,p),Tx*Mod(1,p)),Mod(y*Mod(1,p),Ty*Mod(1,p))))
 * if Tz is NULL it is computed
 * As we do not return it, and the caller will frequently need it,
 * it must compute it and pass it.
 */
GEN
FpX_chinese_coprime(GEN x,GEN y,GEN Tx,GEN Ty,GEN Tz,GEN p)
{
  pari_sp av = avma;
  GEN ax,p1;
  ax = FpX_mul(FpXQ_inv(Tx,Ty,p), Tx,p);
  p1=FpX_mul(ax, FpX_sub(y,x,p),p);
  p1 = FpX_add(x,p1,p);
  if (!Tz) Tz=FpX_mul(Tx,Ty,p);
  p1 = FpX_rem(p1,Tz,p);
  return gerepileupto(av,p1);
}

typedef struct {
  GEN pol, p;
} FpX_muldata;

static GEN
_sqr(void *data, GEN x)
{
  FpX_muldata *D = (FpX_muldata*)data;
  return FpXQ_sqr(x, D->pol, D->p);
}
static GEN
_mul(void *data, GEN x, GEN y)
{
  FpX_muldata *D = (FpX_muldata*)data;
  return FpXQ_mul(x,y, D->pol, D->p);
}

/* x,pol in Z[X], p in Z, n in Z, compute lift(x^n mod (p, pol)) */
GEN
FpXQ_pow(GEN x, GEN n, GEN pol, GEN p)
{
  FpX_muldata D;
  pari_sp av;
  GEN y;

  if (!signe(n)) return pol_1[ varn(x) ];
  if (is_pm1(n)) /* +/- 1 */
    return (signe(n) < 0)? FpXQ_inv(x,pol,p): gcopy(x);
  av = avma;
  if (!is_bigint(p))
  {
    ulong pp = p[2];
    pol = ZX_to_Flx(pol, pp);
    x   = ZX_to_Flx(x,   pp);
    y = Flx_to_ZX( Flxq_pow(x, n, pol, pp) );
  }
  else
  {
    D.pol = pol;
    D.p   = p;
    if (signe(n) < 0) x = FpXQ_inv(x,pol,p);
    y = leftright_pow(x, n, (void*)&D, &_sqr, &_mul);
  }
  return gerepileupto(av, y);
}

static GEN _FpX_mul(void *p,GEN a,GEN b){return FpX_mul(a,b,(GEN)p);}
GEN 
FpXV_prod(GEN V, GEN p)
{
  return divide_conquer_assoc(V, &_FpX_mul,(void *)p);
}

GEN
FpV_roots_to_pol(GEN V, GEN p, long v)
{
  pari_sp ltop=avma;
  long i;
  GEN g=cgetg(lg(V),t_VEC);
  for(i=1;i<lg(V);i++)
    gel(g,i) = deg1pol_i(gen_1,modii(negi(gel(V,i)),p),v);
  return gerepileupto(ltop,FpXV_prod(g,p));
}

/*******************************************************************/
/*                                                                 */
/*                             FpXX                                */
/*                                                                 */
/*******************************************************************/
/*Polynomials whose coefficients are either polynomials or integers*/
GEN
FpXX_red(GEN z, GEN p)
{
  GEN res;
  long i;
  res = cgetg(lg(z),t_POL); res[1] = z[1];
  for(i=2;i<lg(res);i++)
    if (typ(z[i])==t_INT)
      gel(res,i) = modii(gel(z,i),p);
    else
    {
      pari_sp av=avma;
      gel(res,i) = FpX_red(gel(z,i),p);
      if (lg(res[i])<=3)
      {
        if (lg(res[i])==2) {avma=av;gel(res,i) = gen_0;}
        else gel(res,i) = gerepilecopy(av,gmael(res,i,2));
      }
    }
  return FpXX_renormalize(res,lg(res));
}
GEN
FpXX_add(GEN x, GEN y, GEN p)
{
  long i,lz;
  GEN z; 
  long lx=lg(x);
  long ly=lg(y);
  if (ly>lx) swapspec(x,y, lx,ly);
  lz = lx; z = cgetg(lz, t_POL); z[1]=x[1];
  for (i=2; i<ly; i++) gel(z,i) = Fq_add(gel(x,i), gel(y,i), NULL, p);
  for (   ; i<lx; i++) gel(z,i) = gcopy(gel(x,i));
  return FpXX_renormalize(z, lz);
}

/*******************************************************************/
/*                                                                 */
/*                             (Fp[X]/(Q))[Y]                      */
/*                                                                 */
/*******************************************************************/
/*Not malloc nor warn-clean.*/
GEN
FpXQX_from_Kronecker(GEN Z, GEN T, GEN p)
{
  long i,j,lx,l, N = (degpol(T)<<1) + 1;
  GEN x, t = cgetg(N,t_POL), z = FpX_red(Z, p);
  t[1] = T[1] & VARNBITS;
  l = lg(z); lx = (l-2) / (N-2);
  x = cgetg(lx+3,t_POL);
  for (i=2; i<lx+2; i++)
  {
    for (j=2; j<N; j++) t[j] = z[j];
    z += (N-2);
    gel(x,i) = FpX_rem(FpX_renormalize(t,N), T,p);
  }
  N = (l-2) % (N-2) + 2;
  for (j=2; j<N; j++) t[j] = z[j];
  gel(x,i) = FpX_rem(FpX_renormalize(t,N), T,p);
  return FpXQX_renormalize(x, i+1);
}

GEN
FqX_red(GEN z, GEN T, GEN p) { return T? FpXQX_red(z, T, p): FpXX_red(z, p); }
GEN
FpXQX_red(GEN z, GEN T, GEN p)
{
  long i, l = lg(z);
  GEN res = cgetg(l,t_POL); res[1] = z[1];
  for(i=2;i<l;i++)
    if (typ(z[i]) == t_INT)
      gel(res,i) = modii(gel(z,i),p);
    else
      gel(res,i) = FpX_rem(gel(z,i),T,p);
  return FpXQX_renormalize(res,lg(res));
}

GEN
FpXQX_mul(GEN x, GEN y, GEN T, GEN p)
{
  pari_sp ltop=avma;
  GEN z,kx,ky;
  long vx = min(varn(x),varn(y));
  kx= to_Kronecker(x,T);
  ky= to_Kronecker(y,T);
  z = RgX_mulspec(ky+2, kx+2, lgpol(ky), lgpol(kx));
  z = FpXQX_from_Kronecker(z,T,p);
  setvarn(z,vx);/*RgX_mulspec and FpXQX_from_Kronecker are not varn-clean*/
  return gerepileupto(ltop,z);
}
GEN
FpXQX_sqr(GEN x, GEN T, GEN p)
{
  pari_sp ltop=avma;
  GEN z,kx;
  long vx=varn(x);
  kx= to_Kronecker(x,T);
  z = RgX_sqrspec(kx+2, lgpol(kx));
  z = FpXQX_from_Kronecker(z,T,p);
  setvarn(z,vx);/*RgX_mulspec and FpXQX_from_Kronecker are nor varn-clean*/
  return gerepileupto(ltop,z);
}

GEN
FqX_Fq_mul(GEN P, GEN U, GEN T, GEN p)
{
  long i, lP = lg(P);
  GEN res = cgetg(lP,t_POL); res[1] = P[1];
  for(i=2; i<lP; i++) gel(res,i) = Fq_mul(U,gel(P,i), T,p);
  return FpXQX_renormalize(res,lg(res));
}

/* a X^d */
GEN
monomial(GEN a, long d, long v)
{
  long i, lP = d+3;
  GEN P;
  if (d < 0) {
    P = cgetg(3, t_RFRAC);
    gel(P,1) = a;
    gel(P,2) = monomial(gen_1, -d, v);
  } else {
    P = cgetg(lP, t_POL);
    if (gcmp0(a)) P[1] = evalsigne(0) | evalvarn(v);
    else          P[1] = evalsigne(1) | evalvarn(v);
    lP--; gel(P,lP) = a;
    for (i=2; i<lP; i++) gel(P,i) = gen_0;
  }
  return P;
}
GEN
monomialcopy(GEN a, long d, long v)
{
  long i, lP = d+3;
  GEN P;
  if (d < 0) {
    P = cgetg(3, t_RFRAC);
    gel(P,1) = gcopy(a);
    gel(P,2) = monomial(gen_1, -d, v);
  } else {
    P = cgetg(lP, t_POL);
    if (gcmp0(a)) P[1] = evalsigne(0) | evalvarn(v);
    else          P[1] = evalsigne(1) | evalvarn(v);
    lP--; gel(P,lP) = gcopy(a);
    for (i=2; i<lP; i++) gel(P,i) = gen_0;
  }
  return P;
}

GEN
FpXQX_gcd(GEN P, GEN Q, GEN T, GEN p)
{
  pari_sp av2, av = avma, st_lim;
  long dg;
  GEN U, q;
  if (lgefint(p) == 3)
  {
    ulong pp = (ulong)p[2];
    GEN Pl, Ql, Tl;
    Pl = ZXX_to_FlxX(P, pp, varn(T));
    if (!signe(Pl)) { avma = av; return gcopy(Q); }
    Ql = ZXX_to_FlxX(Q, pp, varn(T));
    if (!signe(Ql)) { avma = av; return gcopy(P); }
    Tl = ZX_to_Flx(T, pp);
    U = FlxqX_safegcd(Pl, Ql, Tl, pp);
    if (!U) pari_err(talker, "non-invertible polynomial in FpXQX_gcd");
    return gerepileupto(av, FlxX_to_ZXX(U));
  }
  P = FpXX_red(P, p); av2 = avma;
  Q = FpXX_red(Q, p);
  if (!signe(P)) return gerepileupto(av, Q);
  if (!signe(Q)) { avma = av2; return P; }
  T = FpX_red(T, p);

  av2 = avma; st_lim = stack_lim(av2, 1);
  dg = lg(P)-lg(Q);
  if (dg < 0) { swap(P, Q); dg = -dg; }
  for(;;)
  {
    U = Fq_inv(leading_term(Q), T, p);
    do /* set P := P % Q */
    {
      q = Fq_mul(U, Fq_neg(leading_term(P), T, p), T, p);
      P = FpXX_add(P, FqX_Fq_mul(RgX_shift_shallow(Q, dg), q, T, p), p);
      dg = lg(P)-lg(Q);
    } while (dg >= 0);
    if (!signe(P)) break;

    if (low_stack(st_lim, stack_lim(av2, 1)))
    {
      if (DEBUGMEM>1) pari_warn(warnmem,"FpXQX_gcd");
      gerepileall(av2, 2, &P,&Q);
    }
    swap(P, Q); dg = -dg;
  }
  Q = FqX_Fq_mul(Q, U, T, p); /* normalize GCD */
  return gerepileupto(av, Q);
}

/*******************************************************************/
/*                                                                 */
/*                       (Fp[X]/T(X))[Y] / S(Y)                    */
/*                                                                 */
/*******************************************************************/

/*Preliminary implementation to speed up FpX_ffisom*/
typedef struct {
  GEN S, T, p;
} FpXYQQ_muldata;

/* reduce x in Fp[X, Y] in the algebra Fp[X, Y]/ (P(X),Q(Y)) */
static GEN
FpXYQQ_redswap(GEN x, GEN S, GEN T, GEN p)
{
  pari_sp ltop=avma;
  long n=degpol(S);
  long m=degpol(T);
  long v=varn(T),w=varn(S);
  GEN V = RgXY_swap(x,n,w);
  setvarn(T,w);
  V = FpXQX_red(V,T,p);
  setvarn(T,v);
  V = RgXY_swap(V,m,w);
  return gerepilecopy(ltop,V); 
}
static GEN
FpXYQQ_sqr(void *data, GEN x)
{
  FpXYQQ_muldata *D = (FpXYQQ_muldata*)data;
  return FpXYQQ_redswap(FpXQX_sqr(x, D->S, D->p),D->S,D->T,D->p);
  
}
static GEN
FpXYQQ_mul(void *data, GEN x, GEN y)
{
  FpXYQQ_muldata *D = (FpXYQQ_muldata*)data;
  return FpXYQQ_redswap(FpXQX_mul(x,y, D->S, D->p),D->S,D->T,D->p);
}

/* x in Z[X,Y], S in Z[X] over Fq = Z[Y]/(p,T); compute lift(x^n mod (S,T,p)) */
GEN
FpXYQQ_pow(GEN x, GEN n, GEN S, GEN T, GEN p)
{
  pari_sp av = avma;
  FpXYQQ_muldata D;
  GEN y;
  if (OK_ULONG(p))
  {
    ulong pp = p[2];
    x = ZXX_to_FlxX(x, pp, varn(T));
    S = ZX_to_Flx(S, pp);
    T = ZX_to_Flx(T, pp);
    y = FlxX_to_ZXX( FlxYqQ_pow(x, n, S, T, pp) );
  }
  else
  {
    D.S = S;
    D.T = T;
    D.p = p;
    y = leftright_pow(x, n, (void*)&D, &FpXYQQ_sqr, &FpXYQQ_mul);
  }
  return gerepileupto(av, y);
}

typedef struct {
  GEN T, p, S;
  long v;
} kronecker_muldata;

static GEN
FpXQYQ_red(void *data, GEN x)
{
  kronecker_muldata *D = (kronecker_muldata*)data;
  GEN t = FpXQX_from_Kronecker(x, D->T,D->p);
  setvarn(t, D->v);
  t = FpXQX_divrem(t, D->S,D->T,D->p, ONLY_REM);
  return to_Kronecker(t,D->T);
}
static GEN
FpXQYQ_mul(void *data, GEN x, GEN y) {
  return FpXQYQ_red(data, ZX_mul(x,y));
}
static GEN
FpXQYQ_sqr(void *data, GEN x) {
  return FpXQYQ_red(data, ZX_sqr(x));
}

/* x over Fq, return lift(x^n) mod S */
GEN
FpXQYQ_pow(GEN x, GEN n, GEN S, GEN T, GEN p)
{
  pari_sp ltop = avma;
  GEN y;
  kronecker_muldata D;
  if (OK_ULONG(p))
  {
    ulong pp = p[2];
    GEN z;
    long v = varn(T);
    T = ZX_to_Flx(T, pp);
    x = ZXX_to_FlxX(x, pp, v);
    S = ZXX_to_FlxX(S, pp, v);
    z = FlxqXQ_pow(x, n, S, T, pp);
    y = FlxX_to_ZXX(z);
  }
  else
  {
    long v = varn(x);
    D.S = S;
    D.T = T;
    D.p = p;
    D.v = v;
    y = leftright_pow(to_Kronecker(x,T), n, (void*)&D, &FpXQYQ_sqr, &FpXQYQ_mul);
    y = FpXQX_from_Kronecker(y, T,p);
    setvarn(y, v); 
  }
  return gerepileupto(ltop, y);
}
/*******************************************************************/
/*                                                                 */
/*                             Fq                                  */
/*                                                                 */
/*******************************************************************/

GEN
Fq_add(GEN x, GEN y, GEN T/*unused*/, GEN p)
{
  (void)T;
  switch((typ(x)==t_POL)|((typ(y)==t_POL)<<1))
  {
    case 0: return modii(addii(x,y),p);
    case 1: return FpX_Fp_add(x,y,p);
    case 2: return FpX_Fp_add(y,x,p);
    case 3: return FpX_add(x,y,p);
  }
  return NULL;
}

GEN
Fq_sub(GEN x, GEN y, GEN T/*unused*/, GEN p)
{
  (void)T;
  switch((typ(x)==t_POL)|((typ(y)==t_POL)<<1))
  {
    case 0: return modii(subii(x,y),p);
    case 1: return FpX_Fp_add(x,negi(y),p);
    case 2: return FpX_Fp_add(FpX_neg(y,p),x,p);
    case 3: return FpX_sub(x,y,p);
  }
  return NULL;
}

GEN
Fq_neg(GEN x, GEN T/*unused*/, GEN p)
{
  (void)T;
  switch(typ(x)==t_POL)
  {
    case 0: return signe(x)?subii(p,x):gen_0;
    case 1: return FpX_neg(x,p);
  }
  return NULL;
}

/* If T==NULL do not reduce*/
GEN
Fq_mul(GEN x, GEN y, GEN T, GEN p)
{
  switch((typ(x)==t_POL)|((typ(y)==t_POL)<<1))
  {
    case 0: return modii(mulii(x,y),p);
    case 1: return FpX_Fp_mul(x,y,p);
    case 2: return FpX_Fp_mul(y,x,p);
    case 3: if (T) return FpXQ_mul(x,y,T,p);
            else return FpX_mul(x,y,p);
  }
  return NULL;
}

GEN
Fq_neg_inv(GEN x, GEN T, GEN p)
{
  if (typ(x) == t_INT) return Fp_inv(negi(x),p);
  return FpXQ_inv(FpX_neg(x,p),T,p);
}

GEN
Fq_invsafe(GEN x, GEN pol, GEN p)
{
  if (typ(x) == t_INT) return Fp_invsafe(x,p);
  return FpXQ_invsafe(x,pol,p);
}

GEN
Fq_inv(GEN x, GEN pol, GEN p)
{
  if (typ(x) == t_INT) return Fp_inv(x,p);
  return FpXQ_inv(x,pol,p);
}

GEN
Fq_pow(GEN x, GEN n, GEN pol, GEN p)
{
  if (typ(x) == t_INT) return Fp_pow(x,n,p);
  return FpXQ_pow(x,n,pol,p);
}

GEN
Fq_red(GEN x, GEN T, GEN p)
{
  pari_sp ltop=avma;
  switch(typ(x)==t_POL)
  {
    case 0: return modii(x,p);
    case 1: return gerepileupto(ltop,FpX_rem(FpX_red(x,p),T,p));
  }
  return NULL;
}


/*******************************************************************/
/*                                                                 */
/*                             Fq[X]                               */
/*                                                                 */
/*******************************************************************/

GEN
FqX_mul(GEN x, GEN y, GEN T, GEN p)
{
  return T? FpXQX_mul(x, y, T, p): FpX_mul(x, y, p);
}
GEN
FqX_sqr(GEN x, GEN T, GEN p)
{
  return T? FpXQX_sqr(x, T, p): FpX_sqr(x, p);
}
GEN
FqX_div(GEN x, GEN y, GEN T, GEN p)
{
  return T? FpXQX_divrem(x,y,T,p,NULL): FpX_divrem(x,y,p,NULL);
}
GEN
FqX_rem(GEN x, GEN y, GEN T, GEN p)
{
  return T? FpXQX_divrem(x,y,T,p,ONLY_REM): FpX_divrem(x,y,p,ONLY_REM);
}
GEN
FqX_divrem(GEN x, GEN y, GEN T, GEN p, GEN *z)
{
  return T? FpXQX_divrem(x,y,T,p,z): FpX_divrem(x,y,p,z);
}

struct _FpXQX { GEN T,p; };
static GEN _FpXQX_mul(void *data, GEN a,GEN b)
{
  struct _FpXQX *d=(struct _FpXQX*)data;
  return FpXQX_mul(a,b,d->T,d->p);
}
GEN 
FpXQXV_prod(GEN V, GEN T, GEN p)
{
  if (lgefint(p) == 3)
  {
    pari_sp av = avma;
    ulong pp = p[2];
    GEN Tl = ZX_to_Flx(T, pp);
    GEN Vl = ZXXV_to_FlxXV(V, pp, varn(T));
    Tl = FlxqXV_prod(Vl, Tl, pp);
    return gerepileupto(av, FlxX_to_ZXX(Tl));
  }
  else
  {
    struct _FpXQX d;
    d.p=p; 
    d.T=T;
    return divide_conquer_assoc(V, &_FpXQX_mul,(void*)&d);
  }
}

GEN
FqV_roots_to_pol(GEN V, GEN T, GEN p, long v)
{
  pari_sp ltop = avma;
  long k;
  GEN W;
  if (lgefint(p) == 3)
  {
    ulong pp = p[2];
    GEN Tl = ZX_to_Flx(T, pp);
    GEN Vl = FqV_to_FlxV(V, T, p);
    Tl = FlxqV_roots_to_pol(Vl, Tl, pp, v);
    return gerepileupto(ltop, FlxX_to_ZXX(Tl));
  }
  W = cgetg(lg(V),t_VEC);
  for(k=1; k < lg(V); k++)
    gel(W,k) = deg1pol_i(gen_1,Fq_neg(gel(V,k),T,p),v);
  return gerepileupto(ltop, FpXQXV_prod(W, T, p));
}

GEN
FqV_red(GEN z, GEN T, GEN p)
{
  long i, l = lg(z);
  GEN res = cgetg(l, typ(z));
  for(i=1;i<l;i++)
    if (typ(z[i]) == t_INT)
      gel(res,i) = modii(gel(z,i),p);
    else if (T)
      gel(res,i) = FpX_rem(gel(z,i),T,p);
    else
      gel(res,i) = FpX_red(gel(z,i),p);
  return res;
}

GEN
FqV_to_FlxV(GEN v, GEN T, GEN pp)
{
  long j, N = lg(v);
  long vT = varn(T);
  ulong p = pp[2];
  GEN y = cgetg(N, t_VEC);
  for (j=1; j<N; j++) 
    gel(y,j) = (typ(v[j])==t_INT?  Z_to_Flx(gel(v,j), p, vT)
                                  : ZX_to_Flx(gel(v,j), p));
  return y;
}

GEN
FqC_to_FlxC(GEN v, GEN T, GEN pp)
{
  long j, N = lg(v);
  long vT = varn(T);
  ulong p = pp[2];
  GEN y = cgetg(N, t_COL);
  for (j=1; j<N; j++) 
    gel(y,j) = (typ(v[j])==t_INT?  Z_to_Flx(gel(v,j), p, vT)
                                  : ZX_to_Flx(gel(v,j), p));
  return y;
}

GEN
FqM_to_FlxM(GEN x, GEN T, GEN pp)
{
  long j, n = lg(x);
  GEN y = cgetg(n,t_MAT);
  if (n == 1) return y;
  for (j=1; j<n; j++) 
    gel(y,j) = FqC_to_FlxC(gel(x,j), T, pp);
  return y;
}

/*******************************************************************/
/*                                                                 */
/*                       n-th ROOT in Fq                           */
/*                                                                 */
/*******************************************************************/
/*NO clean malloc*/
static GEN fflgen(GEN l, long e, GEN r, GEN T ,GEN p, GEN *zeta)
{
  pari_sp av1 = avma;
  GEN z, m, m1;
  const long pp = is_bigint(p)? VERYBIGINT: itos(p);
  long x=varn(T),k,u,v,i;

  for (k=0; ; k++)
  {
    z = (degpol(T)==1)? pol_1[x]: pol_x[x];
    u = k/pp; v=2; /* FpX_Fp_add modify y */
    z = gaddgs(z, k%pp);
    while(u)
    {
      z = ZX_add(z, monomial(utoipos(u%pp),v,x));
      u /= pp; v++;
    }
    if ( DEBUGLEVEL>=6 ) fprintferr("FF l-Gen:next %Z\n",z);
    m1 = m = FpXQ_pow(z,r,T,p);
    if (gcmp1(m)) { avma = av1; continue; }

    for (i=1; i<e; i++)
      if (gcmp1(m = FpXQ_pow(m,l,T,p))) break;
    if (i==e) break;
    avma = av1;
  }
  *zeta = m; return m1;
}
/* Solve x^l = a mod (p,T)
 * l must be prime
 * q = p^degpol(T)-1 = (l^e)*r, with e>=1 and pgcd(r,l)=1
 * m = y^(q/l)
 * y not an l-th power [ m != 1 ] */
GEN
FpXQ_sqrtl(GEN a, GEN l, GEN T ,GEN p , GEN q, long e, GEN r, GEN y, GEN m)
{
  pari_sp av = avma, lim;
  long i,k;
  GEN p1,p2,u1,u2,v,w,z;

  if (gcmp1(a)) return gcopy(a);

  (void)bezout(r,l,&u1,&u2); /* result is 1 */
  v = FpXQ_pow(a,u2,T,p);
  w = FpXQ_pow(a, modii(mulii(negi(u1),r),q), T,p);
  lim = stack_lim(av,1);
  while (!gcmp1(w))
  {
    k = 0;
    p1 = w;
    do { /* if p is not prime, loop will not end */
      z = p1;
      p1 = FpXQ_pow(p1,l,T,p);
      k++;
    } while (!gcmp1(p1));
    if (k==e) { avma=av; return NULL; }
    p2 = FpXQ_mul(z,m,T,p);
    for (i=1; !gcmp1(p2); i++) p2 = FpXQ_mul(p2,m,T,p);/*TODO: BS/GS instead*/
    p1= FpXQ_pow(y, modii(mulsi(i,powiu(l,e-k-1)), q), T,p);
    m = FpXQ_pow(m,utoipos(i),T,p);
    e = k;
    v = FpXQ_mul(p1,v,T,p);
    y = FpXQ_pow(p1,l,T,p);
    w = FpXQ_mul(y,w,T,p);
    if (low_stack(lim, stack_lim(av,1)))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"FpXQ_sqrtl");
      gerepileall(av,4, &y,&v,&w,&m);
    }
  }
  return gerepilecopy(av, v);
}
/* Solve x^n = a mod p: n integer, a in Fp[X]/(T) [ p prime, T irred. mod p ]
 *
 * 1) if no solution, return NULL and (if zetan != NULL) set zetan to gen_0.
 *
 * 2) If there is a solution, there are exactly  m=gcd(p-1,n) of them.
 * If zetan != NULL, it is set to a primitive mth root of unity so that the set
 * of solutions is {x*zetan^k;k=0 to m-1}
 *
 * If a = 0, return 0 and (if zetan != NULL) set zetan = gen_1 */
GEN FpXQ_sqrtn(GEN a, GEN n, GEN T, GEN p, GEN *zetan)
{
  pari_sp ltop=avma, av1, lim;
  long i,j,e;
  GEN m,u1,u2;
  GEN q,r,zeta,y,l,z;

  if (typ(a) != t_POL || typ(n) != t_INT || typ(T) != t_POL || typ(p)!=t_INT)
    pari_err(typeer,"FpXQ_sqrtn");
  if (!degpol(T)) pari_err(constpoler,"FpXQ_sqrtn");
  if (!signe(n)) pari_err(talker,"1/0 exponent in FpXQ_sqrtn");
  if (gcmp1(n)) {if (zetan) *zetan=gen_1;return gcopy(a);}
  if (gcmp0(a)) {if (zetan) *zetan=gen_1;return gen_0;}

  q = addsi(-1, powiu(p,degpol(T)));
  m = bezout(n,q,&u1,&u2);
  if (!equalii(m,n)) a = FpXQ_pow(a, modii(u1,q), T,p);
  if (zetan) z = pol_1[varn(T)];
  lim = stack_lim(ltop,1);
  if (!gcmp1(m))
  {
    m = Z_factor(m); av1 = avma;
    for (i = lg(m[1])-1; i; i--)
    {
      l = gcoeff(m,i,1);
      j = itos(gcoeff(m,i,2));
      e = Z_pvalrem(q,l,&r);
      if(DEBUGLEVEL>=6) (void)timer2();
      y = fflgen(l,e,r,T,p,&zeta);
      if(DEBUGLEVEL>=6) msgtimer("fflgen");
      if (zetan) z = FpXQ_mul(z, FpXQ_pow(y,powiu(l,e-j),T,p), T,p);
      for (; j; j--)
      {
	a = FpXQ_sqrtl(a,l,T,p,q,e,r,y,zeta);
	if (!a) {avma=ltop; return NULL;}
      }
      if (low_stack(lim, stack_lim(ltop,1)))
      { /* n can have lots of prime factors */
	if(DEBUGMEM>1) pari_warn(warnmem,"FpXQ_sqrtn");
        gerepileall(av1,zetan? 2: 1, &a,&z);
      }
    }
  }
  if (zetan)
  {
    *zetan = z;
    gerepileall(ltop,2,&a,zetan);
  }
  else
    a = gerepileupto(ltop, a);
  return a;
}
/*******************************************************************/
/*  Isomorphisms between finite fields                             */
/*******************************************************************/
GEN
FpXQ_matrix_pow(GEN y, long n, long m, GEN P, GEN l)
{
  return RgXV_to_RgM(FpXQ_powers(y,m-1,P,l),n);
}

GEN
Flxq_matrix_pow(GEN y, long n, long m, GEN P, ulong l)
{
  return FlxV_to_Flm(Flxq_powers(y,m-1,P,l),n);
}
/* compute the reciprocical isomorphism of S mod T,p, i.e. V such that
   V(S)=X  mod T,p*/
GEN
FpXQ_ffisom_inv(GEN S,GEN T, GEN p)
{
  pari_sp ltop = avma;
  long n = degpol(T);
  GEN V, M = FpXQ_matrix_pow(S,n,n,T,p);
  V = FpM_invimage(M, col_ei(n, 2), p);
  return gerepileupto(ltop, gtopolyrev(V, varn(T)));
}

/* Let M the matrix of the x^p Frobenius automorphism.
 * Compute x^(p^i) for i=0..r */
static GEN
FpM_Frobenius(GEN M, long r, GEN p, long v)
{
  GEN W, V = cgetg(r+2,t_VEC);
  long i;
  gel(V,1) = pol_x[v]; if (!r) return V;
  gel(V,2) = RgV_to_RgX(gel(M,2),v);
  W = gel(M,2);
  for (i = 3; i <= r+1; ++i)
  {
    W = FpM_FpC_mul(M,W,p);
    gel(V,i) = RgV_to_RgX(W,v);
  }
  return V;
}

/* Let M the matrix of the x^p Frobenius automorphism.
 * Compute x^(p^i) for i=0..r */
static GEN
Flm_Frobenius(GEN M, long r, ulong p, long v)
{
  GEN W, V = cgetg(r+2,t_VEC);
  long i;
  gel(V,1) = polx_Flx(v); if (!r) return V;
  gel(V,2) = Flv_to_Flx(gel(M,2),v);
  W = gel(M,2);
  for (i = 3; i <= r+1; ++i)
  {
    W = Flm_Flc_mul(M,W,p);
    gel(V,i) = Flv_to_Flx(W,v);
  }
  return V;
}

/* Let P a polynomial != 0 and M the matrix of the x^p Frobenius automorphism in
 * FFp[X]/T. Compute P(M)
 * V=FpX_Frobenius(M, p, degpol(P), v)
 * not stack clean
 */

static GEN
FpXQV_FpX_Frobenius(GEN V, GEN P, GEN T, GEN p)
{
  pari_sp btop;
  long i;
  long l=degpol(T);
  long v=varn(T);
  GEN M,W,Mi;
  GEN *gptr[2];
  long lV=lg(V);
  GEN  PV=RgX_to_RgV(P, lgpol(P));
  M=cgetg(l+1,t_VEC);
  gel(M,1) = scalarpol(poleval(P,gen_1),v);
  gel(M,2) = FpXV_FpC_mul(V,PV,p);
  btop=avma;
  gptr[0]=&Mi;
  gptr[1]=&W;
  W=shallowcopy(V);
  for(i=3;i<=l;i++)
  {
    long j;
    pari_sp bbot;
    GEN W2=cgetg(lV,t_VEC);
    for(j=1;j<lV;j++)
      gel(W2,j) = FpXQ_mul(gel(W,j),gel(V,j),T,p);
    bbot=avma;
    Mi=FpXV_FpC_mul(W2,PV,p);
    W=gcopy(W2);
    gerepilemanysp(btop,bbot,gptr,2);
    btop=(pari_sp)W;
    gel(M,i) = Mi;
  }
  return RgXV_to_RgM(M,l);
}

static GEN
FlxqV_Flx_Frobenius(GEN V, GEN P, GEN T, ulong p)
{
  pari_sp btop;
  long i;
  long l=degpol(T);
  long v=varn(T);
  GEN M,W,Mi;
  GEN PV=Flx_to_Flv(P, lgpol(P));
  GEN *gptr[2];
  long lV=lg(V);
  M=cgetg(l+1,t_VEC);
  gel(M,1) = Fl_to_Flx(Flx_eval(P,1,p),v);
  gel(M,2) = FlxV_Flc_mul(V,PV,p);
  btop=avma;
  gptr[0]=&Mi;
  gptr[1]=&W;
  W=gcopy(V);
  for(i=3;i<=l;i++)
  {
    long j;
    pari_sp bbot;
    GEN W2=cgetg(lV,t_VEC);
    for(j=1;j<lV;j++)
      gel(W2,j) = Flxq_mul(gel(W,j),gel(V,j),T,p);
    bbot=avma;
    Mi=FlxV_Flc_mul(W2,PV,p);
    W=gcopy(W2);
    gerepilemanysp(btop,bbot,gptr,2);
    btop=(pari_sp)W;
    gel(M,i) = Mi;
  }
  return FlxV_to_Flm(M,l);
}

/* Let M the matrix of the Frobenius automorphism of Fp[X]/(T).
 * Compute M^d
 * TODO: use left-right binary (tricky!)
 */
GEN
Flm_Frobenius_pow(GEN M, long d, GEN T, ulong p)
{
  pari_sp ltop=avma;
  long i,l=degpol(T);
  GEN R, W = gel(M,2);
  for (i = 2; i <= d; ++i) W = Flm_Flc_mul(M,W,p);
  R=Flxq_matrix_pow(Flv_to_Flx(W,T[2]),l,l,T,p);
  return gerepileupto(ltop,R);
}

GEN
FpM_Frobenius_pow(GEN M, long d, GEN T, GEN p)
{
  pari_sp ltop=avma;
  long i,l=degpol(T);
  GEN R, W = gel(M,2);
  for (i = 2; i <= d; ++i) W = FpM_FpC_mul(M,W,p);
  R=FpXQ_matrix_pow(RgV_to_RgX(W,varn(T)),l,l,T,p);
  return gerepilecopy(ltop,R);
}

/* Essentially we want to compute
 * FqM_ker(MA-pol_x[MAXVARN],U,l)
 * To avoid use of matrix in Fq we procede as follows:
 * We compute FpM_ker(U(MA),l) and then we recover
 * the eigen value by Galois action, see formula.
 */
static GEN
intersect_ker(GEN P, GEN MA, GEN U, GEN l)
{
  pari_sp ltop=avma;
  long vp=varn(P);
  long vu=varn(U), r=degpol(U);
  long i;
  GEN A, R, ib0;
  if (DEBUGLEVEL>=4) (void)timer2();
  if (lgefint(l)==3)
  {
    ulong p=l[2];
    GEN M, V=Flm_Frobenius(ZM_to_Flm(MA, p), r, p, evalvarn(vu));
    if (DEBUGLEVEL>=4) msgtimer("pol[Frobenius]");
    M=FlxqV_Flx_Frobenius(V, ZX_to_Flx(U, p), ZX_to_Flx(P, p), p);
    A=Flm_to_ZM(Flm_ker(M,p));
  }
  else
  {
    GEN V=FpM_Frobenius(MA,r,l,vu);
    if (DEBUGLEVEL>=4) msgtimer("pol[Frobenius]");
    A=FpM_ker(FpXQV_FpX_Frobenius(V, U, P, l), l);
  }
  if (DEBUGLEVEL>=4) msgtimer("matrix cyclo");
  if (lg(A)!=r+1)
    pari_err(talker,"ZZ_%Z[%Z]/(%Z) is not a field in FpX_ffintersect"
        ,l,pol_x[vp],P);
  A=gerepileupto(ltop,A);
  /*The formula is 
   * a_{r-1}=-\phi(a_0)/b_0
   * a_{i-1}=\phi(a_i)+b_ia_{r-1}  i=r-1 to 1
   * Where a_0=A[1] and b_i=U[i+2]
   */
  ib0=negi(Fp_inv(gel(U,2),l));
  R=cgetg(r+1,t_MAT);
  R[1]=A[1];
  gel(R,r) = FpM_FpC_mul(MA,gmul(gel(A,1),ib0),l);
  for(i=r-1;i>1;i--)
    gel(R,i) = FpC_red(gadd(FpM_FpC_mul(MA,gel(R,i+1),l),
         gmul(gel(U,i+2),gel(R,r))),l);
  R=shallowtrans(R);
  for(i=1;i<lg(R);i++)
    gel(R,i) = RgV_to_RgX(gel(R,i),vu);
  A=gtopolyrev(R,vp);
  return gerepileupto(ltop,A);
}

/* n must divide both the degree of P and Q.  Compute SP and SQ such
  that the subfield of FF_l[X]/(P) generated by SP and the subfield of
  FF_l[X]/(Q) generated by SQ are isomorphic of degree n.  P and Q do
  not need to be of the same variable.  if MA (resp. MB) is not NULL,
  must be the matrix of the Frobenius map in FF_l[X]/(P) (resp.
  FF_l[X]/(Q) ).  */
/* Note on the implementation choice:
 * We assume the prime p is very large
 * so we handle Frobenius as matrices.
 */
void
FpX_ffintersect(GEN P, GEN Q, long n, GEN l,GEN *SP, GEN *SQ, GEN MA, GEN MB)
{
  pari_sp lbot, ltop = avma;
  long vp, vq, np, nq, e;
  ulong pg;
  GEN A, B, Ap, Bp;
  GEN *gptr[2];
  vp = varn(P); np = degpol(P);
  vq = varn(Q); nq = degpol(Q);
  if (np<=0 || nq<=0 || n<=0 || np%n!=0 || nq%n!=0)
    pari_err(talker,"bad degrees in FpX_ffintersect: %d,%d,%d",n,np,nq);
  e = u_pvalrem(n, l, &pg);
  if(!MA) MA = FpXQ_matrix_pow(FpXQ_pow(pol_x[vp],l,P,l),np,np,P,l);
  if(!MB) MB = FpXQ_matrix_pow(FpXQ_pow(pol_x[vq],l,Q,l),nq,nq,Q,l);
  A = Ap = zeropol(vp);
  B = Bp = zeropol(vq);
  if (pg > 1)
  {
    GEN ipg = utoipos(pg);
    if (umodiu(l,pg) == 1)
    /* No need to use relative extension, so don't. (Well, now we don't
     * in the other case either, but this special case is more efficient) */
    {
      GEN L, An, Bn, z;
      z = gener_Fp_local(l, gel(Z_factor(ipg), 1));
      z = Fp_pow(z, diviuexact(subis(l,1), pg), l); /* prim. pg-th root of 1 */
      z = negi(z);
      if (DEBUGLEVEL>=4) (void)timer2();
      A = FpM_ker(gaddmat(z, MA),l);
      if (lg(A)!=2)
	pari_err(talker,"ZZ_%Z[%Z]/(%Z) is not a field in FpX_ffintersect"
	    ,l,pol_x[vp],P);
      A = RgV_to_RgX(gel(A,1),vp);

      B = FpM_ker(gaddmat(z, MB),l);
      if (lg(B)!=2)
	pari_err(talker,"ZZ_%Z[%Z]/(%Z) is not a field in FpX_ffintersect"
	    ,l,pol_x[vq],Q);
      B = RgV_to_RgX(gel(B,1),vq);

      if (DEBUGLEVEL>=4) msgtimer("FpM_ker");
      An = (GEN) FpXQ_pow(A,ipg,P,l)[2];
      Bn = (GEN) FpXQ_pow(B,ipg,Q,l)[2];
      if (!invmod(Bn,l,&z))
        pari_err(talker,"Polynomials not irreducible in FpX_ffintersect");
      z = modii(mulii(An,z),l);
      L = Fp_sqrtn(z,ipg,l,NULL);
      if ( !L )
        pari_err(talker,"Polynomials not irreducible in FpX_ffintersect");
      if (DEBUGLEVEL>=4) msgtimer("Fp_sqrtn");
      B = FpX_Fp_mul(B,L,l);
    }
    else
    {
      GEN L, An, Bn, z, U;
      U = gmael(FpX_factor(cyclo(pg,MAXVARN),l),1,1);
      A = intersect_ker(P, MA, U, l); 
      B = intersect_ker(Q, MB, U, l);
      if (DEBUGLEVEL>=4) (void)timer2();
      An = (GEN) FpXYQQ_pow(A,ipg,U,P,l)[2];
      Bn = (GEN) FpXYQQ_pow(B,ipg,U,Q,l)[2];
      if (DEBUGLEVEL>=4) msgtimer("pows [P,Q]");
      z = FpXQ_inv(Bn,U,l);
      z = FpXQ_mul(An,z,U,l);
      L = FpXQ_sqrtn(z,ipg,U,l,NULL);
      if (DEBUGLEVEL>=4) msgtimer("FpXQ_sqrtn");
      if (!L) pari_err(talker,"Polynomials not irreducible in FpX_ffintersect");
      B = FqX_Fq_mul(B,L,U,l);
      B = gsubst(B,MAXVARN,gen_0);
      A = gsubst(A,MAXVARN,gen_0);
    }
  }
  if (e)
  {
    GEN VP, VQ, Ay, By, lmun = addis(l,-1);
    long i, j;
    MA = gaddmat(gen_m1,MA);
    MB = gaddmat(gen_m1,MB);
    Ay = pol_1[vp];
    By = pol_1[vq];
    VP = col_ei(np, 1);
    VQ = np == nq? VP: col_ei(nq, 1); /* save memory */
    for(j=0;j<e;j++)
    {
      if (j)
      {
	Ay = FpXQ_mul(Ay,FpXQ_pow(Ap,lmun,P,l),P,l);
	for(i=1;i<lg(Ay)-1;i++) VP[i] = Ay[i+1];
	for(;i<=np;i++) gel(VP,i) = gen_0;
      }
      Ap = FpM_invimage(MA,VP,l);
      Ap = RgV_to_RgX(Ap,vp);

      if (j)
      {
	By = FpXQ_mul(By,FpXQ_pow(Bp,lmun,Q,l),Q,l);
	for(i=1;i<lg(By)-1;i++) VQ[i] = By[i+1];
	for(;i<=nq;i++) gel(VQ,i) = gen_0;
      }
      Bp = FpM_invimage(MB,VQ,l);
      Bp = RgV_to_RgX(Bp,vq);
      if (DEBUGLEVEL>=4) msgtimer("FpM_invimage");
    }
  }
  A = ZX_add(A,Ap);
  B = ZX_add(B,Bp);
  lbot = avma;
  *SP = FpX_red(A,l);
  *SQ = FpX_red(B,l);
  gptr[0] = SP; 
  gptr[1] = SQ; gerepilemanysp(ltop,lbot,gptr,2);
}
/* Let l be a prime number, P, Q in ZZ[X].  P and Q are both
 * irreducible modulo l and degree(P) divides degree(Q).  Output a
 * monomorphism between FF_l[X]/(P) and FF_l[X]/(Q) as a polynomial R such
 * that Q | P(R) mod l.  If P and Q have the same degree, it is of course an
 * isomorphism.  */
GEN
FpX_ffisom(GEN P,GEN Q,GEN l)
{
  pari_sp av = avma;
  GEN SP, SQ, R;
  FpX_ffintersect(P,Q,degpol(P),l,&SP,&SQ,NULL,NULL);
  R = FpXQ_ffisom_inv(SP,P,l);
  return gerepileupto(av, FpX_FpXQ_compo(R,SQ,Q,l));
}

/* Let l be a prime number, P a ZX irreducible modulo l, MP the matrix of the
 * Frobenius automorphism of F_l[X]/(P).
 * Factor P over the subfield of F_l[X]/(P) of index d. */
static GEN
FpX_factorgalois(GEN P, GEN l, long d, long w, GEN MP)
{
  pari_sp ltop = avma;
  GEN R, V, Tl, z, M;
  long k, n = degpol(P), m = n/d;
  long v = varn(P);

  /* x - y */
  if (m == 1) return deg1pol_i(gen_1, deg1pol_i(subis(l,1), gen_0, w), v);
  M = FpM_Frobenius_pow(MP,d,P,l);

  Tl = gcopy(P); setvarn(Tl,w);
  V = cgetg(m+1,t_VEC);
  gel(V,1) = pol_x[w];
  z = gel(M,2);
  gel(V,2) = RgV_to_RgX(z,w);
  for(k=3;k<=m;k++)
  {
    z = FpM_FpC_mul(M,z,l);
    gel(V,k) = RgV_to_RgX(z,w);
  }
  R = FqV_roots_to_pol(V,Tl,l,v);
  return gerepileupto(ltop,R);
}
/* same: P is an Flx, MP an Flm */
static GEN
Flx_factorgalois(GEN P, ulong l, long d, long w, GEN MP)
{
  pari_sp ltop = avma;
  GEN R, V, Tl, z, M;
  long k, n = degpol(P), m = n/d;
  long v = varn(P);

  if (m == 1) {
    R = polx_Flx(v);
    gel(R,2) = z = polx_Flx(w); z[3] = l - 1; /* - y */
    gel(R,3) = Fl_to_Flx(1, w);
    return R; /* x - y */
  }
  M = Flm_Frobenius_pow(MP,d,P,l);

  Tl = gcopy(P); setvarn(Tl,w);
  V = cgetg(m+1,t_VEC);
  gel(V,1) = polx_Flx(w);
  z = gel(M,2);
  gel(V,2) = Flv_to_Flx(z,w);
  for(k=3;k<=m;k++)
  {
    z = Flm_Flc_mul(M,z,l);
    gel(V,k) = Flv_to_Flx(z,w);
  }
  R = FlxqV_roots_to_pol(V,Tl,l,v);
  return gerepileupto(ltop,R);
}

/* P,Q irreducible over F_l. Factor P over FF_l[X] / Q  [factors are ordered as
 * a Frobenius cycle] */
GEN
FpX_factorff_irred(GEN P, GEN Q, GEN l)
{
  pari_sp ltop = avma, av;
  GEN SP, SQ, MP, MQ, M, FP, FQ, E, V, IR, res;
  long np = degpol(P), nq = degpol(Q), d = cgcd(np,nq);
  long i, vp = varn(P), vq = varn(Q);

  if (d==1) return mkcolcopy(P);
  if (DEBUGLEVEL>=4) (void)timer2();
  if (lgefint(l)==3)
  {
    ulong p = l[2];
    GEN Px = ZX_to_Flx(P,p), Qx = ZX_to_Flx(Q,p);
    FQ = Flxq_matrix_pow(Flxq_pow(polx_Flx(vq),l,Qx,p),nq,nq,Qx,p);
    av = avma;
    FP = Flxq_matrix_pow(Flxq_pow(polx_Flx(vp),l,Px,p),np,np,Px,p);
    if (DEBUGLEVEL>=4) msgtimer("FpXQ_matrix_pows");
    FpX_ffintersect(P,Q,d,l,&SP,&SQ, Flm_to_ZM(FP), Flm_to_ZM(FQ));

    E = Flx_factorgalois(Px,p,d,vq, FP);
    E = FlxX_to_Flm(E,np);
    MP= Flxq_matrix_pow(ZX_to_Flx(SP,p),np,d,Px,p);
    IR= (GEN)Flm_indexrank(MP,p)[1];
    E = rowpermute(E, IR);
    M = rowpermute(MP,IR);
    M = Flm_inv(M,p);
    MQ= Flxq_matrix_pow(ZX_to_Flx(SQ,p),nq,d,Qx,p);
    M = Flm_mul(MQ,M,p);
    M = Flm_mul(M,E,p);
    if (DEBUGLEVEL>=4) msgtimer("factor_irred_mat");
    M = gerepileupto(av,M);
    V = cgetg(d+1,t_VEC);
    gel(V,1) = M;
    for(i=2;i<=d;i++)
      gel(V,i) = Flm_mul(FQ,gel(V,i-1),p);
    res=cgetg(d+1,t_COL);
    for(i=1;i<=d;i++)
      gel(res,i) = FlxX_to_ZXX(Flm_to_FlxX(gel(V,i),evalvarn(vp),evalvarn(vq)));
  }
  else
  {
    FQ = FpXQ_matrix_pow(FpXQ_pow(pol_x[vq],l,Q,l),nq,nq,Q,l);
    av = avma;
    FP = FpXQ_matrix_pow(FpXQ_pow(pol_x[vp],l,P,l),np,np,P,l);
    if (DEBUGLEVEL>=4) msgtimer("FpXQ_matrix_pows");
    FpX_ffintersect(P,Q,d,l,&SP,&SQ,FP,FQ);

    E = FpX_factorgalois(P,l,d,vq,FP);
    E = RgXX_to_RgM(E,np);
    MP= FpXQ_matrix_pow(SP,np,d,P,l);
    IR= (GEN)FpM_indexrank(MP,l)[1];
    E = rowpermute(E, IR);
    M = rowpermute(MP,IR);
    M = FpM_inv(M,l);
    MQ= FpXQ_matrix_pow(SQ,nq,d,Q,l);
    M = FpM_mul(MQ,M,l);
    M = FpM_mul(M,E,l);
    M = gerepileupto(av,M);
    if (DEBUGLEVEL>=4) msgtimer("factor_irred_mat");
    V = cgetg(d+1,t_VEC);
    gel(V,1) = M;
    for(i=2;i<=d;i++)
      gel(V,i) = FpM_mul(FQ,gel(V,i-1),l);
    res = cgetg(d+1,t_COL);
    for(i=1;i<=d;i++)
      gel(res,i) = RgM_to_RgXX(gel(V,i),vp,vq);
  }
  if (DEBUGLEVEL>=4) msgtimer("factor_irred");
  return gerepilecopy(ltop,res);
}
/*******************************************************************/
static GEN
to_intmod(GEN x, GEN p) { return mkintmod(modii(x, p), p); }

/* z in Z[X], return z * Mod(1,p), normalized*/
GEN
FpX_to_mod(GEN z, GEN p)
{
  long i,l = lg(z);
  GEN x = cgetg(l,t_POL); p = icopy(p);
  for (i=2; i<l; i++) gel(x,i) = to_intmod(gel(z,i), p);
  x[1] = z[1]; return normalizepol_i(x,l);
}

/* z in Z^n, return z * Mod(1,p), normalized*/
GEN
FpV_to_mod(GEN z, GEN p)
{
  long i,l = lg(z);
  GEN x = cgetg(l, t_VEC); p = icopy(p);
  for (i=1; i<l; i++) gel(x,i) = to_intmod(gel(z,i), p);
  return x;
}
/* z in Z^n, return z * Mod(1,p), normalized*/
GEN
FpC_to_mod(GEN z, GEN p)
{
  long i,l = lg(z);
  GEN x = cgetg(l, t_COL); p = icopy(p);
  for (i=1; i<l; i++) gel(x,i) = to_intmod(gel(z,i), p);
  return x;
}
/* z in Mat m,n(Z), return z * Mod(1,p), normalized*/
GEN
FpM_to_mod(GEN z, GEN p)
{
  long i,j,l = lg(z), m = lg(gel(z,1));
  GEN  x = cgetg(l,t_MAT), y, zi;
  p = icopy(p);
  for (i=1; i<l; i++)
  {
    gel(x,i) = cgetg(m,t_COL);
    y = gel(x,i); zi= gel(z,i);
    for (j=1; j<m; j++) gel(y,j) = to_intmod(gel(zi,j), p);
  }
  return x;
}
/* z in Z[X], return lift(z * Mod(1,p)), normalized*/
GEN
FpX_red(GEN z, GEN p)
{
  long i, l = lg(z); 
  GEN x = cgetg(l, t_POL);
  for (i=2; i<l; i++) gel(x,i) = modii(gel(z,i),p);
  x[1] = z[1]; return FpX_renormalize(x,l);
}

GEN
FpXV_red(GEN z, GEN p)
{
  long i,l = lg(z);
  GEN x = cgetg(l, t_VEC);
  for (i=1; i<l; i++) gel(x,i) = FpX_red(gel(z,i), p);
  return x;
}

/* z in Z^n, return lift(Col(z) * Mod(1,p)) */
GEN
FpC_red(GEN z, GEN p)
{
  long i,l = lg(z);
  GEN x = cgetg(l, t_COL);
  for (i=1; i<l; i++) gel(x,i) = modii(gel(z,i),p);
  return x;
}

/* z in Z^n, return lift(Vec(z) * Mod(1,p)) */
GEN
FpV_red(GEN z, GEN p)
{
  long i,l = lg(z);
  GEN x = cgetg(l, t_VEC);
  for (i=1; i<l; i++) gel(x,i) = modii(gel(z,i),p);
  return x;
}

/* z in Mat m,n(Z), return lift(z * Mod(1,p)) */
GEN
FpM_red(GEN z, GEN p)
{
  long i, l = lg(z);
  GEN x = cgetg(l,t_MAT);
  for (i=1; i<l; i++) gel(x,i) = FpC_red(gel(z,i), p);
  return x;
}

/* no garbage collection, divide by leading coeff, mod p */
GEN
FpX_normalize(GEN z, GEN p)
{
  GEN p1 = leading_term(z);
  if (lg(z) == 2 || gcmp1(p1)) return z;
  return FpX_Fp_mul(z, Fp_inv(p1,p), p);
}

GEN
FqX_normalize(GEN z, GEN T, GEN p)
{
  GEN p1 = leading_term(z);
  if (lg(z) == 2 || gcmp1(p1)) return z;
  if (!T) return FpX_normalize(z,p);
  return FqX_Fq_mul(z, Fq_inv(p1,T,p), T,p);
}

/* z in R[X,Y] representing an elt in R[X,Y] mod T(Y) in Kronecker form,
 * i.e subst(lift(z), x, y^(2deg(z)-1)). Recover the "real" z, with
 * normalized coefficients */
GEN
from_Kronecker(GEN z, GEN T)
{
  long i,j,lx,l = lg(z), N = (degpol(T)<<1) + 1;
  GEN a,x, t = cgetg(N,t_POL);
  t[1] = T[1] & VARNBITS;
  lx = (l-2) / (N-2); x = cgetg(lx+3,t_POL);
  T = gcopy(T);
  for (i=2; i<lx+2; i++)
  {
    a = cgetg(3,t_POLMOD); gel(x,i) = a;
    gel(a,1) = T;
    for (j=2; j<N; j++) t[j] = z[j];
    z += (N-2);
    gel(a,2) = grem(normalizepol_i(t,N), T);
  }
  a = cgetg(3,t_POLMOD); gel(x,i) = a;
  gel(a,1) = T;
  N = (l-2) % (N-2) + 2;
  for (j=2; j<N; j++) t[j] = z[j];
  gel(a,2) = grem(normalizepol_i(t,N), T);
  return normalizepol_i(x, i+1);
}

GEN
to_Kronecker(GEN P, GEN Q)
{
  /* P(X) = sum Mod(.,Q(Y)) * X^i, lift then set X := Y^(2n-1) */
  long i,j,k,l, lx = lg(P), N = (degpol(Q)<<1) + 1, vQ = varn(Q);
  GEN p1, y = cgetg((N-2)*(lx-2) + 2, t_POL);
  for (k=i=2; i<lx; i++)
  {
    p1 = gel(P,i); l = typ(p1);
    if (l == t_POLMOD) { p1 = gel(p1,2); l = typ(p1); }
    if (is_scalar_t(l) || varncmp(varn(p1), vQ) > 0)
    {
      gel(y,k++) = p1; j = 3;
    }
    else
    {
      l = lg(p1);
      for (j=2; j < l; j++) y[k++] = p1[j];
    }
    if (i == lx-1) break;
    for (   ; j < N; j++) gel(y,k++) = gen_0;
  }
  y[1] = Q[1]; setlg(y, k); return y;
}

/*******************************************************************/
/*                                                                 */
/*                          MODULAR GCD                            */
/*                                                                 */
/*******************************************************************/
/*FIXME: Unify the following 3 divrem routines. Treat the case x,y (lifted) in
 * R[X], y non constant. Given: (lifted) [inv(), mul()], (delayed) red() in R */

/* x and y in Z[X].*/
GEN
FpX_divrem(GEN x, GEN y, GEN p, GEN *pr)
{
  long vx, dx, dy, dz, i, j, sx, lr;
  pari_sp av0, av, tetpil;
  GEN z,p1,rem,lead;

  if (!signe(y)) pari_err(gdiver);
  vx = varn(x);
  dy = degpol(y);
  dx = degpol(x);
  if (dx < dy)
  {
    if (pr)
    {
      av0 = avma; x = FpX_red(x, p);
      if (pr == ONLY_DIVIDES) { avma=av0; return signe(x)? NULL: zeropol(vx); }
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
    av0 = avma; z = FpX_normalize(x, p); 
    if (z==x) return gcopy(z);
    else return gerepileupto(av0, z);
  }
  av0 = avma; dz = dx-dy;
  if (OK_ULONG(p))
  { /* assume ab != 0 mod p */
    ulong pp = (ulong)p[2];
    GEN a = ZX_to_Flx(x, pp);
    GEN b = ZX_to_Flx(y, pp);
    z = Flx_divrem(a,b,pp, pr);
    avma = av0; /* HACK: assume pr last on stack, then z */
    z = shallowcopy(z);
    if (pr && pr != ONLY_DIVIDES && pr != ONLY_REM)
    {
      *pr = shallowcopy(*pr);
      *pr = Flx_to_ZX_inplace(*pr);
    }
    return Flx_to_ZX_inplace(z);
  }
  lead = gcmp1(lead)? NULL: gclone(Fp_inv(lead,p));
  avma = av0;
  z=cgetg(dz+3,t_POL); z[1] = x[1];
  x += 2; y += 2; z += 2;

  p1 = gel(x,dx); av = avma;
  gel(z,dz) = lead? gerepileupto(av, modii(mulii(p1,lead), p)): icopy(p1);
  for (i=dx-1; i>=dy; i--)
  {
    av=avma; p1=gel(x,i);
    for (j=i-dy+1; j<=i && j<=dz; j++)
      p1 = subii(p1, mulii(gel(z,j),gel(y,i-j)));
    if (lead) p1 = mulii(p1,lead);
    tetpil=avma; gel(z,i-dy) = gerepile(av,tetpil,modii(p1, p));
  }
  if (!pr) { if (lead) gunclone(lead); return z-2; }

  rem = (GEN)avma; av = (pari_sp)new_chunk(dx+3);
  for (sx=0; ; i--)
  {
    p1 = gel(x,i);
    for (j=0; j<=i && j<=dz; j++)
      p1 = subii(p1, mulii(gel(z,j),gel(y,i-j)));
    tetpil=avma; p1 = modii(p1,p); if (signe(p1)) { sx = 1; break; }
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
      p1 = subii(p1, mulii(gel(z,j),gel(y,i-j)));
    tetpil=avma; gel(rem,i) = gerepile(av,tetpil, modii(p1,p));
  }
  rem -= 2;
  if (lead) gunclone(lead);
  if (!sx) (void)FpX_renormalize(rem, lr);
  if (pr == ONLY_REM) return gerepileupto(av0,rem);
  *pr = rem; return z-2;
}

/* x and y in Z[Y][X]. Assume T irreducible mod p */
GEN
FpXQX_divrem(GEN x, GEN y, GEN T, GEN p, GEN *pr)
{
  long vx, dx, dy, dz, i, j, sx, lr;
  pari_sp av0, av, tetpil;
  GEN z,p1,rem,lead;

  if (!T) return FpX_divrem(x,y,p,pr);
  if (!signe(y)) pari_err(gdiver);
  vx=varn(x); dy=degpol(y); dx=degpol(x);
  if (dx < dy)
  {
    if (pr)
    {
      av0 = avma; x = FpXQX_red(x, T, p);
      if (pr == ONLY_DIVIDES) { avma=av0; return signe(x)? NULL: zeropol(vx); }
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
    av0 = avma; x = FqX_normalize(x, T,p); tetpil = avma;
    return gerepile(av0,tetpil,FpXQX_red(x,T,p));
  }
  av0 = avma; dz = dx-dy;
  if (OK_ULONG(p))
  { /* assume ab != 0 mod p */
    {
      GEN *gptr[2];
      ulong pp = (ulong)p[2];
      long v = varn(T);
      GEN a = ZXX_to_FlxX(x, pp, v);
      GEN b = ZXX_to_FlxX(y, pp, v);
      GEN t = ZX_to_Flx(T, pp);
      z = FlxqX_divrem(a,b,t,pp,pr);
      tetpil=avma;
      z = FlxX_to_ZXX(z); 
      if (pr && pr != ONLY_DIVIDES && pr != ONLY_REM)
        *pr = FlxX_to_ZXX(*pr);
      else return gerepile(av0,tetpil,z);
      gptr[0]=pr; gptr[1]=&z;
      gerepilemanysp(av0,tetpil,gptr,2);
      return z;
    }
  }
  lead = gcmp1(lead)? NULL: gclone(Fq_inv(lead,T,p));
  avma = av0;
  z = cgetg(dz+3,t_POL); z[1] = x[1];
  x += 2; y += 2; z += 2;

  p1 = gel(x,dx); av = avma;
  gel(z,dz) = lead? gerepileupto(av, Fq_mul(p1,lead, T, p)): gcopy(p1);
  for (i=dx-1; i>=dy; i--)
  {
    av=avma; p1=gel(x,i);
    for (j=i-dy+1; j<=i && j<=dz; j++)
      p1 = Fq_sub(p1, Fq_mul(gel(z,j),gel(y,i-j),NULL,p),NULL,p);
    if (lead) p1 = Fq_mul(p1, lead, NULL,p);
    tetpil=avma; gel(z,i-dy) = gerepile(av,tetpil,Fq_red(p1,T,p));
  }
  if (!pr) { if (lead) gunclone(lead); return z-2; }

  rem = (GEN)avma; av = (pari_sp)new_chunk(dx+3);
  for (sx=0; ; i--)
  {
    p1 = gel(x,i);
    for (j=0; j<=i && j<=dz; j++)
      p1 = Fq_sub(p1, Fq_mul(gel(z,j),gel(y,i-j),NULL,p),NULL,p);
    tetpil=avma; p1 = Fq_red(p1, T, p); if (signe(p1)) { sx = 1; break; }
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
      p1 = Fq_sub(p1, Fq_mul(gel(z,j),gel(y,i-j), NULL,p), NULL,p);
    tetpil=avma; gel(rem,i) = gerepile(av,tetpil, Fq_red(p1, T, p));
  }
  rem -= 2;
  if (lead) gunclone(lead);
  if (!sx) (void)FpXQX_renormalize(rem, lr);
  if (pr == ONLY_REM) return gerepileupto(av0,rem);
  *pr = rem; return z-2;
}

/* x and y in Z[X], return lift(gcd(x mod p, y mod p)) */
GEN
FpX_gcd(GEN x, GEN y, GEN p)
{
  GEN a,b,c;
  pari_sp av0, av;

  if (OK_ULONG(p))
  {
    ulong pp=p[2];
    av = avma;
    (void)new_chunk((lg(x) + lg(y)) << 2); /* scratch space */
    a = ZX_to_Flx(x, pp);
    b = ZX_to_Flx(y, pp);
    a = Flx_gcd_i(a,b, pp);
    avma = av; return Flx_to_ZX(a);
  }
  av0=avma;
  a = FpX_red(x, p); av = avma;
  b = FpX_red(y, p);
  while (signe(b))
  {
    av = avma; c = FpX_rem(a,b,p); a=b; b=c;
  }
  avma = av; return gerepileupto(av0, a);
}

/*Return 1 if gcd can be computed
 * else return a factor of p*/

GEN
FpX_gcd_check(GEN x, GEN y, GEN p)
{
  GEN a,b,c;
  pari_sp av=avma;

  a = FpX_red(x, p); 
  b = FpX_red(y, p);
  while (signe(b))
  {
    GEN lead = leading_term(b);
    GEN g = gcdii(lead,p);
    if (!is_pm1(g)) return gerepileupto(av,g);
    c = FpX_rem(a,b,p); a=b; b=c;
  }
  avma = av; return gen_1;
}

/* x and y in Z[X], return lift(gcd(x mod p, y mod p)). Set u and v st
 * ux + vy = gcd (mod p) */
/*TODO: Document the typ() of *ptu and *ptv*/
GEN
FpX_extgcd(GEN x, GEN y, GEN p, GEN *ptu, GEN *ptv)
{
  GEN a,b,q,r,u,v,d,d1,v1;
  pari_sp ltop=avma, lbot;
  GEN *gptr[3]; 
  if (OK_ULONG(p))
  {
    ulong pp=p[2];
    a = ZX_to_Flx(x, pp);
    b = ZX_to_Flx(y, pp);
    d = Flx_extgcd(a,b, pp, &u,&v);
    lbot=avma;
    d=Flx_to_ZX(d);
    u=Flx_to_ZX(u);
    v=Flx_to_ZX(v);
  }
  else
  {
    a = FpX_red(x, p);
    b = FpX_red(y, p);
    d = a; d1 = b; v = gen_0; v1 = gen_1;
    while (signe(d1))
    {
      q = FpX_divrem(d,d1,p, &r);
      v = gadd(v, gneg_i(gmul(q,v1)));
      v = FpX_red(v,p);
      u=v; v=v1; v1=u;
      u=r; d=d1; d1=u;
    }
    u = gadd(d, gneg_i(gmul(b,v)));
    u = FpX_red(u, p);
    lbot = avma;
    u = FpX_div(u,a,p);
    d = gcopy(d);
    v = gcopy(v);
  }
  gptr[0] = &d; gptr[1] = &u; gptr[2] = &v;
  gerepilemanysp(ltop,lbot,gptr,3);
  *ptu = u; *ptv = v; return d;
}

/* x and y in Z[Y][X], return lift(gcd(x mod T,p, y mod T,p)). Set u and v st
 * ux + vy = gcd (mod T,p) */
GEN
FpXQX_extgcd(GEN x, GEN y, GEN T, GEN p, GEN *ptu, GEN *ptv)
{
  GEN a,b,q,r,u,v,d,d1,v1;
  pari_sp ltop=avma, lbot;

  a = FpXQX_red(x, T, p);
  b = FpXQX_red(y, T, p);
  d = a; d1 = b; v = gen_0; v1 = gen_1;
  while (signe(d1))
  {
    q = FpXQX_divrem(d,d1,T,p, &r);
    v = gadd(v, gneg_i(gmul(q,v1)));
    v = FpXQX_red(v,T,p);
    u=v; v=v1; v1=u;
    u=r; d=d1; d1=u;
  }
  u = gadd(d, gneg_i(gmul(b,v)));
  u = FpXQX_red(u,T, p);
  lbot = avma;
  u = FpXQX_divrem(u,a,T,p,NULL);
  d = gcopy(d);
  v = gcopy(v);
  {
    GEN *gptr[3]; gptr[0] = &d; gptr[1] = &u; gptr[2] = &v;
    gerepilemanysp(ltop,lbot,gptr,3);
  }
  *ptu = u; *ptv = v; return d;
}

/*x must be reduced*/
GEN
FpXQ_charpoly(GEN x, GEN T, GEN p)
{
  pari_sp ltop=avma;
  long v=varn(T);
  GEN R;
  T = gcopy(T); setvarn(T, MAXVARN);
  x = gcopy(x); setvarn(x, MAXVARN);
  R = FpY_FpXY_resultant(T, deg1pol_i(gen_1,FpX_neg(x,p),v),p);
  return gerepileupto(ltop,R);
}

GEN 
FpXQ_minpoly(GEN x, GEN T, GEN p)
{
  pari_sp ltop=avma;
  GEN R=FpXQ_charpoly(x, T, p);
  GEN G=FpX_gcd(R,derivpol(R),p);
  G=FpX_normalize(G,p);
  G=FpX_div(R,G,p);
  return gerepileupto(ltop,G);
}

/* return z = a mod q, b mod p (p,q) = 1. qinv = 1/q mod p */
static GEN
Fl_chinese_coprime(GEN a, ulong b, GEN q, ulong p, ulong qinv, GEN pq)
{
  ulong d, amod = umodiu(a, p);
  pari_sp av = avma;
  GEN ax;

  if (b == amod) return NULL;
  d = (b > amod)? b - amod: p - (amod - b); /* (b - a) mod p */
  (void)new_chunk(lgefint(pq)<<1); /* HACK */
  ax = mului(Fl_mul(d,qinv,p), q); /* d mod p, 0 mod q */
  avma = av; return addii(a, ax); /* in ]-q, pq[ assuming a in -]-q,q[ */
}

GEN
ZX_init_CRT(GEN Hp, ulong p, long v)
{
  long i, l = lg(Hp), lim = (long)(p>>1);
  GEN H = cgetg(l, t_POL);
  H[1] = evalsigne(1) | evalvarn(v);
  for (i=2; i<l; i++)
    gel(H,i) = stoi(Fl_center(Hp[i], p, lim));
  return H;
}

/* assume lg(Hp) > 1 */
GEN 
ZM_init_CRT(GEN Hp, ulong p)
{
  long i,j, m = lg(Hp[1]), l = lg(Hp), lim = (long)(p>>1);
  GEN c,cp,H = cgetg(l, t_MAT);
  for (j=1; j<l; j++)
  {
    cp = gel(Hp,j);
    c = cgetg(m, t_COL);
    gel(H,j) = c;
    for (i=1; i<l; i++) gel(c,i) = stoi(Fl_center(cp[i],p, lim));
  }   
  return H;
}

int
Z_incremental_CRT(GEN *H, ulong Hp, GEN q, GEN qp, ulong p)
{
  GEN h, lim = shifti(qp,-1);
  ulong qinv = Fl_inv(umodiu(q,p), p);
  int stable = 1;
  h = Fl_chinese_coprime(*H,Hp,q,p,qinv,qp);
  if (h)
  {
    if (cmpii(h,lim) > 0) h = subii(h,qp);
    *H = h; stable = 0;
  }
  return stable;
}

int
ZX_incremental_CRT(GEN *ptH, GEN Hp, GEN q, GEN qp, ulong p)
{
  GEN H = *ptH, h, lim = shifti(qp,-1);
  ulong qinv = Fl_inv(umodiu(q,p), p);
  long i, l = lg(H), lp = lg(Hp);
  int stable = 1;

  if (l < lp)
  { /* degree increases */
    GEN x = cgetg(lp, t_POL);
    for (i=1; i<l; i++)  x[i] = H[i];
    for (   ; i<lp; i++) gel(x,i) = gen_0;
    *ptH = H = x;
    stable = 0;
  } else if (l > lp)
  { /* degree decreases */
    GEN x = cgetg(l, t_VECSMALL);
    for (i=1; i<lp; i++)  x[i] = Hp[i];
    for (   ; i<l; i++) x[i] = 0;
    Hp = x; lp = l;
  }
  for (i=2; i<lp; i++)
  {
    h = Fl_chinese_coprime(gel(H,i),Hp[i],q,p,qinv,qp);
    if (h)
    {
      if (cmpii(h,lim) > 0) h = subii(h,qp);
      gel(H,i) = h; stable = 0;
    }
  }
  return stable;
}

int
ZM_incremental_CRT(GEN H, GEN Hp, GEN q, GEN qp, ulong p)
{
  GEN h, lim = shifti(qp,-1);
  ulong qinv = Fl_inv(umodiu(q,p), p);
  long i,j, l = lg(H), m = lg(H[1]);
  int stable = 1;
  for (j=1; j<l; j++)
    for (i=1; i<m; i++)
    {
      h = Fl_chinese_coprime(gcoeff(H,i,j), coeff(Hp,i,j),q,p,qinv,qp);
      if (h)
      {
        if (cmpii(h,lim) > 0) h = subii(h,qp);
        gcoeff(H,i,j) = h; stable = 0;
      }
    }
  return stable;
}

/* returns a polynomial in variable v, whose coeffs correspond to the digits
 * of m (in base p) */
GEN
stopoly(ulong m, ulong p, long v)
{
  GEN y = new_chunk(BITS_IN_LONG + 2);
  long l = 2;
  do { ulong q = m/p; gel(y,l++) = utoi(m - q*p); m=q; } while (m);
  y[1] = evalsigne(1) | evalvarn(v); 
  y[0] = evaltyp(t_POL) | evallg(l); return y;
}

GEN
stopoly_gen(GEN m, GEN p, long v)
{
  GEN y = new_chunk(bit_accuracy(lgefint(m))+2);
  long l = 2;
  do { m = dvmdii(m, p, &gel(y,l)); l++; } while (signe(m));
  y[1] = evalsigne(1) | evalvarn(v);
  y[0] = evaltyp(t_POL) | evallg(l); return y;
}

static GEN
muliimod(GEN x, GEN y, GEN p)
{
  return modii(mulii(x,y), p);
}

/* Res(A,B) = Res(B,R) * lc(B)^(a-r) * (-1)^(ab), with R=A%B, a=deg(A) ...*/
GEN
FpX_resultant(GEN a, GEN b, GEN p)
{
  long da,db,dc;
  pari_sp av, lim;
  GEN c,lb, res = gen_1;

  if (!signe(a) || !signe(b)) return gen_0;
  da = degpol(a);
  db = degpol(b);
  if (db > da)
  {
    swapspec(a,b, da,db);
    if (both_odd(da,db)) res = subii(p, res);
  }
  if (!da) return gen_1; /* = res * a[2] ^ db, since 0 <= db <= da = 0 */
  av = avma; lim = stack_lim(av,2);
  while (db)
  {
    lb = gel(b,db+2);
    c = FpX_rem(a,b, p);
    a = b; b = c; dc = degpol(c);
    if (dc < 0) { avma = av; return 0; }

    if (both_odd(da,db)) res = subii(p, res);
    if (!gcmp1(lb)) res = muliimod(res, Fp_powu(lb, da - dc, p), p);
    if (low_stack(lim,stack_lim(av,2)))
    {
      if (DEBUGMEM>1) pari_warn(warnmem,"FpX_resultant (da = %ld)",da);
      gerepileall(av,3, &a,&b,&res);
    }
    da = db; /* = degpol(a) */
    db = dc; /* = degpol(b) */
  }
  res = muliimod(res, Fp_powu(gel(b,2), da, p), p);
  return gerepileuptoint(av, res);
}

/* assuming the PRS finishes on a degree 1 polynomial C0 + C1X, with "generic"
 * degree sequence as given by dglist, set *Ci and return resultant(a,b) */
static ulong
Flx_resultant_all(GEN a, GEN b, long *C0, long *C1, GEN dglist, ulong p)
{
  long da,db,dc,cnt,ind;
  ulong lb, cx = 1, res = 1UL;
  pari_sp av;
  GEN c;

  if (C0) { *C0 = 1; *C1 = 0; }
  if (lgpol(a)==0 || lgpol(b)==0) return 0;
  da = degpol(a);
  db = degpol(b);
  if (db > da)
  {
    swapspec(a,b, da,db);
    if (both_odd(da,db)) res = p-res;
  }
  /* = res * a[2] ^ db, since 0 <= db <= da = 0 */
  if (!da) return 1;
  cnt = ind = 0; av = avma;
  while (db)
  {
    lb = b[db+2];
    c = Flx_rem(a,b, p);
    a = b; b = c; dc = degpol(c);
    if (dc < 0) { avma = av; return 0; }

    ind++;
    if (C0)
    { /* check that Euclidean remainder sequence doesn't degenerate */
      if (dc != dglist[ind]) { avma = av; return 0; }
      /* update resultant */
      if (both_odd(da,db)) res = p-res;
      if (lb != 1)
      {
        ulong t = Fl_pow(lb, da - dc, p);
        res = Fl_mul(res, t, p);
        if (dc) cx = Fl_mul(cx, t, p);
      }
    }
    else
    {
      if (dc > dglist[ind]) dglist[ind] = dc;
    }
    if (++cnt == 4) { cnt = 0; avma = av; }
    da = db; /* = degpol(a) */
    db = dc; /* = degpol(b) */
  }
  if (!C0)
  {
    if (ind+1 > lg(dglist)) setlg(dglist,ind+1);
    return 0;
  }

  if (da == 1) /* last non-constant polynomial has degree 1 */
  {
    *C0 = Fl_mul(cx, a[2], p);
    *C1 = Fl_mul(cx, a[3], p);
    lb = b[2];
  } else lb = Fl_pow(b[2], da, p);
  avma = av; return Fl_mul(res, lb, p);
}

/* u P(X) + v P(-X) */
static GEN
pol_comp(GEN P, GEN u, GEN v)
{
  long i, l = lg(P);
  GEN y = cgetg(l, t_POL);
  for (i=2; i<l; i++)
  {
    GEN t = gel(P,i);
    gel(y,i) = gcmp0(t)? gen_0:
                         (i&1)? gmul(t, gsub(u,v)) /*  odd degree */
                              : gmul(t, gadd(u,v));/* even degree */
  }
  y[1] = P[1]; return normalizepol_i(y,l);
}

GEN
polint_triv(GEN xa, GEN ya)
{
  GEN P = NULL, Q = roots_to_pol(xa,0);
  long i, n = lg(xa);
  pari_sp av = avma, lim = stack_lim(av, 2);
  for (i=1; i<n; i++)
  {
    GEN T, dP, r;
    if (gcmp0(gel(ya,i))) continue;
    T = RgX_div_by_X_x(Q, gel(xa,i), NULL);
    r = poleval(T, gel(xa,i));
    if (i < n-1 && absi_equal(gel(xa,i), gel(xa,i+1)))
    { /* x_i = -x_{i+1} */
      dP = pol_comp(gdiv(T, r), gel(ya,i), gel(ya,i+1));
      i++;
    }
    else
      dP = gdiv(gmul(gel(ya,i), T), r);
    P = P? gadd(P, dP): dP;
    if (low_stack(lim,stack_lim(av,2)))
    {
      if (DEBUGMEM>1) pari_warn(warnmem,"polint_triv2 (i = %ld)",i);
      P = gerepileupto(av, P);
    }
  }
  return P? P: zeropol(0);
}

GEN
FpX_div_by_X_x(GEN a, GEN x, GEN p, GEN *r)
{
  long l = lg(a), i;
  GEN z = cgetg(l-1, t_POL), a0, z0;
  z[1] = evalsigne(1) | evalvarn(0);
  a0 = a + l-1;
  z0 = z + l-2; *z0 = *a0--;
  for (i=l-3; i>1; i--) /* z[i] = (a[i+1] + x*z[i+1]) % p */
  {
    GEN t = addii(gel(a0--,0), muliimod(x, gel(z0--,0), p));
    *z0 = (long)t;
  }
  if (r) *r = addii(gel(a0,0), muliimod(x, gel(z0,0), p));
  return z;
}

GEN
FpV_polint(GEN xa, GEN ya, GEN p)
{
  GEN inv,T,dP, P = NULL, Q = FpV_roots_to_pol(xa, p, 0);
  long i, n = lg(xa);
  pari_sp av, lim;
  av = avma; lim = stack_lim(av,2);
  for (i=1; i<n; i++)
  {
    if (!signe(ya[i])) continue;
    T = FpX_div_by_X_x(Q, gel(xa,i), p, NULL);
    inv = Fp_inv(FpX_eval(T,gel(xa,i), p), p);
    if (i < n-1 && equalii(addii(gel(xa,i), gel(xa,i+1)), p))
    {
      dP = pol_comp(T, muliimod(gel(ya,i),  inv,p),
                       muliimod(gel(ya,i+1),inv,p));
      i++; /* x_i = -x_{i+1} */
    }
    else
      dP = FpX_Fp_mul(T, muliimod(gel(ya,i),inv,p), p);
    P = P? FpX_add(P, dP, p): dP;
    if (low_stack(lim, stack_lim(av,2)))
    {
      if (DEBUGMEM>1) pari_warn(warnmem,"FpV_polint");
      if (!P) avma = av; else P = gerepileupto(av, P);
    }
  }
  return P? P: zeropol(0);
}

static void
FlV_polint_all(GEN xa, GEN ya, GEN C0, GEN C1, ulong p)
{
  GEN T,Q = Flv_roots_to_pol(xa, p, 0);
  GEN dP  = NULL,  P = NULL;
  GEN dP0 = NULL, P0= NULL;
  GEN dP1 = NULL, P1= NULL;
  long i, n = lg(xa);
  ulong inv;
  for (i=1; i<n; i++)
  {
    T = Flx_div_by_X_x(Q, xa[i], p, NULL);
    inv = Fl_inv(Flx_eval(T,xa[i], p), p);

    if (ya[i])
    {
      dP = Flx_Fl_mul(T, Fl_mul(ya[i],inv,p), p);
      P = P ? Flx_add(P , dP , p): dP;
    }
    if (C0[i])
    {
      dP0= Flx_Fl_mul(T, Fl_mul(C0[i],inv,p), p);
      P0= P0? Flx_add(P0, dP0, p): dP0;
    }
    if (C1[i])
    {
      dP1= Flx_Fl_mul(T, Fl_mul(C1[i],inv,p), p);
      P1= P1? Flx_add(P1, dP1, p): dP1;
    }
  }
  gel(ya,1) = (P ? P : zero_Flx(0));
  gel(C0,1) = (P0? P0: zero_Flx(0));
  gel(C1,1) = (P1? P1: zero_Flx(0));
}

/* b a vector of polynomials representing B in Fp[X][Y], evaluate at X = x,
 * Return 0 in case of degree drop. */
static GEN
FlxV_eval(GEN b, ulong x, ulong p)
{
  GEN z;
  long i, lb = lg(b);
  ulong leadz = Flx_eval(leading_term(b), x, p);
  long vs=mael(b,2,1);
  if (!leadz) return zero_Flx(vs);

  z = cgetg(lb, t_VECSMALL); z[1] = vs;
  for (i=2; i<lb-1; i++) z[i] = Flx_eval(gel(b,i), x, p);
  z[i] = leadz; return z;
}

/* as above, but don't care about degree drop */
static GEN
FlxV_eval_gen(GEN b, ulong x, ulong p, long *drop)
{
  GEN z;
  long i, lb = lg(b);
  z = cgetg(lb,t_VECSMALL); z[1]=mael(b,2,1);

  for (i=2; i<lb; i++) z[i] = Flx_eval(gel(b,i), x, p);
  z = Flx_renormalize(z, lb);
  *drop = lb - lg(z);
  return z;
}

static GEN
vec_FpX_eval_gen(GEN b, GEN x, GEN p, long *drop)
{
  GEN z;
  long i, lb = lg(b);
  z = cgetg(lb, t_POL);
  z[1] = b[1];
  for (i=2; i<lb; i++)
    gel(z,i) = FpX_eval(gel(b,i), x, p);
  z = FpX_renormalize(z, lb);
  *drop = lb - lg(z);
  return z;
}

/* Interpolate at roots of 1 and use Hadamard bound for univariate resultant:
 *   bound = N_2(A)^degpol B N_2(B)^degpol(A),  where
 *     N_2(A) = sqrt(sum (N_1(Ai))^2)
 * Return e such that Res(A, B) < 2^e */
ulong
ZY_ZXY_ResBound(GEN A, GEN B, GEN dB)
{
  pari_sp av = avma;
  GEN a = gen_0, b = gen_0;
  long i , lA = lg(A), lB = lg(B);
  double loga, logb;
  for (i=2; i<lA; i++) a = addii(a, sqri(gel(A,i)));
  for (i=2; i<lB; i++)
  {
    GEN t = gel(B,i);
    if (typ(t) == t_POL) t = gnorml1(t, 0);
    b = addii(b, sqri(t));
  }
  loga = dbllog2(a);
  logb = dbllog2(b); if (dB) logb -= 2 * dbllog2(dB);
  i = (long)((degpol(B) * loga + degpol(A) * logb) / 2);
  avma = av; return (i <= 0)? 1: 1 + (ulong)i;
}

/* return Res(a(Y), b(n,Y)) over Fp. la = leading_term(a) [for efficiency] */
static ulong
FlX_eval_resultant(GEN a, GEN b, ulong n, ulong p, ulong la)
{
  long drop;
  GEN ev = FlxV_eval_gen(b, n, p, &drop);
  ulong r = Flx_resultant(a, ev, p);
  if (drop && la != 1) r = Fl_mul(r, Fl_pow(la, drop,p),p);
  return r;
}
static GEN
FpX_eval_resultant(GEN a, GEN b, GEN n, GEN p, GEN la)
{
  long drop;
  GEN ev = vec_FpX_eval_gen(b, n, p, &drop);
  GEN r = FpX_resultant(a, ev, p);
  if (drop && !gcmp1(la)) r = muliimod(r, Fp_powu(la, drop,p),p);
  return r;
}

/* assume dres := deg(Res_Y(a,b), X) <= deg(a,Y) * deg(b,X) < p */
static GEN
Fly_Flxy_resultant_polint(GEN a, GEN b, ulong p, ulong dres)
{
  ulong i, n, la = (ulong)leading_term(a);
  GEN  x = cgetg(dres+2, t_VECSMALL);
  GEN  y = cgetg(dres+2, t_VECSMALL);
 /* Evaluate at dres+ 1 points: 0 (if dres even) and +/- n, so that P_n(X) =
  * P_{-n}(-X), where P_i is Lagrange polynomial: P_i(j) = delta_{i,j} */
  for (i=0,n = 1; i < dres; n++)
  {
    x[++i] = n;   y[i] = FlX_eval_resultant(a,b, x[i], p,la);
    x[++i] = p-n; y[i] = FlX_eval_resultant(a,b, x[i], p,la);
  }
  if (i == dres)
  {
    x[++i] = 0;   y[i] = FlX_eval_resultant(a,b, x[i], p,la);
  }
  return Flv_polint(x,y, p, evalvarn(varn(b)));
}

static GEN
FlxX_pseudorem(GEN x, GEN y, ulong p)
{
  long vx = varn(x), dx, dy, dz, i, lx, dp;
  pari_sp av = avma, av2, lim;

  if (!signe(y)) pari_err(gdiver);
  (void)new_chunk(2);
  dx=degpol(x); x = revpol(x);
  dy=degpol(y); y = revpol(y); dz=dx-dy; dp = dz+1;
  av2 = avma; lim = stack_lim(av2,1);
  for (;;)
  {
    gel(x,0) = Flx_neg(gel(x,0), p); dp--;
    for (i=1; i<=dy; i++)
      gel(x,i) = Flx_add( Flx_mul(gel(y,0), gel(x,i), p),
                              Flx_mul(gel(x,0), gel(y,i), p), p );
    for (   ; i<=dx; i++)
      gel(x,i) = Flx_mul(gel(y,0), gel(x,i), p);
    do { x++; dx--; } while (dx >= 0 && lg(gel(x,0))==2);
    if (dx < dy) break;
    if (low_stack(lim,stack_lim(av2,1)))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"pseudorem dx = %ld >= %ld",dx,dy);
      gerepilecoeffs(av2,x,dx+1);
    }
  }
  if (dx < 0) return zero_Flx(0);
  lx = dx+3; x -= 2;
  x[0]=evaltyp(t_POL) | evallg(lx);
  x[1]=evalsigne(1) | evalvarn(vx);
  x = revpol(x) - 2;
  if (dp)
  { /* multiply by y[0]^dp   [beware dummy vars from FpY_FpXY_resultant] */
    GEN t = Flx_pow(gel(y,0), dp, p);
    for (i=2; i<lx; i++)
      gel(x,i) = Flx_mul(gel(x,i), t, p);
  }
  return gerepilecopy(av, x);
}

static GEN
FlxX_Flx_div(GEN x, GEN y, ulong p)
{
  long i, l;
  GEN z;
  if (degpol(y) == 0)
  {
    ulong t = (ulong)y[2];
    if (t == 1) return x;
    t = Fl_inv(t, p);
    l = lg(x); z = cgetg(l, t_POL); z[1] = x[1];
    for (i=2; i<l; i++) gel(z,i) = Flx_Fl_mul(gel(x,i),t,p);
  }
  else
  {
    l = lg(x); z = cgetg(l, t_POL); z[1] = x[1];
    for (i=2; i<l; i++) gel(z,i) = Flx_div(gel(x,i),y,p);
  }
  return z;
}

static GEN
FlxX_subres(GEN u, GEN v, ulong p)
{
  pari_sp av = avma, av2, lim;
  long degq,dx,dy,du,dv,dr,signh;
  GEN z,g,h,r,p1;

  dx=degpol(u); dy=degpol(v); signh=1;
  if (dx < dy)
  {
    swap(u,v); lswap(dx,dy);
    if (both_odd(dx, dy)) signh = -signh;
  }
  if (dy < 0) return gen_0;
  if (dy==0) return gerepileupto(av, Flx_pow(gel(v,2),dx,p));

  g = h = Fl_to_Flx(1,0); av2 = avma; lim = stack_lim(av2,1);
  for(;;)
  {
    r = FlxX_pseudorem(u,v,p); dr = lg(r);
    if (dr == 2) { avma = av; return gen_0; }
    du = degpol(u); dv = degpol(v); degq = du-dv;
    u = v; p1 = g; g = leading_term(u);
    switch(degq)
    {
      case 0: break;
      case 1:
        p1 = Flx_mul(h,p1, p); h = g; break;
      default:
        p1 = Flx_mul(Flx_pow(h,degq,p), p1, p);
        h = Flx_div(Flx_pow(g,degq,p), Flx_pow(h,degq-1,p), p);
    }
    if (both_odd(du,dv)) signh = -signh;
    v = FlxX_Flx_div(r, p1, p);
    if (dr==3) break;
    if (low_stack(lim,stack_lim(av2,1)))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"subresall, dr = %ld",dr);
      gerepileall(av2,4, &u, &v, &g, &h);
    }
  }
  z = gel(v,2);
  if (dv > 1) z = Flx_div(Flx_pow(z,dv,p), Flx_pow(h,dv-1,p), p);
  if (signh < 0) z = Flx_neg(z,p);
  return gerepileupto(av, z);
}

/* return a t_POL (in variable v) whose coeffs are the coeffs of b,
 * in variable v. This is an incorrect PARI object if initially varn(b) << v.
 * We could return a vector of coeffs, but it is convenient to have degpol()
 * and friends available. Even in that case, it will behave nicely with all
 * functions treating a polynomial as a vector of coeffs (eg poleval). 
 * FOR INTERNAL USE! */
GEN
swap_vars(GEN b0, long v)
{
  long i, n = poldegree(b0, v);
  GEN b, x;
  if (n < 0) return zeropol(v);
  b = cgetg(n+3, t_POL); x = b + 2;
  b[1] = evalsigne(1) | evalvarn(v);
  for (i=0; i<=n; i++) gel(x,i) = polcoeff_i(b0, i, v);
  return b;
}

/* assume varn(b) << varn(a) */
GEN
FpY_FpXY_resultant(GEN a, GEN b0, GEN p)
{
  long i,n,dres, vX = varn(b0), vY = varn(a);
  GEN la,x,y,b = swap_vars(b0, vY);
 
  dres = degpol(a)*degpol(b0);
  if (OK_ULONG(p))
  {
    ulong pp = (ulong)p[2];
    b = ZXX_to_FlxX(b, pp, vX);
    if ((ulong)dres >= pp)
    {
      a = ZXX_to_FlxX(a, pp, vX);
      x = FlxX_subres(a, b, pp);
    }
    else
    {
      a = ZX_to_Flx(a, pp);
      x = Fly_Flxy_resultant_polint(a, b, pp, (ulong)dres);
      setvarn(x, vX);
    }
    return Flx_to_ZX(x);
  }
 
  la = leading_term(a);
  x = cgetg(dres+2, t_VEC);
  y = cgetg(dres+2, t_VEC);
 /* Evaluate at dres+ 1 points: 0 (if dres even) and +/- n, so that P_n(X) =
  * P_{-n}(-X), where P_i is Lagrange polynomial: P_i(j) = delta_{i,j} */
  for (i=0,n = 1; i < dres; n++)
  {
    gel(x,++i) = utoipos(n); gel(y,i) = FpX_eval_resultant(a,b,gel(x,i),p,la);
    gel(x,++i) = subis(p,n); gel(y,i) = FpX_eval_resultant(a,b,gel(x,i),p,la);
  }
  if (i == dres)
  {
    gel(x,++i) = gen_0;        gel(y,i) = FpX_eval_resultant(a,b, gel(x,i), p,la);
  }
  x = FpV_polint(x,y, p);
  setvarn(x, vX); return x;
}

/* check that theta(maxprime) - theta(27448) >= 2^bound */
/* NB: theta(27449) ~ 27225.387, theta(x) > 0.98 x for x>7481
 * (Schoenfeld, 1976 for x > 1155901 + direct calculations) */
static void
check_theta(ulong bound) {
  maxprime_check( (ulong)ceil((bound * LOG2 + 27225.388) / 0.98) );
}
/* 27449 = prime(3000) */
byteptr
init_modular(ulong *p) { *p = 27449; return diffptr + 3000; }

/* 0, 1, -1, 2, -2, ... */
#define next_lambda(a) (a>0 ? -a : 1-a)

/* Assume A in Z[Y], B in Q[Y][X], and Res_Y(A, B) in Z[X].
 * If lambda = NULL, return Res_Y(A,B).
 * Otherwise, find a small lambda (start from *lambda, use the sequence above)
 * such that R(X) = Res_Y(A, B(X + lambda Y)) is squarefree, reset *lambda to
 * the chosen value and return R
 *
 * If LERS is non-NULL, set it to the Last non-constant polynomial in the
 * Euclidean Remainder Sequence */
GEN
ZY_ZXY_resultant_all(GEN A, GEN B0, long *lambda, GEN *LERS)
{
  int checksqfree = lambda? 1: 0, delvar = 0, stable;
  ulong bound, p, dp;
  pari_sp av = avma, av2 = 0, lim;
  long i,n, lb, degA = degpol(A), dres = degA*degpol(B0);
  long vX = varn(B0), vY = varn(A); /* assume vX << vY */
  GEN x, y, dglist, dB, B, q, a, b, ev, H, H0, H1, Hp, H0p, H1p, C0, C1, L;
  byteptr d = init_modular(&p);

  dglist = Hp = H0p = H1p = C0 = C1 = NULL; /* gcc -Wall */
  if (LERS)
  {
    if (!lambda) pari_err(talker,"ZY_ZXY_resultant_all: LERS needs lambda");
    C0 = cgetg(dres+2, t_VECSMALL);
    C1 = cgetg(dres+2, t_VECSMALL);
    dglist = cgetg(dres+1, t_VECSMALL);
  }
  x = cgetg(dres+2, t_VECSMALL);
  y = cgetg(dres+2, t_VECSMALL);
  if (vY == MAXVARN)
  {
    vY = fetch_var(); delvar = 1;
    B0 = gsubst(B0, MAXVARN, pol_x[vY]);
    A = shallowcopy(A); setvarn(A, vY);
  }
  L = pol_x[MAXVARN];
  B0 = Q_remove_denom(B0, &dB);
  lim = stack_lim(av,2);

INIT:
  if (av2) { avma = av2; *lambda = next_lambda(*lambda); } 
  if (lambda)
  {
    L = gadd(pol_x[MAXVARN], gmulsg(*lambda, pol_x[vY]));
    if (DEBUGLEVEL>4) fprintferr("Trying lambda = %ld\n",*lambda);
  }
  B = poleval(B0, L); av2 = avma;

  if (degA <= 3)
  { /* sub-resultant faster for small degrees */
    if (LERS)
    {
      H = subresall(A,B,&q);
      if (typ(q) != t_POL || degpol(q)!=1 || !ZX_is_squarefree(H)) goto INIT;
      H0 = gel(q,2); if (typ(H0) == t_POL) setvarn(H0,vX);
      H1 = gel(q,3); if (typ(H1) == t_POL) setvarn(H1,vX);
    }
    else
    {
      H = subres(A,B);
      if (checksqfree && !ZX_is_squarefree(H)) goto INIT;
    }
    goto END;
  }

  /* make sure p large enough */
  while (p < (ulong)(dres<<1)) NEXT_PRIME_VIADIFF(p,d);

  H = H0 = H1 = NULL;
  lb = lg(B); 
  bound = ZY_ZXY_ResBound(A, B, dB);
  if (DEBUGLEVEL>4) fprintferr("bound for resultant coeffs: 2^%ld\n",bound);
  check_theta(bound);

  dp = 1;
  for(;;)
  {
    NEXT_PRIME_VIADIFF_CHECK(p,d);
    if (dB) { dp = smodis(dB, p); if (!dp) continue; }

    a = ZX_to_Flx(A, p);
    b = ZXX_to_FlxX(B, p, varn(A));
    if (LERS)
    {
      if (!b[lb-1] || degpol(a) < degA) continue; /* p | lc(A)lc(B) */
      if (checksqfree)
      { /* find degree list for generic Euclidean Remainder Sequence */
        long goal = min(degpol(a), degpol(b)); /* longest possible */
        for (n=1; n <= goal; n++) dglist[n] = 0;
        setlg(dglist, 1);
        for (n=0; n <= dres; n++)
        {
          ev = FlxV_eval(b, n, p);
          (void)Flx_resultant_all(a, ev, NULL, NULL, dglist, p);
          if (lg(dglist)-1 == goal) break;
        }
        /* last pol in ERS has degree > 1 ? */
        goal = lg(dglist)-1;
        if (degpol(B) == 1) { if (!goal) goto INIT; }
        else
        {
          if (goal <= 1) goto INIT;
          if (dglist[goal] != 0 || dglist[goal-1] != 1) goto INIT;
        }
        if (DEBUGLEVEL>4)
          fprintferr("Degree list for ERS (trials: %ld) = %Z\n",n+1,dglist);
      }

      for (i=0,n = 0; i <= dres; n++)
      {
        ev = FlxV_eval(b, n, p);
        x[++i] = n; y[i] = Flx_resultant_all(a, ev, C0+i, C1+i, dglist, p);
        if (!C1[i]) i--; /* C1(i) = 0. No way to recover C0(i) */
      }
      FlV_polint_all(x,y,C0,C1, p);
      Hp = gel(y,1);
      H0p= gel(C0,1);
      H1p= gel(C1,1);
    }
    else
      Hp = Fly_Flxy_resultant_polint(a, b, p, (ulong)dres);
    if (!H && degpol(Hp) != dres) continue;
    if (dp != 1) Hp = Flx_Fl_mul(Hp, Fl_pow(Fl_inv(dp,p), degA, p), p);
    if (checksqfree) {
      if (!Flx_is_squarefree(Hp, p)) goto INIT;
      if (DEBUGLEVEL>4) fprintferr("Final lambda = %ld\n",*lambda);
      checksqfree = 0;
    }

    if (!H)
    { /* initialize */
      q = utoipos(p); stable = 0;
      H = ZX_init_CRT(Hp, p,vX);
      if (LERS) {
        H0= ZX_init_CRT(H0p, p,vX);
        H1= ZX_init_CRT(H1p, p,vX);
      }
    }
    else
    {
      GEN qp = muliu(q,p);
      stable = ZX_incremental_CRT(&H, Hp, q,qp, p);
      if (LERS) {
        stable &= ZX_incremental_CRT(&H0,H0p, q,qp, p);
        stable &= ZX_incremental_CRT(&H1,H1p, q,qp, p);
      }
      q = qp;
    }
    /* could make it probabilistic for H ? [e.g if stable twice, etc]
     * Probabilistic anyway for H0, H1 */
    if (DEBUGLEVEL>5)
      msgtimer("resultant mod %ld (bound 2^%ld, stable=%ld)", p,expi(q),stable);
    if (stable && (ulong)expi(q) >= bound) break; /* DONE */
    if (low_stack(lim, stack_lim(av,2)))
    {
      GEN *gptr[4]; gptr[0] = &H; gptr[1] = &q; gptr[2] = &H0; gptr[3] = &H1;
      if (DEBUGMEM>1) pari_warn(warnmem,"ZY_ZXY_rnfequation");
      gerepilemany(av2,gptr,LERS? 4: 2); 
    }
  }
END:
  setvarn(H, vX); if (delvar) (void)delete_var();
  if (LERS)
  {
    *LERS = mkvec2(H0,H1);
    gerepileall(av, 2, &H, LERS);
    return H;
  }
  return gerepilecopy(av, H);
}

GEN
ZY_ZXY_rnfequation(GEN A, GEN B, long *lambda)
{
  return ZY_ZXY_resultant_all(A, B, lambda, NULL);
}

/* If lambda = NULL, return caract(Mod(B, A)), A,B in Z[X].
 * Otherwise find a small lambda such that caract (Mod(B + lambda X, A)) is
 * squarefree */
GEN
ZX_caract_sqf(GEN A, GEN B, long *lambda, long v)
{
  pari_sp av = avma;
  GEN B0, R, a;
  long dB;
  int delvar;

  if (v < 0) v = 0;
  switch (typ(B))
  {
    case t_POL: dB = degpol(B); if (dB > 0) break;
      B = dB? gel(B,2): gen_0; /* fall through */
    default:
      if (lambda) { B = scalarpol(B,varn(A)); dB = 0; break;}
      return gerepileupto(av, gpowgs(gsub(pol_x[v], B), degpol(A)));
  }
  delvar = 0;
  if (varn(A) == 0)
  {
    long v0 = fetch_var(); delvar = 1;
    A = shallowcopy(A); setvarn(A,v0);
    B = shallowcopy(B); setvarn(B,v0);
  }
  B0 = cgetg(4, t_POL);
  B0[1] = evalsigne(1);
  gel(B0,2) = gneg_i(B);
  gel(B0,3) = gen_1;
  R = ZY_ZXY_rnfequation(A, B0, lambda);
  if (delvar) (void)delete_var();
  setvarn(R, v); a = leading_term(A);
  if (!gcmp1(a)) R = gdiv(R, powiu(a, dB));
  return gerepileupto(av, R);
}


/* B may be in Q[X], but assume A and result are integral */
GEN
ZX_caract(GEN A, GEN B, long v)
{
  return (degpol(A) < 16) ? caractducos(A,B,v): ZX_caract_sqf(A,B, NULL, v);
}

static GEN
trivial_case(GEN A, GEN B)
{
  long d;
  if (typ(A) == t_INT) return powiu(A, degpol(B));
  d = degpol(A);
  if (d == 0) return trivial_case(gel(A,2),B);
  if (d < 0) return gen_0;
  return NULL;
}

/* Res(A, B/dB), assuming the A,B in Z[X] and result is integer */
GEN
ZX_resultant_all(GEN A, GEN B, GEN dB, ulong bound)
{
  ulong Hp, dp, p;
  pari_sp av = avma, av2, lim;
  long degA;
  int stable;
  GEN q, a, b, H;
  byteptr d;

  if ((H = trivial_case(A,B)) || (H = trivial_case(B,A))) return H;
  q = H = NULL;
  av2 = avma; lim = stack_lim(av,2);
  degA = degpol(A);
  if (!bound)
  {
    bound = ZY_ZXY_ResBound(A, B, dB);
    if (bound > 50000)
    {
      long eA = gexpo(A), eB = gexpo(B), prec = nbits2prec(max(eA,eB));
      for(;; prec = (prec-1)<<1)
      {
        GEN run = real_1(prec);
        GEN R = subres(gmul(A, run), gmul(B, run));
        bound = gexpo(R) + 1;
        if (!gcmp0(R)) break;
      }
      if (dB) bound -= (long)(dbllog2(dB)*degA);
    }
  }
  if (DEBUGLEVEL>4) fprintferr("bound for resultant: 2^%ld\n",bound);
  d = init_modular(&p);
  check_theta(bound);

  dp = 1; /* denominator mod p */
  for(;;)
  {
    NEXT_PRIME_VIADIFF_CHECK(p,d);
    if (dB) { dp = smodis(dB, p); if (!dp) continue; }

    a = ZX_to_Flx(A, p);
    b = ZX_to_Flx(B, p);
    Hp= Flx_resultant(a, b, p);
    if (dp != 1) Hp = Fl_mul(Hp, Fl_pow(Fl_inv(dp,p), degA, p), p);
    if (!H)
    {
      stable = 0; q = utoipos(p);
      H = stoi(Fl_center(Hp, p, p>>1));
    }
    else /* could make it probabilistic ??? [e.g if stable twice, etc] */
    {
      GEN qp = muliu(q,p);
      stable = Z_incremental_CRT(&H, Hp, q,qp, p);
      q = qp;
    }
    if (DEBUGLEVEL>5)
      msgtimer("resultant mod %ld (bound 2^%ld, stable = %d)",p,expi(q),stable);
    if (stable && (ulong)expi(q) >= bound) break; /* DONE */
    if (low_stack(lim, stack_lim(av,2)))
    {
      GEN *gptr[2]; gptr[0] = &H; gptr[1] = &q;
      if (DEBUGMEM>1) pari_warn(warnmem,"ZX_resultant");
      gerepilemany(av2,gptr, 2);
    }
  }
  return gerepileuptoint(av, icopy(H));
}

GEN
ZX_resultant(GEN A, GEN B) { return ZX_resultant_all(A,B,NULL,0); }

GEN
ZX_QX_resultant(GEN A, GEN B)
{
  GEN c, d, n, R;
  pari_sp av = avma;
  B = Q_primitive_part(B, &c);
  if (!c) return ZX_resultant(A,B);
  n = numer(c);
  d = denom(c); if (is_pm1(d)) d = NULL;
  R = ZX_resultant_all(A, B, d, 0);
  if (!is_pm1(n)) R = mulii(R, powiu(n, degpol(A)));
  return gerepileuptoint(av, R);
}

/* assume x has integral coefficients */
GEN
ZX_disc_all(GEN x, ulong bound)
{
  pari_sp av = avma;
  GEN l, d = ZX_resultant_all(x, derivpol(x), NULL, bound);
  l = leading_term(x); if (!gcmp1(l)) d = diviiexact(d,l);
  if (degpol(x) & 2) d = negi(d);
  return gerepileuptoint(av,d);
}
GEN ZX_disc(GEN x) { return ZX_disc_all(x,0); }

int
ZX_is_squarefree(GEN x)
{
  pari_sp av = avma;
  GEN d = modulargcd(x,derivpol(x));
  int r = (lg(d) == 3); avma = av; return r;
}

static GEN
_gcd(GEN a, GEN b)
{
  if (!a) a = gen_1;
  if (!b) b = gen_1;
  return ggcd(a,b);
}

/* ceil( || p ||_oo / lc(p) ) */
static GEN
maxnorm(GEN p)
{
  long i, n = degpol(p), av = avma;
  GEN x, m = gen_0;

  p += 2;
  for (i=0; i<n; i++)
  {
    x = gel(p,i);
    if (absi_cmp(x,m) > 0) m = x;
  }
  m = divii(m, gel(p,n));
  return gerepileuptoint(av, addis(absi(m),1));
}

/* A0 and B0 in Q[X] */
GEN
modulargcd(GEN A0, GEN B0)
{
  GEN a, b, Hp, D, A, B, q, qp, H, g, bound = NULL;
  long m, n;
  ulong p;
  pari_sp av2, av = avma, avlim = stack_lim(av, 1);
  byteptr d;

  if ((typ(A0) | typ(B0)) !=t_POL) pari_err(notpoler,"modulargcd");
  if (!signe(A0)) return gcopy(B0);
  if (!signe(B0)) return gcopy(A0);
  A = primitive_part(A0, &a); check_ZX(A, "modulargcd");
  B = primitive_part(B0, &b); check_ZX(B, "modulargcd");
  D = _gcd(a,b);
  if (varn(A) != varn(B)) pari_err(talker,"different variables in modulargcd");
 
  /* A, B in Z[X] */
  g = gcdii(leading_term(A), leading_term(B)); /* multiple of lead(gcd) */
  if (is_pm1(g)) g = NULL;
  if (degpol(A) < degpol(B)) swap(A, B);
  n = 1 + degpol(B); /* > degree(gcd) */

  av2 = avma; H = NULL;
  d = init_modular(&p);
  for(;;)
  {
    NEXT_PRIME_VIADIFF_CHECK(p,d);
    if (g && !umodiu(g,p)) continue;
    a = ZX_to_Flx(A, p);
    b = ZX_to_Flx(B, p); Hp = Flx_gcd_i(a,b, p);
    m = degpol(Hp);
    if (m == 0) { H = pol_1[varn(A0)]; break; } /* coprime. DONE */
    if (m > n) continue; /* p | Res(A/G, B/G). Discard */

    if (!g) /* make sure lead(H) = g mod p */
      Hp = Flx_normalize(Hp, p);
    else
    {
      ulong t = Fl_mul(umodiu(g, p), Fl_inv(Hp[m+2],p), p);
      Hp = Flx_Fl_mul(Hp, t, p);
    }
    if (m < n)
    { /* First time or degree drop [all previous p were as above; restart]. */
      H = ZX_init_CRT(Hp,p,varn(A0));
      q = utoipos(p); n = m; continue;
    }
    if (DEBUGLEVEL>5)
      msgtimer("gcd mod %lu (bound 2^%ld)", p,expi(q));

    qp = muliu(q,p);
    if (ZX_incremental_CRT(&H, Hp, q, qp, p))
    { /* H stable: check divisibility */
      if (g) {
        if (!bound)
        {
          GEN mA = maxnorm(A), mB = maxnorm(B);
          if (cmpii(mA, mB) > 0) mA = mB;
          bound = gclone( shifti(mulii(mA, g), n+1) );
          if (DEBUGLEVEL>5)
            msgtimer("bound 2^%ld. Goal 2^%ld", expi(q),expi(bound));
        }
        if (cmpii(qp, bound) < 0) goto next;
        H = primpart(H);
        gunclone(bound); break;
      }
      if (gcmp0(RgX_rem(A,H)) && gcmp0(RgX_rem(B,H))) break; /* DONE */
      if (DEBUGLEVEL) fprintferr("modulargcd: trial division failed");
    }
next:
    q = qp;
    if (low_stack(avlim, stack_lim(av,1)))
    {
      if (DEBUGMEM>1) pari_warn(warnmem,"modulargcd");
      gerepileall(av2, 2, &H, &q);
    }
  }
  return gerepileupto(av, gmul(D,H));
}

/* lift(1 / Mod(A,B)). B0 a t_POL, A0 a scalar or a t_POL. Rational coeffs */
GEN
QXQ_inv(GEN A0, GEN B0)
{
  GEN a,b,D,A,B,q,qp,Up,Vp,U,V,res;
  long stable;
  ulong p;
  pari_sp av2, av = avma, avlim = stack_lim(av, 1);
  byteptr d;

  if (typ(B0) != t_POL) pari_err(notpoler,"QXQ_inv");
  if (typ(A0) != t_POL)
  {
    if (is_scalar_t(typ(A0))) return ginv(A0);
    pari_err(notpoler,"QXQ_inv");
  }
  if (degpol(A0) < 15) return ginvmod(A0,B0);
  A = Q_primitive_part(A0, &D);
  B = Q_primpart(B0);
  /* A, B in Z[X] */
  av2 = avma; U = NULL;
  d = init_modular(&p);
  for(;;)
  {
    NEXT_PRIME_VIADIFF_CHECK(p,d);
    a = ZX_to_Flx(A, p);
    b = ZX_to_Flx(B, p);
    /* if p | Res(A/G, B/G), discard */
    if (!Flx_extresultant(b,a,p, &Vp,&Up)) continue;

    if (!U)
    { /* First time */
      U = ZX_init_CRT(Up,p,varn(A0));
      V = ZX_init_CRT(Vp,p,varn(A0));
      q = utoipos(p); continue;
    }
    if (DEBUGLEVEL>5) msgtimer("QXQ_inv: mod %ld (bound 2^%ld)", p,expi(q));
    qp = muliu(q,p);
    stable  = ZX_incremental_CRT(&U, Up, q,qp, p);
    stable &= ZX_incremental_CRT(&V, Vp, q,qp, p);
    if (stable)
    { /* all stable: check divisibility */
      res = gadd(gmul(A,U), gmul(B,V));
      if (degpol(res) == 0) break; /* DONE */
      if (DEBUGLEVEL) fprintferr("QXQ_inv: char 0 check failed");
    }
    q = qp;
    if (low_stack(avlim, stack_lim(av,1)))
    {
      if (DEBUGMEM>1) pari_warn(warnmem,"QXQ_inv");
      gerepileall(av2, 3, &q,&U,&V);
    }
  }
  D = D? gmul(D,res): res;
  return gerepileupto(av, gdiv(U,D));
}

/* irreducible (unitary) polynomial of degree n over Fp */
GEN
ffinit_rand(GEN p,long n)
{
  pari_sp av = avma;
  GEN pol;

  for(;; avma = av)
  {
    pol = gadd(monomial(gen_1, n, 0), FpX_rand(n-1,0, p));
    if (FpX_is_irred(pol, p)) break;
  }
  return pol;
}

GEN
FpX_direct_compositum(GEN A, GEN B, GEN p)
{
  GEN C,a,b,x;
  a = shallowcopy(A); setvarn(a, MAXVARN);
  b = shallowcopy(B); setvarn(b, MAXVARN);
  x = gadd(pol_x[0], pol_x[MAXVARN]);
  C = FpY_FpXY_resultant(a, poleval(b,x),p);
  return C;
}

GEN
FpX_compositum(GEN A, GEN B, GEN p)
{
  GEN C, a,b;
  long k;

  a = shallowcopy(A); setvarn(a, MAXVARN);
  b = shallowcopy(B); setvarn(b, MAXVARN);
  for (k = 1;; k = next_lambda(k))
  {
    GEN x = gadd(pol_x[0], gmulsg(k, pol_x[MAXVARN]));
    C = FpY_FpXY_resultant(a, poleval(b,x),p);
    if (FpX_is_squarefree(C, p)) break;
  }
  return C;
}

/* return an extension of degree 2^l of F_2, assume l > 0 
 * using Adleman-Lenstra Algorithm.
 * Not stack clean. */
static GEN
f2init(long l)
{
  long i;
  GEN q, T, S;

  if (l == 1) return cyclo(3, MAXVARN);

  S = mkpoln(4, gen_1,gen_1,gen_0,gen_0); /* #(#^2 + #) */
  setvarn(S, MAXVARN);
  q = mkpoln(3, gen_1,gen_1, S); /* X^2 + X + #(#^2+#) */

  /* x^4+x+1, irred over F_2, minimal polynomial of a root of q */
  T = mkpoln(5, gen_1,gen_0,gen_0,gen_1,gen_1);
  for (i=2; i<l; i++)
  { /* q = X^2 + X + a(#) irred. over K = F2[#] / (T(#))
     * ==> X^2 + X + a(#) b irred. over K for any root b of q
     * ==> X^2 + X + (b^2+b)b */
    setvarn(T, MAXVARN);
    T = FpY_FpXY_resultant(T, q, gen_2);
    /* T = minimal polynomial of b over F2 */
  }
  return T;
}

/* return an extension of degree p^l of F_p, assume l > 0 
 * using Adleman-Lenstra Algorithm, see below.
 * Not stack clean. */
GEN
ffinit_Artin_Shreier(GEN ip, long l)
{
  long i;
  long p=itos(ip);
  GEN xp,yp,y2pm1;
  GEN P, Q;
  xp=monomial(gen_1,p,0);
  P = ZX_sub(xp, deg1pol_i(gen_1,gen_1,0));
  if (l == 1) return P;
  
  yp=monomial(gen_1,p,MAXVARN);
  y2pm1=monomial(gen_1,2*p-1,MAXVARN);
  Q = gsub(ZX_sub(xp, pol_x[0]), ZX_sub(y2pm1, yp));
  for (i = 2; i <= l; ++i)
  {
    setvarn(P,MAXVARN);
    P = FpY_FpXY_resultant(P, Q, ip);
  }
  return P;
}


/*Check if subcyclo(n,l,0) is irreducible modulo p*/
static long
fpinit_check(GEN p, long n, long l)
{
  pari_sp ltop=avma;
  long q,o;
  if (!uisprime(n)) {avma=ltop; return 0;}
  q = smodis(p,n);
  if (!q) {avma=ltop; return 0;}
  o = itos(order(mkintmodu(q,n)));
  avma = ltop;
  return ( cgcd((n-1)/o,l) == 1 );
}

/* let k=2 if p%4==1, and k=4 else and assume k*p does not divide l.
 * Return an irreducible polynomial of degree l over F_p.
 * This a variant of an algorithm of Adleman and Lenstra
 * "Finding irreducible polynomials over finite fields",
 * ACM, 1986 (5)  350--355
 * Not stack clean.
 */
static GEN
fpinit(GEN p, long l)
{
  ulong n = 1+l, k = 1;
  while (!fpinit_check(p,n,l)) { n += l; k++; }
  if (DEBUGLEVEL>=4)
    fprintferr("FFInit: using subcyclo(%ld, %ld)\n",n,l);
  return FpX_red(subcyclo(n,l,0),p);
}

static GEN
ffinit_fact(GEN p, long n)
{
  GEN F = (GEN)factoru_pow(n)[3];
  GEN P; /* pol */
  long i;
  /* If n is even, then F[1] is 2^bfffo(n)*/
  if (!odd(n) && equaliu(p, 2))
    P = f2init(vals(n));
  else
    P = fpinit(p, F[1]);
  for (i = 2; i < lg(F); ++i)
    P = FpX_direct_compositum(fpinit(p, F[i]), P, p);
  return P;
}

static GEN
ffinit_nofact(GEN p, long n)
{
  GEN P, Q = NULL;
  if (lgefint(p)==3)
  {
    ulong lp = p[2], q;
    long v = u_lvalrem(n,lp,&q);
    if (v>0)
    {
      if (lp==2) Q = f2init(v);
      else       Q = fpinit(p,n/q);
      n = q;
    }
  }
  if (n==1) P = Q;
  else
  {
    P = fpinit(p, n);
    if (Q) P = FpX_direct_compositum(P, Q, p);
  }
  return P;
}

static GEN
init_Fq_i(GEN p, long n, long v)
{
  GEN P;
  if (n <= 0) pari_err(talker,"non positive degree in ffinit");
  if (typ(p) != t_INT) pari_err(typeer, "ffinit");
  if (v < 0) v = 0;
  if (n == 1) return pol_x[v];
  /*If easy case, use cyclo*/
  if (fpinit_check(p, n + 1, n)) return cyclo(n + 1, v);
  if (lgefint(p)-2 < BITS_IN_LONG-(long)bfffo(n))
    P = ffinit_fact(p,n);
  else
    P = ffinit_nofact(p,n);
  setvarn(P, v); return P;
}
GEN
init_Fq(GEN p, long n, long v)
{
  pari_sp av = avma;
  return gerepileupto(av, init_Fq_i(p, n, v));
}
GEN
ffinit(GEN p, long n, long v)
{
  pari_sp av = avma;
  return gerepileupto(av, FpX_to_mod(init_Fq_i(p, n, v), p));
}
