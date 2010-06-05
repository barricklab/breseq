/* $Id: gen3.c 11395 2008-12-08 10:46:23Z kb $

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
/**                         (third part)                           **/
/**                                                                **/
/********************************************************************/
#include "pari.h"
#include "paripriv.h"

/********************************************************************/
/**                                                                **/
/**                 PRINCIPAL VARIABLE NUMBER                      **/
/**                                                                **/
/********************************************************************/
long
gvar(GEN x)
{
  long i, v, w;
  switch(typ(x))
  {
    case t_POL: case t_SER: return varn(x);
    case t_POLMOD: return varn(gel(x,1));
    case t_RFRAC:  return varn(gel(x,2));
    case t_VEC: case t_COL: case t_MAT:
      v = BIGINT;
      for (i=1; i < lg(x); i++) { w=gvar(gel(x,i)); if (w<v) v=w; }
      return v;
    case t_VECSMALL:
    case t_STR: 
    case t_LIST: 
      pari_err(typeer, "gvar");
  }
  return BIGINT;
}
/* T main polynomial in R[X], A auxiliary in R[X] (possibly degree 0).
 * Guess and return the main variable of R */
static long
var2_aux(GEN T, GEN A)
{
  long a = gvar2(T);
  long b = (typ(A) == t_POL && varn(A) == varn(T))? gvar2(A): gvar(A);
  if (a < b) a = b;
  return a;
}
static long
var2_rfrac(GEN x)  { return var2_aux(gel(x,2), gel(x,1)); }
static long
var2_polmod(GEN x) { return var2_aux(gel(x,1), gel(x,2)); }

/* main variable of x, with the convention that the "natural" main
 * variable of a POLMOD is mute, so we want the next one. */
static long
gvar9(GEN x)
{
  return (typ(x) == t_POLMOD)? var2_polmod(x): gvar(x);
}

/* main variable of the ring over wich x is defined */
long
gvar2(GEN x)
{
  long i, v, w;
  switch(typ(x))
  {
    case t_POLMOD:
      return var2_polmod(x);
    case t_POL: case t_SER:
      v = BIGINT;
      for (i=2; i < lg(x); i++) { w = gvar9(gel(x,i)); if (w<v) v=w; }
      return v;
    case t_RFRAC:
      return var2_rfrac(x);
    case t_VEC: case t_COL: case t_MAT:
      v = BIGINT;
      for (i=1; i < lg(x); i++) { w = gvar2(gel(x,i)); if (w<v) v=w; }
      return v;
  }
  return BIGINT;
}

GEN
gpolvar(GEN x)
{
  long v;
  if (typ(x)==t_PADIC) return gcopy( gel(x,2) );
  v = gvar(x);
  if (v==BIGINT) pari_err(typeer,"gpolvar");
  return pol_x[v];
}

/*******************************************************************/
/*                                                                 */
/*                    PRECISION OF SCALAR OBJECTS                  */
/*                                                                 */
/*******************************************************************/
static long
prec0(long e) { return (e < 0)? 2 - (e >> TWOPOTBITS_IN_LONG): 2; }
static long
precREAL(GEN x) { return signe(x) ? lg(x): prec0(expo(x)); }
/* t t_REAL, s an exact non-complex type. Return precision(|t| + |s|) */
static long
precrealexact(GEN t, GEN s) {
  long l, e = gexpo(s);
  if (e == -(long)HIGHEXPOBIT) return precREAL(t);
  if (e < 0) e = 0;
  e -= expo(t);
  if (!signe(t)) return prec0(-e);
  l = lg(t);
  return (e > 0)? l + (e >> TWOPOTBITS_IN_LONG): l;
}
long
precision(GEN z)
{
  long tx = typ(z), e, ex, ey, lz, lx, ly;

  if (tx == t_REAL) return precREAL(z);
  if (tx == t_COMPLEX)
  { /* ~ precision(|x| + |y|) */
    GEN x = gel(z,1), y = gel(z,2);
    if (typ(x) != t_REAL) {
      if (typ(y) != t_REAL) return 0;
      return precrealexact(y, x);
    }
    if (typ(y) != t_REAL) return precrealexact(x, y);
    /* x, y are t_REALs, cf addrr_sign */
    ex = expo(x);
    ey = expo(y);
    e = ey - ex;
    if (!signe(x)) {
      if (!signe(y)) return prec0( min(ex,ey) );
      if (e < 0) return prec0(ex);
      lz = 3 + (e >> TWOPOTBITS_IN_LONG);
      ly = lg(y); if (lz > ly) lz = ly;
      return lz;
    }
    if (!signe(y)) {
      if (e > 0) return prec0(ey);
      lz = 3 + ((-e)>>TWOPOTBITS_IN_LONG);
      lx = lg(x); if (lz > lx) lz = lx;
      return lz;
    }
    if (e < 0) { swap(x, y); e = -e; }
    lx = lg(x);
    ly = lg(y);
    if (e) {
      long d = e >> TWOPOTBITS_IN_LONG, l = ly-d;
      return (l > lx)? lx + d: ly;
    }
    return min(lx, ly);
  }
  return 0;
}

long
gprecision(GEN x)
{
  long tx=typ(x),lx=lg(x),i,k,l;

  if (is_scalar_t(tx)) return precision(x);
  switch(tx)
  {
    case t_POL: case t_VEC: case t_COL: case t_MAT:
      k=VERYBIGINT;
      for (i=lontyp[tx]; i<lx; i++)
      {
        l = gprecision(gel(x,i));
	if (l && l<k) k = l;
      }
      return (k==VERYBIGINT)? 0: k;

    case t_RFRAC:
    {
      k=gprecision(gel(x,1));
      l=gprecision(gel(x,2)); if (l && (!k || l<k)) k=l;
      return k;
    }
    case t_QFR:
      return gprecision(gel(x,4));
  }
  return 0;
}

GEN
ggprecision(GEN x)
{
  long a = gprecision(x);
  return utoipos(a ? prec2ndec(a): VERYBIGINT);
}

GEN
precision0(GEN x, long n) { return n? gprec(x,n): ggprecision(x); }

/* attention: precision p-adique absolue */
long
padicprec(GEN x, GEN p)
{
  long i,s,t,lx=lg(x),tx=typ(x);

  switch(tx)
  {
    case t_INT: case t_FRAC:
      return VERYBIGINT;

    case t_INTMOD:
      return Z_pval(gel(x,1),p);

    case t_PADIC:
      if (!gequal(gel(x,2),p))
	pari_err(talker,"not the same prime in padicprec");
      return precp(x)+valp(x);

    case t_POL:
    case t_COMPLEX: case t_QUAD: case t_POLMOD: case t_SER: case t_RFRAC:
    case t_VEC: case t_COL: case t_MAT:
      for (s=VERYBIGINT, i=lontyp[tx]; i<lx; i++)
      {
        t = padicprec(gel(x,i),p); if (t<s) s = t;
      }
      return s;
  }
  pari_err(typeer,"padicprec");
  return 0; /* not reached */
}

#define DEGREE0 -VERYBIGINT
/* Degree of x (scalar, t_POL, t_RFRAC) wrt variable v if v >= 0,
 * wrt to main variable if v < 0.
 */
long
poldegree(GEN x, long v)
{
  long tx = typ(x), lx,w,i,d;

  if (is_scalar_t(tx)) return gcmp0(x)? DEGREE0: 0;
  switch(tx)
  {
    case t_POL:
      if (!signe(x)) return DEGREE0;
      w = varn(x);
      if (v < 0 || v == w) return degpol(x);
      if (v < w) return 0;
      lx = lg(x); d = DEGREE0;
      for (i=2; i<lx; i++)
      {
        long e = poldegree(gel(x,i), v);
        if (e > d) d = e;
      }
      return d;

    case t_RFRAC:
      if (gcmp0(gel(x,1))) return DEGREE0;
      return poldegree(gel(x,1),v) - poldegree(gel(x,2),v);
  }
  pari_err(typeer,"degree");
  return 0; /* not reached  */
}

long
degree(GEN x)
{
  return poldegree(x,-1);
}

/* If v<0, leading coeff with respect to the main variable, otherwise wrt v. */
GEN
pollead(GEN x, long v)
{
  long l, tx = typ(x), w;
  pari_sp av;
  GEN xinit;

  if (is_scalar_t(tx)) return gcopy(x);
  w = varn(x);
  switch(tx)
  {
    case t_POL:
      if (v < 0 || v == w)
      {
	l=lg(x);
	return (l==2)? gen_0: gcopy(gel(x,l-1));
      }
      break;

    case t_SER:
      if (v < 0 || v == w) return signe(x)? gcopy(gel(x,2)): gen_0;
      break;

    default:
      pari_err(typeer,"pollead");
      return NULL; /* not reached */
  }
  if (v < w) return gcopy(x);
  av = avma; xinit = x;
  x = gsubst(gsubst(x,w,pol_x[MAXVARN]),v,pol_x[0]);
  if (gvar(x)) { avma = av; return gcopy(xinit);}
  tx = typ(x);
  if (tx == t_POL) {
    l = lg(x); if (l == 2) { avma = av; return gen_0; }
    x = gel(x,l-1);
  }
  else if (tx == t_SER) {
    if (!signe(x)) { avma = av; return gen_0;}
    x = gel(x,2);
  } else pari_err(typeer,"pollead");
  return gerepileupto(av, gsubst(x,MAXVARN,pol_x[w]));
}

/* returns 1 if there's a real component in the structure, 0 otherwise */
int
isinexactreal(GEN x)
{
  long tx=typ(x),lx,i;

  if (is_const_t(tx))
  {
    switch(tx)
    {
      case t_REAL:
        return 1;

      case t_COMPLEX: case t_QUAD:
        return (typ(x[1])==t_REAL || typ(x[2])==t_REAL);
    }
    return 0;
  }
  switch(tx)
  {
    case t_QFR: case t_QFI:
      return 0;

    case t_RFRAC: case t_POLMOD:
      return isinexactreal(gel(x,1)) || isinexactreal(gel(x,2));
  }
  if (is_noncalc_t(tx)) return 0;
  lx = lg(x);
  for (i=lontyp[tx]; i<lx; i++)
    if (isinexactreal(gel(x,i))) return 1;
  return 0;
}

/* returns 1 if there's an inexact component in the structure, and
 * 0 otherwise. */
int
isinexact(GEN x)
{
  long tx = typ(x), lx, i;

  switch(tx)
  {
    case t_REAL: case t_PADIC: case t_SER:
      return 1;
    case t_INT: case t_INTMOD: case t_FRAC: case t_QFR: case t_QFI:
      return 0;
    case t_COMPLEX: case t_QUAD: case t_RFRAC: case t_POLMOD:
      return isinexact(gel(x,1)) || isinexact(gel(x,2));
    case t_POL: case t_VEC: case t_COL: case t_MAT:
      lx = lg(x);
      for (i=lontyp[tx]; i<lx; i++)
        if (isinexact(gel(x,i))) return 1;
      return 0;
    case t_LIST:
      lx = lgeflist(x);
      for (i=lontyp[tx]; i<lx; i++)
        if (isinexact(gel(x,i))) return 1;
      return 0;
  }
  return 0;
}

int
isexactzeroscalar(GEN g)
{
  switch (typ(g))
  {
    case t_INT:
      return !signe(g);
    case t_INTMOD: case t_POLMOD:
      return isexactzeroscalar(gel(g,2));
    case t_FRAC: case t_RFRAC:
      return isexactzeroscalar(gel(g,1));
    case t_COMPLEX:
      return isexactzeroscalar(gel(g,1)) && isexactzeroscalar(gel(g,2));
    case t_QUAD:
      return isexactzeroscalar(gel(g,2)) && isexactzeroscalar(gel(g,3));
    case t_POL: return lg(g) == 2;
  }
  return 0;
}

int
isexactzero(GEN g)
{
  long i;
  switch (typ(g))
  {
    case t_INT:
      return !signe(g);
    case t_INTMOD: case t_POLMOD:
      return isexactzero(gel(g,2));
    case t_COMPLEX:
      return isexactzero(gel(g,1)) && isexactzero(gel(g,2));
    case t_QUAD:
      return isexactzero(gel(g,2)) && isexactzero(gel(g,3));
    case t_POL: return lg(g) == 2;
    case t_VEC: case t_COL: case t_MAT:
      for (i=lg(g)-1; i; i--)
	if (!isexactzero(gel(g,i))) return 0;
      return 1;
  }
  return 0;
}

int
iscomplex(GEN x)
{
  switch(typ(x))
  {
    case t_INT: case t_REAL: case t_FRAC:
      return 0;
    case t_COMPLEX:
      return !gcmp0(gel(x,2));
    case t_QUAD:
      return signe(gmael(x,1,2)) > 0;
  }
  pari_err(typeer,"iscomplex");
  return 0; /* not reached */
}

int
ismonome(GEN x)
{
  long i;
  if (typ(x)!=t_POL || !signe(x)) return 0;
  for (i=lg(x)-2; i>1; i--)
    if (!isexactzero(gel(x,i))) return 0;
  return 1;
}

/*******************************************************************/
/*                                                                 */
/*                    GENERIC REMAINDER                            */
/*                                                                 */
/*******************************************************************/
/* euclidean quotient for scalars of admissible types */
static GEN
_quot(GEN x, GEN y)
{
  GEN q = gdiv(x,y), f = gfloor(q);
  if (gsigne(y) < 0 && !gequal(f,q)) f = gadd(f,gen_1);
  return f;
}
static GEN
_quotgs(GEN x, long y)
{
  GEN q = gdivgs(x,y), f = gfloor(q);
  if (y < 0 && !gequal(f,q)) f = gadd(f,gen_1);
  return f;
}
static GEN
quot(GEN x, GEN y)
{
  pari_sp av = avma;
  return gerepileupto(av, _quot(x, y));
}
static GEN
quotgs(GEN x, long y)
{
  pari_sp av = avma;
  return gerepileupto(av, _quotgs(x, y));
}

GEN
gmod(GEN x, GEN y)
{
  pari_sp av, tetpil;
  long i,lx,ty, tx = typ(x);
  GEN z,p1;

  if (is_matvec_t(tx))
  {
    lx=lg(x); z=cgetg(lx,tx);
    for (i=1; i<lx; i++) gel(z,i) = gmod(gel(x,i),y);
    return z;
  }
  ty = typ(y);
  switch(ty)
  {
    case t_INT:
      switch(tx)
      {
	case t_INT:
	  return modii(x,y);

	case t_INTMOD: z=cgetg(3,tx);
          gel(z,1) = gcdii(gel(x,1),y);
	  gel(z,2) = modii(gel(x,2),gel(z,1)); return z;

	case t_FRAC:
	  av=avma;
	  p1=mulii(gel(x,1),Fp_inv(gel(x,2),y));
	  tetpil=avma; return gerepile(av,tetpil,modii(p1,y));

	case t_QUAD: z=cgetg(4,tx);
          gel(z,1) = gcopy(gel(x,1));
	  gel(z,2) = gmod(gel(x,2),y);
          gel(z,3) = gmod(gel(x,3),y); return z;

	case t_PADIC: return padic_to_Fp(x, y);
	case t_POLMOD: case t_POL:
	  return gen_0;
	case t_REAL: /* NB: conflicting semantic with lift(x * Mod(1,y)). */
	  av = avma;
          return gerepileupto(av, gadd(x, gneg(gmul(_quot(x,y),y))));

	default: pari_err(operf,"%",x,y);
      }

    case t_REAL: case t_FRAC:
      switch(tx)
      {
	case t_INT: case t_REAL: case t_FRAC:
	  av = avma;
          return gerepileupto(av, gsub(x, gmul(_quot(x,y),y)));

	case t_POLMOD: case t_POL:
	  return gen_0;

	default: pari_err(operf,"%",x,y);
      }

    case t_POL:
      if (is_scalar_t(tx))
      {
        if (tx!=t_POLMOD || varncmp(varn(x[1]), varn(y)) > 0)
          return degpol(y)? gcopy(x): gen_0;
	if (varn(x[1])!=varn(y)) return gen_0;
        z=cgetg(3,t_POLMOD);
        gel(z,1) = ggcd(gel(x,1),y);
        gel(z,2) = grem(gel(x,2),gel(z,1)); return z;
      }
      switch(tx)
      {
	case t_POL:
	  return grem(x,y);

	case t_RFRAC:
	  av=avma;
	  p1=gmul(gel(x,1),ginvmod(gel(x,2),y)); tetpil=avma;
          return gerepile(av,tetpil,grem(p1,y));

        case t_SER:
          if (ismonome(y) && varn(x) == varn(y))
          {
            long d = degpol(y);
            if (lg(x)-2 + valp(x) < d) pari_err(operi,"%",x,y);
            av = avma; 
            return gerepileupto(av, gmod(ser2rfrac_i(x), y));
          }
	default: pari_err(operf,"%",x,y);
      }
  }
  pari_err(operf,"%",x,y);
  return NULL; /* not reached */
}

GEN
gmodgs(GEN x, long y)
{
  ulong u;
  long i, lx, tx = typ(x);
  GEN z;
  if (is_matvec_t(tx))
  {
    lx=lg(x); z=cgetg(lx,tx);
    for (i=1; i<lx; i++) gel(z,i) = gmodgs(gel(x,i),y);
    return z;
  }
  switch(tx)
  {
    case t_INT:
      return modis(x,y);

    case t_INTMOD: z=cgetg(3,tx);
      i = cgcd(smodis(gel(x,1), y), y);
      gel(z,1) = utoi(i);
      gel(z,2) = modis(gel(x,2), i); return z;

    case t_FRAC:
      u = (ulong)labs(y);
      return utoi( Fl_div(umodiu(gel(x,1), u),
                          umodiu(gel(x,2), u), u) );

    case t_QUAD: z=cgetg(4,tx);
      gel(z,1) = gcopy(gel(x,1));
      gel(z,2) = gmodgs(gel(x,2),y);
      gel(z,3) = gmodgs(gel(x,3),y); return z;

    case t_PADIC: return padic_to_Fp(x, stoi(y));
    case t_POLMOD: case t_POL:
      return gen_0;
  }
  pari_err(operf,"%",x,stoi(y));
  return NULL; /* not reached */
}

GEN
gmodulsg(long x, GEN y)
{
  GEN z;

  switch(typ(y))
  {
    case t_INT: z = cgetg(3,t_INTMOD);
      gel(z,1) = absi(y); 
      gel(z,2) = modsi(x,y); return z;

    case t_POL: z = cgetg(3,t_POLMOD);
      gel(z,1) = gcopy(y);
      gel(z,2) = stoi(x); return z;
  }
  pari_err(operf,"%",stoi(x),y); return NULL; /* not reached */
}

GEN
gmodulss(long x, long y)
{
  GEN z = cgetg(3,t_INTMOD);
  y = labs(y);
  gel(z,1) = stoi(y);
  gel(z,2) = modss(x, y); return z;
}

static GEN 
specialmod(GEN x, GEN y)
{
  GEN z = gmod(x,y);
  if (varncmp(gvar(z), varn(y)) < 0) z = gen_0;
  return z;
}

GEN
gmodulo(GEN x,GEN y)
{
  long tx=typ(x),l,i;
  GEN z;

  if (is_matvec_t(tx))
  {
    l=lg(x); z=cgetg(l,tx);
    for (i=1; i<l; i++) gel(z,i) = gmodulo(gel(x,i),y);
    return z;
  }
  switch(typ(y))
  {
    case t_INT: z = cgetg(3,t_INTMOD);
      gel(z,1) = absi(y);
      gel(z,2) = Rg_to_Fp(x,y); return z;

    case t_POL: z = cgetg(3,t_POLMOD);
      gel(z,1) = gcopy(y);
      if (is_scalar_t(tx))
      {
        gel(z,2) = (lg(y) > 3)? gcopy(x): gmod(x,y);
        return z;
      }
      if (tx!=t_POL && tx != t_RFRAC && tx!=t_SER) break;
      gel(z,2) = specialmod(x,y); return z;
  }
  pari_err(operf,"%",x,y); return NULL; /* not reached */
}

GEN
Mod0(GEN x,GEN y,long flag)
{
  switch(flag)
  {
    case 0:
    case 1: return gmodulo(x,y);
    default: pari_err(flagerr,"Mod");
  }
  return NULL; /* not reached */
}

/*******************************************************************/
/*                                                                 */
/*                 GENERIC EUCLIDEAN DIVISION                      */
/*                                                                 */
/*******************************************************************/

GEN
gdivent(GEN x, GEN y)
{
  long tx = typ(x);

  if (is_matvec_t(tx))
  {
    long i, lx = lg(x);
    GEN z = cgetg(lx,tx);
    for (i=1; i<lx; i++) gel(z,i) = gdivent(gel(x,i),y);
    return z;
  }
  switch(typ(y))
  {
    case t_INT:
      switch(tx)
      { /* equal to, but more efficient than, quot(x,y) */
        case t_INT: return truedivii(x,y);
        case t_REAL: case t_FRAC: return quot(x,y);
        case t_POL: return gdiv(x,y);
      }
      break;
    case t_REAL: case t_FRAC: return quot(x,y);
    case t_POL:
      if (is_scalar_t(tx))
      {
        if (tx == t_POLMOD) break;
        return degpol(y)? gen_0: gdiv(x,y);
      }
      if (tx == t_POL) return gdeuc(x,y);
  }
  pari_err(operf,"\\",x,y);
  return NULL; /* not reached */
}

GEN
gdiventgs(GEN x, long y)
{
  long tx = typ(x);

  if (is_matvec_t(tx))
  {
    long i, lx = lg(x);
    GEN z = cgetg(lx,tx);
    for (i=1; i<lx; i++) gel(z,i) = gdiventgs(gel(x,i),y);
    return z;
  }
  switch(tx)
  { /* equal to, but more efficient than, quotgs(x,y) */
    case t_INT: return truedivis(x,y);
    case t_REAL: case t_FRAC: return quotgs(x,y);
    case t_POL: return gdivgs(x,y);
  }
  pari_err(operf,"\\",x,stoi(y));
  return NULL; /* not reached */
}

/* with remainder */
static GEN
quotrem(GEN x, GEN y, GEN *r)
{
  pari_sp av;
  GEN q = quot(x,y);
  av = avma;
  *r = gerepileupto(av, gsub(x, gmul(q,y)));
  return q;
}

GEN
gdiventres(GEN x, GEN y)
{
  long tx = typ(x);
  GEN z,q,r;

  if (is_matvec_t(tx))
  {
    long i, lx = lg(x);
    z = cgetg(lx,tx);
    for (i=1; i<lx; i++) gel(z,i) = gdiventres(gel(x,i),y);
    return z;
  }
  z = cgetg(3,t_COL);
  switch(typ(y))
  {
    case t_INT:
      switch(tx)
      { /* equal to, but more efficient than next case */
        case t_INT:
          gel(z,1) = truedvmdii(x,y,(GEN*)(z+2));
          return z;
        case t_REAL: case t_FRAC:
          q = quotrem(x,y,&r);
          gel(z,1) = q;
          gel(z,2) = r; return z;
        case t_POL:
          gel(z,1) = gdiv(x,y);
          gel(z,2) = gen_0; return z;
      }
      break;
    case t_REAL: case t_FRAC:
          q = quotrem(x,y,&r);
          gel(z,1) = q;
          gel(z,2) = r; return z;
    case t_POL:
      if (is_scalar_t(tx))
      {
        if (tx == t_POLMOD) break;
        if (degpol(y))
        {
          q = gen_0;
          r = gcopy(x);
        }
        else
        {
          q = gdiv(x,y);
          r = gen_0;
        }
        gel(z,1) = q;
        gel(z,2) = r; return z;
      }
      if (tx == t_POL)
      {
        gel(z,1) = poldivrem(x,y,(GEN *)(z+2));
        return z;
      }
  }
  pari_err(operf,"\\",x,y);
  return NULL; /* not reached */
}

GEN
divrem(GEN x, GEN y, long v)
{
  pari_sp av = avma;
  long vx, vy;
  GEN q, r;
  if (v < 0 || typ(y) != t_POL || typ(x) != t_POL) return gdiventres(x,y);
  vx = varn(x); if (vx != v) x = swap_vars(x,v);
  vy = varn(y); if (vy != v) y = swap_vars(y,v);
  q = poldivrem(x,y, &r);
  if (v && (vx != v || vy != v))
  {
    q = gsubst(q, v, pol_x[v]); /* poleval broken for t_RFRAC, subst is safe */
    r = gsubst(r, v, pol_x[v]);
  }
  return gerepilecopy(av, mkcol2(q, r));
}

static int
is_scal(long t) { return t==t_INT || t==t_FRAC; }

GEN
diviiround(GEN x, GEN y)
{
  pari_sp av1, av = avma;
  GEN q,r;
  int fl;

  q = dvmdii(x,y,&r); /* q = x/y rounded towards 0, sgn(r)=sgn(x) */
  if (r==gen_0) return q;
  av1 = avma;
  fl = absi_cmp(shifti(r,1),y);
  avma = av1; cgiv(r);
  if (fl >= 0) /* If 2*|r| >= |y| */
  {
    long sz = signe(x)*signe(y);
    if (fl || sz > 0) q = gerepileuptoint(av, addis(q,sz));
  }
  return q;
}

/* If x and y are not both scalars, same as gdivent.
 * Otherwise, compute the quotient x/y, rounded to the nearest integer
 * (towards +oo in case of tie). */
GEN
gdivround(GEN x, GEN y)
{
  pari_sp av1, av;
  long tx=typ(x),ty=typ(y);
  GEN q,r;
  int fl;

  if (tx==t_INT && ty==t_INT) return diviiround(x,y);
  av = avma;
  if (is_scal(tx) && is_scal(ty))
  { /* same as diviiround but less efficient */
    q = quotrem(x,y,&r);
    av1 = avma;
    fl = gcmp(gmul2n(gabs(r,0),1), gabs(y,0));
    avma = av1; cgiv(r);
    if (fl >= 0) /* If 2*|r| >= |y| */
    {
      long sz = gsigne(y);
      if (fl || sz > 0) q = gerepileupto(av, gaddgs(q, sz));
    }
    return q;
  }
  if (is_matvec_t(tx))
  {
    long i,lx = lg(x);
    GEN z = cgetg(lx,tx);
    for (i=1; i<lx; i++) gel(z,i) = gdivround(gel(x,i),y);
    return z;
  }
  return gdivent(x,y);
}

GEN
gdivmod(GEN x, GEN y, GEN *pr)
{
  long ty,tx=typ(x);

  if (tx==t_INT)
  {
    ty=typ(y);
    if (ty==t_INT) return dvmdii(x,y,pr);
    if (ty==t_POL) { *pr=gcopy(x); return gen_0; }
    pari_err(typeer,"gdivmod");
  }
  if (tx!=t_POL) pari_err(typeer,"gdivmod");
  return poldivrem(x,y,pr);
}

/*******************************************************************/
/*                                                                 */
/*                               SHIFT                             */
/*                                                                 */
/*******************************************************************/

/* Shift tronque si n<0 (multiplication tronquee par 2^n)  */

GEN
gshift(GEN x, long n)
{
  long i, lx, tx = typ(x);
  GEN y;

  switch(tx)
  {
    case t_INT:
      return shifti(x,n);
    case t_REAL:
      return shiftr(x,n);

    case t_VEC: case t_COL: case t_MAT:
      lx = lg(x); y = cgetg(lx,tx);
      for (i=1; i<lx; i++) gel(y,i) = gshift(gel(x,i),n);
      return y;
  }
  return gmul2n(x,n);
}

/*******************************************************************/
/*                                                                 */
/*                              INVERSE                            */
/*                                                                 */
/*******************************************************************/
GEN
mpinv(GEN b)
{
  long i, l1, l = lg(b), e = expo(b), s = signe(b);
  GEN x = cgetr(l), a = mpcopy(b);
  double t;

  a[1] = evalexpo(0) | evalsigne(1);
  for (i = 3; i < l; i++) x[i] = 0;
  t = (((double)HIGHBIT) * HIGHBIT) / (double)(ulong)a[2];
  if (((ulong)t) & HIGHBIT)
    x[1] = evalexpo(0) | evalsigne(1);
  else {
    t *= 2;
    x[1] = evalexpo(-1) | evalsigne(1);
  }
  x[2] = (ulong)t;
  l1 = 1; l -= 2;
  while (l1 < l)
  {
    l1 <<= 1; if (l1 > l) l1 = l;
    setlg(a, l1 + 2);
    setlg(x, l1 + 2);
    /* TODO: mulrr(a,x) should be a half product (the higher half is known).
     * mulrr(x, ) already is */
    affrr(addrr(x, mulrr(x, subsr(1, mulrr(a,x)))), x);
    avma = (pari_sp)a;
  }
  x[1] = evalexpo(expo(x)-e) | evalsigne(s);
  avma = (pari_sp)x; return x;
}

GEN
inv_ser(GEN b)
{
  pari_sp av = avma, av2, lim;
  long i, j, le, l = lg(b), e = valp(b), v = varn(b);
  GEN E, y, x = cgetg(l, t_SER), a = shallowcopy(b);

  if (!signe(b)) pari_err(gdiver);

  for (j = 3; j < l; j++) gel(x,j) = gen_0;
  gel(x,2) = ginv(gel(b,2));
  a[1] = x[1] = evalvalp(0) | evalvarn(v) | evalsigne(1);
  E = Newton_exponents(l - 2);
  av2 = avma; lim = stack_lim(av2, 2);
  le = lg(E)-1;
  for (i = le; i > 1; i--) {
    long l1 = E[i-1], l2 = E[i];
    setlg(a, l1 + 2);
    setlg(x, l1 + 2);
    /* TODO: gmul(a,x) should be a half product (the higher half is known) */
    y = gmul(x, gsubsg(1, gmul(a,x))) - l2;
    for (j = l2+2; j < l1+2; j++) x[j] = y[j];
    if (low_stack(lim, stack_lim(av2,2)))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"inv_ser");
      y = gerepilecopy(av2, x);
      for (j = 2; j < l1+2; j++) x[j] = y[j];
    }
  }
  x[1] = evalvalp(valp(x)-e) | evalvarn(v) | evalsigne(1);
  return gerepilecopy(av, x);
}

GEN
ginv(GEN x)
{
  long s;
  pari_sp av, tetpil;
  GEN X, z, y, p1, p2;

  switch(typ(x))
  {
    case t_INT:
      if (is_pm1(x)) return icopy(x);
      s = signe(x); if (!s) pari_err(gdiver);
      z = cgetg(3,t_FRAC);
      gel(z,1) = s<0? gen_m1: gen_1;
      gel(z,2) = absi(x); return z;

    case t_REAL:
      return divsr(1,x);

    case t_INTMOD: z=cgetg(3,t_INTMOD);
      gel(z,1) = icopy(gel(x,1));
      gel(z,2) = Fp_inv(gel(x,2),gel(x,1)); return z;

    case t_FRAC:
      s = signe(x[1]); if (!s) pari_err(gdiver);
      if (is_pm1(x[1]))
        return s>0? icopy(gel(x,2)): negi(gel(x,2));
      z = cgetg(3,t_FRAC);
      gel(z,1) = icopy(gel(x,2));
      gel(z,2) = icopy(gel(x,1));
      if (s < 0)
      {
	setsigne(z[1],-signe(z[1]));
	setsigne(z[2],1);
      }
      return z;

    case t_COMPLEX: case t_QUAD:
      av=avma; p1=gnorm(x); p2=gconj(x); tetpil=avma;
      return gerepile(av,tetpil,gdiv(p2,p1));

    case t_PADIC: z = cgetg(5,t_PADIC);
      if (!signe(x[4])) pari_err(gdiver);
      z[1] = evalprecp(precp(x)) | evalvalp(-valp(x));
      gel(z,2) = icopy(gel(x,2));
      gel(z,3) = icopy(gel(x,3));
      gel(z,4) = Fp_inv(gel(x,4),gel(z,3)); return z;

    case t_POLMOD: z = cgetg(3,t_POLMOD);
      X = gel(x,1); gel(z,1) = gcopy(X);
      if (degpol(X) == 2) { /* optimized for speed */
        av = avma;
        gel(z,2) = gerepileupto(av, gdiv(quad_polmod_conj(gel(x,2), X),
                                  quad_polmod_norm(gel(x,2), X)) );
      }
      else gel(z,2) = ginvmod(gel(x,2), X);
      return z;

    case t_POL: return gred_rfrac_simple(gen_1,x);
    case t_SER: return gdiv(gen_1,x);

    case t_RFRAC:
    {
      GEN n = gel(x,1), d = gel(x,2);
      pari_sp av = avma, ltop;
      if (gcmp0(n)) pari_err(gdiver);

      n = simplify_i(n);
      if (typ(n) != t_POL || varn(n) != varn(d))
      {
        if (gcmp1(n)) { avma = av; return gcopy(d); }
        ltop = avma;
        z = RgX_Rg_div(d,n);
      } else {
        ltop = avma;
        z = cgetg(3,t_RFRAC);
        gel(z,1) = gcopy(d);
        gel(z,2) = gcopy(n);
      }
      stackdummy(av, ltop);
      return z;
    }

    case t_QFR:
      av = avma; z = cgetg(5, t_QFR);
      gel(z,1) = gel(x,1);
      gel(z,2) = negi( gel(x,2) );
      gel(z,3) = gel(x,3);
      gel(z,4) = negr( gel(x,4) );
      return gerepileupto(av, redreal(z));

    case t_QFI:
      y = gcopy(x);
      if (!equalii(gel(x,1),gel(x,2)) && !equalii(gel(x,1),gel(x,3)))
	setsigne(y[2],-signe(y[2]));
      return y;
    case t_MAT:
      return (lg(x)==1)? cgetg(1,t_MAT): invmat(x);
    case t_VECSMALL:
    {
      long i,lx = lg(x);
      y = cgetg(lx,t_VECSMALL);
      for (i=1; i<lx; i++)
      {
        long xi=x[i];
	  if (xi<1 || xi>=lx) pari_err(talker,"incorrect permtuation to inverse");
        y[xi] = i;
      }
      return y;
    }
  }
  pari_err(typeer,"inverse");
  return NULL; /* not reached */
}

/*******************************************************************/
/*                                                                 */
/*           SUBSTITUTION DANS UN POLYNOME OU UNE SERIE            */
/*                                                                 */
/*******************************************************************/

/* Convert t_SER --> t_POL, ignoring valp. INTERNAL ! */
GEN
ser2pol_i(GEN x, long lx)
{
  long i = lx-1;
  GEN y;
  while (i > 1 && isexactzero(gel(x,i))) i--;
  y = cgetg(i+1, t_POL); y[1] = x[1] & ~VALPBITS;
  for ( ; i > 1; i--) y[i] = x[i];
  return y;
}

/*
   subst_poly(pol, from, to) =
   { local(t='subst_poly_t, M);

     \\ if fraction
     M = numerator(from) - t * denominator(from);
     \\ else
     M = from - t;
     subst(pol % M, t, to)
   }
 */
GEN
gsubst_expr(GEN pol, GEN from, GEN to)
{
  pari_sp av = avma;
  long v = fetch_var();		/* XXX Need fetch_var_low_priority() */
  GEN tmp;

  from = simplify_i(from);
  switch (typ(from)) {
  case t_RFRAC: /* M= numerator(from) - t * denominator(from) */
    tmp = gsub(gel(from,1), gmul(pol_x[v], gel(from,2)));
    break;
  default:
    tmp = gsub(from, pol_x[v]);	/* M = from - t */
  }

  if (v <= gvar(from)) pari_err(talker, "subst: unexpected variable precedence");
  tmp = gmul(pol, mkpolmod(gen_1, tmp));
  if (typ(tmp) == t_POLMOD)
    tmp = gel(tmp,2);			/* optimize lift */
  else					/* Vector? */
    tmp = lift0(tmp, gvar(from));
  tmp = gsubst(tmp, v, to);
  (void)delete_var();
  return gerepilecopy(av, tmp);
}

GEN
gsubstpol(GEN x, GEN T, GEN y)
{
  pari_sp av;
  long d, v;
  GEN deflated;

  if (typ(T) != t_POL || !ismonome(T) || !gcmp1(leading_term(T)))
    return gsubst_expr(x,T,y);
  d = degpol(T); v = varn(T);
  if (d == 1) return gsubst(x, v, y);
  av = avma;
  CATCH(cant_deflate) {
    avma = av;
    return gsubst_expr(x,T,y);      
  } TRY {
    deflated = gdeflate(x, v, d);
  } ENDCATCH
  return gerepilecopy(av, gsubst(deflated, v, y));
}

GEN
gsubst(GEN x, long v, GEN y)
{
  long tx = typ(x), ty = typ(y), lx = lg(x), ly = lg(y);
  long l, vx, vy, e, ex, ey, i, j, k, jb;
  pari_sp av, lim;
  GEN t,p1,p2,z;

  if (ty==t_MAT)
  {
    if (ly==1) return cgetg(1,t_MAT);
    if (ly != lg(y[1]))
      pari_err(talker,"forbidden substitution by a non square matrix");
  } else if (is_graphicvec_t(ty))
    pari_err(talker,"forbidden substitution by a vector");

  if (is_scalar_t(tx))
  {
    if (tx!=t_POLMOD || v <= varn(x[1]))
    {
      if (ty==t_MAT) return gscalmat(x,ly-1);
      return gcopy(x);
    }
    av=avma;
    p1=gsubst(gel(x,1),v,y); vx=varn(p1);
    p2=gsubst(gel(x,2),v,y); vy=gvar(p2);
    if (typ(p1)!=t_POL)
      pari_err(talker,"forbidden substitution in a scalar type");
    if (varncmp(vy, vx) >= 0) return gerepileupto(av, gmodulo(p2,p1));
    lx = lg(p2);
    z = cgetg(lx,t_POL); z[1] = p2[1];
    for (i=2; i<lx; i++) gel(z,i) = gmodulo(gel(p2,i),p1);
    return gerepileupto(av, normalizepol_i(z,lx));
  }

  switch(tx)
  {
    case t_POL:
      if (lx==2)
        return (ty==t_MAT)? gscalmat(gen_0,ly-1): gen_0;

      vx = varn(x);
      if (varncmp(vx, v) < 0)
      {
        if (varncmp(gvar(y), vx) > 0)
        { /* easy special case */
          z = cgetg(lx, t_POL); z[1] = x[1];
          for (i=2; i<lx; i++) gel(z,i) = gsubst(gel(x,i),v,y);
          return normalizepol_i(z,lx);
        }
        /* general case */
	av=avma; p1=pol_x[vx]; z= gsubst(gel(x,lx-1),v,y);
	for (i=lx-1; i>2; i--) z=gadd(gmul(z,p1),gsubst(gel(x,i-1),v,y));
	return gerepileupto(av,z);
      }
      /* v <= vx */
      if (ty!=t_MAT)
        return varncmp(vx,v) > 0 ? gcopy(x): poleval(x,y);

      if (varncmp(vx, v) > 0) return gscalmat(x,ly-1);
      if (lx==3) return gscalmat(gel(x,2),ly-1);
      av=avma; z=gel(x,lx-1);
      for (i=lx-1; i>2; i--) z=gaddmat(gel(x,i-1),gmul(z,y));
      return gerepileupto(av,z);

    case t_SER:
      vx = varn(x);
      if (varncmp(vx, v) > 0) return (ty==t_MAT)? gscalmat(x,ly-1): gcopy(x);
      ex = valp(x);
      if (varncmp(vx, v) < 0)
      { /* FIXME: improve this */
        av = avma; p1 = ser2pol_i(x, lx);
        z = tayl(gsubst(p1,v,y), vx, lx-2);
        if (ex) z = gmul(z, monomial(gen_1,ex,vx));
        return gerepileupto(av, z);
      }
      switch(ty) /* here vx == v */
      {
        case t_SER:
	  ey = valp(y);
          vy = varn(y);
	  if (ey < 1) return zeroser(vy, ey*(ex+lx-2));
          if (lg(x) == 2) return zeroser(vy, ey*ex);
	  if (vy != vx)
	  {
	    av = avma; z = zeroser(vy,0);
	    for (i=lx-1; i>=2; i--) z = gadd(gel(x,i), gmul(y,z));
	    if (ex) z = gmul(z, gpowgs(y,ex));
	    return gerepileupto(av,z);
	  }
	  l = (lx-2)*ey+2;
	  if (ex) { if (l>ly) l = ly; }
	  else if (lx != 3)
          {
            long l2;
            for (i = 3; i < lx; i++)
              if (!isexactzero(gel(x,i))) break;
            l2 = (i-2)*ey + (gcmp0(y)? 2 : ly);
            if (l > l2) l = l2;
          }
          p2 = ex? gpowgs(y, ex): NULL;

	  av = avma; lim=stack_lim(av,1);
          t = shallowcopy(y);
          if (l < ly) setlg(t, l);
          z = scalarser(gel(x,2),varn(y),l-2);
	  for (i=3,jb=ey; jb<=l-2; i++,jb+=ey)
	  {
            if (i < lx) {
              for (j=jb+2; j<min(l, jb+ly); j++)
                gel(z,j) = gadd(gel(z,j), gmul(gel(x,i),gel(t,j-jb)));
            }
	    for (j=l-1-jb-ey; j>1; j--)
	    {
	      p1 = gen_0;
	      for (k=2; k<j; k++)
		p1 = gadd(p1, gmul(gel(t,j-k+2),gel(y,k)));
	      gel(t,j) = gadd(p1, gmul(gel(t,2),gel(y,j)));
	    }
            if (low_stack(lim, stack_lim(av,1)))
	    {
	      if(DEBUGMEM>1) pari_warn(warnmem,"gsubst");
	      gerepileall(av,2, &z,&t);
	    }
	  }
	  if (!p2) return gerepilecopy(av,z);
          return gerepileupto(av, gmul(z,p2));

        case t_POL: case t_RFRAC:
          if (isexactzero(y)) return scalarser(gel(x,2),v,lx-2);
          vy = gvar(y); e = gval(y,vy);
          if (e <= 0)
            pari_err(talker,"non positive valuation in a series substitution");
	  av = avma; p1 = gsubst(ser2pol_i(x, lg(x)), v, y);
          z = gmul(gpowgs(y, ex), tayl(p1, vy, e*(lx-2)));
	  return gerepileupto(av, z);

        default:
          pari_err(talker,"non polynomial or series type substituted in a series");
      }
      break;

    case t_RFRAC: av=avma;
      p1=gsubst(gel(x,1),v,y);
      p2=gsubst(gel(x,2),v,y); return gerepileupto(av, gdiv(p1,p2));

    case t_VEC: case t_COL: case t_MAT: z=cgetg(lx,tx);
      for (i=1; i<lx; i++) gel(z,i) = gsubst(gel(x,i),v,y);
      return z;
  }
  return gcopy(x);
}

GEN
gsubstvec(GEN e, GEN v, GEN r)
{
  pari_sp ltop=avma;
  long i,l=lg(v);
  GEN w,z;
  if ( !is_vec_t(typ(v)) || !is_vec_t(typ(r)) )
    pari_err(typeer,"substvec");
  if (lg(r)!=l) 
    pari_err(talker,"different number of variables and values in substvec");
  w=cgetg(l,t_VECSMALL);
  z=cgetg(l,t_VECSMALL);
  for(i=1;i<l;i++)
  { 
    GEN T=(GEN)v[i];
    if (typ(T) != t_POL || !ismonome(T) || !gcmp1(leading_term(T)))
      pari_err(talker,"not a variable in substvec");
    w[i]=varn(T);
    z[i]=fetch_var();
  }
  for(i=1;i<l;i++) e = gsubst(e,w[i],pol_x[z[i]]);
  for(i=1;i<l;i++) e = gsubst(e,z[i],(GEN)r[i]);
  for(i=1;i<l;i++) (void)delete_var();
  return gerepileupto(ltop,e);
}

/*******************************************************************/
/*                                                                 */
/*                SERIE RECIPROQUE D'UNE SERIE                     */
/*                                                                 */
/*******************************************************************/

GEN
recip(GEN x)
{
  long v=varn(x), lx = lg(x);
  pari_sp tetpil, av=avma;
  GEN p1,p2,a,y,u;

  if (typ(x)!=t_SER) pari_err(talker,"not a series in serreverse");
  if (valp(x)!=1 || lx < 3)
    pari_err(talker,"valuation not equal to 1 in serreverse");

  a=gel(x,2);
  if (gcmp1(a))
  {
    long i, j, k, mi;
    pari_sp lim=stack_lim(av, 2);

    mi = lx-1; while (mi>=3 && gcmp0(gel(x,mi))) mi--;
    u = cgetg(lx,t_SER);
    y = cgetg(lx,t_SER);
    u[1] = y[1] = evalsigne(1) | evalvalp(1) | evalvarn(v);
    gel(u,2) = gel(y,2) = gen_1;
    if (lx > 3)
    {
      gel(u,3) = gmulsg(-2,gel(x,3));
      gel(y,3) = gneg(gel(x,3));
    }
    for (i=3; i<lx-1; )
    {
      pari_sp av2;
      for (j=3; j<i+1; j++)
      {
        av2 = avma; p1 = gel(x,j);
        for (k = max(3,j+2-mi); k < j; k++)
          p1 = gadd(p1, gmul(gel(u,k),gel(x,j-k+2)));
        p1 = gneg(p1);
        gel(u,j) = gerepileupto(av2, gadd(gel(u,j), p1));
      }
      av2 = avma;
      p1 = gmulsg(i,gel(x,i+1));
      for (k = 2; k < min(i,mi); k++)
      {
        p2 = gmul(gel(x,k+1),gel(u,i-k+2));
        p1 = gadd(p1, gmulsg(k,p2));
      }
      i++;
      gel(u,i) = gerepileupto(av2, gneg(p1));
      gel(y,i) = gdivgs(gel(u,i), i-1);
      if (low_stack(lim, stack_lim(av,2)))
      {
	if(DEBUGMEM>1) pari_warn(warnmem,"recip");
	for(k=i+1; k<lx; k++) gel(u,k) = gel(y,k) = gen_0; /* dummy */
	gerepileall(av,2, &u,&y);
      }
    }
    return gerepilecopy(av,y);
  }
  y = gdiv(x,a); gel(y,2) = gen_1; y = recip(y);
  a = gdiv(pol_x[v],a); tetpil = avma;
  return gerepile(av,tetpil, gsubst(y,v,a));
}

/*******************************************************************/
/*                                                                 */
/*                    DERIVATION ET INTEGRATION                    */
/*                                                                 */
/*******************************************************************/
GEN
derivpol(GEN x)
{
  long i,lx = lg(x)-1;
  GEN y;

  if (lx<3) return zeropol(varn(x));
  y = cgetg(lx,t_POL);
  for (i=2; i<lx ; i++) gel(y,i) = gmulsg(i-1,gel(x,i+1));
  y[1] = x[1]; return normalizepol_i(y,i);
}

GEN
derivser(GEN x)
{
  long i, vx = varn(x), e = valp(x), lx = lg(x);
  GEN y;
  if (lx == 2) return zeroser(vx,e? e-1: 0);
  if (e)
  {
    y = cgetg(lx,t_SER); y[1] = evalvalp(e-1) | evalvarn(vx);
    for (i=2; i<lx; i++) gel(y,i) = gmulsg(i+e-2,gel(x,i));
  } else {
    if (lx == 3) return zeroser(vx, 0);
    lx--;
    y = cgetg(lx,t_SER); y[1] = evalvalp(0) | evalvarn(vx);
    for (i=2; i<lx; i++) gel(y,i) = gmulsg(i-1,gel(x,i+1));
  }
  return normalize(y);
}

GEN
deriv(GEN x, long v)
{
  long lx, vx, tx, i, j;
  pari_sp av;
  GEN y;

  tx = typ(x); if (is_const_t(tx)) return gen_0;
  if (v < 0) v = gvar9(x);
  switch(tx)
  {
    case t_POLMOD:
      if (v<=varn(x[1])) return gen_0;
      y = cgetg(3,t_POLMOD);
      gel(y,1) = gcopy(gel(x,1));
      gel(y,2) = deriv(gel(x,2),v); return y;

    case t_POL:
      vx = varn(x);
      if (varncmp(vx, v) > 0) return gen_0;
      if (varncmp(vx, v) == 0) return derivpol(x);
      lx = lg(x); y = cgetg(lx,t_POL);
      y[1] = x[1];
      for (i=2; i<lx; i++) gel(y,i) = deriv(gel(x,i),v);
      return normalizepol_i(y,i);

    case t_SER:
      vx = varn(x);
      if (varncmp(vx, v) > 0) return gen_0;
      if (varncmp(vx, v) == 0) return derivser(x);
      lx = lg(x); y = cgetg(lx, t_SER);
      y[1] = x[1];
      for (j=2; j<lx; j++) gel(y,j) = deriv(gel(x,j),v);
      return normalize(y);

    case t_RFRAC: {
      GEN a = gel(x,1), b = gel(x,2), bp, b0, d, t;
      y = cgetg(3,t_RFRAC); av = avma;

      bp = deriv(b, v);
      d = ggcd(bp, b);
      if (gcmp1(d)) {
        d = gadd(gmul(b, deriv(a,v)), gmul(gneg_i(a), bp));
        if (isexactzero(d)) return gerepileupto((pari_sp)(y+3), d);
        gel(y,1) = gerepileupto(av, d);
        gel(y,2) = gsqr(b); return y;
      }
      b0 = gdivexact(b, d);
      bp = gdivexact(bp,d);
      a = gadd(gmul(b0, deriv(a,v)), gmul(gneg_i(a), bp));
      if (isexactzero(a)) return gerepileupto((pari_sp)(y+3), a);
      t = ggcd(a, d);
      if (!gcmp1(t)) { a = gdivexact(a, t); d = gdivexact(d, t); }
      gel(y,1) = a;
      gel(y,2) = gmul(d, gsqr(b0));
      return gerepilecopy((pari_sp)(y+3), y);
    }

    case t_VEC: case t_COL: case t_MAT: lx=lg(x); y=cgetg(lx,tx);
      for (i=1; i<lx; i++) gel(y,i) = deriv(gel(x,i),v);
      return y;
  }
  pari_err(typeer,"deriv");
  return NULL; /* not reached */
}

/********************************************************************/
/**                                                                **/
/**                         TAYLOR SERIES                          **/
/**                                                                **/
/********************************************************************/
static GEN
tayl_vec(long v, long vx) {
  GEN y = cgetg(v+2,t_VEC);
  long i;
  for (i=0; i<v; i++) gel(y,i+1) = pol_x[i];
  gel(y,vx+1) = pol_x[v];
  gel(y,v+1)  = pol_x[vx]; return y;
}

GEN
tayl(GEN x, long v, long precS)
{
  long vx = gvar9(x);
  pari_sp av = avma;
  GEN y, t;

  if (v <= vx) return gadd(zeroser(v,precS),x);
  y = tayl_vec(v, vx);
  t = tayl(changevar(x,y), vx,precS);
  return gerepileupto(av, changevar(t,y));
}

GEN
ggrando(GEN x, long n)
{
  long m, v;

  switch(typ(x))
  {
  case t_INT:/* bug 3 + O(1). We suppose x is a truc() */
    if (signe(x) <= 0) pari_err(talker,"non-positive argument in O()");
    if (!is_pm1(x)) return zeropadic(x,n);
    /* +/-1 = x^0 */
    v = m = 0; break;
  case t_POL:
    if (!signe(x)) pari_err(talker,"zero argument in O()");
    v = varn(x); if ((ulong)v > MAXVARN) pari_err(talker,"incorrect object in O()");
    m = n * polvaluation(x, NULL); break;
  case t_RFRAC:
    if (!gcmp0((GEN)x[1])) pari_err(talker,"zero argument in O()");
    v = gvar(x); if ((ulong)v > MAXVARN) pari_err(talker,"incorrect object in O()");
    m = n * gval(x,v); break;
    default: pari_err(talker,"incorrect argument in O()");
      v = m = 0; /* not reached */
  }
  return zeroser(v,m);
}

/*******************************************************************/
/*                                                                 */
/*                    FORMAL INTEGRATION                           */
/*                                                                 */
/*******************************************************************/

static GEN
triv_integ(GEN x, long v, long tx, long lx)
{
  GEN y = cgetg(lx,tx);
  long i;

  y[1] = x[1];
  for (i=2; i<lx; i++) gel(y,i) = integ(gel(x,i),v);
  return y;
}

GEN
integ(GEN x, long v)
{
  long lx, tx, e, i, vx, n;
  pari_sp av = avma;
  GEN y,p1;

  tx = typ(x);
  if (v < 0) v = gvar(x);
  if (is_scalar_t(tx))
  {
    if (tx == t_POLMOD && v>varn(x[1]))
    {
      y=cgetg(3,t_POLMOD);
      gel(y,1) = gcopy(gel(x,1));
      gel(y,2) = integ(gel(x,2),v); return y;
    }
    if (gcmp0(x)) return gen_0;

    y = cgetg(4,t_POL);
    y[1] = evalsigne(1) | evalvarn(v);
    gel(y,2) = gen_0;
    gel(y,3) = gcopy(x); return y;
  }

  switch(tx)
  {
    case t_POL:
      vx = varn(x); lx = lg(x);
      if (lx == 2) {
        if (varncmp(vx, v) < 0) v = vx;
        return zeropol(v);
      }
      if (varncmp(vx, v) > 0)
      {
        y = cgetg(4,t_POL);
	y[1] = evalsigne(1) | evalvarn(v);
        gel(y,2) = gen_0;
        gel(y,3) = gcopy(x); return y;
      }
      if (varncmp(vx, v) < 0) return triv_integ(x,v,tx,lx);
      y = cgetg(lx+1,tx); y[1] = x[1]; gel(y,2) = gen_0;
      for (i=3; i<=lx; i++) gel(y,i) = gdivgs(gel(x,i-1),i-2);
      return y;

    case t_SER:
      lx = lg(x); vx = varn(x); e = valp(x);
      if (lx == 2)
      {
        if (vx == v) e++; else if (varncmp(vx, v) < 0) v = vx;
        return zeroser(v, e);
      }
      if (varncmp(vx, v) > 0)
      {
        y = cgetg(4,t_POL);
        y[1] = evalvarn(v) | evalsigne(1);
        gel(y,2) = gen_0;
        gel(y,3) = gcopy(x); return y;
      }
      if (varncmp(vx, v) < 0) return triv_integ(x,v,tx,lx);
      y = cgetg(lx,tx);
      for (i=2; i<lx; i++)
      {
	long j = i+e-1;
        if (!j)
	{ /* should be isexactzero, but try to avoid error */
	  if (gcmp0(gel(x,i))) { gel(y,i) = gen_0; continue; }
          pari_err(talker, "a log appears in intformal");
	}
	else gel(y,i) = gdivgs(gel(x,i),j);
      }
      y[1] = evalsigne(1) | evalvarn(vx) | evalvalp(e+1); return y;

    case t_RFRAC:
      vx = gvar(x);
      if (varncmp(vx, v) > 0)
      {
	y=cgetg(4,t_POL);
	y[1] = signe(x[1])? evalvarn(v) | evalsigne(1)
	                  : evalvarn(v);
        gel(y,2) = gen_0;
        gel(y,3) = gcopy(x); return y;
      }
      if (varncmp(vx, v) < 0)
      {
        p1 = tayl_vec(v, vx);
	y = integ(changevar(x,p1),vx);
	return gerepileupto(av, changevar(y,p1));
      }

      tx = typ(x[1]); i = is_scalar_t(tx)? 0: degpol(x[1]);
      tx = typ(x[2]); n = is_scalar_t(tx)? 0: degpol(x[2]);
      n = i+n + 2;
      y = gdiv(gtrunc(gmul(gel(x,2), integ(tayl(x,v,n),v))), gel(x,2));
      if (!gequal(deriv(y,v),x)) pari_err(talker,"a log/atan appears in intformal");
      if (typ(y)==t_RFRAC && lg(y[1]) == lg(y[2]))
      {
        GEN p2;
	tx=typ(y[1]); p1=is_scalar_t(tx)? gel(y,1): leading_term(gel(y,1));
	tx=typ(y[2]); p2=is_scalar_t(tx)? gel(y,2): leading_term(gel(y,2));
	y=gsub(y, gdiv(p1,p2));
      }
      return gerepileupto(av,y);

    case t_VEC: case t_COL: case t_MAT:
      lx=lg(x); y=cgetg(lx,tx);
      for (i=1; i<lg(x); i++) gel(y,i) = integ(gel(x,i),v);
      return y;
  }
  pari_err(typeer,"integ");
  return NULL; /* not reached */
}

/*******************************************************************/
/*                                                                 */
/*                    PARTIES ENTIERES                             */
/*                                                                 */
/*******************************************************************/

GEN
gfloor(GEN x)
{
  GEN y;
  long i,lx, tx = typ(x);

  switch(tx)
  {
    case t_INT:
    case t_POL: return gcopy(x);
    case t_REAL: return floorr(x);
    case t_FRAC: return truedivii(gel(x,1),gel(x,2));
    case t_RFRAC: return gdeuc(gel(x,1),gel(x,2));
    case t_VEC: case t_COL: case t_MAT:
      lx = lg(x); y = cgetg(lx,tx);
      for (i=1; i<lx; i++) gel(y,i) = gfloor(gel(x,i));
      return y;
  }
  pari_err(typeer,"gfloor");
  return NULL; /* not reached */
}

GEN
gfrac(GEN x)
{
  pari_sp av = avma, tetpil;
  GEN p1 = gneg_i(gfloor(x));
  tetpil = avma; return gerepile(av,tetpil,gadd(x,p1));
}

/* assume x t_REAL */
GEN
ceilr(GEN x) {
  pari_sp av = avma;
  GEN y = floorr(x);
  if (cmpri(x, y)) return gerepileuptoint(av, addsi(1,y));
  return y;
}

GEN
gceil(GEN x)
{
  GEN y, p1;
  long i, lx, tx = typ(x);
  pari_sp av;

  switch(tx)
  {
    case t_INT: case t_POL: return gcopy(x);
    case t_REAL: return ceilr(x);
    case t_FRAC:
      av = avma; y = dvmdii(gel(x,1),gel(x,2),&p1);
      if (p1 != gen_0 && gsigne(x) > 0)
      {
        cgiv(p1);
        return gerepileuptoint(av, addsi(1,y));
      }
      return y;

    case t_RFRAC:
      return gdeuc(gel(x,1),gel(x,2));

    case t_VEC: case t_COL: case t_MAT:
      lx = lg(x); y = cgetg(lx,tx);
      for (i=1; i<lx; i++) gel(y,i) = gceil(gel(x,i));
      return y;
  }
  pari_err(typeer,"gceil");
  return NULL; /* not reached */
}

GEN
round0(GEN x, GEN *pte)
{
  if (pte) { long e; x = grndtoi(x,&e); *pte = stoi(e); }
  return ground(x);
}

/* assume x a t_REAL */
GEN
roundr(GEN x)
{
  long ex, s = signe(x);
  pari_sp av;
  GEN t;
  if (!s || (ex=expo(x)) < -1) return gen_0;
  if (ex == -1) return s>0? gen_1:
                            absrnz_egal2n(x)? gen_0: gen_m1;
  av = avma;
  t = real2n(-1, nbits2prec(ex+1)); /* = 0.5 */
  return gerepileuptoint(av, floorr( addrr(x,t) ));
}

GEN
ground(GEN x)
{
  GEN y;
  long i, lx, tx=typ(x);
  pari_sp av;

  switch(tx)
  {
    case t_INT: case t_INTMOD: case t_QUAD: return gcopy(x);
    case t_REAL: return roundr(x);
    case t_FRAC: return diviiround(gel(x,1), gel(x,2));
    case t_POLMOD: y=cgetg(3,t_POLMOD);
      gel(y,1) = gcopy(gel(x,1));
      gel(y,2) = ground(gel(x,2)); return y;

    case t_COMPLEX:
      av = avma; y = cgetg(3, t_COMPLEX);
      gel(y,2) = ground(gel(x,2));
      if (!signe(y[2])) { avma = av; return ground(gel(x,1)); }
      gel(y,1) = ground(gel(x,1)); return y;

    case t_POL:
      y = init_gen_op(x, tx, &lx, &i);
      for (; i<lx; i++) gel(y,i) = ground(gel(x,i));
      return normalizepol_i(y, lx);
    case t_SER:
      y = init_gen_op(x, tx, &lx, &i);
      for (; i<lx; i++) gel(y,i) = ground(gel(x,i));
      return normalize(y);
    case t_RFRAC: case t_VEC: case t_COL: case t_MAT:
      y = init_gen_op(x, tx, &lx, &i);
      for (; i<lx; i++) gel(y,i) = ground(gel(x,i));
      return y;
  }
  pari_err(typeer,"ground");
  return NULL; /* not reached */
}

/* e = number of error bits on integral part */
GEN
grndtoi(GEN x, long *e)
{
  GEN y, p1;
  long i, tx=typ(x), lx, ex, e1;
  pari_sp av;

  *e = -(long)HIGHEXPOBIT;
  switch(tx)
  {
    case t_INT: case t_INTMOD: case t_QUAD: return gcopy(x);
    case t_FRAC: return diviiround(gel(x,1), gel(x,2));
    case t_REAL:
      ex = expo(x);
      if (!signe(x) || ex < -1) { *e = ex; return gen_0; }
      av = avma;
      p1 = addrr(real2n(-1,nbits2prec(ex+2)), x); e1 = expo(p1);
      if (e1 < 0)
      {
	if (signe(p1) >= 0) { *e = ex; avma = av; return gen_0; }
        *e = expo(addsr(1,x)); avma = av; return gen_m1;
      }
      lx = lg(x); 
      e1 = e1 - bit_accuracy(lx) + 1;
      y = ishiftr_lg(p1, lx, e1);
      if (signe(x) < 0) y = addsi(-1,y);
      y = gerepileuptoint(av,y);

      if (e1 <= 0) { av = avma; e1 = expo(subri(x,y)); avma = av; }
      *e = e1; return y;

    case t_POLMOD: y = cgetg(3,t_POLMOD);
      gel(y,1) = gcopy(gel(x,1));
      gel(y,2) = grndtoi(gel(x,2), e); return y;

    case t_COMPLEX:
      av = avma; y = cgetg(3, t_COMPLEX);
      gel(y,2) = grndtoi(gel(x,2), e);
      if (!signe(y[2])) {
        avma = av;
        y = grndtoi(gel(x,1), &e1);
      }
      else
        gel(y,1) = grndtoi(gel(x,1), &e1);
      if (e1 > *e) *e = e1;
      return y;

    case t_POL:
      y = init_gen_op(x, tx, &lx, &i);
      for (; i<lx; i++)
      {
        gel(y,i) = grndtoi(gel(x,i),&e1);
        if (e1 > *e) *e = e1;
      }
      return normalizepol_i(y, lx);
    case t_SER:
      y = init_gen_op(x, tx, &lx, &i);
      for (; i<lx; i++)
      {
        gel(y,i) = grndtoi(gel(x,i),&e1);
        if (e1 > *e) *e = e1;
      }
      return normalize(y);
    case t_RFRAC: case t_VEC: case t_COL: case t_MAT:
      y = init_gen_op(x, tx, &lx, &i);
      for (; i<lx; i++)
      {
        gel(y,i) = grndtoi(gel(x,i),&e1);
        if (e1 > *e) *e = e1;
      }
      return y;
  }
  pari_err(typeer,"grndtoi");
  return NULL; /* not reached */
}

/* floor(x * 2^s) */
GEN
gfloor2n(GEN x, long s)
{
  GEN z;
  switch(typ(x))
  {
    case t_INT:
      return shifti(x, s);
    case t_REAL:
      return ishiftr(x, s);
    case t_COMPLEX:
      z = cgetg(3, t_COMPLEX);
      gel(z,1) = gfloor2n(gel(x,1), s);
      gel(z,2) = gfloor2n(gel(x,2), s);
      return z;
    default: pari_err(typeer,"gfloor2n");
      return NULL; /* not reached */
  }
}

/* e = number of error bits on integral part */
GEN
gcvtoi(GEN x, long *e)
{
  long tx = typ(x), lx, i, ex, e1;
  GEN y;

  if (tx == t_REAL)
  {
    ex = expo(x); if (ex < 0) { *e = ex; return gen_0; }
    lx = lg(x); e1 = ex - bit_accuracy(lx) + 1;
    y = ishiftr_lg(x, lx, e1);
    if (e1 <= 0) { pari_sp av = avma; e1 = expo(subri(x,y)); avma = av; }
    *e = e1; return y;
  }
  *e = -(long)HIGHEXPOBIT;
  if (is_matvec_t(tx))
  {
    lx = lg(x); y = cgetg(lx,tx);
    for (i=1; i<lx; i++)
    {
      gel(y,i) = gcvtoi(gel(x,i),&e1);
      if (e1 > *e) *e = e1;
    }
    return y;
  }
  return gtrunc(x);
}

int
isint(GEN n, GEN *ptk)
{
  switch(typ(n))
  {
    case t_INT: *ptk = n; return 1;
    case t_REAL: {
      pari_sp av0 = avma;
      GEN z = floorr(n);
      pari_sp av = avma;
      long s = signe(subri(n, z));
      if (s) { avma = av0; return 0; }
      *ptk = z; avma = av; return 1;
    }
    case t_FRAC:    return 0;
    case t_COMPLEX: return gcmp0(gel(n,2)) && isint(gel(n,1),ptk);
    case t_QUAD:    return gcmp0(gel(n,3)) && isint(gel(n,2),ptk);
    default: pari_err(typeer,"isint"); return 0; /* not reached */
  }
}

int
issmall(GEN n, long *ptk)
{
  pari_sp av = avma;
  GEN z;
  long k;
  if (!isint(n, &z)) return 0;
  k = itos_or_0(z); avma = av;
  if (k || lgefint(z) == 2) { *ptk = k; return k; }
  return 0;
}

/* smallest integer greater than any incarnations of the real x
 * [avoid mpfloor() and "precision loss in truncation"] */
GEN
ceil_safe(GEN x)
{
  pari_sp av = avma;
  long e;
  GEN y = gcvtoi(x,&e);
  if (e < 0) e = 0;
  y = addii(y, int2n(e));
  return gerepileuptoint(av, y);
}

GEN
ser2rfrac_i(GEN x)
{
  long e = valp(x);
  GEN a = ser2pol_i(x, lg(x));
  if (e) {
    if (e > 0) a = RgX_shift_shallow(a, e);
    else a = gred_rfrac_simple(a, monomial(gen_1, -e, varn(a))); 
  }
  return a;
}

static GEN
ser2rfrac(GEN x)
{
  pari_sp av = avma;
  return gerepilecopy(av, ser2rfrac_i(x));
}

GEN
gtrunc(GEN x)
{
  long tx=typ(x), i, v;
  pari_sp av;
  GEN y;

  switch(tx)
  {
    case t_INT: case t_POL:
      return gcopy(x);

    case t_REAL:
      return mptrunc(x);

    case t_FRAC:
      return divii(gel(x,1),gel(x,2));

    case t_PADIC:
      if (!signe(x[4])) return gen_0;
      v = valp(x);
      if (!v) return gcopy(gel(x,4));
      if (v>0)
      { /* here p^v is an integer */
        av = avma; y = powiu(gel(x,2),v);
        return gerepileuptoint(av, mulii(y,gel(x,4)));
      }
      y=cgetg(3,t_FRAC);
      gel(y,1) = icopy(gel(x,4));
      gel(y,2) = gpowgs(gel(x,2),-v);
      return y;

    case t_RFRAC:
      return gdeuc(gel(x,1),gel(x,2));

    case t_SER:
      return ser2rfrac(x);

    case t_VEC: case t_COL: case t_MAT:
    {
      long lx = lg(x); y = cgetg(lx,tx);
      for (i=1; i<lx; i++) gel(y,i) = gtrunc(gel(x,i));
      return y;
    }
  }
  pari_err(typeer,"gtrunc");
  return NULL; /* not reached */
}

GEN
trunc0(GEN x, GEN *pte)
{
  if (pte) { long e; x = gcvtoi(x,&e); *pte = stoi(e); }
  return gtrunc(x);
}
/*******************************************************************/
/*                                                                 */
/*                  CONVERSIONS -->  INT, POL & SER                */
/*                                                                 */
/*******************************************************************/

/* return a_(n-1) B^(n-1) + ... + a_0, where B = 2^32. 
 * The a_i are 32bits integers */
GEN
mkintn(long n, ...)
{
  va_list ap;
  GEN x, y;
  long i;
#ifdef LONG_IS_64BIT
  long e = (n&1);
  n = (n+1) >> 1;
#endif
  va_start(ap,n);
  x = cgeti(n+2); 
  x[1] = evallgefint(n+2) | evalsigne(1);
  y = int_MSW(x);
  for (i=0; i <n; i++)
  {
#ifdef LONG_IS_64BIT
    ulong a = (e && !i)? 0: va_arg(ap, long);
    ulong b = va_arg(ap, long);
    *y = (a << 32) | b;
#else
    *y = va_arg(ap, long);
#endif
    y = int_precW(y);
  }
  va_end(ap);
  return int_normalize(x, 0);
}

/* 2^32 a + b */
GEN
u2toi(ulong a, ulong b)
{
  GEN x;
  if (!a && !b) return gen_0;
#ifdef LONG_IS_64BIT
  x = cgeti(3);
  x[1] = evallgefint(3)|evalsigne(1);
  x[2] = ((a << 32) | b);
#else
  if (a) {
    x = cgeti(4);
    x[1] = evallgefint(4)|evalsigne(1);
    *(int_MSW(x)) = (long)a;
    *(int_LSW(x)) = (long)b;
  } else {
    x = cgeti(3);
    x[1] = evallgefint(3)|evalsigne(1);
    x[2] = (long)b;
  }
#endif
  return x;
}

/* return a_(n-1) x^(n-1) + ... + a_0 */
GEN
mkpoln(long n, ...)
{
  va_list ap;
  GEN x, y;
  long i;
  va_start(ap,n);
  x = cgetg(n+2, t_POL); y = x + 2;
  x[1] = evalvarn(0);
  for (i=n-1; i >= 0; i--) gel(y,i) = va_arg(ap, GEN);
  va_end(ap); return normalizepol(x);
}

/* return [a_1, ..., a_n] */
GEN
mkvecn(long n, ...)
{
  va_list ap;
  GEN x;
  long i;
  va_start(ap,n);
  x = cgetg(n+1, t_VEC);
  for (i=1; i <= n; i++) gel(x,i) = va_arg(ap, GEN);
  va_end(ap); return x;
}

GEN
mkcoln(long n, ...)
{
  va_list ap;
  GEN x;
  long i;
  va_start(ap,n);
  x = cgetg(n+1, t_COL);
  for (i=1; i <= n; i++) gel(x,i) = va_arg(ap, GEN);
  va_end(ap); return x;
}

GEN
scalarpol(GEN x, long v)
{
  GEN y;
  if (isexactzero(x)) return zeropol(v);
  y = cgetg(3,t_POL);
  y[1] = gcmp0(x)? evalvarn(v)
                 : evalvarn(v) | evalsigne(1);
  gel(y,2) = gcopy(x); return y;
}

/* deg1pol(a,b,x)=a*x+b, assumes a != 0 */
GEN
deg1pol(GEN x1, GEN x0,long v)
{
  GEN x = cgetg(4,t_POL);
  x[1] = evalsigne(1) | evalvarn(v);
  gel(x,2) = gcopy(x0);
  gel(x,3) = gcopy(x1); return normalizepol_i(x,4);
}

/* same, no copy */
GEN
deg1pol_i(GEN x1, GEN x0,long v)
{
  GEN x = cgetg(4,t_POL);
  x[1] = evalsigne(1) | evalvarn(v);
  gel(x,2) = x0;
  gel(x,3) = x1; return normalizepol_i(x,4);
}

static GEN
_gtopoly(GEN x, long v, int reverse)
{
  long tx=typ(x),lx,i,j;
  GEN y;

  if (v<0) v = 0;
  if (isexactzero(x)) return zeropol(v);
  if (is_scalar_t(tx)) return scalarpol(x,v);
  switch(tx)
  {
    case t_POL:
      if (varncmp(varn(x), v) < 0)
        pari_err(talker,"variable must have higher priority in gtopoly");
      y=gcopy(x); break;
    case t_SER:
      if (varncmp(varn(x), v) < 0)
        pari_err(talker,"variable must have higher priority in gtopoly");
      y = ser2rfrac(x);
      if (typ(y) != t_POL)
        pari_err(talker,"t_SER with negative valuation in gtopoly");
      break;
    case t_RFRAC:
      if (varncmp(varn(gel(x,2)), v) < 0)
        pari_err(talker,"variable must have higher priority in gtopoly");
      y=gdeuc(gel(x,1),gel(x,2)); break;
    case t_QFR: case t_QFI: case t_VEC: case t_COL: case t_MAT:
      lx = lg(x); if (tx == t_QFR) lx--;
      if (varncmp(gvar(x), v) <= 0)
        pari_err(talker,"variable must have higher priority in gtopoly");
      if (reverse)
      {
	while (lx-- && isexactzero(gel(x,lx)));
	i = lx+2; y = cgetg(i,t_POL);
	y[1] = gcmp0(x)? 0: evalsigne(1);
	for (j=2; j<i; j++) gel(y,j) = gcopy(gel(x,j-1));
      }
      else
      {
	i=1; j=lx; while (lx-- && isexactzero(gel(x,i++)));
	i = lx+2; y = cgetg(i,t_POL);
	y[1] = gcmp0(x)? 0: evalsigne(1);
	lx = j-1;
	for (j=2; j<i; j++) gel(y,j) = gcopy(gel(x,lx--));
      }
      break;
    default: pari_err(typeer,"gtopoly");
      return NULL; /* not reached */
  }
  setvarn(y,v); return y;
}

GEN
gtopolyrev(GEN x, long v) { return _gtopoly(x,v,1); }

GEN
gtopoly(GEN x, long v) { return _gtopoly(x,v,0); }

GEN
scalarser(GEN x, long v, long prec)
{
  long i, l;
  GEN y;                                                                       

  if (isexactzero(x)) return zeroser(v,0);                                     
  l = prec + 2; y = cgetg(l, t_SER);         
  y[1] = evalsigne(1) | evalvalp(0) | evalvarn(v);
  gel(y,2) = gcopy(x); for (i=3; i<l; i++) gel(y,i) = gen_0;
  return y;
}

static GEN _gtoser(GEN x, long v, long prec);

/* assume x a t_[SER|POL], apply gtoser to all coeffs */
static GEN
coefstoser(GEN x, long v, long prec)
{
  long i, tx = typ(x), lx = lg(x);
  GEN y = cgetg(lx, tx);
  y[1] = x[1];
  for (i=2; i<lx; i++) gel(y,i) = _gtoser(gel(x,i), v, prec);
  return y;
}

/* assume x a scalar or t_POL. Not stack-clean */
GEN
poltoser(GEN x, long v, long prec)
{
  long tx = typ(x), vx = varn(x);
  GEN y;

  if (is_scalar_t(tx) || varncmp(vx, v) > 0) return scalarser(x, v, prec);
  if (varncmp(vx, v) < 0) return coefstoser(x, v, prec);
  if (!lgpol(x)) return zeroser(v, prec);
  y = greffe(x, prec+2, 1);
  setvarn(y, v); return y;
}

/* x a t_RFRAC[N]. Not stack-clean */
GEN
rfractoser(GEN x, long v, long prec)
{
  return gdiv(poltoser(gel(x,1), v, prec), 
              poltoser(gel(x,2), v, prec));
}

GEN
toser_i(GEN x)
{
  switch(typ(x))
  {
    case t_SER: return x;
    case t_POL: return poltoser(x, varn(x), precdl);
    case t_RFRAC: return rfractoser(x, gvar(x), precdl);
  }
  return NULL;
}

static GEN
_gtoser(GEN x, long v, long prec)
{
  long tx=typ(x), lx, i, j, l;
  pari_sp av;
  GEN y;

  if (v < 0) v = 0;
  if (tx == t_SER)
  {
    long vx = varn(x);
    if      (varncmp(vx, v) < 0) y = coefstoser(x, v, prec);
    else if (varncmp(vx, v) > 0) y = scalarser(x, v, prec);
    else y = gcopy(x);
    return y;
  }
  if (is_scalar_t(tx)) return scalarser(x,v,prec);
  switch(tx)
  {
    case t_POL:
      if (varncmp(varn(x), v) < 0)
        pari_err(talker,"main variable must have higher priority in gtoser");
      y = poltoser(x, v, prec); l = lg(y);
      for (i=2; i<l; i++)
        if (gel(y,i) != gen_0) gel(y,i) = gcopy(gel(y,i));
      break;

    case t_RFRAC:
      if (varncmp(varn(gel(x,2)), v) < 0)
        pari_err(talker,"main variable must have higher priority in gtoser");
      av = avma;
      return gerepileupto(av, rfractoser(x, v, prec));

    case t_QFR: case t_QFI: case t_VEC: case t_COL:
      if (varncmp(gvar(x), v) < 0)
        pari_err(talker,"main variable must have higher priority in gtoser");
      lx = lg(x); if (tx == t_QFR) lx--;
      i = 1; while (i<lx && isexactzero(gel(x,i))) i++;
      if (i == lx) return zeroser(v, lx-1);
      lx -= i-2; x += i-2;
      y = cgetg(lx,t_SER);
      y[1] = evalsigne(1) | evalvalp(i-1) | evalvarn(v);
      for (j=2; j<lx; j++) gel(y,j) = gcopy(gel(x,j));
      break;

    default: pari_err(typeer,"gtoser");
      return NULL; /* not reached */
  }
  return y;
}

GEN
gtoser(GEN x, long v) { return _gtoser(x,v,precdl); }

/* assume typ(x) = t_STR */
static GEN
str_to_vecsmall(GEN x)
{
  char *s = GSTR(x);
  long i, l = strlen(s);
  GEN y = cgetg(l+1, t_VECSMALL);
  s--;
  for (i=1; i<=l; i++) y[i] = (long)s[i];
  return y;
}

GEN
gtovec(GEN x)
{
  long tx,lx,i;
  GEN y;

  if (!x) return cgetg(1,t_VEC);
  tx = typ(x);
  if (is_scalar_t(tx) || tx == t_RFRAC) return mkveccopy(x);
  if (tx == t_STR)
  {
    char t[2] = {0,0};
    y = str_to_vecsmall(x);
    lx = lg(y);
    for (i=1; i<lx; i++)
    {
      t[0] = y[i];
      gel(y,i) = strtoGENstr(t);
    }
    settyp(y,t_VEC); return y;
  }
  if (is_graphicvec_t(tx))
  {
    lx=lg(x); y=cgetg(lx,t_VEC);
    for (i=1; i<lx; i++) gel(y,i) = gcopy(gel(x,i));
    return y;
  }
  if (tx==t_POL)
  {
    lx=lg(x); y=cgetg(lx-1,t_VEC);
    for (i=1; i<=lx-2; i++) gel(y,i) = gcopy(gel(x,lx-i));
    return y;
  }
  if (tx==t_LIST)
  {
    lx=lgeflist(x); y=cgetg(lx-1,t_VEC); x++;
    for (i=1; i<=lx-2; i++) gel(y,i) = gcopy(gel(x,i));
    return y;
  }
  if (tx == t_VECSMALL) return vecsmall_to_vec(x);
  if (!signe(x)) return mkvec(gen_0);
  lx=lg(x); y=cgetg(lx-1,t_VEC); x++;
  for (i=1; i<=lx-2; i++) gel(y,i) = gcopy(gel(x,i));
  return y;
}

GEN
gtocol(GEN x)
{
  long lx, tx, i, j, h;
  GEN y;
  if (!x) return cgetg(1,t_COL);
  tx = typ(x);
  if (tx != t_MAT) { y = gtovec(x); settyp(y, t_COL); return y; }
  lx = lg(x); if (lx == 1) return cgetg(1, t_COL);
  h = lg(x[1]); y = cgetg(h, t_COL);
  for (i = 1 ; i < h; i++) {
    gel(y,i) = cgetg(lx, t_VEC);
    for (j = 1; j < lx; j++) gmael(y,i,j) = gcopy(gcoeff(x,i,j));
  }
  return y;
}

GEN
gtovecsmall(GEN x)
{
  GEN V;
  long tx, l,i;
  
  if (!x) return cgetg(1,t_VECSMALL);
  tx = typ(x);
  if (tx == t_VECSMALL) return gcopy(x);
  if (tx == t_INT) return mkvecsmall(itos(x));
  if (tx == t_STR) return str_to_vecsmall(x);
  if (!is_vec_t(tx)) pari_err(typeer,"vectosmall");
  l = lg(x);
  V = cgetg(l,t_VECSMALL);
  for(i=1;i<l;i++) V[i] = itos(gel(x,i));
  return V;
}

GEN
compo(GEN x, long n)
{
  long tx = typ(x);
  ulong l, lx = (ulong)lg(x);

  if (!is_recursive_t(tx))
    pari_err(talker, "this object is a leaf. It has no components");
  if (n < 1) pari_err(talker,"nonexistent component");
  if (tx == t_POL && (ulong)n+1 >= lx) return gen_0;
  if (tx == t_LIST) lx = (ulong)lgeflist(x);
  l = (ulong)lontyp[tx] + (ulong)n-1; /* beware overflow */
  if (l >= lx) pari_err(talker,"nonexistent component");
  return gcopy(gel(x,l));
}

/* assume v > varn(x), extract coeff of pol_x[v]^n */
static GEN
multi_coeff(GEN x, long n, long v, long dx)
{
  long i, lx = dx+3;
  GEN z = cgetg(lx, t_POL); z[1] = x[1];
  for (i = 2; i < lx; i++) gel(z,i) = polcoeff_i(gel(x,i), n, v);
  return normalizepol_i(z, lx);
}

/* assume x a t_POL */
static GEN
_polcoeff(GEN x, long n, long v)
{
  long w, dx;
  dx = degpol(x);
  if (dx < 0) return gen_0;
  if (v < 0 || v == (w=varn(x)))
    return (n < 0 || n > dx)? gen_0: gel(x,n+2);
  if (w > v) return n? gen_0: x;
  /* w < v */
  return multi_coeff(x, n, v, dx);
}

/* assume x a t_SER */
static GEN
_sercoeff(GEN x, long n, long v)
{
  long w, dx = degpol(x), ex = valp(x), N = n - ex;
  GEN z;
  if (dx < 0)
  {
    if (N >= 0) pari_err(talker,"non existent component in truecoeff");
    return gen_0;
  }
  if (v < 0 || v == (w=varn(x)))
  {
    if (N > dx) pari_err(talker,"non existent component in truecoeff");
    return (N < 0)? gen_0: gel(x,N+2);
  }
  if (w > v) return N? gen_0: x;
  /* w < v */
  z = multi_coeff(x, n, v, dx);
  if (ex) z = gmul(z, monomial(gen_1,ex, w));
  return z;
}

/* assume x a t_RFRAC(n) */
static GEN
_rfraccoeff(GEN x, long n, long v)
{
  GEN P,Q, p = gel(x,1), q = gel(x,2);
  long vp = gvar(p), vq = gvar(q);
  if (v < 0) v = min(vp, vq);
  P = (vp == v)? p: swap_vars(p, v);
  Q = (vq == v)? q: swap_vars(q, v);
  if (!ismonome(Q)) pari_err(typeer, "polcoeff");
  n += degpol(Q);
  return gdiv(_polcoeff(P, n, v), leading_term(Q));
}

GEN
polcoeff_i(GEN x, long n, long v)
{
  switch(typ(x))
  {
    case t_POL: return _polcoeff(x,n,v);
    case t_SER: return _sercoeff(x,n,v);
    case t_RFRAC: return _rfraccoeff(x,n,v);
    default: return n? gen_0: x;
  }
}

/* with respect to the main variable if v<0, with respect to the variable v
   otherwise. v ignored if x is not a polynomial/series. */
GEN
polcoeff0(GEN x, long n, long v)
{
  long tx=typ(x);
  pari_sp av;

  if (is_scalar_t(tx)) return n? gen_0: gcopy(x);

  av = avma;
  switch(tx)
  {
    case t_POL: x = _polcoeff(x,n,v); break;
    case t_SER: x = _sercoeff(x,n,v); break;
    case t_RFRAC: x = _rfraccoeff(x,n,v); break;
   
    case t_QFR: case t_QFI: case t_VEC: case t_COL: case t_MAT:
      if (n>=1 && n<lg(x)) return gcopy(gel(x,n));
    /* fall through */

    default: pari_err(talker,"nonexistent component in truecoeff");
  }
  if (x == gen_0) return x;
  if (avma == av) return gcopy(x);
  return gerepilecopy(av, x);
}

GEN
truecoeff(GEN x, long n)
{
  return polcoeff0(x,n,-1);
}

GEN
denom(GEN x)
{
  long lx, i;
  pari_sp av, tetpil;
  GEN s,t;

  switch(typ(x))
  {
    case t_INT: case t_REAL: case t_INTMOD: case t_PADIC: case t_SER:
      return gen_1;

    case t_FRAC:
      return icopy(gel(x,2));

    case t_COMPLEX:
      av=avma; t=denom(gel(x,1)); s=denom(gel(x,2)); tetpil=avma;
      return gerepile(av,tetpil,glcm(s,t));

    case t_QUAD:
      av=avma; t=denom(gel(x,2)); s=denom(gel(x,3)); tetpil=avma;
      return gerepile(av,tetpil,glcm(s,t));

    case t_POLMOD:
      return denom(gel(x,2));

    case t_RFRAC:
      return gcopy(gel(x,2));

    case t_POL:
      return pol_1[varn(x)];

    case t_VEC: case t_COL: case t_MAT:
      lx=lg(x); if (lx==1) return gen_1;
      av = tetpil = avma; s = denom(gel(x,1));
      for (i=2; i<lx; i++)
      {
        t = denom(gel(x,i));
        if (t != gen_1) { tetpil=avma; s=glcm(s,t); }
      }
      return gerepile(av,tetpil,s);
  }
  pari_err(typeer,"denom");
  return NULL; /* not reached */
}

GEN
numer(GEN x)
{
  pari_sp av, tetpil;
  GEN s;

  switch(typ(x))
  {
    case t_INT: case t_REAL: case t_INTMOD:
    case t_PADIC: case t_POL: case t_SER:
      return gcopy(x);

    case t_FRAC:
      return (signe(x[2])>0)? icopy(gel(x,1)): negi(gel(x,1));

    case t_POLMOD:
      av=avma; s=numer(gel(x,2)); tetpil=avma;
      return gerepile(av,tetpil,gmodulo(s,gel(x,1)));

    case t_RFRAC:
      return gcopy(gel(x,1));

    case t_COMPLEX: case t_QUAD: case t_VEC: case t_COL: case t_MAT:
      av=avma; s=denom(x); tetpil=avma;
      return gerepile(av,tetpil,gmul(s,x));
  }
  pari_err(typeer,"numer");
  return NULL; /* not reached */
}

/* Lift only intmods if v does not occur in x, lift with respect to main
 * variable of x if v < 0, with respect to variable v otherwise.
 */
GEN
lift0(GEN x, long v)
{
  long lx,tx=typ(x),i;
  GEN y;

  switch(tx)
  {
    case t_INT: case t_REAL:
      return gcopy(x);

    case t_INTMOD:
      return gcopy(gel(x,2));

    case t_POLMOD:
      if (v < 0 || v == varn(gel(x,1))) return gcopy(gel(x,2));
      y = cgetg(3,tx);
      gel(y,1) = lift0(gel(x,1),v);
      gel(y,2) = lift0(gel(x,2),v); return y;

    case t_PADIC:
      return gtrunc(x);

    case t_FRAC: case t_COMPLEX: case t_RFRAC:
    case t_POL: case t_SER: case t_VEC: case t_COL: case t_MAT:
      y = init_gen_op(x, tx, &lx, &i);
      for (; i<lx; i++) gel(y,i) = lift0(gel(x,i), v);
      return y;

    case t_QUAD:
      y=cgetg(4,t_QUAD); gel(y,1) = gcopy(gel(x,1));
      for (i=2; i<4; i++) gel(y,i) = lift0(gel(x,i), v);
      return y;
  }
  pari_err(typeer,"lift");
  return NULL; /* not reached */
}

GEN
lift(GEN x)
{
  return lift0(x,-1);
}

/* same as lift, without copy. May DESTROY x. For internal use only.
   Conventions on v as for lift. */
GEN
lift_intern0(GEN x, long v)
{
  long i,lx,tx=typ(x);

  switch(tx)
  {
    case t_INT: case t_REAL:
      return x;

    case t_INTMOD:
      return gel(x,2);

    case t_POLMOD:
      if (v < 0 || v == varn(gel(x,1))) return gel(x,2);
      gel(x,1) = lift_intern0(gel(x,1),v);
      gel(x,2) = lift_intern0(gel(x,2),v);
      return x;

    case t_SER:
    case t_FRAC: case t_COMPLEX: case t_QUAD: case t_POL:
    case t_RFRAC: case t_VEC: case t_COL: case t_MAT:
      lx = lg(x);
      for (i = lx-1; i>=lontyp[tx]; i--)
        gel(x,i) = lift_intern0(gel(x,i),v);
      return x;
  }
  pari_err(typeer,"lift_intern");
  return NULL; /* not reached */
}

/* memes conventions pour v que lift */
GEN
centerlift0(GEN x, long v)
{
  long i, lx, tx = typ(x);
  pari_sp av;
  GEN y;

  switch(tx)
  {
    case t_INT:
      return icopy(x);

    case t_INTMOD:
      av = avma; i = cmpii(shifti(gel(x,2),1), gel(x,1)); avma = av;
      return (i > 0)? subii(gel(x,2),gel(x,1)): icopy(gel(x,2));

    case t_POLMOD:
      if (v < 0 || v == varn(gel(x,1))) return gcopy(gel(x,2));
      y = cgetg(3, t_POLMOD);
      gel(y,1) = centerlift0(gel(x,1),v);
      gel(y,2) = centerlift0(gel(x,2),v); return y;

    case t_POL: case t_SER:
    case t_FRAC: case t_COMPLEX: case t_RFRAC:
    case t_VEC: case t_COL: case t_MAT:
      y = init_gen_op(x, tx, &lx, &i);
      for (; i<lx; i++) gel(y,i) = centerlift0(gel(x,i),v);
      return y;

    case t_QUAD:
      y=cgetg(4,t_QUAD); gel(y,1) = gcopy(gel(x,1));
      gel(y,2) = centerlift0(gel(x,2),v);
      gel(y,3) = centerlift0(gel(x,3),v); return y;
  }
  pari_err(typeer,"centerlift");
  return NULL; /* not reached */
}

GEN
centerlift(GEN x)
{
  return centerlift0(x,-1);
}

/*******************************************************************/
/*                                                                 */
/*                  PARTIES REELLE ET IMAGINAIRES                  */
/*                                                                 */
/*******************************************************************/

static GEN
op_ReIm(GEN f(GEN), GEN x)
{
  long lx, i, j, tx = typ(x);
  pari_sp av;
  GEN z;

  switch(tx)
  {
    case t_POL:
      lx = lg(x); z = cgetg(lx,t_POL); z[1] = x[1];
      for (j=2; j<lx; j++) gel(z,j) = f(gel(x,j));
      return normalizepol_i(z, lx);

    case t_SER:
      lx = lg(x); z = cgetg(lx,t_SER); z[1] = x[1];
      for (j=2; j<lx; j++) gel(z,j) = f(gel(x,j));
      return normalize(z);

    case t_RFRAC: {
      GEN dxb, n, d;
      av = avma; dxb = gconj(gel(x,2));
      n = gmul(gel(x,1), dxb);
      d = gmul(gel(x,2), dxb);
      return gerepileupto(av, gdiv(f(n), d));
    }

    case t_VEC: case t_COL: case t_MAT:
      lx=lg(x); z=cgetg(lx,tx);
      for (i=1; i<lx; i++) gel(z,i) = f(gel(x,i));
      return z;
  }
  pari_err(typeer,"greal/gimag");
  return NULL; /* not reached */
}

GEN
real_i(GEN x)
{
  switch(typ(x))
  {
    case t_INT: case t_REAL: case t_FRAC:
      return x;
    case t_COMPLEX:
      return gel(x,1);
    case t_QUAD:
      return gel(x,2);
  }
  return op_ReIm(real_i,x);
}
GEN
imag_i(GEN x)
{
  switch(typ(x))
  {
    case t_INT: case t_REAL: case t_FRAC:
      return gen_0;
    case t_COMPLEX:
      return gel(x,2);
    case t_QUAD:
      return gel(x,3);
  }
  return op_ReIm(imag_i,x);
}
GEN
greal(GEN x)
{
  switch(typ(x))
  {
    case t_INT: case t_REAL: case t_FRAC:
      return gcopy(x);

    case t_COMPLEX:
      return gcopy(gel(x,1));

    case t_QUAD:
      return gcopy(gel(x,2));
  }
  return op_ReIm(greal,x);
}
GEN
gimag(GEN x)
{
  switch(typ(x))
  {
    case t_INT: case t_REAL: case t_FRAC:
      return gen_0;

    case t_COMPLEX:
      return gcopy(gel(x,2));

    case t_QUAD:
      return gcopy(gel(x,3));
  }
  return op_ReIm(gimag,x);
}

/*******************************************************************/
/*                                                                 */
/*                     LOGICAL OPERATIONS                          */
/*                                                                 */
/*******************************************************************/
static long
_egal(GEN x, GEN y)
{
  pari_sp av = avma;
  long r = gequal(simplify_i(x), simplify_i(y));
  avma = av; return r;
}

GEN
glt(GEN x, GEN y) { return gcmp(x,y)<0? gen_1: gen_0; }

GEN
gle(GEN x, GEN y) { return gcmp(x,y)<=0? gen_1: gen_0; }

GEN
gge(GEN x, GEN y) { return gcmp(x,y)>=0? gen_1: gen_0; }

GEN
ggt(GEN x, GEN y) { return gcmp(x,y)>0? gen_1: gen_0; }

GEN
geq(GEN x, GEN y) { return _egal(x,y)? gen_1: gen_0; }

GEN
gne(GEN x, GEN y) { return _egal(x,y)? gen_0: gen_1; }

GEN
gand(GEN x, GEN y) { return gcmp0(x)? gen_0: (gcmp0(y)? gen_0: gen_1); }

GEN
gor(GEN x, GEN y) { return gcmp0(x)? (gcmp0(y)? gen_0: gen_1): gen_1; }

GEN
gnot(GEN x) { return gcmp0(x)? gen_1: gen_0; }

/*******************************************************************/
/*                                                                 */
/*                      FORMAL SIMPLIFICATIONS                     */
/*                                                                 */
/*******************************************************************/

GEN
geval(GEN x)
{
  long lx, i, tx = typ(x);
  pari_sp av, tetpil;
  GEN y,z;

  if (is_const_t(tx)) return gcopy(x);
  if (is_graphicvec_t(tx))
  {
    lx=lg(x); y=cgetg(lx, tx);
    for (i=1; i<lx; i++) gel(y,i) = geval(gel(x,i));
    return y;
  }

  switch(tx)
  {
    case t_STR:
      return gp_read_str(GSTR(x));

    case t_POLMOD: y=cgetg(3,tx);
      gel(y,1) = geval(gel(x,1));
      av=avma; z=geval(gel(x,2)); tetpil=avma;
      gel(y,2) = gerepile(av,tetpil,gmod(z,gel(y,1)));
      return y;

    case t_POL:
      lx=lg(x); if (lx==2) return gen_0;
      {
        long vx = varn(x);
        entree *ep = varentries[vx];
        if (!ep) return gcopy(x);
        z = (GEN)ep->value;
        if (gequal(x, pol_x[vx])) return gcopy(z);
      }
      y=gen_0; av=avma;
      for (i=lx-1; i>1; i--)
        y = gadd(geval(gel(x,i)), gmul(z,y));
      return gerepileupto(av, y);

    case t_SER:
      pari_err(impl, "evaluation of a power series");

    case t_RFRAC:
      return gdiv(geval(gel(x,1)),geval(gel(x,2)));
  }
  pari_err(typeer,"geval");
  return NULL; /* not reached */
}

GEN
simplify_i(GEN x)
{
  long tx=typ(x),i,lx;
  GEN y;

  switch(tx)
  {
    case t_INT: case t_REAL: case t_FRAC:
    case t_INTMOD: case t_PADIC: case t_QFR: case t_QFI:
    case t_LIST: case t_STR: case t_VECSMALL:
      return x;

    case t_COMPLEX:
      if (isexactzero(gel(x,2))) return simplify_i(gel(x,1));
      y=cgetg(3,t_COMPLEX);
      gel(y,1) = simplify_i(gel(x,1));
      gel(y,2) = simplify_i(gel(x,2)); return y;

    case t_QUAD:
      if (isexactzero(gel(x,3))) return simplify_i(gel(x,2));
      y=cgetg(4,t_QUAD);
      y[1]=x[1];
      gel(y,2) = simplify_i(gel(x,2));
      gel(y,3) = simplify_i(gel(x,3)); return y;

    case t_POLMOD: y=cgetg(3,t_POLMOD);
      gel(y,1) = simplify_i(gel(x,1));
      if (typ(y[1]) != t_POL) y[1] = x[1]; /* invalid object otherwise */
      gel(y,2) = simplify_i(gel(x,2)); return y;

    case t_POL:
      lx = lg(x); if (lx==2) return gen_0;
      if (lx==3) return simplify_i(gel(x,2));
      y = cgetg(lx,t_POL); y[1] = x[1];
      for (i=2; i<lx; i++) gel(y,i) = simplify_i(gel(x,i));
      return y;

    case t_SER:
      lx = lg(x); y = cgetg(lx,t_SER); y[1] = x[1];
      for (i=2; i<lx; i++) gel(y,i) = simplify_i(gel(x,i));
      return y;

    case t_RFRAC: y=cgetg(3,t_RFRAC);
      gel(y,1) = simplify_i(gel(x,1));
      gel(y,2) = simplify_i(gel(x,2));
      if (typ(gel(y,2)) != t_POL) return gdiv(gel(y,1), gel(y,2));
      return y;

    case t_VEC: case t_COL: case t_MAT:
      lx=lg(x); y=cgetg(lx,tx);
      for (i=1; i<lx; i++) gel(y,i) = simplify_i(gel(x,i));
      return y;
  }
  pari_err(typeer,"simplify_i");
  return NULL; /* not reached */
}

GEN
simplify(GEN x)
{
  pari_sp av = avma;
  return gerepilecopy(av, simplify_i(x));
}

/*******************************************************************/
/*                                                                 */
/*                EVALUATION OF SOME SIMPLE OBJECTS                */
/*                                                                 */
/*******************************************************************/
static GEN
qfeval0_i(GEN q, GEN x, long n)
{
  long i, j;
  pari_sp av=avma;
  GEN res=gen_0;

  for (i=2;i<n;i++)
    for (j=1;j<i;j++)
      res = gadd(res, gmul(gcoeff(q,i,j), mulii(gel(x,i),gel(x,j))) );
  res=gshift(res,1);
  for (i=1;i<n;i++)
    res = gadd(res, gmul(gcoeff(q,i,i), sqri(gel(x,i))) );
  return gerepileupto(av,res);
}

#if 0
static GEN
qfeval0(GEN q, GEN x, long n)
{
  long i, j;
  pari_sp av=avma;
  GEN res=gen_0;

  for (i=2;i<n;i++)
    for (j=1;j<i;j++)
      res = gadd(res, gmul(gcoeff(q,i,j), gmul(gel(x,i),gel(x,j))) );
  res=gshift(res,1);
  for (i=1;i<n;i++)
    res = gadd(res, gmul(gcoeff(q,i,i), gsqr(gel(x,i))) );
  return gerepileupto(av,res);
}
#else
static GEN
qfeval0(GEN q, GEN x, long n)
{
  long i, j;
  pari_sp av=avma;
  GEN res = gmul(gcoeff(q,1,1), gsqr(gel(x,1)));

  for (i=2; i<n; i++)
  {
    GEN l = gel(q,i);
    GEN sx = gmul(gel(l,1), gel(x,1));
    for (j=2; j<i; j++)
      sx = gadd(sx, gmul(gel(l,j),gel(x,j)));
    sx = gadd(gshift(sx,1), gmul(gel(l,i),gel(x,i)));

    res = gadd(res, gmul(gel(x,i), sx));
  }
  return gerepileupto(av,res);
}
#endif

/* We assume q is a real symetric matrix */
GEN
qfeval(GEN q, GEN x)
{
  long n=lg(q);

  if (n==1)
  {
    if (typ(q) != t_MAT || lg(x) != 1)
      pari_err(talker,"invalid data in qfeval");
    return gen_0;
  }
  if (typ(q) != t_MAT || lg(q[1]) != n)
    pari_err(talker,"invalid quadratic form in qfeval");
  if (typ(x) != t_COL || lg(x) != n)
    pari_err(talker,"invalid vector in qfeval");

  return qfeval0(q,x,n);
}

/* the Horner-type evaluation (mul x 2/3) would force us to use gmul and not
 * mulii (more than 4 x slower for small entries). Not worth it.
 */
static GEN
qfbeval0_i(GEN q, GEN x, GEN y, long n)
{
  long i, j;
  pari_sp av=avma;
  GEN res = gmul(gcoeff(q,1,1), mulii(gel(x,1),gel(y,1)));

  for (i=2;i<n;i++)
  {
    if (!signe(x[i]))
    {
      if (!signe(y[i])) continue;
      for (j=1;j<i;j++)
        res = gadd(res, gmul(gcoeff(q,i,j), mulii(gel(x,j),gel(y,i))));
    }
    else if (!signe(y[i]))
    {
      for (j=1;j<i;j++)
        res = gadd(res, gmul(gcoeff(q,i,j), mulii(gel(x,i),gel(y,j))));
    }
    else
    {
      for (j=1;j<i;j++)
      {
        GEN p1 = addii(mulii(gel(x,i),gel(y,j)), mulii(gel(x,j),gel(y,i)));
        res = gadd(res, gmul(gcoeff(q,i,j),p1));
      }
      res = gadd(res, gmul(gcoeff(q,i,i), mulii(gel(x,i),gel(y,i))));
    }
  }
  return gerepileupto(av,res);
}

#if 0
static GEN
qfbeval0(GEN q, GEN x, GEN y, long n)
{
  long i, j;
  pari_sp av=avma;
  GEN res = gmul(gcoeff(q,1,1), gmul(gel(x,1),gel(y,1)));

  for (i=2;i<n;i++)
  {
    for (j=1;j<i;j++)
    {
      GEN p1 = gadd(gmul(gel(x,i),gel(y,j)), gmul(gel(x,j),gel(y,i)));
      res = gadd(res, gmul(gcoeff(q,i,j),p1));
    }
    res = gadd(res, gmul(gcoeff(q,i,i), gmul(gel(x,i),gel(y,i))));
  }
  return gerepileupto(av,res);
}
#else
static GEN
qfbeval0(GEN q, GEN x, GEN y, long n)
{
  long i, j;
  pari_sp av=avma;
  GEN res = gmul(gcoeff(q,1,1), gmul(gel(x,1), gel(y,1)));

  for (i=2; i<n; i++)
  {
    GEN l = gel(q,i);
    GEN sx = gmul(gel(l,1), gel(y,1));
    GEN sy = gmul(gel(l,1), gel(x,1));
    for (j=2; j<i; j++)
    {
      sx = gadd(sx, gmul(gel(l,j),gel(y,j)));
      sy = gadd(sy, gmul(gel(l,j),gel(x,j)));
    }
    sx = gadd(sx, gmul(gel(l,i),gel(y,i)));

    sx = gmul(gel(x,i), sx);
    sy = gmul(gel(y,i), sy);
    res = gadd(res, gadd(sx,sy));
  }
  return gerepileupto(av,res);
}
#endif

/* We assume q is a real symetric matrix */
GEN
qfbeval(GEN q, GEN x, GEN y)
{
  long n=lg(q);

  if (n==1)
  {
    if (typ(q) != t_MAT || lg(x) != 1 || lg(y) != 1)
      pari_err(talker,"invalid data in qfbeval");
    return gen_0;
  }
  if (typ(q) != t_MAT || lg(q[1]) != n)
    pari_err(talker,"invalid quadratic form in qfbeval");
  if (typ(x) != t_COL || lg(x) != n || typ(y) != t_COL || lg(y) != n)
    pari_err(talker,"invalid vector in qfbeval");

  return qfbeval0(q,x,y,n);
}

/* yield X = M'.q.M, assuming q is symetric.
 * X_ij are X_ji identical, not copies
 * if flag is set, M has integer entries
 */
GEN
qf_base_change(GEN q, GEN M, int flag)
{
  long i,j, k = lg(M), n=lg(q);
  GEN res = cgetg(k,t_MAT);
  GEN (*qf)(GEN,GEN,long)  = flag? &qfeval0_i:  &qfeval0;
  GEN (*qfb)(GEN,GEN,GEN,long) = flag? &qfbeval0_i: &qfbeval0;

  if (n==1)
  {
    if (typ(q) != t_MAT || k != 1)
      pari_err(talker,"invalid data in qf_base_change");
    return res;
  }
  if (typ(M) != t_MAT || k == 1 || lg(M[1]) != n)
    pari_err(talker,"invalid base change matrix in qf_base_change");

  for (i=1;i<k;i++)
  {
    gel(res,i) = cgetg(k,t_COL);
    gcoeff(res,i,i) = qf(q,gel(M,i),n);
  }
  for (i=2;i<k;i++)
    for (j=1;j<i;j++)
      gcoeff(res,i,j)=gcoeff(res,j,i) = qfb(q,gel(M,i),gel(M,j),n);
  return res;
}

/* return Re(x * y), x and y scalars */
GEN
mul_real(GEN x, GEN y)
{
  if (typ(x) == t_COMPLEX)
  {
    if (typ(y) == t_COMPLEX)
    {
      pari_sp av=avma, tetpil;
      GEN p1 = gmul(gel(x,1), gel(y,1));
      GEN p2 = gneg(gmul(gel(x,2), gel(y,2)));
      tetpil=avma; return gerepile(av,tetpil,gadd(p1,p2));
    }
    x = gel(x,1);
  }
  else if (typ(y) == t_COMPLEX) y = gel(y,1);
  return gmul(x,y);
}

/* Compute Re(x * y), x and y matrices of compatible dimensions
 * assume lx, ly > 1, and scalar entries */
GEN
mulmat_real(GEN x, GEN y)
{
  long i, j, k, lx = lg(x), ly = lg(y), l = lg(x[1]);
  pari_sp av;
  GEN p1, z = cgetg(ly,t_MAT);

  for (j=1; j<ly; j++)
  {
    gel(z,j) = cgetg(l,t_COL);
    for (i=1; i<l; i++)
    {
      p1 = gen_0; av=avma;
      for (k=1; k<lx; k++)
        p1 = gadd(p1, mul_real(gcoeff(x,i,k),gcoeff(y,k,j)));
      gcoeff(z,i,j) = gerepileupto(av, p1);
    }
  }
  return z;
}

static GEN
hqfeval0(GEN q, GEN x, long n)
{
  long i, j;
  pari_sp av=avma;
  GEN res=gen_0;

  for (i=2;i<n;i++)
    for (j=1;j<i;j++)
    {
      GEN p1 = gmul(gel(x,i),gconj(gel(x,j)));
      res = gadd(res, mul_real(gcoeff(q,i,j),p1));
    }
  res=gshift(res,1);
  for (i=1;i<n;i++)
    res = gadd(res, gmul(gcoeff(q,i,i), gnorm(gel(x,i))) );
  return gerepileupto(av,res);
}

/* We assume q is a hermitian complex matrix */
GEN
hqfeval(GEN q, GEN x)
{
  long n=lg(q);

  if (n==1)
  {
    if (typ(q) != t_MAT || lg(x) != 1)
      pari_err(talker,"invalid data in hqfeval");
    return gen_0;
  }
  if (typ(q) != t_MAT || lg(q[1]) != n)
    pari_err(talker,"invalid quadratic form in hqfeval");
  if (typ(x) != t_COL || lg(x) != n)
    pari_err(talker,"invalid vector in hqfeval");

  return hqfeval0(q,x,n);
}

GEN
poleval(GEN x, GEN y)
{
  long i, j, imin, tx = typ(x);
  pari_sp av0 = avma, av, lim;
  GEN p1, p2, r, s;

  if (is_scalar_t(tx)) return gcopy(x);
  switch(tx)
  {
    case t_POL:
      i = lg(x)-1; imin = 2; break;

    case t_RFRAC:
      p1 = poleval(gel(x,1),y);
      p2 = poleval(gel(x,2),y);
      return gerepileupto(av0, gdiv(p1,p2));

    case t_VEC: case t_COL:
      i = lg(x)-1; imin = 1; break;
    default: pari_err(typeer,"poleval");
      return NULL; /* not reached */
  }
  if (i<=imin)
    return (i==imin)? gcopy(gel(x,imin)): gen_0;

  lim = stack_lim(av0,2);
  p1 = gel(x,i); i--;
  if (typ(y)!=t_COMPLEX)
  {
#if 0 /* standard Horner's rule */
    for ( ; i>=imin; i--)
      p1 = gadd(gmul(p1,y),gel(x,i));
#endif
    /* specific attention to sparse polynomials */
    for ( ; i>=imin; i=j-1)
    {
      for (j=i; isexactzero(gel(x,j)); j--)
        if (j==imin)
        {
          if (i!=j) y = gpowgs(y, i-j+1);
          return gerepileupto(av0, gmul(p1,y));
        }
      r = (i==j)? y: gpowgs(y, i-j+1);
      p1 = gadd(gmul(p1,r), gel(x,j));
      if (low_stack(lim, stack_lim(av0,2)))
      {
        if (DEBUGMEM>1) pari_warn(warnmem,"poleval: i = %ld",i);
        p1 = gerepileupto(av0, p1);
      }
    }
    return gerepileupto(av0,p1);
  }

  p2 = gel(x,i); i--; r = gtrace(y); s = gneg_i(gnorm(y));
  av = avma;
  for ( ; i>=imin; i--)
  {
    GEN p3 = gadd(p2, gmul(r, p1));
    p2 = gadd(gel(x,i), gmul(s, p1)); p1 = p3;
    if (low_stack(lim, stack_lim(av0,2)))
    {
      if (DEBUGMEM>1) pari_warn(warnmem,"poleval: i = %ld",i);
      gerepileall(av, 2, &p1, &p2);
    }
  }
  return gerepileupto(av0, gadd(p2, gmul(y,p1)));
}
