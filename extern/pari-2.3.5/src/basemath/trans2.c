/* $Id: trans2.c 10281 2008-06-09 11:10:14Z kb $

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
/**                   TRANSCENDENTAL FUNCTIONS                     **/
/**                          (part 2)                              **/
/**                                                                **/
/********************************************************************/
#include "pari.h"
#include "paripriv.h"

/********************************************************************/
/**                                                                **/
/**                          ARCTANGENT                            **/
/**                                                                **/
/********************************************************************/
static GEN
mpatan(GEN x)
{
  long l, l1, l2, n, m, i, lp, e, s, sx = signe(x);
  pari_sp av0, av;
  double alpha, beta, delta;
  GEN y, p1, p2, p3, p4, p5, unr;
  int inv;

  if (!sx) return real_0_bit(expo(x));
  l = lp = lg(x);
  if (absrnz_egal1(x)) { /* |x| = 1 */
    y = Pi2n(-2, l+1); if (sx < 0) setsigne(y,-1);
    return y;
  }
  if (l > AGM_ATAN_LIMIT)
  {
    av = avma;
    return gerepileuptoleaf(av, (GEN)logagmcx(mkcomplex(gen_1, x), l)[2]);
  }

  e = expo(x); inv = (e >= 0); /* = (|x| > 1 ) */
  if (e > 0) lp += (e>>TWOPOTBITS_IN_LONG);

  y = cgetr(lp); av0 = avma;
  p1 = cgetr(l+1); affrr(x,p1); setsigne(p1, 1); /* p1 = |x| */
  if (inv) p1 = divsr(1, p1);
  e = expo(p1);
  if (e < -100)
    alpha = 1.65149612947 - e; /* log_2(Pi) - e */
  else
    alpha = log2(PI / atan(rtodbl(p1)));
  beta = (double)(bit_accuracy(l)>>1);
  delta = 1 + beta - alpha/2;
  if (delta <= 0) { n = 1; m = 0; }
  else
  {
    double fi = alpha-2;
#if 0
    const double gama = 1.; /* optimize this */
    if (delta >= gama*fi*fi)
    {
      n = (long)(1+sqrt(gama*delta));
      m = (long)(1+sqrt(delta/gama) - fi);
    }
#else
    if (delta >= fi*fi)
    {
      double t = 1 + sqrt(delta);
      n = (long)t;
      m = (long)(t - fi);
    }
#endif
    else
    {
      n = (long)(1+beta/fi);
      m = 0;
    }
  }
  l2 = l+1+(m>>TWOPOTBITS_IN_LONG);
  p2 = cgetr(l2); affrr(p1,p2); av = avma;
  for (i=1; i<=m; i++)
  {
    p5 = addsr(1, mulrr(p2,p2)); setlg(p5,l2);
    p5 = addsr(1, sqrtr_abs(p5)); setlg(p5,l2);
    affrr(divrr(p2,p5), p2); avma = av;
  }
  p3 = mulrr(p2,p2); l1 = 4;
  unr = real_1(l2); setlg(unr,4);
  p4 = cgetr(l2); setlg(p4,4);
  affrr(divrs(unr,2*n+1), p4);
  s = 0; e = expo(p3); av = avma;
  for (i = n; i > 1; i--) /* n >= 1. i = 1 done outside for efficiency */
  {
    setlg(p3,l1); p5 = mulrr(p4,p3);
    s -= e; l1 += (s>>TWOPOTBITS_IN_LONG);
    s %= BITS_IN_LONG;
    if (l1 > l2) l1 = l2;
    setlg(unr,l1); p5 = subrr(divrs(unr,2*i-1), p5);
    setlg(p4,l1); affrr(p5,p4); avma = av;
  }
  setlg(p3, l2); p5 = mulrr(p4,p3); /* i = 1 */
  setlg(unr,l2); p4 = subrr(unr, p5);

  p4 = mulrr(p2,p4); setexpo(p4, expo(p4)+m);
  if (inv) p4 = subrr(Pi2n(-1, lp), p4);
  if (sx < 0) setsigne(p4,-signe(p4));
  affr_fixlg(p4,y); avma = av0; return y;
}

GEN
gatan(GEN x, long prec)
{
  pari_sp av;
  GEN a, y;

  switch(typ(x))
  {
    case t_REAL:
      return mpatan(x);

    case t_COMPLEX:
      av = avma; return gerepilecopy(av, mulcxmI(gath(mulcxI(x),prec)));

    case t_INTMOD: case t_PADIC: pari_err(typeer,"gatan");

    default:
      av = avma; if (!(y = toser_i(x))) break;
      if (valp(y) < 0) pari_err(negexper,"gatan");
      if (lg(y)==2) return gcopy(y);
      /* lg(y) > 2 */
      a = integ(gdiv(derivser(y), gaddsg(1,gsqr(y))), varn(y));
      if (!valp(y)) a = gadd(a, gatan(gel(y,2),prec));
      return gerepileupto(av, a);
  }
  return transc(gatan,x,prec);
}
/********************************************************************/
/**                                                                **/
/**                             ARCSINE                            **/
/**                                                                **/
/********************************************************************/
/* |x| < 1, x != 0 */
static GEN
mpasin(GEN x) {
  pari_sp av = avma;
  GEN z, a = sqrtr(subsr(1, mulrr(x,x)));
  if (lg(x) > AGM_ATAN_LIMIT)
    z = (GEN)logagmcx(mkcomplex(a,x), lg(x))[2];
  else
    z = mpatan(divrr(x, a));
  return gerepileuptoleaf(av, z);
}

static GEN mpach(GEN x);
GEN
gasin(GEN x, long prec)
{
  long sx;
  pari_sp av;
  GEN a, y, p1;

  switch(typ(x))
  {
    case t_REAL: sx = signe(x);
      if (!sx) return real_0_bit(expo(x));
      if (absrnz_egal1(x)) { /* |x| = 1 */
        if (sx > 0) return Pi2n(-1, lg(x)); /* 1 */
        y = Pi2n(-1, lg(x)); setsigne(y, -1); return y; /* -1 */
      }
      if (expo(x) < 0) return mpasin(x);
      y = cgetg(3,t_COMPLEX);
      gel(y,1) = Pi2n(-1, lg(x));
      gel(y,2) = mpach(x);
      if (sx < 0)
      {
        setsigne(y[1],-signe(y[1]));
        setsigne(y[2],-signe(y[2]));
      }
      return y;

    case t_COMPLEX:
      av = avma;
      return gerepilecopy(av, mulcxmI(gash(mulcxI(x), prec)));

    case t_INTMOD: case t_PADIC: pari_err(typeer,"gasin");
    default:
      av = avma; if (!(y = toser_i(x))) break;
      if (gcmp0(y)) return gcopy(y);
      /* lg(y) > 2*/
      if (valp(y) < 0) pari_err(negexper,"gasin");
      p1 = gsubsg(1,gsqr(y));
      if (gcmp0(p1))
      {
        GEN t = Pi2n(-1,prec);
        if (gsigne(gel(y,2)) < 0) setsigne(t, -1);
        return gerepileupto(av, scalarser(t, varn(y), valp(p1)>>1));
      }
      p1 = gdiv(derivser(y), gsqrt(p1,prec));
      a = integ(p1,varn(y));
      if (!valp(y)) a = gadd(a, gasin(gel(y,2),prec));
      return gerepileupto(av, a);
  }
  return transc(gasin,x,prec);
}
/********************************************************************/
/**                                                                **/
/**                             ARCCOSINE                          **/
/**                                                                **/
/********************************************************************/
static GEN
acos0(long e) {
  long l = e >> TWOPOTBITS_IN_LONG; if (l >= 0) l = -1;
  return Pi2n(-1, 2-l);
}

/* |x| < 1, x != 0 */
static GEN
mpacos(GEN x)
{
  pari_sp av = avma;
  GEN z, a = sqrtr(subsr(1, mulrr(x,x)));
  if (lg(x) > AGM_ATAN_LIMIT)
    z = (GEN)logagmcx(mkcomplex(x,a), lg(x))[2];
  else {
    z = mpatan(divrr(a, x));
    if (signe(x) < 0) z = addrr(mppi(lg(z)), z);
  }
  return gerepileuptoleaf(av, z);
}

GEN
gacos(GEN x, long prec)
{
  long sx;
  pari_sp av;
  GEN a, y, p1;

  switch(typ(x))
  {
    case t_REAL: sx = signe(x);
      if (!sx) return acos0(expo(x));
      if (absrnz_egal1(x)) /* |x| = 1 */
        return sx > 0? real_0_bit( -(bit_accuracy(lg(x))>>1) ) : mppi(lg(x));
      if (expo(x) < 0) return mpacos(x);

      y = cgetg(3,t_COMPLEX); p1 = mpach(x);
      if (sx < 0) gel(y,1) = mppi(lg(x));
      else {
	gel(y,1) = gen_0;
        setsigne(p1,-signe(p1));
      }
      gel(y,2) = p1; return y;

    case t_COMPLEX: av = avma;
      return gerepilecopy(av, mulcxmI(gach(x,prec)));

    case t_INTMOD: case t_PADIC: pari_err(typeer,"gacos");
    case t_SER:
      av = avma; if (!(y = toser_i(x))) break;
      if (valp(y) < 0) pari_err(negexper,"gacos");
      if (lg(y) > 2)
      {
	p1 = gsubsg(1,gsqr(y));
	if (gcmp0(p1)) return zeroser(varn(y), valp(p1)>>1);
	p1 = integ(gdiv(gneg(derivser(y)), gsqrt(p1,prec)), varn(y));
	if (gcmp1(gel(y,2)) && !valp(y)) /*y = 1+O(y^k), k>=1*/
	  return gerepileupto(av, p1);
      }
      else p1 = y;
      a = (lg(y)==2 || valp(y))? Pi2n(-1, prec): gacos(gel(y,2),prec);
      return gerepileupto(av, gadd(a,p1));
  }
  return transc(gacos,x,prec);
}
/********************************************************************/
/**                                                                **/
/**                            ARGUMENT                            **/
/**                                                                **/
/********************************************************************/

/* we know that x and y are not both 0 */
static GEN
mparg(GEN x, GEN y)
{
  long prec, sx = signe(x), sy = signe(y);
  GEN z;

  if (!sy)
  {
    if (sx > 0) return real_0_bit(expo(y) - expo(x));
    return mppi(lg(x));
  }
  prec = lg(y); if (prec < lg(x)) prec = lg(x);
  if (!sx)
  {
    z = Pi2n(-1, prec); if (sy < 0) setsigne(z,-1);
    return z;
  }

  if (expo(x)-expo(y) > -2)
  {
    z = mpatan(divrr(y,x)); if (sx > 0) return z;
    return addrr_sign(z, signe(z), mppi(prec), sy);
  }
  z = mpatan(divrr(x,y));
  return addrr_sign(z, -signe(z), Pi2n(-1, prec), sy);
}

static GEN
rfix(GEN x,long prec)
{
  switch(typ(x))
  {
    case t_INT: return itor(x, prec);
    case t_FRAC: return rdivii(gel(x,1),gel(x,2), prec);
    case t_REAL: break;
    default: pari_err(typeer,"rfix (conversion to t_REAL)");
  }
  return x;
}

static GEN
cxarg(GEN x, GEN y, long prec)
{
  pari_sp av = avma;
  x = rfix(x,prec);
  y = rfix(y,prec); return gerepileuptoleaf(av, mparg(x,y));
}

GEN
garg(GEN x, long prec)
{
  long tx = typ(x);
  pari_sp av;

  if (gcmp0(x)) pari_err(talker,"zero argument in garg");
  switch(tx)
  {
    case t_REAL: prec = lg(x); /* fall through */
    case t_INT: case t_FRAC:
      return (gsigne(x)>0)? real_0(prec): mppi(prec);

    case t_QUAD:
      av = avma;
      return gerepileuptoleaf(av, garg(quadtoc(x, prec), prec));

    case t_COMPLEX:
      return cxarg(gel(x,1),gel(x,2),prec);

    case t_VEC: case t_COL: case t_MAT:
      return transc(garg,x,prec);
  }
  pari_err(typeer,"garg");
  return NULL; /* not reached */
}

/********************************************************************/
/**                                                                **/
/**                      HYPERBOLIC COSINE                         **/
/**                                                                **/
/********************************************************************/

static GEN
mpch(GEN x)
{
  pari_sp av;
  GEN z;

  if (gcmp0(x)) { /* 1 + x */
    long e = expo(x);
    if (e > 0) return real_0_bit(e);
    return real_1(3 + ((-e)>>TWOPOTBITS_IN_LONG));
  }
  av = avma;
  z = mpexp(x); z = addrr(z, ginv(z)); setexpo(z, expo(z)-1);
  return gerepileuptoleaf(av, z);
}

GEN
gch(GEN x, long prec)
{
  pari_sp av;
  GEN y, p1;

  switch(typ(x))
  {
    case t_REAL: return mpch(x);
    case t_COMPLEX: case t_PADIC: 
      av = avma; p1 = gexp(x,prec); p1 = gadd(p1, ginv(p1));
      return gerepileupto(av, gmul2n(p1,-1));
    case t_INTMOD: pari_err(typeer,"gch");
    default:
      av = avma; if (!(y = toser_i(x))) break;
      if (gcmp0(y) && valp(y) == 0) return gcopy(y);
      p1 = gexp(y,prec); p1 = gadd(p1, ginv(p1));
      return gerepileupto(av, gmul2n(p1,-1));
  }
  return transc(gch,x,prec);
}
/********************************************************************/
/**                                                                **/
/**                       HYPERBOLIC SINE                          **/
/**                                                                **/
/********************************************************************/

static GEN
mpsh(GEN x)
{
  pari_sp av;
  long ex = expo(x), lx;
  GEN z, res;

  if (!signe(x)) return real_0_bit(ex);
  lx = lg(x); res = cgetr(lx); av = avma;
  if (ex < 1 - BITS_IN_LONG) x = rtor(x, lx + nbits2nlong(-ex)-1);
  z = mpexp(x); z = addrr(z, divsr(-1,z)); setexpo(z, expo(z)-1);
  affrr(z, res); avma = av; return res;
}

GEN
gsh(GEN x, long prec)
{
  pari_sp av;
  GEN y, p1;

  switch(typ(x))
  {
    case t_REAL: return mpsh(x);
    case t_COMPLEX: case t_PADIC:
      av = avma; p1 = gexp(x,prec); p1 = gsub(p1, ginv(p1));
      return gerepileupto(av, gmul2n(p1,-1));
    case t_INTMOD: 
    default:
      av = avma; if (!(y = toser_i(x))) break;
      if (gcmp0(y) && valp(y) == 0) return gcopy(y);
      p1 = gexp(y, prec); p1 = gsub(p1, ginv(p1));
      return gerepileupto(av, gmul2n(p1,-1));
  }
  return transc(gsh,x,prec);
}
/********************************************************************/
/**                                                                **/
/**                      HYPERBOLIC TANGENT                        **/
/**                                                                **/
/********************************************************************/

static GEN
mpth(GEN x)
{
  long lx, s = signe(x);
  GEN y;

  if (!s) return real_0_bit(expo(x));
  lx = lg(x);
  if (absr_cmp(x, stor(bit_accuracy(lx), 3)) >= 0) {
    y = real_1(lx);
  } else {
    pari_sp av = avma;
    long ex = expo(x);
    GEN t;
    if (ex < 1 - BITS_IN_LONG) x = rtor(x, lx + nbits2nlong(-ex)-1);
    t = exp1r_abs(gmul2n(x,1)); /* exp(|2x|) - 1 */
    y = gerepileuptoleaf(av, divrr(t, addsr(2,t)));
  }
  if (s < 0) togglesign(y); /* tanh is odd */
  return y;
}

GEN
gth(GEN x, long prec)
{
  pari_sp av;
  GEN y, t;

  switch(typ(x))
  {
    case t_REAL: return mpth(x);
    case t_COMPLEX: case t_PADIC:
      av = avma;
      t = gexp(gmul2n(x,1),prec);
      t = gdivsg(-2, gaddgs(t,1));
      return gerepileupto(av, gaddsg(1,t));
    case t_INTMOD: pari_err(typeer,"gth");
    default:
      av = avma; if (!(y = toser_i(x))) break;
      if (gcmp0(y)) return gcopy(y);
      t = gexp(gmul2n(y, 1),prec);
      t = gdivsg(-2, gaddgs(t,1));
      return gerepileupto(av, gaddsg(1,t));
  }
  return transc(gth,x,prec);
}
/********************************************************************/
/**                                                                **/
/**                     ARG-HYPERBOLIC SINE                        **/
/**                                                                **/
/********************************************************************/

/* x != 0 */
static GEN
mpash(GEN x)
{
  GEN z, res;
  pari_sp av;
  long lx = lg(x), ex = expo(x);
  
  res = cgetr(lx); av = avma;
  if (ex < 1 - BITS_IN_LONG) x = rtor(x, lx + nbits2nlong(-ex)-1);
  z = logr_abs( addrr_sign(x,1, sqrtr( addrs(mulrr(x,x), 1) ), 1) );
  if (signe(x) < 0) togglesign(z);
  affrr(z, res); avma = av; return res;
}

GEN
gash(GEN x, long prec)
{
  long sx, sy, sz;
  pari_sp av;
  GEN a, y, p1;

  if (gcmp0(x)) return gcopy(x);
  switch(typ(x))
  {
    case t_REAL:
      return mpash(x);

    case t_COMPLEX: av = avma; 
      p1 = gadd(x, gsqrt(gaddsg(1,gsqr(x)), prec));
      y = glog(p1,prec);
      sz = (typ(y)==t_COMPLEX)? gsigne(gel(y,1)): gsigne(y);
      if (typ(p1) == t_COMPLEX) {
        sx = gsigne(gel(p1,1));
        sy = gsigne(gel(p1,2));
      } else {
        sx = gsigne(p1);
        sy = 0;
      }
      if (sx > 0 || (!sx && sy*sz<=0)) return gerepileupto(av, y);

      p1 = mppi(prec); if (sy<0) setsigne(p1,-1);
      return gerepileupto(av, gadd(gneg_i(y), pureimag(p1)));
    case t_INTMOD: case t_PADIC: pari_err(typeer,"gash");
    default:
      av = avma; if (!(y = toser_i(x))) break;
      if (gcmp0(y)) return gcopy(y);
      if (valp(y) < 0) pari_err(negexper,"gash");
      p1 = gaddsg(1,gsqr(y));
      if (gcmp0(p1))
      {
        GEN t = PiI2n(-1,prec);
        if ( gsigne(imag_i(gel(y,2))) < 0 ) setsigne(gel(t,2), -1);
        return gerepileupto(av, scalarser(t, varn(y), valp(p1)>>1));
      }
      p1 = gdiv(derivser(y), gsqrt(p1,prec));
      a = integ(p1,varn(y));
      if (!valp(y)) a = gadd(a, gash(gel(y,2),prec));
      return gerepileupto(av, a);
  }
  return transc(gash,x,prec);
}
/********************************************************************/
/**                                                                **/
/**                     ARG-HYPERBOLIC COSINE                      **/
/**                                                                **/
/********************************************************************/

/* |x| >= 1, return ach(|x|) */
static GEN
mpach(GEN x)
{
  pari_sp av = avma;
  GEN z = logr_abs( addrr_sign(x, 1, sqrtr( subrs(mulrr(x,x), 1) ), 1) );
  return gerepileuptoleaf(av, z);
}

GEN
gach(GEN x, long prec)
{
  pari_sp av;
  GEN a, y, p1;
  long v;

  switch(typ(x))
  {
    case t_REAL:
      if (signe(x) == 0) { y=cgetimag(); gel(y,2) = acos0(expo(x)); return y; }
      if (signe(x) > 0 && expo(x) >= 0) return mpach(x); /* x >= 1 */
      /* -1 < x < 1 */
      if (expo(x) < 0) { y = cgetimag(); gel(y,2) = mpacos(x); return y; }
      /* x <= -1 */
      if (absrnz_egal1(x)) { y = cgetimag(); gel(y,2) = mppi(lg(x)); return y; }
      y = cgetg(3,t_COMPLEX);
      av = avma; p1 = mpach(x);
      setsigne(p1, -signe(p1));
      gel(y,1) = p1;
      gel(y,2) = mppi(lg(x)); return y;

    case t_COMPLEX:
      av = avma; 
      p1 = gadd(x, gsqrt(gaddsg(-1,gsqr(x)), prec)); /* x + sqrt(x^2-1) */
      y = glog(p1,prec);
      if (typ(y) == t_COMPLEX && signe(y[2]) < 0) y = gneg(y);
      return gerepileupto(av, y);

    case t_INTMOD: case t_PADIC: pari_err(typeer,"gach");

    default:
      av = avma; if (!(y = toser_i(x))) break;
      v = valp(y);
      if (v < 0) pari_err(negexper,"gach");
      if (gcmp0(y))
      {
        if (!v) return gcopy(y);
        return gerepileupto(av, gadd(y, PiI2n(-1, prec)));
      }
      p1 = gsubgs(gsqr(y),1);
      if (gcmp0(p1)) { avma = av; return zeroser(varn(y), valp(p1)>>1); }
      p1 = gdiv(derivser(y), gsqrt(p1,prec));
      a = integ(p1, varn(y));
      if (v)
        p1 = PiI2n(-1, prec); /* I Pi/2 */
      else
      {
        p1 = gel(y,2); if (gcmp1(p1)) return gerepileupto(av,a);
        p1 = gach(p1, prec);
      }
      return gerepileupto(av, gadd(p1,a));
  }
  return transc(gach,x,prec);
}
/********************************************************************/
/**                                                                **/
/**                     ARG-HYPERBOLIC TANGENT                     **/
/**                                                                **/
/********************************************************************/

/* |x| < 1, x != 0 */
static GEN
mpath(GEN x)
{
  pari_sp av = avma;
  long ex = expo(x);
  GEN z;
  if (ex < 1 - BITS_IN_LONG) x = rtor(x, lg(x) + nbits2nlong(-ex)-1);
  z = logr_abs( addrs(divsr(2,subsr(1,x)), -1) );
  setexpo(z, expo(z)-1); return gerepileuptoleaf(av, z);
}

GEN
gath(GEN x, long prec)
{
  pari_sp av;
  GEN a, y, p1;

  switch(typ(x))
  {
    case t_REAL:
      if (!signe(x)) return real_0_bit(expo(x));
      if (expo(x) < 0) return mpath(x);

      y = cgetg(3,t_COMPLEX);
      av = avma;
      p1 = addrs(divsr(2,addsr(-1,x)),1);
      if (!signe(p1)) pari_err(talker,"singular argument in atanh");
      p1 = logr_abs(p1);
      setexpo(p1, expo(p1)-1);
      gel(y,1) = gerepileuptoleaf(av, p1);
      gel(y,2) = Pi2n(-1, lg(x)); return y;

    case t_COMPLEX:
      av = avma; p1 = glog( gaddgs(gdivsg(2,gsubsg(1,x)),-1), prec );
      return gerepileupto(av, gmul2n(p1,-1));

    case t_INTMOD: case t_PADIC: pari_err(typeer,"gath");
    default:
      av = avma; if (!(y = toser_i(x))) break;
      if (valp(y) < 0) pari_err(negexper,"gath");
      p1 = gdiv(derivser(y), gsubsg(1,gsqr(y)));
      a = integ(p1, varn(y));
      if (!valp(y)) a = gadd(a, gath(gel(y,2),prec));
      return gerepileupto(av, a);
  }
  return transc(gath,x,prec);
}
/********************************************************************/
/**                                                                **/
/**               CACHE BERNOULLI NUMBERS B_2k                     **/
/**                                                                **/
/********************************************************************/
/* is B_{2k} precomputed at precision >= prec ? */
int
OK_bern(long k, long prec)
{
  return (bernzone && bernzone[1] >= k && bernzone[2] >= prec);
}

#define BERN(i)       (B + 3 + (i)*B[2])
#define set_bern(c0, i, B) STMT_START { \
  *(BERN(i)) = c0; affrr(B, BERN(i)); } STMT_END
/* compute B_0,B_2,...,B_2*nb */
void
mpbern(long nb, long prec)
{
  long i, l, c0;
  pari_sp av;
  GEN B;
  pari_timer T;

  prec++; /* compute one more word of accuracy than required */
  if (OK_bern(nb, prec)) return;
  if (nb < 0) nb = 0;
  l = 3 + prec*(nb+1);
  B = newbloc(l);
  B[0] = evaltyp(t_STR) | evallg(l); /* dummy non-recursive type */
  B[1] = nb;
  B[2] = prec;
  av = avma;

  c0 = evaltyp(t_REAL) | evallg(prec);
  *(BERN(0)) = c0; affsr(1, BERN(0));
  if (bernzone && bernzone[2] >= prec)
  { /* don't recompute known Bernoulli */
    for (i = 1; i <= bernzone[1]; i++) set_bern(c0, i, bern(i));
  }
  else i = 1;
  if (DEBUGLEVEL) {
    fprintferr("caching Bernoulli numbers 2*%ld to 2*%ld, prec = %ld\n",
               i,nb,prec);
    TIMERstart(&T);
  }

  if (i == 1 && nb > 0)
  {
    set_bern(c0, 1, divrs(real_1(prec), 6)); /* B2 = 1/6 */
    i = 2;
  }
  for (   ; i <= nb; i++, avma = av)
  { /* i > 1 */
    long n = 8, m = 5, d1 = i-1, d2 = 2*i-3;
    GEN S = BERN(d1);

    for (;;)
    {
      S = divrs(mulrs(S, n*m), d1*d2);
      if (d1 == 1) break;
      n += 4; m += 2; d1--; d2 -= 2;
      S = addrr(BERN(d1), S);
      if ((d1 & 127) == 0) { set_bern(c0, i, S); S = BERN(i); avma = av; }
    }
    S = divrs(subsr(2*i, S), 2*i+1);
    setexpo(S, expo(S) - 2*i);
    set_bern(c0, i, S); /* S = B_2i */
  }
  if (DEBUGLEVEL) msgTIMER(&T, "Bernoulli");
  if (bernzone) gunclone(bernzone);
  avma = av; bernzone = B;
}
#undef BERN

GEN
bernreal(long n, long prec)
{
  GEN B;

  if (n==1) { B = stor(-1, prec); setexpo(B,-1); return B; }
  if (n<0 || n&1) return gen_0;
  n >>= 1; mpbern(n+1,prec); B=cgetr(prec);
  affrr(bern(n),B); return B;
}

#if 0
/* k > 0 */
static GEN
bernfracspec(long k)
{
  ulong n, K = k+1;
  pari_sp av, lim;
  GEN s, c, N, b;

  c = N = utoipos(K); s = gen_1; b = gen_0;
  av = avma; lim = stack_lim(av,2);
  for (n=2; ; n++) /* n <= k+1 */
  {
    c = diviiexact(muliu(c,k+2-n), utoipos(n));
    if (n & 1) setsigne(c, 1); else setsigne(c, -1);
    /* c = (-1)^(n-1) binomial(k+1, n),  s = 1^k + ... + (n-1)^k */

    b = gadd(b, gdivgs(mulii(c,s), n));
    if (n == K) return gerepileupto(av, b);

    gel(N,2) = n; s = addii(s, powiu(N,k));
    if (low_stack(lim, stack_lim(av,2)))
    {
      if (DEBUGMEM>1) pari_warn(warnmem,"bernfrac");
      gerepileall(av,3, &c,&b,&s);
    }
  }
}
#endif

static GEN
B2(void){ GEN z = cgetg(3, t_FRAC);
  gel(z,1) = gen_1;
  gel(z,2) = utoipos(6); return z;
}
static GEN
B4(void) { GEN z = cgetg(3, t_FRAC);
  gel(z,1) = gen_m1;
  gel(z,2) = utoipos(30); return z;
}

GEN
bernfrac(long k)
{
  if (k < 6) switch(k)
  {
    case 0: return gen_1;
    case 1: return gneg(ghalf);
    case 2: return B2();
    case 4: return B4();
    default: return gen_0;
  }
  if (k & 1) return gen_0;
  return bernfrac_using_zeta(k);
}

/* mpbern as exact fractions */
static GEN
bernvec_old(long nb)
{
  long n, i;
  GEN y;

  if (nb < 0) return cgetg(1, t_VEC);
  if (nb > 46340 && BITS_IN_LONG == 32) pari_err(impl, "bernvec for n > 46340");

  y = cgetg(nb+2, t_VEC); gel(y,1) = gen_1;
  for (n = 1; n <= nb; n++)
  { /* compute y[n+1] = B_{2n} */
    pari_sp av = avma;
    GEN b = gmul2n(utoineg(2*n - 1), -1); /* 1 + (2n+1)B_1 = -(2n-1) /2 */
    GEN c = gen_1;
    ulong u1 = 2*n + 1, u2 = n, d1 = 1, d2 = 1;

    for (i = 1; i < n; i++)
    {
      c = diviiexact(muliu(c, u1*u2), utoipos(d1*d2));/*= binomial(2n+1, 2*i) */
      b = gadd(b, gmul(c, gel(y,i+1)));
      u1 -= 2; u2--; d1++; d2 += 2;
    }
    gel(y,n+1) = gerepileupto(av, gdivgs(b, -(1+2*n)));
  }
  return y;
}
GEN
bernvec(long nb)
{
  GEN y = cgetg(nb+2, t_VEC), z = y + 1;
  long i;
  if (nb < 20) return bernvec_old(nb);
  for (i = nb; i > 2; i--) gel(z,i) = bernfrac_using_zeta(i << 1);
  gel(y,3) = B4();
  gel(y,2) = B2();
  gel(y,1) = gen_1; return y;
}

/********************************************************************/
/**                                                                **/
/**                         EULER'S GAMMA                          **/
/**                                                                **/
/********************************************************************/

/* x / (i*(i+1)) */
GEN
divrsns(GEN x, long i)
{
#ifdef LONG_IS_64BIT
  if (i < 3037000500) /* i(i+1) < 2^63 */
#else
  if (i < 46341) /* i(i+1) < 2^31 */
#endif
    return divrs(x, i*(i+1));
  else
    return divrs(divrs(x, i), i+1);
}
/* x / (i*(i+1)) */
GEN
divgsns(GEN x, long i)
{
#ifdef LONG_IS_64BIT
  if (i < 3037000500) /* i(i+1) < 2^63 */
#else
  if (i < 46341) /* i(i+1) < 2^31 */
#endif
    return gdivgs(x, i*(i+1));
  else
    return gdivgs(gdivgs(x, i), i+1);
}

/* arg(s+it) */
double
darg(double s, double t)
{
  double x;
  if (!t) return (s>0)? 0.: PI;
  if (!s) return (t>0)? PI/2: -PI/2;
  x = atan(t/s);
  return (s>0)? x
              : ((t>0)? x+PI : x-PI);
}

void
dcxlog(double s, double t, double *a, double *b)
{
  *a = log(s*s + t*t) / 2; /* log |s| = Re(log(s)) */
  *b = darg(s,t);          /* Im(log(s)) */
}

double
dabs(double s, double t) { return sqrt( s*s + t*t ); }
double
dnorm(double s, double t) { return s*s + t*t; }

GEN
trans_fix_arg(long *prec, GEN *s0, GEN *sig, pari_sp *av, GEN *res)
{
  GEN s, p1;
  long l;
  if (typ(*s0)==t_COMPLEX && gcmp0(gel(*s0,2))) *s0 = gel(*s0,1);
  s = *s0;
  l = precision(s); if (!l) l = *prec;
  if (l < 3) l = 3;

  if (typ(s) == t_COMPLEX)
  { /* s = sig + i t */
    *res = cgetc(l); *av = avma;
    s = ctofp(s, l+1); *sig = gel(s,1);
  }
  else /* real number */
  {
    *res = cgetr(l); *av = avma;
    *sig = s = gtofp(s, l+1);
    p1 = floorr(s);
    if (!signe(subri(s,p1))) *s0 = p1;
  }
  *prec = l; return s;
}

#if 0
/* x, z t_REAL. Compute unique x in ]-z,z] congruent to x mod 2z */
static GEN
red_mod_2z(GEN x, GEN z)
{
  GEN Z = gmul2n(z, 1), d = subrr(z, x);
  /* require little accuracy */
  if (!signe(d)) return x;
  setlg(d, 3 + ((expo(d) - expo(Z)) >> TWOPOTBITS_IN_LONG));
  return addrr(mulir(floorr(divrr(d, Z)), Z), x);
}
#endif

/* update lg(z) before affrr(y, z)  [ to cater for precision loss ]*/
void
affr_fixlg(GEN y, GEN z) {
  long ly = lg(y), lz = lg(z);
  if (ly < lz)
  {
    setlg(z, ly);
    stackdummy((pari_sp)(z + lz), (pari_sp)(z + ly));
  }
  /* lz <= ly */
  affrr(y, z);
}

static GEN
cxgamma(GEN s0, int dolog, long prec)
{
  GEN s, u, a, y, res, tes, sig, invn2, p1, nnx, pi, pi2, sqrtpi2;
  long i, lim, nn;
  pari_sp av, av2, avlim;
  int funeq = 0;

  if (DEBUGLEVEL>5) (void)timer2();
  s = trans_fix_arg(&prec,&s0,&sig,&av,&res);

  if ((signe(sig) <= 0 || expo(sig) < -1)
    && (typ(s) == t_REAL || gexpo(gel(s,2)) <= 16))
  { /* s <--> 1-s */
    funeq = 1; s = gsub(gen_1, s); sig = real_i(s);
  }

  { /* find "optimal" parameters [lim, nn] */
    double ssig = rtodbl(sig);
    double st = rtodbl(imag_i(s));
    double la, l,l2,u,v, rlogs, ilogs;

    dcxlog(ssig,st, &rlogs,&ilogs);
    /* Re (s - 1/2) log(s) */
    u = (ssig - 0.5)*rlogs - st * ilogs;
    /* Im (s - 1/2) log(s) */
    v = (ssig - 0.5)*ilogs + st * rlogs;
    /* l2 = | (s - 1/2) log(s) - s + log(2Pi)/2 |^2 ~ |lngamma(s))|^2 */
    u = u - ssig + log(2.*PI)/2;
    v = v - st;
    l2 = u*u + v*v;
    if (l2 < 0.000001) l2 = 0.000001;
    l = (bit_accuracy_mul(prec, LOG2) - log(l2)/2) / 2.;
    if (l < 0) l = 0.;

    la = 3.; /* FIXME: heuristic... */
    if (st > 1 && l > 0)
    {
      double t = st * PI / l;
      la = t * log(t);
      if (la < 3) la = 3.;
      if (la > 150) la = t;
    }
    lim = (long)ceil(l / (1.+ log(la)));
    if (lim == 0) lim = 1;

    u = (lim-0.5) * la / PI;
    l2 = u*u - st*st;
    if (l2 > 0)
    {
      nn = (long)ceil(sqrt(l2) - ssig);
      if (nn < 1) nn = 1;
    }
    else
      nn = 1;
#if 0
#define pariK2 (1.1239968) /* 1/(1-(log(2)/(2*pi))) */
#define pariK4 (17.079468445347/BITS_IN_LONG) /* 2*e*pi/BIL */
    {/* same: old method */
      long e = gexpo(s);
      double beta;
      if (e > 1000)
      {
        nn = 0;
        beta = log(pariK4 / (prec-2)) / LOG2 + e;
        if (beta > 1.) beta += log(beta)/LOG2;
        lim = (long)((bit_accuracy(prec)>>1)/beta + 1);
      }
      else
      {
        double alpha = sqrt( dnorm(ssig, st) );
        beta = bit_accuracy_mul(prec,LOG2/(2*PI)) - alpha;
        if (beta >= 0) nn = (long)(1+pariK2*beta); else nn = 0;
        if (nn)
          lim = (long)(1+PI*(alpha+nn));
        else
        {
          beta = log( pariK4 * alpha / (prec-2) ) / LOG2;
          if (beta > 1.) beta += log(beta)/LOG2;
          lim = (long)((bit_accuracy(prec)>>1)/beta + 1);
        }
      }
      nn++;
    }
#endif
    if (DEBUGLEVEL>5) fprintferr("lim, nn: [%ld, %ld], la = %lf\n",lim,nn,la);
  }
  prec++;

  av2 = avma; avlim = stack_lim(av2,3);
  y = s;
  if (typ(s0) == t_INT)
  {
    if (signe(s0) <= 0) pari_err(talker,"non-positive integer argument in cxgamma");
    if (is_bigint(s0)) {
      for (i=1; i < nn; i++)
      {
        y = mulri(y, addis(s0, i));
        if (low_stack(avlim,stack_lim(av2,3)))
        {
          if(DEBUGMEM>1) pari_warn(warnmem,"gamma");
          y = gerepileuptoleaf(av2, y);
        }
      }
    } else {
      ulong ss = itou(s0);
      for (i=1; i < nn; i++)
      {
        y = mulru(y, ss + i);
        if (low_stack(avlim,stack_lim(av2,3)))
        {
          if(DEBUGMEM>1) pari_warn(warnmem,"gamma");
          y = gerepileuptoleaf(av2, y);
        }
      }
    }
    if (dolog) y = logr_abs(y);
  }
  else if (!dolog || typ(s) == t_REAL)
  { /* Compute lngamma mod 2 I Pi */
    for (i=1; i < nn; i++)
    {
      y = gmul(y, gaddgs(s,i));
      if (low_stack(avlim,stack_lim(av2,3)))
      {
        if(DEBUGMEM>1) pari_warn(warnmem,"gamma");
        y = gerepileupto(av2, y);
      }
    }
    if (dolog) y = logr_abs(y);
  }
  else
  { /* dolog && complex s: be careful with imaginary part */
    y = glog(y, prec);
    for (i=1; i < nn; i++)
    {
      y = gadd(y, glog(gaddgs(s,i), prec));
      if (low_stack(avlim,stack_lim(av2,3)))
      {
        if(DEBUGMEM>1) pari_warn(warnmem,"gamma");
        y = gerepileupto(av2, y);
      }
    }
  }
  if (DEBUGLEVEL>5) msgtimer("product from 0 to N-1");

  nnx = gaddgs(s, nn);
  a = ginv(nnx); invn2 = gsqr(a);
  tes = divrsns(bernreal(2*lim,prec), 2*lim-1); /* B2l / (2l-1) 2l*/
  if (DEBUGLEVEL>5) msgtimer("Bernoullis");
  for (i = 2*lim-2; i > 1; i -= 2)
  {
    u = divrsns(bernreal(i,prec), i-1); /* Bi / i(i-1) */
    tes = gadd(u, gmul(invn2,tes));
  }
  if (DEBUGLEVEL>5) msgtimer("Bernoulli sum");

  p1 = gsub(gmul(gsub(nnx, ghalf), glog(nnx,prec)), nnx);
  p1 = gadd(p1, gmul(tes, a));

  pi = mppi(prec); pi2 = shiftr(pi, 1); sqrtpi2 = sqrtr(pi2);

  if (dolog)
  {
    if (funeq)
    { /* 2 Pi ceil( (2Re(s) - 3)/4 ) */
      GEN z = mulri(pi2, ceilr(shiftr(subrs(shiftr(sig,1), 3), -2)));
      /* y --> y + log Pi - log sqrt(2Pi) - log sin(Pi s)
       *     = y - log( sin(Pi s) / (sqrt(2Pi)/2) ) */
      y = gsub(y, glog(gdiv(gsin(gmul(pi,s0),prec), shiftr(sqrtpi2,-1)), prec));
      if (signe(z)) {
        if (gsigne(imag_i(s)) < 0) setsigne(z, -signe(z));
        if (typ(y) == t_COMPLEX)
          gel(y,2) = gadd(gel(y,2), z);
        else
          y = gadd(y, pureimag(z));
      }
      p1 = gneg(p1);
    }
    else /* y --> sqrt(2Pi) / y */
      y = gsub(logr_abs(sqrtpi2), y);
    y = gadd(p1, y);
  }
  else
  {
    if (funeq)
    { /* y --> y Pi/(sin(Pi s) * sqrt(2Pi)) = y sqrt(Pi/2)/sin(Pi s) */
      y = gdiv(gmul(shiftr(sqrtpi2,-1),y), gsin(gmul(pi,s0), prec));
      /* don't use s above: sin(pi s0) = sin(pi s) and the former is
       * more accurate, esp. if s0 ~ 0 */
      p1 = gneg(p1);
    }
    else /* y --> sqrt(2Pi) / y */
      y = gdiv(sqrtpi2, y);
    y = gmul(gexp(p1, prec), y);
  }
  if (typ(y) == t_REAL) affr_fixlg(y, res);
  else
  {
    if (typ(res) == t_REAL) return gerepileupto(av, y);
    affr_fixlg(gel(y,1), gel(res,1));
    affr_fixlg(gel(y,2), gel(res,2));
  }
  avma = av; return res;
}

/* Gamma((m+1) / 2) */
static GEN
gammahs(long m, long prec)
{
  GEN y = cgetr(prec), z;
  pari_sp av = avma;
  long ma = labs(m);

  if (ma > 200 + 50*(prec-2)) /* heuristic */
  {
    z = stor(m + 1, prec); setexpo(z, expo(z)-1);
    affrr(cxgamma(z,0,prec), y);
    avma = av; return y;
  }
  z = sqrtr( mppi(prec) );
  if (m)
  {
    GEN p1 = seq_umul(ma/2 + 1, ma);
    long v = vali(p1);
    p1 = shifti(p1, -v); v -= ma;
    if (m >= 0) z = mulri(z,p1);
    else
    {
      z = divri(z,p1); v = -v;
      if ((m&3) == 2) setsigne(z,-1);
    }
    setexpo(z, expo(z) + v);
  }
  affrr(z, y); avma = av; return y;
}
GEN
ggamd(GEN x, long prec)
{
  pari_sp av, tetpil;

  switch(typ(x))
  {
    case t_INT:
    {
      long k = itos(x);
      if (labs(k) > 962353) pari_err(talker, "argument too large in ggamd");
      return gammahs(k<<1, prec);
    }
    case t_REAL: case t_FRAC: case t_COMPLEX: case t_QUAD: case t_PADIC:
      av=avma; x = gadd(x,ghalf); tetpil=avma;
      return gerepile(av,tetpil,ggamma(x,prec));

    case t_INTMOD: pari_err(typeer,"ggamd");
    case t_SER: pari_err(impl,"gamd of a power series");
  }
  return transc(ggamd,x,prec);
}

/* find n such that n-v_p(n!)>=k */
static long nboft(long k, long p)
{
  long s,n;
  for (s=0,n=0; n-s<k; s += u_lval(++n, p));
  return n; 
}

/*
 * Using Dwork's expansion, compute \Gamma(px+1)=-\Gamma(px) with x a
 * unit.
 * See p$-Adic Gamma Functions and Dwork Cohomology,
 * Maurizio Boyarsky
 * Transactions of the American Mathematical Society,
 * Vol. 257, No. 2. (Feb., 1980), pp. 359-369.
 * Inspired by a GP script by Fernando Rodriguez-Villegas 
 */

static GEN
gadw(GEN x, long p)
{
  pari_sp ltop=avma;
  GEN s, t;
  long j, k;
  long n = nboft(precp(x)+valp(x)+1,p);
  GEN  u = cgetg(p+1, t_VEC);
  s = gaddsg(1, zeropadic(gel(x,2), n));
  t = s;
  gel(u, 1) = s;
  for (j = 1; j < p; ++j)
    gel(u, j + 1) = gdivgs(gel(u, j), j);
  for (k = 1; k < n; ++k)
  {
    gel(u, 1) = gdivgs(gdivgs(gadd(gel(u, 1), gel(u, p)), k), p);
    for (j = 1; j < p; ++j)
      gel(u, j + 1) = gdivgs(gadd(gel(u, j), gel(u, j + 1)), (k*p) + j);
    
    t = gmul(t, gaddgs(x, k-1));
    s = gadd(s, gmul(gmul(gel(u, 1), gpowgs(gel(x,2), k)), t));
    if ((k&0xFL)==0) gerepileall(ltop, 3, &u,&s,&t);
  }
  
  return gneg(s);
}

/*Use Dwork expansion*/
/*This is a O(p*e*log(pe)) algorithm, should be used when p small
 * If p==2 this is a O(pe) algorithm...
 */
static GEN
gammap_Dwork(GEN x, long p)
{
  pari_sp ltop = avma;
  long k = itos(gmodgs(x, p));
  GEN p1;
  long j;
  if (k)
  {
    x = gdivgs(gsubgs(x, k), p);
    k--;
    p1 = gadw(x, p);
    if (k%2==1) p1 = gneg(p1);
    for (j = 1; j <= k; ++j)
      p1 = gmul(p1, gaddgs(gmulsg(p, x), j));
  }
  else
    p1 = gneg(gadw(gdivgs(x, p), p));
  return gerepileupto(ltop, p1);
}

/* 
 * Compute gammap using the definition. This is a O(x*M(log(pe))) algorithm.
 * This should be used if x is very small.
 */
static GEN
gammap_Morita(long n, GEN p, long e)
{
  pari_sp ltop=avma;
  GEN p2 = gaddsg((n&1)?-1:1, zeropadic(p, e+1));
  long i;
  long pp=is_bigint(p)? 0: itos(p);
  for (i = 2; i < n; i++)
    if (!pp || i%pp)
    {
      p2 = gmulgs(p2, i);
      if ((i&0xFL) == 0xFL)
        p2 = gerepileupto(ltop, p2);
    }
  return gerepileupto(ltop, p2);
}

/*
 * x\in\N: Gamma(-x)=(-1)^(1+x+x\p)*Gamma(1+x)
 */

static GEN
gammap_neg_Morita(long n, GEN p, long e)
{
  GEN g = ginv(gammap_Morita(n+1,p,e));
  return ((n^sdivsi(n,p)) & 1)? g :gneg(g);
}

/* p-adic Gamma function for x a p-adic integer */
/*
  There are three cases:
  n is small            : we use Morita definition.
  n is large, p is small: we use Dwork expansion.
  n is large, p is large: we don't know how to proceed.
  TODO: handle p=2 better (gammap_Dwork is very slow for p=2).
*/
#define GAMMAP_DWORK_LIMIT 50000UL
static GEN
gammap(GEN x)
{
  GEN p = gel(x,2);
  long e= precp(x);
  GEN n,m,nm;
  if (valp(x)<0) 
    pari_err(talker,"Gamma not defined for non-integral p-adic number");
  n = gtrunc(x);
  m = gtrunc(gneg(x));
  nm= cmpii(n,m)<=0?n:m;
  if (lgefint(nm)==3 && (is_bigint(p) || (ulong)nm[2]<GAMMAP_DWORK_LIMIT))
  {
    if(nm==n)
      return gammap_Morita(itos(n),p,e);
    else 
      return gammap_neg_Morita(itos(m),p,e);
  }
  return gammap_Dwork(x, itos(p));
}

GEN
ggamma(GEN x, long prec)
{
  pari_sp av;
  long m;
  GEN y, z;

  switch(typ(x))
  {
    case t_INT:
      if (signe(x) <= 0) pari_err(talker,"non-positive integer argument in ggamma");
      if (cmpiu(x,481177) > 0) pari_err(talker,"argument too large in ggamma");
      return mpfactr(itos(x) - 1, prec);

    case t_REAL: case t_COMPLEX:
      return cxgamma(x, 0, prec);

    case t_FRAC:
      if (!equaliu(gel(x,2),2)) break;
      z = gel(x,1); /* true argument is z/2 */
      if (is_bigint(z) || labs(m = itos(z)) > 962354)
      {
        pari_err(talker, "argument too large in ggamma");
        return NULL; /* not reached */
      }
      return gammahs(m-1, prec);

    case t_PADIC: return gammap(x);
    case t_INTMOD: pari_err(typeer,"ggamma");
    default:
      av = avma; if (!(y = toser_i(x))) break;
      return gerepileupto(av, gexp(glngamma(y,prec),prec));
  }
  return transc(ggamma,x,prec);
}

GEN
mpfactr(long n, long prec)
{
  GEN f = cgetr(prec);
  pari_sp av = avma;

  if (n+1 > 350 + 70*(prec-2)) /* heuristic */
    affrr(cxgamma(stor(n+1, prec), 0, prec), f);
  else
    affir(mpfact(n), f);
  avma = av; return f;
}

GEN
glngamma(GEN x, long prec)
{
  long i, n;
  pari_sp av;
  GEN a, y, p1;

  switch(typ(x))
  {
    case t_INT:
      if (signe(x) <= 0) pari_err(talker,"non-positive integer in glngamma");
      if (cmpiu(x,200 + 50*(prec-2)) > 0) /* heuristic */
	return cxgamma(x, 1, prec);
      av = avma;
      return gerepileuptoleaf(av, logr_abs( itor(mpfact(itos(x) - 1), prec) ));

    case t_REAL: case t_COMPLEX:
      return cxgamma(x, 1, prec);

    default:
      av = avma; if (!(y = toser_i(x))) break;
      if (valp(y)) pari_err(negexper,"glngamma");
      p1 = gsubsg(1,y);
      if (!valp(p1)) pari_err(impl,"lngamma around a!=1");
      n = (lg(y)-3) / valp(p1);
      a = zeroser(varn(y), lg(y)-2);
      for (i=n; i>=2; i--) a = gmul(p1, gadd(a, gdivgs(szeta(i, prec),i)));
      a = gadd(a, mpeuler(prec));
      return gerepileupto(av, gmul(a, p1));

    case t_PADIC:  pari_err(impl,"p-adic lngamma function");
    case t_INTMOD: pari_err(typeer,"glngamma");
  }
  return transc(glngamma,x,prec);
}
/********************************************************************/
/**                                                                **/
/**                  PSI(x) = GAMMA'(x)/GAMMA(x)                   **/
/**                                                                **/
/********************************************************************/

GEN
cxpsi(GEN s0, long prec)
{
  pari_sp av, av2;
  GEN sum,z,a,res,tes,in2,sig,s,unr;
  long lim,nn,k;
  const long la = 3;
  int funeq = 0;

  if (DEBUGLEVEL>2) (void)timer2();
  s = trans_fix_arg(&prec,&s0,&sig,&av,&res);
  if (signe(sig) <= 0) { funeq = 1; s = gsub(gen_1, s); sig = real_i(s); }
  if (typ(s0) == t_INT && signe(s0) <= 0)
    pari_err(talker,"non-positive integer argument in cxpsi");

  {
    double ssig = rtodbl(sig);
    double st = rtodbl(imag_i(s));
    double l;
    {
      double rlog, ilog; /* log (s - Euler) */
      dcxlog(ssig - 0.57721566, st, &rlog,&ilog);
      l = dnorm(rlog,ilog);
    }
    if (l < 0.000001) l = 0.000001;
    l = log(l) / 2.;
    lim = 2 + (long)ceil((bit_accuracy_mul(prec, LOG2) - l) / (2*(1+log((double)la))));
    if (lim < 2) lim = 2;

    l = (2*lim-1)*la / (2.*PI);
    l = l*l - st*st;
    if (l < 0.) l = 0.;
    nn = (long)ceil( sqrt(l) - ssig );
    if (nn < 1) nn = 1;
    if (DEBUGLEVEL>2) fprintferr("lim, nn: [%ld, %ld]\n",lim,nn);
  }
  prec++; unr = real_1(prec); /* one extra word of precision */

  a = gdiv(unr, gaddgs(s, nn)); /* 1 / (s+n) */
  av2 = avma; sum = gmul2n(a,-1);
  for (k = 0; k < nn; k++)
  {
    sum = gadd(sum, gdiv(unr, gaddgs(s, k)));
    if ((k & 127) == 0) sum = gerepileupto(av2, sum);
  }
  z = gsub(glog(gaddgs(s, nn), prec), sum);
  if (DEBUGLEVEL>2) msgtimer("sum from 0 to N-1");

  in2 = gsqr(a);
  av2 = avma; tes = divrs(bernreal(2*lim, prec), 2*lim);
  for (k=2*lim-2; k>=2; k-=2)
  {
    tes = gadd(gmul(in2,tes), divrs(bernreal(k, prec), k));
    if ((k & 255) == 0) tes = gerepileupto(av2, tes);
  }
  if (DEBUGLEVEL>2) msgtimer("Bernoulli sum");
  z = gsub(z, gmul(in2,tes));
  if (funeq)
  {
    GEN pi = mppi(prec);
    z = gadd(z, gmul(pi, gcotan(gmul(pi,s), prec)));
  }
  if (typ(z) == t_REAL) affr_fixlg(z, res);
  else
  {
    affr_fixlg(gel(z,1), gel(res,1));
    affr_fixlg(gel(z,2), gel(res,2));
  }
  avma = av; return res;
}

GEN
gpsi(GEN x, long prec)
{
  switch(typ(x))
  {
    case t_REAL: case t_COMPLEX: return cxpsi(x,prec);
    case t_INTMOD: case t_PADIC: pari_err(typeer,"gpsi");
    case t_SER: pari_err(impl,"psi of power series");
  }
  return transc(gpsi,x,prec);
}
