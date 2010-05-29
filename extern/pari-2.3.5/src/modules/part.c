/* $Id: part.c 8641 2007-05-02 21:03:43Z kb $

Copyright (C) 2002  The PARI group.

This file is part of the PARI/GP package.

PARI/GP is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation. It is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY WHATSOEVER.

Check the License for details. You should have received a copy of it, along
with the package; see the file 'COPYING'. If not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA. */

/* Original code contributed by: Ralf Stephan (ralf@ark.in-berlin.de).
 * 
 * This program is basically the implementation of the script
 *
 * Psi(n, q) = local(a, b, c); a=sqrt(2/3)*Pi/q; b=n-1/24; c=sqrt(b);
 *             (sqrt(q)/(2*sqrt(2)*b*Pi))*(a*cosh(a*c)-(sinh(a*c)/c))
 * L(n, q) = if(q==1,1,sum(h=1,q-1,if(gcd(h,q)>1,0,cos((g(h,q)-2*h*n)*Pi/q))))
 * g(h, q) = if(q<3,0,sum(k=1,q-1,k*(frac(h*k/q)-1/2)))
 * part(n) = round(sum(q=1,5 + 0.24*sqrt(n),L(n,q)*Psi(n,q)))
 *
 * only faster. It is a translation of the C/mpfr version at
 * http://www.ark.in-berlin.de/part.c
 *
 * ------------------------------------------------------------------
 *   The first restriction depends on Pari's maximum precision of floating
 * point reals, which is 268435454 bits in 2.2.4, since the algorithm needs
 * high precision exponentials. For that engine, the maximum possible argument
 * would be in [5*10^15,10^16], the computation of which would need days on 
 * a ~1-GHz computer. */

#include "pari.h"
#include "paripriv.h"

/****************************************************************/

/* Given:  b = N-1/24;
 *   c = sqrt(2/3)*Pi*sqrt(b)
 *   d = 1 / ((2*b)^(3/2) * Pi);
 *
 * Psi(N, q) = local(a = c/q); sqrt(q) * (a*cosh(a) - sinh(a)) */
static GEN
psi(GEN c, ulong q, long prec)
{
  GEN a = divrs(c, q), ea = mpexp(a), invea = ginv(ea);
  GEN cha = shiftr(mpadd(ea, invea), -1);  /* ch(a) */
  GEN sha = shiftr(mpsub(ea, invea), -1);  /* sh(a) */
  return mulrr(sqrtr(stor(q,prec)), subrr(mulrr(a,cha), sha));
}

/* g(h, q) = if (q<3, 0, sum(k=1,q-1,k*(frac(h*k/q)-1/2))) */
/* assume h < q and (h,q) = 1. Not memory clean. */
static GEN
g(ulong q, ulong h)
{
  ulong k, kh;
  GEN i2;

  if (q < 3)  return gen_0;
  if (h == 1) return gdivgs(mulss(q-1,q-2), 12);
  if (h == 2) return gdivgs(mulss(q-1,q-5), 24); /* q odd since (h,q)=1 */

  k = q % h;
  if (k == 1)
    return gdivgs(mulsi((q-1)/h, subsi(q-1, mulss(h,h))), 12);
  if (k == 2)
    return gdivgs(mulsi((q-2)/h, subsi(q<<1, addsi(1, mulss(h,h)))), 24);
  /* TODO: expr for h-1 mod h  +  gcd-style computation */
  
  kh = h;
  if (MAXULONG/h > q)
  {
    long l3 = 0;
    for (k = 1; k < q; k++)
    {
      l3 += k * ((kh << 1) - q);
      kh += h; if (kh >= q) kh -= q;
    }
    i2 = stoi(l3);
  }
  else
  {
    pari_sp av = avma;
    i2 = gen_0;
    for (k = 1; k < q; k++)
    {
      i2 = addii(i2, mulss(k, (kh << 1) - q));
      if ((k & 31) == 0) i2 = gerepileuptoint(av, i2);
      kh += h; if (kh >= q) kh -= q;
    }
  }
  return gdivgs(i2, q<<1);
}

/* L(n, q) = if(q==1,1,sum(h=1,q-1,if(gcd(h,q)>1,0,cos((g(h,q)-2*h*n)*Pi/q))) */
static GEN
L(GEN n, ulong q, long bitprec)
{
  long pr = nbits2prec(bitprec / q + q);
  ulong h, nmodq = umodiu(n, q), hn;
  GEN r, res = stor(0, pr), pi_q = divrs(mppi(pr), q);
  pari_sp av = avma;
  for (h = 1, hn = 0; h < q; h++, avma = av)
  {
    GEN t;
    hn += nmodq; if (hn >= q) hn -= q;
    if (cgcd((long)q, (long)h) > 1) continue;
    r = gsubgs(g(q,h), hn << 1);
    t = gcmp0(r)? addrs(res, 1): addrr(res, mpcos(gmul(pi_q,r)));
    affrr(t, res);
  }
  return res;
}

/* Return a low precision raw estimate of log p(n). */
static GEN
estim(GEN n)
{
  pari_sp av = avma;
  GEN p1, pi = mppi (DEFAULTPREC);

  p1 = divrs( itor(shifti(n,1), DEFAULTPREC), 3 );
  p1 = mpexp( mulrr(pi, sqrtr(p1)) ); /* exp(Pi * sqrt(2N/3)) */
  p1 = divri (shiftr(p1,-2), n);
  p1 = divrr(p1, sqrtr( stor(3,DEFAULTPREC) ));
  return gerepileupto(av, mplog(p1));
}

static void
pinit(GEN n, GEN *c, GEN *d, ulong prec)
{
  GEN b = divrs( itor( subis(mulis(n,24), 1), prec ), 24 ); /* n - 1/24 */
  GEN sqrtb = sqrtr(b), Pi = mppi (prec), pi2sqrt2, pisqrt2d3;

  pisqrt2d3 = mulrr(Pi, sqrtr( divrs(stor(2, prec), 3) ));
  pi2sqrt2  = mulrr(Pi, sqrtr( stor(8, prec) ));
  *c = mulrr(pisqrt2d3, sqrtb);
  *d = ginv( mulrr(pi2sqrt2, mulrr(b,sqrtb)) );
}

/* part(n) = round(sum(q=1,5 + 0.24*sqrt(n), L(n,q)*Psi(n,q))) */
GEN
numbpart(GEN n)
{
  pari_sp ltop = avma, av;
  GEN sum, est, C, D, p1, p2;
  long prec, bitprec;
  ulong q;

  if (typ(n) != t_INT) pari_err(typeer, "partition function");
  if (signe(n) < 0) return gen_0;
  if (cmpiu(n, 2) < 0) return gen_1;
  if (cmpii(n, u2toi(0x38d7e, 0xa4c68000)) >= 0)
    pari_err(talker, "arg to partition function must be < 10^15");
  est = estim(n);
  bitprec = (long)(rtodbl(est)/LOG2) + 32;
  prec = nbits2prec(bitprec);
  pinit(n, &C, &D, prec);
  sum = cgetr (prec);
  /* Because N < 10^16 and q < sqrt(N), q fits into a long */
  av = avma; togglesign(est);
  for (q = (ulong)(sqrt(gtodouble(n))*0.24 + 5); q >= 3; q--, avma=av)
  {
    GEN t = L(n, q, bitprec);
    if (absr_cmp(t, mpexp(divrs(est,q))) < 0) continue;

    t = mulrr(t, psi(gprec_w(C, nbits2prec(bitprec / q + 32)), q, prec));
    affrr(addrr(sum, t), sum);
  }
  p1 = addrr(sum, psi(C, 1, prec));
  p2 = psi(C, 2, prec);
  affrr(mod2(n)? subrr(p1,p2): addrr(p1,p2), sum);
  return gerepileuptoint (ltop, roundr(mulrr(D,sum)));
}
