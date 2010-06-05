#line 2 "../src/kernel/gmp/gcd.c"
/* $Id: gcd.c 7522 2005-12-09 18:14:24Z kb $

Copyright (C) 2000-2003  The PARI group.

This file is part of the PARI/GP package.

PARI/GP is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation. It is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY WHATSOEVER.

Check the License for details. You should have received a copy of it, along
with the package; see the file 'COPYING'. If not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA. */

GEN
gcdii(GEN a, GEN b)
{
  long v, w;
  pari_sp av;
  GEN t, p1;

  switch (absi_cmp(a,b))
  {
    case 0: return absi(a);
    case -1: t=b; b=a; a=t;
  }
  if (!signe(b)) return absi(a);
  /* here |a|>|b|>0. Try single precision first */
  if (lgefint(a)==3)
    return gcduu((ulong)a[2], (ulong)b[2]);
  if (lgefint(b)==3)
  {
    ulong u = resiu(a,(ulong)b[2]);
    if (!u) return absi(b);
    return gcduu((ulong)b[2], u);
  }
  /* larger than gcd: "avma=av" gerepile (erasing t) is valid */
  av = avma; (void)new_chunk(lgefint(b)+1); /* HACK */
  t = remii(a,b);
  if (!signe(t)) { avma=av; return absi(b); }

  a = b; b = t;
  v = vali(a); a = shifti(a,-v); setsigne(a,1);
  w = vali(b); b = shifti(b,-w); setsigne(b,1);
  if (w < v) v = w;
  switch(absi_cmp(a,b))
  {
    case  0: avma=av; a=shifti(a,v); return a;
    case -1: p1=b; b=a; a=p1;
  }
  if (is_pm1(b)) { avma=av; return int2n(v); }
 {
  /* general case */
  /*This serve two purposes: 1) mpn_gcd destroy its input and need an extra
   * limb 2) this allows us to use icopy instead of gerepile later.  NOTE: we
   * must put u before d else the final icopy could fail.
   */
  GEN res= cgeti(lgefint(a)+1);
  GEN ca = icopy_ef(a,lgefint(a)+1);
  GEN cb = icopy_ef(b,lgefint(b)+1);
  long l = mpn_gcd(LIMBS(res), LIMBS(ca), NLIMBS(ca), LIMBS(cb), NLIMBS(cb));
  res[1] = evalsigne(1)|evallgefint(l+2);
  avma=av;
  return shifti(res,v);
  }
}

/*==================================
 * invmod(a,b,res)
 *==================================
 *    If a is invertible, return 1, and set res  = a^{ -1 }
 *    Otherwise, return 0, and set res = gcd(a,b)
 */
int
invmod(GEN a, GEN b, GEN *res)
{
  if (typ(a) != t_INT || typ(b) != t_INT) pari_err(arither1);
  if (!signe(b)) { *res=absi(a); return 0; }
  if (NLIMBS(b) < INVMOD_GMP_LIMIT) 
    return invmod_pari(a,b,res);
  { /* General case: use gcdext(a+b, b) since mpn_gcdext require S1>=S2 */
    pari_sp av = avma;
    GEN  ca, cb, u, d;
    long lu, l, su, sa = signe(a), lb,lna;
    GEN na;
    if (!sa) { avma = av; *res = absi(b); return 0; }
    if (signe(b) < 0) b = negi(b);
    if (absi_cmp(a, b) < 0)
      na = sa > 0? addii(a, b): subii(a, b);
    else
      na = a;
    /* Copy serves two purposes:
     * 1) mpn_gcdext destroys its input and needs an extra limb
     * 2) allows us to use icopy instead of gerepile later. */
    lb = lgefint(b); lna = lgefint(na);
    ca = icopy_ef(na,lna+1);
    cb = icopy_ef( b,lb+1);
    /* must create u first else final icopy could fail. */
    u = cgeti(lna+1);
    d = cgeti(lna+1);
    /* na >= b */
    l = mpn_gcdext(LIMBS(d), LIMBS(u), &lu, LIMBS(ca), NLIMBS(ca), LIMBS(cb), NLIMBS(cb));
    d[1] = evalsigne(1)|evallgefint(l+2);
    if (!is_pm1(d)) {avma=av; *res=icopy(d); return 0;}
    su = lu?((sa ^ lu) < 0)? -1: 1: 0;
    u[1] = evalsigne(su) | evallgefint(labs(lu)+2);
    if (su < 0) u = addii(u, b);
    avma=av; *res=icopy(u); return 1;
  }
}

/*==================================
 * bezout(a,b,pu,pv)
 *==================================
 *    Return g = gcd(a,b) >= 0, and assign GENs u,v through pointers pu,pv
 *    such that g = u*a + v*b.
 * Special cases:
 *    a == b == 0 ==> pick u=1, v=0
 *    a != 0 == b ==> keep v=0
 *    a == 0 != b ==> keep u=0
 *    |a| == |b| != 0 ==> keep u=0, set v=+-1
 * Assignments through pu,pv will be suppressed when the corresponding
 * pointer is NULL  (but the computations will happen nonetheless).
 */

GEN
bezout(GEN a, GEN b, GEN *pu, GEN *pv)
{
  GEN t,r;
  GEN *pt;
  long s, sa, sb;
  ulong g;
  ulong xu,xu1,xv,xv1;		/* Lehmer stage recurrence matrix */

  if (typ(a) != t_INT || typ(b) != t_INT) pari_err(arither1);
  s = absi_cmp(a,b);
  if (s < 0)
  {
    t=b; b=a; a=t;
    pt=pu; pu=pv; pv=pt;
  }
  /* now |a| >= |b| */

  sa = signe(a); sb = signe(b);
  if (!sb)
  {
    if (pv) *pv = gen_0;
    switch(sa)
    {
    case  0: if (pu) *pu = gen_0; return gen_0;
    case  1: if (pu) *pu = gen_1; return icopy(a);
    case -1: if (pu) *pu = gen_m1; return(negi(a));
    }
  }
  if (s == 0)			/* |a| == |b| != 0 */
  {
    if (pu) *pu = gen_0;
    if (sb > 0)
    { if (pv) *pv = gen_1; return icopy(b); }
    else
    { if (pv) *pv = gen_m1; return(negi(b)); }
  }
  /* now |a| > |b| > 0 */

  if (lgefint(a) == 3)		/* single-word affair */
  {
    g = xxgcduu((ulong)a[2], (ulong)b[2], 0, &xu, &xu1, &xv, &xv1, &s);
    sa = s > 0 ? sa : -sa;
    sb = s > 0 ? -sb : sb;
    if (pu)
    {
      if (xu == 0) *pu = gen_0; /* can happen when b divides a */
      else if (xu == 1) *pu = sa < 0 ? gen_m1 : gen_1;
      else if (xu == 2) *pu = sa < 0 ? negi(gen_2) : gen_2;
      else
      {
	*pu = cgeti(3);
	(*pu)[1] = evalsigne(sa)|evallgefint(3);
	(*pu)[2] = xu;
      }
    }
    if (pv)
    {
      if (xv == 1) *pv = sb < 0 ? gen_m1 : gen_1;
      else if (xv == 2) *pv = sb < 0 ? negi(gen_2) : gen_2;
      else
      {
	*pv = cgeti(3);
	(*pv)[1] = evalsigne(sb)|evallgefint(3);
	(*pv)[2] = xv;
      }
    }
    if (g == 1) return gen_1;
    else if (g == 2) return gen_2;
    else
    {
      r = cgeti(3);
      r[1] = evalsigne(1)|evallgefint(3);
      r[2] = g;
      return r;
    }
  }
  else
  { /* general case */
    pari_sp av = avma;
    /*Copy serves two purposes:
     * 1) mpn_gcdext destroys its input and needs an extra limb
     * 2) allows us to use icopy instead of gerepile later.
     * NOTE: we must put u before d else the final icopy could fail. */
    GEN ca = icopy_ef(a,lgefint(a)+1);
    GEN cb = icopy_ef(b,lgefint(b)+1);
    GEN u = cgeti(lgefint(a)+1), v = NULL;
    GEN d = cgeti(lgefint(a)+1);
    long su,l,lu;
    l = mpn_gcdext(LIMBS(d), LIMBS(u), &lu, LIMBS(ca), NLIMBS(ca), LIMBS(cb), NLIMBS(cb));
    if (lu<=0)
    {
      if (lu==0) su=0;
      else {su=-1;lu=-lu;}
    }
    else
      su=1;
    if (sa<0) su=-su;
    d[1] = evalsigne(1)|evallgefint(l+2);
    u[1] = evalsigne(su)|evallgefint(lu+2);
    if (pv) v=diviiexact(subii(d,mulii(u,a)),b);
    avma = av;
    if (pu) *pu=icopy(u);
    if (pv) *pv=icopy(v);
    return icopy(d);
  }
}
