#line 2 "../src/kernel/none/add.c"
/* $Id: add.c 8602 2007-04-18 20:06:59Z kb $

Copyright (C) 2002-2003  The PARI group.

This file is part of the PARI/GP package.

PARI/GP is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation. It is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY WHATSOEVER.

Check the License for details. You should have received a copy of it, along
with the package; see the file 'COPYING'. If not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA. */

/* prototype of positive small ints */
static long pos_s[] = {
  evaltyp(t_INT) | _evallg(3), evalsigne(1) | evallgefint(3), 0 };

/* prototype of negative small ints */
static long neg_s[] = {
  evaltyp(t_INT) | _evallg(3), (long)evalsigne(-1) | evallgefint(3), 0 };

GEN
addss(long x, long y)
{
  if (!x) return stoi(y);
  if (x>0) { pos_s[2] = x; return addsi(y,pos_s); }
  neg_s[2] = -x; return addsi(y,neg_s);
}

INLINE GEN
icopy_sign(GEN x, long sx)
{
  GEN y=icopy(x);
  setsigne(y,sx);
  return y;
}

GEN
addsi_sign(long x, GEN y, long sy)
{
  long sx,ly;
  GEN z;

  if (!x) return icopy_sign(y, sy);
  if (!sy) return stoi(x);
  if (x<0) { sx=-1; x=-x; } else sx=1;
  if (sx==sy)
  {
    z = addsispec(x,y+2, lgefint(y)-2);
    setsigne(z,sy); return z;
  }
  ly=lgefint(y);
  if (ly==3)
  {
    const long d = y[2] - x;
    if (!d) return gen_0;
    z=cgeti(3);
    if (y[2] < 0 || d > 0) {
      z[1] = evalsigne(sy) | evallgefint(3);
      z[2] = d;
    }
    else {
      z[1] = evalsigne(-sy) | evallgefint(3);
      z[2] =-d;
    }
    return z;
  }
  z = subisspec(y+2,x, ly-2);
  setsigne(z,sy); return z;
}

GEN
addii_sign(GEN x, long sx, GEN y, long sy)
{
  long lx,ly;
  GEN z;

  if (!sx) return sy? icopy_sign(y, sy): gen_0;
  if (!sy) return icopy_sign(x, sx);
  lx=lgefint(x); ly=lgefint(y);

  if (sx==sy)
    z = addiispec(x+2,y+2,lx-2,ly-2);
  else
  { /* sx != sy */
    long i = lx - ly;
    if (i==0) /* lx == ly */
    {
      i = absi_cmp_lg(x,y,lx);
      if (!i) return gen_0;
    }
    if (i<0) { sx=sy; swapspec(x,y, lx,ly); } /* ensure |x| >= |y| */
    z = subiispec(x+2,y+2,lx-2,ly-2);
  }
  setsigne(z,sx); return z;
}

INLINE GEN
rcopy_sign(GEN x, long sx) { GEN y = rcopy(x); setsigne(y,sx); return y; }

GEN
addir_sign(GEN x, long sx, GEN y, long sy)
{
  long e,l,ly;
  GEN z;

  if (!sx) return rcopy_sign(y, sy);
  e = expo(y) - expi(x);
  if (!sy)
  {
    if (e > 0) return rcopy_sign(y, sy);
    z = itor(x, 3 + ((-e)>>TWOPOTBITS_IN_LONG));
    setsigne(z, sx); return z;
  }

  ly = lg(y);
  if (e > 0)
  {
    l = ly - (e>>TWOPOTBITS_IN_LONG);
    if (l < 3) return rcopy_sign(y, sy);
  }
  else l = ly + ((-e)>>TWOPOTBITS_IN_LONG)+1;
  z = (GEN)avma;
  y = addrr_sign(itor(x,l), sx, y, sy);
  ly = lg(y); while (ly--) *--z = y[ly];
  avma = (pari_sp)z; return z;
}

GEN
addsr(long x, GEN y)
{
  if (!x) return rcopy(y);
  if (x>0) { pos_s[2]=x; return addir_sign(pos_s, 1, y, signe(y)); }
  neg_s[2] = -x; return addir_sign(neg_s, -1, y, signe(y));
}

GEN
subsr(long x, GEN y)
{
  if (!x) return rcopy_sign(y, -signe(y));
  if (x>0) { pos_s[2]=x; return addir_sign(pos_s, 1, y, -signe(y)); }
  neg_s[2] = -x; return addir_sign(neg_s, -1, y, -signe(y));
}

/* return x + 1, assuming x > 0 is a normalized t_REAL of exponent 0 */
GEN
addrex01(GEN x)
{
  long l = lg(x);
  GEN y = cgetr(l);
  y[1] = evalsigne(1) | evalexpo(1);
  y[2] = HIGHBIT | (((ulong)x[2] & ~HIGHBIT) >> 1);
  shift_right(y, x, 3,l, x[2], 1); 
  return y;
}
/* return x - 1 to same accuracy as x, assuming x > 1 is a normalized t_REAL
 * of exponent 0 */
GEN
subrex01(GEN x)
{
  long i, sh, k, l = lg(x);
  ulong u;
  GEN y = cgetr(l);
  k = 2;
  u = (ulong)x[2] & (~HIGHBIT);
  while (!u) u = x[++k]; /* terminates: x not a power of 2 */
  sh = bfffo(u);
  if (sh)
  { shift_left(y+2, x+k, 0, l-k-1, 0, sh); }
  else
  { for (i = 2; i < l-k+2; i++) y[i] = x[k-2 + i]; }
  for (i = l-k+2; i < l; i++) y[i] = 0;
  y[1] = evalsigne(1) | evalexpo(- ((k-2)*BITS_IN_LONG + sh));
  return y;
}

GEN
addrr_sign(GEN x, long sx, GEN y, long sy)
{
  long lx, ex = expo(x);
  long ly, ey = expo(y), e = ey - ex;
  long i, j, lz, ez, m;
  int extend, f2;
  GEN z;
  LOCAL_OVERFLOW;

  if (!sy)
  {
    if (!sx)
    {
      if (e > 0) ex = ey;
      return real_0_bit(ex);
    }
    if (e > 0) return real_0_bit(ey);
    lz = 3 + ((-e)>>TWOPOTBITS_IN_LONG);
    lx = lg(x); if (lz > lx) lz = lx;
    z = cgetr(lz); while(--lz) z[lz] = x[lz];
    setsigne(z,sx); return z;
  }
  if (!sx)
  {
    if (e < 0) return real_0_bit(ex);
    lz = 3 + (e>>TWOPOTBITS_IN_LONG);
    ly = lg(y); if (lz > ly) lz = ly;
    z = cgetr(lz); while (--lz) z[lz] = y[lz];
    setsigne(z,sy); return z;
  }

  if (e < 0) { z=x; x=y; y=z; ey=ex; i=sx; sx=sy; sy=i; e=-e; }
  /* now ey >= ex */
  lx = lg(x);
  ly = lg(y);
  /* If exponents differ, need to shift one argument, here x. If
   * extend = 1: extension of x,z by m < BIL bits (round to 1 word) */
  /* in this case, lz = lx + d + 1, otherwise lx + d */
  extend = 0;
  if (e)
  {
    long d = e >> TWOPOTBITS_IN_LONG, l = ly-d;
    if (l <= 2) return rcopy_sign(y, sy);
    m = e & (BITS_IN_LONG-1);
    if (l > lx) { lz = lx + d + 1; extend = 1; }
    else        { lz = ly; lx = l; }
    if (m)
    { /* shift x right m bits */
      const pari_sp av = avma;
      const ulong sh = BITS_IN_LONG-m;
      GEN p1 = x; x = new_chunk(lx + lz + 1);
      shift_right2(x,p1,2,lx, 0,m,sh);
      if (extend) x[lx] = p1[lx-1] << sh;
      avma = av; /* HACK: cgetr(lz) will not overwrite x */
    }
  }
  else
  { /* d = 0 */
    m = 0;
    if (lx > ly) lx = ly;
    lz = lx;
  }

  if (sx == sy)
  { /* addition */
    i = lz-1;
    j = lx-1;
    if (extend) {
      ulong garde = addll(x[lx], y[i]);
      if (m < 4) /* don't extend for few correct bits */
        z = cgetr(--lz);
      else
      {
        z = cgetr(lz);
        z[i] = garde;
      }
    }
    else
    {
      z = cgetr(lz);
      z[i] = addll(x[j], y[i]); j--;
    }
    i--;
    for (; j>=2; i--,j--) z[i] = addllx(x[j],y[i]);
    if (overflow)
    {
      z[1] = 1; /* stops since z[1] != 0 */
      for (;;) { z[i] = y[i]+1; if (z[i--]) break; }
      if (i <= 0)
      {
        shift_right(z,z, 2,lz, 1,1);
        z[1] = evalsigne(sx) | evalexpo(ey+1); return z;
      }
    }
    for (; i>=2; i--) z[i] = y[i];
    z[1] = evalsigne(sx) | evalexpo(ey); return z;
  }

  /* subtraction */
  if (e) f2 = 1;
  else
  {
    i = 2; while (i < lx && x[i] == y[i]) i++;
    if (i==lx) return real_0_bit(ey - bit_accuracy(lx));
    f2 = ((ulong)y[i] > (ulong)x[i]);
  }
  /* result is non-zero. f2 = (y > x) */
  i = lz-1; z = cgetr(lz);
  if (f2)
  {
    j = lx-1;
    if (extend) z[i] = subll(y[i], x[lx]);
    else        z[i] = subll(y[i], x[j--]);
    for (i--; j>=2; i--) z[i] = subllx(y[i], x[j--]);
    if (overflow) /* stops since y[1] != 0 */
      for (;;) { z[i] = y[i]-1; if (y[i--]) break; }
    for (; i>=2; i--) z[i] = y[i];
    sx = sy;
  }
  else
  {
    if (extend) z[i] = subll(x[lx], y[i]);
    else        z[i] = subll(x[i],  y[i]);
    for (i--; i>=2; i--) z[i] = subllx(x[i], y[i]);
  }

  x = z+2; i = 0; while (!x[i]) i++;
  lz -= i; z += i;
  j = bfffo(z[2]); /* need to shift left by j bits to normalize mantissa */
  ez = ey - (j | (i<<TWOPOTBITS_IN_LONG));
  if (extend)
  { /* z was extended by d+1 words [should be e bits = d words + m bits] */
    /* not worth keeping extra word if less than 5 significant bits in there */
    if (m - j < 5 && lz > 3)
    { /* shorten z */
      ulong last = (ulong)z[--lz]; /* cancelled word */

      /* if we need to shift anyway, shorten from left
       * If not, shorten from right, neutralizing last word of z */
      if (j == 0)
        /* stackdummy((pari_sp)(z + lz+1), (pari_sp)(z + lz)); */
        z[lz] = evaltyp(t_VECSMALL) | _evallg(1);
      else
      {
        GEN t = z;
        z++; shift_left(z,t,2,lz-1, last,j);
      }
      if ((last<<j) & HIGHBIT)
      { /* round up */
        i = lz-1;
        while (++z[i] == 0 && i > 1) i--;
        if (i == 1) { ez++; z[2] = HIGHBIT; }
      }
    }
    else if (j) shift_left(z,z,2,lz-1, 0,j);
  }
  else if (j) shift_left(z,z,2,lz-1, 0,j);
  z[1] = evalsigne(sx) | evalexpo(ez);
  z[0] = evaltyp(t_REAL) | evallg(lz);
  avma = (pari_sp)z; return z;
}
