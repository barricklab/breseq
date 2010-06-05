#line 2 "../src/kernel/none/mp.c"
/* $Id: mp.c 8000 2006-08-24 21:51:50Z kb $

Copyright (C) 2000-2003 The PARI group.

This file is part of the PARI/GP package.

PARI/GP is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation. It is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY WHATSOEVER.

Check the License for details. You should have received a copy of it, along
with the package; see the file 'COPYING'. If not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA. */

/***********************************************************************/
/**								      **/
/**		         MULTIPRECISION KERNEL           	      **/
/**                                                                   **/
/***********************************************************************/
#include "pari.h"
#include "paripriv.h"
#include "../src/kernel/none/tune-gen.h"

int pari_kernel_init(void) { return 0; } /*nothing to do*/

/* NOTE: arguments of "spec" routines (muliispec, addiispec, etc.) aren't
 * GENs but pairs (long *a, long na) representing a list of digits (in basis
 * BITS_IN_LONG) : a[0], ..., a[na-1]. [ In ordre to facilitate splitting: no
 * need to reintroduce codewords ] */

/* Normalize a non-negative integer */
GEN
int_normalize(GEN x, long known_zero_words)
{
  long lx = lgefint(x);
  long i = 2 + known_zero_words;
  for ( ; i < lx; i++)
    if (x[i])
    {
      if (i != 2)
      {
        GEN x0 = x;
        i -= 2; x += i;
        if (x0 == (GEN)avma) avma = (pari_sp)x;
        else stackdummy((pari_sp)(x0+i), (pari_sp)x0);
        lx -= i;
        x[0] = evaltyp(t_INT) | evallg(lx);
        x[1] = evalsigne(1) | evallgefint(lx);
      }
      return x;
    }
  x[1] = evalsigne(0) | evallgefint(2); return x;
}

/***********************************************************************/
/**								      **/
/**		         ADDITION / SUBTRACTION          	      **/
/**                                                                   **/
/***********************************************************************/

GEN
setloop(GEN a)
{
  GEN z0 = (GEN)avma; (void)cgetg(lgefint(a) + 3, t_VECSMALL);
  return icopy_av(a, z0); /* two cells of extra space before a */
}

/* we had a = setloop(?), then some incloops. Reset a to b */
GEN
resetloop(GEN a, GEN b) {
  long lb = lgefint(b);
  a += lgefint(a) - lb;
  a[0] = evaltyp(t_INT) | evallg(lb);
  affii(b, a); return a;
}

/* assume a > 0, initialized by setloop. Do a++ */
static GEN
incpos(GEN a)
{
  long i, l = lgefint(a);
  for (i=l-1; i>1; i--)
    if (++a[i]) return a;
  l++; a--; /* use extra cell */
  a[0]=evaltyp(t_INT) | _evallg(l);
  a[1]=evalsigne(1) | evallgefint(l);
  a[2]=1; return a;
}

/* assume a < 0, initialized by setloop. Do a++ */
static GEN
incneg(GEN a)
{
  long l = lgefint(a)-1;
  if (a[l]--)
  {
    if (l == 2 && !a[2])
    {
      a++; /* save one cell */
      a[0] = evaltyp(t_INT) | _evallg(2);
      a[1] = evalsigne(0) | evallgefint(2);
    }
    return a;
  }
  for (l--; l>1; l--)
    if (a[l]--) break;
  l++; a++; /* save one cell */
  a[0] = evaltyp(t_INT) | _evallg(l);
  a[1] = evalsigne(-1) | evallgefint(l);
  return a;
}

/* assume a initialized by setloop. Do a++ */
GEN
incloop(GEN a)
{
  switch(signe(a))
  {
    case 0: a--; /* use extra cell */
      a[0]=evaltyp(t_INT) | _evallg(3);
      a[1]=evalsigne(1) | evallgefint(3);
      a[2]=1; return a;
    case -1: return incneg(a);
    default: return incpos(a);
  }
}

INLINE GEN
addsispec(long s, GEN x, long nx)
{
  GEN xd, zd = (GEN)avma;
  long lz;

  lz = nx+3; (void)new_chunk(lz);
  xd = x + nx;
  *--zd = *--xd + s;
  if ((ulong)*zd < (ulong)s)
    for(;;)
    {
      if (xd == x) { *--zd = 1; break; } /* enlarge z */
      *--zd = ((ulong)*--xd) + 1;
      if (*zd) { lz--; break; }
    }
  else lz--;
  while (xd > x) *--zd = *--xd;
  *--zd = evalsigne(1) | evallgefint(lz);
  *--zd = evaltyp(t_INT) | evallg(lz);
  avma=(pari_sp)zd; return zd;
}

static GEN
addiispec(GEN x, GEN y, long nx, long ny)
{
  GEN xd,yd,zd;
  long lz;
  LOCAL_OVERFLOW;

  if (nx < ny) swapspec(x,y, nx,ny);
  if (ny == 1) return addsispec(*y,x,nx);
  zd = (GEN)avma;
  lz = nx+3; (void)new_chunk(lz);
  xd = x + nx;
  yd = y + ny;
  *--zd = addll(*--xd, *--yd);
  while (yd > y) *--zd = addllx(*--xd, *--yd);
  if (overflow)
    for(;;)
    {
      if (xd == x) { *--zd = 1; break; } /* enlarge z */
      *--zd = ((ulong)*--xd) + 1;
      if (*zd) { lz--; break; }
    }
  else lz--;
  while (xd > x) *--zd = *--xd;
  *--zd = evalsigne(1) | evallgefint(lz);
  *--zd = evaltyp(t_INT) | evallg(lz);
  avma=(pari_sp)zd; return zd;
}

/* assume x >= y */
INLINE GEN
subisspec(GEN x, long s, long nx)
{
  GEN xd, zd = (GEN)avma;
  long lz;
  LOCAL_OVERFLOW;

  lz = nx+2; (void)new_chunk(lz);
  xd = x + nx;
  *--zd = subll(*--xd, s);
  if (overflow)
    for(;;)
    {
      *--zd = ((ulong)*--xd) - 1;
      if (*xd) break;
    }
  if (xd == x)
    while (*zd == 0) { zd++; lz--; } /* shorten z */
  else
    do  *--zd = *--xd; while (xd > x);
  *--zd = evalsigne(1) | evallgefint(lz);
  *--zd = evaltyp(t_INT) | evallg(lz);
  avma=(pari_sp)zd; return zd;
}

/* assume x > y */
static GEN
subiispec(GEN x, GEN y, long nx, long ny)
{
  GEN xd,yd,zd;
  long lz;
  LOCAL_OVERFLOW;

  if (ny==1) return subisspec(x,*y,nx);
  zd = (GEN)avma;
  lz = nx+2; (void)new_chunk(lz);
  xd = x + nx;
  yd = y + ny;
  *--zd = subll(*--xd, *--yd);
  while (yd > y) *--zd = subllx(*--xd, *--yd);
  if (overflow)
    for(;;)
    {
      *--zd = ((ulong)*--xd) - 1;
      if (*xd) break;
    }
  if (xd == x)
    while (*zd == 0) { zd++; lz--; } /* shorten z */
  else
    do  *--zd = *--xd; while (xd > x);
  *--zd = evalsigne(1) | evallgefint(lz);
  *--zd = evaltyp(t_INT) | evallg(lz);
  avma=(pari_sp)zd; return zd;
}

static void
roundr_up_ip(GEN x, long l)
{
  long i = l;
  for(;;)
  {
    if (++x[--i]) break;
    if (i == 2) { x[2] = HIGHBIT; setexpo(x, expo(x)+1); break; }
  }
}

void
affir(GEN x, GEN y)
{
  const long s = signe(x), ly = lg(y);
  long lx, sh, i;

  if (!s)
  {
    y[1] = evalexpo(-bit_accuracy(ly));
    return;
  }

  lx = lgefint(x); sh = bfffo(x[2]);
  y[1] = evalsigne(s) | evalexpo(bit_accuracy(lx)-sh-1);
  if (sh) {
    if (lx <= ly)
    {
      for (i=lx; i<ly; i++) y[i]=0;
      shift_left(y,x,2,lx-1, 0,sh);
      return;
    }
    shift_left(y,x,2,ly-1, x[ly],sh);
    /* lx > ly: round properly */
    if ((x[ly]<<sh) & HIGHBIT) roundr_up_ip(y, ly);
  }
  else {
    if (lx <= ly)
    {
      for (i=2; i<lx; i++) y[i]=x[i];
      for (   ; i<ly; i++) y[i]=0;
      return;
    }
    for (i=2; i<ly; i++) y[i]=x[i];
    /* lx > ly: round properly */
    if (x[ly] & HIGHBIT) roundr_up_ip(y, ly);
  }
}

static GEN
shifti_spec(GEN x, long lx, long n)
{
  long ly, i, m, s = signe(x);
  GEN y;
  if (!s) return gen_0;
  if (!n)
  {
    y = cgeti(lx);
    y[1] = evalsigne(s) | evallgefint(lx);
    while (--lx > 1) y[lx]=x[lx];
    return y;
  }
  if (n > 0)
  {
    GEN z = (GEN)avma;
    long d = n>>TWOPOTBITS_IN_LONG;

    ly = lx+d; y = new_chunk(ly);
    for ( ; d; d--) *--z = 0;
    m = n & (BITS_IN_LONG-1);
    if (!m) for (i=2; i<lx; i++) y[i]=x[i];
    else
    {
      register const ulong sh = BITS_IN_LONG - m;
      shift_left2(y,x, 2,lx-1, 0,m,sh);
      i = ((ulong)x[2]) >> sh;
      /* Extend y on the left? */
      if (i) { ly++; y = new_chunk(1); y[2] = i; }
    }
  }
  else
  {
    n = -n;
    ly = lx - (n>>TWOPOTBITS_IN_LONG);
    if (ly<3) return gen_0;
    y = new_chunk(ly);
    m = n & (BITS_IN_LONG-1);
    if (m) {
      shift_right(y,x, 2,ly, 0,m);
      if (y[2] == 0)
      {
        if (ly==3) { avma = (pari_sp)(y+3); return gen_0; }
        ly--; avma = (pari_sp)(++y);
      }
    } else {
      for (i=2; i<ly; i++) y[i]=x[i];
    }
  }
  y[1] = evalsigne(s)|evallgefint(ly);
  y[0] = evaltyp(t_INT)|evallg(ly); return y;
}

GEN
shifti(GEN x, long n)
{
  return shifti_spec(x, lgefint(x), n);
}

GEN
ishiftr_lg(GEN x, long lx, long n)
{ /*This is a kludge since x is not an integer*/
  return shifti_spec(x, lx, n);
}

GEN
truncr(GEN x)
{
  long d,e,i,s,m;
  GEN y;

  if ((s=signe(x)) == 0 || (e=expo(x)) < 0) return gen_0;
  d = (e>>TWOPOTBITS_IN_LONG) + 3;
  m = e & (BITS_IN_LONG-1);
  if (d > lg(x)) pari_err(precer, "truncr (precision loss in truncation)");

  y=cgeti(d); y[1] = evalsigne(s) | evallgefint(d);
  if (++m == BITS_IN_LONG)
    for (i=2; i<d; i++) y[i]=x[i];
  else
  {
    register const ulong sh = BITS_IN_LONG - m;
    shift_right2(y,x, 2,d,0, sh,m);
  }
  return y;
}

/* integral part */
GEN
floorr(GEN x)
{
  long d,e,i,lx,m;
  GEN y;

  if (signe(x) >= 0) return truncr(x);
  if ((e=expo(x)) < 0) return gen_m1;
  d = (e>>TWOPOTBITS_IN_LONG) + 3;
  m = e & (BITS_IN_LONG-1);
  lx=lg(x); if (d>lx) pari_err(precer, "floorr (precision loss in truncation)");
  y = new_chunk(d);
  if (++m == BITS_IN_LONG)
  {
    for (i=2; i<d; i++) y[i]=x[i];
    i=d; while (i<lx && !x[i]) i++;
    if (i==lx) goto END;
  }
  else
  {
    register const ulong sh = BITS_IN_LONG - m;
    shift_right2(y,x, 2,d,0, sh,m);
    if (x[d-1]<<m == 0)
    {
      i=d; while (i<lx && !x[i]) i++;
      if (i==lx) goto END;
    }
  }
  /* set y:=y+1 */
  for (i=d-1; i>=2; i--) { y[i]++; if (y[i]) goto END; }
  y=new_chunk(1); y[2]=1; d++;
END:
  y[1] = evalsigne(-1) | evallgefint(d);
  y[0] = evaltyp(t_INT) | evallg(d); return y;
}

INLINE int
absi_cmp_lg(GEN x, GEN y, long l)
{
  long i=2;
  while (i<l && x[i]==y[i]) i++;
  if (i==l) return 0;
  return ((ulong)x[i] > (ulong)y[i])? 1: -1;
}

INLINE int
absi_equal_lg(GEN x, GEN y, long l)
{
  long i = l-1; while (i>1 && x[i]==y[i]) i--;
  return i==1;
}

/***********************************************************************/
/**								      **/
/**		          MULTIPLICATION                 	      **/
/**                                                                   **/
/***********************************************************************/
GEN
mulss(long x, long y)
{
  long s,p1;
  GEN z;
  LOCAL_HIREMAINDER;

  if (!x || !y) return gen_0;
  if (x<0) { s = -1; x = -x; } else s=1;
  if (y<0) { s = -s; y = -y; }
  p1 = mulll(x,y);
  if (hiremainder)
  {
    z=cgeti(4); z[1] = evalsigne(s) | evallgefint(4);
    z[2]=hiremainder; z[3]=p1; return z;
  }
  z=cgeti(3); z[1] = evalsigne(s) | evallgefint(3);
  z[2]=p1; return z;
}

GEN
muluu(ulong x, ulong y)
{
  long p1;
  GEN z;
  LOCAL_HIREMAINDER;

  if (!x || !y) return gen_0;
  p1 = mulll(x,y);
  if (hiremainder)
  {
    z=cgeti(4); z[1] = evalsigne(1) | evallgefint(4);
    z[2]=hiremainder; z[3]=p1; return z;
  }
  z=cgeti(3); z[1] = evalsigne(1) | evallgefint(3);
  z[2]=p1; return z;
}

/* assume ny > 0 */
INLINE GEN
muluispec(ulong x, GEN y, long ny)
{
  GEN yd, z = (GEN)avma;
  long lz = ny+3;
  LOCAL_HIREMAINDER;

  (void)new_chunk(lz);
  yd = y + ny; *--z = mulll(x, *--yd);
  while (yd > y) *--z = addmul(x,*--yd);
  if (hiremainder) *--z = hiremainder; else lz--;
  *--z = evalsigne(1) | evallgefint(lz);
  *--z = evaltyp(t_INT) | evallg(lz);
  avma=(pari_sp)z; return z;
}

/* a + b*|Y| */
GEN
addumului(ulong a, ulong b, GEN Y)
{
  GEN yd,y,z;
  long ny,lz;
  LOCAL_HIREMAINDER;
  LOCAL_OVERFLOW;

  if (!signe(Y)) return utoi(a);

  y = Y+2; z = (GEN)avma;
  ny = lgefint(Y)-2;
  lz = ny+3;

  (void)new_chunk(lz);
  yd = y + ny; *--z = addll(a, mulll(b, *--yd));
  if (overflow) hiremainder++; /* can't overflow */
  while (yd > y) *--z = addmul(b,*--yd);
  if (hiremainder) *--z = hiremainder; else lz--;
  *--z = evalsigne(1) | evallgefint(lz);
  *--z = evaltyp(t_INT) | evallg(lz);
  avma=(pari_sp)z; return z;
}

GEN muliispec(GEN a, GEN b, long na, long nb);
/*#define KARAMULR_VARIANT*/
#define muliispec_mirror muliispec

/***********************************************************************/
/**								      **/
/**		          DIVISION                       	      **/
/**                                                                   **/
/***********************************************************************/

ulong
umodiu(GEN y, ulong x)
{
  long sy=signe(y),ly,i;
  LOCAL_HIREMAINDER;

  if (!x) pari_err(gdiver);
  if (!sy) return 0;
  ly = lgefint(y);
  if (x <= (ulong)y[2]) hiremainder=0;
  else
  {
    if (ly==3) return (sy > 0)? (ulong)y[2]: x - (ulong)y[2];
    hiremainder=y[2]; ly--; y++;
  }
  for (i=2; i<ly; i++) (void)divll(y[i],x);
  if (!hiremainder) return 0;
  return (sy > 0)? hiremainder: x - hiremainder;
}

/* return |y| \/ x */
GEN
diviu_rem(GEN y, ulong x, ulong *rem)
{
  long ly,i;
  GEN z;
  LOCAL_HIREMAINDER;

  if (!x) pari_err(gdiver);
  if (!signe(y)) { *rem = 0; return gen_0; }

  ly = lgefint(y);
  if (x <= (ulong)y[2]) hiremainder=0;
  else
  {
    if (ly==3) { *rem = (ulong)y[2]; return gen_0; }
    hiremainder=y[2]; ly--; y++;
  }
  z = cgeti(ly); z[1] = evallgefint(ly) | evalsigne(1);
  for (i=2; i<ly; i++) z[i]=divll(y[i],x);
  *rem = hiremainder; return z;
}

GEN
divis_rem(GEN y, long x, long *rem)
{
  long sy=signe(y),ly,s,i;
  GEN z;
  LOCAL_HIREMAINDER;

  if (!x) pari_err(gdiver);
  if (!sy) { *rem=0; return gen_0; }
  if (x<0) { s = -sy; x = -x; } else s = sy;

  ly = lgefint(y);
  if ((ulong)x <= (ulong)y[2]) hiremainder=0;
  else
  {
    if (ly==3) { *rem = itos(y); return gen_0; }
    hiremainder=y[2]; ly--; y++;
  }
  z = cgeti(ly); z[1] = evallgefint(ly) | evalsigne(s);
  for (i=2; i<ly; i++) z[i]=divll(y[i],x);
  if (sy<0) hiremainder = - ((long)hiremainder);
  *rem = (long)hiremainder; return z;
}

GEN
divis(GEN y, long x)
{
  long sy=signe(y),ly,s,i;
  GEN z;
  LOCAL_HIREMAINDER;

  if (!x) pari_err(gdiver);
  if (!sy) return gen_0;
  if (x<0) { s = -sy; x = -x; } else s = sy;

  ly = lgefint(y);
  if ((ulong)x <= (ulong)y[2]) hiremainder=0;
  else
  {
    if (ly==3) return gen_0;
    hiremainder=y[2]; ly--; y++;
  }
  z = cgeti(ly); z[1] = evallgefint(ly) | evalsigne(s);
  for (i=2; i<ly; i++) z[i]=divll(y[i],x);
  return z;
}

GEN
divrr(GEN x, GEN y)
{
  long sx=signe(x), sy=signe(y), lx,ly,lr,e,i,j;
  ulong y0,y1;
  GEN r, r1;

  if (!sy) pari_err(gdiver);
  e = expo(x) - expo(y);
  if (!sx) return real_0_bit(e);
  if (sy<0) sx = -sx;

  lx=lg(x); ly=lg(y);
  if (ly==3)
  {
    ulong k = x[2], l = (lx>3)? x[3]: 0;
    LOCAL_HIREMAINDER;
    if (k < (ulong)y[2]) e--;
    else
    {
      l >>= 1; if (k&1) l |= HIGHBIT;
      k >>= 1;
    }
    r = cgetr(3); r[1] = evalsigne(sx) | evalexpo(e);
    hiremainder=k; r[2]=divll(l,y[2]); return r;
  }

  lr = min(lx,ly); r = new_chunk(lr);
  r1 = r-1;
  r1[1] = 0; for (i=2; i<lr; i++) r1[i]=x[i];
  r1[lr] = (lx>ly)? x[lr]: 0;
  y0 = y[2]; y1 = y[3];
  for (i=0; i<lr-1; i++)
  { /* r1 = r + (i-1) */
    ulong k, qp;
    LOCAL_HIREMAINDER;
    LOCAL_OVERFLOW;

    if ((ulong)r1[1] == y0)
    {
      qp = MAXULONG; k = addll(y0,r1[2]);
    }
    else
    {
      if ((ulong)r1[1] > y0) /* can't happen if i=0 */
      {
        GEN y1 = y+1;
        j = lr-i; r1[j] = subll(r1[j],y1[j]);
	for (j--; j>0; j--) r1[j] = subllx(r1[j],y1[j]);
	j=i; do r[--j]++; while (j && !r[j]);
      }
      hiremainder = r1[1]; overflow = 0;
      qp = divll(r1[2],y0); k = hiremainder;
    }
    if (!overflow)
    {
      long k3 = subll(mulll(qp,y1), r1[3]);
      long k4 = subllx(hiremainder,k);
      while (!overflow && k4) { qp--; k3 = subll(k3,y1); k4 = subllx(k4,y0); }
    }
    j = lr-i+1;
    if (j<ly) (void)mulll(qp,y[j]); else { hiremainder = 0 ; j = ly; }
    for (j--; j>1; j--)
    {
      r1[j] = subll(r1[j], addmul(qp,y[j]));
      hiremainder += overflow;
    }
    if ((ulong)r1[1] != hiremainder)
    {
      if ((ulong)r1[1] < hiremainder)
      {
        qp--;
        j = lr-i-(lr-i>=ly); r1[j] = addll(r1[j], y[j]);
        for (j--; j>1; j--) r1[j] = addllx(r1[j], y[j]);
      }
      else
      {
	r1[1] -= hiremainder;
	while (r1[1])
	{
	  qp++; if (!qp) { j=i; do r[--j]++; while (j && !r[j]); }
          j = lr-i-(lr-i>=ly); r1[j] = subll(r1[j],y[j]);
          for (j--; j>1; j--) r1[j] = subllx(r1[j],y[j]);
	  r1[1] -= overflow;
	}
      }
    }
    *++r1 = qp;
  }
  /* i = lr-1 */
  /* round correctly */
  if ((ulong)r1[1] > (y0>>1))
  {
    j=i; do r[--j]++; while (j && !r[j]);
  }
  r1 = r-1; for (j=i; j>=2; j--) r[j]=r1[j];
  if (r[0] == 0) e--;
  else if (r[0] == 1) { shift_right(r,r, 2,lr, 1,1); }
  else { /* possible only when rounding up to 0x2 0x0 ... */
    r[2] = HIGHBIT; e++;
  }
  r[0] = evaltyp(t_REAL)|evallg(lr);
  r[1] = evalsigne(sx) | evalexpo(e);
  return r;
}

GEN
divri(GEN x, GEN y)
{
  long lx, s = signe(y);
  pari_sp av;
  GEN z;

  if (!s) pari_err(gdiver);
  if (!signe(x)) return real_0_bit(expo(x) - expi(y));
  if (!is_bigint(y)) return divrs(x, s>0? y[2]: -y[2]);

  lx = lg(x); z = cgetr(lx); av = avma;
  affrr(divrr(x, itor(y, lx+1)), z);
  avma = av; return z;
}

/* Integer division x / y: such that sign(r) = sign(x)
 *   if z = ONLY_REM return remainder, otherwise return quotient
 *   if z != NULL set *z to remainder
 *   *z is the last object on stack (and thus can be disposed of with cgiv
 *   instead of gerepile)
 * If *z is zero, we put gen_0 here and no copy.
 * space needed: lx + ly */
GEN
dvmdii(GEN x, GEN y, GEN *z)
{
  long sx=signe(x),sy=signe(y);
  long lx, ly, lz, i, j, sh, lq, lr;
  pari_sp av;
  ulong y0,y1, *xd,*rd,*qd;
  GEN q, r, r1;

  if (!sy) { if (z == ONLY_REM && !sx) return gen_0; pari_err(gdiver); }
  if (!sx)
  {
    if (!z || z == ONLY_REM) return gen_0;
    *z=gen_0; return gen_0;
  }
  lx=lgefint(x);
  ly=lgefint(y); lz=lx-ly;
  if (lz <= 0)
  {
    if (lz == 0)
    {
      for (i=2; i<lx; i++)
        if (x[i] != y[i])
        {
          if ((ulong)x[i] > (ulong)y[i]) goto DIVIDE;
          goto TRIVIAL;
        }
      if (z == ONLY_REM) return gen_0;
      if (z) *z = gen_0;
      if (sx < 0) sy = -sy;
      return stoi(sy);
    }
TRIVIAL:
    if (z == ONLY_REM) return icopy(x);
    if (z) *z = icopy(x);
    return gen_0;
  }
DIVIDE: /* quotient is non-zero */
  av=avma; if (sx<0) sy = -sy;
  if (ly==3)
  {
    LOCAL_HIREMAINDER;
    y0 = y[2];
    if (y0 <= (ulong)x[2]) hiremainder=0;
    else
    {
      hiremainder = x[2]; lx--; x++;
    }
    q = new_chunk(lx); for (i=2; i<lx; i++) q[i]=divll(x[i],y0);
    if (z == ONLY_REM)
    {
      avma=av; if (!hiremainder) return gen_0;
      r=cgeti(3);
      r[1] = evalsigne(sx) | evallgefint(3);
      r[2]=hiremainder; return r;
    }
    q[1] = evalsigne(sy) | evallgefint(lx);
    q[0] = evaltyp(t_INT) | evallg(lx);
    if (!z) return q;
    if (!hiremainder) { *z=gen_0; return q; }
    r=cgeti(3);
    r[1] = evalsigne(sx) | evallgefint(3);
    r[2] = hiremainder; *z=r; return q;
  }

  r1 = new_chunk(lx); sh = bfffo(y[2]);
  if (sh)
  { /* normalize so that highbit(y) = 1 (shift left x and y by sh bits)*/
    register const ulong m = BITS_IN_LONG - sh;
    r = new_chunk(ly);
    shift_left2(r, y,2,ly-1, 0,sh,m); y = r;
    shift_left2(r1,x,2,lx-1, 0,sh,m);
    r1[1] = ((ulong)x[2]) >> m;
  }
  else
  {
    r1[1] = 0; for (j=2; j<lx; j++) r1[j] = x[j];
  }
  x = r1;
  y0 = y[2]; y1 = y[3];
  for (i=0; i<=lz; i++)
  { /* r1 = x + i */
    ulong k, qp;
    LOCAL_HIREMAINDER;
    LOCAL_OVERFLOW;

    if ((ulong)r1[1] == y0)
    {
      qp = MAXULONG; k = addll(y0,r1[2]);
    }
    else
    {
      hiremainder = r1[1]; overflow = 0;
      qp = divll(r1[2],y0); k = hiremainder;
    }
    if (!overflow)
    {
      long k3 = subll(mulll(qp,y1), r1[3]);
      long k4 = subllx(hiremainder,k);
      while (!overflow && k4) { qp--; k3 = subll(k3,y1); k4 = subllx(k4,y0); }
    }
    hiremainder = 0; j = ly;
    for (j--; j>1; j--)
    {
      r1[j] = subll(r1[j], addmul(qp,y[j]));
      hiremainder += overflow;
    }
    if ((ulong)r1[1] < hiremainder)
    {
      qp--;
      j = ly-1; r1[j] = addll(r1[j],y[j]);
      for (j--; j>1; j--) r1[j] = addllx(r1[j],y[j]);
    }
    *++r1 = qp;
  }

  lq = lz+2;
  if (!z)
  {
    qd = (ulong*)av;
    xd = (ulong*)(x + lq);
    if (x[1]) { lz++; lq++; }
    while (lz--) *--qd = *--xd;
    *--qd = evalsigne(sy) | evallgefint(lq);
    *--qd = evaltyp(t_INT) | evallg(lq);
    avma = (pari_sp)qd; return (GEN)qd;
  }

  j=lq; while (j<lx && !x[j]) j++;
  lz = lx-j;
  if (z == ONLY_REM)
  {
    if (lz==0) { avma = av; return gen_0; }
    rd = (ulong*)av; lr = lz+2;
    xd = (ulong*)(x + lx);
    if (!sh) while (lz--) *--rd = *--xd;
    else
    { /* shift remainder right by sh bits */
      const ulong shl = BITS_IN_LONG - sh;
      ulong l;
      xd--;
      while (--lz) /* fill r[3..] */
      {
        l = *xd >> sh;
        *--rd = l | (*--xd << shl);
      }
      l = *xd >> sh;
      if (l) *--rd = l; else lr--;
    }
    *--rd = evalsigne(sx) | evallgefint(lr);
    *--rd = evaltyp(t_INT) | evallg(lr);
    avma = (pari_sp)rd; return (GEN)rd;
  }

  lr = lz+2;
  rd = NULL; /* gcc -Wall */
  if (lz)
  { /* non zero remainder: initialize rd */
    xd = (ulong*)(x + lx);
    if (!sh)
    {
      rd = (ulong*)avma; (void)new_chunk(lr);
      while (lz--) *--rd = *--xd;
    }
    else
    { /* shift remainder right by sh bits */
      const ulong shl = BITS_IN_LONG - sh;
      ulong l;
      rd = (ulong*)x; /* overwrite shifted y */
      xd--;
      while (--lz)
      {
        l = *xd >> sh;
        *--rd = l | (*--xd << shl);
      }
      l = *xd >> sh;
      if (l) *--rd = l; else lr--;
    }
    *--rd = evalsigne(sx) | evallgefint(lr);
    *--rd = evaltyp(t_INT) | evallg(lr);
    rd += lr;
  }
  qd = (ulong*)av;
  xd = (ulong*)(x + lq);
  if (x[1]) lq++;
  j = lq-2; while (j--) *--qd = *--xd;
  *--qd = evalsigne(sy) | evallgefint(lq);
  *--qd = evaltyp(t_INT) | evallg(lq);
  q = (GEN)qd;
  if (lr==2) *z = gen_0;
  else
  { /* rd has been properly initialized: we had lz > 0 */
    while (lr--) *--qd = *--rd;
    *z = (GEN)qd;
  }
  avma = (pari_sp)qd; return q;
}

/* Montgomery reduction.
 * N has k words, assume T >= 0 has less than 2k.
 * Return res := T / B^k mod N, where B = 2^BIL
 * such that 0 <= res < T/B^k + N  and  res has less than k words */
GEN
red_montgomery(GEN T, GEN N, ulong inv)
{
  pari_sp av;
  GEN Te, Td, Ne, Nd, scratch;
  ulong i, j, m, t, d, k = lgefint(N)-2;
  int carry;
  LOCAL_HIREMAINDER;
  LOCAL_OVERFLOW;

  if (k == 0) return gen_0;
  d = lgefint(T)-2; /* <= 2*k */
#ifdef DEBUG
  if (d > 2*k) pari_err(bugparier,"red_montgomery");
#endif
  if (k == 1)
  { /* as below, special cased for efficiency */
    ulong n = (ulong)N[2];
    t = (ulong)T[d+1];
    m = t * inv;
    (void)addll(mulll(m, n), t); /* = 0 */
    t = hiremainder + overflow;
    if (d == 2)
    {
      t = addll((ulong)T[2], t);
      if (overflow) t -= n; /* t > n doesn't fit in 1 word */
    }
    return utoi(t);
  }
  /* assume k >= 2 */
  av = avma; scratch = new_chunk(k<<1); /* >= k + 2: result fits */

  /* copy T to scratch space (pad with zeroes to 2k words) */
  Td = (GEN)av;
  Te = T + (d+2);
  for (i=0; i < d     ; i++) *--Td = *--Te;
  for (   ; i < (k<<1); i++) *--Td = 0;

  Te = (GEN)av; /* 1 beyond end of T mantissa */
  Ne = N + k+2; /* 1 beyond end of N mantissa */

  carry = 0;
  for (i=0; i<k; i++) /* set T := T/B nod N, k times */
  {
    Td = Te; /* one beyond end of (new) T mantissa */
    Nd = Ne;
    m = *--Td * inv; /* solve T + m N = O(B) */

    /* set T := (T + mN) / B */
    Te = Td;
    (void)addll(mulll(m, *--Nd), *Td); /* = 0 */
    for (j=1; j<k; j++)
    {
      hiremainder += overflow;
      t = addll(addmul(m, *--Nd), *--Td); *Td = t;
    }
    overflow += hiremainder;
    t = addll(overflow, *--Td); *Td = t + carry;
    carry = (overflow || (carry && *Td == 0));
  }
  if (carry)
  { /* Td > N overflows (k+1 words), set Td := Td - N */
    Td = Te;
    Nd = Ne;
    t = subll(*--Td, *--Nd); *Td = t;
    while (Td > scratch) { t = subllx(*--Td, *--Nd); *Td = t; }
  }

  /* copy result */
  Td = (GEN)av;
  while (! *scratch && Te > scratch) scratch++; /* strip leading 0s */
  while (Te > scratch) *--Td = *--Te;
  k = (GEN)av - Td; if (!k) return gen_0;
  k += 2;
  *--Td = evalsigne(1) | evallgefint(k);
  *--Td = evaltyp(t_INT) | evallg(k);
#ifdef DEBUG
{
  long l = lgefint(N)-2, s = BITS_IN_LONG*l;
  GEN R = int2n(s);
  GEN res = remii(mulii(T, Fp_inv(R, N)), N);
  if (k > lgefint(N)
    || !equalii(remii(Td,N),res)
    || cmpii(Td, addii(shifti(T, -s), N)) >= 0) pari_err(bugparier,"red_montgomery");
}
#endif
  avma = (pari_sp)Td; return Td;
}

/* EXACT INTEGER DIVISION */

/* assume xy>0, the division is exact and y is odd. Destroy x */
static GEN
diviuexact_i(GEN x, ulong y)
{
  long i, lz, lx;
  ulong q, yinv;
  GEN z, z0, x0, x0min;

  if (y == 1) return icopy(x);
  lx = lgefint(x);
  if (lx == 3) return utoipos((ulong)x[2] / y);
  yinv = invrev(y);
  lz = (y <= (ulong)x[2]) ? lx : lx-1;
  z = new_chunk(lz);
  z0 = z + lz;
  x0 = x + lx; x0min = x + lx-lz+2;

  while (x0 > x0min)
  {
    *--z0 = q = yinv*((ulong)*--x0); /* i-th quotient */
    if (!q) continue;
    /* x := x - q * y */
    { /* update neither lowest word (could set it to 0) nor highest ones */
      register GEN x1 = x0 - 1;
      LOCAL_HIREMAINDER;
      (void)mulll(q,y);
      if (hiremainder)
      {
        if ((ulong)*x1 < hiremainder)
        {
          *x1 -= hiremainder;
          do (*--x1)--; while ((ulong)*x1 == MAXULONG);
        }
        else
          *x1 -= hiremainder;
      }
    }
  }
  i=2; while(!z[i]) i++;
  z += i-2; lz -= i-2;
  z[0] = evaltyp(t_INT)|evallg(lz);
  z[1] = evalsigne(1)|evallg(lz);
  avma = (pari_sp)z; return z;
}

/* assume y != 0 and the division is exact */
GEN
diviuexact(GEN x, ulong y)
{
  pari_sp av;
  long lx, vy, s = signe(x);
  GEN z;
  
  if (!s) return gen_0;
  if (y == 1) return icopy(x);
  lx = lgefint(x);
  if (lx == 3) {
    ulong q = (ulong)x[2] / y;
    return (s > 0)? utoipos(q): utoineg(q);
  }
  av = avma; (void)new_chunk(lx); vy = vals(y);
  if (vy) { 
    y >>= vy;
    if (y == 1) { avma = av; return shifti(x, -vy); }
    x = shifti(x, -vy);
    if (lx == 3) {
      ulong q = (ulong)x[2] / y;
      avma = av;
      return (s > 0)? utoipos(q): utoineg(q);
    }
  } else x = icopy(x);
  avma = av;
  z = diviuexact_i(x, y);
  setsigne(z, s); return z;
}

/* Find z such that x=y*z, knowing that y | x (unchecked)
 * Method: y0 z0 = x0 mod B = 2^BITS_IN_LONG ==> z0 = 1/y0 mod B.
 *    Set x := (x - z0 y) / B, updating only relevant words, and repeat */
GEN
diviiexact(GEN x, GEN y)
{
  long lx, ly, lz, vy, i, ii, sx = signe(x), sy = signe(y);
  pari_sp av;
  ulong y0inv,q;
  GEN z;

  if (!sy) pari_err(gdiver);
  if (!sx) return gen_0;
  lx = lgefint(x);
  if (lx == 3) {
    q = (ulong)x[2] / (ulong)y[2];
    return (sx+sy) ? utoipos(q): utoineg(q);
  }
  vy = vali(y); av = avma;
  (void)new_chunk(lx); /* enough room for z */
  if (vy)
  { /* make y odd */
    y = shifti(y,-vy);
    x = shifti(x,-vy); lx = lgefint(x);
  }
  else x = icopy(x); /* necessary because we destroy x */
  avma = av; /* will erase our x,y when exiting */
  /* now y is odd */
  ly = lgefint(y);
  if (ly == 3)
  {
    x = diviuexact_i(x,(ulong)y[2]); /* x != 0 */
    setsigne(x, (sx+sy)? 1: -1); return x;
  }
  y0inv = invrev(y[ly-1]);
  i=2; while (i<ly && y[i]==x[i]) i++;
  lz = (i==ly || (ulong)y[i] < (ulong)x[i]) ? lx-ly+3 : lx-ly+2;
  z = new_chunk(lz);

  y += ly - 1; /* now y[-i] = i-th word of y */
  for (ii=lx-1,i=lz-1; i>=2; i--,ii--)
  {
    long limj;
    LOCAL_HIREMAINDER;
    LOCAL_OVERFLOW;

    z[i] = q = y0inv*((ulong)x[ii]); /* i-th quotient */
    if (!q) continue;

    /* x := x - q * y */
    (void)mulll(q,y[0]); limj = max(lx - lz, ii+3-ly);
    { /* update neither lowest word (could set it to 0) nor highest ones */
      register GEN x0 = x + (ii - 1), y0 = y - 1, xlim = x + limj;
      for (; x0 >= xlim; x0--, y0--)
      {
        *x0 = subll(*x0, addmul(q,*y0));
        hiremainder += overflow;
      }
      if (hiremainder && limj != lx - lz)
      {
        if ((ulong)*x0 < hiremainder)
        {
          *x0 -= hiremainder;
          do (*--x0)--; while ((ulong)*x0 == MAXULONG);
        }
        else
          *x0 -= hiremainder;
      }
    }
  }
  i=2; while(!z[i]) i++;
  z += i-2; lz -= (i-2);
  z[0] = evaltyp(t_INT)|evallg(lz);
  z[1] = evalsigne((sx+sy)? 1: -1) | evallg(lz);
  avma = (pari_sp)z; return z;
}


/********************************************************************/
/**                                                                **/
/**               INTEGER MULTIPLICATION (KARATSUBA)               **/
/**                                                                **/
/********************************************************************/
/* nx >= ny = num. of digits of x, y (not GEN, see mulii) */
INLINE GEN
muliispec_basecase(GEN x, GEN y, long nx, long ny)
{
  GEN z2e,z2d,yd,xd,ye,zd;
  long p1,lz;
  LOCAL_HIREMAINDER;

  if (!ny) return gen_0;
  zd = (GEN)avma; lz = nx+ny+2;
  (void)new_chunk(lz);
  xd = x + nx;
  yd = y + ny;
  ye = yd; p1 = *--xd;

  *--zd = mulll(p1, *--yd); z2e = zd;
  while (yd > y) *--zd = addmul(p1, *--yd);
  *--zd = hiremainder;

  while (xd > x)
  {
    LOCAL_OVERFLOW;
    yd = ye; p1 = *--xd;

    z2d = --z2e;
    *z2d = addll(mulll(p1, *--yd), *z2d); z2d--;
    while (yd > y)
    {
      hiremainder += overflow;
      *z2d = addll(addmul(p1, *--yd), *z2d); z2d--;
    }
    *--zd = hiremainder + overflow;
  }
  if (*zd == 0) { zd++; lz--; } /* normalize */
  *--zd = evalsigne(1) | evallgefint(lz);
  *--zd = evaltyp(t_INT) | evallg(lz);
  avma=(pari_sp)zd; return zd;
}

INLINE GEN
sqrispec_basecase(GEN x, long nx)
{
  GEN z2e,z2d,yd,xd,zd,x0,z0;
  long p1,lz;
  LOCAL_HIREMAINDER;
  LOCAL_OVERFLOW;

  if (!nx) return gen_0;
  zd = (GEN)avma; lz = (nx+1) << 1;
  z0 = new_chunk(lz);
  if (nx == 1)
  {
    *--zd = mulll(*x, *x);
    *--zd = hiremainder; goto END;
  }
  xd = x + nx;

  /* compute double products --> zd */
  p1 = *--xd; yd = xd; --zd;
  *--zd = mulll(p1, *--yd); z2e = zd;
  while (yd > x) *--zd = addmul(p1, *--yd);
  *--zd = hiremainder;

  x0 = x+1;
  while (xd > x0)
  {
    LOCAL_OVERFLOW;
    p1 = *--xd; yd = xd;

    z2e -= 2; z2d = z2e;
    *z2d = addll(mulll(p1, *--yd), *z2d); z2d--;
    while (yd > x)
    {
      hiremainder += overflow;
      *z2d = addll(addmul(p1, *--yd), *z2d); z2d--;
    }
    *--zd = hiremainder + overflow;
  }
  /* multiply zd by 2 (put result in zd - 1) */
  zd[-1] = ((*zd & HIGHBIT) != 0);
  shift_left(zd, zd, 0, (nx<<1)-3, 0, 1);

  /* add the squares */
  xd = x + nx; zd = z0 + lz;
  p1 = *--xd;
  zd--; *zd = mulll(p1,p1);
  zd--; *zd = addll(hiremainder, *zd);
  while (xd > x)
  {
    p1 = *--xd;
    zd--; *zd = addll(mulll(p1,p1)+ overflow, *zd);
    zd--; *zd = addll(hiremainder + overflow, *zd);
  }

END:
  if (*zd == 0) { zd++; lz--; } /* normalize */
  *--zd = evalsigne(1) | evallgefint(lz);
  *--zd = evaltyp(t_INT) | evallg(lz);
  avma=(pari_sp)zd; return zd;
}

/* return (x shifted left d words) + y. Assume d > 0, x > 0 and y >= 0 */
static GEN
addshiftw(GEN x, GEN y, long d)
{
  GEN z,z0,y0,yd, zd = (GEN)avma;
  long a,lz,ly = lgefint(y);

  z0 = new_chunk(d);
  a = ly-2; yd = y+ly;
  if (a >= d)
  {
    y0 = yd-d; while (yd > y0) *--zd = *--yd; /* copy last d words of y */
    a -= d;
    if (a)
      z = addiispec(x+2, y+2, lgefint(x)-2, a);
    else
      z = icopy(x);
  }
  else
  {
    y0 = yd-a; while (yd > y0) *--zd = *--yd; /* copy last a words of y */
    while (zd >= z0) *--zd = 0;    /* complete with 0s */
    z = icopy(x);
  }
  lz = lgefint(z)+d;
  z[1] = evalsigne(1) | evallgefint(lz);
  z[0] = evaltyp(t_INT) | evallg(lz); return z;
}

/* Fast product (Karatsuba) of integers. a and b are "special" GENs
 * c,c0,c1,c2 are genuine GENs.
 */
GEN
muliispec(GEN a, GEN b, long na, long nb)
{
  GEN a0,c,c0;
  long n0, n0a, i;
  pari_sp av;

  if (na < nb) swapspec(a,b, na,nb);
  if (nb == 1) return muluispec((ulong)*b, a, na);
  if (nb == 0) return gen_0;
  if (nb < KARATSUBA_MULI_LIMIT) return muliispec_basecase(a,b,na,nb);
  i=(na>>1); n0=na-i; na=i;
  av=avma; a0=a+na; n0a=n0;
  while (!*a0 && n0a) { a0++; n0a--; }

  if (n0a && nb > n0)
  { /* nb <= na <= n0 */
    GEN b0,c1,c2;
    long n0b;

    nb -= n0;
    c = muliispec(a,b,na,nb);
    b0 = b+nb; n0b = n0;
    while (!*b0 && n0b) { b0++; n0b--; }
    if (n0b)
    {
      c0 = muliispec(a0,b0, n0a,n0b);

      c2 = addiispec(a0,a, n0a,na);
      c1 = addiispec(b0,b, n0b,nb);
      c1 = muliispec(c1+2,c2+2, lgefint(c1)-2,lgefint(c2)-2);
      c2 = addiispec(c0+2, c+2, lgefint(c0)-2,lgefint(c) -2);

      c1 = subiispec(c1+2,c2+2, lgefint(c1)-2,lgefint(c2)-2);
    }
    else
    {
      c0 = gen_0;
      c1 = muliispec(a0,b, n0a,nb);
    }
    c = addshiftw(c,c1, n0);
  }
  else
  {
    c = muliispec(a,b,na,nb);
    c0 = muliispec(a0,b,n0a,nb);
  }
  return gerepileuptoint(av, addshiftw(c,c0, n0));
}

/* x % (2^n), assuming x, n >= 0 */
GEN
resmod2n(GEN x, long n)
{
  long hi,l,k,lx,ly;
  GEN z, xd, zd;

  if (!signe(x) || !n) return gen_0;

  l = n & (BITS_IN_LONG-1);    /* n % BITS_IN_LONG */
  k = n >> TWOPOTBITS_IN_LONG; /* n / BITS_IN_LONG */
  lx = lgefint(x);
  if (lx < k+3) return icopy(x);

  xd = x + (lx-k-1);
  /* x = |_|...|#|1|...|k| : copy the last l bits of # and the last k words
   *            ^--- initial xd  */
  hi = ((ulong)*xd) & ((1UL<<l)-1); /* last l bits of # = top bits of result */
  if (!hi)
  { /* strip leading zeroes from result */
    xd++; while (k && !*xd) { k--; xd++; }
    if (!k) return gen_0;
    ly = k+2; xd--;
  }
  else
    ly = k+3;

  zd = z = cgeti(ly);
  *++zd = evalsigne(1) | evallgefint(ly);
  if (hi) *++zd = hi;
  for ( ;k; k--) *++zd = *++xd;
  return z;
}

GEN
sqrispec(GEN a, long na)
{
  GEN a0,c;
  long n0, n0a, i;
  pari_sp av;

  if (na < KARATSUBA_SQRI_LIMIT) return sqrispec_basecase(a,na);
  i=(na>>1); n0=na-i; na=i;
  av=avma; a0=a+na; n0a=n0;
  while (!*a0 && n0a) { a0++; n0a--; }
  c = sqrispec(a,na);
  if (n0a)
  {
    GEN t, c1, c0 = sqrispec(a0,n0a);
#if 0
    c1 = shifti(muliispec(a0,a, n0a,na),1);
#else /* faster */
    t = addiispec(a0,a,n0a,na);
    t = sqrispec(t+2,lgefint(t)-2);
    c1= addiispec(c0+2,c+2, lgefint(c0)-2, lgefint(c)-2);
    c1= subiispec(t+2, c1+2, lgefint(t)-2, lgefint(c1)-2);
#endif
    c = addshiftw(c,c1, n0);
    c = addshiftw(c,c0, n0);
  }
  else
    c = addshiftw(c,gen_0,n0<<1);
  return gerepileuptoint(av, c);
}

/********************************************************************/
/**                                                                **/
/**                    KARATSUBA SQUARE ROOT                       **/
/**      adapted from Paul Zimmermann's implementation of          **/
/**      his algorithm in GMP (mpn_sqrtrem)                        **/
/**                                                                **/
/********************************************************************/

/* Square roots table */
static const unsigned char approx_tab[192] = {
  128,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,
  143,144,144,145,146,147,148,149,150,150,151,152,153,154,155,155,
  156,157,158,159,160,160,161,162,163,163,164,165,166,167,167,168,
  169,170,170,171,172,173,173,174,175,176,176,177,178,178,179,180,
  181,181,182,183,183,184,185,185,186,187,187,188,189,189,190,191,
  192,192,193,193,194,195,195,196,197,197,198,199,199,200,201,201,
  202,203,203,204,204,205,206,206,207,208,208,209,209,210,211,211,
  212,212,213,214,214,215,215,216,217,217,218,218,219,219,220,221,
  221,222,222,223,224,224,225,225,226,226,227,227,228,229,229,230,
  230,231,231,232,232,233,234,234,235,235,236,236,237,237,238,238,
  239,240,240,241,241,242,242,243,243,244,244,245,245,246,246,247,
  247,248,248,249,249,250,250,251,251,252,252,253,253,254,254,255
};

/* N[0], assume N[0] >= 2^(BIL-2).
 * Return r,s such that s^2 + r = N, 0 <= r <= 2s */
static void
p_sqrtu1(ulong *N, ulong *ps, ulong *pr)
{
  ulong prec, r, s, q, u, n0 = N[0];

  q = n0 >> (BITS_IN_LONG - 8);
  /* 2^6 = 64 <= q < 256 = 2^8 */
  s = approx_tab[q - 64];				/* 128 <= s < 255 */
  r = (n0 >> (BITS_IN_LONG - 16)) - s * s;		/* r <= 2*s */
  if (r > (s << 1)) { r -= (s << 1) | 1; s++; }

  /* 8-bit approximation from the high 8-bits of N[0] */
  prec = 8;
  n0 <<= 2 * prec;
  while (2 * prec < BITS_IN_LONG)
  { /* invariant: s has prec bits, and r <= 2*s */
    r = (r << prec) + (n0 >> (BITS_IN_LONG - prec));
    n0 <<= prec;
    u = 2 * s;
    q = r / u; u = r - q * u;
    s = (s << prec) + q;
    u = (u << prec) + (n0 >> (BITS_IN_LONG - prec));
    q = q * q;
    r = u - q;
    if (u < q) { s--; r += (s << 1) | 1; }
    n0 <<= prec;
    prec = 2 * prec;
  }
  *ps = s;
  *pr = r;
}

/* N[0..1], assume N[0] >= 2^(BIL-2).
 * Return 1 if remainder overflows, 0 otherwise */
static int
p_sqrtu2(ulong *N, ulong *ps, ulong *pr)
{
  ulong cc, qhl, r, s, q, u, n1 = N[1];
  LOCAL_OVERFLOW;

  p_sqrtu1(N, &s, &r); /* r <= 2s */
  qhl = 0; while (r >= s) { qhl++; r -= s; }
  /* now r < s < 2^(BIL/2) */
  r = (r << BITS_IN_HALFULONG) | (n1 >> BITS_IN_HALFULONG);
  u = s << 1;
  q = r / u; u = r - q * u;
  q += (qhl & 1) << (BITS_IN_HALFULONG - 1);
  qhl >>= 1;
  /* (initial r)<<(BIL/2) + n1>>(BIL/2) = (qhl<<(BIL/2) + q) * 2s + u */
  s = ((s + qhl) << BITS_IN_HALFULONG) + q;
  cc = u >> BITS_IN_HALFULONG;
  r = (u << BITS_IN_HALFULONG) | (n1 & LOWMASK);
  r = subll(r, q * q);
  cc -= overflow + qhl;
  /* now subtract 2*q*2^(BIL/2) + 2^BIL if qhl is set */
  if ((long)cc < 0)
  {
    if (s) {
      r = addll(r, s);
      cc += overflow;
      s--;
    } else {
      cc++;
      s = ~0UL;
    }
    r = addll(r, s);
    cc += overflow;
  }
  *ps = s;
  *pr = r; return cc;
}

static void
xmpn_zero(GEN x, long n)
{
  while (--n >= 0) x[n]=0;
}
static void
xmpn_copy(GEN z, GEN x, long n)
{
  long k = n;
  while (--k >= 0) z[k] = x[k];
}
static GEN
cat1u(ulong d)
{
  GEN R = cgeti(4);
  R[1] = evalsigne(1)|evallgefint(4);
  R[2] = 1;
  R[3] = d; return R;
}
/* a[0..la-1] * 2^(lb BIL) | b[0..lb-1] */
static GEN
catii(GEN a, long la, GEN b, long lb)
{
  long l = la + lb + 2;
  GEN z = cgeti(l);
  z[1] = evalsigne(1) | evallgefint(l);
  xmpn_copy(z + 2, a, la);
  xmpn_copy(z + 2 + la, b, lb);
  return int_normalize(z, 0);
}

/* sqrt n[0..1], assume n normalized */
static GEN
sqrtispec2(GEN n, GEN *pr)
{
  ulong s, r;
  int hi = p_sqrtu2((ulong*)n, &s, &r);
  GEN S = utoi(s);
  *pr = hi? cat1u(r): utoi(r);
  return S;
}

/* sqrt n[0], _dont_ assume n normalized */
static GEN
sqrtispec1_sh(GEN n, GEN *pr)
{
  GEN S;
  ulong r, s, u0 = (ulong)n[0];
  int sh = bfffo(u0) & ~1UL;
  if (sh) u0 <<= sh;
  p_sqrtu1(&u0, &s, &r);
  /* s^2 + r = u0, s < 2^(BIL/2). Rescale back:
   * 2^(2k) n = S^2 + R
   * so 2^(2k) n = (S - s0)^2 + (2*S*s0 - s0^2 + R), s0 = S mod 2^k. */
  if (sh) {
    int k = sh >> 1;
    ulong s0 = s & ((1<<k) - 1);
    r += s * (s0<<1);
    s >>= k;
    r >>= sh;
  }
  S = utoi(s);
  if (pr) *pr = utoi(r);
  return S;
}

/* sqrt n[0..1], _dont_ assume n normalized */
static GEN
sqrtispec2_sh(GEN n, GEN *pr)
{
  GEN S;
  ulong U[2], r, s, u0 = (ulong)n[0], u1 = (ulong)n[1];
  int hi, sh = bfffo(u0) & ~1UL;
  if (sh) {
    u0 = (u0 << sh) | (u1 >> (BITS_IN_LONG-sh));
    u1 <<= sh;
  }
  U[0] = u0;
  U[1] = u1; hi = p_sqrtu2(U, &s, &r);
  /* s^2 + R = u0|u1. Rescale back:
   * 2^(2k) n = S^2 + R
   * so 2^(2k) n = (S - s0)^2 + (2*S*s0 - s0^2 + R), s0 = S mod 2^k. */
  if (sh) {
    int k = sh >> 1;
    ulong s0 = s & ((1<<k) - 1);
    LOCAL_HIREMAINDER;
    LOCAL_OVERFLOW;
    r = addll(r, mulll(s, (s0<<1)));
    if (overflow) hiremainder++;
    hiremainder += hi; /* + 0 or 1 */
    s >>= k;
    r = (r>>sh) | (hiremainder << (BITS_IN_LONG-sh));
    hi = (hiremainder & (1<<sh));
  }
  S = utoi(s);
  if (pr) *pr = hi? cat1u(r): utoi(r);
  return S;
}

/* Let N = N[0..2n-1]. Return S (and set R) s.t S^2 + R = N, 0 <= R <= 2S
 * Assume N normalized */
static GEN
sqrtispec(GEN N, long n, GEN *r)
{
  GEN S, R, q, z, u;
  long l, h;

  if (n == 1) return sqrtispec2(N, r);
  l = n >> 1;
  h = n - l; /* N = a3(h) | a2(h) | a1(l) | a0(l words) */
  S = sqrtispec(N, h, &R); /* S^2 + R = a3|a2 */

  z = catii(R+2, lgefint(R)-2, N + 2*h, l); /* = R | a1(l) */
  q = dvmdii(z, shifti(S,1), &u);
  z = catii(u+2, lgefint(u)-2, N + n + h, l); /* = u | a0(l) */

  S = addshiftw(S, q, l);
  R = subii(z, sqri(q));
  if (signe(R) < 0)
  {
    GEN S2 = shifti(S,1);
    R = addis(subiispec(S2+2, R+2, lgefint(S2)-2,lgefint(R)-2), -1);
    S = addis(S, -1);
  }
  *r = R; return S;
}

/* Return S (and set R) s.t S^2 + R = N, 0 <= R <= 2S.
 * As for dvmdii, R is last on stack and guaranteed to be gen_0 in case the
 * remainder is 0. R = NULL is allowed. */
GEN
sqrtremi(GEN N, GEN *r)
{
  pari_sp av;
  GEN S, R, n = N+2;
  long k, l2, ln = lgefint(N) - 2;
  int sh;

  if (ln <= 2)
  {
    if (ln == 2) return sqrtispec2_sh(n, r);
    if (ln == 1) return sqrtispec1_sh(n, r);
    if (r) *r = gen_0;
    return gen_0;
  }
  av = avma;
  sh = bfffo(n[0]) >> 1;
  l2 = (ln + 1) >> 1;
  if (sh || (ln & 1)) { /* normalize n, so that n[0] >= 2^BIL / 4 */
    GEN s0, t = new_chunk(ln + 1);
    t[ln] = 0;
    if (sh)
    { shift_left(t, n, 0,ln-1, 0, (sh << 1)); }
    else
      xmpn_copy(t, n, ln);
    S = sqrtispec(t, l2, &R); /* t normalized, 2 * l2 words */
    /* Rescale back:
     * 2^(2k) n = S^2 + R, k = sh + (ln & 1)*BIL/2
     * so 2^(2k) n = (S - s0)^2 + (2*S*s0 - s0^2 + R), s0 = S mod 2^k. */
    k = sh + (ln & 1) * (BITS_IN_LONG/2);
    s0 = resmod2n(S, k);
    R = addii(shifti(R,-1), mulii(s0, S));
    R = shifti(R, 1 - (k<<1));
    S = shifti(S, -k);
  }
  else
    S = sqrtispec(n, l2, &R);

  if (!r) { avma = (pari_sp)S; return gerepileuptoint(av, S); }
  gerepileall(av, 2, &S, &R); *r = R; return S;
}

/* compute sqrt(|a|), assuming a != 0 */

#if 1
GEN
sqrtr_abs(GEN x)
{
  long l = lg(x) - 2, e = expo(x), er = e>>1;
  GEN b, c, res = cgetr(2 + l);
  res[1] = evalsigne(1) | evalexpo(er);
  if (e&1) {
    b = new_chunk(l << 1);
    xmpn_copy(b, x+2, l);
    xmpn_zero(b + l,l);
    b = sqrtispec(b, l, &c);
    xmpn_copy(res+2, b+2, l);
    if (cmpii(c, b) > 0) roundr_up_ip(res, l+2);
  } else {
    ulong u;
    b = new_chunk(2 + (l << 1));
    shift_left(b+1, x+2, 0,l-1, 0, BITS_IN_LONG-1);
    b[0] = ((ulong)x[2])>>1;
    xmpn_zero(b + l+1,l+1);
    b = sqrtispec(b, l+1, &c);
    xmpn_copy(res+2, b+2, l);
    u = (ulong)b[l+2];
    if ( u&HIGHBIT || (u == ~HIGHBIT && cmpii(c,b) > 0))
      roundr_up_ip(res, l+2);
  }
  avma = (pari_sp)res; return res;
}

#else /* use t_REAL: currently much slower (quadratic division) */

#ifdef LONG_IS_64BIT
/* 64 bits of b = sqrt(a[0] * 2^64 + a[1])  [ up to 1ulp ] */
static ulong
sqrtu2(ulong *a)
{
  ulong c, b = dblmantissa( sqrt((double)a[0]) );
  LOCAL_HIREMAINDER;
  LOCAL_OVERFLOW;

  /* > 32 correct bits, 1 Newton iteration to reach 64 */
  if (b <= a[0]) return HIGHBIT | (a[0] >> 1);
  hiremainder = a[0]; c = divll(a[1], b);
  return (addll(c, b) >> 1) | HIGHBIT;
}
/* 64 bits of sqrt(a[0] * 2^63) */
static ulong
sqrtu2_1(ulong *a)
{
  ulong t[2];
  t[0] = (a[0] >> 1);
  t[1] = (a[0] << (BITS_IN_LONG-1)) | (a[1] >> 1);
  return sqrtu2(t);
}
#else
/* 32 bits of sqrt(a[0] * 2^32) */
static ulong
sqrtu2(ulong *a)   { return dblmantissa( sqrt((double)a[0]) ); }
/* 32 bits of sqrt(a[0] * 2^31) */
static ulong
sqrtu2_1(ulong *a) { return dblmantissa( sqrt(2. * a[0]) ); }
#endif

GEN
sqrtr_abs(GEN x)
{
  long l1, i, l = lg(x), ex = expo(x);
  GEN a, t, y = cgetr(l);
  pari_sp av, av0 = avma;

  a = cgetr(l+1); affrr(x,a);
  t = cgetr(l+1);
  if (ex & 1) { /* odd exponent */
    a[1] = evalsigne(1) | evalexpo(1);
    t[2] = (long)sqrtu2((ulong*)a + 2);
  } else { /* even exponent */
    a[1] = evalsigne(1) | evalexpo(0);
    t[2] = (long)sqrtu2_1((ulong*)a + 2);
  }
  t[1] = evalsigne(1) | evalexpo(0);
  for (i = 3; i <= l; i++) t[i] = 0;

  /* |x| = 2^(ex/2) a, t ~ sqrt(a) */
  l--; l1 = 1; av = avma;
  while (l1 < l) { /* let t := (t + a/t)/2 */
    l1 <<= 1; if (l1 > l) l1 = l;
    setlg(a, l1 + 2);
    setlg(t, l1 + 2);
    affrr(addrr(t, divrr(a,t)), t); setexpo(t, expo(t)-1);
    avma = av;
  }
  affrr(t,y); setexpo(y, expo(y) + (ex>>1));
  avma = av0; return y;
}

#endif
