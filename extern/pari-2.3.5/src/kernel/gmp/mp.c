#line 2 "../src/kernel/gmp/mp.c"
/* $Id: mp.c 10295 2008-06-10 16:03:56Z kb $

Copyright (C) 2002-2003  The PARI group.

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
/**		                GMP KERNEL                     	      **/
/** BA2002Sep24                                                       **/
/***********************************************************************/
/* GMP t_INT as just like normal t_INT, just the mantissa is the other way
 * round
 *
 *   `How would you like to live in Looking-glass House, Kitty?  I
 *   wonder if they'd give you milk in there?  Perhaps Looking-glass
 *   milk isn't good to drink--But oh, Kitty! now we come to the
 *   passage.  You can just see a little PEEP of the passage in
 *   Looking-glass House, if you leave the door of our drawing-room
 *   wide open:  and it's very like our passage as far as you can see,
 *   only you know it may be quite different on beyond.  Oh, Kitty!
 *   how nice it would be if we could only get through into Looking-
 *   glass House!  I'm sure it's got, oh! such beautiful things in it!
 *                                                                             
 *  Through the Looking Glass,  Lewis Carrol
 *  
 *  (pityful attempt to beat GN code/comments rate)
 *  */

#include <gmp.h>
#include "pari.h"
#include "paripriv.h"
#include "../src/kernel/none/tune-gen.h"

/*We need PARI invmod renamed to invmod_pari*/
#define INVMOD_PARI

static void *wrap_gprealloc(void *ptr, size_t old_size, size_t new_size) {
  (void)old_size; return (void *) gprealloc(ptr,new_size);
}

int pari_kernel_init(void)
{
  /* Use gpmalloc instead of malloc */
  mp_set_memory_functions((void *(*)(size_t)) gpmalloc, wrap_gprealloc, NULL);
  return 0;
}

#define LIMBS(x)  ((mp_limb_t *)((x)+2))
#define NLIMBS(x) (lgefint(x)-2)
/*This one is for t_REAL to emphasize there are not t_INT*/
#define RLIMBS(x)  ((mp_limb_t *)((x)+2))
#define RNLIMBS(x) (lg(x)-2)

INLINE void
xmpn_copy(mp_limb_t *x, mp_limb_t *y, long n)
{
  while (--n >= 0) x[n]=y[n];
}

INLINE void
xmpn_mirror(mp_limb_t *x, long n)
{
  long i;
  for(i=0;i<(n>>1);i++)
  {
    ulong m=x[i];
    x[i]=x[n-1-i];
    x[n-1-i]=m;
  }
}

INLINE void
xmpn_mirrorcopy(mp_limb_t *z, mp_limb_t *x, long n)
{
  long i;
  for(i=0;i<n;i++)
    z[i]=x[n-1-i];
}

INLINE void
xmpn_zero(mp_ptr x, mp_size_t n)
{
  while (--n >= 0) x[n]=0;
}

INLINE GEN
icopy_ef(GEN x, long l)
{
  register long lx = lgefint(x);
  const GEN y = cgeti(l);

  while (--lx > 0) y[lx]=x[lx];
  return y;
}

/* NOTE: arguments of "spec" routines (muliispec, addiispec, etc.) aren't
 * GENs but pairs (long *a, long na) representing a list of digits (in basis
 * BITS_IN_LONG) : a[0], ..., a[na-1]. [ In ordre to facilitate splitting: no
 * need to reintroduce codewords ]
 * Use speci(a,na) to visualize the corresponding GEN.
 */

/***********************************************************************/
/**								      **/
/**		         ADDITION / SUBTRACTION          	      **/
/**                                                                   **/
/***********************************************************************/

GEN
setloop(GEN a)
{
  GEN z0 = (GEN)avma; (void)cgetg(lgefint(a) + 3, t_VECSMALL);
  return icopy_av(a, z0 - 2); /* two cells of extra space after a */
}

/* we had a = setloop(?), then some incloops. Reset a to b */
GEN
resetloop(GEN a, GEN b) {
  a[0] = evaltyp(t_INT) | evallg(lgefint(b));
  affii(b, a); return a;
}

/* assume a > 0, initialized by setloop. Do a++ */
static GEN
incpos(GEN a)
{
  long i, l = lgefint(a);
  for (i=2; i<l; i++)
    if (++a[i]) return a;
  a[l] = 1; l++;
  a[0]=evaltyp(t_INT) | _evallg(l);
  a[1]=evalsigne(1) | evallgefint(l);
  return a;
}

/* assume a < 0, initialized by setloop. Do a++ */
static GEN
incneg(GEN a)
{
  long i, l = lgefint(a);
  if (a[2]--)
  {
    if (l == 3 && !a[2])
    {
      a[0] = evaltyp(t_INT) | _evallg(2);
      a[1] = evalsigne(0) | evallgefint(2);
    }
    return a;
  }
  for (i=3; i<l; i++)
    if (a[i]--) break;
  l -= i - 2;
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
    case 0:
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
  GEN  zd;
  long lz;

  lz = nx+3; zd = cgeti(lz);
  if (mpn_add_1(LIMBS(zd),(mp_limb_t *)x,nx,s))
    zd[lz-1]=1;
  else
    lz--;
  zd[1] = evalsigne(1) | evallgefint(lz);
  return zd;
}

INLINE GEN
addiispec(GEN x, GEN y, long nx, long ny)
{
  GEN zd;
  long lz;

  if (nx < ny) swapspec(x,y, nx,ny);
  if (ny == 1) return addsispec(*y,x,nx);
  lz = nx+3; zd = cgeti(lz);

  if (mpn_add(LIMBS(zd),(mp_limb_t *)x,nx,(mp_limb_t *)y,ny))
    zd[lz-1]=1;
  else
    lz--;

  zd[1] = evalsigne(1) | evallgefint(lz);
  return zd;
}

/* assume x >= y */
INLINE GEN
subisspec(GEN x, long s, long nx)
{
  GEN zd;
  long lz;
  lz = nx + 2; zd = cgeti(lz);
  
  mpn_sub_1 (LIMBS(zd), (mp_limb_t *)x, nx, s);
  if (! zd[lz - 1]) { --lz; }

  zd[1] = evalsigne(1) | evallgefint(lz);
  return zd;
}

/* assume x > y */
INLINE GEN
subiispec(GEN x, GEN y, long nx, long ny)
{
  GEN zd;
  long lz;
  if (ny==1) return subisspec(x,*y,nx);
  lz = nx+2; zd = cgeti(lz);
  
  mpn_sub (LIMBS(zd), (mp_limb_t *)x, nx, (mp_limb_t *)y, ny);
  while (lz >= 3 && zd[lz - 1] == 0) { lz--; }
  
  zd[1] = evalsigne(1) | evallgefint(lz);
  return zd;
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
  lx = lgefint(x); sh = bfffo(*int_MSW(x));
  y[1] = evalsigne(s) | evalexpo(bit_accuracy(lx)-sh-1);
  if (sh) {
    if (lx <= ly)
    {
      for (i=lx; i<ly; i++) y[i]=0;
      mpn_lshift(LIMBS(y),LIMBS(x),lx-2,sh);
      xmpn_mirror(LIMBS(y),lx-2);
      return;
    }
    mpn_lshift(LIMBS(y),LIMBS(x)+lx-ly,ly-2,sh);
    y[2]|=((ulong) x[lx-ly+1])>>(BITS_IN_LONG-sh);
    xmpn_mirror(LIMBS(y),ly-2);
    /* lx > ly: round properly */
    if ((x[lx-ly+1]<<sh) & HIGHBIT) roundr_up_ip(y, ly);
  }
  else {
    GEN xd=int_MSW(x);
    if (lx <= ly)
    {
      for (i=2; i<lx; i++,xd=int_precW(xd)) y[i]=*xd;
      for (   ; i<ly; i++) y[i]=0;
      return;
    }
    for (i=2; i<ly; i++,xd=int_precW(xd)) y[i]=*xd;
    /* lx > ly: round properly */
    if (x[lx-ly+1] & HIGHBIT) roundr_up_ip(y, ly);
  }
}

GEN
shifti(GEN x, long n)
{
  const long s=signe(x);
  long lz,lx,m;
  GEN z;

  if (!s) return gen_0;
  if (!n) return icopy(x);
  lx = lgefint(x);
  if (n > 0)
  {
    long d = n>>TWOPOTBITS_IN_LONG;
    long i;

    m = n & (BITS_IN_LONG-1);

    lz = lx + d + (m!=0);  
    z = cgeti(lz); 
    for (i=0; i<d; i++) LIMBS(z)[i] = 0;

    if (!m) xmpn_copy(LIMBS(z)+d, LIMBS(x), NLIMBS(x)); 
    else
    {
      ulong carry = mpn_lshift(LIMBS(z)+d, LIMBS(x), NLIMBS(x), m); 
      if (carry) z[lz - 1] = carry; 
      else lz--; 
    }
  }
  else
  {
    long d = (-n)>>TWOPOTBITS_IN_LONG;

    n = -n;
    lz = lx - d;
    if (lz<3) return gen_0;
    z = cgeti(lz);
    m = n & (BITS_IN_LONG-1);

    if (!m) xmpn_copy(LIMBS(z), LIMBS(x) + d, NLIMBS(x) - d);
    else
    {
      mpn_rshift(LIMBS(z), LIMBS(x) + d, NLIMBS(x) - d, m); 
      if (z[lz - 1] == 0)
      {
        if (lz == 3) { avma = (pari_sp)(z+3); return gen_0; }
        lz--; 
      }
    }
  }
  z[1] = evalsigne(s)|evallgefint(lz);
  return z;
}

GEN
ishiftr_lg(GEN x, long lx, long n)
{
  long ly, i, m, s = signe(x);
  GEN y;
  if (!s) return gen_0;
  if (!n) 
  {
    y = cgeti(lx);
    y[1] = evalsigne(s) | evallgefint(lx);
    xmpn_mirrorcopy(LIMBS(y),RLIMBS(x),lx-2);
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
  xmpn_mirror(LIMBS(y),ly-2);
  y[1] = evalsigne(s)|evallgefint(ly);
  y[0] = evaltyp(t_INT)|evallg(ly); return y;
}

GEN
truncr(GEN x)
{
  long s, e, d, m, i;
  GEN y;
  if ((s=signe(x)) == 0 || (e=expo(x)) < 0) return gen_0;
  d = (e>>TWOPOTBITS_IN_LONG) + 3;
  m = e & (BITS_IN_LONG-1);
  if (d > lg(x)) pari_err(precer, "truncr (precision loss in truncation)");

  y=cgeti(d); y[1] = evalsigne(s) | evallgefint(d);
  if (++m == BITS_IN_LONG)
    for (i=2; i<d; i++) y[d-i+1]=x[i];
  else
  {
    GEN z=cgeti(d);
    for (i=2; i<d; i++) z[d-i+1]=x[i];
    mpn_rshift(LIMBS(y),LIMBS(z),d-2,BITS_IN_LONG-m);
    avma=(pari_sp)y;
  }
  return y;
}

/* integral part */
GEN
floorr(GEN x)
{
  long e, d, m, i, lx;
  GEN y;
  if (signe(x) >= 0) return truncr(x);
  if ((e=expo(x)) < 0) return gen_m1;
  d = (e>>TWOPOTBITS_IN_LONG) + 3;
  m = e & (BITS_IN_LONG-1);
  lx=lg(x); if (d>lx) pari_err(precer, "floorr (precision loss in truncation)");
  y = cgeti(d+1);
  if (++m == BITS_IN_LONG)
  {
    for (i=2; i<d; i++) y[d-i+1]=x[i];
    i=d; while (i<lx && !x[i]) i++;
    if (i==lx) goto END;
  }
  else
  {
    GEN z=cgeti(d);
    for (i=2; i<d; i++) z[d-i+1]=x[i];
    mpn_rshift(LIMBS(y),LIMBS(z),d-2,BITS_IN_LONG-m);
    if (x[d-1]<<m == 0)
    {
      i=d; while (i<lx && !x[i]) i++;
      if (i==lx) goto END;
    }
  }
  if (mpn_add_1(LIMBS(y),LIMBS(y),d-2,1))
    y[d++]=1; 
END:
  y[1] = evalsigne(-1) | evallgefint(d);
  return y;
}

INLINE int
absi_cmp_lg(GEN x, GEN y, long l)
{
  return mpn_cmp(LIMBS(x),LIMBS(y),l-2);
}

INLINE int
absi_equal_lg(GEN x, GEN y, long l)
{
  return !mpn_cmp(LIMBS(x),LIMBS(y),l-2);
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
    z[3]=hiremainder; z[2]=p1; return z;
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
    z[3]=hiremainder; z[2]=p1; return z;
  }
  z=cgeti(3); z[1] = evalsigne(1) | evallgefint(3);
  z[2]=p1; return z;
}

/* assume ny > 0 */
INLINE GEN
muluispec(ulong x, GEN y, long ny)
{
  long lz = ny+3;
  GEN z=cgeti(lz);
  ulong hi = mpn_mul_1 (LIMBS(z), (mp_limb_t *)y, ny, x);
  if (hi) { z[lz - 1] = hi; } else lz--;
  z[1] = evalsigne(1) | evallgefint(lz);
  return z;
}

/* a + b*|y| */
GEN
addumului(ulong a, ulong b, GEN y)
{
  GEN z;
  long i, lz;
  ulong hi;
  if (!signe(y)) return utoi(a);
  lz = lgefint(y)+1;
  z = cgeti(lz);
  z[2]=a;
  for(i=3;i<lz;i++) z[i]=0;
  hi=mpn_addmul_1(LIMBS(z), LIMBS(y), NLIMBS(y), b);
  if (hi) z[lz-1]=hi; else lz--;
  z[1] = evalsigne(1) | evallgefint(lz);
  avma=(pari_sp)z; return z;
}

GEN muliispec(GEN x, GEN y, long nx, long ny);

/* We must have nx>=ny. This lets garbage on the stack.
   This handle squares correctly (mpn_mul is optimized
   for squares).
*/

INLINE GEN
muliispec_mirror(GEN x, GEN y, long nx, long ny)
{
  GEN cx=new_chunk(nx),cy;
  GEN z;
  xmpn_mirrorcopy((mp_limb_t *)cx,(mp_limb_t *)x,nx);
  if (x==y) cy=cx; /*If nx<ny cy will be too short*/
  else
  {
    cy=new_chunk(ny);
    xmpn_mirrorcopy((mp_limb_t *)cy,(mp_limb_t *)y,ny);
  }
  z=muliispec(cx, cy, nx, ny);
  xmpn_mirror(LIMBS(z), NLIMBS(z));
  return z;
}

/*#define KARAMULR_VARIANT*/

/***********************************************************************/
/**								      **/
/**		          DIVISION                       	      **/
/**                                                                   **/
/***********************************************************************/

ulong
umodiu(GEN y, ulong x)
{
  long sy=signe(y);
  ulong hi;
  if (!x) pari_err(gdiver);
  if (!sy) return 0;
  hi = mpn_mod_1(LIMBS(y),NLIMBS(y),x);
  if (!hi) return 0;
  return (sy > 0)? hi: x - hi;
}

/* return |y| \/ x */
GEN
diviu_rem(GEN y, ulong x, ulong *rem)
{
  long ly;
  GEN z;

  if (!x) pari_err(gdiver);
  if (!signe(y)) { *rem = 0; return gen_0; }

  ly = lgefint(y);
  if (ly == 3 && (ulong)x > (ulong)y[2]) { *rem = (ulong)y[2]; return gen_0; }

  z = cgeti(ly); 
  *rem = mpn_divrem_1(LIMBS(z), 0, LIMBS(y), NLIMBS(y), x);
  if (z [ly - 1] == 0) ly--;
  z[1] = evallgefint(ly) | evalsigne(1);
  return z;
}

GEN
divis_rem(GEN y, long x, long *rem)
{
  long sy=signe(y),ly,s;
  GEN z;

  if (!x) pari_err(gdiver);
  if (!sy) { *rem = 0; return gen_0; }
  if (x<0) { s = -sy; x = -x; } else s = sy;

  ly = lgefint(y);
  if (ly == 3 && (ulong)x > (ulong)y[2]) { *rem = itos(y); return gen_0; }

  z = cgeti(ly); 
  *rem = mpn_divrem_1(LIMBS(z), 0, LIMBS(y), NLIMBS(y), x);
  if (sy<0) *rem = -  *rem;
  if (z[ly - 1] == 0) ly--;
  z[1] = evallgefint(ly) | evalsigne(s);
  return z;
}

GEN
divis(GEN y, long x)
{
  long sy=signe(y),ly,s;
  GEN z;

  if (!x) pari_err(gdiver);
  if (!sy) return gen_0;
  if (x<0) { s = -sy; x = -x; } else s=sy;

  ly = lgefint(y);
  if (ly == 3 && (ulong)x > (ulong)y[2]) return gen_0;

  z = cgeti(ly); 
  (void)mpn_divrem_1(LIMBS(z), 0, LIMBS(y), NLIMBS(y), x);
  if (z[ly - 1] == 0) ly--;
  z[1] = evallgefint(ly) | evalsigne(s);
  return z;
}

/* We keep llx bits of x and lly bits of y*/
static GEN
divrr_with_gmp(GEN x, GEN y)
{
  long lx=RNLIMBS(x),ly=RNLIMBS(y);
  long lw=min(lx,ly);
  long lly=min(lw+1,ly);
  GEN  w=cgetr(lw+2);
  long lu=lw+lly;
  long llx=min(lu,lx);
  mp_limb_t *u=(mp_limb_t *)new_chunk(lu);
  mp_limb_t *z=(mp_limb_t *)new_chunk(lly);
  mp_limb_t *q,*r;
  long e=expo(x)-expo(y);
  long sx=signe(x),sy=signe(y);
  xmpn_mirrorcopy(z,RLIMBS(y),lly);
  xmpn_mirrorcopy(u+lu-llx,RLIMBS(x),llx);
  xmpn_zero(u,lu-llx);
  q = (mp_limb_t *)new_chunk(lw+1);
  r = (mp_limb_t *)new_chunk(lly);

  mpn_tdiv_qr(q,r,0,u,lu,z,lly);
  
  /*Round up: This is not exactly correct we should test 2*r>z*/
  if ((ulong)r[lly-1] > ((ulong)z[lly-1]>>1))
    mpn_add_1(q,q,lw+1,1);
  
  xmpn_mirrorcopy(RLIMBS(w),q,lw);

  if (q[lw] == 0) e--;
  else if (q[lw] == 1) { shift_right(w,w, 2,lw+2, 1,1); }
  else { w[2] = HIGHBIT; e++; }
  if (sy < 0) sx = -sx;
  w[1] = evalsigne(sx) | evalexpo(e);
  avma=(pari_sp) w;
  return w;
}

/* We keep llx bits of x and lly bits of y*/
static GEN
divri_with_gmp(GEN x, GEN y)
{
  long llx=RNLIMBS(x),ly=NLIMBS(y);
  long lly=min(llx+1,ly);
  GEN  w=cgetr(llx+2);
  long lu=llx+lly, ld=ly-lly;
  mp_limb_t *u=(mp_limb_t *)new_chunk(lu);
  mp_limb_t *z=(mp_limb_t *)new_chunk(lly);
  mp_limb_t *q,*r;
  long sh=bfffo(y[ly+1]);
  long e=expo(x)-expi(y);
  long sx=signe(x),sy=signe(y);
  if (sh) mpn_lshift(z,LIMBS(y)+ld,lly,sh);
  else xmpn_copy(z,LIMBS(y)+ld,lly);
  xmpn_mirrorcopy(u+lu-llx,RLIMBS(x),llx);
  xmpn_zero(u,lu-llx);
  q = (mp_limb_t *)new_chunk(llx+1);
  r = (mp_limb_t *)new_chunk(lly);

  mpn_tdiv_qr(q,r,0,u,lu,z,lly);
  
  /*Round up: This is not exactly correct we should test 2*r>z*/
  if ((ulong)r[lly-1] > ((ulong)z[lly-1]>>1))
    mpn_add_1(q,q,llx+1,1);
  
  xmpn_mirrorcopy(RLIMBS(w),q,llx);

  if (q[llx] == 0) e--;
  else if (q[llx] == 1) { shift_right(w,w, 2,llx+2, 1,1); }
  else { w[2] = HIGHBIT; e++; }
  if (sy < 0) sx = -sx;
  w[1] = evalsigne(sx) | evalexpo(e);
  avma=(pari_sp) w;
  return w;
}

GEN
divri(GEN x, GEN y)
{
  long  s = signe(y);

  if (!s) pari_err(gdiver);
  if (!signe(x)) return real_0_bit(expo(x) - expi(y));
  if (!is_bigint(y)) return divrs(x, s>0? y[2]: -y[2]);
  return divri_with_gmp(x,y);
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

  if (ly>=DIVRR_GMP_LIMIT)
    return divrr_with_gmp(x,y);

  lr = min(lx,ly); r = new_chunk(lr);
  r1 = r-1;
  r1[1] = 0; for (i=2; i<lr; i++) r1[i]=x[i];
  r1[lr] = (lx>ly)? x[lr]: 0;
  y0 = y[2]; y1 = y[3];
  for (i=0; i<lr-1; i++)
  { /* r1 = r + (i-1) */
    ulong k,qp;
    LOCAL_HIREMAINDER;
    LOCAL_OVERFLOW;
    
    if ((ulong)r1[1] == y0)
    {
      qp = MAXULONG; k=addll(y0,r1[2]);
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
      while (!overflow && k4) { qp--; k3=subll(k3,y1); k4=subllx(k4,y0); }
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

/* Integer division x / y: such that sign(r) = sign(x)
 *   if z = ONLY_REM return remainder, otherwise return quotient
 *   if z != NULL set *z to remainder
 *   *z is the last object on stack (and thus can be disposed of with cgiv
 *   instead of gerepile)
 * If *z is zero, we put gen_0 here and no copy.
 * space needed: lx + ly
 */
GEN
dvmdii(GEN x, GEN y, GEN *z)
{
  long sx=signe(x),sy=signe(y);
  long lx, ly, lq;
  pari_sp av;
  GEN r,q;

  if (!sy) { if (z == ONLY_REM && !sx) return gen_0; pari_err(gdiver); }
  if (!sx)
  {
    if (!z || z == ONLY_REM) return gen_0;
    *z=gen_0; return gen_0;
  }
  lx=lgefint(x);
  ly=lgefint(y); lq=lx-ly;
  if (lq <= 0)
  {
    if (lq == 0)
    {
      long s=mpn_cmp(LIMBS(x),LIMBS(y),NLIMBS(x));
      if (s>0) goto DIVIDE;
      if (s==0) 
      {
        if (z == ONLY_REM) return gen_0;
        if (z) *z = gen_0;
        if (sx < 0) sy = -sy;
        return stoi(sy);
      }
    }
    if (z == ONLY_REM) return icopy(x);
    if (z) *z = icopy(x);
    return gen_0;
  }
DIVIDE: /* quotient is non-zero */
  av=avma; if (sx<0) sy = -sy;
  if (ly==3)
  {
    ulong lq = lx;
    ulong si;
    q  = cgeti(lq);
    si = mpn_divrem_1(LIMBS(q), 0, LIMBS(x), NLIMBS(x), y[2]);
    if (q[lq - 1] == 0) lq--;
    if (z == ONLY_REM)
    {
      avma=av; if (!si) return gen_0;
      r=cgeti(3);
      r[1] = evalsigne(sx) | evallgefint(3);
      r[2] = si; return r;
    }
    q[1] = evalsigne(sy) | evallgefint(lq);
    if (!z) return q;
    if (!si) { *z=gen_0; return q; }
    r=cgeti(3);
    r[1] = evalsigne(sx) | evallgefint(3);
    r[2] = si; *z=r; return q;
  }
  if (z == ONLY_REM)
  {
    ulong lr = lgefint(y);
    ulong lq = lgefint(x)-lgefint(y)+3;
    GEN r = cgeti(lr);
    GEN q = cgeti(lq);
    mpn_tdiv_qr(LIMBS(q), LIMBS(r),0, LIMBS(x), NLIMBS(x), LIMBS(y), NLIMBS(y));
    if (!r[lr - 1])
    {
      while(lr>2 && !r[lr - 1]) lr--;
      if (lr == 2) {avma=av; return gen_0;} /* exact division */
    }
    r[1] = evalsigne(sx) | evallgefint(lr);
    avma = (pari_sp) r; return r; 
  }
  else
  {
    ulong lq = lgefint(x)-lgefint(y)+3;
    ulong lr = lgefint(y);
    GEN q = cgeti(lq);
    GEN r = cgeti(lr);
    mpn_tdiv_qr(LIMBS(q), LIMBS(r),0, LIMBS(x), NLIMBS(x), LIMBS(y), NLIMBS(y));
    if (q[lq - 1] == 0) lq--;
    q[1] = evalsigne(sy) | evallgefint(lq);
    if (!z) { avma = (pari_sp)q; return q; }
    if (!r[lr - 1])
    {
      while(lr>2 && !r[lr - 1]) lr--;
      if (lr == 2) {avma=(pari_sp) q; *z=gen_0; return q;} /* exact division */
    }
    r[1] = evalsigne(sx) | evallgefint(lr);
    avma = (pari_sp) r; *z = r; return q; 
  }
}

/* assume y > x > 0. return y mod x */

static ulong
resiu(GEN y, ulong x)
{
  return mpn_mod_1(LIMBS(y), NLIMBS(y), x);
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
  while (! *scratch) scratch++; /* strip leading zeroes */
  while (Te > scratch) *--Td = *--Te;
  k = ((GEN)av - Td) + 2;
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

/* assume y != 0 and the division is exact */
GEN
diviuexact(GEN x, ulong y)
{
  /*TODO: implement true exact division.*/
  return divii(x,utoi(y));
}

/* Find z such that x=y*z, knowing that y | x (unchecked)
 * Method: y0 z0 = x0 mod B = 2^BITS_IN_LONG ==> z0 = 1/y0 mod B.
 *    Set x := (x - z0 y) / B, updating only relevant words, and repeat */
GEN
diviiexact(GEN x, GEN y)
{
  /*TODO: use mpn_bdivmod instead*/
  return divii(x,y);
}

/********************************************************************/
/**                                                                **/
/**               INTEGER MULTIPLICATION                           **/
/**                                                                **/
/********************************************************************/

/* nx >= ny = num. of digits of x, y (not GEN, see mulii) */
GEN
muliispec(GEN x, GEN y, long nx, long ny)
{
  GEN zd;
  long lz;
  ulong hi;

  if (nx < ny) swapspec(x,y, nx,ny);
  if (!ny) return gen_0;
  if (ny == 1) return muluispec((ulong)*y, x, nx);
    
  lz = nx+ny+2;
  zd = cgeti(lz);
  hi = mpn_mul(LIMBS(zd), (mp_limb_t *)x, nx, (mp_limb_t *)y, ny);
  if (!hi) lz--;
  /*else zd[lz-1]=hi; GH tell me it is not necessary.*/

  zd[1] = evalsigne(1) | evallgefint(lz);
  return zd;
}

GEN
sqrispec(GEN x, long nx)
{
  GEN zd;
  long lz;

  if (!nx) return gen_0;
  if (nx==1) return muluu(*x,*x);
    
  lz = (nx<<1)+2;
  zd = cgeti(lz);
  mpn_mul_n(LIMBS(zd), (mp_limb_t *)x, (mp_limb_t *)x, nx);
  if (zd[lz-1]==0) lz--;

  zd[1] = evalsigne(1) | evallgefint(lz);
  return zd;
}

/* x % (2^n), assuming x, n >= 0 */
GEN
resmod2n(GEN x, long n)
{
  ulong hi;
  long l,k,lx,ly;
  GEN z, xd, zd;

  if (!signe(x) || !n) return gen_0;

  l = n & (BITS_IN_LONG-1);    /* n % BITS_IN_LONG */
  k = n >> TWOPOTBITS_IN_LONG; /* n / BITS_IN_LONG */
  lx = lgefint(x);
  if (lx < k+3) return icopy(x);

  xd = x + (2 + k);
  /* x = |k|...|1|#|... : copy the last l bits of # and the first k words
   *              ^--- initial xd  */
  hi = ((ulong)*xd) & ((1UL<<l)-1); /* last l bits of # = top bits of result */
  if (!hi)
  { /* strip leading zeroes from result */
    xd--; while (k && !*xd) { k--; xd--; }
    if (!k) return gen_0;
    ly = k+2;
  }
  else
    ly = k+3;

  zd = z = cgeti(ly);
  *++zd = evalsigne(1) | evallgefint(ly);
  xd = x+1; for ( ;k; k--) *++zd = *++xd;
  if (hi) *++zd = hi;
  return z;
}

/********************************************************************/
/**                                                                **/
/**                      INTEGER SQUARE ROOT                       **/
/**                                                                **/
/********************************************************************/

/* Return S (and set R) s.t S^2 + R = N, 0 <= R <= 2S.
 * As for dvmdii, R is last on stack and guaranteed to be gen_0 in case the 
 * remainder is 0. R = NULL is allowed. */
GEN
sqrtremi(GEN a, GEN *r)
{
  long l, na = NLIMBS(a);
  mp_size_t nr;
  GEN S;
  if (!na) {
    if (r) *r = gen_0;
    return gen_0;
  }
  l = (na + 5) >> 1; /* 2 + ceil(na/2) */
  S = cgeti(l); S[1] = evalsigne(1) | evallgefint(l);
  if (r) {
    GEN R = cgeti(2 + na);
    nr = mpn_sqrtrem(LIMBS(S), LIMBS(R), LIMBS(a), na);
    if (nr) R[1] = evalsigne(1) | evallgefint(nr+2);
    else    { avma = (pari_sp)S; R = gen_0; }
    *r = R;
  }
  else
    (void)mpn_sqrtrem(LIMBS(S), NULL, LIMBS(a), na);
  return S;
}

/* compute sqrt(|a|), assuming a != 0 */
GEN
sqrtr_abs(GEN a)
{
  GEN res;
  mp_limb_t *b, *c;
  long l = RNLIMBS(a), e = expo(a), er = e>>1;
  long n;
  res = cgetr(2 + l);
  res[1] = evalsigne(1) | evalexpo(er);
  if (e&1)
  {
    b = (mp_limb_t *) new_chunk(l<<1);
    xmpn_zero(b,l);
    xmpn_mirrorcopy(b+l, RLIMBS(a), l);
    c = (mp_limb_t *) new_chunk(l);
    n = mpn_sqrtrem(c,b,b,l<<1); /* c <- sqrt; b <- rem */
    if (n>l || (n==l && mpn_cmp(b,c,l) > 0)) mpn_add_1(c,c,l,1);
  }
  else
  {
    ulong u;
    b = (mp_limb_t *) ishiftr_lg(a,l+2,-1);
    b[1] = ((ulong)a[l+1])<<(BITS_IN_LONG-1);
    b = (mp_limb_t *) new_chunk(l);
    xmpn_zero(b,l+1); /* overwrites the former b[0] */
    c = (mp_limb_t *) new_chunk(l + 1);
    n = mpn_sqrtrem(c,b,b,(l<<1)+2); /* c <- sqrt; b <- rem */
    u = (ulong)*c++;
    if ( u&HIGHBIT || (u == ~HIGHBIT && 
             (n>l || (n==l && mpn_cmp(b,c,l)>0))))
      mpn_add_1(c,c,l,1);
  }
  xmpn_mirrorcopy(RLIMBS(res),c,l);
  avma = (pari_sp)res; return res;
}

/* Normalize a non-negative integer.  */
GEN
int_normalize(GEN x, long known_zero_words)
{
  long i =  lgefint(x) - 1 - known_zero_words;
  for ( ; i > 1; i--)
    if (x[i]) { setlgefint(x, i+1); return x; }
  x[1] = evalsigne(0) | evallgefint(2); return x;
}
