/* $Id: base3.c 8269 2006-12-11 14:29:52Z kb $

Copyright (C) 2000  The PARI group.

This file is part of the PARI/GP package.

PARI/GP is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation. It is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY WHATSOEVER.

Check the License for details. You should have received a copy of it, along
with the package; see the file 'COPYING'. If not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA. */

/*******************************************************************/
/*                                                                 */
/*                       BASIC NF OPERATIONS                       */
/*                                                                 */
/*******************************************************************/
#include "pari.h"
#include "paripriv.h"

/*******************************************************************/
/*                                                                 */
/*                OPERATIONS OVER NUMBER FIELD ELEMENTS.           */
/*  represented as column vectors over the integral basis nf[7]    */
/*                                                                 */
/*******************************************************************/

int
RgV_isscalar(GEN x)
{
  long lx = lg(x),i;
  for (i=2; i<lx; i++)
    if (!gcmp0(gel(x, i))) return 0;
  return 1;
}
int
isnfscalar(GEN x) { return typ(x) == t_COL? RgV_isscalar(x): 0; }

static GEN
scal_mul(GEN nf, GEN x, GEN y, long ty)
{
  pari_sp av=avma, tetpil;
  GEN p1;

  if (!is_extscalar_t(ty))
  {
    if (ty!=t_COL) pari_err(typeer,"nfmul");
    y = coltoliftalg(nf, y);
  }
  p1 = gmul(x,y); tetpil=avma;
  return gerepile(av,tetpil,algtobasis(nf,p1));
}

static GEN
get_tab(GEN nf, long *N)
{
  GEN tab = (typ(nf) == t_MAT)? nf: gel(nf,9);
  *N = lg(tab[1])-1; return tab;
}

/* x != 0 t_INT. Return x * y (not memory clean) */
static GEN
_mulix(GEN x, GEN y) {
  return is_pm1(x)? (signe(x) < 0)? gneg(y): y
                  : gmul(x, y);
}
/* x != 0, y t_INT. Return x * y (not memory clean) */
static GEN
_mulii(GEN x, GEN y) {
  return is_pm1(x)? (signe(x) < 0)? negi(y): y
                  : mulii(x, y);
}

/* compute xy as ( sum_i x_i sum_j y_j m^{i,j}_k )_k.
 * Assume tab in M_{N x N^2}(Z), with x, y R^N (R arbitrary) */
static GEN
mul_by_tabi(GEN tab, GEN x, GEN y)
{
  long i, j, k, N = lg(x)-1;
  GEN s, v = cgetg(N+1,t_COL);

  for (k=1; k<=N; k++)
  {
    pari_sp av = avma;
    if (k == 1)
      s = gmul(gel(x,1),gel(y,1));
    else
      s = gadd(gmul(gel(x,1),gel(y,k)),
               gmul(gel(x,k),gel(y,1)));
    for (i=2; i<=N; i++)
    {
      GEN t, xi = gel(x,i);
      long base;
      if (gcmp0(xi)) continue;

      base = (i-1)*N;
      t = NULL;
      for (j=2; j<=N; j++)
      {
	GEN p1, c = gcoeff(tab,k,base+j); /* m^{i,j}_k */
	if (!signe(c)) continue;
        p1 = _mulix(c, gel(y,j));
        t = t? gadd(t, p1): p1;
      }
      if (t) s = gadd(s, gmul(xi, t));
    }
    gel(v,k) = gerepileupto(av,s);
  }
  return v;
}

/* product of x and y in nf */
GEN
element_mul(GEN nf, GEN x, GEN y)
{
  long N,tx,ty;
  GEN tab;

  if (x == y) return element_sqr(nf,x);

  tx=typ(x); ty=typ(y);
  nf=checknf(nf);
  if (tx==t_POLMOD) x=checknfelt_mod(nf,x,"element_mul");
  if (ty==t_POLMOD) y=checknfelt_mod(nf,y,"element_mul");
  if (is_extscalar_t(tx)) return scal_mul(nf,x,y,ty);
  if (is_extscalar_t(ty)) return scal_mul(nf,y,x,tx);
  if (tx != t_COL || ty != t_COL) pari_err(typeer,"element_mul");
  if (RgV_isscalar(x)) return gmul(gel(x,1),y);
  if (RgV_isscalar(y)) return gmul(gel(y,1),x);

  tab = get_tab(nf, &N);
  return mul_by_tabi(tab,x,y);
}

/* inverse of x in nf */
GEN
element_inv(GEN nf, GEN x)
{
  pari_sp av=avma;
  long i,N,tx=typ(x);
  GEN p1;

  nf=checknf(nf); N=degpol(nf[1]);
  if (is_extscalar_t(tx))
  {
    if (tx==t_POLMOD) (void)checknfelt_mod(nf,x,"element_inv");
    else if (tx==t_POL) x=gmodulo(x,gel(nf,1));
    return gerepileupto(av, algtobasis(nf, ginv(x)));
  }
  if (tx != t_COL) pari_err(typeer,"element_inv");
  if (RgV_isscalar(x))
  {
    p1=cgetg(N+1,t_COL); gel(p1,1) = ginv(gel(x,1));
    for (i=2; i<=N; i++) gel(p1,i) = gcopy(gel(x,i));
    return p1;
  }
  p1 = QXQ_inv(coltoliftalg(nf,x), gel(nf,1));
  return gerepileupto(av, poltobasis(nf,p1));
}

/* quotient of x and y in nf */
GEN
element_div(GEN nf, GEN x, GEN y)
{
  pari_sp av=avma;
  long tx = typ(x), ty = typ(y);
  GEN p1;

  nf=checknf(nf);
  if (tx==t_POLMOD) (void)checknfelt_mod(nf,x,"element_div");
  else if (tx==t_POL) x=gmodulo(x,gel(nf,1));

  if (ty==t_POLMOD) (void)checknfelt_mod(nf,y,"element_div");
  else if (ty==t_POL) y=gmodulo(y,gel(nf,1));

  if (is_extscalar_t(tx))
  {
    if (is_extscalar_t(ty)) p1=gdiv(x,y);
    else
    {
      if (ty!=t_COL) pari_err(typeer,"nfdiv");
      p1 = gdiv(x, coltoalg(nf, y));
    }
    return gerepileupto(av, algtobasis(nf,p1));
  }
  if (is_extscalar_t(ty))
  {
    if (tx!=t_COL) pari_err(typeer,"nfdiv");
    p1 = gdiv(coltoalg(nf,x), y);
    return gerepileupto(av, algtobasis(nf,p1));
  }
  if (tx != t_COL || ty != t_COL) pari_err(typeer,"element_div");

  if (RgV_isscalar(y)) return gdiv(x,gel(y,1));
  if (RgV_isscalar(x))
  {
    p1=element_inv(nf,y);
    return gerepileupto(av, gmul(gel(x,1),p1));
  }

  p1 = gmul(coltoliftalg(nf,x), QXQ_inv(coltoliftalg(nf,y), gel(nf,1)));
  p1 = RgX_rem(p1, gel(nf,1));
  return gerepileupto(av, poltobasis(nf,p1));
}

/* product of INTEGERS (i.e vectors with integral coeffs) x and y in nf */
GEN
element_muli(GEN nf, GEN x, GEN y)
{
  long i, j, k, N, tx = typ(x), ty = typ(y);
  GEN s, v, tab = get_tab(nf, &N);

  if (tx == t_INT) { return ty == t_INT? gscalcol(mulii(x,y), N): gmul(x, y); }
  if (tx != t_COL || lg(x) != N+1
   || ty != t_COL || lg(y) != N+1) pari_err(typeer,"element_muli");
  v = cgetg(N+1,t_COL);
  for (k=1; k<=N; k++)
  {
    pari_sp av = avma;
    if (k == 1)
      s = mulii(gel(x,1),gel(y,1));
    else
      s = addii(mulii(gel(x,1),gel(y,k)),
                mulii(gel(x,k),gel(y,1)));
    for (i=2; i<=N; i++)
    {
      GEN t, xi = gel(x,i);
      long base;
      if (!signe(xi)) continue;

      base = (i-1)*N;
      t = NULL;
      for (j=2; j<=N; j++)
      {
	GEN p1, c = gcoeff(tab,k,base+j); /* m^{i,j}_k */
	if (!signe(c)) continue;
        p1 = _mulii(c, gel(y,j));
        t = t? addii(t, p1): p1;
      }
      if (t) s = addii(s, mulii(xi, t));
    }
    gel(v,k) = gerepileuptoint(av,s);
  }
  return v;
}

/* product of INTEGERS (i.e vectors with integral coeffs) x and y in nf */
GEN
element_sqri(GEN nf, GEN x)
{
  long i, j, k, N;
  GEN s, v, tab = get_tab(nf, &N);

  v = cgetg(N+1,t_COL);
  for (k=1; k<=N; k++)
  {
    pari_sp av = avma;
    if (k == 1)
      s = sqri(gel(x,1));
    else
      s = shifti(mulii(gel(x,1),gel(x,k)), 1);
    for (i=2; i<=N; i++)
    {
      GEN p1, c, t, xi = gel(x,i);
      long base;
      if (!signe(xi)) continue;

      base = (i-1)*N;
      c = gcoeff(tab,k,base+i);
      t = signe(c)? _mulii(c,xi): NULL;
      for (j=i+1; j<=N; j++)
      {
	c = gcoeff(tab,k,base+j);
	if (!signe(c)) continue;
        p1 = _mulii(c, shifti(gel(x,j),1));
        t = t? addii(t, p1): p1;
      }
      if (t) s = addii(s, mulii(xi, t));
    }
    gel(v,k) = gerepileuptoint(av,s);
  }
  return v;
}

/* cf mul_by_tabi */
static GEN
sqr_by_tabi(GEN tab, GEN x)
{
  long i, j, k, N = lg(x)-1;
  GEN s, v = cgetg(N+1,t_COL);

  for (k=1; k<=N; k++)
  {
    pari_sp av = avma;
    if (k == 1)
      s = gsqr(gel(x,1));
    else
      s = gmul2n(gmul(gel(x,1),gel(x,k)), 1);
    for (i=2; i<=N; i++)
    {
      GEN p1, c, t, xi = gel(x,i);
      long base;
      if (gcmp0(xi)) continue;

      base = (i-1)*N;
      c = gcoeff(tab,k,base+i);
      t = signe(c)? _mulix(c,xi): NULL;
      for (j=i+1; j<=N; j++)
      {
	c = gcoeff(tab,k,(i-1)*N+j);
	if (!signe(c)) continue;
        p1 = gmul(shifti(c,1), gel(x,j));
        t = t? gadd(t, p1): p1;
      }
      if (t) s = gadd(s, gmul(xi, t));
    }
    gel(v,k) = gerepileupto(av,s);
  }
  return v;
}

/* cf sqr_by_tab. Assume nothing about tab */
GEN
sqr_by_tab(GEN tab, GEN x)
{
  long i, j, k, N = lg(x)-1;
  GEN s, v = cgetg(N+1,t_COL);

  for (k=1; k<=N; k++)
  {
    pari_sp av = avma;
    if (k == 1)
      s = gsqr(gel(x,1));
    else
      s = gmul2n(gmul(gel(x,1),gel(x,k)), 1);
    for (i=2; i<=N; i++)
    {
      GEN p1, c, t, xi = gel(x,i);
      long base;
      if (gcmp0(xi)) continue;

      base = (i-1)*N;
      c = gcoeff(tab,k,base+i);
      t = !gcmp0(c)? gmul(c,xi): NULL;
      for (j=i+1; j<=N; j++)
      {
	c = gcoeff(tab,k,(i-1)*N+j);
	if (gcmp0(c)) continue;
        p1 = gmul(gmul2n(c,1), gel(x,j));
        t = t? gadd(t, p1): p1;
      }
      if (t) s = gadd(s, gmul(xi, t));
    }
    gel(v,k) = gerepileupto(av,s);
  }
  return v;
}

/* square of x in nf */
GEN
element_sqr(GEN nf, GEN x)
{
  long N, tx = typ(x);
  GEN tab;

  nf = checknf(nf);
  if (tx==t_POLMOD) x=checknfelt_mod(nf,x,"element_sqr");
  if (is_extscalar_t(tx))
  {
    pari_sp av = avma;
    return gerepileupto(av, algtobasis(nf, gsqr(x)));
  }
  if (tx != t_COL) pari_err(typeer,"element_sqr");

  tab = get_tab(nf, &N);
  return sqr_by_tabi(tab,x);
}

static GEN
_mul(void *data, GEN x, GEN y) { return element_mul((GEN)data,x,y); }
static GEN
_sqr(void *data, GEN x) { return element_sqr((GEN)data,x); }

/* Compute x^n in nf, left-shift binary powering */
GEN
element_pow(GEN nf, GEN x, GEN n)
{
  pari_sp av = avma;
  long s,N;
  GEN y, cx;

  if (typ(n)!=t_INT) pari_err(talker,"not an integer exponent in nfpow");
  nf=checknf(nf); N=degpol(nf[1]);
  s=signe(n); if (!s) return gscalcol_i(gen_1,N);
  if (typ(x) != t_COL)
  {
    x = algtobasis(nf,x);
    if (typ(x) != t_COL) pari_err(typeer,"element_pow");
  }

  if (RgV_isscalar(x))
  {
    y = gscalcol_i(gen_1,N);
    gel(y,1) = powgi(gel(x,1),n); return y;
  }
  x = primitive_part(x, &cx);
  y = leftright_pow(x, n, (void*)nf, _sqr, _mul);
  if (s<0) y = element_inv(nf, y);
  if (cx) y = gmul(y, powgi(cx, n));
  return av==avma? gcopy(y): gerepileupto(av,y);
}

typedef struct {
  GEN nf, p;
  long I;
} eltmod_muldata;

static GEN
_mulidmod(void *data, GEN x, GEN y)
{
  eltmod_muldata *D = (eltmod_muldata*)data;
  (void)y; /* ignored */
  return FpC_red(element_mulid(D->nf, x, D->I), D->p);
}
static GEN
_sqrmod(void *data, GEN x)
{
  eltmod_muldata *D = (eltmod_muldata*)data;
  return FpC_red(element_sqri(D->nf, x), D->p);
}

/* x = I-th vector of the Z-basis of Z_K, in Z^n, compute lift(x^n mod p) */
GEN
element_powid_mod_p(GEN nf, long I, GEN n, GEN p)
{
  pari_sp av = avma;
  eltmod_muldata D;
  long s,N;
  GEN y;

  if (typ(n) != t_INT) pari_err(talker,"not an integer exponent in nfpow");
  nf = checknf(nf); N = degpol(nf[1]);
  s = signe(n);
  if (s < 0) pari_err(talker,"negative power in element_powid_mod_p");
  if (!s || I == 1) return gscalcol_i(gen_1,N);
  D.nf = nf;
  D.p = p;
  D.I = I;
  y = leftright_pow(col_ei(N, I), n, (void*)&D, &_sqrmod, &_mulidmod);
  return gerepileupto(av,y);
}

GEN
element_mulidid(GEN nf, long i, long j)
{
  long N;
  GEN tab = get_tab(nf, &N);
  tab += (i-1)*N; return gel(tab,j);
}

/* Outputs x.w_i, where w_i is the i-th elt of the integral basis */
GEN
element_mulid(GEN nf, GEN x, long i)
{
  long j, k, N;
  GEN v, tab;

  if (i==1) return gcopy(x);
  tab = get_tab(nf, &N);
  if (typ(x) != t_COL || lg(x) != N+1) pari_err(typeer,"element_mulid");
  tab += (i-1)*N;
  v = cgetg(N+1,t_COL);
  for (k=1; k<=N; k++)
  {
    pari_sp av = avma;
    GEN s = gen_0;
    for (j=1; j<=N; j++)
    {
      GEN c = gcoeff(tab,k,j);
      if (signe(c)) s = gadd(s, _mulix(c, gel(x,j)));
    }
    gel(v,k) = gerepileupto(av,s);
  }
  return v;
}

/* table of multiplication by wi in ZK = Z[w1,..., wN] */
GEN
eltmulid_get_table(GEN nf, long i)
{
  long k,N;
  GEN m, tab = get_tab(nf, &N);
  tab += (i-1)*N;
  m = cgetg(N+1,t_COL);
  for (k=1; k<=N; k++) m[k] = tab[k];
  return m;
}

/* table of multiplication by x in ZK = Z[w1,..., wN] */
GEN
eltmul_get_table(GEN nf, GEN x)
{
  if (typ(x) == t_MAT) return x;
  else
  {
    long i, N = degpol(nf[1]);
    GEN mul = cgetg(N+1,t_MAT);
    x = algtobasis_i(nf, x);
    if (RgV_isscalar(x)) return gscalmat(gel(x,1), N);
    gel(mul,1) = x; /* assume w_1 = 1 */
    for (i=2; i<=N; i++) gel(mul,i) = element_mulid(nf,x,i);
    return mul;
  }
}

/* valuation of integer x, with resp. to prime ideal P above p.
 * p.P^(-1) = b Z_K, v >= val_p(norm(x)), and N = deg(nf)
 * [b may be given as the 'multiplication by b' matrix]
 */
long
int_elt_val(GEN nf, GEN x, GEN p, GEN b, GEN *newx)
{
  long i,k,w, N = degpol(nf[1]);
  GEN r,a,y,mul;

  mul = (typ(b) == t_MAT)? b: eltmul_get_table(nf, b);
  y = cgetg(N+1, t_COL); /* will hold the new x */
  x = shallowcopy(x);
  for(w=0;; w++)
  {
    for (i=1; i<=N; i++)
    { /* compute (x.b)_i */
      a = mulii(gel(x,1), gcoeff(mul,i,1));
      for (k=2; k<=N; k++) a = addii(a, mulii(gel(x,k), gcoeff(mul,i,k)));
      /* is it divisible by p ? */
      gel(y,i) = dvmdii(a,p,&r);
      if (signe(r))
      {
        if (newx) *newx = x;
        return w;
      }
    }
    r=x; x=y; y=r;
  }
}

long
element_val(GEN nf, GEN x, GEN vp)
{
  pari_sp av = avma;
  long w, vcx, e;
  GEN cx,p;

  if (gcmp0(x)) return VERYBIGINT;
  nf = checknf(nf);
  checkprimeid(vp);
  p = gel(vp,1);
  e = itos(gel(vp,3));
  switch(typ(x))
  {
    case t_INT: return e*Z_pval(x,p);
    case t_FRAC:return e*(Z_pval(gel(x,1),p) - Z_pval(gel(x,2),p));
    default: x = algtobasis_i(nf,x); break;
  }
  if (RgV_isscalar(x)) return e*ggval(gel(x,1),p);

  cx = content(x);
  if (gcmp1(cx)) vcx=0; else { x = gdiv(x,cx); vcx = ggval(cx,p); }
  w = int_elt_val(nf,x,p,gel(vp,5),NULL);
  avma = av; return w + vcx*e;
}

/* polegal without comparing variables */
long
polegal_spec(GEN x, GEN y)
{
  long i = lg(x);

  if (i != lg(y)) return 0;
  for (i--; i > 1; i--)
    if (!gequal(gel(x,i),gel(y,i))) return 0;
  return 1;
}

GEN
coltoalg(GEN nf, GEN x)
{
  return mkpolmod( coltoliftalg(nf, x), gel(nf,1) );
}

GEN
basistoalg(GEN nf, GEN x)
{
  long tx=typ(x),lx=lg(x),i,j,l;
  GEN z;

  nf = checknf(nf);
  switch(tx)
  {
    case t_COL:
      for (i=1; i<lx; i++)
      {
        long t = typ(x[i]);
	if (is_matvec_t(t)) break;
      }
      if (i==lx) {
        pari_sp av = avma;
        return gerepilecopy(av, coltoalg(nf, x));
      }
      /* fall through */

    case t_VEC: z=cgetg(lx,tx);
      for (i=1; i<lx; i++) gel(z,i) = basistoalg(nf,gel(x,i));
      return z;
    case t_MAT: z=cgetg(lx,t_MAT);
      if (lx == 1) return z;
      l = lg(x[1]);
      for (j=1; j<lx; j++)
      {
        gel(z,j) = cgetg(l,t_COL);
        for (i=1; i<l; i++) gcoeff(z,i,j) = basistoalg(nf,gcoeff(x,i,j));
      }
      return z;

    case t_POLMOD:
      if (!polegal_spec(gel(nf,1),gel(x,1)))
	pari_err(talker,"not the same number field in basistoalg");
      return gcopy(x);
    default: z=cgetg(3,t_POLMOD);
      gel(z,1) = gcopy(gel(nf,1));
      gel(z,2) = gtopoly(x, varn(nf[1])); return z;
  }
}

/* FIXME: basistoalg and algtobasis should not treat recursive objects
 * since t_COLs are ambiguous, and the functionnality is almost useless.
 * (Even then, matbasistoalg and matalgtobasis can be used instead.)
 * The following shallow functions do what the public functions should do,
 * without sanity checks. 
 * Assume nf is a genuine nf. */
GEN
basistoalg_i(GEN nf, GEN x)
{ return typ(x) == t_COL? coltoalg(nf, x): x; }
GEN
algtobasis_i(GEN nf, GEN x)
{
  switch(typ(x))
  {
    case t_INT: case t_FRAC:
      return gscalcol_i(x, degpol( gel(nf,1) ));
    case t_POLMOD:
      x = gel(x,2);
      if (typ(x) != t_POL) return gscalcol_i(x, degpol( gel(nf,1) ));
      /* fall through */
    case t_POL:
      return poltobasis(nf,x);
    case t_COL:
      if (lg(x) == lg(gel(nf,7))) break;
    default: pari_err(typeer,"algtobasis_i");
  }
  return x;
}
GEN
algtobasis_cp(GEN nf, GEN x)
{ return typ(x) == t_COL? gcopy(x): algtobasis(nf, x); }

/* gmul(A, RgX_to_RgV(x)), A t_MAT (or t_VEC) of compatible dimensions */
GEN
mulmat_pol(GEN A, GEN x)
{
  long i,l;
  GEN z;
  if (typ(x) != t_POL) return gmul(x,gel(A,1)); /* scalar */
  l=lg(x)-1; if (l == 1) return typ(A)==t_VEC? gen_0: zerocol(lg(A[1])-1);
  x++; z = gmul(gel(x,1), gel(A,1));
  for (i=2; i<l ; i++)
    if (!gcmp0(gel(x,i))) z = gadd(z, gmul(gel(x,i), gel(A,i)));
  return z;
}

/* valid for t_POL, nf a genuine nf or an rnf !
 * No garbage collecting. No check.  */
GEN
poltobasis(GEN nf, GEN x)
{
  GEN P = gel(nf,1);
  if (degpol(x) >= degpol(P)) x = RgX_rem(x,P);
  return mulmat_pol(gel(nf,8), x);
}

GEN
algtobasis(GEN nf, GEN x)
{
  long tx=typ(x),lx=lg(x),i,N;
  pari_sp av=avma;
  GEN z;

  nf = checknf(nf);
  switch(tx)
  {
    case t_VEC: case t_COL: case t_MAT:
      z=cgetg(lx,tx);
      for (i=1; i<lx; i++) gel(z,i) = algtobasis(nf,gel(x,i));
      return z;
    case t_POLMOD:
      if (!polegal_spec(gel(nf,1),gel(x,1)))
	pari_err(talker,"not the same number field in algtobasis");
      x = gel(x,2);
      if (typ(x) != t_POL) break;
      /* fall through */
    case t_POL:
      if (varn(x) != varn(gel(nf,1)))
        pari_err(talker,"incompatible variables in algtobasis");
      return gerepileupto(av,poltobasis(nf,x));

  }
  N=degpol(nf[1]); return gscalcol(x,N);
}

GEN
rnfbasistoalg(GEN rnf,GEN x)
{
  long tx = typ(x), lx = lg(x), i;
  pari_sp av = avma;
  GEN p1, z, nf;

  checkrnf(rnf);
  switch(tx)
  {
    case t_VEC: case t_COL:
      p1 = cgetg(lx,t_COL); nf = gel(rnf,10);
      for (i=1; i<lx; i++) gel(p1,i) = basistoalg_i(nf, gel(x,i));
      p1 = gmul(gmael(rnf,7,1), p1);
      return gerepileupto(av, gmodulo(p1,gel(rnf,1)));

    case t_MAT:
      z = cgetg(lx,tx);
      for (i=1; i<lx; i++) gel(z,i) = rnfbasistoalg(rnf,gel(x,i));
      return z;

    case t_POLMOD:
      return gcopy(x);

    default: z = cgetg(3,t_POLMOD);
      gel(z,1) = gcopy(gel(rnf,1));
      gel(z,2) = gmul(x,pol_1[varn(rnf[1])]); return z;
  }
}

GEN
matbasistoalg(GEN nf,GEN x)
{
  long i, j, li, lx = lg(x);
  GEN c, z = cgetg(lx,t_MAT);

  if (typ(x) != t_MAT) pari_err(talker,"not a matrix in matbasistoalg");
  if (lx == 1) return z;
  li = lg(x[1]);
  for (j=1; j<lx; j++)
  {
    c = cgetg(li,t_COL); gel(z,j) = c;
    for (i=1; i<li; i++) gel(c,i) = basistoalg(nf,gcoeff(x,i,j));
  }
  return z;
}

GEN
matalgtobasis(GEN nf,GEN x)
{
  long i, j, li, lx = lg(x);
  GEN c, z = cgetg(lx, t_MAT);

  if (typ(x) != t_MAT) pari_err(talker,"not a matrix in matalgtobasis");
  if (lx == 1) return z;
  li = lg(x[1]);
  for (j=1; j<lx; j++)
  {
    c = cgetg(li,t_COL); gel(z,j) = c;
    for (i=1; i<li; i++) gel(c,i) = algtobasis_cp(nf, gcoeff(x,i,j));
  }
  return z;
}

/* assume x is a t_POLMOD */
GEN
lift_to_pol(GEN x)
{
  GEN y = gel(x,2);
  return (typ(y) != t_POL)? scalarpol(y,varn(x[1])): y;
}

GEN
rnfalgtobasis(GEN rnf,GEN x)
{
  long tx = typ(x), i, lx;
  GEN z;

  checkrnf(rnf);
  switch(tx)
  {
    case t_VEC: case t_COL: case t_MAT:
      lx = lg(x); z = cgetg(lx,tx);
      for (i=1; i<lx; i++) gel(z,i) = rnfalgtobasis(rnf,gel(x,i));
      return z;

    case t_POLMOD:
      if (!polegal_spec(gel(rnf,1),gel(x,1)))
	pari_err(talker,"not the same number field in rnfalgtobasis");
      x = lift_to_pol(x); /* fall through */
    case t_POL: {
      pari_sp av = avma;
      return gerepileupto(av, poltobasis(rnf, x));
    }
  }
  return gscalcol(x, degpol(rnf[1]));
}

/* Given a and b in nf, gives an algebraic integer y in nf such that a-b.y
 * is "small"
 */
GEN
nfdiveuc(GEN nf, GEN a, GEN b)
{
  pari_sp av=avma, tetpil;
  a = element_div(nf,a,b); tetpil=avma;
  return gerepile(av,tetpil,ground(a));
}

/* Given a and b in nf, gives a "small" algebraic integer r in nf
 * of the form a-b.y
 */
GEN
nfmod(GEN nf, GEN a, GEN b)
{
  pari_sp av=avma,tetpil;
  GEN p1=gneg_i(element_mul(nf,b,ground(element_div(nf,a,b))));
  tetpil=avma; return gerepile(av,tetpil,gadd(a,p1));
}

/* Given a and b in nf, gives a two-component vector [y,r] in nf such
 * that r=a-b.y is "small". */
GEN
nfdivrem(GEN nf, GEN a, GEN b)
{
  pari_sp av = avma;
  GEN p1,z, y = ground(element_div(nf,a,b));

  p1 = gneg_i(element_mul(nf,b,y));
  z = cgetg(3,t_VEC);
  gel(z,1) = gcopy(y);
  gel(z,2) = gadd(a,p1); return gerepileupto(av, z);
}

/*************************************************************************/
/**									**/
/**			      (Z_K/I)^*					**/
/**									**/
/*************************************************************************/
/* return sign(sigma_k(x)), x t_COL (integral, primitive) */
static long
eval_sign(GEN M, GEN x, long k)
{
  long i, l = lg(x);
  GEN z = mpmul(gcoeff(M,k,1), gel(x,1));
  for (i = 2; i < l; i++)
    z = mpadd(z, mpmul(gcoeff(M,k,i), gel(x,i)));
  if (lg(z) < DEFAULTPREC) pari_err(precer,"zsigne");
  return signe(z);
}

GEN
arch_to_perm(GEN arch)
{
  long i, k, l;
  GEN perm;

  if (!arch) return cgetg(1, t_VECSMALL);
  switch (typ(arch))
  {
   case t_VECSMALL: return arch;
   case t_VEC: break;
   default: pari_err(typeer,"arch_to_perm");
  }
  l = lg(arch);
  perm = cgetg(l, t_VECSMALL);
  for (k=1, i=1; i < l; i++)
    if (signe(arch[i])) perm[k++] = i;
  setlg(perm, k); return perm;
}
GEN
perm_to_arch(GEN nf, GEN archp)
{
  long i, l;
  GEN v;
  if (typ(archp) == t_VEC) return archp;

  l = lg(archp); nf = checknf(nf);
  v = zerovec( nf_get_r1(nf) );
  for (i = 1; i < l; i++) gel(v, archp[i]) = gen_1;
  return v;
}

/* reduce mod 2 in place */
GEN
F2V_red_ip(GEN v)
{
  long i, l = lg(v);
  for (i = 1; i < l; i++) gel(v,i) = mpodd(gel(v,i))? gen_1: gen_0;
  return v;
}

/* return (column) vector of R1 signatures of x (0 or 1)
 * if arch = NULL, assume arch = [0,..0] */
GEN
zsigne(GEN nf,GEN x,GEN arch)
{
  GEN M, V, archp = arch_to_perm(arch);
  long i, s, l = lg(archp);
  pari_sp av;

  if (l == 1) return cgetg(1,t_COL);
  V = cgetg(l,t_COL); av = avma;
  nf = checknf(nf);
  switch(typ(x))
  {
    case t_MAT: /* factorisation */
    {
      GEN g = gel(x,1), e = gel(x,2), z = vec_setconst(V, gen_0);
      for (i=1; i<lg(g); i++)
        if (mpodd(gel(e,i))) z = gadd(z, zsigne(nf,gel(g,i),archp));
      for (i=1; i<l; i++) gel(V,i) = mpodd(gel(z,i))? gen_1: gen_0;
      avma = av; return V;
    }
    case t_POLMOD: x = gel(x,2);      /* fall through */
    case t_POL: x = algtobasis(nf, x); /* fall through */
    case t_COL: if (!RgV_isscalar(x)) break;
                x = gel(x,1);         /* fall through */
    case t_INT: case t_FRAC:
      s = gsigne(x); if (!s) pari_err(talker,"zero element in zsigne");
      return vec_setconst(V, (s < 0)? gen_1: gen_0);
  }
  x = Q_primpart(x); M = gmael(nf,5,1);
  for (i = 1; i < l; i++)
    gel(V,i) = (eval_sign(M, x, archp[i]) > 0)? gen_0: gen_1;
  avma = av; return V;
}

/* return the t_COL vector of signs of x; the matrix of such if x is a vector
 * of nf elements */
GEN
zsigns(GEN nf, GEN x)
{
  long r1, i, l;
  GEN arch, S;

  nf = checknf(nf); r1 = nf_get_r1(nf);
  arch = cgetg(r1+1, t_VECSMALL); for (i=1; i<=r1; i++) arch[i] = i;
  if (typ(x) != t_VEC) return zsigne(nf, x, arch);
  l = lg(x); S = cgetg(l, t_MAT);
  for (i=1; i<l; i++) gel(S,i) = zsigne(nf, gel(x,i), arch);
  return S;
}

/* For internal use. Reduce x modulo (invertible) y */
GEN
close_modinvertible(GEN x, GEN y)
{
  return gmul(y, ground(gauss(y,x)));
}
GEN
reducemodinvertible(GEN x, GEN y)
{
  return gadd(x, gneg(close_modinvertible(x,y)));
}
GEN
lllreducemodmatrix(GEN x,GEN y)
{
  return reducemodinvertible(x, lllint_ip(y,4));
}

/* Reduce column x modulo y in HNF */
GEN
colreducemodHNF(GEN x, GEN y, GEN *Q)
{
  long i, l = lg(x);
  GEN q;

  if (Q) *Q = cgetg(l,t_COL);
  if (l == 1) return cgetg(1,t_COL);
  for (i = l-1; i>0; i--)
  {
    q = negi(diviiround(gel(x,i), gcoeff(y,i,i)));
    if (Q) gel(*Q, i) = q;
    if (signe(q)) x = gadd(x, gmul(q, gel(y,i)));
  }
  return x;
}

/* for internal use...Reduce matrix x modulo matrix y */
GEN
reducemodmatrix(GEN x, GEN y)
{
  return reducemodHNF(x, hnfmod(y,detint(y)), NULL);
}

/* x = y Q + R */
GEN
reducemodHNF(GEN x, GEN y, GEN *Q)
{
  long lx = lg(x), i;
  GEN R = cgetg(lx, t_MAT);
  if (Q)
  {
    GEN q = cgetg(lx, t_MAT); *Q = q;
    for (i=1; i<lx; i++) gel(R,i) = colreducemodHNF(gel(x,i),y,(GEN*)(q+i));
  }
  else
    for (i=1; i<lx; i++)
    {
      pari_sp av = avma;
      gel(R,i) = gerepileupto(av, colreducemodHNF(gel(x,i),y,NULL));
    }
  return R;
}

/* For internal use. Reduce x modulo ideal (assumed non-zero, in HNF). We
 * want a non-zero result */
GEN
nfreducemodideal_i(GEN x0,GEN ideal)
{
  GEN x = colreducemodHNF(x0, ideal, NULL);
  if (gcmp0(x)) return gel(ideal,1);
  return x == x0? gcopy(x) : x;
}
GEN
nfreducemodideal(GEN nf,GEN x0,GEN ideal)
{
  return nfreducemodideal_i(x0, idealhermite(nf,ideal));
}

/* multiply y by t = 1 mod^* f such that sign(x) = sign(y) at arch = idele[2].
 * If x == NULL, make y >> 0 at sarch */
GEN
set_sign_mod_idele(GEN nf, GEN x, GEN y, GEN idele, GEN sarch)
{
  GEN s, archp, gen;
  long nba,i;
  if (!sarch) return y;
  gen = gel(sarch,2); nba = lg(gen);
  if (nba == 1) return y;

  archp = arch_to_perm(gel(idele,2));
  s = zsigne(nf, y, archp);
  if (x) s = gadd(s, zsigne(nf, x, archp));
  s = gmul(gel(sarch,3), s);
  for (i=1; i<nba; i++)
    if (mpodd(gel(s,i))) y = element_mul(nf,y,gel(gen,i));
  return y;
}

/* compute elt = x mod idele, with sign(elt) = sign(x) at arch */
GEN
nfreducemodidele(GEN nf,GEN x,GEN idele,GEN sarch)
{
  GEN y;

  if (gcmp0(x)) return gcopy(x);
  if (!sarch || typ(idele)!=t_VEC || lg(idele)!=3)
    return nfreducemodideal(nf,x,idele);

  y = nfreducemodideal(nf,x,gel(idele,1));
  y = set_sign_mod_idele(nf, x, y, idele, sarch);
  return (gexpo(y) > gexpo(x))? x: y;
}

/* given an element x in Z_K and an integral ideal y with x, y coprime,
   outputs an element inverse of x modulo y */
GEN
element_invmodideal(GEN nf, GEN x, GEN y)
{
  pari_sp av = avma;
  GEN a, xh, yh;

  nf = checknf(nf);
  if (gcmp1(gcoeff(y,1,1))) return zerocol( degpol(nf[1]) );

  yh = get_hnfid(nf, y);
  switch (typ(x))
  {
    case t_POL: case t_POLMOD: case t_COL:
      xh = idealhermite_aux(nf,x); break;
    default: pari_err(typeer,"element_invmodideal");
      return NULL; /* not reached */
  }
  a = element_div(nf, hnfmerge_get_1(xh, yh), x);
  return gerepilecopy(av, nfreducemodideal_i(a, yh));
}

static GEN
element_sqrmodideal(GEN nf, GEN x, GEN id) {
  return nfreducemodideal_i(element_sqr(nf,x),id);
}
static GEN
element_mulmodideal(GEN nf, GEN x, GEN y, GEN id) {
  return x? nfreducemodideal_i(element_mul(nf,x,y),id): algtobasis_i(nf, y);
}
/* assume k >= 0, ideal in HNF */
GEN
element_powmodideal(GEN nf,GEN x,GEN k,GEN ideal)
{
  GEN y = NULL;
  for(;;)
  {
    if (mpodd(k)) y = element_mulmodideal(nf,y,x,ideal);
    k = shifti(k,-1); if (!signe(k)) break;
    x = element_sqrmodideal(nf,x,ideal);
  }
  return y? y: gscalcol_i(gen_1,degpol(nf[1]));
}

/* assume k >= 0, assume idele = [HNFideal, arch] */
GEN
element_powmodidele(GEN nf,GEN x,GEN k,GEN idele,GEN sarch)
{
  GEN y = element_powmodideal(nf,x,k,gel(idele,1));
  if (mpodd(k))
  { if (!gcmp1(k)) y = set_sign_mod_idele(nf, x, y, idele, sarch); }
  else
  { if (!gcmp0(k)) y = set_sign_mod_idele(nf, NULL, y, idele, sarch); }
  return y;
}

/* a * g^n mod id */
static GEN
elt_mulpow_modideal(GEN nf, GEN a, GEN g, GEN n, GEN id)
{
  return element_mulmodideal(nf, a, element_powmodideal(nf,g,n,id), id);
}

/* assume (num(g[i]), id) = 1 for all i. Return prod g[i]^e[i] mod id.
 * EX = multiple of exponent of (O_K/id)^* */
GEN
famat_to_nf_modideal_coprime(GEN nf, GEN g, GEN e, GEN id, GEN EX)
{
  GEN dh, h, n, plus = NULL, minus = NULL, idZ = gcoeff(id,1,1);
  long i, lx = lg(g);
  GEN EXo2 = (expi(EX) > 10)? shifti(EX,-1): NULL;

  if (is_pm1(idZ)) lx = 1; /* id = Z_K */
  for (i=1; i<lx; i++)
  {
    long sn;
    n = centermodii(gel(e,i), EX, EXo2);
    sn = signe(n); if (!sn) continue;

    h = Q_remove_denom(gel(g,i), &dh);
    if (dh) h = FpC_Fp_mul(h, Fp_inv(dh,idZ), idZ);
    if (sn > 0)
      plus = elt_mulpow_modideal(nf, plus, h, n, id);
    else /* sn < 0 */
      minus = elt_mulpow_modideal(nf, minus, h, negi(n), id);
  }
  if (minus)
    plus = element_mulmodideal(nf, plus, element_invmodideal(nf,minus,id), id);
  return plus? plus: gscalcol_i(gen_1, lg(id)-1);
}

/* given 2 integral ideals x, y in HNF s.t x | y | x^2, compute the quotient
   (1+x)/(1+y) in the form [[cyc],[gen]], if U != NULL, set *U := ux^-1 */
static GEN
zidealij(GEN x, GEN y, GEN *U)
{
  GEN G, p1, cyc;
  long j, N;

  /* x^(-1) y = relations between the 1 + x_i (HNF) */
  cyc = smithrel(hnf_gauss(x, y), U, &G);
  N = lg(cyc); G = gmul(x,G); settyp(G, t_VEC); /* new generators */
  for (j=1; j<N; j++)
  {
    p1 = gel(G,j);
    gel(p1,1) = addsi(1, gel(p1,1)); /* 1 + g_j */
  }
  if (U) *U = gmul(*U, ginv(x));
  return mkvec2(cyc, G);
}

/* smallest integer n such that g0^n=x modulo p prime. Assume g0 has order q
 * (p-1 if q = NULL) */
GEN
Fp_shanks(GEN x,GEN g0,GEN p, GEN q)
{
  pari_sp av=avma,av1,lim;
  long lbaby,i,k;
  GEN p1,smalltable,giant,perm,v,g0inv;

  x = modii(x,p);
  if (is_pm1(x) || equaliu(p,2)) { avma = av; return gen_0; }
  p1 = addsi(-1, p); if (!q) q = p1;
  if (equalii(p1,x)) { avma = av; return shifti(q,-1); }
  p1 = sqrti(q);
  if (cmpiu(p1,LGBITS) >= 0) pari_err(talker,"module too large in Fp_shanks");
  lbaby = itos(p1)+1; smalltable = cgetg(lbaby+1,t_VEC);
  g0inv = Fp_inv(g0, p); p1 = x;

  for (i=1;;i++)
  {
    av1 = avma;
    if (is_pm1(p1)) { avma = av; return stoi(i-1); }
    gel(smalltable,i) = p1; if (i==lbaby) break;
    p1 = gerepileuptoint(av1, remii(mulii(p1,g0inv), p));
  }
  giant = remii(mulii(x, Fp_inv(p1,p)), p);
  p1=cgetg(lbaby+1,t_VEC);
  perm = gen_sort(smalltable, cmp_IND | cmp_C, cmpii);
  for (i=1; i<=lbaby; i++) p1[i]=smalltable[perm[i]];
  smalltable=p1; p1=giant;

  av1 = avma; lim=stack_lim(av1,2);
  for (k=1;;k++)
  {
    i=tablesearch(smalltable,p1,cmpii);
    if (i)
    {
      v=addis(mulss(lbaby-1,k),perm[i]);
      return gerepileuptoint(av,addsi(-1,v));
    }
    p1 = remii(mulii(p1,giant), p);

    if (low_stack(lim, stack_lim(av1,2)))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"Fp_shanks, k = %ld", k);
      p1 = gerepileuptoint(av1, p1);
    }
  }
}

/* Pohlig-Hellman */
GEN
Fp_PHlog(GEN a, GEN g, GEN p, GEN ord)
{
  pari_sp av = avma;
  GEN v,t0,a0,b,q,g_q,n_q,ginv0,qj,ginv;
  GEN fa, ex;
  long e,i,j,l;

  if (equalii(g, a)) return gen_1; /* frequent special case */
  if (!ord) ord = subis(p,1);
  if (typ(ord) == t_MAT)
  {
    fa = ord;
    ord= factorback(fa,NULL);
  }
  else
    fa = Z_factor(ord);
  ex = gel(fa,2);
  fa = gel(fa,1);
  l = lg(fa);
  ginv = Fp_inv(g,p);
  v = cgetg(l, t_VEC);
  for (i=1; i<l; i++)
  {
    q = gel(fa,i);
    e = itos(gel(ex,i));
    if (DEBUGLEVEL>5)
      fprintferr("Pohlig-Hellman: DL mod %Z^%ld\n",q,e);
    qj = new_chunk(e+1); gel(qj,0) = gen_1;
    for (j=1; j<=e; j++) gel(qj,j) = mulii(gel(qj,j-1), q);
    t0 = diviiexact(ord, gel(qj,e));
    a0 = Fp_pow(a, t0, p);
    ginv0 = Fp_pow(ginv, t0, p); /* order q^e */
    g_q = Fp_pow(g, diviiexact(ord,q), p); /* order q */
    n_q = gen_0;
    for (j=0; j<e; j++)
    {
      b = modii(mulii(a0, Fp_pow(ginv0, n_q, p)), p);
      b = Fp_pow(b, gel(qj,e-1-j), p);
      b = Fp_shanks(b, g_q, p, q);
      n_q = addii(n_q, mulii(b, gel(qj,j)));
    }
    gel(v,i) = gmodulo(n_q, gel(qj,e));
  }
  return gerepileuptoint(av, lift(chinese1(v)));
}

/* discrete log in Fq for a in Fp^*, g primitive root in Fq^* */
static GEN
ff_PHlog_Fp(GEN a, GEN g, GEN T, GEN p)
{
  pari_sp av = avma;
  GEN q,n_q,ord,ordp;

  if (gcmp1(a)) { avma = av; return gen_0; }
  if (equaliu(p,2)) {
    if (!signe(a)) pari_err(talker,"a not invertible in ff_PHlog_Fp");
    avma = av; return gen_0;
  }
  ordp = subis(p, 1);
  ord = T? subis(powiu(p,degpol(T)), 1): p;
  if (equalii(a, ordp)) /* -1 */
    return gerepileuptoint(av, shifti(ord,-1));

  if (!T) q = NULL;
  else
  { /* we want < g > = Fp^* */
    q = diviiexact(ord,ordp);
    g = FpXQ_pow(g,q,T,p);
    if (typ(g) == t_POL) g = constant_term(g);
  }
  n_q = Fp_PHlog(a,g,p,NULL);
  if (q) n_q = mulii(q, n_q);
  return gerepileuptoint(av, n_q);
}

/* smallest n >= 0 such that g0^n=x modulo pr, assume g0 reduced
 * N = order of g0 is prime (and != p) */
static GEN
ffshanks(GEN x, GEN g0, GEN N, GEN T, GEN p)
{
  pari_sp av = avma, av1, lim;
  long lbaby,i,k;
  GEN p1,smalltable,giant,perm,v,g0inv;

  if (!degpol(x)) x = constant_term(x);
  if (typ(x) == t_INT)
  {
    if (!gcmp1(modii(p,N))) return gen_0;
    /* g0 in Fp^*, order N | (p-1) */
    if (typ(g0) == t_POL) g0 = constant_term(g0);
    return Fp_PHlog(x,g0,p,N);
  }

  p1 = sqrti(N);
  if (cmpiu(p1,LGBITS) >= 0) pari_err(talker,"module too large in ffshanks");
  lbaby = itos(p1)+1; smalltable = cgetg(lbaby+1,t_VEC);
  g0inv = Fq_inv(g0,T,p);
  p1 = x;

  for (i=1;;i++)
  {
    if (gcmp1(p1)) { avma = av; return stoi(i-1); }
    gel(smalltable,i) = p1; if (i==lbaby) break;
    av1 = avma;
    p1 = gerepileupto(av1, FpXQ_mul(p1,g0inv, T,p));
  }
  giant = FpXQ_div(x,p1,T,p);
  perm = gen_sort(smalltable, cmp_IND | cmp_C, cmp_pol);
  smalltable = vecpermute(smalltable, perm);
  p1 = giant;

  av1 = avma; lim=stack_lim(av1,2);
  for (k=1;;k++)
  {
    i = tablesearch(smalltable,p1,cmp_pol);
    if (i)
    {
      v = addis(mulss(lbaby-1,k), perm[i]);
      return gerepileuptoint(av, addsi(-1,v));
    }
    p1 = FpXQ_mul(p1, giant, T,p);

    if (low_stack(lim, stack_lim(av1,2)))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"ffshanks");
      p1 = gerepileupto(av1, p1);
    }
  }
}

/* same in Fp[X]/T */
GEN
ff_PHlog(GEN a, GEN g, GEN T, GEN p)
{
  pari_sp av = avma;
  GEN v,t0,a0,b,q,g_q,n_q,ginv0,qj,ginv,ord,fa,ex;
  long e,i,j,l;

  if (typ(a) == t_INT)
    return gerepileuptoint(av, ff_PHlog_Fp(a,g,T,p));
  /* f > 1 ==> T != NULL */
  ord = subis(powiu(p, degpol(T)), 1);
  fa = factor(ord);
  ex = gel(fa,2);
  fa = gel(fa,1);
  l = lg(fa);
  ginv = Fq_inv(g,T,p);
  v = cgetg(l, t_VEC);
  for (i=1; i<l; i++)
  {
    q = gel(fa,i);
    e = itos(gel(ex,i));
    if (DEBUGLEVEL>5) fprintferr("nf_Pohlig-Hellman: DL mod %Z^%ld\n",q,e);
    qj = new_chunk(e+1); gel(qj,0) = gen_1;
    for (j=1; j<=e; j++) gel(qj,j) = mulii(gel(qj,j-1), q);
    t0 = diviiexact(ord, gel(qj,e));
    a0 = FpXQ_pow(a, t0, T,p);
    ginv0 = FpXQ_pow(ginv, t0, T,p); /* order q^e */
    g_q = FpXQ_pow(g, diviiexact(ord,q), T,p); /* order q */
    n_q = gen_0;
    for (j=0; j<e; j++)
    {
      b = FpXQ_mul(a0, FpXQ_pow(ginv0, n_q, T,p), T,p);
      b = FpXQ_pow(b, gel(qj,e-1-j), T,p);
      b = ffshanks(b, g_q, q, T,p);
      n_q = addii(n_q, mulii(b, gel(qj,j)));
    }
    gel(v,i) = gmodulo(n_q, gel(qj,e));
  }
  return gerepileuptoint(av, lift(chinese1(v)));
}

/* same in nf.zk / pr */
GEN
nf_PHlog(GEN nf, GEN a, GEN g, GEN pr)
{
  pari_sp av = avma;
  GEN T,p, modpr = nf_to_ff_init(nf, &pr, &T, &p);
  GEN A = nf_to_ff(nf,a,modpr);
  GEN G = nf_to_ff(nf,g,modpr);
  return gerepileuptoint(av, ff_PHlog(A,G,T,p));
}

GEN
znlog(GEN x, GEN g0)
{
  pari_sp av = avma;
  GEN p = gel(g0,1);
  if (typ(g0) != t_INTMOD) pari_err(talker,"not an element of (Z/pZ)* in znlog");
  return gerepileuptoint(av, Fp_PHlog(Rg_to_Fp(x,p),gel(g0,2),p,NULL));
}

GEN
dethnf(GEN mat)
{
  long i,l = lg(mat);
  pari_sp av;
  GEN s;

  if (l<3) return l<2? gen_1: icopy(gcoeff(mat,1,1));
  av = avma; s = gcoeff(mat,1,1);
  for (i=2; i<l; i++) s = gmul(s,gcoeff(mat,i,i));
  return av==avma? gcopy(s): gerepileupto(av,s);
}

GEN
dethnf_i(GEN mat)
{
  pari_sp av;
  long i,l = lg(mat);
  GEN s;

  if (l<3) return l<2? gen_1: icopy(gcoeff(mat,1,1));
  av = avma; s = gcoeff(mat,1,1);
  for (i=2; i<l; i++) s = mulii(s,gcoeff(mat,i,i));
  return gerepileuptoint(av,s);
}

/* as above with cyc = diagonal(Smith Normal Form) */
GEN
detcyc(GEN cyc, long *L)
{
  pari_sp av = avma;
  long i, l = lg(cyc);
  GEN s = gel(cyc,1);

  if (l == 1) { *L = 1; return gen_1; }
  for (i=2; i<l; i++)
  {
    GEN t = gel(cyc,i);
    if (is_pm1(t)) break;
    s = mulii(s, t);
  }
  *L = i; return i <= 2? icopy(s): gerepileuptoint(av,s);
}

/* (U,V) = 1. Return y = x mod U, = 1 mod V (uv: u + v = 1, u in U, v in V).
 * mv = multiplication table by v */
static GEN
makeprimetoideal(GEN UV, GEN u,GEN mv, GEN x)
{
  return nfreducemodideal_i(gadd(u, gmul(mv,x)), UV);
}

static GEN
makeprimetoidealvec(GEN nf, GEN UV, GEN u,GEN v, GEN gen)
{
  long i, lx = lg(gen);
  GEN y = cgetg(lx,t_VEC), mv = eltmul_get_table(nf, v);
  for (i=1; i<lx; i++) gel(y,i) = makeprimetoideal(UV,u,mv, gel(gen,i));
  return y;
}

GEN
FpXQ_gener(GEN T, GEN p)
{
  long i,j, k, vT = varn(T), f = degpol(T);
  GEN g, L, pf_1 = subis(powiu(p, f), 1);
  pari_sp av0 = avma, av;

  L = (GEN)Z_factor(pf_1)[1];
  k = lg(L)-1;

  for (i=1; i<=k; i++) gel(L,i) = diviiexact(pf_1, gel(L,i));
  for (av = avma;; avma = av)
  {
    g = FpX_rand(f, vT, p);
    if (degpol(g) < 1) continue;
    for (j=1; j<=k; j++)
      if (gcmp1(FpXQ_pow(g, gel(L,j), T, p))) break;
    if (j > k) return gerepilecopy(av0, g);
  }
}

/* Given an ideal pr^ep, and an integral ideal x (in HNF form) compute a list
 * of vectors,corresponding to the abelian groups (O_K/pr)^*, and
 * 1 + pr^i/ 1 + pr^min(2i, ep), i = 1,...
 * Each vector has 5 components as follows :
 * [[cyc],[g],[g'],[sign],U.X^-1].
 * cyc   = type as abelian group
 * g, g' = generators. (g',x) = 1, not necessarily so for g
 * sign  = vector of the sign(g') at arch.
 * If x = NULL, the original ideal was a prime power */
static GEN
zprimestar(GEN nf, GEN pr, GEN ep, GEN x, GEN arch)
{
  pari_sp av = avma;
  long a, b, e = itos(ep), f = itos(gel(pr,4));
  GEN p = gel(pr,1), list, g, g0, y, u,v, prh, prb, pre;

  if(DEBUGLEVEL>3) fprintferr("treating pr^%ld, pr = %Z\n",e,pr);
  if (f == 1)
    g = gscalcol_i(gener_Fp(p), degpol(nf[1]));
  else
  {
    GEN T, modpr = zk_to_ff_init(nf, &pr, &T, &p);
    g = ff_to_nf(FpXQ_gener(T,p), modpr);
    g = poltobasis(nf, g);
  }
  /* g generates  (Z_K / pr)^* */
  prh = prime_to_ideal(nf,pr);
  pre = (e==1)? prh: idealpow(nf,pr,ep);
  g0 = g;
  u = v = NULL; /* gcc -Wall */
  if (x)
  {
    GEN uv = idealaddtoone(nf,pre, idealdivpowprime(nf,x,pr,ep));
    u = gel(uv,1);
    v = gel(uv,2); v = eltmul_get_table(nf, v);
    g0 = makeprimetoideal(x,u,v,g);
  }

  list = cget1(e+1, t_VEC);
  y = cgetg(6,t_VEC); appendL(list, y);
  gel(y,1) = mkvec(addis(powiu(p,f), -1));
  gel(y,2) = mkvec(g);
  gel(y,3) = mkvec(g0);
  gel(y,4) = mkvec(zsigne(nf,g0,arch));
  gel(y,5) = gen_1;
  prb = prh;
  for (a = b = 1; a < e; a = b)
  {
    GEN pra = prb, gen, z, s, U;
    long i, l;

    b <<= 1;
    /* compute 1 + pr^a / 1 + pr^b, 2a <= b */
    if(DEBUGLEVEL>3) fprintferr("  treating a = %ld, b = %ld\n",a,b);
    prb = (b >= e)? pre: idealpows(nf,pr,b);
    z = zidealij(pra, prb, &U);
    gen = shallowcopy(gel(z,2));
    l = lg(gen); s = cgetg(l, t_VEC);
    if(DEBUGLEVEL>3) fprintferr("zidealij done\n");
    for (i = 1; i < l; i++)
    {
      if (x) gel(gen,i) = makeprimetoideal(x,u,v,gel(gen,i));
      gel(s,i) = zsigne(nf,gel(gen,i),arch);
    }
    y = cgetg(6,t_VEC); appendL(list, y);
    y[1] = z[1];
    y[2] = z[2];
    gel(y,3) = gen;
    gel(y,4) = s;
    gel(y,5) = U;
  }
  return gerepilecopy(av, list);
}

/* increment y, which runs through [-d,d]^k. Return 0 when done. */
static int
increment(GEN y, long k, long d)
{
  long i = 0, j;
  do
  {
    if (++i > k) return 0;
    y[i]++;
  } while (y[i] > d);
  for (j = 1; j < i; j++) y[j] = -d;
  return 1;
}

GEN
archstar_full_rk(GEN x, GEN bas, GEN v, GEN gen)
{
  long i, r, lgmat, N = lg(bas)-1, nba = lg(v[1]) - 1;
  GEN lambda = cgetg(N+1, t_VECSMALL), mat = cgetg(nba+1,t_MAT);

  lgmat = lg(v); setlg(mat, lgmat+1);
  for (i = 1; i < lgmat; i++) mat[i] = v[i];
  for (     ; i <= nba; i++)  gel(mat,i) = cgetg(nba+1, t_VECSMALL);

  if (x) { x = lllint_ip(x,4); bas = gmul(bas, x); }

  for (r=1;; r++)
  { /* reset */
    (void)vec_setconst(lambda, (GEN)0);
    if (!x) lambda[1] = r;
    while (increment(lambda, N, r))
    {
      pari_sp av1 = avma;
      GEN a = RgM_zc_mul(bas, lambda), c = gel(mat,lgmat);
      for (i = 1; i <= nba; i++)
      {
        GEN t = x? gadd(gel(a,i), gen_1): gel(a,i);
        c[i] = (gsigne(t) < 0)? 1: 0;
      }
      avma = av1; if (Flm_deplin(mat, 2)) continue;

      /* c independant of previous sign vectors */
      if (!x) a = zc_to_ZC(lambda);
      else
      {
        a = ZM_zc_mul(x, lambda);
        gel(a,1) = addis(gel(a,1), 1);
      }
      gel(gen,lgmat) = a;
      if (lgmat++ == nba) return Flm_to_ZM( Flm_inv(mat,2) ); /* full rank */
      setlg(mat,lgmat+1);
    }
  }
}

/* x integral ideal, compute elements in 1+x (in x, if x = zk) whose sign
 * matrix is invertible */
GEN
zarchstar(GEN nf, GEN x, GEN archp)
{
  long i, nba;
  pari_sp av;
  GEN p1, y, bas, gen, mat, gZ, v;

  archp = arch_to_perm(archp);
  nba = lg(archp) - 1;
  y = cgetg(4,t_VEC);
  if (!nba)
  {
    gel(y,1) = cgetg(1,t_VEC);
    gel(y,2) = cgetg(1,t_VEC);
    gel(y,3) = cgetg(1,t_MAT); return y;
  }
  p1 = cgetg(nba+1,t_VEC); for (i=1; i<=nba; i++) gel(p1,i) = gen_2;
  gel(y,1) = p1; av = avma;
  if (gcmp1(gcoeff(x,1,1))) x = NULL; /* x = O_K */
  gZ = x? subsi(1, gcoeff(x,1,1)): gen_m1; /* gZ << 0, gZ = 1 mod x */
  if (nba == 1)
  {
    gel(y,2) = mkvec(gZ);
    gel(y,3) = gscalmat(gen_1,1); return y;
  }
  bas = gmael(nf,5,1);
  if (lg(bas[1]) > lg(archp)) bas = rowpermute(bas, archp);
  gen = cgetg(nba+1,t_VEC);
  v = mkmat( const_vecsmall(nba, 1) );
  gel(gen,1) = gZ;

  mat = archstar_full_rk(x, bas, v, gen);
  gerepileall(av,2,&gen,&mat);

  gel(y,2) = gen;
  gel(y,3) = mat; return y;
}

static GEN
zlog_pk(GEN nf, GEN a0, GEN y, GEN pr, GEN prk, GEN list, GEN *psigne)
{
  GEN a = a0, L, e, cyc, gen, s, U;
  long i,j, llist = lg(list)-1;

  for (j = 1; j <= llist; j++)
  {
    L = gel(list,j);
    cyc = gel(L,1);
    gen = gel(L,2);
    s   = gel(L,4);
    U   = gel(L,5);
    if (j == 1)
      e = mkcol( nf_PHlog(nf, a, gel(gen,1), pr) );
    else if (typ(a) == t_INT)
      e = gmul(subis(a, 1), gel(U,1));
    else
    { /* t_COL */
      GEN t = gel(a,1);
      gel(a,1) = addsi(-1, gel(a,1)); /* a -= 1 */
      e = gmul(U, a);
      gel(a,1) = t; /* restore */
    }
    /* here lg(e) == lg(cyc) */
    for (i = 1; i < lg(cyc); i++)
    {
      GEN t = modii(negi(gel(e,i)), gel(cyc,i));
      gel(++y,0) = negi(t); if (!signe(t)) continue;

      if (mod2(t)) *psigne = *psigne? gadd(*psigne, gel(s,i)): gel(s,i);
      if (j != llist) a = elt_mulpow_modideal(nf, a, gel(gen,i), t, prk);
    }
  }
  return y;
}

static void
zlog_add_sign(GEN y0, GEN sgn, GEN lists)
{
  GEN y, s;
  long i;
  if (!sgn) return;
  y = y0 + lg(y0);
  s = gmul(gmael(lists, lg(lists)-1, 3), sgn);
  for (i = lg(s)-1; i > 0; i--) gel(--y,0) = mpodd(gel(s,i))? gen_1: gen_0;
}

static GEN
famat_zlog(GEN nf, GEN g, GEN e, GEN sgn, GEN bid)
{
  GEN vp = gmael(bid, 3,1), ep = gmael(bid, 3,2), arch = gmael(bid,1,2);
  GEN cyc = gmael(bid,2,2), lists = gel(bid,4), U = gel(bid,5);
  GEN y0, x, y, EX = gel(cyc,1);
  long i, l;

  y0 = y = cgetg(lg(U), t_COL);
  if (!sgn) sgn = zsigne(nf, to_famat(g,e), arch);
  l = lg(vp);
  for (i=1; i < l; i++)
  {
    GEN pr = gel(vp,i), prk;
    prk = (l==2)? gmael(bid,1,1): idealpow(nf, pr, gel(ep,i));
    /* FIXME: FIX group exponent (should be mod prk, not f !) */
    x = famat_makecoprime(nf, g, e, pr, prk, EX);
    y = zlog_pk(nf, x, y, pr, prk, gel(lists,i), &sgn);
  }
  zlog_add_sign(y0, sgn, lists);
  return y0;
}

static GEN
get_index(GEN lists)
{
  long t = 0, j, k, l = lg(lists)-1;
  GEN L, ind = cgetg(l+1, t_VECSMALL);

  for (k = 1; k < l; k++)
  {
    L = gel(lists,k);
    ind[k] = t;
    for (j=1; j<lg(L); j++) t += lg(gmael(L,j,1)) - 1;
  }
  /* for arch. components */
  ind[k] = t; return ind;
}

static void
init_zlog(zlog_S *S, long n, GEN P, GEN e, GEN arch, GEN lists, GEN U)
{
  S->n = n;
  S->U = U;
  S->P = P;
  S->e = e;
  S->archp = arch_to_perm(arch);
  S->lists = lists;
  S->ind = get_index(lists);
}
void
init_zlog_bid(zlog_S *S, GEN bid)
{
  GEN ideal = gel(bid,1), fa = gel(bid,3), fa2 = gel(bid,4), U = gel(bid,5);
  GEN arch = (typ(ideal)==t_VEC && lg(ideal)==3)? gel(ideal,2): NULL;
  init_zlog(S, lg(U)-1, gel(fa,1), gel(fa,2), arch, fa2, U);
}

/* Return decomposition of a on the S->nf successive generators contained in
 * S->lists. If index !=0, do the computation for the corresponding prime
 * ideal and set to 0 the other components. */
static GEN
zlog_ind(GEN nf, GEN a, zlog_S *S, GEN sgn, long index)
{
  GEN y0 = zerocol(S->n), y, list, pr, prk;
  pari_sp av = avma;
  long i, k, kmin, kmax;

  if (typ(a) != t_INT) a = algtobasis_i(nf,a);
  if (DEBUGLEVEL>3)
  {
    fprintferr("entering zlog, "); flusherr();
    if (DEBUGLEVEL>5) fprintferr("with a = %Z\n",a);
  }
  if (index)
  {
    kmin = kmax = index;
    y = y0 + S->ind[index];
  }
  else
  {
    kmin = 1; kmax = lg(S->P)-1;
    y = y0;
  }
  if (!sgn) sgn = zsigne(nf, a, S->archp);
  for (k = kmin; k <= kmax; k++)
  {
    list= (GEN)S->lists[k];
    pr  = (GEN)S->P[k];
    prk = idealpow(nf, pr, (GEN)S->e[k]);
    y = zlog_pk(nf, a, y, pr, prk, list, &sgn);
  }
  zlog_add_sign(y0, sgn, S->lists);
  if (DEBUGLEVEL>3) { fprintferr("leaving\n"); flusherr(); }
  avma = av;
  for (i=1; i <= S->n; i++) gel(y0,i) = icopy(gel(y0,i));
  return y0;
}
/* sgn = sign(a, S->arch) or NULL if unknown */
GEN
zlog(GEN nf, GEN a, GEN sgn, zlog_S *S) { return zlog_ind(nf, a, S, sgn, 0); }

/* Log on bid.gen of generators of P_{1,I pr^{e-1}} / P_{1,I pr^e} (I,pr) = 1,
 * defined implicitly via CRT. 'index' is the index of pr in modulus
 * factorization */
GEN
log_gen_pr(zlog_S *S, long index, GEN nf, long e)
{
  long i, l, yind = S->ind[index];
  GEN y, A, L, L2 = (GEN)S->lists[index];

  if (e == 1)
  {
    L = gel(L2,1);
    y = zerocol(S->n); gel(y, yind+1) = gen_1;
    zlog_add_sign(y, gmael(L,4,1), S->lists);
    A = mkmat(y);
  }
  else
  {
    GEN pr = (GEN)S->P[index], prk, g;

    if (e == 2)
      L = gel(L2,2);
    else
      L = zidealij(idealpows(nf,pr,e-1), idealpows(nf,pr,e), NULL);
    g = gel(L,2);
    l = lg(g);
    A = cgetg(l, t_MAT);
    prk = idealpow(nf, pr, (GEN)S->e[index]);
    for (i = 1; i < l; i++)
    {
      GEN G = gel(g,i), sgn = NULL; /* positive at f_oo */
      y = zerocol(S->n);
      (void)zlog_pk(nf, G, y + yind, pr, prk, L2, &sgn);
      zlog_add_sign(y, sgn, S->lists);
      gel(A,i) = y;
    }
  }
  return gmul(S->U, A);
}
/* Log on bid.gen of generator of P_{1,f} / P_{1,f v[index]}
 * v = vector of r1 real places */
GEN
log_gen_arch(zlog_S *S, long index)
{
  GEN y = zerocol(S->n);
  zlog_add_sign(y, col_ei(lg(S->archp)-1, index), S->lists);
  return gmul(S->U, y);
}

static GEN
compute_gen(GEN nf, GEN u1, GEN gen, GEN bid)
{
  long i, c = lg(u1);
  GEN L = cgetg(c,t_VEC);
  for (i=1; i<c; i++)
    gel(L,i) = famat_to_nf_modidele(nf, gen, gel(u1,i), bid);
  return L;
}
static void
add_clgp(GEN nf, GEN u1, GEN cyc, GEN gen, GEN bid)
{
  GEN c = cgetg(u1? 4: 3, t_VEC);
  long L;
  gel(bid,2) = c;
  gel(c,1) = detcyc(cyc, &L);
  gel(c,2) = cyc;
  if (u1) gel(c,3) = (u1 == gen_1? gen: compute_gen(nf, u1, gen, bid));
}

/* Compute [[ideal,arch], [h,[cyc],[gen]], idealfact, [liste], U]
   gen not computed unless add_gen != 0 */
GEN
Idealstar(GEN nf, GEN ideal,long add_gen)
{
  pari_sp av = avma;
  long i, j, k, nbp, R1, nbgen;
  GEN t, y, cyc, U, u1 = NULL, fa, lists, x, arch, archp, E, P, sarch, gen;

  nf = checknf(nf);
  R1 = nf_get_r1(nf);
  if (typ(ideal) == t_VEC && lg(ideal) == 3)
  {
    arch = gel(ideal,2); ideal = gel(ideal,1);
    i = typ(arch);
    if (!is_vec_t(i) || lg(arch) != R1+1)
      pari_err(talker,"incorrect archimedean component in Idealstar");
    archp = arch_to_perm(arch);
  }
  else
  {
    arch = zerovec(R1);
    archp = cgetg(1, t_VECSMALL);
  }
  x = idealhermite_aux(nf, ideal);
  if (lg(x) == 1 || !gcmp1(denom(gcoeff(x,1,1))))
    pari_err(talker,"Idealstar needs an integral non-zero ideal: %Z",x);
  fa = idealfactor(nf, ideal);
  P = gel(fa,1);
  E = gel(fa,2); nbp = lg(P)-1;
  lists = cgetg(nbp+2,t_VEC);

  gen = cgetg(1,t_VEC);
  t = (nbp==1)? (GEN)NULL: x;
  for (i=1; i<=nbp; i++)
  {
    GEN L = zprimestar(nf, gel(P,i), gel(E,i), t, archp);
    gel(lists,i) = L;
    for (j=1; j<lg(L); j++) gen = shallowconcat(gen,gmael(L,j,3));
  }
  sarch = zarchstar(nf, x, archp);
  gel(lists,i) = sarch;
  gen = shallowconcat(gen, gel(sarch,2));

  nbgen = lg(gen)-1;
  if (nbp)
  {
    GEN h = cgetg(nbgen+1,t_MAT);
    long cp = 0;
    zlog_S S; init_zlog(&S, nbgen, P, E, archp, lists, NULL);
    for (i=1; i<=nbp; i++)
    {
      GEN L2 = gel(lists,i);
      for (j=1; j < lg(L2); j++)
      {
        GEN L = gel(L2,j), F = gel(L,1), G = gel(L,3);
        for (k=1; k<lg(G); k++)
        { /* log(g^f) mod idele */
          GEN g = gel(G,k), f = gel(F,k), a = element_powmodideal(nf,g,f,x);
          GEN sgn = mpodd(f)? zsigne(nf, g, S.archp): zerocol(lg(S.archp)-1);
          gel(h,++cp) = gneg(zlog_ind(nf, a, &S, sgn, i));
          coeff(h,cp,cp) = F[k];
        }
      }
    }
    for (j=1; j<lg(archp); j++)
    {
      gel(h,++cp) = zerocol(nbgen);
      gcoeff(h,cp,cp) = gen_2;
    }
    /* assert(cp == nbgen) */
    h = hnfall_i(h,NULL,0);
    cyc = smithrel(h, &U, add_gen? &u1: NULL);
  }
  else
  {
    cyc = cgetg(nbgen+1, t_VEC);
    for (j=1; j<=nbgen; j++) gel(cyc,j) = gen_2;
    U = matid(nbgen);
    if (add_gen) u1 = gen_1;
  }

  y = cgetg(6,t_VEC);
  gel(y,1) = mkvec2(x, arch);
  gel(y,3) = fa;
  gel(y,4) = lists;
  gel(y,5) = U;
  add_clgp(nf, u1, cyc, gen, y);
  return gerepilecopy(av, y);
}

GEN
zidealstarinitgen(GEN nf, GEN ideal)
{ return Idealstar(nf,ideal,1); }
GEN
zidealstarinit(GEN nf, GEN ideal)
{ return Idealstar(nf,ideal,0); }
/* FIXME: obsolete */
GEN
zidealstar(GEN nf, GEN ideal)
{
  pari_sp av = avma;
  GEN y = Idealstar(nf,ideal,1);
  return gerepilecopy(av,gel(y,2));
}

GEN
idealstar0(GEN nf, GEN ideal,long flag)
{
  switch(flag)
  {
    case 0: return zidealstar(nf,ideal);
    case 1: return zidealstarinit(nf,ideal);
    case 2: return zidealstarinitgen(nf,ideal);
    default: pari_err(flagerr,"idealstar");
  }
  return NULL; /* not reached */
}

void
check_nfelt(GEN x, GEN *den)
{
  long l = lg(x), i;
  GEN t, d = NULL;
  if (typ(x) != t_COL) pari_err(talker,"%Z not a nfelt", x);
  for (i=1; i<l; i++)
  {
    t = gel(x,i);
    switch (typ(t))
    {
      case t_INT: break;
      case t_FRAC:
        if (!d) d = gel(t,2); else d = lcmii(d, gel(t,2));
        break;
      default: pari_err(talker,"%Z not a nfelt", x);
    }
  }
  *den = d;
}

GEN
vecmodii(GEN a, GEN b)
{
  long i, l = lg(a);
  GEN c = cgetg(l, typ(a));
  for (i = 1; i < l; i++) gel(c,i) = modii(gel(a,i), gel(b,i));
  return c;
}

/* Given x (not necessarily integral), and bid as output by zidealstarinit,
 * compute the vector of components on the generators bid[2].
 * Assume (x,bid) = 1 and sgn is either NULL or zsigne(x, bid) */
GEN
zideallog_sgn(GEN nf, GEN x, GEN sgn, GEN bid)
{
  pari_sp av;
  long N, lcyc;
  GEN den, cyc, y;
  int ok = 0;

  nf = checknf(nf); checkbid(bid);
  cyc = gmael(bid,2,2);
  lcyc = lg(cyc); if (lcyc == 1) return cgetg(1, t_COL);
  av = avma;
  N = degpol(nf[1]);
  switch(typ(x))
  {
    case t_INT: case t_FRAC:
      ok = 1; den = denom(x);
      break;
    case t_POLMOD: case t_POL:
      x = algtobasis(nf,x); break;
    case t_COL: break;
    case t_MAT:
      if (lg(x) == 1) return zerocol(lcyc-1);
      y = famat_zlog(nf, gel(x,1), gel(x,2), sgn, bid);
      goto END;

    default: pari_err(talker,"not an element in zideallog");
  }
  if (!ok)
  {
    if (lg(x) != N+1) pari_err(talker,"not an element in zideallog");
    check_nfelt(x, &den);
  }
  if (den)
  {
    GEN g = mkcol2(Q_muli_to_int(x,den), den);
    GEN e = mkcol2(gen_1, gen_m1);
    y = famat_zlog(nf, g, e, sgn, bid);
  }
  else
  {
    zlog_S S; init_zlog_bid(&S, bid);
    y = zlog(nf, x, sgn, &S);
  }
END:
  y = gmul(gel(bid,5), y);
  return gerepileupto(av, vecmodii(y, cyc));
}
GEN
zideallog(GEN nf, GEN x, GEN bid) { return zideallog_sgn(nf, x, NULL, bid); }

/*************************************************************************/
/**									**/
/**		  JOIN BID STRUCTURES, IDEAL LISTS                      **/
/**									**/
/*************************************************************************/

/* bid1, bid2: for coprime modules m1 and m2 (without arch. part).
 * Output: bid [[m1 m2,arch],[h,[cyc],[gen]],idealfact,[liste],U] for m1 m2 */
static GEN
join_bid(GEN nf, GEN bid1, GEN bid2)
{
  pari_sp av = avma;
  long i, nbgen, lx, lx1,lx2, l1,l2;
  GEN I1,I2, f1,f2, G1,G2, fa1,fa2, lists1,lists2, cyc1,cyc2;
  GEN lists, fa, U, cyc, y, u1 = NULL, x, gen;

  f1 = gel(bid1,1); I1 = gel(f1,1);
  f2 = gel(bid2,1); I2 = gel(f2,1);
  if (gcmp1(gcoeff(I1,1,1))) return bid2; /* frequent trivial case */
  G1 = gel(bid1,2); G2 = gel(bid2,2);
  fa1= gel(bid1,3); fa2= gel(bid2,3); x = idealmul(nf, I1,I2);
  fa = concat_factor(fa1, fa2);

  lists1 = gel(bid1,4); lx1 = lg(lists1);
  lists2 = gel(bid2,4); lx2 = lg(lists2);
  /* concat (lists1 - last elt [zarchstar]) + lists2 */
  lx = lx1+lx2-2; lists = cgetg(lx,t_VEC);
  for (i=1; i<lx1-1; i++) lists[i] = lists1[i];
  for (   ; i<lx; i++)    lists[i] = lists2[i-lx1+2];

  cyc1 = gel(G1,2); l1 = lg(cyc1);
  cyc2 = gel(G2,2); l2 = lg(cyc2);
  gen = (lg(G1)>3 && lg(G2)>3)? gen_1: NULL;
  nbgen = l1+l2-2;
  cyc = smithrel(diagonal_i(shallowconcat(cyc1,cyc2)),
                 &U, gen? &u1: NULL);
  if (nbgen) {
    GEN U1 = gel(bid1,5), U2 = gel(bid2,5);
    U1 = l1 == 1? zeromat(nbgen,lg(U1)-1): gmul(vecslice(U, 1, l1-1),   U1);
    U2 = l2 == 1? zeromat(nbgen,lg(U2)-1): gmul(vecslice(U, l1, nbgen), U2);
    U = shallowconcat(U1, U2);
  }
  else
    U = zeromat(0, lx-2);

  if (gen)
  {
    GEN u, v, uv = idealaddtoone(nf,gel(f1,1),gel(f2,1));
    u = gel(uv,1);
    v = gel(uv,2);
    gen = shallowconcat(makeprimetoidealvec(nf,x,u,v, gel(G1,3)),
                   makeprimetoidealvec(nf,x,v,u, gel(G2,3)));
  }
  y = cgetg(6,t_VEC);
  gel(y,1) = mkvec2(x, gel(f1,2));
  gel(y,3) = fa;
  gel(y,4) = lists;
  gel(y,5) = U;
  add_clgp(nf, u1, cyc, gen, y);
  return gerepilecopy(av,y);
}

/* bid1 = for module m1 (without arch. part), arch = archimedean part.
 * Output: bid [[m1,arch],[h,[cyc],[gen]],idealfact,[liste],U] for m1.arch */
static GEN
join_bid_arch(GEN nf, GEN bid1, GEN arch)
{
  pari_sp av = avma;
  long i, lx1;
  GEN f1, G1, fa1, lists1, U;
  GEN lists, cyc, y, u1 = NULL, x, sarch, gen;

  checkbid(bid1);
  f1 = gel(bid1,1); G1 = gel(bid1,2); fa1 = gel(bid1,3);
  x = gel(f1,1);
  sarch = zarchstar(nf, x, arch);
  lists1 = gel(bid1,4); lx1 = lg(lists1);
  lists = cgetg(lx1,t_VEC);
  for (i=1; i<lx1-1; i++) lists[i] = lists1[i];
  gel(lists,i) = sarch;

  gen = (lg(G1)>3)? gen_1: NULL;
  cyc = diagonal_i(shallowconcat(gel(G1,2), gel(sarch,1)));
  cyc = smithrel(cyc, &U, gen? &u1: NULL);
  if (gen) gen = shallowconcat(gel(G1,3), gel(sarch,2));
  y = cgetg(6,t_VEC);
  gel(y,1) = mkvec2(x, arch);
  gel(y,3) = fa1;
  gel(y,4) = lists;
  gel(y,5) = U;
  add_clgp(nf, u1, cyc, gen, y);
  return gerepilecopy(av,y);
}

#if 0 /* could be useful somewhere else */
/* z <- ( z | f(v[i])_{i=1..#v} )*/
void
concatmap(GEN *pz, GEN v, GEN (*f)(void*,GEN), void *data)
{
  long i, nz, lv = lg(v);
  GEN z, Z, Zt;
  if (lv == 1) return;
  z = *pz; nz = lg(z)-1;
  Z = cgetg(lv + nz, typ(z));
  for (i = 1; i <=nz; i++) Z[i] = z[i];
  Zt = Z + nz;
  for (i = 1; i < lv; i++) gel(Zt,i) = f(data, gel(v,i));
  *pz = Z;
}
#endif

typedef struct _ideal_data {
  GEN nf, emb, L, pr, prL, arch, sgnU;
} ideal_data;
static void
concat_join(GEN *pz, GEN v, GEN (*f)(ideal_data*,GEN), ideal_data *data)
{
  long i, nz, lv = lg(v);
  GEN z, Z, Zt;
  if (lv == 1) return;
  z = *pz; nz = lg(z)-1;
  Z = cgetg(lv + nz, typ(z));
  for (i = 1; i <=nz; i++) Z[i] = z[i];
  Zt = Z + nz;
  for (i = 1; i < lv; i++) gel(Zt,i) = f(data, gel(v,i));
  *pz = Z;
}
static GEN
join_idealinit(ideal_data *D, GEN x) {
  return join_bid(D->nf, x, D->prL);
}
static GEN
join_ideal(ideal_data *D, GEN x) {
  return idealmulpowprime(D->nf, x, D->pr, D->L);
}
static GEN
join_unit(ideal_data *D, GEN x) {
  return mkvec2(join_idealinit(D, gel(x,1)), vconcat(gel(x,2), D->emb));
}

/* compute matrix of zlogs of units */
GEN
zlog_units(GEN nf, GEN U, GEN sgnU, GEN bid)
{
  long j, l = lg(U);
  GEN m = cgetg(l, t_MAT);
  zlog_S S; init_zlog_bid(&S, bid);
  for (j = 1; j < l; j++)
    gel(m,j) = zlog(nf, gel(U,j), vecpermute(gel(sgnU,j), S.archp), &S);
  return m;
}
/* compute matrix of zlogs of units, assuming S.archp = [] */
GEN
zlog_units_noarch(GEN nf, GEN U, GEN bid)
{
  long j, l = lg(U);
  GEN m = cgetg(l, t_MAT), empty = cgetg(1, t_COL);
  zlog_S S; init_zlog_bid(&S, bid);
  for (j = 1; j < l; j++) gel(m,j) = zlog(nf, gel(U,j), empty, &S);
  return m;
}

/* calcule la matrice des zlog des unites */
static GEN
zlog_unitsarch(GEN sgnU, GEN bid)
{
  GEN U, liste = gel(bid,4), arch = gmael(bid,1,2);
  long i;
  U = gmul(gmael(liste, lg(liste)-1, 3),
           rowpermute(sgnU, arch_to_perm(arch)));
  for (i = 1; i < lg(U); i++) (void)F2V_red_ip(gel(U,i));
  return U;
}

/*  flag &1 : generators, otherwise no
 *  flag &2 : units, otherwise no
 *  flag &4 : ideals in HNF, otherwise bid */
static GEN
Ideallist(GEN bnf, ulong bound, long flag)
{
  const long do_gen = flag & 1, do_units = flag & 2, big_id = !(flag & 4);
  byteptr ptdif = diffptr;
  pari_sp lim, av, av0 = avma;
  long i, j, l;
  GEN nf, z, p, fa, id, U, empty = cgetg(1,t_VEC);
  ideal_data ID;
  GEN (*join_z)(ideal_data*, GEN) =
    do_units? &join_unit
              : (big_id? &join_idealinit: &join_ideal);

  nf = checknf(bnf);
  if ((long)bound <= 0) return empty;
  id = matid(degpol(nf[1]));
  if (big_id) id = Idealstar(nf,id,do_gen);

  /* z[i] will contain all "objects" of norm i. Depending on flag, this means
   * an ideal, a bid, or a couple [bid, log(units)]. Such objects are stored
   * in vectors, computed one primary component at a time; join_z
   * reconstructs the global object */
  z = cgetg(bound+1,t_VEC);
  if (do_units) {
    U = init_units(bnf);
    gel(z,1) = mkvec( mkvec2(id, zlog_units_noarch(nf, U, id)) );
  } else {
    U = NULL; /* -Wall */
    gel(z,1) = mkvec(id);
  }
  for (i=2; i<=(long)bound; i++) gel(z,i) = empty;
  ID.nf = nf;

  p = cgeti(3); p[1] = evalsigne(1) | evallgefint(3);
  av = avma; lim = stack_lim(av,1);
  maxprime_check(bound);
  for (p[2] = 0; (ulong)p[2] <= bound; )
  {
    NEXT_PRIME_VIADIFF(p[2], ptdif);
    if (DEBUGLEVEL>1) { fprintferr("%ld ",p[2]); flusherr(); }
    fa = primedec(nf, p);
    for (j=1; j<lg(fa); j++)
    {
      GEN pr = gel(fa,j), z2;
      ulong q, iQ, Q = itou_or_0(pr_norm(pr));
      if (!Q || Q > bound) break;

      z2 = shallowcopy(z);
      q = Q;
      ID.pr = ID.prL = pr;
      for (l=1; Q <= bound; l++, Q *= q) /* add pr^l */
      {
        ID.L = utoipos(l);
        if (big_id) {
          if (l > 1) ID.prL = idealpow(nf,pr,ID.L);
          ID.prL = Idealstar(nf,ID.prL,do_gen);
          if (do_units) ID.emb = zlog_units_noarch(nf, U, ID.prL);
        }
        for (iQ = Q,i = 1; iQ <= bound; iQ += Q,i++)
          concat_join(&gel(z,iQ), gel(z2,i), join_z, &ID);
      }
    }
    if (low_stack(lim, stack_lim(av,1)))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"Ideallist");
      z = gerepilecopy(av, z);
    }
  }
  if (do_units) for (i = 1; i < lg(z); i++)
  {
    GEN s = gel(z,i);
    long l = lg(s);
    for (j = 1; j < l; j++) {
      GEN v = gel(s,j), bid = gel(v,1);
      gel(v,2) = gmul(gel(bid,5), gel(v,2));
    }
  }
  return gerepilecopy(av0, z);
}
GEN
ideallist0(GEN bnf,long bound, long flag) {
  if (flag<0 || flag>4) pari_err(flagerr,"ideallist");
  return Ideallist(bnf,bound,flag);
}
GEN
ideallistzstar(GEN nf,long bound) { return Ideallist(nf,bound,0); }
GEN
ideallistzstargen(GEN nf,long bound) { return Ideallist(nf,bound,1); }
GEN
ideallistunit(GEN nf,long bound) { return Ideallist(nf,bound,2); }
GEN
ideallistunitgen(GEN nf,long bound) { return Ideallist(nf,bound,3); }
GEN
ideallist(GEN bnf,long bound) { return Ideallist(bnf,bound,4); }

static GEN
join_arch(ideal_data *D, GEN x) {
  return join_bid_arch(D->nf, x, D->arch);
}
static GEN
join_archunit(ideal_data *D, GEN x) {
  GEN bid = join_arch(D, gel(x,1)), U = gel(x,2);
  U = gmul(gel(bid,5), vconcat(U, zlog_unitsarch(D->sgnU, bid)));
  return mkvec2(bid, U);
}

/* L from ideallist, add archimedean part */
GEN
ideallistarch(GEN bnf, GEN L, GEN arch)
{
  pari_sp av;
  long i, j, l = lg(L), lz;
  GEN v, z, V;
  ideal_data ID;
  GEN (*join_z)(ideal_data*, GEN) = &join_arch;

  if (typ(L) != t_VEC) pari_err(typeer, "ideallistarch");
  if (l == 1) return cgetg(1,t_VEC);
  z = gel(L,1);
  if (typ(z) != t_VEC) pari_err(typeer, "ideallistarch");
  z = gel(z,1); /* either a bid or [bid,U] */
  if (lg(z) == 3) { /* the latter: do units */
    if (typ(z) != t_VEC) pari_err(typeer,"ideallistarch");
    join_z = &join_archunit;
    ID.sgnU = zsignunits(bnf, NULL, 1);
  }
  ID.nf = checknf(bnf);
  arch = arch_to_perm(arch);
  av = avma; V = cgetg(l, t_VEC);
  for (i = 1; i < l; i++)
  {
    z = gel(L,i); lz = lg(z);
    gel(V,i) = v = cgetg(lz,t_VEC);
    for (j=1; j<lz; j++) gel(v,j) = join_z(&ID, gel(z,j));
  }
  return gerepilecopy(av,V);
}
