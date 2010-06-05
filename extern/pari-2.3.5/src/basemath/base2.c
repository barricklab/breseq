/* $Id: base2.c 10291 2008-06-10 15:50:00Z kb $

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
/*                       MAXIMAL ORDERS                            */
/*                                                                 */
/*******************************************************************/
#include "pari.h"
#include "paripriv.h"

/* FIXME: backward compatibility. Should use the proper nf_* equivalents */
#define compat_PARTIAL 1
#define compat_ROUND2  2

static void
allbase_check_args(GEN f, long flag, GEN *dx, GEN *ptw)
{
  GEN w = *ptw;

  if (DEBUGLEVEL) (void)timer2();
  if (typ(f)!=t_POL) pari_err(notpoler,"allbase");
  if (degpol(f) <= 0) pari_err(constpoler,"allbase");

  *dx = w? factorback(w, NULL): ZX_disc(f);
  if (!signe(*dx)) pari_err(talker,"reducible polynomial in allbase");
  if (!w) *ptw = auxdecomp(absi(*dx), (flag & nf_PARTIALFACT)? 0: 1);
  if (DEBUGLEVEL) msgtimer("disc. factorisation");
}

/*******************************************************************/
/*                                                                 */
/*                            ROUND 2                              */
/*                                                                 */
/*******************************************************************/
/* companion matrix of unitary polynomial x */
static GEN
companion(GEN x) /* cf assmat */
{
  long i,j,l;
  GEN y;

  l=degpol(x)+1; y=cgetg(l,t_MAT);
  for (j=1; j<l; j++)
  {
    gel(y,j) = cgetg(l,t_COL);
    for (i=1; i<l-1; i++) gcoeff(y,i,j)=(i+1==j)? gen_1: gen_0;
    gcoeff(y,i,j) = gneg(gel(x,j+1));
  }
  return y;
}

/* assume x, y are square integer matrices of same dim. Multiply them */
static GEN
mulmati(GEN x, GEN y)
{
  pari_sp av;
  long n = lg(x),i,j,k;
  GEN z = cgetg(n,t_MAT),p1,p2;

  for (j=1; j<n; j++)
  {
    gel(z,j) = cgetg(n,t_COL);
    for (i=1; i<n; i++)
    {
      p1=gen_0; av=avma;
      for (k=1; k<n; k++)
      {
        p2=mulii(gcoeff(x,i,k),gcoeff(y,k,j));
        if (p2 != gen_0) p1=addii(p1,p2);
      }
      gcoeff(z,i,j) = gerepileupto(av,p1);
    }
  }
  return z;
}

static GEN
_mulmati(void *data /*ignored*/, GEN x, GEN y) {
  (void)data; return mulmati(x,y);
}
static GEN
_sqrmati(void *data /*ignored*/, GEN x) {
  (void)data; return mulmati(x,x);
}

static GEN
powmati(GEN x, GEN n)
{
  pari_sp av = avma;
  GEN y = leftright_pow(x, n, NULL, &_sqrmati, &_mulmati);
  return gerepileupto(av,y);
}

static GEN
rtran(GEN v, GEN w, GEN q)
{
  pari_sp av,tetpil;
  GEN p1;

  if (signe(q))
  {
    av=avma; p1=gneg(gmul(q,w)); tetpil=avma;
    return gerepile(av,tetpil,gadd(v,p1));
  }
  return v;
}

/* return (v - qw) mod m (only compute entries k0,..,n)
 * v and w are expected to have entries smaller than m */
static GEN
mtran(GEN v, GEN w, GEN q, GEN m, GEN mo2, long k0)
{
  long k;
  GEN p1;

  if (signe(q))
    for (k=lg(v)-1; k >= k0; k--)
    {
      pari_sp av = avma;
      p1 = subii(gel(v,k), mulii(q,gel(w,k)));
      p1 = centermodii(p1, m, mo2);
      gel(v,k) = gerepileuptoint(av, p1);
    }
  return v;
}

/* entries of v and w are C small integers */
static GEN
mtran_long(GEN v, GEN w, long q, long m, long k0)
{
  long k, p1;

  if (q)
  {
    for (k=lg(v)-1; k>= k0; k--)
    {
      p1 = v[k] - q * w[k];
      v[k] = p1 % m;
    }
  }
  return v;
}

/* coeffs of a are C-long integers */
static void
rowred_long(GEN a, long rmod)
{
  long q,j,k,pro, c = lg(a), r = lg(a[1]);

  for (j=1; j<r; j++)
  {
    for (k=j+1; k<c; k++)
      while (coeff(a,j,k))
      {
	q = coeff(a,j,j) / coeff(a,j,k);
	pro=(long)mtran_long(gel(a,j),gel(a,k),q,rmod, j);
	a[j]=a[k]; a[k]=pro;
      }
    if (coeff(a,j,j) < 0)
      for (k=j; k<r; k++) coeff(a,k,j)=-coeff(a,k,j);
    for (k=1; k<j; k++)
    {
      q = coeff(a,j,k) / coeff(a,j,j);
      gel(a,k) = mtran_long(gel(a,k),gel(a,j),q,rmod, k);
    }
  }
  /* don't update the 0s in the last columns */
  for (j=1; j<r; j++)
    for (k=1; k<r; k++) gcoeff(a,j,k) = stoi(coeff(a,j,k));
}

static void
rowred(GEN a, GEN rmod)
{
  long j,k,pro, c = lg(a), r = lg(a[1]);
  pari_sp av=avma, lim=stack_lim(av,1);
  GEN q, rmodo2 = shifti(rmod,-1);

  for (j=1; j<r; j++)
  {
    for (k=j+1; k<c; k++)
      while (signe(gcoeff(a,j,k)))
      {
	q=diviiround(gcoeff(a,j,j),gcoeff(a,j,k));
	pro=(long)mtran(gel(a,j),gel(a,k),q,rmod,rmodo2, j);
	a[j]=a[k]; a[k]=pro;
      }
    if (signe(gcoeff(a,j,j)) < 0)
      for (k=j; k<r; k++) gcoeff(a,k,j) = negi(gcoeff(a,k,j));
    for (k=1; k<j; k++)
    {
      q=diviiround(gcoeff(a,j,k),gcoeff(a,j,j));
      gel(a,k) = mtran(gel(a,k),gel(a,j),q,rmod,rmodo2, k);
    }
    if (low_stack(lim, stack_lim(av,1)))
    {
      long j1,k1;
      GEN p1 = a;
      if(DEBUGMEM>1) pari_warn(warnmem,"rowred j=%ld", j);
      p1 = gerepilecopy(av,a);
      for (j1=1; j1<r; j1++)
        for (k1=1; k1<c; k1++) gcoeff(a,j1,k1) = gcoeff(p1,j1,k1);
    }
  }
}

/* Compute d/x where d is t_INT, x lower triangular t_MAT with t_INT coeffs
 * whose diagonal coeffs divide d (lower triangular ZM result). */
static GEN
matinv(GEN x, GEN d)
{
  pari_sp av,av1;
  long i,j,k, n = lg(x[1]); /* Warning: lg(x) from ordmax is bogus */
  GEN y,h;

  y = matid(n-1);
  for (i=1; i<n; i++)
    gcoeff(y,i,i) = diviiexact(d,gcoeff(x,i,i));
  av=avma;
  for (i=2; i<n; i++)
    for (j=i-1; j; j--)
    {
      for (h=gen_0,k=j+1; k<=i; k++)
      {
        GEN p1 = mulii(gcoeff(y,i,k),gcoeff(x,k,j));
        if (p1 != gen_0) h=addii(h,p1);
      }
      setsigne(h,-signe(h)); av1=avma;
      gcoeff(y,i,j) = gerepile(av,av1,diviiexact(h,gcoeff(x,j,j)));
      av = avma;
    }
  return y;
}

static GEN
ordmax(GEN *cf, GEN p, long epsilon, GEN *ptdelta)
{
  long sp,i,n=lg(cf)-1;
  pari_sp av=avma, av2,limit;
  GEN T,T2,Tn,m,v,delta,hard_case_exponent, *w;
  const GEN pp = sqri(p);
  const GEN ppo2 = shifti(pp,-1);
  const long pps = (2*expi(pp)+2 < (long)BITS_IN_LONG)? pp[2]: 0;

  if (cmpiu(p,n) > 0)
  {
    hard_case_exponent = NULL;
    sp = 0; /* gcc -Wall */
  }
  else
  {
    long k;
    k = sp = itos(p);
    i=1; while (k < n) { k *= sp; i++; }
    hard_case_exponent = utoipos(i);
  }
  T=cgetg(n+1,t_MAT); for (i=1; i<=n; i++) gel(T,i) = cgetg(n+1,t_COL);
  T2=cgetg(2*n+1,t_MAT); for (i=1; i<=2*n; i++) gel(T2,i) = cgetg(n+1,t_COL);
  Tn=cgetg(n*n+1,t_MAT); for (i=1; i<=n*n; i++) gel(Tn,i) = cgetg(n+1,t_COL);
  v = new_chunk(n+1);
  w = (GEN*)new_chunk(n+1);

  av2 = avma; limit = stack_lim(av2,1);
  delta=gen_1; m=matid(n);

  for(;;)
  {
    long j, k, h;
    pari_sp av0 = avma;
    GEN t,b,jp,hh,index,p1, dd = sqri(delta), ppdd = mulii(dd,pp);
    GEN ppddo2 = shifti(ppdd,-1);

    if (DEBUGLEVEL > 3)
      fprintferr("ROUND2: epsilon = %ld\tavma = %ld\n",epsilon,avma);

    b=matinv(m,delta);
    for (i=1; i<=n; i++)
    {
      for (j=1; j<=n; j++)
        for (k=1; k<=n; k++)
        {
          p1 = j==k? gcoeff(m,i,1): gen_0;
          for (h=2; h<=n; h++)
          {
	    GEN p2 = mulii(gcoeff(m,i,h),gcoeff(cf[h],j,k));
            if (p2!=gen_0) p1 = addii(p1,p2);
          }
          gcoeff(T,j,k) = centermodii(p1, ppdd, ppddo2);
        }
      p1 = mulmati(m, mulmati(T,b));
      for (j=1; j<=n; j++)
	for (k=1; k<=n; k++)
	  gcoeff(p1,j,k) = centermodii(diviiexact(gcoeff(p1,j,k),dd),pp,ppo2);
      w[i] = p1;
    }

    if (hard_case_exponent)
    {
      for (j=1; j<=n; j++)
      {
	for (i=1; i<=n; i++) gcoeff(T,i,j) = gcoeff(w[j],1,i);
	/* ici la boucle en k calcule la puissance p mod p de w[j] */
	for (k=1; k<sp; k++)
	{
	  for (i=1; i<=n; i++)
	  {
	    p1 = gen_0;
	    for (h=1; h<=n; h++)
            {
              GEN p2=mulii(gcoeff(T,h,j),gcoeff(w[j],h,i));
	      if (p2!=gen_0) p1 = addii(p1,p2);
            }
            gel(v,i) = modii(p1, p);
	  }
	  for (i=1; i<=n; i++) gcoeff(T,i,j) = gel(v,i);
	}
      }
      t = powmati(T, hard_case_exponent);
    }
    else
    {
      for (i=1; i<=n; i++)
	for (j=1; j<=n; j++)
	{
          pari_sp av1 = avma;
          p1 = gen_0;
	  for (k=1; k<=n; k++)
	    for (h=1; h<=n; h++)
	    {
	      const GEN r=modii(gcoeff(w[i],k,h),p);
	      const GEN s=modii(gcoeff(w[j],h,k),p);
              const GEN p2 = mulii(r,s);
	      if (p2!=gen_0) p1 = addii(p1,p2);
	    }
	  gcoeff(T,i,j) = gerepileupto(av1,p1);
	}
      t = T;
    }

    if (pps)
    {
      long ps = p[2];
      for (i=1; i<=n; i++)
        for (j=1; j<=n; j++)
        {
          coeff(T2,j,i)=(i==j)? ps: 0;
          coeff(T2,j,n+i)=smodis(gcoeff(t,i,j),ps);
        }
      rowred_long(T2,pps);
    }
    else
    {
      for (i=1; i<=n; i++)
        for (j=1; j<=n; j++)
        {
          gcoeff(T2,j,i)=(i==j)? p: gen_0;
          gcoeff(T2,j,n+i) = modii(gcoeff(t,i,j),p);
        }
      rowred(T2,pp);
    }
    jp=matinv(T2,p);
    if (pps)
    {
      for (k=1; k<=n; k++)
      {
        pari_sp av1=avma;
        t = mulmati(mulmati(jp,w[k]), T2);
        for (h=i=1; i<=n; i++)
          for (j=1; j<=n; j++)
            { coeff(Tn,k,h) = itos(diviiexact(gcoeff(t,i,j), p)) % pps; h++; }
        avma=av1;
      }
      avma = av0;
      rowred_long(Tn,pps);
    }
    else
    {
      for (k=1; k<=n; k++)
      {
        t = mulmati(mulmati(jp,w[k]), T2);
        for (h=i=1; i<=n; i++)
          for (j=1; j<=n; j++)
            { gcoeff(Tn,k,h) = diviiexact(gcoeff(t,i,j), p); h++; }
      }
      rowred(Tn,pp);
    }
    for (index=gen_1,i=1; i<=n; i++)
      index = mulii(index,gcoeff(Tn,i,i));
    if (gcmp1(index)) break;

    m = mulmati(matinv(Tn,index), m);
    hh = delta = mulii(index,delta);
    for (i=1; i<=n; i++)
      for (j=1; j<=n; j++) hh = gcdii(gcoeff(m,i,j),hh);
    if (!is_pm1(hh))
    {
      m = gdiv(m,hh);
      delta = diviiexact(delta,hh);
    }
    epsilon -= 2 * Z_pval(index,p);
    if (epsilon < 2) break;
    if (low_stack(limit,stack_lim(av2,1)))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"ordmax");
      gerepileall(av2, 2, &m, &delta);
    }
  }
  gerepileall(av, 2, &m, &delta);
  *ptdelta = delta; return m;
}

/* Input:
 *  x normalized integral polynomial of degree n, defining K=Q(theta).
 *
 *  code 0, 1 or (long)p if we want base, smallbase ou factoredbase (resp.).
 *  y is GEN *, which will receive the discriminant of K.
 *
 * Output
 *  1) A t_COL whose n components are rationnal polynomials (with degree
 *     0,1...n-1) : integral basis for K (putting x=theta).
 *     Rem: common denominator is in da.
 *
 *  2) discriminant of K (in *y).
 */
static GEN
allbase2(GEN f, long flag, GEN *dx, GEN *dK, GEN *ptw)
{
  GEN w,w1,w2,a,pro,at,bt,b,da,db,q, *cf,*gptr[2];
  pari_sp av=avma,tetpil;
  long n,h,j,i,k,r,s,t,mf;

  w = ptw? *ptw: NULL;
  allbase_check_args(f,flag,dx, &w);
  w1 = gel(w,1);
  w2 = gel(w,2);
  n = degpol(f); h = lg(w1)-1;
  cf = (GEN*)cgetg(n+1,t_VEC);
  cf[2]=companion(f);
  for (i=3; i<=n; i++) cf[i]=mulmati(cf[2],cf[i-1]);

  a=matid(n); da=gen_1;
  for (i=1; i<=h; i++)
  {
    pari_sp av1 = avma;
    mf=itos(gel(w2,i)); if (mf==1) continue;
    if (DEBUGLEVEL) fprintferr("Treating p^k = %Z^%ld\n",w1[i],mf);

    b=ordmax(cf,gel(w1,i),mf,&db);
    a=gmul(db,a); b=gmul(da,b);
    da=mulii(db,da);
    at=gtrans(a); bt=gtrans(b);
    for (r=n; r; r--)
      for (s=r; s; s--)
        while (signe(gcoeff(bt,s,r)))
        {
          q=diviiround(gcoeff(at,s,s),gcoeff(bt,s,r));
          pro=rtran(gel(at,s),gel(bt,r),q);
          for (t=s-1; t; t--)
          {
            q=diviiround(gcoeff(at,t,s),gcoeff(at,t,t));
            pro=rtran(pro,gel(at,t),q);
          }
          at[s]=bt[r]; gel(bt,r) = pro;
        }
    for (j=n; j; j--)
    {
      for (k=1; k<j; k++)
      {
        while (signe(gcoeff(at,j,k)))
        {
          q=diviiround(gcoeff(at,j,j),gcoeff(at,j,k));
          pro=rtran(gel(at,j),gel(at,k),q);
          at[j]=at[k]; gel(at,k) = pro;
        }
      }
      if (signe(gcoeff(at,j,j))<0)
        for (k=1; k<=j; k++) gcoeff(at,k,j) = negi(gcoeff(at,k,j));
      for (k=j+1; k<=n; k++)
      {
        q=diviiround(gcoeff(at,j,k),gcoeff(at,j,j));
        gel(at,k) = rtran(gel(at,k),gel(at,j),q);
      }
    }
    for (j=2; j<=n; j++)
      if (equalii(gcoeff(at,j,j), gcoeff(at,j-1,j-1)))
      {
        gcoeff(at,1,j)= gen_0;
        for (k=2; k<=j; k++) gcoeff(at,k,j) = gcoeff(at,k-1,j-1);
      }
    tetpil=avma; a=gtrans(at);
    {
      GEN *gptr[2];
      da = icopy(da); gptr[0]=&a; gptr[1]=&da;
      gerepilemanysp(av1,tetpil,gptr,2);
    }
  }
  *dK = *dx;
  for (j=1; j<=n; j++)
    *dK = diviiexact(mulii(*dK,sqri(gcoeff(a,j,j))), sqri(da));
  tetpil=avma; *dK = icopy(*dK);
  at=cgetg(n+1,t_VEC);
  for (k=1; k<=n; k++)
  {
    q=cgetg(k+2,t_POL); q[1] = f[1]; gel(at,k) = q;
    for (j=1; j<=k; j++) gel(q,j+1) = gdiv(gcoeff(a,k,j),da);
  }
  gptr[0] = &at; gptr[1] = dK;
  gerepilemanysp(av,tetpil,gptr,2);
  return at;
}

GEN
base2(GEN x, GEN *pdK) { return nfbasis(x, pdK, compat_ROUND2, NULL); }

GEN
discf2(GEN x) { return nfdiscf0(x, compat_ROUND2, NULL); }

/*******************************************************************/
/*                                                                 */
/*                            ROUND 4                              */
/*                                                                 */
/*******************************************************************/
GEN maxord_i(GEN p, GEN f, long mf, GEN w, long flag);
static GEN dbasis(GEN p, GEN f, long mf, GEN alpha, GEN U);
static GEN maxord(GEN p,GEN f,long mf);

static int
fnz(GEN x,long j)
{
  long i;
  for (i=1; i<j; i++)
    if (signe(x[i])) return 0;
  return 1;
}

/* return list u[i], 2 by 2 coprime with the same prime divisors as ab */
static GEN
get_coprimes(GEN a, GEN b)
{
  long i, k = 1;
  GEN u = cgetg(3, t_COL);
  gel(u,1) = a;
  gel(u,2) = b;
  /* u1,..., uk 2 by 2 coprime */
  while (k+1 < lg(u))
  {
    GEN d, c = gel(u,k+1);
    if (is_pm1(c)) { k++; continue; }
    for (i=1; i<=k; i++)
    {
      if (is_pm1(u[i])) continue;
      d = gcdii(c, gel(u,i));
      if (d == gen_1) continue;
      c = diviiexact(c, d);
      gel(u,i) = diviiexact(gel(u,i), d);
      u = shallowconcat(u, d);
    }
    gel(u,++k) = c;
  }
  for (i = k = 1; i < lg(u); i++)
    if (!is_pm1(u[i])) u[k++] = u[i];
  setlg(u, k); return u;
}

/* allow p = -1 from factorizations */
static long
safe_Z_pvalrem(GEN x, GEN p, GEN *z)
{
  if (signe(p) < 0) { *z = absi(x); return 1; }
  return Z_pvalrem(x, p, z);
}

/* return integer basis. Set dK = disc(K), dx = disc(f), w (possibly partial)
 * factorization of dK. *ptw can be set by the caller, in which case it is
 * taken to be the factorization of disc(f), then overwritten
 * [no consistency check] */
GEN
allbase(GEN f, long flag, GEN *dx, GEN *dK, GEN *index, GEN *ptw)
{
  VOLATILE GEN w1, w2, a, da, ordmax;
  VOLATILE long n, lw, i, j, k, l;
  GEN w;

  if (flag & nf_ROUND2) return allbase2(f,flag,dx,dK,ptw);
  w = ptw? *ptw: NULL;
  allbase_check_args(f, flag, dx, &w);
  w1 = gel(w,1);
  w2 = vec_to_vecsmall(gel(w,2));
  n = degpol(f); lw = lg(w1);
  ordmax = cgetg(1, t_VEC);
  /* "complete" factorization first */
  for (i=1; i<lw; i++)
  {
    if (w2[i] == 1) { ordmax = shallowconcat(ordmax, gen_1); continue; }

    CATCH(invmoder) { /* caught false prime, update factorization */
      GEN x = (GEN)global_err_data;
      GEN p = gcdii(gel(x,1), gel(x,2));
      GEN N, u;
      if (DEBUGLEVEL) pari_warn(warner,"impossible inverse: %Z", x);

      u = get_coprimes(p, diviiexact(gel(x,1),p));
      l = lg(u);
      /* no small factors, but often a prime power */
      for (k = 1; k < l; k++) gel(u,k) = gcoeff(auxdecomp(gel(u,k), 2),1,1);

      w1[i] = u[1];
      w1 = shallowconcat(w1, vecslice(u, 2, l-1));
      N = *dx;
      w2[i] = Z_pvalrem(N, gel(w1,i), &N);
      k  = lw;
      lw = lg(w1);
      for ( ; k < lw; k++) w2[k] = Z_pvalrem(N, gel(w1,k), &N);
    } RETRY {
      if (DEBUGLEVEL) fprintferr("Treating p^k = %Z^%ld\n",w1[i],w2[i]);
      ordmax = shallowconcat(ordmax, mkvec( maxord(gel(w1,i),f,w2[i]) ));
    } ENDCATCH;
  }

  a = NULL; /* gcc -Wall */
  da = NULL;
  for (i=1; i<lw; i++)
  {
    GEN M, db, b = gel(ordmax,i);
    if (b == gen_1) continue;
    db = gen_1;
    for (j=1; j<=n; j++)
    {
      GEN t = denom(gcoeff(b,j,j));
      if (absi_cmp(t,db) > 0) db = t;
    }
    if (db == gen_1) continue;

    /* db = denom(b), (da,db) = 1. Compute da Im(b) + db Im(a) */
    b = Q_muli_to_int(b,db);
    if (!da) { da = db; a = b; }
    else
    { /* optimization: easy as long as both matrix are diagonal */
      j=2; while (j<=n && fnz(gel(a,j),j) && fnz(gel(b,j),j)) j++;
      k = j-1; M = cgetg(2*n-k+1,t_MAT);
      for (j=1; j<=k; j++)
      {
        gel(M,j) = gel(a,j);
        gcoeff(M,j,j) = mulii(gcoeff(a,j,j),gcoeff(b,j,j));
      }
      /* could reduce mod M(j,j) but not worth it: usually close to da*db */
      for (  ; j<=n;     j++) gel(M,j) = gmul(db, gel(a,j));
      for (  ; j<=2*n-k; j++) gel(M,j) = gmul(da, gel(b,j+k-n));
      da = mulii(da,db);
      a = hnfmodid(M, da);
    }
    if (DEBUGLEVEL>5) fprintferr("Result for prime %Z is:\n%Z\n",w1[i],b);
  }
  if (da)
  {
    *index = diviiexact(da, gcoeff(a,1,1));
    for (j=2; j<=n; j++) *index = mulii(*index, diviiexact(da, gcoeff(a,j,j)));
    a = gdiv(hnfcenter_ip(a), da);
  }
  else
  {
    *index = gen_1;
    a = matid(n);
  }
  *dK = diviiexact(*dx, sqri(*index));

  if (ptw)
  {
    long lfa = 1;
    GEN W1, W2, D = *dK;
    W1 = cgetg(lw, t_COL);
    W2 = cgetg(lw, t_COL);
    for (j=1; j<lw; j++)
    {
      k = safe_Z_pvalrem(D, gel(w1,j), &D);
      if (k) { gel(W1,lfa) = gel(w1,j); gel(W2,lfa) = utoipos(k); lfa++; }
    }
    setlg(W1, lfa);
    setlg(W2, lfa); *ptw = mkmat2(W1,W2);
  }
  return RgM_to_RgXV(a, varn(f));
}

static GEN
update_fact(GEN x, GEN f)
{
  GEN e, q, d = ZX_disc(x), p = gel(f,1);
  long iq, i, k, l;
  if (typ(f)!=t_MAT || lg(f)!=3) pari_err(talker,"not a factorisation in nfbasis");
  l = lg(p);
  q = cgetg(l,t_COL); 
  e = cgetg(l,t_COL); iq = 1;
  for (i=1; i<l; i++)
  {
    k = safe_Z_pvalrem(d, gel(p,i), &d);
    if (k) { q[iq] = p[i]; gel(e,iq) = utoipos(k); iq++; }
  }
  setlg(q,iq); setlg(e,iq);
  return merge_factor_i(Z_factor(d), mkmat2(q,e));
}

/* FIXME: have to deal with compatibility flags */
static void
_nfbasis(GEN x0, long flag, GEN fa, GEN *pbas, GEN *pdK)
{
  GEN x, dx, dK, basis, lead, index;
  long fl = 0;

  if (typ(x0)!=t_POL) pari_err(typeer,"nfbasis");
  if (!degpol(x0)) pari_err(zeropoler,"nfbasis");
  check_ZX(x0, "nfbasis");

  x = pol_to_monic(x0, &lead);
  if (fa && gcmp0(fa)) fa = NULL; /* compatibility. NULL is the proper arg */
  if (fa && lead) fa = update_fact(x, fa);
  if (flag & compat_PARTIAL) fl |= nf_PARTIALFACT;
  if (flag & compat_ROUND2)  fl |= nf_ROUND2;
  basis = allbase(x, fl, &dx, &dK, &index, &fa);
  if (pbas) *pbas = RgXV_unscale(basis, lead);
  if (pdK)  *pdK = dK;
}

GEN
nfbasis(GEN x, GEN *pdK, long flag, GEN fa)
{
  pari_sp av = avma;
  GEN bas; _nfbasis(x, flag, fa, &bas, pdK);
  gerepileall(av, pdK? 2: 1, &bas, pdK); return bas;
}

GEN
nfbasis0(GEN x, long flag, GEN fa)
{
  pari_sp av = avma;
  GEN bas; _nfbasis(x, flag, fa, &bas, NULL);
  return gerepilecopy(av, bas);
}

GEN
nfdiscf0(GEN x, long flag, GEN fa)
{
  pari_sp av = avma;
  GEN dK; _nfbasis(x, flag, fa, NULL, &dK);
  return gerepilecopy(av, dK);
}

GEN
base(GEN x, GEN *pdK) { return nfbasis(x, pdK, 0, NULL); }

GEN
smallbase(GEN x, GEN *pdK) { return nfbasis(x, pdK, compat_PARTIAL, NULL); }

GEN
factoredbase(GEN x, GEN fa, GEN *pdK) { return nfbasis(x, pdK, 0, fa); }

GEN
discf(GEN x) { return nfdiscf0(x, 0, NULL); }

GEN
smalldiscf(GEN x) { return nfdiscf0(x, nf_PARTIALFACT, NULL); }

GEN
factoreddiscf(GEN x, GEN fa) { return nfdiscf0(x, 0, fa); }


/* return U if Z[alpha] is not maximal or 2*dU < m-1; else return NULL */
static GEN
dedek(GEN f, long mf, GEN p,GEN g)
{
  GEN k,h;
  long dk;

  h = FpX_div(f,g,p);
  k = gdivexact(gadd(f, gneg_i(gmul(g,h))), p);
  k = FpX_gcd(k, FpX_gcd(g,h, p), p);

  dk = degpol(k);
  if (DEBUGLEVEL>2)
  {
    fprintferr("  dedek: gcd has degree %ld\n", dk);
    if (DEBUGLEVEL>5) fprintferr("initial parameters p=%Z,\n  f=%Z\n",p,f);
  }
  if (2*dk >= mf-1) return FpX_div(f,k,p);
  return dk? (GEN)NULL: f;
}

/* p-maximal order of Af; mf = v_p(Disc(f)) */
static GEN
maxord(GEN p,GEN f,long mf)
{
  const pari_sp av = avma;
  GEN w = NULL, g, res, fp = FpX_red(f, p);

  if (cmpui(degpol(f),p) < 0)
    g = FpX_div(fp, FpX_gcd(fp,derivpol(fp), p), p);
  else
  {
    w = (GEN)FpX_factor(fp,p)[1];
    g = FpXV_prod(w, p);
  }
  res = dedek(f, mf, p, g);
  if (res)
    res = dbasis(p, f, mf, pol_x[varn(f)], res);
  else
  {
    if (!w) w = (GEN)FpX_factor(fp,p)[1];
    res = maxord_i(p, f, mf, w, 0);
  }
  return gerepileupto(av,res);
}

/* Sylvester's matrix, mod p^m (assumes f1 monic) */
static GEN
sylpm(GEN f1, GEN f2, GEN pm)
{
  long j, n = degpol(f1);
  GEN h, a = cgetg(n+1,t_MAT);
  h = FpX_rem(f2,f1,pm);
  for (j=1;; j++)
  {
    gel(a,j) = RgX_to_RgV(h, n);
    if (j == n) break;
    h = FpX_rem(RgX_shift_shallow(h, 1), f1, pm);
  }
  return hnfmodidpart(a, pm);
}

/* polynomial gcd mod p^m (assumes f1 monic) */
GEN
gcdpm(GEN f1, GEN f2, GEN pm)
{
  pari_sp av = avma;
  long c, n = degpol(f1), v = varn(f1);
  GEN col, a = sylpm(f1,f2,pm);
  for (c = 1; c <= n; c++)
    if (!equalii(gcoeff(a,c,c), pm))
    {
      col = gdiv(gel(a,c), gcoeff(a,c,c));
      return gerepilecopy(av, RgV_to_RgX(col,v));
    }
  avma = av; return zeropol(v);
}

/* reduced resultant mod p^m (assumes x monic) */
GEN
respm(GEN x, GEN y, GEN pm)
{
  pari_sp av = avma;
  GEN z = sylpm(x,y,pm);
  z = gcoeff(z,1,1);
  if (equalii(z,pm)) { avma = av; return gen_0; }
  return gerepileuptoint(av, icopy(z));
}

static void
update_den(GEN *e, GEN *de, GEN *pp)
{
  GEN ce = Q_content(*e);
  if (ce != gen_1) {
    ce = gcdii(*de, ce);
    *de = diviiexact(*de, ce);
    *e  = gdivexact(*e, ce);
    if (pp) *pp =diviiexact(*pp, ce);
  }
}

/* f o g mod (T,p) */
static GEN
compmod(GEN f, GEN g, GEN T, GEN p)
{
  GEN D = NULL, z, df, dg, q;
  f = Q_remove_denom(f, &df);
  g = Q_remove_denom(g, &dg);
  if (df) D = df;
  if (dg) D = mul_content(D, powiu(dg, degpol(f)));
  q = D ? mulii(p, D): p;
  if (dg) f = FpX_rescale(f, dg, q);
  z = FpX_FpXQ_compo(f, g, T, q);
  if (!D) return z;
  update_den(&z, &D, NULL);
  return gdiv( FpX_center(z, mulii(D,p)), D );
}

static GEN
dbasis(GEN p, GEN f, long mf, GEN a, GEN U)
{
  long n = degpol(f), dU, i;
  GEN b, ha, pd, pdp;

  if (n == 1) return gscalmat(gen_1, 1);
  if (DEBUGLEVEL>5)
  {
    fprintferr("  entering Dedekind Basis with parameters p=%Z\n",p);
    fprintferr("  f = %Z,\n  a = %Z\n",f,a);
  }
  ha = pd = powiu(p,mf/2); pdp = mulii(pd,p);
  dU = U ? degpol(U): 0;
  b = cgetg(n, t_MAT); /* Z[a] + U/p Z[a] is maximal */
  /* skip first column = [pd, 0,...,0] */
  for (i=1; i<n; i++)
  {
    if (i == dU)
      ha = gmul(diviiexact(pd, p), compmod(U, a, f, pdp));
    else
    {
      GEN d, mod;
      ha = Q_remove_denom(gmul(ha,a), &d);
      mod = d? mulii(pdp, d): pdp;
      ha = FpX_rem(ha, f, mod);
      if (d) ha = gdivexact(ha,d);
    }
    gel(b,i) = RgX_to_RgV(ha,n);
  }
  b = hnfmodid(b,pd);
  if (DEBUGLEVEL>5) fprintferr("  new order: %Z\n",b);
  return gdiv(b, pd);
}

static GEN
get_partial_order_as_pols(GEN p, GEN f, GEN *d)
{
  GEN b = maxord(p,f, Z_pval(ZX_disc(f),p));
  GEN z = Q_remove_denom( RgM_to_RgXV(b, varn(f)), d );
  if (!*d) *d = gen_1;
  return z;
}

long
FpX_val(GEN x0, GEN t, GEN p, GEN *py)
{
  long k;
  GEN r, y, x = x0;

  for (k=0; ; k++)
  {
    y = FpX_divrem(x, t, p, &r);
    if (signe(r)) break;
    x = y;
  }
  *py = x; return k;
}

#if 0
/* e in Qp, f i Zp. Return r = e mod (f, pk) */
static GEN
QpX_mod(GEN e, GEN f, GEN pk)
{
  GEN mod, d;
  e = Q_remove_denom(e, &d);
  mod = d? mulii(pk,d): pk;
  e = FpX_rem(e, centermod(f, mod), mod);
  e = FpX_center(e, mod);
  if (d) e = gdiv(e, d);
  return e;
}
#endif

typedef struct __decomp {
/* constants */
  GEN p, f; /* goal: factor f p-adically */
  long df; /* p^df = reduced discriminant of f */
/* these are updated along the way */
  GEN phi; /* a p-integer, in Q[X] */
  GEN phi0; /* a p-integer, in Q[X] from testb2 / testc2, to be composed with
             * phi when correct precision is known */
  GEN chi; /* characteristic polynomial of phi (mod p^*), in Z[X] */
  GEN nu; /* irreducible divisor of chi mod p, in Z[X] */
  GEN invnu; /* numerator ( 1/ Mod(nu, chi) mod pmr ) */
  GEN Dinvnu;/* denominator ( ... ) */
  GEN pdr, pmr, pmf;
} decomp_t;

/* if flag = 0, maximal order, else factorization to precision r = flag */
static GEN
Decomp(decomp_t *S, long flag)
{
  GEN fred, res, pr, pk, ph, b1, b2, a, e, de, f1, f2, dt, th;
  GEN p = S->p;
  long k, r = flag? flag: 2*S->df + 1;

  if (DEBUGLEVEL>2)
  {
    fprintferr("  entering Decomp");
    if (DEBUGLEVEL>5) fprintferr(", parameters: %Z^%ld\n  f = %Z",p, r, S->f);
    fprintferr("\n");
  }
  if (!FpX_val(S->chi, S->nu, p, &b1))
    pari_err(talker, "bug in Decomp (not a factor), is p = %Z a prime?", p);
  b2 = FpX_div(S->chi, b1, p);
  a = FpX_mul(FpXQ_inv(b2, b1, p), b2, p);
  /* E = e / de, e in Z[X], de in Z,  E = a(phi) mod (f, p) */
  th = Q_remove_denom(S->phi, &dt);
  if (!dt) dt = gen_1;
  de = powiu(dt, degpol(a));
  pr = mulii(p, de);
  e = FpX_FpXQ_compo(FpX_rescale(a, dt, pr), th, S->f, pr);
  update_den(&e, &de, NULL);

  pk = p; k = 1;
  /* E, (1 - E) tend to orthogonal idempotents in Zp[X]/(f) */
  while (k < r + Z_pval(de, p))
  { /* E <-- E^2(3-2E) mod p^2k, with E = e/de */
    GEN D;
    pk = sqri(pk); k <<= 1;
    e = gmul(gsqr(e), gsub(mulsi(3,de), gmul2n(e,1)));
    de= mulii(de, sqri(de));
    D = mulii(pk, de);
    e = FpX_rem(e, centermod(S->f, D), D); /* e/de defined mod pk */
    update_den(&e, &de, NULL);
  }
  pr = powiu(p, r); /* required precision of the factors */
  ph = mulii(de, pr);
  fred = centermod(S->f, ph);
  e    = centermod(e, ph);

  f1 = gcdpm(fred, gsub(de, e), ph); /* p-adic gcd(f, 1-e) */
  fred = centermod(fred, pr);
  f1   = centermod(f1,   pr);
  f2 = FpX_div(fred,f1, pr);
  f2 = FpX_center(f2, pr);

  if (DEBUGLEVEL>5)
    fprintferr("  leaving Decomp: f1 = %Z\nf2 = %Z\ne = %Z\nde= %Z\n", f1,f2,e,de);

  if (flag)
    return concat_factor(ZX_monic_factorpadic(f1, p, flag),
                         ZX_monic_factorpadic(f2, p, flag));
  else
  {
    GEN D = de, d1, d2, ib1, ib2;
    long n, n1, n2, i;
    ib1 = get_partial_order_as_pols(p,f1, &d1); n1 = lg(ib1)-1;
    ib2 = get_partial_order_as_pols(p,f2, &d2); n2 = lg(ib2)-1; n = n1+n2;
    i = cmpii(d1, d2);
    if (i < 0) {
      ib1 = gmul(ib1, diviiexact(d2, d1));
      d1 = d2;
    }
    else if (i > 0) {
      ib2 = gmul(ib2, diviiexact(d1, d2));
    }
    D = mulii(d1, D);
    fred = centermod(S->f, D);
    res = cgetg(n+1, t_VEC);
    for (i=1; i<=n1; i++)
      gel(res,i) = FpX_center(FpX_rem(gmul(gel(ib1,i),e), fred, D), D);
    e = gsub(de, e); ib2 -= n1;
    for (   ; i<=n; i++)
      gel(res,i) = FpX_center(FpX_rem(gmul(gel(ib2,i),e), fred, D), D);
    res = RgXV_to_RgM(res, n);
    return gdiv(hnfmodid(res,D), D); /* normalized integral basis */
  }
}

/* minimum extension valuation: L/E */
static void
vstar(GEN p,GEN h, long *L, long *E)
{
  long first, j, k, v, w, m = degpol(h);

  first = 1; k = 1; v = 0;
  for (j=1; j<=m; j++)
    if (! gcmp0(gel(h,m-j+2)))
    {
      w = Z_pval(gel(h,m-j+2),p);
      if (first || w*k < v*j) { v = w; k = j; }
      first = 0;
    }
  w = cgcd(v,k);
  *L = v/w;
  *E = k/w;
}

static GEN
redelt_i(GEN a, GEN N, GEN p, GEN *pda)
{
  GEN z;
  a = Q_remove_denom(a, pda);
  if (*pda)
  {
    long v = Z_pvalrem(*pda, p, &z);
    if (v) {
      *pda = powiu(p, v);
      N  = mulii(*pda, N);
    }
    else
      *pda = NULL;
    if (!is_pm1(z)) a = gmul(a, Fp_inv(z, N));
  }
  return centermod(a, N);
}
/* reduce the element a modulo N [ a power of p ], taking first care of the
 * denominators */
static GEN
redelt(GEN a, GEN N, GEN p)
{
  GEN da; a = redelt_i(a, N, p, &da);
  if (da) a = gdiv(a, da);
  return a;
}

/* compute the Newton sums of g(x) mod p, assume deg g > 0 */
GEN
polsymmodp(GEN g, GEN p)
{
  pari_sp av1, av2;
  long d = degpol(g), i, k;
  GEN s , y;

  y = cgetg(d + 1, t_COL);
  gel(y,1) = utoipos(d);
  for (k = 1; k < d; k++)
  {
    av1 = avma;
    s = centermod(mulsi(k, polcoeff0(g,d-k,-1)), p);
    for (i = 1; i < k; i++)
      s = addii(s, mulii(gel(y,k-i+1), polcoeff0(g,d-i,-1)));
    av2 = avma;
    gel(y,k+1) = gerepile(av1, av2, centermod(negi(s), p));
  }

  return y;
}

/* no GC */
static GEN
manage_cache(GEN chi, GEN pp, GEN ns)
{
  if (lgefint(pp) > lg(ns[1]))
  {
    if (DEBUGLEVEL > 4) fprintferr("newtonsums: result doesn't fit in cache\n");
    return polsymmodp(chi, pp);
  }

  if (!signe(ns[1]))
  {
    GEN ns2 = polsymmodp(chi, pp);
    long j, l = lg(ns);
    for (j = 1; j < l; j++) affii(gel(ns2,j), gel(ns,j));
  }
  return ns;
}

/* compute the c first Newton sums modulo pp of the
   characteristic polynomial of a/d mod chi, d > 0 power of p (NULL = gen_1),
   a, chi in Zp[X]
   ns = Newton sums of chi */
static GEN
newtonsums(GEN a, GEN da, GEN chi, long c, GEN pp, GEN ns)
{
  GEN va, pa, dpa, s;
  long j, k, n = degpol(chi);
  pari_sp av, lim;

  a = centermod(a, pp); av = avma; lim = stack_lim(av, 1);
  pa = pol_1[varn(a)]; dpa = gen_1;
  va = zerovec(c);
  for (j = 1; j <= c; j++)
  { /* pa/dpa = (a/d)^(j-1) mod (chi, pp) */
    pa = FpX_rem(FpX_mul(pa, a, pp), chi, pp);
    s  = gen_0;
    for (k = 0; k < n; k++)
      s = addii(s, mulii(polcoeff0(pa, k, -1), gel(ns,k+1)));
    if (da) {
      dpa = mulii(dpa, da);
      s = gdiv(s, dpa);
      if (typ(s) != t_INT) return NULL;
      update_den(&pa, &dpa, &pp);
    }

    gel(va,j) = centermod(s, pp);

    if (low_stack(lim, stack_lim(av, 1)))
    {
      if(DEBUGMEM>1) pari_warn(warnmem, "newtonsums");
      gerepileall(av, 4, &pa, &va, &pp, &dpa);
    }
  }
  return va;
}

/* compute the characteristic polynomial of a/da mod chi (a in Z[X]), given
 * by its Newton sums to a precision of pp using Newton sums */
static GEN
newtoncharpoly(GEN pp, GEN p, GEN NS)
{
  long n = lg(NS)-1, j, k;
  GEN c = cgetg(n + 2, t_VEC);

  if (!NS) return NULL;
  gel(c,1) = (n & 1 ? gen_m1: gen_1);
  for (k = 2; k <= n+1; k++) gel(c,k) = gen_0;
  for (k = 2; k <= n+1; k++)
  {
    pari_sp av2 = avma;
    GEN s = gen_0;
    ulong z;
    long v = u_pvalrem(k - 1, p, &z);
    for (j = 1; j < k; j++)
    {
      GEN t = mulii(gel(NS,j), gel(c,k-j));
      if (!odd(j)) t = negi(t);
      s = addii(s, t);
    }
    if (v) {
      s = gdiv(s, powiu(p, v));
      if (typ(s) != t_INT) return NULL;
    }
    s = mulii(s, Fp_inv(utoipos(z), pp));
    gel(c,k) = gerepileuptoint(av2, centermod(s, pp));
  }
  for (k = odd(n)? 1: 2; k <= n+1; k += 2) gel(c,k) = negi(gel(c,k));
  return gtopoly(c, 0);
}

/* guess if a mod chi has positive valuation
   by looking at the newton sums */
static long
fastvalpos(GEN a, GEN chi, GEN p, GEN ns, long E)
{
  GEN v, d, pp;
  long m, n = degpol(chi), j, c;

  c = equaliu(p, 2)? 2*n/3 : min(2*E, n);
  if (c < 2) c = 2;
  a = Q_remove_denom(a, &d);
  m = d? Z_pval(d, p): 0; /* >= 0 */
  pp = powiu(p, (m+1)*c+1);
  ns = manage_cache(chi, pp, ns);
  v = newtonsums(a, d, chi, c, pp, ns);
  if (!v) return 0;
  for (j = 1; j <= c; j++)
    if (signe(gel(v,j)) && E*Z_pval(gel(v,j), p) - j*(E*m+1) < 0) return 0;
  return 1;
}

/* return v_p(n!) */
long
val_fact(ulong n, ulong p)
{
  ulong q = p, v = 0;
  do { v += n/q; q *= p; } while (n >= q);
  return (long)v;
}

/* return NULL if a mod f is not an integer
 * if dr >= 0, a mod f is an integer and the denominator of any integer
 * in Zp[X]/(f) divides p^dr */
static GEN
mycaract(GEN f, GEN a, GEN p, GEN pp, long dr, GEN ns)
{
  pari_sp av = avma;
  GEN d, chi, npp, NPP, nspp;
  long n = degpol(f);

  if (gcmp0(a)) return zeropol(varn(f));

  a = Q_remove_denom(a, &d);
  npp = pp;
  if (lgefint(p) == 3) npp = mulii(npp, powiu(p, val_fact(n, itou(p))));
  nspp = NPP = npp;
  if (d) {
    NPP = mulii(NPP, powiu(d, n));
    nspp = (dr < 0)? NPP: mulii(nspp, powiu(p, dr));
  }
  ns = newtonsums(a, d, f, n, NPP, manage_cache(f, nspp, ns));
  if (!ns) return NULL;
  chi = newtoncharpoly(npp, p, ns);
  if (!chi) return NULL;
  setvarn(chi, varn(f));
  return gerepileupto(av, centermod(chi, pp));
}

static GEN
get_nu(GEN chi, GEN p, long *ptl)
{
  GEN P = (GEN)FpX_factor(chi, p)[1];
  *ptl = lg(P) - 1;
  return gel(P,*ptl);
}

/* Factor characteristic polynomial of S->phi mod (p, S->chi) */
static long
factcp(decomp_t *S, GEN ns)
{
  GEN chi = mycaract(S->chi, S->phi, S->p, S->pmf, -1, ns);
  long l;
  S->chi= chi;
  S->nu = get_nu(chi, S->p, &l); return l;
}

/* Compute nu_beta in Fp[X]. If something unexpected happens, return NULL */
static GEN
fastnu(GEN p, GEN f, GEN beta, GEN pdr)
{
  long j, k, l, n = degpol(f), v = varn(f), N = 2*n+1, av = avma;
  GEN p1, p2, c, d, B, G, V, nu, h;

  G   = cgetg(N+1, t_MAT);
  c  = gen_0;
  d  = mulii(pdr, sqri(p));

  beta = gmul(pdr, beta);
  B    = beta;
  for (k = 1; k <= n; k++)
  {
    V = zerocol(N); gel(G,N-k) = V;
    gel(V,n+1-k) = gen_1;
    for (j = n+1; j <= N; j++)
    {
      p2 = polcoeff0(B, N-j, -1);
      if (signe(p2)) c = gcdii(c, p2);
      gel(V,j) = p2;
    }
    if (k < n)
    {
      B = gdiv(grem(gmul(B, beta), f), pdr);
      if (!gcmp1(Q_denom(B))) { avma = av; return NULL; }
      B = centermod(B, d);
    }
  }

  if (DEBUGLEVEL >= 6)
    fprintferr(" content in fastnu is %Z\n", c);

  for (k = 1; k <= n; k++)
  {
    p1 = gel(G,N-k);
    for (j = n+1; j <= N; j++)
    {
      p2 = gel(p1,j);
      if (signe(p2)) { p2 = diviiexact(p2, c); gel(p1,j) = p2; }
    }
  }
  pdr = diviiexact(pdr, c);
  d   = diviiexact(d, c);

  V = zerocol(N); gel(G,N) = V;
  gel(V,N) = pdr; gel(V,n+1) = gen_1;

  p1 = mulii(pdr, p);
  for (k = 1; k <= n; k++)
  {
    V = zerocol(N); gel(G,k) = V;
    gel(V,n+k+1) = p1;
  }
  if (DEBUGLEVEL >= 6) fprintferr("  fastnu: G is computed\n");

  G = hnfmodidpart(G, d);
  if (DEBUGLEVEL >= 6) fprintferr("  fastnu: HNF(G) is computed\n");

  setlg(G, n+2);
  G = rowslice(G, 1, n+1);
  h = gtopoly(gel(G,n+1), v);
  for (j = 1; j <= n; j++)
    h = FpX_gcd(h, gtopoly(gel(G,j), v), p);

  if (!degpol(h)) { avma = av; return NULL; }
  nu = get_nu(h, p, &l);
  if (l > 1) { avma = av; return NULL; }
  return gerepilecopy(av, nu);
}

/* return the prime element in Zp[phi], nup, chip in Z[X]
 * if *Ep < oE or Ep divides Ediv (!=0) return NULL (not interesting)
 * */
static GEN
getprime(decomp_t *S, GEN phi, GEN chip, GEN nup, long *Lp, long *Ep, 
         long oE, long Ediv)
{
  GEN chin, q;
  long r, s;

  if (degpol(nup) == 1)
  {
    GEN c = gel(nup,2); /* nup = X + c */
    chin = signe(c)? translate_pol(chip, negi(c)): chip;
  }
  else
    chin = ZX_caract(chip, nup, varn(chip));

  vstar(S->p, chin, Lp, Ep);
  if (*Ep < oE || (Ediv && Ediv % *Ep == 0)) return NULL;

  if (*Ep == 1) return S->p;
  (void)cbezout(*Lp, -*Ep, &r, &s); /* = 1 */
  if (r <= 0)
  {
    long t = 1 + ((-r) / *Ep);
    r += t * *Ep;
    s += t * *Lp;
  }
  /* r > 0 minimal such that r L/E - s = 1/E
   * pi = nu^r / p^s is an element of valuation 1/E,
   * so is pi + O(p) since 1/E < 1. May compute nu^r mod p^(s+1) */

  q = powiu(S->p, s+1);
  nup = FpXQ_pow(nup, utoipos(r), S->chi, q);
  return gdiv(compmod(nup, phi, S->chi, q), powiu(S->p, s));
}

static void
kill_cache(GEN ns) { setsigne(ns[1], 0); }

/* S->phi := T o T0 mod (p, f) */
static void
composemod(decomp_t *S, GEN T, GEN T0) { S->phi = compmod(T, T0, S->f, S->p); }

static int 
update_phi(decomp_t *S, GEN ns, long *ptl, long flag)
{
  GEN PHI = NULL, pdr, pmr = S->pmr, X = pol_x[ varn(S->f) ];
  long k;

  if (!S->chi) 
  {
    kill_cache(ns);
    S->chi = mycaract(S->f, S->phi, S->p, pmr, S->df, ns);
    S->nu = get_nu(S->chi, S->p, ptl);
    if (*ptl > 1) return 0; /* we can get a decomposition */
  }
    
  for (k = 1;; k++)
  {
    kill_cache(ns);
    pdr = respm(S->chi, derivpol(S->chi), pmr);
    if (signe(pdr)) break;
    
    pmr = sqri(pmr); /* try a larger precision */

    PHI = S->phi0? compmod(S->phi, S->phi0, S->f, pmr): S->phi;
    PHI = gadd(PHI, gmul(mulsi(k, S->p), X));
    S->chi = mycaract(S->f, PHI, S->p, pmr, S->df, ns);
  }
  pmr = mulii(sqri(pdr), S->p);
  S->chi = FpX_red(S->chi, pmr);
  if (!PHI) /* ok above for k = 0 */
    PHI = S->phi0? compmod(S->phi, S->phi0, S->f, pmr): S->phi;
  S->phi = PHI;

  if (is_pm1(pdr))
  { /* may happen if p is unramified */
    if (!flag) { *ptl = 1; return 0; }
    S->nu = get_nu(S->chi, S->p, ptl);
    return 0;
  }
  S->pmr = pmr;
  S->pdr = mulii(pdr, S->p); return 1;
}

/* return 1 if at least 2 factors mod p ==> chi can be split
 * Replace S->phi such that F increases (to D) */
static int
testb2(decomp_t *S, long D, GEN theta, GEN ns)
{
  long v = varn(S->chi), dlim = degpol(S->chi)-1;
  GEN T0 = S->phi, chi0 = S->chi;

  if (DEBUGLEVEL>4) fprintferr("  Increasing Fa\n");
  for (;;)
  {
    S->phi = gadd(theta, FpX_rand(dlim, v, S->p));
    /* phi non-primary ? */
    if (factcp(S, ns) > 1) { composemod(S, S->phi, T0); return 1; }
    if (degpol(S->nu) == D) break;
    S->chi = chi0;
  }
  S->phi0 = T0; return 0; /* F_phi=lcm(F_alpha, F_theta)=D and E_phi=E_alpha */
}

/* return 1 if at least 2 factors mod p ==> chi can be split. 
 * compute a new S->phi such that E = lcm(Ea, Et) */
static int
testc2(decomp_t *S, GEN A, long Ea, GEN T, long Et, GEN ns)
{
  GEN c, T0 = S->phi;
  long r, s, t; 

  if (DEBUGLEVEL>4) fprintferr("  Increasing Ea\n");
  (void)cbezout(Ea, Et, &r, &s); t = 0;
  while (r < 0) { r = r + Et; t++; }
  while (s < 0) { s = s + Ea; t++; }

  c = RgX_mul(RgXQ_u_pow(A, s, S->chi), RgXQ_u_pow(T, r, S->chi));
  c = gdiv(RgX_rem(c, S->chi), powiu(S->p, t));
  S->phi = gadd( pol_x[ varn(S->chi) ], redelt(c, S->pmr, S->p) );
  if (factcp(S, ns) > 1) { composemod(S, S->phi, T0); return 1; }
  S->phi0 = T0; return 0; /* E_phi = lcm(E_alpha,E_theta) */
}

/* used to cache the newton sums of chi */
static GEN
init_NS(long N, GEN pp, GEN pmf, GEN pmr)
{
  GEN q, ns = cgetg(N+1, t_COL);
  long i, l, p = itos_or_0(pp);
  q = p? powiu(pp, (ulong)ceil( N * 1. / (p * (p-1)) )): pp;
  q = sqri(mulii(q, mulii(pmf, powiu(pmr, N))));
  l  = lgefint(q); /* should be more than enough ... */
  for (i = 1; i <= N; i++) gel(ns,i) = cgeti(l);
  kill_cache(ns); return ns;
}

static GEN
ch_var(GEN x, long v)
{
  if (typ(x) == t_POL) { x = shallowcopy(x); setvarn(x, v); }
  return x;
}

/* x p^-eq nu^-er mod p */
static GEN
get_gamma(decomp_t *S, GEN x, long eq, long er)
{
  GEN q, g = x, Dg = powiu(S->p, eq);
  if (er)
  {
    if (!S->invnu)
    {
      while (gdvd(S->chi, S->nu)) S->nu = gadd(S->nu, S->p);
      S->invnu = QXQ_inv(S->nu, S->chi);
      S->invnu = redelt_i(S->invnu, S->pmr, S->p, &(S->Dinvnu));
    }
    if (S->Dinvnu) Dg = mulii(Dg, powiu(S->Dinvnu, er));
    q = mulii(S->p, Dg);
    g = gmul(g, FpXQ_pow(S->invnu, stoi(er), S->chi, q));
    g = FpX_rem(g, S->chi, q);
    update_den(&g, &Dg, NULL);
    g = centermod(g, mulii(S->p, Dg));
  }
  if (!is_pm1(Dg)) g = gdiv(g, Dg);
  return g;
}

/* return 1 if at least 2 factors mod p ==> chi can be split */
static int
loop(decomp_t *S, long nv, long Ea, long Fa, GEN ns)
{
  pari_sp av2 = avma, limit = stack_lim(av2, 1);
  GEN w, chib, beta, gamm, chig, nug, delt = NULL;
  long i, l, Fg, fm = 0, go_fm = 2, eq = 0, er = 0;
  long N = degpol(S->f), v = varn(S->f);

  beta  = FpXQ_pow(S->nu, stoi(Ea), S->chi, S->p);
  S->invnu = NULL;
  chib = chig = NULL; /* -Wall */
  for (;;)
  { /* beta tends to a factor of chi */
    if (DEBUGLEVEL>4) fprintferr("  beta = %Z\n", beta);
    
    if (fm == -1) {
      if (DEBUGLEVEL>4) fprintferr("  ** switching to normal mode\n");
      fm = 0;
      go_fm = eq + 2;
    } else if (!fm && eq > go_fm && !er) {
      if (DEBUGLEVEL>4) fprintferr("  ** switching to fast mode\n");
      fm = 1;
    }

    if (fm)
    {
      er++;
      if (er % Ea == 0)  { er = 0; eq++; }
      gamm = get_gamma(S, beta, eq, er); /* = beta p^-eq  nu^-er (a unit) */
      nug = fastnu(S->p, S->chi, gamm, S->pdr);
      if (!nug) { fm = -1; continue; }
    }
    else
    {
      GEN R = modii(ZX_resultant(beta, S->chi), S->pmf);
      long L, E;
      if (signe(R))
      {
        chib = NULL;
        L = Z_pval(R, S->p);
        E = N;
      }
      else
      { /* pmf | norm(beta) ==> useless */
        chib = ZX_caract(S->chi, beta, v);
        vstar(S->p, chib, &L, &E);
      }
      eq = (long)(L / E);
      er = (long)(L*Ea / E - eq*Ea);
      if (DEBUGLEVEL>4) fprintferr("  (eq,er) = (%ld,%ld)\n", eq,er);
      if (er || !chib)
      { /* gamm might not be an integer ==> chig = NULL */
        gamm = get_gamma(S, beta, eq, er); /* = beta p^-eq  nu^-er (a unit) */
        chig = mycaract(S->chi, gamm, S->p, S->pmr, -1, ns);
      }
      else
      { /* gamm = beta/p^eq, special case of the above */
        GEN h = powiu(S->p, eq);
        gamm = gdiv(beta, h);
        chig = gdiv(RgX_unscale(chib, h), powiu(h, N));
        chig = gcmp1(Q_denom(chig))? FpX_red(chig, S->pmf): NULL;
      }

      if (!chig)
      { /* Valuation of beta was wrong ==> gamma fails the v*-test */
        chib = ZX_caract(S->chi, beta, v);
        vstar(S->p, chib, &L, &E);
        eq = (long)(L / E);
        er = (long)(L*Ea / E - eq*Ea);
      
        gamm = get_gamma(S, beta, eq, er); /* an integer */
        chig = mycaract(S->chi, gamm, S->p, S->pmf, -1, ns);
      }
      
      nug = get_nu(chig, S->p, &l);
      if (l > 1) {
        S->chi = chig;
        S->nu  = nug; composemod(S, gamm, S->phi); return 1;
      }
      
      Fg = degpol(nug);
      if (Fa % Fg) return testb2(S, clcm(Fa,Fg), gamm, ns);
    }

    /* nug irreducible mod p */
    w = FpX_factorff_irred(nug, ch_var(S->nu, nv), S->p);
    if (degpol(w[1]) != 1)
    {
      if (fm) { fm = -1; continue; }
      pari_err(talker, "no root in nilord. Is p = %Z a prime?", S->p);
    }

    for (i = 1; i < lg(w); i++)
    { /* Look for a root delt of nug in Fp[phi] such that vp(gamma - delta) > 0
         Can be used to improve beta */
      GEN eta, chie, nue, W = gel(w,i); /* monic linear polynomial */
      delt = gneg_i( ch_var(gel(W,2), v) );
      eta  = gsub(gamm, delt);	
      if (fm)
      {
        if (fastvalpos(eta, S->chi, S->p, ns, Ea)) break;
        continue;
      }
      
      if (typ(delt) == t_INT)
        chie = translate_pol(chig, delt); /* frequent special case */
      else
      {
        if (!dvdii(ZX_QX_resultant(S->chi, eta), S->p)) continue;
        chie = mycaract(S->chi, eta, S->p, S->pmr, -1, ns);
      }
      nue = get_nu(chie, S->p, &l);
      if (l > 1) { 
        S->nu = nue;
        S->chi= chie; composemod(S, eta, S->phi); return 1;
      }

      if (ismonome(nue))
      { /* vp(eta) = vp(gamma - delta) > 0 */
        long Le, Ee;
        GEN pie;
        if (dvdii(constant_term(chie), S->pmr))
          chie = mycaract(S->chi, eta, S->p, S->pmf, -1, ns);
        
        pie = getprime(S, eta, chie, nue, &Le, &Ee,  0,Ea);
        if (pie) return testc2(S, S->nu, Ea, pie, Ee, ns);
        break;
      }
    }
    if (i == lg(w))
    {
      if (fm) { fm = -1; continue; }
      pari_err(talker, "no root in nilord. Is p = %Z a prime?", S->p);
    }

    if (eq) delt = gmul(delt, powiu(S->p,  eq));
    if (er) delt = gmul(delt, gpowgs(S->nu, er));
    beta = gsub(beta, delt);

    if (low_stack(limit,stack_lim(av2,1)))
    {
      if (DEBUGMEM > 1) pari_warn(warnmem, "nilord");
      gerepileall(av2, S->invnu? 3: 1, &beta, &(S->invnu), &(S->Dinvnu));
    }
  }
}

/* flag != 0 iff we're looking for the p-adic factorization,
   in which case it is the p-adic precision we want */
static GEN
nilord(decomp_t *S, GEN dred, long mf, long flag)
{
  GEN p = S->p;
  long Fa, La, Ea, oE, l, N  = degpol(S->f), v = varn(S->f), nv = fetch_var();
  GEN ns, pia, opa;

  if (DEBUGLEVEL>2)
  {
    fprintferr("  entering Nilord");
    if (DEBUGLEVEL>4)
    {
      fprintferr(" with parameters: %Z^%ld\n", p, S->df);
      fprintferr("  fx = %Z, gx = %Z", S->f, S->nu);
    }
    fprintferr("\n");
  }

  S->pmr = mulii(sqri(dred), p);
  S->pdr = mulii(dred, p);
  S->chi = centermod(S->f, S->pmr);
  S->pmf = powiu(p, mf + 1);
  ns = init_NS(N, p, S->pmf, S->pmr);
  oE = 0;
  opa = NULL;

  for(;;)
  {
    l = 2; /* Decomp by default */
    S->phi0 = NULL; /* no delayed composition */
    Fa   = degpol(S->nu);
    for(;;)
    {
      pia  = getprime(S, pol_x[v], S->chi, S->nu, &La, &Ea, oE,0);
      if (pia) break;
      S->phi = gadd(S->phi, opa);
      S->chi = NULL;
      if (!update_phi(S, ns, &l, flag)) break;
    }
    if (!pia) break;
    oE = Ea; opa = RgX_RgXQ_compo(pia, S->phi, S->f);
    if (La > 1)
    { /* change phi such that nu = pia */
      S->phi = gadd(S->phi, opa);
      S->chi = NULL;
      if (!update_phi(S, ns, &l, flag)) break;
    }

    if (DEBUGLEVEL>5) fprintferr("  (Fa, Ea) = (%ld,%ld)\n", Fa, Ea);
    if (Ea*Fa == N)
    { /* O = Zp[phi] */
      if (!flag) S->phi = redelt(S->phi, sqri(p), p);
      S->chi = NULL; l = 1; break;
    }
    l = 2;
    if (loop(S, nv, Ea, Fa, ns)) break;
    if (!update_phi(S, ns, &l, flag)) break;
  }
  (void)delete_var();
  if (l == 1) return flag? NULL: dbasis(p, S->f, mf, S->phi, S->chi);
  return Decomp(S, flag);
}

/* Assume respm(f,g) divides p^M. Return respm(f, g), using dynamic p-adic
 * precision (until result is non-zero or p^M). */
GEN
fast_respm(GEN f, GEN g, GEN p, long M)
{
  long m = 32 / expi(p); /* p^m ~ 2^32 for initial value of m */
  GEN R, q = NULL;
  if (!m) m = 1;
  for(;; m <<= 1) {
    if (M < 2*m) break;
    q = q? sqri(q): powiu(p, m); /* p^m */
    R = respm(f,g, q); if (signe(R)) return R;
  }
  q = powiu(p, M);
  R = respm(f,g, q); return signe(R)? R: q;
}

GEN
maxord_i(GEN p, GEN f, long mf, GEN w, long flag)
{
  long l = lg(w)-1;
  GEN h = gel(w,l); /* largest factor */
  GEN D = fast_respm(f, derivpol(f), p, mf);
  decomp_t S;

  S.f = f;
  S.p = p;
  S.nu = h;
  S.df = Z_pval(D, p);
  S.phi = pol_x[varn(f)];
  if (l == 1) return nilord(&S, D, mf, flag);
  if (flag && flag <= mf) flag = mf + 1;
  S.chi = f; return Decomp(&S, flag);
}

/* DP = multiple of disc(P) or NULL
 * Return a multiple of the denominator of an algebraic integer (in Q[X]/(P))
 * when expressed in terms of the power basis */
GEN
indexpartial(GEN P, GEN DP)
{
  pari_sp av = avma;
  long i, nb;
  GEN fa, res = gen_1, dP = derivpol(P);
  pari_timer T;

  if(DEBUGLEVEL>=5) (void)TIMER(&T);
  if (!DP) DP = ZX_disc(P);
  DP = mpabs(DP);
  if(DEBUGLEVEL>=5) msgTIMER(&T,"IndexPartial: discriminant");
  fa = auxdecomp(DP, 0);
  if(DEBUGLEVEL>=5) msgTIMER(&T,"IndexPartial: factorization");
  nb = lg(fa[1])-1;
  for (i = 1; i <= nb; i++)
  {
    long E = itos(gmael(fa,2,i)), e = E >> 1;
    GEN p = gmael(fa,1,i), q = p;
    if (i == nb)
      q = powiu(p, (odd(E) && !BSW_psp(p))? e+1: e);
    else if (e >= 2)
    {
      if(DEBUGLEVEL>=5) fprintferr("IndexPartial: factor %Z^%ld ",p,E);
      q = fast_respm(P, dP, p, e);
      if(DEBUGLEVEL>=5) { fprintferr("--> %Z : ",q); msgTIMER(&T,""); }
    }
    res = mulii(res, q);
  }
  return gerepileuptoint(av,res);
}

/*******************************************************************/
/*                                                                 */
/*    2-ELT REPRESENTATION FOR PRIME IDEALS (dividing index)       */
/*                                                                 */
/*******************************************************************/
/* to compute norm of elt in algtobasis form */
typedef struct {
  long r1;
  GEN M;  /* via norm_by_embed */

  GEN D, w, T; /* via resultant if M = NULL */
} norm_S;

static GEN
get_norm(norm_S *S, GEN a)
{
  if (S->M)
  {
    long e;
    GEN N = grndtoi( norm_by_embed(S->r1, gmul(S->M, a)), &e );
    if (e > -5) pari_err(precer, "get_norm");
    return N;
  }
  if (S->w) a = gmul(S->w, a);
  return ZX_resultant_all(S->T, a, S->D, 0);
}

/* q = p^(f+1). a/D in pr | p, norm(pr) = pf.
 * Return 1 if (a/D,p) = pr, and 0 otherwise */
static int
is_uniformizer(GEN a, GEN q, norm_S *S)
{
  return (remii(get_norm(S,a), q) != gen_0);
}

/* return x * y, x, y are t_MAT (Fp-basis of in O_K/p), assume (x,y)=1.
 * x or y may be NULL (= ok), not both */
static GEN
mul_intersect(GEN x, GEN y, GEN p)
{
  if (!x) return y;
  if (!y) return x;
  return FpM_intersect(x, y, p);
}

static GEN
Fp_basis(GEN nf, GEN pr)
{
  long i, j, l;
  GEN x, y;
  if (typ(pr) == t_MAT) return pr;
  x = prime_to_ideal(nf, pr);
  l = lg(x);
  y = cgetg(l, t_MAT);
  for (i=j=1; i<l; i++)
    if (gcmp1(gcoeff(x,i,i))) y[j++] = x[i];
  setlg(y, j); return y;
}

static GEN
get_LV(GEN nf, GEN L, GEN p, long N)
{
  long i, l = lg(L)-1;
  GEN LV, LW, A, B;

  LV = cgetg(l+1, t_VEC);
  if (l == 1) { gel(LV,1) = matid(N); return LV; }
  LW = cgetg(l+1, t_VEC);
  for (i=1; i<=l; i++) gel(LW,i) = Fp_basis(nf, gel(L,i));

  /* A[i] = L[1]...L[i-1], i = 2..l */
  A = cgetg(l+1, t_VEC); A[1] = 0;
  for (i=1; i < l; i++) gel(A,i+1) = mul_intersect(gel(A,i), gel(LW,i), p);
  /* B[i] = L[i+1]...L[l], i = 1..(l-1) */
  B = cgetg(l+1, t_VEC); B[l] = 0;
  for (i=l; i>=2; i--)  gel(B,i-1) = mul_intersect(gel(B,i), gel(LW,i), p);
  for (i=1; i<=l; i++) gel(LV,i) = mul_intersect(gel(A,i), gel(B,i), p);
  return LV;
}

static void
errprime(GEN p) { pari_err(talker, "primedec: %Z is not prime", p); }

/* P = Fp-basis (over O_K/p) for pr.
 * V = Z-basis for I_p/pr. ramif != 0 iff some pr|p is ramified.
 * Return a p-uniformizer for pr. */
static GEN
uniformizer(GEN nf, norm_S *S, GEN P, GEN V, GEN p, int ramif)
{
  long i, l, f, m = lg(P)-1, N = degpol(nf[1]);
  GEN u, Mv, x, q;

  if (!m) return gscalcol_i(p,N);
  /* we want v_p(Norm(x)) = p^f, f = N-m */
  f = N - m;
  q = powiu(p,f+1);

  u = FpM_invimage(shallowconcat(P, V), col_ei(N,1), p);
  setlg(u, lg(P));
  u = centermod(gmul(P, u), p);
  if (is_uniformizer(u, q, S)) return u;
  if (signe(u[1]) <= 0) /* make sure u[1] in ]-p,p] */
    gel(u,1) = addii(gel(u,1), p); /* try u + p */
  else
    gel(u,1) = subii(gel(u,1), p); /* try u - p */
  if (!ramif || is_uniformizer(u, q, S)) return u;

  /* P/p ramified, u in P^2, not in Q for all other Q|p */
  Mv = eltmul_get_table(nf, unnf_minus_x(u));
  l = lg(P);
  for (i=1; i<l; i++)
  {
    x = centermod(gadd(u, gmul(Mv, gel(P,i))), p);
    if (is_uniformizer(x, q, S)) return x;
  }
  errprime(p);
  return NULL; /* not reached */
}

static void
init_norm(norm_S *S, GEN nf, GEN p)
{
  GEN T = gel(nf,1);
  long N = degpol(T);

  S->M = NULL;
  if (typ(nf[5]) == t_VEC) /* beware dummy nf from padicff */
  {
    GEN M = gmael(nf,5,1);
    long ex = gexpo(M) + gexpo(mulsi(8 * N, p));
    if (N * ex <= bit_accuracy(gprecision(M)))
    { /* enough prec to use norm_by_embed */
      S->M = M;
      S->r1 = nf_get_r1(nf);
    }
  }
  if (!S->M)
  {
    GEN D, w = Q_remove_denom(gel(nf,7), &D), Dp = sqri(p);
    long i;
    if (!D) w = shallowcopy(w);
    else {
      GEN w1 = D;
      long v = Z_pval(D, p);
      D = powiu(p, v);
      Dp = mulii(D, Dp);
      gel(w, 1) = remii(w1, Dp);
    }
    for (i=2; i<=N; i++) gel(w,i) = FpX_red(gel(w,i), Dp);
    S->D = D;
    S->w = w;
    S->T = T;
  }
}

/* Assuming P = (p,u) prime, return tau such that p Z + tau Z = p P^(-1)*/
static GEN
anti_uniformizer(GEN nf, GEN p, GEN u)
{
  pari_sp av = avma;
  GEN mat = eltmul_get_table(nf, u);
  return gerepileupto(av, FpM_deplin(mat,p));
}

/*******************************************************************/
/*                                                                 */
/*                   BUCHMANN-LENSTRA ALGORITHM                    */
/*                                                                 */
/*******************************************************************/
static GEN
mk_pr(GEN p, GEN u, long e, long f, GEN t)
{
  GEN pr = cgetg(6, t_VEC);
  gel(pr,1) = p;
  gel(pr,2) = u;
  gel(pr,3) = utoipos(e);
  gel(pr,4) = utoipos(f);
  gel(pr,5) = t; return pr;
}

/* pr = (p,u) of ramification index e */
GEN
primedec_apply_kummer(GEN nf,GEN u,long e,GEN p)
{
  GEN t, T = gel(nf,1);
  long f = degpol(u), N = degpol(T);

  if (f == N) /* inert */
  {
    u = gscalcol_i(p,N);
    t = gscalcol_i(gen_1,N);
  }
  else
  { /* make sure v_pr(u) = 1 (automatic if e>1) */
    t = poltobasis(nf, FpX_div(T,u,p));
    t = centermod(t, p);
    u = FpX_center(u, p);
    if (e == 1)
    {
      norm_S S;
      S.D = S.w = S.M = NULL; S.T = T;
      if (!is_uniformizer(u, powiu(p,f+1), &S)) gel(u,2) = addii(gel(u,2), p);
    }
    u = poltobasis(nf,u);
  }
  return mk_pr(p,u,e,f,t);
}

/* return a Z basis of Z_K's p-radical, phi = x--> x^p-x */
static GEN
pradical(GEN nf, GEN p, GEN *phi)
{
  long i,N = degpol(nf[1]);
  GEN q,m,frob,rad;

  /* matrix of Frob: x->x^p over Z_K/p */
  frob = cgetg(N+1,t_MAT);
  for (i=1; i<=N; i++)
    gel(frob,i) = element_powid_mod_p(nf,i,p, p);

  m = frob; q = p;
  while (cmpiu(q,N) < 0) { q = mulii(q,p); m = FpM_mul(m, frob, p); }
  rad = FpM_ker(m, p); /* m = Frob^k, s.t p^k >= N */
  for (i=1; i<=N; i++)
    gcoeff(frob,i,i) = subis(gcoeff(frob,i,i), 1);
  *phi = frob; return rad;
}

/* return powers of a: a^0, ... , a^d,  d = dim A */
static GEN
get_powers(GEN mul, GEN p)
{
  long i, d = lg(mul[1]);
  GEN z, pow = cgetg(d+2,t_MAT), P = pow+1;

  gel(P,0) = gscalcol_i(gen_1, d-1);
  z = gel(mul,1);
  for (i=1; i<=d; i++)
  {
    gel(P,i) = z; /* a^i */
    if (i!=d) z = FpM_FpC_mul(mul, z, p);
  }
  return pow;
}

/* minimal polynomial of a in A (dim A = d).
 * mul = multiplication table by a in A */
static GEN
pol_min(GEN mul, GEN p)
{
  pari_sp av = avma;
  GEN z, pow = get_powers(mul, p);
  z = FpM_deplin(pow, p);
  if (!z) errprime(p);
  return gerepilecopy(av, RgV_to_RgX(z,0));
}

static GEN
get_pr(GEN nf, norm_S *S, GEN p, GEN P, GEN V, int ramif)
{
  GEN u, t;
  long e, f;

  if (typ(P) == t_VEC) return P; /* already done (Kummer) */

  u = uniformizer(nf, S, P, V, p, ramif);
  t = anti_uniformizer(nf,p,u); if (!t) errprime(p);
  e = ramif? 1 + int_elt_val(nf,t,p,t,NULL): 1;
  f = degpol(nf[1]) - (lg(P)-1);
  return mk_pr(p,u,e,f,t);
}

/* prime ideal decomposition of p */
static GEN
_primedec(GEN nf, GEN p)
{
  GEN E, F, L, Ip, H, phi, mat1, f, g, h, p1, UN, T = gel(nf,1);
  long i, k, c, iL, N;

  F = FpX_factor(T, p);
  E = gel(F,2);
  F = gel(F,1);

  k = lg(F); if (k == 1) errprime(p);
  if (signe(modii(gel(nf,4),p))) /* p doesn't divide index */
  {
    L = cgetg(k,t_VEC);
    for (i=1; i<k; i++)
      gel(L,i) = primedec_apply_kummer(nf,gel(F,i), E[i],p);
    return L;
  }

  g = FpXV_prod(F, p);
  h = FpX_div(T,g,p);
  f = FpX_red(gdivexact(gsub(gmul(g,h), T), p), p);

  N = degpol(T);
  L = cgetg(N+1,t_VEC); iL = 1;
  for (i=1; i<k; i++)
    if (E[i] == 1 || signe(FpX_rem(f,gel(F,i),p)))
      gel(L,iL++) = primedec_apply_kummer(nf,gel(F,i), E[i],p);
    else /* F[i] | (f,g,h), happens at least once by Dedekind criterion */
      E[i] = 0;

  /* phi matrix of x -> x^p - x in algebra Z_K/p */
  Ip = pradical(nf,p,&phi);

  /* split etale algebra Z_K / (p,Ip) */
  h = cgetg(N+1,t_VEC);
  if (iL > 1)
  { /* split off Kummer factors */
    GEN mulbeta, beta = NULL;
    for (i=1; i<k; i++)
      if (!E[i]) beta = beta? FpX_mul(beta, gel(F,i), p): gel(F,i);
    if (!beta) errprime(p);
    beta = FpC_red(poltobasis(nf,beta), p);

    mulbeta = FpM_red(eltmul_get_table(nf, beta), p);
    p1 = shallowconcat(mulbeta, Ip);
    /* Fp-base of ideal (Ip, beta) in ZK/p */
    gel(h,1) = FpM_image(p1, p);
  }
  else
    gel(h,1) = Ip;

  UN = col_ei(N, 1);
  for (c=1; c; c--)
  { /* Let A:= (Z_K/p) / Ip; try to split A2 := A / Im H ~ Im M2
       H * ? + M2 * Mi2 = Id_N ==> M2 * Mi2 projector A --> A2 */
    GEN M, Mi, M2, Mi2, phi2;
    long dim;

    H = gel(h,c); k = lg(H)-1;
    M   = FpM_suppl(shallowconcat(H,UN), p);
    Mi  = FpM_inv(M, p);
    M2  = vecslice(M, k+1,N); /* M = (H|M2) invertible */
    Mi2 = rowslice(Mi,k+1,N);
    /* FIXME: FpM_mul(,M2) could be done with vecpermute */
    phi2 = FpM_mul(Mi2, FpM_mul(phi,M2, p), p);
    mat1 = FpM_ker(phi2, p);
    dim = lg(mat1)-1; /* A2 product of 'dim' fields */
    if (dim > 1)
    { /* phi2 v = 0 <==> a = M2 v in Ker phi */
      GEN I, R, r, a, mula, mul2, v = gel(mat1,2);
      long n;

      a = FpM_FpC_mul(M2,v, p);
      mula = FpM_red(eltmul_get_table(nf, a), p);
      mul2 = FpM_mul(Mi2, FpM_mul(mula,M2, p), p);
      R = FpX_roots(pol_min(mul2,p), p); /* totally split mod p */

      n = lg(R)-1;
      for (i=1; i<=n; i++)
      {
        r = lift_intern(gel(R,i));
        I = gaddmat_i(negi(r), mula);
	gel(h,c++) = FpM_image(shallowconcat(H, I), p);
      }
      if (n == dim)
        for (i=1; i<=n; i++) { H = gel(h,--c); gel(L,iL++) = H; }
    }
    else /* A2 field ==> H maximal, f = N-k = dim(A2) */
      gel(L,iL++) = H;
  }
  setlg(L, iL);
{
  GEN Lpr = cgetg(iL, t_VEC);
  GEN LV = get_LV(nf, L,p,N);
  int ramif = dvdii(gel(nf,3), p);
  norm_S S; init_norm(&S, nf, p);
  for (i=1; i<iL; i++)
    gel(Lpr,i) = get_pr(nf, &S, p, gel(L,i), gel(LV,i), ramif);
  return Lpr;
}
}

GEN
primedec(GEN nf, GEN p)
{
  pari_sp av = avma;
  if (typ(p) != t_INT) pari_err(typeer, "primedec");
  return gerepileupto(av, gen_sort(_primedec(checknf(nf),p),
                                   0, cmp_prime_over_p));
}

/* return [Fp[x]: Fp] */
static long
ffdegree(GEN x, GEN frob, GEN p)
{
  pari_sp av = avma;
  long d, f = lg(frob)-1;
  GEN y = x;

  for (d=1; d < f; d++)
  {
    y = FpM_FpC_mul(frob, y, p);
    if (gequal(y, x)) break;
  }
  avma = av; return d;
}

static GEN
lift_to_zk(GEN v, GEN c, long N)
{
  GEN w = zerocol(N);
  long i, l = lg(c);
  for (i=1; i<l; i++) w[c[i]] = v[i];
  return w;
}

/* return integral x = 0 mod p/pr^e, (x,pr) = 1.
 * Don't reduce mod p here: caller may need result mod pr^k */
GEN
special_anti_uniformizer(GEN nf, GEN pr)
{
  GEN p = gel(pr,1), e = gel(pr,3);
  return gdivexact(element_pow(nf,gel(pr,5),e), powiu(p, e[2]-1));
}

/* return t = 1 mod pr, t = 0 mod p / pr^e(pr/p) */
static GEN
anti_uniformizer2(GEN nf, GEN pr)
{
  GEN p = gel(pr,1), z;
  z = FpC_red(special_anti_uniformizer(nf, pr), p);
  z = hnfmodid(eltmul_get_table(nf, z), p);
  z = idealaddtoone_i(nf, pr, z);
  return unnf_minus_x(z);
}

#define mpr_TAU 1
#define mpr_FFP 2
#define mpr_PR  3
#define mpr_T   4
#define mpr_NFP 5
#define SMALLMODPR 4
#define LARGEMODPR 6
static GEN
modpr_TAU(GEN modpr)
{
  GEN tau = gel(modpr,mpr_TAU);
  if (typ(tau) == t_INT && signe(tau) == 0) tau = NULL;
  return tau;
}

/* prh = HNF matrix, which is identity but for the first line. Return a
 * projector to ZK / prh ~ Z/prh[1,1] */
GEN
dim1proj(GEN prh)
{
  long i, N = lg(prh)-1;
  GEN ffproj = cgetg(N+1, t_VEC);
  GEN x, q = gcoeff(prh,1,1);
  gel(ffproj,1) = gen_1;
  for (i=2; i<=N; i++)
  {
    x = gcoeff(prh,1,i);
    if (signe(x)) x = subii(q,x);
    gel(ffproj,i) = x;
  }
  return ffproj;
}

/* p not necessarily prime, but coprime to denom(basis) */
GEN
get_proj_modT(GEN basis, GEN T, GEN p)
{
  long i, l = lg(basis), f = degpol(T);
  GEN z = cgetg(l, t_MAT);
  for (i = 1; i < l; i++)
  {
    GEN cx, w = gel(basis,i);
    if (typ(w) != t_INT)
    {
      w = Q_primitive_part(w, &cx);
      w = FpX_rem(w, T, p);
      if (cx) w = FpX_Fp_mul(w, Rg_to_Fp(cx, p), p);
    }
    gel(z,i) = RgX_to_RgV(w, f); /* w_i mod (T,p) */
  }
  return z;
}

static GEN
modprinit(GEN nf, GEN pr, int zk)
{
  pari_sp av = avma;
  GEN res, tau, mul, x, p, T, pow, ffproj, nfproj, prh, c;
  long N, i, k, f;

  nf = checknf(nf); checkprimeid(pr);
  f = itos( gel(pr,4) );
  N = degpol(nf[1]);
  prh = prime_to_ideal(nf, pr);
  tau = zk? gen_0: anti_uniformizer2(nf, pr);
  p = gel(pr,1);

  if (f == 1)
  {
    res = cgetg(SMALLMODPR, t_COL);
    gel(res,mpr_TAU) = tau;
    gel(res,mpr_FFP) = dim1proj(prh);
    gel(res,mpr_PR) = pr; return gerepilecopy(av, res);
  }

  c = cgetg(f+1, t_VECSMALL);
  ffproj = cgetg(N+1, t_MAT);
  for (k=i=1; i<=N; i++)
  {
    x = gcoeff(prh, i,i);
    if (!is_pm1(x)) { c[k] = i; gel(ffproj,i) = col_ei(N, i); k++; }
    else
      gel(ffproj,i) = gneg(gel(prh,i));
  }
  ffproj = rowpermute(ffproj, c);
  if (! dvdii(gel(nf,4), p))
  {
    GEN basis = gel(nf,7);
    if (N == f) T = gel(nf,1); /* pr inert */
    else
    {
      T = Q_primpart(gmul(basis, gel(pr,2)));
      basis = vecpermute(basis, c);
    }
    T = FpX_red(T, p);
    ffproj = FpM_mul(get_proj_modT(basis, T, p), ffproj, p);

    res = cgetg(SMALLMODPR+1, t_COL);
    gel(res,mpr_TAU) = tau;
    gel(res,mpr_FFP) = ffproj;
    gel(res,mpr_PR) = pr;
    gel(res,mpr_T) = T; return gerepilecopy(av, res);
  }

  if (uisprime(f))
  {
    mul = eltmulid_get_table(nf, c[2]);
    mul = vecpermute(mul, c);
  }
  else
  {
    GEN v, u, u2, frob;
    long deg,deg1,deg2;

    /* matrix of Frob: x->x^p over Z_K/pr = < w[c1], ..., w[cf] > over Fp */
    frob = cgetg(f+1, t_MAT);
    for (i=1; i<=f; i++)
    {
      x = element_powid_mod_p(nf,c[i],p, p);
      gel(frob,i) = FpM_FpC_mul(ffproj, x, p);
    }
    u = col_ei(f,2); k = 2;
    deg1 = ffdegree(u, frob, p);
    while (deg1 < f)
    {
      k++; u2 = col_ei(f, k);
      deg2 = ffdegree(u2, frob, p);
      deg = clcm(deg1,deg2);
      if (deg == deg1) continue;
      if (deg == deg2) { deg1 = deg2; u = u2; continue; }
      u = gadd(u, u2);
      while (ffdegree(u, frob, p) < deg) u = gadd(u, u2);
      deg1 = deg;
    }
    v = lift_to_zk(u,c,N);

    mul = cgetg(f+1,t_MAT);
    gel(mul,1) = v; /* assume w_1 = 1 */
    for (i=2; i<=f; i++) gel(mul,i) = element_mulid(nf,v,c[i]);
  }

  /* Z_K/pr = Fp(v), mul = mul by v */
  mul = FpM_red(mul, p);
  mul = FpM_mul(ffproj, mul, p);

  pow = get_powers(mul, p);
  T = RgV_to_RgX(FpM_deplin(pow, p), varn(nf[1]));
  nfproj = cgetg(f+1, t_MAT);
  for (i=1; i<=f; i++) gel(nfproj,i) = lift_to_zk(gel(pow,i), c, N);
  nfproj = coltoliftalg(nf, nfproj);

  setlg(pow, f+1);
  ffproj = FpM_mul(FpM_inv(pow, p), ffproj, p);

  res = cgetg(LARGEMODPR, t_COL);
  gel(res,mpr_TAU) = tau;
  gel(res,mpr_FFP) = ffproj;
  gel(res,mpr_PR) = pr;
  gel(res,mpr_T) = T;
  gel(res,mpr_NFP) = nfproj; return gerepilecopy(av, res);
}

GEN
nfmodprinit(GEN nf, GEN pr) { return modprinit(nf, pr, 0); }
GEN
zkmodprinit(GEN nf, GEN pr) { return modprinit(nf, pr, 1); }

void
checkmodpr(GEN modpr)
{
  if (typ(modpr) != t_COL || lg(modpr) < SMALLMODPR)
    pari_err(talker,"incorrect modpr format");
  checkprimeid(gel(modpr,mpr_PR));
}


static GEN
to_ff_init(GEN nf, GEN *pr, GEN *T, GEN *p, int zk)
{
  GEN modpr = (typ(*pr) == t_COL)? *pr: modprinit(nf, *pr, zk);
  *T = lg(modpr)==SMALLMODPR? NULL: gel(modpr,mpr_T);
  *pr = gel(modpr,mpr_PR);
  *p = gel(*pr,1); return modpr;
}
GEN
nf_to_ff_init(GEN nf, GEN *pr, GEN *T, GEN *p) {
  GEN modpr = to_ff_init(nf,pr,T,p,0);
  GEN tau = modpr_TAU(modpr);
  if (!tau) gel(modpr,mpr_TAU) = anti_uniformizer2(nf, *pr);
  return modpr;
}
GEN
zk_to_ff_init(GEN nf, GEN *pr, GEN *T, GEN *p) {
  return to_ff_init(nf,pr,T,p,1);
}

/* assume x in 'basis' form (t_COL) */
GEN
zk_to_ff(GEN x, GEN modpr)
{
  GEN pr = gel(modpr,mpr_PR);
  GEN p = gel(pr,1);
  GEN y = gmul(gel(modpr,mpr_FFP), x);
  if (lg(modpr) == SMALLMODPR) return modii(y,p);
  y = FpC_red(y, p);
  return col_to_ff(y, varn(modpr[mpr_T]));
}

/* REDUCTION Modulo a prime ideal */

/* assume x in t_COL form, v_pr(x) >= 0 */
static GEN
kill_denom(GEN x, GEN nf, GEN p, GEN modpr)
{
  GEN cx, den = denom(x);
  long v;
  if (gcmp1(den)) return x;

  v = Z_pval(den,p);
  if (v)
  {
    GEN tau = modpr_TAU(modpr);
    if (!tau) pari_err(talker,"modpr initialized for integers only!");
    x = element_mul(nf,x, element_pow(nf, tau, utoipos(v)));
  }
  x = Q_primitive_part(x, &cx);
  if (cx) x = gmul(Rg_to_Fp(cx, p), x);
  return FpC_red(x, p);
}

/* x integral, reduce mod prh in HNF */
GEN
nfreducemodpr_i(GEN x, GEN prh)
{
  GEN p = gcoeff(prh,1,1);
  long i,j;

  x = shallowcopy(x);
  for (i=lg(x)-1; i>=2; i--)
  {
    GEN t = gel(prh,i), p1 = remii(gel(x,i), p);
    gel(x,i) = p1;
    if (signe(p1) && is_pm1(t[i]))
    {
      for (j=1; j<i; j++)
        gel(x,j) = subii(gel(x,j), mulii(p1, gel(t,j)));
      gel(x,i) = gen_0;
    }
  }
  gel(x,1) = remii(gel(x,1), p); return x;
}

GEN
nfreducemodpr(GEN nf, GEN x, GEN modpr)
{
  pari_sp av = avma;
  GEN pr, p;

  nf = checknf(nf);
  checkmodpr(modpr);
  pr = gel(modpr,mpr_PR);
  p = gel(pr,1);
  x = algtobasis_i(nf,x);
  x = kill_denom(x, nf, p, modpr);
  x = ff_to_nf(zk_to_ff(x,modpr), modpr);
  return gerepilecopy(av, algtobasis_i(nf,x));
}

GEN
nf_to_ff(GEN nf, GEN x, GEN modpr)
{
  pari_sp av = avma;
  GEN pr = gel(modpr,mpr_PR);
  GEN p = gel(pr,1);
  long t = typ(x);

  if (t == t_POLMOD) { x = gel(x,2); t = typ(x); }
  nf = checknf(nf);
  switch(t)
  {
    case t_INT: return modii(x, p);
    case t_FRAC: return Rg_to_Fp(x, p);
    case t_POL: x = poltobasis(nf, x); break;
    case t_COL: break;
    default: pari_err(typeer,"nf_to_ff");
  }
  x = kill_denom(x, nf, p, modpr);
  return gerepilecopy(av, zk_to_ff(x, modpr));
}

GEN
ff_to_nf(GEN x, GEN modpr)
{
  if (lg(modpr) < LARGEMODPR) return x;
  return mulmat_pol(gel(modpr,mpr_NFP), x);
}
GEN
modprM_lift(GEN x, GEN modpr)
{
  long i,j,h,l = lg(x);
  GEN y = cgetg(l, t_MAT);
  if (l == 1) return y;

  h = lg(x[1]);
  for (j=1; j<l; j++)
  {
    GEN p1 = cgetg(h,t_COL); gel(y,j) = p1;
    for (i=1; i<h; i++) gel(p1,i) = ff_to_nf(gcoeff(x,i,j), modpr);
  }
  return y;
}
GEN
modprX_lift(GEN x, GEN modpr)
{
  long i, l;
  GEN z;

  if (typ(x)!=t_POL) return gcopy(x); /* scalar */
  l = lg(x); z = cgetg(l, t_POL); z[1] = x[1];
  for (i=2; i<l; i++) gel(z,i) = ff_to_nf(gel(x,i), modpr);
  return z;
}

/* reduce the coefficients of pol modulo modpr */
GEN
modprX(GEN x, GEN nf,GEN modpr)
{
  long i, l;
  GEN z;

  if (typ(x)!=t_POL) return nf_to_ff(nf,x,modpr);
  l = lg(x); z = cgetg(l,t_POL); z[1] = x[1];
  for (i=2; i<l; i++) gel(z,i) = nf_to_ff(nf,gel(x,i),modpr);
  return normalizepol(z);
}
GEN
modprV(GEN z, GEN nf,GEN modpr)
{
  long i,l = lg(z);
  GEN x;
  x = cgetg(l,typ(z));
  for (i=1; i<l; i++) gel(x,i) = nf_to_ff(nf,gel(z,i), modpr);
  return x;
}
/* assume z a t_VEC/t_COL/t_MAT */
GEN
modprM(GEN z, GEN nf,GEN modpr)
{
  long i,l = lg(z);
  GEN x;

  if (typ(z) != t_MAT) return modprV(z,nf,modpr);
  x = cgetg(l,t_MAT); if (l==1) return x;
  for (i=1; i<l; i++) gel(x,i) = modprV(gel(z,i),nf,modpr);
  return x;
}

/*******************************************************************/
/*                                                                 */
/*                       RELATIVE ROUND 2                          */
/*                                                                 */
/*******************************************************************/

static void
fill(long l, GEN H, GEN Hx, GEN I, GEN Ix)
{
  long i;
  if (typ(Ix) == t_VEC) /* standard */
    for (i=1; i<l; i++) { H[i] = Hx[i]; I[i] = Ix[i]; }
  else /* constant ideal */
    for (i=1; i<l; i++) { H[i] = Hx[i]; gel(I,i) = Ix; }
}

/* given MODULES x and y by their pseudo-bases, returns a pseudo-basis of the
 * module generated by x and y. */
static GEN
rnfjoinmodules_i(GEN nf, GEN Hx, GEN Ix, GEN Hy, GEN Iy)
{
  long lx = lg(Hx), ly = lg(Hy), l = lx+ly-1;
  GEN H, I;

  H = cgetg(l, t_MAT);
  I = cgetg(l, t_VEC);
  fill(lx, H     , Hx, I     , Ix);
  fill(ly, H+lx-1, Hy, I+lx-1, Iy); return nfhermite(nf, mkvec2(H, I));
}
static GEN
rnfjoinmodules(GEN nf, GEN x, GEN y)
{
  if (!x) return y;
  if (!y) return x;
  return rnfjoinmodules_i(nf, gel(x,1), gel(x,2), gel(y,1), gel(y,2));
}

typedef struct {
  GEN nf, multab, modpr,T,p;
  long h;
} rnfeltmod_muldata;

static GEN
_mul(void *data, GEN x, GEN y/* base; ignored */)
{
  rnfeltmod_muldata *D = (rnfeltmod_muldata *) data;
  GEN z = x? element_mulid(D->multab,x,D->h)
           : element_mulidid(D->multab,D->h,D->h);
  (void)y;
  return FqV_red(z,D->T,D->p);
}
static GEN
_sqr(void *data, GEN x)
{
  rnfeltmod_muldata *D = (rnfeltmod_muldata *) data;
  GEN z = x? sqr_by_tab(D->multab,x)
           : element_mulidid(D->multab,D->h,D->h);
  return FqV_red(z,D->T,D->p);
}

/* Compute W[h]^n mod pr in the extension, assume n >= 0 */
static GEN
rnfelementid_powmod(GEN multab, long h, GEN n, GEN T, GEN p)
{
  pari_sp av = avma;
  GEN y;
  rnfeltmod_muldata D;

  if (!signe(n)) return gen_1;

  D.multab = multab;
  D.h = h;
  D.T = T;
  D.p = p;
  y = leftright_pow(NULL, n, (void*)&D, &_sqr, &_mul);
  return gerepilecopy(av, y);
}

/* Relative Dedekind criterion over nf, applied to the order defined by a
 * root of irreducible polynomial P, modulo the prime ideal pr. Assume
 * vdisc = v_pr( disc(P) ).
 * Return NULL if nf[X]/P is pr-maximal. Otherwise, return [flag, O, v]:
 *   O = enlarged order, given by a pseudo-basis
 *   flag = 1 iff O is pr-maximal
 *   v = v_pr(disc(O)). */
static GEN
rnfdedekind_i(GEN nf, GEN P, GEN pr, long vdisc)
{
  long vt, r, d, n, m, i, j;
  pari_sp av = avma;
  GEN Prd, A, I, p, tau, g, matid;
  GEN modpr, h, k, base, nfT, T, gzk, hzk, prinvp, X, pal;

  P = lift(P);
  nf = checknf(nf); nfT = gel(nf,1);
  modpr = nf_to_ff_init(nf,&pr, &T, &p);
  tau = coltoliftalg(nf, gel(pr,5));
  n = degpol(nfT);
  m = degpol(P);

  Prd = modprX(P, nf, modpr);
  A = (GEN)FqX_factor(Prd,T,p)[1];
  r = lg(A); if (r < 2) pari_err(constpoler,"rnfdedekind");
  g = gel(A,1);
  for (i=2; i<r; i++) g = FqX_mul(g, gel(A,i), T, p);
  h = FqX_div(Prd,g, T, p);
  gzk = modprX_lift(g, modpr);
  hzk = modprX_lift(h, modpr);

  k = gsub(P, RgXQX_mul(gzk,hzk, nfT));
  k = gdiv(RgXQX_RgXQ_mul(k, tau, nfT), p);
  k = modprX(k, nf, modpr);
  k  = FqX_gcd(FqX_gcd(g,h,  T,p), k, T,p);
  d = degpol(k);  /* <= m */
  if (!d) return NULL; /* pr-maximal */

  A = cgetg(m+d+1,t_MAT);
  I = cgetg(m+d+1,t_VEC); base = mkvec2(A, I);
 /* base[2] temporarily multiplied by p, for the final nfhermitemod,
  * which requires integral ideals */
  matid = gscalmat(p, n);
  prinvp = pidealprimeinv(nf,pr); /* again multiplied by p */
  for (j=1; j<=m; j++)
  {
    gel(A,j) = col_ei(m, j);
    gel(I,j) = matid;
  }
  X = pol_x[varn(P)];
  pal = FqX_div(Prd,k, T,p);
  pal = modprX_lift(pal, modpr);
  for (   ; j<=m+d; j++)
  {
    gel(A,j) = RgX_to_RgV(pal,m);
    gel(I,j) = prinvp;
    pal = RgXQX_rem(RgXQX_mul(pal,X,nfT),P,nfT);
  }
  /* the modulus is integral */
  base = nfhermitemod(nf,base, gmul(powiu(p, m-d),
                                    idealpows(nf, prinvp, d)));
  gel(base,2) = gdiv(gel(base,2), p); /* cancel the factor p */
  vt = vdisc - 2*d;
  return gerepilecopy(av, mkvec3(vt < 2? gen_1: gen_0, base, stoi(vt)));
}

/* [L:K] = n, [K:Q] = m */
static GEN
triv_order(long n, long m)
{
  GEN I, z = cgetg(3, t_VEC), id = matid(m);
  long i;
  I = cgetg(n+1,t_VEC); for (i=1; i<=n; i++) gel(I,i) = id;
  gel(z,1) = matid(n);
  gel(z,2) = I; return z;
}

/* FIXME: is it really necessary to export this routine ? */
GEN
rnfdedekind(GEN nf, GEN P, GEN pr)
{
  pari_sp av = avma;
  long v = element_val(nf, discsr(P), pr);
  GEN z;
  avma = av; z = rnfdedekind_i(nf, P, pr, v);
  if (!z) {
    z = cgetg(4, t_VEC);
    gel(z,1) = gen_1;
    gel(z,2) = triv_order(degpol(P), degpol(nf[1]));
    gel(z,3) = stoi(v);
  }
  return z;
}

/* return NULL if power order if pr-maximal */
static GEN
rnfordmax(GEN nf, GEN pol, GEN pr, long vdisc)
{
  pari_sp av = avma, av1, lim;
  long i, j, k, n, vpol, cmpt, sep;
  GEN q, q1, p, T, modpr, W, I, MW, C, p1;
  GEN Tauinv, Tau, prhinv, pip, nfT, id, rnfId;

  if (DEBUGLEVEL>1) fprintferr(" treating %Z\n",pr);
  modpr = nf_to_ff_init(nf,&pr,&T,&p);
  p1 = rnfdedekind_i(nf, pol, modpr, vdisc);
  if (!p1) { avma = av; return NULL; }
  if (gcmp1(gel(p1,1))) return gerepilecopy(av,gel(p1,2));
  sep = itos(gel(p1,3));
  W = gmael(p1,2,1);
  I = gmael(p1,2,2);

  pip = coltoalg(nf, gel(pr,2));
  nfT = gel(nf,1);
  n = degpol(pol); vpol = varn(pol);
  q = T? powiu(p,degpol(T)): p;
  q1 = q; while (cmpiu(q1,n) < 0) q1 = mulii(q1,q);
  rnfId = matid(n);
  id    = matid(degpol(nfT));

  prhinv = idealinv(nf, pr);
  C = cgetg(n+1, t_MAT);
  for (j=1; j<=n; j++) gel(C,j) = cgetg(n*n+1, t_COL);
  MW = cgetg(n*n+1, t_MAT);
  for (j=1; j<=n*n; j++) gel(MW,j) = cgetg(n+1, t_COL);
  Tauinv = cgetg(n+1, t_VEC);
  Tau    = cgetg(n+1, t_VEC);
  av1 = avma; lim = stack_lim(av1,1);
  for(cmpt=1; ; cmpt++)
  {
    GEN I0 = shallowcopy(I), W0 = shallowcopy(W);
    GEN Wa, Wainv, Waa, Ip, A, Ainv, MWmod, F, pseudo, G;

    if (DEBUGLEVEL>1) fprintferr("    pass no %ld\n",cmpt);
    for (j=1; j<=n; j++)
    {
      GEN tau, tauinv;
      long v1, v2;
      if (gequal(gel(I,j),id)) { gel(Tau,j) = gel(Tauinv,j) = gen_1; continue; }

      p1 = ideal_two_elt(nf,gel(I,j));
      v1 = element_val(nf,gel(p1,1),pr);
      v2 = element_val(nf,gel(p1,2),pr);
      tau = (v1 > v2)? gel(p1,2): gel(p1,1);
      tauinv = element_inv(nf, tau);
      gel(Tau,j) = tau;
      gel(Tauinv,j) = tauinv;
      gel(W,j) = element_mulvec(nf, tau, gel(W,j));
      gel(I,j) = idealmul(nf,    tauinv, gel(I,j));
    }
    /* W = (Z_K/pr)-basis of O/pr. O = (W0,I0) ~ (W, I) */
    Wa = matbasistoalg(nf,W);

   /* compute MW: W_i*W_j = sum MW_k,(i,j) W_k */
    Waa = lift_intern(RgM_to_RgXV(Wa,vpol));
    Wainv = lift_intern(ginv(Wa));
    for (i=1; i<=n; i++)
      for (j=i; j<=n; j++)
      {
        GEN z = RgXQX_rem(gmul(gel(Waa,i),gel(Waa,j)), pol, nfT);
        long tz = typ(z);
          if (is_scalar_t(tz) || (tz == t_POL && varncmp(varn(z), vpol) > 0))
          z = gmul(z, gel(Wainv,1));
        else
          z = mulmat_pol(Wainv, z);
        for (k=1; k<=n; k++)
        {
          GEN c = grem(gel(z,k), nfT);
          gcoeff(MW, k, (i-1)*n+j) = c;
          gcoeff(MW, k, (j-1)*n+i) = c;
        }
      }

    /* compute Ip =  pr-radical [ could use Ker(trace) if q large ] */
    MWmod = modprM(MW,nf,modpr);
    F = cgetg(n+1, t_MAT); F[1] = rnfId[1];
    for (j=2; j<=n; j++) gel(F,j) = rnfelementid_powmod(MWmod, j, q1, T,p);
    Ip = FqM_ker(F,T,p);
    if (lg(Ip) == 1) { W = W0; I = I0; break; }

    /* Fill C: W_k A_j = sum_i C_(i,j),k A_i */
    A = modprM_lift(FqM_suppl(Ip,T,p), modpr);
    for (j=1; j<lg(Ip); j++)
    {
      p1 = gel(A,j);
      for (i=1; i<=n; i++) gel(p1,i) = mkpolmod(gel(p1,i), nfT);
    }
    for (   ; j<=n; j++)
    {
      p1 = gel(A,j);
      for (i=1; i<=n; i++) gel(p1,i) = gmul(pip, gel(p1,i));
    }
    Ainv = lift_intern(ginv(A));
    A    = lift_intern(A);
    for (k=1; k<=n; k++)
      for (j=1; j<=n; j++)
      {
        GEN z = gmul(Ainv, gmod(element_mulid(MW, gel(A,j),k), nfT));
        for (i=1; i<=n; i++)
        {
          GEN c = grem(gel(z,i), nfT);
          gcoeff(C, (j-1)*n+i, k) = nf_to_ff(nf,c,modpr);
        }
      }
    G = modprM_lift(FqM_ker(C,T,p), modpr);

    pseudo = rnfjoinmodules_i(nf, G,prhinv, rnfId,I);
    /* express W in terms of the power basis */
    W = algtobasis(nf, gmul(Wa, matbasistoalg(nf,gel(pseudo,1))));
    I = gel(pseudo,2);
    /* restore the HNF property W[i,i] = 1. NB: Wa upper triangular, with
     * Wa[i,i] = Tau[i] */
    for (j=1; j<=n; j++)
      if (gel(Tau,j) != gen_1)
      {
        gel(W,j) = element_mulvec(nf, gel(Tauinv,j), gel(W,j));
        gel(I,j) = idealmul(nf,       gel(Tau,j),    gel(I,j) );
      }
    if (DEBUGLEVEL>3) fprintferr(" new order:\n%Z\n%Z\n", W, I);
    if (sep <= 3 || gequal(I,I0)) break;

    if (low_stack(lim, stack_lim(av1,1)) || (cmpt & 3) == 0)
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"rnfordmax");
      gerepileall(av1,2, &W,&I);
    }
  }
  return gerepilecopy(av, mkvec2(W, I));
}

static void
check_pol(GEN *px, long v)
{
  GEN x = *px;
  long i, lx = lg(x);
  for (i=2; i<lx; i++)
  {
    long tx = typ(x[i]);
    if (!is_rational_t(tx)) pari_err(talker,"incorrect coeff in rnf function");
  }
  if (lx == 2) *px = gen_0;
  if (lx == 3) *px = gel(x,2);
}

GEN
fix_relative_pol(GEN nf, GEN x, int chk_lead)
{
  GEN xnf = (typ(nf) == t_POL)? nf: gel(nf,1);
  long i, vnf = varn(xnf), lx = lg(x);
  if (typ(x) != t_POL || varncmp(varn(x), vnf) >= 0)
    pari_err(talker,"incorrect polynomial in rnf function");
  x = shallowcopy(x);
  for (i=2; i<lx; i++)
  {
    GEN c = gel(x,i);
    switch(typ(c))
    {
      case t_INT: case t_FRAC: break;
      case t_POL:
        check_pol(&c, vnf);
        if (typ(c) == t_POL) c = gmodulo(c, xnf);
        break;
      case t_POLMOD:
        if (!gequal(gel(c,1), xnf)) pari_err(consister,"rnf function");
        break;
      default: pari_err(typeer, "rnf function");
    }
    gel(x,i) = c;
  }

  if (chk_lead && !gcmp1(leading_term(x)))
    pari_err(impl,"non-monic relative polynomials");
  return normalizepol_i(x, lx);
}

/* determinant of the trace pairing */
static GEN
get_d(GEN nf, GEN pol, GEN A)
{
  long i, j, n = degpol(pol);
  GEN W = RgM_to_RgXV(lift_intern(matbasistoalg(nf,A)), varn(pol));
  GEN T, nfT = gel(nf,1), sym = polsym_gen(pol, NULL, n-1, nfT, NULL);
  T = cgetg(n+1,t_MAT);
  for (j=1; j<=n; j++) gel(T,j) = cgetg(n+1,t_COL);
  for (j=1; j<=n; j++)
    for (i=j; i<=n; i++)
    {
      GEN c = RgXQX_mul(gel(W,i),gel(W,j), nfT);
      c = RgXQX_rem(c, pol, nfT);
      c = simplify_i(quicktrace(c,sym));
      gcoeff(T,j,i) = gcoeff(T,i,j) = c;
    }
  return algtobasis_i(nf, det(T));
}

/* nf = base field K
 * pol= monic polynomial, coefficients in Z_K, defining a relative
 *   extension L = K[X]/(pol). One MUST have varn(pol) << varn(nf[1]).
 * Returns a pseudo-basis [A,I] of Z_L, set (D,d) to the relative
 * discriminant, and f to the index-ideal */
GEN
rnfallbase(GEN nf, GEN pol, GEN *pD, GEN *pd, GEN *pf)
{
  long i, n, N, l, *ep;
  GEN p1, A, nfT, P, id, I, z, d, D, disc;

  nf = checknf(nf); nfT = gel(nf,1);
  pol = fix_relative_pol(nf,pol,1);
  N = degpol(nfT);
  n = degpol(pol);
  disc = discsr(pol); pol = lift(pol);
  P = idealfactor(nf, disc);
  ep= (long*)P[2];
  P = gel(P,1); l = lg(P);
  for (i=1; i < l; i++) ep[i] = itos(gel(ep,i));
  if (DEBUGLEVEL>1)
  {
    fprintferr("Ideals to consider:\n");
    for (i=1; i < l; i++)
      if (ep[i]>1) fprintferr("%Z^%ld\n",P[i],ep[i]);
    flusherr();
  }
  id = matid(N); z = NULL;
  for (i=1; i < l; i++)
    if (ep[i] > 1)
    {
      GEN y = rnfordmax(nf, pol, gel(P,i), ep[i]);
      z = rnfjoinmodules(nf, z, y);
    }
  if (!z) z = triv_order(n, N);
  A = gel(z,1);
  I = gel(z,2);
  d = get_d(nf, pol, A);

  i=1; while (i<=n && gequal(gel(I,i), id)) i++;
  if (i > n) { D = gen_1; if (pf) *pf = gen_1; }
  else
  {
    D = gel(I,i);
    for (i++; i<=n; i++) D = idealmul(nf,D,gel(I,i));
    if (pf) *pf = idealinv(nf, D);
    D = idealpow(nf,D,gen_2);
  }
  p1 = core2partial(Q_content(d), 0);
  *pd = gdiv(d, sqri(gel(p1,2)));
  *pD = idealmul(nf,D,d); return z;
}

GEN
rnfpseudobasis(GEN nf, GEN pol)
{
  pari_sp av = avma;
  GEN D, d, y = cgetg(5, t_VEC), z = rnfallbase(nf,pol, &D, &d, NULL);
  y[1] = z[1];
  y[2] = z[2];
  gel(y,3) = D;
  gel(y,4) = d; return gerepilecopy(av, y);
}

GEN
rnfdiscf(GEN nf, GEN pol)
{
  pari_sp av = avma;
  GEN D, d; (void)rnfallbase(nf,pol, &D, &d, NULL);
  return gerepilecopy(av, mkvec2(D,d));
}

GEN
gen_if_principal(GEN bnf, GEN x)
{
  pari_sp av = avma;
  GEN z = isprincipalall(bnf,x, nf_GEN_IF_PRINCIPAL | nf_FORCE);
  if (typ(z) == t_INT) { avma = av; return NULL; }
  return z;
}

/* given bnf and a pseudo-basis of an order in HNF [A,I] (or [A,I,D,d] it
 * does not matter), tries to simplify the HNF as much as possible. The
 * resulting matrix will be upper triangular but the diagonal coefficients
 * will not be equal to 1. The ideals are guaranteed to be integral and
 * primitive. */
GEN
rnfsimplifybasis(GEN bnf, GEN x)
{
  pari_sp av = avma;
  long i, l;
  GEN p1, id, Az, Iz, nf, A, I;

  bnf = checkbnf(bnf); nf = gel(bnf,7);
  if (typ(x)!=t_VEC || lg(x)<3)
    pari_err(talker,"not a pseudo-basis in nfsimplifybasis");
  x = shallowcopy(x);
  A = gel(x,1);
  I = gel(x,2); l = lg(I);
  id = matid(degpol(nf[1]));
  Az = cgetg(l, t_MAT); gel(x,1) = Az;
  Iz = cgetg(l, t_VEC); gel(x,2) = Iz;
  for (i = 1; i < l; i++)
  {
    if (gequal(gel(I,i),id)) { gel(Iz,i) = id; Az[i] = A[i]; continue; }

    gel(Iz,i) = Q_primitive_part(gel(I,i), &p1);
    gel(Az,i) = p1? gmul(gel(A,i),p1): gel(A,i);
    if (p1 && gequal(gel(Iz,i), id)) continue;

    p1 = gen_if_principal(bnf, gel(Iz,i));
    if (p1)
    {
      gel(Iz,i) = id;
      gel(Az,i) = element_mulvec(nf,p1,gel(Az,i));
    }
  }
  return gerepilecopy(av, x);
}

GEN
rnfdet2(GEN nf, GEN A, GEN I)
{
  pari_sp av = avma;
  GEN p1;
  nf = checknf(nf);
  if (typ(I) != t_VEC) pari_err(typeer,"rnfdet2");
  p1 = idealmul(nf, det(matbasistoalg(nf, A)), prodid(nf, I));
  return gerepileupto(av, p1);
}

GEN
rnfdet(GEN nf, GEN order)
{
  if (typ(order)!=t_VEC || lg(order)<3)
    pari_err(talker,"not a pseudo-matrix in rnfdet");
  return rnfdet2(nf,gel(order,1),gel(order,2));
}

GEN
rnfdet0(GEN nf, GEN x, GEN y)
{
  return y? rnfdet2(nf,x,y): rnfdet(nf,x);
}

/* Given two fractional ideals a and b, gives x in a, y in b, z in b^-1,
   t in a^-1 such that xt-yz=1. In the present version, z is in Z. */
static GEN
nfidealdet1(GEN nf, GEN a, GEN b)
{
  pari_sp av = avma;
  GEN x,p1,res,u,v,da,db;

  a = idealinv(nf,a);
  da = denom(a); if (!gcmp1(da)) a = gmul(da,a);
  db = denom(b); if (!gcmp1(db)) b = gmul(db,b);
  x = idealcoprime(nf,a,b);
  p1 = idealaddtoone(nf, idealmul(nf,x,a), b);
  u = gel(p1,1);
  v = gel(p1,2);

  res = cgetg(5,t_VEC);
  gel(res,1) = gmul(x,da);
  gel(res,2) = gdiv(v,db);
  gel(res,3) = negi(db);
  gel(res,4) = element_div(nf, u, gel(res,1));
  return gerepilecopy(av,res);
}

static GEN
get_order(GEN nf, GEN O, char *s)
{
  if (typ(O) == t_POL)
    return rnfpseudobasis(nf, O);
  if (typ(O)!=t_VEC || lg(O) < 3 || typ(O[1]) != t_MAT || typ(O[2]) != t_VEC
      || lg(O[1]) != lg(O[2]))
    pari_err(talker,"not a pseudo-matrix in %s", s);
  return O;
}

/* given a pseudo-basis of an order in HNF [A,I] (or [A,I,D,d]), gives an
 * n x n matrix (not in HNF) of a pseudo-basis and an ideal vector
 * [id,id,...,id,I] such that order = nf[7]^(n-1) x I.
 * Uses the approximation theorem ==> slow. */
GEN
rnfsteinitz(GEN nf, GEN order)
{
  pari_sp av = avma;
  long i,n,l;
  GEN Id,A,I,p1,a,b;

  nf = checknf(nf);
  Id = matid(degpol(nf[1]));
  order = get_order(nf, order, "rnfsteinitz");
  A = matalgtobasis(nf, gel(order,1));
  I = shallowcopy(gel(order,2)); n=lg(A)-1;
  for (i=1; i<n; i++)
  {
    GEN c1,c2;
    a = gel(I,i); if (gequal(a,Id)) continue;

    c1 = gel(A,i);
    c2 = gel(A,i+1);
    b = gel(I,i+1);
    if (gequal(b,Id))
    {
      gel(A,i) = c2;
      gel(A,i+1) = gneg(c1);
      gel(I,i) = b;
      gel(I,i+1) = a;
    }
    else
    {
      p1 = nfidealdet1(nf,a,b);
      gel(A,i) = gadd(element_mulvec(nf, gel(p1,1), c1),
                      element_mulvec(nf, gel(p1,2), c2));
      gel(A,i+1) = gadd(element_mulvec(nf, gel(p1,3), c1),
                        element_mulvec(nf, gel(p1,4), c2));
      gel(I,i) = Id;
      gel(I,i+1) = Q_primitive_part(idealmul(nf,a,b), &p1);
      if (p1) gel(A,i+1) = element_mulvec(nf, p1,gel(A,i+1));
    }
  }
  l = lg(order);
  p1 = cgetg(l,t_VEC);
  gel(p1,1) = A;
  gel(p1,2) = I; for (i=3; i<l; i++) p1[i]=order[i];
  return gerepilecopy(av, p1);
}

/* Given bnf and either an order as output by rnfpseudobasis or a polynomial,
 * and outputs a basis if it is free, an n+1-generating set if it is not */
GEN
rnfbasis(GEN bnf, GEN order)
{
  pari_sp av = avma;
  long j, n;
  GEN nf, A, I, cl, col, a, id;

  bnf = checkbnf(bnf); nf = gel(bnf,7);
  id = matid(degpol(nf[1]));
  order = get_order(nf, order, "rnfbasis");
  I = gel(order,2); n = lg(I)-1;
  j=1; while (j<n && gequal(gel(I,j),id)) j++;
  if (j<n)
  {
    order = rnfsteinitz(nf,order);
    I = gel(order,2);
  }
  A = gel(order,1);
  col= gel(A,n); A = vecslice(A, 1, n-1);
  cl = gel(I,n);
  a = gen_if_principal(bnf, cl);
  if (!a)
  {
    GEN p1 = ideal_two_elt(nf, cl);
    A = shallowconcat(A, gmul(gel(p1,1), col));
    a = gel(p1,2);
  }
  A = shallowconcat(A, element_mulvec(nf, a, col));
  return gerepilecopy(av, A);
}

/* Given bnf and either an order as output by rnfpseudobasis or a polynomial,
 * and outputs a basis (not pseudo) in Hermite Normal Form if it exists, zero
 * if not
 */
GEN
rnfhnfbasis(GEN bnf, GEN order)
{
  pari_sp av = avma;
  long j, n;
  GEN nf, A, I, a, id;

  bnf = checkbnf(bnf); nf = gel(bnf,7);
  id = matid(degpol(nf[1]));
  order = get_order(nf, order, "rnfbasis");
  A = gel(order,1); A = shallowcopy(A);
  I = gel(order,2); n = lg(A)-1;
  for (j=1; j<=n; j++)
  {
    if (gequal(gel(I,j),id)) continue;

    a = gen_if_principal(bnf, gel(I,j));
    if (!a) { avma = av; return gen_0; }
    gel(A,j) = element_mulvec(nf, a, gel(A,j));
  }
  return gerepilecopy(av,A);
}

static long
_rnfisfree(GEN bnf, GEN order)
{
  long n, j;
  GEN nf, p1, id, I;

  bnf = checkbnf(bnf);
  if (gcmp1(gmael3(bnf,8,1,1))) return 1;

  nf = gel(bnf,7); id = matid(degpol(nf[1]));
  order = get_order(nf, order, "rnfisfree");
  I = gel(order,2); n = lg(I)-1;
  j=1; while (j<=n && gequal(gel(I,j),id)) j++;
  if (j>n) return 1;

  p1 = gel(I,j);
  for (j++; j<=n; j++)
    if (!gequal(gel(I,j),id)) p1 = idealmul(nf,p1,gel(I,j));
  return gcmp0( isprincipal(bnf,p1) );
}

long
rnfisfree(GEN bnf, GEN order)
{
  pari_sp av = avma;
  long n = _rnfisfree(bnf, order);
  avma = av; return n;
}

/**********************************************************************/
/**								     **/
/**		      COMPOSITUM OF TWO NUMBER FIELDS                **/
/**								     **/
/**********************************************************************/
/* modular version */
GEN
polcompositum0(GEN A, GEN B, long flall)
{
  pari_sp av = avma;
  int same = (A == B || gequal(A,B));
  long v, k;
  GEN C, D, LPRS;

  if (typ(A)!=t_POL || typ(B)!=t_POL) pari_err(typeer,"polcompositum0");
  if (degpol(A)<=0 || degpol(B)<=0) pari_err(constpoler,"compositum");
  v = varn(A);
  if (varn(B) != v) pari_err(talker,"not the same variable in compositum");
  A = Q_primpart(A); check_ZX(A,"compositum");
  if (!ZX_is_squarefree(A)) pari_err(talker,"compositum: %Z inseparable", A);
  if (!same) {
    B = Q_primpart(B); check_ZX(B,"compositum");
    if (!ZX_is_squarefree(B)) pari_err(talker,"compositum: %Z inseparable", B);
  }

  D = NULL; /* -Wall */
  k = same? -1: 1;
  C = ZY_ZXY_resultant_all(A, B, &k, flall? &LPRS: NULL);
  if (same)
  {
    D = RgX_rescale(A, stoi(1 - k));
    C = gdivexact(C, D);
    if (degpol(C) <= 0) C = mkvec(D); else C = shallowconcat(ZX_DDF(C, 0), D);
  }
  else
    C = ZX_DDF(C, 0); /* C = Res_Y (A, B(X + kY)) guaranteed squarefree */
  C = sort_vecpol(C, &cmpii);
  if (flall)
  {
    long i,l = lg(C);
    GEN w,a,b; /* a,b,c root of A,B,C = compositum, c = b - k a */
    for (i=1; i<l; i++)
    { /* invmod possibly very costly */
      a = gmul(gel(LPRS,1), QXQ_inv(gel(LPRS,2), gel(C,i)));
      a = gneg_i(RgX_rem(a, gel(C,i)));
      b = gadd(pol_x[v], gmulsg(k,a));
      w = cgetg(5,t_VEC); /* [C, a, b, n ] */
      w[1] = C[i];
      gel(w,2) = mkpolmod(a, gel(w,1));
      gel(w,3) = mkpolmod(b, gel(w,1));
      gel(w,4) = stoi(-k); gel(C,i) = w;
    }
  }
  settyp(C, t_VEC); return gerepilecopy(av, C);
}

GEN
compositum(GEN pol1,GEN pol2)
{
  return polcompositum0(pol1,pol2,0);
}

GEN
compositum2(GEN pol1,GEN pol2)
{
  return polcompositum0(pol1,pol2,1);
}

int
nfissquarefree(GEN nf, GEN x)
{
  pari_sp av = avma;
  GEN g, y = derivpol(x);
  if (RgX_is_rational(x))
    g = modulargcd(x, y);
  else
    g = nfgcd(x, y, nf, NULL);
  avma = av; return (degpol(g) == 0);
}

GEN
rnfequation_i(GEN A, GEN B, long *pk, GEN *pLPRS)
{
  long k, lA, lB;
  GEN nf, C;

  A = get_nfpol(A, &nf);       lA = lg(A);
  B = fix_relative_pol(A,B,1); lB = lg(B);
  if (lA<=3 || lB<=3) pari_err(constpoler,"rnfequation");

  check_ZX(A,"rnfequation");
  B = primpart(lift_intern(B));
  check_ZXY(B,"rnfequation");
  for (k=2; k<lB; k++)
    if (lg(B[k]) >= lA) gel(B,k) = grem(gel(B,k),A);

  if (!nfissquarefree(A,B))
    pari_err(talker,"inseparable relative equation in rnfequation");

  *pk = 0; C = ZY_ZXY_resultant_all(A, B, pk, pLPRS);
  if (gsigne(leading_term(C)) < 0) C = gneg_i(C);
  *pk = -*pk; return primpart(C);
}

GEN
rnfequation0(GEN A, GEN B, long flall)
{
  pari_sp av = avma;
  long k;
  GEN LPRS, nf, C;

  A = get_nfpol(A, &nf);
  C = rnfequation_i(A, B, &k, flall? &LPRS: NULL);
  if (flall)
  { /* a,b,c root of A,B,C = compositum, c = b + k a */
    GEN a = gmul(gel(LPRS,1), QXQ_inv(gel(LPRS,2), C));/* inv is costly !*/
    a = gneg_i(RgX_rem(a, C));
    /* b = gadd(pol_x[varn(A)], gmulsg(k,a)); */
    C = mkvec3(C, mkpolmod(a, C), stoi(k));
  }
  return gerepilecopy(av, C);
}

GEN
rnfequation(GEN nf,GEN pol2)
{
  return rnfequation0(nf,pol2,0);
}

GEN
rnfequation2(GEN nf,GEN pol2)
{
  return rnfequation0(nf,pol2,1);
}

static GEN
nftau(long r1, GEN x)
{
  long i, l = lg(x);
  GEN s = r1? gel(x,1): gmul2n(real_i(gel(x,1)),1);
  for (i=2; i<=r1; i++) s = gadd(s, gel(x,i));
  for (   ; i < l; i++) s = gadd(s, gmul2n(real_i(gel(x,i)),1));
  return s;
}

static GEN
initmat(long l)
{
  GEN x = cgetg(l, t_MAT);
  long i;
  for (i = 1; i < l; i++) gel(x,i) = cgetg(l, t_COL);
  return x;
}

static GEN
nftocomplex(GEN nf, GEN x)
{
  return gmul(gmael(nf,5,1), algtobasis_i(nf, x));
}
/* assume x a square t_MAT, return a t_VEC of embeddings of its columns */
static GEN
mattocomplex(GEN nf, GEN x)
{
  long i,j, l = lg(x);
  GEN v = cgetg(l, t_VEC);
  for (j=1; j<l; j++)
  {
    GEN c = gel(x,j), b = cgetg(l, t_MAT);
    for (i=1; i<l; i++) gel(b,i) = nftocomplex(nf, gel(c,i));
    b = shallowtrans(b); settyp(b, t_COL);
    gel(v,j) = b;
  }
  return v;
}

static GEN
nf_all_roots(GEN nf, GEN x, long prec)
{
  long i, j, l = lg(x), ru = lg(nf[6]);
  GEN y = cgetg(l, t_POL), v, z;

  x = unifpol(nf, x, t_COL);
  y[1] = x[1];
  for (i=2; i<l; i++) gel(y,i) = nftocomplex(nf, gel(x,i));
  i = gprecision(y); if (i && i <= 3) return NULL;

  v = cgetg(ru, t_VEC);
  z = cgetg(l, t_POL); z[1] = x[1];
  for (i=1; i<ru; i++)
  {
    for (j = 2; j < l; j++) gel(z,j) = gmael(y,j,i);
    gel(v,i) = cleanroots(z, prec);
  }
  return v;
}

static GEN
rnfscal(GEN m, GEN x, GEN y)
{
  long i, l = lg(m);
  GEN z = cgetg(l, t_COL);
  for (i = 1; i < l; i++)
    gel(z,i) = gmul(gconj(shallowtrans(gel(x,i))), gmul(gel(m,i), gel(y,i)));
  return z;
}

static GEN
findmin(GEN nf, GEN x, GEN muf)
{
  pari_sp av = avma;
  long e;
  GEN cx, y, m, M = gmael(nf,5,1);

  x = Q_primitive_part(x, &cx);
  if (gcmp1(gcoeff(x,1,1))) y = M;
  else
  {
    GEN G = gmael(nf,5,2);
    m = lllintern(gmul(G, x), 4, 1, 0);
    if (!m)
    {
      x = lllint_ip(x,4);
      m = lllintern(gmul(G, x), 4, 1, 0);
      if (!m) pari_err(precer,"rnflllgram");
    }
    x = gmul(x, m);
    y = gmul(M, x);
  }
  m = gauss_realimag(y, muf);
  if (cx) m = gdiv(m, cx);
  m = grndtoi(m, &e);
  if (e >= 0) return NULL; /* precision problem */
  if (cx) m = gmul(m, cx);
  return gerepileupto(av, gmul(x,m));
}

static int
RED(long k, long l, GEN U, GEN mu, GEN MC, GEN nf, GEN I, GEN *Ik_inv)
{
  GEN x, xc, ideal;
  long i;

  if (!*Ik_inv) *Ik_inv = idealinv(nf, gel(I,k));
  ideal = idealmul(nf,gel(I,l), *Ik_inv);
  x = findmin(nf, ideal, gcoeff(mu,k,l));
  if (!x) return 0;
  if (gcmp0(x)) return 1;

  xc = nftocomplex(nf,x);
  gel(MC,k) = gsub(gel(MC,k), vecmul(xc,gel(MC,l)));
  gel(U,k) = gsub(gel(U,k), gmul(coltoalg(nf,x), gel(U,l)));
  gcoeff(mu,k,l) = gsub(gcoeff(mu,k,l), xc);
  for (i=1; i<l; i++)
    gcoeff(mu,k,i) = gsub(gcoeff(mu,k,i), vecmul(xc,gcoeff(mu,l,i)));
  return 1;
}

static int
check_0(GEN B)
{
  long i, l = lg(B);
  for (i = 1; i < l; i++)
    if (gsigne(gel(B,i)) <= 0) return 1;
  return 0;
}

static int
do_SWAP(GEN I, GEN MC, GEN MCS, GEN h, GEN mu, GEN B, long kmax, long k,
        const long alpha, long r1)
{
  GEN p1, p2, muf, mufc, Bf, temp;
  long i, j;

  p1 = nftau(r1, gadd(gel(B,k),
                      gmul(gnorml2(gcoeff(mu,k,k-1)), gel(B,k-1))));
  p2 = nftau(r1, gel(B,k-1));
  if (gcmp(gmulsg(alpha,p1), gmulsg(alpha-1,p2)) > 0) return 0;

  lswap(MC[k-1],MC[k]);
  lswap(h[k-1],  h[k]);
  lswap(I[k-1],  I[k]);
  for (j=1; j<=k-2; j++) swap(gcoeff(mu,k-1,j),gcoeff(mu,k,j));
  muf = gcoeff(mu,k,k-1);
  mufc = gconj(muf);
  Bf = gadd(gel(B,k), vecmul(real_i(vecmul(muf,mufc)), gel(B,k-1)));
  if (check_0(Bf)) return 1; /* precision problem */

  p1 = vecdiv(gel(B,k-1),Bf);
  gcoeff(mu,k,k-1) = vecmul(mufc,p1);
  temp = gel(MCS,k-1);
  gel(MCS,k-1) = gadd(gel(MCS,k), vecmul(muf,gel(MCS,k-1)));
  gel(MCS,k) = gsub(vecmul(vecdiv(gel(B,k),Bf), temp),
                    vecmul(gcoeff(mu,k,k-1), gel(MCS,k)));
  gel(B,k) = vecmul(gel(B,k),p1);
  gel(B,k-1) = Bf;
  for (i=k+1; i<=kmax; i++)
  {
    temp = gcoeff(mu,i,k);
    gcoeff(mu,i,k) = gsub(gcoeff(mu,i,k-1), vecmul(muf, gcoeff(mu,i,k)));
    gcoeff(mu,i,k-1) = gadd(temp, vecmul(gcoeff(mu,k,k-1),gcoeff(mu,i,k)));
  }
  return 1;
}

static GEN
rel_T2(GEN nf, GEN pol, long lx, long prec)
{
  long ru, i, j, k, l;
  GEN T2, s, unro, roorder, powreorder;

  roorder = nf_all_roots(nf, pol, prec);
  if (!roorder) return NULL;
  ru = lg(roorder);
  unro = cgetg(lx,t_COL); for (i=1; i<lx; i++) gel(unro,i) = gen_1;
  powreorder = cgetg(lx,t_MAT); gel(powreorder,1) = unro;
  T2 = cgetg(ru, t_VEC);
  for (i = 1; i < ru; i++)
  {
    GEN ro = gel(roorder,i);
    GEN m = initmat(lx);
    for (k=2; k<lx; k++)
    {
      GEN c = cgetg(lx, t_COL); gel(powreorder,k) = c;
      for (j=1; j < lx; j++)
	gel(c,j) = gmul(gel(ro,j), gmael(powreorder,k-1,j));
    }
    for (l = 1; l < lx; l++)
      for (k = 1; k <= l; k++)
      {
        s = gen_0;
        for (j = 1; j < lx; j++)
          s = gadd(s, gmul(gconj(gmael(powreorder,k,j)),
                                 gmael(powreorder,l,j)));
        if (l == k)
          gcoeff(m, l, l) = real_i(s);
        else
        {
          gcoeff(m, k, l) = s;
          gcoeff(m, l, k) = gconj(s);
        }
      }
    gel(T2,i) = m;
  }
  return T2;
}

/* given a base field nf (e.g main variable y), a polynomial pol with
 * coefficients in nf    (e.g main variable x), and an order as output
 * by rnfpseudobasis, outputs a reduced order. */
GEN
rnflllgram(GEN nf, GEN pol, GEN order,long prec)
{
  pari_sp av = avma, lim = stack_lim(av,2);
  long j, k, l, kmax, r1, lx, count = 0;
  GEN M, I, h, H, mth, MC, MPOL, MCS, B, mu, y;
  const long alpha = 10, MAX_COUNT = 4;

  nf = checknf(nf); r1 = nf_get_r1(nf);
  check_ZKmodule(order, "rnflllgram");
  M = gel(order,1);
  I = gel(order,2); lx = lg(I);
  if (lx < 3) return gcopy(order);
  if (lx-1 != degpol(pol)) pari_err(consister,"rnflllgram");
  I = shallowcopy(I);
  H = NULL;
  MPOL = matbasistoalg(nf, M);
  MCS = matid(lx-1); /* dummy for gerepile */
PRECNF:
  if (count == MAX_COUNT)
  {
    prec = (prec<<1)-2; count = 0;
    if (DEBUGLEVEL) pari_warn(warnprec,"rnflllgram",prec);
    nf = nfnewprec(nf,prec);
  }
  mth = rel_T2(nf, pol, lx, prec);
  if (!mth) { count = MAX_COUNT; goto PRECNF; }
  h = NULL;
PRECPB:
  if (h)
  { /* precision problem, recompute. If no progress, increase nf precision */
    if (++count == MAX_COUNT || isidentity(h)) {count = MAX_COUNT; goto PRECNF;}
    H = H? gmul(H, h): h;
    MPOL = gmul(MPOL, h);
  }
  h = matid(lx-1);
  MC = mattocomplex(nf, MPOL);
  mu = cgetg(lx,t_MAT);
  B  = cgetg(lx,t_COL);
  for (j=1; j<lx; j++)
  {
    gel(mu,j) = zerocol(lx - 1);
    gel(B,j) = gen_0;
  }
  if (DEBUGLEVEL) fprintferr("k = ");
  gel(B,1) = real_i(rnfscal(mth,gel(MC,1),gel(MC,1)));
  MCS[1] = MC[1];
  kmax = 1; k = 2;
  do
  {
    GEN Ik_inv = NULL;
    if (DEBUGLEVEL) fprintferr("%ld ",k);
    if (k > kmax)
    { /* Incremental Gram-Schmidt */
      kmax = k; MCS[k] = MC[k];
      for (j=1; j<k; j++)
      {
	gcoeff(mu,k,j) = vecdiv(rnfscal(mth,gel(MCS,j),gel(MC,k)),
	                        gel(B,j));
	gel(MCS,k) = gsub(gel(MCS,k), vecmul(gcoeff(mu,k,j),gel(MCS,j)));
      }
      gel(B,k) = real_i(rnfscal(mth,gel(MCS,k),gel(MCS,k)));
      if (check_0(gel(B,k))) goto PRECPB;
    }
    if (!RED(k, k-1, h, mu, MC, nf, I, &Ik_inv)) goto PRECPB;
    if (do_SWAP(I,MC,MCS,h,mu,B,kmax,k,alpha, r1))
    {
      if (!B[k]) goto PRECPB;
      if (k > 2) k--;
    }
    else
    {
      for (l=k-2; l; l--)
        if (!RED(k, l, h, mu, MC, nf, I, &Ik_inv)) goto PRECPB;
      k++;
    }
    if (low_stack(lim, stack_lim(av,2)))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"rnflllgram");
      gerepileall(av, H?10:9, &nf,&mth,&h,&MPOL,&B,&MC,&MCS,&mu,&I,&H);
    }
  }
  while (k < lx);
  MPOL = gmul(MPOL,h);
  if (H) h = gmul(H, h);
  if (DEBUGLEVEL) fprintferr("\n");
  y = cgetg(3,t_VEC);
  gel(y,1) = mkvec2(algtobasis(nf,MPOL), gcopy(I));
  gel(y,2) = algtobasis(nf,h); return gerepileupto(av, y);
}

GEN
rnfpolred(GEN nf, GEN pol, long prec)
{
  pari_sp av = avma;
  long i, j, n, v = varn(pol);
  GEN id, al, w, I, O, bnf, nfpol;

  if (typ(pol)!=t_POL) pari_err(typeer,"rnfpolred");
  bnf = nf; nf = checknf(bnf);
  bnf = (nf == bnf)? NULL: checkbnf(bnf);
  if (degpol(pol) <= 1) return mkvec(pol_x[v]);
  nfpol = gel(nf,1);

  id = rnfpseudobasis(nf,pol);
  if (bnf && gcmp1(gmael3(bnf,8,1,1))) /* if bnf is principal */
  {
    GEN newI, newO, zk = matid(degpol(nfpol));
    O = gel(id,1);
    I = gel(id,2); n = lg(I)-1;
    newI = cgetg(n+1,t_VEC);
    newO = cgetg(n+1,t_MAT);
    for (j=1; j<=n; j++)
    {
      gel(newI,j) = zk; al = gen_if_principal(bnf,gel(I,j));
      gel(newO,j) = element_mulvec(nf, al, gel(O,j));
    }
    id = mkvec2(newO, newI);
  }

  id = (GEN)rnflllgram(nf,pol,id,prec)[1];
  O = gel(id,1);
  I = gel(id,2); n = lg(I)-1;
  w = cgetg(n+1,t_VEC);
  pol = lift(pol);
  for (j=1; j<=n; j++)
  {
    GEN p1, newpol;

    p1 = gel(I,j); al = gmul(gcoeff(p1,1,1),gel(O,j));
    p1 = coltoalg(nf,gel(al,n));
    for (i=n-1; i; i--)
      p1 = gadd(coltoalg(nf,gel(al,i)), gmul(pol_x[v],p1));
    newpol = RgXQX_red(caract2(pol,lift(p1),v), nfpol);
    newpol = Q_primpart(newpol);

    p1 = nfgcd(newpol, derivpol(newpol), nfpol, gel(nf,4));
    if (degpol(p1) > 0) newpol = RgXQX_div(newpol, p1, nfpol);
    p1 = leading_term(newpol);
    if (typ(p1) != t_POL) p1 = scalarpol(p1, varn(nfpol));
    newpol = RgXQX_div(newpol, p1, nfpol);
    gel(w,j) = newpol;
  }
  return gerepilecopy(av,w);
}

/* given a relative polynomial pol over nf, compute a pseudo-basis for the
 * extension, then an absolute basis */
static GEN
makebasis(GEN nf, GEN pol, GEN rnfeq)
{
  GEN elts,ids,polabs,plg,plg0,B,bs,p1,den,vbs,d,vpro;
  pari_sp av = avma;
  long n,N,m,i,j,k, v = varn(pol);

  polabs= gel(rnfeq,1);
  plg   = gel(rnfeq,2); plg = lift_intern(plg);
  p1 = rnfpseudobasis(nf,pol);
  elts= gel(p1,1);
  ids = gel(p1,2);
  if (DEBUGLEVEL>1) fprintferr("relative basis computed\n");
  N = degpol(pol);
  n = degpol(nf[1]); m = n*N;

  plg0= Q_remove_denom(plg, &den); /* plg = plg0/den */
  /* nf = K = Q(a), vbs[i+1] = a^i as an elt of L = Q[X] / polabs */
  vbs = RgX_powers(plg0, polabs, n-1);
  if (den)
  { /* restore denominators */
    gel(vbs,2) = plg; d = den;
    for (i=3; i<=n; i++) { d = mulii(d,den); gel(vbs,i) = gdiv(gel(vbs,i), d); }
  }

  /* bs = integer basis of K, as elements of L */
  bs = gmul(vbs, RgXV_to_RgM(gel(nf,7),n));

  vpro = cgetg(N+1,t_VEC);
  for (i=1; i<=N; i++) gel(vpro,i) = monomial(gen_1, i-1, v);
  vpro = gmul(vpro, elts);
  B = cgetg(m+1, t_MAT);
  for(i=k=1; i<=N; i++)
  {
    GEN w = element_mulvec(nf, gel(vpro,i), gel(ids,i));
    for(j=1; j<=n; j++)
    {
      p1 = grem(gmul(bs, gel(w,j)), polabs);
      gel(B,k++) = RgX_to_RgV(p1, m);
    }
  }
  B = Q_remove_denom(B, &den);
  if (den) { B = hnfmodid(B, den); B = gdiv(B, den); } else B = matid(m);
  return gerepilecopy(av, mkvec2(polabs, B));
}

/* relative polredabs. Returns relative polynomial by default (flag = 0)
 * flag & nf_ORIG: + element (base change)
 * flag & nf_ADDZK: + integer basis
 * flag & nf_ABSOLUTE: absolute polynomial */
GEN
rnfpolredabs(GEN nf, GEN relpol, long flag)
{
  GEN red, bas, elt, POL, pol, T, a;
  long v, fl = (flag & nf_ADDZK)? nf_ADDZK: nf_RAW;
  pari_sp av = avma;

  if (typ(relpol)!=t_POL) pari_err(typeer,"rnfpolredabs");
  nf = checknf(nf); v = varn(relpol);
  if (DEBUGLEVEL>1) (void)timer2();
  relpol = unifpol(nf, relpol, t_POLMOD);
  T = gel(nf,1);
  if ((flag & nf_ADDZK) && !(flag & nf_ABSOLUTE))
    pari_err(impl,"this combination of flags in rnfpolredabs");
  if (flag & nf_PARTIALFACT)
  {
    long sa;
    fl |= nf_PARTIALFACT;
    POL = rnfequation_i(nf, relpol, &sa, NULL);
    bas = POL;
    a = stoi(sa);
  }
  else
  {
    GEN rel, eq = rnfequation2(nf,relpol);
    POL = gel(eq,1);
    a   = gel(eq,3);
    rel = poleval(relpol, gsub(pol_x[v],
                               gmul(a, gmodulo(pol_x[varn(T)],T))));
    bas = makebasis(nf, rel, eq);
    if (DEBUGLEVEL>1)
    {
      msgtimer("absolute basis");
      fprintferr("original absolute generator: %Z\n", POL);
    }
  }
  red = polredabs0(bas, fl);
  pol = gel(red,1);
  if (DEBUGLEVEL>1) fprintferr("reduced absolute generator: %Z\n",pol);
  if (flag & nf_ABSOLUTE)
  {
    if (flag & nf_ADDZK) pol = mkvec2(pol, gel(red,2));
    return gerepilecopy(av, pol);
  }

  elt = eltabstorel(gel(red,2), T, relpol, a);
  pol = rnfcharpoly(nf,relpol,elt,v);
  if (!(flag & nf_ORIG)) return gerepileupto(av, pol);
  elt = mkpolmod(modreverse_i(gel(elt,2),gel(elt,1)), pol);
  return gerepilecopy(av, mkvec2(pol, elt));
}
