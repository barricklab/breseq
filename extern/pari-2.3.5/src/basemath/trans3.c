/* $Id: trans3.c 10287 2008-06-09 22:03:21Z kb $

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
/**                   TRANSCENDENTAL FONCTIONS                     **/
/**                          (part 3)                              **/
/**                                                                **/
/********************************************************************/
#include "pari.h"
#include "paripriv.h"

/***********************************************************************/
/**                                                                   **/
/**                       BESSEL FUNCTIONS                            **/
/**                                                                   **/
/***********************************************************************/

/* n! sum_{k=0}^m Z^k / (k!*(k+n)!), with Z := (-1)^flag*z^2/4 */
static GEN
_jbessel(GEN n, GEN z, long flag, long m)
{
  long k;
  pari_sp av, lim;
  GEN Z,s;

  Z = gmul2n(gsqr(z),-2); if (flag & 1) Z = gneg(Z);
  if (typ(z) == t_SER)
  {
    long v = valp(z);
    k = lg(Z)-2 - v;
    if (v < 0) pari_err(negexper,"jbessel");
    if (v == 0) pari_err(impl,"jbessel around a!=0");
    if (k <= 0) return gadd(gen_1, zeroser(varn(z), 2*v));
    Z = gprec(Z, k);
  }
  s = gen_1;
  av = avma; lim = stack_lim(av,1);
  for (k=m; k>=1; k--)
  {
    s = gaddsg(1, gdiv(gdivgs(gmul(Z,s),k),gaddsg(k,n)));
    if (low_stack(lim, stack_lim(av,1)))
    {
      if (DEBUGMEM>1) pari_warn(warnmem,"jbessel");
      s = gerepilecopy(av, s);
    }
  }
  return s;
}

/* return L * approximate solution to x log x = B */
static long
bessel_get_lim(double B, double L)
{
  long lim;
  double x = 1 + B; /* 3 iterations are enough except in pathological cases */
  x = (x + B)/(log(x)+1);
  x = (x + B)/(log(x)+1);
  x = (x + B)/(log(x)+1); x = L*x;
  lim = (long)x; if (lim < 2) lim = 2;
  return lim;
}

static GEN
jbesselintern(GEN n, GEN z, long flag, long prec)
{
  long i, lz, lim, k, ki, precnew;
  pari_sp av = avma;
  double B, L;
  GEN p1, p2, y;

  switch(typ(z))
  {
    case t_INT: case t_FRAC: case t_QUAD:
    case t_REAL: case t_COMPLEX:
      i = precision(z); if (i) prec = i;
      p2 = gdiv(gpow(gmul2n(z,-1),n,prec), ggamma(gaddgs(n,1),prec));
      if (gcmp0(z)) return gerepilecopy(av, p2);

      L = 1.3591409 * gtodouble(gabs(z,prec));
      precnew = prec;
      if (L >= 1.0) precnew += 1 + (long)(L/(1.3591409*LOG2*BITS_IN_LONG));
      if (issmall(n,&ki))
      {
	k = labs(ki);
        n = stoi(k);
      } else {
        i = precision(n);
        if (i && i < precnew) n = gtofp(n,precnew);
      }
      z = gtofp(z,precnew);
      B = bit_accuracy_mul(prec, LOG2/2) / L;
      lim = bessel_get_lim(B, L);
      p1 = gprec_wtrunc(_jbessel(n,z,flag,lim), prec);
      return gerepileupto(av, gmul(p2,p1));

    case t_VEC: case t_COL: case t_MAT:
      lz=lg(z); y=cgetg(lz,typ(z));
      for (i=1; i<lz; i++)
	gel(y,i) = jbesselintern(n,gel(z,i),flag,prec);
      return y;

    case t_POLMOD:
      y = cleanroots(gel(z,1), prec); lz = lg(y);
      for (i=1; i<lz; i++) {
        GEN t = poleval(gel(z,2), gel(y,i));
        gel(y,i) = jbesselintern(n,t,flag,prec);
      }
      return gerepilecopy(av,y);

    case t_PADIC: pari_err(impl,"p-adic jbessel function");
    default:
      if (!(y = toser_i(z))) break;
      if (issmall(n,&ki)) n = stoi(labs(ki));
      return gerepilecopy(av, _jbessel(n,y,flag,lg(y)-2));
  }
  pari_err(typeer,"jbessel");
  return NULL; /* not reached */
}
GEN
jbessel(GEN n, GEN z, long prec) { return jbesselintern(n,z,1,prec); }
GEN
ibessel(GEN n, GEN z, long prec) { return jbesselintern(n,z,0,prec); }

static GEN
_jbesselh(long k, GEN z, long prec)
{
  GEN s,c,p0,p1,p2, zinv = ginv(z);
  long i;

  gsincos(z,&s,&c,prec);
  p1 = gmul(zinv,s);
  if (k)
  {
    p0 = p1; p1 = gmul(zinv,gsub(p0,c));
    for (i=2; i<=k; i++)
    {
      p2 = gsub(gmul(gmulsg(2*i-1,zinv),p1), p0);
      p0 = p1; p1 = p2;
    }
  }
  return p1;
}

GEN
jbesselh(GEN n, GEN z, long prec)
{
  long gz, k, l, linit, i, lz;
  pari_sp av;
  GEN res, y, p1;

  if (typ(n)!=t_INT) pari_err(talker,"not an integer index in jbesselh");
  k = itos(n);
  if (k < 0) return jbessel(gadd(ghalf,n), z, prec);

  switch(typ(z))
  {
    case t_INT: case t_FRAC: case t_QUAD:
    case t_REAL: case t_COMPLEX:
      if (gcmp0(z))
      {
        av = avma;
	p1 = gmul(gsqrt(gdiv(z,mppi(prec)),prec),gpowgs(z,k));
	p1 = gdiv(p1, seq_umul(k+1, 2*k+1)); /* x k! / (2k+1)! */
	return gerepileupto(av, gmul2n(p1,2*k));
      }
      gz = gexpo(z);
      linit = precision(z); if (!linit) linit = prec;
      res = cgetc(linit);
      av = avma;
      if (gz>=0) l = linit;
      else l = linit - 1 + ((-2*k*gz)>>TWOPOTBITS_IN_LONG);
      if (l>prec) prec = l;
      prec += (-gz)>>TWOPOTBITS_IN_LONG;
      if (prec < 3) prec = 3;
      z = gadd(z, real_0(prec));
      if (typ(z) == t_COMPLEX) gel(z,2) = gadd(gel(z,2), real_0(prec));
      p1 = gmul(_jbesselh(k,z,prec), gsqrt(gdiv(z,Pi2n(-1,prec)),prec));
      avma = av;
      if (typ(p1) == t_COMPLEX)
      {
        affr_fixlg(gel(p1,1), gel(res,1));
        affr_fixlg(gel(p1,2), gel(res,2));
      }
      else
      {
        res = cgetr(linit);
        affr_fixlg(p1, res);
      }
      return res;

    case t_VEC: case t_COL: case t_MAT:
      lz=lg(z); y=cgetg(lz,typ(z));
      for (i=1; i<lz; i++) gel(y,i) = jbesselh(n,gel(z,i),prec);
      return y;

    case t_POLMOD:
      av = avma;
      y = cleanroots(gel(z,1), prec); lz = lg(y);
      for (i=1; i<lz; i++) {
        GEN t = poleval(gel(z,2), gel(y,i));
        gel(y,i) = jbesselh(n,t,prec);
      }
      return gerepilecopy(av, y);

    case t_PADIC: pari_err(impl,"p-adic jbesselh function");
    default:
      av = avma;
      if (!(y = toser_i(z))) break;
      if (gcmp0(y)) return gpowgs(y,k);
      if (valp(y) < 0) pari_err(negexper,"jbesselh");
      y = gprec(y, lg(y)-2 + (2*k+1)*valp(y));
      p1 = gdiv(_jbesselh(k,y,prec),gpowgs(y,k));
      for (i=1; i<=k; i++) p1 = gmulgs(p1,2*i+1);
      return gerepilecopy(av,p1);
  }
  pari_err(typeer,"jbesselh");
  return NULL; /* not reached */
}

GEN
kbessel2(GEN nu, GEN x, long prec)
{
  pari_sp av = avma;
  GEN p1, x2, a;

  if (typ(x)==t_REAL) prec = lg(x);
  x2 = gshift(x,1);
  a = gcmp0(imag_i(nu))? cgetr(prec): cgetc(prec);
  gaddz(gen_1,gshift(nu,1), a);
  p1 = hyperu(gshift(a,-1),a,x2,prec);
  p1 = gmul(gmul(p1,gpow(x2,nu,prec)), sqrtr(mppi(prec)));
  return gerepileupto(av, gdiv(p1,gexp(x,prec)));
}

GEN
kbessel(GEN nu, GEN gx, long prec)
{
  GEN x,y,yfin,p1,p2,zf,zz,s,t,q,r,u,v,e,f,c,d,ak,pitemp,nu2,w;
  long l, lnew, k, k2, l1, n2, n, ex, rab=0;
  pari_sp av, av1;

  if (typ(nu)==t_COMPLEX) return kbessel2(nu,gx,prec);
  l = (typ(gx)==t_REAL)? lg(gx): prec;
  ex = gexpo(gx);
  if (ex < 0)
  {
    rab = 1 + ((-ex)>>TWOPOTBITS_IN_LONG);
    lnew = l + rab; prec += rab;
  }
  else lnew = l;
  yfin=cgetr(l); l1=lnew+1;
  av=avma; x = gtofp(gx, lnew); y=cgetr(lnew);
  u=cgetr(l1); v=cgetr(l1); c=cgetr(l1); d=cgetr(l1);
  e=cgetr(l1); f=cgetr(l1);
  nu2=gmulgs(gsqr(nu),-4);
  n = (long) (bit_accuracy_mul(l,LOG2) + PI*sqrt(gtodouble(gnorm(nu)))) / 2;
  n2=(n<<1); pitemp=mppi(l1);
  av1=avma;
  if (cmprs(x, n) < 0)
  {
    zf = gsqrt(gdivgs(pitemp,n2),prec);
    zz = ginv(stor(n2<<2, prec));
    s=gen_1; t=gen_0;
    for (k=n2,k2=2*n2-1; k > 0; k--,k2-=2)
    {
      if (k2 & ~(MAXHALFULONG>>1))
        p1 = gadd(mulss(k2,k2),nu2);
      else
        p1 = gaddsg(k2*k2,nu2);
      ak = gdivgs(gmul(p1,zz),-k);
      s = gaddsg(1, gmul(ak,s));
      t = gaddsg(k2,gmul(ak,t));
    }
    gmulz(s,zf,u); t=gmul2n(t,-1);
    gdivgsz(gadd(gmul(t,zf),gmul(u,nu)),-n2,v); avma=av1;
    q = stor(n2, l1);
    r=gmul2n(x,1); av1=avma;
    for(;;)
    {
      p1=divsr(5,q); if (expo(p1) >= -1) p1=ghalf;
      p2=subsr(1,divrr(r,q));
      if (gcmp(p1,p2)>0) p1=p2;
      gnegz(p1,c); gaffsg(1,d); affrr(u,e); affrr(v,f);
      for (k=1; ; k++)
      {
	w=gadd(gmul(gsubsg(k,ghalf),u), gmul(subrs(q,k),v));
	w=gadd(w, gmul(nu,gsub(u,gmul2n(v,1))));
	gdivgsz(gmul(q,v),k,u);
	gdivgsz(w,k,v);
	gmulz(d,c,d);
	gaddz(e,gmul(d,u),e); p1=gmul(d,v);
	gaddz(f,p1,f);
	if (gcmp0(p1) || gexpo(p1) - gexpo(f) <= 1-bit_accuracy(precision(p1)))
	  break;
	avma=av1;
      }
      p1=u; u=e; e=p1;
      p1=v; v=f; f=p1;
      gmulz(q,gaddsg(1,c),q);
      if (expo(subrr(q,r)) - expo(r) <= 1-bit_accuracy(lnew)) break;
      avma=av1;
    }
    gmulz(u,gpow(divrs(x,n),nu,prec),y);
  }
  else
  {
    p2=gmul2n(x,1);
    zf=gsqrt(gdiv(pitemp,p2),prec);
    zz=ginv(gmul2n(p2,2)); s=gen_1;
    for (k=n2,k2=2*n2-1; k > 0; k--,k2-=2)
    {
      if (k2 & ~(MAXHALFULONG>>1))
        p1=gadd(mulss(k2,k2),nu2);
      else
        p1=gaddsg(k2*k2,nu2);
      ak = gdivgs(gmul(p1,zz),k);
      s = gsub(gen_1,gmul(ak,s));
    }
    gmulz(s,zf,y);
  }
  gdivz(y,mpexp(x),yfin);
  avma=av; return yfin;
}

/*   sum_{k=0}^m Z^k (H(k)+H(k+n)) / (k! (k+n)!)
 * + sum_{k=0}^{n-1} (-Z)^(k-n) (n-k-1)!/k!   with Z := (-1)^flag*z^2/4.
 * Warning: contrary to _jbessel, no n! in front.
 * When flag > 1, compute exactly the H(k) and factorials (slow) */
static GEN
_kbessel(long n, GEN z, long flag, long m, long prec)
{
  long k, limit;
  pari_sp av;
  GEN Z, p1, p2, s, H;

  Z = gmul2n(gsqr(z),-2); if (flag & 1) Z = gneg(Z);
  if (typ(z) == t_SER)
  {
    long v = valp(z);
    k = lg(Z)-2 - v;
    if (v < 0) pari_err(negexper,"kbessel");
    if (v == 0) pari_err(impl,"kbessel around a!=0");
    if (k <= 0) return gadd(gen_1, zeroser(varn(z), 2*v));
    Z = gprec(Z, k);
  }
  H = cgetg(m+n+2,t_VEC); gel(H,1) = gen_0;
  if (flag <= 1)
  {
    gel(H,2) = s = real_1(prec);
    for (k=2; k<=m+n; k++) gel(H,k+1) = s = divrs(addsr(1,mulsr(k,s)),k);
  }
  else
  {
    gel(H,2) = s = gen_1;
    for (k=2; k<=m+n; k++) gel(H,k+1) = s = gdivgs(gaddsg(1,gmulsg(k,s)),k);
  }
  s = gadd(gel(H,m+1), gel(H,m+n+1));
  av = avma; limit = stack_lim(av,1);
  for (k=m; k>0; k--)
  {
    s = gadd(gadd(gel(H,k),gel(H,k+n)),gdiv(gmul(Z,s),mulss(k,k+n)));
    if (low_stack(limit,stack_lim(av,1)))
    {
      if (DEBUGMEM>1) pari_warn(warnmem,"kbessel");
      s = gerepilecopy(av, s);
    }
  }
  p1 = (flag <= 1) ? mpfactr(n,prec) : mpfact(n);
  s = gdiv(s,p1);
  if (n)
  {
    Z = gneg(ginv(Z));
    p2 = gmulsg(n, gdiv(Z,p1));
    s = gadd(s,p2);
    for (k=n-1; k>0; k--)
    {
      p2 = gmul(p2, gmul(mulss(k,n-k),Z));
      s = gadd(s,p2);
    }
  }
  return s;
}

static GEN
kbesselintern(GEN n, GEN z, long flag, long prec)
{
  long i, k, ki, lz, lim, precnew, fl, fl2, ex;
  pari_sp av = avma;
  GEN p1, p2, y, p3, pp, pm, s, c;
  double B, L;

  fl = (flag & 1) == 0;
  switch(typ(z))
  {
    case t_INT: case t_FRAC: case t_QUAD:
    case t_REAL: case t_COMPLEX:
      if (gcmp0(z)) pari_err(talker,"zero argument in a k/n bessel function");
      i = precision(z); if (i) prec = i;
      i = precision(n); if (i && prec > i) prec = i;
      ex = gexpo(z);
      /* experimental */
      if (!flag && ex > bit_accuracy(prec)/16 + gexpo(n))
        return kbessel(n,z,prec);
      L = 1.3591409 * gtodouble(gabs(z,prec));
      precnew = prec;
      if (L >= 1.3591409) {
        long rab = (long)(L/(1.3591409*LOG2*BITS_IN_LONG));
        if (fl) rab *= 2;
         precnew += 1 + rab;
      }
      z = gtofp(z, precnew);
      if (issmall(n,&ki))
      {
        GEN z2 = gmul2n(z, -1);
	k = labs(ki);
	B = bit_accuracy_mul(prec,LOG2/2) / L;
	if (fl) B += 0.367879;
        lim = bessel_get_lim(B, L);
	p1 = gmul(gpowgs(z2,k), _kbessel(k,z,flag,lim,precnew));
	p2 = gadd(mpeuler(precnew), glog(z2,precnew));
	p3 = jbesselintern(stoi(k),z,flag,precnew);
	p2 = gsub(gmul2n(p1,-1),gmul(p2,p3));
	p2 = gprec_wtrunc(p2, prec);
        if (fl) {
          if (k & 1) p2 = gneg(p2);
        }
        else
	{
          p2 = gdiv(p2, Pi2n(-1,prec));
          if (ki >= 0 || (k&1)==0) p2 = gneg(p2);
        }
	return gerepilecopy(av, p2);
      }

      n = gtofp(n, precnew);
      gsincos(gmul(n,mppi(precnew)), &s,&c,precnew);
      ex = gexpo(s);
      if (ex < 0)
      {
        long rab = (-ex) >> TWOPOTBITS_IN_LONG;
        if (fl) rab *= 2;
        precnew += 1 + rab;
      }
      if (i && i < precnew) {
        n = gtofp(n,precnew);
        z = gtofp(z,precnew);
        gsincos(gmul(n,mppi(precnew)), &s,&c,precnew);
      }

      pp = jbesselintern(n,      z,flag,precnew);
      pm = jbesselintern(gneg(n),z,flag,precnew);
      if (fl)
        p1 = gmul(gsub(pm,pp), Pi2n(-1,precnew));
      else
        p1 = gsub(gmul(c,pp),pm);
      p1 = gdiv(p1, s);
      return gerepilecopy(av, gprec_wtrunc(p1,prec));

    case t_VEC: case t_COL: case t_MAT:
      lz=lg(z); y=cgetg(lz,typ(z));
      for (i=1; i<lz; i++) gel(y,i) = kbesselintern(n,gel(z,i),flag,prec);
      return y;

    case t_POLMOD:
      y = cleanroots(gel(z,1), prec); lz = lg(y);
      for (i=1; i<lz; i++) {
        GEN t = poleval(gel(z,2), gel(y,i));
        gel(y,i) = kbesselintern(n,t,flag,prec);
      }
      return gerepilecopy(av, y);

    case t_PADIC: pari_err(impl,"p-adic kbessel function");
    default:
      if (!(y = toser_i(z))) break;
      if (issmall(n,&ki))
      {
	k = labs(ki);
	return gerepilecopy(av, _kbessel(k,y,flag+2,lg(y)-2,prec));
      }
      if (!issmall(gmul2n(n,1),&ki))
        pari_err(talker,"cannot give a power series result in k/n bessel function");
      k = labs(ki); n = gmul2n(stoi(k),-1);
      fl2 = (k&3)==1;
      pm = jbesselintern(gneg(n),y,flag,prec);
      if (fl)
      {
        pp = jbesselintern(n,y,flag,prec);
        p2 = gpowgs(y,-k); if (fl2 == 0) p2 = gneg(p2);
        p3 = gmul2n(diviiexact(mpfact(k + 1),mpfact((k + 1) >> 1)),-(k + 1));
        p3 = gdivgs(gmul2n(gsqr(p3),1),k);
        p2 = gmul(p2,p3);
        p1 = gsub(pp,gmul(p2,pm));
      }
      else p1 = pm;
      return gerepileupto(av, fl2? gneg(p1): gcopy(p1));
  }
  pari_err(typeer,"kbessel");
  return NULL; /* not reached */
}

GEN
kbesselnew(GEN n, GEN z, long prec) { return kbesselintern(n,z,0,prec); }
GEN
nbessel(GEN n, GEN z, long prec) { return kbesselintern(n,z,1,prec); }
/* J + iN */
GEN
hbessel1(GEN n, GEN z, long prec)
{
  pari_sp av = avma;
  return gerepileupto(av, gadd(jbessel(n,z,prec), mulcxI(nbessel(n,z,prec))));
}
/* J - iN */
GEN
hbessel2(GEN n, GEN z, long prec)
{
  pari_sp av = avma;
  return gerepileupto(av, gadd(jbessel(n,z,prec), mulcxmI(nbessel(n,z,prec))));
}

GEN
kbessel0(GEN nu, GEN gx, long flag, long prec)
{
  switch(flag)
  {
    case 0: return kbesselnew(nu,gx,prec);
    case 1: return kbessel(nu,gx,prec);
    case 2: return kbessel2(nu,gx,prec);
    default: pari_err(flagerr,"besselk");
  }
  return NULL; /* not reached */
}

/***********************************************************************/
/*                                                                    **/
/**                    FONCTION U(a,b,z) GENERALE                     **/
/**                         ET CAS PARTICULIERS                       **/
/**                                                                   **/
/***********************************************************************/
/* Assume gx > 0 and a,b complex */
/* This might one day be extended to handle complex gx */
/* see Temme, N. M. "The numerical computation of the confluent        */
/* hypergeometric function U(a,b,z)" in Numer. Math. 41 (1983),        */
/* no. 1, 63--82.                                                      */
GEN
hyperu(GEN a, GEN b, GEN gx, long prec)
{
  GEN x,y,p1,p2,p3,zf,zz,s,t,q,r,u,v,e,f,c,d,w,a1,gn;
  long l, k, l1, n, ex;
  pari_sp av, av1, av2;

  if(gsigne(gx) <= 0) pari_err(talker,"hyperu's third argument must be positive");
  ex = (iscomplex(a) || iscomplex(b));

  l = (typ(gx)==t_REAL)? lg(gx): prec;
  if (ex) y=cgetc(l); else y=cgetr(l);
  l1=l+1; av=avma; x = gtofp(gx, l);
  a1=gaddsg(1,gsub(a,b));
  if (ex)
  {
    u=cgetc(l1); v=cgetc(l1); c=cgetc(l1);
    d=cgetc(l1); e=cgetc(l1); f=cgetc(l1);
  }
  else
  {
    u=cgetr(l1); v=cgetr(l1); c=cgetr(l1);
    d=cgetr(l1); e=cgetr(l1); f=cgetr(l1);
  }
  n=(long)(bit_accuracy_mul(l, LOG2) + PI*sqrt(gtodouble(gabs(gmul(a,a1),l1))));
  av1=avma;
  if (cmprs(x,n)<0)
  {
    gn=stoi(n); zf=gpow(gn,gneg_i(a),l1);
    zz=gdivsg(-1,gn); s=gen_1; t=gen_0;
    for (k=n-1; k>=0; k--)
    {
      p1=gdivgs(gmul(gmul(gaddgs(a,k),gaddgs(a1,k)),zz),k+1);
      s=gaddsg(1,gmul(p1,s));
      t=gadd(gaddsg(k,a),gmul(p1,t));
    }
    gmulz(s,zf,u); t=gmul(t,zz); gmulz(t,zf,v); avma=av1;
    q = stor(n, l1); r=x; av1=avma;
    for(;;)
    {
      p1=divsr(5,q); if (expo(p1)>= -1) p1=ghalf;
      p2=subsr(1,divrr(r,q)); if (gcmp(p1,p2)>0) p1=p2;
      gnegz(p1,c); avma=av1;
      gaffsg(1,d); gaffect(u,e); gaffect(v,f);
      p3=gsub(q,b); av2 = avma;
      for(k=1;;k++)
      {
	w=gadd(gmul(gaddsg(k-1,a),u),gmul(gaddsg(1-k,p3),v));
	gdivgsz(gmul(q,v),k,u);
	gdivgsz(w,k,v);
	gmulz(d,c,d);
	gaddz(e,gmul(d,u),e); p1=gmul(d,v);
	gaddz(f,p1,f);
	if (gcmp0(p1) || gexpo(p1) - gexpo(f) <= 1-bit_accuracy(precision(p1))) break;
	avma=av2;
      }
      p1=u; u=e; e=p1;
      p1=v; v=f; f=p1;
      gmulz(q,gadd(gen_1,c),q);
      if (expo(subrr(q,r)) - expo(r) <= 1-bit_accuracy(l)) break;
      avma=av1;
    }
  }
  else
  {
    zf=gpow(x,gneg_i(a),l1);
    zz=divsr(-1,x); s=gen_1;
    for (k=n-1; k>=0; k--)
    {
      p1=gdivgs(gmul(gmul(gaddgs(a,k),gaddgs(a1,k)),zz),k+1);
      s=gadd(gen_1,gmul(p1,s));
    }
    u = gmul(s,zf);
  }
  gaffect(u,y); avma=av; return y;
}

/* = incgam2(0, x, prec). typ(x) = t_REAL. Optimized for eint1 */
static GEN
incgam2_0(GEN x, GEN expx)
{
  long l = lg(x), n, i;
  GEN z;

  if (expo(x) >= 4)
  {
    double m, mx = rtodbl(x);
    m = (bit_accuracy_mul(l,LOG2) + mx)/4;
    n = (long)(1+m*m/mx);
    z = divsr(-n, addsr(n<<1,x));
    for (i=n-1; i >= 1; i--)
      z = divsr(-i, addrr(addsr(i<<1,x), mulsr(i,z))); /* -1 / (2 + z + x/i) */
    return divrr(addrr(real_1(l),z), mulrr(expx, x));
  }
  else
  {
    GEN S, t, H, run = real_1(l+1);
    n = -bit_accuracy(l)-1;
    x = rtor(x, l+1);
    S = z = t = H = run;
    for (i = 2; expo(t) - expo(S) >= n; i++)
    {
      H = addrr(H, divrs(run,i)); /* H = sum_{i=1} 1/i */
      z = divrs(mulrr(x,z), i);   /* z = sum_{i=1} x^(i-1)/i */
      t = mulrr(z, H); S = addrr(S, t);
    }
    return subrr(mulrr(x, divrr(S,expx)), addrr(mplog(x), mpeuler(l)));
  }
}

/* assume x != 0 */
GEN
incgam2(GEN s, GEN x, long prec)
{
  GEN b, x_s, S, y;
  long l, n, i;
  pari_sp av = avma, av2, avlim;
  double m,mx;

  if (typ(x) != t_REAL) x = gtofp(x, prec);
  if (gcmp0(s)) {
    if (typ(x) == t_REAL && signe(x) > 0)
      return gerepileuptoleaf(av, incgam2_0(x, mpexp(x)));
  }
  if (typ(x) == t_COMPLEX)
  {
    double a = rtodbl(gel(x,1));
    double b = rtodbl(gel(x,2));
    l = precision(x);
    mx = sqrt(a*a + b*b);
  }
  else
  {
    l = lg(x);
    mx = fabs(rtodbl(x));
  }
  m = (bit_accuracy_mul(l,LOG2) + mx)/4;
  n = (long)(1+m*m/mx);
  i = typ(s);
  if (i == t_REAL) b = addsr(-1,s);
  else
  { /* keep b integral : final powering more efficient */
    GEN z = gtofp(s, prec);
    b = (i == t_INT)? addsi(-1,s): gaddsg(-1,z);
    s = z;
  }
  y = gmul(gexp(gneg(x), prec), gpow(x,b,prec));
  x_s = gsub(x, s);
  av2 = avma; avlim = stack_lim(av2,3);
  S = gdiv(gaddsg(-n,s), gaddgs(x_s,n<<1));
  for (i=n-1; i>=1; i--)
  {
    S = gdiv(gaddsg(-i,s), gadd(gaddgs(x_s,i<<1),gmulsg(i,S)));
    if (low_stack(avlim,stack_lim(av2,3)))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"incgam2");
      S = gerepileupto(av2, S);
    }
  }
  return gerepileupto(av, gmul(y, gaddsg(1,S)));
}

/* use exp(-x) * (x^s/s) * sum_{k >= 0} x^k / prod(i=1,k, s+i)
 * ( =  exp(-x) * x^s * sum_{k >= 0} x^k / k!(k+s) but the above is more
 * efficient ) */
GEN
incgamc(GEN s, GEN x, long prec)
{
  GEN b, S, t, y;
  long l, n, i;
  pari_sp av = avma, av2, avlim;

  if (typ(x) != t_REAL) x = gtofp(x, prec);
  if (gcmp0(x)) return rcopy(x);

  l = precision(x); n = -bit_accuracy(l)-1;
  i = typ(s); b = s;
  if (i != t_REAL)
  {
    s = gtofp(s, prec);
    if (i != t_INT) b = s;
  }
  av2 = avma; avlim = stack_lim(av2,3);
  S = t = real_1(l);
  for (i=1; gexpo(S) >= n; i++)
  {
    S = gdiv(gmul(x,S), gaddsg(i,s)); /* x^i / ((s+1)...(s+i)) */
    t = gadd(S,t);
    if (low_stack(avlim,stack_lim(av2,3)))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"incgamc");
      gerepileall(av2, 2, &S, &t);
    }
  }
  y = gdiv(gmul(gexp(gneg(x),prec), gpow(x,b,prec)), s);
  return gerepileupto(av, gmul(y,t));
}

/* If g != NULL, assume that g=gamma(s,prec). */
GEN
incgam0(GEN s, GEN x, GEN g, long prec)
{
  pari_sp av = avma;
  long es, e;
  GEN z;

  if (gcmp0(x)) { avma = av; return g? gcopy(g): ggamma(s,prec); }
  es = gexpo(s); e = max(es, 0);
  if (gsigne(real_i(s)) <= 0 || gexpo(x) > e)
    z = incgam2(s,x,prec);
  else
  {
    if (es < 0) {
      long l = precision(s);
      if (!l) l = prec;
      prec = l + nbits2nlong(-es) + 1;
      s = gtofp(s, prec);
      x = gtofp(x, prec);
    }
    z = gadd(g? g: ggamma(s,prec), gneg(incgamc(s,x,prec)));
  }
  return gerepileupto(av, z);
}

GEN
incgam(GEN s, GEN x, long prec) { return incgam0(s, x, NULL, prec); }

GEN
eint1(GEN x, long prec)
{
  long l, n, i;
  pari_sp av = avma;
  GEN p1, t, S, y;

  if (typ(x) != t_REAL) {
    x = gtofp(x, prec);
    if (typ(x) != t_REAL) pari_err(impl,"non-real argument in eint1");
  }
  if (signe(x) >= 0) return gerepileuptoleaf(av, incgam2_0(x, mpexp(x)));
  /* rewritten from code contributed by Manfred Radimersky */
  l  = lg(x);
  n  = bit_accuracy(l);
  y  = negr(x);
  if (cmprs(y, (3*n)/4) < 0) {
    p1 = t = S = y;
    for (i = 2; expo(t) - expo(S) >= -n; i++) {
      p1 = mulrr(y, divrs(p1, i));
      t = divrs(p1, i); S = addrr(S, t);
    }
    y  = addrr(S, addrr(mplog(y), mpeuler(l)));
  } else {
    p1 = divsr(1, y);
    t = S = real_1(l);
    for (i = 1; expo(t) - expo(S) >= -n; i++) {
      t = mulrr(p1, mulrs(t, i)); S = addrr(S, t);
    }
    y  = mulrr(S, mulrr(p1, mpexp(y)));
  }
  return gerepileuptoleaf(av, negr(y));
}

GEN
veceint1(GEN C, GEN nmax, long prec)
{
  long i, n, nstop, nmin, G, chkpoint;
  pari_sp av, av1;
  GEN y, e1, e2, eC, F0, unr;

  if (!nmax) return eint1(C,prec);
  if (typ(nmax) != t_INT) pari_err(typeer,"veceint1");

  if (signe(nmax)<=0) return cgetg(1,t_VEC);
  if (DEBUGLEVEL>1) fprintferr("Entering veceint1:\n");
  if (typ(C) != t_REAL || lg(C) > prec) {
    C = gtofp(C, prec);
    if (typ(C) != t_REAL) pari_err(typeer,"veceint1");
  }
  if (signe(C) <= 0) pari_err(talker,"negative or zero constant in veceint1");

  n = itos(nmax); y = cgetg(n+1,t_VEC);
  for(i=1; i<=n; i++) gel(y,i) = cgetr(prec);
  av = avma; G = expo(C);
  if (G >= 0) nstop = n;
  else
  {
    nstop = itos(gceil(divsr(4,C))); /* >= 4 ~ 4 / C */
    if (nstop > n) nstop = n;
  }

  eC = mpexp(C);
  e1 = gpowgs(eC, -n);
  e2 = gpowgs(eC, 10);
  unr = real_1(prec);
  av1 = avma;
  if(DEBUGLEVEL>1) fprintferr("nstop = %ld\n",nstop);
  if (nstop == n) goto END;

  G = -bit_accuracy(prec);
  F0 = gel(y,n); chkpoint = n;
  affrr(eint1(mulsr(n,C),prec), F0);
  nmin = n;
  for(;;)
  {
    GEN minvn = divrs(unr,-n), My = subrr(minvn,C);
    GEN mcn   = divrs(C,  -n), Mx = mcn;
    GEN t = divrs(e1,-n), D = mkvec2( t, mulrr(My,t) );
    long a, k, cD = 2; /* cD = #D */

    /* D = [ e1/-n, (-1/n-C) * (e1/-n) ] */
    nmin -= 10; if (nmin < nstop) nmin = nstop;
    My = addrr(My, minvn);
    if (DEBUGLEVEL>1 && n < chkpoint)
      { fprintferr("%ld ",n) ; chkpoint -= nstop/20; }
    for (a=1,n--; n>=nmin; n--,a++)
    {
      GEN F = F0, den = stor(-a, prec);
      for (k=1;;)
      {
        GEN add;
	if (k > cD)
	{
	  GEN z = addrr(mulrr(My, gel(D,cD)), mulrr(Mx,gel(D,cD-1)));
          Mx = addrr(Mx,mcn);
	  My = addrr(My,minvn);
          D = shallowconcat(D, z); cD = k;
          /* My = -C - k/n,  Mx = -C k/n */
	}
	add = mulrr(den, gel(D,k));
	if (expo(add) < G) { affrr(F,gel(y,n)); break; }
	F = addrr(F,add); k++;
	den = mulrs(divrs(den, k), -a);
        /* den = prod(i=1,k, -a/i)*/
      }
    }
    avma = av1; F0 = gel(y, ++n);
    if (n <= nstop) break;
    affrr(mulrr(e1,e2), e1);
  }
END:
  affrr(eC, e1);
  for(i=1;; i++)
  { /* e1 = exp(iC) */
    affrr(incgam2_0(mulsr(i,C), e1), gel(y,i));
    if (i == nstop) break;
    affrr(mulrr(e1, eC), e1); avma = av1;
  }
  if (DEBUGLEVEL>1) fprintferr("\n");
  avma = av; return y;
}

GEN
gerfc(GEN x, long prec)
{
  pari_sp av;
  GEN z, sqrtpi;

  if (typ(x) != t_REAL) {
    x = gtofp(x, prec);
    if (typ(x) != t_REAL) pari_err(typeer,"erfc");
  }
  if (!signe(x)) return real_1(prec);
  av = avma; sqrtpi = sqrtr(mppi(lg(x)));
  z = incgam0(ghalf, gsqr(x), sqrtpi, prec);
  z = divrr(z, sqrtpi);
  if (signe(x) < 0) z = subsr(2,z);
  return gerepileupto(av,z);
}

/***********************************************************************/
/**								      **/
/**		        FONCTION ZETA DE RIEMANN                      **/
/**								      **/
/***********************************************************************/
static const double log2PI = 1.83787706641;

static double
get_xinf(double beta)
{
  const double maxbeta = 0.06415003; /* 3^(-2.5) */
  double x0, y0, x1;

  if (beta < maxbeta) return beta + pow(3*beta, 1.0/3.0);
  x0 = beta + PI/2.0;
  for(;;)
  {
    y0 = x0*x0;
    x1 = (beta+atan(x0)) * (1+y0) / y0 - 1/x0;
    if (0.99*x0 < x1) return x1;
    x0 = x1;
  }
}
/* optimize for zeta( s + it, prec ) */
static void
optim_zeta(GEN S, long prec, long *pp, long *pn)
{
  double s, t, alpha, beta, n, B;
  long p;
  if (typ(S) == t_REAL) {
    s = rtodbl(S);
    t = 0.;
  } else {
    s = rtodbl(gel(S,1));
    t = fabs( rtodbl(gel(S,2)) );
  }

  B = bit_accuracy_mul(prec, LOG2);
  if (s <= 0) /* may occur if S ~ 0, and we don't use the func. eq. */
  { /* TODO: the crude bounds below are generally valid. Optimize ? */
    double l,l2, la = 1.; /* heuristic */
    if (dnorm(s-1,t) < 0.1) /* |S - 1|^2 < 0.1 */
      l2 = -(s - 0.5);
    else
    {
      double rlog, ilog; dcxlog(s-1,t, &rlog,&ilog);
      l2 = (s - 0.5)*rlog - t*ilog; /* = Re( (S - 1/2) log (S-1) ) */
    }
    l = (B - l2 + s*log2PI) / (2. * (1.+ log((double)la)));
    l2 = dabs(s, t)/2;
    if (l < l2) l = l2;
    p = (long) ceil(l); if (p < 2) p = 2;
    l2 = (p + s/2. - .25);
    n = 1 + dabs(l2, t/2) * la / PI;
  }
  else if (t != 0)
  {
    double sn = dabs(s, t), L = log(sn/s);
    alpha = B - 0.39 + L + s*(log2PI - log(sn));
    beta = (alpha+s)/t - atan(s/t);
    if (beta <= 0)
    {
      if (s >= 1.0)
      {
	p = 0;
	n = exp((B - LOG2 + L) / s);
      }
      else
      {
	p = 1;
        n = dabs(s + 1, t) / (2*PI);
      }
    }
    else
    {
      beta = 1.0 - s + t * get_xinf(beta);
      if (beta > 0)
      {
	p = (long)ceil(beta / 2.0);
	n = dabs(s + 2*p-1, t) / (2*PI);
      }
      else
      {
	p = 0;
	n = exp((B - LOG2 + L) / s);
      }
    }
  }
  else
  {
    double sn = fabs(s);
    beta = B + 0.61 + s*(log2PI - log(s));
    if (beta > 0)
    {
      p = (long)ceil(beta / 2.0);
      n = fabs(s + 2*p-1)/(2*PI);
    }
    else
    {
      p = 0;
      n = exp((B - LOG2 + log(sn/s)) / s);
    }
  }
  *pp = p;
  *pn = (long)ceil(n);
  if (DEBUGLEVEL) fprintferr("lim, nn: [%ld, %ld]\n", *pp, *pn);
}

#if 0
static GEN
czeta(GEN s0, long prec)
{
  int funeq = 0;
  long n, p, n1, i, i2;
  pari_sp av;
  GEN s, y, z, res, sig, ms, p1, p2, p3, p31, pitemp;

  s = trans_fix_arg(&prec,&s0,&sig,&av,&res);

  if (signe(sig) <= 0 || expo(sig) <= -2)
  {
    if (typ(s0) == t_INT)
    {
      gaffect(szeta(itos(s0), prec), res);
      avma = av; return res;
    }
    funeq = 1; s = gsub(gen_1,s);
  }
  optim_zeta(s, prec, &p, &n);

  n1 = (n < 46340)? n*n: 0;
  y=gen_1; ms=gneg_i(s); p1=cgetr(prec+1); p2=gen_1;
  for (i=2; i<=n; i++)
  {
    affsr(i,p1); p2 = gexp(gmul(ms,mplog(p1)), prec+1);
    y = gadd(y,p2);
  }
  mpbern(p,prec); p31=cgetr(prec+1); z=gen_0;
  for (i=p; i>=1; i--)
  {
    i2 = i<<1;
    p1 = gmul(gaddsg(i2-1,s),gaddsg(i2,s));
    p1 = divgsns(p1, i2);
    p1 = n1? gdivgs(p1,n1): gdivgs(gdivgs(p1,n),n);
    p3 = bern(i);
    if (bernzone[2]>prec+1) { affrr(p3,p31); p3=p31; }
    z = gadd(divrs(p3,i),gmul(p1,z));
  }
  p1 = gsub(gdivsg(n,gsubgs(s,1)),ghalf);
  z = gmul(gadd(p1,gmul(s,gdivgs(z,n<<1))),p2);
  y = gadd(y,z);
  if (funeq)
  {
    pitemp=mppi(prec+1); setexpo(pitemp,2);
    y = gmul(gmul(y,ggamma(s,prec+1)), gpow(pitemp,ms,prec+1));
    setexpo(pitemp, 0);
    y = gmul2n(gmul(y, gcos(gmul(pitemp,s),prec+1)),1);
  }
  gaffect(y,res); avma = av; return res;
}
#endif

/* 1/zeta(n) using Euler product. Assume n > 0.
 * if (lba != 0) it is log(bit_accuracy) we _really_ require */
GEN
inv_szeta_euler(long n, double lba, long prec)
{
  GEN z, res = cgetr(prec);
  pari_sp av = avma, avlim = stack_lim(av, 1);
  byteptr d =  diffptr + 2;
  double A = n / (LOG2*BITS_IN_LONG), D;
  ulong p, lim;

  if (n > bit_accuracy(prec)) return real_1(prec);
  if (!lba) lba = bit_accuracy_mul(prec, LOG2);
  D = exp((lba - log(n-1)) / (n-1));
  lim = 1 + (ulong)ceil(D);
  maxprime_check(lim);

  prec++;
  z = gsub(gen_1, real2n(-n, prec));
  for (p = 3; p <= lim;)
  {
    long l = prec + 1 - (long)floor(A * log(p));
    GEN h;

    if (l < 3)         l = 3;
    else if (l > prec) l = prec;
    h = divrr(z, rpowuu(p, (ulong)n, l));
    z = subrr(z, h);
    if (low_stack(avlim, stack_lim(av,1)))
    {
      if (DEBUGMEM>1) pari_warn(warnmem,"inv_szeta_euler, p = %lu/%lu", p,lim);
      affrr(z, res); avma = av;
    }
    NEXT_PRIME_VIADIFF(p,d);
  }
  affrr(z, res); avma = av; return res;
}

/* assume n even > 0, if iz != NULL, assume iz = 1/zeta(n) */
GEN
bernreal_using_zeta(long n, GEN iz, long prec)
{
  long l = prec + 1;
  GEN z;

  if (!iz) iz = inv_szeta_euler(n, 0., l);
  z = divrr(mpfactr(n, l), mulrr(gpowgs(Pi2n(1, l), n), iz));
  setexpo(z, expo(z) + 1); /* 2 * n! * zeta(n) / (2Pi)^n */
  if ((n & 3) == 0) setsigne(z, -1);
  return z;
}

/* assume n even > 0. Faster than standard bernfrac for n >= 6 */
GEN
bernfrac_using_zeta(long n)
{
  pari_sp av = avma;
  GEN iz, a, d, D = divisors(utoipos( n/2 ));
  long i, prec, l = lg(D);
  double t, u;

  d = utoipos(6); /* 2 * 3 */
  for (i = 2; i < l; i++) /* skip 1 */
  { /* Clausen - von Staudt */
    ulong p = 2*itou(gel(D,i)) + 1;
    if (uisprime(p)) d = muliu(d, p);
  }
  /* 1.712086 = ??? */
  t = log( gtodouble(d) ) + (n + 0.5) * log(n) - n*(1+log2PI) + 1.712086;
  u = t / (LOG2*BITS_IN_LONG); prec = (long)ceil(u);
  prec += 3;
  iz = inv_szeta_euler(n, t, prec);
  a = roundr( mulir(d, bernreal_using_zeta(n, iz, prec)) );
  return gerepilecopy(av, mkfrac(a, d));
}

/* y = binomial(n,k-2). Return binomial(n,k) */
static GEN
next_bin(GEN y, long n, long k)
{
  y = divrs(mulrs(y, n-k+2), k-1);
  return divrs(mulrs(y, n-k+1), k);
}

/* assume k > 1 odd */
static GEN
szeta_odd(long k, long prec)
{
  long kk, n, li = -(1+bit_accuracy(prec));
  pari_sp av = avma, av2, limit;
  GEN y, p1, qn, z, q, pi2 = Pi2n(1, prec), binom= real_1(prec+1);

  q = mpexp(pi2); kk = k+1; /* >= 4 */
  y = NULL; /* gcc -Wall */
  if ((k&3)==3)
  {
    for (n=0; n <= kk>>1; n+=2)
    {
      p1 = mulrr(bernreal(kk-n,prec),bernreal(n,prec));
      if (n) { binom = next_bin(binom,kk,n); setlg(binom,prec+1); }
      p1 = mulrr(binom,p1);
      if (n == kk>>1) setexpo(p1, expo(p1)-1);
      if ((n>>1)&1) setsigne(p1,-signe(p1));
      y = n? addrr(y,p1): p1;
    }
    y = mulrr(divrr(gpowgs(pi2,k),mpfactr(kk,prec)),y);

    av2 = avma; limit = stack_lim(av2,1);
    qn = gsqr(q); z = ginv( addrs(q,-1) );
    for (n=2; ; n++)
    {
      p1 = ginv( mulir(powuu(n,k),addrs(qn,-1)) );

      z = addrr(z,p1); if (expo(p1)< li) break;
      qn = mulrr(qn,q);
      if (low_stack(limit,stack_lim(av2,1)))
      {
        if (DEBUGMEM>1) pari_warn(warnmem,"szeta");
        gerepileall(av2,2, &z, &qn);
      }
    }
    setexpo(z, expo(z)+1);
    y = addrr(y,z); setsigne(y,-signe(y));
  }
  else
  {
    GEN p2 = divrs(pi2, k-1);
    for (n=0; n <= k>>1; n+=2)
    {
      p1 = mulrr(bernreal(kk-n,prec),bernreal(n,prec));
      if (n) binom = next_bin(binom,kk,n);
      p1 = mulrr(binom,p1);
      p1 = mulsr(kk-(n<<1),p1);
      if ((n>>1)&1) setsigne(p1,-signe(p1));
      y = n? addrr(y,p1): p1;
    }
    y = mulrr(divrr(gpowgs(pi2,k),mpfactr(kk,prec)),y);
    y = divrs(y,k-1);
    av2 = avma; limit = stack_lim(av2,1);
    qn = q; z=gen_0;
    for (n=1; ; n++)
    {
      p1=mulir(powuu(n,k),gsqr(addrs(qn,-1)));
      p1=divrr(addrs(mulrr(qn,addsr(1,mulsr(n<<1,p2))),-1),p1);

      z=addrr(z,p1); if (expo(p1) < li) break;
      qn=mulrr(qn,q);
      if (low_stack(limit,stack_lim(av2,1)))
      {
        if (DEBUGMEM>1) pari_warn(warnmem,"szeta");
        gerepileall(av2,2, &z, &qn);
      }
    }
    setexpo(z, expo(z)+1);
    y = subrr(y,z);
  }
  return gerepileuptoleaf(av, y);
}

/* assume k > 0 even. Return B_k */
static GEN
single_bern(long k, long prec)
{
  pari_sp av;
  GEN B;
  if (OK_bern(k >> 1, prec)) B = bernreal(k, prec);
  else if (k * (log(k) - 2.83) > bit_accuracy_mul(prec, LOG2))
    B = bernreal_using_zeta(k, NULL, prec);
  else
  {
    B = cgetr(prec);
    av = avma; gaffect(bernfrac(k), B);
    avma = av;
  }
  return B;
}

/* assume k != 1 */
GEN
szeta(long k, long prec)
{
  pari_sp av = avma;
  GEN y;

  /* treat trivial cases */
  if (!k) { y = real2n(-1, prec); setsigne(y,-1); return y; }
  if (k < 0)
  {
    if ((k&1) == 0) return gen_0;
    /* the one value such that k < 0 and 1 - k < 0, due to overflow */
    if ((ulong)k == (HIGHBIT | 1))
      pari_err(talker, "too large negative arg %ld in zeta", k);
    k = 1-k;
    y = single_bern(k, prec);
    return gerepileuptoleaf(av, divrs(y, -k));
  }
  if (k > bit_accuracy(prec)+1) return real_1(prec);
  if ((k&1) == 0)
  {
    if (!OK_bern(k >> 1, prec)
        && (k * (log(k) - 2.83) > bit_accuracy_mul(prec, LOG2)))
      y = ginv( inv_szeta_euler(k, 0, prec) ); /* would use zeta above */
    else
    {
      y = mulrr(gpowgs(Pi2n(1, prec), k), single_bern(k, prec));
      y = divrr(y, mpfactr(k,prec));
      y[1] = evalsigne(1) | evalexpo(expo(y)-1);
    }
    return gerepileuptoleaf(av, y);
  }
  /* k > 1 odd */
  if (k * log(k) > bit_accuracy_mul(prec, LOG2)) /* heuristic */
    return gerepileuptoleaf(av, ginv( inv_szeta_euler(k, 0, prec) ));
  return szeta_odd(k, prec);
}

/* return x^n, assume n > 0 */
static long
pows(long x, long n)
{
  long i, y = x;
  for (i=1; i<n; i++) y *= x;
  return y;
}

/* return n^-s, n > 1 odd. tab[q] := q^-s, q prime power */
static GEN
n_s(ulong n, GEN *tab)
{
  byteptr d =  diffptr + 2;
  GEN x = NULL;
  long p, e;

  for (p = 3; n > 1; )
  {
    e = u_lvalrem(n, p, &n);
    if (e)
    {
      GEN y = tab[pows(p,e)];
      if (!x) x = y; else x = gmul(x,y);
    }
    NEXT_PRIME_VIADIFF_CHECK(p,d);
  }
  return x;
}

/* s0 a t_INT, t_REAL or t_COMPLEX.
 * If a t_INT, assume it's not a trivial case (i.e we have s0 > 1, odd)
 * */
GEN
czeta(GEN s0, long prec)
{
  GEN s, u, a, y, res, tes, sig, invn2, unr;
  GEN sim, *tab, tabn;
  ulong p, sqn;
  long i, nn, lim, lim2, ct;
  pari_sp av, av2 = avma, avlim;
  int funeq = 0;
  byteptr d;

  if (DEBUGLEVEL>2) (void)timer2();
  s = trans_fix_arg(&prec,&s0,&sig,&av,&res);
  if (gcmp0(s)) { y = gneg(ghalf); goto END; }
  if (gexpo(gsub(s, gen_1)) < -5 ||
      (gexpo(s) > -5 && (signe(sig) <= 0 || expo(sig) < -1)))
  { /* s <--> 1-s */
    if (typ(s0) == t_INT)
    {
      i = itos(s0); avma = av2;
      return szeta(i, prec);
    }
    funeq = 1; s = gsub(gen_1, s); sig = real_i(s);
  }
  if (gcmpgs(sig, bit_accuracy(prec) + 1) > 0) { y = gen_1; goto END; }
  optim_zeta(s, prec, &lim, &nn);
  maxprime_check((ulong)nn);
  prec++; unr = real_1(prec); /* one extra word of precision */

  tab = (GEN*)cgetg(nn, t_VEC); /* table of q^(-s), q = p^e */
  d = diffptr + 1;
  if (typ(s0) == t_INT)
  { /* no explog for 1/p^s */
    ulong k = itou(s0);
    for (p=2; p < (ulong)nn;) {
      tab[p] = divrr(unr, rpowuu(p, k, prec));
      NEXT_PRIME_VIADIFF(p,d);
    }
    a = divrr(unr, rpowuu((ulong)nn, k, prec));
  }
  else
  { /* general case */
    GEN ms = gneg(s), rp = cgetr(prec);
    for (p=2; p < (ulong)nn;)
    {
      affur(p, rp);
      tab[p] = gexp(gmul(ms, mplog(rp)), prec);
      NEXT_PRIME_VIADIFF(p,d);
    }
    affsr(nn, rp);
    a = gexp(gmul(ms, mplog(rp)), prec);
  }
  sqn = (ulong)sqrt(nn-1.); maxprime_check(sqn);
  d = diffptr + 2; /* fill in odd prime powers */
  for (p=3; p <= sqn; )
  {
    ulong oldq = p, q = p*p;
    while (q<(ulong)nn) { tab[q] = gmul(tab[p], tab[oldq]); oldq = q; q *= p; }
    NEXT_PRIME_VIADIFF(p,d);
  }
  if (DEBUGLEVEL>2) msgtimer("tab[q^-s] from 1 to N-1");

  tabn = cgetg(nn, t_VECSMALL); ct = 0;
  for (i = nn-1; i; i>>=1) tabn[++ct] = (i-1)>>1;
  sim = y = unr;
  for (i=ct; i > 1; i--)
  {
    long j;
    pari_sp av2 = avma;
    for (j=tabn[i]+1; j<=tabn[i-1]; j++)
      sim = gadd(sim, n_s(2*j+1, tab));
    sim = gerepileupto(av2, sim);
    y = gadd(sim, gmul(tab[2],y));
  }
  y = gadd(y, gmul2n(a,-1));
  if (DEBUGLEVEL>2) msgtimer("sum from 1 to N-1");

  invn2 = divri(unr, mulss(nn,nn)); lim2 = lim<<1;
  tes = bernreal(lim2, prec);
  if (typ(s0) == t_INT)
  {
    av2 = avma; avlim = stack_lim(av2,3);
    for (i=lim2-2; i>=2; i-=2)
    { /* using single prec (when (s0 + i) < 2^31) not faster (even at \p28) */
      u = mulri(mulrr(tes,invn2), mulii(addsi(i,s0), addsi(i-1,s0)));
      tes = addrr(bernreal(i,prec), divrsns(u, i+1)); /* u / (i+1)(i+2) */
      if (low_stack(avlim,stack_lim(av2,3)))
      {
        if(DEBUGMEM>1) pari_warn(warnmem,"czeta");
        tes = gerepileuptoleaf(av2, tes);
      }
    }
    u = gmul(gmul(tes,invn2), gmul2n(mulii(s0, addsi(-1,s0)), -1));
    tes = gmulsg(nn, gaddsg(1, u));
  }
  else /* typ(s0) != t_INT */
  {
    GEN s1, s2, s3, s4, s5;
    s1 = gsub(gmul2n(s,1), unr);
    s2 = gmul(s, gsub(s,unr));
    s3 = gmul2n(invn2,3);
    av2 = avma; avlim = stack_lim(av2,3);
    s4 = gmul(invn2, gmul2n(gaddsg(4*lim-2,s1),1));
    s5 = gmul(invn2, gadd(s2, gmulsg(lim2, gaddgs(s1, lim2))));
    for (i = lim2-2; i>=2; i -= 2)
    {
      s5 = gsub(s5, s4);
      s4 = gsub(s4, s3);
      tes = gadd(bernreal(i,prec), divgsns(gmul(s5,tes), i+1));
      if (low_stack(avlim,stack_lim(av2,3)))
      {
        if(DEBUGMEM>1) pari_warn(warnmem,"czeta");
        gerepileall(av2,3, &tes,&s5,&s4);
      }
    }
    u = gmul(gmul(tes,invn2), gmul2n(s2, -1));
    tes = gmulsg(nn, gaddsg(1, u));
  }
  if (DEBUGLEVEL>2) msgtimer("Bernoulli sum");
  /* y += tes n^(-s) / (s-1) */
  y = gadd(y, gmul(tes, gdiv(a, gsub(s, unr))));

END:
  if (funeq)
  {
    y = gmul(gmul(y, ggamma(gprec_w(s,prec),prec)),
             gpow(Pi2n(1,prec), gneg(s), prec));
    y = gmul2n(gmul(y, gcos(gmul(Pi2n(-1,prec),s), prec)), 1);
  }
  gaffect(y,res); avma = av; return res;
}

/* return P mod x^n where P is polynomial in x */
static GEN
pol_mod_xn(GEN P, long n)
{
  long j;
  GEN R = cgetg(n+2, t_POL);
  R[1] = evalvarn(0);
  for (j = 0; j < n; j++)
    R[j+2] = (long)polcoeff0(P, j, 0);
  return normalizepol_i(R, n+2);
}

/* compute the values of the twisted partial
   zeta function Z_f(a, c, s) for a in va */
GEN
twistpartialzeta(GEN p, GEN q, long f, long c, GEN va, GEN cff)
{
  long j, k, lva = lg(va)-1, N = lg(cff)-1;
  pari_sp av, av2, lim;
  GEN Ax, Cx, Bx, Dx, x = pol_x[0], y = pol_x[fetch_user_var("y")], eta;
  GEN cyc, psm, rep;

  cyc = gdiv(gsubgs(gpowgs(y, c), 1), gsubgs(y, 1));
  psm = polsym(cyc, degpol(cyc) - 1);
  eta = gmodulo(y, cyc);
  av = avma;
  Ax  = gsubgs(gpowgs(gaddgs(x, 1), f), 1);
  Ax  = gdiv(gmul(Ax, gpowgs(eta, f)), gsubsg(1, gpowgs(eta, f)));
  Ax  = gerepileupto(av, RgX_to_FqX(Ax, cyc, q));
  Cx  = Ax;
  Bx  = gen_1;
  av  = avma; lim = stack_lim(av, 1);
  for (j = 2; j <= N; j++) 
  {
    Bx = gadd(Bx, Cx);
    Bx = FpXQX_red(Bx, cyc, q);
    Cx = FpXQX_mul(Cx, Ax, cyc, q);
    Cx = pol_mod_xn(Cx, N);
    if (gcmp0(Cx)) break; 
    if (low_stack(lim, stack_lim(av, 1)))
    {
      if(DEBUGMEM>1) pari_warn(warnmem, "twistpartialzeta (1), j = %ld/%ld", j, N);
      gerepileall(av, 2, &Cx, &Bx);
    }
  }
  Bx  = lift(gmul(ginv(gsubsg(1, gpowgs(eta, f))), Bx));
  Bx  = gerepileupto(av, RgX_to_FqX(Bx, cyc, q));
  Cx = lift(gmul(eta, gaddsg(1, x)));
  Dx = pol_1[varn(x)];
  av2 = avma; lim = stack_lim(av2, 1);
  for (j = lva; j > 1; j--)
  {  
    GEN Ex;
    long e = va[j] - va[j-1];
    if (e == 1)
      Ex = Cx;
    else
      /* e is very small in general and actually very rarely different
	 to 1, it is always 1 for zetap (so it should be OK not to store 
	 them or to compute them in a smart way) */ 
      Ex = gpowgs(Cx, e);
    Dx = gaddsg(1, FpXQX_mul(Dx, Ex, cyc, q));
    if (low_stack(lim, stack_lim(av2, 1)))
    {
      if(DEBUGMEM>1) 
	pari_warn(warnmem, "twistpartialzeta (2), j = %ld/%ld", lva-j, lva);
      Dx = gerepileupto(av2, FpXQX_red(Dx, cyc, q));
    }
  }
  Dx = FpXQX_mul(Dx, Cx, cyc, q); /* va[1] = 1 */
  Bx = gerepileupto(av, FpXQX_mul(Dx, Bx, cyc, q));  
  rep = gen_0;
  av2 = avma; lim = stack_lim(av2, 1);
  for (k = 1; k <= N; k++)
  {
    GEN p2, ak = polcoeff_i(Bx, k, 0); 
    p2  = quicktrace(ak, psm);
    rep = modii(addii(rep, mulii(gel(cff, k), p2)), q);
    if (low_stack(lim, stack_lim(av2, 1)))
    {
      if(DEBUGMEM>1) pari_warn(warnmem, "twistpartialzeta (3), j = %ld/%ld", k, N);
      rep = gerepileupto(av2, rep);
    }
  }
  return rep;
}

#if 0
/* initialize the roots of unity for the computation
   of the Teichmuller character (also the values of f and c) */
GEN
init_teich(ulong p, GEN q, long prec)
{
  GEN vz, gp = utoipos(p);
  pari_sp av = avma;
  long j;

  if (p == 2UL)
    return NULL;
  else
  { /* primitive (p-1)-th root of 1 */
    GEN z, z0 = padicsqrtnlift(gen_1, utoipos(p-1), gener_Fp(gp), gp, prec);
    z = z0;
    vz = cgetg(p, t_VEC);
    for (j = 1; j < (long)p-2; j++)
    {
      gel(vz, umodiu(z, p)) = z; /* z = z0^i */
      z = modii(mulii(z, z0), q);
    }
    gel(vz, umodiu(z, p)) = z; /* z = z0^(p-2) */
    gel(vz,1) = gen_1; /* z0^(p-1) */
  }
  return gerepileupto(av, gcopy(vz));
}

/* compute phi^(m)_s(x); s must be an integer */
GEN
phi_ms(ulong p, GEN q, long m, GEN s, long x, GEN vz)
{
  long xp = x % p;
  GEN p1, p2;

  if (!xp) return gen_0;
  if (vz)
    p1 =gel(vz,xp); /* vz[x] = Teichmuller(x) */
  else 
    p1 = (x & 2)? gen_m1: gen_1;
  p1 = Fp_pow(p1, addis(s, m), q);
  p2 = Fp_pow(stoi(x), negi(s), q);
  return modii(mulii(p1,p2), q);
}

/* compute the first N coefficients of the Mahler expansion 
   of phi^m_s skipping the first one (which is zero) */
GEN
coeff_of_phi_ms(ulong p, GEN q, long m, GEN s, long N, GEN vz)
{
  GEN qs2 = shifti(q, -1), cff = zerovec(N);
  pari_sp av, lim;
  long k, j;

  av = avma; lim = stack_lim(av, 2);
  for (k = 1; k <= N; k++) 
  {
    gel(cff, k) = phi_ms(p, q, m, s, k, vz); 
    if (low_stack(lim, stack_lim(av, 2)))
    {
      if(DEBUGMEM>1) 
	pari_warn(warnmem, "coeff_of_phi_ms (1), k = %ld/%ld", N-k, N);
      cff = gerepileupto(av, gcopy(cff));
    }
  }
  for (j = N; j > 1; j--)
  {
    GEN b = subii(gel(cff, j), gel(cff, j-1));
    gel(cff, j) = centermodii(b, q, qs2);
    if (low_stack(lim, stack_lim(av, 2)))
    {
      if(DEBUGMEM>1) 
	pari_warn(warnmem, "coeff_of_phi_ms (2), j = %ld/%ld", N-j, N);
      cff = gerepileupto(av, gcopy(cff));
    }
  }
  for (k = 1; k < N; k++)
    for (j = N; j > k; j--)
    {
      GEN b = subii(gel(cff, j), gel(cff, j-1));
      gel(cff, j) = centermodii(b, q, qs2);
      if (low_stack(lim, stack_lim(av, 2)))
      {
	if(DEBUGMEM>1) 
	  pari_warn(warnmem, "coeff_of_phi_ms (3), (k,j) = (%ld,%ld)/%ld", 
	      k, N-j, N);
	cff = gerepileupto(av, gcopy(cff));
      }
    }
  k = N; while(gcmp0(gel(cff, k))) k--;
  setlg(cff, k+1);
  if (DEBUGLEVEL > 2) 
    fprintferr("  coeff_of_phi_ms: %ld coefficients kept out of %ld\n", 
	       k, N);
  return gerepileupto(av, cff);
}

static long
valfact(long N, ulong p)
{
  long f = 0;
  while (N > 1) 
  {
    N /= p;
    f += N;
  }
  return f;
}

static long
number_of_terms(ulong p, long prec)
{
  long N, f;

  if (prec == 0) return p;
  N = (long)((p-1)*prec + (p>>1)*(log2(prec)/log2(p)));
  N = p*(N/p);
  f = valfact(N, p);
  while (f > prec) 
  {
    N = p*(N/p) - 1;
    f -= u_lval(N+1, p);
  }
  while (f < prec) 
  {
    N = p*(N/p+1);
    f += u_lval(N, p);
  }
  return N;
}

static GEN
zetap(GEN s)
{
  ulong p;
  long N, f, c, prec = precp(s);
  pari_sp av = avma;
  GEN gp, q, vz, is, cff, val, va, cft;

  if (valp(s) < 0)
    pari_err(talker, "argument must be a p-adic integer");

  gp = gel(s,2); p = itou(gp);
  is = gtrunc(s);  /* make s an integer */

  N  = number_of_terms(p, prec? prec:1);
  q  = gpowgs(gp, prec? prec:1);

  /* initialize the roots of unity for the computation
     of the Teichmuller character (also the values of f and c) */
  if (DEBUGLEVEL > 1) fprintferr("zetap: computing (p-1)th roots of 1\n");
  vz = init_teich(p, q, prec? prec:1);
  if (p == 2UL) {  f = 4; c = 3; } else { f = (long)p; c = 2; }

  /* compute the first N coefficients of the Mahler expansion
     of phi^(-1)_s skipping the first one (which is zero) */
  if (DEBUGLEVEL > 1)
    fprintferr("zetap: computing Mahler expansion of phi^(-1)_s\n");
  cff = coeff_of_phi_ms(p, q, -1, is, N, vz);

  /* compute the coefficients of the power series corresponding
     to the twisted partial zeta function Z_f(a, c, s) for a in va */
  /* The line below looks a bit stupid but it is to keep the
     possibility of later adding p-adic Dirichlet L-functions */
  va = perm_identity(f - 1);
  if (DEBUGLEVEL > 1)
    fprintferr("zetap: computing values of twisted partial zeta functions\n");
  val = twistpartialzeta(gp, q, f, c, va, cff);

  /* sum over all a's the coefficients of the twisted
     partial zeta functions and integrate */
  if (DEBUGLEVEL > 1)
    fprintferr("zetap: multiplying by correcting factor\n");

  /* multiply by the corrective factor */
  cft = gsubgs(gmulsg(c, phi_ms(p, q, -1, is, c, vz)), 1);
  val = gdiv(val, cft);

  /* adjust the precision and return */
  return gerepileupto(av, cvtop(val, gp, prec));
}
#else
static GEN
hurwitz_p(GEN cache, GEN s, GEN x, GEN p, long prec)
{
  GEN S, x2, x2j, s_1 = gsubgs(s,1);
  long j, J = lg(cache)-2;
  x = ginv(gadd(x, zeropadic(p, prec)));
  x2 = gsqr(x);
  S = gmul2n(gmul(s_1, x), -1);
  x2j = gen_1;
  for (j = 0;; j++)
  {
    S = gadd(S, gmul(gel(cache, j+1), x2j));
    if (j == J) break;
    x2j = gmul(x2, x2j);
  }
  return gmul(gdiv(S, s_1), gexp(gmul(s_1, glog(x, 0)), 0));
}

static GEN
init_cache(long J, GEN s)
{
  GEN C = gen_1, cache = bernvec(J);
  long j;

  for (j = 1; j <= J; j++)
  { /* B_{2j} * binomial(1-s, 2j) */
    GEN t = gmul(gaddgs(s, 2*j-3), gaddgs(s, 2*j-2));
    C = gdiv(gmul(C, t), mulss(2*j, 2*j-1));
    gel(cache, j+1) = gmul(gel(cache, j+1), C);
  }
  return cache;
}

static GEN
zetap(GEN s)
{
  pari_sp av = avma;
  GEN cache, S, gp = gel(s,2);
  ulong a, p = itou(gp);
  long J, prec = valp(s) + precp(s);

  if (prec <= 0) prec = 1;
  if (p == 2) {
    J = ((long)(1+ceil((prec+1.)/2))) >> 1;
    cache = init_cache(J, s);
    S = gmul2n(hurwitz_p(cache, s, gmul2n(gen_1, -2), gen_2, prec), -1);
  } else {
    J = (prec+2) >> 1;
    cache = init_cache(J, s);
    S = gen_0;
    for (a = 1; a <= (p-1)>>1; a++)
      S = gadd(S, hurwitz_p(cache, s, gdivsg(a, gp), gp, prec));
    S = gdiv(gmul2n(S, 1), gp);
  }
  return gerepileupto(av, S);
}
#endif

GEN
gzeta(GEN x, long prec)
{
  if (gcmp1(x)) pari_err(talker, "argument equal to one in zeta");
  switch(typ(x))
  {
    case t_INT:
      if (is_bigint(x))
      {
        if (signe(x) > 0) return real_1(prec);
        if (signe(x) < 0 && mod2(x) == 0) return real_0(prec);
      }
      return szeta(itos(x),prec);

    case t_REAL: case t_COMPLEX:
      return czeta(x,prec);

    case t_INTMOD: pari_err(typeer,"gzeta");

    case t_PADIC:
      return zetap(x);

    case t_SER: pari_err(impl,"zeta of power series");
  }
  return transc(gzeta,x,prec);
}

/***********************************************************************/
/**                                                                   **/
/**                    FONCTIONS POLYLOGARITHME                       **/
/**                                                                   **/
/***********************************************************************/

/* validity domain contains .005 < |x| < 230
 * Li_m(x = e^z) = sum_n=0 zeta(m-n) z^n / n!
 *    with zeta(1) := H_m - log(-z) */
static GEN
cxpolylog(long m, GEN x, long prec)
{
  long li, i, n, bern_upto;
  pari_sp av=avma;
  GEN p1,z,h,q,s;
  int real;

  if (gcmp1(x)) return szeta(m,prec);
  real = typ(x) == t_REAL;

  z = glog(x,prec); h = gen_1;
  for (i=2; i<m; i++) h = gadd(h, ginv(utoipos(i)));
  h = gadd(h, gneg_i( glog(gneg_i(z),prec) ));

  bern_upto = m+50; mpbern(bern_upto,prec);
  q = gen_1; s = szeta(m,prec);
  for (n=1; n<=m+1; n++)
  {
    q = gdivgs(gmul(q,z),n);
    if (n == m-1) {
      p1 = gmul(h, q);
      if (real) p1 = real_i(p1);
    }
    else
      p1 = gmul(szeta(m-n,prec), real? real_i(q): q);
    s = gadd(s, p1);
  }

  z = gsqr(z); li = -(bit_accuracy(prec)+1);
  for(n = m+3;; n += 2)
  {
    GEN zet = szeta(m-n,prec);
    q = divgsns(gmul(q,z), n-1);
    s = gadd(s, gmul(zet, real? real_i(q): q));
    if (gexpo(q) + expo(zet) < li) break;
    if (n >= bern_upto) { bern_upto += 50; mpbern(bern_upto,prec); }
  }
  return gerepileupto(av,s);
}

GEN
polylog(long m, GEN x, long prec)
{
  long l, e, i, G, sx;
  pari_sp av, av1, limpile;
  GEN X, Xn, z, p1, p2, y;

  if (m<0) pari_err(talker,"negative index in polylog");
  if (!m) return gneg(ghalf);
  if (gcmp0(x)) return gcopy(x);
  av = avma;
  if (m==1)
    return gerepileupto(av, gneg(glog(gsub(gen_1,x), prec)));

  l = precision(x);
  if (!l) { l=prec; x=gmul(x, real_1(l)); }
  e = gexpo(gnorm(x)); if (!e || e== -1) return cxpolylog(m,x,prec);
  X = (e > 0)? ginv(x): x;
  G = -bit_accuracy(l);
  av1=avma; limpile=stack_lim(av1,1);
  y = Xn = X;
  for (i=2; ; i++)
  {
    Xn = gmul(X,Xn); p2 = gdiv(Xn,powuu(i,m));
    y = gadd(y,p2);
    if (gexpo(p2) <= G) break;

    if (low_stack(limpile, stack_lim(av1,1)))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"polylog");
      gerepileall(av1,2, &y, &Xn);
    }
  }
  if (e < 0) return gerepileupto(av, y);

  sx = gsigne(imag_i(x));
  if (!sx)
  {
    if (m&1) sx = gsigne(gsub(gen_1, real_i(x)));
    else     sx = - gsigne(real_i(x));
  }
  z = pureimag( divri(mppi(l), mpfact(m-1)) );
  setsigne(z[2], sx);

  if (m == 2)
  { /* same formula as below, written more efficiently */
    y = gneg_i(y);
    if (typ(x) == t_REAL && signe(x) < 0)
      p1 = logr_abs(x);
    else
      p1 = gsub(glog(x,l), z);
    p1 = gmul2n(gsqr(p1), -1); /* = (log(-x))^2 / 2 */
    
    p1 = gadd(p1, divrs(gsqr(mppi(l)), 6));
    p1 = gneg_i(p1);
  }
  else
  {
    GEN logx = glog(x,l), logx2 = gsqr(logx); p1 = gneg_i(ghalf);
    for (i=m-2; i>=0; i-=2)
      p1 = gadd(szeta(m-i,l), gmul(p1,gdivgs(logx2,(i+1)*(i+2))));
    if (m&1) p1 = gmul(logx,p1); else y = gneg_i(y);
    p1 = gadd(gmul2n(p1,1), gmul(z,gpowgs(logx,m-1)));
    if (typ(x) == t_REAL && signe(x) < 0) p1 = real_i(p1);
  }
  return gerepileupto(av, gadd(y,p1));
}

GEN
dilog(GEN x, long prec)
{
  return gpolylog(2, x, prec);
}

GEN
polylogd0(long m, GEN x, long flag, long prec)
{
  long k, l, fl, m2;
  pari_sp av;
  GEN p1,p2,p3,y;

  m2=m&1; av=avma;
  if (gcmp0(x)) return gcopy(x);
  if (gcmp1(x) && m>=2) return m2?szeta(m,prec):gen_0;
  l = precision(x);
  if (!l) { l=prec; x=gmul(x,real_1(l)); }
  p1 = gabs(x,prec); fl=0;
  if (expo(p1) >= 0) { x=ginv(x); p1=gabs(x,prec); fl=!m2; }

  p1=gneg_i(glog(p1,prec)); p2=gen_1;
  y=polylog(m,x,prec); y = m2? real_i(y): imag_i(y);
  for (k=1; k<m; k++)
  {
    p2=gdivgs(gmul(p2,p1),k);
    p3=m2? real_i(polylog(m-k,x,prec)): imag_i(polylog(m-k,x,prec));
    y=gadd(y,gmul(p2,p3));
  }
  if (m2)
  {
    if (flag)
      p2 = gdivgs(gmul(p2,p1),-2*m);
    else
      p2 = gdivgs(gmul(glog(gnorm(gsub(gen_1,x)),prec),p2),2*m);
    y=gadd(y,p2);
  }
  if (fl) y = gneg(y);
  return gerepileupto(av, y);
}

GEN
polylogd(long m, GEN x, long prec)
{
  return polylogd0(m,x,0,prec);
}

GEN
polylogdold(long m, GEN x, long prec)
{
  return polylogd0(m,x,1,prec);
}

GEN
polylogp(long m, GEN x, long prec)
{
  long k, l, fl, m2;
  pari_sp av;
  GEN p1,y;

  m2=m&1; av=avma;
  if (gcmp0(x)) return gcopy(x);
  if (gcmp1(x) && m>=2) return m2?szeta(m,prec):gen_0;
  l=precision(x);
  if (!l) { l=prec; x=gmul(x,real_1(l)); }
  p1=gabs(x,prec); fl=0;
  if (expo(p1) >= 0) { x=ginv(x); p1=gabs(x,prec); fl=!m2; }

  p1=gmul2n(glog(p1,prec),1); mpbern(m>>1,prec);
  y=polylog(m,x,prec); y=m2?real_i(y):imag_i(y);

  if (m==1)
  {
    p1=gmul2n(p1,-2); y=gadd(y,p1);
  }
  else
  {
    GEN p2=gen_1, p3, p4, p5, p51=cgetr(prec);

    for (k=1; k<m; k++)
    {
      p2=gdivgs(gmul(p2,p1),k);
      if (!(k&1) || k==1)
      {
	if (k!=1)
	{
	  p5=(GEN)bern(k>>1);
	  if (bernzone[2]>prec) { affrr(p5,p51); p5=p51; }
	  p4=gmul(p2,p5);
	}
	else p4=gneg_i(gmul2n(p2,-1));
	p3=polylog(m-k,x,prec); p3=m2?real_i(p3):imag_i(p3);
	y=gadd(y,gmul(p4,p3));
      }
    }
  }
  if (fl) y = gneg(y);
  return gerepileupto(av, y);
}

GEN
gpolylog(long m, GEN x, long prec)
{
  long i, lx, n, v;
  pari_sp av = avma;
  GEN a, y, p1;

  if (m <= 0)
  {
    GEN t = mkpoln(2, gen_m1, gen_1); /* 1 - X */
    p1 = pol_x[0];
    for (i=2; i <= -m; i++)
      p1 = gmul(pol_x[0], gadd(gmul(t,derivpol(p1)), gmulsg(i,p1)));
    p1 = gdiv(p1, gpowgs(t,1-m));
    return gerepileupto(av, poleval(p1,x));
  }

  switch(typ(x))
  {
    case t_INT: case t_REAL: case t_FRAC: case t_COMPLEX: case t_QUAD:
      return polylog(m,x,prec);

    case t_POLMOD:
      p1=cleanroots(gel(x,1),prec); lx=lg(p1);
      for (i=1; i<lx; i++) gel(p1,i) = poleval(gel(x,2),gel(p1,i));
      y=cgetg(lx,t_COL);
      for (i=1; i<lx; i++) gel(y,i) = polylog(m,gel(p1,i),prec);
      return gerepileupto(av, y);

    case t_INTMOD: case t_PADIC: pari_err(impl, "padic polylogarithm");
    default:
      av = avma; if (!(y = toser_i(x))) break;
      if (!m) { avma = av; return gneg(ghalf); }
      if (m==1) return gerepileupto(av, gneg( glog(gsub(gen_1,y),prec) ));
      if (gcmp0(y)) return gcopy(y);
      v = valp(y);
      if (v <= 0) pari_err(impl,"polylog around a!=0");
      n = (lg(y)-3 + v) / v;
      a = zeroser(varn(y), lg(y)-2);
      for (i=n; i>=1; i--)
	a = gmul(y, gadd(a, gpowgs(utoipos(i),-m)));
      return gerepileupto(av, a);

    case t_VEC: case t_COL: case t_MAT:
      lx=lg(x); y=cgetg(lx,typ(x));
      for (i=1; i<lx; i++)
	gel(y,i) = gpolylog(m,gel(x,i),prec);
      return y;
  }
  pari_err(typeer,"gpolylog");
  return NULL; /* not reached */
}

void
gpolylogz(long m, GEN x, GEN y)
{
  long prec = precision(y);
  pari_sp av=avma;

  if (!prec) pari_err(infprecer,"gpolylogz");
  gaffect(gpolylog(m,x,prec),y); avma=av;
}

GEN
polylog0(long m, GEN x, long flag, long prec)
{
  switch(flag)
  {
    case 0: return gpolylog(m,x,prec);
    case 1: return polylogd(m,x,prec);
    case 2: return polylogdold(m,x,prec);
    case 3: return polylogp(m,x,prec);
    default: pari_err(flagerr,"polylog");
  }
  return NULL; /* not reached */
}

static GEN
upper_half(GEN x, long *prec)
{
  long tx = typ(x), l;
  if (tx == t_QUAD) { x = quadtoc(x, *prec); tx = typ(x); }
  if (tx != t_COMPLEX || gsigne(gel(x,2)) <= 0)
    pari_err(talker,"argument must belong to upper half-plane");
  l = precision(x); if (l) *prec = l;
  return x;
}

static GEN
qq(GEN x, long prec)
{
  long tx = typ(x);

  if (is_scalar_t(tx))
  {
    if (tx == t_PADIC) return x;
    x = upper_half(x, &prec);
    return gexp(gmul(mulcxI(x), Pi2n(1,prec)), prec); /* e(x) */
  }
  if (! ( x = toser_i(x)) ) pari_err(talker,"bad argument for modular function");
  return x;
}

static GEN
inteta(GEN q)
{
  long tx=typ(q);
  GEN p1,ps,qn,y;

  y=gen_1; qn=gen_1; ps=gen_1;
  if (tx==t_PADIC)
  {
    if (valp(q) <= 0) pari_err(talker,"non-positive valuation in eta");
    for(;;)
    {
      p1 = gneg_i(gmul(ps,gmul(q,gsqr(qn))));
      y = gadd(y,p1); qn = gmul(qn,q); ps = gmul(p1,qn);
      p1 = y;
      y = gadd(y,ps); if (gequal(p1,y)) return y;
    }
  }
  else
  {
    long l, v = 0; /* gcc -Wall */
    pari_sp av = avma, lim = stack_lim(av, 3);

    if (is_scalar_t(tx)) l = -bit_accuracy(precision(q));
    else
    {
      v = gvar(q); l = lg(q)-2; tx = 0;
      if (valp(q) <= 0) pari_err(talker,"non-positive valuation in eta");
    }
    for(;;)
    {
      p1 = gneg_i(gmul(ps,gmul(q,gsqr(qn))));
      /* qn = q^n
       * ps = (-1)^n q^(n(3n+1)/2)
       * p1 = (-1)^(n+1) q^(n(3n+1)/2 + 2n+1) */
      y = gadd(y,p1); qn = gmul(qn,q); ps = gmul(p1,qn);
      y = gadd(y,ps);
      if (tx)
        { if (gexpo(ps)-gexpo(y) < l) return y; }
      else
        { if (gval(ps,v) >= l) return y; }
      if (low_stack(lim, stack_lim(av,3)))
      {
        if(DEBUGMEM>1) pari_warn(warnmem,"eta");
        gerepileall(av, 3, &y, &qn, &ps);
      }
    }
  }
}

GEN
eta(GEN x, long prec)
{
  pari_sp av = avma;
  GEN q = qq(x,prec);
  return gerepileupto(av,inteta(q));
}

/* sqrt(3)/2 */
static GEN
sqrt32(long prec) { GEN z = sqrtr(stor(3, prec)); setexpo(z, -1); return z; }

/* exp(i x), x = k pi/12, assume 0 <= k < 24 */
static GEN
e12(ulong k, long prec)
{
  int s, sPi, sPiov2;
  GEN z, t;
  if (k >12) { s = 1; k = 24 - k; } else s = 0; /* x -> 2pi - x */
  if (k > 6) { sPi = 1; k = 12 - k; } else sPi = 0; /* x -> pi  - x */
  if (k > 3) { sPiov2 = 1; k = 6 - k; } else sPiov2 = 0; /* x -> pi/2 - x */
  z = cgetg(3, t_COMPLEX);
  switch(k)
  {
    case 0: gel(z,1) = icopy(gen_1); gel(z,2) = gen_0; break;
    case 1: t = gmul2n(addrs(sqrt32(prec), 1), -1);
      gel(z,1) = sqrtr(t);
      gel(z,2) = gmul2n(ginv(gel(z,1)), -2); break;

    case 2: gel(z,1) = sqrt32(prec);
            gel(z,2) = real2n(-1, prec); break;

    case 3: gel(z,1) = ginv( gsqrt(gen_2, prec) );
            gel(z,2) = mpcopy(gel(z,1)); break;
  }
  if (sPiov2) lswap(z[1], z[2]);
  if (sPi) setsigne(z[1], -signe(z[1]));
  if (s)   setsigne(z[2], -signe(z[2]));
  return z;
}

/* returns the true value of eta(x) for Im(x) > 0, using reduction to
 * standard fundamental domain */
GEN
trueeta(GEN x, long prec)
{
  long tx = typ(x);
  ulong Nmod24;
  pari_sp av = avma;
  GEN q24, N, n, m, run;

  if (!is_scalar_t(tx)) pari_err(typeer,"trueeta");
  x = upper_half(x, &prec);
  run = dbltor(1 - 1e-8);
  m = gen_1;
  N = gen_0;
  for(;;)
  {
    n = ground( real_i(x) );
    if (signe(n)) { x = gsub(x,n); N = addii(N, n); }
    if (gcmp(cxnorm(x), run) > 0) break;
    x = gdivsg(-1,x);
    m = gmul(m, gsqrt(mulcxmI(x), prec));
  }
  Nmod24 = umodiu(N, 24);
  if (Nmod24) m = gmul(m, e12(Nmod24, prec));
  q24 = gexp(gdivgs(gmul(Pi2n(1,prec), mulcxI(x)), 24),prec); /* e(x/24) */
  m = gmul(q24, m);
  if (24 * gexpo(q24) >= -bit_accuracy(prec))
    m = gmul(m, inteta( gpowgs(q24,24) ));
  return gerepileupto(av, m);
}

GEN
eta0(GEN x, long flag,long prec)
{
  return flag? trueeta(x,prec): eta(x,prec);
}

GEN
jell(GEN x, long prec)
{
  long tx = typ(x);
  pari_sp av = avma;
  GEN p1;

  if (!is_scalar_t(tx) || tx == t_PADIC)
  {
    GEN p2, q = qq(x,prec);
    p1 = gdiv(inteta(gsqr(q)), inteta(q));
    p1 = gmul2n(gsqr(p1),1);
    p1 = gmul(q,gpowgs(p1,12));
    p2 = gaddsg(768,gadd(gsqr(p1),gdivsg(4096,p1)));
    p1 = gmulsg(48,p1);
    return gerepileupto(av, gadd(p2,p1));
  }
  p1 = gdiv(trueeta(gmul2n(x,1),prec), trueeta(x,prec));
  p1 = gsqr(gsqr(gsqr(p1)));
  p1 = gadd(gmul2n(gsqr(p1),8), ginv(p1));
  return gerepileupto(av, gpowgs(p1,3));
}

GEN
weberf2(GEN x, long prec)
{
  pari_sp av=avma, tetpil;
  GEN p1,p2;

  p1=gsqrt(gen_2,prec);
  p2=gdiv(trueeta(gmul2n(x,1),prec),trueeta(x,prec));
  tetpil=avma;
  return gerepile(av,tetpil,gmul(p1,p2));
}

GEN
weberf1(GEN x, long prec)
{
  pari_sp av=avma, tetpil;
  GEN p1,p2;

  p1=trueeta(gmul2n(x,-1),prec); p2=trueeta(x,prec);
  tetpil=avma;
  return gerepile(av,tetpil,gdiv(p1,p2));
}

GEN
weberf(GEN x, long prec)
{
  pari_sp av = avma, tetpil;
  GEN p1, p2;

  p1 = gdiv(trueeta(gmul2n(gaddgs(x,1),-1),prec),trueeta(x,prec));
  p2 = exp_Ir(divrs(mppi(prec),-24));
  tetpil = avma;
  return gerepile(av,tetpil, gmul(p1,p2));
}

GEN
weber0(GEN x, long flag,long prec)
{
  switch(flag)
  {
    case 0: return weberf(x,prec);
    case 1: return weberf1(x,prec);
    case 2: return weberf2(x,prec);
    default: pari_err(flagerr,"weber");
  }
  return NULL; /* not reached */
}

GEN
theta(GEN q, GEN z, long prec)
{
  long l, n;
  pari_sp av = avma;
  GEN ps, qn, y, zy, ps2, k, zold;

  l = precision(q);
  n = precision(z); if (n && n < l) l = n;
  if (l) prec = l;
  z = gtofp(z, prec);
  q = gtofp(q, prec); if (gexpo(q) >= 0) pari_err(talker,"q >= 1 in theta");
  zold = NULL; /* gcc -Wall */
  zy = imag_i(z);
  if (gcmp0(zy)) k = gen_0;
  else
  {
    GEN lq = glog(q,prec); k = roundr(divrr(zy, real_i(lq)));
    if (signe(k)) { zold = z; z = gadd(z, mulcxmI(gmul(lq,k))); }
  }
  qn = gen_1;
  ps2 = gsqr(q);
  ps = gneg_i(ps2);
  y = gsin(z,prec);
  for (n = 1;; n++)
  {
    GEN t;
    qn = gmul(qn,ps);
    ps = gmul(ps,ps2);
    t = gmul(qn, gsin(gmulsg(2*n+1,z),prec)); y = gadd(y, t);
    if (gexpo(t) < -bit_accuracy(prec)) break;
  }
  if (signe(k))
  {
    y = gmul(y, gmul(powgi(q,sqri(k)),
                     gexp(gmul(mulcxI(zold),shifti(k,1)), prec)));
    if (mod2(k)) y = gneg_i(y);
  }
  return gerepileupto(av, gmul(y, gmul2n(gsqrt(gsqrt(q,prec),prec),1)));
}

GEN
thetanullk(GEN q, long k, long prec)
{
  long l, n;
  pari_sp av = avma;
  GEN p1, ps, qn, y, ps2;

  if (k < 0) pari_err(talker,"k < 0 in thetanullk");
  l = precision(q);
  if (l) prec = l;
  q = gtofp(q, prec); if (gexpo(q) >= 0) pari_err(talker,"q >= 1 in theta");

  if (!(k&1)) { avma = av; return gen_0; }
  qn = gen_1;
  ps2 = gsqr(q);
  ps = gneg_i(ps2);
  y = gen_1;
  for (n = 1;; n++)
  {
    GEN t;
    qn = gmul(qn,ps);
    ps = gmul(ps,ps2);
    t = gmul(qn, powuu(2*n+1, k)); y = gadd(y, t);
    if (gexpo(t) < -bit_accuracy(prec)) break;
  }
  p1 = gmul2n(gsqrt(gsqrt(q,prec),prec),1);
  if (k&2) y = gneg_i(y);
  return gerepileupto(av, gmul(p1, y));
}

/* [d^i theta/dz^i(q, 0), i = 1, 3, .., 2*k - 1] */
GEN
vecthetanullk(GEN q, long k, long prec)
{
  long i, l, n;
  pari_sp av = avma;
  GEN p1, ps, qn, y, ps2;

  l = precision(q); if (l) prec = l;
  q = gtofp(q, prec); if (gexpo(q) >= 0) pari_err(talker,"q >= 1 in theta");

  qn = gen_1;
  ps2 = gsqr(q);
  ps = gneg_i(ps2);
  y = cgetg(k+1, t_VEC); for (i = 1; i <= k; i++) gel(y,i) = gen_1;
  for (n = 1;; n++)
  {
    ulong N = 2*n + 1;
    GEN t = NULL/*-Wall*/, P = utoipos(N), N2 = sqru(N);
    qn = gmul(qn,ps);
    ps = gmul(ps,ps2);
    for (i = 1; i <= k; i++)
    {
      t = gmul(qn, P); gel(y,i) = gadd(gel(y,i), t);
      P = mulii(P, N2);
    }
    if (gexpo(t) < -bit_accuracy(prec)) break;
  }
  p1 = gmul2n(gsqrt(gsqrt(q,prec),prec),1);
  for (i = 2; i <= k; i+=2) gel(y,i) = gneg_i(gel(y,i));
  return gerepileupto(av, gmul(p1, y));
}
