/* $Id: buch1.c 12064 2010-01-09 18:30:42Z bill $

Copyright (C) 2000  The PARI group.

This file is part of the PARI/GP package.

PARI/GP is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation. It is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY WHATSOEVER.

Check the License for details. You should have received a copy of it, along
with the package; see the file 'COPYING'. If not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA. */
#include "pari.h"
#include "paripriv.h"

/*******************************************************************/
/*                                                                 */
/*       Hilbert and Ray Class field using CM (Schertz)            */
/*                                                                 */
/*******************************************************************/
static int
isoforder2(GEN form)
{
  GEN a = gel(form,1), b = gel(form,2), c = gel(form,3);
  return !signe(b) || absi_equal(a,b) || equalii(a,c);
}

GEN
getallforms(GEN D, long *pth, GEN *ptz)
{
  ulong d = itou(D), dover3 = d/3, t, b2, a, b, c, h;
  GEN z, L = cgetg((long)(sqrt(d) * log2(d)), t_VEC);
  b2 = b = (d&1); h = 0; z = gen_1;
  if (!b) /* b = 0 treated separately to avoid special cases */
  {
    t = d >> 2; /* (b^2 - D) / 4*/
    for (a=1; a*a<=t; a++)
      if (c = t/a, t == c*a)
      {
	z = mului(a,z);
        gel(L,++h) = mkvecsmall3(a,b,c);
      }
    b = 2; b2 = 4;
  }
  /* now b > 0 */
  for ( ; b2 <= dover3; b += 2, b2 = b*b)
  {
    t = (b2 + d) >> 2; /* (b^2 - D) / 4*/
    /* b = a */
    if (c = t/b, t == c*b)
    {
      z = mului(b,z);
      gel(L,++h) = mkvecsmall3(b,b,c);
    }
    /* b < a < c */
    for (a = b+1; a*a < t; a++)
      if (c = t/a, t == c*a)
      {
	z = mului(a,z);
        gel(L,++h) = mkvecsmall3(a, b,c);
	gel(L,++h) = mkvecsmall3(a,-b,c);
      }
    /* a = c */
    if (a * a == t) { z = mului(a,z); gel(L,++h) = mkvecsmall3(a,b,c); }
  }
  *pth = h; *ptz = z; setlg(L,h+1); return L;
}

static ulong
check_pq(GEN gp, GEN z, long d, GEN D)
{
  ulong p = itou(gp);
  if (!umodiu(z,p) || kross(d,(long)p) <= 0 || 
    gcmp1((GEN)redimag(primeform_u(D,p))[1]))
      pari_err(talker,"[quadhilbert] incorrect values in pq: %lu", p);
  return p;
}
#define MOD4(x) ((x)&3)
/* find P and Q two non principal prime ideals (above p,q) such that
 *   (pq, 2z) = 1  [coprime to all class group representatives]
 *   cl(P) = cl(Q) if P has order 2 in Cl(K)
 * Try to have e = 24 / gcd(24, (p-1)(q-1)) as small as possible */
static void
get_pq(GEN D, GEN z, GEN pq, ulong *ptp, ulong *ptq)
{
  const long MAXL = 50;
  GEN form, wp = cgetg(MAXL,t_VECSMALL), wlf = cgetg(MAXL,t_VEC);
  long i, ell, p, l = 1, d = itos(D);
  byteptr diffell = diffptr + 2;

  if (pq && typ(pq)==t_VEC)
  {
    if (lg(pq) != 3) pari_err(typeer, "quadhilbert (pq)");
    *ptp = check_pq(gel(pq,1),z,d,D);
    *ptq = check_pq(gel(pq,2),z,d,D); return;
  }

  ell = 3;
  while (l < MAXL)
  {
    NEXT_PRIME_VIADIFF_CHECK(ell, diffell);
    if (umodiu(z,ell) && kross(d,ell) > 0)
    {
      form = redimag(primeform_u(D,ell));
      if (gcmp1(gel(form,1))) continue;
      gel(wlf,l) = form;
      wp[l]  = ell; l++;
    }
  }
  setlg(wp,l); setlg(wlf,l);

  for (i=1; i<l; i++)
    if (wp[i] % 3 == 1) break;
  if (i==l) i = 1;
  p = wp[i]; form = gel(wlf,i);
  i = l;
  if (isoforder2(form))
  {
    long oki = 0;
    for (i=1; i<l; i++)
      if (gequal(gel(wlf,i),form))
      {
        if (MOD4(p) == 1 || MOD4(wp[i]) == 1) break;
        if (!oki) oki = i; /* not too bad, still e = 2 */
      }
    if (i==l) i = oki;
    if (!i) pari_err(bugparier,"quadhilbertimag (can't find p,q)");
  }
  else
  {
    if (MOD4(p) == 3)
    {
      for (i=1; i<l; i++)
        if (MOD4(wp[i]) == 1) break;
    }
    if (i==l) i = 1;
  }
  *ptp = p;
  *ptq = wp[i];
}

static GEN
gpq(GEN form, ulong p, ulong q, long e, GEN sqd, GEN u, long prec)
{
  long a = form[1], a2 = a << 1; /* gcd(a2, u) = 2 */
  GEN p1,p2,p3,p4;
  GEN w = lift(chinese(gmodulss(-form[2], a2), u));
  GEN al = mkcomplex(gdivgs(w, -a2), gdivgs(sqd, a2));
  p1 = trueeta(gdivgs(al,p), prec);
  p2 = p == q? p1: trueeta(gdivgs(al,q), prec);
  p3 = trueeta(gdiv(al,muluu(p,q)), prec);
  p4 = trueeta(al, prec);
  return gpowgs(gdiv(gmul(p1,p2),gmul(p3,p4)), e);
}

/* returns an equation for the Hilbert class field of Q(sqrt(D)), D < 0 */
static GEN
quadhilbertimag(GEN D, GEN pq)
{
  GEN z, L, P, qfp, u;
  pari_sp av = avma;
  long h, i, e, prec;
  ulong p, q;

  if (DEBUGLEVEL>1) (void)timer2();
  if (cmpiu(D,11) <= 0) return pol_x[0];
  L = getallforms(D,&h,&z);
  if (DEBUGLEVEL>1) msgtimer("class number = %ld",h);
  if (h == 1) { avma=av; return pol_x[0]; }

  get_pq(D, z, pq, &p, &q);
  e = 24 / cgcd((p%24 - 1) * (q%24 - 1), 24);
  if(DEBUGLEVEL>1) fprintferr("p = %lu, q = %lu, e = %ld\n",p,q,e);
  qfp = primeform_u(D, p);
  if (p == q)
  {
    u = (GEN)compimagraw(qfp, qfp)[2];
    u = gmodulo(u, shifti(sqru(p),1));
  }
  else
  {
    GEN qfq = primeform_u(D, q);
    GEN up = mkintmodu(itou(gel(qfp,2)), p << 1);
    GEN uq = mkintmodu(itou(gel(qfq,2)), q << 1);
    u = chinese(up,uq);
  }
  /* u modulo 2pq */
  prec = 3;
  for(;;)
  {
    long ex, exmax = 0;
    pari_sp av0 = avma;
    GEN lead, sqd = sqrtr_abs(itor(D, prec));
    P = cgetg(h+1,t_VEC);
    for (i=1; i<=h; i++)
    {
      GEN s = gpq(gel(L,i), p, q, e, sqd, u, prec);
      if (DEBUGLEVEL>3) fprintferr("%ld ", i);
      gel(P,i) = s; ex = gexpo(s); if (ex > 0) exmax += ex;
    }
    if (DEBUGLEVEL>1) msgtimer("roots");
    /* to avoid integer overflow (1 + 0.) */
    lead = (exmax < bit_accuracy(prec))? gen_1: real_1(prec);

    P = real_i( roots_to_pol_intern(lead,P,0,0) );
    P = grndtoi(P,&exmax);
    if (DEBUGLEVEL>1) msgtimer("product, error bits = %ld",exmax);
    if (exmax <= -10)
    {
      if (pq && degpol(srgcd(P, derivpol(P)))) { avma = av; return gen_0; }
      break;
    }
    avma = av0; prec += (DEFAULTPREC-2) + (1 + (exmax >> TWOPOTBITS_IN_LONG));
    if (DEBUGLEVEL) pari_warn(warnprec,"quadhilbertimag",prec);
  }
  return gerepileupto(av,P);
}

GEN
quadhilbert(GEN D, GEN flag, long prec)
{
  if (typ(D) != t_INT)
  {
    D = checkbnf(D);
    if (degpol(gmael(D,7,1)) != 2)
      pari_err(talker,"not a polynomial of degree 2 in quadhilbert");
    D = gmael(D,7,3);
  }
  else if (!isfundamental(D))
    pari_err(talker,"quadhilbert needs a fundamental discriminant");
  return (signe(D)>0)? quadhilbertreal(D,prec)
                     : quadhilbertimag(D,flag);
}

#define to_approx(nf,a) ((GEN)gmul(gmael((nf),5,1), (a))[1])
/* Z-basis for a (over C) */
static GEN
get_om(GEN nf, GEN a) {
  return mkvec2(to_approx(nf,gel(a,2)),
                to_approx(nf,gel(a,1)));
}

/* Compute all elts in class group G = [|G|,c,g], c=cyclic factors, g=gens.
 * Set list[j + 1] = g1^e1...gk^ek where j is the integer
 *   ek + ck [ e(k-1) + c(k-1) [... + c2 [e1]]...] */
static GEN
getallelts(GEN bnr)
{
  GEN nf,G,C,c,g, *list, **pows, *gk;
  long lc,i,j,k,no;

  nf = checknf(bnr);
  G = gel(bnr,5);

  no = itos(gel(G,1));
  c = gel(G,2);
  g = gel(G,3); lc = lg(c)-1;
  list = (GEN*) cgetg(no+1,t_VEC);
  if (!lc)
  {
    list[1] = idealhermite(nf,gen_1);
    return (GEN)list;
  }
  pows = (GEN**)cgetg(lc+1,t_VEC);
  c = shallowcopy(c); settyp(c, t_VECSMALL);
  for (i=1; i<=lc; i++)
  {
    c[i] = k = itos(gel(c,i));
    gk = (GEN*)cgetg(k, t_VEC); gk[1] = gel(g,i);
    for (j=2; j<k; j++)
      gk[j] = idealmodidele(bnr, idealmul(nf, gk[j-1], gk[1]));
    pows[i] = gk; /* powers of g[i] */
  }

  C = cgetg(lc+1, t_VECSMALL); C[1] = c[lc];
  for (i=2; i<=lc; i++) C[i] = C[i-1] * c[lc-i+1];
  /* C[i] = c(k-i+1) * ... * ck */
  /* j < C[i+1] <==> j only involves g(k-i)...gk */
  i = 1; list[1] = 0; /* dummy */
  for (j=1; j < C[1]; j++)
    list[j + 1] = pows[lc][j];
  for (   ; j<no; j++)
  {
    GEN p1,p2;
    if (j == C[i+1]) i++;
    p2 = pows[lc-i][j/C[i]];
    p1 = list[j%C[i] + 1];
    if (p1) p2 = idealmodidele(bnr, idealmul(nf,p2,p1));
    list[j + 1] = p2;
  }
  list[1] = idealhermite(nf,gen_1);
  return (GEN)list;
}

/* x quadratic integer (approximate), recognize it. If error return NULL */
static GEN
findbezk(GEN nf, GEN x)
{
  GEN a,b, M = gmael(nf,5,1), u = gcoeff(M,1,2);
  long ea, eb;

  b = grndtoi(gdiv(imag_i(x), imag_i(u)), &eb);
  a = grndtoi(real_i(gsub(x, gmul(b,u))), &ea);
  if (ea > -20 || eb > -20) return NULL;
  if (!signe(b)) return a;
  return coltoalg(nf, mkcol2(a,b));
}

static GEN
findbezk_pol(GEN nf, GEN x)
{
  long i, lx = lg(x);
  GEN y = cgetg(lx,t_POL);
  for (i=2; i<lx; i++)
    if (! (gel(y,i) = findbezk(nf,gel(x,i))) ) return NULL;
  y[1] = x[1]; return y;
}

/* allow t_QF[IR], and t_VEC/t_COL with 3 components */
GEN
form_to_ideal(GEN x)
{
  long tx = typ(x);
  GEN b;
  if ((is_vec_t(tx) || lg(x) != 4)
       && tx != t_QFR && tx != t_QFI) pari_err(typeer,"form_to_ideal");
  b = negi(gel(x,2)); if (mpodd(b)) b = addis(b,1);
  return mkmat2( mkcol2(gel(x,1), gen_0),
                 mkcol2((GEN)shifti(b,-1), gen_1) );
}

/* P approximation computed at initial precision prec. Compute needed prec
 * to know P with 1 word worth of trailing decimals */
static long
get_prec(GEN P, long prec)
{
  long k = gprecision(P);
  if (k == 3) return (prec<<1)-2; /* approximation not trustworthy */
  k = prec - k; /* lost precision when computing P */
  if (k < 0) k = 0;
  k += MEDDEFAULTPREC + (gexpo(P) >> TWOPOTBITS_IN_LONG);
  if (k <= prec) k = (prec<<1)-2; /* dubious: old prec should have worked */
  return k;
}

/* Compute data for ellphist */
static GEN
ellphistinit(GEN om, long prec)
{
  GEN res,om1b,om2b, om1 = gel(om,1), om2 = gel(om,2);

  if (gsigne(imag_i(gdiv(om1,om2))) < 0) { swap(om1,om2); om = mkvec2(om1,om2); }
  om1b = gconj(om1);
  om2b = gconj(om2); res = cgetg(4,t_VEC);
  gel(res,1) = gdivgs(elleisnum(om,2,0,prec),12);
  gel(res,2) = gdiv(PiI2(prec), gmul(om2, imag_i(gmul(om1b,om2))));
  gel(res,3) = om2b; return res;
}

/* Computes log(phi^*(z,om)), using res computed by ellphistinit */
static GEN
ellphist(GEN om, GEN res, GEN z, long prec)
{
  GEN u = imag_i(gmul(z, gel(res,3)));
  GEN zst = gsub(gmul(u, gel(res,2)), gmul(z,gel(res,1)));
  return gsub(ellsigma(om,z,1,prec),gmul2n(gmul(z,zst),-1));
}

/* Computes phi^*(la,om)/phi^*(1,om) where (1,om) is an oriented basis of the
   ideal gf*gc^{-1} */
static GEN
computeth2(GEN om, GEN la, long prec)
{
  GEN p1,p2,res = ellphistinit(om,prec);

  p1 = gsub(ellphist(om,res,la,prec), ellphist(om,res,gen_1,prec));
  p2 = imag_i(p1);
  if (gexpo(real_i(p1))>20 || gexpo(p2)> bit_accuracy(min(prec,lg(p2)))-10)
    return NULL;
  return gexp(p1,prec);
}

/* Computes P_2(X)=polynomial in Z_K[X] closest to prod_gc(X-th2(gc)) where
   the product is over the ray class group bnr.*/
static GEN
computeP2(GEN bnr, GEN la, long prec)
{
  long clrayno, i, first = 1;
  pari_sp av=avma, av2;
  GEN listray,P0,P,f,lanum, nf = checknf(bnr);

  f = gmael3(bnr,2,1,1);
  la = algtobasis_i(nf,la);
  listray = getallelts(bnr);
  clrayno = lg(listray)-1; av2 = avma;
PRECPB:
  if (!first)
  {
    if (DEBUGLEVEL) pari_warn(warnprec,"computeP2",prec);
    nf = gerepileupto(av2, nfnewprec(checknf(bnr),prec));
  }
  first = 0; lanum = to_approx(nf,la);
  P = cgetg(clrayno+1,t_VEC);
  for (i=1; i<=clrayno; i++)
  {
    GEN om = get_om(nf, idealdiv(nf,f,gel(listray,i)));
    GEN s = computeth2(om,lanum,prec);
    if (!s) { prec = (prec<<1)-2; goto PRECPB; }
    gel(P,i) = s;
  }
  P0 = roots_to_pol(P, 0);
  P = findbezk_pol(nf, P0);
  if (!P) { prec = get_prec(P0, prec); goto PRECPB; }
  return gerepilecopy(av, P);
}

#define nexta(a) (a>0 ? -a : 1-a)
static GEN
do_compo(GEN x, GEN y)
{
  long a, i, l = lg(y);
  GEN z;
  y = shallowcopy(y); /* y := t^deg(y) y(#/t) */
  for (i = 2; i < l; i++)
    if (signe(y[i])) gel(y,i) = monomial(gel(y,i), l-i-1, MAXVARN);
  for  (a = 0;; a = nexta(a))
  {
    if (a) x = gsubst(x, 0, gaddsg(a, pol_x[0]));
    z = gsubst(subres(x,y), MAXVARN, pol_x[0]);
    if (issquarefree(z)) return z;
  }
}
#undef nexta

static GEN
galoisapplypol(GEN nf, GEN s, GEN x)
{
  long i, lx = lg(x);
  GEN y = cgetg(lx,t_POL);

  for (i=2; i<lx; i++) gel(y,i) = galoisapply(nf,s,gel(x,i));
  y[1] = x[1]; return y;
}

/* x quadratic, write it as ua + v, u,v rational */
static GEN
findquad(GEN a, GEN x, GEN p)
{
  long tu, tv;
  pari_sp av = avma;
  GEN u,v;
  if (typ(x) == t_POLMOD) x = gel(x,2);
  if (typ(a) == t_POLMOD) a = gel(a,2);
  u = poldivrem(x, a, &v);
  u = simplify(u); tu = typ(u);
  v = simplify(v); tv = typ(v);
  if (!is_scalar_t(tu) || !is_scalar_t(tv))
    pari_err(talker, "incorrect data in findquad");
  x = v;
  if (!gcmp0(u)) x = gadd(gmul(u, pol_x[varn(a)]), x);
  if (typ(x) == t_POL) x = gmodulo(x,p);
  return gerepileupto(av, x);
}

static GEN
findquad_pol(GEN p, GEN a, GEN x)
{
  long i, lx = lg(x);
  GEN y = cgetg(lx,t_POL);
  for (i=2; i<lx; i++) gel(y,i) = findquad(a, gel(x,i), p);
  y[1] = x[1]; return y;
}

static GEN
compocyclo(GEN nf, long m, long d)
{
  GEN sb,a,b,s,p1,p2,p3,res,polL,polLK,nfL, D = gel(nf,3);
  long ell,vx;

  p1 = quadhilbertimag(D, gen_0);
  p2 = cyclo(m,0);
  if (d==1) return do_compo(p1,p2);

  ell = m&1 ? m : (m>>2);
  if (equalui(ell,D)) /* ell = |D| */
  {
    p2 = gcoeff(nffactor(nf,p2),1,1);
    return do_compo(p1,p2);
  }
  if (ell%4 == 3) ell = -ell;
  /* nf = K = Q(a), L = K(b) quadratic extension = Q(t) */
  polLK = quadpoly(stoi(ell)); /* relative polynomial */
  res = rnfequation2(nf, polLK);
  vx = varn(nf[1]);
  polL = gsubst(gel(res,1),0,pol_x[vx]); /* = charpoly(t) */
  a = gsubst(lift(gel(res,2)), 0,pol_x[vx]);
  b = gsub(pol_x[vx], gmul(gel(res,3), a));
  nfL = initalg(polL, DEFAULTPREC);
  p1 = gcoeff(nffactor(nfL,p1),1,1);
  p2 = gcoeff(nffactor(nfL,p2),1,1);
  p3 = do_compo(p1,p2); /* relative equation over L */
  /* compute non trivial s in Gal(L / K) */
  sb= gneg(gadd(b, truecoeff(polLK,1))); /* s(b) = Tr(b) - b */
  s = gadd(pol_x[vx], gsub(sb, b)); /* s(t) = t + s(b) - b */
  p3 = gmul(p3, galoisapplypol(nfL, s, p3));
  return findquad_pol(gel(nf,1), a, p3);
}

/* I integral ideal in HNF. (x) = I, x small in Z ? */
static long
isZ(GEN I)
{
  GEN x = gcoeff(I,1,1);
  if (signe(gcoeff(I,1,2)) || !equalii(x, gcoeff(I,2,2))) return 0;
  return is_bigint(x)? -1: itos(x);
}

/* Treat special cases directly. return NULL if not special case */
static GEN
treatspecialsigma(GEN nf, GEN gf)
{
  GEN p1, p2, tryf, D = gel(nf,3);
  long Ds, fl, i = isZ(gf);

  if (i == 1) return quadhilbertimag(gel(nf,3), NULL); /* f = 1 ? */

  if (equaliu(D,3))
  {
    if (i == 4 || i == 5 || i == 7) return cyclo(i,0);
    if (!equaliu(gcoeff(gf,1,1),9) || !equaliu(content(gf),3)) return NULL;
    p1 = (GEN)nfroots(nf,cyclo(3,0))[1]; /* f = P_3^3 */
    return gadd(monomial(gen_1,3,0), p1); /* x^3+j */
  }
  if (equaliu(D,4))
  {
    if (i == 3 || i == 5) return cyclo(i,0);
    if (i != 4) return NULL;
    p1 = (GEN)nfroots(nf,cyclo(4,0))[1];
    return gadd(monomial(gen_1,2,0), p1); /* x^2+i */
  }
  Ds = smodis(D,48);
  if (i)
  {
    if (i==2 && Ds%16== 8) return compocyclo(nf, 4,1);
    if (i==3 && Ds% 3== 1) return compocyclo(nf, 3,1);
    if (i==4 && Ds% 8== 1) return compocyclo(nf, 4,1);
    if (i==6 && Ds   ==40) return compocyclo(nf,12,1);
    return NULL;
  }

  p1 = gcoeff(gf,1,1); /* integer > 0 */
  p2 = gcoeff(gf,2,2);
  if (gcmp1(p2)) { fl = 0; tryf = p1; }
  else {
    if (Ds % 16 != 8 || !equaliu(Q_content(gf),2)) return NULL;
    fl = 1; tryf = shifti(p1,-1);
  }
  /* tryf integer > 0 */
  if (cmpiu(tryf, 3) <= 0 || signe(remii(D, tryf)) || !isprime(tryf))
    return NULL;

  i = itos(tryf); if (fl) i *= 4;
  return compocyclo(nf,i,2);
}

/* return a vector of all roots of 1 in bnf [not necessarily quadratic] */
static GEN
getallrootsof1(GEN bnf)
{
  GEN T, u, nf = checknf(bnf), tu;
  long i, n = itos(gmael3(bnf,8,4,1));

  if (n == 2) {
    long N = degpol(gel(nf,1));
    return mkvec2(gscalcol_i(gen_m1, N),
                  gscalcol_i(gen_1, N));
  }
  tu = poltobasis(nf, gmael3(bnf,8,4,2));
  T = eltmul_get_table(nf, tu);
  u = cgetg(n+1, t_VEC); gel(u,1) = tu;
  for (i=2; i <= n; i++) gel(u,i) = gmul(T, gel(u,i-1));
  return u;
}

static GEN
get_lambda(GEN bnr)
{
  GEN allf, bnf, nf, pol, f, la, P, labas, lamodf, u;
  long a, b, f2, i, lu, v;

  allf = conductor(bnr,NULL,2);
  bnr = gel(allf,2);
  f = gmael(allf,1,1);
  bnf= gel(bnr,1);
  nf = gel(bnf,7);
  pol= gel(nf,1); v = varn(pol);
  P = treatspecialsigma(nf,f);
  if (P) return P;

  f2 = 2 * itos(gcoeff(f,1,1));
  u = getallrootsof1(bnf); lu = lg(u);
  for (i=1; i<lu; i++)
    gel(u,i) = colreducemodHNF(gel(u,i), f, NULL); /* roots of 1, mod f */
  if (DEBUGLEVEL>1)
    fprintferr("quadray: looking for [a,b] != unit mod 2f\n[a,b] = ");
  for (a=0; a<f2; a++)
    for (b=0; b<f2; b++)
    {
      la = deg1pol_i(stoi(a), stoi(b), v); /* ax + b */
      if (umodiu(gnorm(mkpolmod(la, pol)), f2) != 1) continue;
      if (DEBUGLEVEL>1) fprintferr("[%ld,%ld] ",a,b);

      labas = poltobasis(nf, la);
      lamodf = colreducemodHNF(labas, f, NULL);
      for (i=1; i<lu; i++)
        if (gequal(lamodf, gel(u,i))) break;
      if (i < lu) continue; /* la = unit mod f */
      if (DEBUGLEVEL)
      {
        if (DEBUGLEVEL>1) fprintferr("\n");
        fprintferr("lambda = %Z\n",la);
      }
      return labas;
    }
  pari_err(bugparier,"get_lambda");
  return NULL;
}

GEN
quadray(GEN D, GEN f, GEN flag, long prec)
{
  GEN bnr, y, pol, bnf;
  pari_sp av = avma;

  if (flag)
  {
    if (typ(flag) != t_VEC || lg(flag)!=3) pari_err(flagerr,"quadray");
  }
  if (typ(D) != t_INT)
  {
    bnf = checkbnf(D);
    if (degpol(gmael(bnf,7,1)) != 2)
      pari_err(talker,"not a polynomial of degree 2 in quadray");
    D = gmael(bnf,7,3);
  }
  else
  {
    if (!isfundamental(D))
      pari_err(talker,"quadray needs a fundamental discriminant");
    pol = quadpoly0(D, fetch_user_var("y"));
    bnf = bnfinit0(pol, signe(D)>0?1:0, NULL, prec);
  }
  bnr = bnrinit0(bnf,f,1);
  if (gcmp1(gmael(bnr,5,1))) { avma = av; return pol_x[0]; }
  if (signe(D) > 0)
    y = bnrstark(bnr,NULL,prec);
  else
  {
    bnr = gel(conductor(bnr,NULL,2), 2);
    if (!flag) flag = get_lambda(bnr);
    if (typ(flag) == t_POL) y = flag; /* special case */
    else
      y = computeP2(bnr,flag,prec);
  }
  return gerepileupto(av, y);
}

/*******************************************************************/
/*                                                                 */
/*         CLASS GROUP AND REGULATOR (McCURLEY, BUCHMANN)          */
/*                   QUADRATIC FIELDS                              */
/*                                                                 */
/*******************************************************************/
/* For largeprime() hashtable. Note that hashed pseudoprimes are odd (unless
 * 2 | index), hence the low order bit is not useful. So we hash
 * HASHBITS bits starting at bit 1, not bit 0 */
#define HASHBITS 10
static const long HASHT = 1 << HASHBITS;

static long
hash(long q) { return (q & ((1 << (HASHBITS+1)) - 1)) >> 1; }
#undef HASHBITS

/* See buch2.c:
 * subFB contains split p such that \prod p > sqrt(Disc)
 * powsubFB contains powers of forms in subFB */
#define RANDOM_BITS 4
static const long CBUCH = (1<<RANDOM_BITS)-1;

static ulong limhash;
static long KC, KC2, PRECREG;
static long *primfact, *exprimfact, *FB, *numFB, **hashtab;
static GEN powsubFB, vperm, subFB, Disc, sqrtD, isqrtD, badprim;

/*******************************************************************/
/*                                                                 */
/*  Routines related to binary quadratic forms (for internal use)  */
/*                                                                 */
/*******************************************************************/
/* output canonical representative wrt projection Cl^+ --> Cl (a > 0) */
static GEN
qfr3_canon(GEN x)
{
  GEN a = gel(x,1), c = gel(x,3);
  if (signe(a) < 0) {
    if (absi_equal(a,c)) return qfr3_rho(x,Disc,isqrtD);
    setsigne(a, 1);
    setsigne(c,-1);
  }
  return x;
}
static GEN
qfr5_canon(GEN x)
{
  GEN a = gel(x,1), c = gel(x,3);
  if (signe(a) < 0) {
    if (absi_equal(a,c)) return qfr5_rho(x,Disc,sqrtD,isqrtD);
    setsigne(a, 1);
    setsigne(c,-1);
  }
  return x;
}
static GEN
QFR5_comp(GEN x,GEN y) { return qfr5_canon(qfr5_comp(x,y,Disc,sqrtD,isqrtD)); }
static GEN
QFR3_comp(GEN x, GEN y) { return qfr3_canon(qfr3_comp(x,y,Disc,isqrtD)); }

/* compute rho^n(x) */
static GEN
qrf5_rho_pow(GEN x, long n)
{
  long i;
  pari_sp av = avma, lim = stack_lim(av, 1);
  for (i=1; i<=n; i++)
  {
    x = qfr5_rho(x,Disc,sqrtD,isqrtD);
    if (low_stack(lim, stack_lim(av,1)))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"qrf5_rho_pow");
      x = gerepilecopy(av, x);
    }
  }
  return gerepilecopy(av, x);
}

static GEN
qfr5_pf(GEN D, long p)
{
  GEN y = primeform_u(D,p);
  return qfr5_canon(qfr5_red(qfr_to_qfr5(y,PRECREG), Disc, sqrtD, isqrtD));
}

static GEN
qfr3_pf(GEN D, long p)
{
  GEN y = primeform_u(D,p);
  return qfr3_canon(qfr3_red(y, Disc, isqrtD));
}

#define qfi_pf primeform_u

/* Warning: ex[0] not set in general */
static GEN
init_form(long *ex, GEN (*comp)(GEN,GEN))
{
  long i, l = lg(powsubFB);
  GEN F = NULL;
  for (i=1; i<l; i++)
    if (ex[i])
    {
      GEN t = gmael(powsubFB,i,ex[i]);
      F = F? comp(F,t): t;
    }
  return F;
}
static GEN
qfr5_factorback(long *ex) { return init_form(ex, &QFR5_comp); }
static GEN
qfi_factorback(long *ex) { return init_form(ex, &compimag); }

static GEN
random_form(GEN ex, GEN (*comp)(GEN,GEN))
{
  long i, l = lg(ex);
  pari_sp av = avma;
  GEN F;
  for(;;)
  {
    for (i=1; i<l; i++) ex[i] = random_bits(RANDOM_BITS);
    if ((F = init_form(ex, comp))) return F;
    avma = av;
  }
}
static GEN
qfr3_random(GEN ex){ return random_form(ex, &QFR3_comp); }
static GEN
qfi_random(GEN ex) { return random_form(ex, &compimag); }

/*******************************************************************/
/*                                                                 */
/*                     Common subroutines                          */
/*                                                                 */
/*******************************************************************/
double
check_bach(double cbach, double B)
{
  if (cbach >= B)
   pari_err(talker,"sorry, couldn't deal with this field. PLEASE REPORT");
  cbach *= 2; if (cbach > B) cbach = B;
  if (DEBUGLEVEL) fprintferr("\n*** Bach constant: %f\n", cbach);
  return cbach;
}

#if 0
static long
factorquad(GEN f, long kcz, ulong limp)
{
  ulong p;
  long i, k, lo;
  pari_sp av;
  GEN x = gel(f,1);

  if (is_pm1(x)) { primfact[0] = 0; return 1; }
  av = avma; lo = 0;
  x = absi(x);
  for (i=1; ; i++)
  {
    int stop;
    k = Z_lvalrem_stop(x, (ulong)FB[i], &stop);
    if (k) { primfact[++lo] = i; exprimfact[lo] = k; }
    if (stop) break;
    if (i == kcz) { avma = av; return 0; }
  }
  avma = av;
  if (lgefint(x) != 3 || (p=(ulong)x[2]) > limhash) return 0;

  if (p != 1 && p <= limp)
  {
    if (badprim && cgcd(p, umodiu(badprim,p)) > 1) return 0;
    primfact[++lo] = numFB[p]; exprimfact[lo] = 1;
    p = 1;
  }
  primfact[0] = lo; return p;
}

#else /* Same, Z_lvalrem_stop unrolled. Ugly but more than 30% faster :-( */

/* Is |q| <= p ? */
static int
isless_iu(GEN q, ulong p) {
  long l = lgefint(q);
  return l==2 || (l == 3 && (ulong)q[2] <= p);
}

static long
factorquad(GEN f, long kcz, ulong limp)
{
  ulong X;
  long i, lo;
  pari_sp av;
  GEN x = gel(f,1);

  if (is_pm1(x)) { primfact[0] = 0; return 1; }
  av = avma; lo = 0;
  for (i=1; lgefint(x) > 3; i++)
  {
    ulong p = (ulong)FB[i], r;
    GEN q = diviu_rem(x, p, &r);
    if (!r)
    {
      long k = 0;
      do { k++; x = q; q = diviu_rem(x, p, &r); } while (!r);
      primfact[++lo] = i; exprimfact[lo] = k; 
    }
    if (isless_iu(q,p)) {
      avma = av;
      if (lgefint(x) == 3) { X = (ulong)x[2]; goto END; }
      return 0;
    }
    if (i == kcz) { avma = av; return 0; }
  }
  avma = av; X = (ulong)x[2];
  for (;; i++)
  { /* single precision affair, split for efficiency */
    ulong p = (ulong)FB[i];
    ulong q = X / p, r = X % p; /* gcc makes a single div */
    if (!r)
    {
      long k = 0;
      do { k++; X = q; q = X / p; r = X % p; } while (!r);
      primfact[++lo] = i; exprimfact[lo] = k; 
    }
    if (q <= p) break;
    if (i == kcz) return 0;
  }
END:
  if (X > limhash) return 0;
  if (X != 1 && X <= limp)
  {
    if (badprim && cgcd(X, umodiu(badprim,X)) > 1) return 0;
    primfact[++lo] = numFB[X]; exprimfact[lo] = 1;
    X = 1;
  }
  primfact[0] = lo; return X;
}
#endif

/* Check for a "large prime relation" involving q; q may not be prime */
static long *
largeprime(long q, long *ex, long np, long nrho)
{
  const long hashv = hash(q);
  long *pt, i, l = lg(subFB);

  for (pt = hashtab[hashv]; ; pt = (long*) pt[0])
  {
    if (!pt)
    {
      pt = (long*) gpmalloc((l+3) * sizeof(long));
      *pt++ = nrho; /* nrho = pt[-3] */
      *pt++ = np;   /* np   = pt[-2] */
      *pt++ = q;    /* q    = pt[-1] */
      pt[0] = (long)hashtab[hashv];
      for (i=1; i<l; i++) pt[i]=ex[i];
      hashtab[hashv]=pt; return NULL;
    }
    if (pt[-1] == q) break;
  }
  for(i=1; i<l; i++)
    if (pt[i] != ex[i]) return pt;
  return (pt[-2]==np)? (GEN)NULL: pt;
}

static void
clearhash(long **hash)
{
  long *pt;
  long i;
  for (i=0; i<HASHT; i++) {
    for (pt = hash[i]; pt; ) {
      void *z = (void*)(pt - 3);
      pt = (long*) pt[0]; free(z);
    }
    hash[i] = NULL;
  }
}

/* p | conductor of order of disc D ? */
static int
is_bad(GEN D, ulong p)
{
  pari_sp av = avma;
  int r;
  if (p == 2)
  {
    r = mod16(D) >> 1;
    if (r && signe(D) < 0) r = 8-r;
    return (r < 4);
  }
  r = (remii(D, muluu(p,p)) == gen_0); /* p^2 | D ? */
  avma = av; return r;
}

/* create FB, numFB; set badprim. Return L(kro_D, 1) */
static GEN
FBquad(GEN Disc, long n2, long n)
{
  GEN Res = real_1(DEFAULTPREC);
  long i, p, s, LIM;
  pari_sp av;
  byteptr d = diffptr;

  numFB = cgetg(n2+1, t_VECSMALL);
  FB    = cgetg(n2+1, t_VECSMALL);
  av = avma;
  KC = 0; i = 0;
  maxprime_check((ulong)n2);
  badprim = gen_1;
  for (p = 0;;) /* p <= n2 */
  {
    NEXT_PRIME_VIADIFF(p, d);
    if (!KC && p > n) KC = i;
    if (p > n2) break;
    s = krois(Disc,p);
    Res = mulur(p, divrs(Res, p - s));
    switch (s)
    {
      case -1: break; /* inert */
      case  0: /* ramified */
        if (is_bad(Disc, (ulong)p)) { badprim = muliu(badprim, p); break; }
        /* fall through */
      default:  /* split */
        i++; numFB[p] = i; FB[i] = p; break;
    }
  }
  if (!KC) return NULL;
  KC2 = i;
  setlg(FB, KC2+1);
  if (DEBUGLEVEL)
  {
    msgtimer("factor base");
    if (DEBUGLEVEL>7) fprintferr("FB = %Z\n", FB);
  }
  LIM = (expi(Disc) < 16)? 100: 1000;
  while (p < LIM)
  {
    s = krois(Disc,p);
    Res = mulur(p, divrs(Res, p - s));
    NEXT_PRIME_VIADIFF(p, d);
  }
  if (badprim != gen_1)
    gerepileall(av, 2, &Res, &badprim);
  else
  {
    badprim = NULL;
    Res = gerepileuptoleaf(av, Res);
  }
  return Res;
}

/* create vperm, return subFB */
static GEN
subFBquad(GEN D, double PROD, long KC)
{
  long i, j, minSFB, lgsub = 1, ino = 1, lv = KC+1;
  double prod = 1.;
  pari_sp av;
  GEN no;

  minSFB = (expi(D) > 15)? 3: 2;
  vperm = cgetg(lv, t_VECSMALL);
  av = avma;
  no    = cgetg(lv, t_VECSMALL);
  for (j = 1; j < lv; j++)
  {
    ulong p = FB[j];
    if (!umodiu(D, p)) no[ino++] = j; /* ramified */
    else
    {
      vperm[lgsub++] = j;
      prod *= p;
      if (lgsub > minSFB && prod > PROD) break;
    }
  }
  if (j == lv) return NULL;
  i = lgsub;
  for (j = 1; j < ino;i++,j++) vperm[i] = no[j];
  for (     ; i < lv; i++)     vperm[i] = i;
  if (DEBUGLEVEL) msgtimer("subFBquad (%ld elt.)", lgsub-1);
  no = gclone(vecslice(vperm, 1, lgsub-1));
  avma = av; return no;
}

/* assume n >= 1, x[i][j] = subFB[i]^j, for j = 1..n */
static GEN
powsubFBquad(long n)
{
  pari_sp av = avma;
  long i,j, l = lg(subFB);
  GEN F, y, x = cgetg(l, t_VEC);

  if (PRECREG) /* real */
  {
    for (i=1; i<l; i++)
    {
      F = qfr5_pf(Disc, FB[subFB[i]]);
      y = cgetg(n+1, t_VEC); gel(x,i) = y;
      gel(y,1) = F;
      for (j=2; j<=n; j++) gel(y,j) = QFR5_comp(gel(y,j-1), F);
    }
  }
  else /* imaginary */
  {
    for (i=1; i<l; i++)
    {
      F = qfi_pf(Disc, FB[subFB[i]]);
      y = cgetg(n+1, t_VEC); gel(x,i) = y;
      gel(y,1) = F;
      for (j=2; j<=n; j++) gel(y,j) = compimag(gel(y,j-1), F);
    }
  }
  if (DEBUGLEVEL) msgtimer("powsubFBquad");
  x = gclone(x); avma = av; return x;
}

static void
sub_fact(GEN col, GEN F)
{
  GEN b = gel(F,2);
  long i;
  for (i=1; i<=primfact[0]; i++)
  {
    ulong k = primfact[i], p = FB[k];
    long e = exprimfact[i];
    if (umodiu(b, p<<1) > p) e = -e;
    col[k] -= e;
  }
}
static void
add_fact(GEN col, GEN F)
{
  GEN b = gel(F,2);
  long i;
  for (i=1; i<=primfact[0]; i++)
  {
    ulong k = primfact[i], p = FB[k];
    long e = exprimfact[i];
    if (umodiu(b, p<<1) > p) e = -e;
    col[k] += e;
  }
}

static GEN
get_clgp(GEN Disc, GEN W, GEN *ptD, long prec)
{
  GEN res, *init, u1, D = smithrel(W,NULL,&u1), Z = prec? real_0(prec): NULL;
  long i, j, l = lg(W), c = lg(D);

  if (DEBUGLEVEL) msgtimer("smith/class group");
  res=cgetg(c,t_VEC); init = (GEN*)cgetg(l,t_VEC);
  for (i=1; i<l; i++) init[i] = primeform_u(Disc, FB[vperm[i]]);
  for (j=1; j<c; j++)
  {
    GEN g = NULL;
    if (prec)
    {
      for (i=1; i<l; i++)
      {
        GEN t, u = gcoeff(u1,i,j);
        if (!signe(u)) continue;
        t = qfr3_pow(init[i], u, Disc, isqrtD);
        g = g? qfr3_comp(g, t, Disc, isqrtD): t;
      }
      g = qfr3_to_qfr(qfr3_canon(qfr3_red(g, Disc, isqrtD)), Z);
    }
    else
    {
      for (i=1; i<l; i++)
      {
        GEN t, u = gcoeff(u1,i,j);
        if (!signe(u)) continue;
        t = powgi(init[i], u);
        g = g? compimag(g, t): t;
      }
    }
    gel(res,j) = g;
  }
  if (DEBUGLEVEL) msgtimer("generators");
  *ptD = D; return res;
}

static long
trivial_relations(GEN mat, long KC, GEN C, GEN Disc)
{
  long i, j = 0;
  GEN col;
  for (i = 1; i <= KC; i++) 
  { /* ramified prime ==> trivial relation */
    if (umodiu(Disc, FB[i])) continue;
    col = const_vecsmall(KC, 0);
    col[i] = 2; j++;
    gel(mat,j) = col;
    gel(C,j) = gen_0;
  }
  return j;
}

static void
dbg_all(char *phase, long s, long n)
{
  fprintferr("\nTime %s rel [#rel/#test = %ld/%ld]: %ld\n", phase,s,n,timer2());
}

void
wr_rel(GEN col)
{
  long i, l = lg(col);
  fprintferr("\nrel = ");
  for (i=1; i<l; i++)
    if (col[i]) fprintferr("%ld^%ld ",i,col[i]);
  fprintferr("\n");
}

void
dbg_rel(long s, GEN col)
{
  if (DEBUGLEVEL == 1) fprintferr("%ld ",s);
  else { fprintferr("cglob = %ld. ", s); wr_rel(col); }
  flusherr(); 
}
/* Imaginary Quadratic fields */

static void
imag_relations(long need, long *pc, long lim, ulong LIMC, GEN mat)
{
  long lgsub = lg(subFB), current = *pc, nbtest = 0, s = 0;
  long i, fpc;
  pari_sp av;
  GEN col, form, ex = cgetg(lgsub, t_VECSMALL);

  if (!current) current = 1;
  av = avma;
  for(;;)
  {
    if (s >= need) break;
    avma = av;
    form = qfi_random(ex);
    form = compimag(form, qfi_pf(Disc, FB[current]));
    nbtest++; fpc = factorquad(form,KC,LIMC);
    if (!fpc)
    {
      if (DEBUGLEVEL>1) fprintferr(".");
      continue;
    }
    if (fpc > 1)
    {
      long *fpd = largeprime(fpc,ex,current,0);
      ulong b1, b2, p;
      GEN form2;
      if (!fpd)
      {
        if (DEBUGLEVEL>1) fprintferr(".");
        continue;
      }
      form2 = compimag(qfi_factorback(fpd), qfi_pf(Disc, FB[fpd[-2]]));
      p = fpc << 1;
      b1 = umodiu(gel(form2,2), p);
      b2 = umodiu(gel(form,2),  p);
      if (b1 != b2 && b1+b2 != p) continue;

      col = gel(mat,++s);
      add_fact(col, form);
      (void)factorquad(form2,KC,LIMC);
      if (b1==b2)
      {
        for (i=1; i<lgsub; i++) col[subFB[i]] += fpd[i]-ex[i];
        sub_fact(col, form2); col[fpd[-2]]++;
      }
      else
      {
        for (i=1; i<lgsub; i++) col[subFB[i]] += -fpd[i]-ex[i];
        add_fact(col, form2); col[fpd[-2]]--;
      }
    }
    else
    {
      col = gel(mat,++s);
      for (i=1; i<lgsub; i++) col[subFB[i]] = -ex[i];
      add_fact(col, form);
    }
    col[current]--;
    if (++current > KC) current = 1;
  }
  if (DEBUGLEVEL) dbg_all("random", s, nbtest);
  *pc = current;
}

static int
imag_be_honest()
{
  long p, fpc, s = KC, nbtest = 0;
  GEN F, ex = cgetg(lg(subFB), t_VECSMALL);
  pari_sp av = avma;

  while (s<KC2)
  {
    p = FB[s+1]; if (DEBUGLEVEL) fprintferr(" %ld",p);
    F = compimag(qfi_pf(Disc, p), qfi_random(ex));
    fpc = factorquad(F,s,p-1);
    if (fpc == 1) { nbtest=0; s++; }
    else
      if (++nbtest > 40) return 0;
    avma = av;
  }
  return 1;
}

/* Real Quadratic fields */

static void
real_relations(long need, long *pc, long lim, ulong LIMC, GEN mat, GEN C)
{
  long lgsub = lg(subFB), current = *pc, nbtest = 0, s = 0;
  long i, fpc, endcycle, rhoacc, rho;
  /* in a 2nd phase, don't include FB[current] but run along the cyle 
   * ==> get more units */
  int first = (current == 0);
  pari_sp av, av1, limstack;
  GEN d, col, form, form0, form1, ex = cgetg(lgsub, t_VECSMALL);

  if (!current) current = 1;
  if (lim > need) lim = need;
  av = avma; limstack = stack_lim(av,1);
  for(;;)
  {
    if (s >= need) break;
    if (first && s >= lim) {
      first = 0;
      if (DEBUGLEVEL) dbg_all("initial", s, nbtest);
    }
    avma = av; form = qfr3_random(ex);
    if (!first) form = QFR3_comp(form, qfr3_pf(Disc, FB[current]));
    av1 = avma;
    form0 = form; form1 = NULL;
    endcycle = rhoacc = 0;
    rho = -1;

CYCLE:
    if (endcycle || rho > 5000) continue;
    if (low_stack(limstack, stack_lim(av,1)))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"real_relations");
      gerepileall(av1, form1? 2: 1, &form, &form1);
    }
    if (rho < 0) rho = 0; /* first time in */
    else
    {
      form = qfr3_rho(form, Disc, isqrtD); rho++;
      rhoacc++;
      if (first)
        endcycle = (absi_equal(gel(form,1),gel(form0,1))
             && equalii(gel(form,2),gel(form0,2)));
      else
      {
        if (absi_equal(gel(form,1), gel(form,3))) /* a = -c */
        {
          if (absi_equal(gel(form,1),gel(form0,1)) &&
                  equalii(gel(form,2),gel(form0,2))) continue;
          form = qfr3_rho(form, Disc, isqrtD); rho++;
        }
        else
          { setsigne(form[1],1); setsigne(form[3],-1); }
        if (equalii(gel(form,1),gel(form0,1)) &&
            equalii(gel(form,2),gel(form0,2))) continue;
      }
    }
    nbtest++; fpc = factorquad(form,KC,LIMC);
    if (!fpc)
    {
      if (DEBUGLEVEL>1) fprintferr(".");
      goto CYCLE;
    }
    if (fpc > 1)
    { /* look for Large Prime relation */
      long *fpd = largeprime(fpc,ex,first? 0: current,rhoacc);
      ulong b1, b2, p;
      GEN form2;
      if (!fpd)
      {
        if (DEBUGLEVEL>1) fprintferr(".");
        goto CYCLE;
      }
      if (!form1)
      {
        form1 = qfr5_factorback(ex);
        if (!first) form1 = QFR5_comp(form1, qfr5_pf(Disc, FB[current]));
      }
      form1 = qrf5_rho_pow(form1, rho);
      rho = 0;

      form2 = qfr5_factorback(fpd);
      if (fpd[-2]) form2 = QFR5_comp(form2, qfr5_pf(Disc, FB[fpd[-2]]));
      form2 = qrf5_rho_pow(form2, fpd[-3]);
      if (!absi_equal(gel(form2,1),gel(form2,3)))
      {
        setsigne(form2[1], 1);
        setsigne(form2[3],-1);
      }
      p = fpc << 1;
      b1 = umodiu(gel(form2,2), p);
      b2 = umodiu(gel(form1,2), p);
      if (b1 != b2 && b1+b2 != p) goto CYCLE;

      col = gel(mat,++s);
      add_fact(col, form1);
      (void)factorquad(form2,KC,LIMC);
      if (b1==b2)
      {
        for (i=1; i<lgsub; i++) col[subFB[i]] += fpd[i]-ex[i];
        sub_fact(col, form2);
        if (fpd[-2]) col[fpd[-2]]++;
        d = qfr5_dist(subii(gel(form1,4),gel(form2,4)),
                      divrr(gel(form1,5),gel(form2,5)), PRECREG);
      }
      else
      {
        for (i=1; i<lgsub; i++) col[subFB[i]] += -fpd[i]-ex[i];
        add_fact(col, form2);
        if (fpd[-2]) col[fpd[-2]]--;
        d = qfr5_dist(addii(gel(form1,4),gel(form2,4)),
                      mulrr(gel(form1,5),gel(form2,5)), PRECREG);
      }
      if (DEBUGLEVEL) fprintferr(" %ldP",s);
    }
    else
    { /* standard relation */
      if (!form1)
      {
        form1 = qfr5_factorback(ex);
        if (!first) form1 = QFR5_comp(form1, qfr5_pf(Disc, FB[current]));
      }
      form1 = qrf5_rho_pow(form1,rho);
      rho = 0;

      col = gel(mat,++s);
      for (i=1; i<lgsub; i++) col[subFB[i]] = -ex[i];
      add_fact(col, form1);
      d = qfr5_dist(gel(form1,4), gel(form1,5), PRECREG);
      if (DEBUGLEVEL) fprintferr(" %ld",s);
    }
    affrr(d, gel(C,s));
    if (first)
    {
      if (s >= lim) continue;
      goto CYCLE;
    }
    else
    {
      col[current]--;
      if (++current > KC) current = 1;
    }
  }
  if (DEBUGLEVEL) dbg_all("random", s, nbtest);
  *pc = current;
}

static int
real_be_honest()
{
  long p, fpc, s = KC, nbtest = 0;
  GEN F,F0, ex = cgetg(lg(subFB), t_VECSMALL);
  pari_sp av = avma;

  while (s<KC2)
  {
    p = FB[s+1]; if (DEBUGLEVEL) fprintferr(" %ld",p);
    F = QFR3_comp(qfr3_random(ex), qfr3_pf(Disc, p));
    for (F0 = F;;)
    {
      fpc = factorquad(F,s,p-1);
      if (fpc == 1) { nbtest=0; s++; break; }
      if (++nbtest > 40) return 0;
      F = qfr3_canon(qfr3_rho(F, Disc, isqrtD));
      if (equalii(gel(F,1),gel(F0,1))
       && equalii(gel(F,2),gel(F0,2))) break;
    }
    avma = av;
  }
  return 1;
}

static GEN
gcdreal(GEN a,GEN b)
{
  long e;
  GEN k1,r;

  if (!signe(a)) return mpabs(b);
  if (!signe(b)) return mpabs(a);

  if (typ(a)==t_INT)
  {
    if (typ(b)==t_INT) return gcdii(a,b);
    a = itor(a, lg(b));
  }
  else if (typ(b)==t_INT)
  {
    b = itor(b, lg(a));
  }
  if (expo(a)<-5) return absr(b);
  if (expo(b)<-5) return absr(a);
  a=absr(a); b=absr(b);
  while (expo(b) >= -5  && signe(b))
  {
    k1 = gcvtoi(divrr(a,b),&e);
    if (e > 0) return NULL;
    r=subrr(a,mulir(k1,b)); a=b; b=r;
  }
  return absr(a);
}

static int
get_R(GEN C, long sreg, GEN z, GEN *ptR)
{
  GEN R = gen_1;
  double c;
  long i;

  if (PRECREG)
  {
    R = mpabs(gel(C,1));
    for (i=2; i<=sreg; i++)
    {
      R = gcdreal(gel(C,i), R);
      if (!R) return fupb_PRECI;
    }
    if (gexpo(R) <= -3)
    {
      if (DEBUGLEVEL) fprintferr("regulator is zero.\n");
      return fupb_RELAT;
    }
    if (DEBUGLEVEL) fprintferr("#### Tentative regulator: %Z\n",R);
  }
  c = gtodouble(gmul(z, R));
  if (c < 0.8 || c > 1.3) return fupb_RELAT;
  *ptR = R; return fupb_NONE;
}

static int
quad_be_honest()
{
  int r;
  if (KC2 <= KC) return 1;
  if (DEBUGLEVEL)
    fprintferr("be honest for primes from %ld to %ld\n", FB[KC+1],FB[KC2]);
  r = PRECREG? real_be_honest(): imag_be_honest();
  if (DEBUGLEVEL) { fprintferr("\n"); msgtimer("be honest"); }
  return r;
}

GEN
buchquad(GEN D, double cbach, double cbach2, long RELSUP, long prec)
{
  pari_sp av0 = avma, av, av2;
  long i, s, current, triv, nrelsup, nreldep, need, nsubFB;
  ulong LIMC, LIMC2, cp;
  GEN h, W, cyc, res, gen, dep, mat, C, extraC, B, R, resc, Res, z;
  double drc, lim, LOGD, LOGD2;

  check_quaddisc(D, &s, /*junk*/&i, "buchquad");
  Disc = D;
  if (s < 0)
  {
    if (cmpiu(Disc,4) <= 0)
    {
      GEN z = cgetg(5,t_VEC);
      gel(z,1) = gel(z,4) = gen_1; gel(z,2) = gel(z,3) = cgetg(1,t_VEC);
      return z;
    }
    PRECREG = 0;
  } else {
    PRECREG = max(prec+1, MEDDEFAULTPREC + 2*(expi(Disc)>>TWOPOTBITS_IN_LONG));
  }
  if (DEBUGLEVEL) (void)timer2();
  primfact   = new_chunk(100);
  exprimfact = new_chunk(100);
  hashtab = (long**) new_chunk(HASHT);
  for (i=0; i<HASHT; i++) hashtab[i] = NULL;

  drc = fabs(gtodouble(Disc));
  LOGD = log(drc);
  LOGD2 = LOGD * LOGD;

  lim = sqrt(drc);
  /* resc = sqrt(D) w / 2^r1 (2pi)^r2 ~ hR / L(chi,1) */
  if (PRECREG) resc = dbltor(lim / 2.);
  else         resc = dbltor(lim / PI);
  if (!PRECREG) lim /= sqrt(3.);
  cp = (ulong)exp(sqrt(LOGD*log(LOGD)/8.0));
  if (cp < 20) cp = 20;
  if (cbach > 6.) {
    if (cbach2 < cbach) cbach2 = cbach;
    cbach = 6.;
  }
  if (cbach <= 0.) pari_err(talker,"Bach constant <= 0 in buchquad");
  av = avma; cbach /= 2;
  powsubFB = subFB = NULL;

/* LIMC = Max(cbach*(log D)^2, exp(sqrt(log D loglog D) / 8)) */
START: avma = av; cbach = check_bach(cbach,6.);
  if (subFB) gunclone(subFB);
  if (powsubFB) gunclone(powsubFB);
  clearhash(hashtab);
  nreldep = nrelsup = 0;
  LIMC = (ulong)(cbach*LOGD2);
  if (LIMC < cp) { LIMC = cp; cbach = (double)LIMC / LOGD2; }
  LIMC2 = (ulong)(max(cbach,cbach2)*LOGD2);
  if (LIMC2 < LIMC) LIMC2 = LIMC;
  if (PRECREG)
  {
    sqrtD  = sqrtr(itor(Disc,PRECREG));
    isqrtD = truncr(sqrtD);
  }

  Res = FBquad(Disc, LIMC2, LIMC);
  if (!Res) goto START;
  subFB = subFBquad(Disc, lim + 0.5, KC);
  if (!subFB) goto START;
  nsubFB = lg(subFB) - 1;
  powsubFB = powsubFBquad(CBUCH+1);
  limhash = (LIMC & HIGHMASK)? (HIGHBIT>>1): LIMC*LIMC;

  need = KC + RELSUP - 2;
  current = 0;
  W = NULL;
  s = nsubFB + RELSUP;
  av2 = avma;

MORE:
  if ((nreldep & 3) == 1 || (nrelsup & 7) == 1) {
    if (DEBUGLEVEL) fprintferr("*** Changing sub factor base\n");
    gunclone(subFB);
    gunclone(powsubFB);
    subFB = gclone(vecslice(vperm, 1, nsubFB));
    powsubFB = powsubFBquad(CBUCH+1);
    clearhash(hashtab);
  }
  need += 2;
  mat    = cgetg(need+1, t_MAT);
  extraC = cgetg(need+1, t_VEC);
  if (!W) { /* first time */
    C = extraC;
    triv = trivial_relations(mat, KC, C, Disc);
    if (DEBUGLEVEL) fprintferr("KC = %ld, need %ld relations\n", KC, need);
  } else {
    triv = 0;
    if (DEBUGLEVEL) fprintferr("...need %ld more relations\n", need);
  }
  if (PRECREG) {
    for (i = triv+1; i<=need; i++) {
      gel(mat,i) = const_vecsmall(KC, 0);
      gel(extraC,i) = cgetr(PRECREG);
    }
    real_relations(need - triv, &current, s,LIMC,mat + triv,extraC + triv);
  } else {
    for (i = triv+1; i<=need; i++) {
      gel(mat,i) = const_vecsmall(KC, 0);
      gel(extraC,i) = gen_0;
    }
    imag_relations(need - triv, &current, s,LIMC,mat + triv);
  }

  if (!W)
    W = hnfspec_i((long**)mat,vperm,&dep,&B,&C,nsubFB);
  else
    W = hnfadd_i(W,vperm,&dep,&B,&C, mat,extraC);
  gerepileall(av2, 4, &W,&C,&B,&dep);
  need = lg(dep)>1? lg(dep[1])-1: lg(B[1])-1;
  if (need)
  {
    if (++nreldep > 15 && cbach < 1) goto START;
    goto MORE;
  }

  h = dethnf_i(W);
  if (DEBUGLEVEL) fprintferr("\n#### Tentative class number: %Z\n", h);

  z = mulrr(Res, resc); /* ~ hR if enough relations, a multiple otherwise */
  switch(get_R(C, (lg(C)-1) - (lg(B)-1) - (lg(W)-1), divir(h,z), &R))
  {
    case fupb_PRECI:
      PRECREG = (PRECREG<<1)-2;
      cbach /= 2; goto START;

    case fupb_RELAT:
      if (++nrelsup <= 7 || cbach > 1) {
        need = min(KC, nrelsup); 
        if (cbach > 1 && nsubFB < 3 && lg(vperm) > 3) nsubFB++;
        goto MORE;
      }
      goto START;
  }
  /* DONE */
  if (!quad_be_honest()) goto START;
  clearhash(hashtab);

  gen = get_clgp(Disc,W,&cyc,PRECREG);
  gunclone(subFB);
  gunclone(powsubFB);
  res = cgetg(5,t_VEC);
  gel(res,1) = h;
  gel(res,2) = cyc;
  gel(res,3) = gen;
  gel(res,4) = R; return gerepilecopy(av0,res);
}

GEN
buchimag(GEN D, GEN c, GEN c2, GEN REL)
{ return buchquad(D,gtodouble(c),gtodouble(c2),itos(REL), 0); }

GEN
buchreal(GEN D, GEN flag, GEN c, GEN c2, GEN REL, long prec) {
  if (signe(flag)) pari_err(impl,"narrow class group");
  return buchquad(D,gtodouble(c),gtodouble(c2),itos(REL), prec);
}

GEN
quadclassunit0(GEN x, long flag, GEN data, long prec)
{
  long lx, RELSUP;
  double cbach, cbach2;

  if (!data) lx=1;
  else
  {
    lx = lg(data);
    if (typ(data)!=t_VEC || lx > 7)
      pari_err(talker,"incorrect parameters in quadclassunit");
    if (lx > 4) lx = 4;
  }
  cbach = cbach2 = 0.2; /* was 0.1, but slower on average for 20 digits disc */
  RELSUP = 5;
  switch(lx)
  {
    case 4: RELSUP = itos(gel(data,3));
    case 3: cbach2 = gtodouble(gel(data,2));
    case 2: cbach  = gtodouble(gel(data,1));
  }
  if (flag) pari_err(impl,"narrow class group");
  return buchquad(x,cbach,cbach2,RELSUP,prec);
}
