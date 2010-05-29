/* $Id: galconj.c 7857 2006-04-11 17:28:55Z kb $

Copyright (C) 2000-2003  The PARI group.

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
/*************************************************************************/
/**									**/
/**                           GALOIS CONJUGATES        		        **/
/**									**/
/*************************************************************************/

GEN
galoisconj(GEN nf)
{
  GEN     x, y, z;
  long i, lz, v;
  pari_sp av = avma;
  nf = checknf(nf);
  x = gel(nf,1);
  v = varn(x);
  if (v == 0)
    nf = gsubst(nf, 0, pol_x[MAXVARN]);
  else
  {
    x = shallowcopy(x);
    setvarn(x, 0);
  }
  z = nfroots(nf, x);
  lz = lg(z);
  y = cgetg(lz, t_COL);
  for (i = 1; i < lz; i++)
  {
    GEN     p1 = lift(gel(z,i));
    setvarn(p1, v);
    gel(y,i) = p1;
  }
  return gerepileupto(av, y);
}

/* nbmax: maximum number of possible conjugates */
GEN
galoisconj2pol(GEN x, long nbmax, long prec)
{
  long i, n, v, nbauto;
  pari_sp av = avma;
  GEN     y, w, polr, p1, p2;
  n = degpol(x);
  if (n <= 0)
    return cgetg(1, t_VEC);
  if (gisirreducible(x) == gen_0)
    pari_err(redpoler, "galoisconj2pol");
  polr = roots(x, prec);
  p1 = gel(polr,1);
  nbauto = 1;
  prec = (long)bit_accuracy_mul(prec, L2SL10 * 0.75);
  w = cgetg(n + 2, t_VEC);
  gel(w,1) = gen_1;
  for (i = 2; i <= n; i++)
    gel(w,i) = gmul(p1, gel(w,i - 1));
  v = varn(x);
  y = cgetg(nbmax + 1, t_COL);
  gel(y,1) = pol_x[v];
  for (i = 2; i <= n && nbauto < nbmax; i++)
  {
    w[n + 1] = polr[i];
    p1 = lindep2(w, prec);
    if (signe(p1[n + 1]))
    {
      setlg(p1, n + 1);
      p2 = gdiv(gtopolyrev(p1, v), negi(gel(p1,n + 1)));
      if (gdvd(poleval(x, p2), x))
      {
	gel(y,++nbauto) = p2;
	if (DEBUGLEVEL > 1)
	  fprintferr("conjugate %ld: %Z\n", i, y[nbauto]);
      }
    }
  }
  setlg(y, 1 + nbauto);
  return gerepileupto(av, gen_sort(y, 0, cmp_pol));
}

GEN
galoisconj2(GEN nf, long nbmax, long prec)
{
  long i, j, n, r1, ru, nbauto;
  pari_sp av = avma;
  GEN     x, y, w, polr, p1, p2;
  if (typ(nf) == t_POL)
    return galoisconj2pol(nf, nbmax, prec);
  nf = checknf(nf);
  x = gel(nf,1);
  n = degpol(x);
  if (n <= 0)
    return cgetg(1, t_VEC);
  r1 = nf_get_r1(nf);
  p1 = gel(nf,6);
  prec = precision(gel(p1,1));
  /* accuracy in decimal digits */
  prec = (long)bit_accuracy_mul(prec, L2SL10 * 0.75);
  ru = (n + r1) >> 1;
  nbauto = 1;
  polr = cgetg(n + 1, t_VEC);
  for (i = 1; i <= r1; i++)
    polr[i] = p1[i];
  for (j = i; i <= ru; i++)
  {
    GEN z = gel(p1,i);
    gel(polr,j++) = z;
    gel(polr,j++) = gconj(z);
  }
  p2 = gmael(nf, 5, 1);
  w = cgetg(n + 2, t_VEC);
  for (i = 1; i <= n; i++)
    w[i] = coeff(p2, 1, i);
  y = cgetg(nbmax + 1, t_COL);
  gel(y,1) = pol_x[varn(x)];
  for (i = 2; i <= n && nbauto < nbmax; i++)
  {
    w[n + 1] = polr[i];
    p1 = lindep2(w, prec);
    if (signe(p1[n + 1]))
    {
      setlg(p1, n + 1);
      settyp(p1, t_COL);
      p2 = gdiv(coltoliftalg(nf, p1), negi(gel(p1,n + 1)));
      if (gdvd(poleval(x, p2), x))
      {
	gel(y,++nbauto) = p2;
	if (DEBUGLEVEL > 1)
	  fprintferr("conjugate %ld: %Z\n", i, y[nbauto]);
      }
    }
  }
  setlg(y, 1 + nbauto);
  return gerepileupto(av, gen_sort(y, 0, cmp_pol));
}
/*************************************************************************/
/**									**/
/**                           GALOISCONJ4             		        **/
/**									**/
/**                                                                     **/
/*************************************************************************/
/*DEBUGLEVEL:
  1: timing
  2: outline
  4: complete outline
  6: detail
  7: memory
  9: complete detail
*/

GEN
vandermondeinverseprep(GEN L)
{
  long i, j, n = lg(L);
  GEN V;
  V = cgetg(n, t_VEC);
  for (i = 1; i < n; i++)
  {
    pari_sp ltop=avma;
    GEN W=cgetg(n,t_VEC);
    for (j = 1; j < n; j++)
      if (i==j)
	gel(W,j) = gen_1;
      else
	gel(W,j) = gsub(gel(L,i),gel(L,j));
    gel(V,i) = gerepileupto(ltop,divide_conquer_prod(W,&gmul));
  }
  return V;
}

/* Calcule l'inverse de la matrice de van der Monde de T multiplie par den */
GEN
vandermondeinverse(GEN L, GEN T, GEN den, GEN prep)
{
  pari_sp ltop = avma;
  long i, n = lg(L)-1;
  GEN M, P;
  if (!prep)
    prep = vandermondeinverseprep(L);
  M = cgetg(n+1, t_MAT);
  for (i = 1; i <= n; i++)
  {
    P = gdiv(RgX_div_by_X_x(T, gel(L,i), NULL), gel(prep,i));
    gel(M,i) = RgX_to_RgV(P,n);
  }
  return gerepileupto(ltop, gmul(den, M));
}

/* Calcule les bornes sur les coefficients a chercher */
struct galois_borne
{
  GEN     l;
  long    valsol;
  long    valabs;
  GEN     bornesol;
  GEN     ladicsol;
  GEN     ladicabs;
  GEN     lbornesol;
};


GEN
initgaloisborne(GEN T, GEN dn, long prec, GEN *ptL, GEN *ptprep, GEN *ptdis)
{
  long i, n = degpol(T);
  GEN L, z, prep, den;
  pari_timer ti;

  if (DEBUGLEVEL>=4) (void)TIMER(&ti);
  L = roots(T, prec);
  if (DEBUGLEVEL>=4) msgTIMER(&ti,"roots");
  for (i = 1; i <= n; i++)
  {
    z = gel(L,i);
    if (signe(z[2])) break;
    L[i] = z[1];
  }
  if (DEBUGLEVEL>=4) (void)TIMER(&ti);
  prep = vandermondeinverseprep(L);
  if (!dn)
  {
    GEN dis, res = divide_conquer_prod(gabs(prep,prec), mpmul);
    disable_dbg(0);
    dis = ZX_disc_all(T, 1+logint(res,gen_2,NULL));
    disable_dbg(-1);
    den = indexpartial(T,dis);
    if (ptdis) *ptdis = dis;
  }
  else
  {
    if (typ(dn) != t_INT || signe(dn) <= 0)
      pari_err(talker, "incorrect denominator in initgaloisborne: %Z", dn);
    den = dn;
  }
  if (ptprep) *ptprep = prep;
  *ptL = L; return den;
}

/* ||| M ||| with respect to || x ||_oo. Assume M square t_MAT */
GEN
matrixnorm(GEN M, long prec)
{
  long i,j, n = lg(M);
  GEN B = real_0(prec);

  for (i = 1; i < n; i++)
  {
    GEN z = gabs(gcoeff(M,i,1), prec);
    for (j = 2; j < n; j++)
      z = gadd(z, gabs(gcoeff(M,i,j), prec));
    if (gcmp(z, B) > 0) B = z;
  }
  return B;
}

/* L a t_VEC/t_COL, return ||L||_oo */
GEN
supnorm(GEN L, long prec)
{
  long i, n = lg(L);
  GEN z, B;

  if (n == 1) return real_0(prec);
  B = gabs(gel(L,1), prec);
  for (i = 2; i < n; i++)
  {
    z = gabs(gel(L,i), prec);
    if (gcmp(z, B) > 0) B = z;
  }
  return B;
}

static GEN
galoisborne(GEN T, GEN dn, struct galois_borne *gb)
{
  pari_sp ltop = avma, av2;
  GEN borne, borneroots, borneabs;
  long n, prec;
  GEN L, M, prep, den;
  pari_timer ti;

  prec = ZX_get_prec(T);
  den = initgaloisborne(T,dn,prec, &L,&prep,NULL);
  if (!dn) den = gclone(den);
  if (DEBUGLEVEL>=4) TIMERstart(&ti);
  M = vandermondeinverse(L, gmul(T, real_1(prec)), den, prep);
  if (DEBUGLEVEL>=4) msgTIMER(&ti,"vandermondeinverse");
  borne = matrixnorm(M, prec);
  borneroots = supnorm(L, prec);
  n = degpol(T);
  borneabs = addsr(1, gmulsg(n, gpowgs(borneroots, n)));
  borneroots = addsr(1, gmul(borne, borneroots));
  av2 = avma;
  /*We use d-1 test, so we must overlift to 2^BITS_IN_LONG*/
  gb->valsol = logint(gmul2n(borneroots,2+BITS_IN_LONG), gb->l,NULL);
  gb->valabs = logint(gmul2n(borneabs,2), gb->l,NULL);
  gb->valabs = max(gb->valsol, gb->valabs);
  if (DEBUGLEVEL >= 4)
    fprintferr("GaloisConj:val1=%ld val2=%ld\n", gb->valsol, gb->valabs);
  avma = av2;
  gb->bornesol = gerepileupto(ltop, ceil_safe(mulrs(borneroots,2)));
  if (DEBUGLEVEL >= 9)
    fprintferr("GaloisConj: Bound %Z\n",borneroots);
  gb->ladicsol = powiu(gb->l, gb->valsol);
  gb->ladicabs = powiu(gb->l, gb->valabs);
  gb->lbornesol = subii(gb->ladicsol,gb->bornesol);
  if (!dn) { dn = icopy(den); gunclone(den); }
  return dn;
}

struct galois_lift
{
  GEN     T;
  GEN     den;
  GEN     p;
  GEN     L;
  GEN     Lden;
  long    e;
  GEN     Q;
  GEN     TQ;
  struct galois_borne *gb;
};

static GEN 
makeLden(GEN L,GEN den, struct galois_borne *gb)
{
  pari_sp ltop=avma;
  long i,l=lg(L);
  GEN Lden=cgetg(l,t_VEC);
  for (i=1;i<l;i++)
    gel(Lden,i) = mulii(gel(L,i),den);
  for (i=1;i<l;i++)
    gel(Lden,i) = modii(gel(Lden,i),gb->ladicsol);
  return gerepileupto(ltop,Lden);
}

/* Initialize the structure galois_lift */

static void
initlift(GEN T, GEN den, GEN p, GEN L, GEN Lden, struct galois_borne *gb, struct galois_lift *gl)
{
  pari_sp ltop, lbot;
  gl->gb=gb;
  gl->T = T;
  gl->den = den;
  gl->p = p;
  gl->L = L;
  gl->Lden = Lden;
  ltop = avma;
  gl->e = logint(gmul2n(gb->bornesol, 2+BITS_IN_LONG),p,NULL);
  gl->e = max(2,gl->e);
  lbot = avma;
  gl->Q = powiu(p, gl->e);
  gl->Q = gerepile(ltop, lbot, gl->Q);
  gl->TQ = FpX_red(T,gl->Q);
}

/*
 * Verifie que f est une solution presque surement et calcule sa permutation
 */
static int
poltopermtest(GEN f, struct galois_lift *gl, GEN pf)
{
  pari_sp ltop;
  GEN     fx, fp;
  long     i, j,ll;
  for (i = 2; i< lg(f); i++)
    if (cmpii(gel(f,i),gl->gb->bornesol)>0 
	&& cmpii(gel(f,i),gl->gb->lbornesol)<0)
    {
      if (DEBUGLEVEL>=4)
	fprintferr("GaloisConj: Solution too large, discard it.\n");
      if (DEBUGLEVEL>=8)
	fprintferr("f=%Z\n borne=%Z\n l-borne=%Z\n",f,gl->gb->bornesol,gl->gb->lbornesol);
      return 0;
    }
  ll=lg(gl->L);
  fp = cgetg(ll, t_VECSMALL);
  ltop = avma;
  for (i = 1; i < ll; i++)
    fp[i] = 1;
  for (i = 1; i < ll; i++)
  {
    fx = FpX_eval(f, gel(gl->L,i), gl->gb->ladicsol);
    for (j = 1; j < ll; j++)
    {
      if (fp[j] && equalii(fx, gel(gl->Lden,j)))
      {
	pf[i] = j;
	fp[j] = 0;
	break;
      }
    }
    if (j == ll)
      return 0;
    avma = ltop;
  }
  return 1;
}

/*
 * Soit P one polynome de \ZZ[X] , p one nombre premier , S\in\FF_p[X]/(Q) tel
 * que P(S)=0 [p,Q] Relever S en S_0 tel que P(S_0)=0 [p^e,Q]
 * Unclean stack.
 */
static long
monoratlift(GEN S, GEN q, GEN qm1old,struct galois_lift *gl, GEN frob)
{
  GEN tlift = polratlift(S,q,qm1old,qm1old,gl->den);
  if (tlift)
  {
    pari_sp ltop = avma;
    if(DEBUGLEVEL>=4)
      fprintferr("MonomorphismLift: trying early solution %Z\n",tlift);
    /*Rationals coefficients*/
    tlift = FpX_red(Q_muli_to_int(tlift, gl->den), gl->gb->ladicsol);
    if (poltopermtest(tlift, gl, frob))
    {
      if(DEBUGLEVEL>=4) fprintferr("MonomorphismLift: true early solution.\n");
      avma = ltop; return 1;
    }
    avma = ltop; 
    if(DEBUGLEVEL>=4) fprintferr("MonomorphismLift: false early solution.\n");
  }
  return 0;
}

static GEN
monomorphismratlift(GEN P, GEN S, struct galois_lift *gl, GEN frob)
{
  pari_sp ltop, lbot;
  long rt;
  GEN     Q=gl->T, p=gl->p;
  long    e=gl->e, level=1;
  GEN     q, qold, qm1, qm1old;
  GEN     W, Pr, Qr, Sr, Wr = gen_0, Qrold, Spow;
  long    i,nb,mask;
  GEN    *gptr[2];
  if (DEBUGLEVEL == 1) (void)timer2();
  rt = brent_kung_optpow(degpol(Q),1);
  q = p; qm1 = gen_1; /*during the run, we have p*qm1=q*/
  nb=hensel_lift_accel(e, &mask);
  Pr = FpX_red(P,q);
  Qr = (P==Q)?Pr:FpX_red(Q, q);/*A little speed up for automorphismlift*/
  W=FpX_FpXQ_compo(ZX_deriv(Pr),S,Qr,q);
  W=FpXQ_inv(W,Qr,q);
  qold = p; qm1old=gen_1;
  Qrold = Qr;
  gptr[0] = &S;
  gptr[1] = &Wr;
  for (i=0; i<nb;i++)
  {
    if (DEBUGLEVEL>=2)
    {
      level=(level<<1)-((mask>>i)&1);
      (void)timer2();
    }
    qm1 = (mask&(1<<i))?sqri(qm1):mulii(qm1, q);
    q   =  mulii(qm1, p);
    Pr = FpX_red(P, q);
    Qr = (P==Q)?Pr:FpX_red(Q, q);/*A little speed up for automorphismlift*/
    ltop = avma;
    Sr = S;
    Spow = FpXQ_powers(Sr, rt, Qr, q);

    if (i)
    {
      W = FpXQ_mul(Wr, FpX_FpXQV_compo(ZX_deriv(Pr),FpXV_red(Spow,qold),Qrold,qold), Qrold, qold);
      W = FpX_neg(W, qold);
      W = FpX_Fp_add(W, gen_2, qold);
      W = FpXQ_mul(Wr, W, Qrold, qold);
    }
    Wr = W;
    S = FpXQ_mul(Wr, FpX_FpXQV_compo(Pr, Spow, Qr, q),Qr,q);
    S = ZX_sub(Sr, S);
    lbot = avma;
    Wr = gcopy(Wr);
    S = FpX_red(S, q);
    gerepilemanysp(ltop, lbot, gptr, 2);
    if (i && i<nb-1 && frob && monoratlift(S,q,qm1old,gl,frob))
      return NULL;
    qold = q; qm1old=qm1; Qrold = Qr;
    if (DEBUGLEVEL >= 2)
      msgtimer("MonomorphismLift: lift to prec %d",level);
  }
  if (DEBUGLEVEL == 1)
    msgtimer("monomorphismlift()");
  return S;
}
/*
 * Let T be a polynomial in \ZZ[X] , p a prime number, S\in\FF_p[X]/(T) so
 * that T(S)=0 [p,T] Lift S in S_0 so that T(S_0)=0 [T,p^e]
 * Unclean stack.
 */
static GEN
automorphismlift(GEN S, struct galois_lift *gl, GEN frob)
{
  return  monomorphismratlift(gl->T, S, gl, frob);
}

GEN
monomorphismlift(GEN P, GEN S, GEN Q, GEN p, long e)
{
  struct galois_lift gl;
  gl.T=Q;
  gl.p=p;
  gl.e=e;
  return monomorphismratlift(P,S,&gl,NULL);
}

struct galois_testlift
{
  long    n;
  long    f;
  long    g;
  GEN     bezoutcoeff;
  GEN     pauto;
  GEN     C;
  GEN     Cd;
};
static GEN
galoisdolift(struct galois_lift *gl, GEN frob)
{
  pari_sp ltop=avma;
  long v = varn(gl->T);
  GEN Tp = FpX_red(gl->T, gl->p);
  GEN S = FpXQ_pow(pol_x[v],gl->p, Tp,gl->p);
  GEN plift = automorphismlift(S, gl, frob);
  return gerepileupto(ltop,plift);
}

static void
inittestlift( GEN plift, GEN Tmod, struct galois_lift *gl, 
    struct galois_testlift *gt)
{
  long v = varn(gl->T);
  gt->n = lg(gl->L) - 1;
  gt->g = lg(Tmod) - 1;
  gt->f = gt->n / gt->g;
  gt->bezoutcoeff = bezout_lift_fact(gl->T, Tmod, gl->p, gl->e);
  gt->pauto = cgetg(gt->f + 1, t_VEC);
  gel(gt->pauto,1) = pol_x[v];
  gel(gt->pauto,2) = gcopy(plift);
  if (gt->f > 2)
  {
    pari_sp ltop=avma;
    long i;
    long nautpow=brent_kung_optpow(gt->n-1,gt->f-2);
    GEN autpow;
    if (DEBUGLEVEL >= 1) (void)timer2();
    autpow = FpXQ_powers(plift,nautpow,gl->TQ,gl->Q);
    for (i = 3; i <= gt->f; i++)
      gel(gt->pauto,i) = FpX_FpXQV_compo(gel(gt->pauto,i-1),autpow,gl->TQ,gl->Q);
    /*Somewhat paranoid with memory, but this function use a lot of stack.*/
    gt->pauto=gerepileupto(ltop, gt->pauto);
    if (DEBUGLEVEL >= 1) msgtimer("frobenius power");
  }
}

/* We should have 0<=x<mod. */
/* Explanation of the intheadlong technique:
 * Let B be a bound B, M a modulo M>B*2^BITS_IN_LONG,
 * 0<=a_i<M for i=1,...,n.
 * We want to test if it exists k,l, |k| < B, such that sum a_i = k +l*M
 * The trick is to write a_i*2^BITS_IN_LONG/M=b_i+c_i with b_i integer
 * and 0<=c_i<1. We get sum b_i+c_i = k*2^BITS_IN_LONG/M +l*2^BITS_IN_LONG
 * so sum b_i -l*2^BITS_IN_LONG=k*2^BITS_IN_LONG/M -sum c_i
 * We have -1<k*2^BITS_IN_LONG/M<1, 0<=c_i<1 so
 * so -n-1 < sum b_i -l*2^BITS_IN_LONG<1 so
 * n<=sum b_i -l*2^BITS_IN_LONG<=0
 * So we compute z=sum b_i [2^BITS_IN_LONG] and check if 0<=-z<=n.
 */

long intheadlong(GEN x, GEN mod)
{
  pari_sp ltop=avma;
  long res= (long) itou(divii(shifti(x,BITS_IN_LONG),mod));
  avma=ltop;
  return res;
}

GEN matheadlong(GEN W, GEN mod)
{
  long i,j;
  GEN V=cgetg(lg(W),t_MAT);
  for(i=1;i<lg(W);i++)
  {
    gel(V,i) = cgetg(lg(W[i]),t_VECSMALL);
    for(j=1;j<lg(W[i]);j++)
      mael(V,i,j)=intheadlong(gmael(W,i,j),mod);
  }
  return V;
}

long polheadlong(GEN P, long n, GEN mod)
{
  return (lg(P)>n+2)?intheadlong(gel(P,n+2),mod):0;
}
/*
 * 
 */
static long
frobeniusliftall(GEN sg, long el, GEN *psi, struct galois_lift *gl,
		 struct galois_testlift *gt, GEN frob)
{
  pari_sp av, ltop2, ltop = avma;
  long d, z, m, c, n, ord, i, j, k;
  GEN pf, u, v, C, Cd, SG, cache;
  long N1, N2, R1, Ni, Z, c_idx = gt->g-1;
  long stop = 0, hop = 0;
  GEN NN, NQ;
  m = gt->g;
  ord = gt->f;
  n = lg(gl->L) - 1;
  c = lg(sg) - 1;
  d = m / c;
  pf = cgetg(m, t_VECSMALL);
  *psi = pf;
  ltop2 = avma;
  NN = diviiexact(mpfact(m), mulsi(c, gpowgs(mpfact(d), c)));
  if (DEBUGLEVEL >= 4)
    fprintferr("GaloisConj:I will try %Z permutations\n", NN);
  N1=10000000;
  NQ=divis_rem(NN,N1,&R1);
  if (cmpiu(NQ,1000000000)>0)
  {
    pari_warn(warner,"Combinatorics too hard : would need %Z tests!\n"
	"I will skip it, but it may induce an infinite loop",NN);
    avma = ltop; *psi = NULL; return 0;
  }
  N2=itos(NQ); if(!N2) N1=R1;
  if (DEBUGLEVEL>=4)
  {
    stop=N1/20;
    (void)timer2();
  }
  avma = ltop2;
  C=gt->C;
  Cd=gt->Cd;
  v = FpXQ_mul(gel(gt->pauto, 1+el%ord), gel(gt->bezoutcoeff, m),gl->TQ,gl->Q);
  v = FpX_Fp_mul(v,gl->den,gl->Q);
  SG=cgetg(lg(sg),t_VECSMALL);
  for(i=1;i<lg(SG);i++)
    SG[i]=(el*sg[i])%ord + 1;
  cache=cgetg(m+1,t_VECSMALL);
  cache[m]=polheadlong(v,1,gl->Q);
  Z=polheadlong(v,2,gl->Q);
  for (i = 1; i < m; i++)
    pf[i] = 1 + i / d;
  av = avma;
  for (Ni = 0, i = 0; ;i++)
  {
    for (j = c_idx ; j > 0; j--)
    {
      long h;
      h=SG[pf[j]];
      if (!mael(C,h,j))
      {
	pari_sp av3=avma;
	GEN r;
	r=FpX_Fp_mul(FpXQ_mul(gel(gt->pauto,h), 
	      gel(gt->bezoutcoeff,j),gl->TQ,gl->Q),gl->den,gl->Q);
	gmael(C,h,j) = gclone(r);
	mael(Cd,h,j) = polheadlong(r,1,gl->Q);
	avma=av3;
      }
      cache[j]=cache[j+1]+mael(Cd,h,j);
    }
    if (-(ulong)cache[1]<=(ulong)n)
    {
      long ZZ=Z;
      for (j = 1; j < m; j++)
	ZZ += polheadlong(gmael(C,SG[pf[j]],j),2,gl->Q);
      if (-(ulong)ZZ<=(ulong)n)
      {
	u = v;
	for (j = 1; j < m; j++)
	  u = ZX_add(u, gmael(C,SG[pf[j]],j));
	u = FpX_center(FpX_red(u, gl->Q), gl->Q);
	if (poltopermtest(u, gl, frob))
	{
	  if (DEBUGLEVEL >= 4 )
	  {
	    msgtimer("");
	    fprintferr("GaloisConj: %d hops on %Z tests\n",hop,addis(mulss(Ni,N1),i));
	  }
	  avma = ltop2;
	  return 1;
	}
	else if (DEBUGLEVEL >= 4 )
	  fprintferr("M");
      }
      else hop++;
    }
    if (DEBUGLEVEL >= 4 && i==stop)
    {
      stop+=N1/20;
      msgtimer("GaloisConj:Testing %Z", addis(mulss(Ni,N1),i));
    }
    avma = av;
    if (i == N1 - 1)
    {
      if (Ni==N2-1)
	N1=R1;
      if (Ni==N2)
	break;
      Ni++;
      i=0;
      if (DEBUGLEVEL>=4)
      {
	stop=N1/20;
	(void)timer2();
      }
    }
    for (j = 2; j < m && pf[j - 1] >= pf[j]; j++);
    for (k = 1; k < j - k && pf[k] != pf[j - k]; k++)
    {
      z = pf[k];
      pf[k] = pf[j - k];
      pf[j - k] = z;
    }
    for (k = j - 1; pf[k] >= pf[j]; k--);
    z = pf[j];
    pf[j] = pf[k];
    pf[k] = z;
    c_idx=j;
  }
  if (DEBUGLEVEL>=4)
    fprintferr("GaloisConj: not found, %d hops \n",hop);
  *psi = NULL;
  avma = ltop; return 0;
}

/* structure containing all data for permutation test:
 * 
 * order :ordre des tests pour galois_test_perm order[lg(ordre)]: numero du test
 * principal borne : borne sur les coefficients a trouver ladic: modulo
 * l-adique des racines lborne:ladic-borne TM:vecteur des ligne de M
 * PV:vecteur des clones des matrices de test (Vmatrix) (ou NULL si non
 * calcule) L,M comme d'habitude (voir plus bas)
 */
struct galois_test
{
  GEN     order;
  GEN     borne, lborne, ladic;
  GEN     PV, TM;
  GEN     L, M;
};
/* Calcule la matrice de tests correspondant a la n-ieme ligne de V */
static GEN
Vmatrix(long n, struct galois_test *td)
{
  pari_sp ltop = avma;
  GEN     V;
  long    i;
  V = cgetg(lg(td->L), t_VEC);
  for (i = 1; i < lg(V); i++) gel(V,i) = gmael(td->M,i,n);
  V = FpC_FpV_mul(td->L, V, td->ladic);
  return gerepileupto(ltop, V);
}

/*
 * Initialise la structure galois_test
 */
static void
inittest(GEN L, GEN M, GEN borne, GEN ladic, struct galois_test *td)
{
  pari_sp ltop;
  long i, n = lg(L) - 1;
  if (DEBUGLEVEL >= 8)
    fprintferr("GaloisConj:Entree Init Test\n");
  td->order = cgetg(n + 1, t_VECSMALL);
  for (i = 1; i <= n - 2; i++)
    td->order[i] = i + 2;
  for (; i <= n; i++)
    td->order[i] = i - n + 2;
  td->borne = borne;ltop = avma;
  td->lborne = subii(ladic, borne);
  td->ladic = ladic;
  td->L = L;
  td->M = M;
  td->PV = cgetg(n + 1, t_VECSMALL);
  for (i = 1; i <= n; i++)
    td->PV[i] = 0;
  ltop = avma;
  gel(td->PV, td->order[n]) = gclone(Vmatrix(td->order[n], td));
  avma = ltop;
  td->TM = shallowtrans(M);
  settyp(td->TM, t_VEC);
  for (i = 1; i < lg(td->TM); i++)
    settyp(td->TM[i], t_VEC);
  if (DEBUGLEVEL >= 8)
    fprintferr("GaloisConj:Sortie Init Test\n");
}

/* liberer les clones de la structure galois_test */
static void
freetest(struct galois_test *td)
{
  long i;
  for (i = 1; i < lg(td->PV); i++)
    if (td->PV[i])
    {
      gunclone(gel(td->PV,i));
      td->PV[i] = 0;
    }
}

/*
 * Check if the integer P seen as a p-adic number is close from an integer less
 * than td->borne in absolute value. 
 */
static long
padicisint(GEN P, struct galois_test *td)
{
  pari_sp ltop = avma;
  GEN U  = modii(P, td->ladic);
  long r = cmpii(U, td->borne) <= 0 || cmpii(U, td->lborne) >= 0;
  avma = ltop;
  return r;
}

/*
 * Check if the permutation pf is valid according to td.
 * If not, update td to make subsequent test faster (hopefully).
 */
static long
galois_test_perm(struct galois_test *td, GEN pf)
{
  pari_sp av = avma;
  GEN     P, V;
  long i, j, n = lg(td->L) - 1;
  P = perm_mul(td->L, pf);
  for (i = 1; i < n; i++)
  {
    long    ord;
    GEN     PW;
    ord = td->order[i];
    PW = gel(td->PV, ord);
    if (PW)
    {
      V = gmael(PW,1,pf[1]);
      for (j = 2; j <= n; j++)
	V = addii(V, gmael(PW,j,pf[j]));
    }
    else
      V = centermod(FpV_FpC_mul(gel(td->TM,ord), P, td->ladic), td->ladic);
    if (!padicisint(V, td))
      break;
  }
  if (i == n)
  {
    avma = av;
    return 1;
  }
  if (!td->PV[td->order[i]])
  {
    gel(td->PV, td->order[i]) = gclone(Vmatrix(td->order[i], td));
    if (DEBUGLEVEL >= 4)
      fprintferr("M");
  }
  if (DEBUGLEVEL >= 4)
    fprintferr("%d.", i);
  if (i > 1)
  {
    long    z;
    z = td->order[i];
    for (j = i; j > 1; j--)
      td->order[j] = td->order[j - 1];
    td->order[1] = z;
    if (DEBUGLEVEL >= 8)
      fprintferr("%Z", td->order);
  }
  avma = av;
  return 0;
}
/*Compute a*b/c when a*b will overflow*/
static long muldiv(long a,long b,long c)
{
  return (long)((double)a*(double)b/c);
}

/* F = cycle decomposition of sigma, B = cycle decomposition of cl(tau).
 * Check all permutations pf who can possibly correspond to tau, such that
 * tau*sigma*tau^-1 = sigma^s and tau^d = sigma^t, where d = ord cl(tau)
 * x: vector of choices, G: vector allowing linear access to elts of F. 
 * Choices multiple of e are not changed.
 * */

static GEN
testpermutation(GEN F, GEN B, GEN x, long s, long e, long cut,
		struct galois_test *td)
{
  pari_sp av, avm = avma;
  long a, b, c, d, n, p1, p2, p3, p4, p5, p6, l1, l2, N1, N2, R1;
  long V, i, j, cx, hop = 0, start = 0;
  GEN pf, ar, G, W, NN, NQ;
  if (DEBUGLEVEL >= 1) (void)timer2();
  a = lg(F) - 1;
  b = lg(F[1]) - 1;
  c = lg(B) - 1;
  d = lg(B[1]) - 1;
  n = a * b;
  s = (b + s) % b;
  pf = cgetg(n + 1, t_VECSMALL);
  av = avma;
  ar = cgetg(a + 2, t_VECSMALL); ar[a+1]=0;
  G  = cgetg(a + 1, t_VECSMALL);
  W  = matheadlong(gel(td->PV, td->order[n]), td->ladic);
  for (cx = 1, i = 1, j = 1; cx <= a; cx++, i++)
  {
    gel(G,cx) = gel(F, coeff(B,i,j));
    if (i == d)
    {
      i = 0;
      j++;
    }
  }
  NN = divis(powuu(b, c * (d - d/e)),cut);
  if (DEBUGLEVEL >= 4)
    fprintferr("GaloisConj:I will try %Z permutations\n", NN);
  N1=1000000;
  NQ=divis_rem(NN,N1,&R1);
  if (cmpiu(NQ,100000000)>0)
  {
    avma=avm;
    pari_warn(warner,"Combinatorics too hard : would need %Z tests!\n I'll skip it but you will get a partial result...",NN);
    return perm_identity(n);
  }
  N2=itos(NQ);
  for (l2 = 0; l2 <= N2; l2++)
  {
    long nbiter = (l2<N2) ? N1: R1;
    if (DEBUGLEVEL >= 2 && N2)
      fprintferr("%d%% ", muldiv(l2,100,N2));
    for (l1 = 0; l1 < nbiter; l1++)
    {
      if (start)
      {
	for (i = 1, j = e; i < a;)
	{
	  if ((++(x[i])) != b)
	    break;
	  x[i++] = 0;
	  if (i == j) { i++; j += e; }
	}
      }
      else {start=1; i = a - 1;}
      /* p5 = (p1 % d) - 1 */
      for (p1 = i + 1, p5 = p1 % d - 1 ; p1 >= 1; p1--, p5--)
      { 
        if (p5 == - 1)
        {
          p5 = d - 1;
          p6 = p1 + 1 - d;
        }
        else
          p6 = p1 + 1;
        p4 = p5 ? x[p1 - 1] : 0;
        V = 0;
        for (p2 = 1 + p4, p3 = 1 + x[p1]; p2 <= b; p2++)
        {
          V += mael(W,mael(G,p6,p3),mael(G,p1,p2));
          p3 += s;
          if (p3 > b)
            p3 -= b;
        }
        p3 = 1 + x[p1] - s;
        if (p3 <= 0)
          p3 += b;
        for (p2 = p4; p2 >= 1; p2--)
        {
          V += mael(W,mael(G,p6,p3),mael(G,p1,p2));
          p3 -= s;
          if (p3 <= 0)
            p3 += b;
        }
        ar[p1] = ar[p1+1] + V;
      }

      if (-(ulong)ar[1]<=(ulong)n)
      {
        for (p1 = 1, p5 = d; p1 <= a; p1++, p5++)
        {
          if (p5 == d)
          {
            p5 = 0;
            p4 = 0;
          }
          else
            p4 = x[p1 - 1];
          if (p5 == d - 1)
            p6 = p1 + 1 - d;
          else
            p6 = p1 + 1;
          for (p2 = 1 + p4, p3 = 1 + x[p1]; p2 <= b; p2++)
          {
            pf[mael(G,p1,p2)] = mael(G,p6,p3);
            p3 += s;
            if (p3 > b)
              p3 -= b;
          }
          p3 = 1 + x[p1] - s;
          if (p3 <= 0)
            p3 += b;
          for (p2 = p4; p2 >= 1; p2--)
          {
            pf[mael(G,p1,p2)] = mael(G,p6,p3);
            p3 -= s;
            if (p3 <= 0)
              p3 += b;
          }
        }
        if (galois_test_perm(td, pf))
        {
          if (DEBUGLEVEL >= 1)
          {
            GEN nb=addis(mulss(l2,N1),l1);
            msgtimer("testpermutation(%Z)", nb);
            if (DEBUGLEVEL >= 2 && hop)
              fprintferr("GaloisConj:%d hop sur %Z iterations\n", hop, nb);
          }
          avma = av;
          return pf;
        }
        else
          hop++;
      }
    }
  }
  if (DEBUGLEVEL >= 1)
  {
    msgtimer("testpermutation(%Z)", NN);
    if (DEBUGLEVEL >= 2 && hop)
      fprintferr("GaloisConj:%d hop sur %Z iterations\n", hop, NN);
  }
  avma = avm;
  return NULL;
}

/* List of subgroups of (\ZZ/m\ZZ)^* whose order divide p, and return the list
 * of their elements */
GEN
listznstarelts(long m, long p)
{
  pari_sp ltop = avma;
  GEN zn, zns, lss, res;
  long k, card, i, phi;
  if (m == 2)
  {
    res = cgetg(2, t_VEC);
    gel(res,1) = mkvecsmall(1);
    return res;
  }
  zn = znstar(stoi(m));
  phi = itos(gel(zn,1));
  zns = znstar_small(zn);
  lss = subgrouplist(gel(zn,2), NULL);
  res = cgetg(lg(lss), t_VEC);
  for (k = 1, i = lg(lss) - 1; i >= 1; i--)
  {
    pari_sp av;
    av = avma;
    card = phi / itos(dethnf_i(gel(lss,i)));
    avma = av;
    if (p % card == 0)
      gel(res,k++) = znstar_hnf_elts(zns,gel(lss,i));
  }
  setlg(res,k);
  return gerepileupto(ltop, gen_sort(res,0,&pari_compare_lg));
}
/* A sympol is a symmetric polynomial
 *
 * Currently sympol are couple of t_VECSMALL [v,w]
 * v[1]...v[k], w[1]...w[k]  represent the polynomial
 * sum(i=1,k,v[i]*s_w[i]) where s_i(X_1,...,X_n)=sum(j=1,n,X_j^i)
 */

/*Return s_e*/

static GEN
sympol_eval_newtonsum(long e, GEN O, GEN mod)
{
  long f,g;
  long i,j;
  GEN PL;
  f=lg(O)-1;
  g=lg(O[1])-1;
  PL=cgetg(lg(O), t_COL);
  for(i=1; i<=f; i++)
  {
    pari_sp ltop=avma;
    GEN s=gen_0;
    for(j=1; j<=g; j++)
      s=addii(s,Fp_powu(gmael(O,i,j),(ulong) e,mod));
    gel(PL,i) = gerepileupto(ltop,modii(s,mod));
  }
  return PL;
}

GEN
sympol_eval(GEN v, GEN NS)
{
  pari_sp ltop=avma;
  long i;
  GEN S=gen_0;
  for(i=1;i<lg(v);i++)
    if (v[i]) S=gadd(S,gmulsg(v[i],gel(NS,i)));
  return gerepileupto(ltop, S);
}

/*
 * Let sigma be an automorphism of L (as a polynomial with rational coefs)
 * Let 'sym' be a symmetric polynomial defining alpha in L.
 * We have alpha=sym(x,sigma(x),,,sigma^(g-1)(x))
 * Compute alpha mod p.
 */

GEN
sympol_aut_evalmod(GEN sym, long g, GEN sigma, GEN Tp, GEN p)
{
  pari_sp ltop=avma;
  long i, j, npows;
  GEN  s, f, pows;
  GEN v=gel(sym,1), w=gel(sym,2);
  sigma = RgX_to_FpX(sigma, p);
  f=pol_x[varn(sigma)];
  s=zeropol(varn(sigma));
  for(j=1; j<lg(v); j++)
    s=FpX_add(s,FpX_Fp_mul(FpXQ_pow(f,stoi(w[j]),Tp,p),stoi(v[j]),p),p);
  npows = brent_kung_optpow(lg(Tp)-4,g-1);
  pows  = FpXQ_powers(sigma,npows,Tp,p);
  for(i=2; i<=g;i++)
  {
    f=FpX_FpXQV_compo(f,pows,Tp,p);
    for(j=1; j<lg(v); j++)
      s=FpX_add(s,FpX_Fp_mul(FpXQ_pow(f,stoi(w[j]),Tp,p),stoi(v[j]),p),p);
  }
  return gerepileupto(ltop, s);
}

/* Let Sp be as computed with sympol_aut_evalmod
 * Let Tmod be the factorisation of T mod p.
 * Return the factorisation of the minimal polynomial of S
 * mod p w.r.t. Tmod.
 */

GEN
fixedfieldfactmod(GEN Sp, GEN p, GEN Tmod)
{
  long i;
  long l=lg(Tmod);
  GEN F=cgetg(l,t_VEC);
  for(i=1;i<l;i++)
    gel(F,i) = FpXQ_minpoly(FpX_rem(Sp,gel(Tmod,i),p), gel(Tmod,i),p);
  return F;
}

static GEN
fixedfieldsurmer(GEN O, GEN mod, GEN l, GEN p, long v, GEN NS, GEN W)
{
  long i,j;
  const long step=3;
  long n=lg(W)-1;
  long m=1<<((n-1)<<1);
  GEN sym=cgetg(n+1,t_VECSMALL);
  for (j=1;j<n;j++) sym[j]=step;
  sym[n]=0;
  if (DEBUGLEVEL>=4) fprintferr("FixedField: Weight: %Z\n",W);
  for (i=0;i<m;i++)
  {
    pari_sp av=avma;
    GEN L,P;
    for (j=1;sym[j]==step;j++)
      sym[j]=0;
    sym[j]++;
    if (DEBUGLEVEL>=6) fprintferr("FixedField: Sym: %Z\n",sym);
    L=sympol_eval(sym,NS);
    if (!vec_is1to1(FpC_red(L,l))) continue;
    P=FpX_center(FpV_roots_to_pol(L,mod,v),mod);
    if (!p || FpX_is_squarefree(P,p))
      return mkvec3(mkvec2(sym,W),L,P);
    avma=av;
  }
  return NULL;
}

/*Check whether the line of NS are pair-wise distinct.*/

static long
sympol_is1to1_lg(GEN NS, long n)
{
  long i,j,k;
  long l=lg(NS[1]);
  for (i=1;i<l;i++)
    for(j=i+1;j<l;j++)
    {
      for(k=1;k<n;k++)
        if (!equalii(gmael(NS,k,j),gmael(NS,k,i)))
          break;
      if (k>=n)
        return 0;
    }
  return 1;
}

/* Let O a set of orbits of roots (see fixedfieldorbits) modulo mod,
 * l and p two prime number, l dividing mod
 * Return a vector [sym,s,P] where:
 * s is a sympol, s is the set of images of S on O and 
 * P is the polynomial with roots s.
 */

GEN
fixedfieldsympol(GEN O, GEN mod, GEN l, GEN p, long v)
{
  pari_sp ltop=avma;
  const long n=(BITS_IN_LONG>>1)-1;
  GEN NS=cgetg(n+1,t_MAT);
  GEN sym=NULL, W=cgetg(n+1,t_VECSMALL);
  long i, e=1;
  if (DEBUGLEVEL>=4) 
    fprintferr("FixedField: Size: %ldx%ld\n",lg(O)-1,lg(O[1])-1);
  for(i=1;!sym && i<=n; i++)
  {
    GEN L = sympol_eval_newtonsum(e++, O, mod);
    if (lg(O)>2)
      while (vec_isconst(L))
        L = sympol_eval_newtonsum(e++, O, mod);
    W[i] = e-1; gel(NS,i) = L;
    if (sympol_is1to1_lg(NS,i+1))
      sym=fixedfieldsurmer(O,mod,l,p,v,NS,vecsmall_shorten(W,i));
  }
  if (!sym) pari_err(talker,"p too small in fixedfieldsympol");
  if (DEBUGLEVEL>=2) fprintferr("FixedField: Found: %Z\n",gel(sym,1));
  return gerepilecopy(ltop,sym);
}

/* Let O a set of orbits as indices and L the corresponding roots.
 * Return the set of orbits as roots.
 */

GEN
fixedfieldorbits(GEN O, GEN L)
{
  GEN S = cgetg(lg(O), t_MAT);
  long i, j;
  for (i = 1; i < lg(O); i++)
  {
    GEN z = cgetg(lg(O[i]),t_COL);
    gel(S,i) = z;
    for (j = 1; j < lg(O[i]); j++)
      z[j] = L[mael(O,i,j)];
  }
  return S;
}

GEN
fixedfieldinclusion(GEN O, GEN PL)
{
  GEN S = cgetg((lg(O) - 1) * (lg(O[1]) - 1) + 1, t_COL);
  long i, j;
  for (i = 1; i < lg(O); i++)
    for (j = 1; j < lg(O[i]); j++)
      S[mael(O,i,j)] = PL[i];
  return S;
}

/*Usually mod is bigger than than den so there is no need to reduce it.*/
GEN
vandermondeinversemod(GEN L, GEN T, GEN den, GEN mod)
{
  pari_sp av;
  long i, j, n = lg(L);
  long x = varn(T);
  GEN M, P, Tp;
  M = cgetg(n, t_MAT);
  av=avma;
  Tp = gclone(FpX_deriv(T,mod)); /*clone*/
  avma=av;
  for (i = 1; i < n; i++)
  {
    GEN z;
    av = avma;
    z = Fp_inv(FpX_eval(Tp, gel(L,i),mod),mod);
    z = modii(mulii(den,z),mod);
    P = FpX_Fp_mul(FpX_div(T, deg1pol_i(gen_1,negi(gel(L,i)),x),mod), z, mod); 
    gel(M,i) = cgetg(n, t_COL);
    for (j = 1; j < n; j++)
      gmael(M,i,j) = gcopy(gel(P,1 + j));
    gel(M,i) = gerepileupto(av,gel(M,i));
  }
  gunclone(Tp); /*unclone*/
  return M;
}
/* Calcule le polynome associe a one vecteur conjugue v*/
static GEN
vectopol(GEN v, GEN M, GEN den , GEN mod, long x)
{
  long n = lg(v), i, k;
  pari_sp av;
  GEN z = cgetg(n+1,t_POL),p1,mod2;
  av=avma;
  mod2=gclone(shifti(mod,-1));/*clone*/
  avma=av;
  z[1] = evalsigne(1)|evalvarn(x);
  for (i=2; i<=n; i++)
  {
    p1=gen_0; av=avma;
    for (k=1; k<n; k++)
      p1 = addii(p1, mulii(gcoeff(M,i-1,k),gel(v,k)));
    p1=modii(p1,mod);
    if (cmpii(p1,mod2)>0) p1=subii(p1,mod);
    gel(z,i) = gerepileupto(av, gdiv(p1,den));
  }
  gunclone(mod2);/*unclone*/
  return normalizepol_i(z,n+1);
}
/* Calcule le polynome associe a une permutation des racines*/
static GEN
permtopol(GEN p, GEN L, GEN M, GEN den, GEN mod, long x)
{
  long n = lg(L), i, k;
  pari_sp av;
  GEN z = cgetg(n+1,t_POL),p1,mod2;
  if (lg(p) != n) pari_err(talker,"incorrect permutation in permtopol");
  av=avma;
  mod2=gclone(shifti(mod,-1)); /*clone*/
  avma=av;
  z[1] = evalsigne(1)|evalvarn(x);
  for (i=2; i<=n; i++)
  {
    p1=gen_0; av=avma;
    for (k=1; k<n; k++)
      p1 = addii(p1, mulii(gcoeff(M,i-1,k), gel(L,p[k])));
    p1=modii(p1,mod);
    if (cmpii(p1,mod2)>0) p1=subii(p1,mod);
    gel(z,i) = gerepileupto(av, gdiv(p1,den));
  }
  gunclone(mod2); /*unclone*/
  return normalizepol_i(z,n+1);
}

static GEN
galoisgrouptopol( GEN res, GEN L, GEN M, GEN den, GEN mod, long v)
{
  GEN aut = cgetg(lg(res), t_COL);
  long i;
  for (i = 1; i < lg(res); i++)
  {
    if (DEBUGLEVEL>=6)
      fprintferr("%d ",i);
    gel(aut,i) = permtopol(gel(res,i), L, M, den, mod, v);
  }
  return aut;
}

/* contains the result of the study of Frobenius degrees */
enum ga_code {ga_all_normal=1,ga_ext_2=2,ga_non_wss=4};
struct galois_analysis
{
  long    p; /* prime to be lifted */
  long    deg; /* degree of the lift */
  long    ord;
  long    l; /* l: prime number such that T is totally split mod l */
  long    p4;
  enum ga_code group;
  byteptr primepointer; /* allow computing the primes following p */
};

static void
galoisanalysis(GEN T, struct galois_analysis *ga, long calcul_l)
{
  pari_sp ltop=avma;
  long n,p;
  long i;
  long karma=0;
  long group,linf;
  /*TODO: complete the table to at least 200*/
  const long prim_nonss_orders[]={36,48,56,60,72,75,80,96,108,120,132,0};
  GEN F,Fp,Fe,Fpe,O;
  long min_prime,np;
  long order,phi_order;
  long plift,nbmax,nbtest,deg;
  byteptr primepointer,pp;

  if (!ZX_is_squarefree(T))
    pari_err(talker, "Polynomial not squarefree in galoisinit");
  if (DEBUGLEVEL >= 1) (void)timer2();
  n = degpol(T);
  O = cgetg(n+1,t_VECSMALL);
  for(i=1;i<=n;i++) O[i]=0;
  F = factoru_pow(n);
  Fp =gel(F,1);
  Fe =gel(F,2);
  Fpe=gel(F,3);
  np=lg(Fp)-1;
  /*In this part, we study the cardinal of the group to have an information
    about the orders, so if we are unlucky we can continue.*/

  /*Are there non WSS groups of this order ?*/
  group=0;
  for(i=0;prim_nonss_orders[i];i++)
    if (n%prim_nonss_orders[i] == 0)
    {
      group |= ga_non_wss;
      break;
    }
  if ( n>12 && n%12 == 0 )
  {
    /*We need to know the greatest prime dividing n/12*/
    if ( Fp[np] == 3 && Fe[np] == 1 )
      group |= ga_ext_2;
  }
  phi_order = 1;
  order = 1;
  for (i = np; i > 0; i--)
  {
    p = Fp[i];
    if (phi_order % p != 0)
    {
      order *= p;
      phi_order *= p - 1;
    }
    else
    {
      group |= ga_all_normal;
      break;
    }
    if (Fe[i]>1)
      break;
  }
  /*Now, we study the orders of the Frobenius elements*/
  min_prime=n*max((long)(BITS_IN_LONG-bfffo(n)-4),2);
  plift = 0;
  nbmax = 8+(n>>1);
  nbtest = 0; 
  deg = Fp[np];
  for (p = 0, pp = primepointer = diffptr;
       (plift == 0 
	|| (nbtest < nbmax && (nbtest <=8 || order < (n>>1)))
	|| (n == 24 && O[6] == 0 && O[4] == 0)
        || ((group&ga_non_wss) && order == Fp[np]))
         && (nbtest < 3 * nbmax || !(group&ga_non_wss)) ;)
  {
    pari_sp av;
    GEN ip, FS;
    long d, o, norm_o = 1;

    NEXT_PRIME_VIADIFF_CHECK(p,primepointer);
    /*discard small primes*/
    if (p <= min_prime)
      continue;
    ip = utoipos(p);
    if (!FpX_is_squarefree(T,ip))
      continue;
    nbtest++;
    av=avma;
    FS=(GEN)FpX_degfact(T,ip)[1];
    d = FS[1];
    for(i=2;i<lg(FS);i++)
      if (d != FS[i]) break;
    if (i<lg(FS))
    {
      if (DEBUGLEVEL >= 2)
	fprintferr("GaloisAnalysis:non Galois for p=%ld\n", p);
      ga->p = p;
      ga->deg = 0;
      avma = ltop; return; /* Not a Galois polynomial */
    }
    o=n/(lg(FS)-1);
    avma=av;
    if (!O[o]) O[o]=p;
    if (o % deg == 0)
    {
      /*We try to find a power of the Frobenius which generate
	a normal subgroup just by looking at the order.*/
      if (o * Fp[1] >= n)
	/*Subgroup of smallest index are normal*/
	norm_o = o;
      else		
      {
	norm_o = 1;
	for (i = np; i > 0; i--)
	{
	  if (o % Fpe[i] == 0)
	    norm_o *= Fpe[i];
	  else
	    break;
	}
      }
      if (norm_o != 1)
      {
	if (!(group&ga_all_normal) || o > order || 
	    (o == order && (plift == 0 || norm_o > deg 
			    || (norm_o == deg && cgcd(p-1,n) > (long)karma ))))
	{
	  deg = norm_o;
	  order = o;
	  plift = p;
	  karma=cgcd(p-1,n);
	  pp = primepointer;
	  group |= ga_all_normal;
	}
      }
      else if (!(group&ga_all_normal) && (plift == 0 || o > order 
	    || ( o == order && cgcd(p-1,n) > (long)karma )))
      {
	order = o;
	plift = p;
	karma=cgcd(p-1,n);
	pp = primepointer;
      }
    }
    if (DEBUGLEVEL >= 5)
      fprintferr("GaloisAnalysis:Nbtest=%ld,p=%ld,o=%ld,n_o=%d,best p=%ld,ord=%ld,k=%ld\n",
		 nbtest, p, o, norm_o, plift, order,karma);
  }
  /* This is to avoid looping on non-wss group. 
     To be checked for large groups.  */
  /* Would it be better to disable this check if we are in a good case ?
   * (ga_all_normal and !(ga_ext_2) (e.g. 60)) ?*/
  if (plift == 0 || ((group&ga_non_wss) && order == Fp[np]))
  {
    deg = 0;
    pari_warn(warner, "Galois group almost certainly not weakly super solvable");
  }
  /*linf=(n*(n-1))>>2;*/
  linf=n;
  if (calcul_l && O[1]<=linf)
  {
    pari_sp av;
    long    l=0;
    /*we need a totally split prime l*/
    av = avma;
    while (l == 0)
    {
      long nb;
      GEN Tp;

      NEXT_PRIME_VIADIFF_CHECK(p,primepointer);
      if (p <= linf) continue;
      Tp = ZX_to_Flx(T, p);
      nb = Flx_nbroots(Tp, p);
      if (nb == n)
	l = p;
      else if (nb && Flx_is_squarefree(Tp, p))
      {
	avma = ltop;
	if (DEBUGLEVEL >= 2)
	  fprintferr("GaloisAnalysis:non Galois for p=%ld\n", p);
	ga->p = p;
	ga->deg = 0;
	return;	/* Not a Galois polynomial */
      }
      avma = av;
    }
    O[1]=l;
  }  
  ga->p = plift;
  ga->group = (enum ga_code)group;
  ga->deg = deg;
  ga->ord = order;
  ga->l = O[1];
  ga->primepointer = pp;
  ga->p4 = O[4];
  if (DEBUGLEVEL >= 4)
    fprintferr("GaloisAnalysis:p=%ld l=%ld group=%ld deg=%ld ord=%ld\n",
	       plift, O[1], group, deg, order);
  if (DEBUGLEVEL >= 1)
    msgtimer("galoisanalysis()");
  avma = ltop;
}

/* Groupe A4 */
static GEN
a4galoisgen(GEN T, struct galois_test *td)
{
  pari_sp ltop = avma, av, av2;
  long    i, j, k;
  long    n;
  long    N, hop = 0;
  GEN     O, ar, mt;
  GEN     t, u;
  GEN     res, orb, ry;
  GEN     pft, pfu, pfv;
  n = degpol(T);
  res = cgetg(3, t_VEC);
  ry = cgetg(4, t_VEC);
  gel(res,1) = ry;
  pft = cgetg(n + 1, t_VECSMALL);
  pfu = cgetg(n + 1, t_VECSMALL);
  pfv = cgetg(n + 1, t_VECSMALL);
  gel(ry,1) = pft;
  gel(ry,2) = pfu;
  gel(ry,3) = pfv;
  ry = cgetg(4, t_VECSMALL);
  ry[1] = 2;
  ry[2] = 2;
  ry[3] = 3;
  gel(res,2) = ry;
  av = avma;
  ar = cgetg(n+1, t_VEC);
  for (i = 1; i <= n; i++) gel(ar,i) = cgeti(1 + lg(td->ladic));
  mt = gel(td->PV,td->order[n]);
  t = cgetg(n + 1, t_VECSMALL) + 1;	/* Sorry for this hack */
  u = cgetg(n + 1, t_VECSMALL) + 1;	/* too lazy to correct */
  av2 = avma;
  N = itos(gdiv(mpfact(n), mpfact(n >> 1))) >> (n >> 1);
  if (DEBUGLEVEL >= 4)
    fprintferr("A4GaloisConj:I will test %ld permutations\n", N);
  avma = av2;
  for (i = 0; i < n; i++)
    t[i] = i + 1;
  for (i = 0; i < N; i++)
  {
    GEN     g;
    long     a, x, y;
    if (i == 0)
    {
      affsi(0, gel(ar,(n - 2) >> 1));
      for (k = n - 2; k > 2; k -= 2)
	addiiz(gel(ar,k >> 1), addii(gmael(mt,k + 1,k + 2), gmael(mt,k + 2,k + 1)),
	       gel(ar,(k >> 1) - 1));
    }
    else
    {
      x = i;
      y = 1;
      do
      {
	y += 2;
	a = x%y;
	x = x/y;
      }
      while (!a);
      switch (y)
      {
      case 3:
	x = t[2];
	if (a == 1)
	{
	  t[2] = t[1];
	  t[1] = x;
	}
	else
	{
	  t[2] = t[0];
	  t[0] = x;
	}
	break;
      case 5:
	x = t[0];
	t[0] = t[2];
	t[2] = t[1];
	t[1] = x;
	x = t[4];
	t[4] = t[4 - a];
	t[4 - a] = x;
	addiiz(gel(ar,2), addii(gmael(mt,t[4],t[5]), gmael(mt,t[5],t[4])), gel(ar,1));
	break;
      case 7:
	x = t[0];
	t[0] = t[4];
	t[4] = t[3];
	t[3] = t[1];
	t[1] = t[2];
	t[2] = x;
	x = t[6];
	t[6] = t[6 - a];
	t[6 - a] = x;
	addiiz(gel(ar,3), addii(gmael(mt,t[6],t[7]), gmael(mt,t[7],t[6])), gel(ar,2));
	addiiz(gel(ar,2), addii(gmael(mt,t[4],t[5]), gmael(mt,t[5],t[4])), gel(ar,1));
	break;
      case 9:
	x = t[0];
	t[0] = t[6];
	t[6] = t[5];
	t[5] = t[3];
	t[3] = x;
	x = t[4];
	t[4] = t[1];
	t[1] = x;
	x = t[8];
	t[8] = t[8 - a];
	t[8 - a] = x;
	addiiz(gel(ar,4), addii(gmael(mt,t[8],t[9]), gmael(mt,t[9],t[8])), gel(ar,3));
	addiiz(gel(ar,3), addii(gmael(mt,t[6],t[7]), gmael(mt,t[7],t[6])), gel(ar,2));
	addiiz(gel(ar,2), addii(gmael(mt,t[4],t[5]), gmael(mt,t[5],t[4])), gel(ar,1));
	break;
      default:
	y--;
	x = t[0];
	t[0] = t[2];
	t[2] = t[1];
	t[1] = x;
	for (k = 4; k < y; k += 2)
	{
	  long j;
	  x = t[k];
	  for (j = k; j > 0; j--)
	    t[j] = t[j - 1];
	  t[0] = x;
	}
	x = t[y];
	t[y] = t[y - a];
	t[y - a] = x;
	for (k = y; k > 2; k -= 2)
	  addiiz(gel(ar,k >> 1),
		addii(gmael(mt,t[k],t[k + 1]), gmael(mt,t[k + 1],t[k])),
		gel(ar,(k >> 1) - 1));
      }
    }
    g = addii(gel(ar,1), addii(addii(gmael(mt,t[0],t[1]), gmael(mt,t[1],t[0])),
			 addii(gmael(mt,t[2],t[3]), gmael(mt,t[3],t[2]))));
    if (padicisint(g, td))
    {
      for (k = 0; k < n; k += 2)
      {
	pft[t[k]] = t[k + 1];
	pft[t[k + 1]] = t[k];
      }
      if (galois_test_perm(td, pft))
	break;
      else
	hop++;
    }
    avma = av2;
  }
  if (i == N)
  {
    avma = ltop;
    if (DEBUGLEVEL >= 1 && hop)
      fprintferr("A4GaloisConj: %ld hop sur %ld iterations\n", hop, N);
    return gen_0;
  }
  if (DEBUGLEVEL >= 1 && hop)
    fprintferr("A4GaloisConj: %ld hop sur %ld iterations\n", hop, N);
  N = itos(gdiv(mpfact(n >> 1), mpfact(n >> 2))) >> 1;
  avma = av2;
  if (DEBUGLEVEL >= 4)
    fprintferr("A4GaloisConj:sigma=%Z \n", pft);
  for (i = 0; i < N; i++)
  {
    GEN g;
    long a, x, y;
    if (i == 0)
    {
      for (k = 0; k < n; k += 4)
      {
	u[k + 3] = t[k + 3];
	u[k + 2] = t[k + 1];
	u[k + 1] = t[k + 2];
	u[k] = t[k];
      }
    }
    else
    {
      x = i;
      y = -2;
      do
      {
	y += 4;
	a = x%y;
	x = x/y;
      }
      while (!a);
      x = u[2];
      u[2] = u[0];
      u[0] = x;
      switch (y)
      {
      case 2:
	break;
      case 6:
	x = u[4];
	u[4] = u[6];
	u[6] = x;
	if (!(a & 1))
	{
	  a = 4 - (a >> 1);
	  x = u[6];
	  u[6] = u[a];
	  u[a] = x;
	  x = u[4];
	  u[4] = u[a - 2];
	  u[a - 2] = x;
	}
	break;
      case 10:
	x = u[6];
	u[6] = u[3];
	u[3] = u[2];
	u[2] = u[4];
	u[4] = u[1];
	u[1] = u[0];
	u[0] = x;
	if (a >= 3)
	  a += 2;
	a = 8 - a;
	x = u[10];
	u[10] = u[a];
	u[a] = x;
	x = u[8];
	u[8] = u[a - 2];
	u[a - 2] = x;
	break;
      }
    }
    g = gen_0;
    for (k = 0; k < n; k += 2)
      g = addii(g, addii(gmael(mt,u[k],u[k + 1]), gmael(mt,u[k + 1],u[k])));
    if (padicisint(g, td))
    {
      for (k = 0; k < n; k += 2)
      {
	pfu[u[k]] = u[k + 1];
	pfu[u[k + 1]] = u[k];
      }
      if (galois_test_perm(td, pfu))
	break;
      else
	hop++;
    }
    avma = av2;
  }
  if (i == N)
  {
    avma = ltop;
    return gen_0;
  }
  if (DEBUGLEVEL >= 1 && hop)
    fprintferr("A4GaloisConj: %ld hop sur %ld iterations\n", hop, N);
  if (DEBUGLEVEL >= 4)
    fprintferr("A4GaloisConj:tau=%Z \n", pfu);
  avma = av2;
  orb = cgetg(3, t_VEC);
  gel(orb,1) = pft;
  gel(orb,2) = pfu;
  if (DEBUGLEVEL >= 4)
    fprintferr("A4GaloisConj:orb=%Z \n", orb);
  O = vecperm_orbits(orb, 12);
  if (DEBUGLEVEL >= 4)
    fprintferr("A4GaloisConj:O=%Z \n", O);
  av2 = avma;
  for (j = 0; j < 2; j++)
  {
    pfv[mael(O,1,1)] = mael(O,2,1);
    pfv[mael(O,1,2)] = mael(O,2,3 + j);
    pfv[mael(O,1,3)] = mael(O,2,4 - (j << 1));
    pfv[mael(O,1,4)] = mael(O,2,2 + j);
    for (i = 0; i < 4; i++)
    {
      long    x;
      GEN     g;
      switch (i)
      {
      case 0:
	break;
      case 1:
	x = mael(O,3,1);
	mael(O,3,1) = mael(O,3,2);
	mael(O,3,2) = x;
	x = mael(O,3,3);
	mael(O,3,3) = mael(O,3,4);
	mael(O,3,4) = x;
	break;
      case 2:
	x = mael(O,3,1);
	mael(O,3,1) = mael(O,3,4);
	mael(O,3,4) = x;
	x = mael(O,3,2);
	mael(O,3,2) = mael(O,3,3);
	mael(O,3,3) = x;
	break;
      case 3:
	x = mael(O,3,1);
	mael(O,3,1) = mael(O,3,2);
	mael(O,3,2) = x;
	x = mael(O,3,3);
	mael(O,3,3) = mael(O,3,4);
	mael(O,3,4) = x;
      }
      pfv[mael(O,2,1)] = mael(O,3,1);
      pfv[mael(O,2,3 + j)] = mael(O,3,4 - j);
      pfv[mael(O,2,4 - (j << 1))] = mael(O,3,2 + (j << 1));
      pfv[mael(O,2,2 + j)] = mael(O,3,3 - j);
      pfv[mael(O,3,1)] = mael(O,1,1);
      pfv[mael(O,3,4 - j)] = mael(O,1,2);
      pfv[mael(O,3,2 + (j << 1))] = mael(O,1,3);
      pfv[mael(O,3,3 - j)] = mael(O,1,4);
      g = gen_0;
      for (k = 1; k <= n; k++)
	g = addii(g, gmael(mt,k,pfv[k]));
      if (padicisint(g, td) && galois_test_perm(td, pfv))
      {
	avma = av;
	if (DEBUGLEVEL >= 1)
	  fprintferr("A4GaloisConj:%ld hop sur %d iterations max\n",
		     hop, 10395 + 68);
	return res;
      }
      else
	hop++;
      avma = av2;
    }
  }
  /* Echec? */
  avma = ltop;
  return gen_0;
}

/* Groupe S4 */
static void
s4makelift(GEN u, struct galois_lift *gl, GEN liftpow)
{
  long i;
  gel(liftpow,1) = automorphismlift(u, gl, NULL);
  for (i = 2; i < lg(liftpow); i++)
    gel(liftpow,i) = FpXQ_mul(gel(liftpow,i - 1), gel(liftpow,1),gl->TQ,gl->Q);
}
static long
s4test(GEN u, GEN liftpow, struct galois_lift *gl, GEN phi)
{
  pari_sp ltop = avma;
  GEN res;
  long bl,i,d = lg(u)-2;
  if (DEBUGLEVEL >= 6) (void)timer2();
  if ( !d ) return 0;
  res=gel(u,2);
  for (i = 1; i < d; i++)
  {
    if (lg(liftpow[i])>2)
      res=addii(res,mulii(gmael(liftpow,i,2), gel(u,i + 2))); 
  }
  res=modii(mulii(res,gl->den),gl->Q);
  if (cmpii(res,gl->gb->bornesol)>0 
      && cmpii(res,subii(gl->Q,gl->gb->bornesol))<0)
  {
    avma=ltop;
    return 0;
  }
  res = scalarpol(gel(u,2),varn(u));
  for (i = 1; i < d ; i++)
  {
    GEN z = ZX_Z_mul(gel(liftpow,i), gel(u,i + 2));
    res = FpX_add(res,z ,gl->Q);
  }
  res = FpX_center(FpX_Fp_mul(res,gl->den,gl->Q), gl->Q);
  if (DEBUGLEVEL >= 6)
    msgtimer("s4test()");
  bl = poltopermtest(res, gl, phi);
  avma=ltop;
  return bl;
}
static GEN
s4releveauto(GEN misom,GEN Tmod,GEN Tp,GEN p,long a1,long a2,long a3,long a4,long a5,long a6)
{
  pari_sp ltop=avma;
  GEN u1,u2,u3,u4,u5;
  GEN pu1,pu2,pu3,pu4;
  pu1=FpX_mul( gel(Tmod,a2), gel(Tmod,a1),p);
  u1 = FpX_chinese_coprime(gmael(misom,a1,a2),gmael(misom,a2,a1),
			 gel(Tmod,a2), gel(Tmod,a1),pu1,p);
  pu2=FpX_mul( gel(Tmod,a4), gel(Tmod,a3),p);
  u2 = FpX_chinese_coprime(gmael(misom,a3,a4),gmael(misom,a4,a3),
			 gel(Tmod,a4), gel(Tmod,a3),pu2,p);
  pu3=FpX_mul( gel(Tmod,a6), gel(Tmod,a5),p);
  u3 = FpX_chinese_coprime(gmael(misom,a5,a6),gmael(misom,a6,a5),
			 gel(Tmod,a6), gel(Tmod,a5),pu3,p);
  pu4=FpX_mul(pu1,pu2,p);
  u4 = FpX_chinese_coprime(u1,u2,pu1,pu2,pu4,p);
  u5 = FpX_chinese_coprime(u4,u3,pu4,pu3,Tp,p);
  return gerepileupto(ltop,u5);
}
static GEN
s4galoisgen(struct galois_lift *gl)
{
  struct galois_testlift gt;
  pari_sp av, ltop2, ltop = avma;
  GEN     Tmod, isom, isominv, misom;
  long i, j;
  GEN     sg;
  GEN     sigma, tau, phi;
  GEN     res, ry;
  GEN     pj;
  GEN     p,Q,TQ,Tp;
  GEN     bezoutcoeff, pauto, liftpow, aut;

  p = gl->p;
  Q = gl->Q;
  res = cgetg(3, t_VEC);
  ry  = cgetg(5, t_VEC);
  gel(res,1) = ry;
  for (i = 1; i < lg(ry); i++)
    gel(ry,i) = cgetg(lg(gl->L), t_VECSMALL);
  ry = cgetg(5, t_VECSMALL);
  gel(res,2) = ry;
  ry[1] = 2;
  ry[2] = 2;
  ry[3] = 3;
  ry[4] = 2;
  ltop2 = avma;
  sg = cgetg(7, t_VECSMALL);
  pj = cgetg(7, t_VECSMALL);
  sigma = cgetg(lg(gl->L), t_VECSMALL);
  tau = cgetg(lg(gl->L), t_VECSMALL);
  phi = cgetg(lg(gl->L), t_VECSMALL);
  for (i = 1; i < lg(sg); i++)
    sg[i] = i;
  Tp = FpX_red(gl->T,p);
  TQ = gl->TQ;
  Tmod = lift((GEN) factmod(gl->T, p)[1]);
  isom = cgetg(lg(Tmod), t_VEC);
  isominv = cgetg(lg(Tmod), t_VEC);
  misom = cgetg(lg(Tmod), t_MAT);
  aut=galoisdolift(gl, NULL);
  inittestlift(aut,Tmod, gl, &gt);
  bezoutcoeff = gt.bezoutcoeff;
  pauto = gt.pauto;
  for (i = 1; i < lg(pj); i++)
    pj[i] = 0;
  for (i = 1; i < lg(isom); i++)
  {
    gel(misom,i) = cgetg(lg(Tmod), t_COL);
    gel(isom,i) = FpX_ffisom(gel(Tmod,1), gel(Tmod,i), p);
    if (DEBUGLEVEL >= 6)
      fprintferr("S4GaloisConj:Computing isomorphisms %d:%Z\n", i,
		 gel(isom,i));
    gel(isominv,i) = FpXQ_ffisom_inv(gel(isom,i), gel(Tmod,i),p);
  }
  for (i = 1; i < lg(isom); i++)
    for (j = 1; j < lg(isom); j++)
      gmael(misom,i,j) = FpX_FpXQ_compo(gel(isominv,i),gel(isom,j),
 				        gel(Tmod,j),p);
  liftpow = cgetg(24, t_VEC);
  av = avma;
  for (i = 0; i < 3; i++)
  {
    pari_sp av2, avm1, avm2;
    GEN u;
    long j1, j2, j3;
    GEN u1, u2, u3;
    if (i)
    {
      long x;
      x = sg[3];
      if (i == 1)
      {
	sg[3] = sg[2];
	sg[2] = x;
      }
      else
      {
	sg[3] = sg[1];
	sg[1] = x;
      }
    }
    u=s4releveauto(misom,Tmod,Tp,p,sg[1],sg[2],sg[3],sg[4],sg[5],sg[6]);
    s4makelift(u, gl, liftpow);
    av2 = avma;
    for (j1 = 0; j1 < 4; j1++)
    {
      u1 = FpX_add(FpXQ_mul(gel(bezoutcoeff, sg[5]),
			    gel(pauto,1 + j1),TQ,Q),
                   FpXQ_mul(gel(bezoutcoeff, sg[6]),
                            gel(pauto, ((-j1) & 3) + 1),TQ,Q),Q);
      avm1 = avma;
      for (j2 = 0; j2 < 4; j2++)
      {
	u2 = ZX_add(u1, FpXQ_mul(gel(bezoutcoeff, sg[3]), 
				 gel(pauto,1 + j2),TQ,Q));
	u2 = FpX_add(u2, FpXQ_mul(gel(bezoutcoeff,sg[4]),
				  gel(pauto,((-j2) & 3) + 1), TQ,Q),Q);
	avm2 = avma;
	for (j3 = 0; j3 < 4; j3++)
	{
	  u3 = ZX_add(u2, FpXQ_mul(gel(bezoutcoeff, sg[1]),
				   gel(pauto,1 + j3),TQ,Q));
	  u3 = FpX_add(u3, FpXQ_mul(gel(bezoutcoeff, sg[2]),
				    gel(pauto,((-j3) & 3) + 1), TQ,Q),Q);
	  if (DEBUGLEVEL >= 4)
	    fprintferr("S4GaloisConj:Testing %d/3:%d/4:%d/4:%d/4:%Z\n",
		       i, j1,j2, j3, sg);
	  if (s4test(u3, liftpow, gl, sigma))
	  {
	    pj[1] = j3;
	    pj[2] = j2;
	    pj[3] = j1;
	    goto suites4;
	  }
	  avma = avm2;
	}
	avma = avm1;
      }
      avma = av2;
    }
    avma = av;
  }
  avma = ltop;
  return gen_0;
suites4:
  if (DEBUGLEVEL >= 4)
    fprintferr("S4GaloisConj:sigma=%Z\n", sigma);
  if (DEBUGLEVEL >= 4)
    fprintferr("S4GaloisConj:pj=%Z\n", pj);
  avma = av;
  for (j = 1; j <= 3; j++)
  {
    pari_sp av2;
    GEN     u;
    long w, l, z;
    z = sg[1]; sg[1] = sg[3]; sg[3] = sg[5]; sg[5] = z;
    z = sg[2]; sg[2] = sg[4]; sg[4] = sg[6]; sg[6] = z;
    z = pj[1]; pj[1] = pj[2]; pj[2] = pj[3]; pj[3] = z;
    for (l = 0; l < 2; l++)
    {
      u=s4releveauto(misom,Tmod,Tp,p,sg[1],sg[3],sg[2],sg[4],sg[5],sg[6]);
      s4makelift(u, gl, liftpow);
      av2 = avma;
      for (w = 0; w < 4; w += 2)
      {
	pari_sp av3;
	GEN     uu;
	pj[6] = (w + pj[3]) & 3;
	uu =FpX_add(FpXQ_mul(gel(bezoutcoeff,sg[5]),
			     gel(pauto,(pj[6] & 3) + 1), TQ,Q),
                    FpXQ_mul(gel(bezoutcoeff,sg[6]),
			     gel(pauto,((-pj[6]) & 3) + 1), TQ,Q),Q);
	av3 = avma;
	for (i = 0; i < 4; i++)
	{
	  GEN     u;
	  pj[4] = i;
	  pj[5] = (i + pj[2] - pj[1]) & 3;
	  if (DEBUGLEVEL >= 4)
	    fprintferr("S4GaloisConj:Testing %d/3:%d/2:%d/2:%d/4:%Z:%Z\n",
		       j - 1, w >> 1, l, i, sg, pj);
	  u = FpX_add(uu, FpXQ_mul(gel(pauto,(pj[4] & 3) + 1),
				   gel(bezoutcoeff,sg[1]),TQ,Q),Q);
	  u = FpX_add(u,  FpXQ_mul(gel(pauto,((-pj[4]) & 3) + 1),
				   gel(bezoutcoeff,sg[3]), TQ,Q),Q);
	  u = FpX_add(u,  FpXQ_mul(gel(pauto,(pj[5] & 3) + 1),
				   gel(bezoutcoeff,sg[2]), TQ,Q),Q);
	  u = FpX_add(u,  FpXQ_mul(gel(pauto,((-pj[5]) & 3) + 1),
				   gel(bezoutcoeff,sg[4]), TQ,Q),Q);
	  if (s4test(u, liftpow, gl, tau))
	    goto suites4_2;
	  avma = av3;
	}
	avma = av2;
      }
      z = sg[4];
      sg[4] = sg[3];
      sg[3] = z;
      pj[2] = (-pj[2]) & 3;
      avma = av;
    }
  }
  avma = ltop;
  return gen_0;
suites4_2:
  avma = av;
  {
    long abc, abcdef;
    GEN     u;
    pari_sp av2;
    abc = (pj[1] + pj[2] + pj[3]) & 3;
    abcdef = (((abc + pj[4] + pj[5] - pj[6]) & 3) >> 1);
    u = s4releveauto(misom,Tmod,Tp,p,sg[1],sg[4],sg[2],sg[5],sg[3],sg[6]);
    s4makelift(u, gl, liftpow);
    av2 = avma;
    for (j = 0; j < 8; j++)
    {
      long h, g, i;
      h = j & 3;
      g = abcdef + ((j & 4) >> 1);
      i = h + abc - g;
      u = FpXQ_mul(gel(pauto,(g & 3) + 1),
		   gel(bezoutcoeff,sg[1]),TQ,Q);
      u = ZX_add(u, FpXQ_mul(gel(pauto,((-g) & 3) + 1),
			     gel(bezoutcoeff,sg[4]),TQ,Q));
      u = ZX_add(u, FpXQ_mul(gel(pauto,(h & 3) + 1),
			     gel(bezoutcoeff,sg[2]),TQ,Q));
      u = ZX_add(u, FpXQ_mul(gel(pauto,((-h) & 3) + 1),
			     gel(bezoutcoeff,sg[5]), TQ,Q));
      u = ZX_add(u, FpXQ_mul(gel(pauto,(i & 3) + 1),
			     gel(bezoutcoeff,sg[3]), TQ,Q));
      u = FpX_add(u, FpXQ_mul(gel(pauto,((-i) & 3) + 1),
			      gel(bezoutcoeff,sg[6]), TQ,Q),Q);
      if (DEBUGLEVEL >= 4)
	fprintferr("S4GaloisConj:Testing %d/8 %d:%d:%d\n",
		   j, g & 3, h & 3, i & 3);
      if (s4test(u, liftpow, gl, phi))
	break;
      avma = av2;
    }
  }
  if (j == 8)
  {
    avma = ltop;
    return gen_0;
  }
  for (i = 1; i < lg(gl->L); i++)
  {
    mael3(res,1,1,i) = sigma[tau[i]];
    mael3(res,1,2,i) = phi[sigma[tau[phi[i]]]];
    mael3(res,1,3,i) = phi[sigma[i]];
    mael3(res,1,4,i) = sigma[i];
  }
  avma = ltop2;
  return res;
}
struct galois_frobenius
{
  long p;
  long fp;
  long deg;
  GEN Tmod;
  GEN psi;
};

/*Warning : the output of this function is not gerepileupto
 * compatible...*/
static GEN
galoisfindgroups(GEN lo, GEN sg, long f)
{
  pari_sp ltop=avma;
  GEN V,W;
  long i,j,k;
  V=cgetg(lg(lo),t_VEC);
  for(j=1,i=1;i<lg(lo);i++)
  {
    pari_sp av=avma;
    GEN U;
    W=cgetg(lg(lo[i]),t_VECSMALL);
    for(k=1;k<lg(lo[i]);k++)
      W[k]=mael(lo,i,k)%f;
    vecsmall_sort(W); 
    U=vecsmall_uniq(W);
    if (gequal(U, sg))
    {
      cgiv(U);
      V[j++]=lo[i];
    }
    else
      avma=av;
  }
  setlg(V,j);
  /*warning components of V point to W*/
  return gerepileupto(ltop,V);
}

static long
galoisfrobeniustest(GEN aut, struct galois_lift *gl, GEN frob)
{
  pari_sp ltop = avma;
  GEN tlift = FpX_center(FpX_Fp_mul(aut,gl->den,gl->Q), gl->Q);
  long res = poltopermtest(tlift, gl, frob);
  avma = ltop;
  return res;
}

static GEN
galoismakepsi(long g, GEN sg, GEN pf)
{
  GEN psi=cgetg(g+1,t_VECSMALL);
  long i;
  for (i = 1; i < g; i++)
    psi[i] = sg[pf[i]];
  psi[g]=sg[1];
  return psi;
}

static GEN
galoisfrobeniuslift(GEN T, GEN den, GEN L,  GEN Lden, 
    struct galois_frobenius *gf,  struct galois_borne *gb) 
{
  pari_sp ltop=avma, av2;
  struct galois_testlift gt;
  struct galois_lift gl;
  GEN res;
  long i,j,k;
  long n=lg(L)-1, deg=1, g=lg(gf->Tmod)-1;
  GEN F,Fp,Fe;
  GEN ip = utoipos(gf->p), aut, frob;
  if (DEBUGLEVEL >= 4)
    fprintferr("GaloisConj:p=%ld deg=%ld fp=%ld\n", gf->p, deg, gf->fp);
  res = cgetg(lg(L), t_VECSMALL);
  gf->psi = const_vecsmall(g,1);
  av2=avma;
  initlift(T, den, ip, L, Lden, gb, &gl);
  aut = galoisdolift(&gl, res);
  if (!aut || galoisfrobeniustest(aut,&gl,res))
  {
    avma=av2;
    gf->deg = gf->fp;
    return res;
  }
  inittestlift(aut,gf->Tmod, &gl, &gt);
  gt.C=cgetg(gf->fp+1,t_VEC);
  for (i = 1; i <= gf->fp; i++)
  {
    gel(gt.C,i) = cgetg(gt.g+1,t_VECSMALL);
    for(j = 1; j <= gt.g; j++) mael(gt.C,i,j) = 0;
  }
  gt.Cd=gcopy(gt.C);

  F =factoru(gf->fp);
  Fp=gel(F,1);
  Fe=gel(F,2);
  frob = cgetg(lg(L), t_VECSMALL);
  for(k=lg(Fp)-1;k>=1;k--)
  {
    pari_sp btop=avma;
    GEN psi=NULL,fres=NULL,sg;
    long el=gf->fp, dg=1, dgf=1;
    long e,pr;
    sg=perm_identity(1);
    for(e=1;e<=Fe[k];e++)
    {
      long l;
      GEN lo;
      GEN pf;
      dg *= Fp[k]; el /= Fp[k];
      if ( DEBUGLEVEL>=4 )
	fprintferr("Trying degre %d.\n",dg);
      if (galoisfrobeniustest(gel(gt.pauto,el+1),&gl,frob))
      {
	dgf = dg; 
	psi = const_vecsmall(g,1);
	fres= gcopy(frob);
	continue;
      }
      disable_dbg(0);
      lo = listznstarelts(dg, n / gf->fp);
      disable_dbg(-1);
      if (e!=1)
	lo = galoisfindgroups(lo, sg, dgf);
      if (DEBUGLEVEL >= 4)
	fprintferr("Galoisconj:Subgroups list:%Z\n", lo);
      for (l = 1; l < lg(lo); l++)
	if ( lg(lo[l])>2 && 
	    frobeniusliftall(gel(lo,l), el, &pf, &gl, &gt, frob))
	{
	  sg  = gcopy(gel(lo,l));
	  psi = galoismakepsi(g,sg,pf);
	  dgf = dg;
	  fres=gcopy(frob);
	  break;
	}
      if ( l == lg(lo) )
	break;
    }
    if (dgf==1) { avma=btop; continue; }
    pr=deg*dgf;
    if (deg==1)
    {
      for(i=1;i<lg(res);i++) res[i]=fres[i];
      for(i=1;i<lg(psi);i++) gf->psi[i]=psi[i];
    } 
    else
    {
      GEN cp=perm_mul(res,fres);
      for(i=1;i<lg(res);i++) res[i]=cp[i];
      for(i=1;i<lg(psi);i++) gf->psi[i]=(dgf*gf->psi[i]+deg*psi[i])%pr;
    }
    deg=pr;
    avma=btop;
  }
  for (i = 1; i <= gf->fp; i++)
    for (j = 1; j <= gt.g; j++)
      if (mael(gt.C,i,j)) gunclone(gmael(gt.C,i,j));
  if (DEBUGLEVEL>=4 && res)
    fprintferr("Best lift: %d\n",deg);
  if (deg==1)
  {
    avma=ltop;
    return NULL;
  }
  else
  {
    /*We need to normalise result so that psi[g]=1*/
    long im=Fl_inv(gf->psi[g],deg);
    GEN cp=perm_pow(res, im);
    for(i=1;i<lg(res);i++) res[i]=cp[i];
    for(i=1;i<lg(gf->psi);i++) gf->psi[i] = Fl_mul(im,gf->psi[i],deg);
    avma=av2;
    gf->deg=deg;
    return res;
  }
}

static GEN
galoisfindfrobenius(GEN T, GEN L, GEN den, struct galois_frobenius *gf,
    struct galois_borne *gb, const struct galois_analysis *ga)
{
  pari_sp lbot, ltop=avma;
  long Try=0;
  long n = degpol(T), deg, gmask;
  byteptr primepointer = ga->primepointer;
  GEN Lden,frob;
  Lden=makeLden(L,den,gb);
  gf->deg=ga->deg;gf->p=ga->p; deg=ga->deg;
  gmask=(ga->group&ga_ext_2)?3:1;
  for (;;)
  {
    pari_sp av = avma;
    long    isram;
    long    i;
    GEN     ip,Tmod;
    ip = utoipos(gf->p);
    Tmod = lift_intern(factmod(T, ip));
    isram = 0;
    for (i = 1; i < lg(Tmod[2]) && !isram; i++)
      if (!gcmp1(gmael(Tmod,2,i)))
	isram = 1;
    if (isram == 0)
    {
      gf->fp = degpol(gmael(Tmod,1,1));
      for (i = 2; i < lg(Tmod[1]); i++)
	if (degpol(gmael(Tmod,1,i)) != gf->fp)
	{
	  avma = ltop;
	  return NULL;		/* Not Galois polynomial */
	}
      lbot=avma;
      gf->Tmod=gcopy(gel(Tmod,1));
      if ( ((gmask&1) && gf->fp % deg == 0) || ((gmask&2) && gf->fp % 2== 0) )
      {
	frob=galoisfrobeniuslift(T, den, L, Lden, gf, gb);
	if (frob)
	{
	  GEN *gptr[3];
	  gptr[0]=&gf->Tmod;
	  gptr[1]=&gf->psi;
	  gptr[2]=&frob;
	  gerepilemanysp(ltop,lbot,gptr,3);
	  return frob;
	}
	if ((ga->group&ga_all_normal) && gf->fp % deg == 0)
	  gmask&=~1;
	/*The first prime degree is always divisible by deg, so we don't
	 * have to worry about ext_2 being used before regular supersolvable*/
	if (!gmask)
	{
	  avma = ltop;
	  return NULL;
	}
	Try++;
	if ( (ga->group&ga_non_wss) && Try > n )
	  pari_warn(warner, "galoisconj _may_ hang up for this polynomial");
      }
    }
    NEXT_PRIME_VIADIFF_CHECK(gf->p, primepointer);
    if (DEBUGLEVEL >= 4)
      fprintferr("GaloisConj:next p=%ld\n", gf->p);
    avma = av;
  }
}

static GEN
galoisgen(GEN T, GEN L, GEN M, GEN den, struct galois_borne *gb,
	  const struct galois_analysis *ga);
static GEN
galoisgenfixedfield(GEN Tp, GEN Pmod, GEN V, GEN ip, struct galois_borne *gb, GEN Pg)
{
  pari_sp ltop=avma;
  GEN     P, PL, Pden, PM, Pp, Pladicabs;
  GEN     tau, PG;
  long    g,gp;
  long    x=varn(Tp);
  P=gel(V,3);
  PL=gel(V,2);
  gp=lg(Pmod)-1;
  Pp = FpX_red(P,ip);
  if (DEBUGLEVEL>=6)
    fprintferr("GaloisConj: Fixed field %Z\n",P);
  if (degpol(P)==2)
  {
    PG=cgetg(3,t_VEC);
    gel(PG,1) = mkvec( mkvecsmall2(2,1) );
    gel(PG,2) = mkvecsmall(2);
    tau = deg1pol_i(gen_m1, negi(gel(P,3)), x);
    tau = RgX_to_FpX(tau, ip);
    tau = FpX_FpXQ_compo(gel(Pmod,gp), tau,Pp,ip);
    tau = FpX_gcd(Pp, tau,ip);
    tau = FpX_normalize(tau, ip);
    for (g = 1; g <= gp; g++)
      if (gequal(tau, gel(Pmod,g)))
	break;
    if (g == lg(Pmod))
      return NULL;
    Pg[1]=g;
  }
  else
  {
    struct galois_analysis Pga;
    struct galois_borne Pgb;
    long j;
    galoisanalysis(P, &Pga, 0);
    if (Pga.deg == 0)
      return NULL;		/* Avoid computing the discriminant */
    Pgb.l = gb->l;
    Pden = galoisborne(P, NULL, &Pgb);
    Pladicabs=Pgb.ladicabs;
    if (Pgb.valabs > gb->valabs)
    {
      if (DEBUGLEVEL>=4)
	fprintferr("GaloisConj:increase prec of p-adic roots of %ld.\n"
	    ,Pgb.valabs-gb->valabs);
      PL = ZpX_liftroots(P,PL,gb->l,Pgb.valabs);
    }
    else if (Pgb.valabs < gb->valabs)
      PL = FpC_red(PL, Pgb.ladicabs);
    PM = vandermondeinversemod(PL, P, Pden, Pgb.ladicabs);
    PG = galoisgen(P, PL, PM, Pden, &Pgb, &Pga);
    if (PG == gen_0) return NULL;
    for (j = 1; j < lg(PG[1]); j++)
    {
      pari_sp btop=avma;
      tau = permtopol(gmael(PG,1,j), PL, PM, Pden, Pladicabs, x);
      tau = RgX_to_FpX(tau, ip);
      tau = FpX_FpXQ_compo(gel(Pmod,gp), tau,Pp,ip);
      tau = FpX_gcd(Pp, tau,ip);
      tau = FpX_normalize(tau, ip);
      for (g = 1; g < lg(Pmod); g++)
	if (gequal(tau, gel(Pmod,g)))
	  break;
      if (g == lg(Pmod))
	return NULL;
      avma=btop;
      Pg[j]=g;
    }
  }
  return gerepilecopy(ltop,PG);
}

/* Let 
 * sigma^m=1
 * tau*sigma*tau^-1=sigma^s.
 * Compute n so that 
 * (sigma*tau)^e=sigma^n*tau^e
 * We have n=sum_{k=0}^{e-1} s^k mod m.
 * so n*(1-s) = 1-s^e mod m
 * Unfortunately (1-s) might not invertible mod m.
 */

static long 
stpow(long s, long e, long m)
{
  long i;
  long n = 1;
  for (i = 1; i < e; i++)
    n = (1 + n * s) % m;
  return n;
}

static GEN
wpow(long s, long m, long e, long n)
{
  GEN   w = cgetg(n+1,t_VECSMALL);
  long si = s;
  long i;
  w[1] = 1;
  for(i=2; i<=n; i++)
    w[i] = w[i-1]*e;
  for(i=n; i>=1; i--)
  {
    si = Fl_pow(si,e,m);
    w[i] = Fl_mul(s-1, stpow(si, w[i], m), m);
  }
  return w;
}

static GEN
galoisgen(GEN T, GEN L, GEN M, GEN den, struct galois_borne *gb,
	  const struct galois_analysis *ga)
{
  struct galois_test td;
  struct galois_frobenius gf;
  pari_sp lbot, ltop2, ltop = avma;
  long    n, p, deg, x;
  long    i, j;
  GEN     Lden, sigma;
  GEN     Tmod, res, pf = gen_0, ip;
  GEN     frob;
  GEN     O;
  GEN     PG, Pg;
  n = degpol(T);
  if (!ga->deg)
    return gen_0;
  x = varn(T);
  if (DEBUGLEVEL >= 9)
    fprintferr("GaloisConj:denominator:%Z\n", den);
  if (n == 12 && ga->ord==3)	/* A4 is very probable,so test it first */
  {
    pari_sp av = avma;
    if (DEBUGLEVEL >= 4)
      fprintferr("GaloisConj:Testing A4 first\n");
    inittest(L, M, gb->bornesol, gb->ladicsol, &td);
    lbot = avma;
    PG = a4galoisgen(T, &td);
    freetest(&td);
    if (PG != gen_0)
      return gerepile(ltop, lbot, PG);
    avma = av;
  }
  if (n == 24 && ga->ord==3)	/* S4 is very probable,so test it first */
  {
    pari_sp av = avma;
    struct galois_lift gl;
    if (DEBUGLEVEL >= 4)
      fprintferr("GaloisConj:Testing S4 first\n");
    lbot = avma;
    Lden=makeLden(L,den,gb);
    initlift(T, den, stoi(ga->p4), L, Lden, gb, &gl);
    PG = s4galoisgen(&gl);
    if (PG != gen_0)
      return gerepile(ltop, lbot, PG);
    avma = av;
  }
  frob=galoisfindfrobenius(T, L, den, &gf, gb, ga);
  if (!frob)
  {
    ltop=avma;
    return gen_0;
  }
  p=gf.p;
  ip = utoipos(p);
  Tmod=gf.Tmod;
  O = perm_cycles(frob);
  deg=lg(O[1])-1;
  sigma = permtopol(frob, L, M, den, gb->ladicabs, x);
  if (DEBUGLEVEL >= 9)
    fprintferr("GaloisConj:Orbite:%Z\n", O);
  if (deg == n)			/* Cyclique */
  {
    lbot = avma;
    res = cgetg(3, t_VEC);
    gel(res,1) = mkvec( cyc_pow_perm(O,1) );
    gel(res,2) = mkvecsmall(deg);
    return gerepile(ltop, lbot, res);
  }
  if (DEBUGLEVEL >= 9)
    fprintferr("GaloisConj:Frobenius:%Z\n", sigma);
  Pg=cgetg(lg(O),t_VECSMALL);
  {
    pari_sp btop=avma;
    GEN     V, Tp, Sp, Pmod;
    GEN OL = fixedfieldorbits(O,L);
    V  = fixedfieldsympol(OL, gb->ladicabs, gb->l, ip, x);
    Tp = FpX_red(T,ip);
    Sp = sympol_aut_evalmod(gel(V,1),deg,sigma,Tp,ip);
    Pmod = fixedfieldfactmod(Sp,ip,Tmod);
    PG=galoisgenfixedfield(Tp, Pmod, V, ip, gb, Pg);
    if (PG == NULL)
    {
      avma = ltop;
      return gen_0;
    }
    if (DEBUGLEVEL >= 4)
      fprintferr("GaloisConj:Back to Earth:%Z\n", PG);
    PG=gerepileupto(btop, PG);
  }
  inittest(L, M, gb->bornesol, gb->ladicsol, &td);
  lbot = avma;
  res = cgetg(3, t_VEC);
  gel(res,1) = cgetg(lg(PG[1]) + 1, t_VEC);
  gel(res,2) = cgetg(lg(PG[1]) + 1, t_VECSMALL);
  gmael(res,1,1) = cyc_pow_perm(O,1);
  mael(res,2,1) = deg;
  for (i = 2; i < lg(res[1]); i++)
    gmael(res,1,i) = cgetg(n + 1, t_VECSMALL);
  ltop2 = avma;
  for (j = 1; j < lg(PG[1]); j++)
  {
    long sr;
    long k;
    GEN  X  = cgetg(lg(O), t_VECSMALL);
    GEN  oX = cgetg(lg(O), t_VECSMALL);
    GEN  gj = gmael(PG,1,j);
    long s  = gf.psi[Pg[j]];
    GEN  B  = perm_cycles(gj);
    long oj = lg(B[1]) - 1;
    GEN  F  = factoru_pow(oj);
    GEN  Fp = gel(F,1);
    GEN  Fe = gel(F,2);
    GEN  Fc = gel(F,3);
    if (DEBUGLEVEL >= 6)
      fprintferr("GaloisConj: G[%d]=%Z of relative order %d\n", j, gj, oj);
    B = perm_cycles(gmael(PG,1,j));
    pf = perm_identity(n);
    for (k=lg(Fp)-1; k>=1; k--)
    {
      long f;
      long dg = 1, el = oj;
      long osel = 1;
      long p  = Fp[k];
      long e  = Fe[k];
      long op = oj/Fc[k];
      long a=0;
      GEN  pf1 = NULL, w, wg;
      GEN  Be = cgetg(e+1,t_VEC);
      gel(Be,e) = cyc_pow(B, op);
      for(i=e-1; i>=1; i--)
        gel(Be,i) = cyc_pow(gel(Be,i+1), p);
      w = wpow(Fl_pow(s,op,deg),deg,p,e);
      wg = cgetg(e+2,t_VECSMALL);
      wg[e+1] = deg;
      for (i=e; i>=1; i--)
        wg[i]=cgcd(wg[i+1], w[i]);
      for (i=1; i<lg(O); i++) oX[i] = 0;
      for (f=1; f<=e; f++)
      { 
        long sel;
        GEN Bel = gel(Be,f);
        long t;
        dg *= p; el /= p;
        sel = Fl_pow(s,el,deg); 
        if (DEBUGLEVEL >= 6)
          fprintferr("GaloisConj: B=%Z\n", Bel);
        sr  = cgcd(stpow(sel,p,deg),deg);
        if (DEBUGLEVEL >= 6)
          fprintferr("GaloisConj: exp %d: s=%ld [%ld] a=%ld w=%ld wg=%ld sr=%ld\n",
              dg, sel, deg, a, w[f], wg[f+1], sr);
        for (t = 0; t < sr; t++)
          if ((a+t*w[f])%wg[f+1]==0)
          {
            long i, j, k, st;
            for (i = 1; i < lg(X); i++) X[i] = 0;
            for (i = 0; i < lg(X); i+=dg)
              for (j = 1, k = p, st = t; k <= dg; j++, k += p)
              {
                X[k+i] = (oX[j+i] + st)%deg;
                st = (t + st*osel)%deg;
              }
            pf1 = testpermutation(O, Bel, X, sel, p, sr, &td);
            if (pf1)
              break;
          }
        if (!pf1) { freetest(&td); avma = ltop; return gen_0; }
        for (i=1; i<lg(O); i++) oX[i] = X[i];
        osel = sel; a = (a+t*w[f])%deg;
      }
      pf = perm_mul(pf, perm_pow(pf1, el));
    }
    for (k = 1; k <= n; k++) gmael3(res,1,j + 1,k) = gel(pf,k);
    gmael(res,2,j+1) = gmael(PG,2,j);
    avma = ltop2;
  }
  if (DEBUGLEVEL >= 4)
    fprintferr("GaloisConj:Fini!\n");
  freetest(&td);
  return gerepile(ltop, lbot, res);
}

/* T: polynomial or nf, den multiple of common denominator of solutions or
 * NULL (unknown). If T is nf, and den unknown, use den = denom(nf.zk) */
GEN
galoisconj4(GEN T, GEN den, long flag)
{
  pari_sp ltop = avma;
  GEN     G, L, M, res, aut, grp=NULL;/*keep gcc happy on the wall*/
  struct galois_analysis ga;
  struct galois_borne gb;
  long    n, i, j, k;
  if (typ(T) != t_POL)
  {
    T = checknf(T);
    if (!den) den = Q_denom(gel(T,7));
    T = gel(T,1);
  }
  n = degpol(T);
  if (n <= 0)
    pari_err(constpoler, "galoisconj4");
  for (k = 2; k <= n + 2; k++)
    if (typ(T[k]) != t_INT)
      pari_err(talker, "polynomial not in Z[X] in galoisconj4");
  if (!gcmp1(gel(T,n + 2)))
    pari_err(talker, "non-monic polynomial in galoisconj4");
  n = degpol(T);
  if (n == 1)			/* Too easy! */
  {
    if (!flag) return mkcol( pol_x[varn(T)] );
    ga.l = 3;
    ga.deg = 1;
    den = gen_1;
  }
  else
    galoisanalysis(T, &ga, 1);
  if (ga.deg == 0)
  {
    avma = ltop;
    return utoipos(ga.p); /* Avoid computing the discriminant */
  }
  if (den)
  {
    if (typ(den) != t_INT)
      pari_err(talker, "Second arg. must be integer in galoisconj4");
    den = absi(den);
  }
  gb.l = utoipos(ga.l);
  if (DEBUGLEVEL >= 1) (void)timer2();
  den = galoisborne(T, den, &gb);
  if (DEBUGLEVEL >= 1)
    msgtimer("galoisborne()");
  L = rootpadicfast(T, gb.l, gb.valabs);
  if (DEBUGLEVEL >= 1)
    msgtimer("rootpadicfast()");
  M = vandermondeinversemod(L, T, den, gb.ladicabs);
  if (DEBUGLEVEL >= 1)
    msgtimer("vandermondeinversemod()");
  if (n == 1)
  {
    G = cgetg(3, t_VEC);
    gel(G,1) = cgetg(1, t_VECSMALL);
    gel(G,2) = cgetg(1, t_VECSMALL);
  }
  else
    G = galoisgen(T, L, M, den, &gb, &ga);
  if (DEBUGLEVEL >= 6)
    fprintferr("GaloisConj:%Z\n", G);
  if (G == gen_0)
  {
    avma = ltop;
    return gen_0;
  }
  if (DEBUGLEVEL >= 1) (void)timer2();
  if (flag)
  {
    grp = cgetg(9, t_VEC);
    gel(grp,1) = gcopy(T);
    gel(grp,2) = cgetg(4,t_VEC); /*Make K.B. happy(8 components)*/
    gmael(grp,2,1) = stoi(ga.l);
    gmael(grp,2,2) = stoi(gb.valabs);
    gmael(grp,2,3) = icopy(gb.ladicabs);
    gel(grp,3) = gcopy(L);
    gel(grp,4) = gcopy(M);
    gel(grp,5) = gcopy(den);
    gel(grp,7) = gcopy(gel(G,1));
    gel(grp,8) = gcopy(gel(G,2));
  }
  res = cgetg(n + 1, t_VEC);
  gel(res,1) = perm_identity(n);
  k = 1;
  for (i = 1; i < lg(G[1]); i++)
  {
    long  c = k * (mael(G,2,i) - 1);
    for (j = 1; j <= c; j++)	/* I like it */
      gel(res,++k) = perm_mul(gel(res,j), gmael(G,1,i));
  }
  if (flag)
  {
    gel(grp,6) = res;
    return gerepileupto(ltop, grp);
  }
  aut = galoisgrouptopol(res,L,M,den,gb.ladicsol, varn(T));
  if (DEBUGLEVEL >= 1)
    msgtimer("Calcul polynomes");
  return gerepileupto(ltop, gen_sort(aut, 0, cmp_pol));
}

/* Heuristic computation of #Aut(T), pdepart first prime to be tested */
long
numberofconjugates(GEN T, long pdepart)
{
  pari_sp ltop2, ltop = avma;
  long    n, p, nbmax, nbtest;
  long    card;
  byteptr primepointer;
  long    i;
  GEN     L;
  n = degpol(T);
  card = sturm(T);
  card = cgcd(card, n - card);
  nbmax = (n<<1) + 1;
  if (nbmax < 20) nbmax=20;
  nbtest = 0;
  L = cgetg(n + 1, t_VECSMALL);
  ltop2 = avma;
  for (p = 0, primepointer = diffptr; nbtest < nbmax && card > 1;)
  {
    long    s;
    long    isram;
    GEN     S;

    NEXT_PRIME_VIADIFF_CHECK(p,primepointer);
    if (p < pdepart) continue;
    S = FpX_degfact(T, utoipos(p));
    isram = 0;
    for (i = 1; i < lg(S[2]) && !isram; i++)
      if (mael(S,2,i) != 1) isram = 1;
    if (!isram)
    {
      for (i = 1; i <= n; i++) L[i] = 0;
      for (i = 1; i < lg(S[1]) && !isram; i++) L[ mael(S,1,i) ]++;
      s = L[1];
      for (i = 2; i <= n; i++) s = cgcd(s, L[i] * i);
      card = cgcd(s, card);
    }
    if (DEBUGLEVEL >= 6)
      fprintferr("NumberOfConjugates:Nbtest=%ld,card=%ld,p=%ld\n", nbtest,
		 card, p);
    nbtest++;
    avma = ltop2;
  }
  if (DEBUGLEVEL >= 2)
    fprintferr("NumberOfConjugates:card=%ld,p=%ld\n", card, p);
  avma = ltop;
  return card;
}

GEN
galoisconj0(GEN nf, long flag, GEN d, long prec)
{
  pari_sp ltop;
  GEN     G, T;
  long    card;
  if (typ(nf) != t_POL)
  {
    nf = checknf(nf);
    T = gel(nf,1);
  }
  else
    T = nf;
  switch (flag)
  {
  case 0:
    ltop = avma;
    G = galoisconj4(nf, d, 0);
    if (typ(G) != t_INT)	/* Success */
      return G;
    else
    {
      card = numberofconjugates(T, G == gen_0 ? 2 : itos(G));
      avma = ltop;
      if (card != 1)
      {
	if (typ(nf) == t_POL)
	{
	  G = galoisconj2pol(nf, card, prec);
	  if (lg(G) <= card)
	    pari_warn(warner, "conjugates list may be incomplete in nfgaloisconj");
	  return G;
	}
	else
	  return galoisconj(nf);
      }
    }
    break;			/* Failure */
  case 1:
    return galoisconj(nf);
  case 2:
    return galoisconj2(nf, degpol(T), prec);
  case 4:
    G = galoisconj4(nf, d, 0);
    if (typ(G) != t_INT) return G;
    break;			/* Failure */
  default:
    pari_err(flagerr, "nfgaloisconj");
  }
  return mkcol( pol_x[varn(T)] );	/* Failure */
}



/******************************************************************************/
/* Isomorphism between number fields                                          */
/******************************************************************************/
long
isomborne(GEN P, GEN den, GEN p)
{
  pari_sp ltop=avma;
  struct galois_borne gb;
  gb.l=p;
  (void)galoisborne(P,den,&gb);
  avma=ltop;
  return gb.valsol;
}



/******************************************************************************/
/* Galois theory related algorithms                                           */
/******************************************************************************/
GEN
checkgal(GEN gal)
{
  if (typ(gal) == t_POL) pari_err(talker, "please apply galoisinit first");
  if (typ(gal) != t_VEC || lg(gal) != 9)
    pari_err(talker, "Not a Galois field in a Galois related function");
  return gal;
}

GEN
galoisinit(GEN nf, GEN den)
{
  GEN G = galoisconj4(nf, den, 1);
  if (typ(G) == t_INT)
    pari_err(talker, "field not Galois or Galois group not weakly super solvable");
  return G;
}

GEN
galoispermtopol(GEN gal, GEN perm)
{
  GEN     v;
  long    t = typ(perm), i;
  gal = checkgal(gal);
  switch (t)
  {
  case t_VECSMALL:
    return permtopol(perm, gel(gal,3), gel(gal,4), gel(gal,5),
		     gmael(gal,2,3), varn(gel(gal,1)));
  case t_VEC:
  case t_COL:
  case t_MAT:
    v = cgetg(lg(perm), t);
    for (i = 1; i < lg(v); i++)
      gel(v,i) = galoispermtopol(gal, gel(perm,i));
    return v;
  }
  pari_err(typeer, "galoispermtopol");
  return NULL;			/* not reached */
}


GEN 
galoiscosets(GEN O, GEN perm)
{
  long i,j,k,u;
  long o = lg(O)-1, f = lg(O[1])-1;
  GEN C = cgetg(lg(O),t_VECSMALL);
  pari_sp av=avma;
  GEN RC=cgetg(lg(perm),t_VECSMALL);
  for(i=1;i<lg(RC);i++)
    RC[i]=0;
  u=mael(O,1,1);
  for(i=1,j=1;j<=o;i++)
  {
    if (RC[mael(perm,i,u)])
      continue;
    for(k=1;k<=f;k++)
      RC[mael(perm,i,mael(O,1,k))]=1;
    C[j++]=i;
  }
  avma=av;
  return C;
}

GEN
fixedfieldfactor(GEN L, GEN O, GEN perm, GEN M, GEN den, GEN mod,
                 long x,long y)
{
  pari_sp ltop=avma;
  GEN     F,G,V,res,cosets;
  long    i, j, k;
  F=cgetg(lg(O[1])+1,t_COL);
  gel(F, lg(O[1])) = gen_1;
  G=cgetg(lg(O),t_VEC);
  for (i = 1; i < lg(O); i++)
  {
    GEN Li = cgetg(lg(O[i]),t_VEC);
    for (j = 1; j < lg(O[i]); j++)
      Li[j] = L[mael(O,i,j)];
    gel(G,i) = FpV_roots_to_pol(Li,mod,x);
  }  
  
  cosets=galoiscosets(O,perm);
  if (DEBUGLEVEL>=4) fprintferr("GaloisFixedField:cosets=%Z \n",cosets);
  V=cgetg(lg(O),t_COL);
  if (DEBUGLEVEL>=6) fprintferr("GaloisFixedField:den=%Z mod=%Z \n",den,mod);
  res=cgetg(lg(O),t_VEC);
  for (i = 1; i < lg(O); i++)
  {
    pari_sp av=avma;
    long ii,jj;
    GEN G=cgetg(lg(O),t_VEC);
    for (ii = 1; ii < lg(O); ii++)
    {
      GEN Li = cgetg(lg(O[ii]),t_VEC);
      for (jj = 1; jj < lg(O[ii]); jj++)
	Li[jj] = L[mael(perm,cosets[i],mael(O,ii,jj))];
      gel(G,ii) = FpV_roots_to_pol(Li,mod,x);
    }  
    for (j = 1; j < lg(O[1]); j++)
    {
      for(k = 1; k < lg(O); k++)
	V[k]=mael(G,k,j+1);
      gel(F,j) = vectopol(V, M, den, mod, y);
    }
    gel(res,i) = gerepileupto(av,gtopolyrev(F,x));
  }
  return gerepileupto(ltop,res);
}

GEN
galoisfixedfield(GEN gal, GEN perm, long flag, long y)
{
  pari_sp lbot, ltop = avma;
  GEN     L, P, S, PL, O, res, mod;
  long    x, n, i;
  gal = checkgal(gal);
  x = varn(gel(gal,1));
  L = gel(gal,3); n=lg(L)-1;
  mod = gmael(gal,2,3);
  if (flag<0 || flag>2)
    pari_err(flagerr, "galoisfixedfield");
  if (typ(perm) == t_VEC)
  {
    for (i = 1; i < lg(perm); i++)
      if (typ(perm[i]) != t_VECSMALL || lg(perm[i])!=n+1)
        pari_err(typeer, "galoisfixedfield");
    O = vecperm_orbits(perm, n);
  }
  else if (typ(perm) != t_VECSMALL || lg(perm)!=n+1 )
  {
    pari_err(typeer, "galoisfixedfield");
    return NULL; /* not reached */
  }
  else
    O = perm_cycles(perm);
  {
    GEN OL= fixedfieldorbits(O,L);
    GEN V = fixedfieldsympol(OL, mod, gmael(gal,2,1), NULL, x);
    P=gel(V,3);
    PL=gel(V,2);
  }
  if (flag==1)
    return gerepileupto(ltop,P);
  S = fixedfieldinclusion(O, PL);
  S = vectopol(S, gel(gal,4), gel(gal,5), mod, x);
  if (flag==0)
  {
    lbot = avma;
    res = cgetg(3, t_VEC);
    gel(res,1) = gcopy(P);
    gel(res,2) = gmodulo(S, gel(gal,1));
    return gerepile(ltop, lbot, res);
  }
  else
  {
    GEN PM,Pden;
    {
      struct galois_borne Pgb;
      long val=itos(gmael(gal,2,2));
      Pgb.l = gmael(gal,2,1);
      Pden = galoisborne(P, NULL, &Pgb);
      if (Pgb.valabs > val)
      {
        if (DEBUGLEVEL>=4)
          fprintferr("GaloisConj:increase prec of p-adic roots of %ld.\n"
              ,Pgb.valabs-val);
        PL = ZpX_liftroots(P,PL,Pgb.l,Pgb.valabs);
        L = ZpX_liftroots(gel(gal,1),L,Pgb.l,Pgb.valabs);
        mod = Pgb.ladicabs;
      }
    }
    PM = vandermondeinversemod(PL, P, Pden, mod);
    lbot = avma;
    if (y==-1)
      y = fetch_user_var("y");
    if (y<=x)
      pari_err(talker,"priority of optional variable too high in galoisfixedfield");
    res = cgetg(4, t_VEC);
    gel(res,1) = gcopy(P);
    gel(res,2) = gmodulo(S, gel(gal,1));
    gel(res,3) = fixedfieldfactor(L,O,gel(gal,6),
        PM,Pden,mod,x,y);
    return gerepile(ltop, lbot, res);
  }
}
/* gal being a galois group output the underlying wss group.
 */

GEN
galois_group(GEN gal) { return mkvec2(gel(gal,7), gel(gal,8)); }

GEN
checkgroup(GEN g, GEN *S)
{
  if (typ(g)==t_VEC && lg(g)==3 && typ(g[1])==t_VEC && typ(g[2])==t_VECSMALL) 
  {
    *S = NULL;
    return g;
  }
  g  = checkgal(g); 
  *S = gel(g,6);
  return galois_group(g); 
}

GEN
galoisisabelian(GEN gal, long flag)
{
  pari_sp ltop = avma;
  GEN S, G = checkgroup(gal,&S);
  if (!group_isabelian(G)) {avma=ltop;return gen_0;}
  if (flag==1) {avma=ltop;return gen_1;}
  if (flag==2) return gerepileupto(ltop,group_abelianSNF(G,S));
  if (flag) pari_err(flagerr,"galoisisabelian");
  return gerepileupto(ltop, group_abelianHNF(G,S));
}

GEN
galoissubgroups(GEN gal)
{
  pari_sp ltop=avma;
  GEN S, G = checkgroup(gal,&S);
  return gerepileupto(ltop, group_subgroups(G));
}

GEN
galoissubfields(GEN G, long flag, long v)
{
  pari_sp ltop=avma;
  long i;
  GEN L = galoissubgroups(G);
  long l2 = lg(L);
  GEN p3 = cgetg(l2, t_VEC);
  for (i = 1; i < l2; ++i)
    gel(p3,i) = galoisfixedfield(G, gmael(L,i,1), flag, v);
  return gerepileupto(ltop,p3);
}

GEN
galoisexport(GEN gal, long format)
{
  pari_sp ltop = avma;
  GEN S, G = checkgroup(gal,&S);
  return gerepileupto(ltop, group_export(G,format));
}

GEN
galoisidentify(GEN gal)
{
  pari_sp ltop=avma;
  GEN S, G = checkgroup(gal,&S);
  long idx = group_ident(G,S);
  long card = group_order(G);
  avma = ltop; return mkvec2s(card, idx);
}
