/* $Id: kummer.c 10282 2008-06-09 11:15:00Z kb $

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
/*                      KUMMER EXTENSIONS                          */
/*                                                                 */
/*******************************************************************/
#include "pari.h"
#include "paripriv.h"

typedef struct {
  GEN x;  /* tau ( Mod(x, nf.pol) ) */
  GEN zk; /* action of tau on nf.zk (as t_MAT) */
} tau_s;

typedef struct {
  GEN polnf, invexpoteta1;
  tau_s *tau;
  long m;
} toK_s;

static long
prank(GEN cyc, long ell)
{
  long i;
  for (i=1; i<lg(cyc); i++)
    if (smodis(gel(cyc,i),ell)) break;
  return i-1;
}

/* increment y, which runs through [0,d-1]^(k-1). Return 0 when done. */
static int
increment(GEN y, long k, long d)
{
  long i = k, j;
  do
  {
    if (--i == 0) return 0;
    y[i]++;
  } while (y[i] >= d);
  for (j = i+1; j < k; j++) y[j] = 0;
  return 1;
}

static int
ok_congruence(GEN X, GEN ell, long lW, GEN vecMsup)
{
  long i, l;
  if (gcmp0(X)) return 0;
  l = lg(X);
  for (i=lW; i<l; i++)
    if (gcmp0(gel(X,i))) return 0;
  l = lg(vecMsup);
  for (i=1; i<l; i++)
    if (gcmp0(FpM_FpC_mul(gel(vecMsup,i),X, ell))) return 0;
  return 1;
}

static int
ok_sign(GEN X, GEN msign, GEN arch)
{
  GEN p1 = F2V_red_ip( gmul(msign, X) );
  settyp(p1,t_VEC); return gequal(p1, arch);
}

/* REDUCTION MOD ell-TH POWERS */

static GEN
fix_be(GEN bnfz, GEN be, GEN u)
{
  GEN nf = checknf(bnfz), fu = gmael(bnfz,8,5);
  return element_mul(nf, be, factorbackelt(fu, u, nf));
}

static GEN
logarch2arch(GEN x, long r1, long prec)
{
  long i, lx = lg(x), tx = typ(x);
  GEN y = cgetg(lx, tx);

  if (tx == t_MAT)
  {
    for (i=1; i<lx; i++) gel(y,i) = logarch2arch(gel(x,i), r1, prec);
    return y;
  }
  for (i=1; i<=r1;i++) gel(y,i) = gexp(gel(x,i),prec);
  for (   ; i<lx; i++) gel(y,i) = gexp(gmul2n(gel(x,i),-1),prec);
  return y;
}

/* multiply be by ell-th powers of units as to find small L2-norm for new be */
static GEN
reducebetanaive(GEN bnfz, GEN be, GEN b, GEN ell)
{
  long i,k,n,ru,r1, prec = nfgetprec(bnfz);
  GEN z,p1,p2,nmax,c, nf = checknf(bnfz);

  r1 = nf_get_r1(nf);
  if (!b)
  {
    be = algtobasis_i(nf, be);
    b = gmul(gmael(nf,5,1), be);
  }
  n = max((itos(ell)>>1), 3);
  z = cgetg(n+1, t_VEC);
  c = gmul(real_i(gel(bnfz,3)), ell);
  c = logarch2arch(c, r1, prec); /* = embeddings of fu^ell */
  c = gprec_w(gnorm(c), DEFAULTPREC);
  b = gprec_w(gnorm(b), DEFAULTPREC); /* need little precision */
  gel(z,1) = shallowconcat(c, vecinv(c));
  for (k=2; k<=n; k++) gel(z,k) = vecmul(gel(z,1), gel(z,k-1));
  nmax = T2_from_embed_norm(b, r1);
  ru = lg(c)-1; c = zerovec(ru);
  for(;;)
  {
    GEN B = NULL;
    long besti = 0, bestk = 0;
    for (k=1; k<=n; k++)
      for (i=1; i<=ru; i++)
      {
        p1 = vecmul(b, gmael(z,k,i));    p2 = T2_from_embed_norm(p1,r1);
        if (gcmp(p2,nmax) < 0) { B=p1; nmax=p2; besti=i; bestk = k; continue; }
        p1 = vecmul(b, gmael(z,k,i+ru)); p2 = T2_from_embed_norm(p1,r1);
        if (gcmp(p2,nmax) < 0) { B=p1; nmax=p2; besti=i; bestk =-k; }
      }
    if (!B) break;
    b = B; gel(c,besti) = addis(gel(c,besti), bestk);
  }
  if (DEBUGLEVEL) fprintferr("naive reduction mod U^l: unit exp. = %Z\n",c);
  return fix_be(bnfz, be, gmul(ell,c));
}

static GEN
reduce_mod_Qell(GEN bnfz, GEN be, GEN gell)
{
  GEN c, fa;
  be = algtobasis_i(bnfz, be);
  be = primitive_part(be, &c);
  if (c)
  {
    fa = factor(c);
    gel(fa,2) = FpC_red(gel(fa,2), gell);
    c = factorback(fa, NULL);
    be = gmul(be, c);
  }
  return be;
}

/* return q, q^n r = x, v_pr(r) < n for all pr. Insist q is a genuine n-th
 * root (i.e r = 1) if strict != 0. */
GEN
idealsqrtn(GEN nf, GEN x, GEN gn, int strict)
{
  long i, l, n = itos(gn);
  GEN fa, q, Ex, Pr;

  fa = idealfactor(nf, x);
  Pr = gel(fa,1); l = lg(Pr);
  Ex = gel(fa,2); q = NULL;
  for (i=1; i<l; i++)
  {
    long ex = itos(gel(Ex,i));
    GEN e = stoi(ex / n);
    if (strict && ex % n) pari_err(talker,"not an n-th power in idealsqrtn");
    if (q) q = idealmulpowprime(nf, q, gel(Pr,i), e);
    else   q = idealpow(nf, gel(Pr,i), e);
  }
  return q? q: gen_1;
}

static GEN
reducebeta(GEN bnfz, GEN be, GEN ell)
{
  long j,ru, prec = nfgetprec(bnfz);
  GEN emb,z,u,matunit, nf = checknf(bnfz);

  if (DEBUGLEVEL>1) fprintferr("reducing beta = %Z\n",be);
  /* reduce mod Q^ell */
  be = reduce_mod_Qell(nf, be, ell);
  /* reduce l-th root */
  z = idealsqrtn(nf, be, ell, 0);
  if (typ(z) == t_MAT && !gcmp1(gcoeff(z,1,1)))
  {
    z = idealred_elt(nf, z);
    be = element_div(nf, be, element_pow(nf, z, ell));
    /* make be integral */
    be = reduce_mod_Qell(nf, be, ell);
  }
  if (DEBUGLEVEL>1) fprintferr("beta reduced via ell-th root = %Z\n",be);

  matunit = gmul(real_i(gel(bnfz,3)), ell); /* log. embeddings of fu^ell */
  for (;;)
  {
    z = get_arch_real(nf, be, &emb, prec);
    if (z) break;
    prec = (prec-1)<<1;
    if (DEBUGLEVEL) pari_warn(warnprec,"reducebeta",prec);
    nf = nfnewprec(nf,prec);
  }
  z = shallowconcat(matunit, z);
  u = lllintern(z, 100, 1, prec);
  if (u)
  {
    ru = lg(u);
    for (j=1; j < ru; j++)
      if (gcmp1(gcoeff(u,ru-1,j))) break;
    if (j < ru)
    {
      u = gel(u,j); /* coords on (fu^ell, be) of a small generator */
      ru--; setlg(u, ru);
      be = fix_be(bnfz, be, gmul(ell,u));
    }
  }
  if (DEBUGLEVEL>1) fprintferr("beta LLL-reduced mod U^l = %Z\n",be);
  return reducebetanaive(bnfz, be, NULL, ell);
}

static GEN
tauofalg(GEN x, GEN U) {
  return gsubst(lift(x), varn(U[1]), U);
}

static tau_s *
get_tau(tau_s *tau, GEN nf, GEN U)
{
  GEN bas = gel(nf,7), Uzk;
  long i, l = lg(bas);
  Uzk = cgetg(l, t_MAT);
  for (i=1; i<l; i++)
    gel(Uzk,i) = algtobasis(nf, tauofalg(gel(bas,i), U));
  tau->zk = Uzk;
  tau->x  = U; return tau;
}

static GEN tauoffamat(GEN x, tau_s *tau);

static GEN
tauofelt(GEN x, tau_s *tau)
{
  switch(typ(x))
  {
    case t_COL: return gmul(tau->zk, x);
    case t_MAT: return tauoffamat(x, tau);
    default: return tauofalg(x, tau->x);
  }
}
static GEN
tauofvec(GEN x, tau_s *tau)
{
  long i, l = lg(x);
  GEN y = cgetg(l, typ(x));

  for (i=1; i<l; i++) gel(y,i) = tauofelt(gel(x,i), tau);
  return y;
}
/* [x, tau(x), ..., tau^m(x)] */
static GEN
powtau(GEN x, long m, tau_s *tau)
{
  GEN y = cgetg(m+1, t_VEC);
  long i;
  gel(y,1) = x;
  for (i=2; i<=m; i++) gel(y,i) = tauofelt(gel(y,i-1), tau);
  return y;
}

static GEN
tauoffamat(GEN x, tau_s *tau)
{
  return mkmat2(tauofvec(gel(x,1), tau), gel(x,2));
}

static GEN
tauofideal(GEN nfz, GEN id, tau_s *tau)
{
  return hnfmodid(gmul(tau->zk, id), gcoeff(id, 1,1));
}

static int
isprimeidealconj(GEN nfz, GEN pr1, GEN pr2, tau_s *tau)
{
  GEN p = gel(pr1,1);
  GEN x = gel(pr1,2);
  GEN b1= gel(pr1,5);
  GEN b2= gel(pr2,5);
  if (!equalii(p, gel(pr2,1))
   || !equalii(gel(pr1,3), gel(pr2,3))
   || !equalii(gel(pr1,4), gel(pr2,4))) return 0;
  if (gequal(x,gel(pr2,2))) return 1;
  for(;;)
  {
    if (int_elt_val(nfz,x,p,b2,NULL)) return 1;
    x = FpC_red(tauofelt(x, tau), p);
    if (int_elt_val(nfz,x,p,b1,NULL)) return 0;
  }
}

static int
isconjinprimelist(GEN nfz, GEN S, GEN pr, tau_s *tau)
{
  long i, l;
  
  if (!tau) return 0;
  l = lg(S);
  for (i=1; i<l; i++)
    if (isprimeidealconj(nfz, gel(S,i),pr,tau)) return 1;
  return 0;
}

/* assume x in basistoalg form */
static GEN
downtoK(toK_s *T, GEN x)
{
  long degKz = lg(T->invexpoteta1) - 1;
  GEN y = gmul(T->invexpoteta1, RgX_to_RgV(lift_intern(x), degKz));
  return gmodulo(gtopolyrev(y,varn(T->polnf)), T->polnf);
}

static GEN
no_sol(long all, long i)
{
  if (!all) pari_err(talker,"bug%d in kummer",i);
  return cgetg(1,t_VEC);
}

static GEN
get_gell(GEN bnr, GEN subgp, long all)
{
  GEN gell;
  if (all)        gell = stoi(all);
  else if (subgp) gell = det(subgp);
  else            gell = det(diagonal_i(gmael(bnr,5,2)));
  if (typ(gell) != t_INT) pari_err(arither1);
  return gell;
}

typedef struct {
  GEN Sm, Sml1, Sml2, Sl, ESml2;
} primlist;

static int
build_list_Hecke(primlist *L, GEN nfz, GEN fa, GEN gothf, GEN gell, tau_s *tau)
{
  GEN listpr, listex, pr, p, factell;
  long vd, vp, e, i, l, ell = itos(gell), degKz = degpol(nfz[1]);

  if (!fa) fa = idealfactor(nfz, gothf);
  listpr = gel(fa,1);
  listex = gel(fa,2); l = lg(listpr);
  L->Sm  = cget1(l,t_VEC);
  L->Sml1= cget1(l,t_VEC);
  L->Sml2= cget1(l,t_VEC);
  L->Sl  = cget1(l+degKz,t_VEC);
  L->ESml2=cget1(l,t_VECSMALL);
  for (i=1; i<l; i++)
  {
    pr = gel(listpr,i); p = gel(pr,1); e = itos(gel(pr,3));
    vp = itos(gel(listex,i));
    if (!equalii(p,gell))
    {
      if (vp != 1) return 1;
      if (!isconjinprimelist(nfz, L->Sm,pr,tau)) appendL(L->Sm,pr);
    }
    else
    {
      vd = (vp-1)*(ell-1)-ell*e;
      if (vd > 0) return 4;
      if (vd==0)
      {
	if (!isconjinprimelist(nfz, L->Sml1,pr,tau)) appendL(L->Sml1, pr);
      }
      else
      {
	if (vp==1) return 2;
        if (!isconjinprimelist(nfz, L->Sml2,pr,tau))
        {
          appendL(L->Sml2, pr);
          appendL(L->ESml2,(GEN)vp);
        }
      }
    }
  }
  factell = primedec(nfz,gell); l = lg(factell);
  for (i=1; i<l; i++)
  {
    pr = gel(factell,i);
    if (!idealval(nfz,gothf,pr))
      if (!isconjinprimelist(nfz, L->Sl,pr,tau)) appendL(L->Sl, pr);
  }
  return 0; /* OK */
}

static GEN
logall(GEN nf, GEN vec, long lW, long mginv, long ell, GEN pr, long ex)
{
  GEN m, M, bid = zidealstarinitgen(nf, idealpows(nf, pr, ex));
  long ellrank, i, l = lg(vec);

  ellrank = prank(gmael(bid,2,2), ell);
  M = cgetg(l,t_MAT);
  for (i=1; i<l; i++)
  {
    m = zideallog(nf, gel(vec,i), bid);
    setlg(m, ellrank+1);
    if (i < lW) m = gmulsg(mginv, m);
    gel(M,i) = m;
  }
  return M;
}

/* compute the u_j (see remark 5.2.15.) */
static GEN
get_u(GEN cyc, long rc, GEN gell)
{
  long i, l = lg(cyc);
  GEN u = cgetg(l,t_VEC);
  for (i=1; i<=rc; i++) gel(u,i) = gen_0;
  for (   ; i<  l; i++) gel(u,i) = Fp_inv(gel(cyc,i), gell);
  return u;
}

/* alg. 5.2.15. with remark */
static GEN
isprincipalell(GEN bnfz, GEN id, GEN cycgen, GEN u, GEN gell, long rc)
{
  long i, l = lg(cycgen);
  GEN logdisc, b, y = quick_isprincipalgen(bnfz, id);

  logdisc = FpC_red(gel(y,1), gell);
  b = gel(y,2);
  for (i=rc+1; i<l; i++)
  {
    GEN e = modii(mulii(gel(logdisc,i),gel(u,i)), gell);
    if (signe(e)) b = famat_mul(b, famat_pow(gel(cycgen,i), e));
  }
  setlg(logdisc,rc+1); return mkvec2(logdisc, b);
}

static GEN
famat_factorback(GEN v, GEN e)
{
  long i, l = lg(e);
  GEN V = cgetg(1, t_MAT);
  for (i=1; i<l; i++) 
    if (signe(e[i])) V = famat_mul(V, famat_pow(gel(v,i), gel(e,i)));
  return V;
}

static GEN
compute_beta(GEN X, GEN vecWB, GEN ell, GEN bnfz)
{
  GEN BE, be;
  BE = famat_reduce(famat_factorback(vecWB, X));
  gel(BE,2) = centermod(gel(BE,2), ell);
  be = factorbackelt(BE, bnfz, NULL);
  be = reducebeta(bnfz, be, ell);
  if (DEBUGLEVEL>1) fprintferr("beta reduced = %Z\n",be);
  return be;
}

static GEN
get_Selmer(GEN bnf, GEN cycgen, long rc)
{
  GEN fu = check_units(bnf,"rnfkummer");
  GEN tu = gmael3(bnf,8,4,2);
  return shallowconcat(algtobasis(bnf,shallowconcat(fu,tu)), vecslice(cycgen,1,rc));
}


GEN
lift_if_rational(GEN x)
{
  long lx, i;
  GEN y;

  switch(typ(x))
  {
    default: break;

    case t_POLMOD:
      y = gel(x,2);
      if (typ(y) == t_POL)
      {
        long d = degpol(y);
        if (d > 0) return x;
        return (d < 0)? gen_0: gel(y,2);
      }
      return y;
  
    case t_POL: lx = lg(x);
      for (i=2; i<lx; i++) gel(x,i) = lift_if_rational(gel(x,i));
      break;
    case t_VEC: case t_COL: case t_MAT: lx = lg(x);
      for (i=1; i<lx; i++) gel(x,i) = lift_if_rational(gel(x,i));
  }
  return x;
}


/* if all!=0, give all equations of degree 'all'. Assume bnr modulus is the
 * conductor */
static GEN
rnfkummersimple(GEN bnr, GEN subgroup, GEN gell, long all)
{
  long ell, i, j, degK, dK;
  long lSml2, lSl2, lSp, rc, lW;
  long prec;

  GEN bnf,nf,bid,ideal,arch,cycgen;
  GEN clgp,cyc;
  GEN Sp,listprSp,matP;
  GEN res,u,M,K,y,vecMsup,vecW,vecWB,vecBp,msign;
  primlist L;

  bnf = gel(bnr,1);
  nf  = gel(bnf,7);
  degK = degpol(nf[1]);
  
  bid = gel(bnr,2);
  ideal= gmael(bid,1,1);
  arch = gmael(bid,1,2); /* this is the conductor */
  ell = itos(gell);
  i = build_list_Hecke(&L, nf, gel(bid,3), ideal, gell, NULL);
  if (i) return no_sol(all,i);

  lSml2 = lg(L.Sml2)-1;
  Sp = shallowconcat(L.Sm, L.Sml1); lSp = lg(Sp)-1;
  listprSp = shallowconcat(L.Sml2, L.Sl); lSl2 = lg(listprSp)-1;

  cycgen = check_and_build_cycgen(bnf);
  clgp = gmael(bnf,8,1);
  cyc = gel(clgp,2); rc = prank(cyc, ell);

  vecW = get_Selmer(bnf, cycgen, rc);
  u = get_u(cyc, rc, gell);

  vecBp = cgetg(lSp+1, t_VEC);
  matP  = cgetg(lSp+1, t_MAT);
  for (j=1; j<=lSp; j++)
  {
    GEN e, a, L;
    L = isprincipalell(bnf,gel(Sp,j), cycgen,u,gell,rc);
    e = gel(L,1); gel(matP,j) = e;
    a = gel(L,2); gel(vecBp,j) = a;
  }
  vecWB = shallowconcat(vecW, vecBp);

  prec = DEFAULTPREC +
    nbits2nlong(((degK-1) * (gexpo(vecWB) + gexpo(gmael(nf,5,1)))));
  if (nfgetprec(nf) < prec) nf = nfnewprec(nf, prec);
  msign = zsigns(nf, vecWB);

  vecMsup = cgetg(lSml2+1,t_VEC);
  M = NULL;
  for (i=1; i<=lSl2; i++)
  {
    GEN pr = gel(listprSp,i);
    long e = itos(gel(pr,3)), z = ell * (e / (ell-1));

    if (i <= lSml2)
    {
      z += 1 - L.ESml2[i];
      gel(vecMsup,i) = logall(nf, vecWB, 0,0, ell, pr,z+1);
    }
    M = vconcat(M, logall(nf, vecWB, 0,0, ell, pr,z));
  }
  lW = lg(vecW);
  M = vconcat(M, shallowconcat(zeromat(rc,lW-1), matP));

  K = FpM_ker(M, gell);
  dK = lg(K)-1;
  y = cgetg(dK+1,t_VECSMALL);
  res = cgetg(1,t_VEC); /* in case all = 1 */
  while (dK)
  {
    for (i=1; i<dK; i++) y[i] = 0;
    y[i] = 1; /* y = [0,...,0,1,0,...,0], 1 at dK'th position */
    do
    {
      pari_sp av = avma;
      GEN be, P, X = FpC_red(ZM_zc_mul(K, y), gell);
      if (ok_congruence(X, gell, lW, vecMsup) && ok_sign(X, msign, arch)) 
      {/* be satisfies all congruences, x^ell - be is irreducible, signature
        * and relative discriminant are correct */
        be = compute_beta(X, vecWB, gell, bnf);
        be = lift_if_rational(coltoalg(nf, be));
        P = gsub(monomial(gen_1, ell, 0), be);
        if (all) res = shallowconcat(res, gerepileupto(av, P));
        else
        {
          if (gequal(rnfnormgroup(bnr,P),subgroup)) return P; /*DONE*/
          avma = av;
        }
      }
      else avma = av;
    } while (increment(y, dK, ell));
    y[dK--] = 0;
  }
  if (all) return res;
  return gen_0;
}

/* alg. 5.3.11 (return only discrete log mod ell) */
static GEN
isvirtualunit(GEN bnf, GEN v, GEN cycgen, GEN cyc, GEN gell, long rc)
{
  GEN L, b, eps, y, q, nf = checknf(bnf);
  long i, l = lg(cycgen);

  L = quick_isprincipalgen(bnf, idealsqrtn(nf, v, gell, 1));
  q = gel(L,1);
  if (gcmp0(q)) { eps = v; y = q; }
  else
  {
    b = gel(L,2);
    y = cgetg(l,t_COL);
    for (i=1; i<l; i++)
      gel(y,i) = diviiexact(mulii(gell,gel(q,i)), gel(cyc,i));
    eps = famat_mul(famat_factorback(cycgen, y), famat_pow(b, gell));
    eps = famat_mul(famat_inv(eps), v);
  }
  setlg(y, rc+1);
  b = isunit(bnf,eps);
  if (lg(b) == 1) pari_err(bugparier,"isvirtualunit");
  return shallowconcat(lift_intern(b), y);
}

/* J a vector of elements in nfz = relative extension of nf by polrel,
 * return the Steinitz element associated to the module generated by J */
static GEN
Stelt(GEN nf, GEN J, GEN polrel)
{
  long i, l = lg(J);
  GEN x, A, I, id = matid(degpol(nf[1]));

  A = cgetg(l, t_VEC);
  I = cgetg(l, t_VEC);
  for (i = 1; i < l; i++)
  {
    GEN v = gel(J,i);
    gel(A,i) = (typ(v) != t_POL)? v: RgX_rem(v, polrel);
    gel(I,i) = id;
  }
  x = cgetg(3,t_VEC);
  gel(x,1) = RgXV_to_RgM(A, degpol(polrel));
  gel(x,2) = I;
  return prodid(nf, (GEN)nfhermite(nf,x)[2]);
}

static GEN
polrelKzK(toK_s *T, GEN x)
{
  GEN P = roots_to_pol(powtau(x, T->m, T->tau), 0);
  long i, l = lg(P);
  for (i=2; i<l; i++) gel(P,i) = downtoK(T, gel(P,i));
  return P;
}

/* N: Cl_m(Kz) --> Cl_m(K), lift subgroup from bnr to bnrz using Algo 4.1.11 */
static GEN
invimsubgroup(GEN bnrz, GEN bnr, GEN subgroup, toK_s *T)
{
  long l, j;
  GEN P,raycycz,rayclgpz,raygenz,U,polrel,StZk;
  GEN nf = checknf(bnr), nfz = checknf(bnrz), polz = gel(nfz,1);

  polrel = polrelKzK(T, pol_x[varn(polz)]);
  StZk = Stelt(nf, gel(nfz,7), polrel); 
  rayclgpz = gel(bnrz,5);
  raycycz = gel(rayclgpz,2); l = lg(raycycz);
  raygenz = gel(rayclgpz,3);
  P = cgetg(l,t_MAT);
  for (j=1; j<l; j++)
  {
    GEN g, id = idealhermite(nfz, gel(raygenz,j));
    g = Stelt(nf, gmul(gel(nfz,7), id), polrel);
    g = idealdiv(nf, g, StZk); /* N_{Kz/K}(gen[j]) */
    gel(P,j) = isprincipalray(bnr, g);
  }
  (void)hnfall_i(shallowconcat(P, subgroup), &U, 1);
  setlg(U, l); for (j=1; j<l; j++) setlg(U[j], l);
  return hnfmodid(shallowconcat(U, diagonal_i(raycycz)), gel(raycycz,1));
}

static GEN
pol_from_Newton(GEN S)
{
  long i, k, l = lg(S);
  GEN C = cgetg(l+1, t_VEC), c = C + 1;
  gel(c,0) = gen_1;
  c[1] = S[1];
  for (k = 2; k < l; k++)
  {
    GEN s = gel(S,k);
    for (i = 1; i < k; i++) s = gadd(s, gmul(gel(S,i), gel(c,k-i)));
    gel(c,k) = gdivgs(s, -k);
  }
  return gtopoly(C, 0);
}

/* - mu_b = sum_{0 <= i < m} floor(r_b r_{d-1-i} / ell) tau^i */
static GEN
get_mmu(long b, GEN r, long ell)
{
  long i, m = lg(r)-1;
  GEN M = cgetg(m+1, t_VEC);
  for (i = 0; i < m; i++) gel(M,i+1) = stoi((r[b + 1] * r[m - i]) / ell);
  return M;
}
/* theta^ell = be ^ ( sum tau^a r_{d-1-a} ) */
static GEN
get_reverse(GEN r)
{
  long i, m = lg(r)-1;
  GEN M = cgetg(m+1, t_VEC);
  for (i = 0; i < m; i++) gel(M,i+1) = stoi(r[m - i]);
  return M;
}

/* coeffs(x, a..b) in variable v >= varn(x) */
static GEN
split_pol(GEN x, long v, long a, long b)
{
  long i, l = degpol(x);
  GEN y = x + a, z;

  if (l < b) b = l;
  if (a > b || varn(x) != v) return zeropol(v);
  l = b-a + 3;
  z = cgetg(l, t_POL); z[1] = x[1];
  for (i = 2; i < l; i++) z[i] = y[i];
  return normalizepol_i(z, l);
}

/* return (den_a * z) mod (v^ell - num_a/den_a), assuming deg(z) < 2*ell
 * allow either num/den to be NULL (= 1) */
static GEN
mod_Xell_a(GEN z, long v, long ell, GEN num_a, GEN den_a)
{
  GEN z1 = split_pol(z, v, ell, degpol(z));
  GEN z0 = split_pol(z, v, 0,   ell-1); /* z = v^ell z1 + z0*/
  if (den_a) z0 = gmul(den_a, z0);
  if (num_a) z1 = gmul(num_a, z1);
  return gadd(z0, z1);
}
static GEN
to_alg(GEN nfz, GEN v)
{
  GEN z;
  if (typ(v) != t_COL) return v;
  z = gmul(gel(nfz,7), v);
  if (typ(z) == t_POL) setvarn(z, MAXVARN);
  return z;
}

/* th. 5.3.5. and prop. 5.3.9. */
static GEN
compute_polrel(GEN nfz, toK_s *T, GEN be, long g, long ell)
{
  long i, k, m = T->m, vT = fetch_var();
  GEN r, powtaubet, S, p1, root, num_t, den_t, nfzpol, powtau_prim_invbe;
  GEN prim_Rk, C_Rk, prim_root, C_root, prim_invbe, C_invbe;
  pari_timer ti;
    
  r = cgetg(m+1,t_VECSMALL); /* r[i+1] = g^i mod ell */
  r[1] = 1;
  for (i=2; i<=m; i++) r[i] = (r[i-1] * g) % ell;
  powtaubet = powtau(be, m, T->tau);
  if (DEBUGLEVEL>1) { fprintferr("Computing Newton sums: "); TIMERstart(&ti); }
  prim_invbe = Q_primitive_part(element_inv(nfz, be), &C_invbe);
  powtau_prim_invbe = powtau(prim_invbe, m, T->tau);

  root = cgetg(ell + 2, t_POL);
  root[1] = evalsigne(1) | evalvarn(0);
  for (i = 0; i < ell; i++) gel(root,2+i) = gen_0;
  for (i = 0; i < m; i++)
  { /* compute (1/be) ^ (-mu) instead of be^mu [mu << 0].
     * 1/be = C_invbe * prim_invbe */
    GEN mmu = get_mmu(i, r, ell);
    /* p1 = prim_invbe ^ -mu */
    p1 = to_alg(nfz, factorbackelt(powtau_prim_invbe, mmu, nfz));
    if (C_invbe) p1 = gmul(p1, powgi(C_invbe, sum(mmu,1,m)));
    /* root += zeta_ell^{r_i} T^{r_i} be^mu_i */
    gel(root, 2 + r[i+1]) = monomial(p1, r[i+1], vT);
  }
  /* Other roots are as above with z_ell --> z_ell^j. 
   * Treat all contents (C_*) and principal parts (prim_*) separately */
  prim_Rk = prim_root = Q_primitive_part(root, &C_root);
  C_Rk = C_root;

  /* Compute modulo X^ell - 1, T^ell - t, nfzpol(MAXVARN) */
  p1 = to_alg(nfz, factorbackelt(powtaubet, get_reverse(r), nfz));
  num_t = Q_remove_denom(p1, &den_t);

  nfzpol = shallowcopy(gel(nfz,1));
  setvarn(nfzpol, MAXVARN);
  S = cgetg(ell+1, t_VEC); /* Newton sums */
  gel(S,1) = gen_0;
  for (k = 2; k <= ell; k++)
  { /* compute the k-th Newton sum */
    pari_sp av = avma;
    GEN z, D, Rk = gmul(prim_Rk, prim_root);
    C_Rk = mul_content(C_Rk, C_root);
    Rk = mod_Xell_a(Rk, 0, ell, NULL, NULL); /* mod X^ell - 1 */
    for (i = 2; i < lg(Rk); i++)
    {
      z = mod_Xell_a(gel(Rk,i), vT, ell, num_t,den_t); /* mod T^ell - t */
      gel(Rk,i) = RgXQX_red(z, nfzpol); /* mod nfz.pol */
    }
    if (den_t) C_Rk = mul_content(C_Rk, ginv(den_t));
    prim_Rk = Q_primitive_part(Rk, &D);
    C_Rk = mul_content(C_Rk, D); /* root^k = prim_Rk * C_Rk */

    /* Newton sum is ell * constant coeff (in X), which has degree 0 in T */
    z = polcoeff_i(prim_Rk, 0, 0);
    z = polcoeff_i(z      , 0,vT);
    z = downtoK(T, gmulgs(z, ell));
    if (C_Rk) z = gmul(z, C_Rk);
    gerepileall(av, C_Rk? 3: 2, &z, &prim_Rk, &C_Rk);
    if (DEBUGLEVEL>1) { fprintferr("%ld(%ld) ", k, TIMER(&ti)); flusherr(); }
    gel(S,k) = z;
  }
  if (DEBUGLEVEL>1) fprintferr("\n");
  (void)delete_var(); return pol_from_Newton(S);
}

typedef struct {
  GEN R; /* compositum(P,Q) */
  GEN p; /* Mod(p,R) root of P */
  GEN q; /* Mod(q,R) root of Q */
  GEN k; /* Q[X]/R generated by q + k p */
  GEN rev;
} compo_s;

static GEN
lifttoKz(GEN nfz, GEN nf, GEN id, compo_s *C)
{
  GEN I = ideal_two_elt(nf,id);
  GEN x = coltoliftalg(nf, gel(I,2));
  gel(I,2) = algtobasis(nfz, RgX_RgXQ_compo(x, C->p, C->R));
  return prime_to_ideal(nfz,I);
}
  
static void
compositum_red(compo_s *C, GEN P, GEN Q)
{
  GEN p, q, a, z = (GEN)compositum2(P, Q)[1];
  long v;
  a = gel(z,1); v = varn(a);
  p = lift_intern(gel(z,2));
  q = lift_intern(gel(z,3));
  C->k = gel(z,4);
  /* reduce R */
  z = polredabs0(a, nf_ORIG|nf_PARTIALFACT);
  C->R = gel(z,1);
  a    = gel(z,2);
  C->p = poleval(p, a);
  if (C->p == gen_0) C->p = mkpolmod(zeropol(v),gel(a,1));
  C->q = poleval(q, a); 
  C->rev = modreverse_i(gel(a,2), gel(a,1));
  if (DEBUGLEVEL>1) fprintferr("polred(compositum) = %Z\n",C->R);
}

static GEN
_rnfkummer(GEN bnr, GEN subgroup, long all, long prec)
{
  long ell, i, j, m, d, dK, dc, rc, ru, rv, g, mginv, degK, degKz, vnf;
  long lSp, lSml2, lSl2, lW;
  GEN polnf,bnf,nf,bnfz,nfz,bid,ideal,cycgen,gell,p1,p2,wk,U,vselmer;
  GEN clgp,cyc,gen;
  GEN Q,idealz,gothf;
  GEN res,u,M,K,y,vecMsup,vecW,vecWA,vecWB,vecB,vecC,vecAp,vecBp;
  GEN matP,Sp,listprSp,Tc,Tv,P;
  primlist L;
  toK_s T;
  tau_s _tau, *tau;
  compo_s COMPO;
  pari_timer t;

  if (DEBUGLEVEL) TIMERstart(&t);
  checkbnrgen(bnr);
  bnf = gel(bnr,1);
  nf  = gel(bnf,7);
  polnf = gel(nf,1); vnf = varn(polnf);
  if (!vnf) pari_err(talker,"main variable in kummer must not be x");
  wk = gmael3(bnf,8,4,1);
  /* step 7 */
  p1 = conductor(bnr, subgroup, 2);
  if (DEBUGLEVEL) msgTIMER(&t, "[rnfkummer] conductor");
  bnr      = gel(p1,2); 
  subgroup = gel(p1,3);
  gell = get_gell(bnr,subgroup,all);
  ell = itos(gell);
  if (ell == 1) return pol_x[vnf];
  if (!uisprime(ell)) pari_err(impl,"kummer for composite relative degree");
  if (!umodiu(wk,ell)) return rnfkummersimple(bnr, subgroup, gell, all);

  bid = gel(bnr,2);
  ideal = gmael(bid,1,1);
  /* step 1 of alg 5.3.5. */
  if (DEBUGLEVEL>2) fprintferr("Step 1\n");
  compositum_red(&COMPO, polnf, cyclo(ell,vnf));
  /* step 2 */
  if (DEBUGLEVEL>2) fprintferr("Step 2\n");
  if (DEBUGLEVEL) msgTIMER(&t, "[rnfkummer] compositum");
  degK  = degpol(polnf);
  degKz = degpol(COMPO.R);
  m = degKz / degK;
  d = (ell-1) / m;
  g = (long)Fl_pow(gener_Fl(ell), d, ell);
  if (Fl_pow((ulong)g, m, ell*ell) == 1) g += ell;
  /* ord(g) = m in all (Z/ell^k)^* */
  /* step 3 */
  if (DEBUGLEVEL>2) fprintferr("Step 3\n");
  /* could factor disc(R) using th. 2.1.6. */
  bnfz = bnfinit0(COMPO.R,1,NULL,prec);
  if (DEBUGLEVEL) msgTIMER(&t, "[rnfkummer] bnfinit(Kz)");
  cycgen = check_and_build_cycgen(bnfz);
  nfz = gel(bnfz,7);
  clgp = gmael(bnfz,8,1);
  cyc = gel(clgp,2); rc = prank(cyc,ell);
  gen = gel(clgp,3);
  u = get_u(cyc, rc, gell);

  vselmer = get_Selmer(bnfz, cycgen, rc);
  if (DEBUGLEVEL) msgTIMER(&t, "[rnfkummer] Selmer group");
  ru = (degKz>>1)-1;
  rv = rc+ru+1;

  /* compute action of tau */
  U = gadd(gpowgs(COMPO.q, g), gmul(COMPO.k, COMPO.p));
  U = poleval(COMPO.rev, U);
  tau = get_tau(&_tau, nfz, U);

  /* step 4 */
  if (DEBUGLEVEL>2) fprintferr("Step 4\n");
  vecB=cgetg(rc+1,t_VEC);
  Tc=cgetg(rc+1,t_MAT);
  for (j=1; j<=rc; j++)
  {
    p1 = tauofideal(nfz, gel(gen,j), tau);
    p1 = isprincipalell(bnfz, p1, cycgen,u,gell,rc);
    Tc[j]  = p1[1];
    vecB[j]= p1[2];
  }

  vecC = cgetg(rc+1,t_VEC);
  if (rc)
  {
    for (j=1; j<=rc; j++) gel(vecC,j) = cgetg(1, t_MAT);
    p1 = cgetg(m,t_VEC);
    gel(p1,1) = matid(rc);
    for (j=2; j<=m-1; j++) gel(p1,j) = gmul(gel(p1,j-1),Tc);
    p2 = vecB;
    for (j=1; j<=m-1; j++)
    {
      GEN z = FpM_red(gmulsg((j*d)%ell,gel(p1,m-j)), gell);
      p2 = tauofvec(p2, tau);
      for (i=1; i<=rc; i++)
        gel(vecC,i) = famat_mul(gel(vecC,i), famat_factorback(p2,gel(z,i)));
    }
    for (i=1; i<=rc; i++) gel(vecC,i) = famat_reduce(gel(vecC,i));
  }
  /* step 5 */
  if (DEBUGLEVEL>2) fprintferr("Step 5\n");
  Tv = cgetg(rv+1,t_MAT);
  for (j=1; j<=rv; j++)
  {
    p1 = tauofelt(gel(vselmer,j), tau);
    if (typ(p1) == t_MAT) /* famat */
      p1 = factorbackelt(gel(p1,1), FpC_red(gel(p1,2),gell), nfz);
    gel(Tv,j) = isvirtualunit(bnfz, p1, cycgen,cyc,gell,rc);
  }
  P = FpM_ker(gsubgs(Tv, g), gell);
  lW = lg(P); vecW = cgetg(lW,t_VEC);
  for (j=1; j<lW; j++) gel(vecW,j) = famat_factorback(vselmer, gel(P,j));
  /* step 6 */
  if (DEBUGLEVEL>2) fprintferr("Step 6\n");
  Q = FpM_ker(gsubgs(shallowtrans(Tc), g), gell);
  /* step 8 */
  if (DEBUGLEVEL>2) fprintferr("Step 8\n");
  p1 = RgX_powers(lift_intern(COMPO.p), COMPO.R, degK-1);
  p1 = RgXV_to_RgM(p1, degKz);
  T.invexpoteta1 = invmat(p1); /* left inverse */
  T.polnf = polnf;
  T.tau = tau;
  T.m = m;
 
  idealz = lifttoKz(nfz, nf, ideal, &COMPO);
  if (smodis(gcoeff(ideal,1,1), ell)) gothf = idealz;
  else
  { /* ell | N(ideal) */
    GEN bnrz = buchrayinitgen(bnfz, idealz);
    GEN subgroupz = invimsubgroup(bnrz, bnr, subgroup, &T);
    gothf = conductor(bnrz,subgroupz,0);
  }
  /* step 9, 10, 11 */
  if (DEBUGLEVEL>2) fprintferr("Step 9, 10 and 11\n");
  i = build_list_Hecke(&L, nfz, NULL, gothf, gell, tau);
  if (i) return no_sol(all,i);

  lSml2 = lg(L.Sml2)-1;
  Sp = shallowconcat(L.Sm, L.Sml1); lSp = lg(Sp)-1;
  listprSp = shallowconcat(L.Sml2, L.Sl); lSl2 = lg(listprSp)-1;

  /* step 12 */
  if (DEBUGLEVEL>2) fprintferr("Step 12\n");
  vecAp = cgetg(lSp+1, t_VEC);
  vecBp = cgetg(lSp+1, t_VEC);
  matP  = cgetg(lSp+1, t_MAT);
  for (j=1; j<=lSp; j++)
  {
    GEN e, a, ap;
    p1 = isprincipalell(bnfz, gel(Sp,j), cycgen,u,gell,rc);
    e = gel(p1,1); gel(matP,j) = e;
    a = gel(p1,2);
    p2 = famat_mul(famat_factorback(vecC, gneg(e)), a);
    gel(vecBp,j) = p2;
    ap = cgetg(1, t_MAT);
    for (i=0; i<m; i++)
    {
      ap = famat_mul(ap, famat_pow(p2, utoi(Fl_pow(g,m-1-i,ell))));
      if (i < m-1) p2 = tauofelt(p2, tau);
    }
    gel(vecAp,j) = ap;
  }
  /* step 13 */
  if (DEBUGLEVEL>2) fprintferr("Step 13\n");
  vecWA = shallowconcat(vecW, vecAp);
  vecWB = shallowconcat(vecW, vecBp);

  /* step 14, 15, and 17 */
  if (DEBUGLEVEL>2) fprintferr("Step 14, 15 and 17\n");
  mginv = (m * Fl_inv(g,ell)) % ell;
  vecMsup = cgetg(lSml2+1,t_VEC);
  M = NULL;
  for (i=1; i<=lSl2; i++)
  {
    GEN pr = gel(listprSp,i);
    long e = itos(gel(pr,3)), z = ell * (e / (ell-1));

    if (i <= lSml2)
    {
      z += 1 - L.ESml2[i];
      gel(vecMsup,i) = logall(nfz, vecWA,lW,mginv,ell, pr,z+1);
    }
    M = vconcat(M, logall(nfz, vecWA,lW,mginv,ell, pr,z));
  }
  dc = lg(Q)-1;
  if (dc)
  {
    GEN QtP = gmul(shallowtrans(Q), matP);
    M = vconcat(M, shallowconcat(zeromat(dc,lW-1), QtP));
  }
  if (!M) M = zeromat(1, lSp + lW - 1);
  /* step 16 */
  if (DEBUGLEVEL>2) fprintferr("Step 16\n");
  K = FpM_ker(M, gell);
  /* step 18 & ff */
  if (DEBUGLEVEL>2) fprintferr("Step 18\n");
  dK = lg(K)-1;
  y = cgetg(dK+1,t_VECSMALL);
  res = cgetg(1, t_VEC); /* in case all = 1 */
  if (DEBUGLEVEL) msgTIMER(&t, "[rnfkummer] candidate list");
  while (dK)
  {
    for (i=1; i<dK; i++) y[i] = 0;
    y[i] = 1; /* y = [0,...,0,1,0,...,0], 1 at dK'th position */
    do
    { /* cf. algo 5.3.18 */
      GEN H, be, P, X = FpC_red(ZM_zc_mul(K, y), gell);
      if (ok_congruence(X, gell, lW, vecMsup))
      {
        be = compute_beta(X, vecWB, gell, bnfz);
        P = compute_polrel(nfz, &T, be, g, ell);
        P = lift_if_rational(P);
        if (DEBUGLEVEL>1) fprintferr("polrel(beta) = %Z\n", P);
        if (!all) {
          H = rnfnormgroup(bnr, P);
          if (gequal(subgroup, H)) return P; /* DONE */
        } else {
          GEN P0 = lift(P);
          GEN g = nfgcd(P0, derivpol(P0), polnf, gel(nf,4));
          if (degpol(g)) continue;
          H = rnfnormgroup(bnr, P);
          if (!gequal(subgroup,H) && conductor(bnr, H, -1) == gen_0) continue;
        }
        res = shallowconcat(res, P);
      }
    } while (increment(y, dK, ell));
    y[dK--] = 0;
  }
  if (all) return res;
  return gen_0; /* FAIL */
}

GEN
rnfkummer(GEN bnr, GEN subgroup, long all, long prec)
{
  pari_sp av = avma;
  return gerepilecopy(av, _rnfkummer(bnr, subgroup, all, prec));
}
