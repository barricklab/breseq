/* $Id: buch3.c 7622 2006-01-23 18:46:57Z kb $

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
/*                       RAY CLASS FIELDS                          */
/*                                                                 */
/*******************************************************************/
#include "pari.h"
#include "paripriv.h"

/* Faster than Buchray (because it can use zsignunits: easier zarchstar) */
GEN
buchnarrow(GEN bnf)
{
  GEN nf, cyc, gen, GD, v, invpi, logs, p1, p2, R, basecl, met, u1, archp, clgp;
  long r1, i, j, ngen, t, lo, c;
  pari_sp av = avma;

  bnf = checkbnf(bnf);
  nf = checknf(bnf); r1 = nf_get_r1(nf);
  clgp = gmael(bnf,8,1);
  if (!r1) return gcopy(clgp);

  cyc = gel(clgp,2);
  gen = gel(clgp,3);
  v = FpM_image(zsignunits(bnf, NULL, 1), gen_2);
  t = lg(v)-1;
  if (t == r1) { avma = av; return gcopy(clgp); }

  ngen = lg(gen)-1;
  p1 = cgetg(ngen+r1-t + 1,t_COL);
  for (i=1; i<=ngen; i++) p1[i] = gen[i];
  gen = p1;
  v = archstar_full_rk(NULL, gmael(nf,5,1), ZM_to_Flm(v, 2), gen + (ngen - t));
  v = rowslice(v, t+1, r1);

  logs = cgetg(ngen+1,t_MAT);
  GD = gmael(bnf,9,3); invpi = ginv( mppi(DEFAULTPREC) );
  archp = perm_identity(r1);
  for (j=1; j<=ngen; j++)
  {
    GEN z = zsign_from_logarch(gel(GD,j), invpi, archp);
    gel(logs,j) = F2V_red_ip( gmul(v, z) );
  }
  /* [ cyc  0 ]
   * [ logs 2 ] = relation matrix for Cl_f */
  R = shallowconcat(
    vconcat(diagonal_i(cyc), logs),
    vconcat(zeromat(ngen, r1-t), gscalmat(gen_2,r1-t))
  );
 
  met = smithrel(R,NULL,&u1);
  lo = lg(R); c = lg(met);
  if (DEBUGLEVEL>3) msgtimer("smith/class group");

  basecl = cgetg(c,t_VEC);
  for (j=1; j<c; j++)
  {
    p1 = gcoeff(u1,1,j);
    p2 = idealpow(nf,gel(gen,1),p1);
    if (signe(p1) < 0) p2 = Q_primpart(p2);
    for (i=2; i<lo; i++)
    {
      p1 = gcoeff(u1,i,j);
      if (signe(p1))
      {
	p2 = idealmul(nf,p2, idealpow(nf,gel(gen,i),p1));
        p2 = Q_primpart(p2);
      }
    }
    gel(basecl,j) = p2;
  }
  return gerepilecopy(av, mkvec3(shifti(gel(clgp,1), r1-t), met,basecl));
}

/********************************************************************/
/**                                                                **/
/**                  REDUCTION MOD IDELE                           **/
/**                                                                **/
/********************************************************************/

static GEN
compute_fact(GEN nf, GEN u1, GEN gen)
{
  GEN G, basecl;
  long prec,i,j, l = lg(u1), h = lg(u1[1]); /* l > 1 */

  basecl = cgetg(l,t_VEC);
  prec = nfgetprec(nf);
  G = cgetg(3,t_VEC);
  gel(G,2) = cgetg(1,t_MAT);

  for (j=1; j<l; j++)
  {
    GEN g,e, z = NULL;
    for (i=1; i<h; i++)
    {
      e = gcoeff(u1,i,j); if (!signe(e)) continue;

      g = gel(gen,i);
      if (typ(g) != t_MAT)
      {
        if (z) 
          gel(z,2) = arch_mul(gel(z,2), to_famat_all(g, e));
        else
          z = mkvec2(NULL, to_famat_all(g, e));
        continue;
      }

      gel(G,1) = g;
      g = idealpowred(nf,G,e,prec);
      z = z? idealmulred(nf,z,g,prec): g;
    }
    gel(z,2) = famat_reduce(gel(z,2));
    gel(basecl,j) = z;
  }
  return basecl;
}

/* given two coprime integral ideals x and f (f HNF), compute "small"
 * non-zero a in x, such that a = 1 mod (f). GTM 193: Algo 4.3.3 */
static GEN
redideal(GEN nf,GEN x,GEN f)
{
  if (gcmp1(gcoeff(f,1,1))) return idealred_elt(nf, x); /* f = 1 */
  return idealaddtoone_i(nf,x,f); /* a = b mod (x f), != 0 since 1 mod f */
}

static int
too_big(GEN nf, GEN bet)
{
  GEN x = gnorm(coltoalg(nf,bet));
  switch (typ(x))
  {
    case t_INT: return absi_cmp(x, gen_1);
    case t_FRAC: return absi_cmp(gel(x,1), gel(x,2));
  }
  pari_err(bugparier, "wrong type in too_big");
  return 0; /* not reached */
}

/* GTM 193: Algo 4.3.4. Reduce x mod idele */
static GEN
_idealmodidele(GEN nf, GEN x, GEN idele, GEN sarch)
{
  pari_sp av = avma;
  GEN a,A,D,G, f = gel(idele,1);

  G = redideal(nf, x, f);
  D = redideal(nf, idealdiv(nf,G,x), f);
  A = element_div(nf,D,G);
  if (too_big(nf,A) > 0) { avma = av; return x; }
  a = set_sign_mod_idele(nf, NULL, A, idele, sarch);
  if (a != A && too_big(nf,A) > 0) { avma = av; return x; }
  return idealmul(nf, a, x);
}

GEN
idealmodidele(GEN bnr, GEN x)
{
  GEN bid = gel(bnr,2), fa2 = gel(bid,4);
  GEN idele = gel(bid,1);
  GEN sarch = (GEN)fa2[lg(fa2)-1];
  return _idealmodidele(checknf(bnr), x, idele, sarch);
}

/* v_pr(L0 * cx). tau = pr[5] or (more efficient) mult. table for pr[5] */
static long
fast_val(GEN nf,GEN L0,GEN cx,GEN pr,GEN tau)
{
  pari_sp av = avma;
  GEN p = gel(pr,1);
  long v = int_elt_val(nf,L0,p,tau,NULL);
  if (cx)
  {
    long w = ggval(cx, p);
    if (w) v += w * itos(gel(pr,3));
  }
  avma = av; return v;
}

/* x coprime to fZ, return y = x mod fZ, y integral */
static GEN
make_integral_Z(GEN x, GEN fZ)
{
  GEN d, y = Q_remove_denom(x, &d);
  if (d) y = FpC_Fp_mul(y, Fp_inv(d, fZ), fZ);
  return y;
}

/* p pi^(-1) mod f */
static GEN
get_pinvpi(GEN nf, GEN fZ, GEN p, GEN pi, GEN *v)
{
  if (!*v) {
    GEN invpi = element_inv(nf, pi);
    *v = make_integral_Z(gmul(p, invpi), mulii(p, fZ));
  }
  return *v; 
}
/* p pi^(-1) mod f */
static GEN
get_pi(GEN F, GEN pr, GEN *v)
{
  if (!*v) *v = unif_mod_fZ(pr, F);
  return *v; 
}

static GEN
compute_raygen(GEN nf, GEN u1, GEN gen, GEN bid)
{
  GEN f, fZ, basecl, module, fa, fa2, pr, t, EX, sarch, cyc, F;
  GEN *listpr, *vecpi, *vecpinvpi, *vectau;
  long i,j,l,lp;

  if (lg(u1) == 1) return cgetg(1, t_VEC);

  /* basecl = generators in factored form */
  basecl = compute_fact(nf,u1,gen);

  module = gel(bid,1);
  cyc = gmael(bid,2,2); EX = gel(cyc,1); /* exponent of (O/f)^* */
  f   = gel(module,1); fZ = gcoeff(f,1,1);
  fa  = gel(bid,3);
  fa2 = gel(bid,4); sarch = (GEN)fa2[lg(fa2)-1];
  listpr = (GEN*)fa[1]; F = init_unif_mod_fZ((GEN)listpr);

  lp = lg(listpr);
  vecpinvpi = (GEN*)cgetg(lp, t_VEC);
  vecpi  = (GEN*)cgetg(lp, t_VEC);
  vectau = (GEN*)cgetg(lp, t_VEC);
  for (i=1; i<lp; i++) 
  {
    pr = listpr[i];
    vecpi[i]    = NULL; /* to be computed if needed */
    vecpinvpi[i] = NULL; /* to be computed if needed */
    vectau[i] = eltmul_get_table(nf, gel(pr,5));
  }

  l = lg(basecl);
  for (i=1; i<l; i++)
  {
    GEN p, pi, pinvpi, dmulI, mulI, G, I, A, e, L, newL;
    long la, v, k;
    pari_sp av;
    /* G = [I, A=famat(L,e)] is a generator, I integral */
    G = gel(basecl,i);
    I = gel(G,1);
    A = gel(G,2);
      L = gel(A,1);
      e = gel(A,2);
    /* if no reduction took place in compute_fact, everybody is still coprime
     * to f + no denominators */
    if (!I)
    {
      gel(basecl,i) = famat_to_nf_modidele(nf, L, e, bid);
      continue;
    }
    if (lg(A) == 1)
    {
      gel(basecl,i) = I;
      continue;
    }

    /* compute mulI so that mulI * I coprime to f
     * FIXME: use idealcoprime ??? (Less efficient. Fix idealcoprime!) */
    dmulI = mulI = NULL;
    for (j=1; j<lp; j++)
    {
      pr = listpr[j]; 
      v  = idealval(nf, I, pr);
      if (!v) continue;
      p  = gel(pr,1);
      pi = get_pi(F, pr, &vecpi[j]);
      pinvpi = get_pinvpi(nf, fZ, p, pi, &vecpinvpi[j]);
      t = element_pow(nf, pinvpi, stoi(v));
      mulI = mulI? element_mul(nf, mulI, t): t;
      t = powiu(gel(pr,1), v);
      dmulI = dmulI? mulii(dmulI, t): t;
    }

    /* make all components of L coprime to f. 
     * Assuming (L^e * I, f) = 1, then newL^e * mulI = L^e */
    la = lg(e); newL = cgetg(la, t_VEC);
    for (k=1; k<la; k++)
    {
      GEN L0, cx, LL = algtobasis_i(nf, gel(L,k));
      L0 = Q_primitive_part(LL, &cx); /* LL = L0*cx (faster element_val) */
      for (j=1; j<lp; j++)
      {
        pr = listpr[j];
        v  = fast_val(nf, L0,cx, pr,vectau[j]); /* = val_pr(LL) */
        if (!v) continue;
        p  = gel(pr,1);
        pi = get_pi(F, pr, &vecpi[j]);
        if (v > 0)
        {
          pinvpi = get_pinvpi(nf, fZ, p, pi, &vecpinvpi[j]);
          t = element_pow(nf,pinvpi,stoi(v));
          LL = element_mul(nf, LL, t);
          LL = gdiv(LL, powiu(p, v));
        }
        else
        {
          t = element_pow(nf,pi,stoi(-v));
          LL = element_mul(nf, LL, t);
        }
      }
      gel(newL,k) = FpC_red(make_integral(nf,LL,f,listpr), fZ);
    }

    av = avma;
    /* G in nf, = L^e mod f */
    G = famat_to_nf_modideal_coprime(nf, newL, e, f, EX);
    if (mulI)
    {
      G = element_muli(nf, G, mulI);
      G = colreducemodHNF(G, gmul(f, dmulI), NULL);
    }
    G = set_sign_mod_idele(nf,A,G,module,sarch);
    I = idealmul(nf,I,G);
    if (dmulI) I = gdivexact(I, dmulI);
    /* more or less useless, but cheap at this point */
    I = _idealmodidele(nf,I,module,sarch);
    gel(basecl,i) = gerepilecopy(av, I);
  }
  return basecl;
}

/********************************************************************/
/**                                                                **/
/**                   INIT RAY CLASS GROUP                         **/
/**                                                                **/
/********************************************************************/
static GEN
get_dataunit(GEN bnf, GEN bid)
{
  GEN D, cyc = gmael(bid,2,2), U = init_units(bnf), nf = gel(bnf,7);
  long i, l;
  zlog_S S; init_zlog_bid(&S, bid);
  D = zsignunits(bnf, S.archp, 1); l = lg(D);
  for (i = 1; i < l; i++)
    gel(D,i) = vecmodii(gmul(S.U, zlog(nf, gel(U,i),gel(D,i), &S)), cyc);
  return shallowconcat(D, diagonal_i(cyc));
}

static GEN
Buchray(GEN bnf,GEN module,long flag)
{
  GEN nf, cyc, gen, Gen, u, clg, logs, p1, h, met, u1, u2, U, cycgen;
  GEN bigres, bid, cycbid, genbid, x, y, funits, H, El;
  long RU, Ri, j, ngen, lh;
  const long add_gen = flag & nf_GEN;
  const long do_init = flag & nf_INIT;
  pari_sp av = avma;

  bnf = checkbnf(bnf); nf = checknf(bnf);
  funits = check_units(bnf, "Buchray"); RU = lg(funits);
  El = Gen = NULL; /* gcc -Wall */
  bigres = gel(bnf,8);
  cyc = gmael(bigres,1,2);
  gen = gmael(bigres,1,3); ngen = lg(cyc)-1;

  bid = Idealstar(nf,module,1);
  cycbid = gmael(bid,2,2);
  genbid = gmael(bid,2,3);
  Ri = lg(cycbid)-1; lh = ngen+Ri;

  x = gmael(bid,1,1);
  if (Ri || add_gen || do_init)
  {
    GEN fx = gel(bid,3);
    El = cgetg(ngen+1,t_VEC);
    for (j=1; j<=ngen; j++)
    {
      p1 = idealcoprime_fact(nf, gel(gen,j), fx);
      if (RgV_isscalar(p1)) p1 = gel(p1,1);
      gel(El,j) = p1;
    }
  }
  if (add_gen)
  {
    Gen = cgetg(lh+1,t_VEC);
    for (j=1; j<=ngen; j++) gel(Gen,j) = idealmul(nf,gel(El,j),gel(gen,j));
    for (   ; j<=lh; j++)   Gen[j] = genbid[j - ngen];
  }
  if (!Ri)
  {
    clg = cgetg(add_gen? 4: 3,t_VEC);
    if (add_gen) gel(clg,3) = Gen;
    gel(clg,1) = gmael(bigres,1,1);
    gel(clg,2) = cyc;
    if (!do_init) return gerepilecopy(av,clg);
    y = cgetg(7,t_VEC);
    gel(y,1) = bnf;
    gel(y,2) = bid;
    gel(y,3) = El;
    gel(y,4) = matid(ngen);
    gel(y,5) = clg;
    gel(y,6) = mkvec2(cgetg(1,t_MAT), matid(RU));
    return gerepilecopy(av,y);
  }

  cycgen = check_and_build_cycgen(bnf);
  /* (log(Units)|D) * u = (0 | H) */
  H = hnfall_i( get_dataunit(bnf, bid), do_init? &u: NULL, 1);
  logs = cgetg(ngen+1, t_MAT);
  /* FIXME: cycgen[j] is not necessarily coprime to bid, but it is made coprime
   * in famat_zlog using canonical uniformizers [from bid data]: no need to
   * correct it here. The same ones will be used in bnrisprincipal. Hence
   * modification by El is useless. */
  for (j=1; j<=ngen; j++)
  {
    p1 = gel(cycgen,j);
    if (typ(El[j]) != t_INT) /* <==> != 1 */
    {
      GEN F = to_famat_all(gel(El,j), gel(cyc,j));
      p1 = arch_mul(F, p1);
    }
    gel(logs,j) = zideallog(nf, p1, bid); /* = log(Gen[j]) */
  }
  /* [ cyc  0 ]
   * [-logs H ] = relation matrix for Cl_f */
  h = shallowconcat(
    vconcat(diagonal_i(cyc), gneg_i(logs)),
    vconcat(zeromat(ngen, Ri), H)
  );
  met = smithrel(hnf(h), &U, add_gen? &u1: NULL);
  clg = cgetg(add_gen? 4: 3, t_VEC);
  gel(clg,1) = detcyc(met, &j);
  gel(clg,2) = met;
  if (add_gen) gel(clg,3) = compute_raygen(nf,u1,Gen,bid);
  if (!do_init) return gerepilecopy(av, clg);

  u2 = cgetg(Ri+1,t_MAT);
  u1 = cgetg(RU+1,t_MAT);
  for (j=1; j<=RU; j++) { u1[j]=u[j]; setlg(u[j],RU+1); }
  u += RU;
  for (j=1; j<=Ri; j++) { u2[j]=u[j]; setlg(u[j],RU+1); }
 
  /* log(Units) U2 = H (mod D)
   * log(Units) U1 = 0 (mod D) */
  u1 = lllint_ip(u1,100);
  u2 = gmul(reducemodinvertible(u2,u1), ginv(H));
  y = cgetg(7,t_VEC);
  gel(y,1) = bnf;
  gel(y,2) = bid;
  gel(y,3) = El;
  gel(y,4) = U;
  gel(y,5) = clg;
  gel(y,6) = mkvec2(u2,u1);
  return gerepilecopy(av,y);
}

GEN
buchrayinitgen(GEN bnf, GEN ideal)
{ return Buchray(bnf,ideal, nf_INIT | nf_GEN); }
GEN
buchrayinit(GEN bnf, GEN ideal)
{ return Buchray(bnf,ideal, nf_INIT); }
GEN
buchray(GEN bnf, GEN ideal)
{ return Buchray(bnf,ideal, nf_GEN); }

GEN
bnrclass0(GEN bnf, GEN ideal, long flag)
{
  switch(flag)
  {
    case 0: flag = nf_GEN; break;
    case 1: flag = nf_INIT; break;
    case 2: flag = nf_INIT | nf_GEN; break;
    default: pari_err(flagerr,"bnrclass");
  }
  return Buchray(bnf,ideal,flag);
}

GEN
bnrinit0(GEN bnf, GEN ideal, long flag)
{
  switch(flag)
  {
    case 0: flag = nf_INIT; break;
    case 1: flag = nf_INIT | nf_GEN; break;
    default: pari_err(flagerr,"bnrinit");
  }
  return Buchray(bnf,ideal,flag);
}

GEN
bnrclassno(GEN bnf,GEN ideal)
{
  GEN nf, h, D, bigres, bid, cycbid;
  pari_sp av = avma;

  bnf = checkbnf(bnf); nf = gel(bnf,7);
  bigres = gel(bnf,8); h = gmael(bigres,1,1); /* class number */
  bid = Idealstar(nf,ideal,0);
  cycbid = gmael(bid,2,2);
  if (lg(cycbid) == 1) { avma = av; return icopy(h); }
  D = get_dataunit(bnf, bid); /* (Z_K/f)^* / units ~ Z^n / D */
  return gerepileuptoint(av, mulii(h, dethnf_i(hnf(D))));
}

GEN
quick_isprincipalgen(GEN bnf, GEN x)
{ /* x \prod g[i]^(-ep[i]) = factorisation of principal ideal */
  GEN idep, gen = gmael3(bnf,8,1,3), ep = isprincipal(bnf,x);
  idep = isprincipalfact(bnf, gen, gneg(ep), x, nf_GENMAT|nf_FORCE);
  return mkvec2(ep, gel(idep,2));
}

GEN
bnrisprincipal(GEN bnr, GEN x, long flag)
{
  long i, j, c;
  pari_sp av = avma;
  GEN bnf, nf, bid, U, El, ep, p1, beta, idep, ex, rayclass, divray, genray;
  GEN alpha;

  checkbnr(bnr); rayclass = gel(bnr,5);
  divray = gel(rayclass,2); c = lg(divray);
  ex = cgetg(c,t_COL);
  if (c == 1 && !(flag & nf_GEN)) return ex;

  bnf = gel(bnr,1); nf = gel(bnf,7);
  bid = gel(bnr,2);
  El  = gel(bnr,3);
  U   = gel(bnr,4);

  if (typ(x) == t_VEC && lg(x) == 3)
  { idep = gel(x,2); x = gel(x,1); }  /* precomputed */
  else
    idep = quick_isprincipalgen(bnf, x);
  ep  = gel(idep,1);
  beta= gel(idep,2);
  j = lg(ep);
  for (i=1; i<j; i++) /* modify beta as if gen -> El.gen (coprime to bid) */
    if (typ(El[i]) != t_INT && signe(ep[i])) /* <==> != 1 */
      beta = arch_mul(to_famat_all(gel(El,i), negi(gel(ep,i))), beta);
  p1 = gmul(U, shallowconcat(ep, zideallog(nf,beta,bid)));
  ex = vecmodii(p1, divray);
  if (!(flag & nf_GEN)) return gerepileupto(av, ex);

  /* compute generator */
  if (lg(rayclass)<=3)
    pari_err(talker,"please apply bnrinit(,,1) and not bnrinit(,,0)");

  genray = gel(rayclass,3);
  p1 = isprincipalfact(bnf, genray, gneg(ex), x, nf_GENMAT | nf_FORCE);
  if (!gcmp0(gel(p1,1))) pari_err(bugparier,"isprincipalray");
  p1 = gel(p1,2); alpha = factorbackelt(p1, nf, NULL);
  if (lg(bid[5]) > 1 && lg(gmael(bid,5,1)) > 1)
  {
    GEN u = gel(bnr,6), y = gmul(gel(u,1), zideallog(nf, p1, bid));
    y = reducemodinvertible(y, gel(u,2));
    alpha = element_div(nf, alpha, factorbackelt(init_units(bnf), y, nf));
  }
  return gerepilecopy(av, mkvec2(ex,alpha));
}

GEN
isprincipalray(GEN bnr, GEN x)
{
  return bnrisprincipal(bnr,x,0);
}

GEN
isprincipalraygen(GEN bnr, GEN x)
{
  return bnrisprincipal(bnr,x,nf_GEN);
}

/* N! / N^N * (4/pi)^r2 * sqrt(|D|) */
GEN
minkowski_bound(GEN D, long N, long r2, long prec)
{
  pari_sp av = avma;
  GEN p1;
  p1 = gdiv(mpfactr(N,prec), powuu(N,N));
  p1 = gmul(p1, gpowgs(gdivsg(4,mppi(prec)), r2));
  p1 = gmul(p1, gsqrt(absi(D),prec));
  return gerepileupto(av, p1);
}

/* DK = |dK| */
static long
zimmertbound(long N,long R2,GEN DK)
{
  pari_sp av = avma;
  GEN w;
  long n;

  if (N < 2) return 1;
  if (N < 21)
  {
    static double c[19][11] = {
{/*2*/  0.6931,     0.45158},
{/*3*/  1.71733859, 1.37420604},
{/*4*/  2.91799837, 2.50091538, 2.11943331},
{/*5*/  4.22701425, 3.75471588, 3.31196660},
{/*6*/  5.61209925, 5.09730381, 4.60693851, 4.14303665},
{/*7*/  7.05406203, 6.50550021, 5.97735406, 5.47145968},
{/*8*/  8.54052636, 7.96438858, 7.40555445, 6.86558259, 6.34608077},
{/*9*/ 10.0630022,  9.46382812, 8.87952524, 8.31139202, 7.76081149},
{/*10*/11.6153797, 10.9966020, 10.3907654,  9.79895170, 9.22232770, 8.66213267},
{/*11*/13.1930961, 12.5573772, 11.9330458, 11.3210061, 10.7222412, 10.1378082},
{/*12*/14.7926394, 14.1420915, 13.5016616, 12.8721114, 12.2542699, 11.6490374,
       11.0573775},
{/*13*/16.4112395, 15.7475710, 15.0929680, 14.4480777, 13.8136054, 13.1903162,
       12.5790381},
{/*14*/18.0466672, 17.3712806, 16.7040780, 16.0456127, 15.3964878, 14.7573587,
       14.1289364, 13.5119848},
{/*15*/19.6970961, 19.0111606, 18.3326615, 17.6620757, 16.9999233, 16.3467686,
       15.7032228, 15.0699480},
{/*16*/21.3610081, 20.6655103, 19.9768082, 19.2953176, 18.6214885, 17.9558093,
       17.2988108, 16.6510652, 16.0131906},

{/*17*/23.0371259, 22.3329066, 21.6349299, 20.9435607, 20.2591899, 19.5822454,
       18.9131878, 18.2525157, 17.6007672},

{/*18*/24.7243611, 24.0121449, 23.3056902, 22.6053167, 21.9113705, 21.2242247,
       20.5442836, 19.8719830, 19.2077941, 18.5522234},

{/*19*/26.4217792, 25.7021950, 24.9879497, 24.2793271, 23.5766321, 22.8801952,
       22.1903709, 21.5075437, 20.8321263, 20.1645647},
{/*20*/28.1285704, 27.4021674, 26.6807314, 25.9645140, 25.2537867, 24.5488420,
       23.8499943, 23.1575823, 22.4719720, 21.7935548, 21.1227537}
    };
    w = gmul(dbltor(exp(-c[N-2][R2])), gsqrt(DK,DEFAULTPREC));
  }
  else
  {
    w = minkowski_bound(DK, N, R2, DEFAULTPREC);
  }
  n = itos_or_0( gceil(w) );
  if (!n) pari_err(talker,"Minkowski bound is too large");
  if (n > 500000)
      pari_warn(warner,"large Minkowski bound: certification will be VERY long");
  avma = av; return n;
}

/* return \gamma_n^n if known, an upper bound otherwise */
static GEN
hermiteconstant(long n)
{
  GEN h,h1;
  pari_sp av;

  switch(n)
  {
    case 1: return gen_1;
    case 2: return mkfrac(utoipos(4), utoipos(3));
    case 3: return gen_2;
    case 4: return utoipos(4);
    case 5: return utoipos(8);
    case 6: return mkfrac(utoipos(64), utoipos(3));
    case 7: return utoipos(64);
    case 8: return utoipos(256);
  }
  av = avma;
  h  = gpowgs(divsr(2,mppi(DEFAULTPREC)), n);
  h1 = gsqr(ggamma(gdivgs(utoipos(n+4),2),DEFAULTPREC));
  return gerepileupto(av, gmul(h,h1));
}

/* 1 if L (= nf != Q) primitive for sure, 0 if MAYBE imprimitive (may have a
 * subfield K) */
static long
isprimitive(GEN nf)
{
  long p, i, l, ep, N = degpol(nf[1]);
  GEN d,fa;

  fa = (GEN)factor(utoipos(N))[1]; /* primes | N */
  p = itos(gel(fa,1)); if (p == N) return 1; /* prime degree */

  /* N = [L:Q] = product of primes >= p, same is true for [L:K]
   * d_L = t d_K^[L:K] --> check that some q^p divides d_L */
  d = absi(gel(nf,3));
  fa = (GEN)auxdecomp(d,0)[2]; /* list of v_q(d_L). Don't check large primes */
  if (mod2(d)) i = 1;
  else
  { /* q = 2 */
    ep = itos(gel(fa,1));
    if ((ep>>1) >= p) return 0; /* 2 | d_K ==> 4 | d_K */
    i = 2;
  }
  l = lg(fa);
  for ( ; i < l; i++)
  {
    ep = itos(gel(fa,i));
    if (ep >= p) return 0;
  }
  return 1;
}

static GEN
dft_bound()
{
  if (DEBUGLEVEL>1) fprintferr("Default bound for regulator: 0.2\n");
  return dbltor(0.2);
}

static GEN
regulatorbound(GEN bnf)
{
  long N, R1, R2, R;
  GEN nf, dK, p1, c1;

  nf = gel(bnf,7); N = degpol(nf[1]);
  if (!isprimitive(nf)) return dft_bound();

  dK = absi(gel(nf,3));
  nf_get_sign(nf, &R1, &R2); R = R1+R2-1;
  c1 = (!R2 && N<12)? int2n(N & (~1UL)): powuu(N,N);
  if (cmpii(dK,c1) <= 0) return dft_bound();

  p1 = gsqr(glog(gdiv(dK,c1),DEFAULTPREC));
  p1 = divrs(gmul2n(gpowgs(divrs(mulrs(p1,3),N*(N*N-1)-6*R2),R),R2), N);
  p1 = sqrtr(gdiv(p1, hermiteconstant(R)));
  if (DEBUGLEVEL>1) fprintferr("Mahler bound for regulator: %Z\n",p1);
  return gmax(p1, dbltor(0.2));
}

/* x given by its embeddings */
GEN
norm_by_embed(long r1, GEN x)
{
  long i, ru = lg(x)-1;
  GEN p = gel(x,ru);
  if (r1 == ru)
  {
    for (i=ru-1; i>0; i--) p = gmul(p, gel(x,i));
    return p;
  }
  p = gnorm(p);
  for (i=ru-1; i>r1; i--) p = gmul(p, gnorm(gel(x,i)));
  for (      ; i>0 ; i--) p = gmul(p, gel(x,i));
  return p;
}

static int
is_unit(GEN M, long r1, GEN x)
{
  pari_sp av = avma;
  GEN Nx = ground( norm_by_embed(r1, RgM_zc_mul(M,x)) );
  int ok = is_pm1(Nx);
  avma = av; return ok;
}

/* FIXME: should use smallvectors */
static GEN
minimforunits(GEN nf, long BORNE, GEN w)
{
  const long prec = MEDDEFAULTPREC;
  long n, i, j, k, s, *x, r1, cnt = 0;
  pari_sp av = avma;
  GEN u,r,a,M;
  double p, norme, normin, normax;
  double **q,*v,*y,*z;
  double eps=0.000001, BOUND = BORNE * 1.00001;

  if (DEBUGLEVEL>=2)
  {
    fprintferr("Searching minimum of T2-form on units:\n");
    if (DEBUGLEVEL>2) fprintferr("   BOUND = %ld\n",BORNE);
    flusherr();
  }
  r1 = nf_get_r1(nf); n = degpol(nf[1]);
  minim_alloc(n+1, &q, &x, &y, &z, &v);
  M = gprec_w(gmael(nf,5,1), prec);
  a = gmul(gmael(nf,5,2), real_1(prec));
  r = sqred1_from_QR(a, prec);
  for (j=1; j<=n; j++)
  {
    v[j] = rtodbl(gcoeff(r,j,j));
    for (i=1; i<j; i++) q[i][j] = rtodbl(gcoeff(r,i,j));
  }
  normax = 0.; normin = (double)BOUND;
  s=0; k=n; y[n]=z[n]=0;
  x[n] = (long)(sqrt(BOUND/v[n]));

  for(;;)
  {
    do
    {
      if (k>1)
      {
        long l = k-1;
	z[l] = 0;
	for (j=k; j<=n; j++) z[l] = z[l]+q[l][j]*x[j];
	p = (double)x[k] + z[k];
	y[l] = y[k]+p*p*v[k];
	x[l] = (long)floor(sqrt((BOUND-y[l])/v[l])-z[l]);
        k = l;
      }
      for(;;)
      {
	p = (double)x[k] + z[k];
        if (y[k] + p*p*v[k] <= BOUND) break;
	k++; x[k]--;
      }
    }
    while (k>1);
    if (!x[1] && y[1]<=eps) break;

    if (DEBUGLEVEL>8){ fprintferr("."); flusherr(); }
    if (++cnt == 5000) return NULL; /* too expensive */

    p = (double)x[1] + z[1]; norme = y[1] + p*p*v[1] + eps;
    if (norme > normax) normax = norme;
    if (is_unit(M,r1, x)
    && (norme > 2*n  /* exclude roots of unity */
        || !RgV_isscalar(element_pow(nf, zc_to_ZC(x), w))))
    {
      if (norme < normin) normin = norme;
      if (DEBUGLEVEL>=2) { fprintferr("*"); flusherr(); }
    }
    x[k]--;
  }
  if (DEBUGLEVEL>=2){ fprintferr("\n"); flusherr(); }
  avma = av; u = cgetg(4,t_VEC);
  gel(u,1) = stoi(s<<1);
  gel(u,2) = dbltor(normax);
  gel(u,3) = dbltor(normin);
  return u;
}

#undef NBMAX
static int
is_zero(GEN x, long bitprec) { return (gexpo(x) < -bitprec); }

static int
is_complex(GEN x, long bitprec) { return !is_zero(imag_i(x), bitprec); }

/* assume M_star t_REAL
 * FIXME: what does this do ? To be rewritten */
static GEN
compute_M0(GEN M_star,long N)
{
  long m1,m2,n1,n2,n3,lr,lr1,lr2,i,j,l,vx,vy,vz,vM;
  GEN pol,p1,p2,p3,p4,p5,p6,p7,p8,p9,u,v,w,r,r1,r2,M0,M0_pro,S,P,M;
  GEN f1,f2,f3,g1,g2,g3,pg1,pg2,pg3,pf1,pf2,pf3,X,Y,Z;
  long bitprec = 24;

  if (N == 2) return gmul2n(gsqr(gach(gmul2n(M_star,-1),0)), -1);
  vM = fetch_var(); M = pol_x[vM];
  vz = fetch_var(); Z = pol_x[vz];
  vy = fetch_var(); Y = pol_x[vy];
  vx = fetch_var(); X = pol_x[vx];

  M0 = NULL; m1 = N/3;
  for (n1=1; n1<=m1; n1++)
  {
    m2 = (N-n1)>>1;
    for (n2=n1; n2<=m2; n2++)
    {
      pari_sp av = avma; n3=N-n1-n2;
      if (n1==n2 && n1==n3) /* n1 = n2 = n3 = m1 = N/3 */
      {
	p1 = divrs(M_star, m1);
	p4 = sqrtr_abs( mulrr(addsr(1,p1),subrs(p1,3)) );
        p5 = subrs(p1,1);
	u = gen_1;
        v = gmul2n(addrr(p5,p4),-1);
        w = gmul2n(subrr(p5,p4),-1);
	M0_pro=gmul2n(mulsr(m1,addrr(gsqr(logr_abs(v)),gsqr(logr_abs(w)))), -2);
	if (DEBUGLEVEL>2)
	{
	  fprintferr("[ %ld, %ld, %ld ]: %Z\n",n1,n2,n3,gprec_w(M0_pro,3));
	  flusherr();
	}
	if (!M0 || gcmp(M0_pro,M0) < 0) M0 = M0_pro;
      }
      else if (n1==n2 || n2==n3)
      { /* n3 > N/3 >= n1 */
	long k = N - 2*n2;
	p2 = gsub(M_star, gmulgs(X,n2));
	p3 = gmul(powuu(k,k), 
                  gpowgs(gsubgs(gmul(M_star,p2),k*k),n2));
	pol = gsub(p3, gmul(gmul(powuu(n2,n2),gpowgs(X,n2)),
                            gpowgs(p2, N-n2)));
	r = roots(pol, DEFAULTPREC); lr = lg(r);
	for (i=1; i<lr; i++)
	{
          S = real_i(gel(r,i));
	  if (is_complex(gel(r,i), bitprec) || signe(S) <= 0) continue;

          p4 = subrr(M_star, mulsr(n2,S));
          P = divrr(mulrr(mulsr(n2,S),p4), subrs(mulrr(M_star,p4),k*k));
          p5 = subrr(gsqr(S), gmul2n(P,2));
          if (gsigne(p5) < 0) continue;

          p6 = sqrtr(p5);
          v = gmul2n(subrr(S,p6),-1);
          if (gsigne(v) <= 0) continue;

          u = gmul2n(addrr(S,p6),-1);
          w = gpow(P, gdivgs(utoineg(n2),k), 0);
          p6 = mulsr(n2, addrr(gsqr(logr_abs(u)), gsqr(logr_abs(v))));
          M0_pro = gmul2n(addrr(p6, mulsr(k, gsqr(logr_abs(w)))),-2);
          if (DEBUGLEVEL>2)
          {
            fprintferr("[ %ld, %ld, %ld ]: %Z\n",n1,n2,n3,gprec_w(M0_pro,3));
            flusherr();
          }
          if (!M0 || gcmp(M0_pro,M0) < 0) M0 = M0_pro;
	}
      }
      else
      {
	f1 = gsub(gadd(gmulsg(n1,X),gadd(gmulsg(n2,Y),gmulsg(n3,Z))), M);
	f2 =         gmulsg(n1,gmul(Y,Z));
	f2 = gadd(f2,gmulsg(n2,gmul(X,Z)));
	f2 = gadd(f2,gmulsg(n3,gmul(X,Y)));
	f2 = gsub(f2,gmul(M,gmul(X,gmul(Y,Z))));
	f3 = gsub(gmul(gpowgs(X,n1),gmul(gpowgs(Y,n2),gpowgs(Z,n3))), gen_1);
        /* f1 = n1 X + n2 Y + n3 Z - M */
        /* f2 = n1 YZ + n2 XZ + n3 XY */
        /* f3 = X^n1 Y^n2 Z^n3 - 1*/
	g1=subres(f1,f2); g1=gdiv(g1,content(g1));
	g2=subres(f1,f3); g2=gdiv(g2,content(g2));
	g3=subres(g1,g2); g3=gdiv(g3,content(g3));
	pf1=gsubst(f1,vM,M_star); pg1=gsubst(g1,vM,M_star);
	pf2=gsubst(f2,vM,M_star); pg2=gsubst(g2,vM,M_star);
	pf3=gsubst(f3,vM,M_star); pg3=gsubst(g3,vM,M_star);
        /* g3 = Res_Y,Z(f1,f2,f3) */
	r = roots(pg3,DEFAULTPREC); lr = lg(r);
	for (i=1; i<lr; i++)
	{
          w = real_i(gel(r,i));
	  if (is_complex(gel(r,i), bitprec) || signe(w) <= 0) continue;
          p1=gsubst(pg1,vz,w);
          p2=gsubst(pg2,vz,w);
          p3=gsubst(pf1,vz,w);
          p4=gsubst(pf2,vz,w);
          p5=gsubst(pf3,vz,w);
          r1 = roots(p1, DEFAULTPREC); lr1 = lg(r1);
          for (j=1; j<lr1; j++)
          {
            v = real_i(gel(r1,j));
            if (is_complex(gel(r1,j), bitprec) || signe(v) <= 0
             || !is_zero(gsubst(p2,vy,v), bitprec)) continue;

            p7=gsubst(p3,vy,v);
            p8=gsubst(p4,vy,v);
            p9=gsubst(p5,vy,v);
            r2 = roots(p7, DEFAULTPREC); lr2 = lg(r2);
            for (l=1; l<lr2; l++)
            {
              u = real_i(gel(r2,l));
              if (is_complex(gel(r2,l), bitprec) || signe(u) <= 0
               || !is_zero(gsubst(p8,vx,u), bitprec)
               || !is_zero(gsubst(p9,vx,u), bitprec)) continue;

              M0_pro =              mulsr(n1, gsqr(logr_abs(u)));
              M0_pro = gadd(M0_pro, mulsr(n2, gsqr(logr_abs(v))));
              M0_pro = gadd(M0_pro, mulsr(n3, gsqr(logr_abs(w))));
              M0_pro = gmul2n(M0_pro,-2);
              if (DEBUGLEVEL>2)
              {
               fprintferr("[ %ld, %ld, %ld ]: %Z\n",n1,n2,n3,gprec_w(M0_pro,3));
               flusherr();
              }
              if (!M0 || gcmp(M0_pro,M0) < 0) M0 = M0_pro;
            }
          }
	}
      }
      if (!M0) avma = av; else M0 = gerepilecopy(av, M0);
    }
  }
  for (i=1;i<=4;i++) (void)delete_var();
  return M0? M0: gen_0;
}

static GEN
lowerboundforregulator_i(GEN bnf)
{
  long N,R1,R2,RU,i;
  GEN nf,M0,M,G,bound,minunit,newminunit;
  GEN vecminim,p1,pol,y;
  GEN units = check_units(bnf,"bnfcertify");

  nf = gel(bnf,7); N = degpol(nf[1]);
  nf_get_sign(nf, &R1, &R2); RU = R1+R2-1;
  if (!RU) return gen_1;

  G = gmael(nf,5,2);
  units = algtobasis(bnf,units);
  minunit = gnorml2(gmul(G, gel(units,1))); /* T2(units[1]) */
  for (i=2; i<=RU; i++)
  {
    newminunit = gnorml2(gmul(G, gel(units,i)));
    if (gcmp(newminunit,minunit) < 0) minunit = newminunit;
  }
  if (gexpo(minunit) > 30) return NULL;

  vecminim = minimforunits(nf, itos(gceil(minunit)), gmael3(bnf,8,4,1));
  if (!vecminim) return NULL;
  bound = gel(vecminim,3);
  if (DEBUGLEVEL>1)
  {
    fprintferr("M* = %Z\n", bound);
    if (DEBUGLEVEL>2)
    {
      pol = gaddgs(gsub(monomial(gen_1,N,0),monomial(bound,1,0)),N-1);
      p1 = roots(pol,DEFAULTPREC);
      y= real_i((GEN)p1[ N&1? 3: 2]);
      M0 = gmul2n(gmulsg(N*(N-1),gsqr(glog(y,DEFAULTPREC))),-2);
      fprintferr("pol = %Z\n",pol);
      fprintferr("old method: y = %Z, M0 = %Z\n",y,gprec_w(M0,3));
    }
  }
  M0 = compute_M0(bound, N);
  if (DEBUGLEVEL>1) { fprintferr("M0 = %Z\n",gprec_w(M0,3)); flusherr(); }
  M = gmul2n(gdivgs(gdiv(gpowgs(M0,RU),hermiteconstant(RU)),N),R2);
  if (gcmp(M, dbltor(0.04)) < 0) return NULL;
  M = gsqrt(M,DEFAULTPREC);
  if (DEBUGLEVEL>1)
    fprintferr("(lower bound for regulator) M = %Z\n",gprec_w(M,3));
  return M;
}

static GEN
lowerboundforregulator(GEN bnf)
{
  pari_sp av = avma;
  GEN x = lowerboundforregulator_i(bnf);
  if (!x) { avma = av; x = regulatorbound(bnf); }
  return x;
}

/* Compute a square matrix of rank length(beta) associated to a family
 * (P_i), 1<=i<=length(beta), of primes s.t. N(P_i) = 1 mod p, and
 * (P_i,beta[j]) = 1 for all i,j */
static void
primecertify(GEN bnf, GEN beta, ulong p, GEN bad)
{
  long i, j, nbcol, lb, nbqq, ra;
  GEN nf,mat,gq,LQ,newcol,g,ord,modpr;
  ulong q;

  ord = NULL; /* gcc -Wall */
  nbcol = 0; nf = gel(bnf,7);
  lb = lg(beta)-1; mat = cgetg(1,t_MAT); q = 1UL;
  for(;;)
  {
    q += 2*p;
    if (!umodiu(bad,q) || !uisprime(q)) continue;

    gq = utoipos(q);
    LQ = primedec(bnf,gq); nbqq = lg(LQ)-1;
    g = NULL;
    for (i=1; i<=nbqq; i++)
    {
      GEN mat1, Q = gel(LQ,i); if (!gcmp1(gel(Q,4))) break;
      /* Q has degree 1 */
      if (!g)
      {
        ord = Z_factor( utoipos(q-1) );
        g = gener_Fp_local(gq, gel(ord,1)); /* primitive root */
      }
      modpr = zkmodprinit(nf, Q);
      newcol = cgetg(lb+1,t_COL);
      for (j=1; j<=lb; j++)
      {
        GEN t = to_Fp_simple(nf, gel(beta,j), modpr);
        gel(newcol,j) = Fp_PHlog(t,g,gq,ord);
      }
      if (DEBUGLEVEL>3)
      {
        if (i==1) fprintferr("       generator of (Zk/Q)^*: %Z\n", g);
        fprintferr("       prime ideal Q: %Z\n",Q);
        fprintferr("       column #%ld of the matrix log(b_j/Q): %Z\n",
                   nbcol, newcol);
      }
      mat1 = shallowconcat(mat,newcol); ra = rank(mat1);
      if (ra==nbcol) continue;

      if (DEBUGLEVEL>2) fprintferr("       new rank: %ld\n",ra);
      if (++nbcol == lb) return;
      mat = mat1;
    }
  }
}

static void
check_prime(ulong p, GEN bnf, GEN cyc, GEN cycgen, GEN fu, GEN mu, GEN bad)
{
  pari_sp av = avma;
  long i,b, lc = lg(cyc), w = itos(gel(mu,1)), lf = lg(fu);
  GEN beta = cgetg(lf+lc, t_VEC);

  if (DEBUGLEVEL>1) fprintferr("  *** testing p = %lu\n",p);
  for (b=1; b<lc; b++)
  {
    if (umodiu(gel(cyc,b), p)) break; /* p \nmid cyc[b] */
    if (b==1 && DEBUGLEVEL>2) fprintferr("     p divides h(K)\n");
    beta[b] = cycgen[b];
  }
  if (w % p == 0)
  {
    if (DEBUGLEVEL>2) fprintferr("     p divides w(K)\n");
    beta[b++] = mu[2];
  }
  for (i=1; i<lf; i++) beta[b++] = fu[i];
  setlg(beta, b); /* beta = [cycgen[i] if p|cyc[i], tu if p|w, fu] */
  if (DEBUGLEVEL>3) {fprintferr("     Beta list = %Z\n",beta); flusherr();}
  primecertify(bnf,beta,p,bad); avma = av;
}

long
certifybuchall(GEN bnf)
{
  pari_sp av = avma;
  long nbgen, i, N, R1, R2;
  GEN bad, nf, reg, zu, funits, gen, cycgen, cyc;
  byteptr delta = diffptr;
  ulong bound, p;

  bnf = checkbnf(bnf); nf = gel(bnf,7);
  N=degpol(nf[1]); if (N==1) return 1;
  nf_get_sign(nf, &R1, &R2);
  funits = check_units(bnf,"bnfcertify");
  testprimes(bnf, zimmertbound(N,R2,absi(gel(nf,3))));
  reg = gmael(bnf,8,2);
  cyc = gmael3(bnf,8,1,2); nbgen = lg(cyc)-1;
  gen = gmael3(bnf,8,1,3); zu = gmael(bnf,8,4);
  bound = itou_or_0( ground(gdiv(reg, lowerboundforregulator(bnf))) );
  if (!bound) pari_err(talker,"sorry, too many primes to check");
  maxprime_check(bound);
  if (DEBUGLEVEL>1)
  {
    fprintferr("\nPHASE 2: are all primes good ?\n\n");
    fprintferr("  Testing primes <= B (= %lu)\n\n",bound); flusherr();
  }
  cycgen = check_and_build_cycgen(bnf);
  for (bad=gen_1,i=1; i<=nbgen; i++)
    bad = lcmii(bad, gcoeff(gen[i],1,1));
  for (i=1; i<=nbgen; i++)
  {
    GEN p1 = gel(cycgen,i);
    long j;
    if (typ(p1) == t_MAT)
    {
      GEN h, g = gel(p1,1);
      for (j = 1; j < lg(g); j++)
      {
        h = idealhermite(nf, gel(g,j));
        bad = lcmii(bad, gcoeff(h,1,1));
      }
    }
  }
  /* p | bad <--> p | some element occurring in cycgen[i]  */

  funits = algtobasis(nf, funits);
  zu = mkvec2(gel(zu,1), algtobasis(nf, gel(zu,2)));

  for (p = *delta++; p <= bound; ) {  
    check_prime(p,bnf,cyc,cycgen,funits,zu,bad);
    NEXT_PRIME_VIADIFF(p, delta);
  }

  if (nbgen)
  {
    GEN f = factor(gel(cyc,1)), f1 = gel(f,1);
    long nbf1 = lg(f1);
    if (DEBUGLEVEL>1) { fprintferr("  Testing primes | h(K)\n\n"); flusherr(); }
    for (i=1; i<nbf1; i++)
    {
      p = itou(gel(f1,i));
      if (p > bound) check_prime(p,bnf,cyc,cycgen,funits,zu,bad);
    }
  }
  avma = av; return 1;
}

/*******************************************************************/
/*                                                                 */
/*        RAY CLASS FIELDS: CONDUCTORS AND DISCRIMINANTS           */
/*                                                                 */
/*******************************************************************/
/* Let bnr1, bnr2 be such that mod(bnr2) | mod(bnr1), compute the
   matrix of the surjective map Cl(bnr1) ->> Cl(bnr2) */
GEN
bnrGetSurj(GEN bnr1, GEN bnr2)
{
  long l, i;
  GEN M, gen = gmael(bnr1, 5, 3);

  l = lg(gen); M = cgetg(l, t_MAT);
  for (i = 1; i < l; i++)
    gel(M,i) = isprincipalray(bnr2, gel(gen,i));
  return M;
}

/* s: <gen> = Cl_f --> Cl_f2 --> 0, H subgroup of Cl_f (generators given as
 * HNF on [gen]). Return subgroup s(H) in Cl_f2 */
static GEN
imageofgroup(GEN bnr, GEN bnr2, GEN H)
{
  GEN H2, Delta = diagonal_i(gmael(bnr2,5,2)); /* SNF structure of Cl_n */

  if (!H) return Delta;
  H2 = gmul(bnrGetSurj(bnr, bnr2), H);
  return hnf( shallowconcat(H2, Delta) ); /* s(H) in Cl_n */
}

static GEN
args_to_bnr(GEN arg0, GEN arg1, GEN arg2, GEN *subgroup, int gen)
{
  GEN bnr,bnf;

  if (typ(arg0)!=t_VEC)
    pari_err(talker,"neither bnf nor bnr in conductor or discray");
  if (!arg1) arg1 = gen_0;
  if (!arg2) arg2 = gen_0;

  switch(lg(arg0))
  {
    case 7:  /* bnr */
      bnr = arg0; (void)checkbnf(gel(bnr,1));
      *subgroup = arg1; break;

    case 11: /* bnf */
      bnf = checkbnf(arg0);
      bnr = Buchray(bnf,arg1, gen? nf_INIT | nf_GEN: nf_INIT);
      *subgroup = arg2; break;

    default: pari_err(talker,"neither bnf nor bnr in conductor or discray");
      return NULL; /* not reached */
  }
  if (!gcmp0(*subgroup))
  {
    long tx = typ(*subgroup);
    if (!is_matvec_t(tx))
      pari_err(talker,"bad subgroup in conductor or discray");
  }
  return bnr;
}

GEN
bnrconductor(GEN arg0,GEN arg1,GEN arg2,GEN all)
{
  long flag = all? itos(all): 0;
  GEN sub = arg1, bnr = args_to_bnr(arg0,arg1,arg2,&sub, flag > 0);
  return conductor(bnr,sub, flag);
}

long
bnrisconductor(GEN arg0,GEN arg1,GEN arg2)
{
  GEN sub = arg1, bnr = args_to_bnr(arg0,arg1,arg2,&sub, 0);
  return itos(conductor(bnr,sub,-1));
}

static GEN
check_subgroup(GEN bnr, GEN H, GEN *clhray, int triv_is_NULL, char *s)
{
  GEN h, D = NULL;
  if (H && gcmp0(H)) H = NULL;
  if (H)
  {
    D = diagonal_i(gmael(bnr,5,2));
    H = hnf(H);
    if (!hnfdivide(H, D)) pari_err(talker,"incorrect subgroup in %s", s);
    h = dethnf_i(H);
    if (equalii(h, *clhray)) H = NULL; else *clhray = h;
  }
  if (!H && !triv_is_NULL) H = D? D: diagonal_i(gmael(bnr,5,2));
  return H;
}

/* return bnrisprincipal(bnr, (x)), assuming z = ideallog(x) */
static GEN
ideallog_to_bnr(GEN bnr, GEN z)
{
  GEN rayclass = gel(bnr,5), U = gel(bnr,4), divray = gel(rayclass,2);
  long j, l, lU, lz;
  int col;

  if (lg(z) == 1) return z;
  col = (typ(z) == t_COL); /* else t_MAT */
  lz = col? lg(z): lg(z[1]);
  lU = lg(U);
  if (lz != lU)
  {
    if (lz == 1) return zerocol(lg(U[1]) - 1); /* lU != 1 */
    U = vecslice(U, lU-lz+1, lU-1); /* remove Cl(K) part */
  }
  z = gmul(U, z);
  if (col)
    z = vecmodii(z, divray);
  else
  {
    l = lg(z);
    for (j = 1; j < l; j++) gel(z,j) = vecmodii(gel(z,j), divray);
  }
  return z;
}
static GEN
bnr_log_gen_pr(GEN bnr, zlog_S *S, GEN nf, long e, long index)
{ return ideallog_to_bnr(bnr, log_gen_pr(S, index, nf, e)); }
static GEN
bnr_log_gen_arch(GEN bnr, zlog_S *S, long index)
{ return ideallog_to_bnr(bnr, log_gen_arch(S, index)); }

/* A \subset H ? Allow H = NULL = trivial subgroup */
static int
contains(GEN H, GEN A)
{ return H? (hnf_gauss(H, A) != NULL): gcmp0(A); }

/* (see also Discrayrel). Given a number field bnf=bnr[1], a ray class
 * group structure bnr (with generators if all > 0), and a subgroup H of the
 * ray class group, compute the conductor of H if all=0. If all > 0, compute
 * furthermore the corresponding H' and output
 * if all = 1: [[ideal,arch],[hm,cyc,gen],H']
 * if all = 2: [[ideal,arch],newbnr,H']
 * if all < 0, answer only 1 is module is the conductor, 0 otherwise. */
GEN
conductor(GEN bnr, GEN H0, long all)
{
  pari_sp av = avma;
  long j, k, l;
  GEN bnf, nf, bid, ideal, archp, clhray, bnr2, e2, e, mod, H;
  int iscond = 1;
  zlog_S S;

  if (all > 0) checkbnrgen(bnr); else checkbnr(bnr);
  bnf = gel(bnr,1);
  bid = gel(bnr,2); init_zlog_bid(&S, bid);
  clhray = gmael(bnr,5,1);
  nf = gel(bnf,7);
  H = check_subgroup(bnr, H0, &clhray, 1, "conductor");

  archp = S.archp;
  e     = S.e; l = lg(e);
  e2 = cgetg(l, t_COL);
  for (k = 1; k < l; k++)
  {
    for (j = itos(gel(e,k)); j > 0; j--)
    {
      if (!contains(H, bnr_log_gen_pr(bnr, &S, nf, j, k))) break;
      if (all < 0) { avma = av; return gen_0; }
      iscond = 0;
    }
    gel(e2,k) = stoi(j);
  }
  l = lg(archp);
  for (k = 1; k < l; k++)
  {
    if (!contains(H, bnr_log_gen_arch(bnr, &S, k))) continue;
    if (all < 0) { avma = av; return gen_0; }
    archp[k] = 0;
    iscond = 0;
  }
  if (all < 0) { avma = av; return gen_1; }
  for (j = k = 1; k < l; k++)
    if (archp[k]) archp[j++] = archp[k];
  setlg(archp, j);
  ideal = gequal(e2, e)? gmael(bid,1,1): factorbackprime(nf, S.P, e2);
  mod = mkvec2(ideal, perm_to_arch(nf, archp));
  if (!all) return gerepilecopy(av, mod);

  if (iscond)
  {
    bnr2 = bnr;
    if (!H) H = diagonal_i(gmael(bnr,5,2));
  }
  else
  {
    bnr2 = Buchray(bnf, mod, nf_INIT | nf_GEN);
    H = imageofgroup(bnr, bnr2, H);
  }
  return gerepilecopy(av, mkvec3(mod, (all == 1)? gel(bnr2,5): bnr2, H));
}

/* return the norm group corresponding to the relative extension given by
 * polrel over bnr.bnf, assuming it is abelian and the modulus of bnr is a
 * multiple of the conductor */
GEN
rnfnormgroup(GEN bnr, GEN polrel)
{
  long i, j, reldeg, nfac, k;
  pari_sp av = avma;
  GEN bnf,index,discnf,nf,raycl,group,detgroup,fa,greldeg;
  GEN famo, fac, col;
  byteptr d = diffptr;
  ulong p;

  checkbnr(bnr); bnf=gel(bnr,1); raycl=gel(bnr,5);
  nf=gel(bnf,7);
  polrel = fix_relative_pol(nf,polrel,1);
  if (typ(polrel)!=t_POL) pari_err(typeer,"rnfnormgroup");
  reldeg = degpol(polrel);
  /* reldeg-th powers are in norm group */
  greldeg = utoipos(reldeg);
  group = diagonal_i(FpC_red(gel(raycl,2), greldeg));
  for (i=1; i<lg(group); i++)
    if (!signe(gcoeff(group,i,i))) gcoeff(group,i,i) = greldeg;
  detgroup = dethnf_i(group);
  k = cmpiu(detgroup,reldeg);
  if (k < 0)
    pari_err(talker,"not an Abelian extension in rnfnormgroup?");
  if (!k) return gerepilecopy(av, group);

  discnf = gel(nf,3);
  index  = gel(nf,4);
  for (p=0 ;;)
  {
    long oldf = -1, lfa;
    /* If all pr are unramified and have the same residue degree, p =prod pr
     * and including last pr^f or p^f is the same, but the last isprincipal
     * is much easier! oldf is used to track this */

    NEXT_PRIME_VIADIFF_CHECK(p,d);
    if (!umodiu(index, p)) continue; /* can't be treated efficiently */

    fa = primedec(nf, utoipos(p)); lfa = lg(fa)-1;
    for (i=1; i<=lfa; i++)
    {
      GEN pr = gel(fa,i), pp, T, polr, modpr;
      long f;

      /* primes of degree 1 are enough, and simpler */
      if (itos(gel(pr,4)) > 1) break;

      modpr = nf_to_ff_init(nf, &pr, &T, &pp);
      polr = modprX(polrel, nf, modpr);
      /* if pr (probably) ramified, we have to use all (non-ram) P | pr */
      if (!FqX_is_squarefree(polr, T,pp)) { oldf = 0; continue; }

      famo = FqX_factor(polr, T, pp);
      fac = gel(famo,1); f = degpol(gel(fac,1));
      nfac = lg(fac)-1;
      /* check decomposition of pr has Galois type */
      for (j=2; j<=nfac; j++)
        if (degpol(fac[j]) != f)
          pari_err(talker,"non Galois extension in rnfnormgroup");
      if (oldf < 0) oldf = f; else if (oldf != f) oldf = 0;
      if (f == reldeg) continue; /* reldeg-th powers already included */

      if (oldf && i == lfa && !umodiu(discnf, p)) pr = utoipos(p);

      /* pr^f = N P, P | pr, hence is in norm group */
      col = gmulsg(f, bnrisprincipal(bnr,pr,0));
      group = hnf(shallowconcat(group, col));
      detgroup = dethnf_i(group);
      k = cmpiu(detgroup,reldeg);
      if (k < 0) pari_err(talker,"not an Abelian extension in rnfnormgroup");
      if (!k) { cgiv(detgroup); return gerepileupto(av,group); }
    }
  }
}

static GEN
liftpol(GEN pol, GEN q)
{
  long i, l = lg(pol);
  GEN y = cgetg(l, t_POL); y[1] = pol[1];
  for (i = 2; i < l; i++)
    gel(y,i) = lift_intern(poleval(lift_intern(gel(pol,i)), q));
  return y;
}

static int
rnf_is_abelian(GEN nf, GEN pol)
{
  GEN modpr, pr, T, pp, ro, nfL, eq, C, z, a, sig;
  long i, j, l, v = varn(nf[1]);
  ulong p, k, ka;

  eq = rnfequation2(nf,pol);
  C =   shallowcopy(gel(eq,1)); setvarn(C, v);
  a = lift_intern(gel(eq,2)); setvarn(a, v); /* root of nf[1] */
  nfL = initalg_i(C, nf_PARTIALFACT, DEFAULTPREC);
  z = nfrootsall_and_pr(nfL, liftpol(pol, a));
  if (!z) return 0;
  ro = gel(z,1); l = lg(ro)-1;
  /* small groups are abelian, as are groups of prime order */
  if (l < 6 || uisprime(l)) return 1;

  pr = gel(z,2);
  modpr = nf_to_ff_init(nfL, &pr, &T, &pp);
  p = itou(pp);
  k = umodiu(gel(eq,3), p);
  ka = (k * itou(nf_to_ff(nfL, a, modpr))) % p;
  sig= cgetg(l+1, t_VECSMALL);
  /* image of c = ro[1] + k a [distinguished root of C] by the l automorphisms
   * sig[i]: ro[1] -> ro[i] */
  ro = lift_intern(ro);
  for (i = 1; i <= l; i++)
    sig[i] = Fl_add(ka, itou(nf_to_ff(nfL, gel(ro,i), modpr)), p);
  ro = Q_primpart(ro);
  for (i=2; i<=l; i++) { /* start at 2, since sig[1] = identity */
    gel(ro,i) = ZX_to_Flx(gel(ro,i), p);
    for (j=2; j<i; j++)
      if (Flx_eval(gel(ro,j), sig[i], p)
       != Flx_eval(gel(ro,i), sig[j], p)) return 0;
  }
  return 1;
}

/* Given bnf and polrel defining an abelian relative extension, compute the
 * corresponding conductor and congruence subgroup. Return
 * [[ideal,arch],[hm,cyc,gen],group] where [ideal,arch] is the conductor, and
 * [hm,cyc,gen] is the corresponding ray class group.
 * If flag != 0, check that the extension is abelian */
GEN
rnfconductor(GEN bnf, GEN polrel, long flag)
{
  pari_sp av = avma;
  GEN nf, module, bnr, group, p1, pol2;

  bnf = checkbnf(bnf); nf = gel(bnf,7);
  if (typ(polrel) != t_POL) pari_err(typeer,"rnfconductor");
  p1 = unifpol(nf, polrel, t_COL);
  pol2 = RgX_rescale(polrel, Q_denom(p1));
  if (flag && !rnf_is_abelian(nf, pol2)) { avma = av; return gen_0; }

  pol2 = fix_relative_pol(nf, pol2, 1);
  module = mkvec2((GEN)rnfdiscf(nf,pol2)[1],
                  const_vec(nf_get_r1(nf), gen_1));
  bnr   = Buchray(bnf,module,nf_INIT | nf_GEN);
  group = rnfnormgroup(bnr,pol2);
  if (!group) { avma = av; return gen_0; }
  return gerepileupto(av, conductor(bnr,group,1));
}

/* Given a number field bnf=bnr[1], a ray class group structure bnr, and a
 * subgroup H (HNF form) of the ray class group, compute [n, r1, dk]
 * associated to H (cf. discrayall). If flcond = 1, abort (return gen_0) if
 * module is not the conductor If flrel = 0, compute only N(dk) instead of
 * the ideal dk proper */
static GEN
Discrayrel(GEN bnr, GEN H0, long flag)
{
  pari_sp av = avma;
  long j, k, l, nz, flrel = flag & nf_REL, flcond = flag & nf_COND;
  GEN bnf, nf, bid, ideal, archp, clhray, clhss, P, e, dlk, H;
  zlog_S S;

  checkbnr(bnr);
  bnf = gel(bnr,1);
  bid = gel(bnr,2); init_zlog_bid(&S, bid);
  clhray = gmael(bnr,5,1);
  nf = gel(bnf,7);
  ideal= gmael(bid,1,1);
  H0 = H = check_subgroup(bnr, H0, &clhray, 0, "bnrdiscray");
  archp = S.archp;
  P     = S.P;
  e     = S.e; l = lg(e);
  dlk = flrel? idealpow(nf,ideal,clhray)
             : powgi(dethnf_i(ideal),clhray);
  for (k = 1; k < l; k++)
  {
    GEN pr = gel(P,k), sum = gen_0;
    long ep = itos(gel(e,k));
    for (j = ep; j > 0; j--)
    {
      GEN z = bnr_log_gen_pr(bnr, &S, nf, j, k);
      H = hnf(shallowconcat(H, z));
      clhss = dethnf_i(H);
      if (flcond && j==ep && equalii(clhss,clhray)) { avma = av; return gen_0; }
      if (is_pm1(clhss)) { sum = addis(sum, j); break; }
      sum = addii(sum, clhss);
    }
    dlk = flrel? idealdivpowprime(nf, dlk, pr, sum)
               : diviiexact(dlk, powgi(pr_norm(pr),sum));
  }
  l = lg(archp); nz = nf_get_r1(nf) - (l-1);
  for (k = 1; k < l; k++)
  {
    if (!contains(H0, bnr_log_gen_arch(bnr, &S, k))) continue;
    if (flcond) { avma = av; return gen_0; }
    nz++;
  }
  return gerepilecopy(av, mkvec3(clhray, stoi(nz), dlk));
}

static GEN
Discrayabs(GEN bnr, GEN subgroup, long flag)
{
  pari_sp av = avma;
  long clhray, n, R1;
  GEN z, p1, D, dk, nf, dkabs;

  D = Discrayrel(bnr, subgroup, flag);
  if ((flag & nf_REL) || D == gen_0) return D;

  nf = checknf(bnr);
  dkabs = absi(gel(nf,3));
  clhray = itos(gel(D,1)); p1 = powiu(dkabs, clhray);
  n = clhray * degpol(nf[1]);
  R1= clhray * itos(gel(D,2));
  dk = gel(D,3);
  if (((n-R1)&3) == 2) dk = negi(dk); /* (2r2) mod 4 = 2 : r2(relext) is odd */
  z = cgetg(4,t_VEC);
  gel(z,1) = utoipos(n);
  gel(z,2) = stoi(R1);
  gel(z,3) = mulii(dk,p1); return gerepileupto(av, z);
}

GEN
bnrdisc0(GEN arg0, GEN arg1, GEN arg2, long flag)
{
  GEN H, bnr = args_to_bnr(arg0,arg1,arg2,&H, 0);
  return Discrayabs(bnr,H,flag);
}
GEN
discrayrel(GEN bnr, GEN H)
{ return Discrayrel(bnr,H,nf_REL); }
GEN
discrayrelcond(GEN bnr, GEN H)
{ return Discrayrel(bnr,H,nf_REL | nf_COND); }
GEN
discrayabs(GEN bnr, GEN H)
{ return Discrayabs(bnr,H,0); }
GEN
discrayabscond(GEN bnr, GEN H)
{ return Discrayabs(bnr,H,nf_COND); }

/* chi character of abelian G: chi[i] = chi(z_i), where G = \oplus Z/cyc[i] z_i.
 * Return Ker chi [ NULL = trivial subgroup of G ] */
static GEN
KerChar(GEN chi, GEN cyc)
{
  long i, l = lg(cyc);
  GEN m, U, d1;

  if (lg(chi) != l) pari_err(talker,"incorrect character length in KerChar");
  if (l == 1) return NULL; /* trivial subgroup */
  d1 = gel(cyc,1); m = cgetg(l+1,t_MAT);
  for (i=1; i<l; i++)
  {
    if (typ(chi[i]) != t_INT) pari_err(typeer,"conductorofchar");
    gel(m,i) = mkcol(mulii(gel(chi,i), diviiexact(d1, gel(cyc,i))));
  }
  gel(m,i) = mkcol(d1);
  (void)hnfall_i(m, &U, 1);
  for (i = 1; i < l; i++) setlg(U[i], l);
  setlg(U,l); return U;
}

/* Given a number field bnf=bnr[1], a ray class group structure bnr and a
 * vector chi representing a character on the generators bnr[2][3], compute
 * the conductor of chi. */
GEN
bnrconductorofchar(GEN bnr, GEN chi)
{
  pari_sp av = avma; checkbnr(bnr);
  return gerepileupto(av, conductor(bnr, KerChar(chi, gmael(bnr,5,2)), 0));
}

/* t = [bid,U], h = #Cl(K) */
static GEN
get_classno(GEN t, GEN h)
{
  GEN bid = gel(t,1), cyc = gmael(bid,2,2);
  GEN m = shallowconcat(gel(t,2), diagonal_i(cyc));
  return mulii(h, dethnf_i(hnf(m)));
}

static void
chk_listBU(GEN L, char *s) {
  if (typ(L) != t_VEC) pari_err(typeer,s);
  if (lg(L) > 1) {
    GEN z = gel(L,1);
    if (typ(z) != t_VEC) pari_err(typeer, s);
    if (lg(z) == 1) return;
    z = gel(z,1); /* [bid,U] */
    if (typ(z) != t_VEC || lg(z) != 3) pari_err(typeer, s);
    checkbid(gel(z,1));
  }
}

/* Given lists of [zidealstarinit, unit ideallogs], return lists of ray class
 * numbers */
GEN
bnrclassnolist(GEN bnf,GEN L)
{
  pari_sp av = avma;
  long i, j, lz, l = lg(L);
  GEN v, z, V, h;

  chk_listBU(L, "bnrclassnolist");
  if (l == 1) return cgetg(1, t_VEC);
  bnf = checkbnf(bnf); h = gmael3(bnf,8,1,1);
  V = cgetg(l,t_VEC);
  for (i = 1; i < l; i++)
  {
    z = gel(L,i); lz = lg(z);
    gel(V,i) = v = cgetg(lz,t_VEC);
    for (j=1; j<lz; j++) gel(v,j) = get_classno(gel(z,j), h);
  }
  return gerepilecopy(av, V);
}

static GEN
Lbnrclassno(GEN L, GEN fac)
{
  long i, l = lg(L);
  for (i=1; i<l; i++)
    if (gequal(gmael(L,i,1),fac)) return gmael(L,i,2);
  pari_err(bugparier,"Lbnrclassno");
  return NULL; /* not reached */
}

/* returns the first index i<=n such that x=v[i] if it exits, 0 otherwise */
long
isinvector(GEN v, GEN x)
{
  long i, l = lg(v);
  for (i = 1; i < l; i++)
    if (gequal(gel(v,i), x)) return i;
  return 0;
}

static GEN
factordivexact(GEN fa1,GEN fa2)
{
  long i, j, k, c, l;
  GEN P, E, P1, E1, P2, E2, p1;

  P1 = gel(fa1,1); E1 = gel(fa1,2); l = lg(P1);
  P2 = gel(fa2,1); E2 = gel(fa2,2);
  P = cgetg(l,t_COL);
  E = cgetg(l,t_COL);
  for (c = i = 1; i < l; i++)
  {
    j = isinvector(P2,gel(P1,i));
    if (!j) { P[c] = P1[i]; E[c] = E1[i]; c++; }
    else
    {
      p1 = subii(gel(E1,i), gel(E2,j)); k = signe(p1);
      if (k < 0) pari_err(talker,"factordivexact is not exact!");
      if (k > 0) { P[c] = P1[i]; gel(E,c) = p1; c++; }
    }
  }
  setlg(P, c);
  setlg(E, c); return mkmat2(P, E);
}
/* remove index k */
static GEN
factorsplice(GEN fa, long k)
{
  GEN p = gel(fa,1), e = gel(fa,2), P, E;
  long i, l = lg(p) - 1;
  P = cgetg(l, typ(p));
  E = cgetg(l, typ(e));
  for (i=1; i<k; i++) { P[i] = p[i]; E[i] = e[i]; }
  p++; e++;
  for (   ; i<l; i++) { P[i] = p[i]; E[i] = e[i]; }
  return mkmat2(P,E);
}
static GEN
factorpow(GEN fa, long n)
{
  if (!n) return trivfact();
  return mkmat2(gel(fa,1), gmulsg(n, gel(fa,2)));
}
static GEN
factormul(GEN fa1,GEN fa2)
{
  GEN p, pnew, e, enew, v, P, y = concat_factor(fa1,fa2);
  long i, c, lx;

  p = gel(y,1); v = sindexsort(p); lx = lg(p);
  e = gel(y,2);
  pnew = vecpermute(p, v);
  enew = vecpermute(e, v);
  P = gen_0; c = 0;
  for (i=1; i<lx; i++)
  {
    if (gequal(gel(pnew,i),P))
      gel(e,c) = addii(gel(e,c),gel(enew,i));
    else
    {
      c++; P = gel(pnew,i);
      gel(p,c) = P;
      e[c] = enew[i];
    }
  }
  setlg(p, c+1);
  setlg(e, c+1); return y;
}


static long
get_nz(GEN bnf, GEN ideal, GEN arch, long clhray)
{
  GEN arch2 = shallowcopy(arch), mod = mkvec2(ideal, arch2);
  long nz = 0, l = lg(arch), k, clhss;
  for (k = 1; k < l; k++)
  { /* FIXME: this is wasteful. Use the same algorithm as conductor */
    if (signe(arch[k]))
    {
      gel(arch2,k) = gen_0; clhss = itos(bnrclassno(bnf,mod));
      gel(arch2,k) = gen_1;
      if (clhss == clhray) return -1;
    }
    else nz++;
  }
  return nz;
}

static GEN
get_NR1D(long Nf, long clhray, long degk, long nz, GEN fadkabs, GEN idealrel)
{
  long n, R1;
  GEN dlk;
  if (nz < 0) return NULL;
  n  = clhray * degk;
  R1 = clhray * nz;
  dlk = factordivexact(factorpow(factor(utoipos(Nf)),clhray), idealrel);
  /* r2 odd, set dlk = -dlk */
  if (((n-R1)&3)==2) dlk = factormul(to_famat_all(gen_m1,gen_1), dlk);
  return mkvec3(utoipos(n),
                stoi(R1),
                factormul(dlk,factorpow(fadkabs,clhray)));
}

/* t = [bid,U], h = #Cl(K) */
static GEN
get_discdata(GEN t, GEN h)
{
  GEN bid = gel(t,1), fa = gel(bid,3);
  return mkvec3(mkmat2(gel(fa,1), vec_to_vecsmall(gel(fa,2))),
                (GEN)itou(get_classno(t, h)),
                gel(bid,1));
}
typedef struct _disc_data {
  long degk;
  GEN bnf, fadk, idealrelinit, V;
} disc_data;

static GEN
get_discray(disc_data *D, GEN V, GEN x, GEN z, long N)
{
  GEN idealrel = D->idealrelinit; 
  GEN mod = gel(z,3), Fa = gel(z,1);
  GEN P = gel(Fa,1), E = gel(Fa,2);
  long k, nz, clhray = z[2], lP = lg(P);
  for (k=1; k<lP; k++)
  {
    GEN pr = gel(P,k), p = gel(pr,1);
    long e, ep = E[k], f = itos(gel(pr,4));
    long S = 0, norm = N, Npr, clhss;
    Npr = itos(powiu(p,f));
    for (e=1; e<=ep; e++)
    {
      GEN fad;
      if (e < ep) { E[k] = ep-e; fad = Fa; }
      else fad = factorsplice(Fa, k);
      norm /= Npr;
      clhss = (long)Lbnrclassno(gel(V,norm), fad);
      if (e==1 && clhss==clhray) { E[k] = ep; return cgetg(1, t_VEC); }
      if (clhss == 1) { S += ep-e+1; break; }
      S += clhss;
    }
    E[k] = ep;
    idealrel = factormul(idealrel, to_famat_all(p, utoi(f * S)));
  }
  nz = get_nz(D->bnf, gel(mod,1), gel(mod,2), clhray);
  return get_NR1D(N, clhray, D->degk, nz, D->fadk, idealrel);
}

/* Given a list of bids and associated unit log matrices, return the
 * list of discrayabs. Only keep moduli which are conductors. */
GEN
discrayabslist(GEN bnf, GEN L)
{
  pari_sp av = avma;
  long i, j, lz, l = lg(L);
  GEN nf, v, z, V, D, d, h;
  disc_data ID;

  chk_listBU(L, "discrayabslist");
  if (l == 1) return cgetg(1, t_VEC);
  ID.bnf = bnf = checkbnf(bnf);
  nf = gel(bnf,7);
  h = gmael3(bnf,8,1,1);
  ID.degk = degpol(nf[1]);
  ID.fadk = factor(absi(gel(nf,3)));
  ID.idealrelinit = trivfact();
  V = cgetg(l, t_VEC);
  D = cgetg(l, t_VEC);
  for (i = 1; i < l; i++)
  {
    z = gel(L,i); lz = lg(z);
    gel(V,i) = v = cgetg(lz,t_VEC);
    gel(D,i) = d = cgetg(lz,t_VEC);
    for (j=1; j<lz; j++) {
      gel(d,j) = get_discdata(gel(z,j), h);
      gel(v,j) = get_discray(&ID, D, gel(z,j), gel(d,j), i);
    }
  }
  return gerepilecopy(av, V);
}

/* BIG VECTOR:
 * Interface: a container v whose length is arbitrary (< 2^30), bigel(v,i)
 * refers to the i-th component. It is an lvalue.
 *
 * Implementation: a vector v whose components have exactly 2^LGVINT entries
 * but for the last one which is allowed to be shorter. v[i][j]
 * (where j<=2^LGVINT) is understood as component number I = (i-1)*2^LGVINT+j 
 * in a unique huge vector V. */

#define SHLGVINT 15
#define LGVINT (1L << SHLGVINT)
#define vext0(i) ((((i)-1)>>SHLGVINT)+1)
#define vext1(i) ((i)&(LGVINT-1))
#define bigel(v,i) gmael((v), vext0(i), vext1(i))

/* allocate an extended vector (t_VEC of t_VEC) for N _true_ components */
static GEN
bigcgetvec(long N)
{
  long i, nv = vext0(N);
  GEN v = cgetg(nv+1,t_VEC);
  for (i=1; i<nv; i++) gel(v,i) = cgetg(LGVINT+1,t_VEC);
  gel(v,nv) = cgetg(vext1(N)+1,t_VEC); return v;
}

static GEN
zsimp(GEN bid, GEN embunit)
{
  GEN empty = cgetg(1, t_VECSMALL);
  return mkvec4(mkmat2(empty,empty), gmael(bid,2,2),
                gel(bid,5), embunit);
}

static GEN
zsimpjoin(GEN b, GEN bid, GEN embunit, long prcode, long e)
{
  long i, l1, l2, nbgen, c;
  pari_sp av = avma;
  GEN fa, U, U1, U2, cyc1, cyc2, u1u2, D;

  fa = gel(b,1);
  U1 = gel(b,3);   cyc1 = gel(b,2);      l1 = lg(cyc1);
  U2 = gel(bid,5); cyc2 = gmael(bid,2,2); l2 = lg(cyc2);
  nbgen = l1+l2-2;
  if (nbgen)
  {
    u1u2 = matsnf0(diagonal_i(shallowconcat(cyc1,cyc2)), 1 | 4); /* all && clean */
    U = gel(u1u2,1);
    D = gel(u1u2,3);
    U = shallowconcat(
      l1==1   ? zeromat(nbgen, lg(U1)-1): gmul(vecslice(U, 1,   l1-1), U1),
      l1>nbgen? zeromat(nbgen, lg(U2)-1): gmul(vecslice(U, l1, nbgen), U2)
    );
  }
  else
  {
    c = lg(U1)+lg(U2)-1; U = cgetg(c,t_MAT);
    for (i=1; i<c; i++) gel(U,i) = cgetg(1,t_COL);
    D = cgetg(1,t_MAT);
  }
  return gerepilecopy(av, mkvec4(
    mkmat2(vecsmall_append(gel(fa,1), prcode),
           vecsmall_append(gel(fa,2), e)),
    mattodiagonal_i(D),
    U,
    vconcat(gel(b,4),embunit)
  ));
}

static GEN
bnrclassnointern(GEN B, GEN h)
{
  long lx = lg(B), j;
  GEN b, m, qm, L = cgetg(lx,t_VEC);
  for (j=1; j<lx; j++)
  {
    b = gel(B,j); qm = gmul(gel(b,3),gel(b,4));
    m = shallowconcat(qm, diagonal_i(gel(b,2)));
    gel(L,j) = mkvec2(gel(b,1),
                      mkvecsmall( itou( mulii(h, dethnf_i(hnf(m))) ) ));
  }
  return L;
}

static GEN
bnrclassnointernarch(GEN B, GEN h, GEN matU)
{
  long lx, nc, k, kk, j, r1, jj, nba, nbarch;
  GEN _2, b, qm, L, cyc, m, H, mm, rowsel;

  if (!matU) return bnrclassnointern(B,h);
  lx = lg(B); if (lx == 1) return B;

  r1 = lg(matU[1])-1; _2 = const_vec(r1, gen_2);
  L = cgetg(lx,t_VEC); nbarch = 1<<r1;
  for (j=1; j<lx; j++)
  {
    b = gel(B,j); qm = gmul(gel(b,3),gel(b,4));
    cyc = gel(b,2); nc = lg(cyc)-1;
    /* [ qm   cyc 0 ]
     * [ matU  0  2 ] */
    m = shallowconcat(vconcat(qm, matU),
                 diagonal_i(shallowconcat(cyc, _2)));
    m = hnf(m); mm = shallowcopy(m);
    H = cgetg(nbarch+1,t_VECSMALL);
    rowsel = cgetg(nc+r1+1,t_VECSMALL);
    for (k = 0; k < nbarch; k++)
    {
      nba = nc+1;
      for (kk=k,jj=1; jj<=r1; jj++,kk>>=1)
	if (kk&1) rowsel[nba++] = nc + jj;
      setlg(rowsel, nba);
      rowselect_p(m, mm, rowsel, nc+1);
      H[k+1] = itou( mulii(h, dethnf_i(hnf(mm))) );
    }
    gel(L,j) = mkvec2(gel(b,1), H);
  }
  return L;
}

GEN
decodemodule(GEN nf, GEN fa)
{
  long n, nn, k;
  pari_sp av = avma;
  GEN G, E, id, pr;

  nf = checknf(nf);
  if (typ(fa)!=t_MAT || lg(fa)!=3)
    pari_err(talker,"not a factorisation in decodemodule");
  n = degpol(nf[1]); nn = n*n; id = NULL;
  G = gel(fa,1);
  E = gel(fa,2);
  for (k=1; k<lg(G); k++)
  {
    long code = itos(gel(G,k)), p = code / nn, j = (code%n)+1;
    GEN P = primedec(nf, utoipos(p)), e = gel(E,k);
    if (lg(P) <= j) pari_err(talker, "incorrect hash code in decodemodule");
    pr = gel(P,j);
    id = id? idealmulpowprime(nf,id, pr,e)
           : idealpow(nf, pr,e);
  }
  if (!id) { avma = av; return matid(n); }
  return gerepileupto(av,id);
}

/* List of ray class fields. Do all from scratch, bound < 2^30. No subgroups.
 *
 * Output: a "big vector" V (cf bigcgetvec). V[k] is a vector indexed by
 * the ideals of norm k. Given such an ideal m, the component is as follows:
 *
 * + if arch = NULL, run through all possible archimedean parts; archs are
 * ordered using inverse lexicographic order, [0,..,0], [1,0,..,0], [0,1,..,0],
 * Component is [m,V] where V is a vector with 2^r1 entries, giving for each
 * arch the triple [N,R1,D], with N, R1, D as in discrayabs; D is in factored
 * form.
 *
 * + otherwise [m,N,R1,D] */
GEN
discrayabslistarch(GEN bnf, GEN arch, long bound)
{
  byteptr dif = diffptr + 1;
  int allarch = (arch==NULL), flbou = 0;
  long degk, i, j, k, sqbou, l, nba, nbarch, ii, r1, c;
  pari_sp av0 = avma,  av,  av1,  lim;
  GEN nf, p, Z, fa, ideal, bidp, matarchunit, Disc, U, sgnU, EMPTY;
  GEN res, embunit, h, Ray, discall, idealrel, idealrelinit, fadkabs;

  if (bound <= 0) pari_err(talker,"non-positive bound in Discrayabslist");
  res = discall = NULL; /* -Wall */

  bnf = checkbnf(bnf);
  nf = gel(bnf,7); r1 = nf_get_r1(nf);
  degk = degpol(nf[1]);
  fadkabs = factor(absi(gel(nf,3)));
  h = gmael3(bnf,8,1,1);
  U = init_units(bnf);
  sgnU = zsignunits(bnf, NULL, 1);

  if (allarch) arch = const_vec(r1, gen_1);
  bidp = Idealstar(nf, mkvec2(gen_1, arch), 0);
  if (allarch) {
    matarchunit = zlog_units(nf, U, sgnU, bidp);
    bidp = Idealstar(nf,matid(degk),0);
    if (r1>15) pari_err(talker,"r1>15 in discrayabslistarch");
    nba = r1;
  } else {
    matarchunit = (GEN)NULL;
    for (nba=0,i=1; i<=r1; i++) if (signe(arch[i])) nba++;
  }

  /* what follows was rewritten from Ideallist */
  p = utoipos(2);
  av = avma; lim = stack_lim(av,1);
  sqbou = (long)sqrt((double)bound) + 1;
  Z = bigcgetvec(bound);
  for (i=2; i<=bound; i++) bigel(Z,i) = cgetg(1,t_VEC);
  embunit = zlog_units(nf, U, sgnU, bidp);
  bigel(Z,1) = mkvec(zsimp(bidp,embunit)); 
  if (DEBUGLEVEL>1) fprintferr("Starting zidealstarunits computations\n");
  maxprime_check((ulong)bound);
  /* The goal is to compute Ray (lists of bnrclassno). Z contains "zsimps",
   * simplified zidealstarinit, from which bnrclassno is easy to compute.
   * Once p > sqbou, delete Z[i] for i > sqbou and compute directly Ray */
  Ray = Z;
  while (p[2] <= bound)
  {
    if (!flbou && p[2] > sqbou)
    {
      GEN z;
      flbou = 1;
      if (DEBUGLEVEL>1) fprintferr("\nStarting bnrclassno computations\n");
      Z = gerepilecopy(av,Z); av1 = avma;
      Ray = bigcgetvec(bound);
      for (i=1; i<=bound; i++)
	bigel(Ray,i) = bnrclassnointernarch(bigel(Z,i),h,matarchunit);
      Ray = gerepilecopy(av1,Ray);
      z = bigcgetvec(sqbou);
      for (i=1; i<=sqbou; i++) bigel(z,i) = bigel(Z,i);
      Z = z;
    }
    fa = primedec(nf,p);
    for (j=1; j<lg(fa); j++)
    {
      GEN pr = gel(fa,j);
      long prcode, q, f = itos(gel(pr,4)), Q = itos_or_0(powiu(p,f));
      if (!Q || Q > bound) continue;

      /* p, f-1, j-1 as a single integer in "base degk" (f,j <= degk)*/
      prcode = (p[2]*degk + f-1)*degk + j-1;
      q = Q; ideal = pr;
      for (l=1;; l++) /* Q <= bound */
      {
        bidp = Idealstar(nf,ideal,0);
        embunit = zlog_units_noarch(nf, U, bidp);
        for (i=Q; i<=bound; i+=Q)
        {
          GEN pz, p2, p1 = bigel(Z,i/Q);
          long lz = lg(p1);
          if (lz == 1) continue;

          p2 = cgetg(lz,t_VEC); c = 0;
          for (k=1; k<lz; k++)
          {
            GEN z = gel(p1,k), v = gmael(z,1,1); /* primes in zsimp's fact. */
            long lv = lg(v);
            /* If z has a power of pr in its modulus, skip it */
            if (Q != i && lv > 1 && v[lv-1] == prcode) break;
            gel(p2,++c) = zsimpjoin(z,bidp,embunit,prcode,l);
          }

          setlg(p2, c+1);
          pz = bigel(Ray,i);
          if (flbou) p2 = bnrclassnointernarch(p2,h,matarchunit);
          if (lg(pz) > 1) p2 = shallowconcat(pz,p2);
          bigel(Ray,i) = p2;
        }
        Q = itos_or_0( mulss(Q, q) );
        if (!Q || Q > bound) break;

        ideal = idealmul(nf,ideal,pr);
      }
    }
    if (low_stack(lim, stack_lim(av,1)))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"[1]: discrayabslistarch");
      gerepileall(av, flbou? 2: 1, &Z, &Ray);
    }
    NEXT_PRIME_VIADIFF(p[2], dif);
  }
  if (!flbou) /* occurs iff bound = 1,2,4 */
  {
    if (DEBUGLEVEL>1) fprintferr("\nStarting bnrclassno computations\n");
    Ray = bigcgetvec(bound);
    for (i=1; i<=bound; i++)
      bigel(Ray,i) = bnrclassnointernarch(bigel(Z,i),h,matarchunit);
  }
  Ray = gerepilecopy(av, Ray);
  
  if (DEBUGLEVEL>1) fprintferr("Starting discrayabs computations\n");
  if (allarch) nbarch = 1<<r1;
  else
  {
    nbarch = 1;
    discall = cgetg(2,t_VEC);
  }
  EMPTY = mkvec3(gen_0,gen_0,gen_0);
  idealrelinit = trivfact();
  av1 = avma; lim = stack_lim(av1,1);
  Disc = bigcgetvec(bound);
  for (i=1; i<=bound; i++) bigel(Disc,i) = cgetg(1,t_VEC);
  for (ii=1; ii<=bound; ii++)
  {
    GEN sous, sousdisc;
    long ls;
    i = ii;
    sous = bigel(Ray,i);
    ls = lg(sous); bigel(Disc,ii) = sousdisc = cgetg(ls,t_VEC);
    for (j=1; j<ls; j++)
    {
      GEN b = gel(sous,j), clhrayall = gel(b,2), Fa = gel(b,1);
      GEN P = gel(Fa,1), E = gel(Fa,2);
      long lP = lg(P), karch;

      if (allarch) discall = cgetg(nbarch+1,t_VEC);
      for (karch=0; karch<nbarch; karch++)
      {
        long nz, clhray = clhrayall[karch+1];
        if (allarch)
        {
          long ka, k2;
          nba = 0;
          for (ka=karch,k=1; k<=r1; k++,ka>>=1)
            if (ka & 1) nba++;
          for (k2=1,k=1; k<=r1; k++,k2<<=1)
            if (karch&k2 && clhrayall[karch-k2+1] == clhray)
              { res = EMPTY; goto STORE; }
        }
        idealrel = idealrelinit;
        for (k=1; k<lP; k++) /* cf get_discray */
        {
          long e, ep = E[k], pf = P[k] / degk, f = (pf%degk) + 1;
          long S = 0, normi = i, Npr, clhss;
          p = utoipos(pf / degk);
          Npr = itos(powiu(p,f));
          for (e=1; e<=ep; e++)
          {
            GEN fad;
            if (e < ep) { E[k] = ep-e; fad = Fa; }
            else fad = factorsplice(Fa, k);
            normi /= Npr;
            clhss = Lbnrclassno(bigel(Ray,normi),fad)[karch+1];
            if (e==1 && clhss==clhray) { E[k] = ep; res = EMPTY; goto STORE; }
            if (clhss == 1) { S += ep-e+1; break; }
            S += clhss;
          }
          E[k] = ep;
          idealrel = factormul(idealrel, to_famat_all(p, utoi(f * S)));
        }
        if (!allarch && nba)
          nz = get_nz(bnf, decodemodule(nf,Fa), arch, clhray);
        else
          nz = r1 - nba;
        res = get_NR1D(i, clhray, degk, nz, fadkabs, idealrel);
STORE:  gel(discall,karch+1) = res;
      }
      res = allarch? mkvec2(Fa, discall)
                   : mkvec4(Fa, gel(res,1), gel(res,2), gel(res,3));
      gel(sousdisc,j) = res;
      if (low_stack(lim, stack_lim(av1,1)))
      {
        if(DEBUGMEM>1) pari_warn(warnmem,"[2]: discrayabslistarch");
        Disc = gerepilecopy(av1, Disc);
      }
    }
  }
  return gerepilecopy(av0, Disc);
}
GEN
discrayabslistlong(GEN bnf, long bound) {
  GEN nf = checknf(bnf);
  long r1 = nf_get_r1(nf);
  return discrayabslistarch(bnf,zerovec(r1),bound);
}

static GEN
subgroupcond(GEN bnr, GEN indexbound)
{
  pari_sp av = avma;
  long i, k, l, le, la;
  GEN e, li, p1, lidet, perm, L, nf = checknf(bnr);
  zlog_S S;

  checkbnr(bnr); init_zlog_bid(&S, gel(bnr,2));
  e = S.e; le = lg(e); la = lg(S.archp);
  L = cgetg(le + la - 1, t_VEC);
  i = 1;
  for (k = 1; k < le; k++)
    gel(L,i++) = bnr_log_gen_pr(bnr, &S, nf, itos(gel(e,k)), k);
  for (k = 1; k < la; k++)
    gel(L,i++) = bnr_log_gen_arch(bnr, &S, k);
  li = subgroupcondlist(gmael(bnr,5,2), indexbound, L);
  l = lg(li);
  /* sort by increasing index */
  lidet = cgetg(l,t_VEC);
  for (i=1; i<l; i++) gel(lidet,i) = dethnf_i(gel(li,i));
  perm = sindexsort(lidet); p1 = li; li = cgetg(l,t_VEC);
  for (i=1; i<l; i++) li[i] = p1[perm[l-i]];
  return gerepilecopy(av,li);
}

GEN
subgrouplist0(GEN bnr, GEN indexbound, long all)
{
  if (typ(bnr)!=t_VEC) pari_err(typeer,"subgrouplist");
  if (lg(bnr)!=1 && typ(bnr[1])!=t_INT)
  {
    if (!all) return subgroupcond(bnr,indexbound);
    checkbnr(bnr); bnr = gmael(bnr,5,2);
  }
  return subgrouplist(bnr,indexbound);
}

GEN
bnrdisclist0(GEN bnf, GEN L, GEN arch)
{
  if (typ(L)!=t_INT) return discrayabslist(bnf,L);
  return discrayabslistarch(bnf,arch,itos(L));
}
