/* $Id: stark.c 7857 2006-04-11 17:28:55Z kb $

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
/*        COMPUTATION OF STARK UNITS OF TOTALLY REAL FIELDS        */
/*                                                                 */
/*******************************************************************/
#include "pari.h"
#include "paripriv.h"

#define EXTRA_PREC (DEFAULTPREC-1)
#define ADD_PREC   (DEFAULTPREC-2)*3

/* ComputeCoeff */
typedef struct {
  GEN L0, L1, L11, L2; /* VECSMALL of p */
  GEN *L1ray, *L11ray; /* precomputed isprincipalray(pr), pr | p */
  GEN *rayZ; /* precomputed isprincipalray(i), i < condZ */
  long condZ; /* generates cond(bnr) \cap Z, assumed small */
} LISTray;

/* Char evaluation */
typedef struct {
  long ord;
  GEN *val, chi;
} CHI_t;

/* RecCoeff */
typedef struct {
  GEN M, beta, B, U, nB;
  long v, G, N;
} RC_data;

/********************************************************************/
/*                    Miscellaneous functions                       */
/********************************************************************/
/* exp(2iPi/den), assume den a t_INT */
static GEN
InitRU(GEN den, long prec)
{
  GEN c, s;
  if (equaliu(den, 2)) return gen_m1;
  gsincos(divri(Pi2n(1, prec), den), &s, &c, prec);
  return mkcomplex(c, s);
}
/* Compute the image of logelt by character chi, as a complex number */
static GEN
ComputeImagebyChar(GEN chi, GEN logelt)
{
  GEN gn = gmul(gel(chi,1), logelt), x = gel(chi,2);
  long d = itos(gel(chi,3)), n = smodis(gn, d);
  /* x^d = 1 and, if d even, x^(d/2) = -1 */
  if ((d & 1) == 0)
  {
    d /= 2;
    if (n >= d) return gneg(gpowgs(x, n-d));
  }
  return gpowgs(x, n);
}

/* return n such that C(elt) = z^n */
static ulong
EvalChar_n(CHI_t *C, GEN logelt)
{
  GEN n = gmul(C->chi, logelt);
  return umodiu(n, C->ord);
}
/* return C(elt) */
static GEN
EvalChar(CHI_t *C, GEN logelt)
{
  return C->val[EvalChar_n(C, logelt)];
}

static void
init_CHI(CHI_t *c, GEN CHI, GEN z)
{
  long i, d = itos(gel(CHI,3));
  GEN *v = (GEN*)new_chunk(d);
  v[0] = gen_1;
  v[1] = z;
  for (i=2; i<d; i++) v[i] = gmul(v[i-1], z);
  c->chi = gel(CHI,1);
  c->ord = d;
  c->val = v;
}
/* as t_POLMOD */
static void
init_CHI_alg(CHI_t *c, GEN CHI) {
  long d = itos(gel(CHI,3));
  GEN z;
  switch(d)
  {
    case 1: z = gen_1; break;
    case 2: z = gen_m1; break;
    default: z = mkpolmod(pol_x[0], cyclo(d,0));
  }
  init_CHI(c,CHI, z);
}
/* as t_COMPLEX */
static void
init_CHI_C(CHI_t *c, GEN CHI) {
  init_CHI(c,CHI, gel(CHI,2));
}

/* Compute the conjugate character */
static GEN
ConjChar(GEN chi, GEN cyc)
{
  long i, l = lg(chi);
  GEN z = cgetg(l, t_VEC);
  for (i = 1; i < l; i++)
    gel(z,i) = signe(chi[i])? subii(gel(cyc,i), gel(chi,i)): gen_0;
  return z;
}

typedef struct {
  long r; /* rank = lg(gen) */
  GEN j; /* current elt is gen[1]^j[1] ... gen[r]^j[r] */
  GEN cyc; /* t_VECSMALL of elementary divisors */
} GROUP_t;

static int
NextElt(GROUP_t *G)
{
  long i = 1;
  if (G->r == 0) return 0; /* no more elt */
  while (++G->j[i] == G->cyc[i]) /* from 0 to cyc[i]-1 */
  {
    G->j[i] = 0;
    if (++i > G->r) return 0; /* no more elt */
  }
  return i; /* we have multiplied by gen[i] */
}

/* Compute all the elements of a group given by its SNF */
static GEN
EltsOfGroup(long order, GEN cyc)
{
  long i;
  GEN rep;
  GROUP_t G;

  G.cyc = gtovecsmall(cyc);
  G.r = lg(cyc)-1;
  G.j = const_vecsmall(G.r, 0);

  rep = cgetg(order + 1, t_VEC);
  gel(rep,order) = vecsmall_to_col(G.j);

  for  (i = 1; i < order; i++)
  {
    (void)NextElt(&G);
    gel(rep,i) = vecsmall_to_col(G.j);
  }
  return rep;
}

/* Let dataC as given by InitQuotient, compute a system of
   representatives of the quotient */
static GEN
ComputeLift(GEN dataC)
{
  long order, i;
  pari_sp av = avma;
  GEN cyc, surj, eltq, elt;

  order = itos(gel(dataC,1));
  cyc   = gel(dataC,2);
  surj  = gel(dataC,3);

  eltq = EltsOfGroup(order, cyc);
  elt = cgetg(order + 1, t_VEC);
  for (i = 1; i <= order; i++) gel(elt,i) = inverseimage(surj, gel(eltq,i));

  return gerepileupto(av, elt);
}

/* Return c[1],  [c[1]/c[1] = 1,...,c[n]/c[1]] */
static GEN
init_get_chic(GEN c)
{
  long i, l = lg(c); /* > 1 */
  GEN C, D = cgetg(l, t_VEC);
  if (l == 1) C = gen_1;
  else
  {
    C = gel(c,1); gel(D,1) = gen_1;
    for (i = 2; i < l; i++) gel(D,i) = diviiexact(C, gel(c,i));
  }
  return mkvec2(C, D);
}

static GEN
get_chic(GEN chi, GEN D)
{
  long i, l = lg(chi);
  GEN chic = cgetg(l, t_VEC);
  gel(chic,1) = gel(chi,1);
  for (i = 2; i < l; i++) gel(chic,i) = mulii(gel(chi,i), gel(D,i));
  return chic;
}

/* A character is given by a vector [(c_i), z, d] such that
   chi(id) = z ^ sum(c_i * a_i) where
     a_i= log(id) on the generators of bnr
     z  = exp(2i * Pi / d) */
static GEN
get_Char(GEN chi, GEN initc, GEN U, long prec)
{
  GEN d, ch = cgetg(4, t_VEC), chic = get_chic(chi, gel(initc,2));
  if (U) chic = gmul(chic, U);
  chic = Q_primitive_part(chic, &d);
  if (d) {
    GEN t = gdiv(gel(initc,1), d);
    d = denom(t);
    if (!is_pm1(d)) chic = gmul(d, chic);
    d = numer(t);
  } else
    d = gel(initc,1);

  gel(ch,1) = chic;
  gel(ch,2) = InitRU(d, prec);
  gel(ch,3) = d; return ch;
}

/* prime divisors of conductor */
static GEN
divcond(GEN bnr) { GEN bid = gel(bnr,2); return gmael(bid,3,1); }

/* vector of prime ideals dividing bnr but not bnrc */
static GEN
get_prdiff(GEN bnr, GEN condc)
{
  GEN prdiff, M = gel(condc,1), D = divcond(bnr), nf = gmael(bnr, 1, 7);
  long nd, i, l  = lg(D);
  prdiff = cgetg(l, t_COL);
  for (nd=1, i=1; i < l; i++)
    if (!idealval(nf, M, gel(D,i))) prdiff[nd++] = D[i];
  setlg(prdiff, nd); return prdiff;
}

/* Let chi a character defined over bnr and primitive over bnrc, compute the
 * corresponding primitive character. Returns NULL if bnr = bnrc */
static GEN
GetPrimChar(GEN chi, GEN bnr, GEN bnrc, long prec)
{
  long l;
  pari_sp av = avma;
  GEN U, M, cond, condc, initc, Mrc;

  cond  = gmael(bnr,  2, 1);
  condc = gmael(bnrc, 2, 1); if (gequal(cond, condc)) return NULL;

  initc = init_get_chic(gmael(bnr, 5, 2));
  Mrc   = diagonal_i(gmael(bnrc, 5, 2));
  M = bnrGetSurj(bnr, bnrc);
  (void)hnfall_i(shallowconcat(M, Mrc), &U, 1);
  l = lg(M);
  U = rowslice(vecslice(U, l, lg(U)-1), 1, l-1);
  return gerepilecopy(av, get_Char(chi, initc, U, prec));
}

#define ch_chi(x)  gel(x,1)
#define ch_C(x)    gel(x,2)
#define ch_bnr(x)  gel(x,3)
#define ch_4(x)    gel(x,4)
#define ch_CHI(x)  gel(x,5)
#define ch_diff(x) gel(x,6)
#define ch_cond(x) gel(x,7)
#define ch_CHI0(x) gel(x,8)

static GEN
GetDeg(GEN dataCR)
{
  long i, l = lg(dataCR);
  GEN degs = cgetg(l, t_VECSMALL);

  for (i = 1; i < l; i++) degs[i] = itou(phi(gel(ch_CHI(gel(dataCR, i)), 3)));
  return degs;
}

/********************************************************************/
/*                    1rst part: find the field K                   */
/********************************************************************/
static GEN AllStark(GEN data, GEN nf, long flag, long prec);

/* Columns of C [HNF] give the generators of a subgroup of the finite abelian
 * group A [ in terms of implicit generators ], compute data to work in A/C:
 * 1) order
 * 2) structure
 * 3) the matrix A ->> A/C
 * 4) the group C */
static GEN
InitQuotient(GEN C)
{
  GEN z, U, D = smithall(C, &U, NULL);
  z = cgetg(5, t_VEC);
  gel(z,1) = dethnf_i(D);
  gel(z,2) = mattodiagonal_i(D);
  gel(z,3) = U;
  gel(z,4) = C; return z;
}

/* Let s: A -> B given by P, and let DA, DB be resp. the matrix of the
   relations of A and B, compute the kernel of s. If DA = 0 then A is free */
static GEN
ComputeKernel0(GEN P, GEN DA, GEN DB)
{
  pari_sp av = avma;
  long nbA = lg(DA)-1, rk;
  GEN U;

  rk = nbA + lg(DB) - lg(hnfall_i(shallowconcat(P, DB), &U, 1));
  U = vecslice(U, 1,rk);
  U = rowslice(U, 1,nbA);
  if (!gcmp0(DA)) U = shallowconcat(U, DA);
  return gerepileupto(av, hnf(U));
}

/* Let m and n be two moduli such that n|m and let C be a congruence
   group modulo n, compute the corresponding congruence group modulo m
   ie the kernel of the map Clk(m) ->> Clk(n)/C */
static GEN
ComputeKernel(GEN bnrm, GEN bnrn, GEN dtQ)
{
  long i, nbm;
  pari_sp av = avma;
  GEN Mrm, genm, Mrq, mgq, P;

  Mrm  = diagonal_i(gmael(bnrm, 5, 2));
  Mrq  = diagonal_i(gel(dtQ,2));
  genm = gmael(bnrm, 5, 3);
  nbm  = lg(genm) - 1;
  mgq  = gel(dtQ,3);

  P = cgetg(nbm + 1, t_MAT);
  for (i = 1; i <= nbm; i++)
    gel(P,i) = gmul(mgq, isprincipalray(bnrn, gel(genm,i)));

  return gerepileupto(av, ComputeKernel0(P, Mrm, Mrq));
}

/* Let C a congruence group in bnr, compute its subgroups of index 2 as
   subgroups of Clk(bnr) */
static GEN
ComputeIndex2Subgroup(GEN bnr, GEN C)
{
  pari_sp av = avma;
  long nb, i;
  GEN D, Mr, U, T, subgrp;

  disable_dbg(0);

  Mr = diagonal_i(gmael(bnr, 5, 2));
  D = smithall(hnf_gauss(C, Mr), &U, NULL);
  T = gmul(C,ginv(U));
  subgrp  = subgrouplist(D, mkvec(gen_2));
  nb = lg(subgrp);
  for (i = 1; i < nb; i++)
    gel(subgrp,i) = hnf(shallowconcat(gmul(T, gel(subgrp,i)), Mr));

  disable_dbg(-1);
  return gerepilecopy(av, subgrp);
}

static GEN
Order(GEN cyc, GEN x)
{
  pari_sp av = avma;
  long i, l = lg(cyc);
  GEN c,o,f = gen_1;
  for (i = 1; i < l; i++)
  {
    o = gel(cyc,i);
    c = gcdii(o, gel(x,i));
    if (!is_pm1(c)) o = diviiexact(o,c);
    f = lcmii(f, o);
  }
  return gerepileuptoint(av, f);
}

/* Let pr be a prime (pr may divide mod(bnr)), compute the indexes
   e,f of the splitting of pr in the class field nf(bnr/subgroup) */
static GEN
GetIndex(GEN pr, GEN bnr, GEN subgroup)
{
  long v, e, f;
  pari_sp av = avma;
  GEN bnf, mod, mod0, bnrpr, subpr, M, dtQ, p1;
  GEN rep, cycpr, cycQ;

  bnf  = gel(bnr,1);
  mod  = gmael(bnr, 2, 1);
  mod0 = gel(mod,1);

  v = idealval(bnf, mod0, pr);
  if (v == 0)
  {
    bnrpr = bnr;
    subpr = subgroup;
    e = 1;
  }
  else
  {
    GEN mpr = cgetg(3, t_VEC);
    GEN mpr0 = idealdivpowprime(bnf, mod0, pr, utoipos(v));
    gel(mpr,1) = mpr0; /* part of mod coprime to pr */
    mpr[2] = mod[2];
    bnrpr = buchrayinitgen(bnf, mpr);
    cycpr = gmael(bnrpr, 5, 2);
    M = gmul(bnrGetSurj(bnr, bnrpr), subgroup);
    subpr = hnf(shallowconcat(M, diagonal_i(cycpr)));
    /* e = #(bnr/subgroup) / #(bnrpr/subpr) */
    e = itos( diviiexact(dethnf_i(subgroup), dethnf_i(subpr)) );
  }

  /* f = order of [pr] in bnrpr/subpr */
  dtQ  = InitQuotient(subpr);
  p1   = gmul(gel(dtQ,3), isprincipalray(bnrpr, pr));
  cycQ = gel(dtQ,2);
  f  = itos( Order(cycQ, p1) );
  avma = av;
  rep = cgetg(3, t_VECSMALL);
  rep[1] = e;
  rep[2] = f; return rep;
}

static GEN get_listCR(GEN bnr, GEN dtQ);
static GEN InitChar(GEN bnr, GEN listCR, long prec);

/* Given a conductor and a subgroups, return the corresponding
   complexity and precision required using quickpol. Fill data[5] with
   listCR */
static long
CplxModulus(GEN data, long *newprec, long prec)
{
  long pr, ex, dprec = DEFAULTPREC;
  pari_sp av;
  GEN pol, listCR, cpl, bnr = gel(data,1), nf = checknf(bnr);

  listCR = get_listCR(bnr, gel(data,3));
  for (av = avma;; avma = av)
  {
    gel(data,5) = InitChar(bnr, listCR, dprec);
    pol = AllStark(data, nf, -1, dprec);
    pr = 1 + (gexpo(pol)>>TWOPOTBITS_IN_LONG);
    if (pr < 0) pr = 0;
    dprec = max(dprec, pr) + EXTRA_PREC;
    if (!gcmp0(leading_term(pol)))
    {
      cpl = QuickNormL2(pol, DEFAULTPREC);
      if (!gcmp0(cpl)) break;
    }
    if (DEBUGLEVEL>1) pari_warn(warnprec, "CplxModulus", dprec);
  }
  ex = gexpo(cpl); avma = av;
  if (DEBUGLEVEL>1) fprintferr("cpl = 2^%ld\n", ex);

  gel(data,5) = listCR;
  *newprec = dprec; return ex;
}

/* Let f be a conductor without infinite part and let C be a
   congruence group modulo f, compute (m,D) such that D is a
   congruence group of conductor m where m is a multiple of f
   divisible by all the infinite places but one, D is a subgroup of
   index 2 of Im(C) in Clk(m), no prime dividing f splits in the
   corresponding quadratic extension and m is of minimal norm. Return
   bnr(m), D, quotient Ck(m) / D and Clk(m) / C */
static GEN
FindModulus(GEN bnr, GEN dtQ, long *newprec, long prec)
{
  const long limnorm = 400;
  long n, i, narch, nbp, maxnorm, minnorm, N, nbidnn, s, c, j, nbcand;
  long first = 1, pr, rb, oldcpl = -1, iscyc = 0;
  pari_sp av = avma, av1;
  GEN rep, bnf, nf, f, arch, m, listid, idnormn, bnrm, ImC;
  GEN candD, bpr, indpr, sgp, p1, p2;

  sgp = gel(dtQ,4);
  bnf = gel(bnr,1);
  nf  = gel(bnf,7);
  N   = degpol(nf[1]);
  f   = gmael3(bnr, 2, 1, 1);

  rep = NULL;

  /* if cpl < rb, it is not necessary to try another modulus */
  rb = expi( powgi(gmul(gel(nf,3), det(f)), gmul2n(gmael(bnr, 5, 1), 3)) );

  bpr = divcond(bnr);
  nbp = lg(bpr) - 1;

  indpr = cgetg(nbp + 1,t_VECSMALL);
  for (i = 1; i <= nbp; i++)
  {
    p1 = GetIndex(gel(bpr,i), bnr, sgp);
    indpr[i] = p1[1] * p1[2];
  }

  /* Initialization of the possible infinite part */
  arch = const_vec(N, gen_1);

  /* narch = (N == 2)? 1: N; -- if N=2, only one case is necessary */
  narch = N;
  m = mkvec2(NULL, arch);

  /* go from minnorm up to maxnorm. If necessary, increase these values.
   * If we cannot find a suitable conductor of norm < limnorm, stop */
  maxnorm = 50;
  minnorm = 1;

  /* if the extension is cyclic then we _must_ find a suitable conductor */
  if (lg(dtQ[2]) == 2) iscyc = 1;

  if (DEBUGLEVEL>1)
    fprintferr("Looking for a modulus of norm: ");

  for(;;)
  {
    disable_dbg(0);
    listid = ideallist(nf, maxnorm); /* all ideals of norm <= maxnorm */
    disable_dbg(-1);

    av1 = avma;
    for (n = minnorm; n <= maxnorm; n++)
    {
      if (DEBUGLEVEL>1) fprintferr(" %ld", n);
      avma = av1;

      idnormn = gel(listid,n);
      nbidnn  = lg(idnormn) - 1;
      for (i = 1; i <= nbidnn; i++)
      { /* finite part of the conductor */
	gel(m,1) = idealmul(nf, f, gel(idnormn,i));

	for (s = 1; s <= narch; s++)
	{ /* infinite part */
	  gel(arch,N+1-s) = gen_0;

          /* compute Clk(m), check if m is a conductor */
	  disable_dbg(0);
	  bnrm = buchrayinitgen(bnf, m);
	  p1   = conductor(bnrm, NULL, -1);
	  disable_dbg(-1);
          gel(arch,N+1-s) = gen_1;
	  if (!signe(p1)) continue;

          /* compute Im(C) in Clk(m)... */
          ImC = ComputeKernel(bnrm, bnr, dtQ);

          /* ... and its subgroups of index 2 */
          candD  = ComputeIndex2Subgroup(bnrm, ImC);
          nbcand = lg(candD) - 1;
          for (c = 1; c <= nbcand; c++)
          {
            GEN D  = gel(candD,c);
            long cpl;

            /* check if m is the conductor */
            p1 = conductor(bnrm, D, -1);
            if (!signe(p1)) continue;

            /* check the splitting of primes */
            for (j = 1; j <= nbp; j++)
            {
              p1 = GetIndex(gel(bpr,j), bnrm, D);
              if (p1[1] * p1[2] == indpr[j]) break; /* no good */
            }
            if (j <= nbp) continue;

            p2 = cgetg(6, t_VEC); /* p2[5] filled in CplxModulus */
            gel(p2,1) = bnrm;
            gel(p2,2) = D;
            gel(p2,3) = InitQuotient(D);
            gel(p2,4) = InitQuotient(ImC);
            if (DEBUGLEVEL>1)
              fprintferr("\nTrying modulus = %Z and subgroup = %Z\n",
	                 gmael(bnrm, 2, 1), D);
            cpl = CplxModulus(p2, &pr, prec);
            if (oldcpl < 0 || cpl < oldcpl)
            {
              *newprec = pr;
              if (rep) gunclone(rep);
              rep    = gclone(p2);
              oldcpl = cpl;
            }
            if (oldcpl < rb) goto END; /* OK */

            if (DEBUGLEVEL>1) fprintferr("Trying to find another modulus...");
            first = 0;
          }
	}
        if (!first) goto END; /* OK */
      }
    }
    /* if necessary compute more ideals */
    minnorm = maxnorm;
    maxnorm <<= 1;
    if (!iscyc && maxnorm > limnorm) return NULL;

  }
END:
  if (DEBUGLEVEL>1)
    fprintferr("No, we're done!\nModulus = %Z and subgroup = %Z\n",
               gmael3(rep, 1, 2, 1), gel(rep,2));
  gel(rep,5) = InitChar(gel(rep,1), gel(rep,5), *newprec);
  return gerepilecopy(av, rep);
}

/********************************************************************/
/*                      2nd part: compute W(X)                      */
/********************************************************************/

/* compute the list of W(chi) such that Ld(s,chi) = W(chi) Ld(1 - s, chi*),
 * for all chi in LCHI. All chi have the same conductor (= cond(bnr)).
 * if check == 0 do not check the result */
static GEN
ArtinNumber(GEN bnr, GEN LCHI, long check, long prec)
{
  long ic, i, j, nz, N, nChar = lg(LCHI)-1;
  pari_sp av = avma, av2, lim;
  GEN sqrtnc, dc, cond, condZ, cond0, cond1, lambda, nf, T;
  GEN cyc, vN, vB, diff, vt, idg, mu, idh, zid, gen, z, nchi;
  GEN indW, W, classe, s0, s, den, muslambda, beta2, sarch;
  CHI_t **lC;
  GROUP_t G;

  lC = (CHI_t**)new_chunk(nChar + 1);
  indW = cgetg(nChar + 1, t_VECSMALL);
  W = cgetg(nChar + 1, t_VEC);
  for (ic = 0, i = 1; i <= nChar; i++)
  {
    GEN CHI = gel(LCHI,i);
    if (cmpui(2, gel(CHI,3)) >= 0) { gel(W,i) = gen_1; continue; } /* trivial case */
    ic++; indW[ic] = i;
    lC[ic] = (CHI_t*)new_chunk(sizeof(CHI_t));
    init_CHI_C(lC[ic], CHI);
  }
  if (!ic) return W;
  nChar = ic;

  nf    = gmael(bnr, 1, 7);
  diff  = gmael(nf, 5, 5);
  T     = gmael(nf, 5, 4);
  cond  = gmael(bnr, 2, 1);
  cond0 = gel(cond,1); condZ = gcoeff(cond0,1,1);
  cond1 = arch_to_perm(gel(cond,2));
  N     = degpol(nf[1]);

  sqrtnc  = gsqrt(idealnorm(nf, cond0), prec);
  dc  = idealmul(nf, diff, cond0);
  den = idealnorm(nf, dc);
  z   = InitRU(den, prec);

  /* compute a system of elements congru to 1 mod cond0 and giving all
     possible signatures for cond1 */
  sarch = zarchstar(nf, cond0, cond1);

  /* find lambda in diff.cond such that gcd(lambda.(diff.cond)^-1,cond0) = 1
     and lambda >> 0 at cond1 */
  lambda = idealappr(nf, dc);
  lambda = set_sign_mod_idele(nf, NULL, lambda, cond,sarch);
  idg = idealdivexact(nf, lambda, dc);

  /* find mu in idg such that idh=(mu) / idg is coprime with cond0 and
     mu >> 0 at cond1 */
  if (!gcmp1(gcoeff(idg, 1, 1)))
  {
    GEN P = divcond(bnr);
    GEN f = concat_factor(idealfactor(nf, idg),
                          mkmat2(P, zerocol(lg(P)-1)));

    mu = set_sign_mod_idele(nf, NULL, idealapprfact(nf, f), cond,sarch);
    idh = idealdivexact(nf, mu, idg);
  }
  else
  {
    mu  = gen_1;
    idh = idg;
  }

  muslambda = gmul(den, element_div(nf, mu, lambda));

  /* compute a system of generators of (Ok/cond)^* cond1-positive */
  zid = zidealstarinitgen(nf, cond0);
  cyc = gmael(zid, 2, 2);
  gen = gmael(zid, 2, 3);
  nz = lg(gen) - 1;

  nchi = cgetg(nChar+1, t_VEC);
  for (ic = 1; ic <= nChar; ic++) gel(nchi,ic) = cgetg(nz + 1, t_VECSMALL);

  for (i = 1; i <= nz; i++)
  {
    if (is_bigint(cyc[i]))
      pari_err(talker,"conductor too large in ArtinNumber");
    gel(gen,i) = set_sign_mod_idele(nf, NULL, gel(gen,i), cond,sarch);
    classe = isprincipalray(bnr, gel(gen,i));
    for (ic = 1; ic <= nChar; ic++) {
      GEN n = gel(nchi,ic);
      n[i] = (long)EvalChar_n(lC[ic], classe);
    }
  }

  /* Sum chi(beta) * exp(2i * Pi * Tr(beta * mu / lambda) where beta
     runs through the classes of (Ok/cond0)^* and beta cond1-positive */

  vt = cgetg(N + 1, t_VEC); /* Tr(w_i) */
  for (i = 1; i <= N; i++) gel(vt,i) = gcoeff(T,i,1);

  G.cyc = gtovecsmall(cyc);
  G.r = nz;
  G.j = const_vecsmall(nz, 0);

  vN = cgetg(nChar+1, t_VEC);
  for (ic = 1; ic <= nChar; ic++) gel(vN,ic) = const_vecsmall(nz, 0);

  av2 = avma; lim = stack_lim(av2, 1);
  vB = const_vec(nz, gen_1);

  s0 = powgi(z, Rg_to_Fp(gmul(vt, muslambda), den)); /* for beta = 1 */
  s = const_vec(nChar, s0);

  while ( (i = NextElt(&G)) )
  {
    gel(vB,i) = FpC_red(element_muli(nf, gel(vB,i), gel(gen,i)), condZ);
    for (j=1; j<i; j++) vB[j] = vB[i];

    for (ic = 1; ic <= nChar; ic++)
    {
      GEN v = gel(vN,ic), n = gel(nchi,ic);
      v[i] = (long)Fl_add(v[i], n[i], lC[ic]->ord);
      for (j=1; j<i; j++) v[j] = v[i];
    }

    gel(vB,i) = set_sign_mod_idele(nf, NULL, gel(vB,i), cond,sarch);
    beta2 = element_mul(nf, gel(vB,i), muslambda);
    
    s0 = powgi(z, Rg_to_Fp(gmul(vt, beta2), den));
    for (ic = 1; ic <= nChar; ic++)
    {
      GEN n = gel(vN,ic), val = lC[ic]->val[ n[i] ];
      gel(s,ic) = gadd(gel(s,ic), gmul(val, s0));
    }

    if (low_stack(lim, stack_lim(av2, 1)))
    {
      if (DEBUGMEM > 1) pari_warn(warnmem,"ArtinNumber");
      gerepileall(av2, 2, &s, &vB);
    }
  }

  classe = isprincipalray(bnr, idh);
  z = gpowgs(gneg_i(gi), lg(cond1)-1);

  for (ic = 1; ic <= nChar; ic++)
  {
    s0 = gmul(gel(s,ic), EvalChar(lC[ic], classe));
    s0 = gdiv(s0, sqrtnc);
    if (check && - expo(subrs(gnorm(s0), 1)) < bit_accuracy(prec) >> 1)
      pari_err(bugparier, "ArtinNumber");
    gel(W, indW[ic]) = gmul(s0, z);
  }
  return gerepilecopy(av, W);
}

static GEN
ComputeAllArtinNumbers(GEN dataCR, GEN vChar, int check, long prec)
{
  long j, k, cl = lg(dataCR) - 1, J = lg(vChar)-1;
  GEN W = cgetg(cl+1,t_VEC), WbyCond, LCHI;

  for (j = 1; j <= J; j++)
  {
    GEN LChar = gel(vChar,j), ldata = vecpermute(dataCR, LChar);
    GEN dtcr = gel(ldata,1), bnr = ch_bnr(dtcr);
    long l = lg(LChar);

    if (DEBUGLEVEL>1)
      fprintferr("* Root Number: cond. no %ld/%ld (%ld chars)\n", j, J, l-1);
    LCHI = cgetg(l, t_VEC);
    for (k = 1; k < l; k++) gel(LCHI,k) = ch_CHI0(gel(ldata,k));
    WbyCond = ArtinNumber(bnr, LCHI, check, prec);
    for (k = 1; k < l; k++) W[LChar[k]] = WbyCond[k];
  }
  return W;
}
static GEN
SingleArtinNumber(GEN bnr, GEN chi, long prec)
{ return (GEN)ArtinNumber(bnr, mkvec(chi), 1, prec)[1]; }

/* compute the constant W of the functional equation of
   Lambda(chi). If flag = 1 then chi is assumed to be primitive */
GEN
bnrrootnumber(GEN bnr, GEN chi, long flag, long prec)
{
  long l;
  pari_sp av = avma;
  GEN cond, condc, bnrc, CHI, cyc;

  if (flag < 0 || flag > 1) pari_err(flagerr,"bnrrootnumber");

  checkbnr(bnr);
  cyc = gmael(bnr, 5, 2);
  cond = gmael(bnr, 2, 1);
  l    = lg(cyc);

  if (typ(chi) != t_VEC || lg(chi) != l)
    pari_err(talker, "incorrect character in bnrrootnumber");

  if (flag) condc = NULL;
  else
  {
    condc = bnrconductorofchar(bnr, chi);
    if (gequal(cond, condc)) flag = 1;
  }

  if (flag)
  {
    GEN initc = init_get_chic(cyc);
    bnrc = bnr;
    CHI = get_Char(chi, initc, NULL, prec);
  }
  else
  {
    bnrc = buchrayinitgen(gel(bnr,1), condc);
    CHI = GetPrimChar(chi, bnr, bnrc, prec);
  }
  return gerepilecopy(av, SingleArtinNumber(bnrc, CHI, prec));
}

/********************************************************************/
/*               3rd part: initialize the characters                */
/********************************************************************/

static GEN
LiftChar(GEN cyc, GEN Mat, GEN chi, GEN D)
{
  long lm = lg(cyc), l  = lg(chi), i, j;
  GEN lchi = cgetg(lm, t_VEC);
  for (i = 1; i < lm; i++)
  {
    pari_sp av = avma;
    GEN t, s  = mulii(gel(chi,1), gcoeff(Mat, 1, i));
    for (j = 2; j < l; j++)
    { /* rarely exercised: D[1]/D[j] could be precomputed */
      t = mulii(gel(chi,j), diviiexact(gel(D,1), gel(D,j)));
      s = addii(s, mulii(t, gcoeff(Mat, j, i)));
    }
    t = diviiexact(mulii(s, gel(cyc,i)), gel(D,1));
    gel(lchi,i) = gerepileuptoint(av, modii(t, gel(cyc,i)));
  }
  return lchi;
}

/* Let chi be a character, A(chi) corresponding to the primes dividing diff
   at s = flag. If s = 0, returns [r, A] where r is the order of vanishing
   at s = 0 corresponding to diff. No GC */
static GEN
ComputeAChi(GEN dtcr, long *r, long flag, long prec)
{
  long l, i;
  GEN p1, A, diff, chi, bnrc;

  bnrc = ch_bnr(dtcr);
  diff = ch_diff(dtcr); l = lg(diff);
  chi  = ch_CHI0(dtcr);

  A = gen_1;
  *r = 0;
  for (i = 1; i < l; i++)
  {
    GEN pr = gel(diff,i), B;
    p1  = ComputeImagebyChar(chi, isprincipalray(bnrc, pr));

    if (flag)
      B = gsub(gen_1, gdiv(p1, pr_norm(pr)));
    else if (gcmp1(p1))
    {
      B = glog(pr_norm(pr), prec);
      (*r)++;
    }
    else
      B = gsub(gen_1, p1);
    A = gmul(A, B);
  }
  return A;
}

static GEN
_data4(GEN arch, long r1, long r2)
{
  GEN z = cgetg(5, t_VECSMALL);
  long i, b, q = 0;

  for (i=1; i<=r1; i++) if (signe(arch[i])) q++;
  z[1] = q; b = r1 - q;
  z[2] = b;
  z[3] = r2;
  z[4] = max(b+r2+1, r2+q);
  return z;
}

/* Given a list [chi, F = cond(chi)] of characters over Cl(bnr), compute a
   vector dataCR containing for each character:
   2: the constant C(F) [t_REAL]
   3: bnr(F)
   4: [q, r1 - q, r2, rc] where
        q = number of real places in F
        rc = max{r1 + r2 - q + 1, r2 + q}
   6: diff(chi) primes dividing m but not F
   7: finite part of F

   1: chi
   5: [(c_i), z, d] in bnr(m)
   8: [(c_i), z, d] in bnr(F) */
static GEN
InitChar(GEN bnr, GEN listCR, long prec)
{
  GEN bnf = checkbnf(bnr), nf = checknf(bnf);
  GEN modul, dk, C, dataCR, chi, cond, Mr, initc;
  long N, r1, r2, prec2, i, j, l;
  pari_sp av = avma;

  modul = gmael(bnr, 2, 1);
  Mr    = gmael(bnr, 5, 2);
  dk    = gel(nf,3);
  N     = degpol(nf[1]);
  nf_get_sign(nf, &r1,&r2);
  prec2 = ((prec-2) << 1) + EXTRA_PREC;
  C     = gmul2n(sqrtr_abs(divir(dk, gpowgs(mppi(prec2),N))), -r2);
  initc = init_get_chic(Mr);

  disable_dbg(0);

  l = lg(listCR); dataCR = cgetg(l, t_VEC);
  for (i = 1; i < l; i++)
  {
    GEN olddtcr, dtcr = cgetg(9, t_VEC);
    gel(dataCR,i) = dtcr;

    chi  = gmael(listCR, i, 1);
    cond = gmael(listCR, i, 2);

    /* do we already know the invariants of chi? */
    olddtcr = NULL;
    for (j = 1; j < i; j++)
      if (gequal(cond, gmael(listCR,j,2))) { olddtcr = gel(dataCR,j); break; }

    if (!olddtcr)
    {
      ch_C(dtcr) = gmul(C, gsqrt(det(gel(cond,1)), prec2));
      ch_4(dtcr) = _data4(gel(cond,2),r1,r2);
      ch_cond(dtcr) = gel(cond,1);
      if (gequal(cond,modul))
      {
        ch_bnr(dtcr) = bnr;
        ch_diff(dtcr) = cgetg(1, t_VEC);
      }
      else
      {
        ch_bnr(dtcr) = buchrayinitgen(bnf, cond);
        ch_diff(dtcr) = get_prdiff(bnr, cond);
      }
    }
    else
    {
      ch_C(dtcr) = ch_C(olddtcr);
      ch_bnr(dtcr) = ch_bnr(olddtcr);
      ch_4(dtcr) = ch_4(olddtcr);
      ch_diff(dtcr) = ch_diff(olddtcr);
      ch_cond(dtcr) = ch_cond(olddtcr);
    }

    ch_chi(dtcr) = chi; /* the character */
    ch_CHI(dtcr) = get_Char(chi,initc,NULL,prec); /* associated to bnr(m) */
    chi = GetPrimChar(chi, bnr, ch_bnr(dtcr), prec2);
    if (!chi) chi = ch_CHI(dtcr);
    ch_CHI0(dtcr) = chi;
  }

  disable_dbg(-1);
  return gerepilecopy(av, dataCR);
}

/* compute the list of characters to consider for AllStark and
   initialize precision-independent data to compute with them */
static GEN
get_listCR(GEN bnr, GEN dtQ)
{
  GEN MrD, listCR, vecchi, lchi, Surj, cond, Mr, d, allCR;
  long hD, h, nc, i, j, tnc;

  Surj = gel(dtQ,3);
  MrD  = gel(dtQ,2);
  Mr   = gmael(bnr, 5, 2);
  hD   = itos(gel(dtQ,1));
  h    = hD >> 1;

  disable_dbg(0);

  listCR = cgetg(h + 1, t_VEC); /* non-conjugate characters */
  nc  = 1;
  allCR  = cgetg(h + 1, t_VEC); /* all characters, including conjugates */
  tnc = 1;

  vecchi = EltsOfGroup(hD, MrD);

  for (i = 1; tnc <= h; i++)
  {
    /* lift a character of D in Clk(m) */
    lchi = LiftChar(Mr, Surj, gel(vecchi,i), MrD);

    for (j = 1; j < tnc; j++)
      if (gequal(lchi, gel(allCR,j))) break;
    if (j != tnc) continue;

    cond = bnrconductorofchar(bnr, lchi);
    if (gcmp0(gel(cond,2))) continue;

    /* the infinite part of chi is non trivial */
    gel(listCR,nc++) = mkvec2(lchi, cond);
    gel(allCR,tnc++) = lchi;

    /* if chi is not real, add its conjugate character to allCR */
    d = Order(Mr, lchi);
    if (!equaliu(d, 2))
      gel(allCR,tnc++) = ConjChar(lchi, Mr);
  }
  disable_dbg(-1);
  setlg(listCR, nc); return listCR;
}

/* recompute dataCR with the new precision */
static GEN
CharNewPrec(GEN dataCR, GEN nf, long prec)
{
  GEN dk, C, p1;
  long N, l, j, prec2;

  dk    =  gel(nf,3);
  N     =  degpol(gel(nf,1));
  prec2 = ((prec - 2)<<1) + EXTRA_PREC;

  C = sqrtr(gdiv(absi(dk), gpowgs(mppi(prec2), N)));

  l = lg(dataCR);
  for (j = 1; j < l; j++)
  {
    GEN dtcr = gel(dataCR,j);
    ch_C(dtcr) = gmul(C, gsqrt(dethnf_i(ch_cond(dtcr)), prec2));

    gmael(ch_bnr(dtcr), 1, 7) = nf;

    p1 = ch_CHI( dtcr); gel(p1,2) = InitRU(gel(p1,3), prec2);
    p1 = ch_CHI0(dtcr); gel(p1,2) = InitRU(gel(p1,3), prec2);
  }

  return dataCR;
}

/********************************************************************/
/*             4th part: compute the coefficients an(chi)           */
/*                                                                  */
/* matan entries are arrays of ints containing the coefficients of  */
/* an(chi) as a polmod modulo cyclo(order(chi))                     */
/********************************************************************/

static void
_0toCoeff(int *rep, long deg)
{
  long i;
  for (i=0; i<deg; i++) rep[i] = 0;
}

/* transform a polmod into Coeff */
static void
Polmod2Coeff(int *rep, GEN polmod, long deg)
{
  long i;
  if (typ(polmod) == t_POLMOD)
  {
    GEN pol = gel(polmod,2);
    long d = degpol(pol);

    pol += 2;
    for (i=0; i<=d; i++) rep[i] = itos(gel(pol,i));
    for (   ; i<deg; i++) rep[i] = 0;
  }
  else
  {
    rep[0] = itos(polmod);
    for (i=1; i<deg; i++) rep[i] = 0;
  }
}

/* initialize a deg * n matrix of ints */
static int**
InitMatAn(long n, long deg, long flag)
{
  long i, j;
  int *a, **A = (int**)gpmalloc((n+1)*sizeof(int*));
  A[0] = NULL;
  for (i = 1; i <= n; i++)
  {
    a = (int*)gpmalloc(deg*sizeof(int));
    A[i] = a; a[0] = (i == 1 || flag);
    for (j = 1; j < deg; j++) a[j] = 0;
  }
  return A;
}

static void
FreeMat(int **A, long n)
{
  long i;
  for (i = 0; i <= n; i++)
    if (A[i]) free((void*)A[i]);
  free((void*)A);
}

/* initialize Coeff reduction */
static int**
InitReduction(GEN CHI, long deg)
{
  long j;
  pari_sp av = avma;
  int **A;
  GEN d, polmod, pol;

  A   = (int**)gpmalloc(deg*sizeof(int*));
  d   = gel(CHI,3);
  pol = cyclo(itos(d), 0);
  for (j = 0; j < deg; j++)
  {
    A[j] = (int*)gpmalloc(deg*sizeof(int));
    polmod = gmodulo(monomial(gen_1, deg+j, 0), pol);
    Polmod2Coeff(A[j], polmod, deg);
  }

  avma = av; return A;
}

#if 0
void
pan(int **an, long n, long deg)
{
  long i,j;
  for (i = 1; i <= n; i++)
  {
    fprintferr("n = %ld: ",i);
    for (j = 0; j < deg; j++) fprintferr("%d ",an[i][j]);
    fprintferr("\n");
  }
}
#endif

/* returns 0 if c is zero, 1 otherwise. */
static int
IsZero(int* c, long deg)
{
  long i;
  for (i = 0; i < deg; i++)
    if (c[i]) return 0;
  return 1;
}

/* set c0 <-- c0 * c1 */
static void
MulCoeff(int *c0, int* c1, int** reduc, long deg)
{
  long i,j;
  int c, *T;

  if (IsZero(c0,deg)) return;

  T = (int*)new_chunk(2*deg);
  for (i = 0; i < 2*deg; i++)
  {
    c = 0;
    for (j = 0; j <= i; j++)
      if (j < deg && j > i - deg) c += c0[j] * c1[i-j];
    T[i] = c;
  }
  for (i = 0; i < deg; i++)
  {
    c = T[i];
    for (j = 0; j < deg; j++) c += reduc[j][i] * T[deg+j];
    c0[i] = c;
  }
}

/* c0 <- c0 + c1 * c2 */
static void
AddMulCoeff(int *c0, int *c1, int* c2, int** reduc, long deg)
{
  long i, j;
  pari_sp av;
  int c, *t;

  if (IsZero(c2,deg)) return;
  if (!c1) /* c1 == 1 */
  {
    for (i = 0; i < deg; i++) c0[i] += c2[i];
    return;
  }
  av = avma;
  t = (int*)new_chunk(2*deg); /* = c1 * c2, not reduced */
  for (i = 0; i < 2*deg; i++)
  {
    c = 0;
    for (j = 0; j <= i; j++)
      if (j < deg && j > i - deg) c += c1[j] * c2[i-j];
    t[i] = c;
  }
  for (i = 0; i < deg; i++)
  {
    c = t[i];
    for (j = 0; j < deg; j++) c += reduc[j][i] * t[deg+j];
    c0[i] += c;
  }
  avma = av;
}

/* evaluate the Coeff. No Garbage collector */
static GEN
EvalCoeff(GEN z, int* c, long deg)
{
  long i,j;
  GEN e, r;

  if (!c) return gen_0;
#if 0
  /* standard Horner */
  e = stoi(c[deg - 1]);
  for (i = deg - 2; i >= 0; i--)
    e = gadd(stoi(c[i]), gmul(z, e));
#else
  /* specific attention to sparse polynomials */
  e = NULL;
  for (i = deg-1; i >=0; i=j-1)
  {
    for (j=i; c[j] == 0; j--)
      if (j==0)
      {
        if (!e) return NULL;
        if (i!=j) z = gpowgs(z,i-j+1);
        return gmul(e,z);
      }
    if (e)
    {
      r = (i==j)? z: gpowgs(z,i-j+1);
      e = gadd(gmul(e,r), stoi(c[j]));
    }
    else
      e = stoi(c[j]);
  }
#endif
  return e;
}

/* copy the n * (m+1) array matan */
static void
CopyCoeff(int** a, int** a2, long n, long m)
{
  long i,j;

  for (i = 1; i <= n; i++)
  {
    int *b = a[i], *b2 = a2[i];
    for (j = 0; j < m; j++) b2[j] = b[j];
  }
}

/* return q*p if <= n. Beware overflow */
static long
next_pow(long q, long p, long n)
{
  const GEN x = muluu((ulong)q, (ulong)p);
  const ulong qp = (ulong)x[2];
  return (lgefint(x) > 3 || qp > (ulong)n)? 0: qp;
}

static void
an_AddMul(int **an,int **an2, long np, long n, long deg, GEN chi, int **reduc)
{
  GEN chi2 = chi;
  long q, qk, k;
  int *c, *c2 = (int*)new_chunk(deg);

  CopyCoeff(an, an2, n/np, deg);
  for (q=np;;)
  {
    if (gcmp1(chi2)) c = NULL; else { Polmod2Coeff(c2, chi2, deg); c = c2; }
    for(k = 1, qk = q; qk <= n; k++, qk += q)
      AddMulCoeff(an[qk], c, an2[k], reduc, deg);
    if (! (q = next_pow(q,np, n)) ) break;

    chi2 = gmul(chi2, chi);
  }
}

/* correct the coefficients an(chi) according with diff(chi) in place */
static void
CorrectCoeff(GEN dtcr, int** an, int** reduc, long n, long deg)
{
  pari_sp av = avma;
  long lg, j, np;
  pari_sp av1;
  int **an2;
  GEN bnrc, diff, chi, pr;
  CHI_t C;

  diff = ch_diff(dtcr); lg = lg(diff) - 1;
  if (!lg) return;

  if (DEBUGLEVEL>2) fprintferr("diff(CHI) = %Z", diff);
  bnrc =  ch_bnr(dtcr);
  init_CHI_alg(&C, ch_CHI0(dtcr));

  an2 = InitMatAn(n, deg, 0);
  av1 = avma;
  for (j = 1; j <= lg; j++)
  {
    pr = gel(diff,j);
    np = itos( pr_norm(pr) );

    chi  = EvalChar(&C, isprincipalray(bnrc, pr));

    an_AddMul(an,an2,np,n,deg,chi,reduc);
    avma = av1;
  }
  FreeMat(an2, n); avma = av;
}

/* compute the coefficients an in the general case */
static int**
ComputeCoeff(GEN dtcr, LISTray *R, long n, long deg)
{
  pari_sp av = avma, av2;
  long i, l, np;
  int **an, **reduc, **an2;
  GEN L, CHI, chi;
  CHI_t C;

  CHI = ch_CHI(dtcr); init_CHI_alg(&C, CHI);

  an  = InitMatAn(n, deg, 0);
  an2 = InitMatAn(n, deg, 0);
  reduc  = InitReduction(CHI, deg);
  av2 = avma;

  L = R->L1; l = lg(L);
  for (i=1; i<l; i++, avma = av2)
  {
    np = L[i]; 
    chi  = EvalChar(&C, R->L1ray[i]);
    an_AddMul(an,an2,np,n,deg,chi,reduc);
  }
  FreeMat(an2, n);

  CorrectCoeff(dtcr, an, reduc, n, deg);
  FreeMat(reduc, deg-1);
  avma = av; return an;
}

/********************************************************************/
/*              5th part: compute L-functions at s=1                */
/********************************************************************/
static void
deg11(LISTray *R, long p, GEN bnr, GEN pr) {
  GEN z = isprincipalray(bnr, pr);
  appendL(R->L1, (GEN)p);
  appendL((GEN)R->L1ray, z);
}
static void
deg12(LISTray *R, long p, GEN bnr, GEN Lpr) {
  GEN z = isprincipalray(bnr, gel(Lpr,1));
  appendL(R->L11, (GEN)p);
  appendL((GEN)R->L11ray, z);
}
static void
deg0(LISTray *R, long p) {
  appendL(R->L0, (GEN)p);
}
static void
deg2(LISTray *R, long p) {
  appendL(R->L2, (GEN)p);
}

/* pi(x) <= ?? */
static long
PiBound(long x)
{
  double lx = log((double)x);
  return 1 + (long) (x/lx * (1 + 3/(2*lx)));
}

static void
InitPrimesQuad(GEN bnr, long N0, LISTray *R)
{
  pari_sp av = avma;
  GEN bnf = gel(bnr,1), cond = gmael3(bnr,2,1,1);
  long p,i,l, condZ = itos(gcoeff(cond,1,1)), contZ = itos(content(cond));
  GEN prime, pr, nf = checknf(bnf), dk = gel(nf,3);
  byteptr d = diffptr + 1;
  GEN *gptr[7];

  l = 1 + PiBound(N0);
  R->L0 = cget1(l, t_VECSMALL);
  R->L2 = cget1(l, t_VECSMALL); R->condZ = condZ;
  R->L1 = cget1(l, t_VECSMALL); R->L1ray = (GEN*)cget1(l, t_VEC);
  R->L11= cget1(l, t_VECSMALL); R->L11ray= (GEN*)cget1(l, t_VEC);
  prime = utoipos(2);
  for (p = 2; p <= N0; prime[2] = p) {
    switch (krois(dk, p))
    {
    case -1: /* inert */
      if (condZ % p == 0) deg0(R,p); else deg2(R,p);
      break;
    case 1: /* split */
      pr = primedec(nf, prime);
      if      (condZ % p != 0) deg12(R, p, bnr, pr);
      else if (contZ % p == 0) deg0(R,p);
      else
      {
        pr = idealval(nf, cond, gel(pr,1))? gel(pr,2): gel(pr,1);
        deg11(R, p, bnr, pr);
      }
      break;
    default: /* ramified */
      if (condZ % p == 0) deg0(R,p);
      else
      {
        pr = (GEN)primedec(nf, prime)[1];
        deg11(R, p, bnr, pr);
      }
      break;
    }
    NEXT_PRIME_VIADIFF(p,d);
  }
  /* precompute isprincipalray(x), x in Z */
  R->rayZ = (GEN*)cgetg(condZ, t_VEC);
  for (i=1; i<condZ; i++)
    R->rayZ[i] = (cgcd(i,condZ) == 1)? isprincipalray(bnr, utoipos(i)): gen_0;

  gptr[0] = &(R->L0);
  gptr[1] = &(R->L2);  gptr[2] = (GEN*)&(R->rayZ);
  gptr[3] = &(R->L1);  gptr[5] = (GEN*)&(R->L1ray);
  gptr[4] = &(R->L11); gptr[6] = (GEN*)&(R->L11ray);
  gerepilemany(av,gptr,7);
}

static void
InitPrimes(GEN bnr, long N0, LISTray *R)
{
  GEN bnf = gel(bnr,1), cond = gmael3(bnr,2,1,1);
  long np,p,j,k,l, condZ = itos(gcoeff(cond,1,1)), N = lg(cond)-1;
  GEN *tmpray, tabpr, prime, pr, nf = checknf(bnf);
  byteptr d = diffptr + 1;

  R->condZ = condZ; l = PiBound(N0) * N;
  tmpray = (GEN*)cgetg(N+1, t_VEC);
  R->L1 = cget1(l, t_VECSMALL);
  R->L1ray = (GEN*)cget1(l, t_VEC);
  prime = utoipos(2);
  for (p = 2; p <= N0; prime[2] = p)
  {
    pari_sp av = avma;
    if (DEBUGLEVEL>1 && (p & 2047) == 1) fprintferr("%ld ", p);
    tabpr = primedec(nf, prime);
    for (j = 1; j < lg(tabpr); j++)
    {
      pr  = gel(tabpr,j);
      np = itos_or_0( pr_norm(pr) );
      if (!np || np > N0) break;
      if (condZ % p == 0 && idealval(nf, cond, pr))
      {
        tmpray[j] = NULL; continue;
      }

      appendL(R->L1, (GEN)np);
      tmpray[j] = gclone( isprincipalray(bnr, pr) );
    }
    avma = av;
    for (k = 1; k < j; k++)
    {
      if (!tmpray[k]) continue;
      appendL((GEN)R->L1ray, gcopy(tmpray[k]));
      gunclone(tmpray[k]);
    }
    NEXT_PRIME_VIADIFF(p,d);
  }
}

static GEN /* cf polcoeff */
_sercoeff(GEN x, long n)
{
  long i = n - valp(x);
  return (i < 0)? gen_0: gel(x,i+2);
}

static void
affect_coeff(GEN q, long n, GEN y)
{
  GEN x = _sercoeff(q,-n);
  if (x == gen_0) gel(y,n) = gen_0; else gaffect(x, gel(y,n));
}

typedef struct {
  GEN c1, *aij, *bij, *powracpi, *cS, *cT;
  long i0, a,b,c, r, rc1, rc2;
} ST_t;

/* compute the principal part at the integers s = 0, -1, -2, ..., -i0
   of Gamma((s+1)/2)^a Gamma(s/2)^b Gamma(s)^c / (s - z) with z = 0 and 1 */
/* NOTE: surely not the best way to do this, but it's fast enough! */
static void
ppgamma(ST_t *T, long prec)
{
  GEN eul, gam,gamun,gamdm, an,bn,cn_evn,cn_odd, x,x2,X,Y, cf, sqpi;
  GEN p1, p2, *aij, *bij;
  long a = T->a;
  long b = T->b;
  long c = T->c, r = T->r, i0 = T->i0;
  long i,j, s,t;
  pari_sp av;

  aij = (GEN*)cgetg(i0+1, t_VEC);
  bij = (GEN*)cgetg(i0+1, t_VEC);
  for (i = 1; i <= i0; i++)
  {
    aij[i] = p1 = cgetg(r+1, t_VEC);
    bij[i] = p2 = cgetg(r+1, t_VEC);
    for (j=1; j<=r; j++) { gel(p1,j) = cgetr(prec); gel(p2,j) = cgetr(prec); }
  }
  av = avma;

  x   = pol_x[0];
  x2  = gmul2n(x, -1); /* x/2 */
  eul = mpeuler(prec);
  sqpi= sqrtr_abs(mppi(prec)); /* Gamma(1/2) */

  /* expansion of log(Gamma(u)) at u = 1 */
  gamun = cgetg(r+3, t_SER);
  gamun[1] = evalsigne(1) | evalvalp(0) | evalvarn(0);
  gel(gamun,2) = gen_0;
  gel(gamun,3) = gneg(eul);
  for (i = 2; i <= r; i++)
    gel(gamun,i+2) = divrs(szeta(i,prec), odd(i)? -i: i);
  gamun = gexp(gamun, prec); /* Gamma(1 + x) */
  gam = gdiv(gamun,x); /* Gamma(x) */

  /* expansion of log(Gamma(u) / Gamma(1/2)) at u = 1/2 */
  gamdm = cgetg(r+3, t_SER);
  gamdm[1] = evalsigne(1) | evalvalp(0) | evalvarn(0);
  gel(gamdm,2) = gen_0;
  gel(gamdm,3) = gneg(gadd(gmul2n(mplog2(prec), 1), eul));
  for (i = 2; i <= r; i++)
    gel(gamdm,i+2) = mulri(gel(gamun,i+2), subis(int2n(i), 1));
  gamdm = gmul(sqpi, gexp(gamdm, prec)); /* Gamma(1/2 + x) */

 /* We simplify to get one of the following two expressions
  * if (b > a) : sqrt{Pi}^a 2^{a-au} Gamma(u)^{a+c} Gamma(  u/2  )^{|b-a|}
  * if (b <= a): sqrt{Pi}^b 2^{b-bu} Gamma(u)^{b+c) Gamma((u+1)/2)^{|b-a|} */
  if (b > a)
  {
    t = a; s = b; X = x2; Y = gsub(x2,ghalf);
    p1 = gsubst(gam,0,x2);
    p2 = gdiv(gsubst(gamdm,0,x2), Y); /* Gamma((x-1)/2) */
  }
  else
  {
    t = b; s = a; X = gadd(x2,ghalf); Y = x2;
    p1 = gsubst(gamdm,0,x2);
    p2 = gsubst(gam,0,x2);
  }
  cf = gpowgs(sqpi, t);
  an = gpowgs(gpow(gen_2, gsubsg(1,x), prec), t); /* 2^{t-tx} */
  bn = gpowgs(gam, t+c); /* Gamma(x)^{t+c} */
  cn_evn = gpowgs(p1, s-t); /* Gamma(X)^{s-t} */
  cn_odd = gpowgs(p2, s-t); /* Gamma(Y)^{s-t} */
  for (i = 0; i < i0/2; i++)
  {
    GEN C1,q1, A1 = aij[2*i+1], B1 = bij[2*i+1];
    GEN C2,q2, A2 = aij[2*i+2], B2 = bij[2*i+2];

    C1 = gmul(cf, gmul(bn, gmul(an, cn_evn)));
    p1 = gdiv(C1, gsubgs(x, 2*i));
    q1 = gdiv(C1, gsubgs(x, 2*i+1));

    /* an(x-u-1) = 2^t an(x-u) */
    an = gmul2n(an, t);
    /* bn(x-u-1) = bn(x-u) / (x-u-1)^{t+c} */
    bn = gdiv(bn, gpowgs(gsubgs(x, 2*i+1), t+c));

    C2 = gmul(cf, gmul(bn, gmul(an, cn_odd)));
    p2 = gdiv(C2, gsubgs(x, 2*i+1));
    q2 = gdiv(C2, gsubgs(x, 2*i+2));
    for (j = 1; j <= r; j++)
    {
      affect_coeff(p1, j, A1); affect_coeff(q1, j, B1);
      affect_coeff(p2, j, A2); affect_coeff(q2, j, B2);
    }

    an = gmul2n(an, t);
    bn = gdiv(bn, gpowgs(gsubgs(x, 2*i+2), t+c));
    /* cn_evn(x-2i-2) = cn_evn(x-2i)  / (X - (i+1))^{s-t} */
    /* cn_odd(x-2i-3) = cn_odd(x-2i-1)/ (Y - (i+1))^{s-t} */
    cn_evn = gdiv(cn_evn, gpowgs(gsubgs(X,i+1), s-t));
    cn_odd = gdiv(cn_odd, gpowgs(gsubgs(Y,i+1), s-t));
  }
  T->aij = aij;
  T->bij = bij; avma = av;
}

static GEN
_cond(GEN dtcr) { return mkvec2(ch_cond(dtcr), ch_4(dtcr)); }

/* sort chars according to conductor */
static GEN
sortChars(GEN dataCR)
{
  const long cl = lg(dataCR) - 1;
  GEN vCond  = cgetg(cl+1, t_VEC);
  GEN CC     = cgetg(cl+1, t_VECSMALL);
  GEN nvCond = cgetg(cl+1, t_VECSMALL);
  long j,k, ncond;
  GEN vChar;

  for (j = 1; j <= cl; j++) nvCond[j] = 0;

  ncond = 0;
  for (j = 1; j <= cl; j++)
  {
    GEN cond = _cond(gel(dataCR,j));
    for (k = 1; k <= ncond; k++)
      if (gequal(cond, gel(vCond,k))) break;
    if (k > ncond) gel(vCond,++ncond) = cond;
    nvCond[k]++; CC[j] = k; /* char j has conductor number k */
  }
  vChar = cgetg(ncond+1, t_VEC);
  for (k = 1; k <= ncond; k++)
  {
    gel(vChar,k) = cgetg(nvCond[k]+1, t_VECSMALL);
    nvCond[k] = 0;
  }
  for (j = 1; j <= cl; j++)
  {
    k = CC[j]; nvCond[k]++;
    mael(vChar,k,nvCond[k]) = j;
  }
  return vChar;
}

/* Given W(chi), S(chi) and T(chi), return L(1, chi) if fl & 1, else
   [r(chi), c(chi)] where L(s, chi) ~ c(chi) s^r(chi) at s = 0.
   If fl & 2, adjust the value to get L_S(s, chi). */
static GEN
GetValue(GEN dtcr, GEN W, GEN S, GEN T, long fl, long prec)
{
  pari_sp av = avma;
  GEN cf, z, p1;
  long q, b, c, r;
  int isreal = (itos(gel(ch_CHI0(dtcr), 3)) <= 2);

  p1 = ch_4(dtcr);
  q = p1[1];
  b = p1[2];
  c = p1[3];

  if (fl & 1)
  { /* S(chi) + W(chi).T(chi)) / (C(chi) sqrt(Pi)^{r1 - q}) */
    cf = gmul(ch_C(dtcr), powrshalf(mppi(prec), b));

    z = gadd(S, gmul(W, T));
    if (isreal) z = real_i(z);
    z = gdiv(z, cf);
    if (fl & 2) z = gmul(z, ComputeAChi(dtcr, &r, 1, prec));
  }
  else
  { /* (W(chi).S(conj(chi)) + T(chi)) / (sqrt(Pi)^q 2^{r1 - q}) */
    cf = gmul2n(powrshalf(mppi(prec), q), b);

    z = gadd(gmul(W, gconj(S)), gconj(T));
    if (isreal) z = real_i(z);
    z = gdiv(z, cf); r = 0;
    if (fl & 2) z = gmul(z, ComputeAChi(dtcr, &r, 0, prec));
    z = mkvec2(utoi(b + c + r), z);
  }
  return gerepilecopy(av, z);
}

/* return the order and the first non-zero term of L(s, chi0)
   at s = 0. If flag != 0, adjust the value to get L_S(s, chi0). */
static GEN
GetValue1(GEN bnr, long flag, long prec)
{
  GEN bnf = checkbnf(bnr), nf = checknf(bnf);
  GEN h, R, w, c, diff;
  long i, l, r, r1, r2;
  pari_sp av = avma;

  nf_get_sign(nf, &r1,&r2);
  h = gmael3(bnf, 8, 1, 1);
  R = gmael(bnf, 8, 2);
  w = gmael3(bnf, 8, 4, 1);

  c = gneg_i(gdiv(gmul(h, R), w));
  r = r1 + r2 - 1;

  if (flag)
  {
    diff = divcond(bnr);
    l = lg(diff) - 1; r += l;
    for (i = 1; i <= l; i++)
      c = gmul(c, glog(pr_norm(gel(diff,i)), prec));
  }
  return gerepilecopy(av, mkvec2(stoi(r), c));
}

/********************************************************************/
/*                6th part: recover the coefficients                */
/********************************************************************/
static long
TestOne(GEN plg, RC_data *d)
{
  long j, v = d->v;
  GEN z = gsub(d->beta, gel(plg,v));
  if (expo(z) >= d->G) return 0;
  for (j = 1; j < lg(plg); j++)
    if (j != v && mpcmp(d->B, mpabs(gel(plg,j))) < 0) return 0;
  return 1;
}

static GEN
chk_reccoeff_init(FP_chk_fun *chk, GEN r, GEN mat)
{
  RC_data *d = (RC_data*)chk->data;
  (void)r; d->U = mat; return d->nB;
}

static GEN
chk_reccoeff(void *data, GEN x)
{
  RC_data *d = (RC_data*)data;
  GEN v = gmul(d->U, x), z = gel(v,1);

  if (!gcmp1(z)) return NULL;
  *++v = evaltyp(t_COL) | evallg( lg(d->M) );
  if (TestOne(gmul(d->M, v), d)) return v;
  return NULL;
}

/* Using Cohen's method */
static GEN
RecCoeff3(GEN nf, RC_data *d, long prec)
{
  GEN A, M, nB, cand, p1, B2, C2, tB, beta2, BIG, nf2, Bd;
  GEN beta = d->beta, B = d->B;
  long N = d->N, v = d->v, e;
  long i, j, k, l, ct = 0, prec2;
  FP_chk_fun chk = { &chk_reccoeff, &chk_reccoeff_init, NULL, NULL, 0 };
  chk.data = (void*)d;

  d->G = min(-10, -bit_accuracy(prec) >> 4);
  /* FIXME: why a power of 10 ? Use a power of 2. */
  BIG = powuu(10, max(8, -(d->G >> 1)));
  tB  = sqrtnr(gmul2n(itor(BIG,DEFAULTPREC), -N), N-1);
  Bd    = grndtoi(gmin(B, tB), &e);
  if (e > 0) return NULL; /* failure */
  Bd = addis(Bd, 1);
  prec2 = BIGDEFAULTPREC + (expi(Bd)>>TWOPOTBITS_IN_LONG);
  prec2 = max((prec << 1) - 2, prec2);
  nf2   = nfnewprec(nf, prec2);
  beta2 = gprec_w(beta, prec2);

LABrcf: ct++;
  B2 = sqri(Bd);
  C2 = mulii(B2, sqri(BIG));

  M = gmael(nf2, 5, 1);
  d->M = M;

  A = cgetg(N+2, t_MAT);
  for (i = 1; i <= N+1; i++)
    gel(A,i) = cgetg(N+2, t_COL);

  gcoeff(A, 1, 1) = gadd(gmul(C2, gsqr(beta2)), B2);
  for (j = 2; j <= N+1; j++)
  {
    p1 = gmul(C2, gmul(gneg_i(beta2), gcoeff(M, v, j-1)));
    gcoeff(A, 1, j) = gcoeff(A, j, 1) = p1;
  }
  for (i = 2; i <= N+1; i++)
    for (j = 2; j <= N+1; j++)
    {
      p1 = gen_0;
      for (k = 1; k <= N; k++)
      {
        GEN p2 = gmul(gcoeff(M, k, j-1), gcoeff(M, k, i-1));
        if (k == v) p2 = gmul(C2, p2);
        p1 = gadd(p1,p2);
      }
      gcoeff(A, i, j) = gcoeff(A, j, i) = p1;
    }

  nB = mulsi(N+1, B2);
  d->nB = nB;
  cand = fincke_pohst(A, nB, -1, prec2, &chk);

  if (!cand)
  {
    if (ct > 3) return NULL;

    prec2 = (prec2 << 1) - 2;
    if (DEBUGLEVEL>1) pari_warn(warnprec,"RecCoeff", prec2);
    nf2 = nfnewprec(nf2, prec2);
    beta2 = gprec_w(beta2, prec2);
    goto LABrcf;
  }

  cand = gel(cand,1);
  l = lg(cand) - 1;

  if (l == 1) return coltoalg(nf, gel(cand,1));

  if (DEBUGLEVEL>1) fprintferr("RecCoeff3: no solution found!\n");
  return NULL;
}

/* Using linear dependance relations */
static GEN
RecCoeff2(GEN nf,  RC_data *d,  long prec)
{
  pari_sp av;
  GEN vec, M = gmael(nf, 5, 1), beta = d->beta;
  long i, imin, imax, lM = lg(M);

  d->G = min(-20, -bit_accuracy(prec) >> 4);

  vec  = shallowconcat(mkvec(gneg(beta)), row(M, d->v));
  imin = (long)bit_accuracy_mul(prec, .225);
  imax = (long)bit_accuracy_mul(prec, .315);

  av = avma;
  for (i = imax; i >= imin; i-=16, avma = av)
  {
    long e;
    GEN v = lindep2(vec, i), z = gel(v,1);
    if (!signe(z)) continue;
    *++v = evaltyp(t_COL) | evallg(lM);
    v = grndtoi(gdiv(v, z), &e);
    if (e > 0) break;
    if (TestOne(gmul(M, v), d)) return coltoalg(nf, v);
  }
  /* failure */
  return RecCoeff3(nf,d,prec);
}

/* Attempts to find a polynomial with coefficients in nf such that
   its coefficients are close to those of pol at the place v and
   less than B at all the other places */
static GEN
RecCoeff(GEN nf,  GEN pol,  long v, long prec)
{
  long j, md, cl = degpol(pol);
  pari_sp av = avma;
  RC_data d;

  /* if precision(pol) is too low, abort */
  for (j = 2; j <= cl+1; j++)
  {
    GEN t = gel(pol, j);
    if (bit_accuracy(gprecision(t)) - gexpo(t) < 34) return NULL;
  }

  md = cl/2;
  pol = shallowcopy(pol);

  d.N = degpol(nf[1]);
  d.v = v;

  for (j = 1; j <= cl; j++)
  { /* start with the coefficients in the middle,
       since they are the harder to recognize! */
    long cf = md + (j%2? j/2: -j/2);
    GEN t, bound = shifti(binomial(utoipos(cl), cf), cl-cf);

    if (DEBUGLEVEL>1) fprintferr("RecCoeff (cf = %ld, B = %Z)\n", cf, bound);
    d.beta = real_i( gel(pol,cf+2) );
    d.B    = bound;
    if (! (t = RecCoeff2(nf, &d, prec)) ) return NULL;
    gel(pol, cf+2) = t;
  }
  gel(pol,cl+2) = gen_1;
  return gerepilecopy(av, pol);
}

/* an[q * i] *= chi for all (i,p)=1 */
static void
an_mul(int **an, long p, long q, long n, long deg, GEN chi, int **reduc)
{
  pari_sp av;
  long c,i;
  int *T;

  if (gcmp1(chi)) return;
  av = avma;
  T = (int*)new_chunk(deg); Polmod2Coeff(T,chi, deg);
  for (c = 1, i = q; i <= n; i += q, c++)
    if (c == p) c = 0; else MulCoeff(an[i], T, reduc, deg);
  avma = av;
}
/* an[q * i] = 0 for all (i,p)=1 */
static void
an_set0_coprime(int **an, long p, long q, long n, long deg)
{
  long c,i;
  for (c = 1, i = q; i <= n; i += q, c++)
    if (c == p) c = 0; else _0toCoeff(an[i], deg);
}
/* an[q * i] = 0 for all i */
static void
an_set0(int **an, long p, long n, long deg)
{
  long i;
  for (i = p; i <= n; i += p) _0toCoeff(an[i], deg);
}

/* compute the coefficients an for the quadratic case */
static int**
computean(GEN dtcr, LISTray *R, long n, long deg)
{
  pari_sp av = avma, av2;
  long i, p, q, condZ, l;
  int **an, **reduc;
  GEN L, CHI, chi, chi1;
  CHI_t C;

  CHI = ch_CHI(dtcr); init_CHI_alg(&C, CHI);
  condZ= R->condZ;

  an = InitMatAn(n, deg, 1);
  reduc = InitReduction(CHI, deg);
  av2 = avma;

  /* all pr | p divide cond */
  L = R->L0; l = lg(L);
  for (i=1; i<l; i++) an_set0(an,L[i],n,deg);

  /* 1 prime of degree 2 */
  L = R->L2; l = lg(L);
  for (i=1; i<l; i++, avma = av2)
  {
    p = L[i];
    if (condZ == 1) chi = C.val[0]; /* 1 */
    else            chi = EvalChar(&C, R->rayZ[p % condZ]);
    chi1 = chi;
    for (q=p;;)
    {
      an_set0_coprime(an, p,q,n,deg); /* v_p(q) odd */
      if (! (q = next_pow(q,p, n)) ) break;

      an_mul(an,p,q,n,deg,chi,reduc);
      if (! (q = next_pow(q,p, n)) ) break;
      chi = gmul(chi, chi1);
    }
  }

  /* 1 prime of degree 1 */
  L = R->L1; l = lg(L);
  for (i=1; i<l; i++, avma = av2)
  {
    p = L[i];
    chi = EvalChar(&C, R->L1ray[i]);
    chi1 = chi;
    for(q=p;;)
    {
      an_mul(an,p,q,n,deg,chi,reduc);
      if (! (q = next_pow(q,p, n)) ) break;
      chi = gmul(chi, chi1);
    }
  }

  /* 2 primes of degree 1 */
  L = R->L11; l = lg(L);
  for (i=1; i<l; i++, avma = av2)
  {
    GEN ray1, ray2, chi11, chi12, chi2;

    p = L[i]; ray1 = R->L11ray[i]; /* use pr1 pr2 = (p) */
    if (condZ == 1)
      ray2 = gneg(ray1);
    else
      ray2 = gsub(R->rayZ[p % condZ],  ray1);
    chi11 = EvalChar(&C, ray1);
    chi12 = EvalChar(&C, ray2);

    chi1 = gadd(chi11, chi12);
    chi2 = chi12;
    for(q=p;;)
    {
      an_mul(an,p,q,n,deg,chi1,reduc);
      if (! (q = next_pow(q,p, n)) ) break;
      chi2 = gmul(chi2, chi12);
      chi1 = gadd(chi2, gmul(chi1, chi11));
    }
  }

  CorrectCoeff(dtcr, an, reduc, n, deg);
  FreeMat(reduc, deg-1);
  avma = av; return an;
}

/* compute S and T for the quadratic case where
   the extension is of the type used to construct
   abelian extensions using Stark units */
static void
QuadGetST(GEN bnr, GEN *pS, GEN *pT, GEN dataCR, GEN vChar, long prec)
{
  const long cl  = lg(dataCR) - 1;
  pari_sp av, av1, av2;
  long ncond, n, j, k, n0;
  GEN N0, C, T, S, cf, cfh, an, degs;
  LISTray LIST;

  /* allocate memory for answer */
  *pS = S = cgetg(cl+1, t_VEC);
  *pT = T = cgetg(cl+1, t_VEC);
  for (j = 1; j <= cl; j++)
  {
    gel(S,j) = cgetc(prec);
    gel(T,j) = cgetc(prec);
  }
  av = avma;

  /* initializations */
  degs = GetDeg(dataCR);
  ncond = lg(vChar)-1;
  C    = cgetg(ncond+1, t_VEC);
  N0   = cgetg(ncond+1, t_VECSMALL);
  n0 = 0;
  for (j = 1; j <= ncond; j++)
  {
    GEN dtcr = gel(dataCR, mael(vChar,j,1)), c = ch_C(dtcr);
    gel(C,j) = c;
    N0[j] = (long)bit_accuracy_mul(prec, 0.35 * gtodouble(c));
    if (n0 < N0[j]) n0 = N0[j];
  }
  if ((ulong)n0 > maxprime())
    pari_err(talker, "Not enough precomputed primes (need all p <= %ld)", n0);
  if (DEBUGLEVEL>1) fprintferr("N0 = %ld\n", n0);
  InitPrimesQuad(bnr, n0, &LIST);

  cfh = sqrtr(mppi(prec));
  cf  = gmul2n(cfh, 1);
  av1 = avma;
  /* loop over conductors */
  for (j = 1; j <= ncond; j++)
  {
    const GEN c1 = gel(C,j), c2 = divsr(2,c1), cexp = mpexp(gneg(c2));
    const GEN LChar = gel(vChar,j);
    const long nChar = lg(LChar)-1, NN = N0[j];
    GEN veint1, vcn = cgetg(NN+1, t_VEC);

    if (DEBUGLEVEL>1)
      fprintferr("* conductor no %ld/%ld (N = %ld)\n\tInit: ", j,ncond,NN);
    veint1 = veceint1(c2, stoi(NN), prec);
    gel(vcn,1) = cexp;
    for (n=2; n<=NN; n++) gel(vcn,n) = mulrr(gel(vcn,n-1), cexp);
    av2 = avma;
    for (n=2; n<=NN; n++, avma = av2)
      affrr(divrs(gel(vcn,n),n), gel(vcn,n));

    for (k = 1; k <= nChar; k++)
    {
      const long t = LChar[k], d = degs[t];
      const GEN dtcr = gel(dataCR, t), z = gel(ch_CHI(dtcr), 2);
      GEN p1 = gen_0, p2 = gen_0;
      int **matan;
      long c = 0;

      if (DEBUGLEVEL>1)
        fprintferr("\tcharacter no: %ld (%ld/%ld)\n", t,k,nChar);
      matan = computean(gel(dataCR,t), &LIST, NN, d);
      for (n = 1; n <= NN; n++)
	if ((an = EvalCoeff(z, matan[n], d)))
        {
          p1 = gadd(p1, gmul(an, gel(vcn,n)));
	  p2 = gadd(p2, gmul(an, gel(veint1,n)));
          if (++c == 256) { gerepileall(av2,2, &p1,&p2); c = 0; }
        }
      gaffect(gmul(cfh, gmul(p1,c1)), gel(S,t));
      gaffect(gmul(cf,  gconj(p2)),   gel(T,t));
      FreeMat(matan,NN); avma = av2;
    }
    if (DEBUGLEVEL>1) fprintferr("\n");
    avma = av1;
  }
  if (DEBUGLEVEL) msgtimer("S & T");
  avma = av;
}

/* S & T for the general case */
static void
get_cS_cT(ST_t *T, long n)
{
  pari_sp av;
  GEN csurn, nsurc, lncsurn;
  GEN A,B,s,t,Z,*aij,*bij;
  long i,j,r,i0;

  if (T->cS[n]) return;

  av = avma;
  aij = T->aij; i0= T->i0;
  bij = T->bij; r = T->r;
  Z = cgetg(r+1, t_VEC);
  Z[1] = 0; /* unused */

  csurn = divrs(T->c1, n);
  nsurc = ginv(csurn);
  lncsurn = logr_abs(csurn);

  gel(Z,2) = lncsurn; /* r >= 2 */
  for (i = 3; i <= r; i++)
    gel(Z,i) = divrs(mulrr(gel(Z,i-1), lncsurn), i-1);
  /* Z[i] = ln^(i-1)(c1/n) / (i-1)! */

  /* i = i0 */
    A = aij[i0]; t = gel(A,1);
    B = bij[i0]; s = gel(B,1);
    for (j = 2; j <= r; j++)
    {
      if (signe(B[j])) s = mpadd(s, mulrr(gel(Z,j), gel(B,j)));
      if (signe(A[j])) t = mpadd(t, mulrr(gel(Z,j), gel(A,j)));
    }
  for (i = i0 - 1; i > 1; i--)
  {
    A = aij[i]; if (signe(t)) t = mulrr(t, nsurc);
    B = bij[i]; if (signe(s)) s = mulrr(s, nsurc);
    for (j = odd(i)? T->rc2: T->rc1; j > 1; j--)
    {
      if (signe(B[j])) s = addrr(s, mulrr(gel(Z,j), gel(B,j)));
      if (signe(A[j])) t = addrr(t, mulrr(gel(Z,j), gel(A,j)));
    }
    if (signe(B[1])) s = addrr(s, gel(B,1));
    if (signe(A[1])) t = addrr(t, gel(A,1));
  }
  /* i = 1 */
    A = aij[1]; if (signe(t)) t = mulrr(t, nsurc);
    B = bij[1]; if (signe(s)) s = mulrr(s, nsurc);
    if (signe(B[1])) s = addrr(s, gel(B,1));
    if (signe(A[1])) t = addrr(t, gel(A,1));
    for (j = 2; j <= r; j++)
    {
      if (signe(B[j])) s = addrr(s, mulrr(gel(Z,j), gel(B,j)));
      if (signe(A[j])) t = addrr(t, mulrr(gel(Z,j), gel(A,j)));
    }
  s = mpadd(s, mpmul(csurn, T->powracpi[T->b]));
  T->cS[n] = gclone(s);
  T->cT[n] = gclone(t); avma = av;
}

static void
clear_cScT(ST_t *T, long N)
{
  GEN *cS = T->cS, *cT = T->cT;
  long i;
  for (i=1; i<=N; i++)
    if (cS[i]) { gunclone(cS[i]); gunclone(cT[i]); cS[i] = cT[i] = NULL; }
}

static void
init_cScT(ST_t *T, GEN dtcr, long N, long prec)
{
  GEN p1 = ch_4(dtcr);
  T->a = p1[1];
  T->b = p1[2];
  T->c = p1[3];
  T->rc1 = T->a + T->c;
  T->rc2 = T->b + T->c;
  T->r   = max(T->rc2+1, T->rc1); /* >= 2 */
  ppgamma(T, prec);
  clear_cScT(T, N);
}

static void
GetST(GEN bnr, GEN *pS, GEN *pT, GEN dataCR, GEN vChar, long prec)
{
  const long cl = lg(dataCR) - 1;
  pari_sp av, av1, av2;
  long ncond, n, j, k, jc, n0, prec2, i0, r1, r2;
  GEN nf, racpi, *powracpi;
  GEN N0, C, T, S, an, degs, limx;
  LISTray LIST;
  ST_t cScT;

  nf  = checknf(bnr);

  if (DEBUGLEVEL) (void)timer2();
  /* allocate memory for answer */
  *pS = S = cgetg(cl+1, t_VEC);
  *pT = T = cgetg(cl+1, t_VEC);
  for (j = 1; j <= cl; j++)
  {
    gel(S,j) = cgetc(prec);
    gel(T,j) = cgetc(prec);
  }
  av = avma;

  /* initializations */
  degs = GetDeg(dataCR);
  ncond = lg(vChar)-1;
  nf_get_sign(nf,&r1,&r2);

  C  = cgetg(ncond+1, t_VEC);
  N0 = cgetg(ncond+1, t_VECSMALL);
  n0 = 0;
  limx = zeta_get_limx(r1, r2, bit_accuracy(prec));
  for (j = 1; j <= ncond; j++)
  {
    GEN dtcr = gel(dataCR, mael(vChar,j,1)), c = ch_C(dtcr);
    gel(C,j) = c;
    N0[j] = zeta_get_N0(c, limx);
    if (n0 < N0[j]) n0  = N0[j];
  }
  if ((ulong)n0 > maxprime())
    pari_err(talker, "Not enough precomputed primes (need all p <= %ld)", n0);
  i0 = zeta_get_i0(r1, r2, bit_accuracy(prec), limx);
  InitPrimes(bnr, n0, &LIST);

  prec2 = ((prec-2) << 1) + EXTRA_PREC;
  racpi = sqrtr(mppi(prec2));
  powracpi = (GEN*)cgetg(r1+2,t_VEC);
  powracpi++; powracpi[0] = gen_1; powracpi[1] = racpi;
  for (j=2; j<=r1; j++) powracpi[j] = mulrr(powracpi[j-1], racpi);
  cScT.powracpi = powracpi;

  cScT.cS = (GEN*)cgetg(n0+1, t_VEC);
  cScT.cT = (GEN*)cgetg(n0+1, t_VEC);
  for (j=1; j<=n0; j++) cScT.cS[j] = cScT.cT[j] = NULL;

  cScT.i0 = i0;

  av1 = avma;
  for (jc = 1; jc <= ncond; jc++)
  {
    const GEN LChar = gel(vChar,jc);
    const long nChar = lg(LChar)-1, NN = N0[jc];

    if (DEBUGLEVEL>1)
      fprintferr("* conductor no %ld/%ld (N = %ld)\n\tInit: ", jc,ncond,NN);

    cScT.c1 = gel(C,jc);
    init_cScT(&cScT, gel(dataCR, LChar[1]), NN, prec2);
    av2 = avma;
    for (k = 1; k <= nChar; k++)
    {
      const long t = LChar[k], d = degs[t];
      const GEN dtcr = gel(dataCR, t), z = gel(ch_CHI(dtcr), 2);
      GEN p1 = gen_0, p2 = gen_0;
      long c = 0;
      int **matan;

      if (DEBUGLEVEL>1)
        fprintferr("\tcharacter no: %ld (%ld/%ld)\n", t,k,nChar);
      matan = ComputeCoeff(gel(dataCR,t), &LIST, NN, d);
      for (n = 1; n <= NN; n++)
        if ((an = EvalCoeff(z, matan[n], d)))
        {
          get_cS_cT(&cScT, n);
          p1 = gadd(p1, gmul(an, cScT.cS[n]));
          p2 = gadd(p2, gmul(an, cScT.cT[n]));
          if (++c == 256) { gerepileall(av2,2, &p1,&p2); c = 0; }
        }
      gaffect(p1,        gel(S,t));
      gaffect(gconj(p2), gel(T,t));
      FreeMat(matan, NN); avma = av2;
    }
    if (DEBUGLEVEL>1) fprintferr("\n");
    avma = av1;
  }
  if (DEBUGLEVEL) msgtimer("S&T");
  clear_cScT(&cScT, n0);
  avma = av;
}

/*******************************************************************/
/*                                                                 */
/*     Class fields of real quadratic fields using Stark units     */
/*                                                                 */
/*******************************************************************/

#if 0
typedef struct {
  long cl;
  GEN dkpow;
} DH_t;

/* return 1 if the absolute polynomial pol (over Q) defines the
   Hilbert class field of the real quadratic field bnf */
static GEN
define_hilbert(void *S, GEN pol)
{
  DH_t *T = (DH_t*)S;
  GEN d = modulargcd(derivpol(pol), pol);

  if (degpol(pol) != T->cl + degpol(d)) return NULL;
  pol = gdivexact(pol, d);
  return (T->cl & 1 || !equalii(smalldiscf(pol), T->dkpow))? pol: NULL;
}

/* let polrel define Hk/k,  find L/Q such that Hk=Lk and L and k are
   disjoint */
static GEN
makescindold(GEN nf, GEN polrel, long cl)
{
  long i, l;
  pari_sp av = avma;
  GEN pol, L, BAS, nf2, dk;
  DH_t T;
  FP_chk_fun CHECK;

  BAS = rnfpolredabs(nf, polrel, nf_ABSOLUTE|nf_ADDZK|nf_PARTIALFACT);

  /* attempt to find the subfields using polred */
  T.cl = cl; dk  = gel(nf,3);
  T.dkpow = (cl & 1) ? NULL: powiu(dk, cl>>1);
  CHECK.f = &define_hilbert;
  CHECK.data = (void*)&T;
  pol = polredfirstpol(BAS, 0, &CHECK);
  if (DEBUGLEVEL) msgtimer("polred");

  if (!pol)
  {
    nf2 = nfinit0(BAS, 0, DEFAULTPREC);
    L  = subfields(nf2, utoipos(cl));
    l = lg(L);
    if (DEBUGLEVEL) msgtimer("subfields");

    for (i = 1; i < l; i++)
    {
      pol = gmael(L, i, 1);
      if (cl & 1 || !equalii(smalldiscf(pol), T.dkpow)) break;
    }
    if (i == l)
      for (i = 1; i < l; i++)
      {
        pol = gmael(L, i, 1);
        if (degpol(gcoeff(nffactor(nf, pol), 1, 1)) == cl) break;
      }
    if (i == l)
      pari_err(bugparier, "makescindold (no polynomial found)");
  }
  pol = polredabs0(pol, nf_PARTIALFACT);
  return gerepileupto(av, pol);
}
#endif

/* compute the Hilbert class field using genus class field theory when
   the exponent of the class group is 2 */
static GEN
GenusField(GEN bnf)
{
  long hk, c, l;
  pari_sp av = avma;
  GEN disc, x2, pol, div, d;

  hk   = itos(gmael3(bnf, 8, 1, 1));
  disc = gmael(bnf, 7, 3);
  x2   = gsqr(pol_x[0]);

  if (mod4(disc) == 0) disc = divis(disc, 4);
  div = divisors(disc);
  l = 0;
  pol = NULL;

  for (c = 2; l < hk; c++) /* skip c = 1 : div[1] = gen_1 */
  {
    d = gel(div,c);
    if (mod4(d) == 1)
    {
      GEN t = gsub(x2, d);
      if (!pol)
	pol = t;
      else
	pol = (GEN)compositum(pol, t)[1];

      l = degpol(pol);
    }
  }

  return gerepileupto(av, polredabs0(pol, nf_PARTIALFACT));
}

/* if flag != 0, computes a fast and crude approximation of the result */
static GEN
AllStark(GEN data,  GEN nf,  long flag,  long newprec)
{
  const long BND = 300;
  long cl, i, j, cpt = 0, N, h, v, n, r1, r2, den;
  pari_sp av, av2;
  int **matan;
  GEN bnr, p1, p2, S, T, polrelnum, polrel, Lp, W, veczeta, sig;
  GEN vChar, degs, C, dataCR, cond1, L1, an;
  LISTray LIST;

  bnr = gel(data,1);
  nf_get_sign(nf, &r1,&r2);
  N     = degpol(nf[1]);
  cond1 = gmael3(bnr, 2, 1, 2);
  dataCR = gel(data,5);
  vChar = sortChars(dataCR);

  v = 1;
  while (gcmp1(gel(cond1,v))) v++;

  cl = lg(dataCR)-1;
  degs = GetDeg(dataCR);
  h  = itos(det(gel(data,2))) >> 1;

LABDOUB:
  av = avma;
  W = ComputeAllArtinNumbers(dataCR, vChar, (flag >= 0), newprec);
  if (DEBUGLEVEL) msgtimer("Compute W");
  Lp = cgetg(cl + 1, t_VEC);
  if (!flag)
  {
    if (degpol(nf[1]) == 2) 
      QuadGetST(bnr, &S,&T,dataCR,vChar,newprec); 
    else
      GetST(bnr, &S, &T, dataCR, vChar, newprec);
    for (i = 1; i <= cl; i++)
      Lp[i] = GetValue(gel(dataCR,i), gel(W,i), gel(S,i), gel(T,i),
                       2, newprec)[2];
  }
  else
  { /* compute a crude approximation of the result */
    C = cgetg(cl + 1, t_VEC);
    for (i = 1; i <= cl; i++) gel(C,i) = ch_C(gel(dataCR, i));
    n = zeta_get_N0(vecmax(C), zeta_get_limx(r1, r2, bit_accuracy(newprec)));
    if (n > BND) n = BND;
    if (DEBUGLEVEL) fprintferr("N0 in QuickPol: %ld \n", n);
    InitPrimes(bnr, n, &LIST);

    L1 = cgetg(cl+1, t_VEC);
    /* use L(1) = sum (an / n) */
    for (i = 1; i <= cl; i++)
    {
      GEN dtcr = gel(dataCR,i);
      matan = ComputeCoeff(dtcr, &LIST, n, degs[i]);
      av2 = avma;
      p1 = real_0(newprec); p2 = gel(ch_CHI(dtcr), 2);
      for (j = 1; j <= n; j++)
	if ( (an = EvalCoeff(p2, matan[j], degs[i])) )
          p1 = gadd(p1, gdivgs(an, j));
      gel(L1,i) = gerepileupto(av2, p1);
      FreeMat(matan, n);
    }
    p1 = gmul2n(powrshalf(mppi(newprec), N-2), 1);

    for (i = 1; i <= cl; i++)
    {
      long r;
      GEN WW, A = ComputeAChi(gel(dataCR,i), &r, 0, newprec);
      WW = gmul(gel(C,i), gmul(A, gel(W,i)));
      gel(Lp,i) = gdiv(gmul(WW, gconj(gel(L1,i))), p1);
    }
  }

  p1 = ComputeLift(gel(data,4));

  den = flag ? h: 2*h;
  veczeta = cgetg(h + 1, t_VEC);
  for (i = 1; i <= h; i++)
  {
    GEN z = gen_0;

    sig = gel(p1,i);
    for (j = 1; j <= cl; j++)
    {
      GEN dtcr = gel(dataCR,j), CHI = ch_CHI(dtcr);
      GEN val = ComputeImagebyChar(CHI, sig);
      GEN p2 = real_i(gmul(gel(Lp,j), val));
      if (itos(gel(CHI,3)) != 2) p2 = gmul2n(p2, 1); /* character not real */
      z = gadd(z, p2);
    }
    gel(veczeta,i) = gdivgs(z, den);
  }
  if (DEBUGLEVEL>1) fprintferr("zetavalues = %Z\n", veczeta);

  if (DEBUGLEVEL>1 && !flag)
    fprintferr("Checking the square-root of the Stark unit...\n");

  for (j = 1; j <= h; j++)
    gel(veczeta,j) = gmul2n(gch(gel(veczeta,j), newprec), 1);
  polrelnum = roots_to_pol_intern(real_1(newprec),veczeta, 0,0);
  if (DEBUGLEVEL)
  {
    if (DEBUGLEVEL>1) fprintferr("polrelnum = %Z\n", polrelnum);
    msgtimer("Compute %s", (flag)? "quickpol": "polrelnum");
  }

  if (flag)
    return gerepilecopy(av, polrelnum);

  /* try to recognize this polynomial */
  polrel = RecCoeff(nf, polrelnum, v, newprec);

  if (!polrel)
  {
    if (DEBUGLEVEL>1)
    fprintferr("It's not a square...\n");
    for (j = 1; j <= h; j++)
      gel(veczeta,j) = gsubgs(gsqr(gel(veczeta,j)), 2);
    polrelnum = roots_to_pol_intern(real_1(newprec),veczeta, 0,0);
    if (DEBUGLEVEL)
    {
      if (DEBUGLEVEL>1) fprintferr("polrelnum = %Z\n", polrelnum);
      msgtimer("Compute polrelnum");
    }
    /* try to recognize this polynomial */
    polrel = RecCoeff(nf, polrelnum, v, newprec);
  }

  if (!polrel) /* if it fails... */
  {
    long incr_pr;
    if (++cpt >= 3) pari_err(precer, "stark (computation impossible)");

    /* compute the precision, we need 
          a) get at least EXTRA_PREC fractional digits if there is none;
       or b) double the fractional digits.
    */
    incr_pr = (bit_accuracy(gprecision(polrelnum)) - 
	       gexpo(polrelnum))>>TWOPOTBITS_IN_LONG;
    if (incr_pr < 0) 
      incr_pr = -incr_pr + EXTRA_PREC;
    newprec = newprec + max(ADD_PREC, cpt*incr_pr);

    if (DEBUGLEVEL) pari_warn(warnprec, "AllStark", newprec);

    nf = nfnewprec(nf, newprec);
    dataCR = CharNewPrec(dataCR, nf, newprec);

    gerepileall(av, 2, &nf, &dataCR);
    goto LABDOUB;
  }

  if (DEBUGLEVEL>1) fprintferr("polrel = %Z\n", polrel);
  if (DEBUGLEVEL) msgtimer("Recpolnum");

  return gerepilecopy(av, polrel);
}

/********************************************************************/
/*                        Main functions                            */
/********************************************************************/
/* conj(Mod(x,y)), assume y normalized */
static GEN
quad_conj(GEN x, GEN y)
{
  GEN z, u, v, b;
  if (typ(x) != t_POL || degpol(x) <= 0) return x;
  u = gel(x,3); /*Mod(ux + v, x^2 + bx + c)*/
  v = gel(x,2); b = gel(y,3);
  z = cgetg(4, t_POL); z[1] = x[1];
  gel(z,2) = gadd(v, gmul(u,negi(b)));
  gel(z,3) = gneg(u); return z;
}
static GEN
pol_quad_conj(GEN x, GEN y)
{
  long i, l = lg(x);
  GEN z = cgetg(l, t_POL); z[1] = x[1];
  for (i = 2; i < l; i++) gel(z,i) = quad_conj(gel(x,i), y);
  return z;
}
/* k = nf quadratic field, P relative equation of H_k (Hilbert class field)
 * return T in Z[X], such that H_k / Q is the compositum of Q[X]/(T) and k */
static GEN
makescind(GEN nf, GEN P, long cl)
{
  GEN Pp, p, perm, pol, G, L, a, roo, nfpol = gel(nf,1);
  long i, k, l, is_P;

  P = lift_intern(P);
  pol = RgX_mul(P, pol_quad_conj(P, nfpol)); /* Norm_{k/Q}(P), irreducible/Q */
  for (i = 2; i < lg(pol); i++)
  {
    GEN c = gel(pol,i);
    if (typ(c) != t_POL) continue;
    c = RgX_rem(c, nfpol); /* ZX, degree <= 0 */
    c = signe(c)? gel(c,2) : gen_0;
    gel(pol,i) = c;
  }
  /* pol = rnfequation(nf, P); */
  G = galoisinit(pol, NULL);
  L = gel(G,6);
  p = gmael(G,2,1);
  a = FpX_quad_root(nfpol, p, 0);
  Pp = gsubst(P, varn(nfpol), a);
  Pp = FpX_red(Pp, p); /* P mod a prime \wp above p (which splits) */
  roo = gel(G,3);
  is_P = gcmp0( FpX_eval(Pp, remii(gel(roo,1),p), p) );
  /* each roo[i] mod p is a root of P or (exclusive) tau(P) mod \wp */
  /* record whether roo[1] is a root of P or tau(P) */
  
  perm = NULL; /*-Wall*/
  for (i = 1; lg(L); i++)
  {
    perm = gel(L,i);
    k = perm[1]; if (k == 1) continue;
    k = gcmp0( FpX_eval(Pp, remii(gel(roo,k),p), p) );
    /* roo[k] is a root of the other polynomial */
    if (k != is_P) break;
  }

  l = perm_order(perm);
  if (l != 2) perm = gpowgs(perm, l >> 1);
  /* perm has order two and doesn't belong to Gal(H_k/k) */
  return galoisfixedfield(G, perm, 1, varn(P));
}

/* compute the polynomial over Q of the Hilbert class field of
   Q(sqrt(D)) where D is a positive fundamental discriminant */
GEN
quadhilbertreal(GEN D, long prec)
{
  pari_sp av = avma;
  long newprec;
  VOLATILE long cl;
  VOLATILE GEN pol, bnf, bnr, dtQ, data, nf, exp, M;

  (void)&prec; /* prevent longjmp clobbering it */
  if (DEBUGLEVEL) (void)timer2();
  disable_dbg(0);

  /* quick computation of the class number */
  cl = itos((GEN)quadclassunit0(D, 0, NULL, prec)[1]);
  if (cl == 1) { disable_dbg(-1); avma = av; return pol_x[0]; }

START:
  pol = quadpoly0(D, fetch_user_var("y"));
  bnf = bnfinit0(pol, 1, NULL, prec);
  nf  = gel(bnf,7);
  disable_dbg(-1);
  if (DEBUGLEVEL) msgtimer("Compute Cl(k)");

  /* if the exponent of the class group is 2, use Genus Theory */
  exp = gmael4(bnf, 8, 1, 2, 1);
  if (equaliu(exp,2)) return gerepileupto(av, GenusField(bnf));

  CATCH(precer) {
    prec += EXTRA_PREC; pol = NULL;
    pari_err (warnprec, "quadhilbertreal", prec);
  } TRY {
    /* find the modulus defining N */
    bnr  = buchrayinitgen(bnf, gen_1);
    M = diagonal_i(gmael(bnr,5,2));
    dtQ  = InitQuotient(M);
    data = FindModulus(bnr, dtQ, &newprec, prec);
    if (DEBUGLEVEL) msgtimer("FindModulus");

    if (!data)
    {
      long i, l = lg(M);
      GEN vec = cgetg(l, t_VEC);
      for (i = 1; i < l; i++)
	{
          GEN t = gcoeff(M,i,i);
	  gcoeff(M,i,i) = gen_1;
	  gel(vec,i) = bnrstark(bnr, M, prec);
	  gcoeff(M,i,i) = t;
	}
      CATCH_RELEASE();
      return vec;
    }

    if (newprec > prec)
    {
      if (DEBUGLEVEL>1) fprintferr("new precision: %ld\n", newprec);
      nf = nfnewprec(nf, newprec);
    }
    pol = AllStark(data, nf, 0, newprec);
  } ENDCATCH;
  if (!pol) goto START;

  return gerepileupto(av, makescind(nf, pol, cl));
}

static GEN
get_subgroup(GEN subgp, GEN cyc)
{
  if (!subgp || gcmp0(subgp)) return cyc;
  subgp = hnf(subgp);
  return hnfdivide(subgp, cyc)? subgp: NULL;
}

GEN
bnrstark(GEN bnr, GEN subgrp, long prec)
{
  long N, newprec;
  pari_sp av = avma;
  GEN bnf, p1, Mcyc, nf, data, dtQ;

  /* check the bnr */
  checkbnrgen(bnr);
  bnf = checkbnf(bnr);
  nf  = checknf(bnf);
  N   = degpol(nf[1]);
  if (N == 1) return galoissubcyclo(bnr, subgrp, 0, 0);

  /* check the bnf */
  if (!varn(nf[1])) pari_err(talker, "main variable in bnrstark must not be x");
  if (nf_get_r2(nf)) pari_err(talker, "base field not totally real in bnrstark");

  /* check the subgrp */
  Mcyc = diagonal_i(gmael(bnr, 5, 2));
  if (! (subgrp = get_subgroup(subgrp,Mcyc)) )
    pari_err(talker, "incorrect subgrp in bnrstark");

  /* compute bnr(conductor) */
  p1     = conductor(bnr, subgrp, 2);
  bnr    = gel(p1,2); Mcyc = diagonal_i(gmael(bnr, 5, 2));
  subgrp = gel(p1,3);
  if (gcmp1( dethnf_i(subgrp) )) { avma = av; return pol_x[0]; }

  /* check the class field */
  if (!gcmp0(gmael3(bnr, 2, 1, 2)))
    pari_err(talker, "class field not totally real in bnrstark");

  if (DEBUGLEVEL) (void)timer2();

  /* find a suitable extension N */
  dtQ = InitQuotient(subgrp);
  data  = FindModulus(bnr, dtQ, &newprec, prec);

  if (!data)
  {
    GEN vec, H, cyc = gel(dtQ,2), U = gel(dtQ,3), M = ginv(U);
    long i, j = 1, l = lg(M);

    /* M = indep. generators of Cl_f/subgp, restrict to cyclic components */
    vec = cgetg(l, t_VEC);
    for (i = 1; i < l; i++)
    {
      GEN t = gel(M,i);
      if (is_pm1(cyc[i])) continue;
      M[i] = Mcyc[i]; H = hnf(shallowconcat(M, Mcyc));
      gel(M,i) = t;
      gel(vec,j++) = bnrstark(bnr, H, prec);
    }
    setlg(vec, j); return gerepilecopy(av, vec);
  }

  if (newprec > prec)
  {
    if (DEBUGLEVEL>1) fprintferr("new precision: %ld\n", newprec);
    nf = nfnewprec(nf, newprec);
  }

  return gerepileupto(av, AllStark(data, nf, 0, newprec));
}

/* For each character of Cl(bnr)/subgp, compute L(1, chi) (or equivalently
 * the first non-zero term c(chi) of the expansion at s = 0).
 * If flag & 1: compute the value at s = 1 (for non-trivial characters),
 * else compute the term c(chi) and return [r(chi), c(chi)] where r(chi) is
 *   the order of L(s, chi) at s = 0.
 * If flag & 2: compute the value of the L-function L_S(s, chi) where S is the
 *   set of places dividing the modulus of bnr (and the infinite places),
 * else
 *   compute the value of the primitive L-function associated to chi, 
 * If flag & 4: return also the character */
GEN
bnrL1(GEN bnr, GEN subgp, long flag, long prec)
{
  GEN bnf, nf, cyc, Mcyc, L1, lchi, clchi, allCR, listCR, dataCR;
  GEN W, S, T, indCR, invCR, Qt, vChar;
  long N, cl, i, j, nc, a;
  pari_sp av = avma;

  checkbnrgen(bnr);
  bnf  = gel(bnr,1);
  nf   = gel(bnf,7);
  N    = degpol(nf[1]);

  if (N == 1) pari_err(talker, "the ground field must be distinct from Q");
  if (flag < 0 || flag > 8) pari_err(flagerr,"bnrL1");

  /* compute bnr(conductor) */
  if (!(flag & 2)) bnr = (GEN)conductor(bnr, NULL, 2)[2];
  cyc  = gmael(bnr, 5, 2);
  Mcyc = diagonal_i(cyc);

  /* check the subgroup */
  if (! (subgp = get_subgroup(subgp,Mcyc)) )
    pari_err(talker, "incorrect subgroup in bnrL1");

  cl = itou( dethnf_i(subgp) );
  Qt = InitQuotient(subgp);

  /* compute all characters */
  allCR = EltsOfGroup(cl, gel(Qt,2));

  /* make a list of all non-trivial characters modulo conjugation */
  listCR = cgetg(cl, t_VEC);
  indCR = new_chunk(cl);
  invCR = new_chunk(cl); nc = 0;
  for (i = 1; i < cl; i++)
  {
    /* lift to a character on Cl(bnr) */
    lchi = LiftChar(cyc, gel(Qt,3), gel(allCR,i), gel(Qt,2));
    clchi = ConjChar(lchi, cyc);

    a = i;
    for (j = 1; j <= nc; j++)
      if (gequal(gmael(listCR, j, 1), clchi)) { a = -j; break; }

    if (a > 0)
    {
      nc++;
      gel(listCR,nc) = mkvec2(lchi, bnrconductorofchar(bnr, lchi));
      indCR[i]  = nc;
      invCR[nc] = i;
    }
    else
      indCR[i] = -invCR[-a];

    gel(allCR,i) = lchi;
  }

  /* the trivial character has to be a row vector too! */
  settyp(allCR[cl], t_VEC);

  setlg(listCR, nc + 1);
  if (nc == 0) pari_err(talker, "no non-trivial character in bnrL1");

  /* compute the data for these characters */
  dataCR = InitChar(bnr, listCR, prec);

  vChar = sortChars(dataCR);
  GetST(bnr, &S, &T, dataCR, vChar, prec);
  W = ComputeAllArtinNumbers(dataCR, vChar, 1, prec);

  L1 = cgetg((flag&1)? cl: cl+1, t_VEC);
  for (i = 1; i < cl; i++)
  {
    a = indCR[i];
    if (a > 0)
      gel(L1,i) = GetValue(gel(dataCR,a), gel(W,a), gel(S,a), gel(T,a),
                           flag, prec);
    else
      gel(L1,i) = gconj(gel(L1,-a));
  }
  if (!(flag & 1))
    gel(L1,cl) = GetValue1(bnr, flag & 2, prec);
  else
    cl--;

  if (flag & 4) {
    for (i = 1; i <= cl; i++) gel(L1,i) = mkvec2(gel(allCR,i), gel(L1,i));
  }
  return gerepilecopy(av, L1);
}
