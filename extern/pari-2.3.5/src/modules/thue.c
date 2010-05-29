/* $Id: thue.c 11397 2008-12-08 17:08:13Z kb $

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

/* Thue equation solver. In all the forthcoming remarks, "paper"
 * designs the paper "Thue Equations of High Degree", by Yu. Bilu and
 * G. Hanrot, J. Number Theory (1996). Note that numbering of the constants
 * is that of Hanrot's thesis rather than that of the paper.
 * The last part of the program (bnfisintnorm) was written by K. Belabas. */

/* Check whether tnf is a valid structure */
static int
checktnf(GEN tnf)
{
  long n, R, S, T;
  if (typ(tnf)!=t_VEC || (lg(tnf)!=8 && lg(tnf)!=3)) return 0;
  if (typ(tnf[1]) != t_POL) return 0;
  if (lg(tnf) != 8) return 1; /* S=0 */

  n = degpol(tnf[1]);
  if (n <= 2) pari_err(talker,"invalid polynomial in thue (need n>2)");
  S = sturm(gel(tnf,1)); T = (n-S)>>1; R = S+T-1;
  (void)checkbnf(gel(tnf,2));
  if (typ(tnf[3]) != t_COL || lg(tnf[3]) != n+1) return 0;
  if (typ(tnf[4]) != t_COL || lg(tnf[4]) != R+1) return 0;
  if (typ(tnf[5]) != t_MAT || lg(tnf[5])!=R+1 || lg(gmael(tnf,5,1)) != n+1) return 0;
  if (typ(tnf[6]) != t_MAT || lg(tnf[6])!=R+1 || lg(gmael(tnf,6,1)) != R+1) return 0;
  if (typ(tnf[7]) != t_VEC || lg(tnf[7]) != 8) return 0;
  return 1;
}

static GEN
distoZ(GEN z)
{
  GEN p1 = gfrac(z);
  return gmin(p1, gsub(gen_1,p1));
}

/* Compensates rounding errors for computation/display of the constants.
 * Round up if dir > 0, down otherwise */
static GEN
myround(GEN x, long dir)
{
  GEN eps = gpowgs(stoi(dir > 0? 10: -10), -10);
  return gmul(x, gadd(gen_1, eps));
}

/* Returns the index of the largest element in a vector */
static GEN
_Vecmax(GEN Vec, long *ind)
{
  long k, tno = 1, l = lg(Vec);
  GEN tmax = gel(Vec,1);
  for (k = 2; k < l; k++)
    if (gcmp(gel(Vec,k),tmax) > 0) { tmax = gel(Vec,k); tno = k; }
  if (ind) *ind = tno;
  return tmax;
}

static GEN
Vecmax(GEN v) { return _Vecmax(v, NULL); }

static long
Vecmaxind(GEN v) { long i; (void)_Vecmax(v, &i); return i; }

static GEN
tnf_get_roots(GEN poly, long prec, long S, long T)
{
  GEN R0 = roots(poly, prec), R = cgetg(lg(R0), t_COL);
  long k;

  for (k=1; k<=S; k++) gel(R,k) = real_i(gel(R0,k));
  /* swap roots to get the usual order */
  for (k=1; k<=T; k++)
  {
    R[k+S]  = R0[2*k+S-1];
    R[k+S+T]= R0[2*k+S];
  }
  return R;
}

/* Computation of the logarithmic height of x (given by embeddings) */
static GEN
LogHeight(GEN x, long prec)
{
  long i, n = lg(x)-1;
  GEN LH = gen_1;
  for (i=1; i<=n; i++) LH = gmul(LH, gmax(gen_1, gabs(gel(x,i), prec)));
  return gdivgs(glog(LH,prec), n);
}

/* |x|^(1/n), x t_INT */
static GEN
absisqrtn(GEN x, long n, long prec) {
  GEN r = itor(x,prec);
  if (signe(r) < 0) setsigne(r,1);
  return sqrtnr(r, n);
}

static GEN
get_emb(GEN x, GEN r, long prec)
{
  long l = lg(r), i, tx;
  GEN e, y = cgetg(l, t_COL);

  tx = typ(x);
  if (tx != t_INT && tx != t_POL) pari_err(typeer,"get_emb");
  for (i=1; i<l; i++)
  {
    e = poleval(x, gel(r,i));
    if (gcmp0(e) || (typ(e) != t_INT && precision(e) < prec)) return NULL;
    gel(y,i) = e;
  }
  return y;
}

/* Computation of the conjugates (sigma_i(v_j)), and log. heights, of elts of v */
static GEN
Conj_LH(GEN v, GEN *H, GEN r, long prec)
{
  long j, l = lg(v);
  GEN e, M = (GEN)cgetg(l,t_MAT);

  (*H) = cgetg(l,t_COL);
  for (j = 1; j < l; j++)
  {
    if (! (e = get_emb(gel(v,j), r, prec)) ) return NULL; /* FAIL */
    gel(M,j) = e;
    gel(*H,j) = LogHeight(e, prec);
  }
  return M;
}

static GEN abslog(GEN x, long prec) { return gabs(glog(x,prec), prec); }
static GEN logabs(GEN x, long prec) { return glog(gabs(x,prec), prec); }

/* Computation of M, its inverse A and precision check (see paper) */
static GEN
T_A_Matrices(GEN MatFU, long r, GEN *eps5, long prec)
{
  GEN A, p1, m1, IntM, nia, eps3, eps2;
  long e = bit_accuracy(prec);

  m1 = rowslice(vecslice(MatFU, 1,r), 1,r); /* minor order r */
  m1 = logabs(m1,prec);

  A = invmat(m1);
  IntM = gsub(gmul(A,m1), matid(r));

  eps2 = gadd(vecmax(gabs(IntM,prec)), real2n(-e, prec));
  nia = vecmax(gabs(A,prec));

  /* Check for the precision in matrix inversion. See paper, Lemma 2.4.2. */
  p1 = gadd(gmulsg(r, gmul2n(nia, e)), eps2);
  if (gexpo(p1) < -2*r) pari_err(precer,"thue");

  p1 = gadd(gmulsg(r, gmul2n(nia,-e)), eps2);
  eps3 = gmul(gmulsg(2*r*r,nia), p1);
  eps3 = myround(eps3, 1);

  if (DEBUGLEVEL>1) fprintferr("epsilon_3 -> %Z\n",eps3);
  *eps5 = mulsr(r, eps3); return A;
}

/* Performs basic computations concerning the equation.
 * Returns a "tnf" structure containing
 *  1) the polynomial
 *  2) the bnf (used to solve the norm equation)
 *  3) roots, with presumably enough precision
 *  4) The logarithmic heights of units
 *  5) The matrix of conjugates of units
 *  6) its inverse
 *  7) a few technical constants */
static GEN
inithue(GEN P, GEN bnf, long flag, long prec)
{
  GEN MatFU, x0, tnf, tmp, gpmin, dP, csts, ALH, eps5, ro, c1, c2, Ind = gen_1;
  long k,j, n = degpol(P);
  long s,t, prec_roots;

  if (!bnf)
  {
    if (!gcmp1(leading_term(P))) pari_err(talker,"non-monic polynomial in thue");
    bnf = bnfinit0(P,1,NULL,DEFAULTPREC);
    if (flag) (void)certifybuchall(bnf);
    else 
      Ind = gfloor(mulrs(gmael(bnf, 8, 2), 5)); 
  }

  nf_get_sign(checknf(bnf), &s, &t);
  prec_roots = prec; 
  for(;;)
  {
    ro = tnf_get_roots(P, prec_roots, s, t);
    MatFU = Conj_LH(gmael(bnf,8,5), &ALH, ro, prec);
    if (MatFU) break;
    prec_roots = (prec_roots << 1) - 2;
    if (DEBUGLEVEL>1) pari_warn(warnprec, "inithue", prec_roots); 
  }

  dP = derivpol(P);
  c1 = NULL; /* min |P'(r_i)|, i <= s */
  for (k=1; k<=s; k++)
  {
    tmp = gabs(poleval(dP,gel(ro,k)),prec);
    if (!c1 || gcmp(tmp,c1) < 0) c1 = tmp;
  }
  c1 = gdiv(int2n(n-1), c1);
  c1 = gprec_w(myround(c1, 1), DEFAULTPREC);

  c2 = NULL; /* max |r_i - r_j|, i!=j */
  for (k=1; k<=n; k++)
    for (j=k+1; j<=n; j++)
    {
      tmp = gabs(gsub(gel(ro,j),gel(ro,k)), prec);
      if (!c2 || gcmp(c2,tmp) > 0) c2 = tmp;
    }
  c2 = gprec_w(myround(c2, -1), DEFAULTPREC);

  if (t==0) x0 = gen_1;
  else
  {
    gpmin = NULL; /* min |P'(r_i)|, i > s */
    for (k=1; k<=t; k++)
    {
      tmp = gabs(poleval(dP,gel(ro,s+k)), prec);
      if (!gpmin || gcmp(tmp,gpmin) < 0) gpmin = tmp;
    }
    gpmin = gprec_w(gpmin, DEFAULTPREC);

    /* Compute x0. See paper, Prop. 2.2.1 */
    x0 = gmul(gpmin, Vecmax(gabs(imag_i(ro), prec)));
    x0 = sqrtnr(gdiv(int2n(n-1), x0), n);
  }
  if (DEBUGLEVEL>1) 
    fprintferr("c1 = %Z\nc2 = %Z\nIndice <= %Z\n", c1, c2, Ind);

  ALH = gmul2n(ALH, 1);
  tnf = cgetg(8,t_VEC); csts = cgetg(8,t_VEC);
  gel(tnf,1) = P;
  gel(tnf,2) = bnf;
  gel(tnf,3) = ro;
  gel(tnf,4) = ALH;
  gel(tnf,5) = MatFU;
  gel(tnf,6) = T_A_Matrices(MatFU, s+t-1, &eps5, prec);
  gel(tnf,7) = csts;
  gel(csts,1) = c1; gel(csts,2) = c2;   gel(csts,3) = LogHeight(ro, prec);
  gel(csts,4) = x0; gel(csts,5) = eps5; gel(csts,6) = utoipos(prec);
  gel(csts,7) = Ind; 
  return tnf;
}

typedef struct {
  GEN c10, c11, c13, c15, bak, NE, ALH, Ind, hal, MatFU, ro, Hmu;
  GEN delta, lambda, errdelta;
  long r, iroot, deg;
} baker_s;

/* Compute Baker's bound c9 and B_0, the bound for the b_i's. See Thm 2.3.1 */
static GEN
Baker(baker_s *BS)
{
  const long prec = DEFAULTPREC;
  GEN tmp, B0, hb0, c9 = gen_1, ro = BS->ro, ro0 = (GEN)ro[BS->iroot];
  long k, i1, i2, r = BS->r;

  switch (BS->iroot) {
    case 1: i1=2; i2=3; break;
    case 2: i1=1; i2=3; break;
   default: i1=1; i2=2; break;
  }

  /* Compute h_1....h_r */
  for (k=1; k<=r; k++)
  {
    tmp = gdiv(gcoeff(BS->MatFU,i1,k), gcoeff(BS->MatFU,i2,k));
    tmp = gmax(gen_1, abslog(tmp,prec));
    c9 = gmul(c9, gmax((GEN)BS->ALH[k], gdiv(tmp, BS->bak)));
  }

  /* Compute a bound for the h_0 */
  hb0 = gadd(gmul2n(BS->hal,2), gmul2n(gadd(BS->Hmu,mplog2(prec)), 1));
  tmp = gdiv(gmul(gsub(ro0, gel(ro,i2)), (GEN)BS->NE[i1]),
             gmul(gsub(ro0, gel(ro,i1)), (GEN)BS->NE[i2]));
  tmp = gmax(gen_1, abslog(tmp, prec));
  hb0 = gmax(hb0, gdiv(tmp, BS->bak));
  c9 = gmul(c9,hb0);
  /* Multiply c9 by the "constant" factor */
  c9 = gmul(c9, gmul(mulri(mulsr(18,mppi(prec)), int2n(5*(4+r))),
                     gmul(gmul(mpfact(r+3), powiu(mulis(BS->bak,r+2), r+3)),
                          glog(mulis(BS->bak,2*(r+2)),prec))));
  c9 = gprec_w(myround(c9, 1), DEFAULTPREC);
  /* Compute B0 according to Lemma 2.3.3 */
  B0 = mulir(mulsi(2, BS->Ind), 
	     divrr(addrr(mulrr(c9,mplog(divrr(mulir(BS->Ind, c9),BS->c10))), 
			 mplog(mulir(BS->Ind, BS->c11))),
		   BS->c10));
  B0 = gmax(B0, dbltor(2.71828183));
  B0 = gmax(B0, mulrr(divir(BS->Ind, BS->c10), 
		      mplog(divrr(mulir(BS->Ind, BS->c11), 
				  Pi2n(1, prec))))); 

  if (DEBUGLEVEL>1) {
    fprintferr("  B0  = %Z\n",B0);
    fprintferr("  Baker = %Z\n",c9);
  }
  return B0;
}

/* || x d ||, x t_REAL, d t_INT */
static GEN
errnum(GEN x, GEN d)
{
  GEN dx = mulir(d, x);
  return mpabs(subri(dx, ground(dx)));
}

/* Try to reduce the bound through continued fractions; see paper. */
static int
CF_1stPass(GEN *B0, GEN kappa, baker_s *BS)
{
  GEN q, ql, qd, l0, denbound = mulri(*B0, kappa);
  
  if (gcmp(gmul(dbltor(0.1),gsqr(denbound)), ginv(BS->errdelta)) > 0) 
    return -1;

  q = denom( bestappr(BS->delta, denbound) );
  qd = errnum(BS->delta, q);
  ql = errnum(BS->lambda,q);

  l0 = subrr(ql, addrr(mulrr(qd, *B0), divri(dbltor(0.1),kappa)));
  if (signe(l0) <= 0) return 0;

  if (BS->r > 1)
    *B0 = divrr(mplog(divrr(mulir(q,BS->c15), l0)), BS->c13);
  else
  {
    l0 = mulrr(l0,Pi2n(1, DEFAULTPREC));
    *B0 = divrr(mplog(divrr(mulir(q,BS->c11), l0)), BS->c10);
  }
  if (DEBUGLEVEL>1) fprintferr("    B0 -> %Z\n",*B0);
  return 1;
}

static int
LLL_1stPass(GEN *pB0, GEN kappa, baker_s *BS, GEN *pBx)
{
  GEN B0 = *pB0, Bx = *pBx, lllmat, C, l0, l1, delta, lambda, triv; 
  long e;

  delta = BS -> delta; 
  lambda = BS -> lambda; 
  
  C = grndtoi(mulir(mulii(BS->Ind, kappa), 
		    gpow(B0, dbltor(2.2), DEFAULTPREC)), &e);

  if (DEBUGLEVEL > 1) fprintferr("C (bitsize) : %d\n", expi(C)); 
  lllmat = matid(3);
  if (gcmp(B0, BS -> Ind) > 0) 
  {
    gcoeff(lllmat, 1, 1) = grndtoi(divri(B0, BS -> Ind), &e); 
    triv = mulsr(2, gsqr(B0)); 
  }
  else 
    triv = addir(gsqr(BS -> Ind), gsqr(B0)); 

  gcoeff(lllmat, 3, 1) = ground(gneg(gmul(C, lambda)));
  gcoeff(lllmat, 3, 2) = ground(gneg(gmul(C, delta)));
  gcoeff(lllmat, 3,3) = C;
  lllmat = gmul(lllmat, lllint(lllmat));

  l0 = gnorml2(gel(lllmat,1));
  l0 = subrr(divir(l0, dbltor(1.8262)), triv); /* delta = 0.99 */
  if (signe(l0) <= 0) return 0; 

  l1 = divrs(addri(mulsr(2, B0), BS->Ind), 2);
  l0 = divri(subrr(sqrtr(l0), l1), C);
 
  if (signe(l0) <= 0) return 0;
  B0 = gmul(gdiv(BS->Ind, BS->c13), 
            mplog(gdiv(gmul(BS->Ind, BS->c15), l0)));
  Bx = gpow(gdiv(mulsr(2, gmul(BS->Ind, BS->c15)), l0), 
            ginv(utoipos(BS->deg)), DEFAULTPREC);

  if (DEBUGLEVEL>=2)
    {
      fprintferr("LLL_First_Pass successful !!\n");
      fprintferr("B0 -> %Z\n", B0);
      fprintferr("x <= %Z\n", Bx);
    }

  *pB0 = B0; *pBx = Bx; return 1;
}


/* Check whether a solution has already been found */
static int
new_sol(GEN z, GEN S)
{
  long i, l = lg(S);
  for (i=1; i<l; i++)
    if (gequal(z,gel(S,i))) return 0;
  return 1;
}

static void
add_sol(GEN *pS, GEN x, GEN y)
{
  GEN u = mkvec2(x,y);
  if (new_sol(u, *pS)) *pS = shallowconcat(*pS, mkvec(u));
}

static void
check_sol(GEN x, GEN y, GEN P, GEN rhs, GEN *pS)
{
  if (gcmp0(y))
  { if (equalii(powiu(x,degpol(P)), rhs)) add_sol(pS, x, gen_0); }
  else
  { if (gequal(poleval(RgX_rescale(P,y),x), rhs)) add_sol(pS, x, y); }
}

/* Check whether a potential solution is a true solution. Return 0 if
 * truncation error (increase precision) */
static int
CheckSol(GEN *pS, GEN z1, GEN z2, GEN P, GEN rhs, GEN ro)
{
  GEN x, y, ro1 = gel(ro,1), ro2 = gel(ro,2);
  long e;

  y = grndtoi(real_i(gdiv(gsub(z2,z1), gsub(ro1,ro2))), &e);
  if (e > 0) return 0;
  x = gadd(z1, gmul(ro1, y));
  x = grndtoi(real_i(x), &e);
  if (e > 0) return 0;
  if (e <= -13)
  {
    check_sol(     x ,      y,  P,rhs,pS);
    check_sol(negi(x), negi(y), P,rhs,pS);
  }
  return 1;
}

/* find q1,q2,q3 st q1 + b q2 + c q3 ~ 0 */
static GEN
GuessQi(GEN b, GEN c, GEN *eps)
{
  GEN Q, Lat, C = int2n(33);

  Lat = matid(3);
  gmael(Lat,1,3) = C;
  gmael(Lat,2,3) = ground(gmul(C,b));
  gmael(Lat,3,3) = ground(gmul(C,c));

  Q = (GEN)lllint(Lat)[1];
  if (gcmp0(gel(Q,3))) return NULL; /* FAIL */

  *eps = gadd(gadd(gel(Q,1), gmul(gel(Q,2),b)), gmul(gel(Q,3),c));
  *eps = mpabs(*eps); return Q;
}

/* Check for not-so-small solutions */
static GEN
MiddleSols(GEN *pS, GEN bound, GEN roo, GEN poly, GEN rhs, long s, GEN c1) 
{
  long j, k, nmax, d = degpol(poly); 
  GEN bndcf = sqrtnr(shiftr(c1,1), d - 2); 

  if (cmprr(bound, bndcf) == -1) return bound; 
  /* divide by log((1+sqrt(5))/2) 
   * 1 + ==> ceil
   * 2 + ==> continued fraction is normalized if last entry is 1 */
  nmax = 2 + (long)(gtodouble(mplog(bound)) / 0.4812118250596);
  if (nmax < 3) nmax = 3;

  for (k = 1; k <= s; k++) 
  {
    GEN t = contfrac0(real_i(gel(roo,k)), NULL, nmax); 
    GEN p, q, pm1, qm1, p0, q0, z; 

    pm1 = gen_0; p0 = gen_1; 
    qm1 = gen_1; q0 = gen_0; 

    for (j = 1; j < lg(t); j++) 
    {
      p = addii(mulii(p0, gel(t,j)), pm1); 
      pm1 = p0; p0 = p; 

      q = addii(mulii(q0, gel(t,j)), qm1); 
      qm1 = q0; q0 = q; 
      
      if (cmpir(q, bound) > 0) break; 
      if (DEBUGLEVEL >= 2)
        fprintferr("Checking (\\pm %Z, \\pm %Z)\n",p, q);

      z = poleval(RgX_rescale(poly,q), p); /* = P(p/q) q^dep(P) */
      if (absi_equal(z, rhs)) 	    
      {
        if (signe(z) == signe(rhs))
        {
          add_sol(pS, p, q); 
          if (d % 2 == 0) add_sol(pS, negi(p), negi(q)); 
        }
        else
          if (d % 2 == 1) add_sol(pS, negi(p), negi(q)); 
      }
    }
    if (j == lg(t)) pari_err(talker, "Not enough precision in thue"); 
  }
  return bndcf;
}

/* Check for solutions under a small bound (see paper) */
static GEN
SmallSols(GEN S, long Bx, GEN poly, GEN rhs, GEN ro)
{
  pari_sp av = avma, lim = stack_lim(av, 1);
  const long prec = DEFAULTPREC;
  GEN X, Y, P, r;
  long x, j, n = degpol(poly);

  if (DEBUGLEVEL>1) fprintferr("* Checking for small solutions\n");
  /* x = 0 first */
  Y = ground( absisqrtn(rhs, n, prec) );
  if (gequal(powiu(Y,n), rhs)) add_sol(&S, Y, gen_0);
  Y = negi(Y);
  if (gequal(powiu(Y,n), rhs)) add_sol(&S, Y, gen_0);

  /* x != 0 */
  P = cgetg(lg(poly), t_POL); P[1] = poly[1]; 
  for (x = -Bx; x <= Bx; x++)
  {
    if (!x) continue;

    X = stoi(x); 
    P[lg(poly) - 1] = poly[lg(poly) - 1]; 
    for (j = lg(poly) - 2; j >= 2; j--) 
    {
      gel(P,j) = mulii(X, gel(poly,j));
      X = mulis(X, x); 
    }
    gel(P,2) = subii(gel(P,2), rhs); 
    r = nfrootsQ(P); 

    for (j = 1; j < lg(r); j++) 
      if (typ(gel(r,j)) == t_INT) add_sol(&S, gel(r,j), stoi(x)); 
    if (low_stack(lim,stack_lim(av,1)))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"SmallSols");
      S = gerepilecopy(av, S);
      P = cgetg(lg(poly), t_POL); P[1] = poly[1]; 
    }
  }
  return S;
}

/* Computes [x]! */
static double
fact(double x)
{
  double ft = 1.0;
  x = floor(x); while (x>1) { ft *= x; x--; }
  return ft ;
}

/* From a polynomial and optionally a system of fundamental units for the
 * field defined by pol, computes all relevant constants needed to solve
 * the equation P(x,y)=a given the solutions of N_{K/Q}(x)=a (see inithue). */
GEN
thueinit(GEN pol, long flag, long prec)
{
  GEN tnf, bnf = NULL;
  pari_sp av = avma;
  long k, s;

  if (checktnf(pol)) { bnf = checkbnf(gel(pol,2)); pol = gel(pol,1); }
  if (typ(pol)!=t_POL) pari_err(notpoler,"thueinit");
  if (degpol(pol) <= 2) pari_err(talker,"invalid polynomial in thue (need deg>2)");

  s = sturm(pol);
  if (s)
  {
    long PREC, n = degpol(pol);
    double d, dr, dn = (double)n;

    dr = (double)((s+n-2)>>1); /* s+t-1 */
    d = dn*(dn-1)*(dn-2);
    /* Guess precision by approximating Baker's bound. The guess is most of
     * the time not sharp, ie 10 to 30 decimal digits above what is _really_
     * necessary. Note that the limiting step is the reduction. See paper. */
    PREC = 3 + (long)((5.83 + (dr+4)*5 + log(fact(dr+3)) + (dr+3)*log(dr+2) +
		     (dr+3)*log(d) + log(log(2*d*(dr+2))) + (dr+1))
                     / ((BYTES_IN_LONG/4)* 10.));

    if (flag == 0) PREC = (long)(2.2 * PREC); /* Lazy, to be improved */
    if (PREC < prec) PREC = prec;
    if (DEBUGLEVEL >=2) fprintferr("prec = %d\n", PREC); 

    for (;;)
    {
      if (( tnf = inithue(pol, bnf, flag, PREC) )) break;
      PREC = (PREC<<1)-2;
      if (DEBUGLEVEL>1) pari_warn(warnprec,"thueinit",PREC);
      bnf = NULL; avma = av;
    }
  }
  else
  {
    GEN c0 = gen_1, ro = roots(pol, DEFAULTPREC);
    if (!gisirreducible(pol)) pari_err(redpoler,"thueinit");
    for (k=1; k<lg(ro); k++) c0 = gmul(c0, imag_i(gel(ro,k)));
    c0 = ginv( mpabs(c0) );
    tnf = mkvec2(pol, c0);
  }
  return gerepilecopy(av,tnf);
}

static void
init_get_B(long i1, long i2, GEN Delta, GEN Lambda, GEN eps5, baker_s *BS,
           long prec)
{
  GEN delta, lambda, errdelta;
  if (BS->r > 1)
  {
    delta = divrr(gel(Delta,i2),gel(Delta,i1));
    lambda = gdiv(gsub(gmul(gel(Delta,i2),gel(Lambda,i1)),
                       gmul(gel(Delta,i1),gel(Lambda,i2))),
                  gel(Delta,i1));
    errdelta = mulrr(addsr(1,delta),
                     divrr(eps5, subrr(mpabs(gel(Delta,i1)),eps5)));
  }
  else
  { /* r == 1, single fundamental unit (i1 = s = t = 1) */
    GEN p1, Pi2 = Pi2n(1, prec);
    GEN fu = (GEN)BS->MatFU[1], ro = BS->ro;

    p1 = gdiv(gel(fu,2), gel(fu,3));
    delta = divrr(garg(p1,prec), Pi2);

    p1 = gmul(gdiv(gsub(gel(ro,1), gel(ro,2)),
                   gsub(gel(ro,1), gel(ro,3))),
              gdiv((GEN)BS->NE[3], (GEN)BS->NE[2]));
    lambda = divrr(garg(p1,prec), Pi2);

    errdelta = ginv(gmul2n(gabs(gel(fu,2),prec), bit_accuracy(prec)-1));
  }
  if (DEBUGLEVEL>1) fprintferr("  errdelta = %Z\n",errdelta);
  BS->delta = delta;
  BS->lambda = lambda;
  BS->errdelta = errdelta;
}

static GEN
get_B0(long i1, GEN Delta, GEN Lambda, GEN eps5, long prec, baker_s *BS)
{
  GEN B0 = Baker(BS);
  long i2 = (i1 == 1)? 2: 1;
  for(;;) /* i2 from 1 to r unless r = 1 [then i2 = 2] */
  {
    init_get_B(i1,i2, Delta,Lambda,eps5, BS, prec);
    if (DEBUGLEVEL>1) fprintferr("  Entering CF...\n");
    /* Reduce B0 as long as we make progress: newB0 < oldB0 - 0.1 */
    for (;;)
    {
      GEN oldB0 = B0, kappa = utoipos(10);
      long cf;

      for (cf = 0; cf < 10; cf++, kappa = mulis(kappa,10))
      {
        int res = CF_1stPass(&B0, kappa, BS);
        if (res < 0) return NULL; /* prec problem */
        if (res) break;
        if (DEBUGLEVEL>1) fprintferr("CF failed. Increasing kappa\n");
      }
      if (cf == 10)
      { /* Semirational or totally rational case */
        GEN Q, ep, q, l0, denbound;

        if (! (Q = GuessQi(BS->delta, BS->lambda, &ep)) ) break;

        denbound = gadd(B0, absi(gel(Q,2)));
        q = denom( bestappr(BS->delta, denbound) );
        l0 = subrr(errnum(BS->delta, q), ep);
        if (signe(l0) <= 0) break;

        B0 = divrr(mplog(divrr(mulir(gel(Q,3), BS->c15), l0)),  BS->c13);
        if (DEBUGLEVEL>1) fprintferr("Semirat. reduction: B0 -> %Z\n",B0);
      }
      /* if no progress, stop */
      if (gcmp(oldB0, gadd(B0,dbltor(0.1))) <= 0) return gmin(oldB0, B0);
    }
    i2++; if (i2 == i1) i2++;
    if (i2 > BS->r) break;
  }
  pari_err(bugparier,"thue (totally rational case)");
  return NULL; /* not reached */
}

static GEN
get_Bx_LLL(long i1, GEN Delta, GEN Lambda, GEN eps5, long prec, baker_s *BS)
{
  GEN B0 = Baker(BS), Bx = NULL;
  long i2 = (i1 == 1)? 2: 1;
  for(;;) /* i2 from 1 to r unless r = 1 [then i2 = 2] */
  {
    init_get_B(i1,i2, Delta,Lambda,eps5, BS, prec);
    if (DEBUGLEVEL>1) fprintferr("  Entering LLL...\n");
    /* Reduce B0 as long as we make progress: newB0 < oldB0 - 0.1 */
    for (;;)
    {
      GEN oldBx = Bx, kappa = utoipos(10);
      long cf;

      for (cf = 0; cf < 10; cf++, kappa = mulis(kappa,10))
      {
        int res = LLL_1stPass(&B0, kappa, BS, &Bx);
        if (res) break;
        if (DEBUGLEVEL>1) fprintferr("LLL failed. Increasing kappa\n");
      }

      /* TO COMPLETE */
      if (cf == 10)
      { /* Semirational or totally rational case */
        GEN Q, ep, q, l0, denbound;

        if (! (Q = GuessQi(BS->delta, BS->lambda, &ep)) ) break;

        /* Beware Q[2]] = gen_0 */
        denbound = gadd(mulri(B0, absi(gel(Q,2))), 
			mulii(BS->Ind, absi(gel(Q,3))));
        q = denom( bestappr(BS->delta, denbound) );
        l0 = divri(subrr(errnum(BS->delta, q), ep), absi(gel(Q,3))); 
        if (signe(l0) <= 0) break;

	B0 = divrr(mulir(BS->Ind, mplog(divrr(mulir(BS->Ind, BS->c15), l0))),
		   BS->c13);
	Bx = gpow(gdiv(mulsr(2, gmul(BS->Ind, BS->c15)), l0), 
		  ginv(utoipos(BS->deg)), DEFAULTPREC);

        if (DEBUGLEVEL>1) 
	  fprintferr("Semirat. reduction: B0 -> %Z x <= %Z\n",B0, Bx);
      }
      /* if no progress, stop */
      if (oldBx && gcmp(oldBx, Bx) <= 0) return oldBx;
    }
    i2++; if (i2 == i1) i2++;
    if (i2 > BS->r) break;
  }
  pari_err(bugparier,"thue (totally rational case)");
  return NULL; /* not reached */
}

static GEN
LargeSols(GEN tnf, GEN rhs, GEN ne, GEN *pro, GEN *pS)
{
  GEN Vect, P, ro, bnf, MatFU, A, csts, dP, vecdP, Bx;
  GEN c1,c2,c3,c4,c10,c11,c13,c14,c15, x0, x1, x2, x3, b, zp1, tmp, eps5, Ind;
  long iroot, ine, n, i, r, upb, bi1, Prec, prec, s,t;
  baker_s BS;
  pari_sp av = avma;

  bnf  = gel(tnf,2);
  if (!ne)
  { 
    ne = bnfisintnorm(bnf, rhs);
    if (!gcmp1(gmael(tnf, 7, 7)) &&
        !gcmp1(gmael3(bnf, 8, 1, 1)) && !is_pm1(rhs))
      pari_warn(warner, "Non trivial conditional class group.\n  *** May miss solutions of the norm equation"); 
  }
  if (lg(ne)==1) return NULL;

  nf_get_sign(checknf(bnf), &s, &t);
  BS.r = r = s+t-1;
  P      = gel(tnf,1); n = degpol(P);
  ro     = gel(tnf,3);
  BS.ALH = gel(tnf,4);
  MatFU  = gel(tnf,5);
  A      = gel(tnf,6);
  csts   = gel(tnf,7);
  c1     = gel(csts,1); c1 = gmul(absi(rhs), c1);
  c2     = gel(csts,2);
  BS.hal = gel(csts,3);
  x0     = gel(csts,4);
  eps5   = gel(csts,5);
  Prec = gtolong(gel(csts,6));
  Ind    = gel(csts,7); 
  BS.MatFU = MatFU;
  BS.bak = mulss(n, (n-1)*(n-2)); /* safe */
  BS.deg = n; 
  *pS = cgetg(1, t_VEC);

  if (t) x0 = gmul(x0, absisqrtn(rhs, n, Prec));
  tmp = divrr(c1,c2);
  c3 = mulrr(dbltor(1.39), tmp);
  c4 = mulsr(n-1, c3);
  x1 = gmax(x0, sqrtnr(mulsr(2,tmp),n));

  Vect = cgetg(r+1,t_COL); for (i=1; i<=r; i++) gel(Vect,i) = gen_1;
  Vect = gmul(gabs(A,DEFAULTPREC), Vect);
  c14 = mulrr(c4, Vecmax(Vect));
  x2 = gmax(x1, sqrtnr(mulsr(10,c14), n));
  if (DEBUGLEVEL>1) {
    fprintferr("x1 -> %Z\n",x1);
    fprintferr("x2 -> %Z\n",x2);
    fprintferr("c14 = %Z\n",c14);
  }

  dP = derivpol(P);
  vecdP = cgetg(s+1, t_VEC);
  for (i=1; i<=s; i++) gel(vecdP,i) = poleval(dP, gel(ro,i));

  zp1 = dbltor(0.01);
  x3 = gmax(x2, sqrtnr(mulsr(2,divrr(c14,zp1)),n));

  b = cgetg(r+1,t_COL);
  for (iroot=1; iroot<=s; iroot++)
  {
    GEN Delta, MatNE, Hmu, c5, c7;

    Vect = cgetg(r+1,t_COL); for (i=1; i<=r; i++) gel(Vect,i) = gen_1;
    if (iroot <= r) gel(Vect,iroot) = stoi(1-n);
    Delta = gmul(A,Vect);

    c5 = Vecmax(gabs(Delta,Prec));
    c5  = myround(gprec_w(c5,DEFAULTPREC), 1);
    c7  = mulsr(r,c5);
    c10 = divsr(n,c7); BS.c10 = c10;
    c13 = divsr(n,c5); BS.c13 = c13;
    if (DEBUGLEVEL>1) {
      fprintferr("* real root no %ld/%ld\n", iroot,s);
      fprintferr("  c10 = %Z\n",c10);
      fprintferr("  c13 = %Z\n",c13);
    }

    prec = Prec;
    for (;;)
    {
      if (( MatNE = Conj_LH(ne, &Hmu, ro, prec) )) break;
      prec = (prec<<1)-2;
      if (DEBUGLEVEL>1) pari_warn(warnprec,"thue",prec);
      ro = tnf_get_roots(P, prec, s, t);
    }
    BS.ro    = ro;
    BS.iroot = iroot;

    for (ine=1; ine<lg(ne); ine++)
    {
      GEN Lambda, B0, c6, c8;
      GEN NE = gel(MatNE,ine), Vect2 = cgetg(r+1,t_COL);
      long k, i1;

      if (DEBUGLEVEL>1) fprintferr("  - norm sol. no %ld/%ld\n",ine,lg(ne)-1);
      for (k=1; k<=r; k++)
      {
        if (k == iroot)
          tmp = gdiv(rhs, gmul(gel(vecdP,k), gel(NE,k)));
        else
          tmp = gdiv(gsub(gel(ro,iroot),gel(ro,k)), gel(NE,k));
        gel(Vect2,k) = glog(gabs(tmp,prec), prec);
      }
      Lambda = gmul(A,Vect2);

      c6 = addrr(dbltor(0.1), Vecmax(gabs(Lambda,DEFAULTPREC)));
      c6 = myround(c6, 1);
      c8 = addrr(dbltor(1.23), mulsr(r,c6));
      c11= mulrr(mulsr(2,c3) , mpexp(divrr(mulsr(n,c8),c7)));
      c15= mulrr(mulsr(2,c14), mpexp(divrr(mulsr(n,c6),c5)));

      if (DEBUGLEVEL>1) {
        fprintferr("  c6  = %Z\n",c6);
        fprintferr("  c8  = %Z\n",c8);
        fprintferr("  c11 = %Z\n",c11);
        fprintferr("  c15 = %Z\n",c15);
      }
      BS.c11 = c11;
      BS.c15 = c15;
      BS.NE = NE;
      BS.Hmu = gel(Hmu,ine);
      BS.Ind = Ind; 

      i1 = Vecmaxind(gabs(Delta,prec));
      if (gcmp1(Ind)) 
	{ 
	  if (! (B0 = get_B0(i1, Delta, Lambda, eps5, prec, &BS)) ) 
	    goto PRECPB;
	}
      else
	{
	  if (! (Bx = get_Bx_LLL(i1, Delta, Lambda, eps5, prec, &BS)) )
	    goto PRECPB; 
	  x3 = gmax(Bx, x3); 
	  continue;
	}
     /* For each possible value of b_i1, compute the b_i's
      * and 2 conjugates of z = x - alpha y. Then check. */
      upb = gtolong(gceil(B0));
      for (bi1=-upb; bi1<=upb; bi1++)
      {
        GEN z1, z2;
        for (i=1; i<=r; i++)
        {
          gel(b,i) = gdiv(gsub(gmul(gel(Delta,i), stoi(bi1)),
                           gsub(gmul(gel(Delta,i),gel(Lambda,i1)),
                                gmul(gel(Delta,i1),gel(Lambda,i)))),
                      gel(Delta,i1));
          if (gcmp(distoZ(gel(b,i)), zp1) > 0) break;
        }
        if (i <= r) continue;

        z1 = z2 = gen_1;
        for(i=1; i<=r; i++)
        {
          GEN c = ground(gel(b,i));
          z1 = gmul(z1, powgi(gcoeff(MatFU,1,i), c));
          z2 = gmul(z2, powgi(gcoeff(MatFU,2,i), c));
        }
        z1 = gmul(z1, gel(NE,1));
        z2 = gmul(z2, gel(NE,2));
        if (!CheckSol(pS, z1,z2,P,rhs,ro)) goto PRECPB;
      }
    }
  }
  *pro = ro; return MiddleSols(pS, x3, ro, P, rhs, s, c1); 
  
PRECPB:
  ne = gerepilecopy(av, ne);
  prec += 5 * (DEFAULTPREC-2);
  if (DEBUGLEVEL>1) pari_warn(warnprec,"thue",prec);
  tnf = inithue(P, bnf, 0, prec);
  return LargeSols(tnf, rhs, ne, pro, pS);
}

/* Given a tnf structure as returned by thueinit, a RHS and
 * optionally the solutions to the norm equation, returns the solutions to
 * the Thue equation F(x,y)=a
 */
GEN
thue(GEN tnf, GEN rhs, GEN ne)
{
  pari_sp av = avma;
  GEN P, ro, x3, S;

  if (!checktnf(tnf)) pari_err(talker,"not a tnf in thue");
  if (typ(rhs) != t_INT) pari_err(typeer,"thue");

  P = gel(tnf,1);
  if (lg(tnf) == 8)
  {
    x3 = LargeSols(tnf, rhs, ne, &ro, &S);
    if (!x3) { avma = av; return cgetg(1,t_VEC); }
  }
  else
  { /* Case s=0. All solutions are "small". */
    GEN c0   = gel(tnf,2); /* t_REAL */
    S = cgetg(1,t_VEC);
    ro = roots(P, DEFAULTPREC);
    x3 = sqrtnr(mulir(absi(rhs),c0), degpol(P));
    x3 = addrr(x3, dbltor(0.1)); /* guard from round-off errors */
  }

  if (DEBUGLEVEL>=2)
    fprintferr("All solutions are <= %Z\n", x3);
  
  return gerepilecopy(av, SmallSols(S, itos(gfloor(x3)), P, rhs, ro));
}

static GEN *Relations; /* primes above a, expressed on generators of Cl(K) */
static GEN *Partial;   /* list of vvectors, Partial[i] = Relations[1..i] * u[1..i] */
static GEN *gen_ord;   /* orders of generators of Cl(K) given in bnf */

static long *f;        /* f[i] = f(Primes[i]/p), inertial degree */
static long *n;        /* a = prod p^{ n_p }. n[i]=n_p if Primes[i] divides p */
static long *inext;    /* index of first P above next p, 0 if p is last */
static long *S;        /* S[i] = n[i] - sum_{ 1<=k<=i } f[k].u[k] */
static long *u;        /* We want principal ideals I = prod Primes[i]^u[i] */
static GEN  *normsol; /* lists of copies of the u[] which are solutions */

static long Nprimes; /* length(Relations) = #{max ideal above divisors of a} */
static long sindex, smax; /* current index in normsol; max. index */

/* u[1..i] has been filled. Norm(u) is correct.
 * Check relations in class group then save it.
 */
static void
test_sol(long i)
{
  long k,*sol;

  if (Partial)
  {
    pari_sp av=avma;
    for (k=1; k<lg(Partial[1]); k++)
      if ( signe(modii( gmael(Partial,i,k), gen_ord[k] )) )
        { avma=av; return; }
    avma=av;
  }
  if (sindex == smax)
  {
    long new_smax = smax << 1;
    GEN *new_normsol = (GEN*)new_chunk(new_smax+1);

    for (k=1; k<=smax; k++) new_normsol[k] = normsol[k];
    normsol = new_normsol; smax = new_smax;
  }
  sol = cgetg(Nprimes+1,t_VECSMALL);
  normsol[++sindex] = sol;

  for (k=1; k<=i; k++)       sol[k] = u[k];
  for (   ; k<=Nprimes; k++) sol[k] = 0;
  if (DEBUGLEVEL>2)
  {
    fprintferr("sol = %Z\n",sol);
    if (Partial) fprintferr("Partial = %Z\n",Partial);
    flusherr();
  }
}
static void
fix_Partial(long i)
{
  long k;
  pari_sp av = avma;
  for (k=1; k<lg(Partial[1]); k++)
    affii(addii(gmael(Partial,i-1,k), mulis(gmael(Relations,i,k), u[i])),
          gmael(Partial,i,k));
  avma = av;
}

/* Recursive loop. Suppose u[1..i] has been filled
 * Find possible solutions u such that, Norm(prod Prime[i]^u[i]) = a, taking
 * into account:
 *  1) the relations in the class group if need be.
 *  2) the factorization of a.
 */
static void
isintnorm_loop(long i)
{
  if (S[i] == 0) /* sum u[i].f[i] = n[i], do another prime */
  {
    long k;
    if (inext[i] == 0) { test_sol(i); return; }

    /* there are some primes left */
    if (Partial) gaffect(Partial[i], Partial[inext[i]-1]);
    for (k=i+1; k < inext[i]; k++) u[k]=0;
    i=inext[i]-1;
  }
  else if (i == inext[i]-2 || i == Nprimes-1)
  {
    /* only one Prime left above prime; change prime, fix u[i+1] */
    if (S[i] % f[i+1]) return;
    i++; u[i] = S[i-1] / f[i];
    if (Partial) fix_Partial(i);
    if (inext[i]==0) { test_sol(i); return; }
  }

  i++; u[i]=0;
  if (Partial) gaffect(Partial[i-1], Partial[i]);
  if (i == inext[i-1])
  { /* change prime */
    if (inext[i] == i+1 || i == Nprimes) /* only one Prime above p */
    {
      S[i]=0; u[i] = n[i] / f[i]; /* we already know this is exact */
      if (Partial) fix_Partial(i);
    }
    else S[i] = n[i];
  }
  else S[i] = S[i-1]; /* same prime, different Prime */
  for(;;)
  {
    isintnorm_loop(i); S[i]-=f[i]; if (S[i]<0) break;
    if (Partial)
      gaddz(Partial[i], Relations[i], Partial[i]);
    u[i]++;
  }
}

static void
get_sol_abs(GEN bnf, GEN a, GEN *ptPrimes)
{
  GEN dec, fact, primes, Primes, *Fact;
  long *gcdlist, gcd,nprimes,Ngen,i,j;

  *ptPrimes = NULL;
  if (gcmp1(a))
  {
    GEN sol = cgetg(Nprimes+1, t_VECSMALL);
    sindex = 1; normsol = (GEN*) new_chunk(2);
    normsol[1] = sol; for (i=1; i<=Nprimes; i++) sol[i] = 0;
    return;
  }

  fact=factor(a); primes=gel(fact,1);
  nprimes=lg(primes)-1; sindex = 0;
  gen_ord = (GEN *) mael3(bnf,8,1,2); Ngen = lg(gen_ord)-1;

  Fact = (GEN*) new_chunk(1+nprimes);
  gcdlist = new_chunk(1+nprimes);

  Nprimes=0;
  for (i=1; i<=nprimes; i++)
  {
    long ldec;

    dec = primedec(bnf,gel(primes,i)); ldec = lg(dec)-1;
    gcd = itos(gmael(dec,1,4));

    /* check that gcd_{P|p} f_P | n_p */
    for (j=2; gcd>1 && j<=ldec; j++)
      gcd = cgcd(gcd,itos(gmael(dec,j,4)));

    gcdlist[i]=gcd;

    if (gcd != 1 && smodis(gmael(fact,2,i),gcd))
    {
      if (DEBUGLEVEL>2)
        { fprintferr("gcd f_P  does not divide n_p\n"); flusherr(); }
      return;
    }
    Fact[i] = dec; Nprimes += ldec;
  }

  f = new_chunk(1+Nprimes); u = new_chunk(1+Nprimes);
  n = new_chunk(1+Nprimes); S = new_chunk(1+Nprimes);
  inext = new_chunk(1+Nprimes);
  Primes = cgetg(1+Nprimes, t_VEC);
  *ptPrimes = Primes;

  if (Ngen)
  {
    Partial   = (GEN*) new_chunk(1+Nprimes);
    Relations = (GEN*) new_chunk(1+Nprimes);
  }
  else /* trivial class group, no relations to check */
    Relations = Partial = NULL;

  j=0;
  for (i=1; i<=nprimes; i++)
  {
    long k,lim,v, vn=itos(gmael(fact,2,i));

    gcd = gcdlist[i];
    dec = Fact[i]; lim = lg(dec);
    v = (i==nprimes)? 0: j + lim;
    for (k=1; k < lim; k++)
    {
      j++; Primes[j] = dec[k];
      f[j] = itos(gmael(dec,k,4)) / gcd;
      n[j] = vn / gcd; inext[j] = v;
      if (Partial)
	Relations[j] = isprincipal(bnf, gel(Primes,j));
    }
  }
  if (Partial)
  {
    for (i=1; i <= Nprimes; i++)
      if (!gcmp0(Relations[i])) break;
    if (i > Nprimes) Partial = NULL; /* all ideals dividing a are principal */
  }
  if (Partial)
    for (i=0; i<=Nprimes; i++) /* Partial[0] needs to be initialized */
    {
      Partial[i]=cgetg(Ngen+1,t_COL);
      for (j=1; j<=Ngen; j++)
      {
        GEN z = cgeti(4); z[1] = evalsigne(0)|evallgefint(4);
	Partial[i][j]=(long)z;
      }
    }
  smax=511; normsol = (GEN*) new_chunk(smax+1);
  S[0]=n[1]; inext[0]=1; isintnorm_loop(0);
}

/* Look for unit of norm -1. Return 1 if it exists and set *unit, 0 otherwise */
static long
get_unit_1(GEN bnf, GEN *unit)
{
  GEN v, nf = checknf(bnf);
  long i, n = degpol(nf[7]);

  if (DEBUGLEVEL > 2) fprintferr("looking for a fundamental unit of norm -1\n");
  if (odd(n)) { *unit = gen_m1; return 1; }
  v = zsignunits(bnf, NULL, 0);
  for (i = 1; i < lg(v); i++)
  {
    GEN s = sum(gel(v,i), 1, lg(v[i])-1);
    if (mpodd(s)) {
      GEN fu = check_units(bnf, "bnfisintnorm");
      *unit = gel(fu,i); return 1;
    }
  }
  return 0;
}

GEN
bnfisintnormabs(GEN bnf, GEN a)
{
  GEN nf, res, x, Primes;
  long i;

  bnf = checkbnf(bnf); nf = gel(bnf,7);
  if (typ(a)!=t_INT)
    pari_err(talker,"expected an integer in bnfisintnorm");
  if (!signe(a))  return mkvec(gen_0);
  if (gcmp1(a)) return mkvec(gen_1);

  get_sol_abs(bnf, absi(a), &Primes);

  res = cget1(sindex+1, t_VEC);
  for (i=1; i<=sindex; i++)
  {
    x = normsol[i];
    if (!Nprimes) x = gen_1;
    else
    {
      x = isprincipalfact(bnf, Primes, vecsmall_to_col(x), NULL,
                          nf_FORCE | nf_GEN_IF_PRINCIPAL);
      x = coltoliftalg(nf, x);
    }
    /* x solution, up to sign */
    appendL(res, x);
  }
  return res;
}

GEN
bnfisintnorm(GEN bnf, GEN a)
{
  pari_sp av = avma;
  GEN unit = NULL, z = bnfisintnormabs(bnf, a);
  GEN nf = checknf(bnf), T = gel(nf,1);
  long sNx, i, j, N = degpol(T), l = lg(z), sa = signe(a);
  long norm_1 = 0; /* gcc -Wall */

  for (i = j = 1; i<l; i++)
  {
    GEN x = gel(z,i);
    int xpol = (typ(x) == t_POL);

    if (xpol) sNx = signe(ZX_resultant(T, Q_primpart(x)));
    else      sNx = gsigne(x) < 0 && odd(N) ? -1 : 1;
    if (sNx != sa)
    {
      if (! unit) norm_1 = get_unit_1(bnf, &unit);
      if (!norm_1)
      {
        if (DEBUGLEVEL > 2) fprintferr("%Z eliminated because of sign\n",x);
        continue;
      }
      if (xpol) x = (unit == gen_m1)? RgX_neg(x): RgXQ_mul(unit,x,T);
      else      x = (unit == gen_m1)? gneg(x): RgX_Rg_mul(unit,x);
    }
    gel(z,j++) = x;
  }
  setlg(z, j);
  return gerepilecopy(av, z);
}
