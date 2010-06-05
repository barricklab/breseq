/* $Id: subfield.c 10283 2008-06-09 11:31:42Z kb $

Copyright (C) 2000-2004  The PARI group.

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
/*               SUBFIELDS OF A NUMBER FIELD                       */
/*   J. Klueners and M. Pohst, J. Symb. Comp. (1996), vol. 11      */
/*                                                                 */
/*******************************************************************/
#include "pari.h"
#include "paripriv.h"

typedef struct _poldata {
  GEN pol;
  GEN dis; /* |disc(pol)| */
  GEN roo; /* roots(pol) */
  GEN den; /* multiple of index(pol) */
} poldata;
typedef struct _primedata {
  GEN p;  /* prime */
  GEN pol; /* pol mod p, squarefree */
  GEN ff; /* factorization of pol mod p */
  GEN Z; /* cycle structure of the above [ Frobenius orbits ] */
  long lcm; /* lcm of the above */
  GEN T;  /* ffinit(p, lcm) */

  GEN fk;      /* factorization of pol over F_(p^lcm) */
  GEN firstroot; /* *[i] = index of first root of fk[i] */
  GEN interp;    /* *[i] = interpolation polynomial for fk[i]
                  * [= 1 on the first root firstroot[i], 0 on the others] */
  GEN bezoutC; /* Bezout coefficients associated to the ff[i] */
  GEN Trk;     /* used to compute traces (cf poltrace) */
} primedata;
typedef struct _blockdata {
  poldata *PD; /* data depending from pol */
  primedata *S;/* data depending from pol, p */
  GEN DATA; /* data depending from pol, p, degree, # translations [updated] */
  long N; /* deg(PD.pol) */
  long d; /* subfield degree */
  long size;/* block degree = N/d */
} blockdata;

static GEN print_block_system(blockdata *B, GEN Y, GEN BS);
static GEN test_block(blockdata *B, GEN L, GEN D);

/* return P(X + c) using destructive Horner, optimize for c = 1,-1 */
GEN
translate_pol(GEN P, GEN c)
{
  pari_sp av = avma, lim;
  GEN Q, *R;
  long i, k, n;

  if (!signe(P) || gcmp0(c)) return gcopy(P);
  Q = shallowcopy(P);
  R = (GEN*)(Q+2); n = degpol(P);
  lim = stack_lim(av, 2);
  if (gcmp1(c))
  {
    for (i=1; i<=n; i++)
    {
      for (k=n-i; k<n; k++) R[k] = gadd(R[k], R[k+1]);
      if (low_stack(lim, stack_lim(av,2)))
      {
        if(DEBUGMEM>1) pari_warn(warnmem,"TR_POL(1), i = %ld/%ld", i,n);
        Q = gerepilecopy(av, Q); R = (GEN*)Q+2;
      }
    }
  }
  else if (gcmp_1(c))
  {
    for (i=1; i<=n; i++)
    {
      for (k=n-i; k<n; k++) R[k] = gsub(R[k], R[k+1]);
      if (low_stack(lim, stack_lim(av,2)))
      {
        if(DEBUGMEM>1) pari_warn(warnmem,"TR_POL(-1), i = %ld/%ld", i,n);
        Q = gerepilecopy(av, Q); R = (GEN*)Q+2;
      }
    }
  }
  else
  {
    for (i=1; i<=n; i++)
    {
      for (k=n-i; k<n; k++) R[k] = gadd(R[k], gmul(c, R[k+1]));
      if (low_stack(lim, stack_lim(av,2)))
      {
        if(DEBUGMEM>1) pari_warn(warnmem,"TR_POL, i = %ld/%ld", i,n);
        Q = gerepilecopy(av, Q); R = (GEN*)Q+2;
      }
    }
  }
  return gerepilecopy(av, Q);
}

/* COMBINATORIAL PART: generate potential block systems */

#define BIL 32 /* for 64bit machines also */
/* Computation of potential block systems of given size d associated to a
 * rational prime p: give a row vector of row vectors containing the
 * potential block systems of imprimitivity; a potential block system is a
 * vector of row vectors (enumeration of the roots). */
static GEN
calc_block(blockdata *B, GEN Z, GEN Y, GEN SB)
{
  long r = lg(Z), lK, i, j, t, tp, T, u, nn, lnon, lY;
  GEN K, n, non, pn, pnon, e, Yp, Zp, Zpp;
  pari_sp av0 = avma;

  if (DEBUGLEVEL>3)
  {
    fprintferr("lg(Z) = %ld, lg(Y) = %ld\n", r,lg(Y));
    if (DEBUGLEVEL > 5)
    {
      fprintferr("Z = %Z\n",Z);
      fprintferr("Y = %Z\n",Y);
    }
  }
  lnon = min(BIL, r);
  e    = new_chunk(BIL);
  n    = new_chunk(r);
  non  = new_chunk(lnon);
  pnon = new_chunk(lnon);
  pn   = new_chunk(lnon);

  Zp   = cgetg(lnon,t_VEC);
  Zpp  = cgetg(lnon,t_VEC); nn = 0;
  for (i=1; i<r; i++) { n[i] = lg(Z[i])-1; nn += n[i]; }
  lY = lg(Y); Yp = cgetg(lY+1,t_VEC);
  for (j=1; j<lY; j++) Yp[j] = Y[j];

  {
    pari_sp av = avma;
    long k = nn / B->size;
    for (j = 1; j < r; j++) 
      if (n[j] % k) break;
    if (j == r)
    {
      gel(Yp,lY) = Z;
      SB = print_block_system(B, Yp, SB);
      avma = av;
    }
  }
  gel(Yp,lY) = Zp;

  K = divisors(utoipos(n[1])); lK = lg(K);
  for (i=1; i<lK; i++)
  {
    long ngcd = n[1], k = itos(gel(K,i)), dk = B->size*k, lpn = 0;
    for (j=2; j<r; j++)
      if (n[j]%k == 0)
      {
        if (++lpn >= BIL) pari_err(talker,"overflow in calc_block");
        pn[lpn] = n[j]; pnon[lpn] = j;
        ngcd = cgcd(ngcd, n[j]);
      }
    if (dk % ngcd) continue;
    T = 1<<lpn;
    if (lpn == r-2)
    {
      T--; /* done already above --> print_block_system */
      if (!T) continue;
    }

    if (dk == n[1])
    { /* empty subset, t = 0. Split out for clarity */
      Zp[1] = Z[1]; setlg(Zp, 2);
      for (u=1,j=2; j<r; j++) Zpp[u++] = Z[j];
      setlg(Zpp, u);
      SB = calc_block(B, Zpp, Yp, SB);
    }

    for (t = 1; t < T; t++)
    { /* loop through all non-empty subsets of [1..lpn] */
      for (nn=n[1],tp=t, u=1; u<=lpn; u++,tp>>=1)
      {
        if (tp&1) { nn += pn[u]; e[u] = 1; } else e[u] = 0;
      }
      if (dk != nn) continue;

      for (j=1; j<r; j++) non[j]=0;
      Zp[1] = Z[1];
      for (u=2,j=1; j<=lpn; j++)
        if (e[j]) { Zp[u] = Z[pnon[j]]; non[pnon[j]] = 1; u++; }
      setlg(Zp, u);
      for (u=1,j=2; j<r; j++)
        if (!non[j]) Zpp[u++] = Z[j];
      setlg(Zpp, u);
      SB = calc_block(B, Zpp, Yp, SB);
    }
  }
  avma = av0; return SB;
}

/* product of permutations. Put the result in perm1. */
static void
perm_mul_i(GEN perm1, GEN perm2)
{
  long i, N = lg(perm1);
  pari_sp av = avma;
  GEN perm = new_chunk(N);
  for (i=1; i<N; i++) perm[i] = perm1[perm2[i]];
  for (i=1; i<N; i++) perm1[i]= perm[i];
  avma = av;
}

/* cy is a cycle; compute cy^l as a permutation */
static GEN
cycle_power_to_perm(GEN perm,GEN cy,long l)
{
  long lp,i,j,b, N = lg(perm), lcy = lg(cy)-1;

  lp = l % lcy;
  for (i=1; i<N; i++) perm[i] = i;
  if (lp)
  {
    pari_sp av = avma;
    GEN p1 = new_chunk(N);
    b = cy[1];
    for (i=1; i<lcy; i++) b = (perm[b] = cy[i+1]);
    perm[b] = cy[1];
    for (i=1; i<N; i++) p1[i] = perm[i];

    for (j=2; j<=lp; j++) perm_mul_i(perm,p1);
    avma = av;
  }
  return perm;
}

/* image du block system D par la permutation perm */
static GEN
im_block_by_perm(GEN D,GEN perm)
{
  long i,j,lb,lcy;
  GEN Dn,cy,p1;

  lb=lg(D); Dn=cgetg(lb,t_VEC);
  for (i=1; i<lb; i++)
  {
    cy=gel(D,i); lcy=lg(cy);
    gel(Dn,i) = cgetg(lcy,t_VECSMALL); p1=gel(Dn,i);
    for (j=1; j<lcy; j++) p1[j] = perm[cy[j]];
  }
  return Dn;
}

static void
append(GEN D, GEN a)
{
  long i,l = lg(D), m = lg(a);
  GEN x = D + (l-1);
  for (i=1; i<m; i++) x[i] = a[i];
  setlg(D, l+m-1);
}

static GEN
print_block_system(blockdata *B, GEN Y, GEN SB)
{
  long i, j, l, ll, lp, u, v, ns, r = lg(Y), N = B->N;
  long *k, *n, **e, *t;
  GEN D, De, Z, cyperm, perm, VOID = cgetg(1, t_VECSMALL);

  if (DEBUGLEVEL>5) fprintferr("Y = %Z\n",Y);
  n = new_chunk(N+1);
  D = cget1(N+1, t_VEC);
  t = new_chunk(r+1);
  k = new_chunk(r+1);
  Z = cgetg(r+1, t_VEC);
  for (ns=0,i=1; i<r; i++)
  {
    GEN Yi = gel(Y,i);
    long ki = 0, si = lg(Yi)-1;

    for (j=1; j<=si; j++) { n[j] = lg(Yi[j])-1; ki += n[j]; }
    ki /= B->size;
    De = cgetg(ki+1,t_VEC);
    for (j=1; j<=ki; j++) gel(De,j) = VOID;
    for (j=1; j<=si; j++)
    {
      GEN cy = gel(Yi,j);
      for (l=1,lp=0; l<=n[j]; l++)
      {
        lp++; if (lp > ki) lp = 1;
        gel(De,lp) = vecsmall_append(gel(De,lp), cy[l]);
      }
    }
    append(D, De);
    if (si>1 && ki>1)
    {
      GEN p1 = cgetg(si,t_VEC);
      for (j=2; j<=si; j++) p1[j-1] = Yi[j];
      ns++;
      t[ns] = si-1;
      k[ns] = ki-1;
      gel(Z,ns) = p1;
    }
  }
  if (DEBUGLEVEL>2) fprintferr("\nns = %ld\n",ns);
  if (!ns) return test_block(B, SB, D);

  setlg(Z, ns+1);
  e = (long**)new_chunk(ns+1);
  for (i=1; i<=ns; i++)
  {
    e[i] = new_chunk(t[i]+1);
    for (j=1; j<=t[i]; j++) e[i][j] = 0;
  }
  cyperm= cgetg(N+1,t_VECSMALL);
  perm  = cgetg(N+1,t_VECSMALL); i = ns;
  do
  {
    pari_sp av = avma;
    for (u=1; u<=N; u++) perm[u] = u;
    for (u=1; u<=ns; u++)
      for (v=1; v<=t[u]; v++)
	perm_mul_i(perm, cycle_power_to_perm(cyperm, gmael(Z,u,v), e[u][v]));
    SB = test_block(B, SB, im_block_by_perm(D,perm));
    avma = av;

    /* i = 1..ns, j = 1..t[i], e[i][j] loop through 0..k[i].
     * TODO: flatten to 1-dimensional loop */
    if (++e[ns][t[ns]] > k[ns])
    {
      j = t[ns]-1;
      while (j>=1 && e[ns][j] == k[ns]) j--;
      if (j >= 1) { e[ns][j]++; for (l=j+1; l<=t[ns]; l++) e[ns][l] = 0; }
      else
      {
	i = ns-1;
	while (i>=1)
	{
	  j = t[i];
	  while (j>=1 && e[i][j] == k[i]) j--;
	  if (j<1) i--;
          else
	  {
	    e[i][j]++;
	    for (l=j+1; l<=t[i]; l++) e[i][l] = 0;
	    for (ll=i+1; ll<=ns; ll++)
              for (l=1; l<=t[ll]; l++) e[ll][l] = 0;
            break;
	  }
	}
      }
    }
  }
  while (i > 0);
  return SB;
}

/* ALGEBRAIC PART: test potential block systems */

static GEN
polsimplify(GEN x)
{
  long i,lx = lg(x);
  for (i=2; i<lx; i++)
    if (typ(x[i]) == t_POL) gel(x,i) = constant_term(gel(x,i));
  return x;
}

/* return 0 if |g[i]| > M[i] for some i; 1 otherwise */
static long
ok_coeffs(GEN g,GEN M)
{
  long i, lg = lg(g)-1; /* g is monic, and cst term is ok */
  for (i=3; i<lg; i++)
    if (absi_cmp(gel(g,i), gel(M,i)) > 0) return 0;
  return 1;
}

/* assume x in Fq, return Tr_{Fq/Fp}(x) as a t_INT */
static GEN
trace(GEN x, GEN Trq, GEN p)
{
  long i, l;
  GEN s;
  if (typ(x) == t_INT) return modii(mulii(x, gel(Trq,1)), p);
  l = lg(x)-1; if (l == 1) return gen_0;
  x++; s = mulii(gel(x,1), gel(Trq,1));
  for (i=2; i<l; i++)
    s = addii(s, mulii(gel(x,i), gel(Trq,i)));
  return modii(s, p);
}

/* assume x in Fq[X], return Tr_{Fq[X]/Fp[X]}(x), varn(X) = 0 */
static GEN
poltrace(GEN x, GEN Trq, GEN p)
{
  long i,l;
  GEN y;
  if (typ(x) == t_INT || varn(x) != 0) return trace(x, Trq, p);
  l = lg(x); y = cgetg(l,t_POL); y[1]=x[1];
  for (i=2; i<l; i++) gel(y,i) = trace(gel(x,i),Trq,p);
  return y;
}

/* Find h in Fp[X] such that h(a[i]) = listdelta[i] for all modular factors
 * ff[i], where a[i] is a fixed root of ff[i] in Fq = Z[Y]/(p,T) [namely the
 * first one in FpX_factorff_irred output]. Let f = ff[i], A the given root,
 * then h mod f is Tr_Fq/Fp ( h(A) f(X)/(X-A)f'(A) ), most of the expression
 * being precomputed. The complete h is recovered via chinese remaindering */
static GEN
chinese_retrieve_pol(GEN DATA, primedata *S, GEN listdelta)
{
  GEN interp, bezoutC, h, pol = FpX_red(gel(DATA,1), S->p);
  long i, l;
  interp = gel(DATA,9);
  bezoutC= gel(DATA,6);

  h = NULL; l = lg(interp);
  for (i=1; i<l; i++)
  { /* h(firstroot[i]) = listdelta[i] */
    GEN t = FqX_Fq_mul(gel(interp,i), gel(listdelta,i), S->T,S->p);
    t = poltrace(t, (GEN)S->Trk[i], S->p);
    t = gmul(t, gel(bezoutC,i));
    h = h? gadd(h,t): t;
  }
  return FpX_rem(FpX_red(h, S->p), pol, S->p);
}

/* g in Z[X] potentially defines a subfield of Q[X]/f. It is a subfield iff A
 * (cf subfield) was a block system; then there
 * exists h in Q[X] such that f | g o h. listdelta determines h s.t f | g o h
 * in Fp[X] (cf chinese_retrieve_pol). Try to lift it; den is a
 * multiplicative bound for denominator of lift. */
static GEN
embedding(GEN g, GEN DATA, primedata *S, GEN den, GEN listdelta)
{
  GEN TR, w0_Q, w0, w1_Q, w1, wpow, h0, gp, T, q2, q, maxp, a, p = S->p;
  long rt;
  pari_sp av;

  T   = gel(DATA,1); rt = brent_kung_optpow(degpol(T), 2);
  maxp= gel(DATA,7);
  gp = derivpol(g); av = avma;
  w0 = chinese_retrieve_pol(DATA, S, listdelta);
  w0_Q = centermod(gmul(w0,den), p);
  h0 = FpXQ_inv(FpX_FpXQ_compo(gp,w0, T,p), T,p); /* = 1/g'(w0) mod (T,p) */
  wpow = NULL; q = sqri(p);
  for(;;)
  {/* Given g,w0,h0 in Z[x], s.t. h0.g'(w0) = 1 and g(w0) = 0 mod (T,p), find
    * [w1,h1] satisfying the same conditions mod p^2, [w1,h1] = [w0,h0] (mod p)
    * (cf. Dixon: J. Austral. Math. Soc., Series A, vol.49, 1990, p.445) */
    if (DEBUGLEVEL>1)
      fprintferr("lifting embedding mod p^k = %Z^%ld\n",p, Z_pval(q,p));

    /* w1 := w0 - h0 g(w0) mod (T,q) */
    if (wpow)
      a = FpX_FpXQV_compo(g,wpow, T,q);
    else
      a = FpX_FpXQ_compo(g,w0, T,q); /* first time */
    /* now, a = 0 (p) */
    a = gmul(gneg(h0), gdivexact(a, p));
    w1 = gadd(w0, gmul(p, FpX_rem(a, T,p)));

    w1_Q = centermod(gmul(w1, remii(den,q)), q);
    if (gequal(w1_Q, w0_Q) || cmpii(q,maxp) > 0)
    {
      GEN G = gcmp1(den)? g: RgX_rescale(g,den);
      if (gcmp0(RgX_RgXQ_compo(G, w1_Q, T))) break;
    }
    if (cmpii(q, maxp) > 0)
    {
      if (DEBUGLEVEL) fprintferr("coeff too big for embedding\n");
      return NULL;
    }
    gerepileall(av, 5, &w1,&h0,&w1_Q,&q,&p);
    q2 = sqri(q);
    wpow = FpXQ_powers(w1, rt, T, q2);
    /* h0 := h0 * (2 - h0 g'(w1)) mod (T,q)
     *     = h0 + h0 * (1 - h0 g'(w1)) */
    a = gmul(gneg(h0), FpX_FpXQV_compo(gp, FpXV_red(wpow,q),T,q));
    a = ZX_Z_add(FpX_rem(a, T,q), gen_1); /* 1 - h0 g'(w1) = 0 (p) */
    a = gmul(h0, gdivexact(a, p));
    h0 = gadd(h0, gmul(p, FpX_rem(a, T,p)));
    w0 = w1; w0_Q = w1_Q; p = q; q = q2;
  }
  TR = gel(DATA,5);
  if (!gcmp0(TR)) w1_Q = translate_pol(w1_Q, TR);
  return gdiv(w1_Q,den);
}

/* return U list of polynomials s.t U[i] = 1 mod fk[i] and 0 mod fk[j] for all
 * other j */
static GEN
get_bezout(GEN pol, GEN fk, GEN p)
{
  long i, l = lg(fk);
  GEN A, B, d, u, v, U = cgetg(l, t_VEC);
  for (i=1; i<l; i++)
  {
    A = gel(fk,i);
    B = FpX_div(pol, A, p);
    d = FpX_extgcd(A,B,p, &u, &v);
    if (degpol(d) > 0) pari_err(talker, "relatively prime polynomials expected");
    d = constant_term(d);
    if (!gcmp1(d)) v = FpX_Fp_mul(v, Fp_inv(d, p), p);
    gel(U,i) = FpX_mul(B,v, p);
  }
  return U;
}

static GEN
init_traces(GEN ff, GEN T, GEN p)
{
  long N = degpol(T),i,j,k, r = lg(ff);
  GEN Frob = FpXQ_matrix_pow(FpXQ_pow(pol_x[varn(T)],p, T,p), N,N, T,p);
  GEN y,p1,p2,Trk,pow,pow1;

  k = degpol(ff[r-1]); /* largest degree in modular factorization */
  pow = cgetg(k+1, t_VEC);
  gel(pow,1) = gen_0; /* dummy */
  gel(pow,2) = Frob;
  pow1= cgetg(k+1, t_VEC); /* 1st line */
  for (i=3; i<=k; i++)
    gel(pow,i) = FpM_mul(gel(pow,i-1), Frob, p);
  gel(pow1,1) = gen_0; /* dummy */
  for (i=2; i<=k; i++)
  {
    p1 = cgetg(N+1, t_VEC);
    gel(pow1,i) = p1; p2 = gel(pow,i);
    for (j=1; j<=N; j++) gel(p1,j) = gcoeff(p2,1,j);
  }
  
  /* Trk[i] = line 1 of x -> x + x^p + ... + x^{p^(i-1)} */
  Trk = pow; /* re-use (destroy) pow */
  gel(Trk,1) = vec_ei(N,1);
  for (i=2; i<=k; i++)
    gel(Trk,i) = gadd(gel(Trk,i-1), gel(pow1,i));
  y = cgetg(r, t_VEC);
  for (i=1; i<r; i++) y[i] = Trk[degpol(ff[i])];
  return y;
}

/* return C in Z[X]/(p,T), C[ D[1] ] = 1, C[ D[i] ] = 0 otherwise. H is the
 * list of degree 1 polynomials X - D[i]  (come directly from factorization) */
static GEN
interpol(GEN H, GEN T, GEN p)
{
  long i, m = lg(H);
  GEN X = pol_x[0],d,p1,p2,a;

  p1=pol_1[0]; p2=gen_1; a = gneg(constant_term(gel(H,1))); /* = D[1] */
  for (i=2; i<m; i++)
  {
    d = constant_term(gel(H,i)); /* -D[i] */
    p1 = FpXQX_mul(p1,gadd(X,d), T,p);
    p2 = Fq_mul(p2, gadd(a,d), T,p);
  }
  return FqX_Fq_mul(p1,Fq_inv(p2, T,p), T,p);
}

static void
init_primedata(primedata *S)
{
  long i, j, l, lff = lg(S->ff), v = fetch_var(), N = degpol(S->pol);
  GEN T, p = S->p;

  if (S->lcm == degpol(S->ff[lff-1]))
  {
    T = shallowcopy((GEN)S->ff[lff-1]);
    setvarn(T, v);
  }
  else
    T = init_Fq(p, S->lcm, v);
  name_var(v,"y");
  S->T = T;
  S->firstroot = cgetg(lff, t_VECSMALL);
  S->interp = cgetg(lff, t_VEC);
  S->fk = cgetg(N+1, t_VEC);
  for (l=1,j=1; j<lff; j++)
  { /* compute roots and fix ordering (Frobenius cycles) */
    GEN F = FpX_factorff_irred(gel(S->ff,j), T, p);
    gel(S->interp,j) = interpol(F,T,p);
    S->firstroot[j] = l;
    for (i=1; i<lg(F); i++,l++) S->fk[l] = F[i];
  }
  S->Trk     = init_traces(S->ff, T,p);
  S->bezoutC = get_bezout(S->pol, S->ff, p);
}

static void
choose_prime(primedata *S, GEN pol, GEN dpol)
{
  byteptr di = diffptr + 1;
  long i, j, k, r, lcm, oldlcm, pp, N = degpol(pol), minp = N*N / 4;
  GEN Z, p, ff, oldff, n, oldn;
  pari_sp av;

  if (DEBUGLEVEL) (void)timer2();
  p = utoipos(2);
  while (p[2] <= minp) NEXT_PRIME_VIADIFF(p[2], di);
  oldlcm = 0;
  oldff = oldn = NULL; pp = 0; /* gcc -Wall */
  av = avma;
  for(k = 1; k < 11 || !oldlcm; k++,avma = av)
  {
    do NEXT_PRIME_VIADIFF(p[2], di); while (!smodis(dpol, p[2]));
    if (k > 5 * N) pari_err(talker,"sorry, too many block systems in nfsubfields");
    ff = (GEN)FpX_factor(pol, p)[1];
    r = lg(ff)-1;
    if (r == N || r >= BIL) continue;

    n = cgetg(r+1, t_VECSMALL); lcm = n[1] = degpol(ff[1]);
    for (j=2; j<=r; j++) { n[j] = degpol(ff[j]); lcm = clcm(lcm, n[j]); }
    if (lcm <= oldlcm) continue; /* false when oldlcm = 0 */

    if (DEBUGLEVEL) fprintferr("p = %ld,\tlcm = %ld,\torbits: %Z\n",p[2],lcm,n);
    pp = p[2];
    oldn = n;
    oldff = ff;
    oldlcm = lcm; if (r == 1) break;
    av = avma;
  }
  if (DEBUGLEVEL) fprintferr("Chosen prime: p = %ld\n", pp);
  S->ff = oldff;
  S->lcm= oldlcm;
  S->p  = utoipos(pp);
  S->pol = FpX_red(pol, S->p); init_primedata(S);

  n = oldn; r = lg(n); Z = cgetg(r,t_VEC);
  for (k=0,i=1; i<r; i++)
  {
    GEN t = cgetg(n[i]+1, t_VECSMALL); gel(Z,i) = t;
    for (j=1; j<=n[i]; j++) t[j] = ++k;
  }
  S->Z = Z;
}

static GEN
bound_for_coeff(long m, GEN rr, GEN *maxroot)
{
  long i,r1, lrr=lg(rr);
  GEN p1,b1,b2,B,M, C = matpascal(m-1);

  for (r1=1; r1 < lrr; r1++)
    if (typ(rr[r1]) != t_REAL) break;
  r1--;

  rr = gabs(rr,0); *maxroot = vecmax(rr);
  for (i=1; i<lrr; i++)
    if (gcmp(gel(rr,i), gen_1) < 0) gel(rr,i) = gen_1;
  for (b1=gen_1,i=1; i<=r1; i++) b1 = gmul(b1, gel(rr,i));
  for (b2=gen_1    ; i<lrr; i++) b2 = gmul(b2, gel(rr,i));
  B = gmul(b1, gsqr(b2)); /* Mahler measure */
  M = cgetg(m+2, t_VEC); gel(M,1) = gel(M,2) = gen_0; /* unused */
  for (i=1; i<m; i++)
  {
    p1 = gadd(gmul(gcoeff(C, m, i+1), B),/* binom(m-1, i)   */
              gcoeff(C, m, i));          /* binom(m-1, i-1) */
    gel(M,i+2) = ceil_safe(p1);
  }
  return M;
}

/* d = requested degree for subfield. Return DATA, valid for given pol, S and d
 * If DATA != NULL, translate pol [ --> pol(X+1) ] and update DATA
 * 1: polynomial pol
 * 2: p^e (for Hensel lifts) such that p^e > max(M),
 * 3: Hensel lift to precision p^e of DATA[4]
 * 4: roots of pol in F_(p^S->lcm),
 * 5: number of polynomial changes (translations)
 * 6: Bezout coefficients associated to the S->ff[i]
 * 7: Hadamard bound for coefficients of h(x) such that g o h = 0 mod pol.
 * 8: bound M for polynomials defining subfields x PD->den
 * 9: *[i] = interpolation polynomial for S->ff[i] [= 1 on the first root
      S->firstroot[i], 0 on the others] */
static void
compute_data(blockdata *B)
{
  GEN ffL, roo, pe, p1, p2, fk, fhk, MM, maxroot, pol;
  primedata *S = B->S;
  GEN p = S->p, T = S->T, ff = S->ff, DATA = B->DATA;
  long i, j, l, e, N, lff = lg(ff);

  if (DEBUGLEVEL>1) fprintferr("Entering compute_data()\n\n");
  pol = B->PD->pol; N = degpol(pol);
  roo = B->PD->roo;
  if (DATA) /* update (translate) an existing DATA */
  {
    GEN Xm1 = gsub(pol_x[varn(pol)], gen_1);
    GEN TR = addis(gel(DATA,5), 1);
    GEN mTR = negi(TR), interp, bezoutC;

    gel(DATA,5) = TR;
    pol = translate_pol(gel(DATA,1), gen_m1);
    l = lg(roo); p1 = cgetg(l, t_VEC);
    for (i=1; i<l; i++) gel(p1,i) = gadd(TR, gel(roo,i));
    roo = p1;

    fk = gel(DATA,4); l = lg(fk);
    for (i=1; i<l; i++) gel(fk,i) = gsub(Xm1, gel(fk,i));

    bezoutC = gel(DATA,6); l = lg(bezoutC);
    interp  = gel(DATA,9);
    for (i=1; i<l; i++)
    {
      if (degpol(interp[i]) > 0) /* do not turn pol_1[0] into gen_1 */
      {
        p1 = translate_pol(gel(interp,i), gen_m1);
        gel(interp,i) = FpXX_red(p1, p);
      }
      if (degpol(bezoutC[i]) > 0)
      {
        p1 = translate_pol(gel(bezoutC,i), gen_m1);
        gel(bezoutC,i) = FpXX_red(p1, p);
      }
    }
    ff = cgetg(lff, t_VEC); /* copy, don't overwrite! */
    for (i=1; i<lff; i++)
      gel(ff,i) = FpX_red(translate_pol((GEN)S->ff[i], mTR), p);
  }
  else
  {
    DATA = cgetg(10,t_VEC);
    fk = S->fk;
    gel(DATA,5) = gen_0;
    gel(DATA,6) = shallowcopy(S->bezoutC);
    gel(DATA,9) = shallowcopy(S->interp);
  }
  gel(DATA,1) = pol;
  MM = gmul2n(bound_for_coeff(B->d, roo, &maxroot), 1);
  gel(DATA,8) = MM;
  e = logint(shifti(vecmax(MM),20), p, &pe); /* overlift 2^20 [for d-1 test] */
  gel(DATA,2) = pe;
  gel(DATA,4) = roots_from_deg1(fk);

  /* compute fhk = hensel_lift_fact(pol,fk,T,p,pe,e) in 2 steps
   * 1) lift in Zp to precision p^e */
  ffL = hensel_lift_fact(pol, ff, NULL, p, pe, e);
  fhk = NULL;
  for (l=i=1; i<lff; i++)
  { /* 2) lift factorization of ff[i] in Qp[X] / T */
    GEN F, L = gel(ffL,i);
    long di = degpol(L);
    F = cgetg(di+1, t_VEC);
    for (j=1; j<=di; j++) F[j] = fk[l++];
    L = hensel_lift_fact(L, F, T, p, pe, e);
    fhk = fhk? shallowconcat(fhk, L): L;
  }
  gel(DATA,3) = roots_from_deg1(fhk);

  p1 = mulsr(N, gsqrt(gpowgs(utoipos(N-1),N-1),DEFAULTPREC));
  p2 = gpowgs(maxroot, B->size + N*(N-1)/2);
  p1 = gdiv(gmul(p1,p2), gsqrt(B->PD->dis,DEFAULTPREC));
  gel(DATA,7) = mulii(shifti(ceil_safe(p1), 1), B->PD->den);

  if (DEBUGLEVEL>1) {
    fprintferr("f = %Z\n",DATA[1]);
    fprintferr("p = %Z, lift to p^%ld\n", p, e);
    fprintferr("2 * Hadamard bound * ind = %Z\n",DATA[7]);
    fprintferr("2 * M = %Z\n",DATA[8]);
  }
  if (B->DATA) {
    DATA = gclone(DATA);
    if (isclone(B->DATA)) gunclone(B->DATA);
  }
  B->DATA = DATA;
}

/* g = polynomial, h = embedding. Return [[g,h]] */
static GEN
_subfield(GEN g, GEN h) { return mkvec(mkvec2(g,h)); }

/* Return a subfield, gen_0 [ change p ] or NULL [ not a subfield ] */
static GEN
subfield(GEN A, blockdata *B)
{
  long N, i, j, d, lf, m = lg(A)-1;
  GEN M, pe, pol, fhk, g, e, d_1_term, delta, listdelta, whichdelta;
  GEN T = B->S->T, p = B->S->p, firstroot = B->S->firstroot;

  pol= (GEN)B->DATA[1]; N = degpol(pol); d = N/m; /* m | N */
  pe = (GEN)B->DATA[2];
  fhk= (GEN)B->DATA[3];
  M  = (GEN)B->DATA[8];

  delta = cgetg(m+1,t_VEC);
  whichdelta = cgetg(N+1, t_VECSMALL);
  d_1_term = gen_0;
  for (i=1; i<=m; i++)
  {
    GEN Ai = gel(A,i), p1 = (GEN)fhk[Ai[1]];
    for (j=2; j<=d; j++)
      p1 = Fq_mul(p1, (GEN)fhk[Ai[j]], T, pe);
    gel(delta,i) = p1;
    if (DEBUGLEVEL>2) fprintferr("delta[%ld] = %Z\n",i,p1);
    /* g = prod (X - delta[i])
     * if g o h = 0 (pol), we'll have h(Ai[j]) = delta[i] for all j */
    /* fk[k] belongs to block number whichdelta[k] */
    for (j=1; j<=d; j++) whichdelta[Ai[j]] = i;
    if (typ(p1) == t_POL) p1 = constant_term(p1);
    d_1_term = addii(d_1_term, p1);
  }
  d_1_term = centermod(d_1_term, pe); /* Tr(g) */
  if (absi_cmp(d_1_term, gel(M,3)) > 0) {
    if (DEBUGLEVEL>1) fprintferr("d-1 test failed\n");
    return NULL;
  }
  g = FqV_roots_to_pol(delta, T, pe, 0);
  g = centermod(polsimplify(g), pe); /* assume g in Z[X] */
  if (DEBUGLEVEL>2) fprintferr("pol. found = %Z\n",g);
  if (!ok_coeffs(g,M)) {
    if (DEBUGLEVEL>1) fprintferr("coeff too big for pol g(x)\n");
    return NULL;
  }
  if (!FpX_is_squarefree(g, p)) {
    if (DEBUGLEVEL>1) fprintferr("changing f(x): p divides disc(g)\n");
    compute_data(B);
    return subfield(A, B);
  }

  lf = lg(firstroot); listdelta = cgetg(lf, t_VEC);
  for (i=1; i<lf; i++) listdelta[i] = delta[whichdelta[firstroot[i]]];
  if (DEBUGLEVEL) fprintferr("candidate = %Z\n", g);
  e = embedding(g, B->DATA, B->S, B->PD->den, listdelta);
  if (!e) return NULL;
  if (DEBUGLEVEL) fprintferr("embedding = %Z\n", e);
  return _subfield(g, e);
}

/* L list of current subfields, test whether potential block D is a block,
 * if so, append corresponding subfield */
static GEN
test_block(blockdata *B, GEN L, GEN D)
{
  pari_sp av = avma;
  GEN sub = subfield(D, B);
  if (sub) {
    GEN old = L;
    L = gclone( L? shallowconcat(L, sub): sub );
    if (old) gunclone(old);
  }
  avma = av; return L;
}

/* subfields of degree d */
static GEN
subfields_of_given_degree(blockdata *B)
{
  pari_sp av = avma;
  GEN L;

  if (DEBUGLEVEL) fprintferr("\n* Look for subfields of degree %ld\n\n", B->d);
  B->DATA = NULL; compute_data(B);
  L = calc_block(B, B->S->Z, cgetg(1,t_VEC), NULL);
  if (DEBUGLEVEL) fprintferr("\nSubfields of degree %ld: %Z\n", B->d, L);
  if (isclone(B->DATA)) gunclone(B->DATA);
  avma = av; return L;
}

static GEN
fix_var(GEN x, long v)
{
  long i, l = lg(x);
  if (!v) return x;
  for (i=1; i<l; i++) { GEN t = gel(x,i); setvarn(t[1],v); setvarn(t[2],v); }
  return x;
}

static void
subfields_poldata(GEN T, poldata *PD)
{
  GEN  nf,L,dis;

  T = shallowcopy(get_nfpol(T, &nf));
  PD->pol = T; setvarn(T, 0);
  if (nf)
  {
    PD->den = Q_denom(gel(nf,7));
    PD->roo = gel(nf,6);
    PD->dis = mulii(absi(gel(nf,3)), sqri(gel(nf,4)));
  }
  else
  {
    PD->den = initgaloisborne(T,NULL,ZX_get_prec(T), &L,NULL,&dis);
    PD->roo = L;
    PD->dis = absi(dis);
  }
}

GEN
subfields(GEN nf, GEN d0)
{
  pari_sp av = avma;
  long N, v0, d = itos(d0);
  GEN LSB, pol, G;
  poldata PD;
  primedata S;
  blockdata B;

  pol = get_nfpol(nf, &nf); /* in order to treat trivial cases */
  v0 = varn(pol); N = degpol(pol);
  if (d == N) return gerepilecopy(av, _subfield(pol, pol_x[v0]));
  if (d == 1) return gerepilecopy(av, _subfield(pol_x[v0], pol));
  if (d < 1 || d > N || N % d) return cgetg(1,t_VEC);

  /* much easier if nf is Galois (WSS) */
  G = galoisconj4(nf? nf: pol, NULL, 1);
  if (typ(G) != t_INT)
  { /* Bingo */
    GEN L = galoissubgroups(G), F;
    long k,i, l = lg(L), o = N/d;
    F = cgetg(l, t_VEC);
    k = 1;
    for (i=1; i<l; i++)
    {
      GEN H = gel(L,i);
      if (group_order(H) == o)
        gel(F,k++) = lift_intern(galoisfixedfield(G, gel(H,1), 0, v0));
    }
    setlg(F, k);
    return gerepilecopy(av, F);
  }

  subfields_poldata(nf? nf: pol, &PD);

  B.PD = &PD;
  B.S  = &S;
  B.N  = N;
  B.d  = d;
  B.size = N/d;

  choose_prime(&S, PD.pol, PD.dis);
  LSB = subfields_of_given_degree(&B);
  (void)delete_var(); /* from choose_prime */
  avma = av;
  if (!LSB) return cgetg(1, t_VEC);
  G = gcopy(LSB); gunclone(LSB);
  return fix_var(G, v0);
}

static GEN
subfieldsall(GEN nf)
{
  pari_sp av = avma;
  long N, ld, i, v0;
  GEN G, pol, dg, LSB, NLSB;
  poldata PD;
  primedata S;
  blockdata B;

  /* much easier if nf is Galois (WSS) */
  G = galoisconj4(nf, NULL, 1);
  if (typ(G) != t_INT)
  {
    GEN L, S, p;
    long l;

    pol = get_nfpol(nf, &nf);
    L = lift_intern( galoissubfields(G, 0, varn(pol)) );
    l = lg(L);
    S = cgetg(l, t_VECSMALL);
    for (i=1; i<l; i++) S[i] = lg(gmael(L,i,1));
    p = vecsmall_indexsort(S);
    return gerepilecopy(av,  vecpermute(L, p));
  }

  subfields_poldata(nf, &PD);
  pol = PD.pol;

  v0 = varn(pol); N = degpol(pol);
  dg = divisors(utoipos(N)); ld = lg(dg)-1;
  if (DEBUGLEVEL) fprintferr("\n***** Entering subfields\n\npol = %Z\n",pol);

  LSB = _subfield(pol, pol_x[0]);
  if (ld > 2)
  {
    B.PD = &PD;
    B.S  = &S;
    B.N  = N;
    choose_prime(&S, PD.pol, PD.dis);
    for (i=2; i<ld; i++)
    {
      B.size  = itos(gel(dg,i));
      B.d = N / B.size;
      NLSB = subfields_of_given_degree(&B);
      if (NLSB) { LSB = concat(LSB, NLSB); gunclone(NLSB); }
    }
    (void)delete_var(); /* from choose_prime */
  }
  LSB = shallowconcat(LSB, _subfield(pol_x[0], pol));
  if (DEBUGLEVEL) fprintferr("\n***** Leaving subfields\n\n");
  return fix_var(gerepilecopy(av, LSB), v0);
}

GEN
subfields0(GEN nf,GEN d)
{
  return d? subfields(nf,d): subfieldsall(nf);
}
