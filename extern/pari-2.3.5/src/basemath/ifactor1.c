/* $Id: ifactor1.c 11457 2008-12-16 18:18:47Z bill $

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
/**                     INTEGER FACTORIZATION                      **/
/**                                                                **/
/********************************************************************/
#include "pari.h"
#include "paripriv.h"

int factor_add_primes = 0;

/*********************************************************************/
/**                                                                 **/
/**               PSEUDO PRIMALITY (MILLER-RABIN)                   **/
/**                                                                 **/
/*********************************************************************/
typedef struct {
  ulong n, sqrt1, sqrt2, t1, t;
  long r1;
} Fl_miller_t;

typedef struct {
  GEN n, sqrt1, sqrt2, t1, t;
  long r1;
} miller_t;

static void
init_miller(miller_t *S, GEN n)
{
  if (signe(n) < 0) n = absi(n);
  S->n = n;
  S->t = addsi(-1,n);
  S->r1 = vali(S->t);
  S->t1 = shifti(S->t, -S->r1);
  S->sqrt1 = cgeti(lg(n)); S->sqrt1[1] = evalsigne(0)|evallgefint(2);
  S->sqrt2 = cgeti(lg(n)); S->sqrt2[1] = evalsigne(0)|evallgefint(2);
}
static void
Fl_init_miller(Fl_miller_t *S, ulong n)
{
  S->n = n;
  S->t = n-1;
  S->r1 = vals(S->t);
  S->t1 = S->t >> S->r1;
  S->sqrt1 = 0;
  S->sqrt2 = 0;
}

/* c = sqrt(-1) seen in bad_for_base. End-matching: compare or remember
 * If ends do mismatch, then we have factored n, and this information
 * should somehow be made available to the factoring machinery. But so
 * exceedingly rare... besides we use BSPW now. */
static int
miller_ok(miller_t *S, GEN c)
{
  if (signe(S->sqrt1))
  { /* saw one earlier: compare */
    if (!equalii(c, S->sqrt1) && !equalii(c, S->sqrt2))
    { /* too many sqrt(-1)s mod n */
      if (DEBUGLEVEL) {
        GEN z = gcdii(addii(c, S->sqrt1), S->n);
        pari_warn(warner,"found factor\n\t%Z\ncurrently lost to the factoring machinery", z);
      }
      return 1;
    }
  } else { /* remember */
    affii(c, S->sqrt1);
    affii(subii(S->n, c), S->sqrt2);
  }
  return 0;
}
static int
Fl_miller_ok(Fl_miller_t *S, ulong c)
{
  if (S->sqrt1)
  { /* saw one earlier: compare */
    if (c != S->sqrt1 && c != S->sqrt2) return 1;
  } else { /* remember */
    S->sqrt1 = c;
    S->sqrt2 = S->n - c;
  }
  return 0;
}

/* is n strong pseudo-prime for base a ? `End matching' (check for square
 * roots of -1) added by GN */
static int
bad_for_base(miller_t *S, GEN a)
{
  long r, lim, av = avma;
  GEN c2, c = Fp_pow(a, S->t1, S->n);

  if (is_pm1(c) || equalii(S->t, c)) return 0;

  lim = stack_lim(av,1);
  /* go fishing for -1, not for 1 (saves one squaring) */
  for (r = S->r1 - 1; r; r--) /* r1 - 1 squarings */
  {
    c2 = c; c = remii(sqri(c), S->n);
    if (equalii(S->t, c)) return miller_ok(S, c2);
    if (low_stack(lim, stack_lim(av,1)))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"miller(rabin)");
      c = gerepileuptoint(av, c);
    }
  }
  return 1;
}
static int
Fl_bad_for_base(Fl_miller_t *S, ulong a)
{
  long r;
  ulong c2, c = Fl_pow(a, S->t1, S->n);

  if (c == 1 || c == S->t) return 0;

  /* go fishing for -1, not for 1 (saves one squaring) */
  for (r = S->r1 - 1; r; r--) /* r1 - 1 squarings */
  {
    c2 = c; c = Fl_sqr(c, S->n);
    if (c == S->t) return Fl_miller_ok(S, c2);
  }
  return 1;
}

/* Miller-Rabin test for k random bases */
long
millerrabin(GEN n, long k)
{
  pari_sp av2, av = avma;
  ulong r;
  long i;
  miller_t S;

  if (!signe(n)) return 0;
  /* If |n| <= 3, check if n = +- 1 */
  if (lgefint(n)==3 && (ulong)(n[2])<=3) return (n[2] != 1);

  if (!mod2(n)) return 0;
  init_miller(&S, n); av2 = avma;
  for (i=1; i<=k; i++)
  {
    do r = umodui((ulong)pari_rand31(), n); while (!r);
    if (DEBUGLEVEL > 4) fprintferr("Miller-Rabin: testing base %ld\n", r);
    if (bad_for_base(&S, utoipos(r))) { avma = av; return 0; }
    avma = av2;
  }
  avma = av; return 1;
}

/* As above for k bases taken in pr (i.e not random). We must have |n|>2 and
 * 1<=k<=11 (not checked) or k in {16,17} to select some special sets of bases.
 *
 * From Jaeschke, `On strong pseudoprimes to several bases', Math.Comp. 61
 * (1993), 915--926  (see also http://www.utm.edu/research/primes/prove2.html),
 * we have:
 *
 * k == 4  (bases 2,3,5,7)  detects all composites
 *    n <     118 670 087 467 == 172243 * 688969  with the single exception of
 *    n ==      3 215 031 751 == 151 * 751 * 28351,
 *
 * k == 5  (bases 2,3,5,7,11)  detects all composites
 *    n <   2 152 302 898 747 == 6763 * 10627 * 29947,
 *
 * k == 6  (bases 2,3,...,13)  detects all composites
 *    n <   3 474 749 660 383 == 1303 * 16927 * 157543,
 *
 * k == 7  (bases 2,3,...,17)  detects all composites
 *    n < 341 550 071 728 321 == 10670053 * 32010157,
 * Even this limiting value is caught by an end mismatch between bases 5 and 17
 *
 * Moreover, the four bases chosen at
 *
 * k == 16  (2,13,23,1662803)  detects all composites up
 * to at least 10^12, and the combination at
 *
 * k == 17  (31,73)  detects most odd composites without prime factors > 100
 * in the range  n < 2^36  (with less than 250 exceptions, indeed with fewer
 * than 1400 exceptions up to 2^42). --GN */
static int
Fl_miller(ulong n, long k)
{
  static ulong pr[] =
    { 0, 2,3,5,7,11,13,17,19,23,29, 31,73, 2,13,23,1662803UL, };
  ulong *p;
  ulong r;
  long i;
  Fl_miller_t S;

  if (!(n & 1)) return 0;
  if (k == 16)
  { /* use smaller (faster) bases if possible */
    p = (n < 3215031751UL)? pr: pr+13;
    k = 4;
  }
  else if (k == 17)
  {
    p = (n < 1373653UL)? pr: pr+11;
    k = 2;
  }
  else p = pr; /* 2,3,5,... */
  Fl_init_miller(&S, n);
  for (i=1; i<=k; i++)
  {
    r = p[i] % n; if (!r) break;
    if (Fl_bad_for_base(&S, r)) return 0;
  }
  return 1;
}

int
miller(GEN n, long k)
{
  pari_sp av2, av = avma;
  static ulong pr[] =
    { 0, 2,3,5,7,11,13,17,19,23,29, 31,73, 2,13,23,1662803UL, };
  ulong *p;
  long i;
  miller_t S;

  if (lgefint(n) == 3) return Fl_miller((ulong)n[2], k);

  if (!mod2(n)) return 0;
  if      (k == 16) { p = pr+13; k = 4; } /* 2,13,23,1662803 */
  else if (k == 17) { p = pr+11; k = 2; } /* 31,73 */
  else p = pr; /* 2,3,5,... */
  init_miller(&S, n); av2 = avma;
  for (i=1; i<=k; i++)
  {
    if (bad_for_base(&S, utoipos(p[i]))) { avma = av; return 0; }
    avma = av2;
  }
  avma = av; return 1;
}

/*********************************************************************/
/**                                                                 **/
/**                      PSEUDO PRIMALITY (LUCAS)                   **/
/**                                                                 **/
/*********************************************************************/
/* compute n-th term of Lucas sequence modulo N.
 * v_{k+2} = P v_{k+1} - v_k, v_0 = 2, v_1 = P.
 * Assume n > 0 */
static GEN
LucasMod(GEN n, ulong P, GEN N)
{
  pari_sp av = avma, lim = stack_lim(av, 1);
  GEN nd = int_MSW(n);
  long i, m = *nd, j = 1+bfffo((ulong)m);
  GEN v = utoipos(P), v1 = utoipos(P*P - 2);

  m <<= j; j = BITS_IN_LONG - j;
  for (i=lgefint(n)-2;;) /* cf. leftright_pow */
  {
    for (; j; m<<=1,j--)
    { /* v = v_k, v1 = v_{k+1} */
      if (m < 0)
      { /* set v = v_{2k+1}, v1 = v_{2k+2} */
        v = subis(mulii(v,v1), (long)P);
        v1= subis(sqri(v1), 2);
      }
      else
      {/* set v = v_{2k}, v1 = v_{2k+1} */
        v1= subis(mulii(v,v1), (long)P);
        v = subis(sqri(v), 2);
      }
      v = modii(v, N);
      v1= modii(v1,N);
      if (low_stack(lim,stack_lim(av,1)))
      {
        if(DEBUGMEM>1) pari_warn(warnmem,"LucasMod");
        gerepileall(av, 2, &v,&v1);
      }
    }
    if (--i == 0) return v;
    j = BITS_IN_LONG;
    nd=int_precW(nd);
    m = *nd;
  }
}
/* compute n-th term of Lucas sequence modulo N.
 * v_{k+2} = P v_{k+1} - v_k, v_0 = 2, v_1 = P.
 * Assume n > 0 */
static ulong
u_LucasMod(ulong n, ulong P, ulong N)
{
  long j = 1 + bfffo(n);
  ulong v = P, v1 = P*P - 2, mP = N - P, m2 = N - 2, m = n << j;

  j = BITS_IN_LONG - j;
  for (; j; m<<=1,j--)
  { /* v = v_k, v1 = v_{k+1} */
    if (((long)m) < 0)
    { /* set v = v_{2k+1}, v1 = v_{2k+2} */
      v = Fl_add(Fl_mul(v,v1,N), mP, N);
      v1= Fl_add(Fl_mul(v1,v1,N),m2, N);
    }
    else
    {/* set v = v_{2k}, v1 = v_{2k+1} */
      v1= Fl_add(Fl_mul(v,v1,N),mP, N);
      v = Fl_add(Fl_mul(v,v,N), m2, N);
    }
  }
  return v;
}

static int
u_IsLucasPsP(ulong n)
{
  long i, v;
  ulong b, z, m2, m = n + 1;

  for (b=3, i=0;; b+=2, i++)
  {
    ulong c = b*b - 4; /* = 1 mod 4 */
    if (krouu(n % c, c) < 0) break;
    if (i == 64 && uissquarerem(n, &c)) return 0; /* oo loop if N = m^2 */
  }
  if (!m) return 0; /* neither 2^32-1 nor 2^64-1 are Lucas-pp */
  v = vals(m); m >>= v;
  z = u_LucasMod(m, b, n);
  if (z == 2) return 1;
  m2 = n - 2;
  if (z == m2) return 1;
  for (i=1; i<v; i++)
  {
    if (!z) return 1;
    z = Fl_add(Fl_mul(z,z, n), m2, n);
    if (z == 2) return 0;
  }
  return 0;
}
/* check that N not a square first (taken care of here, but inefficient)
 * assume N > 3 */
static int
IsLucasPsP(GEN N)
{
  GEN N_2, m, z;
  long i, v;
  ulong b;

  for (b=3, i=0;; b+=2, i++)
  {
    ulong c = b*b - 4; /* = 1 mod 4 */
    if (i == 64 && Z_issquare(N)) return 0; /* avoid oo loop if N = m^2 */
    if (krouu(umodiu(N,c), c) < 0) break;
  }
  m = addis(N,1); v = vali(m); m = shifti(m,-v);
  z = LucasMod(m, b, N);
  if (equaliu(z, 2)) return 1;
  N_2 = subis(N,2);
  if (equalii(z, N_2)) return 1;
  for (i=1; i<v; i++)
  {
    if (!signe(z)) return 1;
    z = modii(subis(sqri(z), 2), N);
    if (equaliu(z, 2)) return 0;
  }
  return 0;
}

/* assume u odd, u > 1 */
static int
iu_coprime(GEN N, ulong u)
{
  const ulong n = umodiu(N, u);
  return (n == 1 || ugcd(n, u) == 1);
}
/* assume u odd, u > 1 */
static int
uu_coprime(ulong n, ulong u)
{
  return (n == 1 || ugcd(n, u) == 1);
}

/* Fl_BSW_psp */
int
uisprime(ulong n)
{
  Fl_miller_t S;
  if (n < 103)
    switch(n)
    {
      case 2:
      case 3:
      case 5:
      case 7:
      case 11:
      case 13:
      case 17:
      case 19:
      case 23:
      case 29:
      case 31:
      case 37:
      case 41:
      case 43:
      case 47:
      case 53:
      case 59:
      case 61:
      case 67:
      case 71:
      case 73:
      case 79:
      case 83:
      case 89:
      case 97:
      case 101: return 1;
      default: return 0;
    }
  if (!(n & 1)) return 0;
#ifdef LONG_IS_64BIT
  /* 16294579238595022365 = 3*5*7*11*13*17*19*23*29*31*37*41*43*47*53
   *  7145393598349078859 = 59*61*67*71*73*79*83*89*97*101 */
  if (!uu_coprime(n, 16294579238595022365UL) ||
      !uu_coprime(n,  7145393598349078859UL)) return 0;
#else
  /* 4127218095 = 3*5*7*11*13*17*19*23*37
   * 3948078067 = 29*31*41*43*47*53
   * 4269855901 = 59*83*89*97*101
   * 1673450759 = 61*67*71*73*79 */
  if (!uu_coprime(n, 4127218095UL) ||
      !uu_coprime(n, 3948078067UL) ||
      !uu_coprime(n, 1673450759UL) ||
      !uu_coprime(n, 4269855901UL)) return 0;
#endif
  if (n < 10427) return 1;
  Fl_init_miller(&S, n);
  if (Fl_bad_for_base(&S, 2)) return 0;
  if (n < 1016801) switch(n) {
    case 42799: /* strong 2-pseudoprimes without prime divisors < 103 */
    case 49141:
    case 88357:
    case 90751:
    case 104653:
    case 130561:
    case 196093:
    case 220729:
    case 253241:
    case 256999:
    case 271951:
    case 280601:
    case 357761:
    case 390937:
    case 458989:
    case 486737:
    case 489997:
    case 514447:
    case 580337:
    case 741751:
    case 838861:
    case 873181:
    case 877099:
    case 916327:
    case 976873:
    case 983401: return 0;
    default: return 1;
  }
  return u_IsLucasPsP(n);
}

long
BSW_psp(GEN N)
{
  pari_sp av = avma;
  miller_t S;
  int k;

  if (typ(N) != t_INT) pari_err(arither1);
  if (signe(N) <= 0) return 0;
  if (lgefint(N) == 3) return uisprime((ulong)N[2]);
  if (!mod2(N)) return 0;
#ifdef LONG_IS_64BIT
  /* 16294579238595022365 = 3*5*7*11*13*17*19*23*29*31*37*41*43*47*53
   *  7145393598349078859 = 59*61*67*71*73*79*83*89*97*101 */
  if (!iu_coprime(N, 16294579238595022365UL) ||
      !iu_coprime(N,  7145393598349078859UL)) return 0;
#else
  /* 4127218095 = 3*5*7*11*13*17*19*23*37
   * 3948078067 = 29*31*41*43*47*53
   * 4269855901 = 59*83*89*97*101
   * 1673450759 = 61*67*71*73*79 */
  if (!iu_coprime(N, 4127218095UL) ||
      !iu_coprime(N, 3948078067UL) ||
      !iu_coprime(N, 1673450759UL) ||
      !iu_coprime(N, 4269855901UL)) return 0;
#endif
  /* no prime divisor < 103 */
  av = avma;
  init_miller(&S, N); 
  k = (!bad_for_base(&S, gen_2) && IsLucasPsP(N));
  avma = av; return k;
}

/***********************************************************************/
/**                                                                   **/
/**                       Pocklington-Lehmer                          **/
/**                        P-1 primality test                         **/
/** Crude implementation  BA 2000Apr21                                **/
/***********************************************************************/

/*assume n>=2*/
static ulong
pl831(GEN N, GEN p)
{
  pari_sp ltop = avma, av;
  ulong a;
  GEN Nmunp = diviiexact(addis(N,-1), p);
  av = avma;
  for(a = 2;; a++, avma = av)
  {
    GEN b = Fp_pow(utoipos(a), Nmunp, N);
    GEN c = Fp_pow(b,p,N), g = gcdii(addis(b,-1), N);
    if (!is_pm1(c)) return 0;
    if (is_pm1(g)) { avma=ltop; return a; }
    if (!equalii(g,N)) return 0;
  }
}
/* Assume N is a strong BSW pseudoprime
 *
 * flag 0: return gen_1 (prime), gen_0 (composite)
 * flag 1: return gen_0 (composite), gen_1 (small prime), matrix (large prime)
 *
 * The matrix has 3 columns, [a,b,c] with
 * a[i] prime factor of N-1,
 * b[i] witness for a[i] as in pl831
 * c[i] plisprime(a[i]) */
GEN
plisprime(GEN N, long flag)
{
  pari_sp ltop = avma;
  long i, l, t = typ(N);
  int eps;
  GEN C, F = NULL;

  if (t == t_VEC)
  { /* [ N, [p1,...,pk] ], pi list of pseudoprime divisors of N */
    F = gel(N,2);
    N = gel(N,1); t = typ(N);
  }
  if (t != t_INT) pari_err(arither1);
  if (DEBUGLEVEL>3) fprintferr("PL: proving primality of N = %Z\n", N);

  eps = cmpis(N,2);
  if (eps<=0) return eps? gen_0: gen_1;

  N = absi(N);
  if (!F)
  {
    F = (GEN)Z_factor_limit(addis(N,-1), sqrti(N))[1];
    if (DEBUGLEVEL>3) fprintferr("PL: N-1 factored!\n");
  }

  C = cgetg(4,t_MAT); l = lg(F);
  gel(C,1) = cgetg(l,t_COL);
  gel(C,2) = cgetg(l,t_COL);
  gel(C,3) = cgetg(l,t_COL);
  for(i=1; i<l; i++)
  {
    GEN p = gel(F,i), r;
    ulong witness = pl831(N,p);

    if (!witness) { avma = ltop; return gen_0; }
    gmael(C,1,i) = icopy(p);
    gmael(C,2,i) = utoi(witness);
    if (!flag) r = BSW_isprime(p)? gen_1: gen_0;
    else
    {
      if (BSW_isprime_small(p)) r = gen_1;
      else if (expi(p) > 250)   r = isprimeAPRCL(p)? gen_2: gen_0;
      else                      r = plisprime(p,flag);
    }
    gmael(C,3,i) = r;
    if (r == gen_0) pari_err(talker,"False prime number %Z in plisprime", p);
  }
  if (!flag) { avma = ltop; return gen_1; }
  return gerepileupto(ltop,C);
}

/***********************************************************************/
/**                                                                   **/
/**                       PRIMES IN SUCCESSION                        **/
/** (abstracted by GN 1998Aug21 mainly for use in ellfacteur() below) **/
/**                                                                   **/
/***********************************************************************/

/* map from prime residue classes mod 210 to their numbers in {0...47}.
 * Subscripts into this array take the form ((k-1)%210)/2, ranging from
 * 0 to 104.  Unused entries are */
#define NPRC 128		/* non-prime residue class */

static unsigned char prc210_no[] = {
  0, NPRC, NPRC, NPRC, NPRC, 1, 2, NPRC, 3, 4, NPRC, /* 21 */
  5, NPRC, NPRC, 6, 7, NPRC, NPRC, 8, NPRC, 9, /* 41 */
  10, NPRC, 11, NPRC, NPRC, 12, NPRC, NPRC, 13, 14, NPRC, /* 63 */
  NPRC, 15, NPRC, 16, 17, NPRC, NPRC, 18, NPRC, 19, /* 83 */
  NPRC, NPRC, 20, NPRC, NPRC, NPRC, 21, NPRC, 22, 23, NPRC, /* 105 */
  24, 25, NPRC, 26, NPRC, NPRC, NPRC, 27, NPRC, NPRC, /* 125 */
  28, NPRC, 29, NPRC, NPRC, 30, 31, NPRC, 32, NPRC, NPRC, /* 147 */
  33, 34, NPRC, NPRC, 35, NPRC, NPRC, 36, NPRC, 37, /* 167 */
  38, NPRC, 39, NPRC, NPRC, 40, 41, NPRC, NPRC, 42, NPRC, /* 189 */
  43, 44, NPRC, 45, 46, NPRC, NPRC, NPRC, NPRC, 47, /* 209 */
};

#if 0
/* map from prime residue classes mod 210 (by number) to their smallest
 * positive representatives */
static unsigned char prc210_rp[] = { /* 19 + 15 + 14 = [0..47] */
  1, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79,
  83, 89, 97, 101, 103, 107, 109, 113, 121, 127, 131, 137, 139, 143, 149,
  151, 157, 163, 167, 169, 173, 179, 181, 187, 191, 193, 197, 199, 209,
};
#endif

/* first differences of the preceding */
static unsigned char prc210_d1[] = {
  10, 2, 4, 2, 4, 6, 2, 6, 4, 2, 4, 6, 6, 2, 6, 4, 2, 6,
  4, 6, 8, 4, 2, 4, 2, 4, 8, 6, 4, 6, 2, 4, 6,
  2, 6, 6, 4, 2, 4, 6, 2, 6, 4, 2, 4, 2, 10, 2,
};

GEN
nextprime(GEN n)
{
  long rc, rc0, rcd, rcn;
  pari_sp av = avma;

  if (typ(n) != t_INT) n = gceil(n);
  if (typ(n) != t_INT) pari_err(arither1);
  if (signe(n) <= 0) { avma = av; return gen_2; }
  if (lgefint(n) <= 3)
  { /* check if n <= 7 */
    ulong k = n[2];
    if (k <= 2) { avma = av; return gen_2; }
    if (k == 3) { avma = av; return utoipos(3); }
    if (k <= 5) { avma = av; return utoipos(5); }
    if (k <= 7) { avma = av; return utoipos(7); }
  }
  /* here n > 7 */
  if (!mod2(n)) n = addsi(1,n);
  rc = rc0 = smodis(n, 210);
  /* find next prime residue class mod 210 */
  for(;;)
  {
    rcn = (long)(prc210_no[rc>>1]);
    if (rcn != NPRC) break;
    rc += 2; /* cannot wrap since 209 is coprime and rc odd */
  }
  if (rc > rc0) n = addsi(rc - rc0, n);
  /* now find an actual (pseudo)prime */
  for(;;)
  {
    if (BSW_psp(n)) break;
    rcd = prc210_d1[rcn];
    if (++rcn > 47) rcn = 0;
    n = addsi(rcd, n);
  }
  if (avma == av) return icopy(n);
  return gerepileuptoint(av, n);
}

GEN
precprime(GEN n)
{
  long rc, rc0, rcd, rcn;
  pari_sp av = avma;

  if (typ(n) != t_INT) n = gfloor(n);
  if (typ(n) != t_INT) pari_err(arither1);
  if (signe(n) <= 0) { avma = av; return gen_0; }
  if (lgefint(n) <= 3)
  { /* check if n <= 10 */
    ulong k = n[2];
    if (k <= 1)  { avma = av; return gen_0; }
    if (k == 2)  { avma = av; return gen_2; }
    if (k <= 4)  { avma = av; return utoipos(3); }
    if (k <= 6)  { avma = av; return utoipos(5); }
    if (k <= 10) { avma = av; return utoipos(7); }
  }
  /* here n >= 11 */
  if (!mod2(n)) n = addsi(-1,n);
  rc = rc0 = smodis(n, 210);
  /* find previous prime residue class mod 210 */
  for(;;)
  {
    rcn = (long)(prc210_no[rc>>1]);
    if (rcn != NPRC) break;
    rc -= 2; /* cannot wrap since 1 is coprime and rc odd */
  }
  if (rc < rc0) n = addsi(rc - rc0, n);
  /* now find an actual (pseudo)prime */
  for(;;)
  {
    if (BSW_psp(n)) break;
    if (--rcn < 0) rcn = 47;
    rcd = prc210_d1[rcn];
    n = addsi(-rcd, n);
  }
  if (avma == av) return icopy(n);
  return gerepileuptoint(av, n);
}

/* Find next single-word prime strictly larger than p.
 * If **d is non-NULL (somewhere in a diffptr), this is p + *(*d)++.
 * Otherwise imitate nextprime().
 * *rcn = NPRC or the correct residue class for the current p;  we'll use this
 * to track the current prime residue class mod 210 once we're out of range of
 * the diffptr table, and we'll update it before that if it isn't NPRC.
 * *q is incremented whenever q!=NULL and we wrap from 209 mod 210 to
 * 1 mod 210;  this makes sense
 * k =  second argument for miller(). --GN1998Aug22 */
ulong
snextpr(ulong p, byteptr *d, long *rcn, long *q, long k)
{
  ulong n;
  if (**d)
  {
    byteptr dd = *d;
    long d1 = 0;

    NEXT_PRIME_VIADIFF(d1,dd);
    if (*rcn != NPRC)
    {
      long rcn0 = *rcn;
      while (d1 > 0)
      {
	d1 -= prc210_d1[*rcn];
	if (++*rcn > 47) { *rcn = 0; if (q) (*q)++; }
      }
      if (d1 < 0)
      {
	fprintferr("snextpr: %lu != prc210_rp[%ld] mod 210\n", p, rcn0);
	pari_err(bugparier, "[caller of] snextpr");
      }
    }
    NEXT_PRIME_VIADIFF(p,*d);
    return p;
  }
  /* we are beyond the diffptr table */
  if (*rcn == NPRC)
  { /* initialize */
    *rcn = prc210_no[(p % 210) >> 1];
    if (*rcn == NPRC)
    {
      fprintferr("snextpr: %lu should have been prime but isn\'t\n", p);
      pari_err(bugparier, "[caller of] snextpr");
    }
  }
  /* look for the next one */
  n = p + prc210_d1[*rcn];
  if (++*rcn > 47) *rcn = 0;
  while (!Fl_miller(n, k))
  {
    n += prc210_d1[*rcn];
    if (++*rcn > 47) { *rcn = 0; if (q) (*q)++; }
    if (n <= 11)		/* wraparound mod 2^BITS_IN_LONG */
    {
      fprintferr("snextpr: integer wraparound after prime %lu\n", p);
      pari_err(bugparier, "[caller of] snextpr");
    }
  }
  return n;
}

/***********************************************************************/
/**                                                                   **/
/**                 FACTORIZATION (ECM) -- GN Jul-Aug 1998            **/
/**   Integer factorization using the elliptic curves method (ECM).   **/
/**   ellfacteur() returns a non trivial factor of N, assuming N>0,   **/
/**   is composite, and has no prime divisor below 2^14 or so.        **/
/**   Thanks to Paul Zimmermann for much helpful advice and to        **/
/**   Guillaume Hanrot and Igor Schein for intensive testing          **/
/**                                                                   **/
/***********************************************************************/

static GEN N, gl;
#define nbcmax 64		/* max number of simultaneous curves */
#define bstpmax 1024		/* max number of baby step table entries */

/* addition/doubling/multiplication of a point on an `elliptic curve
 * mod N' may result in one of three things:
 * - a new bona fide point
 * - a point at infinity  (denominator divisible by N)
 * - a point at infinity mod some nontrivial factor of N but finite mod some
 *   other factor  (betraying itself by a denominator which has nontrivial gcd
 *   with N, and this is of course what we want).
 */
/* (In the second case, addition/doubling will simply abort, copying one
 * of the summands to the destination array of points unless they coincide.
 * Multiplication will stop at some unpredictable intermediate stage:  The
 * destination will contain _some_ multiple of the input point, but not
 * necessarily the desired one, which doesn't matter.  As long as we're
 * multiplying (B1 phase) we simply carry on with the next multiplier.
 * During the B2 phase, the only additions are the giant steps, and the
 * worst that can happen here is that we lose one residue class mod 210
 * of prime multipliers on 4 of the curves, so again, we ignore the problem
 * and just carry on.) */
/* The idea is:  Select a handful of curves mod N and one point P on each of
 * them.  Try to compute, for each such point, the multiple [M]P = Q where
 * M is the product of all powers <= B2 of primes <= nextprime(B1), for some
 * suitable B1 and B2.  Then check whether multiplying Q by one of the
 * primes < nextprime(B2) would betray a factor.  This second stage proceeds
 * by looking separately at the primes in each residue class mod 210, four
 * curves at a time, and stepping additively to ever larger multipliers,
 * by comparing X coordinates of points which we would need to add in order
 * to reach another prime multiplier in the same residue class.  `Comparing'
 * means that we accumulate a product of differences of X coordinates, and
 * from time to time take a gcd of this product with N.
 */
/* Montgomery's trick (hiding the cost of computing inverses mod N at a
 * price of three extra multiplications mod N, by working on up to 64 or
 * even 128 points in parallel) is used heavily. */

/* *** auxiliary functions for ellfacteur: *** */

/* Parallel addition on nbc curves, assigning the result to locations at and
 * following *X3, *Y3.  Safe to be called with X3,Y3 equal to X2,Y2  (_not_
 * to X1,Y1).  It is also safe to overwrite Y2 with X3.  (If Y coords of
 * result not desired, set Y3=NULL.)  If nbc1 < nbc, the first summand is
 * assumed to hold only nbc1 distinct points, which are repeated as often
 * as we need them  (useful for adding one point on each of a few curves
 * to several other points on the same curves).
 * Return 0 when successful, 1 when we hit a denominator divisible by N,
 * and 2 when gcd(denominator, N) is a nontrivial factor of N, which will
 * be preserved in gl.
 * Stack space is bounded by a constant multiple of lgefint(N)*nbc.
 * (Phase 2 creates 12 items on the stack, per iteration, of which
 * four are twice as long and one is thrice as long as N -- makes 18 units
 * per iteration.  Phase  1 creates 4 units.  Total can be as large as
 * about 4*nbcmax + 18*8 units.  And elladd2() is just as bad, and
 * elldouble() comes to about 3*nbcmax + 29*8 units.  A few strategic garbage
 * collections every 8 iterations may help when nbc is large.) */

static int
elladd0(long nbc, long nbc1,
	GEN *X1, GEN *Y1, GEN *X2, GEN *Y2, GEN *X3, GEN *Y3)
{
  GEN W[2*nbcmax], *A = W+nbc; /* W[0],A[0] unused */
  long i;
  pari_sp av = avma, tetpil;
  ulong mask = ~0UL;

  /* actually, this is only ever called with nbc1==nbc or nbc1==4, so: */
  if (nbc1 == 4) mask = 3;
  else if (nbc1 < nbc) pari_err(bugparier, "[caller of] elladd0");

  W[1] = subii(X1[0], X2[0]);
  for (i=1; i<nbc; i++)
  {
    A[i] = subii(X1[i&mask], X2[i]); /* don't waste time reducing mod N here */
    W[i+1] = modii(mulii(A[i], W[i]), N);
  }
  tetpil = avma;

  if (!invmod(W[nbc], N, &gl))
  { /* hit infinity */
    if (!equalii(N,gl)) return 2;
    if (X2 != X3)
    { /* cannot add on one of the curves mod N:  make sure X3 contains
       * something useful before letting caller proceed */
      long k;
      for (k = 2*nbc; k--; ) affii(X2[k],X3[k]);
    }
    avma = av; return 1;
  }

  while (i--) /* nbc times, actually */
  {
    pari_sp av2 = avma;
    GEN t, L = modii(mulii(subii(Y1[i&mask], Y2[i]),
	 		   i?mulii(gl, W[i]):gl), N);
    t = subii(sqri(L), addii(X2[i], X1[i&mask]));
    affii(modii(t, N), X3[i]);
    if (Y3) {
      t = subii(mulii(L, subii(X1[i&mask], X3[i])), Y1[i&mask]);
      affii(modii(t, N), Y3[i]);
    }
    if (!i) break;
    avma = av2; gl = modii(mulii(gl, A[i]), N);
    if (!(i&7)) gl = gerepileuptoint(tetpil, gl);
  }
  avma = av; return 0;
}

/* Shortcut, for use in cases where Y coordinates follow their corresponding
 * X coordinates, and first summand doesn't need to be repeated */
static int
elladd(long nbc, GEN *X1, GEN *X2, GEN *X3) {
  return elladd0(nbc, nbc, X1, X1+nbc, X2, X2+nbc, X3, X3+nbc);
}

/* As elladd except it does twice as many additions (and thus hides even more
 * of the cost of the modular inverse); the net effect is the same as
 * elladd(nbc,X1,X2,X3) followed by elladd(nbc,X4,X5,X6).  Safe to have
 * X2==X3, X5==X6, or X1 or X2 coincide with X4 or X5, in any order. */
static int
elladd2(long nbc, GEN *X1, GEN *X2, GEN *X3, GEN *X4, GEN *X5, GEN *X6)
{
  GEN *Y1 = X1+nbc, *Y2 = X2+nbc, *Y3 = X3+nbc;
  GEN *Y4 = X4+nbc, *Y5 = X5+nbc, *Y6 = X6+nbc;
  GEN W[4*nbcmax], *A = W+2*nbc; /* W[0],A[0] unused */
  long i, j;
  pari_sp av=avma, tetpil;

  W[1] = subii(X1[0], X2[0]);
  for (i=1; i<nbc; i++)
  {
    A[i] = subii(X1[i], X2[i]);	/* don't waste time reducing mod N here */
    W[i+1] = modii(mulii(A[i], W[i]), N);
  }
  for (j=0; j<nbc; i++,j++)
  {
    A[i] = subii(X4[j], X5[j]);
    W[i+1] = modii(mulii(A[i], W[i]), N);
  }
  tetpil = avma;

  /* if gl != N we have a factor */
  if (!invmod(W[2*nbc], N, &gl))
  { /* hit infinity */
    if (!equalii(N,gl)) return 2;
    if (X2 != X3)
    { /* cannot add on one of the curves mod N:  make sure X3 contains
       * something useful before letting caller proceed */
      long k;
      for (k = 2*nbc; k--; ) affii(X2[k],X3[k]);
    }
    if (X5 != X6)
    { /* same for X6 */
      long k;
      for (k = 2*nbc; k--; ) affii(X5[k],X6[k]);
    }
    avma = av; return 1;
  }

  while (j--) /* nbc times, actually */
  {
    pari_sp av2 = avma;
    GEN t, L;
    i--;
    L = modii(mulii(subii(Y4[j], Y5[j]), mulii(gl, W[i])), N);
    t = subii(sqri(L), addii(X5[j], X4[j]));         affii(modii(t,N), X6[j]);
    t = subii(mulii(L, subii(X4[j], X6[j])), Y4[j]); affii(modii(t,N), Y6[j]);
    avma = av2; gl = modii(mulii(gl, A[i]), N);
    if (!(i&7)) gl = gerepileuptoint(tetpil, gl);
  }
  while (i--) /* nbc times */
  {
    pari_sp av2 = avma;
    GEN t, L = modii(mulii(subii(Y1[i], Y2[i]),
			   i?mulii(gl, W[i]):gl), N);
    t = subii(sqri(L), addii(X2[i], X1[i]));         affii(modii(t,N), X3[i]);
    t = subii(mulii(L, subii(X1[i], X3[i])), Y1[i]); affii(modii(t,N), Y3[i]);
    if (!i) break;
    avma = av2; gl = modii(mulii(gl, A[i]), N);
    if (!(i&7)) gl = gerepileuptoint(tetpil, gl);
  }
  avma = av; return 0;
}

/* Parallel doubling on nbc curves, assigning the result to locations at
 * and following *X2.  Safe to be called with X2 equal to X1.  Return
 * value as for elladd.  If we find a point at infinity mod N,
 * and if X1 != X2, we copy the points at X1 to X2. */
static int
elldouble(long nbc, GEN *X1, GEN *X2)
{
  GEN *Y1 = X1+nbc, *Y2 = X2+nbc;
  GEN W[nbcmax+1]; /* W[0] unused */
  long i;
  pari_sp av = avma, tetpil;
  /*W[0] = gen_1;*/ W[1] = Y1[0];
  for (i=1; i<nbc; i++) W[i+1] = modii(mulii(Y1[i], W[i]), N);
  tetpil = avma;

  if (!invmod(W[nbc], N, &gl))
  {
    if (!equalii(N,gl)) return 2;
    if (X1 != X2)
    {
      long k;
      /* cannot double on one of the curves mod N:  make sure X2 contains
       * something useful before letting the caller proceed
       */
      for (k = 2*nbc; k--; ) affii(X1[k],X2[k]);
    }
    avma = av; return 1;
  }

  while (i--) /* nbc times, actually */
  {
    pari_sp av2;
    GEN v, w, L, GL = gl;

    if (i) gl = modii(mulii(gl, Y1[i]), N);
    av2 = avma;
    L = modii(mulii(addsi(1, mulsi(3, sqri(X1[i]))),
                    i? mulii(GL,W[i]): GL), N);
    if (signe(L)) /* half of zero is still zero */
      L = shifti(mod2(L)? addii(L, N): L, -1);
    v = modii(subii(sqri(L), shifti(X1[i],1)), N);
    w = modii(subii(mulii(L, subii(X1[i], v)), Y1[i]), N);
    affii(v, X2[i]);
    affii(w, Y2[i]); avma = av2;
    if (!(i&7) && i) gl = gerepileuptoint(tetpil, gl);
  }
  avma = av; return 0;
}

/* Parallel multiplication by an odd prime k on nbc curves, storing the
 * result to locations at and following *X2.  Safe to be called with X2 = X1.
 * Return values as elladd. Uses (a simplified variant of) Peter Montgomery's
 * PRAC (PRactical Addition Chain) algorithm;
 * see ftp://ftp.cwi.nl/pub/pmontgom/Lucas.ps.gz .
 * With thanks to Paul Zimmermann for the reference.  --GN1998Aug13 */

/* k>2 assumed prime, XAUX = scratchpad */
static int
ellmult(long nbc, ulong k, GEN *X1, GEN *X2, GEN *XAUX)
{
  ulong r, d, e, e1;
  long i;
  int res;
  GEN *A = X2, *B = XAUX, *S, *T = XAUX + 2*nbc;

  for (i = 2*nbc; i--; ) affii(X1[i], XAUX[i]);

  /* first doubling picks up X1;  after this we'll be working in XAUX and
   * X2 only, mostly via A and B and T */
  if ((res = elldouble(nbc, X1, X2)) != 0) return res;

  /* split the work at the golden ratio */
  r = (ulong)(k*0.61803398875 + .5);
  d = k - r;
  e = r - d; /* d+e == r, so no danger of ofl below */

  while (d != e)
  { /* apply one of the nine transformations from PM's Table 4. First
     * figure out which, and then go into an eight-way switch, because
     * some of the transformations are similar enough to share code. */
    if (d <= e + (e>>2))	/* floor(1.25*e) */
    {
      if ((d+e)%3 == 0) { i = 0; goto apply; } /* rule 1 */
      if ((d-e)%6 == 0) { i = 1; goto apply; } /* rule 2 */
    } /* else fall through */
    /* d <= 4*e but no ofl */
    if ((d+3)>>2 <= e){ i = 2; goto apply; }	/* rule 3, common case */
    if ((d&1)==(e&1)) { i = 1; goto apply; }	/* rule 4 = rule 2 */
    if (!(d&1))       { i = 3; goto apply; }	/* rule 5 */
    if (d%3 == 0)     { i = 4; goto apply; }	/* rule 6 */
    if ((d+e)%3 == 0) { i = 5; goto apply; }	/* rule 7 */
    if ((d-e)%3 == 0) { i = 6; goto apply; }	/* rule 8 */
    /* when we get here, e must be even, for otherwise one of rules 4,5
     * would have applied */
    i = 7;			/* rule 9 */

  apply:
    switch(i)			/* i takes values in {0,...,7} here */
    {
    case 0:			/* rule 1 */
      e1 = d - e; d = (d + e1)/3; e = (e - e1)/3;
      if ( (res = elladd(nbc, A, B, T)) ) return res;
      if ( (res = elladd2(nbc, T, A, A, T, B, B)) != 0) return res;
      break;			/* end of rule 1 */
    case 1:			/* rules 2 and 4, part 1 */
      d -= e;
      if ( (res = elladd(nbc, A, B, B)) ) return res;
      /* FALL THROUGH */
    case 3:			/* rule 5, and 2nd part of rules 2 and 4 */
      d >>= 1;
      if ( (res = elldouble(nbc, A, A)) ) return res;
      break;			/* end of rules 2, 4, and 5 */
    case 4:			/* rule 6 */
      d /= 3;
      if ( (res = elldouble(nbc, A, T)) ) return res;
      if ( (res = elladd(nbc, T, A, A)) ) return res;
      /* FALL THROUGH */
    case 2:			/* rule 3, and 2nd part of rule 6 */
      d -= e;
      if ( (res = elladd(nbc, A, B, B)) ) return res;
      break;			/* end of rules 3 and 6 */
    case 5:			/* rule 7 */
      d = (d - e - e)/3;
      if ( (res = elldouble(nbc, A, T)) ) return res;
      if ( (res = elladd2(nbc, T, A, A, T, B, B)) != 0) return res;
      break;			/* end of rule 7 */
    case 6:			/* rule 8 */
      d = (d - e)/3;
      if ( (res = elladd(nbc, A, B, B)) ) return res;
      if ( (res = elldouble(nbc, A, T)) ) return res;
      if ( (res = elladd(nbc, T, A, A)) ) return res;
      break;			/* end of rule 8 */
    case 7:			/* rule 9 */
      e >>= 1;
      if ( (res = elldouble(nbc, B, B)) ) return res;
      break;			/* end of rule 9 */
    default: break; /* notreached */
    }
    /* end of Table 4 processing */

    /* swap d <-> e and A <-> B if necessary */
    if (d < e) { r = d; d = e; e = r; S = A; A = B; B = S; }
  } /* while */
  return elladd(nbc, XAUX, X2, X2);
}

/* Auxiliary routines need < (3*nbc+240)*tf words on the PARI stack, in
 * addition to the spc*(tf+1) words occupied by our main table.
 * If stack space is already tight, use the heap & newbloc(). */
static GEN*
alloc_scratch(long nbc, long spc, long tf)
{
  long i, tw = evallg(tf) | evaltyp(t_INT), len = spc + 385 + spc*tf;
  GEN *X, w;
  if ((long)((GEN)avma - (GEN)bot) < len + (3*nbc + 240)*tf)
  {
    if (DEBUGLEVEL>4) fprintferr("ECM: stack tight, using heap space\n");
    X = (GEN*)newbloc(len);
  } else
    X = (GEN*)new_chunk(len);
  /* hack for X[i] = cgeti(tf). X = current point in B1 phase */
  w = (GEN)(X + spc + 385);
  for (i = spc-1; i >= 0; i--) { X[i] = w; *w = tw; w += tf; }
  return X;
}

/* PRAC implementation notes - main changes against the paper version:
 * (1) The general function  [m+n]P = f([m]P,[n]P,[m-n]P)  collapses  (for
 * m!=n)  to an elladd() which does not depend on the third argument;  and
 * thus all references to the third variable (C in the paper) can be elimi-
 * nated. (2) Since our multipliers are prime, the outer loop of the paper
 * version executes only once, and thus is invisible above. (3) The first
 * step in the inner loop of the paper version will always be rule 3, but
 * the addition requested by this rule amounts to a doubling, and it will
 * always be followed by a swap, so we have unrolled this first iteration.
 * (4) Some simplifications in rules 6 and 7 are possible given the above,
 * and we can save one addition in each of the two cases.  NB one can show
 * that none of the other elladd()s in the loop can ever turn out to de-
 * generate into an elldouble. (5) I tried to optimize for rule 3, which
 * is used far more frequently than all others together, but it didn't
 * improve things, so I removed the nested tight loop again.  --GN
 */

/* The main loop body of ellfacteur() runs slightly _slower_  under PRAC than
 * under a straightforward left-shift binary multiplication algorithm when
 * N has <30 digits and B1 is small;  PRAC wins when N and B1 get larger.
 * Weird. --GN
 */

/* memory layout in ellfacteur():  We'll have a large-ish array of GEN
 * pointers, and one huge chunk of memory containing all the actual GEN
 * (t_INT) objects.
 * nbc is constant throughout the invocation.
 */
/* The B1 stage of each iteration through the main loop needs little
 * space:  enough for the X and Y coordinates of the current points,
 * and twice as much again as scratchpad for ellmult().
 */
/* The B2 stage, starting from some current set of points Q, needs, in
 * succession:
 * - space for [2]Q, [4]Q, ..., [10]Q, and [p]Q for building the helix;
 * - space for 48*nbc X and Y coordinates to hold the helix.  Now this
 * could re-use [2]Q,...,[8]Q, but only with difficulty, since we don't
 * know in advance which residue class mod 210 our p is going to be in.
 * It can and should re-use [p]Q, though;
 * - space for (temporarily [30]Q and then) [210]Q, [420]Q, and several
 * further doublings until the giant step multiplier is reached.  This
 * _can_ re-use the remaining cells from above.  The computation of [210]Q
 * will have been the last call to ellmult() within this iteration of the
 * main loop, so the scratchpad is now also free to be re-used.  We also
 * compute [630]Q by a parallel addition;  we'll need it later to get the
 * baby-step table bootstrapped a little faster.
 */
/* Finally, for no more than 4 curves at a time, room for up to 1024 X
 * coordinates only  (the Y coordinates needed whilst setting up this baby
 * step table are temporarily stored in the upper half, and overwritten
 * during the last series of additions).
 */
/* Graphically:  after end of B1 stage  (X,Y are the coords of Q):
 * +--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--
 * | X Y |  scratch  | [2]Q| [4]Q| [6]Q| [8]Q|[10]Q|    ...    | ...
 * +--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--
 * *X    *XAUX *XT   *XD                                       *XB
 *
 * [30]Q is computed from [10]Q.  [210]Q can go into XY, etc:
 * +--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--
 * |[210]|[420]|[630]|[840]|[1680,3360,6720,...,2048*210]      |bstp table...
 * +--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--
 * *X    *XAUX *XT   *XD      [*XG, somewhere here]            *XB .... *XH
 *
 * So we need (13 + 48) * 2 * nbc slots here, and another 4096 slots for
 * the baby step table (not all of which will be used when we start with a
 * small B1, but it's better to allocate and initialize ahead of time all
 * the slots that might be needed later).
 */
/* Note on memory locality:  During the B2 phase, accesses to the helix
 * (once it is set up)  will be clustered by curves  (4 out of nbc at a time).
 * Accesses to the baby steps table will wander from one end of
 * the array to the other and back, one such cycle per giant step, and
 * during a full cycle we would expect on the order of 2E4 accesses when
 * using the largest giant step size.  Thus we shouldn't be doing too bad
 * with respect to thrashing a (512KBy) L2 cache.  However, we don't want
 * the baby step table to grow larger than this, even if it would reduce
 * the number of E.C. operations by a few more per cent for very large B2,
 * lest cache thrashing slow down everything disproportionally. --GN
 */

/* parameters for miller() via snextpr(), for use by ellfacteur() */
#define miller_k1 16		/* B1 phase, foolproof below 10^12 */
#define miller_k2 1		/* B2 phase, not foolproof, much faster */
/* (miller_k2 will let thousands of composites slip through, which doesn't
 * harm ECM, but ellmult() during the B1 phase should only be fed primes
 * which really are prime)
 */
/* ellfacteur() has been re-tuned to be useful as a first stage before
 * MPQS, especially for _large_ arguments, when insist is false, and now
 * also for the case when insist is true, vaguely following suggestions
 * by Paul Zimmermann  (see http://www.loria.fr/~zimmerma/ and especially
 * http://www.loria.fr/~zimmerma/records/ecmnet.html). --GN 1998Jul,Aug */
GEN
ellfacteur(GEN n, int insist)
{
  static ulong TB1[] = { /* table revised, cf. below 1998Aug15 --GN */
    142,172,208,252,305,370,450,545,661,801,972,1180,1430,
    1735,2100,2550,3090,3745,4540,5505,6675,8090,9810,11900,
    14420,17490,21200,25700,31160,37780UL,45810UL,55550UL,67350UL,
    81660UL,99010UL,120050UL,145550UL,176475UL,213970UL,259430UL,
    314550UL,381380UL,462415UL,560660UL,679780UL,824220UL,999340UL,
    1211670UL,1469110UL,1781250UL,2159700UL,2618600UL,3175000UL,
    3849600UL,4667500UL,5659200UL,6861600UL,8319500UL,10087100UL,
    12230300UL,14828900UL,17979600UL,21799700UL,26431500UL,
    32047300UL,38856400UL,	/* 110 times that still fits into 32bits */
#ifdef LONG_IS_64BIT
    47112200UL,57122100UL,69258800UL,83974200UL,101816200UL,
    123449000UL,149678200UL,181480300UL,220039400UL,266791100UL,
    323476100UL,392204900UL,475536500UL,576573500UL,699077800UL,
    847610500UL,1027701900UL,1246057200UL,1510806400UL,1831806700UL,
    2221009800UL,2692906700UL,3265067200UL,3958794400UL,4799917500UL,
    /* the only reason to stop here is that I got bored  (and that users will
     * get bored watching their 64bit machines churning on such large numbers
     * for month after month).  Someone can extend this table when the hardware
     * has gotten 100 times faster than now --GN */
#endif
    };
  static ulong TB1_for_stage[] = { /* table revised 1998Aug11 --GN.
    * Start a little below the optimal B1 for finding factors which would just
    * have been missed by pollardbrent(), and escalate gradually, changing
    * curves sufficiently frequently to give good coverage of the small factor
    * ranges.  Entries grow a bit faster than what Paul says would be optimal
    * but a table instead of a 2D array keeps the code simple */
    500,520,560,620,700,800,900,1000,1150,1300,1450,1600,1800,2000,
    2200,2450,2700,2950,3250,3600,4000,4400,4850,5300,5800,6400,
    7100,7850,8700,9600,10600,11700,12900,14200,15700,17300,
    19000,21000,23200,25500,28000,31000,34500UL,38500UL,43000UL,
    48000UL,53800UL,60400UL,67750UL,76000UL,85300UL,95700UL,
    107400UL,120500UL,135400UL,152000UL,170800UL,191800UL,215400UL,
    241800UL,271400UL,304500UL,341500UL,383100UL,429700UL,481900UL,
    540400UL,606000UL,679500UL,761800UL,854100UL,957500UL,1073500UL,
  };
  long nbc,nbc2,dsn,dsnmax,rep,spc,gse,gss,rcn,rcn0,bstp,bstp0;
  long a, i, j, k, size = expi(n) + 1, tf = lgefint(n);
  ulong B1,B2,B2_p,B2_rt,m,p,p0,dp;
  GEN *X,*XAUX,*XT,*XD,*XG,*YG,*XH,*XB,*XB2,*Xh,*Yh,*Xb;
  GEN res = cgeti(tf);
  pari_sp av1, avtmp, av = avma;
  int rflag;

  N = n; /* make n known to auxiliary functions */
  /* determine where we'll start, how long we'll persist, and how many
   * curves we'll use in parallel */
  if (insist)
  {
    dsnmax = (size >> 2) - 10;
    if (dsnmax < 0) dsnmax = 0;
#ifdef LONG_IS_64BIT
    else if (dsnmax > 90) dsnmax = 90;
#else
    else if (dsnmax > 65) dsnmax = 65;
#endif
    dsn = (size >> 3) - 5;
    if (dsn < 0) dsn = 0;
    else if (dsn > 47) dsn = 47;
    /* pick up the torch where non-insistent stage would have given up */
    nbc = dsn + (dsn >> 2) + 9;	/* 8 or more curves in parallel */
    nbc &= ~3; /* 4 | nbc */
    if (nbc > nbcmax) nbc = nbcmax;
    a = 1 + (nbcmax<<7)*(size&0xffff); /* seed for choice of curves */
    rep = 0; /* gcc -Wall */
  }
  else
  {
    dsn = (size - 140) >> 3;
    if (dsn > 12) dsn = 12;
    dsnmax = 72;
    if (dsn < 0) /* < 140 bits: decline the task */
    {
#ifdef __EMX__
      /* MPQS's disk access under DOS/EMX would be abysmally slow, so... */
      dsn = 0;
      rep = 20;
      nbc = 8;
#else
      if (DEBUGLEVEL >= 4)
	fprintferr("ECM: number too small to justify this stage\n");
      avma = av; return NULL;
#endif
    }
    else
    {
      rep = (size <= 248 ?
	     (size <= 176 ? (size - 124) >> 4 : (size - 148) >> 3) :
	     (size - 224) >> 1);
      nbc = ((size >> 3) << 2) - 80;
      if (nbc < 8) nbc = 8;
      else if (nbc > nbcmax) nbc = nbcmax;
#ifdef __EMX__
      rep += 20;
#endif
    }

    /* it may be convenient to use disjoint sets of curves for the non-insist
     * and insist phases;  moreover, repeated calls acting on factors of the
     * same original number should try to use fresh curves.
     * The following achieves this */
    a = 1 + (nbcmax<<3)*(size & 0xf);
  }
  if (dsn > dsnmax) dsn = dsnmax;

  if (DEBUGLEVEL >= 4)
  {
    (void)timer2();
    fprintferr("ECM: working on %ld curves at a time; initializing", nbc);
    if (!insist)
    {
      if (rep == 1) fprintferr(" for one round");
      else          fprintferr(" for up to %ld rounds", rep);
    }
    fprintferr("...\n");
  }

  nbc2 = nbc << 1;
  spc = (13 + 48) * nbc2 + bstpmax * 4;
  X = alloc_scratch(nbc, spc, tf);
  XAUX = X    + nbc2;	 /* scratchpad for ellmult() */
  XT   = XAUX + nbc2;	 /* ditto, will later hold [3*210]Q */
  XD   = XT   + nbc2;	 /* room for various multiples */
  XB   = XD   + 10*nbc2; /* start of baby steps table */
  XB2  = XB   + 2 * bstpmax; /* middle of baby steps table */
  XH   = XB2  + 2 * bstpmax; /* end of bstps table, start of helix */
  Xh   = XH   + 48*nbc2; /* little helix, X coords */
  Yh   = XH   + 192;	 /* ditto, Y coords */
  /* XG will be set inside the main loop, since it depends on B2 */

  /* Xh range of 384 pointers not set; these will later duplicate the pointers
   * in the XH range, 4 curves at a time.  Some of the cells reserved here for
   * the XB range will never be used, instead, we'll warp the pointers to
   * connect to (read-only) GENs in the X/XD range; it would be complicated to
   * skip them here to conserve merely a few KBy of stack or heap space. */

  /* ECM MAIN LOOP */
  for(;;)
  {
    byteptr d0, d = diffptr;
    
    rcn = NPRC; /* multipliers begin at the beginning */
    /* pick curves & bounds */
    for (i = nbc2; i--; ) affui(a++, X[i]);
    B1 = insist ? TB1[dsn] : TB1_for_stage[dsn];
    B2 = 110*B1;
    B2_rt = (ulong)(sqrt((double)B2));
    /* pick giant step exponent and size.
     * With 32 baby steps, a giant step corresponds to 32*420 = 13440, appro-
     * priate for the smallest B2s.  With 1024, a giant step will be 430080;
     * this will be appropriate for B1 >~ 42000, where 512 baby steps would
     * imply roughly the same number of E.C. additions.
     */
    gse = B1 < 656
            ? (B1 < 200? 5: 6)
            : (B1 < 10500
              ? (B1 < 2625? 7: 8)
              : (B1 < 42000? 9: 10));
    gss = 1UL << gse;
    XG = XT + gse*nbc2;	/* will later hold [2^(gse+1)*210]Q */
    YG = XG + nbc;

    if (DEBUGLEVEL >= 4) {
      fprintferr("ECM: time = %6ld ms\nECM: dsn = %2ld,\tB1 = %4lu,",
                 timer2(), dsn, B1);
      fprintferr("\tB2 = %6lu,\tgss = %4ld*420\n", B2, gss);
    }
    p = 0;
    NEXT_PRIME_VIADIFF(p,d);

    /* ---B1 PHASE--- */
    /* treat p=2 separately */
    B2_p = B2 >> 1;
    for (m=1; m<=B2_p; m<<=1)
    {
      if ((rflag = elldouble(nbc, X, X)) > 1) goto fin;
      else if (rflag) break;
    }
    /* p=3,...,nextprime(B1) */
    while (p < B1 && p <= B2_rt)
    {
      pari_sp av = avma;
      p = snextpr(p, &d, &rcn, NULL, miller_k1);
      B2_p = B2/p;		/* beware integer overflow on 32-bit CPUs */
      for (m=1; m<=B2_p; m*=p)
      {
	if ((rflag = ellmult(nbc, p, X, X, XAUX)) > 1) goto fin;
	else if (rflag) break;
        avma = av;
      }
      avma = av;
    }
    /* primes p larger than sqrt(B2) appear only to the 1st power */
    while (p < B1)
    {
      pari_sp av = avma;
      p = snextpr(p, &d, &rcn, NULL, miller_k1);
      if (ellmult(nbc, p, X, X, XAUX) > 1) goto fin; /* p^2 > B2: no loop */
      avma = av;
    }
    if (DEBUGLEVEL >= 4) {
      fprintferr("ECM: time = %6ld ms, B1 phase done, ", timer2());
      fprintferr("p = %lu, setting up for B2\n", p);
    }

    /* ---B2 PHASE--- */
    /* compute [2]Q,...,[10]Q, which we need to build the helix */
    if (elldouble(nbc, X, XD) > 1)
      goto fin;	/* [2]Q */
    if (elldouble(nbc, XD, XD + nbc2) > 1)
      goto fin;	/* [4]Q */
    if (elladd(nbc, XD, XD + nbc2, XD + (nbc<<2)) > 1)
      goto fin;	/* [6]Q */
    if (elladd2(nbc,
		XD, XD + (nbc<<2), XT + (nbc<<3),
		XD + nbc2, XD + (nbc<<2), XD + (nbc<<3)) > 1)
      goto fin;	/* [8]Q and [10]Q */
    if (DEBUGLEVEL >= 7) fprintferr("\t(got [2]Q...[10]Q)\n");

    /* get next prime (still using the foolproof test) */
    p = snextpr(p, &d, &rcn, NULL, miller_k1);
    /* make sure we have the residue class number (mod 210) */
    if (rcn == NPRC)
    {
      rcn = prc210_no[(p % 210) >> 1];
      if (rcn == NPRC)
      {
	fprintferr("ECM: %lu should have been prime but isn\'t\n", p);
	pari_err(bugparier, "ellfacteur");
      }
    }

    /* compute [p]Q and put it into its place in the helix */
    if (ellmult(nbc, p, X, XH + rcn*nbc2, XAUX) > 1) goto fin;
    if (DEBUGLEVEL >= 7)
      fprintferr("\t(got [p]Q, p = %lu = prc210_rp[%ld] mod 210)\n", p, rcn);

    /* save current p, d, and rcn;  we'll need them more than once below */
    p0 = p;
    d0 = d;
    rcn0 = rcn;	/* remember where the helix wraps */
    bstp0 = 0;	/* p is at baby-step offset 0 from itself */

    /* fill up the helix, stepping forward through the prime residue classes
     * mod 210 until we're back at the r'class of p0.  Keep updating p so
     * that we can print meaningful diagnostics if a factor shows up;  but
     * don't bother checking which of these p's are in fact prime */
    for (i = 47; i; i--) /* 47 iterations */
    {
      p += (dp = (ulong)prc210_d1[rcn]);
      if (rcn == 47)
      {	/* wrap mod 210 */
	if (elladd(nbc, XT + dp*nbc, XH + rcn*nbc2, XH) > 1) goto fin;
	rcn = 0; continue;
      }
      if (elladd(nbc, XT + dp*nbc, XH + rcn*nbc2, XH + rcn*nbc2 + nbc2) > 1)
	goto fin;
      rcn++;
    }
    if (DEBUGLEVEL >= 7) fprintferr("\t(got initial helix)\n");

    /* compute [210]Q etc, which will be needed for the baby step table */
    if (ellmult(nbc, 3, XD + (nbc<<3), X, XAUX) > 1) goto fin;
    if (ellmult(nbc, 7, X, X, XAUX) > 1) goto fin; /* [210]Q */
    /* this was the last call to ellmult() in the main loop body; may now
     * overwrite XAUX and slots XD and following */
    if (elldouble(nbc, X, XAUX) > 1) goto fin; /* [420]Q */
    if (elladd(nbc, X, XAUX, XT) > 1) goto fin;/* [630]Q */
    if (elladd(nbc, X, XT, XD) > 1) goto fin;  /* [840]Q */
    for (i=1; i <= gse; i++)
      if (elldouble(nbc, XT + i*nbc2, XD + i*nbc2) > 1) goto fin;
    /* (the last iteration has initialized XG to [210*2^(gse+1)]Q) */

    if (DEBUGLEVEL >= 4)
      fprintferr("ECM: time = %6ld ms, entering B2 phase, p = %lu\n",
		 timer2(), p);

    /* inner loop over small sets of 4 curves at a time */
    for (i = nbc - 4; i >= 0; i -= 4)
    {
      if (DEBUGLEVEL >= 6)
	fprintferr("ECM: finishing curves %ld...%ld\n", i, i+3);
      /* copy relevant pointers from XH to Xh. Recall memory layout in XH is
       * nbc X coordinates followed by nbc Y coordinates for residue class
       * 1 mod 210, then the same for r.c. 11 mod 210, etc. Memory layout for
       * Xh is: four X coords for 1 mod 210, four for 11 mod 210, ..., four
       * for 209 mod 210, then the corresponding Y coordinates in the same
       * order.  This will allow us to do a giant step on Xh using just three
       * calls to elladd0() each acting on 64 points in parallel */
      for (j = 48; j--; )
      {
	k = nbc2*j + i;
	m = j << 2;		/* X coordinates */
	Xh[m]   = XH[k];   Xh[m+1] = XH[k+1];
	Xh[m+2] = XH[k+2]; Xh[m+3] = XH[k+3];
	k += nbc;		/* Y coordinates */
	Yh[m]   = XH[k];   Yh[m+1] = XH[k+1];
	Yh[m+2] = XH[k+2]; Yh[m+3] = XH[k+3];
      }
      /* build baby step table of X coords of multiples of [210]Q.  XB[4*j]
       * will point at X coords on four curves from [(j+1)*210]Q.  Until
       * we're done, we need some Y coords as well, which we keep in the
       * second half of the table, overwriting them at the end when gse==10.
       * Multiples which we already have  (by 1,2,3,4,8,16,...,2^gse) are
       * entered simply by copying the pointers, ignoring the few slots in w
       * that were initially reserved for them. Here are the initial entries */
      for (Xb=XB,k=2,j=i; k--; Xb=XB2,j+=nbc) /* do first X, then Y coords */
      {
	Xb[0]  = X[j];      Xb[1]  = X[j+1]; /* [210]Q */
	Xb[2]  = X[j+2];    Xb[3]  = X[j+3];
	Xb[4]  = XAUX[j];   Xb[5]  = XAUX[j+1]; /* [420]Q */
	Xb[6]  = XAUX[j+2]; Xb[7]  = XAUX[j+3];
	Xb[8]  = XT[j];     Xb[9]  = XT[j+1]; /* [630]Q */
	Xb[10] = XT[j+2];   Xb[11] = XT[j+3];
	Xb += 4; /* points at [420]Q */
	/* ... entries at powers of 2 times 210 .... */
	for (m = 2; m < (ulong)gse+k; m++) /* omit Y coords of [2^gse*210]Q */
	{
	  long m2 = m*nbc2 + j;
	  Xb += (2UL<<m); /* points at [2^m*210]Q */
	  Xb[0] = XAUX[m2];   Xb[1] = XAUX[m2+1];
	  Xb[2] = XAUX[m2+2]; Xb[3] = XAUX[m2+3];
	}
      }
      if (DEBUGLEVEL >= 7)
	fprintferr("\t(extracted precomputed helix / baby step entries)\n");
      /* ... glue in between, up to 16*210 ... */
      if (elladd0(12, 4,	/* 12 pts + (4 pts replicated thrice) */
		  XB + 12, XB2 + 12,
		  XB,      XB2,
		  XB + 16, XB2 + 16)
	  > 1) goto fin;	/* 4 + {1,2,3} = {5,6,7} */
      if (elladd0(28, 4,	/* 28 pts + (4 pts replicated 7fold) */
		  XB + 28, XB2 + 28,
		  XB,      XB2,
		  XB + 32, XB2 + 32)
	  > 1) goto fin;	/* 8 + {1,...,7} = {9,...,15} */
      /* ... and the remainder of the lot */
      for (m = 5; m <= (ulong)gse; m++)
      {
	/* fill in from 2^(m-1)+1 to 2^m-1 in chunks of 64 and 60 points */
	ulong m2 = 2UL << m;	/* will point at 2^(m-1)+1 */
	for (j = 0; (ulong)j < m2-64; j+=64) /* executed 0 times when m == 5 */
	{
	  if (elladd0(64, 4,
		      XB + m2 - 4, XB2 + m2 - 4,
		      XB + j,      XB2 + j,
		      XB + m2 + j,
		      (m<(ulong)gse ? XB2 + m2 + j : NULL))
	      > 1) goto fin;
	} /* j == m2-64 here, 60 points left */
	if (elladd0(60, 4,
		    XB + m2 - 4, XB2 + m2 - 4,
		    XB + j,      XB2 + j,
		    XB + m2 + j,
		    (m<(ulong)gse ? XB2 + m2 + j : NULL))
	    > 1) goto fin;
	/* when m==gse, drop Y coords of result, and when both equal 1024,
	 * overwrite Y coords of second argument with X coords of result */
      }
      if (DEBUGLEVEL >= 7) fprintferr("\t(baby step table complete)\n");
      /* initialize a few other things */
      bstp = bstp0;
      p = p0; d = d0; rcn = rcn0;
      gl = gen_1; av1 = avma;
      /* scratchspace for prod (x_i-x_j) */
      avtmp = (pari_sp)new_chunk(8 * lgefint(n));
      /* the correct entry in XB to use depends on bstp and on where we are
       * on the helix.  As we skip from prime to prime, bstp will be incre-
       * mented by snextpr() each time we wrap around through residue class
       * number 0 (1 mod 210),  but the baby step should not be taken until
       * rcn>=rcn0  (i.e. until we pass again the residue class of p0).
       * The correct signed multiplier is thus k = bstp - (rcn < rcn0),
       * and the offset from XB is four times (|k| - 1).  When k==0, we may
       * ignore the current prime  (if it had led to a factorization, this
       * would have been noted during the last giant step, or -- when we
       * first get here -- whilst initializing the helix).  When k > gss,
       * we must do a giant step and bump bstp back by -2*gss.
       * The gcd of the product of X coord differences against N is taken just
       * before we do a giant step.
       */
      /* loop over probable primes p0 < p <= nextprime(B2), inserting giant
       * steps as necessary */
      while (p < B2)
      {
	ulong p2 = p; /* save current p for diagnostics */
	/* get next probable prime */
	p = snextpr(p, &d, &rcn, &bstp, miller_k2);
	/* work out the corresponding baby-step multiplier */
	k = bstp - (rcn < rcn0 ? 1 : 0);
	/* check whether it's giant-step time */
	if (k > gss)
	{ /* take gcd */
	  gl = gcdii(gl, n);
	  if (!is_pm1(gl) && !equalii(gl, n)) { p = p2; goto fin; }
	  gl = gen_1; avma = av1;
	  while (k > gss) /* hm, just how large are those prime gaps? */
	  { /* giant step */
	    if (DEBUGLEVEL >= 7) fprintferr("\t(giant step at p = %lu)\n", p);
	    if (elladd0(64, 4,
			XG + i, YG + i,
			Xh, Yh, Xh, Yh) > 1) goto fin;
	    if (elladd0(64, 4,
			XG + i, YG + i,
			Xh + 64, Yh + 64, Xh + 64, Yh + 64) > 1) goto fin;
	    if (elladd0(64, 4,
			XG + i, YG + i,
			Xh + 128, Yh + 128, Xh + 128, Yh + 128) > 1) goto fin;
	    bstp -= (gss << 1);
	    k = bstp - (rcn < rcn0 ? 1 : 0); /* recompute multiplier */
	  }
	}
	if (!k) continue; /* point of interest is already in Xh */
	if (k < 0) k = -k;
	m = ((ulong)k - 1) << 2;
	/* accumulate product of differences of X coordinates */
	j = rcn<<2;
        avma = avtmp; /* go to garbage zone */
	gl = modii(mulii(gl, subii(XB[m],   Xh[j])), n);
	gl = modii(mulii(gl, subii(XB[m+1], Xh[j+1])), n);
	gl = modii(mulii(gl, subii(XB[m+2], Xh[j+2])), n);
	gl = mulii(gl, subii(XB[m+3], Xh[j+3]));
        avma = av1;
        gl = modii(gl, n);
      }	/* loop over p */
      avma = av1;
    } /* for i (loop over sets of 4 curves) */

    /* continuation part of main loop */
    if (dsn < dsnmax)
    {
      dsn += insist ? 1 : 2;
      if (dsn > dsnmax) dsn = dsnmax;
    }

    if (!insist && !--rep)
    {
      if (DEBUGLEVEL >= 4) {
	fprintferr("ECM: time = %6ld ms,\tellfacteur giving up.\n",
		   timer2());
	flusherr();
      }
      res = NULL; goto ret;
    }
  }
  /* END OF ECM MAIN LOOP */
fin:
  affii(gl, res);
  if (DEBUGLEVEL >= 4) {
    fprintferr("ECM: time = %6ld ms,\tp <= %6lu,\n\tfound factor = %Z\n",
	       timer2(), p, res);
    flusherr();
  }
ret:
  if (!isonstack(X)) gunclone((GEN)X);
  avma = av; return res;
}

/***********************************************************************/
/**                                                                   **/
/**                FACTORIZATION (Pollard-Brent rho) --GN1998Jun18-26 **/
/**  pollardbrent() returns a nontrivial factor of n, assuming n is   **/
/**  composite and has no small prime divisor, or NULL if going on    **/
/**  would take more time than we want to spend.  Sometimes it finds  **/
/**  more than one factor, and returns a structure suitable for       **/
/**  interpretation by ifac_crack(). (Cf Algo 8.5.2 in ACiCNT)        **/
/**                                                                   **/
/***********************************************************************/

static void
rho_dbg(long c, long msg_mask)
{
  if (c & msg_mask) return;
  fprintferr("Rho: time = %6ld ms,\t%3ld round%s\n",
             timer2(), c, (c==1?"":"s"));
  flusherr();
}

/* Tuning parameter:  for input up to 64 bits long, we must not spend more
 * than a very short time, for fear of slowing things down on average.
 * With the current tuning formula, increase our efforts somewhat at 49 bit
 * input (an extra round for each bit at first),  and go up more and more
 * rapidly after we pass 80 bits.-- Changed this to adjust for the presence of
 * squfof, which will finish input up to 59 bits quickly. */

#define tune_pb_min 14		/* even 15 seems too much. */

/* Return NULL when we run out of time, or a single t_INT containing a
 * nontrivial factor of n, or a vector of t_INTs, each triple of successive
 * entries containing a factor, an exponent (equal to one),  and a factor
 * class (NULL for unknown or zero for known composite),  matching the
 * internal representation used by the ifac_*() routines below.  Repeated
 * factors may arise;  the caller will sort the factors anyway. */
GEN
pollardbrent(GEN n)
{
  long tf = lgefint(n), size = 0, delta, retries = 0, msg_mask;
  long c0, c, k, k1, l;
  pari_sp GGG, avP, avx, av = avma;
  GEN x, x1, y, P, g, g1, res;

  if (DEBUGLEVEL >= 4) (void)timer2(); /* clear timer */

  if (tf >= 4)
    size = expi(n) + 1;
  else if (tf == 3)		/* try to keep purify happy...  */
    size = BITS_IN_LONG - bfffo((ulong)n[2]);

  if (size <= 28)
    c0 = 32;/* amounts very nearly to `insist'. Now that we have squfof(), we
             * don't insist any more when input is 2^29 ... 2^32 */
  else if (size <= 42)
    c0 = tune_pb_min;
  else if (size <= 59) /* match squfof() cutoff point */
    c0 = tune_pb_min + ((size - 42)<<1);
  else if (size <= 72)
    c0 = tune_pb_min + size - 24;
  else if (size <= 301)
    /* nonlinear increase in effort, kicking in around 80 bits */
    /* 301 gives 48121 + tune_pb_min */
    c0 = tune_pb_min + size - 60 +
      ((size-73)>>1)*((size-70)>>3)*((size-56)>>4);
  else
    c0 = 49152;	/* ECM is faster when it'd take longer */

  c = c0 << 5; /* 2^5 iterations per round */
  msg_mask = (size >= 448? 0x1fff:
                           (size >= 192? (256L<<((size-128)>>6))-1: 0xff));
PB_RETRY:
 /* trick to make a `random' choice determined by n.  Don't use x^2+0 or
  * x^2-2, ever.  Don't use x^2-3 or x^2-7 with a starting value of 2.
  * x^2+4, x^2+9 are affine conjugate to x^2+1, so don't use them either.
  *
  * (the point being that when we get called again on a composite cofactor
  * of something we've already seen, we had better avoid the same delta) */
  switch ((size + retries) & 7)
  {
    case 0:  delta=  1; break;
    case 1:  delta= -1; break;
    case 2:  delta=  3; break;
    case 3:  delta=  5; break;
    case 4:  delta= -5; break;
    case 5:  delta=  7; break;
    case 6:  delta= 11; break;
    /* case 7: */
    default: delta=-11; break;
  }
  if (DEBUGLEVEL >= 4)
  {
    if (!retries)
      fprintferr("Rho: searching small factor of %ld-bit integer\n", size);
    else
      fprintferr("Rho: restarting for remaining rounds...\n");
    fprintferr("Rho: using X^2%+1ld for up to %ld rounds of 32 iterations\n",
               delta, c >> 5);
    flusherr();
  }
  x = gen_2; P = gen_1; g1 = NULL; k = 1; l = 1;
  (void)new_chunk(10 + 6 * tf); /* enough for cgetg(10) + 3 modii */
  y = cgeti(tf); affsi(2, y);
  x1= cgeti(tf); affsi(2, x1);
  avx = avma;
  avP = (pari_sp)new_chunk(2 * tf); /* enough for x = addsi(tf+1) */
  GGG = (pari_sp)new_chunk(4 * tf); /* enough for P = modii(2tf+1, tf) */

  for (;;)			/* terminated under the control of c */
  {
    /* use the polynomial  x^2 + delta */
#define one_iter() {\
    avma = GGG; x = remii(sqri(x), n); /* to garbage zone */\
    avma = avx; x = addsi(delta,x);    /* erase garbage */\
    avma = GGG; P = mulii(P, subii(x1, x));\
    avma = avP; P = modii(P,n); }

    one_iter();

    if ((--c & 0x1f)==0)
    { /* one round complete */
      g = gcdii(n, P); if (!is_pm1(g)) goto fin;
      if (c <= 0)
      {	/* getting bored */
        if (DEBUGLEVEL >= 4)
        {
          fprintferr("Rho: time = %6ld ms,\tPollard-Brent giving up.\n",
                     timer2());
          flusherr();
        }
        avma = av; return NULL;
      }
      P = gen_1;			/* not necessary, but saves 1 mulii/round */
      if (DEBUGLEVEL >= 4) rho_dbg(c0-(c>>5), msg_mask);
      affii(x,y);
    }

    if (--k) continue;		/* normal end of loop body */

    if (c & 0x1f) /* otherwise, we already checked */
    {
      g = gcdii(n, P); if (!is_pm1(g)) goto fin;
      P = gen_1;
    }

   /* Fast forward phase, doing l inner iterations without computing gcds.
    * Check first whether it would take us beyond the alloted time.
    * Fast forward rounds count only half (although they're taking
    * more like 2/3 the time of normal rounds).  This to counteract the
    * nuisance that all c0 between 4096 and 6144 would act exactly as
    * 4096;  with the halving trick only the range 4096..5120 collapses
    * (similarly for all other powers of two)
    */
    if ((c -= (l>>1)) <= 0)
    {				/* got bored */
      if (DEBUGLEVEL >= 4)
      {
	fprintferr("Rho: time = %6ld ms,\tPollard-Brent giving up.\n",
		   timer2());
	flusherr();
      }
      avma = av; return NULL;
    }
    c &= ~0x1f;			/* keep it on multiples of 32 */

    /* Fast forward loop */
    affii(x, x1); k = l; l <<= 1;
    /* don't show this for the first several (short) fast forward phases. */
    if (DEBUGLEVEL >= 4 && (l>>7) > msg_mask)
    {
      fprintferr("Rho: fast forward phase (%ld rounds of 64)...\n", l>>7);
      flusherr();
    }
    for (k1=k; k1; k1--) one_iter();
    if (DEBUGLEVEL >= 4 && (l>>7) > msg_mask)
    {
      fprintferr("Rho: time = %6ld ms,\t%3ld rounds, back to normal mode\n",
		 timer2(), c0-(c>>5));
      flusherr();
    }

    affii(x,y);
  } /* forever */

fin:
  /* An accumulated gcd was > 1 */
  /* if it isn't n, and looks prime, return it */
  if  (!equalii(g,n))
  {
    if (miller(g,17))
    {
      if (DEBUGLEVEL >= 4)
      {
        rho_dbg(c0-(c>>5), 0);
	fprintferr("\tfound factor = %Z\n",g);
	flusherr();
      }
      avma = av; return icopy(g);
    }
    avma = avx; g1 = icopy(g);  /* known composite, keep it safe */
    avx = avma;
  }
  else g1 = n;			/* and work modulo g1 for backtracking */

  /* Here g1 is known composite */
  if (DEBUGLEVEL >= 4 && size > 192)
  {
    fprintferr("Rho: hang on a second, we got something here...\n");
    flusherr();
  }
  for(;;) /* backtrack until period recovered. Must terminate */
  {
    avma = GGG; y = remii(sqri(y), g1);
    avma = avx; y = addsi(delta,y);
    g = gcdii(subii(x1, y), g1); if (!is_pm1(g)) break;

    if (DEBUGLEVEL >= 4 && (--c & 0x1f) == 0) rho_dbg(c0-(c>>5), msg_mask);
  }

  avma = av; /* safe */
  if (g1 == n || equalii(g,g1))
  {
    if (g1 == n && equalii(g,g1))
    { /* out of luck */
      if (DEBUGLEVEL >= 4)
      {
        rho_dbg(c0-(c>>5), 0);
        fprintferr("\tPollard-Brent failed.\n"); flusherr();
      }
      if (++retries >= 4) return NULL;
      goto PB_RETRY;
    }
    /* half lucky: we've split n, but g1 equals either g or n */
    if (DEBUGLEVEL >= 4)
    {
      rho_dbg(c0-(c>>5), 0);
      fprintferr("\tfound %sfactor = %Z\n", (g1!=n ? "composite " : ""), g);
      flusherr();
    }
    res = cgetg(7, t_VEC);
    gel(res,1) = icopy(g);         /* factor */
    gel(res,2) = gen_1;		/* exponent 1 */
    gel(res,3) = (g1!=n? gen_0: NULL); /* known composite when g1!=n */

    gel(res,4) = diviiexact(n,g);       /* cofactor */
    gel(res,5) = gen_1;		/* exponent 1 */
    gel(res,6) = NULL;	/* unknown */
    return res;
  }
  /* g < g1 < n : our lucky day -- we've split g1, too */
  res = cgetg(10, t_VEC);
  /* unknown status for all three factors */
  gel(res,1) = icopy(g);              gel(res,2) = gen_1; gel(res,3) = NULL;
  gel(res,4) = diviiexact(g1,g); gel(res,5) = gen_1; gel(res,6) = NULL;
  gel(res,7) = diviiexact(n,g1); gel(res,8) = gen_1; gel(res,9) = NULL;
  if (DEBUGLEVEL >= 4)
  {
    rho_dbg(c0-(c>>5), 0);
    fprintferr("\tfound factors = %Z, %Z,\n\tand %Z\n", res[1], res[4], res[7]);
    flusherr();
  }
  return res;
}

/***********************************************************************/
/**                                                                   **/
/**              FACTORIZATION (Shanks' SQUFOF) --GN2000Sep30-Oct01   **/
/**  squfof() returns a nontrivial factor of n, assuming n is odd,    **/
/**  composite, not a pure square, and has no small prime divisor,    **/
/**  or NULL if it fails to find one.  It works on two discriminants  **/
/**  simultaneously  (n and 5n for n=1(4), 3n and 4n for n=3(4)).     **/
/**  Present implementation is limited to input <2^59, and works most **/
/**  of the time in signed arithmetic on integers <2^31 in absolute   **/
/**  size. (Cf. Algo 8.7.2 in ACiCNT)                                 **/
/**                                                                   **/
/***********************************************************************/

/* The following is invoked to walk back along the ambiguous cycle* until we
 * hit an ambiguous form and thus the desired factor, which it returns.  If it
 * fails for any reason, it returns 0.  It doesn't interfere with timing and
 * diagnostics, which it leaves to squfof().
 *
 * Before we invoke this, we've found a form (A, B, -C) with A = a^2, where a
 * isn't blacklisted and where gcd(a, B) = 1.  According to ACiCANT, we should
 * now proceed reducing the form (a, -B, -aC), but it is easy to show that the
 * first reduction step always sends this to (-aC, B, a), and the next one,
 * with q computed as usual from B and a (occupying the c position), gives a
 * reduced form, whose third member is easiest to recover by going back to D.
 * From this point onwards, we're once again working with single-word numbers.
 * No need to track signs, just work with the abs values of the coefficients. */
static long
squfof_ambig(long a, long B, long dd, GEN D)
{
  long b, c, q, qc, qcb, a0, b0, b1, c0;
  long cnt = 0; /* count reduction steps on the cycle */

  q = (dd + (B>>1)) / a;
  b = ((q*a) << 1) - B;
  {
    pari_sp av = avma;
    c = itos(divis(shifti(subii(D, sqrs(b)), -2), a));
    avma = av;
  }
#ifdef DEBUG_SQUFOF
  fprintferr("SQUFOF: ambigous cycle of discriminant %Z\n", D);
  fprintferr("SQUFOF: Form on ambigous cycle (%ld, %ld, %ld)\n", a, b, c);
#endif

  a0 = a; b0 = b1 = b;	/* end of loop detection and safeguard */

  for (;;) /* reduced cycles are finite */
  { /* reduction step */
    c0 = c;
    if (c0 > dd)
      q = 1;
    else
      q = (dd + (b>>1)) / c0;
    if (q == 1)
    {
      qcb = c0 - b; b = c0 + qcb; c = a - qcb;
    }
    else
    {
      qc = q*c0; qcb = qc - b; b = qc + qcb; c = a - q*qcb;
    }
    a = c0;

    cnt++; if (b == b1) break;

    /* safeguard against infinite loop: recognize when we've walked the entire
     * cycle in vain. (I don't think this can actually happen -- exercise.) */
    if (b == b0 && a == a0) return 0;

    b1 = b;
  }
  q = a&1 ? a : a>>1;
  if (DEBUGLEVEL >= 4)
  {
    if (q > 1)
      fprintferr("SQUFOF: found factor %ld from ambiguous form\n"
                 "\tafter %ld steps on the ambiguous cycle, time = %ld ms\n",
                 q / cgcd(q,15), cnt, timer2());
    else
      fprintferr("SQUFOF: ...found nothing on the ambiguous cycle\n"
	         "\tafter %ld steps there, time = %ld ms\n", cnt, timer2());
    if (DEBUGLEVEL >= 6) fprintferr("SQUFOF: squfof_ambig returned %ld\n", q);
  }
  return q;
}

#define SQUFOF_BLACKLIST_SZ 64

/* assume 2,3,5 do not divide n */
GEN
squfof(GEN n)
{
  ulong d1, d2;
  long tf = lgefint(n), nm4, cnt = 0;
  long a1, b1, c1, dd1, L1, a2, b2, c2, dd2, L2, a, q, c, qc, qcb;
  GEN D1, D2;
  pari_sp av = avma;
  long blacklist1[SQUFOF_BLACKLIST_SZ], blacklist2[SQUFOF_BLACKLIST_SZ];
  long blp1 = 0, blp2 = 0;
  int act1 = 1, act2 = 1;

#ifdef LONG_IS_64BIT
  if (tf > 3 || (tf == 3 && (ulong)n[2]          >= (1UL << (BITS_IN_LONG-5))))
#else  /* 32 bits */
  if (tf > 4 || (tf == 4 && (ulong)(*int_MSW(n)) >= (1UL << (BITS_IN_LONG-5))))
#endif
    return NULL; /* n too large */

  /* now we have 5 < n < 2^59 */
  nm4 = mod4(n);
  if (nm4 == 1)
  { /* n = 1 (mod4):  run one iteration on D1 = n, another on D2 = 5n */
    D1 = n;
    D2 = mulsi(5,n); d2 = itou(sqrti(D2)); dd2 = (long)((d2>>1) + (d2&1));
    b2 = (long)((d2-1) | 1);	/* b1, b2 will always stay odd */
  }
  else
  { /* n = 3 (mod4):  run one iteration on D1 = 3n, another on D2 = 4n */
    D1 = mulsi(3,n);
    D2 = shifti(n,2); dd2 = itou(sqrti(n)); d2 =  dd2 << 1;
    b2 = (long)(d2 & (~1UL)); /* largest even below d2, will stay even */
  }
  d1 = itou(sqrti(D1));
  b1 = (long)((d1-1) | 1); /* largest odd number not exceeding d1 */
  c1 = itos(shifti(subii(D1, sqru((ulong)b1)), -2));
  if (!c1) pari_err(bugparier,"squfof [caller of] (n or 3n is a square)");
  c2 = itos(shifti(subii(D2, sqru((ulong)b2)), -2));
  if (!c2) pari_err(bugparier,"squfof [caller of] (5n is a square)");
  L1 = (long)usqrtsafe(d1);
  L2 = (long)usqrtsafe(d2);
  /* dd1 used to compute floor((d1+b1)/2) as dd1+floor(b1/2), without
   * overflowing the 31bit signed integer size limit. Same for dd2. */
  dd1 = (long) ((d1>>1) + (d1&1));
  a1 = a2 = 1;

  /* The two (identity) forms (a1,b1,-c1) and (a2,b2,-c2) are now set up.
   *
   * a1 and c1 represent the absolute values of the a,c coefficients; we keep
   * track of the sign separately, via the iteration counter cnt: when cnt is
   * even, c is understood to be negative, else c is positive and a < 0.
   *
   * L1, L2 are the limits for blacklisting small leading coefficients
   * on the principal cycle, to guarantee that when we find a square form,
   * its square root will belong to an ambiguous cycle  (i.e. won't be an
   * earlier form on the principal cycle).
   *
   * When n = 3(mod 4), D2 = 12(mod 16), and b^2 is always 0 or 4 mod 16.
   * It follows that 4*a*c must be 4 or 8 mod 16, respectively, so at most
   * one of a,c can be divisible by 2 at most to the first power.  This fact
   * is used a couple of times below.
   *
   * The flags act1, act2 remain true while the respective cycle is still
   * active;  we drop them to false when we return to the identity form with-
   * out having found a square form  (or when the blacklist overflows, which
   * shouldn't happen). */
  if (DEBUGLEVEL >= 4)
  {
    fprintferr("SQUFOF: entering main loop with forms\n"
	       "\t(1, %ld, %ld) and (1, %ld, %ld)\n\tof discriminants\n"
	       "\t%Z and %Z, respectively\n", b1, -c1, b2, -c2, D1, D2);
    (void)timer2();
  }

  /* MAIN LOOP: walk around the principal cycle looking for a square form.
   * Blacklist small leading coefficients.
   *
   * The reduction operator can be computed entirely in 32-bit arithmetic:
   * Let q = floor(floor((d1+b1)/2)/c1)  (when c1>dd1, q=1, which happens
   * often enough to special-case it).  Then the new b1 = (q*c1-b1) + q*c1,
   * which does not overflow, and the new c1 = a1 - q*(q*c1-b1), which is
   * bounded by d1 in abs size since both the old and the new a1 are positive
   * and bounded by d1. */
  while (act1 || act2)
  {
    if (act1)
    { /* send first form through reduction operator if active */
      c = c1;
      q = (c > dd1)? 1: (dd1 + (b1>>1)) / c;
      if (q == 1)
      { qcb = c - b1; b1 = c + qcb; c1 = a1 - qcb; }
      else
      { qc = q*c; qcb = qc - b1; b1 = qc + qcb; c1 = a1 - q*qcb; }
      a1 = c;

      if (a1 <= L1)
      { /* blacklist this */
	if (blp1 >= SQUFOF_BLACKLIST_SZ) /* overflows: shouldn't happen */
	  act1 = 0;		/* silently */
	else
	{
	  if (DEBUGLEVEL >= 6)
	    fprintferr("SQUFOF: blacklisting a = %ld on first cycle\n", a1);
	  blacklist1[blp1++] = a1;
	}
      }
    }
    if (act2)
    { /* send second form through reduction operator if active */
      c = c2;
      q = (c > dd2)? 1: (dd2 + (b2>>1)) / c;
      if (q == 1)
      { qcb = c - b2; b2 = c + qcb; c2 = a2 - qcb; }
      else
      { qc = q*c; qcb = qc - b2; b2 = qc + qcb; c2 = a2 - q*qcb; }
      a2 = c;

      if (a2 <= L2)
      { /* blacklist this */
	if (blp2 >= SQUFOF_BLACKLIST_SZ) /* overflows: shouldn't happen */
	  act2 = 0;		/* silently */
	else
	{
	  if (DEBUGLEVEL >= 6)
	    fprintferr("SQUFOF: blacklisting a = %ld on second cycle\n", a2);
	  blacklist2[blp2++] = a2;
	}
      }
    }

    /* bump counter, loop if this is an odd iteration (i.e. if the real
     * leading coefficients are negative) */
    if (++cnt & 1) continue;

    /* second half of main loop entered only when the leading coefficients
     * are positive (i.e., during even-numbered iterations) */

    /* examine first form if active */
    if (act1 && a1 == 1) /* back to identity */
    { /* drop this discriminant */
      act1 = 0;
      if (DEBUGLEVEL >= 4)
	fprintferr("SQUFOF: first cycle exhausted after %ld iterations,\n"
		   "\tdropping it\n", cnt);
    }
    if (act1)
    {
      if (uissquarerem((ulong)a1, (ulong*)&a))
      { /* square form */
	if (DEBUGLEVEL >= 4)
	  fprintferr("SQUFOF: square form (%ld^2, %ld, %ld) on first cycle\n"
		     "\tafter %ld iterations, time = %ld ms\n",
		     a, b1, -c1, cnt, timer2());
	if (a <= L1)
	{ /* blacklisted? */
	  long j;
	  for (j = 0; j < blp1; j++)
	    if (a == blacklist1[j]) { a = 0; break; }
	}
	if (a > 0)
	{ /* not blacklisted */
	  q = cgcd(a, b1); /* imprimitive form? */
	  if (q > 1)
	  { /* q^2 divides D1 hence n [ assuming n % 3 != 0 ] */
	    avma = av;
	    if (DEBUGLEVEL >= 4) fprintferr("SQUFOF: found factor %ld^2\n", q);
	    return mkvec3(utoipos(q), gen_2, NULL);/* exponent 2, unknown status */
	  }
	  /* chase the inverse root form back along the ambiguous cycle */
	  q = squfof_ambig(a, b1, dd1, D1);
	  if (nm4 == 3 && q % 3 == 0) q /= 3;
	  if (q > 1) { avma = av; return utoipos(q); } /* SUCCESS! */
	}
	else if (DEBUGLEVEL >= 4) /* blacklisted */
	  fprintferr("SQUFOF: ...but the root form seems to be on the "
		     "principal cycle\n");
      }
    }

    /* examine second form if active */
    if (act2 && a2 == 1) /* back to identity form */
    { /* drop this discriminant */
      act2 = 0;	
      if (DEBUGLEVEL >= 4)
	fprintferr("SQUFOF: second cycle exhausted after %ld iterations,\n"
		   "\tdropping it\n", cnt);
    }
    if (act2)
    {
      if (uissquarerem((ulong)a2, (ulong*)&a))
      { /* square form */
	if (DEBUGLEVEL >= 4)
	  fprintferr("SQUFOF: square form (%ld^2, %ld, %ld) on second cycle\n"
		     "\tafter %ld iterations, time = %ld ms\n",
		     a, b2, -c2, cnt, timer2());
	if (a <= L2)
	{ /* blacklisted? */
	  long j;
	  for (j = 0; j < blp2; j++)
	    if (a == blacklist2[j]) { a = 0; break; }
	}
	if (a > 0)
	{ /* not blacklisted */
	  q = cgcd(a, b2); /* imprimitive form? */
	  /* NB if b2 is even, a is odd, so the gcd is always odd */
	  if (q > 1)
	  { /* q^2 divides D2 hence n [ assuming n % 5 != 0 ] */
	    avma = av;
	    if (DEBUGLEVEL >= 4) fprintferr("SQUFOF: found factor %ld^2\n", q);
	    return mkvec3(utoipos(q), gen_2, NULL);/* exponent 2, unknown status */
	  }
	  /* chase the inverse root form along the ambiguous cycle */
	  q = squfof_ambig(a, b2, dd2, D2);
	  if (nm4 == 1 && q % 5 == 0) q /= 5;
	  if (q > 1) { avma = av; return utoipos(q); } /* SUCCESS! */
	}
	else if (DEBUGLEVEL >= 4)	/* blacklisted */
	  fprintferr("SQUFOF: ...but the root form seems to be on the "
		     "principal cycle\n");
      }
    }
  } /* end main loop */

  /* both discriminants turned out to be useless. */
  if (DEBUGLEVEL>=4) fprintferr("SQUFOF: giving up, time = %ld ms\n", timer2());
  avma = av; return NULL;
}

/***********************************************************************/
/*                                                                     */
/*                    DETECTING ODD POWERS  --GN1998Jun28              */
/*   Factoring engines like MPQS which ultimately rely on computing    */
/*   gcd(N, x^2-y^2) to find a nontrivial factor of N can't split      */
/*   N = p^k for an odd prime p, since (Z/p^k)^* is then cyclic. Here  */
/*   is an analogue of Z_issquarerem() for 3rd, 5th and 7th powers.    */
/*   The general case is handled by is_kth_power                       */
/*                                                                     */
/***********************************************************************/

/* Multistage sieve. First stages work mod 211, 209, 61, 203 in this order
 * (first reduce mod the product of these and then take the remainder apart).
 * Second stages use 117, 31, 43, 71. Moduli which are no longer interesting
 * are skipped. Everything is encoded in a table of 106 24-bit masks. We only
 * need the first half of the residues.  Three bits per modulus indicate which
 * residues are 7th (bit 2), 5th (bit 1) or 3rd (bit 0) powers; the eight
 * moduli above are assigned right-to-left. The table was generated using: */

#if 0
L = [71, 43, 31, [O(3^2),O(13)], [O(7),O(29)], 61, [O(11),O(19)], 211];
ispow(x, N, k)=
{
  if (type(N) == "t_INT", return (ispower(Mod(x,N), k)));
  for (i = 1, #N, if (!ispower(x + N[i], k), return (0))); 1
}
check(r) =
{
  print1("  0");
  for (i=1,#L,
    N = 0;
    if (ispow(r, L[i], 3), N += 1);
    if (ispow(r, L[i], 5), N += 2);
    if (ispow(r, L[i], 7), N += 4);
    print1(N);
  ); print("ul,  /* ", r, " */")
}
for (r = 0, 105, check(r))
#endif
static ulong powersmod[106] = {
  077777777ul,  /* 0 */
  077777777ul,  /* 1 */
  013562440ul,  /* 2 */
  012402540ul,  /* 3 */
  013562440ul,  /* 4 */
  052662441ul,  /* 5 */
  016603440ul,  /* 6 */
  016463450ul,  /* 7 */
  013573551ul,  /* 8 */
  012462540ul,  /* 9 */
  012462464ul,  /* 10 */
  013462771ul,  /* 11 */
  012406473ul,  /* 12 */
  012463641ul,  /* 13 */
  052463646ul,  /* 14 */
  012503446ul,  /* 15 */
  013562440ul,  /* 16 */
  052466440ul,  /* 17 */
  012472451ul,  /* 18 */
  012462454ul,  /* 19 */
  032463550ul,  /* 20 */
  013403664ul,  /* 21 */
  013463460ul,  /* 22 */
  032562565ul,  /* 23 */
  012402540ul,  /* 24 */
  052662441ul,  /* 25 */
  032672452ul,  /* 26 */
  013573551ul,  /* 27 */
  012467541ul,  /* 28 */
  012567640ul,  /* 29 */
  032706450ul,  /* 30 */
  012762452ul,  /* 31 */
  033762662ul,  /* 32 */
  012502562ul,  /* 33 */
  032463562ul,  /* 34 */
  013563440ul,  /* 35 */
  016663440ul,  /* 36 */
  036662550ul,  /* 37 */
  012462552ul,  /* 38 */
  033502450ul,  /* 39 */
  012462643ul,  /* 40 */
  033467540ul,  /* 41 */
  017403441ul,  /* 42 */
  017463462ul,  /* 43 */
  017472460ul,  /* 44 */
  033462470ul,  /* 45 */
  052566450ul,  /* 46 */
  013562640ul,  /* 47 */
  032403640ul,  /* 48 */
  016463450ul,  /* 49 */
  016463752ul,  /* 50 */
  033402440ul,  /* 51 */
  012462540ul,  /* 52 */
  012472540ul,  /* 53 */
  053562462ul,  /* 54 */
  012463465ul,  /* 55 */
  012663470ul,  /* 56 */
  052607450ul,  /* 57 */
  012566553ul,  /* 58 */
  013466440ul,  /* 59 */
  012502741ul,  /* 60 */
  012762744ul,  /* 61 */
  012763740ul,  /* 62 */
  012763443ul,  /* 63 */
  013573551ul,  /* 64 */
  013462471ul,  /* 65 */
  052502460ul,  /* 66 */
  012662463ul,  /* 67 */
  012662451ul,  /* 68 */
  012403550ul,  /* 69 */
  073567540ul,  /* 70 */
  072463445ul,  /* 71 */
  072462740ul,  /* 72 */
  012472442ul,  /* 73 */
  012462644ul,  /* 74 */
  013406650ul,  /* 75 */
  052463471ul,  /* 76 */
  012563474ul,  /* 77 */
  013503460ul,  /* 78 */
  016462441ul,  /* 79 */
  016462440ul,  /* 80 */
  012462540ul,  /* 81 */
  013462641ul,  /* 82 */
  012463454ul,  /* 83 */
  013403550ul,  /* 84 */
  057563540ul,  /* 85 */
  017466441ul,  /* 86 */
  017606471ul,  /* 87 */
  053666573ul,  /* 88 */
  012562561ul,  /* 89 */
  013473641ul,  /* 90 */
  032573440ul,  /* 91 */
  016763440ul,  /* 92 */
  016702640ul,  /* 93 */
  033762552ul,  /* 94 */
  012562550ul,  /* 95 */
  052402451ul,  /* 96 */
  033563441ul,  /* 97 */
  012663561ul,  /* 98 */
  012677560ul,  /* 99 */
  012462464ul,  /* 100 */
  032562642ul,  /* 101 */
  013402551ul,  /* 102 */
  032462450ul,  /* 103 */
  012467445ul,  /* 104 */
  032403440ul,  /* 105 */
};

/* Returns 3, 5, or 7 if x is a cube (but not a 5th or 7th power),  a 5th
 * power (but not a 7th),  or a 7th power, and in this case creates the
 * base on the stack and assigns its address to *pt.  Otherwise returns 0.
 * x must be of type t_INT and positive;  this is not checked.  The *mask
 * argument tells us which things to check -- bit 0: 3rd, bit 1: 5th,
 * bit 2: 7th pwr;  set a bit to have the corresponding power examined --
 * and is updated appropriately for a possible follow-up call */
int
is_357_power(GEN x, GEN *pt, ulong *mask)
{
  long lx = lgefint(x), resbyte;
  ulong residue;
  pari_sp av;
  GEN y;

  *mask &= 7;		/* paranoia */
  if (!*mask) return 0;	/* useful when running in a loop */

  if (DEBUGLEVEL >= 5)
  {
    fprintferr("OddPwrs: is %Z\n\t...a", x);
    if (*mask&1) fprintferr(" 3rd%s", (*mask==7?",":(*mask!=1?" or":"")));
    if (*mask&2) fprintferr(" 5th%s", (*mask==7?", or":(*mask&4?" or":"")));
    if (*mask&4) fprintferr(" 7th");
    fprintferr(" power?\n\tmodulo: resid. (remaining possibilities)\n");
  }
  residue = (lx == 3)? x[2]: umodiu(x, 211*209*61*203);

#define check_res(N, shift) {\
  resbyte = residue%N; if ((ulong)resbyte > (N>>1)) resbyte = N - resbyte;\
  *mask &= (powersmod[resbyte] >> shift); \
  if (DEBUGLEVEL >= 5)\
    fprintferr("\t   %3ld:  %3ld   (3rd %ld, 5th %ld, 7th %ld)\n",\
               N, resbyte, *mask&1, (*mask>>1)&1, (*mask>>2)&1);\
  if (!*mask) return 0;\
}
  check_res(211, 0);
  if (*mask & 3) check_res(209UL, 3);
  if (*mask & 3) check_res( 61UL, 6);
  if (*mask & 5) check_res(203UL, 9);
  residue = (lx == 3)? x[2]: umodiu(x, 117*31*43*71);
  if (*mask & 1) check_res(117UL,12);
  if (*mask & 3) check_res( 31UL,15);
  if (*mask & 5) check_res( 43UL,18);
  if (*mask & 6) check_res( 71UL,21);

  av = avma;
  while (*mask)
  {
    long e, b;
    /* priority to higher powers: if we have a 21st, it is easier to rediscover
     * that its 7th root is a cube than that its cube root is a 7th power */
         if (*mask & 4) { b = 4; e = 7; }
    else if (*mask & 2) { b = 2; e = 5; }
    else                { b = 1; e = 3; }
    y = mpround( sqrtnr(itor(x, DEFAULTPREC + (lx-2) / e), e) );
    if (equalii(powiu(y,e), x))
    {
      if (!pt) { avma = av; return e; }
      avma = (pari_sp)y; *pt = gerepileuptoint(av, y);
      return e;
    }
    if (DEBUGLEVEL >= 5)
      fprintferr("\tBut it nevertheless wasn't a %ld%s power.\n", e,eng_ord(e));
    *mask &= ~b; /* turn the bit off */
    avma = av;
  }
  return 0;
}

/* p not necessarily prime */
ulong
is_kth_power(GEN x, ulong p, GEN *pt, byteptr d)
{
  int init = 0;
  long j, k;
  ulong q, prkmodq, residue, elt;
  GEN y;
  byteptr d0;
  pari_sp av = avma;

  if (d)
  {
    q = p; d0 = d;
  }
  else
  {
    q = 0; d0 = diffptr;
    maxprime_check(p);
    while (q < p) NEXT_PRIME_VIADIFF(q,d0);
  }
  /* for modular checks, use small primes q congruent 1 mod curexp */
  /* #checks is tunable, for small p we can afford to do more than 5 */
  for (j = (p<40 ? 7 : p<80 ? 5 : p<250 ? 4 : 3); j > 0; j--)
  {
    do
    {
      if (*d0) NEXT_PRIME_VIADIFF(q,d0);
      else {
        if (init) q += p; else { init = 1; q += (p + 1 - q % p); }
        while (!uisprime(q)) { q += p; }
        break;
      }
    } while (q % p != 1);
    if (DEBUGLEVEL>4) fprintferr("\tchecking modulo %ld\n", q);
    /* XXX give up if q is too large, huh? */
    residue = umodiu(x, q);
    if (residue == 0) continue;
    /* find a generator of the subgroup of index curexp in (Z/qZ)^* */
    prkmodq = elt = Fl_pow(gener_Fl(q), p, q);
    /* see whether our residue is in the subgroup */
    for (k = (q - 1)/p; k > 0; k--)
    {
      if (elt == residue) break;
      elt = Fl_mul(elt, prkmodq, q);
    }
    if (!k) /* not found */
    {
      if (DEBUGLEVEL>5) fprintferr("\t- ruled out\n");
      avma = av; return 0;
    }
  }
  avma = av;

  if (DEBUGLEVEL>4) fprintferr("OddPwrs: passed modular checks\n");
  /* go to the horse's mouth... */
  y = mpround( sqrtnr(itor(x, nbits2prec((expi(x)+16*p)/p)), p) );
  if (!equalii(powiu(y, p), x)) {
    if (DEBUGLEVEL>4) fprintferr("\tBut it wasn't a pure power.\n");
    avma = av; return 0;
  }
  if (!pt) avma = av; else { avma = (pari_sp)y; *pt = gerepileuptoint(av, y); }
  return 1;
}

/* generic version for exponents 11 or larger.  Cut off when x^(1/k) fits
 * into (say) 14 bits, since we would have found it by trial division.
 * Interface is similar to is_357_power(), but instead of the mask, we keep
 * the current test exponent around.  Everything needed here (primitive roots
 * etc.) is computed from scratch on the fly; compared to the size of numbers
 * under consideration, these word-sized computations take negligible time.
 * Experimentally making the cutoff point caller-configurable... */
int
is_odd_power(GEN x, GEN *pt, ulong *curexp, ulong cutoffbits)
{
  long size = expi(x) /* not +1 */;
  ulong p = 0;
  pari_sp av = avma;
  byteptr d = diffptr;

  /* cutting off at 1 is safe, but direct root extraction attempts will be
   * faster when trial division has been used to discover very small
   * bases.  We become competitive at about cutoffbits = 4 */
  if (!cutoffbits) cutoffbits = 1;
  /* prepare for iterating curexp over primes */
  if (*curexp < 11) *curexp = 11;
  while (p < *curexp) { NEXT_PRIME_VIADIFF(p,d); if (!*d) break; }
  while (p < *curexp) {  p = itou( nextprime(utoipos(p + 1)) ); }
  *curexp = p;

  if (DEBUGLEVEL>4) fprintferr("OddPwrs: examining %Z\n", x);
  /* check size of x vs. curexp */
  /* tunable cutoff was 18 initially, but cheap enough to go further (given
   * the 2^14 minimal cutoff point for trial division) */
  while (size/p >= cutoffbits /*tunable*/ ) {
    if (DEBUGLEVEL>4) fprintferr("OddPwrs: testing for exponent %ld\n", p);
    /* if found, caller should call us again without changing *curexp */
    if (is_kth_power(x, p, pt, d)) return p;
    NEXT_PRIME_VIADIFF(p,d);
    *curexp = p;
  }
  avma = av; return 0; /* give up */
}

/***********************************************************************/
/**                                                                   **/
/**                FACTORIZATION  (master iteration)                  **/
/**      Driver for the various methods of finding large factors      **/
/**      (after trial division has cast out the very small ones).     **/
/**                        GN1998Jun24--30                            **/
/**                                                                   **/
/***********************************************************************/

/*  Direct use:
 *  ifac_start()  registers a number (without prime factors < 100) with the
 *    iterative factorizer, and also registers whether or not we should
 *    terminate early if we find a square factor, and a hint about which
 *    method(s) to use. This must always be called first. If input is not
 *    composite, we will have an oo loop later. The routine decomposes it
 *    nontrivially into a product of two factors except in squarefreeness
 *    (`Moebius') mode.
 *  ifac_primary_factor()  returns a prime divisor (not necessarily the
 *    smallest) and the corresponding exponent. */

/*  Encapsulated user interface:
 *  ifac_decomp()  does the right thing for auxdecomp()  (put a succession of
 *  prime divisor / exponent pairs onto the stack, not necessarily sorted.
 *
 *  For each of the arithmetic functions, there is a `contributor' ifac_xxx,
 *  to be called on any large composite cofactor left over after trial division
 *  by small primes, whose result can then be added to or multiplied with
 *  whatever we already have: xxx is one of moebius, issquarefree,  totient,
 *  omega, bigomega,  numdiv,  sumdiv, sumdivk, */

/* We never test whether the input number is prime or composite, since
 * presumably it will have come out of the small factors finder stage
 * (which doesn't really exist yet but which will test the left-over
 * cofactor for primality once it does).
 */
/* The data structure in which we preserve whatever we know at any given
 * time about our number N is kept on the PARI stack, and updated as needed.
 * This makes the machinery re-entrant, and avoids memory leaks when a lengthy
 * factorization is interrupted. We also make an effort to keep the whole
 * affair connected, and the parent object will always be older than its
 * children.  This may in rare cases lead to some extra copying around, and
 * knowing what is garbage at any given time is not entirely trivial. See below
 * for examples how to do it right.  (Connectedness is destroyed if callers of
 * ifac_main() create stuff on the stack in between calls. This is harmless
 * as long as ifac_realloc() is used to re-create a connected object at
 * the head of the stack just before collecting garbage.)
 */
/* Note that a PARI integer can have hundreds of millions of distinct prime
 * factors larger than 2^16, given enough memory.  And since there's no
 * guarantee that we will find factors in order of increasing size, we must
 * be prepared to drag a very large amount of data around  (although this
 * will _very_ rarely happen for random input!).  We start with a small
 * structure and extend it when necessary.
 */
/* The idea of data structure and algorithm is:
 * Let N0 be whatever is currently left of N after dividing off all the
 * prime powers we have already returned to the caller.  Then we maintain
 * N0 as a product
 * (1)   N0 = \prod_i P_i^{e_i} * \prod_j Q_j^{f_j} * \prod_k C_k^{g_k}
 * where the P_i and Q_j are distinct primes, each C_k is known composite,
 * none of the P_i divides any C_k, and we also know the total ordering
 * of all the P_i, Q_j and C_k  (in particular, we will never try to divide
 * a C_k by a larger Q_j).  Some of the C_k may have common factors, although
 * this will not often be the case.
 */
/* Caveat implementor:  Taking gcds among C_k's is very likely to cost at
 * least as much time as dividing off any primes as we find them, and book-
 * keeping would be a nightmare  (since D=gcd(C_1,C_2) can still have common
 * factors with both C_1/D and C_2/D, and so on...).
 */
/* At startup, we just initialize the structure to
 * (2)        N = C_1^1   (composite).
 */
/* Whenever ifac_primary_factor() or ifac_decomp()  (or, mutatis mutandis,
 * one of the arithmetic user interface routines) needs a primary factor, and
 * the smallest thing in our list is P_1, we return that and its exponent, and
 * remove it from our list. (When nothing is left, we return a sentinel value
 * -- gen_1.  And in Moebius mode, when we see something with exponent > 1,
 * whether prime or composite,  we return gen_0 or 0, depending on the
 * function). In all other cases, ifac_main() iterates the following steps
 * until we have a P_1 in the smallest position.
 */
/* When the smallest item is C_1  (as it is initially):
 * (3.1) Crack C_1 into a nontrivial product  U_1 * U_2  by whatever method
 * comes to mind for this size. (U for `unknown'.)  Cracking will detect
 * squares  (and biquadrates etc), and it may detect odd powers, so we may
 * instead see a power of some U_1 here, or even something of the form
 * U_1^k*U_2^k. (Of course the exponent already attached to C_1 is taken into
 * account in the following.)
 * (3.2) If we have U_1*U_2, sort the two factors; convert to U_1^2 if they
 * happen to be equal  (which they shouldn't -- squares should have been
 * caught in stage 3.1). Note that U_1 and (if it exists) U_2 are smaller than
 * anything else in our list.
 * (3.3) Check U_1 (and U_2) for primality, and flag them accordingly.
 * (3.4) Iterate.
 */
/* When the smallest item is Q_1:
 * This is the potentially unpleasant case.  The idea is to go through the
 * entire list and try to divide Q_1 off each of the current C_k's, which
 * will usually fail, but may succeed several times.  When a division was
 * successful, the corresponding C_k is removed from our list, and the co-
 * factor becomes a U_l for the moment unless it is 1 (which happens when C_k
 * was a power of Q_1).  When we're through we upgrade Q_1 to P_1 status,
 * then do a primality check on each U_l and sort it back into the list
 * either as a Q_j or as a C_k.  If during the insertion sort we discover
 * that some U_l equals some P_i or Q_j or C_k we already have, we just add
 * U_l's exponent to that of its twin.  (The sorting should therefore happen
 * before the primality test).
 * Note that this may produce one or more elements smaller than the P_1
 * we just confirmed, so we may have to repeat the iteration.
 */
/* There's a little trick that avoids some Q_1 instances.  Just after we do
 * a sweep to classify all current unknowns as either composites or primes,
 * we do another downward sweep beginning with the largest current factor
 * and stopping just above the largest current composite.  Every Q_j we
 * pass is turned into a P_i.  (Different primes are automatically coprime
 * among each other, and primes tend not to divide smaller composites.)
 */
/* (We have no use for comparing the square of a prime to N0.  Normally
 * we will get called after casting out only the smallest primes, and
 * since we cannot guarantee that we see the large prime factors in as-
 * cending order, we cannot stop when we find one larger than sqrt(N0).)
 */
/* Data structure:  We keep everything in a single t_VEC of t_INTs.  The
 * first component records whether we're doing full (NULL) or Moebius (one)
 * factorization;  in the latter case many subroutines return a sentinel
 * value as soon as they spot an exponent > 1.  The second component records
 * the hint from factorint()'s optional flag, for use by ifac_crack().
 * The remaining components  (initially 15)  are used in groups of three:
 * a GEN pointer at the t_INT value of the factor, a pointer at the t_INT
 * exponent (usually gen_1 or gen_2 so we don't clutter up the stack too
 * much),  and another t_INT GEN pointer to record the class of the factor:
 * NULL for unknown, zero for known composite C_k, one for known prime Q_j
 * awaiting trial division, and two for finished prime P_i.
 */
/* When during the division stage we re-sort a C_k-turned-U_l to a lower
 * position, we rotate any intervening material upward towards its old
 * slot.  When a C_k was divided down to 1, its slot is left empty at
 * first;  similarly when the re-sorting detects a repeated factor.
 * After the sorting phase, we de-fragment the list and squeeze all the
 * occupied slots together to the high end, so that ifac_crack() has room
 * for new factors.  When this doesn't suffice, we abandon the current
 * vector and allocate a somewhat larger one, defragmenting again during
 * copying.
 */
/* (For internal use, note that all exponents will fit into C longs, given
 * PARI's lgefint field size.  When we work with them, we sometimes read
 * out the GEN pointer, and sometimes do an itos, whatever is more con-
 * venient for the task at hand.)
 */

/*** Overview and forward declarations: ***/

/* The `*where' argument in the following points into *partial at the first of
 * the three fields of the first occupied slot.  It's there because the caller
 * would already know where `here' is, so we don't want to search for it again.
 * We do not preserve this from one user-interface call to the next. */

static GEN ifac_find(GEN *partial, GEN *where);
/* Return GEN pointing at the first nonempty slot strictly behind the current
 * *where, or NULL if such doesn't exist.  Can be used to skip a range of
 * vacant slots, or to initialize *where in the first place (pass partial in
 * both args).  Does not modify its argument pointers. */

void ifac_realloc(GEN *partial, GEN *where, long new_lg);
/* Move to a larger main vector, updating *where if it points into it, and
 * *partial in any case. Can be used as a specialized gcopy before
 * a gerepileupto() (pass 0 as the new length). Normally, one would pass
 * new_lg=1 to let this function guess the new size.  To be used sparingly. */

static long ifac_crack(GEN *partial, GEN *where);
/* Split the first (composite) entry.  There _must_ already be room for another
 * factor below *where, and *where will be updated.  Factor and cofactor are
 * inserted in the correct order, updating *where, or factor^k will be inserted
 * if such should be the case  (leaving *where unchanged). The factor or
 * factors will be set to unknown, and inherit the exponent  (or a multiple
 * thereof) of its/their ancestor.  Returns number of factors written into the
 * structure  (normally 2, but 1 if a factor equalled its cofactor, and may be
 * more than 1 if a factoring engine returned a vector of factors instead of a
 * single factor).  Can reallocate the data structure in the vector-of-factors
 * case  (but not in the most common single-factor case) */

static long ifac_insert_multiplet(GEN *partial, GEN *where, GEN facvec);
/* Gets called to complete ifac_crack()'s job when a factoring engine splits
 * the current factor into a product of three or more new factors. Makes room
 * for them if necessary, sorts them, gives them the right exponents and class.
 * Also returns the number of factors actually written, which may be less than
 * the number of components in facvec if there are duplicates.--- Vectors of
 * factors  (cf pollardbrent()) actually contain `slots' of three GENs per
 * factor with the three fields interpreted as in our partial factorization
 * data structure.  Thus `engines' can tell us what they already happen to
 * know about factors being prime or composite and/or appearing to a power
 * larger than the first */

static long ifac_divide(GEN *partial, GEN *where);
/* Divide all current composites by first (prime, class Q) entry, updating its
 * exponent, and turning it into a finished prime (class P).  Return 1 if any
 * such divisions succeeded  (in Moebius mode, the update may then not have
 * been completed), or 0 if none of them succeeded.  Doesn't modify *where. */

static long ifac_sort_one(GEN *partial, GEN *where, GEN washere);
/* re-sort one (typically unknown) entry from washere to a new position,
 * rotating intervening entries upward to fill the vacant space. It may happen
 * that the new position is the same as the old one, or that the new value of
 * the entry coincides with a value already occupying a lower slot, in which
 * latter case we just add exponents  (and use the `more known' class, and
 * return 1 immediately when in Moebius mode). The slots between *where and
 * washere must be in sorted order, so a sweep using this to re-sort several
 * unknowns must proceed upward  (see ifac_resort()).  Return 1 if we see an
 * exponent > 1  (in Moebius mode without completing the update), 0 otherwise.
 */

static long ifac_resort(GEN *partial, GEN *where);
/* sort all current unknowns downward to where they belong.  Sweeps in the
 * upward direction.  Not needed after ifac_crack(), only when ifac_divide()
 * returned true.  May update *where.  Returns 1 when an ifac_sort_one() call
 * does so to indicate a repeated factor, or 0 if all such calls returned 0 */

static void ifac_defrag(GEN *partial, GEN *where);
/* defragment: collect and squeeze out any unoccupied slots above *where
 * during a downward sweep.  Unoccupied slots arise when a composite factor
 * dissolves completely whilst dividing off a prime, or when ifac_resort()
 * spots a coincidence and merges two factors.  *where will be updated */

static void ifac_whoiswho(GEN *partial, GEN *where, long after_crack);
/* determine primality or compositeness of all current unknowns, and set
 * class Q primes to finished (class P) if everything larger is already
 * known to be prime.  When after_crack is nonnegative, only look at the
 * first after_crack things in the list (do nothing when it's zero) */

static GEN ifac_main(GEN *partial);
/* main loop:  iterate until smallest entry is a finished prime;  returns
 * a `where' pointer, or gen_1 if nothing left, or gen_0 in Moebius mode if
 * we aren't squarefree */

/* In the most common cases, control flows from the user interface to
 * ifac_main() and then to a succession of ifac_crack()s and ifac_divide()s,
 * with (typically) none of the latter finding anything. */

/** user interface: **/
/* return initial data structure, see ifac_crack() for the hint argument */
GEN ifac_start(GEN n, long moebius, long hint);

/* run main loop until primary factor is found, return the prime and assign the
 * exponent. If nothing left, return gen_1 and set exponent to 0; if in Moebius
 * mode and a square factor is discovered, return gen_0 and set exponent to 0 */
GEN ifac_primary_factor(GEN *partial, long *exponent);

/* call ifac_start() and run main loop until factorization is complete,
 * accumulating prime / exponent pairs on the PARI stack to be picked up
 * by aux_end().  Return number of distinct primes found */
long ifac_decomp(GEN n, long hint);

/* encapsulated functions; these call ifac_start() themselves, ensure stack
 * housekeeping etc. Call them on any large composite left over after trial
 * division, and multiply/add the result onto whatever you already have from
 * the small factors.  On large primes, they will run into trouble */
long ifac_moebius(GEN n, long hint);
long ifac_issquarefree(GEN n, long hint);
long ifac_omega(GEN n, long hint);
long ifac_bigomega(GEN n, long hint);
GEN ifac_totient(GEN n, long hint); /* for gp's eulerphi() */
GEN ifac_numdiv(GEN n, long hint);
GEN ifac_sumdiv(GEN n, long hint);
GEN ifac_sumdivk(GEN n, long k, long hint);

/*** implementation ***/

#define ifac_initial_length 24	/* codeword, moebius flag, hint, 7 slots */
/* (more than enough in most cases -- a 512-bit product of distinct 8-bit
 * primes needs at most 7 slots at a time) */

GEN
ifac_start(GEN n, long moebius, long hint)
{
  GEN part, here;

  if (typ(n) != t_INT) pari_err(typeer, "ifac_start");
  if (!signe(n)) pari_err(talker, "factoring 0 in ifac_start");

  part = cgetg(ifac_initial_length, t_VEC);
  here = part + ifac_initial_length;
  gel(part,1) = moebius? gen_1 : NULL;
  gel(part,2) = stoi(hint);
  if (isonstack(n)) n = absi(n);
  /* make copy, because we'll later want to replace it in place.
   * If it's not on stack, then we assume it is a clone made for us by
   * auxdecomp1, and we assume the sign has already been set positive */
  /* fill first slot at the top end */
  gel(--here,0) = gen_0;/* initially composite */
  gel(--here,0) = gen_1;/* initial exponent 1 */
  gel(--here,0) = n;
  while (here > part + 3) gel(--here,0) = NULL;
  return part;
}

static GEN
ifac_find(GEN *partial, GEN *where)
{
  long lgp = lg(*partial);
  GEN end = *partial + lgp;
  GEN scan = *where + 3;

#ifdef IFAC_DEBUG
  if (!*partial || typ(*partial) != t_VEC) pari_err(typeer, "ifac_find");
  if (lg(*partial) < ifac_initial_length)
    pari_err(talker, "partial impossibly short in ifac_find");
  if (!(*where) || *where > *partial + lgp - 3 || *where < *partial)
    pari_err(talker, "`*where\' out of bounds in ifac_find");
#endif
  while (scan < end && !*scan) scan += 3;
  /* paranoia -- check completely NULLed ? nope -- we never inspect the
   * exponent field for deciding whether a slot is empty or occupied */
  if (scan < end)
  {
    if (DEBUGLEVEL >= 5 && !scan[1])
      pari_err(talker, "factor has NULL exponent in ifac_find");
    return scan;
  }
  return NULL;
}

/* simple defragmenter */
static void
ifac_defrag(GEN *partial, GEN *where)
{
  long lgp = lg(*partial);
  GEN scan_new = *partial + lgp - 3, scan_old = scan_new;

  while (scan_old >= *where)
  {
    if (*scan_old)
    { /* slot occupied */
      if (scan_old < scan_new)
      {
	scan_new[2] = scan_old[2];
	scan_new[1] = scan_old[1];
	*scan_new = *scan_old;
      }
      scan_new -= 3; /* point at next slot to be written */
    }
    scan_old -= 3;
  }
  scan_new += 3; /* back up to last slot written */
  *where = scan_new;
  while (scan_new > *partial + 3) gel(--scan_new,0) = NULL; /* erase junk */
}

/* and complex version combined with reallocation.  If new_lg is 0, use the old
 * length, so this acts just like gcopy except that the where pointer is
 * carried along; if it is 1, we make an educated guess. Exception:  If new_lg
 * is 0, the vector is full to the brim, and the first entry is composite, we
 * make it longer to avoid being called again a microsecond later. It is safe
 * to call this with NULL for the where argument;  if it doesn't point anywhere
 * within the old structure, it is left alone */
void
ifac_realloc(GEN *partial, GEN *where, long new_lg)
{
  long old_lg = lg(*partial);
  GEN newpart, scan_new, scan_old;

  if (new_lg == 1)
    new_lg = 2*old_lg - 6;	/* from 7 slots to 13 to 25... */
  else if (new_lg <= old_lg)	/* includes case new_lg == 0 */
  {
    new_lg = old_lg;
    if ((*partial)[3] &&	/* structure full */
	(gel(*partial,5) == gen_0 || gel(*partial,5)==NULL))
				/* and first entry composite or unknown */
      new_lg += 6; /* give it a little more breathing space */
  }
  newpart = cgetg(new_lg, t_VEC);
  if (DEBUGMEM >= 3)
    fprintferr("IFAC: new partial factorization structure (%ld slots)\n",
	       (new_lg - 3)/3);
  newpart[1] = (*partial)[1];	/* moebius */
  icopyifstack((*partial)[2], newpart[2]); /* hint */
  /* downward sweep through the old *partial, picking up where1 and carrying it
   * over if and when we pass it.  (This will only be useful if it pointed at a
   * non-empty slot.)  Factors are icopy()d so that we again have a nice object
   * (parent older than children, connected), except the one factor that may
   * still be living in a clone where n originally was; exponents are similarly
   * copied if they aren't global constants; class-of-factor fields are always
   * global constants so we need only copy them as pointers. Caller may then do
   * a gerepileupto() or a gerepilemanysp() */
  scan_new = newpart + new_lg - 3;
  scan_old = *partial + old_lg - 3;
  for (; scan_old > *partial + 2; scan_old -= 3)
  {
    if (*where == scan_old) *where = scan_new;
    if (!*scan_old) continue; /* skip empty slots */

    icopyifstack(*scan_old, *scan_new);
    icopyifstack(scan_old[1], scan_new[1]);
    scan_new[2] = scan_old[2]; scan_new -= 3;
  }
  scan_new += 3; /* back up to last slot written */
  while (scan_new > newpart + 3) gel(--scan_new,0) = NULL;
  *partial = newpart;
}

#define moebius_mode ((*partial)[1])

/* Bubble-sort-of-thing sort.  Won't be exercised frequently, so this is ok */
static long
ifac_sort_one(GEN *partial, GEN *where, GEN washere)
{
  GEN scan = washere - 3;
  GEN value, exponent, class0, class1;
  long cmp_res;

  if (DEBUGLEVEL >= 5)		/* none of these should ever happen */
  {
    long lgp;
    if (!*partial || typ(*partial) != t_VEC)
      pari_err(typeer, "ifac_sort_one");
    if ((lgp = lg(*partial)) < ifac_initial_length)
      pari_err(talker, "partial impossibly short in ifac_sort_one");
    if (!(*where) || *where < *partial + 3 || *where > *partial + lgp - 3)
      pari_err(talker, "`*where\' out of bounds in ifac_sort_one");
    if (!washere || washere < *where || washere > *partial + lgp - 3)
      pari_err(talker, "`washere\' out of bounds in ifac_sort_one");
  }
  value    = gel(washere,0);
  exponent = gel(washere,1);
  if (exponent != gen_1 && moebius_mode && cmpui(1,exponent) < 0)
    return 1;			/* should have been detected by caller */
  class0 = gel(washere,2);

  if (scan < *where) return 0;	/* nothing to do, washere==*where */

  cmp_res = -1;			/* sentinel */
  while (scan >= *where)	/* therefore at least once */
  {
    if (*scan)
    { /* current slot nonempty, check against where */
      cmp_res = cmpii(value, gel(scan,0));
      if (cmp_res >= 0) break;	/* have found where to stop */
    }
    /* copy current slot upward by one position and move pointers down */
    scan[5] = scan[2];
    scan[4] = scan[1];
    scan[3] = *scan;
    scan -= 3;
  }
  scan += 3;
  /* at this point there are the following possibilities:
   * (*) cmp_res == -1.  Either value is less than that at *where, or for some
   * reason *where was pointing at one or more vacant slots and any factors we
   * saw en route were larger than value.  At any rate, scan == *where now, and
   * scan is pointing at an empty slot, into which we'll stash our entry.
   * (*) cmp_res == 0.  The entry at scan-3 is the one, we compare class0
   * fields and add exponents, and put it all into the vacated scan slot,
   * NULLing the one at scan-3  (and possibly updating *where).
   * (*) cmp_res == 1.  The slot at scan is the one to store our entry into. */
  if (cmp_res)
  {
    if (cmp_res < 0 && scan != *where)
      pari_err(talker, "misaligned partial detected in ifac_sort_one");
    gel(scan,0) = value;
    gel(scan,1) = exponent;
    gel(scan,2) = class0; return 0;
  }
  /* case cmp_res == 0: repeated factor detected */
  if (DEBUGLEVEL >= 4)
    fprintferr("IFAC: repeated factor %Z\n\tdetected in ifac_sort_one\n",
	       value);
  if (moebius_mode) return 1;	/* not squarefree */
  /* if old class0 was composite and new is prime, or vice versa, complain
   * (and if one class0 was unknown and the other wasn't, use the known one) */
  class1 = gel(scan,-1);
  if (class0) /* should never be used */
  {
    if (class1)
    {
      if (class0 == gen_0 && class1 != gen_0)
	pari_err(talker, "composite equals prime in ifac_sort_one");
      else if (class0 != gen_0 && class1 == gen_0)
	pari_err(talker, "prime equals composite in ifac_sort_one");
      else if (class0 == gen_2)	/* should happen even less */
	gel(scan,2) = class0;	/* use it */
    }
    else			/* shouldn't happen either */
      gel(scan,2) = class0;	/* use it */
  }
  /* else stay with the existing known class0 */
  gel(scan,2) = class1;
  /* in any case, add exponents */
  if (gel(scan,-2) == gen_1 && exponent == gen_1)
    gel(scan,1) = gen_2;
  else
    gel(scan,1) = addii(gel(scan,-2), exponent);
  /* move the value over and null out the vacated slot below */
  *scan = scan[-3];
  gel(--scan,0) = NULL;
  gel(--scan,0) = NULL;
  gel(--scan,0) = NULL;
  /* finally, see whether *where should be pulled in */
  if (scan == *where) *where += 3;
  return 0;
}

/* the following loop around the former doesn't need to check moebius_mode
 * because ifac_sort_one() never returns 1 in normal mode */
static long
ifac_resort(GEN *partial, GEN *where)
{
  long lgp = lg(*partial), res;
  GEN scan = *where;

  for (; scan < *partial + lgp; scan += 3)
    if (*scan && !scan[2] /* slot occupied with an unknown */
        && (res = ifac_sort_one(partial, where, scan)) ) return res;
  return 0;
}

/* sweep downward so we can with luck turn some Qs into Ps */
static void
ifac_whoiswho(GEN *partial, GEN *where, long after_crack)
{
  long lgp = lg(*partial), larger_compos = 0;
  GEN scan, scan_end = *partial + lgp - 3;

#ifdef IFAC_DEBUG
  if (!*partial || typ(*partial) != t_VEC)
    pari_err(typeer, "ifac_whoiswho");
  if (lg(*partial) < ifac_initial_length)
    pari_err(talker, "partial impossibly short in ifac_whoiswho");
  if (!(*where) || *where > scan_end || *where < *partial + 3)
    pari_err(talker, "`*where\' out of bounds in ifac_whoiswho");
#endif
  if (after_crack == 0) return;
  if (after_crack > 0)
  {
    larger_compos = 1;		/* disable Q-to-P trick */
    scan = *where + 3*(after_crack - 1);
				/* check at most after_crack entries */
    if (scan > scan_end)	/* ooops... */
    {
      pari_warn(warner, "avoiding nonexistent factors in ifac_whoiswho");
      scan = scan_end;
    }
  }
  else { larger_compos = 0; scan = scan_end; }

  for (; scan >= *where; scan -= 3)
  {
    if (scan[2])
    { /* known class of factor */
      if (gel(scan,2) == gen_0) larger_compos = 1;
      else if (!larger_compos && gel(scan,2) == gen_1)
      {
	if (DEBUGLEVEL >= 3)
	{
	  fprintferr("IFAC: factor %Z\n\tis prime (no larger composite)\n",
		     gel(*where,0));
	  fprintferr("IFAC: prime %Z\n\tappears with exponent = %ld\n",
		     gel(*where,0), itos(gel(*where,1)));
	}
	gel(scan,2) = gen_2;
      }
      continue;
    }
    gel(scan,2) = BSW_psp(gel(scan,0)) ?
       (larger_compos ? gen_1 : gen_2) : /* one- or finished prime */
       gen_0;			     /* composite */

    if (gel(scan,2) == gen_0) larger_compos = 1;
    if (DEBUGLEVEL >= 3)
      fprintferr("IFAC: factor %Z\n\tis %s\n", *scan,
		 (gel(scan,2) == gen_0 ? "composite" : "prime"));
  }
}

/* Here we normally do not check that the first entry is a not-finished
 * prime.  Stack management: we may allocate a new exponent */
static long
ifac_divide(GEN *partial, GEN *where)
{
  long lgp = lg(*partial);
  GEN scan = *where + 3;
  long res = 0, exponent, newexp, otherexp;

#ifdef IFAC_DEBUG
  if (!*partial || typ(*partial) != t_VEC)
    pari_err(typeer, "ifac_divide");
  if (lg(*partial) < ifac_initial_length)
    pari_err(talker, "partial impossibly short in ifac_divide");
  if (!(*where) || *where > *partial + lgp - 3 || *where < *partial + 3)
    pari_err(talker, "`*where\' out of bounds in ifac_divide");
  if (gel(*where,2) != gen_1)
    pari_err(talker, "division by composite or finished prime in ifac_divide");
  if (!(**where))
    pari_err(talker, "division by nothing in ifac_divide");
#endif
  newexp = exponent = itos(gel(*where,1));
  if (exponent > 1 && moebius_mode) return 1;
  /* should've been caught by caller */

  for (; scan < *partial + lgp; scan += 3)
  {
    if (gel(scan,2) != gen_0) continue; /* the other thing ain't composite */
    otherexp = 0;
    /* divide in place to keep stack clutter minimal */
    while (dvdiiz(gel(scan,0), gel(*where,0), gel(scan,0)))
    {
      if (moebius_mode) return 1; /* immediately */
      if (!otherexp) otherexp = itos(gel(scan,1));
      newexp += otherexp;
    }
    if (newexp > exponent)	/* did anything happen? */
    {
      gel(*where,1) = (newexp == 2 ? gen_2 : utoipos(newexp));
      exponent = newexp;
      if (is_pm1(*scan)) /* factor dissolved completely */
      {
	gel(scan,0) = gel(scan,1) = NULL;
	if (DEBUGLEVEL >= 4)
	  fprintferr("IFAC: a factor was a power of another prime factor\n");
      }
      else if (DEBUGLEVEL >= 4)
	fprintferr("IFAC: a factor was divisible by another prime factor,\n"
	           "\tleaving a cofactor = %Z\n", *scan);
      gel(scan,2) = NULL;	/* at any rate it's Unknown now */
      res = 1;
      if (DEBUGLEVEL >= 5)
	fprintferr("IFAC: prime %Z\n\tappears at least to the power %ld\n",
		   gel(*where,0), newexp);
    }
  } /* for */
  gel(*where,2) = gen_2; /* make it a finished prime */
  if (DEBUGLEVEL >= 3)
    fprintferr("IFAC: prime %Z\n\tappears with exponent = %ld\n",
	       gel(*where,0), newexp);
  return res;
}

/* hint == 0 : Use a default strategy
 * hint & 1  : Avoid mpqs(), use ellfacteur() after pollardbrent()
 * hint & 2  : Avoid first-stage ellfacteur() in favour of mpqs()
 * (may still fall back to ellfacteur() if mpqs() is not installed or gives up)
 * hint & 4  : Avoid even the pollardbrent() and squfof() stages. Put under
 *  the same governing  bit, for no good reason other than avoiding a
 *  proliferation of bits.
 * hint & 8  : Avoid final ellfacteur(); this may declare a composite to be
 *  prime.  */

/* stack housekeeping:  this routine may create one or more objects  (a new
 * factor, or possibly several, and perhaps one or more new exponents > 2) */
static long
ifac_crack(GEN *partial, GEN *where)
{
  long hint, cmp_res, exp1 = 1, exp2 = 1;
  pari_sp av;
  GEN factor = NULL, exponent;

#ifdef IFAC_DEBUG
  long lgp;
  if (!*partial || typ(*partial) != t_VEC)
    pari_err(typeer, "ifac_crack");
  if ((lgp = lg(*partial)) < ifac_initial_length)
    pari_err(talker, "partial impossibly short in ifac_crack");
  if (!(*where) || *where < *partial + 6 || *where > *partial + lgp - 3)
    pari_err(talker, "`*where\' out of bounds in ifac_crack");
  if (!(**where) || typ(**where) != t_INT)
    pari_err(typeer, "ifac_crack");
  if (gel(*where,2) != gen_0)
    pari_err(talker, "operand not known composite in ifac_crack");
#endif
  hint = itos(gel(*partial,2)) & 15;
  exponent = gel(*where,1);

  if (DEBUGLEVEL >= 3) fprintferr("IFAC: cracking composite\n\t%Z\n", **where);

  /* crack squares.  Quite fast due to the initial square residue test */
  if (DEBUGLEVEL >= 4) fprintferr("IFAC: checking for pure square\n");
  av = avma;
  while (Z_issquarerem(gel(*where,0), &factor))
  {
    if (DEBUGLEVEL >= 4)
      fprintferr("IFAC: found %Z =\n\t%Z ^2\n", **where, factor);
    affii(factor, gel(*where,0)); avma = av; factor = NULL;
    if (exponent == gen_1)
      gel(*where,1) = gen_2;
    else if (exponent == gen_2)
    { gel(*where,1) = utoipos(4); av = avma; }
    else
      affsi(itos(exponent) << 1, gel(*where,1));
    exponent = gel(*where,1);
    if (moebius_mode) return 0;	/* no need to carry on */
    exp1 = 2;
  } /* while Z_issquarerem */

  /* check whether our composite hasn't become prime */
  if (exp1 > 1 && hint != 15 && BSW_psp(gel(*where,0)))
  {
    gel(*where,2) = gen_1;
    if (DEBUGLEVEL >= 4) fprintferr("IFAC: factor %Z\n\tis prime\n",**where);
    return 0; /* bypass subsequent ifac_whoiswho() call */
  }
  /* still composite -- carry on */

  /* MPQS cannot factor prime powers. Do this even if MPQS is blocked by hint:
   * it is useful in bounded factorization */
  {
    ulong exp0 = 0, mask = 7;
    if (DEBUGLEVEL == 4) fprintferr("IFAC: checking for odd power\n");
    /* At debug levels > 4, is_357_power() prints something more informative */
    av = avma;
    while ( (exp1 = is_357_power(gel(*where,0), &factor, &mask)) )
    {
      if (exp2 == 1) exp2 = exp1; /* remember this after the loop */
      if (DEBUGLEVEL >= 4)
	fprintferr("IFAC: found %Z =\n\t%Z ^%ld\n", **where, factor, exp1);
      affii(factor, gel(*where,0)); avma = av; factor = NULL;
      if (exponent == gen_1)
      { gel(*where,1) = utoipos(exp1); av = avma; }
      else if (exponent == gen_2)
      { gel(*where,1) = utoipos(exp1<<1); av = avma; }
      else
        affsi(exp1 * itos(exponent), gel(*where,1));
      exponent = gel(*where,1);
      if (moebius_mode) return 0; /* no need to carry on */
    } /* while is_357_power */

    /* cutoff at 14 bits as trial division must have found everything below */
    while ( (exp1 = is_odd_power(gel(*where,0), &factor, &exp0, 15)) )
    {
      if (exp2 == 1) exp2 = exp1; /* remember this after the loop */
      if (DEBUGLEVEL >= 4)
	fprintferr("IFAC: found %Z =\n\t%Z ^%ld\n", **where, factor, exp1);
      affii(factor, gel(*where,0)); avma = av; factor = NULL;
      if (exponent == gen_1)
      { gel(*where,1) = utoipos(exp1); av = avma; }
      else if (exponent == gen_2)
      { gel(*where,1) = utoipos(exp1<<1); av = avma; }
      else
        affsi(exp1 * itos(exponent), gel(*where,1));
      exponent = gel(*where,1);
      if (moebius_mode) return 0; /* no need to carry on */
    } /* while is_odd_power */

    if (exp2 > 1 && hint != 15 && BSW_psp(gel(*where,0)))
    { /* Something nice has happened and our composite has become prime */
      gel(*where,2) = gen_1;
      if (DEBUGLEVEL >= 4)
        fprintferr("IFAC: factor %Z\n\tis prime\n", **where);
      return 0;	/* bypass subsequent ifac_whoiswho() call */
    }
  } /* odd power stage */

  if (!(hint & 4))
  { /* pollardbrent() Rho usually gets a first chance */
    if (DEBUGLEVEL >= 4) fprintferr("IFAC: trying Pollard-Brent rho method\n");
    factor = pollardbrent(gel(*where,0));
  }
  if (!factor && !(hint & 4))
  { /* Shanks' squfof() */
    if (DEBUGLEVEL >= 4)
      fprintferr("IFAC: trying Shanks' SQUFOF, will fail silently if input\n"
		 "      is too large for it.\n");
    factor = squfof(gel(*where,0));
  }
  if (!factor && !(hint & 2))
  { /* First ECM stage */
    if (DEBUGLEVEL >= 4) fprintferr("IFAC: trying Lenstra-Montgomery ECM\n");
    factor = ellfacteur(gel(*where,0), 0); /* do not insist */
  }
  if (!factor && !(hint & 1))
  { /* MPQS stage */
    if (DEBUGLEVEL >= 4) fprintferr("IFAC: trying MPQS\n");
    factor = mpqs(gel(*where,0));
  }
  if (!factor)
  {
    if (!(hint & 8))
    { /* still no luck? Final ECM stage, guaranteed to succeed */
      if (DEBUGLEVEL >= 4)
	fprintferr("IFAC: forcing ECM, may take some time\n");
      factor = ellfacteur(gel(*where,0), 1);
    }
    else
    { /* limited factorization */
      if (DEBUGLEVEL >= 2)
      {
        if (hint != 15)
          pari_warn(warner, "IFAC: unfactored composite declared prime");
        else
          pari_warn(warner, "IFAC: untested integer declared prime");

	/* don't print it out at level 3 or above, where it would appear
	 * several times before and after this message already */
	if (DEBUGLEVEL == 2) fprintferr("\t%Z\n",**where);
      }
      gel(*where,2) = gen_1;	/* might as well trial-divide by it... */
      return 1;
    }
  }
  if (typ(factor) == t_VEC) /* delegate this case */
    return ifac_insert_multiplet(partial, where, factor);
#ifdef IFAC_DEBUG
  else if (typ(factor) != t_INT)
  {
    fprintferr("IFAC: factorizer returned strange object to ifac_crack\n");
    outerr(factor);
    pari_err(bugparier, "factoring");
  }
#endif
  /* got single integer back:  work out the cofactor (in place) */
  if (!dvdiiz(gel(*where,0), factor, gel(*where,0)))
  {
    fprintferr("IFAC: factoring %Z\n", **where);
    fprintferr("\tyielded `factor\' %Z\n\twhich isn't!\n", factor);
    pari_err(bugparier, "factoring");
  }

  /* the factoring engines report the factor found when DEBUGLEVEL is
   * large enough;  let's tell about the cofactor */
  if (DEBUGLEVEL >= 4) fprintferr("IFAC: cofactor = %Z\n", **where);

  /* ok, now `factor' is one factor and **where is the other, find out which
   * is larger */
  cmp_res = cmpii(factor, gel(*where,0));
  /* mark factor /cofactor `unknown' */
  gel(*where,2) = NULL;
  gel(*where,-1) = NULL;
  gel(*where,-2) = isonstack(exponent) ? icopy(exponent) : exponent;
  *where -= 3;
  if (cmp_res < 0)
    gel(*where,0) = factor; /* common case */
  else if (cmp_res > 0)
  { /* factor > cofactor, rearrange */
    gel(*where,0) = gel(*where,3); /* move cofactor pointer to lowest slot */
    gel(*where,3) = factor;	/* save factor */
  }
  else pari_err(bugparier,"ifac_crack [Z_issquarerem miss]");
  return 2;
}

/* Don't collect garbage.  No diagnostics: the factoring engine should have
 * printed what it found. facvec contains slots of three components per factor;
 * repeated factors are allowed  (and their classes shouldn't contradict each
 * other whereas their exponents will be added up) */
static long
ifac_insert_multiplet(GEN *partial, GEN *where, GEN facvec)
{
  long j,k=1, lfv=lg(facvec)-1, nf=lfv/3, room=(long)(*where-*partial);
  /* one of the factors will go into the *where slot, so room is now 3 times
   * the number of slots we can use */
  long needroom = lfv - room;
  GEN sorted, auxvec = cgetg(nf+1, t_VEC), factor;
  long exponent = itos(gel(*where,1)); /* the old exponent */
  GEN newexp;

  if (DEBUGLEVEL >= 5) /* squfof may return a single squared factor as a set */
    fprintferr("IFAC: incorporating set of %ld factor(s)\n", nf);
  if (needroom > 0) /* one extra slot for paranoia, errm, future use */
    ifac_realloc(partial, where, lg(*partial) + needroom + 3);

  /* create sort permutation from the values of the factors */
  for (j=nf; j; j--) auxvec[j] = facvec[3*j-2]; /* just the pointers */
  sorted = sindexsort(auxvec);
  /* and readjust the result for the triple spacing */
  for (j=nf; j; j--) sorted[j] = 3*sorted[j]-2;

  /* store factors, beginning at *where, and catching any duplicates */
  **where = facvec[sorted[nf]];
  if ((newexp = gel(facvec,sorted[nf]+1)) != gen_1) /* new exponent > 1 */
  {
    if (exponent == 1)
      gel(*where,1) = isonstack(newexp) ? icopy(newexp) : newexp;
    else
      gel(*where,1) = mulsi(exponent, newexp);
  } /* if new exponent is 1, the old exponent already in place will do */
  (*where)[2] = facvec[sorted[nf]+2]; /* copy class */
  if (DEBUGLEVEL >= 6) fprintferr("\tstored (largest) factor no. %ld...\n", nf);

  for (j=nf-1; j; j--)
  {
    factor = gel(facvec,sorted[j]);
    if (equalii(factor, gel(*where,0)))
    {
      if (DEBUGLEVEL >= 6)
	fprintferr("\tfactor no. %ld is a duplicate%s\n", j, (j>1? "...": ""));
      /* update exponent, ignore class which would already have been set,
       * then forget current factor */
      if ((newexp = gel(facvec,sorted[j]+1)) != gen_1) /* new exp > 1 */
	gel(*where,1) = addis(gel(*where, 1), exponent * itos(newexp));
      else if (gel(*where,1) == gen_1 && exponent == 1)
        gel(*where,1) = gen_2;
      else
        gel(*where,1) = addis(gel(*where,1), exponent);
      if (moebius_mode) return 0; /* stop now, but with exponent updated */
      continue;
    }
    (*where)[-1] = facvec[sorted[j]+2];	/* class as given */
    if ((newexp = gel(facvec,sorted[j]+1)) != gen_1) /* new exp > 1 */
    {
      if (exponent == 1 && newexp == gen_2)
	gel(*where,-2) = gen_2;
      else /* exponent*newexp > 2 */
	gel(*where,-2) = mulsi(exponent, newexp);
    }
    else
    {
      gel(*where,-2) = (exponent == 1 ? gen_1 :
		      (exponent == 2 ? gen_2 :
		       utoipos(exponent)));/* inherit parent's exponent */
    }
    gel(*where,-3) = isonstack(factor) ? icopy(factor) : factor;
				/* keep components younger than *partial */
    *where -= 3;
    k++;
    if (DEBUGLEVEL >= 6)
      fprintferr("\tfactor no. %ld was unique%s\n",
		 j, (j>1 ? " (so far)..." : ""));
  }
  /* make the `sorted' object safe for garbage collection (it should be in the
   * garbage zone from everybody's perspective, but it's easy to do it) */
  *sorted = evaltyp(t_INT) | evallg(nf+1);
  return k;
}

static GEN
ifac_main(GEN *partial)
{
  GEN here = ifac_find(partial, partial);
  long nf, hint = itos(gel(*partial,2)) & 15;

  /* if nothing left, return gen_1 */
  if (!here) return gen_1;

  /* if we are in Moebius mode and have already detected a repeated factor,
   * stop right here.  Shouldn't happen */
  if (moebius_mode && gel(here,1) != gen_1)
  {
    if (DEBUGLEVEL >= 3)
      fprintferr("IFAC: main loop: repeated old factor\n\t%Z\n", *here);
    return gen_0;
  }

  /* loop until first entry is a finished prime.  May involve reallocations,
   * thus updates of *partial */
  while (gel(here,2) != gen_2)
  {
#if IFAC_DEBUG
    if (gel(here,2) == NULL)
    { /* unknown: something is wrong. Try to recover */
      pari_warn(warner, "IFAC: unknown factor seen in main loop");
      /* can only happen in Moebius mode */
      if (ifac_resort(partial, &here)) return gen_0;

      ifac_whoiswho(partial, &here, -1);
      ifac_defrag(partial, &here);
      continue;
    }
#endif
    if (gel(here,2) == gen_0)
    { /* composite: crack it */
      /* make sure there's room for another factor */
      if (here < *partial + 6)
      {	/* try defrag first */
	ifac_defrag(partial, &here);
	if (here < *partial + 6) /* no luck */
	  ifac_realloc(partial, &here, 1); /* guaranteed to work */
	  /* can't do a garbage collection here: we know too little about where
           * in the stack the old components were.*/
      }
      nf = ifac_crack(partial, &here);
      if (moebius_mode && gel(here,1) != gen_1) /* that was a power */
      {
	if (DEBUGLEVEL >= 3)
	  fprintferr("IFAC: main loop: repeated new factor\n\t%Z\n", *here);
	return gen_0;
      }
      /* deal with the new unknowns.  No sort: ifac_crack did it */
      ifac_whoiswho(partial, &here, nf);
      continue;
    }
    if (gel(here,2) == gen_1)
    { /* prime but not yet finished: finish it */
      if (ifac_divide(partial, &here))
      {
	if (moebius_mode)
	{
	  if (DEBUGLEVEL >= 3)
	    fprintferr("IFAC: main loop: another factor was divisible by\n"
	               "\t%Z\n", *here);
	  return gen_0;
	}
	ifac_defrag(partial, &here);
	(void)ifac_resort(partial, &here); /* sort new cofactors down */
        /* it doesn't matter whether this finds a repeated factor: we never
         * get to this point in Moebius mode */
	ifac_defrag(partial, &here); /* resort may have created new gaps */
	ifac_whoiswho(partial, &here, -1);
      }
      continue;
    }
    pari_err(talker, "non-existent factor class in ifac_main");
  } /* while */
  if (moebius_mode && gel(here,1) != gen_1)
  {
    if (DEBUGLEVEL >= 3)
      fprintferr("IFAC: after main loop: repeated old factor\n\t%Z\n", *here);
    return gen_0;
  }
  if (DEBUGLEVEL >= 4)
  {
    nf = (*partial + lg(*partial) - here - 3)/3;
    if (nf)
      fprintferr("IFAC: main loop: %ld factor%s left\n",
		 nf, (nf>1 ? "s" : ""));
    else
      fprintferr("IFAC: main loop: this was the last factor\n");
  }
  if (factor_add_primes && !(hint & 8))
  {
    GEN p = gel(here,0);
    if (lgefint(p)>3 || (ulong)p[2] > 0x1000000UL) (void)addprimes(p);
  }
  return here;
}

#define INIT_HERE(here) gel(here,2) = gel(here,1) = gel(here,0) = NULL;

/* Caller of the following should worry about stack management */
GEN
ifac_primary_factor(GEN *partial, long *exponent)
{
  GEN res, here = ifac_main(partial);

  if (here == gen_1) { *exponent = 0; return gen_1; }
  else if (here == gen_0) { *exponent = 0; return gen_0; }

  res = icopy(gel(here,0));
  *exponent = itos(gel(here,1));
  INIT_HERE(here); return res;
}

/* encapsulated routines */

/* prime/exponent pairs need to appear contiguously on the stack, but we also
 * need our data structure somewhere, and we don't know in advance how many
 * primes will turn up.  The following discipline achieves this:  When
 * ifac_decomp() is called, n should point at an object older than the oldest
 * small prime/exponent pair  (auxdecomp1() guarantees this).
 * We allocate sufficient space to accommodate several pairs -- eleven pairs
 * ought to fit in a space not much larger than n itself -- before calling
 * ifac_start().  If we manage to complete the factorization before we run out
 * of space, we free the data structure and cull the excess reserved space
 * before returning.  When we do run out, we have to leapfrog to generate more
 * (guesstimating the requirements from what is left in the partial
 * factorization structure);  room for fresh pairs is allocated at the head of
 * the stack, followed by an ifac_realloc() to reconnect the data structure and
 * move it out of the way, followed by a few pointer tweaks to connect the new
 * pairs space to the old one. This whole affair translates into a surprisingly
 * compact routine. */

/* ifac_decomp_break:
 *
 * Find primary factors of n until ifac_break return true, or n is factored if
 * ifac_break is NULL.
 */
/* ifac_break: return 1: stop factoring, 0 continue.
 *
 * state is private data for ifac_breack(), which must not leave anything on
 * the stack (except in state).
 * ifac_break is called in auxdecomp1 with here = NULL to register n, then
 * whenever a new factor is found. */

long
ifac_decomp_break(GEN n, long (*ifac_break)(GEN n,GEN pairs,GEN here,GEN state),
		  GEN state, long hint)
{
  pari_sp av = avma, lim = stack_lim(av, 1);
  long nb = 0;
  GEN part, here, workspc, pairs = (GEN)av;

  /* workspc will be doled out in pairs of smaller t_INTs. For n = prod p^{e_p}
   * (p not necessarily prime), need room to store all p and e_p [ cgeti(3) ],
   * bounded by
   *    sum_{p | n} ( log_{2^BIL} (p) + 6 ) <= log_{2^BIL} n + 6 log_2 n */
  workspc = new_chunk((expi(n) + 1) * 7);

  if (!n || typ(n) != t_INT) pari_err(typeer, "ifac_decomp");
  if (!signe(n)) pari_err(talker, "factoring 0 in ifac_decomp");

  part = ifac_start(n, 0, hint);
  here = ifac_main(&part);

  while (here != gen_1)
  {
    long lf = lgefint(*here);
    nb++;
    pairs -= lf;
    *pairs = evaltyp(t_INT) | evallg(lf);
    affii(gel(here,0), pairs);
    pairs -= 3;
    *pairs = evaltyp(t_INT) | evallg(3);
    affii(gel(here,1), pairs);
    if (ifac_break && (*ifac_break)(n,pairs,here,state))
    {
      if (DEBUGLEVEL >= 3) fprintferr("IFAC: (Partial fact.)Stop requested.\n");
      break;
    }
    INIT_HERE(here);
    here = ifac_main(&part);
    if (low_stack(lim, stack_lim(av,1)))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"[2] ifac_decomp");
      ifac_realloc(&part, &here, 0);
      part = gerepileupto((pari_sp)workspc, part);
    }
  }
  avma = (pari_sp)pairs;
  if (DEBUGLEVEL >= 3)
    fprintferr("IFAC: found %ld large prime (power) factor%s.\n",
	       nb, (nb>1? "s": ""));
  return nb;
}

long
ifac_decomp(GEN n, long hint)
{
  return ifac_decomp_break(n, NULL, gen_0, hint);
}

long
ifac_moebius(GEN n, long hint)
{
  long mu = 1;
  pari_sp av = avma, lim = stack_lim(av, 1);
  GEN part = ifac_start(n, 1, hint);
  GEN here = ifac_main(&part);

  while (here != gen_1 && here != gen_0)
  {
    if (itos(gel(here,1)) > 1) { here = gen_0; break; } /* won't happen */
    mu = -mu;
    INIT_HERE(here);
    here = ifac_main(&part);
    if (low_stack(lim, stack_lim(av,1)))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"ifac_moebius");
      ifac_realloc(&part, &here, 0);
      part = gerepileupto(av, part);
    }
  }
  avma = av; return here == gen_1? mu: 0;
}

long
ifac_issquarefree(GEN n, long hint)
{
  pari_sp av=avma, lim=stack_lim(av, 1);
  GEN part = ifac_start(n, 1, hint);
  GEN here = ifac_main(&part);

  while (here != gen_1 && here != gen_0)
  {
    if (itos(gel(here,1)) > 1) { here = gen_0; break; } /* won't happen */
    INIT_HERE(here);
    here = ifac_main(&part);
    if (low_stack(lim, stack_lim(av,1)))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"ifac_issquarefree");
      ifac_realloc(&part, &here, 0);
      part = gerepileupto(av, part);
    }
  }
  avma = av; return here == gen_1? 1: 0;
}

long
ifac_omega(GEN n, long hint)
{
  long omega=0;
  pari_sp av=avma, lim=stack_lim(av, 1);
  GEN part = ifac_start(n, 0, hint);
  GEN here = ifac_main(&part);

  while (here != gen_1)
  {
    omega++;
    INIT_HERE(here);
    here = ifac_main(&part);
    if (low_stack(lim, stack_lim(av,1)))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"ifac_omega");
      ifac_realloc(&part, &here, 0);
      part = gerepileupto(av, part);
    }
  }
  avma = av; return omega;
}

long
ifac_bigomega(GEN n, long hint)
{
  long Omega=0;
  pari_sp av=avma, lim=stack_lim(av, 1);
  GEN part = ifac_start(n, 0, hint);
  GEN here = ifac_main(&part);

  while (here != gen_1)
  {
    Omega += itos(gel(here,1));
    INIT_HERE(here);
    here = ifac_main(&part);
    if (low_stack(lim, stack_lim(av,1)))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"ifac_bigomega");
      ifac_realloc(&part, &here, 0);
      part = gerepileupto(av, part);
    }
  }
  avma = av; return Omega;
}

GEN
ifac_totient(GEN n, long hint)
{
  GEN res = cgeti(lgefint(n));
  pari_sp av=avma, tetpil, lim=stack_lim(av, 1);
  GEN phi = gen_1;
  GEN part = ifac_start(n, 0, hint);
  GEN here = ifac_main(&part);

  while (here != gen_1)
  {
    phi = mulii(phi, addsi(-1, gel(here,0)));
    if (gel(here,1) != gen_1)
    {
      if (gel(here,1) == gen_2)
	phi = mulii(phi, gel(here,0));
      else
	phi = mulii(phi, powiu(gel(here,0), itou(gel(here,1)) - 1));
    }
    INIT_HERE(here);
    here = ifac_main(&part);
    if (low_stack(lim, stack_lim(av,1)))
    {
      GEN *gsav[2];
      if(DEBUGMEM>1) pari_warn(warnmem,"ifac_totient");
      tetpil = avma;
      ifac_realloc(&part, &here, 0);
      phi = icopy(phi);
      gsav[0] = &phi; gsav[1] = &part;
      gerepilemanysp(av, tetpil, gsav, 2);
      /* don't preserve 'here', safer to pick it up again */
      here = ifac_find(&part, &part);
    }
  }
  affii(phi, res); avma = av; return res;
}

GEN
ifac_numdiv(GEN n, long hint)
{
  pari_sp av=avma, lim=stack_lim(av, 1);
  GEN tau = gen_1;
  GEN part = ifac_start(n, 0, hint);
  GEN here = ifac_main(&part);

  while (here != gen_1)
  {
    tau = mulis(tau, 1 + itos(gel(here,1)));
    INIT_HERE(here);
    here = ifac_main(&part);
    if (low_stack(lim, stack_lim(av,1)))
    {
      pari_sp tetpil = avma;
      GEN *gsav[2];
      if(DEBUGMEM>1) pari_warn(warnmem,"ifac_numdiv");
      tetpil = avma;
      ifac_realloc(&part, &here, 0);
      tau = icopy(tau);
      gsav[0] = &tau; gsav[1] = &part;
      gerepilemanysp(av, tetpil, gsav, 2);
      /* (see ifac_totient()) */
      here = ifac_find(&part, &part);
    }
  }
  return gerepileuptoint(av, tau);
}

GEN
ifac_sumdiv(GEN n, long hint)
{
  long exponent;
  pari_sp av=avma, lim=stack_lim(av, 1);
  GEN contrib, sigma = gen_1;
  GEN part = ifac_start(n, 0, hint);
  GEN here = ifac_main(&part);

  while (here != gen_1)
  {
    exponent = itos(gel(here,1));
    contrib = addsi(1, gel(here,0));
    for (; exponent > 1; exponent--)
      contrib = addsi(1, mulii(gel(here,0), contrib));
    sigma = mulii(sigma, contrib);
    INIT_HERE(here);
    here = ifac_main(&part);
    if (low_stack(lim, stack_lim(av,1)))
    {
      pari_sp tetpil = avma;
      GEN *gsav[2];
      if(DEBUGMEM>1) pari_warn(warnmem,"ifac_sumdiv");
      ifac_realloc(&part, &here, 0);
      sigma = icopy(sigma);
      gsav[0] = &sigma; gsav[1] = &part;
      gerepilemanysp(av, tetpil, gsav, 2);
      /* see ifac_totient() */
      here = ifac_find(&part, &part);
    }
  }
  return gerepileuptoint(av, sigma);
}

/* Assume k > 1. The calling function knows what to do with the other cases. */
GEN
ifac_sumdivk(GEN n, long k, long hint)
{
  long exponent;
  pari_sp av=avma, lim=stack_lim(av, 1);
  GEN contrib, q, sigma = gen_1;
  GEN part = ifac_start(n, 0, hint);
  GEN here = ifac_main(&part);

  while (here != gen_1)
  {
    exponent = itos(gel(here,1));
    q = powiu(gel(here,0), k);
    contrib = addsi(1, q);
    for (; exponent > 1; exponent--)
      contrib = addsi(1, mulii(q, contrib));
    sigma = mulii(sigma, contrib);
    INIT_HERE(here);
    here = ifac_main(&part);
    if (low_stack(lim, stack_lim(av,1)))
    {
      pari_sp tetpil = avma;
      GEN *gsav[2];
      if(DEBUGMEM>1) pari_warn(warnmem,"ifac_sumdivk");
      ifac_realloc(&part, &here, 0);
      sigma = icopy(sigma);
      gsav[0] = &sigma; gsav[1] = &part;
      gerepilemanysp(av, tetpil, gsav, 2);
      /* see ifac_totient() */
      here = ifac_find(&part, &part);
    }
  }
  return gerepileuptoint(av, sigma);
}
