/* $Id: mpqs.c 12099 2010-01-29 15:38:23Z bill $

Copyright (C) 2000  The PARI group.

This file is part of the PARI/GP package.

PARI/GP is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation. It is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY WHATSOEVER.

Check the License for details. You should have received a copy of it, along
with the package; see the file 'COPYING'. If not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA. */

/* Written by Thomas Papanikolaou and Xavier Roblot
 *
 * Implementation of the Self-Initializing Multi-Polynomial Quadratic Sieve
 * based on code developed as part of the LiDIA project
 * (http://www.informatik.tu-darmstadt.de/TI/LiDIA/)
 *
 * Extensively modified by The PARI group.
 */
/* Notation commonly used in this file, and sketch of algorithm:
 *
 * Given an odd integer N > 1 to be factored, we throw in a small odd and
 * squarefree multiplier k so as to make kN congruent 1 mod 4 and to have
 * many small primes over which X^2 - kN splits.  We compute a factor base
 * FB of such primes, and then essentially look for values x0 such that
 * Q0(x0) = x0^2 - kN can be decomposed over this factor base, up to a
 * possible factor dividing k and a possible "large prime".  Relations
 * involving the latter can be combined into full relations (working mod N
 * now) which don't, and full relations, by Gaussian elimination over the
 * 2-element field for the exponent vectors, will finally lead us to an
 * expression X^2 - Y^2 divisible by N and hopefully to a nontrivial
 * splitting when we compute gcd(X +- Y, N).  Note that this can never
 * split prime powers.  (Any odd prime dividing the X - Y factor, say, will
 * divide it to the same power as it divides N.)
 *
 * Candidates x0 are found by sieving along arithmetic progressions modulo
 * the small primes in the factor base, and evaluation of candidates picks
 * out those x0 where many of these APs happen to coincide, resulting in a
 * highly divisible Q0(x0).
 *
 * The Multi-Polynomial version improves this by choosing a modest subset of
 * factor base primes and forcing these to divide Q0(x).  Let A be the product
 * of the chosen primes, and write Q(x) = Q0(2Ax + B) = (2Ax + B)^2 - kN =
 * 4A(Ax^2 + Bx + C), where B has been suitably chosen.  For each A, there
 * are 2^omega_A possible values for B when A is the product of omega_A
 * distinct primes, but we'll use only half of these, since the other half
 * is more easily covered by exploiting the symmetry x <-> -x of the original
 * quadratic.  The "Self-Initializating" bit refers to the fact that switching
 * from one B to the next can be done very fast, whereas switching to the
 * next A involves some recomputation from scratch.  (C is never needed
 * explicitly except for debug diagnostics.)  Thus we can very quickly run
 * through an entire "cohort" of polynomials sharing the same A.
 *
 * The sieve now ranges over values x0 such that |x0| < M  (we use x = x0 + M
 * as the non-negative array subscript).  The coefficients A are chosen so
 * that A*M is approximately sqrt(kN).  Then |B| will be bounded by about
 * (j+4)*A, and |C| = -C will be about (M/4)*sqrt(kN), so Q(x0)/(4A) will
 * take values roughly between -|C| and 3|C|.
 *
 * There are numerous refinements to this basic outline  (e.g. it is more
 * efficient to _not_ use the very smallest primes in the FB for sieving,
 * incorporating them only after selecting candidates).  The substition of
 * 2Ax+B into X^2 - kN, with odd B, forces 2 to occur;  when kN is 1 mod 8,
 * it will always occur at least to the 3rd power;  when kN = 5 mod 8, it
 * will always occur exactly to the 2nd power.  We never sieve on 2 and we
 * always pull out the power of 2 directly  (which is easy, of course).
 * The prime factor(s) of k will show up whenever 2Ax + B has a factor in
 * common with k;  we don't sieve on these either but can easily recognize
 * them in a candidate.
 */

#include "pari.h"
#include "paripriv.h"

/* #define MPQS_DEBUG_VERBOSE 1 */
/* histograms are pretty, but don't help performance after all (see below) */
/* #define MPQS_USE_HISTOGRAMS */

#include "mpqs.h"

/*********************************************************************/
/**                                                                 **/
/**                         INITIAL SIZING                          **/
/**                                                                 **/
/*********************************************************************/

static char *i2str(GEN x) { return itostr(x, signe(x) < 0); }

/* number of decimal digits of argument - for parameter choosing and for
 * diagnostics */
static long
decimal_len(GEN N)
{
  pari_sp av = avma;
  long d = strlen(i2str(N));
  avma = av; return d;
}

/* To be called after choosing k and putting kN into the handle:
 * Pick up the requested parameter set for the given size of kN in decimal
 * digits and fill in various fields in the handle.  Return 0 when kN is
 * too large, 1 when we're ok. */
static int
mpqs_set_parameters(mpqs_handle_t *h)
{
  long i;
  mpqs_parameterset_t *P;
  double mb;

  h->digit_size_kN = decimal_len(h->kN);
  if (h->digit_size_kN <= 9)
    i = 0;
  else if (h->digit_size_kN > MPQS_MAX_DIGIT_SIZE_KN)
    return 0;
  else
    i = h->digit_size_kN - 9;

  /* (cf PARI bug#235) the following has always been, and will remain,
   * a moving target... increased thresholds from 64, 80 to 79, 86
   * respectively --GN20050601.  Note that the new values correspond to
   * kN having >= 86 or >= 95 decimal digits, respectively.  Note also
   * that the current sizing parameters for 90 or more digits are based
   * on 100% theory and 0% practice. */
  if (i >= 79)
    pari_warn(warner, "MPQS: factoring this number will take %s hours:\nN = %Z",
        i >= 86 ? "many": "several", h->N);

  if (DEBUGLEVEL >= 5)
  {
    fprintferr("MPQS: kN = %Z\n", h->kN);
    fprintferr("MPQS: kN has %ld decimal digits\n", h->digit_size_kN);
  }

  P = &(mpqs_parameters[i]);
  h->tolerance        = P->tolerance;
  h->lp_scale         = P->lp_scale;
  /* make room for prime factors of k if any: */
  h->size_of_FB       = P->size_of_FB + h->_k.omega_k;
  /* for the purpose of Gauss elimination etc., prime factors of k behave
   * like real FB primes, so take them into account when setting the goal: */
  h->target_no_rels   = (h->size_of_FB >= 200 ?
                         h->size_of_FB + 70 :
                         (mpqs_int32_t)(h->size_of_FB * 1.35));
  h->M                = P->M;
  h->omega_A          = P->omega_A;
  h->no_B             = 1UL << (P->omega_A - 1);
  h->pmin_index1      = P->pmin_index1;
  /* certain subscripts into h->FB should also be offset by omega_k: */
  h->index0_FB        = 3 + h->_k.omega_k;
  /* following are converted from % to parts per thousand: */
  h->first_sort_point = 10 * P->first_sort_point;
  h->sort_pt_interval = 10 * P->sort_pt_interval;

  mb = (h->size_of_FB + 1)/(8.*1048576.) * h->target_no_rels;
  if (mb > 128.)
  {
    pari_warn(warner,
        "MPQS: Gauss elimination will require more than\n\t128MBy of memory");
    if (DEBUGLEVEL >= 1)
      fprintferr("\t(estimated memory needed: %4.1fMBy)\n", mb);
  }

  return 1;
}

/*********************************************************************/
/**                                                                 **/
/**                       OBJECT HOUSEKEEPING                       **/
/**                                                                 **/
/*********************************************************************/

/* The sub-constructors for the pieces of the handle will be called in the
 * same order as their appearance here, and the later ones in part rely on
 * the earlier ones having filled in some fields.
 * There's a single destructor to handle all cleanup at the end  (except
 * for mpqs() itself resetting avma). */

/* main handle constructor */
static mpqs_handle_t *
mpqs_handle_ctor(GEN N)
{
  mpqs_handle_t *h = (mpqs_handle_t *) gpmalloc(sizeof(mpqs_handle_t));
  memset((void *)h, 0, sizeof(mpqs_handle_t));
  h->N = N;
#ifdef MPQS_DEBUG_VERBOSE
  fprintferr("MPQS DEBUG: created handle @0x%p\n", (void *)h);
#endif
  return h;
}

/* factor base constructor. Really a home-grown memalign(3c) underneath.
 * We don't want FB entries to straddle L1 cache line boundaries, and
 * malloc(3c) only guarantees alignment adequate for all primitive data
 * types of the platform ABI - typically to 8 or 16 byte boundaries.
 * Also allocate the inv_A_H array.
 * The FB array pointer is returned for convenience */
static mpqs_FB_entry_t *
mpqs_FB_ctor(mpqs_handle_t *h)
{
  /* leave room for slots 0, 1, and sentinel slot at the end of the array */
  long size_FB_chunk = (h->size_of_FB + 3) * sizeof(mpqs_FB_entry_t);
  /* like FB, except this one does not have a sentinel slot at the end */
  long size_IAH_chunk = (h->size_of_FB + 2) * sizeof(mpqs_inv_A_H_t);
  char *fbp = gpmalloc(size_FB_chunk + 64);
  char *iahp = gpmalloc(size_IAH_chunk + 64);
  long fbl, iahl;

  h->FB_chunk = (void *)fbp;
  h->invAH_chunk = (void *)iahp;
  /* round up to next higher 64-bytes-aligned address */
  fbl = (((long)fbp) + 64) & ~0x3FL;
  /* and put the actual array there */
  h->FB = (mpqs_FB_entry_t *)fbl;
#ifdef MPQS_DEBUG_VERBOSE
  fprintferr("MPQS DEBUG: FB chunk allocated @0x%p\n", (void *)fbp);
  fprintferr("MPQS DEBUG: FB aligned to 0x%p\n", (void *)fbl);
#endif

  iahl = (((long)iahp) + 64) & ~0x3FL;
  h->inv_A_H = (mpqs_inv_A_H_t *)iahl;
#ifdef MPQS_DEBUG_VERBOSE
  fprintferr("MPQS DEBUG: inv_A_H chunk allocated @0x%p\n", (void *)iahp);
  fprintferr("MPQS DEBUG: inv_A_H aligned to 0x%p\n", (void *)iahl);
#endif

  return (mpqs_FB_entry_t *)fbl;
}

/* sieve array constructor;  also allocates the candidates array, the
 * histograms, and temporary storage for relations under construction */
static void
mpqs_sieve_array_ctor(mpqs_handle_t *h)
{
  long size = (h->M << 1) + 1;
  mpqs_int32_t size_of_FB = h->size_of_FB;

  h->sieve_array =
    (unsigned char *) gpmalloc(size * sizeof(unsigned char));
#ifdef MPQS_DEBUG_VERBOSE
  fprintferr("MPQS DEBUG: sieve_array allocated @0x%p\n",
             (void *)h->sieve_array);
#endif
  h->sieve_array_end =
    h->sieve_array + size - 2;
  h->sieve_array_end[1] = 255; /* sentinel */

  h->candidates =
    (long *) gpmalloc(MPQS_CANDIDATE_ARRAY_SIZE * sizeof(long));
#ifdef MPQS_DEBUG_VERBOSE
  fprintferr("MPQS DEBUG: candidates table allocated @0x%p\n",
             (void *)h->candidates);
#endif

  /* Room needed for string representation of a relation - worst case:
   * + leading " 1 1"
   * + trailing " 0\n" with final NUL character
   * + in between up to size_of_FB pairs each consisting of an exponent, a
   *   subscript into FB, and two spaces.
   * Subscripts into FB fit into 5 digits, and exponents fit into 3 digits
   * with room to spare -- anything needing 3 or more digits for the
   * subscript must come with an exponent of at most 2 digits. Moreover the
   * product of the first 58 primes is larger than 10^110  (and the righthand
   * sides of proto-relations are much smaller than kN: on the order of
   * M*sqrt(kN)),  so there cannot be more than 60 pairs in all, even if
   * size_of_FB > 10^5. --GN */

  /* whereas mpqs_self_init() uses size_of_FB+1, we just use the size as
   * it is, not counting FB[1], to start off the following estimate */
  if (size_of_FB > 60) size_of_FB = 60;
  h->relations = (char *) gpmalloc((8 + size_of_FB * 9) * sizeof(char));
  /* and for tracking which primes occur in the current relation: */
  h->relaprimes = (long *) gpmalloc((size_of_FB << 1) * sizeof(long));

#ifdef MPQS_USE_HISTOGRAMS
  /* histograms to be used only when kN isn't very small */
  if (h->size_of_FB > MPQS_MIN_SIZE_FB_FOR_HISTO) {
    h->do_histograms = 1;
    h->histo_full = (long *) gpmalloc(128 * sizeof(long));
    h->histo_lprl = (long *) gpmalloc(128 * sizeof(long));
    h->histo_drop = (long *) gpmalloc(128 * sizeof(long));
    memset((void *)(h->histo_full), 0, 128 * sizeof(long));
    memset((void *)(h->histo_lprl), 0, 128 * sizeof(long));
    memset((void *)(h->histo_drop), 0, 128 * sizeof(long));
  }
#endif
}

/* mpqs() calls the following (after recording avma) to allocate GENs for
 * the current polynomial and self-initialization scratchpad data on the
 * PARI stack.  This space is released by mpqs() itself at the end. */
static void
mpqs_poly_ctor(mpqs_handle_t *h)
{
  mpqs_int32_t i;
  long size_per = h->omega_A * sizeof(mpqs_per_A_prime_t);

  h->per_A_pr = (mpqs_per_A_prime_t *) gpmalloc(size_per);
  memset((void *)(h->per_A_pr), 0, size_per);
  /* Sizing:  A is the product of omega_A primes, each well below word
   * size.
   * |B| is bounded by (omega_A + 4) * A, so can have at most one word
   * more, and that's generous.
   * |C| is less than A*M^2, so can take at most two words more than A.
   * The array H holds residues modulo A, so the same size as used for A
   * is sufficient. */
  h->A = cgeti(h->omega_A + 2);
  h->B = cgeti(h->omega_A + 3);
#ifdef MPQS_DEBUG
  h->C = cgeti(h->omega_A + 4);
#endif
  for (i = 0; i < h->omega_A; i++)
    h->per_A_pr[i]._H = cgeti(h->omega_A + 2);
  /* the handle starts out all zero, so in particular bin_index and index_i
   * are initially 0.
   * XX index_j currently initialized in mqps() but this is going to
   * change. */
}

/* main handle destructor, also cleans up all other allocated pieces
 * (except for stuff created on the PARI stack which the caller should
 * deal with by resetting avma) */
static void
mpqs_handle_dtor(mpqs_handle_t *h)
{
#define myfree(x) if(x) free((void*)x)
  myfree((h->per_A_pr));
  myfree((h->relaprimes));
  myfree(h->relations);

#ifdef MPQS_USE_HISTOGRAMS
  myfree((h->histo_drop));
  myfree((h->histo_lprl));
  myfree((h->histo_full));
#endif

  myfree((h->candidates));
  myfree((h->sieve_array));
  myfree((h->invAH_chunk));
  myfree((h->FB_chunk));
  myfree(h);
}

/* XX todo - relationsdb handle */

/*********************************************************************/
/**                                                                 **/
/**                        FACTOR BASE SETUP                        **/
/**                                                                 **/
/*********************************************************************/

/* our own pointer to PARI's or to our own prime diffs table.
 * NB the latter is never freed, unless we need to replace it with
 * an even larger one. */
static byteptr mpqs_diffptr = NULL;
static long mpqs_prime_count = 0;
static int mpqs_use_our_diffptr = 0;

/* return next prime larger than p, using *primes_ptr on the diffptr table
 * first and pari's other wits after that */
static byteptr
mpqs_iterate_primes(ulong *p, byteptr primes_ptr)
{
  ulong prime = *p;
  if (*primes_ptr)
    NEXT_PRIME_VIADIFF(prime,primes_ptr);
  else
  {
    pari_sp av = avma;
    prime = itou(nextprime(utoipos(prime + 1)));
    avma = av;
  }
  *p = prime; return primes_ptr;
}

/* fill in the best-guess multiplier k for N. We force kN = 1 mod 4.
 * Caller should proceed to fill in kN */
static void
mpqs_find_k(mpqs_handle_t *h)
{
  pari_sp av = avma;
  mpqs_multiplier_t *cand_k;
  long best_i = -1 /* never best */, k, N_mod_4 = mod4(h->N);
  ulong p;
  GEN kN;
  double best_value = -1000. /* essentially -infinity */, value, dp;
  long i, j;
  byteptr primes_ptr;

  for (i = 0; i < MPQS_POSSIBLE_MULTIPLIERS; i++)
  {
    cand_k = &cand_multipliers[i];
    k = cand_k->k;
    if ((k & 3) == N_mod_4) /* kN = 1 (mod 4) */
    {
      value = -0.7 * log2 ((double) k);
      kN = mulis(h->N, k);
      if (mod8(kN) == 1) value += 1.38629;

      j = 0; p = 0;
      primes_ptr = diffptr; /* that's PARI's, not our private one */
      while (j <= MPQS_MULTIPLIER_SEARCH_DEPTH)
      {
        primes_ptr = mpqs_iterate_primes(&p, primes_ptr);
        if (krouu(umodiu(kN, p), p) == 1)
        {
          j++;
          dp = log2((double) p) / p;
          value += (k % p == 0) ? dp : 2 * dp;
        }
      }
      if (value > best_value) { best_value = value; best_i = i; }
    }
  }
  avma = av; h->_k = cand_multipliers[best_i];
}

/******************************/

/* guesstimate up to what size we're going to need precomputed small primes */
static long
mpqs_find_maxprime(long size)
{
  double x;

  if (size < 16000) return 176000;
  x  = log((double)size);
  x += log(x) - 0.9427;
  return (long)(x * size);
}

/* return the number of primes in mpqs_diffptr */
static long
mpqs_count_primes(void)
{
  byteptr p = mpqs_diffptr;
  long gaps = 0;

  for ( ; *p; p++)
    if (*p == DIFFPTR_SKIP) gaps++;
  return (p - mpqs_diffptr - gaps);
}

/* Create a factor base of size primes p_i such that legendre(k*N, p_i) != -1
 * We could have shifted subscripts down from their historical arrangement,
 * but this seems too risky for the tiny potential gain in memory economy.
 * The real constraint is that the subscripts of anything which later shows
 * up at the Gauss stage must be nonnegative, because the exponent vectors
 * there use the same subscripts to refer to the same FB entries.  Thus in
 * particular, the entry representing -1 could be put into FB[0], but could
 * not be moved to FB[-1] (although mpqs_FB_ctor() could be easily adapted
 * to support negative subscripts).-- The historically grown layout is:
 * FB[0] is unused.
 * FB[1] is not explicitly used but stands for -1.
 * FB[2] contains 2 (always).
 * Before we are called, the size_of_FB field in the handle will already have
 * been adjusted by _k.omega_k, so there's room for the primes dividing k,
 * which when present will occupy FB[3] and following.
 * The "real" odd FB primes begin at FB[h->index0_FB].
 * FB[size_of_FB+1] is the last prime p_i.
 * FB[size_of_FB+2] is a sentinel to simplify some of our loops.
 * Thus we allocate size_of_FB+3 slots for FB.
 * If a prime factor of N is found during the construction, it is returned
 * in f, otherwise f = 0  (such a factor necessarily fits into a C long).
 * We use PARI's normal diffptr array if it's large enough, or otherwise
 * our own because we don't want to pull out PARI's array under our caller
 * if someone was in the middle of using it while calling us. (But note
 * that this makes ourselves non-reentrant.) */
/* XX convert (everything!) to new FB entry format */
/* XX also fill in a number of other FB entry fields */
/* returns the FB array pointer for convenience */
static mpqs_FB_entry_t *
mpqs_create_FB(mpqs_handle_t *h, ulong *f)
{
  ulong p = 0;
  mpqs_int32_t size = h->size_of_FB;
  long i, kr /* , *FB */;
  mpqs_uint32_t k = h->_k.k;
  mpqs_FB_entry_t *FB;          /* for ease of reference */
  byteptr primes_ptr;

  FB = mpqs_FB_ctor(h);

  /* tentatively pick up PARI's current differences-of-small-primes array
   * unless we already have our own */
  if (!mpqs_use_our_diffptr) mpqs_diffptr = diffptr;

  if ((mpqs_prime_count? mpqs_prime_count: mpqs_count_primes()) < 3 * size)
  {
    /* not large enough - must use our own then */
    long newsize = 3 * mpqs_find_maxprime(size);
    if (mpqs_use_our_diffptr) free((void *) mpqs_diffptr);
    if (DEBUGLEVEL >= 2)
      fprintferr("MPQS: precomputing auxiliary primes up to %ld\n", newsize);
    /* the following three assignments must happen in this order, to
     * safeguard against corruption when we are being interrupted at
     * the wrong moment: */
    mpqs_diffptr = initprimes(newsize);
    mpqs_use_our_diffptr = 1;   /* and will remain true forever */
    mpqs_prime_count = mpqs_count_primes(); /* count once and remember */
  }

  if (MPQS_DEBUGLEVEL >= 7) fprintferr("MPQS: FB [-1,2");
  FB[2].fbe_p = 2;
  /* the fbe_logval and the fbe_sqrt_kN for 2 are never used */
  FB[2].fbe_flags = MPQS_FBE_CLEAR;
  primes_ptr = mpqs_diffptr;
  primes_ptr = mpqs_iterate_primes(&p, primes_ptr); /* move past 2 */

  /* the first loop executes h->_k.omega_k = 0, 1, or 2 times */
  for (i = 3; i < h->index0_FB; i++)
  {
    mpqs_uint32_t kp = (ulong)h->_k.kp[i-3];
    if (MPQS_DEBUGLEVEL >= 7) fprintferr(",<%lu>", (ulong)kp);
    FB[i].fbe_p = kp;
    /* we *could* flag divisors of k here, but so far I see no need,
     * and no flags bit has been assigned for the purpose */
    FB[i].fbe_flags = MPQS_FBE_CLEAR;
    FB[i].fbe_flogp = (float) log2((double) kp);
    FB[i].fbe_sqrt_kN = 0;
  }
  /* now i == h->index0_FB */
  while (i < size + 2)
  {
    primes_ptr = mpqs_iterate_primes(&p, primes_ptr);

    if ((p > k) || (k%p != 0))
    {
      ulong kN_mod_p = umodiu(h->kN, p);
      kr = krouu(kN_mod_p, p);
      if (kr != -1)
      {
        if (kr == 0)
        {
          if (MPQS_DEBUGLEVEL >= 7)
            fprintferr(",%lu...] Wait a second --\n", p);
          *f = p;
          return FB;
        }
        if (MPQS_DEBUGLEVEL >= 7) fprintferr(",%lu", p);
        FB[i].fbe_p = (mpqs_uint32_t) p;
        FB[i].fbe_flags = MPQS_FBE_CLEAR;
        /* dyadic logarithm of p; single precision suffices */
        FB[i].fbe_flogp = (float) log2((double)p);
        /* cannot yet fill in fbe_logval because the scaling multiplier for
         * it depends on the largest prime in FB, which we're still in the
         * process of finding here! */

        /* x_i such that x_i^2 = (kN % p_i) mod p_i */
        FB[i++].fbe_sqrt_kN =
          (mpqs_uint32_t) Fl_sqrt(kN_mod_p, p);
      }
    }
  }

  if (MPQS_DEBUGLEVEL >= 7) fprintferr("]\n");

  FB[i].fbe_p = 0;              /* sentinel */
  h->largest_FB_p = FB[i-1].fbe_p; /* at subscript size_of_FB + 1 */

#ifdef MPQS_EXPERIMENTAL_COUNT_SIEVE_HITS
  /* temporary addition: ------------8<---------------------- */
  {
    double hitsum;
    long j;
    for (j = 300; j <= 1000; j += 100)
    {
      hitsum = 0.;
      for (i = j; i < size + 2; i++)
      {
        p = FB[i].fbe_p; if (p == 0) break;
        hitsum += ((4.0 * h->M) / p);
      }
      fprintferr("MPQS DEBUG: hits from FB[%d]=%ld up: %10.2f\n",
                 j, (long)(FB[j].fbe_p), hitsum);
    }
  }
#endif
  /* locate the smallest prime that will be used for sieving */
  for (i = h->index0_FB; FB[i].fbe_p != 0; i++)
    if (FB[i].fbe_p >= h->pmin_index1) break;
  h->index1_FB = i;
  /* assert: with our parameters this will never fall of the end of the FB */

  *f = 0;
  return FB;
}

/*********************************************************************/
/**                                                                 **/
/**                      MISC HELPER FUNCTIONS                      **/
/**                                                                 **/
/*********************************************************************/

/* Effect of the following:  multiplying the base-2 logarithm of some
 * quantity by log_multiplier will rescale something of size
 *    log2 ( sqrt(kN) * M / (largest_FB_prime)^tolerance )
 * to 232.  Note that sqrt(kN) * M is just A*M^2, the value our polynomials
 * take at the outer edges of the sieve interval.  The scale here leaves
 * a little wiggle room for accumulated rounding errors from the approximate
 * byte-sized scaled logarithms for the factor base primes which we add up
 * in the sieving phase.-- The threshold is then chosen so that a point in
 * the sieve has to reach a result which, under the same scaling, represents
 *    log2 ( sqrt(kN) * M / (largest_FB_prime)^tolerance )
 * in order to be accepted as a candidate. */
/* The old formula was...
 *   log_multiplier =
 *      127.0 / (0.5 * log2 (handle->dkN)
 *               + log2((double)M)
 *               - tolerance * log2((double)handle->largest_FB_p)
 *               );
 * and we used to use this with a constant threshold of 128. */

/* XX We also used to divide log_multiplier by an extra factor 2, and in
 * XX compensation we were multiplying by 2 when the fbe_logp fields were
 * XX being filled in, making all those bytes even.  Tradeoff: the extra
 * XX bit of precision is helpful, but interferes with a possible sieving
 * XX optimization  (artifically shift right the logp's of primes in A,
 * XX and just run over both arithmetical progressions  (which coincide in
 * XX this case)  instead of skipping the second one, to avoid the conditional
 * XX branch in the mpqs_sieve() loops).  We could still do this, but might
 * XX lose a little bit accuracy for those primes.  Probably no big deal... */
static void
mpqs_set_sieve_threshold(mpqs_handle_t *h)
{
  mpqs_FB_entry_t *FB = h->FB;
  long i;
  double log_maxval;
  double log_multiplier;

  h->l2sqrtkN = 0.5 * log2(h->dkN);
  h->l2M = log2((double)h->M);
  log_maxval = h->l2sqrtkN + h->l2M - MPQS_A_FUDGE;
  log_multiplier = 232.0 / log_maxval;
  h->sieve_threshold =
    (unsigned char) (log_multiplier *
                     (log_maxval
                      - h->tolerance * log2((double)h->largest_FB_p)
                      )
                     ) + 1;
  /* That "+ 1" really helps - we may want to tune towards somewhat smaller
   * tolerances  (or introduce self-tuning one day)... */

  /* If this turns out to be <128, scream loudly.
   * That means that the FB or the tolerance or both are way too
   * large for the size of kN.  (Normally, the threshold should end
   * up in the 150...170 range.) */
  if (h->sieve_threshold < 128) {
    h->sieve_threshold = 128;
    pari_warn(warner,
        "MPQS: sizing out of tune, FB size or tolerance\n\ttoo large");
  }

  /* Now fill in the byte-sized approximate scaled logarithms of p_i */
  if (DEBUGLEVEL >= 5)
  {
    fprintferr("MPQS: computing logarithm approximations for p_i in FB\n");
  }
  for (i = h->index0_FB; i < h->size_of_FB + 2; i++)
  {
    FB[i].fbe_logval =
      (unsigned char) (log_multiplier * FB[i].fbe_flogp);
  }
}

/* Given the partially populated handle, find the optimum place in the FB
 * to pick prime factors for A from.  The lowest admissible subscript is
 * index0_FB, but unless kN is very small, we'll stay away a bit from that.
 * The highest admissible, in a pinch, is size_of_FB + 1, where the largest
 * FB prime resides.  The ideal corner is about (sqrt(kN)/M) ^ (1/omega_A),
 * so that A will end up of size comparable to sqrt(kN)/M;  experimentally
 * it seems desirable to stay very slightly below this.  Moreover, the
 * selection of the individual primes happens to pari_err on the large side, for
 * which it is wise to compensate a bit.  This is what the (small positive)
 * quantity MPQS_A_FUDGE is for.
 * We rely on a few auxiliary fields in the handle to be already set by
 * mqps_set_sieve_threshold() before we are called.
 * This function may fail under highly unfortunate circumstances, so we
 * report back about our success to the caller (mpqs()), allowing it to
 * bail out after emitting a warning. */
static int
mpqs_locate_A_range(mpqs_handle_t *h)
{
  /* i will be counted up to the desirable index2_FB + 1, and omega_A is never
   * less than 3, and we want
   *   index2_FB - (omega_A - 1) + 1 >= index0_FB + omega_A - 3,
   * so: */
  long i = h->index0_FB + 2*(h->omega_A) - 4;
  double l2_target_pA;
  mpqs_FB_entry_t *FB = h->FB;

  h->l2_target_A = (h->l2sqrtkN - h->l2M - MPQS_A_FUDGE);
  l2_target_pA = h->l2_target_A / h->omega_A;

  /* find the sweet spot, normally shouldn't take long */
  while ((FB[i].fbe_p != 0) && (FB[i].fbe_flogp <= l2_target_pA)) i++;

#ifdef MPQS_DEBUG_LOCATE_A_RANGE
  fprintferr("MPQS DEBUG: omega_A=%ld, index0=%ld, i=%ld\n",
             (long) h->omega_A, (long) h->index0_FB, i);
#endif

  /* check whether this hasn't walked off the top end... */
  /* The following should actually NEVER happen. */
  if (i > h->size_of_FB - 3)
  { /* this isn't going to work at all. */
    pari_warn(warner,
        "MPQS: sizing out of tune, FB too small or\n\tway too few primes in A");
    return 0;
  }

  /* GN 20050723 - comparison against index1_FB removed. */
  h->index2_FB = i - 1;
#ifdef MPQS_DEBUG_LOCATE_A_RANGE
  fprintferr("MPQS DEBUG: index2_FB = %ld\n", i - 1);
#endif
  /* GN20050723
   * assert: index0_FB + (omega_A - 3) [the lowest FB subscript eligible to
   * be used in picking primes for A]  plus  (omega_A - 2)  does not exceed
   * index2_FB  [the subscript from which the choice of primes for A starts,
   * putting omega_A - 1 of them at or below index2_FB, and the last and
   * largest one above, cf. mpqs_si_choose_primes() below].
   * Moreover, index2_FB indicates the last prime below the ideal size, unless
   * (when kN is very small) the ideal size was too small to use. */

  return 1;
}


/*********************************************************************/
/**                                                                 **/
/**                HISTOGRAMS AND THRESHOLD FEEDBACK                **/
/**                                                                 **/
/*********************************************************************/

#ifdef MPQS_USE_HISTOGRAMS

/* The histogram-related code is left in this file, but all under the
 * above #ifdef, and disabled by default.  I'm finding that:
 * - merely keeping the numbers updated in mpqs_eval_cand() below  (and
 *   keeping the "negligible" 1.5 or 3KBys' worth of extra arrays in use)
 *   causes us to run quite noticeably slower: 8-10% for a 73-digit number,
 * - mpqs_eval_cand() has already become so much faster than it used to be
 *   that raising the threshold to get rid of many low-valued unpromising
 *   candidates does not save any significant time, and even losing a pretty
 *   small number of additional LP relations actually harms us by lowering
 *   the efficiency of LP relation combining.
 * (The first point might be due merely to code bloat and less effective
 * compiler optimizations - I'm not sure about that.)
 * Just for getting a visual impression of how the sieve is performing,
 * however, this is nice to have available.  Just turn on the #define at
 * the top of the file and recompile with it. --GN2005-02-03
 */

/* histogram evaluation happens very infrequently if at all.  So we'll do
 * all the adding and putting-into-relation-with here, while mpqs_eval_cand()
 * merely bumps one cell at a time and doesn't keep running totals. */

static void
mpqs_print_histo(mpqs_handle_t *h)
{
  long i, tot = 0;

  if (!h->do_histograms) return;

  fprintferr("\nMPQS: values from sieve vs. distribution of evaluated candidates:\n");
  fprintferr("   val  ___full __lprel ___none ___total\n");
  for (i = 127; i >= 0; i--)
  {
    long rowtot = h->histo_full[i] + h->histo_lprl[i] + h->histo_drop[i];
    tot += rowtot;
    if ((rowtot > 0) || (i == h->sieve_threshold))
      fprintferr("%s[%3d] %7ld %7ld %7ld %8ld\n",
                 i + 128 == h->sieve_threshold ? "^-" : "  ", i + 128,
                 h->histo_full[i], h->histo_lprl[i], h->histo_drop[i],
                 rowtot);
  }
  fprintferr("        (total evaluated candidates: %ld)\n", tot);
}

/* evaluation/feedback heuristics:
 * First of all, refuse to draw any conclusions unless and until there's
 * enough material to be statistically significant.
 * Second, after sifting through the histo arrays, the new threshold is
 * set to the minimum of the following three quantities:
 * - the position where going down from the top value the histo_full
 *   totals first exceed 100 - MPQS_HISTO_FREL_QUANTILE percent of all
 *   full relations found so far by direct sieving.  (I.e. if the quantile
 *   is 4, we want to keep all rows which account for 96% of all frels
 *   obtained from the sieve.  Note that once we increase the threshold,
 *   further counts will be biased against smaller values;  but we normally
 *   don't expect to do many adjustments.)
 * - the position where, going down from the top towards smaller values,
 *   the cumulative number of useless candidates in histo_drop first exceeds
 *   MPQS_HISTO_DROP_LIMIT times the number of useful ones.  I.e. when that
 *   limit is 2, we're aiming for at least about 1/3 of all candidates coming
 *   from the sieve to result in usable relations.
 * - one less than the position where the histo_lprl count first falls below
 *   MPQS_HISTO_LPREL_BASEFLOW times the number of useless candidates.  This
 *   one will be capable of lowering the current threshold  (but never below
 *   128).
 * XXX For Double Large Prime mode, this will need to be seriously reworked.
 */

/* This function returns 1 when it actually did something and decided on a
 * good threshold  (although possibly the same as it was before),  -1 when
 * there was nothing to do and never will be ("don't call us again"), 0
 * when the caller should retry somewhat later.  Note that mpqs() already
 * knows the total number of candidates generated so far  (from the return
 * values of mpqs_eval_sieve()),  and won't call us too early;  but we also
 * insist on minimal standards for the column sums.  Conversely when we ever
 * lower the threshold, we ask for a re-evaluation later on.
 * NB With the present accounting, once the threshold has been raised, it
 * won't ever be lowered again, since the increasing counts above it will
 * totally swamp the few earlier measurements below which can no longer
 * grow.  So we might chop off those accumulating loops at the current sieve
 * threshold. */
static int
mpqs_eval_histograms(mpqs_handle_t *h)
{
  long tot_full = 0, tot_lprl = 0, tot_drop = 0, total = 0;
  long target_full, i;
  int th_full, th_base, th_drop;
  int th = h->sieve_threshold - 128;

  if (!h->do_histograms) return -1;

  /* first compute column sums */
  for (i = 127; i >= 0; i--)
  {
    tot_full += h->histo_full[i];
    tot_lprl += h->histo_lprl[i];
    tot_drop += h->histo_drop[i];
  }
  total = tot_full + tot_lprl + tot_drop;
  if ((total < MPQS_MIN_CANDS_FOR_HISTO) ||
      (tot_full < MPQS_MIN_FRELS_FOR_HISTO))
    return 0;                   /* too early to call the race */

  th_full = th_drop = th_base = -1;
  /* find the full relations quantile point */
  target_full = tot_full - (tot_full * MPQS_HISTO_FREL_QUANTILE) / 100.;

  tot_full = 0;
  for (i = 127; i >= th; i--)
  {
    if ((tot_full += h->histo_full[i]) >= target_full)
    {
      th_full = i; break;
    }
  }

  /* find the "lp relations baseflow" point */
  for (i = 127; i >= th; i--)
  {
    if (h->histo_lprl[i] + 1 <
        MPQS_HISTO_LPREL_BASEFLOW * h->histo_drop[i])
    {
      th_base = i; break;
    }
  }

  /* find the wastefulness point */
  tot_lprl = 0; tot_drop = 0;
  for (i = 127; i >= th; i--)
  {
    tot_lprl += h->histo_full[i] + h->histo_lprl[i];
    tot_drop += h->histo_drop[i];
    if (tot_drop >
        MPQS_HISTO_DROP_LIMIT * (tot_lprl + 1))
    {
      th_drop = i; break;
    }
  }
  /* if these loops found nothing, then th_(full|base|drop) will still be -1.
   * We won't tighten the sieve, but th_base would tell us we should loosen
   * it  (reluctantly). */

  if (MPQS_DEBUGLEVEL >= 5)
  {
    /* XX in the long run, histo printing to be limited to *very* high
     * XX debug levels, not 5 */
    mpqs_print_histo(h);
    if (th_full >= 0)
      fprintferr("MPQS: threshold estimate for full rels: %d\n",
                 th_full + 128);
    if (th_drop >= 0)
      fprintferr("MPQS: threshold estimate for useful candidates: %d\n",
                 th_drop + 128);
  }

  /* any reason to open up the sieve?  wait until a minimal number of lprels
   * at the threshold has been seen before going down... */
  if ((th > 0) && (th_base <= th) &&
      (h->histo_lprl[th] > (MPQS_MIN_FRELS_FOR_HISTO * 3.5)) )
  {
    h->sieve_threshold = th + 127;
    if (MPQS_DEBUGLEVEL >= 4)
      fprintferr("MPQS: loosening sieve tolerance, new threshold %d\n",
                 h->sieve_threshold);
    return 0;                   /* this should be re-examined after a bit */
  }
  /* otherwise, any reason to tighten it? */
  th = (th_full < th_drop ? th_full : th_drop) + 128;
  if (th > h->sieve_threshold)
  {
    h->sieve_threshold = th;
    if (MPQS_DEBUGLEVEL >= 4)
      fprintferr("MPQS: tightening sieve tolerance, new threshold %d\n",
                 h->sieve_threshold);
  }
  /* maybe also loosen it if th_drop persistently stays below th... */
  return 1;                     /* wait a good while before rechecking */
}
#endif

/*********************************************************************/
/**                                                                 **/
/**           RELATIONS AS STRINGS AND RELATIONS DATABASE           **/
/**                                                                 **/
/*********************************************************************/

/* determines a unique name for a file based on a short nickname
 * name is allocated on the stack */
static char *
mpqs_get_filename(char *dir, char *s)
{
  char *buf = stackmalloc(strlen(dir) + strlen(s) + 2);
#if defined(__EMX__) || defined(WINCE)
  sprintf(buf, "%s\\%s", dir,s);
#else
  sprintf(buf, "%s/%s", dir,s);
#endif
  return buf;
}

/* compares two `large prime' relations according to their first element
 * (the large prime itself). */
static int
mpqs_relations_cmp(const void *a, const void *b)
{
  char **sa = (char**) a;
  char **sb = (char**) b;
  long qa = strtol(*sa, NULL, 10);
  long qb = strtol(*sb, NULL, 10);
  /* atol() isn't entirely portable for the Full Relations case where the
     strings of digits are too long to fit into a long --GN */
  if (qa < qb) return -1;
  else if (qa > qb) return 1;
  else return strcmp(*sa, *sb);
}

static void
pari_fputs(char *s, pariFILE *f)
{
  if (fputs(s, f->file) < 0)
    pari_err(talker, "error whilst writing to file %s", f->name);
}
#define min_bufspace 120UL /* use new buffer when < min_bufspace left */
#define buflist_size 1024  /* size of list-of-buffers blocks */

/* Given a file "filename" containing full or `large prime' relations,
 * rearrange the file so that relations are sorted by their first elements.
 * Works also for sorting full relations. Works in memory, discards duplicate
 * lines, and overwrites the original file. */
static long
mpqs_sort_lp_file(char *filename)
{
  pariFILE *pTMP;
  FILE *TMP;
  char *old_s, *buf, *cur_line;
  char **sort_table, **buflist, **next_buflist, **buflist_head;
  long i, j, count;
  size_t length, bufspace;
  pari_sp av=avma;

  buflist_head = (char**) stackmalloc(buflist_size * sizeof(char*));
  buflist = buflist_head;
  *buflist++ = NULL; /* flag this as last and only buflist block */
  /* extra blocks may be allocated as needed and linked ahead of 
   * buflist_head.  NB: whilst extra buflist blocks might have been
   * needed when we were still sorting entire FREL files (more than 1023
   * buffers, corresponding to about 20000 lines of ~200 characters), they
   * should never be touched now that we only sort LPNEW and FNEW files, which
   * are rather shorter. But the code might as well stay around for future
   * upgrades to handling even larger numbers (and factor bases and thus
   * relations files).  It costs one comparison per buffer allocation. --GN */

  pTMP = pari_fopen(filename, READ);
  TMP = pTMP->file;
  /* get first buffer and read first line, if any, into it */
  buf = (char*) gpmalloc(MPQS_STRING_LENGTH * sizeof(char));
  cur_line = buf;
  bufspace = MPQS_STRING_LENGTH;

  if (fgets(cur_line, bufspace, TMP) == NULL)
  { /* file empty */
    free(buf); pari_fclose(pTMP);
    avma = av; return 0;
  }
  /* enter first buffer into buflist */
  *buflist++ = buf; /* can't overflow the buflist block */
  length = strlen(cur_line) + 1; /* count the \0 byte as well */
  bufspace -= length;

  sort_table = (char**)avma;
  /* at start of loop, one line from the file is sitting in cur_line inside buf,
   * the next will go into cur_line + length, and there's room for bufspace
   * further characters in buf. The loop reads another line if one exists, and
   * if this overruns the current buffer, it allocates a fresh one --GN */
  for (i=0, sort_table--; /* until end of file */; i++, sort_table--)
  { /* sort_table is allocated on the stack, 0x100 cells at a time. Hence the
     * stack must be left alone in the rest of the loop to keep the array
     * connected. In particular, buffers can't be new_chunk'ed */
    if ((i & 0xff) == 0) (void)new_chunk(0x100);
    *sort_table = cur_line;
    cur_line += length;

    /* if little room is left, allocate a fresh buffer before attempting to
     * read a line, and remember to free it if no further line is forthcoming.
     * This avoids some copying of partial lines --GN */
    if (bufspace < min_bufspace)
    {
      if (MPQS_DEBUGLEVEL >= 7)
        fprintferr("MQPS: short of space -- another buffer for sorting\n");
      buf = (char*) gpmalloc(MPQS_STRING_LENGTH * sizeof(char));
      cur_line = buf;
      bufspace = MPQS_STRING_LENGTH;
      if (fgets(cur_line, bufspace, TMP) == NULL) { free(buf); break; }

      /* remember buffer for later deallocation */
      if (buflist - buflist_head >= buflist_size)
      { /* need another buflist block */
        next_buflist = (char**) gpmalloc(buflist_size * sizeof(char*));
        *next_buflist = (char*)buflist_head; /* link */
        buflist_head = next_buflist;
        buflist = buflist_head + 1;
      }
      *buflist++ = buf;
      length = strlen(cur_line) + 1;
      bufspace -= length; continue;
    }

    /* normal case:  try fitting another line into the current buffer */
    if (fgets(cur_line, bufspace, TMP) == NULL) break; /* none exists */
    length = strlen(cur_line) + 1;
    bufspace -= length;

    /* check whether we got the entire line or only part of it */
    if (bufspace == 0 && cur_line[length-2] != '\n')
    {
      size_t lg1;
      if (MPQS_DEBUGLEVEL >= 7)
        fprintferr("MQPS: line wrap -- another buffer for sorting\n");
      buf = (char*) gpmalloc(MPQS_STRING_LENGTH * sizeof(char));
      /* remember buffer for later deallocation */
      if (buflist - buflist_head >= buflist_size)
      { /* need another buflist block */
        next_buflist = (char**)gpmalloc(buflist_size * sizeof(char*));
        *next_buflist = (char*)buflist_head; /* link */
        buflist_head = next_buflist;
        buflist = buflist_head + 1;
      }
      *buflist++ = buf;

      /* copy what we've got to the new buffer */
      (void)strcpy(buf, cur_line); /* cannot overflow */
      cur_line = buf + length - 1; /* point at the \0 byte */
      bufspace = MPQS_STRING_LENGTH - length + 1;
      /* read remainder of line */
      if (fgets(cur_line, bufspace, TMP) == NULL)
        pari_err(talker,"MPQS: relations file truncated?!\n");
      lg1 = strlen(cur_line);
      length += lg1; /* we already counted the \0 once */
      bufspace -= (lg1 + 1); /* but here we must take it into account */
      cur_line = buf; /* back up to the beginning of the line */
    }
  } /* for */

  pari_fclose(pTMP);

  /* sort the whole lot in place by swapping pointers */
  qsort(sort_table, i, sizeof(char*), mpqs_relations_cmp);

  /* copy results back to the original file, skipping exact duplicates */
  pTMP = pari_fopen(filename, WRITE);
  old_s = sort_table[0];
  pari_fputs(sort_table[0], pTMP);
  count = 1;
  for(j = 1; j < i; j++)
  {
    if (strcmp(old_s, sort_table[j]))
    {
      pari_fputs(sort_table[j], pTMP);
      count++;
    }
    old_s = sort_table[j];
  }
  pari_fclose(pTMP);
  if (MPQS_DEBUGLEVEL >= 6) fprintferr("MPQS: done sorting one file.\n");

  /* deallocate buffers and any extraneous buflist blocks except the first */  
  while (*--buflist)
  {
    if (buflist != buflist_head) /* not a linkage pointer */
      free((void*) *buflist);   /* free a buffer */
    else
    { /* linkage pointer */
      next_buflist = (char**)(*buflist);
      free((void*)buflist_head); /* free a buflist block */
      buflist_head = next_buflist;
      buflist = buflist_head + buflist_size;
    }
  }
  avma = av; return count;
}

/* appends contents of file fp1 to f (auxiliary routine for merge sort) and
 * returns number of lines copied. Close f afterwards */
static long
mpqs_append_file(pariFILE *f, FILE *fp1)
{
  FILE *fp = f->file;
  char line[MPQS_STRING_LENGTH];
  long c = 0;
  while (fgets(line, MPQS_STRING_LENGTH, fp1)) { pari_fputs(line, f); c++; }
  if (fflush(fp)) pari_warn(warner, "error whilst flushing file %s", f->name);
  pari_fclose(f); return c;
}

/* Merge-sort on the files LPREL and LPNEW; assumes that LPREL and LPNEW are
 * already sorted. Creates/truncates the TMP file, writes result to it and
 * closes it (via mpqs_append_file()). Instead of LPREL, LPNEW we may also call
 * this with FREL, FNEW. In the latter case pCOMB should be NULL (and we
 * return the count of all full relations), in the former non-NULL (and we
 * return the count of frels we expect to be able to combine out of the
 * present lprels). If pCOMB is non-NULL, the combinable lprels are written
 * out to this separate file.
 * We keep only one occurrence of each `large prime' in TMP (i.e. in the
 * future LPREL file). --GN */

#define swap_lines() { char *line_tmp;\
  line_tmp = line_new_old; \
  line_new_old = line_new; \
  line_new = line_tmp; }

static long
mpqs_mergesort_lp_file0(FILE *LPREL, FILE *LPNEW, pariFILE *pCOMB,
                        pariFILE *pTMP)
{
  char line1[MPQS_STRING_LENGTH], line2[MPQS_STRING_LENGTH];
  char line[MPQS_STRING_LENGTH];
  char *line_new = line1, *line_new_old = line2;
  long q_new, q_new_old = -1, q, i = 0, c = 0;
  long comb_in_progress;

  if ( !fgets(line_new, MPQS_STRING_LENGTH, LPNEW) )
  { /* LPNEW is empty: copy LPREL to TMP. Could be done by a rename if we
     * didn't want to count the lines (again)... however, this case will not
     * normally happen */
    i = mpqs_append_file(pTMP, LPREL);
    return pCOMB ? 0 : i;
  }
  /* we now have a line_new from LPNEW */

  if (!fgets(line, MPQS_STRING_LENGTH, LPREL))
  { /* LPREL is empty: copy LPNEW to TMP... almost. */
    pari_fputs(line_new, pTMP);
    if (!pCOMB)
    { /* full relations mode */
      i = mpqs_append_file(pTMP, LPNEW);
      return i + 1;
    }

    /* LP mode:  check for combinable relations */
    q_new_old = atol(line_new);
    /* we need to retain a copy of the old line just for a moment, because we
     * may yet have to write it to pCOMB. Do this by swapping the two buffers */
    swap_lines();
    comb_in_progress = 0;
    i = 0;

    while (fgets(line_new, MPQS_STRING_LENGTH, LPNEW))
    {
      q_new = atol(line_new);
      if (q_new_old == q_new)
      { /* found combinables, check whether we're already busy on this
           particular `large prime' */
        if (!comb_in_progress)
        { /* if not, write first line to pCOMB, creating and opening the
           * file first if it isn't open yet */
          pari_fputs(line_new_old, pCOMB);
          comb_in_progress = 1;
        }
        /* in any case, write the current line, and count it */
        pari_fputs(line_new, pCOMB);
        i++;
      }
      else
      { /* not combinable */
        q_new_old = q_new;
        comb_in_progress = 0;
        /* and dump it to the TMP file */
        pari_fputs(line_new, pTMP);
        /* and stash it away for a moment */
        swap_lines();
        comb_in_progress = 0;
      }
    } /* while */
    pari_fclose(pTMP); return i;
  }

  /* normal case: both LPNEW and LPREL are not empty */
  q_new = atol(line_new);
  q = atol(line);

  for(;;)
  { /* main merging loop */
    i = comb_in_progress = 0;

    /* first the harder case:  let LPNEW catch up with LPREL, and possibly
       overtake it, checking for combinables coming from LPNEW alone */
    while (q > q_new)
    {
      if (!pCOMB || !comb_in_progress) pari_fputs(line_new, pTMP);
      if (!pCOMB) c++; /* in FREL mode, count lines written */
      else if (!comb_in_progress)
      {
        q_new_old = q_new;
        swap_lines();
      }
      if (!fgets(line_new, MPQS_STRING_LENGTH, LPNEW))
      {
        pari_fputs(line, pTMP);
        if (!pCOMB) c++; else c += i;
        i = mpqs_append_file(pTMP, LPREL);
        return pCOMB? c: c + i;
      }
      q_new = atol(line_new);
      if (!pCOMB) continue;

      /* LP mode only: */
      if (q_new_old != q_new) /* not combinable */
        comb_in_progress = 0; /* next loop will deal with it, or loop may end */
      else
      { /* found combinables, check whether we're already busy on this
           `large prime' */
        if (!comb_in_progress)
        {
          pari_fputs(line_new_old, pCOMB);
          comb_in_progress = 1;
        }
        /* in any case, write the current line, and count it */
        pari_fputs(line_new, pCOMB);
        i++;
      }
    } /* while q > q_new */

    /* q <= q_new */

    if (pCOMB) c += i;   /* accumulate count of combinables */
    i = 0;               /* and clear it */
    comb_in_progress = 0;/* redundant */

    /* now let LPREL catch up with LPNEW, and possibly overtake it */
    while (q < q_new)
    {
      pari_fputs(line, pTMP);
      if (!pCOMB) c++;
      if (!fgets(line, MPQS_STRING_LENGTH, LPREL))
      {
        pari_fputs(line_new, pTMP);
        i = mpqs_append_file(pTMP, LPNEW);
        return pCOMB? c: c + i + 1;
      }
      else
        q = atol(line);
    }

    /* q >= q_new */

    /* Finally, it may happen that q == q_new, indicating combinables whose
     * `large prime' is already in LPREL, and appears now once or more often in
     * LPNEW. Thus in this sub-loop we advance LPNEW. The `line' from LPREL is
     * left alone, and will be written to TMP the next time around the main for
     * loop; we only write it to pCOMB here -- unless all we find is an exact
     * duplicate of the line we already have, that is. (There can be at most
     * one such, and if so it is simply discarded.) */
    while (q == q_new)
    {
      if (!strcmp(line_new, line))
      { /* duplicate -- move right ahead to the next LPNEW line */
        ;/* do nothing here */
      }
      else if (!pCOMB)
      { /* full relations mode: write line_new out first, keep line */
        pari_fputs(line_new, pTMP);
        c++;
      }
      else
      { /* LP mode, and combinable relation */
        if (!comb_in_progress)
        {
          pari_fputs(line, pCOMB);
          comb_in_progress = 1;
        }
        pari_fputs(line_new, pCOMB);
        i++;
      }
      /* NB comb_in_progress is cleared by q_new becoming bigger than q, thus
       * the current while loop terminating, the next time through the main for
       * loop */

      /* common ending: get another line_new, if any */
      if (!fgets(line_new, MPQS_STRING_LENGTH, LPNEW))
      {
        pari_fputs(line, pTMP);
        if (!pCOMB) c++; else c += i;
        i = mpqs_append_file(pTMP, LPREL);
        return pCOMB? c: c + i;
      }
      else
        q_new = atol(line_new);
    } /* while */

    if (pCOMB) c += i; /* accumulate count of combinables */
  }
}

static long
mpqs_mergesort_lp_file(char *REL_str, char *NEW_str, char *TMP_str, pariFILE *pCOMB)
{
  pariFILE *pREL = pari_fopen(REL_str, READ);
  pariFILE *pNEW = pari_fopen(NEW_str, READ);
  pariFILE *pTMP = pari_fopen(TMP_str, WRITE);
  long tp;
  
  tp = mpqs_mergesort_lp_file0(pREL->file, pNEW->file, pCOMB, pTMP);
  pari_fclose(pREL);
  pari_fclose(pNEW);
  pari_unlink(REL_str);
  if (rename(TMP_str,REL_str))
    pari_err(talker, "cannot rename file %s to %s", TMP_str, REL_str);
  if (MPQS_DEBUGLEVEL >= 6)
    fprintferr("MPQS: renamed file %s to %s\n", TMP_str, REL_str);
  return tp;
}


/*********************************************************************/
/**                                                                 **/
/**                       SELF-INITIALIZATION                       **/
/**                                                                 **/
/*********************************************************************/

#ifdef MPQS_DEBUG
/* Debug-only helper routine: check correctness of the root z mod p_i
 * by evaluting A * z^2 + B * z + C mod p_i  (which should be 0).
 * C is written as (B*B-kN)/(4*A) */
static void
check_root(mpqs_handle_t *h, long p, long start)
{
  long z = start - ((long)(h->M) % p);
  if (smodis(addii(h->C, mulsi(z, addii(h->B, mulsi(z, h->A)))), p))
  {
    fprintferr("MPQS: p = %ld\n", p);
    fprintferr("MPQS: A = %Z\n", h->A);
    fprintferr("MPQS: B = %Z\n", h->B);
    fprintferr("MPQS: C = %Z\n", h->C);
    fprintferr("MPQS: z = %ld\n", z);
    pari_err(bugparier, "MPQS: self_init: found wrong polynomial");
  }
}
#endif

/* There are four parts to self-initialization, which are exercised at
 * somewhat different times:
 * - choosing a new coefficient A  (selecting the prime factors to go into it,
 *   and doing the required bookkeeping in the FB entries, including clearing
 *   out the flags from the previous cohort), together with:
 * - doing the actual computations associated with a new A
 * - choosing a new B keeping the same A (very much simpler and quicker)
 * - and a small common bit that needs to happen in both cases.
 * As to the first item, the new scheme works as follows:
 * We pick omega_A - 1 prime factors for A below the index2_FB point which
 * marks their ideal size, and one prime above this point, choosing the
 * latter so as to get log2(A) as close as possible to l2_target_A.
 * The lower prime factors are chosen using bit patterns of constant weight,
 * gradually moving away from index2_FB towards smaller FB subscripts.
 * If this bumps into index0_FB  (might happen for very small input),  we
 * back up by increasing index2_FB by two, and from then on choosing only
 * bit patterns with either or both of their bottom bits set, so at least
 * one of the omega_A - 1 smaller prime factor will be beyond the original
 * index2_FB point.  In this way we avoid re-using A's which had already
 * been done.
 * (The choice of the upper "flyer" prime is of course constrained by the
 * size of the FB, which normally should never be anywhere close to becoming
 * a problem.  In unfortunate cases, e.g. for very small kN, we might have
 * to live with a rather non-optimal choice.
 * Then again, MPQS as such is surprisingly robust.  One day, I had got the
 * order of entries in mpqs_parameterset_t mixed up, and was running on a
 * smallish N with a huuuuge factor base and extremely tiny sieve interval,
 * and it still managed to factor it fairly quickly...)
 *
 * Mathematically, there isn't much more to this than the usual formula for
 * solving a quadratic  (over the field of p elements for each prime p in
 * the FB which doesn't divide A),  solving a linear equation for each of
 * the primes which do divide A, and precomputing differences between roots
 * mod p so we can adjust the roots quickly when we change B.
 * See Thomas Sosnowski's diploma thesis.
 */

/* Helper function:
 * Increment *x (!=0) to a larger value which has the same number of 1s in its
 * binary representation.  Wraparound can be detected by the caller as long as
 * we keep total_no_of_primes_for_A strictly less than BITS_IN_LONG.
 *
 * Changed switch to increment *x in all cases to the next larger number
 * which (a) has the same count of 1 bits and (b) does not arise from the
 * old value by moving a single 1 bit one position to the left  (which was
 * undesirable for the sieve). --GN based on discussion with TP */
INLINE void
mpqs_increment(mpqs_uint32_t *x)
{
  mpqs_uint32_t r1_mask, r01_mask, slider=1UL;

  /* 32-way computed jump handles 22 out of 32 cases */
  switch (*x & 0x1F)
  {
  case 29:
    (*x)++; break; /* shifts a single bit, but we postprocess this case */
  case 26:
    (*x) += 2; break; /* again */
  case 1: case 3: case 6: case 9: case 11:
  case 17: case 19: case 22: case 25: case 27:
    (*x) += 3; return;
  case 20:
    (*x) += 4; break; /* again */
  case 5: case 12: case 14: case 21:
    (*x) += 5; return;
  case 2: case 7: case 13: case 18: case 23:
    (*x) += 6; return;
  case 10:
    (*x) += 7; return;
  case 8:
    (*x) += 8; break; /* and again */
  case 4: case 15:
    (*x) += 12; return;
  default: /* 0, 16, 24, 28, 30, 31 */
    /* isolate rightmost 1 */
    r1_mask = ((*x ^ (*x - 1)) + 1) >> 1;
    /* isolate rightmost 1 which has a 0 to its left */
    r01_mask = ((*x ^ (*x + r1_mask)) + r1_mask) >> 2;
    /* simple cases.  Both of these shift a single bit one position to the
       left, and will need postprocessing */
    if (r1_mask == r01_mask) { *x += r1_mask; break; }
    if (r1_mask == 1) { *x += r01_mask; break; }
    /* general case -- idea: add r01_mask, kill off as many 1 bits as possible
     * to its right while at the same time filling in 1 bits from the LSB. */
    if (r1_mask == 2) { *x += (r01_mask>>1) + 1; return; }
    while (r01_mask > r1_mask && slider < r1_mask)
    {
      r01_mask >>= 1; slider <<= 1;
    }
    *x += r01_mask + slider - 1;
    return;
  }
  /* post-process all cases which couldn't be finalized above.  If we get
     here, slider still has its original value. */
  r1_mask = ((*x ^ (*x - 1)) + 1) >> 1;
  r01_mask = ((*x ^ (*x + r1_mask)) + r1_mask) >> 2;
  if (r1_mask == r01_mask) { *x += r1_mask; return; }
  if (r1_mask == 1) { *x += r01_mask; return; }
  if (r1_mask == 2) { *x += (r01_mask>>1) + 1; return; }
  while (r01_mask > r1_mask && slider < r1_mask)
  {
    r01_mask >>= 1; slider <<= 1;
  }
  *x += r01_mask + slider - 1;
  return;
}

/* self-init (1): advancing the bit pattern, and choice of primes for A.
 * The first time this is called, it finds h->bin_index == 0, which tells us
 * to initialize things from scratch.  On later occasions, we need to begin
 * by clearing the MPQS_FBE_DIVIDES_A bit in the fbe_flags of the former
 * prime factors of A.  We use, of course, the per_A_pr array for finding
 * them.  Upon successful return, that array will have been filled in, and
 * the flag bits will have been turned on again in the right places.
 * We return 1 (true) when we could set things up successfully, and 0 when
 * we found we'd be using more bits to the left in bin_index than we have
 * matching primes for in the FB.  In the latter case, bin_index will be
 * zeroed out, index2_FB will be incremented by 2, index2_moved will be
 * turned on, and the caller, after checking that index2_FB has not become
 * too large, should just call us again, which then is guaranteed to succeed:
 * we'll start again with a right-justified sequence of 1 bits in bin_index,
 * now interpreted as selecting primes relative to the new index2_FB. */
#ifndef MPQS_DEBUG_SI_CHOOSE_PRIMES
#  define MPQS_DEBUG_SI_CHOOSE_PRIMES 0
#endif
INLINE int
mpqs_si_choose_primes(mpqs_handle_t *h)
{
  mpqs_FB_entry_t *FB = h->FB;
  mpqs_per_A_prime_t *per_A_pr = h->per_A_pr;
  double l2_last_p = h->l2_target_A;
  mpqs_int32_t omega_A = h->omega_A;
  int i, j, v2, prev_last_p_idx;
  int room = h->index2_FB - h->index0_FB - omega_A + 4;
  /* GN 20050723:  I.e., index2_FB minus (index0_FB + omega_A - 3) plus 1
   * The notion of room here (cf mpqs_locate_A_range() above) is the number
   * of primes at or below index2_FB which are eligible for A.
   * At the very least, we need omega_A - 1 of them, and it is guaranteed
   * by mpqs_locate_A_range() that at least this many are available when we
   * first get here.  The logic here ensures that the lowest FB slot used
   * for A is never less than index0_FB + omega_A - 3.  In other words, when
   * omega_A == 3 (very small kN), we allow ourselves to reach all the way
   * down to index0_FB;  otherwise, we keep away from it by at least one
   * position.  For omega_A >= 4 this avoids situations where the selection
   * of the smaller primes here has advanced to a lot of very small ones, and
   * the single last larger one has soared away to bump into the top end of
   * the FB. */
  mpqs_uint32_t room_mask;
  mpqs_int32_t p;
  ulong bits;

  /* XXX also clear the index_j field here? */
  if (h->bin_index == 0)
  {
    /* first time here, or after increasing index2_FB, initialize to a pattern
     * of omega_A - 1 consecutive right-justified 1 bits.
     * Caller will have ensured that there are enough primes for this in the
     * FB below index2_FB. */
    h->bin_index = (1UL << (omega_A - 1)) - 1;
    prev_last_p_idx = 0;
  }
  else
  {
    /* clear out the old flags */
    for (i = 0; i < omega_A; i++)
      MPQS_FLG(i) &= ~MPQS_FBE_DIVIDES_A;
    prev_last_p_idx = MPQS_I(omega_A-1);

    /* find out how much maneuvering room we have before we're up against
     * the index0_FB wall */
    if (room > 30) room = 30;
    room_mask = ~((1UL << room) - 1);

    /* bump bin_index to the next acceptable value.  If index2_moved is off,
     * call mpqs_increment() just once;  otherwise, repeat until there's
     * something in the least significant 2 bits - this to ensure that we
     * never re-use an A which we'd used before increasing index2_FB - but
     * also stop if something shows up in the forbidden bits on the left
     * where we'd run out of bits or out of subscripts  (i.e. walk beyond
     * index0_FB + omega_A - 3). */
    mpqs_increment(&h->bin_index);
    if (h->index2_moved)
    {
      while ((h->bin_index & (room_mask | 0x3)) == 0)
        mpqs_increment(&h->bin_index);
      /* XX a slightly simpler increment algorithm would suffice here */
    }
    /* ok so did we fall off the edge on the left? */
    if ((h->bin_index & room_mask) != 0)
    {
      /* Yes.  Turn on the index2_moved flag in the handle */
      h->index2_FB += 2;        /* caller to check this isn't too large!!! */
      h->index2_moved = 1;
      h->bin_index = 0;
      if (MPQS_DEBUG_SI_CHOOSE_PRIMES || (MPQS_DEBUGLEVEL >= 5))
        fprintferr("MPQS: wrapping, more primes for A now chosen near FB[%ld] = %ld\n",
                   (long)h->index2_FB,
                   (long)FB[h->index2_FB].fbe_p);
      return 0;                 /* back off - caller should retry */
    }
  }
  /* assert: we aren't occupying any of the room_mask bits now, and if
   * index2_moved had already been on, at least one of the two LSBs is on */
  bits = h->bin_index;
  if (MPQS_DEBUG_SI_CHOOSE_PRIMES || (MPQS_DEBUGLEVEL >= 6))
    fprintferr("MPQS: new bit pattern for primes for A: 0x%lX\n", bits);

  /* map bits to FB subscripts, counting downward with bit 0 corresponding
   * to index2_FB, and accumulate logarithms against l2_last_p */
  j = h->index2_FB;
  v2 = vals((long)bits);
  if (v2) { j -= v2; bits >>= v2; }
  for (i = omega_A - 2; i >= 0; i--)
  {
    MPQS_I(i) = j;
    l2_last_p -= MPQS_LP(i);
    MPQS_FLG(i) |= MPQS_FBE_DIVIDES_A;
    bits &= ~1UL;
    if (!bits) break;           /* that was the i=0 iteration */
    v2 = vals((long)bits);
    j -= v2;
    bits >>= v2;
  }
  /* Choose the larger prime.  Note we keep index2_FB <= size_of_FB - 3 */
  for (j = h->index2_FB + 1; (p = FB[j].fbe_p) != 0; j++)
  {
    if (FB[j].fbe_flogp > l2_last_p) break;
  }
  /* GN 20050724: The following trick avoids generating a relatively large
   * proportion of duplicate relations when the last prime happens to fall
   * into an area where there are large gaps from one FB prime to the next,
   * and would otherwise often be repeated  (so that successive A's would
   * wind up too similar to each other).  While this trick isn't perfect,
   * it seems to get rid of a major part of the potential duplication. */
  if ((p != 0) && (j == prev_last_p_idx))
  {
    j++; p = FB[j].fbe_p;
  }
  MPQS_I(omega_A - 1) = (p == 0 ? /* did we fall off the end of the FB? */
                         h->size_of_FB + 1 : /* then improvise */
                         j);
  MPQS_FLG(omega_A - 1) |= MPQS_FBE_DIVIDES_A;

  if (MPQS_DEBUG_SI_CHOOSE_PRIMES || (MPQS_DEBUGLEVEL >= 6))
  {
    fprintferr("MPQS: chose primes for A");
    for (i = 0; i < omega_A; i++)
    {
      fprintferr(" FB[%ld]=%ld%s",
                 (long) MPQS_I(i),
                 (long) MPQS_AP(i),
                 i < omega_A - 1 ? "," : "\n");
    }
  }
  return 1;
}



/* compute coefficients of the sieving polynomial for self initializing
 * variant. Coefficients A and B are returned and several tables are
 * updated.   */
/* A and B are kept on the PARI stack in preallocated GENs.  So is C when
 * we're compiled for debugging. */
/* XX total_no_of_primes_for_A to go away */
static void
mpqs_self_init(mpqs_handle_t *h)
{
  const ulong size_of_FB = h->size_of_FB + 1;
  mpqs_FB_entry_t *FB = h->FB;
  mpqs_inv_A_H_t *inv_A_H = h->inv_A_H;
  /* mpqs_uint32_t *bin_index = &(h->bin_index); ---unused */
  const pari_sp av = avma;
  GEN p1, p2;
  GEN A = h->A;
  GEN B = h->B;
#ifdef MPQS_DEBUG
  GEN C = h->C;
#endif
  mpqs_per_A_prime_t *per_A_pr = h->per_A_pr;
  long i, j;
  long inv_A2;
  ulong p;

#ifdef MPQS_DEBUG_AVMA
  fprintferr("MPQS DEBUG: enter self init, avma = 0x%lX\n", (ulong)avma);
#endif
  if (h->index_j == 0)
  /* "mpqs_self_init_A()" */
  { /* compute first polynomial with new A */
    if (!mpqs_si_choose_primes(h))
    {
      /* this means we ran out of room towards small primes, and index2_FB
       * was raised.  Check that we're still ok in that direction before
       * re-trying the operation, which then is guaranteed to succeed.
       * The invariant we maintain towards the top end is that
       * h->size_of_FB - h->index2_FB >= 3, but note that our size_of_FB
       * is one larger. */
      if (size_of_FB - h->index2_FB < 4)
      {
        /* "throw the exception up to caller." We set bin_index to an
         * impossible value and return.  Errm... it has already been set
         * by mpqs_si_choose_primes() to the impossible value 0. */
        return;
      }
      (void) mpqs_si_choose_primes(h);
    }
    /* assert: bin_index and per_A_pr are now populated with consistent
     * values */

    /* compute A = product of omega_A primes given by bin_index */
    p1 = NULL;
    for (i = 0; i < h->omega_A; i++)
    {
      p = (ulong) MPQS_AP(i);
      p1 = p1 ? muliu(p1, p): utoipos(p);
    }
    affii(p1, A); avma = av;

    /* Compute H[i], 0 <= i < omega_A.  Also compute the initial
     * B = sum(v_i*H[i]), by taking all v_i = +1 */
    /* XX following needs to be changed later for segmented FB and sieve
     * XX interval, where we'll want to precompute several B's. */
    p2 = NULL;
    for (i = 0; i < h->omega_A; i++)
    {
      p = (ulong) MPQS_AP(i);
      p1 = divis(A, (long)p);
      p1 = muliu(p1, Fl_inv(umodiu(p1, p), p));
      p1 = muliu(p1, MPQS_SQRT(i));
      affii(remii(p1, A), MPQS_H(i));
      p2 = p2 ? addii(p2, MPQS_H(i)) : MPQS_H(i);
    }
    affii(p2, B);
    avma = av;                  /* must happen outside the loop */

    /* ensure B = 1 mod 4 */
    if (mod2(B) == 0)
      affii(addii(B, mulsi(mod4(A), A)), B); /* B += (A % 4) * A; */

    p1 = shifti(A, 1);
    /* compute the roots z1, z2, of the polynomial Q(x) mod p_j and
     * initialize start1[i] with the first value p_i | Q(z1 + i p_j)
     * initialize start2[i] with the first value p_i | Q(z2 + i p_j)
     * The following loop "miraculously" does The Right Thing for the
     * primes dividing k (where sqrt_kN is 0 mod p).  Primes dividing A
     * are skipped here, and are handled further down in the common part
     * of SI. */
    for (j = 3; (ulong)j <= size_of_FB; j++)
    {
      ulong mb, tmp1, tmp2, m;
      if (FB[j].fbe_flags & MPQS_FBE_DIVIDES_A) continue;
      p = (ulong)FB[j].fbe_p; m = h->M % p;
      inv_A2 = Fl_inv(umodiu(p1, p), p); /* = 1/(2*A) mod p_j */
      mb = umodiu(B, p); if (mb) mb = p - mb;
      /* mb = -B mod p */
      tmp1 = Fl_sub(mb, FB[j].fbe_sqrt_kN, p);
      tmp1 = Fl_mul(tmp1, inv_A2, p);
      FB[j].fbe_start1 = (mpqs_int32_t)Fl_add(tmp1, m, p);

      tmp2 = Fl_add(mb, FB[j].fbe_sqrt_kN, p);
      tmp2 = Fl_mul(tmp2, inv_A2, p);
      FB[j].fbe_start2 = (mpqs_int32_t)Fl_add(tmp2, m, p);
      for (i = 0; i < h->omega_A - 1; i++)
      {
        ulong h = umodiu(MPQS_H(i), p) << 1; if (h > p) h -= p;
        MPQS_INV_A_H(i,j) = Fl_mul(h, inv_A2, p); /* 1/A * H[i] mod p_j */
      }
    }
  }
  else
  /* "mpqs_self_init_B()" */
  { /* no "real" computation -- use recursive formula */
    /* The following exploits that B is the sum of omega_A terms +-H[i].
     * Each time we switch to a new B, we choose a new pattern of signs;
     * the precomputation of the inv_A_H array allows us to change the
     * two arithmetic progressions equally fast.  The choice of sign
     * patterns does *not* follow the bit pattern of the ordinal number
     * of B in the current cohort;  rather, we use a Gray code, changing
     * only one sign each time.  When the i-th rightmost bit of the new
     * ordinal number index_j of B is 1, the sign of H[i] is changed;
     * the next bit to the left tells us whether we should be adding or
     * subtracting the difference term.  We never need to change the sign
     * of H[omega_A-1] (the topmost one), because that would just give us
     * the same sieve items Q(x) again with the opposite sign of x.  This
     * is why we only precomputed inv_A_H up to i = omega_A - 2. */

    ulong v2 = 0;               /* 2-valuation of h->index_j */

    j = h->index_j;
    /* could use vals() here, but we need to right shift the bit pattern
     * anyway in order to find out which inv_A_H entries must be added to or
     * subtracted from the modular roots */
    while ((j & 1) == 0) { v2++; j >>= 1; }

    /* v2 = v_2(index_j), determine new starting positions for sieving */
    p1 = shifti(MPQS_H(v2), 1);
    if (j & 2)
    { /* j = 3 mod 4 */
      for (j = 3; (ulong)j <= size_of_FB; j++)
      {
	if (FB[j].fbe_flags & MPQS_FBE_DIVIDES_A) continue;
        p = (ulong)FB[j].fbe_p;
        FB[j].fbe_start1 = Fl_sub(FB[j].fbe_start1, MPQS_INV_A_H(v2,j), p);
        FB[j].fbe_start2 = Fl_sub(FB[j].fbe_start2, MPQS_INV_A_H(v2,j), p);
      }
      p1 = addii(B, p1);
    }
    else
    { /* j = 1 mod 4 */
      for (j = 3; (ulong)j <= size_of_FB; j++)
      {
	if (FB[j].fbe_flags & MPQS_FBE_DIVIDES_A) continue;
        p = (ulong)FB[j].fbe_p;
        FB[j].fbe_start1 = Fl_add(FB[j].fbe_start1, MPQS_INV_A_H(v2,j), p);
        FB[j].fbe_start2 = Fl_add(FB[j].fbe_start2, MPQS_INV_A_H(v2,j), p);
      }
      p1 = subii(B, p1);
    }
    affii(p1, B);
  }
  avma = av;

  /* p=2 is a special case.  start1[2], start2[2] are never looked at,
   * so don't bother setting them. */

  /* "mpqs_self_init_common()" */

  /* now compute zeros of polynomials that have only one zero mod p
     because p divides the coefficient A */

  /* compute coefficient -C */
  p1 = diviiexact(subii(h->kN, sqri(B)), shifti(A, 2));

  for (i = 0; i < h->omega_A; i++)
  {
    ulong tmp, s;
    p = (ulong) MPQS_AP(i);
    tmp = Fl_div(umodiu(p1, p), umodiu(B, p), p); s = (tmp + h->M) % p;
    FB[MPQS_I(i)].fbe_start1 = (mpqs_int32_t)s;
    FB[MPQS_I(i)].fbe_start2 = (mpqs_int32_t)s;
  }

  if (MPQS_DEBUGLEVEL >= 6)
  {
    /* must happen before resetting avma, because of the absi() */
    fprintferr("MPQS: chose Q_%ld(x) = %Z x^2 %c %Z x + C\n",
               (long) h->index_j, h->A,
               signe(h->B) < 0? '-': '+', absi(h->B));
  }

  avma = av;

#ifdef MPQS_DEBUG
  /* stash C into the handle.  Since check_root() is the only thing which
   * uses it, and only for debugging, C doesn't even exist as a field in
   * the handle unless we're built with MPQS_DEBUG. */
  affii(negi(p1), h->C);
  for (j = 3; j <= size_of_FB; j++)
  {
    check_root(h, FB[j].fbe_p, FB[j].fbe_start1);
    check_root(h, FB[j].fbe_p, FB[j].fbe_start2); avma = av;
  }
  if (DEBUGLEVEL >= 6)
    PRINT_IF_VERBOSE("MPQS: checking of roots of Q(x) was successful\n");
#endif

#ifdef MPQS_DEBUG_AVMA
  fprintferr("MPQS DEBUG: leave self init, avma = 0x%lX\n", (ulong)avma);
#endif
}

/*********************************************************************/
/**                                                                 **/
/**                           THE SIEVE                             **/
/**                                                                 **/
/*********************************************************************/

/* Main sieving routine:
 *
 * XXX purge the following obsolescent comment...
 * FB: pointer to an array which holds the factor base
 * log_FB: pointer to an array which holds the approximations for
 *        the logarithms of the factor base
 * start_1, start_2: arrays for starting points for different FB elements
 * sieve_array: points to a sieve array
 * sieve_array_end: points to the end of the sieve array
 * M: size of the sieving interval
 * starting_sieving_index: marks the first FB element used for sieving */
INLINE void
mpqs_sieve_p(unsigned char *begin, unsigned char *end,
             long p4, long p, unsigned char log_p)
{
  register unsigned char *e = end - p4;
  /* XX I unrolled this loop some 5 or 6 years ago, but these days I tend
   * XX to believe it's better to let the compiler worry about *this* kind
   * XX of optimization, based on its built-in knowledge of whatever useful
   * XX tricks the machine instruction set architecture may be offering
   * XX ("speculative loads" being the buzzword). --GN */
  while (e - begin >= 0) /* signed comparison */
  {
    (*begin) += log_p, begin += p;
    (*begin) += log_p, begin += p;
    (*begin) += log_p, begin += p;
    (*begin) += log_p, begin += p;
  }
  while (end - begin >= 0) /* again */
    (*begin) += log_p, begin += p;
}

static void
mpqs_sieve(mpqs_handle_t *h)
{
  long p, start1, start2, l = h->index1_FB;
  unsigned char log_p;
  mpqs_FB_entry_t *ptr_FB;
  unsigned char *sieve_array = h->sieve_array;
  unsigned char *sieve_array_end = h->sieve_array_end;

  for (ptr_FB = &(h->FB[l]); (p = ptr_FB->fbe_p) != 0; ptr_FB++, l++)
  {
    /* XX get rid of redundant local variables */
    log_p = ptr_FB->fbe_logval;
    start1 = ptr_FB->fbe_start1;
    start2 = ptr_FB->fbe_start2;

    /* sieve with FB[l] from start_1[l] */
    /* if start1 != start2 sieve with FB[l] from start_2[l] */
    /* XX maybe it is more efficient not to have a conditional branch in
     * XX the present loop body, and instead to right-shift log_p one bit
     * XX based on a flag bit telling us that we're on a one-root prime?
     * XX And instead roll the two invocations of mpqs_sieve_p into one. */
    mpqs_sieve_p(sieve_array + start1, sieve_array_end, p << 2, p, log_p);
    if (start1 != start2)
      mpqs_sieve_p(sieve_array + start2, sieve_array_end, p << 2, p, log_p);
  }
}

/******************************/

/* Could make shameless use of the fact that M is divisible by 4, but
 * let the compiler worry about loop unrolling.  Indeed I wonder whether
 * modern compilers woudln't have an easier time optimizing this if it
 * were written as array accesses.  Doing so. */
static long
mpqs_eval_sieve(mpqs_handle_t *h)
{
  long x = 0, count = 0, M_2 = h->M << 1;
  /* XX Todo: replace the following by an auto-adjusting threshold driven
   * XX by histogram yield measurements */
  unsigned char th = h->sieve_threshold;
  unsigned char *sieve_array = h->sieve_array;
  long *candidates = h->candidates;

  /* The following variation on the original is marginally faster with a
   * good optimizing compiler.  Exploiting the sentinel, we don't need to
   * check for x < M_2 in the inner while loop - this more than makes up
   * for the "lack" of explicit unrolling.  Optimizations like loop
   * unrolling are best left to the compiler anyway... */
  while (count < MPQS_CANDIDATE_ARRAY_SIZE - 1)
  {
    while (sieve_array[x] < th) x++;
    if (x >= M_2) break;
    candidates[count++] = x++;
  }

  candidates[count] = 0; return count;
}

/*********************************************************************/
/**                                                                 **/
/**                     CONSTRUCTING RELATIONS                      **/
/**                                                                 **/
/*********************************************************************/



/* Main relation routine:
 *
 * XXX fix obsolescent comment...
 * FB: pointer to an array which holds the factor base
 * start_1, start_2: arrays for starting points for different FB elements
 * M: size of the sieving interval */

static void
mpqs_add_factor(char **last, ulong ei, ulong pi) {
  sprintf(*last, " %lu %lu", ei, pi);
  *last += strlen(*last);
}

/* concatenate " 0" */
static void
mpqs_add_0(char **last) {
  char *s = *last;
  *s++ = ' ';
  *s++ = '0';
  *s++ = 0; *last = s;
}

#ifdef MPQS_DEBUG
/* XXX use this in *all* places which factor something back...? */
static GEN
mpqs_factorback(mpqs_handle_t *h, char *relations)
{
  char *s, *t = pari_strdup(relations);
  GEN p_e, N = h->N, prod = gen_1;
  long e;
  mpqs_FB_entry_t *FB = h->FB;

  s = strtok(t, " \n");
  while (s != NULL)
  {
    e = atol(s); if (!e) break;
    s = strtok(NULL, " \n");
    if (atol(s) == 1)           /* special case -1... */
    { prod = subii(N, prod); s = strtok(NULL, " \n"); continue; }
    p_e = Fp_powu(utoipos(FB[atol(s)].fbe_p), (ulong)e, N);
    prod = remii(mulii(prod, p_e), N);
    s = strtok(NULL, " \n");
  }
  free(t); return prod;
}
#endif

/* NB FREL, LPREL are actually FNEW, LPNEW when we get called */
static long
mpqs_eval_cand(mpqs_handle_t *h, long number_of_cand,
               FILE *FREL, FILE *LPREL)
{
  pari_sp av;
  long number_of_relations = 0;
  char *relations = h->relations;
  long *relaprimes = h->relaprimes;
  ulong i, pi;
  mpqs_FB_entry_t *FB = h->FB;
  GEN A = h->A;
  GEN B = h->B;                 /* we don't need coefficient C here */
  int pii;
  long *candidates = h->candidates;

  av = avma;
#ifdef MPQS_DEBUG_AVMA
  fprintferr("MPQS DEBUG: enter eval cand, avma = 0x%lX\n", (ulong)avma);
#endif
  for (i = 0; i < (ulong)number_of_cand; i++, avma = av)
  {
    GEN Qx, Qx_part, A_2x_plus_B, Y;
    long powers_of_2, p;
    long x = candidates[i];
    long x_minus_M = x - h->M;
    char *relations_end = relations;
    int relaprpos = 0;

#ifdef MPQS_DEBUG_AVMA
    fprintferr("MPQS DEBUG: eval loop 1, avma = 0x%lX\n", (ulong)avma);
#endif

    *relations_end = 0;
#ifdef MPQS_DEBUG_VERYVERBOSE
    fprintferr("%c", (char)('0' + i%10));
#endif

    /* A_2x_plus_B = (A*(2x)+B), Qx = (A*(2x)+B)^2/(4*A) = Q(x) */
    A_2x_plus_B = addii(mulis(A, x_minus_M << 1), B);
    Y = absi(A_2x_plus_B);

    Qx = subii(sqri(A_2x_plus_B), h->kN);

#ifdef MPQS_DEBUG_AVMA
    fprintferr("MPQS DEBUG: eval loop 2, avma = 0x%lX\n", (ulong)avma);
#endif
    /* When N is relatively small, it may happen that Qx is outright
     * divisible by N at this point.  In any case, when no extensive prior
     * trial division / Rho / ECM had been attempted, gcd(Qx,N) may turn
     * out to be a nontrivial factor of N  (larger than what the FB contains
     * or we'd have found it already, but possibly smaller than the large-
     * prime bound).  This is too rare to check for here in the inner loop,
     * but it will be caught if such an LP relation is ever combined with
     * another. */

    /* XXX Qx cannot possibly vanish here. */
    if (!signe(Qx)) { PRINT_IF_VERBOSE("<+>"); continue; }
    else if (signe(Qx) < 0) {
      setsigne(Qx, 1);
      mpqs_add_factor(&relations_end, 1, 1); /* i = 1, ei = 1, pi */
    }

    /* divide by powers of 2;  we're really dealing with 4*A*Q(x), so we
     * always have at least 2^2 here, and at least 2^3 when kN is 1 mod 4 */
    powers_of_2 = vali(Qx);
    Qx = shifti(Qx, -powers_of_2);
    mpqs_add_factor(&relations_end, powers_of_2, 2);

    /* That has dealt with a possible -1 and the power of 2.  First pass
     * over odd primes in FB: pick up all possible divisors of Qx including
     * those sitting in k or in A, and remember them in relaprimes. Do not
     * yet worry about possible repeated factors, these will be found in the
     * second pass. */
    Qx_part = A;

    /* The first pass recognizes divisors of A by their corresponding flags
     * bit in the FB entry.  (Divisors of k require no special treatment at
     * this stage.)  We construct a preliminary table of FB subscripts and
     * "exponents" of the FB primes which divide Qx.  (We store subscripts
     * rather than the primes themselves because the string representation
     * of a relation is in terms of the subscripts.)
     * We must distinguish three cases so we can do the right thing in the
     * 2nd pass: prime not in A which divides Qx, prime in A which does not
     * divide Qx/A, prime in A which does divide Qx/A. The first and third
     * kinds need checking for repeated factors, the second kind doesn't. The
     * first and second kinds contribute 1 to the exponent in the relation,
     * the 3rd kind contributes 2. We store 1,0,2 respectively in these three
     * cases.
     * Factors in common with k are much simpler - if they occur, they occur
     * exactly to the first power, and this makes no difference in the first
     * pass - here they behave just like every normal odd factor base prime.
     */

    for (pi = 3; (p = FB[pi].fbe_p); pi++)
    {
      long tmp_p = x % p;
      ulong ei = 0;

      /* Here we use that MPQS_FBE_DIVIDES_A equals 1. */
      ei = FB[pi].fbe_flags & MPQS_FBE_DIVIDES_A;

      if (tmp_p == FB[pi].fbe_start1 || tmp_p == FB[pi].fbe_start2)
      { /* p divides Q(x)/A (and possibly A), 1st or 3rd case */
        relaprimes[relaprpos++] = pi;
        relaprimes[relaprpos++] = 1 + ei;
        Qx_part = mulis(Qx_part, p);
      }
      else if (ei)
      { /* p divides A but does not divide Q(x)/A, 2nd case */
        relaprimes[relaprpos++] = pi;
        relaprimes[relaprpos++] = 0;
      }
    }

    /* We have now accumulated the known factors of Qx except for possible
     * repeated factors and for possible large primes.  Divide off what we
     * have.  (This is faster than dividing off A and each prime separately.)
     */
    Qx = diviiexact(Qx, Qx_part);
    /* (ToDo: MPQS_DEBUG sanity check...) */

#ifdef MPQS_DEBUG_AVMA
    fprintferr("MPQS DEBUG: eval loop 3, avma = 0x%lX\n", (ulong)avma);
#endif

    /* second pass - deal with any repeated factors, and write out the string
     * representation of the tentative relation. At this point, the only
     * primes which can occur again in the adjusted Qx are those in relaprimes
     * which are followed by 1 or 2.  We must pick up those followed by a 0,
     * too, though. */
    PRINT_IF_VERBOSE("a");
    for (pii = 0; pii < relaprpos; pii+=2)
    {
      long remd_p;
      ulong ei = relaprimes[pii+1];
      GEN Qx_div_p;
      pi = relaprimes[pii];

      /* Here, prime factors of k go their separate way.  We could have
       * introduced another FB entry flag for marking them, but it is easier
       * to identify them just by their position before index0_FB. */
      if ((mpqs_int32_t)pi < h->index0_FB) {
#ifdef MPQS_DEBUG
        PRINT_IF_VERBOSE("\bk!");
#endif
        mpqs_add_factor(&relations_end, 1, pi);
        continue;
      }
      if (ei == 0) /* p divides A and that was it */
      {
        mpqs_add_factor(&relations_end, 1, pi);
        continue;
      }
      p = FB[pi].fbe_p;
#ifdef MPQS_DEBUG_CANDIDATE_EVALUATION
      fprintferr("MPQS DEBUG: Qx=%Z p=%ld\n", Qx, (long)p);
#endif
      /* otherwise p might still divide the current adjusted Qx. Try it... */
      /* XXX break out of loop when remaining Qx is 1.  Or rather, suppress
       * the trial divisions, since we still need to write our string.
       * Actually instead of testing for 1, test whether Qx is smaller than
       * p;  cf Karim's mail from 20050124.  If it is, without being 1,
       * then it has a common factor with k.  But those factors are soon
       * going to have disappeared before we get here.  However, inserting
       * an explicit if (!is_pm1(Qx)) here did not help any. */
      Qx_div_p = divis_rem(Qx, p, &remd_p);
      while (remd_p == 0) {
        ei++; Qx = Qx_div_p;
        Qx_div_p = divis_rem(Qx, p, &remd_p);
      }
      mpqs_add_factor(&relations_end, ei, pi);
    }

#ifdef MPQS_DEBUG_AVMA
    fprintferr("MPQS DEBUG: eval loop 4, avma = 0x%lX\n", (ulong)avma);
#endif
    PRINT_IF_VERBOSE("\bb");
    if (is_pm1(Qx))
    {
      mpqs_add_0(&relations_end);
      fprintf(FREL, "%s :%s\n", i2str(Y), relations);
      number_of_relations++;
#ifdef MPQS_USE_HISTOGRAMS
      /* bump full relations counter at candidate's value */
      if (h->do_histograms) h->histo_full[sa[x]-128]++;
#endif

#ifdef MPQS_DEBUG
      {
        pari_sp av1 = avma;
        GEN rhs = mpqs_factorback(h, relations);
        GEN Qx_2 = remii(sqri(Y), h->N);
        if (!equalii(Qx_2, rhs))
        {
          PRINT_IF_VERBOSE("\b(!)\n");
          fprintferr("MPQS: %Z @ %Z :%s\n", Y, Qx, relations);
          fprintferr("\tQx_2 = %Z\n", Qx_2);
          fprintferr("\t rhs = %Z\n", rhs);
          pari_err(talker, "MPQS: wrong full relation found!!");
        }
        else
          PRINT_IF_VERBOSE("\b(:)");
        avma = av1;
      }
#endif
    }
    else if (cmpis(Qx, h->lp_bound) > 0)
    { /* TODO: check for double large prime */
#ifdef MPQS_USE_HISTOGRAMS
      /* bump useless-candidates counter at candidate's value */
      if (h->do_histograms) h->histo_drop[sa[x]-128]++;
#endif
      PRINT_IF_VERBOSE("\b.");
    }
    else
    { /* if (mpqs_isprime(itos(Qx))) */
      mpqs_add_0(&relations_end);
      fprintf(LPREL, "%s @ %s :%s\n", i2str(Qx), i2str(Y), relations);
#ifdef MPQS_USE_HISTOGRAMS
      /* bump LP relations counter at candidate's value */
      if (h->do_histograms) h->histo_lprl[sa[x]-128]++;
#endif
#ifdef MPQS_DEBUG
      {
        pari_sp av1 = avma;
        GEN rhs = mpqs_factorback(h, relations);
        GEN Qx_2 = remii(sqri(Y), h->N);

        rhs = modii(mulii(rhs, Qx), h->N);
        if (!equalii(Qx_2, rhs))
        {
          PRINT_IF_VERBOSE("\b(!)\n");
          fprintferr("MPQS: %Z @ %Z :%s\n", Y, Qx, relations);
          fprintferr("\tQx_2 = %Z\n", Qx_2);
          fprintferr("\t rhs = %Z\n", rhs);
          pari_err(talker, "MPQS: wrong large prime relation found!!");
        }
        else
          PRINT_IF_VERBOSE("\b(;)");
        avma = av1;
      }
#endif
    }

#ifdef MPQS_DEBUG_AVMA
    fprintferr("MPQS DEBUG: eval loop end, avma = 0x%lX\n", (ulong)avma);
#endif
  } /* for */
  PRINT_IF_VERBOSE("\n");
#ifdef MPQS_DEBUG_AVMA
  fprintferr("MPQS DEBUG: leave eval cand, avma = 0x%lX\n", (ulong)avma);
#endif
  return number_of_relations;
}

/*********************************************************************/
/**                                                                 **/
/**                      COMBINING RELATIONS                        **/
/**                                                                 **/
/*********************************************************************/

/* combines the large prime relations in COMB to full relations in FNEW.
 * FNEW is assumed to be open for writing / appending. */

typedef struct {
  long q;
  char Y[MPQS_STRING_LENGTH];
  char E[MPQS_STRING_LENGTH];
} mpqs_lp_entry;

static void
mpqs_set_exponents(long *ei, char *r)
{
  char *s, b[MPQS_STRING_LENGTH];
  long e;

  strcpy(b, r);
  s = strtok(b, " \n");
  while (s != NULL)
  {
    e = atol(s); if (!e) break;
    s = strtok(NULL, " \n");
    ei[ atol(s) ] += e;
    s = strtok(NULL, " \n");
  }
}

static void
set_lp_entry(mpqs_lp_entry *e, char *buf)
{
  char *s1, *s2;
  s1 = buf; s2 = strchr(s1, ' '); *s2 = '\0';
  e->q = atol(s1);
  s1 = s2 + 3; s2 = strchr(s1, ' '); *s2 = '\0';
  strcpy(e->Y, s1);
  s1 = s2 + 3; s2 = strchr(s1, '\n'); *s2 = '\0';
  strcpy(e->E, s1);
}

static long
mpqs_combine_large_primes(mpqs_handle_t *h,
                          FILE *COMB, pariFILE *pFNEW, GEN *f)
{
  pari_sp av0 = avma, av, av2;
  char new_relation[MPQS_STRING_LENGTH], buf[MPQS_STRING_LENGTH];
  mpqs_lp_entry e[2]; /* we'll use the two alternatingly */
  long *ei, ei_size = h->size_of_FB + 2;
  long old_q;
  GEN inv_q, Y1, Y2, new_Y, new_Y1;
  long i, l, c = 0;
#ifdef MPQS_DEBUG
  mpqs_FB_entry_t *FB = h->FB;  /* for factoring back */
#endif

  *f = NULL;
  if (!fgets(buf, MPQS_STRING_LENGTH, COMB)) return 0; /* should not happen */

  ei = (long *) new_chunk(ei_size);
  av = avma;
  /* put first lp relation in row 0 of e */
  set_lp_entry(&e[0], buf);

  i = 1; /* second relation will go into row 1 */
  old_q = e[0].q;
  while (!invmod(utoipos(old_q), h->N, &inv_q)) /* can happen */
  {
    inv_q = gcdii(inv_q, h->N);
    /* inv_q can no longer be 1 here (it could while we were doing this mod
     * kN instead of mod N), but never mind - we're not in the fast path
     * at this point.  It could be N when N is quite small;  or we might
     * just have found a divisor by sheer luck. */
    if (is_pm1(inv_q) || equalii(inv_q, h->N)) /* pity */
    {
#ifdef MPQS_DEBUG
      fprintferr("MPQS: skipping relation with non-invertible q\n");
#endif
      if (!fgets(buf, MPQS_STRING_LENGTH, COMB)) { avma = av0; return 0; }
      avma = av;
      set_lp_entry(&e[0], buf);
      old_q = e[0].q; continue;
    }
    *f = gerepileuptoint(av0, inv_q);
    return c;
  }
  Y1 = strtoi(e[0].Y);
  av2 = avma; /* preserve inv_q and Y1 */

  while (fgets(buf, MPQS_STRING_LENGTH, COMB))
  {
    set_lp_entry(&e[i], buf);
    if (e[i].q != old_q)
    {
      /* switch to combining a new bunch, swapping the rows */
      old_q = e[i].q;
      avma = av; /* discard old inv_q and Y1 */
      if (!invmod(utoipos(old_q), h->N, &inv_q)) /* can happen --GN */
      {
        inv_q = gcdii(inv_q, h->N);
        if (is_pm1(inv_q) || equalii(inv_q, h->N)) /* pity */
        {
#ifdef MPQS_DEBUG
          fprintferr("MPQS: skipping relation with non-invertible q\n");
#endif
          old_q = -1; /* sentinel */
          av2 = avma = av;
          continue; /* discard this combination */
        }
        *f = gerepileuptoint(av0, inv_q);
        return c;
      }
      Y1 = strtoi(e[i].Y);
      i = 1 - i; /* subsequent relations go to other row */
      av2 = avma; /* preserve inv_q and Y1 */
      continue;
    }
    /* count and combine the two we've got, and continue in the same row */
    c++;
    memset((void *)ei, 0, ei_size * sizeof(long));
    mpqs_set_exponents(ei, e[0].E);
    mpqs_set_exponents(ei, e[1].E);
    Y2 = strtoi(e[i].Y);
    new_Y = modii(mulii(mulii(Y1, Y2), inv_q), h->N);
    new_Y1 = subii(h->N, new_Y);
    if (absi_cmp(new_Y1, new_Y) < 0) new_Y = new_Y1;
    strcpy(new_relation, i2str(new_Y));
    strcat(new_relation, " :");
    if (ei[1] & 1) strcat(new_relation, " 1 1");
    for (l = 2; l < ei_size; l++)
      if (ei[l])
      {
        sprintf(buf, " %ld %ld", ei[l], l);
        strcat(new_relation, buf);
      }
    strcat(new_relation, " 0");
    if (DEBUGLEVEL >= 6)
    {
      fprintferr("MPQS: combining\n");
      fprintferr("    {%ld @ %s : %s}\n", old_q, e[1-i].Y, e[1-i].E);
      fprintferr("  * {%ld @ %s : %s}\n", e[i].q, e[i].Y, e[i].E);
      fprintferr(" == {%s}\n", new_relation);
    }
    strcat(new_relation, "\n");

#ifdef MPQS_DEBUG
    {
      char ejk [MPQS_STRING_LENGTH];
      GEN Qx_2, prod_pi_ei, pi_ei;
      char *s;
      long pi, exi;
      pari_sp av1 = avma;
      Qx_2 = modii(sqri(new_Y), h->N);

      strcpy(ejk, new_relation);
      s = strchr(ejk, ':') + 2;
      s = strtok(s, " \n");

      /* XXX Why isn't this re-using mpqs_factorback()? */
      /* XXX And de-mess the messy ad-hoc identifiers... */
      prod_pi_ei = gen_1;
      while (s != NULL)
      {
        exi = atol(s); if (!exi) break;
        s = strtok(NULL, " \n");
        pi = atol(s);
        if (pi == 1)            /* special case -1... */
        {
          prod_pi_ei = subii(h->N, prod_pi_ei); s = strtok(NULL, " \n");
          continue;
        }
        pi_ei = Fp_pow(utoipos(FB[pi].fbe_p), utoipos(exi), h->N);
        prod_pi_ei = modii(mulii(prod_pi_ei, pi_ei), h->N);
        s = strtok(NULL, " \n");
      }
      avma = av1;

      if (!equalii(Qx_2, prod_pi_ei))
        pari_err(talker, "MPQS: combined large prime relation is false");
    }
#endif

    pari_fputs(new_relation, pFNEW);
    avma = av2;
  } /* while */

  if (DEBUGLEVEL >= 4)
    fprintferr("MPQS: combined %ld full relation%s\n", c, (c!=1 ? "s" : ""));
  avma = av0; return c;
}

/*********************************************************************/
/**                                                                 **/
/**                    GAUSS-LANCZOS ELIMINATION                    **/
/**                                                                 **/
/*********************************************************************/

#ifdef LONG_IS_64BIT

#define MPQS_GAUSS_BITS 64
static unsigned long mpqs_mask_bit[]  =
{
  0x8000000000000000UL, 0x4000000000000000UL,
  0x2000000000000000UL, 0x1000000000000000UL,
  0x0800000000000000UL, 0x0400000000000000UL,
  0x0200000000000000UL, 0x0100000000000000UL,
  0x0080000000000000UL, 0x0040000000000000UL,
  0x0020000000000000UL, 0x0010000000000000UL,
  0x0008000000000000UL, 0x0004000000000000UL,
  0x0002000000000000UL, 0x0001000000000000UL,
  0x0000800000000000UL, 0x0000400000000000UL,
  0x0000200000000000UL, 0x0000100000000000UL,
  0x0000080000000000UL, 0x0000040000000000UL,
  0x0000020000000000UL, 0x0000010000000000UL,
  0x0000008000000000UL, 0x0000004000000000UL,
  0x0000002000000000UL, 0x0000001000000000UL,
  0x0000000800000000UL, 0x0000000400000000UL,
  0x0000000200000000UL, 0x0000000100000000UL,
  0x0000000080000000UL, 0x0000000040000000UL,
  0x0000000020000000UL, 0x0000000010000000UL,
  0x0000000008000000UL, 0x0000000004000000UL,
  0x0000000002000000UL, 0x0000000001000000UL,
  0x0000000000800000UL, 0x0000000000400000UL,
  0x0000000000200000UL, 0x0000000000100000UL,
  0x0000000000080000UL, 0x0000000000040000UL,
  0x0000000000020000UL, 0x0000000000010000UL,
  0x0000000000008000UL, 0x0000000000004000UL,
  0x0000000000002000UL, 0x0000000000001000UL,
  0x0000000000000800UL, 0x0000000000000400UL,
  0x0000000000000200UL, 0x0000000000000100UL,
  0x0000000000000080UL, 0x0000000000000040UL,
  0x0000000000000020UL, 0x0000000000000010UL,
  0x0000000000000008UL, 0x0000000000000004UL,
  0x0000000000000002UL, 0x0000000000000001UL
};

#else

#define MPQS_GAUSS_BITS 32
static unsigned long mpqs_mask_bit[]  =
{
  0x80000000UL, 0x40000000UL, 0x20000000UL, 0x10000000UL,
  0x08000000UL, 0x04000000UL, 0x02000000UL, 0x01000000UL,
  0x00800000UL, 0x00400000UL, 0x00200000UL, 0x00100000UL,
  0x00080000UL, 0x00040000UL, 0x00020000UL, 0x00010000UL,
  0x00008000UL, 0x00004000UL, 0x00002000UL, 0x00001000UL,
  0x00000800UL, 0x00000400UL, 0x00000200UL, 0x00000100UL,
  0x00000080UL, 0x00000040UL, 0x00000020UL, 0x00000010UL,
  0x00000008UL, 0x00000004UL, 0x00000002UL, 0x00000001UL
};

#endif

static F2_matrix
F2_create_matrix(long rows, long cols)
{
  F2_matrix m;
  long i, j, words = cols / MPQS_GAUSS_BITS;
  if (cols % MPQS_GAUSS_BITS) words++;
  m = (F2_row *) gpmalloc(rows * sizeof(F2_row));
  for (i = 0; i < rows; i++)
  {
    m[i] = (ulong *) gpmalloc(words * sizeof(ulong));
    for (j = 0; j < words; j++) m[i][j] = 0UL;
  }
  return m;
}

static void
F2_destroy_matrix(F2_matrix m, long rows)
{
  long i;
  for (i = 0; i < rows; i++) free(m[i]);
  free(m);
}

static ulong
F2_get_bit(F2_matrix m, long i, long j)
{
  return m[i][j / MPQS_GAUSS_BITS] & mpqs_mask_bit[j % MPQS_GAUSS_BITS];
}
static void
F2_set_bit(F2_matrix m, long i, long j)
{
  m[i][j / MPQS_GAUSS_BITS] |= mpqs_mask_bit[j % MPQS_GAUSS_BITS];
}
#if 0
static void
F2_clear_bit(F2_matrix m, long i, long j)
{
  m[i][j / MPQS_GAUSS_BITS] &= ~mpqs_mask_bit[j % MPQS_GAUSS_BITS];
}
#endif

/* output an F2_matrix in PARI format */
static void
F2_print_matrix(F2_matrix m, long rows, long cols)
{
  long i, j;
  fprintferr("\n[");
  for (i = 0; i < rows; i++)
  {
    for (j = 0; j < cols - 1; j++)
      fprintferr( F2_get_bit(m, i, j)? "1, ": "0, " );
    fprintferr( F2_get_bit(m, i, j)? "1": "0" );
    if (i != rows - 1) fprintferr("; ");
  }
  fprintferr("]\n");
}

/* x ^= y : row addition over F_2 */
static void
F2_add_rows(F2_row y, F2_row x, long k, long n)
{
  long i, q, r;
  n = n - k;
  r = n % 8; q = n - r + k; i = 0 + k;
  for (; i < q; i += 8)
  {
    x[  i] ^= y[  i]; x[1+i] ^= y[1+i]; x[2+i] ^= y[2+i]; x[3+i] ^= y[3+i];
    x[4+i] ^= y[4+i]; x[5+i] ^= y[5+i]; x[6+i] ^= y[6+i]; x[7+i] ^= y[7+i];
  }
  switch (r)
  {
    case 7: x[i] ^= y[i]; i++;
    case 6: x[i] ^= y[i]; i++;
    case 5: x[i] ^= y[i]; i++;
    case 4: x[i] ^= y[i]; i++;
    case 3: x[i] ^= y[i]; i++;
    case 2: x[i] ^= y[i]; i++;
    case 1: x[i] ^= y[i]; i++;
  }
}

/* create and read an F2_matrix from a relation file FREL (opened by caller).
 * Also record the position of each relation in the file for later use
 * rows = size_of_FB+1, cols = rel */
static F2_matrix
F2_read_matrix(FILE *FREL, long rows, long cols, long *fpos)
{
  F2_matrix m = F2_create_matrix(rows, cols);
  long i = 0, e, p;
  char buf[MPQS_STRING_LENGTH], *s;

  if ((fpos[0] = ftell(FREL)) < 0)
    pari_err(talker, "ftell error on full relations file");
  while (fgets(buf, MPQS_STRING_LENGTH, FREL))
  {
    s = strchr(buf, ':') + 2;
    s = strtok(s, " \n");
    while (s != NULL)
    {
      e = atol(s); if (!e) break;
      s = strtok(NULL, " \n");
      p = atol(s);
      if (e & 1) F2_set_bit(m, p - 1, i);
      s = strtok(NULL, " \n");
    }
    i++;
    if (i < cols && (fpos[i] = ftell(FREL)) < 0)
      pari_err(talker, "ftell error on full relations file");
  }
  if (i != cols)
  {
    fprintferr("MPQS: full relations file %s than expected",
               i > cols ? "longer" : "shorter");
    pari_err(talker, "MPQS panicking");
  }
  return m;
}

/* compute the kernel of an F2_matrix over F_2, m = rows, n = cols */
static F2_matrix
mpqs_kernel(F2_matrix x, long m, long n, long *rank)
{
  pari_sp av = avma;
  GEN c, d;
  long k, i, j, t, r = 0;
  F2_matrix ker_x;
  long words = n / MPQS_GAUSS_BITS;
  if (n % MPQS_GAUSS_BITS) words++;

  d = new_chunk(n);
  c = new_chunk(m);
  for (k = 0; k < m; k++) c[k] = -1;

  for (k = 0; k < n; k++)
  {
    j = 0;
    while (j < m && (c[j] >= 0 || F2_get_bit(x, j, k) == 0)) j++;

    if (j == m)
    { /* no pivot found in column k: it's a kernel vector */
      d[k] = -1; r++;
    }
    else
    {
      d[k] = j; /* the pivot in column k is at row j */
      c[j] = k; /* a pivot at position j was found in column k */
      for (t = 0; t < m; t++)
        if (t != j && F2_get_bit(x, t, k))
          F2_add_rows(x[j], x[t], k / MPQS_GAUSS_BITS, words);
    }
  }

  ker_x = F2_create_matrix(n, r);
  for (j = k = 0; j < r; j++, k++)
  {
    while (d[k] != -1) k++;
    for (i = 0; i < k; i++)
      if (d[i] != -1 && F2_get_bit(x, d[i], k)) F2_set_bit(ker_x, i, j);
    F2_set_bit(ker_x, k, j);
  }
  *rank = r; avma = av; return ker_x;
}

/*********************************************************************/
/**                                                                 **/
/**                    FROM RELATIONS TO DIVISORS                   **/
/**                                                                 **/
/*********************************************************************/

/* NB: overwrites rel */
static GEN
mpqs_add_relation(GEN Y_prod, GEN N, long *ei, char *rel)
{
  pari_sp av = avma;
  GEN res;
  char *s;

  s = strchr(rel, ':') - 1;
  *s = '\0';
  
  res = remii(mulii(Y_prod, strtoi(rel)), N);

  s = strtok(s + 3, " \n");
  while (s != NULL)
  {
    long e = atol(s), i;
    if (!e) break;
    s = strtok(NULL, " \n");
    i = atol(s); /* bug in g++-3.4.1: miscompiles ei[ atol(s) ] */
    ei[i] += e;
    s = strtok(NULL, " \n");
  }
  return gerepileuptoint(av, res);
}

static char*
mpqs_get_relation(char *buf, long pos, FILE *FREL)
{
  if (fseek(FREL, pos, SEEK_SET)) pari_err(talker, "cannot seek FREL file");
  if (!fgets(buf, MPQS_STRING_LENGTH, FREL))
    pari_err(talker, "FREL file truncated?!");
  return buf;
}

#define isprobableprime(n) (miller((n),17))
static int
split(GEN N, GEN *e, GEN *res)
{
  ulong mask;
  long flag;
  GEN base;
  if (isprobableprime(N)) { *e = gen_1; return 1; }
  if (Z_issquarerem(N, &base))
  { /* squares could cost us a lot of time */
    /* GN20050707: as used now, this is always called with res!=NULL */
    *res = base;
    *e = gen_2;
    if (DEBUGLEVEL >= 5) fprintferr("MPQS: decomposed a square\n");
    return 1;
  }
  mask = 7;
  /* 5th/7th powers aren't worth the trouble. OTOH once we have the hooks for
   * dealing with cubes, higher powers can be handled essentially for free) */
  if ( (flag = is_357_power(N, &base, &mask)) )
  {
    *res = base;
    *e = utoipos(flag);
    if (DEBUGLEVEL >= 5)
      fprintferr("MPQS: decomposed a %s\n",
                 (flag == 3 ? "cube" :
                  (flag == 5 ? "5th power" : "7th power")));
    return 1;
  }
  *e = gen_0; return 0; /* known composite */
}

static GEN
mpqs_solve_linear_system(mpqs_handle_t *h, pariFILE *pFREL, long rel)
{
  FILE *FREL = pFREL->file;
  GEN N = h->N, X, Y_prod, X_plus_Y, D1, res, new_res;
  mpqs_FB_entry_t *FB = h->FB;
  pari_sp av=avma, av2, av3, lim, lim3;
  long *fpos, *ei;
  long i, j, H_cols, H_rows;
  long res_last, res_next, res_size, res_max;
  F2_matrix m, ker_m;
  long done, rank;
  char buf[MPQS_STRING_LENGTH];

  fpos = (long *) gpmalloc(rel * sizeof(long));

  m = F2_read_matrix(FREL, h->size_of_FB+1, rel, fpos);
  if (DEBUGLEVEL >= 7)
  {
    fprintferr("\\\\ MATRIX READ BY MPQS\nFREL=");
    F2_print_matrix(m, h->size_of_FB+1, rel);
    fprintferr("\n");
  }

  ker_m = mpqs_kernel(m, h->size_of_FB+1, rel, &rank);
  if (DEBUGLEVEL >= 4)
  {
    if (DEBUGLEVEL >= 7)
    {
      fprintferr("\\\\ KERNEL COMPUTED BY MPQS\nKERNEL=");
      F2_print_matrix(ker_m, rel, rank);
      fprintferr("\n");
    }
    fprintferr("MPQS: Gauss done: kernel has rank %ld, taking gcds...\n", rank);
  }

  H_rows = rel;
  H_cols = rank;

  if (!H_cols)
  { /* trivial kernel. Fail gracefully: main loop may look for more relations */
    if (DEBUGLEVEL >= 3)
      pari_warn(warner, "MPQS: no solutions found from linear system solver");
    F2_destroy_matrix(m, h->size_of_FB+1);
    F2_destroy_matrix(ker_m, rel);
    free(fpos); /* ei not yet allocated */
    avma = av; return NULL; /* no factors found */
  }

  /* If the rank is r, we can expect up to 2^r pairwise coprime factors,
   * but it may happen that a kernel basis vector contributes nothing new to
   * the decomposition.  We allocate room for up to eight factors initially
   * (certainly adequate when one or two basis vectors work), adjusting this
   * down at the end to what we actually found, or up if we are very lucky and
   * find more factors.  In the upper half of our vector, we store information
   * about which factors we know to be composite (zero) or believe to be
   * composite ((long)NULL) or suspect to be prime (one), or an exponent (two
   * or some t_INT) if it is a proper power */
  av2 = avma; lim = stack_lim(av2,1);
  if (rank > (long)BITS_IN_LONG - 2)
    res_max = VERYBIGINT; /* the common case, unfortunately */
  else
    res_max = 1L<<rank; /* max number of factors we can hope for */
  res_size = 8; /* no. of factors we can store in this res */
  res = cgetg(2*res_size+1, t_VEC);
  for (i=2*res_size; i; i--) res[i] = 0;
  res_next = res_last = 1;

  ei = (long *) gpmalloc((h->size_of_FB + 2) * sizeof(long));

  for (i = 0; i < H_cols; i++)
  { /* loop over kernel basis */
    X = Y_prod = gen_1;
    memset((void *)ei, 0, (h->size_of_FB + 2) * sizeof(long));

    av3 = avma; lim3 = stack_lim(av3,1);
    for (j = 0; j < H_rows; j++)
    {
      if (F2_get_bit(ker_m, j, i))
        Y_prod = mpqs_add_relation(Y_prod, N, ei,
                                   mpqs_get_relation(buf, fpos[j], FREL));
      if (low_stack(lim3, stack_lim(av3,1)))
      {
        if(DEBUGMEM>1) pari_warn(warnmem,"[1]: mpqs_solve_linear_system");
        Y_prod = gerepileuptoint(av3, Y_prod);
      }
    }
    Y_prod = gerepileuptoint(av3, Y_prod);

    av3 = avma; lim3 = stack_lim(av3,1);
    for (j = 2; j <= h->size_of_FB + 1; j++)
      if (ei[j])
      {
        if (ei[j] & 1) pari_err(bugparier, "MPQS (relation is a nonsquare)");
        X = remii(mulii(X,
                        Fp_powu(utoipos(FB[j].fbe_p), (ulong)ei[j]>>1, N)),
                  N);
        if (low_stack(lim3, stack_lim(av3,1)))
        {
          if(DEBUGMEM>1) pari_warn(warnmem,"[2]: mpqs_solve_linear_system");
          X = gerepileupto(av3, X);
        }
      }
    X = gerepileuptoint(av3, X);
    if (MPQS_DEBUGLEVEL >= 1)
    {
      if (signe(remii(subii(sqri(X), sqri(Y_prod)), N)))
      { /* shouldn't happen */
        fprintferr("MPQS: X^2 - Y^2 != 0 mod N\n");
        fprintferr("\tindex i = %ld\n", i);
        pari_warn(warner, "MPQS: wrong relation found after Gauss");
      }
    }
    /* At this point, X^2 == Y^2 mod N.  Indeed, something stronger is true:
     * We have gcd(X-Y, N) * gcd(X+Y, N) == N.  Why?
     * First, N divides X^2 - Y^2, so it divides the lefthand side.
     * Second, let P be any prime factor of N.  If P were to divide both
     * X-Y and X+Y, then it would divide their sum 2X.  But X (in the present
     * backwards notation!) is a product of powers of FB primes, and no FB
     * prime is a divisor of N, or we would have found out about it already
     * whilst constructing the FB.
     * Therefore in the following it is sufficient to work with gcd(X+Y, N)
     * alone, and never look at gcd(X-Y, N).
     */
    done = 0; /* (re-)counts probably-prime factors (or powers whose bases we
               * don't want to handle any further) */
    X_plus_Y = addii(X, Y_prod);
    if (res_next < 3)
    { /* we still haven't decomposed the original N, and want both a gcd and
       * its cofactor. */
      D1 = gcdii(X_plus_Y, N);
      if (is_pm1(D1) || equalii(D1,N)) { avma = av3; continue; }
      /* got something that works */
      if (DEBUGLEVEL >= 5)
        fprintferr("MPQS: splitting N after %ld kernel vector%s\n",
                   i+1, (i? "s" : ""));
      /* GN20050707 Fixed:
       * Don't divide N in place. We still need it for future X and Y_prod
       * computations! */
      gel(res,1) = diviiexact(N, D1);
      gel(res,2) = D1;
      res_last = res_next = 3;

      if ( split(gel(res,1),  &gel(res,res_size+1), &gel(res,1)) ) done++;
      if ( split(D1, &gel(res,res_size+2), &gel(res,2)) ) done++;
      if (done == 2) break;     /* both factors look prime or were powers */
      /* GN20050707: moved following line down to here, was before the
       * two split() invocations.  Very rare case anyway. */
      if (res_max == 2) break; /* two out of two possible factors seen */
      if (DEBUGLEVEL >= 5)
        fprintferr("MPQS: got two factors, looking for more...\n");
    }
    else
    { /* we already have factors */
      for (j=1; j < res_next; j++)
      { /* loop over known-composite factors */
        if (gel(res,res_size+j) && gel(res,res_size+j) != gen_0)
        {
          done++; continue; /* skip probable primes etc */
        }
        /* actually, also skip square roots of squares etc.  They are a lot
         * smaller than the original N, and should be easy to deal with later */
        av3 = avma;
        D1 = gcdii(X_plus_Y, gel(res,j));
	if (is_pm1(D1) || equalii(D1, gel(res,j))) { avma = av3; continue; }
	/* got one which splits this factor */
	if (DEBUGLEVEL >= 5)
	  fprintferr("MPQS: resplitting a factor after %ld kernel vectors\n",
		     i+1);      /* always plural */
	/* first make sure there's room for another factor */
	if (res_next > res_size)
	{ /* need to reallocate (_very_ rare case) */
	  long i1, new_size = 2*res_size;
	  GEN new_res;
	  if (new_size > res_max) new_size = res_max;
	  new_res = cgetg(2*new_size+1, t_VEC);
	  for (i1=2*new_size; i1>=res_next; i1--) new_res[i1] = 0;
	  for (i1=1; i1<res_next; i1++)
	  {
	    /* GN20050707:
	     * on-stack contents of new_res must be rejuvenated */
	    icopyifstack(res[i1], new_res[i1]); /* factors */
	    if (res[res_size+i1])
	      icopyifstack(res[res_size+i1], new_res[new_size+i1]);
	    /* primality tags */
	  }
	  res = new_res; res_size = new_size;   /* res_next unchanged */
	}
	/* now there is room; divide into existing factor and store the
	   new gcd */
	(void)dvdiiz(gel(res,j), D1, gel(res,j));
	gel(res,res_next) = D1;

	/* following overwrites the old known-composite indication at j */
	if (split( gel(res,j), &gel(res,res_size+j), &gel(res,j)) ) done++;
	/* GN20050707 Fixed:
	 * Don't increment done when the newly stored factor seems to be
	 * prime or otherwise devoid of interest - this happens later
	 * when we routinely revisit it during the present inner loop. */
	(void)split(D1, &gel(res,res_size+res_next), &gel(res,res_next));

	/* GN20050707: Following line moved down to here, was before the
	 * two split() invocations. */
	if (++res_next > res_max)
	{
	  /* all possible factors seen, outer loop postprocessing will
	   * proceed to break out of the outer loop below. */
	  break;
	}
      }       /* loop over known composite factors */

      if (res_next > res_last)
      {
	res_last = res_next - 1; /* we might have resplit more than one */
        if (DEBUGLEVEL >= 5)
          fprintferr("MPQS: got %ld factors%s\n", res_last,
                     (done < res_last ? ", looking for more..." : ""));
        res_last = res_next;
      }
      /* break out of the outer loop when we have seen res_max factors, and
       * also when all current factors are probable primes */
      if (res_next > res_max || done == res_next - 1) break;
    } /* end case of further splitting of existing factors */
    if (low_stack(lim, stack_lim(av2,1)))
    {
      long i1;
      if(DEBUGMEM>1) pari_warn(warnmem,"[3]: mpqs_solve_linear_system");
      /* gcopy would have a problem with our NULL pointers... */
      new_res = cgetg(lg(res), t_VEC);
      for (i1=2*res_size; i1>=res_next; i1--) new_res[i1] = 0;
      for (i1=1; i1<res_next; i1++)
      {
        icopyifstack(res[i1], new_res[i1]);
        /* GN20050707: the marker GENs might need rejuvenating, too */
        if (res[res_size+i1])
          icopyifstack(res[res_size+i1], new_res[res_size+i1]);
      }
      res = gerepileupto(av2, new_res);
    }
  } /* for (loop over kernel basis) */

  F2_destroy_matrix(m, h->size_of_FB+1);
  F2_destroy_matrix(ker_m, rel);
  free(ei); free(fpos);
  if (res_next < 3) { avma = av; return NULL; } /* no factors found */

  /* normal case:  convert internal format to ifac format as described in
   * src/basemath/ifactor1.c  (triples of components: value, exponent, class
   * [unknown, known composite, known prime]),  clean up and return result */
  res_last = res_next - 1; /* number of distinct factors found */
  new_res = cgetg(3*res_last + 1, t_VEC);
  if (DEBUGLEVEL >= 6)
    fprintferr("MPQS: wrapping up vector of %ld factors\n", res_last);
  for (i=1,j=1; i <= res_last; i++)
  {
    GEN F = gel(res, res_size+i);
    icopyifstack(res[i], new_res[j++]); /* factor */
    gel(new_res,j++) = /* exponent */
      F ? (F == gen_0 ? gen_1
                      : (isonstack(F) ? icopy(F) : F))
        : gen_1; /* F was NULL */
    gel(new_res,j++) = /* class */
      F == gen_0 ? gen_0 :      /* known composite */
        NULL;           /* base of power or suspected prime --
                                   mark as `unknown' */
    if (DEBUGLEVEL >= 6)
      fprintferr("\tpackaging %ld: %Z ^%ld (%s)\n", i, res[i],
                 itos(gel(new_res,j-2)), (F == gen_0 ? "comp." : "unknown"));
  }
  return gerepileupto(av, new_res);
}

/*********************************************************************/
/**                                                                 **/
/**               MAIN ENTRY POINT AND DRIVER ROUTINE               **/
/**                                                                 **/
/*********************************************************************/


/* All percentages below are actually fixed-point quantities scaled by 10
 * (value of 1 means 0.1%, 1000 means 100%) */

/* Factors N using the self-initializing multipolynomial quadratic sieve
 * (SIMPQS).  Returns one of the two factors, or (usually) a vector of factors
 * and exponents and information about which ones are still composite, or NULL
 * when something goes wrong or when we can't seem to make any headway. */

/* TODO: this function to be renamed mpqs_main() with several extra parameters,
 * with mpqs() as a wrapper for the standard case, so we can do partial runs
 * across several machines etc.  (from gp or a dedicated C program). --GN */
static GEN
mpqs_i(mpqs_handle_t *handle)
{
  GEN N = handle->N, fact; /* will in the end hold our factor(s) */

  /* XXX fix type! */
  mpqs_int32_t size_of_FB; /* size of the factor base */
  mpqs_FB_entry_t *FB; /* factor base */

  mpqs_int32_t M;               /* sieve interval size [-M, M] */


  /* local loop / auxiliary vars */
  ulong p;

  /* XX already exists in the handle, keep for convenience */
  long lp_bound;                /* size limit for large primes */
  long lp_scale;                /* ...relative to largest FB prime */

  /* bookkeeping */
  long tc;                      /* # of candidates found in one iteration */
  long tp;                      /* # of recently sorted LP rels */
  long tff = 0;                 /* # recently found full rels from sieving */
  long tfc;                     /* # full rels recently combined from LPs */
  double tfc_ratio = 0;         /* recent (tfc + tff) / tff */
  ulong sort_interval;          /* determine when to sort and merge */
  ulong followup_sort_interval; /* temporary files (scaled percentages) */
  long percentage = 0;          /* scaled by 10, see comment above */
  double net_yield;
  long total_full_relations = 0, total_partial_relations = 0, total_no_cand = 0;
  long vain_iterations = 0, good_iterations = 0, iterations = 0;
#ifdef MPQS_USE_HISTOGRAMS
  long histo_checkpoint = MPQS_MIN_CANDS_FOR_HISTO;
#endif

  pariFILE *pFNEW, *pLPNEW, *pCOMB, *pFREL, *pLPREL;
  char *dir, *COMB_str, *FREL_str, *FNEW_str, *LPREL_str, *LPNEW_str, *TMP_str;

/* END: global variables to disappear as soon as possible */

/******************************/


  pari_sp av = avma;

  if (DEBUGLEVEL >= 4)
  {
    (void)timer2();
    fprintferr("MPQS: number to factor N = %Z\n", N);
  }

  handle->digit_size_N = decimal_len(N);
  if (handle->digit_size_N > MPQS_MAX_DIGIT_SIZE_KN)
  {
    pari_warn(warner, "MPQS: number too big to be factored with MPQS,\n\tgiving up");
    return NULL;
  }

  if (DEBUGLEVEL >= 4)
    fprintferr("MPQS: factoring number of %ld decimal digits\n",
               handle->digit_size_N);

  mpqs_find_k(handle);
  if (DEBUGLEVEL >= 5) fprintferr("MPQS: found multiplier %ld for N\n",
                                  handle->_k.k);
  handle->kN = mulis(N, handle->_k.k);

  if (!mpqs_set_parameters(handle))
  {
    pari_warn(warner,
        "MPQS: number too big to be factored with MPQS,\n\tgiving up");
    return NULL;
  }

  size_of_FB = handle->size_of_FB;
  M = handle->M;
  sort_interval = handle->first_sort_point;
  followup_sort_interval = handle->sort_pt_interval;

  if (DEBUGLEVEL >= 5)
    fprintferr("MPQS: creating factor base and allocating arrays...\n");
  FB = mpqs_create_FB(handle, &p);
  if (p)
  {
    /* free(FB); */
    if (DEBUGLEVEL >= 4)
      fprintferr("\nMPQS: found factor = %ld whilst creating factor base\n", p);
    avma = av; return utoipos(p);
  }
  mpqs_sieve_array_ctor(handle);
  mpqs_poly_ctor(handle);

  /* XXX commentary above said that we were taking the max here, but we have
   * always been taking the min - huh?! Thus cutting off LPs at 10^7 even
   * when the largest FB prime was bigger than 1.5E6. */
  lp_bound = handle->largest_FB_p;
  if (lp_bound > MPQS_LP_BOUND) lp_bound = MPQS_LP_BOUND;
  /* don't allow large primes to have room for two factors both bigger than
   * what the FB contains (...yet!) */
  lp_scale = handle->lp_scale;
  if (lp_scale >= handle->largest_FB_p)
    lp_scale = handle->largest_FB_p - 1;
  lp_bound *= lp_scale;
  handle->lp_bound = lp_bound;

  handle->dkN = gtodouble(handle->kN);
  /* compute the threshold and fill in the byte-sized scaled logarithms */
  mpqs_set_sieve_threshold(handle);

  if (!mpqs_locate_A_range(handle)) return NULL;

  if (DEBUGLEVEL >= 4)
  {
    fprintferr("MPQS: sieving interval = [%ld, %ld]\n", -(long)M, (long)M);
    /* that was a little white lie, we stop one position short at the top */
    fprintferr("MPQS: size of factor base = %ld\n",
               (long)size_of_FB);
    fprintferr("MPQS: striving for %ld relations\n",
               (long)handle->target_no_rels);
    fprintferr("MPQS: coefficients A will be built from %ld primes each\n",
               (long)handle->omega_A);
    fprintferr("MPQS: primes for A to be chosen near FB[%ld] = %ld\n",
               (long)handle->index2_FB,
               (long)FB[handle->index2_FB].fbe_p);
    fprintferr("MPQS: smallest prime used for sieving FB[%ld] = %ld\n",
               (long)handle->index1_FB,
               (long)FB[handle->index1_FB].fbe_p);
    fprintferr("MPQS: largest prime in FB = %ld\n",
               (long)handle->largest_FB_p);
    fprintferr("MPQS: bound for `large primes' = %ld\n", (long)lp_bound);
  }

  if (DEBUGLEVEL >= 5)
  {
    fprintferr("MPQS: sieve threshold = %u\n",
               (unsigned int)handle->sieve_threshold);
  }

  if (DEBUGLEVEL >= 4)
  {
    fprintferr("MPQS: first sorting at %ld%%, then every %3.1f%% / %3.1f%%\n",
               sort_interval/10, followup_sort_interval/10.,
               followup_sort_interval/20.);
  }


  /* main loop which
   * - computes polynomials and their zeros (SI)
   * - does the sieving
   * - tests candidates of the sieve array */

  /* Let (A, B_i) the current pair of coeffs. If i == 0 a new A is generated */
  handle->index_j = (mpqs_uint32_t)-1;  /* increment below will have it start at 0 */

  if (DEBUGLEVEL >= 5) fprintferr("MPQS: starting main loop\n");
  /* compute names for the temp files we'll need */
  dir = pari_unique_dir("MPQS");
  TMP_str   = mpqs_get_filename(dir, "LPTMP");
  FREL_str  = mpqs_get_filename(dir, "FREL");
  FNEW_str  = mpqs_get_filename(dir, "FNEW");
  LPREL_str = mpqs_get_filename(dir, "LPREL");
  LPNEW_str = mpqs_get_filename(dir, "LPNEW");
  COMB_str  = mpqs_get_filename(dir, "COMB");
#define unlink_all()\
      pari_unlink(FREL_str);\
      pari_unlink(FNEW_str);\
      pari_unlink(LPREL_str);\
      pari_unlink(LPNEW_str);\
      if (pCOMB) pari_unlink(COMB_str);\
      rmdir(dir); free(dir);

  pFREL = pari_fopen(FREL_str,  WRITE); pari_fclose(pFREL);
  pLPREL = pari_fopen(LPREL_str,  WRITE); pari_fclose(pLPREL);
  pFNEW = pari_fopen(FNEW_str,  WRITE);
  pLPNEW= pari_fopen(LPNEW_str, WRITE);
  pCOMB = NULL;

  for(;;)
  { /* FNEW and LPNEW are open for writing */
    /* XXX the incrementing and stepping logic here really belongs into
     * XXX mpqs_self_init(), too */
    iterations++;
    /* when all of the B's have already been used, choose new A ;
       this is indicated by setting index_j to 0 */
    if (handle->index_j == (mpqs_uint32_t)(handle->no_B - 1))
    {
      handle->index_j = 0;
      handle->index_i++;        /* count finished A's */
    }
    else
      handle->index_j++;

    /* self initialization: compute polynomial and its zeros */
    mpqs_self_init(handle);
    if (handle->bin_index == 0)
    { /* have run out of primes for A */
      /* We might change some parameters.  For the moment, simply give up */
      if (DEBUGLEVEL >= 2)
        fprintferr("MPQS: Ran out of primes for A, giving up.\n");
      pari_fclose(pFNEW);
      pari_fclose(pLPNEW);
      /* FREL, LPREL are closed at this point */
      unlink_all(); avma = av; return NULL;
    }

    memset((void*)(handle->sieve_array), 0, (M << 1) * sizeof(unsigned char));
    mpqs_sieve(handle);

    tc = mpqs_eval_sieve(handle);
    total_no_cand += tc;
    if (DEBUGLEVEL >= 6)
      fprintferr("MPQS: found %lu candidate%s\n", tc, (tc==1? "" : "s"));

    if (tc)
    {
      long t = mpqs_eval_cand(handle, tc, pFNEW->file, pLPNEW->file);
      total_full_relations += t;
      tff += t;
      good_iterations++;
    }

#ifdef MPQS_USE_HISTOGRAMS
    if (handle->do_histograms && !handle->done_histograms &&
        total_no_cand >= histo_checkpoint)
    {
      int res = mpqs_eval_histograms(handle);
      if (res >= 0)
      {
        /* retry later */
        if (res > 0)
          /* histo_checkpoint *= 2.6; */
          handle->do_histograms = 0; /* no, don't retry later */
        else
          histo_checkpoint += (MPQS_MIN_CANDS_FOR_HISTO /* >> 1 */);
      }
      else
        handle->done_histograms = 1;
    }
#endif

    percentage =
      (long)((1000.0 * total_full_relations) / handle->target_no_rels);

    if ((ulong)percentage < sort_interval) continue;
    /* most main loops continue here! */

    /* Extra processing when we have completed a sort interval: */
    if (DEBUGLEVEL >= 3)
    {
      if (DEBUGLEVEL >= 4)
        fprintferr("\nMPQS: passing the %3.1f%% sort point, time = %ld ms\n",
                   sort_interval/10., timer2());
      else
        fprintferr("\nMPQS: passing the %3.1f%% sort point\n",
                   sort_interval/10.);
      flusherr();
    }

#ifdef MPQS_USE_HISTOGRAMS
    /* XXX following to go away when we're done debugging... */
    if (DEBUGLEVEL >= 6)
      mpqs_print_histo(handle);
#endif

    /* sort LPNEW and merge it into LPREL, diverting combinables into COMB */
    pari_fclose(pLPNEW);
    (void)mpqs_sort_lp_file(LPNEW_str);
    pCOMB = pari_fopen(COMB_str, WRITE);
    tp = mpqs_mergesort_lp_file(LPREL_str, LPNEW_str, TMP_str, pCOMB);
    pari_fclose(pCOMB);
    pLPNEW = pari_fopen(LPNEW_str, WRITE);

    /* combine whatever there is to be combined */
    tfc = 0;
    if (tp > 0)
    {
      /* build full relations out of large prime relations */
      pCOMB = pari_fopen(COMB_str, READ);
      tfc = mpqs_combine_large_primes(handle, pCOMB->file, pFNEW, &fact);
      pari_fclose(pCOMB);
      /* now FREL, LPREL are closed and FNEW, LPNEW are still open */
      if (fact)
      { /* factor found during combining */
        if (DEBUGLEVEL >= 4)
        {
          fprintferr("\nMPQS: split N whilst combining, time = %ld ms\n",
                     timer2());
          fprintferr("MPQS: found factor = %Z\n", fact);
        }
        pari_fclose(pLPNEW);
        pari_fclose(pFNEW);
        unlink_all();
        return gerepileupto(av, fact);
      }
      total_partial_relations += tp;
    }

    /* sort FNEW and merge it into FREL */
    pari_fclose(pFNEW);
    (void)mpqs_sort_lp_file(FNEW_str);
    /* definitive count (combinables combined, and duplicates removed) */
    total_full_relations = mpqs_mergesort_lp_file(FREL_str, FNEW_str, TMP_str, NULL);
    /* FNEW stays closed until we need to reopen it for another iteration */

    /* Due to the removal of duplicates, percentage may actually decrease at
     * this point.  Looks funny in the diagnostics but is nothing to worry
     * about: we _are_ making progress. */
    percentage =
      (long)((1000.0 * total_full_relations) / handle->target_no_rels);
    net_yield =
      (total_full_relations * 100.) / (total_no_cand ? total_no_cand : 1);
    vain_iterations =
      (long)((1000.0 * (iterations - good_iterations)) / iterations);

    /* Now estimate the current full relations yield rate:  we directly see
     * each time through the main loop how many full relations we're getting
     * as such from the sieve  (tff since the previous checkpoint),  but
     * only at checkpoints do we see how many we're typically combining
     * (tfc).  So we're really producing (tfc+tff)/tff as many full rels,
     * and when we get close to 100%, we should bias the next interval by
     * the inverse ratio.
     * Avoid drawing conclusions from too-small samples during very short
     * follow-on intervals  (in this case we'll just re-use an earlier
     * estimated ratio). */
    if ((tfc >= 16) && (tff >= 20))
      tfc_ratio = (tfc + tff + 0.) / tff; /* floating-point division */
    tff = 0;                    /* reset this count (tfc is always fresh) */

    if (percentage >= 1000)     /* when Gauss had failed */
      sort_interval = percentage + 2;
    else if (percentage >= 820)
    {
      if (tfc_ratio > 1.)
      {
        if (percentage + (followup_sort_interval >> 1) * tfc_ratio > 994)
        {
          /* aim for a _slight_ overshoot */
          sort_interval = (ulong)(percentage + 2 +
            (1000 - percentage) / tfc_ratio);
        }
        else if (percentage >= 980)
          sort_interval = percentage + 8;
        else
          sort_interval = percentage + (followup_sort_interval >> 1);
      }
      else
      {
        if (percentage >= 980)
          sort_interval = percentage + 10;
        else
          sort_interval = percentage + (followup_sort_interval >> 1);
        if (sort_interval >= 1000 && percentage < 1000)
          sort_interval = 1000;
      }
    }
    else
      sort_interval = percentage + followup_sort_interval;

    if (DEBUGLEVEL >= 4)
    {
      fprintferr("MPQS: done sorting%s, time = %ld ms\n",
                 tp > 0 ? " and combining" : "", timer2());
      fprintferr("MPQS: found %3.1f%% of the required relations\n",
                 percentage/10.);
      if (DEBUGLEVEL >= 5)
      { /* total_full_relations are always plural */
        /* GN20050708: present code doesn't succeed in discarding all
         * dups, so don't lie about it... */
        fprintferr("MPQS: found %ld full relations\n",
                   total_full_relations);
	if (lp_scale > 1)
        fprintferr("MPQS:   (%ld of these from partial relations)\n",
                   total_partial_relations);
        fprintferr("MPQS: Net yield: %4.3g full relations per 100 candidates\n",
                   net_yield);
        fprintferr("MPQS:            %4.3g full relations per 100 polynomials\n",
                   (total_full_relations * 100.) / iterations);
        /* XX also say something about avg cands per poly, huh?
         * XX And for that matter, about the number of polynomials used
         * XX so far. */
        fprintferr("MPQS: %4.1f%% of the polynomials yielded no candidates\n",
                   vain_iterations/10.);
        fprintferr("MPQS: next sort point at %3.1f%%\n", sort_interval/10.);
      }
    }

    if (percentage < 1000)
    {
      pFNEW = pari_fopen(FNEW_str, WRITE);
      /* LPNEW and FNEW are again open for writing */
      continue; /* main loop */
    }

    /* percentage >= 1000, which implies total_full_relations > size_of_FB:
       try finishing it off */

    /* solve the system over F_2 */
    /* present code does NOT in fact guarantee absence of dup FRELs,
     * therefore removing the adjective "distinct" for the time being */
    if (DEBUGLEVEL >= 4)
      fprintferr("\nMPQS: starting Gauss over F_2 on %ld relations\n",
                 total_full_relations);
    pFREL = pari_fopen(FREL_str, READ);
    fact = mpqs_solve_linear_system(handle, pFREL, total_full_relations);
    pari_fclose(pFREL);

    if (fact)
    { /* solution found */
      if (DEBUGLEVEL >= 4)
      {
        fprintferr("\nMPQS: time in Gauss and gcds = %ld ms\n", timer2());
        if (typ(fact) == t_INT) fprintferr("MPQS: found factor = %Z\n", fact);
        else
        {
          long j, nf = (lg(fact)-1)/3;
          if (nf == 2)
            /* GN20050707: Changed the arrangement of the two factors,
             * to match the debug diagnostics in mpqs_solve_linear_system()
             * above */
            fprintferr("MPQS: found factors = %Z\n\tand %Z\n",
                        fact[1], fact[4]);
          else
          {
            /* GN20050707: Changed loop to scan upwards instead of downwards,
             * to match the debug diagnostics in mpqs_solve_linear_system()
             * above */
            fprintferr("MPQS: found %ld factors =\n", nf);
            for (j=1; j<=nf; j++)
              fprintferr("\t%Z%s\n", fact[3*j-2], (j<nf ? "," : ""));
          }
        }
      }
      pari_fclose(pLPNEW);
      unlink_all();
      /* fact not safe for a gerepilecopy(): segfaults on one of the NULL
       * markers. However, it is a nice connected object, and it resides
       * already the top of the stack, so... --GN */
      return gerepileupto(av, fact);
    }
    else
    {
      if (DEBUGLEVEL >= 4)
      {
        fprintferr("\nMPQS: time in Gauss and gcds = %ld ms\n", timer2());
        fprintferr("MPQS: no factors found.\n");
        if (percentage <= MPQS_ADMIT_DEFEAT)
          fprintferr("\nMPQS: restarting sieving ...\n");
        else
          fprintferr("\nMPQS: giving up.\n");
      }
      if (percentage > MPQS_ADMIT_DEFEAT)
      {
        pari_fclose(pLPNEW);
        unlink_all(); avma = av; return NULL;
      }
      pFNEW = pari_fopen(FNEW_str, WRITE);
    }
  } /* main loop */
}
GEN
mpqs(GEN N)
{
  mpqs_handle_t *handle = mpqs_handle_ctor(N);
  GEN fact = mpqs_i(handle);
  mpqs_handle_dtor(handle); return fact;
}
