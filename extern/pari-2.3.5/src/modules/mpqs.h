/* - debug support */

#ifdef MPQS_DEBUG_VERYVERBOSE
#  ifndef MPQS_DEBUG_VERBOSE
#  define MPQS_DEBUG_VERBOSE
#  endif
#endif

#ifdef MPQS_DEBUG_VERBOSE
#  ifndef MPQS_DEBUG
#  define MPQS_DEBUG
#  endif
#  define PRINT_IF_VERBOSE(x) fprintferr(x)
#else
#  define PRINT_IF_VERBOSE(x)
#endif

#ifdef MPQS_DEBUG
#  define MPQS_DEBUGLEVEL 1000  /* infinity */
#else
#  define MPQS_DEBUGLEVEL DEBUGLEVEL
#endif

/* - string and external file stuff for the relations "database" */

#ifndef SEEK_SET
#  define SEEK_SET 0
#endif

#ifdef __CYGWIN32__
/* otherwise fseek() goes crazy due to silent \n <--> LF translations */
#  define WRITE "wb"
#  define READ "rb"
#else
#  define WRITE "w"
#  define READ "r"
#endif

#define MPQS_STRING_LENGTH         (4 * 1024UL)

/* - non-configurable sizing parameters */

#define MPQS_POSSIBLE_MULTIPLIERS  5 /* how many values for k we'll try */
/* following must be in range of the cand_multipliers table below */
#define MPQS_MULTIPLIER_SEARCH_DEPTH 5 /* how many primes to inspect per k */

/* XXX this comment is now totally off base...
 * `large primes' must be smaller than
 *   max(MPQS_LP_BOUND, largest_FB_p) * MPQS_LP_FACTOR
 * - increased this with the idea of capping it at about 2^30
 * (XX but see comment in mpqs() where we actually use this: we take the
 * min, not the max)
 */
#define MPQS_LP_BOUND              12500000 /* works for 32 and 64bit */

/* see mpqs_locate_A_range() for an explanation of the following.  I had
 * some good results with about -log2(0.85) but in the range I was testing,
 * this shifts the primes for A only by one position in the FB.  Don't go
 * over the top with this one... */
#define MPQS_A_FUDGE               0.15 /* ~ -log2(0.9) */

#define MPQS_CANDIDATE_ARRAY_SIZE  2000 /* max. this many cand's per poly */

#ifdef MPQS_USE_HISTOGRAMS
/* histogram evaluation/feedback available when size_of_FB exceeds this: */
#  define MPQS_MIN_SIZE_FB_FOR_HISTO 600
/* min number of candidates to look at before evaluating histograms */
#  define MPQS_MIN_CANDS_FOR_HISTO   4000
/* min number of full relations to have been created before etc. */
#  define MPQS_MIN_FRELS_FOR_HISTO   110

/* see mpqs_eval_histograms() for explanation of the following */
#  define MPQS_HISTO_FREL_QUANTILE   2.4
#  define MPQS_HISTO_DROP_LIMIT      3.6
#  define MPQS_HISTO_LPREL_BASEFLOW  1.4
#endif

/* give up when nothing found after ~1.5 times the required number of
 * relations has been computed  (N might be a prime power or the
 * parameters might be exceptionally unfortunate for it) */
#define MPQS_ADMIT_DEFEAT 1500


/* - structures, types, and constants */

/* -- reasonably-sized integers */
#ifdef LONG_IS_64BIT
typedef int  mpqs_int32_t;
typedef unsigned int  mpqs_uint32_t;
typedef unsigned long mpqs_uint64_t;
#else
typedef long mpqs_int32_t;
typedef unsigned long mpqs_uint32_t;
typedef struct {
  ulong _w0;
  ulong _w1;
} mpqs_uint64_t;
#endif

/* -- we'll sometimes want to use the machine's native size here, and
 * sometimes  (for future double-large-primes)  force 64 bits - thus: */
typedef union mpqs_invp {
  ulong _ul;
  mpqs_uint64_t _u64;
} mpqs_invp_t;

/* -- factor base entries should occupy 32 bytes  (and we'll keep them
 * aligned, for good L1 cache hit rates).  Some of the entries will be
 * abused for e.g. -1 and (factors of) k instead for real factor base
 * primes, and for a sentinel at the end.  This is why __p is a signed
 * field.-- The two start fields depend on the current polynomial and
 * keep changing during sieving, the flags will also change depending on
 * the current A. */
/* Let (z1, z2) be the roots of Q(x) = A x^2 + Bx + C mod p_i; then
 * Q(z1 + p_i Z) == 0 mod p_i and Q(z2 + p_i Z) == 0 mod p_i;
 * start_1, start_2 are the positions where p_i divides Q(x) for the
 * first time, already adjusted for the fact that the sieving array,
 * nominally [-M, M], is represented by a 0-based C array of length
 * 2M + 1.  For the prime factors of A and those of k, the two roots
 * are equal mod p_i. */

#define MPQS_FB_ENTRY_PAD 32

typedef union mpqs_FB_entry {
  char __pad[MPQS_FB_ENTRY_PAD];
  struct {
    mpqs_int32_t __p;           /* the prime p */
    /* XX following two are not yet used: */
    float __flogp;              /* its logarithm as a 4-byte float */
    mpqs_invp_t __invp;         /* 1/p mod 2^64 or 2^BITS_IN_LONG */
    mpqs_int32_t __start1;      /* representatives of the two APs mod p */
    mpqs_int32_t __start2;
    mpqs_uint32_t __sqrt_kN;    /* sqrt(kN) mod p */
    unsigned char __val;        /* 8-bit approx. scaled log for sieving */
    unsigned char __flags;
  } __entry;
} mpqs_FB_entry_t;

/* --- convenience accessor macros for the preceding: */
#define fbe_p           __entry.__p
#define fbe_flogp       __entry.__flogp
#define fbe_invp        __entry.__invp
#define fbe_start1      __entry.__start1
#define fbe_start2      __entry.__start2
#define fbe_sqrt_kN     __entry.__sqrt_kN
#define fbe_logval      __entry.__val
#define fbe_flags       __entry.__flags

/* --- flag bits for fbe_flags: */
/* XXX TODO! */

#define MPQS_FBE_CLEAR       0x0 /* no flags */

/* following used for odd FB primes, and applies to the divisors of A but not
 * those of k.  Must occupy the rightmost bit because we also use it as a
 * shift count after extracting it from the byte. */
#define MPQS_FBE_DIVIDES_A   0x1ul /* and Q(x) mod p only has degree 1 */

/* XX tentative: one bit to mark normal FB primes,
 * XX one to mark the factors of k,
 * XX one to mark primes used in sieving,
 * XX later maybe one to mark primes of which we'll be tracking the square,
 * XX one to mark primes currently in use for A;
 * XX once we segment the FB, one bit marking the members of the first segment
 */

/* -- multiplier k and associated quantities: More than two prime factors
* for k will be pointless in practice, thus capping them at two. */
#define MPQS_MAX_OMEGA_K 2
typedef struct mpqs_multiplier {
  mpqs_uint32_t k;              /* the multiplier (odd, squarefree) */
  mpqs_uint32_t omega_k;        /* number (>=0) of primes dividing k */
  mpqs_uint32_t kp[MPQS_MAX_OMEGA_K]; /* prime factors of k, if any */
} mpqs_multiplier_t;

static mpqs_multiplier_t cand_multipliers[] = {
  {  1, 0, {  0,  0} },
  {  3, 1, {  3,  0} },
  {  5, 1, {  5,  0} },
  {  7, 1, {  7,  0} },
  { 11, 1, { 11,  0} },
  { 13, 1, { 13,  0} },
  { 15, 2, {  3,  5} },
  { 17, 1, { 17,  0} },
  { 19, 1, { 19,  0} },
  { 21, 2, {  3,  7} },
  { 23, 1, { 23,  0} },
  { 29, 1, { 29,  0} },
  { 31, 1, { 31,  0} },
  { 33, 2, {  3, 11} },
  { 35, 2, {  5,  7} },
  { 37, 1, { 37,  0} },
  { 39, 2, {  3, 13} },
  { 41, 1, { 41,  0} },
  { 43, 1, { 43,  0} },
  { 47, 1, { 47,  0} },
  { 51, 2, {  3, 17} },
  { 53, 1, { 53,  0} },
  { 55, 2, {  5, 11} },
  { 57, 2, {  3, 19} },
  { 59, 1, { 59,  0} },
  { 61, 1, { 61,  0} },
  { 65, 2, {  5, 13} },
  { 67, 1, { 67,  0} },
  { 69, 2, {  3, 23} },
  { 71, 1, { 71,  0} },
  { 73, 1, { 73,  0} },
  { 77, 2, {  7, 11} },
  { 79, 1, { 79,  0} },
  { 83, 1, { 83,  0} },
  { 85, 2, {  5, 17} },
  { 87, 2, {  3, 29} },
  { 89, 1, { 89,  0} },
  { 91, 2, {  7, 13} },
  { 93, 2, {  3, 31} },
  { 95, 2, {  5, 19} },
  { 97, 1, { 97,  0} }
};

/* -- the array of (Chinese remainder) idempotents which add/subtract up to
 * the middle coefficient B, and for convenience, the FB subscripts of the
 * primes in current use for A.  We keep these together since both arrays
 * are of the same size and are used at the same times. */
typedef struct mqps_per_A_prime {
  GEN _H;                       /* summand for B */
  mpqs_int32_t _i;              /* subscript into FB */
} mpqs_per_A_prime_t;

/* following cooperate with names of local variables in the self_init fcns.
 * per_A_pr must exist and be an alias for the eponymous handle pointer for
 * all of these, and FB must exist and correspond to the handle FB pointer
 * for all but the first two of them. */
#define MPQS_H(i) (per_A_pr[i]._H)
#define MPQS_I(i) (per_A_pr[i]._i)
#define MPQS_AP(i) (FB[per_A_pr[i]._i].fbe_p)
#define MPQS_LP(i) (FB[per_A_pr[i]._i].fbe_flogp)
#define MPQS_SQRT(i) (FB[per_A_pr[i]._i].fbe_sqrt_kN)
#define MPQS_FLG(i) (FB[per_A_pr[i]._i].fbe_flags)


/* -- the array of addends / subtrahends for changing polynomials during
 * self-initialization: (1/A) H[i] mod p_j, with i subscripting the inner
 * array in each entry, and j choosing the entry in an outer array.
 * Entries will occupy 64 bytes each no matter what  (which imposes one
 * sizing restriction: at most 17 prime factors for A;  thus i will range
 * from 0 to at most 15.)  This wastes a little memory for smaller N but
 * makes it easier for compilers to generate efficient code. */
/* XX At present, memory locality vis-a-vis accesses to this array is good
 * XX in the slow (new A) branch of mpqs_self_init(), but poor in the fast
 * XX (same A, new B) branch, which now loops over the outer array index,
 * XX reading just one field of each inner array each time through the FB
 * XX loop.  This doesn't really harm, but will improve one day when we do
 * XX segmented sieve arrays with the associated segmented FB-range accesses.
 */
#define MPQS_MAX_OMEGA_A 17
typedef struct mpqs_inv_A_H {
  mpqs_uint32_t _i[MPQS_MAX_OMEGA_A - 1];
} mpqs_inv_A_H_t;

#define MPQS_INV_A_H(i,j) (inv_A_H[j]._i[i])

/* -- global handle for keeping track of everything used throughout any one
 * factorization attempt.  The order of the fields is roughly determined by
 * wanting to keep the most frequently used stuff near the beginning. */

typedef struct mpqs_handle {
  /* pointers into gpmalloc()d memory which must be freed at the end: */
  unsigned char *sieve_array;   /* 0-based, representing [-M,M-1] */
  unsigned char *sieve_array_end; /* points at sieve_array[M-1] */
  mpqs_FB_entry_t *FB;          /* (aligned) FB array itself */
  long *candidates;             /* collects promising sieve subscripts */
  char *relations;              /* freshly found relations (strings) */
  long *relaprimes;             /* prime/exponent pairs in a relation */
  mpqs_inv_A_H_t *inv_A_H;      /* self-init: (aligned) stepping array, and */
  mpqs_per_A_prime_t *per_A_pr; /* FB subscripts of primes in A etc. */

  /* other stuff that's being used all the time */
  mpqs_int32_t M;               /* sieving over |x| <= M */
  mpqs_int32_t size_of_FB;      /* # primes in FB (or dividing k) */
  /* the following three are in non-descending order, and the first two must
   * be adjusted for omega_k at the beginning */
  mpqs_int32_t index0_FB;       /* lowest subscript into FB of a "real" prime
                                 * (i.e. other than -1, 2, factors of k) */
  mpqs_int32_t index1_FB;       /* lowest subscript into FB for sieving */
  mpqs_int32_t index2_FB;       /* primes for A are chosen relative to this */
  unsigned char index2_moved;   /* true when we're starved for small A's */
  unsigned char sieve_threshold; /* distinguishes candidates in sieve */
#ifdef MPQS_USE_HISTOGRAMS
  /* histogram feedback */
  unsigned char do_histograms;  /* (boolean) enable histogram updating */
  unsigned char done_histograms; /* histos have been eval'd for feedback */
  /* more gpmalloc()d memory here: */
  long *histo_full;             /* distribution of full rels from sieve */
  long *histo_lprl;             /* - of LP rels from sieve */
  long *histo_drop;             /* - of useless candidates */
#endif
  GEN N;                        /* given number to be factored */
  GEN kN;                       /* N with multiplier (on PARI stack) */
  /* quantities associated with the current polynomial; all these also
   * live in preallocated slots on the PARI stack: */
  GEN A;                        /* leading coefficient */
  GEN B;                        /* middle coefficient */
#ifdef MPQS_DEBUG
  GEN C;                        /* and constant coefficient */
#endif
  mpqs_int32_t omega_A;         /* number of primes going into each A */
  mpqs_int32_t no_B;            /* number of B's for each A: 2^(omega_A-1) */
  double l2_target_A;           /* ~log2 of desired typical A */
  /* counters and bit pattern determining and numbering the current
   * polynomial: */
  mpqs_uint32_t bin_index;      /* bit pattern for selecting primes for A */
  mpqs_uint32_t index_i;        /* running count of A's */
  mpqs_uint32_t index_j;        /* B's ordinal number in A's cohort */
  /* XXX one more to follow here... */

  /* further sizing parameters: */
  mpqs_int32_t target_no_rels;  /* target number of full relations */
  mpqs_int32_t largest_FB_p;    /* largest prime in the FB */
  mpqs_int32_t pmin_index1;	/* lower bound for primes used for sieving */
  mpqs_int32_t lp_scale;        /* factor by which LPs may exceed FB primes */

  mpqs_int32_t first_sort_point; /* when to sort and combine */
  mpqs_int32_t sort_pt_interval; /* (in units of 1/1000) */

  /* XXX subscripts determining where to pick primes for A... */

  /* XX lp_bound might have to be mpqs_int64_t ? or mpqs_invp_t ? */
  long lp_bound;                /* cutoff for Large Primes */
  long digit_size_N;
  long digit_size_kN;
  mpqs_multiplier_t _k;         /* multiplier k and associated quantities */
  double tolerance;             /* controls the tightness of the sieve */
  double dkN;                   /* - double prec. approximation of kN */
  double l2sqrtkN;              /* ~log2(sqrt(kN)) */
  double l2M;                   /* ~log2(M) (cf. below) */
  /* XX need an index2_FB here to remember where to start picking primes */

  /* XX put statistics here or keep them as local variables in mpqs() ? */

  /* bookkeeping pointers to containers of aligned memory chunks: */
  void *FB_chunk;               /* (unaligned) chunk containing the FB */
  void *invAH_chunk;            /* (unaligned) chunk for self-init array */

} mpqs_handle_t;

/* -- for Gauss/Lanczos */
typedef ulong *F2_row;
typedef F2_row *F2_matrix;

/* -- sizing table entries */

/* The "tolerance" is explained below apropos of mpqs_set_sieve_threshold().
 * The LP scale, for very large kN, prevents us from accumulating vast amounts
 * of LP relations with little chance of hitting any particular large prime
 * a second time and being able to combine a full relation from two LP ones;
 * however, the sieve threshold (determined by the tolerance) already works
 * against very large LPs being produced.-- The present relations "database"
 * can detect duplicate full relations only during the sort/combine phases,
 * so we must do some sort points even for tiny kN where we do not admit
 * large primes at all.
 * Some constraints imposed by the present implementation:
 * + omega_A should be at least 3, and no more than MPQS_MAX_OMEGA_A
 * + The size of the FB must be large enough compared to omega_A
 *   (about 2*omega_A + 3, but this is always true below) */
/* XXX Changes needed for segmented mode:
 * XXX When using it (kN large enough),
 * XXX - M must become a multiple of the (cache block) segment size
 * XXX   (or to keep things simple: a multiple of 32K)
 * XXX - we need index3_FB to seperate (smaller) primes used for normal
 * XXX   sieving from larger ones used with transaction buffers
 * XXX   (and the locate_A_range and associated logic must be changed to
 * XXX   cap index2_FB below index3_FB instead of below size_of_FB)
 */
typedef struct mpqs_parameterset {
  float tolerance;              /* "mesh width" of the sieve */
  /* XX following subject to further change */
  mpqs_int32_t lp_scale;        /* factor by which LPs may exceed FB primes */
  mpqs_int32_t M;               /* size of half the sieving interval */
  mpqs_int32_t size_of_FB;      /* #primes to use for FB (including 2) */
  mpqs_int32_t omega_A;         /* #primes to go into each A */
  /* following is auto-adjusted to account for prime factors of k inserted
   * near the start of the FB.  NB never ever sieve on the prime 2  (which
   * would just contribute a constant at each sieve point). */
  mpqs_int32_t pmin_index1;     /* lower bound for primes used for sieving */
  /* the remaining two are expressed in percent  (of the target number of full
   * relations),  and determine when we stop sieving to review the contents
   * of the relations DB and sort them and combine full relations from LP
   * ones.  Note that the handle has these in parts per thousand instead. */
  mpqs_int32_t first_sort_point;
  mpqs_int32_t sort_pt_interval;
} mpqs_parameterset_t;

/* - the table of sizing parameters itself */

/* indexed by size of kN in decimal digits, subscript 0 corresponding to
 * 9 (or fewer) digits */
static mpqs_parameterset_t mpqs_parameters[] =
{ /*       tol lp_scl     M   szFB  oA pmx1 1st  sti */
  {  /*9*/ 0.8,   1,    900,    20,  3,   5, 70,  8},
  { /*10*/ 0.8,   1,    900,    21,  3,   5, 70,  8},
  { /*11*/ 0.8,   1,    920,    22,  3,   5, 70,  6},
  { /*12*/ 0.8,   1,    960,    24,  3,   5, 70,  6},
  { /*13*/ 0.8,   1,   1020,    26,  3,   5, 70,  6},
  { /*14*/ 0.8,   1,   1100,    29,  3,   5, 70,  6},
  { /*15*/ 0.8,   1,   1200,    32,  3,   5, 60,  8},
  { /*16*/ 0.8,   1,   1500,    35,  3,   5, 60,  8},
  { /*17*/ 0.8,   1,   1900,    40,  3,   5, 60,  8},
  { /*18*/ 0.8,   1,   2500,    60,  3,   5, 50, 10},
  { /*19*/ 0.8,   1,   3200,    80,  3,   5, 50, 10},
  { /*20*/ 0.8,   1,   4000,   100,  3,   5, 40, 10},
  { /*21*/ 0.8,   1,   4300,   100,  3,   5, 40, 10},
  { /*22*/ 0.8,   1,   4500,   120,  3,   5, 40, 10},
  { /*23*/ 0.8,   1,   4800,   140,  3,   5, 30, 10},
  { /*24*/ 0.8,   1,   5100,   160,  4,   7, 30, 10},
  { /*25*/ 0.8,   1,   5400,   180,  4,   7, 30, 10},
  { /*26*/ 0.9,   1,   5700,   200,  4,   7, 30, 10},
  { /*27*/ 1.12,  1,   6000,   220,  4,   7, 30, 10},
  { /*28*/ 1.17,  1,   6300,   240,  4,  11, 30, 10},
  { /*29*/ 1.22,  1,   6500,   260,  4,  11, 30, 10},
  { /*30*/ 1.30,  1,   6800,   325,  4,  11, 20, 10},
  { /*31*/ 1.33,  1,   7000,   355,  4,  13, 20, 10},
  { /*32*/ 1.36,  1,   7200,   375,  5,  13, 20, 10},
  { /*33*/ 1.40,  1,   7400,   400,  5,  13, 20, 10},
  { /*34*/ 1.43,  1,   7600,   425,  5,  17, 20, 10},
  /* around here, sieving takes long enough to make it worthwhile recording
   * LP relations into their separate output files, although they tend not
   * yet to contribute a lot to the full relations until we get up to around
   * 47 digits or so. */
  { /*35*/ 1.48, 30,   7800,   550,  5,  17, 20, 10},
  { /*36*/ 1.53, 45,   8100,   650,  5,  17, 20, 10},
  { /*37*/ 1.60, 60,   9000,   750,  6,  19, 20, 10},
  { /*38*/ 1.66, 70,  10000,   850,  6,  19, 20, 10},
  { /*39*/ 1.69, 80,  11000,   950,  6,  23, 20, 10},
  /* around here, the largest prime in FB becomes comparable to M in size */
  { /*40*/ 1.69, 80,  12500,  1000,  6,  23, 20, 10},
  { /*41*/ 1.69, 80,  14000,  1150,  6,  23, 10, 10},
  { /*42*/ 1.69, 80,  15000,  1300,  6,  29, 10, 10},
  { /*43*/ 1.69, 80,  16000,  1500,  6,  29, 10, 10},
  { /*44*/ 1.69, 80,  17000,  1700,  7,  31, 10, 10},
  { /*45*/ 1.69, 80,  18000,  1900,  7,  31, 10, 10},
  { /*46*/ 1.69, 80,  20000,  2100,  7,  37, 10, 10},
  { /*47*/ 1.69, 80,  25000,  2300,  7,  37, 10, 10},
  { /*48*/ 1.69, 80,  27500,  2500,  7,  37, 10, 10},
  { /*49*/ 1.72, 80,  30000,  2700,  7,  41, 10, 10},
  { /*50*/ 1.75, 80,  35000,  2900,  7,  41, 10, 10},
  { /*51*/ 1.80, 80,  40000,  3000,  7,  43, 10, 10},
  { /*52*/ 1.85, 80,  50000,  3200,  7,  43, 10, 10},
  { /*53*/ 1.90, 80,  60000,  3500,  7,  47, 10, 10},
  { /*54*/ 1.95, 80,  70000,  3800,  7,  47, 10, 10},
  { /*55*/ 1.95, 80,  80000,  4100,  7,  53, 10, 10},
  { /*56*/ 1.95, 80,  90000,  4400,  7,  53, 10,  8},
  { /*57*/ 2.00, 80, 100000,  4700,  8,  53, 10,  8},
  { /*58*/ 2.05, 80, 110000,  5000,  8,  59, 10,  8},
  { /*59*/ 2.10, 80, 120000,  5400,  8,  59, 10,  8},
  { /*60*/ 2.15, 80, 130000,  5800,  8,  61, 10,  8},
  { /*61*/ 2.20, 80, 140000,  6100,  8,  61, 10,  8},
  { /*62*/ 2.25, 80, 150000,  6400,  8,  67, 10,  6},
  { /*63*/ 2.39, 80, 160000,  6700,  8,  67, 10,  6},
  { /*64*/ 2.30, 80, 165000,  7000,  8,  67, 10,  6},
  { /*65*/ 2.31, 80, 170000,  7300,  8,  71, 10,  6},
  { /*66*/ 2.32, 80, 175000,  7600,  8,  71, 10,  6},
  { /*67*/ 2.33, 80, 180000,  7900,  8,  73, 10,  6},
  { /*68*/ 2.34, 80, 185000,  8200,  8,  73, 10,  6},
  { /*69*/ 2.35, 80, 190000,  8600,  8,  79,  8,  6},
  { /*70*/ 2.36, 80, 195000,  8800,  8,  79,  8,  6},
  { /*71*/ 2.37, 80, 200000,  9000,  9,  79,  8,  6},
  { /*72*/ 2.38, 80, 205000,  9250,  9,  83,  5,  5},
  { /*73*/ 2.41, 80, 210000,  9500,  9,  83,  5,  5},
  { /*74*/ 2.46, 80, 220000,  9750,  9,  83,  5,  5},
  { /*75*/ 2.51, 80, 230000, 10000,  9,  89,  5,  5},
  { /*76*/ 2.56, 80, 240000, 10500,  9,  89,  5,  5},
  { /*77*/ 2.58, 80, 250000, 11200,  9,  89,  5,  5},
  { /*78*/ 2.60, 80, 260000, 12500,  9,  89,  5,  5},
  { /*79*/ 2.63, 80, 270000, 14000,  9,  97,  5,  4},
  { /*80*/ 2.65, 80, 280000, 15500,  9,  97,  5,  4},
  { /*81*/ 2.72, 80, 300000, 17000,  9,  97,  4,  4},
  { /*82*/ 2.77, 80, 320000, 18500,  9, 101,  4,  4},
  { /*83*/ 2.82, 80, 340000, 20000, 10, 101,  4,  4},
  { /*84*/ 2.84, 80, 360000, 21500, 10, 103,  4,  4},
  { /*85*/ 2.86, 80, 400000, 23000, 10, 103,  4,  3},
  { /*86*/ 2.88, 80, 460000, 24500, 10, 107,  4,  3},
  /* architectures with 1MBy L2 cache will become noticeably slower here
   * as 2*M exceeds that mark - to be addressed in a future version by
   * segmenting the sieve interval */
  { /*87*/ 2.90, 80, 520000, 26000, 10, 107,  4,  3},
  { /*88*/ 2.91, 80, 580000, 27500, 10, 109,  4,  3},
  { /*89*/ 2.92, 80, 640000, 29000, 10, 109,  4,  3},
  { /*90*/ 2.93, 80, 700000, 30500, 10, 113,  2,  2},
  { /*91*/ 2.94, 80, 770000, 32200, 10, 113,  2,  2},
  /* entries below due to Thomas Denny, never tested */
  { /*92*/ 3.6, 90, 2000000, 35000,  9, 113,  2,  2},
  { /*93*/ 3.7, 90, 2000000, 37000,  9, 113,  2,  2},
  { /*94*/ 3.7, 90, 2000000, 39500,  9, 127,  2,  2},
  { /*95*/ 3.7, 90, 2500000, 41500,  9, 127,  2,  2},
  { /*96*/ 3.8, 90, 2500000, 45000, 10, 127,  2,  2},
  { /*97*/ 3.8, 90, 2500000, 47500, 10, 131,  2,  2},
  { /*98*/ 3.7, 90, 3000000, 51000, 10, 131,  2,  2},
  { /*99*/ 3.8, 90, 3000000, 53000, 10, 133,  2,  2},
  {/*100*/ 3.8, 90, 3500000, 51000, 10, 133,  2,  2},
  {/*101*/ 3.8, 90, 3500000, 54000, 10, 139,  2,  2},
  {/*102*/ 3.8, 90, 3500000, 57000, 10, 139,  2,  2},
  {/*103*/ 3.9, 90, 4000000, 61000, 10, 139,  2,  2},
  {/*104*/ 3.9, 90, 4000000, 66000, 10, 149,  2,  2},
  {/*105*/ 3.9, 90, 4000000, 70000, 10, 149,  2,  2},
  {/*106*/ 3.9, 90, 4000000, 75000, 10, 151,  2,  2},
  {/*107*/ 3.9, 90, 4000000, 80000, 10, 151,  2,  2},
};

#define MPQS_MAX_DIGIT_SIZE_KN 107
