/* $Id: parinf.h 7615 2006-01-21 13:38:03Z kb $

Copyright (C) 2000  The PARI group.

This file is part of the PARI/GP package.

PARI/GP is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation. It is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY WHATSOEVER.

Check the License for details. You should have received a copy of it, along
with the package; see the file 'COPYING'. If not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA. */

typedef struct {
  GEN x; /* defining polynomial (monic, integral) */
  GEN dK; /* disc(K) */
  GEN index; /* [O_K : Z[X]/(x)] */
  GEN bas;  /* Z-basis of O_K (t_VEC of t_POL) */
  long r1; /* number of real places of K */
/* possibly NULL = irrelevant or not computed */
  GEN lead; /* leading coeff of initial polynomial defining K if non monic */
  GEN dx;   /* disc(x) */
  GEN basden; /* [nums(bas), dens(bas)] */
} nfbasic_t;

GEN nfbasic_to_nf(nfbasic_t *T, GEN ro, long prec);

typedef struct {
  GEN x;
  GEN ro;   /* roots of x */
  long r1;
  GEN basden;
  long prec;
/* possibly -1 = irrelevant or not computed */
  long extraprec;
/* possibly NULL = irrelevant or not computed */
  GEN M;
  GEN G;
} nffp_t;

void remake_GM(GEN nf, nffp_t *F, long prec);

#define id_PRINCIPAL 0
#define id_PRIME     1
#define id_MAT       2

/* for initalg_i */
#define nf_ROUND2      64
#define nf_NOROOTS     32
#define nf_PARTIALFACT 16 /* and allbase */
#define nf_RED          8
#define nf_PARTRED      2
#define nf_ORIG         1

/* for rnfpolredabs */
#define nf_ABSOLUTE     2
#define nf_ADDZK      256

/* for isprincipal */
#define nf_GEN   1
#define nf_FORCE 2
#define nf_GIVEPREC 4
#define nf_GENMAT 8
#define nf_GEN_IF_PRINCIPAL 512

/* for buchray */
#define nf_INIT  4

/* for buchall */
#define nf_ROOT1 512
#define nf_UNITS 1024
enum { fupb_NONE, fupb_RELAT, fupb_LARGE, fupb_PRECI, fupb_BACH };

/* for discray */
#define nf_REL  1
#define nf_COND 2

/* for polredabs */
#define nf_ALL   4
#define nf_RAW   8

/* for lllgramall[gen] */
#define lll_ALL 0
#define lll_KER 1
#define lll_IM  2
#define lll_GRAM 0x100

/* for minim */
#define min_ALL   0
#define min_FIRST 1
#define min_PERF  2
#define min_VECSMALL 3
#define min_VECSMALL2 4

/* for fincke_pohst() */
typedef struct FP_chk_fun {
  GEN (*f)(void *,GEN);
  /* f_init allowed to permute the columns of u and r */
  GEN (*f_init)(struct FP_chk_fun*,GEN,GEN);
  GEN (*f_post)(struct FP_chk_fun*,GEN,GEN);
  void *data;
  long skipfirst;
} FP_chk_fun;

GEN initalg_i(GEN x, long flag, long prec);
GEN fincke_pohst(GEN a,GEN BOUND,long stockmax,long PREC, FP_chk_fun *CHECK);
GEN polredfirstpol(GEN x, long flag, FP_chk_fun *CHECK);

/* for ideallog / zlog */
typedef struct {
  GEN lists; /* lists[i] = */
  GEN ind;  /* ind[i] = start of vector */
  GEN P, e; /* finit part of conductor = prod P^e */
  GEN archp; /* archimedean part of conductor, in permutation form */
  long n;  /* total number of generators for all (O_K/P^e)^* and (O_K/f_oo) */
  GEN U; /* base change matrix from generators to bid.gen */
} zlog_S;

void init_zlog_bid(zlog_S *S, GEN bid);
GEN  log_gen_arch(zlog_S *S, long index);
GEN  log_gen_pr(zlog_S *S, long index, GEN nf, long e);
GEN  zlog(GEN nf, GEN a, GEN sgn, zlog_S *S);

/* conversions basis / alg */

/* nf a genuine NF, x an nfelt (t_COL) or t_MAT whose columns represent nfelts.
 * Return the corresponding elements as t_POLs (implicitly mod nf.pol) */
#define coltoliftalg(nf,x) (gmul(gel((nf),7), (x)))
GEN    algtobasis_i(GEN nf, GEN x);
GEN    algtobasis_cp(GEN nf, GEN x);
GEN    basistoalg_i(GEN nf, GEN x);
GEN    poltobasis(GEN nf,GEN x);
GEN    coltoalg(GEN nf,GEN x);

/* Other number fields routines */
GEN    arch_mul(GEN x, GEN y);
GEN    archstar_full_rk(GEN x, GEN bas, GEN v, GEN gen);
GEN    bnrGetSurj(GEN bnr1, GEN bnr2);
GEN    check_and_build_cycgen(GEN bnf);
double check_bach(double cbach, double B);
GEN    checkbnf_i(GEN bnf);
GEN    checknf_i(GEN nf);
void   check_ZKmodule(GEN x, char *s);
void   dbg_rel(long s, GEN col);
GEN    element_mulidid(GEN nf, long i, long j);
GEN    element_powid_mod_p(GEN nf, long I, GEN n, GEN p);
GEN    eltabstorel(GEN x, GEN T, GEN pol, GEN k);
GEN    eltmulid_get_table(GEN nf, long i);
GEN    eltreltoabs(GEN rnfeq, GEN x);
GEN    galoisbig(GEN x, long prec);
GEN    get_arch(GEN nf,GEN x,long prec);
GEN    get_arch_real(GEN nf,GEN x,GEN *emb,long prec);
GEN    get_bas_den(GEN bas);
GEN    get_hnfid(GEN nf, GEN x);
GEN    get_mul_table(GEN x,GEN bas,GEN invbas);
GEN    get_nfindex(GEN bas);
GEN    get_proj_modT(GEN basis, GEN T, GEN p);
GEN    get_roots(GEN x,long r1,long prec);
GEN    get_theta_abstorel(GEN T, GEN pol, GEN k);
GEN    idealaddtoone_i(GEN nf, GEN x, GEN y);
GEN    idealcoprime_fact(GEN nf, GEN x, GEN fy);
GEN    idealhermite_aux(GEN nf, GEN x);
GEN    idealsqrtn(GEN nf, GEN x, GEN gn, int strict);
GEN    init_unif_mod_fZ(GEN L);
GEN    init_units(GEN BNF);
long   int_elt_val(GEN nf, GEN x, GEN p, GEN bp, GEN *t);
GEN    make_integral(GEN nf, GEN L0, GEN f, GEN *listpr);
GEN    maxord_i(GEN p, GEN f, long mf, GEN w, long flag);
GEN    modprV(GEN z, GEN nf,GEN modpr);
GEN    nfpol_to_Flx(GEN nf, GEN pol, ulong *ptp);
GEN    nfreducemodideal_i(GEN x0,GEN ideal);
GEN    nfrootsall_and_pr(GEN nf, GEN pol);
GEN    norm_by_embed(long r1, GEN x);
GEN    perm_to_arch(GEN nf, GEN archp);
GEN    pidealprimeinv(GEN nf, GEN x);
GEN    primedec_apply_kummer(GEN nf,GEN pol,long e,GEN p);
GEN    prodid(GEN nf, GEN I);
GEN    pr_norm(GEN pr);
GEN    quadhilbertreal(GEN D, long prec);
GEN    rnfallbase(GEN nf, GEN pol, GEN *pD, GEN *pd, GEN *pfi);
GEN    rnfequation_i(GEN A, GEN B, long *pk, GEN *pLPRS);
GEN    special_anti_uniformizer(GEN nf, GEN pr);
GEN    sqr_by_tab(GEN tab, GEN x);
GEN    subgroupcondlist(GEN cyc, GEN bound, GEN listKer);
GEN    T2_from_embed_norm(GEN x, long r1);
void   testprimes(GEN bnf, ulong bound);
GEN    to_Fp_simple(GEN nf, GEN x, GEN ffproj);
GEN    unif_mod_fZ(GEN pr, GEN F);
GEN    unnf_minus_x(GEN x);
void   wr_rel(GEN col);
GEN    zideallog_sgn(GEN nf, GEN x, GEN sgn, GEN bid);
GEN    zlog_units(GEN nf, GEN U, GEN sgnU, GEN bid);
GEN    zlog_units_noarch(GEN nf, GEN U, GEN bid);
GEN    zsign_from_logarch(GEN Larch, GEN invpi, GEN archp);

/* Dedekind zeta */
GEN  zeta_get_limx(long r1, long r2, long bit);
long zeta_get_i0(long r1, long r2, long bit, GEN limx);
long zeta_get_N0(GEN C,  GEN limx);
