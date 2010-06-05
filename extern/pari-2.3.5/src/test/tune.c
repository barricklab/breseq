/* $Id: tune.c 7522 2005-12-09 18:14:24Z kb $

Copyright (C) 2001  The PARI group.

This file is part of the PARI/GP package.

PARI/GP is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation. It is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY WHATSOEVER.

Check the License for details. You should have received a copy of it, along
with the package; see the file 'COPYING'. If not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA. */

/* This file is a quick hack adapted from gmp-4.0 tuning utilities
 * (T. Granlund et al.)
 *
 * (GMU MP Library is Copyright Free Software Foundation, Inc.) */
#define PARI_TUNE
#include <pari.h>
#include <paripriv.h>

#define numberof(x) (sizeof(x) / sizeof((x)[0]))

int option_trace = 0;
double Step_Factor = .01; /* small steps by default */
ulong DFLT_mod = 17UL;

typedef struct {
  ulong reps, type;
  long *var, size;
  GEN x, y;
} speed_param;

typedef double (*speed_function_t)(speed_param *s);

typedef struct {
  int               kernel;
  const char        *name;
  long              *var;
  int               type; /* t_INT or t_REAL */
  long              min_size;
  long              max_size;
  speed_function_t  fun1;
  speed_function_t  fun2;
  double            step_factor; /* how much to step sizes (rounded down) */
  double            stop_factor;
} tune_param;

/* ========================================================== */
/* To use GMP cycle counting functions, look for GMP in Oxxx/Makefile */
#ifdef GMP_TIMER
/* needed to link with gmp-4.0/tune/{time,freq}.o */
int speed_option_verbose = 0;
extern double speed_unittime;
extern int    speed_precision;
void speed_starttime(void);
double speed_endtime(void);
#else
static pari_timer __T;
static double speed_unittime = 1e-4;
static int    speed_precision= 1000;
static void speed_starttime() { TIMERstart(&__T); }
static double speed_endtime() { return (double)TIMER(&__T)/1000.; }
#endif

/* ========================================================== */
/* int, n words */
static GEN
rand_INT(long n)
{
  pari_sp av = avma;
  GEN x, N = int2n(n*BITS_IN_LONG);
  do x = genrand(N); while (lgefint(x) != n+2);
  return gerepileuptoint(av, x);
}
/* real, n words */
static GEN
rand_REAL(long n) { return gmul2n(itor(rand_INT(n), n+2),-BITS_IN_LONG*n); }

/* Flx, degree n */
static GEN
rand_Flx(long n)
{
  pari_sp av = avma;
  GEN x;
  do x = FpX_rand(n+1, 0, utoipos(DFLT_mod)); while (degpol(x) < n);
  return gerepileuptoleaf(av, ZX_to_Flx(x, DFLT_mod));
}

/* normalized Flx, degree n */
static GEN
rand_NFlx(long n)
{
  pari_sp av = avma;
  GEN x = gadd(monomial(gen_1,n,0), FpX_rand(n, 0, utoipos(DFLT_mod)));
  return gerepileuptoleaf(av, ZX_to_Flx(x, DFLT_mod));
}

#define t_Flx  100
#define t_NFlx 101

static GEN
rand_g(long n, long type)
{
  switch (type) {
    case t_INT:  return rand_INT(n);
    case t_REAL: return rand_REAL(n);
    case t_Flx:  return rand_Flx(n);
    case t_NFlx: return rand_NFlx(n);
  }
  return NULL;
}

/* ========================================================== */
#define TIME_FUN(call) {\
  {                                      \
    pari_sp av = avma;                   \
    int i;                               \
    speed_starttime();                   \
    i = (s)->reps;                       \
    do { call; avma = av; } while (--i); \
  }                                      \
  return speed_endtime();                \
}

#define  enable(s) (*(s->var)=lg(s->x)-2)/* enable  asymptotically fastest */
#define disable(s) (*(s->var)=lg(s->x)+1)/* disable asymptotically fastest */

static double speed_mulrr(speed_param *s)
{ disable(s); TIME_FUN(mulrr(s->x, s->y)); }
static double speed_karamulrr(speed_param *s)
{ enable(s);  TIME_FUN(mulrr(s->x, s->y)); }

static double speed_mulii(speed_param *s)
{ disable(s); TIME_FUN(mulii(s->x, s->y)); }
static double speed_karamulii(speed_param *s)
{ enable(s); TIME_FUN(mulii(s->x, s->y)); }

static double speed_exp(speed_param *s)
{ disable(s); TIME_FUN(mpexp(s->x)); }
static double speed_expnewton(speed_param *s)
{ enable(s);  mplog2(lg(s->x)+2); TIME_FUN(mpexp(s->x)); }

static double speed_log(speed_param *s)
{ disable(s); TIME_FUN(mplog(s->x)); }
static double speed_logagm(speed_param *s)
{ enable(s);  TIME_FUN(mplog(s->x)); }

static double speed_logcx(speed_param *s)
{ GEN z; setexpo(s->x,0); z = gadd(gen_1, gmul(gi, s->x));
  glog(z,s->size);
  disable(s); TIME_FUN(glog(z,s->size)); }
static double speed_logcxagm(speed_param *s)
{ GEN z; setexpo(s->x,0); z = gadd(gen_1, gmul(gi, s->x));
  glog(z,s->size);
  enable(s); TIME_FUN(glog(z,s->size)); }

static double speed_atan(speed_param *s)
{ setexpo(s->x, 0); disable(s);
  gatan(s->x, 0);
  TIME_FUN(gatan(s->x, 0)); }
static double speed_atanagm(speed_param *s)
{ setexpo(s->x, 0); enable(s);
  gatan(s->x, 0);
  TIME_FUN(gatan(s->x, 0)); }

static double speed_sqri (speed_param *s)
{ disable(s); TIME_FUN(sqri(s->x)); }
static double speed_karasqri (speed_param *s)
{ enable(s);  TIME_FUN(sqri(s->x)); }

static double speed_divrr(speed_param *s)
{ disable(s); TIME_FUN(divrr(s->x, s->y)); }
static double speed_divrrgmp(speed_param *s)
{ enable(s); TIME_FUN(divrr(s->x, s->y)); }

static double speed_invmod(speed_param *s)
{ GEN T; disable(s); TIME_FUN(invmod(s->x, s->y, &T)); }
static double speed_invmodgmp(speed_param *s)
{ GEN T; enable(s); TIME_FUN(invmod(s->x, s->y, &T)); }

static double speed_Flx_sqr(speed_param *s)
{ ulong p = DFLT_mod; disable(s); TIME_FUN(Flx_sqr(s->x, p)); }
static double speed_Flx_karasqr(speed_param *s)
{ ulong p = DFLT_mod; enable(s); TIME_FUN(Flx_sqr(s->x, p)); }

static double speed_Flx_inv(speed_param *s)
{ ulong p = DFLT_mod; disable(s);
  TIME_FUN(Flx_invmontgomery(s->x, p)); }
static double speed_Flx_invnewton(speed_param *s)
{ ulong p = DFLT_mod; enable(s);
  TIME_FUN(Flx_invmontgomery(s->x, p)); }

static double speed_Flx_mul(speed_param *s)
{ ulong p = DFLT_mod; disable(s); TIME_FUN(Flx_mul(s->x, s->y, p)); }
static double speed_Flx_karamul(speed_param *s)
{ ulong p = DFLT_mod; enable(s); TIME_FUN(Flx_mul(s->x, s->y, p)); }

#define INIT_RED(s, op)                                 \
  long i, lx = lg(s->x);                                \
  op = cgeti(2*lx - 2);                                 \
  op[1] = evallgefint(2*lx - 2) | evalsigne(1);         \
  for (i=2; i<lx; i++) op[i]      = s->x[i];            \
  for (i=2; i<lx; i++) op[lx-2+i] = s->x[i];            \
  modBIL(s->y) |= 1; /* make sure modulus is odd */
static double speed_redc(speed_param *s) {
  ulong inv = (ulong)-invrev(modBIL(s->y));
  GEN op; INIT_RED(s, op);
  TIME_FUN( red_montgomery(op, s->y, inv) ); };
static double speed_modii(speed_param *s) {
  GEN op; INIT_RED(s, op);
  TIME_FUN( remii(op, s->y) ); };
static double speed_remiimul(speed_param *s) {
  GEN sM = init_remiimul(s->y);
  GEN op; INIT_RED(s, op);
  TIME_FUN( remiimul(op, sM) ); }

static double speed_Flxq_pow_redc(speed_param *s) {
  ulong p = DFLT_mod;
  enable(s); TIME_FUN( Flxq_pow(polx_Flx(0), utoipos(p), s->y, p) );
}
static double speed_Flxq_pow_mod(speed_param *s) {
  ulong p = DFLT_mod;
  disable(s); TIME_FUN( Flxq_pow(polx_Flx(0), utoipos(p), s->y, p) );
}

enum { PARI = 1, GMP = 2 };
#ifdef PARI_KERNEL_GMP
#  define AVOID PARI
#else
#  define AVOID GMP
#endif

/* Thresholds are set in this order. If f() depends on g(), g() should
 * occur first */
#define var(a) # a, &a
static tune_param param[] = {
{PARI,var(KARATSUBA_MULI_LIMIT),   t_INT, 4,0, speed_mulii,speed_karamulii},
{PARI,var(KARATSUBA_SQRI_LIMIT),   t_INT, 4,0, speed_sqri,speed_karasqri},
{0,   var(KARATSUBA_MULR_LIMIT),   t_REAL,4,0, speed_mulrr,speed_karamulrr},
{PARI,var(MONTGOMERY_LIMIT),       t_INT, 3,0, speed_redc,speed_modii},
{0,   var(REMIIMUL_LIMIT),         t_INT, 3,0, speed_modii,speed_remiimul},
{GMP, var(DIVRR_GMP_LIMIT),        t_REAL,4,0, speed_divrr,speed_divrrgmp},
{0,   var(EXPNEWTON_LIMIT),        t_REAL,64,0, speed_exp,speed_expnewton},
{0,   var(LOGAGM_LIMIT),           t_REAL,4,0, speed_log,speed_logagm},
{0,   var(LOGAGMCX_LIMIT),         t_REAL,3,0, speed_logcx,speed_logcxagm,0.05},
{0,   var(AGM_ATAN_LIMIT),         t_REAL,20,0, speed_atan,speed_atanagm,0.05},
{GMP, var(INVMOD_GMP_LIMIT),       t_INT, 3,0, speed_invmod,speed_invmodgmp},
{0,   var(Flx_MUL_LIMIT),          t_Flx, 4,0, speed_Flx_mul,speed_Flx_karamul},
{0,   var(Flx_SQR_LIMIT),          t_Flx, 4,0, speed_Flx_sqr,speed_Flx_karasqr},
{0,   var(Flx_INVMONTGOMERY_LIMIT),t_NFlx,10,30000,
                                   speed_Flx_inv,speed_Flx_invnewton,0.3},
{0,  var(Flx_POW_MONTGOMERY_LIMIT),t_NFlx,1,0,
                                   speed_Flxq_pow_redc,speed_Flxq_pow_mod}
};

/* ========================================================== */
int ndat = 0, allocdat = 0;
struct dat_t {
  long size;
  double d;
} *dat = NULL;

int
double_cmp_ptr(double *x, double *y) { return (int)(*x - *y); }

double
time_fun(speed_function_t fun, speed_param *s)
{
  const double TOLERANCE = 1.005; /* 0.5% */
  pari_sp av = avma;
  double t[30];
  ulong i, j, e;

  s->x = rand_g(s->size, s->type);
  s->y = rand_g(s->size, s->type); s->reps = 1;
  for (i = 0; i < numberof(t); i++)
  {
    for (;;)
    {
      double reps_d;
      t[i] = fun(s);
      if (!t[i]) { s->reps *= 10; continue; }
      if (t[i] >= speed_unittime * speed_precision) break;

      /* go to a value of reps to make t[i] >= precision */
      reps_d = ceil (1.1 * s->reps
                     * speed_unittime * speed_precision
                     / max(t[i], speed_unittime));
      if (reps_d > 2e9 || reps_d < 1.0)
        pari_err(talker, "Fatal error: new reps bad: %.2f\n", reps_d);

      s->reps = (ulong)reps_d;
    }
    t[i] /= s->reps;

    /* require 3 values within TOLERANCE when >= 2 secs, 4 when below */
    e = (t[0] >= 2.0)? 3: 4;

   /* Look for e many t[]'s within TOLERANCE of each other to consider a
      valid measurement.  Return smallest among them.  */
    if (i >= e)
    {
      qsort (t, i+1, sizeof(t[0]), (QSCOMP)double_cmp_ptr);
      for (j = e-1; j < i; j++)
        if (t[j] <= t[j-e+1] * TOLERANCE) { avma = av; return t[j-e+1]; }
    }
  }
  pari_err(talker,"couldn't measure time");
  return -1.0; /* not reached */
}

void
add_dat(long size, double d)
{
  if (ndat == allocdat)
  {
    allocdat += max(allocdat, 100);
    dat = (struct dat_t*) gprealloc((void*)dat, allocdat * sizeof(dat[0]));
  }
  dat[ndat].size = size;
  dat[ndat].d    = d; ndat++;
}

long
analyze_dat(int final)
{
  double  x, min_x;
  int     j, min_j;

  /* If the threshold is set at dat[0].size, any positive values are bad. */
  x = 0.0;
  for (j = 0; j < ndat; j++)
    if (dat[j].d > 0.0) x += dat[j].d;

  if (final && option_trace >= 3)
  {
    printf("\n");
    printf("x is the sum of the badness from setting thresh at given size\n");
    printf("  (minimum x is sought)\n");
    printf("size=%ld  first x=%.4f\n", dat[j].size, x);
  }

  min_x = x;
  min_j = 0;

  /* When stepping to the next dat[j].size, positive values are no longer
     bad (so subtracted), negative values become bad (so add the absolute
     value, meaning subtract). */
  for (j = 0; j < ndat; j++)
  {
    if (final && option_trace >= 3)
      printf ("size=%ld  x=%.4f\n", dat[j].size, x);

    if (x < min_x) { min_x = x; min_j = j; }
    x -= dat[j].d;
  }
  return min_j;
}

void
print_define(const char *name, long value)
{ printf("#define __%-25s  %5ld\n\n", name, value); }

void
Test(tune_param *param)
{
  int since_positive, since_change, thresh, new_thresh;
  speed_param s;

  if (param->kernel == AVOID) { print_define(param->name, -1); return; }

#define DEFAULT(x,n)  if (! (param->x))  param->x = (n);
  DEFAULT(fun2, param->fun1);
  DEFAULT(step_factor, Step_Factor);
  DEFAULT(stop_factor, 1.2);
  DEFAULT(max_size, 10000);

  s.type = param->type;
  s.size = param->min_size;
  s.var  = param->var;
  ndat = since_positive = since_change = thresh = 0;
  if (option_trace >= 1)
    printf("Setting %s... (default %ld)\n", param->name, *(param->var));
  if (option_trace >= 2)
  {
    printf("              algorithm-A  algorithm-B   ratio  possible\n");
    printf("               (seconds)    (seconds)    diff    thresh\n");
  }

  for(;;)
  {
    double t1, t2, d;
    t1 = time_fun(param->fun1, &s);
    t2 = time_fun(param->fun2, &s);
    if (t2 >= t1) d = (t2 - t1) / t2;
    else          d = (t2 - t1) / t1;

    add_dat(s.size, d);
    new_thresh = analyze_dat(0);

    if (option_trace >= 2)
      printf ("size =%4ld     %.8f   %.8f  % .4f %c  %ld\n",
               s.size, t1,t2, d, d < 0? '#': ' ', dat[new_thresh].size);

#define SINCE_POSITIVE 20
#define SINCE_CHANGE 50
    /* Stop if method B has been consistently faster for a while */
    if (d >= 0)
      since_positive = 0;
    else
      if (++since_positive > SINCE_POSITIVE)
      {
        if (option_trace >= 1)
          printf ("Stop: since_positive (%d)\n", SINCE_POSITIVE);
        break;
      }
    /* Stop if method A has become slower by a certain factor */
    if (t1 >= t2 * param->stop_factor)
    {
      if (option_trace >= 1)
        printf ("Stop: t1 >= t2 * factor (%.1f)\n", param->stop_factor);
      break;
    }
    /* Stop if threshold implied hasn't changed for a while */
    if (thresh != new_thresh)
      since_change = 0, thresh = new_thresh;
    else
      if (++since_change > SINCE_CHANGE)
      {
        if (option_trace >= 1)
          printf ("Stop: since_change (%d)\n", SINCE_CHANGE);
        break;
      }
    s.size += max((long)floor(s.size * param->step_factor), 1);
    if (s.size >= param->max_size)
    {
      if (option_trace >= 1)
        printf ("Stop: max_size (%ld). Disable Algorithm B?\n",param->max_size);
      break;
    }
  }
  thresh = dat[analyze_dat(1)].size;
  print_define(param->name, thresh);
  *(param->var) = thresh; /* set to optimal value for next tests */
}

void error(char **argv) {
  long i;
  printf("usage: tune [-t] [-s step_factor] [-p mod] [-u unittime] var1 var2 ...\n");
  printf("Tunable variables: (omitting variable indices tunes everybody)\n");
  for (i = 0; i < (long)numberof(param); i++)
    printf("  %2ld: %-25s (default %4ld)\n", i, param[i].name, *(param[i].var));
  exit(1);
}

int
main(int argc, char **argv)
{
  int i, r, n = 0;
  GEN v;
  pari_init(4000000, 2);
  v = new_chunk(argc);
  for (i = 1; i < argc; i++)
  {
    char *s = argv[i];
    if (*s == '-') {
      switch(*++s) {
        case 't': option_trace += 2; break;
        case 'p':
          if (!*++s)
          {
            if (++i == argc) error(argv);
            s = argv[i];
          }
          DFLT_mod = itou(gp_read_str(s)); break;
        case 's':
          if (!*++s)
          {
            if (++i == argc) error(argv);
            s = argv[i];
          }
          Step_Factor = atof(s); break;
        case 'u': s++;
          if (!*++s)
          {
            if (++i == argc) error(argv);
            s = argv[i];
          }
          speed_unittime = atof(s); break;
        default: error(argv);
      }
    } else {
      if (!isdigit((int)*s)) error(argv);
      r = atol(s); if (r >= (long)numberof(param) || r < 0) error(argv);
      v[n++] = r;
    }
  }
  if (n) { for (i = 0; i < n; i++) Test(&param[ v[i] ]); return 0; }
  n = numberof(param);
  for (i = 0; i < n; i++) Test(&param[i]);
  return 0;
}
