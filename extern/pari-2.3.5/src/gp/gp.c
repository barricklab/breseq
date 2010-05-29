/* $Id: gp.c 8781 2007-06-13 20:50:50Z kb $

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
/**                                                               **/
/**                        PARI CALCULATOR                        **/
/**                                                               **/
/*******************************************************************/
#include "pari.h"
#include "paripriv.h"
#include "../language/anal.h"
#include "gp.h"
#include "../graph/rect.h"

#ifdef _WIN32
#  include <windows.h>
#  ifndef WINCE
#    include <process.h>
#  endif
#endif

#ifdef READLINE
BEGINEXTERN
#  ifdef READLINE_LIBRARY
#    include <readline.h>
#  else
#    include <readline/readline.h>
#  endif
ENDEXTERN
#endif

#define skip_space(s) while (isspace((int)*s)) s++
#define skip_alpha(s) while (isalpha((int)*s)) s++

/*******************************************************************/
/**                                                               **/
/**                    TEXMACS-SPECIFIC STUFF                     **/
/**                                                               **/
/*******************************************************************/
static int tm_is_waiting = 0, tm_did_complete = 0;

/* tell TeXmacs GP will start outputing data */
static void
tm_start_output(void)
{
  if (!tm_is_waiting) { printf("%cverbatim:",DATA_BEGIN); fflush(stdout); }
  tm_is_waiting = 1;
}
/* tell TeXmacs GP is done and is waiting for new data */
static void
tm_end_output(void)
{
  if (tm_is_waiting) { printf("%c", DATA_END); fflush(stdout); }
  tm_is_waiting = 0;
}
static char *
fgets_texmacs(char *s, int n, FILE *f)
{
  if (!tm_did_complete)
  {
    tm_start_output(); tm_end_output(); /* tell TeXmacs we need input */
  }
  return fgets(s,n,f);
}

#ifdef READLINE
typedef struct {
  char *cmd;
  long n; /* number of args */
  char **v; /* args */
} tm_cmd;

static void
parse_texmacs_command(tm_cmd *c, const char *ch)
{
  long l = strlen(ch);
  char *t, *s = (char*)ch, *send = s+l-1;
  growarray A;

  if (*s != DATA_BEGIN || *send-- != DATA_END)
    pari_err(talker, "missing DATA_[BEGIN | END] in TeXmacs command");
  s++;
  if (strncmp(s, "special:", 8)) pari_err(talker, "unrecognized TeXmacs command");
  s += 8;
  if (*s != '(' || *send-- != ')')
    pari_err(talker, "missing enclosing parentheses for TeXmacs command");
  s++; t = s;
  skip_alpha(s);
  c->cmd = pari_strndup(t, s - t);
  grow_init(A);
  for (c->n = 0; s <= send; c->n++)
  {
    char *u = gpmalloc(strlen(s) + 1);
    skip_space(s);
    if (*s == '"') s = readstring(s, u);
    else
    { /* read integer */
      t = s;
      while (isdigit((int)*s)) s++;
      strncpy(u, t, s - t); u[s-t] = 0;
    }
    grow_append(A, (void*)u);
  }
  c->v = (char**)A->v;
}

static void
free_cmd(tm_cmd *c)
{
  while (c->n--) free((void*)c->v[c->n]);
  free((void*)c->v);
}

static void
handle_texmacs_command(const char *s)
{
  tm_cmd c;
  parse_texmacs_command(&c, s);
  if (strcmp(c.cmd, "complete"))
    pari_err(talker,"Texmacs_stdin command %s not implemented", c.cmd);
  if (c.n != 2)
    pari_err(talker,"was expecting 2 arguments for Texmacs_stdin command");
  texmacs_completion(c.v[0], atol(c.v[1]));
  free_cmd(&c);
  tm_did_complete = 1;
}
#else
static void
handle_texmacs_command(const char *s) { pari_err(talker, "readline not available"); }
#endif

/*******************************************************************/
/**                                                               **/
/**                          BUFFERS                              **/
/**                                                               **/
/*******************************************************************/
#define current_buffer (bufstack?((Buffer*)(bufstack->value)):NULL)
static stack *bufstack = NULL;

static void
pop_buffer(void)
{
  Buffer *b = (Buffer*) pop_stack(&bufstack);
  delete_buffer(b);
}

/* kill all buffers until B is met or nothing is left */
static void
kill_all_buffers(Buffer *B)
{
  for(;;) {
    Buffer *b = current_buffer;
    if (b == B || !b) break;
    pop_buffer();
  }
}

static void
jump_to_given_buffer(Buffer *buf)
{
  Buffer *b;
  while ( (b = current_buffer) )
  {
    if (b == buf) break;
    pop_buffer();
  }
  if (!b || !b->env) longjmp(GP_DATA->env, 0);
  longjmp(b->env, 0);
}

/********************************************************************/
/**                                                                **/
/**                             HELP                               **/
/**                                                                **/
/********************************************************************/
static int
has_ext_help(void)
{
  if (GP_DATA->help)
  {
    char *buf = pari_strdup(GP_DATA->help), *s, *t;
    FILE *file;

    for (t = s = buf; *s; *t++ = *s++)
    {
      if (*s == '\\') s++; else if (*s == ' ') break;
    }
    *t = 0; file = fopen(buf,"r");
    free(buf);
    if (file) { fclose(file); return 1; }
  }
  return 0;
}

static int
compare_str(char **s1, char **s2) { return strcmp(*s1, *s2); }

/* Print all elements of list in columns, pausing every nbli lines
 * if nbli is non-zero.
 * list is a NULL terminated list of function names
 */
void
print_fun_list(char **list, long nbli)
{
  long i=0, j=0, maxlen=0, nbcol,len, w = term_width();
  char **l;

  while (list[i]) i++;
  qsort (list, i, sizeof(char *), (QSCOMP)compare_str);

  for (l=list; *l; l++)
  {
    len = strlen(*l);
    if (len > maxlen) maxlen=len;
  }
  maxlen++; nbcol= w / maxlen;
  if (nbcol * maxlen == w) nbcol--;
  if (!nbcol) nbcol = 1;

  pariputc('\n'); i=0;
  for (l=list; *l; l++)
  {
    pariputs(*l); i++;
    if (i >= nbcol)
    {
      i=0; pariputc('\n');
      if (nbli && j++ > nbli) { j = 0; hit_return(); }
      continue;
    }
    len = maxlen - strlen(*l);
    while (len--) pariputc(' ');
  }
  if (i) pariputc('\n');
}

static void
commands(long n)
{
  const size_t LIST_LEN = 1023UL;
  size_t size = LIST_LEN, s = 0;
  long i;
  entree *ep;
  char **list = (char **) gpmalloc((size+1)*sizeof(char *));

  for (i = 0; i < functions_tblsz; i++)
    for (ep = functions_hash[i]; ep; ep = ep->next)
      if ((n < 0 && ep->menu) || ep->menu == n)
      {
        list[s] = ep->name;
        if (++s >= size)
        {
	  size += (LIST_LEN+1);
          list = (char**) gprealloc(list, size*sizeof(char *));
        }
      }
  list[s] = NULL;
  print_fun_list(list, term_height()-4);
  free(list);
}

static void
center(char *s)
{
  long i, l = strlen(s), pad = term_width() - l;
  char *buf, *u;

  if (pad<0) pad=0; else pad >>= 1;
  u = buf = (char*)gpmalloc(l + pad + 2);
  for (i=0; i<pad; i++) *u++ = ' ';
  while (*s) *u++ = *s++;
  *u++ = '\n'; *u = 0;
  pariputs(buf); free(buf);
}

static void
usage(char *s)
{
  printf("### Usage: %s [options] [GP files]\n", s);
  printf("Options are:\n");
  printf("\t[-f,--fast]\tFaststart: do not read .gprc\n");
  printf("\t[-q,--quiet]\tQuiet mode: do not print banner and history numbers\n");
  printf("\t[-p,--primelimit primelimit]\n\t\t\tPrecalculate primes up to the limit\n");
  printf("\t[-s,--stacksize stacksize]\n\t\t\tStart with the PARI stack of given size (in bytes)\n");
  printf("\t[--emacs]\tRun as if in Emacs shell\n");
  printf("\t[--help]\tPrint this message\n");
  printf("\t[--test]\tTest mode. No history, wrap long lines (bench only)\n");
  printf("\t[--texmacs]\tRun as if using TeXmacs frontend\n");
  printf("\t[--version]\tOutput version info and exit\n");
  printf("\t[--version-short]\tOutput version number and exit\n\n");
  exit(0);
}

static void
community(void)
{
  pari_sp av = avma;
  char *s = stackmalloc(strlen(GPDATADIR) + 1024);

  (void)sprintf(s, "The standard distribution of GP/PARI includes a \
reference manual, a tutorial, a reference card and quite a few examples. They \
should have been installed in the directory '%s'. If not, ask the person \
who installed PARI on your system where they can be found. You can also \
download them from the PARI WWW site 'http://pari.math.u-bordeaux.fr/'",
GPDATADIR);
  print_text(s); avma = av;

  pariputs("\nThree mailing lists are devoted to PARI:\n\
  - pari-announce (moderated) to announce major version changes.\n\
  - pari-dev for everything related to the development of PARI, including\n\
    suggestions, technical questions, bug reports and patch submissions.\n\
  - pari-users for everything else!\n");
  print_text("\
To subscribe, send an empty message to <listname>-subscribe@list.cr.yp.to. \
An archive is kept at the WWW site mentioned above. You can also reach the \
authors directly by email: pari@math.u-bordeaux.fr (answer not guaranteed)."); }

static void
gentypes(void)
{
  pariputs("List of the PARI types:\n\
  t_INT    : long integers     [ cod1 ] [ cod2 ] [ man_1 ] ... [ man_k ]\n\
  t_REAL   : long real numbers [ cod1 ] [ cod2 ] [ man_1 ] ... [ man_k ]\n\
  t_INTMOD : integermods       [ code ] [ mod  ] [ integer ]\n\
  t_FRAC   : irred. rationals  [ code ] [ num. ] [ den. ]\n\
  t_COMPLEX: complex numbers   [ code ] [ real ] [ imag ]\n\
  t_PADIC  : p-adic numbers    [ cod1 ] [ cod2 ] [ p ] [ p^r ] [ int ]\n\
  t_QUAD   : quadratic numbers [ cod1 ] [ mod  ] [ real ] [ imag ]\n\
  t_POLMOD : poly mod          [ code ] [ mod  ] [ polynomial ]\n\
  -------------------------------------------------------------\n\
  t_POL    : polynomials       [ cod1 ] [ cod2 ] [ man_1 ] ... [ man_k ]\n\
  t_SER    : power series      [ cod1 ] [ cod2 ] [ man_1 ] ... [ man_k ]\n\
  t_RFRAC  : irred. rat. func. [ code ] [ num. ] [ den. ]\n\
  t_QFR    : real qfb          [ code ] [ a ] [ b ] [ c ] [ del ]\n\
  t_QFI    : imaginary qfb     [ code ] [ a ] [ b ] [ c ]\n\
  t_VEC    : row vector        [ code ] [  x_1  ] ... [  x_k  ]\n\
  t_COL    : column vector     [ code ] [  x_1  ] ... [  x_k  ]\n\
  t_MAT    : matrix            [ code ] [ col_1 ] ... [ col_k ]\n\
  t_LIST   : list              [ code ] [ cod2 ] [ x_1 ] ... [ x_k ]\n\
  t_STR    : string            [ code ] [ man_1 ] ... [ man_k ]\n\
  t_VECSMALL: vec. small ints  [ code ] [ x_1 ] ... [ x_k ]\n\
\n");
}

static void
menu_commands(void)
{
  pariputs("Help topics: for a list of relevant subtopics, type ?n for n in\n\
  0: user-defined identifiers (variable, alias, function)\n\
  1: Standard monadic or dyadic OPERATORS\n\
  2: CONVERSIONS and similar elementary functions\n\
  3: TRANSCENDENTAL functions\n\
  4: NUMBER THEORETICAL functions\n\
  5: Functions related to ELLIPTIC CURVES\n\
  6: Functions related to general NUMBER FIELDS\n\
  7: POLYNOMIALS and power series\n\
  8: Vectors, matrices, LINEAR ALGEBRA and sets\n\
  9: SUMS, products, integrals and similar functions\n\
 10: GRAPHIC functions\n\
 11: PROGRAMMING under GP\n\
 12: The PARI community\n\
\n\
Also:\n\
  ? functionname (short on-line help)\n\
  ?\\             (keyboard shortcuts)\n\
  ?.             (member functions)\n");
  if (has_ext_help()) pariputs("\
Extended help looks available:\n\
  ??             (opens the full user's manual in a dvi previewer)\n\
  ??  tutorial / refcard / libpari (tutorial/reference card/libpari manual)\n\
  ??  keyword    (long help text about \"keyword\" from the user's manual)\n\
  ??? keyword    (a propos: list of related functions).");
}

static void
slash_commands(void)
{
  pariputs("#       : enable/disable timer\n\
##      : print time for last result\n\
\\\\      : comment up to end of line\n\
\\a {n}  : print result in raw format (readable by PARI)\n\
\\b {n}  : print result in beautified format\n\
\\c      : list all commands (same effect as ?*)\n\
\\d      : print all defaults\n\
\\e {n}  : enable/disable echo (set echo=n)\n\
\\g {n}  : set debugging level\n\
\\gf{n}  : set file debugging level\n\
\\gm{n}  : set memory debugging level\n\
\\h {m-n}: hashtable information\n\
\\l {f}  : enable/disable logfile (set logfile=f)\n\
\\m {n}  : print result in prettymatrix format\n\
\\o {n}  : change output method (0=raw, 1=prettymatrix, 2=prettyprint, 3=2-dim)\n\
\\p {n}  : change real precision\n\
\\ps{n}  : change series precision\n\
\\q      : quit completely this GP session\n\
\\r {f}  : read in a file\n\
\\s {n}  : print stack information\n\
\\t      : print the list of PARI types\n\
\\u      : print the list of user-defined functions\n\
\\um     : print the list of user-defined member functions\n\
\\v      : print current version of GP\n\
\\w {nf} : write to a file\n\
\\x {n}  : print complete inner structure of result\n\
\\y {n}  : disable/enable automatic simplification (set simplify=n)\n\
\n\
{f}=optional filename. {n}=optional integer\n");
}

static void
member_commands(void)
{
  pariputs("\
Member functions, followed by relevant objects\n\n\
a1-a6, b2-b8, c4-c6 : coeff. of the curve.            ell\n\
area : area                                           ell\n\
bid  : big ideal                                                    bnr\n\
bnf  : big number field                                        bnf, bnr\n\
clgp : class group                   bid,                      bnf, bnr\n\
cyc  : cyclic decomposition (SNF)    bid,       clgp,          bnf, bnr\n\
diff, codiff: different and codifferent                    nf, bnf, bnr\n\
disc : discriminant                                   ell, nf, bnf, bnr\n\
e, f : inertia/residue  degree            prid\n\
fu   : fundamental units                                       bnf, bnr\n\
gen  : generators                    bid, prid, clgp,          bnf, bnr\n\
index: index                                               nf, bnf, bnr\n\
j    : j-invariant                                    ell\n");
/* split: some compilers can't handle long constant strings */
  pariputs("\
mod  : modulus                       bid,                           bnr\n\
nf   : number field                                            bnf, bnr\n\
no   : number of elements            bid,       clgp,          bnf, bnr\n\
omega, eta: [omega1,omega2] and [eta1, eta2]          ell\n\
p    : rational prime below prid          prid\n\
pol  : defining polynomial                                 nf, bnf, bnr\n\
reg  : regulator                                               bnf, bnr\n\
roots: roots                                          ell  nf, bnf, bnr\n\
sign,r1,r2 : signature                                     nf, bnf, bnr\n\
t2   : t2 matrix                                           nf, bnf, bnr\n\
tate : Tate's [u^2, u, q]                             ell\n\
tu   : torsion unit and its order                              bnf, bnr\n\
w    : Mestre's w                                     ell\n\
zk   : integral basis                                      nf, bnf, bnr\n");
}

#define QUOTE "_QUOTE"
#define DOUBQUOTE "_DOUBQUOTE"
#define BACKQUOTE "_BACKQUOTE"

static char *
_cat(char *s, char *t)
{
  *s = 0; strcat(s,t); return s + strlen(t);
}

static char *
filter_quotes(char *s)
{
  int i, l = strlen(s);
  int quote = 0;
  int backquote = 0;
  int doubquote = 0;
  char *str, *t;

  for (i=0; i < l; i++)
    switch(s[i])
    {
      case '\'': quote++; break;
      case '`' : backquote++; break;
      case '"' : doubquote++;
    }
  str = (char*)gpmalloc(l + quote * (strlen(QUOTE)-1)
                          + doubquote * (strlen(DOUBQUOTE)-1)
                          + backquote * (strlen(BACKQUOTE)-1) + 1);
  t = str;
  for (i=0; i < l; i++)
    switch(s[i])
    {
      case '\'': t = _cat(t, QUOTE); break;
      case '`' : t = _cat(t, BACKQUOTE); break;
      case '"' : t = _cat(t, DOUBQUOTE); break;
      default: *t++ = s[i];
    }
  *t = 0; return str;
}

static int
nl_read(char *s) { size_t l = strlen(s); return s[l-1] == '\n'; }

#define nbof(a) sizeof(a) / sizeof(a[0])
/* query external help program for s. num < 0 [keyword] or chapter number */
static void
external_help(char *s, int num)
{
  long nbli = term_height()-3, li = 0;
  char buf[256], ar[32], *str, *opt = "";
  pariFILE *z;
  FILE *f;

  if (!GP_DATA->help) pari_err(talker,"no external help program");
  s = filter_quotes(s);
  str = gpmalloc(strlen(GP_DATA->help) + strlen(s) + 64);
  *ar = 0;
  if (num < 0)
    opt = "-k";
  else if (s[strlen(s)-1] != '@')
    sprintf(ar,"@%d",num);
  sprintf(str,"%s -fromgp %s %c%s%s%c",GP_DATA->help,opt, SHELL_Q,s,ar,SHELL_Q);
  z = try_pipe(str,0); f = z->file;
  free(str);
  free(s);
  while (fgets(buf, nbof(buf), f))
  {
    if (!strncmp("ugly_kludge_done",buf,16)) break;
    pariputs(buf);
    if (nl_read(buf) && ++li > nbli) { hit_return(); li = 0; }
  }
  pari_fclose(z);
}

char *keyword_list[]={
  "operator",
  "libpari",
  "member",
  "integer",
  "real",
  "readline",
  "refcard",
  "tutorial",
  "nf",
  "bnf",
  "bnr",
  "ell",
  "rnf",
  "bid",
  "modulus",
  NULL
};

static int
ok_external_help(char **s)
{
  long n;
  if (!**s) return 1;
  if (!isalpha((int)**s)) return 3; /* operator or section number */
  if (!strncmp(*s,"t_",2)) { *s += 2; return 2; } /* type name */

  for (n=0; keyword_list[n]; n++)
    if (!strcmp(*s,keyword_list[n])) return 3;
  return 0;
}

/* don't mess readline display */
static void
aide_print(char *s1, char *s2) { pariprintf("%s: %s\n", s1, s2); }

static void
aide0(char *s, int flag)
{
  long n, long_help = flag & h_LONG;
  entree *ep,*ep1;
  char *s1;

  s = get_sep(s);
  if (isdigit((int)*s))
  {
    n = atoi(s);
    if (n == 12) { community(); return; }
    if (n < 0 || n > 12)
      pari_err(talker2,"no such section in help: ?",s,s);
    if (long_help) external_help(s,3); else commands(n);
    return;
  }
  /* Get meaningful entry on \ps 5 */
  if (*s == '\\') { s1 = s+1; skip_alpha(s1); *s1 = '\0';}

  if (flag & h_APROPOS) { external_help(s,-1); return; }
  if (long_help && (n = ok_external_help(&s))) { external_help(s,n); return; }
  switch (*s)
  {
    case '*' : commands(-1); return;
    case '\0': menu_commands(); return;
    case '\\': slash_commands(); return;
    case '.' : member_commands(); return;
  }
  ep = is_entry(s);
  if (ep && long_help)
  {
    if (!strcmp(ep->name, "default"))
    {
      char *t = s+7, *e;
      skip_space(t);
      if (*t == '(') {
	t++; skip_space(t);
        e = t; skip_alpha(e); *e = '\0'; /* safe: get_sep() made it a copy */
	if (is_default(t)) { external_help(t, 2); return; }
      }
    }
  }
  if (!ep)
  {
    n = is_default(s)? 2: 3;
    if (long_help)
      external_help(s,n);
    else
    {
      if (n == 2) { aide_print(s,"default"); return; }
      n = whatnow(s,1);
      if (!n) { aide_print(s,"unknown identifier"); return; }
      aide_print(s, "obsolete function");
      whatnow_new_syntax(s, n);
    }
    return;
  }

  ep1 = ep;  ep = do_alias(ep);
  if (ep1 != ep) pariprintf("%s is aliased to:\n\n",s);

  switch(EpVALENCE(ep))
  {
    case EpUSER:
      if (!ep->help || long_help) print_user_fun(ep);
      if (!ep->help) return;
      if (long_help) { pariputs("\n\n"); long_help=0; }
      break;

    case EpGVAR:
    case EpVAR:
      if (!ep->help) { aide_print(s, "user defined variable"); return; }
      long_help=0; break;

    case EpINSTALL:
      if (!ep->help) { aide_print(s, "installed function"); return; }
      long_help=0; break;

    case EpNEW: aide_print(s, "new identifier (no valence assigned)"); return;
  }
  if (long_help) { external_help(ep->name,3); return; }
  if (ep->help) { print_text(ep->help); return; }

  pari_err(bugparier,"aide (no help found)");
}

void
aide(char *s, long flag)
{
  if ((flag & h_RL) == 0)
  {
    if (*s == '?') { flag |= h_LONG; s++; }
    if (*s == '?') { flag |= h_APROPOS; s++; }
  }
  term_color(c_HELP); aide0(s,flag); term_color(c_NONE);
  if ((flag & h_RL) == 0) pariputc('\n');
}

/********************************************************************/
/**                                                                **/
/**                         GP HEADER                              **/
/**                                                                **/
/********************************************************************/
static char *
what_readline()
{
  char *s;
#ifdef READLINE
  char *ver, *extra = stackmalloc(strlen(READLINE) + 32);
#  if defined(HAS_RL_LIBRARY_VERSION) || defined(FAKE_RL_LIBRARY_VERSION)
#    ifdef FAKE_RL_LIBRARY_VERSION
  extern char *rl_library_version;
#    endif

  if (strcmp(READLINE, rl_library_version))
  {
    ver = (char*)rl_library_version;
    (void)sprintf(extra, " [was v%s in Configure]", READLINE);
  }
  else
#  endif
  {
    ver = READLINE;
    extra[0] = 0;
  }
  s = stackmalloc(3 + strlen(ver) + 8 + strlen(extra));
  (void)sprintf(s, "v%s %s%s", ver,
            (GP_DATA->flags & USE_READLINE)? "enabled": "disabled",
            extra);
#else
  s = "not compiled in";
#endif
  return s;
}

static void
print_shortversion(void)
{
  const ulong mask = (1<<PARI_VERSION_SHIFT) - 1;
  ulong n = PARI_VERSION_CODE, major, minor, patch;

  patch = n & mask; n >>= PARI_VERSION_SHIFT;
  minor = n & mask; n >>= PARI_VERSION_SHIFT;
  major = n;
  pariprintf("%lu.%lu.%lu\n", major,minor,patch); exit(0);
}

static char *
what_cc()
{
  char *s;
#ifdef GCC_VERSION
#  ifdef __cplusplus
#    define Format "g++-%s"
#  else
#    define Format "gcc-%s"
#  endif
  s = stackmalloc(4 + strlen(GCC_VERSION) + 1);
  (void)sprintf(s, Format, GCC_VERSION);
#else
#  ifdef _MSC_VER
  s = stackmalloc(32);
  (void)sprintf(s, "MSVC-%i", _MSC_VER);
#  else
  s = NULL;
#  endif
#endif
  return s;
}

static void
print_version(void)
{
  pari_sp av = avma;
  char *buf, *ver;

  center(PARIVERSION);
  center(PARIINFO);
  ver = what_cc();
  buf = stackmalloc(strlen(__DATE__) +  32 + (ver? strlen(ver): 0));
  if (ver) (void)sprintf(buf, "compiled: %s, %s", __DATE__, ver);
  else     (void)sprintf(buf, "compiled: %s", __DATE__);
  center(buf);
  ver = what_readline();
  buf = stackmalloc(strlen(ver) + 64);
  (void)sprintf(buf, "(readline %s, extended help%s available)", ver,
                has_ext_help()? "": " not");
  center(buf); avma = av;
}

static void
gp_head(void)
{
#ifdef READLINE
  if (GP_DATA->flags & TEXMACS)
    printf("%ccommand:(cas-supports-completions-set! \"pari\")%c\n",
           DATA_BEGIN, DATA_END);
#endif
  print_version(); pariputs("\n");
  center("Copyright (C) 2000-2006 The PARI Group");
  print_text("\nPARI/GP is free software, covered by the GNU General Public \
License, and comes WITHOUT ANY WARRANTY WHATSOEVER");
  pariputs("\n\
Type ? for help, \\q to quit.\n\
Type ?12 for how to get moral (and possibly technical) support.\n\n");
  pariprintf("parisize = %lu, primelimit = %lu\n", top-bot, GP_DATA->primelimit);
}

/********************************************************************/
/**                                                                **/
/**                         METACOMMANDS                           **/
/**                                                                **/
/********************************************************************/
#define pariputs_opt(s) if (!(GP_DATA->flags & QUIET)) pariputs(s)

void
gp_quit(void)
{
  free_graph();
  pari_close();
  kill_all_buffers(NULL);
  term_color(c_NONE);
  pariputs_opt("Goodbye!\n");
  if (GP_DATA->flags & TEXMACS) tm_end_output();
  exit(0);
}

static GEN
gpreadbin(const char *s, int *vector)
{
  GEN x = readbin(s,infile,vector);
  popinfile(); return x;
}

static void
escape0(char *tch)
{
  const char *s;
  char c;

  if (compatible != NONE)
  {
    s = tch;
    while (*s)
      if (*s++ == '=')
      {
        char *f = NULL;
	long len = (s-tch) - 1;
	if      (!strncmp(tch,"precision",len))    f = "realprecision";
	else if (!strncmp(tch,"serieslength",len)) f = "seriesprecision";
	else if (!strncmp(tch,"format",len))       f = "format";
	else if (!strncmp(tch,"prompt",len))       f = "prompt";
	if (f) { (void)setdefault(f, s, d_ACKNOWLEDGE); return; }
	break;
      }
  }
  s = tch;
  switch ((c = *s++))
  {
    case 'w': case 'x': case 'a': case 'b': case 'B': case 'm':
    { /* history things */
      long d;
      GEN x;
      if (c != 'w' && c != 'x') d = get_int(s,0);
      else
      {
	d = atol(s); if (*s == '-') s++;
	while (isdigit((int)*s)) s++;
      }
      x = gp_history(GP_DATA->hist, d, tch+1,tch-1);
      switch (c)
      {
	case 'B':
        { /* prettyprinter */
          gp_data G = *GP_DATA; /* copy */
          gp_hist   h = *(G.hist); /* copy */
          pariout_t f = *(G.fmt);  /* copy */

          G.hist = &h; h.total = 0; /* no hist number */
          G.fmt  = &f; f.prettyp = f_PRETTY;
          G.flags &= ~(TEST|TEXMACS);
          G.lim_lines = 0;
          gp_output(x, &G); break;
        }
	case 'a': brute   (x, GP_DATA->fmt->format, -1); break;
	case 'm': matbrute(x, GP_DATA->fmt->format, -1); break;
	case 'b': sor(x, GP_DATA->fmt->format, -1, GP_DATA->fmt->fieldw); break;
	case 'x': voir(x, get_int(s, -1)); break;
        case 'w':
	  s = get_sep(s); if (!*s) s = current_logfile;
	  write0(s, mkvec(x)); return;
      }
      pariputc('\n'); return;
    }

    case 'c': commands(-1); break;
    case 'd': (void)setdefault("",NULL,0); break;
    case 'e':
      s = get_sep(s);
      if (!*s) s = (GP_DATA->flags & ECHO)? "0": "1";
      (void)sd_echo(s,d_ACKNOWLEDGE); break;
    case 'g':
      switch (*s)
      {
        case 'm': (void)sd_debugmem(++s,d_ACKNOWLEDGE); break;
        case 'f': (void)sd_debugfiles(++s,d_ACKNOWLEDGE); break;
        default : (void)sd_debug(s,d_ACKNOWLEDGE); break;
      }
      break;
    case 'h': print_functions_hash(s); break;
    case 'l':
      s = get_sep(s);
      if (*s)
      {
        (void)sd_logfile(s,d_ACKNOWLEDGE);
        if (logfile) break;
      }
      (void)sd_log(logfile?"0":"1",d_ACKNOWLEDGE);
      break;
    case 'o': (void)sd_output(s,d_ACKNOWLEDGE); break;
    case 'p':
      switch (*s)
      {
        case 's': (void)sd_seriesprecision(++s,d_ACKNOWLEDGE); break;
        default : (void)sd_realprecision(s,d_ACKNOWLEDGE); break;
      }
      break;
    case 'q': gp_quit(); break;
    case 'r':
      s = get_sep(s);
      switchin(s);
      if (file_is_binary(infile))
      {
        int vector;
        GEN x = gpreadbin(s, &vector);
        if (vector) /* many BIN_GEN */
        {
          long i, l = lg(x);
          pari_warn(warner,"setting %ld history entries", l-1);
          for (i=1; i<l; i++) (void)set_hist_entry(GP_DATA->hist, (GEN)x[i]);
        }
      }
      break;
    case 's': etatpile(); break;
    case 't': gentypes(); break;
    case 'u':
      switch (*s)
      {
        case 'm': print_all_user_member(); break;
        default: print_all_user_fun();
      }
      break;
    case 'v': print_version(); break;
    case 'y':
      s = get_sep(s);
      if (!*s) s = (GP_DATA->flags & SIMPLIFY)? "0": "1";
      (void)sd_simplify(s,d_ACKNOWLEDGE); break;
    default: pari_err(caracer1,tch,tch-1);
  }
}

static void
escape(char *tch)
{
  char *old = get_analyseur();
  set_analyseur(tch); /* for error messages */
  escape0(tch);
  set_analyseur(old);
}

enum { ti_NOPRINT, ti_REGULAR, ti_LAST, ti_INTERRUPT };
/* flag:
 *   ti_NOPRINT   don't print
 *   ti_REGULAR   print elapsed time (flags & CHRONO)
 *   ti_LAST      print last elapsed time (##)
 *   ti_INTERRUPT received a SIGINT
 */
static char *
gp_format_time(long flag)
{
  static char buf[64];
  static long last = 0;
  long delay = (flag == ti_LAST)? last: TIMER(GP_DATA->T);
  char *s;

  last = delay;
  switch(flag)
  {
    case ti_REGULAR:   s = "time = "; break;
    case ti_INTERRUPT: s = "user interrupt after "; break;
    case ti_LAST:      s = "  ***   last result computed in "; break;
    default: return NULL;
  }
  strcpy(buf,s); s = buf+strlen(s);
  strcpy(s, term_get_color(c_TIME)); s+=strlen(s);
  if (delay >= 3600000)
  {
    sprintf(s, "%ldh, ", delay / 3600000); s+=strlen(s);
    delay %= 3600000;
  }
  if (delay >= 60000)
  {
    sprintf(s, "%ldmn, ", delay / 60000); s+=strlen(s);
    delay %= 60000;
  }
  if (delay >= 1000)
  {
    sprintf(s, "%ld,", delay / 1000); s+=strlen(s);
    delay %= 1000;
    if (delay < 100)
    {
      sprintf(s, "%s", (delay<10)? "00": "0");
      s+=strlen(s);
    }
  }
  sprintf(s, "%ld ms", delay); s+=strlen(s);
  strcpy(s, term_get_color(c_NONE));
  if (flag != ti_INTERRUPT) { s+=strlen(s); *s++='.'; *s++='\n'; *s=0; }
  return buf;
}

static int
chron(char *s)
{
  if (*s)
  { /* if "#" or "##" timer metacommand. Otherwise let the parser get it */
    if (*s == '#') s++;
    if (*s) return 0;
    pariputs(gp_format_time(ti_LAST));
  }
  else { GP_DATA->flags ^= CHRONO; (void)sd_timer("",d_ACKNOWLEDGE); }
  return 1;
}

/* return 0: can't interpret *buf as a metacommand
 *        1: did interpret *buf as a metacommand or empty command */
static int
check_meta(char *buf)
{
  switch(*buf++)
  {
    case '?': aide(buf, h_REGULAR); break;
    case '#': return chron(buf);
    case '\\': escape(buf); break;
    case '\0': break;
    default: return 0;
  }
  return 1;
}

/********************************************************************/
/*                                                                  */
/*                              GPRC                                */
/*                                                                  */
/********************************************************************/
#if defined(UNIX) || defined(__EMX__)
#  include <pwd.h>
#endif

static int get_line_from_file(char *prompt, filtre_t *F, FILE *file);
#define err_gprc(s,t,u) { fprintferr("\n"); pari_err(talker2,s,t,u); }

/* LOCATE GPRC */

/* return $HOME or the closest we can find */
static char *
get_home(int *free_it)
{
  char *drv, *pth = os_getenv("HOME");
  if (pth) return pth;
  if ((drv = os_getenv("HOMEDRIVE"))
   && (pth = os_getenv("HOMEPATH")))
  { /* looks like WinNT */
    char *buf = gpmalloc(strlen(pth) + strlen(drv) + 1);
    sprintf(buf, "%s%s",drv,pth);
    *free_it = 1; return buf;
  }
#if defined(__EMX__) || defined(UNIX)
  {
    struct passwd *p = getpwuid(geteuid());
    if (p) return p->pw_dir;
  }
#endif
  return ".";
}

static FILE *
gprc_chk(char *s)
{
  FILE *f = fopen(s, "r");
  if (f && !(GP_DATA->flags & QUIET)) fprintferr("Reading GPRC: %s ...", s);
  return f;
}

/* Look for [._]gprc: $GPRC, then in $HOME, ., /etc, path [ to gp binary ] */
static FILE *
gprc_get(char *path)
{
  FILE *f = NULL;
  char *str, *s, c;
  long l;
  s = os_getenv("GPRC");
  if (s) f = gprc_chk(s);
  if (!f)
  {
    int free_it = 0;
    s = get_home(&free_it); l = strlen(s); c = s[l-1];
    str = strcpy(gpmalloc(l+7), s);
    if (free_it) free(s);
    s = str + l;
    if (c != '/' && c != '\\') *s++ = '/';
#ifdef UNIX
    *s = '.'; /* .gprc */
#else
    *s = '_'; /* _gprc */
#endif
    strcpy(s+1, "gprc");
    f = gprc_chk(str); /* in $HOME */
    if (!f) f = gprc_chk(s); /* in . */
    if (!f) f = gprc_chk("/etc/gprc");
    if (!f) f = gprc_chk("C:/_gprc");
    if (!f)
    { /* in 'gp' directory */
      char *t = path + strlen(path);
      while (t > path && *t != '/') t--;
      if (*t == '/')
      {
        long l = t - path + 1;
        t = gpmalloc(l + 6);
        strncpy(t, path, l);
        strcpy(t+l, s); f = gprc_chk(t);
        free(t);
      }
    }
    free(str);
  }
  return f;
}

/* PREPROCESSOR */

static ulong
read_uint(char **s)
{
  long v = atol(*s);
  if (!isdigit((int)**s)) err_gprc("not an integer", *s, *s);
  while (isdigit((int)**s)) (*s)++;
  return v;
}
static ulong
read_dot_uint(char **s)
{
  if (**s != '.') return 0;
  (*s)++; return read_uint(s);
}
/* read a.b.c */
static long
read_version(char **s)
{
  long a, b, c;
  a = read_uint(s);
  b = read_dot_uint(s);
  c = read_dot_uint(s);
  return PARI_VERSION(a,b,c);
}

static int
get_preproc_value(char **s)
{
  if (!strncmp(*s,"EMACS",5)) {
    *s += 5;
    return GP_DATA->flags & (EMACS|TEXMACS);
  }
  if (!strncmp(*s,"READL",5)) {
    *s += 5;
    return GP_DATA->flags & USE_READLINE;
  }
  if (!strncmp(*s,"VERSION",7)) {
    int less = 0, orequal = 0;
    long d;
    *s += 7;
    switch(**s)
    {
      case '<': (*s)++; less = 1; break;
      case '>': (*s)++; less = 0; break;
      default: return -1;
    }
    if (**s == '=') { (*s)++; orequal = 1; }
    d = PARI_VERSION_CODE - read_version(s);
    if (!d) return orequal;
    return less? (d < 0): (d > 0);
  }
  return -1;
}

/* PARSE GPRC */

/* 1) replace next separator by '\0' (t must be writeable)
 * 2) return the next expression ("" if none)
 * see get_sep() */
static char *
next_expr(char *t)
{
  int outer = 1;
  char *s = t;

  for(;;)
  {
    char c;
    switch ((c = *s++))
    {
      case '"':
        if (outer || (s >= t+2 && s[-2] != '\\')) outer = !outer;
        break;
      case '\0':
        return "";
      default:
        if (outer && c == ';') { s[-1] = 0; return s; }
    }
  }
}

static void
gp_initrc(growarray A, char *path)
{
  char *nexts,*s,*t;
  FILE *file = gprc_get(path);
  Buffer *b;
  filtre_t F;
  VOLATILE long c = 0;

  if (!file) return;
  b = new_buffer();
  init_filtre(&F, b);
  for(;;)
  {
    if (setjmp(GP_DATA->env)) fprintferr("...skipping line %ld.\n", c);
    c++;
    if (!get_line_from_file(NULL,&F,file)) break;
    s = b->buf;
    if (*s == '#')
    { /* preprocessor directive */
      int z, NOT = 0;
      s++;
      if (strncmp(s,"if",2)) err_gprc("unknown directive",s,b->buf);
      s += 2;
      if (!strncmp(s,"not",3)) { NOT = !NOT; s += 3; }
      if (*s == '!')           { NOT = !NOT; s++; }
      t = s;
      z = get_preproc_value(&s);
      if (z < 0) err_gprc("unknown preprocessor variable",t,b->buf);
      if (NOT) z = !z;
      if (!*s)
      { /* make sure at least an expr follows the directive */
        if (!get_line_from_file(NULL,&F,file)) break;
        s = b->buf;
      }
      if (!z) continue; /* dump current line */
    }
    /* parse line */
    for ( ; *s; s = nexts)
    {
      nexts = next_expr(s);
      if (!strncmp(s,"read",4) && (s[4] == ' ' || s[4] == '\t' || s[4] == '"'))
      { /* read file */
	s += 4;
	t = gpmalloc(strlen(s) + 1);
	if (*s == '"') (void)readstring(s, t); else strcpy(t,s);
        grow_append(A, t);
      }
      else
      { /* set default */
	t = s; while (*t && *t != '=') t++;
	if (*t != '=') err_gprc("missing '='",t,b->buf);
	*t++ = 0;
	if (*t == '"') (void)readstring(t, t);
	(void)setdefault(s,t,d_INITRC);
      }
    }
  }
  delete_buffer(b);
  if (!(GP_DATA->flags & QUIET)) fprintferr("Done.\n\n");
  fclose(file);
}

/********************************************************************/
/*                                                                  */
/*                           GP MAIN LOOP                           */
/*                                                                  */
/********************************************************************/
static void
brace_color(char *s, int c, int force)
{
  if (disable_color || (gp_colors[c] == c_NONE && !force)) return;
#ifdef RL_PROMPT_START_IGNORE
  if (GP_DATA->flags & USE_READLINE)
    *s++ = RL_PROMPT_START_IGNORE;
#endif
  strcpy(s, term_get_color(c));
#ifdef RL_PROMPT_START_IGNORE
  if (GP_DATA->flags & USE_READLINE)
  {
    s+=strlen(s);
    *s++ = RL_PROMPT_END_IGNORE;
    *s = 0;
  }
#endif
}

char *
color_prompt(char *prompt)
{
  static char buf[MAX_PROMPT_LEN + 24]; /* + room for color codes */
  char *s;

  if (GP_DATA->flags & TEST) return prompt;
  s = buf; *s = 0;
  /* escape sequences bug readline, so use special bracing (if available) */
  brace_color(s, c_PROMPT, 0);
  s += strlen(s); strcpy(s, prompt);
  s += strlen(s);
  brace_color(s, c_INPUT, 1); return buf;
}

void
update_logfile(const char *prompt, const char *s)
{
  switch (logstyle) {
  case logstyle_TeX:
    fprintf(logfile,
            "\\PARIpromptSTART|%s\\PARIpromptEND|%s\\PARIinputEND|%%\n",
            prompt,s);
    break;
  case logstyle_plain:
    fprintf(logfile, "%s%s\n",prompt,s);
    break;
  case logstyle_color:
    /* Can't do in one pass, since term_get_color() returns a static */
    fprintf(logfile, "%s%s", term_get_color(c_PROMPT), prompt);
    fprintf(logfile, "%s%s", term_get_color(c_INPUT), s);
    fprintf(logfile, "%s\n", term_get_color(c_NONE));
    break;
  }
}

/* prompt = NULL --> from gprc. Return 1 if new input, and 0 if EOF */
static int
get_line_from_file(char *PROMPT, filtre_t *F, FILE *file)
{
  const int Texmacs_stdin = ((GP_DATA->flags & TEXMACS) && file == stdin);
  char *s;
  input_method IM;

  IM.file = file;
  IM.fgets= Texmacs_stdin? &fgets_texmacs: &fgets;
  IM.prompt = NULL;
  IM.getline= &file_input;
  IM.free = 0;
  if (! input_loop(F,&IM))
  {
    if (Texmacs_stdin) tm_start_output();
    return 0;
  }

  s = ((Buffer*)F->buf)->buf;
  if (*s && PROMPT) /* don't echo if from gprc */
  {
    if (GP_DATA->flags & ECHO)
      { pariputs(PROMPT); pariputs(s); pariputc('\n'); }
    else if (logfile)
      update_logfile(PROMPT, s);
    pariflush();
  }
  if (GP_DATA->flags & TEXMACS)
  {
    tm_did_complete = 0;
    if (Texmacs_stdin && *s == DATA_BEGIN)
    { handle_texmacs_command(s); *s = 0; }
    else
      tm_start_output();
  }
  return 1;
}

static int
is_interactive(void)
{
  ulong f = GP_DATA->flags;
#if defined(UNIX) || defined(__EMX__) || defined(_WIN32)
  return (infile == stdin && !(f & TEXMACS)
                          && (f & EMACS || isatty(fileno(stdin))));
#else
  return (infile == stdin && !(f & TEXMACS));
#endif
}

/* return 0 if no line could be read (EOF) */
static int
gp_read_line(filtre_t *F, char *PROMPT)
{
  Buffer *b = (Buffer*)F->buf;
  int res;
  F->downcase = (compatible == OLDALL);
  if (b->len > 100000) fix_buffer(b, 100000);
  if (is_interactive())
  {
#ifdef READLINE
    if (GP_DATA->flags & USE_READLINE)
      res = get_line_from_readline(PROMPT? PROMPT: GP_DATA->prompt,
                                   GP_DATA->prompt_cont, F);
    else
#endif
    {
      if (!PROMPT) PROMPT = color_prompt( expand_prompt(GP_DATA->prompt, F) );
      pariputs(PROMPT);
      res = get_line_from_file(PROMPT,F,infile);
    }
    if (!disable_color) { term_color(c_NONE); pariflush(); }
  }
  else
    res = get_line_from_file(DFT_PROMPT,F,infile);
  return res;
}

/* kill all history entries since loc */
static void
prune_history(gp_hist *H, long loc)
{
  long i, j;
  i = (H->total-1) % H->size;
  j = H->total - loc;
  for ( ; j > 0; i--,j--)
  {
    if (H->res[i])
    {
      gunclone(H->res[i]);
      H->res[i] = NULL;
    }
    if (!i) i = H->size;
  }
  H->total = loc;
}

static int
is_silent(char *s) { return s[strlen(s) - 1] == ';'; }

/* If there are other buffers open (bufstack != NULL), we are doing an
 * immediate read (with read, extern...) */
static GEN
gp_main_loop(int ismain)
{
  gp_hist *H  = GP_DATA->hist;
  VOLATILE GEN z = gnil;
  VOLATILE pari_sp av = top - avma;
  Buffer *b = new_buffer();
  filtre_t F;

  init_filtre(&F, b);
  push_stack(&bufstack, (void*)b);
  for(;;)
  {
    if (ismain)
    {
      static long tloc, outtyp;
      tloc = H->total;
      outtyp = GP_DATA->fmt->prettyp;
      recover(0);
      if (setjmp(GP_DATA->env))
      { /* recover from error */
        char *s = (char*)global_err_data;
        if (s && *s) outerr(readseq(s));
	avma = top; av = 0;
        prune_history(H, tloc);
        GP_DATA->fmt->prettyp = outtyp;
        kill_all_buffers(b);
      }
    }

    if (! gp_read_line(&F, NULL))
    {
      if (popinfile()) gp_quit();
      if (ismain) continue;
      pop_buffer(); return z;
    }
    if (check_meta(b->buf)) continue;

    if (ismain)
    {
#if defined(_WIN32) || defined(__CYGWIN32__)
      win32ctrlc = 0;
#endif
      TIMERstart(GP_DATA->T);
    }
    avma = top - av;
    pari_set_last_newline(1);
    z = gpreadseq(b->buf, GP_DATA->flags & STRICTMATCH);
    if (! ismain) continue;

    if (!pari_last_was_newline()) pariputc('\n');

    if (GP_DATA->flags & CHRONO)
      pariputs(gp_format_time(ti_REGULAR));
    else
      (void)gp_format_time(ti_NOPRINT);
    if (z == gnil) continue;

    if (GP_DATA->flags & SIMPLIFY) z = simplify_i(z);
    z = set_hist_entry(H, z);
    if (! is_silent(b->buf) ) gp_output(z, GP_DATA);
  }
}

/********************************************************************/
/*                                                                  */
/*                      EXCEPTION HANDLER                           */
/*                                                                  */
/********************************************************************/
void
gp_sigint_fun(void) { pari_err(siginter, gp_format_time(ti_INTERRUPT)); }

static void
gp_handle_SIGINT(void)
{
#if defined(_WIN32) || defined(__CYGWIN32__)
  win32ctrlc++;
#else
  if (GP_DATA->flags & TEXMACS) tm_start_output();
  gp_sigint_fun();
#endif
}

static void
gp_sighandler(int sig)
{
  char *msg;
#ifndef HAS_SIGACTION
  /*SYSV reset the signal handler in the handler*/
  (void)os_signal(sig,gp_sighandler);
#endif
  switch(sig)
  {
#ifdef SIGBREAK
    case SIGBREAK: gp_handle_SIGINT(); return;
#endif
#ifdef SIGINT
    case SIGINT:   gp_handle_SIGINT(); return;
#endif

#ifdef SIGSEGV
    case SIGSEGV: msg = "GP (Segmentation Fault)"; break;
#endif
#ifdef SIGBUS
    case SIGBUS:  msg = "GP (Bus Error)"; break;
#endif
#ifdef SIGFPE
    case SIGFPE:  msg = "GP (Floating Point Exception)"; break;
#endif

#ifdef SIGPIPE
    case SIGPIPE:
    {
      pariFILE *f = GP_DATA->pp->file;
      if (f && pari_outfile == f->file)
      {
        pari_err(talker, "Broken Pipe, resetting file stack...");
        GP_DATA->pp->file = NULL; /* to avoid oo recursion on error */
        pari_outfile = stdout; pari_fclose(f);
      }
      /*Do not attempt to write to stdout in case it triggered the SIGPIPE*/
      return; /* not reached */
    }
#endif
    default: msg = "signal handling"; break;
  }
  pari_err(bugparier, msg);
}

int
break_loop(long numerr)
{
  static FILE *oldinfile = NULL;
  static char *old = NULL;
  static Buffer *b = NULL;
  VOLATILE int go_on = 0;
  char *s, *t;
  filtre_t F;

  if (b) jump_to_given_buffer(b);
  b = new_buffer();
  push_stack(&bufstack, (void*)b);

  (void)&s; /* emulate volatile */
  old = s = get_analyseur();
  t = NULL;
  if (bufstack->prev)
  {
    Buffer *oldb = (Buffer*)bufstack->prev->value;
    t = oldb->buf;
    /* something fishy, probably a ^C, or we overran analyseur */
    if (!s || !s[-1] || s < t || s >= t + oldb->len) s = NULL;
  }
  oldinfile = infile;
  init_filtre(&F, b);

  term_color(c_ERR); pariputc('\n');
  errcontext("Break loop (type 'break' or Control-d to go back to GP)", s, t);
  if (s) pariputc('\n');
  term_color(c_NONE);
  if (numerr == siginter)
    pariputs("[type <Return> in empty line to continue]\n");
  infile = stdin;
  for(;;)
  {
    GEN x;
    if (setjmp(b->env)) pariputc('\n');
    if (! gp_read_line(&F, BREAK_LOOP_PROMPT))
    {
      if (popinfile()) break;
      continue;
    }
#if defined(_WIN32) || defined(__CYGWIN32__)
    win32ctrlc = 0;
#endif
    if (check_meta(b->buf))
    { /* break loop initiated by ^C? Empty input --> continue computation */
      if (numerr == siginter && *(b->buf) == 0) { go_on=1; break; }
      continue;
    }
    x = readseq(b->buf);
    if (did_break()) break;
    if (x == gnil || is_silent(b->buf)) continue;

    term_color(c_OUTPUT); gen_output(x, GP_DATA->fmt);
    term_color(c_NONE); pariputc('\n');
  }
  if (old && !s) set_analyseur(old);
  b = NULL; infile = oldinfile;
  pop_buffer(); return go_on;
}

int
gp_exception_handler(long numerr)
{
  char *s = (char*)global_err_data;
  if (!s) return 0;
  if (*s) {
    /* prevent infinite recursion in case s raises an exception */
    static int recovering = 0;
    if (recovering)
      recovering = 0;
    else
    {
      recovering = 1;
      fprintferr("\n"); outerr(readseq(s));
      recovering = 0; return 0;
    }
  }
  if (numerr == errpile) { var_make_safe(); avma = top; }
  return break_loop(numerr);
}

/********************************************************************/
/*                                                                  */
/*                      GP-SPECIFIC ROUTINES                        */
/*                                                                  */
/********************************************************************/
static void
check_secure(char *s)
{
  if (GP_DATA->flags & SECURE)
    pari_err(talker, "[secure mode]: system commands not allowed\nTried to run '%s'",s);
}

GEN
read0(char *s)
{
  switchin(s);
  if (file_is_binary(infile)) {
    int junk;
    return gpreadbin(s, &junk);
  }
  return gp_main_loop(0);
}

GEN
extern0(char *s)
{
  check_secure(s);
  infile = try_pipe(s, mf_IN)->file;
  return gp_main_loop(0);
}

GEN
input0(void)
{
  Buffer *b = new_buffer();
  filtre_t F;
  GEN x;

  init_filtre(&F, b);
  push_stack(&bufstack, (void*)b);
  while (! get_line_from_file(DFT_INPROMPT,&F,infile))
    if (popinfile()) { fprintferr("no input ???"); gp_quit(); }
  x = readseq(b->buf);
  pop_buffer(); return x;
}

void
system0(char *s)
{
#if defined(UNIX) || defined(__EMX__) || defined(_WIN32)
  check_secure(s); system(s);
#else
  pari_err(archer);
#endif
}

/*******************************************************************/
/**                                                               **/
/**                        INITIALIZATION                         **/
/**                                                               **/
/*******************************************************************/
static void
testuint(char *s, ulong *d) { if (s) *d = get_uint(s); }

static char *
read_arg(long *nread, char *t, long argc, char **argv)
{
  long i = *nread;
  if (isdigit((int)*t)) return t;
  if (*t || i==argc) usage(argv[0]);
  *nread = i+1; return argv[i];
}

static char *
read_arg_equal(long *nread, char *t, long argc, char **argv)
{
  long i = *nread;
  if (*t=='=' && isdigit((int)t[1])) return t+1;
  if (*t || i==argc) usage(argv[0]);
  *nread = i+1; return argv[i];
}

static void
init_trivial_stack()
{
  const size_t s = 2048;
  bot = (pari_sp)gpmalloc(s);
  avma = top = bot + s;
}

static void
read_opt(growarray A, long argc, char **argv)
{
  char *b = NULL, *p = NULL, *s = NULL;
  long i = 1, initrc = 1;

  (void)&p; (void)&b; (void)&s; /* -Wall gcc-2.95 */

  pari_outfile = stderr;
  while (i < argc)
  {
    char *t = argv[i];

    if (*t++ != '-') break;
    i++;
    switch(*t++)
    {
      case 'b': b = read_arg(&i,t,argc,argv);
        pari_warn(warner, "buffersize is no longer used. -b ignored");
        break;
      case 'p': p = read_arg(&i,t,argc,argv); break;
      case 's': s = read_arg(&i,t,argc,argv); break;

      case 'e':
	if (strncmp(t,"macs",4)) usage(argv[0]); /* obsolete */
        GP_DATA->flags |= EMACS; break;
      case 'q':
        GP_DATA->flags |= QUIET; break;
      case 't':
	if (strncmp(t,"est",3)) usage(argv[0]); /* obsolete */
        GP_DATA->flags |= TEST; /* fall through */
      case 'f':
	initrc = 0; break;
      case '-':
        if (strcmp(t, "version-short") == 0) { print_shortversion(); exit(0); }
        if (strcmp(t, "version") == 0) {
          init_trivial_stack(); print_version();
          free((void*)bot); exit(0);
        }
        if (strcmp(t, "texmacs") == 0) { GP_DATA->flags |= TEXMACS; break; }
        if (strcmp(t, "emacs") == 0) { GP_DATA->flags |= EMACS; break; }
        if (strcmp(t, "test") == 0) { GP_DATA->flags |= TEST; initrc = 0; break; }
        if (strcmp(t, "quiet") == 0) { GP_DATA->flags |= QUIET; break; }
        if (strcmp(t, "fast") == 0) { initrc = 0; break; }
        if (strncmp(t, "primelimit",10) == 0) {p = read_arg_equal(&i,t+10,argc,argv); break; }
        if (strncmp(t, "stacksize",9) == 0) {s = read_arg_equal(&i,t+9,argc,argv); break; }
       /* fall through */
      default:
	usage(argv[0]);
    }
  }
  if (GP_DATA->flags & TEXMACS) tm_start_output();
  if (GP_DATA->flags & TEST) init80col(0);
  if (initrc)
  {
    gp_initrc(A, argv[0]);
    if (setjmp(GP_DATA->env))
    {
      pariputs("### Errors on startup, exiting...\n\n");
      exit(1);
    }
  }
  for ( ; i < argc; i++) grow_append(A, pari_strdup(argv[i]));

  /* override the values from gprc */
  testuint(p, &(GP_DATA->primelimit));
  testuint(s, (ulong*)&top);
  if (GP_DATA->flags & (EMACS|TEXMACS|TEST)) disable_color = 1;
  pari_outfile = stdout;
}

#ifdef WINCE
int
WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance,
        LPWSTR lpCmdLine, int nShowCmd)
{
  char **argv = NULL;
  int argc = 1;
#else
int
main(int argc, char **argv)
{
#endif
  growarray A, *newfun, *oldfun;
  long i;

  GP_DATA = default_gp_data();
  initout(1);

  for (i=0; i<c_LAST; i++) gp_colors[i] = c_NONE;
  bot = (pari_sp)0;
  top = (pari_sp)(1000000*sizeof(long));

  if (setjmp(GP_DATA->env))
  {
    pariputs("### Errors on startup, exiting...\n\n");
    exit(1);
  }
#ifdef __MWERKS__
  argc = ccommand(&argv);
#endif
  grow_init(A);
  pari_init_defaults();
  read_opt(A, argc,argv);

  pari_init_opts(top-bot, GP_DATA->primelimit, 0);
  newfun = pari_get_modules();
  grow_append(*newfun, functions_gp);
  grow_append(*newfun, functions_highlevel);
  oldfun = pari_get_oldmodules();
  grow_append(*oldfun, functions_oldgp);
  if (new_fun_set)
  {
    pari_add_module(functions_gp);
    pari_add_module(functions_highlevel);
  }
  else
    pari_add_module(functions_oldgp);

  pari_sig_init(gp_sighandler);

  init_graph();
#ifdef READLINE
  init_readline();
#endif
  whatnow_fun = whatnow;
  sigint_fun = gp_sigint_fun;
  default_exception_handler = gp_exception_handler;
  gp_expand_path(GP_DATA->path);

  if (!(GP_DATA->flags & QUIET)) gp_head();
  if (A->n)
  {
    VOLATILE ulong f = GP_DATA->flags;
    FILE *l = logfile;
    VOLATILE long i;
    GP_DATA->flags &= ~(CHRONO|ECHO); logfile = NULL;
    for (i = 0; i < A->n; i++) {
      if (setjmp(GP_DATA->env))
      {
        fprintferr("... skipping file '%s'\n", A->v[i]);
        i++; if (i == A->n) break;
      }
      (void)read0((char*)A->v[i]); free(A->v[i]);
    }
    GP_DATA->flags = f; logfile = l;
  }
  grow_kill(A);
  (void)gp_main_loop(1);
  gp_quit(); return 0; /* not reached */
}

/*******************************************************************/
/**                                                               **/
/**                          GP OUTPUT                            **/
/**                                                               **/
/*******************************************************************/
    /* EXTERNAL PRETTYPRINTER */
/* Wait for prettinprinter to finish, to prevent new prompt from overwriting
 * the output.  Fill the output buffer, wait until it is read.
 * Better than sleep(2): give possibility to print */
static void
prettyp_wait(void)
{
  char *s = "                                                     \n";
  long i = 2000;

  pariputs("\n\n"); pariflush(); /* start translation */
  while (--i) pariputs(s);
  pariputs("\n"); pariflush();
}

/* initialise external prettyprinter (tex2mail) */
static int
prettyp_init(void)
{
  gp_pp *pp = GP_DATA->pp;
  if (!pp->cmd) return 0;
  if (pp->file || (pp->file = try_pipe(pp->cmd, mf_OUT))) return 1;

  pari_warn(warner,"broken prettyprinter: '%s'",pp->cmd);
  free(pp->cmd); pp->cmd = NULL; return 0;
}

/* n = history number. if n = 0 no history */
static int
tex2mail_output(GEN z, long n)
{
  pariout_t T = *(GP_DATA->fmt); /* copy */
  FILE *o_out, *o_logfile = logfile;

  if (!prettyp_init()) return 0;
  o_out = pari_outfile; /* save state */

  /* Emit first: there may be lines before the prompt */
  if (n) term_color(c_OUTPUT);
  pariflush();
  pari_outfile = GP_DATA->pp->file->file;
  T.prettyp = f_TEX;
  logfile = NULL;

  /* history number */
  if (n)
  {
    char s[128], c_hist[16], c_out[16];

    strcpy(c_hist, term_get_color(c_HIST));
    strcpy(c_out , term_get_color(c_OUTPUT));
    if (*c_hist || *c_out)
      sprintf(s, "\\LITERALnoLENGTH{%s}\\%%%ld =\\LITERALnoLENGTH{%s} ",
              c_hist, n, c_out);
    else
      sprintf(s, "\\%%%ld = ", n);
    pariputs_opt(s);
    if (o_logfile) {
      switch (logstyle) {
      case logstyle_plain:
        fprintf(o_logfile, "%%%ld = ", n);
        break;
      case logstyle_color:
        fprintf(o_logfile, "%s%%%ld = %s", c_hist, n, c_out);
        break;
      case logstyle_TeX:
        fprintf(o_logfile, "\\PARIout{%ld}", n);
        break;
      }
    }
  }
  /* output */
  gen_output(z, &T);

  /* flush and restore */
  prettyp_wait();
  if (o_logfile) {
    pari_outfile = o_logfile;
    /* XXXX Maybe it is better to output in another format? */
    if (logstyle == logstyle_TeX) {
      T.TeXstyle |= TEXSTYLE_BREAK;
      gen_output(z, &T);
      pariputc('%');
    } else
	outbrute(z);
    pariputc('\n'); pariflush();
  }
  logfile = o_logfile;
  pari_outfile = o_out;
  if (n) term_color(c_NONE);
  return 1;
}

    /* TEXMACS */
static void
texmacs_output(GEN z, long n)
{
  pariout_t T = *(GP_DATA->fmt); /* copy */
  char *sz;

  T.prettyp = f_TEX;
  T.fieldw = 0;
  sz = GENtostr0(z, &T, &gen_output);
  printf("%clatex:", DATA_BEGIN);
  if (n) printf("\\magenta\\%%%ld = ", n);
  printf("$\\blue %s$%c", sz,DATA_END);
  free(sz); fflush(stdout);
}

    /* REGULAR */
static void
normal_output(GEN z, long n)
{
  long l = 0;
  /* history number */
  if (n)
  {
    char s[64];
    term_color(c_HIST);
    sprintf(s, "%%%ld = ", n);
    pariputs_opt(s);
    l = strlen(s);
  }
  /* output */
  term_color(c_OUTPUT);
  if (GP_DATA->lim_lines)
    lim_lines_output(z, GP_DATA->fmt, l, GP_DATA->lim_lines);
  else
    gen_output(z, GP_DATA->fmt);
  term_color(c_NONE); pariputc('\n');
}

void
gp_output(GEN z, gp_data *G)
{
  if (G->flags & TEST) {
    init80col(0);
    gen_output(z, G->fmt); pariputc('\n');
  }
  else if (G->flags & TEXMACS)
    texmacs_output(z, G->hist->total);
  else if (G->fmt->prettyp != f_PRETTY || !tex2mail_output(z, G->hist->total))
    normal_output(z, G->hist->total);
  pariflush();
}
