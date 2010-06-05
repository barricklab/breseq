/* $Id: anal.c 10414 2008-07-08 15:55:11Z kb $

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
/*                  SYNTACTICAL ANALYZER FOR GP                    */
/*                                                                 */
/*******************************************************************/
#include "pari.h"
#include "paripriv.h"
#include "anal.h"

typedef GEN (*PFGEN)(ANYARG);
typedef GEN (*F2GEN)(GEN,GEN);
typedef GEN (*F1GEN)(GEN);

/* for user functions */
typedef struct gp_args {
  long nloc, narg; /* # of local variables, # of format parameters */
  GEN *arg; /* values for these nloc + narg variables */
} gp_args;


char *gp_function_name=NULL;

static GEN    constante();
static GEN    expr();
static GEN    facteur();
static GEN    identifier();
static GEN    read_member(GEN x);
static GEN    seq();
static GEN    truc();
static ulong  number(int *pn, char **ps);
static void   skipconstante();
static void   skipexpr();
static void   skipfacteur();
static void   skipidentifier();
static void   skipseq();
static void   skipstring();
static void   skiptruc();
static entree *entry();
static entree *skipentry(void);

static entree *installep(void *f,char *name,long l,long v,long add,entree **table);
#define VAR_POLS_LONGS	7 /* 4 words for pol_x, 3 for pol_1 */
#define SIZEOF_VAR_POLS	(VAR_POLS_LONGS*sizeof(long))

/* last time we began parsing an object of specified type */
static struct
{
  char *identifier, *symbol, *raw, *member, *start;
} mark;

/* points to the part of the string that remains to be parsed */
static char *analyseur = NULL;

/* when non-0, we are checking the syntax of a new function body */
static long skipping_fun_def;

/* when non-NULL, points to the entree of a new user function (currently
 * being checked). Used by the compatibility engine in the following way:
 *   when user types in a function whose name has changed, it is understood
 *   as EpNEW; first syntax error (missing = after function definition
 *   usually) triggers err_new_fun() if check_new_fun is set. */
static entree *check_new_fun;
#define NOT_CREATED_YET ((entree *)0x1L)
#define initial_value(ep) ((ep)+1)

/* for control statements */
enum { br_NONE = 0, br_BREAK, br_NEXT, br_MULTINEXT, br_RETURN, br_ALLOCMEM };
static long br_status, br_count;
static GEN br_res = NULL;

/* Mnemonic codes parser:
 *
 * TEMPLATE is assumed to be ";"-separated list of items.  Each item
 * may have one of the following forms: id=value id==value id|value id&~value.
 * Each id consists of alphanum characters, dashes and underscores.
 * IDs are case-sensitive.

 * ARG consists of several IDs separated by punctuation (and optional
 * whitespace).  Each modifies the return value in a "natural" way: an
 * ID from id=value should be the first in the sequence and sets RETVAL to
 * VALUE (and cannot be negated), ID from id|value bit-ORs RETVAL with
 * VALUE (and bit-ANDs RETVAL with ~VALUE if negated), ID from
 * id&~value behaves as if it were noid|value, ID from
 * id==value behaves the same as id=value, but should come alone.

 * For items of the form id|value and id&~value negated forms are
 * allowed: either when arg looks like no[-_]id, or when id looks like
 * this, and arg is not-negated. */

enum { A_ACTION_ASSIGN, A_ACTION_SET, A_ACTION_UNSET };
enum { PARSEMNU_TEMPL_TERM_NL, PARSEMNU_ARG_WHITESP };
#define IS_ID(c)	(isalnum((int)c) || ((c) == '_') || ((c) == '-'))
#define ERR(reason)	STMT_START {	\
    if (failure && first) {		\
	*failure = reason; *failure_arg = NULL; return 0;		\
    } else pari_err(talker,reason); } STMT_END
#define ERR2(reason,s)	STMT_START {	\
    if (failure && first) {		\
	*failure = reason; *failure_arg = s; return 0;		\
    } else pari_err(talker,reason,s); } STMT_END

unsigned long
parse_option_string(char *arg, char *tmplate, long flag, char **failure, char **failure_arg)
{
    unsigned long retval = 0;
    char *etmplate = NULL;

    if (flag & PARSEMNU_TEMPL_TERM_NL)
	etmplate = strchr(tmplate, '\n');
    if (!etmplate)
	etmplate = tmplate + strlen(tmplate);

    if (failure)
	*failure = NULL;
    while (1) {
	long numarg;
	char *e, *id;
	char *negated;			/* action found with 'no'-ID */
	int negate;			/* Arg has 'no' prefix removed */
	unsigned long l, action = 0, first = 1, singleton = 0;
	char b[80], *buf, *inibuf;

	if (flag & PARSEMNU_ARG_WHITESP)
	    while (isspace((int)*arg)) arg++;
	if (!*arg)
	    break;
	e = arg;
	while (IS_ID(*e)) e++;
	/* Now the ID is whatever is between arg and e. */
	l = e - arg;
	if (l >= sizeof(b))
	    ERR("id too long in a stringified flag");
	if (!l)				/* Garbage after whitespace? */
	    ERR("a stringified flag does not start with an id");
	strncpy(b, arg, l);
	b[l] = 0;
	arg = e;
	e = inibuf = buf = b;
	while (('0' <= *e) && (*e <= '9'))
	    e++;
	if (*e == 0)
	    ERR("numeric id in a stringified flag");	
	negate = 0;
	negated = NULL;
      find:
	id = tmplate;
	while ((id = strstr(id, buf)) && id < etmplate) {
	    if (IS_ID(id[l])) {		/* We do not allow abbreviations yet */
		id += l;		/* False positive */
		continue;
	    }
	    if ((id >= tmplate + 2) && (IS_ID(id[-1]))) {
		char *s = id;

		if ( !negate && s >= tmplate+3
		     && ((id[-1] == '_') || (id[-1] == '-')) )
		    s--;
		/* Check whether we are preceeded by "no" */
		if ( negate		/* buf initially started with "no" */
		     || (s < tmplate+2) || (s[-1] != 'o') || (s[-2] != 'n')
		     || (s >= tmplate+3 && IS_ID(s[-3]))) {
		    id += l;		/* False positive */
		    continue;
		}
		/* Found noID in the template! */
		id += l;
		negated = id;
		continue;		/* Try to find without 'no'. */
	    }
	    /* Found as is */
	    id += l;
	    break;
	}
	if ( !id && !negated && !negate
	     && (l > 2) && buf[0] == 'n' && buf[1] == 'o' ) {
	    /* Try to find the flag without the prefix "no". */
	    buf += 2; l -= 2;
	    if ((buf[0] == '_') || (buf[0] == '-')) { buf++; l--; }
	    negate = 1;
	    if (buf[0])
		goto find;
	}
	if (!id && negated) {	/* Negated and AS_IS forms, prefer AS_IS */
	    id = negated;	/* Otherwise, use negated form */
	    negate = 1;
	}
	if (!id)
	    ERR2("Unrecognized id '%s' in a stringified flag", inibuf);
	if (singleton && !first)
	    ERR("Singleton id non-single in a stringified flag");
	if (id[0] == '=') {
	    if (negate)
		ERR("Cannot negate id=value in a stringified flag");
	    if (!first)
		ERR("Assign action should be first in a stringified flag");
	    action = A_ACTION_ASSIGN;
	    id++;
	    if (id[0] == '=') {
		singleton = 1;
		id++;
	    }
	} else if (id[0] == '^') {
	    if (id[1] != '~')
		pari_err(talker, "Unrecognized action in a template");
	    id += 2;
	    if (negate)
		action = A_ACTION_SET;
	    else
		action = A_ACTION_UNSET;
	} else if (id[0] == '|') {
	    id++;
	    if (negate)
		action = A_ACTION_UNSET;
	    else
		action = A_ACTION_SET;
	}

	e = id;

	while ((*e >= '0' && *e <= '9')) e++;
	while (isspace((int)*e))
	    e++;
	if (*e && (*e != ';') && (*e != ','))
	    pari_err(talker, "Non-numeric argument of an action in a template");
	numarg = atol(id);		/* Now it is safe to get it... */
	switch (action) {
	case A_ACTION_SET:
	    retval |= numarg;
	    break;
	case A_ACTION_UNSET:
	    retval &= ~numarg;
	    break;
	case A_ACTION_ASSIGN:
	    retval = numarg;
	    break;
	default:
	    ERR("error in parse_option_string");
	}
	first = 0;
	if (flag & PARSEMNU_ARG_WHITESP)
	    while (isspace((int)*arg))
		arg++;
	if (*arg && !(ispunct((int)*arg) && *arg != '-'))
	    ERR("Junk after an id in a stringified flag");
	/* Skip punctuation */
	if (*arg)
	    arg++;
    }
    return retval;
}

/*  Special characters:
 *     ' ', '\t', '\n', '\\' are forbidden internally (suppressed by filtre).
 *     { } are forbidden everywhere and will be used to denote optional
 *     lexemes in the sequel.
 *
 *  Definitions: The sequence
 *    { a }* means any number (possibly 0) of object a.
 *
 *  seq: only this one can be empty.
 *    expr { [:;] expr }* { [:;] }
 *
 *  expr:
 *     expression = sequence of "facteurs" separated by binary operators
 *     whose priority are:
 *      1: *, /, \, \/, %, >>, <<                (highest)
 *      2: +, -
 *      3: <, <=, >, >=, !=, ==, <>
 *      4: &, &&, |, ||                  (lowest)
 *     read from left to right.
 *
 *  facteur:
 *    { [+-] } followed by a "truc", then by any succession of the
 *  following:
 *
 *        ~, ', !
 *  or    ^ facteur
 *  or    matrix_index { matrix_index }*
 *  or    .entry
 *
 *  truc:
 *      ! facteur
 *  or  # facteur
 *  or  ' entry
 *  or  identifier
 *  or  constante
 *  or  string {string}*
 *  or  matrix
 *  or  ( expr )
 *  or  % { ` }*  or %number
 *
 *  identifier:
 *    entry  followed by optional
 *
 *      matrix_assignment_block
 *  or  .entry { = seq }
 *  or  {'} ( arg_list )
 *  or  ( arg_list ) = seq
 *
 *  arg_list
 *    { arg } { , arg }*
 *
 *  arg:
 *    expr  or  &entry
 *    Note: &entry (pointer) not yet implemented for user functions
 *
 *  matrix
 *      [ A { ; A}* ] where A = { expr } { , { expr } }*
 *      All A must share the same length.
 *
 *  matrix_index:
 *      [ expr {,} ]
 *   or [ { expr } , expr ]
 *
 *  matrix_assignment_block:
 *     { matrix_index }  followed by
 *        = expr
 *     or ++  or --
 *     or op= expr  where op is one of the operators in expr 1: and 2:
 *
 *  entry:
 *      [A-Za-z][A-Za-z0-9_]*
 *
 *  string:
 *      " any succession of characters [^\]"
 *
 *  constante:
 *      number { . [0-9]* }  { expo }
 *   or .{number} { expo }
 *
 *  expo:
 *      [eE] {[+-]} { number }
 *
 *  number:
 *      [0-9]+
 */
char*
get_analyseur(void) { return analyseur; }

void
set_analyseur(char *s) { analyseur = s; }

INLINE void
seq_init(char *t)
{
  check_new_fun = NULL;
  skipping_fun_def = 0;
  mark.start = analyseur = t;
  br_status = br_NONE;
  if (br_res) { killbloc(br_res); br_res = NULL; }
}

#define HANDLE_FOREIGN(t)\
  if (foreignExprHandler && *t == foreignExprSwitch)\
    return (*foreignExprHandler)(t);

/* Do not modify (analyseur,mark.start) */
static GEN
readseq0(char *t, GEN (*f)(void))
{
  pari_sp av = top - avma;
  char *olds = analyseur, *olde = mark.start;
  GEN z;

  HANDLE_FOREIGN(t);

  seq_init(t); z = f();
  analyseur = olds; mark.start = olde;
  av = top - av; /* safer than recording av = avma: f() may call allocatemem */
  if (br_status)
  {
    if (br_res) return gerepilecopy(av, br_res);
    if (!z) { avma = av; return gnil; }
  }
  /* ep->value, beware: it may be killed anytime.  */
  if (isclone(z)) { avma = av; return gcopy(z); }
  return gerepileupto(av, z);
}
static GEN
readseq0_nobreak(char *t, GEN (*f)(void))
{
  pari_sp av = top - avma;
  char *olds = analyseur, *olde = mark.start;
  GEN z;

  HANDLE_FOREIGN(t);

  seq_init(t); z = f();
  analyseur = olds; mark.start = olde;
  if (br_status) pari_err(talker,"break not allowed");
  av = top - av; /* safer than recording av = avma: f() may call allocatemem */
  /* ep->value, beware: it may be killed anytime.  */
  if (isclone(z)) { avma = av; return gcopy(z); }
  return gerepileupto(av, z);
}

/* for sumiter: (void)readseq(t) */
void
readseq_void(char *t)
{
  const pari_sp av = top - avma;
  char *olds = analyseur, *olde = mark.start;

  if (foreignExprHandler && *t == foreignExprSwitch)
  { (void)(*foreignExprHandler)(t); return; }

  seq_init(t); (void)seq();
  analyseur = olds; mark.start = olde;
  /* safer than recording av = avma: f() may call allocatemem */
  avma = top - av;
}

GEN readseq_nobreak(char *t)  { return readseq0_nobreak(t, seq);  }
GEN readexpr_nobreak(char *t) { return readseq0_nobreak(t, expr); }
GEN readseq(char *t)  { return readseq0(t, seq);  }
GEN readexpr(char *t) { return readseq0(t, expr); }
/* filtered readseq = remove blanks and comments */
GEN
gp_read_str(char *s)
{
  char *t = filtre(s, (compatible == OLDALL));
  GEN x = readseq0(t, seq);
  free(t); return x;
}

static void
unused_chars(char *c, int strict)
{
  long n = 2 * term_width() - (17+19+1); /* Warning + unused... + . */
  char *s;
  if (strict) pari_err(talker2,"unused characters", analyseur, c);
  if ((long)strlen(analyseur) > n)
  {
    s = gpmalloc(n + 1);
    n -= 5;
    (void)strncpy(s,analyseur, n);
    strcpy(s + n, "[+++]");
  }
  else
    s = pari_strdup(analyseur);
  pari_warn(warner, "unused characters: %s", s);
  free(s);
}

/* check syntax, then execute */
GEN
gpreadseq(char *c, int strict)
{
  char *olds = analyseur, *olde = mark.start;
  GEN z;

  seq_init(c); skipseq(); if (*analyseur) unused_chars(c, strict);
  seq_init(c); z = seq();
  analyseur = olds; mark.start = olde;
  if (br_status)
  {
    if (br_res) return br_res;
    if (!z) return gnil;
  }
  return z;
}

static void
check_proto(char *code)
{
  char *s = code, *old;
  if (*s == 'l' || *s == 'v' || *s == 'i') s++;
  while (*s && *s != '\n') switch (*s++)
  {
    case '&':
    case '=':
    case 'E':
    case 'G':
    case 'I':
    case 'L':
    case 'M':
    case 'P':
    case 'S':
    case 'V':
    case 'f':
    case 'n':
    case 'p':
    case 'r':
    case 'x': break;
    case 's':
      if (*s == '*') s++;
      break;
    case 'D':
      if (*s == 'G' || *s == '&' || *s == 'n' || *s == 'I' || *s == 'V')
      {
        s++; break;
      }
      old = s; while (*s != ',') s++;
      if (*s != ',') pari_err(talker2, "missing comma", old, code);
      break;
    case ',': break;
    case '\n': return; /* Before the mnemonic */

    case 'l':
    case 'i':
    case 'v': pari_err(talker2, "this code has to come first", s-1, code);
    default: pari_err(talker2, "unknown parser code", s-1, code);
  }
}

entree *
install(void *f, char *name, char *code)
{
  long hash;
  entree *ep = is_entry_intern(name, functions_hash, &hash);

  check_proto(code);
  if (ep)
  {
    if (ep->valence != EpINSTALL)
      pari_err(talker,"[install] identifier '%s' already in use", name);
    pari_warn(warner, "[install] updating '%s' prototype; module not reloaded", name);
    if (ep->code) free(ep->code);
  }
  else
  {
    char *s = name;
    if (isalpha((int)*s))
      while (is_keyword_char(*++s)) /* empty */;
    if (*s) pari_err(talker2,"not a valid identifier", s, name);
    ep = installep(f, name, strlen(name), EpINSTALL, 0, functions_hash + hash);
  }
  ep->code = pari_strdup(code);
  return ep;
}

/*******************************************************************/
/*                                                                 */
/*                            VARIABLES                            */
/*                                                                 */
/*******************************************************************/
/* As a rule, ep->value is a clone (COPY). push_val and pop_val are private
 * functions for use in sumiter: we want a temporary ep->value, which is NOT
 * a clone (PUSH), to avoid unnecessary copies. */

/* ep->args is the stack of old values (INITIAL if initial value, from
 * installep) */
typedef struct var_cell {
  struct var_cell *prev; /* cell associated to previous value on stack */
  GEN value; /* last value (not including current one, in ep->value) */
  char flag; /* status of _current_ ep->value: PUSH or COPY ? */
} var_cell;
#define INITIAL NULL
enum {PUSH_VAL = 0, COPY_VAL = 1};
#define killvalue(v) pop_val_full(get_ep(v))
/* Push x on value stack associated to ep. Assume EpVALENCE(ep)=EpVAR/EpGVAR */
static void
new_val_cell(entree *ep, GEN x, char flag)
{
  var_cell *v = (var_cell*) gpmalloc(sizeof(var_cell));
  v->value  = (GEN)ep->value;
  v->prev   = (var_cell*) ep->args;
  v->flag   = flag;

  /* beware: f(p) = Nv = 0
   *         Nv = p; f(Nv) --> this call would destroy p [ isclone ] */
  ep->value = (flag == COPY_VAL)? gclone(x):
                                  (x && isclone(x))? gcopy(x): x;
  /* Do this last. In case the clone is <C-C>'ed before completion ! */
  ep->args  = (void*)v;
}

void
free_ep_args(entree *ep)
{
  long i;
  gp_args *f = (gp_args*)ep->args;
  GEN *y = f->arg;
  for (i = f->narg + f->nloc - 1; i>=0; i--)
    if (isclone(y[i])) gunclone(y[i]);
  ep->args = INITIAL;
}

void
freeep(entree *ep)
{
  if (foreignFuncFree && ep->code && (*ep->code == 'x'))
    (*foreignFuncFree)(ep); /* function created by foreign interpreter */

  if (EpSTATIC(ep)) return; /* gp function loaded at init time */
  if (ep->help) free(ep->help);
  if (ep->code) free(ep->code);
  switch(EpVALENCE(ep))
  {
    case EpVAR:
    case EpGVAR:
      while (ep->args) pop_val(ep);
      break;
    case EpUSER:
      free_ep_args(ep); /* fall through */
    case EpALIAS:
      gunclone((GEN)ep->value); break;
  }
  free(ep);
}

static entree*
get_ep(long v)
{
  entree *ep = varentries[v];
  if (!ep) pari_err(talker2,"this function uses a killed variable",
               mark.identifier, mark.start);
  return ep;
}

INLINE void
copyvalue(long v, GEN x) {
  new_val_cell(get_ep(v), x, typ(x) >= t_VEC ? COPY_VAL: PUSH_VAL);
}
INLINE void
pushvalue(long v, GEN x) { new_val_cell(get_ep(v), x, PUSH_VAL); }

void
push_val(entree *ep, GEN a) { new_val_cell(ep, a, PUSH_VAL); }

/* kill ep->value and replace by preceding one, poped from value stack */
void
pop_val_full(entree *ep)
{
  var_cell *v = (var_cell*) ep->args;

  if (v == INITIAL) return;
  killbloc((GEN)ep->value); /* inspect components to avoid memory leaks */
  ep->value = v->value;
  ep->args  = (void*) v->prev;
  free((void*)v);
}
void
pop_val(entree *ep)
{
  var_cell *v = (var_cell*) ep->args;

  if (v == INITIAL) return;
  /* otherwise often not a valid GEN due to GC */
  if (v->flag == COPY_VAL) killbloc((GEN)ep->value);
  ep->value = v->value;
  ep->args  = (void*) v->prev;
  free((void*)v);
}

/* as above IF ep->value was PUSHed, or was created after block number 'loc'
   return 0 if not deleted, 1 otherwise [for recover()] */
int
pop_val_if_newer(entree *ep, long loc)
{
  var_cell *v = (var_cell*) ep->args;

  if (v == INITIAL) return 0;
  if (v->flag == COPY_VAL && !pop_entree_bloc(ep, loc)) return 0;
  ep->value = v->value;
  ep->args  = (void*) v->prev;
  free((void*)v); return 1;
}

/* set new value of ep directly to val (COPY), do not save last value unless
 * it's INITIAL. */
void
changevalue(entree *ep, GEN x)
{
  var_cell *v = (var_cell*) ep->args;
  if (v == INITIAL) new_val_cell(ep, x, COPY_VAL);
  else
  {
    x = gclone(x); /* beware: killbloc may destroy old x */
    if (v->flag == COPY_VAL) killbloc((GEN)ep->value); else v->flag = COPY_VAL;
    ep->value = (void*)x;
  }
}

/* as above, but PUSH, notCOPY */
void
changevalue_p(entree *ep, GEN x)
{
  var_cell *v = (var_cell*) ep->args;
  if (v == INITIAL) new_val_cell(ep,x, PUSH_VAL);
  else
  {
    if (v->flag == COPY_VAL) { killbloc((GEN)ep->value); v->flag = PUSH_VAL; }
    ep->value = (void*)x;
  }
}

/* make GP variables safe for avma = top */
void
var_make_safe()
{
  long n;
  entree *ep;
  for (n = 0; n < functions_tblsz; n++)
    for (ep = functions_hash[n]; ep; ep = ep->next)
      if (EpVALENCE(ep) == EpVAR)
      { /* make sure ep->value is a COPY */
        var_cell *v = (var_cell*)ep->args;
        if (v && v->flag == PUSH_VAL) {
          GEN x = (GEN)ep->value;
          if (x) changevalue(ep, (GEN)ep->value); else pop_val(ep);
        }
      }
}

/* n = hashvalue(ep->name), delete ep */
void
kill_from_hashlist(entree *ep, long n)
{
  entree *e = functions_hash[n];
  if (e == ep) { functions_hash[n] = ep->next; return; }
  for (; e; e = e->next)
    if (e->next == ep) { e->next = ep->next; return; }
}

/* kill aliases pointing to EP */
static void 
kill_alias(entree *EP)
{
  entree *ep, *epnext;
  long n;
  for (n = 0; n < functions_tblsz; n++)
    for (ep = functions_hash[n]; ep; ep = epnext)
    {
      epnext = ep->next;
      if (EpVALENCE(ep) == EpALIAS &&
          EP == (entree *) ((GEN)ep->value)[1]) kill0(ep);
    }
}

/* Kill entree ep, i.e free all memory it occupies, remove it from hashtable.
 * If it's a variable set a "black hole" in pol_x[v], etc. x = 0-th variable
 * can NOT be killed (only the value): we often use explicitly pol_x[0] */
void
kill0(entree *ep)
{
  char *s = ep->name;
  long v;

  if (EpSTATIC(ep))
    pari_err(talker2,"can't kill that",mark.symbol,mark.start);
  switch(EpVALENCE(ep))
  {
    case EpVAR:
    case EpGVAR:
      while (ep->args) pop_val(ep);
      v = varn(ep->value); if (!v) return; /* never kill x */

      gel(polvar,v+1) = pol_x[v] = pol_1[v] = gnil;
      varentries[v] = NULL; break;
    case EpUSER: kill_alias(ep);
      break;
  }
  kill_from_hashlist(ep, hashvalue(&s));
  freeep(ep);
}

void
addhelp(entree *ep, char *s)
{
  if (ep->help && !EpSTATIC(ep)) free(ep->help);
  ep->help = pari_strdup(s);
}

GEN
type0(GEN x)
{
  const char *s = type_name(typ(x));
  return strtoGENstr(s);
}

/*******************************************************************/
/*                                                                 */
/*                              PARSER                             */
/*                                                                 */
/*******************************************************************/
#define separator(c)  ((c)==';' || (compatible && (c)==':'))

static void
allocate_loop_err(void) {
  pari_err(talker2,"can't allow allocatemem() in loops", analyseur, mark.start);
}

static GEN
seq(void)
{
  const pari_sp av = top - avma;
  GEN res = gnil;
  int allocmem = 0;

  for(;;)
  {
    while (separator(*analyseur)) analyseur++;
    if (!*analyseur || *analyseur == ')' || *analyseur == ',') break;
    res = expr();
    if (br_status) {
      if (br_status != br_ALLOCMEM) break;
      br_status = br_NONE;
      allocmem = 1;
    }
    if (!separator(*analyseur)) break; else analyseur++;

    if (top - avma > ((top - av)>>1))
    {
      if(DEBUGMEM>1) pari_warn(warnmem,"seq");
      if (is_universal_constant(res)) avma = top - av;
      else
	res = gerepilecopy(top - av, res);
    }
  }
  if (allocmem) {
    if (br_status) allocate_loop_err();
    br_status = br_ALLOCMEM;
  }
  return res;
}

static GEN
gshift_l(GEN x, GEN n)  {
  if (is_bigint(n)) pari_err(talker2,"integer too big",analyseur,mark.start);
  return gshift(x, itos(n));
}
static GEN
gshift_r(GEN x, GEN n) {
  if (is_bigint(n)) pari_err(talker2,"integer too big",analyseur,mark.start);
  return gshift(x,-itos(n));
}

#define UNDEF (GEN)0x1
static GEN
expr(void)
{
  pari_sp av = top - avma;
  GEN aux,e,e1,e2,e3;
  F2GEN F1,F2,F3;
  int F0 = 0;
 
  F1 = F2 = F3 = (F2GEN)NULL;
  e1 = e2 = e3 = UNDEF;
L3:
#define act(fun) \
  aux = facteur(); if (br_status) return aux;\
  e3 = fun(e3,aux); goto L

  aux = facteur(); if (br_status) return aux;
  e3 = F3? F3(e3,aux): aux;
L:
  switch(*analyseur)
  {
    case '*': analyseur++; act(gmul);
    case '/': analyseur++; act(gdiv);
    case '%': analyseur++; act(gmod);
    case '\\':
      if (analyseur[1] != '/') { analyseur++; act(gdivent); }
      analyseur += 2; act(gdivround);

    case '<':
      if (analyseur[1] != '<') break;
      analyseur += 2;
      { char *old = analyseur; /* act(shift_l) + error checks */
        aux = facteur(); if (br_status) return aux;
        if (typ(aux) != t_INT) pari_err(talker2,"not an integer",old,mark.start);
        if (is_bigint(aux)) pari_err(talker2,"shift operand too big",old,mark.start);
        e3 = gshift(e3, itos(aux)); goto L;
      }
    case '>':
      if (analyseur[1] != '>') break;
      analyseur += 2;
      { char *old = analyseur; /* act(shift_r) + error checks */
        aux = facteur(); if (br_status) return aux;
        if (typ(aux) != t_INT) pari_err(talker2,"not an integer",old,mark.start);
        if (is_bigint(aux)) pari_err(talker2,"shift operand too big",old,mark.start);
        e3 = gshift(e3,-itos(aux)); goto L;
      }
  }
  F3 = (F2GEN)NULL;

L2:
  if (e3 == UNDEF) goto L3;
  e2 = F2? F2(e2,e3): e3;
  e3 = UNDEF;
  if (top - avma > ((top - bot)>>1))
  {
    if(DEBUGMEM>1) pari_warn(warnmem,"expr");
    gerepileall(top - av, (e1==UNDEF)?1: 2, &e2, &e1);
  }

  switch(*analyseur)
  {
    case '+': analyseur++; F2=&gadd; goto L3;
    case '-': analyseur++; F2=&gsub; goto L3;
  }
  F2 = (F2GEN)NULL;

L1:
  if (e2 == UNDEF) goto L2;
  e1 = F1? F1(e1,e2): e2;
  e2 = UNDEF;
  switch(*analyseur)
  {
    case '<':
      switch(*++analyseur)
      {
        case '=': analyseur++; F1=&gle; goto L2;
        case '>': analyseur++; F1=&gne; goto L2;
      }
      F1=&glt; goto L2;

    case '>':
      if (*++analyseur == '=') { analyseur++; F1=&gge; goto L2; }
      F1=&ggt; goto L2;

    case '=':
      if (analyseur[1] == '=') { analyseur+=2; F1=&geq; goto L2; }
      goto L1;

    case '!':
      if (analyseur[1] == '=') { analyseur+=2; F1=&gne; goto L2; }
      goto L1;
  }
  F1 = (F2GEN)NULL;

/* L0: */
  if (e1 == UNDEF) goto L1;
  e = F0? (gcmp0(e1)? gen_0: gen_1): e1;
  e1 = UNDEF;
  switch(*analyseur)
  {
    case '&':
      if (*++analyseur == '&') analyseur++;
      if (gcmp0(e)) { skipexpr(); return gen_0; }
      F0=1; goto L1;

    case '|':
      if (*++analyseur == '|') analyseur++;
      if (!gcmp0(e)) { skipexpr(); return gen_1; }
      F0=1; goto L1;
  }
  return e;
}
#undef UNDEF

/********************************************************************/
/**                                                                **/
/**                        CHECK FUNCTIONS                         **/
/**                                                                **/
/********************************************************************/
/* Should raise an error. If neighbouring identifier was a function in
 * 1.39.15, raise "obsolete" error instead. If check_new_fun doesn't help,
 * guess offending function was last identifier */
#define LEN 127
static void
err_new_fun()
{
  char s[LEN+1], *t;
  long n;

  if (check_new_fun == NOT_CREATED_YET) check_new_fun = NULL;
  t = check_new_fun? check_new_fun->name: mark.identifier;
  for (n=0; n < LEN; n++)
    if (!is_keyword_char(t[n])) break;
  (void)strncpy(s,t, n); s[n] = 0;
  if (check_new_fun) { kill0(check_new_fun); check_new_fun = NULL ; }
  if (compatible != NONE) return;

  if (whatnow_fun && (n = whatnow_fun(s,1)))
    pari_err(obsoler,mark.identifier,mark.start, s,n);
}
#undef LEN

static void
err_match(char *s, char c)
{
  char str[64];
  if (check_new_fun && (c == '(' || c == '=' || c == ',')) err_new_fun();
  sprintf(str,"expected character: '%c' instead of",c);
  pari_err(talker2,str,s,mark.start);
}

#define match2(s,c) if (*s != c) err_match(s,c);
#define match(c) \
  STMT_START { match2(analyseur, c); analyseur++; } STMT_END
#define NO_BREAK(s, old)\
  if (br_status) pari_err(talker2, "break not allowed "/**/s, (old), mark.start)

static long
readlong()
{
  const pari_sp av = avma;
  const char *old = analyseur;
  long m;
  GEN x = expr();

  NO_BREAK("here (reading long)", old);
  if (typ(x) != t_INT) pari_err(talker2,"this should be an integer", old,mark.start);
  if (is_bigint(x)) pari_err(talker2,"integer too big",old,mark.start);
  m = itos(x); avma=av; return m;
}

static long
check_array_index(long max)
{
  const char *old = analyseur;
  const long c = readlong();

  if (c < 1 || c >= max)
  {
    char s[80];
    sprintf(s,"array index (%ld) out of allowed range ",c);
    if (max == 1) strcat(s, "[none]");
    else if (max == 2) strcat(s, "[1]");
    else sprintf(s,"%s[1-%ld]",s,max-1);
    pari_err(talker2,s,old,mark.start);
  }
  return c;
}

static long
readvar()
{
  const char *old = analyseur;
  const GEN x = expr();

  if (typ(x) != t_POL || lg(x) != 4 ||
    !gcmp0(gel(x,2)) || !gcmp1(gel(x,3))) pari_err(varer1,old,mark.start);
  return varn(x);
}

/* noparen = 1 means function was called without (). Do we need to insert a
 * default argument ? */
static int
do_switch(int noparen, int matchcomma)
{
  const char *s = analyseur;
  if (noparen || *s == ')') return 1;
  if (*s == ',') /* we just read an arg, or first arg */
  {
    if (!matchcomma && s[-1] == '(') return 1; /* first arg */
    if (s[1] == ',' || s[1] == ')') { analyseur++; return 1; }
  }
  return 0;
}

/********************************************************************/
/**                                                                **/
/**                            STRINGS                             **/
/**                                                                **/
/********************************************************************/

static char*
init_buf(long len, char **ptbuf, char **ptlim)
{
  char *buf = (char *)new_chunk(2 + len / sizeof(long));
  *ptbuf = buf; *ptlim = buf + len; return buf;
}

static char*
realloc_buf(char *bp, long len, char **ptbuf,char **ptlimit)
{
  char *buf = *ptbuf;
  long newlen = ((*ptlimit - buf) + len) << 1;
  long oldlen = bp - buf;

  (void)init_buf(newlen, ptbuf, ptlimit);
  memcpy(*ptbuf, buf, oldlen);
  return *ptbuf + oldlen;
}

static char *
expand_string(char *bp, char **ptbuf, char **ptlimit)
{
  char *tmp = NULL; /* -Wall */
  long len = 0; /* -Wall */
  int alloc = 1;

  if (is_keyword_char(*analyseur))
  {
    char *s = analyseur;
    do s++; while (is_keyword_char(*s));

    if ((*s == '"' || *s == ',' || *s == ')') && !is_entry(analyseur))
    { /* Do not create new user variable. Consider as a literal */
      tmp = analyseur;
      len = s - analyseur;
      analyseur = s;
      alloc = 0;
    }
  }

  if (alloc)
  {
    pari_sp av = avma;
    char *old = analyseur;
    GEN z = expr();
    NO_BREAK("here (expanding string)", old);
    tmp = GENtostr(z);
    len = strlen(tmp); avma = av;
  }
  if (ptlimit && bp + len > *ptlimit)
    bp = realloc_buf(bp, len, ptbuf,ptlimit);
  memcpy(bp,tmp,len); /* ignore trailing \0 */
  if (alloc) free(tmp);
  return bp + len;
}

static char *
translate(char **src, char *s, char **ptbuf, char **ptlim)
{
  char *t = *src;
  while (*t)
  {
    while (*t == '\\')
    {
      switch(*++t)
      {
	case 'e':  *s='\033'; break; /* escape */
	case 'n':  *s='\n'; break;
	case 't':  *s='\t'; break;
	default:   *s=*t; if (!*t) pari_err(talker,"unfinished string");
      }
      t++; s++;
    }
    if (*t == '"')
    {
      if (t[1] != '"') break;
      t += 2; continue;
    }
    if (ptlim && s >= *ptlim)
      s = realloc_buf(s,1, ptbuf,ptlim);
    *s++ = *t++;
  }
  *s=0; *src=t; return s;
}

static char *
readstring_i(char *s, char **ptbuf, char **ptlim)
{
  match('"'); s = translate(&analyseur,s, ptbuf,ptlim); match('"');
  return s;
}

static GEN
any_string()
{
  long n = 1, len = 16;
  GEN res = cget1(len + 1, t_VEC);

  while (*analyseur)
  {
    if (*analyseur == ')' || *analyseur == ';') break;
    if (*analyseur == ',')
      analyseur++;
    else
    {
      char *old = analyseur;
      gel(res,n++) = expr();
      NO_BREAK("in print()", old);
    }
    if (n == len)
    {
      long newlen = len << 1;
      GEN p1 = cget1(newlen + 1, t_VEC);
      for (n = 1; n < len; n++) p1[n] = res[n];
      res = p1; len = newlen;
    }
  }
  setlg(res, n); return res;
}

/*  Read a "string" from src. Format then copy it, starting at s. Return
 *  pointer to char following the end of the input string */
char *
readstring(char *src, char *s)
{
  match2(src, '"'); src++; s = translate(&src, s, NULL,NULL);
  match2(src, '"'); return src+1;
}

/* return the first n0 chars of s as a GEN [s may not be 0-terminated] */
static GEN
_strtoGENstr(const char *s, long n0)
{
  long n = nchar2nlong(n0+1);
  GEN x = cgetg(n+1, t_STR);
  char *t = GSTR(x);
  strncpy(t, s, n0); t[n0] = 0; return x;
}

GEN
strtoGENstr(const char *s) { return _strtoGENstr(s, strlen(s)); }
/********************************************************************/
/**                                                                **/
/**                          READ FUNCTIONS                        **/
/**                                                                **/
/********************************************************************/
typedef struct matcomp
{
  GEN *ptcell;
  GEN parent;
  int full_col, full_row;
} matcomp;
typedef struct gp_pointer
{
  matcomp c;
  GEN x;
  entree *ep;
} gp_pointer;

/* Return the content of the matrix cell and sets members of corresponding
 * matrix component 'c'.  Assume *analyseur = '[' */
static GEN
matcell(GEN p, matcomp *C)
{
  GEN *pt = &p;
  long c, r;
  C->full_col = C->full_row = 0;
  do {
    analyseur++; p = *pt;
    switch(typ(p))
    {
      case t_VEC: case t_COL:
        c = check_array_index(lg(p));
        pt = (GEN*)(p + c); match(']'); break;

      case t_LIST:
        c = check_array_index(lgeflist(p)-1) + 1;
        pt = (GEN*)(p + c); match(']'); break;

      case t_VECSMALL:
        c = check_array_index(lg(p));
        pt = (GEN*)(p + c); match(']');
        if (*analyseur == '[') pari_err(caracer1,analyseur,mark.start);
        C->parent = p;
        C->ptcell = pt; return stoi((long)*pt);

      case t_MAT:
        if (lg(p)==1) pari_err(talker2,"a 0x0 matrix has no elements",
                                  analyseur,mark.start);
        C->full_col = C->full_row = 0;
        if (*analyseur==',') /* whole column */
        {
          analyseur++;
          c = check_array_index(lg(p));
          match(']');
          if (*analyseur == '[')
          { /* collapse [,c][r] into [r,c] */
            analyseur++;
            r = check_array_index(lg(p[c]));
            pt = (GEN*)((gel(p,c)) + r); /* &coeff(p,r,c) */
            match(']');
          }
          else
          {
            C->full_col = 1;
            pt = (GEN*)(p + c);
          }
          break;
        }

        r = check_array_index(lg(p[1]));
        match(',');
        if (*analyseur == ']') /* whole row */
        {
          analyseur++;
          if (*analyseur == '[')
          { /* collapse [r,][c] into [r,c] */
            analyseur++;
            c = check_array_index(lg(p));
            pt = (GEN*)((gel(p,c)) + r); /* &coeff(p,r,c) */
            match(']');
          }
          else
          {
            GEN p2 = cgetg(lg(p),t_VEC);
            C->full_row = r; /* record row number */
            for (c=1; c<lg(p); c++) p2[c] = coeff(p,r,c);
            pt = &p2;
          }
        }
        else
        {
          c = check_array_index(lg(p));
          pt = (GEN*)((gel(p,c)) + r); /* &coeff(p,r,c) */
          match(']');
        }
        break;

      default:
        pari_err(caracer1,analyseur-1,mark.start);
    }
  } while (*analyseur == '[');
  C->parent = p;
  C->ptcell = pt; return *pt;
}

static GEN
facteur(void)
{
  const char *old = analyseur;
  GEN x, p1;
  int plus;

  switch(*analyseur)
  {
    case '-': analyseur++; plus = 0; break;
    case '+': analyseur++; plus = 1; break;
    default: plus = 1; break;
  }
  x = truc(); if (br_status) return x;

  for(;;)
    switch(*analyseur)
    {
      case '.':
	analyseur++; x = read_member(x);
        if (!x) pari_err(talker2, "not a proper member definition",
                    mark.member, mark.start);
        break;
      case '^':
	analyseur++; p1 = facteur();
        NO_BREAK("after ^", old);
        x = gpow(x,p1,precreal); break;
      case '\'':
	analyseur++; x = deriv(x, -1); break;
      case '~':
	analyseur++; x = gtrans(x); break;
      case '[':
      {
        matcomp c;
        x = matcell(x, &c);
        if (isonstack(x)) x = gcopy(x);
        break;
      }
      case '!':
	if (analyseur[1] != '=')
	{
	  if (typ(x) != t_INT) pari_err(talker2,"this should be an integer",
                                           old,mark.start);
          if (is_bigint(x)) pari_err(talker2,"integer too big",old,mark.start);
	  analyseur++; x=mpfact(itos(x)); break;
	} /* Fall through */

      default:
        return (plus || x==gnil)? x: gneg(x);
    }
}

/* table array of length N+1, append one expr, growing array if necessary  */
static void
_append(GEN **table, long *n, long *N)
{
  char *old = analyseur;
  if (++(*n) == *N)
  {
    *N <<= 1;
    *table = (GEN*)gprealloc((void*)*table,(*N + 1)*sizeof(GEN));
  }
  (*table)[*n] = expr();
  NO_BREAK("in array context", old);
}

#define check_var_name() \
  if (!isalpha((int)*analyseur)) pari_err(varer1,analyseur,mark.start);

static GEN
truc(void)
{
  long N, i, j, m, n, p;
  GEN *table, z;
  char *old;
  pari_sp av0, av, lim;

  if (isalpha((int)*analyseur)) return identifier();
  if (isdigit((int)*analyseur) || *analyseur == '.') return constante();

  switch(*analyseur)
  {
    case '!': /* NOT */
      analyseur++; old = analyseur;
      z = facteur(); NO_BREAK("after !", old);
      return gcmp0(z)? gen_1: gen_0;

    case ('\''): { /* QUOTE */
      entree *ep;
      analyseur++; check_var_name();
      old = analyseur; ep = entry();
      switch(EpVALENCE(ep))
      {
        case EpVAR: case EpGVAR:
          return (GEN)initial_value(ep);
        default: pari_err(varer1,old,mark.start);
      }
    }
    case '#': /* cardinal */
      analyseur++; old = analyseur;
      z = facteur(); NO_BREAK("after #", old);
      return stoi(glength(z));

    case '"': /* string */
      analyseur++; old = analyseur;
      skipstring();
      n = nchar2nlong(analyseur - old); /* do not count enclosing '"' */
      z = cgetg(n+1, t_STR);
      (void)translate(&old, GSTR(z), NULL,NULL);
      return z;

    case '(':
      analyseur++;
      z = expr(); match(')'); return z;

    case '[': /* constant array/vector */
      analyseur++;
      if (*analyseur == ';' && analyseur[1] == ']')
	{ analyseur += 2; return cgetg(1,t_MAT); } /* [;] */

      n = 0; N = 1024;
      table = (GEN*) gpmalloc((N + 1)*sizeof(GEN));

      av0 = av = avma; lim = stack_lim(av,2);
      if (*analyseur != ']') _append(&table, &n, &N);
      while (*analyseur == ',') {
        analyseur++; _append(&table, &n, &N);
        if (low_stack(lim, stack_lim(av,2))) {
          if(DEBUGMEM>1) pari_warn(warnmem,"truc(): n = %ld", n);
          gerepilecoeffs(av0, (GEN)(table+1), n);
          av = avma; lim = stack_lim(av,2);
        }
      }
      switch (*analyseur++)
      {
	case ']':
        {
          long tx;
          if (*analyseur == '~') { analyseur++; tx=t_COL; } else tx=t_VEC;
	  z = cgetg(n+1,tx);
          if (n < 500) {
            for (i=1; i<=n; i++) gel(z,i) = gcopy(table[i]);
          } else {
            /* huge vector: GC needed */
            for (i=1; i<=n; i++) gel(z,i) = table[i];
            z = gerepilecopy(av0, z);
          }
	  break;
        }

	case ';':
	  m = n;
	  do _append(&table, &n, &N); while (*analyseur++ != ']');
	  z = cgetg(m+1,t_MAT); p = n/m + 1;
	  for (j=1; j<=m; j++)
	  {
            GEN c = cgetg(p,t_COL); gel(z,j) = c;
	    for (i=j; i<=n; i+=m) gel(++c,0) = gcopy(table[i]);
	  }
	  break;

	default: /* can only occur in library mode */
          pari_err(talker,"incorrect vector or matrix");
          return NULL; /* not reached */
      }
      free(table); return z;

    case '%': {
      gp_hist *H = GP_DATA->hist;
      int junk;

      old = analyseur; analyseur++;
      p = 0;
      if (*analyseur == '#') { analyseur++; return utoi(H->total); }
      while (*analyseur == '`') { analyseur++; p++; }
      if (p) return gp_history(H, -p, old, mark.start);
      if (!isdigit((int)*analyseur)) return gp_history(H, 0, old, mark.start);
      p = (long)number(&junk,&analyseur);
      if (!p) pari_err(talker2, "I can't remember before the big bang",
                  old, mark.start);
      return gp_history(H, p, old, mark.start);
    }
  }
  pari_err(caracer1,analyseur,mark.start);
  return NULL; /* not reached */
}

/* valid x opop, e.g x++ */
static GEN
double_op()
{
  char c = *analyseur;
  if (c && c == analyseur[1])
    switch(c)
    {
      case '+': analyseur+=2; return gen_1; /* ++ */
      case '-': analyseur+=2; return gen_m1; /* -- */
    }
  return NULL;
}

/* return op if op= detected */
static F2GEN
get_op_fun()
{
  char c = *analyseur, c1;
  if (c && (c1 = analyseur[1]))
  {
    if (c1 == '=')
    {
      switch(c)
      {
        case '+' : analyseur += 2; return &gadd;
        case '-' : analyseur += 2; return &gsub;
        case '*' : analyseur += 2; return &gmul;
        case '/' : analyseur += 2; return &gdiv;
        case '%' : analyseur += 2; return &gmod;
        case '\\': analyseur += 2; return &gdivent;
      }
    }
    else if (analyseur[2] == '=')
    {
      switch(c)
      {
        case '>' : if (c1=='>') { analyseur += 3; return &gshift_r; }
          break;
        case '<' : if (c1=='<') { analyseur += 3; return &gshift_l; }
          break;
        case '\\': if (c1=='/') { analyseur += 3; return &gdivround; }
          break;
      }
    }
  }
  return (F2GEN)NULL;
}

static GEN
expr_ass()
{
  char *old = analyseur;
  GEN res = expr();
  NO_BREAK("in assignment", old);
  return res;
}

static F2GEN
affect_block(GEN *res)
{
  F2GEN f;
  GEN r;
  if (*analyseur == '=')
  {
    r = NULL; f = NULL;
    if (analyseur[1] != '=') { analyseur++; r = expr_ass(); }
  }
  else if ((r = double_op()))  f = &gadd;
  else if ((f = get_op_fun())) r = expr_ass();
  *res = r; return f;
}

/* assign res at *pt in "simple array object" p and return it, or a copy.
 * In the latter case, reset avma to av */
static GEN
change_compo(pari_sp av, matcomp *c, GEN res)
{
  GEN p = c->parent, *pt = c->ptcell;
  long i;
  char *old = analyseur;

  if (typ(p) == t_VECSMALL)
  {
    if (typ(res) != t_INT || is_bigint(res))
      pari_err(talker2,"not a suitable VECSMALL component",old,mark.start);
    *pt = (GEN)itos(res); return res;
  }
  if (c->full_row)
  {
    if (typ(res) != t_VEC || lg(res) != lg(p))
      pari_err(talker2,"incorrect type or length in matrix assignment",
          old,mark.start);
    for (i=1; i<lg(p); i++)
    {
      GEN p1 = gcoeff(p,c->full_row,i); if (isclone(p1)) killbloc(p1);
      gcoeff(p,c->full_row,i) = gclone(gel(res,i));
    }
    return res;
  }
  if (c->full_col)
    if (typ(res) != t_COL || lg(res) != lg(*pt))
      pari_err(talker2,"incorrect type or length in matrix assignment",
          old,mark.start);

  res = gclone(res); avma = av;
  killbloc(*pt);
  return *pt = res;
}

/* extract from p the needed component, and assign result if needed */
static GEN
matrix_block(GEN p)
{
  matcomp c;
  GEN cpt = matcell(p, &c);

  if (*analyseur != ',' && *analyseur != ')') /* fast shortcut */
  {
    pari_sp av = avma;
    GEN res;
    F2GEN fun = affect_block(&res);
    if (res)
    {
      if (fun) res = fun(cpt, res);
      return change_compo(av, &c,res);
    }
  }
#if 0
  return isonstack(cpt)? gcopy(cpt): cpt; /* no assignment */
#else
  /* too dangerous otherwise: e.g we return x[1] and have x = 0 immediately
   * after, destroying x[1] in changevalue */
  return gcopy(cpt); /* no assignment */
#endif
}

/* x = gen_0: no default value, otherwise a t_STR, formal expression for
 * default argument. Evaluate and return. */
static GEN
make_arg(GEN x) { return (x==gen_0)? x: readseq(GSTR(x)); }

static GEN
fun_seq(char *t) /* readseq0, simplified */
{
  pari_sp av = top - avma;
  char *olds = analyseur, *olde = mark.start;
  GEN z;

  HANDLE_FOREIGN(t);

  seq_init(t); z = seq();
  analyseur = olds; mark.start = olde;
  av = top - av; /* safer than recording av = avma: f() may call allocatemem */
  if (br_status)
  {
    br_status = br_NONE;
    if (br_res) return gerepilecopy(av, br_res);
    if (!z) { avma = av; return gnil; }
  }
  return z == gnil? z: gerepilecopy(av, z);
}

/* p = NULL + array of variable numbers (longs) + function text */
static GEN
call_fun(entree *ep, GEN *arg)
{
  gp_args *f = (gp_args*)ep->args;
  GEN res, p = (GEN)ep->value, bloc = p, *loc = f->arg + f->narg;
  long i;

  gclone_refc(bloc); /* protect bloc while we use it */
  p++; /* skip NULL */
  /* push new values for formal parameters */
  for (i=0; i<f->narg; i++) copyvalue(*p++, *arg++);
  for (i=0; i<f->nloc; i++) pushvalue(*p++, make_arg(*loc++));
  /* dumps arglist from identifier() to the garbage zone */
  res = fun_seq((char *)p);
  /* pop out values of formal parameters */
  for (i=0; i < f->nloc + f->narg; i++) killvalue(*--p);
  gunclone(bloc); return res;
}
/* p = NULL + array of variable numbers (longs) + function text */
static GEN
call_member(GEN p, GEN x)
{
  GEN res;

  p++; /* skip NULL */
  /* push new values for formal parameters */
  pushvalue(*p++, x);
  res = fun_seq((char *)p);
  /* pop out ancient values of formal parameters */
  killvalue(*--p);
  return res;
}

entree *
do_alias(entree *ep)
{
  while (ep->valence == EpALIAS) ep = (entree *) ((GEN)ep->value)[1];
  return ep;
}

static GEN
global0()
{
  GEN res = gnil;
  long i,n;

  for (i=0,n=lg(polvar)-1; n>=0; n--)
  {
    entree *ep = varentries[n];
    if (ep && EpVALENCE(ep) == EpGVAR)
    {
      res = new_chunk(1);
      gel(res,0) = pol_x[n]; i++;
    }
  }
  if (i) { res = cgetg(1,t_VEC); setlg(res, i+1); }
  return res;
}

static void
check_pointers(gp_pointer ptrs[], unsigned int ind)
{
  unsigned int i;
  pari_sp av = avma;
  for (i=0; i<ind; i++)
  {
    gp_pointer *g = &(ptrs[i]);
    if (g->ep) changevalue(g->ep, g->x);
    else (void)change_compo(av, &(g->c), g->x);
  }
}

#define match_comma() \
  STMT_START { if (matchcomma) match(','); else matchcomma = 1; } STMT_END

static void
skipdecl(void)
{
  if (*analyseur == ':') { analyseur++; skipexpr(); }
}

#define skip_arg() STMT_START {\
  if (do_switch(0,matchcomma))\
    matchcomma = 1;\
  else { match_comma(); skipexpr(); skipdecl(); }\
} STMT_END

static void
skip_arg_block(gp_args *f)
{
  int i, matchcomma = 0;
  for (i = f->narg; i; i--) skip_arg();
}

static long
check_args()
{
  long nparam = 0, matchcomma = 0;
  entree *ep;
  char *old;
  GEN cell;

  match('(');
  while (*analyseur != ')')
  {
    old=analyseur; nparam++; match_comma();
    cell = new_chunk(2);
    if (!isalpha((int)*analyseur))
    {
      err_new_fun();
      pari_err(paramer1, mark.identifier, mark.start);
    }
    ep = entry();
    if (EpVALENCE(ep) != EpVAR)
    {
      err_new_fun();
      if (EpVALENCE(ep) == EpGVAR)
        pari_err(talker2,"global variable: ",old , mark.start);
      pari_err(paramer1, old, mark.start);
    }
    cell[0] = varn(initial_value(ep));
    skipdecl();
    if (*analyseur == '=')
    {
      char *old = ++analyseur;
      pari_sp av = avma;
      skipexpr();
      gel(cell,1) = gclone(_strtoGENstr(old, analyseur-old));
      avma = av;
    }
    else gel(cell,1) = gen_0;
  }
  analyseur++; /* match(')') */
  return nparam;
}

/* function is ok. record it */
static void
record_fun(entree *ep, char *start, long len, long narg, long nloc, GEN tmpargs)
{
  long i, NARG = narg + nloc;
  long L1 = NARG + nchar2nlong(len+1) + 1; /* args + fun code + codeword */
  long L2 = NARG + nchar2nlong(sizeof(gp_args)); /* dflt args */
  GEN newfun, *defarg, ptr = (GEN) newbloc(L1 + L2);
  gp_args *f = (gp_args*)(ptr + L1);

  newfun = ptr;
  *newfun++ = evaltyp(t_STR) | evallg(L1 + L2); /* non-recursive dummy */

  ep->args = (void*) f;
  f->nloc = nloc;
  f->narg = narg;
  f->arg = defarg = (GEN*)(f + 1);

  /* record default args and local variables */
  for (i = 1; i <= NARG; i++)
  {
    GEN cell = tmpargs-(i<<1);
    *newfun++ = cell[0];
    *defarg++ = gel(cell,1);
  }
  /* record text */
  strncpy((char *)newfun, start, len);
  ((char *) newfun)[len] = 0;

  if (NARG > 1)
  { /* check for duplicates */
    GEN x = new_chunk(NARG), v = ptr+1;
    long k;
    for (i=0; i<NARG; i++) x[i] = v[i];
    qsort(x,NARG,sizeof(long),(QSCOMP)pari_compare_long);
    for (k=x[0],i=1; i<NARG; k=x[i],i++)
      if (x[i] == k)
        pari_err(talker,"user function %s: variable %Z declared twice",
            ep->name, pol_x[k]);
  }
  ep->value = ptr;
  ep->valence = EpUSER;
}

static GEN
do_call(void *call, GEN x, GEN argvec[])
{
  return ((PFGEN)call)(x, argvec[1], argvec[2], argvec[3], argvec[4],
                          argvec[5], argvec[6], argvec[7], argvec[8]);
}

/* Rationale: (f(2^-e) - f(-2^-e) + O(2^-pr)) / (2 * 2^-e) = f'(0) + O(2^-2e)
 * since 2nd derivatives cancel.
 *   prec(LHS) = pr - e
 *   prec(RHS) = 2e, equal when  pr = 3e = 3/2 fpr (fpr = required final prec)
 *
 * For f'(x), x far from 0: prec(LHS) = pr - e - expo(x)
 * --> pr = 3/2 fpr + expo(x) */
static GEN
num_deriv(void *call, GEN argvec[])
{
  GEN eps,a,b, y, x = argvec[0];
  long fpr, pr, l, e, ex;
  pari_sp av = avma;
  if (!is_const_t(typ(x))) pari_err(impl, "formal derivation");
  fpr = precision(x)-2; /* required final prec (in sig. words) */
  if (fpr == -2) fpr = precreal-2;
  ex = gexpo(x);
  if (ex < 0) ex = 0; /* at 0 */
  pr = (long)ceil(fpr * 1.5 + (ex / BITS_IN_LONG));
  l = 2+pr;
  e = fpr * BITS_IN_HALFULONG; /* 1/2 required prec (in sig. bits) */

  eps = real2n(-e, l);
  y = gtofp(gsub(x, eps), l); a = do_call(call, y, argvec);
  y = gtofp(gadd(x, eps), l); b = do_call(call, y, argvec);
  setexpo(eps, e-1);
  return gerepileupto(av, gmul(gsub(b,a), eps));
}

/* as above, for user functions */
static GEN
num_derivU(entree *ep, GEN *arg)
{
  GEN eps,a,b, x = *arg;
  long fpr, pr, l, e, ex;
  pari_sp av = avma;

  if (!is_const_t(typ(x))) pari_err(impl, "formal derivation");
  fpr = precision(x)-2; /* required final prec (in sig. words) */
  if (fpr == -2) fpr = precreal-2;
  ex = gexpo(x);
  if (ex < 0) ex = 0; /* at 0 */
  pr = (long)ceil(fpr * 1.5 + (ex / BITS_IN_LONG));
  l = 2+pr;
  e = fpr * BITS_IN_HALFULONG; /* 1/2 required prec (in sig. bits) */

  eps = real2n(-e, l);
  *arg = gtofp(gsub(x, eps), l); a = call_fun(ep,arg);
  *arg = gtofp(gadd(x, eps), l); b = call_fun(ep,arg);
  setexpo(eps, e-1);
  return gerepileupto(av, gmul(gsub(b,a), eps));
}

#define DFT_VAR (GEN)-1L
#define DFT_GEN (GEN)NULL
#define _ARGS_ argvec[0], argvec[1], argvec[2], argvec[3],\
               argvec[4], argvec[5], argvec[6], argvec[7], argvec[8]

static GEN
identifier(void)
{
  long m, i, matchcomma, deriv;
  pari_sp av;
  char *ch1;
  entree *ep;
  GEN res, newfun, ptr;

  mark.identifier = analyseur; ep = entry();
  if (EpVALENCE(ep)==EpVAR || EpVALENCE(ep)==EpGVAR)
  { /* optimized for simple variables */
    switch (*analyseur)
    {
      case ')': case ',': return (GEN)ep->value;
      case '.':
      {
        long n, len, v;
        char *name;

        analyseur++; name = analyseur;
        if ((res = read_member((GEN)ep->value)))
        {
          if (*analyseur == '[')
          {
            matcomp c;
            res = matcell(res, &c);
          }
          return res;
        }
        /* define a new member function */
        v = varn(initial_value(ep));
        len = analyseur - name;
        analyseur++; /* skip = */
        ch1 = name;
        ep = installep(NULL,name,len,EpMEMBER,0,
                       members_hash + hashvalue(&ch1));
        ch1 = analyseur; skipseq(); len = analyseur-ch1;

        n = 2 + nchar2nlong(len+1);
        newfun=ptr= (GEN) newbloc(n);
        *newfun++ = evaltyp(t_STR) | evallg(n); /* non-recursive dummy */
        *newfun++ = v;

        /* record text */
        strncpy((char *)newfun, ch1, len);
        ((char *) newfun)[len] = 0;
        ep->value = (void *)ptr; return gnil;
      }
    }
    if (*analyseur != '[')
    { /* whole variable, no component */
      F2GEN fun = affect_block(&res);
      if (res)
      {
        if (fun) res = fun((GEN)ep->value, res);
        changevalue(ep,res);
      }
      return (GEN)ep->value;
    }
    return matrix_block((GEN)ep->value);
  }
  ep = do_alias(ep);
#ifdef STACK_CHECK
  if (PARI_stack_limit && (void*) &ptr <= PARI_stack_limit)
      pari_err(talker2, "deep recursion", mark.identifier, mark.start);
#endif

  if (ep->code)
  {
    char *s = ep->code, *oldanalyseur = NULL, *buf, *limit, *bp;
    unsigned int ret, noparen, ind_pointer=0;
    long fake;
    void *call = ep->value;
    GEN argvec[9];
    gp_pointer ptrs[9];
    char *flags = NULL;

    deriv = (*analyseur == '\'' && analyseur[1] == '(') && analyseur++;
    if (*analyseur == '(')
    {
      analyseur++;
      noparen=0; /* expect matching ')' */
    }
    else
    { /* if no mandatory argument, no () needed */
      if (EpVALENCE(ep)) match('('); /* error */

      if (!*s || (!s[1] && *s == 'p'))
	return ((GEN (*)(long))call)(precreal);
      noparen=1; /* no argument, but valence is ok */
    }
    /* return type */
    if      (*s <  'a')   ret = RET_GEN;
    else if (*s == 'v') { ret = RET_VOID; s++; }
    else if (*s == 'i') { ret = RET_INT;  s++; }
    else if (*s == 'l') { ret = RET_LONG; s++; }
    else                  ret = RET_GEN;
    /* Optimized for G and p. */
    i = 0;
    matchcomma = 0;
    while (*s == 'G')
    {
      s++;
      match_comma(); ch1 = analyseur;
      argvec[i++] = expr();
      NO_BREAK("here (reading arguments)", ch1);
    }
    if (*s == 'p') { argvec[i++] = (GEN) precreal; s++; }

    while (*s && *s != '\n')
      switch (*s++)
      {
	case 'G': /* GEN */
	  match_comma(); ch1 = analyseur;
          argvec[i++] = expr();
          NO_BREAK("here (reading arguments)", ch1);
          break;

	case 'L': /* long */
	  match_comma(); argvec[i++] = (GEN) readlong(); break;

	case 'n': /* var number */
	  match_comma(); argvec[i++] = (GEN) readvar(); break;

        case 'S': /* symbol */
	  match_comma(); mark.symbol=analyseur;
	  argvec[i++] = (GEN)entry(); break;

	case 'V': /* variable */
	  match_comma(); mark.symbol=analyseur;
        {
          entree *e = entry();
          long v = EpVALENCE(e);
          if (v != EpVAR && v != EpGVAR)
            pari_err(talker2,"not a variable:",mark.symbol,mark.start);
	  argvec[i++] = (GEN)e; break;
        }
        case '&': /* *GEN */
	  match_comma(); match('&'); mark.symbol=analyseur;
        {
          entree *ep = entry();
          gp_pointer *g = &ptrs[ind_pointer++];

          if (*analyseur == '[') {
            (void)matcell((GEN)ep->value, &(g->c));
            g->x = *(g->c.ptcell);
            g->ep = NULL;
          } else {
            g->x = (GEN) ep->value;
            g->ep = ep;
          }
	  argvec[i++] = (GEN)&(g->x); break;
        }
        /* Input position */
        case 'E': /* expr */
	case 'I': /* seq */
	  match_comma();
	  argvec[i++] = (GEN) analyseur;
	  skipseq(); break;

	case 'r': /* raw */
	  match_comma(); mark.raw = analyseur;
          bp = init_buf(256, &buf,&limit);
	  while (*analyseur)
	  {
	    if (*analyseur == ',' || *analyseur == ')') break;
	    if (*analyseur == '"')
	      bp = readstring_i(bp, &buf,&limit);
            else
            {
              if (bp > limit)
                bp = realloc_buf(bp,1, &buf,&limit);
              *bp++ = *analyseur++;
            }
	  }
	  *bp++ = 0; argvec[i++] = (GEN) buf;
	  break;

	case 'M': /* Mnemonic flag */
	  match_comma(); ch1 = analyseur;
          argvec[i] = expr();
          NO_BREAK("here (reading arguments)", ch1);
	  if (typ(argvec[i]) == t_STR) {
	    if (!flags) flags = ep->code;
	    flags = strchr(flags, '\n'); /* Skip to the following '\n' */
	    if (!flags)
	        pari_err(talker, "not enough flags in string function signature");
	    flags++;
	    argvec[i] = (GEN) parse_option_string((char*)(argvec[i] + 1),
			flags, PARSEMNU_ARG_WHITESP | PARSEMNU_TEMPL_TERM_NL,
			NULL, NULL);
	  } else
            argvec[i] = (GEN)itos(argvec[i]);
	  i++;
          break;

	case 's': /* expanded string; empty arg yields "" */
	  match_comma();
	  if (*s == '*') /* any number of string objects */
          {
            argvec[i++] = noparen? cgetg(1, t_VEC): any_string();
            s++; break;
          }

          bp = init_buf(256, &buf,&limit);
          while (*analyseur)
          {
            if (*analyseur == ',' || *analyseur == ')') break;
            bp = expand_string(bp, &buf,&limit);
          }
          *bp++ = 0; argvec[i++] = (GEN)buf;
          break;

	case 'p': /* precision */
	  argvec[i++] = (GEN) precreal; break;

	case '=':
	  match('='); matchcomma = 0; break;

	case 'D': /* Has a default value */
	  if (do_switch(noparen,matchcomma))
            switch (*s)
            {
              case 'G':
              case '&':
              case 'I':
              case 'V': matchcomma=1; argvec[i++]=DFT_GEN; s++; break;
              case 'n': matchcomma=1; argvec[i++]=DFT_VAR; s++; break;
              default:
                oldanalyseur = analyseur;
                analyseur = s; matchcomma = 0;
                while (*s++ != ',');
            }
          else
            switch (*s)
            {
              case 'G':
              case '&':
              case 'I':
              case 'V':
              case 'n': break;
              default:
                while (*s++ != ',');
            }
	  break;

	 case 'P': /* series precision */
	   argvec[i++] = (GEN) precdl; break;

	 case 'f': /* Fake *long argument */
	   argvec[i++] = (GEN) &fake; break;

	 case 'x': /* Foreign function */
	   argvec[i++] = (GEN) ep; call = foreignHandler; break;

	 case ',': /* Clean up default */
	   if (oldanalyseur)
	   {
	     analyseur = oldanalyseur;
	     oldanalyseur = NULL; matchcomma=1;
	   }
	   break;
	 default: pari_err(bugparier,"identifier (unknown code)");
      }
#if 0 /* uncomment if using purify: UMR otherwise */
    for ( ; i<9; i++) argvec[i]=NULL;
#endif
{
    char *oldname = gp_function_name;
    gp_function_name = ep->name;
    if (deriv)
    {
      if (!i || (ep->code)[0] != 'G')
        pari_err(talker2, "can't derive this", mark.identifier, mark.start);
      res = num_deriv(call, argvec);
    }
    else switch (ret)
      {
      case RET_GEN:
        res = ((PFGEN)call)(_ARGS_);
        break;

      case RET_INT:
	m = (long)((int (*)(ANYARG))call)(_ARGS_);
	res = stoi(m); break;

      case RET_LONG:
	m = ((long (*)(ANYARG))call)(_ARGS_);
	res = stoi(m); break;

      default: /* RET_VOID */
	((void (*)(ANYARG))call)(_ARGS_);
	res = gnil; break;
      }
    gp_function_name = oldname;
}
    if (ind_pointer) check_pointers(ptrs, ind_pointer);
    if (!noparen) match(')');
    return res;
  }

  if (EpPREDEFINED(ep))
  {
    char *oldname=gp_function_name;
    if (*analyseur != '(')
    {
      if (EpVALENCE(ep) == 88) return global0();
      match('('); /* error */
    }
    analyseur++; ch1 = analyseur;
    switch(EpVALENCE(ep))
    {
      case 50: /* O */
        gp_function_name="O";
        res = truc();
        NO_BREAK("in O()", ch1);
	if (*analyseur=='^') { analyseur++; m = readlong(); } else m = 1;
	res = ggrando(res,m); break;

      case 80: /* if then else */
        gp_function_name="if";
        av = avma; res = expr();
        NO_BREAK("in test expression", ch1);
        m = gcmp0(res); avma = av; match(',');
	if (m) /* false */
	{
	  skipseq();
	  if (*analyseur == ')') res = gnil;
	  else
          {
            match(',');
            res = seq(); if (br_status) { res = NULL; skipseq(); }
          }
	}
	else /* true */
	{
          res = seq(); if (br_status) { res = NULL; skipseq(); }
          if (*analyseur != ')') { match(','); skipseq(); }
	}
	break;

      case 81: /* while do */
        gp_function_name="while";
        av = avma;
	for(;;)
	{
          res = expr();
          NO_BREAK("in test expression", ch1);
	  if (gcmp0(res)) { match(','); break; }

	  avma = av; match(','); (void)seq();
	  if (loop_break()) break;
          analyseur = ch1;
	}
	avma = av; skipseq(); res = gnil; break;

      case 82: /* repeat until */
        gp_function_name="until";
        av = avma; skipexpr();
	for(;;)
	{
	  avma = av; match(','); (void)seq();
	  if (loop_break()) break;
	  analyseur = ch1;
          res = expr();
          NO_BREAK("in test expression", ch1);
	  if (!gcmp0(res)) { match(','); break; }
	}
	avma = av; skipseq(); res = gnil; break;

      case 88: /* global */
        gp_function_name="global";
        if (*analyseur == ')') return global0();
        matchcomma = 0;
        while (*analyseur != ')')
        {
          match_comma(); ch1 = analyseur;
          check_var_name();
          ep = skipentry();
          switch(EpVALENCE(ep))
          {
            case EpGVAR:
            case EpVAR: break;
            default: pari_err(talker2,"symbol already in use",ch1,mark.start);
          }
          analyseur = ch1; ep = entry();
          if (*analyseur == '=')
          {
            pari_sp av = avma; analyseur++;
            ch1 = analyseur;
            res = expr();
            NO_BREAK("here (defining global var)", ch1);
            changevalue(ep, res); avma = av;
          }
          ep->valence = EpGVAR;
        }
        res = gnil; break;

      default: pari_err(valencer1);
        return NULL; /* not reached */
    }
    gp_function_name=oldname;
    match(')'); return res;
  }

  switch (EpVALENCE(ep))
  {
    case EpUSER: /* user-defined functions */
    {
      GEN *arglist;
      gp_args *f = (gp_args*)ep->args;
      deriv = (*analyseur == '\'' && analyseur[1] == '(') && analyseur++;
      arglist = (GEN*) new_chunk(f->narg);
      if (*analyseur != '(') /* no args */
      {
	if (*analyseur != '='  ||  analyseur[1] == '=')
        {
          for (i=0; i<f->narg; i++) arglist[i] = make_arg(f->arg[i]);
	  return call_fun(ep, arglist);
        }
	match('('); /* ==> error */
      }
      analyseur++; /* skip '(' */
      ch1 = analyseur;
      skip_arg_block(f);
      if (*analyseur == ')' && (analyseur[1] != '=' || analyseur[2] == '='))
      {
        matchcomma = 0;
        analyseur = ch1;
        for (i=0; i<f->narg; i++)
        {
          if (do_switch(0,matchcomma))
          { /* default arg */
            arglist[i] = make_arg(f->arg[i]);
            matchcomma = 1;
          }
          else
          { /* user supplied */
            char *old;
            match_comma(); old = analyseur;
            arglist[i] = expr();
            NO_BREAK("here (reading function args)", old);
          }
        }
        analyseur++; /* skip ')' */
        if (deriv)
        {
          if (!f->narg)
            pari_err(talker2, "can't derive this", mark.identifier, mark.start);
          return num_derivU(ep, arglist);
        }
        return call_fun(ep, arglist);
      }
      if (*analyseur != ',' && *analyseur != ')') skipexpr();
      while (*analyseur == ',') { analyseur++; skipexpr(); }
      match(')');
      if (*analyseur != '=' || analyseur[1] == '=')
        pari_err(talker2,"too many parameters in user-defined function call",
            mark.identifier,mark.start);

      analyseur = ch1-1; /* points to '(' */
    } /* Fall through, REDEFINE function */

    case EpNEW: /* new function */
    {
      GEN tokill = NULL, tmpargs = (GEN)avma;
      char *start;
      long narg, nloc;

      if (ep->valence == EpUSER)
      {
        ep->valence = EpNEW;
        tokill = (GEN)ep->value; /* can't kill now */
        free_ep_args(ep);
      }
      check_new_fun = ep;

      /* check arguments */
      narg = check_args(); nloc = 0;
      /* Dirty, but don't want to define a local() function */
      if (*analyseur != '=' && strcmp(ep->name, "local") == 0)
        pari_err(talker2, "local() bloc must appear before any other expression",
                     mark.identifier,mark.start);
      if (*analyseur != '=')
      {
        char *str = stackmalloc(128 + strlen(ep->name));
        sprintf(str,"unknown function '%s', expected '=' instead of", ep->name);
        pari_err(talker2,str, analyseur, mark.start);
      }
      analyseur++;

      /* checking function definition */
      skipping_fun_def++;
      while (strncmp(analyseur,"local(",6) == 0)
      {
        analyseur += 5; /* on '(' */
        nloc += check_args();
        while(separator(*analyseur)) analyseur++;
      }
      start = analyseur; skipseq(); skipping_fun_def--;
      record_fun(ep, start, analyseur-start, narg, nloc, tmpargs);
     /* wait till here for gunclone. In pathological cases, e.g. (f()=f()=x),
      * new text is given by value of old one! */
      if (tokill) gunclone(tokill);
      check_new_fun = NULL;
      avma = (pari_sp)tmpargs; return gnil;
    }
  }
  pari_err(valencer1); return NULL; /* not reached */
}

static ulong
number(int *n, char **s)
{
  ulong m = 0;
  for (*n = 0; *n < 9 && isdigit((int)**s); (*n)++,(*s)++)
    m = 10*m + (**s - '0');
  return m;
}

ulong
u_pow10(int n)
{
  static ulong pw10[] = { 1UL, 10UL, 100UL, 1000UL, 10000UL, 100000UL,
                        1000000UL, 10000000UL, 100000000UL, 1000000000UL };
  return pw10[n];
}

static GEN
int_read_more(GEN y, char **ps)
{
  pari_sp av = avma;
  int i = 0, nb;
  while (isdigit((int)**ps))
  {
    ulong m = number(&nb, ps);
    if (avma != av && ++i > 4) { avma = av; i = 0; } /* HACK gerepile */
    y = addumului(m, u_pow10(nb), y);
  }
  return y;
}

static long
exponent(char **pts)
{
  char *s = *pts;
  long n;
  int nb;
  switch(*++s)
  {
    case '-': s++; n = -(long)number(&nb, &s); break;
    case '+': s++; /* Fall through */
    default: n = (long)number(&nb, &s);
  }
  *pts = s; return n;
}

static GEN
real_0_digits(long n) {
  long b = (n > 0)? (long)(n/L2SL10): (long)-((-n)/L2SL10 + 1);
  return real_0_bit(b);
}

static GEN
real_read(pari_sp av, char **s, GEN y, long PREC)
{
  long l, n = 0;
  switch(**s)
  {
    default: return y; /* integer */
    case '.':
    {
      char *old = ++*s;
      if (isalpha((int)**s))
      {
        if (**s == 'E' || **s == 'e') {
          n = exponent(s);
          if (!signe(y)) { avma = av; return real_0_digits(n); }
          break;
        }
        --*s; return y; /* member */
      }
      y = int_read_more(y, s);
      n = old - *s;
      if (**s != 'E' && **s != 'e')
      {
        if (!signe(y)) { avma = av; return real_0(PREC); }
        break;
      }
    }
    /* Fall through */
    case 'E': case 'e':
      n += exponent(s);
      if (!signe(y)) { avma = av; return real_0_digits(n); }
  }
  l = lgefint(y); if (l < (long)PREC) l = (long)PREC;
  if (!n) return itor(y, l);
  y = itor(y, l+1);
  if (n > 0)
    y = mulrr(y, rpowuu(10UL, (ulong)n, l+1));
  else
    y = divrr(y, rpowuu(10UL, (ulong)-n, l+1));
  return gerepileuptoleaf(av, rtor(y, l));
}

static GEN
int_read(char **s)
{
  int nb;
  GEN y = utoi(number(&nb, s));
  if (nb == 9) y = int_read_more(y, s);
  return y;
}

GEN
strtoi(char *s) { return int_read(&s); }

GEN 
strtor(char *s, long PREC)
{
  pari_sp av = avma;
  GEN y = int_read(&s);
  y = real_read(av, &s, y, PREC);
  if (typ(y) == t_REAL) return y;
  return gerepileuptoleaf(av, itor(y, PREC));
}

static GEN
constante()
{
  pari_sp av = avma;
  GEN y = int_read(&analyseur);
  return real_read(av, &analyseur, y, precreal);
}

/********************************************************************/
/**                                                                **/
/**                   HASH TABLE MANIPULATIONS                     **/
/**                                                                **/
/********************************************************************/
/* inline is_keyword_char(). Not worth a static array. */
#define is_key(c) (isalnum((int)(c)) || (c)=='_')

long
is_keyword_char(char c) { return is_key(c); }

/* return hashing value for identifier s */
long
hashvalue(char **str)
{
  long n = 0;
  char *s = *str;
  while (is_key(*s)) { n = (n<<1) ^ *s; s++; }
  *str = s; if (n < 0) n = -n;
  return n % functions_tblsz;
}

/* Looking for entry in hashtable. ep1 is the cell's first element */
static entree *
findentry(char *name, long len, entree *ep1)
{
  entree *ep;

  for (ep = ep1; ep; ep = ep->next)
    if (!strncmp(ep->name, name, len) && !(ep->name)[len]) return ep;

  if (foreignAutoload) /* Try to autoload. */
    return foreignAutoload(name, len);
  return NULL; /* not found */
}

entree *
is_entry(char *s)
{
  return is_entry_intern(s,functions_hash,NULL);
}

entree *
is_entry_intern(char *s, entree **table, long *pthash)
{
  char *t = s;
  long hash = hashvalue(&t); 
  if (pthash) *pthash = hash;
  return findentry(s, t - s, table[hash]);
}

int
is_identifier(char *s)
{
  while (*s && is_keyword_char(*s)) s++;
  return *s? 0: 1;
}

static entree *
installep(void *f, char *name, long len, long valence, long add, entree **table)
{
  entree *ep = (entree *) gpmalloc(sizeof(entree) + add + len+1);
  const entree *ep1 = initial_value(ep);
  char *u = (char *) ep1 + add;

  ep->name    = u; strncpy(u, name,len); u[len]=0;
  ep->args    = INITIAL; /* necessary, see var_cell definition */
  ep->help    = NULL;
  ep->code    = NULL;
  ep->value   = f? f: (void *) ep1;
  ep->next    = *table;
  ep->valence = valence;
  ep->menu    = 0;
  return *table = ep;
}

long
manage_var(long n, entree *ep)
{
  static long max_avail = MAXVARN; /* max variable not yet used */
  static long nvar; /* first GP free variable */
  long var;
  GEN p;

  switch(n) {
      case manage_var_init: return nvar=0;
      case manage_var_next: return nvar;
      case manage_var_max_avail: return max_avail;
      case manage_var_pop:
      {
        long v = (long)ep;
        if (v != nvar-1) pari_err(talker,"can't pop gp variable");
        setlg(polvar, nvar);
        return --nvar;
      }
      case manage_var_delete:
	/* user wants to delete one of his/her/its variables */
	if (max_avail == MAXVARN-1) return 0; /* nothing to delete */
	free(pol_x[++max_avail]); /* frees both pol_1 and pol_x */
	return max_avail+1;
      case manage_var_create: break;
      default: pari_err(talker, "panic");
  }

  if (nvar == max_avail) pari_err(talker2,"no more variables available",
                             mark.identifier, mark.start);
  if (ep)
  {
    p = (GEN)ep->value;
    var=nvar++;
  }
  else
  {
    p = (GEN) gpmalloc(SIZEOF_VAR_POLS);
    var=max_avail--;
  }

  /* create pol_x[var] */
  p[0] = evaltyp(t_POL) | evallg(4);
  p[1] = evalsigne(1) | evalvarn(var);
  gel(p,2) = gen_0;
  gel(p,3) = gen_1;
  pol_x[var] = p;

  /* create pol_1[nvar] */
  p += 4;
  p[0] = evaltyp(t_POL) | evallg(3);
  p[1] = evalsigne(1) | evalvarn(var);
  gel(p,2) = gen_1;
  pol_1[var] = p;

  varentries[var] = ep;
  if (ep) { gel(polvar,nvar) = (GEN)ep->value; setlg(polvar, nvar+1); }
  return var;
}

long
fetch_var(void)
{
  return manage_var(manage_var_create,NULL);
}

entree *
fetch_named_var(char *s)
{
  char *t = s;
  entree **funhash = functions_hash + hashvalue(&t);
  entree *ep = findentry(s, t - s, *funhash);
  if (ep)
  {
    switch (EpVALENCE(ep))
    {
      case EpVAR: case EpGVAR: break;
      default: pari_err(talker, "%s already exists with incompatible valence", s);
    }
    return ep;
  }
  ep = installep(NULL,s,strlen(s),EpVAR, SIZEOF_VAR_POLS, funhash);
  (void)manage_var(manage_var_create,ep); return ep;
}

long
fetch_user_var(char *s)
{
  return varn( initial_value(fetch_named_var(s)) );
}

void
delete_named_var(entree *ep)
{
  (void)manage_var(manage_var_pop, (entree*)varn(initial_value(ep)));
  kill0(ep);
}

long
delete_var(void)
{
  return manage_var(manage_var_delete,NULL);
}

void
name_var(long n, char *s)
{
  entree *ep;
  char *u;

  if (n < manage_var(manage_var_next,NULL))
    pari_err(talker, "renaming a GP variable is forbidden");
  if (n > (long)MAXVARN)
    pari_err(talker, "variable number too big");

  ep = (entree*)gpmalloc(sizeof(entree) + strlen(s) + 1);
  u = (char *)initial_value(ep);
  ep->valence = EpVAR;
  ep->name = u; strcpy(u,s);
  ep->value = gen_0; /* in case geval is called */
  if (varentries[n]) free(varentries[n]);
  varentries[n] = ep;
}

/* Find entry or create it */
static entree *
entry(void)
{
  char *old = analyseur;
  const long hash = hashvalue(&analyseur), len = analyseur - old;
  entree *ep = findentry(old,len,functions_hash[hash]);
  long val,n;

  if (ep) return ep;
  if (compatible == WARN)
  {
    ep = findentry(old,len,funct_old_hash[hash]);
    if (ep) return ep; /* the warning was done in skipentry() */
  }
  /* ep does not exist. Create it */
  if (*analyseur == '(')
    { n=0; val=EpNEW; }
  else
    { n=SIZEOF_VAR_POLS; val=EpVAR; }
  ep = installep(NULL,old,len,val,n, functions_hash + hash);

  if (n) (void)manage_var(manage_var_create, ep); /* Variable */
  return ep;
}

/********************************************************************/
/**                                                                **/
/**                          SKIP FUNCTIONS                        **/
/**                                                                **/
/********************************************************************/
/* analyseur points on the character following the starting " */
/* skip any number of concatenated strings */
static void
skipstring()
{
  while (*analyseur)
    switch (*analyseur++)
    {
      case '"': if (*analyseur != '"') return;
      /* fall through */
      case '\\': analyseur++;
    }
  match('"');
}

static void
skip_matrix_block()
{
  while (*analyseur == '[')
  {
    analyseur++;
    if (*analyseur == ',') { analyseur++; skipexpr(); }
    else
    {
      skipexpr();
      if (*analyseur == ',')
	if (*++analyseur != ']') skipexpr();
    }
    match(']');
  }
}

/* return 1 if we would be assigning some value after expansion. 0 otherwise.
 * Skip all chars corresponding to the assignment (and assigned value) */
static int
skip_affect_block()
{
  if (*analyseur == '=')
  {
    if (analyseur[1] != '=') { analyseur++; skipexpr(); return 1; }
  }
  else if (double_op()) return 1;
  else if (get_op_fun()) { skipexpr(); return 1; }
  return 0;
}

static void
skipseq(void)
{
  for(;;)
  {
    while (separator(*analyseur)) analyseur++;
    if (*analyseur == ',' || *analyseur == ')' || !*analyseur) return;
    skipexpr(); if (!separator(*analyseur)) return;
  }
}

static void
skipexpr(void)
{
  int e1 = 0, e2 = 0, e3;

L3:
  e3=1; skipfacteur();
  switch(*analyseur)
  {
    case '*': case '/': case '%':
      analyseur++; goto L3;
    case '\\':
      if (*++analyseur == '/') analyseur++;
      goto L3;
    case '<': case '>':
      if (analyseur[1]==*analyseur) { analyseur +=2; goto L3; }
  }

L2:
  if (!e3) goto L3;
  e3=0; e2=1;
  switch(*analyseur)
  {
    case '+': case '-':
      analyseur++; goto L3;
  }

L1:
  if (!e2) goto L2;
  e2=0; e1=1;
  switch(*analyseur)
  {
    case '<':
      switch(*++analyseur)
      {
        case '=': case '>': analyseur++;
      }
      goto L2;

    case '>':
      if (*++analyseur == '=') analyseur++;
      goto L2;

    case '=': case '!':
      if (analyseur[1] == '=') { analyseur+=2; goto L2; }
      goto L1;
  }

/* L0: */
  if (!e1) goto L1;
  e1=0;
  switch(*analyseur)
  {
    case '&':
      if (*++analyseur == '&') analyseur++;
      goto L1;
    case '|':
      if (*++analyseur == '|') analyseur++;
      goto L1;
  }
}

static void
skipmember(void) {
  while (is_key((int)*analyseur)) analyseur++;
}

static void
skipfacteur(void)
{
  if (*analyseur == '+' || *analyseur == '-') analyseur++;
  skiptruc();
  for(;;)
    switch(*analyseur)
    {
      case '.':
	analyseur++; skipmember();
        if (*analyseur == '=' && analyseur[1] != '=')
          { analyseur++; skipseq(); }
        break;
      case '^':
	analyseur++; skipfacteur(); break;
      case '~': case '\'':
	analyseur++; break;
      case '[':
      {
        char *old;
	skip_matrix_block(); old = analyseur;
        if (skip_affect_block()) pari_err(caracer1,old,mark.start);
        break;
      }
      case '!':
	if (analyseur[1] != '=') { analyseur++; break; }
      default: return;
    }
}

/* return the number of elements we need to read if array/matrix */
static void
skiptruc(void)
{
  long i, m, n;
  char *old;

  if (isalpha((int)*analyseur)) { skipidentifier(); return; }
  if (isdigit((int)*analyseur) || *analyseur== '.') { skipconstante(); return; }
  switch(*analyseur)
  {
    case '"':
      analyseur++; skipstring(); return;

    case '!':
    case '#':
      analyseur++; skipfacteur(); return;

    case '&':
    case '\'':
      analyseur++; check_var_name();
      (void)skipentry(); return;

    case '(':
      analyseur++; skipexpr(); match(')'); return;

    case '[':
      old = analyseur; analyseur++;
      if (*analyseur == ';' && analyseur[1] == ']')  /* [;] */
        { analyseur+=2; return; }
      n = 0;
      if (*analyseur != ']')
      {
	do { n++; skipexpr(); old=analyseur; } while (*analyseur++ == ',');
	analyseur--;
      }
      switch (*analyseur++)
      {
	case ']': return;
	case ';':
          for (m=2; ; m++)
          {
            for (i=1; i<n; i++) { skipexpr(); match(','); }
            skipexpr();
            if (*analyseur == ']') break;
            match(';');
          }
          analyseur++; return;
	default:
	  pari_err(talker2,"; or ] expected",old,mark.start);
      }

    case '%':
    {
      int junk;
      analyseur++;
      if (*analyseur == '#') { analyseur++; return; }
      if (*analyseur == '`') { while (*++analyseur == '`') /*empty*/; return; }
      (void)number(&junk, &analyseur); return;
    }
  }
  pari_err(caracer1,analyseur,mark.start);
}

static void
check_var()
{
  char *old = analyseur;
  check_var_name();
  switch(EpVALENCE(skipentry()))
  {
    case EpVAR: break;
    case EpGVAR:
      pari_err(talker2,"global variable not allowed", old,mark.start);
    default: pari_err(varer1,old,mark.start);
  }
}

static void
check_matcell()
{
  char *old = analyseur;
  check_var_name();
  switch(EpVALENCE(skipentry()))
  {
    case EpVAR:
    case EpGVAR: break;
    default: pari_err(varer1,old,mark.start);
  }
  skip_matrix_block();
}

static void
skipidentifier(void)
{
  int matchcomma;
  entree *ep;

  mark.identifier = analyseur; ep = do_alias(skipentry());
  if (ep->code)
  {
    char *s = ep->code;

    if (*analyseur == '\'') analyseur++;
    if (*analyseur != '(')
    {
      if (EpVALENCE(ep) == 0) return; /* no mandatory argument */
      match('('); /* ==> error */
    }
    analyseur++;

    if (*s == 'v' || *s == 'l' || *s == 'i') s++;
    /* Optimized for G and p. */
    matchcomma = 0;
    while (*s == 'G') { match_comma(); skipexpr(); s++; }
    if (*s == 'p') s++;
    while (*s && *s != '\n') switch (*s++)
    {
      case 'G': case 'n': case 'L': case 'M':
        match_comma();
        if (*analyseur == ',' || *analyseur == ')') break;
        skipexpr(); break;
      case 'E':
        match_comma(); skipexpr(); break;
      case 'I':
        match_comma(); skipseq(); break;
      case 'r':
        match_comma();
        while (*analyseur)
        {
          if (*analyseur == '"') { analyseur++; skipstring(); }
          if (*analyseur == ',' || *analyseur == ')') break;
          analyseur++;
        }
        break;
      case 's':
        match_comma();
        if (*s == '*')
        {
          while (*analyseur)
          {
            if (*analyseur == '"') { analyseur++; skipstring(); }
            if (*analyseur == ')') break;
            if (*analyseur == ',') analyseur++;
            else skipexpr();
          }
          s++;
          break;
        }

        while (*analyseur)
        {
          if (*analyseur == '"') { analyseur++; skipstring(); }
          if (*analyseur == ',' || *analyseur == ')') break;
          skipexpr();
        }
        break;

      case 'S': match_comma();
        check_var_name(); (void)skipentry(); break;
      case '&': match_comma(); match('&'); check_matcell(); break;
      case 'V': match_comma(); check_var(); break;

      case 'p': case 'P': case 'f': case 'x':
        break;

      case 'D':
        if ( *analyseur == ')' ) { analyseur++; return; }
        if (*s == 'G' || *s == '&' || *s == 'n' || *s == 'I' || *s == 'V')
          break;
        while (*s++ != ',');
        break;
      case '=':
        match('='); matchcomma = 0; break;
      case ',':
        matchcomma=1; break;
      case '\n':			/* Before the mnemonic */
	break;
      default:
        pari_err(bugparier,"skipidentifier (unknown code)");
    }
    match(')');
    return;
  }
  if (EpPREDEFINED(ep))
  {
    if (*analyseur != '(')
    {
      switch(EpVALENCE(ep))
      {
        case 0:
        case 88: return;
      }
      match('('); /* error */
    }
    analyseur++;
    switch(EpVALENCE(ep))
    {
      case 50: skiptruc();
	if (*analyseur == '^') { analyseur++; skipfacteur(); };
	break;
      case 80: skipexpr(); match(','); skipseq();
          if (*analyseur != ')') { match(','); skipseq(); }
	  break;
      case 81: case 82: skipexpr(); match(','); skipseq(); break;
      case 88:
        matchcomma = 0;
        while (*analyseur != ')') { match_comma(); skipexpr(); };
        break;
      default: pari_err(valencer1);
    }
    match(')'); return;
  }
  switch (EpVALENCE(ep))
  {
    case EpGVAR:
    case EpVAR: /* variables */
      skip_matrix_block(); (void)skip_affect_block(); return;

    case EpUSER: /* fonctions utilisateur */
    {
      char *ch1 = analyseur;

      if (*analyseur == '\'') analyseur++;
      if (*analyseur != '(')
      {
	if ( *analyseur != '=' || analyseur[1] == '=' ) return;
	match('('); /* error */
      }
      analyseur++;  /* skip '(' */
      skip_arg_block((gp_args*)ep->args);
      if (*analyseur == ')' && (analyseur[1] != '=' || analyseur[2] == '='))
	  { analyseur++; return; }

      /* here we are redefining a user function */
      if (*analyseur != ',' && *analyseur != ')') skipexpr();
      while (*analyseur == ',') { analyseur++; skipexpr(); }
      match(')');

      if (*analyseur != '=' || analyseur[1] == '=')
      {
        if (skipping_fun_def) return;
        pari_err(talker2,"too many parameters in user-defined function call",
            mark.identifier,mark.start);
      }
      analyseur = ch1;
    } /* fall through */

    case EpNEW: /* new function */
      if (check_new_fun && ! skipping_fun_def)
      {
	err_new_fun(); /* ep not created yet: no need to kill it */
	pari_err(paramer1, mark.identifier, mark.start);
      }
      check_new_fun = NOT_CREATED_YET; match('(');
      matchcomma = 0;
      while (*analyseur != ')') {
        if (!*analyseur) match(')'); /* error */
        skip_arg();
      }
      analyseur++; /* skip ')' */
      if (*analyseur == '=' && analyseur[1] != '=')
      {
	skipping_fun_def++;
	analyseur++; skipseq();
	skipping_fun_def--;
      }
      check_new_fun=NULL; return;

    default: pari_err(valencer1);
  }
}

static void
skipdigits(void) {
  while (isdigit((int)*analyseur)) analyseur++;
}
static void
skipexponent(void) {
  if (*analyseur=='e' || *analyseur=='E')
  {
    analyseur++;
    if ( *analyseur=='+' || *analyseur=='-' ) analyseur++;
    skipdigits();
  }
}

static void
skipconstante(void)
{
  skipdigits();
  if (*analyseur=='.')
  {
    char *old = ++analyseur;
    if (isalpha((int)*analyseur))
    {
      skipexponent();
      if (analyseur == old) analyseur--; /* member */
      return;
    }
    skipdigits();
  }
  skipexponent();
}

static entree *
skipentry(void)
{
  static entree fakeEpNEW = { "",EpNEW };
  static entree fakeEpVAR = { "",EpVAR };
  char *old = analyseur;
  const long hash = hashvalue(&analyseur), len = analyseur - old;
  entree *ep = findentry(old,len,functions_hash[hash]);

  if (ep) return ep;
  if (compatible == WARN)
  {
    ep = findentry(old,len,funct_old_hash[hash]);
    if (ep)
    {
      pari_warn(warner,"using obsolete function %s",ep->name);
      return ep;
    }
  }
  return (*analyseur == '(') ? &fakeEpNEW : &fakeEpVAR;
}

#include "members.h"

static entree*
find_member()
{
  char *old = analyseur;
  const long hash = hashvalue(&analyseur), len = analyseur - old;
  return findentry(old,len,members_hash[hash]);
}

static GEN
read_member(GEN x)
{
  entree *ep;

  mark.member = analyseur;
  ep = find_member();
  if (ep)
  {
    if (*analyseur == '=' && analyseur[1] != '=')
    {
      if (EpPREDEFINED(ep))
        pari_err(talker2,"can't modify a pre-defined member: ",
            mark.member,mark.start);
      gunclone((GEN)ep->value); return NULL;
    }
    if (EpVALENCE(ep) == EpMEMBER)
      return call_member((GEN)ep->value, x);
    else
    {
      GEN y = ((F1GEN)ep->value)(x);
      return isonstack(y)? gcopy(y): y;
    }
  }
  if (*analyseur != '=' || analyseur[1] == '=')
    pari_err(talker2,"unknown member function",mark.member,mark.start);
  return NULL; /* to be redefined */
}

void
member_err(char *s)
{
  char *str = stackmalloc(strlen(s) + 128);
  sprintf(str, "incorrect type in %s", s);
  pari_err(talker2,str,mark.member,mark.start);
}

/********************************************************************/
/**                                                                **/
/**                        SIMPLE GP FUNCTIONS                     **/
/**                                                                **/
/********************************************************************/
long
loop_break()
{
  switch(br_status)
  {
    case br_MULTINEXT :
      if (! --br_count) br_status = br_NEXT;
      return 1;
    case br_BREAK : if (! --br_count) br_status = br_NONE; /* fall through */
    case br_RETURN: return 1;
    case br_ALLOCMEM: allocate_loop_err();
    case br_NEXT: br_status = br_NONE; /* fall through */
  }
  return 0;
}

long
did_break() { return br_status; }

GEN
return0(GEN x)
{
  GEN y = br_res;
  br_res = (x && x != gnil)? gclone(x): NULL;
  if (y) gunclone(y);
  br_status = br_RETURN; return NULL;
}

GEN
next0(long n)
{
  if (n < 1)
    pari_err(talker2,"positive integer expected",mark.identifier,mark.start);
  if (n == 1) br_status = br_NEXT;
  else
  {
    br_count = n-1;
    br_status = br_MULTINEXT;
  }
  return NULL;
}

GEN
break0(long n)
{
  if (n < 1)
    pari_err(talker2,"positive integer expected",mark.identifier,mark.start);
  br_count = n;
  br_status = br_BREAK; return NULL;
}

void
allocatemem0(size_t newsize)
{
  (void)allocatemoremem(newsize);
  br_status = br_ALLOCMEM;
}

void
alias0(char *s, char *old)
{
  entree *ep, *e;
  long hash;
  GEN x;

  ep = is_entry(old);
  if (!ep) pari_err(talker2,"unknown function",mark.raw,mark.start);
  switch(EpVALENCE(ep))
  {
    case EpVAR: case EpGVAR:
      pari_err(talker2,"only functions can be aliased",mark.raw,mark.start);
  }

  if ( (e = is_entry_intern(s, functions_hash, &hash)) )
  {
    if (EpVALENCE(e) != EpALIAS)
      pari_err(talker2,"can't replace an existing symbol by an alias",
          mark.raw, mark.start);
    kill0(e);
  }
  ep = do_alias(ep); x = newbloc(2);
  x[0] = evaltyp(t_STR)|evallg(2); /* for getheap */
  gel(x,1) = (GEN)ep;
  (void)installep(x, s, strlen(s), EpALIAS, 0, functions_hash + hash);
}

/********************************************************************/
/**                                                                **/
/**            PRINT USER FUNCTION AND MEMBER FUNCTION             **/
/**                                                                **/
/********************************************************************/

static void
print_def_arg(GEN x)
{
  if (x == gen_0) return;
  pariputc('=');
  if (typ(x)==t_STR)
    pariputs(GSTR(x)); /* otherwise it's surrounded by "" */
  else
    brute(x,'g',-1);
}

static void
print_var(long n) {
  entree *ep = varentries[n];
  pariputs(ep? ep->name:"#");
}

void
print_user_fun(entree *ep)
{
  gp_args *f= (gp_args*)ep->args;
  GEN q = (GEN)ep->value, *arg = f->arg;
  int i, narg;

  q++; /* skip initial NULL */
  pariputs(ep->name); pariputc('(');
  narg = f->narg;
  for (i=1; i<=narg; i++, arg++)
  {
    print_var(*q++);
    print_def_arg(*arg);
    if (i == narg) { arg++; break; }
    pariputs(", ");
  }
  pariputs(") = ");
  narg = f->nloc;
  if (narg)
  {
    pariputs("local(");
    for (i=1; i<=narg; i++, arg++)
    {
      print_var(*q++);
      print_def_arg(*arg);
      if (i == narg) break;
      pariputs(", ");
    }
    pariputs("); ");
  }
  pariputs((char*)q);
}

void
print_user_member(entree *ep)
{
  GEN q = (GEN)ep->value;

  q++; /* skip initial NULL */
  print_var(*q++);
  pariprintf(".%s = ", ep->name);
  pariputs((char*)q);
}

static void
brace_print(entree *ep, void print(entree *))
{
  pariputc('{'); print(ep);
  pariputc('}'); pariputs("\n\n");
}

void
print_all_user_fun(void)
{
  entree *ep;
  int i;
  for (i = 0; i < functions_tblsz; i++)
    for (ep = functions_hash[i]; ep; ep = ep->next)
      if (EpVALENCE(ep) == EpUSER) brace_print(ep, &print_user_fun);
}

void
print_all_user_member(void)
{
  entree *ep;
  int i;
  for (i = 0; i < functions_tblsz; i++)
    for (ep = members_hash[i]; ep; ep = ep->next)
      if (EpVALENCE(ep) == EpMEMBER) brace_print(ep, &print_user_member);
}

