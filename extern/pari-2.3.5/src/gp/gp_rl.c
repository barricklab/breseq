/* $Id: gp_rl.c 12098 2010-01-29 13:40:22Z bill $

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
/*                 INTERFACE TO READLINE COMPLETION                */
/*                                                                 */
/*******************************************************************/
#include "pari.h"

#ifdef READLINE

#include "paripriv.h"
#include "../language/anal.h"
#include "gp.h"

typedef int (*RLCI)(int, int); /* rl_complete and rl_insert functions */
typedef char* (*GF)(const char*, int); /* generator function */

BEGINEXTERN
/***** Try to survive broken readline headers and obsolete versions *****/
#ifdef HAS_RL_MESSAGE
#  define USE_VARARGS
#  define PREFER_STDARG
#endif

#ifdef READLINE_LIBRARY
#  include <readline.h>
#  ifdef HAS_HISTORY_H
#    include <history.h>
#  endif
#else
#  include <readline/readline.h>
#  ifdef HAS_HISTORY_H
#    include <readline/history.h>
#  endif
#endif

#ifndef HAS_HISTORY_H
typedef struct HIST_ENTRY__ {char *line; void *data;} HIST_ENTRY;
extern HIST_ENTRY *history_get(int);
extern int history_length;
#endif

#ifdef HAS_RL_SAVE_PROMPT
#  define SAVE_PROMPT() rl_save_prompt()
#  define RESTORE_PROMPT() rl_restore_prompt()
#else
#  ifdef HAS_UNDERSCORE_RL_SAVE_PROMPT
extern void* _rl_restore_prompt(void), _rl_save_prompt(void);
#    define SAVE_PROMPT() _rl_save_prompt()
#    define RESTORE_PROMPT() _rl_restore_prompt()
#  else
#    define SAVE_PROMPT()
#    define RESTORE_PROMPT()
#  endif
#endif

#ifndef HAS_RL_MESSAGE
extern int rl_message (const char*, ...);
extern int rl_clear_message(), rl_begin_undo_group(), rl_end_undo_group();
extern int rl_read_key(), rl_stuff_char();
extern char *filename_completion_function(char *text,int state);
extern char *username_completion_function(char *text,int state);
#endif

extern RLCI rl_last_func;
ENDEXTERN

#ifdef HAS_RL_COMPLETION_MATCHES
#  define COMPLETION_MATCHES(a,b) rl_completion_matches((a),(b))
#  define FILE_COMPLETION rl_filename_completion_function
#  define USER_COMPLETION rl_username_completion_function
#  define DING rl_ding
#else
#  define COMPLETION_MATCHES(a,b) \
      (completion_matches((char*)(a),(CPFunction*)(b)))
#  define FILE_COMPLETION ((GF)filename_completion_function)
#  define USER_COMPLETION ((GF)username_completion_function)
#  define DING ding
#endif
/**************************************************************************/

static int pari_rl_back;
static int did_init_matched = 0;
static entree *current_ep = NULL;

static int
change_state(char *msg, ulong flag, int count)
{
  int c = (readline_state & flag) != 0;
  ulong o_readline_state = readline_state;

  switch(count)
  {
    default: c = 0; break; /* off */
    case -1: c = 1; break; /* on  */
    case -2: c = 1 - c; /* toggle */
  }
  if (c)
    readline_state |= flag;
  else {
    readline_state &= ~flag;
    if (!readline_state && o_readline_state)
	readline_state = 1;
  }
  SAVE_PROMPT();
#if 0				/* Does not work well... */
  rl_message("%s[%s: %s] %s", term_get_color(c_PROMPT),
	     msg, c? "on": "off", term_get_color(c_INPUT));
#else
  rl_message("[%s: %s] ", msg, c? "on": "off");
#endif
  c = rl_read_key();
  RESTORE_PROMPT();
  rl_clear_message();
  rl_stuff_char(c); return 1;
}

/* Wrapper around rl_complete to allow toggling insertion of arguments */
static int
pari_rl_complete(int count, int key)
{
  int ret;

  pari_rl_back = 0;
  if (count <= 0)
    return change_state("complete args", DO_ARGS_COMPLETE, count);

  rl_begin_undo_group();
  if (rl_last_func == pari_rl_complete)
    rl_last_func = (RLCI) rl_complete; /* Make repeated TABs different */
  ret = ((RLCI)rl_complete)(count,key);
  if (pari_rl_back && (pari_rl_back <= rl_point))
    rl_point -= pari_rl_back;
  rl_end_undo_group(); return ret;
}

static int did_matched_insert;

static int
pari_rl_matched_insert_suspend(int count, int key)
{
  ulong o_readline_state = readline_state;
  (void)count; (void)key;

  did_matched_insert = (readline_state & DO_MATCHED_INSERT);
  readline_state &= ~DO_MATCHED_INSERT;
  if (!readline_state && o_readline_state)
    readline_state = 1;
  return 1;
}

static int
pari_rl_matched_insert_restore(int count, int key)
{
  (void)count; (void)key;
  if (did_matched_insert)
    readline_state |= DO_MATCHED_INSERT;
  return 1;
}

static const char paropen[] = "([{";
static const char parclose[] = ")]}";

/* To allow insertion of () with a point in between. */
static int
pari_rl_matched_insert(int count, int key)
{
  int i = 0, ret;

  if (count <= 0)
    return change_state("electric parens", DO_MATCHED_INSERT, count);
  while (paropen[i] && paropen[i] != key) i++;
  if (!paropen[i] || !(readline_state & DO_MATCHED_INSERT) || GP_DATA->flags & EMACS)
    return ((RLCI)rl_insert)(count,key);
  rl_begin_undo_group();
  ((RLCI)rl_insert)(count,key);
  ret = ((RLCI)rl_insert)(count,parclose[i]);
  rl_point -= count;
  rl_end_undo_group(); return ret;
}

static int
pari_rl_default_matched_insert(int count, int key)
{
    if (!did_init_matched) {
      did_init_matched = 1;
      readline_state |= DO_MATCHED_INSERT;
    }
    return pari_rl_matched_insert(count, key);
}

static int
pari_rl_forward_sexp(int count, int key)
{
  int deep = 0, dir = 1, move_point = 0, lfail;

  (void)key;
  if (count < 0)
  {
    count = -count; dir = -1;
    if (!rl_point) goto fail;
    rl_point--;
  }
  while (count || deep)
  {
    move_point = 1;	/* Need to move point if moving left. */
    lfail = 0;		/* Do not need to fail left movement yet. */
    while ( !is_keyword_char(rl_line_buffer[rl_point])
            && !strchr("\"([{}])",rl_line_buffer[rl_point])
            && !( (dir == 1)
                  ? (rl_point >= rl_end)
                  : (rl_point <= 0 && (lfail = 1))))
        rl_point += dir;
    if (lfail || !rl_line_buffer[rl_point]) goto fail;

    if (is_keyword_char(rl_line_buffer[rl_point]))
    {
      while ( is_keyword_char(rl_line_buffer[rl_point])
              && (!((dir == 1) ? (rl_point >= rl_end) : (rl_point <= 0 && (lfail = 1)))
                  || (move_point = 0)))
        rl_point += dir;
      if (deep && lfail) goto fail;
      if (!deep) count--;
    }
    else if (strchr(paropen,rl_line_buffer[rl_point]))
    {
      if (deep == 0 && dir == -1) goto fail; /* We are already out of pars. */
      rl_point += dir;
      deep++; if (!deep) count--;
    }
    else if (strchr(parclose,rl_line_buffer[rl_point]))
    {
      if (deep == 0 && dir == 1)
      {
        rl_point++; goto fail; /* Get out of pars. */
      }
      rl_point += dir;
      deep--; if (!deep) count--;
    }
    else if (rl_line_buffer[rl_point] == '\"')
    {
      int bad = 1;

      rl_point += dir;
      while ( ((rl_line_buffer[rl_point] != '\"') || (bad = 0))
              && (!((dir == 1) ? (rl_point >= rl_end) : (rl_point <= 0))
                  || (move_point = 0)) )
        rl_point += dir;
      if (bad) goto fail;
      rl_point += dir;	/* Skip the other delimiter */
      if (!deep) count--;
    }
    else
    {
      fail: DING(); return 1;
    }
  }
  if (dir != 1 && move_point) rl_point++;
  return 1;
}

static int
pari_rl_backward_sexp(int count, int key)
{
  return pari_rl_forward_sexp(-count, key);
}

/* do we add () at the end of completed word? (is it a function?) */
static int
add_paren(int end)
{
  entree *ep;
  char *s;

  if (end < 0 || rl_line_buffer[end] == '(')
    return 0; /* not from command_generator or already there */
  ep = do_alias(current_ep); /* current_ep set in command_generator */
  if (EpVALENCE(ep) < EpUSER)
  { /* is it a constant masked as a function (e.g Pi)? */
    s = ep->help; if (!s) return 1;
    while (is_keyword_char(*s)) s++;
    return (*s != '=');
  }
  switch(EpVALENCE(ep))
  {
    case EpUSER:
    case EpINSTALL: return 1;
  }
  return 0;
}

static void
match_concat(char **matches, char *s)
{
  matches[0] = gprealloc(matches[0], strlen(matches[0])+strlen(s)+1);
  strcat(matches[0],s);
}

#define add_comma(x) (x==-2) /* from default_generator */

/* a single match, possibly modify matches[0] in place */
static void
treat_single(int code, char **matches)
{
  if (add_paren(code))
  {
    match_concat(matches,"()");
    pari_rl_back = 1;
    if (rl_point == rl_end)
#ifdef HAS_COMPLETION_APPEND_CHAR
      rl_completion_append_character = '\0'; /* Do not append space. */
#else
      pari_rl_back = 2;
#endif
  }
  else if (add_comma(code))
    match_concat(matches,",");
}
#undef add_comma


static char **
matches_for_emacs(const char *text, char **matches)
{
  if (!matches) printf("@");
  else
  {
    int i;
    printf("%s@", matches[0] + strlen(text));
    if (matches[1]) print_fun_list(matches+1,0);

   /* we don't want readline to do anything, but insert some junk
    * which will be erased by emacs.
    */
    for (i=0; matches[i]; i++) free(matches[i]);
    free(matches);
  }
  matches = (char **) gpmalloc(2*sizeof(char *));
  matches[0] = gpmalloc(2); sprintf(matches[0],"_");
  matches[1] = NULL; printf("@E_N_D"); pariflush();
  return matches;
}

/* Attempt to complete on the contents of TEXT. 'code' is used to
 * differentiate between callers when a single match is found.
 * Return the array of matches, NULL if there are none. */
static char **
get_matches(int code, const char *text, GF f)
{
  char **matches = COMPLETION_MATCHES(text, f);
  if (matches && !matches[1]) treat_single(code, matches);
  if (GP_DATA->flags & EMACS) matches = matches_for_emacs(text,matches);
  return matches;
}

#define DFLT 0
#define ENTREE 1
static char *
generator(void *list, const char *text, int *nn, int len, int typ)
{
  char *def = NULL, *name;
  int n = *nn;

  /* Return the next name which partially matches from list.*/
  switch(typ)
  {
    case DFLT :
      do
	def = (((default_type *) list)[n++]).name;
      while (def && strncmp(def,text,len));
      break;

    case ENTREE :
      do
	def = (((entree *) list)[n++]).name;
      while (def && strncmp(def,text,len));
  }

  *nn = n;
  if (def)
  {
    name = strcpy(gpmalloc(strlen(def)+1), def);
    return name;
  }
  return NULL; /* no names matched */
}
static char *
old_generator(const char *text,int state)
{
  static int n,len;
  static char *res;

  if (!state) { res = "a"; n=0; len=strlen(text); }
  if (res)
  {
    res = generator((void *)oldfonctions,text,&n,len,ENTREE);
    if (res) return res;
    n=0;
  }
  return generator((void *)functions_oldgp,text,&n,len,ENTREE);
}
static char *
default_generator(const char *text,int state)
{
  static int n,len;

  if (!state) { n=0; len=strlen(text); }
  return generator(gp_default_list,text,&n,len,DFLT);
}

static char *
add_prefix(char *name, const char *text, long junk)
{
  char *s = strncpy(gpmalloc(strlen(name)+1+junk),text,junk);
  strcpy(s+junk,name); return s;
}
static void
init_prefix(const char *text, int *len, int *junk, char **TEXT)
{
  long l = strlen(text), j = l-1;
  while (j >= 0 && is_keyword_char(text[j])) j--;
  j++;
  *TEXT = (char*)text + j;
  *junk = j;
  *len  = l - j;
}
/* Generator function for command completion.  STATE lets us know whether
 * to start from scratch; without any state (i.e. STATE == 0), then we
 * start at the top of the list. */
static char *
hashtable_generator(const char *text, int state, entree **hash)
{
  static int hashpos, len, junk;
  static entree* ep;
  static char *TEXT;

 /* If this is a new word to complete, initialize now:
  *  + indexes hashpos (GP hash list) and n (keywords specific to long help).
  *  + file completion and keyword completion use different word boundaries,
  *    have TEXT point to the keyword start.
  *  + save the length of TEXT for efficiency.
  */
  if (!state)
  {
    hashpos = 0; ep = hash[hashpos];
    init_prefix(text, &len, &junk, &TEXT);
  }

  /* Return the next name which partially matches from the command list. */
  for(;;)
    if (!ep)
    {
      if (++hashpos >= functions_tblsz) return NULL; /* no names matched */
      ep = hash[hashpos];
    }
    else if (strncmp(ep->name,TEXT,len))
      ep = ep->next;
    else
      break;
  current_ep = ep; ep = ep->next;
  return add_prefix(current_ep->name,text,junk);
}
static char *
command_generator(const char *text, int state)
{ return hashtable_generator(text,state, functions_hash); }
static char *
member_generator(const char *text, int state)
{ return hashtable_generator(text,state, members_hash); }

static char *
ext_help_generator(const char *text, int state)
{
  static int len, junk, n, def, key;
  static char *TEXT;
  if (!state) {
    n = 0;
    def = key = 1;
    init_prefix(text, &len, &junk, &TEXT);
  }
  if (def)
  {
    char *s = default_generator(TEXT, state);
    if (s) return add_prefix(s, text, junk);
    def = 0;
  }
  if (key)
  {
    for ( ; keyword_list[n]; n++)
      if (!strncmp(keyword_list[n],TEXT,len))
        return add_prefix(keyword_list[n++], text, junk);
    key = 0; state = 0;
  }
  return command_generator(text, state);
}

static void
rl_print_aide(char *s, int flag)
{
  int p = rl_point, e = rl_end;
  FILE *save = pari_outfile;

  rl_point = 0; rl_end = 0; pari_outfile = rl_outstream;
  SAVE_PROMPT();
  rl_message("%s",""); /* rl_message("") ==> "zero length format" warning */
  aide(s, flag);
  RESTORE_PROMPT();
  rl_point = p; rl_end = e; pari_outfile = save;
  rl_clear_message();
#ifdef RL_REFRESH_LINE_OLDPROTO
  rl_refresh_line();
#else
  rl_refresh_line(0,0);
#endif
}

/* add a space between \<char> and following text. Attempting completion now
 * would delete char. Hitting <TAB> again will complete properly */
static char **
add_space(int start)
{
  char **m;
  int p = rl_point + 1;
  rl_point = start + 2;
  rl_insert(1, ' '); rl_point = p;
#if 0 /* OK, but rings a bell */
#  ifdef HAS_RL_ATTEMPTED_COMPLETION_OVER
  rl_attempted_completion_over = 1;
#  endif
  return NULL;
#else /* better: fake an empty completion, but don't append ' ' after it! */
#  ifdef HAS_COMPLETION_APPEND_CHAR
  rl_completion_append_character = '\0';
#  endif
  m = (char**)gpmalloc(2 * sizeof(char*));
  m[0] = gpmalloc(1); *(m[0]) = 0;
  m[1] = NULL; return m;
#endif
}

char **
pari_completion(char *text, int START, int END)
{
  int i, first=0, start=START;

#ifdef HAS_COMPLETION_APPEND_CHAR
  rl_completion_append_character = ' ';
#endif
  current_ep = NULL;
/* If the line does not begin by a backslash, then it is:
 * . an old command ( if preceded by "whatnow(" ).
 * . a default ( if preceded by "default(" ).
 * . a member function ( if preceded by "." + keyword_chars )
 * . a file name (in current directory) ( if preceded by 'read' or 'writexx' )
 * . a command */
  if (start >=1 && rl_line_buffer[start] != '~') start--;
  while (start && is_keyword_char(rl_line_buffer[start])) start--;
  if (rl_line_buffer[start] == '~')
  {
    GF f = (GF)USER_COMPLETION;
    for(i=start+1;i<=END;i++)
      if (rl_line_buffer[i] == '/') { f = (GF)FILE_COMPLETION; break; }
    return get_matches(-1, text, f);
  }

  while (rl_line_buffer[first] && isspace((int)rl_line_buffer[first])) first++;
  switch (rl_line_buffer[first])
  {
    case '\\':
      if (first == start) return add_space(start);
      return get_matches(-1, text, FILE_COMPLETION);
    case '?':
      if (rl_line_buffer[first+1] == '?')
        return get_matches(-1, text, ext_help_generator);
      return get_matches(-1, text, command_generator);
  }

  while (start && rl_line_buffer[start] != '('
               && rl_line_buffer[start] != ',') start--;
  if (rl_line_buffer[start] == '(' && start)
  {
    int iend, j,k;
    entree *ep;
    char buf[200];

    i = start;

    while (i && isspace((int)rl_line_buffer[i-1])) i--;
    iend = i;
    while (i && is_keyword_char(rl_line_buffer[i-1])) i--;

    if (strncmp(rl_line_buffer + i,"default",7) == 0)
      return get_matches(-2, text, default_generator);
    if (strncmp(rl_line_buffer + i,"whatnow",7) == 0)
      return get_matches(-1, text, old_generator);
    if ( strncmp(rl_line_buffer + i,"read",4)  == 0
      || strncmp(rl_line_buffer + i,"write",5) == 0)
      return get_matches(-1, text, FILE_COMPLETION);

    j = start + 1;
    while (j <= END && isspace((int)rl_line_buffer[j])) j++;
    k = END;
    while (k > j && isspace((int)rl_line_buffer[k])) k--;
    /* If we are in empty parens, insert the default arguments */
    if ((readline_state & DO_ARGS_COMPLETE) && k == j
         && (rl_line_buffer[j] == ')' || !rl_line_buffer[j])
	 && (iend - i < (long)sizeof(buf))
	 && ( strncpy(buf, rl_line_buffer + i, iend - i),
	      buf[iend - i] = 0, 1)
	 && (ep = is_entry(buf)) && ep->help)
    {
#if 0 /* duplicate F1 */
      rl_print_aide(buf,h_RL);
#  ifdef HAS_RL_ATTEMPTED_COMPLETION_OVER
      rl_attempted_completion_over = 1;
#  endif
      return NULL;
#else
      char *s = ep->help;
      while (is_keyword_char(*s)) s++;
      if (*s++ == '(')
      { /* function call: insert arguments */
        char *e = s; 
        while (*e && *e != ')' && *e != '(') e++;
        if (*e == ')')
        { /* we just skipped over the arguments in short help text */
          char *str = strncpy(gpmalloc(e-s + 1), s, e-s);
          char **ret = (char**)gpmalloc(sizeof(char*)*2);
          str[e-s] = 0;
          ret[0] = str; ret[1] = NULL;
          if (GP_DATA->flags & EMACS) ret = matches_for_emacs("",ret);
          return ret;
        }
      }
#endif
    }
  }
  for(i = END-1; i >= start; i--)
    if (!is_keyword_char(rl_line_buffer[i]))
    {
      if (rl_line_buffer[i] == '.')
        return get_matches(-1, text, member_generator);
      break;
    }
  return get_matches(END, text, command_generator);
}

/* long help if count < 0 */
static int
rl_short_help(int count, int key)
{
  int flag = (count < 0 || rl_last_func == rl_short_help)? h_RL|h_LONG: h_RL;
  int off = rl_point;

  (void)key;
  /* func() with cursor on ')' following completion */
  if (off && rl_line_buffer[off-1] == '('
          && !is_keyword_char(rl_line_buffer[off])) off--;

  while (off && is_keyword_char(rl_line_buffer[off-1])) off--;

  /* Check for \c type situation.  Could check for leading whitespace too... */
  if (off == 1 && rl_line_buffer[off-1] == '\\') off--;
  if (off >= 8) { /* Check for default(whatever) */
    int t = off - 1;

    while (t >= 7 && isspace((int)rl_line_buffer[t])) t--;
    if (rl_line_buffer[t--] == '(') {
      while (t >= 6 && isspace((int)rl_line_buffer[t])) t--;
      t -= 6;
      if (t >= 0
	  && strncmp(rl_line_buffer + t, "default", 7) == 0
	  && (t == 0 || !is_keyword_char(rl_line_buffer[t-1])))
	off = t;
    }
  }
  rl_print_aide(rl_line_buffer + off, flag);
  return 0;
}

static int
rl_long_help(int count, int key)
{
  (void)count;
  return rl_short_help(-1,key);
}

void
init_readline(void)
{
  static int init_done = 0;

  if (init_done) return;
  init_done = 1;

  /* Allow conditional parsing of the ~/.inputrc file. */
  rl_readline_name = "Pari-GP";

  /* added ~, ? and , */
  rl_basic_word_break_characters = " \t\n\"\\'`@$><=;|&{(?,~";
  rl_special_prefixes = "~";

  /* custom completer */
#ifndef HAS_RL_COMPLETION_FUNC_T
# ifndef CPPFunction_defined
#   define CPPFunction Function
# endif
# define rl_completion_func_t CPPFunction
#endif
  rl_attempted_completion_function = (rl_completion_func_t*) pari_completion;

  /* we always want the whole list of completions under emacs */
#ifdef HAS_RL_COMPLETION_QUERY_ITEMS
  if (GP_DATA->flags & EMACS) rl_completion_query_items = 0x8fff;
#endif
#ifdef HAS_RL_BIND_KEY_IN_MAP
#  ifdef _RL_FUNCTION_TYPEDEF
#    define Bind(a,b,c) (rl_bind_key_in_map((a),(b),(c)))
#  else
#    define Bind(a,b,c) (rl_bind_key_in_map((a), (Function*)(b), (c)))
#  endif
#else
#  define Bind(a,b,c)
#endif

#ifdef _RL_FUNCTION_TYPEDEF
#  define Defun(a,b,c) (rl_add_defun((a), (b), (c)))
#else
#  define Defun(a,b,c) (rl_add_defun((char*)(a), (Function*)(b), (c)))
#endif

  Defun("short-help", rl_short_help, -1);
  Defun("long-help", rl_long_help, -1);
  Defun("pari-complete", pari_rl_complete, '\t');
  Defun("pari-matched-insert", pari_rl_default_matched_insert, -1);
  Defun("pari-matched-insert-suspend", pari_rl_matched_insert_suspend, -1);
  Defun("pari-matched-insert-restore", pari_rl_matched_insert_restore, -1);
  Defun("pari-forward-sexp", pari_rl_forward_sexp, -1);
  Defun("pari-backward-sexp", pari_rl_backward_sexp, -1);

  Bind('h', rl_short_help, emacs_meta_keymap);
  Bind('H', rl_long_help,  emacs_meta_keymap);
  Bind('h', rl_short_help, vi_movement_keymap);
  Bind('H', rl_long_help,  vi_movement_keymap);
#  ifdef HAS_RL_GENERIC_BIND
#define KSbind(s,f,k) rl_generic_bind(ISFUNC, (s), (char*)(f), (k))

  KSbind("OP",   rl_short_help,  emacs_meta_keymap); /* f1, vt100 */
  KSbind("[11~", rl_short_help,  emacs_meta_keymap); /* f1, xterm */
  KSbind("OP",   rl_short_help,  vi_movement_keymap); /* f1, vt100 */
  KSbind("[11~", rl_short_help,  vi_movement_keymap); /* f1, xterm */
  /* XTerm may signal start/end of paste by eming F200/F201 */
  /* XXXX For vi mode something more intelligent is needed - to switch to the
     insert mode - and back when restoring. */
  KSbind("[200~", pari_rl_matched_insert_suspend,  emacs_meta_keymap);  /* pre-paste xterm */
  KSbind("[200~", pari_rl_matched_insert_suspend,  vi_movement_keymap); /* pre-paste xterm */
  KSbind("[201~", pari_rl_matched_insert_restore,  emacs_meta_keymap);  /* post-paste xterm */
  KSbind("[201~", pari_rl_matched_insert_restore,  vi_movement_keymap); /* post-paste xterm */
#  endif
  Bind('(', pari_rl_matched_insert, emacs_standard_keymap);
  Bind('[', pari_rl_matched_insert, emacs_standard_keymap);
  Bind(6, pari_rl_forward_sexp,  emacs_meta_keymap); /* M-C-f */
  Bind(2, pari_rl_backward_sexp, emacs_meta_keymap); /* M-C-b */

#ifdef EMACS_DOS_KEYMAP
  Bind(';', rl_short_help, emacs_dos_keymap); /* F1 */
  Bind('T', rl_long_help,  emacs_dos_keymap); /* Shift-F1 */
  Bind(155, pari_rl_backward_sexp, emacs_dos_keymap); /* Alt-Left */
  Bind(157, pari_rl_forward_sexp,  emacs_dos_keymap); /* Alt-Right*/
#endif
}

static void
print_escape_string(char *s)
{
  long l = strlen(s);
  char *t, *t0 = gpmalloc(l * 3 + 3);

  t = t0; *t++ = '"';
  for ( ;*s; *t++ = *s++)
    switch(*s)
    {
      case DATA_BEGIN:
      case DATA_END:
      case DATA_ESCAPE: *t++ = DATA_ESCAPE; continue;

      case '\\':
      case '"': *t++ = '\\'; continue;
    }
  *t++ = '"';
  *t = '\0'; puts(t0); free(t0);
}

static char *
completion_word(long end)
{
  char *s = rl_line_buffer + end, *found_quote = NULL;
  long i;
  /* truncate at cursor position */
  *s = 0;
  /* first look for unclosed string */
  for (i=0; i < end; i++)
  {
    switch(rl_line_buffer[i])
    {
      case '"':
        found_quote = found_quote? NULL: rl_line_buffer + i;
        break;

      case '\\': i++; break;
    }

  }
  if (found_quote) return found_quote + 1; /* return next char after quote */

  /* else find beginning of word */
  while (s >  rl_line_buffer)
  {
    s--;
    if (!is_keyword_char(*s)) { s++; break; }
  }
  return s;
}

/* completion required, cursor on s + pos. Complete wrt strict left prefix */
void
texmacs_completion(char *s, long pos)
{
  char **matches, *text;

  if (rl_line_buffer) free(rl_line_buffer);
  rl_line_buffer = pari_strdup(s);
  text = completion_word(pos);
  /* text = start of expression we complete */
  rl_end = strlen(s)-1;
  rl_point = pos;
  matches = pari_completion(text, text - rl_line_buffer, pos);
  printf("%cscheme:(tuple",DATA_BEGIN);
  if (matches)
  {
    long i, prelen = (rl_line_buffer+pos) - text;
    char *t = gpmalloc(prelen+1);
    strncpy(t, text, prelen); t[prelen] = 0; /* prefix */
    printf(" ");
    print_escape_string(t); free(t);
    for (i = matches[1]? 1: 0; matches[i]; i++)
    {
      printf(" ");
      print_escape_string(matches[i] + prelen);
      free(matches[i]);
    }
    free(matches);
  }
  printf(")%c", DATA_END);
  fflush(stdout);
}

static int
history_is_new(char *s)
{
  HIST_ENTRY *e;
  if (!*s) return 0;
  if (!history_length) return 1;
  e = history_get(history_length);
  /* paranoia: e != NULL, unless readline is in a weird state */
  return e? strcmp(s, e->line): 0;
}

static void
gp_add_history(char *s)
{
  if (history_is_new(s)) add_history(s);
}

/* Read line; returns a malloc()ed string of the user input or NULL on EOF.
   Increments the buffer size appropriately if needed; fix *endp if so. */
static char *
gprl_input(char **endp, int first, input_method *IM, filtre_t *F)
{
  Buffer *b = F->buf;
  ulong used = *endp - b->buf;
  ulong left = b->len - used, l;
  char *s, *t, *prompt = expand_prompt(first? IM->prompt: IM->prompt_cont, F);

  if (! (s = readline(color_prompt(prompt))) ) return NULL; /* EOF */
  gp_add_history(s); /* Makes a copy */
  l = strlen(s) + 1;
  /* put back \n that readline stripped. This is needed for
   * { print("a
   *   b"); }
   * and conforms with the other input methods anyway. */
  t = gpmalloc(l + 1);
  strncpy(t, s, l-1);
  t[l-1] = '\n';
  t[l]   = 0; /* equivalent to sprintf(t,"%s\n", s) */
  if (left < l)
  {
    ulong incr = b->len;
    if (incr < l) incr = l;
    fix_buffer(b, b->len + incr);
    *endp = b->buf + used;
  }
  return t;
}

/* request one line interactively.
 * Return 0: EOF
 *        1: got one line from readline or infile */
int
get_line_from_readline(char *prompt, char *prompt_cont, filtre_t *F)
{
  const int index = history_length;
  char *s;
  input_method IM;

  IM.prompt      = prompt;
  IM.prompt_cont = prompt_cont;
  IM.getline = &gprl_input;
  IM.free = 1;
  if (! input_loop(F,&IM)) { pariputs("\n"); return 0; }

  s = F->buf->buf;
  if (*s)
  {
    if (history_length > index+1)
    { /* Multi-line input. Remove incomplete lines */
      int i = history_length;
      while (i > index) {
        HIST_ENTRY *e = remove_history(--i);
        free(e->line); free(e);
      }
      gp_add_history(s);
    }

    if (logfile) update_logfile(expand_prompt(IM.prompt, F), s);
  }
  return 1;
}
#endif
