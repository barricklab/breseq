/* $Id: gp.h 7476 2005-11-26 16:31:45Z bill $

Copyright (C) 2000  The PARI group.

This file is part of the PARI/GP package.

PARI/GP is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation. It is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY WHATSOEVER.

Check the License for details. You should have received a copy of it, along
with the package; see the file 'COPYING'. If not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA. */

/*************************************************************************/
/*                                                                       */
/*                      GP-SPECIFIC DECLARATIONS                         */
/*                                                                       */
/*************************************************************************/
BEGINEXTERN

void aide(char *s, long flag);
int  get_line_from_readline(char *prompt, char *prompt_cont, filtre_t *F);
void gp_output(GEN z, gp_data *G);
void hit_return(void);
void init_readline(void);
void update_logfile(const char *prompt, const char *s);
void texmacs_completion(char *s, long pos);
char *color_prompt(char *prompt);
void print_fun_list(char **list, long nbli);

/* aide() */
#define h_REGULAR 0
#define h_LONG    1
#define h_APROPOS 2
#define h_RL      4

/* readline completions */
extern char *keyword_list[];

/* TeXmacs */
#define DATA_BEGIN  ((char) 2)
#define DATA_END    ((char) 5)
#define DATA_ESCAPE ((char) 27)

ENDEXTERN
