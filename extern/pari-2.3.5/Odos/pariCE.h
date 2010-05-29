/* only needed for Windows CE */
#ifndef PARICE_H
#define PARICE_H

extern int isspace(int);
extern int isdigit(int);
extern int isalpha(int);
extern int isalnum(int);
extern long strtol(const char *, char **, int );
extern int rename(const char *, const *);
extern int unlink(char *);
extern void *calloc(size_t, size_t);

#endif
