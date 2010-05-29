/* $Id: gp_init.c 7895 2006-04-19 20:34:12Z bill $

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
/*                    GP-SPECIFIC FUNCTIONS                        */
/*                                                                 */
/*******************************************************************/
#include "pari.h"
#include "paripriv.h"
#include "../graph/rect.h"
#include "../language/anal.h"
#include "gp.h"

static void whatnow0(char *s) { whatnow(s,0); }

#include "gp_init.h"

/* Backward Compatibility */

static GEN
gtype(GEN x)
{
  return stoi(typ(x));
}

static GEN
gsettype(GEN x,long t)
{
  x=gcopy(x); settyp(x,t); return x;
}

static long
setserieslength(long n)
{
  long m=precdl;
  if(n>0) precdl=n;
  return m;
}

static long
setprecr(long n)
{
  long m = GP_DATA->fmt->sigd;
  if (n > 0) { GP_DATA->fmt->sigd = n; precreal = ndec2prec(n); }
  return m;
}

entree functions_oldgp[] = {
{"allocatemem",11,(void *)allocatemoremem,2,"vLp","allocatemem(s)=allocates a new stack of s bytes, or doubles the stack if size is 0"},
{"box",35,(void *)rectbox,10,"vLGG","box(w,x2,y2)=if the cursor is at position (x1,y1), draw a box with diagonal (x1,y1) and (x2,y2) in rectwindow w (cursor does not move)"},
{"color",2,(void *)rectcolor,2,"vLL","color(w,c)=set default color to c in rectwindow. Possible values for c are 1=sienna, 2=cornsilk, 3=red, 4=black, 5=grey, 6=blue, 7=gainsborough"},
{"cursor",11,(void*)rectcursor,10,"vLp","cursor(w)=current position of cursor in rectwindow w"},
{"default",0,(void*)default0,11,"D\"\",r,D\"\",s,D0,L,","default({opt},{v},{flag}): set the default opt to v. If v is omitted, print the current default for opt. If no argument is given, print a list of all defaults as well as their values. If flag is non-zero, return the result instead of printing it on screen. See manual for details"},
{"draw",1,(void*)rectdraw,10,"vGp","draw(list)=draw vector of rectwindows list at indicated x,y positions; list is a vector w1,x1,y1,w2,x2,y2,etc..."},
{"initrect",34,(void*)initrect,10,"vLLL","initrect(w,x,y)=initialize rectwindow w to size x,y"},
{"kill",85,(void*)kill0,11,"vS","kill(x)=kills the present value of the variable or function x. Returns new value or 0"},
{"killrect",11,(void *)killrect,10,"vL","killrect(w)=erase the rectwindow w"},
{"line",35,(void *)rectline,10,"vLGG","line(w,x2,y2)=if cursor is at position (x1,y1), draw a line from (x1,y1) to (x2,y2) (and move the cursor) in the rectwindow w"},
{"lines",35,(void *)rectlines,10,"vLGG","lines(w,listx,listy)=draws an open polygon in rectwindow w where listx and listy contain the x (resp. y) coordinates of the vertices"},
{"move",35,(void*)rectmove,10,"vLGG","move(w,x,y)=move cursor to position x,y in rectwindow w"},
{"plot",99,(void *)plot,10,"vV=GGIDGDGp","plot(X=a,b,expr)=crude plot of expression expr, X goes from a to b"},
{"ploth",37,(void *)ploth,10,"V=GGIp","ploth(X=a,b,expr)=plot of expression expr, X goes from a to b in high resolution"},
{"ploth2",37,(void *)ploth2,10,"V=GGIp","ploth2(X=a,b,[expr1,expr2])=plot of points [expr1,expr2], X goes from a to b in high resolution"},
{"plothmult",37,(void *)plothmult,10,"V=GGIp","plothmult(X=a,b,[expr1,...])=plot of expressions expr1,..., X goes from a to b in high resolution"},
{"plothraw",2,(void *)plothraw,10,"GGp","plothraw(listx,listy)=plot in high resolution points whose x (resp. y) coordinates are in listx (resp. listy)"},
{"point",35,(void *)rectpoint,10,"vLGG","point(w,x,y)=draw a point (and move cursor) at position x,y in rectwindow w"},
{"points",35,(void *)rectpoints,10,"vLGG","points(w,listx,listy)=draws in rectwindow w the points whose x (resp y) coordinates are in listx (resp listy)"},
{"postdraw",1,(void *)postdraw,10,"vG","postdraw(list)=same as plotdraw, except that the output is a PostScript program in file \"pari.ps\""},
{"postploth",37,(void *)postploth,10,"V=GGIpD0,L,D0,L,","postploth(X=a,b,expr)=same as ploth, except that the output is a PostScript program in the file \"pari.ps\""},
{"postploth2",37,(void *)postploth2,10,"V=GGIpD0,L,","postploth2(X=a,b,[expr1,expr2])=same as ploth2, except that the output is a PostScript program in the file \"pari.ps\""},
{"postplothraw",2,(void *)postplothraw,10,"GGD0,L,","postplothraw(listx,listy)=same as plothraw, except that the output is a PostScript program in the file \"pari.ps\""},
{"pprint",0,(void*)printp,11,"vs*","pprint(a)=outputs a in beautified format ending with newline"},
{"pprint1",0,(void*)printp1,11,"vs*","pprint1(a)=outputs a in beautified format without ending with newline"},
{"print",0,(void*)print,11,"vs*","print(a)=outputs a in raw format ending with newline"},
{"print1",0,(void*)print1,11,"vs*","print1(a)=outputs a in raw format without ending with newline"},
{"rbox",35,(void *)rectrbox,10,"vLGG","rbox(w,dx,dy)=if the cursor is at (x1,y1), draw a box with diagonal (x1,y1)-(x1+dx,y1+dy) in rectwindow w (cursor does not move)"},
{"read",0,(void *)input0,11,"","read()=read an expression from the input file or standard input"},
{"rline",35,(void *)rectrline,10,"vLGG","rline(w,dx,dy)=if the cursor is at (x1,y1), draw a line from (x1,y1) to (x1+dx,y1+dy) (and move the cursor) in the rectwindow w"},
{"rlines",35,(void *)rectlines,10,"vLGG","rlines(w,dx,dy)=draw in rectwindow w the points given by vector of first coordinates xsand vector of second coordinates, connect them by lines"},
{"rmove",35,(void *)rectrmove,10,"vLGG","rmove(w,dx,dy)=move cursor to position (dx,dy) relative to the present position in the rectwindow w"},
{"rpoint",35,(void *)rectrpoint,10,"vLGG","rpoint(w,dx,dy)=draw a point (and move cursor) at position dx,dy relative to present position of the cursor in rectwindow w"},
{"rpoints",35,(void *)rectpoints,10,"vLGG","rpoints(w,xs,ys)=draw in rectwindow w the points given by vector of first coordinates xs and vector of second coordinates ys"},
{"scale",59,(void *)rectscale,10,"vLGGGG","scale(w,x1,x2,y1,y2)=scale the coordinates in rectwindow w so that x goes from x1 to x2 and y from y1 to y2 (y2<y1 is allowed)"},
{"setprecision",15,(void *)setprecr,2,"lL","setprecision(n)=set the current precision to n decimal digits if n>0, or return the current precision if n<=0"},
{"setserieslength",15,(void *)setserieslength,2,"lL","setserieslength(n)=set the default length of power series to n if n>0, or return the current default length if n<=0"},
{"settype",21,(void *)gsettype,2,"GL","settype(x,t)=make a copy of x with type t (to use with extreme care)"},
{"string",57,(void*)rectstring,10,"vLs","string(w,x)=draw in rectwindow w the string corresponding to x, where x is either a string, or a number in R, written in format 9.3"},
{"system",70,(void*) system0,11,"vs","system(a): a being a string, execute the system command a (not valid on every machine)"},
{"texprint",0,(void*)printtex,11,"vs*","texprint(a)=outputs a in TeX format"},
{"type",1,(void *)gtype,2,"Gp","type(x)=internal type number of the GEN x"},

{NULL,0,NULL,0,NULL,NULL} /* sentinel */
,};
