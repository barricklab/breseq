/* $Id: plotport.c 7601 2006-01-11 11:35:06Z kb $

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
/*                         PLOT ROUTINES                           */
/*                                                                 */
/*******************************************************************/
#include "pari.h"
#include "paripriv.h"
#include "rect.h"
#include "../language/anal.h"

void postdraw0(long *w, long *x, long *y, long lw);
void postdraw00(long *w, long *x, long *y, long lw, long scale);
static void PARI_get_psplot(void);

static long current_color[NUMRECT];
PariRect **rectgraph = NULL;
PARI_plot pari_plot, pari_psplot;
PARI_plot *pari_plot_engine = &pari_plot;
long  rectpoint_itype = 0;
long  rectline_itype  = 0;

#define STRINGRECT (NUMRECT-2)
#define DRAWRECT (NUMRECT-1)

#define PLOTH_NUMPOINTS 1000
#define PARAM_NUMPOINTS 1500
#define RECUR_NUMPOINTS 8

#define RECUR_MAXDEPTH 10
#define RECUR_PREC 0.001
#define PARAMR_MAXDEPTH 10

/********************************************************************/
/**                                                                **/
/**                         LOW-RES PLOT                           **/
/**                                                                **/
/********************************************************************/
#define ISCR 64
#define JSCR 22
#define BLANK ' '
#define ZERO1 ','
#define ZERO2 '-'
#define ZERO3 '`'
#define PICTZERO(j) ((j) % 3 ? ((j) % 3 == 2 ? ZERO3 : ZERO2) : ZERO1)
#define YY '|'
#define XX_UPPER '\''
#define XX_LOWER '.'
#define FF1 '_'
#define FF2 'x'
#define FF3 '"'
#define PICT(j) ((j) % 3 ? ((j) % 3 == 2 ? FF3 : FF2) : FF1)

static char *
dsprintf9(double d, char *buf)
{
  int i = 10;

  while (--i >= 0) {
    sprintf(buf, "%9.*g", i, d);
    if (strlen(buf) <= 9) return buf;
  }
  return buf; /* Should not happen? */
}

typedef unsigned char screen[ISCR+1][JSCR+1];

static void
fill_gap(screen scr, long i, int jnew, int jpre)
{
  int mid, i_up, i_lo, up, lo;

  if (jpre < jnew - 2) {
    up = jnew - 1; i_up = i;
    lo = jpre + 1; i_lo = i - 1;
  } else if (jnew < jpre - 2) {
    up = jpre - 1; i_up = i - 1;
    lo = jnew + 1; i_lo = i;
  } else return; /* if gap < 2, leave it as it is. */

  mid = (jpre+jnew)/2;
  if (mid>JSCR) mid=JSCR; else if (mid<0) mid=0;
  if (lo<0) lo=0;
  if (lo<=JSCR) while (lo <= mid) scr[i_lo][lo++] = ':';
  if (up>JSCR) up=JSCR;
  if (up>=0) while (up > mid) scr[i_up][up--] = ':';
}

static double
todbl(GEN x) { return rtodbl(gtofp(x, 3)); }

#define QUARK  ((char*)NULL) /* Used as a special-case */
static GEN quark_gen;

static GEN
READ_EXPR(char *s, entree *ep, GEN x) {
  if (s == QUARK) return gsubst(quark_gen,0,x);
  ep->value = x; return readseq(s);
}

void
plot(entree *ep, GEN a, GEN b, char *ch,GEN ysmlu,GEN ybigu, long prec)
{
  long jz, j, i, sig;
  pari_sp av = avma, av2, limite;
  int jnew, jpre = 0; /* for lint */
  GEN ysml, ybig, x, diff, dyj, dx, y[ISCR+1];
  screen scr;
  char buf[80], z;

  sig=gcmp(b,a); if (!sig) return;
  if (sig<0) { x=a; a=b; b=x; }
  x = gtofp(a, prec); push_val(ep, x);
  for (i=1; i<=ISCR; i++) y[i]=cgetr(3);
  dx = gtofp(gdivgs(gsub(b,a), ISCR-1), prec);
  ysml=gen_0; ybig=gen_0;
  for (j=1; j<=JSCR; j++) scr[1][j]=scr[ISCR][j]=YY;
  for (i=2; i<ISCR; i++)
  {
    scr[i][1]   = XX_LOWER;
    scr[i][JSCR]= XX_UPPER;
    for (j=2; j<JSCR; j++) scr[i][j] = BLANK;
  }
  av2=avma; limite=stack_lim(av2,1);
  for (i=1; i<=ISCR; i++)
  {
    gaffect(READ_EXPR(ch,ep,x), y[i]);
    if (gcmp(y[i],ysml)<0) ysml=y[i];
    if (gcmp(y[i],ybig)>0) ybig=y[i];
    x = addrr(x,dx);
    if (low_stack(limite, stack_lim(av2,1)))
    {
      pari_sp tetpil=avma;
      if (DEBUGMEM>1) pari_warn(warnmem,"plot");
      x = gerepile(av2,tetpil,rcopy(x));
    }
  }
  if (ysmlu) ysml=ysmlu;
  if (ybigu) ybig=ybigu;
  avma=av2; diff=gsub(ybig,ysml);
  if (gcmp0(diff)) { ybig=gaddsg(1,ybig); diff=gen_1; }
  dyj = gdivsg((JSCR-1)*3+2,diff);
  jz = 3-gtolong(gmul(ysml,dyj));
  av2=avma; z = PICTZERO(jz); jz = jz/3;
  for (i=1; i<=ISCR; i++)
  {
    if (0<=jz && jz<=JSCR) scr[i][jz]=z;
    j = 3+gtolong(gmul(gsub(y[i],ysml),dyj));
    jnew = j/3;
    if (i > 1) fill_gap(scr, i, jnew, jpre);
    if (0<=jnew && jnew<=JSCR) scr[i][jnew] = PICT(j);
    avma = av2;
    jpre = jnew;
  }
  pariputc('\n');
  pariprintf("%s ", dsprintf9(todbl(ybig), buf));
  for (i=1; i<=ISCR; i++) pariputc(scr[i][JSCR]);
  pariputc('\n');
  for (j=(JSCR-1); j>=2; j--)
  {
    pariputs("          ");
    for (i=1; i<=ISCR; i++) pariputc(scr[i][j]);
    pariputc('\n');
  }
  pariprintf("%s ", dsprintf9(todbl(ysml), buf));
  for (i=1; i<=ISCR; i++)  pariputc(scr[i][1]);
  pariputc('\n');
  pariprintf("%10s%-9.7g%*.7g\n"," ",todbl(a),ISCR-9,todbl(b));
  pop_val(ep); avma=av;
}

/********************************************************************/
/**                                                                **/
/**                      RECTPLOT FUNCTIONS                        **/
/**                                                                **/
/********************************************************************/
void
init_graph(void)
{
  int n;

  rectgraph = (PariRect**) gpmalloc(sizeof(PariRect*)*NUMRECT);
  for (n=0; n<NUMRECT; n++)
  {
    PariRect *e = (PariRect*) gpmalloc(sizeof(PariRect));

    e->head = e->tail = NULL;
    e->sizex = e->sizey = 0;
    current_color[n] = DEFAULT_COLOR;
    rectgraph[n] = e;
  }
}

void
free_graph(void)
{
  int i;

  if (!rectgraph)
      return;
  for (i=0; i<NUMRECT; i++)
  {
    PariRect *e=rectgraph[i];

    if (RHead(e)) killrect(i);
    free((void *)e);
  }
  free((void *)rectgraph);
  rectgraph = 0;
}

static PariRect *
check_rect(long ne)
{
  if (!GOODRECT(ne))
    pari_err(talker,
        "incorrect rectwindow number in graphic function (%ld not in [0, %ld])",
        ne, NUMRECT-1);
  return rectgraph[ne];
}

static PariRect *
check_rect_init(long ne)
{
  PariRect *e = check_rect(ne);
  if (!RHead(e)) pari_err(talker,"you must initialize the rectwindow first");
  return e;
}

void
initrect_gen(long ne, GEN x, GEN y, long flag)
{
  long xi, yi;
  if (flag) {
    double xd = gtodouble(x), yd = gtodouble(y);

    PARI_get_plot(0);
    xi = pari_plot.width - 1;
    yi = pari_plot.height - 1;
    if (xd) xi = DTOL(xd*xi);
    if (yd) yi = DTOL(yd*yi);
  } else {
    xi = itos(x);
    yi = itos(y);
    if (!xi || !yi) PARI_get_plot(0);
    if (!xi) xi = pari_plot.width - 1;
    if (!yi) yi = pari_plot.height - 1;
  }
  initrect(ne, xi, yi);
}

void
initrect(long ne, long x, long y)
{
  PariRect *e;
  RectObj *z;

  if (x<=1 || y<=1) pari_err(talker,"incorrect dimensions in initrect");
  e = check_rect(ne); if (RHead(e)) killrect(ne);

  z = (RectObj*) gpmalloc(sizeof(RectObj));
  RoNext(z) = NULL;
  RoType(z) = ROt_NULL;
  RHead(e)=RTail(e)=z;
  RXsize(e)=x; RYsize(e)=y;
  RXcursor(e)=0; RYcursor(e)=0;
  RXscale(e)=1; RXshift(e)=0;
  RYscale(e)=1; RYshift(e)=0;
  RHasGraph(e) = 0;
}

GEN
rectcursor(long ne)
{
  PariRect *e = check_rect_init(ne);
  return mkvec2s((long)RXcursor(e), (long)RYcursor(e));
}

static void
rectscale0(long ne, double x1, double x2, double y1, double y2)
{
  PariRect *e = check_rect_init(ne);
  double p2,p3;

  p2 = RXshift(e) + RXscale(e) * RXcursor(e);
  p3 = RYshift(e) + RYscale(e) * RYcursor(e);
  RXscale(e) = RXsize(e)/(x2-x1); RXshift(e) = -x1*RXscale(e);
  RYscale(e) = RYsize(e)/(y1-y2); RYshift(e) = -y2*RYscale(e);
  RXcursor(e) = (p2 - RXshift(e)) / RXscale(e);
  RYcursor(e) = (p3 - RYshift(e)) / RYscale(e);
}

void
rectscale(long ne, GEN x1, GEN x2, GEN y1, GEN y2)
{
  rectscale0(ne, gtodouble(x1), gtodouble(x2), gtodouble(y1), gtodouble(y2));
}

static void
rectmove0(long ne, double x, double y, long relative)
{
  PariRect *e = check_rect_init(ne);
  RectObj *z = (RectObj*) gpmalloc(sizeof(RectObj1P));

  if (relative)
   { RXcursor(e) += x; RYcursor(e) += y; }
  else
   { RXcursor(e) = x; RYcursor(e) = y; }
  RoNext(z) = 0; RoType(z) = ROt_MV;
  RoMVx(z) = RXcursor(e) * RXscale(e) + RXshift(e);
  RoMVy(z) = RYcursor(e) * RYscale(e) + RYshift(e);
  if (!RHead(e)) RHead(e)=RTail(e)=z;
  else { RoNext(RTail(e))=z; RTail(e)=z; }
}

void
rectmove(long ne, GEN x, GEN y)
{
  rectmove0(ne,gtodouble(x),gtodouble(y),0);
}

void
rectrmove(long ne, GEN x, GEN y)
{
  rectmove0(ne,gtodouble(x),gtodouble(y),1);
}

void
rectpoint0(long ne, double x, double y,long relative) /* code = ROt_MV/ROt_PT */
{
  PariRect *e = check_rect_init(ne);
  RectObj *z = (RectObj*) gpmalloc(sizeof(RectObj1P));

  if (relative)
   { RXcursor(e) += x; RYcursor(e) += y; }
  else
   { RXcursor(e) = x; RYcursor(e) = y; }
  RoNext(z)=0;
  RoPTx(z) = RXcursor(e)*RXscale(e) + RXshift(e);
  RoPTy(z) = RYcursor(e)*RYscale(e) + RYshift(e);
  RoType(z) = ( DTOL(RoPTx(z)) < 0
		|| DTOL(RoPTy(z)) < 0 || DTOL(RoPTx(z)) > RXsize(e)
		|| DTOL(RoPTy(z)) > RYsize(e) ) ? ROt_MV : ROt_PT;
  if (!RHead(e)) RHead(e)=RTail(e)=z;
  else { RoNext(RTail(e))=z; RTail(e)=z; }
  RoCol(z)=current_color[ne];
}

void
rectpoint(long ne, GEN x, GEN y)
{
  rectpoint0(ne,gtodouble(x),gtodouble(y),0);
}

void
rectrpoint(long ne, GEN x, GEN y)
{
  rectpoint0(ne,gtodouble(x),gtodouble(y),1);
}

void
rectcolor(long ne, long color)
{
  check_rect(ne);
  if (!GOODCOLOR(color)) pari_err(talker,"This is not a valid color");
  current_color[ne]=color;
}

void
rectline0(long ne, double gx2, double gy2, long relative) /* code = ROt_MV/ROt_LN */
{
  double dx,dy,dxy,xmin,xmax,ymin,ymax,x1,y1,x2,y2;
  PariRect *e = check_rect_init(ne);
  RectObj *z = (RectObj*) gpmalloc(sizeof(RectObj2P));
  const double c = 1 + 1e-10;

  x1 = RXcursor(e)*RXscale(e) + RXshift(e);
  y1 = RYcursor(e)*RYscale(e) + RYshift(e);
  if (relative)
    { RXcursor(e)+=gx2; RYcursor(e)+=gy2; }
  else
    { RXcursor(e)=gx2; RYcursor(e)=gy2; }
  x2 = RXcursor(e)*RXscale(e) + RXshift(e);
  y2 = RYcursor(e)*RYscale(e) + RYshift(e);
  xmin = max(min(x1,x2),0); xmax = min(max(x1,x2),RXsize(e));
  ymin = max(min(y1,y2),0); ymax = min(max(y1,y2),RYsize(e));
  dxy = x1*y2 - y1*x2; dx = x2-x1; dy = y2-y1;
  if (dy)
  {
    if (dx*dy<0)
      { xmin = max(xmin,(dxy+RYsize(e)*dx)/dy); xmax=min(xmax,dxy/dy); }
    else
      { xmin=max(xmin,dxy/dy); xmax=min(xmax,(dxy+RYsize(e)*dx)/dy); }
  }
  if (dx)
  {
    if (dx*dy<0)
      { ymin=max(ymin,(RXsize(e)*dy-dxy)/dx); ymax=min(ymax,-dxy/dx); }
    else
      { ymin=max(ymin,-dxy/dx); ymax=min(ymax,(RXsize(e)*dy-dxy)/dx); }
  }
  RoNext(z)=0;
  RoLNx1(z) = xmin; RoLNx2(z) = xmax;
  if (dx*dy<0) { RoLNy1(z) = ymax; RoLNy2(z) = ymin; }
  else { RoLNy1(z) = ymin; RoLNy2(z) = ymax; }
  RoType(z) = (xmin>xmax*c || ymin>ymax*c) ? ROt_MV : ROt_LN;
  if (!RHead(e)) RHead(e)=RTail(e)=z;
  else { RoNext(RTail(e))=z; RTail(e)=z; }
  RoCol(z)=current_color[ne];
}

/* Given coordinates of ends of a line, and labels l1 l2 attached to the
   ends, plot ticks where the label coordinate takes "round" values */

static void
rectticks(PARI_plot *WW, long ne,
          double dx1, double dy1, double dx2, double dy2,
          double l1, double l2, long flags)
{
  long dx,dy,dxy,dxy1,x1,y1,x2,y2,nticks,n,n1,dn;
  double minstep, maxstep, step, l_min, l_max, minl, maxl, dl, dtx, dty, x, y;
  double ddx, ddy;
  const double mult[3] = { 2./1., 5./2., 10./5. };
  PariRect *e = check_rect_init(ne);
  int do_double = !(flags & TICKS_NODOUBLE);

  x1 = DTOL(dx1*RXscale(e) + RXshift(e));
  y1 = DTOL(dy1*RYscale(e) + RYshift(e));
  x2 = DTOL(dx2*RXscale(e) + RXshift(e));
  y2 = DTOL(dy2*RYscale(e) + RYshift(e));
  dx = x2 - x1;
  dy = y2 - y1;
  if (dx < 0) dx = -dx;
  if (dy < 0) dy = -dy;
  dxy1 = max(dx, dy);
  dx /= WW->hunit;
  dy /= WW->vunit;
  dxy = (long)sqrt(dx*dx + dy*dy);
  nticks = (long) ((dxy + 2.5)/4);
  if (!nticks) return;

  /* Now we want to find nticks (or less) "round" numbers between l1 and l2.
     For our purpose round numbers have "last significant" digit either
	*) any;
	*) even;
	*) divisible by 5.
     We need to choose which alternative is better.
   */
  if (l1 < l2)
    l_min = l1, l_max = l2;
  else
    l_min = l2, l_max = l1;
  minstep = (l_max - l_min)/(nticks + 1);
  maxstep = 2.5*(l_max - l_min);
  step = exp(log(10) * floor(log10(minstep)));
  if (!(flags & TICKS_ENDSTOO)) {
    double d = 2*(l_max - l_min)/dxy1;	/* Two pixels off */

    l_min += d;
    l_max -= d;
  }
  for (n = 0; ; n++) {
    if (step >= maxstep) return;

    if (step >= minstep) {
      minl = ceil(l_min/step);
      maxl = floor(l_max/step);
      if (minl <= maxl && maxl - minl + 1 <= nticks) {
	nticks = (long) (maxl - minl + 1);
        l_min = minl * step;
        l_max = maxl * step; break;
      }
    }
    step *= mult[ n % 3 ];
  }
  /* Where to position doubleticks, variants:
     small: each 5, double: each 10	(n===2 mod 3)
     small: each 2, double: each 10	(n===1  mod 3)
     small: each 1, double: each  5 */
  dn = (n % 3 == 2)? 2: 5;
  n1 = ((long)minl) % dn; /* unused if do_double = FALSE */

  /* now l_min and l_max keep min/max values of l with ticks, and nticks is
     the number of ticks to draw. */
  if (nticks == 1) ddx = ddy = 0; /* unused: for lint */
  else {
    dl = (l_max - l_min)/(nticks - 1);
    ddx = (dx2 - dx1) * dl / (l2 - l1);
    ddy = (dy2 - dy1) * dl / (l2 - l1);
  }
  x = dx1 + (dx2 - dx1) * (l_min - l1) / (l2 - l1);
  y = dy1 + (dy2 - dy1) * (l_min - l1) / (l2 - l1);
  /* assume hunit and vunit form a square.  For clockwise ticks: */
  dtx = WW->hunit * dy/dxy * (y2 > y1 ? 1 : -1);	/* y-coord runs down */
  dty = WW->vunit * dx/dxy * (x2 > x1 ? 1 : -1);
  for (n = 0; n < nticks; n++) {
    RectObj *z = (RectObj*) gpmalloc(sizeof(RectObj2P));
    double lunit = WW->hunit > 1 ? 1.5 : 2;
    double l = (do_double && (n + n1) % dn == 0) ? lunit: 1;

    RoNext(z) = 0;
    RoLNx1(z) = RoLNx2(z) = x*RXscale(e) + RXshift(e);
    RoLNy1(z) = RoLNy2(z) = y*RYscale(e) + RYshift(e);

    if (flags & TICKS_CLOCKW) {
      RoLNx1(z) += dtx*l;
      RoLNy1(z) -= dty*l;		/* y-coord runs down */
    }
    if (flags & TICKS_ACLOCKW) {
      RoLNx2(z) -= dtx*l;
      RoLNy2(z) += dty*l;		/* y-coord runs down */
    }
    RoType(z) = ROt_LN;

    if (!RHead(e))
      RHead(e)=RTail(e)=z;
    else {
      RoNext(RTail(e))=z; RTail(e)=z;
    }
    RoCol(z)=current_color[ne];
    x += ddx;
    y += ddy;
  }
}

void
rectline(long ne, GEN gx2, GEN gy2)
{
  rectline0(ne, gtodouble(gx2), gtodouble(gy2),0);
}

void
rectrline(long ne, GEN gx2, GEN gy2)
{
  rectline0(ne, gtodouble(gx2), gtodouble(gy2),1);
}

void
rectbox0(long ne, double gx2, double gy2, long relative)
{
  double x1,y1,x2,y2,xmin,ymin,xmax,ymax;
  double xx,yy;
  PariRect *e = check_rect_init(ne);
  RectObj *z = (RectObj*) gpmalloc(sizeof(RectObj2P));

  x1 = RXcursor(e)*RXscale(e) + RXshift(e);
  y1 = RYcursor(e)*RYscale(e) + RYshift(e);
  if (relative)
  { xx = RXcursor(e)+gx2; yy = RYcursor(e)+gy2; }
  else
  {  xx = gx2; yy = gy2; }
  x2 = xx*RXscale(e) + RXshift(e);
  y2 = yy*RYscale(e) + RYshift(e);
  xmin = max(min(x1,x2),0); xmax = min(max(x1,x2),RXsize(e));
  ymin = max(min(y1,y2),0); ymax = min(max(y1,y2),RYsize(e));

  RoNext(z)=0; RoType(z) = ROt_BX;
  RoBXx1(z) = xmin; RoBXy1(z) = ymin;
  RoBXx2(z) = xmax; RoBXy2(z) = ymax;
  if (!RHead(e)) RHead(e)=RTail(e)=z;
  else { RoNext(RTail(e))=z; RTail(e)=z; }
  RoCol(z)=current_color[ne];
}

void
rectbox(long ne, GEN gx2, GEN gy2)
{
  rectbox0(ne, gtodouble(gx2), gtodouble(gy2), 0);
}

void
rectrbox(long ne, GEN gx2, GEN gy2)
{
  rectbox0(ne, gtodouble(gx2), gtodouble(gy2), 1);
}

void
killrect(long ne)
{
  RectObj *p1,*p2;
  PariRect *e = check_rect_init(ne);

  current_color[ne]=DEFAULT_COLOR;
  p1=RHead(e);
  RHead(e) = RTail(e) = NULL;
  RXsize(e) = RYsize(e) = 0;
  RXcursor(e) = RYcursor(e) = 0;
  RXscale(e) = RYscale(e) = 1;
  RXshift(e) = RYshift(e) = 0;
  while (p1)
  {
    if (RoType(p1)==ROt_MP || RoType(p1)==ROt_ML)
    {
      free(RoMPxs(p1)); free(RoMPys(p1));
    }
    if (RoType(p1)==ROt_ST) free(RoSTs(p1));
    p2=RoNext(p1); free(p1); p1=p2;
  }
}

void
rectpoints0(long ne, double *listx, double *listy, long lx) /* code = ROt_MP */
{
  double *ptx, *pty, x, y;
  long i, cp=0;
  PariRect *e = check_rect_init(ne);
  RectObj *z = (RectObj*) gpmalloc(sizeof(RectObjMP));

  ptx=(double*) gpmalloc(lx*sizeof(double));
  pty=(double*) gpmalloc(lx*sizeof(double));
  for (i=0; i<lx; i++)
  {
    x = RXscale(e)*listx[i] + RXshift(e);
    y = RYscale(e)*listy[i] + RYshift(e);
    if ((x>=0)&&(y>=0)&&(x<=RXsize(e))&&(y<=RYsize(e)))
    {
      ptx[cp]=x; pty[cp]=y; cp++;
    }
  }
  RoNext(z)=0; RoType(z) = ROt_MP;
  RoMPcnt(z) = cp; RoMPxs(z) = ptx; RoMPys(z) = pty;
  if (!RHead(e)) RHead(e)=RTail(e)=z;
  else { RoNext(RTail(e))=z; RTail(e)=z; }
  RoCol(z)=current_color[ne];
}

void
rectpoints(long ne, GEN listx, GEN listy)
{
  long i,lx, tx=typ(listx), ty=typ(listy);
  double *px,*py;

  if (!is_matvec_t(tx) || !is_matvec_t(ty)) {
    rectpoint(ne, listx, listy); return;
  }
  lx = lg(listx);
  if (tx == t_MAT || ty == t_MAT || lg(listy) != lx) pari_err(typeer,"rectpoints");
  lx--; if (!lx) return;

  px = (double*) gpmalloc(lx*sizeof(double));
  py = (double*) gpmalloc(lx*sizeof(double));
  for (i=0; i<lx; i++)
  {
    px[i]=gtodouble(gel(listx,i+1)); py[i]=gtodouble(gel(listy,i+1));
  }
  rectpoints0(ne,px,py,lx);
  free(px); free(py);
}

void
rectlines0(long ne, double *x, double *y, long lx, long flag) /* code = ROt_ML */
{
  long i,I;
  double *ptx,*pty;
  PariRect *e = check_rect_init(ne);
  RectObj *z = (RectObj*) gpmalloc(sizeof(RectObj2P));

  I = flag ? lx+1 : lx;
  ptx = (double*) gpmalloc(I*sizeof(double));
  pty = (double*) gpmalloc(I*sizeof(double));
  for (i=0; i<lx; i++)
  {
    ptx[i] = RXscale(e)*x[i] + RXshift(e);
    pty[i] = RYscale(e)*y[i] + RYshift(e);
  }
  if (flag)
  {
    ptx[i] = RXscale(e)*x[0] + RXshift(e);
    pty[i] = RYscale(e)*y[0] + RYshift(e);
  }
  RoNext(z)=0; RoType(z)=ROt_ML;
  RoMLcnt(z)=lx; RoMLxs(z)=ptx; RoMLys(z)=pty;
  if (!RHead(e)) RHead(e)=RTail(e)=z;
  else { RoNext(RTail(e))=z; RTail(e)=z; }
  RoCol(z) = current_color[ne];
}

void
rectlines(long ne, GEN listx, GEN listy, long flag)
{
  long tx=typ(listx), ty=typ(listy), lx=lg(listx), i;
  double *x, *y;

  if (!is_matvec_t(tx) || !is_matvec_t(ty))
  {
    rectline(ne, listx, listy); return;
  }
  if (tx == t_MAT || ty == t_MAT || lg(listy) != lx) pari_err(typeer,"rectlines");
  lx--; if (!lx) return;

  x = (double*) gpmalloc(lx*sizeof(double));
  y = (double*) gpmalloc(lx*sizeof(double));
  for (i=0; i<lx; i++)
  {
    x[i] = gtodouble(gel(listx,i+1));
    y[i] = gtodouble(gel(listy,i+1));
  }
  rectlines0(ne,x,y,lx,flag);
  free(x); free(y);
}

static void
put_string(long win, long x, long y, char *str, long dir)
{
  rectmove0(win,(double)x,(double)y,0); rectstring3(win,str,dir);
}

void
rectstring(long ne, char *str)
{
  rectstring3(ne,str,RoSTdirLEFT);
}

/* Allocate memory, then put string */
void
rectstring3(long ne, char *str, long dir) /* code = ROt_ST */
{
  PariRect *e = check_rect_init(ne);
  RectObj *z = (RectObj*) gpmalloc(sizeof(RectObjST));
  long l = strlen(str);
  char *s = (char *) gpmalloc(l+1);

  strcpy(s,str);
  RoNext(z)=0; RoType(z) = ROt_ST;
  RoSTl(z) = l; RoSTs(z) = s;
  RoSTx(z) = RXscale(e)*RXcursor(e)+RXshift(e);
  RoSTy(z) = RYscale(e)*RYcursor(e)+RYshift(e);
  RoSTdir(z) = dir;
  if (!RHead(e)) RHead(e)=RTail(e)=z;
  else { RoNext(RTail(e))=z; RTail(e)=z; }
  RoCol(z)=current_color[ne];
}

void
rectpointtype(long ne, long type) /* code = ROt_PTT */
{
 if (ne == -1) {
     rectpoint_itype = type;
 } else {
  PariRect *e = check_rect_init(ne);
  RectObj *z = (RectObj*) gpmalloc(sizeof(RectObjPN));

  RoNext(z) = 0; RoType(z) = ROt_PTT;
  RoPTTpen(z) = type;
  if (!RHead(e)) RHead(e)=RTail(e)=z;
  else { RoNext(RTail(e))=z; RTail(e)=z; }
 }
}

/*FIXME: this function is a noop, since no graphic driver implement
 * the ROt_PTS code. 
 * ne==-1 is a legacy, meningless value.
 */
void
rectpointsize(long ne, GEN size) /* code = ROt_PTS */
{
 if (ne == -1) {
     /*do nothing*/
 } else {
     PariRect *e = check_rect_init(ne);
     RectObj *z = (RectObj*) gpmalloc(sizeof(RectObjPS));

     RoNext(z) = 0; RoType(z) = ROt_PTS;
     RoPTSsize(z) = gtodouble(size);
     if (!RHead(e)) RHead(e)=RTail(e)=z;
     else { RoNext(RTail(e))=z; RTail(e)=z; }
 }
}

void
rectlinetype(long ne, long type)
{
 if (ne == -1) {
     rectline_itype = type;
 } else {
  PariRect *e = check_rect_init(ne);
  RectObj *z = (RectObj*) gpmalloc(sizeof(RectObjPN));

  RoNext(z) = 0; RoType(z) = ROt_LNT;
  RoLNTpen(z) = type;
  if (!RHead(e)) RHead(e)=RTail(e)=z;
  else { RoNext(RTail(e))=z; RTail(e)=z; }
 }
}

void
rectcopy_gen(long source, long dest, GEN xoff, GEN yoff, long flag)
{
  long xi, yi;
  if (flag & RECT_CP_RELATIVE) {
    double xd = gtodouble(xoff), yd = gtodouble(yoff);

    PARI_get_plot(0);
    xi = pari_plot.width - 1;
    yi = pari_plot.height - 1;
    xi = DTOL(xd*xi);
    yi = DTOL(yd*yi);
  } else {
    xi = itos(xoff);  yi = itos(yoff);
  }
  if (flag & ~RECT_CP_RELATIVE) {
    PariRect *s = check_rect_init(source), *d = check_rect_init(dest);

    switch (flag & ~RECT_CP_RELATIVE) {
      case RECT_CP_NW:
	break;
      case RECT_CP_SW:
	yi = RYsize(d) - RYsize(s) - yi;
	break;
      case RECT_CP_SE:
	yi = RYsize(d) - RYsize(s) - yi;
	/* FALL THROUGH */
      case RECT_CP_NE:
	xi = RXsize(d) - RXsize(s) - xi;
	break;
    }
  }
  rectcopy(source, dest, xi, yi);
}

void
rectcopy(long source, long dest, long xoff, long yoff)
{
  PariRect *s = check_rect_init(source), *d = check_rect_init(dest);
  RectObj *R = RHead(s);
  RectObj *tail = RTail(d);
  RectObj *next;
  int i;

  while (R)
  {
    switch(RoType(R))
    {
      case ROt_PT:
	next = (RectObj*) gpmalloc(sizeof(RectObj1P));
	memcpy(next,R,sizeof(RectObj1P));
	RoPTx(next) += xoff; RoPTy(next) += yoff;
	RoNext(tail) = next; tail = next;
	break;
      case ROt_LN: case ROt_BX:
	next = (RectObj*) gpmalloc(sizeof(RectObj2P));
	memcpy(next,R,sizeof(RectObj2P));
	RoLNx1(next) += xoff; RoLNy1(next) += yoff;
	RoLNx2(next) += xoff; RoLNy2(next) += yoff;
	RoNext(tail) = next; tail = next;
	break;
      case ROt_MP: case ROt_ML:
	next = (RectObj*) gpmalloc(sizeof(RectObjMP));
	memcpy(next,R,sizeof(RectObjMP));
	RoMPxs(next) = (double*) gpmalloc(sizeof(double)*RoMPcnt(next));
	RoMPys(next) = (double*) gpmalloc(sizeof(double)*RoMPcnt(next));
	memcpy(RoMPxs(next),RoMPxs(R),sizeof(double)*RoMPcnt(next));
	memcpy(RoMPys(next),RoMPys(R),sizeof(double)*RoMPcnt(next));
	for (i=0; i<RoMPcnt(next); i++)
	{
	  RoMPxs(next)[i] += xoff; RoMPys(next)[i] += yoff;
	 }
	RoNext(tail) = next; tail = next;
	break;
      case ROt_ST:
	next = (RectObj*) gpmalloc(sizeof(RectObjST));
	memcpy(next,R,sizeof(RectObjMP));
	RoSTs(next) = (char*) gpmalloc(RoSTl(R)+1);
	memcpy(RoSTs(next),RoSTs(R),RoSTl(R)+1);
	RoSTx(next) += xoff; RoSTy(next) += yoff;
	RoNext(tail) = next; tail = next;
	break;
      case ROt_PTT: case ROt_LNT: case ROt_PTS:
	next = (RectObj*) gpmalloc(sizeof(RectObjPN));
	memcpy(next,R,sizeof(RectObjPN));
	RoNext(tail) = next; tail = next;
	break;
    }
    R=RoNext(R);
  }
  RoNext(tail) = NULL; RTail(d) = tail;
}

#define CLIPLINE_NONEMPTY	1
#define CLIPLINE_CLIP_1		2
#define CLIPLINE_CLIP_2		4

/* A simpler way is to clip by 4 half-planes */
static int
clipline(double xmin, double xmax, double ymin, double ymax,
         double *x1p, double *y1p, double *x2p, double *y2p)
{
    int xy_exch = 0, rc = CLIPLINE_NONEMPTY;
    double t, sl;
    double xi, xmn, xmx;
    double yi, ymn, ymx;
    int x1_is_ymn, x1_is_xmn;
    double x1 = *x1p, x2 = *x2p, y1 = *y1p, y2 = *y2p;

    if ((x1 < xmin &&  x2 < xmin) || (x1 > xmax && x2 > xmax))
	return 0;
    if (fabs(x1 - x2) < fabs(y1 - y2)) { /* Exchange x and y */
	xy_exch = 1;
	t = xmin; xmin = ymin; ymin = t;
	t = xmax; xmax = ymax; ymax = t;
	t = x1; x1 = y1; y1 = t;
	t = x2; x2 = y2; y2 = t;
    }

    /* Build y as a function of x */
    xi = x1, yi = y1;
    sl = ( (x1==x2) ? 0 : (y2 - yi)/(x2 - xi) );

    xmn = x1, xmx = x2;
    if (x1 > x2) {
	x1_is_xmn = 0;
	xmn = x2, xmx = x1;
    } else
	x1_is_xmn = 1;

    if (xmn < xmin) {
	xmn = xmin;
	rc |= (x1_is_xmn ? CLIPLINE_CLIP_1 : CLIPLINE_CLIP_2);
    }

    if (xmx > xmax) {
	xmx = xmax;
	rc |= (x1_is_xmn ? CLIPLINE_CLIP_2 : CLIPLINE_CLIP_1);
    }

    if (xmn > xmx)
	return 0;

    ymn = yi + (xmn - xi)*sl;
    ymx = yi + (xmx - xi)*sl;

    if (sl < 0)
	t = ymn, ymn = ymx, ymx = t;
    if (ymn > ymax || ymx < ymin)
	return 0;

    if (rc & CLIPLINE_CLIP_1)
	x1 = (x1_is_xmn ? xmn : xmx);
    if (rc & CLIPLINE_CLIP_2)
	x2 = (x1_is_xmn ? xmx : xmn);

    /* Now we know there is an intersection, need to move x1 and x2 */
    x1_is_ymn = ((sl >= 0) == (x1 < x2));
    if (ymn < ymin) {
	double x = (ymin - yi)/sl + xi; /* slope != 0  ! */

	if (x1_is_ymn)
	    x1 = x, rc |= CLIPLINE_CLIP_1;
	else
	    x2 = x, rc |= CLIPLINE_CLIP_2;
    }
    if (ymx > ymax) {
	double x = (ymax - yi)/sl + xi; /* slope != 0  ! */

	if (x1_is_ymn)
	    x2 = x, rc |= CLIPLINE_CLIP_2;
	else
	    x1 = x, rc |= CLIPLINE_CLIP_1;
    }
    if (rc & CLIPLINE_CLIP_1)
	y1 = yi + (x1 - xi)*sl;
    if (rc & CLIPLINE_CLIP_2)
	y2 = yi + (x2 - xi)*sl;
    if (xy_exch)			/* Exchange x and y */
	*x1p = y1, *x2p = y2, *y1p = x1, *y2p = x2;
    else
	*x1p = x1, *x2p = x2, *y1p = y1, *y2p = y2;
    return rc;
}

void
rectclip(long rect)
{
  PariRect *s = check_rect_init(rect);
  RectObj *R = RHead(s);
  RectObj **prevp = &RHead(s);
  RectObj *next;
  double xmin = 0;
  double xmax = RXsize(s);
  double ymin = 0;
  double ymax = RYsize(s);

  while (R) {
      int did_clip = 0;

      next = RoNext(R);
      switch(RoType(R)) {
      case ROt_PT:
	  if ( DTOL(RoPTx(R)) < xmin || DTOL(RoPTx(R)) > xmax
	       || DTOL(RoPTy(R)) < ymin || DTOL(RoPTy(R)) > ymax) {
		 remove:
	      *prevp = next;
	      free(R);
	      break;
	  }
	  goto do_next;
      case ROt_BX:
	  if (RoLNx1(R) < xmin)
	      RoLNx1(R) = xmin, did_clip = 1;
	  if (RoLNx2(R) < xmin)
	      RoLNx2(R) = xmin, did_clip = 1;
	  if (RoLNy1(R) < ymin)
	      RoLNy1(R) = ymin, did_clip = 1;
	  if (RoLNy2(R) < ymin)
	      RoLNy2(R) = ymin, did_clip = 1;
	  if (RoLNx1(R) > xmax)
	      RoLNx1(R) = xmax, did_clip = 1;
	  if (RoLNx2(R) > xmax)
	      RoLNx2(R) = xmax, did_clip = 1;
	  if (RoLNy1(R) > ymax)
	      RoLNy1(R) = ymax, did_clip = 1;
	  if (RoLNy2(R) > ymax)
	      RoLNy2(R) = ymax, did_clip = 1;
	  /* Remove zero-size clipped boxes */
	  if ( did_clip
	       && RoLNx1(R) == RoLNx2(R) && RoLNy1(R) == RoLNy2(R) )
	      goto remove;
	  goto do_next;
      case ROt_LN:
	  if (!clipline(xmin, xmax, ymin, ymax,
			&RoLNx1(R), &RoLNy1(R), &RoLNx2(R), &RoLNy2(R)))
	      goto remove;
	  goto do_next;
      case ROt_MP:
      {
	  int c = RoMPcnt(R);
	  int f = 0, t = 0;

	  while (f < c) {
	      if ( DTOL(RoMPxs(R)[f]) >= xmin && DTOL(RoMPxs(R)[f]) <= xmax
		   && DTOL(RoMPys(R)[f]) >= ymin
		   && DTOL(RoMPys(R)[f]) <= ymax) {
		  if (t != f) {
		      RoMPxs(R)[t] = RoMPxs(R)[f];
		      RoMPys(R)[t] = RoMPys(R)[f];
		  }
		  t++;
	      }
	      f++;
	  }
	  if (t == 0)
	      goto remove;
	  RoMPcnt(R) = t;
	  goto do_next;
      }
      case ROt_ML:
      {
	  /* Hard case.  Here we need to break a multiline into
	     several pieces if some part is clipped. */
	  int c = RoMPcnt(R) - 1;
	  int f = 0, t = 0, had_lines = 0, had_hole = 0, rc;
	  double ox = RoMLxs(R)[0], oy = RoMLys(R)[0], oxn, oyn;

	  while (f < c) {
	      /* Endpoint of this segment is the startpoint of the
		 next one, so we need to preserve it if it is clipped. */
	      oxn = RoMLxs(R)[f + 1], oyn = RoMLys(R)[f + 1];
	      rc = clipline(xmin, xmax, ymin, ymax,
			    /* &RoMLxs(R)[f], &RoMLys(R)[f], */
			    &ox, &oy,
			    &RoMLxs(R)[f+1], &RoMLys(R)[f+1]);
	      RoMLxs(R)[f] = ox; RoMLys(R)[f] = oy;
	      ox = oxn; oy = oyn;
	      if (!rc) {
		  if (had_lines)
		      had_hole = 1;
		  f++;
		  continue;
	      }

	      if ( !had_lines || (!(rc & CLIPLINE_CLIP_1) && !had_hole) ) {
		  /* Continuous */
		  had_lines = 1;
		  if (t != f) {
		      if (t == 0) {
			  RoMPxs(R)[t] = RoMPxs(R)[f];
			  RoMPys(R)[t] = RoMPys(R)[f];
		      }
		      RoMPxs(R)[t+1] = RoMPxs(R)[f+1];
		      RoMPys(R)[t+1] = RoMPys(R)[f+1];
		  }
		  t++;
		  f++;
		  if ( rc & CLIPLINE_CLIP_2)
		      had_hole = 1, RoMLcnt(R) = t + 1;
		  continue;
	      }
	      /* Is not continuous, automatically R is not free()ed.  */
	      RoMLcnt(R) = t + 1;
	      if ( rc & CLIPLINE_CLIP_2) { /* Needs separate entry */
		  RectObj *n = (RectObj*) gpmalloc(sizeof(RectObj2P));

		  RoType(n) = ROt_LN;
		  RoCol(n) = RoCol(R);
		  RoLNx1(n) = RoMLxs(R)[f];	RoLNy1(n) = RoMLys(R)[f];
		  RoLNx2(n) = RoMLxs(R)[f+1];	RoLNy2(n) = RoMLys(R)[f+1];
		  RoNext(n) = next;
		  RoNext(R) = n;
		  /* Restore the unclipped value: */
		  RoMLxs(R)[f + 1] = oxn;	RoMLys(R)[f + 1] = oyn;
		  f++;
		  prevp = &RoNext(n);
	      }
	      if (f + 1 < c) {		/* Are other lines */
		  RectObj *n = (RectObj*) gpmalloc(sizeof(RectObjMP));
		  RoType(n) = ROt_ML;
		  RoCol(n) = RoCol(R);
		  RoMLcnt(n) = c - f;
		  RoMLxs(n) = (double*) gpmalloc(sizeof(double)*(c - f));
		  RoMLys(n) = (double*) gpmalloc(sizeof(double)*(c - f));
		  memcpy(RoMPxs(n),RoMPxs(R) + f, sizeof(double)*(c - f));
		  memcpy(RoMPys(n),RoMPys(R) + f, sizeof(double)*(c - f));
		  RoMPxs(n)[0] = oxn; RoMPys(n)[0] = oyn;
		  RoNext(n) = next;
		  RoNext(R) = n;
		  next = n;
	      }
	      goto do_next;
	  }
	  if (t == 0)
	      goto remove;
	  goto do_next;
      }
      default: {
	do_next:
	      prevp = &RoNext(R);
	      break;
	  }
      }
      R = next;
  }
}

/********************************************************************/
/**                                                                **/
/**                        HI-RES PLOT                             **/
/**                                                                **/
/********************************************************************/

static void
Appendx(dblPointList *f, dblPointList *l,double x)
{
  (l->d)[l->nb++]=x;
  if (x < f->xsml) f->xsml=x;
  else if (x > f->xbig) f->xbig=x;
}

static void
Appendy(dblPointList *f, dblPointList *l,double y)
{
  (l->d)[l->nb++]=y;
  if (y < f->ysml) f->ysml=y;
  else if (y > f->ybig) f->ybig=y;
}

/* Convert data from GEN to double before we call rectplothrawin */
static dblPointList*
gtodblList(GEN data, long flags)
{
  dblPointList *l;
  double xsml,xbig,ysml,ybig;
  long tx=typ(data), ty, nl=lg(data)-1, lx1,lx;
  register long i,j,u,v;
  long param=(flags & PLOT_PARAMETRIC);
  GEN x,y;

  if (! is_vec_t(tx)) pari_err(typeer,"gtodblList");
  if (!nl) return NULL;
  lx1 = lg(data[1]);

  if (nl == 1) pari_err(talker,"single vector in gtodblList");
  /* Allocate memory, then convert coord. to double */
  l = (dblPointList*) gpmalloc(nl*sizeof(dblPointList));
  for (i=0; i<nl-1; i+=2)
  {
    u = i+1;
    x = gel(data,u);   tx = typ(x); lx = lg(x);
    y = gel(data,u+1); ty = typ(y);
    if (!is_vec_t(tx) || !is_vec_t(ty) || lg(y) != lx
        || (!param && lx != lx1)) pari_err(typeer,"gtodblList");

    lx--;
    l[i].d = (double*) gpmalloc(lx*sizeof(double));
    l[u].d = (double*) gpmalloc(lx*sizeof(double));
    for (j=0; j<lx; j=v)
    {
      v = j+1;
      l[i].d[j] = gtodouble(gel(x,v));
      l[u].d[j] = gtodouble(gel(y,v));
    }
    l[i].nb = l[u].nb = lx;
  }

  /* Now compute extremas */
  if (param)
  {
    l[0].nb = nl/2;
    for (i=0; i < l[0].nb; i+=2)
      if (l[i+1].nb) break;
    if (i >= l[0].nb) { free(l); return NULL; }
    xsml = xbig = l[i  ].d[0];
    ysml = ybig = l[i+1].d[0];

    for (i=0; i < l[0].nb; i+=2)
    {
      u = i+1;
      for (j=0; j < l[u].nb; j++)
      {
	if      (l[i].d[j] < xsml) xsml = l[i].d[j];
	else if (l[i].d[j] > xbig) xbig = l[i].d[j];

	if      (l[u].d[j] < ysml) ysml = l[u].d[j];
	else if (l[u].d[j] > ybig) ybig = l[u].d[j];
      }
    }
  }
  else
  {
    if (!l[0].nb) { free(l); return NULL; }
    l[0].nb = nl-1;

    xsml = xbig = l[0].d[0];
    ysml = ybig = l[1].d[0];

    for (j=0; j < l[1].nb; j++)
    {
      if      (l[0].d[j] < xsml) xsml = l[0].d[j];
      else if (l[0].d[j] > xbig) xbig = l[0].d[j];
    }
    for (i=1; i <= l[0].nb; i++)
      for (j=0; j < l[i].nb; j++)
      {
	if      (l[i].d[j] < ysml) ysml = l[i].d[j];
	else if (l[i].d[j] > ybig) ybig = l[i].d[j];
      }
  }
  l[0].xsml = xsml; l[0].xbig = xbig;
  l[0].ysml = ysml; l[0].ybig = ybig;
  return l;
}

static void
single_recursion(dblPointList *pl,char *ch,entree *ep,GEN xleft,GEN yleft,
  GEN xright,GEN yright,long depth)
{
  GEN xx,yy;
  pari_sp av=avma;
  double dy=pl[0].ybig - pl[0].ysml;

  if (depth==RECUR_MAXDEPTH) return;

  xx = gmul2n(gadd(xleft,xright),-1);
  yy = gtofp(READ_EXPR(ch,ep,xx), 3);

  if (dy)
  {
    if (fabs(rtodbl(yleft)+rtodbl(yright)-2*rtodbl(yy))/dy < RECUR_PREC)
      return;
  }
  single_recursion(pl,ch,ep, xleft,yleft, xx,yy, depth+1);

  Appendx(&pl[0],&pl[0],rtodbl(xx));
  Appendy(&pl[0],&pl[1],rtodbl(yy));

  single_recursion(pl,ch,ep, xx,yy, xright,yright, depth+1);

  avma=av;
}

static void
param_recursion(dblPointList *pl,char *ch,entree *ep, GEN tleft,GEN xleft,
  GEN yleft, GEN tright,GEN xright,GEN yright, long depth)
{
  GEN tt,xx,yy, p1;
  pari_sp av=avma;
  double dy=pl[0].ybig - pl[0].ysml;
  double dx=pl[0].xbig - pl[0].xsml;

  if (depth==PARAMR_MAXDEPTH) return;

  tt=gmul2n(gadd(tleft,tright),-1);
  p1 = READ_EXPR(ch,ep,tt);
  xx = gtofp(gel(p1,1), 3);
  yy = gtofp(gel(p1,2), 3);

  if (dx && dy)
  {
    if (fabs(rtodbl(xleft)+rtodbl(xright)-2*rtodbl(xx))/dx < RECUR_PREC &&
       fabs(rtodbl(yleft)+rtodbl(yright)-2*rtodbl(yy))/dy < RECUR_PREC)
        return;
  }
  param_recursion(pl,ch,ep, tleft,xleft,yleft, tt,xx,yy, depth+1);

  Appendx(&pl[0],&pl[0],rtodbl(xx));
  Appendy(&pl[0],&pl[1],rtodbl(yy));

  param_recursion(pl,ch,ep, tt,xx,yy, tright,xright,yright, depth+1);

  avma=av;
}

/*
 *  Pure graphing. If testpoints is 0, it is set to the default.
 *  Returns a dblPointList of (absolute) coordinates.
 */
static dblPointList *
rectplothin(entree *ep, GEN a, GEN b, char *ch, long prec, ulong flags,
            long testpoints)
{
  long single_c;
  long param=flags & PLOT_PARAMETRIC;
  long recur=flags & PLOT_RECURSIVE;
  GEN p1,dx,x,xleft,xright,yleft,yright,tleft,tright;
  dblPointList *pl;
  long tx, i, j, sig, nc, nl, nbpoints;
  pari_sp av = avma, av2;
  double xsml,xbig,ysml,ybig,fx,fy;

  if (!testpoints)
  {
    if (recur)
      testpoints = RECUR_NUMPOINTS;
    else
      testpoints = param? PARAM_NUMPOINTS : PLOTH_NUMPOINTS;
  }
  if (recur)
    nbpoints = testpoints << (param? PARAMR_MAXDEPTH : RECUR_MAXDEPTH);
  else
    nbpoints = testpoints;

  sig=gcmp(b,a); if (!sig) return 0;
  if (sig<0) swap(a, b);
  dx=gdivgs(gsub(b,a),testpoints-1);

  x = gtofp(a, prec); push_val(ep, x);
  av2=avma; p1=READ_EXPR(ch,ep,x); tx=typ(p1);
  if (!is_matvec_t(tx))
  {
    xsml = gtodouble(a);
    xbig = gtodouble(b);
    ysml = ybig = gtodouble(p1);
    nc=1; nl=2; /* nc = nb of curves; nl = nb of coord. lists */
    if (param) pari_warn(warner,"flag PLOT_PARAMETRIC ignored");
    single_c=1; param=0;
  }
  else
  {
    if (tx != t_VEC) pari_err(talker,"not a row vector in ploth");
    nl=lg(p1)-1; if (!nl) { avma=av; return 0; }
    single_c=0;
    if (param) nc=nl/2; else { nc=nl; nl++; }
    if (recur && nc > 1) pari_err(talker,"multi-curves cannot be plot recursively");

    if (param)
    {
      xbig=xsml=gtodouble(gel(p1,1));
      ybig=ysml=gtodouble(gel(p1,2));
      for (i=3; i<=nl; i++)
      {
	fx=gtodouble(gel(p1,i)); i++;
        fy=gtodouble(gel(p1,i));
	if (xsml>fx) xsml=fx; else if (xbig<fx) xbig=fx;
	if (ysml>fy) ysml=fy; else if (ybig<fy) ybig=fy;
      }
    }
    else
    {
      xsml=gtodouble(a); xbig=gtodouble(b);
      ysml=gtodouble(vecmin(p1)); ybig=gtodouble(vecmax(p1));
    }
  }

  pl=(dblPointList*) gpmalloc(nl*sizeof(dblPointList));
  for (i = 0; i < nl; i++)
  {
    pl[i].d = (double*) gpmalloc((nbpoints+1)*sizeof(dblPointList));
    pl[i].nb=0;
  }
  pl[0].xsml=xsml; pl[0].ysml=ysml; pl[0].xbig=xbig; pl[0].ybig=ybig;

  if (recur) /* recursive plot */
  {
    xleft=cgetr(3); xright=cgetr(3); yleft=cgetr(3); yright=cgetr(3);
    if (param)
    {
      tleft=cgetr(prec); tright=cgetr(prec);
      av2=avma;
      gaffect(a,tleft); p1=READ_EXPR(ch,ep,tleft);
      gaffect(gel(p1,1),xleft);
      gaffect(gel(p1,2),yleft);
      for (i=0; i<testpoints-1; i++)
      {
	if (i) {
	  gaffect(tright,tleft); gaffect(xright,xleft); gaffect(yright,yleft);
	 }
	gaddz(tleft,dx,tright);
        p1 = READ_EXPR(ch,ep,tright);
        if (lg(p1) != 3) pari_err(talker,"inconsistent data in rectplothin");
        gaffect(gel(p1,1),xright); gaffect(gel(p1,2),yright);

	Appendx(&pl[0],&pl[0],rtodbl(xleft));
	Appendy(&pl[0],&pl[1],rtodbl(yleft));

	param_recursion(pl,ch,ep, tleft,xleft,yleft, tright,xright,yright, 0);
	avma=av2;
      }
      Appendx(&pl[0],&pl[0],rtodbl(xright));
      Appendy(&pl[0],&pl[1],rtodbl(yright));
    }
    else /* single_c */
    {
      av2=avma;
      gaffect(a,xleft);
      gaffect(READ_EXPR(ch,ep,xleft), yleft);
      for (i=0; i<testpoints-1; i++)
      {
        gaddz(xleft,dx,xright);
	gaffect(READ_EXPR(ch,ep,xright),yright);

	Appendx(&pl[0],&pl[0],rtodbl(xleft));
	Appendy(&pl[0],&pl[1],rtodbl(yleft));

        single_recursion(pl,ch,ep,xleft,yleft,xright,yright,0);
        avma=av2;
        gaffect(xright,xleft); gaffect(yright,yleft);
      }
      Appendx(&pl[0],&pl[0],rtodbl(xright));
      Appendy(&pl[0],&pl[1],rtodbl(yright));
    }
  }
  else /* non-recursive plot */
  {
    if (single_c)
      for (i=0; i<testpoints; i++)
      {
	p1 = READ_EXPR(ch,ep,x);
	pl[0].d[i]=gtodouble(x);
	Appendy(&pl[0],&pl[1],gtodouble(p1));
	gaddz(x,dx,x); avma=av2;
      }
    else if (param)
    {
      long k;
      double z;

      for (i=0; i<testpoints; i++)
      {
	p1 = READ_EXPR(ch,ep,x);
        if (lg(p1) != nl+1) pari_err(talker,"inconsistent data in rectplothin");
	for (j=0; j<nl; j=k)
	{
	  k=j+1; z=gtodouble(gel(p1,k));
	  Appendx(&pl[0], &pl[j],z);

	  j=k; k++; z=gtodouble(gel(p1,k));
	  Appendy(&pl[0], &pl[j],z);
	 }
	gaddz(x,dx,x); avma=av2;
      }
    }
    else /* plothmult */
      for (i=0; i<testpoints; i++)
      {
	p1 = READ_EXPR(ch,ep,x);
        if (lg(p1) != nl) pari_err(talker,"inconsistent data in rectplothin");
	pl[0].d[i]=gtodouble(x);
	for (j=1; j<nl; j++) { Appendy(&pl[0],&pl[j],gtodouble(gel(p1,j))); }
	gaddz(x,dx,x); avma=av2;
      }
  }
  pl[0].nb=nc; pop_val(ep); avma = av;
  return pl;
}

/* Uses highlevel plotting functions to implement splines as
   a low-level plotting function.
   Most probably one does not need to make varn==0 into pure variable,
   since plotting functions should take care of this. */
static void
rectsplines(long ne, double *x, double *y, long lx, long flag)
{
  long i, j;
  pari_sp oldavma = avma;
  GEN tas, xa = cgetg(lx+1, t_VEC), ya = cgetg(lx+1, t_VEC);
  entree *var0 = varentries[ordvar[0]];

  if (lx < 4) pari_err(talker, "Too few points (%ld) for spline plot", lx);
  for (i = 1; i <= lx; i++) {
    gel(xa,i) = dbltor(x[i-1]);
    gel(ya,i) = dbltor(y[i-1]);
  }
  if (flag & PLOT_PARAMETRIC) {
    tas = new_chunk(4);
    for (j = 1; j <= 4; j++) gel(tas,j-1) = stoi(j);
    quark_gen = cgetg(3, t_VEC);
  }
  else tas = NULL; /* for lint */
  for (i = 0; i <= lx - 4; i++) {
    pari_sp av = avma;

    xa++; ya++;
    if (flag & PLOT_PARAMETRIC) {
      gel(quark_gen,1) = polint_i(tas, xa, pol_x[0], 4, NULL);
      gel(quark_gen,2) = polint_i(tas, ya, pol_x[0], 4, NULL);
    } else {
      quark_gen = polint_i(xa, ya, pol_x[0], 4, NULL);
      tas = xa;
    }
    rectploth(ne, var0,
               i==0 ? gel(tas,0) : gel(tas,1),
               i==lx-4 ? gel(tas,3) : gel(tas,2),
               QUARK,
               DEFAULTPREC,		/* XXXX precision */
               PLOT_RECURSIVE
                 | PLOT_NO_RESCALE
                 | PLOT_NO_FRAME
                 | PLOT_NO_AXE_Y
                 | PLOT_NO_AXE_X
                 | (flag & PLOT_PARAMETRIC),
               2);			/* Start with 3 points */
    avma = av;
  }
  avma = oldavma;
}

/* Plot a dblPointList. Complete with axes, bounding box, etc.
 * We use two drawing rectangles: one for strings, another for graphs.
 *
 * data is an array of structs. Its meaning depends on flags :
 *
 * + data[0] contains global extremas, the number of curves to plot
 *   (data[0].nb) and a list of doubles (first set of x-coordinates).
 *
 * + data[i].nb (i>0) contains the number of points in the list
 *   data[i].d (hopefully, data[2i].nb=data[2i+1].nb when i>0...)
 *
 * + If flags contain PLOT_PARAMETRIC, the array length should be
 *   even, and successive pairs (data[2i].d, data[2i+1].d) represent
 *   curves to plot.
 *
 * + If there is no such flag, the first element is an array with
 *   x-coordinates and the following ones contain y-coordinates.
 *
 * Additional flags: PLOT_NO_AXE_X, PLOT_NO_AXE_Y, PLOT_NO_FRAME. */
static GEN
rectplothrawin(long stringrect, long drawrect, dblPointList *data,
               long flags, PARI_plot *WW)
{
  GEN res;
  dblPointList y,x;
  double xsml,xbig,ysml,ybig,tmp;
  long ltype;
  pari_sp ltop=avma;
  long i,nc,nbpoints, w[2], wx[2], wy[2];

  w[0]=stringrect; w[1]=drawrect;
  if (!data) return cgetg(1,t_VEC);
  x = data[0]; nc = x.nb;
  xsml = x.xsml; xbig = x.xbig;
  ysml = x.ysml; ybig = x.ybig;
  if (xbig-xsml < 1.e-9)
  {
    tmp=fabs(xsml)/10; if (!tmp) tmp=0.1;
    xbig+=tmp; xsml-=tmp;
  }
  if (ybig-ysml < 1.e-9)
  {
    tmp=fabs(ysml)/10; if (!tmp) tmp=0.1;
    ybig+=tmp; ysml-=tmp;
  }

  if (WW)
  { /* no rectwindow supplied ==> ps or screen output */
    char c1[16],c2[16],c3[16],c4[16];
    PARI_plot W = *WW;
    long lm = W.fwidth*10; /* left margin   */
    long rm = W.hunit-1; /* right margin  */
    long tm = W.vunit-1; /* top margin    */
    long bm = W.vunit+W.fheight-1; /* bottom margin */
    long is = W.width - (lm+rm);
    long js = W.height - (tm+bm);

    wx[0]=wy[0]=0; wx[1]=lm; wy[1]=tm;
   /* Window size (W.width x W.height) is given in pixels, and
    * correct pixels are 0..W.width-1.
    * On the other hand, rect functions work with windows whose pixel
    * range is [0,width]. */
    initrect(stringrect, W.width-1, W.height-1);
    if (drawrect != stringrect) initrect(drawrect, is-1, js-1);

    /* draw labels on stringrect */
    sprintf(c1,"%.5g",ybig); sprintf(c2,"%.5g",ysml);
    sprintf(c3,"%.5g",xsml); sprintf(c4,"%.5g",xbig);

    rectlinetype(stringrect,-2); /* Frame */
    current_color[stringrect]=BLACK;
    put_string( stringrect, lm, 0, c1,
		RoSTdirRIGHT | RoSTdirHGAP | RoSTdirTOP);
    put_string(stringrect, lm, W.height - bm, c2,
		RoSTdirRIGHT | RoSTdirHGAP | RoSTdirVGAP);
    put_string(stringrect, lm, W.height - bm, c3,
		RoSTdirLEFT | RoSTdirTOP);
    put_string(stringrect, W.width - rm - 1, W.height - bm, c4,
		RoSTdirRIGHT | RoSTdirTOP);
  }
  RHasGraph(check_rect(drawrect)) = 1;

  if (!(flags & PLOT_NO_RESCALE))
    rectscale0(drawrect, xsml, xbig, ysml, ybig);

  if (!(flags & PLOT_NO_FRAME))
  {
    int do_double = (flags & PLOT_NODOUBLETICK) ? TICKS_NODOUBLE : 0;
    PARI_plot *pl = WW;
    if (!pl) { PARI_get_plot(0); pl = &pari_plot; }

    rectlinetype(drawrect, -2); 		/* Frame. */
    current_color[drawrect]=BLACK;
    rectmove0(drawrect,xsml,ysml,0);
    rectbox0(drawrect,xbig,ybig,0);
    if (!(flags & PLOT_NO_TICK_X)) {
      rectticks(pl, drawrect, xsml, ysml, xbig, ysml, xsml, xbig,
	TICKS_CLOCKW | do_double);
      rectticks(pl, drawrect, xbig, ybig, xsml, ybig, xbig, xsml,
	TICKS_CLOCKW | do_double);
    }
    if (!(flags & PLOT_NO_TICK_Y)) {
      rectticks(pl, drawrect, xbig, ysml, xbig, ybig, ysml, ybig,
	TICKS_CLOCKW | do_double);
      rectticks(pl, drawrect, xsml, ybig, xsml, ysml, ybig, ysml,
	TICKS_CLOCKW | do_double);
    }
  }

  if (!(flags & PLOT_NO_AXE_Y) && (xsml<=0 && xbig >=0))
  {
    rectlinetype(drawrect, -1); 		/* Axes. */
    current_color[drawrect]=BLUE;
    rectmove0(drawrect,0.0,ysml,0);
    rectline0(drawrect,0.0,ybig,0);
  }

  if (!(flags & PLOT_NO_AXE_X) && (ysml<=0 && ybig >=0))
  {
    rectlinetype(drawrect, -1); 		/* Axes. */
    current_color[drawrect]=BLUE;
    rectmove0(drawrect,xsml,0.0,0);
    rectline0(drawrect,xbig,0.0,0);
  }

  i = (flags & PLOT_PARAMETRIC)? 0: 1;
  for (ltype = 0; ltype < nc; ltype++)
  {
    current_color[drawrect] = ltype&1 ? GREEN: RED;
    if (flags & PLOT_PARAMETRIC) x = data[i++];

    y=data[i++]; nbpoints=y.nb;
    if ((flags & PLOT_POINTS_LINES) || (flags & PLOT_POINTS)) {
	rectlinetype(drawrect, rectpoint_itype + ltype); 	/* Graphs. */
	rectpointtype(drawrect, rectpoint_itype + ltype); 	/* Graphs. */
	rectpoints0(drawrect,x.d,y.d,nbpoints);
    }
    if ((flags & PLOT_POINTS_LINES) || !(flags & PLOT_POINTS)) {
	if (flags & PLOT_SPLINES) {
	    /* rectsplines will call us back with ltype == 0 */
	    int old = rectline_itype;

	    rectline_itype = rectline_itype + ltype;
	    rectsplines(drawrect,x.d,y.d,nbpoints,flags);	
	    rectline_itype = old;
	} else {
	    rectlinetype(drawrect, rectline_itype + ltype); 	/* Graphs. */
	    rectlines0(drawrect,x.d,y.d,nbpoints,0);	
	}
    }
  }
  for (i--; i>=0; i--) free(data[i].d);
  free(data);

  if (WW)
  {
    if (flags & PLOT_POSTSCRIPT)
      postdraw0(w,wx,wy,2);
    else
      rectdraw0(w,wx,wy,2);

    killrect(drawrect); if (stringrect != drawrect) killrect(stringrect);
  }

  avma=ltop;
  res = cgetg(5,t_VEC);
  gel(res,1) = dbltor(xsml); gel(res,2) = dbltor(xbig);
  gel(res,3) = dbltor(ysml); gel(res,4) = dbltor(ybig);
  return res;
}

/*************************************************************************/
/*                                                                       */
/*                          HI-RES FUNCTIONS                             */
/*                                                                       */
/*************************************************************************/

GEN
rectploth(long drawrect,entree *ep,GEN a,GEN b,char *ch,
          long prec,ulong flags,long testpoints)
{
  dblPointList *pl=rectplothin(ep, a,b, ch, prec, flags,testpoints);
  return rectplothrawin(0,drawrect, pl, flags,NULL);
}

GEN
rectplothraw(long drawrect, GEN data, long flags)
{
  dblPointList *pl=gtodblList(data,flags);
  return rectplothrawin(0,drawrect,pl,flags,NULL);
}

static PARI_plot*
init_output(long flags)
{
  if (flags & PLOT_POSTSCRIPT)
    { PARI_get_psplot(); return &pari_psplot; }
  else
    { PARI_get_plot(0); return &pari_plot; }
}

static GEN
ploth0(long stringrect,long drawrect,entree *ep,GEN a,GEN b,char *ch,
             long prec,ulong flags,long testpoints)
{
  PARI_plot *output = init_output(flags);
  dblPointList *pl=rectplothin(ep, a,b, ch, prec, flags,testpoints);
  return rectplothrawin(stringrect,drawrect, pl, flags, output);
}

static GEN
plothraw0(long stringrect, long drawrect, GEN listx, GEN listy, long flags)
{
  PARI_plot *output = init_output(flags);
  long data[] = {evaltyp(t_VEC) | _evallg(3), 0, 0};
  dblPointList *pl;

  gel(data,1) = listx;
  gel(data,2) = listy;
  pl=gtodblList(data,PLOT_PARAMETRIC);
  if (!pl) return cgetg(1,t_VEC);
  return rectplothrawin(stringrect,drawrect,pl,flags | PLOT_PARAMETRIC,output);
}

GEN
plothraw(GEN listx, GEN listy, long flags)
{
  flags = (flags > 1 ? flags : (flags ? 0: PLOT_POINTS));
  return plothraw0(STRINGRECT, DRAWRECT, listx, listy, flags);
}

GEN
ploth(entree *ep, GEN a, GEN b, char *ch, long prec,long flags,long numpoints)
{
  return ploth0(STRINGRECT, DRAWRECT, ep,a,b,ch,prec,flags,numpoints);
}

GEN
ploth2(entree *ep, GEN a, GEN b, char *ch, long prec)
{
  return ploth0(STRINGRECT, DRAWRECT, ep,a,b,ch,prec,PLOT_PARAMETRIC,0);
}

GEN
plothmult(entree *ep, GEN a, GEN b, char *ch, long prec)
{
  return ploth0(STRINGRECT, DRAWRECT, ep,a,b,ch,prec,0,0);
}

GEN
postplothraw(GEN listx, GEN listy, long flags)
{
  flags=(flags)? 0: PLOT_POINTS;
  return plothraw0(STRINGRECT, DRAWRECT, listx, listy, flags|PLOT_POSTSCRIPT);
}

GEN
postploth(entree *ep, GEN a, GEN b, char *ch, long prec,long flags,
           long numpoints)
{
  return ploth0(STRINGRECT,DRAWRECT,ep,a,b,ch,prec,flags|PLOT_POSTSCRIPT,
                numpoints);
}

GEN
postploth2(entree *ep, GEN a, GEN b, char *ch, long prec,
           long numpoints)
{
  return ploth0(STRINGRECT,DRAWRECT,ep,a,b,ch,prec,
		PLOT_PARAMETRIC|PLOT_POSTSCRIPT,numpoints);
}

GEN
plothsizes(void)
{
  return plothsizes_flag(0);
}

GEN
plothsizes_flag(long flag)
{
  GEN vect = cgetg(1+6,t_VEC);

  PARI_get_plot(0);
  gel(vect,1) = stoi(pari_plot.width);
  gel(vect,2) = stoi(pari_plot.height);
  if (flag) {
    gel(vect,3) = dbltor(pari_plot.hunit*1.0/pari_plot.width);
    gel(vect,4) = dbltor(pari_plot.vunit*1.0/pari_plot.height);
    gel(vect,5) = dbltor(pari_plot.fwidth*1.0/pari_plot.width);
    gel(vect,6) = dbltor(pari_plot.fheight*1.0/pari_plot.height);
  } else {
    gel(vect,3) = stoi(pari_plot.hunit);
    gel(vect,4) = stoi(pari_plot.vunit);
    gel(vect,5) = stoi(pari_plot.fwidth);
    gel(vect,6) = stoi(pari_plot.fheight);
  }
  return vect;
}	

void
plot_count(long *w, long lw, col_counter rcolcnt)
{
  RectObj *O;
  long col, i;

  for (col = 1; col < MAX_COLORS; col++)
    for (i = 0; i < ROt_MAX; i++) rcolcnt[col][i] = 0;
  for (i = 0; i < lw; i++)
  {
    PariRect *e = rectgraph[w[i]];
    for (O = RHead(e); O; O=RoNext(O))
      switch(RoType(O))
      {
	case ROt_MP : rcolcnt[RoCol(O)][ROt_PT] += RoMPcnt(O);
	              break;                 /* Multiple Point */
	case ROt_PT :                        /* Point */
	case ROt_LN :                        /* Line */
	case ROt_BX :                        /* Box */
	case ROt_ML :                        /* Multiple lines */
	case ROt_ST : rcolcnt[RoCol(O)][RoType(O)]++;
	              break;                 /* String */
      }
  }
}
/*************************************************************************/
/*                                                                       */
/*                         POSTSCRIPT OUTPUT                             */
/*                                                                       */
/*************************************************************************/

static void
PARI_get_psplot(void)
{
  if (pari_psplot.init) return;
  pari_psplot.init = 1;

  pari_psplot.width = 1120 - 60; /* 1400 - 60 for hi-res */
  pari_psplot.height=  800 - 40; /* 1120 - 60 for hi-res */
  pari_psplot.fheight= 15;
  pari_psplot.fwidth = 6;
  pari_psplot.hunit = 5;
  pari_psplot.vunit = 5;
}

static void
gendraw(GEN list, long ps, long flag)
{
  long i,n,ne,*w,*x,*y;

  if (typ(list) != t_VEC) pari_err(talker,"not a vector in rectdraw");
  n = lg(list)-1; if (!n) return;
  if (n%3) pari_err(talker,"incorrect number of components in rectdraw");
  n = n/3;
  w = (long*)gpmalloc(n*sizeof(long));
  x = (long*)gpmalloc(n*sizeof(long));
  y = (long*)gpmalloc(n*sizeof(long));
  if (flag)
    PARI_get_plot(0);
  for (i=0; i<n; i++)
  {
    GEN win = gel(list,3*i+1), x0 = gel(list,3*i+2), y0 = gel(list,3*i+3);
    long xi, yi;
    if (typ(win)!=t_INT) pari_err(typeer,"rectdraw");
    if (flag) {
      xi = DTOL(gtodouble(x0)*(pari_plot.width - 1));
      yi = DTOL(gtodouble(y0)*(pari_plot.height - 1));
    } else {
      if (typ(x0)!=t_INT || typ(y0)!= t_INT) pari_err(typeer,"rectdraw");
      xi = itos(x0);
      yi = itos(y0);
    }
    x[i] = xi;
    y[i] = yi;
    ne = itos(win); check_rect(ne);
    w[i] = ne;
  }
  if (ps) postdraw00(w,x,y,n,flag); else rectdraw0(w,x,y,n);
  free(x); free(y); free(w);
}

void
postdraw(GEN list) { gendraw(list, 1, 0); }

void
rectdraw(GEN list) { gendraw(list, 0, 0); }

void
postdraw_flag(GEN list, long flag) { gendraw(list, 1, flag); }

void
rectdraw_flag(GEN list, long flag) { gendraw(list, 0, flag); }

void
postdraw0(long *w, long *x, long *y, long lw)
{
  postdraw00(w, x, y, lw, 0);
}

static void
ps_sc(void *data, long col) { (void)data; (void)col; }

static void
ps_point(void *data, long x, long y)
{
  fprintf((FILE*)data,"%ld %ld p\n",y,x);
}

static void
ps_line(void *data, long x1, long y1, long x2, long y2)
{
  fprintf((FILE*)data,"%ld %ld m %ld %ld l\n",y1,x1,y2,x2);
  fprintf((FILE*)data,"stroke\n");
}

static void
ps_rect(void *data, long x, long y, long w, long h)
{
  fprintf((FILE*)data,"%ld %ld m %ld %ld l %ld %ld l %ld %ld l closepath\n",y,x, y,x+w, y+h,x+w, y+h,x);
}

static void
ps_points(void *data, long nb, struct plot_points *p)
{
  long i;
  for (i=0; i<nb; i++) ps_point(data, p[i].x, p[i].y);
}

static void
ps_lines(void *data, long nb, struct plot_points *p)
{
  FILE *psfile = (FILE*)data;
  long i;
  fprintf(psfile,"%ld %ld m\n",p[0].y,p[0].x);
  for (i=1; i<nb; i++) fprintf(psfile, "%ld %ld l\n", p[i].y, p[i].x);
  fprintf(psfile,"stroke\n");
}

static void
ps_string(void *data, long x, long y, char *s, long length)
{
  FILE *psfile = (FILE*)data;
  (void)length;
  if (strpbrk(s, "(\\)")) {
    fprintf(psfile,"(");
    while (*s) {
      if ( *s=='(' || *s==')' || *s=='\\' ) fputc('\\', psfile);
      fputc(*s, psfile);
      s++;
    }
  } else
    fprintf(psfile,"(%s", s);
  fprintf(psfile,") %ld %ld m 90 rotate show -90 rotate\n", y, x);
}

void
postdraw00(long *w, long *x, long *y, long lw, long scale)
{
  struct plot_eng plot;
  FILE *psfile;
  double xscale = 0.65, yscale = 0.65;
  long fontsize = 16;

  PARI_get_psplot();
  if (scale) {
    double psxscale, psyscale;

    PARI_get_plot(0);
    psxscale = pari_psplot.width * 1.0/pari_plot.width ;
    psyscale = pari_psplot.height* 1.0/pari_plot.height;
    fontsize = (long) (fontsize/psxscale);
    xscale *= psxscale;
    yscale *= psyscale;
  }
  psfile = fopen(current_psfile, "a");
  if (!psfile)
    pari_err(openfiler,"postscript",current_psfile);

  /* Definitions taken from post terminal of Gnuplot. */
  fprintf(psfile,"%%!\n\
50 50 translate\n\
/p {moveto 0 2 rlineto 2 0 rlineto 0 -2 rlineto closepath fill} def\n\
/l {lineto} def\n\
/m {moveto} def\n\
/Times-Roman findfont %ld scalefont setfont\n\
%g %g scale\n", fontsize, yscale, xscale);

  plot.sc = &ps_sc;
  plot.pt = &ps_point;
  plot.ln = &ps_line;
  plot.bx = &ps_rect;
  plot.mp = &ps_points;
  plot.ml = &ps_lines;
  plot.st = &ps_string;
  plot.pl = &pari_psplot;

  gen_rectdraw0(&plot, (void*)psfile, w, x, y, lw, 1, 1);
  fprintf(psfile,"stroke showpage\n"); fclose(psfile);
}

void
gen_rectdraw0(struct plot_eng *eng, void *data, long *w, long *x, long *y, long lw, double xs, double ys)
{
  long i, j;
  long hgapsize = eng->pl->hunit, fheight = eng->pl->fheight;
  long vgapsize = eng->pl->vunit,  fwidth = eng->pl->fwidth;
  for(i=0; i<lw; i++)
  {
    PariRect *e = rectgraph[w[i]];
    RectObj *R;
    long x0 = x[i], y0 = y[i];
    for (R = RHead(e); R; R = RoNext(R))
    {
      switch(RoType(R))
      {
      case ROt_PT:
        eng->sc(data,RoCol(R));
        eng->pt(data, DTOL((RoPTx(R)+x0)*xs), DTOL((RoPTy(R)+y0)*ys));
        break;
      case ROt_LN:
        eng->sc(data,RoCol(R));
        eng->ln(data, DTOL((RoLNx1(R)+x0)*xs),
                      DTOL((RoLNy1(R)+y0)*ys),
                      DTOL((RoLNx2(R)+x0)*xs),
                      DTOL((RoLNy2(R)+y0)*ys));
        break;
      case ROt_BX:
        eng->sc(data,RoCol(R));
        eng->bx(data,
                DTOL((RoBXx1(R)+x0)*xs),
                DTOL((RoBXy1(R)+y0)*ys),
                DTOL((RoBXx2(R)-RoBXx1(R))*xs), 
                DTOL((RoBXy2(R)-RoBXy1(R))*ys));
        break;
      case ROt_MP:
        {
          double *ptx = RoMPxs(R); 
          double *pty = RoMPys(R);
          long     nb = RoMPcnt(R);
          struct plot_points *points = 
            (struct plot_points *) gpmalloc(sizeof(*points)*nb);
          for(j=0;j<nb;j++)
          {
            points[j].x = DTOL((ptx[j]+x0)*xs);
            points[j].y = DTOL((pty[j]+y0)*ys);
          }
          eng->sc(data,RoCol(R));
          eng->mp(data, nb, points);
          free(points);
          break;
        }
      case ROt_ML:
        {
          double *ptx = RoMLxs(R); 
          double *pty = RoMLys(R);
          long     nb = RoMLcnt(R);
          struct plot_points *points = 
            (struct plot_points *) gpmalloc(sizeof(*points)*nb);
          for(j=0;j<nb;j++)
          {
            points[j].x = DTOL((ptx[j]+x0)*xs);
            points[j].y = DTOL((pty[j]+y0)*ys);
          }
          eng->sc(data,RoCol(R));
          eng->ml(data, nb, points);
          free(points);
          break;
        }
      case ROt_ST:
        {
          long dir = RoSTdir(R);
          long hjust = dir & RoSTdirHPOS_mask, hgap  = dir & RoSTdirHGAP;
          long vjust = dir & RoSTdirVPOS_mask, vgap  = dir & RoSTdirVGAP;
          char *text = RoSTs(R);
          long l     = RoSTl(R);
          long x, y;
          long shift = (hjust == RoSTdirLEFT ? 0 :
              (hjust == RoSTdirRIGHT ? 2 : 1));
          if (hgap)
            hgap = (hjust == RoSTdirLEFT) ? hgapsize : -hgapsize;
          if (vgap)
            vgap = (vjust == RoSTdirBOTTOM) ? 2*vgapsize : -2*vgapsize;
          if (vjust != RoSTdirBOTTOM)
            vgap -= ((vjust == RoSTdirTOP) ? 2 : 1)*(fheight - 1);
          x = DTOL((RoSTx(R) + x0 + hgap - (l * fwidth * shift)/2)*xs);
          y = DTOL((RoSTy(R) + y0 - vgap/2)*ys);
          eng->sc(data,RoCol(R));
          eng->st(data, x, y, text, l); 
          break;
        }
      default: 
        break;
      }
    }
  }
}
