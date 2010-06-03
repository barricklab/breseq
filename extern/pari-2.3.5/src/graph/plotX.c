/* $Id: plotX.c 7523 2005-12-09 19:34:24Z kb $

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
/*                       HIGH RESOLUTION PLOT                      */
/*                                                                 */
/*******************************************************************/

#include "pari.h"
#include "rect.h"
#include "../language/anal.h"

#ifdef HPPA
#  ifndef __GNUC__
     typedef char *caddr_t;
#  endif
#endif

BEGINEXTERN
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>
ENDEXTERN

static Colormap PARI_Colormap;
static XColor  *PARI_Colors;
static XColor  *PARI_ExactColors;

struct data_x
{
  Display *display;
  Window win;
  GC gc;
};

static void SetForeground(void *data, long col)
{
  struct data_x *dx = (struct data_x *) data;
  XSetForeground(dx->display,dx->gc, PARI_Colors[col].pixel);
}

static void DrawPoint(void *data, long x, long y)
{
  struct data_x *dx = (struct data_x *) data;
  XDrawPoint(dx->display,dx->win,dx->gc, x,y);
}

static void DrawLine(void *data, long x1, long y1, long x2, long y2)
{
  struct data_x *dx = (struct data_x *) data;
  XDrawLine(dx->display,dx->win,dx->gc, x1,y1, x2,y2);
}

static void DrawRectangle(void *data, long x, long y, long w, long h)
{
  struct data_x *dx = (struct data_x *) data;
  XDrawRectangle(dx->display,dx->win,dx->gc, x,y, w,h);
}

static void DrawPoints(void *data, long nb, struct plot_points *p)
{
  struct data_x *dx = (struct data_x *) data;
  XPoint *xp=(XPoint*)gpmalloc(sizeof(xp)*nb);
  long i;
  for (i=0;i<nb;i++)
  {
    xp[i].x=p[i].x;
    xp[i].y=p[i].y;
  }
  XDrawPoints(dx->display,dx->win,dx->gc, xp, nb, 0);
  free(xp);
}

static void DrawLines(void *data, long nb, struct plot_points *p)
{
  struct data_x *dx = (struct data_x *) data;
  XPoint *xp=(XPoint*)gpmalloc(sizeof(xp)*nb);
  long i;
  for (i=0;i<nb;i++)
  {
    xp[i].x=p[i].x;
    xp[i].y=p[i].y;
  }
  XDrawLines(dx->display,dx->win,dx->gc, xp, nb, 0);
  free(xp);
}

static void DrawString(void *data, long x, long y, char *text, long numtext)
{
  struct data_x *dx = (struct data_x *) data;
  XDrawString(dx->display,dx->win,dx->gc, x,y, text, numtext);
}

static char *PARI_DefaultColors[MAX_COLORS] =
{
  " ",
  "black",    /* Default */
  "blue",     /* Axes */
  "violetred",   /* Odd numbered curves in ploth */
  "red",      /* Curves, or Even numbered curves in ploth */
  "green",
  "grey",
  "gainsboro",
};

static void
PARI_ColorSetUp(Display *display, char **Colors, int n)
{
  static int init_done = 0;
  int i;

  if (init_done) return;
  init_done=1;

  PARI_Colormap = DefaultColormap(display, 0);
  PARI_Colors = (XColor *) gpmalloc((n+1) * sizeof(XColor));
  PARI_ExactColors = (XColor *) gpmalloc((n+1) * sizeof(XColor));
  for (i=1; i<n; i++)
    XAllocNamedColor(display, PARI_Colormap, Colors[i],
		     &PARI_ExactColors[i], &PARI_Colors[i]);
}

/* after fork(), we don't want the child to recover but to exit */
static void
exiterr(char *str)
{
  term_color(c_ERR);
  fprintferr("\n  *** X fatal error: %s\n",str);
  term_color(c_NONE); exit(1);
}

#define MAX_BUF 256

static int
Xerror(Display *d, XErrorEvent *pari_err) {
  char buf[MAX_BUF];
  XGetErrorText(d,pari_err->error_code,buf,MAX_BUF);
  exiterr(buf); return 0;
}

static int
IOerror(Display *d) {
  char buf[MAX_BUF];
  sprintf(buf, "lost display on %s", DisplayString(d));
  exiterr(buf); return 0;
}

void
rectdraw0(long *w, long *x, long *y, long lw)
{
  long oldwidth,oldheight;
  struct plot_eng plotX;
  struct data_x dx;
  double xs = 1, ys = 1;
  int screen;
  Display *display;
  GC gc;
  Window win;
  XEvent event;
  XSizeHints size_hints;
  XFontStruct *font_info;
  XSetWindowAttributes attrib;
  Atom wm_delete_window, wm_protocols;

  if (fork()) return;  /* parent process returns */

  pari_close();
  PARI_get_plot(1);

  display = XOpenDisplay(NULL);
  font_info = XLoadQueryFont(display, "9x15");
  if (!font_info) exiterr("cannot open 9x15 font");
  XSetErrorHandler(Xerror);
  XSetIOErrorHandler(IOerror);
  PARI_ColorSetUp(display,PARI_DefaultColors,MAX_COLORS);

  screen = DefaultScreen(display);
  win = XCreateSimpleWindow
    (display, RootWindow(display, screen), 0, 0,
     pari_plot.width, pari_plot.height, 4,
     BlackPixel(display, screen), WhitePixel(display, screen));

  size_hints.flags = PPosition | PSize;
  size_hints.x = 0;
  size_hints.y = 0;
  size_hints.width  = pari_plot.width;
  size_hints.height = pari_plot.height;
  XSetStandardProperties
    (display, win, "rectplot", NULL, None, NULL, 0, &size_hints);

  wm_delete_window = XInternAtom(display, "WM_DELETE_WINDOW", False);
  wm_protocols = XInternAtom(display, "WM_PROTOCOLS", False);
  XSetWMProtocols(display,win,&wm_delete_window, 1);

  XSelectInput (display, win,
    ExposureMask | ButtonPressMask | StructureNotifyMask);

  /* enable backing-store */
  attrib.backing_store = Always;
  attrib.backing_planes = AllPlanes;
  XChangeWindowAttributes(display,win,CWBackingStore|CWBackingPlanes,&attrib);

  gc = XCreateGC(display, win, 0, NULL);
  XSetFont(display, gc, font_info->fid);

  XClearWindow(display, win);
  XMapWindow(display, win);
  oldwidth  = pari_plot.width;
  oldheight = pari_plot.height;
  dx.display= display;
  dx.win = win;
  dx.gc = gc;
  plotX.sc = &SetForeground;
  plotX.pt = &DrawPoint;
  plotX.ln = &DrawLine;
  plotX.bx = &DrawRectangle;
  plotX.mp = &DrawPoints;
  plotX.ml = &DrawLines;
  plotX.st = &DrawString;
  plotX.pl = &pari_plot;

  for(;;)
  {
    XNextEvent(display, &event);
    switch(event.type)
    {
      case ClientMessage:
        if (event.xclient.message_type != wm_protocols ||
            (Atom)event.xclient.data.l[0] != wm_delete_window) break;
      case ButtonPress:
      case DestroyNotify:
	XUnloadFont(display,font_info->fid); XFreeGC(display,gc);
	free_graph();
	XCloseDisplay(display); exit(0);

      case ConfigureNotify:
      {
        int width  = event.xconfigure.width;
        int height = event.xconfigure.height;

        if (width == oldwidth && height == oldheight) break;
        oldwidth  = width;
        oldheight = height; 

        /* recompute scale */
	xs = ((double)width)/pari_plot.width;
        ys = ((double)height)/pari_plot.height;
      }
      case Expose: 
        gen_rectdraw0(&plotX, (void *)&dx, w, x, y,lw,xs,ys);
    }
  }
}

void
PARI_get_plot(long fatal)
{
  Display *display;
  int screen;

  if (pari_plot.init) return;
  if (!(display = XOpenDisplay(NULL)))
  {
    if (fatal) exiterr("no X server");
    pari_err(talker, "no X server");
  }
  screen = DefaultScreen(display);
  pari_plot.width  = DisplayWidth(display, screen) - 40;
  pari_plot.height = DisplayHeight(display, screen) - 60;
  pari_plot.fheight = 15;
  pari_plot.fwidth  = 9;
  pari_plot.hunit   = 5;
  pari_plot.vunit   = 5;
  pari_plot.init = 1;
  XCloseDisplay(display);
}

