/* $Id: plotnull.c 7522 2005-12-09 18:14:24Z kb $

Copyright (C) 2000  The PARI group.

This file is part of the PARI/GP package.

PARI/GP is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation. It is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY WHATSOEVER.

Check the License for details. You should have received a copy of it, along
with the package; see the file 'COPYING'. If not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA. */

#include "pari.h"
#include "rect.h"

void
rectdraw0(long *w, long *x, long *y, long lw)
{
  (void)w;
  (void)x;
  (void)y;
  (void)lw;
}

void
PARI_get_plot(long f)
{
  (void)f;
  pari_err(talker,"high resolution graphics disabled");
}
