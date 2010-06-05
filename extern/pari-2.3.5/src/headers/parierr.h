/* $Id: parierr.h 5674 2004-06-23 17:02:01Z kb $

Copyright (C) 2000  The PARI group.

This file is part of the PARI/GP package.

PARI/GP is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation. It is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY WHATSOEVER.

Check the License for details. You should have received a copy of it, along
with the package; see the file 'COPYING'. If not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA. */

enum {
  cant_deflate = 1,	/* Force errors into non-0 */

/* Always catched up to this point */

  caracer1,
  paramer1, varer1, obsoler, openfiler, talker2,

/* NO CONTEXT now */

  talker, flagerr, warner, warnprec, warnfile,
  accurer, bugparier, impl, archer, warnmem, precer,

  siginter, typeer, consister, user,

/* mp.s ou mp.c */

  affer2, affer3,
  
  errpile, errlg, errexpo, errvalp, rtodber,

/* alglin.c */

  matinv1, mattype1,

/* anal.c   */

  valencer1,

/* arith.c  */

  arither1, primer1, primer2, invmoder,

/* base.c   */
  constpoler, notpoler, redpoler, zeropoler,

/* bibli.c  */

  intger2, lllger3,

/* elliptic.c */

  elliper1,

/* gen.c */
  operi, operf, gdiver, overwriter,

/* init.c */

  memer,

/* trans.c */

  infprecer, negexper, sqrter5,

/* PAS D'ERREUR */

  noer
};
