/* $Id: paritype.h 5793 2004-09-01 23:48:37Z kb $

Copyright (C) 2000  The PARI group.

This file is part of the PARI/GP package.

PARI/GP is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation. It is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY WHATSOEVER.

Check the License for details. You should have received a copy of it, along
with the package; see the file 'COPYING'. If not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA. */

/* DO NOT REORDER THESE
 * actual values can be changed. Adapt lontyp/lontyp2 in gen2.c */
enum {
  t_INT    =  1,
  t_REAL   =  2,
  t_INTMOD =  3,
  t_FRAC   =  4,
  t_COMPLEX=  6,
  t_PADIC  =  7,
  t_QUAD   =  8,
  t_POLMOD =  9,
  t_POL    =  10,
  t_SER    =  11,
  t_RFRAC  =  13,
  t_QFR    =  15,
  t_QFI    =  16,
  t_VEC    =  17,
  t_COL    =  18,
  t_MAT    =  19,
  t_LIST   =  20,
  t_STR    =  21,
  t_VECSMALL= 22
};
#define is_const_t(t) ((t) < t_POLMOD)
#define is_extscalar_t(t) ((t) <= t_POL)
#define is_graphicvec_t(t) ( (t) >= t_QFR && (t) <= t_MAT )
#define is_intreal_t(t) ( (t) <= t_REAL )
#define is_matvec_t(t) ( (t) >= t_VEC && (t) <= t_MAT )
#define is_noncalc_t(tx) ((tx) >= t_LIST)
#define is_rational_t(t) ((t) == t_INT || (t) == t_FRAC)
#define is_recursive_t(t) (lontyp[t])
#define is_scalar_t(t) ((t) < t_POL)
#define is_vec_t(t) ( (t) == t_VEC || (t) == t_COL )

/* backward compatibility */
#define t_FRACN  t_FRAC
#define t_RFRACN t_RFRAC
#define is_frac_t(t) ( (t) == t_FRAC )
#define is_rfrac_t(t) ( (t) == t_RFRAC )
