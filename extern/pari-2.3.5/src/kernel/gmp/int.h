#line 2 "../src/kernel/gmp/int.h"
/* $Id: int.h 5645 2004-06-21 13:40:10Z bill $

Copyright (C) 2000  The PARI group.

This file is part of the PARI/GP package.

PARI/GP is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation. It is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY WHATSOEVER.

Check the License for details. You should have received a copy of it, along
with the package; see the file 'COPYING'. If not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA. */

#define int_MSW(x) ((x)+lgefint((x))-1)
/*x being a t_INT, return a pointer to the most significant word of x.*/

#define int_LSW(x) (x+2)
/*x being a t_INT, return a pointer to the least significant word of x.*/

#define int_precW(x) ((x)-1)
/*x pointing to a mantissa word, return the previous (less significant)
 * mantissa word.*/

#define int_nextW(x) ((x)+1)
/*x pointing to a mantissa word, return the next (more significant) mantissa
 * word.*/

#define int_W(x,l) ((x)+2+(l))
/*x being a t_INT, return a pointer to the l-th least significant word of x.*/

#define PARI_KERNEL_GMP
/*This macro should not be used in libpari itself.*/
