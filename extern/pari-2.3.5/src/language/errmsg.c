/* $Id: errmsg.c 5676 2004-06-23 17:41:11Z kb $

Copyright (C) 2000  The PARI group.

This file is part of the PARI/GP package.

PARI/GP is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation. It is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY WHATSOEVER.

Check the License for details. You should have received a copy of it, along
with the package; see the file 'COPYING'. If not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA. */

char *errmessage[]=
{
/* force error into non-0 */
  "",
/* cant_deflate */
  "can't deflate",

/* Always catched up to this point */

/* caracer1 */
  "unexpected character",
/* paramer1 */
  "unknown function or error in formal parameters",
/* varer1 */
  "variable name expected",
/* obsoler */
  "obsolete function",
/* openfiler */
  "error opening ",
/* talker2 */
  "",

/* NO CONTEXT NOW */

/* talker */
  "",
/* flagerr */
  "invalid flag",
/* warner */
  "Warning:",
/* warnprec */
  "Warning: increasing prec",
/* warnfile */
  "Warning: failed to",
/* accurer */
  "accuracy problems",
/* bugparier */
  "bug in",
/* impl */
  "sorry,",
/* archer */
  "sorry, not yet available on this system",
/* warnmem */
  "collecting garbage in",
/* precer */
  "precision too low",
/* siginter */
  "",
/* typeer */
  "incorrect type",
/* consister */
  "inconsistent data",
/* user */
  "",

    /*  MP.S  */

/* affer2 */
  "impossible assignment I-->S",
/* affer3 */
  "impossible assignment I-->I",
/* errpile */
  "the PARI stack overflows !",
/* errlg */
  "length (lg) overflow",
/* errexpo */
  "exponent (expo) overflow",
/* errvalp */
  "valuation (valp) overflow",
/* rtodber */
  "overflow in R->dbl conversion",

  /*  ALGLIN.C  */

/* matinv1 */
  "non invertible matrix in gauss",
/* mattype1 */
  "not a square matrix",

  /*  ANAL.C  */

/* valencer1 */
  "unknown identifier valence, please report",

  /*  ARITH.C  */

/* arither1 */
  "not an integer argument in an arithmetic function",
/* primer1 */
  "not enough precomputed primes",
/* primer2 */
  "not enough precomputed primes, need primelimit ~ ",
/* invmoder*/
  "impossible inverse modulo: ",

  /*  BASE.C  */

/* constpoler */
  "constant polynomial",
/* notpoler */
  "not a polynomial",
/* redpoler */
  "reducible polynomial",
/* zeropoler */
  "zero polynomial",

  /*  BIBLI.C  */

/* intger2 */
  "too many iterations for desired precision in integration routine",
/* lllger3 */
  "not a definite matrix in lllgram",

  /*  ELLIPTIC.C  */

/* elliper1 */
  "bad argument for an elliptic curve related function",

  /*  GEN.C */

/* operi */
  "impossible",
/* operf */
  "forbidden",
/* gdiver */
  "division by zero",
/* overwriter */
  "trying to overwrite a universal object",

  /*  INIT.C  */

/* memer */
  "not enough memory",

  /*  TRANS.C  */

/* infprecer */
  "infinite precision",
/* negexper */
  "negative exponent",
/* sqrter5 */
  "non quadratic residue in gsqrt",

/* noer */
  "what's going on ?"
};
