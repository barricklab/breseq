#line 2 "../src/kernel/none/invmod.c"
/* $Id: invmod.c 7534 2005-12-12 08:58:13Z kb $

Copyright (C) 2003  The PARI group.

This file is part of the PARI/GP package.

PARI/GP is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation. It is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY WHATSOEVER.

Check the License for details. You should have received a copy of it, along
with the package; see the file 'COPYING'. If not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA. */

/*==================================
 * invmod(a,b,res)
 *==================================
 *    If a is invertible, return 1, and set res  = a^{ -1 }
 *    Otherwise, return 0, and set res = gcd(a,b)
 *
 * This is sufficiently different from bezout() to be implemented separately
 * instead of having a bunch of extra conditionals in a single function body
 * to meet both purposes.
 */

#ifdef INVMOD_PARI
INLINE int
invmod_pari(GEN a, GEN b, GEN *res)
#else
int
invmod(GEN a, GEN b, GEN *res)
#endif
{
  GEN v,v1,d,d1,q,r;
  pari_sp av, av1, lim;
  long s;
  ulong g;
  ulong xu,xu1,xv,xv1;		/* Lehmer stage recurrence matrix */
  int lhmres;			/* Lehmer stage return value */

  if (typ(a) != t_INT || typ(b) != t_INT) pari_err(arither1);
  if (!signe(b)) { *res=absi(a); return 0; }
  av = avma;
  if (lgefint(b) == 3) /* single-word affair */
  {
    ulong d1 = umodiu(a, (ulong)(b[2]));
    if (d1 == 0)
    {
      if (b[2] == 1L)
        { *res = gen_0; return 1; }
      else
        { *res = absi(b); return 0; }
    }
    g = xgcduu((ulong)(b[2]), d1, 1, &xv, &xv1, &s);
#ifdef DEBUG_LEHMER
    fprintferr(" <- %lu,%lu\n", (ulong)(b[2]), (ulong)(d1[2]));
    fprintferr(" -> %lu,%ld,%lu; %lx\n", g,s,xv1,avma);
#endif
    avma = av;
    if (g != 1UL) { *res = utoipos(g); return 0; }
    xv = xv1 % (ulong)(b[2]); if (s < 0) xv = ((ulong)(b[2])) - xv;
    *res = utoipos(xv); return 1;
  }

  (void)new_chunk(lgefint(b));
  d = absi(b); d1 = modii(a,d);

  v=gen_0; v1=gen_1;	/* general case */
#ifdef DEBUG_LEHMER
  fprintferr("INVERT: -------------------------\n");
  output(d1);
#endif
  av1 = avma; lim = stack_lim(av,1);

  while (lgefint(d) > 3 && signe(d1))
  {
#ifdef DEBUG_LEHMER
    fprintferr("Calling Lehmer:\n");
#endif
    lhmres = lgcdii((ulong*)d, (ulong*)d1, &xu, &xu1, &xv, &xv1, MAXULONG);
    if (lhmres != 0)		/* check progress */
    {				/* apply matrix */
#ifdef DEBUG_LEHMER
      fprintferr("Lehmer returned %d [%lu,%lu;%lu,%lu].\n",
	      lhmres, xu, xu1, xv, xv1);
#endif
      if ((lhmres == 1) || (lhmres == -1))
      {
	if (xv1 == 1)
	{
	  r = subii(d,d1); d=d1; d1=r;
	  a = subii(v,v1); v=v1; v1=a;
	}
	else
	{
	  r = subii(d, mului(xv1,d1)); d=d1; d1=r;
	  a = subii(v, mului(xv1,v1)); v=v1; v1=a;
	}
      }
      else
      {
	r  = subii(muliu(d,xu),  muliu(d1,xv));
	a  = subii(muliu(v,xu),  muliu(v1,xv));
	d1 = subii(muliu(d,xu1), muliu(d1,xv1)); d = r;
	v1 = subii(muliu(v,xu1), muliu(v1,xv1)); v = a;
        if (lhmres&1)
	{
          setsigne(d,-signe(d));
          setsigne(v,-signe(v));
        }
        else
	{
          if (signe(d1)) { setsigne(d1,-signe(d1)); }
          setsigne(v1,-signe(v1));
        }
      }
    }
#ifdef DEBUG_LEHMER
    else
      fprintferr("Lehmer returned 0.\n");
    output(d); output(d1); output(v); output(v1);
    sleep(1);
#endif

    if (lhmres <= 0 && signe(d1))
    {
      q = dvmdii(d,d1,&r);
#ifdef DEBUG_LEHMER
      fprintferr("Full division:\n");
      printf("  q = "); output(q); sleep (1);
#endif
      a = subii(v,mulii(q,v1));
      v=v1; v1=a;
      d=d1; d1=r;
    }
    if (low_stack(lim, stack_lim(av,1)))
    {
      GEN *gptr[4]; gptr[0]=&d; gptr[1]=&d1; gptr[2]=&v; gptr[3]=&v1;
      if(DEBUGMEM>1) pari_warn(warnmem,"invmod");
      gerepilemany(av1,gptr,4);
    }
  } /* end while */

  /* Postprocessing - final sprint */
  if (signe(d1))
  {
    /* Assertions: lgefint(d)==lgefint(d1)==3, and
     * gcd(d,d1) is nonzero and fits into one word
     */
    g = xxgcduu((ulong)d[2], (ulong)d1[2], 1, &xu, &xu1, &xv, &xv1, &s);
#ifdef DEBUG_LEHMER
    output(d);output(d1);output(v);output(v1);
    fprintferr(" <- %lu,%lu\n", (ulong)d[2], (ulong)d1[2]);
    fprintferr(" -> %lu,%ld,%lu; %lx\n", g,s,xv1,avma);
#endif
    if (g != 1UL) { avma = av; *res = utoipos(g); return 0; }
    /* (From the xgcduu() blurb:)
     * For finishing the multiword modinv, we now have to multiply the
     * returned matrix  (with properly adjusted signs)  onto the values
     * v' and v1' previously obtained from the multiword division steps.
     * Actually, it is sufficient to take the scalar product of [v',v1']
     * with [u1,-v1], and change the sign if s==1.
     */
    v = subii(muliu(v,xu1),muliu(v1,xv1));
    if (s > 0) setsigne(v,-signe(v));
    avma = av; *res = modii(v,b);
#ifdef DEBUG_LEHMER
    output(*res); fprintfderr("============================Done.\n");
    sleep(1);
#endif
    return 1;
  }
  /* get here when the final sprint was skipped (d1 was zero already) */
  avma = av;
  if (!equalii(d,gen_1)) { *res = icopy(d); return 0; }
  *res = modii(v,b);
#ifdef DEBUG_LEHMER
  output(*res); fprintferr("============================Done.\n");
  sleep(1);
#endif
  return 1;
}

