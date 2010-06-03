/* $Id: members.c 9120 2007-10-29 10:33:14Z kb $

Copyright (C) 2000-2003  The PARI group.

This file is part of the PARI/GP package.

PARI/GP is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation. It is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY WHATSOEVER.

Check the License for details. You should have received a copy of it, along
with the package; see the file 'COPYING'. If not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA. */

#include "pari.h"
#include "paripriv.h"

/********************************************************************/
/**                                                                **/
/**                          MEMBER FUNCTIONS                      **/
/**                                                                **/
/********************************************************************/
#define is_ell(x) (typ(x) == t_VEC && lg(x)>=14)
#define is_bigell(x) (typ(x) == t_VEC && lg(x)>=20)
GEN
member_e(GEN x)
{
  x = get_primeid(x);
  if (!x) member_err("e");
  return gel(x,3);
}

GEN
member_f(GEN x)
{
  x = get_primeid(x);
  if (!x) member_err("f");
  return gel(x,4);
}

GEN
member_p(GEN x)
{
  long t; (void)get_nf(x,&t);
  if (t == typ_GAL)
    return gmael(x,2,1);
  x = get_primeid(x);
  if (!x) member_err("p");
  return gel(x,1);
}

GEN
member_bid(GEN x)
{
  long t; (void)get_nf(x,&t);
  switch(t) {
    case typ_BNR: return gel(x,2);
    case typ_BID: return x;
  }
  member_err("bid");
  return NULL;
}

GEN
member_bnf(GEN x)
{
  long t; x = get_bnf(x,&t);
  if (!x) member_err("bnf");
  return x;
}

GEN
member_nf(GEN x)
{
  long t; x = get_nf(x,&t);
  if (!x) member_err("nf");
  return x;
}

/* integral basis */
GEN
member_zk(GEN x)
{
  long t; GEN y = get_nf(x,&t);
  if (!y)
  {
    switch(t)
    {
      case typ_CLA: return gmael(x,1,4);
      case typ_Q: return mkvec2(gen_1, pol_x[varn(x[1])]);
    }
    member_err("zk");
  }
  return gel(y,7);
}

GEN
member_disc(GEN x) /* discriminant */
{
  long t; GEN y = get_nf(x,&t);
  if (!y)
  {
    switch(t)
    {
      case typ_Q  : return discsr(gel(x,1));
      case typ_CLA:
        x = gmael(x,1,3);
        if (typ(x) != t_VEC || lg(x) != 3) break;
        return gel(x,1);
      case typ_ELL: return gel(x,12);
    }
    member_err("disc");
  }
  return gel(y,3);
}

GEN
member_pol(GEN x) /* polynomial */
{
  long t; GEN y = get_nf(x,&t);
  if (!y)
  {
    switch(t)
    {
      case typ_CLA: return gmael(x,1,1);
      case typ_POL: return x;
      case typ_Q  : return gel(x,1);
      case typ_GAL: return gel(x,1);
    }
    if (typ(x)==t_POLMOD) return gel(x,2);
    if (typ(x)==t_VEC && lg(x) == 13) return gmael(x,11,1);
    member_err("pol");
  }
  return gel(y,1);
}

GEN
member_mod(GEN x) /* modulus */
{
  long t; (void)get_nf(x,&t);
  switch(t) {
    case typ_GAL: return gmael(x,2,3);
    case typ_BNR: x = gel(x,2); /* fall through */
    case typ_BID: return gel(x,1);
  }
  switch(typ(x))
  {
    case t_INTMOD: case t_POLMOD: case t_QUAD: break;
    default: member_err("mod");
  }
  return gel(x,1);
}

GEN
member_sign(GEN x) /* signature */
{
  long t; GEN y = get_nf(x,&t);
  if (!y)
  {
    if (t == typ_CLA) return gmael(x,1,2);
    member_err("sign");
  }
  return gel(y,2);
}
GEN
member_r1(GEN x) { return gel(member_sign(x), 1); }
GEN
member_r2(GEN x) { return gel(member_sign(x), 2); }

GEN
member_index(GEN x)
{
  long t; GEN y = get_nf(x,&t);
  if (!y) member_err("index");
  return gel(y,4);
}

/* x assumed to be output by get_nf: ie a t_VEC with length 11 */
static GEN
nfmats(GEN x)
{
  GEN y;
  if (!x) return NULL;
  y = gel(x,5);
  if (typ(y) == t_VEC && lg(y) != 8) return NULL;
  return y;
}

GEN
member_t2(GEN x) /* T2 matrix */
{
  long t; x = nfmats(get_nf(x,&t));
  if (!x) member_err("t2");
  return gram_matrix(gel(x,2));
}

GEN
member_diff(GEN x) /* different */
{
  long t; x = nfmats(get_nf(x,&t));
  if (!x) member_err("diff");
  return gel(x,5);
}

GEN
member_codiff(GEN x) /* codifferent */
{
  long t; GEN T, D, DinvT, nf = get_nf(x,&t), y = nfmats(nf);
  if (!y) member_err("codiff");
  T = gel(y,4);
  D = absi(gel(nf,3));
  DinvT = ZM_inv(T,D);
  return gdiv(hnfmod(DinvT, D), D);
}

GEN
member_roots(GEN x) /* roots */
{
  long t; GEN y = get_nf(x,&t);
  if (!y)
  {
    if (t == typ_ELL && is_bigell(x)) return gel(x,14);
    if (t == typ_GAL) return gel(x,3);
    member_err("roots");
  }
  return gel(y,6);
}

/* assume x output by get_bnf: ie a t_VEC with length 10 */
static GEN
check_RES(GEN x, char *s)
{
  GEN y = gel(x,8);
  if (typ(y) != t_VEC || lg(y) < 4)
    member_err(s);
  return y;
}

GEN
member_clgp(GEN x) /* class group (3-component row vector) */
{
  long t; GEN y = get_bnf(x,&t);
  if (!y)
  {
    switch(t)
    {
      case typ_QUA: return mkvec3(gel(x,1), gel(x,2), gel(x,3));
      case typ_CLA: return gmael(x,1,5);
      case typ_BID: return x = gel(x,2);
    }
    if (typ(x)==t_VEC)
      switch(lg(x))
      {
        case 3: /* no gen */
        case 4: return x;
      }
    member_err("clgp");
  }
  if (t==typ_BNR) return gel(x,5);
  y = check_RES(y, "clgp");
  return gel(y,1);
}

GEN
member_reg(GEN x) /* regulator */
{
  long t; GEN y = get_bnf(x,&t);
  if (!y)
  {
    switch(t)
    {
      case typ_CLA: return gmael(x,1,6);
      case typ_QUA: return gel(x,4);
    }
    member_err("reg");
  }
  if (t == typ_BNR) pari_err(impl,"ray regulator");
  y = check_RES(y, "reg");
  return gel(y,2);
}

GEN
member_fu(GEN x) /* fundamental units */
{
  long t; GEN y = get_bnf(x,&t);
  if (!y)
  {
    switch(t)
    {
      case typ_CLA: x = gel(x,1); if (lg(x) < 10) break;
        return gel(x,9);
      case typ_Q:
        x = discsr(gel(x,1));
        return (signe(x)<0)? cgetg(1,t_VEC): fundunit(x);
    }
    member_err("fu");
  }
  if (t == typ_BNR) pari_err(impl,"ray units");
  return basistoalg(y, check_units(y,".fu"));
}

/* torsion units. return [w,e] where w is the number of roots of 1, and e a
 * polymod generator */
GEN
member_tu(GEN x)
{
  long t; GEN y, bnf = get_bnf(x,&t), res = cgetg(3,t_VEC);
  if (!bnf)
  {
    switch(t)
    {
      case typ_Q:
        y = discsr(gel(x,1));
        if (signe(y)<0 && cmpiu(y,4)<=0) /* |y| <= 4 */
          y = stoi((itos(y) == -4)? 4: 6);
        else
        { y = gen_2; x = gen_m1; }
        gel(res,1) = y;
        gel(res,2) = x; return res;
      case typ_CLA:
        x = gel(x,1);
        if (lg(x) > 8)
        {
          y = gel(x,8);
          if (typ(y) == t_VEC || lg(y) == 3) { res[2] = y[2]; break; }
        }
      default: member_err("tu");
        return NULL; /* not reached */
    }
  }
  else
  {
    if (t == typ_BNR) pari_err(impl,"ray torsion units");
    x = gel(bnf,7);
    y = gel(bnf,8);
    if (typ(y) == t_VEC && lg(y) > 5) y = gel(y,4);
    else
    {
      y = rootsof1(x);
      gel(y,2) = gmul(gel(x,7), gel(y,2));
    }
    gel(res,2) = basistoalg(bnf, gel(y,2));
  }
  res[1] = y[1]; return res;
}

GEN
member_futu(GEN x) /*  concatenation of fu and tu, w is lost */
{
  GEN fuc = member_fu(x);
  return shallowconcat(fuc, (GEN)member_tu(x)[2]);
}

GEN
member_tufu(GEN x) /*  concatenation of tu and fu, w is lost */
{
  GEN fuc = member_fu(x);
  return shallowconcat((GEN)member_tu(x)[2], fuc);
}

GEN
member_zkst(GEN bid)
/* structure of (Z_K/m)^*, where bid is an idealstarinit (with or without gen)
   or a bnrinit (with or without gen) */
{
  if (typ(bid)==t_VEC)
    switch(lg(bid))
    {
      case 6: return gel(bid,2);   /* idealstarinit */
      case 7: bid = gel(bid,2); /* bnrinit */
        if (typ(bid) == t_VEC && lg(bid) > 2)
          return gel(bid,2);
    }
  member_err("zkst");
  return NULL; /* not reached */
}

GEN
member_no(GEN clg) /* number of elements of a group (of type clgp) */
{
  clg = member_clgp(clg);
  if (typ(clg)!=t_VEC  || (lg(clg)!=3 && lg(clg)!=4))
    member_err("no");
  return gel(clg,1);
}

GEN
member_cyc(GEN clg) /* cyclic decomposition (SNF) of a group (of type clgp) */
{
  clg = member_clgp(clg);
  if (typ(clg)!=t_VEC  || (lg(clg)!=3 && lg(clg)!=4))
    member_err("cyc");
  return gel(clg,2);
}

/* SNF generators of a group (of type clgp), or generators of a prime
 * ideal
 */
GEN
member_gen(GEN x)
{
  long t;
  GEN y = get_primeid(x);
  if (y) return mkvec2(gel(y,1), gel(y,2));
  (void)get_nf(x,&t);
  if (t == typ_GAL)
    return gel(x,7);
  x = member_clgp(x);
  if (typ(x)!=t_VEC || lg(x)!=4)
    member_err("gen");
  if (typ(x[1]) == t_COL) return gel(x,2); /* from bnfisprincipal */
  return gel(x,3);
}
GEN
member_group(GEN x)
{
  long t; (void)get_nf(x,&t);
  if (t == typ_GAL)
    return gel(x,6);
  member_err("group");
  return NULL; /* not reached */
}
GEN
member_orders(GEN x)
{
  long t; (void)get_nf(x,&t);
  if (t == typ_GAL)
    return gel(x,8);
  member_err("orders");
  return NULL; /* not reached */
}

GEN
member_a1(GEN x)
{
  if (!is_ell(x)) member_err("a1");
  return gel(x,1);
}

GEN
member_a2(GEN x)
{
  if (!is_ell(x)) member_err("a2");
  return gel(x,2);
}

GEN
member_a3(GEN x)
{
  if (!is_ell(x)) member_err("a3");
  return gel(x,3);
}

GEN
member_a4(GEN x)
{
  if (!is_ell(x)) member_err("a4");
  return gel(x,4);
}

GEN
member_a6(GEN x)
{
  if (!is_ell(x)) member_err("a6");
  return gel(x,5);
}

GEN
member_b2(GEN x)
{
  if (!is_ell(x)) member_err("b2");
  return gel(x,6);
}

GEN
member_b4(GEN x)
{
  if (!is_ell(x)) member_err("b4");
  return gel(x,7);
}

GEN
member_b6(GEN x)
{
  if (!is_ell(x)) member_err("b6");
  return gel(x,8);
}

GEN
member_b8(GEN x)
{
  if (!is_ell(x)) member_err("b8");
  return gel(x,9);
}

GEN
member_c4(GEN x)
{
  if (!is_ell(x)) member_err("c4");
  return gel(x,10);
}

GEN
member_c6(GEN x)
{
  if (!is_ell(x)) member_err("c6");
  return gel(x,11);
}

GEN
member_j(GEN x)
{
  if (!is_ell(x)) member_err("j");
  return gel(x,13);
}

GEN
member_omega(GEN x)
{
  if (!is_bigell(x)) member_err("omega");
  if (gcmp0(gel(x,19))) pari_err(talker,"curve not defined over R");
  return mkvec2(gel(x,15), gel(x,16));
}

GEN
member_eta(GEN x)
{
  if (!is_bigell(x)) member_err("eta");
  if (gcmp0(gel(x,19))) pari_err(talker,"curve not defined over R");
  return mkvec2(gel(x,17), gel(x,18));
}

GEN
member_area(GEN x)
{
  if (!is_bigell(x)) member_err("area");
  if (gcmp0(gel(x,19))) pari_err(talker,"curve not defined over R");
  return gel(x,19);
}

GEN
member_tate(GEN x)
{
  if (!is_bigell(x)) member_err("tate");
  if (!gcmp0(gel(x,19))) pari_err(talker,"curve not defined over a p-adic field");
  return mkvec3(gel(x,15), gel(x,16), gel(x,17));
}

GEN
member_w(GEN x)
{
  if (!is_bigell(x)) member_err("w");
  if (!gcmp0(gel(x,19))) pari_err(talker,"curve not defined over a p-adic field");
  return gel(x,18);
}
