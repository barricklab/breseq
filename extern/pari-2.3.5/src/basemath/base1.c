/* $Id: base1.c 8363 2007-03-11 22:48:52Z kb $

Copyright (C) 2000  The PARI group.

This file is part of the PARI/GP package.

PARI/GP is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation. It is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY WHATSOEVER.

Check the License for details. You should have received a copy of it, along
with the package; see the file 'COPYING'. If not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA. */

/**************************************************************/
/*                                                            */
/*                        NUMBER FIELDS                       */
/*                                                            */
/**************************************************************/
#include "pari.h"
#include "paripriv.h"

int new_galois_format = 0;

void
checkrnf(GEN rnf)
{
  if (typ(rnf)!=t_VEC || lg(rnf)!=13) pari_err(typeer,"checkrnf");
}

GEN
checkbnf_i(GEN bnf)
{
  if (typ(bnf) == t_VEC)
    switch (lg(bnf))
    {
      case 11: return bnf;
      case 7:  return checkbnf_i(gel(bnf,1));
    }
  return NULL;
}

GEN
checknf_i(GEN nf)
{
  if (typ(nf)==t_VEC)
    switch(lg(nf))
    {
      case 10: return nf;
      case 11: return checknf_i(gel(nf,7));
      case 7:  return checknf_i(gel(nf,1));
      case 3: if (typ(nf[2]) == t_POLMOD) return checknf_i(gel(nf,1));
    }
  return NULL;
}

GEN
checkbnf(GEN x)
{
  GEN bnf = checkbnf_i(x);
  if (!bnf)
  {
    if (checknf_i(x)) pari_err(talker,"please apply bnfinit first");
    pari_err(typeer,"checkbnf");
  }
  return bnf;
}

GEN
checknf(GEN x)
{
  GEN nf = checknf_i(x);
  if (!nf)
  {
    if (typ(x)==t_POL) pari_err(talker,"please apply nfinit first");
    pari_err(typeer,"checknf");
  }
  return nf;
}

void
checkbnr(GEN bnr)
{
  if (typ(bnr)!=t_VEC || lg(bnr)!=7)
    pari_err(talker,"incorrect bigray field");
  (void)checkbnf(gel(bnr,1));
}

void
checkbnrgen(GEN bnr)
{
  checkbnr(bnr);
  if (lg(bnr[5])<=3)
    pari_err(talker,"please apply bnrinit(,,1) and not bnrinit(,)");
}

GEN
check_units(GEN BNF, char *f)
{
  GEN bnf = checkbnf(BNF), x = gel(bnf,8);
  if (lg(x) < 6 || lg(x[5]) != lg(bnf[3])) pari_err(talker,"missing units in %s", f);
  return gel(x,5);
}

void
checkid(GEN x, long N)
{
  if (typ(x)!=t_MAT) pari_err(talker,"incorrect ideal");
  if (lg(x) == 1 || lg(x[1]) != N+1)
    pari_err(talker,"incorrect matrix for ideal");
}

void
checkbid(GEN bid)
{
  if (typ(bid)!=t_VEC || lg(bid)!=6 || lg(bid[1])!=3)
    pari_err(talker,"incorrect bigideal");
}

void
checkprimeid(GEN id)
{
  if (typ(id) != t_VEC || lg(id) != 6 || typ(id[2]) != t_COL)
    pari_err(talker,"incorrect prime ideal");
}

static int
_check_ZX(GEN x)
{
  long k = lg(x)-1;
  for ( ; k>1; k--)
    if (typ(x[k])!=t_INT) return 0;
  return 1;
}

void
check_ZX(GEN x, char *s)
{
  if (! _check_ZX(x)) pari_err(talker,"polynomial not in Z[X] in %s",s);
}
void
check_ZXY(GEN x, char *s)
{
  long k = lg(x)-1;
  for ( ; k>1; k--) {
    GEN t = gel(x,k);
    switch(typ(t)) {
      case t_INT: break;
      case t_POL: if (_check_ZX(t)) break;
      /* fall through */
      default: pari_err(talker,"polynomial not in Z[X,Y] in %s",s);
    }
  }
}

GEN
checknfelt_mod(GEN nf, GEN x, char *s)
{
  if (!gequal(gel(x,1),gel(nf,1)))
    pari_err(talker,"incompatible modulus in %s:\n  mod = %Z,\n  nf  = %Z",
        s, x[1], nf[1]);
  return gel(x,2);
}

GEN
get_primeid(GEN x)
{
  if (typ(x) != t_VEC) return NULL;
  if (lg(x) == 3) x = gel(x,1);
  if (lg(x) != 6 || typ(x[3]) != t_INT) return NULL;
  return x;
}

GEN
get_bnf(GEN x, long *t)
{
  switch(typ(x))
  {
    case t_POL: *t = typ_POL;  return NULL;
    case t_QUAD: *t = typ_Q  ; return NULL;
    case t_VEC:
      switch(lg(x))
      {
        case 3:
          if (typ(x[2]) != t_POLMOD) break;
          return get_bnf(gel(x,1),t);
        case 5 : *t = typ_QUA; return NULL;
        case 6 :
          if (typ(x[1]) != t_VEC || typ(x[3]) != t_MAT) break;
          *t = typ_BID; return NULL;
        case 10: *t = typ_NF; return NULL;
        case 11: *t = typ_BNF; return x;
        case 7 : *t = typ_BNR;
          x = gel(x,1); if (typ(x)!=t_VEC || lg(x)!=11) break;
          return x;
      }
    case t_MAT:
      if (lg(x)==2)
        switch(lg(x[1]))
        {
          case 7: case 10:
            *t = typ_CLA; return NULL;
        }
  }
  *t = typ_NULL; return NULL;
}

GEN
get_nf(GEN x, long *t)
{
  switch(typ(x))
  {
    case t_POL : *t = typ_POL; return NULL;
    case t_QUAD: *t = typ_Q  ; return NULL;
    case t_VEC:
      switch(lg(x))
      {
        case 3:
          if (typ(x[2]) != t_POLMOD) break;
          return get_nf(gel(x,1),t);
        case 10: *t = typ_NF; return x;
        case 11: *t = typ_BNF;
          x = gel(x,7); if (typ(x)!=t_VEC || lg(x)!=10) break;
          return x;
        case 7 : *t = typ_BNR;
          x = gel(x,1); if (typ(x)!=t_VEC || lg(x)!=11) break;
          x = gel(x,7); if (typ(x)!=t_VEC || lg(x)!=10) break;
          return x;
        case 6 :
          if (typ(x[1]) != t_VEC || typ(x[3]) != t_MAT) break;
          *t = typ_BID; return NULL;
        case 9 :
          x = gel(x,2);
          if (typ(x) == t_VEC && lg(x) == 4) *t = typ_GAL;
	  return NULL;
        case 14: case 20:
          *t = typ_ELL; return NULL;
      }break;
    case t_MAT:
      if (lg(x)==2)
        switch(lg(x[1]))
        {
          case 7: case 10:
            *t = typ_CLA; return NULL;
        }
  }
  *t = typ_NULL; return NULL;
}

/*************************************************************************/
/**									**/
/**			       GALOIS GROUP   				**/
/**									**/
/*************************************************************************/

/* exchange elements i and j in vector x */
static GEN
transroot(GEN x, int i, int j)
{
  long k;
  x = shallowcopy(x);
  k=x[i]; x[i]=x[j]; x[j]=k; return x;
}

GEN
tschirnhaus(GEN x)
{
  pari_sp av = avma, av2;
  long a, v = varn(x);
  GEN u, p1 = cgetg(5,t_POL);

  if (typ(x)!=t_POL) pari_err(notpoler,"tschirnhaus");
  if (lg(x) < 4) pari_err(constpoler,"tschirnhaus");
  if (v) { u=shallowcopy(x); setvarn(u,0); x=u; }
  p1[1] = evalsigne(1)|evalvarn(0);
  do
  {
    a = random_bits(2); if (a==0) a  = 1; gel(p1,4) = stoi(a);
    a = random_bits(3); if (a>=4) a -= 8; gel(p1,3) = stoi(a);
    a = random_bits(3); if (a>=4) a -= 8; gel(p1,2) = stoi(a);
    u = caract2(x,p1,v); av2 = avma;
  }
  while (degpol(srgcd(u,derivpol(u)))); /* while u not separable */
  if (DEBUGLEVEL>1)
    fprintferr("Tschirnhaus transform. New pol: %Z",u);
  avma=av2; return gerepileupto(av,u);
}

int
gpolcomp(GEN p1, GEN p2)
{
  long j = lg(p1)-2;
  int s;

  if (lg(p2)-2 != j)
    pari_err(bugparier,"gpolcomp (different degrees)");
  for (; j>=2; j--)
  {
    s = absi_cmp(gel(p1,j), gel(p2,j));
    if (s) return s;
  }
  return 0;
}

/* assume pol in Z[X]. Find C, L in Z such that POL = C pol(x/L) monic in Z[X].
 * Return POL and set *ptlead = L */
GEN
primitive_pol_to_monic(GEN pol, GEN *ptlead)
{
  long i,j,n = degpol(pol);
  GEN lead,fa,e,a, POL = shallowcopy(pol);

  a = POL + 2; lead = gel(a,n);
  if (signe(lead) < 0) { POL = gneg_i(POL); a = POL+2; lead = negi(lead); }
  if (is_pm1(lead)) { if (ptlead) *ptlead = NULL; return POL; }
  fa = auxdecomp(lead,0); lead = gen_1;
  e = gel(fa,2); fa = gel(fa,1);
  for (i=lg(e)-1; i>0;i--) e[i] = itos(gel(e,i));
  for (i=lg(fa)-1; i>0; i--)
  {
    GEN p = gel(fa,i), pk, pku;
    long k = (long)ceil((double)e[i] / n);
    long d = k * n - e[i], v, j0;
    /* find d, k such that  p^d pol(x / p^k) monic */
    for (j=n-1; j>0; j--)
    {
      if (!signe(a[j])) continue;
      v = Z_pval(gel(a,j), p);
      while (v + d < k * j) { k++; d += n; }
    }
    pk = powiu(p,k); j0 = d/k;
    pku = powiu(p,d - k*j0);
    /* a[j] *= p^(d - kj) */
    for (j=j0; j>=0; j--)
    {
      if (j < j0) pku = mulii(pku, pk);
      gel(a,j) = mulii(gel(a,j), pku);
    }
    j0++;
    pku = powiu(p,k*j0 - d);
    /* a[j] /= p^(kj - d) */
    for (j=j0; j<=n; j++)
    {
      if (j > j0) pku = mulii(pku, pk);
      gel(a,j) = diviiexact(gel(a,j), pku);
    }
    lead = mulii(lead, pk);
  }
  if (ptlead) *ptlead = lead; return POL;
}

/* x1*x2^2 + x2*x3^2 + x3*x4^2 + x4*x1^2 */
static GEN
F4(GEN x)
{
  return gadd(
    gmul(gel(x,1), gadd(gsqr(gel(x,2)), gmul(gel(x,4),gel(x,1)))),
    gmul(gel(x,3), gadd(gsqr(gel(x,4)), gmul(gel(x,2),gel(x,3)))));
}

static GEN
roots_to_ZX(GEN z, long *e)
{
  GEN a = roots_to_pol(z,0);
  GEN b = grndtoi(real_i(a),e);
  long e1 = gexpo(imag_i(a));
  if (e1 > *e) *e = e1;
  return b;
}

GEN
polgaloisnames(long a, long b)
{
  const char * const t[]={"S1", "S2", "A3", "S3", 
       "C(4) = 4", "E(4) = 2[x]2", "D(4)", "A4", "S4",
       "C(5) = 5", "D(5) = 5:2", "F(5) = 5:4", "A5", "S5",
       "C(6) = 6 = 3[x]2", "D_6(6) = [3]2", "D(6) = S(3)[x]2",
             "A_4(6) = [2^2]3", "F_18(6) = [3^2]2 = 3 wr 2", 
             "2A_4(6) = [2^3]3 = 2 wr 3", "S_4(6d) = [2^2]S(3)", 
             "S_4(6c) = 1/2[2^3]S(3)", "F_18(6):2 = [1/2.S(3)^2]2", 
             "F_36(6) = 1/2[S(3)^2]2", "2S_4(6) = [2^3]S(3) = 2 wr S(3)", 
             "L(6) = PSL(2,5) = A_5(6)", "F_36(6):2 = [S(3)^2]2 = S(3) wr 2", 
             "L(6):2 = PGL(2,5) = S_5(6)", "A6", "S6",
       "C(7) = 7", "D(7) = 7:2", "F_21(7) = 7:3", "F_42(7) = 7:6", 
             "L(7) = L(3,2)", "A7", "S7"};

   const long idx[]={0,1,2,4,9,14,30};
   return strtoGENstr(t[idx[a-1]+b-1]); 
}

static GEN
galois_res(long d, long n, long s, long k)
{
  long kk = k;
  GEN z = cgetg(5,t_VEC);
  if (!new_galois_format)
  {
    switch (d) {
      case 6:
        kk = (k == 6 || k == 2)? 2: 1;
        break;
      case 3:
        kk = (k == 2)? 1: 2;
        break;
      default:
        kk = 1;
    }
  }
  gel(z,1) = stoi(n);
  gel(z,2) = stoi(s);
  gel(z,3) = stoi(kk);
  gel(z,4) = polgaloisnames(d,k);
  return z;
}

GEN
polgalois(GEN x, long prec)
{
  pari_sp av = avma, av1;
  long i,j,k,n,f,l,l2,e,e1,pr,ind;
  GEN x1,p1,p2,p3,p4,p5,w,z,ee;
  static int ind5[20]={2,5,3,4, 1,3,4,5, 1,5,2,4, 1,2,3,5, 1,4,2,3};
  static int ind6[60]={3,5,4,6, 2,6,4,5, 2,3,5,6, 2,4,3,6, 2,5,3,4,
                       1,4,5,6, 1,5,3,6, 1,6,3,4, 1,3,4,5, 1,6,2,5,
                       1,2,4,6, 1,5,2,4, 1,3,2,6, 1,2,3,5, 1,4,2,3};
  if (typ(x)!=t_POL) pari_err(notpoler,"galois");
  n=degpol(x); if (n<=0) pari_err(constpoler,"galois");
  if (n>11) pari_err(impl,"galois of degree higher than 11");
  x = primpart(x);
  check_ZX(x, "galois");
  if (gisirreducible(x) != gen_1)
    pari_err(impl,"galois of reducible polynomial");

  if (n<4)
  {
    if (n == 1) { avma = av; return galois_res(n,1, 1,1); }
    if (n == 2) { avma = av; return galois_res(n,2,-1,1); }
    /* n = 3 */
    f = Z_issquare(ZX_disc(x));
    avma = av;
    return f? galois_res(n,3,1,1): 
              galois_res(n,6,-1,2);
  }
  x1 = x = primitive_pol_to_monic(x,NULL); av1=avma;
  if (n > 7) return galoisbig(x, prec);
  for(;;)
  {
    double cb = cauchy_bound(x);
    switch(n)
    {
      case 4: z = cgetg(7,t_VEC);
        prec = DEFAULTPREC + (long)(cb*(18./ LOG2 / BITS_IN_LONG));
        for(;;)
	{
	  p1=cleanroots(x,prec);
          gel(z,1) = F4(p1);
          gel(z,2) = F4(transroot(p1,1,2));
          gel(z,3) = F4(transroot(p1,1,3));
          gel(z,4) = F4(transroot(p1,1,4));
          gel(z,5) = F4(transroot(p1,2,3));
          gel(z,6) = F4(transroot(p1,3,4));
          p5 = roots_to_ZX(z, &e); if (e <= -10) break;
	  prec = (prec<<1)-2;
	}
        if (!ZX_is_squarefree(p5)) goto tchi;
	p2 = (GEN)factor(p5)[1];
	switch(lg(p2)-1)
	{
	  case 1: f = Z_issquare(ZX_disc(x)); avma = av;
            return f? galois_res(n,12,1,4): galois_res(n,24,-1,5);

	  case 2: avma = av; return galois_res(n,8,-1,3);
	
	  case 3: avma = av;
	    return (degpol(p2[1])==2)? galois_res(n,4,1,2): galois_res(n,4,-1,1);

	  default: pari_err(bugparier,"galois (bug1)");
	}

      case 5: z = cgetg(7,t_VEC);
        ee= cgetg(7,t_VECSMALL);
        w = cgetg(7,t_VECSMALL);
        prec = DEFAULTPREC + (long)(cb*(21. / LOG2/ BITS_IN_LONG));
        for(;;)
	{
          for(;;)
	  {
	    p1=cleanroots(x,prec);
	    for (l=1; l<=6; l++)
	    {
	      p2=(l==1)?p1: ((l<6)?transroot(p1,1,l): transroot(p1,2,5));
	      p3=gen_0;
              for (k=0,i=1; i<=5; i++,k+=4)
	      {
		p5 = gadd(gmul(gel(p2,ind5[k]),gel(p2,ind5[k+1])),
		          gmul(gel(p2,ind5[k+2]),gel(p2,ind5[k+3])));
		p3 = gadd(p3, gmul(gsqr(gel(p2,i)),p5));
	      }
	      gel(w,l) = grndtoi(real_i(p3),&e);
              e1 = gexpo(imag_i(p3)); if (e1>e) e=e1;
              ee[l]=e; gel(z,l) = p3;
	    }
            p5 = roots_to_ZX(z, &e); if (e <= -10) break;
	    prec = (prec<<1)-2;
	  }
          if (!ZX_is_squarefree(p5)) goto tchi;
	  p3=(GEN)factor(p5)[1];
	  f=Z_issquare(ZX_disc(x));
	  if (lg(p3)-1==1)
	  {
	    avma = av;
            return f? galois_res(n,60,1,4): galois_res(n,120,-1,5);
	  }
	  if (!f) { avma = av; return galois_res(n,20,-1,3); }

          pr = - (bit_accuracy(prec) >> 1);
          for (l=1; l<=6; l++)
	    if (ee[l] <= pr && gcmp0(poleval(p5,gel(w,l)))) break;
	  if (l>6) pari_err(bugparier,"galois (bug4)");
	  p2=(l==6)? transroot(p1,2,5):transroot(p1,1,l);
	  p3=gen_0;
	  for (i=1; i<=5; i++)
	  {
	    j = i+1; if (j>5) j -= 5;
	    p3 = gadd(p3,gmul(gmul(gel(p2,i),gel(p2,j)),
		              gsub(gel(p2,j),gel(p2,i))));
	  }
	  p5=gsqr(p3); p4=grndtoi(real_i(p5),&e);
          e1 = gexpo(imag_i(p5)); if (e1>e) e=e1;
	  if (e <= -10)
	  {
	    if (gcmp0(p4)) goto tchi;
	    f = Z_issquare(p4); avma = av;
            return f? galois_res(n,5,1,1): galois_res(n,10,1,2);
	  }
	  prec=(prec<<1)-2;
	}

      case 6: z = cgetg(7, t_VEC);
        prec = DEFAULTPREC + (long)(cb * (42. / LOG2 / BITS_IN_LONG));
        for(;;)
	{
          for(;;)
	  {
	    p1=cleanroots(x,prec);
	    for (l=1; l<=6; l++)
	    {
	      p2=(l==1)?p1:transroot(p1,1,l);
	      p3=gen_0; k=0;
              for (i=1; i<=5; i++) for (j=i+1; j<=6; j++)
	      {
		p5=gadd(gmul(gel(p2,ind6[k]),gel(p2,ind6[k+1])),
		        gmul(gel(p2,ind6[k+2]),gel(p2,ind6[k+3])));
		p3=gadd(p3,gmul(gsqr(gmul(gel(p2,i),gel(p2,j))),p5));
                k += 4;
	      }
	      gel(z,l) = p3;
	    }
            p5 = roots_to_ZX(z, &e); if (e <= -10) break;
	    prec=(prec<<1)-2;
	  }
	  if (!ZX_is_squarefree(p5)) goto tchi;
	  p2=(GEN)factor(p5)[1];
	  switch(lg(p2)-1)
	  {
	    case 1:
              z = cgetg(11,t_VEC); ind=0;
	      p3=gadd(gmul(gmul(gel(p1,1),gel(p1,2)),gel(p1,3)),
	              gmul(gmul(gel(p1,4),gel(p1,5)),gel(p1,6)));
              gel(z,++ind) = p3;
	      for (i=1; i<=3; i++)
		for (j=4; j<=6; j++)
		{
		  p2=transroot(p1,i,j);
		  p3=gadd(gmul(gmul(gel(p2,1),gel(p2,2)),gel(p2,3)),
		          gmul(gmul(gel(p2,4),gel(p2,5)),gel(p2,6)));
		  gel(z,++ind) = p3;
		}
              p5 = roots_to_ZX(z, &e);
	      if (e <= -10)
	      {
                if (!ZX_is_squarefree(p5)) goto tchi;
		p2 = (GEN)factor(p5)[1];
		f = Z_issquare(ZX_disc(x));
		avma = av;
		if (lg(p2)-1==1)
                  return f? galois_res(n,360,1,15): galois_res(n,720,-1,16);
		else
                  return f? galois_res(n,36,1,10): galois_res(n,72,-1,13);
	      }
	      prec=(prec<<1)-2; break;
		
	    case 2: l2=degpol(p2[1]); if (l2>3) l2=6-l2;
	      switch(l2)
	      {
		case 1: f = Z_issquare(ZX_disc(x));
		  avma = av;
                  return f? galois_res(n,60,1,12): galois_res(n,120,-1,14);
		case 2: f = Z_issquare(ZX_disc(x));
		  if (f) { avma = av; return galois_res(n,24,1,7); }
                  p3 = (degpol(p2[1])==2)? gel(p2,2): gel(p2,1);
                  f = Z_issquare(ZX_disc(p3));
                  avma = av;
                  return f? galois_res(n,24,-1,6): galois_res(n,48,-1,11);
		case 3: f = Z_issquare(ZX_disc(gel(p2,1)))
		         || Z_issquare(ZX_disc(gel(p2,2)));
		  avma = av;
                  return f? galois_res(n,18,-1,5): galois_res(n,36,-1,9);
	      }
	    case 3:
	      for (l2=1; l2<=3; l2++)
		if (degpol(p2[l2]) >= 3) p3 = gel(p2,l2);
	      if (degpol(p3) == 3)
	      {
		f = Z_issquare(ZX_disc(p3)); avma = av;
                return f? galois_res(n,6,-1,1): galois_res(n,12,-1,3);
	      }
	      else
	      {
		f = Z_issquare(ZX_disc(x)); avma = av;
                return f? galois_res(n,12,1,4): galois_res(n,24,-1,8);
	      }
	    case 4: avma = av; return galois_res(n,6,-1,2);
            default: pari_err(bugparier,"galois (bug3)");
	  }
	}
	
      case 7: z = cgetg(36,t_VEC);
        prec = DEFAULTPREC + (long)(cb*(7. / LOG2 / BITS_IN_LONG));
        for(;;)
	{
	  ind = 0; p1=cleanroots(x,prec);
	  for (i=1; i<=5; i++)
	    for (j=i+1; j<=6; j++)
            {
              GEN t = gadd(gel(p1,i),gel(p1,j));
	      for (k=j+1; k<=7; k++) gel(z,++ind) = gadd(t, gel(p1,k));
            }
          p5 = roots_to_ZX(z, &e); if (e <= -10) break;
          prec = (prec<<1)-2;
	}
        if (!ZX_is_squarefree(p5)) goto tchi;
	p2=(GEN)factor(p5)[1];
	switch(lg(p2)-1)
	{
	  case 1: f = Z_issquare(ZX_disc(x)); avma = av;
            return f? galois_res(n,2520,1,6): galois_res(n,5040,-1,7);
	  case 2: f = degpol(p2[1]); avma = av;
	    return (f==7 || f==28)? galois_res(n,168,1,5): galois_res(n,42,-1,4);
	  case 3: avma = av; return galois_res(n,21,1,3);
	  case 4: avma = av; return galois_res(n,14,-1,2);
	  case 5: avma = av; return galois_res(n,7,1,1);
          default: pari_err(bugparier,"galois (bug2)");
	}
    }
    tchi: avma = av1; x = tschirnhaus(x1);
  }
}

#undef _res

/*Deprecated*/
GEN galois(GEN x, long prec) {return polgalois(x,prec);}

GEN
galoisapply(GEN nf, GEN aut, GEN x)
{
  pari_sp av = avma;
  long lx, j, N;
  GEN p, p1, y, pol;

  nf=checknf(nf); pol=gel(nf,1);
  if (typ(aut)==t_POL) aut = gmodulo(aut,pol);
  else
  {
    if (typ(aut)!=t_POLMOD || !gequal(gel(aut,1),pol))
      pari_err(talker,"incorrect galois automorphism in galoisapply");
  }
  switch(typ(x))
  {
    case t_INT: case t_INTMOD: case t_FRAC: case t_PADIC:
      avma=av; return gcopy(x);

    case t_POLMOD: x = gel(x,2); /* fall through */
    case t_POL:
      p1 = gsubst(x,varn(pol),aut);
      if (typ(p1)!=t_POLMOD || !gequal(gel(p1,1),pol)) p1 = gmodulo(p1,pol);
      return gerepileupto(av,p1);

    case t_VEC:
      if (lg(x)==3)
      {
	y=cgetg(3,t_VEC);
	gel(y,1) = galoisapply(nf,aut,gel(x,1));
        gel(y,2) = gcopy(gel(x,2)); return gerepileupto(av, y);
      }
      if (lg(x)!=6) pari_err(typeer,"galoisapply");
      y=cgetg(6,t_VEC); y[1]=x[1]; y[3]=x[3]; y[4]=x[4];
      p = gel(x,1);
      p1=centermod(galoisapply(nf,aut,gel(x,2)), p);
      if (is_pm1(x[3]))
	if (Z_pval(subres(coltoliftalg(nf,p1),pol), p) > itos(gel(x,4)))
	  gel(p1,1) =  (signe(p1[1]) > 0)? subii(gel(p1,1), p)
	                                 : addii(gel(p1,1), p);
      gel(y,2) = p1;
      gel(y,5) = centermod(galoisapply(nf,aut,gel(x,5)), p);
      return gerepilecopy(av,y);

    case t_COL:
      N = degpol(pol);
      if (lg(x)!=N+1) pari_err(typeer,"galoisapply");
      p1 = gsubst(coltoliftalg(nf,x), varn(pol), aut);
      return gerepileupto(av, algtobasis(nf,p1));

    case t_MAT:
      lx=lg(x); if (lx==1) return cgetg(1,t_MAT);
      N=degpol(pol);
      if (lg(x[1])!=N+1) pari_err(typeer,"galoisapply");
      p1=cgetg(lx,t_MAT);
      for (j=1; j<lx; j++) gel(p1,j) = galoisapply(nf,aut,gel(x,j));
      if (lx==N+1) p1 = idealhermite(nf,p1);
      return gerepileupto(av,p1);
  }
  pari_err(typeer,"galoisapply");
  return NULL; /* not reached */
}

GEN pol_to_monic(GEN pol, GEN *lead);
int cmp_pol(GEN x, GEN y);

/* x = relative polynomial nf = absolute nf, bnf = absolute bnf */
GEN
get_bnfpol(GEN x, GEN *bnf, GEN *nf)
{
  *bnf = checkbnf_i(x);
  *nf  = checknf_i(x);
  if (*nf) x = gel(*nf, 1);
  if (typ(x) != t_POL) pari_err(typeer,"get_bnfpol");
  return x;
}

GEN
get_nfpol(GEN x, GEN *nf)
{
  if (typ(x) == t_POL) { *nf = NULL; return x; }
  *nf = checknf(x); return gel(*nf,1);
}

/* if fliso test for isomorphism, for inclusion otherwise. */
static GEN
nfiso0(GEN a, GEN b, long fliso)
{
  pari_sp av = avma;
  long n,m,i,vb,lx;
  GEN nfa, nfb, p1, y, la, lb;

  a = primpart(get_nfpol(a, &nfa)); check_ZX(a, "nsiso0");
  b = primpart(get_nfpol(b, &nfb)); check_ZX(b, "nsiso0");
  if (fliso && nfa && !nfb) { p1=a; a=b; b=p1; p1=nfa; nfa=nfb; nfb=p1; }
  m=degpol(a);
  n=degpol(b); if (m<=0 || n<=0) pari_err(constpoler,"nfiso or nfincl");
  if (fliso)
    { if (n!=m) return gen_0; }
  else
    { if (n%m) return gen_0; }

  if (nfb) lb = NULL; else b = pol_to_monic(b,&lb);
  if (nfa) la = NULL; else a = pol_to_monic(a,&la);
  if (nfa && nfb)
  {
    if (fliso)
    {
      if (!gequal(gel(nfa,2),gel(nfb,2))
       || !gequal(gel(nfa,3),gel(nfb,3))) return gen_0;
    }
    else
      if (!dvdii(gel(nfb,3), powiu(gel(nfa,3),n/m))) return gen_0;
  }
  else
  {
    GEN da = nfa? gel(nfa,3): ZX_disc(a);
    GEN db = nfb? gel(nfb,3): ZX_disc(b);
    if (fliso)
    {
      if (gissquare(gdiv(da,db)) == gen_0) { avma=av; return gen_0; }
    }
    else
    {
      long q=n/m;
      GEN fa=Z_factor(da), ex=gel(fa,2);
      fa=gel(fa,1); lx=lg(fa);
      for (i=1; i<lx; i++)
        if (mod2(gel(ex,i)) && !dvdii(db, powiu(gel(fa,i),q)))
          { avma=av; return gen_0; }
    }
  }
  a = shallowcopy(a); setvarn(a,0);
  b = shallowcopy(b); vb=varn(b);
  if (nfb)
  {
    if (vb == 0) nfb = gsubst(nfb, 0, pol_x[MAXVARN]);
    y = lift_intern(nfroots(nfb,a));
  }
  else
  {
    if (vb == 0) setvarn(b, fetch_var());
    y = (GEN)polfnf(a,b)[1]; lx = lg(y);
    for (i=1; i<lx; i++)
    {
      if (lg(y[i]) != 4) { setlg(y,i); break; }
      gel(y,i) = gneg_i(lift_intern(gmael(y,i,2)));
    }
    y = gen_sort(y, 0, cmp_pol);
    settyp(y, t_VEC);
    if (vb == 0) (void)delete_var();
  }
  lx = lg(y); if (lx==1) { avma=av; return gen_0; }
  for (i=1; i<lx; i++)
  {
    p1 = gel(y,i);
    if (typ(p1) == t_POL) setvarn(p1, vb); else p1 = scalarpol(p1, vb);
    if (lb) p1 = poleval(p1, monomial(lb, 1, vb));
    gel(y,i) = la? gdiv(p1,la): p1;
  }
  return gerepilecopy(av,y);
}

GEN
nfisisom(GEN a, GEN b) { return nfiso0(a,b,1); }

GEN
nfisincl(GEN a, GEN b) { return nfiso0(a,b,0); }

/*************************************************************************/
/**									**/
/**			       INITALG					**/
/**									**/
/*************************************************************************/

GEN
get_roots(GEN x,long r1,long prec)
{
  GEN roo = (typ(x)!=t_POL)? shallowcopy(x): roots(x,prec);
  long i, ru = (lg(roo)-1 + r1) >> 1;

  for (i=1; i<=r1; i++) gel(roo,i) = real_i(gel(roo,i));
  for (   ; i<=ru; i++) roo[i] = roo[(i<<1)-r1];
  roo[0]=evaltyp(t_VEC)|evallg(ru+1); return roo;
}

/* For internal use. compute trace(x mod pol), sym=polsym(pol,deg(pol)-1) */
GEN
quicktrace(GEN x, GEN sym)
{
  GEN p1 = gen_0;
  long i;

  if (typ(x) != t_POL) return gmul(x, gel(sym,1));
  if (signe(x))
  {
    sym--;
    for (i=lg(x)-1; i>1; i--)
      p1 = gadd(p1, gmul(gel(x,i),gel(sym,i)));
  }
  return p1;
}

/* T = trace(w[i]), x = sum x[i] w[i], x[i] integer */
static GEN
trace_col(GEN x, GEN T)
{
  pari_sp av = avma;
  GEN t = gen_0;
  long i, l = lg(x);

  t = mulii(gel(x,1),gel(T,1));
  for (i=2; i<l; i++)
    if (signe(x[i])) t = addii(t, mulii(gel(x,i),gel(T,i)));
  return gerepileuptoint(av, t);
}

/* pol belonging to Z[x], return a monic polynomial generating the same field
 * as pol (x-> ax+b)) set lead = NULL if pol was monic (after dividing
 * by the content), and to to leading coefficient otherwise.
 * No garbage collecting done.
 */
GEN
pol_to_monic(GEN pol, GEN *lead)
{
  long n = lg(pol)-1;

  if (n==1 || gcmp1(gel(pol,n))) { *lead = NULL; return pol; }
  return primitive_pol_to_monic(primpart(pol), lead);
}

/* assume lg(nf) > 3 && typ(nf) = container [hopefully a genuine nf] */
long
nf_get_r1(GEN nf)
{
  GEN x = gel(nf,2);
  if (typ(x) != t_VEC || lg(x) != 3 || typ(x[1]) != t_INT)
    pari_err(talker,"false nf in nf_get_r1");
  return itos(gel(x,1));
}
long
nf_get_r2(GEN nf)
{
  GEN x = gel(nf,2);
  if (typ(x) != t_VEC || lg(x) != 3 || typ(x[2]) != t_INT)
    pari_err(talker,"false nf in nf_get_r2");
  return itos(gel(x,2));
}
void
nf_get_sign(GEN nf, long *r1, long *r2)
{
  GEN x = gel(nf,2);
  if (typ(x) != t_VEC || lg(x) != 3
      || typ(x[1]) != t_INT || typ(x[2]) != t_INT)
    pari_err(talker,"false nf in nf_get_sign");
  *r1 = itos(gel(x,1));
  *r2 = (degpol(nf[1]) - *r1) >> 1;
}

static GEN
get_Tr(GEN mul, GEN x, GEN basden)
{
  GEN tr,T,T1,sym, bas = gel(basden,1), den = gel(basden,2);
  long i,j,n = lg(bas)-1;
  T = cgetg(n+1,t_MAT);
  T1 = cgetg(n+1,t_COL);
  sym = polsym(x, n-1);

  gel(T1,1) = utoipos(n);
  for (i=2; i<=n; i++)
  {
    tr = quicktrace(gel(bas,i), sym);
    if (den && den[i]) tr = diviiexact(tr,gel(den,i));
    gel(T1,i) = tr; /* tr(w[i]) */
  }
  gel(T,1) = T1;
  for (i=2; i<=n; i++)
  {
    gel(T,i) = cgetg(n+1,t_COL); gcoeff(T,1,i) = gel(T1,i);
    for (j=2; j<=i; j++)
      gcoeff(T,i,j) = gcoeff(T,j,i) = trace_col(gel(mul,j+(i-1)*n), T1);
  }
  return T;
}

/* return [bas[i]*denom(bas[i]), denom(bas[i])], denom 1 is given as NULL */
GEN
get_bas_den(GEN bas)
{
  GEN b,d,den, dbas = shallowcopy(bas);
  long i, l = lg(bas);
  int power = 1;
  den = cgetg(l,t_VEC);
  for (i=1; i<l; i++)
  {
    b = Q_remove_denom(gel(bas,i), &d);
    gel(dbas,i) = b;
    gel(den,i) = d; if (d) power = 0;
  }
  if (power) den = NULL; /* power basis */
  return mkvec2(dbas, den);
}

/* allow x or y = NULL (act as 1) */
static GEN
_mulii(GEN x, GEN y)
{
  if (!x) return y;
  if (!y) return x;
  return mulii(x,y);
}

GEN
get_mul_table(GEN x,GEN basden,GEN invbas)
{
  long i,j, n = degpol(x);
  GEN z, d, bas, den, mul = cgetg(n*n+1,t_MAT);
  
  if (typ(basden[1]) != t_VEC) basden = get_bas_den(basden); /*integral basis*/
  bas = gel(basden,1);
  den = gel(basden,2);
  for (i=1; i<=n; i++)
    for (j=i; j<=n; j++)
    {
      pari_sp av = avma;
      z = grem(gmul(gel(bas,j),gel(bas,i)), x);
      z = mulmat_pol(invbas, z); /* integral column */
      if (den)
      {
        d = _mulii(gel(den,i), gel(den,j));
        if (d) z = gdivexact(z, d);
      }
      gel(mul,j+(i-1)*n) = gel(mul,i+(j-1)*n) = gerepileupto(av,z);
    }
  return mul;
}

/* as get_Tr, mul_table not precomputed */
static GEN
make_Tr(GEN x, GEN w)
{
  long i,j, n = degpol(x);
  GEN p1,p2,t,d;
  GEN sym = cgetg(n+2,t_VEC);
  GEN den = cgetg(n+1,t_VEC);
  GEN T = cgetg(n+1,t_MAT);
  pari_sp av;

  sym = polsym(x, n-1);
  p1 = get_bas_den(w);
  w   = gel(p1,1);
  den = gel(p1,2);
  for (i=1; i<=n; i++)
  {
    p1 = cgetg(n+1,t_COL); gel(T,i) = p1;
    for (j=1; j<i ; j++) gel(p1,j) = gcoeff(T,i,j);
    for (   ; j<=n; j++)
    {
      av = avma;
      p2 = grem(gmul(gel(w,i),gel(w,j)), x);
      t = quicktrace(p2, sym);
      if (den)
      {
        d = _mulii(gel(den,i),gel(den,j));
        if (d) t = diviiexact(t, d);
      }
      gel(p1,j) = gerepileuptoint(av, t);
    }
  }
  return T;
}

/* compute roots so that _absolute_ accuracy of M >= prec [also holds for G] */
static void
get_roots_for_M(nffp_t *F)
{
  long n, eBD, PREC;

  if (F->extraprec < 0)
  { /* not initialized yet */
    double er;
    n = degpol(F->x);
    eBD = 1 + gexpo(gel(F->basden,1));
    er  = F->ro? (1+gexpo(F->ro)): cauchy_bound(F->x)/LOG2;
    if (er < 0) er = 0;
    F->extraprec = (long)((n*er + eBD + log2(n)) / BITS_IN_LONG);
  }

  PREC = F->prec + F->extraprec;
  if (F->ro && gprecision(gel(F->ro,1)) >= PREC) return;
  F->ro = get_roots(F->x, F->r1, PREC);
}

/* [bas[i]/den[i]]= integer basis. roo = real part of the roots */
static void
make_M(nffp_t *F, int trunc)
{
  GEN bas = gel(F->basden,1), den = gel(F->basden,2), ro = F->ro;
  GEN m, d, M;
  long i, j, l = lg(ro), n = lg(bas);
  M = cgetg(n,t_MAT);
    m = cgetg(l, t_COL); gel(M,1) = m;
    for (i=1; i<l; i++) gel(m,i) = gen_1; /* bas[1] = 1 */
  for (j=2; j<n; j++)
  {
    m = cgetg(l,t_COL); gel(M,j) = m;
    for (i=1; i<l; i++) gel(m,i) = poleval(gel(bas,j), gel(ro,i));
  }
  if (den)
  {
    GEN invd, rd = cgetr(F->prec + F->extraprec);
    for (j=2; j<n; j++)
    {
      d = gel(den,j); if (!d) continue;
      m = gel(M,j); affir(d,rd); invd = ginv(rd);
      for (i=1; i<l; i++) gel(m,i) = gmul(gel(m,i), invd);
    }
  }

  if (trunc && gprecision(M) > F->prec)
  {
    M     = gprec_w(M, F->prec);
    F->ro = gprec_w(ro,F->prec);
  }
  if (DEBUGLEVEL>4) msgtimer("matrix M");
  F->M = M;
}

/* return G real such that G~ * G = T_2 */
static void
make_G(nffp_t *F)
{
  GEN G, M = F->M;
  long i, j, k, r1 = F->r1, l = lg(M);

  G = cgetg(l, t_MAT);
  for (j=1; j<l; j++)
  {
    GEN g = cgetg(l, t_COL);
    GEN m = gel(M,j);
    gel(G,j) = g;
    for (k=i=1; i<=r1; i++) g[k++] = m[i];
    for (     ; k < l; i++)
    {
      GEN r = gel(m,i);
      if (typ(r) == t_COMPLEX)
      {
        gel(g,k++) = mpadd(gel(r,1), gel(r,2));
        gel(g,k++) = mpsub(gel(r,1), gel(r,2));
      }
      else
      {
        gel(g,k++) = r;
        gel(g,k++) = r;
      }
    }
  }
  F->G = G;
}

static void
make_M_G(nffp_t *F, int trunc)
{
  get_roots_for_M(F);
  make_M(F, trunc);
  make_G(F);
}

void
remake_GM(GEN nf, nffp_t *F, long prec)
{
  F->x  = gel(nf,1);
  F->ro = NULL;
  F->r1 = nf_get_r1(nf);
  F->basden = get_bas_den(gel(nf,7));
  F->extraprec = -1;
  F->prec = prec; make_M_G(F, 1);
}

static void
nffp_init(nffp_t *F, nfbasic_t *T, GEN ro, long prec)
{
  F->x  = T->x;
  F->ro = ro;
  F->r1 = T->r1;
  if (!T->basden) T->basden = get_bas_den(T->bas);
  F->basden = T->basden;
  F->extraprec = -1;
  F->prec = prec;
}

static void
get_nf_fp_compo(nfbasic_t *T, nffp_t *F, GEN ro, long prec)
{
  nffp_init(F,T,ro,prec);
  make_M_G(F, 1);
}

static GEN
get_sign(long r1, long n) { return mkvec2s(r1, (n-r1)>>1); }

GEN
nfbasic_to_nf(nfbasic_t *T, GEN ro, long prec)
{
  GEN nf = cgetg(10,t_VEC);
  GEN x = T->x;
  GEN absdK, invbas, Tr, D, TI, A, dA, MDI, mat = cgetg(8,t_VEC);
  nffp_t F;
  get_nf_fp_compo(T, &F, ro, prec);

  gel(nf,1) = T->x;
  gel(nf,2) = get_sign(T->r1, degpol(T->x));
  gel(nf,3) = T->dK;
  gel(nf,4) = T->index;
  gel(nf,6) = F.ro;
  gel(nf,5) = mat;
  gel(nf,7) = T->bas;

  gel(mat,1) = F.M;
  gel(mat,2) = F.G;

  invbas = QM_inv(RgXV_to_RgM(T->bas, lg(T->bas)-1), gen_1);
  gel(nf,8) = invbas;
  gel(nf,9) = get_mul_table(x, F.basden, invbas);
  if (DEBUGLEVEL) msgtimer("mult. table");

  Tr = get_Tr(gel(nf,9), x, F.basden);
  absdK = T->dK; if (signe(absdK) < 0) absdK = negi(absdK);
  TI = ZM_inv(Tr, absdK); /* dK T^-1 */
  A = Q_primitive_part(TI, &dA);
  gel(mat,6) = A; /* primitive part of codifferent, dA its denominator */
  dA = dA? diviiexact(absdK, dA): absdK;
  A = hnfmodid(A, dA);
  MDI = ideal_two_elt(nf, A);
  gel(MDI,2) = eltmul_get_table(nf, gel(MDI,2));
  gel(mat,7) = MDI;
  if (is_pm1(T->index))
    D = idealhermite_aux(nf, derivpol(x));
  else
    D = gmul(dA, idealinv(nf, A));
  gel(mat,3) = gen_0; /* FIXME: was gram matrix of current mat[2]. Useless */
  gel(mat,4) = Tr;
  gel(mat,5) = D;
  if (DEBUGLEVEL) msgtimer("matrices");
  return nf;
}

static GEN
hnffromLLL(GEN nf)
{
  GEN d, x;
  x = RgXV_to_RgM(gel(nf,7), degpol(nf[1]));
  x = Q_remove_denom(x, &d);
  if (!d) return x; /* power basis */
  return gauss(hnfmodid(x, d), x);
}

static GEN
nfbasechange(GEN u, GEN x)
{
  long i,lx;
  GEN y;
  switch(typ(x))
  {
    case t_COL: /* nfelt */
      return gmul(u, x);

    case t_MAT: /* ideal */
      y = shallowcopy(x);
      lx = lg(x);
      for (i=1; i<lx; i++) gel(y,i) = gmul(u, gel(y,i));
      break;

    case t_VEC: /* pr */
      checkprimeid(x);
      y = shallowcopy(x);
      gel(y,2) = gmul(u, gel(y,2));
      gel(y,5) = gmul(u, gel(y,5));
      break;
    default: y = x;
  }
  return y;
}

GEN
nffromhnfbasis(GEN nf, GEN x)
{
  long tx = typ(x);
  pari_sp av = avma;
  GEN u;
  if (!is_vec_t(tx)) return gcopy(x);
  nf = checknf(nf);
  u = hnffromLLL(nf);
  return gerepilecopy(av, nfbasechange(u, x));
}

GEN
nftohnfbasis(GEN nf, GEN x)
{
  long tx = typ(x);
  pari_sp av = avma;
  GEN u;
  if (!is_vec_t(tx)) return gcopy(x);
  nf = checknf(nf);
  u = ZM_inv(hnffromLLL(nf), gen_1);
  return gerepilecopy(av, nfbasechange(u, x));
}

static GEN
get_red_G(nfbasic_t *T, GEN *pro)
{
  GEN G, u, u0 = NULL;
  pari_sp av;
  long i, prec, extraprec, n = degpol(T->x), MARKED = 1;
  nffp_t F;

  extraprec = (long) (0.25 * n / (sizeof(long) / 4));
  prec = DEFAULTPREC + extraprec;
  nffp_init(&F, T, *pro, prec);
  av = avma;
  for (i=1; ; i++)
  {
    F.prec = prec; make_M_G(&F, 0); G = F.G;
    if (u0) G = gmul(G, u0);
    if (DEBUGLEVEL)
      fprintferr("get_red_G: starting LLL, prec = %ld (%ld + %ld)\n",
                  prec + F.extraprec, prec, F.extraprec);
    if ((u = lllfp_marked(&MARKED, G, 100, 2, prec, 0)))
    {
      if (typ(u) == t_MAT) break;
      u = gel(u,1);
      if (u0) u0 = gerepileupto(av, gmul(u0,u));
      else    u0 = gerepilecopy(av, u);
    }
    prec = (prec<<1)-2 + (gexpo(u0) >> TWOPOTBITS_IN_LONG);
    F.ro = NULL;
    if (DEBUGLEVEL) pari_warn(warnprec,"get_red_G", prec);
  }
  *pro = F.ro; 
  if (u0) u = gmul(u0,u);
  if (MARKED != 1) lswap(u[1], u[MARKED]);
  return u;
}

/* Return the base change matrix giving an LLL-reduced basis for the
 * integer basis of nf(T).
 * *ro = roots of x, computed to precision prec [or NULL] */
static GEN
get_LLL_basis(nfbasic_t *T, GEN *pro)
{
  GEN u;
  if (T->r1 != degpol(T->x)) u = get_red_G(T, pro);
  else
  {
    long MARKED = 1;
    u = lllfp_marked(&MARKED, make_Tr(T->x, T->bas), 100, 0, DEFAULTPREC, 1);
    if (!u) u = matid(1);
    else
      if (MARKED != 1) lswap(u[1], u[MARKED]);
  }
  return u;
}

static void
set_LLL_basis(nfbasic_t *T, GEN *pro)
{
  T->bas = gmul(T->bas, get_LLL_basis(T, pro));
  T->basden = NULL; /* recompute */
  if (DEBUGLEVEL) msgtimer("LLL basis");
}

typedef struct {
  GEN xbest, dxbest;
  long ind, indmax, indbest;
} ok_pol_t;

/* is xn better than x ? */
static int
better_pol(GEN xn, GEN dxn, GEN x, GEN dx)
{
  int fl;
  if (!x) return 1;
  fl = absi_cmp(dxn, dx);
  return (fl < 0 || (!fl && gpolcomp(xn, x) < 0));
}

static GEN
ok_pol(void *TT, GEN xn)
{
  ok_pol_t *T = (ok_pol_t*)TT;
  GEN dxn;

  if (++T->ind > T->indmax && T->xbest) return T->xbest;

  if (!ZX_is_squarefree(xn)) return (T->ind == T->indmax)? T->xbest: NULL;
  if (DEBUGLEVEL>3) outerr(xn);
  dxn = ZX_disc(xn);
  if (better_pol(xn, dxn, T->xbest, T->dxbest))
  {
    T->dxbest = dxn; T->xbest = xn; T->indbest = T->ind;
  }
  if (T->ind >= T->indmax) return T->xbest;
  return NULL;
}

/* z in Z[X] with positive leading coefficient. Set z := z(-X) or z(X) so that
 * second coefficient is > 0. In place. */
static int
canon_pol(GEN z)
{
  long i,s;

  for (i=lg(z)-2; i>=2; i-=2)
  {
    s = signe(z[i]);
    if (s > 0)
    {
      for (; i>=2; i-=2) gel(z,i) = negi(gel(z,i));
      return -1;
    }
    if (s) return 1;
  }
  return 0;
}

static GEN _polred(GEN x, GEN a, GEN *pta, FP_chk_fun *CHECK);

/* Seek a simpler, polynomial pol defining the same number field as
 * x (assumed to be monic at this point) */
static GEN
nfpolred(int part, nfbasic_t *T)
{
  GEN x = T->x, dx = T->dx, a = T->bas;
  GEN phi, xbest, dxbest, mat, d, rev;
  long i, n = lg(a)-1, v = varn(x);
  ok_pol_t O;
  FP_chk_fun chk = { &ok_pol, NULL, NULL, NULL, 0 };

  if (degpol(x) == 1) { T->x = gsub(pol_x[v],gen_1); return gen_1; }

  if (!dx) dx = mulii(T->dK, sqri(T->index));

  O.ind    = 0;
  O.indmax = part? min(n,3): n;
  O.xbest  = NULL;
  chk.data = (void*)&O;
  if (!_polred(x, a, NULL, &chk))
    pari_err(talker,"you found a counter-example to a conjecture, please report!");
  xbest = O.xbest; dxbest = O.dxbest;
  if (!better_pol(xbest, dxbest, x, dx)) return NULL; /* no improvement */

  /* update T */
  phi = gel(a, O.indbest);
  if (canon_pol(xbest) < 0) phi = gneg_i(phi);
  if (DEBUGLEVEL>1) fprintferr("xbest = %Z\n",xbest);
  rev = modreverse_i(phi, x);
  for (i=1; i<=n; i++) gel(a,i) = RgX_RgXQ_compo(gel(a,i), rev, xbest);
  mat = RgXV_to_RgM(Q_remove_denom(a, &d), n);
  if (d) mat = gdiv(hnfmodid(mat,d), d); else mat = matid(n);

  (void)Z_issquarerem(diviiexact(dxbest,T->dK), &(T->index));
  T->bas= RgM_to_RgXV(mat, v);
  T->dx = dxbest;
  T->x  = xbest; return rev;
}

GEN
get_nfindex(GEN bas)
{
  pari_sp av = avma;
  long n = lg(bas)-1;
  GEN d, mat = RgXV_to_RgM(Q_remove_denom(bas, &d), n);
  if (!d) { avma = av; return gen_1; }
  return gerepileuptoint(av, diviiexact(powiu(d, n), det(mat)));
}

void
nfbasic_init(GEN x, long flag, GEN fa, nfbasic_t *T)
{
  GEN bas, dK, dx, index;
  long r1;

  T->basden = NULL;
  T->lead   = NULL;
  if (DEBUGLEVEL) (void)timer2();
  if (typ(x) == t_POL)
  {
    check_ZX(x, "nfinit");
    if (gisirreducible(x) == gen_0) pari_err(redpoler, "nfinit");
    x = pol_to_monic(x, &(T->lead));
    bas = allbase(x, flag, &dx, &dK, &index, &fa);
    if (DEBUGLEVEL) msgtimer("round4");
    r1 = sturm(x);
  }
  else if (typ(x) == t_VEC && lg(x)==3 
           && typ(x[1])==t_POL && lg(x[2])-1 == degpol(x[1]))
  { /* monic integral polynomial + integer basis */
    GEN mat;
    bas = gel(x,2); x = gel(x,1);
    if (typ(bas) == t_MAT)
      { mat = bas; bas = RgM_to_RgXV(mat,varn(x)); }
    else
        mat = RgXV_to_RgM(bas, lg(bas)-1);
    index = get_nfindex(bas);
    dx = ZX_disc(x);
    dK = diviiexact(dx, sqri(index));
    r1 = sturm(x);
  }
  else
  { /* nf, bnf, bnr */
    GEN nf = checknf(x);
    x     = gel(nf,1);
    dK    = gel(nf,3);
    index = gel(nf,4);
    bas   = gel(nf,7);
    r1 = nf_get_r1(nf); dx = NULL;
  }

  T->x     = x;
  T->r1    = r1;
  T->dx    = dx;
  T->dK    = dK;
  T->bas   = bas;
  T->index = index;
}

/* Initialize the number field defined by the polynomial x (in variable v)
 * flag & nf_RED:     try a polred first.
 * flag & nf_PARTRED: as nf_RED but check the first two elements only.
 * flag & nf_ORIG
 *    do a polred and return [nfinit(x), Mod(a,red)], where
 *    Mod(a,red) = Mod(v,x) (i.e return the base change). */
GEN
initalg_i(GEN x, long flag, long prec)
{
  const pari_sp av = avma;
  GEN nf, rev = NULL, ro = NULL;
  nfbasic_t T;

  nfbasic_init(x, flag, NULL, &T);
  set_LLL_basis(&T, &ro);

  if (T.lead && !(flag & (nf_RED|nf_PARTRED)))
  {
    pari_warn(warner,"non-monic polynomial. Result of the form [nf,c]");
    flag |= nf_PARTRED | nf_ORIG;
  }
  if (flag & (nf_RED|nf_PARTRED))
  {
    rev = nfpolred(flag & nf_PARTRED, &T);
    if (DEBUGLEVEL) msgtimer("polred");
    if (rev) { ro = NULL; set_LLL_basis(&T, &ro); } /* changed T.x */
    if (flag & nf_ORIG)
    {
      if (!rev) rev = pol_x[varn(T.x)]; /* no improvement */
      if (T.lead) rev = gdiv(rev, T.lead);
      rev = mkpolmod(rev, T.x);
    }
  }

  nf = nfbasic_to_nf(&T, ro, prec);
  if (flag & nf_ORIG) nf = mkvec2(nf, rev);
  return gerepilecopy(av, nf);
}

GEN
initalgred(GEN x, long prec)  { return initalg_i(x, nf_RED, prec); }
GEN
initalgred2(GEN x, long prec) { return initalg_i(x, nf_RED|nf_ORIG, prec); }
GEN
initalg(GEN x, long prec)     { return initalg_i(x, 0, prec); }

GEN
nfinit0(GEN x, long flag,long prec)
{
  switch(flag)
  {
    case 0:
    case 1: return initalg_i(x,0,prec);
    case 2: return initalg_i(x,nf_RED,prec);
    case 3: return initalg_i(x,nf_RED|nf_ORIG,prec);
    case 4: return initalg_i(x,nf_PARTRED,prec);
    case 5: return initalg_i(x,nf_PARTRED|nf_ORIG,prec);
    default: pari_err(flagerr,"nfinit");
  }
  return NULL; /* not reached */
}

/* assume x a bnr/bnf/nf */
long
nfgetprec(GEN x)
{
  GEN nf = checknf(x), ro = gel(nf,6);
  return (typ(ro)==t_VEC)? precision(gel(ro,1)): DEFAULTPREC;
}

/* assume nf is an nf */
GEN
nfnewprec_i(GEN nf, long prec)
{
  GEN NF = shallowcopy(nf);
  nffp_t F;
  gel(NF,5) = shallowcopy(gel(NF,5));
  remake_GM(NF, &F, prec);
  gel(NF,6) = F.ro;
  gmael(NF,5,1) = F.M;
  gmael(NF,5,2) = F.G;
  return NF;
}
static GEN
_nfnewprec(GEN nf, long prec)
{
  pari_sp av = avma;
  return gerepilecopy(av, nfnewprec_i(nf, prec));
}

GEN
nfnewprec(GEN nf, long prec)
{
  long l = lg(nf);
  GEN z, res = NULL;

  if (typ(nf) != t_VEC) pari_err(talker,"incorrect nf in nfnewprec");
  if (l == 3) {
    res = cgetg(3, t_VEC);
    gel(res,2) = gcopy(gel(nf,2));
    nf = gel(nf,1); l = lg(nf);
  }
  switch(l)
  {
    case 11: z = bnfnewprec(nf,prec); break;
    case  7: z = bnrnewprec(nf,prec); break;
    default: z = _nfnewprec(checknf(nf),prec); break;
  }
  if (res) gel(res,1) = z; else res = z;
  return res;
}

/********************************************************************/
/**                                                                **/
/**                           POLRED                               **/
/**                                                                **/
/********************************************************************/
/* remove duplicate polynomials in y, updating a (same indexes), in place */
static void
remove_duplicates(GEN y, GEN a)
{
  long k, i, l = lg(y);
  pari_sp av = avma;
  GEN z;

  if (l < 2) return;
  z = new_chunk(3);
  gel(z,1) = y;
  gel(z,2) = a; (void)sort_factor(z, gcmp);
  for  (k=1, i=2; i<l; i++)
    if (!gequal(gel(y,i), gel(y,k)))
    {
      k++;
      a[k] = a[i];
      y[k] = y[i];
    }
  l = k+1; setlg(a,l); setlg(y,l);
  avma = av;
}

/* if CHECK != NULL, return the first polynomial pol found such that
 * CHECK->f(CHECK->data, pol) != 0; return NULL if there are no such pol */
static GEN
_polred(GEN x, GEN a, GEN *pta, FP_chk_fun *CHECK)
{
  long i, v = varn(x), l = lg(a);
  GEN ch, d, y = cgetg(l,t_VEC);

  for (i=1; i<l; i++)
  {
    if (DEBUGLEVEL>2) { fprintferr("i = %ld\n",i); flusherr(); }
    ch = ZX_caract(x, gel(a,i), v);
    if (CHECK)
    {
      ch = CHECK->f(CHECK->data, ch);
      if (!ch) continue;
      return ch;
    }
    d = modulargcd(derivpol(ch), ch);
    if (degpol(d)) ch = gdivexact(ch,d);

    if (canon_pol(ch) < 0 && pta) gel(a,i) = gneg_i(gel(a,i));
    if (DEBUGLEVEL > 3) outerr(ch);
    gel(y,i) = ch;
  }
  if (CHECK) return NULL; /* no suitable polynomial found */
  remove_duplicates(y,a);
  if (pta) *pta = a;
  return y;
}

static GEN
allpolred(GEN x, long flag, GEN fa, GEN *pta, FP_chk_fun *CHECK)
{
  GEN ro = NULL;
  nfbasic_t T;
  nfbasic_init(x, flag, fa, &T);
  set_LLL_basis(&T, &ro);
  if (T.lead) pari_err(impl,"polred for non-monic polynomial");
  return _polred(T.x, T.bas, pta, CHECK);
}

/* FIXME: backward compatibility */
#define red_PARTIAL 1
#define red_ORIG    2

GEN
polred0(GEN x, long flag, GEN fa)
{
  pari_sp av = avma;
  GEN y, a;
  long fl = 0;

  if (fa && gcmp0(fa)) fa = NULL; /* compatibility */
  if (flag & red_PARTIAL) fl |= nf_PARTIALFACT;
  if (flag & red_ORIG)    fl |= nf_ORIG;
  y = allpolred(x, fl, fa, &a, NULL);
  if (fl & nf_ORIG) y = mkmat2(a,y);
  return gerepilecopy(av, y);
}

GEN
ordred(GEN x)
{
  pari_sp av = avma;
  GEN y;

  if (typ(x) != t_POL) pari_err(typeer,"ordred");
  if (!gcmp1(leading_term(x))) pari_err(impl,"ordred");
  if (!signe(x)) return gcopy(x);
  y = mkvec2(x, matid(degpol(x)));
  return gerepileupto(av, polred(y));
}

GEN
polredfirstpol(GEN x, long flag, FP_chk_fun *CHECK)
{ return allpolred(x,flag,NULL,NULL,CHECK); }
GEN
polred(GEN x) { return polred0(x, 0, NULL); }
GEN
smallpolred(GEN x) { return polred0(x, red_PARTIAL, NULL); }
GEN
factoredpolred(GEN x, GEN fa) { return polred0(x, 0, fa); }
GEN
polred2(GEN x) { return polred0(x, red_ORIG, NULL); }
GEN
smallpolred2(GEN x) { return polred0(x, red_PARTIAL|red_ORIG, NULL); }
GEN
factoredpolred2(GEN x, GEN fa) { return polred0(x, red_PARTIAL, fa); }

/********************************************************************/
/**                                                                **/
/**                           POLREDABS                            **/
/**                                                                **/
/********************************************************************/

GEN
T2_from_embed_norm(GEN x, long r1)
{
  GEN p = sum(x, 1, r1);
  GEN q = sum(x, r1+1, lg(x)-1);
  if (q != gen_0) p = gadd(p, gmul2n(q,1));
  return p;
}

GEN
T2_from_embed(GEN x, long r1)
{
  return T2_from_embed_norm(gnorm(x), r1);
}

typedef struct {
  long r1, v, prec;
  GEN ZKembed; /* embeddings of fincke-pohst-reduced Zk basis */
  GEN u; /* matrix giving fincke-pohst-reduced Zk basis */
  GEN M; /* embeddings of initial (LLL-reduced) Zk basis */
  GEN bound; /* T2 norm of the polynomial defining nf */
} CG_data;

/* characteristic pol of x (given by embeddings) */
static GEN
get_pol(CG_data *d, GEN x)
{
  long e;
  GEN g = grndtoi(roots_to_pol_r1r2(x, d->r1, d->v), &e);
  if (e > -5) pari_err(precer, "get_pol");
  return g;
}

/* characteristic pol of x (given as vector on (w_i)) */
static GEN
get_polchar(CG_data *d, GEN x)
{
  return get_pol(d, gmul(d->ZKembed,x));
}

/* return a defining polynomial for Q(w_i) */
static GEN
get_polmin_w(CG_data *d, long k)
{
  GEN g = get_pol(d, gel(d->ZKembed,k));
  GEN h = modulargcd(g, derivpol(g));
  if (degpol(h)) g = gdivexact(g,h);
  return g;
}

/* does x generate the correct field ? */
static GEN
chk_gen(void *data, GEN x)
{
  pari_sp av = avma, av1;
  GEN h, g = get_polchar((CG_data*)data,x);
  av1 = avma;
  h = modulargcd(g, derivpol(g));
  if (degpol(h)) { avma = av; return NULL; }
  if (DEBUGLEVEL>3) fprintferr("  generator: %Z\n",g);
  avma = av1; return gerepileupto(av, g);
}

/* set V[k] := matrix of multiplication by nk.zk[k] */
static GEN
set_mulid(GEN V, GEN M, GEN Mi, long r1, long r2, long N, long k)
{
  GEN v, Mk = cgetg(N+1, t_MAT);
  long i, e;
  for (i = 1; i < k; i++) gel(Mk,i) = gmael(V, i, k);
  for (     ; i <=N; i++) 
  {
    v = vecmul(gel(M,k), gel(M,i));
    v = gmul(Mi, split_realimag(v, r1, r2));
    gel(Mk,i) = grndtoi(v, &e);
    if (e > -5) return NULL;
  }
  gel(V,k) = Mk; return Mk;
}

static long
chk_gen_prec(long N, long bit)
{
  return nbits2prec(10 + (long)log2((double)N) + bit);
}

/* U = base change matrix, R = Cholesky form of the quadratic form [matrix
 * Q from algo 2.7.6] */
static GEN
chk_gen_init(FP_chk_fun *chk, GEN R, GEN U)
{
  CG_data *d = (CG_data*)chk->data;
  GEN P, V, S, inv, bound, M;
  long N = lg(U)-1, r1 = d->r1, r2 = (N-r1)>>1;
  long i, j, prec, firstprim = 0, skipfirst = 0;
  pari_sp av;

  d->u = U;
  d->ZKembed = gmul(d->M, U);

  av = avma; bound = d->bound;
  S = cgetg(N+1, t_VECSMALL);
  for (i = 1; i <= N; i++)
  {
    P = get_polmin_w(d, i);
    S[i] = degpol(P);
    if (S[i] == N)
    { /* primitive element */
      GEN B = T2_from_embed(gel(d->ZKembed,i), r1);
      if (gcmp(B,bound) < 0) bound = B ;
      if (!firstprim) firstprim = i; /* index of first primitive element */
      if (DEBUGLEVEL>2) fprintferr("chk_gen_init: generator %Z\n",P);
    }
    else
    {
      if (DEBUGLEVEL>2) fprintferr("chk_gen_init: subfield %Z\n",P);
      if (firstprim)
      { /* cycle basis vectors so that primitive elements come last */
        GEN u = d->u, e = d->ZKembed;
        GEN te = gel(e,i), tu = gel(u,i), tR = gel(R,i);
        long tS = S[i];
        for (j = i; j > firstprim; j--)
        {
          u[j] = u[j-1];
          e[j] = e[j-1];
          R[j] = R[j-1];
          S[j] = S[j-1];
        }
        gel(u,firstprim) = tu;
        gel(e,firstprim) = te;
        gel(R,firstprim) = tR;
        S[firstprim] = tS; firstprim++;
      }
    }
  }
  if (!firstprim)
  { /* try (a little) to find primitive elements to improve bound */
    GEN x = cgetg(N+1, t_COL), e, B;
    if (DEBUGLEVEL>1)
      fprintferr("chk_gen_init: difficult field, trying random elements\n");
    for (i = 0; i < 10; i++)
    {
      for (j = 1; j <= N; j++) gel(x,j) = stoi( (pari_rand31() % 7) - 3 );
      e = gmul(d->ZKembed, x);
      P = get_pol(d, e);
      if (!ZX_is_squarefree(P)) continue;
      if (DEBUGLEVEL>2) fprintferr("chk_gen_init: generator %Z\n",P);
      B = T2_from_embed(e, r1);
      if (gcmp(B,bound) < 0) bound = B ;
    }
  }

  if (firstprim > 1)
  {
    inv = ginv( split_realimag(d->ZKembed, r1, r2) ); /*TODO: use QR?*/
    V = gel(inv,1);
    for (i = 2; i <= r1+r2; i++) V = gadd(V, gel(inv,i)); 
    /* V corresponds to 1_Z */
    V = grndtoi(V, &j);
    if (j > -5) err(bugparier,"precision too low in chk_gen_init");
    M = mkmat(V); /* 1 */

    V = cgetg(N+1, t_VEC);
    for (i = 1; i <= N; i++)
    { /* M = Q-basis of subfield generated by nf.zk[1..i-1] */
      GEN Mx, M2;
      long j, k, h, rkM, dP = S[i];

      if (dP == N) break; /* primitive */
      Mx = set_mulid(V, d->ZKembed, inv, r1, r2, N, i);
      if (!Mx) break; /* prec. problem. Stop */
      if (dP == 1) continue;
      rkM = lg(M)-1;
      M2 = cgetg(N+1, t_MAT); /* we will add to M the elts of M2 */
      gel(M2,1) = col_ei(N, i); /* nf.zk[i] */
      k = 2;
      for (h = 1; h < dP; h++)
      { 
        long r; /* add to M2 the elts of M * nf.zk[i]  */
        for (j = 1; j <= rkM; j++) gel(M2,k++) = gmul(Mx, gel(M,j));
        setlg(M2, k); k = 1;
        M = image(shallowconcat(M, M2));
        r = lg(M) - 1;
        if (r == rkM) break;
        if (r > rkM)
        {
          rkM = r;
          if (rkM == N) break;
        }
      }
      if (rkM == N) break;

      /* Q(w[1],...,w[i-1]) is a strict subfield of nf */
      skipfirst++;
    }
  }
  /* x_1,...,x_skipfirst generate a strict subfield [unless N=skipfirst=1] */
  chk->skipfirst = skipfirst;
  if (DEBUGLEVEL>2) fprintferr("chk_gen_init: skipfirst = %ld\n",skipfirst);

  /* should be DEF + gexpo( max_k C^n_k (bound/k)^(k/2) ) */
  bound = gerepileuptoleaf(av, bound);
  prec = chk_gen_prec(N, (gexpo(bound)*N)/2);
  if (DEBUGLEVEL)
    fprintferr("chk_gen_init: new prec = %ld (initially %ld)\n", prec, d->prec);
  if (prec > d->prec) pari_err(bugparier, "polredabs (precision problem)");
  if (prec < d->prec) d->ZKembed = gprec_w(d->ZKembed, prec);
  return bound;
}

/* store phi(beta mod z). */
static GEN
storeeval(GEN a, GEN x, GEN z, GEN lead)
{
  GEN beta = modreverse_i(a, x);
  if (lead) beta = gdiv(beta, lead);
  return mkvec2(z, mkpolmod(beta,z));
}

static void
findmindisc(GEN *py, GEN *pa)
{
  GEN v, dmin, z, b, discs, y = *py, a = *pa;
  long i,k, l = lg(y);

  if (l == 2) { *py = gel(y,1); *pa = gel(a,1); return; }

  discs = cgetg(l,t_VEC);
  for (i=1; i<l; i++) gel(discs,i) = absi(ZX_disc(gel(y,i)));
  v = sindexsort(discs);
  k = v[1];
  dmin = gel(discs,k);
  z    = gel(y,k);
  b = gel(a,k);
  for (i=2; i<l; i++)
  {
    k = v[i];
    if (!equalii(gel(discs,k), dmin)) break;
    if (gpolcomp(gel(y,k),z) < 0) { z = gel(y,k); b = gel(a,k); }
  }
  *py = z; *pa = b;
}

static GEN
storepol(GEN x, GEN z, GEN a, GEN lead, long flag)
{
  GEN y;
  if (flag & nf_RAW)
    y = mkvec2(z, a);
  else if (flag & nf_ORIG)
    y = storeeval(a, x, z, lead);
  else
    y = z;
  return y;
}

static GEN
storeallpol(GEN x, GEN z, GEN a, GEN lead, long flag)
{
  GEN y;

  if (flag & nf_RAW)
  {
    long i, c = lg(z);
    y = cgetg(c,t_VEC);
    for (i=1; i<c; i++) gel(y,i) = mkvec2(gel(z,i), gel(a,i));
  }
  else if (flag & nf_ORIG)
  {
    long i, c = lg(z);
    y = cgetg(c,t_VEC);
    for (i=1; i<c; i++) gel(y,i) = storeeval(gel(a,i), x, gel(z,i), lead);
  }
  else
    y = z;
  return y;
}

static GEN
_polredabs(nfbasic_t *T, GEN *u)
{
  long prec, e, n = degpol(T->x);
  GEN v, ro = NULL;
  FP_chk_fun chk = { &chk_gen, &chk_gen_init, NULL, NULL, 0 };
  nffp_t F;
  CG_data d; chk.data = (void*)&d;

  set_LLL_basis(T, &ro);

  /* || polchar ||_oo < 2^e */
  e = n * (long)(cauchy_bound(T->x) / LOG2 + log2((double)n)) + 1;
  prec = chk_gen_prec(n, e);
  get_nf_fp_compo(T, &F, ro, prec);

  d.v = varn(T->x);
  d.r1= T->r1;
  d.bound = T2_from_embed(F.ro, T->r1);
  for (;;)
  {
    GEN R = R_from_QR(F.G, prec);
    if (R)
    {
      d.prec = prec;
      d.M    = F.M;
      v = fincke_pohst(mkvec(R),NULL,-1, 0, &chk);
      if (v) break;
    }
    prec = (prec<<1)-2;
    get_nf_fp_compo(T, &F, NULL, prec);
    if (DEBUGLEVEL) pari_warn(warnprec,"polredabs0",prec);
  }
  *u = d.u; return v;
}

GEN
polredabs0(GEN x, long flag)
{
  pari_sp av = avma;
  long i, l, vx;
  GEN y, a, u;
  nfbasic_t T;

  nfbasic_init(x, flag & nf_PARTIALFACT, NULL, &T);
  x = T.x; vx = varn(x);

  if (degpol(x) == 1)
  {
    u = NULL;
    y = mkvec(pol_x[vx]);
    a = mkvec(gsub(gel(y,1), gel(x,2)));
  }
  else
  {
    GEN v = _polredabs(&T, &u);
    y = gel(v,1);
    a = gel(v,2);
  }
  l = lg(a);
  for (i=1; i<l; i++)
    if (canon_pol(gel(y,i)) < 0) gel(a,i) = gneg_i(gel(a,i));
  remove_duplicates(y,a);
  l = lg(a);
  if (l == 1)
  {
    y = mkvec(x);
    a = mkvec(pol_x[vx]);
  }
  if (DEBUGLEVEL) fprintferr("Found %ld minimal polynomials.\n",l-1);
  if (flag & nf_ALL)
  {
    if (u) for (i=1; i<l; i++) gel(a,i) = gmul(T.bas, gmul(u, gel(a,i)));
    y = storeallpol(x, y, a, T.lead, flag);
    if (flag & nf_ADDZK) pari_err(impl,"nf_ADDZK flag when nf_ALL set (polredabs)");
  }
  else
  {
    GEN z;
    findmindisc(&y, &a); z = y;
    if (u && l > 1) a = gmul(T.bas, gmul(u, a));
    y = storepol(x, y, a, T.lead, flag);
    if (flag & nf_ADDZK)
    {
      GEN t, y0 = y, B = RgXV_to_RgM(T.bas, lg(T.bas)-1);
      t = (flag & nf_ORIG)? lift_intern(gel(y,2)): modreverse_i(a, x);
      t = gmul(RgX_powers(t, z, degpol(z)-1), B);
      y = mkvec2(y0, t);
    }
  }
  return gerepilecopy(av, y);
}

GEN
polredabsall(GEN x, long flun) { return polredabs0(x, flun | nf_ALL); }
GEN
polredabs(GEN x) { return polredabs0(x,0); }
GEN
polredabs2(GEN x) { return polredabs0(x,nf_ORIG); }

static long
nf_pm1(GEN y)
{
  long i,l;

  if (!is_pm1(y[1])) return 0;
  l = lg(y);
  for (i=2; i<l; i++)
    if (signe(y[i])) return 0;
  return signe(y[1]);
}

static GEN
is_primitive_root(GEN nf, GEN fa, GEN x, long w)
{
  GEN y, exp = utoipos(2), pp = gel(fa,1);
  long i,p, l = lg(pp);

  for (i=1; i<l; i++)
  {
    p = itos(gel(pp,i));
    exp[2] = w / p; y = element_pow(nf,x,exp);
    if (nf_pm1(y) > 0) /* y = 1 */
    {
      if (p!=2 || !gcmp1(gcoeff(fa,i,2))) return NULL;
      x = gneg_i(x);
    }
  }
  return x;
}

/* only roots of 1 are +/- 1 */
static GEN
trivroots(GEN nf) {
  GEN y = cgetg(3, t_VEC);
  gel(y,1) = gen_2;
  gel(y,2) = gscalcol_i(gen_m1, degpol(nf[1]));
  return y;
}

GEN
rootsof1(GEN nf)
{
  pari_sp av = avma;
  long N, k, i, ws, prec;
  GEN z, y, d, list, w;

  nf = checknf(nf);
  if ( nf_get_r1(nf) ) return trivroots(nf);

  N = degpol(nf[1]); prec = nfgetprec(nf);
  for (;;)
  {
    GEN R = R_from_QR(gmael(nf,5,2), prec);
    if (R)
    {
      y = fincke_pohst(mkvec(R), utoipos(N), 1000, 0, NULL);
      if (y) break;
    }
    prec = (prec<<1)-2;
    if (DEBUGLEVEL) pari_warn(warnprec,"rootsof1",prec);
    nf = nfnewprec(nf,prec);
  }
  if (itos(ground(gel(y,2))) != N) pari_err(bugparier,"rootsof1 (bug1)");
  w = gel(y,1); ws = itos(w);
  if (ws == 2) { avma = av; return trivroots(nf); }

  d = Z_factor(w); list = gel(y,3); k = lg(list);
  for (i=1; i<k; i++)
  {
    z = is_primitive_root(nf, d, gel(list,i), ws);
    if (z) return gerepilecopy(av, mkvec2(w, z));
  }
  pari_err(bugparier,"rootsof1");
  return NULL; /* not reached */
}

/*******************************************************************/
/*                                                                 */
/*                     DEDEKIND ZETA FUNCTION                      */
/*                                                                 */
/*******************************************************************/
static GEN
dirzetak0(GEN nf, long N0)
{
  GEN vect, c, c2, pol = gel(nf,1), index = gel(nf,4);
  pari_sp av = avma;
  long i, k, limk, lx;
  ulong q, p;
  byteptr d = diffptr;
  long court[] = {evaltyp(t_INT)|_evallg(3), evalsigne(1)|evallgefint(3),0};

  (void)evallg(N0+1);
  c  = cgetalloc(t_VECSMALL, N0+1);
  c2 = cgetalloc(t_VECSMALL, N0+1);
  c2[1] = c[1] = 1; for (i=2; i<=N0; i++) c[i] = 0;

  maxprime_check((ulong)N0);
  court[2] = 0;
  while (court[2] <= N0)
  {
    NEXT_PRIME_VIADIFF(court[2], d);
    if (umodiu(index, court[2])) /* court does not divide index */
      { vect = (GEN) FpX_degfact(pol,court)[1]; lx = lg(vect); }
    else
    {
      GEN p1 = primedec(nf,court); lx = lg(p1); vect = cgetg(lx,t_VECSMALL);
      for (i=1; i<lx; i++) vect[i] = itou( gmael(p1,i,4) );
    }
    for (i=1; i<lx; i++)
    {
      GEN N = powiu(court, vect[i]); /* N = court^f */
      if (cmpiu(N, N0) > 0) break;

      q = p = (ulong)N[2]; limk = N0/q;
      for (k=2; k <= N0; k++) c2[k] = c[k];
      while (q <= (ulong)N0)
      {
        LOCAL_HIREMAINDER;
        for (k=1; k<=limk; k++) c2[k*q] += c[k];
        q = mulll(q, p); if (hiremainder) break;
        limk /= p;
      }
      N = c; c = c2; c2 = N;
    }
    avma = av;
  }
  free(c2); return c;
}

GEN
dirzetak(GEN nf, GEN b)
{
  GEN z, c;
  long n;

  if (typ(b) != t_INT) pari_err(talker,"not an integer type in dirzetak");
  if (signe(b) <= 0) return cgetg(1,t_VEC);
  nf = checknf(nf);
  n = itos_or_0(b); if (!n) pari_err(talker,"too many terms in dirzetak");
  c = dirzetak0(nf, n);
  z = vecsmall_to_vec(c); free(c); return z;
}

GEN
zeta_get_limx(long r1, long r2, long bit)
{
  pari_sp av = avma;
  GEN p1, p2, c0, c1, A0;
  long r = r1 + r2, N = r + r2;

  /* c1 = N 2^(-2r2 / N) */
  c1 = mulrs(powrfrac(real2n(1, DEFAULTPREC), -2*r2, N), N);

  p1 = gpowgs(Pi2n(1, DEFAULTPREC), r - 1);
  p2 = gmul2n(mpmul(powuu(N,r), p1), -r2);
  c0 = sqrtr( mpdiv(p2, gpowgs(c1, r+1)) );

  A0 = logr_abs( gmul2n(c0, bit) ); p2 = gdiv(A0, c1);
  p1 = divrr(mulsr(N*(r+1), logr_abs(p2)), addsr(2*(r+1), gmul2n(A0,2)));
  return gerepileuptoleaf(av, divrr(addrs(p1, 1), powrshalf(p2, N)));
}
/* N_0 = floor( C_K / limx ) */
long
zeta_get_N0(GEN C,  GEN limx)
{
  long e;
  pari_sp av = avma;
  GEN z = gcvtoi(gdiv(C, limx), &e); /* avoid truncation error */
  if (e >= 0 || is_bigint(z))
    pari_err(talker, "need %Z coefficients in initzeta: computation impossible", z);
  if (DEBUGLEVEL>1) fprintferr("\ninitzeta: N0 = %Z\n", z);
  avma = av; return itos(z);
}

static long
get_i0(long r1, long r2, GEN B, GEN limx)
{
  long imin = 1, imax   = 1400;
  while(imax - imin >= 4)
  {
    long i = (imax + imin) >> 1;
    GEN t = gpowgs(limx, i);
    t = gmul(t, gpowgs(mpfactr(i/2, DEFAULTPREC), r1));
    t = gmul(t, gpowgs(mpfactr(i  , DEFAULTPREC), r2));
    if (gcmp(t, B) >= 0) imax = i; else imin = i;
  }
  return imax & ~1; /* make it even */
}

/* assume limx = zeta_get_limx(r1, r2, bit) */
long
zeta_get_i0(long r1, long r2, long bit, GEN limx)
{
  pari_sp av = avma;
  GEN B = gmul(sqrtr( gdiv(gpowgs(mppi(DEFAULTPREC), r2-3), limx) ),
               gmul2n(powuu(5, r1), bit + r2));
  long i0 = get_i0(r1, r2, B, limx);
  if (DEBUGLEVEL>1) { fprintferr("i0 = %ld\n",i0); flusherr(); }
  avma = av; return i0;
}

GEN
initzeta(GEN pol, long prec)
{
  GEN nfz, nf, gr1, gr2, gru, p1, p2, cst, coef, bnf = checkbnf_i(pol);
  GEN limx, resi,zet,C,coeflog,racpi,aij,tabj,colzero, *tabcstn, *tabcstni;
  GEN c_even, ck_even, c_odd, ck_odd, serie_even, serie_odd, serie_exp, Pi;
  long N0, i0, r1, r2, r, R, N, i, j, k, n, bit = bit_accuracy(prec) + 6;
  pari_sp av, av2;
  long court[] = {evaltyp(t_INT)|_evallg(3), evalsigne(1)|evallgefint(3),0};
  stackzone *zone, *zone0, *zone1;

  /*************** residue & constants ***************/
  nfz = cgetg(10,t_VEC);
  if (!bnf || nfgetprec(bnf) < prec ) bnf = bnfinit0(pol, 2, NULL, prec);
  prec = (prec<<1) - 2;
  Pi = mppi(prec); racpi = sqrtr(Pi);

  /* class number & regulator */
  nf = gel(bnf,7); N = degpol(nf[1]);
  nf_get_sign(nf, &r1, &r2);
  gr1 = gmael(nf,2,1); gr2 = gmael(nf,2,2);
  r = r1 + r2; R = r+2;
  av = avma; p1 = gel(bnf,8); p2 = gmul(gmul2n(gmael(p1,1,1),r1), gel(p1,2));
  resi = gerepileupto(av, gdiv(p2, gmael(p1,4,1)));

  av = avma;
  p1 = sqrtr_abs(itor(gel(nf,3), prec));
  p2 = gmul2n(gpowgs(racpi,N), r2);
  cst = gerepileuptoleaf(av, divrr(p1,p2));

  /* N0, i0 */
  limx = zeta_get_limx(r1, r2, bit);
  N0 = zeta_get_N0(cst, limx);
  i0 = zeta_get_i0(r1, r2, bit + 4, limx);

  /* Array of i/cst (i=1..N0) */
  av = avma; i = prec*N0;
  zone  = switch_stack(NULL,i + 2*(N0+1) + 6*prec);
  zone1 = switch_stack(NULL,2*i);
  zone0 = switch_stack(NULL,2*i);
  (void)switch_stack(zone,1);
  tabcstn  = (GEN*) cgetg(N0+1,t_VEC);
  tabcstni = (GEN*) cgetg(N0+1,t_VEC);
  p1 = ginv(cst);
  for (i=1; i<=N0; i++) tabcstni[i] = tabcstn[i] = mulsr(i,p1);
  (void)switch_stack(zone,0);

  /********** compute a(i,j) **********/
  zet = cgetg(R,t_VEC); gel(zet,1) = mpeuler(prec);
  for (i=2; i<R; i++) gel(zet,i) = szeta(i, prec);

  aij = cgetg(i0+1,t_VEC);
  for (i=1; i<=i0; i++) gel(aij,i) = cgetg(R,t_VEC);

  c_even = real2n(r1, prec);
  c_odd = gmul(c_even, gpowgs(racpi,r1));
  if (r&1) c_odd = gneg_i(c_odd);
  ck_even = cgetg(R,t_VEC); ck_odd = cgetg(r2+2,t_VEC);
  for (k=1; k<R; k++)
  {
    GEN t = gmul(gel(zet,k), gadd(gr2, gmul2n(gr1,-k)));
    if (k&1) t = gneg(t);
    gel(ck_even,k) = t;
  }
  gru = utoipos(r);
  for (k = 1; k <= r2+1; k++)
  {
    GEN t = gmul(gel(zet,k), gsub(gru, gmul2n(gr1,-k)));
    if (k&1) t = gneg(t);
    gel(ck_odd,k) = gadd(gru, t);
  }
  gel(ck_odd,1) = gsub(gel(ck_odd,1), mulsr(r1, mplog2(prec)));
  serie_even = cgetg(r+3,t_SER);
  serie_odd = cgetg(r2+3,t_SER);
  serie_even[1] = serie_odd[1] = evalsigne(1)+evalvalp(1);
  i = 0;
  while (i < i0/2)
  {
    for (k=1; k<R; k++) gel(serie_even,k+1) = gdivgs(gel(ck_even,k),k);
    serie_exp = gmul(c_even, gexp(serie_even,0));
    p1 = gel(aij,2*i+1);
    for (j=1; j<R; j++) p1[j] = serie_exp[r+3-j];

    for (k=1; k<=r2+1; k++) gel(serie_odd,k+1) = gdivgs(gel(ck_odd,k),k);
    serie_exp = gmul(c_odd, gexp(serie_odd,0));
    p1 = gel(aij,2*i+2);
    for (j=1; j<=r2+1; j++) p1[j] = serie_exp[r2+3-j];
    for (   ; j<R; j++) gel(p1,j) = gen_0;
    i++;

    c_even = gdiv(c_even,gmul(powuu(i,r),powuu(2*i-1,r2)));
    c_odd  = gdiv(c_odd, gmul(powuu(i,r2),powuu(2*i+1,r)));
    c_even = gmul2n(c_even,-r2);
    c_odd  = gmul2n(c_odd,r1-r2);
    if (r1 & 1) { c_even = gneg_i(c_even); c_odd = gneg_i(c_odd); }
    p1 = gr2; p2 = gru;
    for (k=1; k<R; k++)
    {
      p1 = gdivgs(p1,2*i-1); p2 = gdivgs(p2,2*i);
      gel(ck_even,k) = gadd(gel(ck_even,k), gadd(p1,p2));
    }
    p1 = gru; p2 = gr2;
    for (k=1; k<=r2+1; k++)
    {
      p1 = gdivgs(p1,2*i+1); p2 = gdivgs(p2,2*i);
      gel(ck_odd,k) = gadd(gel(ck_odd,k), gadd(p1,p2));
    }
  }
  aij = gerepilecopy(av, aij);
  if (DEBUGLEVEL>1) msgtimer("a(i,j)");
  p1=cgetg(5,t_VEC);
  gel(p1,1) = stoi(r1);
  gel(p1,2) = stoi(r2);
  gel(p1,3) = stoi(i0);
  gel(p1,4) = bnf;
  gel(nfz,1) = p1;
  gel(nfz,2) = resi;
  gel(nfz,5) = cst;
  gel(nfz,6) = glog(cst,prec);
  gel(nfz,7) = aij;

  /************* Calcul du nombre d'ideaux de norme donnee *************/
  coef = dirzetak0(nf,N0); tabj = cgetg(N0+1,t_MAT);
  if (DEBUGLEVEL>1) msgtimer("coef");
  colzero = zerocol(r+1);
  for (i=1; i<=N0; i++)
    if (coef[i])
    {
      GEN t = cgetg(r+2,t_COL);
      p1 = logr_abs(gel(tabcstn,i)); setsigne(p1, -signe(p1));
      gel(t,1) = utoi(coef[i]);
      gel(t,2) = mulur(coef[i], p1);
      for (j=2; j<=r; j++)
      {
        pari_sp av2 = avma;
        gel(t,j+1) = gerepileuptoleaf(av2, divrs(mulrr(gel(t,j), p1), j));
      }
      gel(tabj,i) = t; /* tabj[n,j] = coef(n)*ln(c/n)^(j-1)/(j-1)! */
    }
    else
      gel(tabj,i) = colzero;
  if (DEBUGLEVEL>1) msgtimer("a(n)");

  coeflog=cgetg(N0+1,t_VEC); gel(coeflog,1) = gen_0;
  for (i=2; i<=N0; i++)
    if (coef[i])
    {
      court[2] = i; p1 = glog(court,prec);
      setsigne(p1, -1); gel(coeflog,i) = p1;
    }
    else gel(coeflog,i) = gen_0;
  if (DEBUGLEVEL>1) msgtimer("log(n)");

  gel(nfz,3) = tabj;
  gel(nfz,8) = vecsmall_copy(coef);
  gel(nfz,9) = coeflog;

  /******************** Calcul des coefficients Cik ********************/
  C = cgetg(r+1,t_MAT);
  for (k=1; k<=r; k++) gel(C,k) = cgetg(i0+1,t_COL);
  av2 = avma;
  for (i=1; i<=i0; i++)
  {
    GEN aiji = gel(aij,i);
    stackzone *z;
    for (k=1; k<=r; k++)
    {
      p1 = NULL;
      for (n=1; n<=N0; n++)
        if (coef[n])
        {
          GEN tabjn = gel(tabj,n), p2 = mpmul(gel(aiji,1+k), gel(tabjn,1));
          for (j=2; j<=r-k+1; j++)
            p2 = mpadd(p2, mpmul(gel(aiji,j+k), gel(tabjn,j)));
          if (i > 1) p2 = mpmul(p2, tabcstni[n]);
          p1 = p1? mpadd(p1,p2): p2;
        }
      gcoeff(C,i,k) = gerepileuptoleaf(av2,p1);
      av2 = avma;
    }
    if (i > 1 && i < i0) {
      /* use a parallel stack */
      z = i&1? zone1: zone0;
      (void)switch_stack(z, 1);
      for (n=1; n<=N0; n++)
        if (coef[n]) tabcstni[n] = mpmul(tabcstni[n],tabcstn[n]);
      /* come back */
      (void)switch_stack(z, 0);
    }
  }
  gel(nfz,4) = C;
  if (DEBUGLEVEL>1) msgtimer("Cik");
  free((void*)zone); free((void*)zone1); free((void*)zone0);
  free((void*)coef); return nfz;
}

GEN
gzetakall(GEN nfz, GEN s, long flag, long prec2)
{
  GEN resi,C,cst,cstlog,coeflog,cs,coef;
  GEN lambd,gammas,gammaunmoins,gammas2,gammaunmoins2,var1,var2;
  GEN smoinun,unmoins,gexpro,gar,val,valm,valk,valkm;
  long ts = typ(s), r1,r2,ru,i0,i,j,k,N0,sl,prec,bigprec;
  pari_sp av = avma;

  if (typ(nfz)!=t_VEC || lg(nfz)!=10 || typ(nfz[1]) != t_VEC)
    pari_err(talker,"not a zeta number field in zetakall");
  if (! is_intreal_t(ts) && ts != t_COMPLEX && ts != t_FRAC)
    pari_err(typeer,"gzetakall");
  resi=gel(nfz,2); C=gel(nfz,4); cst=gel(nfz,5);
  cstlog=gel(nfz,6); coef=gel(nfz,8); coeflog=gel(nfz,9);
  r1 = itos(gmael(nfz,1,1));
  r2 = itos(gmael(nfz,1,2)); ru = r1+r2;
  i0 = itos(gmael(nfz,1,3)); N0 = lg(coef)-1;
  bigprec = precision(cst); prec = prec2+1;

  if (ts==t_COMPLEX && gcmp0(imag_i(s))) { s=real_i(s); ts = typ(s); }
  if (ts==t_REAL && !signe(gfrac(s))) { s=mptrunc(s); ts = t_INT; }
  if (ts==t_INT)
  {
    sl = itos(s);
    if (sl==1) pari_err(talker,"s = 1 is a pole (gzetakall)");
    if (sl==0)
    {
      avma = av;
      if (flag) pari_err(talker,"s = 0 is a pole (gzetakall)");
      if (ru == 1) return gneg(r1? ghalf: resi);
      return gen_0;
    }
    if (sl<0 && (r2 || !odd(sl)))
    {
      if (!flag) return gen_0;
      s = subsi(1,s); sl = 1-sl;
    }
    unmoins=subsi(1,s);
    lambd = gdiv(resi, mulis(s,sl-1));
    gammas2=ggamma(gmul2n(s,-1),prec);
    gar=gpowgs(gammas2,r1);
    cs=gexp(gmul(cstlog,s),prec); 	
    val=s; valm=unmoins;
    if (sl < 0) /* r2 = 0 && odd(sl) */
    {
      gammaunmoins2 = ggamma(gmul2n(unmoins,-1),prec);
      var1 = var2 = gen_1;
      for (i=2; i<=N0; i++)
	if (coef[i])
	{
          gexpro=gexp(gmul(gel(coeflog,i),s),bigprec);
	  var1 = gadd(var1,gmulsg(coef[i],gexpro));
	  var2 = gadd(var2,gdivsg(coef[i],gmulsg(i,gexpro)));
	}
      lambd=gadd(lambd,gmul(gmul(var1,cs),gar));
      lambd=gadd(lambd,gmul(gmul(var2,gdiv(cst,cs)),
			    gpowgs(gammaunmoins2,r1)));
      var1 = gen_0;
      for (i=1; i<=i0; i+=2)
      {
	valk  = val;
        valkm = valm;
	for (k=1; k<=ru; k++)	
	{
          GEN c = gcoeff(C,i,k);
	  var1 = gsub(var1,gdiv(c,valk )); valk  = mulii(val,valk);
	  var1 = gsub(var1,gdiv(c,valkm)); valkm = mulii(valm,valkm);
	}
	val  = addis(val, 2);
        valm = addis(valm,2);
      }
    }
    else
    {
      GEN tabj=gel(nfz,3), aij=gel(nfz,7);

      gar = gmul(gar,gpowgs(ggamma(s,prec),r2));
      var1=var2=gen_0;
      for (i=1; i<=N0; i++)
	if (coef[i])
	{
	  gexpro=gexp(gmul(gel(coeflog,i),s),bigprec);
	  var1 = gadd(var1,gmulsg(coef[i],gexpro));
          if (sl <= i0)
          {
            GEN t=gen_0;
            for (j=1; j<=ru+1; j++)
              t = gadd(t, gmul(gmael(aij,sl,j), gmael(tabj,i,j)));
            var2=gadd(var2,gdiv(t,gmulsg(i,gexpro)));
          }
	}
      lambd=gadd(lambd,gmul(gmul(var1,cs),gar));
      lambd=gadd(lambd,gmul(var2,gdiv(cst,cs)));
      var1 = gen_0;
      for (i=1; i<=i0; i++)
      {
	valk  = val;
        valkm = valm;
        if (i == sl)
          for (k=1; k<=ru; k++)
          {	
            GEN c = gcoeff(C,i,k);
            var1 = gsub(var1,gdiv(c,valk)); valk = mulii(val,valk);
          }
        else
	for (k=1; k<=ru; k++)
	{	
            GEN c = gcoeff(C,i,k);
            var1 = gsub(var1,gdiv(c,valk )); valk  = mulii(val,valk);
            var1 = gsub(var1,gdiv(c,valkm)); valkm = mulii(valm,valkm);
	}
	val  = addis(val,1);
        valm = addis(valm,1);
      }
    }
  }
  else
  {
    GEN Pi = mppi(bigprec);
    if (ts == t_FRAC)
      s = gmul(s, real_1(bigprec));
    else
      s = gprec_w(s, bigprec);

    smoinun = gsubgs(s,1);
    unmoins = gneg(smoinun);
    lambd = gdiv(resi,gmul(s,smoinun));
    gammas = ggamma(s,prec);
    gammas2= ggamma(gmul2n(s,-1),prec);
    gar = gmul(gpowgs(gammas,r2),gpowgs(gammas2,r1));
    cs = gexp(gmul(cstlog,s),prec);
    var1 = gmul(Pi,s);
    gammaunmoins = gdiv(Pi, gmul(gsin(var1,prec),gammas));
    gammaunmoins2= gdiv(gmul(gmul(sqrtr(Pi),gpow(gen_2,smoinun,prec)),
                             gammas2),
                        gmul(gcos(gmul2n(var1,-1),prec),gammas));
    var1 = var2 = gen_1;
    for (i=2; i<=N0; i++)
      if (coef[i])
      {
        gexpro = gexp(gmul(gel(coeflog,i),s),bigprec);
	var1 = gadd(var1,gmulsg(coef[i], gexpro));
	var2 = gadd(var2,gdivsg(coef[i], gmulsg(i,gexpro)));
      }
    lambd = gadd(lambd,gmul(gmul(var1,cs),gar));
    lambd = gadd(lambd,gmul(gmul(gmul(var2,gdiv(cst,cs)),
	 		         gpowgs(gammaunmoins,r2)),
                            gpowgs(gammaunmoins2,r1)));
    val  = s;
    valm = unmoins;
    var1 = gen_0;
    for (i=1; i<=i0; i++)
    {
      valk  = val;
      valkm = valm;
      for (k=1; k<=ru; k++)
      {
        GEN c = gcoeff(C,i,k);
	var1 = gsub(var1,gdiv(c,valk )); valk  = gmul(val, valk);
	var1 = gsub(var1,gdiv(c,valkm)); valkm = gmul(valm,valkm);
      }
      if (r2)
      {
	val  = gaddgs(val, 1);
        valm = gaddgs(valm,1);
      }
      else
      {
	val  = gaddgs(val, 2);
        valm = gaddgs(valm,2); i++;
      }
    }
  }
  lambd = gadd(lambd, var1);
  if (!flag) lambd = gdiv(lambd,gmul(gar,cs)); /* zetak */
  if (gprecision(lambd) > prec2) lambd = gprec_w(lambd, prec2);
  return gerepileupto(av, lambd);
}

GEN
gzetak(GEN nfz, GEN s, long prec) { return gzetakall(nfz,s,0,prec); }

GEN
glambdak(GEN nfz, GEN s, long prec) { return gzetakall(nfz,s,1,prec); }
