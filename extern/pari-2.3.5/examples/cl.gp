rgcd(a,b)=
{ local(r);

  a = abs(a);
  b = abs(b);
  while (b > 0.01, r=a%b; a=b; b=r);
  a
}

f(a,b)=
{ local(j,n,l,mv,u,vreg,cp,p);

  n = abs( norm(a + b*t) );
  mv = vectorv(li);
  forprime (p=2, plim,
    if ((l = valuation(n,p)),
      n /= p^l;
      j = ind[p]; cp = v[j][2];
      while((a+b*cp)%p,
	j++; cp = v[j][2]
      );
      mv[j]=l
    )
  );
  if (n!=1, return);

  /* found a relation */
  vreg = vectorv(lireg,j,
    u = a+b*re[j];
    if (j<=r1, abs(u), norm(u))
  );
  mreg = concat(mreg, log(vreg));
  m = concat(m,mv);
  areg = concat(areg, a+b*t);
  print1("(" res++ ": " a "," b ")")
}

global(km,clh,R,nf,areg);

clareg(pol, plim=19, lima=50, extra=5)=
{ local(e,t,r,coreg,lireg,r1,ind,fa,li,co,res,a,b,m,mh,ms,mhs,mreg,mregh);

  nf=nfinit(pol); pol=nf.pol;
  t=Mod(x,pol,1);
  r=nf.roots; r1=nf.sign[1];
  if (nf[4] > 1, /* index: power basis <==> nf[4] = 1 */
    error("sorry, the case f>1 is not implemented")
  );

  print("discriminant = " nf.disc ", signature = " nf.sign);
 
  lireg = (poldegree(pol) + r1) / 2; /* r1 + r2 */
  re=vector(lireg,j,
    if (j<=r1, real(r[j]) , r[j])
  );
  ind=vector(plim); v=[];
  forprime(p=2,plim,
    w = factormod(pol,p);
    e = w[,2];
    find=0;
    for(l=1,#e,
      fa = lift(w[l,1]);
      if (poldegree(fa) == 1,
	if (!find,
	  find=1; ind[p]=#v+1
	);
	v = concat(v, [[p,-polcoeff(fa,0),e[l]]])
      )
    )
  );
  li=#v; co=li+extra;
  res=0; print("need " co " relations");
  areg=[]~; mreg = m = [;];
  a=1; b=1; f(0,1);
  while (res<co,
    if (gcd(a,b)==1,
      f(a,b); f(-a,b)
    );
    a++;
    if (a*b>lima, b++; a=1)
  );
  print(" ");
  mh=mathnf(m); ms=matsize(mh);
  if (ms[1]!=ms[2],
    print("not enough relations for class group: matrix size = ",ms);
    return
  );

  mhs = matsnf(mh,4);
  clh = prod(i=1,#mhs, mhs[i]);
  print("class number = " clh ", class group = " mhs);
  areg=Mat(areg); km=matkerint(m); mregh=mreg*km;
  if (lireg==1,
    R = 1
  ,
    coreg = #mregh;
    if (coreg < lireg-1,
      print("not enough relations for regulator: matsize = " matsize(mregh));
      R = "(not given)";
    ,
      mreg1 = vecextract(mregh, Str(".." lireg-1), "..");
      R = 0;
      for(j=lireg-1,coreg,
        a = matdet(vecextract(mreg1, Str(j-lireg+2 ".." j)));
	R = rgcd(a,R)
      )
    )
  );
  print("regulator = " R)
}

check(lim=200) =
{ local(r1,r2,pol,z,Res,fa);

  r1=nf.sign[1];
  r2=nf.sign[2]; pol=nf.pol;
  z = 2^r1 * (2*Pi)^r2 / sqrt(abs(nf.disc)) / nfrootsof1(nf)[1];
  Res = 1.; \\ Res (Zeta_K,s=1) ~ z * h * R
  forprime (q=2,lim,
    fa = factormod(pol,q,1)[,1];
    Res *= (q-1)/q / prod(i=1, #fa, 1 - q^(-fa[i]))
  );
  z * clh * R / Res 
}

fu() = vector(#km, k, factorback(concat(areg, km[,k])))
