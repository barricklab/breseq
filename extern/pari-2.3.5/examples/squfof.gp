\\ return one non-trivial divisor of n > 1 using Shanks's SQUFOF
squfof(n) =
{
  if (isprime(n) || issquare(n, &n), return(n));

  p = factor(n,0)[1,1];
  if (p != n, return(p));

  if (n%4==1,
    D = n;      d = sqrtint(D); b = (((d-1)\2) << 1) + 1
  ,
    D = n << 2; d = sqrtint(D); b = (d\2) << 1
  );
  f = Qfb(1, b, (b^2-D)>>2);
  l = sqrtint(d);

  q = []; lq = 0; i = 0;
  while (1,
    i++;
    f = qfbred(f, 3, D, d);
    a = component(f, 1);
    if (!(i%2) && issquare(a, &as),
      j = 1;
      while (j<=lq,
        if (as == q[j], break);
        j++
      );
      if (j > lq, break)
    );

    if (abs(a) <= l,
      q = concat(q, abs(a));
      print(q); lq++;
    )
  );

  print("i = ", i); print(f);
  bb = component(f, 2);
  gs = gcd([as, bb, D]);
  if (gs > 1, return (gs));

  f = Qfb(as, -bb, as*component(f,3));
  g = qfbred(f, 3, D, d);
  b = component(g, 2);
  until (b1 == b,
    b1 = b; g = qfbred(g, 3, D, d);
    b = component(g, 2);
  );
  a = abs(component(g, 1));
  if (a % 2, a, a>>1);
}
