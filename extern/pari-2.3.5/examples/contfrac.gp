period(D) =
{ local(u,v,j,r,s);

  if (type(D) != "t_INT" || D < 2, return(-1));
  u = sqrtint(D); v = D-u^2;
  if (!v, return(0));
  s = v;
  r = u; j = 0;
  until (u == r && v == s,
    u = (r+u)\v * v - u;
    v = (D-u^2)\v; j++;
  ); j
}
