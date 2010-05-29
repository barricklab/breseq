lucas(p) =
{
  local(u,q); u=4; q=1<<p - 1;
  for(k=3,p, u = (sqr(u)-2) % q);
  u == 0
}
