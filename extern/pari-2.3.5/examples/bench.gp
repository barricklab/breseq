{
  u=v=p=q=1;
  for (k=1, 1000,
    w = u+v; u = v; v = w;
    p *= w; q = lcm(q,w);
    if (k%50 == 0,
      print(k " " log(p)/log(q))
    )
  )
}
