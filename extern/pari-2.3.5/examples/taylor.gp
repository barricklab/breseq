\\ originally contributed by Ilya Zakharevich

\\ sample function
f(x) = sin(x)

\\ plot Taylor polynomials of f of index first + i*step <= ordlim
\\ for x in [xmin,xmax].
plot_taylor(xmin=-5, xmax=5, ordlim=16, first=1, step=1) = 
{
  local(T,Taylor_array,s,t,w,h,dw,dh,cw,ch,gh,h1, extrasize = 0.6);

  default(seriesprecision,ordlim+1);
  ordlim -= first; ordlim \= step; ordlim += first;
  T = f(tt); Taylor_array = vector(ordlim+1);
  forstep(i=ordlim+1, 1, -1,
    T += O(tt^(1 + first + (i-1)*step));
    Taylor_array[i] = truncate(T)
  );

  t = plothsizes();
  w=floor(t[1]*0.9)-2; dw=floor(t[1]*0.05)+1; cw=t[5];
  h=floor(t[2]*0.9)-2; dh=floor(t[2]*0.05)+1; ch=t[6];
  h1=floor(h/1.2);

  plotinit(2, w+2*dw, h+2*dh);
  plotinit(3, w, h1);
  s = plotrecth(3, x=xmin,xmax, f(x), 2+8+16+32);
  gh=s[4]-s[3];

  plotinit(3, w, h);
  plotscale(3,s[1],s[2],s[3]-gh*extrasize/2,s[4]+gh*extrasize/2);
  plotrecth(3, x=xmin,xmax,
    concat(f(x), subst(Taylor_array, tt, x)), 4);
  plotclip(3);
  plotcopy(3, 2, dw, dh);

  plotmove(2,floor(dw+w/2-15*cw),floor(dh/2));
  plotstring(2,"Multiple Taylor Approximations");
  plotdraw([2, 0, 0])
}
