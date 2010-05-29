\\ ---------------  GP code  ---------------------------------------
\\
\\ Time-stamp: <Fri, Mar 26, 1999 - 14:13:17 - villegas@linux47>
\\
\\ Description: Compute class number of imaginary quadratic field
\\ analytically
\\
\\ File: classno.gp
\\
\\ Original Author: Fernando Rodriguez-Villegas 
\\                  villegas@math.utexas.edu
\\                  University of Texas at Austin
\\
\\ Created:         Fri Mar 26 1999
\\-----------------------------------------------------------------

\\ Class number h(-d), -d fundamental.
\\ Adjust constant cc for accuracy, default at least 9 decimal places.

cl(d,cc) =
{
  local(q0,sd,t,c, s = 0, q = 1);

  if (!isfundamental(-d), error("Discriminant not fundamental"));
  if (!cc, cc = 5);
  sd = sqrt(d);
  q0 = exp(-2*Pi/sd); c = -4*Pi/sd;
  for (n=1, ceil(sd*cc),
    q *= q0; t = 1/(1-q);
    s += kronecker(-d,n) * q * t * (1 + c*t*n)
  );
  if (d==3, s *= 3);
  if (d==4, s *= 2);
  -2*s
}
