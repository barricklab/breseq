#include "libbreseq/common.h"
#include "libbreseq/chisquare.h"

namespace breseq {

  
  /*							gamma.c
   *
   *	Gamma function
   *
   *
   *
   * SYNOPSIS:
   *
   * double x, y, gamma();
   * extern int sgngam;
   *
   * y = gamma( x );
   *
   *
   *
   * DESCRIPTION:
   *
   * Returns gamma function of the argument.  The result is
   * correctly signed, and the sign (+1 or -1) is also
   * returned in a global (extern) variable named sgngam.
   * This variable is also filled in by the logarithmic gamma
   * function lgam().
   *
   * Arguments |x| <= 34 are reduced by recurrence and the function
   * approximated by a rational function of degree 6/7 in the
   * interval (2,3).  Large arguments are handled by Stirling's
   * formula. Large negative arguments are made positive using
   * a reflection formula.  
   *
   *
   * ACCURACY:
   *
   *                      Relative error:
   * arithmetic   domain     # trials      peak         rms
   *    DEC      -34, 34      10000       1.3e-16     2.5e-17
   *    IEEE    -170,-33      20000       2.3e-15     3.3e-16
   *    IEEE     -33,  33     20000       9.4e-16     2.2e-16
   *    IEEE      33, 171.6   20000       2.3e-15     3.2e-16
   *
   * Error for arguments outside the test range will be larger
   * owing to error amplification by the exponential function.
   *
   */
  /*							lgam()
   *
   *	Natural logarithm of gamma function
   *
   *
   *
   * SYNOPSIS:
   *
   * double x, y, lgam();
   * extern int sgngam;
   *
   * y = lgam( x );
   *
   *
   *
   * DESCRIPTION:
   *
   * Returns the base e (2.718...) logarithm of the absolute
   * value of the gamma function of the argument.
   * The sign (+1 or -1) of the gamma function is returned in a
   * global (extern) variable named sgngam.
   *
   * For arguments greater than 13, the logarithm of the gamma
   * function is approximated by the logarithmic version of
   * Stirling's formula using a polynomial approximation of
   * degree 4. Arguments between -33 and +33 are reduced by
   * recurrence to the interval [2,3] of a rational approximation.
   * The cosecant reflection formula is employed for arguments
   * less than -33.
   *
   * Arguments greater than MAXLGM return MAXNUM and an error
   * message.  MAXLGM = 2.035093e36 for DEC
   * arithmetic or 2.556348e305 for IEEE arithmetic.
   *
   *
   *
   * ACCURACY:
   *
   *
   * arithmetic      domain        # trials     peak         rms
   *    DEC     0, 3                  7000     5.2e-17     1.3e-17
   *    DEC     2.718, 2.035e36       5000     3.9e-17     9.9e-18
   *    IEEE    0, 3                 28000     5.4e-16     1.1e-16
   *    IEEE    2.718, 2.556e305     40000     3.5e-16     8.3e-17
   * The error criterion was relative when the function magnitude
   * was greater than one but absolute when it was less than one.
   *
   * The following test used the relative error criterion, though
   * at certain points the relative error could be much higher than
   * indicated.
   *    IEEE    -200, -4             10000     4.8e-16     1.3e-16
   *
   */
  
  /*							gamma.c	*/
  /*	gamma function	*/
  
  /*
   Cephes Math Library Release 2.2:  July, 1992
   Copyright 1984, 1987, 1989, 1992 by Stephen L. Moshier
   Direct inquiries to 30 Frost Street, Cambridge, MA 02140
   */
  
  
//#include "mconf.h"
// USING UNK for DEFINES

  static double PI=
    3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706
    ;
  
  static double P[] = {
    1.60119522476751861407E-4,
    1.19135147006586384913E-3,
    1.04213797561761569935E-2,
    4.76367800457137231464E-2,
    2.07448227648435975150E-1,
    4.94214826801497100753E-1,
    9.99999999999999996796E-1
  };
  static double Q[] = {
    -2.31581873324120129819E-5,
    5.39605580493303397842E-4,
    -4.45641913851797240494E-3,
    1.18139785222060435552E-2,
    3.58236398605498653373E-2,
    -2.34591795718243348568E-1,
    7.14304917030273074085E-2,
    1.00000000000000000320E0
  };
#define MAXGAM 171.624376956302725
  static double LOGPI = 1.14472988584940017414;
  
  /* Stirling's formula for the gamma function */
  static double STIR[5] = {
    7.87311395793093628397E-4,
    -2.29549961613378126380E-4,
    -2.68132617805781232825E-3,
    3.47222221605458667310E-3,
    8.33333333333482257126E-2,
  };
#define MAXSTIR 143.01608
  static double SQTPI = 2.50662827463100050242E0;
  static double LS2PI  =  0.91893853320467274178;
  
#define MAXLGM 2.556348e305
  
  /* A[]: Stirling's formula expansion of log gamma
   * B[], C[]: log gamma function between 2 and 3
   */
  static double A[] = {
    8.11614167470508450300E-4,
    -5.95061904284301438324E-4,
    7.93650340457716943945E-4,
    -2.77777777730099687205E-3,
    8.33333333333331927722E-2
  };
  static double B[] = {
    -1.37825152569120859100E3,
    -3.88016315134637840924E4,
    -3.31612992738871184744E5,
    -1.16237097492762307383E6,
    -1.72173700820839662146E6,
    -8.53555664245765465627E5
  };
  static double C[] = {
    /* 1.00000000000000000000E0, */
    -3.51815701436523470549E2,
    -1.70642106651881159223E4,
    -2.20528590553854454839E5,
    -1.13933444367982507207E6,
    -2.53252307177582951285E6,
    -2.01889141433532773231E6
  };
  static double P0[5] = {
    -5.99633501014107895267E1,
    9.80010754185999661536E1,
    -5.66762857469070293439E1,
    1.39312609387279679503E1,
    -1.23916583867381258016E0,
  };
  static double Q0[8] = {
    /* 1.00000000000000000000E0,*/
    1.95448858338141759834E0,
    4.67627912898881538453E0,
    8.63602421390890590575E1,
    -2.25462687854119370527E2,
    2.00260212380060660359E2,
    -8.20372256168333339912E1,
    1.59056225126211695515E1,
    -1.18331621121330003142E0,
  };
  static double P1[9] = {
    4.05544892305962419923E0,
    3.15251094599893866154E1,
    5.71628192246421288162E1,
    4.40805073893200834700E1,
    1.46849561928858024014E1,
    2.18663306850790267539E0,
    -1.40256079171354495875E-1,
    -3.50424626827848203418E-2,
    -8.57456785154685413611E-4,
  };
  static double Q1[8] = {
    /*  1.00000000000000000000E0,*/
    1.57799883256466749731E1,
    4.53907635128879210584E1,
    4.13172038254672030440E1,
    1.50425385692907503408E1,
    2.50464946208309415979E0,
    -1.42182922854787788574E-1,
    -3.80806407691578277194E-2,
    -9.33259480895457427372E-4,
  };
  static double P2[9] = {
    3.23774891776946035970E0,
    6.91522889068984211695E0,
    3.93881025292474443415E0,
    1.33303460815807542389E0,
    2.01485389549179081538E-1,
    1.23716634817820021358E-2,
    3.01581553508235416007E-4,
    2.65806974686737550832E-6,
    6.23974539184983293730E-9,
  };
  static double Q2[8] = {
    /*  1.00000000000000000000E0,*/
    6.02427039364742014255E0,
    3.67983563856160859403E0,
    1.37702099489081330271E0,
    2.16236993594496635890E-1,
    1.34204006088543189037E-2,
    3.28014464682127739104E-4,
    2.89247864745380683936E-6,
    6.79019408009981274425E-9,
  };
  
  static double MAXLOG =  7.09782712893383996732E2;
  static double MINLOG = -7.451332191019412076235E2;
  static double MACHEP =  1.11022302462515654042E-16;
  static double MAXNUM =  1.79769313486231570815E308;
  
  static int sgngam = 0;  
  static double big = 4.503599627370496e15;
  static double biginv =  2.22044604925031308085e-16;
  static double s2pi = 2.50662827463100050242E0;

  double polevl(double x, double coef[], int N )
  {
    double ans;
    int i;
    double *p;
    
    p = coef;
    ans = *p++;
    i = N;
    
    do
      ans = ans * x  +  *p++;
      while( --i );
    
    return( ans );
  }
  
  /* Gamma function computed by Stirling's formula.
   * The polynomial STIR is valid for 33 <= x <= 172.
   */
  static double stirf(double x)
  {
    double y, w, v;
    
    w = 1.0/x;
    w = 1.0 + w * polevl( w, STIR, 4 );
    y = exp(x);
    if( x > MAXSTIR )
    { /* Avoid overflow in pow() */
      v = pow( x, 0.5 * x - 0.25 );
      y = v * (v / y);
    }
    else
    {
      y = pow( x, x - 0.5 ) / y;
    }
    y = SQTPI * y * w;
    return( y );
  }
  
  /*							p1evl()	*/
  /*                                          N
   * Evaluate polynomial when coefficient of x  is 1.0.
   * Otherwise same as polevl.
   */  
  double p1evl( double x, double coef[], int N )
  {
    double ans;
    double *p;
    int i;
    
    p = coef;
    ans = x + *p++;
    i = N-1;
    
    do
      ans = ans * x  + *p++;
      while( --i );
    
    return( ans );
  }  
  
  double gamma(double x)
  {
    double p, q, z;
    int i;
    
    sgngam = 1;
    
    q = fabs(x);
    
    if( q > 33.0 )
    {
      if( x < 0.0 )
      {
        p = floor(q);
        if( p == q )
        {
          return (NAN);
        }
        i = static_cast<uint32_t>(p);
        if( (i & 1) == 0 )
          sgngam = -1;
          z = q - p;
          if( z > 0.5 )
          {
            p += 1.0;
            z = q - p;
          }
        z = q * sin( PI * z );
        if( z == 0.0 )
        {
          return( sgngam * INFINITY);
        }
        z = fabs(z);
        z = PI/(z * stirf(q) );
      }
      else
      {
        z = stirf(x);
      }
      return( sgngam * z );
    }
    
    z = 1.0;
    while( x >= 3.0 )
    {
      x -= 1.0;
      z *= x;
    }
    
    while( x < 0.0 )
    {
      if( x > -1.E-9 )
        goto small;
      z /= x;
      x += 1.0;
    }
    
    while( x < 2.0 )
    {
      if( x < 1.e-9 )
        goto small;
      z /= x;
      x += 1.0;
    }
    
    if( x == 2.0 )
      return(z);
    
    x -= 2.0;
    p = polevl( x, P, 6 );
    q = polevl( x, Q, 7 );
    return( z * p / q );
    
  small:
    if( x == 0.0 )
    {
      return( INFINITY );
    }
    else
      return( z/((1.0 + 0.5772156649015329 * x) * x) );
  }

  
  double pchisq(double df, double x) {
    
    double pr = 1 / (exp((df/2.0)*log(2.0))*gamma(df/2.0)) * exp(((df/2.0)-1) * log(x)) * exp(-x/2.0);
    return pr;
  }
  
  /* Logarithm of gamma function */  
  double lgam(double x)
  {
    double p, q, u, w, z;
    int i;
    
    sgngam = 1;
    
    if( x < -34.0 )
    {
      q = -x;
      w = lgam(q); /* note this modifies sgngam! */
      p = floor(q);
      if( p == q )
      {
      lgsing:
        return (INFINITY);
      }
      i = static_cast<size_t>(p);
      if( (i & 1) == 0 )
        sgngam = -1;
        else
          sgngam = 1;
          z = q - p;
          if( z > 0.5 )
          {
            p += 1.0;
            z = p - q;
          }
      z = q * sin( PI * z );
      if( z == 0.0 )
        goto lgsing;
      /*	z = log(PI) - log( z ) - w;*/
      z = LOGPI - log( z ) - w;
      return( z );
    }
    
    if( x < 13.0 )
    {
      z = 1.0;
      p = 0.0;
      u = x;
      while( u >= 3.0 )
      {
        p -= 1.0;
        u = x + p;
        z *= u;
      }
      while( u < 2.0 )
      {
        if( u == 0.0 )
          goto lgsing;
        z /= u;
        p += 1.0;
        u = x + p;
      }
      if( z < 0.0 )
      {
        sgngam = -1;
        z = -z;
      }
      else
        sgngam = 1;
        if( u == 2.0 )
          return( log(z) );
      p -= 2.0;
      x = x + p;
      p = x * polevl( x, B, 5 ) / p1evl( x, C, 6);
      return( log(z) + p );
    }
    
    if( x > MAXLGM )
    {
      return( sgngam * INFINITY );
    }
    
    q = ( x - 0.5 ) * log(x) - x + LS2PI;
    if( x > 1.0e8 )
      return( q );
    
    p = 1.0/(x*x);
    if( x >= 1000.0 )
      q += ((   7.9365079365079365079365e-4 * p
             - 2.7777777777777777777778e-3) *p
            + 0.0833333333333333333333) / x;
      else
        q += polevl( p, A, 4 ) / x;
        return( q );
  }

  
  
/////////////////////  
  
	/*************************************************************************
	Natural logarithm of gamma function

	Input parameters:
		X       -   argument

	Result:
		logarithm of the absolute value of the Gamma(X).

	Output parameters:
		SgnGam  -   sign(Gamma(X))

	Domain:
		0 < X < 2.55e305
		-2.55e305 < X < 0, X is not an integer.

	ACCURACY:
	arithmetic      domain        # trials     peak         rms
	   IEEE    0, 3                 28000     5.4e-16     1.1e-16
	   IEEE    2.718, 2.556e305     40000     3.5e-16     8.3e-17
	The error criterion was relative when the function magnitude
	was greater than one but absolute when it was less than one.

	The following test used the relative error criterion, though
	at certain points the relative error could be much higher than
	indicated.
	   IEEE    -200, -4             10000     4.8e-16     1.3e-16

	Cephes Math Library Release 2.8:  June, 2000
	Copyright 1984, 1987, 1989, 1992, 2000 by Stephen L. Moshier
	Translated to AlgoPascal by Bochkanov Sergey (2005, 2006, 2007).
	*************************************************************************/
	double lngamma(double x, double* sgngam)
	{
		double a;
		double b;
		double c;
		double p;
		double q;
		double u;
		double w;
		double z;
		int32_t i;
		double logpi;
		double ls2pi;
		double tmp;
		double result;

		*sgngam = 0;

		*sgngam = 1;
		logpi = 1.14472988584940017414;
		ls2pi = 0.91893853320467274178;
		if( x < -34.0 )
		{
			q = -x;
			w = lngamma(q, &tmp);
			p = (int32_t)floor(q);
			i = (int32_t)floor(p+0.5);
			if( i%2==0 )
			{
				*sgngam = -1;
			}
			else
			{
				*sgngam = 1;
			}
			z = q-p;
			if( z > 0.5 )
			{
				p = p+1;
				z = p-q;
			}
			z = q*sin(M_PI*z);
			result = logpi-log(z)-w;
			return result;
		}
		if( x < 13 )
		{
			z = 1;
			p = 0;
			u = x;
			while(u >= 3)
			{
				p = p-1;
				u = x+p;
				z = z*u;
			}
			while(u < 2)
			{
				z = z/u;
				p = p+1;
				u = x+p;
			}
			if( z < 0 )
			{
				*sgngam = -1;
				z = -z;
			}
			else
			{
				*sgngam = 1;
			}
			if( u == 2 )
			{
				result = log(z);
				return result;
			}
			p = p-2;
			x = x+p;
			b = -1378.25152569120859100;
			b = -38801.6315134637840924+x*b;
			b = -331612.992738871184744+x*b;
			b = -1162370.97492762307383+x*b;
			b = -1721737.00820839662146+x*b;
			b = -853555.664245765465627+x*b;
			c = 1;
			c = -351.815701436523470549+x*c;
			c = -17064.2106651881159223+x*c;
			c = -220528.590553854454839+x*c;
			c = -1139334.44367982507207+x*c;
			c = -2532523.07177582951285+x*c;
			c = -2018891.41433532773231+x*c;
			p = x*b/c;
			result = log(z)+p;
			return result;
		}
		q = (x-0.5)*log(x)-x+ls2pi;
		if( x > 100000000 )
		{
			result = q;
			return result;
		}
		p = 1/(x*x);
		if( x >= 1000.0 )
		{
			q = q+((7.9365079365079365079365*0.0001*p-2.7777777777777777777778*0.001)*p+0.0833333333333333333333)/x;
		}
		else
		{
			a = 8.11614167470508450300*0.0001;
			a = -5.95061904284301438324*0.0001+p*a;
			a = 7.93650340457716943945*0.0001+p*a;
			a = -2.77777777730099687205*0.001+p*a;
			a = 8.33333333333331927722*0.01+p*a;
			q = q+a/x;
		}
		result = q;
		return result;

	}

	/*************************************************************************
	Incomplete gamma integral

	The function is defined by

							  x
							   -
					  1       | |  -t  a-1
	 igam(a,x)  =   -----     |   e   t   dt.
					 -      | |
					| (a)    -
							  0


	In this implementation both arguments must be positive.
	The integral is evaluated by either a power series or
	continued fraction expansion, depending on the relative
	values of a and x.

	ACCURACY:

						 Relative error:
	arithmetic   domain     # trials      peak         rms
	   IEEE      0,30       200000       3.6e-14     2.9e-15
	   IEEE      0,100      300000       9.9e-14     1.5e-14

	Cephes Math Library Release 2.8:  June, 2000
	Copyright 1985, 1987, 2000 by Stephen L. Moshier
	*************************************************************************/
	/*************************************************************************
	Complemented incomplete gamma integral

	The function is defined by


	 igamc(a,x)   =   1 - igam(a,x)

							   inf.
								 -
						1       | |  -t  a-1
				  =   -----     |   e   t   dt.
					   -      | |
					  | (a)    -
								x


	In this implementation both arguments must be positive.
	The integral is evaluated by either a power series or
	continued fraction expansion, depending on the relative
	values of a and x.

	ACCURACY:

	Tested at random a, x.
				   a         x                      Relative error:
	arithmetic   domain   domain     # trials      peak         rms
	   IEEE     0.5,100   0,100      200000       1.9e-14     1.7e-15
	   IEEE     0.01,0.5  0,100      200000       1.4e-13     1.6e-15

	Cephes Math Library Release 2.8:  June, 2000
	Copyright 1985, 1987, 2000 by Stephen L. Moshier
	*************************************************************************/
	double incompletegamma(double a, double x, bool complemented)
	{
		double igammaepsilon;
		double ans;
		double ax;
		double c;
		double r;
		double tmp;
		double result;

		if (!complemented)
		{
			igammaepsilon = 0.000000000000001;
			if( x <= 0||a <= 0 )
			{
				result = 0;
				return result;
			}
			if( x > 1&&x > a )
			{
				result = 1-incompletegamma(a, x, true);
				return result;
			}
			ax = a*log(x)-x-lngamma(a, &tmp);
			if( ax < -709.78271289338399 )
			{
				result = 0;
				return result;
			}
			ax = exp(ax);
			r = a;
			c = 1;
			ans = 1;
			do
			{
				r = r+1;
				c = c*x/r;
				ans = ans+c;
			}
			while(c/ans > igammaepsilon);
			result = ans*ax/a;
		}
		else
		{
			double igammabignumber;
			double igammabignumberinv;
			double yc;
			double t;
			double y;
			double z;
			double pk;
			double pkm1;
			double pkm2;
			double qk;
			double qkm1;
			double qkm2;

			igammaepsilon = 0.000000000000001;
			igammabignumber = 4503599627370496.0;
			igammabignumberinv = 2.22044604925031308085*0.0000000000000001;
			if( x <= 0||a <= 0 )
			{
				result = 1;
				return result;
			}
			if( x < 1||x < a )
			{
				result = 1-incompletegamma(a, x, false);
				return result;
			}
			ax = a*log(x)-x-lngamma(a, &tmp);
			if( ax < -709.78271289338399 )
			{
				result = 0;
				return result;
			}
			ax = exp(ax);
			y = 1-a;
			z = x+y+1;
			c = 0;
			pkm2 = 1;
			qkm2 = x;
			pkm1 = x+1;
			qkm1 = z*x;
			ans = pkm1/qkm1;
			do
			{
				c = c+1;
				y = y+1;
				z = z+2;
				yc = y*c;
				pk = pkm1*z-pkm2*yc;
				qk = qkm1*z-qkm2*yc;
				if( qk != 0 )
				{
					r = pk/qk;
					t = fabs((ans-r)/r);
					ans = r;
				}
				else
				{
					t = 1;
				}
				pkm2 = pkm1;
				pkm1 = pk;
				qkm2 = qkm1;
				qkm1 = qk;
				if( fabs(pk) > igammabignumber )
				{
					pkm2 = pkm2*igammabignumberinv;
					pkm1 = pkm1*igammabignumberinv;
					qkm2 = qkm2*igammabignumberinv;
					qkm1 = qkm1*igammabignumberinv;
				}
			}
			while(t > igammaepsilon);
			result = ans*ax;
		}
		return result;
	}


	/*************************************************************************
	Chi-square distribution

	Returns the area under the left hand tail (from 0 to x)
	of the Chi square probability density function with
	v degrees of freedom.


									  x
									   -
						   1          | |  v/2-1  -t/2
	 P( x | v )   =   -----------     |   t      e     dt
					   v/2  -       | |
					  2    | (v/2)   -
									  0

	where x is the Chi-square variable.

	The incomplete gamma integral is used, according to the
	formula

	y = chdtr( v, x ) = igam( v/2.0, x/2.0 ).

	The arguments must both be positive.

	ACCURACY:

	See incomplete gamma function


	Cephes Math Library Release 2.8:  June, 2000
	Copyright 1984, 1987, 2000 by Stephen L. Moshier
	*************************************************************************/
	double chisquaredistribution(double v, double x)
	{
		double result;

		ASSERT((x >= 0) && (v >= 1), "Domain error in ChiSquareDistribution");
		result = incompletegamma(v/2.0, x/2.0, false);
		return result;
	}
  
  /*							nbdtr.c
   *
   *	Negative binomial distribution
   *
   *
   *
   * SYNOPSIS:
   *
   * int k, n;
   * double p, y, nbdtr();
   *
   * y = nbdtr( k, n, p );
   *
   * DESCRIPTION:
   *
   * Returns the sum of the terms 0 through k of the negative
   * binomial distribution:
   *
   *   k
   *   --  ( n+j-1 )   n      j
   *   >   (       )  p  (1-p)
   *   --  (   j   )
   *  j=0
   *
   * In a sequence of Bernoulli trials, this is the probability
   * that k or fewer failures precede the nth success.
   *
   * The terms are not computed individually; instead the incomplete
   * beta integral is employed, according to the formula
   *
   * y = nbdtr( k, n, p ) = incbet( n, k+1, p ).
   *
   * The arguments must be positive, with p ranging from 0 to 1.
   *
   * ACCURACY:
   *
   * Tested at random points (a,b,p), with p between 0 and 1.
   *
   *               a,b                     Relative error:
   * arithmetic  domain     # trials      peak         rms
   *    IEEE     0,100       100000      1.7e-13     8.8e-15
   * See also incbet.c.
   *
   * @JEB This is equivalent to pnbinom(10, size=5, prob=0.3, lower.tail=F) in R.
   *      
   * Cephes Math Library Release 2.8:  June, 2000
   * Copyright 1984, 1987, 1995, 2000 by Stephen L. Moshier
   */  
  double nbdtrc(int k, double n, double p)
  {
    double dk, dn;
    
    ASSERT( (p >= 0.0) && (p <= 1.0) &&  (k >= 0), "Domain error in nbdtrc" );
    
    dk = k+1;
    dn = n;
    return( incbet( dk, dn, 1.0 - p ) );
  }  
  
  /*							nbdtr.c
   *
   *	Complemented negative binomial distribution
   *
   *
   *
   * SYNOPSIS:
   *
   * int k, n;
   * double p, y, nbdtrc();
   *
   * y = nbdtrc( k, n, p );
   *
   * DESCRIPTION:
   *
   * Returns the sum of the terms k+1 to infinity of the negative
   * binomial distribution:
   *
   *   inf
   *   --  ( n+j-1 )   n      j
   *   >   (       )  p  (1-p)
   *   --  (   j   )
   *  j=k+1
   *
   * The terms are not computed individually; instead the incomplete
   * beta integral is employed, according to the formula
   *
   * y = nbdtrc( k, n, p ) = incbet( k+1, n, 1-p ).
   *
   * The arguments must be positive, with p ranging from 0 to 1.
   *
   * ACCURACY:
   *
   * Tested at random points (a,b,p), with p between 0 and 1.
   *
   *               a,b                     Relative error:
   * arithmetic  domain     # trials      peak         rms
   *    IEEE     0,100       100000      1.7e-13     8.8e-15
   * See also incbet.c.
   *      
   * Cephes Math Library Release 2.8:  June, 2000
   * Copyright 1984, 1987, 1995, 2000 by Stephen L. Moshier
   */  
  double nbdtr(int k, double n, double p)
  {
    double dk, dn;
    
    ASSERT( (p >= 0.0) && (p <= 1.0) &&  (k >= 0), "Domain error in nbdtr" );
    
    dk = k+1;
    dn = n;
    return( incbet( dn, dk, p ) );
  }  
  
  /*							nbdtr.c
   *
   *	Functional inverse of negative binomial distribution
   *
   *
   *
   * SYNOPSIS:
   *
   * int k, n;
   * double p, y, nbdtri();
   *
   * p = nbdtri( k, n, y );
   *
   * DESCRIPTION:
   *
   * Finds the argument p such that nbdtr(k,n,p) is equal to y.
   *
   * ACCURACY:
   *
   * Tested at random points (a,b,y), with y between 0 and 1.
   *
   *               a,b                     Relative error:
   * arithmetic  domain     # trials      peak         rms
   *    IEEE     0,100       100000      1.5e-14     8.5e-16
   * See also incbi.c.
   *      
   * Cephes Math Library Release 2.8:  June, 2000
   * Copyright 1984, 1987, 1995, 2000 by Stephen L. Moshier
   */  
  double nbdtri(int k, double n, double p)
  {
    double dk, dn, w;
        
    ASSERT( (p >= 0.0) && (p <= 1.0) && (k >= 0), "Domain error in nbdtri" );
    
    dk = k+1;
    dn = n;
    w = incbi( dn, dk, p );
    return( w );
  }
  
  /*							incbi()
   *
   *      Inverse of imcomplete beta integral
   *
   *
   *
   * SYNOPSIS:
   *
   * double a, b, x, y, incbi();
   *
   * x = incbi( a, b, y );
   *
   *
   *
   * DESCRIPTION:
   *
   * Given y, the function finds x such that
   *
   *  incbet( a, b, x ) = y .
   *
   * The routine performs interval halving or Newton iterations to find the
   * root of incbet(a,b,x) - y = 0.
   *
   *
   * ACCURACY:
   *
   *                      Relative error:
   *                x     a,b
   * arithmetic   domain  domain  # trials    peak       rms
   *    IEEE      0,1    .5,10000   50000    5.8e-12   1.3e-13
   *    IEEE      0,1   .25,100    100000    1.8e-13   3.9e-15
   *    IEEE      0,1     0,5       50000    1.1e-12   5.5e-15
   *    VAX       0,1    .5,100     25000    3.5e-14   1.1e-15
   * With a and b constrained to half-integer or integer values:
   *    IEEE      0,1    .5,10000   50000    5.8e-12   1.1e-13
   *    IEEE      0,1    .5,100    100000    1.7e-14   7.9e-16
   * With a = .5, b constrained to half-integer or integer values:
   *    IEEE      0,1    .5,10000   10000    8.3e-11   1.0e-11
   *
   * Cephes Math Library Release 2.8:  June, 2000
   * Copyright 1984, 1996, 2000 by Stephen L. Moshier
   */  
  double incbi(double aa, double bb, double yy0)
  {
    double a, b, y0, d, y, x, x0, x1, lgm, yp, di, dithresh, yl, yh, xt;
    int i, rflg, dir, nflg;
    
    
    i = 0;
    if( yy0 <= 0 )
      return(0.0);
    if( yy0 >= 1.0 )
      return(1.0);
    x0 = 0.0;
    yl = 0.0;
    x1 = 1.0;
    yh = 1.0;
    nflg = 0;
    
    if( aa <= 1.0 || bb <= 1.0 )
    {
      dithresh = 1.0e-6;
      rflg = 0;
      a = aa;
      b = bb;
      y0 = yy0;
      x = a/(a+b);
      y = incbet( a, b, x );
      goto ihalve;
    }
    else
    {
      dithresh = 1.0e-4;
    }
    /* approximation to inverse function */
    
    yp = -ndtri(yy0);
    
    if( yy0 > 0.5 )
    {
      rflg = 1;
      a = bb;
      b = aa;
      y0 = 1.0 - yy0;
      yp = -yp;
    }
    else
    {
      rflg = 0;
      a = aa;
      b = bb;
      y0 = yy0;
    }
    
    lgm = (yp * yp - 3.0)/6.0;
    x = 2.0/( 1.0/(2.0*a-1.0)  +  1.0/(2.0*b-1.0) );
    d = yp * sqrt( x + lgm ) / x
    - ( 1.0/(2.0*b-1.0) - 1.0/(2.0*a-1.0) )
    * (lgm + 5.0/6.0 - 2.0/(3.0*x));
    d = 2.0 * d;
    if( d < MINLOG )
    {
      x = 1.0;
      goto under;
    }
    x = a/( a + b * exp(d) );
    y = incbet( a, b, x );
    yp = (y - y0)/y0;
    if( fabs(yp) < 0.2 )
      goto newt;
    
    /* Resort to interval halving if not close enough. */
  ihalve:
    
    dir = 0;
    di = 0.5;
    for( i=0; i<100; i++ )
    {
      if( i != 0 )
      {
        x = x0  +  di * (x1 - x0);
        if( x == 1.0 )
          x = 1.0 - MACHEP;
          if( x == 0.0 )
          {
            di = 0.5;
            x = x0  +  di * (x1 - x0);
            if( x == 0.0 )
              goto under;
          }
        y = incbet( a, b, x );
        yp = (x1 - x0)/(x1 + x0);
        if( fabs(yp) < dithresh )
          goto newt;
        yp = (y-y0)/y0;
        if( fabs(yp) < dithresh )
          goto newt;
      }
      if( y < y0 )
      {
        x0 = x;
        yl = y;
        if( dir < 0 )
        {
          dir = 0;
          di = 0.5;
        }
        else if( dir > 3 )
          di = 1.0 - (1.0 - di) * (1.0 - di);
          else if( dir > 1 )
            di = 0.5 * di + 0.5; 
            else
              di = (y0 - y)/(yh - yl);
              dir += 1;
              if( x0 > 0.75 )
              {
                if( rflg == 1 )
                {
                  rflg = 0;
                  a = aa;
                  b = bb;
                  y0 = yy0;
                }
                else
                {
                  rflg = 1;
                  a = bb;
                  b = aa;
                  y0 = 1.0 - yy0;
                }
                x = 1.0 - x;
                y = incbet( a, b, x );
                x0 = 0.0;
                yl = 0.0;
                x1 = 1.0;
                yh = 1.0;
                goto ihalve;
              }
      }
      else
      {
        x1 = x;
        if( rflg == 1 && x1 < MACHEP )
        {
          x = 0.0;
          goto done;
        }
        yh = y;
        if( dir > 0 )
        {
          dir = 0;
          di = 0.5;
        }
        else if( dir < -3 )
          di = di * di;
          else if( dir < -1 )
            di = 0.5 * di;
            else
              di = (y - y0)/(yh - yl);
              dir -= 1;
              }
    }
    ASSERT(true, "Percision loss in incbi");
    if( x0 >= 1.0 )
    {
      x = 1.0 - MACHEP;
      goto done;
    }
    if( x <= 0.0 )
    {
    under:
      ASSERT(true, "Underflow in incbi");
      x = 0.0;
      goto done;
    }
    
  newt:
    
    if( nflg )
      goto done;
    nflg = 1;
    lgm = lgam(a+b) - lgam(a) - lgam(b);
    
    for( i=0; i<8; i++ )
    {
      /* Compute the function at this point. */
      if( i != 0 )
        y = incbet(a,b,x);
        if( y < yl )
        {
          x = x0;
          y = yl;
        }
        else if( y > yh )
        {
          x = x1;
          y = yh;
        }
        else if( y < y0 )
        {
          x0 = x;
          yl = y;
        }
        else
        {
          x1 = x;
          yh = y;
        }
      if( x == 1.0 || x == 0.0 )
        break;
      /* Compute the derivative of the function at this point. */
      d = (a - 1.0) * log(x) + (b - 1.0) * log(1.0-x) + lgm;
      if( d < MINLOG )
        goto done;
      if( d > MAXLOG )
        break;
      d = exp(d);
      /* Compute the step to the next approximation of x. */
      d = (y - y0)/d;
      xt = x - d;
      if( xt <= x0 )
      {
        y = (x - x0) / (x1 - x0);
        xt = x0 + 0.5 * y * (x - x0);
        if( xt <= 0.0 )
          break;
      }
      if( xt >= x1 )
      {
        y = (x1 - x) / (x1 - x0);
        xt = x1 - 0.5 * y * (x1 - x);
        if( xt >= 1.0 )
          break;
      }
      x = xt;
      if( fabs(d/x) < 128.0 * MACHEP )
        goto done;
    }
    /* Did not converge.  */
    dithresh = 256.0 * MACHEP;
    goto ihalve;
    
  done:
    
    if( rflg )
    {
      if( x <= MACHEP )
        x = 1.0 - MACHEP;
        else
          x = 1.0 - x;
          }
    return( x );
  }
  
  /*							incbet.c
   *
   *	Incomplete beta integral
   *
   *
   * SYNOPSIS:
   *
   * double a, b, x, y, incbet();
   *
   * y = incbet( a, b, x );
   *
   *
   * DESCRIPTION:
   *
   * Returns incomplete beta integral of the arguments, evaluated
   * from zero to x.  The function is defined as
   *
   *                  x
   *     -            -
   *    | (a+b)      | |  a-1     b-1
   *  -----------    |   t   (1-t)   dt.
   *   -     -     | |
   *  | (a) | (b)   -
   *                 0
   *
   * The domain of definition is 0 <= x <= 1.  In this
   * implementation a and b are restricted to positive values.
   * The integral from x to 1 may be obtained by the symmetry
   * relation
   *
   *    1 - incbet( a, b, x )  =  incbet( b, a, 1-x ).
   *
   * The integral is evaluated by a continued fraction expansion
   * or, when b*x is small, by a power series.
   *
   * ACCURACY:
   *
   * Tested at uniformly distributed random points (a,b,x) with a and b
   * in "domain" and x between 0 and 1.
   *                                        Relative error
   * arithmetic   domain     # trials      peak         rms
   *    IEEE      0,5         10000       6.9e-15     4.5e-16
   *    IEEE      0,85       250000       2.2e-13     1.7e-14
   *    IEEE      0,1000      30000       5.3e-12     6.3e-13
   *    IEEE      0,10000    250000       9.3e-11     7.1e-12
   *    IEEE      0,100000    10000       8.7e-10     4.8e-11
   * Outputs smaller than the IEEE gradual underflow threshold
   * were excluded from these statistics.
   *
   * ERROR MESSAGES:
   *   message         condition      value returned
   * incbet domain      x<0, x>1          0.0
   * incbet underflow                     0.0
   *
   * Cephes Math Library, Release 2.8:  June, 2000
   * Copyright 1984, 1995, 2000 by Stephen L. Moshier
   */
  double incbet(double aa, double bb, double xx)
  {
    double a, b, t, x, xc, w, y;
    int flag;
    
    if( aa <= 0.0 || bb <= 0.0 )
      goto domerr;
    
    if( (xx <= 0.0) || ( xx >= 1.0) )
    {
      if( xx == 0.0 )
        return(0.0);
      if( xx == 1.0 )
        return( 1.0 );
    domerr:
      ASSERT(true, "Domain error in incbet");
      return( 0.0 );
    }
    
    flag = 0;
    if( (bb * xx) <= 1.0 && xx <= 0.95)
    {
      t = pseries(aa, bb, xx);
      goto done;
    }
    
    w = 1.0 - xx;
    
    /* Reverse a and b if x is greater than the mean. */
    if( xx > (aa/(aa+bb)) )
    {
      flag = 1;
      a = bb;
      b = aa;
      xc = xx;
      x = w;
    }
    else
    {
      a = aa;
      b = bb;
      xc = w;
      x = xx;
    }
    
    if( flag == 1 && (b * x) <= 1.0 && x <= 0.95)
    {
      t = pseries(a, b, x);
      goto done;
    }
    
    /* Choose expansion for better convergence. */
    y = x * (a+b-2.0) - (a-1.0);
    if( y < 0.0 )
      w = incbcf( a, b, x );
      else
        w = incbd( a, b, x ) / xc;
        
      /* Multiply w by the factor
       a      b   _             _     _
       x  (1-x)   | (a+b) / ( a | (a) | (b) ) .   */
        
        y = a * log(x);
        t = b * log(xc);
        if( (a+b) < MAXGAM && fabs(y) < MAXLOG && fabs(t) < MAXLOG )
        {
          t = pow(xc,b);
          t *= pow(x,a);
          t /= a;
          t *= w;
          t *= gamma(a+b) / (gamma(a) * gamma(b));
          goto done;
        }
    /* Resort to logarithms.  */
    y += t + lgam(a+b) - lgam(a) - lgam(b);
    y += log(w/a);
    if( y < MINLOG )
      t = 0.0;
      else
        t = exp(y);
        
        done:
        
        if( flag == 1 )
        {
          if( t <= MACHEP )
            t = 1.0 - MACHEP;
            else
              t = 1.0 - t;
              }
    return( t );
  }
  
  /* Continued fraction expansion #1
   * for incomplete beta integral
   */  
  double incbcf(double a, double b, double x)
  {
    double xk, pk, pkm1, pkm2, qk, qkm1, qkm2;
    double k1, k2, k3, k4, k5, k6, k7, k8;
    double r, t, ans, thresh;
    int n;
    
    k1 = a;
    k2 = a + b;
    k3 = a;
    k4 = a + 1.0;
    k5 = 1.0;
    k6 = b - 1.0;
    k7 = k4;
    k8 = a + 2.0;
    
    pkm2 = 0.0;
    qkm2 = 1.0;
    pkm1 = 1.0;
    qkm1 = 1.0;
    ans = 1.0;
    r = 1.0;
    n = 0;
    thresh = 3.0 * MACHEP;
    do
    {
      
      xk = -( x * k1 * k2 )/( k3 * k4 );
      pk = pkm1 +  pkm2 * xk;
      qk = qkm1 +  qkm2 * xk;
      pkm2 = pkm1;
      pkm1 = pk;
      qkm2 = qkm1;
      qkm1 = qk;
      
      xk = ( x * k5 * k6 )/( k7 * k8 );
      pk = pkm1 +  pkm2 * xk;
      qk = qkm1 +  qkm2 * xk;
      pkm2 = pkm1;
      pkm1 = pk;
      qkm2 = qkm1;
      qkm1 = qk;
      
      if( qk != 0 )
        r = pk/qk;
        if( r != 0 )
        {
          t = fabs( (ans - r)/r );
          ans = r;
        }
        else
          t = 1.0;
          
          if( t < thresh )
            goto cdone;
      
      k1 += 1.0;
      k2 += 1.0;
      k3 += 2.0;
      k4 += 2.0;
      k5 += 1.0;
      k6 -= 1.0;
      k7 += 2.0;
      k8 += 2.0;
      
      if( (fabs(qk) + fabs(pk)) > big )
      {
        pkm2 *= biginv;
        pkm1 *= biginv;
        qkm2 *= biginv;
        qkm1 *= biginv;
      }
      if( (fabs(qk) < biginv) || (fabs(pk) < biginv) )
      {
        pkm2 *= big;
        pkm1 *= big;
        qkm2 *= big;
        qkm1 *= big;
      }
    }
    while( ++n < 300 );
    
  cdone:
    return(ans);
  }
  
  
  /* Continued fraction expansion #2
   * for incomplete beta integral
   */  
  double incbd(double a, double b, double x)
  {
    double xk, pk, pkm1, pkm2, qk, qkm1, qkm2;
    double k1, k2, k3, k4, k5, k6, k7, k8;
    double r, t, ans, z, thresh;
    int n;
    
    k1 = a;
    k2 = b - 1.0;
    k3 = a;
    k4 = a + 1.0;
    k5 = 1.0;
    k6 = a + b;
    k7 = a + 1.0;;
    k8 = a + 2.0;
    
    pkm2 = 0.0;
    qkm2 = 1.0;
    pkm1 = 1.0;
    qkm1 = 1.0;
    z = x / (1.0-x);
    ans = 1.0;
    r = 1.0;
    n = 0;
    thresh = 3.0 * MACHEP;
    do
    {
      
      xk = -( z * k1 * k2 )/( k3 * k4 );
      pk = pkm1 +  pkm2 * xk;
      qk = qkm1 +  qkm2 * xk;
      pkm2 = pkm1;
      pkm1 = pk;
      qkm2 = qkm1;
      qkm1 = qk;
      
      xk = ( z * k5 * k6 )/( k7 * k8 );
      pk = pkm1 +  pkm2 * xk;
      qk = qkm1 +  qkm2 * xk;
      pkm2 = pkm1;
      pkm1 = pk;
      qkm2 = qkm1;
      qkm1 = qk;
      
      if( qk != 0 )
        r = pk/qk;
        if( r != 0 )
        {
          t = fabs( (ans - r)/r );
          ans = r;
        }
        else
          t = 1.0;
          
          if( t < thresh )
            goto cdone;
      
      k1 += 1.0;
      k2 -= 1.0;
      k3 += 2.0;
      k4 += 2.0;
      k5 += 1.0;
      k6 += 1.0;
      k7 += 2.0;
      k8 += 2.0;
      
      if( (fabs(qk) + fabs(pk)) > big )
      {
        pkm2 *= biginv;
        pkm1 *= biginv;
        qkm2 *= biginv;
        qkm1 *= biginv;
      }
      if( (fabs(qk) < biginv) || (fabs(pk) < biginv) )
      {
        pkm2 *= big;
        pkm1 *= big;
        qkm2 *= big;
        qkm1 *= big;
      }
    }
    while( ++n < 300 );
  cdone:
    return(ans);
  }
  
  /* Power series for incomplete beta integral.
   Use when b*x is small and x not too close to 1.  */  
  double pseries(double a, double b, double x)
  {
    double s, t, u, v, n, t1, z, ai;
    
    ai = 1.0 / a;
    u = (1.0 - b) * x;
    v = u / (a + 1.0);
    t1 = v;
    t = u;
    n = 2.0;
    s = 0.0;
    z = MACHEP * ai;
    while( fabs(v) > z )
    {
      u = (n - b) * x / n;
      t *= u;
      v = t / (a + n);
      s += v; 
      n += 1.0;
    }
    s += t1;
    s += ai;
    
    u = a * log(x);
    if( (a+b) < MAXGAM && fabs(u) < MAXLOG )
    {
      t = gamma(a+b)/(gamma(a)*gamma(b));
      s = s * t * pow(x,a);
    }
    else
    {
      t = lgam(a+b) - lgam(a) - lgam(b) + u + log(s);
      if( t < MINLOG )
        s = 0.0;
        else
          s = exp(t);
          }
    return(s);
  }
  
  /*							ndtri.c
   *
   *	Inverse of Normal distribution function
   *
   *
   *
   * SYNOPSIS:
   *
   * double x, y, ndtri();
   *
   * x = ndtri( y );
   *
   *
   *
   * DESCRIPTION:
   *
   * Returns the argument, x, for which the area under the
   * Gaussian probability density function (integrated from
   * minus infinity to x) is equal to y.
   *
   *
   * For small arguments 0 < y < exp(-2), the program computes
   * z = sqrt( -2.0 * log(y) );  then the approximation is
   * x = z - log(z)/z  - (1/z) P(1/z) / Q(1/z).
   * There are two rational functions P/Q, one for 0 < y < exp(-32)
   * and the other for y up to exp(-2).  For larger arguments,
   * w = y - 0.5, and  x/sqrt(2pi) = w + w**3 R(w**2)/S(w**2)).
   *
   *
   * ACCURACY:
   *
   *                      Relative error:
   * arithmetic   domain        # trials      peak         rms
   *    DEC      0.125, 1         5500       9.5e-17     2.1e-17
   *    DEC      6e-39, 0.135     3500       5.7e-17     1.3e-17
   *    IEEE     0.125, 1        20000       7.2e-16     1.3e-16
   *    IEEE     3e-308, 0.135   50000       4.6e-16     9.8e-17
   *
   *
   * ERROR MESSAGES:
   *
   *   message         condition    value returned
   * ndtri domain       x <= 0        -MAXNUM
   * ndtri domain       x >= 1         MAXNUM
   *
   * Cephes Math Library Release 2.8:  June, 2000
   * Copyright 1984, 1987, 1989, 2000 by Stephen L. Moshier
   */  
  double ndtri(double y0)
  {
    double x, y, z, y2, x0, x1;
    int code;
    
    if( y0 <= 0.0 )
    {
      ASSERT(true, "Domain error in ndtri");
      return( -MAXNUM );
    }
    if( y0 >= 1.0 )
    {
      ASSERT(true, "Domain error in ndtri");
      return( MAXNUM );
    }
    code = 1;
    y = y0;
    if( y > (1.0 - 0.13533528323661269189) ) /* 0.135... = exp(-2) */
    {
      y = 1.0 - y;
      code = 0;
    }
    
    if( y > 0.13533528323661269189 )
    {
      y = y - 0.5;
      y2 = y * y;
      x = y + y * (y2 * polevl( y2, P0, 4)/p1evl( y2, Q0, 8 ));
      x = x * s2pi; 
      return(x);
    }
    
    x = sqrt( -2.0 * log(y) );
    x0 = x - log(x)/x;
    
    z = 1.0/x;
    if( x < 8.0 ) /* y > exp(-32) = 1.2664165549e-14 */
      x1 = z * polevl( z, P1, 8 )/p1evl( z, Q1, 8 );
      else
        x1 = z * polevl( z, P2, 8 )/p1evl( z, Q2, 8 );
        x = x0 - x1;
        if( code != 0 )
          x = -x;
          return( x );
  }
  
  /*							bdtrc()
   *
   *	Complemented binomial distribution
   *
   *
   *
   * SYNOPSIS:
   *
   * int k, n;
   * double p, y, bdtrc();
   *
   * y = bdtrc( k, n, p );
   *
   * DESCRIPTION:
   *
   * Returns the sum of the terms k+1 through n of the Binomial
   * probability density:
   *
   *   n
   *   --  ( n )   j      n-j
   *   >   (   )  p  (1-p)
   *   --  ( j )
   *  j=k+1
   *
   * The terms are not summed directly; instead the incomplete
   * beta integral is employed, according to the formula
   *
   * y = bdtrc( k, n, p ) = incbet( k+1, n-k, p ).
   *
   * The arguments must be positive, with p ranging from 0 to 1.
   *
   * ACCURACY:
   *
   * Tested at random points (a,b,p).
   *
   *               a,b                     Relative error:
   * arithmetic  domain     # trials      peak         rms
   *  For p between 0.001 and 1:
   *    IEEE     0,100       100000      6.7e-15     8.2e-16
   *  For p between 0 and .001:
   *    IEEE     0,100       100000      1.5e-13     2.7e-15
   *
   * ERROR MESSAGES:
   *
   *   message         condition      value returned
   * bdtrc domain      x<0, x>1, n<k       0.0
   *
   * Cephes Math Library Release 2.8:  June, 2000
   * Copyright 1984, 1987, 1995, 2000 by Stephen L. Moshier
   */
  double bdtrc(int k, double n, double p)
  {
    double dk, dn;
    
    if( (p < 0.0) || (p > 1.0) )
      goto domerr;
    if( k < 0 )
      return( 1.0 );
    
    if( n < k )
    {
    domerr:
      ASSERT(true, "Domain error in bdtrc");
      return( 0.0 );
    }
    
    if( k == n )
      return( 0.0 );
    dn = n - k;
    if( k == 0 )
    {
      if( p < .01 )
        dk = -expm1( dn * log1p(-p) );
        else
          dk = 1.0 - pow( 1.0-p, dn );
          }
    else
    {
      dk = k + 1;
      dk = incbet( dk, dn, p );
    }
    return( dk );
  }
  
  
  /*							bdtr()
   *
   *	Binomial distribution
   *
   *
   *
   * SYNOPSIS:
   *
   * int k, n;
   * double p, y, bdtr();
   *
   * y = bdtr( k, n, p );
   *
   * DESCRIPTION:
   *
   * Returns the sum of the terms 0 through k of the Binomial
   * probability density:
   *
   *   k
   *   --  ( n )   j      n-j
   *   >   (   )  p  (1-p)
   *   --  ( j )
   *  j=0
   *
   * The terms are not summed directly; instead the incomplete
   * beta integral is employed, according to the formula
   *
   * y = bdtr( k, n, p ) = incbet( n-k, k+1, 1-p ).
   *
   * The arguments must be positive, with p ranging from 0 to 1.
   *
   * ACCURACY:
   *
   * Tested at random points (a,b,p), with p between 0 and 1.
   *
   *               a,b                     Relative error:
   * arithmetic  domain     # trials      peak         rms
   *  For p between 0.001 and 1:
   *    IEEE     0,100       100000      4.3e-15     2.6e-16
   * See also incbet.c.
   *
   * ERROR MESSAGES:
   *
   *   message         condition      value returned
   * bdtr domain         k < 0            0.0
   *                     n < k
   *                     x < 0, x > 1
   *
   * Cephes Math Library Release 2.8:  June, 2000
   * Copyright 1984, 1987, 1995, 2000 by Stephen L. Moshier
   */
  double bdtr(int k, double n, double p)
  {
    double dk, dn;
    
    if( (p < 0.0) || (p > 1.0) )
      goto domerr;
    if( (k < 0) || (n < k) )
    {
    domerr:
      ASSERT(true, "Domain error in bdtr");
      return( 0.0 );
    }
    
    if( k == n )
      return( 1.0 );
    
    dn = n - k;
    if( k == 0 )
    {
      dk = pow( 1.0-p, dn );
    }
    else
    {
      dk = k + 1;
      dk = incbet( dn, dk, 1.0 - p );
    }
    return( dk );
  }
  
  /*							bdtri()
   *
   *	Inverse binomial distribution
   *
   *
   *
   * SYNOPSIS:
   *
   * int k, n;
   * double p, y, bdtri();
   *
   * p = bdtr( k, n, y );
   *
   * DESCRIPTION:
   *
   * Finds the event probability p such that the sum of the
   * terms 0 through k of the Binomial probability density
   * is equal to the given cumulative probability y.
   *
   * This is accomplished using the inverse beta integral
   * function and the relation
   *
   * 1 - p = incbi( n-k, k+1, y ).
   *
   * ACCURACY:
   *
   * Tested at random points (a,b,p).
   *
   *               a,b                     Relative error:
   * arithmetic  domain     # trials      peak         rms
   *  For p between 0.001 and 1:
   *    IEEE     0,100       100000      2.3e-14     6.4e-16
   *    IEEE     0,10000     100000      6.6e-12     1.2e-13
   *  For p between 10^-6 and 0.001:
   *    IEEE     0,100       100000      2.0e-12     1.3e-14
   *    IEEE     0,10000     100000      1.5e-12     3.2e-14
   * See also incbi.c.
   *
   * ERROR MESSAGES:
   *
   *   message         condition      value returned
   * bdtri domain     k < 0, n <= k         0.0
   *                  x < 0, x > 1
   *
   * Cephes Math Library Release 2.8:  June, 2000
   * Copyright 1984, 1987, 1995, 2000 by Stephen L. Moshier
   */
  double bdtri(int k, double n, double y)
  {
    double dk, dn, p;
    
    if( (y < 0.0) || (y > 1.0) )
      goto domerr;
    if( (k < 0) || (n <= k) )
    {
    domerr:
      ASSERT(true, "Domain error in ndtri");
      return( 0.0 );
    }
    
    dn = n - k;
    if( k == 0 )
    {
      if( y > 0.8 )
        p = -expm1( log1p(y-1.0) / dn );
        else
          p = 1.0 - pow( y, 1.0/dn );
          }
    else
    {
      dk = k + 1;
      p = incbet( dn, dk, 0.5 );
      if( p > 0.5 )
        p = incbi( dk, dn, 1.0-y );
        else
          p = 1.0 - incbi( dn, dk, y );
          }
    return( p );
  }
  
} // namespace breseq
