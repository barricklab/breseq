#include "breseq/common.h"

namespace breseq {

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

		assert(x >= 0 && v >= 1); // Domain error in ChiSquareDistribution
		result = incompletegamma(v/2.0, x/2.0, false);
		return result;
	}

} // namespace breseq
