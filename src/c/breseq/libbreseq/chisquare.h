#include "common.h"

namespace breseq {

  double gamma(double x);
  double pchisq(double df, double x);
  
	double chisquaredistribution(double v, double x);
	double incompletegamma(double a, double x, bool complemented = false);
	double lngamma(double x, double* sgngam);
  
  double nbdtrc(int k, int n, double p);
  double nbdtr(int k, int n, double p);
  double nbdtri(int k, int n, double p);
  
  double incbi(double aa, double bb, double yy0);
  double incbet(double aa, double bb, double xx);
  double incbcf(double a, double b, double x);
  double incbd(double a, double b, double x);
  double pseries(double a, double b, double x);
  
  double ndtri(double y0);

} // namespace breseq
