#include "common.h"

namespace breseq {

  double gamma(double x);
  double pchisq(double df, double x);
  
	double chisquaredistribution(double v, double x);
	double incompletegamma(double a, double x, bool complemented = false);
	double lngamma(double x, double* sgngam);

} // namespace breseq
