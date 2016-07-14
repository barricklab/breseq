#ifndef _CHISQUARE_H_
#define _CHISQUARE_H_

#include "common.h"

namespace breseq {

  double gamma(double x);
  double pchisq(double df, double x);
  
	double chisquaredistribution(double v, double x);
	double incompletegamma(double a, double x, bool complemented = false);
	double lngamma(double x, double* sgngam);
  
  double nbdtrc(double k, double n, double p);
  double nbdtr(double k, double n, double p);
  double nbdtri(double k, double n, double p);
  
  double incbi(double aa, double bb, double yy0);
  double incbet(double aa, double bb, double xx);
  double incbcf(double a, double b, double x);
  double incbd(double a, double b, double x);
  double pseries(double a, double b, double x);
  
  double ndtri(double y0);
  
  double bdtrc(double k, double n, double p);
  double bdtr(double k, double n, double p);
  double bdtri(double k, double n, double y);
  
  inline double combination (int32_t num, int32_t choose)
  {
    
    double log_result = 0.0;
    for (int32_t i=num; i > choose; i--)
    {
      log_result += log(i);
    }
    for (int32_t i=2; i <= num - choose; i++)
    {
      log_result -= log(i);
    }
    
    return exp(log_result);
  }
  
  // Uses fast+accurate approximation for factorials
  inline double fast_combination(uint32_t num, uint32_t choose) {
    double sign;
    double result = exp(lngamma(num+1, &sign) - lngamma(choose+1, &sign) - lngamma(num-choose+1, &sign));
    return result;
  }
  
  // probability of exactly this any successes
  inline double binomial (double pr_success, int32_t num_trials, int32_t num_successes)
  {
    double ret_val = fast_combination(num_trials, num_successes) * pow(pr_success, num_successes) * pow(1-pr_success, num_trials-num_successes);
    return ret_val;
  }
  
  
  // probability of this or fewer successes
  inline uint32_t qbinomial (double tail_value, int32_t num_trials, double pr_success)
  {
    ASSERT((tail_value >= 0) && (tail_value <= 1), "probability out of range");
    double cumulative_pr = 0.0;
    
    int32_t num_successes;
    for (num_successes=0; num_successes < num_trials; num_successes++) {
      cumulative_pr += binomial(pr_success, num_trials, num_successes);
      if (cumulative_pr > tail_value) break;
    }
    
    return num_successes;
  }
  
} // namespace breseq

#endif