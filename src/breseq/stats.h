#ifndef _STATS_H_
#define _STATS_H_

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

  // Negative binomial PMF in R's (size, mu) parametrization.
  // Equivalent to R's dnbinom(k, size=size, mu=mu).
  inline double dnbinom_mu(double k, double size, double mu)
  {
    ASSERT((size > 0) && (mu > 0), "Domain error in dnbinom_mu");
    double p = size / (size + mu);
    double sign;
    double log_pmf = lngamma(k + size, &sign) - lngamma(size, &sign) - lngamma(k + 1, &sign)
                    + size * log(p) + k * log(1.0 - p);
    return exp(log_pmf);
  }

  // Negative binomial quantile function in R's (size, mu) parametrization.
  // Returns the smallest non-negative integer k such that
  // nbdtr(k, size, prob) >= target_pr. Equivalent to R's qnbinom(target_pr, size=size, mu=mu).
  uint32_t qnbinom_mu(double target_pr, double size, double mu);

  // Two-sided Fisher's exact test p-value for a 2x2 contingency table
  //   [ a b ]
  //   [ c d ]
  // Equivalent to R's fisher.test(matrix(c(a,c,b,d), nrow=2), alternative="two.sided")$p.value
  double fisher_exact_test_2x2(uint32_t a, uint32_t b, uint32_t c, uint32_t d);

  // One-sided two-sample Kolmogorov-Smirnov test p-value.
  // Equivalent to R's ks.test(x, y, alternative="less")$p.value
  double ks_test_two_sample_less(const vector<double>& x, const vector<double>& y);

  // Result of a Nelder-Mead simplex minimization.
  struct nelder_mead_result_t {
    vector<double> estimate;
    bool converged;
  };

  // Generic Nelder-Mead simplex minimizer (derivative-free). Replaces R's nlm()
  // for the small, smooth, low-dimensional objective functions used for curve fitting.
  nelder_mead_result_t nelder_mead_minimize(
                                             const function<double(const vector<double>&)>& objective_function,
                                             const vector<double>& initial_guess,
                                             uint32_t max_iterations = 1000,
                                             double tolerance = 1e-8
                                             );

} // namespace breseq

#endif