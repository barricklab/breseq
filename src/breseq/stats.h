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

  // Beta-binomial: an over-dispersed binomial. Parametrized here through the
  // Beta prior's (alpha, beta); for a mean p0 and intra-class correlation rho,
  //   alpha = p0 * (1 - rho) / rho,  beta = (1 - p0) * (1 - rho) / rho
  // which gives Var = n*p0*(1-p0)*[1 + (n-1)*rho], i.e. rho is exactly the
  // variance-inflation-per-extra-trial. rho -> 0 recovers the plain binomial.
  double lbeta(double a, double b);
  double beta_binomial_pmf(double k, double n, double alpha, double beta);
  // Upper tail P(X >= k). Summed directly rather than as 1 - CDF: the lower-tail
  // subtraction cancels catastrophically for exactly the strong signals we care
  // about. At alpha=0.0499, beta=498.95 (p0=1e-4, rho=0.002), k=20 of n=50 is
  // 8.56e-26, but computing it as 1 - CDF yields -3.7e-13 -- a negative probability.
  double beta_binomial_sf(double k, double n, double alpha, double beta);

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

  // ---------------------------------------------------------------------------------------------
  // Exact (Clopper-Pearson) one-sided confidence bounds on a binomial proportion k/n.
  //
  // Used wherever a variant frequency is compared against a cutoff, so the comparison is a
  // CONFIDENCE statement rather than a point estimate: a call is never rejected merely for having
  // small counts, only for being confidently on the wrong side of the cutoff.
  //
  //   lower: L = BetaInv(alpha; k, n-k+1),  L = 0 when k = 0
  //   upper: U = BetaInv(1-alpha; k+1, n-k), U = 1 when k = n
  //
  // Note the identity L(k of k) = alpha^(1/k). Requiring L >= c therefore implies a minimum count
  // k >= ln(alpha)/ln(c) even at 100% frequency -- negligible for a small c (k >= 2 at c = 0.2) but
  // a de facto coverage requirement for a large one (k >= 14 at c = 0.8). Test a HIGH cutoff with
  // the upper bound instead of inverting the lower one.
  //
  // n is a double because some callers average counts from two junction sides; k is a genuine
  // integer count. Both are clamped into range.
  double binomial_frequency_lower_bound(double k, double n, double alpha = 0.05);
  double binomial_frequency_upper_bound(double k, double n, double alpha = 0.05);

  // ---------------------------------------------------------------------------------------------
  // Posterior probability that a read came from a candidate junction rather than from the
  // reference, given the log-likelihood ratio between the two hypotheses.
  //
  //   w = 1 / (1 + prior_odds_against * n_ref_placements * 10^(-log10_odds_for_junction))
  //
  // log10_odds_for_junction is log10 P(read | THIS candidate) - log10 P(read | reference genome),
  // both terms computed by alignment_log10_likelihood() (reference_sequence.cpp) and carried
  // through the BAMs in the X8/X7 tags. Positions where the two references agree contribute the
  // same term to both and cancel, so this measures DISCRIMINATING BASES, not overhang length: a
  // read extending well past a breakpoint into sequence identical to the reference continuation
  // contributes nothing and correctly scores ~0.
  //
  // There is no `base` and no error-rate constant. Each discriminating base is worth
  // log10(3(1-eps)/eps) for its own quality -- ~4.5 at Q40, ~0.35 at Q3, ~0 at Q0 -- rather than a
  // flat 4 score units under an assumed 1% error rate. The earlier score-based form had to convert
  // score units to odds through exactly that assumption; with real qualities the ratio IS the odds.
  //
  // prior_odds_against is the only remaining input that is not read off the data, and callers set
  // it to (1-f)/f from the junction's own frequency (see reweight_window(),
  // resolve_alignments.cpp), which makes it a fitted quantity rather than a tuned one. The
  // half-weight point then sits at log10_odds_for_junction = log10(prior_odds_against).
  double junction_read_weight(double log10_odds_for_junction,
                              uint32_t n_ref_placements,
                              double prior_odds_against);

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