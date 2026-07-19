/*****************************************************************************

 AUTHORS

   Jeffrey E. Barrick <jeffrey.e.barrick@gmail.com> and other contributors

 LICENSE AND COPYRIGHT

   Copyright (c) 2008-2010 Michigan State University
   Copyright (c) 2011-2025 The University of Texas at Austin
   Copyright (c) 2025-     Michigan State University

   breseq is free software; you can redistribute it and/or modify it under the
   terms of the GNU General Public License as published by the Free Software
   Foundation; either version 2, or (at your option) any later version.

   SPDX-License-Identifier: GPL-2.0-or-later

*****************************************************************************/

// @DTF: VERSION: BINARY EDGE SEARCH 9

#include "coverage_distribution.h"

#include "genome_diff.h"
#include "stats.h"

using namespace std;

namespace breseq {

// Externally-visible (declared in coverage_distribution.h) so resolve_alignments/output can build the
// empirical coverage null for the JC pos_hash p-value. Kept out of the anonymous namespace below.
coverage_histogram_t read_coverage_histogram(const string& distribution_file_name)
{
  coverage_histogram_t hist;
  hist.max_coverage = 0;

  ifstream in(distribution_file_name.c_str());
  ASSERT(in.is_open(), "Could not open coverage distribution file: " + distribution_file_name);
  string header;
  getline(in, header);

  vector<pair<uint32_t, double> > rows;
  uint32_t coverage;
  double count;
  while (in >> coverage >> count) {
    rows.push_back(make_pair(coverage, count));
  }

  if (rows.empty()) {
    hist.n.assign(1, 0.0);
    return hist;
  }

  for (size_t i = 0; i < rows.size(); i++) {
    hist.max_coverage = max(hist.max_coverage, rows[i].first);
  }
  hist.n.assign(hist.max_coverage + 1, 0.0);
  for (size_t i = 0; i < rows.size(); i++) {
    hist.n[rows[i].first] = rows[i].second;
  }

  return hist;
}

bool use_empirical_pos_hash_coverage(const coverage_histogram_t& hist, double average,
                                     double neg_log10_p_value_cutoff,
                                     double& out_N, int32_t& out_deletion_floor)
{
  // Treat near-zero coverage as deletions (not valid junction sites). 0.05*average sits well below any
  // real low-coverage shoulder, so it removes only the deletion spike; it must NOT be the nbinom-derived
  // deletion_coverage_propagation_cutoff, which can land mid-shoulder precisely because that fit's lower
  // tail is too thin (the very flaw the empirical null fixes).
  int32_t deletion_floor = max(1, static_cast<int32_t>(round(0.05 * average)));
  out_deletion_floor = deletion_floor;

  double N = 0.0;
  for (int32_t c = deletion_floor + 1; c < static_cast<int32_t>(hist.n.size()); c++)
    N += hist.n[c];
  out_N = N;

  // Empirical resolution ceiling is ~log10(N); require a fixed 1-decade margin above the reject cutoff.
  return (N > 0.0) && (log(N) / log(10.0) >= neg_log10_p_value_cutoff + 1.0);
}

namespace {

// Result of fitting a negative binomial distribution to a censored count histogram via
// fit_censored_negative_binomial. Distribution-agnostic -- reused by both coverage and
// paired-mapping-distance fitting.
struct CensoredNegativeBinomialFitResult {
  double nb_fit_size = 0;
  double nb_fit_mu = 0;
  double average = 0;
  double variance = 0;
  double relative_variance = 0;
  uint32_t censor_start = 0;
  uint32_t censor_end = 0;
  double nb_fit_scale = 0;   // scales a raw dnbinom(...) curve to match the histogram's counts
};

/*! fit_censored_negative_binomial
 @abstract Fits a negative binomial (mu, size) to histogram n (n[i] = count of observations
 with value i, i in [0,N]) via censored least-squares: finds the peak of a 5-point moving
 average, censors the histogram to a window [0.5*peak, 1.5*peak], coarse-grains to ~1000-2000
 bins if larger, and optimizes a log-space-parametrized sum-of-squared-proportions objective
 via multi-restart Nelder-Mead. This is the statistical computation that used to be performed
 by coverage_distribution.r (which is now used only to draw the diagnostic plot, using the
 parameters this function returns -- see CoverageDistribution::fit), extracted here so it can
 be reused for other count histograms (e.g. paired-mapping distance distributions).
 @param n The histogram (n[i] = count of observations with value i); n[0] is ignored.
 @param N The histogram's highest index (n must have at least N+1 elements).
 !*/
CensoredNegativeBinomialFitResult fit_censored_negative_binomial(const vector<double>& n, uint32_t N)
{
  CensoredNegativeBinomialFitResult result;

  double total_positions = 0;
  for (uint32_t i = 1; i <= N; i++) total_positions += n[i];
  if (total_positions == 0) return result;

  double m = 0;
  for (uint32_t i = 1; i <= N; i++) m += i * n[i];
  m /= total_positions;

  // Guard the sample variance: with a single covered position total_positions == 1,
  // so v /= (total_positions - 1) would be 0/0 = NaN (and relative_variance NaN/Inf).
  // A NaN here is later serialized as JSON null and crashes on read-back, so keep
  // these finite (0) for degenerate histograms.
  double v = 0;
  if (total_positions > 1) {
    for (uint32_t i = 1; i <= N; i++) v += n[i] * (i - m) * (i - m);
    v /= (total_positions - 1);
  }

  result.average = m;
  result.variance = v;
  result.relative_variance = (m > 0) ? v / m : 0.0;

  // 5-point centered moving average, to more reliably find the peak.
  // Matches R's filter(X$n, c(1,1,1,1,1)/5): undefined within 2 positions of
  // either end (or everywhere, if there are too few points to smooth at all).
  vector<double> ma(N + 1, 0.0);
  vector<bool> ma_valid(N + 1, false);
  if (N >= 5) {
    for (uint32_t i = 3; i <= N - 2; i++) {
      ma[i] = (n[i - 2] + n[i - 1] + n[i] + n[i + 1] + n[i + 2]) / 5.0;
      ma_valid[i] = true;
    }
  } else {
    for (uint32_t i = 1; i <= N; i++) { ma[i] = n[i]; ma_valid[i] = true; }
  }

  uint32_t min_i = max(static_cast<uint32_t>(m / 4.0), static_cast<uint32_t>(1));
  double max_n = 0;
  uint32_t max_i = 0;
  for (uint32_t i = min_i; i <= N; i++) {
    if (ma_valid[i] && (ma[i] > max_n)) {
      max_n = ma[i];
      max_i = i;
    }
  }

  // Censor data on the right and left of the maximum.
  uint32_t start_i = max(static_cast<uint32_t>(floor(max_i * 0.5)), static_cast<uint32_t>(1));
  uint32_t end_i = min(static_cast<uint32_t>(ceil(max_i * 1.5)), N);
  if (start_i == end_i) return result;

  result.censor_start = start_i;
  result.censor_end = end_i;

  // Coarse-grain so that we are only fitting a number of bins that is
  // 1000-2000. The means of the negative binomial fit are adjusted
  // afterwards by multiplying by num_per_bin (the size parameter doesn't
  // need adjustment).
  uint32_t num_per_bin = (end_i - start_i) / 1000;
  vector<double> x_for_fits;
  uint32_t start_i_for_fits, end_i_for_fits;

  if (num_per_bin > 1) {
    start_i_for_fits = start_i / num_per_bin;
    end_i_for_fits = static_cast<uint32_t>(ceil(static_cast<double>(end_i) / num_per_bin));
    x_for_fits.assign(end_i_for_fits + 1, 0.0);
    for (uint32_t i = start_i_for_fits; i <= end_i_for_fits; i++) {
      for (uint32_t j = 1; j <= num_per_bin; j++) {
        uint32_t idx = i * num_per_bin + j;
        if (idx <= N) x_for_fits[i] += n[idx];
      }
    }
  } else {
    num_per_bin = 1;
    start_i_for_fits = start_i;
    end_i_for_fits = end_i;
    x_for_fits.assign(end_i_for_fits + 1, 0.0);
    for (uint32_t i = 1; i <= end_i_for_fits; i++) x_for_fits[i] = n[i];
  }

  double inner_total = 0;
  for (uint32_t i = start_i_for_fits; i <= end_i_for_fits; i++) inner_total += x_for_fits[i];

  double mean_estimate_num = 0, mean_estimate_den = 0;
  for (uint32_t i = 1; i <= end_i_for_fits; i++) {
    mean_estimate_num += i * x_for_fits[i];
    mean_estimate_den += x_for_fits[i];
  }
  double mean_estimate = mean_estimate_num / mean_estimate_den;

  // Objective: normalized sum-of-squares between the censored histogram and
  // a candidate negative binomial PMF, both expressed as proportions of
  // their respective totals over the fitting window.
  //
  // Parametrized in log-space (par = (log_mu, log_size), exponentiated
  // below) rather than R's f_nb(), which optimizes (mu, size) directly and
  // returns 0 -- not a penalty -- when mu<=0 or size<=0. That makes the
  // invalid region look like a perfect fit to R's nlm(), which turns out to
  // be load-bearing for R's behavior: nlm() readily wanders into that flat
  // "perfect" region, and the only thing rejecting those degenerate results
  // is the caller's explicit negative-size/code!=1 check, which forces
  // another restart. A literal port of that (replacing nlm with Nelder-Mead)
  // instead gets trapped in spurious near-zero-size local minima that are
  // positive (so pass the same check) but practically meaningless --
  // confirmed against R on tests/tmv_plasmid_circular_deletion_start_only's
  // distribution, where this previously produced mu=0/size=0 (fit "failed")
  // against R's mu=826.6/size=19.08. Optimizing in log-space removes the
  // invalid region entirely (mu,size are always positive by construction),
  // so there's no flat "perfect" boundary to get stuck against, and nearly
  // every restart below converges directly to the same true optimum.
  auto f_nb = [&](const vector<double>& par) -> double {
    double mu = exp(par[0]);
    double size = exp(par[1]);
    if (!isfinite(mu) || !isfinite(size)) return 1e10;

    vector<double> dist(end_i_for_fits + 1, 0.0);
    double total = 0;
    for (uint32_t i = start_i_for_fits; i <= end_i_for_fits; i++) {
      dist[i] = dnbinom_mu(static_cast<double>(i), size, mu);
      total += dist[i];
    }
    if (!(total > 0) || !isfinite(total)) return 1e10;

    double l = 0;
    for (uint32_t i = start_i_for_fits; i <= end_i_for_fits; i++) {
      double diff = (x_for_fits[i] / inner_total) - (dist[i] / total);
      l += diff * diff;
    }
    return isfinite(l) ? l : 1e10;
  };

  // Multi-restart strategy: try each of 6 candidate starting means (same
  // candidates coverage_distribution.r used) crossed with a geometrically
  // decreasing sequence of starting sizes (100000 down to 0.001), since the
  // fit is sensitive to its starting point. Unlike R's loop, which stops at
  // the first restart whose result clears its validity check, this scans
  // every restart and keeps the one with the lowest objective value among
  // those that converged -- with the log-space parametrization above nearly
  // all restarts converge to the same true optimum, but a few land in a
  // worse local minimum (also confirmed against this test's distribution),
  // and there's no reason to prefer "first" over "best" once the boundary
  // pathology that motivated R's early-exit no longer applies.
  vector<double> try_means;
  try_means.push_back(mean_estimate);
  try_means.push_back(static_cast<double>(end_i_for_fits));
  try_means.push_back(static_cast<double>(start_i_for_fits));
  try_means.push_back(1.0 * (end_i_for_fits + start_i_for_fits) / 4.0);
  try_means.push_back(2.0 * (end_i_for_fits + start_i_for_fits) / 4.0);
  try_means.push_back(3.0 * (end_i_for_fits + start_i_for_fits) / 4.0);

  double nb_fit_mu = 0;
  double nb_fit_size = 0;
  double best_objective = HUGE_VAL;

  for (size_t try_means_index = 0; try_means_index < try_means.size(); try_means_index++) {
    double try_mean = try_means[try_means_index];
    double try_size = 100000;
    for (int k = 0; k < 8; k++) {
      try_size /= 10.0;

      vector<double> initial_guess;
      initial_guess.push_back(log(try_mean));
      initial_guess.push_back(log(try_size));
      nelder_mead_result_t fit_result = nelder_mead_minimize(f_nb, initial_guess);
      if (!fit_result.converged) continue;

      double objective = f_nb(fit_result.estimate);
      if (objective < best_objective) {
        best_objective = objective;
        nb_fit_mu = exp(fit_result.estimate[0]);
        nb_fit_size = exp(fit_result.estimate[1]);
      }
    }
  }

  // Fit failed -- reset parameters so the caller (and the plotting script)
  // can recognize this.
  if (!(best_objective < HUGE_VAL)) {
    nb_fit_mu = 0;
    nb_fit_size = 0;
  }

  // Validate the fit by checking what fraction of it actually fell inside
  // the fitting window; things can go wrong and still produce a nb_fit_mu
  // that looks superficially plausible.
  double included_fract = 0;
  if (nb_fit_mu > 0) {
    double p = nb_fit_size / (nb_fit_size + nb_fit_mu);
    double end_fract = nbdtr(static_cast<double>(end_i_for_fits), nb_fit_size, p);
    double start_fract = nbdtr(static_cast<double>(start_i_for_fits), nb_fit_size, p);
    included_fract = end_fract - start_fract;

    if (included_fract >= 0.01) {
      // Adjust so that we are back in full coords before reporting the fit.
      if (num_per_bin > 1) nb_fit_mu = nb_fit_mu * num_per_bin;
    }
  }
  if (included_fract < 0.01) {
    nb_fit_mu = 0;
    nb_fit_size = 0;
  }

  result.nb_fit_mu = nb_fit_mu;
  result.nb_fit_size = nb_fit_size;
  if (nb_fit_mu > 0) result.nb_fit_scale = inner_total / included_fract;

  return result;
}

/*! fit_coverage_distribution
 @abstract Fits a negative binomial distribution to a unique-only coverage
 histogram (via fit_censored_negative_binomial) and derives a
 deletion-propagation coverage threshold. Ports the statistical computation
 that used to be performed by coverage_distribution.r (which is now used
 only to draw the diagnostic plot, using the parameters this function
 returns -- see CoverageDistribution::fit).
 @param hist The coverage histogram.
 @param deletion_propagation_pr_cutoff Probability cutoff for the returned
 deletion_coverage_propagation_cutoff.
 @param censor_start_out,censor_end_out The fitting window [start,end] used
 to censor the histogram around its peak, passed through to the plotting
 script so it doesn't have to recompute the same peak-finding logic.
 @param nb_fit_scale_out The factor (inner_total/included_fract) that scales
 a raw dnbinom(...) curve to match the histogram's actual counts, passed
 through to the plotting script for the same reason.
 !*/
CoverageDistributionFitResult fit_coverage_distribution(
                                                          const coverage_histogram_t& hist,
                                                          double deletion_propagation_pr_cutoff,
                                                          uint32_t& censor_start_out,
                                                          uint32_t& censor_end_out,
                                                          double& nb_fit_scale_out
                                                          )
{
  CoverageDistributionFitResult result;

  CensoredNegativeBinomialFitResult fit = fit_censored_negative_binomial(hist.n, hist.max_coverage);

  result.average = fit.average;
  result.variance = fit.variance;
  result.relative_variance = fit.relative_variance;
  result.nb_fit_mu = fit.nb_fit_mu;
  result.nb_fit_size = fit.nb_fit_size;
  censor_start_out = fit.censor_start;
  censor_end_out = fit.censor_end;
  nb_fit_scale_out = fit.nb_fit_scale;

  double m = fit.average;
  double v = fit.variance;
  double nb_fit_mu = fit.nb_fit_mu;
  double nb_fit_size = fit.nb_fit_size;

  // Fit the marginal value used for propagating deletions.
  double deletion_propagation_coverage;
  if (nb_fit_mu > 0) {
    deletion_propagation_coverage = static_cast<double>(qnbinom_mu(deletion_propagation_pr_cutoff, nb_fit_size, nb_fit_mu));
  } else {
    // Fallback to calculating off an estimate of just variance = mu + mu^2/size.
    double size_estimate = (1.0 / (v - m)) * m * m;
    if ((size_estimate > 0) && isfinite(size_estimate)) {
      deletion_propagation_coverage = static_cast<double>(qnbinom_mu(deletion_propagation_pr_cutoff, size_estimate, m));
    } else {
      deletion_propagation_coverage = -1;
    }
    if (!(deletion_propagation_coverage >= 1)) {
      // Double fallback to calculating as just 10% of the mean.
      deletion_propagation_coverage = m * 0.1;
    }
  }

  // Don't allow one read to indicate non-deleted regions.
  if (deletion_propagation_coverage < 1) deletion_propagation_coverage = 1;

  // If we have both low fit coverage and low straight average coverage,
  // then the reference sequence itself is deleted.
  if ((nb_fit_mu <= 3) && (m <= 3)) deletion_propagation_coverage = -1;

  result.deletion_coverage_propagation_cutoff = deletion_propagation_coverage;

  return result;
}

} // anonymous namespace

/*! fit
 @abstract Fits a negative binomial to the unique-only coverage histogram
 natively (see fit_coverage_distribution above), then draws the diagnostic
 plot from the resulting parameters using gnuplot.
 @param distribution_file_name Input histogram, created by error_count() and saved as *.unique_only_coverage_distribution.tab
 @param plot_file Output plot file (PDF).
 @return CoverageDistributionFitResult The fitted parameters and summary statistics.
 !*/

CoverageDistributionFitResult CoverageDistribution::fit(
                                         string distribution_file_name,
                                         string plot_file,
                                         double deletion_propagation_pr_cutoff
                                         )
{
  coverage_histogram_t hist = read_coverage_histogram(distribution_file_name);

  uint32_t censor_start = 0, censor_end = 0;
  double nb_fit_scale = 0;
  CoverageDistributionFitResult result = fit_coverage_distribution(hist, deletion_propagation_pr_cutoff, censor_start, censor_end, nb_fit_scale);

  // Nothing to plot if there's no data or no fit window was found (degenerate distribution).
  if (censor_start >= censor_end) return result;

  const vector<double>& n = hist.n;
  const uint32_t N = hist.max_coverage;

  // Don't graph very high values with very little coverage.
  uint32_t max_i = 1;
  double max_n = 0;
  for (uint32_t i = 1; i <= N; i++) {
    if (n[i] > max_n) { max_n = n[i]; max_i = i; }
  }
  uint32_t graph_end_i = max_i;
  while ((graph_end_i <= N) && (n[graph_end_i] > 0.01 * max_n)) graph_end_i++;
  // Leaves enough room to the right of the peak for the legend.
  graph_end_i = max(static_cast<uint32_t>(floor(2.2 * max_i)), graph_end_i);

  double max_y = max_n;
  for (uint32_t i = 1; i <= N; i++) max_y = max(max_y, n[i]);

  // The negative binomial fit curve has no on-disk representation -- write
  // it to its own small table so gnuplot can plot it like any other series.
  string nb_fit_table_file_name = distribution_file_name + ".nbfit.tab";
  if (result.nb_fit_mu > 0) {
    vector<double> fit_nb(N + 1, 0.0);
    for (uint32_t i = 0; i <= N; i++) {
      fit_nb[i] = dnbinom_mu(static_cast<double>(i), result.nb_fit_size, result.nb_fit_mu) * nb_fit_scale;
      max_y = max(max_y, fit_nb[i]);
    }
    ofstream nb_out(nb_fit_table_file_name.c_str());
    ASSERT(nb_out.is_open(), "Could not write to file: " + nb_fit_table_file_name);
    for (uint32_t i = 0; i <= N; i++) nb_out << i << "\t" << fit_nb[i] << endl;
    nb_out.close();
  }

  ostringstream s;
  s << "set terminal svg size 1320,720 font ',16'" << endl;
  s << "set output " << double_quote(plot_file) << endl;
  s << "set tics out" << endl;
  s << "set border lw 2" << endl;
  s << "set title 'Coverage Distribution at Unique-Only Positions' font ',20'" << endl;
  s << "set xlabel 'Coverage depth (reads)'" << endl;
  s << "set ylabel 'Number of reference positions'" << endl;
  s << "set xrange [0:" << graph_end_i << "]" << endl;
  s << "set yrange [0:" << to_string(max_y * 1.05, 6) << "]" << endl;
  s << "set key top right font ',16' spacing 2" << endl;

  vector<string> plot_clauses;
  plot_clauses.push_back(double_quote(distribution_file_name) + " using ($1>=" + to_string(censor_start) + "&&$1<=" + to_string(censor_end) + "?$1:NaN):2 with points pt 6 lc rgb 'black' title 'Coverage distribution'");
  plot_clauses.push_back(double_quote(distribution_file_name) + " using ($1<" + to_string(censor_start) + "||$1>" + to_string(censor_end) + "?$1:NaN):2 with points pt 6 lc rgb 'red' title 'Censored data'");
  if (result.nb_fit_mu > 0) {
    plot_clauses.push_back(double_quote(nb_fit_table_file_name) + " with lines lw 3 lc rgb 'black' title 'Negative binomial'");
  }
  s << "plot " << join(plot_clauses, string(", \\\n     ")) << endl;

  string script_base_name = plot_file + "." + to_string(getpid());
  string gnuplot_script_name = script_base_name + ".gp";
  string log_file_name = script_base_name + ".gp.log";
  run_gnuplot_script(s.str(), gnuplot_script_name, log_file_name);
  make_svg_responsive(plot_file);
  remove(log_file_name.c_str());
  if (result.nb_fit_mu > 0) remove(nb_fit_table_file_name.c_str());

  return result;
}

// helper functions
/*! analyze_unique_coverage_distribution
 @abstract Fits and plots the unique-only coverage distribution for one coverage group.


 @param settings
 @param summary
 @param ref_seq_info
 !*/

void CoverageDistribution::analyze_unique_coverage_distribution(
                                                                Settings& settings,
                                                                Summary& summary,
                                                                cReferenceSequences& ref_seq_info,
                                                                uint32_t coverage_group_id,
                                                                string plot_key,
                                                                string distribution_file_name,
                                                                string step_key
                                                                )
{
  
  vector<string> seq_ids = settings.refseq_settings.m_seq_ids_by_coverage_group[coverage_group_id];
  
  // We need to know the total sequence length of just this group
  uint32_t sequence_length = 0;
  for (vector<string>::iterator it=seq_ids.begin(); it!=seq_ids.end(); it++) {
    string& seq_id = *it;
    sequence_length += ref_seq_info.get_sequence_length(seq_id);
  }
  
  // Perform no fitting if we are in targeted_sequencing mode.
  if (settings.targeted_sequencing) return;
  
  string unique_only_coverage_plot_file_name = settings.file_name(plot_key, "@", to_string<uint32_t>(coverage_group_id));
  string unique_only_coverage_distribution_file_name = settings.file_name(distribution_file_name, "@", to_string<uint32_t>(coverage_group_id));
  
  // Define various coverage thresholds...
  
  /// DELETION PROPAGATION CUTOFF
  // One-tailed test p=0.05, Bonferroni correction
  //# my del_propagation_pr_cutoff = 0.05 / sequence_length;
  
  // One-tailed test p=0.01, no Bonferroni correction
  //#my del_propagation_pr_cutoff = 0.01;
  
  // We really want somewhere between these two, try this...
  double deletion_propagation_pr_cutoff = 0.05 / sqrt(sequence_length);
  
  /// NEW JUNCTION COVERAGE CUTOFFS
  // Arbitrary value that seems to work....
  //double junction_coverage_pr_cutoff =  sqrt(settings.junction_accept_pr / static_cast<double>(sequence_length));
  double junction_coverage_pr_cutoff = 0.01;
  
  // We really want somewhere between these two, try this...
  double junction_accept_pr_cutoff = 0.01;
  double junction_keep_pr_cutoff = 0.01 / sqrt(sequence_length);
  int32_t junction_max_score = int(2 * summary.sequence_conversion.read_length_avg);
  
  CoverageDistribution dist;
  CoverageDistributionFitResult fit_result = dist.fit(
                                  unique_only_coverage_distribution_file_name,
                                  unique_only_coverage_plot_file_name,
                                  deletion_propagation_pr_cutoff
                                  );
  settings.track_intermediate_file(step_key, unique_only_coverage_plot_file_name);
  settings.track_intermediate_file(step_key, unique_only_coverage_distribution_file_name);

  // Put these into summary

  for (vector<string>::iterator it=seq_ids.begin(); it!=seq_ids.end(); it++) {
    string seq_id = *it;
    summary.unique_coverage[seq_id].nbinom_size_parameter = fit_result.nb_fit_size;
    summary.unique_coverage[seq_id].nbinom_mean_parameter = fit_result.nb_fit_mu;
    // Calculated by formula, prob = size/(size + mu)

    // These remain at their defaults of zero if this is zero = fit failed
    if (summary.unique_coverage[seq_id].nbinom_mean_parameter != 0) {
      summary.unique_coverage[seq_id].nbinom_prob_parameter = summary.unique_coverage[seq_id].nbinom_size_parameter
      / (summary.unique_coverage[seq_id].nbinom_mean_parameter + summary.unique_coverage[seq_id].nbinom_size_parameter);
      // Calculated by formula variance = mu + mu ^ 2 / size
      summary.unique_coverage[seq_id].nbinom_variance = summary.unique_coverage[seq_id].nbinom_mean_parameter
      + pow(summary.unique_coverage[seq_id].nbinom_mean_parameter, 2) / summary.unique_coverage[seq_id].nbinom_size_parameter;
      // Calculated by formula relative_variance = variance / mu
      summary.unique_coverage[seq_id].nbinom_relative_variance = summary.unique_coverage[seq_id].nbinom_variance / summary.unique_coverage[seq_id].nbinom_mean_parameter;
    }

    summary.unique_coverage[seq_id].average = fit_result.average;
    summary.unique_coverage[seq_id].variance = fit_result.variance;
    summary.unique_coverage[seq_id].relative_variance = fit_result.relative_variance;

    summary.unique_coverage[seq_id].deletion_coverage_propagation_cutoff = fit_result.deletion_coverage_propagation_cutoff;

    bool verbose = false;
    if (verbose)
    {
      cout << seq_id << endl;
      cout << "nbinom_size_parameter " << summary.unique_coverage[seq_id].nbinom_size_parameter << endl;
      cout << "nbinom_mean_parameter " << summary.unique_coverage[seq_id].nbinom_mean_parameter << endl;
      cout << "nbinom_prob_parameter " << summary.unique_coverage[seq_id].nbinom_prob_parameter << endl;
      cout << "average " << summary.unique_coverage[seq_id].average << endl;
      cout << "variance " << summary.unique_coverage[seq_id].variance << endl;
      cout << "relative_variance " << summary.unique_coverage[seq_id].relative_variance << endl;
      cout << "deletion_coverage_propagation_cutoff " << summary.unique_coverage[seq_id].deletion_coverage_propagation_cutoff << endl;
    }
  }
}

void CoverageDistribution::analyze_unique_coverage_distributions(
                                                                 Settings& settings,
                                                                 Summary& summary,
                                                                 cReferenceSequences& ref_seq_info,
                                                                 string plot_file_name,
                                                                 string distribution_file_name,
                                                                 string step_key
                                                                 )
{
  vector<vector<string> > seq_ids_by_coverage_group = settings.seq_ids_by_coverage_group();

  for (uint32_t i=0; i< settings.seq_ids_by_coverage_group().size(); i++) {


    analyze_unique_coverage_distribution(
                                         settings,
                                         summary,
                                         ref_seq_info,
                                         i,
                                         plot_file_name,
                                         distribution_file_name,
                                         step_key
                                         );

    //Warning
  }
}

// Smallest index i such that the cumulative histogram count through i reaches half of total --
// a standard "lower median" for grouped/discrete count data.
static uint32_t weighted_median(const vector<double>& n, uint32_t N, double total)
{
  double half = total / 2.0;
  double cumulative = 0;
  for (uint32_t i = 0; i <= N; i++) {
    cumulative += n[i];
    if (cumulative >= half) return i;
  }
  return N;
}

/*! fit
 @abstract Reads the condensed (orientation,distance,count) CSV written by
 PreprocessAlignments::merge_two_sets_of_paired_sam_files/preprocess_one_set_of_paired_sam_files,
 determines the majority orientation by total observation count, computes the median, MAD, and
 a modified-z-score outlier cutoff (Iglewicz & Hoaglin) for the majority-orientation distance
 histogram (ignoring rows for other orientations), then draws a diagnostic plot.
 @param distribution_file_name Input CSV.
 @param plot_file Output plot file (SVG).
 @return PairedMappingDistanceDistributionFitResult The computed statistics.
 !*/
PairedMappingDistanceDistributionFitResult PairedMappingDistanceDistribution::fit(
                                         string distribution_file_name,
                                         string plot_file
                                         )
{
  PairedMappingDistanceDistributionFitResult result;

  ifstream in(distribution_file_name.c_str());
  ASSERT(in.is_open(), "Could not open paired mapping distance distribution file: " + distribution_file_name);
  string header_line;
  getline(in, header_line);

  vector<string> row_orientation;
  vector<int64_t> row_distance;
  vector<uint32_t> row_count;
  map<string, uint64_t> orientation_totals;
  uint64_t total_all_orientations = 0;

  string line;
  while (getline(in, line)) {
    if (line.empty()) continue;
    istringstream iss(line);
    string orientation, distance_str, count_str;
    getline(iss, orientation, ',');
    getline(iss, distance_str, ',');
    getline(iss, count_str, ',');
    int64_t distance = from_string<int64_t>(distance_str);
    uint32_t count = from_string<uint32_t>(count_str);

    row_orientation.push_back(orientation);
    row_distance.push_back(distance);
    row_count.push_back(count);
    orientation_totals[orientation] += count;
    total_all_orientations += count;
  }

  if (orientation_totals.empty()) return result;

  // Majority rule: the orientation with the largest total observation count.
  string majority_orientation;
  uint64_t best_total = 0;
  for (map<string, uint64_t>::iterator it = orientation_totals.begin(); it != orientation_totals.end(); it++) {
    if (it->second > best_total) {
      best_total = it->second;
      majority_orientation = it->first;
    }
  }
  result.majority_orientation = majority_orientation;

  // Build a dense distance->count histogram from only the majority-orientation rows.
  int64_t max_distance = 0;
  for (size_t i = 0; i < row_orientation.size(); i++) {
    if (row_orientation[i] == majority_orientation) max_distance = max(max_distance, row_distance[i]);
  }
  if (max_distance <= 0) return result;

  uint32_t N = static_cast<uint32_t>(max_distance);
  vector<double> n(N + 1, 0.0);
  for (size_t i = 0; i < row_orientation.size(); i++) {
    if (row_orientation[i] == majority_orientation) n[static_cast<uint32_t>(row_distance[i])] += row_count[i];
  }

  double total = 0;
  for (uint32_t i = 0; i <= N; i++) total += n[i];
  if (total == 0) return result;

  uint32_t median_i = weighted_median(n, N, total);
  result.median = median_i;

  // MAD (upper, one-sided): the weighted median of (i - median_i) over ONLY the observations
  // ABOVE the median. The distance distribution is right-skewed (a tight lower bound near the
  // read length, a long tail of larger fragments), so a symmetric MAD is dominated by the narrow
  // lower side and yields an artificially tight cutoff that flags many normal right-tail fragments
  // as discordant. Using only the upper deviations makes the spread estimate reflect the tail that
  // the cutoff actually needs to accommodate.
  uint32_t max_dev = N - median_i;
  vector<double> dev_n(max_dev + 1, 0.0);
  double upper_total = 0;
  for (uint32_t i = median_i + 1; i <= N; i++) {
    if (n[i] == 0) continue;
    dev_n[i - median_i] += n[i];
    upper_total += n[i];
  }
  uint32_t mad_i = (upper_total > 0) ? weighted_median(dev_n, max_dev, upper_total) : 0;
  result.mad = mad_i;

  // Modified z-score rule (0.6745*(x-median)/MAD, per Iglewicz & Hoaglin), but with a
  // Bonferroni-corrected threshold instead of the conventional fixed 3.5: target a family-wise
  // false-positive rate of 0.01 across all total_all_orientations tested pairs, giving
  // alpha = 0.01/N per pair, and take the one-sided standard-normal quantile at that alpha as
  // the z threshold. Flag as an outlier when 0.6745*(x-median)/MAD > z*  =>  x > median + z*MAD/0.6745.
  double alpha = 0.01 / static_cast<double>(total_all_orientations);
  double z_multiplier = ndtri(1.0 - alpha);
  result.distance_cutoff = (result.mad > 0) ? (result.median + z_multiplier * result.mad / 0.6745) : result.median;

  // Preliminary (preprocessing) pair counts: all same-reference mapped pairs, and the concordant
  // subset (majority orientation with distance within the cutoff).
  result.mapped_pairs = static_cast<double>(total_all_orientations);
  double concordant = 0;
  uint32_t cutoff_i = static_cast<uint32_t>(min<double>(result.distance_cutoff, static_cast<double>(N)));
  for (uint32_t i = 0; i <= cutoff_i; i++) concordant += n[i];
  result.concordant_pairs = concordant;

  // Display window: [0, 2*cutoff], falling back to the full histogram range if the cutoff
  // collapsed to 0 (e.g. a degenerate all-identical-value histogram with MAD == 0).
  uint32_t graph_end_i = (result.distance_cutoff > 0) ? static_cast<uint32_t>(ceil(2 * result.distance_cutoff)) : N;
  graph_end_i = min(graph_end_i, N);

  double max_y = 0;
  for (uint32_t i = 0; i <= graph_end_i; i++) max_y = max(max_y, n[i]);

  // The condensed CSV covers all orientations and isn't directly plot-ready (3 comma-separated
  // columns) -- write the majority-orientation histogram to its own small 2-column table so
  // gnuplot can plot it directly, exactly like the cutoff-line sidecar table below.
  string hist_table_file_name = distribution_file_name + ".majority_hist.tab";
  {
    ofstream hist_out(hist_table_file_name.c_str());
    ASSERT(hist_out.is_open(), "Could not write to file: " + hist_table_file_name);
    for (uint32_t i = 0; i <= N; i++) hist_out << i << "\t" << n[i] << endl;
    hist_out.close();
  }

  // The cutoff line has no on-disk representation -- write it to its own small table (a single
  // vertical segment at x=distance_cutoff spanning y=0 to y=max_y) so gnuplot can plot it like
  // any other series.
  string cutoff_table_file_name = distribution_file_name + ".cutoff.tab";
  {
    ofstream cutoff_out(cutoff_table_file_name.c_str());
    ASSERT(cutoff_out.is_open(), "Could not write to file: " + cutoff_table_file_name);
    cutoff_out << to_string(result.distance_cutoff, 6) << "\t0" << endl;
    cutoff_out << to_string(result.distance_cutoff, 6) << "\t" << to_string(max_y, 6) << endl;
    cutoff_out.close();
  }

  // Same sidecar-table convention for a vertical line at the median.
  string median_table_file_name = distribution_file_name + ".median.tab";
  {
    ofstream median_out(median_table_file_name.c_str());
    ASSERT(median_out.is_open(), "Could not write to file: " + median_table_file_name);
    median_out << to_string(result.median, 6) << "\t0" << endl;
    median_out << to_string(result.median, 6) << "\t" << to_string(max_y, 6) << endl;
    median_out.close();
  }

  ostringstream s;
  s << "set terminal svg size 1320,720 font ',16'" << endl;
  s << "set output " << double_quote(plot_file) << endl;
  s << "set tics out" << endl;
  s << "set border lw 2" << endl;
  s << "set title 'Paired-End Mapping Distance Distribution (" << majority_orientation << ")' font ',20'" << endl;
  s << "set xlabel 'Distance between outermost read coordinates (bp)'" << endl;
  s << "set ylabel 'Number of read pairs'" << endl;
  s << "set xrange [0:" << graph_end_i << "]" << endl;
  s << "set yrange [0:" << to_string(max_y * 1.05, 6) << "]" << endl;
  s << "set key top right font ',16' spacing 2" << endl;

  // Light grey band spanning [median-MAD, median+MAD], full plot height, drawn behind the data.
  s << "set object 1 rectangle from " << to_string(result.median - result.mad, 6) << ", graph 0 to "
    << to_string(result.median + result.mad, 6) << ", graph 1 fc rgb 'grey90' fillstyle solid 1.0 noborder behind" << endl;

  vector<string> plot_clauses;
  plot_clauses.push_back(double_quote(hist_table_file_name) + " using 1:2 with points pt 6 lc rgb 'black' title 'Distance distribution'");
  plot_clauses.push_back(double_quote(median_table_file_name) + " with lines lw 2 lc rgb 'black' title 'Median'");
  plot_clauses.push_back(double_quote(cutoff_table_file_name) + " with lines lw 2 lc rgb 'red' title 'Distance cutoff'");
  s << "plot " << join(plot_clauses, string(", \\\n     ")) << endl;

  string script_base_name = plot_file + "." + to_string(getpid());
  string gnuplot_script_name = script_base_name + ".gp";
  string log_file_name = script_base_name + ".gp.log";
  run_gnuplot_script(s.str(), gnuplot_script_name, log_file_name);
  make_svg_responsive(plot_file);
  remove(log_file_name.c_str());
  remove(hist_table_file_name.c_str());
  remove(cutoff_table_file_name.c_str());
  remove(median_table_file_name.c_str());

  return result;
}

/*! fit_paired_mapping_distance_distribution
 @abstract Fits and plots the paired-mapping distance distribution for one paired read file set.
 !*/
void PairedMappingDistanceDistribution::fit_paired_mapping_distance_distribution(
                                                                                 Settings& settings,
                                                                                 Summary& summary,
                                                                                 const cReadFileSet& read_file_set
                                                                                 )
{
  string distribution_file_name = Settings::file_name(settings.paired_mapping_distance_distribution_file_name, "#", read_file_set.m_base_name);
  string plot_file_name = Settings::file_name(settings.paired_mapping_distance_plot_file_name, "#", read_file_set.m_base_name);

  PairedMappingDistanceDistribution dist;
  PairedMappingDistanceDistributionFitResult fit_result = dist.fit(distribution_file_name, plot_file_name);

  // The pair_stats CSV is a transient intermediate: now that the plot is drawn and the stats are
  // extracted, mark it for deletion when this step completes (done_step honors keep-all-intermediates).
  settings.track_intermediate_file(settings.paired_mapping_distance_done_file_name, distribution_file_name);

  // Stage-03 fit is the PRELIMINARY record. The final paired_mapping_distance_distribution is filled
  // later from the BAM pair-flag assignment in resolve_alignments (which reads the cutoff from here).
  PairedMappingDistanceDistributionSummary& prelim = summary.preliminary_paired_mapping_distance_distribution[read_file_set.m_base_name];
  prelim.median = fit_result.median;
  prelim.mad = fit_result.mad;
  prelim.distance_cutoff = fit_result.distance_cutoff;
  prelim.majority_orientation = fit_result.majority_orientation;
  prelim.mapped_pairs = fit_result.mapped_pairs;
  prelim.concordant_pairs = fit_result.concordant_pairs;
}

void PairedMappingDistanceDistribution::fit_paired_mapping_distance_distributions(
                                                                                  Settings& settings,
                                                                                  Summary& summary
                                                                                  )
{
  for (const auto& rfs : settings.read_file_sets) {
    if (!rfs.is_paired()) continue;
    fit_paired_mapping_distance_distribution(settings, summary, rfs);
  }
}


/*
 * Function: calculate_periodicity
 * --------------------------------
 * SUMMARY GOES HERE
 */
void CoverageDistribution::analyze_coverage_bias (
                                                   const string& _fasta_file_name,
                                                   const string& _bam_file_name,
                                                   const string& _output_file_name_prefix,
                                                   const int32_t _read_length
                                                   )
  
{
  
  cReferenceSequences ref_seq_info;
  ref_seq_info.LoadFile(_fasta_file_name);
  
  bam_file bam;
  bam.open_read(_bam_file_name, _fasta_file_name);
  
  
  string read_output_file_name = _output_file_name_prefix + ".read.csv";
  ofstream read_out(read_output_file_name.c_str());
  
  // First line of output is read_length
  read_out << _read_length << endl;
  
  // Second line is lengths of reference sequences
  for (cReferenceSequences::iterator it=ref_seq_info.begin(); it != ref_seq_info.end(); it++) {
    if (it != ref_seq_info.begin()) read_out << ",";
    read_out << it->get_sequence_length() << endl;
  }
  
  // Third line is GC of reference sequences
  for (cReferenceSequences::iterator it=ref_seq_info.begin(); it != ref_seq_info.end(); it++) {
    if (it != ref_seq_info.begin()) read_out << ",";
    read_out << gc_percentage_string(it->m_fasta_sequence.get_sequence()) << endl;
  }
  
  alignment_list alignments;
  while (bam.read_alignments(alignments, false)) {
    double gc(0);

    for(alignment_list::iterator it=alignments.begin(); it!=alignments.end(); it++)
    {
      bam_alignment& a = *(it->get());
      
      int64_t start = (a.strand()==-1) ? a.reference_end_1() - _read_length + 1 : a.reference_start_1();
      string read_string = ref_seq_info[a.reference_target_id()].get_sequence_1_start_size(start, _read_length);
      // Don't need to reverse complement for GC content
      
      gc += gc_percentage_string(read_string);
    }

    gc /= alignments.size();
    
    read_out << gc << endl;
  }
  
  // Now the reference file
  string ref_output_file_name = _output_file_name_prefix + ".ref.csv";
  ofstream ref_out(ref_output_file_name.c_str());
  for (cReferenceSequences::iterator it=ref_seq_info.begin(); it != ref_seq_info.end(); it++) {

    for (size_t i=1; i<= it->get_sequence_length(); i++) {
      string read_string = it->get_sequence_1_start_size(i, _read_length);
      ref_out << gc_percentage_string(read_string) << endl;
    }
  }
}
  

} // namespace breseq
