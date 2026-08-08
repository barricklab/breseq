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

#include "pd_evidence.h"
#include "dp_evidence.h"

#include "genome_diff.h"
#include "pileup.h"
#include "pileup_base.h"
#include "alignment.h"
#include "stats.h"

#include <set>
#include <random>
#include <algorithm>

using namespace std;

namespace breseq {

  // One-sided confidence level for the PD local-frequency (Clopper-Pearson) bounds. Matches DP and MP.
  static const double kPDFrequencyAlpha = 0.05;

  // Likelihood-ratio deadband for sorting a covering pair into "supports the call" / "argues against
  // it" / "cannot tell". See pd_classify_pairs for why a hard threshold at the crossover is wrong.
  static const double kPDSupportOdds = 3.0;

  // Half-width of the chi-square(1) 95% profile-likelihood interval, in log-likelihood units.
  static const double kPDLogLDrop = 1.92;

  // ---------------------------------------------------------------------------------------------
  // The null distribution
  // ---------------------------------------------------------------------------------------------

  bool pd_load_weighted_histogram(const string& file_name, int32_t trunc,
                                  vector<double>& weighted, double& total)
  {
    weighted.clear();
    total = 0.0;
    if (!file_exists(file_name.c_str())) return false;

    ifstream in(file_name.c_str());
    string line;
    while (getline(in, line)) {
      if (line.empty()) continue;
      vector<string> fields = split(line, "\t");
      if (fields.size() < 2) continue;
      int32_t d = from_string<int32_t>(fields[0]);
      double n = from_string<double>(fields[1]);
      if ((d < 0) || (n <= 0.0)) continue;
      double weight = static_cast<double>(d - trunc);
      if (weight <= 0.0) continue;   // overlapping mates: no gap, so never a covering pair
      if (static_cast<int32_t>(weighted.size()) <= d) weighted.resize(d + 1, 0.0);
      weighted[d] = n * weight;
      total += weighted[d];
    }
    if (total <= 0.0) { weighted.clear(); total = 0.0; return false; }
    return true;
  }

  namespace {

  //! The length-bias-corrected pair distance null, smoothed and normalized once up front.
  //
  //  g(d) is the probability that a pair whose gap covers an arbitrary reference point has distance
  //  d, on a genome with no event there. Everything PD decides -- the size shift, the frequency, the
  //  snap -- is a likelihood ratio against this.
  struct pd_null {
    vector<double> g;     // normalized, smoothed density indexed by distance
    double         floor; // small non-zero density outside the observed range
    int32_t        trunc; // 2 * read_length
    int32_t        median;

    pd_null() : floor(0.0), trunc(0), median(0) {}
    bool ok() const { return !g.empty() && (trunc > 0) && (median > 0); }
    int32_t max_distance() const { return static_cast<int32_t>(g.size()) - 1; }

    double at(int32_t d) const {
      if ((d < 0) || (d >= static_cast<int32_t>(g.size()))) return floor;
      return g[d] > floor ? g[d] : floor;
    }

    //! Build from the shared weighted histogram. `smooth` is the half-width of the moving average,
    //  which keeps a sparsely-populated bin from dominating a likelihood ratio.
    bool build(const string& file_name, int32_t read_length, int32_t smooth)
    {
      trunc = 2 * read_length;
      vector<double> w;
      double total = 0.0;
      if (!pd_load_weighted_histogram(file_name, trunc, w, total)) return false;

      g.assign(w.size(), 0.0);
      double norm = 0.0;
      for (size_t d = 0; d < w.size(); d++) {
        double s = 0.0;
        for (int32_t x = static_cast<int32_t>(d) - smooth; x <= static_cast<int32_t>(d) + smooth; x++)
          if ((x >= 0) && (x < static_cast<int32_t>(w.size()))) s += w[x];
        g[d] = s / (2.0 * smooth + 1.0);
        norm += g[d];
      }
      if (norm <= 0.0) { g.clear(); return false; }
      for (size_t d = 0; d < g.size(); d++) g[d] /= norm;
      floor = 1.0 / (static_cast<double>(g.size()) * 1000.0);

      // Median of the CORRECTED distribution, which is what a supporting pair's distance is centred
      // on once the size shift is removed -- so it, not the raw library median, is the reference the
      // initial shift estimate is measured from.
      double run = 0.0;
      median = 0;
      for (size_t d = 0; d < g.size(); d++) {
        run += g[d];
        if (run >= 0.5) { median = static_cast<int32_t>(d); break; }
      }
      return ok();
    }

    //! Density of a SUPPORTING pair's distance given a size shift of delta.
    //
    //  A pair supports the call only if its gap accommodates the whole event, and that selection is
    //  not the same in the two directions -- which is exactly the correction that keeps the estimate
    //  from being dragged toward zero.
    //
    //  Write the pair's true fragment length as L = d - delta and the number of inserted bases as
    //  I = max(0, -delta). The fragment is seen as a supporting pair once for each position at which
    //  its sample-coordinate gap (L - 2r bases) can hold the event, which is (L - 2r - I). For a
    //  deletion I is zero and that is just the weight already folded into g, so h reduces to
    //  g(d - delta). For an insertion it is not: the longer the insertion, the fewer of the fragments
    //  of a given length can bracket it, and ignoring that systematically understates the insertion.
    //
    //  Returns an UNNORMALIZED density; the caller divides by total_supporting(delta).
    //
    //  Note this reads the raw density, NOT the floored at(). The floor exists so that one pair with
    //  an unobserved distance cannot veto an otherwise good fit, but it must not be visible here: a
    //  flat floor is something the shift can be slid onto, and since the normalizer shrinks as more
    //  of the range becomes inadmissible, a shift far outside the observed distribution would score
    //  arbitrarily well. Leaving it at zero makes such a shift produce no supporting density at all,
    //  which is the truth -- no fragment of that length exists in this library.
    double supporting(int32_t d, int32_t delta) const
    {
      int32_t L = d - delta;
      if ((L < 0) || (L >= static_cast<int32_t>(g.size()))) return 0.0;
      double gap = static_cast<double>(L - trunc);
      if (gap <= 0.0) return 0.0;
      double usable = gap - static_cast<double>(delta < 0 ? -delta : 0);
      if (usable <= 0.0) return 0.0;
      return g[L] * (usable / gap);
    }

    double total_supporting(int32_t delta) const
    {
      double s = 0.0;
      int32_t lo = max(0, trunc + delta);
      int32_t hi = max_distance() + delta;
      for (int32_t d = lo; d <= hi; d++) s += supporting(d, delta);
      return s;
    }
  };

  // ---------------------------------------------------------------------------------------------
  // Candidate regions
  // ---------------------------------------------------------------------------------------------

  //! One candidate region parsed from PD_candidate_regions.csv.
  struct pd_region_row {
    string   seq_id;
    uint32_t start;
    uint32_t end;
    char     direction;      // 'L' = pairs mapping farther apart, 'S' = closer
    uint32_t peak_position;
    int32_t  peak_covering;
    int32_t  peak_tail;
    double   peak_z;
  };

  // ---------------------------------------------------------------------------------------------
  // BAM rescan
  // ---------------------------------------------------------------------------------------------

  //! One read pair seen in the rescan, described entirely from its leftmost mate.
  struct pd_pair {
    int32_t read_start;   // leftmost mate's aligned start
    int32_t e1;           // leftmost mate's aligned end   -- low end of the gap, in b-coordinates
    int32_t s2;           // mate's aligned start          -- gap runs to s2 - 1
    int32_t d;            // reference distance (breseq's outer span)

    //! Does this pair's gap reach position b at all? This is the denominator condition: the molecule
    //  was sampled across b, so it has something to say about what is there.
    bool covers(int32_t b) const { return (e1 <= b) && (b <= s2 - 1); }

    //! Can this pair's gap hold an event that removes `deleted` reference bases starting after b?
    //  Stricter than covers(), and it is what a SUPPORTING pair must satisfy: if any part of the
    //  removed span fell under either mate, that mate would have been clipped at the breakpoint and
    //  the pair would not be here in this form.
    bool accommodates(int32_t b, int32_t deleted) const { return (e1 <= b) && (b + deleted <= s2 - 1); }
  };

  //! Collects every pair whose gap could bear on a candidate, over a window around it.
  //
  //  The qualification here is deliberately identical to the seeding step in identify_mutations.cpp,
  //  so a region's own pairs are exactly the ones re-derived here. In particular it does NOT require
  //  proper_pair: a large enough event pushes its spanning pairs past the discordant cutoff, and PD
  //  has to see those too, both to measure them and to be able to supersede the DP they produce.
  class pd_pair_scanner : public pileup_base {
  public:
    pd_pair_scanner(const string& bam, const string& fasta, int32_t max_span)
      : pileup_base(bam, fasta), m_max_span(max_span), m_out(NULL) { set_print_progress(false); }

    int32_t tid_for_seq_id(const string& seq_id) const {
      for (uint32_t t = 0; t < num_targets(); t++)
        if (seq_id == string(target_name(t))) return static_cast<int32_t>(t);
      return -1;
    }
    int32_t seq_length(int32_t tid) const { return static_cast<int32_t>(target_length(tid)); }

    void collect(const string& seq_id, int32_t lo, int32_t hi, vector<pd_pair>& out)
    {
      out.clear();
      int32_t tid = tid_for_seq_id(seq_id);
      if (tid < 0) return;
      lo = max(1, lo);
      hi = min(seq_length(tid), hi);
      if (lo > hi) return;
      m_out = &out;
      do_fetch(seq_id + ":" + to_string(lo) + "-" + to_string(hi));
      m_out = NULL;
    }

    void fetch_callback(const alignment_wrapper& a)
    {
      if (m_out == NULL) return;
      if (!a.is_paired() || a.unmapped()) return;
      if (a.flag() & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY | BAM_FMUNMAP)) return;
      if (a.mate_reference_target_id() != a.reference_target_id()) return;
      if (a.insert_size() <= 0) return;                                   // leftmost mate only
      if (a.reversed() == ((a.flag() & BAM_FMREVERSE) != 0)) return;      // not FR or RF
      if (a.redundancy() > 1) return;                                     // arbitrary placement

      pd_pair p;
      p.read_start = static_cast<int32_t>(a.reference_start_1());
      p.e1 = static_cast<int32_t>(a.reference_end_1());
      p.s2 = a.mate_start_1();
      p.d  = a.insert_size();
      if (p.s2 - 1 < p.e1) return;        // mates overlap: no gap, so no positional information
      if (p.d > m_max_span) return;
      m_out->push_back(p);
    }

  private:
    int32_t m_max_span;
    vector<pd_pair>* m_out;
  };

  // ---------------------------------------------------------------------------------------------
  // Estimating the size shift and the local frequency
  // ---------------------------------------------------------------------------------------------

  //! Log-likelihood of the covering pairs under "a fraction f of the molecules here carry an event
  //  of size shift delta", plus the f that maximizes it.
  //
  //  The mixture is  f * h_delta(d) + (1 - f) * g(d):  a molecule either carries the event, in which
  //  case its pair distance follows the shifted-and-reselected null, or it does not, in which case it
  //  follows the null itself. Both components are normalized, so f is a genuine frequency.
  static double pd_profile_logL(const vector<pd_pair>& pairs, const pd_null& null,
                                int32_t delta, double& f_hat)
  {
    double H = null.total_supporting(delta);
    if (H <= 0.0) { f_hat = 0.0; return -numeric_limits<double>::infinity(); }

    // A few EM rounds are plenty: the likelihood in f alone is concave, and delta -- the parameter
    // actually being profiled -- is held fixed here.
    double f = 0.5;
    for (int iter = 0; iter < 12; iter++) {
      double sum_post = 0.0;
      for (size_t i = 0; i < pairs.size(); i++) {
        double hs = null.supporting(pairs[i].d, delta) / H;
        double hn = null.at(pairs[i].d);
        double num = f * hs;
        double den = num + (1.0 - f) * hn;
        if (den > 0.0) sum_post += num / den;
      }
      double f_new = sum_post / static_cast<double>(pairs.size());
      if (f_new < 1e-9) f_new = 1e-9;
      if (f_new > 1.0 - 1e-9) f_new = 1.0 - 1e-9;
      if (fabs(f_new - f) < 1e-9) { f = f_new; break; }
      f = f_new;
    }

    double logL = 0.0;
    for (size_t i = 0; i < pairs.size(); i++) {
      double hs = null.supporting(pairs[i].d, delta) / H;
      double hn = null.at(pairs[i].d);
      double mix = f * hs + (1.0 - f) * hn;
      logL += log(mix > 0.0 ? mix : 1e-300);
    }
    f_hat = f;
    return logL;
  }

  //! Maximum-likelihood size shift over a bounded grid, with a profile-likelihood interval.
  //
  //  `direction` restricts the search to the sign the seed found, so an insertion candidate cannot
  //  be explained away as a deletion (and the grid is half the size).
  static bool pd_estimate_shift(const vector<pd_pair>& pairs, const pd_null& null, char direction,
                                int32_t max_span,
                                int32_t& delta_hat, double& f_hat, double& logL_hat,
                                int32_t& delta_lower, int32_t& delta_upper)
  {
    if (pairs.empty()) return false;

    // An event can only be seen at all while some fragment's gap still accommodates it, which caps
    // both directions at (longest usable distance - 2 * read length).
    int32_t bound = max_span - null.trunc;
    if (bound < 1) return false;
    int32_t lo = (direction == 'L') ? 1 : -bound;
    int32_t hi = (direction == 'L') ? bound : -1;

    delta_hat = 0;
    f_hat = 0.0;
    logL_hat = -numeric_limits<double>::infinity();
    for (int32_t delta = lo; delta <= hi; delta++) {
      double f = 0.0;
      double logL = pd_profile_logL(pairs, null, delta, f);
      if (logL > logL_hat) { logL_hat = logL; delta_hat = delta; f_hat = f; }
    }
    if (delta_hat == 0) return false;

    // Profile-likelihood interval: walk out either side until the log-likelihood falls by 1.92.
    delta_lower = delta_hat;
    for (int32_t delta = delta_hat - 1; delta >= lo; delta--) {
      double f = 0.0;
      if (pd_profile_logL(pairs, null, delta, f) < logL_hat - kPDLogLDrop) break;
      delta_lower = delta;
    }
    delta_upper = delta_hat;
    for (int32_t delta = delta_hat + 1; delta <= hi; delta++) {
      double f = 0.0;
      if (pd_profile_logL(pairs, null, delta, f) < logL_hat - kPDLogLDrop) break;
      delta_upper = delta;
    }
    return true;
  }

  //! Sort the pairs covering the placed breakpoint into supporting / against / ambiguous.
  //
  //  The obvious rule -- call a pair supporting when its distance is past the crossover between the
  //  two densities -- is badly biased whenever the shift is comparable to the width of the
  //  distribution. At a shift of half a standard deviation, only about 60% of the pairs from a
  //  variant that every molecule carries land on the supporting side, so a clonal event is reported
  //  at frequency 0.6.
  //
  //  A likelihood-ratio deadband fixes that and is self-regulating. A pair only counts when the two
  //  hypotheses differ by odds of kPDSupportOdds or more; otherwise it is reported as ambiguous and
  //  counted for neither. When the shift is large every pair is decisive and the frequency comes out
  //  right; when it is small almost everything is ambiguous, the counts collapse, and the
  //  Clopper-Pearson bound rejects the call on its own without a separate size threshold.
  static void pd_classify_pairs(const vector<pd_pair>& pairs, const pd_null& null,
                                int32_t b, int32_t delta,
                                int& supporting, int& against, int& ambiguous,
                                size_t& distinct, vector<int>& category)
  {
    supporting = 0; against = 0; ambiguous = 0;
    category.assign(pairs.size(), 0);
    set<pair<int32_t, int32_t> > ends;

    double H = null.total_supporting(delta);
    if (H <= 0.0) { distinct = 0; return; }

    int32_t deleted = max(delta, 0);
    for (size_t i = 0; i < pairs.size(); i++) {
      if (!pairs[i].covers(b)) continue;
      double hs = null.supporting(pairs[i].d, delta) / H;
      double hn = null.at(pairs[i].d);
      double lr = (hn > 0.0) ? (hs / hn) : 0.0;
      // A pair whose gap cannot hold the whole removed span is not evidence for the call however
      // well its distance fits, so it is never counted as supporting. Without this a deletion draws
      // support from pairs that demonstrably straddle only part of it.
      if (!pairs[i].accommodates(b, deleted) && (lr >= kPDSupportOdds)) {
        category[i] = 3;
        ambiguous++;
        continue;
      }
      if (lr >= kPDSupportOdds) {
        category[i] = 1;
        supporting++;
        // Distinct FRAGMENT ends, so PCR duplicates of one molecule cannot carry a prediction --
        // the same guard JC's pos_hash score and MP's distinct read count provide.
        ends.insert(make_pair(pairs[i].read_start, pairs[i].s2));
      } else if (lr <= 1.0 / kPDSupportOdds) {
        category[i] = 2;
        against++;
      } else {
        category[i] = 3;
        ambiguous++;
      }
    }
    distinct = ends.size();
  }

  // ---------------------------------------------------------------------------------------------
  // Localization
  // ---------------------------------------------------------------------------------------------

  //! Localize the breakpoint to the positions the supporting pairs agree on.
  //
  //  Each supporting pair licenses an interval of breakpoints. With b in between-base coordinates
  //  (b = "between reference base b and b+1") and `deleted` = max(size shift, 0):
  //
  //      b >= e1   and   b + deleted <= s2 - 1     ->   b in [e1, s2 - 1 - deleted]
  //
  //  Nothing inside that narrows it further: for a fixed set of pairs with fixed distances, the
  //  likelihood does not depend on b at all -- b enters only through whether each pair can
  //  accommodate the event -- so the surface over the interval is flat and there is no maximum to
  //  find. The intervals themselves are the whole of the information.
  //
  //  What the intervals are combined WITH matters, though. Intersecting them all is the obvious move
  //  and it is wrong, because the supporting set is not pure: near an insertion, a naturally short
  //  fragment that happens to cross the region has exactly the distance the insertion hypothesis
  //  predicts, so it is classified as supporting even though its gap stops short of the breakpoint.
  //  One such pair drags the intersection off the true position -- and, worse, drags it off while
  //  reporting a narrow range, so the answer looks confident and is wrong.
  //
  //  So take the point of MAXIMUM OVERLAP instead. Every genuinely supporting pair's interval
  //  contains the true breakpoint, while a misclassified one covers a position only by chance, so
  //  the depth of coverage peaks where the truth is and the plateau at that depth is the interval
  //  worth reporting. This degrades gracefully: with a clean supporting set it returns exactly the
  //  intersection.
  //  Which pairs are fed in matters as much as the sweep. They are selected by DISTANCE alone --
  //  every pair in the window whose distance fits the event hypothesis, whether or not it covers any
  //  particular candidate position. Selecting them by coverage of a starting guess instead would be
  //  circular in the one direction that hurts: it discards exactly the pairs whose gaps begin past
  //  the guess, which are the ones that would have pulled the answer toward the truth, and leaves an
  //  estimate anchored to wherever the seed happened to peak.
  static bool pd_feasible_interval(const vector<pd_pair>& pairs, const pd_null& null,
                                   int32_t delta, int32_t near, int32_t reach,
                                   int32_t& lo, int32_t& hi, int& depth_out)
  {
    int32_t deleted = max(delta, 0);
    double H = null.total_supporting(delta);
    if (H <= 0.0) return false;

    // Sweep endpoints: +1 where a pair's licensed interval opens, -1 just past where it closes.
    vector<pair<int32_t, int> > events;
    for (size_t i = 0; i < pairs.size(); i++) {
      double hs = null.supporting(pairs[i].d, delta) / H;
      double hn = null.at(pairs[i].d);
      if ((hn <= 0.0) || (hs / hn < kPDSupportOdds)) continue;
      int32_t a = pairs[i].e1;
      int32_t z = pairs[i].s2 - 1 - deleted;
      if (z < a) continue;   // cannot hold an event this size anywhere
      // Ignore pairs licensing only positions far from the seeded column: at a few percent of all
      // pairs, naturally short fragments are spread across the whole rescan window, and there is no
      // reason to let ones from the far end of it vote.
      if ((z < near - reach) || (a > near + reach)) continue;
      events.push_back(make_pair(a, +1));
      events.push_back(make_pair(z + 1, -1));
    }
    if (events.empty()) return false;
    sort(events.begin(), events.end());

    // Two passes: find the maximum depth, then take the FULL extent of every run reaching it.
    //
    // Reporting only one such run would be a mistake even when it is the widest. A few contaminating
    // pairs are enough to raise a short plateau to the same depth as the true one a few bases away,
    // and picking between them on width is arbitrary -- the run that wins is then reported with a
    // narrow range, which states confidence the data does not support. Spanning them instead makes
    // the range say what is true: these are the positions consistent with the most pairs. It also
    // centres the estimate between the competing plateaus rather than committing to either.
    int32_t best_lo = 0, best_hi = -1;
    int best_depth = 0;
    for (int pass = 0; pass < 2; pass++) {
      int depth = 0;
      for (size_t i = 0; i < events.size(); ) {
        int32_t at = events[i].first;
        while ((i < events.size()) && (events[i].first == at)) { depth += events[i].second; i++; }
        if (i >= events.size()) break;
        int32_t until = events[i].first - 1;   // depth holds over [at, until]
        if (until < at) continue;
        if (pass == 0) {
          if (depth > best_depth) best_depth = depth;
        } else if (depth == best_depth) {
          if (best_hi < best_lo) { best_lo = at; best_hi = until; }
          else { best_lo = min(best_lo, at); best_hi = max(best_hi, until); }
        }
      }
      if (best_depth <= 0) return false;
    }
    if (best_hi < best_lo) return false;

    lo = best_lo;
    hi = best_hi;
    depth_out = best_depth;
    return true;
  }

  //! One passing JC breakpoint, reduced to the deletion-shaped form PD reports.
  struct pd_jc_breakpoint {
    string  seq_id;
    int32_t position;   // side_1_position: last retained base of the left flank
    int32_t delta;      // side_2_position - side_1_position - 1: bases missing between the flanks
  };

  //! Load the passing (non-rejected) junction evidence, keeping only the junctions shaped like
  //  something PD could have found: both sides on one sequence, left flank kept below and right flank
  //  above. A PD whose feasible interval contains one of these can adopt its split-read coordinates,
  //  which are exact to the base where PD's are only exact to the interval.
  static void pd_load_passing_jcs(const string& file_name, vector<pd_jc_breakpoint>& out)
  {
    out.clear();
    if (!file_exists(file_name.c_str())) return;
    cGenomeDiff jc_gd(file_name);
    diff_entry_list_t jl = jc_gd.get_list(make_vector<gd_entry_type>(JC));
    for (diff_entry_list_t::iterator it = jl.begin(); it != jl.end(); it++) {
      cDiffEntry& j = **it;
      if (j.entry_exists("reject")) continue;
      if (j[SIDE_1_SEQ_ID] != j[SIDE_2_SEQ_ID]) continue;
      if (from_string<int32_t>(j[SIDE_1_STRAND]) != -1) continue;
      if (from_string<int32_t>(j[SIDE_2_STRAND]) != +1) continue;
      pd_jc_breakpoint b;
      b.seq_id = j[SIDE_1_SEQ_ID];
      b.position = from_string<int32_t>(j[SIDE_1_POSITION]);
      b.delta = from_string<int32_t>(j[SIDE_2_POSITION]) - b.position - 1;
      if (b.delta < 0) continue;
      out.push_back(b);
    }
  }

  // ---------------------------------------------------------------------------------------------

  //! A candidate carried from placement through the gates to emission.
  struct pd_call {
    string   seq_id;
    int32_t  position;         // side_1: last retained base of the left flank
    int32_t  delta;            // side_2 - side_1 - 1
    int32_t  delta_lower;
    int32_t  delta_upper;
    int32_t  position_range;   // width of the feasible interval
    int      supporting;
    int      against;
    int      ambiguous;
    size_t   distinct;
    int32_t  candidate_covering;
    double   seed_z;
    bool     snapped;
    bool     contig_end;
    bool     inconsistent;
    bool     suppressed;
  };

  } // anonymous namespace

  // ---------------------------------------------------------------------------------------------

  void predict_pair_distances(const Settings& settings, Summary& summary, cReferenceSequences& ref_seq_info)
  {
    cGenomeDiff pd_gd;

    //
    // Step 0: library geometry and the null distribution.
    //
    bool inner3p = true;
    double D = 0.0;
    double pair_median = 0.0;
    if (!paired_library_params(summary, inner3p, D, pair_median)) {
      WARN("Pair distance (PD) evidence prediction currently supports only FR- and RF-concordant "
           "libraries. No PD evidence will be predicted.");
      pd_gd.write(settings.pd_genome_diff_file_name);
      return;
    }

    // Use the read group with the most mapped pairs, matching how dp_evidence picks its model.
    const PairedMappingDistanceDistributionSummaries& pmdd = summary.preliminary_paired_mapping_distance_distribution;
    string base;
    double best_mapped = -1.0;
    double distance_cutoff = 0.0;
    for (PairedMappingDistanceDistributionSummaries::const_iterator it = pmdd.begin(); it != pmdd.end(); it++) {
      if (it->second.mapped_pairs > best_mapped) {
        best_mapped = it->second.mapped_pairs;
        base = it->first;
        distance_cutoff = it->second.distance_cutoff;
      }
    }

    pd_null null;
    int32_t read_length = static_cast<int32_t>(summary.sequence_conversion.read_length_avg + 0.5);
    if (base.empty()
        || !null.build(Settings::file_name(settings.paired_mapping_distance_histogram_file_name, "#", base),
                       read_length, 5)) {
      WARN("Pair distance (PD) evidence prediction needs the paired mapping distance histogram, which "
           "is not available. No PD evidence will be predicted.");
      pd_gd.write(settings.pd_genome_diff_file_name);
      return;
    }

    int32_t max_span = settings.pair_distance_maximum_span;
    if (max_span <= 0) max_span = static_cast<int32_t>(2.0 * distance_cutoff + 0.5);
    if (max_span <= null.trunc) max_span = null.max_distance();

    //
    // Step 1: re-read the candidate regions CSV.
    //  columns: seq_id,start,end,direction,length,peak_position,peak_covering_count,peak_tail_count,peak_z
    //
    if (!file_exists(settings.pd_candidate_regions_file_name.c_str())) {
      pd_gd.write(settings.pd_genome_diff_file_name);
      return;
    }

    vector<pd_region_row> regions;
    {
      ifstream in(settings.pd_candidate_regions_file_name.c_str());
      string line;
      getline(in, line); // header
      while (getline(in, line)) {
        if (line.empty()) continue;
        vector<string> f = split(line, ",");
        if (f.size() < 9) continue;
        pd_region_row r;
        r.seq_id = f[0];
        r.start = from_string<uint32_t>(f[1]);
        r.end = from_string<uint32_t>(f[2]);
        r.direction = f[3].empty() ? 'L' : f[3][0];
        r.peak_position = from_string<uint32_t>(f[5]);
        r.peak_covering = from_string<int32_t>(f[6]);
        r.peak_tail = from_string<int32_t>(f[7]);
        r.peak_z = from_string<double>(f[8]);
        regions.push_back(r);
      }
    }

    if (regions.empty()
        || !file_exists(settings.reference_bam_file_name.c_str())
        || !file_exists(settings.reference_fasta_file_name.c_str())) {
      pd_gd.write(settings.pd_genome_diff_file_name);
      return;
    }

    //
    // Step 1b: merge same-direction regions that describe one event. The seed test lapses and
    // recovers around a breakpoint as the covering population turns over, which leaves a strong
    // region flanked by short satellites. Merging by proximity and then keeping only the strongest
    // column is simpler and more robust than trying to stop the seed from flickering.
    //
    sort(regions.begin(), regions.end(), [](const pd_region_row& a, const pd_region_row& b) {
      if (a.seq_id != b.seq_id) return a.seq_id < b.seq_id;
      if (a.direction != b.direction) return a.direction < b.direction;
      return a.start < b.start;
    });
    {
      int32_t merge_window = max(1, static_cast<int32_t>(pair_median));
      vector<pd_region_row> merged;
      for (size_t i = 0; i < regions.size(); i++) {
        if (!merged.empty()) {
          pd_region_row& m = merged.back();
          if ((m.seq_id == regions[i].seq_id) && (m.direction == regions[i].direction)
              && (static_cast<int32_t>(regions[i].start) - static_cast<int32_t>(m.end) <= merge_window)) {
            m.end = max(m.end, regions[i].end);
            if (fabs(regions[i].peak_z) > fabs(m.peak_z)) {
              m.peak_position = regions[i].peak_position;
              m.peak_covering = regions[i].peak_covering;
              m.peak_tail = regions[i].peak_tail;
              m.peak_z = regions[i].peak_z;
            }
            continue;
          }
        }
        merged.push_back(regions[i]);
      }
      regions.swap(merged);
    }

    vector<pd_jc_breakpoint> jcs;
    pd_load_passing_jcs(settings.jc_genome_diff_file_name, jcs);

    pd_pair_scanner* scanner =
      new pd_pair_scanner(settings.reference_bam_file_name, settings.reference_fasta_file_name, max_span);

    //
    // Step 2: rescan, estimate, localize.
    //
    vector<pd_call> calls;
    for (size_t ri = 0; ri < regions.size(); ri++) {
      const pd_region_row& r = regions[ri];

      int32_t tid = scanner->tid_for_seq_id(r.seq_id);
      if (tid < 0) continue;
      int32_t seqlen = scanner->seq_length(tid);

      // Every pair whose gap can reach the peak has its leftmost mate within one maximum span of it,
      // so this window is guaranteed to re-derive the full set the seed saw.
      vector<pd_pair> pairs;
      scanner->collect(r.seq_id,
                       static_cast<int32_t>(r.peak_position) - 2 * max_span,
                       static_cast<int32_t>(r.peak_position) + 2 * max_span, pairs);

      // The size shift is estimated from the pairs covering the seeded column, which is where the
      // collective shift is strongest. Coverage of one column is the right condition HERE -- this is
      // a statement about a distance distribution at a point, and it needs an unbiased sample of it.
      vector<pd_pair> covering;
      for (size_t i = 0; i < pairs.size(); i++)
        if (pairs[i].covers(static_cast<int32_t>(r.peak_position))) covering.push_back(pairs[i]);
      if (covering.empty()) continue;

      int32_t delta = 0, delta_lower = 0, delta_upper = 0;
      double f_hat = 0.0, logL = 0.0;
      if (!pd_estimate_shift(covering, null, r.direction, max_span,
                             delta, f_hat, logL, delta_lower, delta_upper)) continue;

      // Localize over the WHOLE window, not the covering subset: see pd_feasible_interval for why
      // restricting to pairs that cover a starting guess pins the answer to that guess.
      int32_t lo = 0, hi = 0;
      int overlap_depth = 0;
      if (!pd_feasible_interval(pairs, null, delta, static_cast<int32_t>(r.peak_position),
                                max(1, static_cast<int32_t>(pair_median)), lo, hi, overlap_depth)) continue;

      int32_t b = (lo + hi) / 2;
      int32_t position_range = hi - lo;

      // Count at the placed position, again over the whole window -- pd_classify_pairs applies the
      // coverage test itself, so the denominator is every molecule sampled across the breakpoint.
      int supporting = 0, against = 0, ambiguous = 0;
      size_t distinct = 0;
      vector<int> category;
      pd_classify_pairs(pairs, null, b, delta, supporting, against, ambiguous, distinct, category);
      if (supporting == 0) continue;

      // The estimated size and the pairs' geometry disagree when most of the pairs whose distance
      // fits the event cannot all be describing one position. Recorded and rejected rather than
      // silently repaired.
      bool inconsistent = (2 * overlap_depth < supporting);

      // Snap onto a validated split-read breakpoint inside the interval whose implied size also
      // agrees. Where split reads exist they locate the event to the base, which PD cannot.
      bool snapped = false;
      if (!inconsistent) {
        for (size_t k = 0; k < jcs.size(); k++) {
          if (jcs[k].seq_id != r.seq_id) continue;
          if ((jcs[k].position < lo) || (jcs[k].position > hi)) continue;
          // The junction's own count of missing reference bases has to be one the shift estimate
          // admits. For a deletion that is the shift itself, and the junction's value -- exact to the
          // base -- replaces the estimate. For an insertion the reference loses nothing, so only a
          // flush junction is compatible, and the inserted length stays as estimated because the
          // junction says nothing about it.
          if (r.direction == 'L') {
            if ((jcs[k].delta < max(delta_lower, 1)) || (jcs[k].delta > delta_upper)) continue;
            delta = jcs[k].delta;
          } else {
            if (jcs[k].delta != 0) continue;
          }
          b = jcs[k].position;
          snapped = true;
          break;
        }
      }

      // Final counts at the placed breakpoint.
      pd_classify_pairs(covering, null, b, delta, supporting, against, ambiguous, distinct, category);

      pd_call c;
      c.seq_id = r.seq_id;
      c.position = b;
      c.delta = delta;
      c.delta_lower = delta_lower;
      c.delta_upper = delta_upper;
      c.position_range = position_range;
      c.supporting = supporting;
      c.against = against;
      c.ambiguous = ambiguous;
      c.distinct = distinct;
      c.candidate_covering = r.peak_covering;
      c.seed_z = r.peak_z;
      c.snapped = snapped;
      c.inconsistent = inconsistent;
      // Within about one fragment length of the end of a LINEAR sequence the covering population is
      // truncated by the end itself, which biases the distance distribution with no event present.
      // Marked IGNORE rather than rejected: an artifact of the reference's shape.
      c.contig_end = !ref_seq_info.is_circular(r.seq_id) &&
                     ((b <= static_cast<int32_t>(pair_median)) || (b >= seqlen - static_cast<int32_t>(pair_median)));
      c.suppressed = false;
      calls.push_back(c);
    }

    delete scanner;

    //
    // Step 3: local non-maximum suppression, in the same order-independent form MP uses.
    //
    {
      int32_t win = max(1, static_cast<int32_t>(pair_median));
      for (size_t i = 0; i < calls.size(); i++) {
        for (size_t j = 0; j < calls.size(); j++) {
          if (i == j) continue;
          if (calls[i].seq_id != calls[j].seq_id) continue;
          if ((calls[i].delta > 0) != (calls[j].delta > 0)) continue;
          if (abs(calls[i].position - calls[j].position) > win) continue;
          bool j_better = (calls[j].supporting > calls[i].supporting) ||
                          ((calls[j].supporting == calls[i].supporting) && (j < i));
          if (j_better) { calls[i].suppressed = true; break; }
        }
      }
    }

    //
    // Step 4: emit.
    //
    for (size_t i = 0; i < calls.size(); i++) {
      const pd_call& c = calls[i];

      // The two sides bracket the reference bases the event removed, so they are separated by
      // max(size_shift, 0) -- and are adjacent when nothing was removed.
      //
      // A negative shift is reported by size_shift alone, not by the sides, because PD genuinely
      // cannot place it there. A pair maps closer together by I bases whether I bases of novel
      // sequence were inserted, or k reference bases were replaced by k + I inserted ones; the
      // distances are identical. The sides therefore state the minimal reading -- nothing removed --
      // which is the same shape a junction reports for an insertion, and size_shift carries the
      // length that the sides cannot.
      int32_t deleted = max(c.delta, 0);

      cDiffEntry pd(PD);
      pd[SIDE_1_SEQ_ID] = c.seq_id;
      pd[SIDE_1_POSITION] = to_string(c.position);
      pd[SIDE_1_STRAND] = to_string(-1);
      pd[SIDE_2_SEQ_ID] = c.seq_id;
      pd[SIDE_2_POSITION] = to_string(c.position + deleted + 1);
      pd[SIDE_2_STRAND] = to_string(+1);

      if (c.contig_end) pd[IGNORE] = "CONTIG_END";

      pd[PD_SIZE_SHIFT] = to_string(c.delta);
      pd[PD_SIZE_SHIFT_LOWER] = to_string(c.delta_lower);
      pd[PD_SIZE_SHIFT_UPPER] = to_string(c.delta_upper);
      pd[PD_POSITION_RANGE] = to_string(c.position_range);
      pd[PD_SUPPORTING_COUNT] = to_string(c.supporting);
      pd[PD_AGAINST_COUNT] = to_string(c.against);
      pd[PD_AMBIGUOUS_COUNT] = to_string(c.ambiguous);
      pd[PD_DISTINCT_COUNT] = to_string(c.distinct);
      pd[PD_TOTAL_COUNT] = to_string(c.supporting + c.against + c.ambiguous);
      pd[PD_CANDIDATE_COUNT] = to_string(c.candidate_covering);
      pd[PD_SEED_Z] = to_string(formatted_double(c.seed_z, 2));
      if (c.snapped) pd["snapped_to_junction"] = "1";

      pd["side_1_annotate_key"] = "gene";
      pd["side_2_annotate_key"] = "gene";

      // The local variant frequency: of the molecules whose gap covers this point and whose distance
      // is decisive either way, what fraction carry the shift.
      double k = static_cast<double>(c.supporting);
      double n = k + static_cast<double>(c.against);
      double freq_lower = 0.0;
      if (n > 0.0) {
        freq_lower = binomial_frequency_lower_bound(k, n, kPDFrequencyAlpha);
        pd[FREQUENCY] = to_string(formatted_double(k / n, 4));
        pd[FREQUENCY_LOWER] = to_string(formatted_double(freq_lower, 4));
        pd[FREQUENCY_UPPER] = to_string(formatted_double(binomial_frequency_upper_bound(k, n, kPDFrequencyAlpha), 4));
      } else {
        pd[FREQUENCY] = "NA";
      }

      // Gates.
      if (c.suppressed)
        pd.add_reject_reason("NEARBY_BETTER_PAIR_DISTANCE");
      if (c.inconsistent)
        pd.add_reject_reason("PAIR_DISTANCE_INCONSISTENT");
      if (c.supporting < settings.pair_distance_minimum_pairs)
        pd.add_reject_reason("PAIR_DISTANCE_COUNT");
      if (static_cast<int32_t>(c.distinct) < settings.pair_distance_minimum_distinct)
        pd.add_reject_reason("PAIR_DISTANCE_DUPLICATES");
      // The size gate is a confidence statement, not a threshold on the point estimate: reject when
      // the interval for the shift still contains zero, so a call is never dropped merely for being
      // small, only for not being established. An explicit --pair-distance-minimum-shift adds an
      // absolute floor on top of that.
      if ((c.delta_lower <= 0) && (c.delta_upper >= 0))
        pd.add_reject_reason("PAIR_DISTANCE_SIZE");
      if ((settings.pair_distance_minimum_shift > 0) && (abs(c.delta) < settings.pair_distance_minimum_shift))
        pd.add_reject_reason("PAIR_DISTANCE_SIZE");
      if ((settings.pair_distance_frequency_cutoff > 0.0) && (freq_lower < settings.pair_distance_frequency_cutoff))
        pd.add_reject_reason("PAIR_DISTANCE_FREQUENCY");

      pd_gd.add(pd);
    }

    pd_gd.write(settings.pd_genome_diff_file_name);
  }

  // ---------------------------------------------------------------------------------------------
  // Evidence plot
  // ---------------------------------------------------------------------------------------------

  namespace {

  //! One read pair drawn on a PD evidence plot.
  struct pd_draw_pair {
    int     category;    // 1 = supporting, 2 = against, 3 = ambiguous
    int32_t read_start, read_end;
    int32_t mate_start, mate_end;
    int32_t d;
    bool    read_reversed;
  };

  //! Gathers the pairs to draw at a PD position, classified by the SAME likelihood ratio the counts
  //  used, so the plot shows exactly the two populations the frequency was computed from.
  //
  //  Both mates are drawn, which means both ends are needed -- and the mate's aligned end is the one
  //  thing an alignment does not carry. Rather than index the whole BAM the way DP does, the second
  //  alignment of each pair is picked up from this same window: for a pair that is being drawn at all,
  //  the mate is by construction inside it.
  class pd_plot_gatherer : public pileup_base {
  public:
    pd_plot_gatherer(const string& bam, const string& fasta)
      : pileup_base(bam, fasta), m_max_span(0) { set_print_progress(false); }

    int32_t tid_for_seq_id(const string& seq_id) const {
      for (uint32_t t = 0; t < num_targets(); t++)
        if (seq_id == string(target_name(t))) return static_cast<int32_t>(t);
      return -1;
    }
    int32_t seq_length(int32_t tid) const { return static_cast<int32_t>(target_length(tid)); }

    void gather(const string& seq_id, int32_t lo, int32_t hi, int32_t max_span)
    {
      m_max_span = max_span;
      m_lefts.clear();
      m_ends.clear();
      int32_t tid = tid_for_seq_id(seq_id);
      if (tid < 0) return;
      lo = max(1, lo);
      hi = min(seq_length(tid), hi);
      if (lo > hi) return;
      do_fetch(seq_id + ":" + to_string(lo) + "-" + to_string(hi));
    }

    void fetch_callback(const alignment_wrapper& a)
    {
      if (a.flag() & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) return;
      if (a.unmapped()) return;

      // Every primary alignment in the window is recorded by (start -> end) so a left mate can look
      // up where its partner actually stopped.
      m_ends[static_cast<int32_t>(a.reference_start_1())] = static_cast<int32_t>(a.reference_end_1());

      if (!a.is_paired()) return;
      if (a.flag() & BAM_FMUNMAP) return;
      if (a.mate_reference_target_id() != a.reference_target_id()) return;
      if (a.insert_size() <= 0) return;
      if (a.reversed() == ((a.flag() & BAM_FMREVERSE) != 0)) return;
      if (a.redundancy() > 1) return;

      pd_draw_pair p;
      p.category = 0;
      p.read_start = static_cast<int32_t>(a.reference_start_1());
      p.read_end = static_cast<int32_t>(a.reference_end_1());
      p.mate_start = a.mate_start_1();
      p.mate_end = 0;
      p.d = a.insert_size();
      p.read_reversed = a.reversed();
      if (p.mate_start - 1 < p.read_end) return;
      if (p.d > m_max_span) return;
      m_lefts.push_back(p);
    }

    //! Fill in each pair's mate end from the window index, dropping any pair whose mate was not seen
    //  (rather than guessing a full-length mate, which at a breakpoint is precisely the mate that
    //  soft-clipped and so would be drawn reaching past it).
    void finish(vector<pd_draw_pair>& out) const
    {
      out.clear();
      for (size_t i = 0; i < m_lefts.size(); i++) {
        map<int32_t, int32_t>::const_iterator e = m_ends.find(m_lefts[i].mate_start);
        if (e == m_ends.end()) continue;
        pd_draw_pair p = m_lefts[i];
        p.mate_end = e->second;
        out.push_back(p);
      }
    }

  private:
    int32_t m_max_span;
    vector<pd_draw_pair> m_lefts;
    map<int32_t, int32_t> m_ends;
  };

  static const string PD_SHIFTED_COLOR   = "'#d62728'";  // distance supports the call  = red
  static const string PD_NORMAL_COLOR    = "'#1f77b4'";  // distance argues against it  = blue
  static const string PD_AMBIGUOUS_COLOR = "'#b0b0b0'";  // cannot separate the two     = gray
  static const string PD_EVENT_COLOR     = "'#ffe0e0'";  // the reference span the event removed

  //! Render the PD evidence plot (SVG) via gnuplot: one read PAIR per lane, both mates drawn as
  //  arrows with a dashed connector across the unsequenced gap between them.
  //
  //  One continuous window holds everything, which is what makes this plot readable where DP's two
  //  separate panels would not be: PD's two sides are adjacent or nearly so, and every pair drawn
  //  spans the same point, so the left mates, the event and the right mates all fit side by side and
  //  the eye can compare connector lengths directly. That comparison IS the evidence -- the red
  //  connectors are systematically longer (a deletion) or shorter (an insertion) than the blue ones,
  //  and the frequency is just how many of each there are.
  static void render_pd_plot(const string& output_svg, vector<pd_draw_pair> pairs,
                             int32_t side_1, int32_t side_2)
  {
    // Lanes are ordered by position alone, NOT grouped by category. The staircase this produces is
    // the structure of the plot: it is what makes the left mates' ends and the right mates' starts
    // converge visibly on the breakpoint. Blocking the categories instead would buy nothing -- colour
    // already separates them -- while pulling every uncounted pair out of the staircase and stacking
    // it at one edge, where it reads as a separate population sitting somewhere else on the genome
    // rather than as one interleaved among its neighbours.
    //
    // The remaining keys only make the order total, so a run stays reproducible when two pairs share
    // a start.
    sort(pairs.begin(), pairs.end(), [](const pd_draw_pair& a, const pd_draw_pair& b) {
      if (a.read_start != b.read_start) return a.read_start < b.read_start;
      if (a.mate_end != b.mate_end)     return a.mate_end < b.mate_end;
      if (a.read_end != b.read_end)     return a.read_end < b.read_end;
      return a.category < b.category;
    });

    int n = static_cast<int>(pairs.size());
    if (n == 0) return;

    int64_t lo = side_1, hi = side_2;
    for (int i = 0; i < n; i++) {
      lo = min(lo, (int64_t)pairs[i].read_start);
      hi = max(hi, (int64_t)pairs[i].mate_end);
    }
    int64_t pad = max((int64_t)10, (int64_t)((hi - lo) / 20));
    int64_t xmin = lo - pad;
    int64_t xmax = hi + pad;

    string f_sr = output_svg + ".sr.tab";  // supporting read arrows (x y dx dy)
    string f_nr = output_svg + ".nr.tab";  // against read arrows
    string f_ar = output_svg + ".ar.tab";  // ambiguous read arrows
    string f_sc = output_svg + ".sc.tab";  // supporting connectors
    string f_nc = output_svg + ".nc.tab";  // against connectors
    string f_ac = output_svg + ".ac.tab";  // ambiguous connectors
    ofstream sr(f_sr.c_str()), nr(f_nr.c_str()), ar(f_ar.c_str());
    ofstream sc(f_sc.c_str()), nc(f_nc.c_str()), ac(f_ac.c_str());
    bool has_s = false, has_n = false, has_a = false;

    for (int i = 0; i < n; i++) {
      const pd_draw_pair& p = pairs[i];
      int y = i + 1;
      ofstream& reads = (p.category == 1) ? sr : ((p.category == 2) ? nr : ar);
      ofstream& conns = (p.category == 1) ? sc : ((p.category == 2) ? nc : ac);
      if      (p.category == 1) has_s = true;
      else if (p.category == 2) has_n = true;
      else                      has_a = true;

      // The left mate points right and the right mate points left (or the mirror under an RF
      // library); drawing each arrow from its own 5' end is what makes that orientation visible.
      int32_t l_tail = p.read_reversed ? p.read_end : p.read_start;
      int32_t l_head = p.read_reversed ? p.read_start : p.read_end;
      reads << l_tail << "\t" << y << "\t" << (l_head - l_tail) << "\t0\n";
      int32_t r_tail = p.read_reversed ? p.mate_start : p.mate_end;
      int32_t r_head = p.read_reversed ? p.mate_end : p.mate_start;
      reads << r_tail << "\t" << y << "\t" << (r_head - r_tail) << "\t0\n";
      conns << p.read_end << "\t" << y << "\n" << p.mate_start << "\t" << y << "\n\n";
    }
    sr.close(); nr.close(); ar.close(); sc.close(); nc.close(); ac.close();

    ostringstream s;
    const int kBelowPlotPx = 145;
    const int kKeyBaselinePx = 30;
    int height = max(300, n * 8 + kBelowPlotPx + 30);
    s << "set terminal svg size 1400," << height << " font ',11' noenhanced" << endl;
    s << "set output " << double_quote(output_svg) << endl;
    s << "unset title" << endl;
    s << "set bmargin at screen " << (static_cast<double>(kBelowPlotPx) / height) << endl;
    s << "set xlabel 'Reference position (bp)' font ',22'" << endl;
    s << "set xrange [" << xmin << ":" << xmax << "]" << endl;
    s << "set x2range [" << xmin << ":" << xmax << "]" << endl;
    s << "set yrange [0:" << (n + 1) << "]" << endl;
    s << "set border 15 lw 1" << endl;
    s << "set tics out" << endl;

    int64_t step = plot_nice_tick(xmax - xmin);
    string tics = plot_xtics_list(xmin, xmax, step);
    s << "set xtics (" << tics << ") nomirror font ',16'" << endl;
    s << "set x2tics (" << tics << ") font ',16'" << endl;
    s << "set lmargin " << plot_lmargin_for_labels(max(plot_commafy(xmin).size(), plot_commafy(xmax).size()), 16.0, 11.0) << endl;
    s << "unset ytics" << endl;

    // The reference span the event removed, shaded. For an insertion the two sides are adjacent and
    // there is nothing to shade, so only the breakpoint lines are drawn.
    if (side_2 > side_1 + 1) {
      s << "set object 1 rectangle from " << (side_1 + 1) << ",0 to " << (side_2 - 1) << "," << (n + 1)
        << " fc rgb " << PD_EVENT_COLOR << " fs solid noborder behind" << endl;
    }
    s << "set arrow from " << side_1 << ",0 to " << side_1 << "," << (n + 1)
      << " nohead lc rgb 'gray50' lw 1 dt 3 back" << endl;
    if (side_2 != side_1 + 1) {
      s << "set arrow from " << side_2 << ",0 to " << side_2 << "," << (n + 1)
        << " nohead lc rgb 'gray50' lw 1 dt 3 back" << endl;
    }

    // Hand-laid-out legend, centred on the canvas rather than on the plot area -- same reasoning as
    // the MP plot, where an explicit lmargin leaves gnuplot's own key visibly off-centre.
    {
      const double kCanvasW = 1400.0;
      const double kCharW   = 8.6;
      const double kSample  = 44.0;
      const double kGap     = 10.0;
      const double kBetween = 70.0;

      vector<pair<string, string> > entries;
      if (has_s) entries.push_back(make_pair(PD_SHIFTED_COLOR,   string("distance shifted (supports PD)")));
      if (has_n) entries.push_back(make_pair(PD_NORMAL_COLOR,    string("distance normal (against)")));
      if (has_a) entries.push_back(make_pair(PD_AMBIGUOUS_COLOR, string("undecidable (not counted)")));

      double total = 0.0;
      for (size_t i = 0; i < entries.size(); i++) {
        if (i) total += kBetween;
        total += kSample + kGap + entries[i].second.size() * kCharW;
      }
      double x = kCanvasW / 2.0 - total / 2.0;
      double ky = static_cast<double>(kKeyBaselinePx) / height;
      for (size_t i = 0; i < entries.size(); i++) {
        s << "set arrow from screen " << (x / kCanvasW) << ", screen " << ky
          << " to screen " << ((x + kSample) / kCanvasW) << ", screen " << ky
          << " head filled size screen 0.008,20,60 lc rgb " << entries[i].first << " lw 2 front" << endl;
        s << "set label " << double_quote(entries[i].second)
          << " at screen " << ((x + kSample + kGap) / kCanvasW) << ", screen " << ky
          << " left font ',18' front" << endl;
        x += kSample + kGap + entries[i].second.size() * kCharW + kBetween;
      }
      s << "unset key" << endl;
    }

    // Ambiguous pairs first, so the counted ones read on top of them -- but at the SAME arrow weight.
    // Every pair drawn here was sampled across the breakpoint; an undecidable one is not a lesser
    // observation, it is one the distance cannot separate, and drawing it thinner would read as less
    // real rather than as uncounted. Colour alone carries that distinction.
    vector<string> clauses;
    if (has_a) clauses.push_back(double_quote(f_ac) + " using 1:2 with lines dt 2 lc rgb " + PD_AMBIGUOUS_COLOR + " lw 1 notitle");
    if (has_n) clauses.push_back(double_quote(f_nc) + " using 1:2 with lines dt 2 lc rgb " + PD_NORMAL_COLOR + " lw 1 notitle");
    if (has_s) clauses.push_back(double_quote(f_sc) + " using 1:2 with lines dt 2 lc rgb " + PD_SHIFTED_COLOR + " lw 1 notitle");
    if (has_a) clauses.push_back(double_quote(f_ar) + " using 1:2:3:4 with vectors head filled size screen 0.008,20,60 lc rgb " + PD_AMBIGUOUS_COLOR + " lw 2 notitle");
    if (has_n) clauses.push_back(double_quote(f_nr) + " using 1:2:3:4 with vectors head filled size screen 0.008,20,60 lc rgb " + PD_NORMAL_COLOR + " lw 2 notitle");
    if (has_s) clauses.push_back(double_quote(f_sr) + " using 1:2:3:4 with vectors head filled size screen 0.008,20,60 lc rgb " + PD_SHIFTED_COLOR + " lw 2 notitle");
    if (clauses.empty()) clauses.push_back("-1 notitle");
    s << "plot " << join(clauses, string(", \\\n     ")) << endl;

    string script_name = output_svg + ".gp";
    string log_name    = output_svg + ".gp.log";
    run_gnuplot_script(s.str(), script_name, log_name);
    remove(log_name.c_str());
    remove(f_sr.c_str()); remove(f_nr.c_str()); remove(f_ar.c_str());
    remove(f_sc.c_str()); remove(f_nc.c_str()); remove(f_ac.c_str());
  }

  //! If v holds more than max_display entries, randomly down-sample it (seeded, so runs reproduce)
  //  and return the ORIGINAL size; otherwise leave v alone and return 0. Mirrors dp/mp_cap_for_display.
  static size_t pd_cap_for_display(vector<pd_draw_pair>& v, uint32_t max_display, uint32_t seed)
  {
    size_t total = v.size();
    if (max_display == 0 || total <= max_display) return 0;
    std::mt19937 rng(seed);
    std::shuffle(v.begin(), v.end(), rng);
    v.resize(max_display);
    return total;
  }

  } // anonymous namespace

  void draw_pair_distance_evidence_plots(const Settings& settings, Summary& summary,
                                         cReferenceSequences& ref_seq_info, cGenomeDiff& gd)
  {
    (void)ref_seq_info;

    bool inner3p = true;
    double D = 0.0, pair_median = 0.0;
    if (!paired_library_params(summary, inner3p, D, pair_median)) return;  // predict already warned
    if (!file_exists(settings.reference_bam_file_name.c_str()) ||
        !file_exists(settings.reference_fasta_file_name.c_str())) return;

    diff_entry_list_t pd_list = gd.get_list(make_vector<gd_entry_type>(PD));
    if (pd_list.empty()) return;

    const PairedMappingDistanceDistributionSummaries& pmdd = summary.preliminary_paired_mapping_distance_distribution;
    string base;
    double best_mapped = -1.0, distance_cutoff = 0.0;
    for (PairedMappingDistanceDistributionSummaries::const_iterator it = pmdd.begin(); it != pmdd.end(); it++) {
      if (it->second.mapped_pairs > best_mapped) {
        best_mapped = it->second.mapped_pairs; base = it->first;
        distance_cutoff = it->second.distance_cutoff;
      }
    }
    pd_null null;
    int32_t read_length = static_cast<int32_t>(summary.sequence_conversion.read_length_avg + 0.5);
    if (base.empty()
        || !null.build(Settings::file_name(settings.paired_mapping_distance_histogram_file_name, "#", base),
                       read_length, 5)) return;

    int32_t max_span = settings.pair_distance_maximum_span;
    if (max_span <= 0) max_span = static_cast<int32_t>(2.0 * distance_cutoff + 0.5);
    if (max_span <= null.trunc) max_span = null.max_distance();

    create_path(settings.evidence_path);
    pd_plot_gatherer g(settings.reference_bam_file_name, settings.reference_fasta_file_name);

    for (diff_entry_list_t::iterator it = pd_list.begin(); it != pd_list.end(); it++) {
      cDiffEntry& pd = **it;

      string  seq_id = pd[SIDE_1_SEQ_ID];
      int32_t side_1 = from_string<int32_t>(pd[SIDE_1_POSITION]);
      int32_t side_2 = from_string<int32_t>(pd[SIDE_2_POSITION]);
      int32_t delta  = pd.entry_exists(PD_SIZE_SHIFT) ? from_string<int32_t>(pd[PD_SIZE_SHIFT]) : 0;
      int32_t deleted = max(delta, 0);

      g.gather(seq_id, side_1 - 2 * max_span, side_2 + 2 * max_span, max_span);
      vector<pd_draw_pair> drawn;
      g.finish(drawn);

      // Classify with the same likelihood ratio and the same geometric requirement the counts used,
      // and keep only the pairs that reach the breakpoint at all -- exactly the set behind the
      // frequency in the table.
      double H = null.total_supporting(delta);
      vector<pd_draw_pair> keep;
      for (size_t i = 0; i < drawn.size(); i++) {
        if ((drawn[i].read_end > side_1) || (side_1 > drawn[i].mate_start - 1)) continue;
        double hs = (H > 0.0) ? (null.supporting(drawn[i].d, delta) / H) : 0.0;
        double hn = null.at(drawn[i].d);
        double lr = (hn > 0.0) ? (hs / hn) : 0.0;
        bool fits = (drawn[i].read_end <= side_1) && (side_1 + deleted <= drawn[i].mate_start - 1);
        if ((lr >= kPDSupportOdds) && fits)        drawn[i].category = 1;
        else if (lr <= 1.0 / kPDSupportOdds)       drawn[i].category = 2;
        else                                       drawn[i].category = 3;
        keep.push_back(drawn[i]);
      }
      if (keep.empty()) continue;

      size_t total = pd_cap_for_display(keep, settings.max_displayed_reads, static_cast<uint32_t>(side_1));
      if (total)
        pd["_pd_plot_message"] = "Only " + to_string(settings.max_displayed_reads) + " of " +
                                 to_string(total) + " read pairs displayed.";

      string svg = settings.evidence_path + "/PD_" + pd._id + ".svg";
      render_pd_plot(svg, keep, side_1, side_2);
      make_svg_responsive(svg);
      pd["_pd_plot_file_name"] = Settings::relative_path(svg, settings.evidence_path);
    }
  }

} // namespace breseq
