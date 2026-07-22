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

#include "dp_evidence.h"
#include "genome_diff.h"
#include "genome_diff_entry.h"
#include "pileup.h"   // pileup_base + alignment_wrapper + BAM flag macros
#include "stats.h"    // nbdtr (negative binomial CDF), incompletegamma (Poisson CDF)
#include "coverage_distribution.h" // fit_negative_binomial_histogram (small-reference DP-skew fallback)

#include <set>
#include <limits>
#include <cmath>
#include <random>
#include <algorithm>

using namespace std;

namespace breseq {

  // Minimum number of shared read pairs for a paired region to be accepted as a DP item.
  // Placeholder threshold; statistically-motivated accept/reject conditions come later.
  static const int kDPMinSharedPairs = 5;

  // ---------------------------------------------------------------------------------------------
  // Rescan of a DP junction side by direct BAM fetch.
  //
  // A DP junction side is (position p, strand s). s=+1 => the retained reference flank lies at coords
  // >= p; s=-1 => flank at coords <= p. The reads that "cross" the side sit on the kept flank with a
  // fixed strand (the same strand as the discordant reads whose region produced this side):
  //   crossing_is_forward = (inner3p == (s == -1))
  // We fetch the kept-flank window out to D = the paired-read distance_cutoff and, for reads on the
  // crossing strand whose body is on the kept side (their junction-facing "anchor" end may extend past
  // p), count three categories:
  //   (1) supporting  : discordant, mate maps to the OTHER side in that side's crossing orientation,
  //                     and (like concordant) the read itself does not extend past p
  //   (2) concordant  : proper-pair that cleanly brackets p (the reference join is still intact)
  //   (3) unpaired    : unpaired or mate-unmapped
  // A separate preliminary pass (refine_outside_median) refines p from the discordant reads' median
  // "outside" position and then nudges it past their junction-facing ends before this counting runs
  // (see predict_discordant_pairs step 4), so the "does not extend past p" guard keeps every supporting
  // read. During that refinement pass the guard is off (it needs the reads that reach past p).
  // ---------------------------------------------------------------------------------------------

  // Geometry of one junction side being scanned (plus the other side, for the supporting test).
  struct dp_side_ctx {
    int32_t p, s;                 // this side's position (1-based) and strand (+/-1)
    bool    cross_fwd;            // crossing-read strand on this side (true = forward)
    int32_t other_tid, other_p;   // the mate/other side
    bool    other_cross_fwd;      // crossing-read strand on the other side
    double  D;                    // distance_cutoff (window half-width / mate-proximity)
    int32_t ovl_p, ovl_other_p;   // reference positions for the overlapping-mate exclusion (the current
                                  // best breakpoint estimate: the initial region estimate during the
                                  // placement/gathering passes, the final placed positions for the count).
                                  // Kept separate from p/other_p -- those are the classification + fetch-
                                  // window position for this pass, which may be a re-anchored, far-from-
                                  // the-reads position that would wrongly drive the overlap test.
  };

  // One supporting read pair, viewed from both sides: outer (away-from-junction) and inner (junction-facing)
  // reference coordinates on each side. The inferred insert at (P1,P2) is reach1 + reach2 (dp_reach).
  struct dp_pair_ends { int32_t o1, i1, o2, i2; };

  // Classify one fetched read at a junction side. Returns 0 = ignore, 1 = supporting/discordant,
  // 2 = concordant-crossing, 3 = unpaired. `anchor` is set (for kept reads) to the junction-facing
  // reference coordinate (reference_end for s=-1, reference_start for s=+1). This is the single source
  // of truth shared by the counting scanner and the plotting gatherer, so the plot shows exactly the
  // reads that are counted.
  static int dp_classify_side_read(const alignment_wrapper& a, const dp_side_ctx& c, int32_t& anchor)
  {
    anchor = 0;
    if (a.flag() & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) return 0;
    if (a.unmapped()) return 0;

    // Only reads on this side's crossing strand.
    bool is_forward = !a.reversed();
    if (is_forward != c.cross_fwd) return 0;

    int32_t rstart = static_cast<int32_t>(a.reference_start_1());
    int32_t rend   = static_cast<int32_t>(a.reference_end_1());

    // The read's body must be on the kept side; the junction-facing "anchor" may extend past p.
    if (c.s == -1) { if (rstart > c.p) return 0; anchor = rend; }
    else           { if (rend   < c.p) return 0; anchor = rstart; }

    bool paired = a.is_paired();
    bool mate_unmapped = (a.flag() & BAM_FMUNMAP) != 0;
    if (!paired || mate_unmapped) return 3;   // unpaired

    if (!a.proper_pair()) {
      // Discordant: the mate must land at the OTHER side in its crossing orientation.
      if (static_cast<int32_t>(a.mate_reference_target_id()) != c.other_tid) return 0;
      int32_t mpos = a.mate_start_1();
      if (mpos < c.other_p - c.D || mpos > c.other_p + c.D) return 0;
      bool mate_forward = (a.flag() & BAM_FMREVERSE) == 0;
      if (mate_forward != c.other_cross_fwd) return 0;
      // A mate pair whose two reads OVERLAP cannot support a discordant junction: the fragment is shorter
      // than the two reads combined, so there is no un-sequenced gap for a breakpoint to sit in. Measure
      // the inner gap when the pair is joined across the junction = (this read's aligned inner end to p) +
      // (the mate's aligned inner end to other_p). A negative sum means the reads overlap -> drop the pair.
      // The mate's junction-facing end follows its orientation (forward -> its right/rend faces the
      // junction, reverse -> its left/rstart), its length approximated by this read's length. The
      // inner gaps are measured against the overlap reference (ovl_p/ovl_other_p = the current best
      // breakpoint estimate), NOT the classification/window position p, which may be re-anchored.
      {
        int32_t matelen = static_cast<int32_t>(a.read_length());
        int32_t g_this  = (c.s == -1) ? (c.ovl_p - rend) : (rstart - c.ovl_p);
        int32_t g_other = mate_forward ? (c.ovl_other_p - (mpos + matelen - 1)) : (mpos - c.ovl_other_p);
        if (g_this + g_other < 0) return 0;
      }
      return 1;
    }
    // Concordant: count only if the pair cleanly BRACKETS p — this read entirely on the kept side and
    // its mate entirely on the other side, with NEITHER read overlapping p. (Discordant reads above may
    // legitimately extend past p; concordant reads may not.)
    //   s=-1: kept side = coords < p, other side = coords > p. This read (forward) must end before p
    //         (rend < p); its reverse mate's leftmost (mate_start) is nearest p, so mate_start > p means
    //         the whole mate is > p.
    //   s=+1: kept side = coords > p, other side = coords < p. This read (reverse) must start after p
    //         (rstart > p); its forward mate's rightmost (mate_end) is nearest p, so mate_start+mate_len-1
    //         < p means the whole mate is < p. The mate CIGAR isn't in this record, so its far extent is
    //         approximated by this read's aligned length (reads are ~equal length).
    bool kept_clear = (c.s == -1) ? (rend < c.p) : (rstart > c.p);
    if (!kept_clear) return 0;
    int32_t mpos = a.mate_start_1();
    int32_t mate_len = static_cast<int32_t>(a.reference_end_1()) - static_cast<int32_t>(a.reference_start_1()) + 1;
    bool completely_other = (c.s == -1) ? (mpos > c.p)
                                        : (mpos + mate_len - 1 < c.p);
    return completely_other ? 2 : 0;
  }

  // BAM target id for a seq_id (-1 if not present).
  static int32_t dp_tid_for_seq_id(const pileup_base& pb, const string& seq_id) {
    for (uint32_t t = 0; t < pb.num_targets(); t++)
      if (seq_id == string(pb.target_name(t))) return static_cast<int32_t>(t);
    return -1;
  }

  // Compute this side's one-sided fetch window [lo, hi] (kept-flank side, out to D).
  static bool dp_side_window(const pileup_base& pb, const string& seq_id, int32_t p, int32_t s, double D,
                             int32_t& lo, int32_t& hi) {
    int32_t tid = dp_tid_for_seq_id(pb, seq_id);
    if (tid < 0) return false;
    int32_t seqlen = static_cast<int32_t>(pb.target_length(tid));
    if (s == -1) { lo = max(1, static_cast<int32_t>(p - D)); hi = min(seqlen, p); }
    else         { lo = max(1, p);                           hi = min(seqlen, static_cast<int32_t>(p + D)); }
    return lo <= hi;
  }

  // A read pair's identifying number, shared by its two mates (the two paired read files write names
  // "<file-prefix>:<read-number>"; the prefix differs between mates, the number is common). Used to
  // pair a read at one junction side with its mate at the other side (both count/plot logic rely on it).
  static string dp_read_num(const string& name) {
    size_t colon = name.find(':');
    return (colon == string::npos) ? name : name.substr(colon + 1);
  }

  // Counts the three read categories at a junction side (used to fill the DP evidence fields), and
  // provides a preliminary refinement pass over the discordant reads at a side.
  class dp_side_scanner : public pileup_base {
  public:
    dp_side_scanner(const string& bam, const string& fasta)
      : pileup_base(bam, fasta), m_collect_outside(false), m_collect_pairs(false),
        m_gother_s(0), m_gread_len(0) { set_print_progress(false); }

    int32_t tid_for_seq_id(const string& seq_id) const { return dp_tid_for_seq_id(*this, seq_id); }
    int32_t seq_length(int32_t tid) const { return static_cast<int32_t>(target_length(tid)); }

    void scan(const string& seq_id, int32_t p, int32_t s, bool crossing_is_forward,
              int32_t other_tid, int32_t other_p, bool other_crossing_is_forward, double D,
              int32_t ovl_p, int32_t ovl_other_p)
    {
      set_ctx(p, s, crossing_is_forward, other_tid, other_p, other_crossing_is_forward, D, ovl_p, ovl_other_p);
      m_supporting = 0; m_concordant = 0; m_unpaired = 0;
      m_supporting_nums.clear();
      m_collect_outside = false;

      int32_t lo, hi;
      if (!dp_side_window(*this, seq_id, p, s, D, lo, hi)) return;
      do_fetch(seq_id + ":" + to_string(lo) + "-" + to_string(hi));
    }

    // Fetch the discordant (supporting) reads at this side and return two things:
    //   median_outer = the median of their "outside" coordinates (the end facing AWAY from the junction:
    //     rstart for s=-1, rend for s=+1). The caller shifts this toward the junction by half the median
    //     pair distance to place the breakpoint.
    //   inner_edge = the furthest a supporting read reaches TOWARD the junction (its junction-facing end
    //     INCLUDING its junction-facing soft-clip: rend+trail for s=-1, rstart-lead for s=+1). The caller
    //     clamps the placed coordinate out to this so no supporting read straddles the seam.
    //   extreme_outer = the OUTERMOST outer coordinate (min for s=-1, max for s=+1) -- the supporting read
    //     whose body reaches furthest AWAY from the junction. Used to re-anchor the wide window for the
    //     refinement pass.
    // Returns false (no placement) if no supporting reads are found. out_count = # supporting reads.
    bool supporting_outer_median(const string& seq_id, int32_t p, int32_t s, bool crossing_is_forward,
                                 int32_t other_tid, int32_t other_p, bool other_crossing_is_forward,
                                 double D, int32_t ovl_p, int32_t ovl_other_p, int32_t& median_outer,
                                 int32_t& inner_edge, int32_t& extreme_outer, size_t& out_count)
    {
      set_ctx(p, s, crossing_is_forward, other_tid, other_p, other_crossing_is_forward, D, ovl_p, ovl_other_p);
      m_outside.clear();
      m_have_inner = false; m_inner_edge = 0;
      m_collect_outside = true;

      int32_t lo = p, hi = p;
      if (dp_side_window(*this, seq_id, p, s, D, lo, hi))
        do_fetch(seq_id + ":" + to_string(lo) + "-" + to_string(hi));

      m_collect_outside = false;
      out_count = m_outside.size();
      if (m_outside.empty()) return false;
      sort(m_outside.begin(), m_outside.end());
      median_outer = m_outside[m_outside.size() / 2];
      // Outermost supporting read WITHIN the fetch window. do_fetch also returns reads that merely OVERLAP
      // the window (starting before its outer boundary); those outliers must not anchor the refinement's
      // re-anchored window, or it would sit too far out and miss the junction-facing cluster. So take the
      // extreme outer coordinate that still lies inside [lo, hi].
      if (s == -1) {
        vector<int32_t>::iterator it = lower_bound(m_outside.begin(), m_outside.end(), lo);
        extreme_outer = (it != m_outside.end()) ? *it : lo;
      } else {
        vector<int32_t>::iterator it = upper_bound(m_outside.begin(), m_outside.end(), hi);
        extreme_outer = (it != m_outside.begin()) ? *(it - 1) : hi;
      }
      inner_edge = m_inner_edge;
      return true;
    }

    // Gather the supporting read pairs at this side over a SYMMETRIC +/-D window (so reads on either side
    // of the region position are seen -- outliers included, to be judged by the Bayes test). Each pair
    // records both mates' outer and inner (junction-facing) reference coordinates. The mate's ends are
    // approximated from its start + read_len (its CIGAR isn't in this record).
    bool gather_pairs(const string& seq_id, int32_t p, int32_t s, bool crossing_is_forward,
                      int32_t other_tid, int32_t other_p, bool other_crossing_is_forward, double D,
                      int32_t ovl_p, int32_t ovl_other_p,
                      int32_t other_s, int32_t read_len, vector<dp_pair_ends>& out) {
      set_ctx(p, s, crossing_is_forward, other_tid, other_p, other_crossing_is_forward, D, ovl_p, ovl_other_p);
      m_pairs.clear(); m_gother_s = other_s; m_gread_len = read_len; m_collect_pairs = true;
      int32_t tid = tid_for_seq_id(seq_id);
      if (tid >= 0) {
        int32_t lo = max(1, static_cast<int32_t>(p - D)), hi = min(seq_length(tid), static_cast<int32_t>(p + D));
        if (lo <= hi) do_fetch(seq_id + ":" + to_string(lo) + "-" + to_string(hi));
      }
      m_collect_pairs = false;
      out = m_pairs;
      return !out.empty();
    }

    void fetch_callback(const alignment_wrapper& a) {
      if (m_collect_pairs) {
        // Minimal supporting classification (no body-side gate): crossing strand, discordant, mate at the
        // other side in its crossing orientation.
        if (a.flag() & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) return;
        if (a.unmapped() || !a.is_paired() || a.proper_pair()) return;
        if ((!a.reversed()) != m_ctx.cross_fwd) return;
        if (static_cast<int32_t>(a.mate_reference_target_id()) != m_ctx.other_tid) return;
        int32_t mpos = a.mate_start_1();
        if (mpos < m_ctx.other_p - m_ctx.D || mpos > m_ctx.other_p + m_ctx.D) return;
        bool mate_forward = (a.flag() & BAM_FMREVERSE) == 0;
        if (mate_forward != m_ctx.other_cross_fwd) return;
        int32_t rs = static_cast<int32_t>(a.reference_start_1()), re = static_cast<int32_t>(a.reference_end_1());
        // Overlapping-mate exclusion, referenced to the current best breakpoint estimate (ovl_p/
        // ovl_other_p), mirroring dp_classify_side_read: drop a pair whose two reads would overlap.
        {
          int32_t matelen = static_cast<int32_t>(a.read_length());
          int32_t g_this  = (m_ctx.s == -1) ? (m_ctx.ovl_p - re) : (rs - m_ctx.ovl_p);
          int32_t g_other = mate_forward ? (m_ctx.ovl_other_p - (mpos + matelen - 1)) : (mpos - m_ctx.ovl_other_p);
          if (g_this + g_other < 0) return;
        }
        dp_pair_ends e;
        e.o1 = (m_ctx.s == -1) ? rs : re;                    e.i1 = (m_ctx.s == -1) ? re : rs;
        e.o2 = (m_gother_s == -1) ? mpos : (mpos + m_gread_len - 1);
        e.i2 = (m_gother_s == -1) ? (mpos + m_gread_len - 1) : mpos;
        m_pairs.push_back(e);
        return;
      }
      int32_t anchor;
      int cat = dp_classify_side_read(a, m_ctx, anchor);
      if (m_collect_outside) {
        if (cat == 1) {
          int32_t rstart = static_cast<int32_t>(a.reference_start_1());
          int32_t rend   = static_cast<int32_t>(a.reference_end_1());
          m_outside.push_back(m_ctx.s == -1 ? rstart : rend);   // outside (away-from-junction) end
          // Junction-facing ALIGNED end (soft-clipped bases excluded -- they map to the other side of the
          // junction); keep the most extreme (furthest toward the junction) as the clamp edge.
          int32_t je = (m_ctx.s == -1) ? rend : rstart;
          if (!m_have_inner || (m_ctx.s == -1 ? je > m_inner_edge : je < m_inner_edge)) {
            m_have_inner = true; m_inner_edge = je;
          }
        }
        return;
      }
      if      (cat == 1) { m_supporting++; m_supporting_nums.insert(dp_read_num(a.read_name())); }
      else if (cat == 2) m_concordant++;
      else if (cat == 3) m_unpaired++;
    }

    int supporting() const { return m_supporting; }
    // Read-pair numbers of the discordant (supporting) reads at the last-scanned side. Intersecting the
    // two sides' sets gives the true count of pairs that bridge THIS junction (a read at a breakpoint
    // shared with a neighboring junction appears on only one side and is excluded).
    const set<string>& supporting_nums() const { return m_supporting_nums; }
    int concordant() const { return m_concordant; }
    int unpaired()   const { return m_unpaired; }

  private:
    void set_ctx(int32_t p, int32_t s, bool crossing_is_forward,
                 int32_t other_tid, int32_t other_p, bool other_crossing_is_forward, double D,
                 int32_t ovl_p, int32_t ovl_other_p) {
      m_ctx.p = p; m_ctx.s = s; m_ctx.cross_fwd = crossing_is_forward;
      m_ctx.other_tid = other_tid; m_ctx.other_p = other_p; m_ctx.other_cross_fwd = other_crossing_is_forward;
      m_ctx.D = D; m_ctx.ovl_p = ovl_p; m_ctx.ovl_other_p = ovl_other_p;
    }
    dp_side_ctx m_ctx;
    int     m_supporting, m_concordant, m_unpaired;
    set<string> m_supporting_nums;
    bool    m_collect_outside;
    // Supporting reads' outside (away-from-junction) coordinates, accumulated during the collect pass,
    // plus the furthest junction-facing edge (soft-clip included) seen among them.
    vector<int32_t> m_outside;
    bool    m_have_inner; int32_t m_inner_edge;
    // Pair-gathering pass (gather_pairs): both mates' ends, plus the other side's strand + read length.
    bool    m_collect_pairs; vector<dp_pair_ends> m_pairs; int32_t m_gother_s, m_gread_len;
  };

  // One read to draw on a per-side plot (its pair anchored at this side).
  struct dp_draw_read {
    int      category;      // 1 = discordant/supporting, 2 = concordant
    int32_t  read_start, read_end;
    bool     read_reversed;
    int32_t  anchor;        // junction-facing coordinate of this read
    int32_t  mate_start;    // concordant only
    bool     mate_reversed; // concordant only
    string   name;          // read (pair) QNAME, shown as the lane label
  };

  // Gathers (a SEPARATE fetch pass from dp_side_scanner) the concordant/discordant reads at a side for
  // plotting, using the SAME classification so the plot matches the counts.
  class dp_side_plot_gatherer : public pileup_base {
  public:
    dp_side_plot_gatherer(const string& bam, const string& fasta)
      : pileup_base(bam, fasta) { set_print_progress(false); }

    int32_t tid_for_seq_id(const string& seq_id) const { return dp_tid_for_seq_id(*this, seq_id); }

    void gather(const string& seq_id, int32_t p, int32_t s, bool crossing_is_forward,
                int32_t other_tid, int32_t other_p, bool other_crossing_is_forward, double D)
    {
      m_ctx.p = p; m_ctx.s = s; m_ctx.cross_fwd = crossing_is_forward;
      m_ctx.other_tid = other_tid; m_ctx.other_p = other_p; m_ctx.other_cross_fwd = other_crossing_is_forward;
      m_ctx.D = D;
      // Same classification as the count, so the plot matches. The plot gathers at the final placed
      // positions (read from the .gd), so the overlap guard references those same positions.
      m_ctx.ovl_p = p; m_ctx.ovl_other_p = other_p;
      m_reads.clear();

      int32_t lo, hi;
      if (!dp_side_window(*this, seq_id, p, s, D, lo, hi)) return;
      do_fetch(seq_id + ":" + to_string(lo) + "-" + to_string(hi));
    }

    void fetch_callback(const alignment_wrapper& a) {
      int32_t anchor;
      int cat = dp_classify_side_read(a, m_ctx, anchor);
      if (cat != 1 && cat != 2) return;
      dp_draw_read r;
      r.category = cat;
      r.read_start = static_cast<int32_t>(a.reference_start_1());
      r.read_end   = static_cast<int32_t>(a.reference_end_1());
      r.read_reversed = a.reversed();
      r.anchor = anchor;
      r.mate_start = a.mate_start_1();
      r.mate_reversed = (a.flag() & BAM_FMREVERSE) != 0;
      r.name = a.read_name();
      m_reads.push_back(r);
    }

    const vector<dp_draw_read>& reads() const { return m_reads; }

  private:
    dp_side_ctx m_ctx;
    vector<dp_draw_read> m_reads;
  };

  // One candidate region parsed from DP_candidate_regions.csv.
  struct dp_region_row {
    string      seq_id;
    uint32_t    start;   // lower coordinate of the region span
    uint32_t    end;     // higher coordinate
    char        strand;  // focal-read strand: 'F' or 'R'
    set<string> keys;    // distinct <read1>__<read2>__<insert_size> pair keys in this region
  };

  // Convert one region into a JC-style junction side (position, strand), given whether the library's
  // reads face the junction with their 3' end (inner3p). See dp_evidence plan / the coordinate
  // convention: strand=+1 means the retained flank lies at coords >= position, -1 at <= position.
  //   inner3p:  F -> (end, -1) ; R -> (begin, +1)
  //  !inner3p:  F -> (begin, +1) ; R -> (end, -1)   [the RF/"outie" mirror]
  static void dp_region_to_side(const dp_region_row& r, bool inner3p, int32_t& position, int32_t& strand)
  {
    bool is_forward = (r.strand == 'F');
    if (inner3p == is_forward) {
      position = static_cast<int32_t>(r.end);
      strand = -1;
    } else {
      position = static_cast<int32_t>(r.start);
      strand = +1;
    }
  }

  // Determine the library's concordant orientation (which fixes which read end faces the junction ->
  // inner3p) and the rescan window half-width D = max distance_cutoff over paired groups. FR is fully
  // supported; RF is the mirror; FF/RR is not supported yet. Returns false (optionally warning) if the
  // orientation is unsupported/unknown.
  static bool dp_library_params(const Summary& summary, bool& inner3p, double& D, double& pair_median, bool warn)
  {
    string majority_orientation;
    D = 0.0;
    pair_median = 0.0;
    map<string, int> votes;
    for (PairedMappingDistanceDistributionSummaries::const_iterator it = summary.preliminary_paired_mapping_distance_distribution.begin();
         it != summary.preliminary_paired_mapping_distance_distribution.end(); it++) {
      if (!it->second.majority_orientation.empty())
        votes[it->second.majority_orientation]++;
      if (it->second.distance_cutoff > D) D = it->second.distance_cutoff;
      if (it->second.median > pair_median) pair_median = it->second.median;
    }
    int best = 0;
    for (map<string, int>::iterator it = votes.begin(); it != votes.end(); it++) {
      if (it->second > best) { best = it->second; majority_orientation = it->first; }
    }

    if (majority_orientation == "FR") { inner3p = true;  return true; }
    if (majority_orientation == "RF") { inner3p = false; return true; }
    if (warn) {
      WARN("Discordant pair (DP) evidence prediction currently supports only FR- and RF-concordant "
           "libraries. The library concordant orientation is '" +
           (majority_orientation.empty() ? string("unknown") : majority_orientation) +
           "'; no DP evidence will be predicted.");
    }
    return false;
  }

  static const double kDPMaxScore = 999999.0;

  // ---------------------------------------------------------------------------------------------
  // Bayesian outlier test for DP position shifts.
  //
  // A supporting pair, with the junction joined at (P1,P2), implies an insert size
  //   d = reach1 + reach2 = (this mate's outer end -> P1) + (its mate's outer end -> P2).
  // Under a mixture the pair's likelihood is L = (1-p_out)*f(d) + p_out*u, where f is the empirical
  // insert PMF (overlapping-pair distances truncated out, as DP does), u is a uniform outlier
  // likelihood, and p_out is the prior that a chance discordant pair lands in this DP's windows.
  // A lone read whose own inferred insert is an outlier (BF_read = (1-p_out)f(d)/(p_out u) < 1/3) at
  // the position set by the next read is dropped, so a one-off read cannot drag side_x_position.
  // ---------------------------------------------------------------------------------------------
  struct dp_insert_model {
    vector<double> counts;   // counts[d] = # concordant pairs with mapping distance d (majority orient.)
    double  total;           // sum of counts
    int32_t trunc;           // 2 * read_length: distances below this (overlapping mates) are excluded
    double  u;               // outlier insert likelihood (uniform over the inferred-insert range)
    double  p_out;           // prior probability of a chance/outlier supporting read
    dp_insert_model() : total(0.0), trunc(0), u(0.0), p_out(0.0) {}
    bool ok() const { return total > 0.0 && trunc > 0; }
    // Empirical PMF f(d): +/-5 bp triangular-free smoothing, a small floor, and 0 below the truncation.
    double f(int32_t d) const {
      double floor = 1.0 / (total * 1000.0);
      if (d < trunc) return floor;
      double s = 0.0;
      for (int32_t x = d - 5; x <= d + 5; x++)
        if (x >= 0 && x < static_cast<int32_t>(counts.size())) s += counts[x];
      double p = s / (total * 11.0);
      return p > floor ? p : floor;
    }
    double BF_read(int32_t d) const { return ((1.0 - p_out) * f(d)) / (p_out * u); }  // inlier:outlier odds
  };

  // Read the persisted majority-orientation distance histogram (distance<TAB>count) into m.counts / m.total.
  static bool dp_read_insert_hist(const string& fn, dp_insert_model& m) {
    m.counts.clear(); m.total = 0.0;
    if (!file_exists(fn.c_str())) return false;
    ifstream in(fn.c_str());
    string line;
    while (getline(in, line)) {
      if (line.empty()) continue;
      vector<string> f = split(line, "\t");
      if (f.size() < 2) continue;
      int d = from_string<int>(f[0]);
      double n = from_string<double>(f[1]);
      if (d < 0) continue;
      if (static_cast<int>(m.counts.size()) <= d) m.counts.resize(d + 1, 0.0);
      m.counts[d] = n; m.total += n;
    }
    return m.total > 0.0;
  }

  static int32_t dp_reach(int32_t P, int32_t outer, int32_t s) { return (s == -1) ? (P - outer) : (outer - P); }

  // Robust innermost edge on ONE side: process that side's reads innermost-first and drop the innermost
  // whenever ITS inferred insert (at the position set by the next distinct read) is an outlier; stop at the
  // first read that fits and return its inner edge. `this_side` selects which side's inner edge varies.
  static int32_t dp_robust_edge(vector<dp_pair_ends> pr, bool this_is_side1,
                                int32_t s_this, int32_t p_other, int32_t s_other, const dp_insert_model& m) {
    // innermost-first sort by this side's inner edge
    if (this_is_side1) sort(pr.begin(), pr.end(), [&](const dp_pair_ends&a,const dp_pair_ends&b){ return s_this==-1 ? a.i1>b.i1 : a.i1<b.i1; });
    else               sort(pr.begin(), pr.end(), [&](const dp_pair_ends&a,const dp_pair_ends&b){ return s_this==-1 ? a.i2>b.i2 : a.i2<b.i2; });
    auto inner = [&](const dp_pair_ends& p){ return this_is_side1 ? p.i1 : p.i2; };
    size_t i = 0;
    while (i + 1 < pr.size()) {
      size_t j = i + 1;
      while (j < pr.size() && inner(pr[j]) == inner(pr[i])) j++;   // next DISTINCT inner edge
      if (j >= pr.size()) break;
      int32_t p_rest = inner(pr[j]);
      int32_t reach_this  = this_is_side1 ? dp_reach(p_rest, pr[i].o1, s_this) : dp_reach(p_rest, pr[i].o2, s_this);
      int32_t reach_other = this_is_side1 ? dp_reach(p_other, pr[i].o2, s_other) : dp_reach(p_other, pr[i].o1, s_other);
      if (m.BF_read(reach_this + reach_other) >= 1.0 / 3.0) break;  // innermost read fits -> keep
      i = j;                                                        // outlier -> peel to next distinct edge
    }
    return inner(pr[i]);
  }

  // Load the empirical insert PMF for the main paired library + derive the mixture's u and p_out.
  //   u     = 1/(2*distance_cutoff): uniform outlier likelihood over the inferred-insert range.
  //   p_out = 1 - (1 - (w/G)^2)^N: prior that >=1 chance discordant pair lands in this DP's windows,
  //           with w = region window (median + 2.42*MAD), G = total reference length, N = genome-wide
  //           discordant pairs (mapped - concordant). Returns false (no Bayes test) if unavailable.
  static bool dp_load_insert_model(const Settings& settings, const Summary& summary,
                                   cReferenceSequences& ref_seq_info, dp_insert_model& m)
  {
    const PairedMappingDistanceDistributionSummaries& pmdd = summary.preliminary_paired_mapping_distance_distribution;
    string base; double best_mapped = -1.0, median = 0, mad = 0, dcut = 0, mapped = 0, concord = 0;
    for (PairedMappingDistanceDistributionSummaries::const_iterator it = pmdd.begin(); it != pmdd.end(); it++) {
      if (it->second.mapped_pairs > best_mapped) {
        best_mapped = it->second.mapped_pairs; base = it->first;
        median = it->second.median; mad = it->second.mad; dcut = it->second.distance_cutoff;
        mapped = it->second.mapped_pairs; concord = it->second.concordant_pairs;
      }
    }
    if (base.empty()) return false;
    string fn = Settings::file_name(settings.paired_mapping_distance_histogram_file_name, "#", base);
    if (!dp_read_insert_hist(fn, m)) return false;
    int32_t readlen = static_cast<int32_t>(summary.sequence_conversion.read_length_avg + 0.5);
    m.trunc = 2 * readlen;
    m.u = (dcut > 0.0) ? 1.0 / (2.0 * dcut) : 0.0;
    double w = median + 2.42 * mad;
    double G = static_cast<double>(ref_seq_info.get_total_length());
    double Nd = mapped - concord;
    double per = (G > 0.0) ? (w / G) * (w / G) : 0.0;
    m.p_out = (Nd > 0.0) ? (1.0 - pow(1.0 - per, Nd)) : 0.0;
    return m.ok() && m.u > 0.0 && m.p_out > 0.0;
  }

  // Read a per-seq_id interior crossing histogram tab (crossing<TAB>count) into `hist` indexed by
  // crossing value. Returns false if absent/empty.
  static bool dp_read_crossing_hist(const string& fn, vector<double>& hist)
  {
    hist.clear();
    if (!file_exists(fn.c_str())) return false;
    ifstream in(fn.c_str());
    string line;
    getline(in, line);  // header "crossing\tcount"
    double total = 0.0;
    while (getline(in, line)) {
      if (line.empty()) continue;
      vector<string> f = split(line, "\t");
      if (f.size() < 2) continue;
      int c = from_string<int>(f[0]);
      double n = from_string<double>(f[1]);
      if (c < 0) continue;
      if (static_cast<int>(hist.size()) <= c) hist.resize(c + 1, 0.0);
      hist[c] = n; total += n;
    }
    return total > 0.0;
  }

  // Average coverage for a seq_id (unique-only fit if available, else the preliminary preprocess value).
  // Only the RATIO between sequences matters for the crossing coverage-projection, so any consistent
  // average-coverage measure works; unique_coverage is preferred (more accurate) and available at the
  // DP scoring and output-plot stages.
  static double dp_seq_coverage(const Summary& summary, const string& seq_id)
  {
    CoverageSummaries::const_iterator u = summary.unique_coverage.find(seq_id);
    if (u != summary.unique_coverage.end() && u->second.average > 0.0) return u->second.average;
    CoverageSummaries::const_iterator p = summary.preprocess_coverage.find(seq_id);
    if (p != summary.preprocess_coverage.end() && p->second.average > 0.0) return p->second.average;
    return 0.0;
  }

  // Censor window [lo, hi] on the crossing distribution to the NORMAL bulk -- excluding the near-zero
  // deletion spike (deletions/no-coverage positions, which otherwise pollute the lower tail) and the
  // high repeat outliers. Peak (mode) is found ignoring the crossing=0 bin; window = [0.25, 2.5]*peak,
  // matching the "normal support" used in offline validation.
  static void dp_crossing_censor(const vector<double>& hist, int32_t& lo, int32_t& hi)
  {
    // Find the BULK peak (mode) on a 5-point moving average, ignoring the crossing=0 bin -- so an
    // isolated low-crossing spike (e.g. deletion-edge positions at crossing 1) doesn't masquerade as
    // the mode.
    int32_t n = static_cast<int32_t>(hist.size());
    double peak = 0.0; int32_t peak_x = 1;
    for (int32_t i = 1; i < n; i++) {
      double sm = 0.0; int32_t cnt = 0;
      for (int32_t j = max(1, i - 2); j <= min(n - 1, i + 2); j++) { sm += hist[j]; cnt++; }
      sm = (cnt > 0) ? sm / cnt : 0.0;
      if (sm > peak) { peak = sm; peak_x = i; }
    }
    lo = max(1, static_cast<int32_t>(0.25 * peak_x));
    hi = min(n - 1, static_cast<int32_t>(2.5 * peak_x));
    if (hi < lo) { lo = 1; hi = n - 1; }
  }

  // The run-wide reference crossing distribution = the histogram of the LONGEST sequence (best-sampled,
  // most interior positions), read from its tab; C_ref = its preliminary average coverage; [lo,hi] =
  // its normal-bulk censor window. Returns false if unavailable.
  static bool dp_load_crossing_reference(const Settings& settings, const Summary& summary,
                                         cReferenceSequences& ref_seq_info,
                                         vector<double>& hist_ref, double& C_ref, string& ref_seq_id,
                                         int32_t& censor_lo, int32_t& censor_hi,
                                         double& N, bool& use_empirical, double& nb_size, double& nb_mu)
  {
    N = 0.0; use_empirical = true; nb_size = 0.0; nb_mu = 0.0;
    ref_seq_id = ""; size_t best_len = 0;
    for (cReferenceSequences::iterator it = ref_seq_info.begin(); it != ref_seq_info.end(); it++)
      if (it->get_sequence_length() > best_len) { best_len = it->get_sequence_length(); ref_seq_id = it->m_seq_id; }
    if (ref_seq_id.empty()) return false;
    C_ref = dp_seq_coverage(summary, ref_seq_id);
    if (!(C_ref > 0.0)) return false;
    string fn = Settings::file_name(settings.concordant_pair_crossing_distribution_file_name, "#", ref_seq_id);
    if (!dp_read_crossing_hist(fn, hist_ref)) return false;
    dp_crossing_censor(hist_ref, censor_lo, censor_hi);

    // N = # non-deletion (censor-window) positions in the reference crossing distribution. The empirical
    // distribution is used only when its resolution ceiling (~log10(N)) clears the DP skew cutoff by a
    // fixed 1-decade margin; smaller references fall back to a negative-binomial fit whose parametric
    // tail is not floored at 0.5/N.
    for (int32_t c = censor_lo; c <= censor_hi && c < static_cast<int32_t>(hist_ref.size()); c++) N += hist_ref[c];
    use_empirical = (N > 0.0) && (log10(N) >= settings.discordant_pair_skew_cutoff + 1.0);
    if (!use_empirical && !hist_ref.empty())
      fit_negative_binomial_histogram(hist_ref, static_cast<uint32_t>(hist_ref.size() - 1), nb_size, nb_mu);
    return true;
  }

  // Negative-binomial pmf P(X=c) parametrized by (size, mu): mean mu, variance mu + mu^2/size.
  static double dp_nbinom_pmf(int32_t c, double size, double mu)
  {
    if (!(size > 0.0) || !(mu > 0.0) || c < 0) return 0.0;
    double prob = size / (size + mu);   // P(success); mean of failures = size*(1-prob)/prob = mu
    return exp(lgamma(static_cast<double>(c) + size) - lgamma(size) - lgamma(static_cast<double>(c) + 1.0)
               + size * log(prob) + static_cast<double>(c) * log(1.0 - prob));
  }

  // The per-position thinning CDF: probability that a reference position with crossing c, at relative
  // coverage r, is spanned by <= k concordant pairs. Binomial thinning for r<=1 (lowering coverage =
  // randomly dropping reads -- exact; validated), Poisson up-scaling for r>1.
  static double dp_thin_cdf(int32_t c, double r, uint32_t k)
  {
    if (r <= 1.0) return (static_cast<double>(k) >= static_cast<double>(c)) ? 1.0 : bdtr(static_cast<double>(k), static_cast<double>(c), r);
    double mu = r * static_cast<double>(c);
    return (mu > 0.0) ? incompletegamma(static_cast<double>(k) + 1.0, mu, /*complemented=*/true) : 1.0;
  }

  // Probability a normal position on a sequence at relative coverage r (=avgcov/C_ref) is spanned by
  // <= k concordant pairs, by PROJECTING the normal reference distribution to that coverage.
  //  - Empirical (use_empirical): weight each reference crossing bin by hist_ref[c] over [lo,hi];
  //    P is floored at 0.5/total_normal (the resolution limit -> skew capped at ~log10(2*total_normal)).
  //  - Negative-binomial fallback (small reference, few positions): weight by the fitted NB pmf over a
  //    broad range that captures the tail, and DO NOT floor -- the parametric tail lets the skew
  //    extrapolate past log10(N). (Conservative: NB overshoots the crossing lower tail, lowering skew.)
  static double dp_crossing_cdf(const vector<double>& hist_ref, int32_t lo, int32_t hi, double r, uint32_t k,
                                bool use_empirical, double nb_size, double nb_mu)
  {
    if (!(r > 0.0)) return std::numeric_limits<double>::quiet_NaN();

    if (use_empirical) {
      double total = 0.0;
      for (int32_t c = lo; c <= hi && c < static_cast<int32_t>(hist_ref.size()); c++) total += hist_ref[c];
      if (!(total > 0.0)) return std::numeric_limits<double>::quiet_NaN();
      double P = 0.0;
      for (int32_t c = lo; c <= hi && c < static_cast<int32_t>(hist_ref.size()); c++) {
        if (hist_ref[c] <= 0.0) continue;
        P += hist_ref[c] * dp_thin_cdf(c, r, k);
      }
      double Pf = P / total;
      double floorP = 0.5 / total;   // resolution limit -> skew capped at ~log10(2*total_normal)
      return (Pf < floorP) ? floorP : Pf;
    }

    // Negative-binomial fallback: sum the NB pmf out to mean + 12 SD (essentially all mass), no floor.
    if (!(nb_mu > 0.0)) return std::numeric_limits<double>::quiet_NaN();
    double nb_var = nb_mu + (nb_size > 0.0 ? nb_mu * nb_mu / nb_size : 0.0);
    int32_t cmax = static_cast<int32_t>(ceil(nb_mu + 12.0 * sqrt(nb_var)));
    if (cmax < 1) cmax = 1;
    if (cmax > 1000000) cmax = 1000000;
    double wsum = 0.0, P = 0.0;
    for (int32_t c = 1; c <= cmax; c++) {
      double w = dp_nbinom_pmf(c, nb_size, nb_mu);
      if (w <= 0.0) continue;
      wsum += w;
      P += w * dp_thin_cdf(c, r, k);
    }
    return (wsum > 0.0) ? (P / wsum) : std::numeric_limits<double>::quiet_NaN();
  }

  // DP "skew" score: -log10 P(crossing <= k), projecting the normal reference to seq X's coverage.
  static double dp_discordance_skew(const vector<double>& hist_ref, int32_t lo, int32_t hi,
                                    double C_ref, double avgcov_X, uint32_t k,
                                    bool use_empirical, double nb_size, double nb_mu)
  {
    if (hist_ref.empty() || !(C_ref > 0.0) || !(avgcov_X > 0.0)) return std::numeric_limits<double>::quiet_NaN();
    double P = dp_crossing_cdf(hist_ref, lo, hi, avgcov_X / C_ref, k, use_empirical, nb_size, nb_mu);
    if (std::isnan(P)) return P;
    double sc = (P > 0.0) ? (-log10(P)) : kDPMaxScore;
    if (sc < 0.0) sc = 0.0;   // P==1 -> -log10 gives -0.0
    return (sc > kDPMaxScore) ? kDPMaxScore : sc;
  }

  // Poisson pmf P(X=v; mu).
  static double dp_poisson_pmf(uint32_t v, double mu)
  {
    if (mu <= 0.0) return (v == 0) ? 1.0 : 0.0;
    return exp(-mu + static_cast<double>(v) * log(mu) - lgamma(static_cast<double>(v) + 1.0));
  }

  // Projected crossing pmf for a sequence at relative coverage r, scaled to `scale` total counts, using
  // only the normal reference bins [clo,chi]: out[v] = scale * Σ_c (hist_ref[c]/total_normal) · PMF(v;c,r)
  // (binomial for r<=1, Poisson for r>1).
  static void dp_project_crossing(const vector<double>& hist_ref, int32_t clo, int32_t chi, double r,
                                  uint32_t maxv, double scale, vector<double>& out)
  {
    out.assign(maxv + 1, 0.0);
    double total = 0.0;
    for (int32_t c = clo; c <= chi && c < static_cast<int32_t>(hist_ref.size()); c++) total += hist_ref[c];
    if (!(total > 0.0) || !(r > 0.0)) return;
    for (int32_t c = clo; c <= chi && c < static_cast<int32_t>(hist_ref.size()); c++) {
      if (hist_ref[c] <= 0.0) continue;
      double w = hist_ref[c] / total;
      if (r <= 1.0) {
        uint32_t vhi = static_cast<uint32_t>(min<int32_t>(c, static_cast<int32_t>(maxv)));
        for (uint32_t v = 0; v <= vhi; v++) out[v] += w * binomial(r, c, static_cast<int32_t>(v));
      } else {
        double mu = r * static_cast<double>(c);
        for (uint32_t v = 0; v <= maxv; v++) out[v] += w * dp_poisson_pmf(v, mu);
      }
    }
    for (uint32_t v = 0; v <= maxv; v++) out[v] *= scale;
  }

  // Render a crossing-distribution SVG in the coverage-plot style: empirical histogram as points,
  // colored black inside the normal window [clo,chi] and red outside (censored deletions/repeats), plus
  // an optional projected line (blue). `emp`/`proj` are indexed by crossing value (counts).
  static void render_crossing_plot(const string& svg, const vector<double>& emp, const vector<double>& proj,
                                   int32_t clo, int32_t chi, const string& title, const string& xlabel,
                                   const string& emp_label, const string& proj_label)
  {
    // Peak / y-scale from the NORMAL window (so the crossing=0 deletion spike and repeats don't dominate
    // the axes and hide the bulk).
    uint32_t maxv = 0; double peak = 0.0, max_y = 0.0; uint32_t peak_i = static_cast<uint32_t>(max(1, clo));
    for (uint32_t i = 0; i < emp.size(); i++) if (emp[i] > 0) maxv = i;
    for (int32_t i = clo; i <= chi && i < static_cast<int32_t>(emp.size()); i++) {
      if (emp[i] > peak) { peak = emp[i]; peak_i = static_cast<uint32_t>(i); }
      if (emp[i] > max_y) max_y = emp[i];
    }
    for (uint32_t i = 0; i < proj.size(); i++) { if (proj[i] > max_y) max_y = proj[i]; if (proj[i] > 0) maxv = max(maxv, i); }
    if (peak <= 0.0) return;
    uint32_t graph_end = static_cast<uint32_t>(chi) + static_cast<uint32_t>((chi - clo) / 4 + 1);
    graph_end = min(graph_end, maxv + 1);
    if (graph_end < static_cast<uint32_t>(1.2 * peak_i)) graph_end = static_cast<uint32_t>(1.2 * peak_i);

    string emp_tab = svg + ".emp.tab";
    { ofstream o(emp_tab.c_str()); for (uint32_t i = 0; i < emp.size(); i++) o << i << "\t" << emp[i] << endl; }
    bool has_proj = !proj.empty();
    string proj_tab = svg + ".proj.tab";
    if (has_proj) { ofstream o(proj_tab.c_str()); for (uint32_t i = 0; i < proj.size(); i++) o << i << "\t" << proj[i] << endl; }

    ostringstream s;
    s << "set terminal svg size 1320,720 font ',16'" << endl;
    s << "set output " << double_quote(svg) << endl;
    s << "set tics out" << endl;
    s << "set border lw 2" << endl;
    s << "set title " << double_quote(title) << " font ',20'" << endl;
    s << "set xlabel " << double_quote(xlabel) << endl;
    s << "set ylabel 'Number of reference positions'" << endl;
    s << "set xrange [0:" << graph_end << "]" << endl;
    s << "set yrange [0:" << to_string(max_y * 1.05, 6) << "]" << endl;
    s << "set key top right font ',16' spacing 2" << endl;
    vector<string> cl;
    cl.push_back(double_quote(emp_tab) + " using ($1>=" + to_string(clo) + "&&$1<=" + to_string(chi) + "?$1:NaN):2 with points pt 6 lc rgb 'black' title " + double_quote(emp_label));
    cl.push_back(double_quote(emp_tab) + " using ($1<" + to_string(clo) + "||$1>" + to_string(chi) + "?$1:NaN):2 with points pt 6 lc rgb 'red' title 'censored'");
    if (has_proj) cl.push_back(double_quote(proj_tab) + " using 1:2 with lines lw 3 lc rgb 'blue' title " + double_quote(proj_label));
    s << "plot " << join(cl, string(", \\\n     ")) << endl;

    string base = svg + "." + to_string(getpid());
    run_gnuplot_script(s.str(), base + ".gp", base + ".gp.log");
    make_svg_responsive(svg);
    remove((base + ".gp.log").c_str());
    remove(emp_tab.c_str());
    if (has_proj) remove(proj_tab.c_str());
  }

  void draw_concordant_pair_crossing_plots(const Settings& settings, Summary& summary, cReferenceSequences& ref_seq_info)
  {
    vector<double> hist_ref; double C_ref = 0.0; string R; int32_t ref_lo = 0, ref_hi = 0;
    double ref_N = 0.0; bool ref_use_empirical = true; double ref_nb_size = 0.0, ref_nb_mu = 0.0;
    if (!dp_load_crossing_reference(settings, summary, ref_seq_info, hist_ref, C_ref, R, ref_lo, ref_hi,
                                    ref_N, ref_use_empirical, ref_nb_size, ref_nb_mu)) return;
    create_path(settings.evidence_path);

    // Run-wide reference distribution (the null actually used; from the longest sequence R).
    render_crossing_plot(settings.concordant_pair_crossing_plot_file_name, hist_ref, vector<double>(),
                         ref_lo, ref_hi, "Concordant Pair Crossing Distribution",
                         "Concordant pairs spanning a position",
                         "Crossing distribution (" + R + ")", "");

    // Per-seq overlays: each sequence's empirical distribution vs the reference projected to its coverage.
    for (cReferenceSequences::iterator it = ref_seq_info.begin(); it != ref_seq_info.end(); it++) {
      string seq_id = it->m_seq_id;
      vector<double> emp;
      string fn = Settings::file_name(settings.concordant_pair_crossing_distribution_file_name, "#", seq_id);
      if (!dp_read_crossing_hist(fn, emp)) continue;
      double avgcov = dp_seq_coverage(summary, seq_id);
      if (!(avgcov > 0.0) || !(C_ref > 0.0)) continue;
      uint32_t maxv = emp.empty() ? 0 : static_cast<uint32_t>(emp.size() - 1);
      int32_t seq_lo = 0, seq_hi = 0; dp_crossing_censor(emp, seq_lo, seq_hi);  // this seq's own normal window
      // Scale the projected pmf to the count of *normal* interior positions (the empirical mass inside the
      // seq's own censor window) -- NOT the grand total, which is dominated by the deletion spike at
      // crossing 0 and the censored tails. At r=1 (seq == reference) this makes projected overlay empirical.
      double total_norm = 0.0;
      for (int32_t v = seq_lo; v <= seq_hi && v < static_cast<int32_t>(emp.size()); v++) total_norm += emp[v];
      vector<double> proj;  // project the normal reference [ref_lo,ref_hi] to this seq's coverage
      dp_project_crossing(hist_ref, ref_lo, ref_hi, avgcov / C_ref, maxv, total_norm, proj);
      string svg = Settings::file_name(settings.concordant_pair_crossing_seq_plot_file_name, "#", seq_id);
      render_crossing_plot(svg, emp, proj, seq_lo, seq_hi, "Concordant Pair Crossing: " + seq_id,
                           "Concordant pairs spanning a position",
                           "empirical (" + seq_id + ")", "projected from " + R);
    }
  }

  string dp_crossing_model_description(const Settings& settings, Summary& summary, cReferenceSequences& ref_seq_info)
  {
    vector<double> hist_ref; double C_ref = 0.0; string R; int32_t ref_lo = 0, ref_hi = 0;
    double N = 0.0; bool use_empirical = true; double nb_size = 0.0, nb_mu = 0.0;
    if (!dp_load_crossing_reference(settings, summary, ref_seq_info, hist_ref, C_ref, R, ref_lo, ref_hi,
                                    N, use_empirical, nb_size, nb_mu)) return "";
    string npos = to_string(static_cast<uint64_t>(N)) + " positions in " + R;
    return use_empirical ? ("empirical (" + npos + ")")
                         : ("negative binomial (fallback; only " + npos + ")");
  }

  void predict_discordant_pairs(const Settings& settings, Summary& summary, cReferenceSequences& ref_seq_info)
  {
    cGenomeDiff dp_gd;

    // Run-wide reference crossing distribution + reference coverage + normal-bulk censor window, plus
    // the empirical-vs-negative-binomial decision (small references fall back to the parametric fit).
    vector<double> crossing_hist_ref; double crossing_C_ref = 0.0; string crossing_ref_seq_id;
    int32_t crossing_lo = 0, crossing_hi = 0;
    double crossing_N = 0.0; bool crossing_use_empirical = true; double crossing_nb_size = 0.0, crossing_nb_mu = 0.0;
    bool have_crossing = dp_load_crossing_reference(settings, summary, ref_seq_info, crossing_hist_ref, crossing_C_ref, crossing_ref_seq_id, crossing_lo, crossing_hi,
                                                    crossing_N, crossing_use_empirical, crossing_nb_size, crossing_nb_mu);

    //
    // Step 0: library orientation (inner3p) + rescan window (distance_cutoff). FF/RR unsupported.
    //
    bool inner3p = true;
    double distance_cutoff = 0.0;
    double pair_median = 0.0;
    if (!dp_library_params(summary, inner3p, distance_cutoff, pair_median, /*warn=*/true)) {
      dp_gd.write(settings.dp_genome_diff_file_name);
      return;
    }

    // Empirical-insert mixture model for the Bayes outlier test on breakpoint placement (optional: if the
    // persisted insert histogram is missing, placement falls back to the raw innermost read edge).
    dp_insert_model insert_model;
    bool have_insert = dp_load_insert_model(settings, summary, ref_seq_info, insert_model);

    //
    // Step 1: Re-read the candidate regions CSV.
    //  columns: seq_id,start,end,strand,orientation,length,max_discordant_count,discordant_pairs
    //  discordant_pairs is the final field: ';'-joined keys, and never contains a comma.
    //
    if (!file_exists(settings.dp_candidate_regions_file_name.c_str())) {
      dp_gd.write(settings.dp_genome_diff_file_name);
      return;
    }

    vector<dp_region_row> regions;
    {
      ifstream in(settings.dp_candidate_regions_file_name.c_str());
      string line;
      getline(in, line); // header
      while (getline(in, line)) {
        if (line.empty()) continue;
        vector<string> f = split(line, ",");
        if (f.size() < 7) continue;

        dp_region_row r;
        r.seq_id = f[0];
        r.start  = from_string<uint32_t>(f[1]);
        r.end    = from_string<uint32_t>(f[2]);
        r.strand = f[3].empty() ? 'F' : f[3][0];

        // The keys field is column index 7 (empty if the region had no descriptors).
        string key_field = (f.size() >= 8) ? f[7] : "";
        vector<string> keys = split(key_field, ";");
        for (size_t i = 0; i < keys.size(); i++) {
          if (!keys[i].empty()) r.keys.insert(keys[i]);
        }
        regions.push_back(r);
      }
    }

    //
    // Step 2: Build the graph. Edge weight between two regions = number of read pairs they share.
    //  Each pair key belongs to at most two distinct regions (a pair's two mates land in <= 2 regions).
    //
    map<string, set<int> > key_to_regions;
    for (size_t ri = 0; ri < regions.size(); ri++) {
      for (set<string>::const_iterator k = regions[ri].keys.begin(); k != regions[ri].keys.end(); k++) {
        key_to_regions[*k].insert(static_cast<int>(ri));
      }
    }
    map<pair<int, int>, int> edge_weight;
    for (map<string, set<int> >::const_iterator it = key_to_regions.begin(); it != key_to_regions.end(); it++) {
      if (it->second.size() == 2) {
        set<int>::const_iterator si = it->second.begin();
        int a = *si; ++si;
        int b = *si;
        edge_weight[make_pair(min(a, b), max(a, b))]++;
      }
    }

    //
    // Step 3: Collect candidate edges (weight = number of shared read pairs), highest weight first.
    //
    vector<pair<int, pair<int, int> > > edges; // (weight, (region_a, region_b))
    for (map<pair<int, int>, int>::const_iterator it = edge_weight.begin(); it != edge_weight.end(); it++) {
      edges.push_back(make_pair(it->second, it->first));
    }
    // Highest weight first; deterministic tiebreak by region indices.
    sort(edges.begin(), edges.end(),
         [](const pair<int, pair<int, int> >& x, const pair<int, pair<int, int> >& y) {
           if (x.first != y.first) return x.first > y.first;
           return x.second < y.second;
         });

    //
    // Step 4: Emit one DP item per edge with >= kDPMinSharedPairs shared read pairs.
    //  No one-to-one matching: a region whose discordant reads jump to several places (a rearrangement
    //  "hub") contributes one DP item per qualifying partner. A read pair's two mates may therefore
    //  appear in more than one DP item (counts across a hub's items can overlap).
    //
    // Optional per-side BAM rescan (fills concordant/unpaired counts, a refined discordant count, and
    // read-based position shifts). Skipped if the reference BAM/FASTA aren't available or D is unknown.
    dp_side_scanner* scanner = NULL;
    if (distance_cutoff > 0.0
        && file_exists(settings.reference_bam_file_name.c_str())
        && file_exists(settings.reference_fasta_file_name.c_str())) {
      scanner = new dp_side_scanner(settings.reference_bam_file_name, settings.reference_fasta_file_name);
    }

    for (size_t e = 0; e < edges.size(); e++) {
      int weight = edges[e].first;
      if (weight < kDPMinSharedPairs) break; // edges are sorted descending; nothing left qualifies
      int a = edges[e].second.first;
      int b = edges[e].second.second;

      int32_t pos_a, strand_a, pos_b, strand_b;
      dp_region_to_side(regions[a], inner3p, pos_a, strand_a);
      dp_region_to_side(regions[b], inner3p, pos_b, strand_b);

      // side_1 = the side with the lower (seq_id, position).
      bool a_is_side_1;
      if (regions[a].seq_id != regions[b].seq_id)
        a_is_side_1 = (regions[a].seq_id < regions[b].seq_id);
      else
        a_is_side_1 = (pos_a <= pos_b);

      string  s1_seq_id, s2_seq_id;
      int32_t s1_pos, s1_strand, s2_pos, s2_strand;
      if (a_is_side_1) {
        s1_seq_id = regions[a].seq_id; s1_pos = pos_a; s1_strand = strand_a;
        s2_seq_id = regions[b].seq_id; s2_pos = pos_b; s2_strand = strand_b;
      } else {
        s1_seq_id = regions[b].seq_id; s1_pos = pos_b; s1_strand = strand_b;
        s2_seq_id = regions[a].seq_id; s2_pos = pos_a; s2_strand = strand_a;
      }

      // Initial (region-derived) breakpoint estimate. The overlapping-mate exclusion in every
      // placement/gathering pass below references these stable positions, not the re-anchored window
      // position it scans at (the final count references the placed positions instead).
      int32_t init1 = s1_pos, init2 = s2_pos;

      // Each side's crossing-read strand (same as the region strand that produced it).
      bool s1_fwd = (inner3p == (s1_strand == -1));
      bool s2_fwd = (inner3p == (s2_strand == -1));
      int32_t s1_tid = scanner ? scanner->tid_for_seq_id(s1_seq_id) : -1;
      int32_t s2_tid = scanner ? scanner->tid_for_seq_id(s2_seq_id) : -1;

      // Place each breakpoint at the supporting reads' median outside end, shifted toward the junction by
      // half the median pair distance (mh1/mh2), then CLAMP the reported/displayed coordinate outward to
      // the furthest supporting read's junction-facing edge so no supporting read straddles the seam (max
      // for s=-1, min for s=+1). Classification (the count below and the plot) uses the UN-clamped mh: the
      // overlap test needs a consistent junction reference for BOTH sides -- clamping each side to its own
      // read edge would make every read's inner gap non-negative and defeat overlap detection.
      int32_t mh1 = s1_pos, mh2 = s2_pos;
      int32_t edge1 = 0, edge2 = 0; bool h1 = false, h2 = false;
      if (scanner && pair_median > 0.0) {
        int32_t H = static_cast<int32_t>(pair_median / 2.0 + 0.5);
        // Side 1: gather with the region other-side position; its shifted median (mh1) then becomes the
        // reference the side-2 gather uses, so both sides' overlap tests share one consistent junction.
        int32_t med1 = 0, med2 = 0, outer1 = 0, outer2 = 0; size_t c1 = 0, c2 = 0;
        h1 = scanner->supporting_outer_median(s1_seq_id, s1_pos, s1_strand, s1_fwd, s2_tid, s2_pos, s2_fwd, distance_cutoff, init1, init2, med1, edge1, outer1, c1);
        if (h1) mh1 = med1 + (s1_strand == -1 ? +H : -H);
        h2 = scanner->supporting_outer_median(s2_seq_id, s2_pos, s2_strand, s2_fwd, s1_tid, mh1, s1_fwd, distance_cutoff, init2, init1, med2, edge2, outer2, c2);
        if (h2) mh2 = med2 + (s2_strand == -1 ? +H : -H);
        // Placement strategy. Conservative (default): put each coordinate at its side's innermost aligned
        // read edge (edge1/edge2) -- aligned reads cannot extend past the breakpoint, so this never
        // overshoots it. The median+half path (max/min of mh and the edge) is retained but disabled -- we
        // will iterate on placement / confidence limits next. Either way, the counting below classifies at
        // the un-clamped mh1/mh2 so the overlap test keeps a consistent junction reference for both sides.
        const bool conservative_edge_placement = true;
        if (h1) {
          int32_t p1 = conservative_edge_placement ? edge1
                     : ((s1_strand == -1) ? max(mh1, edge1) : min(mh1, edge1));
          s1_pos = max(1, min(scanner->seq_length(s1_tid), p1));
        }
        if (h2) {
          int32_t p2 = conservative_edge_placement ? edge2
                     : ((s2_strand == -1) ? max(mh2, edge2) : min(mh2, edge2));
          s2_pos = max(1, min(scanner->seq_length(s2_tid), p2));
        }

        // Refinement. Re-anchor the wide (+/-D) window so the OUTERMOST supporting read sits on its outer
        // boundary -- a [p-D,p] window at p = outer+D is exactly [outer, outer+D] -- and re-scan for
        // discordant pairs whose reads fall in the re-anchored windows on both sides. That window is a
        // superset of the initial one, so the supporting count can only rise; a drop means a logic error.
        // If reads were added, move each coordinate out to the new innermost aligned read edge.
        if (h1 && h2) {
          int32_t Di = static_cast<int32_t>(distance_cutoff);
          int32_t p1r = outer1 + (s1_strand == -1 ? +Di : -Di);
          int32_t p2r = outer2 + (s2_strand == -1 ? +Di : -Di);
          int32_t rmed = 0, router = 0, e1r = edge1, e2r = edge2; size_t rc1 = 0, rc2 = 0;
          bool r1 = scanner->supporting_outer_median(s1_seq_id, p1r, s1_strand, s1_fwd, s2_tid, p2r, s2_fwd, distance_cutoff, init1, init2, rmed, e1r, router, rc1);
          bool r2 = scanner->supporting_outer_median(s2_seq_id, p2r, s2_strand, s2_fwd, s1_tid, p1r, s1_fwd, distance_cutoff, init2, init1, rmed, e2r, router, rc2);
          if (rc1 < c1 || rc2 < c2)
            WARN("DP refine: supporting count dropped at " + s1_seq_id + ":" + to_string(s1_pos) +
                 " (side1 " + to_string(c1) + "->" + to_string(rc1) + ", side2 " + to_string(c2) + "->" + to_string(rc2) + ")");
          // Move each coordinate only when the re-scan added reads, and only TOWARD the junction: take the
          // more-inward of the initial and refined aligned edges (aligned edges never pass the breakpoint,
          // so this can only extend to the true edge, never overshoot or retreat).
          if (r1 && rc1 >= c1)
            s1_pos = max(1, min(scanner->seq_length(s1_tid), s1_strand == -1 ? max(edge1, e1r) : min(edge1, e1r)));
          if (r2 && rc2 >= c2)
            s2_pos = max(1, min(scanner->seq_length(s2_tid), s2_strand == -1 ? max(edge2, e2r) : min(edge2, e2r)));
        }
      }

      // Bayes outlier test on placement: gather the supporting pairs over a symmetric window and move each
      // coordinate to the innermost read edge that is NOT a one-off insert-size outlier -- a lone read
      // whose own inferred insert is anomalous (BF < 1/3) can't drag side_x_position off the cluster.
      if (have_insert && scanner) {
        int32_t readlen = static_cast<int32_t>(summary.sequence_conversion.read_length_avg + 0.5);
        vector<dp_pair_ends> pr;
        if (scanner->gather_pairs(s1_seq_id, s1_pos, s1_strand, s1_fwd, s2_tid, s2_pos, s2_fwd, distance_cutoff, init1, init2, s2_strand, readlen, pr)
            && pr.size() >= 2) {
          int32_t r1 = dp_robust_edge(pr, /*this_is_side1=*/true,  s1_strand, s2_pos, s2_strand, insert_model);
          int32_t r2 = dp_robust_edge(pr, /*this_is_side1=*/false, s2_strand, r1,     s1_strand, insert_model);
          s1_pos = max(1, min(scanner->seq_length(s1_tid), r1));
          s2_pos = max(1, min(scanner->seq_length(s2_tid), r2));
        }
      }

      cDiffEntry dp(DP);
      dp[SIDE_1_SEQ_ID]  = s1_seq_id;
      dp[SIDE_1_POSITION] = to_string(s1_pos);
      dp[SIDE_1_STRAND]  = to_string(s1_strand);
      dp[SIDE_2_SEQ_ID]  = s2_seq_id;
      dp[SIDE_2_POSITION] = to_string(s2_pos);
      dp[SIDE_2_STRAND]  = to_string(s2_strand);

      // Heuristic count from region overlap; kept so we can see the rescan hold steady or increase.
      dp["candidate_discordant_count"] = to_string(weight);

      int k_support = weight;
      if (scanner) {
        // Count the three read categories at each side, classifying at the FINAL placed positions
        // (s1_pos/s2_pos) so the counts describe the reported breakpoint. The overlapping-mate
        // exclusion is referenced to those same placed positions (they are our best breakpoint
        // estimate now that placement is done).
        scanner->scan(s1_seq_id, s1_pos, s1_strand, s1_fwd, s2_tid, s2_pos, s2_fwd, distance_cutoff, s1_pos, s2_pos);
        int c1a = scanner->supporting(), c2a = scanner->concordant(), c3a = scanner->unpaired();
        set<string> support_nums_1 = scanner->supporting_nums();   // copy before the next scan overwrites

        scanner->scan(s2_seq_id, s2_pos, s2_strand, s2_fwd, s1_tid, s1_pos, s1_fwd, distance_cutoff, s2_pos, s1_pos);
        int c1b = scanner->supporting(), c2b = scanner->concordant(), c3b = scanner->unpaired();
        const set<string>& support_nums_2 = scanner->supporting_nums();

        // True support = read pairs whose BOTH mates qualify -- one at each side (intersect the sides'
        // supporting read-pair numbers). A single-side per-side count (c1a/c1b) over-counts at a
        // breakpoint SHARED with a neighboring junction: those reads' mates fall inside this junction's
        // wide (+/-D) partner window but land at the neighbor's breakpoint, so they appear on only one
        // side and drop out of the intersection. This matches exactly the pairs drawn in the joined plot.
        k_support = 0;
        for (set<string>::const_iterator n = support_nums_1.begin(); n != support_nums_1.end(); n++)
          if (support_nums_2.count(*n)) k_support++;

        dp["discordant_count"] = to_string(k_support);
        dp["side_1_discordant_count"] = to_string(c1a);   // raw per-side counts (may exceed the paired
        dp["side_2_discordant_count"] = to_string(c1b);   // count when a side is a shared-breakpoint hub)
        dp["side_1_concordant_count"] = to_string(c2a);
        dp["side_2_concordant_count"] = to_string(c2b);
        dp["side_1_unpaired_count"] = to_string(c3a);
        dp["side_2_unpaired_count"] = to_string(c3b);
      } else {
        // No BAM available: fall back to the heuristic count.
        dp["discordant_count"] = to_string(weight);
      }

      // Discordance "skew" score: -log10 P(a normal position on side_1's seq_id is spanned by <= k
      // concordant pairs), projecting the run-wide reference distribution to this seq_id's coverage.
      double dp_score = std::numeric_limits<double>::quiet_NaN();
      if (have_crossing) {
        double avgcov = dp_seq_coverage(summary, s1_seq_id);
        dp_score = dp_discordance_skew(crossing_hist_ref, crossing_lo, crossing_hi, crossing_C_ref, avgcov, static_cast<uint32_t>(k_support < 0 ? 0 : k_support),
                                       crossing_use_empirical, crossing_nb_size, crossing_nb_mu);
      }
      dp[NEG_LOG10_DISCORDANCE_P_VALUE] = std::isnan(dp_score) ? string("NT") : to_string(dp_score, 1, false);

      // Reject weakly-supported junctions: a skew above the cutoff means the discordant support (k) is
      // anomalously low relative to the concordant crossing expected at this coverage -- i.e. too few
      // discordant pairs where many concordant pairs still span the locus, so it reads as a contiguous
      // (non-junction) position. Rejected DP items are moved to marginal.html (0 = cutoff off).
      if (!std::isnan(dp_score) && settings.discordant_pair_skew_cutoff > 0.0
          && dp_score > settings.discordant_pair_skew_cutoff) {
        dp.add_reject_reason("CONCORDANT_PAIR_SKEW");
      }

      dp_gd.add(dp);
    }

    if (scanner) delete scanner;

    dp_gd.write(settings.dp_genome_diff_file_name);
  }

  // -------------------------------------------------------------------------------------------------
  // Per-side read-pair plots
  // -------------------------------------------------------------------------------------------------

  // Format an integer with thousands separators, e.g. 1234560 -> "1,234,560".
  static string dp_commafy(int64_t v)
  {
    bool neg = v < 0; if (neg) v = -v;
    string digits = to_string(v);
    string out;
    int c = 0;
    for (string::reverse_iterator it = digits.rbegin(); it != digits.rend(); ++it) {
      if (c && c % 3 == 0) out.push_back(',');
      out.push_back(*it);
      c++;
    }
    if (neg) out.push_back('-');
    reverse(out.begin(), out.end());
    return out;
  }

  // A "nice" tick interval (a multiple of 10) giving roughly ~10 labels across a span.
  static int64_t dp_nice_tick(int64_t span)
  {
    double raw = span / 10.0;
    static const int64_t nice[] = {10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000};
    for (size_t i = 0; i < sizeof(nice)/sizeof(nice[0]); i++)
      if (nice[i] >= raw) return nice[i];
    return 100000;
  }

  // Build a gnuplot explicit-tics list ("label" pos, ...) at multiples of `step` within [lo, hi],
  // with comma-formatted whole-number labels.
  static string dp_xtics_list(int64_t lo, int64_t hi, int64_t step)
  {
    string out;
    int64_t first = ((lo + step - 1) / step) * step;  // first multiple of step >= lo
    bool comma = false;
    for (int64_t t = first; t <= hi; t += step) {
      if (comma) out += ", ";
      out += double_quote(dp_commafy(t)) + " " + to_string(t);
      comma = true;
    }
    return out;
  }

  // Render one per-side plot (SVG) via gnuplot: each counted read pair on its own lane, reads drawn as
  // arrows (in mapping direction), connected by a dashed line (mate for concordant; to the side
  // position p for discordant). The x window is symmetric about p and sized ~1.1x the furthest read
  // extent from p. `strand` is the side strand: discordant reads sit on the kept flank (coords <= p for
  // s=-1, >= p for s=+1) and never cross p, so the opposite half of the plot is shaded light gray.
  // `mate_prefix` maps a read-name file prefix to its mate's, to label concordant lanes with both names.
  // `kept_color` shades the kept (read-bearing) half of the plot (light green for side 1, light yellow
  // for side 2 -- matching the joined '*' plot).
  static const string DP_SIDE1_COLOR = "'#e6ffe6'";  // light green
  static const string DP_SIDE2_COLOR = "'#ffffcc'";  // light yellow
  static void render_dp_side_plot(const string& output_svg, vector<dp_draw_read> reads,
                                  int32_t p, int32_t strand, const map<int,int>& mate_prefix,
                                  const string& kept_color)
  {
    // Lay pairs out by leftmost read coordinate for a tidy pileup.
    sort(reads.begin(), reads.end(), [](const dp_draw_read& a, const dp_draw_read& b) {
      if (a.read_start != b.read_start) return a.read_start < b.read_start;
      return a.read_end < b.read_end;
    });

    int n = static_cast<int>(reads.size());

    // Symmetric window about p, sized ~1.1x the furthest extent of any drawn read (or mate) from p.
    int64_t lo = p, hi = p;
    for (int i = 0; i < n; i++) {
      const dp_draw_read& r = reads[i];
      lo = min(lo, (int64_t)r.read_start); hi = max(hi, (int64_t)r.read_end);
      if (r.category == 2) {
        int64_t mlen = r.read_end - r.read_start;
        lo = min(lo, (int64_t)r.mate_start); hi = max(hi, (int64_t)r.mate_start + mlen);
      }
    }
    int64_t reach = max((int64_t)(p - lo), (int64_t)(hi - p));
    if (reach < 10) reach = 10;
    int64_t half = (int64_t)(1.1 * reach + 0.5);
    int64_t xmin = p - half;
    int64_t xmax = p + half;

    const string CONC = "'#1f77b4'";  // concordant = blue
    const string DISC = "'#d62728'";  // discordant = red

    string f_cr = output_svg + ".cr.tab";  // concordant read arrows (x y dx dy)
    string f_dr = output_svg + ".dr.tab";  // discordant read arrows
    string f_cc = output_svg + ".cc.tab";  // concordant connectors (2-point segments)
    string f_dc = output_svg + ".dc.tab";  // discordant connectors
    ofstream cr(f_cr.c_str()), dr(f_dr.c_str()), cc(f_cc.c_str()), dc(f_dc.c_str());
    bool has_cr=false, has_dr=false, has_cc=false, has_dc=false;

    for (int i = 0; i < n; i++) {
      const dp_draw_read& r = reads[i];
      int y = i + 1;

      // The read as an arrow tail->head in its mapping direction, drawn to its ALIGNED extent only
      // (soft-clipped bases are excluded -- they map to the other side of the junction).
      int32_t tail = r.read_reversed ? r.read_end   : r.read_start;
      int32_t head = r.read_reversed ? r.read_start : r.read_end;
      if (r.category == 2) { cr << tail << "\t" << y << "\t" << (head - tail) << "\t0\n"; has_cr = true; }
      else                 { dr << tail << "\t" << y << "\t" << (head - tail) << "\t0\n"; has_dr = true; }

      if (r.category == 2) {
        // Mate arrow (length approximated by this read's aligned length; mate CIGAR isn't available).
        int32_t mlen   = r.read_end - r.read_start;
        int32_t m_tail = r.mate_reversed ? (r.mate_start + mlen) : r.mate_start;
        int32_t m_head = r.mate_reversed ? r.mate_start : (r.mate_start + mlen);
        cr << m_tail << "\t" << y << "\t" << (m_head - m_tail) << "\t0\n"; has_cr = true;
        // Dashed connector between the two reads (anchor end -> mate start).
        cc << r.anchor << "\t" << y << "\n" << r.mate_start << "\t" << y << "\n\n"; has_cc = true;
      } else {
        // Discordant: dashed connector from the read's junction-facing anchor to p (plot center).
        dc << r.anchor << "\t" << y << "\n" << p << "\t" << y << "\n\n"; has_dc = true;
      }
    }
    cr.close(); dr.close(); cc.close(); dc.close();

    ostringstream s;
    int height = max(300, n * 8 + 120);
    s << "set terminal svg size 1400," << height << " font ',11' noenhanced" << endl;
    s << "set output " << double_quote(output_svg) << endl;
    s << "unset title" << endl;
    s << "set xlabel 'Reference position (bp)' font ',22'" << endl;
    s << "set xrange [" << xmin << ":" << xmax << "]" << endl;
    s << "set x2range [" << xmin << ":" << xmax << "]" << endl;
    s << "set yrange [0:" << (n + 1) << "]" << endl;
    s << "set y2range [0:" << (n + 1) << "]" << endl;
    s << "set border 15 lw 1" << endl;   // full box around the plot
    s << "set tics out" << endl;

    // Comma-formatted, whole-number x ticks pinned to round positions, labeled on BOTH bottom and top.
    int64_t step = dp_nice_tick(xmax - xmin);
    string tics = dp_xtics_list(xmin, xmax, step);
    s << "set xtics (" << tics << ") nomirror font ',16'" << endl;
    s << "set x2tics (" << tics << ") font ',16'" << endl;

    // Read-pair names as lane labels down the right side. For a concordant pair (both mates drawn) the
    // label shows both read names; a discordant lane shows the single in-window read.
    s << "unset ytics" << endl;
    {
      ostringstream y2;
      for (int i = 0; i < n; i++) {
        if (i) y2 << ", ";
        string label = reads[i].name;
        if (reads[i].category == 2) {
          size_t colon = label.find(':');
          if (colon != string::npos) {
            int prefix = atoi(label.substr(0, colon).c_str());
            map<int,int>::const_iterator mp = mate_prefix.find(prefix);
            if (mp != mate_prefix.end())
              label += " / " + to_string(mp->second) + ":" + label.substr(colon + 1);
          }
        }
        y2 << double_quote(label) << " " << (i + 1);
      }
      if (n) s << "set y2tics (" << y2.str() << ") font ',9'" << endl;
      else   s << "unset y2tics" << endl;
    }

    // Shade the two halves: light gray over the half the discordant reads never extend into (the
    // non-kept side of p), and kept_color over the kept half where the reads sit.
    int64_t gray_lo, gray_hi, keep_lo, keep_hi;
    if (strand == -1) { gray_lo = p; gray_hi = xmax; keep_lo = xmin; keep_hi = p; }
    else              { gray_lo = xmin; gray_hi = p; keep_lo = p; keep_hi = xmax; }
    s << "set object 1 rectangle from " << gray_lo << ",0 to " << gray_hi << "," << (n + 1)
      << " fc rgb 'gray90' fs solid noborder behind" << endl;
    s << "set object 2 rectangle from " << keep_lo << ",0 to " << keep_hi << "," << (n + 1)
      << " fc rgb " << kept_color << " fs solid noborder behind" << endl;

    // Vertical marker at the side position (plot center).
    s << "set arrow from " << p << ",0 to " << p << "," << (n + 1) << " nohead lc rgb 'gray50' lw 1 dt 3 back" << endl;

    vector<string> clauses;  // connectors first so the read arrows draw on top
    if (has_cc) clauses.push_back(double_quote(f_cc) + " using 1:2 with lines dt 2 lc rgb " + CONC + " lw 1 notitle");
    if (has_dc) clauses.push_back(double_quote(f_dc) + " using 1:2 with lines dt 2 lc rgb " + DISC + " lw 1 notitle");
    if (has_cr) clauses.push_back(double_quote(f_cr) + " using 1:2:3:4 with vectors head filled size screen 0.008,20,60 lc rgb " + CONC + " lw 2 notitle");
    if (has_dr) clauses.push_back(double_quote(f_dr) + " using 1:2:3:4 with vectors head filled size screen 0.008,20,60 lc rgb " + DISC + " lw 2 notitle");
    if (clauses.empty()) clauses.push_back("-1 notitle");
    s << "plot " << join(clauses, string(", \\\n     ")) << endl;

    string script_name = output_svg + ".gp";
    string log_name    = output_svg + ".gp.log";
    run_gnuplot_script(s.str(), script_name, log_name);
    remove(log_name.c_str());
    remove(f_cr.c_str()); remove(f_dr.c_str()); remove(f_cc.c_str()); remove(f_dc.c_str());
  }

  // One discordant read pair, both mates drawn on the joined '*' plot.
  struct dp_joined_pair {
    int32_t s1_start, s1_end, s1_anchor; bool s1_rev;   // side_1 read (genomic)
    int32_t s2_start, s2_end, s2_anchor; bool s2_rev;   // side_2 read (genomic)
    string  s1_name, s2_name;
  };

  // Render the joined '*' plot. The x axis is a custom coordinate that stitches the two junction flanks
  // together at the center: side_1 fills the left half (its kept flank running out to the left, reversed
  // if s1=+1) up to side_1_position at x=0; side_2 fills the right half starting at side_2_position at
  // x=0 (reversed if s2=-1). Each discordant pair is drawn on its own lane -- the side_1 read on the
  // left, its side_2 mate on the right, both pointing toward the center and joined by a dashed line --
  // so a real junction reads as a clean set of pairs bridging the seam. Left half shaded light green,
  // right half light yellow.
  static void render_dp_joined_plot(const string& output_svg, vector<dp_joined_pair> pairs,
                                    int32_t p1, int32_t s1_str, int32_t p2, int32_t s2_str)
  {
    // Genomic -> joined-x transforms: side_1 maps to <=0 (side_1_position at x=0), side_2 to >=0
    // (side_2_position at x=0). The two positions are stitched together at the seam.
    auto j1 = [&](int32_t g) -> int64_t { return (s1_str == -1) ? (int64_t)(g - p1) : (int64_t)(p1 - g); };
    auto j2 = [&](int32_t g) -> int64_t { return (s2_str == -1) ? (int64_t)(p2 - g) : (int64_t)(g - p2); };

    // Lay pairs out by their left (side_1) extent.
    sort(pairs.begin(), pairs.end(), [&](const dp_joined_pair& a, const dp_joined_pair& b) {
      return min(j1(a.s1_start), j1(a.s1_end)) < min(j1(b.s1_start), j1(b.s1_end));
    });

    int n = static_cast<int>(pairs.size());

    // Symmetric window about the seam (x=0), sized ~1.1x the furthest read extent.
    int64_t reach = 10;
    for (int i = 0; i < n; i++) {
      const dp_joined_pair& q = pairs[i];
      reach = max(reach, -min(j1(q.s1_start), j1(q.s1_end)));
      reach = max(reach,  max(j2(q.s2_start), j2(q.s2_end)));
    }
    int64_t half = (int64_t)(1.1 * reach + 0.5);
    int64_t xmin = -half, xmax = half;

    const string DISC = "'#d62728'";  // discordant = red
    string f_r = output_svg + ".r.tab";   // read arrows (x y dx dy)
    string f_c = output_svg + ".c.tab";   // connectors
    ofstream fr(f_r.c_str()), fc(f_c.c_str());
    bool has_r=false, has_c=false;

    for (int i = 0; i < n; i++) {
      const dp_joined_pair& q = pairs[i];
      int y = i + 1;
      // side_1 read arrow (tail->head in mapping direction), drawn to its ALIGNED extent only (soft-clips
      // excluded -- they map across the junction to the other side).
      int64_t a_tail = j1(q.s1_rev ? q.s1_end : q.s1_start);
      int64_t a_head = j1(q.s1_rev ? q.s1_start : q.s1_end);
      fr << a_tail << "\t" << y << "\t" << (a_head - a_tail) << "\t0\n"; has_r = true;
      // side_2 mate arrow.
      int64_t b_tail = j2(q.s2_rev ? q.s2_end : q.s2_start);
      int64_t b_head = j2(q.s2_rev ? q.s2_start : q.s2_end);
      fr << b_tail << "\t" << y << "\t" << (b_head - b_tail) << "\t0\n"; has_r = true;
      // dashed connector between the two reads' junction-facing anchors, across the seam.
      fc << j1(q.s1_anchor) << "\t" << y << "\n" << j2(q.s2_anchor) << "\t" << y << "\n\n"; has_c = true;
    }
    fr.close(); fc.close();

    ostringstream s;
    int height = max(360, n * 8 + 190);
    s << "set terminal svg size 1400," << height << " font ',11' noenhanced" << endl;
    s << "set output " << double_quote(output_svg) << endl;
    s << "unset title" << endl;
    // Reserve room for the two-row tick labels on top and bottom (and the pushed-down xlabel).
    s << "set tmargin 5" << endl;
    s << "set bmargin 7" << endl;
    s << "set xlabel 'Joined reference position (bp)' font ',22' offset 0,-3" << endl;
    s << "set xrange [" << xmin << ":" << xmax << "]" << endl;
    s << "set x2range [" << xmin << ":" << xmax << "]" << endl;
    s << "set yrange [0:" << (n + 1) << "]" << endl;
    s << "set y2range [0:" << (n + 1) << "]" << endl;
    s << "set border 15 lw 1" << endl;
    s << "set tics out" << endl;

    // x ticks: label with each side's genomic coordinate (comma-formatted, reversed where the side is),
    // as TWO-ROW labels -- side_1 (left) coords in the upper row, side_2 (right) coords in the lower row
    // (via an embedded newline). The center tick at x=0 carries both: side_1_position above,
    // side_2_position below (consecutive genomic positions sharing one tick at this scale). The same
    // two-row labels are placed on both the top and bottom axes.
    {
      int64_t step = dp_nice_tick(xmax - xmin);
      ostringstream tks;
      bool first = true;
      auto emit = [&](int64_t x, const string& label) {
        if (!first) tks << ", "; first = false;
        tks << double_quote(label) << " " << x;
      };
      // Seam: side_1_position (upper) over side_2_position (lower).
      emit(0, dp_commafy(p1) + "\\n" + dp_commafy(p2));
      // Left half = side_1, upper row (trailing newline leaves the lower row blank).
      for (int64_t x = -step; x >= xmin; x -= step) {
        int64_t g = (s1_str == -1) ? (p1 + x) : (p1 - x);
        emit(x, dp_commafy(g) + "\\n");
      }
      // Right half = side_2, lower row (leading newline leaves the upper row blank).
      for (int64_t x = step; x <= xmax; x += step) {
        int64_t g = (s2_str == -1) ? (p2 - x) : (p2 + x);
        emit(x, "\\n" + dp_commafy(g));
      }
      s << "set xtics ("  << tks.str() << ") nomirror font ',16'" << endl;  // bottom axis
      s << "set x2tics (" << tks.str() << ") font ',16'" << endl;           // top axis (same two rows)
    }

    // Read-pair names down the right side (both mates).
    s << "unset ytics" << endl;
    if (n) {
      ostringstream y2;
      for (int i = 0; i < n; i++) {
        if (i) y2 << ", ";
        y2 << double_quote(pairs[i].s1_name + " / " + pairs[i].s2_name) << " " << (i + 1);
      }
      s << "set y2tics (" << y2.str() << ") font ',9'" << endl;
    } else {
      s << "unset y2tics" << endl;
    }

    // Shade the two flanks: side_1 (left) light green up to x=0, side_2 (right) light yellow from x=0.
    s << "set object 1 rectangle from " << xmin << ",0 to 0," << (n + 1)
      << " fc rgb " << DP_SIDE1_COLOR << " fs solid noborder behind" << endl;
    s << "set object 2 rectangle from 0,0 to " << xmax << "," << (n + 1)
      << " fc rgb " << DP_SIDE2_COLOR << " fs solid noborder behind" << endl;
    // Seam marker at x=0.
    s << "set arrow from 0,0 to 0," << (n + 1) << " nohead lc rgb 'gray50' lw 1 dt 3 back" << endl;

    vector<string> clauses;  // connectors first so the read arrows draw on top
    if (has_c) clauses.push_back(double_quote(f_c) + " using 1:2 with lines dt 2 lc rgb " + DISC + " lw 1 notitle");
    if (has_r) clauses.push_back(double_quote(f_r) + " using 1:2:3:4 with vectors head filled size screen 0.008,20,60 lc rgb " + DISC + " lw 2 notitle");
    if (clauses.empty()) clauses.push_back("-1 notitle");
    s << "plot " << join(clauses, string(", \\\n     ")) << endl;

    string script_name = output_svg + ".gp";
    string log_name    = output_svg + ".gp.log";
    run_gnuplot_script(s.str(), script_name, log_name);
    remove(log_name.c_str());
    remove(f_r.c_str()); remove(f_c.c_str());
  }

  // Extract the read-number (the part after ':') from a QNAME "<prefix>:<num>".
  // If v holds more than max_display entries, randomly down-sample it to max_display (seeded so runs are
  // reproducible) and return the ORIGINAL size; otherwise leave v alone and return 0. Used to bound the
  // number of read/pair lanes drawn in a DP plot (0 = show all).
  template <typename T>
  static size_t dp_cap_for_display(vector<T>& v, uint32_t max_display, uint32_t seed)
  {
    size_t total = v.size();
    if (max_display == 0 || total <= max_display) return 0;
    std::mt19937 rng(seed);
    std::shuffle(v.begin(), v.end(), rng);
    v.resize(max_display);
    return total;
  }

  static string dp_capped_message(uint32_t shown, size_t total)
  {
    return "Only " + to_string(shown) + " of " + to_string(total) + " mapped read pairs displayed.";
  }

  void draw_discordant_pair_evidence_plots(const Settings& settings, Summary& summary,
                                           cReferenceSequences& ref_seq_info, cGenomeDiff& gd)
  {
    (void)ref_seq_info;

    bool inner3p = true;
    double D = 0.0;
    double pair_median = 0.0;
    if (!dp_library_params(summary, inner3p, D, pair_median, /*warn=*/false)) return;  // predict already warned
    if (D <= 0.0) return;
    if (!file_exists(settings.reference_bam_file_name.c_str()) ||
        !file_exists(settings.reference_fasta_file_name.c_str())) return;

    diff_entry_list_t dp_list = gd.get_list(make_vector<gd_entry_type>(DP));
    if (dp_list.empty()) return;

    create_path(settings.evidence_path);
    dp_side_plot_gatherer g(settings.reference_bam_file_name, settings.reference_fasta_file_name);

    // Map a read-name file prefix (m_id + 1, as written by identify_mutations) to its mate's, so a
    // concordant lane can be labeled with both read names.
    map<int,int> mate_prefix;
    for (cReadFileSets::const_iterator rfs = settings.read_file_sets.begin(); rfs != settings.read_file_sets.end(); rfs++) {
      if (rfs->is_paired()) {
        int a = static_cast<int>(rfs->m_files[0].m_id) + 1;
        int b = static_cast<int>(rfs->m_files[1].m_id) + 1;
        mate_prefix[a] = b; mate_prefix[b] = a;
      }
    }

    for (diff_entry_list_t::iterator it = dp_list.begin(); it != dp_list.end(); it++) {
      cDiffEntry& dp = **it;

      string  s1_seq = dp[SIDE_1_SEQ_ID];
      int32_t s1_pos = from_string<int32_t>(dp[SIDE_1_POSITION]);
      int32_t s1_str = from_string<int32_t>(dp[SIDE_1_STRAND]);
      string  s2_seq = dp[SIDE_2_SEQ_ID];
      int32_t s2_pos = from_string<int32_t>(dp[SIDE_2_POSITION]);
      int32_t s2_str = from_string<int32_t>(dp[SIDE_2_STRAND]);

      bool s1_fwd = (inner3p == (s1_str == -1));
      bool s2_fwd = (inner3p == (s2_str == -1));
      int32_t s1_tid = g.tid_for_seq_id(s1_seq);
      int32_t s2_tid = g.tid_for_seq_id(s2_seq);

      uint32_t max_display = settings.max_displayed_reads;

      // Side 1 -- render the per-side plot and keep the (FULL) discordant reads for the joined plot.
      g.gather(s1_seq, s1_pos, s1_str, s1_fwd, s2_tid, s2_pos, s2_fwd, D);
      vector<dp_draw_read> s1_disc;
      for (vector<dp_draw_read>::const_iterator r = g.reads().begin(); r != g.reads().end(); r++)
        if (r->category == 1) s1_disc.push_back(*r);
      {
        // Cap a COPY for the plot (concordant + discordant lanes); the joined plot below keeps s1_disc full.
        vector<dp_draw_read> plot_reads = g.reads();
        size_t total = dp_cap_for_display(plot_reads, max_display, static_cast<uint32_t>(s1_pos));
        if (total) dp["_side_1_dp_plot_message"] = dp_capped_message(max_display, total);
        string svg = settings.evidence_path + "/DP_SIDE_1_" + dp._id + ".svg";
        render_dp_side_plot(svg, plot_reads, s1_pos, s1_str, mate_prefix, DP_SIDE1_COLOR);
        make_svg_responsive(svg);
        dp["_side_1_dp_plot_file_name"] = Settings::relative_path(svg, settings.evidence_path);
      }

      // Side 2
      g.gather(s2_seq, s2_pos, s2_str, s2_fwd, s1_tid, s1_pos, s1_fwd, D);
      vector<dp_draw_read> s2_disc;
      for (vector<dp_draw_read>::const_iterator r = g.reads().begin(); r != g.reads().end(); r++)
        if (r->category == 1) s2_disc.push_back(*r);
      {
        vector<dp_draw_read> plot_reads = g.reads();
        size_t total = dp_cap_for_display(plot_reads, max_display, static_cast<uint32_t>(s2_pos));
        if (total) dp["_side_2_dp_plot_message"] = dp_capped_message(max_display, total);
        string svg = settings.evidence_path + "/DP_SIDE_2_" + dp._id + ".svg";
        render_dp_side_plot(svg, plot_reads, s2_pos, s2_str, mate_prefix, DP_SIDE2_COLOR);
        make_svg_responsive(svg);
        dp["_side_2_dp_plot_file_name"] = Settings::relative_path(svg, settings.evidence_path);
      }

      // Joined '*' plot: match each side_1 discordant read to its side_2 mate by read-number and draw
      // the pair bridging the seam.
      {
        map<string, const dp_draw_read*> s2_by_num;
        for (size_t i = 0; i < s2_disc.size(); i++) s2_by_num[dp_read_num(s2_disc[i].name)] = &s2_disc[i];
        vector<dp_joined_pair> pairs;
        for (size_t i = 0; i < s1_disc.size(); i++) {
          map<string, const dp_draw_read*>::iterator m = s2_by_num.find(dp_read_num(s1_disc[i].name));
          if (m == s2_by_num.end()) continue;
          const dp_draw_read& a = s1_disc[i];
          const dp_draw_read& b = *m->second;
          dp_joined_pair q;
          q.s1_start = a.read_start; q.s1_end = a.read_end; q.s1_anchor = a.anchor; q.s1_rev = a.read_reversed; q.s1_name = a.name;
          q.s2_start = b.read_start; q.s2_end = b.read_end; q.s2_anchor = b.anchor; q.s2_rev = b.read_reversed; q.s2_name = b.name;
          pairs.push_back(q);
        }
        size_t total = dp_cap_for_display(pairs, max_display, static_cast<uint32_t>(s1_pos) ^ (static_cast<uint32_t>(s2_pos) << 1));
        if (total) dp["_dp_plot_message"] = dp_capped_message(max_display, total);
        string svg = settings.evidence_path + "/DP_" + dp._id + ".svg";
        render_dp_joined_plot(svg, pairs, s1_pos, s1_str, s2_pos, s2_str);
        make_svg_responsive(svg);
        dp["_dp_plot_file_name"] = Settings::relative_path(svg, settings.evidence_path);
      }
    }
  }

} // namespace breseq
