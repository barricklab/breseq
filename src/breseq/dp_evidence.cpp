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

  // Ceiling on the data-derived shared-pair floor: if even this many pairs cannot be distinguished from
  // the spurious-pair background, the run has no usable DP evidence and there is nothing to search past.
  static const int32_t kDPMaxBackgroundFloor = 100;

  // One-sided confidence level for the DP local-frequency (Clopper-Pearson) lower bound.
  static const double kDPFrequencyAlpha = 0.05;

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

  // A read pair's identifying number, shared by its two mates (the two paired read files write names
  // "<file-prefix>:<read-number>"; the prefix differs between mates, the number is common). Used to
  // pair a read at one junction side with its mate at the other side (both count/plot logic rely on it).
  static string dp_read_num(const string& name) {
    size_t colon = name.find(':');
    return (colon == string::npos) ? name : name.substr(colon + 1);
  }

  // Where one discordant alignment sits, so the OTHER mate can look up its true aligned extent.
  struct dp_mate_rec { int32_t tid, pos, end; };
  // pair number -> that pair's discordant alignments (more than two when a mate maps redundantly).
  typedef map<string, vector<dp_mate_rec> > dp_mate_index;

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
    const dp_mate_index* mates;   // discordant-alignment index for exact mate ends; NULL -> estimate
  };

  // One supporting read pair, viewed from both sides: outer (away-from-junction) and inner (junction-facing)
  // reference coordinates on each side. The inferred insert at (P1,P2) is reach1 + reach2 (dp_reach).
  struct dp_pair_ends { int32_t o1, i1, o2, i2; };


  // How far this read's MATE reaches along the reference.
  //
  // SAM stores a mate's position but not its CIGAR, so this used to be estimated as
  // mate_start + read_length_avg, which is wrong for every soft-clipped mate -- and those cluster at
  // breakpoints, exactly where this quantity decides a coordinate. The index resolves it exactly: it
  // holds every discordant alignment in the run (1.2% of the BAM, ~115k records on an LTEE clone),
  // built in one pass, so the mate's record is an O(1) lookup with no further I/O.
  //
  // Falls back to this read's own query length when the mate is not indexed (no index supplied, or a
  // mate that is not itself discordant): mates share a query length in 99.8% of pairs, which is far
  // closer than the run-wide average -- that average is short of the true length for nearly every read
  // once trimming is in play (0.5% of mate ends exact, versus 93.6% for the read's own length).
  static int32_t dp_mate_reference_span(const alignment_wrapper& a, const dp_mate_index* mates)
  {
    if (mates != NULL) {
      dp_mate_index::const_iterator it = mates->find(dp_read_num(a.read_name()));
      if (it != mates->end()) {
        int32_t mtid = static_cast<int32_t>(a.mate_reference_target_id());
        int32_t mpos = a.mate_start_1();
        for (size_t i = 0; i < it->second.size(); i++) {
          const dp_mate_rec& r = it->second[i];
          if ((r.tid == mtid) && (r.pos == mpos)) return r.end - r.pos + 1;
        }
      }
    }
    return static_cast<int32_t>(a.read_length());
  }

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
        int32_t matelen = dp_mate_reference_span(a, c.mates);
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
    //         < p means the whole mate is < p. The mate's far extent comes from dp_mate_reference_span.
    bool kept_clear = (c.s == -1) ? (rend < c.p) : (rstart > c.p);
    if (!kept_clear) return 0;
    int32_t mpos = a.mate_start_1();
    int32_t mate_len = dp_mate_reference_span(a, c.mates);
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

  // Builds the discordant-alignment index in one sequential pass over the BAM. Only paired, primary,
  // non-proper alignments are kept -- 1.2% of records on an LTEE clone (115k of 9.8M), so the index
  // costs one scan and a few MB, and every mate end afterwards is exact with no further fetches. The
  // per-side windows the scanner fetches are far too coarse for this: at a DP side only ~3-4% of the
  // records in a +/-D window are discordant, so per-call window reads would cost ~30x more I/O.
  class dp_mate_indexer : public pileup_base {
  public:
    dp_mate_indexer(const string& bam, const string& fasta)
      : pileup_base(bam, fasta), m_out(NULL) { set_print_progress(false); }

    void build(dp_mate_index& out) {
      m_out = &out;
      for (uint32_t t = 0; t < num_targets(); t++)
        do_fetch(string(target_name(t)) + ":1-" + to_string(target_length(t)));
      m_out = NULL;
    }

    void fetch_callback(const alignment_wrapper& a) {
      if (m_out == NULL) return;
      if (a.flag() & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) return;
      // BAM_FMUNMAP: a singleton is flagged paired (mark_mate_unmapped in resolve_alignments.cpp) with
      // mtid/mpos pointing at ITSELF, so without this guard it would index itself as its own mate.
      if (a.unmapped() || !a.is_paired() || a.proper_pair() || (a.flag() & BAM_FMUNMAP)) return;
      dp_mate_rec r;
      r.tid = static_cast<int32_t>(a.reference_target_id());
      r.pos = static_cast<int32_t>(a.reference_start_1());
      r.end = static_cast<int32_t>(a.reference_end_1());
      (*m_out)[dp_read_num(a.read_name())].push_back(r);
    }

  private:
    dp_mate_index* m_out;
  };

  // Counts the three read categories at a junction side (used to fill the DP evidence fields), and
  // provides a preliminary refinement pass over the discordant reads at a side.
  class dp_side_scanner : public pileup_base {
  public:
    dp_side_scanner(const string& bam, const string& fasta)
      : pileup_base(bam, fasta), m_mates(NULL), m_collect_outside(false), m_collect_pairs(false),
        m_gother_s(0) { set_print_progress(false); }

    // Exact mate ends for every classification and gathering pass from here on (see
    // dp_mate_reference_span). Not owned; must outlive the scanner.
    void set_mate_index(const dp_mate_index* mates) { m_mates = mates; }

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
    // records both mates' outer and inner (junction-facing) reference coordinates. The mate's ends come
    // from its start plus dp_mate_reference_span (the XE tag, or this read's length as a fallback).
    bool gather_pairs(const string& seq_id, int32_t p, int32_t s, bool crossing_is_forward,
                      int32_t other_tid, int32_t other_p, bool other_crossing_is_forward, double D,
                      int32_t ovl_p, int32_t ovl_other_p,
                      int32_t other_s, vector<dp_pair_ends>& out) {
      set_ctx(p, s, crossing_is_forward, other_tid, other_p, other_crossing_is_forward, D, ovl_p, ovl_other_p);
      m_pairs.clear(); m_gother_s = other_s; m_collect_pairs = true;
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
        // BAM_FMUNMAP: a singleton is flagged paired with mtid/mpos pointing at ITSELF
        // (mark_mate_unmapped in resolve_alignments.cpp). Without this guard a singleton on the
        // crossing strand whose own position happens to fall inside [other_p +/- D] would pass every
        // test below -- its "mate" is itself, on the same tid, in range, and forward (FMREVERSE is
        // never set) -- and be counted as a supporting pair. dp_classify_side_read already screens
        // these out (category 3, unpaired); this branch bypasses it, so it needs its own guard.
        if (a.unmapped() || !a.is_paired() || a.proper_pair() || (a.flag() & BAM_FMUNMAP)) return;
        if ((!a.reversed()) != m_ctx.cross_fwd) return;
        if (static_cast<int32_t>(a.mate_reference_target_id()) != m_ctx.other_tid) return;
        int32_t mpos = a.mate_start_1();
        if (mpos < m_ctx.other_p - m_ctx.D || mpos > m_ctx.other_p + m_ctx.D) return;
        bool mate_forward = (a.flag() & BAM_FMREVERSE) == 0;
        if (mate_forward != m_ctx.other_cross_fwd) return;
        int32_t rs = static_cast<int32_t>(a.reference_start_1()), re = static_cast<int32_t>(a.reference_end_1());
        int32_t matelen = dp_mate_reference_span(a, m_ctx.mates);
        // Overlapping-mate exclusion, referenced to the current best breakpoint estimate (ovl_p/
        // ovl_other_p), mirroring dp_classify_side_read: drop a pair whose two reads would overlap.
        {
          int32_t g_this  = (m_ctx.s == -1) ? (m_ctx.ovl_p - re) : (rs - m_ctx.ovl_p);
          int32_t g_other = mate_forward ? (m_ctx.ovl_other_p - (mpos + matelen - 1)) : (mpos - m_ctx.ovl_other_p);
          if (g_this + g_other < 0) return;
        }
        dp_pair_ends e;
        e.o1 = (m_ctx.s == -1) ? rs : re;                    e.i1 = (m_ctx.s == -1) ? re : rs;
        e.o2 = (m_gother_s == -1) ? mpos : (mpos + matelen - 1);
        e.i2 = (m_gother_s == -1) ? (mpos + matelen - 1) : mpos;
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
      if      (cat == 1) {
        m_supporting++;
        // Remember each supporting read's OUTSIDE (away-from-junction) coordinate alongside its pair
        // number. Intersecting the two sides then yields not just how many pairs bridge the junction but
        // how many DISTINCT (outer_1, outer_2) fragment starts they represent -- the DP analogue of JC's
        // pos_hash_score, so PCR duplicates of one molecule cannot inflate the support count.
        int32_t rstart = static_cast<int32_t>(a.reference_start_1());
        int32_t rend   = static_cast<int32_t>(a.reference_end_1());
        m_supporting_nums[dp_read_num(a.read_name())] = (m_ctx.s == -1) ? rstart : rend;
      }
      else if (cat == 2) m_concordant++;
      else if (cat == 3) m_unpaired++;
    }

    int supporting() const { return m_supporting; }
    // Read-pair numbers of the discordant (supporting) reads at the last-scanned side, each mapped to
    // that read's outside (away-from-junction) coordinate. Intersecting the two sides' key sets gives the
    // true count of pairs that bridge THIS junction (a read at a breakpoint shared with a neighboring
    // junction appears on only one side and is excluded); the paired-up values give the distinct-fragment
    // count.
    const map<string, int32_t>& supporting_nums() const { return m_supporting_nums; }
    int concordant() const { return m_concordant; }
    int unpaired()   const { return m_unpaired; }

  private:
    void set_ctx(int32_t p, int32_t s, bool crossing_is_forward,
                 int32_t other_tid, int32_t other_p, bool other_crossing_is_forward, double D,
                 int32_t ovl_p, int32_t ovl_other_p) {
      m_ctx.p = p; m_ctx.s = s; m_ctx.cross_fwd = crossing_is_forward;
      m_ctx.other_tid = other_tid; m_ctx.other_p = other_p; m_ctx.other_cross_fwd = other_crossing_is_forward;
      m_ctx.D = D; m_ctx.ovl_p = ovl_p; m_ctx.ovl_other_p = ovl_other_p;
      m_ctx.mates = m_mates;
    }
    dp_side_ctx m_ctx;
    const dp_mate_index* m_mates;
    int     m_supporting, m_concordant, m_unpaired;
    map<string, int32_t> m_supporting_nums;   // pair number -> that side's outside coordinate
    bool    m_collect_outside;
    // Supporting reads' outside (away-from-junction) coordinates, accumulated during the collect pass,
    // plus the furthest junction-facing edge (soft-clip included) seen among them.
    vector<int32_t> m_outside;
    bool    m_have_inner; int32_t m_inner_edge;
    // Pair-gathering pass (gather_pairs): both mates' ends, plus the other side's strand.
    bool    m_collect_pairs; vector<dp_pair_ends> m_pairs; int32_t m_gother_s;
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
      m_ctx.mates = NULL;   // plotting only needs the same classification, not exact mate ends
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

  // The aligned extent of one region's read of a pair (its reference start/end), as recorded per key
  // in DP_candidate_regions.csv. A key can be listed more than once in one region (the region-open
  // step re-lists every read still in the sliding window), so the two are folded to the widest span.
  struct dp_key_extent {
    int32_t read_start, read_end;
    bool    have;      // false for a pre-extent CSV: fall back to the region span
  };

  // One candidate region parsed from DP_candidate_regions.csv.
  struct dp_region_row {
    string      seq_id;
    uint32_t    start;   // lower coordinate of the region span
    uint32_t    end;     // higher coordinate
    char        strand;  // focal-read strand: 'F' or 'R'
    bool        redundant; // majority of this region's discordant reads mapped redundantly (multicopy side)
    // Distinct <read1>__<read2>__<insert_size> pair keys in this region -> that read's aligned extent.
    map<string, dp_key_extent> keys;
  };

  // Shared with MP; see the declaration in dp_evidence.h for the coordinate convention.
  void paired_region_to_side(char region_strand, uint32_t region_start, uint32_t region_end,
                             bool inner3p, int32_t& position, int32_t& strand)
  {
    bool is_forward = (region_strand == 'F');
    if (inner3p == is_forward) {
      position = static_cast<int32_t>(region_end);
      strand = -1;
    } else {
      position = static_cast<int32_t>(region_start);
      strand = +1;
    }
  }

  // Convert one DP region into a JC-style junction side (position, strand).
  static void dp_region_to_side(const dp_region_row& r, bool inner3p, int32_t& position, int32_t& strand)
  {
    paired_region_to_side(r.strand, r.start, r.end, inner3p, position, strand);
  }

  // Starting coordinate for one side of a candidate junction, taken from the aligned extents of the
  // read pairs that SEED that candidate (the keys the two regions share) instead of the region's span.
  //
  // The region span is a product of the sliding-window detector, not of the reads: a region closes
  // when a read ages out, so its far bound is a read's start plus the window width -- a coordinate no
  // read occupies, which can sit past the true breakpoint. It is also shared by every partner of a hub
  // region. The seeding reads' junction-facing ALIGNED ends cannot pass the breakpoint, so the extreme
  // one is a conservative start that is specific to this candidate:
  //   strand -1 (kept flank <= p, junction to the right) -> max(read_end)
  //   strand +1 (kept flank >= p, junction to the left)  -> min(read_start)
  // Returns false if no seeding key carried an extent (a pre-extent CSV) -- the caller then keeps the
  // region-derived position.
  static bool dp_seed_side_position(const dp_region_row& r, const vector<string>& seed_keys,
                                    int32_t strand, int32_t& position)
  {
    bool found = false;
    int32_t best = 0;
    for (size_t i = 0; i < seed_keys.size(); i++) {
      map<string, dp_key_extent>::const_iterator it = r.keys.find(seed_keys[i]);
      if (it == r.keys.end() || !it->second.have) continue;
      int32_t inner = (strand == -1) ? it->second.read_end : it->second.read_start;
      if (!found) { best = inner; found = true; }
      else        best = (strand == -1) ? max(best, inner) : min(best, inner);
    }
    if (found) position = best;
    return found;
  }

  // Shared with MP; see the declaration in dp_evidence.h.
  bool paired_library_params(const Summary& summary, bool& inner3p, double& D, double& pair_median)
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
    return false;
  }

  // Return the library orientation for the majority of paired read groups, or the (unsupported)
  // orientation name for the warning below.
  static string paired_library_orientation_name(const Summary& summary)
  {
    map<string, int> votes;
    for (PairedMappingDistanceDistributionSummaries::const_iterator it = summary.preliminary_paired_mapping_distance_distribution.begin();
         it != summary.preliminary_paired_mapping_distance_distribution.end(); it++) {
      if (!it->second.majority_orientation.empty())
        votes[it->second.majority_orientation]++;
    }
    string best_name;
    int best = 0;
    for (map<string, int>::iterator it = votes.begin(); it != votes.end(); it++) {
      if (it->second > best) { best = it->second; best_name = it->first; }
    }
    return best_name.empty() ? string("unknown") : best_name;
  }

  // DP's wrapper: the shared geometry plus DP's own warning.
  static bool dp_library_params(const Summary& summary, bool& inner3p, double& D, double& pair_median, bool warn)
  {
    if (paired_library_params(summary, inner3p, D, pair_median)) return true;
    if (warn) {
      WARN("Discordant pair (DP) evidence prediction currently supports only FR- and RF-concordant "
           "libraries. The library concordant orientation is '" +
           paired_library_orientation_name(summary) +
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
    vector<double> counts;   // LENGTH-BIAS-WEIGHTED concordant-range density: counts[d] = (# concordant
                             // pairs at distance d) * (d - trunc), for trunc <= d < cutoff, else 0.
                             // (see dp_load_insert_model for why -- DP pairs are a length-biased sample)
    double  total;           // sum of the weighted counts
    int32_t trunc;           // 2 * read_length: distances below this (overlapping mates) are excluded
    int32_t cutoff;          // distance_cutoff: distances at/above this are outliers (concordant range end)
    double  u;               // outlier insert likelihood (uniform over the inferred-insert range)
    double  p_out;           // prior probability of a chance/outlier supporting read
    dp_insert_model() : total(0.0), trunc(0), cutoff(0), u(0.0), p_out(0.0) {}
    bool ok() const { return total > 0.0 && trunc > 0 && cutoff > trunc; }
    // DP-pair inferred-insert null g(d): +/-5 bp smoothing of the length-bias-weighted concordant
    // density, a small floor, and 0 outside the concordant range [trunc, cutoff).
    double f(int32_t d) const {
      double floor = 1.0 / (total * 1000.0);
      if (d < trunc || d >= cutoff) return floor;
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

  // Total log-likelihood of a set of supporting pairs' inferred inserts given a junction at (P1,P2),
  // under the mixture L = (1-p_out)*g(d) + p_out*u with the length-bias-corrected insert null g (m.f).
  // Used by the JC-snap test: compare this at the DP position vs a candidate JC position.
  static double dp_pairs_logL(const vector<dp_pair_ends>& pr, int32_t P1, int32_t s1,
                              int32_t P2, int32_t s2, const dp_insert_model& m) {
    double lp = 0.0;
    for (vector<dp_pair_ends>::const_iterator e = pr.begin(); e != pr.end(); e++) {
      int32_t d = dp_reach(P1, e->o1, s1) + dp_reach(P2, e->o2, s2);
      lp += log((1.0 - m.p_out) * m.f(d) + m.p_out * m.u);
    }
    return lp;
  }

  // One passing JC breakpoint side pair, for snapping a nearby DP's coordinates onto it. red1/red2 are
  // the JC's side_N_redundant flags: an IS-mediated junction has exactly one redundant side, and the
  // unique-side snap below needs to know which one it is.
  struct dp_jc_sides {
    string s1seq; int32_t p1, st1; bool red1;
    string s2seq; int32_t p2, st2; bool red2;
  };

  // Load the passing (non-rejected) junction (JC) evidence breakpoints from jc_evidence.gd. A DP whose
  // pair-based edges land within the snap window of one of these validated split-read breakpoints is
  // snapped onto it (see the snap block in predict_discordant_pairs).
  static void dp_load_passing_jcs(const string& fn, vector<dp_jc_sides>& out) {
    out.clear();
    if (!file_exists(fn.c_str())) return;
    cGenomeDiff jc_gd(fn);
    diff_entry_list_t jl = jc_gd.get_list(make_vector<gd_entry_type>(JC));
    for (diff_entry_list_t::iterator it = jl.begin(); it != jl.end(); it++) {
      cDiffEntry& j = **it;
      if (j.entry_exists("reject")) continue;
      dp_jc_sides s;
      s.s1seq = j[SIDE_1_SEQ_ID]; s.p1 = from_string<int32_t>(j[SIDE_1_POSITION]); s.st1 = from_string<int32_t>(j[SIDE_1_STRAND]);
      s.s2seq = j[SIDE_2_SEQ_ID]; s.p2 = from_string<int32_t>(j[SIDE_2_POSITION]); s.st2 = from_string<int32_t>(j[SIDE_2_STRAND]);
      s.red1 = j.entry_exists(SIDE_1_REDUNDANT) && (j[SIDE_1_REDUNDANT] == "1");
      s.red2 = j.entry_exists(SIDE_2_REDUNDANT) && (j[SIDE_2_REDUNDANT] == "1");
      out.push_back(s);
    }
  }

  // Is this DP side sitting on (or within 50 bp of) an annotated repeat boundary? Used to tell an
  // IS-mediated DP's IS side from its unique side while placing it -- before the repeat-end snap below
  // has set side_N_repeat, and for a side whose reads were not themselves flagged redundant.
  static bool dp_side_on_repeat(cReferenceSequences& ref_seq_info, const string& seq_id, int32_t pos, int32_t strand)
  {
    int32_t md = 50;
    return cReferenceSequences::find_closest_repeat_region_boundary(pos, ref_seq_info[seq_id].m_repeats, md, strand, true) != NULL;
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
    m.cutoff = static_cast<int32_t>(dcut + 0.5);
    // Length-bias-correct the null. DP supporting pairs are NOT an unbiased sample of the concordant
    // fragment distribution f(L): a spanning fragment only survives as a DP pair if the breakpoint
    // falls in its unsequenced middle gap (else a read crosses it and is soft-clipped into a JC
    // -M1/-M2 read). That selection weights each fragment by L*f(L) [inspection-paradox: fragments
    // crossing a fixed point are length-biased by L] times (L-2r)/L [fraction of crossing positions
    // leaving both reads clean] = f(L)*(L-2r). So reweight the concordant-range histogram by (d-trunc)
    // and renormalize; f() then reads out this corrected null g directly. (Uncorrected f is ~60 bp
    // too short for DP pairs -- validated: median of g = 715 vs measured DP inferred insert 713.)
    if (m.cutoff > m.trunc) {
      m.total = 0.0;
      for (int32_t d = 0; d < static_cast<int32_t>(m.counts.size()); d++) {
        m.counts[d] = (d >= m.trunc && d < m.cutoff) ? m.counts[d] * (d - m.trunc) : 0.0;
        m.total += m.counts[d];
      }
    }
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
    // P==1 gives -log10 = -0.0, which compares EQUAL to 0.0 (so a "< 0.0" test misses it) but still
    // prints with its sign as "-0.0". Use <= so the negative zero is normalized away.
    if (sc <= 0.0) sc = 0.0;
    return (sc > kDPMaxScore) ? kDPMaxScore : sc;
  }

  // ---------------------------------------------------------------------------------------------
  // Local frequency test -- the accept/reject gate that survives a short pair distance distribution.
  //
  // A fragment whose unsequenced middle gap falls across the breakpoint is observed as a DISCORDANT
  // pair when that molecule carries the variant junction and as a CONCORDANT pair when it does not.
  // The two are the variant and reference observations of one sampling event, so k / (k + c) is the
  // local frequency of the structural variant -- self-normalized to local coverage, local repeat
  // content and local deletions, all of which cancel out of the ratio.
  //
  // Why this and not the concordant-pair skew below: the skew asks whether k is low compared with how
  // many concordant pairs span a NORMAL position, an expectation of (coverage / 2*read_length) *
  // E[(insert - 2*read_length)+]. As the insert distribution shortens toward twice the read length
  // that expectation collapses toward zero, P(crossing <= k) goes to 1 for every k, and the skew test
  // accepts everything -- which is exactly the regime where DP evidence over-predicts. The frequency
  // test is indifferent: k and c shrink together and their ratio does not move.
  //
  // Reported as an exact (Clopper-Pearson) lower confidence bound -- binomial_frequency_lower_bound
  // in stats.h, shared with the RA and JC frequency tests -- so a candidate is never rejected merely
  // for having small counts, only for being confidently LOW frequency. k of k observations therefore
  // always passes, however small k is; ruling those out is the background test's job.

  // Mean of a crossing distribution over ALL bins, INCLUDING crossing = 0 -- the expected number of
  // concordant pairs spanning a normal position. Deliberately uncensored: dp_crossing_censor drops the
  // zero bin as a "deletion spike", but on a library whose fragments are barely longer than two reads
  // zero IS the normal state (70% of positions is routine), and hiding it makes the skew null describe
  // a library that was never sequenced. This uncensored mean is what says whether the skew test has any
  // power at all at this coverage.
  static double dp_crossing_mean(const vector<double>& hist)
  {
    double n = 0.0, s = 0.0;
    for (int32_t c = 0; c < static_cast<int32_t>(hist.size()); c++) { n += hist[c]; s += hist[c] * static_cast<double>(c); }
    return (n > 0.0) ? (s / n) : 0.0;
  }

  // ---------------------------------------------------------------------------------------------
  // Background model for the number of read pairs that two UNRELATED candidate regions happen to
  // share -- the null the old hard-coded "at least 5 shared pairs" was standing in for.
  //
  // Fitted from the observed edge-weight spectrum itself, using only the ratios between the first
  // few bins. For a Poisson, n(w+1)/n(w) = mu/(w+1); for a negative binomial with size r and
  // p = mu/(mu+r), n(w+1)/n(w) = p*(w+r)/(w+1). Ratios are unaffected by the fact that weight-0
  // region pairs are never observed, so no zero-truncation correction is needed to estimate the
  // parameters. A negative binomial rather than a Poisson because spurious discordant pairs are not
  // placed uniformly -- they pile up at repeats, rRNA operons and mismapping hotspots -- and a
  // Poisson null is anti-conservative in exactly those loci.
  struct dp_background_model {
    bool   ok;
    double mu;        // background mean shared pairs between two regions
    double size;      // negative binomial size; <= 0 means Poisson
    double n_edges;   // candidate junctions tested (the multiple-testing correction)
    dp_background_model() : ok(false), mu(0.0), size(0.0), n_edges(0.0) {}

    // P(X >= k | X >= 1): the conditional tail, since only region pairs sharing at least one read
    // pair ever become an edge.
    double tail(int32_t k) const {
      if (!ok || k <= 1) return 1.0;
      double p0 = (size > 0.0) ? dp_nbinom_pmf(0, size, mu) : exp(-mu);
      double denom = 1.0 - p0;
      if (!(denom > 0.0)) return 1.0;
      double below = 0.0;                                    // P(X <= k-1)
      for (int32_t i = 0; i < k; i++)
        below += (size > 0.0) ? dp_nbinom_pmf(i, size, mu) : exp(-mu + i * log(mu) - lgamma(i + 1.0));
      double above = 1.0 - below;
      if (above < 0.0) above = 0.0;
      return above / denom;
    }
    double e_value(int32_t k) const { return n_edges * tail(k); }
  };

  // Fit the background from the edge-weight histogram (weights >= 1). `mu_floor` is the closed-form
  // uniform-placement expectation (derived from p_out); the fit is never allowed below it, so a
  // degenerate spectrum can only make the test more conservative, never less.
  static void dp_fit_background(const vector<int>& weights, double mu_floor, dp_background_model& m)
  {
    m.ok = false;
    m.n_edges = static_cast<double>(weights.size());
    double n1 = 0, n2 = 0, n3 = 0;
    for (size_t i = 0; i < weights.size(); i++) {
      if      (weights[i] == 1) n1++;
      else if (weights[i] == 2) n2++;
      else if (weights[i] == 3) n3++;
    }
    double mu = 0.0, size = 0.0;
    if (n1 > 0.0 && n2 > 0.0) {
      double r1 = n2 / n1;
      // Try the negative binomial first; it needs a third bin and must come out over-dispersed.
      if (n3 > 0.0) {
        double r2 = n3 / n2;
        double q = 1.5 * r2 / r1;              // = (2 + size) / (1 + size)
        if (q > 1.0 && q < 2.0) {
          size = (2.0 - q) / (q - 1.0);
          double p = 2.0 * r1 / (1.0 + size);
          if (p > 0.0 && p < 1.0) mu = size * p / (1.0 - p);
          else size = 0.0;
        }
      }
      if (!(size > 0.0 && mu > 0.0)) { size = 0.0; mu = 2.0 * r1; }   // Poisson fallback
    }
    if (mu < mu_floor) { mu = mu_floor; size = 0.0; }
    if (!(mu > 0.0) || m.n_edges <= 0.0) return;
    m.mu = mu; m.size = size; m.ok = true;
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
    string s = use_empirical ? ("empirical (" + npos + ")")
                             : ("negative binomial (fallback; only " + npos + ")");
    // State whether the skew test is actually in force. Its effect size is the expected number of
    // concordant pairs spanning a normal position; when that is small (a paired-mapping distance
    // distribution barely longer than two reads) the test cannot discriminate and only the local
    // frequency test rejects. Saying so is the quickest way to understand a short-insert run.
    double mean_x = dp_crossing_mean(hist_ref);
    s += ". Expected concordant pairs spanning a normal position in " + R + ": " + to_string(mean_x, 1, false);
    if (settings.discordant_pair_minimum_crossing > 0.0 && mean_x < settings.discordant_pair_minimum_crossing)
      s += " -- below --discordant-pair-minimum-crossing, so the skew test is reported but does not reject; "
           "DP evidence is accepted or rejected on its local frequency alone";
    return s;
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
    // Uncensored mean of that reference distribution -- the skew test's effect size (see the power gate
    // in the emit loop). Kept separate from the censored window the skew null itself is built on.
    double crossing_mean_ref = have_crossing ? dp_crossing_mean(crossing_hist_ref) : 0.0;

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

    // Passing split-read junction (JC) breakpoints, for snapping DP coordinates onto a validated
    // junction when the pair evidence is consistent with it (see the snap block below).
    //
    // The snap window is how far a pair-derived coordinate can plausibly sit from the true breakpoint.
    // Placement puts it at the innermost supporting read's aligned edge, and the breakpoint lies
    // somewhere in that fragment's unsequenced middle gap, so the error is bounded by the LARGEST such
    // gap: distance_cutoff - 2*read_length. (It was previously median - 2*read_length, which is not a
    // bound on anything and goes NEGATIVE whenever the paired-mapping distance distribution is shorter
    // than two reads -- exactly the short-insert libraries at issue here. Every `abs(...) > snap_win`
    // test then always fired, so neither the JC snap nor the circular-origin snap could ever run, and
    // circular-origin artifacts were reported as ordinary DP junctions.) Widening this is safe: each
    // snap is separately gated on the supporting pairs' inferred inserts favoring the candidate.
    vector<dp_jc_sides> passing_jcs;
    if (have_insert) dp_load_passing_jcs(settings.jc_genome_diff_file_name, passing_jcs);
    double read_len_avg = summary.sequence_conversion.read_length_avg;
    int32_t snap_win = static_cast<int32_t>(max(read_len_avg, distance_cutoff - 2.0 * read_len_avg) + 0.5);

    //
    // Step 1: Re-read the candidate regions CSV.
    //  columns: seq_id,start,end,strand,orientation,length,max_discordant_count,redundant,discordant_pairs
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

        // 'redundant' is column index 7 (1 = tie-broken multicopy side). Absent in older CSVs.
        r.redundant = (f.size() >= 8 && f[7] == "1");

        // The keys field is column index 8 (empty if the region had no descriptors). Each entry is
        // <read1>__<read2>__<insert_size>__<read_start>__<read_end> (identify_mutations.cpp's
        // dp_descriptor): the identity key is the first three fields, the last two are this region's
        // read's aligned extent. A CSV written before the extents existed has only the three.
        string key_field = (f.size() >= 9) ? f[8] : "";
        vector<string> keys = split(key_field, ";");
        for (size_t i = 0; i < keys.size(); i++) {
          if (keys[i].empty()) continue;
          vector<string> kf = split(keys[i], "__");
          dp_key_extent x; x.read_start = 0; x.read_end = 0; x.have = false;
          string key = keys[i];
          if (kf.size() >= 5) {
            key = kf[0] + "__" + kf[1] + "__" + kf[2];
            x.read_start = from_string<int32_t>(kf[3]);
            x.read_end   = from_string<int32_t>(kf[4]);
            x.have = true;
          }
          // Fold a repeated key (re-listed by a later region-open snapshot) to the widest extent.
          map<string, dp_key_extent>::iterator prev = r.keys.find(key);
          if (prev == r.keys.end()) r.keys[key] = x;
          else if (x.have) {
            if (!prev->second.have) prev->second = x;
            else {
              prev->second.read_start = min(prev->second.read_start, x.read_start);
              prev->second.read_end   = max(prev->second.read_end,   x.read_end);
            }
          }
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
      for (map<string, dp_key_extent>::const_iterator k = regions[ri].keys.begin(); k != regions[ri].keys.end(); k++) {
        key_to_regions[k->first].insert(static_cast<int>(ri));
      }
    }
    // The keys behind each edge's weight are kept alongside it: they are the pairs that SEED that
    // candidate junction, and each side's starting coordinate is taken from their aligned extents
    // (dp_seed_side_position) rather than from the region span.
    map<pair<int, int>, int> edge_weight;
    map<pair<int, int>, vector<string> > edge_keys;
    for (map<string, set<int> >::const_iterator it = key_to_regions.begin(); it != key_to_regions.end(); it++) {
      if (it->second.size() == 2) {
        set<int>::const_iterator si = it->second.begin();
        int a = *si; ++si;
        int b = *si;
        pair<int, int> ab = make_pair(min(a, b), max(a, b));
        edge_weight[ab]++;
        edge_keys[ab].push_back(it->first);
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
    // Step 3b: Fit the spurious-pair background to the edge-weight spectrum, and derive the minimum
    // shared-pair count worth examining. This replaces a fixed "at least 5 shared pairs" threshold with
    // one the run's own data sets: a candidate whose weight is explained by the background across all
    // edges.size() junctions tested carries no information, so nothing below that floor is emitted.
    // Because the edges are sorted by descending weight, this is applied as the loop's break point.
    //
    dp_background_model background;
    {
      vector<int> weights;
      weights.reserve(edges.size());
      for (size_t i = 0; i < edges.size(); i++) weights.push_back(edges[i].first);
      // Closed-form floor: p_out is the probability that at least one chance discordant pair lands in
      // this DP's two windows under uniform placement, so -log(1 - p_out) is the corresponding mean.
      double mu_floor = (have_insert && insert_model.p_out > 0.0 && insert_model.p_out < 1.0)
                        ? -log(1.0 - insert_model.p_out) : 0.0;
      dp_fit_background(weights, mu_floor, background);
    }
    int32_t min_pairs = max(1, settings.discordant_pair_minimum_pairs);
    if (background.ok && settings.discordant_pair_background_e_value_cutoff > 0.0) {
      // Smallest weight whose background e-value clears the cutoff (tail is monotone in k).
      for (int32_t k = min_pairs; k <= kDPMaxBackgroundFloor; k++) {
        if (background.e_value(k) <= settings.discordant_pair_background_e_value_cutoff) { min_pairs = k; break; }
        min_pairs = k;
      }
    }
    // Record the model + gates for summary.html / summary.json (see DiscordantPairSummary).
    {
      DiscordantPairSummary& dps = summary.discordant_pair;
      dps.read_length_avg = summary.sequence_conversion.read_length_avg;
      const PairedMappingDistanceDistributionSummaries& pmdd = summary.preliminary_paired_mapping_distance_distribution;
      double best_mapped = -1.0;
      for (PairedMappingDistanceDistributionSummaries::const_iterator it = pmdd.begin(); it != pmdd.end(); it++) {
        if (it->second.mapped_pairs > best_mapped) {
          best_mapped = it->second.mapped_pairs;
          dps.pair_distance_median = it->second.median;
          dps.pair_distance_mad = it->second.mad;
          dps.pair_distance_cutoff = it->second.distance_cutoff;
          dps.pair_orientation = it->second.majority_orientation;
        }
      }
      dps.crossing_reference_seq_id = crossing_ref_seq_id;
      dps.crossing_reference_coverage = crossing_C_ref;
      dps.expected_concordant_crossing = crossing_mean_ref;
      dps.crossing_use_empirical = crossing_use_empirical;
      dps.crossing_normal_positions = static_cast<uint64_t>(crossing_N);
      dps.frequency_cutoff = settings.discordant_pair_frequency_cutoff;
      dps.skew_cutoff = settings.discordant_pair_skew_cutoff;
      dps.minimum_crossing = settings.discordant_pair_minimum_crossing;
      // Reported at the crossing reference's own coverage; the per-item gate rescales by each side's.
      dps.skew_in_force = (settings.discordant_pair_minimum_crossing <= 0.0)
                          || (crossing_mean_ref >= settings.discordant_pair_minimum_crossing);
      dps.background_e_value_cutoff = settings.discordant_pair_background_e_value_cutoff;
      dps.minimum_pairs_option = settings.discordant_pair_minimum_pairs;
      dps.minimum_pairs_used = min_pairs;
      dps.candidate_junctions = static_cast<uint64_t>(edges.size());
      dps.background_mean = background.ok ? background.mu : 0.0;
      dps.background_size = background.ok ? background.size : 0.0;
    }

    cerr << "  Discordant pair (DP) background: " << edges.size() << " candidate junctions, "
         << (background.ok ? ((background.size > 0.0)
               ? ("negative binomial mean " + to_string(background.mu, 4, true) + ", size " + to_string(background.size, 3, false))
               : ("Poisson mean " + to_string(background.mu, 4, true)))
           : string("not fit"))
         << "; requiring at least " << min_pairs << " shared read pairs." << endl;

    //
    // Step 4: Emit one DP item per edge with >= min_pairs shared read pairs.
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

    // One pass over the BAM to index every discordant alignment, so each mate's aligned extent is
    // exact rather than estimated from a read length (see dp_mate_reference_span). Must outlive the
    // scanner, which only borrows it.
    dp_mate_index mate_index;
    if (scanner) {
      dp_mate_indexer indexer(settings.reference_bam_file_name, settings.reference_fasta_file_name);
      indexer.build(mate_index);
      scanner->set_mate_index(&mate_index);
      cerr << "  Discordant pair (DP): indexed " << mate_index.size()
           << " read pairs with a discordant alignment (exact mate ends)." << endl;
    }

    // Items already emitted, keyed by the six side fields they were placed at. A DP item is identified
    // in a .gd by exactly those fields, so a second item at the same breakpoint is not a near-duplicate
    // to be tolerated -- it is a fatal duplicate on write. See the fold-in block inside the loop.
    map<string, diff_entry_ptr_t> emitted_by_breakpoint;

    for (size_t e = 0; e < edges.size(); e++) {
      int weight = edges[e].first;
      if (weight < min_pairs) break; // edges are sorted descending; nothing left qualifies
      int a = edges[e].second.first;
      int b = edges[e].second.second;

      // Each side's strand comes from its region; its starting coordinate comes from the aligned
      // extents of the pairs that seed THIS edge, falling back to the region span when the CSV
      // carries no extents.
      int32_t pos_a, strand_a, pos_b, strand_b;
      dp_region_to_side(regions[a], inner3p, pos_a, strand_a);
      dp_region_to_side(regions[b], inner3p, pos_b, strand_b);
      {
        const vector<string>& seed_keys = edge_keys[edges[e].second];
        dp_seed_side_position(regions[a], seed_keys, strand_a, pos_a);
        dp_seed_side_position(regions[b], seed_keys, strand_b, pos_b);
      }

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
        vector<dp_pair_ends> pr;
        if (scanner->gather_pairs(s1_seq_id, s1_pos, s1_strand, s1_fwd, s2_tid, s2_pos, s2_fwd, distance_cutoff, init1, init2, s2_strand, pr)
            && pr.size() >= 2) {
          int32_t r1 = dp_robust_edge(pr, /*this_is_side1=*/true,  s1_strand, s2_pos, s2_strand, insert_model);
          int32_t r2 = dp_robust_edge(pr, /*this_is_side1=*/false, s2_strand, r1,     s1_strand, insert_model);
          s1_pos = max(1, min(scanner->seq_length(s1_tid), r1));
          s2_pos = max(1, min(scanner->seq_length(s2_tid), r2));
        }
      }

      // Snap the coarse pair-based edges onto a better-supported breakpoint. Three candidate snaps are
      // tried in turn; each is gated by the SAME Bayesian requirement: accept a candidate only if the
      // supporting pairs' inferred inserts do NOT favor the current position over it, per pair, by more
      // than the floor set below -- i.e. lp(candidate) - lp(current) >= floor under the length-bias-
      // corrected insert model. So a real near-origin junction, or a DP whose pairs contradict a nearby
      // transposable element, is neither moved nor marked.
      bool circular_dp = false, side1_repeat = false, side2_repeat = false;
      bool side1_jc_snapped = false, side2_jc_snapped = false;
      // Which side's supporting reads mapped redundantly (from the DP candidate regions). A redundant
      // IS side that is not pinned to a JC gets a guaranteed IS-end snap in the final pass below.
      bool side1_redundant_reads = a_is_side_1 ? regions[a].redundant : regions[b].redundant;
      bool side2_redundant_reads = a_is_side_1 ? regions[b].redundant : regions[a].redundant;
      if (have_insert && scanner) {
        vector<dp_pair_ends> pr;
        if (scanner->gather_pairs(s1_seq_id, s1_pos, s1_strand, s1_fwd, s2_tid, s2_pos, s2_fwd, distance_cutoff, init1, init2, s2_strand, pr)
            && !pr.empty()) {
          // Accept a snap unless the supporting pairs favor the current position over the candidate.
          // The evidence is a SUM over those pairs, so a fixed total threshold means something
          // different at 10 pairs than at 500: at ~150 pairs the old log(1/3) worked out to -0.007
          // nats per pair -- "the pairs must not disagree at all", about +/-18 bp. Scaling the floor
          // with the pair count makes the criterion depth-independent and gives it a physical reading.
          // Moving a breakpoint by d changes the mean per-pair log-likelihood by roughly d^2/(2*sigma^2)
          // with sigma the insert spread, so a floor of -0.25 nats/pair accepts a move of up to about
          // 0.7*sigma and rejects anything beyond -- on an LTEE clone (sigma ~150) that is a ~100 bp
          // move accepted at -0.22/pair and a 538 bp move rejected at -4.2/pair, where the old total
          // threshold could not tell the two apart (both were simply "far outside").
          const double kSnapPerPairLBF = -0.25;
          const double snap_lbf_floor = kSnapPerPairLBF * static_cast<double>(pr.size());
          double lp_cur = dp_pairs_logL(pr, s1_pos, s1_strand, s2_pos, s2_strand, insert_model);

          // (1) JC snap: a passing split-read junction (JC) gives a base-resolution breakpoint. Snap to
          // the nearest passing JC matching both sides within snap_win (either side order).
          if (!passing_jcs.empty()) {
            const dp_jc_sides* best = NULL; int32_t best_off = 0; int32_t jp1 = 0, jp2 = 0;
            for (vector<dp_jc_sides>::const_iterator j = passing_jcs.begin(); j != passing_jcs.end(); j++) {
              for (int order = 0; order < 2; order++) {
                const string& a_seq = order ? j->s2seq : j->s1seq; int32_t a_pos = order ? j->p2 : j->p1;
                const string& b_seq = order ? j->s1seq : j->s2seq; int32_t b_pos = order ? j->p1 : j->p2;
                if (a_seq != s1_seq_id || b_seq != s2_seq_id) continue;
                if (abs(s1_pos - a_pos) > snap_win || abs(s2_pos - b_pos) > snap_win) continue;
                int32_t off = abs(s1_pos - a_pos) + abs(s2_pos - b_pos);
                if (!best || off < best_off) { best = &(*j); best_off = off; jp1 = a_pos; jp2 = b_pos; }
              }
            }
            if (best && (jp1 != s1_pos || jp2 != s2_pos)
                && dp_pairs_logL(pr, jp1, s1_strand, jp2, s2_strand, insert_model) - lp_cur >= snap_lbf_floor) {
              s1_pos = max(1, min(scanner->seq_length(s1_tid), jp1));
              s2_pos = max(1, min(scanner->seq_length(s2_tid), jp2));
              lp_cur = dp_pairs_logL(pr, s1_pos, s1_strand, s2_pos, s2_strand, insert_model);
              side1_jc_snapped = true; side2_jc_snapped = true;   // both sides pinned to this JC
            }
          }

          // (1b) Unique-side JC snap. An IS-mediated DP and the JC describing the same insertion do NOT
          // agree on their IS side: the DP's is chosen by the per-locus copy vote in resolve_alignments,
          // which has no reason to pick the copy the split reads landed on, so (1) -- which requires both
          // sides -- never fires for them. The DP is then left with a unique-side coordinate that can sit
          // tens of bp off the junction it agrees with, because the pair gathering approximates a mate's
          // junction-facing end as mate_start + read_length, which overshoots whenever that mate is
          // soft-clipped at the breakpoint. Match on the unique side alone and pin ONLY that coordinate;
          // the IS side stays on its voted copy (the reads cannot say which copy it is), and the final
          // redundant-IS pass below puts it exactly on that copy's element end.
          if (!side1_jc_snapped && !side2_jc_snapped) {
            bool s1_is = side1_redundant_reads || dp_side_on_repeat(ref_seq_info, s1_seq_id, s1_pos, s1_strand);
            bool s2_is = side2_redundant_reads || dp_side_on_repeat(ref_seq_info, s2_seq_id, s2_pos, s2_strand);
            if (s1_is != s2_is) {                      // exactly one IS side, one unique side
              bool unique_is_side1 = !s1_is;
              const string& u_seq = unique_is_side1 ? s1_seq_id : s2_seq_id;
              int32_t u_pos = unique_is_side1 ? s1_pos : s2_pos;
              int32_t u_str = unique_is_side1 ? s1_strand : s2_strand;

              int32_t best_pos = 0, best_off = 0; bool found = false;
              for (vector<dp_jc_sides>::const_iterator j = passing_jcs.begin(); j != passing_jcs.end(); j++) {
                // The JC must be IS-mediated too. Classify its IS side the same way as the DP's, by the
                // repeat annotation -- NOT by side_N_redundant alone: a side whose split reads uniquely
                // matched one variant copy is a repeat side that is nonetheless flagged redundant=0
                // (the same trap combine_DP_with_MOB_by_unique_side documents).
                bool j1_is = j->red1 || dp_side_on_repeat(ref_seq_info, j->s1seq, j->p1, j->st1);
                bool j2_is = j->red2 || dp_side_on_repeat(ref_seq_info, j->s2seq, j->p2, j->st2);
                if (j1_is == j2_is) continue;
                const string& ju_seq = j1_is ? j->s2seq : j->s1seq;
                int32_t ju_pos = j1_is ? j->p2 : j->p1;
                int32_t ju_str = j1_is ? j->st2 : j->st1;
                if (ju_seq != u_seq || ju_str != u_str) continue;
                int32_t off = abs(u_pos - ju_pos);
                if (off > snap_win) continue;
                if (!found || off < best_off) { found = true; best_off = off; best_pos = ju_pos; }
              }

              // Gated exactly like the snaps around it: accept unless the supporting pairs' inferred
              // inserts favor the current position by more than 3x.
              if (found && (best_pos != u_pos)) {
                int32_t c1 = unique_is_side1 ? best_pos : s1_pos;
                int32_t c2 = unique_is_side1 ? s2_pos : best_pos;
                if (dp_pairs_logL(pr, c1, s1_strand, c2, s2_strand, insert_model) - lp_cur >= snap_lbf_floor) {
                  if (unique_is_side1) {
                    s1_pos = max(1, min(scanner->seq_length(s1_tid), best_pos)); side1_jc_snapped = true;
                  } else {
                    s2_pos = max(1, min(scanner->seq_length(s2_tid), best_pos)); side2_jc_snapped = true;
                  }
                  lp_cur = dp_pairs_logL(pr, s1_pos, s1_strand, s2_pos, s2_strand, insert_model);
                }
              }
            }
          }

          // (2) Circular-origin snap: a DP that just reconnects the two ends of a circular seq_id is an
          // artifact (like a CIRCULAR_CHROMOSOME JC). Snap onto (1, length) and mark it ignored.
          if (s1_seq_id == s2_seq_id && s1_strand != s2_strand && ref_seq_info.is_circular(s1_seq_id)) {
            int32_t L = static_cast<int32_t>(ref_seq_info.get_sequence_length(s1_seq_id));
            int32_t c1 = s1_pos, c2 = s2_pos; bool cand = false;
            if (s1_pos <= snap_win && s2_pos >= L - snap_win + 1)      { c1 = 1; c2 = L; cand = true; }
            else if (s2_pos <= snap_win && s1_pos >= L - snap_win + 1) { c1 = L; c2 = 1; cand = true; }
            // Mark it circular whenever it reconnects the origin (cand) and the pairs are consistent with
            // (1, L) -- NOT only when the coordinates move. A DP already pinned to the origin (e.g. it
            // snapped onto the circular JC in step 1) is at c1==s1_pos/c2==s2_pos, for which the log-
            // likelihood difference is 0 (>= snap_lbf_floor), so it is still flagged/ignored.
            if (cand && dp_pairs_logL(pr, c1, s1_strand, c2, s2_strand, insert_model) - lp_cur >= snap_lbf_floor) {
              s1_pos = c1; s2_pos = c2; circular_dp = true;
              lp_cur = dp_pairs_logL(pr, s1_pos, s1_strand, s2_pos, s2_strand, insert_model);
            }
          }

          // (3) Transposable-element (repeat) end snap: a side within 50 bp of a repeat/IS boundary is
          // snapped onto that boundary and marked redundant (same flags as JC), so it renders orange.
          if (!circular_dp) {
            int32_t md1 = 50;
            cFeatureLocation* is1 = cReferenceSequences::find_closest_repeat_region_boundary(s1_pos, ref_seq_info[s1_seq_id].m_repeats, md1, s1_strand);
            if (is1) {
              int32_t c1 = (s1_strand == -1) ? is1->get_end_1() : is1->get_start_1();
              if (dp_pairs_logL(pr, c1, s1_strand, s2_pos, s2_strand, insert_model) - lp_cur >= snap_lbf_floor) {
                s1_pos = max(1, min(scanner->seq_length(s1_tid), c1)); side1_repeat = true;
                lp_cur = dp_pairs_logL(pr, s1_pos, s1_strand, s2_pos, s2_strand, insert_model);
              }
            }
            int32_t md2 = 50;
            cFeatureLocation* is2 = cReferenceSequences::find_closest_repeat_region_boundary(s2_pos, ref_seq_info[s2_seq_id].m_repeats, md2, s2_strand);
            if (is2) {
              int32_t c2 = (s2_strand == -1) ? is2->get_end_1() : is2->get_start_1();
              if (dp_pairs_logL(pr, s1_pos, s1_strand, c2, s2_strand, insert_model) - lp_cur >= snap_lbf_floor) {
                s2_pos = max(1, min(scanner->seq_length(s2_tid), c2)); side2_repeat = true;
              }
            }
          }
        }
      }

      // Final pass (always runs, no Bayesian gate): guarantee that a redundant IS side NOT already
      // pinned to a JC lands exactly on its IS element end -- on its OWN specific copy (no cross-copy
      // move). This is the catch-all for cases the JC snap did not cover, so DP always reports a clean
      // IS-end coordinate from which predict_mutations can read the specific copy.
      if (!circular_dp) {
        if (side1_redundant_reads && !side1_jc_snapped) {
          int32_t md1 = 50;
          cFeatureLocation* is1 = cReferenceSequences::find_closest_repeat_region_boundary(s1_pos, ref_seq_info[s1_seq_id].m_repeats, md1, s1_strand, true);
          if (is1) {
            int32_t c1 = (s1_strand == -1) ? is1->get_end_1() : is1->get_start_1();
            s1_pos = max(1, min(static_cast<int32_t>(ref_seq_info.get_sequence_length(s1_seq_id)), c1));
            side1_repeat = true;
          }
        }
        if (side2_redundant_reads && !side2_jc_snapped) {
          int32_t md2 = 50;
          cFeatureLocation* is2 = cReferenceSequences::find_closest_repeat_region_boundary(s2_pos, ref_seq_info[s2_seq_id].m_repeats, md2, s2_strand, true);
          if (is2) {
            int32_t c2 = (s2_strand == -1) ? is2->get_end_1() : is2->get_start_1();
            s2_pos = max(1, min(static_cast<int32_t>(ref_seq_info.get_sequence_length(s2_seq_id)), c2));
            side2_repeat = true;
          }
        }
      }

      // More than one graph edge can end up at ONE breakpoint. A single breakpoint shoulder is split
      // into several candidate regions whenever its (strand x orientation) bin's in-window count dips
      // below --discordant-pair-seed, and each fragment forms its own edge with the same partner; the
      // placement and snapping passes above then pull those edges onto the same coordinates. Every
      // field below this point is recomputed from the BAM at the placed positions, so the resulting
      // items differ only in candidate_discordant_count -- and writing both aborts the run.
      //
      // Fold this edge into the item already emitted here rather than emitting a second one. Done
      // BEFORE the BAM rescan below, so a folded edge costs no scan and does not enter the run tallies
      // as a separately tested item.
      const string breakpoint_key = s1_seq_id + ":" + to_string(s1_pos) + ":" + to_string(s1_strand)
                                  + "|" + s2_seq_id + ":" + to_string(s2_pos) + ":" + to_string(s2_strand);
      {
        map<string, diff_entry_ptr_t>::iterator prev = emitted_by_breakpoint.find(breakpoint_key);
        if (prev != emitted_by_breakpoint.end()) {
          cDiffEntry& kept = *(prev->second);
          // Summing is exact rather than approximate: a pair key is counted toward an edge only when it
          // lands in exactly two regions (see the size() == 2 guard in Step 2), so the edges partition
          // the keys between them and no key is counted twice here.
          const int32_t merged_weight = from_string<int32_t>(kept["candidate_discordant_count"]) + weight;
          kept["candidate_discordant_count"] = to_string(merged_weight);
          // Keep the reported expectation describing the count now shown beside it.
          if (background.ok)
            kept["background_e_value"] = to_string(background.e_value(merged_weight), 3, true);
          // Redundancy is the only other thing that can differ: it comes from the regions, and this
          // edge's regions are not the kept item's. A side seen as redundant by either stays redundant.
          if (side1_repeat || side1_redundant_reads) {
            kept["side_1_annotate_key"] = "repeat";
            kept[SIDE_1_REDUNDANT] = "1";
          }
          if (side2_repeat || side2_redundant_reads) {
            kept["side_2_annotate_key"] = "repeat";
            kept[SIDE_2_REDUNDANT] = "1";
          }
          summary.discordant_pair.items_merged_duplicate++;
          continue;
        }
      }

      cDiffEntry dp(DP);
      dp[SIDE_1_SEQ_ID]  = s1_seq_id;
      dp[SIDE_1_POSITION] = to_string(s1_pos);
      dp[SIDE_1_STRAND]  = to_string(s1_strand);
      dp[SIDE_2_SEQ_ID]  = s2_seq_id;
      dp[SIDE_2_POSITION] = to_string(s2_pos);
      dp[SIDE_2_STRAND]  = to_string(s2_strand);
      // Circular-origin artifact -> ignore (hidden in output like a CIRCULAR_CHROMOSOME JC).
      if (circular_dp) dp[IGNORE] = "CIRCULAR_CHROMOSOME";
      // Per-side annotate_key drives the HTML highlight exactly as for JC ("repeat" -> orange). A
      // side is redundant if EITHER it snapped onto a repeat/IS boundary (side_N_repeat) OR its
      // supporting reads mapped redundantly (side_N_redundant_reads, derived above from
      // region.redundant). Both feed the same side_N_redundant flag JC uses.
      bool side1_red = side1_repeat || side1_redundant_reads;
      bool side2_red = side2_repeat || side2_redundant_reads;
      dp["side_1_annotate_key"] = side1_red ? "repeat" : "gene";
      dp["side_2_annotate_key"] = side2_red ? "repeat" : "gene";
      if (side1_red) dp[SIDE_1_REDUNDANT] = "1";
      if (side2_red) dp[SIDE_2_REDUNDANT] = "1";

      // Heuristic count from region overlap; kept so we can see the rescan hold steady or increase.
      dp["candidate_discordant_count"] = to_string(weight);

      int k_support = weight;
      int k_distinct = weight;
      bool have_local_concordant = false;
      double c_local = 0.0;
      if (scanner) {
        // Count the three read categories at each side, classifying at the FINAL placed positions
        // (s1_pos/s2_pos) so the counts describe the reported breakpoint. The overlapping-mate
        // exclusion is referenced to those same placed positions (they are our best breakpoint
        // estimate now that placement is done).
        scanner->scan(s1_seq_id, s1_pos, s1_strand, s1_fwd, s2_tid, s2_pos, s2_fwd, distance_cutoff, s1_pos, s2_pos);
        int c1a = scanner->supporting(), c2a = scanner->concordant(), c3a = scanner->unpaired();
        map<string, int32_t> support_nums_1 = scanner->supporting_nums();   // copy before the next scan overwrites

        scanner->scan(s2_seq_id, s2_pos, s2_strand, s2_fwd, s1_tid, s1_pos, s1_fwd, distance_cutoff, s2_pos, s1_pos);
        int c1b = scanner->supporting(), c2b = scanner->concordant(), c3b = scanner->unpaired();
        const map<string, int32_t>& support_nums_2 = scanner->supporting_nums();

        // True support = read pairs whose BOTH mates qualify -- one at each side (intersect the sides'
        // supporting read-pair numbers). A single-side per-side count (c1a/c1b) over-counts at a
        // breakpoint SHARED with a neighboring junction: those reads' mates fall inside this junction's
        // wide (+/-D) partner window but land at the neighbor's breakpoint, so they appear on only one
        // side and drop out of the intersection. This matches exactly the pairs drawn in the joined plot.
        //
        // k_distinct counts the DISTINCT (outer_1, outer_2) fragment ends among those bridging pairs --
        // the DP analogue of JC's pos_hash_score. PCR/optical duplicates of a single original molecule
        // share both outer coordinates, so they collapse to one. Every statistical test below uses
        // k_distinct; k_support is kept for reporting and for downstream IS-copy tie-breaking.
        // (Strand adds nothing here: dp_classify_side_read already restricts each side to one strand.)
        set<pair<int32_t, int32_t> > distinct_ends;
        k_support = 0;
        for (map<string, int32_t>::const_iterator n = support_nums_1.begin(); n != support_nums_1.end(); n++) {
          map<string, int32_t>::const_iterator m = support_nums_2.find(n->first);
          if (m == support_nums_2.end()) continue;
          k_support++;
          distinct_ends.insert(make_pair(n->second, m->second));
        }
        k_distinct = static_cast<int>(distinct_ends.size());

        dp["discordant_count"] = to_string(k_support);
        dp["distinct_discordant_count"] = to_string(k_distinct);
        dp["side_1_discordant_count"] = to_string(c1a);   // raw per-side counts (may exceed the paired
        dp["side_2_discordant_count"] = to_string(c1b);   // count when a side is a shared-breakpoint hub)
        dp["side_1_concordant_count"] = to_string(c2a);
        dp["side_2_concordant_count"] = to_string(c2b);
        dp["side_1_unpaired_count"] = to_string(c3a);
        dp["side_2_unpaired_count"] = to_string(c3b);
        // Concordant pairs at a REDUNDANT side cross the intact reference junction of SOME copy of
        // the repeat, not necessarily this one, so they say nothing about whether THIS junction is
        // present -- averaging them in dilutes a junction whose unique side has no concordant support
        // at all. JC excludes such a side from its reference count the same way (resolve_alignments.cpp:
        // a side counts only when side_N_redundant != "1" and its annotate_key is not "repeat", and an
        // excluded side becomes "NA" and drops out of the mean). Average only the usable sides; when
        // both are redundant there is no denominator and the frequency is NA (written below).
        // The raw side_N_concordant_count fields above still report what was seen at each side.
        double c_sum = 0.0;
        int n_usable = 0;
        if (!side1_red) { c_sum += static_cast<double>(c2a); n_usable++; }
        if (!side2_red) { c_sum += static_cast<double>(c2b); n_usable++; }
        have_local_concordant = (n_usable > 0);
        c_local = have_local_concordant ? (c_sum / static_cast<double>(n_usable)) : 0.0;

        // Nothing survived verification: the candidate's shared-pair count came from the coarse
        // region-overlap heuristic, but at the placed breakpoint not one read pair actually bridges the
        // two sides. That is not weak evidence for a junction, it is the absence of any -- the frequency
        // is 0/c, which measures the reference being intact rather than a variant being rare. Emitting it
        // would put a "0 pairs, 0.0%" row in marginal.html and render a read-pair plot with nothing on it.
        // Typical cause is a repeat-driven mismapping pile-up: both sides redundant, huge unpaired counts.
        // Dropped rather than rejected, and counted in the run summary so the drop is not silent.
        if (k_distinct <= 0) {
          summary.discordant_pair.items_dropped_unsupported++;
          continue;
        }
      } else {
        // No BAM available: fall back to the heuristic count.
        dp["discordant_count"] = to_string(weight);
      }

      // How many pairs the spurious-pair background would be expected to place at this junction, across
      // every candidate junction tested. Reported for every item so the floor above can be audited.
      if (background.ok)
        dp["background_e_value"] = to_string(background.e_value(weight), 3, true);

      // Local frequency: discordant pairs vs concordant pairs spanning the SAME breakpoint, i.e. the
      // variant and reference observations of one sampling event. The non-redundant sides are averaged
      // because each measures the same fragment population from one end (see where c_local is built).
      // See also the comment above dp_crossing_mean.
      double f_lcb = std::numeric_limits<double>::quiet_NaN();
      if (have_local_concordant) {
        double k = static_cast<double>(k_distinct < 0 ? 0 : k_distinct);
        f_lcb = binomial_frequency_lower_bound(k, k + c_local, kDPFrequencyAlpha);
        // The shared FREQUENCY key, as every evidence type now uses. This used to carry its own
        // name so an evidence-level frequency could never be mistaken for a mutation's allele
        // frequency; that separation now lives where it belongs -- mutation_predictor attaches DP
        // items purely as supporting evidence matched by breakpoint and never reads a frequency
        // off one, so nothing can propagate this number to a mutation in the first place.
        dp[FREQUENCY] = to_string((k + c_local > 0.0) ? (k / (k + c_local)) : 0.0, 4, false);
        dp["concordant_count"] = to_string(c_local, 1, false);
        dp[FREQUENCY_LOWER] = to_string(f_lcb, 4, false);
        // Display only -- the gate below still tests the LOWER bound alone, since a high discordant
        // fraction is the signal rather than a problem. The upper bound exists so the report can
        // show an interval instead of a naked point estimate.
        dp[FREQUENCY_UPPER] = to_string(binomial_frequency_upper_bound(k, k + c_local, kDPFrequencyAlpha), 4, false);
      } else if (scanner) {
        // Both sides redundant: every concordant pair seen belongs to some copy of a repeat, so there
        // is no reference observation of THIS breakpoint to divide by -- the frequency is unknown, not
        // zero and not 100%. Written as "NA" (freq_to_string and freq_range_to_string already render
        // that literal), and f_lcb stays NaN so the frequency gate below skips the item entirely.
        dp[FREQUENCY] = "NA";
        dp["concordant_count"] = "NA";
        dp[FREQUENCY_LOWER] = "NA";
        dp[FREQUENCY_UPPER] = "NA";
      }

      // Discordance "skew" score: -log10 P(a normal position on side_1's seq_id is spanned by <= k
      // concordant pairs), projecting the run-wide reference distribution to this seq_id's coverage.
      double avgcov = dp_seq_coverage(summary, s1_seq_id);
      double dp_score = std::numeric_limits<double>::quiet_NaN();
      if (have_crossing) {
        dp_score = dp_discordance_skew(crossing_hist_ref, crossing_lo, crossing_hi, crossing_C_ref, avgcov, static_cast<uint32_t>(k_distinct < 0 ? 0 : k_distinct),
                                       crossing_use_empirical, crossing_nb_size, crossing_nb_mu);
      }
      dp[NEG_LOG10_DISCORDANCE_P_VALUE] = std::isnan(dp_score) ? string("NT") : to_string(dp_score, 1, false);

      // Expected concordant pairs spanning a normal position HERE: the uncensored reference mean scaled
      // to this sequence's coverage. This is the skew test's own effect size, and when it is small the
      // test cannot separate a real junction from noise -- for a Poisson null with mean lambda even k=0
      // only reaches -log10 P = lambda/ln(10), so below roughly 10 the skew can never clear a cutoff of
      // 3 on merit. Any rejection it produces there comes from the shape of the censored null rather
      // than from the data, so the test is reported but not allowed to reject.
      double lambda_x = (crossing_C_ref > 0.0) ? crossing_mean_ref * (avgcov / crossing_C_ref) : 0.0;
      bool skew_has_power = (settings.discordant_pair_minimum_crossing <= 0.0)
                            || (lambda_x >= settings.discordant_pair_minimum_crossing);
      dp["expected_concordant_count"] = to_string(lambda_x, 1, false);

      // Reject a junction whose local variant frequency is confidently below the cutoff: many concordant
      // pairs still span the breakpoint, so the reference is intact there and the discordant pairs are
      // not reporting a real rearrangement. Rejected DP items move to marginal.html (0 = cutoff off).
      if (!std::isnan(f_lcb) && settings.discordant_pair_frequency_cutoff > 0.0
          && f_lcb < settings.discordant_pair_frequency_cutoff) {
        // Its own reason string rather than the shared FREQUENCY_CUTOFF: what is compared against the
        // cutoff is the lower CONFIDENCE BOUND, not the frequency shown in the table, and the generic
        // "frequency below cutoff" wording makes a 25% item rejected on a 14% bound look like a bug.
        dp.add_reject_reason("DISCORDANT_PAIR_FREQUENCY");
      }

      // The original skew test, now only applied where it has power (see above).
      if (skew_has_power && !std::isnan(dp_score) && settings.discordant_pair_skew_cutoff > 0.0
          && dp_score > settings.discordant_pair_skew_cutoff) {
        dp.add_reject_reason("CONCORDANT_PAIR_SKEW");
      }

      // Tally the outcome for the run summary (an item can carry both reject reasons).
      {
        DiscordantPairSummary& dps = summary.discordant_pair;
        dps.items_tested++;
        if (circular_dp) dps.items_ignored_circular++;
        vector<string> reasons = dp.get_reject_reasons();
        bool by_freq = false, by_skew = false;
        for (size_t i = 0; i < reasons.size(); i++) {
          if (reasons[i] == "DISCORDANT_PAIR_FREQUENCY") by_freq = true;
          if (reasons[i] == "CONCORDANT_PAIR_SKEW") by_skew = true;
        }
        if (by_freq) dps.items_rejected_frequency++;
        if (by_skew) dps.items_rejected_skew++;
        if (!by_freq && !by_skew && !circular_dp) dps.items_accepted++;
      }

      emitted_by_breakpoint[breakpoint_key] = dp_gd.add(dp);
    }

    if (scanner) delete scanner;

    dp_gd.write(settings.dp_genome_diff_file_name);
  }

  // -------------------------------------------------------------------------------------------------
  // Per-side read-pair plots
  // -------------------------------------------------------------------------------------------------

  // Format an integer with thousands separators, e.g. 1234560 -> "1,234,560".
  string plot_commafy(int64_t v)
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
  int64_t plot_nice_tick(int64_t span)
  {
    double raw = span / 10.0;
    static const int64_t nice[] = {10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000};
    for (size_t i = 0; i < sizeof(nice)/sizeof(nice[0]); i++)
      if (nice[i] >= raw) return nice[i];
    return 100000;
  }

  // Left margin, in base-font character widths, so the leftmost x tick label is not clipped.
  //
  // These plots carry no y-axis labels (read names sit on the RIGHT via y2tics), so gnuplot leaves the
  // left axis almost no margin. But an x tick label is CENTERED on its tick and drawn at a larger font
  // than the terminal's base font, so about half of the leftmost label hangs past the plot edge and is
  // cut off by the canvas. Reserve half a label's width, converted from label-font to base-font
  // character units (the average-character-width factor cancels in the ratio), plus a character of slack.
  int plot_lmargin_for_labels(size_t max_label_chars, double label_font, double base_font)
  {
    if (!(base_font > 0.0)) return 8;
    double chars = ceil((static_cast<double>(max_label_chars) * label_font / 2.0) / base_font);
    return static_cast<int>(chars) + 1;
  }

  // Build a gnuplot explicit-tics list ("label" pos, ...) at multiples of `step` within [lo, hi],
  // with comma-formatted whole-number labels.
  string plot_xtics_list(int64_t lo, int64_t hi, int64_t step)
  {
    string out;
    int64_t first = ((lo + step - 1) / step) * step;  // first multiple of step >= lo
    bool comma = false;
    for (int64_t t = first; t <= hi; t += step) {
      if (comma) out += ", ";
      out += double_quote(plot_commafy(t)) + " " + to_string(t);
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
    int64_t step = plot_nice_tick(xmax - xmin);
    string tics = plot_xtics_list(xmin, xmax, step);
    s << "set xtics (" << tics << ") nomirror font ',16'" << endl;
    s << "set x2tics (" << tics << ") font ',16'" << endl;
    // Keep the widest label (the extremes of the range) inside the canvas.
    s << "set lmargin " << plot_lmargin_for_labels(max(plot_commafy(xmin).size(), plot_commafy(xmax).size()), 16.0, 11.0) << endl;

    // Read-pair names as lane labels down the right side. For a concordant pair (both mates drawn) the
    // label shows both read names; a discordant lane shows the single in-window read.
    // Unlabeled tick marks at each read lane on the left axis; the names stay on the right only.
    s << "set ytics 1 nomirror format ''" << endl;
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
      int64_t step = plot_nice_tick(xmax - xmin);
      ostringstream tks;
      bool first = true;
      // Longest single ROW of a label -- these labels are two-row (embedded "\n"), and only the wider
      // row governs how far the centered label sticks out past the plot edge.
      size_t widest_row = 0;
      auto emit = [&](int64_t x, const string& label) {
        if (!first) tks << ", "; first = false;
        tks << double_quote(label) << " " << x;
        for (size_t b = 0; b <= label.size(); ) {
          size_t e = label.find("\\n", b);
          size_t len = (e == string::npos ? label.size() : e) - b;
          if (len > widest_row) widest_row = len;
          if (e == string::npos) break;
          b = e + 2;
        }
      };
      // Seam: side_1_position (upper) over side_2_position (lower).
      emit(0, plot_commafy(p1) + "\\n" + plot_commafy(p2));
      // Left half = side_1, upper row (trailing newline leaves the lower row blank).
      for (int64_t x = -step; x >= xmin; x -= step) {
        int64_t g = (s1_str == -1) ? (p1 + x) : (p1 - x);
        emit(x, plot_commafy(g) + "\\n");
      }
      // Right half = side_2, lower row (leading newline leaves the upper row blank).
      for (int64_t x = step; x <= xmax; x += step) {
        int64_t g = (s2_str == -1) ? (p2 - x) : (p2 + x);
        emit(x, "\\n" + plot_commafy(g));
      }
      s << "set xtics ("  << tks.str() << ") nomirror font ',16'" << endl;  // bottom axis
      s << "set x2tics (" << tks.str() << ") font ',16'" << endl;           // top axis (same two rows)
      s << "set lmargin " << plot_lmargin_for_labels(widest_row, 16.0, 11.0) << endl;
    }

    // Read-pair names down the right side (both mates).
    // Unlabeled tick marks at each read lane on the left axis; the names stay on the right only.
    s << "set ytics 1 nomirror format ''" << endl;
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
