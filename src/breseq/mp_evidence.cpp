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

#include "mp_evidence.h"
#include "dp_evidence.h"   // paired_library_params / paired_region_to_side (shared pair geometry)
#include "genome_diff.h"
#include "genome_diff_entry.h"
#include "pileup.h"   // pileup_base + alignment_wrapper + BAM flag macros
#include "stats.h"    // binomial_frequency_lower_bound / _upper_bound

#include <set>
#include <random>
#include <algorithm>

using namespace std;

namespace breseq {

  // One-sided confidence level for the MP local-frequency (Clopper-Pearson) bounds. Matches DP.
  static const double kMPFrequencyAlpha = 0.05;

  // ---------------------------------------------------------------------------------------------
  // Candidate regions
  // ---------------------------------------------------------------------------------------------

  // One candidate region parsed from MP_candidate_regions.csv.
  struct mp_region_row {
    string   seq_id;
    uint32_t start;      // lower coordinate of the region span
    uint32_t end;        // higher coordinate
    char     strand;     // focal-read strand: 'F' or 'R'
    uint32_t max_count;  // peak in-window count while the region was open (seed diagnostic)
    uint32_t max_total;  // same-strand in-window total when that peak was reached (seed diagnostic)
    bool     redundant;  // majority of this region's reads mapped redundantly
    int32_t  extreme_read_end;    // max(read_end) over the region's reads
    int32_t  extreme_read_start;  // min(read_start) over the region's reads
    bool     have_extents;
  };

  // ---------------------------------------------------------------------------------------------
  // Exact mate extents
  // ---------------------------------------------------------------------------------------------

  // Read number (the part after ':') of a breseq QNAME "<file_index>:<num>". Both mates of a pair
  // share it, so it is how a read finds its mate's alignment.
  static string mp_read_num(const string& name)
  {
    size_t c = name.find(':');
    return (c == string::npos) ? name : name.substr(c + 1);
  }

  // read number -> the ALIGNED extents of that pair's alignments near the position of interest.
  typedef map<string, vector<pair<int32_t, int32_t> > > mp_extent_index;

  // Indexes aligned read extents over a window, so a read can look up its mate's exact aligned end.
  //
  // Nothing on a read gives that end directly. TLEN is no help: breseq computes the pair distance
  // from UNCLIPPED extents, so a mate that aligns flush to a breakpoint and soft-clips the rest still
  // reports a fragment reaching past it -- which is exactly backwards, since those clipped bases are
  // the novel sequence. Assuming a full-length mate is worse for the same reason. Only the mate's own
  // CIGAR settles it, so it has to be looked up.
  //
  // The window is local (a few fragment lengths around one candidate) rather than genome-wide, which
  // is what keeps this cheap: a proper pair's mate is by definition within a fragment length, so a
  // small window around the position always contains it. DP indexes the whole BAM instead, but it can
  // afford to -- it only needs the 1.2% of alignments that are discordant, whereas MP would have to
  // index every read.
  class mp_extent_indexer : public pileup_base {
  public:
    mp_extent_indexer(const string& bam, const string& fasta)
      : pileup_base(bam, fasta), m_out(NULL) { set_print_progress(false); }

    int32_t tid_for_seq_id(const string& seq_id) const {
      for (uint32_t t = 0; t < num_targets(); t++)
        if (seq_id == string(target_name(t))) return static_cast<int32_t>(t);
      return -1;
    }

    void build(const string& seq_id, int32_t lo, int32_t hi, mp_extent_index& out)
    {
      out.clear();
      int32_t tid = tid_for_seq_id(seq_id);
      if (tid < 0) return;
      lo = max(1, lo);
      hi = min(static_cast<int32_t>(target_length(tid)), hi);
      if (lo > hi) return;
      m_out = &out;
      do_fetch(seq_id + ":" + to_string(lo) + "-" + to_string(hi));
      m_out = NULL;
    }

    void fetch_callback(const alignment_wrapper& a)
    {
      if (m_out == NULL) return;
      if (a.flag() & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) return;
      if (a.unmapped()) return;
      (*m_out)[mp_read_num(a.read_name())].push_back(
        make_pair(static_cast<int32_t>(a.reference_start_1()), static_cast<int32_t>(a.reference_end_1())));
    }

  private:
    mp_extent_index* m_out;
  };

  // Look up the mate's exact aligned extent. Returns false when it is not in the indexed window,
  // in which case the caller treats the pair as uninformative rather than guessing.
  static bool mp_mate_extent(const alignment_wrapper& a, const mp_extent_index* index,
                             int32_t& mlo, int32_t& mhi)
  {
    if (index == NULL) return false;
    mp_extent_index::const_iterator it = index->find(mp_read_num(a.read_name()));
    if (it == index->end()) return false;
    int32_t mpos = a.mate_start_1();
    for (size_t i = 0; i < it->second.size(); i++) {
      if (it->second[i].first == mpos) { mlo = it->second[i].first; mhi = it->second[i].second; return true; }
    }
    return false;
  }

  // ---------------------------------------------------------------------------------------------
  // Side scanner
  // ---------------------------------------------------------------------------------------------

  // Why MP cannot reuse dp_side_scanner: dp_classify_side_read requires a supporting read's MATE to
  // land at a second locus in its own crossing orientation. An MP read has no mapped mate at all, so
  // that classifier rejects every read MP cares about (it folds them into its "unpaired" tally and
  // stops). The window arithmetic below is the same, but the classification is MP's own.
  class mp_side_scanner : public pileup_base {
  public:
    mp_side_scanner(const string& bam, const string& fasta)
      : pileup_base(bam, fasta), m_p(0), m_s(0), m_cross_fwd(true),
        m_window(0.0), m_collect(false), m_ext(NULL), m_supporting(0), m_spanning(0), m_window_total(0),
        m_have_inner(false), m_inner_edge(0)
      { set_print_progress(false); }

    //! Exact mate extents for the position about to be scanned. Not owned; must outlive the scan.
    void set_extent_index(const mp_extent_index* ext) { m_ext = ext; }

    int32_t tid_for_seq_id(const string& seq_id) const {
      for (uint32_t t = 0; t < num_targets(); t++)
        if (seq_id == string(target_name(t))) return static_cast<int32_t>(t);
      return -1;
    }
    int32_t seq_length(int32_t tid) const { return static_cast<int32_t>(target_length(tid)); }

    //! Fetch the supporting (mate-unmapped) reads on this side and return:
    //   median_outer  = median of their "outside" coordinates -- the end facing AWAY from the
    //                   insertion point (reference_start for s=-1, reference_end for s=+1). The
    //                   caller shifts this toward the insertion by half the median pair distance;
    //                   that is the same half-MAD-derived placement DP uses, and it works here for
    //                   the same reason: the unsequenced middle of each fragment is what separates a
    //                   read's outer end from the breakpoint, and its expected length is half the
    //                   median fragment length.
    //   inner_edge    = the furthest a supporting read reaches TOWARD the insertion (its aligned
    //                   junction-facing end). An aligned base cannot lie past the breakpoint, so the
    //                   caller clamps the placed coordinate out to this.
    //  Returns false (no placement possible) when no supporting read is found.
    bool supporting_outer_median(const string& seq_id, int32_t p, int32_t s, bool cross_fwd, double D,
                                 int32_t& median_outer, int32_t& inner_edge)
    {
      set_ctx(p, s, cross_fwd, 0.0);
      m_outside.clear();
      m_have_inner = false; m_inner_edge = 0;
      m_collect = true;

      fetch_side_window(seq_id, p, s, D);

      m_collect = false;
      if (m_outside.empty()) return false;
      sort(m_outside.begin(), m_outside.end());
      median_outer = m_outside[m_outside.size() / 2];
      inner_edge = m_inner_edge;
      return true;
    }

    //! Count the supporting and spanning-pair reads at a (final) position, along with the
    //  number of DISTINCT outer coordinates among the supporting reads -- the MP analogue of JC's
    //  pos_hash_score, so PCR duplicates of one molecule cannot carry a prediction. Counted here, in
    //  the SAME pass as m_supporting, so the two numbers always describe the same set of reads.
    //  W is the counting window: both the fetch half-width AND the bound on how far a read's outer
    //  end may lie from p. It must be the width the null was tabulated at -- see classify().
    void scan(const string& seq_id, int32_t p, int32_t s, bool cross_fwd, double W)
    {
      set_ctx(p, s, cross_fwd, W);
      m_supporting = 0; m_spanning = 0; m_window_total = 0;
      m_distinct.clear();
      m_collect = false;
      fetch_side_window(seq_id, p, s, W);
    }

    int supporting() const { return m_supporting; }
    int spanning() const { return m_spanning; }
    int window_total() const { return m_window_total; }
    size_t distinct() const { return m_distinct.size(); }

    void fetch_callback(const alignment_wrapper& a)
    {
      int32_t outer = 0, inner = 0;
      int cat = classify(a, outer, inner);
      if (cat == 0) return;
      if (m_collect) {
        if (cat == 1) {
          m_outside.push_back(outer);
          if (!m_have_inner || (m_s == -1 ? inner > m_inner_edge : inner < m_inner_edge)) {
            m_have_inner = true; m_inner_edge = inner;
          }
        }
        return;
      }
      // Every read that reached here is on the crossing strand with its body on the kept flank
      // inside the window: the opportunity set. The numerator is a subset of it by construction,
      // and the genome-wide null is tabulated over exactly this definition in the stage-08 pileup,
      // so k/m_window_total is the quantity the score's null describes.
      m_window_total++;
      if      (cat == 1) { m_supporting++; m_distinct.insert(outer); }
      else if (cat == 2) m_spanning++;
    }

  private:
    void set_ctx(int32_t p, int32_t s, bool cross_fwd, double window)
    { m_p = p; m_s = s; m_cross_fwd = cross_fwd; m_window = window; }

    // One-sided fetch window [lo, hi]: the kept flank, out to D.
    void fetch_side_window(const string& seq_id, int32_t p, int32_t s, double D)
    {
      int32_t tid = tid_for_seq_id(seq_id);
      if (tid < 0) return;
      int32_t seqlen = seq_length(tid);
      int32_t lo, hi;
      if (s == -1) { lo = max(1, static_cast<int32_t>(p - D)); hi = min(seqlen, p); }
      else         { lo = max(1, p);                           hi = min(seqlen, static_cast<int32_t>(p + D)); }
      if (lo > hi) return;
      do_fetch(seq_id + ":" + to_string(lo) + "-" + to_string(hi));
    }

    // 0 = ignore, 1 = supporting (mate unmapped), 2 = spanning pair (the denominator),
    // 3 = present but uninformative (a pair that never reached the breakpoint).
    int classify(const alignment_wrapper& a, int32_t& outer, int32_t& inner) const
    {
      outer = 0; inner = 0;
      if (a.flag() & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) return 0;
      if (a.unmapped()) return 0;

      // Only reads on this side's crossing strand -- the ones pointing INTO the insertion.
      if ((!a.reversed()) != m_cross_fwd) return 0;

      int32_t rs = static_cast<int32_t>(a.reference_start_1());
      int32_t re = static_cast<int32_t>(a.reference_end_1());

      // The read's body must be on the kept side.
      if (m_s == -1) { if (rs > m_p) return 0; outer = rs; inner = re; }
      else           { if (re < m_p) return 0; outer = re; inner = rs; }

      // ...and its OUTER end must lie within one window of the placed position. Two reasons, and
      // they agree:
      //
      //  - Geometrically, the outer end is what separates this read from the breakpoint. A read
      //    further out than one fragment length has no mate that could have reached p, so it says
      //    nothing about an insertion there either way.
      //  - Statistically, this is what makes the opportunity set exactly m_window wide. do_fetch
      //    returns everything OVERLAPPING the window, so without this the set spans
      //    m_window + read_length in outer-end coordinates while the null was tabulated over a
      //    sliding window of exactly m_window. Calibrating on one width and testing on another
      //    misstates n, and rho's variance inflation (1 + (n-1)*rho) depends on n.
      if (m_window > 0.0 && abs(outer - m_p) > m_window) { outer = 0; inner = 0; return 0; }

      // Supporting: paired, mapped, and the mate produced NO alignment at all.
      //
      // Both conditions are needed. BAM_FMUNMAP is also set for a mate the aligner placed but that
      // breseq then rejected in test_read_alignment_requirements -- most often for falling short of
      // --require-match-fraction, which --predict-soft-clipping silently lowers from 0.9 to 0.5.
      // Counting those made this numerator a measure of alignment stringency: the same sample
      // called an order of magnitude more MP evidence with soft clipping off than on, and the extra
      // calls were partially-aligning mates (SC's signal), not sequence missing from the reference.
      // A mate that lands wholly inside a novel insert aligns 0% of its length, not 50-90%.
      if (a.is_paired() && (a.flag() & BAM_FMUNMAP) && mate_never_aligned(a)) return 1;

      // Denominator: ONLY a pair whose mate actually maps past the breakpoint, into the region the
      // insertion would occupy. That molecule demonstrably carries reference sequence where this
      // prediction says novel sequence sits, so it is evidence against -- and it is the only kind of
      // pair that is. A pair whose mate lands back on the same flank never reached the breakpoint and
      // cannot speak to it either way; counting those was wrong, and it dragged the frequency of a
      // homozygous insertion down towards the fraction of local fragments that happen to be short.
      //
      // "Reaches past p" is decided from the mate's own ALIGNED extent (see mp_extent_indexer for why
      // neither TLEN nor the read length will do). At a real breakpoint the mate is precisely the read
      // that soft-clips there, so any full-length assumption pushes it past p and turns the strongest
      // supporting reads into evidence against themselves.
      if (!a.is_paired() || !a.proper_pair()) return 3;
      if (static_cast<int32_t>(a.mate_reference_target_id()) != static_cast<int32_t>(a.reference_target_id())) return 3;
      int32_t mlo = 0, mhi = 0;
      if (!mp_mate_extent(a, m_ext, mlo, mhi)) return 3;
      bool spans = (m_s == -1) ? (mhi > m_p) : (mlo < m_p);
      return spans ? 2 : 3;
    }

    int32_t m_p;
    int32_t m_s;
    bool    m_cross_fwd;
    //! Half-width of the opportunity set, in outer-end coordinates. 0 disables the bound, which is
    //! what the placement pass wants: it scans out to D looking for ANY supporting read to place on.
    double  m_window;
    bool    m_collect;
    const mp_extent_index* m_ext;
    int     m_supporting, m_spanning, m_window_total;
    set<int32_t> m_distinct;
    vector<int32_t> m_outside;
    bool    m_have_inner;
    int32_t m_inner_edge;
  };

  // ---------------------------------------------------------------------------------------------
  // Evidence plot
  // ---------------------------------------------------------------------------------------------

  // One read drawn on an MP evidence plot.
  struct mp_draw_read {
    int      category;      // 1 = supporting, 2 = spanning pair, 3 = pair that never reached p
    int32_t  read_start, read_end;
    bool     read_reversed;
    int32_t  anchor;        // this read's insertion-facing coordinate
    int32_t  mate_start;    // categories 2/3 only
    int32_t  mate_end;      // categories 2/3 only (from TLEN where it pins the fragment end)
    bool     mate_reversed; // categories 2/3 only
    string   name;          // QNAME, shown as the lane label
  };

  // Gathers the reads to draw at an MP position in a separate fetch pass, using the SAME
  // classification as mp_side_scanner so the plot shows exactly the reads that were counted.
  class mp_plot_gatherer : public pileup_base {
  public:
    mp_plot_gatherer(const string& bam, const string& fasta)
      : pileup_base(bam, fasta), m_p(0), m_s(0), m_cross_fwd(true), m_window(0.0), m_ext(NULL)
      { set_print_progress(false); }

    //! Exact mate extents for the position about to be gathered. Not owned; must outlive the gather.
    void set_extent_index(const mp_extent_index* ext) { m_ext = ext; }

    int32_t tid_for_seq_id(const string& seq_id) const {
      for (uint32_t t = 0; t < num_targets(); t++)
        if (seq_id == string(target_name(t))) return static_cast<int32_t>(t);
      return -1;
    }

    void gather(const string& seq_id, int32_t p, int32_t s, bool cross_fwd, double W)
    {
      m_p = p; m_s = s; m_cross_fwd = cross_fwd; m_window = W;
      m_reads.clear();
      int32_t tid = tid_for_seq_id(seq_id);
      if (tid < 0) return;
      int32_t seqlen = static_cast<int32_t>(target_length(tid));
      int32_t lo, hi;
      if (s == -1) { lo = max(1, static_cast<int32_t>(p - W)); hi = min(seqlen, p); }
      else         { lo = max(1, p);                           hi = min(seqlen, static_cast<int32_t>(p + W)); }
      if (lo > hi) return;
      do_fetch(seq_id + ":" + to_string(lo) + "-" + to_string(hi));
    }

    const vector<mp_draw_read>& reads() const { return m_reads; }

    void fetch_callback(const alignment_wrapper& a)
    {
      if (a.flag() & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) return;
      if (a.unmapped()) return;
      if ((!a.reversed()) != m_cross_fwd) return;

      int32_t rs = static_cast<int32_t>(a.reference_start_1());
      int32_t re = static_cast<int32_t>(a.reference_end_1());
      if (m_s == -1) { if (rs > m_p) return; } else { if (re < m_p) return; }

      // Same outer-end bound mp_side_scanner::classify applies, so the plot does not show reads that
      // were outside the counting window and therefore counted nowhere.
      int32_t outer = (m_s == -1) ? rs : re;
      if (m_window > 0.0 && abs(outer - m_p) > m_window) return;

      mp_draw_read r;
      r.read_start = rs;
      r.read_end = re;
      r.read_reversed = a.reversed();
      r.anchor = (m_s == -1) ? re : rs;
      r.mate_start = 0;
      r.mate_end = 0;
      r.mate_reversed = false;
      r.name = a.read_name();

      // Same three-way split as mp_side_scanner::classify, so the plot shows exactly what was counted.
      bool mate_unmapped = a.is_paired() && ((a.flag() & BAM_FMUNMAP) != 0) && mate_never_aligned(a);
      if (mate_unmapped) {
        r.category = 1;
      } else if (!a.is_paired() || !a.proper_pair() ||
                 (static_cast<int32_t>(a.mate_reference_target_id()) != static_cast<int32_t>(a.reference_target_id()))) {
        r.category = 3;
        if (a.is_paired()) {
          r.mate_start = a.mate_start_1();
          r.mate_end = r.mate_start + static_cast<int32_t>(a.read_length()) - 1;
          r.mate_reversed = (a.flag() & BAM_FMREVERSE) != 0;
        }
      } else {
        // Exact mate extent, exactly as in mp_side_scanner::classify -- so the arrow is drawn to where
        // the mate really aligns, and the color matches the category it was counted as.
        int32_t mlo = 0, mhi = 0;
        bool have = mp_mate_extent(a, m_ext, mlo, mhi);
        bool spans = have && ((m_s == -1) ? (mhi > m_p) : (mlo < m_p));
        r.category = spans ? 2 : 3;
        r.mate_start = have ? mlo : a.mate_start_1();
        r.mate_end = have ? mhi : (r.mate_start + static_cast<int32_t>(a.read_length()) - 1);
        r.mate_reversed = (a.flag() & BAM_FMREVERSE) != 0;
      }
      m_reads.push_back(r);
    }

  private:
    int32_t m_p;
    int32_t m_s;
    bool    m_cross_fwd;
    double  m_window;
    const mp_extent_index* m_ext;
    vector<mp_draw_read> m_reads;
  };

  static const string MP_UNPAIRED_COLOR = "'#7f3fbf'";  // mate unmapped        = purple (matches the .PM read tag)
  static const string MP_PAIRED_COLOR   = "'#1f77b4'";  // mate maps past p     = blue (matches DP's concordant)
  static const string MP_IGNORED_COLOR  = "'#b0b0b0'";  // pair never reached p = gray, shown but not counted
  static const string MP_KEPT_COLOR     = "'#efe6ff'";  // light purple over the flank the reads sit on

  // Render the MP evidence plot (SVG) via gnuplot: one read per lane, drawn as an arrow in its mapping
  // direction, colored by what it says about the prediction:
  //   purple  a read whose mate did not map anywhere -- the supporting evidence. Drawn as the arrow
  //           ALONE: there is no second alignment, so there is no extent to connect it to, and a
  //           connector would assert a mate position that was never observed.
  //   blue    a pair whose mate maps PAST the placed position, into the region the insertion would
  //           occupy. That molecule carries reference sequence there, so it counts against -- drawn
  //           with its mate and a dashed connector, the intact molecule visible end to end.
  //   gray    a pair that never reached the placed position. Shown for context (it is why the locus
  //           looks covered) but counted neither way, because it cannot speak to the prediction.
  // An insertion therefore reads as bare purple arrows crowding the shaded edge with no blue at all.
  static void render_mp_plot(const string& output_svg, vector<mp_draw_read> reads,
                             int32_t p, int32_t strand, const map<int,int>& mate_prefix)
  {
    sort(reads.begin(), reads.end(), [](const mp_draw_read& a, const mp_draw_read& b) {
      if (a.read_start != b.read_start) return a.read_start < b.read_start;
      return a.read_end < b.read_end;
    });

    int n = static_cast<int>(reads.size());

    // Symmetric window about p, sized ~1.1x the furthest extent of any drawn read (or mate) from p.
    int64_t lo = p, hi = p;
    for (int i = 0; i < n; i++) {
      const mp_draw_read& r = reads[i];
      lo = min(lo, (int64_t)r.read_start); hi = max(hi, (int64_t)r.read_end);
      if (r.category == 2 && r.mate_start > 0) {
        lo = min(lo, (int64_t)r.mate_start); hi = max(hi, (int64_t)r.mate_end);
      }
    }
    int64_t reach = max((int64_t)(p - lo), (int64_t)(hi - p));
    if (reach < 10) reach = 10;
    int64_t half = (int64_t)(1.1 * reach + 0.5);
    int64_t xmin = p - half;
    int64_t xmax = p + half;

    string f_pr = output_svg + ".pr.tab";  // spanning-pair read arrows (x y dx dy)
    string f_ur = output_svg + ".ur.tab";  // mate-unmapped read arrows
    string f_xr = output_svg + ".xr.tab";  // uninformative-pair read arrows
    string f_pc = output_svg + ".pc.tab";  // spanning-pair connectors
    string f_xc = output_svg + ".xc.tab";  // uninformative-pair connectors
    ofstream pr(f_pr.c_str()), ur(f_ur.c_str()), xr(f_xr.c_str()), pc(f_pc.c_str()), xc(f_xc.c_str());
    bool has_pr=false, has_ur=false, has_xr=false, has_pc=false, has_xc=false;

    for (int i = 0; i < n; i++) {
      const mp_draw_read& r = reads[i];
      int y = i + 1;

      int32_t tail = r.read_reversed ? r.read_end   : r.read_start;
      int32_t head = r.read_reversed ? r.read_start : r.read_end;
      if      (r.category == 2) { pr << tail << "\t" << y << "\t" << (head - tail) << "\t0\n"; has_pr = true; }
      else if (r.category == 3) { xr << tail << "\t" << y << "\t" << (head - tail) << "\t0\n"; has_xr = true; }
      else                      { ur << tail << "\t" << y << "\t" << (head - tail) << "\t0\n"; has_ur = true; }

      if ((r.category == 2 || r.category == 3) && r.mate_start > 0) {
        int32_t m_tail = r.mate_reversed ? r.mate_end   : r.mate_start;
        int32_t m_head = r.mate_reversed ? r.mate_start : r.mate_end;
        if (r.category == 2) {
          pr << m_tail << "\t" << y << "\t" << (m_head - m_tail) << "\t0\n"; has_pr = true;
          pc << r.anchor << "\t" << y << "\n" << r.mate_start << "\t" << y << "\n\n"; has_pc = true;
        } else {
          xr << m_tail << "\t" << y << "\t" << (m_head - m_tail) << "\t0\n"; has_xr = true;
          xc << r.anchor << "\t" << y << "\n" << r.mate_start << "\t" << y << "\n\n"; has_xc = true;
        }
      }
      // A supporting read gets no connector: its mate has no alignment, so there is nothing to
      // connect to and no observed coordinate to draw toward.
    }
    pr.close(); ur.close(); xr.close(); pc.close(); xc.close();

    ostringstream s;
    // Space reserved below the plot area, in pixels: x tics, then the x label, then a clear gap, then
    // the key. Pinning the plot bottom and the key to explicit screen fractions (rather than letting
    // "key outside bottom" pack them together) is what keeps that gap from collapsing onto the label,
    // and keeps it constant as the canvas grows with the number of read lanes.
    const int kBelowPlotPx = 145;   // tics + label + gap + key
    const int kKeyBaselinePx = 30;  // key baseline above the bottom edge
    int height = max(300, n * 8 + kBelowPlotPx + 30);
    s << "set terminal svg size 1400," << height << " font ',11' noenhanced" << endl;
    s << "set output " << double_quote(output_svg) << endl;
    s << "unset title" << endl;
    s << "set bmargin at screen " << (static_cast<double>(kBelowPlotPx) / height) << endl;
    s << "set xlabel 'Reference position (bp)' font ',22'" << endl;
    s << "set xrange [" << xmin << ":" << xmax << "]" << endl;
    s << "set x2range [" << xmin << ":" << xmax << "]" << endl;
    s << "set yrange [0:" << (n + 1) << "]" << endl;
    s << "set y2range [0:" << (n + 1) << "]" << endl;
    s << "set border 15 lw 1" << endl;
    s << "set tics out" << endl;

    int64_t step = plot_nice_tick(xmax - xmin);
    string tics = plot_xtics_list(xmin, xmax, step);
    s << "set xtics (" << tics << ") nomirror font ',16'" << endl;
    s << "set x2tics (" << tics << ") font ',16'" << endl;
    s << "set lmargin " << plot_lmargin_for_labels(max(plot_commafy(xmin).size(), plot_commafy(xmax).size()), 16.0, 11.0) << endl;

    s << "set ytics 1 nomirror format ''" << endl;
    {
      ostringstream y2;
      for (int i = 0; i < n; i++) {
        if (i) y2 << ", ";
        string label = reads[i].name;
        if (reads[i].category == 2 && reads[i].mate_start > 0) {
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

    // Shade the half beyond the placed position -- where the missing mates went -- gray, and the flank
    // the reads sit on light purple.
    int64_t gray_lo, gray_hi, keep_lo, keep_hi;
    if (strand == -1) { gray_lo = p; gray_hi = xmax; keep_lo = xmin; keep_hi = p; }
    else              { gray_lo = xmin; gray_hi = p; keep_lo = p; keep_hi = xmax; }
    s << "set object 1 rectangle from " << gray_lo << ",0 to " << gray_hi << "," << (n + 1)
      << " fc rgb 'gray90' fs solid noborder behind" << endl;
    s << "set object 2 rectangle from " << keep_lo << ",0 to " << keep_hi << "," << (n + 1)
      << " fc rgb " << MP_KEPT_COLOR << " fs solid noborder behind" << endl;

    s << "set arrow from " << p << ",0 to " << p << "," << (n + 1) << " nohead lc rgb 'gray50' lw 1 dt 3 back" << endl;

    // Legend, drawn by hand as labels + arrows rather than with "set key". The two colors are the whole
    // point of this plot, so they have to be named -- but gnuplot centers a key on the PLOT AREA, which
    // here is pushed right by the read-name labels down the right edge and an explicit lmargin on the
    // left, leaving the legend visibly off-center whichever way it is anchored. Laying it out in screen
    // pixels centers it on the canvas exactly. Entries are placed only for the categories actually
    // drawn, so a plot with just one kind of read gets a single centered entry.
    {
      const double kCanvasW = 1400.0;
      const double kCharW   = 8.6;   // ~average glyph width at font size 18 in the svg terminal
      const double kSample  = 44.0;  // legend arrow length
      const double kGap     = 10.0;  // arrow -> text
      const double kBetween = 70.0;  // entry -> entry

      vector<pair<string, string> > entries;  // (color, text)
      if (has_ur) entries.push_back(make_pair(MP_UNPAIRED_COLOR, string("mate did not map (supports MP)")));
      if (has_pr) entries.push_back(make_pair(MP_PAIRED_COLOR,   string("mate maps past the breakpoint")));
      if (has_xr) entries.push_back(make_pair(MP_IGNORED_COLOR,  string("pair does not span (not counted)")));

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

    // Uninformative pairs first and thinner, so the counted reads read on top of them.
    vector<string> clauses;
    if (has_xc) clauses.push_back(double_quote(f_xc) + " using 1:2 with lines dt 2 lc rgb " + MP_IGNORED_COLOR + " lw 1 notitle");
    if (has_pc) clauses.push_back(double_quote(f_pc) + " using 1:2 with lines dt 2 lc rgb " + MP_PAIRED_COLOR + " lw 1 notitle");
    if (has_xr) clauses.push_back(double_quote(f_xr) + " using 1:2:3:4 with vectors head filled size screen 0.006,20,60 lc rgb " + MP_IGNORED_COLOR + " lw 1 notitle");
    if (has_pr) clauses.push_back(double_quote(f_pr) + " using 1:2:3:4 with vectors head filled size screen 0.008,20,60 lc rgb " + MP_PAIRED_COLOR + " lw 2 notitle");
    if (has_ur) clauses.push_back(double_quote(f_ur) + " using 1:2:3:4 with vectors head filled size screen 0.008,20,60 lc rgb " + MP_UNPAIRED_COLOR + " lw 2 notitle");
    if (clauses.empty()) clauses.push_back("-1 notitle");
    s << "plot " << join(clauses, string(", \\\n     ")) << endl;

    string script_name = output_svg + ".gp";
    string log_name    = output_svg + ".gp.log";
    run_gnuplot_script(s.str(), script_name, log_name);
    remove(log_name.c_str());
    remove(f_pr.c_str()); remove(f_ur.c_str()); remove(f_xr.c_str());
    remove(f_pc.c_str()); remove(f_xc.c_str());
  }

  // If v holds more than max_display entries, randomly down-sample it (seeded, so runs reproduce) and
  // return the ORIGINAL size; otherwise leave v alone and return 0. Mirrors dp_cap_for_display.
  static size_t mp_cap_for_display(vector<mp_draw_read>& v, uint32_t max_display, uint32_t seed)
  {
    size_t total = v.size();
    if (max_display == 0 || total <= max_display) return 0;
    std::mt19937 rng(seed);
    std::shuffle(v.begin(), v.end(), rng);
    v.resize(max_display);
    return total;
  }

  void draw_missing_pair_evidence_plots(const Settings& settings, Summary& summary,
                                        cReferenceSequences& ref_seq_info, cGenomeDiff& gd)
  {
    (void)ref_seq_info;

    bool inner3p = true;
    double D = 0.0, pair_median = 0.0;
    if (!paired_library_params(summary, inner3p, D, pair_median)) return;  // predict already warned
    if (!file_exists(settings.reference_bam_file_name.c_str()) ||
        !file_exists(settings.reference_fasta_file_name.c_str())) return;

    diff_entry_list_t mp_list = gd.get_list(make_vector<gd_entry_type>(MP));
    if (mp_list.empty()) return;

    create_path(settings.evidence_path);
    mp_plot_gatherer g(settings.reference_bam_file_name, settings.reference_fasta_file_name);
    mp_extent_indexer indexer(settings.reference_bam_file_name, settings.reference_fasta_file_name);

    // Map a read-name file prefix (m_id + 1) to its mate's, so a mate-mapped lane shows both names.
    map<int,int> mate_prefix;
    for (cReadFileSets::const_iterator rfs = settings.read_file_sets.begin(); rfs != settings.read_file_sets.end(); rfs++) {
      if (rfs->is_paired()) {
        int a = static_cast<int>(rfs->m_files[0].m_id) + 1;
        int b = static_cast<int>(rfs->m_files[1].m_id) + 1;
        mate_prefix[a] = b; mate_prefix[b] = a;
      }
    }

    // Outside the loop: the gatherer holds a pointer to this, so it must outlive every use. See the
    // note in predict_missing_pairs -- the same object declared inside the loop there was a
    // use-after-free.
    mp_extent_index ext;

    for (diff_entry_list_t::iterator it = mp_list.begin(); it != mp_list.end(); it++) {
      cDiffEntry& mp = **it;

      string  seq_id = mp[SEQ_ID];
      int32_t pos    = from_string<int32_t>(mp[POSITION]);
      int32_t strand = from_string<int32_t>(mp[STRAND]);
      bool cross_fwd = (inner3p == (strand == -1));

      // Draw over the same window the counts were taken in, with the same exact mate extents, so the
      // plot and the table agree on which reads are involved and how each counts. The width comes
      // from the same place the scorer's did -- the null carried on the candidate-region CSV, which
      // predict_missing_pairs left in summary.missing_pair -- rather than from pair_median directly,
      // so the two cannot drift apart.
      double W = (summary.missing_pair.window_width > 0.0)
               ? summary.missing_pair.window_width : max(1.0, pair_median);
      indexer.build(seq_id, static_cast<int32_t>(pos - 2 * W), static_cast<int32_t>(pos + 2 * W), ext);
      g.set_extent_index(&ext);
      g.gather(seq_id, pos, strand, cross_fwd, W);

      vector<mp_draw_read> plot_reads = g.reads();
      size_t total = mp_cap_for_display(plot_reads, settings.max_displayed_reads, static_cast<uint32_t>(pos));
      if (total)
        mp["_mp_plot_message"] = "Only " + to_string(settings.max_displayed_reads) + " of " +
                                 to_string(total) + " reads displayed.";

      string svg = settings.evidence_path + "/MP_" + mp._id + ".svg";
      render_mp_plot(svg, plot_reads, pos, strand, mate_prefix);
      make_svg_responsive(svg);
      mp["_mp_plot_file_name"] = Settings::relative_path(svg, settings.evidence_path);
    }
  }

  // ---------------------------------------------------------------------------------------------
  // One refined, scored candidate, before the local non-maximum suppression pass.
  // ---------------------------------------------------------------------------------------------
  struct mp_call {
    string   seq_id;
    int32_t  position;
    int32_t  strand;
    size_t   supporting;
    size_t   distinct;
    int      spanning;        // pairs whose mate maps past the breakpoint (evidence against)
    int      window_total;    // every crossing-strand read on the kept flank -- the TEST denominator
    double   score;           // -log10 of the genome-wide e-value
    uint32_t candidate_count;
    bool     redundant;
    bool     contig_end;
    bool     suppressed;
  };

  // Upper tail of the MP null: the chance that k or more of n reads lose their mates when the
  // genome-wide rate is p0 and the background is over-dispersed by rho. Identical in form to
  // add_sc_evidence (identify_mutations.cpp), including the fall back to a plain binomial at rho == 0.
  static double mp_tail_p_value(double k, double n, double p0, double rho)
  {
    if (!(n > 0.0) || !(k > 0.0)) return 1.0;
    if ((rho > 0.0) && (rho < 1.0)) {
      double alpha = p0 * (1.0 - rho) / rho;
      double beta  = (1.0 - p0) * (1.0 - rho) / rho;
      if ((alpha > 0.0) && (beta > 0.0)) return beta_binomial_sf(k, n, alpha, beta);
    }
    return bdtrc(k - 1.0, n, p0);
  }

  // ---------------------------------------------------------------------------------------------

  void predict_missing_pairs(const Settings& settings, Summary& summary, cReferenceSequences& ref_seq_info)
  {
    (void)ref_seq_info;

    cGenomeDiff mp_gd;

    //
    // Step 0: library geometry. Which read end faces the insertion (inner3p), how far out to rescan
    // (D), and the median pair distance whose half sets the placement shift.
    //
    bool inner3p = true;
    double D = 0.0;
    double pair_median = 0.0;
    if (!paired_library_params(summary, inner3p, D, pair_median)) {
      WARN("Missing pair (MP) evidence prediction currently supports only FR- and RF-concordant "
           "libraries. No MP evidence will be predicted.");
      mp_gd.write(settings.mp_genome_diff_file_name);
      return;
    }

    //
    // Step 1: re-read the candidate regions CSV.
    //  columns: seq_id,start,end,strand,length,max_unpaired_count,redundant,unpaired_reads
    //  unpaired_reads is the final field: ';'-joined <read_start>__<read_end>, never contains a comma.
    //
    if (!file_exists(settings.mp_candidate_regions_file_name.c_str())) {
      mp_gd.write(settings.mp_genome_diff_file_name);
      return;
    }

    vector<mp_region_row> regions;
    bool have_format = false;
    {
      ifstream in(settings.mp_candidate_regions_file_name.c_str());
      string line;
      bool past_header = false;
      while (getline(in, line)) {
        if (line.empty()) continue;
        // Calibration measured during the stage-08 pileup, written as `#key=value` lines ahead of
        // the CSV header. See identify_mutations_pileup::write_mp_candidate_regions -- the null is
        // measured over the columns where nothing was seeded, so it cannot be re-derived from the
        // regions below.
        if (line[0] == '#') {
          size_t eq = line.find('=');
          if (eq == string::npos) continue;
          string key = line.substr(1, eq - 1);
          string value = line.substr(eq + 1);
          if      (key == "mp_format")         have_format = true;
          else if (key == "window_width")      summary.missing_pair.window_width = from_string<double>(value);
          else if (key == "n_effective_tests") summary.missing_pair.n_effective_tests = from_string<double>(value);
          else if (key == "null_rate")         summary.missing_pair.null_rate = from_string<double>(value);
          else if (key == "null_rate_raw")     summary.missing_pair.null_rate_raw = from_string<double>(value);
          else if (key == "dispersion")        summary.missing_pair.dispersion = from_string<double>(value);
          else if (key == "dispersion_raw")    summary.missing_pair.dispersion_raw = from_string<double>(value);
          else if (key == "pearson_phi")       summary.missing_pair.pearson_phi = from_string<double>(value);
          else if (key == "tested_columns")    summary.missing_pair.tested_columns = from_string<uint64_t>(value);
          else if (key == "trimmed_columns")   summary.missing_pair.trimmed_columns = from_string<uint64_t>(value);
          else if (key == "mean_tested_reads") summary.missing_pair.mean_tested_reads = from_string<double>(value);
          else if (key == "seed_fraction")     summary.missing_pair.seed_fraction_used = from_string<double>(value);
          continue;
        }
        if (!past_header) { past_header = true; continue; }   // the CSV column header
        vector<string> f = split(line, ",");
        if (f.size() < 8) continue;

        mp_region_row r;
        r.seq_id = f[0];
        r.start = from_string<uint32_t>(f[1]);
        r.end = from_string<uint32_t>(f[2]);
        r.strand = f[3].empty() ? 'F' : f[3][0];
        r.max_count = from_string<uint32_t>(f[5]);
        r.max_total = from_string<uint32_t>(f[6]);
        r.redundant = (f[7] == "1");
        r.extreme_read_end = 0;
        r.extreme_read_start = 0;
        r.have_extents = false;

        string extent_field = (f.size() >= 9) ? f[8] : "";
        vector<string> extents = split(extent_field, ";");
        for (size_t i = 0; i < extents.size(); i++) {
          if (extents[i].empty()) continue;
          vector<string> ef = split(extents[i], "__");
          if (ef.size() < 2) continue;
          int32_t rs = from_string<int32_t>(ef[0]);
          int32_t re = from_string<int32_t>(ef[1]);
          if (!r.have_extents) {
            r.extreme_read_start = rs; r.extreme_read_end = re; r.have_extents = true;
          } else {
            r.extreme_read_start = min(r.extreme_read_start, rs);
            r.extreme_read_end   = max(r.extreme_read_end,   re);
          }
        }
        regions.push_back(r);
      }
    }

    ASSERT(have_format,
           "Missing pair candidate region file is missing its calibration header and is probably "
           "from an older run:\n  " + settings.mp_candidate_regions_file_name +
           "\nDelete 08_mutation_identification/mutation_identification.done and "
           "08_mutation_identification/missing_pair.done and re-run to regenerate it.");
    // A zero null rate is reachable without a stale file: --missing-pair-minimum-rate 0 removes the
    // floor, and a sample in which no mate ever failed to align then measures exactly zero. Refuse
    // rather than proceed, because every p-value against a zero rate is zero and every score is
    // infinite -- the failure mode would be to accept everything, silently.
    if (!(summary.missing_pair.null_rate > 0.0)) {
      WARN("The genome-wide rate at which a read's mate fails to align measured as zero, so there is "
           "no null for missing pair (MP) evidence to be tested against and no MP evidence will be "
           "predicted. Raise --missing-pair-minimum-rate above 0 to set a floor.");
      mp_gd.write(settings.mp_genome_diff_file_name);
      return;
    }

    summary.missing_pair.pair_distance_median = pair_median;
    summary.missing_pair.score_cutoff = settings.missing_pair_log10_e_value_cutoff;
    summary.missing_pair.regions_seeded = regions.size();

    if (regions.empty()) {
      mp_gd.write(settings.mp_genome_diff_file_name);
      return;
    }

    if (!file_exists(settings.reference_bam_file_name.c_str())
        || !file_exists(settings.reference_fasta_file_name.c_str())) {
      mp_gd.write(settings.mp_genome_diff_file_name);
      return;
    }
    mp_side_scanner* scanner =
      new mp_side_scanner(settings.reference_bam_file_name, settings.reference_fasta_file_name);
    mp_extent_indexer indexer(settings.reference_bam_file_name, settings.reference_fasta_file_name);

    //
    // Step 2: place and count each candidate.
    //
    int32_t H = static_cast<int32_t>(pair_median / 2.0 + 0.5);

    // The counting window, and the null the counts are judged against. Both come from the header the
    // pileup wrote: calibrating at one width and testing at another misstates the variance by the
    // ratio of the two. A header from a run whose read groups produced no usable width falls back to
    // the library median, which is what that width is.
    double W_null = (summary.missing_pair.window_width > 0.0)
                  ? summary.missing_pair.window_width : max(1.0, pair_median);
    double p0  = summary.missing_pair.null_rate;
    double rho = summary.missing_pair.dispersion;
    double log10_n_tests = log10(max(1.0, summary.missing_pair.n_effective_tests));

    vector<mp_call> calls;

    // Declared OUTSIDE the loop, and deliberately so: the scanner holds a POINTER to it for the whole
    // of each pass. If this lived inside the loop it would be destroyed at the end of an iteration
    // while the scanner still pointed at it, and the next iteration's placement pass would read freed
    // memory before rebuilding it -- which is exactly the use-after-free this once was. It surfaced as
    // an intermittent segfault inside map::find on macOS and as a hang on Linux (walking a corrupted
    // red-black tree can loop forever), and it survived local testing because freed memory often still
    // looks valid. Keep the lifetime of this object strictly longer than the scanner's use of it.
    mp_extent_index ext;

    for (size_t ri = 0; ri < regions.size(); ri++) {
      const mp_region_row& r = regions[ri];

      int32_t p = 0, s = 0;
      paired_region_to_side(r.strand, r.start, r.end, inner3p, p, s);

      // Prefer the reads' own extents to the region span: a region closes when a read ages out of the
      // sliding window, so its far bound is a read's start plus the window width -- a coordinate no
      // read occupies. The reads' insertion-facing aligned ends cannot pass the breakpoint, so the
      // extreme one is a conservative starting point. (Same reasoning as dp_seed_side_position.)
      if (r.have_extents)
        p = (s == -1) ? r.extreme_read_end : r.extreme_read_start;

      int32_t tid = scanner->tid_for_seq_id(r.seq_id);
      if (tid < 0) continue;
      int32_t seqlen = scanner->seq_length(tid);

      // The crossing strand for this side, matching DP's convention.
      bool cross_fwd = (inner3p == (s == -1));

      // Half-MAD placement: the median of the supporting reads' outside ends, shifted toward the
      // insertion by half the median pair distance, then clamped OUT to the innermost aligned read
      // edge so no supporting read straddles the placed coordinate. (DP computes the same quantity
      // but currently reports the aligned edge alone; MP reports the shifted median, which is the
      // better estimate here because there is no second side to cross-check it against.)
      // The placement pass classifies reads too, so it needs the mate extents just as the counting
      // pass does. Build them for the placement window BEFORE calling it -- leaving a stale index in
      // place here is what caused the use-after-free noted above.
      indexer.build(r.seq_id, static_cast<int32_t>(p - 2 * D), static_cast<int32_t>(p + 2 * D), ext);
      scanner->set_extent_index(&ext);

      int32_t median_outer = 0, inner_edge = 0;
      if (scanner->supporting_outer_median(r.seq_id, p, s, cross_fwd, D, median_outer, inner_edge)) {
        int32_t mh = median_outer + (s == -1 ? +H : -H);
        int32_t placed = (s == -1) ? max(mh, inner_edge) : min(mh, inner_edge);
        p = max(1, min(seqlen, placed));
      } else {
        // No supporting read in the rescan window -- the seed came from reads the window no longer
        // reaches. Nothing to place or count; drop the candidate rather than emit an unsupported item.
        summary.missing_pair.items_dropped_unplaced++;
        continue;
      }

      // Count over ONE fragment length, not the wide placement window D. A read further than that
      // from the breakpoint has no mate that could have reached past it, so it carries no information
      // about whether the insertion is there -- including it would only dilute the frequency by the
      // arbitrary ratio D/pair_median (here about 3x) and make the gate depend on the rescan width.
      //
      // The width comes from the CSV header, not from pair_median directly, because the null the
      // score tests against was tabulated at exactly that width during the pileup and rho's meaning
      // depends on n. They are the same number for a single read group; taking it from the header is
      // what keeps them the same number when they would otherwise drift apart.
      double W = max(1.0, W_null);
      // Rebuild around the PLACED position, which the pass above may have moved: index exact mate
      // extents a couple of fragment lengths either side, so every counted read's mate is present (a
      // proper pair's mate is within one fragment length by definition).
      indexer.build(r.seq_id, static_cast<int32_t>(p - 2 * W), static_cast<int32_t>(p + 2 * W), ext);
      scanner->set_extent_index(&ext);
      scanner->scan(r.seq_id, p, s, cross_fwd, W);

      mp_call c;
      c.seq_id = r.seq_id;
      c.position = p;
      c.strand = s;
      c.supporting = static_cast<size_t>(scanner->supporting());
      c.distinct = scanner->distinct();
      c.spanning = scanner->spanning();
      c.window_total = scanner->window_total();
      c.candidate_count = r.max_count;
      c.redundant = r.redundant;

      // The genome-wide test. Everything else about a candidate is a local sanity check; this is the
      // only thing that asks how often a pile this good arises by chance anywhere in the reference.
      //
      // Note which denominator it uses: window_total, every crossing-strand read on the kept flank,
      // NOT supporting + spanning. That is the quantity the pileup tabulated the null over, and it
      // is the only one whose null is known. supporting + spanning answers a different question --
      // of the molecules that COULD have contradicted this call, how many did not -- which is a
      // frequency, and it stays below as exactly that.
      {
        double kk = static_cast<double>(c.supporting);
        double nn = static_cast<double>(c.window_total);
        double pv = mp_tail_p_value(kk, nn, p0, rho);
        if (!(pv > 0.0) || !std::isfinite(pv)) pv = 1e-300;
        c.score = -(log10(pv) + log10_n_tests);
      }
      summary.missing_pair.items_tested++;
      // Within about one fragment length of the end of a LINEAR sequence, the mate of a flank-facing
      // read legitimately runs off the end and is unmapped, so every such end produces a maximal pile
      // that has nothing to do with an insertion. Marked IGNORE (not rejected): it is an artifact of
      // the reference's shape, not a call that came close to passing.
      c.contig_end = !ref_seq_info.is_circular(r.seq_id) &&
                     ((p <= static_cast<int32_t>(pair_median)) || (p >= seqlen - static_cast<int32_t>(pair_median)));
      c.suppressed = false;

      // Expected to turn up by chance somewhere in the reference: hopeless, and keeping these would
      // swamp the marginal evidence table with the noisiest windows in the sample. Mirrors the
      // early-out in add_sc_evidence and in pd_evidence.
      if ((settings.missing_pair_log10_e_value_cutoff > 0.0) && (c.score < 0.0)) {
        summary.missing_pair.items_dropped_score++;
        continue;
      }

      calls.push_back(c);
    }

    delete scanner;

    //
    // Step 3: local non-maximum suppression. The sliding window smears one insertion point across
    // several columns, which can close and reopen a region as reads age in and out, so one event can
    // yield several candidates that place within a fragment length of each other. Keep the
    // best-SCORING one per (seq_id, strand) neighbourhood. (Same motivation as the SC suppression
    // in identify_mutations.)
    //
    // Ranked by score rather than by raw supporting count: at unequal depth the raw count prefers
    // whichever neighbour happens to sit in more coverage, which is not the same as the one better
    // supported relative to its own background.
    //
    {
      int32_t win = max(1, static_cast<int32_t>(pair_median));
      for (size_t i = 0; i < calls.size(); i++) {
        for (size_t j = 0; j < calls.size(); j++) {
          if (i == j) continue;
          if (calls[i].seq_id != calls[j].seq_id) continue;
          if (calls[i].strand != calls[j].strand) continue;
          if (abs(calls[i].position - calls[j].position) > win) continue;
          // Strictly better, or equal but earlier -- so exactly one of any tied group survives.
          // Deliberately NOT conditioned on calls[j] itself having survived: that would make the
          // outcome depend on the order candidates happen to be visited in. "Is there a better
          // neighbour" is order-independent, at the cost of the usual non-maximum-suppression chain
          // behaviour (a weak call next to a middling call next to a strong one is dropped).
          bool j_better = (calls[j].score > calls[i].score) ||
                          ((calls[j].score == calls[i].score) && (j < i));
          if (j_better) { calls[i].suppressed = true; break; }
        }
      }
    }

    //
    // Step 3b: collapse calls that would emit the SAME MP entry.
    //
    // An MP's identity in a genome diff is (seq_id, position, strand) and nothing else -- every count
    // is a key=value field, which cDiffEntry::compare does not look at when testing for duplicates.
    // Two candidates that place on the same base and strand therefore emit identical entries, and
    // cGenomeDiff::write treats that as a fatal duplicate.
    //
    // The suppression step above does not prevent this: it only marks the weaker call, it does not
    // drop it, and marked calls are still emitted. Two candidates at the SAME position are also a
    // case suppression cannot break, since neither is nearer than the other.
    //
    // Keep exactly one call per emitted identity: the one suppression already chose, falling back to
    // the better score and then to the earlier call. Every tie-break is order-independent, so which
    // call survives does not depend on the order candidates were visited in.
    //
    vector<bool> emit(calls.size(), true);
    {
      map<string,size_t> best_by_identity;
      for (size_t i = 0; i < calls.size(); i++) {
        string identity = calls[i].seq_id + ":" + to_string(calls[i].position)
                        + ":" + to_string(calls[i].strand);
        map<string,size_t>::iterator it = best_by_identity.find(identity);
        if (it == best_by_identity.end()) {
          best_by_identity[identity] = i;
          continue;
        }
        size_t j = it->second;
        bool i_better = (calls[j].suppressed && !calls[i].suppressed) ||
                        ((calls[i].suppressed == calls[j].suppressed) &&
                         (calls[i].score > calls[j].score));
        if (i_better) {
          emit[j] = false;
          it->second = i;
        } else {
          emit[i] = false;
        }
      }
    }

    //
    // Step 4: emit.
    //
    for (size_t i = 0; i < calls.size(); i++) {
      if (!emit[i]) continue;
      const mp_call& c = calls[i];

      cDiffEntry mp(MP);
      mp[SEQ_ID] = c.seq_id;
      mp[POSITION] = to_string(c.position);
      mp[STRAND] = to_string(c.strand);

      if (c.contig_end) mp[IGNORE] = "CONTIG_END";

      double k = static_cast<double>(c.supporting);
      double n = k + static_cast<double>(c.spanning);

      mp[MP_READ_COUNT] = to_string(c.supporting);
      mp[MP_DISTINCT_COUNT] = to_string(c.distinct);
      mp[MP_CONCORDANT_COUNT] = to_string(c.spanning);
      mp[MP_TOTAL_COUNT] = to_string(static_cast<uint32_t>(n));
      mp[MP_WINDOW_COUNT] = to_string(c.window_total);
      mp[MP_CANDIDATE_COUNT] = to_string(c.candidate_count);
      mp[MP_SCORE] = to_string(formatted_double(c.score, 1));

      // The local variant frequency: of the molecules that could have CONTRADICTED this call --
      // pairs whose mate maps past the placed position, carrying reference sequence where the
      // insertion is claimed to be -- what fraction instead lost their mate.
      //
      // This says how much of the sample carries the insertion. It does NOT say whether the
      // insertion is there: that is what MP_SCORE tests, against the measured genome-wide rate at
      // which mates fail to align. The two used to be conflated, with a fixed frequency cutoff
      // standing in for a null of "the background is zero" -- which holds in a clean simulation and
      // in nothing else. On a real 2x150 bacterial library the background is about 2%, and it is
      // spread unevenly enough that 4% of all windows in the genome sit above 10%.
      double freq_lower = 0.0;
      if (n > 0.0) {
        freq_lower = binomial_frequency_lower_bound(k, n, kMPFrequencyAlpha);
        mp[FREQUENCY] = to_string(formatted_double(k / n, 4));
        mp[FREQUENCY_LOWER] = to_string(formatted_double(freq_lower, 4));
        mp[FREQUENCY_UPPER] = to_string(formatted_double(binomial_frequency_upper_bound(k, n, kMPFrequencyAlpha), 4));
      } else {
        mp[FREQUENCY] = "NA";
      }

      mp["annotate_key"] = c.redundant ? "repeat" : "gene";
      if (c.redundant) mp["redundant"] = "1";

      // Gates. The score is the one that decides; the rest are local sanity checks.
      bool rejected_by_score = (settings.missing_pair_log10_e_value_cutoff > 0.0)
                            && (c.score < settings.missing_pair_log10_e_value_cutoff);
      if (rejected_by_score)
        mp.add_reject_reason("MISSING_PAIR_SCORE");
      if (c.suppressed)
        mp.add_reject_reason("NEARBY_BETTER_MISSING_PAIR");
      if (static_cast<int32_t>(c.supporting) < settings.missing_pair_minimum_reads)
        mp.add_reject_reason("MISSING_PAIR_COUNT");
      if (static_cast<int32_t>(c.distinct) < settings.missing_pair_minimum_distinct)
        mp.add_reject_reason("MISSING_PAIR_DUPLICATES");
      if ((settings.missing_pair_frequency_cutoff > 0.0) && (freq_lower < settings.missing_pair_frequency_cutoff))
        mp.add_reject_reason("MISSING_PAIR_FREQUENCY");

      if (!mp.entry_exists(REJECT))          summary.missing_pair.items_accepted++;
      else if (rejected_by_score)            summary.missing_pair.items_rejected_score++;
      else                                   summary.missing_pair.items_rejected_other++;

      mp_gd.add(mp);
    }

    mp_gd.write(settings.mp_genome_diff_file_name);
  }

} // namespace breseq
