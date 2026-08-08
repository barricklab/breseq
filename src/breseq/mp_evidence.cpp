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
        m_collect(false), m_ext(NULL), m_supporting(0), m_spanning(0), m_have_inner(false), m_inner_edge(0)
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
      set_ctx(p, s, cross_fwd);
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
    void scan(const string& seq_id, int32_t p, int32_t s, bool cross_fwd, double D)
    {
      set_ctx(p, s, cross_fwd);
      m_supporting = 0; m_spanning = 0;
      m_distinct.clear();
      m_collect = false;
      fetch_side_window(seq_id, p, s, D);
    }

    int supporting() const { return m_supporting; }
    int spanning() const { return m_spanning; }
    size_t distinct() const { return m_distinct.size(); }

    void fetch_callback(const alignment_wrapper& a)
    {
      int32_t outer = 0, inner = 0;
      int cat = classify(a, outer, inner);
      if (m_collect) {
        if (cat == 1) {
          m_outside.push_back(outer);
          if (!m_have_inner || (m_s == -1 ? inner > m_inner_edge : inner < m_inner_edge)) {
            m_have_inner = true; m_inner_edge = inner;
          }
        }
        return;
      }
      if      (cat == 1) { m_supporting++; m_distinct.insert(outer); }
      else if (cat == 2) m_spanning++;
    }

  private:
    void set_ctx(int32_t p, int32_t s, bool cross_fwd) { m_p = p; m_s = s; m_cross_fwd = cross_fwd; }

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

      // Supporting: paired, mapped, mate unmapped anywhere (mark_mate_unmapped).
      if (a.is_paired() && (a.flag() & BAM_FMUNMAP)) return 1;

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
    bool    m_collect;
    const mp_extent_index* m_ext;
    int     m_supporting, m_spanning;
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
      : pileup_base(bam, fasta), m_p(0), m_s(0), m_cross_fwd(true), m_ext(NULL) { set_print_progress(false); }

    //! Exact mate extents for the position about to be gathered. Not owned; must outlive the gather.
    void set_extent_index(const mp_extent_index* ext) { m_ext = ext; }

    int32_t tid_for_seq_id(const string& seq_id) const {
      for (uint32_t t = 0; t < num_targets(); t++)
        if (seq_id == string(target_name(t))) return static_cast<int32_t>(t);
      return -1;
    }

    void gather(const string& seq_id, int32_t p, int32_t s, bool cross_fwd, double W)
    {
      m_p = p; m_s = s; m_cross_fwd = cross_fwd;
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
      bool mate_unmapped = a.is_paired() && ((a.flag() & BAM_FMUNMAP) != 0);
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

    for (diff_entry_list_t::iterator it = mp_list.begin(); it != mp_list.end(); it++) {
      cDiffEntry& mp = **it;

      string  seq_id = mp[SEQ_ID];
      int32_t pos    = from_string<int32_t>(mp[POSITION]);
      int32_t strand = from_string<int32_t>(mp[STRAND]);
      bool cross_fwd = (inner3p == (strand == -1));

      // Draw over the same one-fragment-length window the counts were taken in, with the same exact
      // mate extents, so the plot and the table agree on which reads are involved and how each counts.
      double W = max(1.0, pair_median);
      mp_extent_index ext;
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
    uint32_t candidate_count;
    bool     redundant;
    bool     contig_end;
    bool     suppressed;
  };

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
    {
      ifstream in(settings.mp_candidate_regions_file_name.c_str());
      string line;
      getline(in, line); // header
      while (getline(in, line)) {
        if (line.empty()) continue;
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
    vector<mp_call> calls;

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
      int32_t median_outer = 0, inner_edge = 0;
      if (scanner->supporting_outer_median(r.seq_id, p, s, cross_fwd, D, median_outer, inner_edge)) {
        int32_t mh = median_outer + (s == -1 ? +H : -H);
        int32_t placed = (s == -1) ? max(mh, inner_edge) : min(mh, inner_edge);
        p = max(1, min(seqlen, placed));
      } else {
        // No supporting read in the rescan window -- the seed came from reads the window no longer
        // reaches. Nothing to place or count; drop the candidate rather than emit an unsupported item.
        continue;
      }

      // Count over ONE fragment length, not the wide placement window D. A read further than that
      // from the breakpoint has no mate that could have reached past it, so it carries no information
      // about whether the insertion is there -- including it would only dilute the frequency by the
      // arbitrary ratio D/pair_median (here about 3x) and make the gate depend on the rescan width.
      double W = max(1.0, pair_median);
      // Index exact mate extents a couple of fragment lengths either side, so every counted read's
      // mate is present (a proper pair's mate is within one fragment length by definition).
      mp_extent_index ext;
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
      c.candidate_count = r.max_count;
      c.redundant = r.redundant;
      // Within about one fragment length of the end of a LINEAR sequence, the mate of a flank-facing
      // read legitimately runs off the end and is unmapped, so every such end produces a maximal pile
      // that has nothing to do with an insertion. Marked IGNORE (not rejected): it is an artifact of
      // the reference's shape, not a call that came close to passing.
      c.contig_end = !ref_seq_info.is_circular(r.seq_id) &&
                     ((p <= static_cast<int32_t>(pair_median)) || (p >= seqlen - static_cast<int32_t>(pair_median)));
      c.suppressed = false;
      calls.push_back(c);
    }

    delete scanner;

    //
    // Step 3: local non-maximum suppression. The sliding window smears one insertion point across
    // several columns, which can close and reopen a region as reads age in and out, so one event can
    // yield several candidates that place within a fragment length of each other. Keep the
    // best-supported one per (seq_id, strand) neighbourhood. (Same motivation as the SC suppression
    // in identify_mutations.)
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
      mp[MP_CANDIDATE_COUNT] = to_string(c.candidate_count);

      // The local variant frequency: of the molecules sampled at this point on this strand, what
      // fraction lost their mate. This is the gate that separates a real insertion (most crossing
      // fragments lose their mate, so the fraction approaches the variant frequency) from the
      // genome-wide background of reads whose mates simply failed to align, which runs at a few
      // percent and is roughly constant along the genome.
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

      // Gates.
      if (c.suppressed)
        mp.add_reject_reason("NEARBY_BETTER_MISSING_PAIR");
      if (static_cast<int32_t>(c.supporting) < settings.missing_pair_minimum_reads)
        mp.add_reject_reason("MISSING_PAIR_COUNT");
      if (static_cast<int32_t>(c.distinct) < settings.missing_pair_minimum_distinct)
        mp.add_reject_reason("MISSING_PAIR_DUPLICATES");
      if ((settings.missing_pair_frequency_cutoff > 0.0) && (freq_lower < settings.missing_pair_frequency_cutoff))
        mp.add_reject_reason("MISSING_PAIR_FREQUENCY");

      mp_gd.add(mp);
    }

    mp_gd.write(settings.mp_genome_diff_file_name);
  }

} // namespace breseq
