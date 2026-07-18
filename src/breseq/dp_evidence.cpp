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

#include <set>

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
    bool    guard_discordant;     // if true, a discordant read that extends past p is not counted (the
                                  // same "kept-clear" guard used for concordant reads). Off during the
                                  // preliminary refinement pass (which needs reads that reach past p to
                                  // find where p should be), on for counting/plotting.
  };

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
      // When counting (guard on), a discordant read must also lie entirely on the kept side and not
      // extend past p -- the same guard as concordant. p has been shifted (predict step 4) to clear all
      // the discordant reads, so this normally excludes nothing; it just enforces the invariant.
      if (c.guard_discordant) {
        bool kept_clear = (c.s == -1) ? (rend < c.p) : (rstart > c.p);
        if (!kept_clear) return 0;
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

  // Counts the three read categories at a junction side (used to fill the DP evidence fields), and
  // provides a preliminary refinement pass over the discordant reads at a side.
  class dp_side_scanner : public pileup_base {
  public:
    dp_side_scanner(const string& bam, const string& fasta)
      : pileup_base(bam, fasta), m_collect_outside(false) { set_print_progress(false); }

    int32_t tid_for_seq_id(const string& seq_id) const { return dp_tid_for_seq_id(*this, seq_id); }
    int32_t seq_length(int32_t tid) const { return static_cast<int32_t>(target_length(tid)); }

    void scan(const string& seq_id, int32_t p, int32_t s, bool crossing_is_forward,
              int32_t other_tid, int32_t other_p, bool other_crossing_is_forward, double D)
    {
      set_ctx(p, s, crossing_is_forward, other_tid, other_p, other_crossing_is_forward, D, /*guard=*/true);
      m_supporting = 0; m_concordant = 0; m_unpaired = 0;
      m_collect_outside = false;

      int32_t lo, hi;
      if (!dp_side_window(*this, seq_id, p, s, D, lo, hi)) return;
      do_fetch(seq_id + ":" + to_string(lo) + "-" + to_string(hi));
    }

    // Preliminary refinement pass: fetch the discordant (supporting) reads at this side (guard OFF, so
    // reads reaching past p are included) and collect, per read: the "outside" coordinate -- the end
    // facing AWAY from the junction (rstart for s=-1, rend for s=+1) -- and the "junction extent" -- the
    // end facing TOWARD the junction (rend for s=-1, rstart for s=+1). Returns true if any discordant
    // reads were found; sets median_outside (median of the outside coords; robust to a stray pair) and
    // junction_extent (max rend for s=-1 / min rstart for s=+1 -- the furthest a read reaches toward the
    // junction, used to shift p so no read overlaps it).
    bool refine_outside_median(const string& seq_id, int32_t p, int32_t s, bool crossing_is_forward,
                               int32_t other_tid, int32_t other_p, bool other_crossing_is_forward,
                               double D, int32_t& median_outside, int32_t& junction_extent)
    {
      set_ctx(p, s, crossing_is_forward, other_tid, other_p, other_crossing_is_forward, D, /*guard=*/false);
      m_outside.clear();
      m_have_extent = false; m_extent = 0;
      m_collect_outside = true;

      int32_t lo, hi;
      if (dp_side_window(*this, seq_id, p, s, D, lo, hi))
        do_fetch(seq_id + ":" + to_string(lo) + "-" + to_string(hi));

      m_collect_outside = false;
      if (m_outside.empty()) return false;
      sort(m_outside.begin(), m_outside.end());
      median_outside = m_outside[m_outside.size() / 2];
      junction_extent = m_extent;
      return true;
    }

    void fetch_callback(const alignment_wrapper& a) {
      int32_t anchor;
      int cat = dp_classify_side_read(a, m_ctx, anchor);
      if (m_collect_outside) {
        if (cat == 1) {
          int32_t rstart = static_cast<int32_t>(a.reference_start_1());
          int32_t rend   = static_cast<int32_t>(a.reference_end_1());
          m_outside.push_back(m_ctx.s == -1 ? rstart : rend);
          int32_t je = (m_ctx.s == -1) ? rend : rstart;   // junction-facing end
          if (!m_have_extent)          { m_extent = je; m_have_extent = true; }
          else if (m_ctx.s == -1)      { if (je > m_extent) m_extent = je; }
          else                         { if (je < m_extent) m_extent = je; }
        }
        return;
      }
      if      (cat == 1) m_supporting++;
      else if (cat == 2) m_concordant++;
      else if (cat == 3) m_unpaired++;
    }

    int supporting() const { return m_supporting; }
    int concordant() const { return m_concordant; }
    int unpaired()   const { return m_unpaired; }

  private:
    void set_ctx(int32_t p, int32_t s, bool crossing_is_forward,
                 int32_t other_tid, int32_t other_p, bool other_crossing_is_forward, double D, bool guard) {
      m_ctx.p = p; m_ctx.s = s; m_ctx.cross_fwd = crossing_is_forward;
      m_ctx.other_tid = other_tid; m_ctx.other_p = other_p; m_ctx.other_cross_fwd = other_crossing_is_forward;
      m_ctx.D = D; m_ctx.guard_discordant = guard;
    }
    dp_side_ctx m_ctx;
    int     m_supporting, m_concordant, m_unpaired;
    bool    m_collect_outside;
    vector<int32_t> m_outside;
    bool    m_have_extent; int32_t m_extent;
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
      m_ctx.D = D; m_ctx.guard_discordant = true;  // plot the reads that are actually counted
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

  void predict_discordant_pairs(const Settings& settings, Summary& summary, cReferenceSequences& ref_seq_info)
  {
    (void)ref_seq_info;

    cGenomeDiff dp_gd;

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

      // Each side's crossing-read strand (same as the region strand that produced it).
      bool s1_fwd = (inner3p == (s1_strand == -1));
      bool s2_fwd = (inner3p == (s2_strand == -1));
      int32_t s1_tid = scanner ? scanner->tid_for_seq_id(s1_seq_id) : -1;
      int32_t s2_tid = scanner ? scanner->tid_for_seq_id(s2_seq_id) : -1;

      // Refine the two breakpoint coordinates from the discordant reads BEFORE counting. For each side
      // we take the median "outside" position of its discordant reads (robust to a stray pair) and move
      // the breakpoint from there toward the junction by half the concordant pair-distance median. Doing
      // this on both sides rebalances the pair so the discordant reads sit centered between the two
      // coordinates. Then we nudge the coordinate past the junction-facing end of every discordant read
      // (max rend for s=-1 / min rstart for s=+1) so NO read overlaps it -- this is what lets counting
      // apply the same "does not extend past p" guard to discordant reads as to concordant reads without
      // discarding the very reads that support the junction. Both refinements use the region-based
      // positions as input, then update together.
      if (scanner && pair_median > 0.0) {
        int32_t med1, ext1, med2, ext2;
        bool r1 = scanner->refine_outside_median(s1_seq_id, s1_pos, s1_strand, s1_fwd, s2_tid, s2_pos, s2_fwd, distance_cutoff, med1, ext1);
        bool r2 = scanner->refine_outside_median(s2_seq_id, s2_pos, s2_strand, s2_fwd, s1_tid, s1_pos, s1_fwd, distance_cutoff, med2, ext2);
        int32_t half = static_cast<int32_t>(pair_median / 2.0 + 0.5);
        if (r1) {
          int32_t p = med1 + (s1_strand == -1 ? +half : -half);   // move toward the junction
          if (s1_strand == -1) { if (p < ext1 + 1) p = ext1 + 1; } // clear all reads (rend < p)
          else                 { if (p > ext1 - 1) p = ext1 - 1; } // clear all reads (rstart > p)
          s1_pos = max(1, min(scanner->seq_length(s1_tid), p));
        }
        if (r2) {
          int32_t p = med2 + (s2_strand == -1 ? +half : -half);
          if (s2_strand == -1) { if (p < ext2 + 1) p = ext2 + 1; }
          else                 { if (p > ext2 - 1) p = ext2 - 1; }
          s2_pos = max(1, min(scanner->seq_length(s2_tid), p));
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

      if (scanner) {
        // Count the three read categories at each (now refined) side.
        scanner->scan(s1_seq_id, s1_pos, s1_strand, s1_fwd, s2_tid, s2_pos, s2_fwd, distance_cutoff);
        int c1a = scanner->supporting(), c2a = scanner->concordant(), c3a = scanner->unpaired();

        scanner->scan(s2_seq_id, s2_pos, s2_strand, s2_fwd, s1_tid, s1_pos, s1_fwd, distance_cutoff);
        int c1b = scanner->supporting(), c2b = scanner->concordant(), c3b = scanner->unpaired();

        // The two sides should agree on the supporting count (one mate of each pair per side).
        if (c1a != c1b) {
          WARN("DP supporting count differs between sides (" + to_string(c1a) + " vs " + to_string(c1b) +
               ") for " + s1_seq_id + ":" + to_string(s1_pos) + " <-> " + s2_seq_id + ":" + to_string(s2_pos));
        }
        dp["discordant_count"] = to_string(max(c1a, c1b));
        dp["side_1_concordant_count"] = to_string(c2a);
        dp["side_2_concordant_count"] = to_string(c2b);
        dp["side_1_unpaired_count"] = to_string(c3a);
        dp["side_2_unpaired_count"] = to_string(c3b);
      } else {
        // No BAM available: fall back to the heuristic count.
        dp["discordant_count"] = to_string(weight);
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

      // The read as an arrow tail->head in its mapping direction.
      int32_t tail = r.read_reversed ? r.read_end : r.read_start;
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

    vector<string> clauses;  // connectors first so the arrows draw on top
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
    // Genomic -> joined-x transforms (side_1 maps to <=0, side_2 to >=0).
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
      // side_1 read arrow (tail->head in mapping direction), transformed.
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

    // Shade the two flanks: side_1 (left) light green, side_2 (right) light yellow.
    s << "set object 1 rectangle from " << xmin << ",0 to 0," << (n + 1)
      << " fc rgb " << DP_SIDE1_COLOR << " fs solid noborder behind" << endl;
    s << "set object 2 rectangle from 0,0 to " << xmax << "," << (n + 1)
      << " fc rgb " << DP_SIDE2_COLOR << " fs solid noborder behind" << endl;
    // Seam marker at x=0.
    s << "set arrow from 0,0 to 0," << (n + 1) << " nohead lc rgb 'gray50' lw 1 dt 3 back" << endl;

    vector<string> clauses;
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
  static string dp_read_num(const string& name) {
    size_t colon = name.find(':');
    return (colon == string::npos) ? name : name.substr(colon + 1);
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

      // Side 1 -- render the per-side plot and keep the discordant reads for the joined plot.
      g.gather(s1_seq, s1_pos, s1_str, s1_fwd, s2_tid, s2_pos, s2_fwd, D);
      vector<dp_draw_read> s1_disc;
      for (vector<dp_draw_read>::const_iterator r = g.reads().begin(); r != g.reads().end(); r++)
        if (r->category == 1) s1_disc.push_back(*r);
      {
        string svg = settings.evidence_path + "/DP_SIDE_1_" + dp._id + ".svg";
        render_dp_side_plot(svg, g.reads(), s1_pos, s1_str, mate_prefix, DP_SIDE1_COLOR);
        make_svg_responsive(svg);
        dp["_side_1_dp_plot_file_name"] = Settings::relative_path(svg, settings.evidence_path);
      }

      // Side 2
      g.gather(s2_seq, s2_pos, s2_str, s2_fwd, s1_tid, s1_pos, s1_fwd, D);
      vector<dp_draw_read> s2_disc;
      for (vector<dp_draw_read>::const_iterator r = g.reads().begin(); r != g.reads().end(); r++)
        if (r->category == 1) s2_disc.push_back(*r);
      {
        string svg = settings.evidence_path + "/DP_SIDE_2_" + dp._id + ".svg";
        render_dp_side_plot(svg, g.reads(), s2_pos, s2_str, mate_prefix, DP_SIDE2_COLOR);
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
        string svg = settings.evidence_path + "/DP_" + dp._id + ".svg";
        render_dp_joined_plot(svg, pairs, s1_pos, s1_str, s2_pos, s2_str);
        make_svg_responsive(svg);
        dp["_dp_plot_file_name"] = Settings::relative_path(svg, settings.evidence_path);
      }
    }
  }

} // namespace breseq
