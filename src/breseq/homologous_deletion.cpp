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

#include "homologous_deletion.h"
#include "genome_diff_entry.h"

#include <algorithm>

using namespace std;

namespace breseq {

  void homology_base_counter::count_region(const string& seq_id, int32_t start_1, int32_t end_1)
  {
    if (end_1 < start_1) return;
    m_region_start = start_1;
    m_region_end = end_1;
    // clip=true so the pileup only hands us the columns we asked for. Without it, do_pileup() calls
    // back for every position from 1 up to the region as well (handle_position only gates on the clip
    // bounds), which for a region a few megabases in would be millions of empty columns.
    do_pileup(seq_id + ":" + to_string(start_1) + "-" + to_string(end_1), true);
  }

  void homology_base_counter::pileup_callback(const pileup& p)
  {
    int32_t position = static_cast<int32_t>(p.position_1());
    // Belt and braces with the clip above: never tally a column outside the requested region.
    if ((position < m_region_start) || (position > m_region_end)) return;

    bool record = !m_record_positions.empty() && (m_record_positions.count(position) > 0);
    base_tally& t = m_tallies[position];

    for (pileup::const_iterator i = p.begin(); i != p.end(); i++) {

      // The read deletes this reference base, so it has no base to contribute.
      if (i->on_base_indel() < 0) continue;

      uint8_t read_base_bam = i->read_base_bam_0(i->query_position_0());
      if (_base_bam_is_N(read_base_bam)) continue;
      if (i->read_base_quality_0(i->query_position_0()) < m_base_quality_cutoff) continue;

      char base_char = basebam2char(read_base_bam);
      uint8_t base_index = basechar2index(base_char);
      if (base_index > 3) continue;

      // Note what is deliberately NOT filtered here -- redundancy, trimming, mapping quality. See
      // the class comment: the reads this predictor needs are exactly the ones identify_mutations
      // discards, and every placement counts once because the two homologous columns are pooled.
      t.n[base_index]++;
      t.total++;

      if (record) {
        // Key per PLACEMENT, not per read: the same read placed at both copies, and the two mates of
        // a pair, must not be merged into one allele series.
        string key = i->read_name() + ":" + to_string(i->reference_start_1());
        m_read_alleles[key].push_back(make_pair(position - m_position_offset, base_char));
      }
    }
  }

  double homology_identity(const cAnnotatedSequence& seq, int32_t d, int32_t lo, int32_t hi)
  {
    if (d <= 0) return 0.0;
    int32_t seq_len = static_cast<int32_t>(seq.m_length);

    // Clamp so that both x and x+d stay inside the sequence.
    int32_t x_lo = max(1, lo);
    int32_t x_hi = min(hi, seq_len - d);
    if (x_hi < x_lo) return 0.0;

    int32_t same = 0;
    for (int32_t x = x_lo; x <= x_hi; x++)
      if (seq.get_sequence_1(x) == seq.get_sequence_1(x + d)) same++;

    return static_cast<double>(same) / static_cast<double>(x_hi - x_lo + 1);
  }

  int32_t find_homologous_register(const cAnnotatedSequence& seq,
                                   int32_t anchor_left, int32_t anchor_right, int32_t d_seed,
                                   int32_t search_half_width, int32_t window,
                                   double min_identity, double min_margin, int32_t margin_exclusion,
                                   double& best_identity, double& runner_up_identity)
  {
    best_identity = 0.0;
    runner_up_identity = 0.0;

    int32_t seq_len = static_cast<int32_t>(seq.m_length);
    int32_t d_min = max(1, d_seed - search_half_width);
    int32_t d_max = min(seq_len - 1, d_seed + search_half_width);
    if (d_max < d_min) return 0;

    // The two anchor windows are fixed; only their homologous partners slide with d. Pull each
    // anchor and the whole span its partner can occupy out of the sequence once, so scanning ten
    // thousand registers is a few million array lookups rather than a few million bounds-checked
    // calls.
    int32_t a_lo = max(1, anchor_left - window);
    int32_t a_hi = min(seq_len, anchor_left + window);
    int32_t b_lo = max(1, anchor_right - window);
    int32_t b_hi = min(seq_len, anchor_right + window);

    string a_win, a_partner, b_win, b_partner;
    int32_t a_partner_lo = 0, b_partner_lo = 0;

    if (a_hi >= a_lo) {
      a_win = seq.get_sequence_1(a_lo, a_hi);
      a_partner_lo = a_lo + d_min;
      int32_t a_partner_hi = min(seq_len, a_hi + d_max);
      if (a_partner_hi >= a_partner_lo) a_partner = seq.get_sequence_1(a_partner_lo, a_partner_hi);
    }
    if (b_hi >= b_lo) {
      b_win = seq.get_sequence_1(b_lo, b_hi);
      b_partner_lo = max(1, b_lo - d_max);
      int32_t b_partner_hi = b_hi - d_min;
      if (b_partner_hi >= b_partner_lo) b_partner = seq.get_sequence_1(b_partner_lo, b_partner_hi);
    }

    int32_t best_d = 0;
    vector<double> identity_by_d(d_max - d_min + 1, 0.0);

    for (int32_t d = d_min; d <= d_max; d++) {

      // Anchor on the left MC boundary: the left copy at x against the right copy at x+d.
      int32_t a_same = 0, a_total = 0;
      for (int32_t x = a_lo; x <= a_hi; x++) {
        int32_t partner = x + d;
        if (partner > seq_len) break;
        size_t pi = static_cast<size_t>(partner - a_partner_lo);
        if (pi >= a_partner.size()) break;
        a_total++;
        if (a_win[x - a_lo] == a_partner[pi]) a_same++;
      }

      // Anchor on the right MC boundary: the right copy at y against the left copy at y-d. Either
      // boundary may have landed outside the homologous block, so take the better of the two rather
      // than letting the worse one veto a real register.
      int32_t b_same = 0, b_total = 0;
      for (int32_t y = b_lo; y <= b_hi; y++) {
        int32_t partner = y - d;
        if (partner < 1) continue;
        size_t pi = static_cast<size_t>(partner - b_partner_lo);
        if (pi >= b_partner.size()) continue;
        b_total++;
        if (b_win[y - b_lo] == b_partner[pi]) b_same++;
      }

      // Near a contig end an anchor window gets truncated, and a handful of positions can agree
      // perfectly by chance. Require most of the window to be measurable before believing the number.
      int32_t min_measured = window;
      double a_identity = (a_total >= min_measured) ? static_cast<double>(a_same) / a_total : 0.0;
      double b_identity = (b_total >= min_measured) ? static_cast<double>(b_same) / b_total : 0.0;
      double identity = max(a_identity, b_identity);

      identity_by_d[d - d_min] = identity;
      if (identity > best_identity) { best_identity = identity; best_d = d; }
    }

    if (best_d == 0) return 0;
    if (best_identity < min_identity) return 0;

    // Between paralogs that are identical in long stretches everywhere, a high score on its own means
    // nothing -- what distinguishes the true register is that nothing else comes close. Measure that
    // against the best register far enough away to be a genuinely different alignment.
    for (int32_t d = d_min; d <= d_max; d++) {
      if (abs(d - best_d) < margin_exclusion) continue;
      if (identity_by_d[d - d_min] > runner_up_identity) runner_up_identity = identity_by_d[d - d_min];
    }
    if (best_identity - runner_up_identity < min_margin) return 0;

    return best_d;
  }

  void homologous_block_extent(const cAnnotatedSequence& seq, int32_t d, int32_t seed,
                               int32_t step, double min_identity, int32_t max_extent,
                               int32_t& block_start, int32_t& block_end)
  {
    int32_t seq_len = static_cast<int32_t>(seq.m_length);
    block_start = seed;
    block_end = seed;
    if ((d <= 0) || (seed < 1) || (seed > seq_len - d)) return;

    while ((block_end - seed) < max_extent) {
      int32_t lo = block_end + 1;
      int32_t hi = min(block_end + step, seq_len - d);
      if (hi < lo) break;
      if (homology_identity(seq, d, lo, hi) < min_identity) break;
      block_end = hi;
    }

    while ((seed - block_start) < max_extent) {
      int32_t hi = block_start - 1;
      int32_t lo = max(1, block_start - step);
      if (hi < lo) break;
      if (homology_identity(seq, d, lo, hi) < min_identity) break;
      block_start = lo;
    }
  }

  void discriminating_positions(const cAnnotatedSequence& seq, int32_t d,
                                int32_t block_start, int32_t block_end,
                                vector<int32_t>& positions)
  {
    positions.clear();
    if (d <= 0) return;
    int32_t seq_len = static_cast<int32_t>(seq.m_length);
    int32_t lo = max(1, block_start);
    int32_t hi = min(block_end, seq_len - d);

    for (int32_t x = lo; x <= hi; x++)
      if (seq.get_sequence_1(x) != seq.get_sequence_1(x + d)) positions.push_back(x);
  }

  //! Fill in position/left_base/right_base and reject columns the reference cannot discriminate.
  static bool init_homology_column(const cAnnotatedSequence& seq, int32_t d, int32_t x,
                                   homology_column& c)
  {
    c = homology_column();
    c.position = x;
    c.left_base = seq.get_sequence_1(x);
    c.right_base = seq.get_sequence_1(x + d);
    if (c.left_base == c.right_base) return false;

    // An N in either copy discriminates nothing.
    if ((c.left_base == 'N') || (c.right_base == 'N')) return false;
    return true;
  }

  void build_homology_columns(const cAnnotatedSequence& seq, int32_t d,
                              const vector<int32_t>& positions,
                              const base_tally_map_t& tallies,
                              uint32_t min_column_coverage, double min_allele_fraction,
                              homology_columns_t& columns)
  {
    columns.clear();

    for (size_t i = 0; i < positions.size(); i++) {
      homology_column c;
      if (!init_homology_column(seq, d, positions[i], c)) continue;

      uint8_t left_index = basechar2index(c.left_base);
      uint8_t right_index = basechar2index(c.right_base);
      if ((left_index > 3) || (right_index > 3)) continue;

      base_tally_map_t::const_iterator t1 = tallies.find(c.position);
      base_tally_map_t::const_iterator t2 = tallies.find(c.position + d);

      if (t1 != tallies.end()) {
        c.left_allele += t1->second.n[left_index];
        c.right_allele += t1->second.n[right_index];
        c.total += t1->second.total;
      }
      if (t2 != tallies.end()) {
        c.left_allele += t2->second.n[left_index];
        c.right_allele += t2->second.n[right_index];
        c.total += t2->second.total;
      }

      // With both copies still present this lands near 0.5 on every column and nothing is called,
      // which is exactly the behaviour we want for a locus that has not been deleted at all.
      if (c.total > 0) {
        if ((c.right_allele >= min_column_coverage)
            && (c.right_allele >= min_allele_fraction * c.total)) c.call = +1;
        else if ((c.left_allele >= min_column_coverage)
                 && (c.left_allele >= min_allele_fraction * c.total)) c.call = -1;
      }

      columns.push_back(c);
    }
  }

  void build_homology_columns_from_ra(const cAnnotatedSequence& seq, int32_t d,
                                      const vector<int32_t>& positions,
                                      const map<int32_t, string>& ra_base,
                                      homology_columns_t& columns)
  {
    columns.clear();

    for (size_t i = 0; i < positions.size(); i++) {
      homology_column c;
      if (!init_homology_column(seq, d, positions[i], c)) continue;

      // Only a called difference says anything here. An RA on the left copy carrying the right
      // copy's base means the hybrid holds the right allele at this column; an RA on the right copy
      // carrying the left copy's base means it holds the left allele. Absence of an RA is not
      // evidence, because we cannot tell a matching base from an uncovered one.
      map<int32_t, string>::const_iterator left_ra = ra_base.find(c.position);
      map<int32_t, string>::const_iterator right_ra = ra_base.find(c.position + d);

      if ((left_ra != ra_base.end()) && (left_ra->second.size() == 1)
          && (left_ra->second[0] == c.right_base)) {
        c.call = +1;
      } else if ((right_ra != ra_base.end()) && (right_ra->second.size() == 1)
                 && (right_ra->second[0] == c.left_base)) {
        c.call = -1;
      }

      columns.push_back(c);
    }
  }

  bool find_crossover_interval(const homology_columns_t& columns,
                               int32_t min_columns_each_side, int32_t max_violations,
                               int32_t& last_left, int32_t& first_right,
                               int32_t& n_left, int32_t& n_right, int32_t& n_violations)
  {
    last_left = 0;
    first_right = 0;
    n_left = 0;
    n_right = 0;
    n_violations = 0;

    vector<const homology_column*> called;
    for (size_t i = 0; i < columns.size(); i++)
      if (columns[i].call != 0) called.push_back(&columns[i]);
    if (called.empty()) return false;

    int32_t n = static_cast<int32_t>(called.size());
    int32_t total_left = 0;
    for (int32_t i = 0; i < n; i++)
      if (called[i]->call == -1) total_left++;

    // A real crossover is a step: left allele, then right allele, once. Score every split by how many
    // called columns it puts on the wrong side and take the best one, so a single mismapped column
    // shifts the count instead of destroying the answer.
    int32_t best_split = 0, best_cost = n + 1;
    int32_t plus_before = 0, minus_before = 0;
    for (int32_t k = 0; k <= n; k++) {
      if (k > 0) {
        if (called[k-1]->call == +1) plus_before++;
        else minus_before++;
      }
      int32_t cost = plus_before + (total_left - minus_before);
      if (cost < best_cost) { best_cost = cost; best_split = k; }
    }
    n_violations = best_cost;

    for (int32_t i = 0; i < best_split; i++) {
      if (called[i]->call != -1) continue;
      n_left++;
      if (called[i]->position > last_left) last_left = called[i]->position;
    }
    for (int32_t i = best_split; i < n; i++) {
      if (called[i]->call != +1) continue;
      n_right++;
      if ((first_right == 0) || (called[i]->position < first_right)) first_right = called[i]->position;
    }

    if (n_left < min_columns_each_side) return false;
    if (n_right < min_columns_each_side) return false;
    if (n_violations > max_violations) return false;
    if ((last_left == 0) || (first_right == 0) || (last_left >= first_right)) return false;

    return true;
  }

  double read_allele_consistency(const cAnnotatedSequence& seq, int32_t d,
                                 const read_alleles_t& read_alleles,
                                 int32_t& informative_reads)
  {
    informative_reads = 0;
    int32_t consistent_reads = 0;
    int32_t seq_len = static_cast<int32_t>(seq.m_length);

    for (read_alleles_t::const_iterator it = read_alleles.begin(); it != read_alleles.end(); it++) {

      read_allele_list_t observed = it->second;
      sort(observed.begin(), observed.end());

      vector<int32_t> calls;
      for (size_t i = 0; i < observed.size(); i++) {
        int32_t x = observed[i].first;
        if ((x < 1) || (x > seq_len - d)) continue;
        char left_base = seq.get_sequence_1(x);
        char right_base = seq.get_sequence_1(x + d);
        // A base matching neither copy is a sequencing error and says nothing about the register.
        if (observed[i].second == left_base) calls.push_back(-1);
        else if (observed[i].second == right_base) calls.push_back(+1);
      }
      if (calls.size() < 2) continue;

      informative_reads++;

      // The read itself must cross the boundary at most once, in the one direction the hybrid allows.
      bool seen_right = false, consistent = true;
      for (size_t i = 0; i < calls.size(); i++) {
        if (calls[i] == +1) seen_right = true;
        else if (seen_right) { consistent = false; break; }
      }
      if (consistent) consistent_reads++;
    }

    if (informative_reads == 0) return 1.0;
    return static_cast<double>(consistent_reads) / static_cast<double>(informative_reads);
  }

  cDiffEntry make_homologous_deletion(cAnnotatedSequence& seq, int32_t position, int32_t size,
                                      const vector<string>& evidence_ids,
                                      bool polymorphism_prediction)
  {
    cDiffEntry mut(DEL);
    mut[SEQ_ID] = seq.m_seq_id;
    mut[POSITION] = to_string(position);
    mut[SIZE] = to_string(size);
    mut._evidence = evidence_ids;
    if (polymorphism_prediction) mut[FREQUENCY] = "1";

    // Name the two homologous features the deletion sits between, following the 'between' convention
    // used for a deletion between two copies of a repeat family. annotate_repeat_hotspots() preserves
    // a value set here unless it finds annotated repeat_regions at the boundaries.
    cFeatureLocation* f1 = cReferenceSequences::get_overlapping_feature(seq.m_gene_locations, position);
    cFeatureLocation* f2 = cReferenceSequences::get_overlapping_feature(seq.m_gene_locations, position + size - 1);
    if ((f1 != NULL) && (f2 != NULL)) {
      string n1 = (*f1->get_feature())["name"], n2 = (*f2->get_feature())["name"];
      if (!n1.empty() && !n2.empty()) mut["between"] = (n1 == n2) ? n1 : (n1 + "," + n2);
    }

    return mut;
  }

} // namespace breseq
