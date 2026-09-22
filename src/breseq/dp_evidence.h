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

#ifndef _BRESEQ_DP_EVIDENCE_H_
#define _BRESEQ_DP_EVIDENCE_H_

#include "common.h"
#include "settings.h"
#include "summary.h"
#include "reference_sequence.h"

namespace breseq {

  //! Paired-library geometry shared by every pair-based evidence type (DP, MP, PD).
  //
  //  Votes the majority orientation across the paired read groups into a pair_geometry (alignment.h),
  //  and also returns the rescan window half-width D (max distance_cutoff over the groups) and the
  //  median paired-mapping distance. Returns false, leaving the outputs as set so far, if the
  //  orientation names no supported geometry; the caller emits its own warning.
  bool paired_library_geometry(const Summary& summary, pair_geometry& geometry, double& D, double& pair_median);

  //! A candidate region of reads that all FACE one way, as a junction side. Reads facing right sit on
  //  the left flank of whatever they reach toward, so the side is the region's END with the retained
  //  flank at <= position (strand -1); reads facing left give the region's START and strand +1. For an
  //  FR library "faces right" is "forward", so this is the rule that was always used there; the
  //  facing is what pair_geometry computes for every supported library.
  void paired_region_facing_to_side(bool faces_right, uint32_t region_start, uint32_t region_end,
                                    int32_t& position, int32_t& strand);

  //! Small gnuplot axis-formatting helpers, shared by the pair-based evidence plots (DP, MP).
  //! plot_commafy: whole number with thousands separators. plot_nice_tick: a round tick step giving
  //! ~10 ticks over a span. plot_xtics_list: an explicit gnuplot tics list over [lo,hi] at `step`.
  //! plot_lmargin_for_labels: left margin (in character units) that keeps a label of the given width
  //! inside the canvas when the tick font differs from the base font.
  string plot_commafy(int64_t v);
  int64_t plot_nice_tick(int64_t span);
  string plot_xtics_list(int64_t lo, int64_t hi, int64_t step);
  int plot_lmargin_for_labels(size_t max_label_chars, double label_font, double base_font);

  //! Predict Discordant Pair (DP) evidence.
  //
  //  Re-reads the DP candidate-region CSV written during identify_mutations
  //  (settings.dp_candidate_regions_file_name), builds a graph of candidate regions weighted by the
  //  number of discordant read pairs they share, and emits one DP GenomeDiff evidence item per region
  //  pair that shares >= 5 read pairs (no one-to-one matching: a region may pair with several others).
  //  Writes settings.dp_genome_diff_file_name.
  //
  //  Re-entrant: reads the CSV fresh each run and overwrites the output gd.
  void predict_discordant_pairs(const Settings& settings, Summary& summary, cReferenceSequences& ref_seq_info);

  //! Draw a per-side read-pair diagnostic plot (SVG) for each DP evidence item.
  //
  //  For each DP item and each side, re-fetches the reads counted as concordant/discordant at that
  //  side (using the same classification as the rescan) and renders a gnuplot SVG (arrows for reads,
  //  dashed connectors) into settings.evidence_path. Stamps the relative plot filename onto the DP
  //  entry (side_1/side_2) so cOutputEvidenceFiles can surface it as a per-side '?' evidence page.
  //  Call during the Output stage, before cOutputEvidenceFiles.
  void draw_discordant_pair_evidence_plots(const Settings& settings, Summary& summary, cReferenceSequences& ref_seq_info, cGenomeDiff& gd);

  //! Draw the run-wide concordant-pair crossing distribution plot and the per-sequence
  //  empirical-vs-projected overlay plots (from the tab files written during resolve_alignments).
  //  Call during the Output stage, before summary.html is generated.
  void draw_concordant_pair_crossing_plots(const Settings& settings, Summary& summary, cReferenceSequences& ref_seq_info);

  //! Human-readable description of which model the DP "skew" score used -- the empirical concordant-pair
  //  crossing distribution or the negative-binomial fallback -- and the number of non-deletion reference
  //  positions (N) that drove the choice. Empty string if unavailable. For the summary.html note.
  string dp_crossing_model_description(const Settings& settings, Summary& summary, cReferenceSequences& ref_seq_info);

} // namespace breseq

#endif
