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

} // namespace breseq

#endif
