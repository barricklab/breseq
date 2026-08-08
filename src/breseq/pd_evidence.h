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

#ifndef _BRESEQ_PD_EVIDENCE_H_
#define _BRESEQ_PD_EVIDENCE_H_

#include "common.h"
#include "settings.h"
#include "summary.h"
#include "reference_sequence.h"

namespace breseq {

  //! Load the null distribution that pair-distance (PD) evidence tests against.
  //
  //  Reads the persisted majority-orientation pair distance histogram
  //  (data/#.pair_distance_histogram.tab) and applies PD's length-bias correction: each bin is
  //  weighted by (d - trunc), with trunc = 2 * read_length. The reason is exact rather than
  //  heuristic. PD only ever looks at a pair through the question "does this pair's unsequenced
  //  middle gap cover position b?", and the pairs whose gap covers a FIXED point are sampled in
  //  proportion to the length of that gap. So the correct null for a covering pair is the raw
  //  fragment histogram reweighted by the gap length, which is what this returns.
  //
  //  Note what is deliberately NOT done: unlike the insert model dp_evidence builds from the same
  //  file, nothing is truncated from above at the discordant distance cutoff. For DP that truncation
  //  is right -- past the cutoff a pair is its own evidence. For PD it would be fatal, because every
  //  pair beyond the cutoff would then be indistinguishable from every other, and PD could never
  //  arbitrate against DP on a breakpoint they both see.
  //
  //  Shared with the seeding step in identify_mutations.cpp, so the null the candidate regions were
  //  found under and the null they are then judged under are the same by construction.
  //
  //  `weighted` is indexed by distance; `total` is its sum. Returns false (leaving both empty/zero)
  //  when the histogram is missing or carries no usable weight.
  bool pd_load_weighted_histogram(const string& file_name, int32_t trunc,
                                  vector<double>& weighted, double& total);

  //! Predict Pair Distance (PD) evidence.
  //
  //  PD finds the events that are invisible to DP because no single read pair is anomalous. A
  //  deletion of a few hundred bases makes every pair whose gap spans it map that much FARTHER apart
  //  than the library says it should; an insertion makes them map that much CLOSER. Neither shift
  //  need be large enough to cross the discordant distance cutoff -- and on the short side there is
  //  no cutoff at all -- so DP tests each pair and finds nothing. PD tests the collective shift of
  //  all the pairs covering one point instead.
  //
  //  Re-reads the PD candidate-region CSV written during identify_mutations
  //  (settings.pd_candidate_regions_file_name), rescans the BAM around each region's strongest
  //  column, estimates the size shift and the local variant frequency by maximum likelihood against
  //  the length-bias-corrected null, localizes the breakpoint to the interval every supporting
  //  pair's gap can accommodate, optionally snaps onto a validated split-read (JC) breakpoint inside
  //  that interval, applies the acceptance gates, and emits one two-sided PD GenomeDiff evidence
  //  item per surviving region. Writes settings.pd_genome_diff_file_name.
  //
  //  The two sides use JC's convention: side_1 = (b, -1) is the last retained base of the left flank
  //  and side_2 = (b + max(size_shift, 0) + 1, +1) the first retained base of the right flank, so
  //  they bracket exactly the reference bases the event removed and sit adjacent when it removed
  //  none. An insertion is reported by a negative size_shift with adjacent sides, because that is
  //  all PD can honestly say: a pair maps I bases closer whether I bases of novel sequence were
  //  inserted or k reference bases were replaced by k + I inserted ones.
  //
  //  Re-entrant: reads the CSV fresh each run and overwrites the output gd.
  void predict_pair_distances(const Settings& settings, Summary& summary, cReferenceSequences& ref_seq_info);

  //! Draw a read-pair diagnostic plot (SVG) for each PD evidence item.
  //
  //  Unlike DP -- whose two sides can be a whole chromosome apart, and which therefore draws each
  //  side separately and stitches them -- every pair PD draws covers the same point, so one
  //  continuous window holds the left mates, the right mates and the event between them. Pairs whose
  //  distance supports the call are drawn in one color, those arguing against it in another, and
  //  those the likelihood ratio cannot separate in a third, so the reader can see the two
  //  populations that produced the frequency. Stamps the relative plot filename onto the PD entry so
  //  cOutputEvidenceFiles can show it on the evidence page.
  //  Call during the Output stage, before cOutputEvidenceFiles.
  void draw_pair_distance_evidence_plots(const Settings& settings, Summary& summary, cReferenceSequences& ref_seq_info, cGenomeDiff& gd);

} // namespace breseq

#endif
