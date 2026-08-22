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

#ifndef _BRESEQ_MP_EVIDENCE_H_
#define _BRESEQ_MP_EVIDENCE_H_

#include "common.h"
#include "settings.h"
#include "summary.h"
#include "reference_sequence.h"

namespace breseq {

  //! Predict Missing Pair (MP) evidence.
  //
  //  MP flags a point where a NOVEL sequence -- one present in neither the reference nor any
  //  candidate junction -- has been inserted into the genome. Fragments spanning such a point put one
  //  mate in reference sequence, where it maps, and the other in the insert, where it maps nowhere.
  //  The mapped mates therefore pile up facing the insertion point: on the forward strand just to its
  //  left and on the reverse strand just to its right (under an FR library; RF is the mirror). This is
  //  precisely the case neither JC nor DP can see -- JC needs a split read whose two halves both map,
  //  DP needs both mates mapped.
  //
  //  A supporting read is one whose mate produced NO alignment at all -- not merely one flagged
  //  BAM_FMUNMAP, which also covers mates the aligner placed and breseq then rejected for falling
  //  short of --require-match-fraction. See kBreseqMateNeverAlignedBAMTag: conflating the two made
  //  the numerator a measure of alignment stringency, and on a real bacterial library 86% of the
  //  mates it counted turn out to align to the reference perfectly well.
  //
  //  Re-reads the MP candidate-region CSV written during identify_mutations
  //  (settings.mp_candidate_regions_file_name), refines each region's boundary with the same
  //  half-median-pair-distance statistics DP uses, scores each against the null carried on that
  //  CSV's `#key=value` header, applies the acceptance gates, and emits one one-sided MP GenomeDiff
  //  evidence item per surviving region. Writes settings.mp_genome_diff_file_name.
  //
  //  WHAT DECIDES A CALL is the score: minus log10 of the expected number of positions anywhere in
  //  the reference where this many of the reads on one strand would lose their mates by chance,
  //  given the genome-wide rate at which mates fail to align and how unevenly that rate is spread
  //  (a beta-binomial upper tail, times the number of independent windows). Everything else -- the
  //  read count, the distinct-start count, the local frequency -- is a local sanity check.
  //
  //  It has to be a measured null rather than an assumed one. The frequency gate this replaced was
  //  a Clopper-Pearson lower bound against a fixed 0.10, which is a precision bound, not a
  //  significance test: it controls sampling error only, so at large n it collapses onto the point
  //  estimate and more coverage makes a false positive MORE confident. Its implicit null was "the
  //  background is zero", true in a clean simulation and in nothing else. Measured on a 2x150
  //  E. coli library the background is 2.2% with Pearson dispersion phi = 4.8, and 3.9% of all
  //  windows in the genome already exceed 0.10 -- so the old cutoff sat below the 96th percentile
  //  of pure noise, and the run accepted 182 items of which essentially none were real.
  //
  //  Note the two denominators, which are deliberately different. MP_WINDOW_COUNT (every
  //  crossing-strand read on the kept flank) is what the score tests, because that is what the null
  //  was tabulated over. MP_TOTAL_COUNT (supporting + spanning) is what FREQUENCY reports, because
  //  "how much of the sample carries this" is a different question from "is it there at all".
  //
  //  Unlike DP there is no snapping step: with no mate coordinates there is nothing to snap onto.
  //
  //  Re-entrant: reads the CSV fresh each run and overwrites the output gd.
  void predict_missing_pairs(const Settings& settings, Summary& summary, cReferenceSequences& ref_seq_info);

  //! Draw a read diagnostic plot (SVG) for each MP evidence item.
  //
  //  For each MP item, re-fetches the reads counted at that position (using the same classification as
  //  the count) and renders a gnuplot SVG: one read per lane as an arrow in its mapping direction. A
  //  read whose mate mapped is drawn in blue with that mate and a dashed connector; a supporting read,
  //  whose mate did not map anywhere, is drawn in purple with its connector running past the placed
  //  position into the shaded region where the mate went missing. Stamps the relative plot filename
  //  onto the MP entry so cOutputEvidenceFiles can show it on the evidence page.
  //  Call during the Output stage, before cOutputEvidenceFiles.
  void draw_missing_pair_evidence_plots(const Settings& settings, Summary& summary, cReferenceSequences& ref_seq_info, cGenomeDiff& gd);

} // namespace breseq

#endif
