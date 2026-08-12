/*****************************************************************************

AUTHORS

  Jeffrey E. Barrick <jeffrey.e.barrick@gmail.com

LICENSE AND COPYRIGHT

  Copyright (c) 2008-2010 Michigan State University
  Copyright (c) 2011-2026 The University of Texas at Austin

  breseq is free software; you can redistribute it and/or modify it under the
  terms the GNU General Public License as published by the Free Software
  Foundation; either version 1, or (at your option) any later version.

*****************************************************************************/

#ifndef _BRESEQ_CN_EVIDENCE_H_
#define _BRESEQ_CN_EVIDENCE_H_

#include "common.h"

#include "genome_diff.h"
#include "reference_sequence.h"
#include "settings.h"

namespace breseq {

  // Runs the external tool CNery (https://github.com/barricklab/CNery) on the
  // current breseq output, then ingests its per-reference-sequence copy
  // number calls into the CN evidence genome diff files that the rest of the
  // pipeline (mutation prediction passthrough, HTML display) already expects.
  class CNEvidence
  {
  public:

    static void predict(
                        Settings& settings,
                        Summary& summary,
                        cReferenceSequences& ref_seq_info
                        );

    // Draws the coverage CNery actually made its calls from: one whole-reference overview per
    // sequence, plus a zoomed plot per CN entry stamped onto that entry as
    // "_cn_corrected_plot_file_name" for cOutputEvidenceFiles to surface.
    //
    // Must run during the Output step, for two independent reasons: the "_"-prefixed key it stamps
    // is never serialized to a .gd file (so it cannot be set back in stage 09 and read here), and
    // 09_copy_number_variation/ -- where CNery's CSVs live -- is deleted once Output finishes.
    // Missing CNery output is therefore a warning and a skip, never an error.
    static void draw_evidence_plots(
                                    const Settings& settings,
                                    cReferenceSequences& ref_seq_info,
                                    cGenomeDiff& gd
                                    );

  private:

    // One row of CNery's <prefix><seq_id>_CNV.csv: a single sliding window.
    struct cnery_window {
      int32_t start;          // win_st                (1-based, inclusive)
      int32_t end;            // win_end
      int32_t length;         // win_len
      double  raw_cov;        // norm_raw_cov          (normalized, uncorrected)
      double  corrected_cov;  // otr_gc_corr_norm_cov  (normalized, GC- and ori-ter-corrected)
      double  otr_fit_cov;    // otr_gc_corr_fact      (the ori-ter ramp that was divided out)
      int32_t copy_number;    // prob_copy_number      (HMM Viterbi state)
    };

    // The origin and terminus of replication CNery inferred, and used to build its OTR correction.
    struct cnery_otr {
      bool    detected;      // false => CNery found no ori-ter bias; the coordinates are meaningless
      int32_t origin;        // 1-based reference coordinate
      int32_t terminus;
      double  origin_cov;    // the fitted ramp's value at the origin ...
      double  terminus_cov;  // ... and at the terminus: the two ends of the straight line
      double  ratio;         // origin_cov / terminus_cov -- the magnitude of the bias
    };

    static void run_cnery(Settings& settings, Summary& summary, const string& cnery_output_prefix);

    // Reads CNery's <prefix><seq_id>_otr_results.json. Returns false (leaving otr.detected false)
    // if the file is missing/unparseable or CNery reported no ori-ter bias.
    static bool read_cnery_otr(const string& otr_file_name, cnery_otr& otr);

    // Parses the per-window CSV by COLUMN NAME (CNery's column order is not a contract). Returns
    // false if the file cannot be opened; asserts if it is there but lacks a needed column.
    // Optional fields get a neutral value when their column is absent -- see the definition.
    static bool read_cnery_windows(const string& cnv_file_name, vector<cnery_window>& windows);

    // Reduces a long window list to at most max_points entries, so a whole-genome overview does not
    // become a multi-megabyte SVG.
    static vector<cnery_window> bin_cnery_windows(const vector<cnery_window>& in, size_t max_points);

    // The ori-ter ramp at one position, evaluated in GENOMIC coordinates from the two endpoints
    // CNery reports. See the definition for why its per-window column cannot be drawn as a line.
    static double otr_ramp_at(const cnery_otr& otr, int32_t position, int32_t seq_length);

    static void ingest_csv_for_seq_id(
                                      const string& seq_id,
                                      const string& cnv_file_name,
                                      const string& break_pts_file_name,
                                      const string& gd_file_name
                                      );

    // Emits one gnuplot SVG over [plot_start, plot_end]. Windows outside that range are ignored.
    // A shaded_start/shaded_end narrower than the plot range greys out the flanks around it.
    // The ori/ter markers are drawn only where they fall inside the plotted range. seq_length is the
    // whole sequence, needed because the ori-ter ramp wraps around its end.
    static void render_cn_plot(
                               const string& output_svg,
                               const string& seq_id,
                               const vector<cnery_window>& windows,
                               int32_t plot_start,
                               int32_t plot_end,
                               int32_t shaded_start,
                               int32_t shaded_end,
                               const cnery_otr& otr,
                               int32_t seq_length
                               );

  }; // class CNEvidence

} // namespace breseq

#endif
