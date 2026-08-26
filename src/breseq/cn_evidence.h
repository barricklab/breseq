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
    //
    // The three coverage values are the SAME measurement after successive corrections --
    // raw_cov -> gc_corrected_cov -> corrected_cov -- which is what makes comparing their spread a
    // readout of how much each correction stage actually accomplished.
    struct cnery_window {
      int32_t start;             // win_st                (1-based, inclusive)
      int32_t end;               // win_end
      int32_t length;            // win_len
      double  gc_percent;        // gc_percent            (a FRACTION, 0-1, despite the name)
      double  raw_cov;           // norm_raw_cov          (normalized, uncorrected)
      double  gc_corrected_cov;  // gc_corr_norm_cov      (normalized, GC-corrected)
      double  gc_corr_fact;      // gc_corr_fact          (the LOWESS GC curve that was divided out)
      double  corrected_cov;     // otr_gc_corr_norm_cov  (normalized, GC- and ori-ter-corrected)
      double  otr_fit_cov;       // otr_gc_corr_fact      (the ori-ter ramp that was divided out)
      int32_t copy_number;       // prob_copy_number      (HMM Viterbi state)
    };

    // One merged run of equal copy number, from CNery's <prefix><seq_id>_break_pts.csv. CNery has
    // already collapsed contiguous windows of the same HMM state, so these are the copy-number
    // segments in genomic coordinates -- real boundaries, with no window granularity left to
    // resolve and no gap where CNery dropped a window. Both the CN evidence entries and the step
    // line on the CN plots are built from these, which is what keeps the two agreeing.
    // Inclusive on both ends, and 1-based like every other reference coordinate in breseq --
    // read_cnery_segments() validates that as it reads and rejects anything else, so nothing
    // downstream has to re-check it.
    //
    // CNery before its "Open the first segment of a sequence at 1, not 0" fix wrote the very first
    // segment of every sequence with Startpos 0, and breseq passed it straight through on the theory
    // that the first segment is always copy number 1 and so always dropped by the evidence ingest.
    // It is not: a run whose first window is called CN != 1 -- a deletion at the start of a contig,
    // or a library so thin that CNery calls the whole genome CN 0 -- wrote start = 0 into the .gd and
    // then died parsing it back at the Output stage.
    struct cnery_segment {
      int32_t start;
      int32_t end;           // Startpos + Segment_Size - 1
      int32_t copy_number;   // CNery's State
    };

    // The origin and terminus of replication CNery inferred, and used to build its OTR correction.
    //
    // The four fields below the fit are read even when there is no fit -- they are what says WHY
    // there is none, which "detected == false" on its own does not.
    struct cnery_otr {
      bool    detected;      // false => CNery found no ori-ter bias; the coordinates are meaningless
      int32_t origin;        // 1-based reference coordinate
      int32_t terminus;
      double  origin_cov;    // the fitted ramp's value at the origin ...
      double  terminus_cov;  // ... and at the terminus: the two ends of the straight line
      double  ratio;         // origin_cov / terminus_cov -- the magnitude of the bias

      // This sequence's coverage relative to the LONGEST sequence of the run, which reads exactly
      // 1.0. Deliberately non-integral -- a plasmid at 2.96x is a measurement, and rounding it to 3
      // would throw away the precision that makes it worth reporting. 0 => CNery did not report it.
      double  relative_copy_number;
      // How the ori/ter above were arrived at ("Ori-ter coordinates fit by coverage", the GC-skew
      // method string, "No usable coverage", or "No OTR correction (--bias gc|none)") ...
      string  correction_type;
      // ... and which arm supplied them: "coverage fit", "GC skew" or "not corrected".
      string  breakpoint_source;
      // Non-empty => the sequence had nothing to measure at all, and says which way: no position
      // rows in the coverage table, or every window at zero coverage. That is a different statement
      // from "no ori-ter bias", which is what an empty value here leaves it as.
      string  no_coverage_reason;
    };

    // The per-reference-sequence scale every plotted coverage series is drawn on.
    //
    // CNery normalizes against ONE pooled median across every reference sequence it was handed
    // (core.py process_multi_genome: norm_raw_cov = read_count_cov / global_median), so on a
    // multi-copy replicon all three coverage columns arrive at a MULTIPLE of single copy -- 2.38x
    // and 5.11x on two real plasmids -- while CNery's HMM refits the single-copy level per sequence
    // and calls both copy number 1. Plotted as they arrive, the coverage traces float far above the
    // copy-number line they exist to be read against, and the GC-bias plot's y-axis clip leaves a
    // multi-copy plasmid's cloud mostly off the top of the figure. CNery has the same problem
    // internally and solves it the same way (core.py censored_median_coverage, _in_band_fractions).
    //
    // TWO divisors, not one per series, and not one for everything.
    //
    // raw_cov, gc_corrected_cov and otr_fit_cov all arrive on CNery's pooled scale, so ONE divisor
    // covers them and every offset BETWEEN them survives -- which matters, because those offsets are
    // the corrections. On one real plasmid the GC correction moves the sequence's level from 3.90x
    // the chromosome to 5.12x: dividing each series by its own median would flatten that 25% out of
    // the picture, when showing it is most of the point of drawing the uncorrected trace at all.
    //
    // corrected_cov is the exception, and CNery's doing: its OTR stage renormalizes to ~1 wherever a
    // tent actually fires (core.py otr_fit: y_corr = y/y_fit) and leaves the value untouched where
    // one does not (that branch returns y). One divisor could not put it at single copy in both
    // cases. Where no tent fired the two divisors are equal anyway -- corrected_cov IS
    // gc_corrected_cov there -- so this splits only where CNery already split it.
    //
    // `gc` is CNery's own censored_median_coverage: the median of the GC-corrected coverage, which
    // is what its "Relative copy number" is built from and the closest thing the run has to a
    // definition of this sequence's single-copy level. 1.0 is the neutral value, so a series with
    // nothing to measure is left exactly as it arrives.
    struct cn_plot_scale {
      double gc;          // divides raw_cov, gc_corrected_cov, otr_fit_cov, ori/ter marker heights
      double corrected;   // divides corrected_cov, and equals gc wherever no ori-ter tent fired
      cn_plot_scale() : gc(1.0), corrected(1.0) {}
    };

    static void run_cnery(Settings& settings, Summary& summary, cReferenceSequences& ref_seq_info, const string& cnery_output_prefix);

    // Reads CNery's <prefix><seq_id>_otr_results.json. Returns false (leaving otr.detected false)
    // if the file is missing/unparseable or CNery reported no ori-ter bias.
    static bool read_cnery_otr(const string& otr_file_name, cnery_otr& otr);

    // Parses the per-window CSV by COLUMN NAME (CNery's column order is not a contract). Returns
    // false if the file cannot be opened; asserts if it is there but lacks a needed column.
    // Optional fields get a neutral value when their column is absent -- see the definition.
    static bool read_cnery_windows(const string& cnv_file_name, vector<cnery_window>& windows);

    // Parses CNery's <prefix><seq_id>_break_pts.csv (Startpos,State,Segment_Size), positionally --
    // unlike the per-window CSV this file has exactly three columns and the assert is fatal.
    // Returns false if the file cannot be opened, so each caller can decide: ingesting the evidence
    // treats that as fatal, drawing the plots treats it as a reason to leave the line off.
    //
    // Coordinates arrive as breseq's: 1-based and inclusive, with anything that cannot be rejected
    // rather than corrected. sequence_length bounds the far end; pass 0 to skip that one check.
    static bool read_cnery_segments(const string& break_pts_file_name, vector<cnery_segment>& segments,
                                    int32_t sequence_length);

    // Reduces a long window list to at most max_points entries, so a whole-genome overview does not
    // become a multi-megabyte SVG.
    static vector<cnery_window> bin_cnery_windows(const vector<cnery_window>& in, size_t max_points);

    // The ori-ter ramp at one position, evaluated in GENOMIC coordinates from the two endpoints
    // CNery reports. See the definition for why its per-window column cannot be drawn as a line.
    static double otr_ramp_at(const cnery_otr& otr, int32_t position, int32_t seq_length);

    // The windows a per-sequence statistic is measured over: every correction stage positive, then
    // restricted to copy number 1 when there are at least 100 of those. Shared by summarize() and
    // compute_plot_scale() so that the spread reported in summary.html and the scale the plots are
    // drawn on cannot drift into describing different parts of the sequence.
    static void select_measured_windows(const vector<cnery_window>& windows,
                                        vector<size_t>& selected,
                                        bool& single_copy_only);

    // The median of each coverage series over those windows -- see cn_plot_scale for why the plots
    // need them. Must be computed over the FULL UNBINNED window list of the whole sequence: binning
    // discards the HMM state this selects on, and a per-CN-item plot covers a subset that is
    // typically inside an amplification, so neither can be allowed to derive a scale of its own.
    static cn_plot_scale compute_plot_scale(const vector<cnery_window>& windows);

    // Distills the fit and the per-window coverage into the numbers summary.html and summary.json
    // report. Must happen in stage 09: everything it reads is deleted when the pipeline finishes.
    static void summarize(
                          const cnery_otr& otr,
                          const vector<cnery_window>& windows,
                          CopyNumberSummary& cns
                          );

    static void ingest_csv_for_seq_id(
                                      const string& seq_id,
                                      const vector<cnery_window>& windows,
                                      const string& break_pts_file_name,
                                      const string& gd_file_name,
                                      int32_t sequence_length
                                      );

    // Emits one gnuplot SVG over [plot_start, plot_end]. Windows outside that range are ignored.
    // A shaded_start/shaded_end narrower than the plot range greys out the flanks around it.
    // The ori/ter markers are drawn only where they fall inside the plotted range. seq_length is the
    // whole sequence, needed because the ori-ter ramp wraps around its end.
    // `segments` drives the copy-number step line and `windows` the coverage traces. They are kept
    // separate on purpose: the coverage is a per-window measurement, while the copy number is a
    // piecewise-constant call whose boundaries CNery already resolved to a base.
    static void render_cn_plot(
                               const string& output_svg,
                               const string& seq_id,
                               const vector<cnery_window>& windows,
                               const vector<cnery_segment>& segments,
                               int32_t plot_start,
                               int32_t plot_end,
                               int32_t shaded_start,
                               int32_t shaded_end,
                               const cnery_otr& otr,
                               int32_t seq_length,
                               const cn_plot_scale& scale
                               );

    // Coverage against window GC content, before and after the GC correction. Drawn from CNery's own
    // gc_corr_fact column rather than from its GC_bias PDF, which is pooled across every reference
    // sequence and whose "LOWESS fit" line is really a degree-2 polyfit through the correction
    // factors rather than the LOWESS curve that was actually divided out.
    //
    // Does nothing if CNery emitted no GC columns (--bias none or --bias otr).
    // relative_copy_number is only used for the label saying what scale the axis is on, and is
    // CNery's own (cnery_otr::relative_copy_number); 0 leaves that half of the label off.
    static void render_gc_bias_plot(
                                    const string& output_svg,
                                    const string& seq_id,
                                    const vector<cnery_window>& windows,
                                    const cn_plot_scale& scale,
                                    double relative_copy_number
                                    );

  }; // class CNEvidence

} // namespace breseq

#endif
