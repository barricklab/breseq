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

#include "cn_evidence.h"

#include "genome_diff.h"

using namespace std;

namespace breseq {

void CNEvidence::predict(Settings& settings, Summary& summary, cReferenceSequences& ref_seq_info)
{
  if (settings.installed.count("cnery") == 0) {
    settings.installed["cnery"] = which("CNery");
  }
  if (settings.installed["cnery"].size() == 0) {
    ERROR("Could not find 'CNery' command in $PATH.\n"
          "Install it (e.g. 'pip install CNery') to use --predict-copy-number.\n"
          "See https://github.com/barricklab/CNery");
  }

  run_cnery(settings, summary, ref_seq_info, settings.cnery_output_path);

  const set<string> analyzed_seq_ids = settings.call_mutations_seq_id_set();

  for (cReferenceSequences::iterator it = ref_seq_info.begin(); it != ref_seq_info.end(); ++it) {

    cAnnotatedSequence& seq = *it;

    // Junction-only references are never pileup'd, so they have no coverage table and CNery was
    // never given one for them. See run_cnery().
    if (analyzed_seq_ids.count(seq.m_seq_id) == 0) continue;

    string cnv_file_name = settings.file_name(settings.cnery_cnv_csv_file_name, "@", seq.m_seq_id);
    string break_pts_file_name = settings.file_name(settings.cnery_break_points_file_name, "@", seq.m_seq_id);
    string gd_file_name = settings.file_name(settings.copy_number_evidence_genome_diff_file_name, "@", seq.m_seq_id);

    vector<cnery_window> windows;
    ASSERT(read_cnery_windows(cnv_file_name, windows),
           "Could not open CNery output file: " + cnv_file_name);

    ingest_csv_for_seq_id(seq.m_seq_id, windows, break_pts_file_name, gd_file_name,
                          static_cast<int32_t>(seq.m_length));

    // Recorded now because none of what it is derived from survives the run: this whole directory is
    // deleted once the pipeline completes, but Output still has to report how the analysis went.
    cnery_otr otr;
    read_cnery_otr(settings.file_name(settings.cnery_otr_results_file_name, "@", seq.m_seq_id), otr);
    summarize(otr, windows, summary.copy_number[seq.m_seq_id]);
  }
}

// CNery takes coverage tables and nothing else: no BAM, no reference FASTA, no breseq output folder.
// It reads the reference sequence out of the table's own ref_base column, and derives each sequence
// id from the table's file name (basename minus ".coverage.tsv") -- which is what makes its output
// land at the cnery_<seq_id>_* paths Settings already expects.
void CNEvidence::run_cnery(Settings& settings, Summary& summary, cReferenceSequences& ref_seq_info, const string& cnery_output_prefix)
{
  (void)summary;
  // Example of how to use read length as the fragment length – but this is not appropriate
  // Use
  //uint32_t fragment_length(0);
  //float total_bases(0);
  //float total_reads(0);
  //for(auto& rf : settings.read_file_sets.flat_files())
  //{
  //  const AnalyzeFastqSummary& s = summary.sequence_conversion.reads[it->m_base_name];
  //  total_bases += static_cast<double>(s.num_bases);
  //  total_reads += static_cast<double>(s.num_reads);
  //}
  //fragment_length = trunc(total_bases/total_reads);

  string command = double_quote(settings.installed["cnery"]);

  // Name every table explicitly rather than handing over 08_mutation_identification/. A directory
  // argument would make the input set whatever happens to match CNery's file endings in there;
  // listing the files means a missing one is an error instead of a silently smaller analysis.
  //
  // The set is the one stage 08 pileups, NOT every reference sequence: a junction-only reference
  // (-s/--junction-only-reference) is deliberately excluded from mutation calling, so no coverage
  // table is written for it. Naming one anyway made CNery raise FileNotFoundError before it read a
  // single table, and SYSTEM() below does not ignore errors -- one junction-only reference took the
  // whole run down. Copy number on a reference used only to call junctions means nothing anyway.
  const set<string> analyzed_seq_ids = settings.call_mutations_seq_id_set();
  for (cReferenceSequences::iterator it = ref_seq_info.begin(); it != ref_seq_info.end(); ++it) {
    if (analyzed_seq_ids.count(it->m_seq_id) == 0) continue;
    command += " " + double_quote(settings.file_name(settings.complete_coverage_text_file_name, "@", it->m_seq_id));
  }

  command += " -o " + double_quote(cnery_output_prefix);
  //command += " --frag-size " + to_string<uint32_t>(fragment_length);

  SYSTEM(command, false, false, false);
}

// CNery's <prefix><seq_id>_CNV.csv has one row per sliding window. Columns are matched by NAME,
// not position: CNery builds the frame by successively adding columns across its GC- and
// ori-ter-correction stages, so the order is an implementation detail rather than a contract.
//
// Only win_st, win_len and otr_gc_corr_norm_cov are required -- those are what turning the
// segments into CN entries needs. The rest are what the plots want, and get a neutral value when
// absent (0 for coverage, -1 = "no call" for copy number) so a CNery that stopped emitting them
// degrades the plot instead of failing the run.
bool CNEvidence::read_cnery_windows(const string& cnv_file_name, vector<cnery_window>& windows)
{
  windows.clear();

  ifstream cnv_file(cnv_file_name.c_str());
  if (!cnv_file.good()) return false;

  string header_line;
  getline(cnv_file, header_line);
  vector<string> header_fields = split(header_line, ",");

  size_t win_st_col = string::npos, win_len_col = string::npos, rel_cov_col = string::npos;
  size_t win_end_col = string::npos, raw_cov_col = string::npos, copy_number_col = string::npos;
  size_t otr_fit_col = string::npos;
  size_t gc_percent_col = string::npos, gc_cov_col = string::npos, gc_fit_col = string::npos;
  for (size_t i = 0; i < header_fields.size(); i++) {
    if (header_fields[i] == "win_st") win_st_col = i;
    else if (header_fields[i] == "win_len") win_len_col = i;
    else if (header_fields[i] == "otr_gc_corr_norm_cov") rel_cov_col = i;
    else if (header_fields[i] == "win_end") win_end_col = i;
    else if (header_fields[i] == "norm_raw_cov") raw_cov_col = i;
    else if (header_fields[i] == "prob_copy_number") copy_number_col = i;
    else if (header_fields[i] == "otr_gc_corr_fact") otr_fit_col = i;
    else if (header_fields[i] == "gc_percent") gc_percent_col = i;
    else if (header_fields[i] == "gc_corr_norm_cov") gc_cov_col = i;
    else if (header_fields[i] == "gc_corr_fact") gc_fit_col = i;
  }
  ASSERT((win_st_col != string::npos) && (win_len_col != string::npos) && (rel_cov_col != string::npos),
         "Unexpected column layout in CNery output file: " + cnv_file_name);

  // Highest column index actually read, so a truncated row can be skipped rather than indexed past.
  size_t max_needed_col = 0;
  const size_t used_cols[] = { win_st_col, win_len_col, rel_cov_col, win_end_col, raw_cov_col,
                               copy_number_col, otr_fit_col, gc_percent_col, gc_cov_col, gc_fit_col };
  for (size_t i = 0; i < sizeof(used_cols) / sizeof(used_cols[0]); i++) {
    if ((used_cols[i] != string::npos) && (used_cols[i] > max_needed_col)) max_needed_col = used_cols[i];
  }

  string line;
  while (getline(cnv_file, line)) {
    if (line.size() == 0) continue;
    vector<string> fields = split(line, ",");
    if (fields.size() <= max_needed_col) continue;

    cnery_window w;
    w.start = from_string<int32_t>(fields[win_st_col]);
    w.length = from_string<int32_t>(fields[win_len_col]);
    // CNery computes win_len as win_end - win_st, so the fallback needs no +/-1 correction.
    w.end = (win_end_col != string::npos) ? from_string<int32_t>(fields[win_end_col]) : (w.start + w.length);
    w.corrected_cov = from_string<double>(fields[rel_cov_col]);
    w.raw_cov = (raw_cov_col != string::npos) ? from_string<double>(fields[raw_cov_col]) : 0.0;
    w.otr_fit_cov = (otr_fit_col != string::npos) ? from_string<double>(fields[otr_fit_col]) : 0.0;
    w.gc_percent = (gc_percent_col != string::npos) ? from_string<double>(fields[gc_percent_col]) : 0.0;
    w.gc_corrected_cov = (gc_cov_col != string::npos) ? from_string<double>(fields[gc_cov_col]) : 0.0;
    w.gc_corr_fact = (gc_fit_col != string::npos) ? from_string<double>(fields[gc_fit_col]) : 0.0;
    w.copy_number = (copy_number_col != string::npos) ? from_string<int32_t>(fields[copy_number_col]) : -1;
    windows.push_back(w);
  }
  cnv_file.close();

  return true;
}

// The headline ratio's key. CNery spelled it "Termius" from the day it was introduced until its
// 1.1.0 release (its commit 98aa6ba), and kept the misspelling on purpose because this reader looks
// the key up by name -- correcting it there would have silently cost breseq all of its ori-ter
// reporting. It is now spelled correctly, so read that first and fall back to the old spelling: a
// user with an older CNery on $PATH still gets their bias reported rather than a silent "none".
static const char* kOTRRatioKey       = "Origin-to-Terminus/Bias Ratio";
static const char* kOTRRatioLegacyKey = "Origin-to-Termius/Bias Ratio";

// CNery's <prefix><seq_id>_otr_results.json records the origin and terminus of replication it fit,
// which is what its ori-ter coverage correction is built from.
//
// Two traps in that file. The keys are named "Origin window"/"Terminus window" but hold 1-based
// GENOMIC COORDINATES, not window indices -- CNery converts them with `df["win_st"].iloc[ori_idx]`
// before writing (core.py, apply_otr_correction). And when it fits no ramp it still writes both
// keys -- with the first and last coordinate of the sequence, or 0 when the sequence has no windows
// at all -- plus a "Not detected" string for the ratio; those placeholders must not be drawn as if
// they were a real ori/ter, which is what the ratio's type is used to distinguish here.
//
// Returns false when there is no fit, but the diagnostic fields are filled in BEFORE that: they are
// exactly what says why there is none, and both callers keep the record either way.
bool CNEvidence::read_cnery_otr(const string& otr_file_name, cnery_otr& otr)
{
  otr.detected = false;
  otr.origin = 0;
  otr.terminus = 0;
  otr.origin_cov = 0.0;
  otr.terminus_cov = 0.0;
  otr.ratio = 0.0;
  otr.relative_copy_number = 0.0;
  otr.correction_type = "";
  otr.breakpoint_source = "";
  otr.no_coverage_reason = "";

  ifstream otr_file(otr_file_name.c_str());
  if (!otr_file.good()) return false;

  json j;
  try {
    otr_file >> j;
  } catch (...) {
    WARN("Could not parse CNery output file: " + otr_file_name);
    return false;
  }

  // Every one of these can be JSON null -- CNery maps its own non-finite values that way so the file
  // stays strict RFC JSON -- so each is guarded on actually having the type it is read as.
  if (j.count("Relative copy number") && j["Relative copy number"].is_number())
    otr.relative_copy_number = j["Relative copy number"].get<double>();
  if (j.count("Correction type") && j["Correction type"].is_string())
    otr.correction_type = j["Correction type"].get<string>();
  if (j.count("Breakpoint source") && j["Breakpoint source"].is_string())
    otr.breakpoint_source = j["Breakpoint source"].get<string>();
  // Written on every record as of CNery 1.1.0, null on a healthy one, so a reader never has to tell
  // "this sequence had coverage" from "this CNery predates the key".
  if (j.count("No usable coverage reason") && j["No usable coverage reason"].is_string())
    otr.no_coverage_reason = j["No usable coverage reason"].get<string>();

  if (!j.count("Origin window") || !j.count("Terminus window")) return false;

  // "Not detected" (a string) rather than a number means CNery found no bias to correct. A null is
  // possible too, when the terminus anchor came back at zero and CNery declined to invent a ratio.
  const char* ratio_key = j.count(kOTRRatioKey) ? kOTRRatioKey
                        : (j.count(kOTRRatioLegacyKey) ? kOTRRatioLegacyKey : NULL);
  if ((ratio_key == NULL) || !j[ratio_key].is_number()) return false;

  otr.origin = j["Origin window"].get<int32_t>();
  otr.terminus = j["Terminus window"].get<int32_t>();
  otr.ratio = j[ratio_key].get<double>();

  // The two ends of the fitted ramp. NaN in CNery becomes JSON null, so these are only usable when
  // they really are numbers -- the markers that use them are guarded on that separately.
  if (j.count("Origin coverage (normalized)") && j["Origin coverage (normalized)"].is_number())
    otr.origin_cov = j["Origin coverage (normalized)"].get<double>();
  if (j.count("Terminus coverage (normalized)") && j["Terminus coverage (normalized)"].is_number())
    otr.terminus_cov = j["Terminus coverage (normalized)"].get<double>();

  otr.detected = true;

  return true;
}

// CNery's <prefix><seq_id>_break_pts.csv already merges contiguous runs of the same predicted copy
// number into ranges for us (Startpos,State,Segment_Size), so nothing here has to re-implement that.
//
// Read positionally rather than by column name, unlike the per-window CSV: this file has exactly
// three columns and a fourth would mean CNery changed something we need to look at, so the assert is
// deliberately fatal. Returns false only when the file cannot be opened at all -- what that means
// differs by caller, so it is theirs to decide.
//
// This is the ONE place CNery's coordinates become breseq's, so it is also where they are checked,
// and nothing here corrects them -- a coordinate breseq cannot use is an error, not something to
// clamp. The alternative is what used to happen: a bad coordinate rides through
// ingest_csv_for_seq_id() into a .gd that breseq itself then refuses to parse, and the run dies in
// cGenomeDiff::read() at the Output stage with a stack trace five stages away from the file that
// actually caused it.
//
// sequence_length is this reference's length, used only to keep a segment from running off the end
// of it; pass 0 to skip that one check.
bool CNEvidence::read_cnery_segments(const string& break_pts_file_name, vector<cnery_segment>& segments, int32_t sequence_length)
{
  segments.clear();

  ifstream break_pts_file(break_pts_file_name.c_str());
  if (!break_pts_file.good()) return false;

  string line;
  getline(break_pts_file, line); // discard header

  bool warned_past_end = false;
  size_t row = 0;

  while (getline(break_pts_file, line)) {
    if (line.size() == 0) continue;
    vector<string> fields = split(line, ",");
    ASSERT(fields.size() == 3, "Unexpected number of columns in CNery output file: " + break_pts_file_name);
    row++;

    // Startpos is where the segment opens and Segment_Size is its length, so the last base it covers
    // is one before the next segment's start.
    int32_t raw_start = from_string<int32_t>(fields[0]);

    cnery_segment s;
    s.copy_number = from_string<int32_t>(fields[1]);
    s.start       = raw_start;
    s.end         = raw_start + from_string<int32_t>(fields[2]) - 1;

    // Reference coordinates are 1-based and inclusive. Both of these produce a CN entry that fails
    // GenomeDiff validation, so they are fatal here rather than at Output: START and END are
    // positive-integer fields, and an entry that ends before it starts is not a range.
    //
    // A start of 0 specifically is what CNery used to write on the first segment of every sequence
    // -- _segments_from_path() in its core.py seeded start_pos = 0 before it had looked at a window,
    // so it was unconditional and had nothing to do with what state that segment was in. breseq
    // passed it through, on the theory that the first segment is always copy number 1 and so always
    // dropped by the evidence ingest; a first window called CN != 1 wrote start = 0 into the .gd and
    // killed the run at Output. Fixed in CNery by its commit "Open the first segment of a sequence
    // at 1, not 0", and NOT worked around here: dev-environment.yml pins cnery-prerelease to an
    // exact commit build, so there is one CNery a given breseq tree is ever run against and a
    // tolerated 0 would only hide a pin that had slipped backwards. If this fires, that is what
    // happened -- check the pin rather than adding a clamp.
    ASSERT(s.start >= 1,
           "Segment starting before the first base of the sequence (" + to_string<int32_t>(s.start) +
           ") on row " + to_string<size_t>(row) + " of CNery output file: " + break_pts_file_name +
           "\nA start of 0 on row 1 means CNery predates its \"Open the first segment of a sequence\n"
           "at 1, not 0\" fix; check the cnery-prerelease pin in dev-environment.yml.");
    ASSERT(s.end >= s.start,
           "Segment ending (" + to_string<int32_t>(s.end) + ") before it starts (" +
           to_string<int32_t>(s.start) + ") on row " + to_string<size_t>(row) +
           " of CNery output file: " + break_pts_file_name);

    // Running past the end of the contig is recoverable in a way the two above are not -- the .gd
    // parses either way and every plot clamps to its own range -- so trim it and say so instead of
    // ending the run. CNery stops its last window short of the sequence end rather than overrunning
    // it, so this is not expected to fire; it is here so that if that ever changes, the CN evidence
    // stays inside the reference it describes. Warn once: only the last row can realistically trip.
    if ((sequence_length > 0) && (s.end > sequence_length)) {
      if (!warned_past_end) {
        WARN("Copy number segment ending at " + to_string<int32_t>(s.end) + " runs past the end of the " +
             to_string<int32_t>(sequence_length) + " bp reference sequence; trimming it.\n"
             "CNery output file: " + break_pts_file_name);
        warned_past_end = true;
      }
      s.end = sequence_length;
      if (s.start > s.end) continue;
    }

    segments.push_back(s);
  }
  break_pts_file.close();

  return true;
}

// Turn each merged range into a CN evidence entry, using the per-window file only to compute a
// representative relative coverage value to display for that range.
void CNEvidence::ingest_csv_for_seq_id(
                                       const string& seq_id,
                                       const vector<cnery_window>& windows,
                                       const string& break_pts_file_name,
                                       const string& gd_file_name,
                                       int32_t sequence_length
                                       )
{
  uint32_t window_size = windows.size() ? static_cast<uint32_t>(windows[0].length) : 0;

  vector<cnery_segment> segments;
  ASSERT(read_cnery_segments(break_pts_file_name, segments, sequence_length),
         "Could not open CNery output file: " + break_pts_file_name);

  cGenomeDiff gd;

  for (size_t seg_i = 0; seg_i < segments.size(); seg_i++) {
    int32_t start_pos = segments[seg_i].start;
    int32_t end_pos = segments[seg_i].end;
    int32_t copy_number = segments[seg_i].copy_number;

    // Copy number 1 is the baseline (haploid, single-copy) state -- only
    // regions CNery calls as different from that are evidence-worthy.
    if (copy_number == 1) continue;

    double coverage_sum = 0;
    uint32_t coverage_n = 0;
    for (size_t i = 0; i < windows.size(); i++) {
      if ((windows[i].start >= start_pos) && (windows[i].start <= end_pos)) {
        coverage_sum += windows[i].corrected_cov;
        coverage_n++;
      }
    }

    cDiffEntry item(CN);
    item[SEQ_ID] = seq_id;
    item[START] = to_string<int32_t>(start_pos);
    item[END] = to_string<int32_t>(end_pos);
    item["tile_size"] = to_string<uint32_t>(window_size);
    item["copy_number"] = to_string<int32_t>(copy_number);
    // Three significant figures, deliberately. This number is CNery's, fitted with scipy and
    // statsmodels, and its trailing digits move whenever that stack is re-resolved -- 0.035776 became
    // 0.0357757 on long_ltee_ara_m3_32k_mp2800 from nothing but a fresh conda env, failing the test on
    // a difference of 2e-7. Nothing consumes that precision: the only reader is the CN evidence page,
    // which renders this at two decimal places (output.cpp), and mutation prediction reads
    // copy_number, never this. Three digits leaves about a thousandfold margin over that drift.
    //
    // Default float format, not fixed or scientific: fixed(3) would crush the small values that
    // matter most here (0.00117571 -> 0.001), and scientific would turn the clean 0 that a fully
    // deleted region gets into 0.00e+00.
    ostringstream relative_coverage;
    relative_coverage << setprecision(3) << ((coverage_n > 0) ? (coverage_sum / coverage_n) : 0.0);
    item["relative_coverage"] = relative_coverage.str();

    gd.add(item);
  }

  gd.write(gd_file_name);
}

// ---------------------------------------------------------------------------------------------
// Distilling the run into the numbers reported in summary.html and summary.json
// ---------------------------------------------------------------------------------------------

// Sorts in place; the callers below have no use for the original order.
static double cnery_median(vector<double>& v)
{
  if (v.empty()) return 0.0;
  sort(v.begin(), v.end());
  size_t n = v.size();
  return (n % 2) ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
}

// Robust coefficient of variation: 1.4826 * MAD / median. The constant makes the MAD a consistent
// estimator of the standard deviation for normally distributed data, so this reads on the same scale
// as an ordinary CV while ignoring the handful of extreme windows that CNery's coverage clipping and
// the ragged ends of a reference sequence leave behind. Takes its argument by value: it reorders it.
static double cnery_robust_cv(vector<double> v)
{
  double med = cnery_median(v);
  if (med <= 0.0) return 0.0;
  for (size_t i = 0; i < v.size(); i++) v[i] = fabs(v[i] - med);
  return 1.4826 * cnery_median(v) / med;
}

void CNEvidence::summarize(const cnery_otr& otr, const vector<cnery_window>& windows, CopyNumberSummary& cns)
{
  cns.otr_detected = otr.detected;
  cns.origin = otr.origin;
  cns.terminus = otr.terminus;
  cns.origin_coverage = otr.origin_cov;
  cns.terminus_coverage = otr.terminus_cov;
  cns.otr_ratio = otr.ratio;
  cns.relative_copy_number = otr.relative_copy_number;
  cns.correction_type = otr.correction_type;
  cns.breakpoint_source = otr.breakpoint_source;
  cns.no_coverage_reason = otr.no_coverage_reason;

  if (windows.empty()) return;

  cns.window_size = windows[0].length;

  // The step, as the most common gap between consecutive window starts rather than the mean gap:
  // CNery drops every window overlapping repeat coverage, and those holes -- up to 5.9 kb on REL606
  // -- would drag a mean well above the step actually used.
  {
    map<int32_t, uint32_t> step_counts;
    for (size_t i = 1; i < windows.size(); i++) {
      int32_t step = windows[i].start - windows[i - 1].start;
      if (step > 0) step_counts[step]++;
    }
    uint32_t best_count = 0;
    for (map<int32_t, uint32_t>::const_iterator it = step_counts.begin(); it != step_counts.end(); it++) {
      if (it->second > best_count) { best_count = it->second; cns.window_step = it->first; }
    }
  }

  // How large the GC correction was. A range that barely straddles 1.0 means there was almost no GC
  // bias to remove, so a dramatic-looking correction would deserve suspicion.
  bool have_gc_factor = false;
  for (size_t i = 0; i < windows.size(); i++) {
    double f = windows[i].gc_corr_fact;
    if (f <= 0.0) continue;
    if (!have_gc_factor) {
      cns.gc_correction_min = cns.gc_correction_max = f;
      have_gc_factor = true;
    } else {
      if (f < cns.gc_correction_min) cns.gc_correction_min = f;
      if (f > cns.gc_correction_max) cns.gc_correction_max = f;
    }
  }

  // The three spreads are only meaningful against each other, so they must come from ONE set of
  // windows. Requiring all three values to be positive is what makes that possible: CNery leaves
  // otr_gc_corr_norm_cov at 0 wherever raw coverage fell below 10% of the median, so scoring each
  // stage over "its own" non-zero windows would compare different parts of the genome.
  vector<size_t> usable, single_copy;
  for (size_t i = 0; i < windows.size(); i++) {
    const cnery_window& w = windows[i];
    if ((w.raw_cov <= 0.0) || (w.gc_corrected_cov <= 0.0) || (w.corrected_cov <= 0.0)) continue;
    usable.push_back(i);
    if (w.copy_number == 1) single_copy.push_back(i);
  }

  // Restricting to single-copy windows is the point: a real amplification or deletion inflates the
  // spread at every stage equally, which would hide exactly the improvement being measured. Fall
  // back to all windows when there are too few calls to select on -- either because the HMM column
  // was missing, or because so much of the sequence is non-single-copy that the subset is not
  // representative -- and record which happened so the report can say so.
  const size_t kMinimumSpreadWindows = 100;
  const vector<size_t>& selected = (single_copy.size() >= kMinimumSpreadWindows) ? single_copy : usable;
  cns.spread_single_copy = (single_copy.size() >= kMinimumSpreadWindows);
  cns.spread_windows = selected.size();
  if (selected.empty()) return;

  vector<double> raw, gc, otr_gc;
  raw.reserve(selected.size());
  gc.reserve(selected.size());
  otr_gc.reserve(selected.size());
  for (size_t k = 0; k < selected.size(); k++) {
    const cnery_window& w = windows[selected[k]];
    raw.push_back(w.raw_cov);
    gc.push_back(w.gc_corrected_cov);
    otr_gc.push_back(w.corrected_cov);
  }

  cns.cv_uncorrected = cnery_robust_cv(raw);
  cns.cv_gc_corrected = cnery_robust_cv(gc);
  cns.cv_otr_gc_corrected = cnery_robust_cv(otr_gc);
}

// ---------------------------------------------------------------------------------------------
// Plots of CNery's corrected coverage
//
// These are deliberately styled to match coverage_output::plot -- same canvas, fonts, tics, grey
// shading and key placement -- because on a CN evidence page this plot is stacked directly beneath
// that one over the same interval. The two are meant to be read against each other, so the axes
// have to land in the same place.
// ---------------------------------------------------------------------------------------------

// The ori-ter ramp at one position, in GENOMIC coordinates.
//
// The ramp is two straight lines: one from the terminus up to the origin, and one from the origin
// back down to the terminus the other way around the chromosome. Both are evaluated here from the
// two endpoints CNery reports, rather than read out of its per-window otr_gc_corr_fact column.
//
// That column cannot be plotted against position as a line. CNery builds the ramp linear in WINDOW
// INDEX, and it drops every window overlapping repeat coverage -- on REL606 about 1,400 windows are
// missing, with holes up to 5.9 kb. Linear in index is therefore NOT linear in position: the measured
// slope wanders (0.0789 to 0.0827 per Mb on one test) and the line visibly bends and breaks at every
// hole, so the two arms never meet. Evaluating the ramp here restores what it actually is -- two
// straight segments that join at the origin, at the terminus, and across the end of the sequence.
double CNEvidence::otr_ramp_at(const cnery_otr& otr, int32_t position, int32_t seq_length)
{
  int32_t ori = otr.origin, ter = otr.terminus;
  double  yori = otr.origin_cov, yter = otr.terminus_cov;

  // Arm 1 runs directly between the two, arm 2 the other way round, through the sequence end.
  int32_t lo = min(ori, ter), hi = max(ori, ter);
  double  ylo = (lo == ori) ? yori : yter;
  double  yhi = (hi == ori) ? yori : yter;

  if ((position >= lo) && (position <= hi)) {
    if (hi == lo) return ylo;
    return ylo + (yhi - ylo) * static_cast<double>(position - lo) / static_cast<double>(hi - lo);
  }

  // Outside [lo, hi]: on the wrapping arm, measured from hi forward through the end and back to lo.
  int32_t span = seq_length - (hi - lo);
  if (span <= 0) return ylo;
  int32_t along = (position > hi) ? (position - hi) : (seq_length - hi + position);
  return yhi + (ylo - yhi) * static_cast<double>(along) / static_cast<double>(span);
}

// Reduce a window list to at most max_points bins, averaging the coverage over each.
//
// Coverage only. A bin spans many windows and so has no single HMM state, and the copy-number line
// is no longer drawn from windows at all -- it comes from CNery's merged segments, at full precision
// whether or not the coverage under it was binned. This used to carry the copy number FURTHEST FROM
// 1 through each bin so that narrow events stayed visible on a whole-genome overview, which cost
// every event up to a bin of false width on each side; a segment drawn at its true extent is still
// visible, as the two vertical transitions of a spike.
vector<CNEvidence::cnery_window> CNEvidence::bin_cnery_windows(const vector<cnery_window>& in,
                                                               size_t max_points)
{
  if ((max_points == 0) || (in.size() <= max_points)) return in;

  vector<cnery_window> out;
  out.reserve(max_points);

  for (size_t bin = 0; bin < max_points; bin++) {
    size_t lo = (in.size() * bin) / max_points;
    size_t hi = (in.size() * (bin + 1)) / max_points;   // exclusive
    if (hi <= lo) continue;

    cnery_window w = in[lo];
    w.end = in[hi - 1].end;
    w.length = w.end - w.start;

    double raw_sum = 0, corrected_sum = 0, otr_fit_sum = 0;
    double gc_percent_sum = 0, gc_corrected_sum = 0, gc_fit_sum = 0;
    for (size_t i = lo; i < hi; i++) {
      raw_sum += in[i].raw_cov;
      corrected_sum += in[i].corrected_cov;
      otr_fit_sum += in[i].otr_fit_cov;
      gc_percent_sum += in[i].gc_percent;
      gc_corrected_sum += in[i].gc_corrected_cov;
      gc_fit_sum += in[i].gc_corr_fact;
    }
    w.raw_cov = raw_sum / static_cast<double>(hi - lo);
    w.corrected_cov = corrected_sum / static_cast<double>(hi - lo);
    // Not drawn on the overview this binning feeds, but averaged rather than left at the first
    // window's value so a binned list stays a faithful summary of the windows it replaced.
    w.gc_percent = gc_percent_sum / static_cast<double>(hi - lo);
    w.gc_corrected_cov = gc_corrected_sum / static_cast<double>(hi - lo);
    w.gc_corr_fact = gc_fit_sum / static_cast<double>(hi - lo);
    // The ramp is piecewise linear, so a bin mean is its value at the bin centre -- averaging does
    // not distort it the way it smooths the coverage traces.
    w.otr_fit_cov = otr_fit_sum / static_cast<double>(hi - lo);
    // The struct's "no call" sentinel: a bin spans many windows and has no one HMM state, and
    // nothing reads a binned copy number now.
    w.copy_number = -1;

    out.push_back(w);
  }

  return out;
}

void CNEvidence::render_cn_plot(
                                const string& output_svg,
                                const string& seq_id,
                                const vector<cnery_window>& windows,
                                const vector<cnery_segment>& segments,
                                int32_t plot_start,
                                int32_t plot_end,
                                int32_t shaded_start,
                                int32_t shaded_end,
                                const cnery_otr& otr,
                                int32_t seq_length
                                )
{
  // CNery drops any window overlapping repeat/redundant coverage rather than correcting it, so the
  // window list has real holes in it. Emitting a blank line across each hole stops gnuplot drawing
  // a straight line through a region where nothing was measured.
  string tab_file_name = output_svg + ".tab";
  ofstream tab(tab_file_name.c_str());
  ASSERT(tab.good(), "Could not write file: " + tab_file_name);
  tab << "position\traw\tcorrected" << endl;

  double max_y = 0;
  double otr_fit_sum = 0;
  size_t otr_fit_n = 0;
  size_t n_drawn = 0;
  int32_t previous_end = 0;

  for (size_t i = 0; i < windows.size(); i++) {
    const cnery_window& w = windows[i];
    if ((w.end < plot_start) || (w.start > plot_end)) continue;

    if (n_drawn && (w.start > previous_end)) tab << endl;
    previous_end = w.end;

    // Both corners of the window, at the same height, so the coverage traces draw as the step
    // functions they are. A window's coverage is ONE statistic over the whole window, not a sample
    // at a point: joining consecutive midpoints slopes between two windows where the data says a
    // step, and puts values on the plot at positions nothing was measured at. CNery tiles without
    // overlapping, so w.end is exactly where the next window opens (win_len is win_end - win_st,
    // which makes w.end exclusive) and every join lands on a real window boundary.
    //
    // The copy-number line does NOT come from here -- see the segment file below.
    tab << w.start << "\t" << w.raw_cov << "\t" << w.corrected_cov << endl;
    tab << w.end   << "\t" << w.raw_cov << "\t" << w.corrected_cov << endl;

    if (w.otr_fit_cov > 0) {
      // Only kept to set the axis and to place a flat line when there is no ori-ter bias; the ramp
      // itself is drawn from its endpoints, not from these per-window values.
      otr_fit_sum += w.otr_fit_cov;
      otr_fit_n++;
      if (w.otr_fit_cov > max_y) max_y = w.otr_fit_cov;
    }

    if (w.raw_cov > max_y) max_y = w.raw_cov;
    if (w.corrected_cov > max_y) max_y = w.corrected_cov;
    n_drawn++;
  }
  tab.close();

  if (n_drawn == 0) {
    remove(tab_file_name.c_str());
    return;
  }

  // The copy-number calls, as their own file of explicit corner vertices.
  //
  // Drawn from CNery's merged segments rather than from the per-window HMM column, and with lines
  // rather than with steps. Both of those are deliberate. gnuplot's `steps` holds each y from its
  // own x to the NEXT point's x, so a series plotted at window centres -- which is where a coverage
  // measurement belongs -- puts every transition half a window late; on lambda a region called
  // 21701-27700 drew its edges at 21751 and 27751. `steps` also draws no horizontal at all for the
  // last point of a block, so the level before every dropped window went missing. Emitting both
  // corners of each segment and joining them with straight lines has neither problem, and it is the
  // same reasoning the ori-ter ramp file below is written on.
  //
  // The segments tile the sequence, so segment i's closing vertex and segment i+1's opening vertex
  // share an x: the connecting line is exactly vertical and lands exactly on the boundary CNery
  // called -- the same boundary the CN evidence entry reports and the shading above is drawn at.
  string cn_file_name = output_svg + ".cn.tab";
  bool have_copy_number = false;
  {
    ofstream cn(cn_file_name.c_str());
    ASSERT(cn.good(), "Could not write file: " + cn_file_name);
    cn << "position\tcopy_number" << endl;

    for (size_t i = 0; i < segments.size(); i++) {
      const cnery_segment& s = segments[i];
      if ((s.end < plot_start) || (s.start > plot_end)) continue;

      // Clipped to the plotted range. The closing vertex is one past the segment's last base, which
      // is where the next segment opens.
      int32_t from = max(s.start, plot_start);
      int32_t to   = min(s.end + 1, plot_end);
      if (to <= from) continue;

      cn << from << "\t" << s.copy_number << endl;
      cn << to   << "\t" << s.copy_number << endl;
      have_copy_number = true;
      if (s.copy_number > max_y) max_y = s.copy_number;
    }
    cn.close();
  }

  // The ori-ter ramp, as its own file with no gaps in it. Only the vertices are written -- the ends
  // of the plotted range plus whichever of the origin and terminus fall inside it -- because between
  // those the ramp IS a straight line, so gnuplot drawing straight between them is exact. Keeping it
  // out of the main table also keeps it clear of the blank lines that break the coverage traces at
  // unmeasured regions: those breaks are right for measurements and wrong for a model.
  string fit_file_name = output_svg + ".fit.tab";
  bool have_otr_fit = false;
  {
    vector<int32_t> vertices;
    vertices.push_back(plot_start);
    if ((otr.terminus > plot_start) && (otr.terminus < plot_end)) vertices.push_back(otr.terminus);
    if ((otr.origin > plot_start) && (otr.origin < plot_end)) vertices.push_back(otr.origin);
    vertices.push_back(plot_end);
    sort(vertices.begin(), vertices.end());
    vertices.erase(unique(vertices.begin(), vertices.end()), vertices.end());

    ofstream fit(fit_file_name.c_str());
    ASSERT(fit.good(), "Could not write file: " + fit_file_name);
    fit << "position\totr_fit" << endl;
    if (otr.detected && (otr.origin_cov > 0) && (otr.terminus_cov > 0) && (seq_length > 0)) {
      for (size_t k = 0; k < vertices.size(); k++) {
        double y = otr_ramp_at(otr, vertices[k], seq_length);
        fit << vertices[k] << "\t" << y << endl;
        if (y > max_y) max_y = y;
      }
      have_otr_fit = true;
    } else if (otr_fit_n > 0) {
      // No ori-ter bias found: CNery divides by a constant, so the honest picture is a flat line.
      double flat = otr_fit_sum / static_cast<double>(otr_fit_n);
      fit << plot_start << "\t" << flat << endl;
      fit << plot_end << "\t" << flat << endl;
      have_otr_fit = true;
    }
    fit.close();
  }

  // Floor the axis at 2.5 so that a region called at copy number 1 still shows the baseline line
  // with room around it, instead of gnuplot zooming into the noise about 1.0.
  max_y = max(2.5, max_y * 1.1);

  ostringstream s;
  s << "set datafile columnheaders" << endl;
  s << "set datafile missing 'NaN'" << endl;
  // noenhanced: reference sequence ids routinely contain underscores, which gnuplot's enhanced text
  // mode would otherwise render as subscripts ("NC_001416.1" -> "NC(0)01416.1").
  s << "set terminal svg size 2200,1200 font ',28' noenhanced" << endl;
  s << "set output " << double_quote(output_svg) << endl;
  s << "set tics out" << endl;
  s << "set border lw 2" << endl;
  s << "set xlabel 'Coordinate in " << seq_id << "'" << endl;
  s << "set ylabel 'Relative Coverage / Copy Number'" << endl;
  s << "set format x '%.0f'" << endl;
  s << "set xrange [" << plot_start << ":" << plot_end << "]" << endl;
  s << "set yrange [0:" << to_string<double>(max_y) << "]" << endl;

  // Grey out whatever lies outside the region this plot is about, matching coverage_output::plot.
  int obj_id = 1;
  if (plot_start < shaded_start) {
    s << "set object " << (obj_id++) << " rect from " << plot_start << ", graph 0 to " << shaded_start
      << ", graph 1 fc rgb 'grey85' fillstyle solid 1.0 noborder behind" << endl;
  }
  if (shaded_end < plot_end) {
    s << "set object " << (obj_id++) << " rect from " << shaded_end << ", graph 0 to " << plot_end
      << ", graph 1 fc rgb 'grey85' fillstyle solid 1.0 noborder behind" << endl;
  }

  // The origin and terminus of replication CNery fit. These are what the "corrected" series is
  // corrected FOR -- the ramp between them is the bias being divided out -- so seeing where they
  // landed is how you judge whether that correction was reasonable. A terminus sitting exactly on a
  // deleted region, for instance, means the fit was pulled by the deletion rather than by
  // replication timing.
  if (otr.detected) {
    const int32_t markers[2] = { otr.origin, otr.terminus };
    const char* names[2] = { "ori", "ter" };
    const char* colors[2] = { "'dark-green'", "'dark-violet'" };
    for (int m = 0; m < 2; m++) {
      if ((markers[m] < plot_start) || (markers[m] > plot_end)) continue;
      s << "set arrow from " << markers[m] << ", graph 0 to " << markers[m]
        << ", graph 1 nohead lc rgb " << colors[m] << " lw 3 dt 4 front" << endl;
      s << "set label " << double_quote(names[m]) << " at " << markers[m] << ", graph 0.97"
        << " left offset 0.5,0 tc rgb " << colors[m] << " font ',24' front" << endl;
    }
    // Mark where the two straight segments of the ramp meet, i.e. the fitted coverage at each end.
    // Guarded on being real numbers: CNery writes NaN -> null for these when it found no bias.
    if ((otr.origin_cov > 0) && (otr.terminus_cov > 0)) {
      const double marker_y[2] = { otr.origin_cov, otr.terminus_cov };
      for (int m = 0; m < 2; m++) {
        if ((markers[m] < plot_start) || (markers[m] > plot_end)) continue;
        s << "set label \"\" at " << markers[m] << "," << to_string<double>(marker_y[m])
          << " point pointtype 7 pointsize 2 lc rgb " << colors[m] << " front" << endl;
      }
    }

    // State the magnitude numerically as well as drawing it: a ratio near 1 means there was hardly
    // any bias to remove, so a dramatic-looking correction would deserve suspicion.
    ostringstream ratio_label;
    ratio_label << fixed << setprecision(2)
                << "ori-ter bias: " << otr.origin_cov << " at ori to " << otr.terminus_cov
                << " at ter (" << otr.ratio << "x)";
    s << "set label " << double_quote(ratio_label.str()) << " at graph 0.01, graph 0.04"
      << " left tc rgb 'gray30' font ',20' front" << endl;
  }

  s << "set bmargin 6" << endl;
  s << "set key below horizontal Left reverse font ',20' width 4 samplen 1.5" << endl;

  string quoted_tab = double_quote(tab_file_name);
  vector<string> clauses;
  // The single-copy baseline, so a step away from it is visible without reading the axis.
  clauses.push_back("1 with lines lc rgb 'dark-grey' lw 3 dt 2 title 'single copy'");
  clauses.push_back(quoted_tab + " using \"position\":\"raw\" with lines lc rgb 'grey60' lw 3 title 'uncorrected'");
  clauses.push_back(quoted_tab + " using \"position\":\"corrected\" with lines lc rgb 'blue' lw 4 title 'GC + ori-ter corrected'");
  if (have_otr_fit) {
    // The straight-line ramp itself, from its own gap-free file. Its distance from the single-copy
    // line at 1.0 IS the correction applied at each position: corrected = GC-corrected / this.
    // Drawing both together is what makes the size and direction of the correction readable rather
    // than inferred.
    clauses.push_back(double_quote(fit_file_name) + " using \"position\":\"otr_fit\" with lines lc rgb 'black' lw 3 title 'ori-ter bias fit (divided out)'");
  }
  if (have_copy_number) {
    // Piecewise-constant by nature, and drawn as the exact corners of each called segment -- see
    // where the file is written for why this is not `with steps`.
    clauses.push_back(double_quote(cn_file_name) + " using \"position\":\"copy_number\" with lines lc rgb 'red' lw 4 title 'copy number (HMM)'");
  }
  s << "plot " << join(clauses, string(", \\\n     ")) << endl;

  string script_name = output_svg + ".gp";
  string log_name    = output_svg + ".gp.log";
  run_gnuplot_script(s.str(), script_name, log_name);

  make_svg_responsive(output_svg);

  remove(tab_file_name.c_str());
  remove(fit_file_name.c_str());
  remove(cn_file_name.c_str());
  remove(log_name.c_str());
}

void CNEvidence::render_gc_bias_plot(
                                     const string& output_svg,
                                     const string& seq_id,
                                     const vector<cnery_window>& windows
                                     )
{
  // The correction curve doubles as the check for whether there is anything to draw: CNery only
  // emits gc_corr_fact when it ran a GC correction at all.
  vector<pair<double, double> > curve;   // (GC%, correction factor)
  for (size_t i = 0; i < windows.size(); i++) {
    if ((windows[i].gc_corr_fact <= 0.0) || (windows[i].gc_percent <= 0.0)) continue;
    curve.push_back(make_pair(windows[i].gc_percent * 100.0, windows[i].gc_corr_fact));
  }
  if (curve.empty()) return;

  // gc_corr_fact is a function of GC alone, so the tens of thousands of windows collapse to one
  // point per distinct GC value. Sorting by GC is what lets gnuplot join them into the curve.
  sort(curve.begin(), curve.end());
  curve.erase(unique(curve.begin(), curve.end()), curve.end());

  // A whole bacterial chromosome is ~46,000 windows, and two scatter series that size make an SVG
  // tens of megabytes. Striding preserves the shape of the cloud; the curve above is left intact
  // because it is small already and is the part that has to be read precisely.
  const size_t kMaximumScatterPoints = 5000;
  size_t stride = 1;
  {
    size_t n_scatter = 0;
    for (size_t i = 0; i < windows.size(); i++) {
      if ((windows[i].gc_percent > 0.0) && (windows[i].raw_cov > 0.0)) n_scatter++;
    }
    if (n_scatter > kMaximumScatterPoints) stride = (n_scatter / kMaximumScatterPoints) + 1;
  }

  string tab_file_name = output_svg + ".tab";
  ofstream tab(tab_file_name.c_str());
  ASSERT(tab.good(), "Could not write file: " + tab_file_name);
  tab << "gc\traw\tgc_corrected" << endl;

  double max_y = 0;
  size_t kept = 0, n_drawn = 0;
  for (size_t i = 0; i < windows.size(); i++) {
    const cnery_window& w = windows[i];
    if ((w.gc_percent <= 0.0) || (w.raw_cov <= 0.0)) continue;
    if ((kept++ % stride) != 0) continue;

    tab << (w.gc_percent * 100.0) << "\t" << w.raw_cov << "\t";
    // Absent under --bias none: keep the window's raw value on the plot rather than dropping the row.
    if (w.gc_corrected_cov > 0.0) {
      tab << w.gc_corrected_cov;
      if (w.gc_corrected_cov > max_y) max_y = w.gc_corrected_cov;
    } else {
      tab << "NaN";
    }
    tab << endl;

    if (w.raw_cov > max_y) max_y = w.raw_cov;
    n_drawn++;
  }
  tab.close();

  if (n_drawn == 0) {
    remove(tab_file_name.c_str());
    return;
  }

  string curve_file_name = output_svg + ".fit.tab";
  {
    ofstream fit(curve_file_name.c_str());
    ASSERT(fit.good(), "Could not write file: " + curve_file_name);
    fit << "gc\tfactor" << endl;
    for (size_t i = 0; i < curve.size(); i++) {
      fit << curve[i].first << "\t" << curve[i].second << endl;
      if (curve[i].second > max_y) max_y = curve[i].second;
    }
  }

  // Clip the axis rather than fitting it to the tail: a handful of windows sit at many times the
  // median (rDNA, prophage), and letting them set the range flattens the bulk of the cloud into a
  // band too thin to read. 3.0 keeps the single-copy line at 1.0 comfortably mid-plot.
  max_y = min(3.0, max(2.0, max_y * 1.1));

  ostringstream s;
  s << "set datafile columnheaders" << endl;
  s << "set datafile missing 'NaN'" << endl;
  // noenhanced, as in render_cn_plot: reference sequence ids routinely contain underscores.
  s << "set terminal svg size 2200,1200 font ',28' noenhanced" << endl;
  s << "set output " << double_quote(output_svg) << endl;
  s << "set tics out" << endl;
  s << "set border lw 2" << endl;
  s << "set xlabel 'GC content of " << seq_id << " window (%)'" << endl;
  s << "set ylabel 'Relative Coverage'" << endl;
  s << "set yrange [0:" << to_string<double>(max_y) << "]" << endl;
  s << "set bmargin 6" << endl;
  s << "set key below horizontal Left reverse font ',20' width 4 samplen 1.5" << endl;

  string quoted_tab = double_quote(tab_file_name);
  vector<string> clauses;
  clauses.push_back("1 with lines lc rgb 'dark-grey' lw 3 dt 2 title 'single copy'");
  clauses.push_back(quoted_tab + " using \"gc\":\"raw\" with points pt 7 ps 0.4 lc rgb 'grey60' title 'uncorrected'");
  clauses.push_back(quoted_tab + " using \"gc\":\"gc_corrected\" with points pt 7 ps 0.4 lc rgb 'blue' title 'GC corrected'");
  // The curve that was divided out. How far it departs from 1.0 IS the size of the correction, and
  // the blue cloud should sit flat about 1.0 exactly where the grey cloud follows this line.
  clauses.push_back(double_quote(curve_file_name) + " using \"gc\":\"factor\" with lines lc rgb 'black' lw 4 title 'GC bias fit (divided out)'");
  s << "plot " << join(clauses, string(", \\\n     ")) << endl;

  string script_name = output_svg + ".gp";
  string log_name    = output_svg + ".gp.log";
  run_gnuplot_script(s.str(), script_name, log_name);

  make_svg_responsive(output_svg);

  remove(tab_file_name.c_str());
  remove(curve_file_name.c_str());
  remove(log_name.c_str());
}

void CNEvidence::draw_evidence_plots(
                                     const Settings& settings,
                                     cReferenceSequences& ref_seq_info,
                                     cGenomeDiff& gd
                                     )
{
  // draw_coverage() creates this too, but it is called for its own reasons and the two are
  // independent; create_path is a no-op when the directory is already there.
  create_path(settings.evidence_path);

  diff_entry_list_t cn_items = gd.show_list(make_vector<gd_entry_type>(CN));

  // Same set predict() ingested: a junction-only reference was never analyzed, so skipping it here
  // is not a missing plot, and warning about one would be noise on every run that uses -s.
  const set<string> analyzed_seq_ids = settings.call_mutations_seq_id_set();

  for (cReferenceSequences::iterator it = ref_seq_info.begin(); it != ref_seq_info.end(); ++it) {

    cAnnotatedSequence& seq = *it;
    if (analyzed_seq_ids.count(seq.m_seq_id) == 0) continue;

    string cnv_file_name = settings.file_name(settings.cnery_cnv_csv_file_name, "@", seq.m_seq_id);

    vector<cnery_window> windows;
    if (!read_cnery_windows(cnv_file_name, windows) || windows.empty()) {
      // Expected whenever Output is re-run on its own: 09_copy_number_variation/ is deleted once
      // the pipeline completes, so the corrected coverage is simply no longer on disk.
      WARN("Skipping copy number coverage plots for " + seq.m_seq_id +
           ": could not read CNery output file " + cnv_file_name);
      continue;
    }

    cnery_otr otr;
    read_cnery_otr(settings.file_name(settings.cnery_otr_results_file_name, "@", seq.m_seq_id), otr);

    // The copy-number calls as CNery merged them. Drawn instead of the per-window HMM column, whose
    // boundaries are only known to a window and which goes missing wherever CNery dropped one.
    // Absent is survivable -- the plot then simply carries no copy-number line.
    vector<cnery_segment> segments;
    string break_pts_file_name = settings.file_name(settings.cnery_break_points_file_name, "@", seq.m_seq_id);
    if (!read_cnery_segments(break_pts_file_name, segments, static_cast<int32_t>(seq.m_length))) {
      WARN("Omitting the copy number line from the plots for " + seq.m_seq_id +
           ": could not read CNery output file " + break_pts_file_name);
    }

    // Whole-reference overview. Drawn even under --brief-html-mode, matching how draw_coverage
    // always produces its own per-reference overviews.
    {
      string overview_svg = settings.file_name(settings.cn_overview_coverage_plot_file_name, "@", seq.m_seq_id);
      // Only the coverage traces are binned; the copy-number line is drawn from the segments at full
      // precision, so a narrow event stays at its true width instead of being painted a bin wide.
      vector<cnery_window> overview_windows = bin_cnery_windows(windows, 4000);
      render_cn_plot(overview_svg, seq.m_seq_id, overview_windows, segments,
                     1, static_cast<int32_t>(seq.m_length),
                     1, static_cast<int32_t>(seq.m_length), otr,
                     static_cast<int32_t>(seq.m_length));
    }

    // The GC correction the coverage above was corrected BY. Wants the unbinned windows: binning
    // averages together windows of unrelated GC content, which is exactly the axis this plots.
    render_gc_bias_plot(settings.file_name(settings.gc_bias_plot_file_name, "@", seq.m_seq_id),
                        seq.m_seq_id, windows);

    if (settings.no_evidence_html) continue;

    for (diff_entry_list_t::iterator item_it = cn_items.begin(); item_it != cn_items.end(); item_it++) {
      cDiffEntry& item = **item_it;
      if (item[SEQ_ID] != seq.m_seq_id) continue;

      int32_t start = from_string<int32_t>(item[START]);
      int32_t end = from_string<int32_t>(item[END]);

      // Same flanking rule as the raw coverage plot in draw_coverage(), so the two plots stacked
      // on the evidence page cover exactly the same interval.
      int32_t flanking = max(static_cast<int32_t>(100), (end - start + 1) / 10);
      int32_t plot_start = max(static_cast<int32_t>(1), start - flanking);
      int32_t plot_end = min(static_cast<int32_t>(seq.m_length), end + flanking);

      string svg = settings.evidence_path + "/CN_" + item._id + ".svg";
      render_cn_plot(svg, seq.m_seq_id, windows, segments, plot_start, plot_end, start, end, otr,
                     static_cast<int32_t>(seq.m_length));

      if (file_exists(svg.c_str()))
        item["_cn_corrected_plot_file_name"] = Settings::relative_path(svg, settings.evidence_path);
    }
  }
}

} // namespace breseq
