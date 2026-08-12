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

  run_cnery(settings, summary, settings.cnery_output_path);

  for (cReferenceSequences::iterator it = ref_seq_info.begin(); it != ref_seq_info.end(); ++it) {

    cAnnotatedSequence& seq = *it;

    string cnv_file_name = settings.file_name(settings.cnery_cnv_csv_file_name, "@", seq.m_seq_id);
    string break_pts_file_name = settings.file_name(settings.cnery_break_points_file_name, "@", seq.m_seq_id);
    string gd_file_name = settings.file_name(settings.copy_number_evidence_genome_diff_file_name, "@", seq.m_seq_id);

    ingest_csv_for_seq_id(seq.m_seq_id, cnv_file_name, break_pts_file_name, gd_file_name);
  }
}

void CNEvidence::run_cnery(Settings& settings, Summary& summary, const string& cnery_output_prefix)
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
  
  // Prepend breseq's own directory to PATH so CNery can call "breseq bam2cov"
  string command = "PATH=" + double_quote(Settings::get_bin_path() + ":$PATH") + " ";
  command += double_quote(settings.installed["cnery"]);
  command += " -i " + double_quote(settings.base_output_path.size() ? settings.base_output_path : string("."));
  command += " -ref " + double_quote(settings.reference_fasta_file_name);
  command += " -o " + double_quote(cnery_output_prefix);
  //command += " -f " + to_string<uint32_t>(fragment_length);
  
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
  for (size_t i = 0; i < header_fields.size(); i++) {
    if (header_fields[i] == "win_st") win_st_col = i;
    else if (header_fields[i] == "win_len") win_len_col = i;
    else if (header_fields[i] == "otr_gc_corr_norm_cov") rel_cov_col = i;
    else if (header_fields[i] == "win_end") win_end_col = i;
    else if (header_fields[i] == "norm_raw_cov") raw_cov_col = i;
    else if (header_fields[i] == "prob_copy_number") copy_number_col = i;
    else if (header_fields[i] == "otr_gc_corr_fact") otr_fit_col = i;
  }
  ASSERT((win_st_col != string::npos) && (win_len_col != string::npos) && (rel_cov_col != string::npos),
         "Unexpected column layout in CNery output file: " + cnv_file_name);

  // Highest column index actually read, so a truncated row can be skipped rather than indexed past.
  size_t max_needed_col = 0;
  const size_t used_cols[] = { win_st_col, win_len_col, rel_cov_col, win_end_col, raw_cov_col,
                               copy_number_col, otr_fit_col };
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
    w.copy_number = (copy_number_col != string::npos) ? from_string<int32_t>(fields[copy_number_col]) : -1;
    windows.push_back(w);
  }
  cnv_file.close();

  return true;
}

// CNery's <prefix><seq_id>_otr_results.json records the origin and terminus of replication it fit,
// which is what its ori-ter coverage correction is built from.
//
// Two traps in that file. The keys are named "Origin window"/"Terminus window" but hold 1-based
// GENOMIC COORDINATES, not window indices -- CNery converts them with
// `df["win_st"].iloc[ori_idx]` before writing (core.py, otr_correction). And when it detects no
// ori-ter bias it still writes both keys, filled with the first and last coordinate of the sequence
// and a "Not detected" string for the ratio; those placeholders must not be drawn as if they were a
// real ori/ter, which is what the ratio's type is used to distinguish here.
bool CNEvidence::read_cnery_otr(const string& otr_file_name, cnery_otr& otr)
{
  otr.detected = false;
  otr.origin = 0;
  otr.terminus = 0;
  otr.origin_cov = 0.0;
  otr.terminus_cov = 0.0;
  otr.ratio = 0.0;

  ifstream otr_file(otr_file_name.c_str());
  if (!otr_file.good()) return false;

  json j;
  try {
    otr_file >> j;
  } catch (...) {
    WARN("Could not parse CNery output file: " + otr_file_name);
    return false;
  }

  if (!j.count("Origin window") || !j.count("Terminus window")) return false;

  // "Not detected" (a string) rather than a number means CNery found no bias to correct.
  if (!j.count("Origin-to-Termius/Bias Ratio") || !j["Origin-to-Termius/Bias Ratio"].is_number())
    return false;

  otr.origin = j["Origin window"].get<int32_t>();
  otr.terminus = j["Terminus window"].get<int32_t>();
  otr.ratio = j["Origin-to-Termius/Bias Ratio"].get<double>();

  // The two ends of the fitted ramp. NaN in CNery becomes JSON null, so these are only usable when
  // they really are numbers -- the markers that use them are guarded on that separately.
  if (j.count("Origin coverage (normalized)") && j["Origin coverage (normalized)"].is_number())
    otr.origin_cov = j["Origin coverage (normalized)"].get<double>();
  if (j.count("Terminus coverage (normalized)") && j["Terminus coverage (normalized)"].is_number())
    otr.terminus_cov = j["Terminus coverage (normalized)"].get<double>();

  otr.detected = true;

  return true;
}

// Its sibling <prefix><seq_id>_break_pts.csv already merges contiguous runs
// of the same predicted copy number into ranges for us
// (Startpos,State,Segment_Size), so we don't need to re-implement that
// merging here -- we just convert each merged range into a CN evidence
// entry, using the per-window file only to compute a representative
// relative coverage value to display for that range.
void CNEvidence::ingest_csv_for_seq_id(
                                       const string& seq_id,
                                       const string& cnv_file_name,
                                       const string& break_pts_file_name,
                                       const string& gd_file_name
                                       )
{
  vector<cnery_window> windows;
  ASSERT(read_cnery_windows(cnv_file_name, windows),
         "Could not open CNery output file: " + cnv_file_name);

  uint32_t window_size = windows.size() ? static_cast<uint32_t>(windows[0].length) : 0;

  ifstream break_pts_file(break_pts_file_name.c_str());
  ASSERT(break_pts_file.good(), "Could not open CNery output file: " + break_pts_file_name);
  string line;
  getline(break_pts_file, line); // discard header

  cGenomeDiff gd;

  while (getline(break_pts_file, line)) {
    if (line.size() == 0) continue;
    vector<string> fields = split(line, ",");
    ASSERT(fields.size() == 3, "Unexpected number of columns in CNery output file: " + break_pts_file_name);

    // CNery reports a 1-based start position for each merged segment.
    int32_t start_pos = from_string<int32_t>(fields[0]);
    int32_t copy_number = from_string<int32_t>(fields[1]);
    int32_t segment_size = from_string<int32_t>(fields[2]);
    int32_t end_pos = start_pos + segment_size - 1;

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
    item["relative_coverage"] = to_string<double>((coverage_n > 0) ? (coverage_sum / coverage_n) : 0.0);

    gd.add(item);
  }
  break_pts_file.close();

  gd.write(gd_file_name);
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

// Coverage is averaged over each bin, but the copy number kept is the one FURTHEST FROM 1 -- taking
// a max would hide deletions and taking a mean would smear every narrow event into the baseline,
// and it is exactly the narrow events an overview is scanned for.
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
    int32_t extreme_cn = in[lo].copy_number;
    for (size_t i = lo; i < hi; i++) {
      raw_sum += in[i].raw_cov;
      corrected_sum += in[i].corrected_cov;
      otr_fit_sum += in[i].otr_fit_cov;
      if (labs(in[i].copy_number - 1) > labs(extreme_cn - 1)) extreme_cn = in[i].copy_number;
    }
    w.raw_cov = raw_sum / static_cast<double>(hi - lo);
    w.corrected_cov = corrected_sum / static_cast<double>(hi - lo);
    // The ramp is piecewise linear, so a bin mean is its value at the bin centre -- averaging does
    // not distort it the way it smooths the coverage traces.
    w.otr_fit_cov = otr_fit_sum / static_cast<double>(hi - lo);
    w.copy_number = extreme_cn;

    out.push_back(w);
  }

  return out;
}

void CNEvidence::render_cn_plot(
                                const string& output_svg,
                                const string& seq_id,
                                const vector<cnery_window>& windows,
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
  tab << "position\traw\tcorrected\tcopy_number" << endl;

  double max_y = 0;
  bool have_copy_number = false;
  double otr_fit_sum = 0;
  size_t otr_fit_n = 0;
  size_t n_drawn = 0;
  int32_t previous_end = 0;

  for (size_t i = 0; i < windows.size(); i++) {
    const cnery_window& w = windows[i];
    if ((w.end < plot_start) || (w.start > plot_end)) continue;

    if (n_drawn && (w.start > previous_end)) tab << endl;
    previous_end = w.end;

    // Windows overlap (default 200 bp wide on a 100 bp step), so each one is plotted at its
    // midpoint rather than at either edge.
    int32_t midpoint = w.start + (w.end - w.start) / 2;
    tab << midpoint << "\t" << w.raw_cov << "\t" << w.corrected_cov << "\t";
    if (w.otr_fit_cov > 0) {
      // Only kept to set the axis and to place a flat line when there is no ori-ter bias; the ramp
      // itself is drawn from its endpoints, not from these per-window values.
      otr_fit_sum += w.otr_fit_cov;
      otr_fit_n++;
      if (w.otr_fit_cov > max_y) max_y = w.otr_fit_cov;
    }
    if (w.copy_number >= 0) {
      tab << w.copy_number;
      have_copy_number = true;
      if (w.copy_number > max_y) max_y = w.copy_number;
    } else {
      tab << "NaN";
    }
    tab << endl;

    if (w.raw_cov > max_y) max_y = w.raw_cov;
    if (w.corrected_cov > max_y) max_y = w.corrected_cov;
    n_drawn++;
  }
  tab.close();

  if (n_drawn == 0) {
    remove(tab_file_name.c_str());
    return;
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
    // Piecewise-constant by nature (one HMM state per window), hence steps rather than lines.
    clauses.push_back(quoted_tab + " using \"position\":\"copy_number\" with steps lc rgb 'red' lw 4 title 'copy number (HMM)'");
  }
  s << "plot " << join(clauses, string(", \\\n     ")) << endl;

  string script_name = output_svg + ".gp";
  string log_name    = output_svg + ".gp.log";
  run_gnuplot_script(s.str(), script_name, log_name);

  make_svg_responsive(output_svg);

  remove(tab_file_name.c_str());
  remove(fit_file_name.c_str());
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

  for (cReferenceSequences::iterator it = ref_seq_info.begin(); it != ref_seq_info.end(); ++it) {

    cAnnotatedSequence& seq = *it;
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

    // Whole-reference overview. Drawn even under --brief-html-mode, matching how draw_coverage
    // always produces its own per-reference overviews.
    {
      string overview_svg = settings.file_name(settings.cn_overview_coverage_plot_file_name, "@", seq.m_seq_id);
      vector<cnery_window> overview_windows = bin_cnery_windows(windows, 4000);
      render_cn_plot(overview_svg, seq.m_seq_id, overview_windows,
                     1, static_cast<int32_t>(seq.m_length),
                     1, static_cast<int32_t>(seq.m_length), otr,
                     static_cast<int32_t>(seq.m_length));
    }

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
      render_cn_plot(svg, seq.m_seq_id, windows, plot_start, plot_end, start, end, otr,
                     static_cast<int32_t>(seq.m_length));

      if (file_exists(svg.c_str()))
        item["_cn_corrected_plot_file_name"] = Settings::relative_path(svg, settings.evidence_path);
    }
  }
}

} // namespace breseq
