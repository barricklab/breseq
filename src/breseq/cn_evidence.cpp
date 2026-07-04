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
    settings.installed["cnery"] = SYSTEM_CAPTURE("which CNery", true);
  }
  if (settings.installed["cnery"].size() == 0) {
    ERROR("Could not find 'CNery' command in $PATH.\n"
          "Install it (e.g. 'pip install CNery') to use --cn-evidence.\n"
          "See https://github.com/barricklab/CNery");
  }

  string cnery_output_prefix = settings.copy_number_variation_path + "/cnery_out";
  run_cnery(settings, summary, cnery_output_prefix);

  string cnery_basename = path_to_filename(cnery_output_prefix);
  string csv_path = cnery_output_prefix + "/CNV_csv";

  for (cReferenceSequences::iterator it = ref_seq_info.begin(); it != ref_seq_info.end(); ++it) {

    cAnnotatedSequence& seq = *it;

    string cnv_file_name = csv_path + "/" + cnery_basename + seq.m_seq_id + "_CNV.csv";
    string break_pts_file_name = csv_path + "/" + cnery_basename + seq.m_seq_id + "_break_pts.csv";
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

// CNery's <prefix><seq_id>_CNV.csv has one row per sliding window, with a
// per-window predicted copy number ("prob_copy_number", the final column).
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
  ifstream cnv_file(cnv_file_name.c_str());
  ASSERT(cnv_file.good(), "Could not open CNery output file: " + cnv_file_name);

  string header_line;
  getline(cnv_file, header_line);
  vector<string> header_fields = split(header_line, ",");

  size_t win_st_col = string::npos, win_len_col = string::npos, rel_cov_col = string::npos;
  for (size_t i = 0; i < header_fields.size(); i++) {
    if (header_fields[i] == "win_st") win_st_col = i;
    else if (header_fields[i] == "win_len") win_len_col = i;
    else if (header_fields[i] == "otr_gc_corr_norm_cov") rel_cov_col = i;
  }
  ASSERT((win_st_col != string::npos) && (win_len_col != string::npos) && (rel_cov_col != string::npos),
         "Unexpected column layout in CNery output file: " + cnv_file_name);

  vector<int32_t> window_start;
  vector<double> window_relative_coverage;
  uint32_t window_size = 0;

  string line;
  while (getline(cnv_file, line)) {
    if (line.size() == 0) continue;
    vector<string> fields = split(line, ",");
    window_start.push_back(from_string<int32_t>(fields[win_st_col]));
    window_relative_coverage.push_back(from_string<double>(fields[rel_cov_col]));
    if (window_size == 0) window_size = from_string<uint32_t>(fields[win_len_col]);
  }
  cnv_file.close();

  ifstream break_pts_file(break_pts_file_name.c_str());
  ASSERT(break_pts_file.good(), "Could not open CNery output file: " + break_pts_file_name);
  getline(break_pts_file, header_line); // discard header

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
    for (size_t i = 0; i < window_start.size(); i++) {
      if ((window_start[i] >= start_pos) && (window_start[i] <= end_pos)) {
        coverage_sum += window_relative_coverage[i];
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

} // namespace breseq
