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

#include "libbreseq/fastq.h"
#include "libbreseq/reference_sequence.h"

using namespace std;

namespace breseq {
  
  /*
   normalize_fastq
   
   correct common errors in input fastq and normalize to standard SANGER format

   */

  AnalyzeFastqSummary normalize_fastq(
                                        const string &file_name,
                                        const string &convert_file_name,
                                        const uint32_t file_index,
                                        const int32_t trim_end_on_base_quality,
                                        const bool filter_reads,
                                        uint64_t current_read_file_bases,
                                        const uint64_t read_file_base_limit,
                                        const uint32_t _read_length_min,
                                        const double _max_same_base_fraction,
                                        const double _max_N_fraction,
                                        const uint32_t _long_read_trigger_length,
                                        const uint32_t _long_read_split_length,
                                        const bool _long_read_distribute_remainder,
                                        const uint32_t num_threads
                                        )
{
    cerr << "    Converting/filtering FASTQ file..." << endl;

    // Set up maps between formats
    map<string,uint8_t> format_to_chr_offset;
    format_to_chr_offset["SANGER"] = 33;
    format_to_chr_offset["SOLEXA"] = 64;
    format_to_chr_offset["ILLUMINA_1.3+"] = 64;

    // Honor that zero means no splitting
    uint32_t long_read_trigger_length = _long_read_trigger_length;
    if (long_read_trigger_length == 0) {
      long_read_trigger_length = numeric_limits<uint32_t>::max();
    }
    uint32_t long_read_split_length = _long_read_split_length;
    if (long_read_split_length == 0) {
      long_read_split_length = numeric_limits<uint32_t>::max();
    }
    bool long_read_distribute_remainder = _long_read_distribute_remainder;

    const uint64_t initial_current_read_file_bases = current_read_file_bases;

    // Same thresholds as cFastqQualityConverter::predict_fastq_file_format(), but
    // operating on an already-known min raw quality score.
    auto predict_format_from_min_quality_score = [&](uint8_t min_q) -> string {
      string format = "SANGER";
      if (min_q >= format_to_chr_offset["SOLEXA"] - 5) {
        format = "SOLEXA";
      }
      if (min_q >= format_to_chr_offset["ILLUMINA_1.3+"]) {
        format = "ILLUMINA_1.3+";
      }
      return format;
    };

    // Per-read processing: quality convert, trim, split, filter, write, and update stats.
    struct AttemptResult {
      uint64_t num_original_reads = 0;
      uint64_t num_original_bases = 0;
      uint8_t  overall_min_quality_score = 255;
      uint8_t  overall_max_quality_score = 0;
      uint64_t num_reads = 0;
      uint64_t num_bases = 0;
      uint64_t num_filtered_too_short_reads = 0;
      uint64_t num_filtered_too_short_bases = 0;
      uint64_t num_filtered_same_base_reads = 0;
      uint64_t num_filtered_same_base_bases = 0;
      uint64_t num_filtered_too_many_N_reads = 0;
      uint64_t num_filtered_too_many_N_bases = 0;
      uint32_t read_length_min = numeric_limits<uint32_t>::max();
      uint32_t read_length_max = 0;
      bool     file_has_split_reads = false;
      bool     reached_eof = true;
    };

    auto run_pass = [&](const string& pass_quality_format,
                         std::function<bool(cFastqSequence&)> get_next_sequence,
                         uint64_t starting_read_file_bases) -> AttemptResult
    {
      AttemptResult r;

      cFastqQualityConverter fqc(pass_quality_format, "SANGER");
      cFastqFile output_fastq_file(convert_file_name.c_str(), fstream::out, num_threads);

      uint64_t local_current_read_file_bases = starting_read_file_bases;
      uint32_t on_read = 1;
      cFastqSequence original_sequence;

      while (get_next_sequence(original_sequence)) {

        r.num_original_reads++;
        r.num_original_bases += original_sequence.length();

        // track raw (pre-conversion) min/max quality score
        for (uint32_t i=0; i<original_sequence.m_qualities.size(); i++) {
          uint8_t this_score = static_cast<uint8_t>(original_sequence.m_qualities[i]);
          r.overall_min_quality_score = min(r.overall_min_quality_score, this_score);
          r.overall_max_quality_score = max(r.overall_max_quality_score, this_score);
        }

        // truncate second name
        original_sequence.m_name_plus = "";

        // fastq quality convert
        fqc.convert_sequence(original_sequence);

        // trim bad quality scores from the end
        if (trim_end_on_base_quality) {
          fastq_sequence_trim_end_on_base_quality(original_sequence, trim_end_on_base_quality);
        }

        // New loop that enables us to split long reads

        // Uniformly name, to prevent problems drawing alignments
        // and allows us to know the input order when merge sorting later
        original_sequence.m_name = to_string(file_index) + ":" + to_string(on_read++);

        // Decide on a per-read basis whether this read is subject to splitting
        uint32_t effective_split_length = (original_sequence.m_sequence.length() < long_read_trigger_length)
                                            ? numeric_limits<uint32_t>::max()
                                            : long_read_split_length;
        if (effective_split_length != numeric_limits<uint32_t>::max()) {
          r.file_has_split_reads = true;
        }

        // Always have at least one piece and only add new pieces if it is over the length (hence minus 1)
        uint32_t num_split_read_pieces = 1 + (original_sequence.m_sequence.length()-1) / effective_split_length;

        size_t chunk_start_0, chunk_end_0;
        double chunk_size = effective_split_length;
        bool read_was_split = num_split_read_pieces > 1;

        if (read_was_split) {
          if (!long_read_distribute_remainder) {
            // One fewer pieces because we ignore the last if not distributing the remainder
            num_split_read_pieces--;
          } else {
            // Update the chunk size (can be fractional) if we are distributing the remainder
            chunk_size = static_cast<double>(original_sequence.m_sequence.length()) / static_cast<double>(num_split_read_pieces);
          }
        }

        for (uint32_t i=0; i<num_split_read_pieces; i++) {

          cFastqSequence on_sequence;
          chunk_start_0 = ceil(i * chunk_size);
          chunk_end_0 = ceil((i+1) * chunk_size) - 1;

          // Suffix name with chunk number
          on_sequence.m_name = original_sequence.m_name;
          if (read_was_split) on_sequence.m_name += "S" + to_string(i+1);

          // Find correct chunk of bases and quals
          on_sequence.m_sequence = original_sequence.m_sequence.substr(chunk_start_0, chunk_end_0 - chunk_start_0 + 1);
          on_sequence.m_qualities = original_sequence.m_qualities.substr(chunk_start_0, chunk_end_0 - chunk_start_0 + 1);

          if ( filter_reads ) {

            // Discard sequences that are too short
            if ( _read_length_min && (on_sequence.length() < _read_length_min) ) {
              r.num_filtered_too_short_reads++;
              r.num_filtered_too_short_bases += on_sequence.length();
              continue;
            }

            // If doing tests that require, copy over info that is normally calculated on reading
            if (_max_N_fraction || _max_same_base_fraction) {

              if (read_was_split) {
                // Have to recalculate these since they are only done on reading the FASTQ right now
                for(size_t j=0; j<on_sequence.length(); j++) {
                  on_sequence.m_base_counts[basechar2index(on_sequence.m_sequence[j])]++;
                }
              } else {
                // Copy over if we didn't split the read
                for (uint8_t b=0; b<base_list_including_N_size; b++) {
                  on_sequence.m_base_counts[b] = original_sequence.m_base_counts[b];
                }
              }
            }

            // Discard sequences that are 50% or more N.
            if ( _max_N_fraction ) {
              if (_max_N_fraction * static_cast<double>(on_sequence.length()) <= on_sequence.m_base_counts[base_list_N_index]) {
                r.num_filtered_too_many_N_reads++;
                r.num_filtered_too_many_N_bases += on_sequence.length();
                continue;
              }
            }

            // Ignore heavily homopolymer reads, as these are a common type of machine error
            // Discard sequences that are 90% or more of a single base or N.
            bool same_base_filtered = false;
            if (_max_same_base_fraction) {
              for (uint8_t b=0; b<base_list_including_N_size; b++) {
                if ((_max_same_base_fraction * static_cast<double>(on_sequence.length())) <=
                    static_cast<double>(on_sequence.m_base_counts[b] + on_sequence.m_base_counts[base_list_N_index]) ) {
                  same_base_filtered = true;
                  break;
                }
              }
            }

            if (same_base_filtered)  {
              r.num_filtered_same_base_reads++;
              r.num_filtered_same_base_bases += on_sequence.length();
              continue;
            }

          } // end filter read block
          r.num_reads++;
          r.num_bases += on_sequence.m_sequence.length();

          if (read_file_base_limit) {
            local_current_read_file_bases += on_sequence.m_sequence.length();
          }

          r.read_length_min = min<size_t>(on_sequence.length(), r.read_length_min);
          r.read_length_max = max<size_t>(on_sequence.length(), r.read_length_max);

          output_fastq_file.write_sequence(on_sequence);
        }

        // check to see if we've reached the limit
        // outside of loop b/c we always process complete input reads
        if (read_file_base_limit) {
          if (local_current_read_file_bases > read_file_base_limit) {
            r.reached_eof = false;
            break;
          }
        }
      }

      return r;
    };

    // ---- Sniff up to FASTQ_FORMAT_SNIFF_READS reads to predict the quality format ----
    const size_t FASTQ_FORMAT_SNIFF_READS = 10000;
    const int32_t prelim_offset = 64;

    cFastqFile input_fastq_file(file_name.c_str(), fstream::in, num_threads);
    cFastqQualityConverter prelim_fqc("ILLUMINA_1.3+", "SANGER");

    vector<cFastqSequence> buffered_reads;
    buffered_reads.reserve(FASTQ_FORMAT_SNIFF_READS);

    uint8_t sniffed_min_quality_score = 255;
    uint8_t sniffed_max_quality_score = 0;

    {
      cFastqSequence seq;
      while ( (buffered_reads.size() < FASTQ_FORMAT_SNIFF_READS) && input_fastq_file.read_sequence(seq, prelim_fqc) ) {
        for (uint32_t i=0; i<seq.m_qualities.size(); i++) {
          uint8_t this_score = static_cast<uint8_t>(seq.m_qualities[i]);
          sniffed_min_quality_score = min(sniffed_min_quality_score, this_score);
          sniffed_max_quality_score = max(sniffed_max_quality_score, this_score);
        }
        buffered_reads.push_back(seq);
        seq = cFastqSequence();
      }
    }

    string quality_format = predict_format_from_min_quality_score(sniffed_min_quality_score);

    // Adjustment for buffered reads that were read in numerical-quality format using
    // the preliminary ILLUMINA_1.3+ offset (64), but should have used the chosen format's offset.
    int32_t numerical_quality_offset_adjustment = prelim_offset - static_cast<int32_t>(format_to_chr_offset[quality_format]);

    size_t buffered_idx = 0;
    auto get_buffered_then_stream = [&](cFastqSequence& seq) -> bool {
      if (buffered_idx < buffered_reads.size()) {
        seq = std::move(buffered_reads[buffered_idx++]);
        if (seq.m_numerical_qualities && numerical_quality_offset_adjustment != 0) {
          for (uint32_t i=0; i<seq.m_qualities.size(); i++) {
            seq.m_qualities[i] = static_cast<char>(static_cast<int32_t>(static_cast<uint8_t>(seq.m_qualities[i])) - numerical_quality_offset_adjustment);
          }
        }
        return true;
      }
      cFastqQualityConverter fqc_for_offset(quality_format, "SANGER");
      return input_fastq_file.read_sequence(seq, fqc_for_offset);
    };

    AttemptResult result = run_pass(quality_format, get_buffered_then_stream, current_read_file_bases);

    // ---- Check whether the format predicted from the sniff matches the format
    // derived from the actual min/max quality scores seen over the whole file ----
    string corrected_quality_format = predict_format_from_min_quality_score(result.overall_min_quality_score);

    cFastqFile* tail_count_stream = &input_fastq_file;
    unique_ptr<cFastqFile> fresh_input_fastq_file;

    if (corrected_quality_format != quality_format) {

      cerr << "    Warning: Quality score format predicted from the first " << FASTQ_FORMAT_SNIFF_READS << " reads (" << quality_format
           << ") did not match the format determined from the entire file (" << corrected_quality_format << ")." << endl;
      cerr << "    Re-converting FASTQ file using corrected quality score format..." << endl;

      quality_format = corrected_quality_format;

      fresh_input_fastq_file.reset(new cFastqFile(file_name.c_str(), fstream::in, num_threads));
      cFastqQualityConverter fresh_read_fqc(quality_format, "SANGER");

      auto get_fresh = [&](cFastqSequence& seq) -> bool {
        return fresh_input_fastq_file->read_sequence(seq, fresh_read_fqc);
      };

      result = run_pass(quality_format, get_fresh, initial_current_read_file_bases);
      tail_count_stream = fresh_input_fastq_file.get();
    }

    // If we stopped early due to read_file_base_limit, count (but don't convert/write)
    // the remaining reads so the "filtered due to coverage limit" stats are accurate.
    uint64_t num_filtered_coverage_limit_reads = 0;
    uint64_t num_filtered_coverage_limit_bases = 0;

    if (!result.reached_eof) {
      cFastqQualityConverter tail_fqc("ILLUMINA_1.3+", "SANGER");
      cFastqSequence seq;
      while (tail_count_stream->read_sequence(seq, tail_fqc)) {
        num_filtered_coverage_limit_reads++;
        num_filtered_coverage_limit_bases += seq.length();
      }
    }

    uint64_t num_original_reads = result.num_original_reads + num_filtered_coverage_limit_reads;
    uint64_t num_original_bases = result.num_original_bases + num_filtered_coverage_limit_bases;

    // Convert the overall min/max raw quality scores to SANGER for reporting
    cFastqQualityConverter final_fqc(quality_format, "SANGER");
    cFastqSequence min_max_sequence;
    min_max_sequence.m_qualities.append(1, (char)result.overall_min_quality_score);
    min_max_sequence.m_qualities.append(1, (char)result.overall_max_quality_score);
    final_fqc.convert_sequence(min_max_sequence);
    uint8_t min_quality_score = (uint8_t)min_max_sequence.m_qualities[0] - format_to_chr_offset["SANGER"];
    uint8_t max_quality_score = (uint8_t)min_max_sequence.m_qualities[1] - format_to_chr_offset["SANGER"];

    uint32_t width_for_reads = to_string(num_original_reads).size();
    uint32_t width_for_bases = to_string(num_original_bases).size();

    cerr << "    Original base quality format: " << quality_format << " New format: SANGER"<< endl;
    cerr << "    Original reads: " << num_original_reads << " bases: "<< num_original_bases << endl;

    if (num_filtered_coverage_limit_reads) {
      cerr << "    Filtered reads: " << setw(width_for_reads) << num_filtered_coverage_limit_reads;
      cerr << " bases: "<< setw(width_for_bases) << num_filtered_coverage_limit_bases;
      cerr << " (coverage limit option)" << endl;
    }

    // Then report splitting because other filters happened on split reads
    if (result.file_has_split_reads) {
      uint64_t total_split_reads = result.num_reads + result.num_filtered_same_base_reads + result.num_filtered_too_many_N_reads;
      uint64_t total_split_bases = result.num_bases + result.num_filtered_same_base_bases + result.num_filtered_too_many_N_bases;

      cerr << "    >> Long reads split to " << (long_read_distribute_remainder ? "≤": "exactly ");
      cout << long_read_split_length << " bases" << (long_read_distribute_remainder ? "" : " (extra bases discarded)") << endl;
      cerr << "    >> Split reads: " << setw(width_for_reads) << total_split_reads;
      cerr << " bases: "<< setw(width_for_bases) << total_split_bases << endl;
    }

    if (filter_reads) {
      if (result.num_filtered_too_many_N_reads + result.num_filtered_same_base_reads + result.num_filtered_too_short_reads + num_filtered_coverage_limit_reads == 0) {
        cerr << "    Filtered reads: none" << endl;
      } else {

        if (result.num_filtered_too_short_reads) {
          cerr << "    Filtered reads: " << setw(width_for_reads) << result.num_filtered_too_short_reads;
          cerr << " bases: "<< setw(width_for_bases) << result.num_filtered_too_short_bases;
          cerr << " (<" << _read_length_min << " bases long)" << endl;
        }
        if (result.num_filtered_too_many_N_reads) {
          string percentage = formatted_double(100 * _max_N_fraction, 1).to_string();
          cerr << "    Filtered reads: " << setw(width_for_reads) << result.num_filtered_too_many_N_reads;
          cerr << " bases: "<< setw(width_for_bases) << result.num_filtered_too_many_N_bases;
          cerr << " (≥" << percentage << "% bases N)" << endl;
        }
        if (result.num_filtered_same_base_reads) {
          string percentage = formatted_double(100 * _max_same_base_fraction, 1).to_string();
          cerr << "    Filtered reads: " << setw(width_for_reads) << result.num_filtered_same_base_reads;
          cerr << " bases: "<< setw(width_for_bases) << result.num_filtered_same_base_bases;
          cerr << " (≥" << percentage << "% same base)" << endl;
        }
      }
    }
    cerr << "    Analyzed reads: " << setw(width_for_reads) << result.num_reads << " bases: " << setw(width_for_bases) << result.num_bases << endl;

    double read_length_avg = (result.num_reads > 0) ? static_cast<double>(result.num_bases) / static_cast<double>(result.num_reads) : 0;

    AnalyzeFastqSummary retval(
                                 result.read_length_min,
                                 result.read_length_max,
                                 read_length_avg,
                                 num_original_reads,
                                 result.num_filtered_too_short_reads,
                                 result.num_filtered_same_base_reads,
                                 result.num_filtered_too_many_N_reads,
                                 num_filtered_coverage_limit_reads,
                                 result.num_reads,
                                 min_quality_score,
                                 max_quality_score,
                                 num_original_bases,
                                 result.num_bases,
                                 result.file_has_split_reads,
                                 quality_format,
                                 "SANGER",
                                 convert_file_name
                                 );
    return retval;
  }
  
  pair<AnalyzeFastqSummary, AnalyzeFastqSummary> normalize_fastq_paired(
                                        const string &r1_file_name,
                                        const string &r1_convert_file_name,
                                        const string &r2_file_name,
                                        const string &r2_convert_file_name,
                                        const uint32_t r1_file_index,
                                        const uint32_t r2_file_index,
                                        const int32_t trim_end_on_base_quality,
                                        const bool filter_reads,
                                        uint64_t current_read_file_bases,
                                        const uint64_t read_file_base_limit,
                                        const uint32_t _read_length_min,
                                        const double _max_same_base_fraction,
                                        const double _max_N_fraction,
                                        const uint32_t _long_read_trigger_length,
                                        const uint32_t _long_read_split_length,
                                        const bool _long_read_distribute_remainder,
                                        const uint32_t num_threads
                                        )
  {
    (void)_long_read_split_length;
    (void)_long_read_distribute_remainder;

    cerr << "    Converting/filtering paired FASTQ files..." << endl;

    map<string,uint8_t> format_to_chr_offset;
    format_to_chr_offset["SANGER"] = 33;
    format_to_chr_offset["SOLEXA"] = 64;
    format_to_chr_offset["ILLUMINA_1.3+"] = 64;

    const uint64_t initial_current_read_file_bases = current_read_file_bases;

    auto predict_format_from_min_quality_score = [&](uint8_t min_q) -> string {
      string fmt = "SANGER";
      if (min_q >= format_to_chr_offset["SOLEXA"] - 5) fmt = "SOLEXA";
      if (min_q >= format_to_chr_offset["ILLUMINA_1.3+"]) fmt = "ILLUMINA_1.3+";
      return fmt;
    };

    // Per-file stats accumulated inside run_pass_paired
    struct FileResult {
      uint64_t num_original_reads = 0;
      uint64_t num_original_bases = 0;
      uint8_t  overall_min_quality_score = 255;
      uint8_t  overall_max_quality_score = 0;
      uint64_t num_reads = 0;
      uint64_t num_bases = 0;
      uint64_t num_filtered_too_short_reads = 0;
      uint64_t num_filtered_too_short_bases = 0;
      uint64_t num_filtered_same_base_reads = 0;
      uint64_t num_filtered_same_base_bases = 0;
      uint64_t num_filtered_too_many_N_reads = 0;
      uint64_t num_filtered_too_many_N_bases = 0;
      uint32_t read_length_min = numeric_limits<uint32_t>::max();
      uint32_t read_length_max = 0;
    };

    struct PairResult {
      FileResult r1, r2;
      bool reached_eof = true;
    };

    auto run_pass_paired = [&](
        const string& pass_quality_format,
        std::function<bool(cFastqSequence&)> get_next_r1,
        std::function<bool(cFastqSequence&)> get_next_r2,
        uint64_t starting_read_file_bases) -> PairResult
    {
      PairResult pr;

      cFastqQualityConverter fqc(pass_quality_format, "SANGER");
      cFastqFile out_r1(r1_convert_file_name.c_str(), fstream::out, num_threads);
      cFastqFile out_r2(r2_convert_file_name.c_str(), fstream::out, num_threads);

      uint64_t local_current_bases = starting_read_file_bases;
      uint32_t on_read = 1;
      bool warned_long_reads = false;
      uint32_t long_read_trigger = (_long_read_trigger_length == 0) ? numeric_limits<uint32_t>::max() : _long_read_trigger_length;

      cFastqSequence orig_r1, orig_r2;

      while (get_next_r1(orig_r1)) {
        if (!get_next_r2(orig_r2)) {
          cerr << "  Warning: R2 file ended before R1 in paired FASTQ processing. Files may have different read counts." << endl;
          break;
        }

        // Count originals
        pr.r1.num_original_reads++;
        pr.r1.num_original_bases += orig_r1.length();
        pr.r2.num_original_reads++;
        pr.r2.num_original_bases += orig_r2.length();

        // Track raw quality scores (pre-conversion) for format detection
        for (uint32_t i = 0; i < orig_r1.m_qualities.size(); i++) {
          uint8_t q = static_cast<uint8_t>(orig_r1.m_qualities[i]);
          pr.r1.overall_min_quality_score = min(pr.r1.overall_min_quality_score, q);
          pr.r1.overall_max_quality_score = max(pr.r1.overall_max_quality_score, q);
        }
        for (uint32_t i = 0; i < orig_r2.m_qualities.size(); i++) {
          uint8_t q = static_cast<uint8_t>(orig_r2.m_qualities[i]);
          pr.r2.overall_min_quality_score = min(pr.r2.overall_min_quality_score, q);
          pr.r2.overall_max_quality_score = max(pr.r2.overall_max_quality_score, q);
        }

        orig_r1.m_name_plus = "";
        orig_r2.m_name_plus = "";

        fqc.convert_sequence(orig_r1);
        fqc.convert_sequence(orig_r2);

        if (trim_end_on_base_quality) {
          fastq_sequence_trim_end_on_base_quality(orig_r1, trim_end_on_base_quality);
          fastq_sequence_trim_end_on_base_quality(orig_r2, trim_end_on_base_quality);
        }

        // Name reads with their respective file indices so names remain unique
        orig_r1.m_name = to_string(r1_file_index) + ":" + to_string(on_read);
        orig_r2.m_name = to_string(r2_file_index) + ":" + to_string(on_read);
        on_read++;

        if (!warned_long_reads && (orig_r1.length() >= long_read_trigger || orig_r2.length() >= long_read_trigger)) {
          cerr << "    Warning: Long read(s) detected in paired FASTQ files. Long-read splitting is not supported in paired mode; reads will be processed without splitting." << endl;
          warned_long_reads = true;
        }

        // Filter: check R1 first, then R2. If either fails, skip both and
        // increment the triggered filter count in both R1 and R2 stats.
        if (filter_reads) {

          if (_read_length_min && (orig_r1.length() < _read_length_min)) {
            pr.r1.num_filtered_too_short_reads++;
            pr.r1.num_filtered_too_short_bases += orig_r1.length();
            pr.r2.num_filtered_too_short_reads++;
            pr.r2.num_filtered_too_short_bases += orig_r2.length();
            continue;
          }
          if (_read_length_min && (orig_r2.length() < _read_length_min)) {
            pr.r1.num_filtered_too_short_reads++;
            pr.r1.num_filtered_too_short_bases += orig_r1.length();
            pr.r2.num_filtered_too_short_reads++;
            pr.r2.num_filtered_too_short_bases += orig_r2.length();
            continue;
          }

          if (_max_N_fraction) {
            if (_max_N_fraction * static_cast<double>(orig_r1.length()) <= orig_r1.m_base_counts[base_list_N_index]) {
              pr.r1.num_filtered_too_many_N_reads++;
              pr.r1.num_filtered_too_many_N_bases += orig_r1.length();
              pr.r2.num_filtered_too_many_N_reads++;
              pr.r2.num_filtered_too_many_N_bases += orig_r2.length();
              continue;
            }
            if (_max_N_fraction * static_cast<double>(orig_r2.length()) <= orig_r2.m_base_counts[base_list_N_index]) {
              pr.r1.num_filtered_too_many_N_reads++;
              pr.r1.num_filtered_too_many_N_bases += orig_r1.length();
              pr.r2.num_filtered_too_many_N_reads++;
              pr.r2.num_filtered_too_many_N_bases += orig_r2.length();
              continue;
            }
          }

          if (_max_same_base_fraction) {
            bool r1_same_base = false;
            for (uint8_t b = 0; b < base_list_including_N_size; b++) {
              if (_max_same_base_fraction * static_cast<double>(orig_r1.length()) <=
                  static_cast<double>(orig_r1.m_base_counts[b] + orig_r1.m_base_counts[base_list_N_index])) {
                r1_same_base = true;
                break;
              }
            }
            if (r1_same_base) {
              pr.r1.num_filtered_same_base_reads++;
              pr.r1.num_filtered_same_base_bases += orig_r1.length();
              pr.r2.num_filtered_same_base_reads++;
              pr.r2.num_filtered_same_base_bases += orig_r2.length();
              continue;
            }

            bool r2_same_base = false;
            for (uint8_t b = 0; b < base_list_including_N_size; b++) {
              if (_max_same_base_fraction * static_cast<double>(orig_r2.length()) <=
                  static_cast<double>(orig_r2.m_base_counts[b] + orig_r2.m_base_counts[base_list_N_index])) {
                r2_same_base = true;
                break;
              }
            }
            if (r2_same_base) {
              pr.r1.num_filtered_same_base_reads++;
              pr.r1.num_filtered_same_base_bases += orig_r1.length();
              pr.r2.num_filtered_same_base_reads++;
              pr.r2.num_filtered_same_base_bases += orig_r2.length();
              continue;
            }
          }

        } // end filter block

        // Both reads pass all filters; write them
        pr.r1.num_reads++;
        pr.r1.num_bases += orig_r1.length();
        pr.r2.num_reads++;
        pr.r2.num_bases += orig_r2.length();

        pr.r1.read_length_min = min<size_t>(orig_r1.length(), pr.r1.read_length_min);
        pr.r1.read_length_max = max<size_t>(orig_r1.length(), pr.r1.read_length_max);
        pr.r2.read_length_min = min<size_t>(orig_r2.length(), pr.r2.read_length_min);
        pr.r2.read_length_max = max<size_t>(orig_r2.length(), pr.r2.read_length_max);

        out_r1.write_sequence(orig_r1);
        out_r2.write_sequence(orig_r2);

        if (read_file_base_limit) {
          local_current_bases += orig_r1.length() + orig_r2.length();
          if (local_current_bases > read_file_base_limit) {
            pr.reached_eof = false;
            break;
          }
        }
      }

      return pr;
    };

    // ---- Sniff up to FASTQ_FORMAT_SNIFF_READS reads from R1 to predict the quality format ----
    const size_t FASTQ_FORMAT_SNIFF_READS = 10000;
    const int32_t prelim_offset = 64;

    cFastqFile r1_input(r1_file_name.c_str(), fstream::in, num_threads);
    cFastqQualityConverter prelim_fqc("ILLUMINA_1.3+", "SANGER");

    vector<cFastqSequence> buffered_reads;
    buffered_reads.reserve(FASTQ_FORMAT_SNIFF_READS);

    uint8_t sniffed_min_quality_score = 255;

    {
      cFastqSequence seq;
      while ((buffered_reads.size() < FASTQ_FORMAT_SNIFF_READS) && r1_input.read_sequence(seq, prelim_fqc)) {
        for (uint32_t i = 0; i < seq.m_qualities.size(); i++) {
          uint8_t q = static_cast<uint8_t>(seq.m_qualities[i]);
          sniffed_min_quality_score = min(sniffed_min_quality_score, q);
        }
        buffered_reads.push_back(seq);
        seq = cFastqSequence();
      }
    }

    string quality_format = predict_format_from_min_quality_score(sniffed_min_quality_score);

    int32_t numerical_quality_offset_adjustment = prelim_offset - static_cast<int32_t>(format_to_chr_offset[quality_format]);

    size_t buffered_idx = 0;
    auto get_buffered_then_stream_r1 = [&](cFastqSequence& seq) -> bool {
      if (buffered_idx < buffered_reads.size()) {
        seq = std::move(buffered_reads[buffered_idx++]);
        if (seq.m_numerical_qualities && numerical_quality_offset_adjustment != 0) {
          for (uint32_t i = 0; i < seq.m_qualities.size(); i++) {
            seq.m_qualities[i] = static_cast<char>(static_cast<int32_t>(static_cast<uint8_t>(seq.m_qualities[i])) - numerical_quality_offset_adjustment);
          }
        }
        return true;
      }
      cFastqQualityConverter fqc_for_offset(quality_format, "SANGER");
      return r1_input.read_sequence(seq, fqc_for_offset);
    };

    cFastqFile r2_input(r2_file_name.c_str(), fstream::in, num_threads);
    cFastqQualityConverter r2_initial_fqc(quality_format, "SANGER");
    auto get_r2_stream = [&](cFastqSequence& seq) -> bool {
      return r2_input.read_sequence(seq, r2_initial_fqc);
    };

    PairResult result = run_pass_paired(quality_format, get_buffered_then_stream_r1, get_r2_stream, current_read_file_bases);

    // ---- Check format from R1's actual min quality score (same logic as normalize_fastq) ----
    string corrected_quality_format = predict_format_from_min_quality_score(result.r1.overall_min_quality_score);

    cFastqFile* tail_count_r1_stream = &r1_input;
    cFastqFile* tail_count_r2_stream = &r2_input;
    unique_ptr<cFastqFile> fresh_r1, fresh_r2;

    if (corrected_quality_format != quality_format) {
      cerr << "    Warning: Quality score format predicted from the first " << FASTQ_FORMAT_SNIFF_READS << " reads (" << quality_format
           << ") did not match the format determined from the entire file (" << corrected_quality_format << ")." << endl;
      cerr << "    Re-converting FASTQ files using corrected quality score format..." << endl;

      quality_format = corrected_quality_format;

      fresh_r1.reset(new cFastqFile(r1_file_name.c_str(), fstream::in, num_threads));
      fresh_r2.reset(new cFastqFile(r2_file_name.c_str(), fstream::in, num_threads));

      cFastqQualityConverter fresh_r1_fqc(quality_format, "SANGER");
      cFastqQualityConverter fresh_r2_fqc(quality_format, "SANGER");

      auto get_fresh_r1 = [&](cFastqSequence& seq) -> bool {
        return fresh_r1->read_sequence(seq, fresh_r1_fqc);
      };
      auto get_fresh_r2 = [&](cFastqSequence& seq) -> bool {
        return fresh_r2->read_sequence(seq, fresh_r2_fqc);
      };

      result = run_pass_paired(quality_format, get_fresh_r1, get_fresh_r2, initial_current_read_file_bases);
      tail_count_r1_stream = fresh_r1.get();
      tail_count_r2_stream = fresh_r2.get();
    }

    // If stopped early, count remaining reads in both files for stats accuracy
    uint64_t num_filtered_coverage_limit_reads_r1 = 0;
    uint64_t num_filtered_coverage_limit_bases_r1 = 0;
    uint64_t num_filtered_coverage_limit_reads_r2 = 0;
    uint64_t num_filtered_coverage_limit_bases_r2 = 0;

    if (!result.reached_eof) {
      cFastqQualityConverter tail_fqc("ILLUMINA_1.3+", "SANGER");
      cFastqSequence seq;
      while (tail_count_r1_stream->read_sequence(seq, tail_fqc)) {
        num_filtered_coverage_limit_reads_r1++;
        num_filtered_coverage_limit_bases_r1 += seq.length();
      }
      while (tail_count_r2_stream->read_sequence(seq, tail_fqc)) {
        num_filtered_coverage_limit_reads_r2++;
        num_filtered_coverage_limit_bases_r2 += seq.length();
      }
    }

    // Convert raw min/max quality scores to SANGER for reporting
    cFastqQualityConverter final_fqc(quality_format, "SANGER");

    auto build_summary = [&](const FileResult& fr,
                              uint64_t ncov_reads, uint64_t ncov_bases,
                              const string& convert_name) -> AnalyzeFastqSummary {
      cFastqSequence mm_seq;
      mm_seq.m_qualities.append(1, (char)fr.overall_min_quality_score);
      mm_seq.m_qualities.append(1, (char)fr.overall_max_quality_score);
      final_fqc.convert_sequence(mm_seq);
      uint8_t min_q = (uint8_t)mm_seq.m_qualities[0] - format_to_chr_offset["SANGER"];
      uint8_t max_q = (uint8_t)mm_seq.m_qualities[1] - format_to_chr_offset["SANGER"];

      uint64_t num_orig_reads = fr.num_original_reads + ncov_reads;
      uint64_t num_orig_bases = fr.num_original_bases + ncov_bases;
      double avg = (fr.num_reads > 0) ? static_cast<double>(fr.num_bases) / static_cast<double>(fr.num_reads) : 0.0;

      return AnalyzeFastqSummary(
          fr.read_length_min,
          fr.read_length_max,
          avg,
          num_orig_reads,
          fr.num_filtered_too_short_reads,
          fr.num_filtered_same_base_reads,
          fr.num_filtered_too_many_N_reads,
          ncov_reads,
          fr.num_reads,
          min_q,
          max_q,
          num_orig_bases,
          fr.num_bases,
          false,   // reads_were_split: not supported in paired mode
          quality_format,
          "SANGER",
          convert_name
      );
    };

    AnalyzeFastqSummary s_r1 = build_summary(result.r1, num_filtered_coverage_limit_reads_r1, num_filtered_coverage_limit_bases_r1, r1_convert_file_name);
    AnalyzeFastqSummary s_r2 = build_summary(result.r2, num_filtered_coverage_limit_reads_r2, num_filtered_coverage_limit_bases_r2, r2_convert_file_name);

    // Print combined summary for paired files
    uint64_t total_orig_reads = s_r1.num_original_reads + s_r2.num_original_reads;
    uint64_t total_orig_bases = s_r1.num_original_bases + s_r2.num_original_bases;
    uint64_t total_reads = s_r1.num_reads + s_r2.num_reads;
    uint64_t total_bases = s_r1.num_bases + s_r2.num_bases;

    uint32_t width_for_reads = to_string(total_orig_reads).size();
    uint32_t width_for_bases = to_string(total_orig_bases).size();

    cerr << "    Original base quality format: " << quality_format << " New format: SANGER" << endl;
    cerr << "    Original reads: " << total_orig_reads << " bases: " << total_orig_bases << endl;

    uint64_t total_cov_filtered = num_filtered_coverage_limit_reads_r1 + num_filtered_coverage_limit_reads_r2;
    if (total_cov_filtered) {
      cerr << "    Filtered reads: " << setw(width_for_reads) << total_cov_filtered;
      cerr << " bases: " << setw(width_for_bases) << (num_filtered_coverage_limit_bases_r1 + num_filtered_coverage_limit_bases_r2);
      cerr << " (coverage limit option)" << endl;
    }

    if (filter_reads) {
      uint64_t total_short = result.r1.num_filtered_too_short_reads + result.r2.num_filtered_too_short_reads;
      uint64_t total_N = result.r1.num_filtered_too_many_N_reads + result.r2.num_filtered_too_many_N_reads;
      uint64_t total_same = result.r1.num_filtered_same_base_reads + result.r2.num_filtered_same_base_reads;

      if (total_short + total_N + total_same + total_cov_filtered == 0) {
        cerr << "    Filtered reads: none" << endl;
      } else {
        if (total_short) {
          cerr << "    Filtered reads: " << setw(width_for_reads) << total_short;
          cerr << " bases: " << setw(width_for_bases) << (result.r1.num_filtered_too_short_bases + result.r2.num_filtered_too_short_bases);
          cerr << " (<" << _read_length_min << " bases long)" << endl;
        }
        if (total_N) {
          string percentage = formatted_double(100 * _max_N_fraction, 1).to_string();
          cerr << "    Filtered reads: " << setw(width_for_reads) << total_N;
          cerr << " bases: " << setw(width_for_bases) << (result.r1.num_filtered_too_many_N_bases + result.r2.num_filtered_too_many_N_bases);
          cerr << " (≥" << percentage << "% bases N)" << endl;
        }
        if (total_same) {
          string percentage = formatted_double(100 * _max_same_base_fraction, 1).to_string();
          cerr << "    Filtered reads: " << setw(width_for_reads) << total_same;
          cerr << " bases: " << setw(width_for_bases) << (result.r1.num_filtered_same_base_bases + result.r2.num_filtered_same_base_bases);
          cerr << " (≥" << percentage << "% same base)" << endl;
        }
      }
    }
    cerr << "    Analyzed reads: " << setw(width_for_reads) << total_reads << " bases: " << setw(width_for_bases) << total_bases << endl;

    return {s_r1, s_r2};
  }

  // converts a sequence file
  void convert_fastq(const string &from_file_name, const string &to_file_name, const string &from_format, const string &to_format, bool _reverse_complement)
  {
    cFastqFile input_fastq_file(from_file_name.c_str(), ios::in);
    cFastqFile output_fastq_file(to_file_name.c_str(), ios::out);

    cFastqQualityConverter fqc(from_format, to_format);

    cFastqSequence on_sequence;
    while (input_fastq_file.read_sequence(on_sequence, fqc)) 
    {
      fqc.convert_sequence(on_sequence);
      if (_reverse_complement) 
        on_sequence = reverse_complement(on_sequence);
      output_fastq_file.write_sequence(on_sequence);
    }
    
  }
  
  void fastq_sequence_trim_end_on_base_quality(cFastqSequence& seq, const uint32_t base_quality)
  {
    for (uint32_t i=0; i<seq.m_qualities.size(); i++) {
      if (static_cast<uint8_t>(seq.m_qualities[i]-33) < base_quality) {
        seq.m_sequence.resize(i);
        seq.m_qualities.resize(i);
        break; 
      }
    }
  }
  
  bool cFastqSequence::identical(cFastqSequence& seq)
  {
    return ( (this->m_sequence == seq.m_sequence) && (this->m_qualities == seq.m_qualities) );
  }


  // constructor
  cFastqQualityConverter::cFastqQualityConverter(const string &_from_quality_format, const string &_to_quality_format)
  {
    // Set up maps between formats
    map<string,uint8_t> format_to_chr_offset;
    format_to_chr_offset["SANGER"] = 33;
    format_to_chr_offset["SOLEXA"] = 64;
    format_to_chr_offset["ILLUMINA_1.3+"] = 64;
    
    map<string,string> format_to_quality_type;
    format_to_quality_type["SANGER"] = "PHRED";
    format_to_quality_type["SOLEXA"] = "SOLEXA";
    format_to_quality_type["ILLUMINA_1.3+"] = "PHRED";
    
    from_quality_format = _from_quality_format;
    to_quality_format = _to_quality_format;
    
    // check what we asked for is valid...
    ASSERT(format_to_chr_offset.count(from_quality_format), 
           "Unknown FASTQ quality score format: " + from_quality_format + "\nValid choices are 'SANGER', 'SOLEXA', 'ILLUMINA_1.3+', 'NUMERICAL'");
    ASSERT(format_to_chr_offset.count(to_quality_format), 
           "Unknown FASTQ quality score format: " + to_quality_format + "\nValid choices are 'SANGER', 'SOLEXA', 'ILLUMINA_1.3+', 'NUMERICAL'");

    
    from_quality_type = format_to_quality_type[from_quality_format];
    to_quality_type = format_to_quality_type[to_quality_format];

    from_chr_offset = format_to_chr_offset[from_quality_format];
    to_chr_offset = format_to_chr_offset[to_quality_format];

    
    this->resize(256);
    for (uint16_t i = 0; i<=255; i++) {
      (*this)[i] = 0;
    }
    
    for (uint16_t from_chr = 0; from_chr<=255; from_chr++) {

      int32_t from_quality = from_chr - from_chr_offset;
      
      // Calculate the probability of error
      double probability_of_error;
      
      if (from_quality_type == "SOLEXA") {
        probability_of_error = 1 / (1+pow(10,(double)from_quality/10));
      } else if (from_quality_type == "PHRED") {
        probability_of_error = pow(10,-(double)from_quality/10);
      } else {
        cerr << "Unknown base quality score type: " << from_quality_type << endl;
        exit(-1);
      }
      
      //Convert back to quality score
      int32_t to_quality;
            
      if (to_quality_type == "SOLEXA") {
        to_quality = static_cast<uint32_t>(round(10 * log((1-probability_of_error)/probability_of_error) / log(10)));
      } else if (to_quality_type == "PHRED") {
        to_quality = static_cast<uint32_t>(round(-10 * log(probability_of_error) / log(10)));
      } else {
        cerr << "Unknown base quality score type: " << to_quality_type << endl;
        exit(-1);
      }
            
      int16_t to_chr = to_quality + to_chr_offset;
      
      // May be out of range
      if ((to_chr < 0) || (to_chr > 255)) continue;

      (*this)[(uint8_t)from_chr] = (uint8_t)to_chr;
      
      // Debug
      //cerr << from_chr << " => " << to_chr << endl;
    }     
    
  }

  void cFastqQualityConverter::convert_sequence(cFastqSequence &seq) {
    
    for(uint32_t i=0; i < seq.m_qualities.size(); i++)
    {
      seq.m_qualities[i] = (*this)[seq.m_qualities[i]];
    }
  }
  
  string cFastqQualityConverter::predict_fastq_file_format(const string& file_name, uint64_t& num_original_reads, uint64_t& num_original_bases, uint32_t& read_length_min, uint32_t& read_length_max, uint8_t& min_quality_score, uint8_t& max_quality_score)
  {
  // Initialize the input variables!
    num_original_reads = 0;
    num_original_bases = 0;
    read_length_min = numeric_limits<uint32_t>::max();
    read_length_max = 0;
    min_quality_score = 255;
    max_quality_score = 0;
    
  // Set up maps between formats
  map<string,uint8_t> format_to_chr_offset;
  format_to_chr_offset["SANGER"] = 33;
  format_to_chr_offset["SOLEXA"] = 64;
  format_to_chr_offset["ILLUMINA_1.3+"] = 64;
    
  cFastqFile input_fastq_file(file_name.c_str(), fstream::in);
  input_fastq_file.m_check_for_repeated_read_names = true;
  
  cFastqSequence on_sequence;
  cFastqQualityConverter prelim_fqc("ILLUMINA_1.3+", "SANGER");
  
  while (input_fastq_file.read_sequence(on_sequence, prelim_fqc)) {
    
    //increment read number
    num_original_reads++;
    
    //check sequence length
    read_length_min = min<uint32_t>(read_length_min, on_sequence.m_sequence.size());
    read_length_max = max<uint32_t>(read_length_max, on_sequence.m_sequence.size());

    
    //add current sequence length to number of bases
    num_original_bases += on_sequence.m_sequence.size();
      
      //iterate through sequence grabbing the associated scores
    for (uint32_t i=0; i<on_sequence.m_qualities.size(); i++) {
      int this_score(uint8_t(on_sequence.m_qualities[i]));
      if( this_score > max_quality_score ) max_quality_score = this_score;
        if( this_score < min_quality_score ) min_quality_score = this_score;
    }
  }
  
  // Default is SANGER
  string quality_format = "SANGER";
  
  // Typical range: (-5, 40) + 64
  if (min_quality_score >= format_to_chr_offset["SOLEXA"] - 5) {
    quality_format = "SOLEXA";
  } 
  // Typical range:  (0, 40) + 64
  if (min_quality_score >= format_to_chr_offset["ILLUMINA_1.3+"]) {
    quality_format = "ILLUMINA_1.3+";
  }
    
    return quality_format;
  }
  
  //constructor
  cFastqFile::cFastqFile()
    : flexgzfstream()
    , m_current_line(0)
    , m_file_name("")
    , m_check_for_repeated_read_names(false)
    , m_last_read_name(""), m_repeated_read_name_count(0)
  {
  }
 
  
  cFastqFile::cFastqFile(const string &file_name, std::ios_base::openmode mode, unsigned int num_threads)
    : flexgzfstream(file_name.c_str(), mode, num_threads)
    , m_current_line(0)
    , m_file_name(file_name)
    , m_check_for_repeated_read_names(false)
    , m_last_read_name("")
    , m_repeated_read_name_count(0)
  {
  }

  // read one sequence record from the file
  bool cFastqFile::read_sequence(cFastqSequence &sequence, cFastqQualityConverter& fqc) {
    
    // We're done, no error
    if (m_stream->eof())
     return false; 

    
    uint32_t count = 0;
    string line;
        
    memset(sequence.m_base_counts, 0, sizeof(sequence.m_base_counts));
    
    // get the next four lines
    while (count < 4) {
      breseq::getline(*m_stream, line);
      
      m_current_line++;
      
      // Didn't get a first line, then we ended correctly
      if (m_stream->eof()) {
        if (count == 0) {
          return false;
        } else {
          uint32_t last_valid_line = static_cast<uint32_t>(floor((m_current_line-1)/4.0) * 4);
          fprintf(stderr, "Incomplete FASTQ sequence record found at end of file.\nFile %s\nLine: %d\n", m_file_name.c_str(), m_current_line-1);
          fprintf(stderr, "You may be able to repair this damage and salvage the reads before this point with the command:\n");
          fprintf(stderr, "  head -n %u %s > new.fastq\n", last_valid_line, m_file_name.c_str());
          fprintf(stderr, "Then use \"new.fastq\" as input.\n");
          exit(-1);
        }
      }
      
      // Skip empty lines
      if (line.size()==0) continue;
      
      switch (count) {
        case 0:
          if( line[0] != '@' ) {
            fprintf(stderr, "FASTQ sequence record does not begin with @NAME line.\nFile %s\nLine: %d\n", m_file_name.c_str(), m_current_line);
            exit(-1);
          }
          sequence.m_name = line.substr(1,string::npos);
          
          // Delete any sequence name information after the first space...
          // Necessary for scrubbing SRA FASTQs, for example.
          { // block to keep inside this switch case
            size_t space_pos = sequence.m_name.find(" ");
            if (space_pos != string::npos) 
            {
              sequence.m_name.erase(space_pos);
            }
          }
          
          // some SRA files have identical read names, we don't like this...
          if (m_check_for_repeated_read_names)
          {
            string original_read_name = sequence.m_name;
            if (m_last_read_name == sequence.m_name)
            {
              m_repeated_read_name_count++;
              sequence.m_name += "r" + to_string(m_repeated_read_name_count);
            }
            else
            {
              m_repeated_read_name_count = 0;
            }
            m_last_read_name = original_read_name;
          }
          
          break;
          
        case 1:
          sequence.m_sequence = line;
          
          for (uint32_t i=0; i<sequence.m_sequence.size(); i++) {
            
            // convert to uppercase and require
            // reformatting if this was necessary
            switch (sequence.m_sequence[i]) {
                
              case 'A':
              case 'T':
              case 'C':
              case 'G':
                break;
                
              case 'N':
                break;
                
              case 'a':
                sequence.m_sequence.replace(i,1,1,'A');
                break;
                
              case 't':
                sequence.m_sequence.replace(i,1,1,'T');
                break;
                
              case 'c':
                sequence.m_sequence.replace(i,1,1,'C');
                break;
                
              case 'g':
                sequence.m_sequence.replace(i,1,1,'G');
                break;

              case 'n':
                sequence.m_sequence.replace(i,1,1,'N');
                break;
              
              // all other characters converted to 'N'
              default :
                sequence.m_sequence.replace(i,1,1,'N');

            }

            // keep a count of the number of each base for detecting homopolymeric reads
            sequence.m_base_counts[basechar2index(sequence.m_sequence[i])]++;
            
            if(sequence.m_sequence[i] != 'A' && 
               sequence.m_sequence[i] != 'T' && 
               sequence.m_sequence[i] != 'G' && 
               sequence.m_sequence[i] != 'C' && 
               sequence.m_sequence[i] != 'N') {
              
              fprintf(stderr, "FASTQ sequence character not allowed %c.\nSequence: %s\nFile %s\nLine: %d\n", 
                      sequence.m_sequence[i], sequence.m_sequence.c_str(), m_file_name.c_str(), m_current_line);
              exit(-1);
            }
          }
          
          break;
        case 2:
          
          //Only need to see if the first character is a +
          if( line[0] != '+' ) {
            fprintf(stderr, "FASTQ sequence record does not contain +NAME line.\nFile %s\nLine: %d\n", m_file_name.c_str(), m_current_line);
            exit(-1);
          }
          // Could optionally check to see if the name after the + was either absent or identical to the earlier name
          sequence.m_name_plus = line.substr(1,string::npos);

          break;
        case 3:
          
          if (sequence.m_sequence.size() == line.size()) {
            sequence.m_qualities = line;
          } else if ((line.find_first_of(" ") != string::npos) && (line.find_first_not_of(" -0123456789\t") == string::npos)) {
            
            sequence.m_numerical_qualities = true;
            vector<string> numerical_qualities(split(line, " "));
            if( sequence.m_sequence.size() != numerical_qualities.size() ) {
              fprintf(stderr, "FASTQ sequence record has different SEQUENCE and numerical QUALITY lengths.\nFile %s\nLine: %d\n", m_file_name.c_str(), m_current_line);
              exit(-1);
            }
            
            // convert the qualities to characters with the Illumina offset (which keeps things from being negative)
            sequence.m_qualities = "";
            for(vector<string>::iterator it = numerical_qualities.begin(); it != numerical_qualities.end(); it++)
            {
              // use of uint16_t is on purpose to force proper conversion @JEB
              sequence.m_qualities += static_cast<char>(from_string<int16_t>(*it)) + fqc.from_chr_offset;
            }
          } else {
            ERROR("FASTQ QUALITY line length does not match SEQUENCE length.\nFile: " + m_file_name + " Line: " + to_string(m_current_line) + "\nSequence:     " + sequence.m_sequence + "\nQuality Line: " + line);
          }

          break;
      }
      
      count++;
  }
    
    return true;
  }

  void cFastqFile::write_sequence(const cFastqSequence &sequence) {
    *m_stream << "@" << sequence.m_name << endl;
    *m_stream << sequence.m_sequence << endl;
    *m_stream << "+" << sequence.m_name_plus << endl;
    *m_stream << sequence.m_qualities << endl;
  }


  int32_t cSimFastqSequence::SEED_VALUE = time(NULL);
  /*!
    qscore_cumulative_probability_table

    Achieved by:

      Step: Gather frequencies of scores in the following DCAMP read files:
        SRR014475.fastq
        SRR014476.fastq
        SRR014477.fastq
        SRR014478.fastq
        SRR014479.fastq
        SRR014480.fastq
        SRR014481.fastq
        SRR014482.fastq
        SRR014483.fastq
        SRR014484.fastq
        SRR014485.fastq

      NOTE: The scores are in the Sanger format and represent a Phred quality
      score from 0 to 93 using ASCII 33 to 126.
      (http://en.wikipedia.org/wiki/FASTQ_format)

      Step: Normalize each frequency to 1 by dividing each score's frequency by
      the sum of all frequencies that occured.

      Step: Calculate the cumulative probability for each score.

      Step: Multiply each probability by 100,000.

      ***@GRC: Multiplied by 100,000 to achieve precision to the 100,000th
      decimal place. Wasn't sure if C's rand() function would be able to roll
      between 0 and 1 by increment of .00001 so opted for this method.
    */

  map<uint32_t, uint32_t>
  cSimFastqSequence::qscore_cumulative_probability_table =
  make_map<uint32_t, uint32_t>
    (33, 57)
    (34, 96)
    (35, 214)
    (36, 475)
    (37, 985)
    (38, 1501)
    (39, 1759)
    (40, 2017)
    (41, 2275)
    (42, 2532)
    (43, 27931)
    (44, 28185)
    (45, 28439)
    (46, 28693)
    (47, 28946)
    (48, 29201)
    (49, 29457)
    (50, 29714)
    (51, 29971)
    (52, 30230)
    (53, 30491)
    (54, 30753)
    (55, 31016)
    (56, 31280)
    (57, 31545)
    (58, 31809)
    (59, 32074)
    (60, 32337)
    (61, 32600)
    (62, 32861)
    (63, 33121)
    (64, 58268)
    (65, 65171)
    (66, 65422)
    (67, 69288)
    (68, 69533)
    (69, 69775)
    (70, 70013)
    (71, 81668)
    (72, 81899)
    (73, 96356)
    (78, 96412)
    (84, 100000);

  map<char, string> cSimFastqSequence::random_snp_base_options =
  make_map<char, string>
    ('A', "TCG")
    ('T', "ACG")
    ('C', "ATG")
    ('G', "ATC")
    ('N', "ACTG");

  char cSimFastqSequence::random_insertion_base_options[] =
  {'A', 'C', 'T', 'G'};

  char cSimFastqSequence::get_random_quality_score(void)
  {
    uint32_t reserved_offset = qscore_cumulative_probability_table[35];
    
    //Roll between 0 and 99,999.
    uint32_t random_die = rand() % 100000 - reserved_offset;
    random_die += reserved_offset;

    map<uint32_t, uint32_t>::const_iterator it =
        qscore_cumulative_probability_table.begin();

    //Iterate through until random probability is greater then a cumulative
    //probability in the table.
    while (random_die >= it->second) {
      ++it;
    }

    //Found it! Return as a character.
    return char(it->first);
  }

  char cSimFastqSequence::get_random_error_base(const char not_this_base)
  {
    ASSERT(random_snp_base_options.count(not_this_base), "Error!");

    uint32_t size = random_snp_base_options[not_this_base].size();
    return random_snp_base_options[not_this_base][rand() % size];
  }

  char cSimFastqSequence::get_random_insertion_base(void)
  {
    //Roll from 0 to 3.
    return base_char_list[rand() % 4];
  }

  //Return if this particular base is an error given a quality score.
  bool cSimFastqSequence::is_random_error_base(char ascii_qscore)
  {
    const int32_t qscore = int32_t(ascii_qscore) - 33;

    double p_value = pow(10, -qscore/10.0);

    //We want precision to the 100,000th decimal place.
    double p_value_max = p_value * 100000;

    return ((rand() % 100000) <= p_value_max);
  }

  const uint32_t deletion_probability = 100000;//10E-5
  bool cSimFastqSequence::is_random_deletion_base(void)
  {
    //Roll from 0 to deletion_probability.
    return (1 == (rand() % deletion_probability));
  }

  const uint32_t insertion_probability = 100000;//10E-5
  bool cSimFastqSequence::is_random_insertion_base(void)
  {
    //Roll from 0 to insertion_probability.
    return (1 == (rand() % insertion_probability));
  }

  void cSimFastqSequence::GaussianRNG::box_muller_transform(float* z0, float* z1) {
  static const float PI =
  3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706f;
    //Two random variables in the interval (0, 1] with a precision of .001
    float u1 = static_cast<float>((rand() % 1000 + 1) / 1000.f);
    float u2 = static_cast<float>((rand() % 1000 + 1) / 1000.f);

    *z0 = sqrtf(-2.f * log(u1)) * std::cos(2.f * PI * u2); 
    *z1 = sqrtf(-2.f * log(u1)) * std::sin(2.f * PI * u2); 

    return;
  }

  cSimFastqSequence::GaussianRNG::GaussianRNG(int mean, int stdev) 
    : m_mean(mean)
    , m_stdev(stdev) {

    srand(cSimFastqSequence::SEED_VALUE);
    GaussianRNG::box_muller_transform(&m_z0, &m_z1);
    m_z0 = (m_z0 * m_stdev) + m_mean; 
    m_z1 = (m_z1 * m_stdev) + m_mean; 

    return;
  }

  int32_t cSimFastqSequence::GaussianRNG::sample() {
    int32_t ret_val = static_cast<int32_t>(round(m_z0));
    m_store = m_z1;

    GaussianRNG::box_muller_transform(&m_z0, &m_z1);
    m_z0 = (m_z0 * m_stdev) + m_mean; 
    m_z1 = (m_z1 * m_stdev) + m_mean; 

    return ret_val;
  }


  vector<int32_t> cSimFastqSequence::GaussianRNG::samples(int32_t size) {
    vector<int32_t> ret_val(size, 0);
    for (int32_t i = 0; i < size; ++i) {
      ret_val[i] = this->sample();
    }
    return ret_val;
  }


  cFastqSequence cSimFastqSequence::simulate(const cAnnotatedSequence& ref_sequence,
                                             uint32_t start_1,
                                             uint32_t read_size,
                                             int8_t   strand,
                                             uint32_t  id,
                                             uint32_t  n_reads,
                                             bool     verbose) {
    cFastqSequence ret_val;
    sprintf(ret_val.m_name, "READ-%i", id);
    sprintf(ret_val.m_name_plus, "[strand]:%i\t[start_1]:%u", strand, start_1);

    if(verbose && id && (!(id % 10000) || (id == n_reads))){
      ostringstream progress_message;
      progress_message << "\tREAD: " << setw(12) << right << id;
      if (n_reads) progress_message << "/" << n_reads;
      print_progress_line(progress_message.str());
    }

    /*! Algorithm in use:

        1) Grab 2 x read_size segment from a random location in the reference
          sequence. Get the reverse complement of this sequence if it is
          a negative strand.

          ***Note: We treat the sequence circularly such that if we index
          bases from past the end of the sequence we will instead start
          indexing them from the front.

          ***Note: Grab 2x read_size for easier access to the following bases
          when deletions occur and for the impossible chance that each base
          is a deletion.

        2) Iterate through this segment.

        3) At each index/base determine if insertion(INS), deletion(DEL),
          if it's a possible error due to quality score or if the base is
          normal.

           If DEL  : Continue on to next base in segment.
           If INS  : Randomly insert one of the four possible bases.
           If Error: Assign random base other then the one given.
           Normal  : Assign base given in segment.

          Assign the base a random quality score which is determine by the
          qscore_cumulative_probability_table.
    */


    string verbose_errors(read_size, ' ');
    string verbose_deletions (read_size, ' ');
    string verbose_insertions(read_size, ' ');

    //! Algorithm Step 1:
    string ref_segment =
        ref_sequence.get_sequence_1_start_size(start_1, 2 * read_size);

    if (strand == -1) {
      ref_segment = reverse_complement(ref_segment);
    }


    //! Algorithm Step 2:
    //Initializations for iterating through bases.
    ret_val.m_sequence.resize(read_size);
    ret_val.m_qualities.resize(read_size);

    //Index for bases in the simulated read.
    size_t index_to_assign = 0;

    //Index for bases in the reference segment.
    size_t index_in_ref_segment = 0;

    //Begin assigning bases to simulated read.
    while (index_to_assign < read_size) {
    //! Algorithm Step 3:

      //! DEL base.
      if (cSimFastqSequence::is_random_deletion_base()) {
       //For verbose.
        verbose_deletions[index_in_ref_segment] =
            ref_segment[index_in_ref_segment];
        //Continue to next index/base in reference segment.
        ++index_in_ref_segment;
        continue;
      }
      //! INS base.
      else if (cSimFastqSequence::is_random_insertion_base()) {
        char base_to_insert =
            cSimFastqSequence::get_random_insertion_base();

        //Assign base.
        ret_val.m_sequence[index_to_assign] = base_to_insert;

        //Assign quality score.
        ret_val.m_qualities[index_to_assign] =
            cSimFastqSequence::get_random_quality_score();

        //For verbose
        verbose_insertions[index_in_ref_segment] = base_to_insert;

        /*Increment index_to_assign but not index_in_ref_segment and
         continue to next iteration.*/
        ++index_to_assign;
        continue;
      }
      //! Determine if error or normal base.
      else {
        char quality_score =
            cSimFastqSequence::get_random_quality_score();

        //! Error base.
        if (cSimFastqSequence::is_random_error_base(quality_score)) {
          /*Since an error occured we want to assign a base different from
           the one given. */
          char not_this_base = ref_segment[index_in_ref_segment];

          //Assign different base.
          char new_base =
              cSimFastqSequence::get_random_error_base(not_this_base);
          ret_val.m_sequence[index_to_assign] = new_base;

          //Assing quality score.
          ret_val.m_qualities[index_to_assign] = quality_score;

          //For verbose.
          verbose_errors[index_to_assign] = new_base;
        }
        //! Normal base.
        else {
          //Assign base.
          ret_val.m_sequence[index_to_assign] = ref_segment[index_in_ref_segment];

          //Assign quality score.
          ret_val.m_qualities[index_to_assign] = quality_score;
        }
        /*Increment both index_to_assign and index_in_ref_segment and
        continue to next iteration. */
        ++index_to_assign;
        ++index_in_ref_segment;
      }
    } //End assigning bases to simulated read.
    if (verbose) {
      if (verbose_deletions.find_first_not_of(' ')  != string::npos ||
          verbose_insertions.find_first_not_of(' ') != string::npos ) {
        printf("\tVerbose output for simulated read    :  %s\n",
               ret_val.m_name.c_str());

        const string &original =
            ref_sequence.get_sequence_1_start_size(start_1, 2 * read_size);
        printf("\tReference Segment(2 x Read Size)     :  %s\n",
               original.c_str());

        if (strand == -1) {
          printf("\tSimulated negative strand            :  %s\n",
                 ref_segment.c_str());
        }

        if (verbose_errors.find_first_not_of(' ')     != string::npos) {
          printf("\tSimulated errors                     :  %s\n",
                 verbose_errors.c_str());
        }

        if (verbose_deletions.find_first_not_of(' ')  != string::npos) {
          printf("\tSimulated DEL                        :  %s\n",
                 verbose_deletions.c_str());
        }

        if (verbose_insertions.find_first_not_of(' ') != string::npos) {
          printf("\tSimulated INS                        :  %s\n",
                 verbose_insertions.c_str());
        }

        printf("\tFinal simulated read sequence        :  %s\n",
                 ret_val.m_sequence.c_str());
        printf("\tFinal simulated read quality scores  :  %s\n",
                 ret_val.m_qualities.c_str());
        printf("\n");
      }
    }//End verbose output.

    return ret_val;
  }
  void cSimFastqSequence::simulate_single_ends(const cAnnotatedSequence& sequence,
                                               uint32_t n_reads,
                                               uint32_t read_size, 
                                               string file_name,
                                               bool verbose) 
{
    cFastqFile out(file_name.c_str(), ios_base::out);
    
    srand(cSimFastqSequence::SEED_VALUE);
    for (uint32_t i = 0; i < n_reads; ++i) {
      uint32_t start_1 = rand() % sequence.get_sequence_length() + 1;
      int8_t strand    = (rand() % 2) == 0 ? 1 : -1;
      cFastqSequence read = cSimFastqSequence::simulate(sequence,
                                                        start_1,
                                                        read_size,
                                                        strand,
                                                        i + 1,
                                                        n_reads,
                                                        verbose);
      out.write_sequence(read);
    }
    
  }

  void cSimFastqSequence::simulate_paired_ends(const cAnnotatedSequence& sequence,
                                               uint32_t n_reads,
                                               uint32_t read_size, 
                                               uint32_t mean,
                                               uint32_t stdev,
                                               string pair_1_file_name,
                                               string pair_2_file_name,
                                               bool verbose)
  {
    cFastqFile pair_1_out(pair_1_file_name.c_str(), ios_base::out);
    cFastqFile pair_2_out(pair_2_file_name.c_str(), ios_base::out);

    cSimFastqSequence::GaussianRNG random_size(mean, stdev);
    srand(cSimFastqSequence::SEED_VALUE);
    for (uint32_t i = 0; i < n_reads; ++i) {

      //Pair 1
      uint32_t start_1 = rand() % sequence.get_sequence_length() + 1;
      int8_t strand    = 1;
      cFastqSequence pair_1 = cSimFastqSequence::simulate(sequence,
                                                          start_1,
                                                          read_size,
                                                          strand,
                                                          i + 1,
                                                          n_reads,
                                                          verbose);

      //Pair 2
      start_1 = (start_1 + random_size.sample() - (read_size * 2)) % sequence.get_sequence_length();
      strand  = -1;
      cFastqSequence pair_2 = cSimFastqSequence::simulate(sequence,
                                                          start_1,
                                                          read_size,
                                                          strand,
                                                          i + 1,
                                                          n_reads,
                                                          verbose);
      if (rand() % 100 < 50) {
        pair_1.m_name += "/1";
        pair_2.m_name += "/2";
        pair_1_out.write_sequence(pair_1);
        pair_2_out.write_sequence(pair_2);
      } else {
        pair_1.m_name += "/2";
        pair_2.m_name += "/1";
        pair_1_out.write_sequence(pair_2);
        pair_2_out.write_sequence(pair_1);
      }

    }

  }

// Simulates perfectly tiled reads - assuming a circular genome
//
// If coverage is 2 x read length (or greater), then every possible read is simulated
// on top and bottom strands.
//
// if it is less, then every ceiling(2 x length / coverage) bases, reads on both strands
// are simulated on each strand for every start position.
//
// This assumes a circular genome
void cSimFastqSequence::simulate_tiled(const cAnnotatedSequence& sequence,
                                             uint32_t read_size,
                                             uint32_t coverage,
                                             string file_name,
                                             bool verbose)
{
  (void) verbose;
  
  cFastqFile out(file_name.c_str(), ios_base::out);
  vector<int8_t>strands = make_vector<int8_t>(-1)(+1);
  
  uint32_t spacing = ceil( (2.0 * read_size) / coverage );
  spacing = max<uint32_t>(spacing, 0);
  cout << "Read length : " << read_size << endl;
  cout << "Coverage    : " << coverage << endl;
  cout << "Spacing     : " << spacing << endl;
  
  uint32_t read_index = 1;
  for (uint32_t start_1 = 1; start_1 <= sequence.get_sequence_length(); start_1+=spacing) {
    for (vector<int8_t>::iterator strand_it = strands.begin(); strand_it != strands.end(); strand_it++) {
      int8_t strand = *strand_it;
      
      cFastqSequence read;
      sprintf(read.m_name, "READ-%i", read_index++);
      sprintf(read.m_name_plus, "[strand]:%i\t[start_1]:%u", strand, start_1);
      
      read.m_sequence = sequence.get_sequence_1_start_size(start_1, read_size);
      if (strand == -1) {
        read.m_sequence = reverse_complement(read.m_sequence);
      }
      
      read.m_qualities = repeat_char(char(73), read_size);
      
      out.write_sequence(read);
    }
  }

return;
}

// Reverse complement and also uppercase
// Convert most characters to 'N'. Might want to give errors on non-printable characters
  char reverse_complement_lookup_table[256] = {
/*  0*/    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
/* 16*/    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
/* 32*/    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
/* 48*/    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
/* 64*/    'N', 'T', 'N', 'G', 'N', 'N', 'N', 'C', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
/* 80*/    'N', 'N', 'N', 'N', 'A', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
/* 96*/    'N', 'T', 'N', 'G', 'N', 'N', 'N', 'C', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
/*112*/    'N', 'N', 'N', 'N', 'A', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
/*128*/    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
/*144*/    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
/*160*/    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
/*176*/    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
/*192*/    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
/*208*/    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
/*224*/    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
/*240*/    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N'
  };
  
} // breseq namespace

