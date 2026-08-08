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

#include "fastq.h"
#include "reference_sequence.h"

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
      // Two output streams are open at once here, each with its own compressor
      // thread pool. Split the -j budget so the two pools sum to ~num_threads
      // rather than 2*num_threads. (At num_threads==1 this rounds up to 1 each,
      // i.e. 2 total; avoiding that last overshoot would require the two streams
      // to share a single pool.)
      unsigned int const per_stream_threads = std::max<unsigned int>(1, num_threads / 2);
      cFastqFile out_r1(r1_convert_file_name.c_str(), fstream::out, per_stream_threads);
      cFastqFile out_r2(r2_convert_file_name.c_str(), fstream::out, per_stream_threads);

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


  // Default seed is the wall clock, so an unseeded run is NOT reproducible. That is deliberate --
  // it preserves the historical behavior -- but do_simulate_reads prints the resolved seed so any
  // run can be reproduced after the fact, and every test passes --seed explicitly.
  int32_t cSimFastqSequence::SEED_VALUE = time(NULL);

  map<char, string> cSimFastqSequence::random_snp_base_options =
  make_map<char, string>
    ('A', "TCG")
    ('T', "ACG")
    ('C', "ATG")
    ('G', "ATC")
    ('N', "ACTG");

  char cSimFastqSequence::random_insertion_base_options[] =
  {'A', 'C', 'T', 'G'};

  // 1e9 * 10^(-k/10) for k = 0..9. Every Phred score's error probability is one of these ten
  // mantissas divided by a power of ten, so the whole scale is exact integer arithmetic. This
  // replaces a pow(10, -q/10.0) call: libm is not required to be correctly rounded, and a last-ulp
  // difference between glibc and Apple libm would change which bases become errors and so change
  // the simulated reads between platforms. See portable_random.h.
  static const uint32_t kPhredMantissaPerBillion[10] = {
    1000000000u, 794328235u, 630957344u, 501187234u, 398107171u,
     316227766u, 251188643u, 199526231u, 158489319u, 125892541u
  };

  static uint32_t phred_error_rate_per_billion(uint32_t phred)
  {
    uint32_t decades = phred / 10;
    uint32_t rate = kPhredMantissaPerBillion[phred % 10];
    for (uint32_t i = 0; i < decades; i++) rate /= 10;
    return rate;
  }

  bool cSimFastqSequence::is_error_base(cPortableRandom& rng, uint32_t phred)
  {
    return rng.next_event_per_billion(phred_error_rate_per_billion(phred));
  }

  char cSimFastqSequence::get_random_error_base(cPortableRandom& rng, const char not_this_base)
  {
    ASSERT(random_snp_base_options.count(not_this_base), "Unexpected base: " + to_string(not_this_base));

    const string& options = random_snp_base_options[not_this_base];
    return options[rng.next_below(static_cast<uint32_t>(options.size()))];
  }

  char cSimFastqSequence::get_random_insertion_base(cPortableRandom& rng)
  {
    return base_char_list[rng.next_below(4)];
  }

  char cSimQualityModel::draw(cPortableRandom& rng, uint32_t cycle, uint32_t read_size) const
  {
    // Real Illumina quality is roughly flat over most of the read and falls off fastest over the
    // last few cycles, so the decline is quadratic in the cycle rather than linear. Integer
    // throughout: the division truncates identically everywhere.
    uint32_t n = (read_size > 1) ? (read_size - 1) : 1;
    uint32_t drop = (phred_start > phred_end) ? (phred_start - phred_end) : 0;
    uint32_t shift = static_cast<uint32_t>(
      (static_cast<uint64_t>(drop) * cycle * cycle) / (static_cast<uint64_t>(n) * n));

    int32_t center = static_cast<int32_t>(phred_start) - static_cast<int32_t>(shift);
    int32_t q = rng.next_gaussian(center, static_cast<int32_t>(phred_stdev));

    if (q < static_cast<int32_t>(phred_min)) q = static_cast<int32_t>(phred_min);
    if (q > static_cast<int32_t>(phred_max)) q = static_cast<int32_t>(phred_max);

    return static_cast<char>(q + 33);   // Sanger offset
  }

  void cSimFastqSequence::warn_if_quality_format_ambiguous(char min_quality_char)
  {
    // normalize_fastq classifies a file from the minimum quality character in it: >= 59 is taken to
    // be SOLEXA and every score is shifted by 31. A library with no low-quality tail is therefore
    // silently reinterpreted, which is very hard to diagnose downstream.
    if (min_quality_char >= 59) {
      WARN("The lowest quality score written was '" + to_string(min_quality_char) +
           "' (Phred " + to_string(static_cast<int32_t>(min_quality_char) - 33) + "). breseq infers a "
           "FASTQ's quality encoding from the lowest score in the file and will misread this one as "
           "SOLEXA rather than SANGER. Lower --quality-end or --quality-min so that the file "
           "contains at least one score below Phred 26.");
    }
  }

  cFastqSequence cSimFastqSequence::simulate_read(const cAnnotatedSequence& ref_sequence,
                                                  int32_t  left_1,
                                                  uint32_t read_size,
                                                  int8_t   strand,
                                                  const cSimQualityModel& quality_model,
                                                  const cSimErrorModel&   error_model,
                                                  cPortableRandom& rng,
                                                  const string& read_name,
                                                  char& min_quality_char,
                                                  bool verbose)
  {
    cFastqSequence ret_val;
    ret_val.m_name = read_name;
    ret_val.m_name_plus = "[strand]:" + to_string(static_cast<int32_t>(strand)) +
                          "\t[left_1]:" + to_string(left_1);

    /*! Algorithm:
        1) Fetch the read's reference footprint [left_1, left_1+read_size-1], plus a tail of up to
           read_size more bases in the direction the read is built, so that simulated deletions have
           somewhere to draw replacement bases from. Reverse-complement the whole buffer for the
           minus strand, which puts the read's 5' end first.
        2) Walk the buffer, deciding per base whether it is a deletion, an insertion, a substitution
           error, or a normal base.
        3) Give every emitted base a quality score from the per-cycle model, and derive the
           substitution probability from that score.
     */

    // Buffer origin: the plus-strand read starts at its footprint; the minus-strand read's buffer
    // must START read_size earlier, because after reverse-complementing, the footprint ends up at
    // the FRONT. Fetching an equal tail in each case keeps the deletion slack symmetric.
    int32_t seq_len = static_cast<int32_t>(ref_sequence.get_sequence_length());
    int32_t fetch_start = (strand == -1) ? (left_1 - static_cast<int32_t>(read_size)) : left_1;

    // On a circular sequence get_sequence_1_start_size wraps; on a linear one it silently CLAMPS and
    // returns a shorter string. The old code indexed the returned buffer blindly and so could read
    // past the end of it for any read near a linear sequence's end -- both undefined behavior and a
    // reproducibility hazard, since the bytes past the end are allocator-dependent. Clamp here
    // instead, and let the loop below cope with a short buffer.
    string ref_segment;
    if (ref_sequence.is_circular()) {
      ref_segment = ref_sequence.get_sequence_1_start_size(fetch_start, 2 * read_size);
    } else {
      int32_t lo = max<int32_t>(1, fetch_start);
      int32_t hi = min<int32_t>(seq_len, fetch_start + 2 * static_cast<int32_t>(read_size) - 1);
      if (lo <= hi) ref_segment = ref_sequence.get_sequence_1(lo, hi);
    }

    if (strand == -1) {
      ref_segment = reverse_complement(ref_segment);
    }

    ASSERT(ref_segment.size() >= read_size,
           "Simulated read at " + to_string(left_1) + " would extend past the end of sequence " +
           ref_sequence.m_seq_id + ". The caller must place reads so their whole footprint fits.");

    string verbose_errors(read_size, ' ');
    string verbose_deletions(ref_segment.size(), ' ');
    string verbose_insertions(ref_segment.size(), ' ');

    ret_val.m_sequence.resize(read_size);
    ret_val.m_qualities.resize(read_size);

    size_t index_to_assign = 0;       // position in the simulated read
    size_t index_in_ref_segment = 0;  // position in the reference buffer

    while (index_to_assign < read_size) {

      // A deletion consumes a reference base without emitting one, so it is only possible while
      // there is still buffer left to consume. Near the end of a linear sequence there may not be.
      bool can_delete = (index_in_ref_segment + 1 < ref_segment.size());

      if (can_delete && rng.next_event_per_billion(error_model.deletion_rate_per_million * 1000)) {
        verbose_deletions[index_in_ref_segment] = ref_segment[index_in_ref_segment];
        ++index_in_ref_segment;
        continue;
      }

      char quality_char = quality_model.draw(rng, static_cast<uint32_t>(index_to_assign), read_size);
      if (quality_char < min_quality_char) min_quality_char = quality_char;
      uint32_t phred = static_cast<uint32_t>(quality_char) - 33;

      if (rng.next_event_per_billion(error_model.insertion_rate_per_million * 1000)) {
        // An insertion emits a base without consuming one from the reference.
        char base_to_insert = cSimFastqSequence::get_random_insertion_base(rng);
        ret_val.m_sequence[index_to_assign] = base_to_insert;
        ret_val.m_qualities[index_to_assign] = quality_char;
        verbose_insertions[index_in_ref_segment] = base_to_insert;
        ++index_to_assign;
        continue;
      }

      char ref_base = ref_segment[index_in_ref_segment];

      if (cSimFastqSequence::is_error_base(rng, phred)) {
        char new_base = cSimFastqSequence::get_random_error_base(rng, ref_base);
        ret_val.m_sequence[index_to_assign] = new_base;
        verbose_errors[index_to_assign] = new_base;
      } else {
        ret_val.m_sequence[index_to_assign] = ref_base;
      }
      ret_val.m_qualities[index_to_assign] = quality_char;

      ++index_to_assign;
      ++index_in_ref_segment;
    }

    if (verbose) {
      if (verbose_deletions.find_first_not_of(' ')  != string::npos ||
          verbose_insertions.find_first_not_of(' ') != string::npos ||
          verbose_errors.find_first_not_of(' ')     != string::npos) {
        printf("\tVerbose output for simulated read    :  %s\n", ret_val.m_name.c_str());
        printf("\tReference segment                    :  %s\n", ref_segment.c_str());
        if (verbose_errors.find_first_not_of(' ')     != string::npos)
          printf("\tSimulated errors                     :  %s\n", verbose_errors.c_str());
        if (verbose_deletions.find_first_not_of(' ')  != string::npos)
          printf("\tSimulated DEL                        :  %s\n", verbose_deletions.c_str());
        if (verbose_insertions.find_first_not_of(' ') != string::npos)
          printf("\tSimulated INS                        :  %s\n", verbose_insertions.c_str());
        printf("\tFinal simulated read sequence        :  %s\n", ret_val.m_sequence.c_str());
        printf("\tFinal simulated read quality scores  :  %s\n", ret_val.m_qualities.c_str());
        printf("\n");
      }
    }

    return ret_val;
  }

  void cSimFastqSequence::simulate_single_ends(const cAnnotatedSequence& sequence,
                                               uint32_t n_reads,
                                               uint32_t read_size,
                                               uint32_t seed,
                                               const cSimQualityModel& quality_model,
                                               const cSimErrorModel& error_model,
                                               string file_name,
                                               bool verbose)
  {
    cFastqFile out(file_name.c_str(), ios_base::out);

    int32_t seq_len = static_cast<int32_t>(sequence.get_sequence_length());
    char min_quality_char = 126;

    for (uint32_t i = 0; i < n_reads; ++i) {
      // One independent stream per read, so a read depends only on its own index. Adding a draw
      // elsewhere then cannot shift every subsequent read.
      cPortableRandom rng(seed, i);

      int8_t strand = rng.next_coin() ? +1 : -1;

      // Place the whole read inside the sequence. The old code drew a start anywhere in the
      // sequence and let the fetch clamp, which produced short reads (and out-of-bounds indexing).
      int32_t left_1 = sequence.is_circular()
                         ? static_cast<int32_t>(rng.next_in_range_1(1, seq_len))
                         : static_cast<int32_t>(rng.next_in_range_1(1, seq_len - read_size + 1));

      cFastqSequence read = simulate_read(sequence, left_1, read_size, strand,
                                          quality_model, error_model, rng,
                                          "READ-" + to_string(i + 1), min_quality_char, verbose);
      out.write_sequence(read);

      if (verbose && ((i + 1) % 10000 == 0 || (i + 1) == n_reads)) {
        ostringstream progress_message;
        progress_message << "\tREAD: " << setw(12) << right << (i + 1) << "/" << n_reads;
        print_progress_line(progress_message.str());
      }
    }

    warn_if_quality_format_ambiguous(min_quality_char);
  }

  void cSimFastqSequence::simulate_paired_ends(const cAnnotatedSequence& sequence,
                                               uint32_t n_pairs,
                                               uint32_t read_size,
                                               uint32_t fragment_mean,
                                               uint32_t fragment_stdev,
                                               uint32_t seed,
                                               const cSimQualityModel& quality_model,
                                               const cSimErrorModel& error_model,
                                               string pair_1_file_name,
                                               string pair_2_file_name,
                                               bool verbose)
  {
    cFastqFile pair_1_out(pair_1_file_name.c_str(), ios_base::out);
    cFastqFile pair_2_out(pair_2_file_name.c_str(), ios_base::out);

    int32_t seq_len = static_cast<int32_t>(sequence.get_sequence_length());
    char min_quality_char = 126;

    for (uint32_t i = 0; i < n_pairs; ++i) {
      cPortableRandom rng(seed, i);

      // Fragment length. Clamping at read_size (rather than at 2*read_size) is what allows the mates
      // to OVERLAP, which is a perfectly normal library -- tests/data/tmv_plasmid is 190 bp
      // fragments read 2 x 151. It also removes the unsigned underflow the old code had when the
      // drawn length came out below 2*read_size.
      int32_t L = rng.next_gaussian(static_cast<int32_t>(fragment_mean),
                                    static_cast<int32_t>(fragment_stdev));
      L = max<int32_t>(L, static_cast<int32_t>(read_size));
      L = min<int32_t>(L, seq_len);

      // Leftmost base of the fragment. On a linear sequence the whole fragment must fit inside, just
      // as a library prepared from linear DNA would; a circular one may wrap.
      int32_t left = sequence.is_circular()
                       ? static_cast<int32_t>(rng.next_in_range_1(1, seq_len))
                       : static_cast<int32_t>(rng.next_in_range_1(1, seq_len - L + 1));

      // FR ("innie"): each mate reads inward from one end of the fragment. Which physical end is
      // read 1 is a coin flip, so both strands are equally represented in R1 -- previously every R1
      // was a plus-strand read. This does NOT change the pair's orientation label, which breseq
      // derives from 5'-end order rather than R1/R2 identity (see
      // best_pair_orientation_and_distance in candidate_junctions.cpp).
      bool r1_is_left = rng.next_coin();
      int32_t left_read_left_1  = left;
      int32_t right_read_left_1 = left + L - static_cast<int32_t>(read_size);

      string name = "READ-" + to_string(i + 1);

      cFastqSequence r1 = simulate_read(sequence,
                                        r1_is_left ? left_read_left_1 : right_read_left_1,
                                        read_size, r1_is_left ? +1 : -1,
                                        quality_model, error_model, rng,
                                        name + "/1", min_quality_char, verbose);
      cFastqSequence r2 = simulate_read(sequence,
                                        r1_is_left ? right_read_left_1 : left_read_left_1,
                                        read_size, r1_is_left ? -1 : +1,
                                        quality_model, error_model, rng,
                                        name + "/2", min_quality_char, verbose);

      // ALWAYS in this order: breseq pairs the two files positionally, reading them in lockstep.
      pair_1_out.write_sequence(r1);
      pair_2_out.write_sequence(r2);

      if (verbose && ((i + 1) % 10000 == 0 || (i + 1) == n_pairs)) {
        ostringstream progress_message;
        progress_message << "\tPAIR: " << setw(12) << right << (i + 1) << "/" << n_pairs;
        print_progress_line(progress_message.str());
      }
    }

    warn_if_quality_format_ambiguous(min_quality_char);
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
                                             const cSimQualityModel& quality_model,
                                             string file_name,
                                             bool verbose)
{
  (void) verbose;

  // Tiled reads are deterministic by construction -- fixed spacing, no RNG, no simulated errors --
  // so this mode reproduces regardless of the seed. Qualities come from the shared model rather than
  // a hard-coded constant, but with the flat settings this mode uses (start == end, stdev 0) the
  // model consumes no randomness and emits the same character for every base, as before.
  cPortableRandom unused_rng(0, 0);

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
      
      read.m_qualities.resize(read_size);
      for (uint32_t c = 0; c < read_size; c++)
        read.m_qualities[c] = quality_model.draw(unused_rng, c, read_size);

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

