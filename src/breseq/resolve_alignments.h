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

#ifndef _BRESEQ_RESOLVE_ALIGNMENTS_H_
#define _BRESEQ_RESOLVE_ALIGNMENTS_H_

#include "common.h"

#include "reference_sequence.h"
#include "genome_diff.h"
#include "calculate_trims.h"
#include "candidate_junctions.h"
#include "pileup_base.h"


using namespace std;

namespace breseq {

  // Pileup class for fetching reads that align across from start to end
  class junction_read_counter : pileup_base {
  public:
    junction_read_counter(const string& bam, const string& fasta, bool verbose, const string& trim_file_name);

    uint32_t count(
                   const string& seq_id,
                   const int32_t start,
                   const int32_t end,
                   const map<string,bool> ignore_read_names,
                   map<string,bool>& counted_read_names
                   );

    // Weighting context for the counts. Each counted read is split between the two hypotheses by
    // junction_read_weight() (stats.h): weight w to "came from the junction", 1-w to "came from
    // the reference". delta_sign is +1 when this counter reads the junction BAM (where a record's
    // own log-likelihood X8 belongs to the junction and the stamped X7 to the reference) and -1
    // when it reads the reference BAM, where the two swap. With weighting off every read scores
    // w = 1 toward whichever BAM it sits in and the sums reproduce the raw counts exactly.
    //
    // The sums produced here use a flat prior; they are a starting point, and the caller refits
    // them against the junction frequency. Only the recorded per-read evidence is needed for that.
    void set_weighting(bool enabled, int32_t delta_sign);

    // Ignore bases at or below this quality when deciding how far a read reaches, so a read is not
    // credited with spanning the counting window on the strength of a miscalled tail. 0 keeps raw
    // alignment bounds.
    void set_base_quality_cutoff(uint32_t cutoff) { _base_quality_cutoff = cutoff; }


    // Valid after count(). Sums run over exactly the reads that count() tallied, and always add
    // to that raw count -- each read's evidence is split between the two hypotheses but stays in
    // this window, so both sums share the window's register normalization.
    //
    // "Weight" is always named from the point of view of the BAM being read: sum_weight() is the
    // evidence supporting whichever hypothesis this BAM represents, sum_complement() the evidence
    // for the other one. So on the junction BAM sum_weight() is junction support, while on the
    // reference BAM it is reference support.
    double sum_weight() const { return _sum_weight; }
    double sum_weight_sq() const { return _sum_weight_sq; }       // for a Kish effective depth
    double sum_complement() const { return _sum_complement; }

    // Per-counted-read evidence, in the order counted: the log10 likelihood ratio in favour of
    // THE JUNCTION (not of this window -- the reference window's sign flip is already folded in
    // when the value is recorded), and how many places in the reference the read maps equally well.
    // have_odds is false when the read carried no competing hypothesis, which is NOT the same as a
    // ratio of zero (that would mean "the two hypotheses explain it equally well").
    // Kept so the caller can re-weight without a second pileup pass -- the mixture below iterates,
    // and re-fetching the same reads once per iteration would multiply the cost of counting.
    struct read_evidence_t { double log10_odds_for_junction; uint32_t n_ref_placements; bool have_odds; };
    const vector<read_evidence_t>& counted_read_evidence() const { return _read_evidence; }

    // Counts how many of the candidate reference start positions for a hypothetical,
    // gapless, average-length read spanning [window_start, window_end] would still
    // span it after being shrunk by the trims at that position (see fetch_callback).
    uint32_t count_confident_overlap_registers(
                   const string& seq_id,
                   const int32_t window_start,
                   const int32_t window_end,
                   const uint32_t read_length_avg
                   ) const;

    virtual void fetch_callback ( const alignment_wrapper& a );

  protected:
    uint32_t _count;
    int32_t _start;
    int32_t _end;
    map<string,bool> _ignore_read_names, _counted_read_names;
    bool _verbose;
    SequenceTrimsList _trims_list;

    bool _weighting_enabled;
    uint32_t _base_quality_cutoff;
    int32_t _delta_sign;
    double _sum_weight;
    double _sum_weight_sq;
    double _sum_complement;
    vector<read_evidence_t> _read_evidence;
  };

  class ResolveJunctionInfo : public JunctionInfo
  {
  public:
    // Extended properties for resolve_alignments.cpp
		string key;
		int32_t overlap;
		uint32_t unique_side;
		uint32_t is_side;
    
    ResolveJunctionInfo() : JunctionInfo() {};
    
    ResolveJunctionInfo(const string& junction_name) 
    : JunctionInfo(junction_name)
    , key(junction_name)
    , overlap(0)
    , is_side(0)
    {
    }
  };
  
  class JunctionTestInfo {
  public:
    int32_t max_left;
    int32_t max_left_minus;
    int32_t max_left_plus;
    int32_t max_right;
    int32_t max_right_minus;
    int32_t max_right_plus;
    int32_t max_min_right;
    int32_t max_min_right_minus;
    int32_t max_min_right_plus;
    int32_t max_min_left;
    int32_t max_min_left_minus;
    int32_t max_min_left_plus;
    uint32_t coverage_minus;
    uint32_t coverage_plus;
    uint32_t total_non_overlap_reads;
    uint32_t pos_hash_score;
    uint32_t max_pos_hash_score;
    uint32_t unique_matches_size;
    uint32_t repeat_matches_size;
    bool has_reads_with_both_different_start_and_end;
    bool redundant_1;
    bool redundant_2;
    string junction_id;
    double neg_log10_pos_hash_p_value;
    uint32_t side_1_continuation;
    uint32_t side_2_continuation;
    vector<string> reject_reasons;

    bool operator <(const JunctionTestInfo& _in) const
    {
      // sort by pos_hash_score, unique_matches_size, repeat_matches_size
      
      if (this->pos_hash_score != _in.pos_hash_score) {
        return (this->pos_hash_score < _in.pos_hash_score);
      }
      
      if (this->unique_matches_size != _in.unique_matches_size) {
        return (this->unique_matches_size < _in.unique_matches_size);
      }
      
      return (this->repeat_matches_size < _in.repeat_matches_size);
    }
  };
  
    
  typedef map<string, JunctionTestInfo> JunctionTestInfoMap;
  
	class JunctionMatch
	{
  public:
    alignment_list reference_alignments;
		alignment_list junction_alignments;
		uint32_t fastq_file_index;
		int32_t mapping_quality_difference;
		uint32_t degenerate_count;
    // The read's best alignment score against the reference genome. Kept alongside
    // mapping_quality_difference (which is best-junction MINUS this) because the difference alone
    // is best over ALL candidate junctions: to score a read against the ONE candidate under
    // consideration we need this baseline plus that candidate's own X5 tag. See score_junction().
    int32_t best_reference_score;

    JunctionMatch() : fastq_file_index(0), mapping_quality_difference(0), degenerate_count(0), best_reference_score(0) {}

    JunctionMatch(
                    const alignment_list& _reference_alignments,
                    const alignment_list&  _junction_alignments,
                    uint32_t _fastq_file_index,
                    int32_t _mapping_quality_difference,
                    uint32_t _degenerate_count,
                    int32_t _best_reference_score = 0
                    )
          :reference_alignments(_reference_alignments)
          ,junction_alignments(_junction_alignments)
          ,fastq_file_index(_fastq_file_index)
          ,mapping_quality_difference(_mapping_quality_difference)
          ,degenerate_count(_degenerate_count)
          ,best_reference_score(_best_reference_score)
    { }
	};

	struct VectorSize {
		string junction_id; uint32_t size; uint32_t size2;
    VectorSize(string _junction_id, uint32_t _size, uint32_t _size2) 
      : junction_id(_junction_id), size(_size), size2(_size2) {}
		static bool sort_reverse_by_size(const VectorSize& lhs, const VectorSize& rhs) {
			if (lhs.size != rhs.size) return lhs.size > rhs.size;
      if (lhs.size2 != rhs.size2) return lhs.size2 > rhs.size2;
      return lhs.junction_id < rhs.junction_id;
		}
	};
  
  typedef counted_ptr<JunctionMatch> JunctionMatchPtr;
  typedef map<string, vector<JunctionMatchPtr> > UniqueJunctionMatchMap;      // Map of junction_id to list of MatchedJunctions
  typedef map<string, map<string, JunctionMatchPtr> > RepeatJunctionMatchMap; // Map of junction_id to read_name to MatchedJunction
  
  
  class PosHashProbabilityTable {
  public:
        
    struct Parameters {
      double negative_binomial_size;
      double negative_binomial_prob;
      double chance_per_pos_strand_no_read_start;
      double average_coverage;
      // Empirical unique-coverage null for the pos_hash (skew) p-value. When use_empirical is true,
      // probability() marginalizes coverage over coverage_hist (the actual per-position histogram)
      // instead of the negative binomial. Small references fall back to the nbinom (use_empirical=false)
      // because the empirical p-value is floored near 1/coverage_hist_total. See
      // use_empirical_pos_hash_coverage() in coverage_distribution.h.
      bool use_empirical;
      int32_t deletion_floor;         // coverage <= this is treated as deletion (excluded from the null)
      double coverage_hist_total;     // N = # non-deletion positions (sum of coverage_hist above the floor)
      vector<double> coverage_hist;   // n[c] = # reference positions at unique coverage c (index 0 unused)
    };

    map<string, Parameters> param;
    uint32_t average_read_length;
    map<string, map<uint32_t, map<uint32_t, double> > > probability_table;
    
    PosHashProbabilityTable(Summary& summary, const Settings& settings);
    
    double probability(string& seq_id, uint32_t pos_hash_score, uint32_t overlap);
  };

  
  class pos_hash_p_value_table : public vector<vector<double> > {
  public:
    pos_hash_p_value_table() {};
    
    pos_hash_p_value_table(const string& in_file_name)
    {      
      // file may not have been created if fitting failed
      if (!file_exists(in_file_name.c_str())) return;
      
      ifstream in_file(in_file_name.c_str());
      ASSERT(!in_file.fail(), "Could not open file: " + in_file_name);
      
      string line; 
      while (!in_file.eof()) {
        getline(in_file, line);
        if (line == "") continue;
        vector<string> line_list = split(line, "\t");
        
        vector<double> converted_line;
        for (vector<string>::iterator it=line_list.begin(); it!=line_list.end(); it++) {
          converted_line.push_back(from_string<double>(*it));
        }
        
        (*this).push_back(converted_line);
      }
    }
  };
  
  void calculate_continuation(
                              ResolveJunctionInfo& rji, 
                              cReferenceSequences& ref_seq_info, 
                              cReferenceSequences& junction_ref_seq_info, 
                              uint32_t& side_1_continuation,
                              uint32_t& side_2_continuation
                              );
  
  void resolve_alignments(
                          Settings& settings,
                          Summary& summary,
                          cReferenceSequences& ref_seq_info,
                          const bool junction_prediction,
                          cReadFileSets &read_files
                          );
  
  // A discordant read pair whose IS (redundant) side maps to several near-identical repeat copies,
  // held aside during the streaming resolve pass so the specific copy can be chosen by a global
  // per-locus vote in the post-streaming merge-back. Every best-score copy of the redundant mate is
  // kept (redundant coverage preserved when written) plus the single unique mate; the vote only
  // decides WHICH copy carries the discordant-pair flags.
  struct HeldDiscordantPair {
    alignment_list unique_alignments;      // the unique mate (exactly one alignment)
    alignment_list redundant_alignments;   // the IS mate (all best-score copies)
    uint32_t unique_fastq_file_index;
    uint32_t redundant_fastq_file_index;
    string   unique_seq_id;
    int32_t  unique_position;              // reference_start_1 of the unique mate (clustering key)
    double   window;                       // per-read-set clustering window (paired distance cutoff)
  };

  void load_junction_alignments(
                                const Settings& settings,
                                Summary& summary,
                                cReadFileSets& read_files,
                                cReferenceSequences& ref_seq_info,
                                cReferenceSequences& junction_ref_seq_info,
                                SequenceTrimsList& trims_list,
                                map<string,uint32_t>& all_junction_ids,
                                bool junction_prediction,
                                const vector<ResolveJunctionInfo>& junction_info_list,
                                UniqueJunctionMatchMap& unique_junction_match_map,
                                RepeatJunctionMatchMap& repeat_junction_match_map,
                                bam_file& resolved_reference_tam,
                                vector<HeldDiscordantPair>& held_discordant_pairs
                                );
  
  void load_sam_only_alignments(
                           const Settings& settings, 
                           Summary& summary, 
                           cReadFileSets& read_files, 
                           cReferenceSequences& ref_seq_info,
                           SequenceTrimsList& trims_list,
                           bam_file& resolved_reference_tam
                                );
  
  bool alignment_overlaps_junction(const vector<ResolveJunctionInfo>& junction_info_list, const alignment_wrapper& in_a);

  
  void score_junction(
                      const Settings& settings, 
                      Summary& summary, 
                      const string& junction_id, 
                      UniqueJunctionMatchMap& matched_junction_ref, 
                      RepeatJunctionMatchMap& degenerate_matches_ref, 
                      bam_file& junction_tam,
                      JunctionTestInfo& junction_test_info, 
                      vector<ResolveJunctionInfo>& junction_info_list,
                      cReferenceSequences& ref_seq_info, 
                      cReferenceSequences& junction_ref_seq_info
                      );
  
  void resolve_junction(
                        const Settings& settings,
                        Summary& summary,
                        cReferenceSequences& ref_seq_info,
                        cReferenceSequences& junction_ref_seq_info,
                        const SequenceTrimsList& trims_list,
                        const SequenceTrimsList& junction_trims_list,
                        const string& junction_id,
                        UniqueJunctionMatchMap& unique_junction_match_map,
                        RepeatJunctionMatchMap& repeat_junction_match_map,
                        bam_file& reference_tam,
                        bam_file& junction_tam,
                        bool failed,
                        bool has_non_overlap_alignment
                        );

	void _write_reference_matches(
                                const Settings& settings, 
                                Summary& summary,
                                cReferenceSequences& ref_seq_info, 
                                const SequenceTrimsList& trim_list, 
                                alignment_list& reference_alignments, 
                                bam_file& reference_tam, 
                                uint32_t fastq_file_index
                                );
                                

	vector<string> get_sorted_junction_ids(
                                         UniqueJunctionMatchMap& unique_map, 
                                         RepeatJunctionMatchMap& degenerate_map, 
                                         const vector<string>& keys
                                         );

  
  cDiffEntry junction_to_diff_entry(
                                    const string& key, 
                                    cReferenceSequences& ref_seq_info, 
                                    JunctionTestInfo& test_info
                                    );

  void  assign_one_junction_read_counts(
                                        const Settings& settings,
                                        Summary& summary,
                                        cDiffEntry& j,
                                        const counted_ptr<junction_read_counter>& reference_jrc,
                                        const counted_ptr<junction_read_counter>& junction_jrc
                                        );
  
  // Builds the junction/reference read counters with identical configuration for every caller.
  void make_junction_read_counters(const Settings& settings,
                                   counted_ptr<junction_read_counter>& reference_jrc,
                                   counted_ptr<junction_read_counter>& junction_jrc);

  void  assign_junction_read_counts(
                                    const Settings& settings,
                                    Summary& summary,
                                    cGenomeDiff& gd
                                    );
  


}

#endif
