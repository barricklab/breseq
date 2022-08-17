/*****************************************************************************

AUTHORS

  Jeffrey E. Barrick <jeffrey.e.barrick@gmail.com
  David B. Knoester

LICENSE AND COPYRIGHT

  Copyright (c) 2008-2010 Michigan State University
  Copyright (c) 2011-2022 The University of Texas at Austin

  breseq is free software; you can redistribute it and/or modify it under the  
  terms the GNU General Public License as published by the Free Software 
  Foundation; either version 1, or (at your option) any later version.

*****************************************************************************/

#ifndef _BRESEQ_CANDIDATE_JUNCTIONS_H_
#define _BRESEQ_CANDIDATE_JUNCTIONS_H_

#include "common.h"

#include "alignment.h"
#include "reference_sequence.h"

using namespace std;

namespace breseq {

  const string junction_name_separator = "__";
  const string junction_user_defined_tag = "UD";
  
  /*! information about one side of a new junction 
   */
  class JunctionSide
  {
  public:
    string seq_id;
    int32_t position;   // 1-indexed
    int32_t strand;     // -1 or +1
    int32_t redundant;
    
    // Extended properties for resolve_alignments.cpp
    cFeatureLocation* is;
    string is_interval;
		string is_side_key;
    
    int32_t overlap;    // This is the amount of junction overlap that was given
                        // to this side when assigning it to each side. Must be >0.
    
    JunctionSide()
    : is(NULL)
    {
      position = 0;
      strand = 0;
      redundant = 0;
      overlap = 0;
    }
    
    JunctionSide(const string& _seq_id, int32_t _position, int32_t _strand, int32_t _redundant = false)
    : is(NULL)
    {
      seq_id = _seq_id;
      position = _position;
      strand = _strand;
      redundant = _redundant;
      overlap = 0;
    }
    
    bool operator ==(const JunctionSide& side) const
    {
      return (this->seq_id == side.seq_id) && (this->position == side.position) && (this->strand == side.strand);
    }
    
  };
  
	class JunctionInfo
	{
  public:
		JunctionSide sides[2];
    
		int32_t alignment_overlap;
		string unique_read_sequence;
		int32_t flanking_left;
		int32_t flanking_right;
    bool user_defined;
    
    //! Create empty JunctionInfo
    JunctionInfo()
    {
      alignment_overlap = 0;
      flanking_left = 0;
      flanking_right = 0;
      user_defined = false;
    }
    
    //! Create a JunctionInfo from supplied parameters
    JunctionInfo(const JunctionSide& _side_1, const JunctionSide& _side_2, int32_t _alignment_overlap, const string& _unique_read_sequence, int32_t _flanking_left=0, int32_t _flanking_right=0, bool _user_defined = false)
    {       
      sides[0] = _side_1;
      sides[1] = _side_2;
      
      alignment_overlap = _alignment_overlap;
      unique_read_sequence = _unique_read_sequence;
      flanking_left = _flanking_left;
      flanking_right = _flanking_right;
      user_defined = _user_defined;
    }
    
    //! Create JunctionInfo from genome diff entry
    JunctionInfo(cDiffEntry & de)
    {
      alignment_overlap = from_string<int32_t>(de["overlap"]);
      unique_read_sequence = "";
      if (de.count("unique_read_sequence")) 
        unique_read_sequence = de["unique_read_sequence"];
      flanking_left = from_string<int32_t>(de["flanking_left"]);
      flanking_right = from_string<int32_t>(de["flanking_right"]);
      
      sides[0].seq_id = de["side_1_seq_id"];
      sides[0].position = from_string<int32_t>(de["side_1_position"]);
      sides[0].strand = from_string<int32_t>(de["side_1_strand"]);
      
      sides[0].redundant = false;
      if (de.count("side_1_redundant")) {
        sides[0].redundant = from_string<int32_t>(de["side_1_redundant"]); 
      }
      
      sides[1].seq_id = de["side_2_seq_id"];
      sides[1].position = from_string<int32_t>(de["side_2_position"]);
      sides[1].strand = from_string<int32_t>(de["side_2_strand"]);
      
      sides[1].redundant = false;
      if (de.count("side_2_redundant")) {
        sides[1].redundant = from_string<int32_t>(de["side_2_redundant"]); 
      }
      
      user_defined = false;
    }
    //! Deserializes JunctionInfo from a key string
    JunctionInfo(const string& junction_name)
    {
      vector<string> s = split(junction_name, junction_name_separator);
      
      JunctionSide side_1(s[0], from_string<int32_t>(s[1]), from_string<int32_t>(s[2]), from_string<int32_t>(s[10]));
      if (side_1.strand == 0) side_1.strand = -1;
      
      JunctionSide side_2(s[3], from_string<int32_t>(s[4]), from_string<int32_t>(s[5]), from_string<int32_t>(s[11]));
      if (side_2.strand == 0) side_2.strand = -1;
      
      bool user_defined = false;
      if (s.size() == 13) user_defined = true;
      
      // last portion only there if user defined
      JunctionInfo retval(
                          side_1, 
                          side_2,
                          from_string<int32_t>(s[6]),
                          s[7],
                          from_string<int32_t>(s[8]),
                          from_string<int32_t>(s[9]),
                          user_defined
                          );
            
      *this = retval;
    }
    
    
    // Serializes a JunctionInfo to a string
    string junction_key(bool include_redundant_tags = true)
    {
      ASSERT( (sides[0].strand == +1) || (sides[0].strand == -1), "side 1 strand uninitialized or wrong: must be -1/+1");
      ASSERT( (sides[1].strand == +1) || (sides[1].strand == -1), "side 2 strand uninitialized or wrong: must be -1/+1");
      
      // allocate vector of correct size
      vector<string> values(12); 
      
      // place values in vector
      values[0] = sides[0].seq_id;
      values[1] = to_string(sides[0].position);
      values[2] = to_string(sides[0].strand);
      
      values[3] = sides[1].seq_id;
      values[4] = to_string(sides[1].position);
      values[5] = to_string(sides[1].strand);
      
      values[6] = to_string(alignment_overlap);
      values[7] = to_string(unique_read_sequence);
      values[8] = to_string(flanking_left);
      values[9] = to_string(flanking_right);
 
      //@JEB: Including these can lead to identical user candidate junctions having the same key
      
      if (include_redundant_tags) {
        values[10] = to_string(sides[0].redundant);
        values[11] = to_string(sides[1].redundant);
      }
      
      if (user_defined) values.push_back(junction_user_defined_tag);
      
      return join(values, junction_name_separator);
    }
    
    bool operator ==(const JunctionInfo& junction) const
    {
      return (this->sides[0] == junction.sides[0]) && (this->sides[1] == junction.sides[1]);
    }
	};

  // Predefinitition
  class JunctionCandidate;
  typedef counted_ptr<JunctionCandidate> JunctionCandidatePtr;
  typedef map<string, JunctionCandidatePtr> KeyToJunctionCandidateMap;
  typedef map<string, KeyToJunctionCandidateMap > SequenceToKeyToJunctionCandidateMap;
  
  class JunctionCandidate : public JunctionInfo {
  public:
    string sequence;
    string reverse_complement_sequence;
		map<uint32_t, uint32_t> read_begin_hash;
    vector<JunctionCandidatePtr> merged_from;    // Keeps track of all junctions that were merged to create this (including itself)
    
    JunctionCandidate() {}
    
    JunctionCandidate(
                      const JunctionInfo& _junction_info, 
                      const string& _sequence
                      ) 
    : JunctionInfo(_junction_info)
    , sequence(_sequence)
    {
      reverse_complement_sequence = reverse_complement(_sequence);
    }
    
    size_t pos_hash_score() const
    {
      return read_begin_hash.size();
    }
    
    size_t num_matching_reads() const
    {
      size_t num(0);
      for(map<uint32_t, uint32_t>::const_iterator it=read_begin_hash.begin(); it!=read_begin_hash.end(); it++ ) {
        num += it->second;
      }
      return num;
    }
    
    // Sort by unique coordinate, then redundant (or second unique) coordinate to get reliable ordering for output
    static bool sort_by_score_unique_coord(const JunctionCandidate& a, const JunctionCandidate &b)
    {
      int32_t a_uc = a.sides[0].position;
      int32_t a_rc = a.sides[1].position;
      if (a.sides[0].redundant != 0) swap(a_uc, a_rc);
      
      int32_t b_uc = b.sides[0].position;
      int32_t b_rc = b.sides[1].position;
      if (b.sides[0].redundant != 0) swap(b_uc, b_rc);
      
      if (b.pos_hash_score() != a.pos_hash_score())
        return (b.pos_hash_score() < a.pos_hash_score());
      else if (a_uc != b_uc)
        return (a_uc < b_uc);
      else
        return (a_rc < b_rc);
    }
    
    static bool sort_by_scores_and_seq_length(const JunctionCandidate& a, const JunctionCandidate &b)
    {
      if (b.pos_hash_score() != a.pos_hash_score())
        return (b.pos_hash_score() < a.pos_hash_score());
      else
        return (a.sequence.size() < b.sequence.size());
    }
    
    static bool sort_by_ref_seq_coord(const JunctionCandidate& a, const JunctionCandidate &b)
    {
      //TODO: Uncomment this code after supplying it with a ref_seq_info with a seq_order field
      /*if (ref_seq_info.seq_order[acj.sides[0].seq_id] != ref_seq_info.seq_order[bcj.sides[0].seq_id])
       return (ref_seq_info.seq_order[acj.sides[0].seq_id] < ref_seq_info.seq_order[bcj.sides[0].seq_id]);
       else*/
      return (a.sides[0].position < b.sides[0].position);
    }
    
    //! Returns true then we merge into cj, otherwise we merge into this
    bool operator <(const JunctionCandidate& cj) const
    {
      // we want to merge into the shorter sequence (which means less overlap)
      if ( this->sequence.size() > cj.sequence.size() ) return true;  // merge into shorter cj
      if ( this->sequence.size() < cj.sequence.size() ) return false;
      
      // sequence lengths are equal
      // we want to merge into the one that is within the same reference sequence fragment
      int32_t equal_seq_1 = (this->sides[0].seq_id == this->sides[1].seq_id) ? 1 : 0;
      int32_t equal_seq_2 = (cj.sides[0].seq_id == cj.sides[1].seq_id) ? 1 : 0;
      if (equal_seq_1 < equal_seq_2) return true; // merge into cj if only it is on the same fragment
      if (equal_seq_1 > equal_seq_2) return false;
      
      // we want to merge into the one with the closest coordinates?
      int32_t distance_1 = abs( this->sides[0].position - this->sides[1].position );
      int32_t distance_2 = abs( cj.sides[0].position - cj.sides[1].position );
      if (distance_1 > distance_2) return true; // merge into cj if it is separated by a smaller distance
      if (distance_1 < distance_2) return false;
      
      return (this->sides[0].position > cj.sides[0].position); // merge into cj if it has a smaller coord
    }
  };
  


  
  
  uint32_t eligible_read_alignments(
                                    const Settings& settings, 
                                    const cReferenceSequences& ref_seq_info, 
                                    alignment_list& alignments,
                                    bool keep_suboptimal_matches = false,
                                    int32_t min_match_score = -1 // OFF
                                    );
  
	bool test_read_alignment_requirements(
                                        const Settings& settings, 
                                        const cReferenceSequences& ref_seq_info, 
                                        const alignment_wrapper& a
                                        );
  
  class PreprocessAlignments
  {
  public:
    /*! Preprocesses alignments
		 */
		static void preprocess_alignments(Settings& settings, Summary& summary, const cReferenceSequences& ref_seq_info);
    static void split_alignments_on_indels(const Settings& settings, Summary& summary, tam_file& PSAM, int32_t min_indel_split_len, const alignment_list& alignments);

    static void split_matched_and_unmatched_alignments(
                                                       uint32_t fastq_file_index,
                                                       string fasta_file_name, 
                                                       string input_sam_file_name, 
                                                       string matched_sam_file_name, 
                                                       string unmatched_fastq_file_name
                                                       );
    static void split_matched_alignments(uint32_t fastq_file_index, 
                                         string fasta_file_name, 
                                         string input_sam_file_name,
                                         string matched_sam_file_name);
    
    static void merge_sort_sam_files(
                                     string input_sam_file_name_1,
                                     string input_sam_file_name_2,
                                     string output_sam_file_name
                                     );
  };
  
  /*! Structure for storing and testing whether alignment pairs support new junctions
   */
  class AlignmentPair
  {
  public:
    bam_alignment	a1;
    bam_alignment a2;
    int32_t hash_coord;
    
    int32_t a1_unique_start;
    int32_t a1_unique_end;
    int32_t a1_unique_length;
    
    int32_t a2_unique_start;
    int32_t a2_unique_end;    
    int32_t a2_unique_length;

    int32_t end_to_end_length;
    int32_t union_length;
    int32_t intersection_length;
    
    bool pass;
    
    // the constructor - also calculates statistics
    AlignmentPair(bam_alignment& _a1, bam_alignment& _a2, const Settings& settings);
    
    static bool reverse_sort_by_overlap (const AlignmentPair& i,const AlignmentPair&  j) 
      { return (i.intersection_length > j.intersection_length); }
    
    string junction_key(const cReferenceSequences& ref_seq_info);
    
  private:
    void calculate_union_and_unique();
    
    //! determines whether a read pair passes
    bool test(const Settings& settings);
    
  };
  
	class CandidateJunctions
	{
	public:

		/*! Predicts candidate junctions
		 */
		static void identify_candidate_junctions(
                                             const Settings& settings, 
                                             Summary& summary, 
                                             const cReferenceSequences& ref_seq_info
                                             );

    static map<string,cDiffEntry> load_user_junctions (
                                     const Settings& settings,
                                     const Summary& summary,
                                     const cReferenceSequences& ref_seq_info
                                     );
    
    static void normalize_junction_overlap (
                                            const cReferenceSequences& ref_seq_info,
                                            cDiffEntry& jc
                                            );
  
    
    static string construct_junction_sequence( 
                                              const cReferenceSequences& ref_seq_info,
                                              cDiffEntry& jc,
                                              int32_t flanking_length,
                                              bool inclusive_overlap = false
                                              );
  
  private:

    
    /*! Sorting functions
     */
    typedef std::map<string, string> map_t;

    static void _by_ref_seq_coord(map_t a, map_t b, map_t ref_seq_info);
		static void _by_score_unique_coord(map_t a, map_t b);

    
    //! Merge two Junction candidates, updating redundancy. Returns the one that was merged into, then the one that was merged from
    static bool merge_candidate_junctions(JunctionCandidatePtr*& jcp1, JunctionCandidatePtr*& jcp2);
    
    static 	bool alignment_pair_to_candidate_junction(
                                                  const Settings& settings, 
                                                  Summary& summary, 
                                                  const cReferenceSequences& ref_seq_info, 
                                                  AlignmentPair& ap,
                                                  JunctionCandidatePtr& returned_junction_candidate
                                                  );
    
		static uint64_t alignments_to_candidate_junctions(
                                                  const Settings& settings, 
                                                  Summary& summary,  
                                                  const cReferenceSequences& ref_seq_info, 
                                                  SequenceToKeyToJunctionCandidateMap& candidate_junctions, 
                                                  alignment_list& alignments
                                                  );
    
    
	}; // class CandidateJunction

} // namespace breseq

#endif
