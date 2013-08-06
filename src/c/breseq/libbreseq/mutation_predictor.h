/*****************************************************************************

AUTHORS

  Jeffrey E. Barrick <jeffrey.e.barrick@gmail.com
  David B. Knoester

LICENSE AND COPYRIGHT

  Copyright (c) 2008-2010 Michigan State University
  Copyright (c) 2011-2012 The University of Texas at Austin

  breseq is free software; you can redistribute it and/or modify it under the  
  terms the GNU General Public License as published by the Free Software 
  Foundation; either version 1, or (at your option) any later version.

*****************************************************************************/

#ifndef _BRESEQ_MUTATION_PREDICTOR_H_
#define _BRESEQ_MUTATION_PREDICTOR_H_

#include "common.h"

#include "alignment.h"
#include "reference_sequence.h"
#include "genome_diff.h"
#include "settings.h"
#include "fastq.h"
#include "candidate_junctions.h"

using namespace std;

namespace breseq {

	class MutationPredictor
	{
	public:

		static cReferenceSequences ref_seq_info;

		MutationPredictor(cReferenceSequences& ref_seq_info);
    
    void prepare_junctions(Settings& settings, Summary& summary, cGenomeDiff& gd);
		void predict(Settings& settings, Summary& summary, cGenomeDiff& gd);

		static bool sort_by_hybrid(const counted_ptr<cDiffEntry>& a, const counted_ptr<cDiffEntry>& b);
		static bool sort_by_reject_score(const counted_ptr<cDiffEntry>& a, const counted_ptr<cDiffEntry>& b);
		static bool sort_by_pos(const counted_ptr<cDiffEntry>& a, const counted_ptr<cDiffEntry>& b);

	private:

		cSequenceFeature* within_repeat(string seq_id, int32_t position);

	}; // class MutationPredictor
  
  
  // effects of base substitutions
  // List is 4 times as long as the genome, 
  // with a slot or all possible nucleotide changes (including no change) at a position
  // Each subset of entries is the effect of a change to ('A')('T')('C')('G'), respectively.
  enum BaseSubstitutionEffect {
    intergenic_base_substitution,
    noncoding_base_substitution,
    pseudogene_base_substitution,
    synonymous_base_substitution,
    nonsynonymous_base_substitution,
    no_change_base_substitution
  };
  typedef vector<BaseSubstitutionEffect> SequenceBaseSubstitutionEffects;
  
  enum BaseType {
    intergenic_base,
    noncoding_base,
    pseudogene_base,
    protein_base
  };
  
  // 
  // List is as long as the genome
  enum BaseCDSStrand {
    no_CDS,
    forward,
    reverse,
    conflict
  };
  typedef vector<BaseCDSStrand> SequenceBaseCDSStrands;
  
  // Organizes information about a single position for easier counting
  struct BaseSubstitutionEffectPositionInfo {
    base_char                       m_base_char;
    BaseType                        m_base_type;
    vector<BaseSubstitutionEffect>  m_base_substitution_effect; // For change to ('A')('T')('C')('G')
    BaseCDSStrand                   m_base_cds_strand;          // For change to ('A')('T')('C')('G')
  };
  
  //<! Class for annotating effects of base substitution mutations
  //   that properly accounts for overlapping features.
  class BaseSubstitutionEffects {    
  public:
        
    static string                   separator;
    static vector<base_char>        base_char_list;
    static map<base_char, uint8_t>  base_char_to_base_index;
    static map<string,uint8_t>      base_change_to_base_pair_change_index;
    static vector<string>           base_change_list;
    static vector<string>           base_pair_change_list;
    static map<string,string>       base_change_to_base_pair_change;
    static vector<string>           base_change_type_list;
    static vector<string>           base_type_list;
    static map<BaseSubstitutionEffect,BaseType>  snp_type_to_base_type;
    
    map<string,SequenceBaseSubstitutionEffects> m_bse;
    map<string,SequenceBaseCDSStrands> m_bcs;
    
    void initialize_from_sequence(cReferenceSequences& ref_seq_info);
    
    BaseSubstitutionEffectPositionInfo position_info_1(
                                                       cReferenceSequences& ref_seq_info, 
                                                       string seq_id, 
                                                       uint32_t pos_1
                                                       );
    
  }; // class BaseSubstitutionEffects

  
  class BaseSubstitutionEffectCounts  {
  public:
    void initialize_possible_totals(cReferenceSequences& ref_seq_info, BaseSubstitutionEffects& bse);
    void change_position_1_possible_totals(cReferenceSequences& ref_seq_info, BaseSubstitutionEffects& bse, string seq_id, uint32_t pos_1, int32_t inc);
    void change_position_1_observed_totals(cReferenceSequences& ref_seq_info, BaseSubstitutionEffects& bse, string seq_id, uint32_t pos_1, string new_base, int32_t inc);

    static vector<string> base_pair_change_count_list; 
    static vector<string> base_change_type_count_list;

    struct BaseTypeCounts : public map<string,int32_t> {
      
      BaseTypeCounts()
      { 
        (*this)["nt"] = 0; (*this)["gc"] = 0; (*this)["at"] = 0;
      }
    };
    
    // map by base pair change
    struct BasePairChangeCounts : public map<string,int32_t> {
      
      BasePairChangeCounts()
      { 
        for (vector<string>::iterator it = BaseSubstitutionEffects::base_pair_change_list.begin();
             it != BaseSubstitutionEffects::base_pair_change_list.end(); ++it) {
          (*this)[*it] = 0;
        }
        (*this)["TOTAL"] = 0;
      }
    }; //class BaseSubstitutionEffectCounts
    
    map<BaseType, BaseTypeCounts > m_base_counts;
    map<string, BasePairChangeCounts > m_possible_base_pair_change_counts;
    map<string, BasePairChangeCounts > m_observed_base_pair_change_counts;

  };
  
  void MutationCountFile(
                         cReferenceSequences& ref_seq_info, 
                         vector<cGenomeDiff>& genome_diffs, 
                         string& output_file_name, 
                         bool base_substitution_statistics,
                         bool verbose
                         );
  
  
} // namespace breseq

#endif
