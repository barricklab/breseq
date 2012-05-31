/*****************************************************************************

AUTHORS

  Jeffrey E. Barrick <jeffrey.e.barrick@gmail.com
  David B. Knoester

LICENSE AND COPYRIGHT

  Copyright (c) 2008-2010 Michigan State University
  Copyright (c) 2011 The University of Texas at Austin

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
  // Each subset of entries is the effect of a change to ('A')('T')('C')('G').
  enum BaseSubstitutionEffect {
    intergenic,
    noncoding,
    synonymous,
    nonsynonymous,
    no_change
  };
  typedef vector<BaseSubstitutionEffect> SequenceBaseSubstitutionEffects;
  
  // 
  // List is as long as the genome
  enum BaseCDSStrand {
    no_CDS,
    forward,
    reverse,
    conflict
  };
  typedef vector<BaseCDSStrand> SequenceBaseCDSStrands;
  
  //<! Class for annotating effects of base substitution mutations
  //   that properly accounts for overlapping features.
  class BaseSubstitutionEffects {    
  public:
        
    static map<string,uint8_t>  base_change_to_index;
    static vector<string>       base_change_list;
    static vector<string>       base_pair_change_list;
    static map<string,string>   base_change_to_base_pair_change;
    static vector<string>       snp_types;
    static map<string,uint8_t>  nt_type_list;
    
    map<string,SequenceBaseSubstitutionEffects> m_bse;
    map<string,SequenceBaseCDSStrands> m_bcs;
    
    void initialize_from_sequence(cReferenceSequences& ref_seq_info);
    
  }; // class BaseSubstitutionEffects

  
  class BaseSubstitutionEffectCounts : public map<string, map<string, uint32_t> > {
  public:
    void initialize_totals(cReferenceSequences& ref_seq_info, BaseSubstitutionEffects& bse);
    
    void add_position_1_to_totals(cReferenceSequences& ref_seq_info, BaseSubstitutionEffects& bse, string seq_id, uint32_t on_pos);
    void subtract_position_1_from_totals(cReferenceSequences& ref_seq_info, BaseSubstitutionEffects& bse, string seq_id, uint32_t on_pos);
    void change_position_1_totals(cReferenceSequences& ref_seq_info, BaseSubstitutionEffects& bse, string seq_id, uint32_t pos_1, int32_t inc);
    
    void add_base_pair_change_to_totals(cReferenceSequences& ref_seq_info, BaseSubstitutionEffects& bse, string seq_id, uint32_t pos_1, string base_pair_change);

  };
  
  void MutationCountFile(
                         cReferenceSequences& ref_seq_info, 
                         vector<cGenomeDiff>& genome_diffs, 
                         string& output_file_name, 
                         bool base_substitution_statistics
                         );
  
  
} // namespace breseq

#endif
