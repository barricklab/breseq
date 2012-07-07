
/*****************************************************************************

AUTHORS

  Jeffrey E. Barrick <jeffrey.e.barrick@gmail.com>
  David B. Knoester

LICENSE AND COPYRIGHT

  Copyright (c) 2008-2010 Michigan State University
  Copyright (c) 2011 The University of Texas at Austin

  breseq is free software; you can redistribute it and/or modify it under the  
  terms the GNU General Public License as published by the Free Software 
  Foundation; either version 1, or (at your option) any later version.

*****************************************************************************/


#include "libbreseq/candidate_junctions.h"

#include "libbreseq/fastq.h"

using namespace std;

namespace breseq {

  /*! Utility for sorting alignments by the number of mismatches
   */
  
  class alignment_score_map_class : public map<bam_alignment*,uint32_t> {
  public:
    inline bool operator() (alignment_list::iterator a1, alignment_list::iterator a2) { return (*this)[a1->get()] > (*this)[a2->get()]; } 
  } alignment_score_map;
  
  
  bool sort_by_alignment_score (counted_ptr<bam_alignment>& a1, counted_ptr<bam_alignment>& a2) { return alignment_score_map[a1.get()] > alignment_score_map[a2.get()]; }  
  
  
  /*! Filter a list of alignments to only those that are eligible for mapping
   *  and return the BEST SCORE, which is currently the length of the aligned region of the read
   *  minus the number of base mismatches and minus the number of indels.
   *
   *  In junction_mode, negaive overlap of the junction sequence is also subtracted from the score.
   *    and returns all alignments above the minimum.
   *  Otherwise, returns only alignments tied for best.
   */
  
 uint32_t eligible_read_alignments(const Settings& settings, const cReferenceSequences& ref_seq_info, alignment_list& alignments, bool junction_mode, uint32_t min_match_score)
  {
    bool verbose = false;
    
    if (alignments.size() == 0) return 0;
    
    // require read to be mapped! -- @JEB maybe this should be checked sooner?
    for (alignment_list::iterator it = alignments.begin(); it != alignments.end();)
    {
      if (it->get()->unmapped())
        it = alignments.erase(it);
      else
        it++;
    }
    if (alignments.size() == 0) return 0;
    
    // require a minimum length of the read to be mapped
    for (alignment_list::iterator it = alignments.begin(); it != alignments.end();)
    {
      if ( !test_read_alignment_requirements(settings, ref_seq_info, *(it->get())) )
        it = alignments.erase(it);
      else
        it++;
    }
    if (alignments.size() == 0) return 0;
    
    // @JEB: We would ideally sort by alignment/mapping score here
    // the number of mismatches is our current proxy for this.
    
    uint32_t read_length = alignments.front()->read_length();
    alignment_score_map.clear();
    for (alignment_list::iterator it = alignments.begin(); it != alignments.end();)
    {  
      

      bam_alignment* ap = it->get(); // we are saving the pointer value as the map key
      
      
      uint32_t i;
      if (settings.bowtie2_align && settings.bowtie2) {
        i = ap->aux_get_i("AS"); // test using alignment score instead of mismatches.
      } else {
        i = alignment_mismatches(*ap, ref_seq_info);
      
        // @JEB may want to revisit this.
        // add in read only bases (negative overlap) as mismatches
        /*
        if (junction_mode) {
          JunctionInfo junction_info( ref_seq_info[ap->reference_target_id()].m_seq_id );
          if (junction_info.alignment_overlap < 0)
            // NOTE: subtraction here is adding b/c alignment_overlap < 0
            i -= junction_info.alignment_overlap;
        }
        */
        ASSERT(read_length >= i, "More mismatches than matches for read alignment. ");
        //ASSERT(read_length >= i, "More mismatches than matches for read alignment. " + ap->read_name());
         
         
        i = read_length - i;
      }
      
      // Only keep ones with a minimum score
      if (i < min_match_score) {
        it = alignments.erase(it);
      }
      else {
        ap->aux_set(kBreseqAlignmentScoreBAMTag, 'I', sizeof(uint32_t), (uint8_t*)&i);
        alignment_score_map[ap] = i;
        it++;
      }
      
    }
    
    alignments.sort(sort_by_alignment_score);
    
    if (verbose)
    {
      for (alignment_list::iterator it = alignments.begin(); it != alignments.end(); it++)
      {
        bam_alignment& a = *(it->get());
        cerr << a.query_start_1() << "-" << a.query_end_1() << " ";
        cerr << alignment_score_map[it->get()] << "\n";
      }
    }
        
    // how many reads share the best score?
    uint32_t last_best(0);
    uint32_t best_score = alignment_score_map[alignments.front().get()];
    
    // In non-junction mode we only return the best
    // In junction mode we return all above or equal to the minimum score (which is for the alignment to the reference genome).
    
    // no scores meet minimum
    for (alignment_list::iterator it = alignments.begin()++; it != alignments.end(); it++)
    {
      bam_alignment* ap = it->get();
      uint8_t is_best_score = (alignment_score_map[ap] == best_score) ? 1 : 0;
      
      ap->aux_set(kBreseqBestAlignmentScoreBAMTag, 'C', sizeof(uint8_t), (uint8_t*)&is_best_score);
      
      if (is_best_score)
        last_best++;
    }
    
    // Truncate only to ties for best if we are not in junction mode
    if (!junction_mode)
      alignments.resize(last_best);
    
    if (verbose)
    {
      cerr << last_best << endl;
      for (alignment_list::iterator it = alignments.begin(); it != alignments.end(); it++)
      {
        bam_alignment& a = *(it->get());
        cerr << a.query_start_1() << "-" << a.query_end_1() << endl;
      }
    }
    
    // Note that the score we return is higher for better matches...
    return best_score;
  }
  
  /*! Test the requirements for an alignment to be counted.
   *
   *  Returns whether the alignment passes.
   */
  bool test_read_alignment_requirements(const Settings& settings, const cReferenceSequences& ref_seq_info, const alignment_wrapper& a)
  {
    bool accept = true;
    
    if (a.unmapped()) return false;
    
    if (a.query_match_length() < settings.require_match_length)
      return false;
    
    if (a.query_match_length() < settings.require_match_fraction * static_cast<double>(a.read_length()) )
      return false;
    
    if (settings.maximum_read_mismatches >= 0)
    {
      int32_t mismatches = alignment_mismatches(a, ref_seq_info);
      if (mismatches > settings.maximum_read_mismatches)
        return false; 
    }
    
    return true;
  }

  // Splits unaligned reads out of a SAM file
  void PreprocessAlignments::split_matched_and_unmatched_alignments(
                                                                    uint32_t fastq_file_index,
                                                                    string fasta_file_name, 
                                                                    string input_sam_file_name, 
                                                                    string matched_sam_file_name, 
                                                                    string unmatched_fastq_file_name
                                                                    ) 
  {

    tam_file matched(matched_sam_file_name, fasta_file_name, ios::out);
    ofstream unmatched(unmatched_fastq_file_name.c_str());
    
    tam_file in(input_sam_file_name, fasta_file_name, ios::in);

    alignment_list al;
    while (in.read_alignments(al, false)) {
      if (al.front()->unmapped()) {
        unmatched << "@" << al.front()->read_name() << endl << al.front()->read_char_sequence() << endl << "+" << endl << al.front()->read_base_quality_char_string() << endl;
      } else {
        matched.write_alignments(fastq_file_index, al);
      }
    }
  }
  void PreprocessAlignments::split_matched_alignments(uint32_t fastq_file_index,
                                                      string fasta_file_name, 
                                                      string input_sam_file_name,
                                                      string matched_sam_file_name) 
  {
    tam_file in(input_sam_file_name, fasta_file_name, ios::in);
    tam_file matched(matched_sam_file_name, fasta_file_name, ios::out);

    alignment_list al;
    while (in.read_alignments(al, false)) {
      if (!al.front()->unmapped()) {
        matched.write_alignments(fastq_file_index, al);
      }
    }

    return;
  }
  
  void line_to_read_index(string s, int32_t& index, bool& mapped) {
    size_t start = s.find_first_of(":");
    size_t end = s.find_first_of(" \t", start);
    index = n(s.substr(start+1, end-(start+1)));
    
    size_t next = s.find_first_not_of(" \t", end);
    next = s.find_first_of(" \t", next);
    next = s.find_first_not_of(" \t", next);

    mapped = (s[next] != '*');
    
  }
  
  // Takes two SAM files and re-sorts read order so that it matches the original FASTQ file.
  // Requires reads to be renamed as we expect from 01_sequence_onversion: file_num/read_num.
  void PreprocessAlignments::merge_sort_sam_files(
                                                  string input_sam_file_name_1,
                                                  string input_sam_file_name_2,
                                                  string output_sam_file_name
                                                  )
  {
    
    ifstream input_sam_file_1(input_sam_file_name_1.c_str());
    ifstream input_sam_file_2(input_sam_file_name_2.c_str());
    
    ofstream out_sam_file(output_sam_file_name.c_str());
    
    string line_1, line_2;
    bool not_done_1 = getline(input_sam_file_1, line_1);
    while (not_done_1 && (line_1[0] == '@')) {
      not_done_1 = getline(input_sam_file_1, line_1);
    }
    
    bool not_done_2 = getline(input_sam_file_2, line_2);
    while (not_done_2 && (line_2[0] == '@')) {
      not_done_2 = getline(input_sam_file_2, line_2);
    }
         
    int32_t index_1 = -1; 
    int32_t index_2 = -1; 
    bool mapped_1 = false;
    bool mapped_2 = false;
    
    if (not_done_1) {
      line_to_read_index(line_1, index_1, mapped_1);
    }
    if (not_done_2) {
      line_to_read_index(line_2, index_2, mapped_2);
    }
    
    while (not_done_1 || not_done_2) {
      
      if (not_done_1 && (index_1 < index_2)) {
        
        if (mapped_1)
          out_sam_file << line_1 << endl;
        
        // read next line
        not_done_1 = getline(input_sam_file_1, line_1);
        if (not_done_1) {
          line_to_read_index(line_1, index_1, mapped_1);
        } else {
          index_1 = numeric_limits<int32_t>::max();
        }
      }
      else if (not_done_2) {
        if (mapped_2)
          out_sam_file << line_2 << endl;
        
        // read next line
        not_done_2 = getline(input_sam_file_2, line_2);
        if (not_done_2) {
          line_to_read_index(line_2, index_2, mapped_2);
        } else {
          index_2 = numeric_limits<int32_t>::max();
        }
      }
    }
  }
  
  /*
  // Takes two SAM files and re-sorts read order so that it matches the original FASTQ file.
  // Requires reads to be renamed as we expect from 01_sequence_onversion: file_num/read_num.
  void PreprocessAlignments::slow_merge_sort_sam_files(
                                                  uint32_t fastq_file_index,
                                                  string fasta_file_name,
                                                  string input_sam_file_name_1,
                                                  string input_sam_file_name_2,
                                                  string output_sam_file_name
                                                 )
  {
    
    tam_file input_sam_file_1(input_sam_file_name_1, fasta_file_name, ios::in);
    tam_file input_sam_file_2(input_sam_file_name_2, fasta_file_name, ios::in);
    
    tam_file out_sam_file(output_sam_file_name, fasta_file_name, ios::out);

    alignment_list alignment_list_1;
    alignment_list alignment_list_2;
    
    bool not_done_1 = input_sam_file_1.read_alignments(alignment_list_1, false);
    bool not_done_2 = input_sam_file_2.read_alignments(alignment_list_2, false);
    
    int32_t index_1 = -1;
    int32_t index_2 = -1;
    
    if (not_done_1) {
      index_1 = n(split(alignment_list_1.front()->read_name(),":")[1]);
    }
    if (not_done_2) {
      index_2 = n(split(alignment_list_2.front()->read_name(),":")[1]);

    }
    
    while (not_done_1 || not_done_2) {
      
      if (not_done_1 && (index_1 < index_2)) {
        
        if (!alignment_list_1.front()->unmapped())
          out_sam_file.write_alignments(fastq_file_index, alignment_list_1);

        // read next line
        not_done_1 = input_sam_file_1.read_alignments(alignment_list_1, false);
        if (not_done_1) {
          index_1 = n(split(alignment_list_1.front()->read_name(),":")[1]);
        } else {
          index_1 = numeric_limits<int32_t>::max();
        }
      }
      else if (not_done_2) {
        
        if (!alignment_list_2.front()->unmapped())
          out_sam_file.write_alignments(fastq_file_index, alignment_list_2);

        // read next line
        not_done_2 = input_sam_file_2.read_alignments(alignment_list_2, false);   
        if (not_done_2) {
          index_2 = n(split(alignment_list_2.front()->read_name(),":")[1]);
        } else {
          index_2 = numeric_limits<int32_t>::max();
        }
      }
    }
  }
  */
  
  /*! PreprocessAlignments::preprocess_alignments
   *
   *  Writes one SAM file of partial read alignments that could support junctions
   *  and another SAM file of the best read matches to the reference genome for
   *  a preliminary analysis of coverage that is necessary for junction prediction.
	 */
	void PreprocessAlignments::preprocess_alignments(const Settings& settings, Summary& summary, const cReferenceSequences& ref_seq_info)
	{
		cout << "Preprocessing alignments." << endl;
    
		// get the cutoff for splitting alignments with indels
		int32_t min_indel_split_len = settings.preprocess_junction_min_indel_split_length;
    
		// BSAM includes best matches as they are merged from all alignment files
		string preprocess_junction_best_sam_file_name = settings.preprocess_junction_best_sam_file_name;
    string reference_fasta_file_name = settings.reference_fasta_file_name;
    tam_file BSAM(preprocess_junction_best_sam_file_name, reference_fasta_file_name, ios_base::out);
    
    uint32_t i = 0;
		for (uint32_t index = 0; index < settings.read_files.size(); index++)
		{
			cReadFile read_file = settings.read_files[index];
			cerr << "  READ FILE::" << read_file.m_base_name << endl;
      
			string reference_sam_file_name = Settings::file_name(settings.reference_sam_file_name, "#", read_file.m_base_name);
      tam_file tam(reference_sam_file_name, reference_fasta_file_name, ios_base::in);
      
			// includes all matches, and splits long indels for this one read file
      string preprocess_junction_split_sam_file_name = Settings::file_name(settings.preprocess_junction_split_sam_file_name, "#", read_file.m_base_name);
      tam_file PSAM(preprocess_junction_split_sam_file_name, reference_fasta_file_name, ios_base::out);
      
			alignment_list alignments;
			while (tam.read_alignments(alignments, false))
			{
				if (++i % 100000 == 0)
					cerr << "    ALIGNED READ:" << i << endl;
        
        summary.preprocess_alignments.aligned_reads++;
        summary.preprocess_alignments.alignments += alignments.size();
        
				// for testing...
				if (settings.candidate_junction_read_limit != 0 && i > settings.candidate_junction_read_limit) break;
        
        if (alignments.front()->unmapped()) continue;
        
				// write split alignments
				if (min_indel_split_len != -1)
        {
					split_alignments_on_indels(settings, summary, PSAM, min_indel_split_len, alignments);
        }
        
				// write best alignments
        int32_t best_score = eligible_read_alignments(settings, ref_seq_info, alignments);
        BSAM.write_alignments(0, alignments, NULL);
      }
    }
    
    cout << "  Summary... " << endl
         << "  Aligned reads:                         " << setw(12) << right << summary.preprocess_alignments.aligned_reads << endl
         << "  Read alignments:                       " << setw(12) << right << summary.preprocess_alignments.alignments << endl
         << "  Alignments split on indels:            " << setw(12) << right << summary.preprocess_alignments.alignments_split_on_indels << endl
         << "  Reads with alignments split on indels: " << setw(12) << right << summary.preprocess_alignments.reads_with_alignments_split_on_indels << endl
         << "  Split alignments:                      " << setw(12) << right << summary.preprocess_alignments.split_alignments << endl
         << "  Reads with split alignments:           " << setw(12) << right << summary.preprocess_alignments.reads_with_split_alignments << endl
    ;
  }

  
  /*! Split alignments interrupted by indels into their segments and write to a SAM file. 
	 */
  
  void PreprocessAlignments::split_alignments_on_indels(const Settings& settings, Summary& summary, tam_file& PSAM, int32_t min_indel_split_len, const alignment_list& alignments)
  {
    //##
    //## @JEB: Note that this may affect the order of alignments in the SAM file. This has 
    //## consequences for the AlignmentCorrection step, which may assume highest scoring
    //## alignments are first in the TAM file.
    //##
    //## So only USE the split alignment file for creating a list of candidate junctions.
    //##
    
    assert(min_indel_split_len >= 0);
    
    // Keeps track of which to split
    alignment_list untouched_alignments;
    alignment_list split_alignments;
    
    untouched_alignments.read_base_quality_char_string = alignments.read_base_quality_char_string;
    untouched_alignments.read_base_quality_char_string_reversed = alignments.read_base_quality_char_string_reversed;
    uint32_t alignments_written = 0;      
    
    for(alignment_list::const_iterator it = alignments.begin(); it != alignments.end(); it++) {
      
      uint32_t* cigar_list = (*it)->cigar_array();
      bool do_split = false;
      
      for(uint32_t i=0; i<(*it)->cigar_array_length(); i++) {
        uint32_t op = cigar_list[i] & BAM_CIGAR_MASK;
        uint32_t len = cigar_list[i] >> BAM_CIGAR_SHIFT;
        if (((op == BAM_CINS) || (op == BAM_CDEL)) && (len >= static_cast<uint32_t>(min_indel_split_len))) {
          do_split = true;
        }
      }
      
      if (do_split) {
        split_alignments.push_back( *it );
        
      }
      else {
        // Check guard showing that the main alignment is good and thus we don't look at other alignments
        if ( (*it)->query_match_length() >= alignments.front()->read_length() - settings.required_both_unique_length_per_side)
          return;
        
        untouched_alignments.push_back( *it );
      }      
    }
        
    // Don't write if there is only one alignment to be written,
    // it takes at least two to make a candidate junction.

    for(alignment_list::iterator it = split_alignments.begin(); it != split_alignments.end(); ++it) {
      PSAM.write_split_alignment(min_indel_split_len, **it, alignments);
      alignments_written += 2;
    }
    
    if (untouched_alignments.size() + alignments_written > 1) {
      PSAM.write_alignments(0, untouched_alignments);
      alignments_written += untouched_alignments.size();
    }

    // record statistics
    if (split_alignments.size()>0) summary.preprocess_alignments.reads_with_alignments_split_on_indels++; 
    summary.preprocess_alignments.alignments_split_on_indels +=split_alignments.size() ;
    
    if (alignments_written > 0) summary.preprocess_alignments.reads_with_split_alignments++;
    summary.preprocess_alignments.split_alignments += alignments_written;
  }
    
	/*! Predicts candidate junctions
	 */
	void CandidateJunctions::identify_candidate_junctions(const Settings& settings, Summary& summary, const cReferenceSequences& ref_seq_info)
	{
		(void)summary; // TODO: save statistics
    bool verbose = false;
    
		// hash by junction sequence
		SequenceToKeyToJunctionCandidateMap candidate_junctions;
    
		// shortcut to summary data for this step
		Summary::CandidateJunctionSummaryData& hcs(summary.candidate_junction);
    
		uint32_t i = 0;
    
		for (uint32_t j = 0; j < settings.read_files.size(); j++)
		{
			cReadFile read_file = settings.read_files[j];
      
			string read_file_name = read_file.m_base_name;
			cerr << "  READ FILE::" << read_file_name << endl;

			// Decide which input SAM file we are using...
      
			string reference_sam_file_name = Settings::file_name(settings.preprocess_junction_split_sam_file_name, "#", settings.read_files[j].m_base_name);
      
      tam_file tam(reference_sam_file_name, settings.reference_fasta_file_name, ios_base::in);
			alignment_list alignments;
      
			while (tam.read_alignments(alignments, false))
			{
				if (alignments.size() == 0)
					break;
        
				if (++i % 10000 == 0)
					cerr << "    ALIGNED READ:" << i << " CANDIDATE JUNCTIONS:" << candidate_junctions.size() << endl;
        
				// for testing...
				if (settings.candidate_junction_read_limit != 0 && i > settings.candidate_junction_read_limit)
					break;
        
				alignments_to_candidate_junctions(settings, summary, ref_seq_info, candidate_junctions, alignments);
			}
    }

    ////
		// Merge all junctions with the same exact sequence 
    //   * They are hashed together for speed in this comparison
		////
    
    // New list consisting of merged junction candidates
    list<JunctionCandidatePtr> junction_candidate_list;
    
		for (SequenceToKeyToJunctionCandidateMap::iterator same_seq_it = candidate_junctions.begin(); same_seq_it != candidate_junctions.end(); same_seq_it++)
    {
      KeyToJunctionCandidateMap& kcjm = same_seq_it->second;
      
      KeyToJunctionCandidateMap::iterator it = kcjm.begin();
      
      // compare each pair 
      JunctionCandidatePtr one_junction = it->second;

      for (++it; it != kcjm.end(); ++it)
      {
        JunctionCandidatePtr *jcp1 = &one_junction;
        JunctionCandidatePtr *jcp2 = &it->second;

        merge_candidate_junctions(jcp1, jcp2);
        one_junction = *jcp1;
      }
      
      if (one_junction->pos_hash_score() >= settings.minimum_candidate_junction_pos_hash_score) junction_candidate_list.push_back(one_junction);
    }
    
    ////
		// Second round of merging, attempt to merge each item in the list with _every_ other one
    //   * This takes care of cases where one is a subsequence of another - we keep the shortest one
		////
    
    for (list<JunctionCandidatePtr>::iterator it1 = junction_candidate_list.begin(); it1 != junction_candidate_list.end();)
    {
      list<JunctionCandidatePtr>::iterator it2 = it1;
      bool deleted_it1 = false;
      
      for (++it2; it2 != junction_candidate_list.end(); )
      {
        bool deleted_it2 = false;
        JunctionCandidatePtr *jcp1 = &(*it1);
        JunctionCandidatePtr *jcp2 = &(*it2);
        
        if (verbose)
        {
          cout << "Comparing:" << endl;
          cout << (*jcp1)->junction_key() << endl;
          cout << (*jcp1)->sequence << endl;
          cout << (*jcp2)->junction_key() << endl;
          cout << (*jcp2)->sequence << endl;
        }
        
        
        bool merged = merge_candidate_junctions(jcp1, jcp2);
        if (merged)
        {
          if (verbose) cout << "MERGED" << endl;
          
          if (*jcp1 == *it1)
          {
            it2 = junction_candidate_list.erase(it2);
            deleted_it2 = true;
          }
          else
          {
            it1 = junction_candidate_list.erase(it1);
            deleted_it1 = true;
            break;
          }
        }
        
        if (!deleted_it2) it2++;
      }
      if (!deleted_it1) it1++;
    }

		////
		//  Combine hash into a list, retaining only one item (the best pos_hash score) for each unique sequence 
    //    (and also its reverse complement)
		////
    
		map<int32_t, int32_t> observed_pos_hash_score_distribution;
    
		vector<JunctionCandidate> combined_candidate_junctions;
    
		map<string, int32_t> handled_seq;
		map<string, JunctionCandidate> ids_to_print;
        
		// Map is pre-sorted by CandidateJunction::Sorter
		for (list<JunctionCandidatePtr>::iterator it = junction_candidate_list.begin(); it != junction_candidate_list.end(); ++it)
		{
      JunctionCandidatePtr& jcp = *it;
      JunctionCandidate& best_candidate_junction = *jcp;
      
			string junction_seq = best_candidate_junction.sequence;
      string rc_junction_seq = best_candidate_junction.reverse_complement_sequence;

      if (verbose) cout << "Handling " << junction_seq << endl;
      
			// These shouldn't be necessary checks, but keeping to detect unintended errors
      
      // We already handled the sequence or its reverse complement, skip.
      if (handled_seq.count(junction_seq)) continue;
      if (handled_seq.count(rc_junction_seq)) continue;

      string junction_id = best_candidate_junction.junction_key();
            
			handled_seq[junction_seq]++;
			handled_seq[rc_junction_seq]++;
            
			// Save the score in the distribution
			add_score_to_distribution(observed_pos_hash_score_distribution, best_candidate_junction.pos_hash_score());
      
			// Check minimum requirements
			if (best_candidate_junction.pos_hash_score() < settings.minimum_candidate_junction_pos_hash_score)
				continue;
      
      if (verbose) cout << "  best candidate junction:" << endl << best_candidate_junction.junction_key() << endl;
      
			// Make sure it isn't a duplicate junction id -- this should NEVER happen and causes downstream problem.
			// ---> Begin sanity check
			if (ids_to_print.count(junction_id))
			{
        JunctionCandidate& ccj = ids_to_print[junction_id];
        
				cout << "Attempt to create junction candidate with duplicate id: " << junction_id << endl;
        
				cout << "==Existing junction==" << endl;      
        cout << "  id: " << ccj.junction_key() << endl;
        cout << "  pos_hash_score: " << ccj.pos_hash_score() << endl;
        cout << "  pos_hash_score: " << ccj.pos_hash_score() << endl;
        cout << "  seq: " << ccj.sequence << endl;
        cout << "  rc_seq: " << ccj.reverse_complement_sequence << endl;   
        
				cout << "==New junction==" << endl;
        cout << "  id: " << best_candidate_junction.junction_key() << endl;
        cout << "  pos_hash_score: " << best_candidate_junction.pos_hash_score() << endl;
        cout << "  pos_hash_score: " << best_candidate_junction.pos_hash_score() << endl;
        cout << "  seq: " << best_candidate_junction.sequence << endl;
        cout << "  rc_seq: " << best_candidate_junction.reverse_complement_sequence << endl;  
        
				assert (best_candidate_junction.sequence == ids_to_print[junction_id].sequence);
				exit(-1);
			}
      // <--- End sanity check

			ids_to_print[best_candidate_junction.junction_key()] = best_candidate_junction;
      
			combined_candidate_junctions.push_back(best_candidate_junction);
		}
    
		sort(combined_candidate_junctions.begin(), combined_candidate_junctions.end(), JunctionCandidate::sort_by_scores_and_seq_length);
    
    if (verbose)
    {
      for (vector<JunctionCandidate>::iterator it=combined_candidate_junctions.begin(); it < combined_candidate_junctions.end(); it++ )
      {
        cout << "ID: " << it->junction_key() << endl;
        cout << "  pos_hash_score: " << it->pos_hash_score() << endl;
        cout << "  seq: " << it->sequence << endl;
        cout << "  rc_seq: " << it->reverse_complement_sequence << endl;
      }
    }
    
		///
		// Limit the number of candidate junctions that we print by:
		//   (1) A maximum number of candidate junctions
		//   (2) A maximum length of the sequences in candidate junctions
    //   (3) But take at least some minimum despite these
		///
    
		cerr << "  Taking top candidate junctions..." << endl;
    
		// adding up the lengths might be too time-consuming to be worth it...
		int32_t total_cumulative_cj_length = 0;
		int32_t total_candidate_junction_number = combined_candidate_junctions.size();
		for (uint32_t j = 0; j < combined_candidate_junctions.size(); j++)
			total_cumulative_cj_length += combined_candidate_junctions[j].sequence.size();
    
		uint32_t cumulative_cj_length = 0;
		int32_t lowest_accepted_pos_hash_score = 0;
    
		// Right now we limit the candidate junctions to have a length no longer than the reference sequence times some factor.
		uint32_t cj_length_limit = static_cast<uint32_t>(summary.sequence_conversion.total_reference_sequence_length * settings.maximum_candidate_junction_length_factor);
		uint32_t maximum_candidate_junctions = settings.maximum_candidate_junctions;
		uint32_t minimum_candidate_junctions = settings.minimum_candidate_junctions;
    
		fprintf(stderr, "  Minimum number to keep: %7d \n", minimum_candidate_junctions);
		fprintf(stderr, "  Maximum number to keep: %7d \n", maximum_candidate_junctions);
		fprintf(stderr, "  Maximum length to keep: %7d bases\n", cj_length_limit);
    
		cerr << "    Initial: Number = " << total_candidate_junction_number << ", Cumulative Length = " << total_cumulative_cj_length << " bases" << endl;
    
		if (combined_candidate_junctions.size() > 0)
		{
			vector<JunctionCandidate> remaining_ids;
			vector<JunctionCandidate> list_in_waiting;
			int32_t add_cj_length = 0;
			int32_t num_duplicates = 0;
      
			i = 0;
			uint32_t current_pos_hash_score = combined_candidate_junctions[i].pos_hash_score();
      
			// Check to make sure that adding everything from the last iteration doesn't put us over any limits...
			uint32_t new_number = remaining_ids.size() + list_in_waiting.size();
			uint32_t new_length = cumulative_cj_length + add_cj_length;
			while (	( new_number <= minimum_candidate_junctions ) || ((new_length <= cj_length_limit) && (new_number <= maximum_candidate_junctions)) )
			{
				// OK, add everything from the last iteration
				cumulative_cj_length += add_cj_length;
				remaining_ids.reserve(remaining_ids.size() + list_in_waiting.size());
				remaining_ids.insert(remaining_ids.end(), list_in_waiting.begin(), list_in_waiting.end());
        
				lowest_accepted_pos_hash_score = current_pos_hash_score;
        
				// Zero out what we will add
				add_cj_length = 0;
				list_in_waiting.clear();
				num_duplicates = 0;
        
				// Check to make sure we haven't exhausted the list
				if (i >= combined_candidate_junctions.size()) break;
        
        // Grab the next chunk with the same score
				current_pos_hash_score = combined_candidate_junctions[i].pos_hash_score();
				while (
               (i < combined_candidate_junctions.size())
               && (combined_candidate_junctions[i].pos_hash_score() == current_pos_hash_score)
               )
				{
					JunctionCandidate c = combined_candidate_junctions[i];
					list_in_waiting.push_back(c);
					add_cj_length += c.sequence.size();
					i++;
				}
        
				new_number = remaining_ids.size() + list_in_waiting.size();
				new_length = cumulative_cj_length + add_cj_length;
        
				fprintf(stderr, "      Testing Pos Hash Score = %4d, Num = %6d, Length = %6d\n", current_pos_hash_score, int32_t(list_in_waiting.size()), add_cj_length);
			}
			combined_candidate_junctions = remaining_ids;
		}
    
		int32_t accepted_candidate_junction_number = combined_candidate_junctions.size();
		cerr << "    Accepted: Number = " << accepted_candidate_junction_number << ", Pos Hash Score >= " << lowest_accepted_pos_hash_score << ", Cumulative Length = " << cumulative_cj_length << " bases" << endl;
    
		// Save summary statistics
		hcs.total.number = total_candidate_junction_number;
		hcs.total.length = total_cumulative_cj_length;
    
		hcs.accepted.number = accepted_candidate_junction_number;
		hcs.accepted.length = cumulative_cj_length;
		hcs.accepted.pos_hash_score_cutoff = lowest_accepted_pos_hash_score;
		hcs.pos_hash_score_distribution = observed_pos_hash_score_distribution;
    
		///
		// Print out the candidate junctions, sorted by the lower coordinate, higher coord, then number
		///
    
		sort(combined_candidate_junctions.begin(), combined_candidate_junctions.end(), JunctionCandidate::sort_by_ref_seq_coord);
    
    cFastaFile out(settings.candidate_junction_fasta_file_name, ios_base::out);
    
		for (uint32_t j = 0; j < combined_candidate_junctions.size(); j++)
		{
			JunctionCandidate junction = combined_candidate_junctions[j];
      cFastaSequence seq; //= { junction.junction_key(), "", junction.sequence };
      seq.m_name = junction.junction_key();
      seq.m_description = "";
      seq.m_sequence = junction.sequence;
			out.write_sequence(seq);
		}
		out.close();
    
		summary.candidate_junction = hcs;
	}
  
  
  bool CandidateJunctions::merge_candidate_junctions(JunctionCandidatePtr*& jcp1, JunctionCandidatePtr*& jcp2)
  {
    bool verbose = false;
    uint32_t merged = 0; // 0 for not merged, 1 for sequence, 2 for reverse complement of sequence
    
    JunctionCandidate& jc1 = **jcp1;
    JunctionCandidate& jc2 = **jcp2;
    
    // Determine whether one is a subsequence of the other (including on the opposite strand)
    if (jc1.sequence.size() > jc2.sequence.size())
    {
      if (jc1.sequence.find(jc2.sequence) != string::npos)
        merged = 1;
      else if (jc1.reverse_complement_sequence.find(jc2.sequence) != string::npos )
        merged = 2;
    }
    else
    {
      if (jc2.sequence.find(jc1.sequence) != string::npos)
        merged = 1;
      else if (jc2.sequence.find(jc1.reverse_complement_sequence) != string::npos)
        merged = 2;
    }
    
    if (!merged) return (merged != 0);
    
    if (verbose) cout << ((merged == 1) ? "Merged same strand" : "Merged reverse complement") << endl;
    
    JunctionCandidatePtr* merge_into_p = NULL;
    JunctionCandidatePtr* merge_from_p = NULL;
    
    // this is a rather complicated compare function to favor longer sequences and those with close coords
  
    // This comparison of junctions favors...
    //   1) the longer sequence 
    //   2) both sides on the same reference sequence 
    //   3) the smallest coordinate on side 1
    
    if ( jc2 < jc1 )
    {
      merge_into_p = jcp1;
      merge_from_p = jcp2; 
    }
    else
    {
      merge_into_p = jcp2;
      merge_from_p = jcp1;
    }
    
    JunctionCandidate& merge_into = **merge_into_p;
    JunctionCandidate& merge_from = **merge_from_p;
    
    
    if (verbose) cout << "Merging into:" << merge_into.junction_key() << endl;
    if (verbose) cout << merge_into.sequence << endl;
    if (verbose) cout << "Merging from:" << merge_from.junction_key() << endl;
    if (verbose) cout << merge_from.sequence << endl;
    
    // Carry over redundancy that was already assigned (based on different coord matching, but sequence being the same)
    if (merged == 2)
    {
      merge_into.sides[0].redundant = merge_from.sides[1].redundant || merge_into.sides[0].redundant;
      merge_into.sides[1].redundant = merge_from.sides[0].redundant || merge_into.sides[1].redundant;
    }
    else // merged == 1
    {
      merge_into.sides[0].redundant = merge_from.sides[0].redundant || merge_into.sides[0].redundant;
      merge_into.sides[1].redundant = merge_from.sides[1].redundant || merge_into.sides[1].redundant;
    }
    
    // If one of the sides is identical (in terms of reference coordinate and strand), then mark it as redundant
    if (verbose) cout << "Merging into carryover:" << merge_into.junction_key() << endl;
    for (uint32_t into_side = 0; into_side < 2; into_side++) 
    {
      bool found = false;
      
      for (uint32_t from_side = 0; from_side < 2; from_side++) 
      {
        if (merge_into.sides[into_side] == merge_from.sides[from_side])
          found = true;
      }
      
      if (!found)
      {
        merge_into.sides[into_side].redundant = true;
        if (verbose) cout << "Marking side " << into_side << " as redundant." << endl;
      }
      
    }
    
    jcp1 = merge_into_p;
    jcp2 = merge_from_p;

    return (merged != 0);
  }


	bool CandidateJunctions::alignment_pair_to_candidate_junction(
                                                            const Settings& settings, 
                                                            Summary& summary, 
                                                            const cReferenceSequences& ref_seq_info, 
                                                            AlignmentPair& ap,
                                                            JunctionCandidatePtr& returned_junction_candidate
                                                            )
	{    
		bool verbose = false;
    
    //if ( (ap.a1.read_name() == "1:907874") || (ap.a1.read_name() == "1:15230"))
    //   verbose = true;
    
    
    // clear the return value
    returned_junction_candidate = JunctionCandidatePtr(NULL);
    

		// set up local settings
		int32_t flanking_length = summary.sequence_conversion.max_read_length;

    bam_alignment& a1 = ap.a1;
    bam_alignment& a2 = ap.a2;
    
		// Method
		//
		// Hash junctions by a key showing the inner coordinate of the read.
		// and the direction propagating from that position in the reference sequence.
		// Prefer the lower coordinate side of the read for main hash.
		//
		// REL606__1__1__REL606__4629812__-1__0
		// means the junction sequence is 36-1 + 4629812-4629777 from the reference sequence
		//
		// On the LEFT side: -1 means this is highest coord of alignment, junction seq begins at lower coord
		//                    1 means this is lowest coord of alignment, junction seq begins at higher coord
		// On the RIGHT side:-1 means this is highest coord of alignment, junction seq continues to lower coord
		//                    1 means this is lowest coord of alignment, junction seq continues to higher coord

		int32_t i = 0;

		string read_id = a1.read_name();
    
		// First, sort matches by their order in the query
		bam_alignment& q1 = a1;
		bam_alignment& q2 = a2;
		uint32_t q1_start, q1_end;
		uint32_t q2_start, q2_end;
		q1.query_stranded_bounds_1(q1_start, q1_end);
		q2.query_stranded_bounds_1(q2_start, q2_end);

    if (verbose)
      cout << q1.read_name() << endl;
		
    if (verbose)
			cout << q1_start << ", " << q1_end << ", " << q2_start << ", " << q2_end << endl;

		// Reverse the coordinates to be consistently such that 1 refers to lowest...
		if (q2_start < q1_start)
		{
			swap(q1_start, q2_start);
			swap(q1_end, q2_end);
			swap(q1, q2);
		}
    
		// create hash key and store information about the location of this hit
		bool hash_strand_1 = q1.reversed();
		string hash_seq_id_1 = ref_seq_info[q1.reference_target_id()].m_seq_id;
		const string& ref_seq_1(ref_seq_info[q1.reference_target_id()].m_fasta_sequence.m_sequence);

		bool hash_strand_2 = !q2.reversed();
		string hash_seq_id_2 = ref_seq_info[q2.reference_target_id()].m_seq_id;
		const string& ref_seq_2(ref_seq_info[q2.reference_target_id()].m_fasta_sequence.m_sequence);
    
		// how much overlap is there between the two matches?
		// positive if the middle sequence can match either side of the read
		// negative if there is sequence in the read NOT matched on either side
		int32_t overlap = -1 * (q2_start - q1_end - 1);
    
		//
		// OVERLAP MISMATCH CORRECTION
		//
		// If there are mismatches in one or the other read in the overlap region
		// then we need to adjust the coordinates. Why? All sequences that we
		// retrieve are from the reference sequence and there are two choices
		// for where to extract this non-necessarily identical sequence!

		// save these as variables, because we may have to adjust them
		int32_t r1_start = a1.reference_start_1();
		int32_t r1_end = a1.reference_end_1();
		int32_t r2_start = a2.reference_start_1();
		int32_t r2_end = a2.reference_end_1();

		if (verbose)
		{
			string ref_seq_matched_1 = ref_seq_1.substr(r1_start - 1, r1_end - r1_start + 1);
			string ref_seq_matched_2 = ref_seq_2.substr(r2_start - 1, r2_end - r2_start + 1);

			cout << "==============> Initial Matches" << endl;
			cout << "Alignment #1" << endl;
			cout << "qpos: " << q1_start << "-" << q1_end << " rpos: " << r1_start << "-" << r1_end << " reversed: " << q1.reversed() << endl;
			cout << q1.read_char_sequence() << endl << ref_seq_matched_1 << endl;

			cout << "Alignment #2" << endl;
			cout << "qpos: " << q2_start << "-" << q2_end << " rpos: " << r2_start << "-" << r2_end << " reversed: " << q2.reversed() << endl;
			cout << q2.read_char_sequence() << endl << ref_seq_matched_2 << endl;
			cout << "<==============" << endl;

			// debug print information
			cout << "=== overlap: " << overlap << endl;
		}

		// Adjust the overlap in cases where there is a mismatch within the overlap region
    int32_t overlap_in_reference = 0;
    
    if (q1.strand() == q2.strand()) {
      
      if ((r2_start >=  r1_start) && (r2_start <=  r1_end)) {
        if (q1.strand() == +1)
          overlap_in_reference = r1_end - r2_start + 1;
        else
          overlap_in_reference = r2_end - r1_start + 1;      
      }
    }

    if (verbose)
      cout << "=== overlap in reference: " << overlap_in_reference << endl;
    
		if ((overlap > 0) || (overlap_in_reference > 0))
		{
      overlap = max(overlap, overlap_in_reference);
      
      if (verbose)
        cout << "=== overlap: " << overlap << endl;
      
			int32_t q1_move, q2_move, r1_move, r2_move;
			q1.num_matches_from_end(ref_seq_info, false, overlap, q1_move, r1_move);
			q2.num_matches_from_end(ref_seq_info, true, overlap, q2_move, r2_move);

			if (q1_move >= 0)
			{
				if (verbose)
					cout << "ALIGNMENT #1 OVERLAP MISMATCH: " << q1_move << ", " << r1_move << endl;
				// change where it ENDS
				q1_end -= q1_move;
				if (q1.reversed())
					r1_start += r1_move;
				else
					r1_end -= r1_move;
			}

			if (q2_move >= 0)
			{
				if (verbose)
					cout << "ALIGNMENT #2 OVERLAP MISMATCH: " << q2_move << ", " << r2_move << endl;
				// change where it STARTS
				q2_start += q2_move;
				if (!q2.reversed())
					r2_start += r2_move;
				else
					r2_end -= r2_move;
			}

			// We need to re-check that they have the required amount of unique length if they moved
			if (q1_move >= 0 || q2_move >= 0)
			{
				if (verbose)
					cout << "Rejecting alignment because there is not enough overlap" << endl;
				int32_t dont_care;
        
        AlignmentPair ap(q1, q2, settings);
				if (!ap.pass) return false;
			}

			//re-calculate the overlap
			overlap = -1 * (q2_start - q1_end - 1);
			if (verbose)
				cout << "=== overlap corrected for mismatches " << overlap << endl;
		}
        
    //Recalculate the overlap in the reference
    if (q1.strand() == q2.strand()) {
      
      if ((r2_start >=  r1_start) && (r2_start <=  r1_end)) {
        if (q1.strand() == +1)
          overlap_in_reference = r1_end - r2_start + 1;
        else
          overlap_in_reference = r2_end - r1_start + 1;      
      }
    }
    
    // create hash coords AFTER overlap adjustment
		int32_t hash_coord_1 = (hash_strand_1) ? r1_start : r1_end;
		int32_t hash_coord_2 = (hash_strand_2) ? r2_start : r2_end;

    // Further correction for zero overlap.
    //
    // If there are multiple ways the two sides could have been aligned...shift them
    // over so as much is included in the lower reference coordinate side as possible
    // 
    // This case can arise for reads matching the same reference bases in SSAHA2
    // or after correcting for mismatches in the overlap (?)
    
    if (overlap == 0)
    {
      
      int32_t lower_coord_side = (hash_coord_1 < hash_coord_2) ? -1 : +1;
      int32_t move_r1_pos = (hash_strand_1) ? -1 : +1;
      int32_t move_r2_pos = (!hash_strand_2) ? -1 : +1;
      move_r1_pos *= lower_coord_side;
      move_r2_pos *= lower_coord_side;

      if (verbose) cout << "ZERO OVERLAP SHIFT TEST" << endl;
      if (verbose) cout << "hash coord 1:" << hash_coord_1 << " hash coord 2: " << hash_coord_2 << endl;
      uint32_t test_r1_pos = (hash_strand_1) ? r1_start : r1_end;
      uint32_t test_r2_pos = (hash_strand_2) ? r2_start : r2_end;
      
      if (lower_coord_side == -1)
        test_r2_pos += move_r2_pos;
      else
        test_r1_pos += move_r1_pos;
      
      if ( (test_r1_pos >= 1) && (test_r1_pos <= ref_seq_1.size() )
        && (test_r2_pos >= 1) && (test_r2_pos <= ref_seq_2.size() ) )
      {
        string test_r1_char;
        string test_r2_char;

        test_r1_char = ref_seq_1.substr(test_r1_pos - 1, 1);
        if (hash_strand_1) test_r1_char = reverse_complement(test_r1_char);
        test_r2_char = ref_seq_2.substr(test_r2_pos - 1, 1);
        if (!hash_strand_2) test_r2_char = reverse_complement(test_r2_char);
        
        while (test_r1_char == test_r2_char)
        {          
          test_r1_pos += move_r1_pos;
          test_r2_pos += move_r2_pos;
          
          if (! (
                 (test_r1_pos >= 1) && (test_r1_pos <= ref_seq_1.size())
                 && (test_r2_pos >= 1) && (test_r2_pos <= ref_seq_2.size()) 
                 ) )
          {
            test_r1_pos -= move_r1_pos;
            test_r2_pos -= move_r2_pos;
            break;
          }

          
          test_r1_char = ref_seq_1.substr(test_r1_pos - 1, 1);
          if (hash_strand_1) test_r1_char = reverse_complement(test_r1_char);
          test_r2_char = ref_seq_2.substr(test_r2_pos - 1, 1);
          if (!hash_strand_2) test_r2_char = reverse_complement(test_r2_char);
        }
        
        // backtrack by one
        if (lower_coord_side == -1)
          test_r2_pos -= move_r2_pos;
        else
          test_r1_pos -= move_r1_pos;
        
        hash_coord_1 = test_r1_pos;
        hash_coord_2 = test_r2_pos;
        if (verbose) cout << "hash coord 1:" << hash_coord_1 << " hash coord 2: " << hash_coord_2 << endl;
        
      }
    }

		// these are the positions of the beginning and end of the read, across the junction
		// query 1 is the start of the read, which is why we hash by this coordinate
		// (it is less likely to be shifted by a nucleotide or two by base errors)
		int32_t read_begin_coord = (hash_strand_1) ? r1_end : r1_start;

		if (verbose)
		{
			string ref_seq_matched_1 = ref_seq_1.substr(r1_start - 1, r1_end - r1_start + 1);
			string ref_seq_matched_2 = ref_seq_2.substr(r2_start - 1, r2_end - r2_start + 1);

			cout << "==============> Initial Matches" << endl;
			cout << "Alignment #1" << endl;
			cout << "qpos: " << q1_start << "-" << q1_end << " rpos: " << r1_start << "-" << r1_end << " reversed: " << q1.reversed() << endl;
			cout << q1.read_char_sequence() << endl << ref_seq_matched_1 << endl;
//			print Dumper($q1->cigar_array);

			cout << "Alignment #2" << endl;
			cout << "qpos: " << q2_start << "-" << q2_end << " rpos: " << r2_start << "-" << r2_end << " reversed: " << q2.reversed() << endl;
			cout << q2.read_char_sequence() << endl << ref_seq_matched_2 << endl;
//			print Dumper($q2->cigar_array);
			cout << "<==============" << endl;
		}

		// Calculate an offset that only applies if the overlap is positive (sequence is shared between the two ends)
		int32_t overlap_offset = (overlap > 0) ? overlap : 0;
		if (verbose)
			cout << "Overlap offset: " << overlap_offset << endl;

    ////
		// Create the sequence of the candidate junction
    ////
		string junction_seq_string = "";

		// first end - contains the overlap for >0 overlap
		int32_t flanking_left = flanking_length;
		if (hash_strand_1 != 1) // alignment is not reversed
		{
			// start_pos is in 1-based coordinates
			int32_t start_pos = hash_coord_1 - (flanking_left - 1) - overlap_offset;
			if (start_pos < 1)
			{
				if (verbose)
					cout << "START POS 1: " << start_pos << " < 0" << endl;
				flanking_left += start_pos - 1;
				start_pos = 1;
			}
			string add_seq = ref_seq_1.substr(start_pos - 1, flanking_left + overlap_offset);
			if (verbose) cout << "1F: " << add_seq << endl;
			junction_seq_string += add_seq;
		}
		else
		{
			// end_pos is in 1-based coordinates
			uint32_t end_pos = hash_coord_1 + (flanking_left - 1) + overlap_offset;
			if (end_pos > ref_seq_1.size())
			{
				if (verbose) cout << "END POS 1: (" << end_pos << " < length" << endl;
				flanking_left -= end_pos - ref_seq_1.size();
				end_pos = ref_seq_1.size();
			}

			string add_seq = ref_seq_1.substr(end_pos - (flanking_left + overlap_offset), flanking_left + overlap_offset);
			add_seq = reverse_complement(add_seq);
			if (verbose) cout << "1R: " << add_seq << endl;
			junction_seq_string += add_seq;
    }

		// Add any unique junction sequence that was only in the read
		// and NOT present in the reference genome (overlap < 0)
		string unique_read_seq_string = "";
		if (overlap < 0)
			unique_read_seq_string = q1.read_char_sequence().substr(q1_end, -1 * overlap);
		junction_seq_string += unique_read_seq_string;

		if (verbose) cout << "+U: " << unique_read_seq_string << endl;

		// second end - added without overlapping sequence
		int32_t flanking_right = flanking_length;
    
		if (hash_strand_2 == 1) //alignment is not reversed
		{
			// end_pos is in 1-based coordinates
			uint32_t end_pos = hash_coord_2 + (flanking_right - 1) + overlap_offset;
			if (end_pos > ref_seq_2.size())
			{
				if (verbose)
					cout << "END POS 2: (" << end_pos << " < length" << endl;
				flanking_right -= (end_pos - ref_seq_2.size());
				end_pos = ref_seq_2.size();
			}
			string add_seq = ref_seq_2.substr(end_pos - flanking_right, flanking_right);
			if (verbose) cout << "2F: " << add_seq << endl;
			junction_seq_string += add_seq;
		}
		else // alignment is reversed
		{
			// start_pos is in 1-based coordinates
			int32_t start_pos = hash_coord_2 - (flanking_right - 1) - overlap_offset;
			if (start_pos < 1)
			{
				if (verbose) cout << "START POS 2: " << start_pos << " < 0" << endl;
				flanking_right += start_pos - 1;
				start_pos = 1;
			}
			string add_seq = ref_seq_2.substr(start_pos - 1, flanking_right);			
      add_seq = reverse_complement(add_seq);
			if (verbose) cout << "2R: " << add_seq << endl;
			junction_seq_string += add_seq;
		}

		if (verbose) cout << "3: " << junction_seq_string << endl;

		// create hash coords after adjustment for overlap
		if (!hash_strand_1)
			r1_end -= overlap_offset;
		else //reversed
			r1_start += overlap_offset;

		if (!hash_strand_2)
			r2_end -= overlap_offset;
		else //reversed
			r2_start += overlap_offset;

    ///
    //  Important: Here is where we choose the strand of the junction sequence
    ///
    
		// want to be sure that lowest ref coord is always first for consistency
		if ( hash_seq_id_1.compare(hash_seq_id_2) > 0 || ((hash_seq_id_1.compare(hash_seq_id_2) == 0) && (hash_coord_2 < hash_coord_1)) )
		{
			swap(hash_coord_1, hash_coord_2);
			swap(hash_strand_1, hash_strand_2);
			swap(hash_seq_id_1, hash_seq_id_2);
			swap(flanking_left, flanking_right);

			junction_seq_string = reverse_complement(junction_seq_string);
			unique_read_seq_string = reverse_complement(unique_read_seq_string);
		}
    
    // create junction candidate pointer
		JunctionCandidate* candidate_junction_ptr = new JunctionCandidate(
          JunctionInfo(
                       JunctionSide(hash_seq_id_1,	hash_coord_1,	hash_strand_1 ? +1 : -1), // note conversion of strand 0/1 to -1/+1
                       JunctionSide(hash_seq_id_2,	hash_coord_2,	hash_strand_2 ? +1 : -1), // note conversion of strand 0/1 to -1/+1
                       overlap, 				
                       unique_read_seq_string,
                       flanking_left,
                       flanking_right
                       ), 
          junction_seq_string
    );
    candidate_junction_ptr->read_begin_hash[read_begin_coord]++;

    // set the return value (which takes control of the allocated pointer)
    returned_junction_candidate = JunctionCandidatePtr(candidate_junction_ptr);

    
    // Bowtie2 specific code: (Not needed for SSAHA2) -->
    // need to rule out short indels that we want to fit by RA methods.
    // would be better to get rid of these at an earlier stage.
    // This code may not get rid of all possible SUB conditions.

    if ( hash_seq_id_1.compare(hash_seq_id_2) == 0 ) {
      
      // Insertion of the same base
      //    ------>
      // ----> 
      if ((overlap > 0) && (hash_strand_1 ==1) && (hash_strand_2 != 1) && (hash_coord_2 - hash_coord_1 <= overlap) && ( hash_coord_1 + overlap - hash_coord_2 + 1  < settings.preprocess_junction_min_indel_split_length)) {
        //cerr << "Rejected/Insertion of same base(s):" << candidate_junction_ptr->junction_key() << endl;
        return false;
      }
      
      // Insertion of unique bases at an existing junction
      if ((overlap < 0) && (hash_strand_1 !=1) && (hash_strand_2 == 1) && ( hash_coord_1 + 1 == hash_coord_2) && (-overlap < settings.preprocess_junction_min_indel_split_length)) {
        //cerr << "Rejected/Insertion of unique base(s):" << candidate_junction_ptr->junction_key() << endl;
        return false;
      }
    }
    
    //
    // <-- End Bowtie2 specific code
    //
    
		if (verbose)
		{
      string junction_id = candidate_junction_ptr->junction_key();
			cout << "READ ID: " << a1.read_name() << endl;
			cout << "JUNCTION ID: " << junction_id << endl;
		}

		ASSERT(junction_seq_string.size() > 0, "Junction sequence not found."); 
		ASSERT(junction_seq_string.size() == flanking_left + flanking_right + static_cast<uint32_t>(abs(overlap)),
           "Incorrect junction sequence length for " + candidate_junction_ptr->junction_key() + "\n" + junction_seq_string);
    
		return true;
	}

	void CandidateJunctions::alignments_to_candidate_junctions(
                                                             const Settings& settings, 
                                                             Summary& summary, const 
                                                             cReferenceSequences& ref_seq_info, 
                                                             SequenceToKeyToJunctionCandidateMap& candidate_junctions, 
                                                             alignment_list& alignments
                                                             )
	{
    (void)summary;
		bool verbose = false;

		if (verbose)
		{
			cout << endl << "###########################" << endl;
			cout << alignments.front()->read_name();
			cout << endl << "###########################" << endl;
		}

		// Must  have multiple matches to support a new junction.
		if (alignments.size() <= 1)
			return;

    ////
    // Split the read alignments into two lists
    //  (1) list1 contains all matches starting at the beginning of the read
    //  (2) list2 contains all other matches
    // Accepted pairs must have a member from each list
    // Requiring (1) saves a number of comparisons and gets rid of a lot of bad matches.
    ////

    alignment_list list1, list2;
    
		if (verbose)
		{
			cout << alignments.front()->read_name() << endl;
			cout << "Total matches: " << alignments.size() << endl;
		}

    const uint32_t unmatched_end_min_coord = settings.required_junction_read_end_min_coordinate(alignments.front()->read_length());
    
    for (alignment_list::iterator it=alignments.begin(); it != alignments.end(); it++)
		{
			bam_alignment_ptr a = *it;
      
			uint32_t a_start, a_end;
      a->query_stranded_bounds_1(a_start, a_end);

			if (verbose) cout << "(" << a_start << ", " << a_end << ")" << endl;
          
			if (a_start == 1)
				list1.push_back(a);
			else if (a_end >= unmatched_end_min_coord)
				list2.push_back(a);
		}
    
		// The first match in this category is the longest
		if (verbose)
		{
			cout << "  List1: " << list1.size() << endl;
			cout << "  List2: " << list2.size() << endl;
		}


    ////
		// Try adding together each pair of matches to make a junction, by looking at read coordinates
    //   Only keep pairs that cover the maximum number of query bases encountered
    ////
    
    int32_t max_union_length = 0;
		vector<AlignmentPair> passed_pair_list;

    for (alignment_list::iterator it1 = list1.begin(); it1 != list1.end(); it1++)
		{
			bam_alignment_ptr& a1 = *it1;

      for (alignment_list::iterator it2 = list2.begin(); it2 != list2.end(); it2++)
			{
				bam_alignment_ptr& a2 = *it2;

				// constructing the alignment pair calculates statistics about their overlap
        // and tests the guards that are in Settings.
        
        AlignmentPair ap(*a1, *a2, settings);

				if (ap.pass)
        {
          // right now we use the union length as a kind of score
          // matches with greater union length trump those with smaller
          if (ap.union_length > max_union_length)
          {
            max_union_length = ap.union_length;
            passed_pair_list.clear();
          }

          passed_pair_list.push_back(ap);
        }
      }
		}
    
    // see if the junction sequence is unique (not contained in or containing any other sequences)    
		for (uint32_t i = 0; i < passed_pair_list.size(); i++)
		{
			AlignmentPair& ap = passed_pair_list[i];
      
      JunctionCandidatePtr new_junction_ptr;
			bool passed = alignment_pair_to_candidate_junction(settings, summary, ref_seq_info, ap, new_junction_ptr);
			if (!passed) continue;
      
      JunctionCandidate& new_junction = *new_junction_ptr;
      if (verbose) cout << "Testing junction: " << new_junction_ptr->junction_key() << endl << new_junction_ptr->sequence << endl;
      
			string junction_id = new_junction.junction_key();
      if (verbose) cout << junction_id << endl;

      if ((candidate_junctions.count(new_junction.sequence) == 0) || (candidate_junctions[new_junction.sequence].count(junction_id) == 0)) {
        // this is a new combination of sequence and id
        if (verbose) cout << "New saved junction " << new_junction.sequence << " " << junction_id << endl;
        candidate_junctions[new_junction.sequence][junction_id] = new_junction_ptr;
      }
      else
      {
        // update score of existing junction
        JunctionCandidate& cj = *candidate_junctions[new_junction.sequence][junction_id];
        cj.read_begin_hash[new_junction.read_begin_hash.begin()->first]++;
        if (verbose) cout << "Updating score of existing " << new_junction.sequence << " " << junction_id << endl 
          << "Pos: " << new_junction.read_begin_hash.begin()->first << " Score: " << cj.pos_hash_score()  << endl;
      }

    } // end passed pair list
	}
  
  /*!
   * Calculates statistics and tests whether it passes required conditions in settings
   */
  
  AlignmentPair::AlignmentPair(bam_alignment& _a1, bam_alignment& _a2, const Settings &settings)
  : a1(_a1), a2(_a2), hash_coord(0)
    , a1_unique_start(0), a1_unique_end(0), a1_unique_length(0)
    , a2_unique_start(0), a2_unique_end(0), a2_unique_length(0)
    , end_to_end_length(0), union_length(0), intersection_length(0)
    , pass(false)
  {
    calculate_union_and_unique();
    pass = test(settings);
  }
  
  /*!
   * Calculates statistics about how the alignments overlap
   */
  
  void AlignmentPair::calculate_union_and_unique()
  {
    bool verbose = false;
    
    uint32_t a1_start, a1_end;
    a1.query_stranded_bounds_1(a1_start, a1_end);
    
    uint32_t a2_start, a2_end;
    a2.query_stranded_bounds_1(a2_start, a2_end);
    
    // move adjust so order is a1-a2 on read sequence
    if (a1_start > a2_start)
    {
      swap(a1, a2);
      swap(a1_start, a2_start);
      swap(a1_end, a2_end);
    }

    // this is the coord of the first base in the read
    // (since we required this to match)
    uint32_t a1_reference_start, a1_reference_end;
    a1.reference_stranded_bounds_1(a1_reference_start, a1_reference_end);
    hash_coord = a1_reference_start;
    ASSERT(a1_start == 1, "Read does not match from beginning.");
    
		int32_t intersection_start = a2_start;
		int32_t intersection_end = a1_end;
		this->intersection_length = intersection_end - intersection_start + 1;
    
    int32_t union_start = a1_start;
		int32_t union_end = a2_end;
    this->end_to_end_length = union_end - union_start + 1;
    // Note: last term subtracts missing bases when the two alignments don't overlap in the middle
		this->union_length = end_to_end_length- max(0, -intersection_length);
    
    a1_unique_start = a1_start;
    a1_unique_end = min(a2_start - 1, a1_end);
    this->a1_unique_length = max(a1_unique_end - a1_unique_start + 1, 0);
    
    a2_unique_start = max(a1_end + 1, a2_start);
    a2_unique_end = a2_end;
    this->a2_unique_length = max(a2_unique_end - a2_unique_start + 1, 0);
    
    if (verbose)
    {
      cout << " Read: " << a1.read_name() << endl;
      cout << "=== Match1: " << a1_start << "-" << a1_end << "   Match2: " << a2_start << "-" << a2_end << endl;
			cout << "    Union: " << this->union_length << "   Intersection: " << this->intersection_length << endl;
      
      cout << "    Unique length 1: " << this->a1_unique_length << " Unique length 2:"<< this->a2_unique_length << endl;
    }
  }
  
  /*!
   * Calculates whether the read pair passes required conditions in settings
   */
  
  bool AlignmentPair::test(const Settings& settings)
  {
    int32_t intersection_length_negative = -min(0, intersection_length);
    int32_t intersection_length_positive = max(0, intersection_length);
    
    int32_t  scaled_maximum_junction_sequence_insertion_overlap_length_fraction
      = static_cast<int32_t>(floor(static_cast<double>(a1.read_length()) * settings.maximum_junction_sequence_insertion_overlap_length_fraction));
    
		//// Require negative overlap (inserted unique sequence length) to be less than some value
		if (intersection_length_negative > scaled_maximum_junction_sequence_insertion_overlap_length_fraction)
			return false;
    
		if (settings.maximum_junction_sequence_insertion_length &&
        (intersection_length_negative > static_cast<int32_t>(settings.maximum_junction_sequence_insertion_length)))
			return false;
    
    //// Require positive overlap (shared by both ends) to be less than some value
    if (intersection_length_positive > scaled_maximum_junction_sequence_insertion_overlap_length_fraction)
			return false;
    
    if (settings.maximum_junction_sequence_overlap_length && 
        (intersection_length_positive > static_cast<int32_t>(settings.maximum_junction_sequence_overlap_length)))
			return false;
    
		//// Require both ends to extend a certain minimum length outside of the overlap

    // This can be an absolute number or a fraction of the read length
    int32_t scaled_required_both_unique_length_per_side 
      = static_cast<int32_t>(ceil(static_cast<double>(a1.read_length()) * settings.required_both_unique_length_per_side_fraction));

		if (a1_unique_length < static_cast<int32_t>(settings.required_both_unique_length_per_side))
			return false;
        
    if (a1_unique_length < scaled_required_both_unique_length_per_side)
      return false;
    
		if (a2_unique_length < static_cast<int32_t>(settings.required_both_unique_length_per_side))
			return false;
    
    if (a2_unique_length < scaled_required_both_unique_length_per_side)
      return false;

		//// Require one end to extend a higher minimum length outside of the overlap
		if ((a1_unique_length <  static_cast<int32_t>(settings.required_one_unique_length_per_side))
        && (a2_unique_length <  static_cast<int32_t>(settings.required_one_unique_length_per_side)))
			return false;
    
		//// Test all of the normal criteria for counting a match to the reference
		if (end_to_end_length < static_cast<int32_t>(settings.require_match_length))
			return false;
    
    if (end_to_end_length < settings.require_match_fraction * static_cast<double>(a1.read_length()) )
      return false;

    return true;
  }

} // namespace breseq
