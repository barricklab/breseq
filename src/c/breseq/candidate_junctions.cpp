
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
  
  class mismatch_map_class : public map<bam_alignment*,double> {
  public:
    inline bool operator() (alignment_list::iterator a1, alignment_list::iterator a2) { return (*this)[a1->get()] < (*this)[a2->get()]; } 
  } mismatch_map;
  
  
  bool sort_by_mismatches (counted_ptr<bam_alignment>& a1, counted_ptr<bam_alignment>& a2) { return mismatch_map[a1.get()] < mismatch_map[a2.get()]; }  
  
  
  /*! Filter a list of alignments to only those that are eligible for mapping
   */
  
 uint32_t eligible_read_alignments(const Settings& settings, const cReferenceSequences& ref_seq_info, alignment_list& alignments)
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
    
    // @JEB: We would ideally sort by alignment score here
    // the number of mismatches is our current proxy for this.
    
    //@JEB This method of sorting may be slower than alternatives
    //     Ideally, the scores should be hashes and only references should be sorted.
    mismatch_map.clear();
    for (alignment_list::iterator it = alignments.begin(); it != alignments.end(); it++)
    {  
      bam_alignment* ap = it->get(); // we are saving the pointer value as the map key
      uint32_t i = alignment_mismatches(*ap, ref_seq_info);
      mismatch_map[ap] = static_cast<double>(i);
    }
    
    alignments.sort(sort_by_mismatches);
    
    if (verbose)
    {
      for (alignment_list::iterator it = alignments.begin(); it != alignments.end(); it++)
      {
        bam_alignment& a = *(it->get());
        cerr << a.query_start_1() << "-" << a.query_end_1() << " ";
        cerr << a.query_match_length()-mismatch_map[it->get()] << "\n";
      }
    }
    
    
    // how many reads share the best score?
    uint32_t last_best(0);
    uint32_t best_score = static_cast<uint32_t>(mismatch_map[(alignments.front().get())]);
    
    // no scores meet minimum
    for (alignment_list::iterator it = alignments.begin()++; it != alignments.end(); it++)
    {
      if (mismatch_map[it->get()] != best_score) break;
      last_best++;
    }
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
    
    // Note that the score we return is higher for better matches so we negative this value...
    return alignments.front()->read_length() - best_score;
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
    
    if (settings.require_complete_match)
    {
      if (!a.beginning_to_end_match())
        return false; 
    }
    if (settings.max_read_mismatches >= 0)
    {
      int32_t mismatches = alignment_mismatches(a, ref_seq_info);
      if (mismatches > settings.max_read_mismatches)
        return false; 
    }
    
    return true;
  }

  
  
  /*! Preprocess read alignments for predicting junctions
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
    
		// includes best matches as they are
		string preprocess_junction_best_sam_file_name = settings.preprocess_junction_best_sam_file_name;
    string reference_fasta_file_name = settings.reference_fasta_file_name;
    
    tam_file BSAM(preprocess_junction_best_sam_file_name, reference_fasta_file_name, ios_base::out);
    
		for (uint32_t index = 0; index < settings.read_files.size(); index++)
		{
			cReadFile read_file = settings.read_files[index];
			cerr << "  READ FILE::" << read_file.m_base_name << endl;
      
			string reference_sam_file_name = Settings::file_name(settings.reference_sam_file_name, "#", read_file.m_base_name);
      tam_file tam(reference_sam_file_name, reference_fasta_file_name, ios_base::in);
      
			// includes all matches, and splits long indels
      string preprocess_junction_split_sam_file_name = Settings::file_name(settings.preprocess_junction_split_sam_file_name, "#", read_file.m_base_name);
      tam_file PSAM(preprocess_junction_split_sam_file_name, reference_fasta_file_name, ios_base::out);
      
			alignment_list alignments;
			uint32_t i = 0;
			while (tam.read_alignments(alignments, false))
			{
				if (++i % 10000 == 0)
					cerr << "    ALIGNED READ:" << i << endl;
        
				// for testing...
				if (settings.candidate_junction_read_limit != 0 && i > settings.candidate_junction_read_limit) break;
        
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
    
    (void)settings; //TODO? unused?
    (void)summary; //TODO: track statistics?
    
    assert(min_indel_split_len >= 0);
    
    alignment_list untouched_alignments;
    uint32_t alignments_written = 0;      
    
    for(alignment_list::const_iterator it = alignments.begin(); it != alignments.end(); it++) 
    {
      uint32_t* cigar_list = (*it)->cigar_array();
      bool do_split = false;
      for(uint32_t i=0; i<(*it)->cigar_array_length(); i++)
      {
        uint32_t op = cigar_list[i] & BAM_CIGAR_MASK;
        uint32_t len = cigar_list[i] >> BAM_CIGAR_SHIFT;
        if (((op == BAM_CINS) || (op == BAM_CDEL)) && (len >= static_cast<uint32_t>(min_indel_split_len)))
        {
          do_split = true;
        }
      }
      
      if (do_split)
      {
        PSAM.write_split_alignment(min_indel_split_len, **it);
        alignments_written += 2;
      }
      else
      {
        untouched_alignments.push_back(*it);
      }
      
    }

    //## @JEB POSSIBLE OPTIMIZATION
    //## We could actually be a little smarter than this and
    //## Not write when we know that this matches so much of the
    //## Read that it cannot be used in a pair to create a candidate junction
    //##
    //## Use $self->{required_both_unique_length_per_side} to rule them out.
    
    //## Don't write in possible junction file when it covers the entire read!
    for(alignment_list::iterator it = untouched_alignments.begin(); it != untouched_alignments.end();) 
    {
      if ((*it)->beginning_to_end_match())
        it = untouched_alignments.erase(it);
      else
        it++;
    }
        
    // Don't write if there is only one alignment to be written,
    // it takes at least two to make a candidate junction.
    if (alignments_written + untouched_alignments.size() > 1)
    {
      PSAM.write_alignments(0, untouched_alignments);
    }
  }
    
	/*! Predicts candidate junctions
	 */
	void CandidateJunctions::identify_candidate_junctions(const Settings& settings, Summary& summary, const cReferenceSequences& ref_seq_info)
	{
		(void)summary; // TODO: save statistics
    bool verbose = true;
    
		// hash by junction sequence concatenated with name of counts
		map<string, map<string, CandidateJunction>, CandidateJunction::Sorter> candidate_junctions;
    
		// summary data for this step
		Summary::CandidateJunctionSummaryData hcs;
    
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
    
		//
		// Calculate pos_hash score for each candidate junction now that we have the complete list of read match support
		//
		for (map<string, map<string, CandidateJunction>, CandidateJunction::Sorter>::iterator outer_it = candidate_junctions.begin(); outer_it != candidate_junctions.end(); outer_it++)
		{
      cout << "Sequence: " << outer_it->first << endl;

			for (map<string, CandidateJunction>::iterator it = (*outer_it).second.begin(); it != (*outer_it).second.end(); it++)
			{
				string junction_id = (*it).first;
				CandidateJunction& cj = (*it).second;
				cj.pos_hash_score = cj.read_begin_hash.size();
        
        // prints out
				if (verbose)
				{
					cout << ">>>" << junction_id << endl;
          cout << "  Pos Hash Score: " << cj.pos_hash_score << endl;
          // this prints out the entire read_begin_hash
          /*
          for (map<uint32_t,uint32_t>::iterator rhbit=cj.read_begin_hash.begin(); rhbit != cj.read_begin_hash.end(); rhbit++)
          {
            cout << "   " << rhbit->first << " " << rhbit->second << endl;
          }
           */
				}
			}
		}
    
		//
		// Combine hash into a list, retaining only one item (the best pos_hash score) for each unique sequence 
    // (and also its reverse complement)
		//
    
		map<int32_t, int32_t> observed_pos_hash_score_distribution;
    
		vector<CombinedCandidateJunction> combined_candidate_junctions;
    
		map<string, int32_t> handled_seq;
		map<string, CombinedCandidateJunction> ids_to_print;
    
		// Not sorted like in perl script to get reproducible ordering, because map is pre-sorted by CandidateJunction::Sorter
		for (map<string, map<string, CandidateJunction>, CandidateJunction::Sorter>::iterator outer_it = candidate_junctions.begin(); outer_it != candidate_junctions.end(); outer_it++)
		{
			string junction_seq = (*outer_it).first;
			if (verbose) cout << "Handling " << junction_seq << endl;
      
			// We may have already done the reverse complement
			if (handled_seq.count(junction_seq) > 0) continue;
      
			// holds all junctions with the same seq
			vector<CombinedCandidateJunction> combined_candidate_junction_list;
      
			string rc_junction_seq = reverse_complement(junction_seq);
      
			for (map<string, CandidateJunction>::iterator it = (*outer_it).second.begin(); it != (*outer_it).second.end(); it++)
			{
				string junction_id = (*it).first;
				// add redundancy to the junction_id
				CandidateJunction cj = (*it).second;
        
				CombinedCandidateJunction ccj = {
					junction_id,			// id
          "",
					cj.pos_hash_score,		// pos_hash_score
					junction_seq,			// seq
					rc_junction_seq			// rc_seq
				};
				combined_candidate_junction_list.push_back(ccj);
			}
			handled_seq[junction_seq]++;
      
			// add the reverse complement
			if (candidate_junctions.count(rc_junction_seq) > 0)
			{
				for (map<string, CandidateJunction>::iterator it = candidate_junctions[rc_junction_seq].begin(); it != candidate_junctions[rc_junction_seq].end(); it++)
				{
					string junction_id = (*it).first;
					// add redundancy to the junction_id (reversed)
					CandidateJunction cj = (*it).second;
          
					CombinedCandidateJunction ccj = {
						junction_id,			// id
            "",
						cj.pos_hash_score,		// pos_hash_score
						rc_junction_seq,		// seq
						junction_seq			// rc_seq
					};
					combined_candidate_junction_list.push_back(ccj);
				}
				handled_seq[rc_junction_seq]++;
			}
      
			sort(combined_candidate_junction_list.begin(), combined_candidate_junction_list.end(), CombinedCandidateJunction::sort_by_score_unique_coord);
			CombinedCandidateJunction& best_candidate_junction = combined_candidate_junction_list[0];
      
			// Save the score in the distribution
			add_score_to_distribution(observed_pos_hash_score_distribution, best_candidate_junction.pos_hash_score);
      
			// Check minimum requirements
			if (best_candidate_junction.pos_hash_score < settings.minimum_candidate_junction_pos_hash_score)
				continue;
      
      // ASSIGN REDUNDANCY
      // =================
      // Assign redundancy to all sides of the BEST junction that are not in EVERY candidate junction with the same sequence
      for (vector<CombinedCandidateJunction>::iterator it=combined_candidate_junction_list.begin()++; it != combined_candidate_junction_list.end(); it++)
      {
        CombinedCandidateJunction& test_candidate_junction(*it);
        
        for (uint32_t best_side_index=0; best_side_index<=1; best_side_index++)
        {
          bool found = false;
          for (uint32_t test_side_index=0; test_side_index<=1; test_side_index++)
          {
            if (best_candidate_junction.junction_info.sides[best_side_index] == test_candidate_junction.junction_info.sides[test_side_index])
              found=true;
          }
          
          // didn't find this side -- mark as redundant
          if (!found)
          {
            best_candidate_junction.junction_info.sides[best_side_index].redundant = true;
            if (verbose) cout << "Marking side " << best_side_index << " as redundant." << endl;
          }
        }
      }
      
      // ONLY NOW can we create the correct junction id key
      best_candidate_junction.junction_key = best_candidate_junction.junction_info.junction_key();
      
			// Make sure it isn't a duplicate junction id -- this should NEVER happen and causes downstream problem.
			// <--- Begin sanity check
			if (ids_to_print.count(best_candidate_junction.junction_key))
			{
        CombinedCandidateJunction& ccj = ids_to_print[best_candidate_junction.junction_key];
        
				cout << "Attempt to create junction candidate with duplicate id: " << best_candidate_junction.junction_key << endl;
        
				cout << "==Existing junction==" << endl;      
        cout << "  id: " << ccj.junction_key << endl;
        cout << "  pos_hash_score: " << ccj.pos_hash_score << endl;
        cout << "  pos_hash_score: " << ccj.pos_hash_score << endl;
        cout << "  seq: " << ccj.seq << endl;
        cout << "  rc_seq: " << ccj.rc_seq << endl;   
        
				cout << "==New junction==" << endl;
        cout << "  id: " << best_candidate_junction.junction_key << endl;
        cout << "  pos_hash_score: " << best_candidate_junction.pos_hash_score << endl;
        cout << "  pos_hash_score: " << best_candidate_junction.pos_hash_score << endl;
        cout << "  seq: " << best_candidate_junction.seq << endl;
        cout << "  rc_seq: " << best_candidate_junction.rc_seq << endl;  
        
				assert (best_candidate_junction.seq == ids_to_print[best_candidate_junction.junction_key].seq);
				exit(-1);
			}
			ids_to_print[best_candidate_junction.junction_key] = best_candidate_junction;
			// <--- End sanity check
      
			combined_candidate_junctions.push_back(best_candidate_junction);
		}
    
		sort(combined_candidate_junctions.begin(), combined_candidate_junctions.end(), CombinedCandidateJunction::sort_by_scores_and_seq_length);
    
    if (verbose)
    {
      for (vector<CombinedCandidateJunction>::iterator it=combined_candidate_junctions.begin(); it < combined_candidate_junctions.end(); it++ )
      {
        cout << "ID: " << it->junction_info.junction_key() << endl;
        cout << "  pos_hash_score: " << it->pos_hash_score << endl;
        cout << "  seq: " << it->seq << endl;
        cout << "  rc_seq: " << it->rc_seq << endl;
      }
    }
    
		///
		// Limit the number of candidate junctions that we print by:
		//   (1) A maximum number of candidate junctions
		//   (2) A maximum length of the sequences in candidate junctions
		///
    
		cerr << "  Taking top candidate junctions..." << endl;
    
		// adding up the lengths might be too time-consuming to be worth it...
		int32_t total_cumulative_cj_length = 0;
		int32_t total_candidate_junction_number = combined_candidate_junctions.size();
		for (uint32_t j = 0; j < combined_candidate_junctions.size(); j++)
			total_cumulative_cj_length += combined_candidate_junctions[j].seq.size();
    
		uint32_t cumulative_cj_length = 0;
		int32_t lowest_accepted_pos_hash_score = 0;
    
		// Right now we limit the candidate junctions to have a length no longer than the reference sequence.
		uint32_t cj_length_limit = static_cast<uint32_t>(summary.sequence_conversion.total_reference_sequence_length * settings.maximum_candidate_junction_length_factor);
		uint32_t maximum_candidate_junctions = settings.maximum_candidate_junctions;
		uint32_t minimum_candidate_junctions = settings.minimum_candidate_junctions;
    
		fprintf(stderr, "  Minimum number to keep: %7d \n", minimum_candidate_junctions);
		fprintf(stderr, "  Maximum number to keep: %7d \n", maximum_candidate_junctions);
		fprintf(stderr, "  Maximum length to keep: %7d bases\n", cj_length_limit);
    
		cerr << "    Initial: Number = " << total_candidate_junction_number << ", Cumulative Length = " << total_cumulative_cj_length << " bases" << endl;
    
		if ((settings.maximum_candidate_junctions > 0) && (combined_candidate_junctions.size() > 0))
		{
			vector<CombinedCandidateJunction> remaining_ids;
			vector<CombinedCandidateJunction> list_in_waiting;
			int32_t add_cj_length = 0;
			int32_t num_duplicates = 0;
      
			i = 0;
			uint32_t current_pos_hash_score = combined_candidate_junctions[i].pos_hash_score;
      
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
        
				current_pos_hash_score = combined_candidate_junctions[i].pos_hash_score;
				while (
               (i < combined_candidate_junctions.size())
               && (combined_candidate_junctions[i].pos_hash_score == current_pos_hash_score)
               )
				{
					CombinedCandidateJunction c = combined_candidate_junctions[i];
					list_in_waiting.push_back(c);
					add_cj_length += c.seq.size();
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
    
		sort(combined_candidate_junctions.begin(), combined_candidate_junctions.end(), CombinedCandidateJunction::sort_by_ref_seq_coord);
    
    cFastaFile out(settings.candidate_junction_fasta_file_name, ios_base::out);
    
		for (uint32_t j = 0; j < combined_candidate_junctions.size(); j++)
		{
			CombinedCandidateJunction junction = combined_candidate_junctions[j];
			cFastaSequence seq = { junction.junction_info.junction_key(), "", junction.seq };
			out.write_sequence(seq);
		}
		out.close();
    
		summary.candidate_junction = hcs;
	}
  

	bool CandidateJunctions::alignment_pair_to_candidate_junction(
                                                            const Settings& settings, 
                                                            Summary& summary, 
                                                            const cReferenceSequences& ref_seq_info, 
                                                            AlignmentPair& ap,
                                                            SingleCandidateJunction& junction
                                                            )
	{    
		bool verbose = false;

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
		// REL606__1__1__REL606__4629812__0__0
		// means the junction sequence is 36-1 + 4629812-4629777 from the reference sequence
		//
		// On the LEFT side:  0 means this is highest coord of alignment, junction seq begins at lower coord
		//                    1 means this is lowest coord of alignment, junction seq begins at higher coord
		// On the RIGHT side: 0 means this is highest coord of alignment, junction seq continues to lower coord
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
		if (overlap > 0)
		{
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
					cout << "Rejecting alignment because there is not not enough overlap" << endl;
				int32_t dont_care;
        
        AlignmentPair ap(q1, q2, settings);
				if (!ap.pass) return false;
			}

			//re-calculate the overlap
			overlap = -1 * (q2_start - q1_end - 1);
			if (verbose)
				cout << "=== overlap corrected for mismatches " << overlap << endl;
		}
    
    
    // create hash coords AFTER overlap adjustment
		int32_t hash_coord_1 = (hash_strand_1) ? r1_start : r1_end;
		int32_t hash_coord_2 = (hash_strand_2) ? r2_start : r2_end;

    // if there are multiple ways the two sides could have been aligned...shift them
    // over so as much is included in the lower reference coordinate side as possible
      
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

		// record what parts of the reference sequence were actually matched on each side
		// this is to determine whether that side was redundant or unique in the reference

		int32_t ref_seq_matched_length_1 = r1_end - r1_start + 1;
		int32_t ref_seq_matched_length_2 = r2_end - r2_start + 1;

		// create the sequence of the candidate junction
		string junction_seq_string = "";

		// first end
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
		// and NOT present in the reference genome
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
		else // ($m_2->{hash_strand} * $m_2->{read_side} == -1)
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

		// matched reference sequence, EXCLUDING OVERLAP because it is removed by offsets above.
		string ref_seq_matched_1 = ref_seq_1.substr(r1_start - 1, r1_end - r1_start + 1);
		string ref_seq_matched_2 = ref_seq_2.substr(r2_start - 1, r2_end - r2_start + 1);

		// want to be sure that lowest ref coord is always first for consistency
		if ( hash_seq_id_1.compare(hash_seq_id_2) > 0 || (hash_seq_id_1.compare(hash_seq_id_2) == 0 && hash_coord_2 < hash_coord_1) )
		{
			swap(hash_coord_1, hash_coord_2);
			swap(hash_strand_1, hash_strand_2);
			swap(hash_seq_id_1, hash_seq_id_2);
			swap(flanking_left, flanking_right);
			swap(ref_seq_matched_1, ref_seq_matched_2);

			junction_seq_string = reverse_complement(junction_seq_string);
			unique_read_seq_string = reverse_complement(unique_read_seq_string);
		}

		JunctionInfo::Side side_1(
			hash_seq_id_1,	
			hash_coord_1,	
			hash_strand_1
		);
    
    JunctionInfo::Side side_2(
			hash_seq_id_2,	
			hash_coord_2,	
			hash_strand_2
		);
    
		JunctionInfo new_junction(
			side_1, side_2,
			overlap, 				
			unique_read_seq_string,
			flanking_left,
			flanking_right
		);
    
    
		string join1[] = { hash_seq_id_1, to_string(hash_coord_1), to_string(hash_strand_1) };
		string junction_coord_1 = join(join1, "::");
		string join2[] = { hash_seq_id_2, to_string(hash_coord_2), to_string(hash_strand_2) };
		string junction_coord_2 = join(join2, "::");

		if (verbose)
		{
      string junction_id = new_junction.junction_key();
			cout << "READ ID: " << a1.read_name() << endl;
			cout << "JUNCTION ID: " << junction_id << endl;
		}

		assert(junction_seq_string.size() > 0); // die "Junction sequence not found: $junction_id " . $q1->qname . " " . $a2->qname  if (!$junction_seq_string);
		assert(junction_seq_string.size() == flanking_left + flanking_right + static_cast<uint32_t>(abs(overlap))); // die "Incorrect length for $junction_seq_string: $junction_id " . $q1->qname . " " . $a2->qname if (length $junction_seq_string != $flanking_left + $flanking_right + abs($overlap));

    // Set return values
    junction.junction_info = new_junction;
    junction.sequence = junction_seq_string;
    junction.read_begin_coord = read_begin_coord;
    junction.side_1_ref_seq = ref_seq_matched_1;
    junction.side_2_ref_seq = ref_seq_matched_2;
    junction.junction_coord_1 = junction_coord_1;
    junction.junction_coord_2 = junction_coord_2;
    
		return true;
	}

	void CandidateJunctions::alignments_to_candidate_junctions(const Settings& settings, Summary& summary, const cReferenceSequences& ref_seq_info, map<string, map<string, CandidateJunction>, CandidateJunction::Sorter>& candidate_junctions, alignment_list& alignments)
	{
		bool verbose = false;

		if (verbose)
		{
			cout << endl << "###########################" << endl;
			cout << alignments.front()->read_name();
			cout << endl << "###########################" << endl;
		}

		// Must still have multiple matches to support a new junction.
		if (alignments.size() <= 1)
			return;

		// Now is our chance to decide which groups of matches are compatible,
		// to change their boundaries and to come up with a final list.

		if (verbose)
		{
			cout << alignments.front()->read_name() << endl;
			cout << "Total matches: " << alignments.size() << endl;
		}

    // list1 contains all matches starting at the beginning of the read
    // list2 contains all other matches
    // you can never have a pair that doesn't have a member from each list
		alignment_list list1, list2;

    // Try only pairs where one match starts at the beginning of the read >>>
		// This saves a number of comparisons and gets rid of a lot of bad matches.
		int32_t max_alignment_length = 0;

    for (alignment_list::iterator it=alignments.begin(); it != alignments.end(); it++)
		{
			bam_alignment_ptr a = *it;
      
			uint32_t a_start, a_end;
      a->query_stranded_bounds_1(a_start, a_end);

			if (verbose) cout << "(" << a_start << ", " << a_end << ")" << endl;

			int32_t length = a_end - a_start + 1;
			if (length > max_alignment_length)
				max_alignment_length = length;

			if (a_start == 1)
				list1.push_back(a);
			else
				list2.push_back(a);
		}

		// Pairs must CLEAR the maximum length of any one alignment by a certain amount
		int32_t required_union_length = max_alignment_length + settings.required_extra_pair_total_length;
    int32_t max_union_length = 0;
    
		// The first match in this category is the longest
		if (verbose)
		{
			cout << "  List1: " << list1.size() << endl;
			cout << "  List2: " << list2.size() << endl;
		}

		vector<AlignmentPair> passed_pair_list;

		// Try adding together each pair of matches to make a junction, by looking at read coordinates
    for (alignment_list::iterator it1 = list1.begin(); it1 != list1.end(); it1++)
		{
			bam_alignment& a1 = *(it1->get());

      for (alignment_list::iterator it2 = list2.begin(); it2 != list2.end(); it2++)
			{
				bam_alignment& a2 = *(it2->get());

				// constructing the alignment pair calculates statistics about their overlap
        // and tests the guards that are in Settings.
        
        AlignmentPair ap(a1, a2, settings);

				if (ap.pass && (ap.union_length >= required_union_length))
				{
					// right now we use the union length as a score
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
    
    // create a list of all the candidate junctions
    vector<SingleCandidateJunction> junctions;
		for (uint32_t i = 0; i < passed_pair_list.size(); i++)
		{
			AlignmentPair& ap = passed_pair_list[i];
      
      SingleCandidateJunction cj;
			bool passed = alignment_pair_to_candidate_junction(settings, summary, ref_seq_info, ap, cj);
			if (!passed) continue;
      
      junctions.push_back(cj);
    }

		if (verbose) cout << "  Junctions: " << junctions.size() << endl;

		// Done if everything already ruled out...
		if (junctions.size() == 0) return;

		if (verbose) cout << alignments.front()->read_name() << endl;

    // Add these to the main hash (by sequence)
		for (uint32_t i = 0; i < junctions.size(); i++)
		{
			SingleCandidateJunction& jct = junctions[i];
			JunctionInfo& junction_info = jct.junction_info;
			string junction_seq_string = jct.sequence;
			int32_t read_begin_coord = jct.read_begin_coord;

			string side_1_ref_seq = jct.side_1_ref_seq;
			string side_2_ref_seq = jct.side_2_ref_seq;

			string junction_id = junction_info.junction_key();
			if (verbose) cout << junction_id << endl;

			// initialize candidate junction if it didn't exist
			// they are redundant by default, until proven otherwise

      if ((candidate_junctions.count(junction_seq_string) == 0) || (candidate_junctions[junction_seq_string].count(junction_id) == 0)) {
        
        // these values get written over immediately
        CandidateJunction new_cj;
        candidate_junctions[junction_seq_string][junction_id] = new_cj;
      }
      
      CandidateJunction& cj(candidate_junctions[junction_seq_string][junction_id]);

			// Update score of junction and the redundancy of each side
      cj.read_begin_hash[read_begin_coord]++;
			candidate_junctions[junction_seq_string][junction_id] = cj;
		}
	}
  
  /*!
   * Calculates statistics and tests whether it passes required conditions in settings
   */
  
  AlignmentPair::AlignmentPair(bam_alignment& _a1, bam_alignment& _a2, const Settings &settings)
    : a1(_a1), a2(_a2), union_length(0), a1_unique_length(0), a2_unique_length(0), pass(false)
    , intersection_length(0), a1_length(0), a2_length(0)
  {
    calculate_union_and_unique();
    pass = test(settings);
  }
  
  /*!
   * Calculates the union and unique lengths of the read pair
   */
  
  void AlignmentPair::calculate_union_and_unique()
  {
    bool verbose = false;
   
    uint32_t a1_start, a1_end;
    a1.query_stranded_bounds_1(a1_start, a1_end);
    
    uint32_t a2_start, a2_end;
    a2.query_stranded_bounds_1(a2_start, a2_end);
    
		if (verbose)
			cout << "=== Match1: " << a1_start << "-" << a1_end << "   Match2: " << a2_start << "-" << a2_end << endl;
    
    //@JEB TODO Is it ok to remove the following section?
    
		// 0. Require one match to start at the beginning of the read
		//if (a1_start != 1 && a2_start != 1)
		//	return false;
		// Already checked when two lists were constructed: TEST AND REMOVE

		a1_length = a1_end - a1_start + 1;
		a2_length = a2_end - a2_start + 1;
    
		int32_t union_start = (a1_start < a2_start) ? a1_start : a2_start;
		int32_t union_end = (a1_end > a2_end) ? a1_end : a2_end;
		union_length = union_end - union_start + 1;
    
		int32_t intersection_start = (a1_start > a2_start) ? a1_start : a2_start;
		int32_t intersection_end = (a1_end < a2_end) ? a1_end : a2_end;
		this->intersection_length = intersection_end - intersection_start + 1;
    
		if (intersection_length < 0)
			this->union_length += this->intersection_length;
    
    if (verbose)
			cout << "    Union: " << this->union_length << "   Intersection: " << this->intersection_length << endl;
  }
  
  /*!
   * Calculates whether the read pair passes required conditions in settings
   */
  
  bool AlignmentPair::test(const Settings& settings)
  {
		int32_t intersection_length_positive = (intersection_length > 0) ? intersection_length : 0;
		int32_t intersection_length_negative = (intersection_length < 0) ? -1 * intersection_length : 0;
    
		//// 1. Require maximum negative overlap (inserted unique sequence length) to be less than some value
		if (intersection_length_negative > static_cast<int32_t>(settings.maximum_inserted_junction_sequence_length))
			return false;
    
		//// 2. Require both ends to extend a certain minimum length outside of the overlap
		if (a1_length < intersection_length_positive + static_cast<int32_t>(settings.required_both_unique_length_per_side))
			return false;
    
		if (a2_length < intersection_length_positive + static_cast<int32_t>(settings.required_both_unique_length_per_side))
			return false;
    
		//// 3. Require one end to extend a higher minimum length outside of the overlap
		if ((a1_length < intersection_length_positive + static_cast<int32_t>(settings.required_one_unique_length_per_side))
        && (a2_length < intersection_length_positive + static_cast<int32_t>(settings.required_one_unique_length_per_side)))
			return false;
    
		//// 4. Require both matches together to cover a minimum part of the read.
		if (union_length < static_cast<int32_t>(settings.require_match_length))
			return false;

    return true;
  }

  
  



} // namespace breseq
