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

#include "libbreseq/resolve_alignments.h"

#include "libbreseq/genome_diff.h"
#include "libbreseq/fastq.h"
#include "libbreseq/fasta.h"
#include "libbreseq/alignment.h"
#include "libbreseq/annotated_sequence.h"

using namespace std;

namespace breseq {
    
// Compares matches to candidate junctions with matches to original genome
void resolve_alignments(
                        Settings& settings,
                        Summary& summary,
                        cReferenceSequences& ref_seq_info,
                        bool junction_prediction,
                        cReadFiles &read_files,
                        const uint32_t max_read_length,
                        const uint32_t alignment_read_limit
                        ) 
{
  (void)max_read_length; // @TODO: needed?
    
	bool verbose = false;
    
  // Load the trims @JEB could speed up by only loading this once in main
  SequenceTrimsList trims_list;
  read_trims(trims_list, ref_seq_info, settings.reference_trim_file_name);
  
  // create junction trims
  // directly from FASTA
  SequenceTrimsList junction_trims_list;
  
  if (junction_prediction)
  {  
    cFastaFile ff(settings.candidate_junction_fasta_file_name, ios::in);
    cFastaSequence fs;
    while (ff.read_sequence(fs))
    {
      SequenceTrims st(fs.m_sequence);
      junction_trims_list.push_back(st);
    }
  }
    
	// ####
	// ##	Junction sequences
	// ####

	//## if there were no candidate junctions (file is empty) then we seg fault if we try to use samtools on it...
  cReferenceSequences junction_ref_seq_info;
	if (junction_prediction
			&& !file_exists(settings.candidate_junction_fasta_file_name.c_str())
			&& !file_empty(settings.candidate_junction_fasta_file_name.c_str())
		)
  {
		junction_prediction = 0;
  }
    
	vector<JunctionInfo> junction_info_list;

	tam_file* junction_tam = NULL;
	tam_file* reference_tam = NULL;
    

  if (junction_prediction)
	{
    LoadFeatureIndexedFastaFile(junction_ref_seq_info, "", settings.candidate_junction_fasta_file_name);
		string junction_sam_file_name = settings.file_name(settings.candidate_junction_sam_file_name, "#", read_files[0].m_base_name);
		tam_file junction_tam(junction_sam_file_name, settings.candidate_junction_fasta_file_name, ios::in);

		//## Preload all of the information about junctions
		//## so that we only have to split the names once
		for (int i = 0; i < junction_tam.bam_header->n_targets; i++)
		{
			JunctionInfo ji = junction_name_split(junction_tam.bam_header->target_name[i]);
			junction_info_list.push_back(ji);
		}

	}

	//####
	//##	Output files
	//####

	genome_diff gd;
    
  tam_file resolved_reference_tam(settings.resolved_reference_sam_file_name, settings.reference_fasta_file_name, ios::out);
  tam_file resolved_junction_tam(settings.resolved_junction_sam_file_name, settings.candidate_junction_fasta_file_name, ios::out);
    
  map<string, vector<MatchedJunction> > matched_junction;
	map<string, map<string, MatchedJunction> > degenerate_matches;
	uint32_t reads_processed = 0;

  // stores all junction ids that we have encountered
  map<string,uint32_t> all_junction_ids;

	for (uint32_t fastq_file_index = 0; fastq_file_index < read_files.size(); fastq_file_index++)
	{
    const cReadFile& rf = read_files[fastq_file_index];
    string fastq_file_name = read_files.base_name_to_read_file_name(rf.m_base_name);
    
		cerr << "  READ FILE:" << rf.m_base_name << endl;

		map<string,int32_t> summary_info = make_map<string,int32_t>("unmatched_reads", 0);

		// Traverse the original fastq files to keep track of order
		// b/c some matches may exist in only one or the other file
    
		cFastqFile in_fastq(fastq_file_name, ios::in);

    string this_unmatched_file_name = settings.data_path + "/unmatched."
        + rf.m_base_name + ".fastq";
    cFastqFile out_unmatched_fastq(this_unmatched_file_name, ios::out);
    assert(!out_unmatched_fastq.fail());

		string reference_sam_file_name = settings.file_name(settings.reference_sam_file_name, "#", rf.m_base_name);
    string reference_fasta = settings.reference_fasta_file_name;

		reference_tam = new tam_file(reference_sam_file_name, settings.reference_fasta_file_name, ios::in); 

		if (junction_prediction)
		{
			string junction_sam_file_name = settings.file_name(settings.candidate_junction_sam_file_name, "#", rf.m_base_name);
			junction_tam = new tam_file(junction_sam_file_name, settings.candidate_junction_fasta_file_name, ios::in); 
		}

		alignment_list junction_alignments;

		//proceed through all of the alignments
		if (junction_prediction)
			junction_tam->read_alignments(junction_alignments, false);

		alignment_list reference_alignments;
		reference_tam->read_alignments(reference_alignments, false);

		///
		//  Test each read for its matches to the reference and candidate junctions
		///

		cFastqSequence seq;
		while (in_fastq.read_sequence(seq)) // READ
		{
			reads_processed++;

			if ((alignment_read_limit) && (reads_processed > alignment_read_limit))
				break;

			if (reads_processed % 10000 == 0)
				cerr << "    READS:" << reads_processed << endl;

			if (verbose)
				cerr << "===> Read: " << seq.m_name << endl;

			uint32_t best_junction_score = 0;
			uint32_t best_reference_score = 0;

			// Does this read have eligible candidate junction matches?
			alignment_list this_junction_alignments;

			if (verbose)
			{
				cerr << " Before Overlap Reference alignments = " << reference_alignments.size() << endl;
				cerr << " Before Overlap Junction alignments = " << junction_alignments.size() << endl;
			}
      
      if (verbose && (junction_alignments.size() > 0))
      {
        cerr << " Junction SAM read name: " << junction_alignments.front()->read_name() <<endl;
      }

			if ((junction_alignments.size() > 0) && (seq.m_name == junction_alignments.front()->read_name()))
			{

				this_junction_alignments = junction_alignments;
				junction_tam->read_alignments(junction_alignments, false);

				///
				// Matches to candidate junctions MUST overlap the junction.
				//
				// Reduce this list to those that overlap ANY PART of the junction.
				// Alignments that extend only into the overlap region, are only additional
				//  evidence for predicted junctions and NOT support for a new junction on
				// their own. (They will also match the original reference genome equally well).
        // ... but this last point only if overlap >=0 for the junction
				///

				for (alignment_list::iterator it = this_junction_alignments.begin(); it != this_junction_alignments.end(); )
        {
					if (!_alignment_overlaps_junction(junction_info_list, it->get()))
          {
						this_junction_alignments.erase(it++);
          }
          else
          {
            it++; 
          }
        }

				best_junction_score = _eligible_read_alignments(settings, junction_ref_seq_info, this_junction_alignments);

				if (verbose)
					cerr << " Best junction score: " << best_junction_score
							<< endl;
			}

      if (verbose && (reference_alignments.size() > 0))
      {
        cerr << " Reference SAM read name: " << reference_alignments.front()->read_name() <<endl;
      }
      
			// Does this read have eligible reference sequence matches?
			alignment_list this_reference_alignments;
			if ((reference_alignments.size() > 0) && (seq.m_name == reference_alignments.front()->read_name()))
			{

				this_reference_alignments = reference_alignments;
				reference_tam->read_alignments(reference_alignments, false);

				best_reference_score = _eligible_read_alignments(settings, ref_seq_info, this_reference_alignments);

				if (verbose)
					cerr << " Best reference score: " << best_reference_score << endl;
			}

			// Nothing to be done if there were no eligible matches to either
			// Record in the unmatched FASTQ data file
			if ((this_junction_alignments.size() == 0) && (this_reference_alignments.size() == 0))
			{
				summary_info["unmatched_reads"]++;
				out_unmatched_fastq.write_sequence(seq);
			}

			///
			// Determine if the read has a better match to a candidate junction
			// or to the reference sequence.
			///

			/// There are three possible kinds of reads at this point
			//
			// 1: Read has a best match to the reference genome
			// --> Write this match and we are done
			// 2: Read has a best match (or multiple best matches) to junctions
			// --> Keep an item that describes these matches
			// 3: Read has an equivalent match to the reference genome
			//      and goes into the overlap part of a junction condidate
			// --> Keep an item that is not used during scoring
			///

			//#			if (@$this_junction_al)
			//#			{
			//#				# This is the length of the match on the query -- NOT the length of the query
			//#				$best_junction_score = $this_junction_al->[0]->query->length;
			//# THESE SCORES ARE NOT CONSISTENT ACROSS STRANDS DUE TO DIFFERENT INDEL MATCHES
			//#				$best_candidate_junction_score = $ca->aux_get("AS");
			//#			}
			//#
			//#			if (@$this_reference_al)
			//#			{
			//#				# This is the length of the match on the query -- NOT the length of the query
			//#				$best_reference_score = $this_reference_al->[0]->query->length;
			//# THESE SCORES ARE NOT CONSISTENT ACROSS STRANDS DUE TO DIFFERENT INDEL MATCHES
			//#				$best_reference_score = $ra->aux_get("AS");
			//#			}

			// if < 0, then the best match is to the reference
			int32_t mapping_quality_difference = best_junction_score - best_reference_score;

			if (verbose)
			{
				cerr << " Best junction score: " << best_junction_score << endl;
				cerr << " Best reference score: " << best_reference_score << endl;
				cerr << " Mapping quality difference: " << best_reference_score << endl;
				cerr << " Final Reference alignments = " << this_reference_alignments.size() << endl;
				cerr << " Final Candidate junction alignments = " << this_junction_alignments.size() << endl;
			}

			if ((this_junction_alignments.size() == 0) && (this_reference_alignments.size() == 0))
				continue;

			///
			// The best match we found to the reference was no better than the best to the
			// candidate junction. This read potentially supports the candidate junction.
			//
			// ONLY allow EQUAL matches through if they match the overlap only, otherwise
			// you can get predictions of new junctions with all reads supporting them
			// actually mapping perfectly to the reference.
			///

			// best match is to the reference, record in that SAM file.
			if (mapping_quality_difference < 0)
			{
				if (verbose)
					cout << "Best alignment to reference. MQD: " << mapping_quality_difference << endl;

				_write_reference_matches(settings, ref_seq_info, trims_list, this_reference_alignments, resolved_reference_tam, fastq_file_index);
			}
			else
			{
				if (verbose)
					cout << "Best alignment is to candidate junction. MQD: " << mapping_quality_difference << endl;

				MatchedJunction item(
					this_reference_alignments, // reference sequence alignments
					this_junction_alignments, // the BEST candidate junction alignments
					fastq_file_index, // index of the fastq file this read came from
					mapping_quality_difference, 0 // degenerate count
				);

				///
				// Just one best hit to candidate junctions, that is better than every match to the reference
				///
				if ((this_junction_alignments.size() == 1) && (mapping_quality_difference > 0))
				{
					bam_alignment& a = *(this_junction_alignments.front().get());
					string junction_id = junction_tam->bam_header->target_name[a.reference_target_id()];
					matched_junction[junction_id].push_back(item);
          all_junction_ids[junction_id]++;
				}
				///
				// Multiple equivalent matches to junctions and reference, ones with most hits later will win these matches
				// If $mapping_quality_difference > 0, then they will count for scoring
				///
				else
				{
          if (verbose)
            cout << "this_junction_alignments: " << this_junction_alignments.size() << endl;
          
          for(alignment_list::iterator it=this_junction_alignments.begin(); it!=this_junction_alignments.end(); it++)
					{
            bam_alignment& a = *(it->get());
						item.degenerate_count = this_junction_alignments.size(); // mark as degenerate
						string junction_id = junction_tam->bam_header->target_name[a.reference_target_id()];
						degenerate_matches[junction_id][seq.m_name] = item;
            all_junction_ids[junction_id]++;
					}
				}
			} // READ

		} // End loop through every $read_struct

		// save statistics
		summary.alignment_correction.read_file[read_files[0].m_base_name] = summary_info;
    
    // safe only because we know they are always or never used
    if (junction_tam != NULL) delete junction_tam;
		if (reference_tam != NULL) delete reference_tam;

	} // End of Read File loop

  verbose = false;

	///
	// Determine which junctions are real, prefer ones with most matches
	///

	map<int32_t, int32_t> accepted_pos_hash_score_distribution;
	map<int32_t, int32_t> observed_pos_hash_score_distribution;
    
	vector<string> passed_junction_ids;
	vector<string> rejected_junction_ids;
	map<string, CandidateJunction> junction_test_info; // scoring information about junctions

	///
	// Candidate junctions with unique matches
	///

	//sort junction ids based on size of vector

	vector<string> sorted_junction_ids = get_sorted_junction_ids(matched_junction, degenerate_matches, get_keys(all_junction_ids));
    
  if (verbose) cout << "Number of unique matches: " << sorted_junction_ids.size() << endl;


	for (uint32_t i = 0; i < sorted_junction_ids.size(); i++)
	{
		string key = sorted_junction_ids[i];

    if (verbose) 
    {
      cout << "Testing Junction with Unique Matches:" << key << endl;
      cout << "  Number of unique matches:" << matched_junction[key].size() << endl;
      size_t num_degenerate_matches = degenerate_matches.count(key) ? degenerate_matches[key].size() : 0;
      cout << "  Number of degenerate matches:" << num_degenerate_matches << endl;
    }
		bool has_non_overlap_alignment = false;
		bool success = _test_junction(settings, summary, key, matched_junction, degenerate_matches, junction_test_info, ref_seq_info, trims_list, resolved_reference_tam, resolved_junction_tam, has_non_overlap_alignment);

		// save the score in the distribution
		add_score_to_distribution(observed_pos_hash_score_distribution, junction_test_info[key].pos_hash_score);

    if (verbose && !has_non_overlap_alignment) cout << "Does not have nonoverlap alignments" << endl;
    
		// only keep matches that span overlap
		if (has_non_overlap_alignment)
		{
			if (success)
      {
// @JEB we might want to re-sort the list here.
// Strategy would be to use a <list> and shorten it each time
// Then re-sort when we successfully added one that might
// Have eaten up degenerate matches shared with others.
				passed_junction_ids.push_back(key);
        if (verbose) cout << "  PASSED" << endl;
			}
      else
			{
        rejected_junction_ids.push_back(key);
        if (verbose) cout << "  REJECTED" << endl;
      }
    }
  }

  if (verbose)
    cout << "Degenerate matches after handling ones with unique matches: " << degenerate_matches.size() << endl;
    
	// print successful ones out
	if (verbose) cout << "Successful hybrids" << endl;

	//Re-sort
	passed_junction_ids = get_sorted_junction_ids(matched_junction, degenerate_matches, passed_junction_ids);
	rejected_junction_ids = get_sorted_junction_ids(matched_junction, degenerate_matches, rejected_junction_ids);

  if (verbose)
  {
    cout << "passed_junction_ids" << endl;
    for(uint32_t i = 0; i < passed_junction_ids.size(); i++)
    {
      string key = passed_junction_ids[i];
      cout << key << endl;
    }
      
    cout << "matched_junction" << endl;
    for(map<string, vector<MatchedJunction> >::const_iterator it = matched_junction.begin(); it != matched_junction.end(); it++)
    {
      cout << it->first << " " << it->second.size() << endl;
    }
  }
    
	for (uint32_t i = 0; i < passed_junction_ids.size(); i++)
	{
		string key = passed_junction_ids[i];
		if (verbose) cout << key << endl;
		diff_entry item = _junction_to_hybrid_list_item(key, ref_seq_info, matched_junction[key].size(), junction_test_info[key]);
		gd.add(item);

		// save the score in the distribution
		add_score_to_distribution(accepted_pos_hash_score_distribution, junction_test_info[key].pos_hash_score);

		// Create matches from UNIQUE sides of each match to reference genome
		// this fixes, for example appearing to not have any coverage at the origin of a circular DNA fragment
		//  Currently, we do not add coverage to redundantly matched sides because we don't know which copy.
		if (!settings.add_split_junction_sides) continue;

		for (uint32_t j = 0; j < matched_junction[key].size(); j++)
		{
			MatchedJunction& match = matched_junction[key][j];
			bam_alignment& a = *(match.junction_alignments.front().get());
			uint32_t fastq_file_index = match.fastq_file_index;

			if (verbose) {
				cout << ">>>>" << a.read_name() << " #" << match.junction_alignments.size() << endl;
				cout << ">>>>Alignment start-end: " << a.reference_start_1() << "  " << a.reference_end_1() << endl;
			}
      
      // at this point, all degeneracy should have been removed!
      assert(match.junction_alignments.size() == 1);

			for (uint32_t side = 1; side <= 2; side++)
			{
				string side_key = "side_" + to_string(side);

				// Do not count for coverage if it is redundant!!
				if (from_string<int32_t>(item[side_key + "_redundant"])) continue;
        
				// Write out match corresponding to this part to SAM file
				// By trimming in the candidate junctions sequence, rather than on each half,
				// this is done properly.
				Trims trims = _trim_ambiguous_ends(a, junction_trims_list);
        
				resolved_reference_tam.write_moved_alignment(
					a,
          resolved_junction_tam.target_name(a),
					fastq_file_index,
					item[side_key + "_seq_id"],
					from_string<int32_t>(item[side_key + "_position"]),
					from_string<int32_t>(item[side_key + "_strand"]),
					from_string<int32_t>(item[side_key + "_overlap"]),
					side,
					from_string<int32_t>(item["flanking_left"]),
					from_string<int32_t>(item["alignment_overlap"]),
					&trims
				);
			}
		}
	}

	// Save summary statistics
	summary.alignment_correction.new_junctions.observed_pos_hash_score_distribution = observed_pos_hash_score_distribution;
	summary.alignment_correction.new_junctions.accepted_pos_hash_score_distribution = accepted_pos_hash_score_distribution;

	for (uint32_t i = 0; i < rejected_junction_ids.size(); i++)
	{
		string key = rejected_junction_ids[i];
		diff_entry item = _junction_to_hybrid_list_item(key, ref_seq_info, matched_junction[key].size(), junction_test_info[key]);
		add_reject_reason(item, "NJ");
		gd.add(item);
	}

	gd.write(settings.jc_genome_diff_file_name);
}
    
//
//=head2 _eligible_alignments
//
//Title   : _eligible_alignments
//Usage   : _eligible_alignments( );
//Function:
//Returns : Best score
//
//=cut
//
  
class mismatch_map_class : public map<bam_alignment*,double> {
public:
  inline bool operator() (alignment_list::iterator a1, alignment_list::iterator a2) { return (*this)[a1->get()] < (*this)[a2->get()]; } 
} mismatch_map;


bool sort_by_mismatches (counted_ptr<bam_alignment>& a1, counted_ptr<bam_alignment>& a2) { return mismatch_map[a1.get()] < mismatch_map[a2.get()]; }  
  
uint32_t _eligible_read_alignments(const Settings& settings, const cReferenceSequences& ref_seq_info, alignment_list& alignments)
{
	bool verbose = false;

	if (alignments.size() <= 0) return 0;

	// require a minimum length of the read to be mapped
	for (alignment_list::iterator it = alignments.begin(); it != alignments.end(); )
  {
		if ( !_test_read_alignment_requirements(settings, ref_seq_info, *(it->get())) )
    {
			alignments.erase(it++);
    }
    else
    {
      it++;
    }
  }
	if (alignments.size() == 0) return 0;

  // require read to be mapped! -- @JEB maybe this should be checked sooner?
	for (alignment_list::iterator it = alignments.begin(); it != alignments.end(); it++)
  {
    if (it->get()->unmapped())
    {
      alignments.erase(it);
    }
  }
  if (alignments.size() == 0) return 0;

	// @JEB v1> Unfortunately sometimes better matches don't get better alignment scores!
	// example is 30K88AAXX_LenskiSet2:1:37:1775:92 in RJW1129
	// Here a read with an extra match at the end doesn't get a higher score!!!

	// This sucks, but we need to re-sort matches and check scores ourselves...
	// for now just count mismatches (which include soft padding!)
  
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
  uint32_t best_score = mismatch_map[(alignments.front().get())];
	
  // no scores meet minimum
  
  for (alignment_list::iterator it = alignments.begin()++; it != alignments.end(); it++)
  {
    if (mismatch_map[it->get()] != best_score) break;
    last_best++;
  }

	//#broken
	//## no scores meet minimum difference between best and next best
	//#if (defined $minimum_best_score_difference && (scalar @al > $last_best+1))
	//#{
	//#	my $second_best_score = $al[$last_best+1]->aux_get("AS");
	//#	return () if ($second_best_score + $minimum_best_score_difference >= $best_score)
	//#}

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

	// Note that the score we return is higher for matches so we negative this value...
  return alignments.front()->read_length() - best_score;
}
  
//
//
//=head2 _read_alignment_passes_requirements
//
//Title   : _test_read_alignment_requirements
//Usage   : _test_read_alignment_requirements( );
//Function: Tests an individual read alignment for required match lengths
//and number of mismatches
//Returns : 
//
//=cut
//
bool _test_read_alignment_requirements(const Settings& settings, const cReferenceSequences& ref_seq_info, const alignment_wrapper& a)
{
	bool accept = true;

	if (a.unmapped()) return false;

	if (settings.required_match_length > 0)
	{
		int32_t alignment_length_on_query = a.query_match_length(); //this is the length of the alignment on the read
		if (alignment_length_on_query < settings.required_match_length)
    {
			return false;
    }
  }

	if (settings.require_complete_match)
	{
    if (!a.beginning_to_end_match())
    {
      return false; 
    }
	}
	if (settings.max_read_mismatches >= 0)
	{
		int32_t mismatches = alignment_mismatches(a, ref_seq_info);
		if (mismatches > settings.max_read_mismatches)
    {
      return false; 
    }
	}

	return true;
}

//=head2 _alignment_overlaps_junction
//
//Title   : _test_read_alignment_requirements
//Usage   : _test_read_alignment_requirements( );
//Function: Tests an individual read alignment for required match lengths
//and number of mismatches
//Returns : 
//
//=cut
//
bool _alignment_overlaps_junction(const vector<JunctionInfo>& junction_info_list, const alignment_wrapper& a)
{
  // unmapped reads don't overlap the junction
  if (a.unmapped()) return false;
  
  int32_t tid = a.reference_target_id();
  assert (tid >= 0);
  
  const JunctionInfo& this_junction_info = junction_info_list[tid];
  int32_t overlap = this_junction_info.alignment_overlap;
  
  uint32_t junction_start = this_junction_info.flanking_left + 1;
  uint32_t junction_end = this_junction_info.flanking_left + abs(overlap);

  //## If it didn't overlap the junction at all
  //## Check coordinates in the "reference" junction sequence
  if (a.reference_start_1() > junction_end) return false;
  if (a.reference_end_1() < junction_start) return false;
  return true;
}


void _write_reference_matches(const Settings& settings, cReferenceSequences& ref_seq_info, const SequenceTrimsList& trims_list, alignment_list& reference_alignments, tam_file& reference_tam, uint32_t fastq_file_index)
{
  (void)settings; //TODO: unused?
	// Nice try, no alignments
	if (reference_alignments.size() == 0) return;

	vector<Trims> trims;

  for(alignment_list::iterator it=reference_alignments.begin(); it!=reference_alignments.end(); it++)
  {
    Trims t = _trim_ambiguous_ends(**it, trims_list);
		trims.push_back(t);
  }
  
	reference_tam.write_alignments((int32_t)fastq_file_index, reference_alignments, &trims, &ref_seq_info, true);
}

bool _test_junction(const Settings& settings, Summary& summary, const string& junction_seq_id, map<string, vector<MatchedJunction> >& matched_junction_ref, map<string, map<string, MatchedJunction> >& degenerate_matches_ref, map<string, CandidateJunction>& junction_test_info_ref, cReferenceSequences& ref_seq_info, const SequenceTrimsList& trims_list, tam_file& reference_tam, tam_file& junction_tam, bool& has_non_overlap_alignment)
{
	bool verbose = false;

	if (verbose) cout << "Testing " << junction_seq_id << endl;

	// There are two kinds of matches to a candidate junction:

	// (1) Reads that uniquely map to one candidate junction (but any number of times to reference)
	vector<MatchedJunction>* unique_matches = NULL;
	if (matched_junction_ref.count(junction_seq_id))
		unique_matches = &(matched_junction_ref[junction_seq_id]);

  if (verbose) cout << "Unique size: " << (unique_matches ? unique_matches->size() : 0) << endl;

	// (2) Reads that uniquely map equally well to more than one candidate junction (and any number of times to reference)
	map<string, MatchedJunction>* degenerate_matches = NULL;
	if (degenerate_matches_ref.count(junction_seq_id))
		degenerate_matches = &(degenerate_matches_ref[junction_seq_id]);

	// FAI target id -- there is no easy way to get this short of loading the entire array and going through them...
	// Debatable about whether we save more string comparisons by doing this here or each time

	// @JEB v1> hash by tid rather than alignment junction names!!
	uint32_t junction_tid;
	for (junction_tid = 0; junction_tid < static_cast<uint32_t>(junction_tam.bam_header->n_targets); junction_tid++)
		if (junction_tam.bam_header->target_name[junction_tid] == junction_seq_id) break;

	assert(junction_tid < static_cast<uint32_t>(junction_tam.bam_header->n_targets));

	if (verbose) {
		cout << "Testing Junction Candidate: " << junction_seq_id << endl;
    size_t unique_matches_size = (unique_matches) ? unique_matches->size() : 0;
    size_t degenerate_matches_size = (degenerate_matches) ? degenerate_matches->size() : 0;
		cout << "Unique Matches: " << unique_matches_size << " Degenerate Matches: " << degenerate_matches_size << endl;
	}

	//// TEST 1: Reads that go a certain number of bp into the nonoverlap sequence on each side of the junction on each strand
	map<bool,uint32_t> max_left_per_strand = make_map<bool,uint32_t>(true,0)(false,0);
	map<bool,uint32_t> max_right_per_strand = make_map<bool,uint32_t>(true,0)(false,0);
	map<bool,uint32_t> max_min_left_per_strand = make_map<bool,uint32_t>(true,0)(false,0);
	map<bool,uint32_t> max_min_right_per_strand = make_map<bool,uint32_t>(true,0)(false,0);
	map<bool,uint32_t> count_per_strand = make_map<bool,uint32_t>(true,0)(false,0);
	uint32_t total_non_overlap_reads = 0;
	map<int32_t,bool> pos_hash[2];
  uint32_t pos_hash_count(0);

	// basic information about the junction
	JunctionInfo scj = junction_name_split(junction_seq_id);
	int32_t overlap = scj.alignment_overlap;
	int32_t flanking_left = scj.flanking_left;

	// Is there at least one read that isn't overlap only?
	// displaying ones where it doesn't as marginals looks really confusing
	has_non_overlap_alignment = false;

	// We also need to count degenerate matches b/c sometimes ambiguity unfairly penalizes real reads...
	vector<MatchedJunction*> items;
  if (unique_matches)
    for (vector<MatchedJunction>::iterator it = unique_matches->begin(); it != unique_matches->end(); it++)
      items.push_back(&(*it));
  if (degenerate_matches)
    for (map<string, MatchedJunction>::iterator it = degenerate_matches->begin(); it != degenerate_matches->end(); it++)
      items.push_back(&(it->second));

	for (uint32_t i = 0; i < items.size(); i++) // READ (loops over unique_matches, degenerate_matches)
	{
		MatchedJunction* item = items[i];
		//!!> Matches that don't extend through the overlap region will have the same quality
		//!!> as a reference match and, therefore, no difference in mapping quality
		//!!> do not count these toward scoring!
		if (item->mapping_quality_difference == 0) {
      if (verbose) cout << "  Degenerate:" << item->junction_alignments.front()->read_name() << endl;
     continue; 
    }

		total_non_overlap_reads++;
		has_non_overlap_alignment = true;

		//If there were no degenerate matches, then we could just take the
		//one and only match in the 'junction_alignments' array
		//#my $a = $item->{junction_alignments}->[0];

		// as it is, we must be sure we are looking at the one that matches
		alignment_wrapper* a = NULL;
    for(alignment_list::iterator it=item->junction_alignments.begin(); it!=item->junction_alignments.end(); it++)
		{
			alignment_wrapper* candidate_a = &(**it);
			if (candidate_a->reference_target_id() == junction_tid) {
				a = candidate_a;
				break;
			}
		}
		assert(a != NULL);

		bool rev_key = a->reversed();
		count_per_strand[rev_key]++;

		// Look at reference coords of aligned part of read sequence
    // and of unaligned ends of query continued straight to where
    // they would have aligned in the reference.
    //
    // All four of these coordinates must have never been seen before to count.
    
		int32_t begin_coord = a->reference_start_1();
    int32_t end_coord   = a->reference_end_1();
    int32_t begin_read_coord = a->reference_start_1() - (a->query_start_1() - 1);
    int32_t end_read_coord   = a->reference_end_1() + (a->read_length() - a->query_end_1());
    
    if (
           !pos_hash[rev_key].count(begin_coord)
        && !pos_hash[rev_key].count(end_coord)
        && !pos_hash[rev_key].count(begin_read_coord)
        && !pos_hash[rev_key].count(end_read_coord)
        ) 
    {
      pos_hash_count++;
    }
    
    pos_hash[rev_key][begin_coord] = true;
    pos_hash[rev_key][end_coord] = true;
    pos_hash[rev_key][begin_coord] = true;
    pos_hash[rev_key][begin_coord] = true;

		if (verbose)
			cout << "  " << item->junction_alignments.front()->read_name() << ' ' << static_cast<int32_t>(rev_key) << ' ' << begin_coord << ' ' << end_coord << ' ' << begin_read_coord << ' ' << end_read_coord << endl;

		// The left side goes exactly up to the flanking length
		uint32_t this_left = flanking_left + 1;
		this_left -= a->reference_start_1();

		// The right side starts after moving past any overlap (negative or positive)
		uint32_t this_right = flanking_left + 1;
		this_right += abs(overlap);
		this_right = a->reference_end_1() - this_right + 1;

		// Update:
		// Score = the minimum unique match length on a side
		// Max_Min = the maximum of the minimum length match sides
		// Max = the maximum match on a side
		// Note that the max and min filtering is really a kind of poor man's KS test
		//   if we implemented that with a certain coverage cutoff it would be a
		//   more principled way of doing things...
		if (this_left < this_right) {
			if (max_min_left_per_strand[rev_key] < this_left)
				max_min_left_per_strand[rev_key] = this_left;
		}
		else
		{
			if (max_min_right_per_strand[rev_key] < this_right)
				max_min_right_per_strand[rev_key] = this_right;
		}

		if (max_left_per_strand[rev_key] < this_left)
			max_left_per_strand[rev_key] = this_left;
		if (max_right_per_strand[rev_key] < this_right)
			max_right_per_strand[rev_key] = this_right;

	}

	uint32_t max_left = max(max_left_per_strand[false], max_left_per_strand[true]);
	uint32_t max_right = max(max_right_per_strand[false], max_right_per_strand[true]);

	uint32_t max_min_left = max(max_min_left_per_strand[false], max_min_left_per_strand[true]);
	uint32_t max_min_right = max(max_min_right_per_strand[false], max_min_right_per_strand[true]);

	// Save the test info about this junction.
	CandidateJunction::TestInfo test_info = {
		max_left,                           //max_left
		max_left_per_strand[false],         //max_left_minus
		max_left_per_strand[true],          //max_left_plus
		max_right,                          //max_right
		max_right_per_strand[false],        //max_right_minus
		max_right_per_strand[true],         //max_right_plus
		max_min_right,                      //max_min_right
		max_min_right_per_strand[false],    //max_min_right_minus
		max_min_right_per_strand[true],     //max_min_right_plus
		max_min_left,                       //max_min_left
		max_min_left_per_strand[false],     //max_min_left_minus
		max_min_left_per_strand[true],      //max_min_left_plus
		count_per_strand[false],            //coverage_minus
		count_per_strand[true],             //coverage_plus
		total_non_overlap_reads,            //total_non_overlap_reads
    pos_hash_count                      //pos_hash_score
	};
  
	junction_test_info_ref[junction_seq_id].test_info = test_info;
	junction_test_info_ref[junction_seq_id].pos_hash_score = pos_hash_count;

	// Old way, requiring certain overlap on each side on each strand
	// @JEB !> Best results may be to combine these methods

	// These parameters still need additional testing
	// and, naturally, they have problems with scaling with the
	// total number of reads...

	uint32_t alignment_on_each_side_cutoff = 14; //16
	uint32_t alignment_on_each_side_cutoff_per_strand = 9; //13
	uint32_t alignment_on_each_side_min_cutoff = 3;

	bool failed = (max_left < alignment_on_each_side_cutoff)
				|| (max_right < alignment_on_each_side_cutoff)
        || (max_left_per_strand[false] < alignment_on_each_side_cutoff_per_strand)
				|| (max_left_per_strand[true] < alignment_on_each_side_cutoff_per_strand)
	      || (max_right_per_strand[false] < alignment_on_each_side_cutoff_per_strand)
				|| (max_right_per_strand[true] < alignment_on_each_side_cutoff_per_strand)
				|| (max_min_left < alignment_on_each_side_min_cutoff)
				|| (max_min_right < alignment_on_each_side_min_cutoff)
	;

	// POS_HASH test
	// New way, but we need to have examined the coverage distribution to calibrate what scores to accept!
	uint32_t junction_accept_score_cutoff_1 = summary.preprocess_coverage[scj.sides[0].seq_id].junction_accept_score_cutoff;
	uint32_t junction_accept_score_cutoff_2 = summary.preprocess_coverage[scj.sides[1].seq_id].junction_accept_score_cutoff;
  
  // both score cutoffs might be zero - indicating these are missing contigs that are basically deleted
  // fail if this is the case. Revisit this logic at a future time. @JEB
  failed = failed || ((junction_accept_score_cutoff_1 == 0) && (junction_accept_score_cutoff_2 == 0)) ;
  
	failed = failed || ( ( test_info.pos_hash_score < junction_accept_score_cutoff_1 ) && ( test_info.pos_hash_score < junction_accept_score_cutoff_2 ) );
	if (verbose) cout << (failed ? "Failed" : "Passed") << endl;

	// TODO:
	// ADD -- NEED TO CORRECT OVERLAP AND ADJUST NUMBER OF READS SUPPORTING HERE, RATHER THAN LATER
	//

	// DEGENERATE JUNCTION MATCHES
	// ===========================
	// Determine the fate of degenerate reads that map to this junction

	if (degenerate_matches)
	{
		for (map<string, MatchedJunction>::iterator it = degenerate_matches->begin(); it != degenerate_matches->end(); it++)
		{      
			MatchedJunction& degenerate_match = it->second;
			uint32_t fastq_file_index = degenerate_match.fastq_file_index;

			// Success for this candidate junction...
			// purge all references to this from the degenerate match hash
			// so that they will not be counted for other junctions
			if (!failed)
			{
        if (verbose) cout << "Before size: " << (unique_matches ? unique_matches->size() : 0) << endl;
        
				// Purge all references to this read from the degenerate match hash
        // so that it cannot be counted for any other junction
        
        counted_ptr<bam_alignment> matched_alignment(NULL);
        for (alignment_list::iterator it2=degenerate_match.junction_alignments.begin(); it2 !=degenerate_match.junction_alignments.end(); )
				{          
          // we make a copy and then increment, in case the current iterator value will be erased
					counted_ptr<bam_alignment> a = *it2; it2++;
          string test_junction_seq_id = junction_tam.target_name(*a);
          
          //this is the one for the current candidate junction
          if (a->reference_target_id() == junction_tid)
          {
            matched_alignment = a;
          }
          else
          {
            size_t deleted = degenerate_matches_ref[test_junction_seq_id].erase(a->read_name());
          }
          
          if (degenerate_matches_ref[test_junction_seq_id].size() == 0)
          {
            degenerate_matches_ref.erase(test_junction_seq_id);
          }
        }

				assert(matched_alignment.get() != NULL);
				degenerate_match.junction_alignments.clear();
        degenerate_match.junction_alignments.push_back(matched_alignment);
        
        // We need to add this degenerately matched read to the other ones supporting this junction
        // Create empty list if necessary...
        if (matched_junction_ref.count(junction_seq_id) == 0) {
          matched_junction_ref.insert( pair<string, vector<MatchedJunction> >(junction_seq_id, vector<MatchedJunction>()) );
          unique_matches = &(matched_junction_ref[junction_seq_id]);
        }
				unique_matches->push_back(degenerate_match);
        
        if (verbose) cout << "After size: " << (unique_matches ? unique_matches->size() : 0) << endl;

			}

			// Failure for this candidate junction...
			// Remove just the degenerate hits to this candidate junction
			// Once all have failed, then we need to add the reference alignments (if any)!
			else
			{
				degenerate_match.degenerate_count--;

				if (verbose) cout << "New Degenerate match count: " << degenerate_match.degenerate_count << endl;

				// This degenerate match missed on all opportunities,
				// we should add it to the reference sequence
				if (degenerate_match.degenerate_count == 0)
				{
					alignment_list& this_reference_al = degenerate_match.reference_alignments;
					_write_reference_matches(settings, ref_seq_info, trims_list, this_reference_al, reference_tam, fastq_file_index);
				}
        
        counted_ptr<bam_alignment> matched_alignment(NULL);
        for (alignment_list::iterator it2=degenerate_match.junction_alignments.begin(); it2 !=degenerate_match.junction_alignments.end(); it2++)
				{
					counted_ptr<bam_alignment>& candidate_a = *it2; //this is the one for the current candidate junction
					if (candidate_a->reference_target_id() == junction_tid)
          {
						matched_alignment = candidate_a;
            break;
          }
        }
        
        // Write alignment to SAM file for candidate junctions regardless of success...
        // Note that successful ones get written below, because they were pushed to the other list
        assert(matched_alignment.get() != NULL);
        if (has_non_overlap_alignment) {
          alignment_list alignments;
          alignments.push_back(matched_alignment);
          junction_tam.write_alignments(fastq_file_index, alignments, NULL, &ref_seq_info, true);
        }
			}
		}

		// We are completely done with degenerate matches to this junction id.
		// Deleting them here means that we will never go through this loop with them again
		// and is necessary for not doubly writing them.
		degenerate_matches_ref.erase(junction_seq_id);
	}

	// UNIQUE JUNCTION MATCHES
	// =======================
  // If there were no unique matches to begin with, we may have created this entry...

  
  if (unique_matches)
  {
    if (verbose) cout << "Printing size:" << (unique_matches ? unique_matches->size() : 0) << endl;

    for (uint32_t i = 0; i< unique_matches->size(); i++)
    {
      MatchedJunction& item = (*unique_matches)[i];
      // Write out the matches to the proper SAM file(s) depending on whether the junction succeeded or failed
      uint32_t fastq_file_index = item.fastq_file_index;

      // ONLY if we failed: write matches to reference sequences
      if (failed)
      {
        alignment_list this_reference_al = item.reference_alignments;
        _write_reference_matches(settings, ref_seq_info, trims_list, this_reference_al, reference_tam, fastq_file_index);
      }

      // REGARDLESS of success: write matches to the candidate junction SAM file
      //if (has_non_overlap_alignment) - since it passed, we want to show them all
        junction_tam.write_alignments(fastq_file_index, item.junction_alignments, NULL, &ref_seq_info, true);
    }
  }
	return !failed;
}

diff_entry _junction_to_hybrid_list_item(const string& key, cReferenceSequences& ref_seq_info, uint32_t total_reads, CandidateJunction& test_info)
{
  
  CandidateJunction::TestInfo& this_test_info = test_info.test_info;
  
	// split the key to an item with information about the junction
	JunctionInfo jc = junction_name_split(key);

	jc.key = key;

	// overlap may be adjusted below... this messes up making the alignment
	// 'alignment_overlap' is the original one that applies to the candidate junction BAM file
	// 'overlap' is a version where overlap has been resolved if possible for adding sides of the
	//    alignment
	jc.overlap = jc.alignment_overlap;
	jc.total_reads = total_reads;

	// Redundancy is loaded from the key, but we doubly enforce it when IS elements are involved.

	// Correct for overlapping IS elements

	///
	// IS insertion overlap correction
	//
	// For these the coordinates may have been offset incorrectly initially (because both sides of the junction may look unique)
	// The goal is to offset through positive overlap to get as close as possible to the ends of the IS
	///

	cSequenceFeature* is = NULL;
	for (int32_t i = 0; i <= 1; i++)
	{
		// Determine IS elements
		// Is it within an IS or near the boundary of an IS in the direction leading up to the junction?
		is = cReferenceSequences::find_closest_repeat_region(jc.sides[i].position, ref_seq_info.repeat_lists[jc.sides[i].seq_id], 200, jc.sides[i].strand);
		if (is != NULL)
		{
			jc.sides[i].is.name = is->SafeGet("name");
			jc.sides[i].is.interval = (is->m_strand == 1) ? to_string(is->m_start) + "-" + to_string(is->m_end) : to_string(is->m_end) + "-" + to_string(is->m_start);
			jc.sides[i].is.product = is->SafeGet("product");
		}
	}

	//_add_is_coords_from_interval($jc->{side_1});
	//_add_is_coords_from_interval($jc->{side_2});
	for (int32_t i = 0; i <= 1; i++)
	{
		JunctionInfo::Side& c = jc.sides[i];
    
		//return if (!defined $c->{is});
    if (c.is.interval.size()!=0)
    {
      vector<string> is_start_end = split(c.is.interval, "-");
      int32_t is_start = from_string<int32_t>(is_start_end[0]);
      int32_t is_end = from_string<int32_t>(is_start_end[1]);
      c.is.strand = (is_start < is_end) ? +1 : -1;
      c.is.start = c.is.strand ? is_start : is_end;
      c.is.end = c.is.strand ? is_end : is_start;
    }
	}

	jc.sides[0].read_side = false;
	jc.sides[1].read_side = true;

	// Determine which side of the junction is the IS and which is unique
	// these point to the correct initial interval...
	jc.is_side = UNDEFINED_UINT32;
	if (jc.sides[0].is.name.size() > 0 && jc.sides[1].is.name.size() == 0)
	{
		if (abs(static_cast<int32_t>(jc.sides[0].is.start) - static_cast<int32_t>(jc.sides[0].position)) <= 20)
		{
			jc.is_side = 0;
			jc.sides[jc.is_side].is.side_key = "start";
		}
		else if (abs(static_cast<int32_t>(jc.sides[0].is.end) - static_cast<int32_t>(jc.sides[0].position)) <= 20 )
		{
			jc.is_side = 0;
			jc.sides[jc.is_side].is.side_key = "end";
		}
		jc.unique_side = 1;
	}

	else if (jc.sides[0].is.name.size() == 0 && jc.sides[1].is.name.size() > 0)
	{
		if (abs(static_cast<int32_t>(jc.sides[1].is.start) - static_cast<int32_t>(jc.sides[1].position)) <= 20)
		{
			jc.is_side = 1;
			jc.sides[jc.is_side].is.side_key = "start";
		}
		else if (abs(static_cast<int32_t>(jc.sides[1].is.end) - static_cast<int32_t>(jc.sides[1].position)) <= 20 )
		{
			jc.is_side = 1;
			jc.sides[jc.is_side].is.side_key = "end";
		}
		jc.unique_side = 0;
	}
	// both were IS! -- define as redundant here
	else if (jc.sides[0].is.name.size() > 0)
		jc.sides[0].redundant = true;
	else if (jc.sides[1].is.name.size() > 0)
		jc.sides[1].redundant = true;

	// By default, overlap is included on both sides of the junction (possibly changed below)
	jc.sides[0].overlap = 0;
	jc.sides[1].overlap = 0;

	// Resolve redundant overlap
	if (jc.overlap > 0)
	{
		jc.sides[0].overlap = jc.overlap;
		jc.sides[1].overlap = jc.overlap;

		// If there was in IS, resolve overlap so it goes to the edge of the IS element
		if (jc.is_side != UNDEFINED_UINT32)
		{
			// first, adjust the repetitive sequence boundary to get as close to the IS as possible
      assert(jc.sides[jc.is_side].is.side_key.size() > 0);
			int32_t move_dist = jc.sides[jc.is_side].strand * (static_cast<int32_t>((jc.sides[jc.is_side].is.side_key == "start" 
          ? jc.sides[jc.is_side].is.start : jc.sides[jc.is_side].is.end)) - jc.sides[jc.is_side].position);

			if (move_dist < 0) move_dist = 0;
			if (move_dist > jc.overlap) move_dist = jc.overlap ;

			jc.sides[jc.is_side].position += jc.sides[jc.is_side].strand * move_dist;
			jc.overlap -= move_dist;
			jc.sides[jc.is_side].overlap -= move_dist;

			// second, adjust the unique sequence side with any remaining overlap
			jc.sides[jc.unique_side].position += jc.sides[jc.unique_side].strand * jc.overlap;
			jc.sides[jc.unique_side].overlap -= jc.overlap;

			jc.overlap = 0;
		}
		/// If there is no IS element and
		//    (1) both sides are unique
		// OR (2) only the second side is redundant,
		// OR (3) both sides are redundant
		/// then give overlap to first side.
		/// This gives proper support for junctions.
		/// and ensures we don't count this coverage twice.
		else if ((!jc.sides[0].redundant) || (jc.sides[0].redundant && jc.sides[1].redundant) )
		{
			uint32_t strand_direction = (jc.sides[1].strand > 0 ? 1 : -1);
			jc.sides[1].position += jc.overlap * strand_direction;
			jc.sides[1].overlap = 0;
			jc.overlap = 0;
		}
		else  // side_1 was redundant, give overlap to side_2
		{
			uint32_t strand_direction = (jc.sides[0].strand > 0 ? -1 : 1);
			jc.sides[0].position += jc.overlap * strand_direction;
			jc.sides[0].overlap = 0;
			jc.overlap = 0;
		}

		// If both sides were redundant, no adjustment because we are not going to count coverage
	}

	// flatten things to only what we want to keep
	//TODO: Are parameters to constructor correct?
	diff_entry item("JC");
	item
		("side_1_seq_id", jc.sides[0].seq_id)
		("side_1_position", to_string(jc.sides[0].position))
		("side_1_redundant", to_string(jc.sides[0].redundant))
		("side_1_strand", to_string(jc.sides[0].strand))
		("side_1_overlap", to_string(jc.sides[0].overlap))

		("side_2_seq_id", jc.sides[1].seq_id)
		("side_2_position", to_string(jc.sides[1].position))
		("side_2_redundant", to_string(jc.sides[1].redundant))
		("side_2_strand", to_string(jc.sides[1].strand))
		("side_2_overlap", to_string(jc.sides[1].overlap))

		("key", jc.key)
		("alignment_overlap", to_string(jc.alignment_overlap))
		("overlap", to_string(jc.overlap))
		("total_reads", to_string(jc.total_reads))
		("flanking_left", to_string(jc.flanking_left))
		("flanking_right", to_string(jc.flanking_right))

		("unique_read_sequence", to_string(jc.unique_read_sequence))
	;
  
//	## may want to take only selected of these fields..
  
  item
  ("max_left", to_string(this_test_info.max_left))
  ("max_left_minus", to_string(this_test_info.max_left_minus))
  ("max_left_plus", to_string(this_test_info.max_left_plus))
  ("max_right", to_string(this_test_info.max_right))
  ("max_right_minus", to_string(this_test_info.max_right_minus))
  ("max_right_plus", to_string(this_test_info.max_right_plus))
  ("max_min_right", to_string(this_test_info.max_min_right))
  ("max_min_right_minus", to_string(this_test_info.max_min_right_minus))
  ("max_min_right_plus", to_string(this_test_info.max_min_right_plus))
  ("max_min_left", to_string(this_test_info.max_min_left))
  ("max_min_left_minus", to_string(this_test_info.max_min_left_minus))
  ("max_min_left_plus", to_string(this_test_info.max_min_left_plus))
  ("coverage_minus", to_string(this_test_info.coverage_minus))
  ("coverage_plus", to_string(this_test_info.coverage_plus))
  ("total_non_overlap_reads", to_string(this_test_info.total_non_overlap_reads))
  ("pos_hash_score", to_string(this_test_info.pos_hash_score))
  ;

	/// Note: Other adjustments to overlap can happen at the later annotation stage
	/// and they will not affect coverage for calling deletions or mutations
	/// because they will be in REDUNDANTLY matched sides of junctions
	return item;
}

  
Trims _trim_ambiguous_ends(const alignment_wrapper& a, const SequenceTrimsList& trims)
{
	bool verbose = false;
  
	// which reference sequence?
	uint32_t tid = a.reference_target_id();
  
  Trims t;
  t.L= trims[tid].left_trim_0(a.reference_start_0());
  t.R = trims[tid].right_trim_0(a.reference_end_0());

  t.L += a.query_start_0();
  t.R += a.read_length() - a.query_end_1();

//  cerr << a.read_name() << endl;
//  cerr << "start: " << a.reference_start_1() << " end: " << a.reference_end_1() << endl;
//  cerr << "left: " << t.L << " right: " << t.R << endl;

  return t;
}
  
void read_trims(SequenceTrimsList& trims, const cReferenceSequences& ref_seqs, const string &in_trims_file_name ) 
{
  trims.resize(ref_seqs.size());
  for(uint32_t i = 0; i < ref_seqs.size(); i++) {
    string this_file_name = Settings::file_name(in_trims_file_name, "@", ref_seqs[i].m_seq_id);
    trims[i].ReadFile(this_file_name, ref_seqs[i].m_length);
  }
}


//sort junction ids based on size of vector contained in map
vector<string> get_sorted_junction_ids(map<string, vector<MatchedJunction> >& unique_map, map<string, map<string, MatchedJunction> >& degenerate_map, const vector<string>& keys)
{
  vector<VectorSize> vector_sizes;
  for (uint32_t i = 0; i < keys.size(); i++)
  {
    // may or may not exist
    uint32_t degenerate_count = 0;
    if (degenerate_map.count(keys[i]))
    {
      degenerate_count = degenerate_map[keys[i]].size();
    }
    
    VectorSize info(keys[i], unique_map[keys[i]].size()+degenerate_count, degenerate_count);
    vector_sizes.push_back(info);
  }
  sort(vector_sizes.begin(), vector_sizes.end(), VectorSize::sort_reverse_by_size);
  
  vector<string> sorted_junction_ids;
  for (uint32_t i = 0; i < keys.size(); i++)
    sorted_junction_ids.push_back(vector_sizes[i].junction_id);
  return sorted_junction_ids;
}
  
} // namespace breseq

