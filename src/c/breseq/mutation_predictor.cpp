
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

#include "breseq/mutation_predictor.h"

#include "breseq/alignment.h"
#include "breseq/annotated_sequence.h"
#include "breseq/genome_diff.h"
#include "breseq/settings.h"
#include "breseq/fastq.h"


using namespace std;

namespace breseq {

  cReferenceSequences MutationPredictor::ref_seq_info;
  
	MutationPredictor::MutationPredictor(cReferenceSequences& _ref_seq_info)
	{
		 ref_seq_info = _ref_seq_info;
	}

	// Private methods

	// Utility function
	string MutationPredictor::get_sequence(string seq_id, uint32_t start, uint32_t end)
	{
		bool verbose = false;
		if (verbose)
			cout << "Get sequence: " << seq_id << ":" << start << "-" << end << endl;
		return ref_seq_info.ref_strings[seq_id].substr(start - 1, end - start + 1);
	}

	cSequenceFeature* MutationPredictor::within_repeat(string seq_id, uint32_t position)
	{
		vector<cSequenceFeature>& r = ref_seq_info.repeat_lists[seq_id];
    
		for (uint32_t i = 0; i < r.size(); i++)
			if ((r[i].m_start <= position) && (position <= r[i].m_end))
				return &(r[i]);

		return NULL;
	}

	bool MutationPredictor::sort_by_hybrid(const counted_ptr<diff_entry>& a, const counted_ptr<diff_entry>& b)
	{
    
		//return (a.field < b.field);
		int32_t a_pos = n(a->entry_exists("_side_1_is") ? (*a)["side_2_position"] : (*a)["side_1_position"]);
		int32_t b_pos = n(b->entry_exists("_side_1_is") ? (*b)["side_2_position"] : (*b)["side_1_position"]);

		int32_t a_seq_order = (a->entry_exists("_side_1_is") ? ref_seq_info.seq_order[(*a)["side_2_seq_id"]] : ref_seq_info.seq_order[(*a)["side_1_seq_id"]]);
		int32_t b_seq_order = (b->entry_exists("_side_1_is") ? ref_seq_info.seq_order[(*b)["side_2_seq_id"]] : ref_seq_info.seq_order[(*b)["side_1_seq_id"]]);

		uint32_t a_reject_order = number_reject_reasons(*a);
		uint32_t b_reject_order = number_reject_reasons(*b);

		// sort by seq_id, position, fewer reject reasons, then score (highest to lowest)
		return (
			(a_seq_order != b_seq_order) ? (a_seq_order < b_seq_order) :
			(a_pos != b_pos) ? (a_pos < b_pos) :
			(a_reject_order != b_reject_order) ? (a_reject_order < b_reject_order) :
			((*a)["pos_hash_score"] != (*b)["pos_hash_score"]) ? (n((*a)["pos_hash_score"]) > n((*b)["pos_hash_score"])) :
			(n((*a)["min_overlap_score"]) > n((*b)["min_overlap_score"]))
		);
	}

	bool MutationPredictor::sort_by_reject_score(const counted_ptr<diff_entry>& a, const counted_ptr<diff_entry>& b)
	{
		uint32_t a_reject_order = number_reject_reasons(*a);
		uint32_t b_reject_order = number_reject_reasons(*b);

		// sort by seq_id, position, fewer reject reasons, then score (highest to lowest)
		return (
			(a_reject_order != b_reject_order) ? (a_reject_order < b_reject_order) :
			(n((*a)["score"]) < n((*b)["score"]))
		);
	}

	// look at SNPs and small indels predicted by read alignments.
	// be sure they are sorted by position
	bool MutationPredictor::sort_by_pos(const counted_ptr<diff_entry>& a, const counted_ptr<diff_entry>& b)
	{
		return (
			((*a)["seq_id"] != (*b)["seq_id"]) ? ((*a)["seq_id"] < (*b)["seq_id"]) :
			((*a)["position"] != (*b)["position"]) ? (n((*a)["position"]) < n((*b)["position"])) :
			(n((*a)["insert_position"]) < n((*b)["insert_position"]))
		);
	}

	/*
	 Title   : predict
	 Usage   : $mp->predict();
	 Function: Predicts mutations from evidence in a GenomeDiff and adds them to it
	 Returns :

	*/
	void MutationPredictor::predict(Settings& settings, genome_diff& gd, uint32_t max_read_length)
	{
		
		bool verbose = false;

		///
		//  Preprocessing of JC evidence
		///

		// For all that follows, we need information about repeat_regions overlapping the sides of junctions
		vector<string> jc_types = make_list<string>("JC");
		genome_diff::entry_list_t jc = gd.list(jc_types);

		for (uint32_t i = 0; i < jc.size(); i++)
		{
			diff_entry& j = *(jc[i].get());
      
			j["_side_1_read_side"] = "false";
			j["_side_2_read_side"] = "true";

			for (uint32_t side = 1; side <= 2; side++)
			{
				string side_key = "side_" + s(side);

				cSequenceFeature* is = ref_seq_info.find_closest_repeat_region(
					n(j[side_key + "_position"]),
					ref_seq_info.repeat_lists[j[side_key + "_seq_id"]],
					50,
					n(j[side_key + "_strand"])
				);
				if (is != NULL)
				{
					j["_" + side_key + "_is"] = "1";
					j["_" + side_key + "_is_start"] = s(is->m_start);
					j["_" + side_key + "_is_end"] = s(is->m_end);
          j["_" + side_key + "_is_name"] = (*is)["name"];
          j["_" + side_key + "_is_strand"] = s(is->m_strand);
				}
				
        bool test_val = b(j[side_key + "_redundant"]);
				j[side_key + "_annotate_key"] = (j.entry_exists("_" + side_key + "_is_start") || b(j[side_key + "_redundant"])) ? "repeat" : "gene";
			}

			// by default, we are sorted by this coord
			j["_unique_interval"] = "side_1";

			// Determine which side of the junction is the IS and which is unique
			// these point to the correct initial interval...
			if (j.entry_exists("_side_1_is"))
			{
        //cout << n(j["_side_1_is_start"]) << " " << n(j["_side_1_is_end"]) << " " << n(j["side_1_position"]) << endl;
        
				if (abs(n(j["_side_1_is_start"]) - n(j["side_1_position"])) <= 20)
				{
					j["_is_interval"] = "side_1";
					j["_is_interval_closest_side_key"] = "start";
					j["_unique_interval"] = "side_2";
				}
				else if (abs(n(j["_side_1_is_end"]) - n(j["side_1_position"])) <= 20)
				{
					j["_is_interval"] = "side_1";
					j["_is_interval_closest_side_key"] = "end";
					j["_unique_interval"] = "side_2";
				}
			}
      
      if (j.entry_exists("_side_2_is"))
			{
        //cout << n(j["_side_2_is_start"]) << " " << n(j["_side_2_is_end"]) << " " << n(j["side_2_position"]) << endl;

				if (abs(n(j["_side_2_is_start"]) - n(j["side_2_position"])) <= 20)
				{
					j["_is_interval"] = "side_2";
					j["_is_interval_closest_side_key"] = "start";
					j["_unique_interval"] = "side_1";
				}
				else if (abs(n(j["_side_2_is_end"]) - n(j["side_2_position"])) <= 20)
				{
					j["_is_interval"] = "side_2";
					j["_is_interval_closest_side_key"] = "end";
					j["_unique_interval"] = "side_1";
				}
			}
      
      // add in the rest of the unique interval information
      j["_unique_interval_seq_id"] = j[j["_unique_interval"] + "_seq_id"];
      j["_unique_interval_position"] = j[j["_unique_interval"] + "_position"];
      j["_unique_interval_strand"] = j[j["_unique_interval"] + "_strand"];
    }
    
    // Don't count rejected ones, this can be relaxed, but it makes MOB prediction much more complicated and prone to errors.
		for(genome_diff::entry_list_t::iterator it = jc.begin(); it < jc.end(); it++)
    {
      diff_entry& de = **it;
			if (de.entry_exists("reject"))
      {
				jc.erase(it);
        it--;
      }
		}
    
		vector<string> mc_types = make_list<string>("MC");
		genome_diff::entry_list_t mc = gd.list(mc_types);

		///
		// evidence MC + JC => DEL mutation
		///

		// DEL prediction:
		// (1) there is a junction that exactly crosses the deletion boundary deletion
		// (2) there is is no junction, but both ends of the deletion are in repeat sequences
		// (3) there is a junction between unique sequence and a repeat element

		for (uint32_t i = 0; i < mc.size(); i++) //MC
		{
      diff_entry& mc_item = *mc[i];
      
			if (mc_item.entry_exists("reject"))
			  continue;

			// set up generic deletion item
			diff_entry mut;
			mut._type = "DEL";
			mut._evidence = make_list<string>(mc_item._id);
			mut
				("seq_id", mc_item["seq_id"])
				("position", mc_item["start"])
				("size", s(n(mc_item["end"]) - n(mc_item["start"]) + 1));
			;

			///
			// (1) there is a junction that exactly crosses the deletion boundary deletion
			///

			for(genome_diff::entry_list_t::iterator it = jc.begin(); it < jc.end(); it++) //JC
			{
				diff_entry& jc_item = **it;

				if (jc_item["side_1_seq_id"] != mut["seq_id"] || jc_item["side_2_seq_id"] != mut["seq_id"])
					continue;

				// We always know that the lower coordinate part of the junction is first, hence these
				// assumptions about the strands hold.
				if (
					   (n(jc_item["side_1_position"]) == n(mut["position"])-1)
					&& (n(jc_item["side_1_strand"]) == -1)
					&& (n(jc_item["side_2_position"]) == n(mut["position"])+n(mut["size"]))
					&& (n(jc_item["side_2_strand"]) == +1)
           )
				{
					mut._evidence.push_back(jc_item._id);
					jc.erase(it);
					gd.add(mut);
					continue;
				}
			}

			///
			// (2) there is is no junction, but both ends of the deletion are in repeat sequences
			///

			cSequenceFeature* r1_pointer = within_repeat(mut["seq_id"], n(mut["position"]));
			cSequenceFeature* r2_pointer = within_repeat(mut["seq_id"], n(mut["position"]) + n(mut["size"]));

			// Then we will adjust the coordinates to remove...
			if (r1_pointer != NULL && r2_pointer != NULL && ((*r1_pointer)["name"] == (*r2_pointer)["name"]))
			{
				cSequenceFeature& r1 = *r1_pointer, r2 = *r2_pointer;

				// there may be more evidence that one or the other is deleted...
				uint32_t r1_overlap_end = n(mc_item["start"]) + n(mc_item["start_range"]);
				if (r1_overlap_end > r1.m_end)
					r1_overlap_end = r1.m_end;
				int32_t r1_overlap = r1_overlap_end - n(mc_item["start"]) + 1;

				uint32_t r2_overlap_start = n(mc_item["end"]) - n(mc_item["end_range"]);
				if (r2_overlap_start < r1.m_start)
					r2_overlap_start = r2.m_start;
				int32_t r2_overlap = n(mc_item["end"]) - r2_overlap_start + 1;

				// it may be really close...defined by read length of genome in which case
				uint32_t slop_distance = max_read_length;

				// prefer to delete the second copy
				if ( (static_cast<uint32_t>(abs(r1_overlap - r2_overlap)) <= slop_distance) || (r2_overlap > r1_overlap ))
				{
					mut["position"] = s(r1.m_end + 1);
					mut["size"] = s(r2.m_end - r1.m_end);
				}
				else // delete the first copy
				{
					mut["position"] = s(r1.m_start);
					mut["size"] = s(r2.m_start - r1.m_start);
				}

				// remember the name of the element
				mut["between"] = r1["name"];
				gd.add(mut);
				continue;
			}

			// Both sides were unique or redundant, nothing more we can do...
			if ( (r1_pointer == NULL && r2_pointer == NULL) || (r1_pointer != NULL && r2_pointer != NULL) )
				continue;

			///
			// (3) there is a junction between unique sequence and a repeat element
			///
			cSequenceFeature& r = (r1_pointer != NULL) ? *r1_pointer : *r2_pointer;
			int32_t redundant_deletion_side = (r1_pointer != NULL) ? -1 : +1;
			int32_t unique_deletion_strand = -redundant_deletion_side;
			int32_t needed_coord = (r1_pointer != NULL)
			  ? n(mut["position"]) + n(mut["size"])
			  : n(mut["position"]) - 1;

			bool ok_were_good = false;
			for(genome_diff::entry_list_t::iterator it = jc.begin(); it < jc.end(); it++) //JUNCTION
			{
				diff_entry& j = **it;

				if (!j.entry_exists("_is_interval")) continue;

				if (verbose)
					cout << "Check 1: " << j[j["_unique_interval"] + "_seq_id"] << " ne " << mut["seq_id"] << endl;
				if (j[j["_unique_interval"] + "_seq_id"] != mut["seq_id"])
					continue;
				if (verbose)
					cout << "Pass 1" << endl;

				// check type of IS
				if (verbose)
					cout << "Check 2: " << r["name"] << " ne " << j["_" + j["_is_interval"] + "_is_name"] << endl;
				if (r["name"] != j["_" + j["_is_interval"] + "_is_name"])
					continue;

				if (verbose) cout << "Pass 2" << endl;

				// check that the unique side matches coordinate
				if (verbose)
					cout << "Check 3: " << j[j["_unique_interval"] + "_position"] << " != " << needed_coord << endl;
				if (n(j[j["_unique_interval"] + "_position"]) != needed_coord)
					continue;
				if (verbose)
					cout << "Pass 3" << endl;

				/*print Dumper($mut) if ($verbose);
				print Dumper($j) if ($verbose);
				print Dumper($r) if ($verbose);*/

				// check that IS is on the right strand
				if (verbose)
					cout << "Check 4: " << redundant_deletion_side << " * " << r.m_strand << " != " << j[j["_is_interval"] + "_strand"] << " * " << j["_" + j["_is_interval"] + "_is_strand"] << endl;
				if ( (redundant_deletion_side * r.m_strand) != (n(j[j["_is_interval"] + "_strand"]) * n(j["_" + j["_is_interval"] + "_is_strand"])) )
					continue;
				if (verbose)
					cout << "Pass 4" << endl;

				// check that the unique side is on the right strand
				if (verbose)
					cout << "Check 5: " << unique_deletion_strand << " != " << j[j["_unique_interval"] + "_strand"] << endl;
				if ( unique_deletion_strand != n(j[j["_unique_interval"] + "_strand"]) )
					continue;
				if (verbose)
					cout << "Pass 5" << endl;

				// need to adjust the non-unique coords
				if (!redundant_deletion_side)
				{
					uint32_t move_dist = r.m_end + 1 - n(mut["position"]);
					mut["position"] = s(n(mut["position"]) + move_dist);
					mut["size"] = s(n(mut["size"]) - move_dist);
				}
				else
				{
					uint32_t move_dist = (n(mut["position"]) + n(mut["size"]) - 1) - (r.m_start - 1);
					mut["size"] = s(n(mut["size"]) - move_dist);
				}

				// OK, we're good!
				mut["mediated"] = r["name"];
				mut._evidence.push_back(j._id);
				jc.erase(it);
				gd.add(mut);
				ok_were_good = true;
			}
			if (ok_were_good) continue;
		}


		///
		// evidence JC + JC = MOB mutation
		///

		for(genome_diff::entry_list_t::iterator it = jc.begin(); it < jc.end(); it++) //JC
		{
			diff_entry& j = **it;
		
			// Ah, we don't have an IS, we are done
			if (!j.entry_exists("_is_interval")) continue;

			// Ah, there is no overlap to play with, we are done
			if (n(j["overlap"]) <= 0) continue;

			// The following code implies $j->{overlap} > 0

			/// first, adjust the repetitive sequence boundary to get as close to the IS as possible
			int32_t is_interval_position = n(j[j["_is_interval"] + "_position"]);
			uint32_t move_dist = abs(is_interval_position - n(j[j["_is_interval"] + "_is_" + j["_is_interval_closest_side_key"]])); //  + "_position" ?
			uint32_t overlap = n(j["overlap"]);
			if (move_dist > overlap) move_dist = overlap;

			is_interval_position += (b(j[j["_is_interval"] + "_strand"]) ? 1 : -1) * move_dist;
			j[j["_is_interval"] + "_position"] = s(is_interval_position);
			overlap -= move_dist;

			/// second, adjust the unique sequence side with any remaining overlap
			int32_t unique_interval_position = n(j[j["_unique_interval"] + "_position"]);
			unique_interval_position += (b(j[j["_unique_interval"] + "_strand"]) ? 1 : -1) * overlap;
			j[j["_unique_interval"] + "_position"] = s(unique_interval_position);

			j["overlap"] = "0";
		}

		sort(jc.begin(), jc.end(), MutationPredictor::sort_by_hybrid);


		for (uint32_t i = 0; i < jc.size(); i++) //JC1
		{
			diff_entry& j1 = *jc[i];

			// Compile a list of the next possibilities within a certain length of bases
      genome_diff::entry_list_t j2_list;
			
			for (uint32_t j = 1; i + j < jc.size(); j++) //JC2
			{
				diff_entry& j2 = *jc[i+j];
        
				// must be close together in real coords
				if ( (j1[j1["_unique_interval"] + "_seq_id"] != j2[j2["_unique_interval"] + "_seq_id"])
					|| (abs(n(j1[j1["_unique_interval"] + "_position"]) - n(j2[j2["_unique_interval"] + "_position"])) > 20 ) )
					break;

				if ( (!j1.entry_exists("_is_interval") || !j2.entry_exists("_is_interval"))
					|| (j1["_" + j1["_is_interval"] + "_is_name"] != j2["_" + j2["_is_interval"] + "_is_name"]) )
					continue;

				j2["_delete_index"] = s(i + j); // for remembering what to delete if this one succeeds
				j2_list.push_back(jc[i+j]);

				j++;
			}

			//sort the $j2_list by reject reason and score

			sort(j2_list.begin(), j2_list.end(), sort_by_reject_score);

			// We need to go through all with the same coordinate (or within a certain coordinate stretch?)
			// because sometimes a failed junction will be in between the successful junctions
			for(genome_diff::entry_list_t::iterator it = j2_list.begin(); it < j2_list.end(); it++) //J2
			{
				diff_entry& j2 = **it;

				// positive overlap should be resolved by now
				assert(n(j1["overlap"]) <= 0);
				assert(n(j2["overlap"]) <= 0);

				// the first unique coords are going into the IS element
				bool uc1_strand = b(j1[j1["_unique_interval"] + "_strand"]);
				bool uc2_strand = b(j2[j2["_unique_interval"] + "_strand"]);
				if (uc1_strand == uc2_strand) continue;

				// What strand is the IS on relative to the top strand of the genome
				bool is1_strand = ! (b(j1[j1["_is_interval"] + "_strand"]) == b(j1["_" + j1["_is_interval"] + "_is_strand"]) == b(j1[j1["_unique_interval"] + "_strand"]));
				bool is2_strand = ! (b(j2[j2["_is_interval"] + "_strand"]) == b(j2["_" + j2["_is_interval"] + "_is_strand"]) == b(j2[j2["_unique_interval"] + "_strand"]));

				// Remove these predictions from the list, $j2 first, so indices don't shift
				jc.erase(jc.begin() + n(j2["_delete_index"]));
				jc.erase(jc.begin() + i);
				i--; // only minus one b/c the current one was deleted

				// Create the mutation, with evidence
				diff_entry mut;
				mut._type = "MOB";
				mut._evidence.push_back(j1._id);
				mut._evidence.push_back(j2._id);
				mut
					("seq_id", j1[j1["_unique_interval"] + "_seq_id"])
				;
				mut["_start"] = (!uc1_strand) ? j2[j2["_unique_interval"] + "_position"] : j1[j1["_unique_interval"] + "_position"];
				mut["_end"] = (!uc1_strand) ? j1[j1["_unique_interval"] + "_position"] : j2[j2["_unique_interval"] + "_position"];
				mut["repeat_name"] = j1["_" + j1["_is_interval"] + "_is_name"];

				//print Dumper($j1, $j2) if ($verbose);

				mut["position"] = mut["_start"]; // - 1; //position is the first duplicated base...
				mut["duplication_size"] = s(n(mut["_end"]) - n(mut["_start"]) + 1);

				// ok, we're actually missing a base of the reference...
				if (n(mut["duplication_size"]) < 0)
				{
					mut["position"] = s(n(mut["position"]) + n(mut["duplication_size"]));
				}

				// get any unique junction sequence
				JunctionInfo j1i = junction_name_split(j1["key"]);
				string j1_unique_read_sequence = j1i.unique_read_sequence;

				JunctionInfo j2i = junction_name_split(j2["key"]);
				string j2_unique_read_sequence = j2i.unique_read_sequence;

				// _gap_left and _gap_right also refer to the top strand of the genome

				mut["_ins_start"] = "";
				mut["_ins_end"] = "";

				mut["_del_start"] = "0";
				mut["_del_end"] = "0";

				// sometimes the ends of the IS are not quite flush
				string j1_not_flush_seq = "";
				if (!b(j1[j1["_is_interval"] + "_strand"]))
				{
					mut["_gap_left"] = s(n(j1[j1["_is_interval"] + "_position"]) - n(j1["_" + j1["_is_interval"] + "_is_end"]));

					if (n(mut["_gap_left"]) > 0)
					{
						j1_not_flush_seq = get_sequence (
							j1[j1["_is_interval"] + "_seq_id"],
							n(j1["_" + j1["_is_interval"] + "_is_end"]) + 1,
							n(j1[j1["_is_interval"] + "_position"])
						);
					}
				}
				else
				{
					mut["_gap_left"] = s(n(j1["_" + j1["_is_interval"] + "_is_start"]) - n(j1[j1["_is_interval"] + "_position"]));
					if (n(mut["_gap_left"]) > 0)
					{
						j1_not_flush_seq = get_sequence (
							j1[j1["_is_interval"] + "_seq_id"],
							n(j1[j1["_is_interval"] + "_position"]),
							n(j1["_" + j1["_is_interval"] + "_is_start"]) - 1
						);
					}
				}

				if (n(mut["_gap_left"]) >= 0)
				{
					if (verbose)
						cout << "J1 NF:" << j1_not_flush_seq << " U:" << j1_unique_read_sequence << endl;

					if (j1["_" + j1["_is_interval"] + "_read_side"] != j1[j1["_is_interval"] + "_strand"])
					{
						j1_not_flush_seq = reverse_complement(j1_not_flush_seq);
					}

					if (!b(j1["_" + j1["_is_interval"] + "_read_side"]))
					{
						mut["_gap_left"] = j1_not_flush_seq + j1_unique_read_sequence;
					}
					else
					{
						mut["_gap_left"] = j1_unique_read_sequence + j1_not_flush_seq;
					}
					mut["_ins_start"] = mut["_gap_left"];
				}
				else if (n(mut["_gap_left"]) < 0)
				{
					mut["_del_start"] = s(abs(n(mut["_gap_left"])));
				}

				string j2_not_flush_seq = "";
				if (!b(j2[j2["_is_interval"] + "_strand"]))
				{
					mut["_gap_right"] = s(n(j2[j2["_is_interval"] + "_position"]) - n(j2["_" + j2["_is_interval"] + "_is_end"]));
					if (n(mut["_gap_right"]) > 0)
					{
						j2_not_flush_seq = get_sequence (
							j1[j2["_is_interval"] + "_seq_id"],
							n(j2["_" + j2["_is_interval"] + "_is_end"]) + 1,
							n(j2[j2["_is_interval"] + "_position"])
						);
					}
				}
				else
				{
					mut["_gap_right"] = s(n(j2["_" + j2["_is_interval"] + "_is_start"]) - n(j2[j2["_is_interval"] + "_position"]));
					if (n(mut["_gap_right"]) > 0)
					{
						j2_not_flush_seq = get_sequence (
							j1[j2["_is_interval"] + "_seq_id"],
							n(j2[j2["_is_interval"] + "_position"]),
							n(j2["_" + j2["_is_interval"] + "_is_start"]) - 1
						);
					}
				}


				if (n(mut["_gap_right"]) >= 0)
				{
					if (verbose)
						cout << "J2 NF:" << j2_not_flush_seq << " U:" << j2_unique_read_sequence << endl;

					if (j2["_" + j2["_is_interval"] + "_read_side"] != j2[j2["_is_interval"] + "_strand"])
					{
						j2_not_flush_seq = reverse_complement(j2_not_flush_seq);
					}

					if (!b(j2["_" + j2["_is_interval"] + "_read_side"]))
					{
						mut["_gap_right"] = j2_not_flush_seq + j2_unique_read_sequence;
					}
					else
					{
						mut["_gap_right"] = j2_unique_read_sequence + j2_not_flush_seq;
					}
					mut["_ins_end"] = mut["_gap_right"];
				}
				else if (n(mut["_gap_right"]) < 0)
				{
					mut["_del_end"] = s(abs(n(mut["_gap_right"])));
				}

				// At this point any added junction sequences are on the strand as you would see them in the alignment.
				// we may need to reverse complement and change sides.


				// At this point any added junction sequences are on the strand as you would see them in the alignment.
				// we may need to reverse complement and change sides.

				if (verbose)
					cout << mut["_gap_left"] << " :: " << mut["_gap_right"] << endl;

				if (j1[j1["_unique_interval"] + "_strand"] != j1["_" + j1["_unique_interval"] + "_read_side"])
				{
					if (verbose) cout << "RC left" << endl;
					mut["_ins_start"] = reverse_complement(mut["_ins_start"]);
				}

				if (j2[j2["_unique_interval"] + "_strand"] != j2["_" + j2["_unique_interval"] + "_read_side"])
				{
					if (verbose) cout << "RC right" << endl;
					mut["_ins_end"] = reverse_complement(mut["_ins_end"]);
				}

				//// Check for ambiguous insertion direction!
				// Sometimes a strand will be assigned just because there is a 50-50 chance of getting the correct sides of the IS.
				// We need to actually check the sequence on each side of the repeat element on the end as far in as the maximum overlap on that side.
				// Use:			{max_left"] {max_right"]
				// Retrieve sequence on unique side and compare to sequence on the other side of a repeat element

				uint32_t j1_not_flush_length = j1_not_flush_seq.size();
				uint32_t j2_not_flush_length = j2_not_flush_seq.size();
				uint32_t max_not_flush_length = max(j1_not_flush_length, j2_not_flush_length);

				if (verbose) {
					cout << "J1 not flush length: " << j1_not_flush_length << endl;
					cout << "J2 not flush length: " << j2_not_flush_length << endl;
					cout << "Max not flush length: " << max_not_flush_length << endl;
				}

				int32_t j1_is_overlap_length = n(!b(j1["_" + j1["_is_interval"] + "_read_side"]) ? j1["max_left"] : j1["max_right"]);
				int32_t j2_is_overlap_length = n(!b(j2["_" + j2["_is_interval"] + "_read_side"]) ? j2["max_left"] : j2["max_right"]);

				if (verbose) {
					cout << "J1 IS overlap length: " << j1_is_overlap_length << endl;
					cout << "J2 IS overlap length: " << j2_is_overlap_length << endl;
				}

				string j1_is_seq_matched = "";
				if (!b(j1[j1["_is_interval"] + "_strand"]))
				{
					j1_is_seq_matched = get_sequence (
						j1[j1["_is_interval"] + "_seq_id"],
						n(j1[j1["_is_interval"] + "_position"]) - (j1_is_overlap_length - 1),
						n(j1[j1["_is_interval"] + "_position"]) - j1_not_flush_length
					);
					j1_is_seq_matched = reverse_complement(j1_is_seq_matched);
				}
				else
				{
					j1_is_seq_matched = get_sequence (
						j1[j1["_is_interval"] + "_seq_id"],
						n(j1[j1["_is_interval"] + "_position"]) + j1_not_flush_length,
						n(j1[j1["_is_interval"] + "_position"]) + j1_is_overlap_length - 1
					);
				}

				string j2_is_seq_matched = "";
				if (!b(j2[j2["_is_interval"] + "_strand"]))
				{
					j2_is_seq_matched = get_sequence (
						j2[j2["_is_interval"] + "_seq_id"],
						n(j2[j2["_is_interval"] + "_position"]) - (j2_is_overlap_length - 1),
						n(j2[j2["_is_interval"] + "_position"]) - j2_not_flush_length
					);
					j2_is_seq_matched = reverse_complement(j2_is_seq_matched);
				}
				else
				{
					j2_is_seq_matched = get_sequence (
						j2[j2["_is_interval"] + "_seq_id"],
						n(j2[j2["_is_interval"] + "_position"]) + j2_not_flush_length,
						n(j2[j2["_is_interval"] + "_position"]) + j2_is_overlap_length - 1
					);
				}

				// what are the actual sequences of this length at the end of the IS elements?

				string j1_left_is_sequence = get_sequence (
					j1[j1["_is_interval"] + "_seq_id"],
					n(j1["_" + j1["_is_interval"] + "_is_start"]),
					n(j1["_" + j1["_is_interval"] + "_is_start"]) + j1_is_overlap_length - 1
				);

				string j1_right_is_sequence = get_sequence (
					j1[j1["_is_interval"] + "_seq_id"],
					n(j1["_" + j1["_is_interval"] + "_is_end"]) - (j1_is_overlap_length - 1),
					n(j1["_" + j1["_is_interval"] + "_is_end"])
				);
				j1_right_is_sequence = reverse_complement(j1_right_is_sequence);

				if (verbose) {
					cout << "J1 LEFT : " << j1_left_is_sequence << endl;
					cout << "J1 RIGHT: " << j1_right_is_sequence << endl;
				}

				// believe the direction if the sequences are different
				bool j1_is_ambiguous = (j1_left_is_sequence == j1_right_is_sequence);

				string j2_left_is_sequence = get_sequence (
					j2[j2["_is_interval"] + "_seq_id"],
					n(j2["_" + j2["_is_interval"] + "_is_start"]),
					n(j2["_" + j2["_is_interval"] + "_is_start"]) + j2_is_overlap_length - 1
				);

				string j2_right_is_sequence = get_sequence (
					j2[j2["_is_interval"] + "_seq_id"],
					n(j2["_" + j2["_is_interval"] + "_is_end"]) - (j2_is_overlap_length - 1),
					n(j2["_" + j2["_is_interval"] + "_is_end"])
				);
				j2_right_is_sequence = reverse_complement(j2_right_is_sequence);

				// believe the direction if the sequences are different
				bool j2_is_ambiguous = (j2_left_is_sequence == j2_right_is_sequence);

				
				if (verbose) {
					cout << "J2 LEFT : " << j2_left_is_sequence << endl;
					cout << "J2 RIGHT: " << j2_right_is_sequence << endl;

					cout << "J1 IS matched length " << j1_is_overlap_length << ": " << j1_is_seq_matched << endl;
					cout << "J2 IS matched length " << j2_is_overlap_length << ": " << j2_is_seq_matched << endl;
				}

				// if the matched IS element sequences are the same then the direction is AMBIGUOUS
				bool is_strand = false;
				if (j1_is_ambiguous && j2_is_ambiguous)
				{
					if (verbose) cout << "AMBIGUOUS strand for mobile element insertion" << endl;
				}
				else if (j1_is_ambiguous)
				{
					is_strand = is2_strand;
				}
				else if (j2_is_ambiguous)
				{
					is_strand = is1_strand;
				}
				else // neither is ambiguous, hopefully the strands agree
				{
					assert (is1_strand == is2_strand);
          is_strand = is1_strand;
				}
				mut["strand"] = s(is_strand);

				////
				//// We are still not checking for a case where one junction side extends far enough to uniquely define the
				//// side of the IS, but the other side does not (giving it the wrong strand).


				// Finally, do this AFTER checking for the IS-matched sequences...
				// $j1 may be the left side, rather than the right side of the insertion, if so...
				if (uc1_strand)
				{
					if (verbose) cout << "reverse right and left" << endl;
					swap(mut["_ins_start"], mut["_ins_end"]);
					swap(mut["_del_start"], mut["_del_end"]);
				}

				// clean up unused keys
				if (mut.entry_exists("_del_start")) mut["del_start"] = mut["_del_start"];
				if (mut.entry_exists("_del_end")) mut["del_end"] = mut["_del_end"];

				 if (mut.entry_exists("_ins_start")) mut["ins_start"] = mut["_ins_start"];
				 if (mut.entry_exists("_ins_end")) mut["ins_end"] = mut["_ins_end"];

				if (verbose)
					cout << mut["_gap_left"] << " :: " << mut["_gap_right"] << endl;

				gd.add(mut);
				break; // next JC1
			}
		}

		///
		// evidence JC => INS, SUB, AMP mutations
		///

		for (uint32_t i = 0; i < jc.size(); i++) //JC
		{
			diff_entry& j = *jc[i];

      //cout << j["side_1_seq_id"] << " " << j["side_2_seq_id"] << endl;
      //cout << j["side_1_strand"] << " " << j["side_2_strand"] << endl;
      
			if (
				// must be on same sequence
				   (j["side_1_seq_id"] != j["side_2_seq_id"])
				// must be in same orientation (implies strands are opposite)
				|| (j["side_1_strand"] == j["side_2_strand"])
			)
				continue;

			string seq_id = j["side_1_seq_id"];

			// We can assume that the lower coordinate will be first since this is NOT a deletion
			// (which would be handled above)
      // By this point any positive overlap should have been resolved.
			assert(n(j["overlap"]) <= 0);
      
			// mutation will always be after this position
			int32_t position = n(j["side_1_position"]);

			// Special case of circular chromosome
      // TODO: change to use ref_seq_info length field
			if ( (j["side_1_position"] == "1") && ( static_cast<uint32_t>(n(j["side_2_position"])) == ref_seq_info[ref_seq_info.seq_id_to_index(j["side_2_seq_id"])].m_length ) )
			{
				j["circular_chromosome"] = "1";
				continue;
			}
			// protection against mistakes
			if (n(j["side_2_position"]) - n(j["side_1_position"]) + 1 > 100000)
				continue;

			// 'AMP'
			if (!b(j["unique_read_sequence"]))
			{
				int32_t size = n(j["side_2_position"]) - n(j["side_1_position"]) + 1;
				if (size < 0) continue; // this is a deletion!
				//#	next if ($size > 100); #spurious duplication, need extra evidence from coverage!

				diff_entry mut;
				mut._type = "AMP";
				mut
					("seq_id", seq_id)
					("position", s(position))
					("size", s(size))
					("new_copy_number", "2")
				;
				mut._evidence = make_list<string>(j._id);
				gd.add(mut);
			}
			// 'SUB'
			else if (n(j["side_1_position"]) >= n(j["side_2_position"]))
			{
				string ref_seq = "";
				string new_seq = j["unique_read_sequence"];
				if (n(j["side_1_position"]) >= n(j["side_2_position"])) //TODO: When would this be false?
				{
					new_seq = get_sequence (
						seq_id,
						n(j["side_2_position"]),
						n(j["side_1_position"])
					);
				}

				diff_entry mut;
				mut._type = "SUB";
				mut
					("seq_id", seq_id)
					("position", s(position))
					("size", s(n(j["side_1_position"]) - n(j["side_2_position"]) + 1))
					("new_seq", new_seq)
				;
				mut._evidence = make_list<string>(j._id);

				gd.add(mut);
			}
			// "INS"
			else if (n(j["side_1_position"]) + 1 == n(j["side_2_position"]))
			{
				diff_entry mut;
				mut._type = "INS";
				mut
					("seq_id", seq_id)
					("position", s(position))
					("new_seq", j["unique_read_sequence"])
				;
				mut._evidence = make_list<string>(j._id);
				gd.add(mut);
			}
		}


		///
		// Read Alignments => SNP, DEL, INS, SUB
		///

		vector<string> ra_types = make_list<string>("RA");
    genome_diff::entry_list_t ra = gd.list(ra_types);

		///
		// Ignore RA that overlap DEL or MC
		// They are due to low spurious coverage in deleted regions!
		///

		{
			vector<string> del_types = make_list<string>("DEL");
			vector<string> mc_types = make_list<string>("MC");
			genome_diff::entry_list_t del = gd.list(del_types);
			genome_diff::entry_list_t mc = gd.list(mc_types);

			for (uint32_t i = 0; i < ra.size(); i++) // RA
			{
				diff_entry& ra_item = *ra[i];

				bool next_ra = false;
				for (uint32_t j = 0; j < del.size(); j++) // DEL
				{
					diff_entry& del_item = *del[j];
					if (ra_item["seq_id"] != del_item["seq_id"])
						continue;

					// there might be a problem here with insert_position > 0
					if ( (n(ra_item["position"]) >= n(del_item["position"])) && (n(ra_item["position"]) <= n(del_item["position"]) + n(del_item["size"]) - 1) )
					{
						ra_item["deleted"] = 1;
						next_ra = true;
						break;
					}
				}
				if (next_ra) continue;

				for (uint32_t j = 0; j < mc.size(); j++) // MC
				{
					diff_entry& mc_item = *mc[j];

					if (ra_item["seq_id"] != mc_item["seq_id"]) continue;

					if ( (n(ra_item["position"]) >= n(mc_item["start"])) && (n(ra_item["position"]) <= n(mc_item["end"])) )
					{
						ra_item["deleted"] = 1;
						break;
					}
				}

			}
		}

		sort(ra.begin(), ra.end(), MutationPredictor::sort_by_pos);

		///
		// Gather together read alignment mutations that occur next to each other
		// ...unless they are polymorphisms
		///

		bool first_time = true;
		diff_entry mut;
		vector<diff_entry> muts;

		for (uint32_t i = 0; i < ra.size(); i++) // RA
		{
			diff_entry& item = *ra[i];

			if ( item.entry_exists("reject") || item.entry_exists("deleted"))
			  continue;

			// Sometimes a SNP might be called in a deleted area because the end was wrong,
			// but it was corrected using a junction. (This catches this case.)

			bool same = false;
			if (!first_time)
			{
				if ( ((mut["end"] == item["position"]) && (n(mut["insert_end"]) + 1 == n(item["insert_position"])))
						|| ((n(mut["end"]) + 1 == n(item["position"])) && (item["insert_position"] == "0")) )
					same = true;
				if ( (item["frequency"] != "1") || (mut["frequency"] != "1") //don't join polymorphisms
						|| (mut["seq_id"] != item["seq_id"]) )
					same = false;
			}

			if (!same)
			{
				if (!first_time) muts.push_back(mut);
				first_time = false;
				diff_entry new_mut;
				new_mut._evidence = make_list<string>(item._id);
				new_mut
					("seq_id", item["seq_id"])
					("position", item["position"])
					("start", item["position"])
					("end", item["position"])
					("insert_start", item["insert_position"])
					("insert_end", item["insert_position"])
					("ref_seq", (item["ref_base"] != ".") ? item["ref_base"] : "")
					("new_seq", (item["new_base"] != ".") ? item["new_base"] : "")
					("frequency", item["frequency"])
				;
				mut = new_mut;
			}
			else
			{
				mut
					("insert_end", item["insert_position"])
					("end", item["position"])
				;
				if (item["ref_base"] != ".") mut["ref_seq"] += item["ref_base"];
				if (item["new_base"] != ".") mut["new_seq"] += item["new_base"];
				mut._evidence.push_back(item._id);
			}
		}
		//don't forget the last one
		if (!first_time) muts.push_back(mut);

		///
		// Finally, convert these items into the fields needed for the various types of mutations
		///

		for (uint32_t i = 0; i < muts.size(); i++)
		{
			mut = muts[i];
			// insertion
			if (mut["ref_seq"].size() == 0)
			{
				mut._type = "INS";
				// unused fields
				mut._fields.erase("ref_seq");
			}
			// deletion
			else if (mut["new_seq"].size() == 0)
			{
				mut._type = "DEL";
				mut["size"] = s(n(mut["end"]) - n(mut["start"]) + 1);

				// unused fields
				mut._fields.erase("new_seq");
				mut._fields.erase("ref_seq");
			}
			// block substitution
			else if ((mut["ref_seq"].size() > 1) || (mut["new_seq"].size() > 1))
			{
				mut._type = "SUB";
				mut["size"] = s(mut["ref_seq"].size());
				mut._fields.erase("ref_seq");
			}
			//snp
			else
			{
				mut._fields.erase("ref_seq");
				mut._type = "SNP";
			}

			// we don't need these fields
			if (mut["frequency"] == "1") mut._fields.erase("frequency");
			mut._fields.erase("start");
			mut._fields.erase("end");
			mut._fields.erase("insert_start");
			mut._fields.erase("insert_end");

			gd.add(mut);
		}

    /* @JEB This appears to have been made redundant with new junction methods!
		// Remove remaining junctions that we didn't pair up with anything that are below a coverage cutoff.
		jc = gd.filter_used_as_evidence(gd.list(jc_types));
		for (uint32_t i = 0; i < jc.size(); i++)
		{
			diff_entry& item = *jc[i];

			int32_t coverage_cutoff_1 = settings.unique_coverage[item["side_1_seq_id"]].junction_coverage_cutoff;
			int32_t coverage_cutoff_2 = settings.unique_coverage[item["side_2_seq_id"]].junction_coverage_cutoff;

			//TODO: What range of values means undefined? Is 0 undefined?
			if ( (coverage_cutoff_1 == 0 || (n(item["total_reads"]) < coverage_cutoff_1) )
				&& (coverage_cutoff_2 == 0 || (n(item["total_reads"]) < coverage_cutoff_2) ) )
			{
				add_reject_reason(item, "COV");
			}
		}
    */
		
	}


} // namespace breseq
