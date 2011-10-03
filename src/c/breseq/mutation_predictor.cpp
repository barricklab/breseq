
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

#include "libbreseq/mutation_predictor.h"

using namespace std;

namespace breseq {

  cReferenceSequences MutationPredictor::ref_seq_info;
  
	MutationPredictor::MutationPredictor(cReferenceSequences& _ref_seq_info)
	{
		 ref_seq_info = _ref_seq_info;
	}

	// Private methods

	cSequenceFeature* MutationPredictor::within_repeat(string seq_id, uint32_t position)
	{
		cSequenceFeatureList& rl = ref_seq_info[seq_id].m_repeats;
    cSequenceFeature* r= NULL;
    
    // by returning the last one we encounter that we are inside, 
    // we get the inner repeat in nested cases
		for (uint32_t i = 0; i < rl.size(); i++)
			if ((rl[i]->m_start <= position) && (position <= rl[i]->m_end))
				r = rl[i].get();

		return r;
	}

	bool MutationPredictor::sort_by_hybrid(const counted_ptr<diff_entry>& a, const counted_ptr<diff_entry>& b)
	{
    int32_t a_pos = n(a->entry_exists("_side_1_is") ? (*a)["side_2_position"] : (*a)["side_1_position"]);
		int32_t b_pos = n(b->entry_exists("_side_1_is") ? (*b)["side_2_position"] : (*b)["side_1_position"]);

		int32_t a_seq_order = (a->entry_exists("_side_1_is") ? ref_seq_info.seq_order[(*a)["side_2_seq_id"]] : ref_seq_info.seq_order[(*a)["side_1_seq_id"]]);
		int32_t b_seq_order = (b->entry_exists("_side_1_is") ? ref_seq_info.seq_order[(*b)["side_2_seq_id"]] : ref_seq_info.seq_order[(*b)["side_1_seq_id"]]);

		uint32_t a_reject_order = a->number_reject_reasons();
		uint32_t b_reject_order = b->number_reject_reasons();

		// sort by seq_id, position, fewer reject reasons, then score (highest to lowest)
    
    if (a_seq_order != b_seq_order) 
      return (a_seq_order < b_seq_order);
    if (a_pos != b_pos) 
      return (a_pos < b_pos);
    if (a_reject_order != b_reject_order) 
      return (a_reject_order < b_reject_order);

    return (n((*a)["pos_hash_score"]) > n((*b)["pos_hash_score"]));
  }

	bool MutationPredictor::sort_by_reject_score(const counted_ptr<diff_entry>& a, const counted_ptr<diff_entry>& b)
	{
		uint32_t a_reject_order = a->number_reject_reasons();
		uint32_t b_reject_order = b->number_reject_reasons();

		// sort by seq_id, position, fewer reject reasons, then score (highest to lowest)
    if (a_reject_order != b_reject_order) 
      return (a_reject_order < b_reject_order);
            
    return (n((*a)["pos_hash_score"]) > n((*b)["pos_hash_score"]));
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
	void MutationPredictor::predict(Settings& settings, genome_diff& gd, uint32_t max_read_length, double avg_read_length)
	{
		(void)settings; //TODO; unused?
    bool verbose = false; // for debugging

    //@JEB This could be replaced by passing summary
    if (avg_read_length == 0.0) avg_read_length = max_read_length;

		///
		//  Preprocessing of JC evidence
		///

		// For all that follows, we need information about repeat_regions overlapping the sides of junctions
    vector<Type> jc_types = make_list<Type>(JC);
		diff_entry_list jc = gd.list(jc_types);
    
    const int32_t max_distance_to_repeat = 50;

		for (diff_entry_list::iterator jc_it=jc.begin(); jc_it!=jc.end(); jc_it++)
		{
			diff_entry& j = **jc_it;
      
			j["_side_1_read_side"] = "-1";
			j["_side_2_read_side"] = "1";

			for (uint32_t side = 1; side <= 2; side++)
			{
				string side_key = "side_" + s(side);

				cSequenceFeature* is = ref_seq_info.find_closest_repeat_region(
					n(j[side_key + "_position"]),
					ref_seq_info[j[side_key + "_seq_id"]].m_repeats,
					max_distance_to_repeat,
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
				
        uint32_t test_val = n(j[side_key + "_redundant"]);
				j[side_key + "_annotate_key"] = (j.entry_exists("_" + side_key + "_is_start") || n(j[side_key + "_redundant"])) ? "repeat" : "gene";
			}

			// by default, we are sorted by this coord
			j["_unique_interval"] = "side_1";

			// Determine which side of the junction is the IS and which is unique
			// these point to the correct initial interval...
			if (j.entry_exists("_side_1_is"))
			{
        //cout << n(j["_side_1_is_start"]) << " " << n(j["_side_1_is_end"]) << " " << n(j["side_1_position"]) << endl;
        
				if (abs(n(j["_side_1_is_start"]) - n(j["side_1_position"])) <= max_distance_to_repeat)
				{
					j["_is_interval"] = "side_1";
					j["_is_interval_closest_side_key"] = "start";
					j["_unique_interval"] = "side_2";
				}
				else if (abs(n(j["_side_1_is_end"]) - n(j["side_1_position"])) <= max_distance_to_repeat)
				{
					j["_is_interval"] = "side_1";
					j["_is_interval_closest_side_key"] = "end";
					j["_unique_interval"] = "side_2";
				}
			}
      
      if (j.entry_exists("_side_2_is"))
			{
        //cout << n(j["_side_2_is_start"]) << " " << n(j["_side_2_is_end"]) << " " << n(j["side_2_position"]) << endl;

				if (abs(n(j["_side_2_is_start"]) - n(j["side_2_position"])) <= max_distance_to_repeat)
				{
					j["_is_interval"] = "side_2";
					j["_is_interval_closest_side_key"] = "start";
					j["_unique_interval"] = "side_1";
				}
				else if (abs(n(j["_side_2_is_end"]) - n(j["side_2_position"])) <= max_distance_to_repeat)
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
		/*for(diff_entry_list::iterator jc_it = jc.begin(); jc_it != jc.end(); jc_it++)
    {
      diff_entry& de = **jc_it;
			if (de.entry_exists("reject"))
      {
				jc.erase(jc_it);
        jc_it--;
      }
		}*/
    jc.remove_if(diff_entry::field_exists("reject"));
    
    vector<Type> mc_types = make_list<Type>(MC);
		diff_entry_list mc = gd.list(mc_types);

		///
		// evidence MC + JC => DEL mutation
		///

		// DEL prediction:
		// (1) there is a junction that exactly crosses the deletion boundary deletion
		// (2) there is is no junction, but both ends of the deletion are in repeat sequences
		// (3) there is a junction between unique sequence and a repeat element

    for(diff_entry_list::iterator mc_it = mc.begin(); mc_it != mc.end(); mc_it++)
    {
      diff_entry& mc_item = **mc_it;
      
			if (mc_item.entry_exists("reject"))
			  continue;

			// set up generic deletion item
			diff_entry mut;
      mut._type = DEL;
			mut._evidence = make_list<string>(mc_item._id);
			mut
				("seq_id", mc_item["seq_id"])
				("position", mc_item["start"])
				("size", s(n(mc_item["end"]) - n(mc_item["start"]) + 1));
			;

			///
			// (0) this is a deletion of an entire fragment
			///     
      
      uint32_t tid = ref_seq_info.seq_id_to_index(mut[SEQ_ID]); 
      if ( (n(mut[POSITION]) == 1) && (n(mut[POSITION]) + n(mut["size"]) - 1 == static_cast<int32_t>(ref_seq_info[tid].m_length)) )
      {
        gd.add(mut);
        continue;
      }
      
      
			///
			// (1) there is a junction that exactly crosses the deletion boundary 
			///

			for(diff_entry_list::iterator jc_it = jc.begin(); jc_it != jc.end(); jc_it++) //JC
			{
				diff_entry& jc_item = **jc_it;

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
          
          //it's possible that one or both sides are in repeat elements
          cSequenceFeature* r1_pointer = within_repeat(jc_item["side_1_seq_id"], n(jc_item["side_1_position"]));
          cSequenceFeature* r2_pointer = within_repeat(jc_item["side_2_seq_id"], n(jc_item["side_2_position"]));
          
          // one repeat cases where the end matches up exactly
          if (r1_pointer) 
          {
            // must match up to an end of the repeat
            if ((n(jc_item["side_1_position"]) == static_cast<int32_t>(r1_pointer->m_start))
             || (n(jc_item["side_1_position"]) == static_cast<int32_t>(r1_pointer->m_end)))
            {
              mut["mediated"] = (*r1_pointer)["name"];
            }
          }
          else if (r2_pointer) 
          {
            // must match up to an end of the repeat
            if ((n(jc_item["side_2_position"]) == static_cast<int32_t>(r2_pointer->m_start))
                || (n(jc_item["side_2_position"]) == static_cast<int32_t>(r2_pointer->m_end)))
            {
              mut["mediated"] = (*r2_pointer)["name"];
            }
          }    
                      
					mut._evidence.push_back(jc_item._id);
					jc.erase(jc_it);
          jc_it--;
					gd.add(mut);
					continue;
				}
			}

      cSequenceFeature* r1_pointer = within_repeat(mut["seq_id"], n(mut["position"]));
      cSequenceFeature* r2_pointer = within_repeat(mut["seq_id"], n(mut["position"]) + n(mut["size"]));
      
			///
			// (2) there is is no junction, but both ends of the deletion are in different copies of the same repeat sequence
			///

			// Then we will adjust the coordinates to remove...
			if  (
           (r1_pointer != r2_pointer) 
           && (r1_pointer != NULL) 
           && (r2_pointer != NULL) 
           && ((*r1_pointer)["name"] == (*r2_pointer)["name"])
          )
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
			for(diff_entry_list::iterator jc_it = jc.begin(); jc_it != jc.end(); jc_it++) //JUNCTION
			{
				diff_entry& j = **jc_it;

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


				// check that IS is on the right strand
				if (verbose)
					cout << "Check 4: " << redundant_deletion_side << " * " << to_string(r.m_strand) << " != " << j[j["_is_interval"] + "_strand"] << " * " << j["_" + j["_is_interval"] + "_is_strand"] << endl;
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
				if (redundant_deletion_side == -1)
				{
					uint32_t move_dist = r.m_end + 1 - n(mut["position"]);
					mut["position"] = s(n(mut["position"]) + move_dist);
					mut["size"] = s(n(mut["size"]) - move_dist);
				}
				else
				{
					int32_t move_dist = (n(mut["position"]) + n(mut["size"]) - 1) - (r.m_start - 1);
					mut["size"] = s(n(mut["size"]) - move_dist);
				}

				// OK, we're good!
				mut["mediated"] = r["name"];
				mut._evidence.push_back(j._id);
				jc.erase(jc_it);
        jc_it--;
				gd.add(mut);
				ok_were_good = true;
			}
			if (ok_were_good) continue;
		}


		///
		// evidence JC + JC = MOB mutation
		///

		for(diff_entry_list::iterator it = jc.begin(); it != jc.end(); it++) //JC
		{
			diff_entry& j = **it;
		
			// Junction isn't near an IS. Move on.
			if (!j.entry_exists("_is_interval")) continue;

			// There is no overlap to correct. Move on.
			if (n(j["overlap"]) <= 0) continue;

			// The following code implies n(j["overlap"]) > 0

			/// first, adjust the repetitive sequence boundary to get as close to the IS as possible
			int32_t is_interval_position = n(j[j["_is_interval"] + "_position"]);
			uint32_t move_dist = abs(is_interval_position - n(j[j["_is_interval"] + "_is_" + j["_is_interval_closest_side_key"]])); //  + "_position" ?
			uint32_t overlap = n(j["overlap"]);
			if (move_dist > overlap) move_dist = overlap;

			is_interval_position += n(j[j["_is_interval"] + "_strand"]) * move_dist;
			j[j["_is_interval"] + "_position"] = s(is_interval_position);
			overlap -= move_dist;

			/// second, adjust the unique sequence side with any remaining overlap
			int32_t unique_interval_position = n(j[j["_unique_interval"] + "_position"]);
			unique_interval_position += n(j[j["_unique_interval"] + "_strand"]) * overlap;
			j[j["_unique_interval"] + "_position"] = s(unique_interval_position);

			j["overlap"] = "0";
		}

		jc.sort(MutationPredictor::sort_by_hybrid);

    for(diff_entry_list::iterator jc1_it = jc.begin(); jc1_it != jc.end(); jc1_it++) //JC1
		{
			diff_entry& j1 = **jc1_it;
      
			// Compile a list of the next possibilities within a certain length of bases
      vector<diff_entry_ptr> j2_list;
			vector<diff_entry_list::iterator> it_delete_list_2;
      
      // start looking at the next JC entry
      diff_entry_list::iterator jc2_it = jc1_it;
      for(jc2_it++ ; jc2_it != jc.end(); jc2_it++) //JC2
			{
				diff_entry& j2 = **jc2_it;
        
				// must be close together in real coords
				if ( (j1[j1["_unique_interval"] + "_seq_id"] != j2[j2["_unique_interval"] + "_seq_id"])
					|| (abs(n(j1[j1["_unique_interval"] + "_position"]) - n(j2[j2["_unique_interval"] + "_position"])) > 20 ) )
					break;

				if ( (!j1.entry_exists("_is_interval") || !j2.entry_exists("_is_interval"))
					|| (j1["_" + j1["_is_interval"] + "_is_name"] != j2["_" + j2["_is_interval"] + "_is_name"]) )
					continue;
        
        it_delete_list_2.push_back(jc2_it);        
				j2_list.push_back(*jc2_it);
			}
      if (verbose)
        cout << "Size of J2 list: " << j2_list.size() << endl;

			//sort the $j2_list by reject reason and score

			sort(j2_list.begin(), j2_list.end(), sort_by_reject_score);

			// We need to go through all with the same coordinate (or within a certain coordinate stretch?)
			// because sometimes a failed junction will be in between the successful junctions
      for(size_t i=0; i<j2_list.size(); i++)
			{
				diff_entry& j2 = *(j2_list[i]);

        if (verbose) 
        {
          cout << "Sorted: == J1 ==" << endl;
          for(map<string,string>::iterator it=j1._fields.begin(); it!=j1._fields.end(); it++)
          {
            cout << it->first << " = " << it->second << endl; 
          }
          cout << "Sorted: == J2 ==" << endl;
          for(map<string,string>::iterator it=j2._fields.begin(); it!=j2._fields.end(); it++)
          {
            cout << it->first << " = " << it->second << endl; 
          }
        }
        
				// positive overlap should be resolved by now
				assert(n(j1["overlap"]) <= 0);
				assert(n(j2["overlap"]) <= 0);

				// the first unique coords are going into the IS element
				int32_t uc1_strand = n(j1[j1["_unique_interval"] + "_strand"]);
				int32_t uc2_strand = n(j2[j2["_unique_interval"] + "_strand"]);
				if (uc1_strand != -uc2_strand) continue;

				// What strand is the IS on relative to the top strand of the genome
				int32_t is1_strand = - (n(j1[j1["_is_interval"] + "_strand"]) * n(j1["_" + j1["_is_interval"] + "_is_strand"]) * n(j1[j1["_unique_interval"] + "_strand"]));
				int32_t is2_strand = - (n(j2[j2["_is_interval"] + "_strand"]) * n(j2["_" + j2["_is_interval"] + "_is_strand"]) * n(j2[j2["_unique_interval"] + "_strand"]));

				// Remove these predictions from the list
				jc.erase(jc1_it);
        jc1_it--;
				jc.erase(it_delete_list_2[i]);

				// Create the mutation, with evidence
				diff_entry mut;
        mut._type = MOB;
				mut._evidence.push_back(j1._id);
				mut._evidence.push_back(j2._id);
				mut
					("seq_id", j1[j1["_unique_interval"] + "_seq_id"])
				;
				mut["_start"] = (uc1_strand == -1) ? j2[j2["_unique_interval"] + "_position"] : j1[j1["_unique_interval"] + "_position"];
				mut["_end"] = (uc1_strand == -1) ? j1[j1["_unique_interval"] + "_position"] : j2[j2["_unique_interval"] + "_position"];
				mut["repeat_name"] = j1["_" + j1["_is_interval"] + "_is_name"];

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
				if (n(j1[j1["_is_interval"] + "_strand"]) == -1)
				{
					mut["_gap_left"] = s(n(j1[j1["_is_interval"] + "_position"]) - n(j1["_" + j1["_is_interval"] + "_is_end"]));

					if (n(mut["_gap_left"]) > 0)
					{
						j1_not_flush_seq = ref_seq_info.get_sequence_1 (
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
						j1_not_flush_seq = ref_seq_info.get_sequence_1 (
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

					if (n(j1["_" + j1["_is_interval"] + "_read_side"]) != n(j1[j1["_is_interval"] + "_strand"]))
					{
						j1_not_flush_seq = reverse_complement(j1_not_flush_seq);
					}

					if (n(j1["_" + j1["_is_interval"] + "_read_side"]) == -1)
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
				if (n(j2[j2["_is_interval"] + "_strand"]) == -1)
				{
					mut["_gap_right"] = s(n(j2[j2["_is_interval"] + "_position"]) - n(j2["_" + j2["_is_interval"] + "_is_end"]));
					if (n(mut["_gap_right"]) > 0)
					{
						j2_not_flush_seq = ref_seq_info.get_sequence_1 (
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
						j2_not_flush_seq = ref_seq_info.get_sequence_1 (
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

					if ( n(j2["_" + j2["_is_interval"] + "_read_side"]) * n(j2[j2["_is_interval"] + "_strand"]) == -1)
					{
						j2_not_flush_seq = reverse_complement(j2_not_flush_seq);
					}

					if (n(j2["_" + j2["_is_interval"] + "_read_side"]) == -1)
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

				if ( n(j1[j1["_unique_interval"] + "_strand"]) != n(j1["_" + j1["_unique_interval"] + "_read_side"]))
				{
					if (verbose) cout << "RC left" << endl;
					mut["_ins_start"] = reverse_complement(mut["_ins_start"]);
				}

				if ( n(j2[j2["_unique_interval"] + "_strand"]) != n(j2["_" + j2["_unique_interval"] + "_read_side"]))
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

				int32_t j1_is_overlap_length = (n(j1["_" + j1["_is_interval"] + "_read_side"]) == -1) ? n(j1["max_left"]) : n(j1["max_right"]);
				int32_t j2_is_overlap_length = (n(j2["_" + j2["_is_interval"] + "_read_side"]) == -1) ? n(j2["max_left"]) : n(j2["max_right"]);

				if (verbose) {
					cout << "J1 IS overlap length: " << j1_is_overlap_length << endl;
					cout << "J2 IS overlap length: " << j2_is_overlap_length << endl;
				}

				string j1_is_seq_matched = "";
				if (n(j1[j1["_is_interval"] + "_strand"]) == -1)
				{
					j1_is_seq_matched = ref_seq_info.get_sequence_1 (
						j1[j1["_is_interval"] + "_seq_id"],
						n(j1[j1["_is_interval"] + "_position"]) - (j1_is_overlap_length - 1),
						n(j1[j1["_is_interval"] + "_position"]) - j1_not_flush_length
					);
					j1_is_seq_matched = reverse_complement(j1_is_seq_matched);
				}
				else
				{
					j1_is_seq_matched = ref_seq_info.get_sequence_1 (
						j1[j1["_is_interval"] + "_seq_id"],
						n(j1[j1["_is_interval"] + "_position"]) + j1_not_flush_length,
						n(j1[j1["_is_interval"] + "_position"]) + j1_is_overlap_length - 1
					);
				}

				string j2_is_seq_matched = "";
				if (n(j2[j2["_is_interval"] + "_strand"]) == -1)
				{
					j2_is_seq_matched = ref_seq_info.get_sequence_1 (
						j2[j2["_is_interval"] + "_seq_id"],
						n(j2[j2["_is_interval"] + "_position"]) - (j2_is_overlap_length - 1),
						n(j2[j2["_is_interval"] + "_position"]) - j2_not_flush_length
					);
					j2_is_seq_matched = reverse_complement(j2_is_seq_matched);
				}
				else
				{
					j2_is_seq_matched = ref_seq_info.get_sequence_1 (
						j2[j2["_is_interval"] + "_seq_id"],
						n(j2[j2["_is_interval"] + "_position"]) + j2_not_flush_length,
						n(j2[j2["_is_interval"] + "_position"]) + j2_is_overlap_length - 1
					);
				}

				// what are the actual sequences of this length at the end of the IS elements?

				string j1_left_is_sequence = ref_seq_info.get_sequence_1 (
					j1[j1["_is_interval"] + "_seq_id"],
					n(j1["_" + j1["_is_interval"] + "_is_start"]),
					n(j1["_" + j1["_is_interval"] + "_is_start"]) + j1_is_overlap_length - 1
				);

				string j1_right_is_sequence = ref_seq_info.get_sequence_1 (
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

				string j2_left_is_sequence = ref_seq_info.get_sequence_1 (
					j2[j2["_is_interval"] + "_seq_id"],
					n(j2["_" + j2["_is_interval"] + "_is_start"]),
					n(j2["_" + j2["_is_interval"] + "_is_start"]) + j2_is_overlap_length - 1
				);

				string j2_right_is_sequence = ref_seq_info.get_sequence_1 (
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
				int32_t is_strand = 0;
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
				else // neither is ambiguous and the strands don't agree (this is not a simple IS insertion)
				{
					if (is1_strand != is2_strand) continue;
          is_strand = is1_strand;
				}
				mut["strand"] = s(is_strand);

				////
				//// We are still not checking for a case where one junction side extends far enough to uniquely define the
				//// side of the IS, but the other side does not (giving it the wrong strand).
				////

				// Finally, do this AFTER checking for the IS-matched sequences...
				// $j1 may be the left side, rather than the right side of the insertion, if so...
				if (uc1_strand)
				{
					if (verbose) cout << "reverse right and left" << endl;
					swap(mut["_ins_start"], mut["_ins_end"]);
					swap(mut["_del_start"], mut["_del_end"]);
				}

				// only transfer the hidden _keys to normal keys that will be printed if they are different from 0
				if (mut.entry_exists("_del_start") && (mut["_del_start"] != "0")) mut["del_start"] = mut["_del_start"];
				if (mut.entry_exists("_del_end")   && (mut["_del_end"] != "0"))   mut["del_end"] = mut["_del_end"];

				 if (mut.entry_exists("_ins_start") && (mut["_del_start"].length() != 0)) mut["ins_start"] = mut["_ins_start"];
				 if (mut.entry_exists("_ins_end")   && (mut["_del_end"].length() != 0))   mut["ins_end"] = mut["_ins_end"];

				if (verbose)
					cout << mut["_gap_left"] << " :: " << mut["_gap_right"] << endl;

        // print out everything
        if (verbose)
        {
          cout << "== J1 ==" << endl;
          for(map<string,string>::iterator it=j1._fields.begin(); it!=j1._fields.end(); it++)
          {
            cout << it->first << " = " << it->second << endl; 
          }
          
          cout << "== J2 ==" << endl;
          for(map<string,string>::iterator it=j2._fields.begin(); it!=j2._fields.end(); it++)
          {
            cout << it->first << " = " << it->second << endl; 
          }
          
          cout << "== Mut ==" << endl;
          for(map<string,string>::iterator it=mut._fields.begin(); it!=mut._fields.end(); it++)
          {
            cout << it->first << " = " << it->second << endl; 
          }
        }
        
				gd.add(mut);
				break; // next JC1
			}
		}

		///
		// evidence JC => INS, SUB, AMP mutations
		///
    for(diff_entry_list::iterator jc_it = jc.begin(); jc_it != jc.end(); jc_it++) //JC
		{
			diff_entry& j = **jc_it;

      //cout << j["side_1_seq_id"] << " " << j["side_2_seq_id"] << endl;
      //cout << j["side_1_strand"] << " " << j["side_2_strand"] << endl;
      
			if (
				// must be on same sequence
				   (j["side_1_seq_id"] != j["side_2_seq_id"])
				// must be in same orientation (implies strands are opposite)
				|| (n(j["side_1_strand"]) == n(j["side_2_strand"]))
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
			if ( (j["side_1_position"] == "1") && ( static_cast<uint32_t>(n(j["side_2_position"])) == ref_seq_info[ref_seq_info.seq_id_to_index(j["side_2_seq_id"])].m_length ) )
			{
				j["circular_chromosome"] = "1";
				continue;
			}
			// protection against mistakes
			if (n(j["side_2_position"]) - n(j["side_1_position"]) + 1 > 100000)
				continue;

			// 'DEL'
			if (!j.entry_exists("unique_read_sequence"))
			{
        int32_t side_1_position = n(j["side_1_position"]);
        int32_t side_2_position = n(j["side_2_position"]);
        
        if (n(j["side_1_strand"]) == -1)
        {
          int32_t t = side_1_position;
          side_1_position = side_2_position;
          side_2_position = t;
        }

				int32_t size = side_2_position - side_1_position + 1;
        
        // @JEB TODO: Large duplications are likely spurious... need to x-ref with coverage differences.
        if (abs(size) > avg_read_length / 3.0) continue;
        
				if (size < 0) // this is a deletion!
        {
          diff_entry mut;
          mut._type = DEL;
          mut
          ("seq_id", seq_id)
          ("position", s(position+1))
          ("size", s(-size))         // note adjustment due to +1 above
          ;
          mut._evidence = make_list<string>(j._id);
          gd.add(mut);
        } 
      }
			// 'SUB'
			else if (n(j["side_1_position"]) >= n(j["side_2_position"]))
			{
				string ref_seq = "";
				string new_seq = j["unique_read_sequence"];
				if (n(j["side_1_position"]) >= n(j["side_2_position"])) //TODO: When would this be false?
				{
					new_seq = ref_seq_info.get_sequence_1 (
						seq_id,
						n(j["side_2_position"]),
						n(j["side_1_position"])
					);
				}

				diff_entry mut;
        mut._type = SUB;
				mut
					("seq_id", seq_id)
					("position", s(position))
					("size", s(n(j["side_1_position"]) - n(j["side_2_position"]) + 1))
					("new_seq", new_seq)
				;
				mut._evidence = make_list<string>(j._id);

				gd.add(mut);
			}
			// "INS" || "AMP"
			else if (n(j["side_1_position"]) + 1 == n(j["side_2_position"]))
			{
        // Check to see if unique sequence matches sequence directly before
        size_t size = j["unique_read_sequence"].size();
        size_t position = n(j["side_2_position"]) - size;
        string dup_check_seq = ref_seq_info.get_sequence_1(j["side_1_seq_id"], position, position + size - 1);
        
        if (j["unique_read_sequence"] == dup_check_seq)
        {
          diff_entry mut;
          mut._type = AMP;
          mut
          ("seq_id", seq_id)
          ("position", s(n(j["side_2_position"]) - j["unique_read_sequence"].size()))
          ("size", s(j["unique_read_sequence"].size()))
          ("new_copy_number", "2")
          ;
          mut._evidence = make_list<string>(j._id);
          gd.add(mut);
        }
        else
        {
          diff_entry mut;
          mut._type = INS;
          mut
            ("seq_id", seq_id)
            ("position", s(position))
            ("new_seq", j["unique_read_sequence"])
          ;
          mut._evidence = make_list<string>(j._id);
          gd.add(mut);
        }
			}
		}


		///
		// Read Alignments => SNP, DEL, INS, SUB
		///

    vector<Type> ra_types = make_list<Type>(RA);
    diff_entry_list ra = gd.list(ra_types);

		///
		// Ignore RA that overlap DEL or MC
		// They are due to low spurious coverage in deleted regions!
		///

		{
      vector<Type> del_types = make_list<Type>(DEL);
      vector<Type> mc_types = make_list<Type>(MC);
			diff_entry_list del = gd.list(del_types);
      diff_entry_list mc = gd.list(mc_types);

      for(diff_entry_list::iterator ra_it = ra.begin(); ra_it != ra.end(); ra_it++) //RA
      {
        diff_entry& ra_item = **ra_it;

				bool next_ra = false;
        
        for(diff_entry_list::iterator del_it = del.begin(); del_it != del.end(); del_it++) //DEL
        {
          diff_entry& del_item = **del_it;

					if (ra_item["seq_id"] != del_item["seq_id"])
						continue;

					// there might be a problem here with insert_position > 0
					if ( (n(ra_item["position"]) >= n(del_item["position"])) && (n(ra_item["position"]) <= n(del_item["position"]) + n(del_item["size"]) - 1) )
					{
						ra_item["deleted"] = "1";
						next_ra = true;
						break;
					}
				}
				if (next_ra) continue;

        for(diff_entry_list::iterator mc_it = mc.begin(); mc_it != mc.end(); mc_it++) //MC
        {
          diff_entry& mc_item = **mc_it;

					if (ra_item["seq_id"] != mc_item["seq_id"]) continue;

					if ( (n(ra_item["position"]) >= n(mc_item["start"])) && (n(ra_item["position"]) <= n(mc_item["end"])) )
					{
						ra_item["deleted"] = "1";
						break;
					}
				}

			}
		}

		ra.sort(MutationPredictor::sort_by_pos);

		///
		// Gather together read alignment mutations that occur next to each other
		// ...unless they are polymorphisms
		///

		bool first_time = true;
		diff_entry mut;
		vector<diff_entry> muts;

    
    for(diff_entry_list::iterator ra_it = ra.begin(); ra_it != ra.end(); ra_it++) //RA
    {
      diff_entry& item = **ra_it;

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
        mut._type = INS;
				// unused fields
				mut._fields.erase("ref_seq");
			}
			// deletion
			else if (mut["new_seq"].size() == 0)
			{
        mut._type = DEL;
				mut["size"] = s(n(mut["end"]) - n(mut["start"]) + 1);

				// unused fields
				mut._fields.erase("new_seq");
				mut._fields.erase("ref_seq");
			}
			// block substitution
			else if ((mut["ref_seq"].size() > 1) || (mut["new_seq"].size() > 1))
			{
        mut._type = SUB;
				mut["size"] = s(mut["ref_seq"].size());
				mut._fields.erase("ref_seq");
			}
			//snp
			else
			{
				mut._fields.erase("ref_seq");
        mut._type = SNP;
			}

			// we don't need these fields
			if (mut["frequency"] == "1") mut._fields.erase("frequency");
			mut._fields.erase("start");
			mut._fields.erase("end");
			mut._fields.erase("insert_start");
			mut._fields.erase("insert_end");

			gd.add(mut);
		}
		
	}


} // namespace breseq
