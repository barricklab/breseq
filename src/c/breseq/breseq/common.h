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

#ifndef _BRESEQ_COMMON_H_
#define _BRESEQ_COMMON_H_

// System headers

// C
#include <assert.h>
#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// C++
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <list>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

// Library specific headers
#include <bam.h>
#include <sam.h>
#include <faidx.h>
#include <boost/any.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/optional.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/variant.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>

// Begin breseq specific --->

#define _base_bam_is_A(x) (x == 0x01)
#define _base_bam_is_C(x) (x == 0x02)
#define _base_bam_is_G(x) (x == 0x04)
#define _base_bam_is_T(x) (x == 0x08)
#define _base_bam_is_N(x) (x == 0x0f)
#define _base_char_is_N(x) (x == 'N')

using namespace std;

namespace breseq {
	
  
  // There are three ways to represent a base.
  // We use typing to prevent coding errors when converting.
  // 
  // bam: SamTools representation, uses four bit codes 
  //      A(0x1), C(0x2), G(0x4), T(0x8), N(0xf), .('.')
  //      Note that we add '.' for convenience to this list.
  typedef uint8_t base_bam;
  
  // char: Normal human-readable representation
  //       'A', 'C', 'G', 'T', 'N', '.'
  typedef char base_char;
  
  // index: Numbered starting at zero, used for array storage and lookups
  //       A(0), C(1), G(2), T(3), .(4)  No 'N' bases allowed.
  typedef uint8_t base_index;

  /*! Definition of all single-base states that are considered */
  static base_char base_char_list[] = {'A', 'C', 'G', 'T', '.'};
  static uint8_t base_list_size = 5;

  
	struct Settings
	{
		// Fields

		map<string, bool> installed;
		string bin_path;

		string candidate_junction_score_method;

		string candidate_junction_fasta_file_name;
		string candidate_junction_faidx_file_name;
		string candidate_junction_sam_file_name;
		string jc_genome_diff_file_name;
		string preprocess_junction_best_sam_file_name;
		string preprocess_junction_split_sam_file_name;
		string reference_fasta_file_name;
		string reference_faidx_file_name;
		string reference_sam_file_name;
		string resolved_reference_sam_file_name;
		string resolved_junction_sam_file_name;
		string unmatched_read_file_name;

		bool no_junction_prediction;
		bool unmatched_reads;
		bool add_split_junction_sides;
		bool require_complete_match;

		int32_t alignment_read_limit;
		int32_t candidate_junction_read_limit;
		int32_t minimum_candidate_junction_pos_hash_score;
		int32_t minimum_candidate_junction_min_overlap_score;
		int32_t minimum_candidate_junctions;
		int32_t maximum_candidate_junctions;
		int32_t maximum_candidate_junction_length_factor;
		int32_t max_read_length;
		int32_t maximum_inserted_junction_sequence_length;
		int32_t maximum_read_mismatches;
		int32_t required_both_unique_length_per_side;
		int32_t required_one_unique_length_per_side;
		int32_t required_extra_pair_total_length;
		int32_t required_match_length;

		boost::optional<int32_t> preprocess_junction_min_indel_split_length;

		struct ReadStructure
		{
			string base_name;
		};
		vector<ReadStructure> read_structures;

		// Utility function to substitute specific details into a generic file name
		static string file_name(string file_name_key)
		{
			return file_name_key;

//			my ($self, $file_name_key, $sub_hash)= @_;
//			my $file_name = $self->{$file_name_key};
//			$file_name or $self->throw("Settings file \"$file_name_key\" not found.");
//
//			return $self->substitute_file_name($file_name, $sub_hash);
		}

		string ctool(string tool_name)
		{
//			my ($self, $tool_name, $allow_fail) = @_;
//
//			if (!$self->{installed}->{$tool_name})
//			{
//				if ($allow_fail)
//				{
//					$self->warn("Executable \"$tool_name\" not found in breseq bin path\"$self->{bin_path}\".");
//					return undef; # couldn't find it, but it's not an error.
//				}
//				else
//				{
//					$self->throw("Executable \"$tool_name\" not found in breseq bin path\"$self->{bin_path}\".");
//				}
//			}

			return bin_path + "/" + tool_name;
		}
	};

	struct Summary
	{
		// Fields

		struct AlignmentCorrection
		{
			string read_file;
			struct NewJunction
			{
				int32_t observed_min_overlap_score_distribution;
				int32_t accepted_min_overlap_score_distribution;
				int32_t observed_pos_hash_score_distribution;
				int32_t accepted_pos_hash_score_distribution;
			};
			list<NewJunction> new_junctions;

		} alignment_correction;

		struct PreprocessCoverage
		{
			int32_t junction_accept_score_cutoff_1;
			int32_t junction_accept_score_cutoff_2;
		} preprocess_coverage;

		struct CandidateJunctionSummaryData
		{
			struct Total
			{
				int32_t number;
				int32_t length;
			} total;

			struct Accepted
			{
				int32_t number;
				int32_t length;
				int32_t pos_hash_score_cutoff;
				int32_t min_overlap_score_cutoff;
			} accepted;

			map<int32_t, int32_t> pos_hash_score_distribution;
			map<int32_t, int32_t> min_overlap_score_distribution;

			map<string, map<string, int32_t> > read_file;
		} candidate_junction;

		struct SequenceConversion
		{
			int32_t total_reference_sequence_length;
		} sequence_conversion;
	};

	/*! Reverse a base.
	 */
	inline base_bam complement_base_bam(base_bam base) {
    // sam-style 4-bit field
    switch(base) {
      case 0x1: return 0x8;
      case 0x2: return 0x4;
      case 0x4: return 0x2;
      case 0x8: return 0x1;        
      case 0xf: return 0xf;
      case '.': return '.';
      default: assert(false);
    }		
	}
  
  inline base_char complement_base_char(base_char base ) {
    switch(base) {
      case 'A': return 'T';
      case 'C': return 'G';
      case 'G': return 'C';
      case 'T': return 'A';
      case '.': return '.';
      case 'N': return 'N';
      default: assert(false);
    }
  }
  
  inline base_index complement_base_index(base_index base) {
    // ascii
    switch(base) {
      case 0: return 3;
      case 1: return 2;
      case 2: return 1;
      case 3: return 0;
      case 4: return 4;
      default: assert(false);
    }
  }

	
	/*! Convert a base to an ASCII character.
	 */
	inline char basebam2char(base_bam base) {
    // sam-style 4-bit field
    switch(base) {
      case 0x01: return 'A';
      case 0x02: return 'C';
      case 0x04: return 'G';
      case 0x08: return 'T';
      case '.': return '.';
      case 0x0f: return 'N'; // might want to assert here 
      default: assert(false);
    }		
	}

  inline char baseindex2char(base_index base) {
    switch(base) {
      case 0: return 'A';
      case 1: return 'C';
      case 2: return 'G';
      case 3: return 'T';        
      case 4: return '.';
      default: assert(false);
    }
  }
 
 	/*! Convert a base to an index
	 */
	inline char basebam2index(base_bam base) {
    switch(base) {
      case 0x1: return 0;
      case 0x2: return 1;
      case 0x4: return 2;
      case 0x8: return 3;  
      case '.': return 4;        
      case 0xf: assert(false);
      default: assert(false);
    }			
	}
  
  inline char basechar2index(base_char base) {
    switch(base) {
      case 'A': return 0;
      case 'C': return 1;
      case 'G': return 2;
      case 'T': return 3;
      case '.': return 4;
      default: assert(false);
    }
  }

  // Utility functions
  inline bool file_exists(const char *filename)
  {
    ifstream ifile(filename);
    return ifile;
  }
  
  inline bool file_empty(const char *filename)
  {
    ifstream ifile(filename);
    string test;
    getline(ifile, test);
    return !ifile.eof();
  }
  
  //!< Split a string on a delimiter into a vector
  inline void
  split( vector<string> & theStringVector,  /* Altered/returned value */
        const  string  & theString,
        const  string  & theDelimiter )
  {
    assert( theDelimiter.size() > 0 ); // My own ASSERT macro.
    
    size_t  start = 0, end = 0;
    
    while ( end != string::npos )
    {
      end = theString.find( theDelimiter, start );
      
      // If at end, use length=maxLength.  Else use length=end-start.
      theStringVector.push_back( theString.substr( start,
                                                  (end == string::npos) ? string::npos : end - start ) );
      
      // If at end, use start=maxSize.  Else use start=end+delimiter.
      start = (   ( end > (string::npos - theDelimiter.size()) )
               ?  string::npos  :  end + theDelimiter.size()    );
    }
  }

	inline vector<bam1_t*> tam_next_read_alignments(tamFile tam, bam_header_t* header, bam1_t* last_alignment, bool paired = false)
	{
		int num_to_slurp = (paired) ? 2 : 1;
		string last_read_name;
		vector<bam1_t*> al_ref;
		if (last_alignment != NULL)
		{
			last_read_name = bam1_qname(last_alignment);
			al_ref.push_back(last_alignment);
			last_alignment = NULL;
		}

		int num_slurped = 0;
		while (true)
		{
			last_alignment = new bam1_t();
			int bytes = sam_read1(tam, header, last_alignment);

			//returns bytes == -1 if EOF reached
			if (bytes < 0)
			{
				last_alignment = NULL;
				return al_ref;
			}

			string read_name = bam1_qname(last_alignment);

			if ( (last_read_name.size() > 0) && (read_name != last_read_name) && (++num_slurped == num_to_slurp) )
				break;

			if (last_read_name.size() == 0)
				last_read_name = read_name;

			al_ref.push_back(last_alignment);
		}

		return al_ref;
	}

	inline uint32_t fix_flags(uint32_t flags)
	{
		flags = ((flags >> 9) << 9) + flags % 128;
		return flags;
	}

	template <class T> inline string to_string (const T& t)
	{
		stringstream ss;
		ss << t;
		return ss.str();
	}

	inline string join(const vector<string>& values, const string& separator)
	{
		if(values.size() == 0)
			return "";

		string::size_type size = separator.length() * values.size();
		for(uint32_t i=0; i < values.size(); i++)
			size += values[i].size();

		string retval;
		retval.reserve(size);
		retval = values[0];
		for(uint32_t i = 1; i < values.size(); i++)
			retval += separator + values[i];

		return retval;
	}

	inline string join(string values[], const string& separator)
	{
		return join(vector<string> (values, values + sizeof(values) / sizeof(*values)), separator);
	}

	struct Trim
	{
		string L;
		string R;
	};

	inline void tam_write_read_alignments(ofstream& fh, bam_header_t* header, int32_t fastq_file_index, vector<bam1_t*> al, boost::optional< vector<Trim> > trims)
	{
		for (int32_t i = 0; i < al.size(); i++)
		{
			bam1_t* a = al[i];

			stringstream aux_tags_ss;
			aux_tags_ss << "AS:i:" << bam_aux_get(a, "AS") << "\t" << "X1:i:" << al.size() << "\t" << "X2:i:" << fastq_file_index;

			if (trims && trims.get().size() > i)
			{
				Trim trim = trims.get()[i];
				aux_tags_ss << "\t" << "XL:i:" << trim.L << "\t" << "XR:i:" << trim.R;
			}

			string aux_tags = aux_tags_ss.str();

			string* qscore = (string*)bam1_qual(a);
			string quality_score_string = *qscore;
			for (int32_t j = 0; j < quality_score_string.size(); j++)
				quality_score_string[j] = quality_score_string[j] + 33;

			uint32_t* cigar_list = bam1_cigar(a);
			stringstream cigar_string_ss;

			for (int32_t j = 0; j <= a->core.n_cigar; j++) //foreach my $c (@$cigar_list)
			{
				uint32_t op = cigar_list[i] & BAM_CIGAR_MASK;
				uint32_t len = cigar_list[i] >> BAM_CIGAR_SHIFT;
				cigar_string_ss << len << op; //$cigar_string += $c->[1] + $c->[0];
			}
			string cigar_string = cigar_string_ss.str();

			vector<string> ll;
			ll.push_back(bam1_qname(a));
			ll.push_back(to_string(fix_flags(a->core.flag)));
			ll.push_back(header->target_name[a->core.tid]);
			ll.push_back(to_string(a->core.pos));
			ll.push_back(to_string(a->core.qual));
			ll.push_back(cigar_string);

			//something strange in new version... such that mate_start sometimes
			//returns 1 even though there is no mate
			if (a->core.flag & BAM_FPROPER_PAIR != 0)
			{
				ll.push_back("*");
				ll.push_back(0);
				ll.push_back(0);
			}
			else
			{
				ll.push_back("=");
				ll.push_back(0);
				int32_t mate_start = a->core.mpos + 1;
				ll.push_back(to_string(mate_start));
				ll.push_back(to_string(a->core.isize));
			}

			ll.push_back(to_string(*bam1_seq(a)));
			ll.push_back(quality_score_string);
			ll.push_back(aux_tags);

			fh << join(ll, "\t") << endl;
		}
	}

	// returns start and end in coordinates of query
	//lowest query base first, regardless of which strand
	// it matched in the reference genome (i.e., reversed alignment).
	inline void alignment_query_start_end(bam1_t* a, int32_t& start, int32_t& end)
	{
		uint32_t* cigar = bam1_cigar(a); // cigar array for this alignment
		start = 1;
		end = bam_cigar2qlen(&a->core,cigar);

		// start:
		for(uint32_t i=0; i<= a->core.n_cigar; i++)
		{
		    uint32_t op = cigar[i] & BAM_CIGAR_MASK;
		    uint32_t len = cigar[i] >> BAM_CIGAR_SHIFT;
		    // if we encounter padding, or a gap in reference then we are done
		    if((op != BAM_CSOFT_CLIP) && (op != BAM_CHARD_CLIP) && (op != BAM_CREF_SKIP)) {
				break;
		    }
		    start += len;
		}

		// end:
		for(uint32_t i=(a->core.n_cigar-1); i>0; --i)
		{
			uint32_t op = cigar[i] & BAM_CIGAR_MASK;
			uint32_t len = cigar[i] >> BAM_CIGAR_SHIFT;
			// if we encounter padding, or a gap in reference then we are done
			if((op != BAM_CSOFT_CLIP) && (op != BAM_CHARD_CLIP) && (op != BAM_CREF_SKIP)) {
			  break;
			}
			end -= len;
		}
	}

	inline uint32_t alignment_query_length(bam1_t* a) {
		uint32_t* cigar = bam1_cigar(a); // cigar array for this alignment
		uint32_t qlen = bam_cigar2qlen(&a->core, cigar); // total length of the query
		return qlen;
	}

	inline bool is_reversed(bam1_t* a)
	{
		return bam1_strand(a);
	}

	inline string bama_qseq(bam1_t* a)
	{
	    string seq(a->core.l_qseq, ' ');
	    for (int32_t i = 0; i < a->core.l_qseq; i++)
			seq[i] = bam_nt16_rev_table[bam1_seqi(bam1_seq(a),i)];
	}

	inline string reverse_complement(string seq)
	{
		char trade['Z'];
		trade['A'] = 'T'; trade['T'] = 'A'; trade['C'] = 'G'; trade['G'] = 'C';
		string retval = seq;
		for (int i = 0; i < seq.size(); i++)
			retval[i] = trade[seq[seq.size() - 1 - i]];
		return retval;
	}

	struct JunctionList
	{
		struct Side
		{
			string seq_id;
			int32_t position;
			bool strand;
			int32_t redundant;
		} side_1, side_2;

		int32_t alignment_overlap;
		string unique_read_sequence;
		int32_t flanking_left;
		int32_t flanking_right;
	};

	const string junction_name_separator = "__";

	// Deserializes a JunctionList from a string
	inline JunctionList junction_name_split(string junction_name)
	{
		vector<string> s;
		split(s, junction_name, junction_name_separator);

		JunctionList::Side side_1 = {
			s[0],
			boost::lexical_cast<int32_t>(s[1]),
			boost::lexical_cast<bool>(s[2]),
			boost::lexical_cast<int32_t>(s[10])
		}, side_2 = {
			s[3],
			boost::lexical_cast<int32_t>(s[4]),
			boost::lexical_cast<bool>(s[5]),
			boost::lexical_cast<int32_t>(s[11])
		};
		JunctionList retval =
		{
			side_1,
			side_2,
			boost::lexical_cast<int32_t>(s[6]),
			s[7],
			boost::lexical_cast<int32_t>(s[8]),
			boost::lexical_cast<int32_t>(s[9])
		};

		return retval;
	}

	// Serializes a JunctionList to a string
	inline string junction_name_join(JunctionList item)
	{
		bool has_redundant = (item.side_1.redundant >= 0 && item.side_2.redundant >= 0);
		vector<string> values(has_redundant ? 12 : 10);

		values.push_back(item.side_1.seq_id);
		values.push_back(boost::lexical_cast<string>(item.side_1.position));
		values.push_back(boost::lexical_cast<string>(item.side_1.strand));

		values.push_back(item.side_2.seq_id);
		values.push_back(boost::lexical_cast<string>(item.side_2.position));
		values.push_back(boost::lexical_cast<string>(item.side_2.strand));

		values.push_back(boost::lexical_cast<string>(item.alignment_overlap));
		values.push_back(boost::lexical_cast<string>(item.unique_read_sequence));
		values.push_back(boost::lexical_cast<string>(item.flanking_left));
		values.push_back(boost::lexical_cast<string>(item.flanking_right));

		if (has_redundant)
		{
			values.push_back(boost::lexical_cast<string>(item.side_1.redundant));
			values.push_back(boost::lexical_cast<string>(item.side_2.redundant));
		}

		return join(values, junction_name_separator);
	}

	inline void add_score_to_distribution(map<int32_t, int32_t> distribution_hash_ref, int32_t score)
	{
		if (distribution_hash_ref.count(score) == 0)
			distribution_hash_ref[score] = 1; // Initialize value
		else
			distribution_hash_ref[score]++;
	}

} // breseq

#endif
