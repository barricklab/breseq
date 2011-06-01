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
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

// Library specific headers
#include <bam.h>
#include <sam.h>
#include <faidx.h>
#include <boost/optional.hpp>

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

	inline vector<bam1_t*> tam_next_read_alignments(tamFile tam, bam_header_t* header, bam1_t* last_alignment, bool paired)
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

			if (read_name != last_read_name && ++num_slurped == num_to_slurp)
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

	inline string str_join(const vector<string>& vec,const string& sep)
	{
		if(vec.size()==0)
			return "";

		string::size_type size=sep.length()*vec.size();
		for(unsigned int i=0;i<vec.size();i++)
			size+=vec[i].size();

		string tmp;
		tmp.reserve(size);
		tmp=vec[0];
		for(unsigned int i=1;i<vec.size();i++)
			tmp += sep + vec[i];

		return tmp;
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
			uint32_t cigar_list_size = a->core.n_cigar;

			for (int32_t j = 0; j < cigar_list_size; j++) //foreach my $c (@$cigar_list)
				cigar_string_ss << cigar_list[j]; //$cigar_string += $c->[1] + $c->[0];
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

			fh << str_join(ll, "\t") << endl;
		}
	}

} // breseq

#endif
