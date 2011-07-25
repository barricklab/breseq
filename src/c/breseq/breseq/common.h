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
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <list>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

// Library specific headers
#include <bam.h>
#include <sam.h>
#include <faidx.h>

// Breseq
#include "breseq/settings.h"

// Begin breseq specific --->

#define _base_bam_is_A(x) (x == 0x01)
#define _base_bam_is_C(x) (x == 0x02)
#define _base_bam_is_G(x) (x == 0x04)
#define _base_bam_is_T(x) (x == 0x08)
#define _base_bam_is_N(x) (x == 0x0f)
#define _base_char_is_N(x) (x == 'N')

#define UNDEFINED UINT_MAX
#define is_defined(x) (x != UINT_MAX)

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
    
    
    inline void donothing(){
        
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

	inline uint32_t fix_flags(uint32_t flags)
	{
		flags = ((flags >> 9) << 9) + flags % 128;
		return flags;
	}

	template <typename T, typename U> struct make_map : public map<T,U>
	{
	public:
		make_map(const T& key, const U& val) { (*this)(key, val); }
		make_map<T, U>& operator()(const T& key, const U& val)
		{
			(*this)[key] = val;
			return *this;
		}
	};

	template <typename T> struct make_list : public vector<T>
	{
	public:
	    make_list(const T& t) { (*this)(t); }
	    make_list& operator()(const T& t) {
	        this->push_back(t);
	        return *this;
	    }
	};


	template <class T, class U> inline vector<T> get_keys(map<T,U> input)
	{
		vector<T> retval;
		for (class map<T,U>::iterator it = input.begin(); it != input.end(); it++)
			retval.push_back(it->first);
		return retval;
	}

	//!< Split a string on a delimiter into a vector
	inline vector<string> split(
					const  string  & theString,
					const  string  & theDelimiter
	) {
		assert(theDelimiter.size() > 0); // My own ASSERT macro.
    assert(theString.size() != 1); //@GRC check that programmer hasn't switched variables in declaration.

		size_t start = 0, end = 0;
		vector<string> theStringVector;

    if (theString.size() == 0) return theStringVector;
    // return empty list if the string is empty
    
		while (end != string::npos)
		{
			end = theString.find( theDelimiter, start );

			// If at end, use length=maxLength.  Else use length=end-start.
			theStringVector.push_back(
					theString.substr(
						start,
						(end == string::npos) ? string::npos : end - start
					  )
				  );

			// If at end, use start=maxSize.  Else use start=end+delimiter.
			start =
					(end > (string::npos - theDelimiter.size()))
					? string::npos
					: end + theDelimiter.size();
		}
		return theStringVector;
	}
  
  //!< Split a string on any char in string of delimiters into a vector
	inline vector<string> split_on_any(
                              const  string  & theString,
                              const  string  & theDelimiters
                              ) {
		assert(theDelimiters.size() > 0); // My own ASSERT macro.
    
		size_t start = 0, end = 0;
		vector<string> theStringVector;
    
		while (end != string::npos)
		{
			end = theString.find_first_of( theDelimiters, start );
      
			// If at end, use length=maxLength.  Else use length=end-start.
			theStringVector.push_back(
        theString.substr(
           start,
           (end == string::npos) ? string::npos : end - start
           )
        );
      
			// If at end, use start=maxSize.  Else use start=end+delimiter.
			start =
      (end > (string::npos - 1))
      ? string::npos
      : end + 1;
		}
		return theStringVector;
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

	inline string chomp(const string& str)
	{
		return str.substr(0, str.find_last_not_of("\n \t")-1);
	}


	inline ostream &operator << (ostream &stream, vector<string> lhs)
	{
		stream << join(lhs, ",");
		return stream;
	}
	/*istream &operator >> (istream &stream, vector<string>& lhs)
	{
		string value;
		stream >> value;
		lhs = split(value, ",");
		return stream;
	}*/
	template <typename T> inline istream &operator >> (istream &stream, vector<T>& rhs)
	{
		rhs.clear();
		string value;
    // the different values are separated by newlines "\n" in the string
    // so they are read out correctly this way! @JEB
    while (!stream.eof())
    {
      T t;
      stream >> boolalpha >> t;
      rhs.push_back(t);
		}
		return stream;
	}


	template <typename T> inline string to_string (const T& t)
	{
		stringstream ss;
		ss << t;
		return ss.str();
	}
	inline string to_string (const pair<int,int>& t)
	{
		return to_string(t.first) + '/' + to_string(t.second);
	}
	inline string to_string (const double& t, const uint32_t precision=1)
	{
		if(isnan(t)) {
			return "NA";
		} else {
			ostringstream interpreter;
			interpreter << fixed << setprecision(precision) << t;
			return interpreter.str();
		}
	}

  // handle bool as either TRUE/FALSE or zero/non-zero number
  // Does not handle single-character T/F correctly 
  // @JEB this is never called
  inline bool from_string(const string& s)
  {
    bool t = false;
    istringstream iss1(s);
    iss1 >> boolalpha >> t;
    
    int32_t t2;
    istringstream iss2(s);
    iss2 >> noboolalpha >> t2;
    t = t || (t2 != 0);
    
    return t;
  }
  
  template <typename T> inline T from_string(const string &s)
	{
    assert(!s.empty());
    T t;
		istringstream iss(s);
		iss >> boolalpha >> t;
		return t;
	}

	inline string to_upper(const string& input)
	{
		string str = input;
		transform(str.begin(), str.end(),str.begin(), ::toupper);
		return str;
	}
	
	inline string to_lower(const string& input)
  {
      string str = input;
      transform(str.begin(), str.end(),str.begin(), ::tolower);
      return str;
  }

  inline string reverse_string(const string &in_string)
  {
    string rev_string("");
    for(string::const_reverse_iterator rit=in_string.rbegin(); rit<in_string.rend(); rit++)
    {
      rev_string+=*rit;
    }
    return rev_string;
  }
  
  inline string repeat_char(const char char_to_append, const uint32_t num_times)
  {
    string s;
    // reserve memory to be a little more efficient
    s.reserve(num_times);
    for (uint32_t i=0; i< num_times; i++)
    {
      s += char_to_append;
    }
    return s;
  }
  
  inline string substitute(const string& in_s, const string& replace_this, const string& with_this)
  {
    string s = in_s;
    size_t pos = s.find(replace_this);
    while (pos != string::npos)
    {
      s.replace(pos, replace_this.size(), with_this);
      pos += with_this.size();
      pos = s.find(replace_this, pos);
    }
    return s;
  }

	inline string reverse_complement(string seq)
	{
		char trade[static_cast<uint8_t>('Z')];
		trade[static_cast<uint8_t>('A')] = 'T'; 
    trade[static_cast<uint8_t>('T')] = 'A'; 
    trade[static_cast<uint8_t>('C')] = 'G'; 
    trade[static_cast<uint8_t>('G')] = 'C';
		string retval = seq;
		for (uint32_t i = 0; i < seq.size(); i++)
			retval[i] = trade[static_cast<uint8_t>(seq[seq.size() - 1 - i])];
		return retval;
	}

	struct CandidateJunction
	{
		int32_t r1;
		int32_t r2;
		int32_t L1;
		int32_t L2;
		int32_t min_overlap_score;
		int32_t pos_hash_score;
		map<uint32_t, uint32_t> read_begin_hash;

		struct TestInfo {
			uint32_t max_left;
			uint32_t max_left_minus;
			uint32_t max_left_plus;
			uint32_t max_right;
			uint32_t max_right_minus;
			uint32_t max_right_plus;
			uint32_t max_min_right;
			uint32_t max_min_right_minus;
			uint32_t max_min_right_plus;
			uint32_t max_min_left;
			uint32_t max_min_left_minus;
			uint32_t max_min_left_plus;
			uint32_t coverage_minus;
			uint32_t coverage_plus;
			uint32_t total_non_overlap_reads;
			uint32_t min_overlap_score;
			uint32_t pos_hash_score;
		} test_info;

		struct Sorter {
			bool operator() (const string& lhs, const string& rhs) const { return (lhs < rhs); }
		};
	};

	struct Feature
	{
		string type;
		string name;
		uint32_t start;
		uint32_t end;
		bool strand;
		string product;
		string pseudogene;
		string cds;
		string note;

		string interval;
		string side_key;
	};

	struct JunctionInfo
	{
		struct Side
		{
			string seq_id;
			int32_t position;
			int32_t strand;
			int32_t redundant;

			// Extended properties for resolve_alignments.cpp
			Feature is;
			bool read_side;
			int32_t overlap;
		};
		Side sides[2];

		int32_t alignment_overlap;
		string unique_read_sequence;
		int32_t flanking_left;
		int32_t flanking_right;

		// Extended properties for resolve_alignments.cpp
		string key;
		int32_t overlap;
		uint32_t total_reads;
		uint32_t unique_side;
		uint32_t is_side;
	};

	const string junction_name_separator = "__";

	// Deserializes a JunctionInfo from a string
	inline JunctionInfo junction_name_split(string junction_name)
	{
		vector<string> s = split(junction_name, junction_name_separator);

		JunctionInfo::Side side_1 = {
			s[0],
			from_string<int32_t>(s[1]),
			from_string<int32_t>(s[2]),
			from_string<int32_t>(s[10])
		};
    if (side_1.strand == 0) side_1.strand = -1;
    
    JunctionInfo::Side side_2 = {
			s[3],
			from_string<int32_t>(s[4]),
			from_string<int32_t>(s[5]),
			from_string<int32_t>(s[11])
		};
    if (side_2.strand == 0) side_2.strand = -1;
    
		JunctionInfo retval =
		{
			{ side_1, side_2 },
			from_string<int32_t>(s[6]),
			s[7],
			from_string<int32_t>(s[8]),
			from_string<int32_t>(s[9])
		};

		return retval;
	}

	// Serializes a JunctionInfo to a string
	inline string junction_name_join(const JunctionInfo& item)
	{
		bool has_redundant = (item.sides[0].redundant >= 0 && item.sides[1].redundant >= 0);
    
    // allocate vector of correct size
		vector<string> values(has_redundant ? 12 : 10); 

    // place values in vector
		values[0] = item.sides[0].seq_id;
		values[1] = to_string(item.sides[0].position);
		values[2] = to_string(item.sides[0].strand);

		values[3] = item.sides[1].seq_id;
		values[4] = to_string(item.sides[1].position);
		values[5] = to_string(item.sides[1].strand);

		values[6] = to_string(item.alignment_overlap);
		values[7] = to_string(item.unique_read_sequence);
		values[8] = to_string(item.flanking_left);
		values[9] = to_string(item.flanking_right);

		if (has_redundant)
		{
			values[10] = to_string(item.sides[0].redundant);
			values[11] = to_string(item.sides[1].redundant);
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
	
	//! Returns true if exp matches anywhere in input
	//! true regex expressions ie !^[\s](\w)+ do not work.
  inline	bool regex_m(string exp, string input)
  {
    if(input.find(exp) != string::npos)
      return true;
    else
      return false;
  }
  
  //! Applies regex_m to a vector of strings, returns a vector of strings
  //! that do/don't (dependent on bool match) contain an exp.
  //! true regex expressions ie !^[\s](\w)+ do not work.
  inline vector<string> grep(bool match, string exp, vector<string> lines)
  {
    typedef vector<string> Lines;
    Lines matching_lines;
    
    for(Lines::iterator line = lines.begin();
        line != lines.end(); line ++)
        if(match == regex_m(exp, (*line)))
          matching_lines.push_back((*line));
        return matching_lines;    
  }
  ///! Returns first element, then removes it from container.
  template <typename T> inline T shift(vector<T>& input)
  {
    assert(!input.empty());
    class vector<T>::iterator first = input.begin();
    T retval = (*first);
    input.erase(first);
    
    return retval;
  }
  
  //! 1       2       3       4       5       6       7
  //! splice(list, 3, 2);
  //! 1       2       3       6       7
  inline vector<string> splice (vector<string> &input, const uint8_t &offset, const uint8_t &length)
  {
    vector<string>::iterator itr = input.begin();
    vector<string>::iterator start = itr + offset;
    vector<string>::iterator end = start + length;
    
    input.erase(start,end);
    
    return input;
  }

  
// counted_ptr keeps track of number of references 
  
  template <class X> class counted_ptr
  {
  public:
    typedef X element_type;
    
    explicit counted_ptr(X* p = 0) // allocate a new counter
    : itsCounter(0) {if (p) itsCounter = new counter(p);}
    ~counted_ptr()
    {release();}
    counted_ptr(const counted_ptr& r) throw()
    {acquire(r.itsCounter);}
    counted_ptr& operator=(const counted_ptr& r)
    {
      if (this != &r) {
        release();
        acquire(r.itsCounter);
      }
      return *this;
    }
    
    X& operator*()  const throw()   {return *itsCounter->ptr;}
    X* operator->() const throw()   {return itsCounter->ptr;}
    X* get()        const throw()   {return itsCounter ? itsCounter->ptr : 0;}
    bool unique()   const throw()
    {return (itsCounter ? itsCounter->count == 1 : true);}
    
  private:
    
    struct counter {
      counter(X* p = 0, unsigned c = 1) : ptr(p), count(c) {}
      X*          ptr;
      unsigned    count;
    }* itsCounter;
    
    void acquire(counter* c) throw()
    { // increment the count
      itsCounter = c;
      if (c) ++c->count;
    }
    
    void release()
    { // decrement the count, delete if it is 0
      if (itsCounter) {
        if (--itsCounter->count == 0) {
          delete itsCounter->ptr;
          delete itsCounter;
          itsCounter = 0; //@JEB edit
        }
//@JEB -- this is an error?        itsCounter = 0;
      }
    }
  };
  

} // breseq

#endif
