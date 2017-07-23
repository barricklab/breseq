/*****************************************************************************

AUTHORS

  Jeffrey E. Barrick <jeffrey.e.barrick@gmail.com>
  David B. Knoester

LICENSE AND COPYRIGHT

  Copyright (c) 2008-2010 Michigan State University
  Copyright (c) 2011-2012 The University of Texas at Austin

  breseq is free software; you can redistribute it and/or modify it under the  
  terms the GNU General Public License as published by the Free Software 
  Foundation; either version 1, or (at your option) any later version.

*****************************************************************************/

#ifndef _BRESEQ_COMMON_H_
#define _BRESEQ_COMMON_H_

// System headers

#include <config.h>

// C
#include <signal.h>
#include <execinfo.h>
#include <assert.h>
#include <libgen.h>
#include <limits.h>
#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <limits.h>
#include <sys/stat.h>
#include <sys/types.h>

#if HAVE_LIBUNWIND
  #include <libunwind.h>
#endif

// C++
// Containers
#include <list>
#include <map>
#include <set>
#include <vector>
#include <string>
// Streams
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <iomanip>
// Other
#include <cerrno>
#include <cmath>
#include <limits>
#include <memory>
#include <utility>
#include <functional>
#include <iterator>

#include "gzstream.h"

// Begin breseq specific --->
// Library specific headers
#include "bam.h"
#include "sam.h"
#include "htslib/faidx.h"

#define _base_bam_is_A(x) (x == 0x01)
#define _base_bam_is_C(x) (x == 0x02)
#define _base_bam_is_G(x) (x == 0x04)
#define _base_bam_is_T(x) (x == 0x08)
#define _base_bam_is_N(x) (x == 0x0f)
#define _base_char_is_N(x) (x == 'N')

#define UNDEFINED_UINT32 UINT_MAX
#define UNDEFINED_INT32 INT_MAX

using namespace std;

namespace breseq {
	
  // These are our own local wrappers for common functions.
  
  inline void  my_error_handler(bool condition, bool fatal, bool include_backtrace, const char *file, const char *base_file, int line, const string& message = "")
  {
    (void)base_file;
    if (condition) return;
    
    
    if (fatal) {
      
      cerr << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!> FATAL ERROR <!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
      if (message.length() > 0) cerr << message << endl;
      if (file && base_file && line) {
        cerr << "FILE: " << file << "   LINE: " << line << endl;
      }

    } else {
      cerr << "----------------------------------> WARNING <-----------------------------------" << endl;
      if (message.length() > 0) cerr << message << endl;
      if (file && base_file && line) {
        cerr << "FILE: " << file << "   LINE: " << line << endl;
      }
      
    }
    
#if HAVE_LIBUNWIND
    if (include_backtrace) {

      if (fatal) {
        cerr << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!> STACK TRACE <!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
      } else {
        cerr << "--------------------------------> STACK TRACE <---------------------------------" << endl;
      }

      unw_cursor_t    cursor;
      unw_context_t   context;
      
      unw_getcontext(&context);
      unw_init_local(&cursor, &context);
      
      while (unw_step(&cursor) > 0)
      {
        unw_word_t  offset, pc;
        char        fname[64];
        
        unw_get_reg(&cursor, UNW_REG_IP, &pc);
        
        fname[0] = '\0';
        (void) unw_get_proc_name(&cursor, fname, sizeof(fname), &offset);
        cerr << "0x" << hex << setw(16) << setfill('0') << pc << " " << fname << endl;
      }
    }
    
#else
    //Alternative fallback version which uses backtrace_symbols
    //This doesn't work with statically linked libraries.
  
    if (include_backtrace) {
      
      if (fatal) {
        cerr << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!> STACK TRACE <!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
      } else {
        cerr << "--------------------------------> STACK TRACE <---------------------------------" << endl;
      }
      
      void *array[20];
      size_t size;
      char **strings;
      size_t i;
      
      size = backtrace(array, 20);
      strings = backtrace_symbols(array, size);
      
      printf ("Backtrace with %zd stack frames.\n", size);
      
      for (i = 0; i < size; i++)
        cerr << strings[i] << endl;
      
      free (strings);
    }
#endif
    
    if (fatal) {
      cerr << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
    } else {
      cerr << "--------------------------------------------------------------------------------" << endl;
    }

    
    if (fatal) {
#ifdef DEBUG
    assert(false);
#else
    exit(0);
#endif
    }
  }
  
  inline void seg_fault_handler(int sig) {
    (void) sig;
    my_error_handler(false, true, true, NULL, NULL, 0 , "Segmentation Fault");
  }
  
  
// Fatal versions, always add backtrace
#define ASSERT(condition, message) { my_error_handler( condition, true, true, __FILE__, __BASE_FILE__, __LINE__, message); }
#define ERROR(message) { my_error_handler( false, true, true, __FILE__, __BASE_FILE__, __LINE__, message); }
  
// Nonfatal versions, add backtrace if requested
#define CHECK(condition, message) { my_error_handler( condition, false, false, __FILE__, __BASE_FILE__, __LINE__, message); }
#define WARN(message) { my_error_handler( false,  false, false, NULL, NULL, 0, message); }
#define WARN_WITH_BACKTRACE(message) { my_error_handler( false, false, true, __FILE__, __BASE_FILE__, __LINE__, message); }
  
  // There are three ways to represent a base.
  // We use typing to prevent coding errors when converting.
  // 
  // bam: SamTools representation, uses four bit codes 
  //      A(0x1), C(0x2), G(0x4), T(0x8), N(0xf), .('.')
  // &     Note that we add '.' for convenience to this list.
  typedef uint8_t base_bam;
  
  // char: Normal human-readable representation
  //       'A', 'C', 'G', 'T', 'N', '.'
  typedef char base_char;
  
  // index: Numbered starting at zero, used for array storage and lookups
  //       A(0), C(1), G(2), T(3), .(4)  No 'N' bases allowed!
  typedef uint8_t base_index;

  /*! Definition of all single-base states that are considered */
  static base_char base_char_list[] = {'A', 'C', 'G', 'T', '.'};
  static const uint8_t base_list_size = 5;
  static const uint8_t base_list_including_N_size = 6;
  static const uint8_t base_list_N_index = 5;

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
      default: ASSERT(false, "Unrecognized BAM base");
    }		
    return ' ';
	}
  
  inline base_char complement_base_char(base_char base ) {
    switch(base) {
      case 'A': return 'T';
      case 'C': return 'G';
      case 'G': return 'C';
      case 'T': return 'A';
      case '.': return '.';
      case 'N': return 'N';
      default: ASSERT(false, "Unrecognized base char");
    }
    return ' ';
  }
    
  inline base_index complement_base_index(base_index base) {
    // ascii
    switch(base) {
      case 0: return 3;
      case 1: return 2;
      case 2: return 1;
      case 3: return 0;
      case 4: return 4;
      default: return 4; //ASSERT(false, "Unrecognized base index");
    }
    return 0;
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
      default: return 'N';   //ASSERT(false, "Unrecognized BAM base");
    }		
    return ' ';
	}

  inline char baseindex2char(base_index base) {
    switch(base) {
      case 0: return 'A';
      case 1: return 'C';
      case 2: return 'G';
      case 3: return 'T';        
      case 4: return '.';
      case 5: return 'N';
      default: ASSERT(false, "Unrecognized base index");
    }
    return ' ';
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
      case 0xf: ASSERT(false, "BAM base not allowed: 0xf");
      default: ASSERT(false, "Unrecognized BAM base");
    }			
    return 0;
	}
  
  inline uint8_t basechar2index(base_char base) {
    switch(base) {
      case 'A': return 0;
      case 'C': return 1;
      case 'G': return 2;
      case 'T': return 3;
      case '.': return 4;
      case 'N': return 5;
      default: ASSERT(false, "Unrecognized base char");
    }
    return ' ';
  }
  
  // takes a string that is supposed to be bases
  inline bool is_base_sequence(const string& base_string) {
    for(size_t i=0; i<base_string.size(); i++) {
      if (!strchr( "ATCGN",  toupper(base_string[i])))
        return false;
    }
    return true;
  }
  
  // this is our version that scrubs windows line endings by default
  inline istream& getline(istream& is, string& str) {
    std::getline(is, str);
    
    if (str.size() == 0) return is;
      
    if( str[str.size()-1] == '\r') {
      str.resize(str.size()-1);
    }
    
    return is;
  }
  
  // this is our version that scrubs windows line endings by default
  inline istream& getline(istream& is, string& str, char delim) {
    std::getline(is, str, delim);
    
    if( str[str.size()-1] == '\r') {
      str.resize(str.size()-1);
    }
    
    return is;
  }
  
  // Utility functions
  inline bool file_exists(const char *filename)
  {
    ifstream ifile(filename);
    return !ifile.fail();
  }
  
  inline bool file_is_gzipped(const char * filename)
  {
    ifstream in(filename, ios::binary);
    ASSERT(in.is_open(), "Could not open file for input: " + string(filename));
    char magic[4] = {0};
    in.read(magic, sizeof(magic));
    uint8_t temp = 0x1f;
    char compare = static_cast<char>(0x1f);
    return ( magic[0] == static_cast<char>(0x1f) && magic[1] == static_cast<char>(0x8b));
  }
  
  // handles file that may or may not be
  class flexgzfstream {
  protected:
    std::iostream* m_stream;
    
  public:
    flexgzfstream() : m_stream(NULL) {}
    
    flexgzfstream(const char * file_name, std::ios::openmode mode)
    {
      if ( (mode & ios::in) && file_is_gzipped(file_name)) {
        m_stream = new iogzstream(file_name, mode);
      } else {
        m_stream = new std::fstream(file_name, mode);
      }
      ASSERT(m_stream && !m_stream->fail(), "Could not open file: " +  std::string(file_name));
    }
    
    std::iostream * get_stream() { return m_stream; }
    
    ~flexgzfstream() { if (m_stream) delete m_stream; }
  };
 
  inline void copy_file(const string& in_fn, const string& out_fn)
  {
    ifstream ifile(in_fn.c_str());
    ASSERT(!ifile.fail(), "Could not open file for input: " + in_fn);
    ofstream ofile(out_fn.c_str());
    ASSERT(!ofile.fail(), "Could not open file for output: " + out_fn);
    ofile << ifile.rdbuf();
  }
  
  inline void replace_file_contents_using_map(
    const string & in_file_name, 
    const string & out_file_name,
    map<string, string> & replacement_map)
  {
    //written by: Aaron Reba
    
    //given a map with keys as strings to look for in a file,
    //replace those strings found with the corresponding values in the map
    
    //this will only replace the first key found!!!
    
    //method for doing this:
    //1: find the largest string to look for. this will be the largest key in
    //the map.
    //2: read the file 1 character at a time, keep track of all characters in
    //a string the size of the largest key in the map.
    //3: if a key is encountered, write the replacement to the output,
    //      write all characters that won't be replaced that have already been read,
    //      and clear all characters in the string
    //4: if the last_read string was full for this iteration, write the
    //      character that was kicked out to make room for the next character
    
    ifstream in_file(in_file_name.c_str());
    ASSERT(!in_file.fail(), "Could not open file for input: " + in_file_name);
    ofstream out_file(out_file_name.c_str());
    ASSERT(!out_file.fail(), "Could not open file for output: " + out_file_name);
    
    size_t largest_key_size = 0;
    
    map<string, string>::iterator it;
    
    //1: find the largest string to look for. this will be the largest key in
    //the map.
    for (it = replacement_map.begin(); it != replacement_map.end(); it++){
      size_t key_size = (*it).first.size();
      if (key_size > largest_key_size){
        largest_key_size = key_size;
      }
    }
    
    string last_read = "";
    char single;
    bool do_write_single;
    bool do_write_last_read;
    
    //2: read the file 1 character at a time, keep track of all characters in
    //a string the size of the largest key in the map.
    while (in_file.good()){
      string not_comparison_string;
      string comparison_string;
      
      single = in_file.get();
      last_read += single;
      single = last_read[0];
      
      do_write_single = false;
      if (last_read.size() == largest_key_size + 1){
        last_read.erase(0, 1);
        do_write_single = true;
      }
      
      do_write_last_read = false;
      for (it = replacement_map.begin(); it != replacement_map.end(); it++){
        size_t key_size = (*it).first.size();
        if (key_size > last_read.size()){
          continue;
        }
        not_comparison_string = last_read.substr(0, last_read.size() - key_size);
        comparison_string = last_read.substr(last_read.size() - key_size, key_size);
        //3: if a key is encountered, write the replacement to the output,
        //      write all characters that won't be replaced that have already been read,
        //      and clear all characters in the string
        if ((*it).first == comparison_string){
          do_write_last_read = true;
          last_read = "";
          break;
        }
      }
      //4: if the last_read string was full for this iteration, write the
      //      character that was kicked out to make room for the next character
      if (do_write_single){
        out_file << single;
      }
      if (do_write_last_read){
        out_file << not_comparison_string;
        out_file << (*it).second;
      }
    }
    
    
    last_read.erase(last_read.size() - 1, 1);
    //On my system, there seems to be an end of file character.
    //This removes it.
    
    out_file << last_read;
    
    in_file.close();
    out_file.close();
  }

  inline bool file_empty(const char *filename)
  {
    struct stat filestatus;
    if (stat( filename, &filestatus ) == 0)
      return (filestatus.st_size == 0);
    
    // error getting stat, file must not exist
    return true;
  }
  inline bool file_empty(const string& filename) 
  {
    return file_empty(filename.c_str());
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

  template <typename T> struct make_vector : public vector<T>
	{
	public:
      make_vector(const T& t) { (*this)(t); }
      make_vector& operator()(const T& t) {
	        this->push_back(t);
	        return *this;
	    }
	};


	template <class T, class U> inline vector<T> get_keys(const map<T,U>& input)
	{
		vector<T> retval;
		for (typename map<T,U>::const_iterator it = input.begin(); it != input.end(); it++)
			retval.push_back(it->first);
		return retval;
	}

	//!< Split a string on a delimiter into a vector
	inline vector<string> split(
					const  string  & theString,
					const  string  & theDelimiter
	) {
    assert(theDelimiter.size() > 0); // My own ASSERT macro.

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
		CHECK(theDelimiters.size() > 0, ""); // My own ASSERT macro.
    
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
  
  //!< Split a string on any char in string of delimiters into a vector
	inline vector<string> split_on_whitespace(
                                     const  string  & theString
                                     ) {
    
    string theDelimiters = "\t \n\r";
		size_t start = 0, end = 0;
		vector<string> theStringVector;
    
    start = theString.find_first_not_of( theDelimiters, start );
    
		while (start != string::npos)
		{
			end = theString.find_first_of( theDelimiters, start );
      
			// If at end, use length=maxLength.  Else use length=end-start.
			theStringVector.push_back(
                                theString.substr(
                                                 start,
                                                 (end == string::npos) ? string::npos : end - start
                                                 )
                                );
      
      
      start = theString.find_first_not_of( theDelimiters, end );
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

	inline string chomp(const string& str)
	{
		return str.substr(0, str.find_last_not_of("\n \t")-1);
	}


	inline ostream &operator << (ostream &stream, vector<string> lhs)
	{
		stream << join(lhs, ",");
		return stream;
	}
	/*inline ostream &operator << (ostream &stream, time_t lhs)
	{
		//stream << sprintf("%04d-%02d-%02d %02d:%02d:%02d", lhs.year+1900, lhs.mon+1, lhs.mday, lhs.hour, lhs.min, lhs.sec);
		stream << ctime(&lhs);
		return stream;
	}*/
	/*istream &operator >> (istream &stream, vector<string>& lhs)
	{
		string value;
		stream >> value;
		lhs = split(value, ",");
		return stream;
	}*/
	inline istream &operator >> (istream &stream, vector<string>& rhs)
  {
    /*
     * NOTE: Can't rely on templated form to properly read in strings,
     * since it will split them on whitespace (not \n) with the operator >>.
     */
    string value = "";

    while(!stream.eof()) {
      breseq::getline(stream, value, '\n');
      rhs.push_back(value);
    }

    return stream;
  }

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
	inline string to_string (const double& t, const uint32_t precision=1, const bool use_scientific=false)
	{
		if(isnan(t)) {
			return "NA";
		} else {
			ostringstream interpreter;
			interpreter << (use_scientific ? scientific : fixed) << setprecision(precision) << t;
      
      // Not a fan of negative zeros...
      if (interpreter.str().substr(0,3) == "-0.")
      {
        interpreter.str().replace(0,3, "0.");
      }
      string return_str = interpreter.str();
      
      // Some C++ implementations may use 3 digits for the scientific notation exponent part
      if (use_scientific && (return_str[return_str.size()-3] == '0' )) {
        return_str.erase(return_str.size()-3, 1);
      }
      
			return return_str;
		}
	}
	inline string to_string (time_t& t)
	{
		return ctime(&t);
	}

  // handle bool as either TRUE/FALSE or zero/non-zero number
  // Does not handle single-character T/F correctly 
  // @JEB this is never called??
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
    ASSERT(!s.empty(), "Attempt to convert empty string");
    T t;
		istringstream iss(s);
		iss >> boolalpha >> t;
		return t;
	}
  
  inline bool is_integer(string& s, int32_t& ret_value)
  {
    char* p;
    ret_value = static_cast<int32_t>(strtol(s.c_str(), &p, 10));
    return !(*p);
  }
  
  
  //! Short aliases to conversions
  inline int32_t n(string input) { return from_string<int32_t>(input); }
  inline uint32_t un(string input) { return from_string<uint32_t>(input); }
	inline bool b(string input) { return from_string(input); }
	inline string s(int32_t input) { return to_string(input); }
  
	inline string to_upper(const string& input)
	{
		string str = input;
		transform(str.begin(), str.end(),str.begin(), (int (*)(int))std::toupper);
		return str;
	}
	
	inline string to_lower(const string& input)
  {
      string str = input;
      transform(str.begin(), str.end(),str.begin(), (int (*)(int))std::tolower);
      return str;
  }
  
  inline string double_quote(const string& input)
  {
    return "\"" + input + "\"";
  }
  
  //  Special handling of NA and INF values in doubles
  //  Used for consensus and polymorphism scores
  //  compatible with C++ limits and R representations
  //  Using this function prevents some incompatibilities
  //  In different compilers (MacOS 10.7 compilation, for example)
  inline double double_from_string(const string& input)
  {
    string compare_string = to_upper(input);
    if ((compare_string == "NA") || (compare_string == "#NA") || (compare_string == "NAN"))
      return numeric_limits<double>::quiet_NaN();
    if (compare_string == "INF")
      return numeric_limits<double>::infinity();
    if (compare_string == "-INF")
      return -numeric_limits<double>::infinity();
    
    return from_string<double>(input);
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


	inline void add_score_to_distribution(map<int32_t, int32_t> distribution_hash_ref, int32_t score)
	{
		if (distribution_hash_ref.count(score) == 0)
			distribution_hash_ref[score] = 1; // Initialize value
		else
			distribution_hash_ref[score]++;
	}
  
  // Return the path of the file without the trailing forward-slash
  inline string path_to_dirname(string file_name)
  {
    size_t found = file_name.rfind("/");
    return ((found != string::npos) ? file_name.substr(0, found) : ".");
  }
  
  inline string path_to_filename(string file_name)
  {
    size_t found = file_name.rfind("/");
    if (found == file_name.length()) return "";
    return ((found != string::npos) ? file_name.substr(found+1, string::npos) : file_name);
  }
  
  // Returns the number of bases overlapping 
  inline int32_t overlap_length(int32_t a1_start, int32_t a1_end, int32_t a2_start, int32_t a2_end)
  {
    // bigger start
    int32_t intersection_start = max(a1_start, a2_start);
    int32_t intersection_end = min(a1_end, a2_end);
    int32_t intersection_length = intersection_end - intersection_start + 1;
    
    return max(0, intersection_length);
  }
  
  
  inline string create_path(string path)
  {
    int status = umask(0);
    
    if (path.find("/") != string::npos) {
      string accumulate_path;
      vector<string> directories = split(path, "/");
      
      for (vector<string>::iterator itr = directories.begin();
           itr != directories.end(); itr ++) {
        accumulate_path.append(*itr);
        status = mkdir(accumulate_path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        accumulate_path.append("/");
      }
    } else {
      status = mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    }
    
    if (status && (errno != EEXIST))
    {
      cerr << "Could not create path: '" << path << "'" << endl;
      exit(-1);
    }
    
		return path;
	}
  
  /* Returns executable path */
  inline string getExecPath (char * argv0)
  {
    uint32_t dest_len = 512;
    char path[dest_len];
    
    char * baseName = NULL;
    char * systemPath = NULL;
    char * candidateDir = NULL;
    
    /* the easiest case: we are in linux */
    if (readlink ("/proc/self/exe", path, dest_len) != -1)
    {
      size_t len = readlink ("/proc/self/exe", path, dest_len);
      path[len] = '\0';
      return path_to_dirname(path);
    }
    
    /* Ups... not in linux, no  guarantee */
    
    /* check if we have something like execve("foobar", NULL, NULL) */
    if (argv0 == NULL)
    {
      /* we surrender and give current path instead */
      if (getcwd (path, dest_len) == NULL) return NULL;
      return path;
    }
    
    
    /* argv[0] */
    /* if dest_len < PATH_MAX may cause buffer overflow */
    if ((realpath (argv0, path)) && (!access (path, F_OK)))
    {
      return path_to_dirname(path);
    }
    
    /* Current path */
    baseName = basename (argv0);
    if (getcwd (path, dest_len - strlen (baseName) - 1) == NULL)
      return "";
    
    strcat (path, "/");
    strcat (path, baseName);
    if (access (path, F_OK) == 0)
    {
      return path_to_dirname(path);
    }
    
    /* Try the PATH. */
    systemPath = getenv ("PATH");
    if (systemPath != NULL)
    {
      dest_len--;
      systemPath = strdup (systemPath);
      for (candidateDir = strtok (systemPath, ":"); candidateDir != NULL; candidateDir = strtok (NULL, ":"))
      {
        strncpy (path, candidateDir, dest_len);
        strncat (path, "/", dest_len);
        strncat (path, baseName, dest_len);
        
        if (access(path, F_OK) == 0)
        {
          free (systemPath);
          return path_to_dirname(path);
        }
      }
      free(systemPath);
      dest_len++;
    }
    
    /* again someone has use execve: we dont know the executable name; we surrender and give instead current path */
    if (getcwd (path, dest_len - 1) == NULL) return NULL;
    return path;
  }
  
  inline vector<string> prefix_each_in_vector(const vector<string>& in_list, string prefix)
  {
    vector<string> return_list;
    for (vector<string>::const_iterator it = in_list.begin(); it != in_list.end(); it++)
    {
      return_list.push_back(prefix + *it);
    }
    return return_list;
  }
  
  
  template <typename T, typename U> vector<T> map_keys_to_list (map<T,U>& the_map)
  {
    vector<T> return_list;
    for (typename map<T,U>::const_iterator it = the_map.begin(); it != the_map.end(); it++ )
    {
      return_list.push_back(it->first);
    }
    return return_list;
  }
  
  
  template <typename T, typename U> vector<U> map_key_list_to_values (map<T,U>& the_map,vector<T> the_keys)
  {
    vector<U> return_list;
    for (typename vector<T>::iterator it = the_keys.begin(); it != the_keys.end(); it++ ) {
      return_list.push_back(the_map[*it]);
    }
    return return_list;
  }
  
  template <typename T, typename U> vector<string> map_key_list_to_values_as_strings(map<T,U>& the_map,vector<T> the_keys, string prefix = "")
  {
    vector<string> return_list;
    for (typename vector<T>::iterator it = the_keys.begin(); it != the_keys.end(); it++ ) {
      return_list.push_back(prefix + to_string(the_map[*it]));
    }
    return return_list;
  }

  
  template <typename T, typename U> inline void print_map(const map<T,U>& the_map)
  {
    for (class map<T,U>::const_iterator it = the_map.begin(); it != the_map.end(); it++ )
      cout << it->first << '=' << it->second << endl;
  }

  template <class inputiterator, class outputiterator, class predicate>
  outputiterator copy_if(inputiterator first, inputiterator last, outputiterator result, predicate pred)
  {
    inputiterator currentin = first;
    outputiterator currentout = result;
    while (currentin != last) {
      if (pred(*currentin)) {
        *currentout = *currentin;
        ++currentout;
        ++currentin;
      } else {
        ++currentin;
      }
    }
    return currentout;
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
    friend bool operator!= (const counted_ptr<X>& lhs, const counted_ptr<X>& rhs)
    {return (lhs.itsCounter->ptr == rhs.itsCounter->ptr);}
    friend bool operator== (const counted_ptr<X>& lhs, const counted_ptr<X>& rhs)
    {return (lhs.itsCounter->ptr == rhs.itsCounter->ptr);}
    X& operator*()  const throw()   {return *itsCounter->ptr;}
    X* operator->() const throw()   {return itsCounter->ptr;}
    X* get()        const throw()   {return itsCounter ? itsCounter->ptr : 0;}
    bool unique()   const throw()
    {return (itsCounter ? itsCounter->count == 1 : true);}
    
    //! compare what we point to, for convenience
    bool operator<(const counted_ptr& _in) const
    { return *(this->get())< *(_in.get()); }
    
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
          itsCounter = 0;
        }
      }
    }
  };
  
  inline string SYSTEM_CAPTURE(string command, bool silent = false)
	{
    if (!silent) cout << "[system] " << command << endl;
    
		// Open the command for reading.
    string piped_command = command + " 2>&1";
		FILE *fp = popen(piped_command.c_str(), "r");
		assert(fp != NULL);
    
		// Read the output a line at a time
		stringstream ss;
		char path[1035];
		while (fgets(path, sizeof (path) - 1, fp) != NULL)
		{
			ss << path;
		}
    
		// Close
		pclose(fp);
    
    // Delete the trailing line ending as a convenience for 'which'
    string s = ss.str();
    size_t line_break_pos = s.rfind("\n");
    if (line_break_pos != string::npos) 
      s.erase(line_break_pos);
		return s;
	}

inline int32_t fprintf(ostream &os, const char *format,...) {
  int32_t ret_val;
  char buffer[32000];

  va_list p_args;

  const size_t &size = sizeof(buffer) - 1;
  va_start(p_args, format);
  ret_val = vsnprintf(buffer, size, format, p_args);
  va_end(p_args);

  buffer[size] = '\0';

  os.write(buffer, string(buffer).size());

  return ret_val;
}

inline int32_t sprintf(string &value, const char *format,...) {
  int32_t ret_val;
  char buffer[32000];

  va_list p_args;

  const size_t &size = sizeof(buffer) - 1;
  va_start(p_args, format);
  ret_val = vsnprintf(buffer, size, format, p_args);
  va_end(p_args);

  buffer[size] = '\0';

  value = buffer;

  return ret_val;
}

template <typename T, typename U>
inline bool map_comp_second(const pair<T, U> & lhs, const pair<T, U> &rhs) { return lhs.second < rhs.second;}

template<typename T> inline void SYSTEM_CAPTURE(T out_itr, string command, bool silent = false)
{
  if (!silent) cout << "[system] " << command << endl;

  // Open the command for reading.
  string piped_command = command + " 2>&1";
  FILE *fp = popen(piped_command.c_str(), "r");
  assert(fp != NULL);

  // Read the output a line at a time
  char str[1035];
  while (fgets(str, sizeof (str) - 1, fp) != NULL)
  {
    const char *pch = strchr(str, '\n');
    if (pch != NULL) {
      str[pch-str] = '\0';
    }
    *out_itr++ = str;
  }
  // Close
  pclose(fp);
}

class cKeyValuePair : public string
{
  public:
    cKeyValuePair(const string &line, const char split_chr);

    bool   valid(void)  const;
    string get_key(void)   const;
    string get_value(void) const;

  private:
    char   _split_chr;
    string _key;
    string _value;
};

inline cKeyValuePair::cKeyValuePair(const string &line, const char split_chr)
  : string(line)
  , _split_chr(split_chr)
  , _key(line.substr(0, line.find(split_chr)))
  , _value(line.substr(line.find(split_chr) + 1))
{}

inline bool cKeyValuePair::valid(void) const
{ return count(this->begin(), this->end(), _split_chr) >= 1 && _key.size() && _value.size(); }

inline string cKeyValuePair::get_key() const
{ return _key; }

inline string cKeyValuePair::get_value() const
{ return _value; }

class cString : public string
{
  public:
    template<class T> cString(const T &val) : string(val) {}
    cString(const char *format,...);
    cString() : string ("") {}

    bool   starts_with(const string &prefix) const;
    bool   ends_with(const string &suffix) const;
    bool   contains(const char chr)const
    ;
    cString remove_ending(const string &suffix);
    cString remove_starting(const string &prefix);
    cString trim_ends_of(const char val);
    cString& escape_shell_chars(void);

    cString get_base_name() const;
    cString get_base_name_no_extension(bool remove_all_extensions = false, bool one_name_for_pair = false) const;
    cString get_base_name_unzipped() const;
    cString get_file_extension() const;
    cString get_directory_path() const;

};

inline cString::cString(const char *format,...)
{
  char buffer[32000];

  va_list p_args;

  const size_t &size = sizeof(buffer) - 1;
  va_start(p_args, format);
  vsnprintf(buffer, size, format, p_args);
  va_end(p_args);

  buffer[size] = '\0';

  *this = buffer;
}

inline bool cString::starts_with(const string &prefix) const
{
  if (this->find(prefix) == 0) {
    return true;
  } else {
    return false;
  }
}

inline bool cString::ends_with(const string &suffix) const
{
  if (this->size() < suffix.size()) {
    return false;
  } else {
    return this->substr(this->size() - suffix.size(), suffix.size()) == suffix;
  }
}

inline cString cString::remove_ending(const string &suffix)
{
  if (this->ends_with(suffix))
    this->erase(this->size() - suffix.size());
  return *this;
}

inline cString cString::remove_starting(const string &prefix)
{
  if (this->starts_with(prefix))
    this->erase(0, prefix.size());
  return *this;
}

inline cString cString::trim_ends_of(const char val)
{
  for (size_t i = 0; (*this)[i] == val; this->erase(i++, 1)) {}
  for (size_t i = this->size() - 1; (*this)[i] == val; this->erase(i--)) {}
  return *this;
}

//! Returns file name and extension, removes any directory path beforehand.
inline cString cString::get_base_name() const
{
  const size_t pos = this->rfind('/');
  if (pos == string::npos)
    return *this;
  else
    return this->substr(pos + 1);
}

//! Returns file name with no extension, removes any directory path beforehand.
//  input.1.gd
//  remove_all_extensions = FALSE removes everything past the first period  => input 
//  remove_all_extensions = TRUE removes only the last period and beyond   => input.1
inline cString cString::get_base_name_no_extension(bool remove_all_extensions, bool one_name_for_pair) const
{
  cString this_return = this->get_base_name();
  
  if (one_name_for_pair) {
    this_return = substitute(this_return, "_R1_", "_RX_");
    this_return = substitute(this_return, "_R2_", "_RX_");
  }
  
  size_t pos;
  
  if (remove_all_extensions) 
    pos = this_return.find('.');
  else
    pos = this_return.rfind('.');
  
  if (pos == string::npos)
    return this_return;
  else
    return this_return.substr(0, pos);
}
  
inline cString cString::get_base_name_unzipped() const
{
  cString this_return = this->get_base_name();
  this_return.remove_ending(".gz");
  this_return.remove_ending(".zip");
  return this_return;
}


inline cString cString::get_file_extension() const
{
  const size_t n = this->rfind('.');
  if (n != string::npos) {
    return this->substr(n);
  } else {
    return "";
  }
}

inline bool cString::contains(const char chr) const
{
  return this->find(chr) != string::npos;
}

inline cString cString::get_directory_path() const
{
  size_t pos = this->rfind('/');
  if (pos == string::npos) {
    return "./";
  } else {
    return this->substr(0, pos);;
  }
}

inline cString& cString::escape_shell_chars(void) {
  cString& value = *this;
  char escapees[] = {'<', '>', '|', '&', '\0', ';'};
  char temp[1000];

  uint32_t k = 0;
  for(uint32_t i = 0; i < value.size(); ++i) {
    char found = '\0';
    for(uint32_t j = 0; escapees[j]; ++j) {
      if (value[i] == escapees[j]) {
        found = escapees[j];
      }
    }
    if (found) {
      temp[k++] = '\\';
    }
    temp[k++] = value[i];
  }
  temp[k] = '\0'; 

  *this = temp;

  return *this;
}

//! Used to add types that will print with a specified precision
struct formatted_double {
  
  double  _value;     //actual value
  uint8_t _precision; //number of digits past zero to print
  bool _use_scientific;
  
  //! Constructor.
  formatted_double(const double v, const uint8_t p=1, const bool use_scientific=false)
  : _value(v), _precision(p), _use_scientific(use_scientific) {}
  
  virtual ~formatted_double() { }
  
  string to_string() const {
    return breseq::to_string(_value, _precision, _use_scientific);
  }
  
};

// For writing formatted_doubles
inline ostream &operator<<( ostream &out, const formatted_double &fd ) {
  out << fd.to_string();
  return out;
}
  
inline int SYSTEM(string command, bool silent = false, bool ignore_errors = false, bool escape_shell_chars = true)
{
  if (escape_shell_chars) {
    command = cString(command).escape_shell_chars();
  }
  if (!silent) cerr << "[system] " << command << endl;
  int return_value = system(command.c_str());
  
  string error_message = "Error running command:\n[system] " + command + "\nResult code: " + to_string(return_value);
  if (ignore_errors)
  {
    if (return_value != 0) cerr << error_message;
  }
  else
  {
    ASSERT(return_value == 0, error_message);
  }
  return return_value;
}
  
inline string remove_file(string path, bool silent = false)
{
  //remove(path.c_str()); // @JEB this does not work with wildcards
  SYSTEM("rm -f " + path, silent);
  return path;
}

template<uint32_t nth_place> double roundp(double value) {
  return floor(value * nth_place + 0.5f) / nth_place;
}

} // breseq

#endif
