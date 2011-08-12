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

#ifndef _BRESEQ_GENOME_DIFF_H_
#define _BRESEQ_GENOME_DIFF_H_

#include "common.h"



namespace breseq {

// Common keywords used for diff entries:
extern const char* TYPE;
extern const char* ID;
extern const char* PID;
extern const char* SEQ_ID;
extern const char* START;
extern const char* END;
extern const char* START_RANGE;
extern const char* END_RANGE;
extern const char* LEFT_OUTSIDE_COV;
extern const char* LEFT_INSIDE_COV;
extern const char* RIGHT_INSIDE_COV;
extern const char* RIGHT_OUTSIDE_COV;
extern const char* POSITION;
extern const char* INSERT_POSITION;
extern const char* QUALITY;
extern const char* POLYMORPHISM_QUALITY;
extern const char* REF_BASE;
extern const char* NEW_BASE;
extern const char* FREQUENCY;
extern const char* REJECT;
extern const char* REF_COV;
extern const char* NEW_COV;
extern const char* TOT_COV;
extern const char* ERROR;
//For JC
extern const char* _SIDE_1_SEQ_ID;
extern const char* _SIDE_1_POSITION;
extern const char* _SIDE_1_STRAND;
extern const char* _SIDE_2_SEQ_ID;
extern const char* _SIDE_2_POSITION;
extern const char* _SIDE_2_STRAND;
extern const char* _SIDE_KEY_JC;



// Types of diff entries:
extern const char* SNP;
extern const char* SUB;
extern const char* DEL;
extern const char* INS;
extern const char* MOB;
extern const char* AMP;
extern const char* INV;
extern const char* CON;

// Types of diff entries:
extern const char* RA;
extern const char* MC;
extern const char* JC;
extern const char* UN;

//!  
static uint8_t kPolymorphismFrequencyPrecision = 4;
static uint8_t kMutationQualityPrecision = 4;

//! Convenience typedef, used during diff entry marshalling.
typedef vector<string> field_list_t;

//! Used to add types that will print with a specified precision
struct formatted_double {

  //! Constructor.
  formatted_double(const double v, const uint8_t p=1)
    : _value(v), _precision(p) {}

  virtual ~formatted_double() { }

  string to_string() {
    return breseq::to_string(_value, _precision);
  }

  double  _value;     //actual value
  uint8_t _precision; //number of digits past zero to print
};


/*! Genome diff entry type.
 
 Instead of trying to define (and maintain!) concrete classes for each different 
 type of diff entry, we instead define a single generic diff entry, and then 
 provide easy-to-use mechanisms to manipulate them.
 
 The key abstraction here is that a diff entry is a map of keys (strings) to
 values.  In this case, values are strings that represent int, double, string,
 formatted_double and pair<,> types.
 
 For convenience, there are also factory methods on genome_diff that create diff
 entries that have been pre-populated according to their type.
 */
struct diff_entry {
  typedef string key_t; //!< Diff entry keys.
  typedef string value_t; //!< Diff entry values.
  typedef map<key_t, value_t> map_t; //!< Diff entry key-value map.
  
  //! Constructor.
  diff_entry(const string& type);
  diff_entry();
  
  //! Copy constructor
  //diff_entry(const diff_entry& rhs) : _fields(rhs._fields), _type(rhs._type), _id(rhs._id), _parents(rhs._parents) {}

  //! Destructor.
  virtual ~diff_entry() { }
  
  //! Ease-of-use accessor for setting fields on this entry.
  value_t& operator[](const key_t& k) { return _fields[k]; }
  bool entry_exists(const key_t& k) { return (_fields.count(k) > 0); }

  //! Marshal this diff entry into an ordered list of fields.
  virtual void marshal(field_list_t& s);

  diff_entry& operator()(const key_t& key, const value_t& value) {
  	(*this)[key] = value;
  	return *this;
  }

  //! Clone this entry.
  //virtual diff_entry* clone() const = 0;
  
  //! Parameters most diff_entrys have in common
  string _type;
  string _id;
  vector<string> _evidence; 
  map_t _fields; //!< Additional information about this diff entry. Look at 
};

void add_reject_reason(diff_entry& de, const string &reason);
uint32_t number_reject_reasons(diff_entry& de);

//! Convert genome diff to GVF
void GDtoGVF( const string& gdfile, const string& gfffile );

//! Convert VCF to genome diff
void VCFtoGD( const string& vcffile, const string& gfffile );

//! Output operator for a diff entry.
ostream& operator<<(ostream& out, diff_entry& de);

//! Genome Diff Sorting
//! For sorting by a number, then by fields to break ties
struct sort_fields_item {

  //! Constructor.
  sort_fields_item() {_f1=0; _f2=""; _f3=""; };

  sort_fields_item(uint8_t f1, string f2, string f3) :
    _f1(f1), _f2(f2), _f3(f3) {};
  
  //! Destructor.
  virtual ~sort_fields_item() {};

  uint8_t _f1;
  string _f2;
  string _f3;
};

//! Sort routine
bool diff_entry_sort(const counted_ptr<diff_entry>& a, const counted_ptr<diff_entry>& b);





/*! Genome diff class.
 
 //  		Genome Diff files are tab delimitted. The first column defines the type of entry on a given line.
 //  		The second and third columns are type-nonspecific (id, parents), followed by type-specific
 //  		columns, then an arbitrary number of columns of more detailed data in a key=value format.
 
 */
class genome_diff {
public:

  typedef vector<counted_ptr<diff_entry> > entry_list_t; //!< Type for a list of diff entries.a
  typedef vector<counted_ptr<diff_entry> > entry_vector_t;
  typedef counted_ptr<diff_entry> diff_entry_ptr;

  typedef string key_t; 
  typedef vector<string> list_t;
  
  //! Constructor.
  genome_diff() : _unique_id_counter(0) { }
  
  //! Constructor that sets a default filename.
  genome_diff(const string& filename);

  //! Constructor that replaces ::merge(1,2) function
  genome_diff(genome_diff& merge1, genome_diff& merge2);

  //! Destructor.
  ~genome_diff() { }

  //! Retrieve a new diff entry id for this genome diff.
  uint32_t new_unique_id();
  
  //! Add evidence to this genome diff.
  void add(const diff_entry& item);

  //! Read a genome diff from the default file to build entry_list_t,
  void read() { read(_default_filename); }

  //! Write the genome diff to the default file.
  void write() { write(_default_filename); }

  //! Read a genome diff from a file.
  //TEST 6/12 Completed @GRC
  void read(const string& filename);
  
  //! Write the genome diff to a file.
  void write(const string& filename);
  
  //! Remove items used as evidence by any mutations out of input list
  ::list<counted_ptr<diff_entry> > filter_used_as_evidence(const entry_list_t& list);
  
  //! Retrieve diff_entrys that match given type(s) 
  entry_list_t list(vector<string> types);
  
  //! Converts a genome_diff(.gd) file's line to a diff_entry
  diff_entry _line_to_item(const string& line);

  //Returns _entry_list with matching item._evidence
  entry_list_t mutation_evidence_list(const diff_entry& item);

  entry_list_t mutation_list();

  diff_entry_ptr parent(const diff_entry& item);

  void strcopy(char* arg1, const char* arg2);

  //Metadata kept in .gd files
  struct Metadata
  {
  string run_id;
  string version; 
  string author;
  string ref_seq;
  vector<string> read_seq;
  };

  Metadata metadata;
  
  

protected:  	
  const string _default_filename; //!< Default filename for this diff.
  entry_list_t _entry_list; //!< All diff entries.
  uint32_t _unique_id_counter; //!< Smallest available id.
  map<uint32_t,bool> unique_id_used;
  string version; //!< .gd file #=GENOME_DIFF version, found on first line
};


}

#endif
