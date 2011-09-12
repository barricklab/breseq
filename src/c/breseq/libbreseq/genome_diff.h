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
enum Type {TYPE_UNKOWN, SNP, SUB, DEL, INS, MOB, AMP, INV, CON, NOT_MUTATION, RA, MC, JC, UN};
inline string to_string(const Type type)
{
  switch(type) {
    case SNP: return "SNP";
    case SUB: return "SUB";
    case DEL: return "DEL";
    case INS: return "INS";
    case MOB: return "MOB";
    case AMP: return "AMP";
    case INV: return "INV";
    case CON: return "CON";
    case RA: return "RA";
    case MC: return "MC";
    case JC: return "JC";
    case UN: return "UN";
    default: return "?";
  }
}

inline Type to_type(const string& type)
{
  if (type == "SNP") return SNP;
  if (type == "SUB") return SUB;
  if (type == "DEL") return DEL;
  if (type == "INS") return INS;
  if (type == "MOB") return MOB;
  if (type == "AMP") return AMP;
  if (type == "INV") return INV;
  if (type == "CON") return CON;
  if (type == "RA") return RA;
  if (type == "MC") return MC;
  if (type == "JC") return JC;
  if (type == "UN") return UN;
return TYPE_UNKOWN;
}


enum Strand {POS_STRAND = 1, NEG_STRAND = -1};
  
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
  diff_entry(const Type type);
  diff_entry();
  
  //! Copy constructor
  //diff_entry(const diff_entry& rhs) : _fields(rhs._fields), _type(rhs._type), _id(rhs._id), _parents(rhs._parents) {}

  //! Destructor.
  virtual ~diff_entry() { }
  
  //! Ease-of-use accessor for setting fields on this entry.
  value_t& operator[](const key_t& k)
  { 
    return _fields[k]; 
  }
  value_t operator[](const key_t& k) const
  { 
    map_t::const_iterator it = _fields.find(k);
    assert(it != _fields.end());
    return it->second; 
  }
  
  bool entry_exists(const key_t& k) const { return (_fields.count(k) > 0); }

  //! Marshal this diff entry into an ordered list of fields.
  virtual void marshal(field_list_t& s);

  //! Returns values for diff_entry["reject"]
  vector<string> get_reject_reasons();
  
  size_t number_reject_reasons();

  struct by_scores;
  struct is_type;
  struct is_not_type;
  struct field_exists;
  struct fields_exist;
  struct frequency_less_than_two_or_no_show;
  struct coverage_or_no_show_is_true;
  struct rejected;


  diff_entry& operator()(const key_t& key, const value_t& value) {
  	(*this)[key] = value;
  	return *this;
  }

  //! Clone this entry.
  //virtual diff_entry* clone() const = 0;
  
  //! Parameters most diff_entrys have in common
  Type _type;
  string _id;
  vector<string> _evidence; 
  map_t _fields; //!< Additional information about this diff entry. Look at 
};

void add_reject_reason(diff_entry& de, const string &reason);

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
  
typedef list<counted_ptr<diff_entry> > diff_entry_list; //!< Type for a list of diff entries.a
typedef counted_ptr<diff_entry> diff_entry_ptr;

//! Sort routine
bool diff_entry_sort(const diff_entry_ptr& a, const diff_entry_ptr& b);


/*! Genome diff class.
 
 //  		Genome Diff files are tab delimitted. The first column defines the type of entry on a given line.
 //  		The second and third columns are type-nonspecific (id, parents), followed by type-specific
 //  		columns, then an arbitrary number of columns of more detailed data in a key=value format.
 
 */
class genome_diff {
public:

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

  //! Read a genome diff from a file.
  void read(const string& filename);
  
  //! Write the genome diff to a file.
  void write(const string& filename);
  
  //! Remove items used as evidence by any mutations out of input list
  diff_entry_list filter_used_as_evidence(const diff_entry_list& list);
  
  //! Retrieve diff_entrys that match given type(s) 
  diff_entry_list list(const vector<Type>& types = vector<Type>());
  //diff_entry_list list() {return _entry_list;}
  
  //! retrieve diff_entrys that match given type(s) and do not have 'no_show' set
  diff_entry_list show_list(const vector<Type>& types = vector<Type>());
  
  //! Converts a genome_diff(.gd) file's line to a diff_entry
  diff_entry _line_to_item(const string& line);

  //! Returns _entry_list with matching item._evidence
  diff_entry_list mutation_evidence_list(const diff_entry& item);

  diff_entry_list mutation_list();
  diff_entry_list evidence_list();

  diff_entry_ptr parent(const diff_entry& item);
    
  //Additional functions that need? adding from GenomeDiff.gm
  void add_reject_reasons(diff_entry item, const string& reason);
  size_t number_reject_reasons(diff_entry item);
  bool mutation_unknown(diff_entry mut);
  bool has_mutation(diff_entry test_item);
  bool interval_un (const uint32_t& start, const uint32_t& end);

  void apply_to_sequences(cReferenceSequences &ref_seq_info);

  void strcopy(char* arg1, const char* arg2);

  //! Metadata kept in .gd files
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
  diff_entry_list _entry_list; //!< All diff entries.
  uint32_t _unique_id_counter; //!< Smallest available id.
  map<uint32_t,bool> unique_id_used;
};

//! Functor. Wraps diff_entry.entry_exists() for use in STL algorithms.
//Returns true if a diff_entry contains the given field_key.
struct diff_entry::field_exists : public unary_function <diff_entry_ptr, bool>
{
  //! Constructor
  explicit field_exists (const diff_entry::key_t& field_key)
    : m_field_key(field_key) {}
  
  //! Predicate
  virtual bool operator() (const diff_entry_ptr& p_diff_entry) 
    const {return (p_diff_entry->entry_exists(m_field_key));}

  protected:
    diff_entry::key_t m_field_key;
};

//! Functor. Wraps diff_entry.entry_exists() for use in STL algorithms.
//Returns true if a diff_entry contains all of the given field_keys.
//ie:  diff_entry[field_key_1] && diff_entry[field_key_2]
struct diff_entry::fields_exist : public unary_function <diff_entry_ptr, bool> 
{
  //! Constructor
  explicit fields_exist (const vector<diff_entry::key_t>& field_keys)
    : m_field_keys(field_keys) {}
  
  //! Predicate
  virtual bool operator() (const diff_entry_ptr& p_diff_entry) const
  {
    for (vector<diff_entry::key_t>::const_iterator itr = m_field_keys.begin();
         itr != m_field_keys.end(); itr++) {
      genome_diff::key_t field_key(*itr);
      if (p_diff_entry->entry_exists(field_key)) 
        continue;
      else 
        return false;
    }
    
    // diff_entry contains all field_keys
    return true;
  }
  protected:
    vector<diff_entry::key_t> m_field_keys;
};

//! Functor. Sorts diff_entrys in decending order depending on given fields that
//can be evaluated as an unsigned integer.
struct diff_entry::by_scores : public binary_function
  <diff_entry_ptr, diff_entry_ptr, bool>
{

  //! Constructor
  explicit by_scores (const vector<key_t>& field_keys)
    : m_field_keys(field_keys) {}
  
  //! Predicate
  virtual bool operator() (const diff_entry_ptr& a, const diff_entry_ptr& b) const 
  {
    for (vector<key_t>::const_iterator itr = m_field_keys.begin();
         itr != m_field_keys.end(); itr++) {
      string key(*itr);

      if (from_string<uint32_t>((*a)[key]) == from_string<uint32_t>((*b)[key]))
        continue;
      else 
        return from_string<uint32_t>((*a)[key]) > from_string<uint32_t>((*b)[key]);
    }
    return false;
  }
  
  protected:
    vector<key_t> m_field_keys;
};


struct diff_entry::is_type: unary_function <diff_entry_ptr, bool>
{
  //! Constructor
  explicit is_type(const string& type)
    : m_type(type) {}

  //! Predicate 
  virtual bool operator() (const diff_entry_ptr& diff_entry)
    const {return (*diff_entry)._type == to_type(m_type);}


  protected:
    string m_type;
};
  
struct diff_entry::is_not_type: unary_function <diff_entry_ptr, bool>
{
  //! Constructor
  explicit is_not_type(const string& type)
  : m_type(type) {}
  
  //! Predicate 
  virtual bool operator() (const diff_entry_ptr& diff_entry)
  const {return (*diff_entry)._type != to_type(m_type);}
  
  
protected:
  string m_type;
};

struct diff_entry::frequency_less_than_two_or_no_show:public unary_function<diff_entry_ptr, bool>
{
  virtual bool operator() (const diff_entry_ptr& diff_entry) const
  {
    return ((*diff_entry)["frequency"] == "0" || (*diff_entry)["frequency"] == "1" || (*diff_entry).entry_exists("no_show"));
  }
};

struct diff_entry::coverage_or_no_show_is_true:public unary_function<diff_entry_ptr,bool>
{
  virtual bool operator() (const diff_entry_ptr& diff_entry) const
  {
    if ( (*diff_entry).entry_exists("no_show") && from_string<bool>((*diff_entry)["no_show"]) )
    {
      return true;
    }
    
    return false;
  }
};

struct diff_entry::rejected:public unary_function<diff_entry_ptr,bool>
{
  virtual bool operator() (diff_entry_ptr diff_entry) 
  {
    return diff_entry->entry_exists(REJECT);
  }

};

}
#endif
