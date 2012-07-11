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

#ifndef _BRESEQ_GENOME_DIFF_H_
#define _BRESEQ_GENOME_DIFF_H_

#include "common.h"

namespace breseq {

class cAnnotatedSequence;
class cDiffEntry;
class cGenomeDiff;
class cReferenceSequences;

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
extern const char* SIDE_1_SEQ_ID;
extern const char* SIDE_1_POSITION;
extern const char* SIDE_1_STRAND;
extern const char* SIDE_1_REDUNDANT;

extern const char* SIDE_2_SEQ_ID;
extern const char* SIDE_2_POSITION;
extern const char* SIDE_2_STRAND;
extern const char* SIDE_2_REDUNDANT;
  
extern const char* SIDE_1_READ_COUNT;
extern const char* SIDE_2_READ_COUNT;
extern const char* NEW_JUNCTION_READ_COUNT;
  
extern const char* SIDE_1_COVERAGE;
extern const char* SIDE_2_COVERAGE;
extern const char* NEW_JUNCTION_COVERAGE;

// Types of diff entries:
enum gd_entry_type {UNKNOWN = 0, SNP, SUB, DEL, INS, MOB, AMP, INV, CON, RA,
                    MC, JC, CN, UN, CURA, FPOS, PHYL, TSEQ, PFLP, RFLP, PFGE, NOTE};

extern const vector<string> gd_entry_type_lookup_table;
  
inline string to_string(const gd_entry_type type)
{
  return gd_entry_type_lookup_table[type];
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
 
 For convenience, there are also factory methods on cGenomeDiff that create diff
 entries that have been pre-populated according to their type.
 */
  
typedef string diff_entry_key_t; //!< Diff entry keys.
typedef string diff_entry_value_t; //!< Diff entry values.
typedef map<diff_entry_key_t, diff_entry_value_t> diff_entry_map_t; //!< Diff entry key-value map.

typedef counted_ptr<cDiffEntry> diff_entry_ptr_t;

class cDiffEntry : public diff_entry_map_t {
public: 

  //! Constructor.
  cDiffEntry(const gd_entry_type type);
  cDiffEntry(const string &line); //For deserialization from gd file.
  cDiffEntry(diff_entry_map_t& de) : diff_entry_map_t(de) {};
  cDiffEntry();
  
  //! Copy constructor
  //cDiffEntry(const cDiffEntry& rhs) : _fields(rhs._fields), _type(rhs._type), _id(rhs._id), _parents(rhs._parents) {}

  //! Destructor.
  virtual ~cDiffEntry() { }

  //! Return if a given key value exists in _fields
  bool entry_exists(const diff_entry_key_t& k) const { return (count(k) > 0); }
  
  //! Const accessor
  diff_entry_value_t get(const diff_entry_key_t key) const
  {
    return this->find(key)->second;
  }
    
  //! Return if this diff entry is a mutation
  bool is_mutation() const;

  //! Return if this diff entry is evidence
  bool is_evidence() const;

  //! Return if this diff entry is a validation
  bool is_validation() const;
  
  bool is_marked_deleted()
  { return (this->entry_exists("deleted") && ((from_string<int32_t>((*this)["deleted"]) != 0))); }
  
  //! Common function for getting start or end of mutation or evidence
  uint32_t get_start();
  uint32_t get_end();

  //! Marshal this diff entry into an ordered list of fields.
  virtual void marshal(field_list_t& s) const;

  //! Serialize this diff entry into a string for output.
  virtual string to_string(void) const;

  //! Returns values for cDiffEntry["reject"]
  vector<string> get_reject_reasons();

  //const cDiffEntry& operator= (diff_entry_ptr_t ptr) const { return *ptr; }
  
  cDiffEntry to_spec(void) const;


  void normalize_to_sequence(const cAnnotatedSequence &seq);

  size_t number_reject_reasons();
  int32_t mutation_size_change(cReferenceSequences& ref_seq_info);

  struct by_scores;
  struct is_type;
  struct is_not_type;
  struct field_exists;
  struct fields_exist;
  struct no_show;
  struct rejected;

  cDiffEntry& operator()(const diff_entry_key_t& key, const diff_entry_value_t& value) {
    (*this)[key] = value;
    return *this;
  }

  bool operator== (const cDiffEntry& de);

  //! Clone this entry.
  //virtual cDiffEntry* clone() const = 0;
  
  //! Parameters most cDiffEntrys have in common
  gd_entry_type _type;
  string _id;
  vector<string> _evidence; 
  
};
typedef list<diff_entry_ptr_t> diff_entry_list_t; //!< Type for a list of diff entries.

void add_reject_reason(cDiffEntry& de, const string &reason);

//! Convert genome diff to GVF
void GDtoGVF( const string& gdfile, const string& gvffile, bool snv_only = false );

//! Convert VCF to genome diff
void VCFtoGD( const string& vcffile, const string& gfffile );

//! Convert GD to Circos files
void GDtoCircos(const vector<string> &gd_file_names,
                const vector<string> &reference_file_names,
                const string &circos_directory,
                double distance_scale,
                double feature_scale);

//! Convert MIRA feature analysis file to GD
void MIRAtoGD(const string &mira_file_name, const string &gd_file_name);
  
//! Output operator for a diff entry.
ostream& operator<<(ostream& out, const cDiffEntry& de);

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
  



//! Sort routines
bool diff_entry_ptr_sort(const diff_entry_ptr_t& a, const diff_entry_ptr_t& b);
bool diff_entry_sort(const cDiffEntry &a, const cDiffEntry &b);

/*! Genome diff class.
 
 //  		Genome Diff files are tab delimitted. The first column defines the type of entry on a given line.
 //  		The second and third columns are type-nonspecific (id, parents), followed by type-specific
 //  		columns, then an arbitrary number of columns of more detailed data in a key=value format.
 
 */
class cGenomeDiff
{
public:

  enum group { MUTATIONS = 0, EVIDENCE, VALIDATION }; 

  typedef string key_t; 
  typedef vector<string> list_t;
  
  //! Constructor.
  cGenomeDiff() : _unique_id_counter(0) { }

  //! Constructor from file
  cGenomeDiff(const string& filename);

  //! Constructor that replaces ::merge(1,2) function
  cGenomeDiff(cGenomeDiff& merge1, cGenomeDiff& merge2, bool unique=true, bool new_id=true, bool verbose=false);

  //! Destructor.
  ~cGenomeDiff() { }

  //! Retrieve a new diff entry id for this genome diff.
  uint32_t new_unique_id();
  
  //! Add evidence to this genome diff, return assigned ID.
  uint32_t add(const cDiffEntry& item, bool lowest_unique=false);
  
  //! Subtract mutations using gd_ref as reference.
  void set_subtract(cGenomeDiff& gd_ref, bool verbose=false);

  void set_intersect(cGenomeDiff& gd_ref, bool verbose=false);
  
  void set_union(cGenomeDiff& gd_ref, bool verbose=false);
  
  //! Merge GenomeDiff information using gd_new as potential new info.
  void merge(cGenomeDiff& gd_new, bool unique=true, bool new_id=false, bool verbose=false);

  //! fast merge, doesn't compare entries, but does renumber
  void fast_merge(const cGenomeDiff& gd);
  
  //! sort
  void sort() { _entry_list.sort(diff_entry_ptr_sort); }
  void unique();

  void compare(cGenomeDiff& gd, bool verbose);

  void assign_unique_ids(void);


  static cGenomeDiff from_vcf(const string &file_name);

  //! Read a genome diff from a file.
  void read(const string& filename);
  
  //! Write the genome diff to a file.
  void write(const string& filename);

  //! Remove mutations, evidence, validation.
  void remove(cGenomeDiff::group group);

  //! Removes all GD entries that aren't used as evidence.
  void filter_not_used_as_evidence(bool verbose=false);

  //! Remove items used as evidence by any mutations out of input list
  diff_entry_list_t filter_used_as_evidence(const diff_entry_list_t& list);
  
  //! Call to check if loaded info seq_ids match supplied reference.
  bool is_valid(cReferenceSequences& ref_seq_info, bool verbose=false);
  
  //! Call to generate random mutations.
  void random_mutations(const string& exclusion_file, const string& type, uint32_t number, uint32_t read_length, cAnnotatedSequence& ref_seq_info, uint32_t rand_seed, bool verbose=false);
  
  //! Retrieve cDiffEntrys that match given type(s) 
  const diff_entry_list_t list() const { return _entry_list; }
  diff_entry_list_t list(const vector<gd_entry_type>& types = vector<gd_entry_type>());
  
  //! retrieve cDiffEntrys that match given type(s) and do not have 'no_show' set
  diff_entry_list_t show_list(const vector<gd_entry_type>& types = vector<gd_entry_type>());
  
  //! Returns _entry_list with matching item._evidence
  diff_entry_list_t mutation_evidence_list(const cDiffEntry& item);

  diff_entry_list_t mutation_list();
  diff_entry_list_t evidence_list();
  diff_entry_list_t validation_list();

  diff_entry_ptr_t parent(const cDiffEntry& evidence);

  void normalize_to_sequence(cReferenceSequences &ref_seq);
  void mutations_to_evidence(cReferenceSequences &ref_seq);
    
  //Additional functions that need? adding from GenomeDiff.gm
  void add_reject_reasons(cDiffEntry item, const string& reason);
  size_t number_reject_reasons(cDiffEntry item);
  
  bool mutation_in_entry_of_type(cDiffEntry mut, const gd_entry_type type);
  bool mutation_unknown(cDiffEntry mut) { return mutation_in_entry_of_type(mut, UN); }
  bool mutation_deleted(cDiffEntry mut) { return mutation_in_entry_of_type(mut, DEL); }

  void apply_to_sequences(cReferenceSequences &ref_seq_info, cReferenceSequences& new_ref_seq_info, bool verbose=false);
  void shift_positions(cDiffEntry& item, cReferenceSequences& ref_seq_info, bool verbose=false);

  void strcopy(char* arg1, const char* arg2);

  void add_breseq_data(const key_t& key, const string& value);

  //! Metadata kept in .gd files
  struct Metadata
  {
    Metadata() : version("1.0") {}

    string run_name;
    string version; 
    string author;
    vector<string> ref_seqs;
    vector<string> read_seqs;
    map<string,string> breseq_data; // Use this to write values from pipeline to gd
  } metadata;


protected:  	
  string _default_filename; //!< Default filename for this diff.
  diff_entry_list_t _entry_list; //!< All diff entries.
  uint32_t _unique_id_counter; //!< Smallest available id.
  map<uint32_t,bool> unique_id_used;
};

//! Functor. Wraps cDiffEntry.entry_exists() for use in STL algorithms.
//Returns true if a cDiffEntry contains the given field_key.
struct cDiffEntry::field_exists : public unary_function <diff_entry_ptr_t, bool>
{
  //! Constructor
  explicit field_exists (const diff_entry_key_t& field_key)
    : m_field_key(field_key) {}
  
  //! Predicate
  virtual bool operator() (const diff_entry_ptr_t& p_diff_entry)
    const {return (p_diff_entry->entry_exists(m_field_key));}

  protected:
    diff_entry_key_t m_field_key;
};

//! Functor. Wraps cDiffEntry.entry_exists() for use in STL algorithms.
//Returns true if a cDiffEntry contains all of the given field_keys.
//ie:  cDiffEntry[field_key_1] && cDiffEntry[field_key_2]
struct cDiffEntry::fields_exist : public unary_function <diff_entry_ptr_t, bool>
{
  //! Constructor
  explicit fields_exist (const vector<diff_entry_key_t>& field_keys)
    : m_field_keys(field_keys) {}
  
  //! Predicate
  virtual bool operator() (const diff_entry_ptr_t& p_diff_entry) const
  {
    for (vector<diff_entry_key_t>::const_iterator itr = m_field_keys.begin();
         itr != m_field_keys.end(); itr++) {
      diff_entry_key_t field_key(*itr);
      if (p_diff_entry->entry_exists(field_key)) {
        continue;
      }
      else { 
        return false;
      }
    }
    
    // cDiffEntry contains all field_keys
    return true;
  }
  protected:
    vector<diff_entry_key_t> m_field_keys;
};

//! Functor. Sorts cDiffEntrys in decending order depending on given fields that
//can be evaluated as an unsigned integer.
struct cDiffEntry::by_scores : public binary_function
  <diff_entry_ptr_t, diff_entry_ptr_t, bool>
{

  //! Constructor
  explicit by_scores (const vector<diff_entry_key_t>& field_keys)
    : m_field_keys(field_keys) {}
  
  //! Predicate
  virtual bool operator() (const diff_entry_ptr_t& a, const diff_entry_ptr_t& b) const
  {
    for (vector<diff_entry_key_t>::const_iterator itr = m_field_keys.begin(); itr != m_field_keys.end(); itr++) {
      string key(*itr);

      if (from_string<uint32_t>((*a)[key]) == from_string<uint32_t>((*b)[key]))
        continue;
      else 
        return from_string<uint32_t>((*a)[key]) > from_string<uint32_t>((*b)[key]);
    }
    return false;
  }
  
  protected:
    vector<diff_entry_key_t> m_field_keys;
};


struct cDiffEntry::is_type: unary_function <diff_entry_ptr_t, bool>
{
  //! Constructor
  explicit is_type(const gd_entry_type type)
    : m_type(type) {}

  //! Predicate 
  virtual bool operator() (const diff_entry_ptr_t& cDiffEntry)
    const {return (*cDiffEntry)._type == m_type;}


  protected:
    gd_entry_type m_type;
};
  
struct cDiffEntry::is_not_type: unary_function <diff_entry_ptr_t, bool>
{
  //! Constructor
  explicit is_not_type(const gd_entry_type type)
  : m_type(type) {}
  
  //! Predicate 
  virtual bool operator() (const diff_entry_ptr_t& cDiffEntry)
  const {return (*cDiffEntry)._type != m_type;}
  
  
protected:
  gd_entry_type m_type;
};

struct cDiffEntry::no_show:public unary_function<diff_entry_ptr_t, bool>
{
  virtual bool operator() (const diff_entry_ptr_t& cDiffEntry) const
  {
    return cDiffEntry->entry_exists("no_show");
  }
};


struct cDiffEntry::rejected:public unary_function<diff_entry_ptr_t,bool>
{
  virtual bool operator() (diff_entry_ptr_t cDiffEntry)
  {
    return cDiffEntry->entry_exists(REJECT);
  }

};



}
#endif
