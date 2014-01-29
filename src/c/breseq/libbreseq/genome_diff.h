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
class cSequenceFeature;
  
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
extern const char* GENOTYPE_QUALITY;
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
extern const char* NEW_JUNCTION_FREQUENCY;
  
extern const char* SIDE_1_COVERAGE;
extern const char* SIDE_2_COVERAGE;
extern const char* NEW_JUNCTION_COVERAGE;

// Types of diff entries:
enum gd_entry_type {UNKNOWN = 0, SNP, SUB, DEL, INS, MOB, AMP, INV, CON, RA,
                    MC, JC, CN, UN, CURA, FPOS, PHYL, TSEQ, PFLP, RFLP, PFGE, NOTE, MASK};

extern const vector<string> gd_entry_type_lookup_table;
  
inline string to_string(const gd_entry_type type)
{ return gd_entry_type_lookup_table[type]; }

  
/*! Parse errors
 
 Helper class for aggregating all of the errors encountered when reading a Genome Diff file
 
 */

class cFileParseErrors {
public:
  
  class sFileParseError {
  public:
    uint32_t _line_number;
    string _line;
    string _error_message;
    
    sFileParseError(uint32_t line_number, const string& line, const string& error_message )
    :_line_number(line_number), _line(line), _error_message(error_message)
    {}
    
    bool operator <(const sFileParseError& compare) const
    { return this->_line_number < compare._line_number; }
    
    void print()
    {
      string tab_line = substitute(_line, "\t", "<tab>");
      
      cerr << ">>ERROR: " << _error_message << endl;
      if ((_line_number != 0) || (_line != "")) {
        cerr << ">>LINE: " << setw(5) <<  _line_number << endl << tab_line << endl;
      }
    }
  };
  
  // List pairs of line number and error message
  list<sFileParseError> _errors;
  string _filename;
  bool _fatal;
  
  cFileParseErrors(const string& filename = "")
  : _filename(filename), _fatal(false)
  { }
  
  void add_line_error(const uint32_t line_number, const string& line, const string& message, bool fatal)
  {
    _fatal = _fatal || fatal;
    _errors.push_back( sFileParseError(line_number, line, message) );
  }
  
  void print_errors(bool print_file = true)
  {
    if (_errors.size() == 0) return;
    
    if (print_file) {
      cerr << endl;
      cerr << ">>> Error(s) in GenomeDiff format. FILE: " << _filename << " <<<" << endl;
    }
    
    _errors.sort();
    
    for(list<sFileParseError>::iterator it = _errors.begin(); it != _errors.end(); it++) {
      it->print();
    }
    cerr << endl;
  }
  
  bool fatal() {return _fatal;}
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
typedef list<diff_entry_ptr_t> diff_entry_list_t; //!< Type for a list of diff entries.

class cDiffEntry : public diff_entry_map_t {
public: 

  //!---- Variables ---- !//
  // @JEB: these should be made protected and given accessors
  gd_entry_type _type;
  string _id;
  vector<string> _evidence; 
  
  //!---- Constructor / Destructor ---- !//

  //! Constructor.
  cDiffEntry();
  cDiffEntry(const gd_entry_type type);
  cDiffEntry(const string &line, uint32_t line_number, cFileParseErrors* file_parse_errors = NULL); //For deserialization from gd file.
  cDiffEntry(diff_entry_map_t& de) : diff_entry_map_t(de) {};
  
  //! Helper function
  static gd_entry_type type_to_enum(string type);
  
  //! Checks fields for expected types (such as integers)
  void valid_field_variable_types(cFileParseErrors& parse_errors);

  //! Copy constructor
  //cDiffEntry(const cDiffEntry& rhs) : _fields(rhs._fields), _type(rhs._type), _id(rhs._id), _parents(rhs._parents) {}

  //! Destructor.
  virtual ~cDiffEntry() { }

  
  //!---- Accessors to generic properties ---- !//

  //! Accessor convenience function
  cDiffEntry& operator()(const diff_entry_key_t& key, const diff_entry_value_t& value) {
    (*this)[key] = value;
    return *this;
  }
  
  //! Const accessor
  diff_entry_value_t get(const diff_entry_key_t key) const { return this->find(key)->second; }
    
  //! Comparison operator
  bool operator== (const cDiffEntry& de);
  
  //! Return if a given key value exists in _fields
  bool entry_exists(const diff_entry_key_t& k) const { return (count(k) > 0); }
  
  //! Return if this diff entry is a mutation
  bool is_mutation() const;

  //! Return if this diff entry is evidence
  bool is_evidence() const;

  //! Return if this diff entry is a validation
  bool is_validation() const;
  
  bool is_marked_deleted()
  { return (this->entry_exists("deleted") && ((from_string<int32_t>((*this)["deleted"]) != 0))); }
  
  //!---- Accessors to calculated properties ---- !//
  
  //! Common function for getting start or end of mutation or evidence
  uint32_t get_start();
  uint32_t get_end();

  //! Common function giving change in size of genome at site of applying entry
  int32_t mutation_size_change(cReferenceSequences& ref_seq_info);
  
  //!---- Output ---- !//
  
  //! Marshal this diff entry into an ordered list of fields.
  virtual void marshal(vector<string>& s) const;
  
  //! Serialize this diff entry into a string for output.
  virtual string as_string(void) const;
  
  //!---- Reject Reasons Field ---- !//
  
  size_t number_reject_reasons();
  
  //! Returns values for cDiffEntry["reject"]
  vector<string> get_reject_reasons();
  
  //! Adds a reject reason to cDiffEntry["reject"] as a list
  void add_reject_reason(const string &reason);

  //!---- Simplifying entries ---- !//
  
  //! Remove all information except required fields
  cDiffEntry to_spec(void) const;

  void normalize_to_sequence(const cAnnotatedSequence &seq, bool verbose = false);

  //!---- Sorting ---- !//
  
  //! Various functors for testing many entries at once for a property

  
  //! Functor. Sorts cDiffEntrys in decending order depending on given fields that
  //can be evaluated as an unsigned integer.
  struct by_scores : public binary_function
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
        
        if (from_string<double>((*a)[key]) == from_string<double>((*b)[key]))
          continue;
        else 
          return from_string<double>((*a)[key]) > from_string<double>((*b)[key]);
      }
      return false;
    }
    
  protected:
    vector<diff_entry_key_t> m_field_keys;
  };
  
  
  struct is_type: unary_function <diff_entry_ptr_t, bool>
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
  
  struct is_not_type: unary_function <diff_entry_ptr_t, bool>
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
  
  struct no_show:public unary_function<diff_entry_ptr_t, bool>
  {
    virtual bool operator() (const diff_entry_ptr_t& cDiffEntry) const
    {
      return cDiffEntry->entry_exists("no_show");
    }
  };
  
  
  struct rejected:public unary_function<diff_entry_ptr_t,bool>
  {
    virtual bool operator() (diff_entry_ptr_t cDiffEntry)
    {
      return cDiffEntry->entry_exists(REJECT);
    }
    
  };
  
  struct rejected_and_not_user_defined:public unary_function<diff_entry_ptr_t,bool> {
    virtual bool operator() (diff_entry_ptr_t cDiffEntry)
    {
      return cDiffEntry->entry_exists(REJECT) && !cDiffEntry->entry_exists("user_defined");
    }
    
  };
  
  //! Functor. Wraps cDiffEntry.entry_exists() for use in STL algorithms.
  //Returns true if a cDiffEntry contains the given field_key.
  struct field_exists : public unary_function <diff_entry_ptr_t, bool>
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
  struct fields_exist : public unary_function <diff_entry_ptr_t, bool>
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
  

};

// Overload to output Genome DIff entries
inline ostream &operator<<( ostream &out, const cDiffEntry &de ) {
  out << de.as_string();
  return out;
}
  

/*! Genome Diff class.
 
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

  //!---- Variables ---- !//
protected:  	
  string _filename;                   //!< File name associated with this diff.
  diff_entry_list_t _entry_list;      //!< All diff entries.
  uint32_t _unique_id_counter;        //!< Smallest available id.
  map<uint32_t,bool> unique_id_used;

public:
  // @JEB should make this protected and add an accessor
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
  
  
  //! ---- Constructor / Destructor ---- !//

  //! Constructor.
  cGenomeDiff() : _unique_id_counter(0) { }

  //! Constructor from file
  cGenomeDiff(const string& filename);

  //! Constructor that replaces ::merge(1,2) function
  cGenomeDiff(cGenomeDiff& merge1, cGenomeDiff& merge2, bool unique=true, bool new_id=true, bool verbose=false);

  //! Destructor.
  ~cGenomeDiff() { }
  
  
  //!---- Accessors ---- !//

  string file_name() const {return _filename;}

  void add_breseq_data(const key_t &key, const string& value)
    { this->metadata.breseq_data.insert(pair<string,string>(key, value)); }

  //! Gets parent of entry, if there is one
  diff_entry_ptr_t parent(const cDiffEntry& evidence);
  
  //!---- Input and Output ---- !//
  
  //! Read a genome diff from a file.
  cFileParseErrors read(const string& filename, bool suppress_errors = false);
  
  //! Check to see if genome diff is valid with reference sequences
  cFileParseErrors valid_with_reference_sequences(cReferenceSequences& ref_seq, bool suppress_errors = false);
    
  //! Write the genome diff to a file.
  void write(const string& filename);
  
  //!---- Adding and Removing Entries ---- !//
  
  //! Helper function to find next unused id
  uint32_t new_unique_id();
  
  //! Add an item to this genome diff.
  diff_entry_ptr_t add(const cDiffEntry& item, bool lowest_unique=false);
  
  //! Remove mutations, evidence, or validation.
  void remove(cGenomeDiff::group group);
  
  //!---- Accessing Entries ---- !//
  
  diff_entry_ptr_t find_by_id(string _id);
  
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
  
  //! Removes all GD entries that aren't used as evidence.
  void filter_not_used_as_evidence(bool verbose=false);
  
  //! Remove items used as evidence by any mutations out of input list
  diff_entry_list_t filter_used_as_evidence(const diff_entry_list_t& list);
  
  //! Helper function for returning subsets below
  bool mutation_in_entry_of_type(cDiffEntry mut, const gd_entry_type type);
  bool mutation_unknown(cDiffEntry mut) { return mutation_in_entry_of_type(mut, UN); }
  bool mutation_deleted(cDiffEntry mut) { return mutation_in_entry_of_type(mut, DEL); }
  
  //!---- Set Operations ---- !//
  
  //! Subtract mutations using gd_ref as reference.
  void set_subtract(cGenomeDiff& gd_ref, bool verbose=false);

  void set_intersect(cGenomeDiff& gd_ref, bool verbose=false);
  
  void set_union(cGenomeDiff& gd_ref, bool verbose=false);
  
  //! Helper function for union
  void unique();
  
  //! Merge GenomeDiff information using gd_new as potential new info.
  void merge(cGenomeDiff& gd_new, bool unique=true, bool new_id=false, bool verbose=false);
  
  //! fast merge, doesn't compare entries, but does renumber
  void fast_merge(const cGenomeDiff& gd);
  
  //! Helper function for fixing IDs after a set operation
  void reassign_unique_ids();
  
  //!---- Sorting Items in Genome Diff ---- !//
  
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
  
  static bool diff_entry_ptr_sort(const diff_entry_ptr_t& a, const diff_entry_ptr_t& b);
  void sort() { _entry_list.sort(diff_entry_ptr_sort); }
  
  
  //!---- Simulating and Applying Mutations ---- !//
  
  //! Call to generate random mutations.
  void random_mutations(string exclusion_file,
                        string type,
                        uint32_t n_muts,
                        uint32_t buffer,
                        cAnnotatedSequence& ref,
                        bool verbose = false);

  void mutations_to_evidence(cReferenceSequences &ref_seq, bool remove_mutations = true);
  
  // Helper function for apply_to_sequences
  void shift_positions(cDiffEntry& item, cReferenceSequences& ref_seq_info, bool verbose=false);

  // For constructing the sequence a MOB replaces things with
  string mob_replace_sequence(cReferenceSequences& ref_seq_info, 
                              cDiffEntry& mut, 
                              string* picked_seq_id = NULL, 
                              cSequenceFeature* picked_sequence_feature = NULL
                              );
  
  //! Call to apply Genome Diff to sequences
  void apply_to_sequences(cReferenceSequences &ref_seq_info, cReferenceSequences& new_ref_seq_info, bool verbose=false);
  
  //! Shift mutations to preferred descriptions
  void normalize_to_sequence(cReferenceSequences &ref_seq);
   
  //!---- Comparing known lists of mutations/evidence to test files ---- !//
  
  static cGenomeDiff check(cGenomeDiff& ctrl, cGenomeDiff& test, bool verbose = false);
  
  static cGenomeDiff check_evidence(cReferenceSequences& sequence,
                                       uint32_t buffer,
                                       uint32_t shorten_length,
                                       cGenomeDiff& ctrl,
                                       cGenomeDiff& test,
                                       bool jc_only_accepted,
                                       bool verbose = false);
  static void write_jc_score_table(cGenomeDiff& compare, string table_file_path, bool verbose = false); 

  //!---- Format Conversion Functions: Member ---- !//

  // ! VCF files
  void read_vcf(const string& filename);
  void write_vcf(const string& filename, cReferenceSequences& ref_seq_info);


  //! GVF files
  void write_gvf(const string& filename, cReferenceSequences& ref_seq_info, bool snv_only = false);
  
  //!---- Format Conversion Functions: Static Convenience ---- !//

  //! Convert genome diff to GVF
  static void GD2GVF( const string& gdfile, const string& gvffile, cReferenceSequences& ref_seq_info, bool snv_only = false )
    { cGenomeDiff gd(gdfile); gd.write_gvf(gvffile, ref_seq_info, snv_only); }
  
  //! Convert VCF to genome diff
  static void GD2VCF( const string &gdfile, const string & vcffile, cReferenceSequences& ref_seq_info)
    { cGenomeDiff gd(gdfile); gd.write_vcf(vcffile, ref_seq_info); }

  static void VCF2GD( const string& vcffile, const string& gdfile )
  { cGenomeDiff gd; gd.read_vcf(vcffile); gd.write(gdfile); }
  
  //! Convert GD to Circos files
  static void GD2Circos(const vector<string> &gd_file_names,
                 const vector<string> &reference_file_names,
                 const string &circos_directory,
                 double distance_scale,
                 double feature_scale);
  
  //! Convert MIRA feature analysis file to GD @JEB deprecated 08-16-2013
  //static void MIRA2GD(const string &mira_file_name, const string &gd_file_name);

};
  
  
/* Helper class for cGenomeDiff::random_mutations, handles sorting of start_1,
 * end_1 positions and merges the positions if they overlap.
 */
class cFlaggedRegions  {
public:
  typedef pair<uint32_t, uint32_t>  region_t;
  typedef set<region_t>             regions_t;
  
  cFlaggedRegions()
  :m_regions() {
    return;
  }
  
  //! I/O.
  cFlaggedRegions& read_mummer(string file_path, cAnnotatedSequence& ref_seq);
  cFlaggedRegions& read_nucmer_tab_coords(string file_path);
  void write(string file_path);
  void print(void);
  
  
  //! Add region to be marked.
  cFlaggedRegions& flag_region(uint32_t start_1, uint32_t end_1 = 0);
  
  //! Remove overlapping regions, adds segments if partial overlapping occurs.
  //cFlaggedRegions& unflag_region(uint32_t start_1, uint32_t end_1 = 0);
  
  //! Tests if start_1 to end_1 spans over a marked region.
  bool is_flagged(uint32_t start_1, uint32_t end_1 = 0);
  
  //! Test if two regions overlap each other.
  bool overlaps(region_t region_1, region_t region_2);
  
  //! Test if a position is within a region.
  bool overlaps(uint32_t pos_1, region_t region);
  
  //! Return overlapping regions, defaults to all regions.
  regions_t regions(uint32_t start_1 = 0, uint32_t end_1 = 0);
  
  //! Remove regions.
  cFlaggedRegions& remove(regions_t regions);
  
protected:
  regions_t m_regions;
};

}
#endif
