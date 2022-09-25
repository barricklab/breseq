/*****************************************************************************
 
 AUTHORS
 
 Jeffrey E. Barrick <jeffrey.e.barrick@gmail.com>
 David B. Knoester
 
 LICENSE AND COPYRIGHT
 
 Copyright (c) 2008-2010 Michigan State University
 Copyright (c) 2011-2022 The University of Texas at Austin
 
 breseq is free software; you can redistribute it and/or modify it under the
 terms the GNU General Public License as published by the Free Software
 Foundation; either version 1, or (at your option) any later version.
 
 *****************************************************************************/

#ifndef _BRESEQ_GENOME_DIFF_ENTRY_H_
#define _BRESEQ_GENOME_DIFF_ENTRY_H_

#include "common.h"
#include "file_parse_errors.h"
#include "counted_ptr.h"

using namespace std;

namespace breseq {

  class cDiffEntry;
  class cReferenceCoordinate;
  class cReferenceSequences;
  class cAnnotatedSequence;
  
  // Defined in other files
  extern const int32_t kBreseq_size_cutoff_AMP_becomes_INS_DEL_mutation;
  extern const int32_t kBreseq_large_mutation_size_cutoff;
  
  
  // Common keywords used for diff entries:
  extern const char* SEQ_ID;
  extern const char* START;
  extern const char* END;
  extern const char* STRAND;
  extern const char* POSITION;
  extern const char* INSERT_POSITION;
  extern const char* PHYLOGENY_ID;
  extern const char* FREQUENCY;

  // REJECT = comma delimited list of why evidence was rejected according to statistical criteria
  // Evidence is not rejected but is masked because it should not be used to predict mutations
  extern const char* REJECT;
  extern const char* POLYMORPHISM_REJECT;
  extern const char* CONSENSUS_REJECT;

  // IGNORE = comma delimited list of why passing evidence was not used to predict mutations
  // This can be b/c it is near a contig_end or is only evidence of a known CIRCULAR_CHROMOSOME
  extern const char* IGNORE;
  extern const char* USER_DEFINED;

  // For relationships to repeats and IS
  extern const char* MEDIATED;
  extern const char* BETWEEN;
  
  // For APPLY
  extern const char* WITHIN;
  extern const char* BEFORE;
  extern const char* APPLY_SIZE_ADJUST;
  
  // For DEL
  extern const char* SIZE;
  
  // For INS
  extern const char* NEW_SEQ;
  
  // For MOB
  extern const char* REPEAT_NAME;
  extern const char* DUPLICATION_SIZE;
  extern const char* INS_START;
  extern const char* INS_END;
  extern const char* DEL_START;
  extern const char* DEL_END;
  extern const char* MOB_REGION;
  
  // For INS/DEL
  extern const char* REPEAT_SEQ;
  extern const char* REPEAT_LENGTH;
  extern const char* REPEAT_REF_COPIES;
  extern const char* REPEAT_NEW_COPIES;
  
  // For AMP
  extern const char* NEW_COPY_NUMBER;
  extern const char* MEDIATED_STRAND;
  
  // For CON
  extern const char* REPLACE_SIZE;
  extern const char* REGION;
  
  //For RA
  // old + new required field
  extern const char* REF_BASE;
  // old fields to maintain RA definition
  extern const char* NEW_BASE;
  extern const char* REF_COV;
  extern const char* NEW_COV;
  
  // new fields
  extern const char* MAJOR_BASE;
  extern const char* MINOR_BASE;
  extern const char* MAJOR_COV;
  extern const char* MINOR_COV;
  extern const char* TOTAL_COV;
  extern const char* PREDICTION;
  extern const char* SCORE;
  extern const char* CONSENSUS_SCORE;
  extern const char* POLYMORPHISM_SCORE;
  extern const char* POLYMORPHISM_FREQUENCY;
  extern const char* MAJOR_FREQUENCY;
  extern const char* POLYMORPHISM_EXISTS;
  
  //For MC
  extern const char* START_RANGE;
  extern const char* END_RANGE;
  extern const char* LEFT_OUTSIDE_COV;
  extern const char* LEFT_INSIDE_COV;
  extern const char* RIGHT_INSIDE_COV;
  extern const char* RIGHT_OUTSIDE_COV;
  
  //For JC
  extern const char* SIDE_1_SEQ_ID;
  extern const char* SIDE_1_POSITION;
  extern const char* SIDE_1_STRAND;
  extern const char* SIDE_1_REDUNDANT;
  
  extern const char* SIDE_2_SEQ_ID;
  extern const char* SIDE_2_POSITION;
  extern const char* SIDE_2_STRAND;
  extern const char* SIDE_2_REDUNDANT;
  
  extern const char* OVERLAP;
  extern const char* UNIQUE_READ_SEQUENCE;
  
  extern const char* SIDE_1_READ_COUNT;
  extern const char* SIDE_2_READ_COUNT;
  extern const char* NEW_JUNCTION_READ_COUNT;
  extern const char* NEW_JUNCTION_FREQUENCY;
  
  extern const char* SIDE_1_COVERAGE;
  extern const char* SIDE_2_COVERAGE;
  extern const char* NEW_JUNCTION_COVERAGE;
  
  //For CN
  extern const char* COPY_NUMBER;
  
  extern const vector<string> gd_entry_type_lookup_table;

  // Types of diff entries:
  enum gd_entry_type {UNKNOWN = 0, SNP, SUB, DEL, INS, MOB, AMP, INV, CON, INT,
    RA, MC, JC, CN, UN, CURA, FPOS, PHYL, TSEQ, PFLP, RFLP, PFGE, NOTE, MASK};
  
  extern const vector<string> gd_keys_with_ids;
  
  inline string to_string(const gd_entry_type type)
  { return gd_entry_type_lookup_table[type]; }
  
  /*! Genome diff entry type.
   
   Instead of trying to define (and maintain!) concrete classes for each different
   type of diff entry, we instead define a single generic diff entry, and then
   provide easy-to-use mechanisms to manipulate them.
   
   The key abstraction here is that a diff entry is a map of keys (strings) to
   values.  In this case, values are strings that represent int, double, string,
   formatted_double and pair<,> types.
   
   
   */
  
  typedef string diff_entry_key_t; //!< Diff entry keys.
  typedef string diff_entry_value_t; //!< Diff entry values.
  typedef map<diff_entry_key_t, diff_entry_value_t> diff_entry_map_t; //!< Diff entry key-value map.
  
  typedef counted_ptr<cDiffEntry> diff_entry_ptr_t;
  typedef list<diff_entry_ptr_t> diff_entry_list_t; //!< Type for a list of diff entries.
  
  class cDiffEntry : public diff_entry_map_t {
  public:
    
    //!---- Sorting Items in Genome Diff ---- !//
    
    //! Genome Diff Sorting
    //! For sorting by a number, then by fields to break ties
    struct sort_fields_item {
      
      //! Constructor.
      sort_fields_item() {_f1=0; _f2=""; _f3=""; _f4="";};
      
      sort_fields_item(uint8_t f1, string f2, string f3, string f4 = "") :
      _f1(f1), _f2(f2), _f3(f3), _f4(f4) {};
      
      //! Destructor.
      virtual ~sort_fields_item() {};
      
      uint8_t _f1;
      string _f2;
      string _f3;
      string _f4;
    };
    
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
    
    //! Returns whether an id string is valid (>1 and an int)
    static bool valid_id(string& test_id);
    
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
    
    //! Comparison operators
    
    static int32_t compare(const cDiffEntry& a, const cDiffEntry& b);
    // Returns -1 if a < b, 0 if a == b, and +1 is a > b
    
    bool operator<(const cDiffEntry& de) const { return (compare(*this, de) < 0); }
    bool operator<=(const cDiffEntry& de) const { return (compare(*this, de) <= 0); }
    bool operator>(const cDiffEntry& de) const { return (compare(*this, de) > 0); }
    bool operator>=(const cDiffEntry& de) const { return (compare(*this, de) >= 0); }
    bool operator==(const cDiffEntry& de) const { return (compare(*this, de) == 0); }
    bool operator!=(const cDiffEntry& de) const { return (compare(*this, de) != 0); }

    //! Return if a given key value exists in _fields
    bool entry_exists(const diff_entry_key_t& k) const { return (count(k) > 0); }
    
    //! Return if this entry should be printed (does not begin with an underscore)
    static const string unprintable_key_prefix;
    static bool is_unprintable_key(const diff_entry_key_t& k)
    { return k.find(unprintable_key_prefix) == 0; }
    static string make_unprintable_key(const diff_entry_key_t & k)
    { return unprintable_key_prefix + k; }
    
    //! Return if this diff entry is a mutation
    bool is_mutation() const;
    
    //! Return if this diff entry is evidence
    bool is_evidence() const;
    
    //! Return if this diff entry is a validation
    bool is_validation() const;
    
    bool is_marked_deleted()
    { return (this->entry_exists("deleted") && ((from_string<int32_t>((*this)["deleted"]) != 0))); }
    
    bool is_polymorphism()
    {
      return (
              this->entry_exists("frequency")
              && ((from_string<double>((*this)["frequency"]) != 1))
              && ((from_string<double>((*this)["frequency"]) != 0))
      );
    }

    
    //!---- Accessors to calculated properties ---- !//
    
    //! Common function for getting start or end of mutation or evidence
    cReferenceCoordinate get_reference_coordinate_start() const;
    cReferenceCoordinate get_reference_coordinate_end() const;
    
    //! Is this item located within another and on same seq_id? Testing by location.
    bool located_within(const cDiffEntry &within) const;
    
    //! Common function giving change in size of genome at site of applying entry
    int32_t mutation_size_change(cReferenceSequences& ref_seq_info) const;
    
    
    //!---- Functions for updating mutations ---- !//
    
    //! Common function for updating mutations based on a mutation occurring in the interval shift_start to shift_end and changing size by shift_size;
    void mutation_shift_position(const string& seq_id, const cReferenceCoordinate& shift_start, const cReferenceCoordinate& shift_end, int32_t shift_size);
    
    void mutation_reverse_complement();
    
    // Updates positions for inversion and reverse-complements mutation
    void mutation_invert_position_sequence(cDiffEntry& inverting_mutation);
    
    bool is_small_mutation(uint32_t large_size_cutoff=kBreseq_large_mutation_size_cutoff);
    
    // Updates
    void annotate_repeat_hotspots(cReferenceSequences& new_ref_seq_info, int32_t slop_distance, int32_t size_cutoff_AMP_becomes_INS_DEL_mutation, bool remove_old_tags, bool warn_mode = false);
    
    //!---- Functions related to complex APPLY situations ---- !//
    
    //! Does this item have a "within" key for a certain mutation? Returns empty string if not.
    const string get_within_mutation_id() const
    {
      if (!this->entry_exists(WITHIN)) return "";
      vector<string> split_within = split(this->get(WITHIN), ":");
      return split_within[0];
    }
    
    //!---- Output ---- !//
    
    //! Marshal this diff entry into an ordered list of fields.
    virtual void marshal(vector<string>& s, bool include_unprintable_fields=false) const;
    
    //! Serialize this diff entry into a string for output.
    virtual string as_string(bool include_unprintable_fields=false) const;
    
    //! Output all keys and values
    string as_key_values() const {
      string s;
      for (diff_entry_map_t::const_iterator it = this->begin(); it != this->end(); it++)
      {
        s+= it->first + " = " + it->second + "\n";
      }
      return s;
    }
    
    //!---- Reject Reasons Field ---- !//
    
    size_t number_reject_reasons();
    
    bool is_rejected() { return (number_reject_reasons() > 0); }
    
    //! Returns values for cDiffEntry["reject"]
    vector<string> get_reject_reasons(const string& field=REJECT);
    
    //! Adds a reject reason to cDiffEntry["reject"] as a list
    void add_reject_reason(const string &reason, const string& field=REJECT);
    
    //! Removes just this reject reason
    void remove_reject_reason(const string &reason, const string& field=REJECT);
    
    //! Removes all of the reject reasons
    void clear_reject_reasons(const string& field=REJECT)
    { this->erase(field); }
    
    bool is_rejected_and_not_user_defined()
    { return entry_exists(REJECT) && !entry_exists(USER_DEFINED); }

    bool is_ignored_and_not_user_defined()
    { return entry_exists(IGNORE) && !entry_exists(USER_DEFINED); }
    
    //!---- Simplifying entries ---- !//
    
    //! Remove all information except required fields
    cDiffEntry to_spec(void) const; //void returns stripped, leaving this unchanged
    void strip_to_spec(); // strips this item
    
    //! @JEB 03-16-2014 should deprecate. Only used in random mutation generator, which should be handled by MutationPredictor
    //  functionality is replaced by cGenomeDiff::normalize_mutations.
    void normalize_to_sequence(const cAnnotatedSequence &seq, bool verbose = false);
    
    //!---- Sorting ---- !//
    
    //! Functor. Sorts cDiffEntrys in descending order depending on given fields that
    //can be evaluated as an unsigned integer.
    struct descending_by_scores : public binary_function
    <diff_entry_ptr_t, diff_entry_ptr_t, bool>
    {
      
      //! Constructor
      explicit descending_by_scores (const vector<diff_entry_key_t>& field_keys)
      : m_field_keys(field_keys) {}
      
      //! Predicate
      virtual bool operator() (const diff_entry_ptr_t& a, const diff_entry_ptr_t& b) const
      {
        for (vector<diff_entry_key_t>::const_iterator itr = m_field_keys.begin(); itr != m_field_keys.end(); itr++) {
          string key(*itr);
          
          double a_val = -9999;
          double b_val = -9999;
          
          if (a->entry_exists(key) && ((*a)[key].size() > 0))
            a_val = from_string<double>((*a)[key]);
          if (b->entry_exists(key) && ((*b)[key].size() > 0))
            b_val = from_string<double>((*b)[key]);
          
          if (a_val == b_val)
            continue;
          else
            return a_val > b_val;
        }
        return false;
      }
      
    protected:
      vector<diff_entry_key_t> m_field_keys;
    };
    
    //!---- Filtering ---- !//
    
    //! Functors for selecting certain types of entries
    
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
    
    struct ignored_and_not_user_defined:public unary_function<diff_entry_ptr_t,bool> {
      virtual bool operator() (diff_entry_ptr_t cDiffEntry)
      {
        return cDiffEntry->entry_exists(IGNORE) && !cDiffEntry->entry_exists("user_defined");
      }
    };
    
    struct is_not_consensus:public unary_function<diff_entry_ptr_t,bool> {
      virtual bool operator() (diff_entry_ptr_t cDiffEntry)
      {
        return from_string<double>((*cDiffEntry)[FREQUENCY]) != 1.0;
      }
    };
    
    struct ignored_but_not_circular:public unary_function<diff_entry_ptr_t,bool> {
      virtual bool operator() (diff_entry_ptr_t cDiffEntry)
      {
        if ( (*cDiffEntry).entry_exists(IGNORE) ) {
          return (*cDiffEntry)[IGNORE] != "CIRCULAR_CHROMOSOME";
        }
        return false;
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
    //Returns true if a cDiffEntry contains the given field_key.
    struct field_equals : public unary_function <diff_entry_ptr_t, bool>
    {
      //! Constructor
      explicit field_equals (const string& field_key, const string& field_value)
      : m_field_key(field_key), m_field_value(field_value) {}
      
      //! Predicate
      virtual bool operator() (const diff_entry_ptr_t& p_diff_entry)
      const {return (*p_diff_entry)[m_field_key] == m_field_value;}
      
    protected:
      diff_entry_key_t m_field_key;
      string m_field_value;
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
  

}
#endif
