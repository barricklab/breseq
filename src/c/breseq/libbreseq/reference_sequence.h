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

#ifndef _BRESEQ_REFERENCE_SEQUENCE_H_
#define _BRESEQ_REFERENCE_SEQUENCE_H_

#include "libbreseq/common.h"

#include "libbreseq/genome_diff.h"
#include "libbreseq/fasta.h"
#include "libbreseq/fastq.h"
#include "libbreseq/alignment.h"
#include "libbreseq/anyoption.h"
#include "libbreseq/summary.h"


namespace breseq {


  // Pre-declaration
	class Settings;
  
  // Convenience structure for passing many options
  class MutationTableOptions {
  public:
    MutationTableOptions(const Settings& _settings);
    
    uint32_t repeat_header;
    bool legend_row;
    bool force_show_sample_headers;
    bool one_ref_seq;
    bool force_frequencies_for_one_reference;
    bool shade_frequencies;
    bool detailed;
    vector<string> gd_name_list_ref;
    string relative_link;
    uint32_t max_nucleotides_to_show_in_tables;
    bool no_javascript;
  };

  /* Position for dealing with things that are relative to reference genom
    and the annoying issues of insert counts
  */
  class cReferenceCoordinate {
  private:
    int32_t m_position;
    int32_t m_insert_position;
  public:
    
    cReferenceCoordinate(): m_position(0), m_insert_position(0) { }
    
    cReferenceCoordinate(const int32_t position, const int32_t insert_position = 0)
    : m_position(position), m_insert_position(insert_position) { }
    
    int32_t get_position() const
    { return m_position; }
    
    int32_t get_insert_position() const
    { return m_insert_position; }
    
    bool operator<(const cReferenceCoordinate& _in) const
    {
      if (this->m_position == _in.m_position) {
          return (this->m_insert_position < _in.m_insert_position);
      }
      return (this->m_position < _in.m_position);
    }
    
    bool operator>(const cReferenceCoordinate& _in) const
    {
      if (this->m_position == _in.m_position) {
        return (this->m_insert_position > _in.m_insert_position);
      }
      return (this->m_position > _in.m_position);
    }
    
    bool operator<=(const cReferenceCoordinate& _in) const
    {
      if (this->m_position == _in.m_position) {
        return (this->m_insert_position <= _in.m_insert_position);
      }
      return (this->m_position <= _in.m_position);
    }
    
    bool operator>=(const cReferenceCoordinate& _in) const
    {
      if (this->m_position == _in.m_position) {
        return (this->m_insert_position >= _in.m_insert_position);
      }
      return (this->m_position >= _in.m_position);
    }
    
    bool operator==(const cReferenceCoordinate& _in) const
    {
      return (this->m_position == _in.m_position) && (this->m_insert_position == _in.m_insert_position);
    }
    
    cReferenceCoordinate& operator+=(const cReferenceCoordinate& rhs)
    {
      this->m_position += rhs.get_position();
      this->m_insert_position += rhs.get_insert_position();
      
      /* addition of rhs to *this takes place here */
      return *this; // return the result by reference
    }
    
    friend cReferenceCoordinate operator+(cReferenceCoordinate lhs,
                                          const cReferenceCoordinate& rhs)
    {
      lhs += rhs; // reuse compound assignment
      return lhs; // return the result by value (uses move constructor)
    }
    
    cReferenceCoordinate& operator-=(const cReferenceCoordinate& rhs)
    {
      this->m_position -= rhs.get_position();
      this->m_insert_position -= rhs.get_insert_position();
      
      /* addition of rhs to *this takes place here */
      return *this; // return the result by reference
    }
    
    friend cReferenceCoordinate operator-(cReferenceCoordinate lhs,
                                          const cReferenceCoordinate& rhs)
    {
      lhs -= rhs; // reuse compound assignment
      return lhs; // return the result by value (uses move constructor)
    }
  };
  
  // Pre-declaration
  class cSequenceFeature;

  class cLocation {
    
  private:
    int32_t m_start_1, m_end_1;
    int8_t m_strand;
    
    // for when coords are marked as extending last here (<1..514)
    bool m_start_is_indeterminate, m_end_is_indeterminate;
    
  public:
    cLocation()
    : m_start_1(0)
    , m_end_1(0)
    , m_strand(0)
    , m_start_is_indeterminate(false)
    , m_end_is_indeterminate(false)
    {};
    
    cLocation(
             int32_t start_1,
             int32_t end_1,
             int8_t strand,
             const cSequenceFeature* feature = NULL,
             bool start_is_indeterminate = false,
             bool end_is_indeterminate = false
             )
    : m_start_1(start_1)
    , m_end_1(end_1)
    , m_strand(strand)
    , m_start_is_indeterminate(start_is_indeterminate)
    , m_end_is_indeterminate(end_is_indeterminate)
    {
      this->check_valid(feature);
    }
    
    cLocation(const cLocation& in)
    : m_start_1(in.m_start_1)
    , m_end_1(in.m_end_1)
    , m_strand(in.m_strand)
    , m_start_is_indeterminate(in.m_start_is_indeterminate)
    , m_end_is_indeterminate(in.m_end_is_indeterminate)
    {
      this->check_valid();
    }
    
    //>! For sorting
    bool operator < (const cLocation& in) const
    {
      if (m_start_1 != in.m_start_1)
        return (m_start_1 < in.m_start_1);
      
      return (m_end_1 < in.m_end_1);
    }
    
    bool operator > (const cLocation& in) const
    {
      return (in < *this);
    }
    
    bool operator != (const cLocation& in) const
    {
      return (*this < in) || (*this > in);
    }
    
    bool operator == (const cLocation& in) const
    {
      return !(*this < in) && !(*this > in);
    }
    
    //>! Get 1-indexed start position
    int32_t get_start_1() const {
      return m_start_1;
    }
    //>! Get 1-indexed end position
    int32_t get_end_1() const {
      return m_end_1;
    }
    
    //>! Get whether start position is indeterminate
    bool start_is_indeterminate() const {
      return m_start_is_indeterminate;
    }
    //>! Get whether end position is indeterminate
    bool end_is_indeterminate() const {
      return m_end_is_indeterminate;
    }
    
    //>! Get whether start position is indeterminate
    bool stranded_start_is_indeterminate() const {
      return (m_strand==1) ? m_start_is_indeterminate : m_end_is_indeterminate;
    }
    //>! Get whether end position is indeterminate
    bool stranded_end_is_indeterminate() const {
      return (m_strand==1) ? m_end_is_indeterminate : m_start_is_indeterminate;
    }
    
    //>! Strand is -1 or +1
    int8_t get_strand() const { return m_strand; }
    
    bool is_top_strand() const { return m_strand==1; }
    
    int32_t get_strand_aware_initial_position_1() const {
      return (m_strand == 1) ? m_start_1 : m_end_1;
    }
    
    string as_string() {
      cString s;
      s = "start = " + to_string(m_start_1) + " end = " + to_string(m_end_1) + " strand " + to_string<int32_t>(m_strand);
      return std::move(s);
    }
    
    bool is_valid() {
      return (m_start_1 <= m_end_1) && (m_strand >=-1) && (m_strand <=1);
    }
    
    void check_valid(const cSequenceFeature* feature = NULL);
    
    void set_start_1(int32_t start_1) {
      m_start_1 = start_1;
      this->check_valid();
    }
    void set_end_1(int32_t end_1) {
      m_end_1 = end_1;
      this->check_valid();
    }
    
    void offset(int32_t shift) {
      m_start_1 += shift;
      m_end_1 += shift;
    }
    
    // This function is provided to not transiently create unallowed start/end combinations
    void set_start_end_1(int32_t start_1, int32_t end_1) {
      m_start_1 = start_1;
      m_end_1 = end_1;
      this->check_valid();
    }
    
    void set_start_is_indeterminate(bool start_is_indeterminate) {
      m_start_is_indeterminate = start_is_indeterminate;
    }
    void set_end_is_indeterminate(bool end_is_indeterminate) {
      m_end_is_indeterminate = end_is_indeterminate;
    }
    
    void set_strand(int8_t _strand) {
      m_strand = _strand;
    }
    
    int32_t distance_to_position(const int32_t pos) const {
      
      if ((pos >= m_start_1) && (pos <= m_end_1))
        return 0;
      
      return min(abs(pos - m_start_1), abs(pos - m_end_1));
    }
    
    void invert_within_region(const int32_t invert_start_1, int32_t invert_end_1) {
      int32_t original_start_1 = get_start_1();
      int32_t original_end_1 = get_end_1();
      set_start_end_1(
                      invert_start_1 +  (invert_end_1 - original_end_1),
                      invert_start_1 + (invert_end_1 - original_start_1)
                      );
      set_strand(-get_strand());
    }
    
    string as_string() const {
      return to_string(get_start_1()) + "-" + to_string(get_end_1()) + " " + to_string<int32_t>(get_strand());
    }
    
  };
  
  
  class cSequenceFeature;
  
  /*! We allow a feeature to consist of an ordered list of regions
   * For most features, there is only one location
   * this enables representing spliced and even trans spliced genes
   */
  
  class cFeatureLocation : public cLocation {
    
  private:
    cSequenceFeature* m_feature;
    
    // additional information for linking together sublocations
    uint32_t m_index;   // index of this exon or piece within the feature
    int32_t m_offset;   // within-gene location of the start of the location
                        // 1 for the first piece
    
  public:
    cFeatureLocation()
      : cLocation()
      , m_feature(NULL)
      , m_index(0)
      , m_offset(1)
    {};
    
    cFeatureLocation(cSequenceFeature* feature, const cLocation& loc);
    
    // This constructor takes into account
    // any previously assigned locations in feature
    // to automatically fill in m_index and m_offset
    cFeatureLocation(
                    cSequenceFeature* feature,
                    int32_t start_1,
                    int32_t end_1,
                    int8_t strand,
                    bool start_is_indeterminate = false,
                    bool end_is_indeterminate = false
                    );

    
    //>! Get sort index
    int32_t get_index() const {
      return m_index;
    }
    
    cSequenceFeature* get_feature() {
      ASSERT(m_feature, "Feature location does not have associated feature.");
      return m_feature;
    }
    
    void set_feature(cSequenceFeature* feature) {
      m_feature = feature;
    }
    
  };
  
class cFeatureLocationList: public list<cFeatureLocation> {
  
public:
  string as_string(string sep=",") {
    vector<string> ret;
    for(list<cFeatureLocation>::iterator it=this->begin(); it!=this->end(); it++) {
      ret.push_back(it->as_string());
    }
    return join(ret, sep);
  }
  
};
  
  extern const vector<string> snp_types;
    
	/*! Sequence Feature class.
	 */ 

  // Currently everything is stored as strings...
  typedef string sequence_feature_key_t;
  typedef string sequence_feature_value_t;
  typedef map<sequence_feature_key_t, sequence_feature_value_t> sequence_feature_map_t; //!< Diff entry key-value map.
  
  class cSequenceFeature : public sequence_feature_map_t {
    
    // Required Fields
    //
    // GFF3: http://www.sequenceontology.org/gff3.shtml
    // GenBank: http://www.ncbi.nlm.nih.gov/collab/FT/
    //
    //   "type"  GenBank (primary key) |  GFF3 (type)
    //      string [CDS, gene, repeat_region, tRNA, rRNA]
    //
    //   m_start, m_end, m_strand GenBank (location) | GFF (start, end, strand)
    //      uint32_t and -1/+1
    //
    //  "name" GenBank (gene, locus_tag) | GFF (name)
    //      string
    //
    //  "alias" GenBank (locus_tag) | GFF (alias)
    //      string
    //  
    
    public:
    
      // Could add accessors that convert strings to numbers...
      cFeatureLocationList m_locations;
      bool m_start_is_indeterminate;
      bool m_end_is_indeterminate;
      bool m_pseudo;
      map<string, vector<string> > m_gff_attributes;
      vector<string> m_original_genbank_tags;
    
      cSequenceFeature() : m_start_is_indeterminate(false), m_end_is_indeterminate(false), m_pseudo(false) {}

      cSequenceFeature(const cSequenceFeature& copy) 
        : sequence_feature_map_t(copy)
        , m_locations(copy.m_locations)
        , m_start_is_indeterminate(copy.m_start_is_indeterminate)
        , m_end_is_indeterminate(copy.m_end_is_indeterminate)
        , m_pseudo(copy.m_pseudo)
        , m_gff_attributes(copy.m_gff_attributes)
        , m_original_genbank_tags(copy.m_original_genbank_tags)
      {
        // need to update locations to point to this feature
        for (list<cFeatureLocation>::iterator it=m_locations.begin(); it!=m_locations.end(); it++) {
          it->set_feature(this);
        }
      }

      bool operator<(const cSequenceFeature& in) const
      {
        list<cFeatureLocation>::const_iterator it1 = this->m_locations.begin();
        list<cFeatureLocation>::const_iterator it2 = in.m_locations.begin();
        
        // Smallest start positions
        while ((it1 != this->m_locations.end()) && (it2 != in.m_locations.end())) {
          if (*it1 != *it2) {
            return it1->get_start_1() < it2->get_start_1();
          }
          it1++; it2++;
        } ;
        
        /*
        // Sort 'source' and 'region' entries to the front of the list
        if ( ((this->SafeGet("type") == "source") || (this->SafeGet("type") == "region")) && (in.SafeGet("type") != "source") && (in.SafeGet("type") != "region") )
          return true;
        if ( (this->SafeGet("type") != "source") && (this->SafeGet("type") != "region") && ((in.SafeGet("type") == "source") || (in.SafeGet("type") == "region")) )
          return false;
        */
        
        // Different number of positions? Fewest first.
        if (! ((it1 == this->m_locations.end()) && (it2 == in.m_locations.end())) ) {
          if (it1 == this->m_locations.end())
            return true;
          if (it2 == in.m_locations.end())
            return false;
        }
        
        // Different length? Longer first.
        int32_t l1 = this->get_length();
        int32_t l2 = in.get_length();
        if (l1 != l2) {
          return l1 > l2;
        }

        
        // Break ties by type. It's important to put genes first before CDS
        if ( (this->SafeGet("type") == "gene") && (in.SafeGet("type") != "gene"))
          return true;
        if ( (this->SafeGet("type") != "gene") && (in.SafeGet("type") == "gene") )
          return false;
        
        return (this->SafeGet("type") < in.SafeGet("type"));
      }
    
      //<! Safe accessor that returns empty string if not defined.
      string SafeGet(sequence_feature_key_t in_key) const 
      { 
        sequence_feature_map_t::const_iterator it = this->find(in_key);
        if (it == this->end()) return std::string("");
        return it->second;
      }
  
      bool is_source() {
        return ((*this)["type"] == "region") || ((*this)["type"] == "source");
      }
    
      bool is_repeat() {
        return (   ((*this)["type"] == "mobile_element")
                || ((*this)["type"] == "repeat_region")
        );
      }
    
      bool is_gene() {
        return (   ((*this)["type"] == "CDS")
                || ((*this)["type"] == "rRNA")
                || ((*this)["type"] == "tRNA")
                || ((*this)["type"] == "ncRNA")
                || ((*this)["type"] == "RNA")
        );
      }
    
      string get_locus_tag() const {
        if (m_gff_attributes.count("Alias")) {
          return m_gff_attributes.at("Alias")[0];
        }
        else if (this->count("locus_tag")) {
          return this->at("locus_tag");
        }
        return "";

      }

    //Mark it as pseudo
      void flag_pseudo(bool verbose=false)
      {
        //If this feature is already pseudo or a region, do nothing.
        if(m_pseudo || (*this)["type"] == "region" || (*this)["type"] == "source" || (*this)["type"] == "repeat_region")
          return;
        
        //Notify the user of the action
        if(verbose) cout << "PSEUDO\t" << (*this)["type"] << "\t" << m_gff_attributes["ID"];
        
        m_pseudo = true;
        m_gff_attributes["Pseudo"].push_back("true");
        
        //Notify the user of the action (cont.)
        if(verbose) cout << endl;
      }
    
      // Adds append_str to all important IDs and gene names
      void append_to_accession(const string& append_str)
      {
        //(*this)["name"] += append_str;
        (*this)["accession"] += append_str;
        if (m_gff_attributes.count("Alias")) m_gff_attributes["Alias"][0] += append_str;
        if (m_gff_attributes.count("ID")) m_gff_attributes["ID"][0] += append_str;
        //if (m_gff_attributes.count("Name")) m_gff_attributes["Name"][0] += append_str;
      }

      int32_t get_length() const {
        int32_t total_length(0);
        for(cFeatureLocationList::const_iterator it=m_locations.begin(); it!=m_locations.end(); ++it) {
          const cFeatureLocation& loc = *it;
          total_length += loc.get_end_1() - loc.get_start_1() + 1;
        }
        return total_length;
      }
    
      void add_location(const cLocation& loc) {
        this->m_locations.push_back(cFeatureLocation(this, loc));
        
        // Update indeterminate start and end
        if (this->m_locations.size()==1) {
          m_start_is_indeterminate = loc.stranded_start_is_indeterminate();
        }
        m_end_is_indeterminate = loc.stranded_end_is_indeterminate();
      }
    
      void ReadGenBankTag(string& tag, string& s, ifstream& in);
    
      string get_nucleotide_sequence(const cAnnotatedSequence& seq) const;
    
      void genomic_position_to_index_strand_1(int32_t pos_1, int32_t& index_1, int8_t& strand) const;
    
      bool start_is_indeterminate() { return m_start_is_indeterminate; }
      bool end_is_indeterminate() { return m_end_is_indeterminate; }
    
      void make_feature_strings_safe();
  };
  
  //!< Subclass of reference features with more information
  class cGeneFeature : public cSequenceFeature {
  public:
    string name;
    string product;
    string type;
    bool pseudogene; 
    uint32_t translation_table;
    
    cGeneFeature() {};
    cGeneFeature(cSequenceFeature& src) : cSequenceFeature(src)
    {
      name = src["name"];
      product = src["product"];
      type = src["type"];
      pseudogene = src.m_pseudo;
      translation_table = 1;
      if (src.count("transl_table")) 
        translation_table = from_string<uint32_t>(src["transl_table"]);
    }

    cGeneFeature(const cGeneFeature& copy) : cSequenceFeature(copy)
    {
      name = copy.name;
      product = copy.product;
      type = copy.type;
      pseudogene = copy.pseudogene;
      translation_table = copy.translation_table;
    }


  };
  
  typedef counted_ptr<cGeneFeature> cGeneFeaturePtr;
  typedef counted_ptr<cSequenceFeature> cSequenceFeaturePtr;
  typedef list<cSequenceFeaturePtr> cSequenceFeatureList;
  
	/*! Sequence class.
	 */   
   
  class cAnnotatedSequence {
    
    public:      
      int32_t m_length;
      bool m_is_circular;
      string m_description; // GenBank (DEFINITION) | GFF (description), from main feature line
      string m_seq_id;      // GenBank (LOCUS)      | GFF (seqid), from ##sequence-region line
      string m_file_name;   // Name of file this sequence was loaded from
      string m_file_format; // Format of file used to load sequence

      cFastaSequence m_fasta_sequence;            //!< Nucleotide sequence
    
      vector<string> m_genbank_raw_header_lines;   //!< Raw GenBank header linex (except LOCUS and FEATURES)
 
      // Features are stored as counted pointers so that we can have ready-made lists
      // of different types of features. 
      cSequenceFeatureList m_features;    //!< Full list of sequence features
    
      // These lists are subsets of the original list
      cSequenceFeatureList m_genes;       //!< Subset of features
      cSequenceFeatureList m_repeats;     //!< Subset of features
    
      // Storage as list of locations is used for identifying features
      // that overlap a certain position
      cFeatureLocationList m_gene_locations;
      cFeatureLocationList m_repeat_locations;
    
      string m_features_loaded_from_file; //!< File name features were loaded from
      string m_sequence_loaded_from_file; //!< File name sequence was loaded from
    
    public:
    
      //Constructor for empty object
      cAnnotatedSequence() : 
        m_length(0),
        m_is_circular(false),
        m_description("na"), 
        m_seq_id("na")
        {} ;
    
      void set_features_loaded_from_file(const string& file_name, bool allow_reload = false) {
        ASSERT(allow_reload || (m_features_loaded_from_file.size() == 0), "Duplicate information for sequence found in file '" + file_name + "'!\nFeatures for sequence '" + m_seq_id + "' were already loaded from file '" + m_features_loaded_from_file + "'.")
        m_features_loaded_from_file = file_name;
      }
    
      void set_sequence_loaded_from_file(const string& file_name, bool allow_reload = false) {
        ASSERT(allow_reload || (m_sequence_loaded_from_file.size() == 0), "Duplicate information for sequence found in file '" + file_name + "'!\nDNA sequence for sequence '" + m_seq_id + "' was already loaded from file '" + m_sequence_loaded_from_file + "'.")
        m_sequence_loaded_from_file = file_name;
      }
    
      void set_file_format(const string& file_format) {
        if (m_file_format.length() > 0) m_file_format += "+";
        m_file_format += file_format;
      }
    
      void sort_features()
      {
        m_features.sort();
        m_genes.sort();
        m_repeats.sort();
      }
    
      bool is_circular() const
      {
        return m_is_circular;
      }
    
      // Utility to get top strand sequence
      string get_sequence_1(int32_t start_1, int32_t end_1) const
      {
        return m_fasta_sequence.get_sequence_1(start_1, end_1);
      }
    
      char get_sequence_1(int32_t pos_1) const
      {
        return m_fasta_sequence.get_sequence_1(pos_1);
      }
    
      string get_stranded_sequence_1(int32_t strand, int32_t start_1, int32_t end_1) const;
      char get_stranded_sequence_1(int32_t strand, int32_t pos_1) const;
  
      string get_sequence_1_start_size(int32_t start_1, uint32_t size) const;
      int32_t get_circular_distance_1(int32_t pos_1, int32_t pos_2) const;

    
      size_t get_sequence_length(void) const
      {
        return m_fasta_sequence.get_sequence_length();
      }

      // Replace Sequence with Input
      void replace_sequence_1(int32_t start_1, int32_t end_1, const string &replacement_seq, string mut_type="", bool verbose=false);
      
      // Insert Input AFTER Position
      void insert_sequence_1(int32_t pos_1, const string &insertion_seq, string mut_type="", bool verbose=false);
    
      // Invert sequence between two positions
      void invert_sequence_1(int32_t start_1, int32_t end_1, string mut_type="", bool verbose=false);
    
      // Repeat Feature at Position
      void repeat_feature_1(int32_t pos, int32_t start_del, int32_t end_del, cReferenceSequences& orig_ref_seq_info, string& orig_seq_id, int8_t strand, const cLocation &repeated_region, bool verbose = false);
      
      // Find Specific Feature
      // Given a cSequenceFeatureList feat_list, iterate through it until
      // finding the cSequenceFeature feat.
      // Iterators passed in should be used to quickly find the features
      bool find_feature(cSequenceFeatureList &feat_list, cSequenceFeature &feat, list<cSequenceFeaturePtr>::iterator &it)
      {        
        //Iterate through the list
        //If a matching feature is found return true
        for (it = feat_list.begin(); it != feat_list.end(); it++)
        {
          //The current feature we're looking at
          cSequenceFeature& temp_feat = **it;
          
          //If the feature in the list matches the the main feature
          //break out of here and return true.
          if(temp_feat == feat){return true;}
        }
        
        //If we got this far and found nothing, return false
        return false;
      }
    
      string get_file_name() const
      {
        return m_file_name;
      }
    
      void make_feature_strings_safe() {
        for (list<cSequenceFeaturePtr>::iterator it = m_features.begin(); it != m_features.end(); it++)
        {
          (*it)->make_feature_strings_safe();
        }
      }

      void update_feature_lists();
    
      static cAnnotatedSequence deep_copy(cAnnotatedSequence in) {
        cAnnotatedSequence copy;

        copy.m_length = in.m_length;
        copy.m_is_circular = in.m_is_circular;
        copy.m_description = in.m_description;
        copy.m_seq_id = in.m_seq_id;
        copy.m_fasta_sequence = in.m_fasta_sequence;
        copy.m_genbank_raw_header_lines = in.m_genbank_raw_header_lines;

        //Features.
        for (cSequenceFeatureList::iterator it = in.m_features.begin(); it != in.m_features.end(); ++it) {
          cSequenceFeaturePtr new_feature(new cSequenceFeature(**it));
          copy.m_features.push_back(new_feature);
        }
        copy.update_feature_lists();
        return copy;
      }
    
      // Read GenBank coords
      list<cLocation> ReadGenBankCoords(const cSequenceFeature& in_feature, string& s, ifstream& in, bool safe_create_feature_locations);
      //Parse portion of GenBank coords string
      list<cLocation> ParseGenBankCoords(const cSequenceFeature& in_feature, string& s, bool safe_create_feature_locations, int8_t in_strand = 1);
    
      list<cLocation> SafeCreateFeatureLocations(
                                          const cSequenceFeature& in_feature,
                                          int32_t in_start_1,
                                          int32_t in_end_1,
                                          int8_t in_strand,
                                          bool in_start_is_indeterminate,
                                          bool in_end_is_indeterminate,
                                          bool safe_create_feature_locations = true
                                          );
  };
  
  /*! Reference Sequences
   
   Holds sequences and features for ALL reference sequences.
	 */ 
  
  typedef map<string,uint32_t> str_uint;
  
  class cReferenceSequences : public vector<cAnnotatedSequence> 
  {
    
  protected:
    map<string,int> m_seq_id_to_index; // for looking up sequences by seq_id
    str_uint m_seq_id_loaded; // for recording how many times we tried to load this seq_id
    uint32_t m_index_id;
    bool m_initialized;
    bool m_use_safe_seq_ids;
    map<string,string> m_seq_id_to_original_file_name;

    //!< Currently supported file types.
    enum FileType {UNKNOWN, GENBANK, FASTA, GFF3, BULL};

  public:
    
    // Used for annotating mutations
    static const string intergenic_separator;
    static const string text_intergenic_separator;
    static const string html_intergenic_separator;
    
    static const string gene_list_separator;
    static const string text_gene_list_separator;
    static const string html_gene_list_separator;
    
    static const string no_gene_name;
    static const string gene_range_separator;
    
    static const string multiple_separator;
    static const string text_multiple_separator;
    static const string html_multiple_separator;
    
    static const string gene_strand_reverse_char;
    static const string gene_strand_forward_char;

    static const double k_inactivate_overlap_fraction;
    static const int32_t k_promoter_distance;
    
    cReferenceSequences(bool _use_safe_seq_ids = true)
      : m_index_id(0)
      , m_initialized(false)
      , m_use_safe_seq_ids(_use_safe_seq_ids)
    {}

    //!< Load all reference files and verify - this is the only public load method!
    void LoadFiles(const vector<string>& file_names, const string& genbank_field_for_seq_id = "AUTOMATIC");
    
    //!: Convenience function to load just one file
    void LoadFile(const string& file_name)
      { LoadFiles(make_vector<string>(file_name)); }
    
    //!< Fixes gene/product names so that our separator character is unique. Called after load.
    void make_feature_strings_safe() {
      for (vector<cAnnotatedSequence>::iterator its= this->begin(); its != this->end(); its++) {
        cAnnotatedSequence& as = *its;
        as.make_feature_strings_safe();
      }
    }
    
    //!< Updates gene lists and other properties. Called after load.
    void update_feature_lists() {
      for (vector<cAnnotatedSequence>::iterator its= this->begin(); its != this->end(); its++) {
        cAnnotatedSequence& as = *its;
        as.update_feature_lists();
      }
    }
    
    //!< Adds new IS_element entries
    void ReadISEScan(const string& isescan_csv_file_name);
    
  protected:
    
    //!< Read reference file - not safe to call on their own = private
    //!< because we need to do some pre- and post-load
    // !< initializations that are handled centrally in LoadFiles
    
    //!< Load reference file into object
    //!< Detect file type and load the information in it appropriately
    //!< Should not be called directly, only through LoadFiles!
    void PrivateLoadFile(const string& file_name, const string& genbank_field_for_seq_id = "AUTOMATIC");
    
    //!< Verify that all seq_id have sequence and that features fit in sequence;
    void VerifySequenceFeatureMatch();
    void VerifyCDSLengthsAreValid();
    bool Initialized() {return m_initialized;}
    
    void ReadFASTA(const std::string &file_name);
    void ReadFASTA(cFastaFile& ff);
    void ReadGFF(const string& file_name);
    void ReadBull(const string& file_name);

    //!< Read GenBank file
    void ReadGenBank(const string& in_file_names, const string& genbank_field_for_seq_id = "AUTOMATIC");
    bool ReadGenBankFileHeader(std::ifstream& in, const string& file_name, const string& genbank_field_for_seq_id = "AUTOMATIC");
    //void ReadGenBankTag(std::string& tag, std::string& s, std::ifstream& in);
    void ReadGenBankFileSequenceFeatures(std::ifstream& in, cAnnotatedSequence& s);
    void ReadGenBankFileSequence(std::ifstream& in, cAnnotatedSequence& s);
    
    //!< Write GenBank file
    void WriteGenBankFileHeader(std::ofstream& out, const cAnnotatedSequence& s);
    void WriteGenBankFileSequenceFeatures(std::ofstream& out, const cAnnotatedSequence& s);
    void WriteGenBankFileSequence(std::ofstream& out, const cAnnotatedSequence& s);
    
  public:
    void WriteFASTA(const string &file_name);
    void WriteFASTA(cFastaFile& ff);
    void WriteGFF(const string &file_name, bool no_sequence = false);
    void WriteGenBank(const string &file_name, bool no_sequence = false);
    void WriteCSV(const string &file_name);
    
    //!< Moves over original file names to the current names
    //!< This is key for the relaod of GFF that happens during a breseq
    //!< run maintaining junction-only and contig reference designations
    void use_original_file_names() {
      
      ASSERT(m_seq_id_to_original_file_name.size() == this->size(), "Number of original file names saved does not match the actual number of sequences.");
      
      for(vector<cAnnotatedSequence>::iterator it = this->begin(); it != this->end(); it++) {
        ASSERT(m_seq_id_to_original_file_name.find(it->m_seq_id) != m_seq_id_to_original_file_name.end(), "Could not find original file name for seq id: " + it->m_seq_id);
        it->m_file_name = m_seq_id_to_original_file_name[it->m_seq_id];
      }
    }
    
    //!< Calculates the total length of all reference sequences together
    uint64_t get_total_length() const
    {
      uint64_t ret_val(0);
      for (cReferenceSequences::const_iterator it = (*this).begin(); it != (*this).end(); it++)
      {
        ret_val += it->m_length;
      }
      return ret_val;
    }
    
    uint64_t get_total_num_genes() const
    {
      uint64_t ret_val(0);
      for (cReferenceSequences::const_iterator it = (*this).begin(); it != (*this).end(); it++)
      {
        ret_val += it->m_genes.size();
      }
      return ret_val;
    }
    
    uint64_t get_total_num_repeats() const
    {
      uint64_t ret_val(0);
      for (cReferenceSequences::const_iterator it = (*this).begin(); it != (*this).end(); it++)
      {
        ret_val += it->m_repeats.size();
      }
      return ret_val;
    }
    
    string get_file_formats() const
    {
      set<string> formats;
      for (cReferenceSequences::const_iterator it = (*this).begin(); it != (*this).end(); it++)
      {
        vector<string> this_formats = split(it->m_file_format,"+");
        for (vector<string>::const_iterator it2 = this_formats.begin(); it2 != this_formats.end(); it2++) {
          formats.insert(*it2);
        }
      }
      
      vector<string> concat_formats;
      std::copy(formats.begin(), formats.end(), std::back_inserter(concat_formats));
      std::sort(concat_formats.begin(), concat_formats.end());
      return join(concat_formats, "+");
    }
    
    double get_total_gc_content() const
    {
      double gc = 0;
      uint64_t len = 0;
      for (cReferenceSequences::const_iterator it = (*this).begin(); it != (*this).end(); it++)
      {
        string s = it->m_fasta_sequence.get_sequence();
        len+= s.size();
        for (size_t i=0; i<s.length(); i++)
        {
          if ( (s[i] == 'G') ||(s[i] == 'C') )
            gc++;
        }
      }
      return gc / len;
    }

    void add_new_seq(const string& seq_id, const string& file_name)
    {
      m_seq_id_loaded[seq_id]++;
      if (m_seq_id_to_index.count(seq_id)) {
        return;
      } else {
        cAnnotatedSequence as;
        as.m_seq_id = seq_id;
        as.m_file_name = file_name;
        this->push_back(as);
        m_seq_id_to_index[seq_id] = m_index_id;
        m_index_id++;
      }
      return;
    }
    
    vector<string> seq_ids() const
    { 
      vector<string> ret_val;
      for(map<string,int>::const_iterator it=m_seq_id_to_index.begin(); it!= m_seq_id_to_index.end(); it++) {
        ret_val.push_back(it->first);
      }
      return ret_val;
    }
    
    //!< These gymnastics allow us to use [] to get a sequence by target_id (uint_32t) or by seq_id (string)
    
    void set_seq_id_to_index(const string& seq_id, int id)
      { m_seq_id_to_index[seq_id] = id; }
    
    uint32_t seq_id_to_index(const string& seq_id)
    { 
        ASSERT(m_seq_id_to_index.count(seq_id), "Reference seq id not found: " + seq_id + "\nValid seq ids: " + join(seq_ids(), ", ")); 
        return m_seq_id_to_index[seq_id]; 
    };

    cAnnotatedSequence& operator[](const size_t target_id)
      { return this->at(target_id); }
 
    const cAnnotatedSequence& operator[](const size_t target_id) const
    { return this->at(target_id); }
    
    cAnnotatedSequence& operator[](const string& seq_id)
      { ASSERT(m_seq_id_to_index.count(seq_id), "SEQ_ID not found: " + seq_id); return this->at(m_seq_id_to_index[seq_id]); }

    const cAnnotatedSequence& operator[](const string& seq_id) const
    { ASSERT(m_seq_id_to_index.count(seq_id), "SEQ_ID not found: " + seq_id); return this->at(m_seq_id_to_index.find(seq_id)->second); }

    bool seq_id_exists(const string& seq_id) const {
      return m_seq_id_to_index.count(seq_id);
    }
    
    bool is_circular(const string& seq_id) const
    {
      return (*this)[seq_id].is_circular();
    }
    
    //!< Utility to get sequences by seq_id
    string get_sequence_1(const string& seq_id, int32_t start_1, int32_t end_1)
    {
      return (*this)[seq_id].get_sequence_1(start_1, end_1);
    }
    
    string get_sequence_1(uint32_t tid, int32_t start_1, int32_t end_1) const
    {
      return (*this)[tid].get_sequence_1(start_1, end_1);
    }
    
    string get_sequence_1_start_size(const string& seq_id, const size_t start_1, const size_t size)
    {
      return (*this)[seq_id].get_sequence_1_start_size(start_1, size);
    }
    
    string get_sequence_1_start_size(uint32_t tid, const size_t start_1, const size_t size) const
    {
      return (*this)[tid].get_sequence_1_start_size(start_1, size);
    }

    void replace_sequence_1(const string& seq_id, int32_t start_1, int32_t end_1, const string& replacement_seq, string mut_type="", bool verbose=false)
    {
      (*this)[seq_id].replace_sequence_1(start_1, end_1, replacement_seq, mut_type, verbose);
    }

    void insert_sequence_1(const string& seq_id, int32_t pos, const string &insertion_seq, string mut_type="", bool verbose=false)
    {
      (*this)[seq_id].insert_sequence_1(pos, insertion_seq, mut_type, verbose);
    }
    
    void invert_sequence_1(const string& seq_id, int32_t start_1, int32_t end_1, string mut_type="", bool verbose=false)
    {
      (*this)[seq_id].invert_sequence_1(start_1, end_1, mut_type, verbose);
    }
    
    void repeat_feature_1(const string& seq_id, int32_t pos, int32_t start_del, int32_t end_del, cReferenceSequences& orig_ref_seq_info, string& orig_seq_id, int8_t strand, const cLocation&repeated_region, bool verbose = false)
    {
      (*this)[seq_id].repeat_feature_1(pos, start_del, end_del, orig_ref_seq_info, orig_seq_id, strand, repeated_region, verbose);
    }

    uint32_t get_sequence_length(const string& seq_id) const
    {
      return (*this)[seq_id].get_sequence_length();
    }

    vector<string> seq_ids()
    {
      vector<string> return_value;
      for(vector<cAnnotatedSequence>::iterator it=this->begin(); it != this->end(); it++)
        return_value.push_back(it->m_seq_id);
      return return_value;
    }
    
    // This is a duplicate of a function in pileup.h for when
    // we don't have a BAM available.
    void parse_region(const string& region, uint32_t& target_id, uint32_t& start_pos_1, uint32_t& end_pos_1)
    {      
      vector<string> split_region = split(region, ":");
      // Must check first split for second to not potentially crash
      ASSERT(split_region.size() == 2, "Unrecognized region: " + region + "\nExpected format [seq_id:start-end]");
      vector<string> split_positions = split(split_region[1], "-");
      ASSERT(split_positions.size() == 2, "Unrecognized region: " + region + "\nExpected format [seq_id:start-end]");

      target_id = seq_id_to_index(split_region[0]);
      start_pos_1 = from_string<uint32_t>(split_positions[0]);
      end_pos_1 = from_string<uint32_t>(split_positions[1]);
    }

    // This parses to a string for target
    static void parse_region(const string& region, string& seq_id, uint32_t& start_pos_1, uint32_t& end_pos_1)
    {
      vector<string> split_region = split(region, ":");
      // Must check first split for second to not potentially crash
      ASSERT(split_region.size() == 2, "Unrecognized region: " + region + "\nExpected format [seq_id:start-end]");
      vector<string> split_positions = split(split_region[1], "-");
      ASSERT(split_positions.size() == 2, "Unrecognized region: " + region + "\nExpected format [seq_id:start-end]");
      
      seq_id = split_region[0];
      start_pos_1 = from_string<uint32_t>(split_positions[0]);
      end_pos_1 = from_string<uint32_t>(split_positions[1]);
    }
    
    // normalize a region description
    // * remove commas from numbers
    // * flip start and end so start is smaller (returning bool whether this was done)
    static bool normalize_region(string& region) {
      region = substitute(region, ",", "");
      
      string seq_id;
      uint32_t start_pos_1;
      uint32_t end_pos_1;
      
      parse_region(region, seq_id, start_pos_1, end_pos_1);
      bool reverse = false;
      
      if (start_pos_1>end_pos_1) {
        reverse = true;
        std::swap(start_pos_1, end_pos_1);
      }
      
      region = seq_id + ":" + to_string<int32_t>(start_pos_1) + "-" + to_string<int32_t>(end_pos_1);
      return reverse;
    }
    
    map<string,int32_t> seq_order;
    map<string,string> trims;
    
    // Translation tables
    static vector<string> translation_tables;
    static vector<string> initiation_codon_translation_tables;
    static map<string,uint16_t> codon_to_aa_index;

    static cFeatureLocation* find_closest_repeat_region_boundary(int32_t position, const cSequenceFeatureList& repeat_list, int32_t& max_distance, int32_t direction, bool include_interior_matches = false);
    static cFeatureLocation* get_overlapping_feature(cFeatureLocationList& feature_list, int32_t pos);
    static char translate_codon(string seq, uint32_t translation_table, uint32_t codon_number_1, const string& gene="");
    static char translate_codon(string seq, string translation_table, string translation_table_1, uint32_t codon_number_1, const string& gene="");
    static string translate_protein(cAnnotatedSequence& seq, cSequenceFeature& loc, string translation_table, string translation_table_1);
    static void find_nearby_genes(
                                  cFeatureLocationList& gene_list,
                                  int32_t pos_1, 
                                  int32_t pos_2, 
                                  vector<cFeatureLocation*>& within_genes,
                                  vector<cFeatureLocation*>& between_genes,
                                  vector<cFeatureLocation*>& inside_left_genes,
                                  vector<cFeatureLocation*>& inside_right_genes,
                                  cFeatureLocation*& prev_gene,
                                  cFeatureLocation*& next_gene);
    
    bool mutation_overlapping_gene_is_inactivating(const cDiffEntry& mut, const string& snp_type, const uint32_t start, const uint32_t end, const cGeneFeature& gene, const double inactivate_overlap_fraction);
    
    static string list_to_entry(const vector<string>& _list, const string& _ignore);
    
    // Cleans out intergenic and list separators
    static string gene_strand_to_string(const bool forward)
    { return (forward ? gene_strand_forward_char : gene_strand_reverse_char); }
    
    void annotate_1_mutation_in_genes(cDiffEntry& mut, vector<cFeatureLocation*>& within_gene_locs, uint32_t start, uint32_t end, bool ignore_pseudogenes);
    void annotate_1_mutation(cDiffEntry& mut, uint32_t start, uint32_t end, bool repeat_override = false, bool ignore_pseudogenes = false);
    void categorize_1_mutation(cDiffEntry& mut, int32_t large_size_cutoff);
    void annotate_mutations(cGenomeDiff& gd, bool only_muts = false, bool ignore_pseudogenes = false, bool compare_mode = false, int32_t large_size_cutoff=kBreseq_large_mutation_size_cutoff, bool verbose = false);
    void polymorphism_statistics(Settings& settings, Summary& summary);
    string repeat_family_sequence(const string& repeat_name, int8_t strand, string* repeat_region = NULL, string* picked_seq_id=NULL, cFeatureLocation* picked_sequence_feature=NULL, bool fatal_error=true);
    
    string safe_seq_id_name(const string& input)
    {
      // return if not using safe seq ids
      if (!m_use_safe_seq_ids) return input;
      
      string s(input);
      
      size_t pos = s.find_first_not_of("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz1234567890._-");
      while(pos != string::npos) {
        //s.replace(pos, 1, ""); // remove
        s[pos] = '_'; pos++;// or replace with "safe" character
        pos = s.find_first_not_of("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz1234567890._-", pos);
      }
      
      // Remove double underscores which conflict with junction dividers
      string new_s = substitute(s,"__","_");
      while (new_s != s) {
        s = new_s;
        new_s = substitute(s,"__","_");
      }
      s = new_s;
      
      // Also remove any leading or trailing underscores
      pos = s.find_first_of("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz1234567890.-");
      s.replace(0, pos, "");

      pos = s.find_last_of("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz1234567890.-");
      s.replace(pos+1, s.size()-pos+1, "");
      
      ASSERT(s.find_first_of("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz1234567890.-") != string::npos, "Seq id [" + input + "] does not contain any valid alphanumeric characters.");
      
      if (s != input) {
        WARN("Reference seq id converted from '" + input + "' to '" + s + "'.\nOnly alphanumeric characters, periods, dashes, and single underscores '_' are allowed in seq ids.");
      }
      
      return s;
    }
    
    static string GFF3EscapeString(const string& s)
    {
      string escaped_s(s);
      
      escaped_s = substitute(escaped_s, "%", "%25");
      escaped_s = substitute(escaped_s, ";", "%3B");
      escaped_s = substitute(escaped_s, "=", "%3D");
      escaped_s = substitute(escaped_s, "&", "%26");
      escaped_s = substitute(escaped_s, ",", "%2C");
      return escaped_s;
    }
    
    //Substitute in reverse order as GFF3EscapeString, in case we introduced inconsistentancies.
    static string GFF3UnescapeString(const string& escaped_s)
    {
      string s(escaped_s);
      
      // Additional step for removing unnecessary outer quotes
      if ( (s[0] == '"') && (s[s.length()-1] == '"') ) {
        s = s.substr(1, s.length()-2);
      }
      
      s = substitute(s, "%2C", ",");
      s = substitute(s, "%26", "&");
      s = substitute(s, "%3D", "=");
      s = substitute(s, "%3B", ";");
      s = substitute(s, "%25", "%");
      return s;
    }    

    static cReferenceSequences deep_copy(const cReferenceSequences& in) {
      cReferenceSequences copy;

      copy.m_seq_id_to_index =  in.m_seq_id_to_index;
      copy.m_seq_id_loaded   =  in.m_seq_id_loaded; 
      copy.m_index_id        =  in.m_index_id;
      copy.m_initialized     =  in.m_initialized;

      copy.resize(in.size());
      for (uint32_t i = 0; i < in.size(); ++i) {
        copy[i] = cAnnotatedSequence::deep_copy(in[i]);
      }

      return copy;
    }


  };
    
  /*! Utility functions.
  */
  string GetWord(string &s);
  void RemoveLeadingWhitespace(string &s);
  void RemoveLeadingTrailingWhitespace(string &s);

  uint32_t alignment_mismatches(const alignment_wrapper& a, const cReferenceSequences& ref_seq_info);
  int32_t alignment_score(const alignment_wrapper& a, const cReferenceSequences& ref_seq_info);

  string shifted_cigar_string(const alignment_wrapper& a, const cReferenceSequences& ref_seq_info);
  struct sort_by_file_name : public binary_function<string, string, bool> {

    explicit sort_by_file_name(const string& name)
      : m_name(name) {}

     bool operator() (string a, string b) {
      if(a.find(m_name) != string::npos && b.find(m_name) == string::npos) {
        return true;
      }
      else if (a.find(m_name) != string::npos && b.find(m_name) != string::npos) {
        return a < b;
      }
       return false;
    }

    private:
      string m_name;
  };
  
  inline double gc_percentage_string(const string& seq)
  {
    double gc = 0;
    for (size_t i=0; i<seq.length(); i++)
    {
      if ( (seq[i] == 'G') ||(seq[i] == 'C') )
        gc++;
    }
    return gc / seq.length();
  }

} // breseq namespace

#endif
