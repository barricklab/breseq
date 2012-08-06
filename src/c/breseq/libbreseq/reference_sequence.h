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

#ifndef _BRESEQ_REFERENCE_SEQUENCE_H_
#define _BRESEQ_REFERENCE_SEQUENCE_H_

#include "libbreseq/common.h"

#include "libbreseq/genome_diff.h"
#include "libbreseq/fasta.h"
#include "libbreseq/fastq.h"
#include "libbreseq/alignment.h"
#include "libbreseq/anyoption.h"
#include "libbreseq/settings.h"
#include "libbreseq/summary.h"

namespace breseq {
	
	/*! Interface for loading sequences and sequence features from GenBank files.
  */
  class cLocation {
  public:
    int32_t m_start, m_end; // 1-indexed
    int8_t m_strand;
    vector<cLocation> m_sub_locations;
     
    cLocation() : m_start(0), m_end(0), m_strand(0) {};
    cLocation(int32_t start, int32_t end, int8_t strand) 
      : m_start(start), m_end(end), m_strand(strand) { 
    }

    cLocation(const cLocation& copy)
      : m_start(copy.m_start)
      , m_end(copy.m_end)
      , m_strand(copy.m_strand)
      , m_sub_locations(copy.m_sub_locations) { }


    cLocation& operator=(const cLocation& assign) {
      m_start         = assign.m_start;
      m_end           = assign.m_end;
      m_strand        = assign.m_strand;
      m_sub_locations = assign.m_sub_locations; 
      return *this;
    }

    //>! Get 1-indexed start position
    int32_t get_start_1() const {
      return m_start;
    }
    //>! Get 1-indexed end position
    int32_t get_end_1() const {
      return m_end;
    }
    int8_t get_strand() const {
      return m_strand;
    }

    vector<cLocation> get_all_sub_locations() {
       
      if (m_sub_locations.size() > 0) {
        vector<cLocation> return_locations;
        for (vector<cLocation>::iterator it=m_sub_locations.begin(); it!=m_sub_locations.end(); ++it) {
          vector<cLocation> this_locations = it->get_all_sub_locations();
          return_locations.insert(return_locations.end(), this_locations.begin(), this_locations.end());
        }
        return return_locations;
      }
      return make_vector<cLocation>(*this);
    }
    
    void set_start_1(int32_t _start) {
      m_start = _start;
    }
    void set_end_1(int32_t _end) {
      m_end = _end;
    }
    void set_strand(int8_t _strand) {
      m_strand = _strand;
    }

    void add_sub_location( const cLocation& value) {
      m_sub_locations.push_back(value);
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
      cLocation m_location;
      bool m_pseudo;
      map<string, vector<string> > m_gff_attributes;
    
      cSequenceFeature() : m_pseudo(0) {}

      cSequenceFeature(const cSequenceFeature& copy) 
        : sequence_feature_map_t(copy)
        , m_location(copy.m_location)
        , m_pseudo(copy.m_pseudo)
        , m_gff_attributes(copy.m_gff_attributes) { }

      bool operator<(const cSequenceFeature& _in) const
      {
        if (this->m_location.get_start_1() == _in.m_location.get_start_1())
          return (this->m_location.get_end_1() > _in.m_location.get_end_1());
        return (this->m_location.get_start_1() < _in.m_location.get_start_1());
      }
    
      //<! Safe accessor that returns empty string if not defined. 
      string SafeGet(sequence_feature_key_t in_key) 
      { 
        sequence_feature_map_t::const_iterator it = this->find(in_key);
        if (it == this->end()) return std::string("");
        return it->second;
      }
    
      //<! accessors
      int32_t get_start_1()     const { return m_location.get_start_1();  }
      int32_t get_end_1()       const { return m_location.get_end_1();    }
      int32_t get_strand()      const { return m_location.get_strand();   }
      bool is_top_strand()      const { return get_strand() == +1;        }
      bool is_bottom_strand()   const { return get_strand() == -1;        }
      string get_locus_tag() {
        if (m_gff_attributes.count("Alias")) {
          return m_gff_attributes["Alias"][0];
        }
        else if (this->count("locus_tag")) {
          return (*this)["locus_tag"];
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
      
      // Read GenBank coords
      void ReadGenBankCoords(string& s, ifstream& in);
      //Parse portion of GenBank coords string
      static cLocation ParseGenBankCoords(string& s, int8_t in_strand = 1);
    
      void ReadGenBankTag(string& tag, string& s, ifstream& in);
  };
  
  //!< Subclass of reference features with more information
  class Gene : public cSequenceFeature {
  public:
    string name;
    string product;
    string type;
    bool pseudogene; 
    uint32_t translation_table;
    
    Gene() {};
    Gene(cSequenceFeature& src) : cSequenceFeature(src)
    {
      name = src["name"];
      product = src["product"];
      type = src["type"];
      pseudogene = src.m_pseudo;
      translation_table = 1;
      if (src.count("transl_table")) 
        translation_table = from_string<uint32_t>(src["transl_table"]);
    }

    Gene(const Gene& copy) : cSequenceFeature(copy)
    {
      name = copy.name;
      product = copy.product;
      type = copy.type;
      pseudogene = copy.pseudogene;
      translation_table = copy.translation_table;
    }


  };

  typedef counted_ptr<Gene> cGenePtr;
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
    
      cFastaSequence m_fasta_sequence;            //!< Nucleotide sequence
    
      // Features are stored as counted pointers so that we can have ready-made lists
      // of different types of features. 
      cSequenceFeatureList m_features;    //!< Full list of sequence features
      cSequenceFeatureList m_genes;       //!< Subset of features
      cSequenceFeatureList m_repeats;     //!< Subset of features
    
    public:
    
      //Constructor for empty object
      cAnnotatedSequence() : 
        m_length(0),
        m_is_circular(false),
        m_description("na"), 
        m_seq_id("na"),
        m_features(0) {} ;

      // Utility to get top strand sequence
      string get_sequence_1(int32_t start_1, int32_t end_1) const
      {
        ASSERT(start_1 <= end_1, "start (" + to_string(start_1) + ") not less than or equal to end (" + to_string(end_1) + ")");
        return m_fasta_sequence.m_sequence.substr(start_1-1, end_1-start_1+1);
      }
    
      char get_sequence_1(int32_t pos_1) const
      {
        return m_fasta_sequence.m_sequence.substr(pos_1-1, 1)[0];
      }

      string get_circular_sequence_1(int32_t start_1, uint32_t size) const
      {
        const string &seq = m_fasta_sequence.m_sequence;

        if (start_1 > static_cast<int32_t>(seq.size())) {
          //If start_1 is too large, set to the beginning of the sequence.
          start_1 = start_1 % seq.size();
        } 
        else
        if (start_1 < 0) {
          //If start_1 is negative, set to the end of the sequence.
          start_1 = seq.size() + start_1; 
        }

        string ret_val = "";
        try {
          ret_val = seq.substr(start_1, size);

          if (size > ret_val.size()) {
            ret_val.append(seq.substr(0, size - ret_val.size()));
          }

        } catch(...) {
          ERROR("[start_1]: " + to_string(start_1) +
                "[size]:    " + to_string(size));
        }

        return ret_val;
      }

      size_t get_sequence_size(void) const
      {
        return m_fasta_sequence.m_sequence.size();
      }

      // Replace Sequence with Input
      void replace_sequence_1(int32_t start_1, int32_t end_1, const string &replacement_seq, string mut_type="", bool verbose=false);
      
      // Insert Input AFTER Position
      void insert_sequence_1(int32_t pos_1, const string &insertion_seq, string mut_type="", bool verbose=false);
    
      // Repeat Feature at Position
      void repeat_feature_1(int32_t pos, int32_t start_del, int32_t end_del, cReferenceSequences& ref_seq_info, const string &repeat_name, int8_t strand, int32_t region_pos, bool verbose=false);
      
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
        
      uint32_t get_sequence_length()
      {
        return m_fasta_sequence.m_sequence.length();
      }
    
      //! Correctly adds features across different lists
      void feature_push_back(cSequenceFeaturePtr& fp)
      {
        m_features.push_back(fp);
        
        if ( ((*fp)["type"] == "repeat_region") || ((*fp)["type"] == "mobile_element") )
        {
          m_repeats.push_back(fp);
        }
        else if ( ((*fp)["type"] == "CDS") || ((*fp)["type"] == "tRNA") || ((*fp)["type"] == "rRNA") || ((*fp)["type"] == "RNA") )
        { 
          m_genes.push_back(fp);
        }
      }
    
      //! Correctly adds features across different lists
      void feature_push_front(cSequenceFeaturePtr& fp)
      {
        m_features.push_front(fp);
        
        if ((*fp)["type"] == "repeat_region")
        {
          m_repeats.push_front(fp);
        }
        else if (((*fp)["type"] == "CDS") || ((*fp)["type"] == "tRNA") || ((*fp)["type"] == "rRNA") || ((*fp)["type"] == "RNA"))
        {
          m_genes.push_front(fp);
        }
      }
      
      static cAnnotatedSequence deep_copy(cAnnotatedSequence in) {
        cAnnotatedSequence copy;


        copy.m_length = in.m_length;
        copy.m_is_circular = in.m_is_circular;
        copy.m_description = in.m_description;
        copy.m_seq_id = in.m_seq_id;
        copy.m_fasta_sequence = in.m_fasta_sequence;

        cSequenceFeatureList::iterator it;
        //Features.
        for (it = in.m_features.begin(); it != in.m_features.end(); ++it) {
          cSequenceFeature* temp = new cSequenceFeature(**it);
          copy.m_features.push_back(cSequenceFeaturePtr(temp));
        }

        //Genes.
        for (it = in.m_genes.begin(); it != in.m_genes.end(); ++it) {
          cSequenceFeature* temp = new cSequenceFeature(**it);
          copy.m_genes.push_back(cSequenceFeaturePtr(temp));
        }

        //Repeats.
        for (it = in.m_repeats.begin(); it != in.m_repeats.end(); ++it) {
          cSequenceFeature* temp = new cSequenceFeature(**it);
          copy.m_repeats.push_back(cSequenceFeaturePtr(temp));
        }

        return copy;
      }
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

    //!< Currently supported file types.
    enum FileType {UNKNOWN, GENBANK, FASTA, GFF3, BULL};


  public:
    
    cReferenceSequences()
      : m_index_id(0)
      , m_initialized(false)
    {}

    //!< Load all reference files and verify
    void LoadFiles(const vector<string>& file_names);
    
    //!< Load reference file into object
    //!< Detect file type and load the information in it appropriately
    void LoadFile(const string& file_name);
    
    //!< Verify that all seq_id have sequence;
    void Verify();
    bool Initialized() {return m_initialized;}
    
    //!< Read/Write FASTA file     
    void ReadFASTA(const std::string &file_name);
    void ReadFASTA(cFastaFile& ff);
    void WriteFASTA(const string &file_name, bool verbose=false);
    void WriteFASTA(cFastaFile& ff, bool verbose=false);
        
    //!< Read/Write a tab delimited GFF3 file
    void ReadGFF(const string& file_name);
    void WriteGFF(const string &file_name, bool verbose=false);

    //!< Read GenBank file
    void ReadGenBank(const string& in_file_names);
    bool ReadGenBankFileHeader(std::ifstream& in);
    void ReadGenBankCoords(string& s, ifstream& in);
    //void ReadGenBankTag(std::string& tag, std::string& s, std::ifstream& in);
    void ReadGenBankFileSequenceFeatures(std::ifstream& in, cAnnotatedSequence& s);
    void ReadGenBankFileSequence(std::ifstream& in, cAnnotatedSequence& s);
    
    //!< Read Bull gene table file
    void ReadBull(const string& file_name);
    
    //!< Calculates the total length of all reference sequences together
    uint32_t total_length()
    {
      uint32_t ret_val(0);
      for (cReferenceSequences::iterator it = (*this).begin(); it != (*this).end(); it++)
      {
        ret_val += it->m_length;
      }
      return ret_val;
    }

    void add_new_seq(const string& seq_id)
    {
      m_seq_id_loaded[seq_id]++;
      if (m_seq_id_to_index.count(seq_id)) {
        return;
      } else {
        cAnnotatedSequence as;
        as.m_seq_id = seq_id;
        this->push_back(as);
        m_seq_id_to_index[seq_id] = m_index_id;
        m_index_id++;
      }
      return;
    }
    
    //!< These gymnastics allow us to use [] to get a sequence by target_id (uint_32t) or by seq_id (string)
    
    void set_seq_id_to_index(const string& seq_id, int id)
      { m_seq_id_to_index[seq_id] = id; }
    
    uint32_t seq_id_to_index(const string& seq_id)
      { ASSERT(m_seq_id_to_index.count(seq_id), "SEQ_ID not found: " + seq_id); return m_seq_id_to_index[seq_id]; };

    cAnnotatedSequence& operator[](const size_t target_id)
      { return this->at(target_id); }
 
    const cAnnotatedSequence& operator[](const size_t target_id) const
    { return this->at(target_id); }
    
    cAnnotatedSequence& operator[](const string& seq_id)
      { ASSERT(m_seq_id_to_index.count(seq_id), "SEQ_ID not found: " + seq_id); return this->at(m_seq_id_to_index[seq_id]); }

    
    //!< Utility to get sequences by seq_id
    string get_sequence_1(const string& seq_id, int32_t start_1, int32_t end_1)
    {
      return (*this)[seq_id].get_sequence_1(start_1, end_1);
    }
    
    string get_sequence_1(uint32_t tid, int32_t start_1, int32_t end_1) const
    {
      // TODO: Handle circular genomes
      return (*this)[tid].get_sequence_1(start_1, end_1);
    }
    
    string get_circular_sequence_1(const string& seq_id, const size_t start_1, const size_t size)
    {
      return (*this)[seq_id].get_circular_sequence_1(start_1, size);
    }
    
    string get_circular_sequence_1(uint32_t tid, const size_t start_1, const size_t size) const
    {
      return (*this)[tid].get_circular_sequence_1(start_1, size);
    }

    void replace_sequence_1(const string& seq_id, int32_t start_1, int32_t end_1, const string& replacement_seq, string mut_type="", bool verbose=false)
    {
      (*this)[seq_id].replace_sequence_1(start_1, end_1, replacement_seq, mut_type, verbose);
    }

    void insert_sequence_1(const string& seq_id, int32_t pos, const string &insertion_seq, string mut_type="", bool verbose=false)
    {
      (*this)[seq_id].insert_sequence_1(pos, insertion_seq, mut_type, verbose);
    }
    
    void repeat_feature_1(const string& seq_id, int32_t pos, int32_t start_del, int32_t end_del, cReferenceSequences& ref_seq_info, const string &repeat_name, int8_t strand, int32_t region_pos, bool verbose=false)
    {
      (*this)[seq_id].repeat_feature_1(pos, start_del, end_del, ref_seq_info, repeat_name, strand, region_pos, verbose);
    }

    uint32_t get_sequence_length(const string& seq_id)
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
      vector<string> split_positions = split(split_region[1], "-");
      ASSERT(split_region.size()+split_positions.size()-1 == 3, "Unrecognized region: " + region + "\n(Expected seq_id:start-end)");

      target_id = seq_id_to_index(split_region[0]);
      start_pos_1 = from_string<uint32_t>(split_positions[0]);
      end_pos_1 = from_string<uint32_t>(split_positions[1]);
    }


    map<string,int32_t> seq_order;
    map<string,string> trims;
    
    // Translation tables
    static vector<string> translation_tables;
    static vector<string> initiation_codon_translation_tables;
    static map<string,uint16_t> codon_to_aa_index;

    static cSequenceFeaturePtr find_closest_repeat_region_boundary(int32_t position, cSequenceFeatureList& repeat_list, int32_t& max_distance, int32_t direction);
    static cSequenceFeaturePtr get_overlapping_feature(cSequenceFeatureList& feature_list, int32_t pos);
    static char translate_codon(string seq, uint32_t translation_table, uint32_t codon_pos_1);
    static void find_nearby_genes(
                                  cSequenceFeatureList& gene_list, 
                                  int32_t pos_1, 
                                  int32_t pos_2, 
                                  vector<Gene>& within_genes, 
                                  vector<Gene>& between_genes, 
                                  vector<Gene>& inside_left_genes, 
                                  vector<Gene>& inside_right_genes,
                                  Gene& prev_gene, 
                                  Gene& next_gene);
    void annotate_1_mutation(cDiffEntry& mut, uint32_t start, uint32_t end, bool repeat_override = false, bool ignore_pseudogenes = false);
    void annotate_mutations(cGenomeDiff& gd, bool only_muts = false, bool ignore_pseudogenes = false);
    void polymorphism_statistics(Settings& settings, Summary& summary);
    string repeat_family_sequence(const string& repeat_name, int8_t strand, int32_t& region_pos, string* seq_id=NULL, cSequenceFeature* sequence_feature=NULL);
    
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

} // breseq namespace

#endif
