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

#ifndef _BRESEQ_ANNOTATED_SEQUENCE_H_
#define _BRESEQ_ANNOTATED_SEQUENCE_H_

#include "libbreseq/common.h"

#include "libbreseq/fasta.h"
#include "libbreseq/fastq.h"
#include "libbreseq/alignment.h"
#include "libbreseq/genome_diff.h"

using namespace std;

namespace breseq {
	
	/*! Interface for loading sequences and sequence features from GenBank files.
  */
  
  
	/*! Sequence Feature class.

	 */ 

  // Currently everything is stored as strings...
  typedef string sequence_feature_key_t;
  typedef string sequence_feature_value_t;
  typedef map<sequence_feature_key_t, sequence_feature_value_t> sequence_feature_map_t; //!< Diff entry key-value map.
  
  class cSequenceFeature : public sequence_feature_map_t {
    
    public:
    
      // Could add accessors that convert strings to numbers...
      uint32_t m_start, m_end;
      int8_t m_strand;
    
      cSequenceFeature() {}
      cSequenceFeature(const cSequenceFeature& _in) : sequence_feature_map_t(_in) {
        m_start = _in.m_start;
        m_end = _in.m_end;
        m_strand = _in.m_strand;
      }
	  cSequenceFeature operator=(const cSequenceFeature& _in) {
        m_start = _in.m_start;
        m_end = _in.m_end;
        m_strand = _in.m_strand;
        sequence_feature_map_t::operator=(_in);
        return *this;
      }
      
      //<! Safe accessor that returns empty string if not defined. 
      std::string SafeGet(sequence_feature_key_t in_key) { 
        sequence_feature_map_t::const_iterator it = this->find(in_key);
        if (it == this->end()) return std::string("");
        return it->second;
      }
      
      void ReadCoords(string& s, ifstream& in);
      void ReadTag(string& tag, string& s, ifstream& in);
  };
  
  //!< Subclass of reference features with more information
  class Gene : public cSequenceFeature {
  public:
    string name;
    string product;
    string type;
    uint32_t start;
    uint32_t end;
    bool strand;
    bool pseudogene; 
    
    Gene() {};
    Gene(cSequenceFeature& src)
    {
      name = src["name"];
      product = src["product"];
      type = src["type"];
      start = src.m_start;
      end = src.m_end;
      strand = (src.m_strand >= 1);
      pseudogene = false;
    }
  };

  typedef counted_ptr<Gene> cGenePtr;
  typedef counted_ptr<cSequenceFeature> cSequenceFeaturePtr;
  typedef vector<cSequenceFeaturePtr> cSequenceFeatureList;
  
	/*! Sequence class.
	 */   
   
  class cAnnotatedSequence {
    
    public:      
      uint32_t m_length;
      string m_definition, m_version, m_seq_id;
    
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
        m_definition("na"), 
        m_version("na"), 
        m_seq_id("na"),
        m_features(0) {} ;
    
      // Utility to get top strand sequence
      string get_sequence_1(uint32_t start_1, uint32_t end_1) 
      {
        ASSERTM(start_1 <= end_1, "start (" + to_string(start_1) + ")not less than or equal to end(" + to_string(end_1) + ")");
        return m_fasta_sequence.m_sequence.substr(start_1-1, end_1-start_1+1);
      }

      void replace_sequence_1(uint32_t start_1, uint32_t end_1, const string &replacement_seq)
      {
        ASSERTM(start_1 <= end_1, "start (" + to_string(start_1) + ")not less than or equal to end(" + to_string(end_1) + ")");
        m_fasta_sequence.m_sequence.replace(start_1-1, end_1-start_1+1, replacement_seq);
      }

      // Inserts AFTER the input position
      void insert_sequence_1(uint32_t pos_1, const string &insertion_seq)
      {
        m_fasta_sequence.m_sequence.insert(pos_1, insertion_seq);
      }

      uint32_t get_sequence_length()
      {
        return m_fasta_sequence.m_sequence.length();
      }
    
      // Correctly adds features across different lists
      void add_feature(cSequenceFeaturePtr& fp)
      {
        m_features.push_back(fp);
        
        if ((*fp)["type"] == "repeat_region")
        {
          m_repeats.push_back(fp);
        }
        else
        {
          m_genes.push_back(fp);
        }
      }
  };

  
  /*! Reference Sequences
   
   Holds sequences and features for ALL reference sequences.
	 */ 
  
  class cReferenceSequences : public vector<cAnnotatedSequence> {
  protected:
    map<string,int> m_seq_id_to_index; // for looking up sequences by seq_id
    
  public:
    
    cReferenceSequences() {};    
    
    //!< Write a tab delimited feature 
    void WriteFeatureTable(const string &file_name);
    
    //!< Read a tab delimited feature
    void ReadFeatureTable(const string &file_name);
    
    //!< Write FASTA file       
    void WriteFASTA(const string &file_name);
    
    //!< Read FASTA file       
    void ReadFASTA(const std::string &file_name);
    
    //!< Write a tab delimited GFF3 file
    void WriteGFF( const string &file_name );  
      
    
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
    
    //!< These gymnastics allow us to use [] to get a sequence by target_id (uint_32t) or by seq_id (string)
    
    void set_seq_id_to_index(const string& seq_id, int id)
      { m_seq_id_to_index[seq_id] = id; }
    
    uint32_t seq_id_to_index(const string& seq_id)
      { ASSERT(m_seq_id_to_index.count(seq_id)); return m_seq_id_to_index[seq_id]; };

    cAnnotatedSequence& operator[](const size_t target_id)
      { return this->at(target_id); }
 
    const cAnnotatedSequence& operator[](const size_t target_id) const
    { return this->at(target_id); }
    
    cAnnotatedSequence& operator[](const string& seq_id)
      { ASSERT(m_seq_id_to_index.count(seq_id)); return this->at(m_seq_id_to_index[seq_id]); }

    
    //!< Utility to get sequences by seq_id
    string get_sequence_1(const string& seq_id, uint32_t start_1, uint32_t end_1) 
    {
      return (*this)[seq_id].get_sequence_1(start_1, end_1);
    }

    void replace_sequence_1(const string& seq_id, uint32_t start_1, uint32_t end_1, const string& replacement_seq)
    {
      (*this)[seq_id].replace_sequence_1(start_1, end_1, replacement_seq);
    }

    void insert_sequence_1(const string& seq_id, uint32_t pos, const string &insertion_seq)
    {
      (*this)[seq_id].insert_sequence_1(pos, insertion_seq);
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
      
      vector<string> split_region = split_on_any(region, ":-");
      ASSERTM(split_region.size() == 3, "Unrecognized region: " + region + "\n(Expected seq_id:start-end)");

      target_id = seq_id_to_index(split_region[0]);
      start_pos_1 = from_string<uint32_t>(split_region[1]);
      end_pos_1 = from_string<uint32_t>(split_region[2]);
    }


    map<string,int32_t> seq_order;
    map<string,string> trims;
    
    static map<string,char> translation_table_11;

    static cSequenceFeature* find_closest_repeat_region(uint32_t position, vector<cSequenceFeaturePtr>& repeat_list_ref, uint32_t max_distance, int32_t direction);
    static cSequenceFeature* get_overlapping_feature(vector<cSequenceFeaturePtr>& feature_list_ref, uint32_t pos);
    static char translate(string seq);
    static void find_nearby_genes(
                                  cSequenceFeatureList& gene_list_ref, 
                                  uint32_t pos_1, 
                                  uint32_t pos_2, 
                                  vector<Gene>& within_genes, 
                                  vector<Gene>& between_genes, 
                                  vector<Gene>& inside_left_genes, 
                                  vector<Gene>& inside_right_genes,
                                  Gene& prev_gene, 
                                  Gene& next_gene);
    void annotate_1_mutation(diff_entry& mut, uint32_t start, uint32_t end, bool repeat_override = false);
    void annotate_mutations(genome_diff& gd, bool only_muts = false);
    void polymorphism_statistics(Settings& settings, Summary& summary);
    string repeat_example(const string& repeat_name, int8_t strand);

    
  };
  
  /*! Helper function for creating cReferenceSequences
   */
  
  void LoadGenBankFile(cReferenceSequences& s, const vector<string>& in_file_names);
  bool LoadGenBankFileHeader(ifstream& in, cReferenceSequences& s);
  void LoadGenBankFileSequenceFeatures(ifstream& in, cAnnotatedSequence& s);
  void LoadGenBankFileSequence(ifstream& in, cAnnotatedSequence& s);
  
  void LoadFeatureIndexedFastaFile(cReferenceSequences& s, const string &in_feature_file_name, const string &in_fasta_file_name);
  
  void LoadBullFile(cReferenceSequences& s, const vector<string>& in_file_names);
  void LoadBullFeatureFile(ifstream& in, cAnnotatedSequence& s);

  /*! Utility functions.
  */
    
  std::string GetWord(string &s);
  void RemoveLeadingWhitespace(string &s);
  void RemoveLeadingTrailingWhitespace(string &s);

  uint32_t alignment_mismatches(const alignment_wrapper& a, const cReferenceSequences& ref_seq_info);
  string shifted_cigar_string(const alignment_wrapper& a, const cReferenceSequences& ref_seq_info);

} // breseq namespace

#endif
