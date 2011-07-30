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

#include "breseq/common.h"

#include "breseq/fasta.h"
#include "breseq/alignment.h"

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
      cSequenceFeature(cSequenceFeature* _in) : sequence_feature_map_t(*_in) {
        m_start = _in->m_start;
        m_end = _in->m_end;
        m_strand = _in->m_strand;
      }
	  cSequenceFeature operator=(cSequenceFeature* _in) {
        m_start = _in->m_start;
        m_end = _in->m_end;
        m_strand = _in->m_strand;
        sequence_feature_map_t::operator=(*_in);
        return this;
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

  

	/*! Sequence class.
	 */   
   
  class cAnnotatedSequence {
    
    public:      
      uint32_t m_length;
      string m_definition, m_version, m_seq_id;
    
      cFastaSequence m_fasta_sequence;            //!< Nucleotide sequence
      vector<cSequenceFeature> m_features;  //!< List of sequence features

    public:
    
      //Constructor for empty object
      cAnnotatedSequence() : 
        m_length(0), 
        m_definition("na"), 
        m_version("na"), 
        m_seq_id("na"),
        m_features(0) {} ; 
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
      
    //!< Convert 
    uint32_t seq_id_to_index(const string& seq_id) { return m_seq_id_to_index[seq_id]; };

    map<string,int32_t> seq_order;
    map<string,string> trims;
    map<string,string> ref_strings;

    map<string, vector<cSequenceFeature> > repeat_lists;
    static cSequenceFeature* find_closest_repeat_region(uint32_t position, vector<cSequenceFeature>& repeat_list_ref, uint32_t max_distance, int32_t direction);
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
