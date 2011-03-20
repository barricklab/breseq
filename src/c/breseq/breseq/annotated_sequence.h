/*****************************************************************************

AUTHORS

  Jeffrey E. Barrick <jeffrey.e.barrick@gmail.com>
  David B. Knoester

LICENSE AND COPYRIGHT

  Copyright (c) 2010 Michigan State University

  breseq is free software; you can redistribute it and/or modify it under the  
  terms the GNU General Public License as published by the Free Software 
  Foundation; either version 1, or (at your option) any later version.

*****************************************************************************/

#ifndef _BRESEQ_ANNOTATED_SEQUENCE_H_
#define _BRESEQ_ANNOTATED_SEQUENCE_H_

#include <map>
#include <string>
#include <vector>

namespace breseq {
	
	/*! Interface for loading sequences and sequence features from GenBank files.
  */
  

	/*! Sequence Feature class.

	 */ 

  // Currently everything is stored as strings...
  typedef std::string sequence_feature_key_t;               //!< Diff entry keys.
  typedef std::string sequence_feature_value_t;            //!< Diff entry values.
  typedef std::map<sequence_feature_key_t, sequence_feature_value_t> sequence_feature_map_t; //!< Diff entry key-value map.
  
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
      
      //<! Safe accessor that returns empty string if not defined. 
      std::string SafeGet(sequence_feature_key_t in_key) { 
        sequence_feature_map_t::const_iterator it = this->find(in_key);
        if (it == this->end()) return std::string("");
        return it->second;
      }
      
      void ReadCoords(std::string& s);
      void ReadTag(std::string& tag, std::string& s, std::ifstream& in);
  };


	/*! Sequence class.
	 */   
   
  class cAnnotatedSequence {
    
    public:      
      uint32_t m_length;
      std::string m_definition, m_version, m_seq_id;
    
      std::string m_sequence;                             //!< Nucleotide sequence
      std::vector<cSequenceFeature> m_features;  //!< List of sequence features

    public:
      cAnnotatedSequence() {} ; //Constructor
            
      //!< Write a tab delimited feature 
      void WriteFeatureTable(const std::string &file_name);
      
      //!< Write FASTA file       
      void WriteFASTA(const std::string &file_name);

  };

  /*! Functions for creating cAnnotatedSequences
  */
  
  void LoadGenBankFile(const std::string &in_file, cAnnotatedSequence& s);
  
  /*! Utility functions.
  */
    
  std::string GetWord(std::string &s);
  void RemoveLeadingWhitespace(std::string &s);
  void RemoveLeadingTrailingWhitespace(std::string &s);
  
  void LoadGenBankFileHeader(std::ifstream& in, cAnnotatedSequence& s);
  void LoadGenBankFileSequenceFeatures(std::ifstream& in, cAnnotatedSequence& s);
  void LoadGenBankFileSequence(std::ifstream& in, cAnnotatedSequence& s);
	
} // breseq namespace

#endif
