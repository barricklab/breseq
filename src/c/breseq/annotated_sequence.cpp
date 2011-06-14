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

#include "breseq/annotated_sequence.h"

using namespace std;

namespace breseq {

  void cReferenceSequences::WriteFeatureTable(const std::string &file_name) {

    ofstream out(file_name.c_str());
    assert(!out.fail()); 

    for(vector<cAnnotatedSequence>::iterator it_as = this->begin(); it_as < this->end(); it_as++) {

      for (vector<cSequenceFeature>::iterator it = it_as->m_features.begin(); it < it_as->m_features.end(); it++) {
        out << it_as->m_seq_id;
        out << "\t" << (*it)["type"];
        out << "\t" << (*it)["accession"];
        out << "\t" << (*it)["name"];
        out << "\t" << (*it).m_start;
        out << "\t" << (*it).m_end;
        out << "\t" << (int)(*it).m_strand;
        out << "\t" << (*it)["product"];
        out << std::endl;
      }
      
    }
  }

  void cReferenceSequences::ReadFeatureTable(const string &file_name) {

    ifstream infile(file_name.c_str());
    assert(!infile.fail());
    
    string line;
    getline(infile,line);
    
    while (!infile.eof())
    {
      // split line on tabs
      char * cstr = new char [line.size()+1];
      strcpy (cstr, line.c_str());
      
      cSequenceFeature feature;
      string seq_id;
      
      char * pch;
      pch = strtok(cstr,"\t");
      uint32_t on_token = 0;
      while (pch != NULL)
      {
        //printf ("%s\n",pch);
        
        switch (on_token) {
          case 0:
            seq_id = pch;
            break;
          
          case 1:
            feature["type"] = pch;
            break;
            
          case 2:
            feature["accession"] = pch;
            break;
            
          case 3:
            feature["name"] = pch;
            break;
            
          case 4:
            feature.m_start = atoi(pch);
            break;
            
          case 5:
            feature.m_end = atoi(pch);
            break;
            
          case 6:
            feature.m_strand = atoi(pch);
            break;
            
          case 7:
            feature["product"] = pch;
            break;
        }
        
        pch = strtok (NULL, "\t");
        on_token++;
      }  
      
      delete[] cstr;
      
      uint32_t seq_idx = seq_id_to_index(seq_id);
      (*this)[seq_idx].m_features.push_back(feature);
            
      getline(infile,line);
    }

  }


  void cReferenceSequences::WriteFASTA(const std::string &file_name) {
    
    cFastaFile ff(file_name, ios_base::out);
    for(vector<cAnnotatedSequence>::iterator it_as = this->begin(); it_as < this->end(); it_as++) {
      ff.write_sequence(it_as->m_fasta_sequence);
    }
  }
  
  void cReferenceSequences::ReadFASTA(const std::string &file_name) {

    cFastaFile ff(file_name, ios_base::in);
    cFastaSequence on_seq;
    uint32_t on_seq_id = 0;
    
    while ( ff.read_sequence(on_seq) ) {
      
      cAnnotatedSequence new_seq;
      new_seq.m_fasta_sequence = on_seq;
      (*this).push_back(new_seq);
      
      m_seq_id_to_index[on_seq.m_name] = on_seq_id++;
    }
  }


  void LoadGenBankFile(cReferenceSequences& rs, const vector<string>& in_file_names) {
    
    for (vector<string>::const_iterator it = in_file_names.begin(); it < in_file_names.end(); it++) {
      
      ifstream in(it->c_str(), ios_base::in);
      assert(!in.fail()); 
      
      while (LoadGenBankFileHeader(in, rs)) {
        LoadGenBankFileSequenceFeatures(in, rs.back());
        LoadGenBankFileSequence(in, rs.back());
      }
    }
  }

  std::string GetWord(std::string &s) {
    RemoveLeadingWhitespace(s);
    int found = s.find_first_of(" =");
    std::string w = s.substr(0,found);
    s.erase(0, found);
    return w;
  }


  void RemoveLeadingWhitespace(std::string &s) {
    int found = s.find_first_not_of(" \t");
    s.erase(0,found);
  }

  void RemoveLeadingTrailingWhitespace(std::string &s) {
    int found = s.find_first_not_of(" \t\n\r");
    s.erase(0,found);  
    found = s.find_last_not_of(" \t\n\r");
    s.erase(found+1,s.length());
  }

  bool LoadGenBankFileHeader(std::ifstream& in, cReferenceSequences& rs) {		

    //std::cout << "header" << std::endl;
    std::string line;
    cAnnotatedSequence * s = NULL;
    bool found_LOCUS_line = false;
    while (!in.eof()) {
      std::getline(in, line);
      std::string first_word = GetWord(line);
      RemoveLeadingTrailingWhitespace(line);
      
      //std::cout << first_word << "::" << line << std::endl;

      // This is the first line -- if we see it we resize the list
      if (first_word == "LOCUS") {
        
        rs.resize(rs.size()+1);
        s = &(rs.back());
        
        std::string w;
        w = GetWord(line);
        s->m_seq_id = w;
        w = GetWord(line);
        s->m_length = atoi(w.c_str());
        assert (!found_LOCUS_line); // Should only be one line like this per record!
        found_LOCUS_line = true;
      }

      if (first_word == "DEFINITION") {
        s->m_definition = line;
      }
      
      if (first_word == "VERSION") {
        s->m_version = line;
      }
      
      if (first_word == "FEATURES") break;
    }
        
    return (found_LOCUS_line);
  }

  void cSequenceFeature::ReadCoords(std::string& s) {

    //std::cerr << "whole: " << s << std::endl;

    m_strand = 1;
    
    size_t start_complement = s.find("complement(");
    
    if (start_complement != string::npos) {
      //std::cerr << "before: " << s << std::endl;
      s.erase(0,11);
      
      size_t end_complement = s.find(")");
      assert(end_complement !=string::npos);
      s.erase(end_complement);
      
      //std::cerr << "after: " << s << std::endl;
      m_strand = -1;
    }
    
    int pos = s.find("..");
    std::string start_s = s.substr(0,pos);
    std::string end_s = s.substr(pos+2,s.length());

    m_start = atoi(start_s.c_str());
    m_end = atoi(end_s.c_str());
  }

  void cSequenceFeature::ReadTag(std::string& tag, std::string& s, std::ifstream& in) {

    // delete leading slash
    tag.erase(0,1);
      
    // erase through the equals on this line
    int pos = s.find("=");
    s.erase(0,pos+1);
    
    // If there is a quote then we need to continue to the last quote
    int first_quote_pos = s.find("\"");

    if (first_quote_pos == -1) {
      (*this)[tag] = s;
      return;
    } 
    
    s.erase(0,first_quote_pos+1);
    
    int second_quote_pos = s.find("\"");
    
    // One liner
    if (second_quote_pos != -1) {
      s.erase(second_quote_pos,s.length());
      (*this)[tag] = s;
      return;
    }

    // If the value is quoted, we have to read additional lines until end quote  
    std::string value = s;
    
    bool found_last_quote = false;
    while (!found_last_quote && !in.eof()) {
      std::getline(in, s);
      RemoveLeadingTrailingWhitespace(s);
      
      second_quote_pos = s.find("\"");
      if (second_quote_pos != -1) {
        s.erase(second_quote_pos,s.length());
        found_last_quote = true;
      }
      
      if (tag != "translation") value += " ";
      value += s;
    
    }
    assert(found_last_quote);
    
    (*this)[tag] = value;
  }

  void LoadGenBankFileSequenceFeatures(std::ifstream& in, cAnnotatedSequence& s) {
    //std::cout << "features" << std::endl;
    cSequenceFeature* current_feature=NULL;
    vector<cSequenceFeature> all_features;
    string line;
    while (!in.eof()) {
      getline(in, line);
      string first_word = GetWord(line);

      //std::cout << first_word << "::" << line << std::endl;
      
      // Done with this section...
      if (first_word == "ORIGIN") break;
      
      // Major tag = new feature or information block
      if (first_word[0] != '/') {
   
        if (first_word != "BASE") {

          all_features.resize(all_features.size()+1);
          current_feature = &(all_features[all_features.size()-1]);
          (*current_feature)["type"] = first_word;
          // parse the rest of the line
          std::string coord_s = GetWord(line);
          current_feature->ReadCoords(coord_s);
        }
      }
      // Minor tag = information about current feature
      else {
        assert(current_feature);
        
        // Remove leading slash
        current_feature->ReadTag(first_word, line, in); // reads multi-line entries
      }
    }
    
    
    for (vector<cSequenceFeature>::iterator it = all_features.begin(); it < all_features.end(); it++) {
    
      if ((*it)["type"] == "repeat_region") {

        // Don't add unnamed ones to the list...
        if (it->SafeGet("mobile_element") == "") continue;
                  
        (*it)["name"] = (*it)["mobile_element"];
        std::string& name = (*it)["name"];
        
        // remove prefix
        int pos = name.find("insertion sequence:");
        if (pos != -1) {
          name.erase(pos,pos+19);
        }
        
        // remove suffix if "IS(\d)}
        pos = name.find("IS");
        if (pos != -1) {
          int found = name.find_first_not_of("0123456789", 3);
          if (found != -1) { 
            name.erase(found, name.length());
          }
        }
        
        //std::cerr << (*it).SafeGet("mobile_element") << " " << (*it).SafeGet("name") << std::endl;
        
        (*it)["product"] = "repeat_region";
      
        s.m_features.push_back(*it);

      }
      else if ( ((*it)["type"] == "CDS") 
             || ((*it)["type"] == "tRNA") 
             || ((*it)["type"] == "rRNA") ) {
             
        // Add information
        if (it->SafeGet("gene") != "") {
          (*it)["name"] = (*it)["gene"];
        }
        if ( (it->SafeGet("name") == "") && (it->SafeGet("locus_tag") != "") ) {
          (*it)["name"] = (*it)["locus_tag"];
        }
        
        //std::cerr << (*it).SafeGet("name") << " " << (*it).SafeGet("gene") << " " << (*it).SafeGet("locus_tag") << std::endl;
              
        if ((*it).SafeGet("type") == "CDS") {
          (*it)["type"] = "protein";
        }
        
        (*it)["accession"] = (*it).SafeGet("locus_tag");
        
        if ((*it).SafeGet("pseudo") != "") {
          (*it)["type"] = "pseudogene";
        }
        
        (*it)["index"] = s.m_features.size();
        s.m_features.push_back(*it);
      } 
    }
  }

  void LoadGenBankFileSequence(std::ifstream& in, cAnnotatedSequence& s) {
    //std::cout << "sequence" << std::endl;

    s.m_fasta_sequence.m_name = s.m_seq_id;
    
    string line;
    while (!in.eof()) {
      getline(in, line);
      string first_word = GetWord(line);
      RemoveLeadingTrailingWhitespace(line);
      
      //cerr << first_word << "::" << line << std::endl;

      if (first_word == "//") break;
      for(string::iterator i=line.begin(); i!=line.end(); ++i) {
        if (*i == ' ') continue;
        *i = toupper(*i);
        
        // Scrub nonstandard characters
        if ((*i != 'A' ) && (*i != 'T' ) && (*i != 'C' ) && (*i != 'G' )) {
          *i = 'N';
        }
        
        s.m_fasta_sequence.m_sequence += *i;
      }
    }
    
    //cout << s.m_sequence << std::endl;
  }
  
  void LoadFeatureIndexedFastaFile(
                                   cReferenceSequences& s, 
                                   const string &in_feature_file_name, 
                                   const string &in_fasta_file_name
                                   ) 
  {
    s.ReadFASTA(in_fasta_file_name);
    if (in_feature_file_name.size() > 0)
    	s.ReadFeatureTable(in_feature_file_name);
  }

	uint32_t alignment_mismatches(alignment a, bam_header_t* header, faidx_t* fai, const cReferenceSequences* ref_seq_info)
	{
		bool verbose = false;
		uint32_t mismatches = 0;

		string seq_id = header->target_name[a.reference_target_id()];

		int32_t a_start, a_end;
		a.query_bounds_0(a_start, a_end);

		string ref_string;
		if (ref_seq_info != NULL)
		{
			cReferenceSequences ref_seq_info_value = *ref_seq_info;
			uint32_t index = ref_seq_info_value.seq_id_to_index(seq_id);
			ref_string = ref_seq_info_value[index].m_fasta_sequence.m_sequence;
			ref_string = ref_string.substr(a_start - 1, a_end - a_start + 1);
		}
		else
		{
			stringstream region_ss;
			region_ss << seq_id << ':' << a_start << '-' << a_end;
			string region = region_ss.str();
			int len = 0;
			ref_string = fai_fetch(fai, region.c_str(), &len);
		}

		vector<string> ref_string_list = split(ref_string, "/");
		uint32_t ref_pos = 0;

		string qseq = a.qseq();
		string read_string = qseq.substr(a_start - 1, a_end - a_start + 1);
		vector<string> read_string_list = split(read_string, "/");
		uint32_t read_pos = 0;

		uint32_t* cigar_list = a.cigar_array(); // cigar array for this alignment

		if (verbose)
		{
			cout << a.query_name() << endl;
			//cout << Dumper(cigar_list)
		}
		//#	my $cigar_string = '';

		for (uint32_t i = 0; i <= a.cigar_array_length(); i++)
		{
			char op = cigar_list[i] & BAM_CIGAR_MASK;
			uint32_t len = cigar_list[i] >> BAM_CIGAR_SHIFT;

			// soft padding counts as a mismatch
			if (op == 'S')
			{
				mismatches += len;
			}
			else if (op == 'D')
			{
				mismatches += len;
				ref_pos += len;
			}
			else if (op == 'I')
			{
				mismatches += len;
				read_pos += len;
			}
			else if (op == 'M')
			{
				for (uint32_t j = 0; j < len; j++)
				{
					if (verbose)
					{
						cout << "REF: " << ref_pos << "  " << ref_string_list[ref_pos] << endl;
						cout << "READ: " << read_pos << "  " << read_string_list[read_pos] << endl;
					}
					if (ref_string_list[ref_pos].compare(read_string_list[read_pos]) != 0);
						mismatches++;
					//#print "$read_pos $ref_pos\n";

					read_pos++;
					ref_pos++;
				}
			}
			//#$cigar_string .= $len . $op;
		}

		//#	print $a->qname . "\n$mismatches\n$cigar_string\n$ref_string\n$read_string\n" if ($mismatches);
		return mismatches;
	}

} // breseq namespace

