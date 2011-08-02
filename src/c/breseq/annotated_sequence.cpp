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
    
    out << "##gff-version\t3" << endl;
    for(vector<cAnnotatedSequence>::iterator it_as = this->begin(); it_as < this->end(); it_as++) {
    
      
      uint32_t m_length;
      string m_definition, m_version, m_seq_id;
      
      out << "##sequence-region " << it_as->m_seq_id << " 1 " << it_as->m_length << endl;
      out << "##description " << it_as->m_seq_id << " " << it_as->m_definition << endl;

    }
    
    for(vector<cAnnotatedSequence>::iterator it_as = this->begin(); it_as < this->end(); it_as++) {

      for (vector<cSequenceFeature>::iterator it = it_as->m_features.begin(); it < it_as->m_features.end(); it++) {
        out << it_as->m_seq_id;
        
        if (it->SafeGet("type") == "") 
          out << "\t" << ".";
        else
          out << "\t" << (*it)["type"];
        
        if( it->SafeGet("accession") == "" )
          out << "\t" << ".";
        else
          out << "\t" <<(*it)["accession"];
        
        if( it->SafeGet("name") == "" )
          out << "\t" << ".";
        else
          out << "\t" <<(*it)["name"];
        
        out << "\t" << (*it).m_start;
        out << "\t" << (*it).m_end;
        out << "\t" << (int)(*it).m_strand;
        
        if( it->SafeGet("product") == "" )
          out << "\t" << ".";
        else
          out << "\t" <<(*it)["product"];
        out << std::endl;
      }
      
    }
  }
  
  void cReferenceSequences::WriteGFF( const string &file_name ){
      
      ofstream out(file_name.c_str());
      assert(!out.fail()); 
      
      for(vector<cAnnotatedSequence>::iterator it_as = this->begin(); it_as < this->end(); it_as++) {
          
          for (vector<cSequenceFeature>::iterator it = it_as->m_features.begin(); it < it_as->m_features.end(); it++) {
            out << it_as->m_seq_id;
            out << "\tGenBank";
          
            if (it->SafeGet("type") == "") 
              out << ".";
            else
              out << "\t" << (*it)["type"];
            
            
            out << "\t" << (*it).m_start;
            out << "\t" << (*it).m_end;
            out << "\t.";
            out << "\t" << (int)(*it).m_strand;
            out << "\t.";
            
            if( it->SafeGet("accession") == "" )
              out << ".";
            else
              out << "\t" << "accession=" <<(*it)["accession"];
            
            if( it->SafeGet("name") == "" )
              out << ";" << "name=" << ".";
            else
              out << ";" << "name=" <<(*it)["name"];
            
            if( it->SafeGet("product") == "" )
              out << ";" << "produce=" << ".";
            else
              out << ";" << "produce=" <<(*it)["product"];
            
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
      if (line[0] != '#') {
        
      
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
        
        if (feature["type"] == "repeat_region")
          (*this).repeat_lists[seq_id].push_back(feature);
      }
      getline(infile,line);
    }

  }

  cSequenceFeature* cReferenceSequences::find_closest_repeat_region(uint32_t position, vector<cSequenceFeature>& repeat_list_ref, uint32_t max_distance, int32_t direction)
  {
	if (repeat_list_ref.size() == 0) return NULL;

	cSequenceFeature* is = NULL;
	int32_t best_distance = 0;

	for (uint32_t i = 0; i < repeat_list_ref.size(); i++) //IS
	{
		cSequenceFeature* test_is = &(repeat_list_ref[i]);

		//count within the IS element as zero distance
		//if this happens then we are immediately done
		if ( (test_is->m_start <= position) && (test_is->m_end >= position) )
			return test_is;

		//otherwise calculate the distance
		//keep if less than max_distance, in the correct direction, and the best found so far
		int32_t test_distance = ((direction == -1) ? position - test_is->m_end : test_is->m_start - position);
		if (test_distance < 0) continue; //wrong direction...

		if ((test_distance <= (int32_t)max_distance) && ((is == NULL) || (test_distance < best_distance)) )
		{
			is = test_is;
			best_distance = test_distance;
		}
	}
	return is;
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
      new_seq.m_seq_id = on_seq.m_name;
      new_seq.m_length = on_seq.m_sequence.size();
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

  void cSequenceFeature::ReadCoords(string& s, ifstream& in) {

    //std::cerr << "whole: " << s << std::endl;

    // Read through all parentheses
    int32_t parentheses_level = 0;
    size_t parenthesis_pos = s.find_first_of("()");
    
    while(parenthesis_pos != string::npos) {
      if (s.at(parenthesis_pos) == '(') {
        parentheses_level++;
      } else {
        parentheses_level--;
      }
      parenthesis_pos = s.find_first_of("()", parenthesis_pos+1);
    }
    
    string value = s;
    
    // Multiline
    while ((parentheses_level > 0) && !in.eof()) {
      std::getline(in, s);
      RemoveLeadingTrailingWhitespace(s);
      
      size_t on_pos = s.find_first_of("()");
      while(on_pos != string::npos) {
        
        if (s.at(on_pos) == '(') {
          parentheses_level++;
        } else {
          parentheses_level--;
        }
        on_pos = s.find_first_of("()", on_pos+1);
      }
      value += s;
    }
    assert(parentheses_level == 0);    
    
    m_strand = 1;
    size_t start_complement = value.find("complement(");
    
    if (start_complement != string::npos) {
      //std::cerr << "before: " << s << std::endl;
      //value.erase(0,11);
      
      //size_t end_complement = value.find(")");
      //assert(end_complement !=string::npos);
      //value.erase(end_complement);
      
      //std::cerr << "after: " << s << std::endl;
      m_strand = -1;
    }
    
    size_t start_start_pos = value.find_first_of("1234567890");
    assert(start_start_pos != string::npos);
    size_t start_end_pos = value.find_first_not_of("1234567890", start_start_pos);
    if (start_end_pos == string::npos) start_end_pos = value.size() ;
    string start_s = value.substr(start_start_pos,start_end_pos-start_start_pos);
    
    size_t end_end_pos = value.find_last_of("1234567890");
    assert(end_end_pos != string::npos);
    size_t end_start_pos = value.find_last_not_of("1234567890", end_end_pos);
    if (end_start_pos == string::npos) start_end_pos = -1;
    
    string end_s = value.substr(end_start_pos+1,end_end_pos+1 - (end_start_pos+1));

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
    size_t first_quote_pos = s.find("\"");
    
    if (first_quote_pos == string::npos) {
      (*this)[tag] = s;
      return;
    }
    
    s.erase(0,first_quote_pos+1);
    
    size_t second_quote_pos = s.find("\"");
    
    // One liner
    if (second_quote_pos != string::npos) {
      s.erase(second_quote_pos,s.length());
      (*this)[tag] = s;
      return;
    }

    // If the value is quoted, we have to read additional lines until end quote  
    string value = s;
    
    bool found_last_quote = false;
    while (!found_last_quote && !in.eof()) {
      std::getline(in, s);
      RemoveLeadingTrailingWhitespace(s);
      
      second_quote_pos = s.find("\"");
      if (second_quote_pos != string::npos) {
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
      //debug
      //cout << line << endl;
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
          current_feature->ReadCoords(coord_s, in);
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
        
        // /pseudo tag doesn't take a value
        if ((*it).count("pseudo") != 0) {
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
  
  void LoadBullFile(cReferenceSequences& rs, const vector<string>& in_file_names) {
    
    for (vector<string>::const_iterator it = in_file_names.begin(); it < in_file_names.end(); it++) {
      
      ifstream in(it->c_str(), ios_base::in);
      assert(!in.fail()); 
      
      LoadBullFeatureFile(in, rs.back());
    }
    
   
    
    /*for (vector<cSequenceFeature>::iterator it = s.m_features.begin(); it < s.m_features.end(); it++) {
      
      if (it->SafeGet("type") == "") 
        cout << "\t.";
      else
        cout << "\t" << (*it)["type"];
      
      
      cout << "\t" << (*it).m_start;
      cout << "\t" << (*it).m_end;
      cout << "\t.";
      cout << "\t" << (int)(*it).m_strand;
      cout << "\t.";
      
      if( it->SafeGet("accession") == "" )
        cout << "\t.";
      else
        cout << "\t" << "accession=" <<(*it)["accession"];
      
      if( it->SafeGet("name") == "" )
        cout << "\t" << ".";
      else
        cout << "\t" << (*it)["name"];
      
      if( it->SafeGet("product") == "" )
        cout << "\t" << ".";
      else
        cout << "\t" <<(*it)["product"];
      
      cout << std::endl;
    }*/
    
  }
  
  void LoadBullFeatureFile(ifstream& in, cAnnotatedSequence& s) {
      
    char line[10];
    
    vector<cSequenceFeature> all_features;
    
    while ( !in.eof() ) {
      
      cSequenceFeature* current_feature(NULL);
      
      all_features.resize(all_features.size()+1);
      current_feature = &(all_features[all_features.size()-1]);
      
      current_feature->m_strand = 1;
      
      in.getline(line, 10, ' ');
      current_feature->m_start = from_string<uint32_t>(line);
      
      in.getline(line, 10, ' ');
      current_feature->m_end = from_string<uint32_t>(line);
      
      in.getline(line, 10, '\n');
      string name( (string) line );
      name.erase(std::remove(name.begin(), name.end(), '\n'), name.end());
      (*current_feature)["name"] = name;
    }
    
    for (vector<cSequenceFeature>::iterator it = all_features.begin(); it < all_features.end(); it++) {
      //s.m_features.push_back(*it);
      cout << "Start: " << (*it).m_start << " Stop: " << (*it).m_end << " Strand: " << (*it).m_strand << endl;
    }
  }

	uint32_t alignment_mismatches(const alignment_wrapper& a, const cReferenceSequences& ref_seq_info)
	{
		bool verbose = false;
		uint32_t mismatches = 0;

    uint32_t index = a.reference_target_id();
    const string& const_ref_string = ref_seq_info[index].m_fasta_sequence.m_sequence;
    string ref_string = const_ref_string.substr(a.reference_start_0(), a.reference_match_length());
		uint32_t ref_pos = 0;

		string read_string = a.read_char_sequence().substr(a.query_start_0(), a.query_match_length());
		uint32_t read_pos = 0;

		uint32_t* cigar_list = a.cigar_array(); // cigar array for this alignment

		if (verbose)
		{
			cout << a.read_name() << endl;
		}

		for (uint32_t i = 0; i < a.cigar_array_length(); i++)
		{
			char op = cigar_list[i] & BAM_CIGAR_MASK;
			uint32_t len = cigar_list[i] >> BAM_CIGAR_SHIFT;

			// soft padding counts as a mismatch
			if (op == BAM_CSOFT_CLIP)
			{
				mismatches += len;
			}
			else if (op == BAM_CDEL)
			{
				mismatches += len;
				ref_pos += len;
			}
			else if (op == BAM_CINS)
			{
				mismatches += len;
				read_pos += len;
			}
			else if (op == BAM_CMATCH)
			{
				for (uint32_t j = 0; j < len; j++)
				{
					if (verbose)
					{
						cout << "REF: " << ref_pos << "  " << ref_string[ref_pos] << endl;
						cout << "READ: " << read_pos << "  " << read_string[read_pos] << endl;
					}
					if (ref_string[ref_pos] != read_string[read_pos])
          {
						mismatches++;
					}
					read_pos++;
					ref_pos++;
				}
			}
		}

		//#	print $a->qname . "\n$mismatches\n$cigar_string\n$ref_string\n$read_string\n" if ($mismatches);
		return mismatches;
  }

	string shifted_cigar_string(const alignment_wrapper& a, const cReferenceSequences& ref_seq_info)
	{
    bool verbose = true;
    
		string ref_seq = ref_seq_info[a.reference_target_id()].m_fasta_sequence.m_sequence;
		uint32_t ref_seq_index = 0;
		string read_seq = a.read_char_sequence();
		uint32_t read_seq_index = 0;
		vector<pair<char,uint16_t> > cigar_pair_array = a.cigar_pair_array();

    if (verbose)
    {
      cout << a.read_name() << endl;
      
      if (a.read_name() == "GW1ULQG02EIUY7")
      {
        cout << "debug" << endl;
      }
      cout << a.cigar_string() << endl;
    }
    
		for (vector<pair<char,uint16_t> >::iterator it = cigar_pair_array.begin(); it != cigar_pair_array.end(); it++)
		{
			char op = it->first;
			uint16_t len = it->second;

			if (op == 'D')
			{
				char base = ref_seq[ref_seq_index];
				bool all_same_base = true;
				for (uint32_t j = 1; j < len; j++)
					all_same_base = all_same_base && (ref_seq[ref_seq_index + j] == base);

				if (all_same_base)
				{
					uint32_t shift_amount = 0;
					while (ref_seq[ref_seq_index + len + shift_amount] == base)
						shift_amount++;

					if (shift_amount > 0)
					{
            if (verbose)
            {
              cout << "Shifting D: " << a.read_name() << endl;
            }
            
						if (it != cigar_pair_array.begin())
							(it - 1)->second += shift_amount;
						else
							cigar_pair_array.insert(it, make_pair<char,uint16_t>('M', shift_amount));

						if (it != cigar_pair_array.end())
							(it + 1)->second -= shift_amount;
						else
							cigar_pair_array.push_back(make_pair<char,uint16_t>('M', shift_amount));
					}

					ref_seq_index += len + shift_amount;
				}
				else
				{
					ref_seq_index += len;
				}
			}
			else if (op == 'I')
			{
				char base = read_seq[read_seq_index];
				bool all_same_base = true;
				for (uint32_t j = 1; j < len; j++)
					all_same_base = all_same_base && (read_seq[read_seq_index + j] == base);

				if (all_same_base)
				{
					uint32_t shift_amount = 0;
					while (read_seq[read_seq_index + len + shift_amount] == base)
						shift_amount++;

					if (shift_amount > 0)
					{
            if (verbose)
            {
              cout << "Shifting I: " << a.read_name() << endl;
            }
            
						if (it != cigar_pair_array.begin())
							(it - 1)->second += shift_amount;
						else
							cigar_pair_array.insert(it, make_pair<char,uint16_t>('M', shift_amount));

						if (it != cigar_pair_array.end())
							(it + 1)->second -= shift_amount;
						else
							cigar_pair_array.push_back(make_pair<char,uint16_t>('M', shift_amount));
					}

					read_seq_index += len + shift_amount;
				}
				else
				{
					read_seq_index += len;
				}
			}
			else
			{
				ref_seq_index += len;
				read_seq_index += len;
			}
		}

		stringstream shifted_cigar_string_ss;
		for (vector<pair<char,uint16_t> >::iterator it = cigar_pair_array.begin(); it != cigar_pair_array.end(); it++)
			shifted_cigar_string_ss << it->second << it->first;

		string shifted_cigar_string = shifted_cigar_string_ss.str();
    if (verbose)
    {
      cout << shifted_cigar_string << endl;
    }
    
		return shifted_cigar_string;
	}

} // breseq namespace

