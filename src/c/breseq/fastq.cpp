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

#include "breseq/fastq.h"

using namespace std;

namespace breseq {
  
  /*
   analyze_fastq
   
   convert if necessary
   print statistics about the fastq file

   */

  AnalyzeFastq analyze_fastq(const string &file_name, const string &convert_file_name) {
    
    // Set up maps between formats
    map<string,uint8_t> format_to_chr_offset;
    format_to_chr_offset["SANGER"] = 33;
    format_to_chr_offset["SOLEXA"] = 64;
    format_to_chr_offset["ILLUMINA_1.3+"] = 64;
    
    // Summary information that will be printed at the end
    uint32_t max_read_length = 0;
    int min_quality_score = 255;
    int max_quality_score = 0;
    uint32_t num_bases = 0; 
    uint32_t num_reads = 0;
    
    // Process the input file, one sequence at a time
    cFastqFile input_fastq_file(file_name.c_str(), fstream::in);
    assert(input_fastq_file.is_open());

    cFastqSequence on_sequence;
    while (input_fastq_file.read_sequence(on_sequence)) {
      
      //increment read number
      num_reads++;
      
      //check sequence length
      if( on_sequence.m_sequence.size() > max_read_length ) max_read_length = on_sequence.m_sequence.size();

      //add current sequence length to number of bases
      num_bases += on_sequence.m_sequence.size();
      
      vector<int> numerical_quality_scores;
      
      if( on_sequence.m_numerical_qualities )
        convert_numeric_scores(on_sequence);
      
      //cout << "I'm here: " << num_reads << endl;
      
      //iterate through sequence grabbing the associated scores
      for (uint32_t i=0; i<on_sequence.m_qualities.size(); i++) {
        int this_score(uint8_t(on_sequence.m_qualities[i]));
        if( this_score > max_quality_score ) max_quality_score = this_score;
        if( this_score < min_quality_score ) min_quality_score = this_score;

      format_to_chr_offset["SANGER"];
      }
    }
    input_fastq_file.close();
    
    // Default is SANGER
    string quality_format = "SANGER";
    
    // Typical range: (-5, 40) + 64
    if (min_quality_score >= format_to_chr_offset["SOLEXA"] - 5) {
      quality_format = "SOLEXA";
    } 
    // Typical range:  (0, 40) + 64
    if (min_quality_score >= format_to_chr_offset["ILLUMINA_1.3+"]) {
      quality_format = "ILLUMINA_1.3+";
    }
    
    //cerr << "min_quality_score "     << (int)min_quality_score  << endl;
    //cerr << "max_quality_score "     << (int)max_quality_score  << endl;

    string converted_fastq_name(file_name);
    
    //std::cout << m_quality_format << std::endl;
    if(( quality_format != "SANGER" ) || (input_fastq_file.needs_conversion()) )  {
      cerr << "  Converting/repairing FASTQ file..." << endl;
      cerr << "  Original format: " << quality_format << " New format: SANGER"<< endl;

      cFastqQualityConverter fqc(quality_format, "SANGER");
      
      // re-open input for another pass
      input_fastq_file.open(file_name.c_str(), fstream::in);
      
      //open output converted file
      cFastqFile output_fastq_file(convert_file_name.c_str(), fstream::out);

      // recalculate min and max quality scores from table
      cFastqSequence min_max_sequence;
      min_max_sequence.m_qualities.append(1,min_quality_score);
      min_max_sequence.m_qualities.append(1,max_quality_score);
      fqc.convert_sequence(min_max_sequence);
      min_quality_score = (uint8_t)min_max_sequence.m_qualities[0];
      max_quality_score = (uint8_t)min_max_sequence.m_qualities[1];
      
      // (much faster than looking through all qualities again)
      
      while (input_fastq_file.read_sequence(on_sequence)) {
        
        // truncate second name name
        on_sequence.m_name_plus = "+";
      
        if( on_sequence.m_numerical_qualities )
          convert_numeric_scores(on_sequence);
        
        // fastq quality convert
        fqc.convert_sequence(on_sequence);
        
        // replace "/" characters
        size_t pos = 0;
        pos = on_sequence.m_name.find("/", pos);
        while (pos != string::npos)
        {
          on_sequence.m_name.replace(pos, 1, "_");
          pos = on_sequence.m_name.find("/", pos+1);
        }
        
        // convert base qualities
        output_fastq_file.write_sequence(on_sequence);
        
      }
      input_fastq_file.close();
      output_fastq_file.close();
      
      converted_fastq_name = convert_file_name;
    }
    
    // quality scores are in SANGER at this point
    min_quality_score = min_quality_score - format_to_chr_offset["SANGER"];
    max_quality_score = max_quality_score - format_to_chr_offset["SANGER"];
        
    //cerr << "min_quality_score "     << (int)min_quality_score  << endl;
    //cerr << "max_quality_score "     << (int)max_quality_score  << endl;
    
    // Output information to stdout
    cout << "max_read_length "       << max_read_length         << endl;
    cout << "num_reads "             << num_reads               << endl;
    cout << "min_quality_score "     << (int)min_quality_score  << endl;
    cout << "max_quality_score "     << (int)max_quality_score  << endl;
    cout << "num_bases "             << num_bases               << endl;
    cout << "original_qual_format "  << quality_format          << endl;
    cout << "qual_format "           << "SANGER"                << endl;
    cout << "converted_fastq_name "  << converted_fastq_name    << endl;
	AnalyzeFastq retval = { max_read_length, num_reads, min_quality_score, max_quality_score, num_bases, quality_format, "SANGER", converted_fastq_name };
	return retval;
  }

  // constructor
  cFastqQualityConverter::cFastqQualityConverter(const string &from_quality_type, const string &to_quality_type)
  {
    // Set up maps between formats
    map<string,uint8_t> format_to_chr_offset;
    format_to_chr_offset["SANGER"] = 33;
    format_to_chr_offset["SOLEXA"] = 64;
    format_to_chr_offset["ILLUMINA_1.3+"] = 64;
    
    map<string,string> format_to_quality_type;
    format_to_quality_type["SANGER"] = "PHRED";
    format_to_quality_type["SOLEXA"] = "SOLEXA";
    format_to_quality_type["ILLUMINA_1.3+"] = "PHRED";
    
    this->resize(255);
    for (uint16_t i = 0; i<=255; i++) {
      (*this)[i] = 0;
    }
    
    for (uint16_t from_chr = 0; from_chr<=255; from_chr++) {

      int16_t from_quality = from_chr - format_to_chr_offset[from_quality_type];
      
      // Calculate the probability of error
      double probability_of_error;
      
      if (format_to_quality_type[from_quality_type] == "SOLEXA") {
        probability_of_error = 1 / (1+pow(10,(double)from_quality/10));
      } else if (format_to_quality_type[from_quality_type] == "PHRED") {
        probability_of_error = pow(10,-(double)from_quality/10);
      } else {
        cerr << "Unknown base quality score type: " << from_quality_type << endl;
        exit(-1);
      }
      
      //Convert back to quality score
      int to_quality;
            
      if (format_to_quality_type[to_quality_type] == "SOLEXA") {
        to_quality = round(-10 * log((1-probability_of_error)/probability_of_error) / log(10));
      } else if (format_to_quality_type[to_quality_type] == "PHRED") {
        double t_quality = round(-10 * log(probability_of_error) / log(10));

        to_quality = round(t_quality);
      } else {
        cerr << "Unknown base quality score type: " << to_quality_type << endl;
        exit(-1);
      }
      
      int16_t to_chr = from_quality + format_to_chr_offset[to_quality_type];

      
      // May be out of range
      if ((to_chr < 0) || (to_chr > 255)) continue;

      (*this)[(uint8_t)from_chr] = (uint8_t)to_chr;
      
      // Debug
      //cerr << from_chr << " => " << to_chr << endl;
    }     
    
  }

  void cFastqQualityConverter::convert_sequence(cFastqSequence &seq) {
    
    for(uint32_t i=0; i < seq.m_qualities.size(); i++)
    {
      //seq.m_qualities.replace(i,1,1,(*this)[seq.m_qualities[i]]);
      seq.m_qualities[i]= (*this)[seq.m_qualities[i]];
    }
  }
  
  //constructor
  cFastqFile::cFastqFile() :
    fstream(), m_current_line(0), m_file_name(""), m_needs_conversion(false)
  {
  }
 
  
  cFastqFile::cFastqFile(const std::string &file_name, std::ios_base::openmode mode) :
    fstream(file_name.c_str(), mode), m_current_line(0), m_file_name(file_name), m_needs_conversion(false)
  { 
    assert(!(*this).fail());
  }

  // read one sequence record from the file
  bool cFastqFile::read_sequence(cFastqSequence &sequence) {
    
    // We're done, no error
    if (this->eof()) return false;
    
    uint32_t count = 0;
    string line;
    
    // get the next four lines
    while (count < 4) {
      std::getline(*this, line);
      m_current_line++;
      
      // Didn't get a first line, then we ended correctly
      if (this->eof()) {
        if (count == 0) {
          return false;
        } else {
          uint32_t last_valid_line = floor((m_current_line-1)/4.0) * 4;
          fprintf(stderr, "Incomplete FASTQ sequence record found at end of file.\nFile %s\nLine: %d\n", m_file_name.c_str(), m_current_line-1);
          fprintf(stderr, "You may be able to repair this damage and salvage the reads before this point with the command:\n");
          fprintf(stderr, "  head -n %u %s > new.fastq\n", last_valid_line, m_file_name.c_str());
          fprintf(stderr, "Then use \"new.fastq\" as input.\n");
          exit(-1);
        }
      }
      
      switch (count) {
        case 0:
          if( line[0] != '@' ) {
            fprintf(stderr, "FASTQ sequence record does not begin with @NAME line.\nFile %s\nLine: %d\n", m_file_name.c_str(), m_current_line);
            exit(-1);
          }
          sequence.m_name = line.substr(1,string::npos);
          if (m_current_line == 1) {
            if (sequence.m_name.find("/") != string::npos) m_needs_conversion = true;
          }
          break;
        case 1:
          sequence.m_sequence = line;
          
          // check for extraneous DOS ending
          if (m_current_line == 1) {
            if( sequence.m_sequence[sequence.m_sequence.size()-1] == '\r') {
              sequence.m_sequence.resize(sequence.m_sequence.size()-1);
              m_needs_conversion = true;
            }
          }
          
          for (uint32_t i=0; i<sequence.m_sequence.size(); i++) {
            
            // convert to uppercase and require
            // reformatting if this was necessary
            switch (sequence.m_sequence[i]) {
              case 'a':
                sequence.m_sequence.replace(i,1,1,'A');
                m_needs_conversion = true;
                break;
                
              case 't':
                sequence.m_sequence.replace(i,1,1,'T');
                m_needs_conversion = true;
                break;
                
              case 'c':
                sequence.m_sequence.replace(i,1,1,'C');
                m_needs_conversion = true;
                break;
                
              case 'g':
                sequence.m_sequence.replace(i,1,1,'G');
                m_needs_conversion = true;
                break;

              case 'n':
                sequence.m_sequence.replace(i,1,1,'N');
                m_needs_conversion = true;
                break;
            }

            
            if(sequence.m_sequence[i] != 'A' && 
               sequence.m_sequence[i] != 'T' && 
               sequence.m_sequence[i] != 'G' && 
               sequence.m_sequence[i] != 'C' && 
               sequence.m_sequence[i] != 'N') {
              
              fprintf(stderr, "FASTQ sequence character not allowed %c.\nSequence: %s\nFile %s\nLine: %d\n", 
                      sequence.m_sequence[i], sequence.m_sequence.c_str(), m_file_name.c_str(), m_current_line);
              exit(-1);
            }
          }
          
          break;
        case 2:
          
          //Only need to see if the first character is a +
          if( line[0] != '+' ) {
            fprintf(stderr, "FASTQ sequence record does not contain +NAME line.\nFile %s\nLine: %d\n", m_file_name.c_str(), m_current_line);
            exit(-1);
          }
          // Could optionally check to see if the name after the + was either absent or identical to the earlier name
          sequence.m_name_plus = line.substr(1,string::npos);

          break;
        case 3:
          sequence.m_qualities = line;
          
          if( sequence.m_sequence.size() != sequence.m_qualities.size() ) {
            vector<string> numerical_qualities(split(sequence.m_qualities, " "));
            
            if( sequence.m_sequence.size() != numerical_qualities.size() ) {
              fprintf(stderr, "FASTQ sequence record has different SEQUENCE and QUALITY lengths.\nFile %s\nLine: %d\n", m_file_name.c_str(), m_current_line);
              exit(-1);
            }
            else if( sequence.m_sequence.size() == numerical_qualities.size() ) {
              sequence.m_numerical_qualities = true; 
            }
          }
          else {
            sequence.m_numerical_qualities = false;
          }

          break;
      }
      
      count++;
  }
    
    return true;
  }

  void cFastqFile::write_sequence(const cFastqSequence &sequence) {
    (*this) << '@' << sequence.m_name << std::endl;
    (*this) << sequence.m_sequence << std::endl;
    (*this) << '+' << sequence.m_name_plus << std::endl;
    (*this) << sequence.m_qualities << std::endl;
  }
  
  void convert_numeric_scores(cFastqSequence &sequence) {
    vector<string> numerical_qualities_as_string( split(sequence.m_qualities, " " ) );
    
    sequence.m_qualities.resize(numerical_qualities_as_string.size());
    
    for (uint32_t this_quality = 0; this_quality < numerical_qualities_as_string.size(); ++this_quality) {
      sequence.m_qualities[this_quality] = from_string<int>(numerical_qualities_as_string[this_quality]) + 64;
    }
  }
  
} // breseq namespace

