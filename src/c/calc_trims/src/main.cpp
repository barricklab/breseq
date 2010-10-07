#include <boost/program_options.hpp>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>

#include <faidx.h>
#include <bam.h>

using namespace std;

/*
calc_trims

Currently each reference sequence is loaded entirely into memory, 
which isn't necessary.
*/


// Note: this doesn't quite get the ends of the genome right when repeats overlap there
bool repeat_match ( const string& _in_seq, const uint32_t pos, const uint32_t repeat_size, const uint32_t repeat_num ) {

  const uint32_t base_offset = repeat_size * repeat_num;
  for (uint32_t i=0; i<repeat_size; i++) {
  
    //cerr << i << " " << _in_seq[pos+i] << " " << _in_seq[pos+i+base_offset] << endl;
    
    if (pos+i >= _in_seq.length()) return false;
    if (_in_seq[pos+i] != _in_seq[pos+i+base_offset]) return false;
  }
  return true;
}

void calculate_trims_1 ( const string& _in_seq, const string& in_output_filename ) {

  //cerr << _in_seq.length() << endl;

  // use one structure to avoid byte alignment issues when writing
  uint8_t * trim = new unsigned char[2*_in_seq.length()];
  const uint32_t left_trim_offset = 0;
  const uint32_t right_trim_offset = _in_seq.length();
  
  uint32_t max_repeat_length = 18;

  for (unsigned int pos=0; pos<_in_seq.length(); pos++)
  {
  
    // always trim at least one bp
    uint32_t max_trim_length = 1;
  
    //compare starting at this nucleotide
    for (uint32_t repeat_size=1; repeat_size<=max_repeat_length; repeat_size++)
    {
    
      unsigned int repeat_num = 1;
      while (repeat_match(_in_seq, pos, repeat_size, repeat_num) ) 
      {
        repeat_num++;
      }
    
      //cerr << (pos+1) << " " << repeat_size << " " << repeat_num << endl;
      
      if (repeat_num > 1)
      {
        uint32_t this_trim = repeat_num * repeat_size;
        if (this_trim > max_trim_length)
        {
          max_trim_length = this_trim;
        }
      }
    }
  
    // currently limited to an unsigned char
    if (max_trim_length > 255) max_trim_length = 255;
    uint8_t add_max_trim_length = max_trim_length;
  
    // update relevant trims
    for (int32_t offset=0; offset<add_max_trim_length; offset++)
    {
      if (pos + offset >= _in_seq.length()) break;
      
      uint8_t this_trim = offset + 1;
   //   cerr << (pos+offset+1) << " " << this_trim << " " << right_trim[pos + offset] << endl;
      trim[right_trim_offset + pos + offset] = max(this_trim, trim[right_trim_offset + pos + offset]);
    }
 
    for (int32_t offset=add_max_trim_length; offset>=0; offset--)
    {
      uint8_t this_trim = add_max_trim_length - offset;      
  //    cerr << (pos+offset+1) << " " << this_trim << " " << right_trim[pos + offset] << endl;
      trim[left_trim_offset + pos + offset] = max(this_trim, trim[left_trim_offset + pos + offset]);
    }
    
    // debug
    
    //if (pos == 10000) break;
  }

  // Write out the trims to abinary file
  
  //cerr << in_output_filename << endl;
  ofstream out(in_output_filename.c_str(), ios::binary);
  out << trim;
  out.close();

//debugging...
//  for (uint32_t i=0; i<_in_seq.length(); i++)
//  {
//    cerr << (i+1) << " " << (int)left_trim[i] << " " << (int)right_trim[i] << endl;
//  }
      
  delete[] trim;
}

void calculate_trims( const string& in_fasta, const string& in_output_path) {

  // Load the sequence index
  string fai_filename(in_fasta);
  fai_filename+=".fai";
  
  //cerr << fai_filename << endl;
  
  bam_header_t* bam_header = sam_header_read2(fai_filename.c_str());
  assert(bam_header);
  int nseq = bam_header->n_targets;

  faidx_t* fasta_index = fai_load(in_fasta.c_str());
  assert(fasta_index);


	// load all the reference sequences:
	for(int i=0; i<nseq; ++i) {

		cerr << "  REFERENCE: " << bam_header->target_name[i] << endl;
		cerr << "  LENGTH: " << bam_header->target_len[i] << endl;
    
    int len;
    const char* cseq = fai_fetch(fasta_index, bam_header->target_name[i], &len);
    
    assert(cseq);
    assert(len > 0);
    assert(static_cast<unsigned int>(len) == bam_header->target_len[i]);
    
    const string seq(cseq);
    
    string output_filename(in_output_path);
    output_filename += "/";
    output_filename += bam_header->target_name[i];
    output_filename += ".trims";

    calculate_trims_1(seq, output_filename);
	}
  
  bam_header_destroy(bam_header);
}



/*! Identify mutations.
 
 This file only does the command-line parsing bit; the real work is over in
 identify_mutations.cpp.
 */
int main(int argc, char* argv[]) {
	using namespace std;
	namespace po = boost::program_options;
	
	// setup and parse configuration options:
	po::options_description cmdline_options("Allowed options");
	cmdline_options.add_options()
	("help,h", "produce this help message")
	("fasta,f", po::value<string>(), "FASTA file of reference sequence")
	("output,o", po::value<string>(), "output directory")
  ;
  
	po::variables_map options;
	po::store(po::parse_command_line(argc, argv, cmdline_options), options);
	po::notify(options);
	
	// make sure that the config options are good:
	if(options.count("help")
		 || !options.count("fasta")
		 || !options.count("output")
		 ) {
		cout << "Usage: identify_mutations --bam <sequences.bam> --fasta <reference.fasta> --error_dir <path> --genome_diff <path> --output <path> --readfiles <filename> --coverage_dir <dirname> [--minimum-quality-score 3]" << endl;
		cout << cmdline_options << endl;
		return -1;
	}                       
  
	// attempt to calculate error calibrations:
	try {
		calculate_trims(options["fasta"].as<string>(),options["output"].as<string>());
		} catch(...) {
		// failed; 
		return -1;
	}
	
	return 0;
}


