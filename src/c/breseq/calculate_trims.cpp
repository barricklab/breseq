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


#include "libbreseq/calculate_trims.h"

#include "libbreseq/reference_sequence.h"
#include "libbreseq/alignment.h"

using namespace std;

namespace breseq {

const uint32_t k_max_repeat_length = 18;
  
/*
calc_trims

Currently each reference sequence is loaded entirely into memory, 
which isn't necessary.
*/


// Note: this doesn't quite get the ends of the genome right when repeats overlap there
// The "easy" way to do this correctly is to paste some sequence from each end of circular fragments on the other end
// and then correct everything for these offsets.
bool repeat_match ( const string& _in_seq, const uint32_t pos, const uint32_t repeat_size, const uint32_t repeat_num ) {

  const uint32_t base_offset = repeat_size * repeat_num;
  for (uint32_t i=0; i<repeat_size; i++) {
  
    //cerr << i << " " << _in_seq[pos+i] << " " << _in_seq[pos+i+base_offset] << endl;
    
    if (pos+i >= _in_seq.length()) return false;
    if (pos+i+base_offset >= _in_seq.length()) return false;
    if (_in_seq[pos+i] != _in_seq[pos+i+base_offset]) return false;
  }
  return true;
}

// writes straight to file
void calculate_trims_1 ( const string& _in_seq, const string& in_output_filename ) {
  SequenceTrims trims(_in_seq);
  trims.WriteFile(in_output_filename);
}

  
  
// Initialize trims for a sequence from the string
SequenceTrims::SequenceTrims(const string& _in_seq) 
: trim_data(NULL), m_length(0) 
{

  //cerr << _in_seq.length() << endl;
    
  // use one structure to avoid byte alignment issues when writing
  trim_data = new unsigned char[2*_in_seq.length()];
  m_length = _in_seq.length();

  assert(trim_data != NULL);
  memset( trim_data, 0, 2*m_length );
  
  const uint32_t left_trim_offset = 0;
  const uint32_t right_trim_offset = m_length;
  

  for (uint32_t pos_0=0; pos_0<m_length; pos_0++) {
    
    // always trim at least one bp
    uint32_t max_trim_length = 1;
  
    //compare starting at this nucleotide
    for (uint32_t repeat_size=1; repeat_size<=k_max_repeat_length; repeat_size++) {
      
      unsigned int repeat_num = 1;
      while (repeat_match(_in_seq, pos_0, repeat_size, repeat_num) ) {
        repeat_num++;
      }
    
      //cerr << (pos+1) << " " << repeat_size << " " << repeat_num << endl;
      
      if (repeat_num > 1) {
        uint32_t this_trim = repeat_num * repeat_size;
        if (this_trim > max_trim_length) {
          max_trim_length = this_trim;
        }
      }
    }
  
    // currently limited to an unsigned char
    if (max_trim_length > 255) max_trim_length = 255;
    uint8_t add_max_trim_length = max_trim_length;
  
    // update relevant trims
    for (int32_t offset=0; offset<add_max_trim_length; offset++) {
      if (pos_0 + offset >= _in_seq.length()) continue;
      
      uint8_t this_trim = offset + 1;
      //cerr << (right_trim_offset + pos + offset) << " " << (int)(this_trim) << " " << (int)(trim[right_trim_offset + pos + offset]) << endl;
      trim_data[right_trim_offset + pos_0 + offset] = max(this_trim, trim_data[right_trim_offset + pos_0 + offset]);
    }
 
    for (int32_t offset=add_max_trim_length; offset>=0; offset--) {
      if (pos_0 + offset >= _in_seq.length()) continue;

      uint8_t this_trim = add_max_trim_length - offset;      
      //cerr << (left_trim_offset + pos + offset) << " " << (int)(this_trim) << " " << int(trim[left_trim_offset + pos + offset]) << endl;
      trim_data[left_trim_offset + pos_0 + offset] = max(this_trim, trim_data[left_trim_offset + pos_0 + offset]);
    }
    
    
    // debug
    //if (pos == 10000) break;
  }

//debugging...
//  for (uint32_t i=0; i<_in_seq.length(); i++)
//  {
//    cerr << (i+1) << " " << (int)left_trim[i] << " " << (int)right_trim[i] << endl;
//  }
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
	for(int32_t i=0; i<nseq; ++i) {

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

  
// Trims are from each end
Trims edge_trims_for_sequence(const string &_in_seq)
{
  Trims this_trims;
  this_trims.L = 1;
  this_trims.R = 1;
  
  // Left edge
  for (uint32_t repeat_size=1; repeat_size<=k_max_repeat_length; repeat_size++) {
    
    uint32_t repeat_num = 1;
    uint32_t pos_0 = 0;
    
    while (repeat_match(_in_seq, pos_0, repeat_size, repeat_num)) {
      repeat_num++;
    }
    
    //cerr << (pos+1) << " " << repeat_size << " " << repeat_num << endl;
    
    if (repeat_num > 1) {
      this_trims.L = repeat_num * repeat_size;
    }
  }
  
  string rc_seq = reverse_complement(_in_seq);
  
  // Right edge
  for (uint32_t repeat_size=1; repeat_size<=k_max_repeat_length; repeat_size++) {
    
    uint32_t repeat_num = 1;
    uint32_t pos_0 = 0;
    
    while (repeat_match(rc_seq, pos_0, repeat_size, repeat_num)) {
      repeat_num++;
    }
    
    //cerr << (pos+1) << " " << repeat_size << " " << repeat_num << endl;
    
    if (repeat_num > 1) {
      this_trims.R = repeat_num * repeat_size;
    }
  }
  
  return this_trims;
}

  
// @JEB 2015-07-29 This code also trims based on ends of mapped read sequence
// Assuming that fixing mismatches near the ends of the reads breaking match
// @JEB 2015-11-12 This is necessary for avoiding certain misalignments

#define TRIM_READ_ENDS

Trims get_alignment_trims(const alignment_wrapper& a, const SequenceTrimsList& trims)
{
  bool verbose = false;
  
  // which reference sequence?
  uint32_t tid = a.reference_target_id();
  
  
#ifdef TRIM_READ_ENDS
  
  // Read-based trims
  Trims et = edge_trims_for_sequence(a.query_char_sequence());
  
  // Genomic-based trims
  Trims gt;
  gt.L= trims[tid].left_trim_0(a.reference_start_0());
  gt.R = trims[tid].right_trim_0(a.reference_end_0());
  
  Trims t;
  t.L = max(gt.L, et.L);
  t.R = max(gt.R, et.R);
  
  // debug output
  //cerr << a.read_name() << endl;
  //  cerr << "start: " << a.reference_start_1() << " end: " << a.reference_end_1() << endl;
  //  cerr << "left read  : " << et.L << " right read  : " << et.R << endl;
  //  cerr << "left genome: " << gt.L << " right genome: " << gt.R << endl;
  //  cerr << "left offset: " << a.query_start_0() << " right offset : " << (a.read_length() - a.query_end_1()) << endl;
  //  cerr << "left final : " << t.L  << " right final : " << t.R << endl;
  
#else
  
  // Begin alternative block alternative
  Trims t;
  t.L= trims[tid].left_trim_0(a.reference_start_0());
  t.R = trims[tid].right_trim_0(a.reference_end_0());
  // End alternative code block
  
#endif
  
  // Offset based on where our match actually was in the read
  
  t.L += a.query_start_0();
  t.R += a.read_length() - a.query_end_1();
  
  return t;
}

void read_trims(SequenceTrimsList& trims, const cReferenceSequences& ref_seqs, const string &in_trims_file_name ) 
{
  trims.resize(ref_seqs.size());
  for(uint32_t i = 0; i < ref_seqs.size(); i++) {
    string this_file_name = Settings::file_name(in_trims_file_name, "@", ref_seqs[i].m_seq_id);
    trims[i].ReadFile(this_file_name, ref_seqs[i].m_length);
  }
}


} // breseq

