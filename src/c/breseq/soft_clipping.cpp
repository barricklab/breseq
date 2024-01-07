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

#include "libbreseq/soft_clipping.h"

#include "libbreseq/alignment.h"
#include "libbreseq/identify_mutations.h"
#include "libbreseq/reference_sequence.h"
#include "libbreseq/output.h"

using namespace std;

namespace breseq {

bool is_read_soft_clipped(
                          Settings& settings,
                          Summary& summary,
                          cReferenceSequences& ref_seq_info,
                          const string&  output_file_name
                          )
{
  
}

void analyze_soft_clipping(
                           const string& bam_file_name,
                           const string& fasta_file_name,
                           const string& output_file_name,
                           const uint32_t minimum_clipped_bases
                           )
{
  ofstream out_file(output_file_name);
  out_file << join(make_vector<string>("read_name")("seq_id")("position")("direction")("num_bases")("strand"), ",") << endl;
  
  bam_file final_bam_file(bam_file_name, fasta_file_name, ios::in);
  
  alignment_list al;
  
  while (final_bam_file.read_alignments(al, false)) {
    
    uint32_t num_alignments = al.size();
    
    // We don't want any with multiple alignments...
    // this should have been resolved or they are uninformative.
    //if (num_alignments > 1) continue;
    
    for (alignment_list::iterator it = al.begin(); it != al.end(); it++) {
      bam_alignment& a = *it->get();
      
      // is it mapped?
      if (a.unmapped()) continue;
      
      // check if it has our mark for being multiply mapped and skip if so
      uint32_t num_equivalent_alignments;
      if (a.aux_get_i("X1", num_equivalent_alignments)) {
        if (num_equivalent_alignments > 1) continue;
      }
      
      // Check to see if this is one side of a junction match
      // in which case we need to ignore the soft trimming on one side
      // But, we don't know which side... so
      int32_t ignore_side = -1;
      uint32_t junction_side;
      // either 1 or 2, this tells us which side to ignore
      if (a.aux_get_i("XJ", junction_side)) {
        ignore_side = junction_side;
      }
      
      // NOTE: a.read_length() can be zero when there are multiple alignments
      //       for minimap2. To avoid these problems, we need to use the CIGAR string to calculate the length
      uint32_t read_length = a.cigar_query_length();
      
      uint32_t query_begin_soft_clipping = a.query_start_1() - 1;
      uint32_t query_end_soft_clipping = read_length - a.query_end_1();
    
      // debugging
      //if (a.read_name() == "e176333e-29ea-4179-be3c-890e0e65dc60") {
      //  cout << "DEBUG!" << endl;
      //}
      
      //cout << a.read_name() << endl;
      //cout << "  Read Length/Strand: " << read_length << "  " << a.strand() << endl;
      //cout << "  Reference: " << a.reference_start_1() << "-" << a.reference_end_1() << endl;
      //cout << "  Query: " << a.query_start_1() << "-" << a.query_end_1() << endl;
      //cout << "  Clipping: " << query_begin_soft_clipping << "   " << query_end_soft_clipping << endl;

              
      if (query_begin_soft_clipping > minimum_clipped_bases) {
        
        uint32_t clipping_coord = a.reference_start_1();
        int32_t clipping_direction = -1;
        uint32_t num_bases = query_begin_soft_clipping;
        
        //cout << "PRINTED BEGIN" << endl;
        out_file << join(
                         make_vector<string>
                         (a.read_name())
                         (final_bam_file.target_name(a))
                         (to_string(clipping_coord))
                         (to_string(clipping_direction))
                         (to_string(num_bases))
                         (to_string(a.strand())),
                         ",") << endl;
      }
      
      if (query_end_soft_clipping > minimum_clipped_bases) {

        uint32_t clipping_coord = a.reference_end_1();
        int32_t clipping_direction = +1;
        uint32_t num_bases = query_end_soft_clipping;
        
        //cout << "PRINTED END" << endl;
        out_file << join(
                         make_vector<string>
                         (a.read_name())
                         (final_bam_file.target_name(a))
                         (to_string(clipping_coord))
                         (to_string(clipping_direction))
                         (to_string(num_bases))
                         (to_string(a.strand())),
                         ",") << endl;
        
      }
      
      } // end for loop
      
    } // end while loop
    
  } // end function
 
} // namespace breseq

