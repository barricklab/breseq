/*****************************************************************************
 *
 * AUTHORS
 *
 *    Jeffrey E. Barrick <jeffrey.e.barrick@gmail.com>
 *    David B. Knoester
 *
 * LICENSE AND COPYRIGHT
 *
 *    Copyright (c) 2008-2010 Michigan State University
 *    Copyright (c) 2011-2017 The University of Texas at Austin
 *
 *    breseq is free software; you can redistribute it and/or modify it under the
 *    terms the GNU General Public License as published by the Free Software
 *    Foundation; either version 1, or (at your option) any later version.
 *
 *****************************************************************************/

#include "libbreseq/coverage_output.h"

using namespace std;

namespace breseq
{
  
void coverage_output::plot(const string& region, const string& output_file_name, const uint32_t resolution)
{
  
  uint32_t target_id, start_pos, end_pos, insert_start, insert_end;
  this->parse_region(region, target_id, start_pos, end_pos, insert_start, insert_end);
  string seq_id = this->target_name(target_id);
  
  ASSERT( (seq_id.length() > 0) || !start_pos || !end_pos, "Invalid region: \"" + region + "\"");

  // Load summary file of fitting coverage from breseq output
  if (m_show_average) {
    ASSERT(file_exists(m_average_file_name.c_str()), "Could not open file with fit coverage: " + m_average_file_name + "\n");
    m_summary.retrieve(m_average_file_name);
  }
  
  // extend the region and re-check endpoints
  int32_t extended_start = start_pos - m_shaded_flanking;
  if (extended_start < 1) extended_start = 1;

  int32_t extended_end = end_pos + m_shaded_flanking;
  if (extended_end > static_cast<int32_t>(this->target_length(target_id)) ) extended_end = this->target_length(target_id);
  	
	string extended_region = seq_id + ":" + to_string(extended_start) + "-" + to_string(extended_end);  
	uint32_t size = extended_end - extended_start + 1;
	
  string _output_file_name(output_file_name);
  
  // auto-construct filename if none was provided
  if (output_file_name.length() == 0) 
  {
    _output_file_name = region;
    _output_file_name = substitute(_output_file_name, ":", "_");
    _output_file_name += "." + m_output_format;
  }
	
  pid_t pid = getpid();
	string tmp_coverage = m_intermediate_path + "/" + to_string(pid) + "_" + to_string(m_thread_id) + ".coverage.tab";
	
  this->table(extended_region, tmp_coverage, resolution);
  
  // default is zero
  m_reference_average_coverage = m_show_average ? m_summary.references.reference[seq_id].coverage_average : 0;
	
  string log_file_name = m_intermediate_path + "/" + to_string(pid) + ".r.log";
  string command = "R --vanilla < " + double_quote(m_r_script_file_name)+ " > " + double_quote(log_file_name);
  command += " --args";
  command += " in_file=" + double_quote(tmp_coverage);
  command += " out_file=" + double_quote(_output_file_name);
  command += " pdf_output=";
  command += ((m_output_format=="pdf") ? "1" : "0");
  command += " total_only=";
  command += ((m_total_only) ? "1" : "0");
  command += " window_start=" + to_string(start_pos);
  command += " window_end=" + to_string(end_pos);
  command += " avg_coverage=" + to_string(m_reference_average_coverage);
  command += " fixed_coverage_scale=" + ( m_reference_average_coverage ? to_string<double>(m_reference_average_coverage * m_fixed_coverage_scale) : to_string<double>(m_fixed_coverage_scale) );


  
	SYSTEM(command, true, false, false); //NOTE: Not escaping shell characters here.
	
	remove(tmp_coverage.c_str());
	remove(log_file_name.c_str());
  remove("Rplots.pdf");
}

void coverage_output::table(const string& region, const string& output_file_name, uint32_t resolution)
{
  uint32_t target_id, start_pos, end_pos, insert_start, insert_end;
  this->parse_region(region, target_id, start_pos, end_pos, insert_start, insert_end);
  string seq_id = this->target_name(target_id);
  
  // Load summary file of fitting coverage from breseq output
  if (m_show_average) {
    ASSERT(file_exists(m_average_file_name.c_str()), "Could not open file with fit coverage: " + m_average_file_name + "\n");
    m_summary.retrieve(m_average_file_name);
  }
  
	uint32_t size = end_pos - start_pos + 1;

  // Resolution of 0 implies no downsampling, otherwise adjust to get close to resolution # of pts
  uint32_t downsample = 1;
  if (resolution != 0)
  {
    downsample = static_cast<uint32_t>(floor(static_cast<double>(size) / static_cast<double>(resolution)));
    if (downsample < 1) downsample = 1;
  }
  
  m_reference_average_coverage = m_show_average ? m_summary.references.reference[seq_id].coverage_average : 0;

  m_output_table.open(output_file_name.c_str());
  
  if (m_read_begin_output_file_name.length() > 0) m_read_begin_output.open(m_read_begin_output_file_name.c_str());
  if (m_gc_output_file_name.length() > 0) m_gc_output.open(m_gc_output_file_name.c_str());
  
  if (m_total_only) {
    m_output_table << "position" << "\t" << "ref_base" << "\t"
      << "unique_cov" << "\t" << "redundant_cov" << "\t" "total_cov"
      << std::endl;
  } else {
    m_output_table << "position" << "\t" << "ref_base" << "\t"
      << "unique_top_cov" << "\t" << "unique_bot_cov" << "\t"
      << "redundant_top_cov" << "\t" << "redundant_bot_cov" << "\t"
      << "raw_redundant_top_cov" << "\t" << "raw_redundant_bot_cov" << "\t"
      << "unique_top_begin" << "\t" << "unique_bot_begin"
      << std::endl;
  }
    
  this->clear();
  
  // pileup handles everything else, including into file
  if (region.length() == 0)
    this->do_pileup();
  else
    this->do_pileup(region, true, downsample);
  
  // Write the averages over the region and the reference average (if available) as a commented addendum at the end
  if (m_show_average) {
    m_output_table << "#\treference_unique_average_cov\t" << m_reference_average_coverage << std::endl;
  }
   
  m_output_table << "#\tregion_unique_average_cov\t" << ( m_region_average_unique_coverage / m_region_num_positions) << std::endl;
  m_output_table << "#\tregion_repeat_average_cov\t" << ( m_region_average_repeat_coverage / m_region_num_positions) << std::endl;
  m_output_table << "#\tregion_average_cov\t" << ( m_region_average_coverage / m_region_num_positions) << std::endl;
  m_output_table << "#\tnumber_of_positions\t" << m_region_num_positions << std::endl;
  
  m_output_table.close();
}
  
/*! Clear all statistics.
 */  
void coverage_output::clear()
{
  m_read_begin_top_bins.clear();
  m_read_begin_bot_bins.clear();
  m_ref_begin_top_bins.clear();
  m_ref_begin_bot_bins.clear();
  m_gc_content_bins.clear();
  m_region_average_unique_coverage = 0;
  m_region_average_repeat_coverage = 0;
  m_region_average_coverage = 0;
  m_region_num_positions = 0;
}

/*! Called for each alignment.
 */
void coverage_output::pileup_callback(const breseq::pileup& p) {
  
  char* refseq = p.reference_sequence(); // reference sequence for this target
  uint32_t pos = p.position_1();
  
  // don't handle indels before first position
  if (pos==0) return;
  
  // print positions not called because there were no reads
  for (uint32_t i=m_last_position_1+1; i<pos; i++) {
    if (m_total_only) {
      m_output_table << i << "\t" << refseq[i-1] << "\t" << 0 << "\t" << 0 << "\t" << 0
      << std::endl;
    } else {
      m_output_table << i << "\t" << refseq[i-1] << "\t" << 0 << "\t" << 0 << "\t"
        << 0 << "\t" << 0 << "\t" << 0 << "\t" << 0 << "\t" << 0 << "\t" << 0 << std::endl;
    }
    
    // Add to the positions, but don't add to the coverage
    m_region_num_positions++;
  }
  
  // catches this position

  if (m_read_begin_output.is_open()) {

    for (uint32_t j=m_last_position_1+1; j<=pos; j++) {
      if (j>=3) {
        std::string read_begin_s;
        {
          for (int i=1; i<=3; i++) {
            read_begin_s += complement_base_char(p.reference_base_char_1(j-i+1));
          }
        }
        m_ref_begin_bot_bins[read_begin_s]++;    
      }
      if (j<p.target_length()-3) {
        std::string read_begin_s;
        for (int i=1; i<=3; i++) {
          read_begin_s += p.reference_base_char_1(j+i-1);
        }        
        m_ref_begin_top_bins[read_begin_s]++;    
      }
    }
    
  }
  
  uint8_t ref_base = refseq[pos-1];
  uint32_t unique_cov[2] = {0,0};
  double redundant_cov[2] = {0.0, 0.0};
  uint32_t raw_redundant_cov[2] = {0,0};
  uint32_t unique_begin_reads[2] = {0,0};
  
  // for each alignment within this pileup:
  for(pileup::const_iterator a=p.begin(); a!=p.end(); ++a) {
    
    // skip deletions entirely, they are handled by adjacent matching positions
    if(a->is_del()) {
      continue;
    }
    
    uint32_t redundancy = a->redundancy();
    uint32_t reversed = a->reversed();
    
    bool first_base_matched;
    bool this_is_first_base;
    if (!reversed) { 
      // top strand
      first_base_matched = (a->query_start_1() == 1);
      this_is_first_base = (a->query_position_1() == 1); 
    } else {
      // bottom strand
      first_base_matched = (a->query_end_1() == a->read_length());
      this_is_first_base = (a->query_position_1() == a->read_length());     
    }
    
    if (redundancy == 1)
    {	    
      unique_cov[reversed]++;
      if (this_is_first_base) {
        unique_begin_reads[reversed]++;
      }
      
      
      if (this_is_first_base && m_read_begin_output.is_open()) {
        
        std::string read_begin_s;
        if (!reversed) {
          
          for (uint32_t i=1; i<=3; i++) {
            if (i > a->read_length()) break;
            base_bam bb = a->read_base_bam_1(i);
            if ( !_base_bam_is_N(bb) ) {
              read_begin_s += basebam2char(bb);
            }
          }
          
          // all must be not N
          if (read_begin_s.length() == 3) {
            m_read_begin_top_bins[read_begin_s]++;    
          }
          
          /*
           std::cout << "Pos " << pos << " Top strand " << read_begin_s << "  " << a->query_char_sequence() << std::endl;
           */
        } else {
          for (uint32_t i=1; i<=3; i++) {
            if (a->read_length()-i+1 < 1) break;
            base_bam bb = a->read_base_bam_1(a->read_length()-i+1);
            if ( !_base_bam_is_N(bb) ) {
              read_begin_s += basebam2char(complement_base_bam(bb));
            }
          }
          // all must be not N
          if (read_begin_s.length() == 3) {
            m_read_begin_bot_bins[read_begin_s]++;
          }
          /*
           std::cout << "Pos " << pos << " Bottom strand " << read_begin_s << "  " << a->query_char_sequence() << std::endl;
           */
        }
      }
    }
    else
    {
      raw_redundant_cov[reversed]++;
      redundant_cov[reversed] += 1.0/redundancy;			
    }		
  }
  
  
  //output
  
  m_region_num_positions++;
  m_region_average_unique_coverage += unique_cov[0] + unique_cov[1];
  m_region_average_repeat_coverage += redundant_cov[0] + redundant_cov[1];
  m_region_average_coverage += unique_cov[0] + unique_cov[1] + redundant_cov[0] + redundant_cov[1];
  
  if (m_total_only) {
    m_output_table << pos << "\t" << ref_base << "\t"
      << (unique_cov[0] + unique_cov[1]) << "\t"
      << (redundant_cov[0] + redundant_cov[1]) << "\t"
      << (unique_cov[0] + unique_cov[1] + redundant_cov[0] + redundant_cov[1])
      << std::endl;
  } else {
    m_output_table << pos << "\t" << ref_base << "\t"
      << unique_cov[0] << "\t" << unique_cov[1] << "\t"
      << redundant_cov[0] << "\t" << redundant_cov[1] << "\t"
      << raw_redundant_cov[0] << "\t" << raw_redundant_cov[1] << "\t"
      << unique_begin_reads[0] << "\t" << unique_begin_reads[1]
      << std::endl;
  }
}
  
  

/*! Called at the end of the pileup.
 */
void coverage_output::at_target_end(const uint32_t tid) {
  
  char* refseq = get_refseq(tid); // reference sequence for this target
  uint32_t pos = target_length(tid)+1;
  
  if (m_read_begin_output.is_open()) {
    
    // catches this position
    for (uint32_t j=m_last_position_1+1; j<pos; j++) {
      if (j>=3) {
        std::string read_begin_s;
        for (int i=1; i<=3; i++) {
          read_begin_s += complement_base_char(reference_base_char_1(tid, j-i+1));
        }
        m_ref_begin_bot_bins[read_begin_s]++;    
      }
      if (j<target_length(tid)-3) {
        std::string read_begin_s;
        for (int i=1; i<=3; i++) {
          read_begin_s += reference_base_char_1(tid, j+i-1);
        }        
        m_ref_begin_top_bins[read_begin_s]++;    
      }
    }
    
    m_read_begin_output << "base_1\tbase_2\tbase_3\tread_top\tread_bot\tref_top\tref_bot" << std::endl;
    for (int b1=0; b1<base_list_size-1; b1++) {
      for (int b2=0; b2<base_list_size-1; b2++) {
        for (int b3=0; b3<base_list_size-1; b3++) {
          
          std::string key_s;
          key_s += base_char_list[b1];
          key_s += base_char_list[b2];
          key_s += base_char_list[b3];
          
          uint32_t read_begin_top_count=0;
          if (m_read_begin_top_bins.find(key_s) != m_read_begin_top_bins.end()) {
            read_begin_top_count = m_read_begin_top_bins[key_s]; 
          }
          
          uint32_t read_begin_bot_count=0;
          if (m_read_begin_bot_bins.find(key_s) != m_read_begin_bot_bins.end()) {
            read_begin_bot_count = m_read_begin_bot_bins[key_s]; 
          }
          
          uint32_t ref_begin_top_count=0;
          if (m_ref_begin_top_bins.find(key_s) != m_ref_begin_top_bins.end()) {
            ref_begin_top_count = m_ref_begin_top_bins[key_s]; 
          }
          
          uint32_t ref_begin_bot_count=0;
          if (m_ref_begin_bot_bins.find(key_s) != m_ref_begin_bot_bins.end()) {
            ref_begin_bot_count = m_ref_begin_bot_bins[key_s]; 
          }
          
          m_read_begin_output << base_char_list[b1] << "\t" << base_char_list[b2] << "\t" << base_char_list[b3]
          << "\t" << read_begin_top_count << "\t" << read_begin_bot_count << "\t" << ref_begin_top_count << "\t" << ref_begin_bot_count<< std::endl;
        }
      }
    }
  }
}
  
} // namespace breseq




