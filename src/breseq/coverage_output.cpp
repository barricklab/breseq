/*****************************************************************************

 AUTHORS

   Jeffrey E. Barrick <jeffrey.e.barrick@gmail.com> and other contributors

 LICENSE AND COPYRIGHT

   Copyright (c) 2008-2010 Michigan State University
   Copyright (c) 2011-2025 The University of Texas at Austin
   Copyright (c) 2025-     Michigan State University

   breseq is free software; you can redistribute it and/or modify it under the
   terms of the GNU General Public License as published by the Free Software
   Foundation; either version 2, or (at your option) any later version.

   SPDX-License-Identifier: GPL-2.0-or-later

*****************************************************************************/

#include "coverage_output.h"

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

  double fixed_coverage_scale_value = m_reference_average_coverage
    ? m_reference_average_coverage * m_fixed_coverage_scale
    : m_fixed_coverage_scale;

  // Expressions for the "total" series: a single raw column when the table
  // only has totals (m_total_only), otherwise the sum of the top/bottom
  // strand columns.
  string unique_tot_expr    = m_total_only ? "(column(\"unique_cov\"))"     : "(column(\"unique_top_cov\")+column(\"unique_bot_cov\"))";
  string redundant_tot_expr = m_total_only ? "(column(\"redundant_cov\"))" : "(column(\"redundant_top_cov\")+column(\"redundant_bot_cov\"))";

  ostringstream s;
  s << "set datafile columnheaders" << endl;
  if (m_output_format == "pdf") {
    s << "set terminal pdfcairo size 11in,6in font ',12'" << endl;
  } else if (m_output_format == "svg") {
    s << "set terminal svg size 2200,1200 font ',28'" << endl;
  } else {
    s << "set terminal pngcairo size 2200,1200 font ',28'" << endl;
  }
  s << "set output " << double_quote(_output_file_name) << endl;
  s << "set tics out" << endl;
  s << "set border lw 2" << endl;
  s << "set xlabel 'Coordinate in Reference Genome'" << endl;
  s << "set ylabel 'Read Coverage Depth'" << endl;
  s << "set format x '%.0f'" << endl;
  s << "set xrange [" << extended_start << ":" << extended_end << "]" << endl;

  // Compute the y-axis maximum the same way plot_coverage.r did, but let
  // gnuplot's own stats command find the peak instead of scanning the table
  // in C++.
  s << "stats " << double_quote(tmp_coverage) << " using (" << unique_tot_expr << " > " << redundant_tot_expr
    << " ? " << unique_tot_expr << " : " << redundant_tot_expr << ") nooutput name 'COV'" << endl;
  s << "maxy = COV_max + 5" << endl;
  s << "if (" << to_string<double>(fixed_coverage_scale_value) << " != 0) {" << endl;
  s << "  maxy = " << to_string<double>(fixed_coverage_scale_value) << endl;
  s << "} else {" << endl;
  s << "  if (" << to_string<double>(m_reference_average_coverage) << "*1.1 > maxy) {" << endl;
  s << "    maxy = " << to_string<double>(m_reference_average_coverage) << "*1.1" << endl;
  s << "  }" << endl;
  s << "}" << endl;
  s << "set yrange [0:maxy]" << endl;

  // Grey out the flanking region outside the originally requested window,
  // drawn behind all plotted data via gnuplot's "behind" layer.
  int obj_id = 1;
  if (extended_start < static_cast<int32_t>(start_pos)) {
    s << "set object " << (obj_id++) << " rect from " << extended_start << ", graph 0 to " << start_pos
      << ", graph 1 fc rgb 'grey85' fillstyle solid 1.0 noborder behind" << endl;
  }
  if (static_cast<int32_t>(end_pos) < extended_end) {
    s << "set object " << (obj_id++) << " rect from " << (end_pos + 1) << ", graph 0 to " << extended_end
      << ", graph 1 fc rgb 'grey85' fillstyle solid 1.0 noborder behind" << endl;
  }

  // Spread out horizontally in a single row, with a filled square as the
  // sample swatch placed to the left of each label ("reverse") instead of
  // gnuplot's default of the swatch trailing the text (set below via dummy
  // NaN "with boxes" plot clauses, since "with steps"/"with lines" series
  // would otherwise show a line segment as their swatch). The slightly
  // enlarged bottom margin gives the key a bit of its own space below the
  // xlabel instead of crowding it (gnuplot's auto-sized margin otherwise
  // places the two right next to each other).
  string key_font = (m_output_format == "pdf") ? "',11'" : "',20'";
  s << "set bmargin 6" << endl;
  s << "set key below horizontal reverse Left font " << key_font << " width 4 samplen 1.5" << endl;

  // Legend swatch (filled square) clauses and the actual data clauses
  // (notitle, since the swatch above already provides the key entry) are
  // built in parallel and emitted as: all swatches, then all data -- this
  // keeps the key in the same left-to-right order as before without
  // affecting the draw order/z-stacking of the real data among itself.
  vector<string> legend_clauses;
  vector<string> data_clauses;
  string quoted_tmp_coverage = double_quote(tmp_coverage);

  if (m_reference_average_coverage != 0) {
    legend_clauses.push_back("NaN with boxes fc rgb 'dark-grey' fs solid title 'unique average'");
    data_clauses.push_back(to_string<double>(m_reference_average_coverage) + " with lines lc rgb 'dark-grey' lw 8 notitle");
  }
  if (m_total_only) {
    legend_clauses.push_back("NaN with boxes fc rgb 'green' fs solid title 'total'");
    data_clauses.push_back(quoted_tmp_coverage + " using \"position\":\"total_cov\" with steps lc rgb 'green' lw 8 notitle");
  }
  legend_clauses.push_back("NaN with boxes fc rgb 'red' fs solid title 'repeat total'");
  data_clauses.push_back(quoted_tmp_coverage + " using \"position\":" + redundant_tot_expr + " with steps lc rgb 'red' lw 3 notitle");
  if (!m_total_only) {
    legend_clauses.push_back("NaN with boxes fc rgb 'yellow' fs solid title 'repeat top'");
    data_clauses.push_back(quoted_tmp_coverage + " using \"position\":\"redundant_top_cov\" with steps lc rgb 'yellow' lw 1.4 notitle");
    legend_clauses.push_back("NaN with boxes fc rgb 'orange' fs solid title 'repeat bottom'");
    data_clauses.push_back(quoted_tmp_coverage + " using \"position\":\"redundant_bot_cov\" with steps lc rgb 'orange' lw 1.4 notitle");
  }
  legend_clauses.push_back("NaN with boxes fc rgb 'blue' fs solid title 'unique total'");
  data_clauses.push_back(quoted_tmp_coverage + " using \"position\":" + unique_tot_expr + " with steps lc rgb 'blue' lw 3 notitle");
  if (!m_total_only) {
    legend_clauses.push_back("NaN with boxes fc rgb 'cyan' fs solid title 'unique top'");
    data_clauses.push_back(quoted_tmp_coverage + " using \"position\":\"unique_top_cov\" with steps lc rgb 'cyan' lw 1.4 notitle");
    legend_clauses.push_back("NaN with boxes fc rgb 'purple' fs solid title 'unique bottom'");
    data_clauses.push_back(quoted_tmp_coverage + " using \"position\":\"unique_bot_cov\" with steps lc rgb 'purple' lw 1.4 notitle");
  }

  vector<string> plot_clauses(legend_clauses);
  plot_clauses.insert(plot_clauses.end(), data_clauses.begin(), data_clauses.end());

  s << "plot " << join(plot_clauses, string(", \\\n     ")) << endl;

  string gnuplot_script_name = m_intermediate_path + "/" + to_string(pid) + "_" + to_string(m_thread_id) + ".coverage.gp";
  string log_file_name = m_intermediate_path + "/" + to_string(pid) + "_" + to_string(m_thread_id) + ".coverage.gp.log";
  run_gnuplot_script(s.str(), gnuplot_script_name, log_file_name);

  if (m_output_format == "svg") make_svg_responsive(_output_file_name);

  remove(tmp_coverage.c_str());
  remove(log_file_name.c_str());
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

  // Select the field separator based on the output format (comma for csv, tab otherwise)
  m_table_delimiter = (m_output_format == "csv") ? "," : "\t";
  const string& d = m_table_delimiter;

  if (m_total_only) {
    m_output_table << "position" << d << "ref_base" << d
      << "unique_cov" << d << "redundant_cov" << d << "total_cov"
      << std::endl;
  } else {
    m_output_table << "position" << d << "ref_base" << d
      << "unique_top_cov" << d << "unique_bot_cov" << d
      << "redundant_top_cov" << d << "redundant_bot_cov" << d
      << "raw_redundant_top_cov" << d << "raw_redundant_bot_cov" << d
      << "unique_top_begin" << d << "unique_bot_begin"
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
    m_output_table << "#" << d << "reference_unique_average_cov" << d << m_reference_average_coverage << std::endl;
  }

  m_output_table << "#" << d << "region_unique_average_cov" << d << ( m_region_average_unique_coverage / m_region_num_positions) << std::endl;
  m_output_table << "#" << d << "region_repeat_average_cov" << d << ( m_region_average_repeat_coverage / m_region_num_positions) << std::endl;
  m_output_table << "#" << d << "region_average_cov" << d << ( m_region_average_coverage / m_region_num_positions) << std::endl;
  m_output_table << "#" << d << "number_of_positions" << d << m_region_num_positions << std::endl;
  
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
  const string& d = m_table_delimiter;
  for (uint32_t i=m_last_position_1+1; i<pos; i++) {
    if (m_total_only) {
      m_output_table << i << d << refseq[i-1] << d << 0 << d << 0 << d << 0
      << std::endl;
    } else {
      m_output_table << i << d << refseq[i-1] << d << 0 << d << 0 << d
        << 0 << d << 0 << d << 0 << d << 0 << d << 0 << d << 0 << std::endl;
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
    m_output_table << pos << d << ref_base << d
      << (unique_cov[0] + unique_cov[1]) << d
      << (redundant_cov[0] + redundant_cov[1]) << d
      << (unique_cov[0] + unique_cov[1] + redundant_cov[0] + redundant_cov[1])
      << std::endl;
  } else {
    m_output_table << pos << d << ref_base << d
      << unique_cov[0] << d << unique_cov[1] << d
      << redundant_cov[0] << d << redundant_cov[1] << d
      << raw_redundant_cov[0] << d << raw_redundant_cov[1] << d
      << unique_begin_reads[0] << d << unique_begin_reads[1]
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




