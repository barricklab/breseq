/*****************************************************************************
 * 
 * AUTHORS
 * 
 *  Jeffrey E. Barrick <jeffrey.e.barrick@gmail.com>
 *  David B. Knoester
 * 
 * LICENSE AND COPYRIGHT
 * 
 *  Copyright (c) 2008-2010 Michigan State University
 *  Copyright (c) 2011 The University of Texas at Austin
 * 
 *  breseq is free software; you can redistribute it and/or modify it under the
 *  terms the GNU General Public License as published by the Free Software
 *  Foundation; either version 1, or (at your option) any later version.
 * 
 *****************************************************************************/

#ifndef _BRESEQ_COVERAGE_OUTPUT_H_
#define _BRESEQ_COVERAGE_OUTPUT_H_

#include "common.h"
#include "pileup_base.h"
#include "pileup.h"

using namespace std;

namespace breseq
{
  
  /*! This class is a FACTORY for generating coverage plots
   */
  class coverage_output : pileup_base
  {
  protected:
    string    m_output_format; 
    uint32_t  m_downsample;
    bool      m_total_only;
    uint32_t  m_shaded_flanking;
    
    string    m_r_script_file_name;
    string    m_intermediate_path;
    string    m_read_begin_output_file_name; // extra output file set as option
    string    m_gc_output_file_name;         // extra output file set as option
    
    // variables used during tabulate pileup to collect and print statistics
    ofstream  m_output_table;
    ofstream  m_read_begin_output; 
    ofstream  m_gc_output;
    
    std::map<std::string,uint32_t> m_read_begin_top_bins;
    std::map<std::string,uint32_t> m_read_begin_bot_bins;
    std::map<std::string,uint32_t> m_ref_begin_top_bins;
    std::map<std::string,uint32_t> m_ref_begin_bot_bins;
    
    std::vector<uint32_t> m_gc_content_bins;
    
    //! Clear saved statistics before beginning a new table
    void clear();
    
    //! Called for each alignment.
		virtual void pileup_callback(const pileup& p);
    
		//! Called at end of fragment.
    void at_target_end(const uint32_t tid);
    
  public:
    
    coverage_output( const string& bam, const string& fasta, const string& r_script_file_name, const string& intermediate_path = "/tmp" )
      : pileup_base(bam, fasta), m_output_format("png"), m_downsample(0), m_total_only(false)
      , m_shaded_flanking(0), m_r_script_file_name(r_script_file_name), m_intermediate_path(intermediate_path) {};
    
    string output_format(const string& _output_format = "") 
    { 
      if (_output_format.length() > 0)
      {
        m_output_format = _output_format;
        _assert( (m_output_format=="png") || (m_output_format=="pdf"), 
                "Unrecognized coverage plot output format: " + m_output_format);
      }
      return m_output_format; 
    }
    bool shaded_flanking(uint32_t _shaded_flanking) { m_shaded_flanking = _shaded_flanking; return m_shaded_flanking; }
    
    void plot(const string& region, const string& output_file_name, uint32_t resolution = 600);
    void tabulate(const string& region, const string& output_file_name, uint32_t downsample = 1);
  };
  
  
}//end namespace breseq
#endif
