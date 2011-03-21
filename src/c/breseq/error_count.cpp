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

#include <iostream>
#include <cmath>
#include <map>
#include <sstream>
#include <string>
#include <fstream>
#include <vector>
#include <sam.h>
#include <faidx.h>
#include <assert.h>
#include <boost/tuple/tuple.hpp>

#include "breseq/common.h"
#include "breseq/pileup.h"
#include "breseq/error_count.h"

namespace breseq {

const std::string cErrorTable::covariate_names[] = {
  "read_set", 
  "ref_base", 
  "prev_ref_base", 
  "obs_base", 
  "quality", 
  "read_pos", 
  "base_repeat"
}; 
const char cErrorTable::m_sep = '\t';

/*! Count errors.
 */
void error_count(const std::string& bam, 
												 const std::string& fasta,
												 const std::string& output_dir,
												 const std::vector<std::string>& readfiles,
                         bool do_coverage,
                         bool do_errors,
                         uint8_t min_qual_score,
                         const std::string& covariates
                      ) 
{
	error_count_pileup ecp(bam, fasta, do_coverage, do_errors, min_qual_score, covariates);
	ecp.do_pileup();
	if (do_coverage) ecp.print_coverage(output_dir);
	if (do_errors) ecp.print_error(output_dir, readfiles);
}


/*! Constructor.
 */
error_count_pileup::error_count_pileup(const std::string& bam, const std::string& fasta, bool do_coverage, bool do_errors, uint8_t min_qual_score, const std::string& covariates)
: pileup_base(bam, fasta), m_do_coverage(do_coverage), m_do_errors(do_errors), m_min_qual_score(min_qual_score), m_error_table(covariates) {
	// reserve enough space for the sequence info:
	_seq_info.resize(_bam->header->n_targets);
  m_use_CErrorTable = (covariates.length() > 0);
}


/*! Destructor.
 */
error_count_pileup::~error_count_pileup() {
}


/*! Called for each alignment.
 */
void error_count_pileup::callback(const pileup& p) {
  
	using namespace std;
	assert(p.target() < _seq_info.size());
	sequence_info& info=_seq_info[p.target()];

	size_t unique_coverage=0; // number of non-deletion, non-redundant alignments at this position.
	bool has_redundant_reads=false; // flag indicating whether this position has any redundant reads.

	// for each alignment within this pileup:
	for(pileup::const_iterator i=p.begin(); i!=p.end(); ++i) {
		
    // skip deletions entirely, they are handled by adjacent matching positions
    if(i->is_del()) {
      continue;
    }
		
		// is this a redundant read?  if so, don't process it - we're all done.
		// also, mark this position as having a redundant read so that we don't update the
		// coverage count when we're done looking at all the alignments.
		if(i->is_redundant()) {
			has_redundant_reads = true;
			continue;
		}
		
		// track the number of non-deletion, non-redundant alignments:
		++unique_coverage;
		
		// if we are only tracking coverage, go to the next alignment
    if (!m_do_errors) {
      continue;
    }		

    // Newer implementation...
    if (m_use_CErrorTable) {
      m_error_table.count_alignment_position(*i, p);
      continue;
    }

				
		uint32_t reversed = i->reversed(); // are we on the reverse strand?
		uint8_t* qseq = i->query_bam_sequence(); // query sequence (read)
		int32_t qpos = i->query_position_0(); // position of the alignment in the query
		int32_t qstart = i->query_start_0(); // @dk: want 0-indexed, so subtract 1
		int32_t qend = i->query_end_0(); // @dk: want 0-indexed, so subtract 1

		uint8_t* qscore = i->quality_scores(); // quality score array
		int32_t fastq_file_index = i->fastq_file_index(); // sequencer-generated read file that this alignment belongs to
		
		uint32_t pos = p.position_0(); // position of this alignment on the reference sequence
		char* refseq = p.reference_sequence(); // reference sequence for this target
		
		// Things to remember in the following:
    // -->1 Reverse the base when the read is on the other strand
    // -->2 Correct for which base is AFTER in error counts when on the other strand
    // -->3 Observations that involve an N base in any way should not be counted

			//# (1) base substitutions
			//#     e.g. 'AG' key for observing a G in read at a place where the reference has an A
			//#     IMPROVE by keeping track of flanking base context and scores?
			//#     this would, for example, penalize low scoring sequences more

    //We already skipped deletions of *this* base so we know we have a (mis)match
    {
      base_bam read_base_bam = bam1_seqi(qseq,qpos);    
      base_char ref_base_char = refseq[pos];
      
      if (!_base_bam_is_N(read_base_bam) && !_base_char_is_N(ref_base_char)) {
        if(reversed) {
          read_base_bam = complement_base_bam(read_base_bam);
          ref_base_char = complement_base_char(ref_base_char);
        }
        string key; key += static_cast<char>(ref_base_char); key += basebam2char(read_base_bam);
        uint8_t quality = qscore[qpos];
        ++_error_hash[fastq_file_index][quality][key];
      }
    }

		//# the next base also matches 
    //# (1) base substitution or match
    //#     e.g. '..' key indicating an observation of a "non-gap, non-gap"
    //#     quality score is of the second non-gap in the pair
    if(i->indel() == 0) {
      //## don't count past last match position
      if (qpos < qend) {	
        int32_t mqpos = qpos + 1 - reversed;
        uint8_t base = bam1_seqi(qseq,mqpos);
        
        int32_t mrpos = pos + 1 - reversed;
        uint8_t ref_base = refseq[mrpos];
        
        if (!_base_bam_is_N(base) && !_base_char_is_N(ref_base)) {     
          string key = ".."; 
          uint8_t quality = qscore[mqpos];
          ++_error_hash[fastq_file_index][quality][key];
        }
      }	
    }
		
		//# there is a deletion of EXACTLY one base in the read relative to the reference before the next read base
    //# (2) deletion in read relative to reference
    //#     e.g. 'A.' key for observing nothing in a read at a position where the reference has an A
    //#     quality score is of the next non-gap base in the read
    else if (i->indel() == -1) {
      //## count the quality of this or next base depending on reversed, and make sure it is not an N
      int32_t mqpos = qpos + 1 - reversed;
      base_bam read_base_bam = bam1_seqi(qseq,mqpos);   
			
      //## the reference base opposite the deletion is really the NEXT base
      int32_t mrpos = pos + 1;
      base_char ref_base_char = refseq[mrpos];      
      
      if (!_base_bam_is_N(read_base_bam) && !_base_char_is_N(ref_base_char)) {
        if(reversed) {
          ref_base_char = complement_base_char(ref_base_char);
        }
        string key; key += static_cast<char>(ref_base_char); key += '.';
        uint8_t quality = qscore[mqpos];
        ++_error_hash[fastq_file_index][quality][key];									
      }
    }

		
    //# there is an insertion of EXACTLY one base in the read relative to the reference before the next reference base
    //# (3) insertion in read relative to reference
    //#     e.g. '.A' key for observing an A in a read at a position where the reference has no base
    //#     quality score is that of the observed inserted base
    else if (i->indel() == +1) {
      int32_t mqpos = qpos + 1;
			
      if ((mqpos <= qend) && (mqpos >= qstart)) {
        uint8_t read_base_bam = bam1_seqi(qseq,mqpos);    
				
        if (!_base_bam_is_N(read_base_bam)) {
          if(reversed) {
            read_base_bam = complement_base_bam(read_base_bam);
          }        
          string key; key += '.'; key += basebam2char(read_base_bam);
          uint8_t quality = qscore[mqpos];
          ++_error_hash[fastq_file_index][quality][key];									
        }	
      }
    }
	}
	
	// NOTE: this can't move inside the for-loop; we're tracking information about
	// the position, not each alignment.
	// record coverage at this position, but only if there were no redundant reads:
	if(!has_redundant_reads) {
		// resize our coverage table as needed:
		if(unique_coverage >= info.unique_only_coverage.size()) {
			info.unique_only_coverage.resize(unique_coverage+1,0); // >= and +1 because of 0-indexing.
		}
		++info.unique_only_coverage[unique_coverage];
	}
}


/*! Print coverage distribution.
 */
void error_count_pileup::print_coverage(const std::string& output_dir) {
	using namespace std;
	for(std::size_t i=0; i<_seq_info.size(); ++i) {
		string filename(output_dir + _bam->header->target_name[i] + ".unique_only_coverage_distribution.tab");
		ofstream out(filename.c_str());					
		
		out << "coverage\tn" << endl;
		for(std::size_t j=1; j<_seq_info[i].unique_only_coverage.size(); ++j) {
			out << j << "\t" << _seq_info[i].unique_only_coverage[j] << endl;
		}	
		out.close();
	}
}


/*! Print error file.
 */
void error_count_pileup::print_error(const std::string& output_dir, const std::vector<std::string>& readfiles) {

  if (m_use_CErrorTable) {
      std::string output_file;
      output_file = output_dir + "error_counts.tab"; 
      m_error_table.write_count_table(output_file);
      
      output_file = output_dir + "error_rates.tab"; 
      m_error_table.counts_to_log10_prob();
      m_error_table.write_log10_prob_table(output_file);
      return;
  }
    
	using namespace std;
	char bases[] = {'A', 'T', 'C', 'G', '.'};
	
	assert(readfiles.size() == _error_hash.size());
	
	for(fastq_map_t::iterator iter=_error_hash.begin(); iter!=_error_hash.end(); ++iter) {
		ostringstream filename;
		filename << output_dir << readfiles[iter->first] << ".error_counts.tab";
		ofstream out(filename.str().c_str());
		
		out << "quality";
		for(int i=0; i<5; ++i) {
			for(int j=0; j<5; ++j) {
				out << "\t" << bases[i] << bases[j];
			}
		}
		out << endl;
		
		qual_map_t& qual_map=iter->second;
		for(qual_map_t::reverse_iterator iter=qual_map.rbegin(); iter!=qual_map.rend(); ++iter) {
			out << static_cast<unsigned int>(iter->first);
			for(int i=0; i<5; ++i) {
				for(int j=0; j<5; ++j) {
					string k; k += bases[i]; k += bases[j];
					out << "\t" << iter->second[k];
				}
			}
			out << endl;
		}
		out.close();
	}
}


/*! Load error rates.
 */
error_count_results::error_count_results(const std::string& input_dir, const std::vector<std::string>& readfiles) {
	using namespace std;
	char bases[] = {'A', 'T', 'C', 'G', '.'}; // order is important!!! must match the header from the error rates file...
	
	_error_rates.resize(readfiles.size());
	_log10_error_rates.resize(readfiles.size());

	// load the error rates:
	for(std::size_t i=0; i<readfiles.size(); ++i) {
		string filename(input_dir + readfiles[i] + ".error_rates.tab");
		ifstream in(filename.c_str());
		in.ignore(1024, '\n'); // get rid of header line.
		
		error_map_t& emt = _error_rates[i];
		error_map_t& log10_emt = _log10_error_rates[i];

		while(!in.eof()) {
			int quality;
			in >> quality;
			
      // @JEB: If we don't ignore trailing line after line ending of last line, 
      // then we get bogus values written over the lowest quality score!
      if (!in) break;
      
			error_map_t::iterator em = emt.insert(make_pair(static_cast<uint8_t>(quality),base_error_t())).first;
			error_map_t::iterator log10_em = log10_emt.insert(make_pair(static_cast<uint8_t>(quality),base_error_t())).first;

			for(int i=0; i<5; ++i) {
				for(int j=0; j<5; ++j) {
					string k; k += bases[i]; k += bases[j];
					double r;
					in >> r;
          if (r == 0) r = 1e-256;
          double one_minus_pr = 1.0-r;
          if (one_minus_pr == 0) one_minus_pr = 1e-256;
          
          double cor = r;
					double err = one_minus_pr;
          em->second[k] = make_pair(err,cor);
          
					double log10_cor = log10(r);
					double log10_err = log10(one_minus_pr);
					log10_em->second[k] = make_pair(log10_err,log10_cor);
          
          //std::cerr << r << std::endl;
				}
			}
		}
		in.close();
	}
}

/*! Return the log10 error rate for the given base pair, quality, and FASTQ file index.
 */
double error_count_results::log10_error_rates(int32_t fastq_file_index, uint8_t quality, const std::string& base_key) {
	//  fail if this doesn't exist?
	assert((std::size_t)fastq_file_index<_log10_error_rates.size());
	assert(_log10_error_rates[fastq_file_index].find(quality)!=_log10_error_rates[fastq_file_index].end());
	assert(_log10_error_rates[fastq_file_index][quality].find(base_key)!=_log10_error_rates[fastq_file_index][quality].end());
	
	return _log10_error_rates[fastq_file_index][quality][base_key].first;
}


/*! Return the log10 correct rate for the given base pair, quality, and FASTQ file index.
 */
double error_count_results::log10_correct_rates(int32_t fastq_file_index, uint8_t quality, const std::string& base_key) {
		//  fail if this doesn't exist?
	assert((std::size_t)fastq_file_index<_log10_error_rates.size());
	assert(_log10_error_rates[fastq_file_index].find(quality)!=_log10_error_rates[fastq_file_index].end());
	assert(_log10_error_rates[fastq_file_index][quality].find(base_key)!=_log10_error_rates[fastq_file_index][quality].end());
	
	return _log10_error_rates[fastq_file_index][quality][base_key].second;
}

/*! Return the log10 correct rate for the given base pair, quality, and FASTQ file index.
 */
double error_count_results::correct_rates(int32_t fastq_file_index, uint8_t quality, const std::string& base_key) {
		//  fail if this doesn't exist?
	assert((std::size_t)fastq_file_index<_error_rates.size());
	assert(_error_rates[fastq_file_index].find(quality)!=_error_rates[fastq_file_index].end());
	assert(_error_rates[fastq_file_index][quality].find(base_key)!=_error_rates[fastq_file_index][quality].end());
	
	return _error_rates[fastq_file_index][quality][base_key].second;
}

/*! Return the pair of (correct,error) rates.
 */
const std::pair<double,double>& error_count_results::log10_rates(int32_t fastq_file_index, uint8_t quality, const std::string& base_key) {
	return _log10_error_rates[fastq_file_index][quality][base_key];
}

/*  cErrorTable::cErrorTable()

    This version of the initializer is for constructing the full
    table that will be used from command line arguments.
*/
cErrorTable::cErrorTable (const std::string& colnames)
{
  if (colnames.length() == 0) return;
  read_covariates(colnames);
  allocate_table();
}

/*  cErrorTable::cErrorTable()

    Allocates an empty table with the given covariates.
*/
cErrorTable::cErrorTable(covariates_used_t covariate_used, covariates_max_t covariate_max, covariates_enforce_max_t covariate_enforce_max) {

  int current_offset = 1;

  for (int i=0; i<k_num_covariates; i++) {
  
    m_covariate_used[i] = covariate_used[i];
    m_covariate_max[i] = covariate_max[i];
    m_covariate_enforce_max[i] = covariate_enforce_max[i];
    
    m_covariate_offset[i] = current_offset;
    
    if (m_covariate_used[i]) {
        current_offset *= m_covariate_max[i];
    }
  }
  
  allocate_table();
}

/*  cErrorTable::cErrorTable()

    This version of the initializer is for constructing sub-tables
    that sum across various dimensions of the full table.
    
    The new list of covariates should be a subset of the ones used by 
    the error_table that is provided.
*/
cErrorTable::cErrorTable(cErrorTable& error_table, covariates_used_t covariates) {

  int current_offset = 1;

  // Set up new covariates
  for (int i=0; i<k_num_covariates; i++) {
    if (covariates[i]) {
      m_covariate_used[i] = true;
      m_covariate_max[i] = error_table.m_covariate_max[i];
      m_covariate_enforce_max[i] = error_table.m_covariate_enforce_max[i];    
      m_covariate_offset[i] = current_offset;
    } else {
      m_covariate_used[i] = false;
      m_covariate_max[i] = 0;
      m_covariate_enforce_max[i] = false;
      m_covariate_offset[i] = 0;
    }
    
    if (m_covariate_used[i]) {
        current_offset *= m_covariate_max[i];
    }
  }
  
  allocate_table();

  // Assign old data bins to new ones, combining different bin totals
  for (uint32_t i=0; i<error_table.m_count_table.size(); i++) {
    covariate_values_t cv;
    error_table.index_to_covariates(i, cv);
    uint32_t j = covariates_to_index(cv);
    m_count_table[j] += error_table.m_count_table[i];
  }
}

/*  cErrorTable::allocate_table()

    Resize the count table.
*/
void cErrorTable::allocate_table() {
  int n = 1;
  for (int i=0; i<k_num_covariates; i++) {
    if (m_covariate_used[i]) n *= m_covariate_max[i];
  }
  m_count_table.resize(n);
}

/*  cErrorTable::covariates_to_index()

    Find which line of the table corresponds to the given covariates.
*/

uint32_t cErrorTable::covariates_to_index(const covariate_values_t& cv) {

  // Calculate the row index in which to record the observation
  uint32_t idx = 0;
  for (int i=0; i<k_num_covariates; i++) {
    if (m_covariate_used[i]) {
    
      uint32_t val = cv[i];
      if(val >= m_covariate_max[i]) {
        if (m_covariate_enforce_max[i] == false) {
          assert(false);
        } else {
          std::cout << "Adjusted down " << val << " to " << m_covariate_max[i]-1 << std::endl;
          val = m_covariate_max[i]-1;
        }
      }
      idx += val * m_covariate_offset[i];
    }
  }
  return idx;
}

/*  cErrorTable::index_to_covariates()

    Fill in values of covariates, given a table line index.
    Note: will not change value if covariate is not used.
*/

void cErrorTable::index_to_covariates(const uint32_t idx, covariate_values_t& cv) {
          
  for (int i=0; i<k_num_covariates; i++) {
    if (!m_covariate_used[i]) continue;
    
    uint32_t j = (idx / m_covariate_offset[i]) % m_covariate_max[i];
   
    // set the corresponding covariate
    cv[i] = j;
  }
}


/*  cErrorTable::read_covariates()

    Convert string representation to covariates
*/
void cErrorTable::read_covariates(const std::string& colnames) {

  // set default values
  for (int i=0; i<k_num_covariates; i++) {
    m_covariate_used[i] = false;
    m_covariate_enforce_max[i] = false;
  }

  std::vector<std::string> columns_to_use;
  cErrorTable::split(colnames, ',', columns_to_use);

  // Turn on and assign max
  for(std::vector<std::string>::iterator i=columns_to_use.begin(); i!=columns_to_use.end(); ++i) {
    std::vector<std::string> columns_parts;
    cErrorTable::split(*i, '=', columns_parts);

		if (columns_parts[0] == "ref_base") {
      m_covariate_used[k_ref_base] = true;
      m_covariate_max[k_ref_base] = 5;
    }
		else if (columns_parts[0] == "prev_base") {
      m_covariate_used[k_prev_base] = true;
      m_covariate_max[k_prev_base] = 5;
    }
		else if (columns_parts[0] == "obs_base") {
      m_covariate_used[k_obs_base] = true;
      m_covariate_max[k_obs_base] = 5;
    }
    else if (columns_parts[0] == "quality") {
      m_covariate_used[k_quality] = true;
      // This is a bit lazy, we could assign a minimum quality offset to save empty bins...
      m_covariate_max[k_quality] = atoi(columns_parts[1].c_str()) + 1;
    }
    else if (columns_parts[0] == "read_set") {
      m_covariate_used[k_read_set] = true;
      m_covariate_max[k_read_set] = atoi(columns_parts[1].c_str());
    }
    else if (columns_parts[0] == "read_pos") {
      m_covariate_used[k_read_pos] = true;
      m_covariate_max[k_read_pos] = atoi(columns_parts[1].c_str());
    }
    else if (columns_parts[0] == "base_repeat") {
      m_covariate_used[k_base_repeat] = true;
      m_covariate_max[k_base_repeat] = atoi(columns_parts[1].c_str())+1;
      // Base repeat will trypically have a cap, above which we combine all values
      m_covariate_enforce_max[k_base_repeat] = true;
    }
    else {
      std::cerr << "Unrecognized covariate: " << columns_parts[0] << std::endl;
      assert(1);
    }
	}
  
  // Then assign offsets, so they are in a consistent order...
  int current_offset = 1;

  for (int i=0; i<k_num_covariates; i++) {
    if (m_covariate_used[i]) {
      m_covariate_offset[i] = current_offset;
      current_offset *= m_covariate_max[i];
    }
  }
}


/*  cErrorTable::print_covariates()

    Convert covariates to string representation
*/
std::string cErrorTable::print_covariates() {

  std::stringstream covariate_string;
  
  for (uint32_t i=0; i<k_num_covariates; i++) {
    if (m_covariate_used[i]) {
      if (covariate_string.str().length() > 0) covariate_string << ",";
      
      // special cases where there are always 5 choices
      if ((i == k_ref_base) || (i == k_obs_base)) {
        covariate_string << covariate_names[i];
      }
      else {
        covariate_string << covariate_names[i];
        covariate_string << "=";
        covariate_string << m_covariate_max[i];
      }
    }
  }
  return covariate_string.str();
}


/*  cErrorTable::split()

    Split the string into a vector on each occurrence of a character.
*/
void cErrorTable::split(const std::string& s, char c, std::vector<std::string>& v) {
  std::string::size_type i = 0;
  std::string::size_type j = s.find(c);

  while (j != std::string::npos) {
      v.push_back(s.substr(i, j-i));
      i = ++j;
      j = s.find(c, j);
  }
  if (j == std::string::npos) v.push_back(s.substr(i, s.length( )));
}



/*  cErrorTable::read_log10_prob_table()

    Read table of log10 probabilities of observations
*/
void cErrorTable::read_log10_prob_table(const std::string& filename) {

  std::ifstream in(filename.c_str());
  
  std::string s;
  
  // First line contains the covariates  
  getline(in, s);
  read_covariates(s);
  allocate_table();
  m_log10_prob_table.resize(m_count_table.size());
  
  // Ignore the neader line, we expect order to be constant
  getline(in, s);
  
  // We don't need the beginning of the line, since we know the
  // order everything was printed in. That's just there to make it
  // human readable
  
  // Read one line for each
  for (uint32_t i=0; i< m_log10_prob_table.size(); i++) {
    getline(in, s);
    std::vector<std::string> split_line;
    split(s, '\t', split_line);
    m_log10_prob_table[i] = strtod(split_line.back().c_str(), NULL);
  }           
}


/*  cErrorTable::print()

    Print out a table of covariates and counts.
*/
void cErrorTable::write_log10_prob_table(const std::string& filename) {

  std::ofstream out(filename.c_str());
  
  // First line contains the covariates
  out << print_covariates() << std::endl;
  
  for (int i=0; i<k_num_covariates; i++) {
    if(m_covariate_used[i]) out << covariate_names[i] << m_sep;
  }
  out << "log10_probability" << std::endl;

  for (uint32_t idx=0; idx<m_log10_prob_table.size(); idx++) {
      
      for (int i=0; i<k_num_covariates; i++) {
        if (!m_covariate_used[i]) continue;
        
        uint32_t j = (idx / m_covariate_offset[i]) % m_covariate_max[i];
       
        if ( (i == k_ref_base) || (i == k_obs_base) ) {
          char base = baseindex2char(j);
          out << base << m_sep;
        }
        else {
          out << j << m_sep;
        }
    }
        
    out << m_log10_prob_table[idx] << std::endl;
  }
}

/*  cErrorTable::print()

    Print out a table of covariates and counts.
*/
void cErrorTable::write_count_table(const std::string& filename) {

  std::ofstream out(filename.c_str());
  
  // First line contains the covariates
  out << print_covariates() << std::endl;
  
  for (int i=0; i<k_num_covariates; i++) {
    if(m_covariate_used[i]) out << covariate_names[i] << m_sep;
  }
  out << "count" << std::endl;

  for (uint32_t idx=0; idx<m_count_table.size(); idx++) {
      
      for (int i=0; i<k_num_covariates; i++) {
        if (!m_covariate_used[i]) continue;
        
        uint32_t j = (idx / m_covariate_offset[i]) % m_covariate_max[i];
       
        if ( (i == k_ref_base) || (i == k_obs_base) ) {
          char base = baseindex2char(j);
          out << base << m_sep;
        }
        else {
          out << j << m_sep;
        }
    }
        
    out << m_count_table[idx] << std::endl;
  }
}


/*  cErrorTable::count_alignment_position()

    Record all observations in an alignment to a column of the reference genome
    in the count table. This function is called by the pileup callback.
*/
void cErrorTable::count_alignment_position(const alignment& i, const pileup& p) {
    
    uint32_t ref_pos = p.position_0();
    const char* ref_seq = p.reference_sequence();
    
		uint32_t reversed = i.reversed(); // are we on the reverse strand?
		uint8_t* qseq = i.query_bam_sequence(); // query sequence (read)
		int32_t q_pos_0 = i.query_position_0(); // 0-indexed
		int32_t q_start_0 = i.query_start_0(); // 0-indexed
		int32_t q_end_0 = i.query_end_0(); // 0-indexed

		uint8_t* qscore = i.quality_scores(); // quality score array
		int32_t fastq_file_index = i.fastq_file_index(); // sequencer-generated read file that this alignment belongs to

    // Fill in all covariates that are used...
    covariate_values_t cv;
    cv.read_set() = fastq_file_index;
    cv.read_pos() = q_pos_0;
		
		// Things to remember in the following:
    // -->1 Reverse the base when the read is on the other strand
    // -->2 Correct for which base is AFTER in error counts when on the other strand
    // -->3 Observations that involve an N base in any way should not be counted

			//# (1) base substitutions
			//#     e.g. 'AG' key for observing a G in read at a place where the reference has an A
			//#     IMPROVE by keeping track of flanking base context and scores?
			//#     this would, for example, penalize low scoring sequences more
      
    //We already skipped deletions of *this* base so we know we have a (mis)match
    {
      uint8_t obs_base_bam = bam1_seqi(qseq,q_pos_0);    
      uint8_t ref_base_char = ref_seq[ref_pos];
      
      if (!_base_bam_is_N(obs_base_bam) && !_base_char_is_N(ref_base_char)) {
        if(reversed) {
          obs_base_bam =complement_base_bam(obs_base_bam);
          ref_base_char = complement_base_char(ref_base_char);
        }
        
        cv.quality() = qscore[q_pos_0];
        cv.obs_base() = basebam2index(obs_base_bam);
        cv.ref_base() = basechar2index(ref_base_char);
        if (m_covariate_used[k_base_repeat]) cv.base_repeat() = i.base_repeat_0(q_pos_0);
        count_covariate(cv);
      }
    }

		//# the next base also matches 
    //# (1) base substitution or match
    //#     e.g. '..' key indicating an observation of a "non-gap, non-gap"
    //#     quality score is of the second non-gap in the pair
    if(i.indel() == 0) {
      //## don't count past last match position
      if (q_pos_0 < q_end_0) {	
        int32_t mqpos = q_pos_0 + 1 - reversed;
        base_bam obs_base_bam = bam1_seqi(qseq,mqpos);
        
        int32_t mrpos = ref_pos + 1 - reversed;
        base_char ref_base_char = ref_seq[mrpos];
        
        if (!_base_bam_is_N(obs_base_bam) && !_base_char_is_N(ref_base_char)) {     
          cv.quality() = qscore[mqpos];
          cv.obs_base() = basechar2index('.');
          cv.ref_base() = basechar2index('.');
          if (m_covariate_used[k_base_repeat]) cv.base_repeat() = i.base_repeat_0(mqpos);
          count_covariate(cv);
        }
      }	
    }
		
		//# there is a deletion of EXACTLY one base in the read relative to the reference before the next read base
    //# (2) deletion in read relative to reference
    //#     e.g. 'A.' key for observing nothing in a read at a position where the reference has an A
    //#     quality score is of the next non-gap base in the read
    else if (i.indel() == -1) {
      //## count the quality of this or next base depending on reversed, and make sure it is not an N
      int32_t mqpos = q_pos_0 + 1 - reversed;
      base_char obs_base_bam = bam1_seqi(qseq,mqpos);   
			
      //## the reference base opposite the deletion is really the NEXT base
      int32_t mrpos = ref_pos + 1;
      base_char ref_base_char = ref_seq[mrpos];      
      
      if (!_base_bam_is_N(obs_base_bam) && !_base_char_is_N(ref_base_char)) {
        
        if(reversed) ref_base_char = complement_base_char(ref_base_char);
        
        cv.quality() = qscore[mqpos];
        cv.obs_base() = basechar2index('.');
        cv.ref_base() = basechar2index(ref_base_char);
        if (m_covariate_used[k_base_repeat]) cv.base_repeat() = i.base_repeat_0(mqpos);
        count_covariate(cv);
      }
    }

		
    //# there is an insertion of EXACTLY one base in the read relative to the reference before the next reference base
    //# (3) insertion in read relative to reference
    //#     e.g. '.A' key for observing an A in a read at a position where the reference has no base
    //#     quality score is that of the observed inserted base
    else if (i.indel() == +1) {
      int32_t mqpos = q_pos_0 + 1;
			
      if ((mqpos <= q_end_0) && (mqpos >= q_start_0)) {
        base_bam obs_base_bam = bam1_seqi(qseq,mqpos);    
				
        if (!_base_bam_is_N(obs_base_bam)) {
        
          if(reversed) obs_base_bam = complement_base_bam(obs_base_bam);       
          
          cv.quality() = qscore[mqpos];
          cv.obs_base() = basebam2index(obs_base_bam);
          cv.ref_base() = basechar2index('.');
          
          if (m_covariate_used[k_base_repeat]) cv.base_repeat() = i.base_repeat_0(mqpos);
          count_covariate(cv);
        }	
      }
    }
    
    
    
}

/*  cErrorTable::count_covariate()

    Record an observation of the specified covariates in table.
*/

void cErrorTable::count_covariate(const covariate_values_t& cv) {

  uint32_t i=covariates_to_index(cv); 
  m_count_table[i]++;
}


/*  cErrorTable::print_empirical_error_rates()

  Calculate empirical error rates using Yates correction (i.e., adding 1 to each bin)
*/

void cErrorTable::counts_to_log10_prob() {

  const uint32_t smoothing_factor = 1;

  m_log10_prob_table.resize(m_count_table.size());

  // Creates new table that has summed across all occurences of observed base
  covariates_used_t covariates;
  memcpy(covariates, m_covariate_used, sizeof(covariates));
  covariates[k_obs_base] = false;
  cErrorTable sum_error_table(*this, covariates);
  
  // Now calculate the ratios of things to these
  for (uint32_t i=0; i<m_count_table.size(); i++) {
    covariate_values_t cv;
    index_to_covariates(i, cv);
    
    uint32_t j = sum_error_table.covariates_to_index(cv);
    m_log10_prob_table[i] = log10((double)(m_count_table[i]+smoothing_factor)) 
      - log10((double)sum_error_table.m_count_table[j]+smoothing_factor*m_covariate_max[k_obs_base]);    
  }
}

/*  cErrorTable::alignment_position_to_covariates()

  Extracts covariates from the alignment and insert index.
  
  Note: Does not fill in the ref_base covariate!
*/

bool cErrorTable::alignment_position_to_covariates(const alignment& a, int32_t insert_count, covariate_values_t& cv) {
  
  // -1 for deletion, otherwise 1-number of bases inserted at this position
  int indel=a.on_base_indel();
  base_bam read_base_bam = a.on_base_bam(insert_count);
        
  //##don't use bases without qualities!!
  if(_base_bam_is_N(read_base_bam)) return false;
  
  //## These are the start and end coordinates of the aligned part of the read
  uint32_t q_start_0,q_end_0;
  boost::tie(q_start_0,q_end_0) = a.query_bounds_0(); // @dk: 1-indexed!
 
  //## (1) Mis(match) in read relative to reference...
  //##       Quality is of the current base in the read, we have ALREADY checked that it is not an N					
  //##       Don't need to do anything.
  //  if (indel == 0)
          
  uint32_t q_pos_0 = a.query_position_0();
  uint32_t quality = 0;
  
  //## (2) Deletion in read relative to reference...
  //##       Quality is of the NEXT base in the read, and check that it is not an N
  if (indel == -1)
  {			
    q_pos_0 += 1 - a.reversed(); 
    base_bam check_base_bam = a.query_base_bam_0(q_pos_0);
    if (_base_bam_is_N(check_base_bam)) return false;
  }
  
  //## (3) Insertion in read relative to reference...
  //##       Quality is of the NEXT base in the read, and check that it is not an N
  //##       Note that it is possible this read base may be a '.' (supporting the non-insert call)
  else if (insert_count > 0) 
  {		
    // Offset as much as asked for by insert_count, but this particular read
    // may nothave inserted bases at all of those positons, so cap at indel.
    int32_t max_offset = insert_count;
    if (indel < max_offset) max_offset = indel;
    q_pos_0 += max_offset + 1 - a.reversed(); 
            
    //## Check bounds: it's possible to go past the end of the read because
    //## this is the last base of this read, but other reads have inserted bases
    if (q_pos_0 > q_end_0) return false;
    base_bam check_base_bam = a.query_base_bam_0(q_pos_0);
    if (_base_bam_is_N(check_base_bam)) return false;
  }

  //eventually include in above...
  cv.obs_base() = basebam2index(read_base_bam);
  cv.quality() = a.quality_base_0(q_pos_0);
  cv.read_set() = a.fastq_file_index();
  cv.read_pos() = q_pos_0;
  
  if (m_covariate_used[k_base_repeat]) cv.base_repeat() = a.base_repeat_0(q_pos_0);

  return true;
}

double cErrorTable::get_log10_prob(covariate_values_t& cv) {

  assert(m_log10_prob_table.size() > 0);
  uint32_t i = covariates_to_index(cv);

  assert(i < m_log10_prob_table.size());
  return m_log10_prob_table[i];
}


} //namespace breseq

