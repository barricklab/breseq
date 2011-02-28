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

#include "breseq/common.h"
#include "breseq/pileup.h"
#include "breseq/error_count.h"

namespace breseq {

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
      m_error_table.record_alignment(*i, p.position(), p.reference_sequence());
      continue;
    }

				
		uint32_t reversed = i->reversed(); // are we on the reverse strand?
		uint8_t* qseq = i->query_sequence(); // query sequence (read)
		int32_t qpos = i->query_position(); // position of the alignment in the query

		int32_t qstart = i->query_start() - 1; // @dk: want 0-indexed, so subtract 1
		int32_t qend = i->query_end() - 1; // @dk: want 0-indexed, so subtract 1

		uint8_t* qscore = i->quality_scores(); // quality score array
		int32_t fastq_file_index = i->fastq_file_index(); // sequencer-generated read file that this alignment belongs to
		
		uint32_t pos = p.position(); // position of this alignment on the reference sequence
		char* refseq = p.reference_sequence(); // reference sequence for this target
		
		// Things to remember in the following:
    // -->1 Reverse the base when the read is on the other strand
    // -->2 Correct for which base is AFTER in error counts when on the other strand
    // -->3 Observations that involve an N base in any way should not be counted

			//# (1) base substitutions
			//#     e.g. 'AG' key for observing a G in read at a place where the reference has an A
			//#     IMPROVE by keeping track of flanking base context and scores?
			//#     this would, for example, penalize low scoring sequences more
		{
      uint8_t base = bam1_seqi(qseq,qpos);    
      uint8_t ref_base = refseq[pos];
      
      if (!is_N(base) && !is_char_N(ref_base)) {
        if(reversed) {
          base = reverse_base(base);
          ref_base = reverse_base(ref_base);
        }
        string key; key += static_cast<char>(ref_base); key += base2char(base);
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
        
        if (!is_N(base) && !is_char_N(ref_base)) {     
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
      uint8_t base = bam1_seqi(qseq,mqpos);   
			
      //## the reference base opposite the deletion is really the NEXT base
      int32_t mrpos = pos + 1;
      uint8_t ref_base = refseq[mrpos];      
      
      if (!is_N(base) && !is_char_N(ref_base)) {
        if(reversed) {
          ref_base = reverse_base(ref_base);
        }
        string key; key += static_cast<char>(ref_base); key += '.';
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
        uint8_t base = bam1_seqi(qseq,mqpos);    
				
        if (!is_N(base)) {
          if(reversed) {
            base = reverse_base(base);
          }        
          string key; key += '.'; key += base2char(base);
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
      m_error_table.print(output_dir);
      m_error_table.recalibrate();
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
cErrorTable::cErrorTable (const std::string& colnames) : m_sep("\t")
{

  for (int i=0; i<k_num_covariates; i++) {
    m_covariate_used[i] = false;
  }

  std::vector<std::string> columns_to_use;
  cErrorTable::split(colnames, ',', columns_to_use);

  int current_offset = 1;

  for(std::vector<std::string>::iterator i=columns_to_use.begin(); i!=columns_to_use.end(); ++i) {
    std::vector<std::string> columns_parts;
    cErrorTable::split(*i, '=', columns_parts);

		if (columns_parts[0] == "ref") {
      m_covariate_used[k_ref_base] = true;
      m_covariate_max[k_ref_base] = 5;
      m_covariate_offset[k_ref_base] = current_offset;
      current_offset *= m_covariate_max[k_ref_base];
    }
		if (columns_parts[0] == "obs") {
      m_covariate_used[k_obs_base] = true;
      m_covariate_max[k_obs_base] = 5;
      m_covariate_offset[k_obs_base] = current_offset;
      current_offset *= m_covariate_max[k_obs_base];
    }
    if (columns_parts[0] == "quality") {
      m_covariate_used[k_quality] = true;
      m_covariate_max[k_quality] = atoi(columns_parts[1].c_str());
      m_covariate_offset[k_quality] = current_offset;
      current_offset *= m_covariate_max[k_quality];
    }
    if (columns_parts[0] == "read_set") {
      m_covariate_used[k_read_set] = true;
      m_covariate_max[k_read_set] = atoi(columns_parts[1].c_str());
      m_covariate_offset[k_read_set] = current_offset;
      current_offset *= m_covariate_max[k_read_set];
    }
    if (columns_parts[0] == "read_pos") {
      m_covariate_used[k_read_pos] = true;
      m_covariate_max[k_read_pos] = atoi(columns_parts[1].c_str());
      m_covariate_offset[k_read_pos] = current_offset;
      current_offset *= m_covariate_max[k_read_pos];
    }
    if (columns_parts[0] == "sep") {
      m_sep = columns_parts[1];
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
cErrorTable::cErrorTable(cErrorTable& error_table, covariates_used_t covariates) : m_sep("\t") {

  for (int i=0; i<k_num_covariates; i++) {
    if (covariates[i]) {
      m_covariate_used[i] = true;
      m_covariate_max[i] = error_table.m_covariate_max[i];
      m_covariate_offset[i] = error_table.m_covariate_offset[i];
    } else {
      m_covariate_used[i] = false;
      m_covariate_max[i] = 0;
      m_covariate_offset[i] = 0;
    }
  }
  
  allocate_table();
  
  char ref_base;
  char obs_base; 
  uint32_t quality; 
  uint32_t read_pos; 
  uint32_t read_set;
    
  for (uint32_t i=0; i<error_table.m_count_table.size(); i++) {
  
    error_table.index_to_covariates(i, ref_base, obs_base, quality, read_pos, read_set);
    uint32_t j = covariates_to_index(ref_base, obs_base, quality, read_pos, read_set);
    m_count_table[j] += error_table.m_count_table[i];
  }
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

/*  cErrorTable::allocate_table()

    Resize the count table.
*/
void cErrorTable::allocate_table () {
  int n = 1;
  for (int i=0; i<k_num_covariates; i++) {
    if (m_covariate_used[i]) n *= m_covariate_max[i];
  }
  m_count_table.resize(n);
}

/*  cErrorTable::print_line()

    Print out one line of the table of covariates and counts.
*/
void cErrorTable::print_line(std::ostream& out, const uint32_t idx) {
  
  for (int i=0; i<k_num_covariates; i++) {
    if (!m_covariate_used[i]) continue;
    
    uint32_t j = (idx / m_covariate_offset[i]) % m_covariate_max[i];
   
    if ( (i == k_ref_base) || (i == k_obs_base) ) {
      char base = index2basechar(j);
      out << base << m_sep;
    }
    else {
      out << j << m_sep;
    }
  }
      
  out << m_count_table[idx] << std::endl;
    
}


/*  cErrorTable::print()

    Print out a table of covariates and counts.
*/
void cErrorTable::print(const std::string& filename) {

  std::ofstream out(filename.c_str());

  std::string covariate_names[] = {"ref_base", "obs_base", "quality", "read_pos", "read_set"};
  
  for (int i=0; i<k_num_covariates; i++) {
    if(m_covariate_used[i]) out << covariate_names[i] << m_sep;
  }
  out << "count" << std::endl;

  for (uint32_t i=0; i<m_count_table.size(); i++) {
    print_line(out, i);
  }
}


/*  cErrorTable::record_observation()

    Record all observations in an alignment to a column of the reference genome
    in the count table. This function is called by the pileup callback.
*/
void cErrorTable::record_alignment(const alignment& i, const uint32_t ref_pos, const char* ref_seq) {

		uint32_t reversed = i.reversed(); // are we on the reverse strand?
		uint8_t* qseq = i.query_sequence(); // query sequence (read)
		int32_t qpos = i.query_position(); // position of the alignment in the query, 0-indexed

		int32_t qstart = i.query_start() - 1; // @dk: want 0-indexed, so subtract 1
		int32_t qend = i.query_end() - 1; // @dk: want 0-indexed, so subtract 1

		uint8_t* qscore = i.quality_scores(); // quality score array
		int32_t fastq_file_index = i.fastq_file_index(); // sequencer-generated read file that this alignment belongs to
		
		// Things to remember in the following:
    // -->1 Reverse the base when the read is on the other strand
    // -->2 Correct for which base is AFTER in error counts when on the other strand
    // -->3 Observations that involve an N base in any way should not be counted

			//# (1) base substitutions
			//#     e.g. 'AG' key for observing a G in read at a place where the reference has an A
			//#     IMPROVE by keeping track of flanking base context and scores?
			//#     this would, for example, penalize low scoring sequences more
		{
      uint8_t obs_base = bam1_seqi(qseq,qpos);    
      uint8_t ref_base = ref_seq[ref_pos];
      
      if (!is_N(obs_base) && !is_char_N(ref_base)) {
        if(reversed) {
          obs_base = reverse_base(obs_base);
          ref_base = reverse_base(ref_base);
        }
        uint8_t quality = qscore[qpos];
        record_observation(ref_base, obs_base, quality, qpos, fastq_file_index);
      }
    }

		//# the next base also matches 
    //# (1) base substitution or match
    //#     e.g. '..' key indicating an observation of a "non-gap, non-gap"
    //#     quality score is of the second non-gap in the pair
    if(i.indel() == 0) {
      //## don't count past last match position
      if (qpos < qend) {	
        int32_t mqpos = qpos + 1 - reversed;
        uint8_t obs_base = bam1_seqi(qseq,mqpos);
        
        int32_t mrpos = ref_pos + 1 - reversed;
        uint8_t ref_base = ref_seq[mrpos];
        
        if (!is_N(obs_base) && !is_char_N(ref_base)) {     
          uint8_t quality = qscore[mqpos];
          record_observation('.', '.', quality, qpos, fastq_file_index);
        }
      }	
    }
		
		//# there is a deletion of EXACTLY one base in the read relative to the reference before the next read base
    //# (2) deletion in read relative to reference
    //#     e.g. 'A.' key for observing nothing in a read at a position where the reference has an A
    //#     quality score is of the next non-gap base in the read
    else if (i.indel() == -1) {
      //## count the quality of this or next base depending on reversed, and make sure it is not an N
      int32_t mqpos = qpos + 1 - reversed;
      uint8_t obs_base = bam1_seqi(qseq,mqpos);   
			
      //## the reference base opposite the deletion is really the NEXT base
      int32_t mrpos = ref_pos + 1;
      uint8_t ref_base = ref_seq[mrpos];      
      
      if (!is_N(obs_base) && !is_char_N(ref_base)) {
        if(reversed) {
          ref_base = reverse_base(ref_base);
        }
        uint8_t quality = qscore[mqpos];
        record_observation(ref_base, '.', quality, qpos, fastq_file_index);
      }
    }

		
    //# there is an insertion of EXACTLY one base in the read relative to the reference before the next reference base
    //# (3) insertion in read relative to reference
    //#     e.g. '.A' key for observing an A in a read at a position where the reference has no base
    //#     quality score is that of the observed inserted base
    else if (i.indel() == +1) {
      int32_t mqpos = qpos + 1;
			
      if ((mqpos <= qend) && (mqpos >= qstart)) {
        uint8_t obs_base = bam1_seqi(qseq,mqpos);    
				
        if (!is_N(obs_base)) {
          if(reversed) {
            obs_base = reverse_base(obs_base);
          }        
          uint8_t quality = qscore[mqpos];
          record_observation('.', obs_base, quality, qpos,fastq_file_index);
        }	
      }
    }
}

/*  cErrorTable::record_observation()

    Record an observation of the specified covariates in table.
*/

void cErrorTable::record_observation(char ref_base, char obs_base, uint32_t quality, uint32_t read_pos, uint32_t read_set ) {

  uint32_t i=covariates_to_index(ref_base, obs_base, quality, read_pos, read_set); 
  m_count_table[i]++;
}

/*  cErrorTable::covariates_to_index()

    Find which line of the table corresponds to the given covariates.
*/

uint32_t cErrorTable::covariates_to_index(char ref_base, char obs_base, uint32_t quality, uint32_t read_pos, uint32_t read_set) {

  // Calculate the row index in which to record the observation

  uint32_t i = 0;
  if (m_covariate_used[k_read_set]) {
    assert(read_set < m_covariate_max[k_read_set]);
    i += read_set * m_covariate_offset[k_read_set];
  }
  if (m_covariate_used[k_read_pos]) {
    assert(read_pos < m_covariate_max[k_read_pos]);
    i += read_pos * m_covariate_offset[k_read_pos];
  }
  if (m_covariate_used[k_quality]) {
    assert(quality < m_covariate_max[k_quality]);
    i += quality * m_covariate_offset[k_quality];
  }
  if (m_covariate_used[k_obs_base]) {
    i += base2index(obs_base) * m_covariate_offset[k_obs_base];
  }  
  if (m_covariate_used[k_ref_base]) {
    i += base2index(ref_base) * m_covariate_offset[k_ref_base];
  }
  return i;
}

/*  cErrorTable::index_to_covariates()

    Fill in values of covariates, given a table line index.
    Note: will not change value if covariate is not used.
*/

void cErrorTable::index_to_covariates(const uint32_t idx, char& ref_base, char& obs_base, uint32_t& quality, uint32_t& read_pos, uint32_t& read_set) {
          
  for (int i=0; i<k_num_covariates; i++) {
    if (!m_covariate_used[i]) continue;
    
    uint32_t j = (idx / m_covariate_offset[i]) % m_covariate_max[i];
   
    switch (i) {
      case k_ref_base:
        ref_base = index2basechar(j);
        break;

      case k_obs_base:
        obs_base = index2basechar(j);
        break;

      case k_quality:
        quality = j;
        break;

      case k_read_pos:
        read_pos = j;
        break;
        
      case k_read_set:
        read_set = j;
        break;    
    }
  }
}

/*  cErrorTable::recalibrate()

  JEB: currently a sandbox for trying out log-odds calculations
*/

void cErrorTable::recalibrate() {
  
  assert(m_covariate_used[k_obs_base]);
  assert(m_covariate_used[k_ref_base]);
  
  covariates_used_t covariates;
  for (int i=0; i<k_num_covariates; i++) {
    covariates[i]=0;
  }
  covariates[k_ref_base] = true;
  covariates[k_obs_base] = true;
  cErrorTable total_table(*this, covariates);
  total_table.print(std::string("total.tab"));

  // Need to calculate sub-tables for each covariate
  // One row for each value of covariates, One column for each observed base
  for (int i=0; i<k_num_covariates; i++) {
    covariates[i]=0;
  }
  covariates[k_ref_base] = true;
  covariates[k_obs_base] = true;
  covariates[k_quality] = true;
  
  cErrorTable subtable(*this, covariates);
  subtable.print(std::string("newoutput.tab"));
  
  subtable.calculate_odds_ratios(total_table);
}

void cErrorTable::calculate_odds_ratios(cErrorTable& total_table) {

  char ref_base;
  char obs_base; 
  uint32_t quality; 
  uint32_t read_pos; 
  uint32_t read_set;

  for (uint32_t i=0; i<m_count_table.size(); i++) {
    index_to_covariates(i, ref_base, obs_base, quality, read_pos, read_set);
    uint32_t j = total_table.covariates_to_index(ref_base, obs_base, quality, read_pos, read_set);
    
    for (uint32_t p=0; p<5; p++) {
      char b1 = index2basechar(p);
      
      // Find the total in this category
      uint32_t total_total = 0;
      uint32_t this_total = 0;
      
      for (uint32_t q=0; q<5; q++) {
        char b2 = index2basechar(q);
        uint32_t t1 = total_table.covariates_to_index(b1, b2, quality, read_pos, read_set);
        uint32_t t2 = covariates_to_index(b1, b2, quality, read_pos, read_set);
        
        total_total += total_table.m_count_table[t1] + 1;
        this_total += m_count_table[t2] + 1;
      }
      
      for (uint32_t q=0; q<5; q++) {
        char b2 = index2basechar(q);

        uint32_t t1 = total_table.covariates_to_index(b1, b2, quality, read_pos, read_set);
        uint32_t t2 = covariates_to_index(b1, b2, quality, read_pos, read_set);
        
        uint32_t total_obs = total_table.m_count_table[t1] + 1;
        uint32_t this_obs = m_count_table[t2] + 1;
        
        double m_log_odds_ratio = log((double)this_total / (total_total - this_total)) - log((double)this_obs / (total_obs - this_obs));
        double adjust = exp(m_log_odds_ratio);
      }
    }
  }
}

} //namespace breseq

