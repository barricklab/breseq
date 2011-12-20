/*****************************************************************************

AUTHORS

  Jeffrey E. Barrick <jeffrey.e.barrick@gmail.com>
  David B. Knoester

LICENSE AND COPYRIGHT

  Copyright (c) 2008-2010 Michigan State University
  Copyright (c) 2011 The University of Texas at Austin

  breseq is free software; you can redistribute it and/or modify it under the  
  terms the GNU General Public License as published by the Free Software 
  Foundation; either version 1, or (at your option) any later version.

*****************************************************************************/

#include "libbreseq/common.h"
#include "libbreseq/pileup.h"
#include "libbreseq/error_count.h"

using namespace std;

namespace breseq {

const string cErrorTable::covariate_names[] = {
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
  @abstract Calls a pileup function to count errors along the reference bases.
  @param bam File name for bam file
  @param fasta File name for fasta file
  @param readfiles Vector of file names for read files
  @param do_coverage Output .unique_only_coverage_distribution.tab
  @param do_errors
  @param min_qual_score
  @param covariates
 */
void error_count(
                 Summary& summary,
                 const string& bam,
                 const string& fasta,
								 const string& output_dir,
								 const vector<string>& readfiles,
                 bool do_coverage,
                 bool do_errors,
                 uint8_t min_qual_score,
                 const string& covariates
                 ) 
{
	error_count_pileup ecp(summary, bam, fasta, output_dir, do_coverage, do_errors, min_qual_score, covariates);
	ecp.do_pileup();
	if (do_coverage) ecp.print_coverage();
	if (do_errors) ecp.print_error(readfiles);
}


/*! Constructor.
 */
error_count_pileup::error_count_pileup(
                                       Summary& _summary,
                                       const string& bam, 
                                       const string& fasta, 
                                       const string& output_dir,
                                       bool do_coverage, 
                                       bool do_errors, 
                                       uint8_t min_qual_score, 
                                       const string& covariates
                                       )
: pileup_base(bam, fasta)
, summary(_summary)
, m_output_dir(output_dir)
, m_do_coverage(do_coverage)
, m_do_errors(do_errors)
, m_min_qual_score(min_qual_score)
, m_error_table(covariates)
{
	// reserve enough space for the sequence info:
	_seq_info.resize(num_targets());
  set_print_progress(true);
  
  // Set up to print out file during pileup
  if (m_error_table.m_per_position)
  {
    string output_file;

    output_file = output_dir + "error_counts.tab";     
    m_per_position_file.open(output_file.c_str());
    assert(!m_per_position_file.fail());
    m_error_table.write_count_table_header(m_per_position_file);
  }
  
  m_read_found_starting_at_pos[0] = 0;
  m_read_found_starting_at_pos[1] = 0;

}


/*! Destructor.
 */
error_count_pileup::~error_count_pileup() {
}


/*! Called for each alignment.
 */
void error_count_pileup::pileup_callback(const pileup& p) {
  
	assert(p.target() < _seq_info.size());
	sequence_info& info=_seq_info[p.target()];

  //cerr << p.target() << " " << p.position_1() << endl;
  
	size_t unique_coverage=0; // number of non-deletion, non-redundant alignments at this position.
	bool has_redundant_reads=false; // flag indicating whether this position has any redundant reads.
  int32_t has_query_start[2] = {0,0};
  
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
        
    if (i->query_position_1() == i->query_stranded_bounds_1().first)
    {
      has_query_start[i->reversed() ? 1 : 0] = 1;
    }
    
		// if we are only tracking coverage, go to the next alignment
    if (!m_do_errors) {
      continue;
    }		

    // count the error
    m_error_table.count_alignment_position(*i, p);
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

    m_read_found_starting_at_pos[has_query_start[0]]++;
    m_read_found_starting_at_pos[has_query_start[1]]++;
	}
  
  // per-position prints each line separately
  if (m_error_table.m_per_position) 
  {
    m_error_table.write_count_table_content(m_per_position_file, p.position_1());
    m_error_table.clear();  
  }
}
  
void error_count_pileup::at_target_end(const uint32_t tid)
{
  double total = static_cast<double>(m_read_found_starting_at_pos[0] + m_read_found_starting_at_pos[1]);
  summary.preprocess_error_count[target_name(tid)].no_pos_hash_per_position_pr = 1.0;
  if (total != 0) {
    summary.preprocess_error_count[target_name(tid)].no_pos_hash_per_position_pr = static_cast<double>(m_read_found_starting_at_pos[0]) / total;
  }

  m_read_found_starting_at_pos[0] = 0;
  m_read_found_starting_at_pos[1] = 0;
}


/*! Print coverage distribution.
  @abstract Creates coverage count table:
       This is a table of non-deletion reads per position to non-redundancy counts.
       For example, given unique_only_coverage[i] = x, for all aligned positions p:
       i is the number of reads that do not indicate a deletion at p
       x is the number of positions that have no redundancies

 */
void error_count_pileup::print_coverage() {
	using namespace std;
	for(size_t i=0; i<_seq_info.size(); ++i) {
		string filename(m_output_dir + "/" + m_bam->header->target_name[i] + ".unique_only_coverage_distribution.tab");
		ofstream out(filename.c_str());					
		
		out << "coverage\tn" << endl;
		for(size_t j=1; j<_seq_info[i].unique_only_coverage.size(); ++j) {
			out << j << "\t" << _seq_info[i].unique_only_coverage[j] << endl;
		}	
		out.close();
	}
}


/*! Print error file.
 */
void error_count_pileup::print_error(const vector<string>& readfiles) {

    string output_file;
  
    // It will be printed during pileup if !m_per_position
    // and we don't want the additional processing.
    if (!m_error_table.m_per_position) 
    {
      output_file = m_output_dir + "/error_counts.tab"; 
      m_error_table.write_count_table(output_file);
      
      output_file = m_output_dir + "/error_rates.tab"; 
      m_error_table.counts_to_log10_prob();
      m_error_table.write_log10_prob_table(output_file);
  
      output_file = m_output_dir + "/base_qual_error_prob.#.tab"; 
      m_error_table.write_base_qual_only_prob_table(output_file, readfiles);
      return;
    }
}


/*  cErrorTable::cErrorTable()

    This version of the initializer is for constructing the full
    table that will be used from command line arguments.
*/
cErrorTable::cErrorTable (const string& colnames)
{
  m_per_position = 0;
  if (colnames.length() == 0) return;
  read_covariates(colnames);
  allocate_table();
}

/*  cErrorTable::cErrorTable()

    Allocates an empty table with the given covariates.
*/
cErrorTable::cErrorTable(covariates_used_t covariate_used, covariates_max_t covariate_max, covariates_enforce_max_t covariate_enforce_max, bool per_position) {

  m_per_position = per_position;

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

  m_per_position = error_table.m_per_position;
  
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
          ASSERT(false, "Covariate \'" + covariate_names[i] + "\'exceeded enforced maximum value of \'" + to_string(m_covariate_max[i]) + "\'.");
        } else {
          cout << "Adjusted down " << val << " to " << m_covariate_max[i]-1 << endl;
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
void cErrorTable::read_covariates(
  const string& colnames, 
  covariates_used_t&         _covariate_used,         // list of covariates that are used by table
  covariates_max_t&          _covariate_max,          // maximum value of each covariate
  covariates_enforce_max_t&  _covariate_enforce_max,  // do not throw an error if max exceeded, reassign value to max
  covariates_offset_t&       _covariate_offset,       // number to multiply this covariate by when constructing row numbers
  bool&                      _per_position
  )
{
  // set default values
  for (int i=0; i<k_num_covariates; i++) {
    _covariate_used[i] = false;
    _covariate_enforce_max[i] = false;
  }

  vector<string> columns_to_use = split(colnames, ",");

  // Turn on and assign max
  for(vector<string>::iterator i=columns_to_use.begin(); i!=columns_to_use.end(); ++i) {
    vector<string> columns_parts = split(*i, "=");

		if (columns_parts[0] == "ref_base") {
      _covariate_used[k_ref_base] = true;
      _covariate_max[k_ref_base] = 5;
    }
		else if (columns_parts[0] == "prev_base") {
      _covariate_used[k_prev_base] = true;
      _covariate_max[k_prev_base] = 5;
    }
		else if (columns_parts[0] == "obs_base") {
      _covariate_used[k_obs_base] = true;
      _covariate_max[k_obs_base] = 5;
    }
    else if (columns_parts[0] == "quality") {
      _covariate_used[k_quality] = true;
      // This is a bit lazy, we could assign a minimum quality offset to save empty bins...
      _covariate_max[k_quality] = atoi(columns_parts[1].c_str());
    }
    else if (columns_parts[0] == "read_set") {
      _covariate_used[k_read_set] = true;
      _covariate_max[k_read_set] = atoi(columns_parts[1].c_str());
    }
    // Note: ref_pos is treated special, because the table would be too big
    // to store, and we can go through the pile_up and print as we can.
    else if (columns_parts[0] == "ref_pos") {
      _per_position = true;
    }
    else if (columns_parts[0] == "read_pos") {
      _covariate_used[k_read_pos] = true;
      _covariate_max[k_read_pos] = atoi(columns_parts[1].c_str());
    }
    else if (columns_parts[0] == "base_repeat") {
      _covariate_used[k_base_repeat] = true;
      _covariate_max[k_base_repeat] = atoi(columns_parts[1].c_str());
      // Base repeat will typically have a cap, above which we combine all values
      _covariate_enforce_max[k_base_repeat] = true;
    }
    else {
      cerr << "Unrecognized covariate: " << columns_parts[0] << endl;
      assert(1);
    }
	}
  
  // Then assign offsets, so they are in a consistent order...
  int current_offset = 1;

  for (int i=0; i<k_num_covariates; i++) {
    if (_covariate_used[i]) {
      _covariate_offset[i] = current_offset;
      current_offset *= _covariate_max[i];
    }
  }
}


/*  cErrorTable::print_covariates()

    Convert covariates to string representation
*/
string cErrorTable::print_covariates() {

  stringstream covariate_string;
  if (m_per_position) covariate_string << "ref_pos";
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



/*  cErrorTable::read_log10_prob_table()

    Read table of log10 probabilities of observations
*/
void cErrorTable::read_log10_prob_table(const string& filename) {

  ifstream in(filename.c_str());
  
  string s;
  
  // First line contains the covariates  
  getline(in, s);
  read_covariates(s);
  allocate_table();
  m_log10_prob_table.resize(m_count_table.size());
  
  // Ignore the header line, we expect order to be constant
  getline(in, s);
  
  // We don't need the beginning of the line, since we know the
  // order everything was printed in. That's just there to make it
  // human readable
  
  // Read one line for each
  for (uint32_t i=0; i< m_log10_prob_table.size(); i++) {
    getline(in, s);
    vector<string> split_line = split(s, "\t");
    m_log10_prob_table[i] = strtod(split_line.back().c_str(), NULL);
  }           
}

/*  cErrorTable::write_log10_prob_table()

    Print out a table of covariates and counts.
*/
void cErrorTable::write_log10_prob_table(const string& filename) {

  ofstream out(filename.c_str());
  
  // First line contains the covariates
  out << print_covariates() << endl;
  
  for (int i=0; i<k_num_covariates; i++) {
    if(m_covariate_used[i]) out << covariate_names[i] << m_sep;
  }
  out << "log10_probability" << endl;

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
        
    out << m_log10_prob_table[idx] << endl;
  }
}
  
  /*  cErrorTable::write_base_qual_only_prob_table()
   
   Write out an error table where only the base and quality covariates are used.
   This is used as input to R for making a plot of error rates.
   */
  void cErrorTable::write_base_qual_only_prob_table(const std::string& filename, const vector<string>& readfiles) {
        
    // we know this is in order...
    //enum {k_read_set, k_ref_base, k_prev_base, k_obs_base, k_quality, k_read_pos, k_base_repeat, k_num_covariates};
    
    // in order baseindex2char(obs_base) * baseindex2char(ref_base); * quality
    vector<double> output_table;
    output_table.resize(m_covariate_max[k_read_set] * m_covariate_max[k_obs_base] * m_covariate_max[k_ref_base] * m_covariate_max[k_quality]);
    
    
    for (uint32_t idx=0; idx<m_count_table.size(); idx++) {
      
      uint32_t add_index = 0;
      
      for (int i=0; i<k_num_covariates; i++) {
        if (!m_covariate_used[i]) continue;
        
        uint32_t j = (idx / m_covariate_offset[i]) % m_covariate_max[i];
        
        if (i == k_read_set) {
          add_index += m_covariate_max[k_obs_base] * m_covariate_max[k_ref_base] * m_covariate_max[k_quality] * j;
        } else if (i == k_ref_base) {
          add_index +=  j;
        } else if (i == k_obs_base) {
          add_index +=  m_covariate_max[k_ref_base] * j;
        } else if (i == k_quality) {
          add_index +=  m_covariate_max[k_obs_base] * m_covariate_max[k_ref_base] * j;
        }
      }
      
      output_table[add_index] += m_count_table[idx];
    }
    
    // Calculate probabilities within each set of m_covariate_max[k_ref_base]
    double running_total = 0;
    
    for (uint32_t r=0; r<m_covariate_max[k_read_set]; r++) {
      
      for (uint32_t k=0; k<m_covariate_max[k_obs_base] * m_covariate_max[k_ref_base] * m_covariate_max[k_quality]; k++) {
        uint32_t i = r * m_covariate_max[k_obs_base] * m_covariate_max[k_ref_base] * m_covariate_max[k_quality] + k;
        running_total += output_table[i];
        if (i % m_covariate_max[k_ref_base] == m_covariate_max[k_ref_base] - 1) {
          
            for (uint32_t j=i-(m_covariate_max[k_ref_base]-1); j<=i; j++) {
              if (running_total > 0) {
                output_table[j] /= running_total;
              }
              if (output_table[j] == 0) {
               output_table[j] = NAN; 
              }
            }
            running_total = 0;
        }
      }
      
      // create the proper filename
      string this_file_name = filename;
      size_t pos = this_file_name.find("#");
      assert(pos);
      this_file_name.replace(pos, 1, readfiles[r]);
      
      std::ofstream out(this_file_name.c_str());
      
      // First line contains the headers: quality and then all possible base changes
      out << "quality";
      
      for (uint32_t b1=0; b1< m_covariate_max[k_ref_base]; b1++) {
        char base1 = baseindex2char(b1);

        for (uint32_t b2=0; b2< m_covariate_max[k_obs_base]; b2++) {
          char base2 = baseindex2char(b2);
          out << "\t" << base1 << base2;
        }
      }
      out << endl;
       
      for (uint32_t q=0; q< m_covariate_max[k_quality]; q++) {
        out << q;
        for (uint32_t b1=0; b1< m_covariate_max[k_ref_base]; b1++) {        
          for (uint32_t b2=0; b2< m_covariate_max[k_obs_base]; b2++) {
            uint32_t idx = r * m_covariate_max[k_obs_base] * m_covariate_max[k_ref_base] * m_covariate_max[k_quality]
                      + q * m_covariate_max[k_obs_base] * m_covariate_max[k_ref_base] + b1 * m_covariate_max[k_ref_base] + b2;
            out << "\t";
            if (isnan(output_table[idx]))
              out << "NA";
            else 
              out << output_table[idx];
          }
        }
        out << endl;
      }
    }
  }

/*  cErrorTable::write_count_table()

    Print out a table of covariates and counts.
*/
void cErrorTable::write_count_table(const string& filename) {

  ofstream out(filename.c_str());
  
  // First line contains the covariates
  write_count_table_header(out);
  write_count_table_content(out);
}
  
void cErrorTable::write_count_table_header(ofstream& out) {
  
  out << print_covariates() << endl;
  
  if (m_per_position) out << "ref_pos" << m_sep;

  for (int i=0; i<k_num_covariates; i++) {
    if(m_covariate_used[i]) out << covariate_names[i] << m_sep;
  }
  out << "count" << endl;

}

  
/*  cErrorTable::write_count_table_content()
 
 Print out a table of covariates and counts.
 */
void cErrorTable::write_count_table_content(ofstream& out, const uint32_t position) {
  
  for (uint32_t idx=0; idx<m_count_table.size(); idx++) {
    
    // if in per_position mode, skip empty lines
    if ((position != 0) && (m_count_table[idx] == 0))
    {
      continue;
    }
    
    if (m_per_position) out << position << m_sep;
    
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
    
    out << m_count_table[idx] << endl;
  }
}


/*  cErrorTable::count_alignment_position()

    Record all observations in an alignment to a column of the reference genome
    in the count table. This function is called by the pileup_callback.
*/
void cErrorTable::count_alignment_position(const pileup_wrapper& i, const pileup& p) {
    
    uint32_t ref_pos = p.position_0();
    const char* ref_seq = p.reference_sequence();
    
		uint32_t reversed = i.reversed(); // are we on the reverse strand?
		uint8_t* qseq = i.read_bam_sequence(); // query sequence (read)
		int32_t q_pos_0 = i.query_position_0(); // 0-indexed
		int32_t q_start_0 = i.query_start_0(); // 0-indexed
		int32_t q_end_0 = i.query_end_0(); // 0-indexed

		uint8_t* qscore = i.read_base_quality_bam_sequence(); // quality score array
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
  
/*  cErrorTable::log10_prob_to_prob()
 
 Change log10 table to probability table 
 */
void cErrorTable::log10_prob_to_prob() {
  
  ASSERT(m_log10_prob_table.size() > 0, "No values loaded into log10 prob table.");
  m_prob_table.resize(m_log10_prob_table.size(), 0.0);
  for (uint32_t i=0; i< m_log10_prob_table.size(); i++) {
    m_prob_table[i] = pow(10, m_log10_prob_table[i]);    
  }  
  
}

/*  cErrorTable::alignment_position_to_covariates()

  Extracts covariates from the alignment and insert index.
  
  Note: Does not fill in the ref_base covariate!
*/

bool cErrorTable::alignment_position_to_covariates(const pileup_wrapper& a, int32_t insert_count, covariate_values_t& cv) {
  
  // -1 for deletion, otherwise 1-number of bases inserted at this position
  int indel=a.on_base_indel();
  base_bam read_base_bam = a.on_base_bam(insert_count);
        
  //##don't use bases without qualities!!
  if(_base_bam_is_N(read_base_bam)) return false;
  
  //## These are the start and end coordinates of the aligned part of the read
  uint32_t q_start_0,q_end_0;
  a.query_bounds_0(q_start_0, q_end_0); // @dk: 1-indexed!
 
  //## (1) Mis(match) in read relative to reference...
  //##       Quality is of the current base in the read, we have ALREADY checked that it is not an N					
  //##       Don't need to do anything.
  //  if (indel == 0)
          
  uint32_t q_pos_0 = a.query_position_0();
  
  //## (2) Deletion in read relative to reference...
  //##       Quality is of the NEXT base in the read, and check that it is not an N
  if (indel == -1)
  {
    q_pos_0 += 1 - a.reversed(); 
    base_bam check_base_bam = a.read_base_bam_0(q_pos_0);
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
    base_bam check_base_bam = a.read_base_bam_0(q_pos_0);
    if (_base_bam_is_N(check_base_bam)) return false;
  }

  //eventually include in above...
  cv.obs_base() = basebam2index(read_base_bam);
  cv.quality() = a.read_base_quality_0(q_pos_0);
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
  
double cErrorTable::get_prob(covariate_values_t& cv) {
  
  assert(m_prob_table.size() > 0);
  uint32_t i = covariates_to_index(cv);
  
  assert(i < m_prob_table.size());
  return m_prob_table[i];
}


} //namespace breseq

