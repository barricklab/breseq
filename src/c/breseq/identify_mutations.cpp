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

#include <iostream>
#include <cmath>
#include <map>
#include <sstream>
#include <string>
#include <vector>
#include <bam.h>
#include <sam.h>
#include <faidx.h>
#include <assert.h>

#include "breseq/common.h"
#include "breseq/pileup.h"
#include "breseq/identify_mutations.h"
#include "breseq/error_count.h"
#include <boost/math/distributions/chi_squared.hpp>


/*! Convenience wrapper around the identify_mutations_pileup class.
 */
void breseq::identify_mutations(
								const string& bam,
								const string& fasta,
								const string& error_dir,
								const string& gd_file,
								const string& output_dir,
								const vector<string>& readfiles,
								const string& coverage_dir,
								const vector<double>& deletion_propagation_cutoff,
								double mutation_cutoff,
								bool predict_deletions,
								bool predict_polymorphisms,
								uint8_t min_qual_score,
								double polymorphism_cutoff,
								double polymorphism_frequency_cutoff,
								const string& error_table_file,
								bool print_per_position_file
 ) {
                                                                                            
	// do the mutation identification:
	identify_mutations_pileup imp(
								bam,
								fasta,
								error_dir,
								gd_file,
								output_dir,
								readfiles,
								coverage_dir,
								deletion_propagation_cutoff,
								mutation_cutoff,
								predict_deletions,
								predict_polymorphisms,
								min_qual_score,
								polymorphism_cutoff,
								polymorphism_frequency_cutoff,
								error_table_file,
								print_per_position_file
							);
	imp.do_pileup();
}


/*! Constructor.
 */
breseq::identify_mutations_pileup::identify_mutations_pileup(
															const string& bam,
															const string& fasta,
															const string& error_dir,
															const string& gd_file,
															const string& output_dir,
															const vector<string>& readfiles,
															const string& coverage_dir,
															const vector<double>& deletion_propagation_cutoff,
															double mutation_cutoff,
															bool predict_deletions,
															bool predict_polymorphisms,
															uint8_t min_qual_score,
															double polymorphism_cutoff,
															double polymorphism_frequency_cutoff,
															const string& error_table_file,
															bool print_per_position_file
                                                            )
: breseq::pileup_base(bam, fasta)
, _ecr(error_dir, readfiles)
, _gd(gd_file)
, _min_qual_score(min_qual_score)
, _deletion_seed_cutoff(0)
, _deletion_propagation_cutoff(deletion_propagation_cutoff)
, _mutation_cutoff(mutation_cutoff)
, _predict_deletions(predict_deletions)
, _predict_polymorphisms(predict_polymorphisms)
, _polymorphism_cutoff(polymorphism_cutoff)
, _polymorphism_frequency_cutoff(polymorphism_frequency_cutoff)
, _coverage_dir(coverage_dir)
, _output_dir(output_dir)
, _log10_ref_length(0)
, _on_deletion_seq_id(UNDEFINED)
, _this_deletion_reaches_seed_value(false)
, _last_position_coverage_printed(0)
, _print_per_position_file(print_per_position_file) {
	  
  assert(m_bam->header->n_targets == (int32_t)_deletion_propagation_cutoff.size());
    
	// reserve enough space for the sequence info:
	_seq_info.resize(m_bam->header->n_targets);
	
  // tally up the reference lengths from the bam file
	for(int i=0; i<m_bam->header->n_targets; ++i) {
		_log10_ref_length += static_cast<double>(m_bam->header->target_len[i]);
	}
	assert(_log10_ref_length != 0);
	_log10_ref_length = log10(_log10_ref_length);
  
  // are we printing detailed coverage information?
  _print_coverage_data = (coverage_dir != "");
  
  // use new error
  m_use_cErrorTable = false;
  if (error_table_file.length() > 0) {
    m_use_cErrorTable = true;
    m_error_table.read_log10_prob_table(error_table_file);
  }
  
  if (_print_per_position_file) {
    string filename(_output_dir);
		filename += "/full_identify_mutation.out";
    _per_position_file.open(filename.c_str());
  }

	_on_deletion_seq_id = UNDEFINED;
	_last_deletion_start_position = UNDEFINED;
	_last_deletion_end_position = UNDEFINED;
	_last_deletion_redundant_start_position = UNDEFINED;
	_last_deletion_redundant_end_position = UNDEFINED;
	_last_start_unknown_interval = UNDEFINED;
}


/*! Destructor.
 */
breseq::identify_mutations_pileup::~identify_mutations_pileup()
{
}


/*! Called for each alignment.
 */
void breseq::identify_mutations_pileup::pileup_callback(const breseq::pileup& p) {
	using namespace std;
	assert(p.target() < _seq_info.size());
  _this_deletion_propagation_cutoff = _deletion_propagation_cutoff[p.target()];
  
  uint32_t position = p.position_1();
  //cout << position << endl;
  
	int insert_count=-1;
	bool next_insert_count_exists=true;
	
	// check to see if we already opened the coverage file:
	if(_print_coverage_data && !_coverage_data.is_open()) {
		//		open COV, ">$coverage_tab_file_name" if (defined $coverage_tab_file_name);
		//		print COV join("\t", 'unique_top_cov', 'unique_bot_cov', 'redundant_top_cov', 'redundant_bot_cov', 'raw_redundant_top_cov', 'raw_redundant_bot_cov', 'e_value', 'position') . "\n";
		string filename(_coverage_dir);
		filename += p.target_name();
		filename += ".coverage.tab";
		_coverage_data.open(filename.c_str());
		_coverage_data << "unique_top_cov" << "\t" << "unique_bot_cov" << "\t" << "redundant_top_cov" << "\t" << "redundant_bot_cov" << "\t" << "raw_redundant_top_cov" << "\t" << "raw_redundant_bot_cov" << "\t" << "e_value" << "\t" << "position" << endl;
	}	
  
  // temporary file for debugging polymorphism prediction
  if(_predict_polymorphisms && !_polymorphism_r_input_file.is_open()) {
		string filename(_output_dir);
		filename += "/polymorphism_statistics_input.tab";
		_polymorphism_r_input_file.open(filename.c_str());
		_polymorphism_r_input_file 
      << "seq_id" << "\t"
      << "position" << "\t" 
      << "insert_position" << "\t" 
      << "ref_base" << "\t" 
      << "best_base" << "\t" 
      << "second_best_base" << "\t" 
      << "frequency" << "\t"
      << "log10_base_likelihood" << "\t"
      << "E-value" << "\t"
      << "ref_top_strand" << "\t" 
      << "ref_bot_strand" << "\t" 
      << "new_top_strand" << "\t" 
      << "new_bot_strand" << "\t"
      << "best_quals" << "\t"
      << "second_best_quals" << "\t"
      << endl;
  }
	
	while(next_insert_count_exists) {
		++insert_count; // we're doing this here because the perl version uses a while-continue.
    next_insert_count_exists = false;
		
		base_char ref_base_char = '.';
		if(!insert_count) {
			ref_base_char = p.reference_base_char_1(position);
		}
		
		//## zero out the info about this position
		map<uint8_t,position_info> pos_info;
		for(size_t j=0; j<base_list_size; ++j) {
			pos_info.insert(make_pair(base_char_list[j],position_info()));
		}
		
		//## keep track of coverage for deletion prediction
		position_coverage this_position_coverage;
		bool this_position_unique_only_coverage=true;
		
    //## SNP caller
    cDiscreteSNPCaller snp(1);
        
		//## polymorphism prediction data
		vector<polymorphism_data> pdata;
    
		//## for each alignment within this pileup:
		for(pileup::const_iterator i=p.begin(); i!=p.end(); ++i) {

      //## After these substitutions...
      //## Indel is -1 if the ref base is deleted in the read,
      //## Zero if the read base is aligned to a ref base, and
      //## Positive if the read base is an insertion relative to the ref base
      int indel=i->indel();
      if(indel < 0) {
        indel = 0;
      }
      if(i->is_del()) {
        indel = -1;
      }
      
      base_bam read_base_bam='.';
      if(indel >= insert_count) {
        read_base_bam = i->query_base_bam_0(i->query_position_0() + insert_count);
      }
            
      //##don't use bases without qualities!!
      if(_base_bam_is_N(read_base_bam)) continue;
      
      //## gather information about the aligned base
      int32_t redundancy = i->redundancy();
      int32_t fastq_file_index = i->fastq_file_index();
      int strand = i->strand();
      bool trimmed = i->is_trimmed();
            
      //### Optionally, only count reads that completely match
      //			my $complete_match = 1;
      //			if ($settings->{require_complete_match})
      //			{
      //				$complete_match = ($q_start == 1) && ($q_end == $a->l_qseq);
      //				next if (!$complete_match);
      //			}
      //### End complete match condition
      
      //##### update coverage if this is not a deletion in read relative to reference
      //### note that we count trimmed reads here, but not when looking for short indel mutations...	
      
      if(redundancy == 1) {
        //## keep track of unique coverage	
        ++this_position_coverage.unique[1+strand];
        ++pos_info[basebam2char(read_base_bam)].unique_cov[1+strand];
        
        //## we don't continue to consider further insertions relative
        //## to the reference unless uniquely aligned reads have them 
        if(indel > insert_count) {
          next_insert_count_exists = true;
        }
        
      } else {		
        //## mark that this position has some non-unique coverage
        this_position_unique_only_coverage = false;
        //## keep track of redundant coverage
        this_position_coverage.redundant[1+strand] += 1.0/redundancy;
        ++this_position_coverage.raw_redundant[1+strand];
      }
      
      //## don't use information from trimmed or redundant reads
      //## to predict base substitutions and short indels!!
      if(trimmed || (redundancy > 1)) {
        continue;
      }

 			
			//##### deal with base calls
      //cerr << "POSITION:" << position << endl;
      
      if (m_use_cErrorTable) {
      
        covariate_values_t cv; 

        bool is_ok = m_error_table.alignment_position_to_covariates(*i, insert_count, cv);
        //cv.obs_base is still not a char here...
        
        if (is_ok)  {
        
          if (cv.quality() < _min_qual_score) {
            continue;
          }
  
          //## this is the coverage for SNP counts, tabulate AFTER skipping trimmed reads
          ++pos_info[baseindex2char(cv.obs_base())].unique_trimmed_cov[1+strand];
          
          //##### this is for polymorphism prediction and making strings
          pdata.push_back(polymorphism_data(baseindex2char(cv.obs_base()),cv.quality(),i->strand(),cv.read_set(), cv));
          
          //cerr << " " << cv.obs_base() << " " << (char)ref_base << endl;

          snp.update(cv, strand == 1, m_error_table);
        }
      }
      else {        
        
        //## These are the start and end coordinates of the aligned part of the read
        //			my ($q_start, $q_end) = Breseq::Shared::alignment_query_start_end($a, {no_reverse=>1});
        int32_t q_start,q_end;
        i->query_bounds_1(q_start, q_end); // @dk: 1-indexed!
        
        uint8_t quality=0;
        
        //## Deletion in read relative to reference...
        //## Quality is of the NEXT base in the read, and check that it is not an N
        //## Note: This is for a deletion when $insert_count == 0 and against an insertion when $insert_count > 0

        if (indel == -1)
        {			
          int32_t mqpos = i->query_position_0() + 1 - i->reversed(); 
          base_bam check_base_bam = i->query_base_bam_0(mqpos);
          if (_base_bam_is_N(check_base_bam)) continue;
          quality = i->quality_base_0(mqpos);
        }

        //## Substitution in read relative to reference...
        //## Quality is of the current base in the read, we have ALREADY checked that it is not an N					
        else if (insert_count == 0)
        {
          quality = i->quality_base_0(i->query_position_0());
        }
        
        //## Insertion in read relative to reference...
        //## Quality is of the NEXT base in the read, and check that it is not an N
        //## Note that it is possible this read base may be a '.' (supporting the non-insert call)
        else //if (insert_count > 0) 
        {		
          int32_t max_offset = insert_count;
          if (indel < max_offset) max_offset = indel;
          int32_t mqpos = i->query_position_0() + max_offset + 1 - i->reversed(); 
                  
          //## Check bounds: it's possible to go past the end of the read because
          //## this is the last base of this read, but other reads have inserted bases
        
          if (mqpos >= q_end) continue;
          //next ALIGNMENT if ($mqpos > $q_end);
        
          base_bam check_base_bam = i->query_base_bam_0(mqpos);
          if (_base_bam_is_N(check_base_bam)) continue;
          
          quality = i->quality_base_0(mqpos);
        }

        //## We may want to ignore all bases below a certain quality when calling mutations and polymorphisms
        //## This is the check for whether the base fails; it should be after coverage counting
       if (quality < _min_qual_score) {
  //        cerr << position << " " << (unsigned int)quality << " " << endl;
          continue;
        }
        
        
        //## this is the coverage for SNP counts, tabulate AFTER skipping trimmed reads
        ++pos_info[basebam2char(read_base_bam)].unique_trimmed_cov[1+strand];
        
        //##### this is for polymorphism prediction and making strings
        pdata.push_back(polymorphism_data(basebam2char(read_base_bam),quality,strand,fastq_file_index));

        snp.update(read_base_bam, quality, strand == 1, fastq_file_index, _ecr);
      }
		} // end for-each read
		
    
    //#############################
		//## PER POSITION/INSERT COUNT
    //#############################
    
		//#sum up coverage observations	
		this_position_coverage.sum();
				
		//#we are trying to find the base with the most support

    pair<uint8_t,double> snp_pred = snp.get_prediction();
    base_char best_base_char = snp_pred.first;
    double e_value_call = snp_pred.second - _log10_ref_length;
    
    //Do we predict a base at this position?
    bool base_predicted=false;
    if(e_value_call >= _mutation_cutoff) {
      base_predicted = true;
    }
    
		//cerr << position << " e:" << e_value_call << " b:" << base_predicted << endl;

		int total_cov[3]={0,0,0}; // triple, same as above
    
    // Don't need to print, but is nice for debug
		ostringstream line;
		if (_print_per_position_file) {
      line << position << " " << insert_count << " " << ref_base_char << " " << e_value_call;
		}
    
		for(size_t j=0; j<base_list_size; ++j) {
			double top_cov = pos_info[base_char_list[j]].unique_trimmed_cov[2];
			double bot_cov = pos_info[base_char_list[j]].unique_trimmed_cov[0];
			total_cov[2] += (int)round(top_cov);
			total_cov[0] += (int)round(bot_cov);
      if (_print_per_position_file) {
        line << " " << base_char_list[j] << " (" << bot_cov << "/" << top_cov << ")";
      }
    }
    
    // Debug: print additional information to file.
    if (_print_per_position_file) {
      _per_position_file << line.str() << endl;
		}
		
		//###
		//## DELETION DELETION DELETION
		//###
		
    
		//#update information on deletions
		if(insert_count == 0) {
			if(_predict_deletions) {
        // @JEB: note change in call so position sent to check_deletion_completion is 1-based
				check_deletion_completion(position, p.target(), this_position_coverage, e_value_call);
				// @dk: skip update_copy_number_variation(pos, this_position_coverage, ref_base);
			}
		}
		
		//###
		//## POLYMORPHISM POLYMORPHISM POLYMORPHISM
		//###								

		bool polymorphism_predicted=false;
    polymorphism_prediction ppred;
    base_char second_best_base_char;
  
		if(_predict_polymorphisms) {
	
        // Debug output
        /* 
        cerr << position << endl;  
         for(size_t j=0; j<5; ++j) {
            cerr << " " << base_list[j] << " " << (pos_info[base_list[j]].unique_trimmed_cov[0] + pos_info[base_list[j]].unique_trimmed_cov[2]);
        }
        cerr << endl;
        */
        
        // Find the bases with the highest and second highest coverage
        // We only predict polymorphisms involving these
        base_index best_base_index;
        int best_base_coverage = 0;
        base_index second_best_base_index;
        int second_best_base_coverage = 0;
        
        vector<double> snp_probs = snp.get_log10_probabilities();
                
        for (uint8_t i=0; i<base_list_size; i++) {
          base_char this_base_char = base_char_list[i];
          base_index this_base_index = i;
          int this_base_coverage = pos_info[this_base_char].unique_trimmed_cov[0] + pos_info[this_base_char].unique_trimmed_cov[2];
          
          // if better coverage or tied in coverage and better probability
          if ((this_base_coverage > best_base_coverage) ||
            ((this_base_coverage == best_base_coverage) && (snp_probs[this_base_index] > snp_probs[best_base_index]))) {
            second_best_base_index = best_base_index;
            second_best_base_coverage = best_base_coverage;
            best_base_index = this_base_index;
            best_base_coverage = this_base_coverage;
          }
          else if ((this_base_coverage > second_best_base_coverage) 
            || ((this_base_coverage == second_best_base_coverage) && (snp_probs[this_base_index] > snp_probs[second_best_base_index]))) {
            second_best_base_index = this_base_index;
            second_best_base_coverage = this_base_coverage;
          }
        }
        
        // Only try mixed SNP model if there is coverage for more than one base!
        if (second_best_base_coverage) {
          best_base_char = base_char_list[best_base_index];
          second_best_base_char = base_char_list[second_best_base_index];

          ppred = predict_polymorphism(best_base_char, second_best_base_char, pdata);
          ppred.log10_e_value = -(log(ppred.p_value)/log(10)) - _log10_ref_length;
   
          if (ppred.log10_e_value >= 0.0) {
            polymorphism_predicted = 1;
          }
          //cerr << ppred.frequency << " " << ppred.log10_base_likelihood << " " << ppred.p_value << endl;
        }		
      
		}				
		
		
		
		//###
		//## UNKNOWN UNKNOWN UNKNOWN
		//###				
		//		if ($insert_count == 0)
		if(insert_count == 0) {
			update_unknown_intervals(position, p.target(), base_predicted, this_position_unique_only_coverage);
		}
		
		//## evaluate whether to call an actual mutation!				
		//### skip if there is not enough evidence for a call or if it agrees with the reference
		if(isnan(e_value_call) || (e_value_call < -_log10_ref_length)) {
			continue;
		}
    
		//cerr << e_value_call << endl;
		
		//## mutation and polymorphism are exclusive predictions.
    
    // ----> potential problem 
		bool mutation_predicted = !polymorphism_predicted && (best_base_char != ref_base_char);
				
		//## bail if it's just the reference base and we aren't interested in polymorphisms...
		if(!mutation_predicted && !polymorphism_predicted) {
			continue; // goes to next insert count...
		}
    		
		//cerr << position << " " << e_value_call << " " << mutation_predicted << " " << polymorphism_predicted << " " << base2char(best_base) << " " << base2char(ref_base) << endl;

    //## Create new base mutation for genome diff
		ra mut(to_string(_gd.new_id()), "");
    
		//## Fields common to consensus mutations and polymorphisms
		mut[SEQ_ID] = p.target_name();
		mut[POSITION] = to_string<uint32_t>(position);
		mut[INSERT_POSITION] = to_string<uint32_t>(insert_count);
		mut[QUALITY] = formatted_double(e_value_call, 1).to_string();
      
    // both should never be true!
    assert( !(mutation_predicted && polymorphism_predicted) );
    
    //## Specific initializations for consensus mutations
		if(mutation_predicted) {
			mut[REF_BASE] = ref_base_char;
			mut[NEW_BASE] = best_base_char;
			mut[FREQUENCY] = "1";
			if(e_value_call < _mutation_cutoff) {
        breseq::add_reject_reason(mut, "EVALUE");
			}
    }
    //## Specific initilizations for polymorphisms
    else if (polymorphism_predicted) {
 
      
/// PROBLEM: We need to more robustly deal with cases of polymorphisms where neither of the two bases 
///          involved is the reference base; Solve by adding two lines to genome diff?
                
			//# the frequency returned is the probability of the FIRST base
			//# we want to quote the probability of the second base (the change from the reference).      
			if (best_base_char == ref_base_char) {
        mut[REF_BASE] = best_base_char;
        mut[NEW_BASE] = second_best_base_char;
        mut[FREQUENCY] = formatted_double(1 - ppred.frequency, kPolymorphismFrequencyPrecision).to_string();
      } else if (second_best_base_char == ref_base_char) {
        mut[REF_BASE] = second_best_base_char;
        mut[NEW_BASE] = best_base_char;
        mut[FREQUENCY] = formatted_double(ppred.frequency, kPolymorphismFrequencyPrecision).to_string();
      } else {
        cerr << "Warning: polymorphism between two bases not including reference base found at position " << position << endl;
        mut[REF_BASE] = best_base_char;
        mut[NEW_BASE] = second_best_base_char;
        mut[FREQUENCY] = formatted_double(1 - ppred.frequency, kPolymorphismFrequencyPrecision).to_string();
        mut[ERROR] = "polymorphic_without_reference_base";
      }
      
			mut[POLYMORPHISM_QUALITY] = formatted_double(ppred.log10_e_value, kMutationQualityPrecision).to_string();
			if (ppred.log10_e_value < _polymorphism_cutoff ) {
        breseq::add_reject_reason(mut, "EVALUE");
      } 

 /// <--- PROBLEM     
      
      // Need to create lists of all quality scores for and against base for R output  
      
      if (ref_base_char != second_best_base_char) {
        _polymorphism_r_input_file 
          << p.target_name() << "\t"
          << position << "\t"
          << insert_count << "\t"
          << ref_base_char << "\t"
          << best_base_char << "\t"
          << second_best_base_char << "\t"
          << ppred.frequency << "\t"
          << ppred.log10_base_likelihood << "\t"
          << (0.1 * round(ppred.log10_e_value*10)) << "\t"
          << pos_info[best_base_char].unique_trimmed_cov[2] << "\t"
          << pos_info[best_base_char].unique_trimmed_cov[0] << "\t"
          << pos_info[second_best_base_char].unique_trimmed_cov[2] << "\t"
          << pos_info[second_best_base_char].unique_trimmed_cov[0] << "\t"
        ;    
      } else {
        _polymorphism_r_input_file 
          << p.target_name() << "\t"
          << position << "\t"
          << insert_count << "\t"
          << ref_base_char << "\t"
          << second_best_base_char << "\t"
          << best_base_char << "\t"
          << (1-ppred.frequency) << "\t"
          << ppred.log10_base_likelihood << "\t"
          << (0.1 * round(ppred.log10_e_value*10)) << "\t"
          << pos_info[second_best_base_char].unique_trimmed_cov[2] << "\t"
          << pos_info[second_best_base_char].unique_trimmed_cov[0] << "\t"
          << pos_info[best_base_char].unique_trimmed_cov[2] << "\t"
          << pos_info[best_base_char].unique_trimmed_cov[0] << "\t"
        ;    
      }

      
      string best_base_qualities;
      string second_best_base_qualities;

      for(vector<polymorphism_data>::iterator it=pdata.begin(); it<pdata.end(); ++it) {
                
        if (it->_base_char == best_base_char) {
          if (best_base_qualities.length() > 0) {
            best_base_qualities += ",";
          }
          stringstream convert_quality;
          convert_quality << (unsigned int)it->_quality;
          best_base_qualities += convert_quality.str();
        }
        
        if (it->_base_char == second_best_base_char) {
          if (second_best_base_qualities.length() > 0) {
            second_best_base_qualities += ",";
          }
          stringstream convert_quality;
          convert_quality << (unsigned int)it->_quality;
          second_best_base_qualities += convert_quality.str();
        }
      }
      
      if (ref_base_char != second_best_base_char) {
        _polymorphism_r_input_file << best_base_qualities << "\t";
        _polymorphism_r_input_file << second_best_base_qualities << "\t";
      } else {
        _polymorphism_r_input_file << second_best_base_qualities << "\t";
        _polymorphism_r_input_file << best_base_qualities << "\t";
      }
      
      _polymorphism_r_input_file << endl;
			
      if ( (ppred.frequency < _polymorphism_frequency_cutoff) || (ppred.frequency > 1 -_polymorphism_frequency_cutoff) ) {
        breseq::add_reject_reason(mut, "FREQ");
      }
		}
    
		//## More fields common to consensus mutations and polymorphisms
		//## ...now that ref_base and new_base are defined
		int* ref_cov = pos_info[from_string<base_char>(mut[REF_BASE])].unique_trimmed_cov;
		mut[REF_COV] = to_string(make_pair(ref_cov[2], ref_cov[0]));
		
		int* new_cov = pos_info[from_string<base_char>(mut[NEW_BASE])].unique_trimmed_cov;
		mut[NEW_COV] = to_string(make_pair(new_cov[2], new_cov[0]));
		
		mut[TOT_COV] = to_string(make_pair(total_cov[2], total_cov[0]));
		
		_gd.add(mut);
	}
}



/*! Called at the end of the pileup.
 */
void breseq::identify_mutations_pileup::at_target_end(const uint32_t tid) {

  // end "open" intervals
	check_deletion_completion(target_length(tid)+1, tid, position_coverage(numeric_limits<double>::quiet_NaN()), numeric_limits<double>::quiet_NaN());
  update_unknown_intervals(target_length(tid)+1, tid, true, false);

  // write genome diff file
	_gd.write();
	
  // close open files
	_coverage_data.close();
  _polymorphism_r_input_file.close();
}


/*! Helper method to track information about putative deleted regions.
 
 Used at each pileup iteration and at the end.
 //## when called at the end of a fragment, the position is fragment_length+1
 //## and $this_position_coverage is undefined
 
 @JEB This function expects 1-indexed positions!!!
 
 */
void breseq::identify_mutations_pileup::check_deletion_completion(uint32_t position, uint32_t seq_id, const position_coverage& this_position_coverage, double e_value_call) {

	//cerr << position << " " << e_value_call << endl;
	
  // special case = beginning of new seq_id
  if (position == 1) _last_position_coverage = position_coverage(numeric_limits<double>::quiet_NaN());
  
	//## called with an undef $this_position_coverage at the end of the genome
  // print to optional output file
  if (!isnan(this_position_coverage.unique[1]) && _coverage_data.is_open()) {
    _coverage_data << this_position_coverage.unique[0] << "\t"
    << this_position_coverage.unique[2] << "\t"
    << this_position_coverage.redundant[0] << "\t"
    << this_position_coverage.redundant[2] << "\t"
    << this_position_coverage.raw_redundant[0] << "\t"
    << this_position_coverage.raw_redundant[2] << "\t"
    << e_value_call << "\t" << position << endl;
  }
  
  
  //## UNIQUE COVERAGE
  //#start a new possible deletion if we fall below the propagation cutoff
  if(this_position_coverage.unique[1] <= _this_deletion_propagation_cutoff) {
    if(_last_deletion_start_position == UNDEFINED) {
      _last_deletion_start_position = position;
      _left_outside_coverage_item = _last_position_coverage;
      _left_inside_coverage_item = this_position_coverage;
    }
  }
		
  //##keep track of whether we've encountered the seed value
  //		if ($this_position_coverage->{total} <= $deletion_seed_cutoff)
  if(!isnan(this_position_coverage.unique[1]) && (this_position_coverage.total <= _deletion_seed_cutoff)) {
    _this_deletion_reaches_seed_value = true;
  }
    
  //## REDUNDANT COVERAGE
  //## updated only if we are currently within a deletion
  if (_last_deletion_start_position != UNDEFINED) {
    
    if (this_position_coverage.redundant[1] == 0) {
      _this_deletion_redundant_reached_zero = true;
      _last_deletion_redundant_end_position = UNDEFINED;
    }
    else if (this_position_coverage.redundant[1] > 0) {
    //## if there is any redundant coverage remember the start (until we find zero redundant coverage)
      if (!_this_deletion_redundant_reached_zero) {
        _last_deletion_redundant_start_position = position;
      }
      else {
        if (_last_deletion_redundant_end_position == UNDEFINED) _last_deletion_redundant_end_position = position;
      }
    }
  }
	
	//## If we are in a deletion and rise back above the propagation cutoff OR we are at the end of this fragment (NAN),
  //## then record the current deletion.
	if( (_last_deletion_start_position != UNDEFINED) && 
     ( isnan(this_position_coverage.unique[1]) || (this_position_coverage.unique[1] > _this_deletion_propagation_cutoff) ) )
  {
		
		if(_this_deletion_reaches_seed_value) {

      _last_deletion_end_position = position-1;
      if (_last_deletion_redundant_end_position == UNDEFINED) _last_deletion_redundant_end_position = _last_deletion_end_position;
      if (_last_deletion_redundant_start_position == UNDEFINED) _last_deletion_redundant_start_position = _last_deletion_start_position;

			mc del(to_string(_gd.new_id()), "");
			del[SEQ_ID] = target_name(seq_id);
			del[START] = to_string<uint32_t>(_last_deletion_start_position);
			del[END] = to_string<uint32_t>(_last_deletion_end_position);
			del[START_RANGE] = to_string<uint32_t>(_last_deletion_redundant_start_position - _last_deletion_start_position);
			del[END_RANGE] = to_string<uint32_t>(_last_deletion_end_position - _last_deletion_redundant_end_position);
			
      del[LEFT_OUTSIDE_COV] = formatted_double(_left_outside_coverage_item.unique[1], 0).to_string();
      del[LEFT_INSIDE_COV] = formatted_double(_left_inside_coverage_item.unique[1], 0).to_string();
      del[RIGHT_INSIDE_COV] = formatted_double(_last_position_coverage.unique[1], 0).to_string();
      del[RIGHT_OUTSIDE_COV] = formatted_double(this_position_coverage.unique[1], 0).to_string();
      
			_gd.add(del);
		}
		
		//#reset the search
		_this_deletion_reaches_seed_value = false;
		_this_deletion_redundant_reached_zero = false;
		_last_deletion_start_position = UNDEFINED;
		_last_deletion_end_position = UNDEFINED;
		_last_deletion_redundant_start_position = UNDEFINED;
		_last_deletion_redundant_end_position = UNDEFINED;
	}
	
  
  _last_position_coverage = this_position_coverage;
}


/*! Helper method to track unknowns.
 */
void breseq::identify_mutations_pileup::update_unknown_intervals(uint32_t position, uint32_t seq_id, bool base_predicted, bool this_position_unique_only_coverage) 
{
  //debug
  /*
	cerr << position << " " << base_predicted << " " << this_position_unique_only_coverage << endl;
	if(_last_start_unknown_interval == DEFINED) {
		cerr << *_last_start_unknown_interval << endl;
	} else {
		cerr << "undef" << endl;
	}
  */
	
	if(!base_predicted) {
		if(this_position_unique_only_coverage) {
			++s.coverage_unique_uncalled;
		}
		if(_last_start_unknown_interval == UNDEFINED) {
			_last_start_unknown_interval = position;
		}
	}	else {
		if(this_position_unique_only_coverage) {
			++s.coverage_unique_called;
		}
			
		//#end interval where we were unable to call mutations
		if(_last_start_unknown_interval != UNDEFINED) {
			un new_interval(to_string(_gd.new_id()), "");
			new_interval[SEQ_ID] = target_name(seq_id);
			new_interval[START] = to_string<uint32_t>(_last_start_unknown_interval);
			new_interval[END] = to_string<uint32_t>(position - 1);
			_gd.add(new_interval);
			
			_last_start_unknown_interval = UNDEFINED;
		}
	}
}

/*! Predict the significance of putative polymorphisms.
 */
breseq::polymorphism_prediction breseq::identify_mutations_pileup::predict_polymorphism (base_char best_base_char, base_char second_best_base_char, vector<polymorphism_data>& pdata ) {
    
  //#calculate the likelihood of observed reads given this position is 100% the best base  
	double log10_likelihood_of_one_base_model = 0;
  for(vector<polymorphism_data>::iterator it=pdata.begin(); it<pdata.end(); ++it) {
  
    double log10_correct_pr;
    if (m_use_cErrorTable) {
      covariate_values_t this_cv = it->_cv;
      
      if(it->_strand == 1) {
        this_cv.ref_base() = basechar2index(best_base_char);
      } else {
        this_cv.ref_base() = basechar2index(complement_base_char(best_base_char)); 
        this_cv.obs_base() = complement_base_index(this_cv.obs_base());
      }
      log10_correct_pr = m_error_table.get_log10_prob(this_cv);
    } else {
      string base_key;
      if(it->_strand == 1) {
        base_key += best_base_char;
        base_key += it->_base_char;
      } else {
        base_key += complement_base_char(best_base_char); 
        base_key += complement_base_char(it->_base_char);
      }
      log10_correct_pr = _ecr.log10_correct_rates(it->_fastq_file_index, it->_quality, base_key);
    }
    
    log10_likelihood_of_one_base_model  += log10_correct_pr;
  }

	vector<uint8_t> best_base_qualities;
	vector<uint8_t> second_best_base_qualities;
	uint32_t best_base_strand_hash[] = {0, 0};
	uint32_t second_best_base_strand_hash[] = {0, 0};
  
  for(vector<polymorphism_data>::iterator it=pdata.begin(); it<pdata.end(); ++it) {
  
  		if (it->_base_char == best_base_char) {
        best_base_qualities.push_back(it->_quality);
        best_base_strand_hash[it->_strand]++;
      }
      else if (it->_base_char == second_best_base_char) {
        second_best_base_qualities.push_back(it->_quality);
        second_best_base_strand_hash[it->_strand]++;
      }
  }
  
  
	//## Maximum likelihood of observing alignment if sequenced bases were a mixture of the top two bases  
  pair<double,double> best_two_base_model = best_two_base_model_log10_likelihood(best_base_char, second_best_base_char, pdata);
  double max_likelihood_fr_first_base = best_two_base_model.first;
  double log10_likelihood_of_two_base_model = best_two_base_model.second;
    
  //## Likelihood ratio test
  double log10_likelihood_difference = log10_likelihood_of_one_base_model - log10_likelihood_of_two_base_model;

  //debug output 
  /*
  cerr  << "ML Best Base Fraction: " << max_likelihood_fr_first_base << endl;
  cerr  << " Log10 Likelihood (one base model): " << log10_likelihood_of_one_base_model << endl;
  cerr  << " Log10 Likelihood (two base model): " << log10_likelihood_of_two_base_model << endl;
  cerr  << " Log10 Likelihood (different): " << log10_likelihood_difference << endl;
  */
    
  long double p_value = 1;
  if (max_likelihood_fr_first_base != 1.0) {
    double likelihood_ratio_test_value = -2*log(10)*log10_likelihood_difference;
    
    //boost::math::chi_squared myChiSquared(1.0L);
    //p_value = boost::math::pdf(myChiSquared, likelihood_ratio_test_value);
    //cerr << "likelihood_ratio_test_value: " << likelihood_ratio_test_value << " p-value: " << p_value << endl;

    p_value = pchisq(1.0L, likelihood_ratio_test_value);
    //cerr << "likelihood_ratio_test_value: " << likelihood_ratio_test_value << " p-value: " << p_value << endl;
  }

  //debug output 
  /*
  cerr
    << " Log10 Likelihood Difference (one vs two base model): " << (log10_likelihood_of_one_base_model - log10_likelihood_of_two_base_model) 
    << " P-value: " << p_value 
    << endl;
  */
   
  polymorphism_prediction p(max_likelihood_fr_first_base, log10_likelihood_of_one_base_model - log10_likelihood_of_two_base_model, p_value);
		
	return p;
}


/*! Find the best fraction for the best base at a polymorphic site.
 */
pair<double,double> breseq::identify_mutations_pileup::best_two_base_model_log10_likelihood(base_char best_base_char, base_char second_best_base_char, vector<polymorphism_data>& pdata)
{	
	double cur_pr_first_base = 1;
	double cur_log_pr = calculate_two_base_model_log10_likelihood(best_base_char, second_best_base_char, pdata, cur_pr_first_base);

	double last_pr_first_base = 1;
	double last_log_pr = cur_log_pr;

	//print "$cur_pr_first_base $cur_log_pr\n" if ($verbose);

	while (cur_log_pr >= last_log_pr)
	{
		last_log_pr = cur_log_pr;
		last_pr_first_base = cur_pr_first_base;
    if (cur_pr_first_base < 0) break;

		cur_pr_first_base -= 0.001;
		cur_log_pr = calculate_two_base_model_log10_likelihood(best_base_char, second_best_base_char, pdata, cur_pr_first_base);
		//print "$cur_pr_first_base $cur_log_pr\n" if ($verbose);
	}
	
	return make_pair(last_pr_first_base, last_log_pr);
}

/*! Calculate the likelihood of a mixture model of two bases leading to the observed read bases.
 */
double breseq::identify_mutations_pileup::calculate_two_base_model_log10_likelihood (base_char best_base_char, base_char second_best_base_char, const vector<polymorphism_data>& pdata, double best_base_freq)
{
	double log10_likelihood = 0;	
	
  for(vector<polymorphism_data>::const_iterator it=pdata.begin(); it<pdata.end(); ++it) {
  
    //## the first value is pr_base, second is pr_not_base
    double best_base_log10pr;
    double second_best_base_log10pr;
    if (m_use_cErrorTable) {
      covariate_values_t this_cv = it->_cv;
      
      if (it->_strand == -1) {
        this_cv.obs_base() = complement_base_index(this_cv.obs_base());
      }
      
      if(it->_strand == 1) {
        this_cv.ref_base() = basechar2index(best_base_char);
      } else {
        this_cv.ref_base() = basechar2index(complement_base_char(best_base_char));
      }
      best_base_log10pr = m_error_table.get_log10_prob(this_cv);
      
      if(it->_strand == 1) {
        this_cv.ref_base() = basechar2index(second_best_base_char);
      } else {
        this_cv.ref_base() = basechar2index(complement_base_char(second_best_base_char));
      }
      second_best_base_log10pr = m_error_table.get_log10_prob(this_cv);

    } else {
    
      string best_base_key;
      string second_best_base_key;
      
      if(it->_strand == 1) {
        best_base_key += best_base_char;
        best_base_key += it->_base_char;
        second_best_base_key += second_best_base_char;
        second_best_base_key += it->_base_char;
      } else {
        best_base_key += complement_base_char(best_base_char); 
        best_base_key += complement_base_char(it->_base_char);
        second_best_base_key += complement_base_char(second_best_base_char);
        second_best_base_key += complement_base_char(it->_base_char);
      }
      best_base_log10pr = _ecr.log10_correct_rates(it->_fastq_file_index, it->_quality, best_base_key);
      second_best_base_log10pr = _ecr.log10_correct_rates(it->_fastq_file_index, it->_quality, second_best_base_key);
    }
    
    //debug output
    //cerr << "Base in Read: " << it->base << " Read Strand: " << it->strand << endl;
    //cerr << "Best Base: " << best_base << " Key: " << best_base_key << " Chance of Observing: " << pow(10,best_base_log10pr) << endl;
    //cerr << "Second Best Base: " << second_best_base << " Key: " << second_best_base_key << " Chance of Observing: " << pow(10,second_best_base_log10pr) << endl;

    double pr_ref_base_given_obs = best_base_freq * pow(10, best_base_log10pr) + (1-best_base_freq) * pow(10, second_best_base_log10pr);

    log10_likelihood += log(pr_ref_base_given_obs);		
  }

	log10_likelihood /= log(10);
  
  //debug output
  /*
  cerr << "Best Base: " << best_base << " Second Best Base: " << second_best_base << " Fraction Best Base: " << best_base_freq << " Log10 Likelihood " << log10_likelihood << endl;
  */
	return log10_likelihood;
}

breseq::cDiscreteSNPCaller::cDiscreteSNPCaller(uint8_t ploidy) 
: _observations(0), _ploidy(ploidy)
{
  assert(ploidy==1); // only the haploid version is implemented

  double _log10_prior_probability = 0;

  //create all of the states and initialize priors
  for (int i=0; i<base_list_size; i++) {
    _log10_priors.push_back(_log10_prior_probability);
    _log10_probabilities.push_back(_log10_prior_probability);
  }
  _normalized = false;
}

// obs_base is a BAM style base when input
void breseq::cDiscreteSNPCaller::update(base_bam obs_base_bam, uint8_t obs_quality, bool obs_top_strand, int32_t fastq_file_index, error_count_results &ecr)
{  
  //update probabilities give observation using Bayes rule
  for (int i=0; i<base_list_size; i++) {
    string base_key;
    if (obs_top_strand) {
      base_key += base_char_list[i]; 
      base_key += basebam2char(obs_base_bam);    
    } else {
      base_key += complement_base_char(base_char_list[i]); 
      base_key += basebam2char(complement_base_bam(obs_base_bam));
    }
    _log10_probabilities[i] += ecr.log10_correct_rates(fastq_file_index, obs_quality, base_key);
    
    //cerr << "  " << fastq_file_index << " " << base_key << " " << (int)obs_quality << " " << ecr.log10_correct_rates(fastq_file_index, obs_quality, base_key) << endl;
  }
  _normalized = false;
  _observations++;
}

void breseq::cDiscreteSNPCaller::update(const covariate_values_t& cv, bool obs_top_strand, cErrorTable& et) {

  covariate_values_t this_cv = cv;
  //update probabilities give observation using Bayes rule
  if (!obs_top_strand) {
    this_cv.obs_base() = complement_base_index(this_cv.obs_base()); 
  }

  for (int i=0; i<base_list_size; i++) {
    this_cv.ref_base() = i;
    if (!obs_top_strand) {
      this_cv.ref_base() = complement_base_index(this_cv.ref_base());
    }    
    _log10_probabilities[i] += et.get_log10_prob(this_cv);
    
    
    //cout << " " << i << " " << obs_top_strand << " ? " << cv.obs_base() << " " << cv.ref_base()<< " " << et.get_log10_prob(cv) << endl;
  }
  
  _normalized = false;
  _observations++;
}


pair<uint8_t,double> breseq::cDiscreteSNPCaller::get_prediction()
{
  
  if (_observations == 0) {
    //Best base is 'N' and E-value is NAN
    return make_pair('N', NAN);
  }


  // which is the largest log probability?
  vector<double>::iterator max = max_element(_log10_probabilities.begin(), _log10_probabilities.end());
  double max_log10_probability = *max;
  int max_index = distance(_log10_probabilities.begin(), max);
  uint8_t best_base = base_char_list[max_index];

  // we want to normalize the probabilities, but avoid floating point errors
  if (!_normalized) {
  
    double total_offset_probability = 0;
    for(vector<double>::iterator i=_log10_probabilities.begin(); i!=_log10_probabilities.end(); ++i) {
      total_offset_probability += pow(10, *i - max_log10_probability);
    }
    double total_log10_offset_probability = log10(total_offset_probability) + max_log10_probability;
    
    for(vector<double>::iterator i=_log10_probabilities.begin(); i!=_log10_probabilities.end(); ++i) {
      *i = *i - total_log10_offset_probability;
    }    
    _normalized = true;
  }

  // ignore the highest probability and combine all others to get the total error rate of this base call
  double second_best_log10_error_probability = numeric_limits<double>::quiet_NaN();
  vector<double>::iterator second_best;
  for(vector<double>::iterator i=_log10_probabilities.begin(); i!=_log10_probabilities.end(); ++i) {
    if (i != max) {
      if (isnan(second_best_log10_error_probability) || (*i > second_best_log10_error_probability) ) {
        second_best_log10_error_probability = *i;
        second_best = i;
      }
    }
  }
  
  double total_error_probability = 0.0;
  for(vector<double>::iterator i=_log10_probabilities.begin(); i!=_log10_probabilities.end(); ++i) {
    if ((i != max) ) {
      total_error_probability += pow(10, *i - second_best_log10_error_probability);
    }
  }
  double total_log10_error_probability = log10(total_error_probability) + second_best_log10_error_probability;

  // this second best base isn't what we want for polymorphism prediction
  // because the second best may have less coverage than the third best...
  //int second_best_index = distance(_log10_probabilities.begin(), second_best);
  //uint8_t second_best_base = base_list[second_best_index];
  
  return make_pair(best_base, -total_log10_error_probability);
}

