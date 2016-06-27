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

#include "libbreseq/common.h"
#include "libbreseq/pileup.h"
#include "libbreseq/identify_mutations.h"
#include "libbreseq/error_count.h"

using namespace std;

namespace breseq {

/*! Convenience wrapper around the identify_mutations_pileup class.
 */
  
void identify_mutations(
                const Settings& settings,
                const Summary& summary,
								const string& bam,
								const string& fasta,
								const string& gd_file,
                const vector<double>& deletion_propagation_cutoff,
                const vector<double>& deletion_seed_cutoffs,
								double mutation_cutoff,
								double polymorphism_cutoff,
								double polymorphism_frequency_cutoff,
                double polymorphism_precision_decimal,
                uint32_t polymorphism_precision_places,
								bool print_per_position_file
 ) {
                                                                                            
	// do the mutation identification:
	identify_mutations_pileup imp(
                settings,
                summary,
								bam,
								fasta,
                deletion_propagation_cutoff,
                deletion_seed_cutoffs,
								mutation_cutoff,
								polymorphism_cutoff,
								polymorphism_frequency_cutoff,
                polymorphism_precision_decimal,
                polymorphism_precision_places,
								print_per_position_file
							);
	imp.do_pileup(settings.call_mutations_seq_id_set());
  imp.write_gd(gd_file);
}

/*! Chooses which RA predictions to reject and for what reasons
 */
enum ePredictionType {unknown, consensus, polymorphism};


bool rejected_RA_polymorphism_bias(cDiffEntry& ra,
                                       cReferenceSequences& ref_seq_info,
                                       const Settings& settings
                                       )
{
  (void)ref_seq_info;
  bool rejected = false;


  ////// Strand and base quality score biases
  
  // Test quality score bias
  if (settings.polymorphism_bias_p_value_cutoff) {
    if (ra.entry_exists("ks_quality_p_value")) {
      if (from_string<double>(ra["ks_quality_p_value"]) < settings.polymorphism_bias_p_value_cutoff) {
        rejected = true;
        ra.add_reject_reason("KS_BASE_QUALITY");
      }
    }
  }
  
  if (settings.polymorphism_bias_p_value_cutoff) {
    if (ra.entry_exists("fisher_strand_p_value")) {
      if (from_string<double>(ra["fisher_strand_p_value"]) < settings.polymorphism_bias_p_value_cutoff) {
        rejected = true;
        ra.add_reject_reason("FISHER_STRAND");
      }
    }
  }
  
  return rejected;
}
  
bool rejected_RA_polymorphism_coverage(cDiffEntry& ra,
                               cReferenceSequences& ref_seq_info,
                               const Settings& settings
                               )
{
  (void)ref_seq_info;
  bool rejected = false;

  // Minimum coverage on both strands for both reference and new allele
  if (settings.polymorphism_minimum_new_coverage_each_strand > 0) {
    double polymorphism_coverage_limit_both_bases = settings.polymorphism_minimum_new_coverage_each_strand;
    bool passed = true;
    vector<string> top_bot;
    double top;
    double bot;
    
    top_bot = split(ra[MAJOR_COV], "/");
    top = from_string<double>(top_bot[0]);
    bot = from_string<double>(top_bot[1]);
    passed = passed && (top >= polymorphism_coverage_limit_both_bases);
    passed = passed && (bot >= polymorphism_coverage_limit_both_bases);
    
    top_bot = split(ra[MINOR_COV], "/");
    top = from_string<double>(top_bot[0]);
    bot = from_string<double>(top_bot[1]);
    passed = passed && (top >= polymorphism_coverage_limit_both_bases);
    passed = passed && (bot >= polymorphism_coverage_limit_both_bases);
    
    if (!passed) {
      rejected = true;
      ra.add_reject_reason("STRAND_COVERAGE");
    }
  }
  
  return rejected;
}
  
bool rejected_RA_indel_homopolymer(cDiffEntry& ra,
                              cReferenceSequences& ref_seq_info,
                              const Settings& settings
                               )
{
  
  bool rejected = false;
  
  ////// Optionally, ignore if in a homopolymer stretch longer than this
  // There are several cases here that may be of interest for removing false-positives
  // 1) indel in what was originally a homopolymer adding a base that is the same
  // 2) A mutation in the middle of a stretch of one base to what is in the rest of that stretch
  //    TTTTTTATTTTTT
  //    TTTTTTTTTTTTT
  //
  // First case,
  // 1) indel in what was originally a homopolymer adding a base that is the same
  if (settings.polymorphism_reject_indel_homopolymer_length && ((ra[MAJOR_BASE] == ".") || (ra[MINOR_BASE] == ".")))
  {
    // Code needs to be robust to the mutation being at the beginning
    // OR end of the homopolymer tract
    
    string seq_id = ra["seq_id"];
    int32_t mut_pos = from_string<uint32_t>(ra["position"]);
    bool is_insertion = (from_string<int32_t>(ra[INSERT_POSITION]) > 0);
    string mut_base = (ra[MAJOR_BASE] == ".") ? ra[MINOR_BASE] : ra[MAJOR_BASE];
    int32_t homopolymer_length = 0;
    
    // If we are an insertion, we need to start on the base before or the base after,
    // depending on which one we match (if we match neither than our length is zero
    bool no_match = false;
    if (is_insertion) {
      
      // This is the base before
      if (ref_seq_info.get_sequence_1(seq_id, mut_pos, mut_pos) != mut_base) {
        
        // No match? -- Now check the base after
        mut_pos++;
        if (mut_pos > static_cast<int32_t>(ref_seq_info.get_sequence_length(seq_id))) {
          no_match = true;
        } else if (ref_seq_info.get_sequence_1(seq_id, mut_pos, mut_pos) != mut_base) {
          no_match = true;
        }
      }
    }
    
    // extend match
    if (!no_match) {
      
      // check before
      int32_t start_pos = mut_pos - 1;
      while ( (ref_seq_info.get_sequence_1(seq_id, start_pos, start_pos) == mut_base) && (start_pos > 0) ) {
        start_pos--;
      }
      start_pos++;
      
      int32_t end_pos = mut_pos + 1;
      while ( (ref_seq_info.get_sequence_1(seq_id, end_pos, end_pos) == mut_base) && (end_pos <= static_cast<int32_t>(ref_seq_info.get_sequence_length(seq_id)) ) ) {
        end_pos++;
      }
      end_pos--;
      homopolymer_length = end_pos - start_pos + 1;
    }
    
    if (homopolymer_length >= static_cast<int32_t>(settings.polymorphism_reject_indel_homopolymer_length))
    {
      rejected = true;
      ra.add_reject_reason("INDEL_HOMOPOLYMER");
    }
  }
  
  // Second case
  // 2) A mutation in the middle of a stretch of one base to what is in the rest of that stretch
  //    TTTTTTATTTTTT
  //    TTTTTTTTTTTTT
  if (settings.polymorphism_reject_surrounding_homopolymer_length && (ra[MAJOR_BASE] != ".") && (ra[MINOR_BASE] != "."))    {
    
    string seq_id = ra["seq_id"];
    int32_t mut_pos = from_string<int32_t>(ra["position"]);
    string ref_base = ref_seq_info.get_sequence_1(seq_id, mut_pos, mut_pos);
    string mut_base = (ra[MAJOR_BASE] != ref_base) ? ra[MAJOR_BASE] : ra[MINOR_BASE];
    
    // determine the start and end of homopolymer created by substitution
    int32_t start_pos = mut_pos - 1;
    while (start_pos >= 1)
    {
      if (ref_seq_info.get_sequence_1(seq_id, start_pos, start_pos) != mut_base) {
        break;
      }
      start_pos--;
    }
    start_pos++;
    
    int32_t end_pos = mut_pos + 1;
    while (end_pos <= static_cast<int32_t>(ref_seq_info.get_sequence_length(seq_id)))
    {
      if (ref_seq_info.get_sequence_1(seq_id, end_pos, end_pos) != mut_base) {
        break;
      }
      end_pos++;
    }
    end_pos--;
    
    // check bounds
    if ( (start_pos < mut_pos) && (end_pos > mut_pos) && (end_pos - start_pos + 1 >= static_cast<int32_t>(settings.polymorphism_reject_surrounding_homopolymer_length)) ) {
      rejected = true;
      ra.add_reject_reason("SURROUNDING_HOMOPOLYMER");
    }
  }
  ////// END checking on homopolymers
  
  if (settings.no_indel_polymorphisms && ((ra[MAJOR_BASE] == ".") || (ra[MINOR_BASE] == ".")))
  {
    rejected = true;
    ra.add_reject_reason("INDEL_POLYMORPHISM");
  }
  
  return rejected;
}
  
// Returns whether it should be deleted
bool test_RA_evidence_CONSENSUS_mode(
                                     cDiffEntry& ra,
                                     cReferenceSequences& ref_seq_info,
                                     const Settings& settings
                                     )
{
  
  // Decide if we are a polymorphism (or mixed base) prediction or a consensus prediction
  ePredictionType prediction(unknown);
  
  double consensus_score = from_string<double>(ra[CONSENSUS_SCORE]);
  double polymorphism_score = (ra[POLYMORPHISM_SCORE]=="NA") ? numeric_limits<double>::quiet_NaN() : from_string<double>(ra[POLYMORPHISM_SCORE]);
  double variant_frequency = from_string<double>(ra[POLYMORPHISM_FREQUENCY]);
  
  // Score must be below one of the thresholds
  if (consensus_score >= settings.mutation_log10_e_value_cutoff) {
    prediction = consensus;
  } else if (polymorphism_score >= settings.polymorphism_log10_e_value_cutoff) {
    prediction = polymorphism;
  } else {
    return true;
  }
  
  // Drop down to a polymorphism if we don't pass consensus frequency criterion
  if ((prediction == consensus) && (variant_frequency < settings.consensus_frequency_cutoff - settings.polymorphism_precision_decimal)) {
    prediction = polymorphism;
  }
  
  // Drop down to a polymorphism if we don't pass the consensus strand criterion
  if ( (prediction == consensus) && (settings.consensus_minimum_new_coverage_each_strand > 0) ) {
    vector<string> top_bot = split(ra[MAJOR_COV], "/");
    if ( (from_string<double>(top_bot[0]) < settings.consensus_minimum_new_coverage_each_strand) ||
        (from_string<double>(top_bot[1]) < settings.consensus_minimum_new_coverage_each_strand) ) {
      prediction = polymorphism;
    }
  }
  
  // Bail now if consensus
  if (prediction == consensus) {
    ra[PREDICTION] = "consensus";
    ra[FREQUENCY] = (ra[REF_BASE] == ra[MAJOR_BASE]) ? "0" : "1";
    
    // and delete if we are just the reference base!
    return (ra[REF_BASE] == ra[MAJOR_BASE]);
  }
  
  // Code only reaaches here for non-consensus predictions
  
  // Set us up as a potential polymorphism
  ra[PREDICTION] = "polymorphism";
  ra[FREQUENCY] = ra[POLYMORPHISM_FREQUENCY];
  
  // Perform further checks for polymorphisms
  bool rejected_RA_polymorphism_frequency = false;
  if (variant_frequency > 1.0 - settings.polymorphism_frequency_cutoff - settings.polymorphism_precision_decimal) {
    rejected_RA_polymorphism_frequency = true;
    ra.add_reject_reason("FREQUENCY_CUTOFF");
  } else if (variant_frequency < settings.polymorphism_frequency_cutoff) {
    rejected_RA_polymorphism_frequency = true;
    ra.add_reject_reason("FREQUENCY_CUTOFF");
  }
  if (rejected_RA_polymorphism_frequency)
    return true;
  
  if (rejected_RA_polymorphism_bias(ra, ref_seq_info, settings))
    return true;
  
  if (rejected_RA_polymorphism_coverage(ra, ref_seq_info, settings))
    return true;
  
  if (rejected_RA_indel_homopolymer(ra, ref_seq_info, settings))
    return true;
  
  return false;
}
  
// Returns whether it should be deleted
bool test_RA_evidence_POLYMORPHISM_mode(
                                     cDiffEntry& ra,
                                     cReferenceSequences& ref_seq_info,
                                     const Settings& settings
                                     )
{
  
  // Decide if we are a polymorphism (or mixed base) prediction or a consensus prediction
  ePredictionType prediction(unknown);
  
  double consensus_score = from_string<double>(ra[CONSENSUS_SCORE]);
  double polymorphism_score = (ra[POLYMORPHISM_SCORE]=="NA") ? numeric_limits<double>::quiet_NaN() : from_string<double>(ra[POLYMORPHISM_SCORE]);
  double variant_frequency = from_string<double>(ra[POLYMORPHISM_FREQUENCY]);
  
  // Set up polymorphism as default
  ra[PREDICTION] = "polymorphism";
  ra[FREQUENCY] = ra[POLYMORPHISM_FREQUENCY];
  
  // Score must be below one of the thresholds
  if (polymorphism_score >= settings.polymorphism_log10_e_value_cutoff) {
    prediction = polymorphism;
  } else if (consensus_score >= settings.mutation_log10_e_value_cutoff) {
    prediction = consensus;
  } else {
    return true;
  }
  
  bool tested_indel_homopolymer = false;
  bool failed_indel_homopolymer_test = false;
  string saved_polymorphism_reject;
  if (prediction == polymorphism) {
    
    // Perform further checks for polymorphisms
    bool rejected_RA_polymorphism_frequency = false;
    if (variant_frequency > 1.0 - settings.polymorphism_frequency_cutoff - settings.polymorphism_precision_decimal) {
      rejected_RA_polymorphism_frequency = true;
      ra.add_reject_reason("FREQUENCY_CUTOFF");
    } else if (variant_frequency < settings.polymorphism_frequency_cutoff) {
      rejected_RA_polymorphism_frequency = true;
      ra.add_reject_reason("FREQUENCY_CUTOFF");
    }
    if (rejected_RA_polymorphism_frequency)
      prediction = consensus;
    
    if (prediction == polymorphism) {
      if (rejected_RA_polymorphism_bias(ra, ref_seq_info, settings))
        prediction = consensus;
    }
    if (prediction == polymorphism) {
      if (rejected_RA_polymorphism_coverage(ra, ref_seq_info, settings))
        prediction = consensus;
    }
    if (prediction == polymorphism) {
      tested_indel_homopolymer = true;
      failed_indel_homopolymer_test = rejected_RA_indel_homopolymer(ra, ref_seq_info, settings);
      if (failed_indel_homopolymer_test)
        prediction = consensus;
    }
    
    // If we're still a good polymorphism -- fill out values and return
    if (prediction == polymorphism) {
      return false;
    }
    
    saved_polymorphism_reject = ra[REJECT];
    ra.clear_reject_reasons();
  }
  
  // If we get here, we are back to testing for a consensus
  
  // Drop back to a polymorphism if we don't pass consensus frequency criterion
  if ((prediction == consensus) && (variant_frequency < settings.consensus_frequency_cutoff - settings.polymorphism_precision_decimal)) {
    prediction = polymorphism;
    ra.add_reject_reason("FREQUENCY_CUTOFF");
  }
  
  // Drop down to a rejected polymorphism if we don't pass the consensus strand criterion
  if ( (prediction == consensus) && (settings.consensus_minimum_new_coverage_each_strand > 0) ) {
    vector<string> top_bot = split(ra[MAJOR_COV], "/");
    if ( (from_string<double>(top_bot[0]) < settings.consensus_minimum_new_coverage_each_strand) ||
        (from_string<double>(top_bot[1]) < settings.consensus_minimum_new_coverage_each_strand) ) {
      prediction = polymorphism;
      ra.add_reject_reason("STRAND_COVERAGE");
    }
  }
  
  // We can only go back to a consensus prediction if we are not user_evidence
  if (!ra.entry_exists("user_defined")) {
  
    // Bail now if still consensus
    if (prediction == consensus) {
      
      ra[PREDICTION] = "consensus";
      ra[FREQUENCY] = "1";
      
      // and delete if we are just the reference base!
      return (ra[REF_BASE] == ra[MAJOR_BASE]);
    }
    
    // We got kicked out of consensus mode, but we don't pass the polymorphism
    // score cutoff, so we are now a rejected consensus mutation
    // This is the only case that inherits the rejections from consensus tests above
    if (polymorphism_score < settings.polymorphism_log10_e_value_cutoff) {
      
      ra[PREDICTION] = "consensus";
      ra[FREQUENCY] = "1";
      
      // and delete if we are just the reference base!
      return (ra[REF_BASE] == ra[MAJOR_BASE]);
      
    }
  }
  
  // We are back to a polymorphism, albeit a rejected one
  ra[REJECT] = saved_polymorphism_reject;
  
  // Must pass this test in order to be retained as a marginal prediction
  // Saving previously calculated value prevents doing a computationally intensive test twice.
  if ( (tested_indel_homopolymer && failed_indel_homopolymer_test)
      || rejected_RA_indel_homopolymer(ra, ref_seq_info, settings)) {
    return true;
  }

  return false;
  
}
  
void test_RA_evidence(
                        cGenomeDiff& gd,
                        cReferenceSequences& ref_seq_info,
                        const Settings& settings
                        )
{
  
  // Assumes these entries are present initially
  //  * REF_BASE
  //  * CONSENSUS_SCORE is present for all entries        = maximumum likelihood frequency for 100% mutated model
  //  * POLYMORPHISM_FREQUENCY is present for all entries = maximumum likelihood frequency for mixed model
  //  * QUALITY is present for all entries                = consensus model
  //  * POLYMORPHISM_SCORE may or may not be present      = mixed model
  //  * NEW_COV, REF_COV present for all entries
  
  diff_entry_list_t list = gd.get_list();
  // Note nonstandard non-increment, since we remove items
  for (diff_entry_list_t::iterator it = list.begin(); it != list.end(); ) {
    
    cDiffEntry& ra = **it;
    if (ra._type != RA) {
      it++;
      continue;
    }
    
    // At this point, both consensus and polymorphism quality may be > 0,
    // we see if either is below the significance + frequency cutoff first.

    if (!ra.entry_exists(REF_BASE)) {
      ERROR("Expected field 'ref_base' in evidence item\n" + ra.as_string());
    }
    if (!ra.entry_exists(CONSENSUS_SCORE)) {
      ERROR("Expected field 'consensus_score' in evidence item\n" + ra.as_string());
    }
    if (!ra.entry_exists(POLYMORPHISM_SCORE)) {
      ERROR("Expected field 'polymorphism_score' in evidence item\n" + ra.as_string());
    }
    if (!ra.entry_exists(POLYMORPHISM_FREQUENCY)) {
      ERROR("Expected field '" + cString(POLYMORPHISM_FREQUENCY) + "' in evidence item\n" + ra.as_string());
    }
    
    bool delete_entry = false;
    

    if (settings.polymorphism_prediction) {
      delete_entry = test_RA_evidence_POLYMORPHISM_mode(ra, ref_seq_info, settings);
    } else {
      delete_entry = test_RA_evidence_CONSENSUS_mode(ra, ref_seq_info, settings);
    }
  
    // User evidence items should be treated normally (to define polymorphism vs. consensus predictions)
    // but they should not be rejected, and they should not be deleted.
    if (ra.entry_exists("user_defined")) {
      delete_entry=false;
    }

    if (delete_entry) {
      it = list.erase(it);
    } else {
      it++;
    }
  
  } // end loop over ra evidence items

  // Resave the list (from which we may have removed items)
  gd.set_list(list);
}


/*! Constructor.
 */
identify_mutations_pileup::identify_mutations_pileup(
                              const Settings& settings,
                              const Summary& summary,
															const string& bam,
															const string& fasta,
                              const vector<double>& deletion_propagation_cutoffs,
                              const vector<double>& deletion_seed_cutoffs,
															double consensus_score_cutoff,
															double polymorphism_score_cutoff,
															double polymorphism_frequency_cutoff,
                              double polymorphism_precision_decimal,
                              uint32_t polymorphism_precision_places,
															bool print_per_position_file
                                                            )
: pileup_base(bam, fasta)
, _settings(settings)
, _gd()
, _deletion_seed_cutoffs(deletion_seed_cutoffs)
, _deletion_propagation_cutoffs(deletion_propagation_cutoffs)
, _consensus_score_cutoff(consensus_score_cutoff)
, _polymorphism_score_cutoff(polymorphism_score_cutoff)
, _polymorphism_frequency_cutoff(polymorphism_frequency_cutoff)
, _polymorphism_precision_decimal(polymorphism_precision_decimal)
, _polymorphism_precision_places(polymorphism_precision_places)
, _log10_ref_length(0)
, _snp_caller("haploid", summary.sequence_conversion.total_reference_sequence_length)
, _this_deletion_reaches_seed_value(false)
, _this_deletion_redundant_reached_zero(false)
, _last_position_coverage_printed(0)
, _print_per_position_file(print_per_position_file)
{
	
  // remove once used
  (void)settings;
  
  set_print_progress(true);
  
  ASSERT(m_bam->header->n_targets == (int32_t)_deletion_propagation_cutoffs.size(), 
         "Number of targets in BAM file [" + to_string(m_bam->header->n_targets) + "] " +
         "does not match + number in cutoff table [" + to_string(_deletion_propagation_cutoffs.size()) + "]."
         );
  ASSERT(m_bam->header->n_targets == (int32_t)_deletion_seed_cutoffs.size(),
         "Number of targets in BAM file [" + to_string(m_bam->header->n_targets) + "] " +
         "does not match + number in cutoff table [" + to_string(_deletion_propagation_cutoffs.size()) + "]."
         );
    
	// reserve enough space for the sequence info:
	_seq_info.resize(m_bam->header->n_targets);
	
  // tally up the reference lengths from the bam file
	for(int i=0; i<m_bam->header->n_targets; ++i) {
		_log10_ref_length += static_cast<double>(m_bam->header->target_len[i]);
	}
	assert(_log10_ref_length != 0);
	_log10_ref_length = log10(_log10_ref_length);
  
  // are we printing detailed coverage information?
  _print_coverage_data = true;
  
  // load the error table file and convert back to probabilities
  _error_table.read_log10_prob_table(settings.error_rates_file_name);
  _error_table.log10_prob_to_prob();
  
  if (_print_per_position_file) {
    _per_position_file.open(_settings.mutation_identification_per_position_file_name.c_str());
  }
  
  // load the list of user RA evidence
  if (_settings.user_evidence_genome_diff_file_name != "") {
    load_user_ra_evidence_from_gd();
  }
  
}

void identify_mutations_pileup::load_user_ra_evidence_from_gd()
{
  cGenomeDiff gd(_settings.user_evidence_genome_diff_file_name);
  gd.sort(); // very important to be sorted!!
  
  // get rid of any previous annotations
  gd.strip_to_spec();
  
  // Return just RA entries
  _user_evidence_ra_list = gd.get_list(make_vector<gd_entry_type>(RA));
}

/*! Destructor.
 */
identify_mutations_pileup::~identify_mutations_pileup()
{
}

/*! Called for each alignment.
 */
void identify_mutations_pileup::pileup_callback(const pileup& p) {
  
  bool verbose = false;
  ASSERT(p.target() < _seq_info.size(), "Unknown target id: " + to_string<uint32_t>(p.target()));
  if (verbose) cout << "Target id: " << p.target() << " position: " << p.position_1() << endl;

  _this_deletion_propagation_cutoff = _deletion_propagation_cutoffs[p.target()];
  _this_deletion_seed_cutoff = _deletion_seed_cutoffs[p.target()];
  
  // if the propagation cutoff is zero then the coverage distribution failed
  if (_this_deletion_propagation_cutoff < 0.0) return;
  
  uint32_t position = p.position_1();
  
	int32_t insert_count=-1;
	bool next_insert_count_exists=true;
	
  // User RA evidence
  // We need to continue looping in some cases where there is user evidence RA
  // but no reads actually have a base at that insert position. Here's how:
  int32_t force_insert_count_max = 0;
  if (_user_evidence_ra_list.size()) {
    
    diff_entry_list_t::iterator user_ra = _user_evidence_ra_list.begin();
    
    while ( (user_ra != _user_evidence_ra_list.end()) && (from_string<uint32_t>((**user_ra)[POSITION]) == position)) {
      force_insert_count_max = from_string<int32_t>((**user_ra)[INSERT_POSITION]);
      user_ra++;
    }
  }
  
  // BIG outer loop to traverse this position AND any positions after this one
  // that are inserted relative to the reference genome
	while(next_insert_count_exists || (insert_count < force_insert_count_max)) {
		++insert_count;
    next_insert_count_exists = false;
		
		base_char ref_base_char = '.';
		if(!insert_count) {
			ref_base_char = p.reference_base_char_1(position);
		}
		
		//## zero out the info about this position
		position_base_info pos_info;
    position_base_info redundant_pos_info;
		
		//## keep track of coverage for deletion prediction
		position_coverage this_position_coverage;
		bool this_position_unique_only_coverage=true;
		
    //## reset SNP caller
    _snp_caller.reset(basechar2index(ref_base_char));
        
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
      bool on_insert_position_past_base = true; 
      if(indel >= insert_count) {
        read_base_bam = i->read_base_bam_0(i->query_position_0() + insert_count);
        on_insert_position_past_base = false; 
      }
            
      //## don't use bases without qualities!!
      if(_base_bam_is_N(read_base_bam)) continue;
      
      //## gather information about the aligned base
      int32_t redundancy = i->redundancy();
      int32_t fastq_file_index = i->fastq_file_index();
      int strand = i->strand();
      bool trimmed = i->is_trimmed(on_insert_position_past_base);
      
      //##### update coverage if this is not a deletion in read relative to reference
      //### note that we count trimmed reads here, but not when looking for short indel mutations...	
      
      if(redundancy == 1) {
        //## keep track of unique coverage	
        ++this_position_coverage.unique[1+strand];
        
        //## we don't continue to consider further insertions relative
        //## to the reference unless uniquely aligned reads have them 
        if(indel > insert_count)
          next_insert_count_exists = true;
        
      } else {		
        //## mark that this position has some non-unique coverage
        this_position_unique_only_coverage = false;
        //## keep track of redundant coverage
        this_position_coverage.redundant[1+strand] += 1.0/redundancy;
        ++this_position_coverage.raw_redundant[1+strand];
        
        if (!trimmed) 
          ++redundant_pos_info[basebam2char(read_base_bam)][1+strand];
      }
      
      //## don't use information from trimmed or redundant reads
      //## to predict base substitutions and short indels!!
      if(trimmed || (redundancy > 1)) {
        continue;
      }

			//##### deal with base calls
      //cerr << "POSITION:" << position << endl;
      
      covariate_values_t cv; 

      bool is_ok = _error_table.alignment_position_to_covariates(*i, insert_count, cv);
      //cv.obs_base is still not a char here...
      
      if (is_ok)  {
      
        if (cv.quality() < _settings.base_quality_cutoff) {
          //cout << cv.quality()  << endl;
          continue;
        }

        //## this is the coverage for SNP counts, tabulate AFTER skipping trimmed reads
        ++pos_info[baseindex2char(cv.obs_base())][1+strand];
        
        //##### this is for polymorphism prediction and making strings
        pdata.push_back(polymorphism_data(baseindex2char(cv.obs_base()),cv.quality(),i->strand(),cv.read_set(), cv));
        
        //cerr << " " << cv.obs_base() << " " << (char)ref_base << endl;

        _snp_caller.update(cv, strand == 1, i->mapping_quality(), _error_table);
      }
		} // end for-each read
		
    
    //#############################
		//## PER POSITION/INSERT COUNT
    //#############################
    
		//#sum up coverage observations	
		this_position_coverage.sum();
				
		//#we are trying to find the base with the most support
    cSNPCall snp_call = _snp_caller.get_prediction();
    
    base_char best_base_char('N');
    double consensus_bonferroni_score(numeric_limits<double>::quiet_NaN());
    double polymorphism_bonferroni_score(numeric_limits<double>::quiet_NaN());

    // SNP caller returns one genotype
    best_base_char = snp_call.genotype[0];
    // Here is where we convert to a Bonferroni type correction for phred score from SNP caller
    consensus_bonferroni_score = snp_call.score - (_log10_ref_length);
    
    //Do we predict a base at this position?
    bool base_predicted=false;
    if(consensus_bonferroni_score >= _consensus_score_cutoff) {
      base_predicted = true;
    }
    
    //cout << position << " " << best_base_char << " " << e_value_call << " " << (base_predicted ? "T" : "F") << endl;
    
		int total_cov[3]={0,0,0}; // triple, same as above
    
    //// BEGIN Print per-position output file
		ostringstream line;
		if (_print_per_position_file) {
      line << position << " " << insert_count << " " << ref_base_char << " " << consensus_bonferroni_score;
		}
    
    for(size_t j=0; j<base_list_size; ++j) {
			double top_cov = pos_info[base_char_list[j]][2];
			double bot_cov = pos_info[base_char_list[j]][0];
			total_cov[2] += (int)round(top_cov);
			total_cov[0] += (int)round(bot_cov);
    }
    
    //// Summing coverage here should be moved to where coverage is updated?
    
    // Debug: print additional information to file.
    if (_print_per_position_file) {
       
      // Print unique bases
      for(size_t j=0; j<base_list_size; ++j) {
        double top_cov = pos_info[base_char_list[j]][2];
        double bot_cov = pos_info[base_char_list[j]][0];
        if (_print_per_position_file) {
          line << " " << base_char_list[j] << " (" << bot_cov << "/" << top_cov << ")";
        }
      }
    
      // Print redundant bases
      for(size_t j=0; j<base_list_size; ++j) {
        double top_cov = redundant_pos_info[base_char_list[j]][2];
        double bot_cov = redundant_pos_info[base_char_list[j]][0];
        if (_print_per_position_file) {
          line << " r" << base_char_list[j] << " (" << bot_cov << "/" << top_cov << ")";
        }
      }
    }
    
    // Finally print line (kept separate from above because the line is
    // added to at various points in the code).
    if (_print_per_position_file) {
      _per_position_file << line.str() << endl;
    }
    //// END Per-position output file
		
		//###
		//## DELETION DELETION DELETION
		//###
		
    if(!_settings.no_deletion_prediction && (insert_count == 0))
      check_deletion_completion(position, p.target(), this_position_coverage, consensus_bonferroni_score);
		
		//###
		//## POLYMORPHISM POLYMORPHISM POLYMORPHISM
		//###								

		bool passed_as_polymorphism_prediction(false);
    polymorphism_prediction ppred;
    base_char second_best_base_char('N');

    cDiffEntry mut(RA);
    
    //## evaluate whether to call an actual mutation!
    // -- note that we accept > 0 and only reject later
    // so that these can potentially make it into the marginal data
    bool passed_as_consensus_prediction = (best_base_char != ref_base_char) && (!isnan(consensus_bonferroni_score) && (consensus_bonferroni_score > 0));
    
    // Find the bases with the highest and second highest coverage
    // We only predict polymorphisms involving these 'major' and 'minor' alleles
    base_index best_base_index(base_list_N_index);
    int best_base_coverage(0);
    base_index second_best_base_index(base_list_N_index);
    int second_best_base_coverage(0);

    vector<double> snp_probs = _snp_caller.get_genotype_log10_probabilities();
            
    for (uint8_t i=0; i<base_list_size; i++) {
      base_char this_base_char = base_char_list[i];
      base_index this_base_index = i;
      int this_base_coverage = pos_info[this_base_char][0] + pos_info[this_base_char][2];
      
      if (this_base_coverage==0) continue;
      
      // if better coverage or tied in coverage and better probability
      if ((this_base_coverage > best_base_coverage) ||
        ((this_base_coverage == best_base_coverage) && ((best_base_index==base_list_N_index) || (snp_probs[this_base_index] > snp_probs[best_base_index])))) {
        second_best_base_index = best_base_index;
        second_best_base_coverage = best_base_coverage;
        best_base_index = this_base_index;
        best_base_coverage = this_base_coverage;
      }
      else if ((this_base_coverage > second_best_base_coverage) 
        || ((this_base_coverage == second_best_base_coverage) && (((second_best_base_index==base_list_N_index) ||snp_probs[this_base_index] > snp_probs[second_best_base_index])))) {
        second_best_base_index = this_base_index;
        second_best_base_coverage = this_base_coverage;
      }
    }
    
    int this_base_coverage = min(pos_info[best_base_char][0], pos_info[best_base_char][2]);
    
    // Only try mixed SNP model if there is coverage for more than one base!
    ppred.frequency = 1;
    if (second_best_base_index != base_list_N_index) {
      best_base_char = base_char_list[best_base_index];
      second_best_base_char = base_char_list[second_best_base_index];

      // tries all frequencies of the best two
      if (_settings.polymorphism_prediction) {
        ppred = predict_polymorphism(best_base_char, second_best_base_char, pdata);
      }
      
      // tries only the raw ML frequency of the best two
      else /* if (_settings.mixed_base_prediction) */ {
        ppred = predict_mixed_base(best_base_char, second_best_base_char, pdata);
      }
      
      // Calculate E-value for polymorphism score (E theta)
      polymorphism_bonferroni_score = (-(log(ppred.likelihood_ratio_test_p_value)/log(10)) - _log10_ref_length);
      
      // Do we accept this as a polymorphism?
      if (polymorphism_bonferroni_score >= _polymorphism_score_cutoff)
        passed_as_polymorphism_prediction = true;
      
      //cerr << ppred.frequency << " " << ppred.log10_base_likelihood << " " << ppred.p_value << endl;
    }
		
		//###
		//## UNKNOWN UNKNOWN UNKNOWN
		//###
		if(insert_count == 0) {
			update_unknown_intervals(position, p.target(), base_predicted, this_position_unique_only_coverage);
		}
    
    //###
		//## Does any RA evidence pass tests?
		//###
				
    // This is a pointer to the added mutation, so that we can change it later for user-RA
    diff_entry_ptr_t added_mut_p(NULL);
    if (passed_as_consensus_prediction || passed_as_polymorphism_prediction) {
      
      //###
      //## Create new RA evidence mutation for genome diff
      //###		
      
      //## Fields common to consensus mutations and polymorphisms
      mut[SEQ_ID] = p.target_name();
      mut[POSITION] = to_string<uint32_t>(position);
      mut[INSERT_POSITION] = to_string<uint32_t>(insert_count);
      // Genotype quality is for the top called genotype
      mut[CONSENSUS_SCORE] = formatted_double(consensus_bonferroni_score, kMutationScorePrecision).to_string();
      mut[POLYMORPHISM_SCORE] = formatted_double(polymorphism_bonferroni_score, kMutationScorePrecision).to_string();
      
      //## Specific initializations for polymorphisms. Must take precedence.

      //# the frequency returned is the probability of the FIRST base
      //# we want to quote the probability of the second base (the minor allele from the reference).
      mut[REF_BASE] = ref_base_char;
      mut[NEW_BASE] = (best_base_char == ref_base_char) ? second_best_base_char : best_base_char;
      
      mut[MAJOR_BASE] = best_base_char;
      mut[MINOR_BASE] = second_best_base_char;
      mut[MAJOR_FREQUENCY] = formatted_double(ppred.frequency, _polymorphism_precision_places, true).to_string();
      
      double variant_frequency = ppred.frequency;
      if (mut[REF_BASE] == mut[MAJOR_BASE]) {
        variant_frequency = 1.0 - variant_frequency;
      }
      mut[POLYMORPHISM_FREQUENCY] = formatted_double(variant_frequency, _polymorphism_precision_places, true).to_string();

      
      // Add line to R input file if we are only a polymorphism
      
      // @JEB TO_DO: deprecate going to R here
      // and in the main PIPELINE. We can now do
      // Fisher's exact test in C++ and the KS
      // test is probably not necessary
      
      if (ppred.frequency != 1) {
        mut[POLYMORPHISM_EXISTS] = "1";
        write_polymorphism_input_file_line(p, insert_count, ref_base_char, best_base_char, second_best_base_char, ppred, pos_info, pdata );
      }

      //## More fields common to consensus mutations and polymorphisms
      //## ...now that ref_base and new_base are defined
      
      vector<uint32_t>& ref_cov = pos_info[from_string<base_char>(mut[REF_BASE])];
      mut[REF_COV] = to_string(make_pair(static_cast<int32_t>(ref_cov[2]), static_cast<int32_t>(ref_cov[0])));
      
      vector<uint32_t>& new_cov = pos_info[from_string<base_char>(mut[NEW_BASE])];
      mut[NEW_COV] = to_string(make_pair(static_cast<int32_t>(new_cov[2]), static_cast<int32_t>(new_cov[0])));
      
      vector<uint32_t>& major_cov = pos_info[from_string<base_char>(mut[MAJOR_BASE])];
      mut[MAJOR_COV] = to_string(make_pair(static_cast<int32_t>(major_cov[2]), static_cast<int32_t>(major_cov[0])));
      
      vector<uint32_t>& minor_cov = pos_info[from_string<base_char>(mut[MINOR_BASE])];
      mut[MINOR_COV] = to_string(make_pair(static_cast<int32_t>(minor_cov[2]), static_cast<int32_t>(minor_cov[0])));
      
      mut[TOTAL_COV] = to_string(make_pair(total_cov[2], total_cov[0]));
      
      _gd.add(mut);
      
      // what's going on here? we may need to change a value latter,
      // and add added a copy not the current one
      added_mut_p = _gd.get_list().back();
    } // END ra_output
    
    // Now we print additional RA items as user= if they have not already been printed.
    // for this position and insert_count
    while (_user_evidence_ra_list.size()
           && (from_string<uint32_t>((*(_user_evidence_ra_list.front()))[POSITION]) == position)
           && (from_string<int32_t>((*(_user_evidence_ra_list.front()))[INSERT_POSITION]) == insert_count)) {
      
      // Note! mut only exists as a valid record for comparison
      // if output_ra_on_merits is true!!
      if ((passed_as_consensus_prediction || passed_as_polymorphism_prediction) && (*added_mut_p == *(_user_evidence_ra_list.front()))) {
        // in this case we have already created a valid junction
        (*added_mut_p)["user_defined"] = "1";
        
        // do not mark it as user_defined or there will be problems:
        // (it will not show up in HTML and it will mess up polymorphism stats)
      } else {
        // copy the input entry and fill in more information
        cDiffEntry mut = *(_user_evidence_ra_list.front().get());
        mut.to_spec(); //remove additional fields that might be left over!!
        mut["user_defined"] = "1";
        mut["user_defined_no_poly"] = "1"; // prevents doing stats on it
        
        // These are already assigned correctly by copy of RA
        //mut[SEQ_ID] = p.target_name();
        //mut[POSITION] = to_string<uint32_t>(position);
        //mut[INSERT_POSITION] = to_string<uint32_t>(insert_count);

        // tries all frequencies of the best two
        
        base_char best_base_char = from_string<base_char>(mut[REF_BASE]);
        base_char second_best_base_char = from_string<base_char>(mut[NEW_BASE]);
        
        if (_settings.polymorphism_prediction) {
          ppred = predict_polymorphism(best_base_char, second_best_base_char, pdata);
        }
        // tries only the raw ML frequency of the best two
        else {
          ppred = predict_mixed_base(best_base_char, second_best_base_char, pdata);
        }
        
        double polymorphism_bonferroni_score = 10 * (-(log(ppred.likelihood_ratio_test_p_value)/log(10)) - _log10_ref_length);
        
        // These defs of major and minor base are not quite always accurate.
        // because they only take from ref_base and new_base of the RA line
        // ---> in the future it may be better to KEEP the major and minor alleles and work with
        //      these instead, but that leads to its own issues
        mut[MAJOR_BASE] = (ppred.frequency > 0.5) ? best_base_char : second_best_base_char;
        mut[MINOR_BASE] = (ppred.frequency > 0.5) ? second_best_base_char : best_base_char;
        mut[MAJOR_FREQUENCY] = formatted_double( ((ppred.frequency > 0.5) ? ppred.frequency : 1 - ppred.frequency), _polymorphism_precision_places, true).to_string();
        
        // The frequency of the variant is always this, due to the way the ref_base is set as best_base_char
        mut[POLYMORPHISM_FREQUENCY] = formatted_double(1 - ppred.frequency, _polymorphism_precision_places, true).to_string();
        
        // Genotype quality is for the top called genotype
        mut[CONSENSUS_SCORE] = formatted_double(consensus_bonferroni_score, kMutationScorePrecision).to_string();
        mut[POLYMORPHISM_SCORE] = formatted_double(polymorphism_bonferroni_score, kMutationScorePrecision).to_string();
        
        vector<uint32_t>& ref_cov = pos_info[from_string<base_char>(mut[REF_BASE])];
        mut[REF_COV] = to_string(make_pair(static_cast<int32_t>(ref_cov[2]), static_cast<int32_t>(ref_cov[0])));
        
        vector<uint32_t>& new_cov = pos_info[from_string<base_char>(mut[NEW_BASE])];
        mut[NEW_COV] = to_string(make_pair(static_cast<int32_t>(new_cov[2]), static_cast<int32_t>(new_cov[0])));
        
        vector<uint32_t>& major_cov = pos_info[from_string<base_char>(mut[MAJOR_BASE])];
        mut[MAJOR_COV] = to_string(make_pair(static_cast<int>(major_cov[2]), static_cast<int>(major_cov[0])));
        
        vector<uint32_t>& minor_cov = pos_info[from_string<base_char>(mut[MINOR_BASE])];
        mut[MINOR_COV] = to_string(make_pair(static_cast<int>(minor_cov[2]), static_cast<int>(minor_cov[0])));
        
        mut[TOTAL_COV] = to_string(make_pair(total_cov[2], total_cov[0]));
        
        _gd.add(mut);
      }
      
      _user_evidence_ra_list.pop_front();
    }
    
	}
}

/*! Called at the beginning of a reference sequence fragment
    Open per-reference files
*/

void identify_mutations_pileup::at_target_start(const uint32_t tid)
{
    
  // Open per-reference coverage file:
	if(_print_coverage_data) {
		string filename = _settings.file_name(_settings.complete_coverage_text_file_name, "@", target_name(tid));
		_coverage_data.open(filename.c_str());
    ASSERT(!_coverage_data.fail(), "Could not open output file:" + filename);
		_coverage_data << "unique_top_cov" << "\t" << "unique_bot_cov" << "\t" << "redundant_top_cov" << "\t" << "redundant_bot_cov" << "\t" << "raw_redundant_top_cov" << "\t" << "raw_redundant_bot_cov" << "\t" << "e_value" << "\t" << "position" << endl;
	}	
  
  // Polymorphism file used as input to R
  // Only one file for all reference sequences
  if((_settings.polymorphism_prediction || _settings.mixed_base_prediction) && !_polymorphism_r_input_file.is_open()) {
		string filename = _settings.polymorphism_statistics_input_file_name;
		_polymorphism_r_input_file.open(filename.c_str());
    ASSERT(!_polymorphism_r_input_file.fail(), "Could not open output file:" + filename);
		_polymorphism_r_input_file 
    << "seq_id" << "\t"
    << "position" << "\t" 
    << "insert_position" << "\t"
    << "reference_base" << "\t"
    << "major_base" << "\t"
    << "minor_base" << "\t"
    << "major_frequency" << "\t"
    << "major_top_strand" << "\t"
    << "major_bot_strand" << "\t"
    << "minor_top_strand" << "\t"
    << "minor_bot_strand" << "\t"
    << "major_quals" << "\t"
    << "minor_quals" << "\t"
    << endl;
  }
  
  // Reset the Missing Coverage evidence variables
  _last_deletion_start_position = UNDEFINED_UINT32;
	_last_deletion_end_position = UNDEFINED_UINT32;
	_last_deletion_redundant_start_position = UNDEFINED_UINT32;
	_last_deletion_redundant_end_position = UNDEFINED_UINT32;
	_last_start_unknown_interval = UNDEFINED_UINT32;
  
}
  
/*! Called at the end of a reference sequence fragment
    Close per-reference files
 */
void identify_mutations_pileup::at_target_end(const uint32_t tid) {

  // end "open" Missing Coverahge and Unknown intervals
  if (!_settings.no_deletion_prediction) {
    check_deletion_completion(target_length(tid)+1, tid, position_coverage(numeric_limits<double>::quiet_NaN()), numeric_limits<double>::quiet_NaN());
  }
  update_unknown_intervals(target_length(tid)+1, tid, true, false);

  // if this target failed to have its coverage fit, mark the entire thing as a deletion
  double _this_deletion_propagation_cutoff = _deletion_propagation_cutoffs[tid];
  // if the propagation cutoff is zero then the coverage distribution failed
  if (!_settings.no_deletion_prediction && (_this_deletion_propagation_cutoff < 0.0))
  {
    cDiffEntry del(MC);
    del[SEQ_ID] = target_name(tid);
    
    del[START] = to_string<uint32_t>(1);
    del[END] = to_string<uint32_t>(target_length(tid));
    del[START_RANGE] = to_string<uint32_t>(0);
    del[END_RANGE] = to_string<uint32_t>(0);
    
    del[LEFT_OUTSIDE_COV] = "NA";
    del[LEFT_INSIDE_COV] = formatted_double(0.0, 0).to_string();
    del[RIGHT_INSIDE_COV] = formatted_double(0.0, 0).to_string();
    del[RIGHT_OUTSIDE_COV] = "NA";

    _gd.add(del);
  }
  
  // Close per-reference coverage file:
	if(_print_coverage_data)
		_coverage_data.close();
}


/*! Helper method to track information about putative deleted regions.
 
 Used at each pileup iteration and at the end.
 //## when called at the end of a fragment, the position is fragment_length+1
 //## and this_position_coverage is undefined
 
 @JEB This function expects 1-indexed positions!!!
 
 */
void identify_mutations_pileup::check_deletion_completion(uint32_t position, uint32_t seq_id, const position_coverage& this_position_coverage, double e_value_call) {

	//cerr << position << " " << e_value_call << endl;
	
  // special case = beginning of new seq_id
  if (position == 1) 
    _last_position_coverage = position_coverage(numeric_limits<double>::quiet_NaN());
  
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
    if(_last_deletion_start_position == UNDEFINED_UINT32) {
      _last_deletion_start_position = position;
      _left_outside_coverage_item = _last_position_coverage;
      _left_inside_coverage_item = this_position_coverage;
    }
  }
		
  //##keep track of whether we've encountered the seed value
  if(!isnan(this_position_coverage.unique[1]) && (this_position_coverage.total <= _this_deletion_seed_cutoff)) {
    _this_deletion_reaches_seed_value = true;
  }
	
	//## If we are in a deletion and rise back above the propagation cutoff OR we are at the end of this fragment (NAN),
  //## then record the current deletion.
	if( (_last_deletion_start_position != UNDEFINED_UINT32) && 
     ( isnan(this_position_coverage.unique[1]) || (this_position_coverage.unique[1] > _this_deletion_propagation_cutoff) ) )
  {
		
		if(_this_deletion_reaches_seed_value) {

      _last_deletion_end_position = position-1;
      if (_last_deletion_redundant_end_position == UNDEFINED_UINT32) 
        _last_deletion_redundant_end_position = _last_deletion_end_position;
      if (_last_deletion_redundant_start_position == UNDEFINED_UINT32) 
        _last_deletion_redundant_start_position = _last_deletion_start_position;

      cDiffEntry del(MC);
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
		_last_deletion_start_position = UNDEFINED_UINT32;
		_last_deletion_end_position = UNDEFINED_UINT32;
		_last_deletion_redundant_start_position = UNDEFINED_UINT32;
		_last_deletion_redundant_end_position = UNDEFINED_UINT32;
	}
  
  //## REDUNDANT COVERAGE
  //## updated only if we are still within a deletion
  if (_last_deletion_start_position != UNDEFINED_UINT32) {
    
    if (this_position_coverage.redundant[1] == 0) {
      _this_deletion_redundant_reached_zero = true;
      _last_deletion_redundant_end_position = UNDEFINED_UINT32;
    }
    else if (this_position_coverage.redundant[1] > 0) {
      //## if there is any redundant coverage remember the start (until we find zero redundant coverage)
      if (!_this_deletion_redundant_reached_zero) {
        _last_deletion_redundant_start_position = position;
      }
      else if (_last_deletion_redundant_end_position == UNDEFINED_UINT32) {
        _last_deletion_redundant_end_position = position;
      }
    }
  }
	
  
  _last_position_coverage = this_position_coverage;
}


/*! Helper method to track unknowns.
 */
void identify_mutations_pileup::update_unknown_intervals(uint32_t position, uint32_t seq_id, bool base_predicted, bool this_position_unique_only_coverage) 
{
  //debug
  /*
	cerr << position << " " << base_predicted << " " << this_position_unique_only_coverage << endl;
	if(_last_start_unknown_interval != UNDEFINED_UINT32) {
		cerr << _last_start_unknown_interval << endl;
	} else {
		cerr << "undef" << endl;
	}
  */
	
	if(!base_predicted) {
		if(this_position_unique_only_coverage) {
			++s.coverage_unique_uncalled;
		}
		if(_last_start_unknown_interval == UNDEFINED_UINT32) {
			_last_start_unknown_interval = position;
		}
	}	else {
		if(this_position_unique_only_coverage) {
			++s.coverage_unique_called;
		}
			
		//#end interval where we were unable to call mutations
		if(_last_start_unknown_interval != UNDEFINED_UINT32) {
      cDiffEntry new_interval(UN);
			new_interval[SEQ_ID] = target_name(seq_id);
			new_interval[START] = to_string<uint32_t>(_last_start_unknown_interval);
			new_interval[END] = to_string<uint32_t>(position - 1);
			_gd.add(new_interval);
			
			_last_start_unknown_interval = UNDEFINED_UINT32;
		}
	}
}
  
void identify_mutations_pileup::write_polymorphism_input_file_line(const pileup& p, const uint32_t insert_count, char ref_base_char, char best_base_char, char second_best_base_char, const polymorphism_prediction& ppred, position_base_info& pos_info, const vector<polymorphism_data>& pdata )
{
  
  _polymorphism_r_input_file
  << p.target_name() << "\t"
  << p.position_1() << "\t"
  << insert_count << "\t"
  << ref_base_char << "\t"
  << best_base_char << "\t"
  << second_best_base_char << "\t"
  << ppred.frequency << "\t"
  << pos_info[best_base_char][2] << "\t"
  << pos_info[best_base_char][0] << "\t"
  << pos_info[second_best_base_char][2] << "\t"
  << pos_info[second_best_base_char][0] << "\t"
  ;
  
  
  string best_base_qualities;
  string second_best_base_qualities;
  
  for(vector<polymorphism_data>::const_iterator it=pdata.begin(); it<pdata.end(); ++it) {
    
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
  
  _polymorphism_r_input_file << best_base_qualities << "\t";
  _polymorphism_r_input_file << second_best_base_qualities << "\t";
  
  _polymorphism_r_input_file << endl;
  
}


/*! Predict the significance of putative polymorphisms.
 */
polymorphism_prediction identify_mutations_pileup::predict_polymorphism(base_char best_base_char, base_char second_best_base_char, vector<polymorphism_data>& pdata ) {
  
  //#calculate the likelihood of observed reads given this position is 100% the best base
	double log10_likelihood_of_one_base_model = 0;
  for(vector<polymorphism_data>::iterator it=pdata.begin(); it<pdata.end(); ++it) {
  
    double log10_correct_pr;

    covariate_values_t this_cv = it->_cv;
      
    if(it->_strand == 1) {
      this_cv.ref_base() = basechar2index(best_base_char);
    } else {
      this_cv.ref_base() = basechar2index(complement_base_char(best_base_char)); 
      this_cv.obs_base() = complement_base_index(this_cv.obs_base());
    }
    log10_correct_pr = _error_table.get_log10_prob(this_cv);
       
    log10_likelihood_of_one_base_model  += log10_correct_pr;
  }

	vector<uint8_t> best_base_qualities;
	vector<uint8_t> second_best_base_qualities;
	uint32_t best_base_strand_hash[] = {0, 0};
	uint32_t second_best_base_strand_hash[] = {0, 0};
  
  for(vector<polymorphism_data>::iterator it=pdata.begin(); it<pdata.end(); ++it) {
  
    int8_t zp_strand = (it->_strand == +1) ? 1 : 0;
    if (it->_base_char == best_base_char) {
      best_base_qualities.push_back(it->_quality);
      best_base_strand_hash[zp_strand]++;
    }
    else if (it->_base_char == second_best_base_char) {
      second_best_base_qualities.push_back(it->_quality);
      second_best_base_strand_hash[zp_strand]++;
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
  
/*! Predict the significance of putative polymorphisms.
 */
polymorphism_prediction identify_mutations_pileup::predict_mixed_base(base_char best_base_char, base_char second_best_base_char, vector<polymorphism_data>& pdata ) {
  
  //#calculate the likelihood of observed reads given this position is 100% the best base  
  double log10_likelihood_of_one_base_model = 0;
  for(vector<polymorphism_data>::iterator it=pdata.begin(); it<pdata.end(); ++it) {
    
    double log10_correct_pr;
    
    covariate_values_t this_cv = it->_cv;
    
    if(it->_strand == 1) {
      this_cv.ref_base() = basechar2index(best_base_char);
    } else {
      this_cv.ref_base() = basechar2index(complement_base_char(best_base_char)); 
      this_cv.obs_base() = complement_base_index(this_cv.obs_base());
    }
    log10_correct_pr = _error_table.get_log10_prob(this_cv);
    
    log10_likelihood_of_one_base_model  += log10_correct_pr;
  }
  
  vector<uint8_t> best_base_qualities;
  vector<uint8_t> second_best_base_qualities;
  uint32_t best_base_strand_hash[] = {0, 0};
  uint32_t second_best_base_strand_hash[] = {0, 0};
  
  for(vector<polymorphism_data>::iterator it=pdata.begin(); it<pdata.end(); ++it) {
    
    int8_t zp_strand = (it->_strand == +1) ? 1 : 0;

    if (it->_base_char == best_base_char) {
      best_base_qualities.push_back(it->_quality);
      best_base_strand_hash[zp_strand]++;
    }
    else if (it->_base_char == second_best_base_char) {
      second_best_base_qualities.push_back(it->_quality);
      second_best_base_strand_hash[zp_strand]++;
    }
  }
  
  
  
  // Unlike full polymorphism prediction, we test just the raw frequency of the two bases
  // and do not check for bias later.
  
  // This would calculate the true best frequency -- slower, but more accurate
  // pair<double,double> best_two_base_model = best_two_base_model_log10_likelihood(best_base_char, second_best_base_char, pdata);
  
  double max_likelihood_fr_first_base = static_cast<double>(best_base_strand_hash[0] + best_base_strand_hash[1]) 
    / static_cast<double>(best_base_strand_hash[0] + best_base_strand_hash[1] + second_best_base_strand_hash[0] + second_best_base_strand_hash[1]);
  
  double log10_likelihood_of_two_base_model = calculate_two_base_model_log10_likelihood(
                                                                                        best_base_char, 
                                                                                        second_best_base_char, 
                                                                                        pdata, 
                                                                                        max_likelihood_fr_first_base
                                                                                        );
  
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
    
    p_value = pchisq(1.0L, likelihood_ratio_test_value);
    //cerr << "likelihood_ratio_test_value: " << likelihood_ratio_test_value << " p-value: " << p_value << endl;
  }
  
  polymorphism_prediction p(max_likelihood_fr_first_base, log10_likelihood_of_one_base_model - log10_likelihood_of_two_base_model, p_value);
  
  return p;
}
 
  
double identify_mutations_pileup::slope_at_percentage_best_base(base_char best_base_char, base_char second_best_base_char, vector<polymorphism_data>& pdata, const double guess, const double precision, double& middle_point_log10_likelihood)  
{
  // precision is a fraction of the value
  double point_1 = max(guess * (1 - precision), 0.0);  
  double point_2 = min(guess * (1 + precision), 1.0);
  
  double point_1_log10_likelihood = calculate_two_base_model_log10_likelihood(best_base_char, second_best_base_char, pdata, point_1);
  double point_2_log10_likelihood = calculate_two_base_model_log10_likelihood(best_base_char, second_best_base_char, pdata, point_2);
  
  middle_point_log10_likelihood = calculate_two_base_model_log10_likelihood(best_base_char, second_best_base_char, pdata, guess);

  // Check for local minimum!
  if ( (middle_point_log10_likelihood < point_1_log10_likelihood) && (middle_point_log10_likelihood < point_2_log10_likelihood) ) {
    return 0;
  }
  
  return (point_2_log10_likelihood - point_1_log10_likelihood) / (point_2 - point_1);
}

/*! Find the best fraction for the best base at a polymorphic site.
 */
  
pair<double,double> identify_mutations_pileup::best_two_base_model_log10_likelihood(base_char best_base_char, base_char second_best_base_char, vector<polymorphism_data>& pdata)
{	
  uint32_t iterations = 0;
  double current_upper_pr_first_base = 1.0;
  double current_lower_pr_first_base = 0.0;
  
  double current_upper_pr_first_base_log10_likelihood = calculate_two_base_model_log10_likelihood(best_base_char, second_best_base_char, pdata, current_upper_pr_first_base);
  double current_lower_pr_first_base_log10_likelihood = calculate_two_base_model_log10_likelihood(best_base_char, second_best_base_char, pdata, current_lower_pr_first_base);
  
  // precision is a fraction of the value...
  while (current_upper_pr_first_base - current_lower_pr_first_base > (current_upper_pr_first_base + current_lower_pr_first_base) / 2 * _polymorphism_precision_decimal) {
    
    iterations++;
    //cout << iterations << " "  << current_lower_pr_first_base << " " << current_upper_pr_first_base << endl;
    
    double current_middle_pr_first_base = (current_upper_pr_first_base + current_lower_pr_first_base) / 2;
    
    double current_middle_pr_first_base_log10_likelihood;
    double middle_slope = slope_at_percentage_best_base(best_base_char, second_best_base_char, pdata, current_middle_pr_first_base, _polymorphism_precision_decimal, current_middle_pr_first_base_log10_likelihood);
    
    // Slope is set to zero if the tested point is better than the ones 
    // on either side, when calculating the slope.
    if ( middle_slope == 0) {
      return make_pair(current_middle_pr_first_base, current_middle_pr_first_base_log10_likelihood);
    
    } else if ( middle_slope < 0) {
      
      current_upper_pr_first_base = current_middle_pr_first_base;
      current_upper_pr_first_base_log10_likelihood = current_middle_pr_first_base_log10_likelihood;
      
    } else {
      
      current_lower_pr_first_base = current_middle_pr_first_base;
      current_lower_pr_first_base_log10_likelihood = current_middle_pr_first_base_log10_likelihood;
    }
  }  
    
  if (current_lower_pr_first_base_log10_likelihood > current_upper_pr_first_base_log10_likelihood) {

    return make_pair(current_lower_pr_first_base, current_lower_pr_first_base_log10_likelihood);
  }
  
  return make_pair(current_upper_pr_first_base, current_upper_pr_first_base_log10_likelihood);
}
  
/*
pair<double,double> identify_mutations_pileup::best_two_base_model_log10_likelihood(base_char best_base_char, base_char second_best_base_char, vector<polymorphism_data>& pdata)
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
	}
	
	return make_pair(last_pr_first_base, last_log_pr);
}
 */

/*! Calculate the likelihood of a mixture model of two bases leading to the observed read bases.
 */
double identify_mutations_pileup::calculate_two_base_model_log10_likelihood(
                                                                            base_char best_base_char, 
                                                                            base_char second_best_base_char, 
                                                                            const vector<polymorphism_data>& pdata, 
                                                                            double best_base_freq
                                                                            )
{
	double log10_likelihood = 0;	
	
  for(vector<polymorphism_data>::const_iterator it=pdata.begin(); it<pdata.end(); ++it) {
  
    //## the first value is pr_base, second is pr_not_base
    double best_base_log10pr;
    double second_best_base_log10pr;
    
    covariate_values_t this_cv = it->_cv;
    
    if (it->_strand == -1) {
      this_cv.obs_base() = complement_base_index(this_cv.obs_base());
    }
    
    if(it->_strand == 1) {
      this_cv.ref_base() = basechar2index(best_base_char);
    } else {
      this_cv.ref_base() = basechar2index(complement_base_char(best_base_char));
    }
    best_base_log10pr = _error_table.get_log10_prob(this_cv);
    
    if(it->_strand == 1) {
      this_cv.ref_base() = basechar2index(second_best_base_char);
    } else {
      this_cv.ref_base() = basechar2index(complement_base_char(second_best_base_char));
    }
    second_best_base_log10pr = _error_table.get_log10_prob(this_cv);

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

cDiscreteSNPCaller::cDiscreteSNPCaller(
                                       const string& type,
                                       uint32_t reference_length
                                       ) 
: _type(type)
{
  

  // Uniform priors across all bases.
  if (_type == "haploid") {
    
    double uniform_probability = 1.0 / static_cast<double>(base_list_size);
    add_genotype("A", uniform_probability);
    add_genotype("C", uniform_probability);
    add_genotype("G", uniform_probability);
    add_genotype("T", uniform_probability);
    add_genotype(".", uniform_probability);
  }
  
  /*
  // Prior that we expect one change from reference
  else if (_type == "haploid-change") {
    
    // recall that first one counts as reference
    double uniform_probability = 1.0 / reference_length;    
    add_genotype("A", 1.0 - 4.0 * uniform_probability);
    add_genotype("C", uniform_probability);
    add_genotype("G", uniform_probability);
    add_genotype("T", uniform_probability);
    add_genotype(".", uniform_probability);
  }
  */
  
  // Extra states and priors for unexpected mixed states... experimental
  else if (_type == "haploid-cnv") {
    
    double mixed_probability = 1.0 / reference_length;    
    double uniform_probability = (1.0 - 10 * mixed_probability) / base_list_size;    
    
    add_genotype("A", uniform_probability);
    add_genotype("C", uniform_probability);
    add_genotype("G", uniform_probability);
    add_genotype("T", uniform_probability);
    add_genotype(".", uniform_probability);
    
    add_genotype("AC", mixed_probability);
    add_genotype("AG", mixed_probability);
    add_genotype("AT", mixed_probability);
    add_genotype("A.", mixed_probability);
    
    add_genotype("CG", mixed_probability);
    add_genotype("CT", mixed_probability);
    add_genotype("C.", mixed_probability);
    
    add_genotype("GT", mixed_probability);
    add_genotype("G.", mixed_probability);
    
    add_genotype("T.", mixed_probability);
  }  
  
  else
  {
    ERROR("Unknown SNP Caller type:" + type);
  }
  
  // Check priors
  double total_probability = 0;
  for(size_t i=0; i<_log10_genotype_prior_probabilities.size(); i++) {
    total_probability += pow(10, _log10_genotype_prior_probabilities[i]);
  }
  ostringstream ss;
  ss << setprecision(5) << total_probability;
  
  ASSERT( from_string<double>(ss.str()) == 1.0, "Prior probabilities do not sum to 1. (" + to_string(total_probability) + ").")
  
  reset(0);
}
  
void cDiscreteSNPCaller::add_genotype(const string& genotype, double probability) {
  
  _log10_genotype_prior_probabilities.push_back(log10(probability));
  
  vector<base_index> gv;
  for(size_t i=0; i<genotype.length(); i++) {
    gv.push_back(basechar2index(genotype[i]));
  }
  _genotype_vector.push_back(gv);
  
}
  
  
void cDiscreteSNPCaller::reset(uint8_t ref_base_index) {
  _best_genotype_index = 0;
  _observations = 0;
  _normalized_observations = 0;
  _log10_genotype_probabilities = _log10_genotype_prior_probabilities;
  
  (void) ref_base_index;
  //this is for where there are unbalanced priors -- haploid-change
  //they do not behave properly when the reference is 'N'
  //swap(_genotype_probability[0], _genotype_probability[ref_base_index]);
}
  
void cDiscreteSNPCaller::update(const covariate_values_t& cv, bool obs_top_strand, int32_t mapping_quality, cErrorTable& et) {

  covariate_values_t this_cv = cv;
  //update probabilities give observation using Bayes rule
  
  if (!obs_top_strand) {
    this_cv.obs_base() = complement_base_index(this_cv.obs_base()); 
  }
  
  double incorrect_mapping_prob = pow(10, -static_cast<double>(mapping_quality) / 10);
  double correct_mapping_prob = 1 - incorrect_mapping_prob;
  this->_normalized_observations += correct_mapping_prob;
  
  double total_prob = 0.0;
  for (uint32_t i=0; i<_genotype_vector.size(); i++) {
  
    vector<base_index>& gv = this->_genotype_vector[i];
    double this_pr = 0.0;
    
    for (uint32_t j=0; j < gv.size(); j++) {
      
      this_cv.ref_base() = gv[j];
      
      if (!obs_top_strand) {
        this_cv.ref_base() = complement_base_index(this_cv.ref_base()); 
      }
      this_pr += (correct_mapping_prob * et.get_prob(this_cv) + incorrect_mapping_prob * 1.0 / _genotype_vector.size()) * pow(10, this->_log10_genotype_probabilities[i]) / gv.size();
      
      // Floating point error can make this a very  negative number
      if (this_pr < 0.0) this_pr = 0.0;
    }
    
    total_prob += this_pr;
  }

  double highest_pr = -numeric_limits<double>::max();
  for (uint32_t i=0; i<this->_genotype_vector.size(); i++) {
    
    vector<base_index>& gv = this->_genotype_vector[i];
    double this_pr = 0.0;
    
    for (uint32_t j=0; j < gv.size(); j++) {
      
      this_cv.ref_base() = gv[j];
      
      if (!obs_top_strand) {
        this_cv.ref_base() = complement_base_index(this_cv.ref_base()); 
      }
      this_pr += (correct_mapping_prob * et.get_prob(this_cv) + incorrect_mapping_prob * 1 / _genotype_vector.size()) / gv.size();
    }
    
    this->_log10_genotype_probabilities[i] += log10(this_pr) - log10(total_prob);
    
    if (this->_log10_genotype_probabilities[i] > highest_pr) {
      this->_best_genotype_index = i;
      highest_pr = this->_log10_genotype_probabilities[i];
    }
  }
  
  //print();
  
  _observations++;
}
  
void cDiscreteSNPCaller::print() {
  
  for (uint32_t i=0; i<_genotype_vector.size(); i++) {
    
    vector<base_index>& gv = _genotype_vector[i];
    
    cout << "Genotype: " ;
    for (uint32_t j=0; j < gv.size(); j++) {
      
      cout << baseindex2char(gv[j]);
    }
    
    cout << " " << _log10_genotype_probabilities[i] << endl;
  }
}


cSNPCall cDiscreteSNPCaller::get_prediction()
{
  cSNPCall snp_call;
  
  if (_observations == 0) {
    //Best base is 'N' and E-value is NAN
    return snp_call;
  }
  
  // need to go through and find the most probable
  snp_call.genotype = "";
  for(size_t i=0; i<_genotype_vector[_best_genotype_index].size(); i++)
    snp_call.genotype += baseindex2char(_genotype_vector[_best_genotype_index][i]);
    
  // we want to normalize the probabilities, but avoid floating point errors    
  double log10_offset_probability = -numeric_limits<double>::max();
  for(size_t i=0; i < _log10_genotype_probabilities.size(); i++) {
    if (i != _best_genotype_index)
      log10_offset_probability = max(log10_offset_probability, _log10_genotype_probabilities[i]);
  }
  
  double total_error_probability = 0;
  for (uint32_t i=0; i<_genotype_vector.size()-1; i++) {
    if (i != _best_genotype_index)
      total_error_probability += pow(10, _log10_genotype_probabilities[i] - log10_offset_probability);
  }
  double log10_total_error_probability = log10(total_error_probability);
  log10_total_error_probability += log10_offset_probability;
  
  snp_call.score = (_log10_genotype_probabilities[_best_genotype_index] - log10_total_error_probability);
  
  return snp_call;
}

} // namespace breseq

