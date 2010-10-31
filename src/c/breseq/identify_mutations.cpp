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
#include <vector>
#include <bam.h>
#include <sam.h>
#include <faidx.h>
#include <assert.h>
#include <boost/tuple/tuple.hpp>
#include <boost/algorithm/string.hpp>

#include "breseq/common.h"
#include "breseq/pileup.h"
#include "breseq/identify_mutations.h"

/*! Convenience wrapper around the identify_mutations_pileup class.
 */
void breseq::identify_mutations(const std::string& bam, 
																const std::string& fasta,
																const std::string& error_dir,
																const std::string& gd_file,
																const std::string& output_dir,
																const std::vector<std::string>& readfiles,
																const std::string& coverage_dir,
																const std::vector<double>& deletion_propagation_cutoff,
																double mutation_cutoff,
																bool predict_deletions,
																bool predict_polymorphisms,
                                uint8_t min_qual_score) {
                                                                                            
	// do the mutation identification:
	identify_mutations_pileup imp(bam, fasta, error_dir, gd_file, output_dir, readfiles, coverage_dir, deletion_propagation_cutoff, mutation_cutoff, predict_deletions, predict_polymorphisms, min_qual_score);
	imp.do_pileup();
}


/*! Constructor.
 */
breseq::identify_mutations_pileup::identify_mutations_pileup(const std::string& bam, 
																														 const std::string& fasta,
																														 const std::string& error_dir,
																														 const std::string& gd_file,
																														 const std::string& output_dir,
																														 const std::vector<std::string>& readfiles,
																														 const std::string& coverage_dir,
																														 const std::vector<double>& deletion_propagation_cutoff,
																														 double mutation_cutoff,
																														 bool predict_deletions,
																														 bool predict_polymorphisms,
                                                             uint8_t min_qual_score)
: breseq::pileup_base(bam, fasta)
, _ecr(error_dir, readfiles)
, _gd(gd_file)
, _min_qual_score(min_qual_score)
, _deletion_seed_cutoff(0)
, _deletion_propagation_cutoff(deletion_propagation_cutoff)
, _mutation_cutoff(mutation_cutoff)
, _predict_deletions(predict_deletions)
, _predict_polymorphisms(predict_polymorphisms)
, _coverage_dir(coverage_dir)
, _log10_ref_length(0)
, _on_deletion_seq_id(boost::none)
, _this_deletion_reaches_seed_value(false)
, _last_position_coverage_printed(0) {
	  
  assert(_bam->header->n_targets == (int32_t)_deletion_propagation_cutoff.size());
    
	// reserve enough space for the sequence info:
	_seq_info.resize(_bam->header->n_targets);
	
	//	my $total_ref_length = 0;
	//	foreach my $seq_id (@seq_ids)
	//	{
	//		$total_ref_length+= $bam->length($seq_id);
	//	}
	//	
	//	my $log10_ref_length = log($total_ref_length) / log(10);	
	for(int i=0; i<_bam->header->n_targets; ++i) {
		_log10_ref_length += static_cast<double>(_bam->header->target_len[i]);
	}
	assert(_log10_ref_length != 0);
	_log10_ref_length = log10(_log10_ref_length);
	
	//	my $snps_all_tab_file_name = $settings->file_name('complete_mutations_text_file_name', {'@'=>$seq_id}); 
	//	my $coverage_tab_file_name = $settings->file_name('complete_coverage_text_file_name', {'@'=>$seq_id}); 
	//	
	//	open MUT, ">$snps_all_tab_file_name" if (defined $snps_all_tab_file_name);
	//	open COV, ">$coverage_tab_file_name" if (defined $coverage_tab_file_name);
	//	
	//	my $cnv_coverage_tab_file_name = $settings->file_name('cnv_coverage_tab_file_name', {'@'=>$seq_id}); 
	//	open CNV_COV, ">$cnv_coverage_tab_file_name" if (defined $cnv_coverage_tab_file_name);
}


/*! Destructor.
 */
breseq::identify_mutations_pileup::~identify_mutations_pileup() {
}


/*! Called for each alignment.
 */
void breseq::identify_mutations_pileup::callback(const breseq::pileup& p) {
	using namespace std;
	assert(p.target() < _seq_info.size());
  _this_deletion_propagation_cutoff = _deletion_propagation_cutoff[p.target()];

	//	our @base_list = ('A', 'T', 'C', 'G', '.');
	static uint8_t base_list[] = {'A', 'T', 'C', 'G', '.'};
	
	//	my ($seqid,$pos,$pileup) = @_;
	//	
	//	print STDERR "    POSITION:$pos\n" if ($pos % 10000 == 0);			
	//	
	//	my $insert_count = 0;
	//	my $next_insert_count_exists = 1;
  
  // @JEB Use 1-indexing throughout!!
  uint32_t position = p.position()+1;
	int insert_count=-1;
	bool next_insert_count_exists=true;
	
	// check to see if we already opened the coverage file:
	if(!_coverage_data.is_open()) {
		//		open COV, ">$coverage_tab_file_name" if (defined $coverage_tab_file_name);
		//		print COV join("\t", 'unique_top_cov', 'unique_bot_cov', 'redundant_top_cov', 'redundant_bot_cov', 'raw_redundant_top_cov', 'raw_redundant_bot_cov', 'e_value', 'position') . "\n";
		std::string filename(_coverage_dir);
		filename += p.target_name();
		filename += ".coverage.tab";
		_coverage_data.open(filename.c_str());
		_coverage_data << "unique_top_cov" << "\t" << "unique_bot_cov" << "\t" << "redundant_top_cov" << "\t" << "redundant_bot_cov" << "\t" << "raw_redundant_top_cov" << "\t" << "raw_redundant_bot_cov" << "\t" << "e_value" << "\t" << "position" << std::endl;
	}	
	
	//INSERT_COUNT: while ($next_insert_count_exists)
	while(next_insert_count_exists) {
		++insert_count; // we're doing this here because the perl version uses a while-continue.
    //cerr << "Position:" << position << " Insert Count:" << insert_count << endl;
		//		@dk: positions are 0-indexed here, while genome diff is 1-indexed; this will have to be fixed (maybe in the genome diff?).
		//		if(insert_count) {
		//			cerr << position << " " << insert_count << endl;
		//		}
		
		//	$next_insert_count_exists = 0;
		next_insert_count_exists = false;
		
		//	my $ref_base = ($insert_count) ? '.' : $bam->segment($seqid,$pos,$pos)->dna;
		uint8_t ref_base = '.';
		if(!insert_count) {
			ref_base = p.reference_sequence()[position-1]; //reference_sequence is 1-indexed
		}
		
		//# zero out the info about this position
		//	my $pos_info;
		//	foreach my $base (@base_list)
		//	{
		//		$pos_info->{$base}->{unique_cov}->{1} = 0;
		//		$pos_info->{$base}->{unique_cov}->{-1} = 0;
		//		$pos_info->{$base}->{unique_trimmed_cov}->{1} = 0;
		//		$pos_info->{$base}->{unique_trimmed_cov}->{-1} = 0;
		//		$pos_info->{$base}->{mutation_cov}->{1} = 0;
		//		$pos_info->{$base}->{mutation_cov}->{-1} = 0;
		//	}
		map<uint8_t,position_info> pos_info;
		for(std::size_t j=0; j<5; ++j) {
			pos_info.insert(make_pair(base_list[j],position_info()));
		}
		
		//	
		//## keep track of coverage for deletion prediction
		//	my $this_position_coverage;
		//	$this_position_coverage->{unique} = {'-1' => 0, '1' => 0, 'total' => 0};
		//	$this_position_coverage->{redundant} = {'-1' => 0, '1' => 0, 'total' => 0};
		//	$this_position_coverage->{raw_redundant} = {'-1' => 0, '1' => 0, 'total' => 0};
		//	$this_position_coverage->{total} = 0;
		//	my $this_position_unique_only_coverage = 1;
		position_coverage this_position_coverage;
		bool this_position_unique_only_coverage=true;
		
		//## calculate the chance of observing alignment given each possible base was 100% of the population
		//	my $pr_base_hash;
		//	my $pr_not_base_hash;
		map<uint8_t,double> pr_base_hash;
		map<uint8_t,double> pr_not_base_hash;
		
		//## polymorphism prediction
		//	my $pdata;
		vector<polymorphism_data> pdata;
		
		
		///*		
		//	ALIGNMENT: foreach my $p (@$pileup) 
		//		{
		
		// for each alignment within this pileup:
		for(pileup::const_iterator i=p.begin(); i!=p.end(); ++i) {
			//			my $a = $p->alignment;
			//		*i;
			
			//## This setup gives expected behavior from indel!
			//			my $indel = $p->indel;       ## insertions relative to the reference have the 
			//## number of this inserted base
			//			$indel = 0 if ($indel < 0);  ## substitute such that
			//			$indel = -1 if ($p->is_del); ## deletions relative to reference have -1 as indel
			int indel=i->indel();
			if(indel < 0) {
				indel = 0;
			}
			if(i->is_del()) {
				indel = -1;
			}
			
			//my $base = ($indel < $insert_count) ? '.' : substr($a->qseq,$p->qpos + $insert_count,1);		
			uint8_t base='.';
			if(indel >= insert_count) {
        base = i->query_base(i->query_position() + insert_count);
			}
			
			//##don't use bases without qualities!!
			//next if ($base =~ /[nN]/);
			if(is_N(base)) {
				continue;
			}
			
			//			my $redundancy = $a->aux_get('X1');
			int32_t redundancy = i->redundancy();
			
			//			my $fastq_file_index = $a->aux_get('X2');
			int32_t fastq_file_index = i->fastq_file_index();
			
			//			my $strand = $a->reversed ? -1 : +1;
			int strand = i->strand();
			
			//## Handle trimming
			//## Note that trimming INCLUDES the unaligned bases on each end
			//			my $trimmed = 0;
			//			my $trim_left = $a->aux_get('XL');  
			//			my $trim_right = $a->aux_get('XR');
			//
			//			$trimmed = 1 if ((defined $trim_left) && ($p->qpos+1 <= $trim_left));
			//			$trimmed = 1 if ((defined $trim_right) && ($a->query->length-$p->qpos <= $trim_right));
			//			
			bool trimmed = i->is_trimmed();
			
			//std::cerr << q_start << " " << q_end << std::endl;			
			
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
				//## this is only used when reporting coverage for within-read indels
				//## NOT for calling deletions...		
				//$this_position_coverage->{unique}->{$strand}++;
				//$pos_info->{$base}->{unique_cov}->{$strand}++;
				++this_position_coverage.unique[1+strand];
				++pos_info[base2char(base)].unique_cov[1+strand];
				
				//if ($indel > $insert_count)
				//{
				//$next_insert_count_exists = 1;
				//}
				if(indel > insert_count) {
					next_insert_count_exists = true;
				}
				
			} else {
				//$this_position_unique_only_coverage = 0;
				//$this_position_coverage->{redundant}->{$strand} += 1/$redundancy;			
				//$this_position_coverage->{raw_redundant}->{$strand}++;			
				this_position_unique_only_coverage = false;
				this_position_coverage.redundant[1+strand] += 1.0/redundancy;
				++this_position_coverage.raw_redundant[1+strand];
			}
			
			
			//## EXPERIMENTAL -- moved above		
			//##don't use information from trimmed reads!!
			//			next if ($trimmed);
			//			
			//##don't use information from redundant reads!!
			//next if ($redundancy > 1);				
			if(trimmed || (redundancy > 1)) {
				continue;
			}
			
			
			//## These are the start and end coordinates of the aligned part of the read
			//			my ($q_start, $q_end) = Breseq::Shared::alignment_query_start_end($a, {no_reverse=>1});
			int32_t q_start,q_end;
			boost::tie(q_start,q_end) = i->query_bounds(); // @dk: 1-indexed!
			
			uint8_t quality=0;
      
      //## Deletion in read relative to reference...
      //## Quality is of the NEXT base in the read, and check that it is not an N
      //## Note: This is for a deletion when $insert_count == 0 and against an insertion when $insert_count > 0

      if (indel == -1)
      {			
        int32_t mqpos = i->query_position() + 1 - i->reversed(); 
        uint8_t check_base = i->query_base(mqpos);
        if (is_N(check_base)) continue;
        quality = i->quality_base(mqpos);
        //	my $mqpos = $qpos + 1 - $reversed;
        //	my $check_base = substr($a->qseq,$mqpos,1);
        //	next ALIGNMENT if ($check_base eq 'N');
        //	$quality = $a->qscore->[$mqpos];
      }

      //## Substitution in read relative to reference...
      //## Quality is of the current base in the read, we have ALREADY checked that it is not an N					
      else if (insert_count == 0)
      {
        quality = i->quality_base(i->query_position());
        //	$quality = $a->qscore->[$qpos];
      }
      
      //## Insertion in read relative to reference...
      //## Quality is of the NEXT base in the read, and check that it is not an N
      //## Note that it is possible this read base may be a '.' (supporting the non-insert call)
      else //if (insert_count > 0) 
      {		
        int32_t max_offset = insert_count;
        if (indel < max_offset) max_offset = indel;
        int32_t mqpos = i->query_position() + max_offset + 1 - i->reversed(); 
      
        //my $max_offset = $insert_count;
        //$max_offset = $indel if ($indel < $max_offset);
        //my $mqpos = $qpos + $max_offset + 1 - $reversed;
          
        //## Check bounds: it's possible to go past the end of the read because
        //## this is the last base of this read, but other reads have inserted bases
      
        if (mqpos >= q_end) continue;  // @JEB unlike Perl, this is comparing 0-indexed to 1-indexed
        //next ALIGNMENT if ($mqpos > $q_end);
      
        uint8_t check_base = i->query_base(mqpos);
        if (is_N(check_base)) continue;
        //my $check_base = substr($a->qseq,$mqpos,1);
        //next ALIGNMENT if ($check_base eq 'N');
        
        quality = i->quality_base(mqpos);
        //$quality = $a->qscore->[$mqpos];
      }

      //## We may want to ignore all bases below a certain quality when calling mutations and polymorphisms
      //## This is the check for whether the base fails; it should be after coverage counting
      //next ALIGNMENT if ( $settings->{base_quality_cutoff} && ($quality < $settings->{base_quality_cutoff}) );

      if (quality < _min_qual_score) {
//        std::cerr << position << " " << (unsigned int)quality << " " << std::endl;
        continue;
      }
			
			
			//## this is the coverage for SNP counts, tabulate AFTER skipping trimmed reads
			//$pos_info->{$base}->{unique_trimmed_cov}->{$strand}++;
			++pos_info[base2char(base)].unique_trimmed_cov[1+strand];
			
			//##### this is for polymorphism prediction and making strings
			//push @$pdata, { base => $base, quality => $quality, strand => $strand, fastq_file_index => $fastq_file_index };
			pdata.push_back(polymorphism_data(base,quality,strand,fastq_file_index));
			
			//##### deal with base calls
			//foreach my $hypothetical_base (@base_list)
			//{				
			
			//std::cerr << "========" << std::endl << position+1 << std::endl;
			
			for(std::size_t j=0; j<5; ++j) {
				// base_list[j] == hypothetical base
				//my $base_key =  ($strand == +1) ? $hypothetical_base . $base : Breseq::Fastq::revcom($hypothetical_base) . Breseq::Fastq::revcom($base);
				string base_key;
				if(strand == 1) {
					base_key += base_list[j];
					base_key += base2char(base);
				} else {
					base_key += reverse_base(base_list[j]); 
					base_key += base2char(reverse_base(base));
				}
				
				//##sanity checks
				//
				//##is the error rate defined?
				//
				//if (!defined $log10_correct_rates->[$fastq_file_index]->{$quality}->{$base_key})
				//{
				//print "$fastq_file_index $quality $base_key\n";
				//print Dumper($log10_correct_rates->[$fastq_file_index]);
				//die;
				//}
				
				//##record evidence for and against this hypothetical base being the reference, given the observation
				//$pr_base_hash->{$hypothetical_base} += $log10_correct_rates->[$fastq_file_index]->{$quality}->{$base_key};
				//$pr_not_base_hash->{$hypothetical_base} += $log10_error_rates->[$fastq_file_index]->{$quality}->{$base_key};
				
				//std::cerr << base2char(base) << " " << base_key << " " << (unsigned int)quality << " " << _ecr.log10_correct_rates(fastq_file_index, quality, base_key) << " " << _ecr.log10_error_rates(fastq_file_index, quality, base_key) << std::endl;
				
				std::pair<double,double> log10rates = _ecr.log10_rates(fastq_file_index, quality, base_key);
				//std::cerr << "pr: " << log10rates.second << " not: " << log10rates.first << std::endl;
				pr_base_hash[base_list[j]] += log10rates.second; //_ecr.log10_correct_rates(fastq_file_index, quality, base_key);
				pr_not_base_hash[base_list[j]] += log10rates.first; //_ecr.log10_error_rates(fastq_file_index, quality, base_key);
			}
		} // end for-each read
		
		//## PER POSITION/INSERT COUNT
		//		
		//#sum up coverage observations
		//		$this_position_coverage->{unique}->{total} = $this_position_coverage->{unique}->{-1} + $this_position_coverage->{unique}->{+1};
		//		$this_position_coverage->{redundant}->{total} = $this_position_coverage->{redundant}->{-1} + $this_position_coverage->{redundant}->{+1};
		//		$this_position_coverage->{raw_redundant}->{total} = $this_position_coverage->{raw_redundant}->{-1} + $this_position_coverage->{raw_redundant}->{+1};
		//		$this_position_coverage->{total} = $this_position_coverage->{unique}->{total} + $this_position_coverage->{redundant}->{total};
		//		
		this_position_coverage.sum();
		
		//		$s->{coverage}->{unique_total}++ if ($this_position_unique_only_coverage && ($insert_count == 0));
		if(this_position_unique_only_coverage && (insert_count==0)) {
			//		$s->{coverage}->{unique_total}++			
			++s.coverage_unique_total; // @dk: this could be skipped for now?
		}
		
		//		std::cerr << this_position_coverage.total << std::endl;
		
		//#we are trying to find the base with the most support;		
		//#calculate ratios of each base to all other bases
		//		my $best_base;
		//		my $pr_call;
		//		foreach my $test_base (keys %$pr_base_hash)
		uint8_t best_base=0;
		double pr_call=0; 
		bool pr_call_defined=false;
		for(map<uint8_t,double>::iterator i=pr_base_hash.begin(); i!=pr_base_hash.end(); ++i) {
			uint8_t test_base=i->first;
			
			//		{			
			//##obsolete code
			//# 'P' is a special key for evidence against indels
			//###next if ($test_base eq 'P');
			//			
			//			my $this_pr_call = $pr_base_hash->{$test_base} - $pr_not_base_hash->{$test_base};
			//			if ( (!defined $pr_call) || ($this_pr_call > $pr_call) )
			//			{
			//				$pr_call = $this_pr_call;
			//				$best_base = $test_base;
			//			}
			//		}
			double this_pr_call = i->second - pr_not_base_hash[test_base];
			//std::cerr << base2char(test_base) << " " << this_pr_call << " " << i->second << " " << pr_not_base_hash[test_base] << std::endl;
			if((!pr_call_defined) || (this_pr_call > pr_call)) {
				pr_call = this_pr_call; 
				pr_call_defined = true;
				best_base = test_base;
				//std::cerr << base2char(test_base) << " " << pr_call << std::endl;
			}
		}
		
		//#for the case where there are no counts
		//		my $e_value_call = 'NA';
		double e_value_call=std::numeric_limits<double>::quiet_NaN();
		
		//#otherwise correct to an e-value based on genome size
		//#could correct only to number of unique positions???
		//		if (defined $pr_call)
		//		{
		//			$e_value_call = $pr_call - $log10_ref_length;
		//			$e_value_call = sprintf "%.1f", $e_value_call; #round immediately
		//		}
		if(pr_call_defined) {
			e_value_call = pr_call - _log10_ref_length;
			// @dk: based on the code above, this value should be truncated, but based on comments above it's rounded...
			e_value_call = 0.1 * round(e_value_call*10);
//			std::cerr << pr_call << " " << e_value_call << std::endl;
		}
		
		//##did we predict a base at this position?
		//		my $base_predicted = ($e_value_call ne 'NA' && ($e_value_call >= $settings->{mutation_log10_e_value_cutoff}));
		bool base_predicted=false;
		if(!std::isnan(e_value_call) && (e_value_call >= _mutation_cutoff)) {
			base_predicted = true;
		}
		
		//std::cerr << position << " e:" << e_value_call << " b:" << base_predicted << std::endl;
		
		//##print out SNP call information
		//		my $total_cov;
		//		$total_cov->{1} = 0;
		//		$total_cov->{-1} = 0;
		//		
		//		$line = "$pos\t$insert_count\t$ref_base\t$e_value_call";
		//		foreach my $base (@base_list)
		//		{			
		//			my $current_base_info = $pos_info->{$base};
		//			my $top_cov = $current_base_info->{unique_trimmed_cov}->{1};
		//			my $bot_cov = $current_base_info->{unique_trimmed_cov}->{-1};
		//			$total_cov->{1} += $top_cov;
		//			$total_cov->{-1} += $bot_cov;
		//			$line .= "\t$base\t" . "\t($bot_cov/$top_cov)";
		//		}
		//		print MUT "$line\n" if (defined $snps_all_tab_file_name);
		
		int total_cov[3]={0,0,0}; // triple, same as above
    
    // Don't need to print, but is nice for debug
		//ostringstream line;
		//line << position << " " << insert_count << " " << ref_base << " " << e_value_call;
		
		for(std::size_t j=0; j<5; ++j) {
			double top_cov = pos_info[base_list[j]].unique_trimmed_cov[2];
			double bot_cov = pos_info[base_list[j]].unique_trimmed_cov[0];
			total_cov[2] += round(top_cov);
			total_cov[0] += round(bot_cov);
		//	line << " " << base_list[j] << " (" << bot_cov << "/" << top_cov << ")";
		}
		//std::cerr << line.str() << endl; // @dk print to file
		

		
		
		//###
		//## DELETION DELETION DELETION
		//###
		//		
		//#print to coverage file
		//#update information on deletions
		//		if ($insert_count == 0)
		//		{
		//			if (!$settings->{no_deletion_prediction})
		//			{
		//				_check_deletion_completion($pos, $this_position_coverage, $e_value_call);
		//				_update_copy_number_variation($pos, $this_position_coverage, $ref_base); 							
		//			}
		//		}			
		if(insert_count == 0) {
			if(_predict_deletions) {
        // @JEB: note change in call so position sent to check_deletion_completion is 1-based
				check_deletion_completion(position, p.target(), &this_position_coverage, e_value_call);
				// @dk: skip update_copy_number_variation(pos, this_position_coverage, ref_base);
			}
		}
		
		//###
		//## POLYMORPHISM POLYMORPHISM POLYMORPHISM
		//###								
		//		my $polymorphism_predicted = 0;
		//		my $polymorphism;
		//		if ($settings->{polymorphism_prediction})
		bool polymorphism_predicted=false;
		if(_predict_polymorphisms) {
			//			$polymorphism = _predict_polymorphism($settings, $pdata, $log10_correct_rates, $error_rates, $ref_base);
			//			
			//			if ($polymorphism)
			//			{
			//				$polymorphism->{log10_e_value} = 'ND'; 
			//				if ($polymorphism->{p_value} ne 'ND')
			//				{
			//					$polymorphism->{log10_e_value} = ($polymorphism->{p_value} == 0) ? "999" : -log($total_ref_length * $polymorphism->{p_value})/log(10);						
			//				}
			//				if ($polymorphism->{log10_e_value} >= 2)
			//				{
			//					$polymorphism_predicted = 1;
			//#print Dumper($polymorphism);
			//				}
			//				$base_predicted = 1 if ($polymorphism_predicted);
			//			}					
		}				
		
		
		
		//###
		//## UNKNOWN UNKNOWN UNKNOWN
		//###				
		//		if ($insert_count == 0)
		if(insert_count == 0) {
			//			_update_unknown_intervals($seq_id, $pos, $base_predicted, $this_position_unique_only_coverage);
			update_unknown_intervals(position, p.target(), base_predicted, this_position_unique_only_coverage);
		}
		
		
		//## evaluate whether to call an actual mutation!				
		//### skip if there is not enough evidence for a call or if it agrees with the reference
		//#	next if (!$base_predicted);	
		//		next if (($e_value_call eq 'NA') || ($e_value_call < -$log10_ref_length));
		if(std::isnan(e_value_call) || (e_value_call < -_log10_ref_length)) {
			continue;
		}
		
		//std::cerr << e_value_call << std::endl;
		
		//## mutation and polymorphism are exclusive predictions.
		// my $mutation_predicted = (!$polymorphism_predicted) && ($best_base ne $ref_base);
		bool mutation_predicted = !polymorphism_predicted && (base2char(best_base) != base2char(ref_base));
				
		//## bail if it's just the reference base and we aren't interested in polymorphisms...
		//		next INSERT_COUNT if (!$mutation_predicted && !$polymorphism_predicted);
		if(!mutation_predicted && !polymorphism_predicted) {
			continue;
		}
		//## bail if we are predicting polymorphisms, but there wasn't one
		
		//std::cerr << position << " " << e_value_call << " " << mutation_predicted << " " << polymorphism_predicted << " " << base2char(best_base) << " " << base2char(ref_base) << std::endl;
		

		
		//## Fields common to consensus mutations and polymorphisms
		//my $mut;				
		//$mut->{type} = 'RA';
		//$mut->{seq_id} = $seq_id;
		//$mut->{position} = $pos;
		//$mut->{insert_position} = $insert_count;
		//$mut->{quality} = $e_value_call;		
		ra mut(boost::lexical_cast<std::string>(_gd.new_id()), "");
		mut[SEQ_ID] = p.target_name();
		mut[POSITION] = position;
		mut[INSERT_POSITION] = insert_count;
		mut[QUALITY] = e_value_call;
		
		//
		//## code that prints out even more information
		//## slow because it sorts things, and not necessary
		//if (0)
		//{
		//my ($base_string, $quality_string, $strand_string) = _pdata_to_strings(@$pdata);
		//$mut->{bases} = $base_string;
		//$mut->{qualities} = $quality_string;
		//$mut->{strands} = $strand_string;
		//}
		
		if(mutation_predicted) {
			//$mut->{ref_base} = $ref_base;
			//$mut->{new_base} = $best_base;		
			//$mut->{frequency} = 1; ## this is not a polymorphism
			//$mut->{reject} = "EVALUE" if ($e_value_call < $settings->{mutation_log10_e_value_cutoff});
			mut[REF_BASE] = ref_base;
			mut[NEW_BASE] = best_base;
			mut[FREQUENCY] = 1;
			if(e_value_call < _mutation_cutoff) {
				mut[REJECT] = "EVALUE";
			}
		}
		
		
		//if ($polymorphism_predicted)
		if(polymorphism_predicted) {
			//$mut->{quality} = $polymorphism->{log10_e_value};		
			//#$mut->{fisher_strand_p_value} = $polymorphism->{fisher_strand_p_value};
			
			//# the frequency returned is the probability of the FIRST base
			//# we want to quote the probability of the second base (the change from the reference).
			//my $polymorphism_coverage_both_bases = 0;
			//if ($polymorphism->{first_base} eq $ref_base)
			//{
			//$mut->{frequency} = 1-$polymorphism->{frequency};
			//$mut->{ref_base} = $polymorphism->{first_base};
			//$mut->{new_base} = $polymorphism->{second_base};
			//$polymorphism_coverage_both_bases = 
			//( ($polymorphism->{second_base_strand_coverage}->{-1} > 0)
			//&& ($polymorphism->{second_base_strand_coverage}->{+1} > 0) );
			//}	
			//elsif ($polymorphism->{second_base} eq $ref_base)
			//{
			//$mut->{frequency} = $polymorphism->{frequency};
			//$mut->{ref_base} = $polymorphism->{second_base};
			//$mut->{new_base} = $polymorphism->{first_base};
			//
			//$polymorphism_coverage_both_bases = 
			//( ($polymorphism->{first_base_strand_coverage}->{-1} > 0)
			//&& ($polymorphism->{first_base_strand_coverage}->{+1} > 0) );					
			//}
			//### NOTE: This neglects the case where neither the first nor second base is the reference base! Should almost never happen					
			//# die if (($polymorphism->{first_base} ne $ref_base) && ($polymorphism->{second_base} ne $ref_base));
			//else
			//{
			//$mut->{frequency} = $polymorphism->{frequency};
			//$mut->{ref_base} = $polymorphism->{first_base};
			//$mut->{new_base} = $polymorphism->{second_base};
			//
			//$mut->{error} = "polymorphic_without_reference_base";
			//}
			
			//$mut->{reject} = "EVALUE" if ($mut->{quality} < $settings->{polymorphism_log10_e_value_cutoff});
			
			//###
			//## Print input file for R
			//###
			//my $ref_cov = $pos_info->{$mut->{ref_base}}->{unique_trimmed_cov};
			//my $new_cov = $pos_info->{$mut->{new_base}}->{unique_trimmed_cov};
			//
			//my @ref_base_qualities;
			//my @new_base_qualities;
			//foreach my $item (@$pdata)
			//{
			//if ($item->{base} eq $mut->{ref_base})
			//{
			//push @ref_base_qualities, $item->{quality};
			//}
			//elsif ($item->{base} eq $mut->{new_base})
			//{
			//push @new_base_qualities, $item->{quality};
			//}
			//}
			//my $ref_quality_string = join ',', @ref_base_qualities;
			//my $new_quality_string = join ',', @new_base_qualities;
			//
			//print $polymorphism_statistics_input_fh +join( "\t",
			//$new_cov->{1}, $new_cov->{-1}, $ref_cov->{1}, $ref_cov->{-1}, $new_quality_string, $ref_quality_string
			//) . "\n";
			//###
			//## End printing input file for R
			//###
			
			//$mut->{reject} = "STRAND" if (!$polymorphism_coverage_both_bases);
			//$mut->{reject} = "FREQ" if ($mut->{frequency} < $settings->{polymorphism_frequency_cutoff});
			//$mut->{reject} = "FREQ" if ($mut->{frequency} > 1-$settings->{polymorphism_frequency_cutoff});		
		}
		
		//## More fields common to consensus mutations and polymorphisms
		//## ...now that ref_base and new_base are defined
		//my $ref_cov = $pos_info->{$mut->{ref_base}}->{unique_trimmed_cov};
		//$mut->{ref_cov} = $ref_cov->{-1} . "/" . $ref_cov->{1};
		int* ref_cov = pos_info[boost::get<uint8_t>(mut[REF_BASE])].unique_trimmed_cov;
		mut[REF_COV] = std::make_pair(ref_cov[2], ref_cov[0]);
		
		//my $new_cov = $pos_info->{$mut->{new_base}}->{unique_trimmed_cov};
		//$mut->{new_cov} = $new_cov->{-1} . "/" . $new_cov->{1};
		int* new_cov = pos_info[boost::get<uint8_t>(mut[NEW_BASE])].unique_trimmed_cov;
		mut[NEW_COV] = std::make_pair(new_cov[2], new_cov[0]);
		
		//$mut->{tot_cov} = $total_cov->{-1} . "/" . $total_cov->{1};
		mut[TOT_COV] = std::make_pair(total_cov[2], total_cov[0]);
		
		//$gd->add($mut);
		_gd.add(mut);
	}
}



/*! Called at the end of the pileup.
 */
void breseq::identify_mutations_pileup::at_end(uint32_t tid, uint32_t seqlen) {
	//	_check_deletion_completion($sequence_length+1); 		
	//	_update_unknown_intervals($seq_id, $sequence_length+1, 1);

  // @JEB changed check_deletion_completion to 1-based
	check_deletion_completion(seqlen+1, tid, 0, std::numeric_limits<double>::quiet_NaN());
	update_unknown_intervals(seqlen+1, tid, true, false);

	//	my $ra_mc_genome_diff_file_name = $settings->file_name('ra_mc_genome_diff_file_name');	
	//	$gd->write($ra_mc_genome_diff_file_name);
	_gd.write();
	
	_coverage_data.close();
}


/*! Helper method to track information about putative deleted regions.
 
 Used at each pileup iteration and at the end.
 //## when called at the end of a fragment, the position is fragment length +1
 //## and $this_position_coverage is undefined
 
 @JEB This function expects 1-indexed positions!!!
 
 */
void breseq::identify_mutations_pileup::check_deletion_completion(uint32_t position, uint32_t seq_id, const position_coverage* this_position_coverage, double e_value_call) {

	//std::cerr << position << " " << e_value_call << std::endl;
		
  // reset if we hit a new fragment  
  if (seq_id != _on_deletion_seq_id)
  {
    _last_position_coverage_printed = 0;
    _on_deletion_seq_id = seq_id;
  }  
    
	//# we need to fill in reference positions with NO reads aligned to them
	//# pileup won't be called at these positions
	//foreach (my $i = $last_position_coverage_printed + 1; $i < $pos; $i++)
	//{
	//std::cerr << _last_position_coverage_printed << " " << position << std::endl;

	for(uint32_t i=_last_position_coverage_printed+1; i<position; ++i) {
		if(_last_deletion_start_position == boost::none) {
			//## special treatment for the beginning of a fragment
			if(_last_position_coverage_printed == 0) {
				//$left_outside_coverage_item =  {
				//unique => {'1'=>'NA', '-1'=>'NA', 'total' => 'NA' },
				//redundant => {'1'=>'NA', '-1'=>'NA', 'total' => 'NA' }				//};
				//$left_inside_coverage_item = { 
				//unique => {'1'=>0, '-1'=>0, 'total' => 0 },
				//redundant => {'1'=>0, '-1'=>0, 'total' => 0 }				//};
				_left_outside_coverage_item.reset(position_coverage(std::numeric_limits<double>::quiet_NaN()));
				_left_inside_coverage_item.reset(position_coverage());				
			} else {
				//## normal treatment is that coverage went to zero
				//$left_outside_coverage_item = { 
				//unique => {'1'=>0, '-1'=>0, 'total' => 0 },
				//redundant => {'1'=>0, '-1'=>0, 'total' => 0 }				//};
				_left_inside_coverage_item.reset(position_coverage());
				
				//$left_inside_coverage_item = $last_position_coverage;	
				_left_outside_coverage_item = _last_position_coverage;
			}
			
			// moved into the ifs above, to handle the special case of the fragment beginning.
			// had to do this because of the change in 0-1 indexing.
			// $last_deletion_start_position = $last_position_coverage_printed+1;			
			_last_deletion_start_position = _last_position_coverage_printed+1;
      if (_last_deletion_redundant_start_position == boost::none) {
        _last_deletion_redundant_start_position = _last_position_coverage_printed+1;
      }
    }
		
		//$this_deletion_reaches_seed_value = 1;
		_this_deletion_reaches_seed_value = true;
    _this_deletion_redundant_reached_zero = true;
		//$last_position_coverage = { 
		//unique => {'1'=>0, '-1'=>0, 'total' => 0 },
		//redundant => {'1'=>0, '-1'=>0, 'total' => 0 }		//};
		_last_position_coverage.reset(position_coverage());
		
		//print COV join("\t", 0, 0, 0, 0, 0, 0, 'NA', $i) . "\n"; 
		//_coverage_data << 0 << "\t" << 0 << "\t" << 0 << "\t" << 0 << "\t" << 0 << "\t" << 0 << "\t" << "NA\t" << i << std::endl;
	}
	_last_position_coverage_printed = position;
	
	//## called with an undef $this_position_coverage at the end of the genome
	//if ((defined $this_position_coverage) && (defined $e_value_call))
	if(this_position_coverage) {
		//my $tu = $this_position_coverage->{unique};
		//my $tr = $this_position_coverage->{redundant};
		//my $trr = $this_position_coverage->{raw_redundant};
		
		//#print this information
		//print COV join("\t", $tu->{-1}, $tu->{1}, $tr->{-1}, $tr->{1}, $trr->{-1}, $trr->{1}, $e_value_call, $pos) . "\n";
		//_coverage_data << this_position_coverage->unique[0] << "\t"
		//<< this_position_coverage->unique[2] << "\t"
		//<< this_position_coverage->redundant[0] << "\t"
		//<< this_position_coverage->redundant[2] << "\t"
		//<< this_position_coverage->raw_redundant[0] << "\t"
		//<< this_position_coverage->raw_redundant[2] << "\t"
		//<< e_value_call << "\t" << position << std::endl;
		
    //## UNIQUE COVERAGE
		//#start a new possible deletion if we fall below the propagation cutoff
		//if ($this_position_coverage->{unique}->{total} <= $deletion_propagation_cutoff)
		if(this_position_coverage->unique[1] <= _this_deletion_propagation_cutoff) {
			//if (!defined $last_deletion_start_position)
			if(_last_deletion_start_position == boost::none) {
				//$last_deletion_start_position = $pos;
				_last_deletion_start_position = position;
				//$left_outside_coverage_item = $last_position_coverage;
				_left_outside_coverage_item = _last_position_coverage;
				//$left_inside_coverage_item = $this_position_coverage;
				_left_inside_coverage_item = *this_position_coverage;
			}
		}
		
		//##keep track of whether we've encountered the seed value
		//		if ($this_position_coverage->{total} <= $deletion_seed_cutoff)
		if(this_position_coverage->total <= _deletion_seed_cutoff) {
			//			$this_deletion_reaches_seed_value = 1;
			_this_deletion_reaches_seed_value = true;
		}
    
    //## REDUNDANT COVERAGE
    //## updated only if we are currently within a deletion
    if (_last_deletion_start_position != boost::none) {
    //if (defined $last_deletion_start_position)
    //{					
      if (this_position_coverage->redundant[1] == 0) {
      //if ($this_position_coverage->{redundant}->{total} == 0)
      //{
          _this_deletion_redundant_reached_zero = true;
        //$redundant_reached_zero = 1; #switch from adjusting start to end
          _last_deletion_redundant_end_position = boost::none;
        //undef $last_deletion_redundant_end_position;
        }
      //}
      //elsif ($this_position_coverage->{redundant}->{total} > 0)
        else if (this_position_coverage->redundant[1] > 0) {
      //{
        //## if there is any redundant coverage remember the start (until we find zero redundant coverage)
        //if (!$redundant_reached_zero)
          if (!_this_deletion_redundant_reached_zero) {
            _last_deletion_redundant_start_position = position;
          }
          else {
            if (_last_deletion_redundant_end_position == boost::none) _last_deletion_redundant_end_position = position;
          }
        //{
          //$last_deletion_redundant_start_position = $pos;
        //}
        //## if we are working on the right side update the end position if it is not already defined.
        //else
        //{
          //$last_deletion_redundant_end_position = $pos if (!defined $last_deletion_redundant_end_position);
        //}
      //}
    //}
        }
    }
    
	}
	
	//##if we are above the propagation cutoff then record the current deletion
	//	if ( (defined $last_deletion_start_position) && ((!defined $this_position_coverage) 
	//																									 || ($this_position_coverage->{unique}->{total} > $deletion_propagation_cutoff)) )
	if(_last_deletion_start_position != boost::none
		 && (!this_position_coverage || (this_position_coverage->unique[1] > _this_deletion_propagation_cutoff))) {
		
		//		if ($this_deletion_reaches_seed_value)
		if(_this_deletion_reaches_seed_value) {
			boost::optional<position_coverage> tmp;
			
			// ### for the end of the genome....
			//			if (!defined $this_position_coverage)
			if(!this_position_coverage) {
				//				$this_position_coverage = {
				//					unique => {'1'=>'NA', '-1'=>'NA', 'total' => 'NA' },
				//					redundant => {'1'=>'NA', '-1'=>'NA', 'total' => 'NA' }				//				};
				tmp.reset(position_coverage(std::numeric_limits<double>::quiet_NaN()));
			} else {
				tmp.reset(*this_position_coverage);
			}

      _last_deletion_end_position = position-1;
      //my $last_deletion_end_position = $pos-1;
      if (_last_deletion_redundant_end_position == boost::none) _last_deletion_redundant_end_position = _last_deletion_end_position;
      if (_last_deletion_redundant_start_position == boost::none) _last_deletion_redundant_start_position = _last_deletion_start_position;

      //$last_deletion_redundant_end_position = $last_deletion_end_position	if (!defined $last_deletion_redundant_end_position);
      //$last_deletion_redundant_start_position = $last_deletion_start_position	if (!defined $last_deletion_redundant_start_position);

			
			//			my $del = {
			//				type => 'MC',
			//				seq_id => $seq_id,
			//				start => $last_deletion_start_position,
			//				end => $pos-1,
			//        start_range => $last_deletion_redundant_start_position - $last_deletion_start_position,
			//        end_range => $last_deletion_end_position - $last_deletion_redundant_end_position,
			//#			size => ($pos-1) - $last_deletion_start_position + 1, #end - start + 1
			//				left_outside_cov => $left_outside_coverage_item->{unique}->{total},
			//				left_inside_cov => $left_inside_coverage_item->{unique}->{total},
			//				right_inside_cov => $last_position_coverage->{unique}->{total},
			//				right_outside_cov => $this_position_coverage->{unique}->{total},
			//			};
			mc del(boost::lexical_cast<std::string>(_gd.new_id()), "");
			del[SEQ_ID] = target_name(seq_id);
			del[START] = *_last_deletion_start_position;
			del[END] = *_last_deletion_end_position;
			del[START_RANGE] = *_last_deletion_redundant_start_position - *_last_deletion_start_position;
			del[END_RANGE] = *_last_deletion_end_position - *_last_deletion_redundant_end_position;
			
			if(!_left_outside_coverage_item || std::isnan(_left_outside_coverage_item->unique[1])) {
				del[LEFT_OUTSIDE_COV] = std::numeric_limits<double>::quiet_NaN();
			} else {
				del[LEFT_OUTSIDE_COV] = static_cast<uint32_t>(_left_outside_coverage_item->unique[1]);
			}

			if(!_left_inside_coverage_item || std::isnan(_left_inside_coverage_item->unique[1])) {
				del[LEFT_INSIDE_COV] = std::numeric_limits<double>::quiet_NaN();
			} else {
				del[LEFT_INSIDE_COV] = static_cast<uint32_t>(_left_inside_coverage_item->unique[1]);
			}
			
			if(!_last_position_coverage || std::isnan(_last_position_coverage->unique[1])) {
				del[RIGHT_INSIDE_COV] = std::numeric_limits<double>::quiet_NaN();
			} else {
				del[RIGHT_INSIDE_COV] = static_cast<uint32_t>(_last_position_coverage->unique[1]);
			}
			
			if(!this_position_coverage || std::isnan(this_position_coverage->unique[1])) {
				del[RIGHT_OUTSIDE_COV] = std::numeric_limits<double>::quiet_NaN();
			} else {
				del[RIGHT_OUTSIDE_COV] = static_cast<uint32_t>(this_position_coverage->unique[1]);
			}
			
//			del[LEFT_OUTSIDE_COV] = _left_outside_coverage_item ? _left_outside_coverage_item->unique[1] : std::numeric_limits<double>::quiet_NaN();
//			del[LEFT_INSIDE_COV] = _left_inside_coverage_item ? _left_inside_coverage_item->unique[1] : std::numeric_limits<double>::quiet_NaN();
//			del[RIGHT_INSIDE_COV] = _last_position_coverage ? _last_position_coverage->unique[1] : std::numeric_limits<double>::quiet_NaN();
//			del[RIGHT_OUTSIDE_COV] = this_position_coverage ? this_position_coverage->unique[1] : std::numeric_limits<double>::quiet_NaN();
			
			//			$del->{left_inside_cov} = 'NA' if (!defined $del->{left_inside_cov});
			//			$del->{right_inside_cov} = 'NA' if (!defined $del->{right_inside_cov});
			//			
			//			$del->{left_outside_cov} = 'NA' if (!defined $del->{left_outside_cov});
			//			$del->{right_outside_cov} = 'NA' if (!defined $del->{right_outside_cov});
			//			
			//			$gd->add($del);					
			
			_gd.add(del);
		}
		
		//#reset the search
		//		$this_deletion_reaches_seed_value = 0;
		_this_deletion_reaches_seed_value = false;
		//		undef $last_deletion_start_position;
    _this_deletion_redundant_reached_zero = false;
		_last_deletion_start_position = boost::none;
    _last_deletion_end_position = boost::none;
		_last_deletion_redundant_start_position = boost::none;
    _last_deletion_redundant_end_position = boost::none;
	}
	
	//	$last_position_coverage = $this_position_coverage;
	if(this_position_coverage) {
		_last_position_coverage = *this_position_coverage;
	}
}


//
//sub _update_unknown_intervals
//{
//	my ($seq_id, $pos, $base_predicted, $this_position_unique_only_coverage) = @_;
//	

/*! Helper method to track unknowns.
 */
void breseq::identify_mutations_pileup::update_unknown_intervals(uint32_t position, uint32_t seq_id, bool base_predicted, bool this_position_unique_only_coverage) {
	//std::cerr << position << " " << base_predicted << " " << this_position_unique_only_coverage << std::endl;
	//if(_last_start_unknown_interval) {
	//	std::cerr << *_last_start_unknown_interval << std::endl;
	//} else {
	//	std::cerr << "undef" << std::endl;
	//}
	
	//	if (!$base_predicted)
	if(!base_predicted) {
		//		$s->{coverage}->{unique_uncalled}++ if (($this_position_unique_only_coverage));
		if(this_position_unique_only_coverage) {
			++s.coverage_unique_uncalled;
		}
		//		if (!defined $last_start_unknown_interval)
		if(!_last_start_unknown_interval) {
			//			$last_start_unknown_interval = $pos;
			_last_start_unknown_interval = position;
		}
	}	else {
		//		$s->{coverage}->{unique_called}++ if (($this_position_unique_only_coverage));
		if(this_position_unique_only_coverage) {
			++s.coverage_unique_called;
		}
			
		//#end interval where we were unable to call mutations
		//		if (defined $last_start_unknown_interval)
		if(_last_start_unknown_interval) {
			//			my $new_interval = { 'type'=>'UN', 'start'=> $last_start_unknown_interval, 'end'=> $pos-1, 'seq_id' => $seq_id };
			un new_interval(boost::lexical_cast<std::string>(_gd.new_id()), "");
			new_interval[SEQ_ID] = target_name(seq_id);
			new_interval[START] = *_last_start_unknown_interval;
			new_interval[END] = position - 1;
			_gd.add(new_interval);
			
			_last_start_unknown_interval = boost::none;
			//#	print Dumper($new_interval); ##DEBUG
		}
	}
}
