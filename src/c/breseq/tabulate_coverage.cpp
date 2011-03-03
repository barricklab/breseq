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
#include <boost/math/distributions.hpp>

#include "breseq/common.h"
#include "breseq/tabulate_coverage.h"


/*! Convenience wrapper around the identify_mutations_pileup class.
 */
void breseq::tabulate_coverage( const std::string& bam, 
																const std::string& fasta,
                                const std::string& output,
                                const std::string& region,
                                const uint32_t downsample )
{

  tabulate_coverage_pileup tcp(bam, fasta, output);
  tcp.do_pileup();
}

/*! Constructor.
 */
breseq::tabulate_coverage_pileup::tabulate_coverage_pileup(const std::string& bam, const std::string& fasta, const std::string& output)
: breseq::pileup_base(bam, fasta), _last_position_processed(0) {

  _output_table.open(output.c_str());
  
  _output_table << "position" << "\t" << "ref_base" << "\t" 
    << "unique_top_cov" << "\t" << "unique_bot_cov" << "\t" 
    << "redundant_top_cov" << "\t" << "redundant_bot_cov" << "\t" 
    << "raw_redundant_top_cov" << "\t" << "raw_redundant_bot_cov" << "\t"
    << "unique_top_begin" << "\t" << "unique_bot_begin"
  << std::endl;
}


/*! Destructor.
 */
breseq::tabulate_coverage_pileup::~tabulate_coverage_pileup() {
}

/*! Called for each alignment.
 */
void breseq::tabulate_coverage_pileup::callback(const breseq::pileup& p) {

  char* refseq = p.reference_sequence(); // reference sequence for this target
  uint32_t pos = p.position();
  
  // don't handle indels before first position
  if (pos==0) return;
  
  // print positions not called because there were no reads
  for (uint32_t i=_last_position_processed+1; i<pos; i++) {
    _output_table << i << "\t" << refseq[i-1] << "\t" << 0 << "\t" << 0 << "\t" 
      << 0 << "\t" << 0 << "\t" << 0 << "\t" << 0 << "\t" << 0 << "\t" << 0 << std::endl;
  }
  
  uint8_t ref_base = refseq[pos-1];
  uint32_t unique_cov[2] = {0,0};
  double redundant_cov[2] = {0.0, 0.0};
  uint32_t raw_redundant_cov[2] = {0,0};
  uint32_t unique_begin_reads[2] = {0,0};
  
	// for each alignment within this pileup:
	for(pileup::const_iterator i=p.begin(); i!=p.end(); ++i) {
    
    // skip deletions entirely, they are handled by adjacent matching positions
    if(i->is_del()) {
      continue;
    }

		uint32_t redundancy = i->redundancy();
    uint32_t reversed = i->reversed();
    bool first_base_matched;
    bool this_is_first_base;
    if (!reversed) { 
      // top strand
      first_base_matched = (i->query_start() == 1);
      this_is_first_base = (i->query_position() == 0); 
    } else {
      // bottom strand
      first_base_matched = (i->query_end() == i->query_length());
      this_is_first_base = (i->query_position()+1 == i->query_length());       
    }

// WHOA -- do we really want this....?    
    if (!first_base_matched) {
      continue;
    }
    
		if (redundancy == 1)
		{	
			unique_cov[reversed]++;
      if (this_is_first_base) {
        unique_begin_reads[reversed]++;
      } 
		}
		else
		{
      raw_redundant_cov[reversed]++;
			redundant_cov[reversed] += 1.0/redundancy;			
		}		
  }
  
  
  //output
  _output_table << pos << "\t" << ref_base << "\t" 
    << unique_cov[0] << "\t" << unique_cov[1] << "\t" 
    << redundant_cov[0] << "\t" << redundant_cov[1] << "\t" 
    << raw_redundant_cov[0] << "\t" << raw_redundant_cov[1] << "\t"
    << unique_begin_reads[0] << "\t" << unique_begin_reads[1] 
    
  << std::endl;
 
  _last_position_processed = pos;
}

/*! Called at the end of the pileup.
 */
void breseq::tabulate_coverage_pileup::at_end(uint32_t tid, uint32_t seqlen) {

  char* refseq = get_refseq(tid); // reference sequence for this target
  uint32_t pos = seqlen+1;
  
  // print positions not called because there were no reads
  for (uint32_t i=_last_position_processed+1; i<pos; i++) {
    _output_table << i << "\t" << refseq[i-1] << "\t" << 0 << "\t" << 0 << "\t" 
      << 0 << "\t" << 0 << "\t" << 0 << "\t" << 0 << "\t" << 0 << "\t" << 0 << std::endl;
  }
  
  _last_position_processed = pos;
}



/*
sub tabulate_coverage
{
	my ($self, $tmp_coverage, $region, $options) = @_;
	my $downsample = $options->{downsample};
	$downsample = 1 if (!defined $downsample);
	my $bam = $self->{bam};
		
	my ($seq_id, $start, $end, $insert_start, $insert_end);
	($seq_id, $start, $end, $insert_start, $insert_end, $region) = Breseq::Shared::region_to_coords($region, $bam);
	
	## Open file for output
	open COV, ">$tmp_coverage";
	COV->autoflush(1);
	print COV join("\t", 'position', 'unique_top_cov', 'unique_bot_cov', 'redundant_top_cov', 'redundant_bot_cov') . "\n";
	our $coverage;

	###
	##  Behold, the dreaded SLOW fetch function...
	###
	
	my $fetch_function = sub {
		my ($a) = @_;
  						
		my $redundancy = $a->aux_get('X1');
		my $reversed = $a->reversed;
		my $strand = $reversed ? -1 : +1;
		
		##### Update coverage if this is not a deletion in read relative to reference
		### Count trimmed reads here, but not when looking for short indel mutations...	
		if ($redundancy == 1)
		{	
			$coverage->{unique}->{$strand}++;
		}
		else
		{
			$coverage->{redundant}->{$strand} += 1/$redundancy;			
#			$this_position_coverage->{raw_redundant}->{$strand}++;			
		}		
	}; 

	$start = $downsample * int($start/$downsample);
	$start = 1 if ($start < 1);
	
	$end = $downsample * int((($end-1)/$downsample) + 1);
	my $reference_length = $bam->length($seq_id);
	$end = $reference_length if ($end > $reference_length);	
				
	for (my $pos = $start; $pos <= $end; $pos += $downsample)
	{		
		#initialize coverage observations
		$coverage->{unique} = {'-1' => 0, '1' => 0, 'total' => 0};
		$coverage->{redundant} = {'-1' => 0, '1' => 0, 'total' => 0};
#		$coverage->{raw_redundant} = {'-1' => 0, '1' => 0, 'total' => 0};
				
		my $fetch_region = $seq_id . ":" . $pos . "-" . $pos;
#		print "$fetch_region\n";
		$bam->fetch($fetch_region, $fetch_function);
		
		#sum up coverage observations
		$coverage->{unique}->{total} = $coverage->{unique}->{-1} + $coverage->{unique}->{+1};
		$coverage->{redundant}->{total} = $coverage->{redundant}->{-1} + $coverage->{redundant}->{+1};
#		$coverage->{raw_redundant}->{total} = $coverage->{raw_redundant}->{-1} + $coverage->{raw_redundant}->{+1};

		my $tu = $coverage->{unique};
		my $tr = $coverage->{redundant};
#		my $trr = $this_position_coverage->{raw_redundant};
		print COV join("\t", $pos, $tu->{-1}, $tu->{1}, $tr->{-1}, $tr->{1}) . "\n";
	}
	
	close COV;
}
*/


