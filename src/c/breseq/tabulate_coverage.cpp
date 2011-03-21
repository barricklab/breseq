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

  std::string output_gc(output);
  output_gc += ".gc.tab";
  std::string output_read_begin(output);
  output_read_begin += ".read_begin.tab";

  tabulate_coverage_pileup tcp(bam, fasta, output, output_read_begin, output_gc);
  tcp.do_pileup();
}

/*! Constructor.
 */
breseq::tabulate_coverage_pileup::tabulate_coverage_pileup(const std::string& bam, const std::string& fasta, const std::string& output,
    const std::string& read_begin_output, const std::string& gc_output)
: breseq::pileup_base(bam, fasta), m_last_position_processed(0) {

  m_output_table.open(output.c_str());
  m_read_begin_output.open(read_begin_output.c_str());
  m_gc_output.open(gc_output.c_str());
  
  
  m_output_table << "position" << "\t" << "ref_base" << "\t" 
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
  uint32_t pos = p.position_1();
  
  // don't handle indels before first position
  if (pos==0) return;
  
  // print positions not called because there were no reads
  for (uint32_t i=m_last_position_processed+1; i<pos; i++) {
    m_output_table << i << "\t" << refseq[i-1] << "\t" << 0 << "\t" << 0 << "\t" 
      << 0 << "\t" << 0 << "\t" << 0 << "\t" << 0 << "\t" << 0 << "\t" << 0 << std::endl;
  }

  // catches this position
  for (uint32_t j=m_last_position_processed+1; j<=pos; j++) {
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
      first_base_matched = (a->query_end_1() == a->query_length());
      this_is_first_base = (a->query_position_1() == a->query_length());     
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
      
      std::string read_begin_s;
      if (!reversed) { 
        for (int i=1; i<=3; i++) {
          read_begin_s += basebam2char(a->query_base_bam_1(i));
        }
        m_read_begin_top_bins[read_begin_s]++;    
      } else {
        for (int i=1; i<=3; i++) {
          read_begin_s += basebam2char(complement_base_bam(a->query_base_bam_1(a->query_length()-i+1)));
        }
        m_read_begin_bot_bins[read_begin_s]++;    
      }

		}
		else
		{
      raw_redundant_cov[reversed]++;
			redundant_cov[reversed] += 1.0/redundancy;			
		}		
  }
  
  
  //output
  m_output_table << pos << "\t" << ref_base << "\t" 
    << unique_cov[0] << "\t" << unique_cov[1] << "\t" 
    << redundant_cov[0] << "\t" << redundant_cov[1] << "\t" 
    << raw_redundant_cov[0] << "\t" << raw_redundant_cov[1] << "\t"
    << unique_begin_reads[0] << "\t" << unique_begin_reads[1] 
    
  << std::endl;
 
  m_last_position_processed = pos;
}

/*! Called at the end of the pileup.
 */
void breseq::tabulate_coverage_pileup::at_end(uint32_t tid, uint32_t seqlen) {

  char* refseq = get_refseq(tid); // reference sequence for this target
  uint32_t pos = seqlen+1;
  
  // print positions not called because there were no reads
  for (uint32_t i=m_last_position_processed+1; i<pos; i++) {
    m_output_table << i << "\t" << refseq[i-1] << "\t" << 0 << "\t" << 0 << "\t" 
      << 0 << "\t" << 0 << "\t" << 0 << "\t" << 0 << "\t" << 0 << "\t" << 0 << std::endl;
  }

  // catches this position
  for (uint32_t j=m_last_position_processed+1; j<pos; j++) {
    if (j>=3) {
      std::string read_begin_s;
      {
        for (int i=1; i<=3; i++) {
          read_begin_s += complement_base_char(reference_base_char_1(tid, j-i+1));
        }
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

  
  m_last_position_processed = pos;
  
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


