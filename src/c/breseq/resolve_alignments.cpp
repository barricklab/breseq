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

#include "breseq/resolve_alignments.h"

#include "breseq/genome_diff.h"
#include "breseq/fastq.h"
#include "breseq/fasta.h"
#include "breseq/alignment.h"
#include "breseq/annotated_sequence.h"

using namespace std;

namespace breseq {

  const string k_junction_name_separator("__");
  
  void junction_name_split(junction_info& ji, const string& junction_name) {
    
    size_t last_pos=0;
    size_t next_pos = junction_name.find(k_junction_name_separator);
    string item;
    int item_num = 0;
    while(next_pos != string::npos) {
      item = junction_name.substr(last_pos, next_pos - last_pos - 1);
      
      switch (item_num) {
        case 0:
          ji.m_side_1.m_seq_id = item;
          break;
          
        case 1:
          ji.m_side_1.m_position = atoi(item.c_str());
          break;      
          
        case 2:
          ji.m_side_1.m_strand = atoi(item.c_str());
           if (ji.m_side_1.m_strand == 0) ji.m_side_1.m_strand = -1;
          break;
          
        case 3:
          ji.m_side_2.m_seq_id = item;
          break;
          
        case 4:
          ji.m_side_2.m_position = atoi(item.c_str());
          break;      
          
        case 5:
          ji.m_side_2.m_strand = atoi(item.c_str());
          if (ji.m_side_2.m_strand == 0) ji.m_side_2.m_strand = -1;
          break;   
          
        case 6:
          ji.m_alignment_overlap = atoi(item.c_str());
          break;
          
        case 7:
          ji.m_unique_read_sequence = item;
          break;      
        
        case 8:
          ji.m_flanking_left = atoi(item.c_str());
          break;      
          
        case 9:
          ji.m_flanking_right = atoi(item.c_str());
          break;
          
        case 10:
          ji.m_side_1.m_redundant = atoi(item.c_str());
          break;
          
        case 11:
          ji.m_side_2.m_redundant = atoi(item.c_str());
          break;
      }
      
      last_pos = next_pos + k_junction_name_separator.length();
      next_pos = junction_name.find(k_junction_name_separator);

    }
  }
  
  
  // Things to FIX
  // 1) Memory allocation and deallocation
  
  //#compare matches to candidate junctions with matches to original genome
  void resolve_alignments( 
                          bool junction_prediction,
                          const string &reference_fasta,
                          const string &junction_fasta,
                          const string &reference_sam_path,
                          const string &junction_sam_path,
                          const string &resolved_path,
                          const string &data_path,
                          const string &features_file,
                          const cReadFiles &read_files,
                          const uint32_t max_read_length,
                          const uint32_t alignment_read_limit
                          ) 
  {
    Settings settings;
    
    int verbose = 1;
    
    // my $verbose = 0;
    // my ($settings, $summary, $ref_seq_info) = @_;
    // my $gene_list_hash_ref = $ref_seq_info->{gene_lists};
    // my $repeat_list_hash_ref = $ref_seq_info->{repeat_lists};
    // my $flanking_length = $settings->{max_read_length};
    
    // ####
    // ##	Reference sequences
    // ####		
        
    // Load the reference sequence info
    cReferenceSequences refseqs;
    LoadFeatureIndexedFastaFile(refseqs, features_file, reference_fasta);    
    
    // ####
    // ##	Junction sequences
    // ####
    
    
    // clean up allocated objects

    // use to get sequences
    // m_seq = fai_fetch(m_ref, target.c_str(), &m_len);

    //## if there were no candidate junctions (file is empty) then we seg fault if we try to use samtools on it...
    //    $settings->{no_junction_prediction} = 1 if ( (!-e $junction_faidx_file_name) || (-s $junction_fasta_file_name == 0) );
    if (junction_prediction && !file_exists(junction_fasta.c_str()) && !file_empty(junction_fasta.c_str())) {
      junction_prediction = 0;
    }
      
    vector<junction_info> junction_info_list;

    //    if (!$settings->{no_junction_prediction})
    if (junction_prediction) {

      // $junction_fai = Bio::DB::Sam::Fai->load($junction_fasta_file_name);	
      ////faidx_t * junction_faidx = fai_load(junction_fasta.c_str());
      ////assert(junction_faidx);

      //      
      //## Load header once at the beginning (but have to peek at TAM file to do this).
      //      my @read_structures = $settings->read_structures;
      //      my $read_file = $read_structures[0]->{base_name};
      //      my $junction_sam_file_name = $settings->file_name('candidate_junction_sam_file_name', {'#'=>$read_file});

      //      my $junction_sam_file_name = $settings->file_name('candidate_junction_sam_file_name', {'#'=>$read_file});

      //      my $junction_tam = Bio::DB::Tam->open($junction_sam_file_name) or die " Could not open junction SAM file\n";
      //      $junction_header = $junction_tam->header_read2($junction_faidx_file_name) or die("Error reading reference fasta index file: $junction_faidx_file_name");				
      //      

      tamFile junction_tam = NULL;
      bam_header_t *junction_header = NULL;
      
      string junction_sam_file_name = junction_sam_path + "/" + read_files[0].m_base_name + ".candidate_junction.sam";

      // $junction_tam = Bio::DB::Tam->open($junction_sam_file_name) or die " Could not open junction SAM file\n";
      junction_tam = sam_open(junction_sam_file_name.c_str());
      assert(junction_tam);
      
      //## junction_header = $junction_tam->header_read2($junction_faidx_file_name) or die("Error reading reference fasta index file: $junction_faidx_file_name");				
      string junction_faidx_file_name = junction_fasta + ".fai";
      junction_header = sam_header_read2(junction_faidx_file_name.c_str()); 

      //## Preload all of the information about junctions
      //## so that we only have to split the names once	
      // my $junction_ids = $junction_header->target_name;
      
      // for (my $i=0; $i< $junction_header->n_targets; $i++)
      for (int i=0; i<junction_header->n_targets; i++) {
        // $junction_info->[$i] = Breseq::Shared::junction_name_split($junction_ids->[$i]);
        junction_info ji;
        junction_name_split(ji, junction_header->target_name[i]);
        junction_info_list.push_back(ji);
      }		
      
      bam_header_destroy(junction_header);
    }
    
    //####
    //##	Output files
    //####

    // our $gd = GenomeDiff->new();
    genome_diff gd;
      
    // my $resolved_reference_sam_file_name = $settings->file_name('resolved_reference_sam_file_name');
    string resolved_reference_sam_file_name = resolved_path + "/reference.sam";
    // my $RREF;
    // open $RREF, ">$resolved_reference_sam_file_name" or die;
    ifstream RREF(resolved_reference_sam_file_name.c_str(), ios_base::in);

    // my $resolved_junction_sam_file_name = $settings->file_name('resolved_junction_sam_file_name');
    string resolved_junction_sam_file_name = resolved_path + "/junction.sam";

    // my $RCJ;
    // open $RCJ, ">$resolved_junction_sam_file_name" or die;
    ifstream RCJ(resolved_junction_sam_file_name.c_str(), ios_base::in);
    
    // my %matched_junction;
    map<string,uint32_t> matched_junction;

    // my %degenerate_matches;
    map<string,uint32_t> degenerate_matches;
        
    // my $reads_processed = 0;
    uint32_t reads_processed = 0;
    
    // keep track of overall index of fastq files 
    uint32_t on_fastq_file_index = 0;
    
    // foreach my $read_struct ($settings->read_structures)
    for (vector<cReadFile>::const_iterator rf = read_files.begin(); rf < read_files.end(); rf++) {
                    
      //print STDERR "  READ FILE::$read_file\n";
      fprintf(stderr, "  READ FILE::%s", rf->m_base_name.c_str());
      cerr << endl;

      // ## setup the summary
      // my $s;
      map<string,uint32_t> summary_info;

      // $s->{unmatched_reads} = 0;
      summary_info["unmatched_reads"] = 0;
      
      //## Traverse the original fastq files to keep track of order
      //## b/c some matches may exist in only one or the other file  
      
      // my @in_fastq;
      // $in_fastq[$i] = Breseq::Fastq->new(-file => $fastq_file_name[$i]);        

      cFastqFile in_fastq(rf->m_fastq_file_name.c_str(), ios::in);
      
      // my @fastq_file_name;
      // my @fastq_file_index;
      
      // my @out_unmatched_fastq;
      cFastqFile out_unmatched_fastq;
            
      // @JEB No longer looping
      // for (my $i=0; $i < scalar @{$read_struct->{base_names}}; $i++)
        
      // my $this_read_file = $read_struct->{base_names}->[$i];        
      // $fastq_file_name[$i] = $settings->read_file_to_fastq_file_name($this_read_file);	
      // $fastq_file_index[$i] = $settings->read_file_to_fastq_file_index($this_read_file);	
      
      // if ($settings->{unmatched_reads})
      if (true)
      {				
        // my $unmatched_file_name = $settings->file_name('unmatched_read_file_name', {'#'=>$this_read_file});
        // $out_unmatched_fastq[$i] = Breseq::Fastq->new(-file => ">$unmatched_file_name");
        
        string this_unmatched_file_name = data_path + "/unmatched." + rf->m_base_name + ".fastq";        
        out_unmatched_fastq.open(this_unmatched_file_name.c_str(), ios::out);
        assert(!out_unmatched_fastq.fail());
      }
      
        
      // my $reference_sam_file_name = $settings->file_name('reference_sam_file_name', {'#'=>$read_file});
      string reference_sam_file_name = reference_sam_path + "/" + rf->m_base_name + ".reference.sam";
      
      // my $reference_tam = Bio::DB::Tam->open($reference_sam_file_name) or die "Could not open $reference_sam_file_name";
      tamFile reference_tam = sam_open(reference_sam_file_name.c_str());
      assert(reference_tam);
      
      // $reference_header = $reference_tam->header_read2($reference_faidx_file_name) or throw("Error reading reference fasta index file: $reference_faidx_file_name");		
      string reference_faidx_file_name = reference_fasta + ".fai";
      bam_header_t *reference_header = sam_header_read2(reference_faidx_file_name.c_str()); 
      assert(reference_header);
      
      // my $junction_tam;
      tamFile junction_tam = NULL;
      bam_header_t *junction_header = NULL;
      
      // if (!$settings->{no_junction_prediction})
      if (junction_prediction)
      {
        // my $junction_sam_file_name = $settings->file_name('candidate_junction_sam_file_name', {'#'=>$read_file});
        string junction_sam_file_name = junction_sam_path + "/" + rf->m_base_name + ".candidate_junction.sam";

        // $junction_tam = Bio::DB::Tam->open($junction_sam_file_name) or die " Could not open junction SAM file\n";
        junction_tam = sam_open(junction_sam_file_name.c_str());
        assert(junction_tam);
        
        //## junction_header = $junction_tam->header_read2($junction_faidx_file_name) or die("Error reading reference fasta index file: $junction_faidx_file_name");				
        string junction_faidx_file_name = junction_fasta + ".fai";
        junction_header = sam_header_read2(junction_faidx_file_name.c_str()); 
        assert(junction_header);
      }
      
        
      // my $reference_al;
      vector<alignment> reference_alignments;
      // my $last_reference_alignment;
      alignment* last_reference_alignment = NULL;
      
      // my $junction_al;
      vector<alignment> junction_alignments;
      // my $last_junction_alignment;
      alignment* last_junction_alignment = NULL;

      //#proceed through all of the alignments
      // if (!$settings->{no_junction_prediction})
      if (junction_prediction) {
        
        // ($junction_al, $last_junction_alignment) 
        //   = Breseq::Shared::tam_next_read_alignments($junction_tam, $junction_header, $last_junction_alignment);		
  
//@JEB: Rewrite to use tam_file
//        junction_alignments = alignment::tam_next_read_alignments(junction_tam, junction_header, last_junction_alignment);
      }
      
      //      
      // ($reference_al, $last_reference_alignment) 
      // = Breseq::Shared::tam_next_read_alignments($reference_tam, $reference_header, $last_reference_alignment);		
      //      
//@JEB: Rewrite to use tam_file      
//      reference_alignments = alignment::tam_next_read_alignments(reference_tam, reference_header, last_reference_alignment);

      
      //###
      //##  Test each read for its matches to the reference and candidate junctions
      //###
      
      // my $f = 0;
      cFastqSequence seq;
      // READ: while (my $seq = $in_fastq[$f]->next_seq)
      
      while (in_fastq.read_sequence(seq))
      {			
        // $reads_processed++;
        reads_processed++;
      
        // last if ($settings->{alignment_read_limit} && ($reads_processed > $settings->{alignment_read_limit}));
        // @JEB - PORT This is unsafe because memory dealloaction has not occurred.
        if ((alignment_read_limit) && (reads_processed > alignment_read_limit)) break;

      
        // print STDERR "    READS:$reads_processed\n" if ($reads_processed % 10000 == 0);
        if (reads_processed % 10000 == 0) fprintf(stderr, "    READS:%u\n", reads_processed);

        //			print "===> Read: $seq->{id}\n" if ($verbose);
        if (verbose) fprintf(stderr, "===> Read: %s\n", seq.m_name.c_str());

        // my $best_junction_score = 0;
        // my $best_reference_score = 0;
        uint32_t best_junction_score = 0;
        uint32_t best_reference_score = 0;
			
        //## Does this read have eligible candidate junction matches?
        // my $this_junction_al = [];
        vector<alignment> this_junction_alignments;

        // print " Before Overlap Reference alignments = "  . (scalar @$this_reference_al) . "\n" if ($verbose);
        if (verbose) fprintf(stderr, " Before Overlap Reference alignments = %u\n", (uint32_t)reference_alignments.size());

        
        // print " Before Overlap Junction alignments = " . (scalar @$this_junction_al) . "\n" if ($verbose);
        if (verbose) fprintf(stderr, " Before Overlap Junction alignments = %u\n", (uint32_t)junction_alignments.size());

        
        // if (($junction_al) && ($junction_al->[0]->qname =~ m/^$seq->{id}/))
        if ((junction_alignments.size() > 0) && (seq.m_name == junction_alignments[0].read_name())) {
          
          //				$this_junction_al = $junction_al;
          //				($junction_al, $last_junction_alignment) 
          //        = Breseq::Shared::tam_next_read_alignments($junction_tam, $junction_header, $last_junction_alignment);
          //       

          this_junction_alignments = junction_alignments;
//          junction_alignments = alignment::tam_next_read_alignments(junction_tam, junction_header, last_junction_alignment);

      
          //###			
          //## Matches to candidate junctions MUST overlap the junction.
          //##
          //## Reduce this list to those that overlap ANY PART of the junction.
          //## Alignments that extend only into the overlap region, are only additional
          //##  evidence for predicted junctions and NOT support for a new junction on 
          //## their own. (They will also match the original reference genome equally well).
          //###
          
          // @$this_junction_al = grep {_alignment_overlaps_junction($junction_info, $_) } @$this_junction_al;
          for (vector<alignment>::iterator it=this_junction_alignments.begin(); it<this_junction_alignments.end(); it++) {
            if (!_alignment_overlaps_junction(junction_info_list, *it)) {
                this_junction_alignments.erase(it);
            }
          }  
        
          // ($best_junction_score, @$this_junction_al) = _eligible_read_alignments($settings, $junction_header, $junction_fai, undef, @$this_junction_al);
//          best_reference_score = _eligible_read_alignments(settings, junction_header, junction_fai, refseqs, this_junction_alignments);

          // print " Best junction score: $best_junction_score\n" if ($verbose);
          if (verbose) fprintf(stderr, " Best junction score: %u\n", best_junction_score);
        }
          
          
        //## Does this read have eligible reference sequence matches?
        // my $this_reference_al = [];
        vector<alignment> this_reference_alignments;
        
        
        // if (($reference_al) && ($reference_al->[0]->qname =~ m/^$seq->{id}/))
        if ((junction_alignments.size() > 0) && (seq.m_name == junction_alignments[0].read_name())) {
          
          // $this_reference_al = $reference_al;
          // ($reference_al, $last_reference_alignment) 
          //    = Breseq::Shared::tam_next_read_alignments($reference_tam, $reference_header, $last_reference_alignment);

          this_reference_alignments = reference_alignments;
//          reference_alignments = alignment::tam_next_read_alignments(reference_tam, reference_header, last_reference_alignment);

          //($best_reference_score, @$this_reference_al) = _eligible_read_alignments($settings, $reference_header, $reference_fai, $ref_seq_info, @$this_reference_al);	
//          best_reference_score = _eligible_read_alignments(settings, reference_header, reference_fai, refseqs, this_reference_alignments);
          
          // print " Best reference score: $best_reference_score\n" if ($verbose);
          if (verbose) fprintf(stderr, " Best reference score: %u\n", best_reference_score);
        }		
      
        //## Nothing to be done if there were no eligible matches to either
        //## Record in the unmatched FASTQ data file
        
        // if ((@$this_junction_al == 0) && (@$this_reference_al == 0))
        if ( (this_junction_alignments.size() == 0) && (this_reference_alignments.size() == 0) ) {
          // $s->{unmatched_reads}++;
          // $out_unmatched_fastq[$f]->write_seq($seq);
          
          summary_info["unmatched_reads"]++;
          out_unmatched_fastq.write_sequence(seq);
        }

        
      
      
            
      
  
//###			
//## Determine if the read has a better match to a candidate junction
//## or to the reference sequence.
//###
//      
//### There are three possible kinds of reads at this point
//##
//## 1: Read has a best match to the reference genome
//## --> Write this match and we are done
//## 2: Read has a best match (or multiple best matches) to junctions
//## --> Keep an item that describes these matches
//## 3: Read has an equivalent match to the reference genome
//##      and goes into the overlap part of a junction condidate
//## --> Keep an item that is not used during scoring
//###
//      
//# Moved to above		
//#			my $best_junction_score = 0;
//#			my $best_reference_score = 0;
//#
//#			if (@$this_junction_al)
//#			{
//#				# This is the length of the match on the query -- NOT the length of the query
//#				$best_junction_score = $this_junction_al->[0]->query->length;
//# THESE SCORES ARE NOT CONSISTENT ACROSS STRANDS DUE TO DIFFERENT INDEL MATCHES
//#				$best_candidate_junction_score = $ca->aux_get("AS");
//#			}
//#			
//#			if (@$this_reference_al)
//#			{
//#				# This is the length of the match on the query -- NOT the length of the query
//#				$best_reference_score = $this_reference_al->[0]->query->length;
//# THESE SCORES ARE NOT CONSISTENT ACROSS STRANDS DUE TO DIFFERENT INDEL MATCHES
//#				$best_reference_score = $ra->aux_get("AS");
//#			}
//			
//# if < 0, then the best match is to the reference
      // my $mapping_quality_difference = $best_junction_score - $best_reference_score;
      int32_t mapping_quality_difference = best_junction_score - best_reference_score;

      // print " Best junction score: $best_junction_score\n" if ($verbose);
        if (verbose) fprintf(stderr, " Best junction score: %i\n", best_junction_score);
        
        // print " Best reference score: $best_reference_score\n" if ($verbose);
        if (verbose) fprintf(stderr, " Best reference score: %i\n", best_reference_score);

        // print " Mapping quality difference: $mapping_quality_difference\n" if ($verbose);
        if (verbose) fprintf(stderr, " Mapping quality difference:  %i\n", best_reference_score);
        
      // print " Final Reference alignments = "  . (scalar @$this_reference_al) . "\n" if ($verbose);
      // print " Final Candidate junction alignments = " . (scalar @$this_junction_al) . "\n" if ($verbose);
//      
//      
//			next READ if ( (scalar @$this_junction_al == 0) && (scalar @$this_reference_al == 0));
//      
//### The best match we found to the reference was no better than the best to the
//### candidate junction. This read potentially supports the candidate junction.
//###
//### ONLY allow EQUAL matches through if they match the overlap only, otherwise
//### you can get predictions of new junctions with all reads supporting them
//### actually mapping perfectly to the reference.
//			
//### best match is to the reference, record in that SAM file.
//			if ($mapping_quality_difference < 0)
//			{
//				print "Best alignment to reference.\n" if ($verbose);
//				_write_reference_matches($settings, $reference_fai, $ref_seq_info, $RREF, $reference_header, $fastq_file_index[$f], @$this_reference_al);
//			}
//			else
//			{	
//				print "Best alignment is to candidate junction. MQD: $mapping_quality_difference\n" if ($verbose);
//				
//				my $item = {
//					reference_alignments => $this_reference_al, 					# reference sequence alignments
//					junction_alignments => $this_junction_al, 						#the BEST candidate junction alignments
//					fastq_file_index => $fastq_file_index[$f],						#index of the fastq file this read came from
//					mapping_quality_difference => $mapping_quality_difference,
//				};
//#print Dumper($item) if ($verbose);
//        
//				
//###
//## Just one best hit to candidate junctions, that is better than every match to the reference
//## #
//				if ((scalar @$this_junction_al == 1) && (scalar $mapping_quality_difference > 0))
//				{
//					my $a = $this_junction_al->[0];
//					my $junction_id = $junction_header->target_name()->[$a->tid];
//#print "$junction_id\n";
//					push @{$matched_junction{$junction_id}}, $item;
//				}	
//###				
//## Multiple equivalent matches to junctions and refernece, ones with most hits later will win these matches
//## If $mapping_quality_difference > 0, then they will count for scoring
//###
//				else
//				{					
//					foreach my $a (@$this_junction_al)
//					{	
//						$item->{degenerate_count} = scalar @$this_junction_al; #mark as degenerate
//						my $junction_id = $junction_header->target_name()->[$a->tid];
//						$degenerate_matches{$junction_id}->{$seq->{id}} = $item;
//					}
//				}
//			}
        
      } // End loop through every read

        
//      
//		} continue {
//			$f++;
//			$f %= scalar @in_fastq;
//		}
        

//      
//## save statistics
//      $summary->{alignment_correction}->{read_file}->{$read_file} = $s;
//    }	
    

  } // End of Read File loop
    
    
//    
//###			
//## Determine which junctions are real, prefer ones with most matches
//###
//    
//    my %accepted_pos_hash_score_distribution;
//    my %observed_pos_hash_score_distribution;
//    
//    my %accepted_min_overlap_score_distribution;
//    my %observed_min_overlap_score_distribution;
//    
//    my @passed_junction_ids = (); 
//    my @rejected_junction_ids = ();
//    my %junction_test_info;	#scoring information about junctions
//		
//###
//## Candidate junctions with unique matches
//###
//    my @sorted_junction_ids = sort {-(scalar @{$matched_junction{$a}} <=> scalar @{$matched_junction{$b}})} keys %matched_junction;
//		
//    print "Degenerate matches before handling ones with unique matches: " . (scalar keys %degenerate_matches) . "\n" if ($verbose);
//		
//    foreach my $key (@sorted_junction_ids)
//    {
//      my ($failed, $has_non_overlap_only) = _test_junction($settings, $summary, $key, \%matched_junction, \%degenerate_matches, \%junction_test_info, $reference_fai, $ref_seq_info, $RREF, $reference_header, $RCJ, $junction_header);
//## save the score in the distribution
//      Breseq::Shared::add_score_to_distribution(\%observed_pos_hash_score_distribution, $junction_test_info{$key}->{pos_hash_score});
//      Breseq::Shared::add_score_to_distribution(\%observed_min_overlap_score_distribution, $junction_test_info{$key}->{min_overlap_score});
//      
//## only count matches that span overlap
//      if (!$has_non_overlap_only)
//      {
//        if (!$failed)
//        {				
//          push @passed_junction_ids, $key; 
//        }
//        else
//        {				
//          push @rejected_junction_ids, $key;
//        }
//      }
//    }
//    
//    print "Degenerate matches after handling ones with unique matches: " . (scalar keys %degenerate_matches) . "\n" if ($verbose);
//    
//###
//## Candidate junctions with ONLY degenerate matches
//##
//###
//    @sorted_junction_ids = sort {-(scalar keys %{$degenerate_matches{$a}} <=> scalar keys %{$degenerate_matches{$b}})} keys %degenerate_matches;
//    while (@sorted_junction_ids)
//    {
//      my $key = shift @sorted_junction_ids;
//      
//      print "Trying degenerate $key...\n" if ($verbose);
//      
//      my ($failed, $has_non_overlap_only) = _test_junction($settings, $summary, $key, \%matched_junction, \%degenerate_matches, \%junction_test_info, $reference_fai, $ref_seq_info, $RREF, $reference_header, $RCJ, $junction_header);
//## save the score in the distribution
//      Breseq::Shared::add_score_to_distribution(\%observed_pos_hash_score_distribution, $junction_test_info{$key}->{pos_hash_score});
//      Breseq::Shared::add_score_to_distribution(\%observed_min_overlap_score_distribution, $junction_test_info{$key}->{min_overlap_score});
//      
//## if it succeeded, then it may have changed the order of the remaining ones by removing some reads...
//      if (!$failed)
//      {
//        @sorted_junction_ids = sort {-(scalar keys %{$degenerate_matches{$a}} <=> scalar keys %{$degenerate_matches{$b}})} keys %degenerate_matches;
//      }
//      
//## only count matches that span overlap
//      if (!$has_non_overlap_only)
//      {
//        if (!$failed)
//        {
//          push @passed_junction_ids, $key;
//        }
//## Failed ones are not kept in the rejected list (but they could be?).
//      }
//    }
//    
//#print successful ones out
//    print "Successful hybrids\n" if ($verbose);
//    
//#Re-sort
//    @passed_junction_ids = sort {-(scalar @{$matched_junction{$a}} <=> scalar @{$matched_junction{$b}})} @passed_junction_ids;
//    @rejected_junction_ids = sort {-(scalar @{$matched_junction{$a}} <=> scalar @{$matched_junction{$b}})} @rejected_junction_ids;
//    
//    foreach my $key (@passed_junction_ids)
//    {	  
//      print "$key\n" if ($verbose);
//      my $item = _junction_to_hybrid_list_item($key, $ref_seq_info, scalar @{$matched_junction{$key}}, $junction_test_info{$key});
//      $gd->add($item);
//      
//## save the score in the distribution
//      Breseq::Shared::add_score_to_distribution(\%accepted_pos_hash_score_distribution, $junction_test_info{$key}->{pos_hash_score});
//      Breseq::Shared::add_score_to_distribution(\%accepted_min_overlap_score_distribution, $junction_test_info{$key}->{min_overlap_score});
//      
//## Create matches from UNIQUE sides of each match to reference genome
//## this fixes, for example appearing to not have any coverage at the origin of a circular DNA fragment
//### Currently, we do not add coverage to redundantly matched sides because we don't know which copy.
//
//next if ( !$settings->{add_split_junction_sides} );
//
//foreach my $match (@{$matched_junction{$key}})
//{
//my $a = $match->{junction_alignments}->[0];
//my $fastq_file_index = $match->{fastq_file_index};
//
//print ">>>>" . $a->qname . "\n" if ($verbose);
//print ">>>>Alignment start-end: " . $a->start . "  " . $a->end . "\n" if ($verbose);
//
//foreach my $side (1, 2)
//{
//my $side_key = 'side_' . $side;
//## Do not count for coverage if it is redundant!!
//
//next if ($item->{"$side_key\_redundant"});
//
//## Write out match corresponding to this part to SAM file
//## By trimming in the candidate junctions sequence, rather than on each half,
//## this is done properly.						
//my $trim = _trim_ambiguous_ends($a, $junction_header, $junction_fai);						
//Breseq::Shared::tam_write_moved_alignment(
//$RREF, 
//$junction_header,
//$a, 
//$fastq_file_index, 
//$item->{"$side_key\_seq_id"}, 
//$item->{"$side_key\_position"}, 
//$item->{"$side_key\_strand"},
//$item->{"$side_key\_overlap"}, 
//$side, 
//$item->{flanking_left}, 
//$item->{alignment_overlap}, 
//$trim
//);		
//}
//}
//}
//
//## Save summary statistics
//$summary->{alignment_correction}->{new_junctions}->{observed_min_overlap_score_distribution} = \%observed_min_overlap_score_distribution;
//$summary->{alignment_correction}->{new_junctions}->{accepted_min_overlap_score_distribution} = \%accepted_min_overlap_score_distribution;
//
//$summary->{alignment_correction}->{new_junctions}->{observed_pos_hash_score_distribution} = \%observed_pos_hash_score_distribution;
//$summary->{alignment_correction}->{new_junctions}->{accepted_pos_hash_score_distribution} = \%accepted_pos_hash_score_distribution;
//
//my @rejected_hybrid_predictions = ();
//foreach my $key (@rejected_junction_ids)
//{
//print "$key\n" if ($verbose);
//my $item = _junction_to_hybrid_list_item($key, $ref_seq_info, scalar @{$matched_junction{$key}}, $junction_test_info{$key});
//GenomeDiff::add_reject_reason($item, "NJ");
//$gd->add($item);
//}
//
//my $jc_genome_diff_file_name = $settings->file_name('jc_genome_diff_file_name');
//$gd->write($jc_genome_diff_file_name);	
//}
//
    
  }
    
//
//=head2 _eligible_alignments
//
//Title   : _eligible_alignments
//Usage   : _eligible_alignments( );
//Function:
//Returns : Best score
//
//=cut
//
  
class mismatch_map_class : public map<alignment*,double> {
public:
  inline bool operator() (vector<alignment>::iterator a1, vector<alignment>::iterator a2) { return (*this)[&(*a1)] < (*this)[&(*a2)]; } 
} mismatch_map;


bool sort_by_mismatches (alignment a1, alignment a2) { return mismatch_map[&a1] < mismatch_map[&a2]; }  
  
uint32_t _eligible_read_alignments(const Settings& settings, const cReferenceSequences& ref_seq_info, vector<alignment>& alignments)
{
	bool verbose = false;

	// These settings are currently not used
	uint32_t minimum_best_score = 0;
	uint32_t minimum_best_score_difference = 0;
	// but the code below works if they are set

	if (alignments.size() <= 0) return false;

	// require a minimum length of the read to be mapped
	for (vector<alignment>::iterator it = alignments.end() - 1; it >= alignments.begin(); it--)
  {
		if (_test_read_alignment_requirements(settings, ref_seq_info, (*it)))
    {
			alignments.erase(it);
    }
  }
	if (alignments.size() == 0) return 0;

  // require read to be mapped! -- @JEB maybe this should be checked sooner?
	for (vector<alignment>::iterator it = alignments.end() - 1; it >= alignments.begin(); it--)
  {
    if (it->unmapped())
    {
      alignments.erase(it);
    }
  }
  if (alignments.size() == 0) return 0;

	// @JEB v1> Unfortunately sometimes better matches don't get better alignment scores!
	// example is 30K88AAXX_LenskiSet2:1:37:1775:92 in RJW1129
	// Here a read with an extra match at the end doesn't get a higher score!!!

	// This sucks, but we need to re-sort matches and check scores ourselves...
	// for now just count mismatches (which include soft padding!)
	//my %mismatch_hash;

  //foreach my $a (@al)
	//{
	//	$mismatch_hash{$a} = Breseq::Shared::alignment_mismatches($a, $reference_header, $reference_fai, $ref_seq_info);
	//}
	//@al = sort { $mismatch_hash{$a} <=> $mismatch_hash{$b} } @al;
  
  //@JEB This method of sorting may be slower than alternatives
  //     Ideally, the scores should be hashes and only references should be sorted.
  mismatch_map.clear();
  for (vector<alignment>::iterator it = alignments.begin(); it < alignments.end(); it++)
  {  
    alignment* ap = &(*it); // we are saving the pointer value as the map key
    mismatch_map[ap] = static_cast<double>(alignment_mismatches(*it, ref_seq_info));
  }
  sort(alignments.begin(), alignments.end(), sort_by_mismatches ); 

	if (verbose)
	{
    for (vector<alignment>::iterator it = alignments.begin(); it < alignments.end(); it++)
    {
      cerr << it->query_start_1() << "-" << it->query_end_1() << " ";
      cerr << it->query_match_length()-mismatch_map[static_cast<alignment*>(&(*it))] << "\n";
    }
	}

	// how many reads share the best score?
  uint32_t last_best(0);
	//my $last_best = 0;
  uint32_t best_score = mismatch_map[&(alignments[0])];
	//my $best_score = $mismatch_hash{$al[0]};
	
  // no scores meet minimum
  if ((minimum_best_score > 0) && (best_score < minimum_best_score)) return 0; 
	//return (0) if (defined $minimum_best_score && ($best_score < $minimum_best_score));

  
	while ((last_best+1 < alignments.size()) && (mismatch_map[&(alignments[last_best+1])] == best_score))
	{
    last_best++;
	}

	//#broken
	//## no scores meet minimum difference between best and next best
	//#if (defined $minimum_best_score_difference && (scalar @al > $last_best+1))
	//#{
	//#	my $second_best_score = $al[$last_best+1]->aux_get("AS");
	//#	return () if ($second_best_score + $minimum_best_score_difference >= $best_score)
	//#}

	//@al = splice @al, 0, $last_best+1;
  
  alignments.resize(last_best+1);

  
	if (verbose)
	{
    cerr << last_best << endl;
    for (vector<alignment>::iterator it = alignments.begin(); it < alignments.end(); it++)
    {
      cerr << it->query_start_1() << "-" << it->query_end_1() << endl;
    }
	}

	// Note that the score we return is higher for matches so we negative this value...
	//return $al[0]->l_qseq-$best_score;
  return best_score;
}
  
//
//
//=head2 _read_alignment_passes_requirements
//
//Title   : _test_read_alignment_requirements
//Usage   : _test_read_alignment_requirements( );
//Function: Tests an individual read alignment for required match lengths
//and number of mismatches
//Returns : 
//
//=cut
//
bool _test_read_alignment_requirements(const Settings& settings, const cReferenceSequences& ref_seq_info, const alignment& a)
{
	bool accept = true;

	bool unmapped = ((a.flag() & BAM_FUNMAP) != 0);
	if (unmapped) return false;

	if (settings.required_match_length > 0)
	{
		//uint32_t* cigar = a.cigar_array(); // cigar array for this alignment
		int32_t alignment_length_on_query = a.cigar_query_length(); //this is the length of the alignment on the read
		if (alignment_length_on_query < settings.required_match_length)
			return false;
	}

	if (settings.require_complete_match)
	{
    if (a.beginning_to_end_match()) return false;
	}
	if (settings.max_read_mismatches > 0)
	{
		int32_t mismatches = alignment_mismatches(a, ref_seq_info);
		if (mismatches > settings.max_read_mismatches) return false;
	}

	return true;
}

//=head2 _alignment_overlaps_junction
//
//Title   : _test_read_alignment_requirements
//Usage   : _test_read_alignment_requirements( );
//Function: Tests an individual read alignment for required match lengths
//and number of mismatches
//Returns : 
//
//=cut
//
//sub _alignment_overlaps_junction
bool _alignment_overlaps_junction(const vector<junction_info>& junction_info_list, alignment a)
{
  //my ($junction_info, $a) = @_;
  //my $this_junction_info = $junction_info->[$a->tid];
  uint32_t tid = a.reference_target_id();
  const junction_info& this_junction_info = junction_info_list[tid];
  //my $overlap = $this_junction_info->{alignment_overlap};
  int32_t overlap = this_junction_info.m_alignment_overlap;
  
  //my $flanking_left = $this_junction_info->{flanking_left};
  //## find the start and end coordinates of the overlap
  //my ($junction_start, $junction_end);
  //
  //$junction_start = $flanking_left + 1;
  uint32_t junction_start = this_junction_info.m_flanking_left + 1;

  //$junction_end = $flanking_left + abs($overlap);
  uint32_t junction_end = this_junction_info.m_flanking_left + abs(overlap);


  //## If it didn't overlap the junction at all
  //## Check coordinates in the "reference" junction sequence
  // return 0 if ($a->start > $junction_end);
  if (a.reference_start_1() > junction_end) return false;
  // return 0 if ($a->end < $junction_start);
  if (a.reference_end_1() < junction_start) return false;
  // return 1;
  return true;
}

  
//    
//    
//    sub _write_reference_matches
//    {
//      my ($settings, $reference_fai, $ref_seq_info, $RREF, $reference_header, $fastq_file_index, @reference_al) = @_;
//      return if (scalar @reference_al == 0); #nice try, no alignments
//      
//      my @trims;
//      foreach my $a (@reference_al)
//      {
//#	push @trims, _trim_ambiguous_ends($a, $reference_header, $reference_fai);	
//        push @trims, _trim_ambiguous_ends($a, $reference_header, $reference_fai, $ref_seq_info); #slightly faster than using fai	
//      }
//      
//      Breseq::Shared::tam_write_read_alignments($RREF, $reference_header, $fastq_file_index, \@reference_al, \@trims);
//    }
//    
//    sub _test_junction
//    {
//      my $verbose = 0;
//      
//      my ($settings, $summary, $junction_seq_id, $matched_junction_ref, $degenerate_matches_ref, $junction_test_info_ref, $reference_fai, $ref_seq_info, $RREF, $reference_header, $RCJ, $candidate_junction_header) = @_;
//      
//      print "Testing $junction_seq_id\n" if ($verbose);
//      
//## variable initialization
//      my $test_info;
//      my $failed = 0;
//      
//## There are two kinds of matches to a candidate junction:
//## (1) Reads that uniquely map to one candidate junction (but any number of times to reference)
//      my @unique_matches = ();
//      @unique_matches = @{$matched_junction_ref->{$junction_seq_id}} if (defined $matched_junction_ref->{$junction_seq_id});
//      
//## (2) Reads that uniquely map equally well to more than one candidate junction (and any number of times to reference)
//      my	@degenerate_matches = ();
//      @degenerate_matches = map { $degenerate_matches_ref->{$junction_seq_id}->{$_} } sort keys %{$degenerate_matches_ref->{$junction_seq_id}}
//      if (defined $degenerate_matches_ref->{$junction_seq_id});
//      
//## FAI target id -- there is no easy way to get this short of loading the entire array and going through them...
//## Debatable about whether we save more string comparisons by doing this here or each time
//      
//## @JEB v1> hash by tid rather than alignment junction names!!
//      my $junction_tid = 0;
//      foreach my $test_junction_seq_id  (@{$candidate_junction_header->target_name})
//      {
//        last if ($junction_seq_id eq $test_junction_seq_id);
//        $junction_tid++;
//      }	
//      die if ($junction_tid >= $candidate_junction_header->n_targets);
//      
//      
//      print "Testing Junction Candidate: $junction_seq_id\n" if ($verbose);
//      print "Unique Matches: " . (scalar @unique_matches) . " Degenerate Matches: " . (scalar @degenerate_matches) . "\n" if ($verbose);
//      
//#### TEST 1: Reads that go a certain number of bp into the nonoverlap sequence on each side of the junction on each strand
//      my $max_left_per_strand = { '0'=> 0, '1'=>0 };
//      my $max_right_per_strand = { '0'=> 0, '1'=>0 };
//      my $max_min_left_per_strand = { '0'=> 0, '1'=>0 };
//      my $max_min_right_per_strand = { '0'=> 0, '1'=>0 };	
//      my $count_per_strand = { '0'=> 0, '1'=>0 };
//      my $total_non_overlap_reads = 0;
//      my $count_per_coord_per_strand;
//      my $min_overlap_score = 0;
//      
//## basic information about the junction
//      my $scj = Breseq::Shared::junction_name_split($junction_seq_id);
//      my $overlap = $scj->{alignment_overlap};
//      my $flanking_left = $scj->{flanking_left};
//      
//## Is there at least one read that isn't overlap only?
//## displaying ones where it doesn't as marginals looks really confusing
//      my $has_non_overlap_only = 1; 
//## @JEB v1> Rename to $has_non_overlap_alignment to give right sense!!
//      
//### We also need to count degenerate matches b/c sometimes ambiguity unfairly penalizes real reads...
//    READ: foreach my $item (@unique_matches, @degenerate_matches)
//      {
//##!!> Matches that don't extend through the overlap region will have the same quality 
//##!!> as a reference match and, therefore, no difference in mapping quality
//##!!> do not count these toward scoring!
//next READ if ($item->{mapping_quality_difference} == 0);
//
//$total_non_overlap_reads++;
//$has_non_overlap_only = 0;
//
//#If there were no degenerate matches, then we could just take the
//#one and only match in the 'junction_alignments' array
//#my $a = $item->{junction_alignments}->[0]; 
//
//## as it is, we must be sure we are looking at the one that matches
//my $a;
//ALIGNMENT: foreach my $candidate_a (@{$item->{junction_alignments}})
//{
//if ($candidate_a->tid == $junction_tid)
//{
//$a = $candidate_a;
//last ALIGNMENT;
//}
//}
//die if (!defined $a);
//
//my $rev_key = ($a->reversed ? 1 : 0);
//$count_per_strand->{$rev_key}++;
//
//# The start coordinate is less likely to be misaligned due to errors
//# than the end coordinate
//my $begin_coord = $rev_key ? $a->end : $a->start;
//$count_per_coord_per_strand->{"$begin_coord-$rev_key"}++;
//
//print "  " . $item->{junction_alignments}->[0]->qname . " $begin_coord-$rev_key\n" if ($verbose);
//
//
//##The left side goes exactly up to the flanking length
//my $this_left = $flanking_left;
//$this_left = $this_left - $a->start+1;
//
//#The right side starts after moving past any overlap (negative or positive)
//my $this_right = $flanking_left+1;
//$this_right += abs($overlap);
//$this_right = $a->end - $this_right+1;
//
//## Update:
//### Score = the minimum unique match length on a side
//### Max_Min = the maximum of the minimum length match sides
//### Max = the maximum match on a side
//### Note that the max and min filtering is really a kind of poor man's KS test
//###   if we implemented that with a certain coverage cutoff it would be a
//###   more principled way of doing things... 
//        if ($this_left < $this_right)
//        {
//          $min_overlap_score += $this_left;
//          $max_min_left_per_strand->{$rev_key} = $this_left if ($max_min_left_per_strand->{$rev_key} < $this_left);
//        }
//        else
//        {
//          $min_overlap_score += $this_right;
//          $max_min_right_per_strand->{$rev_key} = $this_right if ($max_min_right_per_strand->{$rev_key} < $this_right);
//        }
//        
//        $max_left_per_strand->{$rev_key} = $this_left if ($max_left_per_strand->{$rev_key} < $this_left);
//        $max_right_per_strand->{$rev_key} = $this_right if ($max_right_per_strand->{$rev_key} < $this_right);
//        
//      }				
//      
//      my $max_left = ($max_left_per_strand->{'0'} > $max_left_per_strand->{'1'}) ? $max_left_per_strand->{'0'} : $max_left_per_strand->{'1'};
//      my $max_right = ($max_right_per_strand->{'0'} > $max_right_per_strand->{'1'}) ? $max_right_per_strand->{'0'} : $max_right_per_strand->{'1'};
//      
//      my $max_min_left = ($max_min_left_per_strand->{'0'} > $max_min_left_per_strand->{'1'}) ? $max_min_left_per_strand->{'0'} : $max_min_left_per_strand->{'1'};
//      my $max_min_right = ($max_min_right_per_strand->{'0'} > $max_min_right_per_strand->{'1'}) ? $max_min_right_per_strand->{'0'} : $max_min_right_per_strand->{'1'};
//      
//      
//      $test_info = {
//        max_left => $max_left,
//        max_left_minus => $max_left_per_strand->{0},
//        max_left_plus => $max_left_per_strand->{1},
//        max_right => $max_right,
//        max_right_minus => $max_right_per_strand->{0},
//        max_right_plus =>$max_right_per_strand->{1},
//        max_min_right => $max_min_right,
//        max_min_right_minus => $max_min_right_per_strand->{0},
//        max_min_right_plus =>$max_min_right_per_strand->{1},		
//        max_min_left => $max_min_left,
//        max_min_left_minus => $max_min_left_per_strand->{0},
//        max_min_left_plus =>$max_min_left_per_strand->{1},		
//        coverage_minus => $count_per_strand->{0},
//        coverage_plus => $count_per_strand->{1},
//        total_non_overlap_reads => $total_non_overlap_reads,
//        min_overlap_score => $min_overlap_score,
//        pos_hash_score => scalar keys %$count_per_coord_per_strand,
//      };
//      
//      
//## Old way, requiring certain overlap on each side on each strand
//## @JEB !> Best results may be to combine these methods
//      
//## These parameters still need additional testing
//## and, naturally, they have problems with scaling with the
//## total number of reads...
//      
//      my $alignment_on_each_side_cutoff = 14; #16
//      my $alignment_on_each_side_cutoff_per_strand = 9; #13
//      my $alignment_on_each_side_min_cutoff = 3;
//      
//      $failed = 	   ($max_left < $alignment_on_each_side_cutoff) 
//      || ($max_right < $alignment_on_each_side_cutoff)
//      || ($max_left_per_strand->{'0'} < $alignment_on_each_side_cutoff_per_strand) 
//      || ($max_left_per_strand->{'1'} < $alignment_on_each_side_cutoff_per_strand)
//      || ($max_right_per_strand->{'0'} < $alignment_on_each_side_cutoff_per_strand) 
//      || ($max_right_per_strand->{'1'} < $alignment_on_each_side_cutoff_per_strand)
//      || ($max_min_left < $alignment_on_each_side_min_cutoff)
//      || ($max_min_right < $alignment_on_each_side_min_cutoff)
//      ;
//      
//## POS_HASH test
//## New way, but we need to have examined the coverage distribution to calibrate what scores to accept!
//      my $junction_accept_score_cutoff_1 = $summary->{preprocess_coverage}->{$scj->{side_1}->{seq_id}}->{junction_accept_score_cutoff};
//      my $junction_accept_score_cutoff_2 = $summary->{preprocess_coverage}->{$scj->{side_2}->{seq_id}}->{junction_accept_score_cutoff};
//      $failed ||= ( $test_info->{pos_hash_score} < $junction_accept_score_cutoff_1 ) && ( $test_info->{pos_hash_score} < $junction_accept_score_cutoff_2 );
//      
//      print Dumper($test_info) if ($verbose);
//      print ($failed ? "Failed\n" : "Passed\n") if ($verbose);
//      
//###	
//### ADD -- NEED TO CORRECT OVERLAP AND ADJUST NUMBER OF READS SUPPORTING HERE, RATHER THAN LATER
//###
//      
//## DEGENERATE JUNCTION MATCHES
//## ===========================
//## Determine the fate of degenerate reads that map to this junction
//      
//      if (defined $degenerate_matches_ref->{$junction_seq_id})
//      {
//        foreach my $read_name (keys %{$degenerate_matches_ref->{$junction_seq_id}})
//        {	
//          my $degenerate_match = $degenerate_matches_ref->{$junction_seq_id}->{$read_name};
//          my $fastq_file_index = $degenerate_match->{fastq_file_index};
//          my $matched_alignment;
//          
//## Success for this candidate junction... 
//## purge all references to this from the degenerate match hash
//## so that they will not be counted for other junctions
//          if (!$failed)
//          {
//## We need to add this degenerately matched read to the other ones supporting this junction
//            push @{$matched_junction_ref->{$junction_seq_id}}, $degenerate_match;
//            
//# Purge all references to this read from the degenerate match hash
//# so that it cannot be counted for any other junction
//            foreach my $a (@{$degenerate_match->{junction_alignments}})
//            {
//              my $test_junction_seq_id = $candidate_junction_header->target_name()->[$a->tid];
//              $matched_alignment = $a if ($a->tid eq $junction_tid); #this is the one for the current candidate junction
//                
//                delete $degenerate_matches_ref->{$test_junction_seq_id}->{$read_name};
//              delete $degenerate_matches_ref->{$test_junction_seq_id} if (scalar keys %{$degenerate_matches_ref->{$test_junction_seq_id}} == 0);
//            }
//            
//## Keep only the alignment
//##------>		## WE SHOULD ALSO UPDATE THE MAPPING SCORE!
//            my $a;
//          DOMINANT_ALIGNMENT: foreach my $candidate_a (@{$degenerate_match->{junction_alignments}})
//            {
//              if ($candidate_a->tid == $junction_tid)
//              {
//                $a = $candidate_a;
//                last DOMINANT_ALIGNMENT;
//              }
//            }
//            @{$degenerate_match->{junction_alignments}} = ($a);
//          }
//          
//## Failure for this candidate junction...
//## Remove just the degenerate hits to this candidate junction
//## Once all have failed, then we need to add the reference alignments (if any)!
//          else
//          {
//            $degenerate_match->{degenerate_count}--;
//            
//            print "New Degenerate match count: $degenerate_match->{degenerate_count}\n" if ($verbose);
//            
//## This degenerate match missed on all opportunities,
//## we should add it to the reference sequence
//            if ($degenerate_match->{degenerate_count} == 0)
//            {					
//              my $this_reference_al = $degenerate_match->{reference_alignments};
//              _write_reference_matches($settings, $reference_fai, $ref_seq_info, $RREF, $reference_header, $fastq_file_index, @$this_reference_al);					
//            }
//            
//            foreach my $a (@{$degenerate_match->{junction_alignments}})
//            {
//              $matched_alignment = $a if ($a->tid eq $junction_tid); #this is the one for the current candidate junction
//                }	
//          }
//          
//# Write alignment to SAM file for candidate junctions regardless of success...
//          die if (!$matched_alignment);
//          Breseq::Shared::tam_write_read_alignments($RCJ, $candidate_junction_header, $fastq_file_index, [$matched_alignment]) if (!$has_non_overlap_only);
//        }
//        
//## We are completely done with degenerate matches to this junction id.
//## Deleting them here means that we will never go through this loop with them again
//## and is necessary for not doubly writing them.
//        delete $degenerate_matches_ref->{$junction_seq_id};
//      }		
//      
//## UNIQUE JUNCTION MATCHES
//## =======================	
//    READ: foreach my $item (@unique_matches)
//      {
//## Write out the matches to the proper SAM file(s) depending on whether the junction succeeded or failed
//        my $fastq_file_index = $item->{fastq_file_index};
//        
//## ONLY if we failed: write matches to reference sequences
//        if ($failed)
//        {		
//          my $this_reference_al = $item->{reference_alignments};
//          _write_reference_matches($settings, $reference_fai, $ref_seq_info, $RREF, $reference_header, $fastq_file_index, @$this_reference_al);
//        }
//        
//## REGARDLESS of success: write matches to the candidate junction SAM file 
//        Breseq::Shared::tam_write_read_alignments($RCJ, $candidate_junction_header, $fastq_file_index, $item->{junction_alignments})  if (!$has_non_overlap_only);
//      }	
//      
//# Save the test info about this junction.
//      $junction_test_info_ref->{$junction_seq_id} = $test_info;
//      return ($failed, $has_non_overlap_only);
//    }
//    
//    sub _junction_to_hybrid_list_item
//    {
//      my ($key, $ref_seq_info, $total_reads, $test_info) = @_;
//      
//## split the key to an item with information about the junction
//      my $jc = Breseq::Shared::junction_name_split($key);
//      $jc->{key} = $key;	
//      
//## overlap may be adjusted below... this messes up making the alignment
//## 'alignment_overlap' is the original one that applies to the candidate junction BAM file
//## 'overlap' is a version where overlap has been resolved if possible for adding sides of the
//##    alignment
//      $jc->{overlap} = $jc->{alignment_overlap};			
//      $jc->{total_reads} = $total_reads;
//      
//## Redundancy is loaded from the key, but we doubly enforce it when IS elements are involved.
//      
//## Correct for overlapping IS elements
//      
//###
//# IS insertion overlap correction
//#
//# For these the coordinates may have been offset incorrectly initially (because both sides of the junction may look unique)
//# The goal is to offset through positive overlap to get as close as possible to the ends of the IS
//###
//      
//      my $is;
//      foreach my $side_key ('side_1', 'side_2')
//      {
//## Determine IS elements
//## Is it within an IS or near the boundary of an IS in the direction leading up to the junction?			
//        if (my $is = Breseq::ReferenceSequence::find_closest_repeat_region($jc->{$side_key}->{position}, $ref_seq_info->{repeat_lists}->{$jc->{$side_key}->{seq_id}}, 200, $jc->{$side_key}->{strand}))
//        {
//          $jc->{$side_key}->{is}->{name} = $is->{name};
//          $jc->{$side_key}->{is}->{interval} = ($is->{strand} == +1) ? "$is->{start}-$is->{end}" : "$is->{end}-$is->{start}"; 
//          $jc->{$side_key}->{is}->{product} = $is->{product};
//        }
//      }
//      
//## use of $j is historical due to a move of this part of the function from annotate_rearrangements
//      my $j = $jc;
//      sub _add_is_coords_from_interval
//      {
//        my ($c) = @_;
//        return if (!defined $c->{is}); 
//        
//        my ($is_start, $is_end) = split /-/, $c->{is}->{interval};
//        $c->{is}->{strand} = ($is_start < $is_end) ? +1 : -1; 
//        $c->{is}->{start} = ($is_start < $is_end) ? $is_start : $is_end; 
//        $c->{is}->{end} = ($is_start < $is_end) ? $is_end : $is_start;
//      }
//      
//      _add_is_coords_from_interval($j->{side_1});
//      _add_is_coords_from_interval($j->{side_2});
//      
//      $j->{side_1}->{read_side} = -1;
//      $j->{side_2}->{read_side} = +1;
//      
//## Determine which side of the junction is the IS and which is unique
//## these point to the correct initial interval...
//      if (defined $j->{side_1}->{is} && !defined $j->{side_2}->{is})
//      {
//        if (abs($j->{side_1}->{is}->{start} - $j->{side_1}->{position}) <= 20)
//        {
//          $j->{is_side} = $j->{side_1};
//          $j->{is_side}->{is}->{side_key} = 'start';
//        }
//        elsif (abs($j->{side_1}->{is}->{end} - $j->{side_1}->{position}) <= 20 )
//        {
//          $j->{is_side} = $j->{side_1};
//          $j->{is_side}->{is}->{side_key} = 'end';
//        }
//        $j->{unique_side} = $j->{side_2};
//      }
//      
//      if (!defined $j->{side_1}->{is} && defined $j->{side_2}->{is})
//      {
//        if (abs($j->{side_2}->{is}->{start} - $j->{side_2}->{position}) <= 20)
//        {
//          $j->{is_side} = $j->{side_2};
//          $j->{is_side}->{is}->{side_key} = 'start';
//        }
//        elsif (abs($j->{side_2}->{is}->{end} - $j->{side_2}->{position}) <= 20 )
//        {
//          $j->{is_side} = $j->{side_2};
//          $j->{is_side}->{is}->{side_key} = 'end';
//        }
//        $j->{unique_side} = $j->{side_1};
//      }
//      
//## both were IS! -- define as redundant here
//      $j->{side_1}->{redundant} = 1 if (defined $j->{side_1}->{is});
//      $j->{side_2}->{redundant} = 1 if (defined $j->{side_2}->{is});
//      
//#by default, overlap is included on both sides of the junction (possibly changed below)
//      $jc->{side_1}->{overlap} = 0;
//      $jc->{side_2}->{overlap} = 0;		
//      
//## Resolve redundant overlap
//      if ($jc->{overlap} > 0)
//      {
//        $jc->{side_1}->{overlap} = $jc->{overlap};
//        $jc->{side_2}->{overlap} = $jc->{overlap};
//        
//## If there was in IS, resolve overlap so it goes to the edge of the IS element
//        if (defined $j->{is_side})
//        {			
//### first, adjust the repetitive sequence boundary to get as close to the IS as possible
//          my $move_dist = $j->{is_side}->{strand} * ($j->{is_side}->{is}->{$j->{is_side}->{is}->{side_key}} - $j->{is_side}->{position});
//          
//          $move_dist = 0 if ($move_dist < 0);
//          $move_dist = $j->{overlap} if ($move_dist > $j->{overlap});
//          $j->{is_side}->{position} += $j->{is_side}->{strand} * $move_dist;
//          $j->{overlap} -= $move_dist;
//          $j->{is_side}->{overlap} -= $move_dist;
//          
//### second, adjust the unique sequence side with any remaining overlap
//          $j->{unique_side}->{position} += $j->{unique_side}->{strand} * $j->{overlap};	
//          $j->{unique_side}->{overlap} -= $j->{overlap};
//          
//          $j->{overlap} = 0;
//        }
//        
//### If there is no IS element and 
//##    (1) both sides are unique
//## OR (2) only the second side is redundant,
//## OR (3) both sides are redundant
//### then give overlap to first side. 
//### This gives proper support for junctions.
//### and ensures we don't count this coverage twice.
//elsif ((!$jc->{side_1}->{redundant}) || ($jc->{side_1}->{redundant} && $jc->{side_2}->{redundant}) )
//{
//my $strand_direction = ($jc->{side_2}->{strand} > 0) ? +1 : -1;
//
//$jc->{side_2}->{position} += $jc->{overlap} * $strand_direction;
//$jc->{side_2}->{overlap} = 0;
//$jc->{overlap} = 0;			
//}
//else  ## side_1 was redundant, give overlap to side_2
//{
//my $strand_direction = ($jc->{side_1}->{strand} > 0) ? -1 : +1;
//$jc->{side_1}->{position} += $jc->{overlap} * $strand_direction;
//$jc->{side_1}->{overlap} = 0;
//$jc->{overlap} = 0;
//}
//
//## If both sides were redundant, no adjustment because we are not going to count coverage
//}
//
//##flatten things to only what we want to keep
//my $item = {
//type => 'JC',
//
//side_1_seq_id => $jc->{side_1}->{seq_id},
//side_1_position => $jc->{side_1}->{position},
//side_1_redundant => $jc->{side_1}->{redundant},
//side_1_strand => $jc->{side_1}->{strand},
//side_1_overlap => $jc->{side_1}->{overlap},
//
//side_2_seq_id => $jc->{side_2}->{seq_id},
//side_2_position => $jc->{side_2}->{position},
//side_2_redundant => $jc->{side_2}->{redundant},
//side_2_strand => $jc->{side_2}->{strand},
//side_2_overlap => $jc->{side_2}->{overlap},
//
//key => $jc->{key},
//alignment_overlap => $jc->{alignment_overlap},
//overlap => $jc->{overlap},
//total_reads => $jc->{total_reads},
//flanking_left => $jc->{flanking_left},
//flanking_right => $jc->{flanking_right},
//
//unique_read_sequence => $jc->{unique_read_sequence},
//};	
//
//## may want to take only selected of these fields...
//foreach my $key (keys %$test_info)
//{		
//$item->{$key} = $test_info->{$key};
//}	
//
//### Note: Other adjustments to overlap can happen at the later annotation stage
//### and they will not affect coverage for calling deletions or mutations
//### because they will be in REDUNDANTLY matched sides of junctions
//return $item;
//}
//
//sub _trim_ambiguous_ends
//{	
//my $verbose = 0;
//my ($a, $header, $fai, $ref_seq_info) = @_;
//
//#which reference sequence?
//my $seq_id = $header->target_name->[$a->tid];
//my $ref_seq_length = $header->target_len->[$a->tid];
//
//###############################################
//## NEW version using preacalculated trim files	
//###############################################
//if ((defined $ref_seq_info) && (defined $ref_seq_info->{trims}))
//{			
//### To accelerate this code further we may want to have trims
//### in array by target ID so translation to name and hash lookup
//### can be skipped.
//
//my ($left_trim, $right_trim);
//
//## These give identical results and seem to give identical speeds. 
//## First one is more scalable if we want to allow trims >255		
//$left_trim = vec($ref_seq_info->{trims}->{$seq_id}, $a->start-1, 8);
//$right_trim = vec($ref_seq_info->{trims}->{$seq_id}, $a->end-1+$ref_seq_length, 8);
//#$left_trim = unpack "C", substr($ref_seq_info->{trims}->{$seq_id}, $a->start-1, 1);
//#$right_trim = unpack "C", substr($ref_seq_info->{trims}->{$seq_id}, $a->end-1+$ref_seq_length, 1);
//
//$left_trim += $a->query->start - 1;
//$right_trim += $a->l_qseq - $a->query->end;
//
//#print STDERR $a->query->name . "\n";
//#print STDERR "start: " . ($a->start) . " end: " . ($a->end) . "\n";		
//#print STDERR "left: " . $left_trim . " right: " . $right_trim . "\n";
//
//return ( {'L'=>$left_trim, 'R'=>$right_trim} );
//}
//
//###############################################
//## OLD version using Perl string comparisons	
//###############################################
//
//# Has two keys: 'left' and 'right' which are how far to inset in REFERENCE coords.
//my $trims;
//
//## using $fai is more compatible, and must currently be used for junctions
//## using $ref_seq_info is slightly quicker, and currently used for the reference sequence
//my $ref_strings;
//if (defined $ref_seq_info)
//{		
//$ref_strings = $ref_seq_info->{ref_strings};	
//}
//
//#create sequence snippets that we need to pay attention to ends of sequence
//my $expand_by = 18; #36
//my $expand_left = ($a->start-1 < $expand_by) ? $a->start-1 : $expand_by;
//my $expand_right = ($ref_seq_length - $a->end < $expand_by) ? $ref_seq_length-$a->end : $expand_by;
//
//my $expanded_ref_string = '';	
//if (defined $ref_strings)
//{
//$expanded_ref_string = substr $ref_strings->{$seq_id}, $a->start-$expand_left-1, ($a->end+$expand_right) - ($a->start-$expand_left) + 1;
//}
//## >>> transition to not using ref_seq_info
//else
//{
//my $expanded_ref_range = $seq_id . ':' . ($a->start-$expand_left) . '-' . ($a->end+$expand_right);
//$expanded_ref_string = $fai->fetch($expanded_ref_range);		
//}
//
//my $ref_string;
//if (defined $ref_strings)
//{
//$ref_string = substr $ref_strings->{$seq_id}, $a->start-1, $a->end - $a->start + 1;
//}	
//## >>> transition to not using ref_seq_info	
//else
//{
//my $ref_range = $seq_id . ':' . $a->start . '..' . $a->end;
//$ref_string = $fai->fetch( $seq_id . ':' . $a->start . '-' . $a->end );
//}
//
//my ($q_start, $q_end) = Breseq::Shared::alignment_query_start_end($a);
//my $q_length = $a->l_qseq;
//my $qry_string = substr $a->qseq, $q_start-1, $q_end - $q_start + 1;
//my $full_qry_string = $a->qseq;
//
//#take maximum of inset for query and reference
//#	my ($left_ref_inset, $right_ref_inset) = (0,0); ##TESTING	
//my ($left_ref_inset, $right_ref_inset) = _ambiguous_end_offsets_from_sequence($ref_string);
//
//#add UNALIGNED bases at te end of reads
//$left_ref_inset += $q_start - 1;
//$right_ref_inset += $q_length - $q_end;
//
//## save a little time if qry and ref sequences are identical.
//my ($left_qry_inset, $right_qry_inset) = (0,0);
//if ($ref_string ne $qry_string)
//{
//($left_qry_inset, $right_qry_inset) = _ambiguous_end_offsets_from_sequence($qry_string);
//#add UNALIGNED bases at the end of reads
//$left_qry_inset += $q_start - 1;
//$right_qry_inset += $q_length - $q_end;
//}
//
//
//#	my ($left_full_qry_inset, $right_full_qry_inset) = (0,0); ##TESTING
//#	my ($left_ref_expanded_inset, $right_ref_expanded_inset) = (0,0); ##TESTING
//
//my ($left_full_qry_inset, $right_full_qry_inset) = _ambiguous_end_offsets_from_sequence($full_qry_string);
//my ($left_ref_expanded_inset, $right_ref_expanded_inset) 
//= _ambiguous_end_offsets_from_expanded_sequence($expand_left, $expand_right, $expanded_ref_string);
//
//if ($verbose)
//{
//print "Whole Read: $ref_string\n";
//print "Qry Start, End: $q_start, $q_end\n";	
//print "Ref: $ref_string\n";
//print "Ref insets: $left_ref_inset, $right_ref_inset\n";
//print "Qry: $qry_string\n";
//print "Qry insets: $left_qry_inset, $right_qry_inset\n";
//print "Full Qry: $full_qry_string\n";
//print "Full Qry insets: $left_full_qry_inset, $right_full_qry_inset\n";
//print "Expanded: $expanded_ref_string\n";
//print "Expand: $expand_left, $expand_right\n";
//print "Expanded Ref insets: $left_ref_expanded_inset, $right_ref_expanded_inset\n";		
//}
//
//$left_qry_inset = ($left_qry_inset > $left_full_qry_inset) ? $left_qry_inset : $left_full_qry_inset;
//$right_qry_inset = ($right_qry_inset > $right_full_qry_inset) ? $right_qry_inset : $right_full_qry_inset;
//
//$left_ref_inset = ($left_ref_inset > $left_ref_expanded_inset) ? $left_ref_inset : $left_ref_expanded_inset;
//$right_ref_inset = ($right_ref_inset > $right_ref_expanded_inset) ? $right_ref_inset : $right_ref_expanded_inset;
//
//##
//# Correct insets in the ref sequence to read coordinates, which must pay attention to gaps
//##
//
//#if no gaps, then just need to take the greater one
//my $left_ref_count = ($left_ref_inset > $left_qry_inset) ? $left_ref_inset : $left_qry_inset;
//my $left_qry_count = $left_ref_count;
//my $right_ref_count = ($right_ref_inset > $right_qry_inset) ? $right_ref_inset : $right_qry_inset;
//my $right_qry_count = $right_ref_count;
//
//if ($verbose)
//{
//print "Ref insets: $left_ref_inset, $right_ref_inset\n";
//print "Full Qry insets: $left_full_qry_inset, $right_full_qry_inset\n";
//print "Expanded Ref insets: $left_ref_expanded_inset, $right_ref_expanded_inset\n";	
//print "Qry Count: $left_qry_count, $right_qry_count\n";	
//}
//
//return ( {'L'=>$left_qry_count, 'R'=>$right_qry_count} );
//}
//
//sub _ambiguous_end_offsets_from_expanded_sequence
//{
//my $verbose = 0;
//my ($expand_left, $expand_right, $ref_string) = @_;
//
//my $left_inset = 0;
//my $right_inset = 0;
//
//{ #left side
//
//#maximum size to check is $expand_by
//my $test_left_inset = 0;
//while ($test_left_inset < $expand_left)
//{			
//my $found_left_inset = $test_left_inset;
//my $test_length = $expand_left-$test_left_inset;
//
//my $match_found = 0;
//my $test_end_string;
//while ( ($test_length > 0) && !$match_found)
//{	
//$test_end_string = substr $ref_string, $found_left_inset, $test_length;
//#force removal of end nucleotides if no higher order repeats found
//$found_left_inset-- if ($test_length == 1);			
//my $test_interior_string = substr $ref_string, $found_left_inset+$test_length, $test_length;
//
//while ($test_end_string eq $test_interior_string)
//{
//$match_found = 1;
//$found_left_inset += $test_length;
//$test_interior_string = substr $ref_string, $found_left_inset+$test_length, $test_length;
//}
//
//print "Complete:\n" if ($verbose);
//print "left_inset $test_left_inset :: test_length $test_length\n" if ($verbose);
//print "test_end: $test_end_string\ntest_interior: $test_interior_string\n" if ($verbose);
//
//$found_left_inset += $test_length if ($match_found);
//$test_length--;
//} 
//
//#test partial matches (that continue part of repeat further)
//#note: already starts at one less than actual repeat size
//my $test_partial_size = $test_length; 
//while ($test_partial_size > 0)
//{
//my $test_partial_end_string = substr $test_end_string, 0, $test_partial_size;
//my $test_interior_string = substr $ref_string, $found_left_inset, $test_partial_size;
//if ($test_partial_end_string eq $test_interior_string)
//{
//$found_left_inset += $test_partial_size;
//last;
//}
//
//print "Partial:\n" if ($verbose);
//print "test_partial_size $test_partial_size :: test_length $test_length\n" if ($verbose);
//print "test_end: $test_partial_end_string\ntest_interior: $test_interior_string\n" if ($verbose);
//
//$test_partial_size--;
//}
//
//$left_inset = $found_left_inset-$expand_left if ($found_left_inset-$expand_left > $left_inset);
//$test_left_inset++;
//}
//}
//
//{ #right side
//
//#maximum size to check is $expand_by
//my $test_right_inset = 0;
//while ($test_right_inset < $expand_right)
//{
//my $found_right_inset = $test_right_inset;
//my $test_length = $expand_right-$test_right_inset;
//
//my $match_found = 0;
//my $test_end_string;
//while ( ($test_length > 0) && !$match_found)
//{	
//$test_end_string = substr $ref_string, $found_right_inset, $test_length;
//#force removal of end nucleotides if no higher order repeats found
//$found_right_inset-- if ($test_length == 1);			
//my $test_interior_string = substr $ref_string, $found_right_inset+$test_length, $test_length;
//
//while ($test_end_string eq $test_interior_string)
//{
//$match_found = 1;
//$found_right_inset += $test_length;
//$test_interior_string = substr $ref_string, $found_right_inset+$test_length, $test_length;
//}
//$found_right_inset += $test_length if ($match_found);
//$test_length--;
//} 
//
//#test partial matches (that continue part of repeat further)
//#note: already starts at one less than actual repeat size
//my $test_partial_size = $test_length; 
//while ($test_partial_size > 0)
//{
//my $test_partial_end_string = substr $test_end_string, 0, $test_partial_size;
//my $test_interior_string = substr $ref_string, $found_right_inset, $test_partial_size;
//if ($test_partial_end_string eq $test_interior_string)
//{
//$found_right_inset += $test_partial_size;
//last;
//}
//$test_partial_size--;
//}
//
//$right_inset = $found_right_inset-$expand_right if ($found_right_inset-$expand_right > $right_inset);
//$test_right_inset++;
//}
//}
//
//
//
//return ($left_inset, $right_inset);
//}
//
//
//sub _ambiguous_end_offsets_from_sequence
//{
//my ($ref_string) = @_;
//
//#test longest substrings
//my $left_inset = 0;
//my $right_inset = 0;
//
//{ #left side
//my $test_length = int((length $ref_string) / 2);
//my $match_found = 0;
//my $test_end_string;
//while ( ($test_length > 0) && !$match_found)
//{	
//$test_end_string = substr $ref_string, $left_inset, $test_length;
//#force removal of end nucleotides if no higher order repeats found
//$left_inset-- if ($test_length == 1);			
//my $test_interior_string = substr $ref_string, $left_inset+$test_length, $test_length;
//
//while ($test_end_string eq $test_interior_string)
//{
//$match_found = 1;
//$left_inset += $test_length;
//$test_interior_string = substr $ref_string, $left_inset+$test_length, $test_length;
//}
//$left_inset += $test_length if ($match_found);
//$test_length--;
//} 
//
//#test partial matches (that continue part of repeat further)
//#note: already starts at one less than actual repeat size
//my $test_partial_size = $test_length; 
//while ($test_partial_size > 0)
//{
//my $test_partial_end_string = substr $test_end_string, 0, $test_partial_size;
//my $test_interior_string = substr $ref_string, $left_inset, $test_partial_size;
//if ($test_partial_end_string eq $test_interior_string)
//{
//$left_inset += $test_partial_size;
//last;
//}
//$test_partial_size--;
//}
//
//}
//
//{ #right side
//my $test_length = int((length $ref_string) / 2);
//my $match_found = 0;
//my $test_end_string;
//while ( ($test_length > 0) && !$match_found)
//{	
//$test_end_string = substr $ref_string, -$right_inset-$test_length, $test_length;
//#force removal of end nucleotides if no higher order repeats found
//$right_inset-- if ($test_length == 1);			
//my $test_interior_string = substr $ref_string, -$right_inset-2 * $test_length, $test_length;
//
//while ($test_end_string eq $test_interior_string)
//{
//$match_found = 1;
//$right_inset += $test_length;
//$test_interior_string = substr $ref_string, -$right_inset-2 * $test_length, $test_length;
//}
//$right_inset += $test_length if ($match_found);
//
//$test_length--;
//} 
//
//#test partial matches (that continue part of repeat further)
//#note: already starts at one less than actual repeat size
//my $test_partial_size = $test_length; 
//while ($test_partial_size > 0)
//{
//my $test_partial_end_string = substr $test_end_string, -$test_partial_size, $test_partial_size;
//my $test_interior_string = substr $ref_string, -$right_inset-$test_partial_size, $test_partial_size;
//if ($test_partial_end_string eq $test_interior_string)
//{
//$right_inset += $test_partial_size;
//last;
//}
//$test_partial_size--;
//}
//}
//
//return ($left_inset, $right_inset);
//}
//
//
//
//return 1;

  
} // namespace breseq

