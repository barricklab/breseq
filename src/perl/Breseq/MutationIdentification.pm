###
# Pod Documentation
###

=head1 NAME

Breseq::Fastq.pm

=head1 SYNOPSIS

Module for reading and writing fastq files more rapidly than BioPerl.

=head1 DESCRIPTION

=head1 AUTHOR

Jeffrey Barrick
<jbarrick@msu.edu>

=head1 COPYRIGHT

Copyright 2009.  All rights reserved.

=cut

###
# End Pod Documentation
###

use strict;

use Math::CDF;

use Bio::Root::Root;
use Bio::DB::Sam;

use Breseq::GenomeDiff;
use Breseq::Shared;


package Breseq::MutationIdentification;
use vars qw(@ISA);
@ISA = qw( Bio::Root::Root );

use Data::Dumper;

use Breseq::GenomeDiff;

our @base_list = ('A', 'T', 'C', 'G', '.');

#does both within-read SNPs+indels and missing coverage deletions
sub identify_mutations
{
	our ($settings, $summary, $ref_seq_info, $error_rates) = @_;
		
	my $reference_fasta_file_name = $settings->file_name('reference_fasta_file_name');
	my $reference_bam_file_name = $settings->file_name('reference_bam_file_name');
	my $bam = Bio::DB::Sam->new(-fasta => $reference_fasta_file_name, -bam => $reference_bam_file_name);
	##my @seq_ids = $bam->seq_ids;
	my @seq_ids = @{$ref_seq_info->{seq_ids}};
	our $gd = Breseq::GenomeDiff->new();
	
	## some local variable lookups for convenience
	my $total_ref_length = 0;
	foreach my $seq_id (@seq_ids)
	{
		$total_ref_length+= $bam->length($seq_id);
	}
	
	my $log10_ref_length = log($total_ref_length) / log(10);	
	our $s;

	### 
    ## Data Preparation
	###

	my ($log10_correct_rates, $log10_error_rates) = Breseq::ErrorCalibration::log10_error_rates($error_rates);
   
	
	REFERENCE: foreach our $seq_id (@seq_ids)
	{				
		my $this_predicted_mutation_file_name = $settings->file_name('predicted_mutation_file_name', {'@'=>$seq_id});
		my $this_tiled_coverage_tab_file_name = $settings->file_name('tiled_coverage_text_file_name', {'@'=>$seq_id});	
		
		###
		##  Already done with this reference sequence.
		###
#		if (-e $this_mutation_identification_done_file_name)
#		{
#			print " REFERENCE: $seq_id :: Already Complete. \n";			
#			next REFERENCE;
#		}
		
		###
		##  Handle predictions for this reference sequence.
		###
		my $snps_all_tab_file_name = $settings->file_name('complete_mutations_text_file_name', {'@'=>$seq_id}); 
		my $coverage_tab_file_name = $settings->file_name('complete_coverage_text_file_name', {'@'=>$seq_id}); 
		
		my $sequence_length = $bam->length($seq_id);
		
		### deletion settings may differ between reference sequences
		our $deletion_seed_cutoff = 0;
		our $deletion_propagation_cutoff = $settings->{unique_coverage}->{$seq_id}->{deletion_coverage_propagation_cutoff};
		
		print STDERR "  REFERENCE: $seq_id\n";
		print STDERR "  LENGTH: $sequence_length\n";
	
		## open detailed output files
		open MUT, ">$snps_all_tab_file_name" if (defined $snps_all_tab_file_name);
		open COV, ">$coverage_tab_file_name" if (defined $coverage_tab_file_name);
		print COV join("\t", 'unique_top_cov', 'unique_bot_cov', 'redundant_top_cov', 'redundant_bot_cov', 'raw_redundant_top_cov', 'raw_redundant_bot_cov', 'e_value', 'position') . "\n";
	
		### DELETION PREDICTION: variables to keep track of during pileup
		our $last_deletion_start_position = undef;
		our $this_deletion_reaches_seed_value = 0;
		our $left_side_coverage_item = undef;
		our $left_inside_coverage_item = undef;
		our $last_position_coverage = undef;

		### COPY NUMBER VARIATION: variables to keep track of during pileup
		our $cnv_tile_size = 100;
		our $cnv_tile_num = 0;
		our $cnv_encountered_redundant = 0;		
		our $cnv_cumulative_coverage = 0; 
		our $cnv_base_counts = { 'A'=> 0, 'C'=>0, 'G'=>0, 'T'=> 0};
		
		### UNKNOWN INTERVALS: variable to keep track of during pileup
		our $last_start_unknown_interval;

		my $cnv_coverage_tab_file_name = $settings->file_name('cnv_coverage_tab_file_name', {'@'=>$seq_id}); 
		open CNV_COV, ">$cnv_coverage_tab_file_name" if (defined $cnv_coverage_tab_file_name);
		print CNV_COV "pos\tcov\tred\tgc\n";

		my $last_insert_count_was_accepted = 1;
		our $last_position_coverage_printed = 0;

		##
		#  Utility function used at each pileup iteration and at the end
		##
		sub _check_deletion_completion 
		{
			my ($pos, $this_position_coverage, $e_value_call) = @_;
			
			## when called at the end of a fragment, the position is 1+fragment length
			## and 
			
			#we need to fill in positions with NO reads			
			foreach (my $i=$last_position_coverage_printed+1; $i<$pos; $i++)
			{
				my $zero_coverage = { 
					'unique' => {'1'=>0, '-1'=>0, 'total' => 0 },
					'redundant' => {'1'=>0, '-1'=>0, 'total' => 0 }
				};
				
				if (!defined $last_deletion_start_position)
				{
					$last_deletion_start_position = $i;				
					$left_side_coverage_item = $zero_coverage;
					$left_inside_coverage_item = $last_position_coverage;
				}
				$this_deletion_reaches_seed_value = 1;
				
				print COV join("\t", 0, 0, 0, 0, 0, 0, 'NA', $i) . "\n";
			}
			$last_position_coverage_printed = $pos;
			
			## called with an undef $this_position_coverage at the end of the genome
			if ((defined $this_position_coverage) && (defined $e_value_call))
			{				
				my $tu = $this_position_coverage->{unique};
				my $tr = $this_position_coverage->{redundant};
				my $trr = $this_position_coverage->{raw_redundant};
						
				#print this information
				print COV join("\t", $tu->{-1}, $tu->{1}, $tr->{-1}, $tr->{1}, $trr->{-1}, $trr->{1}, $e_value_call, $pos) . "\n";
			
				#start a new possible deletion if we fall below the propagation cutoff
				if ($this_position_coverage->{unique}->{total} <= $deletion_propagation_cutoff)
				{	
					if (!defined $last_deletion_start_position)
					{
						$last_deletion_start_position = $pos;				
						$left_side_coverage_item = $this_position_coverage;
						$left_inside_coverage_item = $last_position_coverage;
					}
				}
				
				##keep track of whether we've encountered the seed value
				if ($this_position_coverage->{total} <= $deletion_seed_cutoff)
				{
					$this_deletion_reaches_seed_value = 1;
				}
			}
			
			##if we are above the propagation cutoff then record the current deletion
			if ( (defined $last_deletion_start_position) && ((!defined $this_position_coverage) 
				|| ($this_position_coverage->{unique}->{total} > $deletion_propagation_cutoff)) )
			{

				if ($this_deletion_reaches_seed_value)
				{				
					### for the end of the genome....
					if (!defined $this_position_coverage)
					{
						$this_position_coverage = {
							unique => {'1'=>'NA', '-1'=>'NA', 'total' => 'NA' },
							raw_redundant => {'1'=>'NA', '-1'=>'NA', 'total' => 'NA' },
							redundant => {'1'=>'NA', '-1'=>'NA', 'total' => 'NA' }
						};
					}
					
					my $del = {
						type => 'MC',
						seq_id => $seq_id,
						start => $last_deletion_start_position,
						end => $pos-1,
						start_range => 0,
						end_range => 0,
			#			size => ($pos-1) - $last_deletion_start_position + 1, #end - start + 1
						left_unique_cov => $left_side_coverage_item->{unique}->{total},
						left_inside_unique_cov => $left_inside_coverage_item->{unique}->{total},
						right_inside_unique_cov => $last_position_coverage->{unique}->{total},
						right_unique_cov => $this_position_coverage->{unique}->{total},
					};
					$del->{left_inside_unique_cov} = 'NA' if (!defined $del->{left_inside_unique_cov});
					$del->{right_inside_unique_cov} = 'NA' if (!defined $del->{right_inside_unique_cov});
					
					$del->{left_unique_cov} = 'NA' if (!defined $del->{left_unique_cov});
					$del->{right_unique_cov} = 'NA' if (!defined $del->{right_unique_cov});
					
					$gd->add($del);					
				}

				#reset the search
				$this_deletion_reaches_seed_value = 0;
				undef $last_deletion_start_position;
			}

			$last_position_coverage = $this_position_coverage;
		} #end _check_deletion_completion

		##
		#  Utility function used at each pileup iteration and at the end
		##
		sub _update_copy_number_variation
		{
			## Do not assume that we will be called for every position
			my ($pos, $this_position_coverage, $ref_base) = @_;
				
			my $this_tile_num = int(($pos+1)/$cnv_tile_size);
			
			## make a prediction for this window
			if ($this_tile_num != $cnv_tile_num)
			{
				## these all give less variation than is actually present
	
				my $cn = $cnv_cumulative_coverage / ($summary->{unique_coverage}->{$seq_id}->{average} * $cnv_tile_size);
				my $start_pos = $cnv_tile_num * $cnv_tile_size + 1;
				my $end_pos = ($cnv_tile_num+1) * $cnv_tile_size;			
				my $GC = 'NA';	
				my $total = ($cnv_base_counts->{C} + $cnv_base_counts->{G} + $cnv_base_counts->{A} + $cnv_base_counts->{T});
				$GC = ($cnv_base_counts->{C} + $cnv_base_counts->{G}) / $total if ($total > 0);
				
				#$print "$cnv_tile_num $start_pos\-$end_pos COV: $cn\n";
	
				print CNV_COV "$cnv_tile_num\t$cnv_cumulative_coverage\t$cnv_encountered_redundant\t$GC\n";
				
				##reset values
				##use Central Limit Theorem?
				# my $var = $summary->{unique_coverage}->{$seq_id}->{variance} * $cnv_tile_size;
				# my $mean = $summary->{unique_coverage}->{$seq_id}->{average} * $cnv_tile_size;
				# my $z = ($cnv_cumulative_coverage-$mean)/sqrt($var);
				# 
				# ## calculate one tailed-probability of observed counts
				# my $p_value = Math::CDF::pnorm($z);
				# print "$cnv_tile_num $start_pos\-$end_pos COV: $cnv_cumulative_coverage MEAN: $mean VAR: $var PROB: $p_value Z: $z\n";
				# 
				# ## use Sum of Negative Binomials, which will have the same prob
				# my $nb_mean = $summary->{unique_coverage}->{$seq_id}->{average} * $cnv_tile_size;
				# my $nb_prob = $summary->{unique_coverage}->{$seq_id}->{nbinom_prob_parameter};
				# my $nb_size = $summary->{unique_coverage}->{$seq_id}->{nbinom_size_parameter} * $cnv_tile_size;
				# 
				# my $nb_p_value = Math::CDF::pnbinom($cnv_cumulative_coverage, $nb_size, $nb_prob);				
				# print "$cnv_tile_num $start_pos\-$end_pos COV: $cnv_cumulative_coverage MEAN: $nb_mean P: $nb_prob NBPROB: $nb_p_value\n";

				$cnv_tile_num = $this_tile_num;
				$cnv_cumulative_coverage = 0;
				$cnv_encountered_redundant = 0;
				$cnv_base_counts = { 'A'=> 0, 'C'=>0, 'G'=>0, 'T'=> 0};
			}

			$cnv_encountered_redundant = 1 if ($this_position_coverage->{redundant}->{total} >= 1);
			$cnv_cumulative_coverage += $this_position_coverage->{unique}->{total};
			$cnv_base_counts->{$ref_base}++;
			
		}
		

		sub _update_unknown_intervals
		{
			my ($seq_id, $pos, $base_predicted, $this_position_unique_only_coverage) = @_;
		
			if (!$base_predicted)
			{
				$s->{coverage}->{unique_uncalled}++ if (($this_position_unique_only_coverage));
				if (!defined $last_start_unknown_interval)
				{
					$last_start_unknown_interval = $pos;
				}
			}	
			else
			{
				$s->{coverage}->{unique_called}++ if (($this_position_unique_only_coverage));

				#end interval where we were unable to call mutations
				if (defined $last_start_unknown_interval)
				{
					my $new_interval = { 'type'=>'UN', 'start'=> $last_start_unknown_interval, 'end'=> $pos-1, 'seq_id' => $seq_id };
					$gd->add($new_interval);
					undef $last_start_unknown_interval;
				#	print Dumper($new_interval); ##DEBUG
				}
			}
		}
		
		###
		##  The pileup function takes care of the bulk of the processing
		###
		my $pileup_function = sub {
			my ($seqid,$pos,$pileup) = @_;
	       	 
			print STDERR "    POSITION:$pos\n" if ($pos % 10000 == 0);			

			my $insert_count = 0;
			my $next_insert_count_exists = 1;
			INSERT_COUNT: while ($next_insert_count_exists)
			{				
				#print STDERR "$pos $insert_count\n" if ($insert_count > 0);

				$next_insert_count_exists = 0;

				my $ref_base = ($insert_count) ? '.' : $bam->segment($seqid,$pos,$pos)->dna;

				# zero out the info about this position
				my $pos_info;
				foreach my $base (@base_list)
				{
					$pos_info->{$base}->{unique_cov}->{1} = 0;
					$pos_info->{$base}->{unique_cov}->{-1} = 0;
					$pos_info->{$base}->{unique_trimmed_cov}->{1} = 0;
					$pos_info->{$base}->{unique_trimmed_cov}->{-1} = 0;
					$pos_info->{$base}->{mutation_cov}->{1} = 0;
					$pos_info->{$base}->{mutation_cov}->{-1} = 0;
				}
			
				## keep track of coverage for deletion prediction
				my $this_position_coverage;
				$this_position_coverage->{unique} = {'-1' => 0, '1' => 0, 'total' => 0};
				$this_position_coverage->{redundant} = {'-1' => 0, '1' => 0, 'total' => 0};
				$this_position_coverage->{raw_redundant} = {'-1' => 0, '1' => 0, 'total' => 0};
				$this_position_coverage->{total} = 0;
				my $this_position_unique_only_coverage = 1;

				my $line = '';

				## calculate the chance of observing alignment given each possible base was 100% of the population
				my $pr_base_hash;
				my $pr_not_base_hash;
				
				## polymorphism prediction
				my $pdata;

				ALIGNMENT: foreach my $p (@$pileup) 
				{
					my $a = $p->alignment;
									
					## This setup gives expected behavior from indel!
					my $indel = $p->indel;       ## insertions relative to the reference have the 
												 ## number of this inserted base
					$indel = 0 if ($indel < 0);  ## substitute such that
					$indel = -1 if ($p->is_del); ## deletions relative to reference have -1 as indel

					my $redundancy = $a->aux_get('X1');
					my $fastq_file_index = $a->aux_get('X2');
					my $strand = $a->reversed ? -1 : +1;

					## Handle trimming
					## Note that trimming INCLUDES the unaligned bases on each end
					my $trimmed = 0;
					my $trim_left = $a->aux_get('XL');  
					my $trim_right = $a->aux_get('XR');
					
					$trimmed = 1 if ((defined $trim_left) && ($p->qpos+1 <= $trim_left));
					$trimmed = 1 if ((defined $trim_right) && ($a->query->length-$p->qpos <= $trim_right));

					##also trimmed if up to next position and this is an indel.
					if ($indel != 0)
					{
						#remember that qpos is 0-indexed
						$trimmed = 1 if ((defined $trim_left) && ($p->qpos+1 <= $trim_left)); #so this is position <1
						$trimmed = 1 if ((defined $trim_right) && ($a->query->length-($p->qpos+1)+1 <= $trim_right)); #this is position+1
					}
					
					## These are the start and end coordinates of the aligned part of the read
					my ($q_start, $q_end) = Breseq::Shared::alignment_query_start_end($a, {no_reverse=>1});
					
					### Optionally, only count reads that completely match
					my $complete_match = 1;
					if ($settings->{require_complete_match})
					{
						$complete_match = ($q_start == 1) && ($q_end == $a->l_qseq);
						next if (!$complete_match);
					}
					### End complete match condition
					
					my $base = ($indel < $insert_count) ? '.' : substr($a->qseq,$p->qpos + $insert_count,1);
					##don't use bases without qualities!!
					next if ($base =~ /[nN]/);

					##### update coverage if this is not a deletion in read relative to reference
					### note that we count trimmed reads here, but not when looking for short indel mutations...	
					if ($redundancy == 1)
					{
						## this is only used when reporting coverage for within-read indels
						## NOT for calling deletions...
						
						$this_position_coverage->{unique}->{$strand}++;
						$pos_info->{$base}->{unique_cov}->{$strand}++;

						if ($indel > $insert_count)
						{
							$next_insert_count_exists = 1;
						}
					}
					else
					{
						$this_position_unique_only_coverage = 0;
						$this_position_coverage->{redundant}->{$strand} += 1/$redundancy;			
						$this_position_coverage->{raw_redundant}->{$strand}++;			
					}
	
					## EXPERIMENTAL -- moved above		
					##don't use information from trimmed reads!!
					next if ($trimmed);
					
					##don't use information from redundant reads!!
					next if ($redundancy > 1);				
					
					my $quality;
					if ($indel < $insert_count) 
					{
						## We would really like for this to not count any of the insertions unless
						## the read makes it out of the inserted region (this will be the case if it was aligned)
						##no information about gap if next position is end of alignment
						## THIS SHOULD NEVER BE TRUE, but it is sometimes...
						next ALIGNMENT if ($p->qpos+$indel+1 >= $q_end);
						#die if ($p->qpos+$indel+1 >= $q_end);

						my $average_quality = POSIX::floor(($a->qscore->[$p->qpos+$indel] + $a->qscore->[$p->qpos+$indel])/2);
						$quality = $average_quality;
					}
					else
					{
						$quality = $a->qscore->[$p->qpos+$insert_count];
					}

					## this is the coverage for SNP counts, tabulate AFTER skipping trimmed reads
					$pos_info->{$base}->{unique_trimmed_cov}->{$strand}++;
					
					##### this is for polymorphism prediction and making strings
					push @$pdata, { base => $base, quality => $quality, strand => $strand, fastq_file_index => $fastq_file_index };
								
					##### deal with base calls
					foreach my $hypothetical_base (@base_list)
					{				
						##sanity checks

						##is the error rate defined?
						my $base_key =  ($strand == +1) 
							? $hypothetical_base . $base
							: Breseq::Fastq::revcom($hypothetical_base) . Breseq::Fastq::revcom($base);
						
						if (!defined $log10_correct_rates->[$fastq_file_index]->{$quality}->{$base_key})
						{
							print "$fastq_file_index $quality $base_key\n";
							print Dumper($log10_correct_rates->[$fastq_file_index]);
							die;
						}

						##record evidence for and against this hypothetical base being the reference, given the observation
						$pr_base_hash->{$hypothetical_base} += $log10_correct_rates->[$fastq_file_index]->{$quality}->{$base_key};
						$pr_not_base_hash->{$hypothetical_base} += $log10_error_rates->[$fastq_file_index]->{$quality}->{$base_key};
					}
				} #end FOREACH READ
					
				## PER POSITION/INSERT COUNT
					
				#sum up coverage observations
				$this_position_coverage->{unique}->{total} = $this_position_coverage->{unique}->{-1} + $this_position_coverage->{unique}->{+1};
				$this_position_coverage->{redundant}->{total} = $this_position_coverage->{redundant}->{-1} + $this_position_coverage->{redundant}->{+1};
				$this_position_coverage->{raw_redundant}->{total} = $this_position_coverage->{raw_redundant}->{-1} + $this_position_coverage->{raw_redundant}->{+1};
				$this_position_coverage->{total} = $this_position_coverage->{unique}->{total} + $this_position_coverage->{redundant}->{total};
				
				$s->{coverage}->{unique_total}++ if ($this_position_unique_only_coverage && ($insert_count == 0));
					
				#we are trying to find the base with the most support;		
				#calculate ratios of each base to all other bases
				my $best_base;
				my $pr_call;
				foreach my $test_base (keys %$pr_base_hash)
				{			
					##obsolete code
					# 'P' is a special key for evidence against indels
					###next if ($test_base eq 'P');

					my $this_pr_call = $pr_base_hash->{$test_base} - $pr_not_base_hash->{$test_base};
					if ( (!defined $pr_call) || ($this_pr_call > $pr_call) )
					{
						$pr_call = $this_pr_call;
						$best_base = $test_base;
					}
				}
				
				#for the case where there are no counts
				my $e_value_call = 'NA';

				#otherwise correct to an e-value based on genome size
				#could correct only to number of unique positions???
				if (defined $pr_call)
				{
					$e_value_call = $pr_call - $log10_ref_length;
					$e_value_call = sprintf "%.1f", $e_value_call; #round immediately
				}
				
				##did we predict a base at this position?
				my $base_predicted = ($e_value_call ne 'NA' && ($e_value_call >= $settings->{mutation_log10_e_value_cutoff}));
				
				##print out SNP call information
				my $total_cov;
				$total_cov->{1} = 0;
				$total_cov->{-1} = 0;

				$line = "$pos\t$insert_count\t$ref_base\t$e_value_call";
				foreach my $base (@base_list)
				{			
					my $current_base_info = $pos_info->{$base};
					my $top_cov = $current_base_info->{unique_trimmed_cov}->{1};
					my $bot_cov = $current_base_info->{unique_trimmed_cov}->{-1};
					$total_cov->{1} += $top_cov;
					$total_cov->{-1} += $bot_cov;
					$line .= "\t$base\t" . "\t($bot_cov/$top_cov)";
				}
				print MUT "$line\n" if (defined $snps_all_tab_file_name);
			
				my ($base_string, $quality_string, $strand_string) = _pdata_to_strings(@$pdata);
			
				###
				## DELETION DELETION DELETION
				###

				#print to coverage file
				#update information on deletions
				if ($insert_count == 0)
				{
					if (!$settings->{no_deletion_prediction})
					{
						_check_deletion_completion($pos, $this_position_coverage, $e_value_call);
						_update_copy_number_variation($pos, $this_position_coverage, $ref_base); 							
					}
				}
				
				###
				## POLYMORPHISM POLYMORPHISM POLYMORPHISM
				###								
				my $polymorphism_predicted = 0;
				my $polymorphism;
				if ($settings->{polymorphism_prediction})
				{
					$polymorphism = _predict_polymorphism($pdata, $log10_correct_rates, $error_rates, $ref_base);
				
					if ($polymorphism)
					{
						$polymorphism->{log10_e_value} = ($polymorphism->{p_value} == 0) ? "999" : -log($total_ref_length * $polymorphism->{p_value})/log(10);						
						if ($polymorphism->{log10_e_value} >= 2)
						{
							$polymorphism_predicted = 1;
							#print Dumper($polymorphism);
						}
						$base_predicted = 1 if ($polymorphism_predicted);
					}					
				}				
				
				###
				## UNKNOWN UNKNOWN UNKNOWN
				###				
				if ($insert_count == 0)
				{
					_update_unknown_intervals($seq_id, $pos, $base_predicted, $this_position_unique_only_coverage);
				}
				
				## evaluate whether to call an actual mutation!				
				### skip if there is not enough evidence for a call or if it agrees with the reference
			#	next if (!$base_predicted);	
				next if (($e_value_call eq 'NA') || ($e_value_call < -$log10_ref_length));

				## mutation and polymorphism are exclusive predictions.
				my $mutation_predicted = (!$polymorphism_predicted) && ($best_base ne $ref_base);

				## bail if it's just the reference base and we aren't interested in polymorphisms...
				next INSERT_COUNT if (!$mutation_predicted && !$polymorphism_predicted);
				## bail if we are predicting polymorphisms, but there wasn't one

				my $mut;
				if ($mutation_predicted || $polymorphism_predicted)
				{
					#only once we are here can we be sure there is a $best_base!
					my $best_cov = $pos_info->{$best_base}->{unique_trimmed_cov};
					$best_cov->{1} = 0 if (!defined $best_cov->{1});
					$best_cov->{-1} = 0 if (!defined $best_cov->{-1});

					$mut->{type} = 'RA';
					$mut->{seq_id} = $seq_id;
					$mut->{position} = $pos;
					$mut->{insert_position} = $insert_count;
					$mut->{quality} = $e_value_call;		
					$mut->{tot_cov} = $total_cov->{-1} . "/" . $total_cov->{1};
					$mut->{new_cov} = $best_cov->{-1} . "/" . $best_cov->{1};
					
					$mut->{bases} = $base_string;
					$mut->{qualities} = $quality_string;
					$mut->{strands} = $strand_string;
				}
				if ($mutation_predicted)
				{
					$mut->{ref_base} = $ref_base;
					$mut->{new_base} = $best_base;		
					$mut->{frequency} = 1; ## this is not a polymorphism
					$mut->{marginal} = 1 if ($e_value_call < $settings->{mutation_log10_e_value_cutoff})
				}
				if ($polymorphism_predicted)
				{						
					
					$mut->{log10_e_value} = $polymorphism->{log10_e_value};
					$mut->{fisher_strand_p_value} = $polymorphism->{fisher_strand_p_value};
											
					if ($polymorphism->{first_base} eq $ref_base)
					{
						$mut->{frequency} = $polymorphism->{frequency};
						$mut->{ref_base} = $polymorphism->{first_base};
						$mut->{new_base} = $polymorphism->{second_base};
					}	
					elsif ($mut->{second_base} eq $ref_base)
					{
						$mut->{frequency} = 1-$polymorphism->{frequency};
						$mut->{ref_base} = $polymorphism->{second_base};
						$mut->{new_base} = $polymorphism->{first_base};
					}

					### NOTE: This neglects the case where neither the first nor second base is the reference base! Should almost never happen					
					# die if (($polymorphism->{first_base} ne $ref_base) && ($polymorphism->{second_base} ne $ref_base));
					else
					{
						$mut->{frequency} = $polymorphism->{frequency};
						$mut->{ref_base} = $polymorphism->{first_base};
						$mut->{new_base} = $polymorphism->{second_base};
						
						$mut->{error} = "polymorphic_without_reference_base";
					}
										
									
					$mut->{marginal} = 1 if ($mut->{log10_e_value} < $settings->{polymorphism_log10_e_value_cutoff});
					$mut->{marginal} = 1 if ($mut->{fisher_strand_p_value} < $settings->{polymorphism_fisher_strand_p_value_cutoff});
					$mut->{marginal} = 1 if ($mut->{frequency} < $settings->{polymorphism_frequency_cutoff});
				 	$mut->{marginal} = 1 if ($mut->{frequency} > 1-$settings->{polymorphism_frequency_cutoff});		
				}
				
				$gd->add($mut);
				
			} continue {
				$insert_count++;
			}#end INSERT COUNT	
		}; 

		$bam->pileup($seq_id, $pileup_function);

		###
		## Need to clean up deletions and unknowns that occur at the end
		###
		_check_deletion_completion($sequence_length+1); 		
		_update_unknown_intervals($seq_id, $sequence_length+1, 1); 		
		
		close COV;
		close MUT;	

# now this is handled after the evidence phase... 		
#		remove_deletions_overlapping_mutations(\@deletions, \@mutations);
	}
	
	##finally, write out the genome diff file
	my $predicted_mutation_genome_diff_file_name = $settings->file_name('predicted_mutation_genome_diff_file_name');	
	$gd->write($predicted_mutation_genome_diff_file_name);
}

sub remove_deletions_overlapping_mutations
{
	my ($deletions_ref, $mutations_ref) = @_;
	
	DEL: for (my $i=0; $i<scalar @$deletions_ref; $i++)
	{
		my $del = $deletions_ref->[$i];	 
		
		MUT: foreach my $mut (@$mutations_ref)
		{
			##remove the deletion if exactly the same as the SNP
			if ( ($del->{start} >= $mut->{start}) && ($del->{end} <= $mut->{end}) )
			{
				splice @$deletions_ref, $i, 1;
				$i--;
				last MUT;
			}
		}
	}
	
	##
	# Remove SNPs that overlap deletions
	##

	### This may now be unnecessary????
	MUT: for (my $i=0; $i<scalar @$mutations_ref; $i++)
	{
		my $mut = $mutations_ref->[$i];
		DEL: foreach my $del (@$deletions_ref)
		{
			##remove the SNP if it is located within the deletion
			if ( ($del->{start} <= $mut->{start}) && ($del->{end} >= $mut->{end}) )
			{
				splice @$mutations_ref, $i, 1;
				$i--;
				last DEL;
			}
		}
	}

}

sub _predict_polymorphism
{
	my $polymorphism;
	my ($info_list, $log10_correct_rates, $error_rates, $ref_base) = @_;
	# pdata is reference to list of observations
	# { base => $base, quality => $quality, strand => $strand, fastq_file_index => $fastq_file_index };							
		
	my $num_ref_base = 0;
	my $num_not_ref_base = 0;	
		
	#tabulate the raw counts of each base
	my $base_counts;
	foreach my $base (@base_list)
	{
		$base_counts->{$base} = 0;
	}	
	
	my $total = 0;
	my $total_minus = 0;
	my $total_plus = 0;
	my $mismatches = 0;
	foreach my $item (@$info_list)
	{
		if ($item->{base} ne $ref_base)
		{
			$mismatches++;
			$num_not_ref_base++;
		}
		else
		{
			$num_ref_base++;
		}	

		if ($item->{strand} == 1)
		{
			$total_plus++;
		}
		else
		{
			$total_minus++;
		}

		$base_counts->{$item->{base}}++;
		$total++;
	}
	
	##obviously no polymorphism! but do we believe consensus?
	return undef if ($mismatches == 0);
	##expand here to have partial polymorphism treatment replace the older mutation caller...
		
	#calculate the likelihood of observed bases given all positions were actually the reference base
	my $log10_likelihood_given_ref_base;
	foreach my $test_ref_base (@base_list)
	{
		$log10_likelihood_given_ref_base->{$test_ref_base} = 0;

		foreach my $item (@$info_list)
		{
			
			my $log10_pr_ref_base_given_obs = $log10_correct_rates->[$item->{fastq_file_index}]->{$item->{quality}}->{$test_ref_base . $item->{base}};
			die "ERROR: $item->{quality}, $test_ref_base$item->{base}" if (!defined $log10_pr_ref_base_given_obs);
			$log10_likelihood_given_ref_base->{$test_ref_base}  += $log10_pr_ref_base_given_obs;
		}
	}	

#	print Dumper($info_list);
#	print Dumper($log10_likelihood_given_ref_base);
	
	#we want the one with the most bases or the largest log likelihood (less negative is better)			
	my @bases_sorted_by_likelihood = sort { -($base_counts->{$a} <=> $base_counts->{$b}) || -($log10_likelihood_given_ref_base->{$a} <=> $log10_likelihood_given_ref_base->{$b}) } @base_list;
#	print Dumper(\@bases_sorted_by_likelihood);
#	print Dumper($base_counts);
	
	my $first_base = $bases_sorted_by_likelihood[0];
	my $second_base = $bases_sorted_by_likelihood[1];

	my @first_base_qualities;
	my @second_base_qualities;
	my $first_base_strand_hash = { '1' => 0, '-1' => 0};
	my $second_base_strand_hash = { '1' => 0, '-1' => 0};
	
	my $verbose = 0;
	
	print STDERR "1st = $first_base\n" if ($verbose);
	print STDERR "2nd = $second_base\n" if ($verbose);

	foreach my $item (@$info_list)
	{
		if ($item->{base} eq $first_base)
		{
			push @first_base_qualities, $item->{quality};
			$first_base_strand_hash->{$item->{strand}}++;
		}
		elsif ($item->{base} eq $second_base)
		{
			push @second_base_qualities, $item->{quality};
			$second_base_strand_hash->{$item->{strand}}++;
		}
	}

 	eval "use Statistics::FishersExactTest";
	my $fisher_strand_p_value = Statistics::FishersExactTest::fishers_exact($first_base_strand_hash->{1},$first_base_strand_hash->{-1},$second_base_strand_hash->{1},$second_base_strand_hash->{-1},1);
	$fisher_strand_p_value = sprintf "%.1e", $fisher_strand_p_value; #round immediately
	
#	print STDERR Dumper(\@first_base_qualities) if ($verbose);
#	print STDERR Dumper(\@second_base_qualities) if ($verbose);
	my $first_base_median_quality = median(\@first_base_qualities);
	my $second_base_median_quality = median(\@second_base_qualities);

	#use likelihood test to estimate probability of mixed model vs. one-base model
	my ($log10_likelihood_of_two_base_model, $max_likelihood_fr_first_base) = _find_best_log_likelihood($info_list, $error_rates, $first_base, $second_base);

	my $likelihood_ratio_test_value = -2*log(10)*($log10_likelihood_given_ref_base->{$first_base} - $log10_likelihood_of_two_base_model);
	
	eval "use Math::GSL::CDF";
	my $chi_squared_pr = Math::GSL::CDF::gsl_cdf_chisq_Q($likelihood_ratio_test_value, 1);

	print STDERR "$num_not_ref_base/$total $max_likelihood_fr_first_base // $log10_likelihood_given_ref_base->{$first_base} / $log10_likelihood_of_two_base_model = $likelihood_ratio_test_value $chi_squared_pr\n" if ($verbose);
	
	# my $print_p_cutoff = 1;
	# if ($chi_squared_pr <= $print_p_cutoff)
	# {	
	# 	print join("\t", 
	# 		'POLY', $ref_base, $first_base, $second_base, 
	# 		$base_counts->{'A'}, $base_counts->{'T'}, $base_counts->{'C'}, $base_counts->{'G'}, 
	# 		$mismatches, $total, $first_base_median_quality, $second_base_median_quality, 
	# 		$first_base_strand_hash->{-1}, $first_base_strand_hash->{1}, $second_base_strand_hash->{-1}, $second_base_strand_hash->{1},
	# 		$max_likelihood_fr_first_base, $chi_squared_pr,
	# 		join(",", @first_base_qualities), join(",", @second_base_qualities)
	# 	) . "\n";
	# }
	# if ($first_base ne $ref_base)
	# {
	# 	print join("\t", 
	# 		'MUT', $ref_base, $first_base, $second_base, 
	# 		$base_counts->{'A'}, $base_counts->{'T'}, $base_counts->{'C'}, $base_counts->{'G'}, 
	# 		$mismatches, $total, $first_base_median_quality, $second_base_median_quality, 
	# 		$first_base_strand_hash->{-1}, $first_base_strand_hash->{1}, $second_base_strand_hash->{-1}, $second_base_strand_hash->{1},
	# 		$max_likelihood_fr_first_base, $chi_squared_pr,
	# 		join(",", @first_base_qualities), join(",", @second_base_qualities)
	# 	) . "\n";		
	# }
	
	$polymorphism = {
		'frequency' => $max_likelihood_fr_first_base,
		'first_base' => $first_base,
		'second_base' => $second_base,
		'p_value' => $chi_squared_pr,
		'fisher_strand_p_value' => $fisher_strand_p_value,
	};
	
	return $polymorphism;
}

sub _find_best_log_likelihood
{
	my $verbose = 0;
	my ($info_list_ref, $error_rates, $first_base, $second_base) = @_;
	
	my $cur_pr_first_base = 1;
	my $cur_log_pr = _calculate_log10_likelihood($info_list_ref, $error_rates, $first_base, $second_base, $cur_pr_first_base);
	my $last_log_pr = -9999;
	my $last_pr_first_base = 1;

	while ($cur_log_pr >= $last_log_pr)
	{
		print "$cur_pr_first_base $cur_log_pr\n" if ($verbose);
		
		$last_log_pr = $cur_log_pr;
		$last_pr_first_base = $cur_pr_first_base;
		die if ($cur_pr_first_base < 0);

		$cur_pr_first_base -= 0.001;
		$cur_log_pr = _calculate_log10_likelihood($info_list_ref, $error_rates, $first_base, $second_base, $cur_pr_first_base);
	}
	
	return ($last_log_pr, $last_pr_first_base);
}

sub _calculate_log10_likelihood
{
	my ($info_list_ref, $error_rates, $first_base, $second_base, $pr_first_base) = @_;
	
	my $log10_likelihood = 0;	
	foreach my $item (@$info_list_ref)
	{
		my $pr_ref_base_given_obs 
			= $pr_first_base * $error_rates->[$item->{fastq_file_index}]->{$item->{quality}}->{$first_base . $item->{base}}
			+ (1-$pr_first_base) * $error_rates->[$item->{fastq_file_index}]->{$item->{quality}}->{$second_base . $item->{base}};
		
		print "Probability of error is zero.\n  Fastq File = $item->{fastq_file_index}\n"
		  .  "  Quality = $item->{quality}\n  Error = $first_base$item->{base} $second_base$item->{base}\n"
		  if ($pr_ref_base_given_obs == 0);
			
		$log10_likelihood += log($pr_ref_base_given_obs);		
	}	
	$log10_likelihood /= log(10);
	return $log10_likelihood;
}

sub _pdata_to_strings
{
	my @pdata = @_;
	
	@pdata = sort { ($a->{base} cmp $b->{base}) || -($a->{quality} <=> $b->{quality})  } @pdata;
	
	my $base_string = join '', map ($_->{base}, @pdata);
	my $quality_string = join '', map ( (chr($_->{quality}+33)), @pdata);
	my $strand_string = join '', map ( (($_->{strand}) == -1) ? 1 : 0, @pdata);

	return ($base_string, $quality_string, $strand_string);
}

sub median {
    my $rpole = shift;
    my @pole = @$rpole;

    my $ret;

	return "NA" if (@pole == 0);

    @pole = sort {$a <=> $b} @pole;

	if (@pole == 1) {
		$ret = $pole[0];
	}
    elsif( (@pole % 2) == 1 ) {
        $ret = $pole[((@pole+1) / 2)-1];
    } else {
        $ret = ($pole[(@pole / 2)-1] + $pole[@pole / 2]) / 2;
    }

    return $ret;
}

sub calculate_confidence_interval
{
	my $trials = 1000;
	my $precision = 1000;
	my $verbose = 0;
	
	my ($error_rates, $ref_base, $new_base, $ref_base_quality_list, $new_base_quality_list) = @_;
	my $ml = 0;
	my $ml_prob = 0;
	my $lower = 0;
	my $upper = 1;
	
	my $num_ref_base = scalar @$ref_base_quality_list;
	my $num_new_base = scalar @$new_base_quality_list;	
	
	#fill in our error chart as ref/new given (ref or new)
	#0 = ref
	#1 = new
	
	my $chance_of_new_hash;
	foreach (my $fastq_file_index=0; $fastq_file_index < scalar @$error_rates; $fastq_file_index++)
	{
		
		foreach my $qual (keys %{$error_rates->[$fastq_file_index]})
		{
			my $key = "$fastq_file_index+$qual";
		
			my @item = (
				$error_rates->[$fastq_file_index]->{$qual}->{$new_base . $new_base} / ($error_rates->[$fastq_file_index]->{$qual}->{$new_base . $ref_base} + $error_rates->[$fastq_file_index]->{$qual}->{$new_base . $new_base}),
				$error_rates->[$fastq_file_index]->{$qual}->{$ref_base . $new_base} / ($error_rates->[$fastq_file_index]->{$qual}->{$ref_base . $ref_base} + $error_rates->[$fastq_file_index]->{$qual}->{$ref_base . $new_base}), 
			);
		
			$chance_of_new_hash->{$key} = \@item; 
		}
	}
				
	#val fraction is for new base
	my $val = 0;
	my $since_last_upper = 0;

	while ($val <= 1)
	{
		my $count_lower = 0;
		my $count_upper = 0;
		my $count_actual = 0;
	
		for my $t (1..$trials)
		{
			my $new_count = 0;
			foreach my $qual (@$ref_base_quality_list, @$new_base_quality_list)
			{
				my $base = (rand 1.0 <= $val) ? 0 : 1;
				$new_count++ if (rand 1.0 < $chance_of_new_hash->{$qual}->[$base]);
			}
		
			$count_lower++ if ($new_count <= $num_new_base);
			$count_upper++ if ($new_count >= $num_new_base);
			$count_actual++ if ($new_count == $num_new_base);
		}
		$val+=$precision;
		
		#upper bound
		if (($count_lower/$trials >= 0.025))
		{
			$upper = $val;
			$since_last_upper = 0;
		}
		else
		{
			$since_last_upper += $precision;
		}
		
		#lower bound
		my $lower_frac = $count_upper/$trials;
		$lower = $val if (($count_upper/$trials <= 0.025) && ($val > $lower));


		if ($count_actual > $ml_prob)
		{
			$ml = $val;
			$ml_prob = $count_actual;
		}
		last if ($since_last_upper >= 0.05);
	
		print "$val, $count_lower, $count_upper, $count_actual, $since_last_upper, $lower_frac\n" if ($verbose);
	}
		
	return ($ml, $lower, $upper);
}

return 1;