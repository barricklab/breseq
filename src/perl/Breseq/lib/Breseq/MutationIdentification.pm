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

use Bio::Root::Root;
use Bio::DB::Sam;

use Breseq::GenomeDiff;
use Breseq::Shared;
use Breseq::ErrorCalibration;


package Breseq::MutationIdentification;
use vars qw(@ISA);
@ISA = qw( Bio::Root::Root );

use Data::Dumper;

use Breseq::GenomeDiff;
use Statistics::Distributions;


our @base_list = ('A', 'T', 'C', 'G', '.');

#does both within-read SNPs+indels and missing coverage deletions
sub identify_mutations
{
	our ($settings, $summary, $ref_seq_info, $error_rates) = @_;
		
	my $reference_fasta_file_name = $settings->file_name('reference_fasta_file_name');
	my $reference_bam_file_name = $settings->file_name('reference_bam_file_name');
		
	## check to see if identify_mutations++ is available:
	my $cident_mut = $settings->ctool('identify_mutations', 1);
	
	## fall back to perl is requested or if predicting polymorphisms
	if ( (!$settings->{perl_identify_mutations}) && (!$settings->{polymorphism_prediction}) && (defined $cident_mut) )
	{
		my $coverage_fn = $settings->file_name('unique_only_coverage_distribution_file_name', {'@'=>""});
		my $error_dir = `dirname $coverage_fn`;
		chomp $error_dir; $error_dir .= '/';
		my $this_predicted_mutation_file_name = $settings->file_name('predicted_mutation_file_name', {'@'=>""});
		my $output_dir = `dirname $this_predicted_mutation_file_name`;
		chomp $output_dir; $output_dir .= '/';
		my $readfiles = join(" --readfile ", $settings->read_files);
		my $cmdline = "$cident_mut --bam $reference_bam_file_name --fasta $reference_fasta_file_name --readfile $readfiles";
		$cmdline .= " --error_dir $error_dir";
		my $ra_mc_genome_diff_file_name = $settings->file_name('ra_mc_genome_diff_file_name');	
		$cmdline .= " --genome_diff $ra_mc_genome_diff_file_name";
		$cmdline .= " --output $output_dir";
		if(defined $settings->{mutation_log10_e_value_cutoff}) {
			$cmdline .= " --mutation_cutoff $settings->{mutation_log10_e_value_cutoff}"; # defaults to 2.0.
		}		
		my $coverage_tab_file_name = $settings->file_name('complete_coverage_text_file_name', {'@'=>""});
		my $coverage_dir = `dirname $coverage_tab_file_name`;
		chomp $coverage_dir; $coverage_dir .= '/';
		$cmdline .= " --coverage_dir $coverage_dir";
		if(defined $settings->{deletion_propagation_cutoff}) {
			$cmdline .= " --deletion_propagation_cutoff $settings->{deletion_propagation_cutoff}"; # defaults to 28.0.
		}
		if((defined $settings->{no_deletion_prediction}) && ($settings->{no_deletion_prediction})) {
			$cmdline .= " --predict_deletions 0";
		} else {
			$cmdline .= " --predict_deletions 1"; # defaults TO predicting deletions.
		}
		if((defined $settings->{polymorphism_prediction}) && ($settings->{polymorphism_prediction})) {
			$cmdline .= " --predict_polymorphisms 1";
		} else {
			$cmdline .= " --predict_polymorphisms 0"; # defaults to NOT predicting polymorphisms.
		}

		Breseq::Shared::system($cmdline);
		return; # identify_mutations++ worked, so we're all done here.
	}
	
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
 
	my $polymorphism_statistics_input_fh;
	if ($settings->{polymorphism_prediction})
	{
		my $polymorphism_statistics_input_file_name = $settings->file_name('polymorphism_statistics_input_file_name');
		open $polymorphism_statistics_input_fh, ">$polymorphism_statistics_input_file_name" or die "Could not open file: $polymorphism_statistics_input_file_name";
		print $polymorphism_statistics_input_fh +join("\t",
			'position', 'insert_position', 'frequency', 'log10_base_likelihood', 'new_top_strand', 'new_bot_strand', 'ref_top_strand', 'ref_bot_strand', 'new_quals', 'ref_quals'
		) . "\n";
	}
	
	REFERENCE: foreach our $seq_id (@seq_ids)
	{				
		my $this_predicted_mutation_file_name = $settings->file_name('predicted_mutation_file_name', {'@'=>$seq_id});
		my $this_tiled_coverage_tab_file_name = $settings->file_name('tiled_coverage_text_file_name', {'@'=>$seq_id});	
		
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
		our $left_outside_coverage_item = undef
		our $left_inside_coverage_item = undef;
		our $last_position_coverage = undef;
		# additions to check whether we are within initial or end redundant positions
		our $last_deletion_redundant_start_position = undef;
		our $last_deletion_redundant_end_position = undef;	
		our $redundant_reached_zero = 0;

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
			
			## when called at the end of a fragment, the position is fragment length +1
			## and $this_position_coverage is undefined
				
			# we need to fill in reference positions with NO reads aligned to them
			# pileup won't be called at these positions
			foreach (my $i = $last_position_coverage_printed + 1; $i < $pos; $i++)
			{
				if (!defined $last_deletion_start_position)
				{
					## special treatment for the beginning of a genome fragment
					if ($last_position_coverage_printed == 0)
					{
						$left_outside_coverage_item =  {
							unique => {'1'=>'NA', '-1'=>'NA', 'total' => 'NA' },
							redundant => {'1'=>'NA', '-1'=>'NA', 'total' => 'NA' }
						};
						$left_inside_coverage_item = { 
							unique => {'1'=>0, '-1'=>0, 'total' => 0 },
							redundant => {'1'=>0, '-1'=>0, 'total' => 0 }
						};
					}
					## normal treatment is that coverage went to zero
					else
					{
						$left_inside_coverage_item = { 
							unique => {'1'=>0, '-1'=>0, 'total' => 0 },
							redundant => {'1'=>0, '-1'=>0, 'total' => 0 }
						};
						$left_outside_coverage_item = $last_position_coverage;	
					}
					$last_deletion_start_position = $last_position_coverage_printed+1;
					$last_deletion_redundant_start_position = $last_position_coverage_printed+1 if (!defined $last_deletion_redundant_start_position);
				}
				
				$this_deletion_reaches_seed_value = 1;
				$redundant_reached_zero = 1;
				$last_position_coverage = { 
					unique => {'1'=>0, '-1'=>0, 'total' => 0 },
					redundant => {'1'=>0, '-1'=>0, 'total' => 0 }
				};

				print COV join("\t", 0, 0, 0, 0, 0, 0, 'NA', $i) . "\n";
			}
			$last_position_coverage_printed = $pos;
			
			## called with an undef $this_position_coverage at the end of the genome
			if (defined $this_position_coverage)
			{
				# Print this information to the coverage file (used as input by R plotting script)
				my $tu = $this_position_coverage->{unique};
				my $tr = $this_position_coverage->{redundant};
				my $trr = $this_position_coverage->{raw_redundant};
				print COV join("\t", $tu->{-1}, $tu->{1}, $tr->{-1}, $tr->{1}, $trr->{-1}, $trr->{1}, $e_value_call, $pos) . "\n";
		
				## UNIQUE COVERAGE
				#start a new possible deletion if we fall below the propagation cutoff
				if ($this_position_coverage->{unique}->{total} <= $deletion_propagation_cutoff)
				{	
					if (!defined $last_deletion_start_position)
					{
						$last_deletion_start_position = $pos;
						$left_outside_coverage_item = $last_position_coverage;
						$left_inside_coverage_item = $this_position_coverage;
					}
				}
				
				##keep track of whether we've encountered the seed value for the current deletion
				if ($this_position_coverage->{total} <= $deletion_seed_cutoff)
				{
					$this_deletion_reaches_seed_value = 1;
				}
				
				## REDUNDANT COVERAGE
				## updated only if we are currently within a deletion
				if (defined $last_deletion_start_position)
				{					
					if ($this_position_coverage->{redundant}->{total} == 0)
					{
						$redundant_reached_zero = 1; #switch from adjusting start to end
						undef $last_deletion_redundant_end_position;
					}
					elsif ($this_position_coverage->{redundant}->{total} > 0)
					{
						## if there is any redundant coverage remember the start (until we find zero redundant coverage)
						if (!$redundant_reached_zero)
						{
							$last_deletion_redundant_start_position = $pos;
						}
						## if we are working on the right side update the end position if it is not already defined.
						else
						{
							$last_deletion_redundant_end_position = $pos if (!defined $last_deletion_redundant_end_position);
						}
					}
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
							redundant => {'1'=>'NA', '-1'=>'NA', 'total' => 'NA' }
						};
					}
					
					my $last_deletion_end_position = $pos-1;
					$last_deletion_redundant_end_position = $last_deletion_end_position	if (!defined $last_deletion_redundant_end_position);
					$last_deletion_redundant_start_position = $last_deletion_start_position	if (!defined $last_deletion_redundant_start_position);

					my $del = {
						type => 'MC',
						seq_id => $seq_id,
						start => $last_deletion_start_position,
						end => $last_deletion_end_position,
						start_range => $last_deletion_redundant_start_position - $last_deletion_start_position,
						end_range => $last_deletion_end_position - $last_deletion_redundant_end_position,
						left_outside_cov => $left_outside_coverage_item->{unique}->{total},
						left_inside_cov => $left_inside_coverage_item->{unique}->{total},
						right_inside_cov => $last_position_coverage->{unique}->{total},
						right_outside_cov => $this_position_coverage->{unique}->{total},
					};
					$del->{left_inside_cov} = 'NA' if (!defined $del->{left_inside_cov});
					$del->{right_inside_cov} = 'NA' if (!defined $del->{right_inside_cov});
					
					$del->{left_outside_cov} = 'NA' if (!defined $del->{left_outside_cov});
					$del->{right_outside_cov} = 'NA' if (!defined $del->{right_outside_cov});
					
					$gd->add($del);					
				}

				#reset the search
				$this_deletion_reaches_seed_value = 0;
				$redundant_reached_zero = 0;
				undef $last_deletion_start_position;
				undef $last_deletion_redundant_start_position;
				undef $last_deletion_redundant_end_position;
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
				## This works even with nonstandard reference bases.

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

					my $qpos = $p->qpos;	#0-indexed
					my $redundancy = $a->aux_get('X1');
					my $fastq_file_index = $a->aux_get('X2');
					my $reversed = $a->reversed;
					my $strand = $reversed ? -1 : +1;

					## Handle trimming
					## Note that trimming INCLUDES the unaligned (padded) bases on each end
					my $trimmed = 0;
					my $trim_left = $a->aux_get('XL');  #1-indexed
					my $trim_right = $a->aux_get('XR');  #1-indexed
					
					$trimmed = 1 if ((defined $trim_left) && ($p->qpos+1 <= $trim_left));
					$trimmed = 1 if ((defined $trim_right) && ($a->query->length-($p->qpos+1) <= $trim_right));
					
					## These are the start and end coordinates of the aligned part of the read
					my ($q_start, $q_end) = ($a->query->start-1, $a->query->end-1); #0-indexed
					
					### This is now handled during AlignmentCorrection so error model is correct.
					### Optionally, only count reads that completely match --->
					#my $complete_match = 1;
					#if ($settings->{require_complete_match})
					#{
					#	$complete_match = ($q_start+1 == 1) && ($q_end+1 == $a->l_qseq);
					#	next if (!$complete_match);
					#}
					###<--- End complete match condition
					
					
					## What is the read base?
					my $base = ($indel < $insert_count) ? '.' : substr($a->qseq,$qpos + $insert_count,1);
					next ALIGNMENT if ($base eq 'N'); ## We don't want to use read N bases for anything 

					##### Update coverage if this is not a deletion in read relative to reference
					### Count trimmed reads here, but not when looking for short indel mutations...	
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
	
					##don't use information from trimmed reads!!
					next ALIGNMENT if ($trimmed);
					
					##don't use information from redundant reads!!
					next ALIGNMENT if ($redundancy > 1);				
										
					my $quality;
					
					## Deletion in read relative to reference...
					## Quality is of the NEXT base in the read, and check that it is not an N
					## Note: This is for a deletion when $insert_count == 0 and against an insertion when $insert_count > 0
					
					if ($indel == -1)  
					{					
						my $mqpos = $qpos + 1 - $reversed;
						my $check_base = substr($a->qseq,$mqpos,1);
						next ALIGNMENT if ($check_base eq 'N');
						$quality = $a->qscore->[$mqpos];
					}

					## Substitution in read relative to reference...
					## Quality is of the current base in the read, we have ALREADY checked that it is not an N					
					elsif ($insert_count == 0)
					{
						$quality = $a->qscore->[$qpos];
					}
					
					## Insertion in read relative to reference...
					## Quality is of the NEXT base in the read, and check that it is not an N
					## Note that it is possible this read base may be a '.' (supporting the non-insert call)
					else ## if ($insert_count > 0) 
					{		
						my $max_offset = $insert_count;
						$max_offset = $indel if ($indel < $max_offset);
						my $mqpos = $qpos + $max_offset + 1 - $reversed;
						
						## Check bounds: it's possible to go past the end of the read because
						## this is the last base of this read, but other reads have inserted bases
						next ALIGNMENT if ($mqpos > $q_end);
												
						my $check_base = substr($a->qseq,$mqpos,1);
						next ALIGNMENT if ($check_base eq 'N');
						
						$quality = $a->qscore->[$mqpos];
					}
					
					## We may want to ignore all bases below a certain quality when calling mutations and polymorphisms
					## This is the check for whether the base fails; it should be after coverage counting
					next ALIGNMENT if ( $settings->{base_quality_cutoff} && ($quality < $settings->{base_quality_cutoff}) );

					## this is the coverage for SNP counts, tabulate AFTER skipping trimmed reads
					$pos_info->{$base}->{unique_trimmed_cov}->{$strand}++;
					
					##### this is for polymorphism prediction and making strings
					push @$pdata, { base => $base, quality => $quality, strand => $strand, fastq_file_index => $fastq_file_index };
								
					#print STDERR "========\n$pos\n";
					
					##### deal with base calls
					foreach my $hypothetical_base (@base_list)
					{				
						##sanity checks

						##is the error rate defined?
						my $base_key =  ($strand == +1) 
							? $hypothetical_base . $base
							: Breseq::Fastq::complement($hypothetical_base . $base);
						
						if (!defined $log10_correct_rates->[$fastq_file_index]->{$quality}->{$base_key})
						{
							print "=== Error Rate not defined === \n";
							print "  Fastq File Index: $fastq_file_index\n";
							print "  Base Quality: $quality\n";
							print "  Base Key: $base_key\n";
							print "  Encountered in read: " . $a->qname . "\n";
							die;
						}

						##record evidence for and against this hypothetical base being the reference, given the observation
						my $pr = $log10_correct_rates->[$fastq_file_index]->{$quality}->{$base_key};
						my $notpr = $log10_error_rates->[$fastq_file_index]->{$quality}->{$base_key};
						#print STDERR "pr: $pr not: $notpr\n";
						$pr_base_hash->{$hypothetical_base} += $pr;#$log10_correct_rates->[$fastq_file_index]->{$quality}->{$base_key};
						$pr_not_base_hash->{$hypothetical_base} += $notpr;#$log10_error_rates->[$fastq_file_index]->{$quality}->{$base_key};
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
					$polymorphism = _predict_polymorphism($settings, $pdata, $log10_correct_rates, $error_rates, $ref_base);
				
					if ($polymorphism)
					{
						#print "Position: $pos\n";
						
						$polymorphism->{log10_e_value} = 'ND'; 
						if ($polymorphism->{p_value} ne 'ND')
						{
							$polymorphism->{log10_e_value} = ($polymorphism->{p_value} == 0) ? "999" : -log($total_ref_length * $polymorphism->{p_value})/log(10);						
						}
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

				### 
				## Code from here to end of loop is all about adding to the genome diff.
				###
				## Fields common to consensus mutations and polymorphisms
				my $mut;				
				$mut->{type} = 'RA';
				$mut->{seq_id} = $seq_id;
				$mut->{position} = $pos;
				$mut->{insert_position} = $insert_count;
				$mut->{quality} = $e_value_call;		
				$mut->{snp_quality} = $e_value_call;		

				## code that prints out even more information
				## slow because it sorts things, and not necessary
				if (0)
				{
					my ($base_string, $quality_string, $strand_string) = _pdata_to_strings(@$pdata);
					$mut->{bases} = $base_string;
					$mut->{qualities} = $quality_string;
					$mut->{strands} = $strand_string;
				}
				
				if ($mutation_predicted)
				{
					$mut->{ref_base} = $ref_base;
					$mut->{new_base} = $best_base;		
					$mut->{frequency} = 1; ## this is not a polymorphism
					
					Breseq::GenomeDiff::add_reject_reason($mut, "EVALUE") if ($e_value_call < $settings->{mutation_log10_e_value_cutoff});					
				}
				if ($polymorphism_predicted)
				{	
					$mut->{quality} = $polymorphism->{log10_e_value};		
					#$mut->{fisher_strand_p_value} = $polymorphism->{fisher_strand_p_value};
					
					# the frequency returned is the probability of the FIRST base
					# we want to quote the probability of the second base (the change from the reference).
					if ($polymorphism->{first_base} eq $ref_base)
					{
						$mut->{frequency} = 1-$polymorphism->{frequency};
						$mut->{ref_base} = $polymorphism->{first_base};
						$mut->{new_base} = $polymorphism->{second_base};
					}	
					elsif ($polymorphism->{second_base} eq $ref_base)
					{
						$mut->{frequency} = $polymorphism->{frequency};
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
					
					#round
					$mut->{frequency} = sprintf "%.4f", $mut->{frequency};
					
					Breseq::GenomeDiff::add_reject_reason($mut, "EVALUE") if ($mut->{quality} < $settings->{polymorphism_log10_e_value_cutoff});
					
					###
					## Print input file for R
					###
					my $ref_cov = $pos_info->{$mut->{ref_base}}->{unique_trimmed_cov};
					my $new_cov = $pos_info->{$mut->{new_base}}->{unique_trimmed_cov};
					
					my @ref_base_qualities;
					my @new_base_qualities;
					foreach my $item (@$pdata)
					{
						if ($item->{base} eq $mut->{ref_base})
						{
							push @ref_base_qualities, $item->{quality};
						}
						elsif ($item->{base} eq $mut->{new_base})
						{
							push @new_base_qualities, $item->{quality};
						}
					}
					
					my $ref_quality_string = join ',', @ref_base_qualities;
					my $new_quality_string = join ',', @new_base_qualities;
					
					print $polymorphism_statistics_input_fh +join( "\t",
						$mut->{position}, $mut->{insert_position}, $mut->{frequency}, $polymorphism->{log10_base_likelihood}, $new_cov->{1}, $new_cov->{-1}, $ref_cov->{1}, $ref_cov->{-1}, $new_quality_string, $ref_quality_string
					) . "\n";
					###
					## End printing input file for R
					###
					
					Breseq::GenomeDiff::add_reject_reason($mut, "FREQ") if ($mut->{frequency} < $settings->{polymorphism_frequency_cutoff});
					Breseq::GenomeDiff::add_reject_reason($mut, "FREQ") if ($mut->{frequency} > 1-$settings->{polymorphism_frequency_cutoff});		
				}
				
				
				## More fields common to consensus mutations and polymorphisms
				## ...now that ref_base and new_base are defined
				my $ref_cov = $pos_info->{$mut->{ref_base}}->{unique_trimmed_cov};
				$mut->{ref_cov} = $ref_cov->{1} . "/" . $ref_cov->{-1};
				
				my $new_cov = $pos_info->{$mut->{new_base}}->{unique_trimmed_cov};
				$mut->{new_cov} = $new_cov->{1} . "/" . $new_cov->{-1};
				
				$mut->{tot_cov} = $total_cov->{1} . "/" . $total_cov->{-1};
				
				$gd->add($mut);
				
			} continue {
				$insert_count++;
			}#end INSERT COUNT	
		}; 

		#for testing
		#$bam->pileup("$seq_id:1-200", $pileup_function);
		$bam->pileup("$seq_id", $pileup_function);

		###
		## Need to clean up deletions and unknowns that occur at the end
		###
		_check_deletion_completion($sequence_length+1); 		
		_update_unknown_intervals($seq_id, $sequence_length+1, 1); 		
		
		close COV;
		close MUT;	
	}
	
	if ($settings->{polymorphism_prediction})
	{
		close $polymorphism_statistics_input_fh;
	}
	
	##finally, write out the genome diff file
	my $ra_mc_genome_diff_file_name = $settings->file_name('ra_mc_genome_diff_file_name');	
	$gd->write($ra_mc_genome_diff_file_name);
}

sub _predict_polymorphism
{
	my $verbose = 0;
	
	my $polymorphism;
	my ($settings, $info_list, $log10_correct_rates, $error_rates, $ref_base) = @_;
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
		($item->{base} ne $ref_base) ? $mismatches++ && $num_not_ref_base++ : $num_ref_base++;
		($item->{strand} == 1) ? $total_plus++ : $total_minus++;

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
			my $base_change_key = $test_ref_base . $item->{base};
		 	$base_change_key = Breseq::Fastq::complement($base_change_key) if ($item->{strand} == -1); 
			my $log10_pr_ref_base_given_obs = $log10_correct_rates->[$item->{fastq_file_index}]->{$item->{quality}}->{$base_change_key};
			die "ERROR: $item->{quality}, $test_ref_base$item->{base}" if (!defined $log10_pr_ref_base_given_obs);
			$log10_likelihood_given_ref_base->{$test_ref_base}  += $log10_pr_ref_base_given_obs;
		}
	}	

#	print Dumper($info_list) if ($verbose);
#	print Dumper($log10_likelihood_given_ref_base) if ($verbose);
	
	#we want the one with the most bases or the largest log likelihood (less negative is better)			
	my @bases_sorted_by_likelihood = sort { -($base_counts->{$a} <=> $base_counts->{$b}) || -($log10_likelihood_given_ref_base->{$a} <=> $log10_likelihood_given_ref_base->{$b}) } @base_list;
#	print Dumper(\@bases_sorted_by_likelihood);
#	print Dumper($base_counts);
	
	my $first_base = $bases_sorted_by_likelihood[0];
	my $second_base = $bases_sorted_by_likelihood[1];

	print "1st = $first_base\n" if ($verbose);
	print "2nd = $second_base\n" if ($verbose);

	## we may not want to deal with indels
	return undef if ($settings->{no_indel_polymorphisms} && (($first_base eq '.') || ($second_base eq '.')) );

	##Don't need to do this...
	## from now on we ignore any observation that was not the first or second base!
	##my $total_observation_num = scalar @$info_list;
	##@$info_list = grep { ($_->{base} eq $first_base) || ($_->{base} eq $second_base) } @$info_list;
	##my $top_two_observation_num = scalar @$info_list;
	##print "Observations of the top two bases $top_two_observation_num / $total_observation_num total\n" if ($verbose);

	my @first_base_qualities;
	my @second_base_qualities;
	my $first_base_strand_hash = { '1' => 0, '-1' => 0};
	my $second_base_strand_hash = { '1' => 0, '-1' => 0};
	my $total_strand_hash = { '1' => 0, '-1' => 0};
	
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
		$total_strand_hash->{$item->{strand}}++;
	}

	##Likelihood of observing alignment if all sequenced bases were actually the first base
	my $log10_likelihood_of_one_base_model = $log10_likelihood_given_ref_base->{$first_base};

	## Maximum likelihood of observing alignment if sequenced bases were a mixture of the top two bases
	my ($log10_likelihood_of_two_base_model, $max_likelihood_fr_first_base) = _best_two_base_model_log10_likelihood($info_list, $error_rates, $first_base, $second_base);
#	my $pr_first_base = ($first_base_strand_hash->{'-1'} + $first_base_strand_hash->{'1'}) / ($first_base_strand_hash->{'-1'} + $first_base_strand_hash->{'1'} + $second_base_strand_hash->{'-1'} + $second_base_strand_hash->{'1'});
#	my ($log10_likelihood_of_two_base_model) = _calculate_two_base_model_log10_likelihood($info_list, $error_rates, $first_base, $second_base);

## this requires evidence to reject two base model that is lacking when there are too few obs.
#	my $strand_1_pr_first_base = ($first_base_strand_hash->{'-1'}) / ($first_base_strand_hash->{'-1'} + $second_base_strand_hash->{'-1'});
#	my $strand_2_pr_first_base = ($first_base_strand_hash->{'1'}) / ($first_base_strand_hash->{'1'} + $second_base_strand_hash->{'1'});
#	my ($log10_likelihood_of_two_base_strand_model) = _calculate_two_base_model_strand_log10_likelihood($info_list, $error_rates, $first_base, $second_base,$strand_1_pr_first_base, $strand_2_pr_first_base);

	## Maximum likelihood of observing alignment if sequenced bases were actually the first base, except
	## that on one strand there was a much higher chance of observing the second base
	my $log10_likelihood_of_strand_bias_model;	

	##test models on both strands
	my $strand_biased_error_rate;
	foreach my $s ('-1', '1')
	{
		my $strand_total = $total_strand_hash->{$s};
		$strand_biased_error_rate->{$s} = 0;
		next if ($strand_total == 0);
		$strand_biased_error_rate->{$s} = $second_base_strand_hash->{$s} / $strand_total;
		
		my $this_log10_likelihood_of_strand_bias_model = _calculate_strand_bias_model_log10_likelihood($info_list, $error_rates, $first_base, $second_base, $s, $strand_biased_error_rate->{$s});
		$log10_likelihood_of_strand_bias_model = $this_log10_likelihood_of_strand_bias_model if (!defined $log10_likelihood_of_strand_bias_model || ($this_log10_likelihood_of_strand_bias_model > $log10_likelihood_of_strand_bias_model));
	}

	print "Num First Base (+/-): $first_base_strand_hash->{1}/$first_base_strand_hash->{-1}\n" if ($verbose);
	print "Num Second Base (+/-): $second_base_strand_hash->{1}/$second_base_strand_hash->{-1}\n" if ($verbose);

	print "== Best Base Model ==\n" if ($verbose);
	print "  Log10 Likelihood Best Base Only = $log10_likelihood_given_ref_base->{$first_base}\n" if ($verbose);

	$Statistics::Distributions::SIGNIFICANT = 50; ## we need many more significant digits than the default 5

	## this is the natural logarithm of the ratio of the probabilities
	my $likelihood_ratio_test_value = -2*log(10)*($log10_likelihood_of_one_base_model - $log10_likelihood_of_two_base_model);
	my $chi_squared_pr = 'ND';
	$chi_squared_pr = Statistics::Distributions::chisqrprob(1, $likelihood_ratio_test_value);
	
	print "== 2 Base Model ==\n" if ($verbose);
	print "  Log10 Likelihood 2 Base Model = $log10_likelihood_of_two_base_model\n" if ($verbose);
	print "  LR test value = $likelihood_ratio_test_value\n" if ($verbose);
	print "  Chi-squared Pr = $chi_squared_pr\n" if ($verbose);

	my $likelihood_ratio_test_value_strand_bias = -2*log(10)*($log10_likelihood_of_one_base_model - $log10_likelihood_of_strand_bias_model);
	my $chi_squared_pr_strand_bias = 'ND';
	$chi_squared_pr_strand_bias = Statistics::Distributions::chisqrprob(1, $likelihood_ratio_test_value);

	print "== Strand Bias Model ==\n" if ($verbose);
	print "  Log10 Likelihood Strand Bias Model= $log10_likelihood_of_strand_bias_model\n" if ($verbose);
	print "  LR test value = $likelihood_ratio_test_value_strand_bias\n" if ($verbose);
	print "  Chi-squared Pr = $chi_squared_pr_strand_bias\n" if ($verbose);

	## we need the two base model to beat the strand bias model by a fair bit
	my $p_value_two_base_vs_strand_bias =  10**($log10_likelihood_of_strand_bias_model - $log10_likelihood_of_two_base_model);
	print "==> P-value two-base model vs. strand bias model = $p_value_two_base_vs_strand_bias\n" if ($verbose);

	$polymorphism = {
		'frequency' => $max_likelihood_fr_first_base,
		'first_base' => $first_base,
		'second_base' => $second_base,
		'first_base_strand_coverage' => $first_base_strand_hash, 
		'second_base_strand_coverage' => $second_base_strand_hash, 
		'p_value' => $chi_squared_pr,
		'log10_base_likelihood' => $log10_likelihood_given_ref_base->{$first_base} - $log10_likelihood_of_two_base_model,
		'p_value_strand_bias' => $chi_squared_pr_strand_bias,
		'log10_base_likelihood_strand_bias' => $log10_likelihood_of_one_base_model - $log10_likelihood_of_strand_bias_model,
		'p_value_two_base_vs_strand_bias' => $p_value_two_base_vs_strand_bias,
	};
	
	return $polymorphism;
}

sub _best_two_base_model_log10_likelihood
{
	my $verbose = 0;
	my ($info_list, $error_rates, $first_base, $second_base) = @_;
	
	my $cur_pr_first_base = 1;
	my $cur_log_pr = _calculate_two_base_model_log10_likelihood($info_list, $error_rates, $first_base, $second_base, $cur_pr_first_base);
	my $last_log_pr = -9999;
	my $last_pr_first_base = 1;

	while ($cur_log_pr >= $last_log_pr)
	{
		print "$cur_pr_first_base $cur_log_pr\n" if ($verbose);
		
		$last_log_pr = $cur_log_pr;
		$last_pr_first_base = $cur_pr_first_base;
		die if ($cur_pr_first_base < 0);

		$cur_pr_first_base -= 0.001;
		$cur_log_pr = _calculate_two_base_model_log10_likelihood($info_list, $error_rates, $first_base, $second_base, $cur_pr_first_base);
	}
	
	return ($last_log_pr, $last_pr_first_base);
}

sub _calculate_two_base_model_log10_likelihood
{
	my ($info_list, $error_rates, $first_base, $second_base, $pr_first_base) = @_;
	
	my $log10_likelihood = 0;	
	
	foreach my $item (@$info_list)
	{
		my $base_change_key_1 = $first_base . $item->{base};
	 	$base_change_key_1 = Breseq::Fastq::complement($base_change_key_1) if ($item->{strand} == -1); 
		my $base_change_key_2 = $second_base . $item->{base};
	 	$base_change_key_2 = Breseq::Fastq::complement($base_change_key_2) if ($item->{strand} == -1); 
		
		my $pr_ref_base_given_obs 
			= $pr_first_base * $error_rates->[$item->{fastq_file_index}]->{$item->{quality}}->{$base_change_key_1}
			+ (1-$pr_first_base) * $error_rates->[$item->{fastq_file_index}]->{$item->{quality}}->{$base_change_key_2};
		
		print "Probability of error is zero.\n  Fastq File = $item->{fastq_file_index}\n"
		  .  "  Quality = $item->{quality}\n  Error = $first_base$item->{base} $second_base$item->{base}\n"
		  if ($pr_ref_base_given_obs == 0);
			
		$log10_likelihood += log($pr_ref_base_given_obs);		
	}	
	$log10_likelihood /= log(10);
	return $log10_likelihood;
}


sub _best_strand_bias_model_log10_likelihood
{
## STUB
}

sub _calculate_strand_bias_model_log10_likelihood
{
	## strand_biased_error_rate is the chance of being in a new category with error probability=1 that only exists on $bias_strand
	my ($info_list, $error_rates, $first_base, $second_base, $bias_strand, $strand_biased_error_rate) = @_;
		
	my $log10_likelihood = 0;	
		
	foreach my $item (@$info_list)
	{		
		my $base_change_key = $first_base . $item->{base};
	 	$base_change_key = Breseq::Fastq::complement($base_change_key) if ($item->{strand} == -1); 

		my $strand_biased_change_key = $second_base . $item->{base};
	 	$strand_biased_change_key = Breseq::Fastq::complement($strand_biased_change_key) if ($item->{strand} == -1); 
		
		my $pr_obs = 0;
		if ($item->{strand} == $bias_strand)
		{						
			$pr_obs = $strand_biased_error_rate * $error_rates->[$item->{fastq_file_index}]->{$item->{quality}}->{$strand_biased_change_key}
				+ (1-$strand_biased_error_rate) * $error_rates->[$item->{fastq_file_index}]->{$item->{quality}}->{$base_change_key};
		}
		else
		{
			$pr_obs = $error_rates->[$item->{fastq_file_index}]->{$item->{quality}}->{$base_change_key};
		}
			
		$log10_likelihood += log($pr_obs);		
	}	
	$log10_likelihood /= log(10);
	return $log10_likelihood;
}


sub _calculate_two_base_strand_model_log10_likelihood
{
	my ($info_list, $error_rates, $first_base, $second_base, $strand_1_pr_first_base, $strand_2_pr_first_base) = @_;
	
	my $log10_likelihood = 0;	
	
	foreach my $item (@$info_list)
	{
		my $base_change_key_1 = $first_base . $item->{base};
	 	$base_change_key_1 = Breseq::Fastq::complement($base_change_key_1) if ($item->{strand} == -1); 
		my $base_change_key_2 = $second_base . $item->{base};
	 	$base_change_key_2 = Breseq::Fastq::complement($base_change_key_2) if ($item->{strand} == -1); 
		
		my $pr_first_base = ($item->{strand} == -1) ? $strand_1_pr_first_base : $strand_2_pr_first_base;
		
		my $pr_ref_base_given_obs 
			= $pr_first_base * $error_rates->[$item->{fastq_file_index}]->{$item->{quality}}->{$base_change_key_1}
			+ (1-$pr_first_base) * $error_rates->[$item->{fastq_file_index}]->{$item->{quality}}->{$base_change_key_2};
		
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

sub polymorphism_statistics
{
	our ($settings, $summary, $ref_seq_info) = @_;

	my $reference_fasta_file_name = $settings->file_name('reference_fasta_file_name');
	my $reference_bam_file_name = $settings->file_name('reference_bam_file_name');
	my $bam = Bio::DB::Sam->new(-fasta => $reference_fasta_file_name, -bam => $reference_bam_file_name);
	##my @seq_ids = $bam->seq_ids;
	my @seq_ids = @{$ref_seq_info->{seq_ids}};
		
	## some local variable lookups for convenience
	my $total_ref_length = 0;
	foreach my $seq_id (@seq_ids)
	{
		$total_ref_length+= $bam->length($seq_id);
	}
	
	my $log10_ref_length = log($total_ref_length) / log(10);	

	#for now, just average together all of the quality counts to get distribution
	my $max_quality = 0;
	my $error_counts_hash;
	foreach my $read_file ($settings->read_files)
	{
		my $fastq_file_index = $settings->read_file_to_fastq_file_index($read_file);
		my $this_error_counts_file_name = $settings->file_name('error_counts_file_name', {'#' => $read_file} );
		
		my $this_error_counts_hash = Breseq::ErrorCalibration::load_error_file($this_error_counts_file_name);
		
		foreach my $q (keys  %$this_error_counts_hash)
		{
			foreach my $b1 ('A', 'T', 'C', 'G')
			{
				foreach my $b2 ('A', 'T', 'C', 'G')
				{
					$error_counts_hash->{$q} += $this_error_counts_hash->{$q}->{"$b1$b2"};
				}
			}
			$max_quality = $q if ($q > $max_quality);
		}
	}

	my $genome_error_counts_file_name = $settings->file_name('genome_error_counts_file_name');

	open GEC, ">$genome_error_counts_file_name";
	my @counts;
	for (my $i=1; $i<=$max_quality; $i++)
	{
		my $val = 0;
		$val = $error_counts_hash->{$i} if (defined $error_counts_hash->{$i});
		print GEC "$val\n";
	}
	close GEC;

	my $polymorphism_statistics_input_file_name = $settings->file_name('polymorphism_statistics_input_file_name');
	my $polymorphism_statistics_output_file_name = $settings->file_name('polymorphism_statistics_output_file_name');

	### Load the older GenomeDiff and add new fields
	my $ra_mc_genome_diff_file_name = $settings->file_name('ra_mc_genome_diff_file_name');
	my $gd = Breseq::GenomeDiff->new(-in => $ra_mc_genome_diff_file_name);
	
	my $polymorphism_statistics_r_script_file_name = $settings->file_name('polymorphism_statistics_r_script_file_name');
	my $polymorphism_statistics_r_script_log_file_name = $settings->file_name('polymorphism_statistics_r_script_log_file_name');
	my $total_reference_length = $summary->{sequence_conversion}->{total_reference_sequence_length};
	
	Breseq::Shared::system("R --vanilla total_length=$total_reference_length in_file=$polymorphism_statistics_input_file_name out_file=$polymorphism_statistics_output_file_name qual_file=$genome_error_counts_file_name < $polymorphism_statistics_r_script_file_name > $polymorphism_statistics_r_script_log_file_name");
	
	## Read R file and add new results corresponding to all columns
	open ROUT, "<$polymorphism_statistics_output_file_name" or die "Could not find file: $polymorphism_statistics_output_file_name";
	my $header = <ROUT>;
	chomp $header;
	my @header_list = split /\t/, $header;
	
	my $new_gd = Breseq::GenomeDiff->new();
	foreach my $mut ($gd->list)
	{
		## lines only exist for RA evidence
		if ($mut->{type} ne 'RA')
		{
			$new_gd->add($mut);
			next;
		}
		
		## lines only exist for polymorphisms
		if (($mut->{frequency} == 1) || ($mut->{frequency} == 0))
		{
			$new_gd->add($mut);
			next;
		}

		my $line = <ROUT>;
		chomp $line;
		my @line_list = split /\t/, $line;
		
		for (my $i=0; $i< scalar @header_list; $i++)
		{
			$mut->{$header_list[$i]} = $line_list[$i];
			die "Incorrect number of items on line:\n$line" if (!defined $line_list[$i]);
		}

#		print Dumper($mut);

		
		# EXPERIMENTAL to not use bias_p_value
##		$mut->{quality} = $mut->{quality_genome_model};		
#		if ($mut->{quality} ne 'Inf')
#		{
#			$mut->{quality} = -(log ($mut->{quality}) / log(10)) - $log10_ref_length;
#			$new_gd->add($mut) if ($mut->{quality} > 2);
#		}	

		my $polymorphism_coverage_limit_both_bases = $settings->{polymorphism_coverage_both_bases};
		my $passed = 1;
		my ($top,$bot) = split /\//, $mut->{ref_cov};
		$passed &&= $top >= $polymorphism_coverage_limit_both_bases;
		$passed &&= $bot >= $polymorphism_coverage_limit_both_bases;		
		($top,$bot) = split /\//, $mut->{new_cov};
		$passed &&= $top >= $polymorphism_coverage_limit_both_bases;
		$passed &&= $bot >= $polymorphism_coverage_limit_both_bases;


		Breseq::GenomeDiff::add_reject_reason($mut, "POLYMORPHISM_STRAND") if (!$passed);
		Breseq::GenomeDiff::add_reject_reason($mut, "KS_QUALITY_P_VALUE_UNUSUAL_POLY") if ($mut->{ks_quality_p_value_unusual_poly} < $settings->{polymorphism_bias_p_value_cutoff});
#		Breseq::GenomeDiff::add_reject_reason($mut, "KS_QUALITY_P_VALUE_UNUSUAL_NEW") if ($mut->{ks_quality_p_value_unusual_new} < $settings->{polymorphism_bias_p_value_cutoff});
#		Breseq::GenomeDiff::add_reject_reason($mut, "KS_QUALITY_P_VALUE_UNUSUAL_REF") if ($mut->{ks_quality_p_value_unusual_ref} < $settings->{polymorphism_bias_p_value_cutoff});

#		Breseq::GenomeDiff::add_reject_reason($mut, "KS_QUALITY_P_VALUE_UNUSUAL_ALL") if ($mut->{ks_quality_p_value_unusual_all} < $settings->{polymorphism_bias_p_value_cutoff});
		Breseq::GenomeDiff::add_reject_reason($mut, "FISHER_STRAND_P_VALUE") if ($mut->{fisher_strand_p_value} < $settings->{polymorphism_bias_p_value_cutoff});

		if ($mut->{reject} && ($mut->{snp_quality} > $settings->{mutation_log10_e_value_cutoff}) && ($mut->{frequency} > 0.5) )
		{
			print Dumper($mut);
			$mut->{frequency} = 1;
			delete $mut->{reject};
		}

		$new_gd->add($mut);

		## END EXPERIMENTAL
	}
	
	### Write out the file which now has much more data
	my $polymorphism_statistics_ra_mc_genome_diff_file_name = $settings->file_name('polymorphism_statistics_ra_mc_genome_diff_file_name');
	$new_gd->write($polymorphism_statistics_ra_mc_genome_diff_file_name);
	
}

return 1;