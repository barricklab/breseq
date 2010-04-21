#!/usr/bin/env perl -w

###
# Pod Documentation
###

=head1 NAME

breseq.pl

=head1 SYNOPSIS

Usage: breseq.pl [-a] [-q solexa|phred] -r reference.gbk reads.fastq

Perform an analysis of bacterial resequencing data.

=head1 DESCRIPTION

=over

=item B<-a>

Print additional intermediate files.

=item B<-q> solexa|phred

Style of quality scores in fastq file. Must be either "solexa" (default) or "phred".
Generally, if the largest ASCII values you see are no greater than 70 (=F) then you have "phred" style scores.
If you have "solexa" style scores, then you probably see lowercase letters (=a). 

=item B<--snp-error-analysis>

Perform error analysis of SNPs. Do not predict mutations. Do not allow indels in read alignments. DEFAULT = OFF

=item B<--require-full-matches>

Matche entire reads. DEFAULT = OFF

=item B<--max-mismatched-per-read>

Disregard reads with more than this number of mismatches. DEFAULT = OFF

=item B<--no-mismatch-prediction>

Do not perform mismatch prediction.

=item B<--no-deletion-prediction>

Do not perform deletion prediction.

=item B<--no-junction-prediction>

Do not perform junction prediction.

=item B<--no-alignment-generation>

Do not generate alignments of reads to predicted mutations.

=item B<--alignment-read-limit> int

Only align this many reads (the first in the input files) to predicted mutations. FOR TESTING.

=back

=head1 AUTHOR

Jeffrey Barrick <jbarrick@msu.edu>

=head1 COPYRIGHT

Copyright 2008-2009.  All rights reserved.

=cut

###
# End Pod Documentation
###

=head1 TO-DO LIST

Need to print out reads that do not match or match unwanted sequences during candidate identification step.

Problem with display of reads overlapping new junction boundaries?

Need to have a step that disallows short read matches (when compiling alignment database)?

Correct reads that match junctions into two halves that can be used for real coverage and
on each half that they overlap. Eventually, get this information from alignment database and trim the remainder of the read
and give it a special tag to show that it matches a junction.


=cut


our $VERSION = '0.00_02';

#### Standard Perl Modules ####
use strict;
use Data::Dumper;
use Storable;
use POSIX qw(ceil floor);

use FindBin;
use lib $FindBin::Bin;
$ENV{PATH} = "$ENV{PATH}:" . $FindBin::Bin;

#### Breseq Perl Modules ####
use Breseq::AlignmentCorrection;
use Breseq::ErrorCalibration;
use Breseq::Output;
use Breseq::Settings;
use Breseq::Shared;
use Breseq::CandidateJunctions;
use Breseq::MutationIdentification;
use Breseq::ReferenceSequence;

#### BioPerl Modules ####
use Bio::SeqIO;

#### Configuration Options ####

our %unwanted_sequences = ( 
	'ILLUMINA_ADAPTOR'      => 'GATCGGAAGAGCTCGTATGCCGCTTCTGCTTCCGATC',  #Solexa Adaptor sequence.
);

## Keep a summary of certain statistics. 
my $summary = {};

###
### Get options from the command line
###    handles all GetOpt and filling in many other settings
### 
my $settings = Breseq::Settings->new;
Breseq::Output::record_time("Start");

##
# Convert the read fastq file to fasta for input into MUMmer
sub sequence_conversion {}
##

my $ref_seq_info;
my $sequence_converson_summary_file_name = $settings->file_name('sequence_conversion_summary_file_name');
my $sequence_conversion_done_file_name = $settings->file_name('sequence_conversion_done_file_name');
my $ref_seq_info_file_name = $settings->file_name('ref_seq_info_file_name');
if (!-e $sequence_conversion_done_file_name) 
{
	my $s;
	$settings->create_path('sequence_conversion_path');
	
	## quality trimming the reads
	if ($settings->{trim_reads})
	{
	 	print STDERR "  Trimming fastq reads for quality\n";
	 	foreach my $read_file ($settings->read_files)
		{
			my $fastq_file_name = $settings->read_file_to_fastq_file_name($read_file);	
			my $trimmed_fastq_file_name = $settings->file_name('trimmed_fastq_file_name', {'#'=>$read_file});
			Breseq::Fastq::fastq_to_trimmed_fastq($fastq_file_name, $trimmed_fastq_file_name, $settings->{quality_score_offset});
		}
	}
	
	##we need to know the maximum read length in each file before constructing junctions
	##and collect some information about the input read files at the same time
	print STDERR "  Analyzing fastq read files...\n";
	my $overall_max_read_length;
	foreach my $read_file ($settings->read_files)
	{
		my $max_read_length;
		my $total_bases = 0;
		my $num_reads = 0;
		my $fastq_file_name = $settings->read_file_to_fastq_file_name($read_file);	
		my $in = Breseq::Fastq->new(-file => $fastq_file_name);
		while (my $seq = $in->next_seq)
		{
			my $read_length = length $seq->{seq};
			$total_bases += $read_length;
			$num_reads++;
			$max_read_length = $read_length if ((!defined $max_read_length) || ($read_length > $max_read_length));
		}
		$s->{reads}->{$read_file}->{total_bases} = $total_bases;
		$s->{reads}->{$read_file}->{num_reads} = $num_reads;
		$s->{reads}->{$read_file}->{max_read_length} = $max_read_length;
		$overall_max_read_length = $max_read_length if ((!defined $overall_max_read_length) || ($max_read_length > $overall_max_read_length));
	}
	$s->{max_read_length} = $overall_max_read_length;
	$summary->{sequence_conversion} = $s;
	

	### The real way to do this is to add it to the FASTA file...
	### Create UNWANTED fasta sequence file.
	# if (!$settings->{no_filter_unwanted}) 
	# {
	# 	my $unwanted_sequences_fasta_file_name = $settings->file_name('unwanted_fasta_file_name');
	# 	open UNWANTED, ">$unwanted_sequences_fasta_file_name" or die "Could not open file $unwanted_sequences_fasta_file_name\n"; 
	# 	foreach my $key (keys %unwanted_sequences) 
	# 	{
	# 		print UNWANTED ">$settings->{unwanted_prefix}$key\n$unwanted_sequences{$key}\n";
	# 	}
	# 	close UNWANTED;
	# }
		
	## convert reference sequence to fasta and store other information so it can be reloaded quickly w/o bioperl
	$ref_seq_info = Breseq::Shared::process_reference_sequences($settings, $summary);
	Storable::store($ref_seq_info, $ref_seq_info_file_name) or die "Can't store data in file $ref_seq_info_file_name!\n";
	
	## store summary information
	Storable::store($summary->{sequence_conversion}, $sequence_converson_summary_file_name) or die "Can't store data in file $sequence_converson_summary_file_name!\n";
	
	open DONE, ">$sequence_conversion_done_file_name";
	close DONE;
	Breseq::Output::record_time("Sequence conversion");
}

#load this info
$ref_seq_info = Storable::retrieve($ref_seq_info_file_name);
die "Can't retrieve data from file $ref_seq_info_file_name!\n" if (!$ref_seq_info);
$summary->{sequence_conversion} = Storable::retrieve($sequence_converson_summary_file_name);
die "Can't retrieve data from file $sequence_converson_summary_file_name!\n" if (!$summary->{sequence_conversion});
(defined $summary->{sequence_conversion}->{max_read_length}) or die "Can't retrieve max read length from file $sequence_converson_summary_file_name!\n";
$settings->{max_read_length} = $summary->{sequence_conversion}->{max_read_length};
$settings->{total_reference_sequence_length} = $summary->{sequence_conversion}->{total_reference_sequence_length};

##
# Match all reads against the reference genome
sub alignment_to_reference {}
##

my $reference_alignment_done_file_name = $settings->file_name('reference_alignment_done_file_name');
if (!-e $reference_alignment_done_file_name) 
{
	print STDERR "Aligning reads to reference genome...\n";	
	$settings->create_path('reference_alignment_path');

#
#nohup ssaha2 -save REL606_kmer1_skip1 -kmer 10 -skip 1 -output sam_soft -outfile JEB559_REL606_kmer1_skip1.sam JEB559.fastq &> JEB559_REL606_kmer1_skip1.out
#nohup ssaha2 -save REL606_solexa -rtype solexa -output sam_soft -outfile JEB559_REL606_solexa.sam JEB559.fastq &> JEB559_REL606_solexa.out
#
#the lower seed value is important for finding split matches
#nohup ssaha2 -save REL606_solexa -kmer 13 -skip 2 -seeds 1 -score 12 -cmatch 9 -ckmer 6 -output sam_soft -outfile JEB559_REL606_solexa_3.sam JEB559.fastq >&JEB559_REL606_solexa_3.out &


	### create ssaha2 hash
	my $reference_hash_file_name = $settings->file_name('reference_hash_file_name');
	my $reference_fasta_file_name = $settings->file_name('reference_fasta_file_name');
	Breseq::Shared::system("ssaha2Build -rtype solexa -skip 1 -save $reference_hash_file_name $reference_fasta_file_name");		
	
	### ssaha2 align reads to reference sequences
	foreach my $read_struct ($settings->read_structures)
	{		
		##reads are paired
		if (defined $read_struct->{min_pair_dist} && defined $read_struct->{max_pair_dist})
		{
			die if (scalar @{$read_struct->{read_fastq_list}} != 2);
			
			my $fastq_1 = $read_struct->{read_fastq_list}->[0];
			my $fastq_2 = $read_struct->{read_fastq_list}->[1];
			my $min = $read_struct->{min_pair_dist};
			my $max = $read_struct->{max_pair_dist};
			
			my $reference_sam_file_name = $settings->file_name('reference_sam_file_name', {'#'=>$read_struct->{base_name}});	
			Breseq::Shared::system("ssaha2 -save $reference_hash_file_name -rtype solexa -ckmer 1  -skip 1 -cut 1000000000 -seeds 1 -output sam_soft -outfile $reference_sam_file_name -multi 1 -mthresh 9 -pair $min,$max $fastq_1 $fastq_2");
			
		}
		
		##reads are not paired
		else
		{			
			die if (scalar @{$read_struct->{base_names}} != 1);
			my $read_name = $read_struct->{base_names}->[0];
			my $read_fastq_file = $settings->read_file_to_fastq_file_name($read_name);
			my $reference_sam_file_name = $settings->file_name('reference_sam_file_name', {'#'=>$read_name});	
	#		Breseq::Shared::system("ssaha2 -save $reference_hash_file_name -disk 1 -score 9 -rtype solexa -skip 1 -cut 1000000000 -seeds 1 -output sam_soft -outfile $reference_sam_file_name $read_fastq_file");		
			Breseq::Shared::system("ssaha2 -save $reference_hash_file_name -rtype solexa -ckmer 1 -skip 1 -cut 1000000000 -seeds 1 -output sam_soft -outfile $reference_sam_file_name $read_fastq_file");		
		}
	}
	

	open DONE, ">$reference_alignment_done_file_name";
	close DONE;
	Breseq::Output::record_time("reference alignment");
}
else
{
	print STDERR "Reference alignment already exists.\n";
}

##
# Identify candidate junctions reads
sub identify_candidate_junctions {}
##


if (!$settings->{no_junction_prediction})
{
	my $candidate_junction_summary_file_name = $settings->file_name('candidate_junction_summary_file_name');
	
	if (!-e $settings->file_name('candidate_junction_done_file_name'))
	{
		print STDERR "Identifying candidate junctions...\n";
		$settings->create_path('candidate_junction_path');	
		Breseq::CandidateJunction::identify_candidate_junctions($settings, $summary, $ref_seq_info);
		Breseq::Output::record_time("Candidate junction identification");
		open DONE, ">" . $settings->file_name('candidate_junction_done_file_name');
		close DONE;

		Storable::store($summary->{candidate_junction}, $candidate_junction_summary_file_name) 
			or die "Can't store data in file $candidate_junction_summary_file_name!\n";

	}
	else
	{
		print STDERR "Identifying candidate junctions already complete.\n";
	}
	
	#load this info
	$summary->{candidate_junction} = Storable::retrieve($candidate_junction_summary_file_name);
	die "Can't retrieve data from file $candidate_junction_summary_file_name!\n" if (!$summary->{candidate_junction});
}
else
{
	print STDERR "Skipping identifying candidate junctions.\n";
}

		
##
# Find matches to new junction candidates
sub candidate_junction_alignment {}
##

if (!$settings->{no_junction_prediction})
{
	if (!-e $settings->file_name('candidate_junction_alignment_done_file_name'))
	{
		print STDERR "Finding matches to candidate junctions...\n";
		$settings->create_path('candidate_junction_alignment_path');	

		### create ssaha2 hash
		my $candidate_junction_hash_file_name = $settings->file_name('candidate_junction_hash_file_name');
		my $candidate_junction_fasta_file_name = $settings->file_name('candidate_junction_fasta_file_name');

		if (-s $candidate_junction_fasta_file_name > 0)
		{
			Breseq::Shared::system("ssaha2Build -rtype solexa -skip 1 -save $candidate_junction_hash_file_name $candidate_junction_fasta_file_name");		
		}
		
	
		### ssaha2 align reads to candidate junction sequences

		foreach my $read_name ($settings->read_files)
		{		
			my $candidate_junction_sam_file_name = $settings->file_name('candidate_junction_sam_file_name', {'#'=>$read_name});	
			my $read_fastq_file = $settings->read_file_to_fastq_file_name($read_name);
	#		Breseq::Shared::system("ssaha2 -save $candidate_junction_hash_file_name -disk 1 -rtype solexa -skip 1 -seeds 1 -output sam_soft -outfile $candidate_junction_sam_file_name $read_fastq_file");		
	# Note: Added -best parameter to try to avoid too many matches to redundant junctions!

			if (-e "$candidate_junction_hash_file_name.base")
			{
				Breseq::Shared::system("ssaha2 -save $candidate_junction_hash_file_name -best 1 -rtype solexa -skip 1 -seeds 1 -output sam_soft -outfile $candidate_junction_sam_file_name $read_fastq_file");		
			}
			else
			{
				open EMPTY, ">$candidate_junction_sam_file_name";
				close EMPTY;
			}
		}
	
	
		open DONE, ">" . $settings->file_name('candidate_junction_alignment_done_file_name');
		close DONE;
	}
	else
	{
		print STDERR "Finding matches to candidate junctions already complete.\n";
	}
}
else
{
	print STDERR "Skipping candidate junction alignment.\n";
}

##
# Resolve matches to new junction candidates
sub alignment_correction {}
##

my @hybrids;
my $alignment_correction_done_file_name = $settings->file_name('alignment_correction_done_file_name');
if (!-e $alignment_correction_done_file_name)
{
	print STDERR "Resolving alignments with candidate junctions and ambiguous ends...\n";
	$settings->create_path('alignment_correction_path');		
	my $predicted_junctions_file_name = $settings->file_name('predicted_junction_file_name');
	@hybrids = Breseq::AlignmentCorrection::correct_alignments($settings, $summary, $ref_seq_info);
	Storable::store(\@hybrids, $predicted_junctions_file_name) or die "Can't store data in file $predicted_junctions_file_name!\n";
	my $alignment_correction_summary_file_name = $settings->file_name('alignment_correction_summary_file_name');	
	Storable::store($summary->{alignment_correction}, $alignment_correction_summary_file_name) 
		or die "Can't store data in file $alignment_correction_summary_file_name!\n";	
	Breseq::Output::record_time("Resolve candidate junctions");
	open DONE, ">$alignment_correction_done_file_name";
	close DONE;
}
else
{
	print STDERR "Resolving candidate junctions/ambiguous ends already complete.\n";
	my $predicted_junctions_file_name = $settings->file_name('predicted_junction_file_name');	
	@hybrids = ();
	@hybrids = @{Storable::retrieve($predicted_junctions_file_name)} if (!$settings->{no_junction_prediction});
	my $alignment_correction_summary_file_name = $settings->file_name('alignment_correction_summary_file_name');	
	$summary->{alignment_correction} = Storable::retrieve($alignment_correction_summary_file_name) if (-e $alignment_correction_summary_file_name);
}	


##
# Create BAM files
sub bam_creation {}
##

my $bam_done_file_name = $settings->file_name('bam_done_file_name');
if (!-e $bam_done_file_name)
{
	print STDERR "Creating BAM files...\n";
	
	$settings->create_path('bam_path');		
	
	my $reference_faidx_file_name = $settings->file_name('reference_faidx_file_name');
	my $candidate_junction_faidx_file_name = $settings->file_name('candidate_junction_faidx_file_name');

	my $resolved_junction_sam_file_name = $settings->file_name('resolved_junction_sam_file_name');
	my $junction_bam_unsorted_file_name = $settings->file_name('junction_bam_unsorted_file_name');
	my $junction_bam_prefix = $settings->file_name('junction_bam_prefix');	
	my $junction_bam_file_name = $settings->file_name('junction_bam_file_name');

	if (!$settings->{no_junction_prediction})
	{
		Breseq::Shared::system("samtools import $candidate_junction_faidx_file_name $resolved_junction_sam_file_name $junction_bam_unsorted_file_name");
		Breseq::Shared::system("samtools sort $junction_bam_unsorted_file_name $junction_bam_prefix");
		Breseq::Shared::system("samtools index $junction_bam_file_name");
	}
	
	my $resolved_reference_sam_file_name = $settings->file_name('resolved_reference_sam_file_name');
	my $reference_bam_unsorted_file_name = $settings->file_name('reference_bam_unsorted_file_name');
	my $reference_bam_prefix = $settings->file_name('reference_bam_prefix');
	my $reference_bam_file_name = $settings->file_name('reference_bam_file_name');
	
	Breseq::Shared::system("samtools import $reference_faidx_file_name $resolved_reference_sam_file_name $reference_bam_unsorted_file_name");
	Breseq::Shared::system("samtools sort $reference_bam_unsorted_file_name $reference_bam_prefix");
	Breseq::Shared::system("samtools index $reference_bam_file_name");
		
	#Might be better to check for existence of these output files manually, as samtools does not always return an error code.
	
	Breseq::Output::record_time("BAM file creation");
	open DONE, ">$bam_done_file_name";
	close DONE;
}
else
{
	print STDERR "BAM file creation already complete.\n";
}



sub tabulate_pairs_in_tam {}

{
	my @rs = $settings->read_structures;
		
	my @min_pair_dist;
	my @max_pair_dist;
	
	my $paired = 0;
	
	my $i=0;
	foreach my $rfi (@{$settings->{read_file_index_to_struct_index}})
	{
		$min_pair_dist[$i] = 0;
		$max_pair_dist[$i] = 0;
		
		if ($rs[$rfi]->{paired})
		{
			$paired = 1;
			$min_pair_dist[$i] = $rs[$rfi]->{min_pair_dist};
			$max_pair_dist[$i] = $rs[$rfi]->{max_pair_dist};
		}
		$i++;
	}
	
	my $long_pairs_file_name = $settings->file_name('long_pairs_file_name');
	
	if ($paired && (!-e $long_pairs_file_name))
	{
	
		my $reference_sam_file_name = $settings->file_name('resolved_reference_sam_file_name');
		my $reference_tam = Bio::DB::Tam->open($reference_sam_file_name) or die "Could not open $reference_sam_file_name";

		my $reference_faidx_file_name = $settings->file_name('reference_faidx_file_name');
		my $reference_header = $reference_tam->header_read2($reference_faidx_file_name) or throw("Error reading reference fasta index file: $reference_faidx_file_name");		
		my $target_names = $reference_header->target_name;

		my $save;
		my $on_alignment = 0;
		my $last;		
	
		while (1)
		{
			$a = Bio::DB::Bam::Alignment->new();
			my $bytes = $reference_tam->read1($reference_header, $a);
			last if ($bytes <= 0);
		
		
			my $start       = $a->start;
		    my $end         = $a->end;
		    my $seqid       = $target_names->[$a->tid];
		
			$on_alignment++;
			print "$on_alignment\n" if ($on_alignment % 10000 == 0);
		
			#last if ($on_alignment > 100000);
			
			#print $a->qname . "\n";		
				
			if (!$a->munmapped)
			{
				my $mate_insert_size = abs($a->isize);
				my $mate_end = $a->mate_end;
				my $mate_start = $a->mate_start;
				my $mate_reversed = 2*$a->mreversed + $a->reversed;
		 		my $mreversed = $a->mreversed;
		 		my $reversed = $a->reversed;	
	
				my $fastq_file_index = $a->aux_get('X2');
				#print "$mate_insert_size $min_pair_dist[$fastq_file_index] $max_pair_dist[$fastq_file_index]\n";
				#if (($mate_insert_size < $min_pair_dist[$fastq_file_index]) || ($mate_insert_size > $max_pair_dist[$fastq_file_index]))
				if ((($mate_insert_size >= 400) && ($mate_insert_size < $min_pair_dist[$fastq_file_index])) || ($mate_insert_size > $max_pair_dist[$fastq_file_index]))
				{
					#correct pair
				
					if ($last && ($last->{start} == $mate_start))
					{					
						$save->{int($start/100)}->{int($mate_start/100)}->{$mate_reversed}++;
						$save->{int($last->{start}/100)}->{int($last->{mate_start}/100)}->{$last->{mate_reversed}}++;
						undef $last;				
					}
					else
					{
						($last->{mate_reversed}, $last->{start}, $last->{mate_start}) = ($mate_reversed, $start, $mate_start);
					}
								
					#$save->{$mate_reversed}->{int($start/100)}->{int($mate_start/100)}++;
				    #print $a->qname," aligns to $seqid:$start..$end, $mate_start $mate_reversed ($mreversed $reversed) $mate_insert_size\n";	
				}

			}
		}
	
		open LP, ">$long_pairs_file_name" or die;

		foreach my $key_1 (sort {$a <=> $b} keys %$save)
		{
			foreach my $key_2 (sort {$a <=> $b} keys %{$save->{$key_1}})
			{
				foreach my $key_reversed (sort {$a <=> $b} keys %{$save->{$key_1}->{$key_2}})
				{
					print LP "$key_1\t$key_2\t$key_reversed\t$save->{$key_1}->{$key_2}->{$key_reversed}\n";
				}
			}
		}
		close LP;
	}
	
	if ($paired)
	{
		open LP, "$long_pairs_file_name" or die;
		while ($_ = <LP>)
		{
			chomp $_;
			my ($start, $end, $key_reversed);
		}
	}
}





=comment


sub tabulate_pairs {}

{
	my $reference_bam_file_name = $settings->file_name('reference_bam_file_name');

	my $index = Bio::DB::Bam->index_open($reference_bam_file_name);
	my $bam = Bio::DB::Bam->open($reference_bam_file_name);
	my $header = $bam->header;
	my $target_names = $header->target_name;
	my @rs = $settings->read_structures;
		
	my @min_pair_dist;
	my @max_pair_dist;
	
	my $i=0;
	foreach my $rfi (@{$settings->{read_file_index_to_struct_index}})
	{
		$min_pair_dist[$i] = $rs[$rfi]->{min_pair_dist};
		$max_pair_dist[$i] = $rs[$rfi]->{max_pair_dist};
		$i++;
	}
		
	my $save;
	
	my $on_alignment;
	
	my $long_pairs_file_name = $settings->file_name('long_pairs_file_name');
	open LP, ">$long_pairs_file_name";
	my $callback = sub {
	    my $a = shift;
	    my $start       = $a->start;
	    my $end         = $a->end;
	    my $seqid       = $target_names->[$a->tid];
		
		$on_alignment++;
		print "$on_alignment\n" if ($on_alignment % 10000 == 0);
		
		if (!$a->munmapped)
		{
			my $mate_insert_size = abs($a->isize);
		#	my $mate_end = $a->mate_end;
			my $mate_start = $a->mate_start;
			my $mate_reversed = ($a->mreversed + $a->reversed) % 2;
	 		my $mreversed = $a->mreversed;
	 		my $reversed = $a->reversed;
	
		    #print $a->qname," aligns to $seqid:$start..$end, $mate_start $mreversed $reversed\n";	
	
	
			my $fastq_file_index = $a->aux_get('X2');
			#print "$mate_insert_size $min_pair_dist[$fastq_file_index] $max_pair_dist[$fastq_file_index]\n";
			#if (($mate_insert_size < $min_pair_dist[$fastq_file_index]) || ($mate_insert_size > $max_pair_dist[$fastq_file_index]))
			
			if ($mate_insert_size > $max_pair_dist[$fastq_file_index])
			{
				$save->{$mate_reversed}->{int($start/100)}->{int($mate_start/100)}++;
			    #print $a->qname," aligns to $seqid:$start..$end, $mate_start $mate_reversed ($mreversed $reversed) $mate_insert_size\n";	
			}
		}
	};
	
	for (my $i=0; $i<$header->n_targets; $i++)
	{
		$index->fetch($bam,$i,1,$header->target_len->[$i],$callback);
	}
	
	
	foreach my $key_reversed (sort {$a <=> $b} keys %$save)
	{	
		foreach my $key_1 (sort {$a <=> $b} keys %{$save->{$key_reversed}})
		{
			foreach my $key_2 (sort {$a <=> $b} keys %{$save->{$key_reversed}->{$key_1}})
			{
				print LP "$key_1\t$key_2\t$key_reversed\t$save->{$key_reversed}->{$key_1}->{$key_2}\n";
			}
		}
	}
	
	
	close LP;
}

=cut


##
# Tabulate error counts and coverage distribution at unique only sites
sub error_count {}
##

my $error_counts_done_file_name = $settings->file_name('error_counts_done_file_name');
if (!-e $error_counts_done_file_name) 
{
	$settings->create_path('error_calibration_path');
	print STDERR "Tabulating error counts...\n";
	Breseq::ErrorCalibration::count($settings, $summary, $ref_seq_info);
	open DONE, ">$error_counts_done_file_name";
	close DONE;
}
else 
{
	print STDERR "Tabulating error counts complete.\n";
}


##
# Calculate error rates
sub error_rates {}
##

my $error_rates_done_file_name = $settings->file_name('error_rates_done_file_name');
my $error_rates_summary_file_name = $settings->file_name('error_rates_summary_file_name');
if (!-e $error_rates_done_file_name) 
{
	$settings->create_path('output_path'); ##need output for plots	
	print STDERR "Calibrating error rates...\n";
	Breseq::ErrorCalibration::error_counts_to_error_rates($settings, $summary, $ref_seq_info);
	Breseq::ErrorCalibration::analyze_unique_coverage_distributions($settings, $summary, $ref_seq_info);
	
	Storable::store($summary->{unique_coverage}, $error_rates_summary_file_name) or die "Can't store data in file $error_rates_summary_file_name!\n";
	open DONE, ">$error_rates_done_file_name";
	close DONE;
}
else 
{
	print STDERR "Calibrating error model complete.\n";
}
$summary->{unique_coverage} = Storable::retrieve($error_rates_summary_file_name);
die "Can't retrieve data from file $error_rates_summary_file_name!\n" if (!$summary->{unique_coverage});
#these are determined by the loaded summary information
$settings->{unique_coverage} = $summary->{unique_coverage};

##
# Make predictions of point mutations, small indels, and large deletions
sub mutation_prediction {}
##

my @mutations;
my @deletions;
my @unknowns;

if (!$settings->{no_mutation_prediction}) #could remove conditional?
{		
	print STDERR "Collecting and evaluating SNPs...\n";
	
	$settings->create_path('mutation_identification_path');
	
	##
	# Load error rates
	##
	my $error_rates = Breseq::ErrorCalibration::load_error_rates($settings, $summary, $ref_seq_info);
	
	##
	# Predict SNPS, indels within reads, and large deletions
	#   this function handles all file creation...
	##
	my $mutation_info =  Breseq::MutationIdentification::identify_mutations($settings, $summary, $error_rates);

	@mutations = @{$mutation_info->{mutations}};
	@deletions = @{$mutation_info->{deletions}};
	@unknowns = @{$mutation_info->{unknowns}};
}





my $output_done_file_name = $settings->file_name('output_done_file_name');	
if (!-e $output_done_file_name)
{
	###
	## Output Genome Diff File
	sub genome_diff_output {}
	###
	print STDERR "Creating genome diff file...\n";

	### filter what hybrids we believe in...
	### this should also annotate and do much more of the filtering that is currently done
	### during the prediction step
	foreach my $hybrid (@hybrids)
	{
		my $coverage_cutoff_1 = $settings->{unique_coverage}->{$hybrid->{interval_1}->{seq_id}}->{deletion_coverage_propagation_cutoff};
		my $coverage_cutoff_2 = $settings->{unique_coverage}->{$hybrid->{interval_2}->{seq_id}}->{deletion_coverage_propagation_cutoff};
		$hybrid->{marginal} = 1 if (($hybrid->{total_reads} < $coverage_cutoff_1) &&  ($hybrid->{total_reads} < $coverage_cutoff_2));
	}
	
	my $full_genome_diff_file_name = $settings->file_name('full_genome_diff_file_name');
	my $filtered_genome_diff_file_name = $settings->file_name('filtered_genome_diff_file_name');
	my $marginal_genome_diff_file_name = $settings->file_name('marginal_genome_diff_file_name');
	
	Breseq::Output::write_genome_diff($full_genome_diff_file_name, $settings, \@mutations, \@deletions, \@hybrids, \@unknowns);
	my $conditions = "marginal!=1";
	Breseq::Shared::system("gd_utils.pl FILTER -i $full_genome_diff_file_name -o $filtered_genome_diff_file_name -r $marginal_genome_diff_file_name $conditions");


	###
	## Below this point is optional and has heftier prerequisites
	## for running, including R and BioPerl
	## 
	###
	if (1)
	{
		my $marginal_html_file_name = $settings->file_name('marginal_html_file_name');
		
		my $res;
		$res = Breseq::Shared::system("gd_utils.pl ANNOTATE -i $filtered_genome_diff_file_name $settings->{arguments} ", 0, 1);
		my $annotated_mutations = ($res == 0);
		$res = Breseq::Shared::system("gd_utils.pl ANNOTATE -i $marginal_genome_diff_file_name --html-mutation-file=$marginal_html_file_name --sort_junctions_by_score $settings->{arguments} ", 0, 1);
		my $annotated_marginal = ($res == 0);

		my $index_html_file_name = $settings->file_name('index_html_file_name');	
		Breseq::Output::html_index($index_html_file_name, $settings, $summary, $ref_seq_info, $annotated_mutations, $annotated_marginal);

		###
		## Temporary debug output using Data::Dumper
		###

		my $summary_text_file_name = $settings->file_name('summary_text_file_name');
		open SUM, ">$summary_text_file_name";
		print SUM Dumper($summary);
		close SUM;
	
		my $settings_text_file_name = $settings->file_name('settings_text_file_name');
		open SETTINGS, ">$settings_text_file_name";
		print SETTINGS Dumper($settings);
		close SETTINGS;
		
		###
		## Temporary debug output using Data::Dumper
		###
	}

	## record the final time and print summary table
	Breseq::Output::record_time("End");
	Breseq::Output::html_summary_table($settings->{summary_html_file_name}, $settings, $summary, \@Breseq::Output::execution_times);

	my $output_done_file_name = $settings->file_name('output_done_file_name');	
	open DONE, ">$output_done_file_name";
	close DONE;
}
