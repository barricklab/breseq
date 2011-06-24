###
# Pod Documentation
###

=head1 NAME

Breseq::Settings

=head1 SYNOPSIS

Perl modules used internally by breseq.

=head1 AUTHOR

Jeffrey E. Barrick <jeffrey.e.barrick@gmail.com>

=head1 LICENSE AND COPYRIGHT

Copyright (c) 2010 Michigan State University

breseq is free software; you can redistribute it and/or modify it under the terms the 
GNU General Public License as published by the Free Software Foundation; either 
version 1, or (at your option) any later version.

=cut

###
# End Pod Documentation
###

package Breseq::Settings;
use strict;
use Getopt::Long;
use Data::Dumper;
use Pod::Usage;
use FindBin;
use Breseq::File::Path;

use Breseq;

use vars qw(@ISA);
@ISA = qw( Bio::Root::Root );

our %unwanted_sequences = ( 
	'UNWANTED::ILLUMINA_ADAPTOR_1'    => 'GATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG',  #Solexa Adaptor sequence.
	'UNWANTED::ILLUMINA_ADAPTOR_2'	  => 'ACACTCTTTCCCTACACGACGCTCTTCCGATCT'
);                              

## Options for turning analysis off ##
## Mainly for development, long names only ##

sub new
{	
	my($caller,@args) = @_;
	my $class = ref($caller) || $caller;
	my $self = new Bio::Root::Root($caller, @args);
	bless ($self, $class);

	$self->pre_option_initialize;
	
	#keep this on by default
	
	my ($help, $man);
	GetOptions(
		'help|?' => \$help, 'man' => \$man,
		'verbose|v' => \$self->{verbose},
	## Options for turning various analysis chunks off or on
		'no-junction-prediction' => \$self->{no_junction_prediction},
		'no-mismatch-prediction' => \$self->{no_mutation_prediction},
		'no-deletion-prediction' => \$self->{no_deletion_prediction},
		'no-alignment-generation' => \$self->{no_alignment_generation},
		'no-filter-unwanted' => \$self->{no_filter_unwanted},
		'no-unmatched-reads' => \$self->{no_unmatched_reads},
		'copy-number-variation' => \$self->{copy_number_variation},
	## Options for only using part of the data		
		'read-limit|l=s' => \$self->{read_limit},
		'candidate-junction-read-limit=s' => \$self->{candidate_junction_read_limit},
		'alignment-read-limit=s' => \$self->{alignment_read_limit},		
	## Options for input
		'name|n=s' => \$self->{run_name},	
		'output-path|o=s' => \$self->{base_output_path},	
		'reference-sequence|r=s' => \@{$self->{reference_genbank_file_names}},
		'junction-sequence|j=s' => \@{$self->{junction_only_reference_genbank_file_names}},	
	## Options for output			
		'keep-all-intermediates' => \$self->{keep_all_intermediates},
		'shade-frequencies' => \$self->{shade_frequencies},
		'max-rejected-polymorphisms-to-show=s' => \$self->{max_rejected_polymorphisms_to_show},
		'max-rejected-junctions-to-show=s' => \$self->{max_rejected_junctions_to_show},	
	## Options for read aligment
		'maximum-mismatches|m=s' => \$self->{maximum_read_mismatches},	
		'required-match-length=s' => \$self->{required_match_length},
	## Options for snp error analysis
		'require-complete-match' => \$self->{require_complete_match},
		'require-no-indel-match' => \$self->{require_no_indel_match},
		'require-max-mismatches=s' => \$self->{require_max_mismatches},
		'do-not-trim-ambiguous-ends' => \$self->{do_not_trim_ambiguous_ends},
		'base-quality-cutoff|b=s' => \$self->{base_quality_cutoff},
		'error-model-method=s' => \$self->{error_model_method},
	## Options for polymorphism analysis
		'polymorphism-prediction' => \$self->{polymorphism_prediction},		
		'polymorphism-log10-e-value-cutoff=s' => \$self->{polymorphism_log10_e_value_cutoff},
		'polymorphism-frequency-cutoff=s' =>  \$self->{polymorphism_frequency_cutoff},	
		'polymorphism-p-value-cutoff=s' => \$self->{polymorphism_bias_p_value_cutoff},
		'polymorphism-coverage-both-strands=s' => \$self->{polymorphism_coverage_both_strands},	
		'polymorphism-reject-homopolymer-length=s' => \$self->{polymorphism_reject_homopolymer_length},
		'no-indel-polymorphisms' => \$self->{no_indel_polymorphisms},	
	## Options for candidate junction identification	
		'maximum-candidate-junctions=s' => \$self->{maximum_candidate_junctions},
	## Options for using deprecated and slow Perl methods		
		'perl-error-count' => \$self->{perl_error_count},
		'perl-identify-mutations' => \$self->{perl_identify_mutations},
		'perl-calc-trims' => \$self->{perl_calc_trims},
		'strict-polymorphism-prediction' => \$self->{strict_polymorphism_prediction},
		'perl-preprocess-alignments' => \$self->{perl_preprocess_alignments},
		'perl-identify-candidate-junctions' => \$self->{perl_identify_candidate_junctions},
##		'smalt' => \$self->{smalt},
		'no_ssaha2' => \$self->{no_ssaha2},
		
	) or pod2usage(2);

	pod2usage(1) if $help;
	pod2usage(-exitstatus => 0, -verbose => 2) if $man;
	pod2usage(-exitstatus => 0, -verbose => 2) if (scalar @ARGV == 0);
	
	## default to using smalt
	$self->{smalt} = 1 if (!$self->{no_ssaha2});
	
	$self->post_option_initialize;
	
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
	my $time_stamp = sprintf "%4d-%02d-%02d %02d:%02d:%02d",$year+1900,$mon+1,$mday,$hour,$min,$sec;
	$self->log($time_stamp);	
	$self->log($self->{full_command_line});
	
	return $self;
}

## called before getting options from command line
sub pre_option_initialize
{
	my ($self) = @_;

	@{$self->{reference_genbank_file_names}} = ();  # files containing reference sequences	
	@{$self->{junction_only_reference_genbank_file_names}} = (); #files to look for junctions to but not align to
	
	## Set up default values for options
	$self->{full_command_line} = "$0 @ARGV"; 
	$self->{arguments} = "@ARGV";	
	$self->{predicted_quality_type} = '';
	$self->{min_quality} = 0;
	$self->{max_quality} = 0;
	$self->{run_name} = 'unnamed';
	$self->{clean} = 0;
	$self->{base_output_path} = '';
	$self->{error_model_method} = 'EMPIRICAL';
	$self->{base_quality_cutoff} = 3;			# avoids problem with Illumina assigning 2 to bad ends of reads!
	
	###
	#   CandidateJunctions.pm
	###
	$self->{no_junction_prediction} = undef;					# don't perform junction prediction steps
	$self->{candidate_junction_score_method} = 'POS_HASH'; 		# Either POS_HASH, or MIN_OVERLAP
	$self->{preprocess_junction_min_indel_split_length} = 3;	# Split the SAM entries on indels of this many bp or more before identifying CJ
																# Undefined = OFF
																
	#### Scoring to decide which pairs of alignments to the same read to consider
	$self->{required_extra_pair_total_length} = 2;				# The union of the pairs must exceed the best of any single match by this length
																# Setting this does penalize some *real* new junctions, so be careful!
	$self->{required_both_unique_length_per_side} = 5;			# Require both of the pair of matches supporting a junction to have this
																	# much of their matches unique in the reference sequence.
	$self->{required_one_unique_length_per_side} = 10;			# Require at least one of the pair of matches supporting a junction to have this
																	# much of its match that is unique in the reference sequence.

	#### Scoring section to choose which ones from list to take
	$self->{minimum_candidate_junction_pos_hash_score} = 0;		# Require at least this many unique start coordinate/strand reads to accept a CJ
																# OFF by default, because a fixed number are taken
	$self->{minimum_candidate_junction_min_overlap_score} = 0;	# Require at least this many unique start coordinate/strand reads to accept a CJ
																# OFF by default, because a fixed number are taken

#	$self->{required_unique_length_per_side} = 10;				# Require at least one of the pair of matches supporting a junction to have this
																# much of its match that is unique in the reference sequence.
	$self->{maximum_inserted_junction_sequence_length} = 20;	# Ignore junctions with negative overlap (unique inserted sequence between reference 
																# matches) greater than this length. Prevents evidence file names that are too long.
	$self->{minimum_candidate_junctions} = 200;					# Minimum number of candidate junctions to keep
	$self->{maximum_candidate_junctions} = 5000;				# Maximum number of candidate junctions to keep
	$self->{maximum_candidate_junction_length_factor} = 0.1;	# Only keep CJ cumulative lengths adding up to this factor times the total reference size 
	$self->{candidate_junction_read_limit} = undef;		    	# FOR TESTING: only process this many reads when creating candidate junctions

	#used by AlignmentCorrection.pm
	$self->{add_split_junction_sides} = 1;			# Add the sides of passed junctions to the SAM file?
	$self->{required_match_length} = 28;			# Match must span this many bases in query to count as a match
	
	$self->{no_mutation_prediction} = undef;		# don't perform read mismatch/indel prediction steps
	$self->{no_deletion_prediction} = undef;		# don't perform deletion prediction steps
	$self->{no_alignment_generation} = undef;		# don't generate alignments
	$self->{alignment_read_limit} = undef;			# only go through this many reads when creating alignments
	$self->{correction_read_limit} = undef;			# only go through this many reads when correcting alignments
	
	## NOT IMPLEMENTED
	$self->{no_filter_unwanted} = undef;			# don't filter out unwanted reads with adaptor matches
	$self->{unwanted_prefix} = "UNWANTED:::";	# prefix on unwanted read names
	
	#used by MutationIdentification.pm
	$self->{mutation_log10_e_value_cutoff} = 2;		   # log10 of evidence required for SNP calls 
	$self->{polymorphism_log10_e_value_cutoff} = 2;
	$self->{polymorphism_bias_p_value_cutoff} = 0.05;
	$self->{polymorphism_frequency_cutoff} = 0;        # cut off if < X or > 1-X
	$self->{polymorphism_coverage_both_strands} = 0;     # require this much coverage on each strand
	$self->{no_indel_polymorphisms} = 0;		
	
	#used by Output.pm
	$self->{max_rejected_polymorphisms_to_show} = 20;
	$self->{max_rejected_junctions_to_show} = 20;
	$self->{hide_circular_genome_junctions} = 1;		
		
	 @{$self->{execution_times}} = ();
}

## called after getting options from command line
sub post_option_initialize
{
	my ($self) = @_;
	
	$self->{version} = $Breseq::VERSION;
	$self->{byline} = "<b><i>breseq</i></b>&nbsp;&nbsp;version $self->{version}";
	$self->{website} = "http://barricklab.org/breseq";
	$self->{bin_path} = $FindBin::Bin;
	$self->{lib_path} = "$self->{bin_path}/../lib/perl5/Breseq";

	#neaten up some settings for later string comparisons
	$self->{error_model_method} = "\U$self->{error_model_method}";

	#on by default
	$self->{unmatched_reads} = ($self->{no_unmatched_reads}) ? 0 : 1;
		
	## block option
	if ($self->{strict_polymorphism_prediction})
	{
		$self->{polymorphism_prediction} = 1;
		$self->{maximum_read_mismatches} = 1;
		$self->{require_complete_match} = 1;
		$self->{no_indel_polymorphisms} = 1;
		$self->{polymorphism_log10_e_value_cutoff} = 5;
	}
	
	## problems if there are spaces b/c shell removes quotes before we know about them
	## thus require run names to only use underscores (but when printing output, remove).
	if ($self->{run_name} =~ m/\s/)
	{
		die("Option <--name|-n> must not contain whitespace characters. Please use underscores '_' in place of spaces.\n");
	}
	$self->{print_run_name} = $self->{run_name};
	$self->{print_run_name} =~ s/_/ /g;
	
	#######  SETUP FILE NAMES  #######
	### '#' replaced with read fastq name
	### '@' replaced by seq_id of reference sequence

	##### sequence conversion #####
	$self->{sequence_conversion_path} = "01_sequence_conversion";
	$self->{sequence_conversion_path} = "$self->{base_output_path}/$self->{sequence_conversion_path}" if ($self->{base_output_path});
	$self->{converted_fastq_file_name} = "$self->{sequence_conversion_path}/#.converted.fastq";
    $self->{ref_seq_info_file_name} = "$self->{sequence_conversion_path}/ref_seq_info.bin";	
	$self->{unwanted_fasta_file_name} = "$self->{sequence_conversion_path}/unwanted.fasta";
	$self->{reference_trim_file_name} = "$self->{sequence_conversion_path}/@.trims";
	$self->{sequence_conversion_summary_file_name} = "$self->{sequence_conversion_path}/summary.bin";
	$self->{sequence_conversion_done_file_name} = "$self->{sequence_conversion_path}/sequence_conversion.done";

	##### reference #####
	$self->{reference_alignment_path} = "02_reference_alignment";
	$self->{reference_alignment_path} = "$self->{base_output_path}/$self->{reference_alignment_path}" if ($self->{base_output_path});
	$self->{reference_hash_file_name} = "$self->{reference_alignment_path}/reference";
	$self->{reference_sam_file_name} = "$self->{reference_alignment_path}/#.reference.sam";
	$self->{reference_alignment_done_file_name} = "$self->{reference_alignment_path}/alignment.done";
	
	##### candidate junction #####
	$self->{candidate_junction_path} = "03_candidate_junctions";
	$self->{candidate_junction_path} = "$self->{base_output_path}/$self->{candidate_junction_path}" if ($self->{base_output_path});

	$self->{preprocess_junction_best_sam_file_name} = "$self->{candidate_junction_path}/best.sam";
	$self->{preprocess_junction_split_sam_file_name} = "$self->{candidate_junction_path}/#.split.sam";
	$self->{preprocess_junction_done_file_name} = "$self->{candidate_junction_path}/preprocess_junction_alignment.done";

	$self->{coverage_junction_best_bam_unsorted_file_name} = "$self->{candidate_junction_path}/best.unsorted.bam";
	$self->{coverage_junction_best_bam_file_name} = "$self->{candidate_junction_path}/best.bam";
	$self->{coverage_junction_best_bam_prefix} = "$self->{candidate_junction_path}/best";
	$self->{coverage_junction_distribution_file_name} = "$self->{candidate_junction_path}/@.unique_only_coverage_distribution.tab";
	$self->{coverage_junction_plot_file_name} = "$self->{candidate_junction_path}/@.coverage.pdf";
	$self->{coverage_junction_summary_file_name} = "$self->{candidate_junction_path}/coverage.summary.bin";
	$self->{coverage_junction_done_file_name} = "$self->{candidate_junction_path}/coverage_junction_alignment.done";

	$self->{candidate_junction_summary_file_name} = "$self->{candidate_junction_path}/candidate_junction_summary.bin";	
	$self->{candidate_junction_fasta_file_name} = "$self->{candidate_junction_path}/candidate_junction.fasta";
	$self->{candidate_junction_faidx_file_name} = "$self->{candidate_junction_path}/candidate_junction.fasta.fai";	
	$self->{candidate_junction_done_file_name} = "$self->{candidate_junction_path}/candidate_junction.done";
	
	##### candidate junction alignment #####
	$self->{candidate_junction_alignment_path} = "04_candidate_junction_alignment";
	$self->{candidate_junction_alignment_path} = "$self->{base_output_path}/$self->{candidate_junction_alignment_path}" if ($self->{base_output_path});
	$self->{candidate_junction_hash_file_name} = "$self->{candidate_junction_alignment_path}/candidate_junction";
	$self->{candidate_junction_sam_file_name} = "$self->{candidate_junction_alignment_path}/#.candidate_junction.sam";
	$self->{candidate_junction_alignment_done_file_name} = "$self->{candidate_junction_alignment_path}/candidate_junction_alignment.done";

	##### alignment correction #####
	$self->{alignment_correction_path} = "05_alignment_correction";
	$self->{alignment_correction_path} = "$self->{base_output_path}/$self->{alignment_correction_path}" if ($self->{base_output_path});
	$self->{resolved_reference_sam_file_name} = "$self->{alignment_correction_path}/reference.sam";
	$self->{resolved_junction_sam_file_name} = "$self->{alignment_correction_path}/junction.sam";
	$self->{alignment_correction_summary_file_name} = "$self->{alignment_correction_path}/summary.bin";
	$self->{alignment_correction_done_file_name} = "$self->{alignment_correction_path}/alignment_resolution.done";
	$self->{jc_genome_diff_file_name} = "$self->{alignment_correction_path}/jc_evidence.gd";

	##### index BAM #####
	$self->{bam_path} = "06_bam";
	$self->{bam_path} = "$self->{base_output_path}/$self->{bam_path}" if ($self->{base_output_path});
	$self->{reference_bam_unsorted_file_name} = "$self->{bam_path}/reference.unsorted.bam";
	$self->{junction_bam_unsorted_file_name} = "$self->{bam_path}/junction.unsorted.bam";
	$self->{junction_bam_prefix} = "$self->{bam_path}/junction";
	$self->{junction_bam_file_name} = "$self->{bam_path}/junction.bam";
	$self->{bam_done_file_name} = "$self->{bam_path}/bam.done";
	
	##### error rates and coverage distribution #####
	$self->{error_calibration_path} = "07_error_calibration";
	$self->{error_calibration_path} = "$self->{base_output_path}/$self->{error_calibration_path}" if ($self->{base_output_path});
	$self->{error_counts_file_name} = "$self->{error_calibration_path}/#.error_counts.tab";
	#FOR TESTING: $self->{complex_error_counts_file_name} = "$self->{error_calibration_path}/#.complex_error_counts.tab";
	$self->{error_rates_file_name} = "$self->{error_calibration_path}/#.error_rates.tab";
	$self->{error_counts_done_file_name} = "$self->{error_calibration_path}/error_counts.done";
	$self->{error_rates_done_file_name} = "$self->{error_calibration_path}/error_rates.done";
	$self->{coverage_file_name} = "$self->{error_calibration_path}/@.coverage.tab";
	$self->{unique_only_coverage_distribution_file_name} = "$self->{error_calibration_path}/@.unique_only_coverage_distribution.tab";
	$self->{error_rates_summary_file_name} = "$self->{error_calibration_path}/summary.bin";
	$self->{error_rates_base_qual_error_prob_file_name} = "$self->{error_calibration_path}/base_qual_error_prob.#.tab";
	$self->{plot_error_rates_r_script_file_name} = "$self->{lib_path}/plot_error_rate.r";
	$self->{plot_error_rates_fit_r_script_file_name} = "$self->{error_calibration_path}/fit.#.r_script";
	$self->{plot_error_rates_r_script_log_file_name} = "$self->{error_calibration_path}/#.plot_error_rate.log";

	##### mutation identification #####
	$self->{mutation_identification_path} = "08_mutation_identification";
	$self->{mutation_identification_path} = "$self->{base_output_path}/$self->{mutation_identification_path}" if ($self->{base_output_path});
	$self->{predicted_mutation_file_name} = "$self->{mutation_identification_path}/@.predicted_mutations.bin";
	$self->{ra_mc_genome_diff_file_name} = "$self->{mutation_identification_path}/ra_mc_evidence.gd";
	$self->{complete_mutations_text_file_name} = "$self->{mutation_identification_path}/@.mutations.tab";
	$self->{complete_coverage_text_file_name} = "$self->{mutation_identification_path}/@.coverage.tab";
	$self->{mutation_identification_done_file_name} = "$self->{mutation_identification_path}/mutation_identification.done";
	$self->{cnv_coverage_tab_file_name} = "$self->{mutation_identification_path}/@.cnv_coverage.tab";
	$self->{genome_error_counts_file_name} = "$self->{mutation_identification_path}/error_counts.tab";
	$self->{polymorphism_statistics_input_file_name} = "$self->{mutation_identification_path}/polymorphism_statistics_input.tab";
	$self->{polymorphism_statistics_output_file_name} = "$self->{mutation_identification_path}/polymorphism_statistics_output.tab";
	$self->{polymorphism_statistics_r_script_file_name} = "$self->{lib_path}/polymorphism_statistics.r";
	$self->{polymorphism_statistics_r_script_log_file_name} = "$self->{mutation_identification_path}/polymorphism_statistics_output.log";
	$self->{polymorphism_statistics_ra_mc_genome_diff_file_name} = "$self->{mutation_identification_path}/ra_mc_evidence_polymorphism_statistics.gd";
	$self->{polymorphism_statistics_done_file_name} = "$self->{mutation_identification_path}/polymorphism_statistics.done";

	##### data #####
	## things in this location are needed for running post-processing steps
	$self->{data_path} = "data";
	$self->{data_path} = "$self->{base_output_path}/$self->{data_path}" if ($self->{base_output_path});
	$self->{reference_bam_prefix} = "$self->{data_path}/reference";
	$self->{reference_bam_file_name} = "$self->{data_path}/reference.bam";
	$self->{reference_fasta_file_name} = "$self->{data_path}/reference.fasta";
	$self->{reference_faidx_file_name} = "$self->{data_path}/reference.fasta.fai";
	$self->{reference_features_file_name} = "$self->{data_path}/reference.features.tab";
	$self->{unmatched_read_file_name} = "$self->{data_path}/#.unmatched.fastq";
	
	##### output #####
	## things in this location are part of the user-readable output
	$self->{output_path} = "output";
	$self->{output_path} = "$self->{base_output_path}/$self->{output_path}" if ($self->{base_output_path});
	$self->{output_done_file_name} = "$self->{output_path}/output.done";
	$self->{log_file_name} = "$self->{output_path}/log.txt";	
	$self->{index_html_file_name} = "$self->{output_path}/index.html";
	$self->{summary_html_file_name} = "$self->{output_path}/summary.html";
	$self->{marginal_html_file_name} = "$self->{output_path}/marginal.html";
	$self->{final_genome_diff_file_name} = "$self->{output_path}/output.gd";	
	$self->{local_evidence_path} = "evidence";
	$self->{evidence_path} = "$self->{output_path}/$self->{local_evidence_path}";
	$self->{evidence_genome_diff_file_name} = "$self->{evidence_path}/evidence.gd";
	$self->{local_coverage_plot_path} = "evidence";
	$self->{coverage_plot_path} = "$self->{output_path}/$self->{local_coverage_plot_path}";
	$self->{deletions_text_file_name} = "$self->{coverage_plot_path}/deletions.tab";
	$self->{coverage_plot_file_name} = "$self->{coverage_plot_path}/@.overview.png";
	$self->{output_calibration_path} = "$self->{output_path}/calibration";
	$self->{unique_only_coverage_plot_file_name} = "$self->{output_calibration_path}/@.unique_coverage.pdf";
	$self->{error_rates_plot_file_name} = "$self->{output_calibration_path}/#.error_rates.pdf";	
	
	## text output files, to be replaced...
	$self->{settings_text_file_name} = "$self->{output_path}/settings.tab";
	$self->{summary_text_file_name} = "$self->{output_path}/summary.tab";
	$self->{tiled_coverage_text_file_name} = "$self->{output_path}/@.tiled_coverage.tab";

	$self->{breseq_small_graphic_from_file_name} = "$self->{lib_path}/breseq_small.png";
	$self->{breseq_small_graphic_to_file_name} = "$self->{output_path}/$self->{local_evidence_path}/breseq_small.png";

	$self->{long_pairs_file_name} = "$self->{output_path}/long_pairs.tab";

	###### tmp ######
	$self->{tmp_path} = "tmp";
	$self->{tmp_path} = "$self->{base_output_path}/$self->{tmp_path}" if ($self->{base_output_path});
	
	#read sequence filenames are given as straight arguments
	@{$self->{read_fastq_list}} = @ARGV;

	## Read sequence file provided?
	if (scalar @{$self->{read_fastq_list}} == 0) 
	{
		print STDERR "No read sequence files provided.";
		pod2usage(1);
	}

	##
	# Order the input fastq files, remove their '.fastq' endings,
	# and give each read sequence a unique id (used in alignment database and file names)
	##
	my $fastq_file_index = -1;
	my @new_fastq_list;
	foreach my $raw_read_fastq_file (@{$self->{read_fastq_list}}) 
	{
		my $read_structure = {};
		
		#paired files have form READFILE1,READFILE2::MIN,MAX
		my @fastq_files = ($raw_read_fastq_file);
		$read_structure->{paired} = 0;
		if ($raw_read_fastq_file =~ m/^(.+),(.+)::(.+)-(.+)$/)
		{
			@fastq_files = ($1,$2);
			$read_structure->{min_pair_dist} = $3;
			$read_structure->{max_pair_dist} = $4;
			$read_structure->{paired} = 1;
		}
		@{$read_structure->{read_fastq_list}} = @fastq_files;
		push @{$self->{read_structures}}, $read_structure;
		
		foreach my $read_fastq_file (@fastq_files) 
		{
			$fastq_file_index++;
			#name without path or fastq ending
			my $read_file = $read_fastq_file;
			$read_file =~ s/\.fastq$//; #no trailing .fastq
			$read_file =~ s/.+\///; #no beginning path
			push @{$self->{read_file_base_names}}, $read_file;				
			push @{$read_structure->{base_names}}, $read_file;	
			push @{$self->{read_file_index_to_struct_index}}, $#{$self->{read_structures}};

			$self->{read_file_to_fastq_file_index}->{$read_file} = $fastq_file_index; 
			$self->{read_file_to_fastq_file}->{$read_file} = $read_fastq_file; 

			#index for keeping track of what file reads came from in alignment database
			#max is 256 b/c stored as unsigned byte in alignment database
			$self->throw("Maximum of 256 input files allowed.") if ($fastq_file_index > 255);
		}
		
		$read_structure->{base_name} = join ("-pair-", @{$read_structure->{base_names}});
	}
	@{$self->{read_fastq_list}} = @new_fastq_list;

	## Reference sequence provided?
	if (scalar @{$self->{reference_genbank_file_names}} == 0) 
	{
		print STDERR "No reference sequences provided (-r).";
		pod2usage(1);
	}
	
	$self->compare_to_saved_settings;

	$self->installed;

	return $self;
}

### Checks to see if there were any changes since the last time we were run.
sub compare_to_saved_settings
{
	my ($self) = @_;
	
	##load saved settings file

	##cycle through all keys of new settings and compare
	
	## warn if 
}

### Utility function to substitute specific details into a generic file name
sub file_name
{
	my ($self, $file_name_key, $sub_hash)= @_;
	my $file_name = $self->{$file_name_key};
	$file_name or $self->throw("Settings file \"$file_name_key\" not found.");

	return $self->substitute_file_name($file_name, $sub_hash);
}

sub substitute_file_name
{
	my ($self, $file_name, $sub_hash)= @_;
	
	foreach my $key (keys %$sub_hash)
	{		
		(!ref($file_name)) or $self->throw("Substituting in file name arrays not supported.");
		defined $sub_hash->{$key} or $self->throw("Value to substitute for $key not defined");
		$file_name =~ s/$key/$sub_hash->{$key}/g or $self->throw("No '$key' found to substitute for in file name name '$file_name'.");
	}
	
	#return as a list if this is a reference to a list, note that replacing things in names does not work
	return @$file_name if (ref($file_name) eq 'ARRAY');
		
	#otherwise return single value
	return $file_name;
}

sub html_path
{
	my ($self, $file_name_key, $sub_hash)= @_;
	my $file_name = $self->file_name($file_name_key,$sub_hash);	
	
	## strip off the leading output path so we have a relative link
	if (substr ($file_name, 0, length($self->{output_path})+1) eq "$self->{output_path}/")
	{
		substr ($file_name, 0, length($self->{output_path})+1) = '';
	}
# regular expression chokes when there are special characters like '+' in the string 
# that are not escaped, it's treated like a double quoted string here.
#	$file_name =~ s/^$self->{output_path}\///; 
	$file_name or $self->throw("Settings file \"$file_name_key\" not found.");
}


sub create_path
{
	my ($self, $path_key) = @_;
	my $path = $self->file_name($path_key);
	(-e $path) or Breseq::File::Path::make_path($path) or $self->throw("Could not create path \'$path\'.");
	return $path;
}

sub remove_path
{
	my ($self, $path_key) = @_;
	my $path = $self->file_name($path_key);
	(-e $path) and Breseq::File::Path::remove_tree($path) or $self->throw("Could not remove path \'$path\'.");
	return $path;
}

sub read_file_to_fastq_file_index
{
	my ($self, $read_file) = @_;
	$self->throw() if (!defined $self->{read_file_to_fastq_file_index}->{$read_file});
	return $self->{read_file_to_fastq_file_index}->{$read_file};
}

#transparent to whether read trimming is on
sub read_file_to_fastq_file_name
{
	my ($self, $read_file) = @_;
	
	if (defined $self->{read_file_to_converted_fastq_file}->{$read_file}) {
		return $self->{read_file_to_converted_fastq_file}->{$read_file};
	}
	$self->throw if (!defined $self->{read_file_to_fastq_file}->{$read_file});
	return $self->{read_file_to_fastq_file}->{$read_file};
}

# same as above but forces original fastq file name if trimming is on
sub read_file_to_original_fastq_file_name
{
	my ($self, $read_file) = @_;
	$self->throw if (!defined $self->{read_file_to_fastq_file}->{$read_file});
	return $self->{read_file_to_fastq_file}->{$read_file};
}

sub read_files
{
	my ($self) = @_;
	$self->throw if (!defined $self->{read_file_base_names});
	return @{$self->{read_file_base_names}};
}

sub read_fastq_file_names
{
	my ($self) = @_;
	$self->throw if (!defined $self->{read_fastq_list});
	return @{$self->{read_fastq_list}};
}

#information about what is paired, etc...
sub read_structures
{
	my ($self) = @_;
	$self->throw if (!defined $self->{read_structures});
	return @{$self->{read_structures}};
}

sub log
{
	my ($self, $message) = @_;
	
	$self->create_path('output_path');
	$self->create_path('output_calibration_path');
	
	open LOG, ">>", $self->file_name('log_file_name');
	print LOG "$message\n";
	close LOG;
}

sub ctool
{
	my ($self, $tool_name, $allow_fail) = @_;
		
	if (!$self->{installed}->{$tool_name})
	{
		if ($allow_fail)
		{
			$self->warn("Executable \"$tool_name\" not found in breseq bin path\"$self->{bin_path}\".");
			return undef; # couldn't find it, but it's not an error.
		}
		else
		{
			$self->throw("Executable \"$tool_name\" not found in breseq bin path\"$self->{bin_path}\".");
		}
	}
	return "$self->{installed}->{$tool_name}";
}

sub installed
{
	my ($self) = @_;
	
	## breseq C++ executables
	$self->{installed}->{cbreseq} = (-x "$self->{bin_path}/cbreseq") ? "$self->{bin_path}/cbreseq" : 0;

	## absolutely required ssaha2 or smalt
	$self->{installed}->{SSAHA2} = (`which ssaha2`) ? "ssaha2" : 0;
	
	## check for default names
	$self->{installed}->{smalt} = (`which smalt`) ? "smalt" : 0;
	if (!$self->{installed}->{smalt}) { $self->{installed}->{smalt} = (`which smalt_i386`) ? "smalt_i386" : 0; }
	if (!$self->{installed}->{smalt}) { $self->{installed}->{smalt} = (`which smalt_ia64`) ? "smalt_ia64" : 0; }
	if (!$self->{installed}->{smalt}) { $self->{installed}->{smalt} = (`which smalt_x86_64`) ? "smalt_x86_64" : 0; }
	if (!$self->{installed}->{smalt}) { $self->{installed}->{smalt} = (`which smalt_MacOSX_i386`) ? "smalt_MacOSX_i386" : 0; }

	$self->{installed}->{R} = (`which R`) ? "R" : 0;	
	if ($self->{installed}->{R}) 
	{
		my $R_version = `R --version`;
		if ($R_version =~ m/R\s+(version\s+|)(\d+)\.(\d+)\.(\d+)/)
		{
			$self->{installed_version}->{R} = $2 * 1000000 +  $3 * 1000 + $4;
		}
		else
		{
			$self->{installed_version}->{R} = 0;
		}
	}
	
	$self->{installed}->{Statistics_Distributions} = (eval 'require Statistics::Distributions');	
	$self->{installed}->{samtools} = (-x "$self->{bin_path}/samtools") ? "$self->{bin_path}/samtools" : 0;
	$self->{installed}->{bioperl} = (eval 'require Bio::Root::Root');	

	## installed locally
	$self->{installed}->{Bio_DB_Sam} = (eval 'require Bio::DB::Sam');	
}

sub check_installed
{
	my ($self) = @_;

	my $good_to_go = 1;
	
	if (!$self->{smalt} && !$self->{installed}->{SSAHA2})
	{
		$good_to_go = 0;
		print STDERR "---> ERROR Required executable \"ssaha2\" not found.\n";
		print STDERR "---> See http://www.sanger.ac.uk/resources/software/ssaha2\n";
	}
	
	if ($self->{smalt} && !$self->{installed}->{smalt})
	{
		$good_to_go = 0;
		print STDERR "---> ERROR Required executable \"smalt\" not found.\n";
		print STDERR "---> See http://www.sanger.ac.uk/resources/software/smalt/\n";
	}

	## R version 2.1 required
	if (!$self->{installed}->{R})
	{
		$good_to_go = 0;
		print STDERR "---> ERROR Required executable \"R\" not found.\n";
		print STDERR "---> See http://www.r-project.org\n";
	}
	elsif ( (!defined $self->{installed_version}->{R}) || ($self->{installed_version}->{R} < 2001000) )
	{
		my $R_version = 'unknown';
		if ($self->{installed_version}->{R})
		{
			$R_version = int($self->{installed_version}->{R}/1000000) 
				. "." . int($self->{installed_version}->{R}%1000000/1000)
				. "." . int($self->{installed_version}->{R}%1000);
		}
		
		$good_to_go = 0;
		print STDERR "---> ERROR Required executable \"R version 2.1.0 or later\" not found.\n";
		print STDERR "---> Your version is $R_version\n";
		print STDERR "---> See http://www.r-project.org\n";
	}
	
	if (!$self->{installed}->{samtools})
	{
		$good_to_go = 0;
		print STDERR "---> ERROR Required executable \"samtools\" not found.\n";
		print STDERR "---> This should have been installed by the breseq installer.\n";
	}
	
	if (!$self->{installed}->{bioperl})
	{
		$good_to_go = 0;
		print STDERR "---> ERROR Required \"Bioperl\" modules not found.\n";
		print STDERR "---> See http://www.bioperl.org\n";
	}

	if (!$self->{installed}->{Bio_DB_Sam})
	{
		$good_to_go = 0;
		print STDERR "---> ERROR Required Perl module \"Bio::DB::Sam\" not found.\n";
		print STDERR "---> This module should have been installed by the breseq installer.\n";
	}
	
	if ( ($self->{perl_identify_mutations} && $self->{polymorphism_prediction} && !$self->{installed}->{Statistics_Distributions}) )
	{
		$good_to_go = 0;
		print STDERR "---> ERROR Perl module Statistics::Distributions not found.\n";
		print STDERR "---> Required for Perl polymorphism prediction.\n";
	}
	
	die "\n" if (!$good_to_go);
	
#	return $good_to_go;
}


sub do_step
{
	my ($self, $done_key, $message) = @_;
	
	my $done_file_name = $self->file_name($done_key);
	$self->{done_key_messages}->{$done_key} = $message;
	if (!-e $done_file_name)
	{
		print STDERR "+++   NOW PROCESSING $message\n";
		$self->record_start_time($message);
		return 1;
	}
	
	print STDERR "--- ALREADY COMPLETE $message\n";
	
	my $time;
	$time = Storable::retrieve($done_file_name) if (-s $done_file_name > 0);
	if (!$time) 
	{	
		$time = {};
		$self->warn("Can't retrieve time data from file $done_file_name");
	}
	push @{$self->{execution_times}}, $time;
	
	return 0;
}

sub done_step
{
	my ($self, $done_key) = @_;
	
	my $done_file_name = $self->file_name($done_key);
	my $message = $self->{done_key_messages}->{$done_key};
	$self->record_end_time($message);
	
	## create the done file with timing information
	Storable::store($self->{execution_times}->[-1], $done_file_name)
		or $self->throw("Can't store time data in file $done_file_name");
}

sub record_start_time
{
	my ($self, $message) = @_;
	
	my $this_time = time;
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime($this_time);	
	my $formatted_time = sprintf "%04d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $mday, $hour, $min, $sec;
	my $new_time = { 
		_time_start => $this_time,
		_formatted_time_start => $formatted_time, 
		_time_end => 0,
		_formatted_time_end => '', 
		_time_elapsed => 0,
		_formatted_time_elapsed => '',
		_message => $message,
	 };
	
	push @{$self->{execution_times}}, $new_time;
	return $formatted_time;
}

sub record_end_time
{
	my ($self, $message) = @_;
	
	my $i = 0;
	while ($i < scalar @{$self->{execution_times}})
	{
		last if ($self->{execution_times}->[$i]->{_message} eq $message);
		$i++;
	}
	
	if ($i >= scalar @{$self->{execution_times}})
	{
		$self->warn("Did not find matching start time for:\n$message");
	}
	
	my $ex_time = $self->{execution_times}->[$i];
	$ex_time->{_message} = $message;

	my $this_time = time;
	$ex_time->{_time_end} = $this_time;
	$ex_time->{_formatted_time_end} = time2string($this_time);

	##if we had a previous time, calculate elapsed
	if ($i < scalar @{$self->{execution_times}})
	{
		my $time_interval = $ex_time->{_time_end} - $ex_time->{_time_start};
		$ex_time->{_time_elapsed} = $time_interval;
		$ex_time->{_formatted_time_elapsed} = time2string($time_interval, 1);
	}
}

sub time2string
{
	my ($this_time, $relative) = @_;
	my $s = '';

	## time is absolute
	if (!$relative)
	{
		my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime($this_time);	
		$s = sprintf "%04d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $mday, $hour, $min, $sec;
	}
	
	##time is an elapsed interval
	else 
	{
		my $rem_time = $this_time;
		
		
		my $sec = $rem_time % 60;
		$rem_time = int($rem_time / 60);
		$s = "$sec sec";
		
		if ($rem_time)
		{
			my $min = $rem_time % 60;
			$rem_time = int($rem_time / 60);
			$s = "$min min " . $s;
		}
		
		if ($rem_time)
		{
			my $hours = $rem_time % 24;
			$rem_time = int($rem_time / 24);
			$s = "$hours hr " . $s;
		}
		
		if ($rem_time)
		{
			my $days = $rem_time;
			
			if ($days == 1)
			{
				$s = "$days day " . $s;
			}
			else
			{
				$s = "$days days " . $s;
			}
		}
	}
	
	return $s;
}


return 1;