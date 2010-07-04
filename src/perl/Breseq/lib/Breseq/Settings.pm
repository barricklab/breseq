###
# Pod Documentation
###

=head1 NAME

Breseq::Shared.pm

=head1 SYNOPSIS

Various utility functions.

=head1 DESCRIPTION

=head1 AUTHOR

Jeffrey Barrick
<jbarrick@msu.edu>

=head1 COPYRIGHT

Copyright 2008.  All rights reserved.

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

	$self->initialize_1;
	
	#keep this on by default
	
	my ($help, $man);
	my ($resume_run, $continue_run);
	GetOptions(
		'help|?' => \$help, 'man' => \$man,
		'verbose|v' => \$self->{verbose},
	## Options for input and output files
		'name|n=s' => \$self->{run_name},	
		'output-path|o=s' => \$self->{base_output_path},	
		'reference-sequence|r=s' => \@{$self->{reference_genbank_file_names}},
		'junction-sequence|j=s' => \@{$self->{junction_only_reference_genbank_file_names}},		
		'quality-style|q=s' => \$self->{quality_type},
		'clean=s' => \$self->{clean},
	## Options for snp error analysis
		'require-complete-match' => \$self->{require_complete_match},
		'require-no-indel-match' => \$self->{require_no_indel_match},
		'require-unique-match' => \$self->{require_unique_match},
		'require-max-mismatches=s' => \$self->{require_max_mismatches},
		'do-not-trim-ambiguous-ends' => \$self->{do_not_trim_ambiguous_ends},
		'polymorphism-prediction' => \$self->{polymorphism_prediction},		
	## Options for turning various analysis chunks off or on
		'no-junction-prediction' => \$self->{no_junction_prediction},
		'no-mismatch-prediction' => \$self->{no_mutation_prediction},
		'no-deletion-prediction' => \$self->{no_deletion_prediction},
		'no-alignment-generation' => \$self->{no_alignment_generation},
		'no-filter-unwanted' => \$self->{no_filter_unwanted},
		'copy-number-variation' => \$self->{copy_number_variation},		
		'read-limit|l=s' => \$self->{read_limit},
		'candidate-junction-read-limit=s' => \$self->{candidate_junction_read_limit},
		'alignment-read-limit=s' => \$self->{alignment_read_limit},
		'polymorphism-log10-e-value-cutoff=s' => \$self->{polymorphism_log10_e_value_cutoff},
		'polymorphism-frequency-cutoff=s' =>  \$self->{polymorphism_frequency_cutoff},	
		'polymorphism-p-value-cutoff=s' => \$self->{polymorphism_bias_p_value_cutoff},
		'no-unmatched-reads' => \$self->{no_unmatched_reads},
		'maximum-candidate-junctions=s' => \$self->{maximum_candidate_junctions},
		'error-model=s' => \$self->{error_model_method},
		'trim-read-ends' => \$self->{trim_read_ends},
		'trim-read-base-quality=s' => \$self->{trim_read_base_quality},
		'base-quality-cutoff|b=s' => \$self->{base_quality_cutoff},
		'maximum-mismatches|m=s' => \$self->{maximum_read_mismatches},		
		'perl-error-count' => \$self->{perl_error_count},
	) or pod2usage(2);

	pod2usage(1) if $help;
	pod2usage(-exitstatus => 0, -verbose => 2) if $man;
	pod2usage(-exitstatus => 0, -verbose => 2) if (scalar @ARGV == 0);
	
	$self->initialize_2;
	
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
	my $time_stamp = sprintf "%4d-%02d-%02d %02d:%02d:%02d",$year+1900,$mon+1,$mday,$hour,$min,$sec;
	$self->log($time_stamp);	
	$self->log($self->{full_command_line});
	
	return $self;
}


sub new_annotate
{	
	my($caller,@args) = @_;
	my $class = ref($caller) || $caller;
	my $self = {}; # = new Bio::Root::Root($caller, @args);
	bless ($self, $class);

	$self->initialize_1;
	
	#keep this on by default
	my ($help, $man);
	my ($resume_run, $continue_run);
	GetOptions(
		'help|?' => \$help, 'man' => \$man,
		'verbose|v' => \$self->{verbose},
	## Options for input and output files
		'name|n=s' => \$self->{run_name},	
		'output-path|o=s' => \$self->{base_output_path},	
		'reference-sequence|r=s' => \@{$self->{reference_genbank_file_names}},
		'junction-sequence|j=s' => \@{$self->{junction_only_reference_genbank_file_names}},
		'quality-style|q=s' => \$self->{quality_type},
		'clean=s' => \$self->{clean},
	## Options for what results are printed
		'snp-quality-cutoff|c=s' => \$self->{snp_log10_prob_cutoff},
	## Options for snp error analysis
		'require-complete-match' => \$self->{require_complete_match},
		'require-no-indel-match' => \$self->{require_no_indel_match},
		'require-unique-match' => \$self->{require_unique_match},
		'require-max-mismatches=s' => \$self->{require_max_mismatches},
		'do-not-trim-ambiguous-ends' => \$self->{do_not_trim_ambiguous_ends},
		'polymorphism-prediction' => \$self->{polymorphism_prediction},		
	## Options for turning various analysis chunks off or on
		'no-junction-prediction' => \$self->{no_junction_prediction},
		'no-mismatch-prediction' => \$self->{no_mutation_prediction},
		'no-deletion-prediction' => \$self->{no_deletion_prediction},
		'no-alignment-generation' => \$self->{no_alignment_generation},
		'no-filter-unwanted' => \$self->{no_filter_unwanted},
		'copy-number-variation' => \$self->{copy_number_variation},		
		'read-limit|l=s' => \$self->{read_limit},
		'candidate-junction-read-limit=s' => \$self->{candidate_junction_read_limit},
		'alignment-read-limit=s' => \$self->{alignment_read_limit},
		'polymorphism-log10-e-value-cutoff=s' => \$self->{polymorphism_log10_e_value_cutoff},
		'polymorphism-frequency-cutoff=s' =>  \$self->{polymorphism_frequency_cutoff},	
		'polymorphism-p-value-cutoff=s' => \$self->{polymorphism_bias_p_value_cutoff},
		'no-unmatched-reads' => \$self->{no_unmatched_reads},
		'maximum-candidate-junctions=s' => \$self->{maximum_candidate_junctions},
		'trim-read-ends' => \$self->{trim_read_ends},
		'trim-read-base-quality=s' => \$self->{trim_read_base_quality},
		'base-quality-cutoff|b=s' => \$self->{base_quality_cutoff},
		'maximum-mismatches|m=s' => \$self->{maximum_read_mismatches},
		'perl-error-count' => \$self->{perl_error_count},
		
	## This is the only different option	
		'input-genome-diff|i=s' => \@{$self->{input_genome_diffs}},	
	) or pod2usage(2);

	pod2usage(1) if $help;
	pod2usage(-exitstatus => 0, -verbose => 2) if $man;
	pod2usage(-exitstatus => 0, -verbose => 2) if (scalar @ARGV == 0);
	
	$self->{max_junctions_to_print} = 100  if ($self->{marginal_mode});
	$self->{sort_junctions_by_score} = 1  if ($self->{marginal_mode});
		
	$self->initialize_2;
		
	return $self;
}

## called before getting options from command line
sub initialize_1
{
	my ($self) = @_;

	@{$self->{reference_genbank_file_names}} = ();  # files containing reference sequences	
	@{$self->{junction_only_reference_genbank_file_names}} = (); #files to look for junctions to but not align to
	
	## Set up default values for options
	$self->{full_command_line} = "$0 @ARGV"; 
	$self->{arguments} = "@ARGV";	
	$self->{quality_type} = '';					# quality score format
	$self->{predicted_quality_type} = '';
	$self->{min_quality} = 0;
	$self->{max_quality} = 0;
	$self->{min_match_length} = 24;
	$self->{run_name} = 'unnamed';
	$self->{clean} = 0;
	$self->{base_output_path} = '';
	$self->{error_model_method} = 'FIT';
	$self->{trim_reads} = 0;
	
	#used by CandidateJunctions.pm
	$self->{required_unique_length_per_side} = 10;				# Require at least one of the pair of matches supporting a junction to have this
																# much of its match that is unique in the reference sequence.
	$self->{maximum_inserted_junction_sequence_length} = 20;	# Ignore junctions with negative overlap (unique inserted sequence between reference 
																# matches) greater than this length. Prevents evidence file names that are too long.
	$self->{minimum_candidate_junctions} = 200;					# Minimum number of candidate junctions to keep
	$self->{maximum_candidate_junctions} = 5000;				# Maximum number of candidate junctions to keep
	$self->{maximum_candidate_junction_length_factor} = 0.1;	# Only keep CJ lengths adding up to this factor times the total reference size 
	$self->{candidate_junction_read_limit} = undef;		    	# FOR TESTING: only process this many reads when creating candidate junctions

	#used by AlignmentCorrection.pm
	$self->{add_split_junction_sides} = 1;

	$self->{no_junction_prediction} = undef;		# don't perform junction prediction steps
	$self->{no_mutation_prediction} = undef;		# don't perform read mismatch/indel prediction steps
	$self->{no_deletion_prediction} = undef;		# don't perform deletion prediction steps
	$self->{no_alignment_generation} = undef;		# don't generate alignments
	$self->{alignment_read_limit} = undef;			# only go through this many reads when creating alignments
	$self->{correction_read_limit} = undef;			# only go through this many reads when correcting alignments
	$self->{no_filter_unwanted} = undef;			# don't filter out unwanted reads with adaptor matches
	$self->{unwanted_prefix} = "UNWANTED:::";	# prefix on unwanted read names
	
	#used by MutationIdentification.pm
	$self->{mutation_log10_e_value_cutoff} = 2;			# log10 of evidence required for SNP calls 
	$self->{polymorphism_log10_e_value_cutoff} = 2;
	$self->{polymorphism_bias_p_value_cutoff} = 0.05;
	$self->{polymorphism_frequency_cutoff} = 0;   # cut off if < this or > 1-this
	$self->{polymorphism_coverage_both_bases} = 2;
	
	#used by Output.pm
	$self->{max_rejected_polymorphisms_to_show} = 100;
	$self->{max_rejected_junctions_to_show} = 25;
}

## called after getting options from command line
sub initialize_2
{
	my ($self) = @_;
	
	$self->{version} = $Breseq::VERSION;
	$self->{byline} = "<b><i>breseq</i></b> Version $self->{version} | Developed by Barrick JE and Knoester DB";
	$self->{bin_path} = $FindBin::Bin;
	$self->{lib_path} = "$self->{bin_path}/../lib/perl5/Breseq";

	#neaten up some settings for later string comparisons
	$self->{quality_type} = "\L$self->{quality_type}";
	$self->{error_model_method} = "\U$self->{error_model_method}";

	#on by default
	$self->{unmatched_reads} = ($self->{no_unmatched_reads}) ? 0 : 1;
	
	$self->{trim_reads} = $self->{trim_read_ends} || $self->{trim_read_base_quality};
	
	#######  SETUP FILE NAMES  #######
	### '#' replaced with read fastq name
	### '@' replaced by seq_id of reference sequence

	##### sequence conversion #####
	$self->{sequence_conversion_path} = "01_sequence_conversion";
	$self->{sequence_conversion_path} = "$self->{base_output_path}/$self->{sequence_conversion_path}" if ($self->{base_output_path});
	$self->{trimmed_fastq_file_name} = "$self->{sequence_conversion_path}/#.trimmed.fastq";
	$self->{reference_fasta_file_name} = "$self->{sequence_conversion_path}/reference.fasta";
	$self->{reference_faidx_file_name} = "$self->{sequence_conversion_path}/reference.fasta.fai";
    $self->{ref_seq_info_file_name} = "$self->{sequence_conversion_path}/ref_seq_info.bin";	
	$self->{unwanted_fasta_file_name} = "$self->{sequence_conversion_path}/unwanted.fasta";
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
	$self->{predicted_junction_file_name} = "$self->{alignment_correction_path}/jc_evidence.gd";
	$self->{unmatched_read_file_name} = "$self->{alignment_correction_path}/#.unmatched.fastq";
	$self->{alignment_correction_summary_file_name} = "$self->{alignment_correction_path}/summary.bin";
	$self->{alignment_correction_done_file_name} = "$self->{alignment_correction_path}/alignment_resolution.done";
	$self->{jc_genome_diff_file_name} = "$self->{alignment_correction_path}/jc_evidence.gd";

	##### index BAM #####
	$self->{bam_path} = "06_bam";
	$self->{bam_path} = "$self->{base_output_path}/$self->{bam_path}" if ($self->{base_output_path});
	$self->{reference_bam_unsorted_file_name} = "$self->{bam_path}/reference.unsorted.bam";
	$self->{reference_bam_prefix} = "$self->{bam_path}/reference";
	$self->{reference_bam_file_name} = "$self->{bam_path}/reference.bam";
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
	$self->{unique_coverage_distribution_r_script_file_name} = "$self->{error_calibration_path}/coverage.@.r_script";
	$self->{error_rates_r_script_file_name} = "$self->{error_calibration_path}/error_rates.#.r_script";
	
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

	##### output #####
	$self->{output_path} = "output";
	$self->{output_path} = "$self->{base_output_path}/$self->{output_path}" if ($self->{base_output_path});
	$self->{output_done_file_name} = "$self->{output_path}/output.done";
	$self->{log_file_name} = "$self->{output_path}/log.txt";	
	$self->{index_html_file_name} = "$self->{output_path}/index.html";
	$self->{summary_html_file_name} = "$self->{output_path}/summary.html";
	$self->{final_genome_diff_file_name} = "$self->{output_path}/$self->{run_name}.gd";
	$self->{evidence_genome_diff_file_name} = "$self->{output_path}/evidence.gd";
	$self->{local_evidence_path} = "evidence";
	$self->{evidence_path} = "$self->{output_path}/$self->{local_evidence_path}";
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
			$self->{read_file_to_trimmed_fastq_file}->{$read_file} = $self->file_name('trimmed_fastq_file_name', {'#' => $read_file});

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

	## if the resume option is chosen, all values will be read from saved settings file
	if ($self->{resume_run} || $self->{continue_run})
	{
		$self->load_from_file;
		print "Resuming run that was interrupted with same settings\n".
		return;
	}
	
	$self->compare_to_saved_settings;

	##
	## Clean old results if requested
	##

	if ($self->{clean} > 0)
	{
		rmtree([$self->{output_path}, $self->{hybrid_read_path}, $self->{alignment_correction_path}]);
	}
	if ($self->{clean} > 1)
	{
		rmtree([$self->{mummer_path}, $self->{fasta_conversion_path}]);
	}

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
	use File::Path qw(make_path);
	my ($self, $path_key) = @_;
	my $path = $self->file_name($path_key);
	(-e $path) or make_path($path) or $self->throw("Could not create path \'$path\'.");
	return $path;
}

sub remove_path
{
	use File::Path qw(remove_tree);
	my ($self, $path_key) = @_;
	my $path = $self->file_name($path_key);
	(-e $path) and remove_tree($path) or $self->throw("Could not remove path \'$path\'.");
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
	if ($self->{trim_reads})
	{
		$self->throw if (!defined $self->{read_file_to_trimmed_fastq_file}->{$read_file});
		return $self->{read_file_to_trimmed_fastq_file}->{$read_file};
	}
	else
	{
		$self->throw if (!defined $self->{read_file_to_fastq_file}->{$read_file});
		return $self->{read_file_to_fastq_file}->{$read_file};
	}
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
	my ($self, $tool_name) = @_;
		
	if (!$self->{installed}->{$tool_name})
	{
		$self->throw("Executable \"$tool_name\" not found in breseq bin path\"$self->{bin_path}\".");
	}
	return "$self->{bin_path}/$tool_name";
}

sub installed
{
	my ($self) = @_;
	
	## breseq C++ executables
	$self->{installed}->{error_count} = (-x "$self->{bin_path}/error_count") ? 1 : 0;

	## absolutely required
	$self->{installed}->{SSAHA2} = (`which ssaha2`) ? 1 : 0;
	$self->{installed}->{R} = (`which R`) ? 1 : 0;	
	$self->{installed}->{samtools} = (-x "$self->{bin_path}/samtools") ? 1 : 0;
	$self->{installed}->{bioperl} = (eval 'require Bio::Root::Root');	

	## installed locally
	$self->{installed}->{Statistics_Distributions} = (eval 'require Statistics::Distributions');	
	$self->{installed}->{Bio_DB_Sam} = (eval 'require Bio::DB::Sam');
		

	## optional
#	$self->{installed}->{Math_Pari} = (eval 'require Math::Pari');

	## optional, being phased out, seems to make no difference compared to Statistics::Distributions
	#$self->{installed}->{Math_GSL_CDF} = (eval 'require Math::GSL::CDF');	
}

sub check_installed
{
	my ($self) = @_;

	my $good_to_go = 1;
	
	if (!$self->{installed}->{SSAHA2})
	{
		$good_to_go = 0;
		print STDERR "---> ERROR Required executable \"ssaha2\" is not installed.\n";
		print STDERR "---> See http://www.sanger.ac.uk/resources/software/ssaha2/\n";
	}

	if (!$self->{installed}->{R})
	{
		$good_to_go = 0;
		print STDERR "---> ERROR Required executable \"R\" is not installed.\n";
		print STDERR "---> See http://www.r-project.org/\n";
	}
	
	if (!$self->{installed}->{samtools})
	{
		$good_to_go = 0;
		print STDERR "---> ERROR Required executable \"samtools\" is not installed.\n";
		print STDERR "---> This should have been installed by the breseq installer.\n";
	}
	
	if (!$self->{installed}->{bioperl})
	{
		$good_to_go = 0;
		print STDERR "---> ERROR Required \"Bioperl\" modules not installed.\n";
		print STDERR "---> They should have been installed by the breseq installer.\n";
	}

	if (!$self->{installed}->{Bio_DB_Sam})
	{
		$good_to_go = 0;
		print STDERR "---> ERROR Required Perl module \"Bio::DB::Sam\" not installed.\n";
		print STDERR "---> This should have been installed by the breseq installer.\n";
	}
	
	if ( ($self->{polymorphism_prediction} && !$self->{installed}->{Statistics_Distributions}) )
	{
		print STDERR "---> WARNING! Perl module Statistics::Distributions not installed. Fisher's Exact\n";
		print STDERR "---> Test will not be used to check strand bias in polymorphism predictions.\n";
		print STDERR "---> This should have been installed by the breseq installer.\n";		
	}
	
	$self->throw if (!$good_to_go);
	
#	return $good_to_go;
}

return 1;