#!/usr/bin/env perl -w

###
# Pod Documentation
###

=head1 NAME

subsequence.pl

=head1 SYNOPSIS

Usage: merge_paired_end_reads.pl seq_file_1 seq_file_2 output_file

Merge two files containing paired-end reads.

=head1 DESCRIPTION

=over

=item B<-i> <file path> 

Input sequence file.

=item B<-1|--start> <int> 

Start coordinate.

=item B<-2|--end> <int> 

End coordinate.

=back

=head1 AUTHOR

Jeffrey Barrick

=head1 COPYRIGHT

Copyright 2007.  All rights reserved.

=cut

###
# End Pod Documentation
###


use strict;
use Data::Dumper;

use FindBin;
use lib $FindBin::Bin;

my $command = shift @ARGV;
$command = "\U$command";

my $output_gd;
if ($command eq 'NUMBER_PAIRED')
{
	use Breseq::Fastq;
	use Getopt::Long;
	use Pod::Usage;
	my ($help, $man);
	GetOptions(
		'help|?' => \$help, 'man' => \$man,
	) or pod2usage(2);
	pod2usage(1) if $help;
	pod2usage(-exitstatus => 0, -verbose => 2) if $man;	
		
	my ( $input_file_2_1, $input_file_2_2 ) = @ARGV;
	
	my @in_list = ($input_file_2_1, $input_file_2_2); 
	
	my $input_file_2_1_base = $input_file_2_1;
	$input_file_2_1_base =~ s/\.(\S+?)$//;
	my $input_file_2_1_suffix = $1 ? $1 : "";
	
	my $input_file_2_2_base = $input_file_2_2;
	$input_file_2_2_base =~ s/\.(\S+?)$//;
	my $input_file_2_2_suffix = $1 ? $1 : "";	
	
	my @out_list = ("$input_file_2_1_base.rn.$input_file_2_1_suffix", "$input_file_2_2_base.rn.$input_file_2_2_suffix");
	
	for (my $i=1; $i<=scalar @in_list; $i++)
	{
		my $in_fastq = FastqLite->new(-file => $in_list[$i-1]);
		my $out_fastq = FastqLite->new(-file => ">" . $out_list[$i-1]);
		my $n = 1;
		
		while ($_ = $in_fastq->next_seq)
		{
			$_->{id} .= "_$i\_$n";
			$out_fastq->write_seq($_);
			$n++;
		}
	}
}

## Converts a FASTQ to a series of FASTA files
## --> was originally used for MUMmer input.
if ($command eq 'FASTA')
{
	use Breseq::Fastq;
	use Getopt::Long;
	use Pod::Usage;
	
	my ( $input, $output ) = @ARGV;
	my $options;
	my ($help, $man);
	GetOptions(
		'help|?' => \$help, 'man' => \$man,
		'strip|s' => \$options->{strip_names},
		'paired|p' => \$options->{paired},
		'number-per-file|n=s' => \$options->{number_per_file},
	) or pod2usage(2);
	pod2usage(1) if $help;
	pod2usage(-exitstatus => 0, -verbose => 2) if $man;
	
	my $reads = Breseq::Fastq::split_fastq_to_fasta($input, $output, $options);
	print "$reads reads converted.\n"
}


## Converts two paired FASTQ files to one
## with the FASTQ entries interspersed.

if ($command eq 'MERGE_PAIRED')
{
	
	#Get options
	use Getopt::Long;
	use Pod::Usage;
	my ($help, $man);
	my ($input_file, $start, $end);
	GetOptions(
		'help|?' => \$help, 'man' => \$man,
	) or pod2usage(2);
	pod2usage(1) if $help;
	pod2usage(-exitstatus => 0, -verbose => 2) if $man;

	my ( $input_1, $input_2, $output ) = @ARGV;

	#Open Input/output
	my $in_fastq_1 = FastqLite->new(-file => $input_1);
	my $in_fastq_2 = FastqLite->new(-file => $input_2);
	my $out_fastq = FastqLite->new(-file => ">$output");

	my $line_1 = <INPUT_1>;
	my $line_2 = <INPUT_2>;
	while (($seq1 = $in_fastq->next_seq) && ($seq2 = $in_fastq->next_seq))
	{
		$out_fastq->write_seq($seq1);
		$out_fastq->write_seq($seq2);
	}
}

## This seems to convert the quality scores...
## Not really sure...

if ($command eq "QUALS_FASTQ")
{
	#Get options
	use Getopt::Long;
	use Pod::Usage;
	my ($help, $man);
	my ($input_file, $start, $end);
	my $strip_names = 0;
	GetOptions(
		'help|?' => \$help, 'man' => \$man,
		'strip-names|s' => \$strip_names,
	) or pod2usage(2);
	pod2usage(1) if $help;
	pod2usage(-exitstatus => 0, -verbose => 2) if $man;
	pod2usage(-exitstatus => 0, -verbose => 2) if (scalar @ARGV != 2);
	my ( $input_1, $input_2 ) = @ARGV;

	# Solexa quality scores are -10 log10(P/(1-P))
	# $Q = ord($q) - 64;

	# PHRED (454) quality scores are -10 log10(P) 
	# $Q = ord($q) - 33;

	#Open Input/output

	open INPUT_1, "<$input_1";
	open INPUT_2, "<$input_2";

	#my $in_1 = Bio::SeqIO->new( -file => "<$input_1", -format => "fastq");
	#my $in_2 = Bio::SeqIO->new( -file => "<$input_2", -format => "fastq");
	#my $out  = Bio::SeqIO->new( -file => ">$output", -format => "fastq");

	#my $seq_1 = $in_1->next_seq;
	#my $seq_2 = $in_2->next_seq;

	my $line_1 = <INPUT_1>;
	my $line_2 = <INPUT_2>;

	my $num_seq = 0;
	while ($line_1)
	{
		#
		$num_seq++;
		chomp $line_1;
		my $read_name = $line_1;
		$read_name =~ s/^>/@/;
		if ($strip_names)
		{
			$read_name = "\@$num_seq\n";
		}

		##collect sequence
		my $sequence;

		while ($line_1 = <INPUT_1>)
		{
			last if ($line_1 =~ m/^>/);
			chomp $line_1;
			$sequence .= $line_1;
		}

		##collect qualities
		my @quality_list;

		while ($line_2 = <INPUT_2>)
		{
			last if ($line_2 =~ m/^>/);
			chomp $line_2;
			push @quality_list, split /\s+/, $line_2;
		}

		#output...
		print "$read_name\n$sequence\n+\n";
	#	print +(phred_quality_to_char(0));
		foreach my $quality (@quality_list)
		{
			print +(phred_quality_to_char($quality));
		}
		print "\n";


		##print in fastq format

		#$seq_1 = $in_1->next_seq;
		#$seq_2 = $in_2->next_seq;

		#$out->write_seq($seq_1);
		#$out->write_seq($seq_2);
	}

	#my $index_m1 = $index-1;
	#system "cat $input_1\_$index_m1 $input_1\_$index > temp; mv temp > $input_1\_$index_m1";
	#system "cat $input_2\_$index_m1 $input_2\_$index > temp; mv temp > $input_2\_$index_m1";



	sub phred_to_solexa_quality
	{
		my ($pq) = @_;
		my $e = exp(-$pq * log(10) / 10);

		my $sq = -10 * log ($e / (1-$e)) / log (10);
		return int floor $sq;
	}

	sub phred_quality_to_char
	{
		my ($q) = @_;
		return chr($q+33);
	}
}

#converts a fastq list file (where the quality scores are whitespace separated integers
#to a standard format where they are characters.

if ($command eq "QUALS_FASTQ")
{
	use Getopt::Long;
	use Pod::Usage;
	use Breseq::Fastq;
	my ($help, $man);
	my $quality_score_format = 'phred';
	GetOptions(
		'help|?' => \$help, 'man' => \$man,
		'quality-score-format|q=s' => \$quality_score_format,
	) or pod2usage(2);
	pod2usage(1) if $help;
	pod2usage(-exitstatus => 0, -verbose => 2) if $man;
	pod2usage(-exitstatus => 0, -verbose => 2) if (scalar @ARGV != 2);
	my ( $input_file, $output_file) = @ARGV;

	print "Using quality score format: $quality_score_format\n";
	my $in_fastq = Breseq::Fastq->new( -file => "<$input_file", -list_format => 1, -quality_format => $quality_score_format);
	my $out_fastq = Breseq::Fastq->new( -file => ">$output_file");
	
	while (my $read_seq = $in_fastq->next_seq) 
	{
		$out_fastq->write_seq(, $read_seq);
	}

}