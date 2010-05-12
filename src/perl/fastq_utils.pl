#!/usr/bin/env perl

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

if ($command eq 'SPLIT')
{
	use Breseq::Fastq;
	use Getopt::Long;
	use Pod::Usage;
	my ($help, $man);
	my $num_seq_per_file;
	GetOptions(
		'help|?' => \$help, 'man' => \$man,
		'num|m=s' => \$num_seq_per_file,
	) or pod2usage(2);
	pod2usage(1) if $help;
	pod2usage(-exitstatus => 0, -verbose => 2) if $man;	
	die "Must provide number of sequences per output file (-n)\n" if (!$num_seq_per_file);

	my ( $input_file ) = @ARGV;
		
	my $input_file_base = $input_file;
	$input_file_base =~ s/\.(\S+?)$//;
	my $input_file_suffix = $1 ? $1 : "";	
	
	my $i = 0;
	my $read_num;
	my $output_file;
	
	my $in_fastq = Breseq::Fastq->new(-file => $input_file);
	my $out_fastq;

	while ($_ = $in_fastq->next_seq)
	{	
		if (!$read_num || $read_num >= $num_seq_per_file)
		{
			$read_num = 0;
			$i++;
			$output_file = "$input_file_base.$i.$input_file_suffix";
			$out_fastq = Breseq::Fastq->new(-file => ">" . $output_file);	
		}
		$out_fastq->write_seq($_);
		$read_num++;	
	}
}


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
		my $in_fastq = Breseq::Fastq->new(-file => $in_list[$i-1]);
		my $out_fastq = Breseq::Fastq->new(-file => ">" . $out_list[$i-1]);
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
	my $in_fastq_1 = Breseq::Fastq->new(-file => $input_1);
	my $in_fastq_2 = Breseq::Fastq->new(-file => $input_2);
	my $out_fastq = Breseq::Fastq->new(-file => ">$output");

	my $line_1 = <INPUT_1>;
	my $line_2 = <INPUT_2>;
	
	my $seq1;
	my $seq2;
	while (($seq1 = $in_fastq_1->next_seq) && ($seq2 = $in_fastq_2->next_seq))
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
		$out_fastq->write_seq($read_seq);
	}

}

if ($command eq "SIMULATE")
{
	do_simulate();
}

sub do_simulate
{
	use Getopt::Long;
	use Pod::Usage;
	use Bio::SeqIO;
	use Breseq::Fastq;
	my ($help, $man);

	my $output = 'output.fastq';
	my $coverage = 50;
	my $read_length = 36;
	my $circular = 1;
	
	my $base_error_rate = 0.001;
	my $indel_error_rate = 0.0001;
	
	GetOptions(
		'help|?' => \$help, 'man' => \$man,
		'read-length|l=s' => \$read_length,
		'coverage|c=s' => \$coverage,
		'output|o=s' => \$output,
	) or pod2usage(2);
	pod2usage(1) if $help;
	pod2usage(-exitstatus => 0, -verbose => 2) if $man;
	my ( $input ) = @ARGV;
	pod2usage(1) if (!$input);

	my $ref_i = Bio::SeqIO->new( -file => "<$input");
	my $ref_seq = $ref_i->next_seq;
	my $ref_sequence_string = $ref_seq->seq;
	my $len = $ref_seq->length;

	## add to the end to so that we can just subseq enven though it is circular
#	$ref_sequence .= substr $seq, 0, $read_length;

	our @bases = ("A", "T", "C", "G");

	sub rand_base
	{
		return $bases[int rand 4];
	}
	
	my @ref_seq;
	$ref_seq[0] = $ref_sequence_string;
	$ref_seq[1] = Breseq::Fastq::revcom($ref_seq[0]);
	
	my $out_fastq = Breseq::Fastq->new(-file => ">$output");
	
	my $num_reads = $coverage * $len / $read_length;
	
	print "NUM READS: $num_reads\n";
	for (my $i=1; $i<=$num_reads; $i++)
	{	
		print STDERR "READ: $i\n" if ($i % 10000 == 0);
			
		my $strand = int rand(2);
		my $seq = '';
		my $pos = int rand($len);
		POS: while (length $seq < $read_length)
		{			
			$pos+=1;
			$pos = 0 if ($pos > $len);
			
			##deletion
			next if (rand(1) < $indel_error_rate);
			
			##mutation
			if (rand(1) < $base_error_rate)
			{
				$seq .= rand_base();
			}
			else
			{				
				$seq .= substr $ref_seq[$strand], $pos, 1;
			}
			
			##insertion
			while (rand(1) < $indel_error_rate)
			{
				last if (length $seq == $read_length);
				$seq .= rand_base();
			}
			
		}
		
		my $read_seq;
		$read_seq->{id} = "READ-$i";
		$read_seq->{seq} = $seq;				
		$read_seq->{qual_chars} = 'A' x int(length($read_seq->{seq})/2);
		$read_seq->{qual_chars} .= 'B' x (length($seq) - length($read_seq->{qual_chars}));
		$out_fastq->write_seq($read_seq);
	}
}