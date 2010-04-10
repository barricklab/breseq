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

use Bio::Location::Simple;
use Bio::DB::Flat;
use Bio::SeqIO;

my ( $input_1, $input_2, $output ) = @ARGV;

#print "$input_1 $input_2 $output\n";

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



#Open Input/output

open INPUT_1, "<$input_1";
open INPUT_2, "<$input_2";
open OUTPUT, ">$output";


#my $in_1 = Bio::SeqIO->new( -file => "<$input_1", -format => "fastq");
#my $in_2 = Bio::SeqIO->new( -file => "<$input_2", -format => "fastq");
#my $out  = Bio::SeqIO->new( -file => ">$output", -format => "fastq");

#my $seq_1 = $in_1->next_seq;
#my $seq_2 = $in_2->next_seq;

my $line_1 = <INPUT_1>;
my $line_2 = <INPUT_2>;
while ($line_1)
{
	foreach my $count (1..4)
	{
		print OUTPUT $line_1;
		$line_1 = <INPUT_1>;
	}

	foreach my $count (1..4)
	{
		print OUTPUT $line_2;
		$line_2 = <INPUT_2>;
	}

	#$seq_1 = $in_1->next_seq;
	#$seq_2 = $in_2->next_seq;
	
	#$out->write_seq($seq_1);
	#$out->write_seq($seq_2);
}
