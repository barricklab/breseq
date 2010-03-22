#!/usr/bin/perl -w

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

use FastqLite;

my $command = shift @ARGV;

#print "$input_1 $input_2 $output\n";

#Get options
use Getopt::Long;
use Pod::Usage;
my ($help, $man);
my ($input_file_1, $input_file_2_1, $input_file_2_2, $output);
GetOptions(
	'help|?' => \$help, 'man' => \$man,
	'input|i=s' => \$input_file_1,
	'input_1|1=s' => \$input_file_2_1,
	'input_2|2=s' => \$input_file_2_2,
) or pod2usage(2);
pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;

print Dumper(@ARGV);

my $output_gd;
if ("\U$command" eq 'NUMBER_PAIRED')
{
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
