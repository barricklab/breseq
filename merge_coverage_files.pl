#!/usr/bin/perl -w

###
# Pod Documentation
###

=head1 NAME

graph_coverage.pl

=head1 SYNOPSIS

Usage: graph_coverage.pl -i deletions.tab -c coverage.tab

Use R to create graphs of the read coverage at specified positions, or to 
show the boundaries of all predicted deletions in a genome.

=head1 DESCRIPTION

=over

=item B<-i> <file path> 

Input "deletions.tab" file produced by breseq.pl. These intervals will be graphed
in the output file "deletions.pdf".

=item B<-c> <file path> 

Path to "coverage.tab" file produced by breseq.pl.

=item B<--interval> <start>-<end>

Instead of using the -i option to input a file specifying positions to graph, you can
specify start-end coordinates (e.g. --interval=12345-12456). Multiple intervals can be 
specified, and the graphs will be made as separate pages in the output file "specific_coverage.pdf".

=back

=head1 AUTHOR

Jeffrey Barrick
<jbarrick@msu.edu>

=head1 COPYRIGHT

Copyright 2008.  All rights reserved.

=cut

##
## Note that really we should use the average unique coverage as a poisson mean
## and fit this form of generalized linear model to the observation in each file
## simultaneously to get a per-position offset to the counts. Then we can calculate
## a true expected 95% confidence interval. 
##

###
# End Pod Documentation
###
use strict;

use File::Path;
use FindBin;
use lib $FindBin::Bin;
use Data::Dumper;

#Get options
use Getopt::Long;
use Pod::Usage;

my ($help, $man);
my $verbose;
my $input;
my $output = "output.tab";
my $max_files = 99999;
my $window_size = 100;

GetOptions(
	'help|?' => \$help, 'man' => \$man,
	'input-path|p=s' => \$input,
	'output-path|o=s' => \$output,
	'window-size|w=s' => \$window_size,
	'verbose|v=s' => \$verbose,
);
pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;

open OUTPUT, ">$output" or die;
print OUTPUT "#window_size=$window_size\n";

opendir INPUT, $input;
my @complete_file_list = readdir INPUT;
@complete_file_list = grep m/\.coverage\.tab$/, @complete_file_list;

my @deletion_inputs;
my @coverage_inputs;
foreach my $f (@complete_file_list)
{
	$f =~ m/(.+?)\./;
	push @coverage_inputs, $f;
	push @deletion_inputs, "$1.deletions.tab"; 
}

#read complete length of all deletion files
my @deletion_lists;
for (my $i=0; $i<scalar @deletion_inputs; $i++)
{
	my $input = $deletion_inputs[$i];
	
	open IN, "$input";
	my @lines = <IN>;
	shift @lines;
	
	print "File $input\n";
	
	$deletion_lists[$i] = [];
	foreach my $l (@lines)
	{
		my @sl = split /\t/, $l; 
		push @{$deletion_lists[$i]}, { start=> $sl[1], end=>$sl[2] };
		
		print "  Deletion $sl[1]-$sl[2]\n";
	}
	last if ($i >= $max_files-1);
}

my @coverage_files;
my @unique_averages;
for (my $i=0; $i<scalar @coverage_inputs; $i++)
{
	my $input = $coverage_inputs[$i];
	
	my $unique_positions = 0;
	my $unique_total = 0;
	
	open COV, "$input";
	my $dontcare = <COV>;
	while (my $_ = <COV>)
	{
		##debug
		last if ($unique_positions > 200000);
		
		chomp $_;
		my @sl = split "\t", $_;
		
		if ($sl[2] + $sl[3] == 0)
		{
			$unique_total+= $sl[0] + $sl[1];
			$unique_positions++;
		}
	}
	my $unique_average = $unique_total/$unique_positions;
	push @unique_averages, $unique_average;
	close COV;
	
	print "Input $input, Unique Average $unique_average\n";
		
	last if ($i >= $max_files-1);
}


#open all coverage files
for (my $i=0; $i<scalar @coverage_inputs; $i++)
{
	my $input = $coverage_inputs[$i];
	
	open $coverage_files[$i], "$input";
	my $coverage_file = $coverage_files[$i];
	my $line = readline $coverage_file;	
	last if ($i >= $max_files-1);
}



#go through every position
#assumes every position exists in file
my $position = 0;
POS: while (1) {
	my @unique_coverages;
	
	my @top;
	my @bot;
	my @tot;
	my @count;
	
	for (my $i=0; $i<scalar @coverage_inputs; $i++)
	{
		$top[$i] = 0;
		$bot[$i] = 0;
		$tot[$i] = 0;
		$count[$i] = 0;
	}
	
	W: for (my $w=0; $w < $window_size; $w++)
	{	
		#each file
		$position++;
		my $i=0;
		FILE: while ($i<scalar @coverage_inputs)
		{		
			my $coverage_file = $coverage_files[$i];
			my $line = readline $coverage_file;
		
			last W if (!$line);
			#SHOULD check here to see if we are actually at the right position from the split line.
	
			#check the corresponding deletion we are on
			if (scalar @{$deletion_lists[$i]} > 0)
			{
				#remove first deletion if we are past it
				shift @{$deletion_lists[$i]} if ($position > $deletion_lists[$i]->[0]->{end});
			
				if (scalar @{$deletion_lists[$i]} > 0)
				{
					#ignore position if we are inside a deletion
					next FILE if ($position >= $deletion_lists[$i]->[0]->{start} && $position <= $deletion_lists[$i]->[0]->{end})
				}
			}
	
			#print "$line\n";
			chomp $line;
			my @sl = split "\t", $line;
		
			#note that we average out IS elements...
		
			if ($sl[2] + $sl[3] == 0)
			{
				$top[$i] += $sl[0] / $unique_averages[$i];
				$bot[$i] += $sl[1] / $unique_averages[$i];
				$tot[$i] += ($sl[0]+$sl[1]) / $unique_averages[$i];
				$count[$i] += 1;
			}
		
		} continue {	
			last if ($i >= $max_files-1);
			$i++;
		}
	}
	
	#calculate average and standard deviation
	my $mean = 'NA';
	my $stddev = 'NA';
	my $num = 0;
	
	foreach my $c (@count)
	{
		$num++ if ($c > 0);
	}
		
	if ($num > 0)
	{
		$mean = 0;
		for (my $i=0; $i< scalar @tot; $i++)
		{
			next if ($count[$i] == 0);
			$mean += $tot[$i]/$count[$i];
		}
		$mean /= $num;
		$mean = sprintf "%.5f", $mean;
	
		$stddev = 0;
		for (my $i=0; $i< scalar @tot; $i++)
		{
			next if ($count[$i] == 0);
			$stddev += ($tot[$i]/$count[$i]-$mean)**2;
		}
		$stddev /= $num;
		$stddev = $stddev**0.5;
		$stddev = sprintf "%.5f", $stddev;
		
	}

	print OUTPUT "$mean\t$stddev\t$num\n";
}