#!/usr/bin/env perl

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

###
# End Pod Documentation
###
use strict;

use FindBin;
use lib $FindBin::Bin;
use Data::Dumper;

#Get options
use Getopt::Long;
use Pod::Usage;

my ($help, $man);
my $verbose;
my ($input_file, $coverage_file);
my $output_deletions_file = 'deletions';
my $output_overview_file = 'coverage';
my $output_file = 'specific_coverage';
my @specific_intervals;
my $max_coverage;
my $downsample;
our $output_path;
my $tmp_path = ".";
my $drawing_options = "height=450, width=900";
#$drawing_options = "";

my $drawing_format = 'png';
my $seq_id;

GetOptions(
	'help|?' => \$help, 'man' => \$man,
	'input|i=s' => \$input_file,
	'coverage-file|c=s' => \$coverage_file,
	'output-deletions-file|o=s' => \$output_file,
	'output-overview-file|w=s' => \$output_overview_file,
	'interval=s' => \@specific_intervals,
	'downsample|d=s' => \$downsample,
	'max_coverage|m=s' => \$max_coverage, #applies only to overview graphs!
	'seq_id=s' => \$seq_id,
	'output-path|p=s' => \$output_path,
	'tmp-path|t=s' => \$tmp_path,
	'drawing-format=s' => \$drawing_format,
	'verbose|v=s' => \$verbose,
);

$output_deletions_file .= ".$drawing_format";
$output_overview_file .= ".$drawing_format";
$output_file .= ".$drawing_format";


pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;
pod2usage(1) if (!defined $coverage_file);
die "Coverage file does not exist: $coverage_file" if (!-e $coverage_file);
if (!defined $downsample)
{
	if (@specific_intervals)  
	{
		$downsample = 1;
	}
	elsif ($input_file)
	{
#		$downsample = 1000;
		$downsample = -1; #AUTO
	}
	else
	{
		$downsample = 1000;
#		$downsample = 'AUTO';
	}
}

mkdir $output_path if (defined $output_path);

print STDERR "Downsampling by $downsample\n" if ($verbose);

if (@specific_intervals)
{
	foreach my $interval (@specific_intervals)
	{
		my $interval_item;
		($interval_item->{start}, $interval_item->{end}) = split /-/, $interval;
		$interval = $interval_item;
	}
	
	#print Dumper(@specific_intervals);
	graph_coverage(\@specific_intervals, $coverage_file, $output_file, '', $downsample, $max_coverage, '');
}
elsif ($input_file)
{
	my @intervals = read_deletion_tab($input_file);
	#print Dumper(@intervals);
	@intervals = grep {$_->{seq_id} eq $seq_id} @intervals if ($seq_id);
	#print Dumper(@intervals);
	graph_coverage(\@intervals, $coverage_file, $output_deletions_file, $output_overview_file, $downsample, $max_coverage, $seq_id);
}
else #graph whole genome only
{
	graph_coverage([], $coverage_file, $output_deletions_file, $output_overview_file, $downsample, $max_coverage, $seq_id);
}

sub read_deletion_tab
{
	my ($file_name) = @_;
	
	open INFILE, "<$file_name";
	my $col_hash;
	
	#gather headers
	my $header_line = <INFILE>;
	chomp $header_line;
	my @header_list = split /\t/, $header_line;
	
	my @intervals;
	while (my $line = <INFILE>)
	{
		chomp $line;
		my @line_list = split /\t/, $line;
		
		my $new_interval;
		for (my $i=0; $i<scalar @line_list; $i++)
		{
			$new_interval->{$header_list[$i]} = $line_list[$i];
		}
		push @intervals, $new_interval;
	}
	
	close INFILE;
	
	return @intervals;
}


sub graph_coverage
{
	my ($intervals, $coverage_file, $output_deletions_file, $output_overview_file, $downsample, $max_coverage, $seq_id) = @_;
	
	$max_coverage = "max(cov\$redundant_top_cov[start:end]+cov\$redundant_bot_cov[start:end], cov\$unique_top_cov[start:end]+cov\$unique_bot_cov[start:end])" if (!defined $max_coverage);
	my $max_unique_coverage = "max(cov\$unique_top_cov[start:end]+cov\$unique_bot_cov[start:end])";

	open R_SCRIPT, ">$tmp_path/$$.coverage_graph.r_script";
	print R_SCRIPT <<END;
plot_coverage <- function(cov, del_start, del_end, my_title, genes)
{
	pos = c(1:dim(cov[1]));
	
	len<- del_end - del_start + 1;
	offset = len / 5;
	if (offset < 100)
	{
		offset <- 100; 
	}
	start <- del_start - offset;
	end <- del_end + offset;
	
	if (end > length(pos))
	{
		end = length(pos);
	}
	
	if (start < 1)
	{
		start = 1;
	}	
	
	maxy=$max_coverage;

#	options(scipen = 15)

	par(mar=c(6.1,4.1,1.2,0.2));
#	par(mar=c(8.1,4.1,4.1,2.1));
	plot(0:10, 0:10, type="n", lty="solid", ylim=c(0, maxy), xlim=c(pos[start], pos[end]), lwd=2, xaxs="i", yaxs="i", xlab="Coordinate in Reference Genome", ylab="Read Coverage Depth") #, main=my_title )
#	arrows(pos[del_start],maxy,pos[del_start], 0, length=0, col="black", lwd=3)
#	arrows(pos[del_end],maxy,pos[del_end], 0, length=0, col="black", lwd=3)
	rect(pos[start], 0, pos[del_start], maxy, col="grey80", lty=0)
	rect(pos[del_end], 0, pos[end], maxy, col="grey80", lty=0)

	lines(pos[start:end],cov\$redundant_top_cov[start:end]+cov\$redundant_bot_cov[start:end], type="s", col="red", lty="solid", lwd=2 )
	lines(pos[start:end],cov\$redundant_top_cov[start:end], type="s", col="yellow", lty="solid")
	lines(pos[start:end],cov\$redundant_bot_cov[start:end], type="s", col="orange", lty="solid")
	
	lines(pos[start:end],cov\$unique_top_cov[start:end]+cov\$unique_bot_cov[start:end], type="s", col="blue", lty="solid", lwd=2 )
	lines(pos[start:end],cov\$unique_top_cov[start:end], type="s", col="cyan", lty="solid")
	lines(pos[start:end],cov\$unique_bot_cov[start:end], type="s", col="purple", lty="solid")
		
	mtext(genes, side=1, col="blue", cex=0.7, outer=TRUE, line=-2)
}
plot_overview <- function(cov, start, end, my_title, genes)
{
	maxy=$max_unique_coverage;
	par(mar=c(6.1,4.1,1.2,0.2));
	#	par(mar=c(8.1,4.1,4.1,2.1));	
	plot(cov\$pos[start:end],cov\$redundant_top_cov[start:end]+cov\$redundant_bot_cov[start:end], type="s", col="red", lty="solid", lwd=1, xaxs="i", yaxs="i", ylim=c(0,maxy), xlab="Coordinate in Reference Genome", ylab="Read Coverage Depth") #, main=my_title )
	mtext(genes, side=1, col="blue", cex=0.7, outer=TRUE, line=-2)	
	lines(cov\$pos[start:end],cov\$unique_top_cov[start:end]+cov\$unique_bot_cov[start:end], type="s", col="blue", lty="solid", lwd=1 )
}

downsample <- function(m, factor)
{	
	new_entries = trunc(dim(m)[1] / factor);	
	out_m = data.frame(pos=c(1:new_entries));
	out_m\$pos = out_m\$pos * factor;
	
	#foreach column
	for (j in c('unique_top_cov', 'unique_bot_cov', 'redundant_top_cov', 'redundant_bot_cov')) 
	{
		#print(j)
		#print (m[[j]]);
		
		new_c = c(1:new_entries);
		for (i in 1:new_entries) 
		{ 
			start_coord <- 1 + (i-1)*factor
			end_coord <- i*factor
				
			if (end_coord > dim(m)[1])
			{
				end_coord = dim(m)[1]
			}
			new_c[i] = mean(m[[j]][start_coord:end_coord]);
		}
		out_m[[j]] = new_c;
	}
	return(out_m);
}

cov<-read.table("$coverage_file", header=T);

END
	
	if (!$output_path)
	{
		print R_SCRIPT "$drawing_format(\"$output_deletions_file\", $drawing_options)\n";
	}
		
	my $i=1;
	foreach my $interval (@$intervals)
	{
		# $interval->{genes} =~ s/\[/\\\[/;
		# $interval->{genes} =~ s/\]/\\\]/;
		$interval->{genes} = "" if (!defined $interval->{genes});
		$interval->{genes} =~ s/(.{80,}?)\s+/$1\\n/g;
		
		my $deletion_title = "$output_path/$i\.$drawing_format";
		$deletion_title = "$output_path/$seq_id.$i\.$drawing_format" if ($seq_id);
		
		print R_SCRIPT "$drawing_format(\"$deletion_title\", $drawing_options)\n" if ($output_path);
		my $title = "Predicted Deletion $interval->{start}-$interval->{end}";
		$title = "Predicted Deletion $seq_id:$interval->{start}-$interval->{end}" if ($seq_id);		
		print R_SCRIPT "plot_coverage(cov,$interval->{start},$interval->{end}, \"$title\", \"$interval->{genes}\")\n";
		print R_SCRIPT "dev.off()\n" if ($output_path);
		$i++;
	}	
		
	if (!$output_path)
	{
		print R_SCRIPT "dev.off()\n";
	}	
		
	if ($output_overview_file)
	{		
		print R_SCRIPT <<END;
		
downsample_by <- $downsample
if (downsample_by == -1)
{
	downsample_by = trunc(log10(dim(cov)[1]));
	downsample_by = downsample_by - 3;
	if (downsample_by < 0)
	{
		downsample_by = 0
	}
	downsample_by = 10**downsample_by
}
#print (downsample_by);
	
cov <- downsample(cov, downsample_by)
#print(cov)

END

		my $coverage_title = "$output_path/overview.$drawing_format";
		$coverage_title = "$output_path/$seq_id\.overview.$drawing_format" if ($seq_id);

		print R_SCRIPT "$drawing_format(\"$coverage_title\", $drawing_options)\n" if ($output_path);
		print R_SCRIPT "plot_overview(cov,1,length(cov\$pos), \"Genome Overview\", \"\")\n";
		print R_SCRIPT "dev.off()\n" if ($output_path);

#		print R_SCRIPT "$drawing_format(\"$output_path/cnv.$drawing_format\", $drawing_options)\n" if ($output_path);
#		print R_SCRIPT "plot_cnv(pos100,cov100,1,dim(pos100)[1], \"Copy Number Predictions\", \"\")\n";
#		print R_SCRIPT "dev.off()\n" if ($output_path);

	}

	close R_SCRIPT;
	
	my $command = "R --vanilla < $tmp_path/$$.coverage_graph.r_script &> $tmp_path/$$.coverage_graph.r_output";
	my $res = system $command;
	die "Error running command: $command" if ($res);
	
	unlink "$$.coverage_graph.r_script" if (!$verbose);
	unlink "$$.coverage_graph.r_output" if (!$verbose);
}
