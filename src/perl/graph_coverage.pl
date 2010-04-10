#!/usr/bin/env perl -w

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

use File::Path;
use FindBin;
use lib $FindBin::Bin;
use Data::Dumper;

#Get options
use Getopt::Long;
use Pod::Usage;

my ($help, $man);
my $verbose;
my ($input_file, $coverage_file);
my $output_deletions_file = 'deletions.pdf';
my $output_overview_file = 'coverage.pdf';
my $output_file = 'specific_coverage.pdf';
my @specific_intervals;
my $max_coverage;
my $downsample;
our $output_path;
my $pdf_options = "height=6, width=10";
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
	'verbose|v=s' => \$verbose,
);
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
		$downsample = 1000;
	}
	else
	{
		$downsample = 1000;
	}
}

mkdir $output_path if (defined $output_path);

print STDERR "Downsampling by $downsample\n";

if (@specific_intervals)
{
	foreach my $interval (@specific_intervals)
	{
		my $interval_item;
		($interval_item->{start}, $interval_item->{end}) = split /-/, $interval;
		$interval = $interval_item;
	}
	
	#print Dumper(@specific_intervals);
	new_graph_coverage(\@specific_intervals, $coverage_file, $output_file, $downsample, $max_coverage);
}
elsif ($input_file)
{
	my @intervals = read_deletion_tab($input_file);
	#print Dumper(@intervals);
	@intervals = grep {$_->{seq_id} eq $seq_id} @intervals if ($seq_id);
	#print Dumper(@intervals);
	graph_coverage(\@intervals, $coverage_file, $output_deletions_file, $output_overview_file, $downsample, $max_coverage);
}
else #graph whole genome only
{
	graph_coverage([], $coverage_file, $output_deletions_file, $output_overview_file, $downsample, $max_coverage);
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
	
	#my %header_hash;
	#my $header_list;
	#for (my $i=0; $i<scalar @header_list; $i++)
	#{
	#	$header_hash{$header_list[$i]} = $i;
	#}
	
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

sub new_graph_coverage
{
	my ($intervals, $coverage_file, $output_file, $downsample, $max_coverage) = @_;
	
	$max_coverage = "max(cov\$redundant_top_cov[start:end]+cov\$redundant_bot_cov[start:end], cov\$unique_top_cov[start:end]+cov\$unique_bot_cov[start:end])" if (!defined $max_coverage);

	open R_SCRIPT, ">$$.coverage_graph.r_script";
	print R_SCRIPT <<END;
plot_coverage <- function(pos, cov, del_start, del_end, downsample, my_title, genes)
{
	start <- trunc(del_start/downsample);
	end <- trunc(del_end/downsample);

	if (end > length(pos))
	{
		end = length(pos);
	}
	
	if (start < 1)
	{
		start = 1;
	}
	
	maxy=$max_coverage;
	options(scipen = 15)

	par(mar=c(8.1,4.1,4.1,2.1));
	plot(0:10, 0:10, type="n", lty="solid", ylim=c(0, maxy), xlim=c(pos[start], pos[end]), lwd=2, xlab="Coordinate in Reference Genome (Mb)", ylab="Read Coverage Depth", main=my_title )
#	arrows(pos[del_start],maxy,pos[del_start], 0, length=0, col="black", lwd=3)
#	arrows(pos[del_end],maxy,pos[del_end], 0, length=0, col="black", lwd=3)
#	rect(pos[start], 0, pos[del_start], maxy, col="grey80", lty=0)
#	rect(pos[del_end], 0, pos[end], maxy, col="grey80", lty=0)

	lines(pos[start:end],cov\$redundant_top_cov[start:end]+cov\$redundant_bot_cov[start:end], type="s", col="red", lty="solid", lwd=2 )
	lines(pos[start:end],cov\$redundant_top_cov[start:end], type="s", col="yellow", lty="solid")
	lines(pos[start:end],cov\$redundant_bot_cov[start:end], type="s", col="orange", lty="solid")
	
	lines(pos[start:end],cov\$unique_top_cov[start:end]+cov\$unique_bot_cov[start:end], type="s", col="blue", lty="solid", lwd=2 )
	lines(pos[start:end],cov\$unique_top_cov[start:end], type="s", col="cyan", lty="solid")
	lines(pos[start:end],cov\$unique_bot_cov[start:end], type="s", col="purple", lty="solid")
		
	mtext(genes, side=1, col="blue", cex=0.7, outer=TRUE, line=-2)
}

downsample <- function(pos, factor)
{
	new_pos <- 1:trunc(length(pos)/factor)
	for (i in 1:length(new_pos)) { 
			new_pos[i] = pos[i*factor]
	}
	return(new_pos);
}

cov<-read.table("$coverage_file", header=T);
pos<-1:dim(cov)[1];
pos = pos-0.5;
pos = pos / 1000000;
pdf("$output_file", $pdf_options)
pos_down <- downsample(pos, $downsample)
cov_down = c();
cov_down\$unique_top_cov <- downsample(cov\$unique_top_cov, $downsample)
cov_down\$unique_bot_cov <- downsample(cov\$unique_bot_cov, $downsample)
cov_down\$redundant_top_cov <- downsample(cov\$redundant_top_cov, $downsample)
cov_down\$redundant_bot_cov <- downsample(cov\$redundant_bot_cov, $downsample)
END
	
	foreach my $interval (@$intervals)
	{
		# $interval->{genes} =~ s/\[/\\\[/;
		# $interval->{genes} =~ s/\]/\\\]/;
		$interval->{genes} = "" if (!defined $interval->{genes});
		$interval->{genes} =~ s/(.{80,}?)\s+/$1\\n/g;
		print R_SCRIPT "plot_coverage(pos_down,cov_down,$interval->{start},$interval->{end}, $downsample, \"Genome Coordinates $interval->{start}-$interval->{end}\", \"$interval->{genes}\")\n";
	}	
	print R_SCRIPT "dev.off()\n";

	close R_SCRIPT;
	
	my $command = "R --vanilla < $$.coverage_graph.r_script &> $$.coverage_graph.r_output";
	my $res = system $command;
	die "Error running command: $command" if ($res);
	
	unlink "$$.coverage_graph.r_script" if (!$verbose);
	unlink "$$.coverage_graph.r_output" if (!$verbose);
}


sub graph_coverage
{
	my ($intervals, $coverage_file, $output_deletions_file, $output_overview_file, $downsample, $max_coverage) = @_;
	
	$max_coverage = "max(cov\$redundant_top_cov[start:end]+cov\$redundant_bot_cov[start:end], cov\$unique_top_cov[start:end]+cov\$unique_bot_cov[start:end])" if (!defined $max_coverage);
	my $max_unique_coverage = "max(cov\$unique_top_cov[start:end]+cov\$unique_bot_cov[start:end])";

	open R_SCRIPT, ">$$.coverage_graph.r_script";
	print R_SCRIPT <<END;
plot_coverage <- function(pos, cov, del_start, del_end, my_title, genes)
{
	
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

	options(scipen = 15)

	par(mar=c(8.1,4.1,4.1,2.1));
	plot(0:10, 0:10, type="n", lty="solid", ylim=c(0, maxy), xlim=c(pos[start], pos[end]), lwd=2, xlab="Coordinate in Reference Genome", ylab="Read Coverage Depth", main=my_title )
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
plot_overview <- function(pos, cov, start, end, my_title, genes)
{
	maxy=$max_unique_coverage;
	par(mar=c(8.1,4.1,4.1,2.1));
	plot(pos[start:end],cov\$redundant_top_cov[start:end]+cov\$redundant_bot_cov[start:end], type="s", col="red", lty="solid", lwd=1, ylim=c(0,maxy), xlab="Coordinate in Reference Genome", ylab="Read Coverage Depth", main=my_title )
	mtext(genes, side=1, col="blue", cex=0.7, outer=TRUE, line=-2)	
	lines(pos[start:end],cov\$unique_top_cov[start:end]+cov\$unique_bot_cov[start:end], type="s", col="blue", lty="solid", lwd=1 )
}
plot_cnv <- function(pos, cov, start, end, my_title, genes)
{
# Need to re-implement
#	maxy=max(cov\$V5[start:end]);
#	par(mar=c(8.1,4.1,4.1,2.1));
#	plot(pos[start:end],cov\$V5[start:end], type="s", col="green", lty="solid", lwd=1, ylim=c(0,maxy), xlab="Coordinate in Reference Genome", ylab="Predicted Copy Number", main=my_title )
#	mtext(genes, side=1, col="blue", cex=0.7, outer=TRUE, line=-2)	
}

downsample <- function(pos, factor)
{
	new_pos <- pos[1]
	for (i in 2:dim(cov)[1]) { 
		if ((i-1) %% factor == 0) { 
			new_pos = rbind(new_pos,pos[i])
		} 
	}
	return(new_pos);
}

cov<-read.table("$coverage_file", header=T);
pos<-1:dim(cov)[1];
pos = pos-0.5;
pos = pos / 1000000;
END
	
	if (!$output_path)
	{
		print R_SCRIPT "pdf(\"$output_deletions_file\", $pdf_options)\n";
	}
		
	my $i=1;
	foreach my $interval (@$intervals)
	{
		# $interval->{genes} =~ s/\[/\\\[/;
		# $interval->{genes} =~ s/\]/\\\]/;
		$interval->{genes} = "" if (!defined $interval->{genes});
		$interval->{genes} =~ s/(.{80,}?)\s+/$1\\n/g;
		
		print R_SCRIPT "pdf(\"$output_path/$i.pdf\")\n" if ($output_path);
		print R_SCRIPT "plot_coverage(pos,cov,$interval->{start},$interval->{end}, \"Predicted Deletion $interval->{start}-$interval->{end}\", \"$interval->{genes}\")\n";
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
pos100 <- downsample(pos, $downsample)
cov100 = c()
cov100\$unique_top_cov <- downsample(cov\$unique_top_cov, $downsample)
cov100\$unique_bot_cov <- downsample(cov\$unique_bot_cov, $downsample)
cov100\$redundant_top_cov <- downsample(cov\$redundant_top_cov, $downsample)
cov100\$redundant_bot_cov <- downsample(cov\$redundant_bot_cov, $downsample)
#re-implement
#cov100\$V5 <- filter(cov\$V5, rep(1/5,5), method="convolution", sides=2, circular=TRUE)
#cov100\$V5 <- downsample(cov100\$V5, $downsample)
END

		print R_SCRIPT "pdf(\"$output_path/overview.pdf\", $pdf_options)\n" if ($output_path);
		print R_SCRIPT "plot_overview(pos100,cov100,1,dim(pos100)[1], \"Genome Overview\", \"\")\n";
		print R_SCRIPT "dev.off()\n" if ($output_path);

		print R_SCRIPT "pdf(\"$output_path/cnv.pdf\", $pdf_options)\n" if ($output_path);
		print R_SCRIPT "plot_cnv(pos100,cov100,1,dim(pos100)[1], \"Copy Number Predictions\", \"\")\n";
		print R_SCRIPT "dev.off()\n" if ($output_path);

	}

	close R_SCRIPT;
	
	my $command = "R --vanilla < $$.coverage_graph.r_script &> $$.coverage_graph.r_output";
	my $res = system $command;
	die "Error running command: $command" if ($res);
	
	unlink "$$.coverage_graph.r_script" if (!$verbose);
	unlink "$$.coverage_graph.r_output" if (!$verbose);
}
