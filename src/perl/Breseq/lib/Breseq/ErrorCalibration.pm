###
# Pod Documentation
###

=head1 NAME

Breseq::Fastq.pm

=head1 SYNOPSIS

Module for reading and writing fastq files more rapidly than BioPerl.

=head1 DESCRIPTION

=head1 AUTHOR

Jeffrey Barrick
<jbarrick@msu.edu>

=head1 COPYRIGHT

Copyright 2009.  All rights reserved.

=cut

###
# End Pod Documentation
###

use strict;

package Breseq::ErrorCalibration;

use Data::Dumper;

###
# Global Variables
###

our @base_list = ('A', 'T', 'C', 'G', '.');
our $minimum_error_probability = 1E-256;

=head2 alignment_database_calibrate

 Title   : alignment_database_calibrate
 Usage   : 
 Function: 
 Returns : 

=cut

sub count
{
	my $verbose = 0;
	my ($settings, $summary, $ref_seq_info) = @_;

	## for this error method, we don't need to iterate through the sequence and count...
	return if ($settings->{error_model_method} eq 'QUALITY');

	my $reference_fasta_file_name = $settings->file_name('reference_fasta_file_name');
	my $reference_bam_file_name = $settings->file_name('reference_bam_file_name');

	## check to see if error_count++ is available:
	my $error_count = $settings->ctool('error_count');
	if ($error_count && !$settings->{perl_error_count}) 
	{
		my $coverage_fn = $settings->file_name('unique_only_coverage_distribution_file_name', {'@'=>""});
		my $outputdir = `dirname $coverage_fn`;
		chomp $outputdir; $outputdir .= '/';
		my $readfiles = join(" --readfile ", $settings->read_files);
		my $cmdline = "$error_count --bam $reference_bam_file_name --fasta $reference_fasta_file_name --output $outputdir --readfile $readfiles";
		print STDERR "$cmdline\n";
		system ($cmdline) == 0
			or die $!;
		return; # error_count++ worked, so we're all done here.
	}

	my $bam = Bio::DB::Sam->new(-fasta => $reference_fasta_file_name, -bam => $reference_bam_file_name);
	##my @seq_ids = $bam->seq_ids;
	my @seq_ids = @{$ref_seq_info->{seq_ids}};
		
	## populated by pileup
	my $error_hash = [];		#list by fastq file index
	#my $complex_error_hash = [];		#list by fastq file index

	foreach my $seq_id (@seq_ids)
	{							
		my $sequence_length = $bam->length($seq_id);
		print STDERR "  REFERENCE: $seq_id\n";
		print STDERR "  LENGTH: $sequence_length\n";
		
		## populated by pileup
		my $unique_only_coverage = [];
		my $ref_seq_string = $ref_seq_info->{ref_strings}->{$seq_id};
				
        my $pileup_function = sub {
        	my ($seq_id,$pos,$pileup) = @_;

			print STDERR "    POSITION:$pos\n" if ($pos % 10000 == 0);

			my @ref_base; #index is 'reversed'
			
			## temporary debug
			my $ref_length = length $ref_seq_string;
			
			if ($pos > $ref_length)
			{
				print $pileup->[0]->alignment->qname . "\n";
				die "Position outside of reference sequence $pos > $ref_length\n";
			}
			
       		#$ref_base[0] = $bam->segment($seq_id,$pos,$pos)->dna;
			#$ref_base[1] = $bam->segment($seq_id,$pos,$pos)->seq->revcom->seq;

			my $has_redundant_reads = 0;
			my $unique_coverage = 0;
         	ALIGNMENT: for my $p (@$pileup) 
			{
				my $a = $p->alignment;
				
				##
				# We don't need to process deleted positions
				# ... all information comes from adjacent positions
				##
				next if ($p->is_del);
				
				my $redundancy = $a->aux_get('X1');
				if ($redundancy > 1)
				{
					$has_redundant_reads = 1;
					next ALIGNMENT;
				}
				
				
				## this position has coverage b/c it is not deleted
				$unique_coverage++;
				
				## gather information about this position in read
				my $qseq = $a->qseq;
				my $qpos = $p->qpos; #0-indexed
				my $qscore = $a->qscore;
				my $indel = $p->indel;
				my $reversed = $a->reversed;
				my $query_end = $a->query->end-1; #0-indexed
				my $query_start = $a->query->start-1; #0-indexed
				my $fastq_file_index = $a->aux_get('X2');
				
				print $a->qname . " " . $pos . " " . $qpos . " " . " " . $indel . " " . $reversed . " " . $query_start . " " . $query_end . "\n" if ($verbose);
								
				##
				# In all that follows, be sure to keep track of strandedness of mutations!
				##
				

				# (1) base substitution or match
				#     e.g. 'AG' key for observing a G in read at a place where the reference has an A
				#     IMPROVE by keeping track of flanking base context and scores?
				#     this would, for example, penalize low scoring sequences more

				{
					my $base  = substr($qseq, $qpos,1);

					my $ref_base = substr $ref_seq_string, $pos-1, 1;

					if (($base ne 'N') && ($ref_base ne 'N'))
					{		
						$base = Breseq::Fastq::revcom($base) if ($reversed);
						$ref_base = Breseq::Fastq::revcom($ref_base) if ($reversed);
						my $key = $ref_base . $base; 
						$error_hash->[$fastq_file_index]->{$qscore->[$qpos]}->{$key}++;
						print "$key $qscore->[$qpos]\n" if ($verbose);
					}
				}
				
				# the next base also matches 
				# (1) base substitution or match
				#     e.g. '..' key indicating an observation of a "non-gap, non-gap"
				#     quality score is of the second non-gap in the pair
				if ($indel == 0)
				{
					## don't count past last match position
					if ($qpos < $query_end)
					{	
						my $mqpos = $qpos + 1 - $reversed;	
						my $base  = substr($qseq, $mqpos,1);
						
						my $mrpos = $pos + 1 - $reversed;	
						my $ref_base = substr $ref_seq_string, $mrpos-1, 1;
						
						if (($base ne 'N') && ($ref_base ne 'N'))
						{
							my $key = '.' . '.'; 
							$error_hash->[$fastq_file_index]->{$qscore->[$mqpos]}->{$key}++;
							print "$key $qscore->[$mqpos]\n" if ($verbose);						
						}
					}	
				}
				
				# there is a deletion of EXACTLY one base in the read relative to the reference before the next read base
				# (2) deletion in read relative to reference
				#     e.g. 'A.' key for observing nothing in a read at a position where the reference has an A
				#     quality score is of the next non-gap base in the read
				elsif ($indel == -1)
				{
					my $mqpos = $qpos + 1 - $reversed;
					my $base  = substr($qseq, $mqpos,1);

					## the reference base opposite the deletion is really the NEXT base
					my $mrpos = $pos + 1;
					my $ref_base = substr $ref_seq_string, $mrpos-1, 1;

					if (($base ne 'N') && ($ref_base ne 'N'))
					{
						$ref_base = Breseq::Fastq::revcom($ref_base) if ($reversed);
						my $key = $ref_base . '.'; 
						$error_hash->[$fastq_file_index]->{$qscore->[$mqpos]}->{$key}++;
						print "$key $qscore->[$mqpos]\n" if ($verbose);												
					}
				}
				
				# there is an insertion of EXACTLY one base in the read relative to the reference before the next reference base
				# (3) insertion in read relative to reference
				#     e.g. '.A' key for observing an A in a read at a position where the reference has no base
				#     quality score is that of the observed inserted base
				elsif ($indel == +1)
				{
					my $mqpos = $qpos + 1;
				
					if (($mqpos <= $query_end) && ($mqpos >= $query_start))
					{
						my $base  = substr($qseq, $mqpos,1);
				
						if ($base ne 'N')
						{	
							$base = Breseq::Fastq::revcom($base) if ($reversed);
							my $key = '.' . $base; 
							$error_hash->[$fastq_file_index]->{$qscore->[$mqpos]}->{$key}++;
							print "$key $qscore->[$mqpos]\n" if ($verbose);												
						}	
					}
				}
			} #end ALIGNMENT

			# record unique only coverage
			$unique_only_coverage->[$unique_coverage]++ if (!$has_redundant_reads);
		}; #end $pileup_function

        $bam->pileup($seq_id,$pileup_function);
        #$bam->pileup("$seq_id:1-10000", $pileup_function);
			
		## save the unique only coverage distribution
		my $this_unique_only_coverage_distribution_file_name = $settings->file_name('unique_only_coverage_distribution_file_name', {'@'=>$seq_id});
		save_unique_coverage_distribution_file($this_unique_only_coverage_distribution_file_name, $unique_only_coverage);
	}
	
	foreach my $read_file ($settings->read_files)
	{
		my $fastq_file_index = $settings->read_file_to_fastq_file_index($read_file);
		my $this_error_counts_file_name = $settings->file_name('error_counts_file_name', {'#'=> $read_file});
		save_error_file($this_error_counts_file_name, $error_hash->[$fastq_file_index]);
	}
	
	# foreach my $read_file ($settings->read_files)
	# {
	# 	my $fastq_file_index = $settings->read_file_to_fastq_file_index($read_file);
	# 	my $this_complex_error_counts_file_name = $settings->file_name('complex_error_counts_file_name', {'#'=> $read_file});
	# 	save_complex_error_file($this_complex_error_counts_file_name, $complex_error_hash->[$fastq_file_index]);
	# }
	
}


sub save_unique_coverage_distribution_file
{
	my ($unique_only_coverage_file_name, $unique_only_coverage_list_ref) = @_;
	open COV, ">$unique_only_coverage_file_name" or die;
	print COV "coverage\tn\n";
	# print out coverage, but ZERO OUT observations of zero
	# because these may be bona fide deletions
	for (my $i=1; $i<scalar @$unique_only_coverage_list_ref; $i++)
	{
		my $cov = $unique_only_coverage_list_ref->[$i];
		$cov = 0 if (!defined $cov); #some may be undefined, they mean zero
		print COV "$i\t$cov\n";
	}		
	close COV;
}

sub analyze_unique_coverage_distributions
{
	my ($settings, $summary, $ref_seq_info) = @_;

	foreach my $seq_id (@{$ref_seq_info->{seq_ids}})
	{
		my $this_unique_only_coverage_plot_file_name = $settings->file_name('unique_only_coverage_plot_file_name', {'@'=>$seq_id});
		my $this_unique_only_coverage_distribution_file_name = $settings->file_name('unique_only_coverage_distribution_file_name', {'@'=>$seq_id});
		
		analyze_unique_coverage_distribution_file_using_R($settings, $this_unique_only_coverage_distribution_file_name, $this_unique_only_coverage_plot_file_name, $seq_id, $summary);
	}
}

## uses R to fit negative binomial distribution to coverage distribution
## and adds this information to the $summary
sub analyze_unique_coverage_distribution_file_using_R
{
	my ($settings, $unique_only_coverage_file_name, $unique_only_coverage_plot_file_name, $seq_id, $summary) = @_;

	##initialize summary information
	$summary->{unique_coverage}->{$seq_id}->{nbinom_size_parameter} = 'NA';
	$summary->{unique_coverage}->{$seq_id}->{nbinom_prob_parameter} = 'NA'; 
	$summary->{unique_coverage}->{$seq_id}->{average} = 1;
	$summary->{unique_coverage}->{$seq_id}->{variance} = 'NA';
	$summary->{unique_coverage}->{$seq_id}->{dispersion} = 'NA';
	$summary->{unique_coverage}->{$seq_id}->{deletion_coverage_propagation_cutoff} = 5;

## method to avoid R code...
## we bail early if in QUALITY mode
	return if ($settings->{error_model_method} eq 'QUALITY');

	my $unique_coverage_distribution_r_script_file_name = $settings->file_name('unique_coverage_distribution_r_script_file_name', {'@' => $seq_id});
	open RSCRIPT, ">$unique_coverage_distribution_r_script_file_name" or die "Could not create file: $unique_coverage_distribution_r_script_file_name\n";

	### Define various coverage thresholds...
	my $sequence_length = $summary->{sequence_conversion}->{reference_sequences}->{$seq_id}->{length};

	### DELETION PROPAGATION CUTOFF
	## One-tailed test p=0.05, Bonferroni correction
	# my $del_propagation_pr_cutoff = 0.05 / $sequence_length;
	
	## One-tailed test p=0.01, no Bonferroni correction
	#my $del_propagation_pr_cutoff = 0.01;
	
	## We really want somewhere between these two, try this...
	my $deletion_propagation_pr_cutoff = 0.05 / sqrt($sequence_length);

	### NEW JUNCTION CUTOFF
	## Arbitrary value that seems to work....
	my $new_junction_pr_cutoff = 1/$sequence_length; # *0.05;
	

print RSCRIPT <<EOF;
#load data
X<-read.table("$unique_only_coverage_file_name", header=T)

#create the distribution vector and fit
Y<-rep(X\$coverage, X\$n)
m<-mean(Y)
v<-var(Y)
D<-v/m

## we find the range to fit using a moving average to avoid noise problems
ma5 = c(1, 1, 1, 1, 1)/5;
X\$ma = filter(X\$n, ma5)

i<-0
max_n <- 0;
max_i <- i;
for (i in trunc(m/4):length(X\$ma))
{		
	if (!is.na(X\$ma[i]) && X\$ma[i] > max_n)
	{
		max_n = X\$ma[i];
		max_i = i;
	}
}

#censor data on the right and left of the biggest maximum
coverage_factor <- 0.25;

i<-max_i
while (i >= 1 && X\$ma[i] && X\$ma[i]>coverage_factor*max_n)
{	
	i <- i-1;
}
start_i = i;

i<-length(X\$ma);
#i<-max_i
#while (i <= length(X\$ma) && X\$ma[i]>coverage_factor*max_n)
#{		
#	i <- i+1;
#}
end_i <-i

#now maximum likelihood find the best dispersion parameter

inner_total<-0;
for (i in start_i:end_i)
{
	inner_total = inner_total + X\$n[i]; 
}
total_total<-sum(X\$n);

i<-max_i
while (i <= length(X\$n) && X\$n[i]>0.01*max_n)
{		
	i <- i+1;
}
graph_end_i <-i

f_nb <- function(par) {

	mu = par[1];
	size = par[2];
	
	dist<-c()
	total <- 0;
	for (i in start_i:end_i)
	{	
		dist[i] <- dnbinom(i, size=size, mu=mu);
		total <- total + dist[i] 
	}
	#print (mu, size)
	
 	l <- 0;
	for (i in start_i:end_i)
	{
		l <- l + ((X\$n[i]/inner_total)-(dist[i]/total))^2;
	}
	l;
}

nb_fit<-nlm(f_nb, c(m,1) )
nb_fit_mu = nb_fit\$estimate[1];
nb_fit_size = nb_fit\$estimate[2];

print(nb_fit_size);
print(nb_fit_mu);

fit_nb = dnbinom(0:max(X\$coverage), mu = nb_fit_mu, size=nb_fit_size)*total_total;


f_p <- function(par) {

	lambda = par[1];
	
	dist<-c()
	total <- 0;
	for (i in start_i:end_i)
	{	
		dist[i] <- dpois(i, lambda=lambda);
		total <- total + dist[i] 
	}
	#print (total)
	
 	l <- 0;
	for (i in start_i:end_i)
	{
		l <- l + ((X\$n[i]/inner_total)-(dist[i]/total))^2;
	}
	l;
}

p_fit<-nlm(f_p, c(m))
p_fit_lambda = nb_fit\$estimate[1];
fit_p<-dpois(0:max(X\$coverage), lambda = p_fit_lambda)*total_total;


my_pch = 21
my_col = "black";

pdf("$unique_only_coverage_plot_file_name", height=6, width=7)
par(mar=c(5,5,3,3));
plot(0:10, 0:10, type="n", lty="solid", ylim=c(0, max(X\$n, fit_p, fit_nb))*1.05, xlim=c(0, graph_end_i), lwd=1, xaxs="i", yaxs="i", axes=F, las=1, xlab="", ylab="", cex.lab=1.2, cex.axis=1.2 )
box()
#axis(2, cex.lab=1.2, las=1, cex.axis=1.2, labels=T, at=(0:6)*50000)
axis(2, cex.lab=1.2, las=1, cex.axis=1.2, labels=T)
axis(1, cex.lab=1.2, cex.axis=1.2, labels=T)

#graph the coverage as points
points(X\$coverage, X\$n, pch=my_pch, col=my_col, bg="white", cex=1.2)

#graph the poisson fit
lines(0:max(X\$coverage), fit_p, lwd=3, lty="22", col="black");

#graph the negative binomial fit
lines(0:max(X\$coverage), fit_nb, lwd=3, col="black");

dev.off()


#print out some statistics
print(m)
print(v)
print(D)

deletion_propagation_coverage = qnbinom($deletion_propagation_pr_cutoff, size = nb_fit_size, mu = nb_fit_mu)
print(deletion_propagation_coverage)

new_junction_propagation = qnbinom($new_junction_pr_cutoff, size = nb_fit_size, mu = nb_fit_mu)
print(new_junction_propagation)
EOF

	close RSCRIPT;
	
	Breseq::Shared::system("R --vanilla < $unique_coverage_distribution_r_script_file_name > $unique_coverage_distribution_r_script_file_name\.out");

	open ROUT, "<$unique_coverage_distribution_r_script_file_name\.out" or die;
	my @lines = <ROUT>;
	close ROUT;
	chomp @lines;
	
	@lines = grep s/^\[1\]\s+//, @lines;
#	print Dumper(@lines);
	
	#First two lines are negative binomial parameters.
	#Next three lines are average, standard deviation, and index of overdispersion

	#Put these into summary
	$summary->{unique_coverage}->{$seq_id}->{nbinom_size_parameter} = $lines[0];
		#prob = size/(size + mu)
	$summary->{unique_coverage}->{$seq_id}->{nbinom_prob_parameter} = $lines[0] / ($lines[0] + $lines[1]); 
	$summary->{unique_coverage}->{$seq_id}->{average} = $lines[2];
	$summary->{unique_coverage}->{$seq_id}->{variance} = $lines[3];
	$summary->{unique_coverage}->{$seq_id}->{dispersion} = $lines[4];

	$summary->{unique_coverage}->{$seq_id}->{deletion_coverage_propagation_cutoff} = $lines[5];
	$summary->{unique_coverage}->{$seq_id}->{junction_coverage_cutoff} = $lines[6];
}

sub save_error_file
{
	my ($file_name, $error_hash_ref) = @_;
	
	#print out a table of errors stratified by quality scores
	open ERR, ">$file_name" or die;
	
	##
	## Note that we are ignoring 'N' bases by NOT PRINTING THEM OUT
	## Later when encountering these they should be SKIPPED
	##
	
	## Create header list which is the same for all files		
	my @bases = ('A', 'T', 'C', 'G', '.');
	my @header_list = ('quality');
	foreach my $base_1 (@bases)
	{
		foreach my $base_2 (@bases)
		{
			push @header_list, "$base_1$base_2";
		}
	}

	print ERR join("\t", @header_list). "\n";
	
	foreach my $quality (sort { -($a<=>$b) } keys %{$error_hash_ref})
	{
		my @line_list;
		push @line_list, $quality;
		
		foreach my $base_1 (@bases)
		{
			foreach my $base_2 (@bases)
			{
				my $val = $error_hash_ref->{$quality}->{"$base_1$base_2"};
				$val = 0 if (!defined $val);
				push @line_list, $val;
			}
		}
		print ERR join("\t", @line_list). "\n";
	}
	close ERR;
}

sub save_complex_error_file
{
	my ($file_name, $error_hash_ref) = @_;
	
	#print out a table of errors stratified by quality scores
	open ERR, ">$file_name" or die;
	
	##
	## Note that we are ignoring 'N' bases by NOT PRINTING THEM OUT
	## Later when encountering these they should be SKIPPED
	##
	
	## Create header list which is the same for all files		
	my @bases = ('A', 'T', 'C', 'G', '.');
	my @header_list = ('quality', 'neighborhood');
	foreach my $base_1 (@bases)
	{
		foreach my $base_2 (@bases)
		{
			push @header_list, "$base_1$base_2";
		}
	}

	print ERR join("\t", @header_list). "\n";
	
	foreach my $neighborhood (sort { -($a<=>$b) } keys %{$error_hash_ref})
	{
		foreach my $quality (sort { -($a<=>$b) } keys %{$error_hash_ref->{$neighborhood}})
		{
			my @line_list;
			push @line_list, $quality, $neighborhood;
		
			foreach my $base_1 (@bases)
			{
				foreach my $base_2 (@bases)
				{
					my $val = $error_hash_ref->{$neighborhood}->{$quality}->{"$base_1$base_2"};
					$val = 0 if (!defined $val);
					push @line_list, $val;
				}
			}
			print ERR join("\t", @line_list). "\n";
		}
	}
	close ERR;
}

sub load_error_file
{
	my ($file_name) = @_;
	
	open ERR, "<$file_name" or die "Could not open file: $file_name";
	my $header_line = <ERR>;
	
	chomp $header_line;
	my @header_list = split "\t", $header_line;
		
	my $error_hash;
	while (my $this_line = <ERR>)
	{
		#print "$this_line";
		chomp $this_line;
		my @line_list = split "\t", $this_line;
		my $this_quality = $line_list[0];
		for (my $i=1; $i<scalar @header_list; $i++)
		{			
			$error_hash->{$this_quality}->{$header_list[$i]} = $line_list[$i];
		}
	}
	close ERR;
	
	return $error_hash;
}

sub load_error_rates
{
	my ($settings, $summary, $ref_seq_info) = @_;
	my @seq_ids = @{$ref_seq_info->{seq_ids}};
		
	my @error_rates_list;
	foreach my $read_file ($settings->read_files)
	{
		my $fastq_file_index = $settings->read_file_to_fastq_file_index($read_file);
		my $this_error_rates_file_name = $settings->file_name('error_rates_file_name', {'#' => $read_file});
		$error_rates_list[$fastq_file_index] = load_error_file($this_error_rates_file_name);
	}
	return \@error_rates_list;
}

sub log10_error_rates
{
	my ($error_rates) = @_;
		
	## precalculate log10 probabilities and not probabilities 
	## results in a significant speed-up of calculations	
	my $log10_correct_rates;
	my $log10_error_rates;
	my $log10_random_rates;
	foreach (my $i=0; $i<scalar @$error_rates; $i++)
	{
		foreach my $q (sort keys %{$error_rates->[$i]})
		{
			foreach my $base_1 (@base_list)
			{				
				foreach my $base_2 (@base_list)
				{
					my $c = $base_1 . $base_2;
					my $pr = $error_rates->[$i]->{$q}->{$c};	
					$pr = 1E-256 if ($pr == 0);
					my $one_minus_pr = 1-$pr;
					$one_minus_pr = 1E-256 if ($one_minus_pr == 0);
					
					$log10_correct_rates->[$i]->{$q}->{$c} = log($pr) / log(10);
					$log10_error_rates->[$i]->{$q}->{$c} = log($one_minus_pr) / log(10);
					$log10_random_rates->[$i]->{$q}->{$c} = log($pr/0.25) / log(10);
				}
			}
		}
	}
	
	return ($log10_correct_rates, $log10_error_rates, $log10_random_rates);
}


sub error_counts_to_error_rates
{
	my ($settings, $summary, $ref_seq_info) = @_;
		
	my @seq_ids = @{$ref_seq_info->{seq_ids}};
		
	#for each error file: load results, convert, and save
	foreach my $read_file ($settings->read_files)
	{
		my $fastq_file_index = $settings->read_file_to_fastq_file_index($read_file);
		my $this_error_counts_file_name = $settings->file_name('error_counts_file_name', {'#' => $read_file} );
		my $this_error_rates_file_name = $settings->file_name('error_rates_file_name', {'#' => $read_file} );
		my $this_error_rates_plot_file_name = $settings->file_name('error_rates_plot_file_name', {'#' => $read_file } );
		
		#two options:
		#(1) Calculate these probabilities directly
	
		if ($settings->{error_model_method} eq 'QUALITY')
		{			
			my $error_rates_hash;
			$error_rates_hash = error_rates_calculated_from_quality();
			save_error_file($this_error_rates_file_name, $error_rates_hash);
		}
		elsif ($settings->{error_model_method} eq 'EMPIRICAL')
		{	
			my $error_rates_hash;
			$error_rates_hash = error_counts_to_error_rates_empirical($this_error_counts_file_name);
			save_error_file($this_error_rates_file_name, $error_rates_hash);			
		}
		#(2) Calculate using log-linear model (assumes base quality calibration is correct!)
		elsif ($settings->{error_model_method} eq 'FIT')
		{			
			error_counts_to_error_rates_using_R($settings, $this_error_counts_file_name, $this_error_rates_file_name, $this_error_rates_plot_file_name, $read_file);
		}
		else
		{
			die "Unknown Error Model Method:$settings->{error_model_method}\n";
		}
	}
}

sub error_rates_calculated_from_quality
{
	my ($quality_style, $relative_indel) = @_;
	
	##currently assumed...
	$quality_style = 'ILLUMINA-1.3' if (!defined $quality_style);
	
	$relative_indel = 0.01 if (!defined $relative_indel);
	
	## expect current Illumina phred format for this to work??
	my $error_rates_hash;
	for (my $s=0; $s<=62; $s++)
	{
		my $total_error_rate = 10**(-$s/10);
		
		my $no_error = 1 - $total_error_rate;
		my $each_rate = $total_error_rate / (3 + $relative_indel);
		
		my @bases = ('A', 'C', 'T', 'G', '.');
		foreach my $b1 (@bases)
		{
			foreach my $b2 (@bases)
			{
				if ($b1 eq $b2)
				{
					$error_rates_hash->{$s}->{"$b1$b2"} = $no_error;
				}
				elsif (($b1 eq '.') || ($b2 eq '.'))
				{
					$error_rates_hash->{$s}->{"$b1$b2"} = $each_rate * $relative_indel;
				}
				else
				{
					$error_rates_hash->{$s}->{"$b1$b2"} = $each_rate;
				}
			}
		}
	}
	
	return $error_rates_hash;
}

sub error_counts_to_error_rates_empirical
{
	my ($this_error_counts_file_name) = @_;
	my $error_counts_hash = load_error_file($this_error_counts_file_name);

	my $error_rates_hash;
	foreach my $q (sort keys %$error_counts_hash)
	{
		#find total with quality
		my $total_bases_of_quality = 0;
		my @quality_list = sort keys %{$error_counts_hash->{$q}};
			
		my @bases = ('A', 'C', 'T', 'G', '.');
		foreach my $b1 (@bases)
		{
			my $total_this_base = 0;
			foreach my $b2 (@bases)
			{
				$total_this_base += $error_counts_hash->{$q}->{"$b1$b2"};
			}
			
			foreach my $b2 (@bases)
			{
				$error_rates_hash->{$q}->{"$b1$b2"} = $error_counts_hash->{$q}->{"$b1$b2"} / $total_this_base;
				$error_rates_hash->{$q}->{"$b1$b2"} = $minimum_error_probability if ($error_rates_hash->{$q}->{"$b1$b2"} < $minimum_error_probability);
			}			
		}
	}
	
	return $error_rates_hash;
}

sub error_counts_to_error_rates_using_R
{
	my ($settings, $input_file_name, $output_file_name, $plot_file_name, $read_file) = @_;
	
	my $error_rates_r_script_file_name = $settings->file_name('error_rates_r_script_file_name', {'#' => $read_file});
	open RSCRIPT, ">$error_rates_r_script_file_name" or die "Could not create file: $error_rates_r_script_file_name\n";
			
	my @bases = ('A', 'C', 'T', 'G', '.');
	
	my $regression_string = "X\$quality";
	print RSCRIPT <<EOF;	
library(nnet);
X <- read.table("$input_file_name", header = TRUE, sep = "\t")
X <- subset(X, quality>=$settings->{base_quality_cutoff});
Y <- data.frame(quality=X\$quality)
EOF

	foreach my $base_1 (@bases)
	{		
		my @base_list;
		foreach my $base_2 (@bases)
		{
			push @base_list, "$base_1$base_2";
		}
		
		my $response_string = "cbind(" . join(",", @base_list) . ")";
		my $subset_string = join("+", @base_list);
		
		
		print RSCRIPT <<EOF;	
myfit <- multinom( $response_string ~ $regression_string, X, na.action=na.fail, subset=($subset_string != 0), maxit=10000, MaxNWts=10000)
Z = predict(myfit, X, type="probs")

min_value <- $minimum_error_probability;
for (i in 1:dim(Z)[1])
{
	for (j in 1:dim(Z)[2])
	{
		if (Z[i,j] < min_value)
		{
			Z[i,j] <- min_value;
		}
	}
}

Y = cbind(Y,Z)
EOF
	}
	print RSCRIPT <<EOF;	
write.table(Y, file = "$output_file_name", sep = "\t", quote=F, row.names=F)

#create empty plot
pdf("$plot_file_name", height=6, width=9)

B<-as.matrix(X)
print(X)

A_tot = X\$AA+X\$AT+X\$AC+X\$AG+X\$A.;
A_tot[A_tot==0] = 1;
X\$AA = X\$AA / A_tot
X\$AT = X\$AT / A_tot
X\$AC = X\$AC / A_tot
X\$AG = X\$AG / A_tot
X\$A. = X\$A. / A_tot

C_tot = X\$CA+X\$CT+X\$CC+X\$CG+X\$C.;
C_tot[C_tot==0] = 1;
X\$CA = X\$CA / C_tot
X\$CT = X\$CT / C_tot
X\$CC = X\$CC / C_tot
X\$CG = X\$CG / C_tot
X\$C. = X\$C. / C_tot

T_tot = X\$TA+X\$TT+X\$TC+X\$TG+X\$T.;
T_tot[T_tot==0] = 1;
X\$TA = X\$TA / T_tot
X\$TT = X\$TT / T_tot
X\$TC = X\$TC / T_tot
X\$TG = X\$TG / T_tot
X\$T. = X\$T. / T_tot

G_tot = X\$GA+X\$GT+X\$GC+X\$GG+X\$G.;
G_tot[G_tot==0] = 1;
X\$GA = X\$GA / G_tot
X\$GT = X\$GT / G_tot
X\$GC = X\$GC / G_tot
X\$GG = X\$GG / G_tot
X\$G. = X\$G. / G_tot

._tot = X\$.A+X\$.T+X\$.C+X\$.G+X\$..;
._tot[._tot==0] = 1;
X\$.A = X\$.A / ._tot
X\$.T = X\$.T / ._tot
X\$.C = X\$.C / ._tot
X\$.G = X\$.G / ._tot
X\$.. = X\$.. / ._tot


X[X==0] = 0.00000001;

print(X)
X<-data.frame(X)
print(X)

new_plot <- function()
{
	#bottom, left, top, right
	par(mar=c(5,5,3,3));
	options(warn=-1)
	plot(0, type="n", lty="solid", log="y", ylim=c(0.00000001, 1), xlim=c(min(X\$quality), max(X\$quality)), lwd=1, axes=F, xlab="", ylab="", las=1, cex.lab=1.2, cex.axis=1.2 )
	options(warn=1)
	box()
	#y-axis
	axis(2, cex.lab=1.2, las=1, cex.axis=1.2, yaxs="i", at = c(0.00000001, 0.0000001, 0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 1), label = c(expression(10^-8), expression(10^-7), expression(10^-6), expression(10^-5), expression(10^-4), expression(10^-3), expression(10^-2), expression(10^-1), expression(1)))
	title(ylab="error rate", mgp = c(3.5, 1, 0), cex.lab=1.2)
	#x-axis
	axis(1, cex.lab=1.2, cex.axis=1.2, xaxs="i")
	title(xlab="base quality score", mgp = c(3, 1, 0), cex.lab=1.2)

}

A_col = "green";
A_pch = 15;
C_col = "red"; 
C_pch = 15;
T_col = "blue"; 
T_pch = 15;
G_col = "black"; 
G_pch = 15;
._col = "orange"; 
._pch = 15;

new_plot()
lines(Y\$quality, Y\$AC, col=C_col)
points(X\$quality, X\$AC, pch=C_pch, col=C_col)

lines(Y\$quality, Y\$AT, col=T_col)
points(X\$quality, X\$AT, pch=T_pch, col=T_col)

lines(Y\$quality, Y\$AG, col=G_col)
points(X\$quality, X\$AG, pch=G_pch, col=G_col)

lines(Y\$quality, Y\$A., col=._col)
points(X\$quality, X\$A., pch=._pch, col=._col)

new_plot()
lines(Y\$quality, Y\$TA, col=A_col)
points(X\$quality, X\$TA, pch=A_pch, col=A_col)

lines(Y\$quality, Y\$TC, col=C_col)
points(X\$quality, X\$TC, pch=C_pch, col=C_col)

lines(Y\$quality, Y\$TG, col=G_col)
points(X\$quality, X\$TG, pch=G_pch, col=G_col)

lines(Y\$quality, Y\$T., col=._col)
points(X\$quality, X\$T., pch=._pch, col=._col)

new_plot()
lines(Y\$quality, Y\$CA, col=A_col)
points(X\$quality, X\$CA, pch=A_pch, col=A_col)

lines(Y\$quality, Y\$CT, col=T_col)
points(X\$quality, X\$CT, pch=T_pch, col=T_col)

lines(Y\$quality, Y\$CG, col=G_col)
points(X\$quality, X\$CG, pch=G_pch, col=G_col)

lines(Y\$quality, Y\$C., col=._col)
points(X\$quality, X\$C., pch=._pch, col=._col)

new_plot()
lines(Y\$quality, Y\$GA, col=A_col)
points(X\$quality, X\$GA, pch=A_pch, col=A_col)

lines(Y\$quality, Y\$GT, col=T_col)
points(X\$quality, X\$GT, pch=T_pch, col=T_col)

lines(Y\$quality, Y\$GC, col=C_col)
points(X\$quality, X\$GC, pch=C_pch, col=C_col)

lines(Y\$quality, Y\$G., col=._col)
points(X\$quality, X\$G., pch=._pch, col=._col)

new_plot()
lines(Y\$quality, Y\$.A, col=A_col)
points(X\$quality, X\$.A, pch=A_pch, col=A_col)

lines(Y\$quality, Y\$.T, col=T_col)
points(X\$quality, X\$.T, pch=T_pch, col=T_col)

lines(Y\$quality, Y\$.C, col=C_col)
points(X\$quality, X\$.C, pch=C_pch, col=C_col)

lines(Y\$quality, Y\$.G, col=G_col)
points(X\$quality, X\$.G, pch=G_pch, col=G_col)

dev.off()

EOF
	
	close RSCRIPT;

	Breseq::Shared::system("R --vanilla < $error_rates_r_script_file_name > $error_rates_r_script_file_name\.out");
}

sub assign_mapping_quality_to_matches
{
	my ($match_list_ref, $read_seq, $ref_strings, $log10_random_rates, $settings) = @_;

#	print Dumper($log10_random_rates);


	##note error rates are specific to file it came from!
	
#	print Dumper($read_seq);

	#quality scores are the same for each match
	my $quals;
	@{$quals->{+1}} = Breseq::Fastq::quals($read_seq);
	@{$quals->{-1}} = reverse @{$quals->{+1}};
	
	## it would make sense to sort by length and mismatch before 
	## doing this and bail once they drop below a certain threshold
	
	## assign each match a score
	foreach my $m (@$match_list_ref)
	{		
#		print Dumper($m);		
#		print Dumper($ref_strings->{$m->{reference}});
#		print Dumper($read_seq);
		
		my $ref_string = substr $ref_strings->{$m->{reference}}, $m->{reference_start}-1, $m->{reference_end}-$m->{reference_start}+1;
		my $qry_string = substr $read_seq->{seq}, $m->{query_start}-1, $m->{query_end} - $m->{query_start}+1;
		$qry_string = Breseq::Fastq::revcom($qry_string) if ($m->{strand} == -1);
		
#		print "$ref_string\n$qry_string\n";
		
		## since we use the quality of the base 'after' the gap it actually matters if we are on the opposite strand for this
		my $qry_offset_for_gaps = ($m->{strand} == -1) ? 0 : 1;
			
		my ($r, $q) = AlignmentCorrection::add_gaps_to_alignment($m, $ref_string, $qry_string);
		my @ref = @$r;
		my @qry = @$q;
		my @qls = @{$quals->{$m->{strand}}};
	
#		print Dumper(\@ref, \@qry, \@qls);
	
		my $qual_pos = $m->{query_start}-1;
		my $read_pos = 0;
	
#		print "@ref\n@qry\n";
		## Loop through comparing each character
			
		my $quality = 0;  #log10 of quality
		while ($read_pos < scalar @qry)
		{	
			## have to check all alternative to ignore N's!!! (they don't count for or against)						
			if ($qry[$read_pos] eq 'N')
			{
			}
			elsif ($qry[$read_pos] eq '.')
			{
				my $base_key = "$ref[$read_pos]$qry[$read_pos+$qry_offset_for_gaps]";
				$quality += $log10_random_rates->{$qls[$qual_pos]}->{$base_key};
			}
			elsif ($ref[$read_pos] eq '.')
			{
				my $base_key = "$ref[$read_pos]$qry[$read_pos]";
				$quality += $log10_random_rates->{$qls[$qual_pos]}->{$base_key};
			}
			elsif ($qry[$read_pos] ne $ref[$read_pos])
			{
				my $base_key = "$ref[$read_pos]$qry[$read_pos]";				
				$quality += $log10_random_rates->{$qls[$qual_pos]}->{$base_key};
			}
			elsif ($qry[$read_pos] eq $ref[$read_pos])
			{
				my $base_key = "$ref[$read_pos]$qry[$read_pos]";
				$quality -= $log10_random_rates->{$qls[$qual_pos]}->{$base_key};
			}		
								
 			#increment positions
			$qual_pos++ if ($qry[$read_pos] ne '.');
			$read_pos++;
		}
		
		$m->{mapping_quality} = $quality;
	}	
}

return 1;
