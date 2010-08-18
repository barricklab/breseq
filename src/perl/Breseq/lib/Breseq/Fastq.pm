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
use Bio::Root::Root;

package Breseq::Fastq;
use vars qw(@ISA);
@ISA = qw( Bio::Root::RootI );

use Data::Dumper;

###
# Global Variables
###

=head2 new

 Title   : new
 Usage   : $diff = Intergenic::454Diff->new( -file_name => 'sequence.HCDiff.txt' );
 Function: Creates a diff object and initializes it from a file.
 Returns : a new 454Diff object

=cut


our $format_to_chr_offset = {
	'sanger' => 33,
	'solexa' => 64,
	'illumina' => 64,
};

our $format_to_quality_type = {
	'sanger' => 'phred',
	'solexa' => 'solexa',
	'illumina' => 'phred', 
};


sub new
{
	my($caller,@args) = @_;
	my $class = ref($caller) || $caller;
	my $self = new Bio::Root::Root($caller, @args);

	bless ($self, $class);
	($self->{file_name}, $self->{format}) = $self->Bio::Root::RootI::_rearrange([qw(FILE FORMAT)], @args);
	$self->{file_name} or $self->throw("Must provide -file parameter");
	#$self->{dont_reverse} = $self->Bio::Root::RootI::_rearrange([qw(dont_reverse)], @args);	
	
	## assume sanger
	(defined $self->{format}) or $self->{format} = 'sanger';
	$self->{chr_offset} = $format_to_chr_offset->{$self->{format}};
	$self->{quality_type} = $format_to_quality_type->{$self->{format}};
	die "invalid format" if (!defined $self->{chr_offset});

	$self->{list_format} = $self->Bio::Root::RootI::_rearrange([qw(LIST_FORMAT)], @args);
	
	#opening for writing
	open $self->{fh}, "$self->{file_name}" or die "Could not open file $self->{file_name}";
	if ($self->{file_name} =~ m/^>/)
	{
	}
	#opening for reading
	else
	{	
		open $self->{fh}, "$self->{file_name}" or die "Could not open file $self->{file_name}";
		$self->{next_line} = readline $self->{fh};
		chomp $self->{next_line} if ($self->{next_line});
	}
	
	return $self;
}

=head2 next_seq

 Title   : next_seq
 Usage   : $read = $fastq_lite->read_seq;
 Function: get the next fastq sequence
 Returns :

=cut

sub next_seq
{
	my ($self) = @_;
	return undef if (!defined $self->{next_line});
		
	my $new_seq;

	#print STDERR "\@LINE:::  $self->{next_line}\n";

	#first line should start with '@'
	if ($self->{next_line} =~ s/^@//)
	{
		$new_seq->{id} = $self->{next_line};
	}
	else
	{
		$self->throw("Sequence record does not begin with \'@\' as expected.");	
	}
	#strip id after spaces!!
	$new_seq->{id} =~ s/\s.*$//;
	
	#next line is "sequence"
	$self->{next_line} = readline $self->{fh};
	chomp $self->{next_line};
	$new_seq->{seq} = "\U$self->{next_line}";
	$new_seq->{seq} =~ s/\r//g; #clean sequence
		
	#next line is "+"
	$self->{next_line} = readline $self->{fh};
	
	#next line is quals
	$self->{next_line} = readline $self->{fh};
	chomp $self->{next_line};
	$new_seq->{qual_chars} = $self->{next_line};
	$new_seq->{qual_chars} =~ s/\r//g; #clean sequence
	
	#next line should be next sequence
	$self->{next_line} = readline $self->{fh};
	chomp $self->{next_line} if (defined $self->{next_line});
	
	### if we are using list format we need to convert quals to ascii
	if ($self->{list_format})
	{
		my @qual_list = split /\s+/, $new_seq->{qual_chars};
		$new_seq->{qual_chars} = '';
				
		foreach my $q (@qual_list)
		{
			die "chr out of range $q, fastq type may be incorrect" if ($q+$self->{chr_offset} < 0);
			die "chr out of range $q, fastq type may be incorrect" if ($q+$self->{chr_offset} > 255);
			$new_seq->{qual_chars}.= chr($q+$self->{chr_offset});
		}
	}
	
	### convert to sanger
	if ($self->{format} ne 'sanger') 
	{
		my $new_quals = '';
		foreach my $q (split //, $new_seq->{qual_chars})
		{
			##convert offset
			$q = ord($q);
			$q -= $self->{chr_offset};
			
			if ($self->{quality_type} ne 'phred')
			{
				$q = convert_quality_solexa_to_phred($q);
			}
			$q += 33;
			$new_quals .= chr($q);
		}
		$new_seq->{qual_chars} = $new_quals;
	}
	
	return $new_seq;
}


our %solexa_to_phred_conversion_table;
our $init_solexa_to_phred_conversion_table = 0;

sub populate_solexa_to_phred_conversion_table
{
	my ($nothing) = @_;
	
	for (my $sq=-5; $sq<62; $sq++)
	{
		my $p = 10**(-$sq/10) / (1+10**(-$sq/10));
		my $pq = -10 * log($p) / log(10);
		$solexa_to_phred_conversion_table{$sq} = int($pq);
	} 
	$init_solexa_to_phred_conversion_table = 1;
#	print STDERR Dumper(\%solexa_to_phred_conversion_table);
}

sub convert_quality_solexa_to_phred
{
	my ($sq) = @_;	
	populate_solexa_to_phred_conversion_table() if (!$init_solexa_to_phred_conversion_table);
	return $solexa_to_phred_conversion_table{$sq};
}

=head2 predict_qual_format

 Title   : predict_qual_format
 Usage   : $quality_type = $fastq_lite->predict_qual_format($n);
 Function: predict quality type from first $n sequences
 Returns :

=cut

sub predict_format
{
	my ($self, $n, $verbose) = @_;

	my %qual_found_hash;
	my $average_qual;
	my $max_qual;
	my $min_qual;
	my $qual_num = 0;
	for (my $i=0; $i<$n; $i++)
	{
		my $seq = $self->next_seq;
		last if (!$seq);
		
		my @quals = Breseq::Fastq::quals($seq);
		foreach my $q (@quals)
		{
			$qual_found_hash{$q}++;
			$average_qual+=$q;
			$qual_num++;
			
			$max_qual = $q if ((!defined $max_qual) || ($q > $max_qual));
			$min_qual = $q if ((!defined $min_qual) || ($q < $min_qual));			
		}
	}
	$average_qual /= $qual_num;
	my @sorted_quals = sort {-($qual_found_hash{$a} <=> $qual_found_hash{$b})} keys (%qual_found_hash);
	my $most_common_qual = $sorted_quals[0];	
	
	## what is the most common quality score?
	## RECALL: these have already been treated as if they were phred style with 33 offset ***
	## So if we believe it is Solexa, we must return values that are further offset;
	
	print STDERR "  Average quality: $average_qual\n" if ($verbose);
	print STDERR "  Most common quality: $most_common_qual\n" if ($verbose);
	print STDERR "  Quality Score Distribution:\n"  if ($verbose);
	
	## record the cumulative distribution of quality scores
	my $cdf = 0;
	my @qual_cdf;
	for (my $i=$min_qual; $i<=$max_qual; $i++)
	{
		my $pr = (defined $qual_found_hash{$i}) ? $qual_found_hash{$i} : 0;
		$pr /= $qual_num;
		$cdf += $pr;
		$qual_cdf[$i-$min_qual] = $cdf;				
		print STDERR "    $i: pr=" . sprintf("%.2f", $pr * 100) . "% cdf=" . sprintf("%.2f", $cdf * 100) . "%\n"  if ($verbose);
	}
	
	##choose fairly conservative cutoffs
	return ('sanger', $min_qual, $max_qual, \@qual_cdf)  if ($average_qual < 50);
	return ('illumina', $min_qual-31, $max_qual-31, \@qual_cdf) if ($average_qual > 50);

	#print STDERR "  Failure to auto-detect quality score format.\n  Specify one manually!\n";
	#print Dumper(\%qual_found_hash);
	#die "\n";
	
	return ('', 0, 0, undef);
}


sub write_seq
{
	my ($self, $seq) = @_;
	print {$self->{fh}} "\@$seq->{id}\n$seq->{seq}\n+\n$seq->{qual_chars}\n";
}

sub complement
{
	my ($seq) = @_;
	$seq =~ tr/ATCG/TAGC/;
	return $seq;
}

sub revcom
{
	my ($seq) = @_;
	$seq =~ tr/ATCG/TAGC/;
	return reverse $seq;
}

##only do conversion if explicitly asked to save time
##correctly removes offsets for phred(+33) vs solexa (+64)
##solexa scores will be adjusted at a later stage...
sub quals
{
	my ($seq) = @_;
	return @{$seq->{qual_list}} if (defined $seq->{qual_list});
	@{$seq->{qual_list}} = split //, $seq->{qual_chars};
	@{$seq->{qual_list}} = map { unpack("C",$_) - $seq->{chr_offset} } (@{$seq->{qual_list}}); 
	return @{$seq->{qual_list}};
}

##
## stripped down utility routine to split fastq file into smaller chunks for mummer
##
sub split_fastq_to_fasta
{
	my ($input, $output, $options) = @_;
	die if (!defined $options->{number_per_file});

	open INPUT, "<$input";

	my $line = <INPUT>;

	my $total_bases = 0;
	my $max_read_length = 0;
	my $num_seq = 0;
	my $num_per = 2;
	my $remainder = 0;

	my $sequence_index = 0;
	my $file_index = 0;

	do {
		my $output_name = $output;
		if ($options->{number_per_file})
		{
			if ($output_name =~ m/=/)
			{
				$output_name =~ s/=/$file_index/;
			}
			else
			{
				$output_name .= "_$file_index";
			}
		}
		open OUTPUT, ">$output_name" or die "Could not open $output_name\n";
		
		while ($line)
		{
			$line =~ s/^\@/>/;
			if ($options->{add_read_index})
			{
				chomp $line;
				$line =~ s/>//;
				$line = ">" . $sequence_index . ":" . $line . "\n";
			}
			elsif ($options->{strip_names})
			{
				$line =~ s/[^\/]+/$num_seq/;
				$line = ">" . $line;
			}
			
			print OUTPUT $line;
			$line = <INPUT>;
			print OUTPUT $line;
			
			#chomp line ending after printing, but before recording length!
			chomp $line;
			$total_bases += length($line);
			$max_read_length = length($line) if (length($line) > $max_read_length);
			
			$line = <INPUT>;
			$line = <INPUT>;
			$line = <INPUT>;

			if ($options->{paired})
			{
				$remainder++;
				if ($remainder % $num_per == 0)
				{
					$num_seq++;
					$remainder = 0;
				}
			}
			else
			{
				$num_seq++;
			}
			$sequence_index++;
			last if (defined ($options->{number_per_file}) && ($sequence_index % $options->{number_per_file} == 0));
		}
		
		#if we get here, then loop exited because there are no more lines
		$file_index++;		
	} while ($line);
	
	##return the number of reads
	return ($sequence_index, $total_bases, $max_read_length);
}


##
## shorten the names of fastq sequences and give them a file index
##
sub fastq_to_encoded_name_fastq
{
	my ($input, $output, $fastq_file_index, $trim) = @_;

	open IN, "<$input" or die "Could not open input file $output";
	open OUT, ">$output" or die "Could not open output file $output";

	my $line = <IN>;
	my $read_index = 0;
	
	while ($line)
	{
		$read_index++;
		($line =~ m/^\@/) or die; #small sanity check
		
		#new name
		
		#then copy all of the other three lines normally
		my $sequence = <IN>;
		my $second_accession = <IN>;
		my $quality = <IN>;
		
		#read must still be 25 nt long
		print OUT "\@$fastq_file_index:$read_index\n";
		print OUT $sequence;
		print OUT $second_accession;
		print OUT $quality;
		
		$line = <IN>;
	}
}

sub base_quality_filter_fastq
{
	my $verbose = 0;
	my ($input, $output, $base_quality_cutoff) = @_;
	my $minimum_read_Length = 25;

	my $total_reads_processed = 0;
	my $total_reads_removed = 0;
	my $total_bases_kept = 0;
	my $total_bases_removed = 0;
	
	open IN, "<$input" or die "Could not open input file $output";
	open OUT, ">$output" or die "Could not open output file $output";

	my $name = <IN>;
	chomp $name;
	my $read_index = 0;

	my $chr_offset = 33; #breseq
	my $base_chr_cutoff = $chr_offset + $base_quality_cutoff;

	while ($name)
	{
		$read_index++;
		$total_reads_processed++;
		($name =~ m/^\@/) or die "Read FASTQ must have sequence on a single line"; #small sanity check

		#then copy all of the other three lines normally
		my $sequence = <IN>;
		chomp $sequence;
		my $second_accession = <IN>;
		chomp $second_accession;
		my $quality = <IN>;
		chomp $quality;
		my @quals = split //, $quality;
		my $pos = $#quals;

		my $original_read_length = length $sequence;

		my $last_high_quality_base_pos_encountered = 0;
		foreach (my $pos = 0; $pos < $original_read_length; $pos++)
		{
#			print "$pos " . ord($quals[$pos]) . "\n" if ($verbose);

			#replace low scoring bases with N's
			if (ord($quals[$pos]) < $base_chr_cutoff)
			{
				substr $sequence, $pos, 1, 'N';
			}
			else
			{
				$last_high_quality_base_pos_encountered = $pos;
			}
		}

		my $high_quality_length = $last_high_quality_base_pos_encountered + 1;

		print "$high_quality_length / $original_read_length\n" if ($verbose);

		#read must still be 25 nt long
		if ($high_quality_length >= $minimum_read_Length)
		{
			$total_bases_removed += $original_read_length - $high_quality_length;
			$total_bases_kept += $high_quality_length;

			print OUT "$name\n";
			print OUT substr($sequence, 0, $high_quality_length) . "\n";
			print OUT "+\n"; ##always truncate this line to save space
			print OUT substr($quality, 0, $high_quality_length) . "\n";

		}
		else
		{
			$total_reads_removed++;
			$total_bases_removed += $original_read_length;
		}


		$name = <IN>;
		chomp $name if (defined $name);
		### should keep track of some statistics about what was lost!
	}

	print "Total reads processed: $total_reads_processed\n";
	print "Total reads removed: $total_reads_removed\n";
	print "Total bases kept: $total_bases_kept\n";
	print "Total bases removed: $total_bases_removed\n";

	return ($total_reads_processed, $total_reads_removed, $total_bases_kept, $total_bases_removed);


}


##
## quality trim fastq
##
sub fastq_to_trimmed_fastq
{
	my $verbose = 0;
	my ($input, $output) = @_;

	my $total_reads_processed = 0;
	my $total_reads_removed = 0;
	my $total_bases_kept = 0;
	my $total_bases_removed = 0;

	open IN, "<$input" or die "Could not open input file $output";
	open OUT, ">$output" or die "Could not open output file $output";

	my $name = <IN>;
	chomp $name;
	my $read_index = 0;

	my $chr_offset = 33; #breseq
	my $quality_score_cutoff = $chr_offset + 10;
	
	while ($name)
	{
		$read_index++;
		$total_reads_processed++;
		($name =~ m/^\@/) or die; #small sanity check
				
		#then copy all of the other three lines normally
		my $sequence = <IN>;
		chomp $sequence;
		my $second_accession = <IN>;
		chomp $second_accession;
		my $quality = <IN>;
		chomp $quality;
		my @quals = split //, $quality;
		my $pos = $#quals;
		
		my $original_read_length = length $sequence;
		
		my $number_low_quality_bases_encountered = 0;
		my $number_low_quality_bases_limit = 3;
		my $last_high_quality_base_pos_encountered = 0;
		foreach (my $pos = 0; $pos < $original_read_length; $pos++)
		{
#			print "$pos " . ord($quals[$pos]) . "\n" if ($verbose);
			
			if (ord($quals[$pos]) < $quality_score_cutoff)
			{
				$number_low_quality_bases_encountered++;
			}
			else
			{
				$last_high_quality_base_pos_encountered = $pos;
			}
			
			last if ($number_low_quality_bases_encountered > $number_low_quality_bases_limit);
		}
		
		my $high_quality_length = $last_high_quality_base_pos_encountered + 1;
		my $minimum_read_Length = 25;
		
		print "$high_quality_length / $original_read_length\n" if ($verbose);
		
		#read must still be 25 nt long
		if ( ($number_low_quality_bases_encountered <= $number_low_quality_bases_limit) && ($high_quality_length >= $minimum_read_Length) )
		{
			$total_bases_removed += $original_read_length - $high_quality_length;
			$total_bases_kept += $high_quality_length;
			
			print OUT "$name\n";
			print OUT substr($sequence, 0, $high_quality_length) . "\n";
			print OUT "+\n"; ##always truncate this line to save space
			print OUT substr($quality, 0, $high_quality_length) . "\n";
			
		}
		else
		{
			$total_reads_removed++;
			$total_bases_removed += $original_read_length;
		}
		
		
		$name = <IN>;
		chomp $name if (defined $name);
		### should keep track of some statistics about what was lost!
	}
	
	print "Total reads processed: $total_reads_processed\n";
	print "Total reads removed: $total_reads_removed\n";
	print "Total bases kept: $total_bases_kept\n";
	print "Total bases removed: $total_bases_removed\n";
	
	return ($total_reads_processed, $total_reads_removed, $total_bases_kept, $total_bases_removed);
}

return 1;
