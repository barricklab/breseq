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

package Breseq::Shared;
use strict;

use Breseq::Fastq;

require Exporter;
our @ISA = qw( Exporter );
our @EXPORT = qw( 	@int_to_base 
					%base_to_int 
					poisson
				);

use Data::Dumper;

##
# Global variables that are exported.
##

our %base_to_int = ( 'A'=>0, 'C'=>1, 'T'=>2, 'G'=>3, '.'=>4, 'N'=>5 );
our @int_to_base = ('A','C','T','G','.','N');

##
# Global variables that are NOT exported.
##	

sub system
{
	my ($command, $silent, $continue) = @_;	
	print STDERR "[system] $command\n" if (!$silent);
	my $res = CORE::system $command;
	die "$res" if ($res && !$continue);
	return $res;
}

sub poisson
{
	my ($x, $lambda) = @_;

	my $pr = exp(-$lambda);

	#divide by factorial
	for (my $i=1; $i<=$x; $i++)
	{
		$pr *= $lambda;
		$pr /= $i;
	}
	
	return $pr;
}

sub tam_next_read_alignments
{
	my ($tam, $header, $last_alignment, $paired) = @_;
	
	my $num_to_slurp = ($paired) ? 2 : 1;
	my $last_read_name;
	my $al_ref;
	if (defined $last_alignment)
	{
		$last_read_name = $last_alignment->qname;
		push @$al_ref, $last_alignment;
		$last_alignment = undef;
	}
	
	
	my $num_slurped = 0;
	ALIGNMENT: while (1)
	{
		
		$last_alignment = Bio::DB::Bam::Alignment->new();
		my $bytes = $tam->read1($header, $last_alignment);
		
		#returns bytes == -1 if EOF reached
		return ($al_ref, undef) if ($bytes < 0);
		
		my $read_name = $last_alignment->qname;
		
		if ((defined $last_read_name) && ($read_name ne $last_read_name))
		{
			$num_slurped++;
			last ALIGNMENT if ($num_slurped == $num_to_slurp);
		}
		
		$last_read_name = $last_alignment->qname if (!defined $last_read_name);
		
		push @$al_ref, $last_alignment;
	}
	
	return ($al_ref, $last_alignment);
}

sub tam_write_read_alignments
{
	my ($fh, $header, $fastq_file_index, $al, $trims) = @_;
#	print Dumper($fh, $header, $fastq_file_index, $al, $trims);
	
	for (my $i=0; $i<scalar @$al; $i++)
	{
		my $a = $al->[$i];
		my $trim = $trims->[$i];
		
		my $aux_tags = 'AS:i:' . $a->aux_get('AS') . "\t" . 'X1:i:' . (scalar @$al) . "\t" . "X2:i:$fastq_file_index";
		$aux_tags .= "\t" . "XL:i:$trim->{L}" . "\t" . "XR:i:$trim->{R}" if (defined $trim);
	
		my @score_array = $a->qscore;
		my $quality_score_string = '';
		foreach my $s (@score_array)
		{
			$quality_score_string .= chr($s+33);
		}
		
		my $cigar_list = $a->cigar_array;
		my $cigar_string = '';
		foreach my $c (@$cigar_list)
		{
			$cigar_string .= $c->[1] . $c->[0];
		}
		
		my @ll;
		push @ll, $a->qname;
		push @ll, $a->flag;
		push @ll, $header->target_name()->[$a->tid];
		push @ll, $a->start;
		push @ll, $a->qual, $cigar_string, ($a->proper_pair ? '=' : '*'), $a->mate_start, $a->isize, $a->qseq, $quality_score_string, $aux_tags;
				
		print $fh join("\t", @ll) . "\n";
	}
}

## We have to project the alignment onto a different reference sequence
## one end marked as a continuation so this can be shown in alignment.

sub tam_write_moved_alignment
{
	my $verbose = 0;
	my ($fh, $seq_id, $fastq_file_index, $a, $junction_side, $flanking_length, $junction_overlap, $junction_pos, $junction_strand, $trim) = @_;
	print STDERR Dumper($seq_id, $fastq_file_index, $a, $junction_side, $flanking_length, $junction_overlap, $junction_pos, $junction_strand, $trim) . "\n" if ($verbose);

	#alter it if we are on the opposite strand
	my $relative_strand = (($junction_side==1) ? -1 : +1) * $junction_strand;
	print STDERR "Rel strand = $relative_strand\n" if ($verbose);
	
	my $aux_tags = 'AS:i:' . $a->aux_get('AS') . "\t" . 'X1:i:' . '1' . "\t" . "X2:i:$fastq_file_index";
	$aux_tags .= "\t" . "XJ:i:1"; #this flag indicates this is a junction match (could also add which direction?)
	
	##handle putting the trims in the right places
	##need to be aware if read is trimmed out of existence??
	
	if (defined $trim)
	{
		my $trim_left = ($junction_side == 1) ? $trim->{L} : 0;
		my $trim_right = ($junction_side == 1) ? 0 : $trim->{R};
		($trim_left, $trim_right) = ($trim_right, $trim_left) if ($relative_strand == -1);
		$aux_tags .= "\t" . "XL:i:$trim_left" . "\t" . "XR:i:$trim_right";
	}

	my ($q_start, $q_end) = Breseq::Shared::alignment_query_start_end($a);

	#setup all of the original read information
	my @qual_scores = $a->qscore;
	my $cigar_list = $a->cigar_array;
	my $seq = $a->qseq;
	my $flags = $a->flag;
	
	##remove soft padding, we will add it back
	my $left_padding = 0;
	my $right_padding = 0;
	if ($cigar_list->[0]->[0] eq 'S')
	{
		$left_padding = $cigar_list->[0]->[1];
		shift @$cigar_list; 
	}
	if ($cigar_list->[-1]->[0] eq 'S')
	{
		$right_padding = $cigar_list->[-1]->[1];
		pop @$cigar_list; 
	}
	
	if ($relative_strand == -1)
	{
		$seq = Breseq::Fastq::revcom($seq);
		@qual_scores = reverse @qual_scores;
		$flags ^= 16; #bitwise XOR to flip strand
	}
		
	#need to adjust read position if insertions or deletions relative to reference
	my $seek_ref_pos = $flanking_length; 
	
	##TESTING need to capture actual place it starts for second half...	
	$seek_ref_pos -= $junction_overlap if ($junction_side == 2);
	
	my $test_ref_pos = $a->start;
	my $read_pos = $q_start;
		
	my $added = 0;
	foreach my $c (@$cigar_list)
	{
		my $n = $c->[1];
		if ($c->[0] eq 'D')
		{
			$test_ref_pos += $n;
		}
		elsif ($c->[0] eq 'I')
		{
			$read_pos += $n;
		}
		else
		{
			$read_pos += $n;
			$test_ref_pos += $n;
		}
		last if ($test_ref_pos >= $seek_ref_pos);
	}
	
	#correct for overshoot
#	my $side_1_read_pos -= ($test_ref_pos - $seek_ref_pos);	
	
	my $side_1_ref_match_length = $seek_ref_pos - $a->start + 1;
	my $side_2_ref_match_length = $a->end - $seek_ref_pos;
	my $reference_match_length = ($junction_side == 1) ? $side_1_ref_match_length : $side_2_ref_match_length;

	print STDERR "SEEK REF POS: $seek_ref_pos\n" if ($verbose);
	print STDERR "SIDE1 REF MATCH LEN: $side_1_ref_match_length\n" if ($verbose);
	print STDERR "SIDE2 REF MATCH LEN: $side_2_ref_match_length\n" if ($verbose);
	
	print STDERR "a->start: ".  $a->start . "\n" if ($verbose);
	print STDERR "a->end: ".  $a->end . "\n" if ($verbose);
	print STDERR "REF MATCH LENGTH: $reference_match_length\n" if ($verbose);

	#recall: strand == 1 means this is the lowest coordinate of that junction side sequence
	#        strand == 0 means this is the highest coordinate
	#start of match in new reference coords
	my $ref_start = ($junction_strand == 1)
		? $junction_pos
		: $junction_pos - $reference_match_length + 1;
	
	###
	## split the cigar list into two halves and keep track of their length in the read
	###
	my (@side_1_cigar_list, @side_2_cigar_list);
	my ($side_1_ref_length, $side_2_ref_length) = (0, 0);
	my ($side_1_read_length, $side_2_read_length) = (0, 0);	
	
	print STDERR Dumper($cigar_list) if ($verbose);
	while (my $c = shift @$cigar_list)
	{
		##gets longer unless this is an insertion relative to reference
		my $n = $c->[1];		
		if ($c->[0] ne 'I')
		{
			if ($side_1_ref_length + $n > $side_1_ref_match_length)
			{
				my $leftover = $side_1_ref_length + $n - $side_1_ref_match_length;
				$n = $side_1_ref_match_length - $side_1_ref_length;
				unshift @$cigar_list, [$c->[0], $leftover];
			}
			$side_1_ref_length += $n;
		}
		if ($c->[0] ne 'D')
		{
			$side_1_read_length += $n;	
		}
		push @side_1_cigar_list, [$c->[0], $n];
		
		last if ($side_1_ref_length >= $side_1_ref_match_length);
	}
	while (my $c = shift @$cigar_list)
	{
		my $n = $c->[1];
		push @side_2_cigar_list, [$c->[0], $n];
		$side_2_ref_length += $n if ($c->[0] ne 'I');
		$side_2_read_length += $n if ($c->[0] ne 'D');
	}
	
	print STDERR Dumper(\@side_1_cigar_list, $side_1_ref_length, $side_1_read_length, \@side_2_cigar_list, $side_2_ref_length, $side_2_read_length) if ($verbose);
	print STDERR "JUNCTION SIDE: $junction_side\n" if ($verbose);
	
	@$cigar_list = ($junction_side == 1) ? @side_1_cigar_list : @side_2_cigar_list;
	
	print STDERR Dumper($cigar_list) if ($verbose);
	#add padding to both sides	
	
	#original padding at the ends
	print STDERR "Left = $left_padding Right = $right_padding\n" if ($verbose);
	
	#additional padding on the end that is blocked	
	$left_padding += $side_1_read_length if ($junction_side == 2);
	$right_padding += $side_2_read_length if ($junction_side == 1);

#	($left_padding, $right_padding) = ($right_padding, $left_padding) if ($relative_strand == -1);
	
	unshift @$cigar_list, ['S', $left_padding] if ($left_padding);
	push @$cigar_list, ['S', $right_padding] if ($right_padding);

	@$cigar_list = reverse @$cigar_list if ($relative_strand == -1);

	print STDERR Dumper($cigar_list) if ($verbose);
	print STDERR "Ref match length: $reference_match_length\n" if ($verbose);
	
	my $cigar_string = '';
	my $cigar_length = 0;
	foreach my $c (@$cigar_list)
	{
		$cigar_string .= $c->[1] . $c->[0];
		$cigar_length += $c->[1] if ($c->[0] ne 'D');
	}
	
			
	## assemble the quality score string
	my $quality_score_string = '';
	foreach my $s (@qual_scores)
	{
		$quality_score_string .= chr($s+33);
	}
	
	my @ll;
	push @ll, $a->qname;
	push @ll, $flags;
	push @ll, $seq_id;
	push @ll, $ref_start;
	push @ll, $a->qual, $cigar_string, ($a->proper_pair ? '=' : '*'), $a->mate_start, $a->isize, $seq, $quality_score_string, $aux_tags;	
	my $l = join("\t", @ll) . "\n";
	print STDERR $l if ($verbose);

	die "CIGAR length does not match calculated length" if ($cigar_length != $a->l_qseq);

	print $fh $l;	
}


sub alignment_query_start_end
{
	my ($a, $options) = @_;
	my $ca = $a->cigar_array;
	
	## The cigar array will be undefined if there was no match for this sequence
	return (0,0) if (!defined $ca || (scalar @$ca < 1));
	
	my $start = 1;
	$start += $ca->[0]->[1] if ($ca->[0]->[0] eq 'S');
	my $end = $a->query->length;
	$end -= $ca->[-1]->[1] if ($ca->[-1]->[0] eq 'S');
	
	if ($a->reversed && !$options->{no_reverse})
	{
		($start, $end) = ($a->query->length - $start + 1, $a->query->length - $end + 1);
		($start, $end) = ($end, $start);
	}
	
	return ($start, $end);
}

our $junction_name_separator = '__';
sub junction_name_join
{
	$_[6] =~ tr/-/n/;	
	return join "$junction_name_separator", @_;
}
sub junction_name_split
{
	my @s = split /$junction_name_separator/, $_[0];
	$s[6] =~ tr/n/-/;
	return @s;
}

## moved here from ReferenceSequence.pm to avoid BioPerl requirement.
sub process_reference_sequences
{
	my ($settings, $summary) = @_;
	$summary->{sequence_conversion}->{total_reference_sequence_length} = 0;
	my $s;
	my $ref_seq_info;
	my $i = 0;
	
	# get two pieces of information from $settings
	my $reference_fasta_file_name = $settings->file_name('reference_fasta_file_name');
	my @genbank_file_names = $settings->file_name('reference_genbank_file_names'); 
		
	print STDERR "Loading reference sequences...\n";

	##open fasta file
	open FASTA, ">$reference_fasta_file_name";

	my %loaded_seq_ids;

	foreach my $genbank_file_name (@genbank_file_names)
	{
		## open this GenBank file		
		open GENBANK, "<$genbank_file_name";
		
		my $seq_id;
		my $ref_seq = '';
		my $specified_length = -1;
		my $actual_length = 0;
		my $seq_definition;
		my $seq_version;
		my $started_on_sequence = 0;
		while ($_ = <GENBANK>)
		{			
			chomp $_;
			
			##end of a record
			if ($_ =~ m/^\s*\/\/\s*$/)
			{				
				## error checking
				($actual_length == $specified_length) or die "Error reading GenBank file entry: $seq_id\nLength in header ($specified_length) does not match length of sequence ($actual_length).\n";  
				(!$loaded_seq_ids{$seq_id}) or die "Duplicate GenBank file entry: $seq_id\n";  
				$loaded_seq_ids{$seq_id}++;
				
				$s->{$seq_id}->{length} = $actual_length;
				$s->{$seq_id}->{definition} = (defined $seq_definition) ? $seq_definition : '';
				
				$s->{$seq_id}->{seq_id} = $seq_id;
				$s->{$seq_id}->{version} = $seq_version;
				
				$s->{$seq_id}->{string} = $seq_id;
				
				
				$summary->{sequence_conversion}->{total_reference_sequence_length} += $actual_length;
				
				## it would be nice to get rid of storing the whole genome in memory
				$ref_seq_info->{ref_strings}->{$seq_id} = $ref_seq;
				$ref_seq_info->{seq_order}->{$seq_id} = $i++;
				push @{$ref_seq_info->{seq_ids}}, $seq_id;

				### re-initialize to nothing for next sequence
				undef $seq_id;
				undef $seq_definition;
				undef $seq_version;
				$ref_seq = '';
				$specified_length = -1;
				$started_on_sequence = 0;
				
			}
			elsif ($_ =~ m/\s*LOCUS\s+(\S+)\s+(\d+)\s+bp/)
			{				
				$seq_id = $1; 
				$specified_length = $2;
			}
			elsif ($_ =~ m/\s*DEFINITION\s+(.+)$/)
			{
				$seq_definition = $1; 
			}
			elsif ($_ =~ m/\s*VERSION\s+(.+)$/)
			{
				$seq_version = $1; 
			}			
			elsif ($_ =~ m/\s*ORIGIN/)
			{
				print FASTA ">$seq_id\n";
				$started_on_sequence = 1; 
			}		
			elsif ($started_on_sequence)
			{
				$_ =~ s/(\d|\s)//g;
				$actual_length += length $_;
				print FASTA "\U$_\n"; 
				$ref_seq .= "\U$_";
			}			
		}
		
		print STDERR "  Loading File::$genbank_file_name\n";
	}
	
	## create SAM faidx
	Breseq::Shared::system("samtools faidx $reference_fasta_file_name", 1);
	
	$summary->{sequence_conversion}->{reference_sequences} = $s;	
	
	return $ref_seq_info;
}



return 1;

