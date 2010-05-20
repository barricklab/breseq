###
# Pod Documentation
###

=head1 NAME

GenomeDiff.pm

=head1 SYNOPSIS

Module for reading and writing Genome Diff (.gd) files that describe 
evidence for mutations and predicted mutations.

=head1 DESCRIPTION

=head1 FORMAT

Genome Diff files are tab delimitted. The first column defines the type of entry on a given line.
The next few columns are type-nonspecific (id, reference, start, end, parent, child) data, followed by
an arbitrary number of columns of more detailed data in a key=value format.

=head1 ENTRY TYPES

=head2 Within read evidence (code:WR)

Fixed columns: entry_type, entry_id, seq_id, start, end, ref_seq, new_seq, 

=back

=head2 Missing coverage evidence (MC)

Fixed columns: entry_type, entry_id, seq_id, start, end 

=back

=head2 New junction (NJ) evidence

Fixed columns: entry_type, entry_id, seq_id_1, position_1, direction_1, seq_id_2, position_2, direction_2, overlap


=back

=head1 AUTHOR

Jeffrey Barrick
<jbarrick@msu.edu>

=head1 COPYRIGHT

Copyright 2010.  All rights reserved.

=cut

###
# End Pod Documentation
###

package Breseq::GenomeDiff;
use vars qw(@ISA);
use strict;

use Bio::Root::Root;
use Data::Dumper;


@ISA = qw( Bio::Root::Root );

=head2 new

 Title   : new
 Usage   : $gd = Intergenic::454Diff->new( -file_name => 'evolved.gd' );
 Function: Creates a GenomeDiff object, loading it from a file if a file_name is provided
 Returns : new GenomeDiff object

=cut

sub new
{
	my($caller,@args) = @_;
	my $class = ref($caller) || $caller;
	my $self = new Bio::Root::Root($caller, @args);

	#initialize
	@{$self->{mutations}} = ();

	bless ($self, $class);
	$self->{file_name} = $self->Bio::Root::RootI::_rearrange([qw(FILE_NAME)], @args);
	$self->{file_name} = $self->Bio::Root::RootI::_rearrange([qw(FILE)], @args) if (!defined $self->{file_name});
	$self->{file_name} = $self->Bio::Root::RootI::_rearrange([qw(IN)], @args) if (!defined $self->{file_name});
	$self->read($self->{file_name}) if ($self->{file_name});
	
	return $self;
}

=head2 get_next

 Title   : get_next
 Usage   : $read = $delta_file->get_next;
 Function: get the next read from a delta file, loads all regions matched within the read
 Returns :

=cut

sub hash_out
{
	my ($self, $l, $key) = @_;
	
	my @list = split /[,|]\s+/, $l;
	my $return_hash;
	foreach my $item (@list)
	{
		next if (!$item);
		my ($item_key, $item_value) = split /=/, $item;
		die "$item_key, $item_value" if (!defined $item_value);
		$return_hash->{$item_key} = $item_value;
		$self->{$key}->{$item_key} = $item_value if (defined $key);
	}
	
	return $return_hash;
}

#convert hash to line
sub hash_to_line
{
	my ($self, $hash) = @_;
	my $line;

	foreach my $key (sort keys %$hash)
	{
		$line .= ", " if ($line);
		$line .= $key . "=" . "$hash->{$key}"
	}
	
	return $line;
}

# this should do some checking
sub add_mutation
{
	my ($self, $mut) = @_;
	push @{$self->{mutations}}, $mut;
}


sub read
{
	my ($self, $file_name) = @_;
	
	open IN, "<$file_name";

	## read version from first line
	my $l;
	$l = <IN>;
	chomp $l;
	$l =~ m/#=GENOME_DIFF\s+(\d+)/ or die "Could not match version line in file $self->{file_name}.";
	$self->{version} = $1;

	while ($l = <IN>)
	{
		chomp $l;
		if ($l =~ s/#=SAMPLE\w*//)
		{
			$self->hash_out($l, 'SAMPLE');
		}
		else
		{
			push @{$self->{'mutations'}}, $self->hash_out($l);
		}
	}
	close IN;
}

sub mutations
{
	my ($self) = @_;
	return @{$self->{'mutations'}};
}

=head2 write

 Title   : write
 Usage   : $gd->write("output.gd");
 Function: writes a genome diff file
 Returns :

=cut

sub write
{
	my ($self, $file_name, $no_sort) = @_;

	## read version from first line
	open OUT, ">$file_name" or $self->throw("Could not write file: $file_name");
	print OUT "#=GENOME_DIFF 1.0\n";
	print OUT "#=SAMPLE " . $self->hash_to_line($self->{'SAMPLE'}) . "\n" if (defined $self->{'SAMPLE'});

	sub by_seq_id_pos 
	{ 		
		my ($a_pos, $b_pos) = ($a->{pos}, $b->{pos});
		$a_pos =~ s/-\d+$//;
		$b_pos =~ s/-\d+$//;
		
		return (($a->{seq_id} cmp $b->{seq_id}) || ($a_pos <=> $b_pos));
	}

# sorting broken
#	@{$self->{mutations}} = sort by_seq_id_pos @{$self->{mutations}} if (!$no_sort);

	foreach my $mut (@{$self->{mutations}})
	{
		print OUT $self->hash_to_line($mut) . "\n";
	}
	
	close OUT;
}

=head2 has_mutation

 Title   : has_mutation
 Usage   : $gd->has_mutation;
 Function: test whether a genome diff has the specified mutation
 Returns : hether a genome diff has the specified mutation

=cut

sub has_mutation
{
	my ($self, $test_mut) = @_;

	foreach my $mut (@{$self->{mutations}})
	{
		if (	($mut->{type} eq $test_mut->{type})
			 &&	($mut->{pos} == $test_mut->{pos})
			 && ($mut->{new} eq $test_mut->{new})
			 && ($mut->{ref} eq $test_mut->{ref}) )
		{
			return 1;
		}
	}

	return 0;
}

sub filter_mutations
{
	my ($self, $filter_function) = @_;
	
	my @new_mutations = ();
	my $removed_gd = Breseq::GenomeDiff->new();
	
	foreach my $mut ($self->mutations)
	{
		if ($filter_function->($mut))
		{
			push @new_mutations, $mut;
		}
		else
		{
			$removed_gd->add_mutation($mut);
		}
	}
	@{$self->{mutations}} = @new_mutations;
	
	return $removed_gd;
}

=head2 find_common_mutations

 Title   : get_next
 Usage   : $read = $delta_file->get_next;
 Function: new genome diff object with only common mutations
 Returns :

=cut

sub last_common_ancestor
{
#	my ($list_1, $list_2) = @_;	

#	my $list_1_intersection = GenomeDiff::intersection($list_1);
#	my $list_1_intersection = GenomeDiff::intersection($list_1);


#	my $gd1 = shift @$list_1;
#	if (scalar $list_1 == 1)
#	{
#		last_common_ancestor;
#	}
#	elseif (scalar $list_1 == 1)
#	{
#		last_common_ancestor;
#	}
#	{
#	}


	
#	my @list = split /[,|]\s+/, $l;
#	my $return_hash;
#	foreach my $item (@list)
#	{
#		next if (!$item);
#		my ($item_key, $item_value) = split /=/, $item;
#		die "$item_key, $item_value" if (!$item_value);
#		$return_hash->{$item_key} = $item_value;
#		$self->{$key}->{$item_key} = $item_value if (defined $key);
#	}
	
#	return $return_hash;
}


=head2 find_common_mutations

 Title   : get_next
 Usage   : $read = $delta_file->get_next;
 Function: new genome diff object with only common mutations
 Returns :

=cut

sub union
{
	my ($list) = @_;
	
	my $new_gd = GenomeDiff->new();
	$new_gd->{'SAMPLE'}->{merged} = join ",", map {$_->{SAMPLE}->{strain}} @$list;
	
	my $snp_hash;
	foreach my $gd (@$list)
	{
		foreach my $snp (@{$gd->{SNPS}})
		{
			if (!$new_gd->has_mut($snp))
			{
				push @{$new_gd->{SNPS}}, $snp;
			}
		}
	}	
	return $new_gd;
}


sub intersection
{
	my ($list) = @_;
	
	my $union_gd = GenomeDiff::union($list);
	
	my $new_gd = GenomeDiff->new();
	$new_gd->{'SAMPLE'}->{intersection } = join ",", map {$_->{SAMPLE}->{strain}} @$list;
	
	SNP: foreach my $test_snp (@{$union_gd->{mutations}})
	{
		foreach my $gd (@$list)
		{
			next SNP if (!$gd->has_mututation($test_snp));
		}
		push @{$new_gd->{mutations}}, $test_snp;
	}
	
	return $new_gd;
}

sub subtract
{
	my ($list_1, $list_2) = @_;
	
	my $union1_gd = GenomeDiff::union($list_1);
	my $union2_gd = GenomeDiff::union($list_2);

	my $new_gd = GenomeDiff->new();
	$new_gd->{'SAMPLE'}->{subtract} = join( ",", map {$_->{SAMPLE}->{strain}} @$list_1) . "-" . join( ",", map {$_->{SAMPLE}->{strain}} @$list_2);
	
	SNP: foreach my $test_snp (@{$union1_gd->{SNPS}})
	{
		next SNP if ($union2_gd->has_mut($test_snp));
		push @{$new_gd->{SNPS}}, $test_snp;
	}
	
	return $new_gd;
}

=head2 find_subtract_mutations

 Title   : get_next
 Usage   : $read = $delta_file->get_next;
 Function: get the next read from a delta file, loads all regions matched within the read
 Returns :

=cut




return 1;
