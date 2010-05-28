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

# Format specification
# 
our $line_specification = {
	## mutational event
	'SNP' => ['position', 'ref_base', 'new_base'],
	'SUB' => [''],
	'DEL' => [''],
	'INS' => [''],
	'MOB' => [''],
	'DUP' => [''],
	'INV' => [''],
	
	## evidence
	'RA' => ['seq_id', 'position', 'insert_position', 'ref_base', 'new_base'],
	'MC' => ['seq_id', 'start', 'end'],
	'JC' => ['side_1_seq_id', 'side_1_position', 'side_1_strand', 'side_2_seq_id', 'side_2_position', 'side_2_strand', 'overlap'],
	'UN' => ['seq_id', 'start', 'end'],
};

our $tag_sort_fields = {
	'SNP' => [1, 'seq_id', 'position'],
	'INS' => [1],
	'DEL' => [1],
	'INS' => [1],
	'MOB' => [1],
	'DUP' => [1],
	'INV' => [1],
	'RA' => [2, 'seq_id', 'position'],
	'MC' => [2, 'seq_id', 'start'],
	'JC' => [2, 'side_1_seq_id', 'side_1_position'],
	'UN' => [3, 'seq_id', 'start'],
};


=head2 new

 Title   : new
 Usage   : $gd = Breseq::GenomeDiff->new();
           $gd = Breseq::GenomeDiff->new( -file => 'evolved.gd' );
 Function: Creates a GenomeDiff object, loading it from a file if one is provided
 Returns : new GenomeDiff object

=cut

sub new
{
	my($caller,@args) = @_;
	my $class = ref($caller) || $caller;
	my $self = new Bio::Root::Root($caller, @args);

	# initialize
	@{$self->{list}} = ();
	$self->{unique_id_counter} = 0;
	$self->{unique_id_used} = {};

	bless ($self, $class);
	
	# load from file if one of these options found...
	$self->{file_name} = $self->Bio::Root::RootI::_rearrange([qw(FILE_NAME)], @args);
	$self->{file_name} = $self->Bio::Root::RootI::_rearrange([qw(FILE)], @args) if (!defined $self->{file_name});
	$self->{file_name} = $self->Bio::Root::RootI::_rearrange([qw(IN)], @args) if (!defined $self->{file_name});
	$self->read($self->{file_name}) if ($self->{file_name});
	
	return $self;
}


=head2 add

 Title   : add
 Usage   : $gd->add( {seq_id=>'', start=>456, end=>546 } );
 Function: Adds a new item to a GenomeDiff from a hash. Auto-populates empty required fields.
 Returns : 

=cut
sub add
{
	my ($self, $item) = @_;
	my @missing_required_columns = ();

	## no ID, give it a new one
	if ( !defined $item->{id} )
	{
		$item->{id} = $self->new_unique_id;
	}
	elsif ($self->used_unique_id( $item->{id}) )
	{
		$self->warn("Ignoring attempt to add item with an existing id: $item->{id}");
		return;
	}
	
	##mark ID as used
	$self->{unique_id_used}->{$item->{id}} = 1;

	sub check_required_field
	{
		my ($item, $field, $missing_ref) = @_;
		push @$missing_ref, $field if (!defined $item->{$field});
	}
	
	## check to be sure the item has required fields, or auto-populate them
	$item->{type} = '' if (!defined $item->{type});
	
	my $spec = $line_specification->{$item->{type}};
	if (!defined $spec)
	{
		$self->warn("Type \'$item->{type}\' is not recognized. Ignoring item.");
		return;
	}
	
	## check for required fields
	foreach my $key (@$spec)
	{
		check_required_field($item, $key, \@missing_required_columns);
	}

	if (scalar @missing_required_columns > 0)
	{
		$self->warn("GenomeDiff::Ignoring item of type \'$item->{type}\' that is missing required field(s):" . join (',', @missing_required_columns));
		return;
	}

	## these are all required columns
	$item->{SORT_1} = $tag_sort_fields->{$item->{type}}->[0];
	$item->{SORT_2} = $item->{$tag_sort_fields->{$item->{type}}->[1]};
	$item->{SORT_3} = $item->{$tag_sort_fields->{$item->{type}}->[2]};
	
	push @{$self->{list}}, $item;
}


## 
sub used_unique_id
{
	my ($self, $id) = @_;
	return $self->{unique_id_used}->{$id};
}

##
sub new_unique_id
{
	my ($self) = @_;
	my $assigned_id = ++$self->{unique_id_counter};
	while ( $self->{unique_id_used}->{$assigned_id} )
	{
		$assigned_id = ++$self->{unique_id_counter};
	}
	return $assigned_id;
}

##
sub mark_unique_id
{
	my ($self, $id) = @_;
	$self->{unique_id_used}->{$id} = 1;
}


sub list
{
	my ($self) = @_;
	return @{$self->{list}};
}

sub list_ref
{
	my ($self) = @_;
	return $self->{list};
}



## Internal function for creating a hash from a line
sub _line_to_item
{
	my ($self, $line) = @_;
	
	my @line_list = split /\t/, $line;
	my $item = {};
	$item->{type} = shift @line_list;
	$item->{id} = shift @line_list;
	
	my $spec = $line_specification->{$item->{type}};
	if (!defined $spec)
	{
		$self->warn("Type \'$item->{type}\' is not recognized for line:\n$line");
		return undef;
	}
	
	foreach my $key (@$spec)
	{
		my $next = shift @line_list;
		if (!defined $next)
		{
			$self->warn("Number of required items is less than expected for type \'$item->{type}\' line:\n$line");
			return undef;
		}
		$item->{$key} = $next;
	}

	## Remainder of the line is comprised of non-required key value pairs
	foreach my $key_value_pair (@line_list)
	{
		next if (!$key_value_pair);
		my $matched = ($key_value_pair =~ m/^(.+)=(.+)$/);
		if (!$matched)
		{
			$self->warn("Not a key value pair \'$key_value_pair\' line:\n$line");
			next;
		}		
		
		my ($item_key, $item_value) = ($1, $2);	
		$item->{$item_key} = $item_value;
	}
	
	$item->{SORT_1} = $item->{$tag_sort_fields->{$item->{type}}->[0]};
	$item->{SORT_2} = $item->{$tag_sort_fields->{$item->{type}}->[1]};
	$item->{SORT_3} = $item->{$tag_sort_fields->{$item->{type}}->[2]};

	return $item;
}

## Internal function for creating a line from a hash
sub _item_to_line
{
	my ($self, $item) = @_;
	my $line = '';
		
	my $spec = $line_specification->{$item->{type}};
	if (!defined $spec)
	{
		$self->warn("Type \'$item->{type}\' not found for item. Ignoring.");
		return '';
	}
	
	my %ignore;
	foreach my $key ('type', 'id', @$spec)
	{
		$line .= "$item->{$key}\t";
		$ignore{$key} = 1;
	}
	
	my @keys = sort keys %$item;
	@keys = grep { !$ignore{$_} } @keys;
	@keys = grep { !($_=~ m/^SORT_/) } @keys;
	
	my @key_value_pairs = map { $_ . "=" . $item->{$_} } @keys;
	$line .= join "\t", @key_value_pairs;
	
	return $line;
}

sub read
{
	my ($self, $file_name) = @_;
	
	open IN, "<$file_name";

	#read lines, skip comment lines, and blank lines
	my @lines = <IN>;
	chomp @lines;
	@lines = grep {!/^\s*#[^=]/} @lines;
	@lines = grep {!/^\s*$/} @lines;
	
	## read version from first line
	my $l = shift @lines;
	$l =~ m/#=GENOME_DIFF\s+(\d+)/ or $self->throw("Could not match version line in file $self->{file_name}.");
	$self->{version} = $1;

	## read header information
	
	## read data
	while ($l = shift @lines)
	{
		my $item = $self->_line_to_item($l);
		push @{$self->{list}}, $item if ($item);
	}
	close IN;
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
	
	sub by_sort_fields
	{ 			
		return ($a->{SORT_1} <=> $b->{SORT_1}) || ($a->{SORT_2} cmp $b->{SORT_2}) || ($a->{SORT_3} <=> $b->{SORT_3});
	}
	

	@{$self->{list}} = sort by_sort_fields @{$self->{list}} if (!$no_sort);

	foreach my $item (@{$self->{list}})
	{		
		my $line = $self->_item_to_line($item);
		print OUT "$line\n" if ($line);
	}
	
	close OUT;
}

=head2 has_mutation

 Title   : has_mutation
 Usage   : $gd->has_mutation;
 Function: test whether a genome diff has the specified mutation
 Returns : hether a genome diff has the specified mutation

=cut

sub exists
{
	my ($self, $test_item) = @_;

	foreach my $item (@{$self->{list}})
	{
		if (	($item->{type} eq $test_item->{type})
			 &&	($item->{pos} == $test_item->{pos})
			 && ($item->{new} eq $test_item->{new})
			 && ($item->{ref} eq $test_item->{ref}) )
		{
			return 1;
		}
	}

	return 0;
}

sub filter
{
	my ($self, $filter_function) = @_;
		
	my @new_list = ();
	my $removed_gd = Breseq::GenomeDiff->new();
	
	foreach my $item ($self->list)
	{
		if ($filter_function->($item))
		{
			push @new_list, $item;
		}
		else
		{
			$removed_gd->add($item);
		}
	}
	@{$self->{list}} = @new_list;
	
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


=head2 merge

 Title   : merge
 Usage   : $gd = Bio::Breseq::GenomeDiff::merge($gd1, $gd2, $gd3, ...);
 Function: merge evidence and predictions from multiple GenomeDiffs into one
 Returns : new GenomeDiff object with items from all input GenomeDiffs

=cut

sub merge
{
	use Storable qw(dclone);
	my (@list) = @_;
	
	my $new_gd = Breseq::GenomeDiff->new();			
	while (my $gd = shift @list)
	{
		## deep copy the list, so we aren't changing the original items
		my @item_list = @{ dclone($gd->list_ref) };
		
		foreach my $item (@item_list)
		{
			## first we need to check to be sure
			## each id we are using isn't present in what we have added so far
			if ($new_gd->used_unique_id($item->{id}))
			{
				$item->{id} = $new_gd->new_unique_id();
			}
			
			$new_gd->add($item);
		}
	}
	return $new_gd;
}


sub intersection
{
	my ($list) = @_;
	
	my $union_gd = Breseq::GenomeDiff::union($list);
	
	my $new_gd = Breseq::GenomeDiff->new();
	$new_gd->{'SAMPLE'}->{intersection } = join ",", map {$_->{SAMPLE}->{strain}} @$list;
	
	SNP: foreach my $test_item (@{$union_gd->{list}})
	{
		foreach my $gd (@$list)
		{
			next SNP if (!$gd->has_item($test_item));
		}
		push @{$new_gd->{list}}, $test_item;
	}
	
	return $new_gd;
}

sub subtract
{
	my ($list_1, $list_2) = @_;
	
	my $union1_gd = Breseq::GenomeDiff::union($list_1);
	my $union2_gd = Breseq::GenomeDiff::union($list_2);

	my $new_gd = Breseq::GenomeDiff->new();
	$new_gd->{'SAMPLE'}->{subtract} = join( ",", map {$_->{SAMPLE}->{strain}} @$list_1) . "-" . join( ",", map {$_->{SAMPLE}->{strain}} @$list_2);
		
	SNP: foreach my $test_item (@{$union1_gd->{list}})
	{
		next SNP if ($union2_gd->has_mutation($test_item));
		push @{$new_gd->{list}}, $test_item;
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
