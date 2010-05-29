###
# Pod Documentation
###

=head1 NAME

MutationPredictor.pm

=head1 SYNOPSIS

Takes evidence in a GenomeDiff and adds predicted mutations

=head1 DESCRIPTION

=head1 AUTHOR

Jeffrey Barrick
<jbarrick@msu.edu>

=head1 COPYRIGHT

Copyright 2010.  All rights reserved.

=cut

###
# End Pod Documentation
###

package Breseq::MutationPredictor;
use vars qw(@ISA);
use strict;

use Bio::Root::Root;
use Data::Dumper;


@ISA = qw( Bio::Root::Root );


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
	bless ($self, $class);
	
	# initialize options
	$self->{ref_seq_info} = $self->Bio::Root::RootI::_rearrange([qw(REF_SEQ_INFO)], @args);
	$self->throw("Must provide -ref_seq_info to constructor.") if (!defined $self->{ref_seq_info});
	
	$self->{mutation_log10_e_value_cutoff} = $self->Bio::Root::RootI::_rearrange([qw(MUTATION_LOG10_E_VALUE_CUTOFF)], @args);
	$self->{mutation_log10_e_value_cutoff} = 2 if (!defined $self->{mutation_log10_e_value_cutoff});

	return $self;
}


=head2 predict

 Title   : predict
 Usage   : $mp->predict();
 Function: Predicts mutations from evidence in a GenomeDiff and adds them to it
 Returns : 

=cut
sub predict
{
	my ($self, $gd) = @_;
	
	print "here!\n";
	
	##first look at SNPs and small indels predicted by read alignments.
	my @ra = $gd->list('RA');

	##be sure they are sorted by position
	sub by_pos
	{
	       ($a->{seq_id} cmp $b->{seq_id})
		|| ($a->{position} <=> $b->{position}) 
		|| ($a->{indel_position} <=> $b->{indel_position})
	}
	@ra = sort by_pos @ra;
	
	
	## Our approach is to create a list of generic items with full properties
	## and then to delete unneeded ones afterward, when we know what kind of
	## mutation this actually was.
	
	my $mut;
	my @muts;
	
	foreach my $item (@ra)
	{
		next if ($item->{marginal});
		
		#decide whether to merge with the last mutation
		my $same = 0;
		if (defined $mut)
		{
			$same = 1 if ($mut->{end} == $item->{position}) && ($mut->{insert_end} + 1 == $item->{insert_position});
			$same = 1 if ($mut->{end} + 1 == $item->{position}) && ($item->{insert_position} == 0);
			$same = 0 if ($item->{frequency} != 1); #don't join polymorphisms
			$same = 0 if ($mut->{seq_id} ne $item->{seq_id});
		}
		
		if (!$same)
		{
			push @muts, $mut if (defined $mut);
			my $new_mut = { 
				seq_id => $item->{seq_id},
				position => $item->{position},
				start => $item->{position},
				end => $item->{position},
				insert_start => $item->{insert_position},
				insert_end => $item->{insert_position},
				ref_seq => $item->{ref_base},
				new_seq => $item->{new_base},
				evidence => [$item->{id}],
				frequency => $item->{frequency}, 
			};			
			$mut = $new_mut;
		}
		else
		{
			$mut->{insert_end} = $item->{insert_position};
			$mut->{end} = $item->{position};
			$mut->{ref_seq} .= $item->{ref_base};
			$mut->{new_seq} .= $item->{new_base};
			push @{$mut->{evidence}}, $item->{id}; 
		}
	}	
	
	foreach my $mut (@muts)
	{
		#insertion
		if ($mut->{ref_seq} =~ m/\./)
		{			
			$mut->{type} = 'INS';
		}
		#deletion
		elsif ($mut->{new_seq} =~ m/\./)
		{
			$mut->{type} = 'DEL';
			$mut->{size} = $mut->{end} - $mut->{start} + 1;
			$mut->{start_range} = 0;
			$mut->{end_range} = 0;
		}
		#block substitution
		elsif ((length $mut->{ref_seq} > 1) || (length $mut->{new_seq} > 1))
		{
			$mut->{type} = 'SUB';
		}
		#snp
		else
		{
			$mut->{type} = 'SNP';
		}
		
		$gd->add($mut);
	}	
	
	my @mc = $gd->list('MC');	
	my @jc = $gd->list('JC');
	
	## add deletions, looking for ones that are also supported by new junctions, or by ending in repeat regions
	foreach my $mc_item (@mc)
	{
		next if ($mc_item->{marginal});

		my $mut = { 
			type => 'DEL',
			seq_id => $mc_item->{seq_id},
			position => $mc_item->{start},
			size => $mc_item->{end} - $mc_item->{start} + 1,
			start_range => $mc_item->{start_range},
			end_range => $mc_item->{end_range},
			evidence => [$mc_item->{id}],
		};			
		
		## search for junctions that give the same (or similar...?) extents of deletion
		for (my $i=0; $i < scalar @jc; $i++)
		{
			my $jc_item = $jc[$i];
			
			next if ($jc_item->{side_1_seq_id} ne $mut->{seq_id});
			next if ($jc_item->{side_2_seq_id} ne $mut->{seq_id});

			if (  ($jc_item->{side_1_position} == $mut->{position}-1) && ($jc_item->{side_1_strand} == -1)
			   && ($jc_item->{side_2_position} == $mut->{position}+$mut->{size}) && ($jc_item->{side_2_strand} == +1) )
			{
				push @{$mut->{evidence}}, $jc_item->{id};
				splice @jc, $i, 1; 
				$i--;
			}
		}
		
		print Dumper($mut);	
				
		$gd->add($mut);
	}
	
	
	
	
	

	
	## first, extract SNPs and order them
	
	
	## PROBLEM: We can't apply the coverage cutoff until AFTER we count errors
	##   (because only then do we have the distribution to fit)
	##   but we have to choose which junctions we believe BEFORE counting
	##   (because we put their split alignments in the BAM file)
	## Ideally we would do this after step 7, then remove the offending read pieces from the BAM file
	## before proceeding to SNP calling.
	## this could be done by reserving these pieces in a separate SAM file
	## then merging them later? But full matches would also have to be kept separate...
	##

=comment	
	foreach my $hybrid ($hybrid_gd->list)
	{
		my $coverage_cutoff_1 = $settings->{unique_coverage}->{$hybrid->{side_1_seq_id}}->{junction_coverage_cutoff};
		my $coverage_cutoff_2 = $settings->{unique_coverage}->{$hybrid->{side_2_seq_id}}->{junction_coverage_cutoff};
		
		if ( (!defined $coverage_cutoff_1 || ($hybrid->{total_reads} < $coverage_cutoff_1) ) 
		  && (!defined $coverage_cutoff_2 || ($hybrid->{total_reads} < $coverage_cutoff_2) ) )
		{
			$hybrid->{marginal} = 1;
		}
	}
	
		
	
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
	
=cut	
}


return 1;
