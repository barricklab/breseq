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
	our ($self, $settings, $summary, $ref_seq_info, $gd) = @_;
		
	##first look at SNPs and small indels predicted by read alignments.
	my @ra = $gd->list('RA');

	##be sure they are sorted by position
	sub by_pos
	{		
	       ($a->{seq_id} cmp $b->{seq_id})
		|| ($a->{position} <=> $b->{position}) 
		|| ($a->{insert_position} <=> $b->{insert_position})
	}
	@ra = sort by_pos @ra;
	
	
	## Our approach is to create a list of generic items with full properties
	## and then to delete unneeded ones afterward, when we know what kind of
	## mutation this actually was.

	###
	## Also, gather together read alignment mutations that occur next to each other
	## ...unless they are polymorphisms
	###
	
	my $mut;
	my @muts;
	
	foreach my $item (@ra)
	{
		next if ($item->{reject});
		
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
	##don't forget the last one
	push @muts, $mut if (defined $mut);
	
	
	
	##add information about repeat_regions overlapping the sides of junctions	
	my @jc = $gd->list('JC');
	JC: foreach my $j (@jc)
	{
		$j->{_side_1_read_side} = -1;
		$j->{_side_2_read_side} = +1;					

		foreach my $side_key ('side_1', 'side_2')
		{						
			$j->{"_$side_key\_is"} = Breseq::ReferenceSequence::find_closest_repeat_region(
				$j->{"$side_key\_position"}, 
				$ref_seq_info->{repeat_lists}->{$j->{"$side_key\_seq_id"}}, 
				50, 
				$j->{"$side_key\_strand"}
			);
			
			$j->{"$side_key\_annotate_key"} = ((defined $j->{"_$side_key\_is"}) || ($j->{"$side_key\_redundant"})) ? 'repeat' : 'gene';				
		}
		
		## Determine which side of the junction is the IS and which is unique
		## these point to the correct initial interval...
		if (defined $j->{_side_1_is})
		{	
			if (abs($j->{_side_1_is}->{start} - $j->{side_1_position}) <= 20)
			{
				$j->{_is_interval} = 'side_1';
				$j->{_is_interval_closest_side_key} = 'start';
				$j->{_unique_interval} = 'side_2';	
			}
			elsif (abs($j->{_side_1_is}->{end} - $j->{side_1_position}) <= 20 )
			{
				$j->{_is_interval} = 'side_1';
				$j->{_is_interval_closest_side_key} = 'end';
				$j->{_unique_interval} = 'side_2';		
			}			
		}
		elsif (defined $j->{_side_2_is})
		{
			if (abs($j->{_side_2_is}->{start} - $j->{side_2_position}) <= 20)
			{
				$j->{_is_interval} = 'side_2';
				$j->{_is_interval_closest_side_key} = 'start';
				$j->{_unique_interval} = 'side_1';
			}
			elsif (abs($j->{_side_2_is}->{end} - $j->{side_2_position}) <= 20 )
			{
				$j->{_is_interval} = 'side_2';
				$j->{_is_interval_closest_side_key} = 'end';
				$j->{_unique_interval} = 'side_1';
			}
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
	
	## add deletions
	## look for ones that are also supported by new junctions, and clean up those that end
	## in repeat regions of the same kind.
	##
	MC: foreach my $mc_item (@mc)
	{
		next if ($mc_item->{reject});

		my $mut = { 
			type => 'DEL',
			seq_id => $mc_item->{seq_id},
			position => $mc_item->{start},
			size => $mc_item->{end} - $mc_item->{start} + 1,
			evidence => [$mc_item->{id}],
		};			


		## Search for scenarios where>
		##
		## (1) junctions that give the exact same extents of deletion
		## (2) there is a junction between unique sequence and a repeat element
		## (3) there is is no junction, but both ends of the deletion are in repeat sequences

		
		JUNCTION: for (my $i=0; $i < scalar @jc; $i++)
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
				$gd->add($mut);
				next MC;
			}
		}
		
		sub within_repeat
		{
			my ($seq_id, $position) = @_;
			foreach my $r (@{$self->{ref_seq_info}->{repeat_lists}->{$seq_id}})
			{
				return $r if ($r->{start} <= $position) && ($position <= $r->{end})
			}
			return undef;
		}
		
		## Are we within two copies of the same repeat region??
		my $r1 = within_repeat($mut->{seq_id}, $mut->{position}); 
		my $r2 = within_repeat($mut->{seq_id}, $mut->{position} + $mut->{size}); 
		
		## Then we will adjust the coordinates to remove...
		if (defined $r1 && defined $r2 && ($r1->{name} eq $r2->{name}))
		{
			#there may be more evidence that one or the other is deleted...
			my $r1_overlap_end = $mc_item->{start} + $mc_item->{start_range};
			$r1_overlap_end = $r1->{end} if ($r1_overlap_end > $r1->{end});
			my $r1_overlap = $r1_overlap_end - $mc_item->{start} + 1;
			
			my $r2_overlap_start = $mc_item->{end} - $mc_item->{end_range};
			$r2_overlap_start = $r2->{start} if ($r2_overlap_start < $r1->{start});
			my $r2_overlap = $mc_item->{end} - $r2_overlap_start + 1;				
			
			# it may be really close...defined by read length of genome in which case
			my $slop_distance = $summary->{sequence_conversion}->{max_read_length};
			
			## prefer to delete the second copy
			if ((abs($r1_overlap - $r2_overlap) <= $slop_distance) || ($r2_overlap > $r1_overlap ))
			{
				$mut->{position} = $r1->{end} + 1;
				$mut->{size} = $r2->{end} - $r1->{end};
			}
			else #delete the first copy
			{
				$mut->{position} = $r1->{start};
				$mut->{size} = $r2->{start} - $r1->{start};
			}				
			
			## remember the name of the element
			$mut->{between} = $r1->{name};
			$gd->add($mut);	
			next MC;			
		}
		
		## Both sides were unique or redundant, nothing more we can do...
		next MC if (!defined $r1 && !defined $r2);
		next MC if (defined $r1 && defined $r2);
		
		## One of the two sides was defined as a repeat
		my $r = (defined $r1) ? $r1 : $r2; 
		my $redundant_deletion_side = (defined $r1) ? -1 : +1; 
		my $needed_coord = (defined $r1) ?  $mut->{position}+$mut->{size} : $mut->{position} - 1;
				
				
		my $verbose = 0;		
		print Dumper($mut) if ($verbose);
		print Dumper($r) if ($verbose);
		

		JUNCTION: for (my $i=0; $i < scalar @jc; $i++)
		{
			my $j = $jc[$i];


			next JUNCTION if (!defined $j->{_is_interval});
			
			print Dumper($j) if ($verbose);
			
			print "Check 1: " . $j->{"$j->{_unique_interval}_seq_id"} . " ne $mut->{seq_id}\n" if ($verbose);
			next JUNCTION if ($j->{"$j->{_unique_interval}_seq_id"} ne $mut->{seq_id});
			print "Pass 1\n" if ($verbose);

			#check type of IS
			print "Check 2: " . $r->{name} . " ne " .  $j->{"_$j->{_is_interval}_is"}->{name} . "\n" if ($verbose);			
			next JUNCTION if ( $r->{name} ne $j->{"_$j->{_is_interval}_is"}->{name} );
			print "Pass 2\n" if ($verbose);
			
			#check that IS is on the right strand
			print "Check 3: " . $redundant_deletion_side . " * " . $r->{strand} . " != " .  $j->{"$j->{_unique_interval}\_strand"}  . " * " . $j->{"_$j->{_is_interval}_is"}->{strand} . " * " . $j->{"_$j->{_is_interval}_read_side"}  . "\n" if ($verbose);							
			next JUNCTION if ( $redundant_deletion_side * $r->{strand} !=  $j->{"$j->{_unique_interval}\_strand"} * $j->{"_$j->{_is_interval}_is"}->{strand} * $j->{"_$j->{_is_interval}_read_side"} );
			print "Pass 3\n" if ($verbose);
			
			#check that the unique side matches coordinate
			print "Check 4: " . $j->{"$j->{_unique_interval}\_position"} . " != " .  $needed_coord . "\n" if ($verbose);			
			next JUNCTION if ( $j->{"$j->{_unique_interval}\_position"} != $needed_coord );
			print "Pass 4\n" if ($verbose);

			#check that the unique side is on the right strand	
			print "Check 5: " . $redundant_deletion_side . " != " .  $j->{"$j->{_unique_interval}\_strand"} . " * " . $j->{"_$j->{_is_interval}_read_side"} . "\n" if ($verbose);				
			next JUNCTION if ( $redundant_deletion_side != $j->{"$j->{_unique_interval}\_strand"} * $j->{"_$j->{_is_interval}_read_side"} );
			print "Pass 5\n" if ($verbose);

			## need to adjust the non-unique coords
			if ($redundant_deletion_side == -1)
			{
				my $move_dist = $r->{end} + 1 - $mut->{position};
				$mut->{position} += $move_dist;
				$mut->{size} -= $move_dist;
			}
			else
			{
				my $move_dist = ($mut->{position} + $mut->{size} - 1) - ($r->{start}-1);
				$mut->{size} -= $move_dist;
			}

			## OK, we're good!
			push @{$mut->{evidence}}, $j->{id};
			splice @jc, $i, 1; 
			$i--;
			$gd->add($mut);
			next MC;
		}
						
	}
	
	##infer IS element insertions from pairs of new junctions	
	JC: foreach my $j (@jc)
	{					
		## Ah, we don't have an IS, we are done
		next JC if (!defined $j->{_is_interval});
		
		## Ah, there is no overlap to play with, we are done
		next JC if ($j->{overlap} <= 0);
		
		## The following code implies $j->{overlap} > 0
				
		### first, adjust the repetitive sequence boundary to get as close to the IS as possible
		my $move_dist = abs($j->{"$j->{_is_interval}\_position"} - $j->{"_$j->{_is_interval}\_is"}->{$j->{_is_interval_closest_side_key}});
		$move_dist = $j->{overlap} if ($move_dist > $j->{overlap});
		$j->{"$j->{_is_interval}\_position"} += $j->{"$j->{_is_interval}\_strand"} * $move_dist;
		$j->{overlap} -= $move_dist;
		
		### second, adjust the unique sequence side with any remaining overlap
		$j->{"$j->{_unique_interval}\_position"} += $j->{"$j->{_unique_interval}\_strand"} * $j->{overlap};						
					
		$j->{overlap} = 0;
		
	}
	

	sub by_hybrid
	{
		my $a_pos = (defined $a->{_side_1_is}) ? $a->{_side_2}->{position} : $a->{_side_1}->{position};
		my $b_pos = (defined $b->{_side_1_is}) ? $b->{_side_2}->{position} : $b->{_side_1}->{position};

		my $a_seq_order = (defined $a->{_side_1_is}) ? $ref_seq_info->{seq_order}->{$a->{_side_2}->{seq_id}} : $ref_seq_info->{seq_order}->{$a->{_side_1}->{seq_id}};
		my $b_seq_order = (defined $b->{_side_1_is}) ? $ref_seq_info->{seq_order}->{$b->{_side_2}->{seq_id}} : $ref_seq_info->{seq_order}->{$b->{_side_1}->{seq_id}};		

		return (($a_seq_order <=> $b_seq_order) || ($a_pos <=> $b_pos));
	}
	@jc = sort by_hybrid @jc;
	
	
	foreach (my $i=0; $i<scalar(@jc)-1; $i++)
	{	
		my $j1 = $jc[$i];
		my $j2 = $jc[$i+1];		

		#must be same IS
		next if (!defined $j1->{_is_interval} || !defined $j2->{_is_interval});
		next if ($j1->{"_$j1->{_is_interval}\_is"}->{name} ne $j2->{"_$j2->{_is_interval}\_is"}->{name});

		## positive overlap should be resolved by now
		die if ($j1->{overlap} > 0);
		die if ($j2->{overlap} > 0);

		#must be close together in real coords
		next if (abs($j1->{"$j1->{_unique_interval}\_position"} - $j2->{"$j2->{_unique_interval}\_position"}) > 20);

		#the first unique coords are going into the IS element
		my $uc1_strand = $j1->{"$j1->{_unique_interval}\_strand"};
		my $uc2_strand = $j2->{"$j2->{_unique_interval}\_strand"};
		next if ($uc1_strand != -$uc2_strand);

		my $is1_strand = - $j1->{"$j1->{_is_interval}\_strand"} * $j1->{"_$j1->{_is_interval}\_is"}->{strand} * $j1->{"$j1->{_unique_interval}\_strand"};
		my $is2_strand = - $j2->{"$j2->{_is_interval}\_strand"} * $j2->{"_$j2->{_is_interval}\_is"}->{strand} * $j2->{"$j2->{_unique_interval}\_strand"};
		my $is_strand = ($is1_strand == $is2_strand) ? $is2_strand : '0';

		### add additional information to the first match, which will 
		### cause a new line to be drawn in the new junction table

		splice @jc, $i, 2; 
		$i-=2;

		my $mut = { 
			type => 'MOB',
			seq_id => $j1->{"$j1->{_unique_interval}\_seq_id"},
			evidence => [ $j1->{id}, $j2->{id} ],
		};

		$mut->{start} = ($uc1_strand == -1) ? $j2->{"$j2->{_unique_interval}\_position"} : $j1->{"$j1->{_unique_interval}\_position"};
		$mut->{end} = ($uc1_strand == -1) ? $j1->{"$j1->{_unique_interval}\_position"} : $j2->{"$j2->{_unique_interval}\_position"};
		$mut->{repeat_name} = $j1->{"_$j1->{_is_interval}\_is"}->{name};
		$mut->{strand} = $is_strand;
		
		##this is an estimate for the size of the entire element
### TODO check the size range for the element across the reference genome!!		
		$mut->{size} = abs($j1->{"_$j1->{_is_interval}\_is"}->{end} - $j1->{"_$j1->{_is_interval}\_is"}->{start} + 1);
		$mut->{position} = $mut->{start} - 1;
		$mut->{duplication_size} = $mut->{end} - $mut->{start} + 1;
		
		#sometimes the ends of the IS are not quite flush		
		if ($j1->{"$j1->{_is_interval}\_strand"} == -1)
		{
			$mut->{gap_left} = $j1->{"$j1->{_is_interval}\_position"} - $j1->{"_$j1->{_is_interval}\_is"}->{end} + abs($j1->{overlap});
		}
		else
		{
			$mut->{gap_left} = $j1->{"_$j1->{_is_interval}\_is"}->{start} - $j1->{"$j1->{_is_interval}\_position"} + abs($j1->{overlap});
		}

		if ($j2->{"$j2->{_is_interval}\_strand"} == -1)
		{
			$mut->{gap_right} = $j2->{"$j2->{_is_interval}\_position"} - $j2->{"_$j2->{_is_interval}\_is"}->{end} + abs($j2->{overlap});
		}
		else
		{
			$mut->{gap_right} = $j2->{"_$j2->{_is_interval}\_is"}->{start} - $j2->{"$j2->{_is_interval}\_position"} + abs($j2->{overlap});
		}

		if ($j1->{"$j1->{_unique_interval}\_strand"} *  $j1->{"_$j1->{_unique_interval}\_read_side"} == +1)
		{
			($mut->{gap_right}, $mut->{gap_left}) = ($mut->{gap_left}, $mut->{gap_right});
		}		
		
		$gd->add($mut);
	}

	

	##RA that overlap deletions should not be shown
	my @del = $gd->list('DEL');	
	RA: foreach my $ra_item (@ra)
	{
		DEL: foreach my $del_item (@del)
		{
			next DEL if ($ra_item->{seq_id} ne $del_item->{seq_id});
			
			## there might be a problem here with insert_position > 0
			if ( ($ra_item->{position} >= $del_item->{position}) && ($ra_item->{position} <= $del_item->{position} + $del_item->{size} - 1) )
			{
				$ra_item->{deleted} = 1;
				next RA;
			}
		}
	}
	
	
	
	## PROBLEM: We can't apply the coverage cutoff until AFTER we count errors
	##   (because only then do we have the distribution to fit)
	##   but we have to choose which junctions we believe BEFORE counting
	##   (because we put their split alignments in the BAM file)
	## Ideally we would do this after step 7, then remove the offending read pieces from the BAM file
	## before proceeding to SNP calling.
	## this could be done by reserving these pieces in a separate SAM file
	## then merging them later? But full matches would also have to be kept separate...
	##

	## Remove remaining junctions that we didn't pair up with anything that are below a coverage cutoff.
	@jc = $gd->filter_used_as_evidence($gd->list('JC'));	
	foreach my $item (@jc)
	{		
		my $coverage_cutoff_1 = $settings->{unique_coverage}->{$item->{side_1_seq_id}}->{junction_coverage_cutoff};
		my $coverage_cutoff_2 = $settings->{unique_coverage}->{$item->{side_2_seq_id}}->{junction_coverage_cutoff};
		
		if ( (!defined $coverage_cutoff_1 || ($item->{total_reads} < $coverage_cutoff_1) ) 
		  && (!defined $coverage_cutoff_2 || ($item->{total_reads} < $coverage_cutoff_2) ) )
		{
			Breseq::GenomeDiff::add_reject_reason($item, "COV");
		}
	}
	
}


return 1;
