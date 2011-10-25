###
# Pod Documentation
###

=head1 NAME

RNAfold::FoldedSequence.pm

=head1 SYNOPSIS

An object for reading an Avida output file with named columns. A list of
hashes with the column names as keys is created.

=head1 AUTHOR

Jeffrey Barrick

=head1 COPYRIGHT

Copyright (C) 2007.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

=cut

###
# End Pod Documentation
###

package GenomeDiff::PreterminationTables;
use vars qw(@ISA);
use strict;

use Bio::Seq;

@ISA = qw();

use Data::Dumper;

=head2 new

 Title   : new
 Usage   : $af = RNAfold::PreterminationTables->new( -filename => 'input.sto' );
 Function: pass in a sequence and number of requested folds
 Returns : a new AvidaFile object

 Options: -separate_codon_blocks 
		  treat occurences of an amino acid that differ at a nt
		  other than the third codon position as separate amino acids
=cut

sub new
{
	my($caller,%options) = @_;
	
	#Default to bacterial genetic code
	$options{translation_table} = 11 if ( !defined $options{translation_table} );	
	$options{transition_to_transversion_ratio} = 4.0 if ( !defined $options{transition_to_transversion_ratio} );
	
	my $normalize_mutation_rate;

	##decide whether we are using HKY85 model
	if (defined($options{codon_3_composition}))
	{
		$normalize_mutation_rate = (2 + $options{transition_to_transversion_ratio});
	}
	##or using a uniform distribution
	else
	{
		$options{codon_3_composition} = { 'A'=> 0.25, 'T'=> 0.25, 'C'=> 0.25, 'G'=> 0.25 };
		$normalize_mutation_rate = (2+$options{transition_to_transversion_ratio});
	}
	
	my $self = { }; #and I
	my $class = ref($caller) || $caller;
	bless ($self, $class);
	
	#Construct the table	
	my @nt_list = ('A', 'G', 'T', 'C');
	my %transition_if_equal = ( 'A' => 1, 'G' => 1, 'T' => 0, 'C' => 0); #parity table
	
	#Information we will be saving
	my %codon_to_aa;
	my %aa_to_codons;
	my %pretermination_codons;
	my %diff_pretermination_aa;
	my %pretermination_aa;
	my %gc3_same_aa_codons_preterm;
	my %gc3_same_aa_codons_no_preterm;
	my %gc3_no_preterm;
	my %gc3_equilibrium;
	
	#calls to bioperl are slow, do only once
	my %basic_codons_to_aa;
	foreach my $nt_1 (@nt_list)
	{
		foreach my $nt_2 (@nt_list)
		{	
			foreach my $nt_3 (@nt_list)
			{
				my $codon = $nt_1 . $nt_2 . $nt_3;
				
				#if it is a stop codon, then we're done
				my $codon_seq = Bio::Seq->new('-seq' => $codon);
				my $aa = $codon_seq->translate( undef, undef, undef, $options{translation_table} )->seq();
				$basic_codons_to_aa{$codon} = $aa;
			}
		}
	}
	
	foreach my $nt_1 (@nt_list)
	{
		foreach my $nt_2 (@nt_list)
		{
			my $all_gc3_are_same_aa = 1;
			my $gc_3_const_aa;
			my $gc_3_includes_preterm = 0;
		
			foreach my $nt_3 (@nt_list)
			{
				my $codon = $nt_1 . $nt_2 . $nt_3;
				
				#if it is a stop codon, then we're done
				#my $codon_seq = Bio::Seq->new('-seq' => $codon);
				#my $aa = $codon_seq->translate( undef, undef, undef, $options{translation_table} )->seq();
				my $aa = $basic_codons_to_aa{$codon};
				
				##Optionally, look at other codons we already have for this amino acid
				##and give a new group if it differs at a nt other than the wobble position.
				if ($options{separate_codon_blocks_3})
				{
					#try all blocks to see if we match one of them.
					my $block_num = 1;
					BLOCK: while (defined $aa_to_codons{$aa . $block_num})
					{						
						foreach my $test_codon (@{$aa_to_codons{$aa . $block_num}})
						{
							my $copy_codon = $test_codon;
							#print STDERR "$codon $test_codon $copy_codon\n";
							$copy_codon =~ s/.$//;
							#print STDERR "$codon $test_codon $copy_codon\n";
							last BLOCK if ($codon =~ m/^$copy_codon/);
						}
						
						$block_num++;
					}
					
					$aa .= $block_num;
				}
				elsif ($options{separate_codon_blocks})
				{
					#try all blocks to see if we match one of them.
					my $block_num = 1;
					BLOCK: while (defined $aa_to_codons{$aa . $block_num})
					{						
						foreach my $test_codon (@{$aa_to_codons{$aa . $block_num}})
						{
							my $copy_codon = $test_codon;
							#print STDERR "$codon $test_codon $copy_codon\n";
							$copy_codon =~ s/.$//;
							#print STDERR "$codon $test_codon $copy_codon\n";
							last BLOCK if ($codon =~ m/^$copy_codon/);
						}
						
						$block_num++;
					}
					
					$aa .= $block_num;
				}
				
				if (!defined $gc_3_const_aa)
				{
					$gc_3_const_aa = $aa;
				}
				elsif ($gc_3_const_aa ne $aa)
				{
					$all_gc3_are_same_aa = 0;
				}

				$codon_to_aa{$codon} = $aa;
				push @{$aa_to_codons{$aa}}, $codon;
				
				#we don't count stop codons
				next if ($aa =~ m/^\*/);
				
				#try all one-neighbors
				for (my $pos = 0; $pos < 3; $pos++)
				{
					NB_NT: foreach my $nb_nt (@nt_list)
					{
						my $nb_codon = $codon;
						my $old_nt = substr ($nb_codon, $pos, 1, $nb_nt);
						
						#print "$codon $nb_codon $old_nt\n";
						
						#bail if same codon
						next NB_NT if ($old_nt eq $nb_nt);
										
						#my $nb_codon_seq = Bio::Seq->new('-seq' => $nb_codon);
						#my $nb_aa = $nb_codon_seq->translate( undef, undef, undef, $options{translation_table} )->seq();
						my $nb_aa = $basic_codons_to_aa{$nb_codon};
										
						#bail if not a stop codon
						next if ($nb_aa ne '*');
						
						my $add_prob = $options{codon_3_composition}->{$nb_nt};
						$add_prob *= ($transition_if_equal{$nb_nt} == $transition_if_equal{$old_nt}) 
							? $options{transition_to_transversion_ratio} : 1;
						
						$pretermination_codons{$codon} = 0 if (!defined $pretermination_codons{$codon});
						$pretermination_codons{$codon} += $add_prob / $normalize_mutation_rate;
					}
					#end all one-neighbors	
				
				}	
				
				#if one of the codons with the same gc3 is a pretermination,
				#then this is no longer a background codon/aa
				#or if the amino acid encoded has changed
				if ($pretermination_codons{$codon} || !$all_gc3_are_same_aa)
				{
					$gc_3_includes_preterm = 1;	
				}
			}	
			
			
			#record gc3 stats
			if ($all_gc3_are_same_aa)
			{
				foreach my $nt_3 (@nt_list)
				{
					my $codon = $nt_1 . $nt_2 . $nt_3;
					
					if ($gc_3_includes_preterm)
					{
						$gc3_same_aa_codons_preterm{$codon} = 1;
					}
					else
					{
						$gc3_same_aa_codons_no_preterm{$codon} = 1;
					}
				}
			}
			
			if (!$gc_3_includes_preterm)
			{
				foreach my $nt_3 (@nt_list)
				{
					my $codon = $nt_1 . $nt_2 . $nt_3;
					$gc3_no_preterm{$codon} = 1;
				}		
			}
			if ($all_gc3_are_same_aa)
			{
				foreach my $nt_3 (@nt_list)
				{
					my $codon = $nt_1 . $nt_2 . $nt_3;
					$gc3_equilibrium{$codon} = 1;
				}
			}
			
		}	
	}

	#Decide which amino acids have codons that are informative when amino acid content is unchanged
	AA: foreach my $key (keys %aa_to_codons)
	{
		#Don't count stops as informative
		next if ($key =~ m/^\*/);

		my $codon_value;
		foreach my $codon (@{$aa_to_codons{$key}})
		{
			my $this_codon_value = $pretermination_codons{$codon};
			$this_codon_value = 0 if (!defined $this_codon_value);
			
			if ($this_codon_value > 0)
			{
				$pretermination_aa{$key} = 1;
			}
			
			#if any are different, then this amino acid is informative. RLSG by default no matter what TTR ratio is.
			if ((defined $codon_value) && ($this_codon_value != $codon_value))
			{
				$diff_pretermination_aa{$key} = 1;

				next AA;
			}
			$codon_value = $this_codon_value;
		}	
		#Remove codons where there is no difference between chance over all codons.
		delete $diff_pretermination_aa{$key};
	}

	#save all of the tables we have calculated
	$self->{codon_to_aa} = \%codon_to_aa;
	$self->{aa_to_codons} = \%aa_to_codons;
	$self->{pretermination_codons} = \%pretermination_codons;
	$self->{diff_pretermination_aa} = \%diff_pretermination_aa;
	$self->{pretermination_aa} = \%pretermination_aa;
	$self->{gc3_same_aa_codons_preterm} = \%gc3_same_aa_codons_preterm;
	$self->{gc3_same_aa_codons_no_preterm} = \%gc3_same_aa_codons_no_preterm;
	$self->{gc3_no_preterm} = \%gc3_no_preterm;
	$self->{gc3_equilibrium} = \%gc3_equilibrium;
	$self->{codon_3_composition} = $options{codon_3_composition};
	$self->{transition_if_equal} = \%transition_if_equal;
	$self->{transition_to_transversion_ratio} = $options{transition_to_transversion_ratio};
	$self->{normalize_mutation_rate} = $normalize_mutation_rate;
	return $self;
}


sub relative_mutation_rate
{
	my ($self, $from, $to) = @_;

	my $rate = $self->{codon_3_composition}->{$to};
	$rate *= ($self->{transition_if_equal}->{$from} == $self->{transition_if_equal}->{$to}) 
		? $self->{transition_to_transversion_ratio} : 1;
	return $rate / $self->{normalize_mutation_rate};
}


return 1;
