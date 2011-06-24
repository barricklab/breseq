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

package GenomeDiff::SynonymousNonsynonymous;
use vars qw(@ISA);
use strict;

use Bio::Seq;

@ISA = qw();

use Data::Dumper;

=head2 new

 Title   : new
 Usage   : $sn = Codon::SynonymousNonsynonymous>new( -filename => 'input.sto' );
 Function: 
 Returns : Options = {count_synonymous_stop_codons}

=cut

sub new
{
	my($caller,%options) = @_;
	
	#Default to bacterial genetic code
	$options{translation_table} = 11 if ( !defined $options{translation_table} );	
	$options{count_synonymous_stop_codons} = $options{count_synonymous_stop_codons};	

	my $self = { }; #and I
	my $class = ref($caller) || $caller;
	bless ($self, $class);
	
	#Construct the table	
	my @nt_list = ('A', 'G', 'T', 'C');
	
	#Information we will be saving
	my %codon_to_aa;
	my %aa_to_codons;
	my %codon_synonymous_changes;
	my %codon_nonsynonymous_changes;
	my %codon_num_synonymous_changes;
	my %codon_num_mutT_synonymous_changes;
	
	## structure: contains keys of form A/C/T/G_A/C/T/G and values are number of such changes per codon.
	
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
				$codon_to_aa{$codon} = $aa;
				push @{$aa_to_codons{$aa}}, $codon;

			}
		}
	}
	
	##zero counts so nothing is undefined
	foreach my $nt_1 (@nt_list)
	{
		foreach my $nt_2 (@nt_list)
		{
			foreach my $codon (keys %codon_to_aa)
			{
				$codon_synonymous_changes{$codon}->{"$nt_1-$nt_2"} = 0;
				$codon_nonsynonymous_changes{$codon}->{"$nt_1-$nt_2"} = 0;
				$codon_num_synonymous_changes{$codon} = 0;
				$codon_num_mutT_synonymous_changes{$codon} = 0;
			}
		}
	}
	
	foreach my $nt_1 (@nt_list)
	{
		foreach my $nt_2 (@nt_list)
		{		
			foreach my $nt_3 (@nt_list)
			{
				my $codon = $nt_1 . $nt_2 . $nt_3;
				
				#if it is a stop codon, then we're done
				#my $codon_seq = Bio::Seq->new('-seq' => $codon);
				#my $aa = $codon_seq->translate( undef, undef, undef, $options{translation_table} )->seq();
				my $aa = $codon_to_aa{$codon};
				
				#we don't count stop codons
				next if (($aa =~ m/^\*/) && (!$options{count_synonymous_stop_codons}));
				
				#try all one-neighbors
				for (my $pos = 0; $pos < 3; $pos++)
				{
					NB_NT: foreach my $nb_nt (@nt_list)
					{
						my $nb_codon = $codon;
						my $old_nt = substr ($nb_codon, $pos, 1, $nb_nt);
											
						#bail if same codon
						next NB_NT if ($old_nt eq $nb_nt);
										
						#my $nb_codon_seq = Bio::Seq->new('-seq' => $nb_codon);
						#my $nb_aa = $nb_codon_seq->translate( undef, undef, undef, $options{translation_table} )->seq();
						my $nb_aa = $codon_to_aa{$nb_codon};
										
						## this is an interesting question, what do we do about stop codons?
						#for now we count them as nonsynonymous
						#next if ($nb_aa ne '*');
						
						if ($aa eq $nb_aa)
						{
							$codon_synonymous_changes{$codon}->{"$old_nt-$nb_nt"}++;
							$codon_num_synonymous_changes{$codon}++;
							
							if ( (($old_nt eq 'A') && ($nb_nt eq 'C'))
							  || (($old_nt eq 'T') && ($nb_nt eq 'G')) )
							{
								$codon_num_mutT_synonymous_changes{$codon}++;
							}
						}
						else
						{
							$codon_nonsynonymous_changes{$codon}->{"$old_nt-$nb_nt"}++;
						}
					}
					#end all one-neighbors	
				}	
			}	
		}	
	}
	
	my %codon_position_mutation_synonymous;

	foreach my $from_codon (sort keys %codon_to_aa)
	{	
		foreach my $pos (1..3)
		{
			my $from_nt = substr($from_codon, $pos-1, 1);
			my $from_aa = $codon_to_aa{$from_codon};
			foreach my $to_nt (@nt_list)
			{				
				my $new_codon = $from_codon;
				substr($new_codon, $pos-1, 1) = $to_nt;
				my $new_aa = $codon_to_aa{$new_codon};
				
				if ($from_aa eq $new_aa)
				{
					$codon_position_mutation_synonymous{$from_codon . "_" . $pos . "_" . $from_nt . "_" . $to_nt} = 1;
				}
			}
		}
	}
	
	
	#save all of the tables we have calculated
	
	$self->{codon_to_aa} = \%codon_to_aa;
	$self->{aa_to_codons} = \%aa_to_codons;
	$self->{codon_synonymous_changes} = \%codon_synonymous_changes;
	$self->{codon_nonsynonymous_changes} = \%codon_nonsynonymous_changes;
	$self->{codon_num_synonymous_changes} = \%codon_num_synonymous_changes;
	$self->{codon_num_mutT_synonymous_changes} = \%codon_num_mutT_synonymous_changes;
	
	$self->{codon_position_mutation_synonymous} = \%codon_position_mutation_synonymous;
	
	return $self;
}


return 1;
