###
# Pod Documentation
###

=head1 NAME

GenomeDiff::BaseSubstitutionFile.pm

=head1 SYNOPSIS

Keeps information about the effects of base substitions in a genome.

=head1 AUTHOR

Jeffrey Barrick

=head1 COPYRIGHT

Copyright (C) 2011.

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

package GenomeDiff::BaseSubstitutionFile;
use vars qw(@ISA);
use strict;

use Bio::Seq;
use Bio::Tools::SeqStats;
use GenomeDiff::SynonymousNonsynonymous;

@ISA = qw();

use Data::Dumper;


# Conventions for base change labels

our @nt_list = ('A', 'G', 'T', 'C');

our $bp_change_to_index = {
	"A-G" => 0, "T-C" => 0,
	"A-C" => 1, "T-G" => 1,
	"A-T" => 2, "T-A" => 2,
	"G-A" => 3, "C-T" => 3,
	"G-T" => 4, "C-A" => 4,
	"G-C" => 5, "C-G" => 5,	
};

our @bp_change_list = (
	["A-G", "T-C"],
	["A-C", "T-G"],
	["A-T", "T-A"],
	["G-A", "C-T"],
	["G-T", "C-A"],
	["G-C", "C-G"],
);

our @bp_change_label_list = (
	"AT.GC",
	"AT.CG",
	"AT.TA",
	"CG.TA",
	"CG.AT",
	"CG.GC",
);

our $base_to_bp_change = {
	"A.G" => "AT.GC", "T.C" => "AT.GC",
	"A.C" => "AT.CG", "T.G" => "AT.CG",
	"A.T" => "AT.TA", "T.A" => "AT.TA",
	"C.T" => "CG.TA", "G.A" => "CG.TA",
	"C.A" => "CG.AT", "G.T" => "CG.AT",
	"C.G" => "CG.GC", "G.C" => "CG.GC",
};

our @bsf_snp_types = ( 'PROTEIN', 'RNA', 'PSEUDOGENE', 'INTERGENIC', 'NONSYNONYMOUS', 'SYNONYMOUS', 'TOTAL');

our $k_nt_type_PROTEIN = 0;
our $k_nt_type_RNA = 1;
our $k_nt_type_PSEUDOGENE = 2;
our $k_nt_type_INTERGENIC = 3;
our $k_nt_type_MAX = 4;



=head2 new

 Title   : new
 Usage   : $sn = Codon::SynonymousNonsynonymous>new( -filename => 'input.sto' );
 Function: 
 Returns : Options = {count_synonymous_stop_codons}

=cut

sub new
{
	my ($caller,%options) = @_;
	my $count_synonymous_stop_codons = 1;
	my $verbose = 0;
	
	my $self = { }; #and I
	my $class = ref($caller) || $caller;
	bless ($self, $class);
	
	my $seq = $options{sequence};
	defined ($seq) or die "Undefined sequence.";
		
	my $nss;
	my %codon_to_aa;
	my %aa_to_codons;
	my %codon_synonymous_changes;
	my %codon_nonsynonymous_changes;
	my %codon_num_synonymous_changes;
	my %codon_position_mutation_synonymous;

	my %nonsynonymous_mutations;
	my %synonymous_mutations;

	my $total_num_synonymous_changes = 0;
	my $total_num_nonsynonymous_changes = 0;
	my $total_codon_nt_positions = 0;
	my $total_nt_position = 0;

	my $translation_table;

	#print STDERR Dumper(\%highly_expressed);

	my $total_bases;

	##Process sequences
	my @fragments = ();
	my @proteins = ();
	my $total_codons = 0;
	my $total_orfs = 0;

	#lower numbers get preference

	#Load sequence
	my $seq_length = $seq->length; #BioPerl is slow...
	$self->{seq_length} = $seq_length;

	my $monomers = Bio::Tools::SeqStats->count_monomers($seq);
	foreach my $base (sort keys %$monomers) 
	{
		$total_bases->{$base} += $monomers->{$base};
    	#print STDERR "Number of bases of type ", $base, " = ", $monomers->{$base},"\n";
    	print STDERR "WARNING: Nonstandard base in sequence: $base\n" if (($base ne 'A') && ($base ne 'T') && ($base ne 'C') && ($base ne 'G'));
	}

	$total_nt_position += $seq->length;

	print STDERR "Allocating genome arrays...\n" if ($verbose);

	my @nt_change_is_nonsynonymous;
	$#nt_change_is_nonsynonymous = $seq_length * 6;  		
	for (my $i=0; $i<$seq_length * 6; $i++) {
		$nt_change_is_nonsynonymous[$i] = 0;
	}

	my @nt_type;
	$#nt_type = $seq_length;  
	for (my $i=0; $i<$seq_length; $i++) {
		$nt_type[$i] = $k_nt_type_INTERGENIC;
	}

	my $fragment;
	$fragment->{accession_version} = $seq->accession() . "." . $seq->version();
	$fragment->{description} = $seq->description();
	print STDERR "$fragment->{description}\n" if ($verbose);		
	push @fragments, $fragment;

	my @feature_list = $seq->get_SeqFeatures(); # just top level features

	##Record proteins
	my $prev_gene;
	FEATURE: foreach my $f (@feature_list) 
	{		
		## figure out position
		my @Location_List = $f->location->each_Location();

		## Record three categories: noncoding (RNA genes), pseudogenes, protein coding genes, intergenic (last is assumed if others not true)
		if (($f->primary_tag eq 'tRNA') || ($f->primary_tag eq 'rRNA') ) {

			foreach my $loc (@Location_List) {
				for (my $i=$loc->start; $i<=$loc->end; $i++) {
					$nt_type[$i-1] = $k_nt_type_RNA if ($nt_type[$i-1] > $k_nt_type_RNA);
				}
			}
			next FEATURE;				
		}


		if ($f->has_tag("pseudo") || ( ($f->location->start_pos_type ne 'EXACT') || ($f->location->end_pos_type ne 'EXACT')) ) {

			foreach my $loc (@Location_List) {
				for (my $i=$loc->start; $i<=$loc->end; $i++) {
					$nt_type[$i-1] = $k_nt_type_PSEUDOGENE if ($nt_type[$i-1] > $k_nt_type_PSEUDOGENE);						
				}
			}
			next FEATURE;
		}
		
		next FEATURE if ($f->primary_tag ne 'CDS');
		

		#no gene fragments!

		$total_orfs++;

		my $p; 
		$p->{start} = $Location_List[0]->start if (!defined $p->{start});
		$p->{end} = $Location_List[$#Location_List]->end;
		$p->{strand} =  $f->location->strand;			$p->{gene} = GetTag($f, "gene");
		$p->{name} = ($p->{gene}) ? $p->{gene} : GetTag($f, "locus_tag");
		$p->{description} = GetTag($f, "product");
		$p->{locus_tag} = GetTag($f, "locus_tag");
		$p->{gene} = '-' if (!$p->{gene});

		#print STDERR "$p->{gene}\n";

		#Load a list of synonyms, format is likely to vary in each file...
		push @{$p->{synonyms}}, $p->{gene} if ($p->{gene});
		push @{$p->{synonyms}}, $p->{locus_tag} if ($p->{locus_tag});	



		## Use the correct genetic code
		## Generate new pretermination table if this one has not been seen yet.
		$p->{translation_table} = GetTag($f, "transl_table");
		die "No translation table for CDS!" if (!$p->{translation_table});
		die if (defined $translation_table && ($translation_table != $p->{translation_table}));

		if (!defined $translation_table)
		{
			$translation_table = $p->{translation_table};
			print STDERR "Using translation table: $p->{translation_table}\n" if ($verbose);

			$nss =  GenomeDiff::SynonymousNonsynonymous->new(	
													'translation_table' => $p->{translation_table}, 
													'count_synonymous_stop_codons' => $count_synonymous_stop_codons,
												);
			%codon_synonymous_changes = %{$nss->{codon_synonymous_changes}};
			%codon_num_synonymous_changes = %{$nss->{codon_num_synonymous_changes}};
			%codon_nonsynonymous_changes = %{$nss->{codon_nonsynonymous_changes}};
			%codon_to_aa = %{$nss->{codon_to_aa}};
			%aa_to_codons = %{$nss->{aa_to_codons}};
			%codon_position_mutation_synonymous = %{$nss->{codon_position_mutation_synonymous}};

			#print Dumper(%codon_position_mutation_synonymous);
			#print STDERR Dumper(%codon_num_synonymous_changes);
			#print STDERR Dumper(%codon_synonymous_changes, %codon_nonsynonymous_changes, %codon_to_aa, %aa_to_codons);

		}
		elsif ($translation_table != $p->{translation_table})
		{
			#no need to die, but downstream programs need to worry about this
			die "Different translation tables used in the same organism! $translation_table and $p->{translation_table}\n";
		}

		# Piece together the gene
		# Code should be savvy to an internal
		# frameshift (bacterial) or intron (eukaryotic).
		$p->{nt_seq} = '';
		my @gene_positions; #nt positions in referencethat correspond to each base of nt_seq
		foreach my $loc (@Location_List)
		{
			#print STDERR $loc->start . " " .  $loc->end . "\n";
			
			my $add_seq = $seq->trunc($loc->start, $loc->end);
			my $on_nt;
			if ($p->{strand} == -1)
			{
				$on_nt = $loc->end;
				$add_seq = $add_seq->revcom;
				$p->{nt_seq} = $add_seq->seq . $p->{nt_seq};
			}
			else
			{
				$on_nt = $loc->start;
				$p->{nt_seq} = $p->{nt_seq} . $add_seq->seq;
			}

			my @new_gene_positions;
			my $on_nt_index = 0;
			for (my $pos = $loc->start-1;  $pos < $loc->end; $pos++)
			{
				$new_gene_positions[$on_nt_index] = $on_nt-1;
				$on_nt_index++;
				$on_nt += $p->{strand};
			}
			if ($p->{strand} == +1) {
				push @gene_positions, @new_gene_positions;
			} else {
				unshift @gene_positions, @new_gene_positions;
			}
		}

		#print STDERR Dumper($p);


		#Split to codons
		my $nt_seq = $p->{nt_seq};
		my @aa_list;
		my $codon_index = -1;
		CODON: while (my $codon = substr $nt_seq, 0, 3, "") {
			$codon_index++;

			#add to codon list and keep track of codons used
			if (length $nt_seq == 0) {
				if ($codon_to_aa{$codon} ne '*') {
					print STDERR "WARNING: Reading frame has no stop codon: $p->{gene}|$p->{locus_tag}\n";
				}
			}

			#check for stop codon
			if ($codon_to_aa{$codon} eq '*') {
				if (length $nt_seq > 0) {
					print STDERR "WARNING: Stop codon ($codon) is within gene: $p->{gene}|$p->{locus_tag}\n";
				}
				## count or not?
				next CODON if (!$count_synonymous_stop_codons);
			}

			push @{$p->{codon_list}}, $codon;
			push @aa_list, $codon_to_aa{$codon};

			#print $nt_seq . "\n";
			#print $codon . "\n";

			if (length $codon != 3) {
				print STDERR "ERROR: Codon that is not of length 3 found!\n";
				print STDERR $nt_seq . "\n";
				print STDERR $codon . "\n";
			 	print STDERR Dumper($p);
			 	die;
			}

			print "  $codon\n" if ($verbose);

			## update genome synonymous calls
			## Note that stop codon has not been removed!
			foreach my $codon_position (1..3)
			{ 
				my $gene_position = $codon_index * 3 + $codon_position - 1;
				my $genome_position = $gene_positions[$gene_position];

				my $from_nt = substr $p->{nt_seq}, $gene_position, 1;
				## note -> this is already on the coding strand

				print "  $codon_position:$from_nt ($gene_position, $genome_position)\n" if ($verbose);

				## count this position as in a gene
				$nt_type[$genome_position] = $k_nt_type_PROTEIN if ($nt_type[$genome_position] > $k_nt_type_PROTEIN);


				TO_NT: foreach my $to_nt (@nt_list) {
					next TO_NT if ($from_nt eq $to_nt);
					#print STDERR "$genome_position $codon_index $gene_position $from_nt-$to_nt" . "\n";

					if (!$codon_position_mutation_synonymous{$codon . "_" . $codon_position . "_" . $from_nt . "_" . $to_nt}) {	
							$nt_change_is_nonsynonymous[$genome_position * 6 + $bp_change_to_index->{"$from_nt-$to_nt"} ] = 1;
						}
				}		
			} # end codon position

		} ## end codon



		#Save length in AA
		$p->{length} = scalar @{$p->{codon_list}};

		#Tabulate all codons used
		CODON: for (my $on_codon = 0; $on_codon < scalar @{$p->{codon_list}}; $on_codon++)
		{
			my $codon = $p->{codon_list}->[$on_codon];

			#Need to catch these with /transl_except=(pos:1546010..1546012,aa:Sec)

			if (!defined $codon_to_aa{$codon})
			{
				print STDERR "Undefined translation for codon ($codon) in the middle of this gene: $p->{gene}|$p->{locus_tag}\n";
				#print STDERR Dumper($p);
				next CODON;
			}
			elsif ( ($codon_to_aa{$codon} eq '*') && ($on_codon != scalar @{$p->{codon_list}}) - 1 )
			{
				print STDERR "Stop codon ($codon) in the middle of this gene: $p->{gene}|$p->{locus_tag}\n";
				#print STDERR Dumper($p);
				next CODON;
			} 
		}
		
		
		## debug code
		#print ">$p->{name}\n";
		#print +(join "", @{$p->{codon_list}}) . "\n";
		#print +(join "", @aa_list) . "\n";
			
	} ## end feature

	###
	#  Per-Sequence Counting Loop
	###

	my $sequence = $seq->seq;
	print STDERR "Tallying per-position information...\n" if ($verbose);

	## Now go through the sequence and normalize expectations for overlapping genes
	## If it is marked as nonsynonymous, then it's nonsynonymous even if others are synonymous
	for (my $genome_position=0; $genome_position<$seq_length; $genome_position++) {
		
		my $from_nt = substr $sequence, $genome_position, 1;
		
		
		## debug code for consistency with C++
		#print +($genome_position+1) . " " . $from_nt;
		#TO_NT: foreach my $to_nt ('A', 'T', 'C', 'G') {
		#	#print +($genome_position+1) . " " . $from_nt . " " . $to_nt . " " . (($from_nt != $to_nt) ? "1" : "0") . " " . (($nt_change_is_nonsynonymous[$genome_position * 6 + $bp_change_to_index->{"$from_nt-$to_nt"} ] == 1) ? "1" : "0") . "\n";
		#	print " " . ((($from_nt ne $to_nt) && ($nt_change_is_nonsynonymous[$genome_position * 6 + $bp_change_to_index->{"$from_nt-$to_nt"} ] == 1)) ? "1" : "0");
		#}
		#print "\n";
		## end consistency code
		
		#my $print_num;
		#if ($nt_type[$genome_position] == $k_nt_type_PROTEIN) {
		#	$print_num = 3;
		#}
		#if ($nt_type[$genome_position] == $k_nt_type_RNA) {
		#	$print_num = 1;
		#}
		#if ($nt_type[$genome_position] == $k_nt_type_PSEUDOGENE) {
		#	$print_num = 2;
		#}
		#if ($nt_type[$genome_position] == $k_nt_type_INTERGENIC) {
		#	$print_num = 0;
		#}
		#print " " . $print_num . "\n";
		
		
		## only count coding positions
		next if ($nt_type[$genome_position] != $k_nt_type_PROTEIN);
		$total_codon_nt_positions++;


		TO_NT: foreach my $to_nt (@nt_list) {

			next TO_NT if ($from_nt eq $to_nt);
			#print STDERR "$genome_position $codon_index $gene_position $from_nt-$to_nt" . "\n";

			my $key = "$from_nt-$to_nt";


			if ($nt_change_is_nonsynonymous[$genome_position * 6 + $bp_change_to_index->{"$from_nt-$to_nt"} ] == 1) {
				$nonsynonymous_mutations{$key}++;
				$total_num_nonsynonymous_changes++;
			} else {
				$synonymous_mutations{$key}++;
				$total_num_synonymous_changes++;
			}
		}
	}

	###
	#  / Per-Sequence Counting Loop
	###

	my @syn_labels = ("syn-AT.GC-CG.TA", "syn-AT.CG-CG.AT", "syn-AT.TA-CG.GC");
	@syn_labels = map { $_ . ".syn" } @syn_labels;
	my @tot_labels = ('Type');

	foreach (my $i=0; $i<$seq_length; $i++)
	{
		my @tot_list;

		my $from_nt = substr $sequence, $i, 1;

		my $is_AT = 0;
		$is_AT = (($from_nt eq 'A') || ($from_nt eq 'T')) ? 1 : 0;

		push @tot_list, $nt_type[$i];		

		my @syn_list = ();
		if ($nt_type[$i] == $k_nt_type_PROTEIN) {

			my $j_offset = 0;
			$j_offset = 3 if (!$is_AT);
			for (my $j=0; $j < 3; $j++) {
				push @syn_list, ($nt_change_is_nonsynonymous[$i*6+$j+$j_offset] == 0) ? 1 : 0;
			}
		}
	}
	
	#save
	$self->{packed_data} = '';
	for (my $i=0; $i < $self->{seq_length}; $i++) 
	{	
		my $from_nt = substr $sequence, $i, 1;
		my $is_AT = 0;
		$is_AT = (($from_nt eq 'A') || ($from_nt eq 'T')) ? 1 : 0;
		my $offset = 0;
		$offset = 3 if (!$is_AT);

		my $start = $i*6 + $offset;
		my $end = $start+2;
		
		my $bit_string = '';
		
		my ($this_nt_type, $this_bit);
		my $this_nt_type = $nt_type[$i];
		
		$this_bit = $this_nt_type % 2;
		$bit_string .= $this_bit;
		$this_nt_type /= 2;
		$this_bit = $this_nt_type % 2;
		$bit_string .= $this_bit;

		#print +($i+1) . " " . $nt_type[$i] . "\n";
		#print "$bit_string\n";
		for (my $j=$start; $j<=$end;$j++) 
		{
			$bit_string .= $nt_change_is_nonsynonymous[$j];
		}
		#print "$bit_string\n";
		$self->{packed_data} .= pack 'b8', $bit_string;
	}
	
	return $self;
}

sub read
{
	my ($caller,%options) = @_;
	my $count_synonymous_stop_codons = 1;
	my $verbose = 0;
	
	my $self = { }; #and I
	my $class = ref($caller) || $caller;
	bless ($self, $class);
	
	$self->{nt_sequence} = $options{nt_sequence};
	defined ($self->{nt_sequence}) or die "Undefined nt_sequence.";
	$self->{nt_sequence} = "\U$self->{nt_sequence}";
	
	my $input_file = $options{input_file};
	defined ($input_file) or die "Undefined input_file.";
	
	open IN, "<$input_file" or die "Could'd open file $input_file";
	$self->{packed_data} = '';
	while (<IN>) {
		$self->{packed_data} .= $_;
	}
	close IN;
	
	return $self;
}
	
=comment

	## print out the synonymous / nonsynonymous table
	print OUT "#Probabilities of synonymous mutations (NCBI Translation Table $translation_table)\n";
	print OUT +join("\t", 'aa', 'codon', 'synonymous_muts', 'total_muts') . "\n";
	foreach my $aa (sort keys %aa_to_codons)
	{
		foreach my $codon (sort @{$aa_to_codons{$aa}})
		{
			print OUT +join("\t", $aa, $codon, $codon_num_synonymous_changes{$codon}, 9) . "\n";
		}
	}
	print OUT "\n";

	print OUT "#total_nt_positions\t$total_nt_position\n";
	print OUT "#total_aa_coding_nt_positions\t$total_codon_nt_positions\n";
	print OUT "\n";

	my $total_possible_mutations = $total_nt_position * 3;
	print OUT "#total_possible_mutations $total_possible_mutations\n";
	my $total_possible_aa_coding_mutations = $total_codon_nt_positions * 3;
	print OUT "#total_possible_aa_coding_mutations $total_possible_aa_coding_mutations\n";

	print OUT "\n";
	print OUT "#total_num_synonymous_changes\t$total_num_synonymous_changes\n";
	print OUT "#total_num_nonsynonymous_changes\t$total_num_nonsynonymous_changes\n";

	print OUT "\n";
	my $chance_of_genome_mutation_synonymous = $total_num_synonymous_changes / $total_possible_mutations;
	print OUT "#chance_of_genome_mutation_synonymous\t$chance_of_genome_mutation_synonymous\n";
	my $chance_of_aa_coding_mutation_synonymous = $total_num_synonymous_changes / $total_possible_aa_coding_mutations;
	print OUT "#chance_of_aa_coding_mutation_synonymous\t$chance_of_aa_coding_mutation_synonymous\n";

	print OUT "\n";
	print OUT "#Base distribution for entire genome\n";
	print OUT +join("\t", 'base', 'num') . "\n";
	foreach my $base (keys %$total_bases)
	{
		print OUT +join("\t", $base, $total_bases->{$base}) . "\n";
	}

	print OUT "\n";
	print OUT "#Probabilities of synonymous mutations given base change over entire genome (NCBI Translation Table $translation_table)\n";
	print OUT +join("\t", "mutation", "fr_synonymous", "fr_nonsynonymous") . "\n";
	foreach my $key (sort keys %nonsynonymous_mutations)
	{
		next if ($nonsynonymous_mutations{$key} + $synonymous_mutations{$key} == 0);
		my $total = $nonsynonymous_mutations{$key} + $synonymous_mutations{$key};
		print OUT +join("\t",  $key, $synonymous_mutations{$key} / $total, $nonsynonymous_mutations{$key} / $total) . "\n";
	}

	## print out probabilities of all mutations GIVEN a synonymous mutation
	print OUT "#Probabilities of base changes GIVEN synonymous mutation (NCBI Translation Table $translation_table)\n";
	print OUT +join("\t", 'from-bp', 'to-bp', 'probability') . "\n";

	#print Dumper(\%nonsynonymous_mutations);
	#print Dumper(\@bp_change_list);

	foreach (my $i = 0; $i < scalar @bp_change_list; $i++)
	{
		my ($from_bp_change, $to_bp_change) = ($bp_change_list[$i]->[0], $bp_change_list[$i]->[1]);
		my $probability = ($synonymous_mutations{$from_bp_change} + $synonymous_mutations{$to_bp_change} ) / $total_num_synonymous_changes;	
		my ($from_bp, $to_bp) = ($bp_change_label_list[$i]->[0], $bp_change_label_list[$i]->[1]);
		print OUT +join("\t", $from_bp, $to_bp, $probability) . "\n";
	}
	print OUT "\n";

=cut

sub write 
{
	my ($self, $output_file) = @_;
	
	open OUT, ">$output_file";
	print OUT "$self->{packed_data}";
	close OUT;
}

sub totals
{
	my ($self) = @_;
	
	my $seq_length = length $self->{nt_sequence};
	my $totals;
	for (my $i=1; $i <= $seq_length; $i++)
	{
		$totals = $self->add_position_1_to_totals($totals, $i);
		#print Dumper($totals);
		print "$i\n" if ($i % 10000 == 0);
	}
	return $totals;
}


sub add_position_1_to_totals
{
	my ($self, $totals, $on_pos) = @_;
	#print "$on_pos\n";
	#print Dumper($totals);
	$totals = $self->change_position_1_totals($totals, $on_pos, +1);
	return $totals;
}

sub subtract_position_1_from_totals
{
	my ($self, $totals, $on_pos) = @_;
	$totals = $self->change_position_1_totals($totals, $on_pos, -1);
	return $totals;
}

sub change_position_1_totals
{
	my ($self, $totals, $on_pos, $inc) = @_;
	
	#print Dumper($totals);
	#print "$on_pos, $inc\n";
	
	my $pos_info = $self->pos_info_1($on_pos);

	#print Dumper($pos_info);

	$totals->[$pos_info->{nt_type}]->{nt} += $inc;
	if ($pos_info->{is_GC}) {
		$totals->[$pos_info->{nt_type}]->{GC} += $inc;
	} else {
		$totals->[$pos_info->{nt_type}]->{AT} += $inc;
	}
	
	foreach my $bp_mutation (@{$pos_info->{bp_mutations}})
	{
		$totals->[$pos_info->{nt_type}]->{$bp_mutation} += $inc;
		$totals->[$pos_info->{nt_type}]->{TOTAL} += $inc;
		
		##keep track of total in each category
		$totals->[6]->{$bp_mutation} += $inc;
		$totals->[6]->{TOTAL} += $inc;
	}
	
	## extra things to keep track of for proteins
	if ($pos_info->{nt_type} == 0)
	{
		for (my $j=0; $j< scalar @{$pos_info->{bp_mutations}}; $j++)
		{
			##NS
			$totals->[($pos_info->{bp_mutation_nonsynonymous}->[$j]) ? 4 : 5]->{$pos_info->{bp_mutations}->[$j]} += $inc;
			$totals->[($pos_info->{bp_mutation_nonsynonymous}->[$j]) ? 4 : 5]->{TOTAL} += $inc;
		}
	}
	
	#print Dumper($totals);
	return $totals;
}


sub add_bp_change_to_totals
{
	my ($self, $totals, $on_pos, $bp_change) = @_;
	my $inc = +1;
	my $pos_info = $self->pos_info_1($on_pos);
	
	$totals->[$pos_info->{nt_type}]->{$bp_change} += $inc;
	$totals->[$pos_info->{nt_type}]->{TOTAL} += $inc;
	
	##keep track of total in each category
	$totals->[6]->{$bp_change} += $inc;
	$totals->[6]->{TOTAL} += $inc;

	## extra things to keep track of for proteins
	if ($pos_info->{nt_type} == 0)
	{
		for (my $j=0; $j< scalar @{$pos_info->{bp_mutations}}; $j++)
		{
			if ($pos_info->{bp_mutations}->[$j] eq $bp_change)
			{
				$totals->[($pos_info->{bp_mutation_nonsynonymous}->[$j]) ? 4 : 5]->{$bp_change} += $inc;
				$totals->[($pos_info->{bp_mutation_nonsynonymous}->[$j]) ? 4 : 5]->{TOTAL} += $inc;
			}
		}
	}
	
	return $totals;				
}

sub pos_info_1
{
	my ($self, $i) = @_;
	my $pos = $i-1;
	my $pos_info;
	
	my $base = substr $self->{nt_sequence}, $pos, 1;
	my $is_GC = (($base eq 'G') || ($base eq 'C')) ? 1 : 0;
	$pos_info->{base} = $base;
	$pos_info->{is_GC} = $is_GC;
	
	my $offset = 3*$is_GC;
	@{$pos_info->{bp_mutations}} = @bp_change_label_list[$offset..$offset+2];
	
	#print "Length of packed data: " . +(length $self->{packed_data}) . "\n";
	
	my $encoded_char = substr $self->{packed_data}, $pos, 1;
	my ($nt_type_1, $nt_type_2, @is_nonsynonymous) =
	      split( //, unpack( 'b8', $encoded_char ) );
	
	#print "$nt_type_1, $nt_type_2, @is_nonsynonymous\n";
	
	$pos_info->{nt_type} = $nt_type_1 + 2*$nt_type_2;
	
	for (my $j=0; $j< scalar @{$pos_info->{bp_mutations}}; $j++)
	{
		$pos_info->{bp_mutation_nonsynonymous}->[$j] = $is_nonsynonymous[$j];
	}
	
	return $pos_info;
}


sub GetTag
{
	my ($Feature, $Tag, $Allow_Array) = @_;	
	$Allow_Array = 0 if (!defined $Allow_Array);
	
	return '' if (!$Feature->has_tag("$Tag"));
	my @Tag_Values = $Feature->get_tag_values("$Tag");
	return '' if (scalar @Tag_Values == 0);
	return $Tag_Values[0] if (!$Allow_Array);
	return @Tag_Values;
}	


return 1;
