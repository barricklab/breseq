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

package Breseq::ReferenceSequence;
use vars qw(@ISA);
@ISA = qw( Bio::Root::Root );

use Bio::Seq::RichSeq;
use Breseq::Settings;
use Breseq::Shared;
use Data::Dumper;

sub process_reference_sequences
{
	my ($settings, $summary) = @_;
	my $s = $summary->{sequence_conversion};
	
	# get two pieces of information from $settings
	my $reference_fasta_file_name = $settings->file_name('reference_fasta_file_name');
	my @genbank_file_names = $settings->file_name('reference_genbank_file_names'); 
		
	print STDERR "Loading reference sequences...\n";

	##list of sequence ids
	my @seq_ids;

	##open fasta file
	my $ref_o = Bio::SeqIO->new( -file => ">$reference_fasta_file_name", -format => 'fasta');

	my %ref_seqs;
	my %ref_strings;
	my %gene_lists;
	my %is_lists;
	my %fasta_file_names;
	my %seq_order;
	my $i = 0;

	foreach my $genbank_file_name (@genbank_file_names)
	{
		## open this GenBank file
		my $ref_i = Bio::SeqIO->new( -file => "<$genbank_file_name");
		print STDERR "  Loading File::$genbank_file_name\n";

		while (my $ref_seq = $ref_i->next_seq)
		{
			my $seq_id = $ref_seq->id;
			push @seq_ids, $seq_id;
			
			print STDERR "    Sequence::$seq_id loaded.\n";
			$ref_seqs{$seq_id} = $ref_seq;

			$ref_o->write_seq($ref_seq);

			##it is much faster to use substr to create the lists for nt comparisons than BioPerl trunc
			$ref_strings{$seq_id} = $ref_seq->seq;
			$ref_strings{$seq_id} = "\U$ref_strings{$seq_id}"; ##uppercase for comparisons

			##load the genbank record
			my @Feature_List = $ref_seq->get_SeqFeatures();
			my @gene_list;
			my @is_list;
			
			FEATURE: foreach my $Feature ( @Feature_List ) 
			{ 	
				if ($Feature->primary_tag eq 'repeat_region')
				{	
					my $is;
					$is->{gene} = get_tag($Feature, "mobile_element");
					$is->{gene} =~ s/insertion sequence:// if ($is->{gene});
					$is->{gene} = "unknown" if (!defined ($is->{gene}));
					$is->{product} = "repeat region";

					$is->{start} = $Feature->start;
					$is->{end} = $Feature->end;
					$is->{strand} = $Feature->strand;
					push @is_list, $is;
					next FEATURE;
				}
				
				## add additional information to the last 
				if (   ($Feature->primary_tag eq 'CDS') 
					|| ($Feature->primary_tag eq 'tRNA') 
					|| ($Feature->primary_tag eq 'rRNA') )

				{
					#Add information to last gene
					my $gene;
					$gene->{gene} = get_tag($Feature, "gene");
					$gene->{gene} = get_tag($Feature, "locus_tag") if (!$gene->{gene});
					$gene->{start} = $Feature->start;
					$gene->{end} = $Feature->end;
					$gene->{strand} = $Feature->strand;
					$gene->{product} = "";
					$gene->{note} = get_tag($Feature, "note");
						
					$gene->{accession} = get_tag($Feature, "protein_id");
					$gene->{translation} = get_tag($Feature, "translation");
					$gene->{product} = get_tag($Feature, "product");
					
					#set a type for the feature
					$gene->{type} = $Feature->primary_tag;
					$gene->{type} = "protein" if ($gene->{type} eq 'CDS');
					$gene->{pseudogene} = get_tag($Feature, "pseudo");
					$gene->{type} = "pseudogene" if ($gene->{pseudogene});
					
					##assume if there is no translation that we have a pseudogene...
					$gene->{type} = "pseudogene" if (($Feature->primary_tag eq 'CDS') && (!$gene->{translation}));
					$gene->{index} = scalar @gene_list;
					
					push @gene_list, $gene;		
				}
			}
			
			$gene_lists{$seq_id} = \@gene_list;
			$is_lists{$seq_id} = \@is_list;	
			
			#add information to summary
			$s->{reference_seqs}->{$seq_id}->{length} = $ref_seq->length;
			$s->{reference_seqs}->{$seq_id}->{display_id} = $ref_seq->display_id;
			$s->{reference_seqs}->{$seq_id}->{accession} = $ref_seq->accession;
			$s->{reference_seqs}->{$seq_id}->{accession} .= "." . $ref_seq->version if ($ref_seq->version);
			$seq_order{$seq_id} = $i++;
			}
	}
	
	## create SAM faidx
	Breseq::Shared::system("samtools faidx $reference_fasta_file_name", 1);
		
	return {'bioperl_ref_seqs' => \%ref_seqs, 'ref_strings' => \%ref_strings, 'gene_lists' =>\%gene_lists, 'is_lists' =>\%is_lists, 'seq_ids' => \@seq_ids, 'seq_order' => \%seq_order };	
}


sub annotate_mutations
{
	my ($settings, $summary, $ref_seq_info, $mutations_list_ref) = @_;

	my $current_line;
	SNP: foreach my $c (@$mutations_list_ref)
	{			
		#set up to work with this sequence
		my $seq_id = $c->{seq_id};
		my $gene_list_ref = $ref_seq_info->{gene_lists}->{$seq_id};
		my $is_list_ref = $ref_seq_info->{is_lists}->{$seq_id};
		my $ref_seq = $ref_seq_info->{bioperl_ref_seqs}->{$seq_id};
								
		$c->{type} = "";
		$c->{aa_position} = "";
		$c->{aa_ref_seq} = "";
		$c->{aa_new_seq} = "";
		$c->{codon_position} = "";
		$c->{codon_ref_seq} = "";
		$c->{codon_new_seq} = "";		
		$c->{gene} = "";
		$c->{gene_position} = "";
		$c->{gene_product} = "";

#		#save original start and end before we possibly shift them around
		$c->{original_start} = $c->{start};
		$c->{original_end} = $c->{end};
		
		if (defined $c->{quality})
		{		
			my @quality_list = split /,/, $c->{quality};
			@quality_list = map { sprintf("%.1f", $_) } @quality_list;
			$c->{display_quality} = join ",", @quality_list;
		}
		
		## Create versions of each sequence that can be substituted into
		## reading frames to determine alternate translations.
		$c->{ref_rep_seq} = $c->{ref_seq};
		$c->{ref_rep_seq} = '' if ($c->{ref_rep_seq} eq '.');
		
		$c->{new_rep_seq} = $c->{new_seq};
		$c->{new_rep_seq} = '' if ($c->{new_rep_seq} eq '.');	
			
		## Determine whether the change is within a gene
		my ($prev_gene, $gene, $next_gene) = find_nearby_genes($c, $gene_list_ref);		
		
		## The change is within a gene on the bottom strand
		## For INDELS, move the change to be as late on the gene's strand as possible
		## i.e. if we deleted a G in a run of GGGGGGG
		## This may actually shift it out of the reading frame! 
		## (Clean up for this possibility happens afterwards.)
		
		if ($gene && ($gene->{strand} == - 1)) 
		{
			## deletion 
			if (length($c->{ref_rep_seq}) > length($c->{new_rep_seq}))
			{
				my $len = length($c->{ref_rep_seq}) - length($c->{new_rep_seq});
				if ($c->{start}-$len >= 1)
				{
					my $test_shift_sequence = $ref_seq->trunc($c->{start}-$len, $c->{start}-1)->seq;
					#print STDERR "$c->{start} $len $test_shift_sequence ?= $c->{ref_rep_seq} -> $c->{new_rep_seq}\n";

					MOVE: while ( "$test_shift_sequence" eq $c->{ref_rep_seq})
					{
						if ($c->{start} - $len < 1)
						{
							last MOVE;
						}
						$c->{start} -= $len;
						$c->{end} -= $len;
						
						$test_shift_sequence = $ref_seq->trunc($c->{start}-$len, $c->{start}-1)->seq;
						#print STDERR "$c->{start} $len $test_shift_sequence ?= $c->{ref_rep_seq} -> $c->{new_rep_seq}\n";
					}
				}
			}
			## insertion
			elsif (length($c->{new_rep_seq}) > length($c->{ref_rep_seq}))
			{
				my $len = length($c->{new_rep_seq}) - length($c->{ref_rep_seq});
				if ($c->{start}-$len >= 1)
				{
					my $test_shift_sequence = $ref_seq->trunc($c->{start}-$len, $c->{start}-1)->seq;
					#print STDERR "$c->{start} $len $test_shift_sequence ?= $c->{ref_rep_seq} -> $c->{new_rep_seq}\n";
					MOVE: while ("$test_shift_sequence" eq $c->{new_rep_seq})
					{
						if ($c->{start} - $len < 1)
						{
							last MOVE;
						}
						$c->{start} -= $len;
						$c->{end} -= $len;

						$test_shift_sequence = $ref_seq->trunc($c->{start}-$len, $c->{start}-1)->seq;
						#print STDERR "$c->{start} $len $test_shift_sequence ?= $c->{ref_rep_seq} -> $c->{new_rep_seq}\n";
					}
				}
			}
			
			## redetermine the nearby genes (it may now be intergenic)
			($prev_gene, $gene, $next_gene) = find_nearby_genes($c, $gene_list_ref);		
		}		
		
		
		
		###
		### Change is noncoding, $prev_gene and/or $next_gene will be defined
		###
		
		if (!$gene)
		{			
			$c->{type} = "noncoding";
	
			$c->{gene} .= $prev_gene->{gene} if (defined $prev_gene);
			$c->{gene} .= "/";
			$c->{gene} .= $next_gene->{gene} if (defined $next_gene);
	
			if (defined $prev_gene)
			{
				$c->{gene_position} .= ($prev_gene->{strand} == +1) ? "+" : "-";
				$c->{gene_position} .= ($c->{start} - $prev_gene->{end});
			}
			$c->{gene_position} .= "/";
			if (defined $next_gene)
			{
				$c->{gene_position} .= ($next_gene->{strand} == +1) ? "-" : "+";
				$c->{gene_position} .= ($next_gene->{start} - $c->{end});
			}
			
			$c->{gene_product} .= $prev_gene->{product} if (defined $prev_gene);
			$c->{gene_product} .= "/";			
			$c->{gene_product} .= $next_gene->{product} if (defined $next_gene);

			#it's a SNP
			if ((length($c->{ref_rep_seq}) == 1) && (length($c->{new_rep_seq}) == 1))
			{
				$summary->{snps}->{noncoding}->{num}++;
				$summary->{snps}->{noncoding}->{base_changes}->{"$c->{ref_rep_seq}$c->{new_rep_seq}"}++;
				$summary->{snps}->{total}->{num}++;
				$summary->{snps}->{total}->{base_changes}->{"$c->{ref_rep_seq}$c->{new_rep_seq}"}++;
			}
			next SNP;
		}
		
		###
		### Change is coding
		###
		
		$c->{gene} = $gene->{gene};
		$c->{gene_product} = $gene->{product};

		#remember to revcom the changes if gene is on opposite strand
		$c->{ref_rep_seq} = Breseq::Fastq::revcom($c->{ref_rep_seq}) if ($gene->{strand} == -1);		
		$c->{new_rep_seq} = Breseq::Fastq::revcom($c->{new_rep_seq}) if ($gene->{strand} == -1);
		
		### Single base substitution -- SNP-type
		if (($gene->{translation}) && (length($c->{ref_rep_seq}) == 1) && (length($c->{new_rep_seq}) == 1))
		{
		    ## determine the old and new translation of this codon  
			my $gene_start = ($gene->{strand} == +1) ? $gene->{start} : $gene->{end};	
			$c->{gene_position} = abs($c->{start}-$gene_start) + 1; 
    		$c->{aa_position} = int(($c->{gene_position}-1)/3) + 1;
			$c->{codon_position} = abs($c->{start}-$gene_start) % 3;
	    	
			my $codon_seq = ($gene->{strand} == +1) ?
				$ref_seq->trunc($gene->{start} + 3 * ($c->{aa_position}-1),$gene->{start} + 3 * $c->{aa_position} - 1) :
				$ref_seq->trunc($gene->{end} - 3 * $c->{aa_position}+1,$gene->{end} - 3 * ($c->{aa_position}-1))->revcom;
					
			$c->{codon_ref_seq} = $codon_seq->seq();
			$c->{aa_ref_seq} = $codon_seq->translate( undef, undef, undef, 11 )->seq();
			
			$c->{codon_new_seq} = $codon_seq->seq();
			substr($c->{codon_new_seq}, $c->{codon_position}, 1) = $c->{new_rep_seq};
			$codon_seq->seq($c->{codon_new_seq});
			$c->{aa_new_seq} = $codon_seq->translate( undef, undef, undef, 11 )->seq();
		    
		    $c->{type} = "substitution";
		    $c->{snp_type} = ($c->{aa_ref_seq} ne $c->{aa_new_seq}) ? "nonsynonymous" : "synonymous";
			
			## record for dN/dS statistics
			$summary->{snps}->{$c->{snp_type}}->{num}++;
			$summary->{snps}->{$c->{snp_type}}->{base_changes}->{"$c->{ref_rep_seq}$c->{new_rep_seq}"}++;
			$summary->{snps}->{total}->{num}++;
			$summary->{snps}->{total}->{base_changes}->{"$c->{ref_rep_seq}$c->{new_rep_seq}"}++;
			next SNP;
		}
	
	
		if (!$gene->{translation})
		{
			$c->{type} = "pseudogene";
		}
		
		#Deletions, insertions, and frameshifts
		
		## Change in length does not result in frameshift, retranslate area
		if (abs(length($c->{new_rep_seq}) - length($c->{ref_rep_seq})) % 3 == 0)
		{
			$c->{type} = "in-frame indel";
			
			#need to round before and after sequences to whole codons
		
#			## Insertion of amino acids
#			my ($ref_codon_seq, $new_codon_seq);		
#			$c->{aa_position} = ($gene->{strand} == +1) ? int(($c->{start}-$gene->{start})/3) + 1 : int(($gene->{end}-$c->{end})/3) + 1;
#			$c->{aa_position} = "";
#			
#			## translate to the inserted or deleted amino acid sequence
#			$ref_codon_seq =  Bio::Seq->new( '-seq' => $c->{ref_rep_seq}) if (length $c->{ref_rep_seq} > 0);
#			$new_codon_seq =  Bio::Seq->new( '-seq' => $c->{new_rep_seq}) if (length $c->new_rep_seq} > 0);
#
#			my ($ref_aa, $new_aa) = ('.','.');
#			$c->{aa_ref_seq} = $ref_codon_seq->translate( undef, undef, undef, 11 )->seq if ($ref_codon_seq);
#			$c->{aa_new_seq} = $new_codon_seq->translate( undef, undef, undef, 11 )->seq if ($new_codon_seq);
#			
#			$c->{type} = (length $new_rep_seq > length $ref_rep_seq) ? "insertion" : "deletion";
#			
#			
#			
			next SNP;
		} 
		## Change results in a frameshift
		my ($codon_num, $ref_reading_frame, $new_reading_frame, $ref_orf_seq, $new_orf_seq);
		$c->{type} = "frameshift";		
		
		my $gene_start = ($gene->{strand} == +1) ? $gene->{start} : $gene->{end};	
		$c->{gene_position} = abs($c->{start}-$gene_start) + 1; 				
						
		if ($gene->{strand} == +1)
		{
			$c->{aa_position} = int(($c->{start}-$gene->{start})/3) + 1;
			
			$ref_reading_frame = $ref_seq->trunc($gene->{start}, $gene->{end});
			$ref_orf_seq = $ref_reading_frame->seq();
			$new_orf_seq = substr($ref_orf_seq, 0, $c->{start}-$gene->{start}+1) . $c->{new_rep_seq} 
				. substr($ref_orf_seq, $c->{start}+length($c->{ref_rep_seq} )-$gene->{start}+1, $gene->{end}-($c->{start}+length($c->{ref_rep_seq} )));
	
			#extend the reading frame so we can find the new protein length (and eventually sequence)
			if ($gene->{end}+1 <= $ref_seq->length)
			{
				my $safe_end = $gene->{end}+10000;
				$safe_end = $ref_seq->length if ($safe_end > $ref_seq->length);
				$new_orf_seq .=  $ref_seq->trunc($gene->{end}+1, $safe_end)->seq;
			}
						
			##ok, finally have the new reading frame
			$new_reading_frame = Bio::Seq->new( '-seq' => $new_orf_seq);
		}
		else #Reverse strand
		###
		### Add code to shift differences over so they show up as early in the reading frame as possible!!!
		### Be sure that this works for insertions with length > 1 as well.
		###
		
		{
			$c->{aa_position} = int(($gene->{end}-$c->{end})/3) + 1;
			
			$ref_reading_frame = $ref_seq->trunc($gene->{start}, $gene->{end});
			$ref_orf_seq = $ref_reading_frame->seq();
			$new_orf_seq = substr($ref_orf_seq, 0, $c->{start}-$gene->{start}+1) . $c->{new_rep_seq}
				. substr($ref_orf_seq, $c->{start}+length($c->{ref_rep_seq})-$gene->{start}+1, $gene->{end}-($c->{start}+length($c->{ref_rep_seq})));

			#extend the reading frame so we can find the new protein length (and eventually sequence)
			if ($gene->{start}-1 >= 1)
			{
				my $safe_start = $gene->{start}-10000;
				$safe_start = 1 if ($safe_start <= 1);
				$new_orf_seq =  $ref_seq->trunc($safe_start, $gene->{start}-1)->seq . $new_orf_seq;
			}
			$new_reading_frame = Bio::Seq->new( '-seq' => $new_orf_seq);

			$ref_reading_frame = $ref_reading_frame->revcom();
			$new_reading_frame = $new_reading_frame->revcom();
		}		
		
		## determine how the length of the ORF changes
		my $ref_protein = $ref_reading_frame->translate( undef, undef, undef, 11 )->seq;
	#	die "Did not find a stop codon in reference protein" if (! ($ref_protein =~ s/\*.+/\*/) );		
		my $ref_length = length($ref_protein) - 1;

		my $new_protein = $new_reading_frame->translate( undef, undef, undef, 11 )->seq;
		die "Did not find a stop codon in new protein" if (! ($new_protein =~ s/\*.+/\*/) );		
		my $new_length = length($new_protein) - 1;

		$c->{aa_ref_seq} = substr $ref_protein, $c->{aa_position}-1, 1;
		$c->{aa_new_seq} = substr $new_protein, $c->{aa_position}-1, 1; 	

		$c->{codon_ref_seq} = substr $ref_reading_frame->seq, ($c->{aa_position}-1) * 3, 3;
		$c->{codon_new_seq} = substr $new_reading_frame->seq, ($c->{aa_position}-1) * 3, 3; 
	} continue {
		$c->{gene_shifted_start} = $c->{start};
		$c->{gene_shifted_end} = $c->{end};
		$c->{start} = $c->{original_start};
		$c->{end} = $c->{original_end};
	}
}

sub annotate_deletions
{
	my ($settings, $summary, $ref_seq_info, $deletions_ref) = @_;	
	
	foreach my $deletion (@$deletions_ref)
	{
		#set up to use correct reference sequence
		my $seq_id = $deletion->{seq_id};
		my $gene_list_ref = $ref_seq_info->{gene_lists}->{$seq_id};
		my $is_list_ref = $ref_seq_info->{is_lists}->{$seq_id};
		
		my $gene_string = '';
		my $feature_start = get_overlapping_feature($deletion->{start}, $is_list_ref);
		$feature_start = get_overlapping_feature($deletion->{start}, $gene_list_ref) if (!$feature_start);
	
		my $feature_end = get_overlapping_feature($deletion->{end}, $is_list_ref);
		$feature_end = get_overlapping_feature($deletion->{end}, $gene_list_ref) if (!$feature_end);
	
		if ($feature_start)
		{
			$gene_string .= "[$feature_start->{gene}]";
		}
	
		my $start_bound = ($feature_start) ? $feature_start->{end}-25 : $deletion->{start};
		my $end_bound = ($feature_end) ? $feature_end->{start} : $deletion->{end};

		foreach my $gene (@$gene_list_ref)
		{
			last if ($gene->{start} >= $end_bound);
		
			if ($gene->{start} > $start_bound)
			{
				$gene_string .= " " if ($gene_string);
				$gene_string .= "$gene->{gene}";
			}
		}
	
		## only add the end gene if it isn't the start gene
		if ($feature_end && (!$feature_start || ($feature_start != $feature_end)))
		{
			$gene_string .= " " if ($gene_string);
			$gene_string .= "[$feature_end->{gene}]";
		}
	
		##no genes found
		$gene_string .= "noncoding" if (!$gene_string);
		
		$deletion->{genes} = $gene_string;
	}
}


sub annotate_rearrangements
{
	my $verbose = 0;
	our ($settings, $summary, $ref_seq_info, $rearrangements_ref) = @_;	
	
	foreach my $item (@$rearrangements_ref)
	{				
		## This additional information is used for the complex reference line
		my $alignment_reference_info_1 = { 
			truncate_end => $item->{flanking_length}, 
			ghost_end => $item->{interval_1}->{start}, 
			ghost_strand => $item->{interval_1}->{strand},
			ghost_seq_id => $item->{interval_1}->{seq_id}
		};

		my $alignment_reference_info_2 = { 
			truncate_start => $item->{flanking_length}+1-$item->{overlap}, 
			ghost_start => $item->{interval_2}->{start}, 
			ghost_strand => $item->{interval_2}->{strand},
			ghost_seq_id => $item->{interval_2}->{seq_id}
		};

		push @{$item->{alignment_reference_info_list}}, $alignment_reference_info_1, $alignment_reference_info_2;
		$item->{alignment_empty_change_line} = 1;

		#add gene information for each end
		$item->{hybrid} = $item;
		foreach my $key ('interval_1', 'interval_2')
		{			
			##create circular reference to self so we can print table at the top of the alignment
			$item->{$key}->{hybrid} = $item;
			
			my ($prev_gene, $gene, $next_gene) 
				= Breseq::ReferenceSequence::find_nearby_genes($item->{$key}, $ref_seq_info->{gene_lists}->{$item->{$key}->{seq_id}});		

			## noncoding
			if (!$gene)
			{
				$item->{$key}->{gene}->{gene} .= ($prev_gene && $prev_gene->{gene}) ? $prev_gene->{gene} : '-';
				$item->{$key}->{gene}->{gene} .= "/";
				$item->{$key}->{gene}->{gene} .= ($next_gene && $next_gene->{gene}) ? $next_gene->{gene} : '-';
	
				if (defined $prev_gene)
				{
					$item->{$key}->{gene}->{position} .= ($prev_gene->{strand} == +1) ? "+" : "-";
					$item->{$key}->{gene}->{position} .= ($item->{$key}->{start} - $prev_gene->{end});
				}
				$item->{$key}->{gene}->{gene_position} .= "/";
				if (defined $next_gene)
				{
					$item->{$key}->{gene}->{position} .= ($next_gene->{strand} == +1) ? "-" : "+";
					$item->{$key}->{gene}->{position} .= ($next_gene->{start} - $item->{$key}->{end});
				}
		
				$item->{$key}->{gene}->{product} .= ($prev_gene && $prev_gene->{product}) ? $prev_gene->{product} : '-';
				$item->{$key}->{gene}->{product} .= "/";			
				$item->{$key}->{gene}->{product} .= ($next_gene && $next_gene->{product}) ? $next_gene->{product} : '-';

				$item->{$key}->{gene}->{interval} .= $prev_gene->{end}+1 if $prev_gene;
				$item->{$key}->{gene}->{interval} .= "/"; 
				$item->{$key}->{gene}->{interval} .= $next_gene->{start}-1 if $next_gene; 
			}
			## coding
			else
			{
				$item->{$key}->{gene}->{gene} = $gene->{gene};
				$item->{$key}->{gene}->{product} = $gene->{product};
				my $gene_start = ($gene->{strand} == +1) ? $gene->{start} : $gene->{end};	
				$item->{$key}->{gene}->{position} = abs($item->{$key}->{start}-$gene_start) + 1;
				$item->{$key}->{gene}->{interval} = ($gene->{strand} == +1) ? "$gene->{start}-$gene->{end}" : "$gene->{end}-$gene->{start}"; 
 
			}
		
			## determine IS elements
			## Is it within an IS or near the boundary of an IS in the direction leading up to the junction?			
			if (my $is = Breseq::ReferenceSequence::find_closest_is_element($item->{$key}, $ref_seq_info->{is_lists}->{$item->{$key}->{seq_id}}, 200, $item->{$key}->{strand}))
			{
				$item->{$key}->{is}->{gene} = $is->{gene};
				$item->{$key}->{is}->{interval} = ($is->{strand} == +1) ? "$is->{start}-$is->{end}" : "$is->{end}-$is->{start}"; 
				$item->{$key}->{is}->{product} = $is->{product};
			}
			$item->{$key}->{annotate_key} = (defined $item->{$key}->{is}) ? 'is' : 'gene';			
		}
		print Dumper($item) if ($verbose);
	}

	sub by_hybrid
	{
		my $a_pos = (defined $a->{interval_1}->{is}) ? $a->{interval_2}->{start} : $a->{interval_1}->{start};
		my $b_pos = (defined $b->{interval_1}->{is}) ? $b->{interval_2}->{start} : $b->{interval_1}->{start};
	
		my $a_seq_order = (defined $a->{interval_1}->{is}) ? $ref_seq_info->{seq_order}->{$a->{interval_2}->{seq_id}} : $ref_seq_info->{seq_order}->{$a->{interval_1}->{seq_id}};
		my $b_seq_order = (defined $b->{interval_1}->{is}) ? $ref_seq_info->{seq_order}->{$b->{interval_2}->{seq_id}} : $ref_seq_info->{seq_order}->{$b->{interval_1}->{seq_id}};		
	
		return (($a_seq_order <=> $b_seq_order) || ($a_pos <=> $b_pos));
	}
	@$rearrangements_ref = sort by_hybrid @$rearrangements_ref;
}


sub get_overlapping_feature
{
	my ($pos, $feature_list_ref) = @_;
	foreach my $f (@$feature_list_ref)
	{
		return $f if ( ($pos >= $f->{start}) && ($pos <= $f->{end}) );
	}
	return undef;
}

sub find_nearby_genes
{
	my ($c, $gene_list_ref) = @_;

	my ($gene, $prev_gene, $next_gene);
	GENE: for (my $i=0; $i < scalar @$gene_list_ref; $i++)
	{
		my $test_gene = $gene_list_ref->[$i];
		#This change is within a gene
		if ( ($test_gene->{start} <= $c->{start}) && ($test_gene->{end} >= $c->{end}) )
		{	
			$gene = $test_gene;
			last GENE;
		}
		#We've passed the changes, so it is in the previous intergenic space
		if ($test_gene->{start} > $c->{end})
		{
			$prev_gene = $gene_list_ref->[$i-1] if ($i > 0);
			$next_gene = $test_gene;
			last GENE;
		}
	}
	$prev_gene = $gene_list_ref->[-1] if (!$gene && !$prev_gene && !$next_gene);
	return ($prev_gene, $gene, $next_gene);
}

sub find_closest_is_element
{
	my ($c, $is_list_ref, $max_distance, $direction) = @_;

	my $is;
	my $best_distance;
	IS: for (my $i=0; $i < scalar @$is_list_ref; $i++)
	{
		my $test_is = $is_list_ref->[$i];
				
		#count within the IS element as zero distance
		#if this happens then we are immediately done
		if ( ($test_is->{start} <= $c->{start}) && ($test_is->{end} >= $c->{end}) )
		{	
			return $test_is;
		}
		
		#otherwise calculate the distance
		#keep if less than max_distance, in the correct direction, and the best found so far 
		my $test_distance = ($direction == -1) ? $c->{start} - $test_is->{end} : $test_is->{start} - $c->{end};
		next if ($test_distance < 0); #wrong direction...
		
		if (($test_distance <= $max_distance) && ((!defined $is) || ($test_distance < $best_distance)) )
		{
			$is = $test_is;
			$best_distance = $test_distance;
		}
	}
	return $is;
}


sub get_tag
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